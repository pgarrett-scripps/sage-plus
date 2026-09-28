use crate::database::{binary_search_slice, IndexedDatabase, PeptideIx};
use crate::mass::{composition, Composition, Tolerance, NEUTRON};
use crate::ml::{matrix::Matrix, retention_alignment::Alignment};
use crate::scoring::Feature;
use crate::spectrum::ProcessedSpectrum;
use dashmap::DashMap;
use fnv::FnvHashSet;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Minimum normalized spectral angle required to integrate a peak
// const MIN_SPECTRAL_ANGLE: f64 = 0.70;
/// Width of gaussian kernel used for smoothing intensities
const K_WIDTH: usize = 10;
/// Mass tolerance, in ppm, to seach for precursor ions
// const PPM_TOL: f32 = 5.0;
/// Number of equally spaced bins that will be used to integrate ions in the
/// configured retention-time window
const GRID_SIZE: usize = 100;
/// Number of isotopes to search for
const N_ISOTOPES: usize = 3;
/// Bins either side of an identification RT searched for isotope-consistent
/// signal before climbing to the apex.
const APEX_SNAP_BINS: usize = 10;
/// Warp search range, in bins, for files without their own identification.
const WARP_SLACK: isize = 75;

fn default_rt_pct_tolerance() -> f32 {
    0.5
}

#[derive(
    Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, schemars::JsonSchema,
)]
pub enum PeakScoringStrategy {
    RetentionTime,
    SpectralAngle,
    Intensity,
    Hybrid,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize, schemars::JsonSchema)]
pub enum IntegrationStrategy {
    Apex,
    Sum,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PrecursorId {
    Combined(PeptideIx),
    Charged((PeptideIx, u8)),
}

#[derive(Copy, Clone, Debug, Deserialize, Serialize)]
pub struct LfqSettings {
    pub peak_scoring: PeakScoringStrategy,
    pub integration: IntegrationStrategy,
    pub spectral_angle: f64,
    pub ppm_tolerance: f32,
    /// Symmetric retention-time tolerance as a percentage of total run length.
    #[serde(default = "default_rt_pct_tolerance")]
    pub rt_pct_tolerance: f32,
    pub mobility_pct_tolerance: f32,
    pub combine_charge_states: bool,
    pub peptide_q_value: f32,
    /// Trace identified precursors into files without direct MS2 evidence.
    #[serde(default = "default_true")]
    pub mbr: bool,
    /// Center extraction on the MS1 elution apex instead of the MS2
    /// identification RT: with MBR the traced window is centered on the median
    /// identification RT across files, each identified file's apex is refined
    /// from its own isotope trace, the cross-run apex is the median of those
    /// apexes, and integration bounds follow the peak shape.
    #[serde(default = "default_recenter_on_apex")]
    pub recenter_on_apex: bool,
}

fn default_recenter_on_apex() -> bool {
    false
}

fn default_true() -> bool {
    true
}

impl Default for LfqSettings {
    fn default() -> Self {
        Self {
            peak_scoring: PeakScoringStrategy::Hybrid,
            integration: IntegrationStrategy::Sum,
            spectral_angle: 0.70,
            ppm_tolerance: 5.0,
            rt_pct_tolerance: default_rt_pct_tolerance(),
            mobility_pct_tolerance: 1.0,
            combine_charge_states: true,
            peptide_q_value: 0.01,
            mbr: true,
            recenter_on_apex: default_recenter_on_apex(),
        }
    }
}

impl LfqSettings {
    fn rt_tolerance(&self) -> f32 {
        self.rt_pct_tolerance / 100.0
    }
}

/// Retention time of the best identification of a precursor in one file.
#[derive(Copy, Clone, Debug)]
pub struct IdentificationRt {
    /// Consensus (aligned) RT.
    pub aligned: f32,
    /// RT in the file's own time units.
    pub observed: f32,
}

#[derive(Copy, Clone, Debug)]
pub struct PrecursorRange {
    pub rt: f32,
    pub mass_lo: f32,
    pub mass_hi: f32,
    pub mobility_lo: f32,
    pub mobility_hi: f32,
    pub charge: u8,
    pub isotope: usize,
    pub peptide: PeptideIx,
    pub file_id: usize,
    pub decoy: bool,
}

/// Create a data structure analogous to [`IndexedDatabase`] - instaed of
/// storing fragment masses binned by precursor mass, store MS1 precursors
/// binned by RT - This should enable rapid quantification as well
pub struct FeatureMap {
    pub ranges: Vec<PrecursorRange>,
    pub min_rts: Vec<f32>,
    pub bin_size: usize,
    pub settings: LfqSettings,
    mass_search_margin: f32,
    /// Direct MS2 evidence for each LFQ precursor in each acquisition file.
    /// LFQ intensities are always obtained through the cross-run feature-tracing
    /// workflow; this set records only whether a matching accepted PSM was
    /// observed in the file itself.
    ms2_confirmed: FnvHashSet<(PrecursorId, usize)>,
    ms2_confirmed_strict: FnvHashSet<(PrecursorId, usize)>,
    /// RT of the most confident accepted PSM of each precursor in each file.
    id_rts: fnv::FnvHashMap<(PrecursorId, usize), IdentificationRt>,
}

/// A quantified LFQ precursor across all acquisition files.
#[derive(Clone, Debug)]
pub struct QuantifiedPeak {
    pub peak: Peak,
    /// Integrated MS1 intensity for each acquisition file. `None` represents
    /// an absent integrated signal rather than a measured zero.
    pub intensities: Vec<Option<f64>>,
    /// Whether the corresponding file contains a matching accepted target PSM.
    pub ms2_confirmed: Vec<bool>,
    /// Direct target evidence passing both the PSM and peptide thresholds.
    pub ms2_confirmed_strict: Vec<bool>,
    /// Diagnostics for individual file signals. These are not calibrated probabilities.
    pub file_evidence: Vec<Option<FileEvidence>>,
    /// For a target, its shifted decoy evaluated at this target's peak (same
    /// warps, apex, window and reference trace) in each file. These are the
    /// competitors for per-file extraction q-values. Empty for decoys.
    pub paired_decoy_evidence: Vec<Option<FileEvidence>>,
}

#[derive(Clone, Debug, Default)]
pub struct FileEvidence {
    /// No strict direct target evidence in this file, including for shifted decoys.
    pub transfer_candidate: bool,
    pub spectral_angle: f64,
    pub trace_cosine: f64,
    pub rt_shift_bins: i32,
    /// Per-file ranking score: isotope agreement, trace similarity and warp
    /// proximity. It also ranks extractions for [`crate::fdr::extraction_q_values`].
    pub score: f64,
    /// Target-decoy q-value of this target (precursor, file) extraction; see
    /// [`crate::fdr::extraction_q_values`]. `None` for decoys and until assigned.
    pub extraction_q_value: Option<f32>,
    /// Apex of the integrated peak in this file, in the file's RT units.
    pub apex_rt: f32,
    /// Integration bounds in this file, in the file's RT units.
    pub peak_start_rt: f32,
    pub peak_end_rt: f32,
    /// Full width at half maximum of this file's trace around the apex, in
    /// the file's RT units; `None` when the trace does not fall to half
    /// height inside the traced window.
    pub fwhm: Option<f32>,
    /// Identification RT minus apex RT, for files with an accepted PSM.
    pub id_apex_offset: Option<f32>,
}

pub fn build_feature_map(
    settings: LfqSettings,
    precursor_charge: (u8, u8),
    features: &[Feature],
    db: &IndexedDatabase,
) -> FeatureMap {
    let rt_tol = settings.rt_tolerance();
    let ms2_confirmed = features
        .iter()
        .filter(|feat| feat.peptide_q <= settings.peptide_q_value && feat.label == 1)
        .map(|feat| {
            let id = if settings.combine_charge_states {
                PrecursorId::Combined(feat.peptide_idx)
            } else {
                PrecursorId::Charged((feat.peptide_idx, feat.charge))
            };
            (id, feat.file_id)
        })
        .collect::<FnvHashSet<_>>();
    let ms2_confirmed_strict = features
        .iter()
        .filter(|feat| {
            feat.label == 1
                && feat.peptide_q <= settings.peptide_q_value
                && feat.spectrum_q <= settings.peptide_q_value
        })
        .map(|feat| {
            let id = if settings.combine_charge_states {
                PrecursorId::Combined(feat.peptide_idx)
            } else {
                PrecursorId::Charged((feat.peptide_idx, feat.charge))
            };
            (id, feat.file_id)
        })
        .collect();
    // `features` is sorted by confidence, so the first accepted PSM of a
    // precursor in a file is its best identification there.
    let mut id_rts = fnv::FnvHashMap::default();
    let mut peptide_rts: fnv::FnvHashMap<PeptideIx, fnv::FnvHashMap<usize, f32>> =
        fnv::FnvHashMap::default();
    for feat in features
        .iter()
        .filter(|feat| feat.peptide_q <= settings.peptide_q_value && feat.label == 1)
    {
        let id = if settings.combine_charge_states {
            PrecursorId::Combined(feat.peptide_idx)
        } else {
            PrecursorId::Charged((feat.peptide_idx, feat.charge))
        };
        id_rts
            .entry((id, feat.file_id))
            .or_insert(IdentificationRt {
                aligned: feat.aligned_rt,
                observed: feat.rt,
            });
        peptide_rts
            .entry(feat.peptide_idx)
            .or_default()
            .entry(feat.file_id)
            .or_insert(feat.aligned_rt);
    }
    // With apex recentering and MBR, the traced window is centered on the
    // median identification RT across files rather than on the single most
    // confident PSM; the apexes themselves are refined from the traces.
    let anchor_rt = |feat: &Feature| -> f32 {
        if !(settings.recenter_on_apex && settings.mbr) {
            return feat.aligned_rt;
        }
        let mut rts = peptide_rts
            .get(&feat.peptide_idx)
            .map(|files| files.values().copied().collect::<Vec<_>>())
            .unwrap_or_default();
        if rts.is_empty() {
            return feat.aligned_rt;
        }
        rts.sort_unstable_by(f32::total_cmp);
        rts[(rts.len() - 1) / 2]
    };
    let map: DashMap<(PeptideIx, usize), PrecursorRange, fnv::FnvBuildHasher> = DashMap::default();
    let label_groups = db
        .peptides
        .iter()
        .enumerate()
        .filter(|(_, peptide)| !peptide.decoy && peptide.label_channel.is_some())
        .fold(
            HashMap::<String, Vec<PeptideIx>, fnv::FnvBuildHasher>::default(),
            |mut groups, (index, peptide)| {
                groups
                    .entry(peptide.label_group())
                    .or_default()
                    .push(PeptideIx(index as u32));
                groups
            },
        );
    features
        .iter()
        .filter(|feat| feat.peptide_q <= settings.peptide_q_value && feat.label == 1)
        .for_each(|feat| {
            let peptide = &db[feat.peptide_idx];
            let members = peptide
                .label_channel
                .as_ref()
                .and_then(|_| label_groups.get(&peptide.label_group()))
                .map(Vec::as_slice)
                .unwrap_or(std::slice::from_ref(&feat.peptide_idx));
            for peptide_idx in members {
                let map_file_id = if settings.mbr {
                    usize::MAX
                } else {
                    feat.file_id
                };
                // `features` is sorted by confidence, so take the first anchor
                // for each exact channel precursor and requested file scope.
                if map.contains_key(&(*peptide_idx, map_file_id)) {
                    continue;
                }
                // let mass = if feat.isotope_error > 0.0 || feat.delta_mass >= settings.ppm_tolerance * 3.0 {
                //    feat.expmass - feat.isotope_error
                // } else {
                //     feat.calcmass
                // };
                let (mobility_lo, mobility_hi) = Tolerance::Pct(
                    -settings.mobility_pct_tolerance,
                    settings.mobility_pct_tolerance,
                )
                .bounds(feat.ims);
                map.insert(
                    (*peptide_idx, map_file_id),
                    PrecursorRange {
                        rt: anchor_rt(feat),
                        mass_lo: db[*peptide_idx].monoisotopic,
                        mass_hi: 0.0,
                        peptide: *peptide_idx,
                        charge: feat.charge,
                        isotope: 0,
                        file_id: feat.file_id,
                        mobility_lo,
                        mobility_hi,
                        decoy: false,
                    },
                );
            }
        });

    // Cartesian product of (observed RT ranges, mass) with charge and isotopes
    let mut ranges = map
        .into_par_iter()
        .flat_map_iter(|(_, range)| {
            (precursor_charge.0..=precursor_charge.1).flat_map(move |charge| {
                (0..N_ISOTOPES).flat_map(move |isotope| {
                    let mass = (range.mass_lo + isotope as f32 * NEUTRON) / charge as f32;
                    let (mass_lo, mass_hi) =
                        Tolerance::Ppm(-settings.ppm_tolerance, settings.ppm_tolerance)
                            .bounds(mass);

                    let fwd = PrecursorRange {
                        mass_lo,
                        mass_hi,
                        charge,
                        isotope,
                        decoy: false,
                        ..range
                    };

                    let (mass_lo, mass_hi) =
                        Tolerance::Ppm(-settings.ppm_tolerance, settings.ppm_tolerance)
                            .bounds(mass + 11.06);

                    // Shift the decoy by a full traced window, so that its
                    // window does not overlap the target's.
                    let rev = PrecursorRange {
                        rt: (fwd.rt - rt_tol * 2.0).max(0.0),
                        mass_lo,
                        mass_hi,
                        decoy: true,
                        ..fwd
                    };

                    [fwd, rev]
                })
            })
        })
        .collect::<Vec<_>>();

    // Essentially the same procedure for binning as the MS2 search engine
    // (see `database.rs` for explanation)
    ranges.par_sort_unstable_by(|a, b| a.rt.total_cmp(&b.rt));
    let min_rts = ranges
        .par_chunks_mut(16 * 1024)
        .map(|chunk| {
            // There should always be at least one item in the chunk!
            //  we know the chunk is already sorted by retention time too, so this is minimum value
            let min = chunk[0].rt;
            chunk.par_sort_unstable_by(|a, b| a.mass_lo.total_cmp(&b.mass_lo));
            min
        })
        .collect::<Vec<_>>();

    let mass_search_margin = ranges
        .iter()
        .map(|range| range.mass_hi - range.mass_lo)
        .fold(0.0_f32, f32::max);

    log::trace!("building feature map");
    FeatureMap {
        ranges,
        min_rts,
        bin_size: 16 * 1024,
        settings,
        mass_search_margin,
        ms2_confirmed,
        ms2_confirmed_strict,
        id_rts,
    }
}

struct Query<'a> {
    ranges: &'a [PrecursorRange],
    page_lo: usize,
    page_hi: usize,
    bin_size: usize,
    min_rt: f32,
    max_rt: f32,
    mass_search_margin: f32,
}

impl FeatureMap {
    fn rt_slice(&self, rt: f32, rt_tol: f32) -> Query<'_> {
        let (page_lo, page_hi) = binary_search_slice(
            &self.min_rts,
            |rt, x| rt.total_cmp(x),
            rt - rt_tol,
            rt + rt_tol,
        );

        Query {
            ranges: &self.ranges,
            page_lo,
            page_hi,
            bin_size: self.bin_size,
            max_rt: rt + rt_tol,
            min_rt: rt - rt_tol,
            mass_search_margin: self.mass_search_margin,
        }
    }
}

impl FeatureMap {
    /// Run label-free quantification module
    pub fn quantify(
        &self,
        db: &IndexedDatabase,
        spectra: &[ProcessedSpectrum],
        alignments: &[Alignment],
    ) -> HashMap<(PrecursorId, bool), QuantifiedPeak, fnv::FnvBuildHasher> {
        // Grids are keyed by the anchor file as well. With MBR every precursor
        // has a single cross-run anchor (`usize::MAX`); without MBR each
        // identified file has its own anchor RT, so it is traced on a private
        // single-file grid centered on that anchor.
        let scores: DashMap<(PrecursorId, bool, usize), Grid, fnv::FnvBuildHasher> =
            DashMap::default();
        let n_files = alignments.len();

        log::info!("tracing MS1 features");

        if spectra.is_empty() {
            log::warn!("no MS1 spectra found for quantification");
        } else {
            let rt_tol = self.settings.rt_tolerance();
            spectra.par_iter().for_each(|spectrum| {
                let rt = alignments[spectrum.file_id].transform(spectrum.scan_start_time);
                let query = self.rt_slice(rt, rt_tol);

                let add_entry = |entry: &PrecursorRange, intensity: f32| {
                    if !self.settings.mbr && entry.file_id != spectrum.file_id {
                        return;
                    }
                    let id = match self.settings.combine_charge_states {
                        true => PrecursorId::Combined(entry.peptide),
                        false => PrecursorId::Charged((entry.peptide, entry.charge)),
                    };

                    let (anchor, row) = match self.settings.mbr {
                        true => (usize::MAX, spectrum.file_id),
                        false => (entry.file_id, 0),
                    };

                    let mut grid = scores.entry((id, entry.decoy, anchor)).or_insert_with(|| {
                        let p = &db[entry.peptide];
                        let composition = p
                            .sequence
                            .iter()
                            .map(|r| composition(*r))
                            .sum::<Composition>();
                        let dist = crate::isotopes::peptide_isotopes(
                            composition.carbon,
                            composition.sulfur,
                        );
                        match self.settings.mbr {
                            true => Grid::new(entry, rt_tol, dist, n_files, GRID_SIZE),
                            false => {
                                let mut grid = Grid::new(entry, rt_tol, dist, 1, GRID_SIZE);
                                grid.reference_file_id = 0;
                                grid
                            }
                        }
                    });

                    grid.add_entry(rt, entry.isotope, row, intensity);
                };

                if spectrum.mobilities.is_empty() {
                    for (&mass, &intensity) in
                        spectrum.masses.iter().zip(spectrum.intensities.iter())
                    {
                        for entry in query.mass_lookup(mass) {
                            add_entry(entry, intensity);
                        }
                    }
                } else {
                    for ((&mass, &intensity), &mobility) in spectrum
                        .masses
                        .iter()
                        .zip(spectrum.intensities.iter())
                        .zip(spectrum.mobilities.iter())
                    {
                        for entry in query.mass_mobility_lookup(mass, mobility) {
                            add_entry(entry, intensity);
                        }
                    }
                }
            });
        }

        log::info!("integrating MS1 features");

        // Pair each target grid with its shifted decoy grid, so the decoy can
        // also be evaluated at the target's peak.
        let mut pairs: HashMap<(PrecursorId, usize), [Option<Grid>; 2], fnv::FnvBuildHasher> =
            HashMap::default();
        for ((id, decoy, anchor), grid) in scores {
            pairs.entry((id, anchor)).or_default()[usize::from(decoy)] = Some(grid);
        }

        let mut quantified = pairs
            .into_par_iter()
            .flat_map_iter(|((id, anchor), [target, decoy])| {
                // Identification RTs as bins of the target grid. The shifted
                // decoy grid has the same geometry, so it is seeded at the
                // same bins and runs through exactly the same apex search.
                let seeds = match (&target, self.settings.recenter_on_apex) {
                    (Some(grid), true) => Some(self.seeds(grid, id, anchor, n_files)),
                    _ => None,
                };
                let seeds = seeds.as_deref();
                let mut decoy_traces = decoy.map(|mut grid| grid.summarize_traces());
                let mut out = Vec::with_capacity(2);
                if let Some(mut grid) = target {
                    let mut traces = grid.summarize_traces();
                    if let Some((peak, areas, evidence, window)) =
                        traces.pick(&self.settings, seeds)
                    {
                        let paired = match &decoy_traces {
                            Some(decoy) => {
                                decoy
                                    .clone()
                                    .paired_evidence(&window, &self.settings, seeds)
                            }
                            None => vec![None; evidence.len()],
                        };
                        out.push(
                            self.expand(n_files, id, false, anchor, peak, areas, evidence, paired),
                        );
                    }
                }
                if let Some(traces) = decoy_traces.as_mut() {
                    if let Some((peak, areas, evidence, _)) = traces.pick(&self.settings, seeds) {
                        out.push(self.expand(
                            n_files,
                            id,
                            true,
                            anchor,
                            peak,
                            areas,
                            evidence,
                            Vec::new(),
                        ));
                    }
                }
                out
            })
            .collect::<Vec<_>>();

        // Merge per-file peaks in a fixed (precursor, anchor file) order so the
        // reported peak does not depend on thread scheduling.
        quantified.par_sort_unstable_by_key(|(key, anchor, _)| (*key, *anchor));
        let mut peaks: HashMap<_, QuantifiedPeak, fnv::FnvBuildHasher> = HashMap::default();
        for (key, _, quantified) in quantified {
            match peaks.entry(key) {
                std::collections::hash_map::Entry::Occupied(mut merged) => {
                    merged.get_mut().merge(quantified)
                }
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(quantified);
                }
            }
        }
        // Report peak positions in each file's own RT units.
        peaks.par_iter_mut().for_each(|((id, decoy), quantified)| {
            for (file, evidence) in quantified.file_evidence.iter_mut().enumerate() {
                let Some(evidence) = evidence else { continue };
                let alignment = &alignments[file];
                let apex = alignment.inverse(evidence.apex_rt);
                evidence.fwhm = evidence.fwhm.map(|width| {
                    let half = width / 2.0;
                    alignment.inverse(evidence.apex_rt + half)
                        - alignment.inverse(evidence.apex_rt - half)
                });
                evidence.apex_rt = apex;
                evidence.peak_start_rt = alignment.inverse(evidence.peak_start_rt);
                evidence.peak_end_rt = alignment.inverse(evidence.peak_end_rt);
                // Decoys have no identification of their own.
                evidence.id_apex_offset = (!*decoy)
                    .then(|| self.id_rts.get(&(*id, file)))
                    .flatten()
                    .map(|rt| rt.observed - apex);
            }
        });
        peaks
    }

    /// Bin of the best identification in each grid row, if it is inside the
    /// grid. Row `r` is file `r` with MBR, or the anchor file without.
    fn seeds(
        &self,
        grid: &Grid,
        id: PrecursorId,
        anchor: usize,
        n_files: usize,
    ) -> Vec<Option<usize>> {
        let files = match anchor {
            usize::MAX => (0..n_files).collect::<Vec<_>>(),
            anchor => vec![anchor],
        };
        files
            .into_iter()
            .map(|file| {
                let rt = self.id_rts.get(&(id, file))?.aligned;
                let bin = ((rt - grid.rt_min) / grid.rt_step).round();
                (bin >= 0.0 && (bin as usize) < grid.matrix.cols).then_some(bin as usize)
            })
            .collect()
    }
}

impl FeatureMap {
    /// Expand a single-file (MBR disabled) grid back to all files and attach
    /// the per-file labels.
    #[allow(clippy::too_many_arguments)]
    fn expand(
        &self,
        n_files: usize,
        id: PrecursorId,
        decoy: bool,
        anchor: usize,
        peak: Peak,
        areas: Vec<Option<f64>>,
        evidence: Vec<Option<FileEvidence>>,
        paired: Vec<Option<FileEvidence>>,
    ) -> ((PrecursorId, bool), usize, QuantifiedPeak) {
        let (intensities, mut file_evidence, mut paired_decoy_evidence) = match anchor {
            usize::MAX => (areas, evidence, paired),
            anchor => {
                let mut intensities = vec![None; n_files];
                let mut file_evidence = vec![None; n_files];
                intensities[anchor] = areas[0];
                file_evidence[anchor] = evidence.into_iter().next().flatten();
                let mut paired_decoy_evidence = Vec::new();
                if !decoy {
                    paired_decoy_evidence = vec![None; n_files];
                    paired_decoy_evidence[anchor] = paired.into_iter().next().flatten();
                }
                (intensities, file_evidence, paired_decoy_evidence)
            }
        };
        for (file_id, evidence) in file_evidence
            .iter_mut()
            .chain(paired_decoy_evidence.iter_mut())
            .enumerate()
        {
            if let Some(evidence) = evidence {
                evidence.transfer_candidate = self.settings.mbr
                    && !self.ms2_confirmed_strict.contains(&(id, file_id % n_files));
            }
        }
        let ms2_confirmed = (0..n_files)
            .map(|file_id| !decoy && self.ms2_confirmed.contains(&(id, file_id)))
            .collect();
        let ms2_confirmed_strict = (0..n_files)
            .map(|file_id| !decoy && self.ms2_confirmed_strict.contains(&(id, file_id)))
            .collect();
        (
            (id, decoy),
            anchor,
            QuantifiedPeak {
                peak,
                intensities,
                ms2_confirmed,
                ms2_confirmed_strict,
                file_evidence,
                paired_decoy_evidence,
            },
        )
    }
}

impl QuantifiedPeak {
    /// Combine the peak of another anchor file into this one: per-file signals
    /// are taken from whichever anchor observed them, and the highest-scoring
    /// peak (earliest anchor on ties) represents the precursor.
    fn merge(&mut self, other: QuantifiedPeak) {
        let mut paired = other.paired_decoy_evidence.into_iter();
        for (file, (area, evidence)) in other
            .intensities
            .into_iter()
            .zip(other.file_evidence)
            .enumerate()
        {
            let decoy = paired.next().flatten();
            if area.is_some() {
                self.intensities[file] = area;
                self.file_evidence[file] = evidence;
                if let Some(slot) = self.paired_decoy_evidence.get_mut(file) {
                    *slot = decoy;
                }
            }
        }
        if other.peak.score > self.peak.score {
            self.peak = other.peak;
        }
    }
}

pub struct Grid {
    rt_min: f32,
    rt_step: f32,
    files: usize,
    /// File with the most confident PSM
    reference_file_id: usize,
    /// Relative theoretical isotopic abundances
    pub distribution: [f32; N_ISOTOPES],
    /// Matrix of summed intensities for each isotopic trace in each file, divided
    /// among equally spaced retention time bins. This is a [N_FILES * N_ISOTOPES, GRID_SIZE]
    /// sized matrix.
    ///
    /// Isotopic summed intensities are arranged in consequtive rows, ordered by file. The first
    /// N_ISOTOPE rows correspond to file 0, then the following N_ISOTOPE rows correspond to file 1, etc.
    /// Indexing into the matrix is done by [(n_file * N_ISOTOPES + isotope, rt_window)]
    pub matrix: Matrix,
}

#[derive(Clone)]
pub struct Traces {
    /// Matrix of dot(MS1 ions, Grid.distribution). This collapses our N_FILES * N_ISOTOPES rows
    /// down to just N_FILES
    pub dot_product: Matrix,
    /// Matrix of spectral angles at each retention time for each file
    pub spectral_angle: Matrix,
    /// File with the most confident PSM
    reference_file_id: usize,
    /// Consensus RT of bin 0 and the bin width.
    pub geometry: GridGeometry,
}

/// Position of a trace grid on the consensus RT axis.
#[derive(Copy, Clone, Debug, Default)]
pub struct GridGeometry {
    pub rt_min: f32,
    pub rt_step: f32,
}

impl GridGeometry {
    fn rt(&self, bin: f64) -> f32 {
        self.rt_min + (bin as f32) * self.rt_step
    }
}

#[derive(Clone, Debug, Default)]
pub struct Peak {
    /// Discretized retention time
    pub rt: usize,
    /// Intensity weighted normalized spectral angle
    pub spectral_angle: f64,
    /// Peak score
    pub score: f64,

    pub q_value: f32,
}

impl Traces {
    /// Calculate and apply time warping factors
    fn warp(&mut self) -> Vec<isize> {
        let time_warps = self.find_time_warps(&self.dot_product, WARP_SLACK);
        Self::apply_time_warps(&mut self.spectral_angle, &time_warps);
        Self::apply_time_warps(&mut self.dot_product, &time_warps);
        time_warps
    }

    /// Choose, align and integrate the peak: around the identification RTs'
    /// apexes when `seeds` are given (apex recentering), otherwise the
    /// best-scoring bin of the traced window.
    #[allow(clippy::type_complexity)]
    pub fn pick(
        &mut self,
        settings: &LfqSettings,
        seeds: Option<&[Option<usize>]>,
    ) -> Option<(
        Peak,
        Vec<Option<f64>>,
        Vec<Option<FileEvidence>>,
        PeakWindow,
    )> {
        // Precursors without an identification inside the grid (e.g. label
        // channels seeded by another channel) keep the window search.
        match seeds.filter(|seeds| seeds.iter().any(Option::is_some)) {
            Some(seeds) => self.integrate_apex(settings, seeds),
            None => self.integrate_window(settings),
        }
    }

    /// Isotope-consistent trace: signal where the isotope pattern matches.
    fn consistent(&self, file: usize, settings: &LfqSettings) -> Vec<f64> {
        self.dot_product
            .row_slice(file)
            .iter()
            .zip(self.spectral_angle.row_slice(file))
            .map(|(&dot, &angle)| {
                if angle >= settings.spectral_angle {
                    dot
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Integrate the peak whose apex is nearest the identification RTs.
    ///
    /// * Each seeded (identified) row climbs from its identification bin to
    ///   the apex of the isotope-consistent trace it sits on.
    /// * The consensus apex is the median of those apexes; each seeded row is
    ///   shifted so that its own apex lands on it.
    /// * Unseeded (transfer) rows get the usual warp search towards the sum of
    ///   the aligned seeded traces.
    /// * Integration bounds follow the summed seeded trace down to a valley
    ///   or [`APEX_BOUND_FRACTION`] of its apex, without fixed caps.
    #[allow(clippy::type_complexity)]
    pub fn integrate_apex(
        &mut self,
        settings: &LfqSettings,
        seeds: &[Option<usize>],
    ) -> Option<(
        Peak,
        Vec<Option<f64>>,
        Vec<Option<FileEvidence>>,
        PeakWindow,
    )> {
        let rows = self.dot_product.rows;
        let cols = self.dot_product.cols;
        let apexes = (0..rows)
            .map(|file| {
                let seed = (*seeds.get(file)?)?;
                climb_to_apex(&self.consistent(file, settings), seed)
            })
            .collect::<Vec<_>>();
        let mut found = apexes.iter().flatten().copied().collect::<Vec<_>>();
        if found.is_empty() {
            return None;
        }
        found.sort_unstable();
        let apex = found[(found.len() - 1) / 2];

        let mut shifts = vec![0isize; rows];
        let mut reference = vec![0.0; cols];
        let mut bounds_trace = vec![0.0; cols];
        for (file, file_apex) in apexes.iter().enumerate() {
            if let Some(file_apex) = file_apex {
                shifts[file] = *file_apex as isize - apex as isize;
                let consistent = self.consistent(file, settings);
                for (i, (r, b)) in reference.iter_mut().zip(&mut bounds_trace).enumerate() {
                    let j = i as isize + shifts[file];
                    if j >= 0 && (j as usize) < cols {
                        *r += self.dot_product[(file, j as usize)];
                        *b += consistent[j as usize];
                    }
                }
            }
        }
        for (file, file_apex) in apexes.iter().enumerate() {
            if file_apex.is_none() {
                shifts[file] =
                    Self::best_shift(&reference, self.dot_product.row_slice(file), WARP_SLACK);
            }
        }
        Self::apply_time_warps(&mut self.spectral_angle, &shifts);
        Self::apply_time_warps(&mut self.dot_product, &shifts);

        let (left, right) = peak_bounds(&bounds_trace, apex, APEX_BOUND_FRACTION);

        // Same form as the hybrid score, with retention-time proximity
        // measured from the identification RTs rather than the window center.
        let mut max = 0.0f64;
        let mut at_apex = (0.0, 1.0);
        for col in 0..cols {
            let mut summed_int = 1.0;
            let mut weighted = 0.0;
            for (sa, dotp) in self.spectral_angle.col(col).zip(self.dot_product.col(col)) {
                weighted += sa * dotp;
                summed_int += dotp;
            }
            max = max.max(summed_int);
            if col == apex {
                at_apex = (weighted / summed_int, summed_int);
            }
        }
        let mut seeded = seeds.iter().flatten().copied().collect::<Vec<_>>();
        seeded.sort_unstable();
        let center = seeded[(seeded.len() - 1) / 2];
        let half = (cols / 2).max(1) as f64;
        let proximity = (1.0 - (apex as f64 - center as f64).abs() / half).max(0.0);
        let (spectral_angle, intensity) = at_apex;
        let score =
            spectral_angle.max(0.0).powi(3) * proximity.powf(0.33) * (intensity / max).sqrt();
        if score <= 0.0 || !score.is_finite() {
            return None;
        }
        let peak = Peak {
            rt: apex,
            spectral_angle,
            score,
            q_value: 0.0,
        };
        let window = PeakWindow {
            rt: apex,
            left,
            right,
            shifts,
            reference,
        };
        let (areas, evidence) = self.window_evidence(&window, settings);
        Some((peak, areas, evidence, window))
    }

    /// Find time warping offsets for each file that maximize the dot product
    /// with the most intense run
    ///
    /// * Use the LC-MS run with the most confident PSM for a peptide as the reference run
    /// * For each LC-MS run, find the time warping shif that maximizes dot product
    ///   with the `reference` run
    pub fn find_time_warps(&self, matrix: &Matrix, slack: isize) -> Vec<isize> {
        Self::time_warps_to(matrix.row_slice(self.reference_file_id), matrix, slack)
    }

    /// Time warping offsets for each row of `matrix` that maximize the dot
    /// product with `reference`
    fn time_warps_to(reference: &[f64], matrix: &Matrix, slack: isize) -> Vec<isize> {
        (0..matrix.rows)
            .map(|row| Self::best_shift(reference, matrix.row_slice(row), slack))
            .collect()
    }

    /// Offset of `run` that maximizes its dot product with `reference`.
    fn best_shift(reference: &[f64], run: &[f64], slack: isize) -> isize {
        let mut best_offset = (0, 0.0);
        for offset in -slack..=slack {
            let mut dot = 0.0;
            for (i, ref_int) in reference.iter().enumerate() {
                let j = i as isize + offset;
                if j >= 0 && j < run.len() as isize {
                    dot += ref_int * run[j as usize];
                }
            }

            if dot >= best_offset.1 {
                best_offset = (offset, dot);
            }
        }
        best_offset.0
    }

    /// Perform local Correlation Optimization Warping
    fn apply_time_warps(matrix: &mut Matrix, time_warps: &[isize]) {
        for (row, warp) in time_warps.iter().enumerate() {
            let run = matrix.row_slice_mut(row);
            let mut shifted = vec![0.0; run.len()];
            for (i, val) in shifted.iter_mut().enumerate() {
                let j = i as isize + warp;
                if j >= 0 && j < run.len() as isize {
                    *val = run[j as usize];
                }
            }
            run.copy_from_slice(&shifted);
        }
    }

    pub fn scores(&self, strategy: PeakScoringStrategy) -> (Vec<f64>, Vec<f64>) {
        let mut spectral = Vec::with_capacity(self.spectral_angle.cols);
        let mut intensity = Vec::with_capacity(self.spectral_angle.cols);
        let mut max = 0.0f64;
        for col in 0..self.spectral_angle.cols {
            let mut summed_int = 1.0;
            let mut weighted = 0.0;
            for (sa, dotp) in self.spectral_angle.col(col).zip(self.dot_product.col(col)) {
                weighted += sa * dotp;
                summed_int += dotp;
            }
            spectral.push(weighted / summed_int);
            intensity.push(summed_int);
            max = max.max(summed_int);
        }

        let center = self.spectral_angle.cols as isize / 2;
        let scores = spectral
            .iter()
            .zip(intensity.iter())
            .enumerate()
            .map(|(rt, (s, i))| match strategy {
                PeakScoringStrategy::RetentionTime => {
                    (1.0 - ((rt as isize - center).abs() as f64 / center as f64)).powf(0.33)
                }
                PeakScoringStrategy::SpectralAngle => *s,
                PeakScoringStrategy::Intensity => (*i / max).sqrt(),
                PeakScoringStrategy::Hybrid => {
                    let rt = 1.0 - ((rt as isize - center).abs() as f64 / center as f64);
                    s.powi(3) * rt.powf(0.33) * (*i / max).sqrt()
                    // s.powi(3) * (*i / max).sqrt()
                }
            })
            .collect();
        (scores, spectral)
    }

    /// Align and integrate MS1 traces across files
    ///
    /// * Calculate time warping factors for each file, aligning them so that
    ///   the correlation between them is maximized (we just do a dot product)
    /// * Locate the retention time window corresponding to the maximum average
    ///   angle observed across all of the files
    /// * Integrate all of the MS1 traces within said window, returning a vector
    ///   of length `n_files` containing the summed MS1 intensities
    pub fn integrate(&mut self, settings: &LfqSettings) -> Option<(Peak, Vec<Option<f64>>)> {
        self.integrate_with_evidence(settings)
            .map(|(peak, areas, _)| (peak, areas))
    }

    #[allow(clippy::type_complexity)]
    pub fn integrate_with_evidence(
        &mut self,
        settings: &LfqSettings,
    ) -> Option<(Peak, Vec<Option<f64>>, Vec<Option<FileEvidence>>)> {
        self.integrate_window(settings)
            .map(|(peak, areas, evidence, _)| (peak, areas, evidence))
    }

    /// Integrate the best peak and also return its window, so that a paired
    /// decoy can be evaluated at exactly the same place (see [`Traces::paired_evidence`]).
    #[allow(clippy::type_complexity)]
    pub fn integrate_window(
        &mut self,
        settings: &LfqSettings,
    ) -> Option<(
        Peak,
        Vec<Option<f64>>,
        Vec<Option<FileEvidence>>,
        PeakWindow,
    )> {
        let shifts = self.warp();

        let (scores, spectral) = self.scores(settings.peak_scoring);
        let mut best = Peak::default();
        for (rt, s) in scores.iter().enumerate() {
            if *s > best.score && spectral[rt] >= settings.spectral_angle {
                best.score = *s;
                best.rt = rt;
            }
        }

        if best.score == 0.0 {
            return None;
        }

        // Find peak boundaries
        let mut left = best.rt.saturating_sub(1);
        let mut right = best.rt.saturating_add(1);

        let threshold = best.score * 0.50;

        // Don't let peaks extend more than GRID_SIZE/5 bins to either side
        while left > best.rt.saturating_sub(scores.len() / 5)
            && scores[left] >= threshold
            && spectral[left] >= settings.spectral_angle
        {
            left -= 1;
        }

        while right < scores.len().saturating_sub(1).min(best.rt + 20)
            && scores[right] >= threshold
            && spectral[right] >= settings.spectral_angle
        {
            right += 1;
        }

        let mut summed_int = 1.0;
        let mut weighted = 0.0;
        for (sa, dotp) in self
            .spectral_angle
            .col(best.rt)
            .zip(self.dot_product.col(best.rt))
        {
            weighted += sa * dotp;
            summed_int += dotp;
        }
        best.spectral_angle = weighted / summed_int;
        let window = PeakWindow {
            rt: best.rt,
            left,
            right,
            shifts,
            reference: self.dot_product.row_slice(self.reference_file_id).to_vec(),
        };
        let (areas, evidence) = self.window_evidence(&window, settings);
        Some((best, areas, evidence, window))
    }

    /// Evaluate these (decoy) traces at a peak chosen from other traces: the
    /// same apex, integration window and reference trace. Each file gets the
    /// same local warp search the target files got (towards that reference),
    /// so a file row answers the same question as the target row it is paired
    /// with: can an isotope envelope be aligned to the expected elution?
    ///
    /// With apex recentering (`seeds`), a seeded row instead climbs from the
    /// same identification bin as the target row and is shifted so that its
    /// own apex lands on the target's apex, exactly as the target row was. A
    /// row whose climb finds no signal keeps the warp search, as a target row
    /// does in [`Traces::integrate_apex`].
    pub fn paired_evidence(
        &mut self,
        window: &PeakWindow,
        settings: &LfqSettings,
        seeds: Option<&[Option<usize>]>,
    ) -> Vec<Option<FileEvidence>> {
        let mut shifts = Self::time_warps_to(&window.reference, &self.dot_product, WARP_SLACK);
        for (file, seed) in seeds.unwrap_or_default().iter().enumerate() {
            if let Some(apex) =
                seed.and_then(|seed| climb_to_apex(&self.consistent(file, settings), seed))
            {
                shifts[file] = apex as isize - window.rt as isize;
            }
        }
        Self::apply_time_warps(&mut self.spectral_angle, &shifts);
        Self::apply_time_warps(&mut self.dot_product, &shifts);
        let window = PeakWindow {
            shifts,
            reference: window.reference.clone(),
            ..*window
        };
        self.window_evidence(&window, settings).1
    }

    fn window_evidence(
        &self,
        window: &PeakWindow,
        settings: &LfqSettings,
    ) -> (Vec<Option<f64>>, Vec<Option<FileEvidence>>) {
        let mut areas = Vec::with_capacity(self.dot_product.rows);
        for file in 0..self.dot_product.rows {
            let area = match settings.integration {
                IntegrationStrategy::Sum => self.dot_product.row_slice(file)
                    [window.left..window.right]
                    .iter()
                    .sum::<f64>(),
                IntegrationStrategy::Apex => self.dot_product.row_slice(file)[window.rt],
            };

            areas.push((area.is_finite() && area > 0.0).then_some(area));
        }

        let reference = &window.reference;
        let reference_norm = reference.iter().map(|x| x * x).sum::<f64>().sqrt();
        let evidence = areas
            .iter()
            .enumerate()
            .map(|(file, area)| {
                area.map(|_| {
                    let trace = self.dot_product.row_slice(file);
                    let norm = trace.iter().map(|x| x * x).sum::<f64>().sqrt();
                    let dot = trace.iter().zip(reference).map(|(x, y)| x * y).sum::<f64>();
                    let trace_cosine = if norm > 0.0 && reference_norm > 0.0 {
                        (dot / (norm * reference_norm)).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    let spectral_angle = self.spectral_angle[(file, window.rt)].clamp(0.0, 1.0);
                    let proximity = (1.0
                        - window.shifts[file].unsigned_abs() as f64 / self.dot_product.cols as f64)
                        .max(0.0);
                    // Positions in the unwarped grid of this file, as consensus RT.
                    let shift = window.shifts[file] as f64;
                    let fwhm = half_maximum_width(trace, window.rt)
                        .map(|bins| bins as f32 * self.geometry.rt_step);
                    FileEvidence {
                        transfer_candidate: false,
                        spectral_angle,
                        trace_cosine,
                        rt_shift_bins: window.shifts[file] as i32,
                        score: spectral_angle.powi(3) * trace_cosine * proximity,
                        extraction_q_value: None,
                        apex_rt: self.geometry.rt(window.rt as f64 + shift),
                        peak_start_rt: self.geometry.rt(window.left as f64 + shift),
                        peak_end_rt: self.geometry.rt(window.right as f64 - 1.0 + shift),
                        fwhm,
                        id_apex_offset: None,
                    }
                })
            })
            .collect();
        (areas, evidence)
    }
}

/// Where a precursor's cross-run peak was found: local warps, apex and
/// integration bounds, and the warped reference trace.
#[derive(Clone, Debug)]
pub struct PeakWindow {
    pub rt: usize,
    pub left: usize,
    pub right: usize,
    pub shifts: Vec<isize>,
    pub reference: Vec<f64>,
}

impl Grid {
    pub fn new(
        entry: &PrecursorRange,
        rt_tol: f32,
        distribution: [f32; N_ISOTOPES],
        files: usize,
        grid_size: usize,
    ) -> Grid {
        let matrix = Matrix::new(
            vec![0.0; grid_size * files * N_ISOTOPES],
            files * N_ISOTOPES,
            grid_size,
        );
        let rt_step = (rt_tol * 2.0) / (grid_size) as f32;

        Grid {
            rt_min: entry.rt - rt_tol,
            rt_step,
            distribution,
            matrix,
            files,
            reference_file_id: entry.file_id,
        }
    }

    /// Add a data point to the integration grid
    pub fn add_entry(&mut self, spectrum_rt: f32, isotope: usize, file_id: usize, intensity: f32) {
        let bin_lo = ((spectrum_rt - self.rt_min) / self.rt_step).floor() as usize;
        let bin_lo = bin_lo.min(self.matrix.cols - 1);
        let bin_hi = (bin_lo + 1).min(self.matrix.cols - 1);

        let bin_lo_rt = bin_lo as f32 * self.rt_step + self.rt_min;
        // what fraction [0.0, 1.0] of the way are we to the higher bin?
        let interp = (spectrum_rt - bin_lo_rt) / self.rt_step;

        self.matrix[(file_id * N_ISOTOPES + isotope, bin_lo)] +=
            ((1.0 - interp) * intensity) as f64;
        self.matrix[(file_id * N_ISOTOPES + isotope, bin_hi)] += (interp * intensity) as f64;
    }

    /// Combine individual isotopic traces across files into aligned, summed
    /// MS1 traces for each file
    ///
    /// * Perform gaussian smoothing on summed intensities
    /// * Calculate normalized spectral angle for observed isotopic distribution
    ///   relative to theoretical distribution
    pub fn summarize_traces(&mut self) -> Traces {
        let k = gaussian_kernel(0.5, K_WIDTH);

        let mut spectral_angle = Matrix::new(
            vec![0.0; self.files * self.matrix.cols],
            self.files,
            self.matrix.cols,
        );

        let mut dot_product = spectral_angle.clone();

        // square root of the summed squared relative abundances of theoretical
        // isotopic distribution
        let ss_dist = self
            .distribution
            .iter()
            .map(|x| x.powi(2))
            .sum::<f32>()
            .sqrt() as f64;

        for file in 0..self.files {
            let mut summed_squared_intensities = vec![0.0; self.matrix.cols];
            for isotope in 0..N_ISOTOPES {
                let convolved = convolve(self.matrix.row_slice(file * N_ISOTOPES + isotope), &k);
                for (col, intensity) in convolved.iter().enumerate() {
                    spectral_angle[(file, col)] += intensity * self.distribution[isotope] as f64;
                    summed_squared_intensities[col] += intensity.powi(2);
                }
                self.matrix
                    .row_slice_mut(file * N_ISOTOPES + isotope)
                    .copy_from_slice(&convolved);
            }

            for (col, ss) in summed_squared_intensities.iter().enumerate() {
                let dot = spectral_angle[(file, col)];
                let similarity = if *ss > 0.0 {
                    dot / (ss.sqrt() * ss_dist)
                } else {
                    0.0
                };

                // Calculate the normalized spectral angle
                spectral_angle[(file, col)] = 1.0 - 2.0 * similarity.acos() / std::f64::consts::PI;
                dot_product[(file, col)] = dot;
            }
        }

        Traces {
            dot_product,
            spectral_angle,
            reference_file_id: self.reference_file_id,
            geometry: GridGeometry {
                rt_min: self.rt_min,
                rt_step: self.rt_step,
            },
        }
    }
}

/// Fraction of the summed apex height at which apex-mode integration stops,
/// unless a valley comes first. Half height matches the window search; 20%
/// took in more interference and worsened the E. coli ratio on PXD028735.
const APEX_BOUND_FRACTION: f64 = 0.5;

/// Climb from `seed` to the apex of the peak it sits on. The seed snaps to
/// the nearest signal within [`APEX_SNAP_BINS`]; the climb then moves to the
/// highest point reachable without the trace falling below half of the
/// current apex, so that small dips in a smoothed peak do not stop it early.
fn climb_to_apex(trace: &[f64], seed: usize) -> Option<usize> {
    let n = trace.len();
    let mut apex = (0..=APEX_SNAP_BINS)
        .flat_map(|d| [seed.checked_sub(d), Some(seed + d).filter(|&i| i < n)])
        .flatten()
        .find(|&i| trace[i] > 0.0)?;
    loop {
        let height = trace[apex];
        let mut next = None;
        for dir in [-1isize, 1] {
            let mut i = apex as isize + dir;
            while i >= 0 && (i as usize) < n && trace[i as usize] >= 0.5 * height {
                if trace[i as usize] > trace[next.unwrap_or(apex)] {
                    next = Some(i as usize);
                }
                i += dir;
            }
        }
        match next {
            Some(higher) => apex = higher,
            None => return Some(apex),
        }
    }
}

/// Integration bounds `[left, right)` around `apex`: extend while the trace
/// stays above `fraction` of the apex height and keeps falling (a valley ends
/// the peak).
fn peak_bounds(trace: &[f64], apex: usize, fraction: f64) -> (usize, usize) {
    let floor = trace[apex] * fraction;
    let mut left = apex;
    while left > 0 && trace[left - 1] >= floor && trace[left - 1] <= trace[left] {
        left -= 1;
    }
    let mut right = apex;
    while right + 1 < trace.len() && trace[right + 1] >= floor && trace[right + 1] <= trace[right] {
        right += 1;
    }
    (left, right + 1)
}

/// Width, in (fractional) bins, of `trace` at half its height at `apex`.
fn half_maximum_width(trace: &[f64], apex: usize) -> Option<f64> {
    let half = trace[apex] / 2.0;
    if half <= 0.0 {
        return None;
    }
    let crossing = |mut i: usize, step: isize| -> Option<f64> {
        loop {
            let j = i as isize + step;
            if j < 0 || j as usize >= trace.len() {
                return None;
            }
            let j = j as usize;
            if trace[j] < half {
                // Interpolate between i (above) and j (below).
                let frac = (trace[i] - half) / (trace[i] - trace[j]);
                return Some(i as f64 + step as f64 * frac);
            }
            i = j;
        }
    };
    Some(crossing(apex, 1)? - crossing(apex, -1)?)
}

/// Create a symmetrical gaussian kernel of given standard deviation and length
fn gaussian_kernel(sigma: f64, len: usize) -> Vec<f64> {
    let step = 2.0 / (len - 1) as f64;
    let constant = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());

    let mut kernel = (0..len)
        .map(|i| {
            let x = i as f64 * step - 1.0;
            constant * (-0.5 * (x / sigma).powi(2)).exp()
        })
        .collect::<Vec<_>>();

    let sum = kernel.iter().sum::<f64>();
    kernel.iter_mut().for_each(|x| *x /= sum);
    kernel
}

/// Convolve a signal with a symmetrical kernel
/// - This should behave the same as `np.convolve(..., mode='same')`
fn convolve(slice: &[f64], kernel: &[f64]) -> Vec<f64> {
    // Middle index of kernel
    let n = kernel.len() - (kernel.len() / 2);

    (0..slice.len())
        .map(|idx| {
            // If idx < kernel.len(), take only some subset of the kernel (at least half)
            let k = &kernel[kernel.len().saturating_sub(n + idx)..];
            // If idx < kernel.len(), then start from 0
            let w = &slice[idx.saturating_sub(n - 1)..];
            // Dot product
            w.iter().zip(k).fold(0.0, |acc, (x, y)| acc + x * y)
        })
        .collect()
}

impl Query<'_> {
    pub fn mass_lookup(&self, mass: f32) -> impl Iterator<Item = &PrecursorRange> {
        (self.page_lo..self.page_hi).flat_map(move |page| {
            let left_idx = page * self.bin_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = (left_idx + self.bin_size).min(self.ranges.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &self.ranges[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |frag, bounds| frag.mass_lo.total_cmp(bounds),
                mass - self.mass_search_margin,
                mass,
            );

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                frag.rt <= self.max_rt
                    && frag.rt >= self.min_rt
                    && mass >= frag.mass_lo
                    && mass <= frag.mass_hi
            })
        })
    }

    pub fn mass_mobility_lookup(
        &self,
        mass: f32,
        mobility: f32,
    ) -> impl Iterator<Item = &PrecursorRange> {
        self.mass_lookup(mass).filter(move |precursor| {
            (precursor.mobility_hi >= mobility) && (precursor.mobility_lo <= mobility)
        })
    }
}

#[cfg(test)]
#[path = "../tests/unit/lfq.rs"]
mod tests;
