//! Label-free quantification: MS1 feature tracing, integration and the reported areas.

use crate::database::{binary_search_slice, IndexedDatabase, PeptideIx};
use crate::mass::{composition, Composition, Tolerance, NEUTRON};
use crate::ml::retention_alignment::Alignment;
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

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
        .filter(|(_, peptide)| !peptide.decoy && peptide.label_channel().is_some())
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
                .label_channel()
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

mod grid;
mod signal;

pub use grid::*;
use signal::*;

#[cfg(test)]
#[path = "../tests/unit/lfq.rs"]
mod tests;
