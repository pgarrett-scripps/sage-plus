//! Spectrum-indexed exact prefilter.
//!
//! The classic prefilter builds a theoretical fragment index for every
//! database chunk and streams all spectra through it. This module inverts the
//! loop: the spectra are indexed once, and generated peptides are streamed
//! through the spectrum index without ever building or sorting a fragment
//! index.
//!
//! A peptide is kept when, for at least one precursor hypothesis (charge,
//! isotope error, and mass offset) of a spectrum whose window contains the
//! peptide mass, at least `min_matched_peaks` (peak, preliminary fragment)
//! pairs fall inside the fragment tolerance. This is the preliminary match
//! count of the search. With `min_matched_peaks = 1` and no peak cap, the
//! retained set is identical to [`crate::scoring::Scorer::exact_prefilter`].
//! Every window is computed with the same floating-point operations as the
//! fragment index query, so the comparison is exact rather than approximate.

use crate::database::{preliminary_fragment_masses, MassOffset, Parameters};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::peptide::Peptide;
use crate::scoring::{max_fragment_charge, offset_query, AtomicBitSet, FragmentMatchIndex};
use crate::spectrum::ProcessedSpectrum;
use rayon::prelude::*;

const NO_PROBE: u32 = u32::MAX;

/// Build the global peak index only when per-probe lookups are expected to
/// cost this many times more than global scans (open and wide-window searches).
const GLOBAL_INDEX_MARGIN: f64 = 2.0;

/// Search settings that determine the precursor hypotheses of a spectrum.
/// These mirror the corresponding [`crate::scoring::Scorer`] fields.
#[derive(Clone, Debug)]
pub struct SpectrumIndexSettings {
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub min_isotope_err: i8,
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub wide_window: bool,
    pub min_peaks: usize,
    /// Preliminary matches a precursor hypothesis needs to keep a peptide.
    pub min_matched_peaks: u16,
    /// Only the most intense peaks of each spectrum are indexed.
    pub max_peaks: Option<usize>,
}

/// Accumulates spectra, possibly across several file batches, before the
/// index is sorted.
pub struct SpectrumIndexBuilder {
    settings: SpectrumIndexSettings,
    /// Precursor mass delta and fragment shift of each configured offset.
    offsets: Vec<(f32, f32)>,
    windows: Vec<PrecursorWindow>,
    /// Neutral fragment masses of every peak set, sorted within each set.
    peaks: Vec<f32>,
    peakset_starts: Vec<u64>,
    probes: Vec<Probe>,
    spectra: usize,
}

#[derive(Copy, Clone)]
struct PrecursorWindow {
    lo: f32,
    hi: f32,
    /// Unshifted probe and, for mass offsets, the shifted probe of the same
    /// peaks.
    probes: [u32; 2],
}

/// Peaks of one spectrum for one fragment-charge limit, searched with one
/// fragment shift. Unshifted probes use a shift of zero, which leaves every
/// bound unchanged.
#[derive(Copy, Clone)]
struct Probe {
    peakset: u32,
    shift: f32,
}

#[derive(Default)]
struct LocalSpectrum {
    peaksets: Vec<Vec<f32>>,
    probes: Vec<Probe>,
    windows: Vec<PrecursorWindow>,
}

/// Fragment window of a peak, computed exactly as the fragment index query
/// does for unshifted and shifted lookups.
#[inline(always)]
fn fragment_window(tolerance: Tolerance, mass: f32, shift: f32) -> (f32, f32) {
    let (lo, hi) = tolerance.bounds(mass);
    (lo - shift, hi - shift)
}

impl SpectrumIndexBuilder {
    pub fn new(settings: SpectrumIndexSettings, mass_offsets: &[MassOffset]) -> Self {
        Self {
            settings,
            offsets: mass_offsets
                .iter()
                .map(|offset| (offset.mass(), offset.fragment_shift()))
                .collect(),
            windows: Vec::new(),
            peaks: Vec::new(),
            peakset_starts: vec![0],
            probes: Vec::new(),
            spectra: 0,
        }
    }

    pub fn spectra(&self) -> usize {
        self.spectra
    }

    /// Approximate size of the finished index.
    pub fn allocated_bytes(&self) -> usize {
        self.peaks.capacity() * 4
            + self.peakset_starts.capacity() * 8
            + self.probes.capacity() * std::mem::size_of::<Probe>()
            + self.windows.capacity() * (std::mem::size_of::<PrecursorWindow>() + 4)
    }

    /// Add a batch of spectra. The batch can be dropped afterwards.
    pub fn add(&mut self, spectra: &[ProcessedSpectrum]) {
        let local = spectra
            .par_iter()
            .filter(|spectrum| spectrum.is_searchable(self.settings.min_peaks))
            .map(|spectrum| self.local_spectrum(spectrum))
            .collect::<Vec<_>>();

        for spectrum in local {
            self.spectra += 1;
            let peakset_base = (self.peakset_starts.len() - 1) as u32;
            let probe_base = self.probes.len() as u32;
            for peaks in spectrum.peaksets {
                self.peaks.extend_from_slice(&peaks);
                self.peakset_starts.push(self.peaks.len() as u64);
            }
            self.probes
                .extend(spectrum.probes.into_iter().map(|probe| Probe {
                    peakset: probe.peakset + peakset_base,
                    ..probe
                }));
            self.windows
                .extend(spectrum.windows.into_iter().map(|mut window| {
                    for probe in &mut window.probes {
                        if *probe != NO_PROBE {
                            *probe += probe_base;
                        }
                    }
                    window
                }));
        }
        assert!(
            self.probes.len() < NO_PROBE as usize,
            "spectrum index exceeds 32-bit probe identifiers"
        );
    }

    /// Enumerate the precursor hypotheses exactly as
    /// [`crate::scoring::Scorer::exact_prefilter`] does.
    fn local_spectrum(&self, query: &ProcessedSpectrum) -> LocalSpectrum {
        let s = &self.settings;
        let precursor = query
            .precursors
            .first()
            .unwrap_or_else(|| panic!("missing MS1 precursor for {}", query.id));
        let mz = precursor.mz - PROTON;

        let mut hypotheses = Vec::new();
        if s.wide_window {
            for charge in s.min_precursor_charge..=s.max_precursor_charge {
                let tolerance = precursor
                    .isolation_window
                    .unwrap_or(Tolerance::Da(-2.4, 2.4))
                    * charge as f32;
                hypotheses.push((mz * charge as f32, charge, tolerance));
            }
        } else if let Some(charge) = precursor.charge.filter(|_| !s.override_precursor_charge) {
            hypotheses.push((mz * charge as f32, charge, s.precursor_tol));
        } else {
            for charge in s.min_precursor_charge..=s.max_precursor_charge {
                hypotheses.push((mz * charge as f32, charge, s.precursor_tol));
            }
        }

        // Peaks among the `max_peaks` most intense, ties broken by position.
        let mut intense = vec![true; query.masses.len()];
        if let Some(max_peaks) = s.max_peaks.filter(|&n| n < query.masses.len()) {
            let mut order = (0..query.masses.len()).collect::<Vec<_>>();
            order.sort_by(|&a, &b| {
                query.intensities[b]
                    .total_cmp(&query.intensities[a])
                    .then(a.cmp(&b))
            });
            for &ix in &order[max_peaks..] {
                intense[ix] = false;
            }
        }

        let mut local = LocalSpectrum::default();
        // Probe ids per fragment-charge limit: unshifted, then one per offset.
        let mut probe_sets: Vec<(u8, Vec<u32>)> = Vec::new();
        for (precursor_mass, charge, tolerance) in hypotheses {
            let fragment_charge = max_fragment_charge(s.max_fragment_charge, charge);
            let probes = match probe_sets.iter().find(|(c, _)| *c == fragment_charge) {
                Some((_, probes)) => probes.clone(),
                None => {
                    let peakset = local.peaksets.len() as u32;
                    // Non-finite peaks never match a fragment in the exact
                    // prefilter, and would break the monotone window order.
                    local.peaksets.push(
                        FragmentMatchIndex::new(query, fragment_charge)
                            .peaks
                            .iter()
                            .filter(|peak| intense[peak.query_index])
                            .map(|peak| peak.neutral_mass)
                            .filter(|mass| mass.is_finite())
                            .collect(),
                    );
                    let shifts = std::iter::once(0.0).chain(self.offsets.iter().map(|o| o.1));
                    let probes = shifts
                        .map(|shift| {
                            local.probes.push(Probe { peakset, shift });
                            local.probes.len() as u32 - 1
                        })
                        .collect::<Vec<_>>();
                    probe_sets.push((fragment_charge, probes.clone()));
                    probes
                }
            };

            for offset in 0..=self.offsets.len() {
                for isotope_error in s.min_isotope_err..=s.max_isotope_err {
                    let mass = precursor_mass - isotope_error as f32 * NEUTRON;
                    let (mass, tolerance) = match offset {
                        0 => (mass, tolerance),
                        _ => offset_query(mass, tolerance, self.offsets[offset - 1].0),
                    };
                    let (lo, hi) = tolerance.bounds(mass);
                    // Windows from non-finite precursor masses contain no
                    // finite peptide mass (NaN bounds compare false).
                    if !(lo < f32::INFINITY && hi > f32::NEG_INFINITY) {
                        continue;
                    }
                    local.windows.push(PrecursorWindow {
                        lo,
                        hi,
                        probes: [
                            probes[0],
                            if offset == 0 {
                                NO_PROBE
                            } else {
                                probes[offset]
                            },
                        ],
                    });
                }
            }
        }
        local
    }

    pub fn finish(self) -> SpectrumIndex {
        self.finish_with(None)
    }

    /// `force_global` overrides the cost model for the global peak index.
    pub(crate) fn finish_with(self, force_global: Option<bool>) -> SpectrumIndex {
        let Self {
            settings,
            mut windows,
            peaks,
            peakset_starts,
            probes,
            spectra,
            ..
        } = self;
        let tolerance = settings.fragment_tol;
        let peakset = |probe: &Probe| {
            &peaks[peakset_starts[probe.peakset as usize] as usize
                ..peakset_starts[probe.peakset as usize + 1] as usize]
        };

        // Fragment windows within a probe share one tolerance and shift, and
        // their bounds are monotone in the peak mass, so peaks sorted by mass
        // have non-decreasing lower and upper bounds.
        assert!(
            probes.par_iter().all(|probe| {
                peakset(probe).windows(2).all(|pair| {
                    let (lo0, hi0) = fragment_window(tolerance, pair[0], probe.shift);
                    let (lo1, hi1) = fragment_window(tolerance, pair[1], probe.shift);
                    lo0 <= lo1 && hi0 <= hi1
                })
            }),
            "internal bug: spectrum index fragment windows are not monotone"
        );

        windows.sort_unstable_by(|a, b| a.lo.total_cmp(&b.lo));
        let window_lo = windows.iter().map(|w| w.lo).collect::<Vec<_>>();
        let window_hi = windows.iter().map(|w| w.hi).collect::<Vec<_>>();
        let window_max_hi = prefix_max(&window_hi);
        let window_probes = windows.iter().map(|w| w.probes).collect::<Vec<_>>();

        let depth = match (window_lo.first(), window_max_hi.last()) {
            (Some(&lo), Some(&hi)) if hi > lo => {
                windows.iter().map(|w| (w.hi - w.lo) as f64).sum::<f64>() / (hi - lo) as f64
            }
            _ => 0.0,
        };

        // Per-probe lookups binary search each candidate spectrum's peaks;
        // global lookups scan every peak window around a fragment mass.
        let windows_total = probes
            .iter()
            .map(|p| peakset(p).len())
            .sum::<usize>()
            .max(1) as f64;
        let probe_cost = (windows_total / probes.len().max(1) as f64).log2().max(1.0) + 1.0;
        let expected_global_cost = {
            let (lo, hi) = (
                peaks.iter().copied().reduce(f32::min),
                peaks.iter().copied().reduce(f32::max),
            );
            let span = match (lo, hi) {
                (Some(lo), Some(hi)) if hi > lo => (hi - lo) as f64,
                _ => 1.0,
            };
            let width = peaks
                .iter()
                .map(|&mass| {
                    let (lo, hi) = tolerance.bounds(mass);
                    (hi - lo) as f64
                })
                .sum::<f64>()
                / peaks.len().max(1) as f64;
            windows_total.log2() + windows_total * width / span
        };
        let build_global =
            force_global.unwrap_or(depth * probe_cost > GLOBAL_INDEX_MARGIN * expected_global_cost);

        let global = build_global.then(|| {
            let mut entries = probes
                .iter()
                .enumerate()
                .flat_map(|(ix, probe)| {
                    peakset(probe).iter().map(move |&mass| {
                        let (lo, hi) = fragment_window(tolerance, mass, probe.shift);
                        (lo, hi, ix as u32)
                    })
                })
                .collect::<Vec<_>>();
            entries.par_sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
            let hi = entries.iter().map(|e| e.1).collect::<Vec<_>>();
            GlobalPeaks {
                lo: entries.iter().map(|e| e.0).collect(),
                max_hi: prefix_max(&hi),
                hi,
                probe: entries.iter().map(|e| e.2).collect(),
            }
        });
        // Measure the actual scan length at a sample of peak masses.
        let global_cost = global.as_ref().map_or(f64::INFINITY, |global| {
            let step = (global.lo.len() / 1024).max(1);
            let samples = (0..global.lo.len()).step_by(step);
            let count = samples.len().max(1) as f64;
            let scanned = samples
                .map(|ix| global.scan_length((global.lo[ix] + global.hi[ix]) / 2.0))
                .sum::<usize>();
            (global.lo.len() as f64).log2() + scanned as f64 / count
        });

        SpectrumIndex {
            fragment_tol: tolerance,
            min_matched_peaks: settings.min_matched_peaks.max(1),
            window_lo,
            window_hi,
            window_max_hi,
            window_probes,
            peaks,
            peakset_starts,
            probes,
            global,
            spectra,
            depth,
            probe_cost,
            global_cost: if force_global == Some(true) {
                0.0
            } else {
                global_cost
            },
        }
    }
}

fn prefix_max(values: &[f32]) -> Vec<f32> {
    let mut max = f32::NEG_INFINITY;
    values
        .iter()
        .map(|&value| {
            max = max.max(value);
            max
        })
        .collect()
}

/// Peak windows of every probe, sorted by lower bound.
struct GlobalPeaks {
    lo: Vec<f32>,
    hi: Vec<f32>,
    /// Running maximum of `hi`, used to stop backward scans.
    max_hi: Vec<f32>,
    probe: Vec<u32>,
}

impl GlobalPeaks {
    /// Number of windows a backward scan visits for `fragment`.
    fn scan_length(&self, fragment: f32) -> usize {
        let end = self
            .lo
            .partition_point(|lo| lo.total_cmp(&fragment).is_le());
        end - self.max_hi[..end].partition_point(|hi| *hi < fragment)
    }

    fn allocated_bytes(&self) -> usize {
        (self.lo.capacity() + self.hi.capacity() + self.max_hi.capacity() + self.probe.capacity())
            * 4
    }
}

pub struct SpectrumIndex {
    fragment_tol: Tolerance,
    min_matched_peaks: u16,
    window_lo: Vec<f32>,
    window_hi: Vec<f32>,
    window_max_hi: Vec<f32>,
    window_probes: Vec<[u32; 2]>,
    peaks: Vec<f32>,
    peakset_starts: Vec<u64>,
    probes: Vec<Probe>,
    global: Option<GlobalPeaks>,
    spectra: usize,
    depth: f64,
    /// Estimated cost of one fragment lookup in one probe.
    probe_cost: f64,
    /// Estimated cost of one fragment lookup in the global peak index.
    global_cost: f64,
}

/// Per-thread scratch space.
struct Scratch {
    stamp: Vec<u32>,
    /// Preliminary matches per probe, valid where `stamp` is current.
    counts: Vec<u16>,
    epoch: u32,
    probes: Vec<u32>,
    /// Precursor windows containing the peptide mass.
    windows: Vec<u32>,
    fragments: Vec<f32>,
}

impl SpectrumIndex {
    pub fn spectra(&self) -> usize {
        self.spectra
    }

    pub fn probes(&self) -> usize {
        self.probes.len()
    }

    pub fn peaks(&self) -> usize {
        self.peaks.len()
    }

    /// Mean number of precursor windows that contain a peptide mass.
    pub fn depth(&self) -> f64 {
        self.depth
    }

    pub fn uses_global_index(&self) -> bool {
        self.global.is_some()
    }

    pub fn allocated_bytes(&self) -> usize {
        (self.window_lo.capacity() + self.window_hi.capacity() + self.window_max_hi.capacity()) * 4
            + self.window_probes.capacity() * 8
            + self.peaks.capacity() * 4
            + self.peakset_starts.capacity() * 8
            + self.probes.capacity() * std::mem::size_of::<Probe>()
            + self.global.as_ref().map_or(0, GlobalPeaks::allocated_bytes)
    }

    pub fn min_matched_peaks(&self) -> u16 {
        self.min_matched_peaks
    }

    fn scratch(&self) -> Scratch {
        Scratch {
            stamp: vec![0; self.probes()],
            counts: vec![0; self.probes()],
            epoch: 0,
            probes: Vec::new(),
            windows: Vec::new(),
            fragments: Vec::new(),
        }
    }

    /// Mark every peptide with at least `min_matched_peaks` preliminary
    /// matches to one precursor hypothesis of any spectrum.
    pub fn filter(&self, parameters: &Parameters, peptides: &[Peptide], keep: &AtomicBitSet) {
        assert_eq!(keep.len(), peptides.len());
        peptides.par_iter().enumerate().for_each_init(
            || self.scratch(),
            |scratch, (index, peptide)| {
                if self.retains(parameters, peptide, scratch) {
                    keep.insert(index);
                }
            },
        );
    }

    fn retains(&self, parameters: &Parameters, peptide: &Peptide, scratch: &mut Scratch) -> bool {
        if !self.candidate_probes(peptide.monoisotopic, scratch) {
            return false;
        }
        scratch.fragments.clear();
        scratch
            .fragments
            .extend(preliminary_fragment_masses(parameters, peptide));
        if self.min_matched_peaks <= 1 {
            return self.any_match(scratch);
        }
        // A window counts the matches of its unshifted and shifted probes,
        // so any single probe reaching the threshold decides early.
        self.count_matches(scratch, self.min_matched_peaks)
            || scratch
                .windows
                .iter()
                .any(|&ix| self.window_count(ix, scratch) >= self.min_matched_peaks)
    }

    /// Stamp the probes of every precursor window containing `mass`, matching
    /// the peptide-range bounds of `IndexedDatabase::query`.
    fn candidate_probes(&self, mass: f32, scratch: &mut Scratch) -> bool {
        scratch.epoch = scratch.epoch.wrapping_add(1);
        if scratch.epoch == 0 {
            scratch.stamp.fill(0);
            scratch.epoch = 1;
        }
        let epoch = scratch.epoch;
        scratch.probes.clear();
        scratch.windows.clear();

        let end = self
            .window_lo
            .partition_point(|lo| lo.total_cmp(&mass).is_le());
        for ix in (0..end).rev() {
            if self.window_max_hi[ix] < mass {
                break;
            }
            if self.window_hi[ix] >= mass {
                scratch.windows.push(ix as u32);
                for &probe in &self.window_probes[ix] {
                    if probe != NO_PROBE && scratch.stamp[probe as usize] != epoch {
                        scratch.stamp[probe as usize] = epoch;
                        scratch.counts[probe as usize] = 0;
                        scratch.probes.push(probe);
                    }
                }
            }
        }
        !scratch.probes.is_empty()
    }

    fn use_global(&self, scratch: &Scratch) -> Option<&GlobalPeaks> {
        self.global
            .as_ref()
            .filter(|_| scratch.probes.len() as f64 * self.probe_cost > self.global_cost)
    }

    fn peaks_of(&self, probe: Probe) -> &[f32] {
        &self.peaks[self.peakset_starts[probe.peakset as usize] as usize
            ..self.peakset_starts[probe.peakset as usize + 1] as usize]
    }

    /// Whether any fragment matches any stamped probe.
    fn any_match(&self, scratch: &Scratch) -> bool {
        let epoch = scratch.epoch;
        match self.use_global(scratch) {
            Some(global) => scratch.fragments.iter().any(|&fragment| {
                let end = global
                    .lo
                    .partition_point(|lo| lo.total_cmp(&fragment).is_le());
                for ix in (0..end).rev() {
                    if global.max_hi[ix] < fragment {
                        return false;
                    }
                    if global.hi[ix] >= fragment
                        && scratch.stamp[global.probe[ix] as usize] == epoch
                    {
                        return true;
                    }
                }
                false
            }),
            None => scratch.probes.iter().any(|&probe| {
                let probe = self.probes[probe as usize];
                let peaks = self.peaks_of(probe);
                let window = |mass| fragment_window(self.fragment_tol, mass, probe.shift);
                scratch.fragments.iter().any(|&fragment| {
                    // Bounds are monotone within a probe, so the last window
                    // starting at or below the fragment has the largest upper
                    // bound among all candidates.
                    let end = peaks.partition_point(|&mass| window(mass).0 <= fragment);
                    end > 0 && window(peaks[end - 1]).1 >= fragment
                })
            }),
        }
    }

    /// Count (peak, fragment) matches per stamped probe. Returns true as soon
    /// as one probe reaches `stop`, leaving other counts partial.
    fn count_matches(&self, scratch: &mut Scratch, stop: u16) -> bool {
        let epoch = scratch.epoch;
        match self.use_global(scratch) {
            Some(global) => {
                for &fragment in &scratch.fragments {
                    let end = global
                        .lo
                        .partition_point(|lo| lo.total_cmp(&fragment).is_le());
                    for ix in (0..end).rev() {
                        if global.max_hi[ix] < fragment {
                            break;
                        }
                        let probe = global.probe[ix] as usize;
                        if global.hi[ix] >= fragment && scratch.stamp[probe] == epoch {
                            let count = scratch.counts[probe].saturating_add(1);
                            scratch.counts[probe] = count;
                            if count >= stop {
                                return true;
                            }
                        }
                    }
                }
                false
            }
            None => {
                for &probe_ix in &scratch.probes {
                    let probe = self.probes[probe_ix as usize];
                    let peaks = self.peaks_of(probe);
                    let window = |mass| fragment_window(self.fragment_tol, mass, probe.shift);
                    // Lower and upper bounds are both monotone within a probe,
                    // so the windows containing a fragment are contiguous.
                    let count = scratch.fragments.iter().fold(0u16, |count, &fragment| {
                        let end = peaks.partition_point(|&mass| window(mass).0 <= fragment);
                        let start = peaks[..end].partition_point(|&mass| window(mass).1 < fragment);
                        count.saturating_add((end - start).min(u16::MAX as usize) as u16)
                    });
                    scratch.counts[probe_ix as usize] = count;
                    if count >= stop {
                        return true;
                    }
                }
                false
            }
        }
    }

    fn window_count(&self, ix: u32, scratch: &Scratch) -> u16 {
        let [unshifted, shifted] = self.window_probes[ix as usize];
        let count = scratch.counts[unshifted as usize];
        match shifted {
            NO_PROBE => count,
            shifted => count.saturating_add(scratch.counts[shifted as usize]),
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/spectrum_index.rs"]
mod test;
