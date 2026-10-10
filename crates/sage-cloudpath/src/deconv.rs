//! SCOPING PROTOTYPE (b18/mzdeisotope-scope): mzdeisotope / mzsignal hooks.
//!
//! Not a supported feature. `apply_from_env` is called after mzML reading and
//! is controlled by `SAGE_PLUS_MZDEISOTOPE`:
//!   unset        -> only profile spectra are centroided with mzsignal
//!   `ms2`        -> MS2 peaks replaced by mzdeisotope monoisotopic peaks + charge array
//!   `ms2+prec`   -> also re-fit precursor monoisotopic m/z from the preceding MS1
//!   `prec`       -> only the precursor re-fit
use mzdeisotope::{
    isotopic_model::{CachingIsotopicModel, IsotopicModels, IsotopicPatternParams, PROTON},
    scorer::{MSDeconvScorer, MaximizingFitFilter, PenalizedMSDeconvScorer},
    DeconvolutionEngine, PeaksAndTargets,
};
use mzpeaks::{CentroidPeak, MZPeakSetType, Tolerance};
use rayon::prelude::*;
use sage_core::spectrum::{RawSpectrum, Representation};

pub type Ms2Engine = DeconvolutionEngine<'static, CentroidPeak, MSDeconvScorer, MaximizingFitFilter>;
pub type Ms1Engine =
    DeconvolutionEngine<'static, CentroidPeak, PenalizedMSDeconvScorer, MaximizingFitFilter>;

/// MSn engine with the parameters of mzdeisotoper's defaults.
pub fn ms2_engine() -> Ms2Engine {
    let mut engine = DeconvolutionEngine::new(
        IsotopicPatternParams::new(0.8, 0.001, None, PROTON),
        CachingIsotopicModel::from(IsotopicModels::Peptide),
        MSDeconvScorer::default(),
        MaximizingFitFilter::new(10.0),
        true,
    );
    engine.populate_isotopic_model_cache(80.0, 3000.0, 1, 8);
    engine
}

pub fn ms1_engine() -> Ms1Engine {
    let mut engine = DeconvolutionEngine::new(
        Default::default(),
        CachingIsotopicModel::from(IsotopicModels::Peptide),
        PenalizedMSDeconvScorer::new(0.02, 2.0),
        MaximizingFitFilter::new(10.0),
        true,
    );
    engine.populate_isotopic_model_cache(80.0, 3000.0, 1, 8);
    engine
}

fn centroids(mz: &[f32], intensity: &[f32]) -> MZPeakSetType<CentroidPeak> {
    mz.iter()
        .zip(intensity)
        .enumerate()
        .map(|(i, (&mz, &int))| CentroidPeak::new(mz as f64, int, i as u32))
        .collect()
}

/// Deconvolved MS2 peaks as (m/z at the assigned charge, intensity, charge), sorted by m/z.
pub fn deconvolve_ms2(
    engine: &mut Ms2Engine,
    spectrum: &RawSpectrum,
    max_charge: i32,
) -> Vec<(f32, f32, u8)> {
    let peaks = centroids(&spectrum.mz, &spectrum.intensity);
    let Ok(out) = engine.deconvolute_peaks(peaks, Tolerance::PPM(20.0), (1, max_charge), 1) else {
        return Vec::new();
    };
    let mut v = out
        .iter()
        .map(|p| {
            let z = p.charge.unsigned_abs().max(1) as f64;
            (
                (p.neutral_mass / z + PROTON) as f32,
                p.intensity,
                p.charge.unsigned_abs() as u8,
            )
        })
        .collect::<Vec<_>>();
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    v
}

/// Centroid a profile spectrum with mzsignal's Gaussian peak fitter.
pub fn pick_profile(spectrum: &mut RawSpectrum) {
    let mz = spectrum.mz.iter().map(|&m| m as f64).collect::<Vec<_>>();
    let picker = mzsignal::PeakPicker::default();
    let mut acc = Vec::new();
    if picker
        .discover_peaks(&mz, &spectrum.intensity, &mut acc)
        .is_err()
    {
        return;
    }
    spectrum.mz = acc.iter().map(|p| p.mz as f32).collect();
    spectrum.intensity = acc.iter().map(|p| p.intensity).collect();
    spectrum.fragment_charges = None;
    spectrum.representation = Representation::Centroid;
}

#[derive(Debug, Default, Clone, Copy)]
pub struct PrecursorStats {
    pub targets: usize,
    pub fitted: usize,
    pub same: usize,
    pub mono_shift: usize,
    pub charge_differs: usize,
    pub charge_filled: usize,
    pub applied: usize,
}

/// Re-fit precursor monoisotopic m/z from the preceding MS1 scan, applying the
/// result only when the charge is missing or agrees (mzdeisotoper's rule).
pub fn refit_precursors(spectra: &mut [RawSpectrum], apply: bool) -> PrecursorStats {
    // group: MS1 index followed by MS2 indices
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, s) in spectra.iter().enumerate() {
        if s.ms_level == 1 {
            groups.push((i, Vec::new()));
        } else if s.ms_level == 2 {
            if let Some(g) = groups.last_mut() {
                g.1.push(i);
            }
        }
    }
    let engine = ms1_engine();
    let results = groups
        .par_iter()
        .map_init(
            || engine.clone(),
            |engine, (ms1, ms2s)| {
                let targets_mz = ms2s
                    .iter()
                    .map(|&i| spectra[i].precursors.first().map(|p| p.mz as f64).unwrap_or(0.0))
                    .collect::<Vec<_>>();
                let peaks = centroids(&spectra[*ms1].mz, &spectra[*ms1].intensity);
                match engine.deconvolute_peaks_with_targets(
                    peaks,
                    Tolerance::PPM(10.0),
                    (1, 8),
                    1,
                    &targets_mz,
                ) {
                    Ok(PeaksAndTargets { targets, .. }) => ms2s
                        .iter()
                        .zip(targets)
                        .map(|(&i, t)| {
                            (
                                i,
                                t.map(|p| {
                                    let z = p.charge.unsigned_abs().max(1);
                                    ((p.neutral_mass / z as f64 + PROTON) as f32, z as u8)
                                }),
                            )
                        })
                        .collect::<Vec<_>>(),
                    Err(_) => ms2s.iter().map(|&i| (i, None)).collect(),
                }
            },
        )
        .flatten()
        .collect::<Vec<_>>();

    let mut stats = PrecursorStats::default();
    for (i, fit) in results {
        stats.targets += 1;
        let Some((mz, z)) = fit else { continue };
        stats.fitted += 1;
        let Some(prec) = spectra[i].precursors.first_mut() else {
            continue;
        };
        let charge_ok = match prec.charge {
            None => {
                stats.charge_filled += 1;
                true
            }
            Some(c) if c == z => true,
            Some(_) => {
                stats.charge_differs += 1;
                false
            }
        };
        if charge_ok {
            let shift = (mz - prec.mz).abs() * f32::from(z);
            if shift < 0.02 {
                stats.same += 1;
            } else {
                stats.mono_shift += 1;
            }
            if apply {
                prec.mz = mz;
                prec.charge = Some(z);
                stats.applied += 1;
            }
        }
    }
    stats
}

/// Replace MS2 peak lists with mzdeisotope output (charge array set).
pub fn deconvolve_all_ms2(spectra: &mut [RawSpectrum]) {
    let engine = ms2_engine();
    spectra.par_iter_mut().for_each_init(
        || engine.clone(),
        |engine, s| {
            if s.ms_level != 2 {
                return;
            }
            let max_z = s
                .precursors
                .first()
                .and_then(|p| p.charge)
                .map(i32::from)
                .unwrap_or(3);
            let out = deconvolve_ms2(engine, s, max_z);
            s.mz = out.iter().map(|p| p.0).collect();
            s.intensity = out.iter().map(|p| p.1).collect();
            s.fragment_charges = Some(out.iter().map(|p| p.2).collect());
        },
    );
}

pub fn apply_from_env(spectra: &mut [RawSpectrum]) {
    spectra
        .par_iter_mut()
        .filter(|s| s.representation == Representation::Profile)
        .for_each(pick_profile);
    let mode = std::env::var("SAGE_PLUS_MZDEISOTOPE").unwrap_or_default();
    if mode.contains("prec") {
        let stats = refit_precursors(spectra, true);
        log::info!("mzdeisotope precursor refit: {stats:?}");
    }
    if mode.contains("ms2") {
        deconvolve_all_ms2(spectra);
        log::info!("mzdeisotope: deconvolved MS2 peak lists");
    }
}
