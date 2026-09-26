use super::*;
use crate::database::{Builder, EnzymeBuilder, IndexedDatabase};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::modification::{
    NeutralLossMode, SearchMode, SiteMode, StaticModEntry, VarModEntry, VariableModification,
};
use crate::scoring::{offset_query, ScoreType, Scorer};
use crate::spectrum::Precursor;
use std::collections::HashMap;

const FASTA: &str = ">P1\nMAGSPEPTSIDEKLLSAYGNRWTTPEGSARMKVLDEFGHIKR\n\
>P2\nGGSTVLAPEDKAAAARNMSTYWPLLKEEGHCTMSPARQQWNNLK\n\
>P3\nLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEK\n\
>P4\nYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCK\n";

/// Deterministic generator so failures are reproducible.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn unit(&mut self) -> f32 {
        (self.next() % 1_000_000) as f32 / 1_000_000.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn parameters(offset: bool) -> Parameters {
    let phospho = VarModEntry::Detailed(VariableModification {
        mass: 79.966_33,
        max_count: Some(1),
        max_total_count: None,
        name: Some("Phospho".into()),
        neutral_losses: vec![97.976_9],
        neutral_loss_mode: NeutralLossMode::Required,
        site_mode: SiteMode::Exhaustive,
        search_mode: if offset {
            SearchMode::MassOffset
        } else {
            SearchMode::Database
        },
        channel_offsets: Default::default(),
    });
    let mut variable_mods = HashMap::new();
    for residue in ["S", "T", "Y"] {
        variable_mods.insert(residue.to_string(), vec![phospho.clone()]);
    }
    variable_mods.insert("M".to_string(), vec![VarModEntry::Mass(15.994_915)]);
    Builder {
        max_variable_mods: Some(2),
        variable_mods: Some(variable_mods),
        static_mods: Some(HashMap::from([(
            "C".to_string(),
            StaticModEntry::Mass(57.021_464),
        )])),
        enzyme: Some(EnzymeBuilder {
            min_len: Some(5),
            missed_cleavages: Some(2),
            ..Default::default()
        }),
        ..Builder::default()
    }
    .make_parameters()
}

fn spectra(database: &IndexedDatabase, count: usize, seed: u64) -> Vec<ProcessedSpectrum> {
    let mut rng = Lcg(seed);
    (0..count)
        .map(|index| {
            let peptide = &database.peptides[rng.below(database.peptides.len())];
            let mut masses = Vec::new();
            let mut charges = Vec::new();
            let mut known = Vec::new();
            for kind in [Kind::B, Kind::Y] {
                for ion in IonSeries::new(peptide, kind) {
                    if rng.unit() < 0.3 {
                        continue;
                    }
                    // Place peaks around the 20 ppm fragment tolerance edge.
                    let ppm = (rng.unit() - 0.5) * 50.0;
                    let charge = 1 + rng.below(2) as u8;
                    masses.push(ion.monoisotopic_mass * (1.0 + ppm / 1e6) / charge as f32);
                    charges.push(charge);
                    known.push(rng.unit() < 0.5);
                }
            }
            for _ in 0..20 {
                masses.push(100.0 + rng.unit() * 1500.0);
                charges.push(1);
                known.push(false);
            }
            let charge = 2 + rng.below(3) as u8;
            let isotope = rng.below(4) as f32 - 1.0;
            let shift = if rng.unit() < 0.3 { 79.966_33 } else { 0.0 };
            let ppm = (rng.unit() - 0.5) * 30.0;
            let mass = (peptide.monoisotopic + shift + isotope * NEUTRON) * (1.0 + ppm / 1e6);
            let mut order = (0..masses.len()).collect::<Vec<_>>();
            order.sort_by(|&a, &b| masses[a].total_cmp(&masses[b]));
            ProcessedSpectrum {
                level: 2,
                id: format!("scan={index}"),
                precursors: vec![Precursor {
                    mz: mass / charge as f32 + PROTON,
                    charge: (rng.unit() < 0.6).then_some(charge),
                    isolation_window: Some(Tolerance::Da(-0.8, 0.8)),
                    ..Precursor::default()
                }],
                intensities: vec![1.0; masses.len()],
                total_ion_current: masses.len() as f32,
                masses: order.iter().map(|&i| masses[i]).collect(),
                charges: order.iter().map(|&i| charges[i]).collect(),
                charge_is_known: order.iter().map(|&i| known[i]).collect(),
                ..ProcessedSpectrum::default()
            }
        })
        .collect()
}

fn settings(
    precursor_tol: Tolerance,
    wide_window: bool,
    override_charge: bool,
) -> SpectrumIndexSettings {
    SpectrumIndexSettings {
        precursor_tol,
        fragment_tol: Tolerance::Ppm(-20.0, 20.0),
        min_isotope_err: -1,
        max_isotope_err: 2,
        min_precursor_charge: 2,
        max_precursor_charge: 4,
        override_precursor_charge: override_charge,
        max_fragment_charge: Some(2),
        wide_window,
        min_peaks: 1,
        min_matched_peaks: 1,
        max_peaks: None,
    }
}

fn scorer<'db>(database: &'db IndexedDatabase, s: &SpectrumIndexSettings) -> Scorer<'db> {
    Scorer {
        db: database,
        precursor_tol: s.precursor_tol,
        fragment_tol: s.fragment_tol,
        min_matched_peaks: 1,
        min_isotope_err: s.min_isotope_err,
        max_isotope_err: s.max_isotope_err,
        min_precursor_charge: s.min_precursor_charge,
        max_precursor_charge: s.max_precursor_charge,
        override_precursor_charge: s.override_precursor_charge,
        max_fragment_charge: s.max_fragment_charge,
        chimera: false,
        report_psms: 2,
        wide_window: s.wide_window,
        annotate_matches: false,
        mass_shift_ppm: crate::ambiguity::DEFAULT_MASS_SHIFT_PPM,
        score_type: ScoreType::SageHyperScore,
        mass_recalibration: None,
    }
}

fn assert_same_survivors(offset: bool, s: SpectrumIndexSettings, open: bool) {
    assert_same_survivors_with(offset, s, open, |_| {});
}

fn assert_same_survivors_with(
    offset: bool,
    s: SpectrumIndexSettings,
    open: bool,
    corrupt: impl Fn(&mut [ProcessedSpectrum]),
) {
    let parameters = parameters(offset);
    let fasta = Fasta::parse(FASTA.into(), "rev_", true).unwrap();
    let database = parameters.clone().build(fasta);
    let mut queries = spectra(&database, 400, 0x5eed ^ offset as u64);
    corrupt(&mut queries);

    let expected = AtomicBitSet::new(database.peptides.len());
    let classic = scorer(&database, &s);
    for query in &queries {
        classic.exact_prefilter(query, &expected);
    }
    let expected = (0..database.peptides.len())
        .filter(|&ix| expected.contains(ix))
        .collect::<Vec<_>>();
    assert!(!expected.is_empty());
    if !open {
        assert!(expected.len() < database.peptides.len());
    }

    // Both lookup paths must agree with the classic prefilter.
    for global in [false, true] {
        let mut builder = SpectrumIndexBuilder::new(s.clone(), &database.mass_offsets);
        // Batches must combine exactly like a single pass.
        let (first, second) = queries.split_at(queries.len() / 3);
        builder.add(first);
        builder.add(second);
        let index = builder.finish_with(Some(global));
        assert_eq!(index.uses_global_index(), global);
        let actual = AtomicBitSet::new(database.peptides.len());
        index.filter(&parameters, &database.peptides, &actual);
        let actual = (0..database.peptides.len())
            .filter(|&ix| actual.contains(ix))
            .collect::<Vec<_>>();
        assert_eq!(expected, actual, "global index: {global}");
    }

    // Separate indexes over spectrum batches accumulate the same survivors.
    let actual = AtomicBitSet::new(database.peptides.len());
    for batch in queries.chunks(150) {
        let mut builder = SpectrumIndexBuilder::new(s.clone(), &database.mass_offsets);
        builder.add(batch);
        builder
            .finish()
            .filter(&parameters, &database.peptides, &actual);
    }
    let actual = (0..database.peptides.len())
        .filter(|&ix| actual.contains(ix))
        .collect::<Vec<_>>();
    assert_eq!(expected, actual, "batched indexes");
}

#[test]
fn closed_ppm_search_matches_classic_prefilter() {
    assert_same_survivors(
        false,
        settings(Tolerance::Ppm(-10.0, 10.0), false, false),
        false,
    );
}

#[test]
fn override_charge_matches_classic_prefilter() {
    assert_same_survivors(
        false,
        settings(Tolerance::Ppm(-10.0, 10.0), false, true),
        false,
    );
}

#[test]
fn mass_offsets_match_classic_prefilter() {
    assert_same_survivors(
        true,
        settings(Tolerance::Ppm(-10.0, 10.0), false, false),
        false,
    );
}

#[test]
fn da_tolerance_matches_classic_prefilter() {
    assert_same_survivors(
        true,
        settings(Tolerance::Da(-0.02, 0.02), false, false),
        false,
    );
}

#[test]
fn wide_window_matches_classic_prefilter() {
    assert_same_survivors(
        false,
        settings(Tolerance::Ppm(-10.0, 10.0), true, false),
        false,
    );
}

#[test]
fn open_search_matches_classic_prefilter() {
    assert_same_survivors(
        true,
        settings(Tolerance::Da(-50.0, 150.0), false, false),
        true,
    );
}

/// Non-finite peaks never match in the exact prefilter; the spectrum index
/// must ignore them too instead of failing its monotonicity check.
fn corrupt_with_non_finite_masses(queries: &mut [ProcessedSpectrum]) {
    for (index, query) in queries.iter_mut().enumerate().step_by(7) {
        // Keep masses sorted under `total_cmp`, which places these last.
        for mass in [f32::INFINITY, f32::NAN] {
            query.masses.push(mass);
            query.intensities.push(1.0);
            query.charges.push(1);
            query.charge_is_known.push(index % 2 == 0);
        }
        if index % 3 == 0 {
            query.precursors[0].mz = if index % 2 == 0 {
                f32::NAN
            } else {
                f32::INFINITY
            };
        }
    }
}

#[test]
fn non_finite_masses_match_classic_prefilter() {
    for tolerance in [Tolerance::Ppm(-10.0, 10.0), Tolerance::Da(-0.02, 0.02)] {
        for offset in [false, true] {
            assert_same_survivors_with(
                offset,
                settings(tolerance, false, false),
                false,
                corrupt_with_non_finite_masses,
            );
        }
    }
    assert_same_survivors_with(
        false,
        settings(Tolerance::Ppm(-10.0, 10.0), true, false),
        false,
        corrupt_with_non_finite_masses,
    );
}

#[test]
fn wide_windows_select_the_global_index() {
    let parameters = parameters(false);
    let fasta = Fasta::parse(FASTA.into(), "rev_", true).unwrap();
    let database = parameters.build(fasta);
    let queries = spectra(&database, 400, 7);
    let build = |tolerance| {
        let mut builder =
            SpectrumIndexBuilder::new(settings(tolerance, false, false), &database.mass_offsets);
        builder.add(&queries);
        builder.finish()
    };
    assert!(!build(Tolerance::Ppm(-10.0, 10.0)).uses_global_index());
    assert!(build(Tolerance::Da(-50.0, 150.0)).uses_global_index());
}

/// Best preliminary match count of every peptide over the precursor
/// hypotheses of `query`, counted with the fragment index as the search does.
fn classic_counts(
    database: &IndexedDatabase,
    s: &SpectrumIndexSettings,
    query: &ProcessedSpectrum,
    best: &mut [u16],
) {
    let precursor = &query.precursors[0];
    let mz = precursor.mz - PROTON;
    let mut hypotheses = Vec::new();
    if s.wide_window {
        for charge in s.min_precursor_charge..=s.max_precursor_charge {
            let tolerance = precursor.isolation_window.unwrap() * charge as f32;
            hypotheses.push((mz * charge as f32, charge, tolerance));
        }
    } else if let Some(charge) = precursor.charge.filter(|_| !s.override_precursor_charge) {
        hypotheses.push((mz * charge as f32, charge, s.precursor_tol));
    } else {
        for charge in s.min_precursor_charge..=s.max_precursor_charge {
            hypotheses.push((mz * charge as f32, charge, s.precursor_tol));
        }
    }

    let mut order = (0..query.masses.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| {
        query.intensities[b]
            .total_cmp(&query.intensities[a])
            .then(a.cmp(&b))
    });
    let mut intense = vec![false; query.masses.len()];
    for &ix in order.iter().take(s.max_peaks.unwrap_or(usize::MAX)) {
        intense[ix] = true;
    }

    for (precursor_mass, charge, tolerance) in hypotheses {
        if !precursor_mass.is_finite() {
            continue;
        }
        let peaks =
            FragmentMatchIndex::new(query, max_fragment_charge(s.max_fragment_charge, charge))
                .peaks
                .into_iter()
                .filter(|peak| intense[peak.query_index])
                .collect::<Vec<_>>();
        for offset in 0..=database.mass_offsets.len() {
            for isotope_error in s.min_isotope_err..=s.max_isotope_err {
                let mass = precursor_mass - isotope_error as f32 * NEUTRON;
                let (mass, tolerance) = match offset {
                    0 => (mass, tolerance),
                    _ => offset_query(mass, tolerance, database.mass_offsets[offset - 1].mass()),
                };
                let candidates = database.query(mass, tolerance, s.fragment_tol);
                let mut counts = HashMap::<usize, u16>::new();
                for peak in &peaks {
                    for fragment in candidates.page_search(peak.neutral_mass) {
                        *counts.entry(fragment.peptide_index.0 as usize).or_default() += 1;
                    }
                    if offset > 0 {
                        let shift = database.mass_offsets[offset - 1].fragment_shift();
                        for fragment in candidates.page_search_shifted(peak.neutral_mass, shift) {
                            *counts.entry(fragment.peptide_index.0 as usize).or_default() += 1;
                        }
                    }
                }
                for (peptide, count) in counts {
                    best[peptide] = best[peptide].max(count);
                }
            }
        }
    }
}

/// Survivors of every lookup path must equal the peptides whose best
/// preliminary count reaches the threshold.
fn assert_threshold_survivors(offset: bool, base: SpectrumIndexSettings) {
    let parameters = parameters(offset);
    let fasta = Fasta::parse(FASTA.into(), "rev_", true).unwrap();
    let database = parameters.clone().build(fasta);
    let mut queries = spectra(&database, 400, 0xc0de ^ offset as u64);
    let mut rng = Lcg(11);
    for query in &mut queries {
        for intensity in &mut query.intensities {
            *intensity = rng.unit();
        }
    }

    for max_peaks in [None, Some(12)] {
        let mut best = vec![0u16; database.peptides.len()];
        let s = SpectrumIndexSettings {
            max_peaks,
            ..base.clone()
        };
        for query in &queries {
            classic_counts(&database, &s, query, &mut best);
        }
        let mut previous = usize::MAX;
        for min_matched_peaks in [1u16, 2, 3, 5] {
            let s = SpectrumIndexSettings {
                min_matched_peaks,
                ..s.clone()
            };
            let expected = (0..best.len())
                .filter(|&ix| best[ix] >= min_matched_peaks)
                .collect::<Vec<_>>();
            assert!(!expected.is_empty(), "threshold {min_matched_peaks}");
            assert!(expected.len() <= previous);
            previous = expected.len();

            let survivors = |keep: &AtomicBitSet| {
                (0..database.peptides.len())
                    .filter(|&ix| keep.contains(ix))
                    .collect::<Vec<_>>()
            };
            for global in [false, true] {
                let mut builder = SpectrumIndexBuilder::new(s.clone(), &database.mass_offsets);
                builder.add(&queries);
                let index = builder.finish_with(Some(global));
                let actual = AtomicBitSet::new(database.peptides.len());
                index.filter(&parameters, &database.peptides, &actual);
                assert_eq!(
                    expected,
                    survivors(&actual),
                    "threshold {min_matched_peaks}, max peaks {max_peaks:?}, global {global}"
                );
            }
            let actual = AtomicBitSet::new(database.peptides.len());
            for batch in queries.chunks(150) {
                let mut builder = SpectrumIndexBuilder::new(s.clone(), &database.mass_offsets);
                builder.add(batch);
                builder
                    .finish()
                    .filter(&parameters, &database.peptides, &actual);
            }
            assert_eq!(expected, survivors(&actual), "batched indexes");
        }
    }
}

#[test]
fn match_threshold_closed_search() {
    assert_threshold_survivors(false, settings(Tolerance::Ppm(-10.0, 10.0), false, false));
}

#[test]
fn match_threshold_mass_offsets() {
    assert_threshold_survivors(true, settings(Tolerance::Da(-0.02, 0.02), false, false));
}

#[test]
fn match_threshold_wide_window() {
    assert_threshold_survivors(false, settings(Tolerance::Ppm(-10.0, 10.0), true, false));
}

#[test]
fn match_threshold_open_search() {
    assert_threshold_survivors(true, settings(Tolerance::Da(-50.0, 150.0), false, false));
}
