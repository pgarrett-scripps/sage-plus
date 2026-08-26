//! Deterministic scorer benchmark for charge-aware fragment matching.
//! Run with `cargo run --release -p sage-core --example charge_matching_benchmark`.

use sage_core::database::Builder;
use sage_core::enzyme::Digest;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::{Tolerance, NEUTRON, PROTON};
use sage_core::peptide::Peptide;
use sage_core::scoring::{ScoreType, Scorer};
use sage_core::spectrum::{
    DeisotopeSettings, Precursor, RawSpectrum, Representation, SpectrumProcessor,
};
use std::hint::black_box;
use std::sync::Arc;
use std::time::Instant;

const PEPTIDE_COUNT: usize = 25_000;
const SPECTRUM_COUNT: usize = 160;
const ROUNDS: usize = 25;
const PREPROCESS_ROUNDS: usize = 100;

fn sequence(mut value: usize) -> String {
    const ALPHABET: &[u8] = b"ACDEFGHIKMNPQRSTVWY";
    let mut state = value as u64 ^ 0x9E37_79B9_7F4A_7C15;
    let mut sequence = vec![b'A'; 15];
    for residue in &mut sequence[..14] {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        *residue = ALPHABET[state as usize % ALPHABET.len()];
    }
    for position in [0, 5, 10, 1] {
        sequence[position] = ALPHABET[value % ALPHABET.len()];
        value /= ALPHABET.len();
    }
    sequence[14] = b'K';
    String::from_utf8(sequence).unwrap()
}

fn peptide(index: usize) -> Peptide {
    Peptide::try_from(Digest {
        sequence: sequence(index),
        protein: Arc::from(format!("protein-{index}")),
        ..Digest::default()
    })
    .unwrap()
}

fn raw_spectrum(peptide: &Peptide, index: usize) -> RawSpectrum {
    let mut peaks = Vec::new();
    for (kind_offset, kind) in [Kind::B, Kind::Y].into_iter().enumerate() {
        for (ordinal, ion) in IonSeries::new(peptide, kind).enumerate() {
            let charge = if (ordinal + kind_offset) % 2 == 0 {
                1
            } else {
                2
            };
            let mz = ion.monoisotopic_mass / charge as f32 + PROTON;
            let intensity = 1_000.0 - ordinal as f32 * 7.0;
            peaks.push((mz, intensity));
            peaks.push((mz + NEUTRON / charge as f32, intensity * 0.35));
        }
    }
    for noise in 0..30 {
        let mz = 175.0 + ((index * 37 + noise * 53) % 1_250) as f32;
        peaks.push((mz, 5.0 + noise as f32));
    }
    peaks.sort_by(|left, right| left.0.total_cmp(&right.0));
    let (mz, intensity) = peaks.into_iter().unzip();

    RawSpectrum {
        ms_level: 2,
        id: format!("synthetic-{index}"),
        representation: Representation::Centroid,
        precursors: vec![Precursor {
            mz: peptide.monoisotopic / 3.0 + PROTON,
            charge: Some(3),
            ..Precursor::default()
        }],
        mz,
        intensity,
        ..RawSpectrum::default()
    }
}

fn main() {
    let peptides = (0..PEPTIDE_COUNT).map(peptide).collect::<Vec<_>>();
    let targets = (0..SPECTRUM_COUNT)
        .map(|index| peptides[(index * 149) % peptides.len()].clone())
        .collect::<Vec<_>>();
    let database = Builder::default()
        .make_parameters()
        .build_from_peptides(peptides);
    let mode = std::env::var("SAGE_DEISOTOPE_MODE").unwrap_or_else(|_| "legacy".into());
    let processor = match mode.as_str() {
        "averagine" => {
            let mut settings = DeisotopeSettings::default();
            if let Ok(value) = std::env::var("SAGE_MIN_ISOTOPE_SCORE") {
                settings.min_score = value.parse().expect("valid SAGE_MIN_ISOTOPE_SCORE");
            }
            if let Ok(value) = std::env::var("SAGE_MAX_ISOTOPE_LOG2_RATIO") {
                settings.max_isotope_log2_ratio =
                    value.parse().expect("valid SAGE_MAX_ISOTOPE_LOG2_RATIO");
            }
            SpectrumProcessor::with_deisotope_settings(150, settings, 0.0)
        }
        "legacy" => SpectrumProcessor::new(150, true, 0.0),
        other => panic!("unknown SAGE_DEISOTOPE_MODE {other}"),
    };
    let raw_spectra = targets
        .iter()
        .enumerate()
        .map(|(index, peptide)| raw_spectrum(peptide, index))
        .collect::<Vec<_>>();
    let spectra = raw_spectra
        .iter()
        .cloned()
        .map(|spectrum| processor.process(spectrum))
        .collect::<Vec<_>>();
    let scorer = Scorer {
        db: &database,
        precursor_tol: Tolerance::Da(-5.0, 5.0),
        fragment_tol: Tolerance::Ppm(-10.0, 10.0),
        min_matched_peaks: 4,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 2,
        max_precursor_charge: 4,
        override_precursor_charge: false,
        max_fragment_charge: None,
        chimera: false,
        report_psms: 1,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: 50.0,
        score_type: ScoreType::SageHyperScore,
    };

    if std::env::var_os("SAGE_BENCH_DETAILS").is_some() {
        for (index, spectrum) in spectra.iter().enumerate() {
            let features = scorer.score(spectrum);
            if let Some(feature) = features.first() {
                println!(
                    "{index}\t{}\t{}\t{}\t{}\t{}",
                    feature.peptide_idx.0,
                    feature.matched_peaks,
                    feature.hyperscore.to_bits(),
                    feature.calcmass.to_bits(),
                    targets[index].monoisotopic.to_bits()
                );
            } else {
                println!("{index}\tNONE");
            }
        }
        return;
    }

    for spectrum in &raw_spectra {
        black_box(processor.process(black_box(spectrum.clone())));
    }
    let preprocess_start = Instant::now();
    let mut processed_peaks = 0usize;
    let mut known_charges = 0usize;
    for _ in 0..PREPROCESS_ROUNDS {
        for spectrum in &raw_spectra {
            let processed = black_box(processor.process(black_box(spectrum.clone())));
            processed_peaks += processed.len();
            known_charges += processed
                .charge_is_known
                .iter()
                .filter(|&&known| known)
                .count();
        }
    }
    let preprocess_elapsed = preprocess_start.elapsed();

    for spectrum in &spectra {
        black_box(scorer.score(black_box(spectrum)));
    }

    let start = Instant::now();
    let mut psms = 0usize;
    let mut matched_peaks = 0u64;
    let mut peptide_checksum = 0u64;
    let mut score_checksum = 0u64;
    for _ in 0..ROUNDS {
        for spectrum in &spectra {
            let features = black_box(scorer.score(black_box(spectrum)));
            psms += features.len();
            for feature in features {
                matched_peaks += u64::from(feature.matched_peaks);
                peptide_checksum = peptide_checksum.wrapping_add(u64::from(feature.peptide_idx.0));
                score_checksum = score_checksum.wrapping_add(feature.hyperscore.to_bits());
            }
        }
    }
    let elapsed = start.elapsed();
    println!("deisotope_mode={mode}");
    println!("preprocess_rounds={PREPROCESS_ROUNDS}");
    println!(
        "preprocess_ns_per_spectrum={:.1}",
        preprocess_elapsed.as_nanos() as f64 / (SPECTRUM_COUNT * PREPROCESS_ROUNDS) as f64
    );
    println!("processed_peaks={processed_peaks}");
    println!("known_charges={known_charges}");
    println!("peptides={PEPTIDE_COUNT}");
    println!("spectra={SPECTRUM_COUNT}");
    println!("rounds={ROUNDS}");
    println!("searches={}", SPECTRUM_COUNT * ROUNDS);
    println!("elapsed_ms={:.3}", elapsed.as_secs_f64() * 1_000.0);
    println!(
        "ns_per_search={:.1}",
        elapsed.as_nanos() as f64 / (SPECTRUM_COUNT * ROUNDS) as f64
    );
    println!("psms={psms}");
    println!("matched_peaks={matched_peaks}");
    println!("peptide_checksum={peptide_checksum}");
    println!("score_checksum={score_checksum}");
}
