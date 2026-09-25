use super::*;
use crate::{database::Builder, scoring::Feature};

#[test]
fn default_rt_tolerance_preserves_existing_window() {
    let mut settings = LfqSettings::default();

    assert_eq!(settings.rt_pct_tolerance, 0.5);
    assert_eq!(settings.rt_tolerance(), 0.005);

    settings.rt_pct_tolerance = 1.25;
    assert_eq!(settings.rt_tolerance(), 0.0125);
    assert!(settings.mbr);
}

#[test]
fn disabling_mbr_keeps_one_anchor_per_identified_file() {
    let parameters: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false
    }))
    .unwrap();
    let parameters = parameters.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDE\n");
    let db = parameters.build_from_peptides(peptides);
    let features = [0, 1].map(|file_id| Feature {
        peptide_idx: crate::database::PeptideIx(0),
        peptide_q: 0.0,
        label: 1,
        file_id,
        aligned_rt: 0.5,
        charge: 2,
        ..Feature::default()
    });

    let with_mbr = build_feature_map(LfqSettings::default(), (2, 2), &features, &db);
    let without_mbr = build_feature_map(
        LfqSettings {
            mbr: false,
            ..LfqSettings::default()
        },
        (2, 2),
        &features,
        &db,
    );
    assert_eq!(with_mbr.ranges.len(), 6);
    assert_eq!(without_mbr.ranges.len(), 12);
    assert_eq!(
        without_mbr
            .ranges
            .iter()
            .map(|range| range.file_id)
            .collect::<std::collections::HashSet<_>>(),
        std::collections::HashSet::from([0, 1])
    );
}

#[test]
fn one_identified_label_channel_seeds_all_channel_precursors() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "static_mods": {
            "Arg10": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269},
                "sites": ["R"]
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    parameters.validate_channels().unwrap();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDER\n");
    let db = parameters.build_from_peptides(peptides);
    let light = db
        .peptides
        .iter()
        .position(|peptide| peptide.label_channel.as_deref() == Some("light"))
        .unwrap();
    let feature = Feature {
        peptide_idx: crate::database::PeptideIx(light as u32),
        peptide_q: 0.0,
        label: 1,
        aligned_rt: 0.5,
        ims: 1.0,
        calcmass: db.peptides[light].monoisotopic,
        charge: 2,
        ..Feature::default()
    };

    let map = build_feature_map(LfqSettings::default(), (2, 2), &[feature], &db);
    let seeded = map
        .ranges
        .iter()
        .filter(|range| !range.decoy && range.isotope == 0)
        .map(|range| range.peptide)
        .collect::<std::collections::HashSet<_>>();
    assert_eq!(seeded.len(), 2);
}

#[test]
fn shared_light_variable_channel_seeds_every_heavy_site_pattern() {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "generate_decoys": false,
        "max_variable_mods": 2,
        "variable_mods": {
            "Optional-Lys8": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 8.014199},
                "sites": ["K"]
            }
        }
    }))
    .unwrap();
    let parameters = builder.make_parameters();
    parameters.validate_channels().unwrap();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDEKK\n");
    let db = parameters.build_from_peptides(peptides);
    let light = db
        .peptides
        .iter()
        .position(|peptide| peptide.label_channel.as_deref() == Some("light"))
        .unwrap();
    let feature = Feature {
        peptide_idx: crate::database::PeptideIx(light as u32),
        peptide_q: 0.0,
        label: 1,
        aligned_rt: 0.5,
        ims: 1.0,
        calcmass: db.peptides[light].monoisotopic,
        charge: 2,
        ..Feature::default()
    };

    let map = build_feature_map(LfqSettings::default(), (2, 2), &[feature], &db);
    let seeded = map
        .ranges
        .iter()
        .filter(|range| !range.decoy && range.isotope == 0)
        .map(|range| range.peptide)
        .collect::<std::collections::HashSet<_>>();
    assert_eq!(seeded.len(), 4);
}

#[test]
fn gaussian_kernel_is_symmetric_and_normalized() {
    let kernel = gaussian_kernel(0.5, 11);

    assert!((kernel.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!(kernel
        .iter()
        .zip(kernel.iter().rev())
        .all(|(left, right)| (left - right).abs() < 1e-12));
    assert_eq!(
        kernel.iter().copied().max_by(f64::total_cmp),
        Some(kernel[5])
    );
}

#[test]
fn convolution_matches_same_mode_at_signal_boundaries() {
    let convolved = convolve(&[1.0, 2.0, 3.0], &[0.25, 0.5, 0.25]);

    assert_eq!(convolved, vec![1.0, 2.0, 2.0]);
}

fn precursor(rt: f32, mass_lo: f32, mass_hi: f32, mobility: (f32, f32)) -> PrecursorRange {
    PrecursorRange {
        rt,
        mass_lo,
        mass_hi,
        mobility_lo: mobility.0,
        mobility_hi: mobility.1,
        charge: 2,
        isotope: 0,
        peptide: PeptideIx(0),
        file_id: 0,
        decoy: false,
    }
}

#[test]
fn query_filters_by_mass_retention_time_and_mobility() {
    let ranges = vec![
        precursor(10.0, 499.9, 500.1, (0.9, 1.1)),
        precursor(20.0, 499.9, 500.1, (0.9, 1.1)),
        precursor(10.0, 599.9, 600.1, (0.9, 1.1)),
    ];
    let query = Query {
        ranges: &ranges,
        page_lo: 0,
        page_hi: 1,
        bin_size: ranges.len(),
        min_rt: 9.0,
        max_rt: 11.0,
        mass_search_margin: 0.2,
    };

    assert_eq!(query.mass_lookup(500.0).count(), 1);
    assert_eq!(query.mass_mobility_lookup(500.0, 1.0).count(), 1);
    assert_eq!(query.mass_mobility_lookup(500.0, 0.9).count(), 1);
    assert_eq!(query.mass_mobility_lookup(500.0, 1.1).count(), 1);
    assert_eq!(query.mass_mobility_lookup(500.0, 0.899).count(), 0);
    assert_eq!(query.mass_mobility_lookup(500.0, 1.2).count(), 0);
    assert_eq!(query.mass_lookup(700.0).count(), 0);
}

#[test]
fn query_uses_the_configured_mass_range_instead_of_a_fixed_margin() {
    let ranges = vec![precursor(10.0, 499.75, 500.25, (0.9, 1.1))];
    let query = Query {
        ranges: &ranges,
        page_lo: 0,
        page_hi: 1,
        bin_size: ranges.len(),
        min_rt: 9.0,
        max_rt: 11.0,
        mass_search_margin: 0.5,
    };

    assert_eq!(query.mass_lookup(500.2).count(), 1);
    assert_eq!(query.mass_lookup(500.3).count(), 0);
}

#[test]
fn grid_interpolation_conserves_intensity() {
    let entry = precursor(10.0, 499.9, 500.1, (0.9, 1.1));
    let mut grid = Grid::new(&entry, 1.0, [1.0, 0.0, 0.0], 2, 10);

    grid.add_entry(9.1, 0, 0, 100.0);
    grid.add_entry(10.0, 1, 1, 50.0);

    assert!((grid.matrix.row_slice(0).iter().sum::<f64>() - 100.0).abs() < 1e-6);
    assert!((grid.matrix.row_slice(4).iter().sum::<f64>() - 50.0).abs() < 1e-6);
    assert_eq!(grid.reference_file_id, 0);
}

fn traces() -> Traces {
    Traces {
        dot_product: Matrix::new(
            [
                0.0, 1.0, 6.0, 10.0, 6.0, 1.0, 0.0, 0.0, 1.0, 6.0, 10.0, 6.0, 1.0, 0.0,
            ],
            2,
            7,
        ),
        spectral_angle: Matrix::new(
            [
                0.0, 0.4, 0.9, 1.0, 0.9, 0.4, 0.0, 0.0, 0.4, 0.9, 1.0, 0.9, 0.4, 0.0,
            ],
            2,
            7,
        ),
        reference_file_id: 0,
    }
}

#[test]
fn time_warp_finds_and_applies_a_shifted_trace() {
    let matrix = Matrix::new([0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0], 2, 5);
    let trace = Traces {
        dot_product: matrix.clone(),
        spectral_angle: matrix,
        reference_file_id: 0,
    };
    let warps = trace.find_time_warps(&trace.dot_product, 2);

    assert_eq!(warps, vec![0, 1]);
    let mut shifted = trace.dot_product.clone();
    Traces::apply_time_warps(&mut shifted, &warps);
    assert_eq!(shifted.row_slice(1), &[0.0, 1.0, 2.0, 1.0, 0.0]);
}

#[test]
fn every_peak_scoring_strategy_prefers_the_centered_matching_peak() {
    let trace = traces();

    for strategy in [
        PeakScoringStrategy::RetentionTime,
        PeakScoringStrategy::SpectralAngle,
        PeakScoringStrategy::Intensity,
        PeakScoringStrategy::Hybrid,
    ] {
        let (scores, spectral) = trace.scores(strategy);
        assert_eq!(
            scores
                .iter()
                .copied()
                .enumerate()
                .max_by(|left, right| left.1.total_cmp(&right.1))
                .unwrap()
                .0,
            3
        );
        assert!(spectral[3] > spectral[2]);
        assert!(scores.iter().all(|score| score.is_finite()));
    }
}

#[test]
fn integration_supports_sum_and_apex_and_rejects_weak_matches() {
    let mut sum_trace = traces();
    let mut settings = LfqSettings {
        spectral_angle: 0.5,
        integration: IntegrationStrategy::Sum,
        ..Default::default()
    };
    let (sum_peak, sum_areas) = sum_trace.integrate(&settings).unwrap();

    let mut apex_trace = traces();
    settings.integration = IntegrationStrategy::Apex;
    let (apex_peak, apex_areas) = apex_trace.integrate(&settings).unwrap();

    assert_eq!(sum_peak.rt, 3);
    assert_eq!(apex_peak.rt, 3);
    assert!(sum_areas[0].unwrap() > apex_areas[0].unwrap());
    assert_eq!(apex_areas, vec![Some(10.0), Some(10.0)]);

    let mut weak = traces();
    settings.spectral_angle = 1.1;
    assert!(weak.integrate(&settings).is_none());
}

#[test]
fn summarized_isotope_traces_reward_theoretical_abundance() {
    let entry = precursor(10.0, 499.9, 500.1, (0.9, 1.1));
    let distribution = [0.8, 0.15, 0.05];
    let mut grid = Grid::new(&entry, 1.0, distribution, 1, 21);
    for (isotope, abundance) in distribution.into_iter().enumerate() {
        grid.add_entry(10.0, isotope, 0, abundance * 1000.0);
    }

    let traces = grid.summarize_traces();
    let center = traces.spectral_angle.cols / 2;

    assert!(traces.dot_product[(0, center)] > 0.0);
    assert!(traces.spectral_angle[(0, center)] > 0.99);
    assert!(traces.spectral_angle[(0, 0)] <= traces.spectral_angle[(0, center)]);
}

#[test]
fn file_evidence_exposes_a_weak_signal_despite_a_strong_shared_peak() {
    let mut trace = traces();
    trace.spectral_angle.row_slice_mut(1).fill(0.1);
    for value in trace.dot_product.row_slice_mut(1) {
        *value *= 0.01;
    }
    let (_, areas, evidence) = trace
        .integrate_with_evidence(&LfqSettings::default())
        .unwrap();
    assert!(areas.iter().all(Option::is_some));
    assert!(evidence[0].as_ref().unwrap().score > 0.9);
    assert!(evidence[1].as_ref().unwrap().score < 0.01);
}

#[test]
fn file_evidence_does_not_invent_a_score_for_a_missing_trace() {
    let mut trace = traces();
    trace.dot_product.row_slice_mut(1).fill(0.0);
    let (_, areas, evidence) = trace
        .integrate_with_evidence(&LfqSettings::default())
        .unwrap();
    assert!(areas[1].is_none());
    assert!(evidence[1].is_none());
}

#[test]
fn strict_ms2_evidence_requires_both_psm_and_peptide_acceptance() {
    let builder: Builder =
        serde_json::from_value(serde_json::json!({"generate_decoys": false})).unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDE\n");
    let db = parameters.build_from_peptides(peptides);
    let features = [0.001, 0.5, f32::NAN]
        .into_iter()
        .enumerate()
        .map(|(file_id, spectrum_q)| Feature {
            peptide_idx: PeptideIx(0),
            peptide_q: 0.001,
            spectrum_q,
            label: 1,
            file_id,
            aligned_rt: 0.5,
            charge: 2,
            ..Default::default()
        })
        .collect::<Vec<_>>();
    let map = build_feature_map(LfqSettings::default(), (2, 2), &features, &db);
    let id = PrecursorId::Combined(PeptideIx(0));
    assert!(map.ms2_confirmed.contains(&(id, 1)));
    assert!(map.ms2_confirmed_strict.contains(&(id, 0)));
    assert!(!map.ms2_confirmed_strict.contains(&(id, 1)));
    assert!(!map.ms2_confirmed_strict.contains(&(id, 2)));
}

fn identity_alignment(file_id: usize) -> Alignment {
    Alignment {
        file_id,
        max_rt: 1.0,
        slope: 1.0,
        intercept: 0.0,
        knots: Vec::new(),
    }
}

/// Synthesize MS1 spectra with an identical isotope envelope eluting at the
/// given `(file_id, rt)` apexes.
fn eluting_ms1_spectra(
    map: &FeatureMap,
    db: &IndexedDatabase,
    apexes: &[(usize, f32)],
) -> Vec<ProcessedSpectrum> {
    let composition = db.peptides[0]
        .sequence
        .iter()
        .map(|r| composition(*r))
        .sum::<Composition>();
    let dist = crate::isotopes::peptide_isotopes(composition.carbon, composition.sulfur);
    let mut masses = (0..N_ISOTOPES)
        .map(|isotope| {
            let range = map
                .ranges
                .iter()
                .find(|range| !range.decoy && range.isotope == isotope)
                .unwrap();
            (isotope, (range.mass_lo + range.mass_hi) / 2.0)
        })
        .collect::<Vec<_>>();
    masses.sort_by(|a, b| a.1.total_cmp(&b.1));

    let rt_tol = map.settings.rt_tolerance();
    apexes
        .iter()
        .flat_map(|&(file_id, anchor)| {
            let masses = masses.clone();
            (-40..=40).map(move |step| {
                let rt = anchor + step as f32 * rt_tol / 50.0;
                let apex = (-0.5 * (step as f32 / 10.0).powi(2)).exp() * 1000.0;
                ProcessedSpectrum {
                    level: 1,
                    file_id,
                    scan_start_time: rt,
                    masses: masses.iter().map(|(_, mass)| *mass).collect(),
                    intensities: masses
                        .iter()
                        // Slightly off the theoretical envelope: a perfect
                        // match rounds the spectral similarity above 1.0
                        .map(|(isotope, _)| dist[*isotope] * apex * (1.0 - 0.05 * *isotope as f32))
                        .collect(),
                    ..Default::default()
                }
            })
        })
        .collect()
}

#[test]
fn disabling_mbr_quantifies_each_file_against_its_own_anchor() {
    let builder: Builder =
        serde_json::from_value(serde_json::json!({"generate_decoys": false})).unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDE\n");
    let db = parameters.build_from_peptides(peptides);
    let features = [(0, 0.2), (1, 0.6)].map(|(file_id, aligned_rt)| Feature {
        peptide_idx: PeptideIx(0),
        peptide_q: 0.0,
        spectrum_q: 0.0,
        label: 1,
        file_id,
        aligned_rt,
        charge: 2,
        ..Feature::default()
    });
    let settings = LfqSettings {
        mbr: false,
        ..LfqSettings::default()
    };
    let map = build_feature_map(settings, (2, 2), &features, &db);
    let spectra = eluting_ms1_spectra(&map, &db, &[(0, 0.2), (1, 0.6)]);
    let alignments = [identity_alignment(0), identity_alignment(1)];

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let key = (PrecursorId::Combined(PeptideIx(0)), false);
    let first = pool.install(|| map.quantify(&db, &spectra, &alignments));
    let peak = &first[&key];
    let areas = peak
        .intensities
        .iter()
        .map(|area| area.expect("each identified file has a quantified peak"))
        .collect::<Vec<_>>();
    assert!(
        (areas[0] - areas[1]).abs() <= areas[0] * 1e-4,
        "identical files quantified differently: {areas:?}"
    );
    assert_eq!(peak.ms2_confirmed, vec![true, true]);
    for evidence in &peak.file_evidence {
        let evidence = evidence.as_ref().unwrap();
        assert!(!evidence.transfer_candidate);
        assert_eq!(evidence.rt_shift_bins, 0);
    }

    for _ in 0..20 {
        let repeat = pool.install(|| map.quantify(&db, &spectra, &alignments));
        let repeat = &repeat[&key];
        assert_eq!(repeat.intensities, peak.intensities);
        assert_eq!(repeat.peak.score.to_bits(), peak.peak.score.to_bits());
        assert_eq!(repeat.peak.rt, peak.peak.rt);
    }
}

#[test]
fn mbr_traces_one_anchor_across_files_deterministically() {
    let builder: Builder =
        serde_json::from_value(serde_json::json!({"generate_decoys": false})).unwrap();
    let parameters = builder.make_parameters();
    let peptides = parameters.peptides_from_tsv("sequence\nPEPTIDE\n");
    let db = parameters.build_from_peptides(peptides);
    let feature = Feature {
        peptide_idx: PeptideIx(0),
        peptide_q: 0.0,
        spectrum_q: 0.0,
        label: 1,
        file_id: 0,
        aligned_rt: 0.4,
        charge: 2,
        ..Feature::default()
    };
    let map = build_feature_map(LfqSettings::default(), (2, 2), &[feature], &db);
    let spectra = eluting_ms1_spectra(&map, &db, &[(0, 0.4), (1, 0.4), (2, 0.4)]);
    let alignments = [0, 1, 2].map(identity_alignment);

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let key = (PrecursorId::Combined(PeptideIx(0)), false);
    let first = pool.install(|| map.quantify(&db, &spectra, &alignments));
    let peak = &first[&key];
    let reference = peak.intensities[0].unwrap();
    assert!(peak
        .intensities
        .iter()
        .all(|area| (area.unwrap() - reference).abs() <= reference * 1e-4));
    assert_eq!(peak.ms2_confirmed, vec![true, false, false]);
    let transfers = peak
        .file_evidence
        .iter()
        .map(|evidence| evidence.as_ref().unwrap().transfer_candidate)
        .collect::<Vec<_>>();
    assert_eq!(transfers, vec![false, true, true]);

    for _ in 0..20 {
        let repeat = pool.install(|| map.quantify(&db, &spectra, &alignments));
        let repeat = &repeat[&key];
        assert_eq!(repeat.intensities, peak.intensities);
        assert_eq!(repeat.peak.score.to_bits(), peak.peak.score.to_bits());
    }
}
