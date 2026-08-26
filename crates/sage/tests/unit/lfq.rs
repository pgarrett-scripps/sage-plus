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
            "R": {
                "mass": 0.0,
                "channel_offsets": {"light": 0.0, "heavy": 10.008269}
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
            "K": [{
                "mass": 0.0,
                "name": "Optional-Lys8",
                "channel_offsets": {"light": 0.0, "heavy": 8.014199}
            }]
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
