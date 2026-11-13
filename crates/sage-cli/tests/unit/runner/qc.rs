use super::*;

fn ms1(mz: Vec<f32>, intensity: Vec<f32>, representation: Representation) -> RawSpectrum {
    RawSpectrum {
        ms_level: 1,
        representation,
        mz,
        intensity,
        ..Default::default()
    }
}

#[test]
fn polymer_scan_counts_centroided_ms1_only() {
    // Singly protonated PEG n = 10..=15, listed out of order, plus one
    // unrelated peak of the same total intensity.
    let mut mz = (10..=15)
        .rev()
        .map(|n| (18.010_565 + n as f64 * 44.026_215 + 1.007_276) as f32)
        .collect::<Vec<_>>();
    mz.push(300.0);
    let mut intensity = vec![1.0; 6];
    intensity.push(6.0);
    let spectra = vec![
        ms1(mz.clone(), intensity.clone(), Representation::Centroid),
        ms1(mz, intensity, Representation::Profile),
        RawSpectrum {
            ms_level: 2,
            ..Default::default()
        },
    ];
    let stats = scan_polymers(&spectra).expect("MS1 present");
    assert_eq!(stats.ms1_spectra, 1);
    assert_eq!(stats.skipped_profile_spectra, 1);
    let peg = &stats.polymers[0];
    assert_eq!(peg.name, "peg");
    assert!((peg.tic_pct - 50.0).abs() < 1e-9, "{peg:?}");
    assert!(stats.polymers[1..].iter().all(|share| share.tic_pct == 0.0));
}

#[test]
fn polymer_scan_needs_ms1() {
    let spectra = vec![RawSpectrum {
        ms_level: 2,
        ..Default::default()
    }];
    assert!(scan_polymers(&spectra).is_none());
}

fn residue_counts(sequences: &[&str]) -> CompositionProfile {
    use sage_core::digestion::TerminusClass;
    let mut profile = CompositionProfile::default();
    for sequence in sequences {
        // One missed-cleavage class, so the weighting cannot move a fraction.
        profile.add_sequence(sequence.as_bytes(), Some(TerminusClass::NonEnzymatic), 0);
    }
    profile
}

fn profiles() -> CompositionProfiles {
    // One file sees half the cysteine the database predicts; the other matches.
    let short = CompositionSummary::new(&residue_counts(&["CAAA"]), &residue_counts(&["CCAA"]));
    let exact = CompositionSummary::new(&residue_counts(&["CAAA"]), &residue_counts(&["CAAA"]));
    CompositionProfiles {
        q_value: 0.01,
        files: vec![
            CompositionFileProfile {
                file: "alkylated.mzML".into(),
                summary: short,
            },
            CompositionFileProfile {
                file: "control.mzML".into(),
                summary: exact.clone(),
            },
        ],
        total: exact,
    }
}

#[test]
fn composition_tsv_has_one_row_per_file_and_residue() {
    let bytes = serialize_composition(&profiles()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "file\tresidue\tobserved\tobserved_fraction\texpected\texpected_fraction\tlog2_ratio"
    );
    let rows = lines.collect::<Vec<_>>();
    // Two residues for each of two files, then the total rows.
    assert_eq!(rows.len(), 6);
    assert_eq!(
        rows[0],
        "alkylated.mzML\tA\t3\t0.750000\t2\t0.500000\t0.5850"
    );
    assert_eq!(
        rows[1],
        "alkylated.mzML\tC\t1\t0.250000\t2\t0.500000\t-1.0000"
    );
    assert!(rows[4].starts_with("total\tA\t"));
    assert!(rows[5].starts_with("total\tC\t"));
}

#[test]
fn composition_tsv_leaves_an_undefined_ratio_empty() {
    let summary = CompositionSummary::new(&residue_counts(&["AAAA"]), &residue_counts(&["AACC"]));
    let profiles = CompositionProfiles {
        q_value: 0.01,
        files: Vec::new(),
        total: summary,
    };
    let text = String::from_utf8(serialize_composition(&profiles).unwrap()).unwrap();
    let carbon = text
        .lines()
        .find(|line| line.starts_with("total\tC\t"))
        .unwrap();
    // Observed none of it, so the ratio is written as nothing at all.
    assert_eq!(carbon, "total\tC\t0\t0.000000\t2\t0.500000\t");
}

#[test]
fn composition_run_stats_carry_the_headline_numbers_only() {
    let stats = profiles().run_stats();
    assert_eq!(stats.q_value, 0.01);
    assert_eq!(stats.files.len(), 2);
    let flagged = &stats.files[0];
    assert_eq!(flagged.file, "alkylated.mzML");
    assert_eq!(flagged.observed_peptides, 1);
    assert_eq!(flagged.observed_residues, 4);
    assert_eq!(flagged.matched_observed_residues, 4);
    assert_eq!(flagged.expected_peptides, 1);
    assert_eq!(flagged.expected_residues, 4);
    assert_eq!(flagged.largest_shift_residue, Some('C'));
    assert_eq!(flagged.largest_shift_log2_ratio, Some(-1.0));
    assert!(flagged.jensen_shannon_distance.unwrap() > 0.0);
    // The file that matches its database has nothing to report.
    assert_eq!(stats.files[1].largest_shift_log2_ratio, Some(0.0));
    assert_eq!(stats.files[1].jensen_shannon_distance, Some(0.0));
    assert_eq!(stats.total.file, "total");
}

#[test]
fn composition_run_stats_omit_absent_numbers_from_json() {
    let profiles = CompositionProfiles {
        q_value: 0.01,
        files: vec![CompositionFileProfile {
            file: "empty.mzML".into(),
            summary: CompositionSummary::new(
                &CompositionProfile::default(),
                &residue_counts(&["AC"]),
            ),
        }],
        total: CompositionSummary::default(),
    };
    let json = serde_json::to_string(&profiles.run_stats()).unwrap();
    // A file with no identifications states its zeros and claims no distance.
    assert!(json.contains("\"observed_peptides\":0"), "{json}");
    assert!(!json.contains("jensen_shannon_distance"), "{json}");
    assert!(!json.contains("largest_shift_residue"), "{json}");
}

fn occupancy(
    site_rule: &str,
    modified_sites: usize,
    eligible_sites: usize,
) -> ModificationOccupancy {
    ModificationOccupancy {
        site_rule: site_rule.into(),
        mass: 15.9949,
        name: None,
        search_mode: "database".into(),
        peptidoforms: 3,
        modified_peptidoforms: 1,
        eligible_sites,
        modified_sites,
        occupancy_pct: 100.0 * modified_sites as f64 / eligible_sites as f64,
    }
}

#[test]
fn modification_occupancy_tsv_has_a_row_per_file_and_rule_then_the_total() {
    let stats = ModificationOccupancyRunStats {
        q_value: 0.01,
        modifications: vec![occupancy("M", 3, 8)],
        files: vec![
            ModificationOccupancyFileStats {
                file: "a.mzML".into(),
                modifications: vec![occupancy("M", 1, 4)],
            },
            ModificationOccupancyFileStats {
                file: "b.mzML".into(),
                modifications: vec![occupancy("M", 2, 4)],
            },
        ],
    };
    let tsv = String::from_utf8(serialize_modification_occupancy(&stats).unwrap()).unwrap();
    let rows = tsv.lines().collect::<Vec<_>>();
    assert_eq!(
        rows[0],
        "file\tsite_rule\tmass\tname\tsearch_mode\tpeptidoforms\tmodified_peptidoforms\teligible_sites\tmodified_sites\toccupancy_pct"
    );
    assert_eq!(
        rows[1],
        "a.mzML\tM\t15.994900\t\tdatabase\t3\t1\t4\t1\t25.000"
    );
    assert_eq!(
        rows[2],
        "b.mzML\tM\t15.994900\t\tdatabase\t3\t1\t4\t2\t50.000"
    );
    assert_eq!(
        rows[3],
        "total\tM\t15.994900\t\tdatabase\t3\t1\t8\t3\t37.500"
    );
    assert_eq!(rows.len(), 4);
}
