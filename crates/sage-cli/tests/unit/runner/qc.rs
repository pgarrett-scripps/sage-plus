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
