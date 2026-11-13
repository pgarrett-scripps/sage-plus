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

fn at(ms_level: u8, scan_start_time: f32, total_ion_current: f32) -> RawSpectrum {
    RawSpectrum {
        ms_level,
        scan_start_time,
        total_ion_current,
        ..Default::default()
    }
}

#[test]
fn census_bins_the_gradient_into_ten_equal_widths() {
    // MS2 every minute from 0 to 9 inclusive, so the span is 0..=9 and each
    // bin is 0.9 minutes wide: bin 0 takes 0.0 and 0.9 would start bin 1.
    let mut spectra = (0..10).map(|n| at(2, n as f32, 0.0)).collect::<Vec<_>>();
    spectra.push(at(1, 0.0, 500.0));
    spectra.push(at(1, 9.0, 1500.0));
    let census = census(&spectra);
    assert_eq!(census.ms2_spectra, 10);
    assert_eq!(census.ms1_spectra, 2);
    assert!(census.ms1_read);
    assert_eq!((census.rt_start, census.rt_end), (0.0, 9.0));
    assert_eq!(census.bins.len(), 10);
    let acquired = census
        .bins
        .iter()
        .map(|bin| bin.ms2_acquired)
        .collect::<Vec<_>>();
    // Ten spectra over a nine-minute span, bins 0.9 minutes wide: minute 8
    // lands in bin 8 (8.0 / 0.9 = 8.9) and minute 9, the span's end, in the
    // last bin, so each bin holds one.
    assert_eq!(acquired, vec![1; 10]);
    assert_eq!(
        census
            .bins
            .iter()
            .map(|bin| bin.ms2_acquired)
            .sum::<usize>(),
        10
    );
    assert_eq!(census.bins[0].ms1_tic, 500.0);
    assert_eq!(census.bins[9].ms1_tic, 1500.0);
}

#[test]
fn census_without_ms1_measures_nothing_about_ms1() {
    let census = census(&[at(2, 1.0, 0.0), at(2, 2.0, 0.0)]);
    assert_eq!(census.ms1_spectra, 0);
    assert!(!census.ms1_read);
    assert!(census.bins.iter().all(|bin| bin.ms1_tic == 0.0));
}

#[test]
fn census_of_no_spectra_has_no_bins() {
    let census = census(&[]);
    assert_eq!(census, SpectrumCensus::default());
    assert!(census.bins.is_empty());
}

#[test]
fn census_tic_falls_back_to_the_summed_peaks() {
    // An mzML without the `total ion current` cvParam states none, so the
    // peaks are summed instead.
    let stated = census(&[at(1, 0.0, 42.0)]);
    let summed = census(&[RawSpectrum {
        ms_level: 1,
        intensity: vec![20.0, 22.0],
        ..Default::default()
    }]);
    assert_eq!(stated.bins[0].ms1_tic, 42.0);
    assert_eq!(summed.bins[0].ms1_tic, 42.0);
}

#[test]
fn a_file_spanning_no_time_puts_everything_in_the_first_bin() {
    let census = census(&[at(2, 5.0, 0.0), at(2, 5.0, 0.0)]);
    assert_eq!((census.rt_start, census.rt_end), (5.0, 5.0));
    assert_eq!(census.bins[0].ms2_acquired, 2);
    assert!(census.bins[1..].iter().all(|bin| bin.ms2_acquired == 0));
}

#[test]
fn gradient_bin_is_inclusive_of_the_last_edge() {
    assert_eq!(gradient_bin(0.0, 0.0, 10.0), 0);
    assert_eq!(gradient_bin(0.999, 0.0, 10.0), 0);
    assert_eq!(gradient_bin(1.0, 0.0, 10.0), 1);
    assert_eq!(gradient_bin(10.0, 0.0, 10.0), 9);
    // Outside the span, and a retention time no file states.
    assert_eq!(gradient_bin(-1.0, 0.0, 10.0), 0);
    assert_eq!(gradient_bin(11.0, 0.0, 10.0), 9);
    assert_eq!(gradient_bin(f32::NAN, 0.0, 10.0), 0);
}

/// Identification stats of one file, built by hand so the serializers can be
/// tested without a search.
fn identification_fixture() -> IdentificationFileStats {
    let bins = (0..10)
        .map(|bin| GradientBinStats {
            bin,
            rt_start: bin as f32,
            rt_end: bin as f32 + 1.0,
            ms2_acquired: 10,
            // The middle of the gradient identifies nothing: a dropout.
            ms2_identified: if (4..6).contains(&bin) { 0 } else { 5 },
            ms1_spectra: 2,
            ms1_tic: 1000.0,
        })
        .collect::<Vec<_>>();
    IdentificationFileStats {
        file: "a.mzML".into(),
        census: SpectrumCensus {
            ms1_spectra: 20,
            ms2_spectra: 100,
            rt_start: 0.0,
            rt_end: 10.0,
            ms1_read: true,
            bins: vec![GradientBinCensus::default(); 10],
        },
        psms: 40,
        ms2_identified: 40,
        peptides: 30,
        proteins: 12,
        ms1_tic: 10_000.0,
        bins,
    }
}

#[test]
fn identification_rates_are_shares_of_acquired_spectra() {
    let stats = identification_fixture();
    assert_eq!(stats.identified_pct(), 40.0);
    assert_eq!(stats.empty_bins(), 2);
    assert_eq!(stats.bin_pct_range(), Some((0.0, 50.0)));
    // A bin that acquired nothing is not a dropout and not a rate.
    let mut quiet = stats.clone();
    quiet.bins.iter_mut().for_each(|bin| {
        bin.ms2_acquired = 0;
        bin.ms2_identified = 0;
    });
    assert_eq!(quiet.empty_bins(), 0);
    assert_eq!(quiet.bin_pct_range(), None);
    assert_eq!(quiet.bins[0].identified_pct(), 0.0);
}

#[test]
fn qc_tables_have_the_published_columns_and_one_row_per_file() {
    let stats = vec![identification_fixture()];
    let qc = String::from_utf8(
        serialize_qc(
            &stats,
            &DigestionRunStats::default(),
            &[],
            &CompositionProfiles::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let lines = qc.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2, "one header and one file row: {qc}");
    let header = lines[0].split('\t').collect::<Vec<_>>();
    let row = lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), row.len());
    let column = |name: &str| row[header.iter().position(|&h| h == name).expect(name)];
    assert_eq!(column("file"), "a.mzML");
    assert_eq!(column("ms2_spectra"), "100");
    assert_eq!(column("ms2_identified_pct"), "40.000");
    assert_eq!(column("peptides"), "30");
    assert_eq!(column("proteins"), "12");
    assert_eq!(column("ms1_read"), "1");
    assert_eq!(column("empty_gradient_bins"), "2");
    assert_eq!(column("min_bin_identified_pct"), "0.000");
    assert_eq!(column("max_bin_identified_pct"), "50.000");
    // No digestion row and no polymer row for this file: the digestion
    // columns are zero and the polymer column says so.
    assert_eq!(column("digestion_peptides"), "0");
    assert_eq!(column("top_polymer"), "-");
    // No composition profile for this file: an absent shift is empty, not a
    // shift of zero, as in qc.composition.tsv.
    assert_eq!(column("largest_shift_residue"), "");
    assert_eq!(column("largest_shift_log2_ratio"), "");
    assert_eq!(column("jensen_shannon_distance"), "");

    let gradient = String::from_utf8(serialize_qc_gradient(&stats).unwrap()).unwrap();
    let lines = gradient.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 11, "one header and ten bins: {gradient}");
    assert_eq!(
        lines[1].split('\t').collect::<Vec<_>>(),
        vec!["a.mzML", "0", "0.0000", "1.0000", "10", "5", "50.000", "2", "1000.0"]
    );
    assert_eq!(lines[5].split('\t').nth(6), Some("0.000"));
}

#[test]
fn qc_table_columns_match_the_published_schema() {
    for (schema, header) in [
        (
            "../../schemas/qc.v1.tsv.schema.json",
            String::from_utf8(
                serialize_qc(
                    &[],
                    &DigestionRunStats::default(),
                    &[],
                    &CompositionProfiles::default(),
                )
                .unwrap(),
            )
            .unwrap(),
        ),
        (
            "../../schemas/qc_gradient.v1.tsv.schema.json",
            String::from_utf8(serialize_qc_gradient(&[]).unwrap()).unwrap(),
        ),
    ] {
        let table: serde_json::Value =
            serde_json::from_slice(&std::fs::read(schema).expect(schema)).expect(schema);
        let columns = table["fields"]
            .as_array()
            .expect("fields")
            .iter()
            .map(|field| field["name"].as_str().expect("name"))
            .collect::<Vec<_>>();
        assert_eq!(
            header
                .lines()
                .next()
                .expect("header")
                .split('\t')
                .collect::<Vec<_>>(),
            columns,
            "{schema}"
        );
    }
}

#[test]
fn qc_table_carries_the_composition_headline_of_each_file() {
    // `profiles()` gives "alkylated.mzML" half the cysteine the database
    // predicts and "control.mzML" a profile that matches it exactly.
    let composition = profiles();
    let stats = composition
        .files
        .iter()
        .map(|profile| IdentificationFileStats {
            file: profile.file.clone(),
            ..identification_fixture()
        })
        .collect::<Vec<_>>();
    let qc = String::from_utf8(
        serialize_qc(&stats, &DigestionRunStats::default(), &[], &composition).unwrap(),
    )
    .unwrap();
    let lines = qc.lines().collect::<Vec<_>>();
    let header = lines[0].split('\t').collect::<Vec<_>>();
    let cell = |row: usize, name: &str| {
        lines[row].split('\t').collect::<Vec<_>>()
            [header.iter().position(|&h| h == name).expect(name)]
        .to_string()
    };
    assert_eq!(cell(1, "file"), "alkylated.mzML");
    assert_eq!(cell(1, "largest_shift_residue"), "C");
    assert_eq!(cell(1, "largest_shift_log2_ratio"), "-1.0000");
    assert_ne!(cell(1, "jensen_shannon_distance"), "");
    // A profile that matches the database has no shift and no distance.
    assert_eq!(cell(2, "file"), "control.mzML");
    assert_eq!(cell(2, "largest_shift_log2_ratio"), "0.0000");
    assert_eq!(cell(2, "jensen_shannon_distance"), "0.0000");
}
