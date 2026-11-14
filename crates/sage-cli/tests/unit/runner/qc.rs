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

fn chromatography(charges: &[u8]) -> ChromatographyStats {
    let mut stats = ChromatographyStats::default();
    for &charge in charges {
        stats.add_charge(charge);
    }
    stats
}

#[test]
fn chromatography_tsv_labels_which_width_each_row_carries() {
    let mut quantified = chromatography(&[2, 2, 3, 5]);
    quantified.peak_width = Some(0.25);
    quantified.peak_width_kind = Some(PeakWidth::LfqMedianFwhm);
    quantified.quantified_precursors = Some(100);
    quantified.transferred_precursors = Some(25);
    let mut unquantified = chromatography(&[2, 0]);
    unquantified.peak_width = Some(0.5);
    unquantified.peak_width_kind = Some(PeakWidth::RepeatPsmSpread);
    let stats = ChromatographyRunStats {
        q_value: 0.01,
        files: vec![
            ChromatographyFileStats {
                file: "lfq.mzML".into(),
                stats: quantified,
            },
            ChromatographyFileStats {
                file: "plain.mzML".into(),
                stats: unquantified,
            },
        ],
        total: chromatography(&[2, 2, 3, 5, 2, 0]),
    };
    let text = String::from_utf8(serialize_chromatography(&stats).unwrap()).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "file\tconfident_precursors\tcharge_1_pct\tcharge_2_pct\tcharge_3_pct\tcharge_4_plus_pct\tcharge_unknown_pct\tpeak_width\tpeak_width_kind\tquantified_precursors\ttransferred_precursors\ttransferred_pct"
    );
    let rows = lines.collect::<Vec<_>>();
    assert_eq!(rows.len(), 3);
    assert_eq!(
        rows[0],
        "lfq.mzML\t4\t0.000\t50.000\t25.000\t25.000\t0.000\t0.2500\tlfq_median_fwhm\t100\t25\t25.000"
    );
    // Without quantification the three LFQ columns are empty, and the width is
    // named as the repeat-PSM spread it is.
    assert_eq!(
        rows[1],
        "plain.mzML\t2\t0.000\t50.000\t0.000\t0.000\t50.000\t0.5000\trepeat_psm_spread\t\t\t"
    );
    // The total pools the charges and leaves the width out.
    assert!(rows[2].starts_with("total\t6\t"));
    assert!(rows[2].ends_with("\t\t\t\t\t"), "{}", rows[2]);
}

#[test]
fn chromatography_run_stats_keep_the_kind_in_json() {
    let mut stats = chromatography(&[2]);
    stats.peak_width = Some(0.3);
    stats.peak_width_kind = Some(PeakWidth::RepeatPsmSpreadFallback);
    let json = serde_json::to_string(&ChromatographyRunStats {
        q_value: 0.01,
        files: vec![ChromatographyFileStats {
            file: "one.mzML".into(),
            stats,
        }],
        total: chromatography(&[2]),
    })
    .unwrap();
    assert!(
        json.contains("\"peak_width_kind\":\"repeat_psm_spread_fallback\""),
        "{json}"
    );
    assert!(json.contains("\"quantified_precursors\":null"), "{json}");
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
    // A bin that acquired nothing is not a dropout and not a rate: it is
    // counted on its own, so the two cannot be confused in qc.tsv.
    let mut quiet = stats.clone();
    quiet.bins.iter_mut().for_each(|bin| {
        bin.ms2_acquired = 0;
        bin.ms2_identified = 0;
    });
    assert_eq!(quiet.no_ms2_bins(), quiet.bins.len());
    assert_eq!(quiet.empty_bins(), 0);
    // Every bin of the fixture acquired MS2 spectra, so the two counts
    // separate the dropouts from the silence.
    assert_eq!(stats.no_ms2_bins(), 0);
    assert_eq!(quiet.bin_pct_range(), None);
    assert_eq!(quiet.bins[0].identified_pct(), 0.0);
}

/// `QcSources` with only the passes a test cares about filled in. Every other
/// pass is its own empty stats object, which is what a run without that
/// feature gives `qc.tsv`.
fn qc_sources<'a>(
    identification: &'a [IdentificationFileStats],
    digestion: &'a DigestionRunStats,
    composition: &'a CompositionProfiles,
    chromatography: &'a ChromatographyRunStats,
) -> QcSources<'a> {
    QcSources {
        identification,
        digestion,
        polymers: &[],
        composition,
        chromatography,
        modifications: None,
        contaminants: None,
        in_source: None,
    }
}

#[test]
fn qc_tables_have_the_published_columns_and_one_row_per_file() {
    let stats = vec![identification_fixture()];
    let qc = String::from_utf8(
        serialize_qc(&qc_sources(
            &stats,
            &DigestionRunStats::default(),
            &CompositionProfiles::default(),
            &ChromatographyRunStats::default(),
        ))
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
    assert_eq!(column("no_ms2_gradient_bins"), "0");
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
                serialize_qc(&qc_sources(
                    &[],
                    &DigestionRunStats::default(),
                    &CompositionProfiles::default(),
                    &ChromatographyRunStats::default(),
                ))
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
        serialize_qc(&qc_sources(
            &stats,
            &DigestionRunStats::default(),
            &composition,
            &ChromatographyRunStats::default(),
        ))
        .unwrap(),
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

#[test]
fn qc_table_repeats_the_headline_of_every_other_pass() {
    // Each column below is the number the pass that owns it already computed;
    // nothing in qc.tsv recomputes anything.
    let identification = vec![IdentificationFileStats {
        file: "a.mzML".into(),
        ..identification_fixture()
    }];
    let mut stats = chromatography(&[2, 2, 3, 4, 5]);
    stats.charge_unknown = 1;
    stats.confident_precursors += 1;
    stats.peak_width = Some(0.3125);
    stats.peak_width_kind = Some(PeakWidth::LfqMedianFwhm);
    stats.quantified_precursors = Some(40);
    stats.transferred_precursors = Some(10);
    let chromatography = ChromatographyRunStats {
        q_value: 0.01,
        files: vec![ChromatographyFileStats {
            file: "a.mzML".into(),
            stats,
        }],
        total: ChromatographyStats::default(),
    };
    let occupancy = |site_rule: &str, name: Option<&str>, pct: f64| ModificationOccupancy {
        site_rule: site_rule.into(),
        mass: 15.9949,
        name: name.map(str::to_string),
        search_mode: "database".into(),
        peptidoforms: 10,
        modified_peptidoforms: 4,
        eligible_sites: 20,
        modified_sites: (pct / 5.0) as usize,
        occupancy_pct: pct,
    };
    let modifications = ModificationOccupancyRunStats {
        q_value: 0.01,
        modifications: Vec::new(),
        files: vec![ModificationOccupancyFileStats {
            file: "a.mzML".into(),
            modifications: vec![
                occupancy("M", Some("Oxidation"), 42.5),
                occupancy("K", None, 7.0),
            ],
        }],
    };
    let contaminants = ContaminantRunStats {
        path: "contaminants.fasta".into(),
        proteins: 381,
        peptides: 5000,
        confident_psms: 100,
        contaminant_psms: 7,
        contaminant_psm_pct: 7.0,
        files: vec![ContaminantFileStats {
            file_id: 0,
            filename: "a.mzML".into(),
            confident_psms: 100,
            contaminant_psms: 7,
            contaminant_psm_pct: 7.0,
        }],
    };
    let qc = String::from_utf8(
        serialize_qc(&QcSources {
            identification: &identification,
            digestion: &DigestionRunStats::default(),
            polymers: &[],
            composition: &CompositionProfiles::default(),
            chromatography: &chromatography,
            modifications: Some(&modifications),
            contaminants: Some(&contaminants),
            in_source: None,
        })
        .unwrap(),
    )
    .unwrap();
    let lines = qc.lines().collect::<Vec<_>>();
    let header = lines[0].split('\t').collect::<Vec<_>>();
    let row = lines[1].split('\t').collect::<Vec<_>>();
    let column = |name: &str| row[header.iter().position(|&h| h == name).expect(name)];
    // Six precursors: two at 2+, one each at 3+ and 4+, one at 5+ counted as
    // 4 or above, and one whose charge the file did not state.
    assert_eq!(column("charge_1_pct"), "0.000");
    assert_eq!(column("charge_2_pct"), "33.333");
    assert_eq!(column("charge_3_pct"), "16.667");
    assert_eq!(column("charge_4_plus_pct"), "33.333");
    assert_eq!(column("charge_unknown_pct"), "16.667");
    assert_eq!(column("peak_width"), "0.3125");
    assert_eq!(column("peak_width_kind"), PeakWidth::LfqMedianFwhm.label());
    assert_eq!(column("quantified_precursors"), "40");
    assert_eq!(column("transferred_pct"), "25.000");
    // The most occupied definition, named as qc.modifications.tsv names it.
    assert_eq!(column("modifications_searched"), "2");
    assert_eq!(column("top_modification"), "Oxidation");
    assert_eq!(column("top_modification_occupancy_pct"), "42.500");
    // The contaminant share keeps the run summary's unit: a fraction.
    assert_eq!(column("contaminant_psms"), "7");
    // A percent, like every other share in the table.
    assert_eq!(column("contaminant_psm_pct"), "7.000");
    // in_source is off here, so the column is empty rather than zero.
    assert_eq!(column("in_source_flagged_psms"), "");
}

#[test]
fn a_definition_without_a_name_is_labelled_by_its_site_rule_and_mass() {
    let bare = ModificationOccupancy {
        site_rule: "K".into(),
        mass: 42.010565,
        name: None,
        search_mode: "mass_offset".into(),
        peptidoforms: 1,
        modified_peptidoforms: 1,
        eligible_sites: 1,
        modified_sites: 1,
        occupancy_pct: 100.0,
    };
    assert_eq!(modification_label(&bare), "K +42.0106");
    assert_eq!(
        modification_label(&ModificationOccupancy {
            name: Some("Acetyl".into()),
            ..bare
        }),
        "Acetyl"
    );
}
