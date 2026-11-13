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
    let qc = String::from_utf8(serialize_qc(&stats, &DigestionRunStats::default(), &[]).unwrap())
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
            String::from_utf8(serialize_qc(&[], &DigestionRunStats::default(), &[]).unwrap())
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
