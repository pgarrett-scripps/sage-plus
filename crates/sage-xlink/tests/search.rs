//! End-to-end crosslink search on a synthetic DSSO spectrum.

use sage_core::database::{Builder, IndexedDatabase};
use sage_core::fasta::Fasta;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::{Tolerance, PROTON};
use sage_core::scoring::{ScoreType, Scorer};
use sage_core::spectrum::{Precursor, ProcessedSpectrum};
use sage_xlink::{Class, CleavableLinker, CrosslinkSearch, CrosslinkSettings};

const FASTA: &str = "
>sp|P00001|ONE
MSTRAGLEKVDFPSRWQTEIKHNNGRLLSPEDAVKR
>sp|P00002|TWO
MKRYTDEAKLNWPGRHHSSDFKGAMEVLPRTTNQEFGKR
";
const ALPHA: &str = "AGLEKVDFPSR";
const BETA: &str = "YTDEAKLNWPGR";

fn database() -> IndexedDatabase {
    let builder: Builder = serde_json::from_value(serde_json::json!({
        "bucket_size": 8192,
        "fasta": "static",
        "enzyme": {"missed_cleavages": 1, "min_len": 5, "max_len": 30},
        "static_mods": {"C": 57.021464},
        "variable_mods": {},
    }))
    .unwrap();
    let fasta = Fasta::parse(FASTA.into(), "rev_", true).unwrap();
    builder.make_parameters().build(fasta)
}

fn scorer(db: &IndexedDatabase) -> Scorer<'_> {
    Scorer {
        db,
        precursor_tol: Tolerance::Ppm(-10.0, 10.0),
        fragment_tol: Tolerance::Ppm(-10.0, 10.0),
        min_matched_peaks: 4,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 2,
        max_precursor_charge: 4,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms: 1,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: 5.0,
        score_type: ScoreType::SageHyperScore,
        mass_recalibration: None,
    }
}

/// Neutral fragment masses of a released chain: b and y ions, with the ions
/// that contain the linked lysine carrying the given stub.
fn chain_fragments(db: &IndexedDatabase, sequence: &str, stub: f32) -> (f32, Vec<f32>) {
    let peptide = db
        .peptides
        .iter()
        .find(|p| !p.decoy && std::str::from_utf8(&p.sequence).unwrap() == sequence)
        .unwrap_or_else(|| panic!("{sequence} not in database"));
    let site = sequence.find('K').unwrap();
    let mut masses = Vec::new();
    for kind in [Kind::B, Kind::Y] {
        for (index, ion) in IonSeries::new(peptide, kind).enumerate() {
            // b_{i+1} covers residues 0..=i; y ions are emitted from the
            // longest down, y_{n-1-i} covers residues i+1..n.
            let carries = match kind {
                Kind::B => index >= site,
                _ => index < site,
            };
            masses.push(ion.monoisotopic_mass + if carries { stub } else { 0.0 });
        }
    }
    (peptide.monoisotopic, masses)
}

fn spectrum(db: &IndexedDatabase, linker: &CleavableLinker) -> ProcessedSpectrum {
    let mut peaks: Vec<(f32, u8, f32)> = Vec::new();
    let mut total = 0.0;
    for (sequence, stub) in [(ALPHA, linker.light_stub), (BETA, linker.heavy_stub)] {
        let (chain, fragments) = chain_fragments(db, sequence, stub);
        total += chain;
        // Signature doublet at charge 2.
        peaks.push(((chain + linker.light_stub) / 2.0, 2, 200.0));
        peaks.push(((chain + linker.heavy_stub) / 2.0, 2, 150.0));
        peaks.extend(fragments.into_iter().map(|mass| (mass, 1, 50.0)));
    }
    // Noise.
    peaks.extend([(333.3, 1, 10.0), (777.7, 1, 10.0), (1111.1, 1, 10.0)]);
    let precursor_mass = total + linker.crosslink_mass;
    // Stored as neutral mass with its charge; doublets are charge 2, so the
    // stored mass is the neutral chain+stub mass.
    let mut stored: Vec<(f32, u8, f32)> = peaks
        .into_iter()
        .map(|(mass, charge, intensity)| (mass * charge as f32, charge, intensity))
        .collect();
    stored.sort_by(|a, b| a.0.total_cmp(&b.0));
    ProcessedSpectrum {
        level: 2,
        id: "controllerType=0 controllerNumber=1 scan=42".into(),
        precursors: vec![Precursor {
            mz: precursor_mass / 4.0 + PROTON,
            charge: Some(4),
            ..Default::default()
        }],
        masses: stored.iter().map(|p| p.0).collect(),
        charges: stored.iter().map(|p| p.1).collect(),
        charge_is_known: vec![true; stored.len()],
        intensities: stored.iter().map(|p| p.2).collect(),
        total_ion_current: stored.iter().map(|p| p.2).sum(),
        ..Default::default()
    }
}

#[test]
fn synthetic_dsso_spectrum_identifies_both_chains() {
    let db = database();
    let settings = CrosslinkSettings {
        isotope_errors: (0, 0),
        ..Default::default()
    };
    let search = CrosslinkSearch::new(settings.clone()).unwrap();
    let query = spectrum(&db, &search.linker);
    let csm = search
        .search(&scorer(&db), &query)
        .expect("a crosslink match");

    let mut chains = [
        std::str::from_utf8(&db[csm.alpha.peptide].sequence).unwrap(),
        std::str::from_utf8(&db[csm.beta.peptide].sequence).unwrap(),
    ];
    chains.sort();
    assert_eq!(chains, [ALPHA, BETA]);
    assert_eq!(csm.class, Class::TT);
    assert!(!csm.intra);
    assert_eq!(csm.charge, 4);
    assert_eq!(csm.isotope_error, 0);
    assert!(csm.precursor_ppm.abs() < 1.0, "{}", csm.precursor_ppm);
    assert!(csm.alpha.doublet && csm.beta.doublet);
    // Both chains' fragments count towards the combined match.
    assert!(csm.matched_peaks as usize >= 2 * (ALPHA.len() + BETA.len()) - 8);
    assert!(csm.hyperscore > csm.alpha.hyperscore);

    // Linked lysines: ALPHA K5 at P00001 position 9, BETA K6 at P00002 position 9.
    for chain in [&csm.alpha, &csm.beta] {
        let (_, position) = sage_xlink::search::protein_position(&db[chain.peptide], chain.site)
            .expect("protein position");
        assert_eq!(position, 9);
    }

    let mut csms = vec![csm];
    sage_xlink::assign_q_values(&mut csms, &db, 0.01);
    let bytes = sage_xlink::output::serialize(
        &csms.iter().collect::<Vec<_>>(),
        &db,
        &["synthetic.mzML".into()],
        &search.linker.name,
    )
    .unwrap();
    assert!(bytes.starts_with(b"PAR1"));
}

#[test]
fn spectrum_without_doublets_has_no_crosslink() {
    let db = database();
    let search = CrosslinkSearch::new(CrosslinkSettings::default()).unwrap();
    let mut query = spectrum(&db, &search.linker);
    let keep: Vec<bool> = query.charges.iter().map(|&c| c == 1).collect();
    let filter = |values: &mut Vec<_>| {
        let mut index = 0;
        values.retain(|_| {
            index += 1;
            keep[index - 1]
        });
    };
    filter(&mut query.masses);
    filter(&mut query.intensities);
    let mut index = 0;
    query.charges.retain(|_| {
        index += 1;
        keep[index - 1]
    });
    query.charge_is_known.truncate(query.masses.len());
    assert!(search.search(&scorer(&db), &query).is_none());
}
