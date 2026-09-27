//! Beta 12 database features combined with older ones in one search: initiator
//! Met clipping, ambiguous-residue expansion and J, a PTM library with a
//! protein N-terminal record on a clipped protein and a residue record in an
//! X-expanded peptide, a `<`-anchored motif, `max_total_count`, semi-enzymatic
//! digestion, generated decoys and the streamed prefilter.

use sage_core::enzyme::Position;
use sage_core::ion_series::{IonSeries, Kind};
use std::collections::BTreeSet;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

const FASTA: &str = ">sp|CLIP|ACETYLATED\nMASPEPTIDEAAKGGLLR\n\
                     >sp|MOTIF|THREONINE\nMTNLEPGYDSLKGGR\n\
                     >sp|AMBIG|MIXED\nMKGWEMLXDFMGRYBDSFPAGHKVLGJAPQTWEK\n";

/// One-based M6 of AMBIG sits in the X-containing digest GWEMLXDFMGR.
const LIBRARY: &str = "protein\tposition\tresidue\tmodification\tattachment\n\
                       sp|CLIP|ACETYLATED\t2\tA\tAcetyl\tprotein_n_term\n\
                       sp|AMBIG|MIXED\t6\tM\tOxidation\tresidue\n";

/// Spectra are written for these database peptides.
const TRUTHS: [&str; 5] = [
    // Clipped, acetylated only through the library's protein N-terminal record.
    "[Acetyl]-ASPEPTIDEAAK",
    // Clipped; the motif's `<` anchor matches residue 2.
    "T[AcetylMotif]NLEPGYDSLK",
    // X as A; one library oxidation plus one new one needs max_total_count 2.
    "GWEM[Oxidation]LADFM[Oxidation]GR",
    // B as N.
    "YNDSFPAGHK",
    // J searched with the I/L mass and reported as J.
    "VLGJAPQTWEK",
];

fn database(fasta: &std::path::Path, library: &std::path::Path) -> serde_json::Value {
    serde_json::json!({
        "fasta": fasta,
        "enzyme": {"min_len": 5, "semi_enzymatic": true},
        "static_mods": {},
        "clip_n_term_met": true,
        "expand_ambiguous_residues": true,
        "generate_decoys": true,
        "prefilter_min_matched_peaks": 1,
        "variable_mods": {
            "Oxidation": {
                "mass": 15.994915, "sites": ["M"], "site_mode": "both",
                "max_count": 1, "max_total_count": 2
            },
            "Acetyl": {
                "mass": 42.010565, "sites": ["protein_n_term"], "site_mode": "library",
                "max_count": 1
            },
            "AcetylMotif": {"mass": 42.010565, "sites": ["motif:<[AGST]*"], "max_count": 1}
        },
        "ptm_library": {"path": library, "strict": true}
    })
}

#[test]
fn combined_database_features_match_with_and_without_prefilter() -> anyhow::Result<()> {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let root = std::env::temp_dir().join(format!("sage-combined-{}-{nonce}", std::process::id()));
    std::fs::create_dir_all(&root)?;
    let fasta_path = root.join("proteins.fasta");
    let library_path = root.join("sites.tsv");
    std::fs::write(&fasta_path, FASTA)?;
    std::fs::write(&library_path, LIBRARY)?;
    let database = database(&fasta_path, &library_path);

    // Build the database in process to write theoretical spectra and to check
    // that every truth has an N-terminal-consistent generated decoy.
    let mut parameters =
        serde_json::from_value::<sage_core::database::Builder>(database.clone())?.make_parameters();
    let library =
        sage_core::ptm_library::PtmLibrary::from_tsv(LIBRARY).map_err(anyhow::Error::msg)?;
    parameters
        .validate_ptm_library(&library)
        .map_err(anyhow::Error::msg)?;
    parameters.loaded_ptm_library = Some(std::sync::Arc::new(library));
    let parsed = sage_core::fasta::Fasta::parse(FASTA.into(), "rev_", true)?;
    let peptides = parameters.modify_digests(parameters.digest_unmodified(&parsed));

    let proton = sage_core::mass::PROTON;
    let mut mgf = String::new();
    for (index, truth) in TRUTHS.iter().enumerate() {
        let peptide = peptides
            .iter()
            .find(|peptide| !peptide.decoy && peptide.to_string() == *truth)
            .unwrap_or_else(|| panic!("{truth} is not in the database"));
        assert!(!peptide.semi_enzymatic, "{truth} is fully enzymatic");
        if index < 2 {
            assert_eq!(peptide.position, Position::Nterm, "{truth} is clipped");
            let decoy = peptides
                .iter()
                .find(|decoy| {
                    decoy.decoy
                        && decoy.monoisotopic == peptide.monoisotopic
                        && decoy.position == Position::Nterm
                        && decoy.sequence.first() == peptide.sequence.first()
                        && decoy.sequence.len() == peptide.sequence.len()
                })
                .unwrap_or_else(|| panic!("{truth} has no N-terminal decoy"));
            assert_ne!(decoy.sequence, peptide.sequence);
        }
        let mut ions = [Kind::B, Kind::Y]
            .into_iter()
            .flat_map(|kind| IonSeries::new(peptide, kind).map(|ion| ion.monoisotopic_mass))
            .collect::<Vec<_>>();
        ions.sort_by(f32::total_cmp);
        mgf.push_str(&format!(
            "BEGIN IONS\nTITLE=truth-{index}\nSCANS={}\nRTINSECONDS={}\nPEPMASS={:.6}\nCHARGE=2+\n",
            index + 1,
            60 * (index + 1),
            (peptide.monoisotopic + 2.0 * proton) / 2.0
        ));
        for (rank, mass) in ions.iter().enumerate() {
            mgf.push_str(&format!("{:.6} {}\n", mass + proton, 1000 - rank));
        }
        mgf.push_str("END IONS\n");
    }
    std::fs::write(root.join("spectra.mgf"), mgf)?;

    let mut outputs = Vec::new();
    for prefilter in [false, true] {
        let mut database = database.clone();
        database["prefilter"] = prefilter.into();
        let config = serde_json::json!({
            "database": database,
            "mzml_paths": [root.join("spectra.mgf")],
            "deisotope": false,
            "precursor_tol": {"ppm": [-10, 10]},
            "fragment_tol": {"ppm": [-10, 10]},
            "min_matched_peaks": 4,
            "report_psms": 3,
            "write_pin": true,
            "diagnostic_ions": true,
            "output_filter": {"psm_q_value": 1.0}
        });
        let run = root.join(format!("prefilter-{prefilter}"));
        std::fs::create_dir_all(&run)?;
        std::fs::write(run.join("config.json"), serde_json::to_vec_pretty(&config)?)?;
        let output = Command::new(env!("CARGO_BIN_EXE_sage"))
            .arg(run.join("config.json"))
            .arg("--output_directory")
            .arg(run.join("output"))
            .arg("--disable-telemetry-i-dont-want-to-improve-sage")
            .output()?;
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        // Every pin column of every PSM, independent of row order.
        let pin = std::fs::read_to_string(run.join("output/results.sage.pin"))?;
        let mut lines = pin.lines();
        let header = lines.next().expect("pin header").to_string();
        let rows = lines.map(str::to_string).collect::<BTreeSet<_>>();
        let digestion = std::fs::read_to_string(run.join("output/digestion.tsv"))?;
        outputs.push((header, rows, digestion));
    }

    let (header, rows, _) = &outputs[0];
    let columns = header.split('\t').collect::<Vec<_>>();
    let column = |name: &str| columns.iter().position(|c| *c == name).expect(name);
    let (label, rank, peptide, semi) = (
        column("Label"),
        column("rank"),
        column("Peptide"),
        column("semi_enzymatic"),
    );
    for truth in TRUTHS {
        let row = rows
            .iter()
            .map(|row| row.split('\t').collect::<Vec<_>>())
            .find(|row| row[rank] == "1" && row[peptide] == truth)
            .unwrap_or_else(|| panic!("{truth} is not a rank-1 PSM:\n{rows:#?}"));
        assert_eq!(row[label], "1", "{truth}");
        assert_eq!(row[semi], "0", "{truth} is not semi-enzymatic");
    }
    // Every pin row and column, and the digestion summary, are unchanged.
    assert_eq!(outputs[0], outputs[1], "prefilter changed the results");

    std::fs::remove_dir_all(root)?;
    Ok(())
}
