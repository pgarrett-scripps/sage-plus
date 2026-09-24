use super::*;

fn ctx<'a>(left: &'a [u8], right: &'a [u8]) -> MotifContext<'a> {
    MotifContext {
        left,
        right,
        ..Default::default()
    }
}

#[test]
fn parses_and_canonicalizes() {
    let motif = SiteMotif::parse("N*-{P}-[TS]").unwrap();
    assert_eq!(motif.canonical(), "motif:N*-{P}-[ST]");
    assert_eq!(motif.reach(), (0, 2));
    assert_eq!(motif.site_residues(), b"N");
    let kinase = SiteMotif::parse("R-x(2)-[ST]*").unwrap();
    assert_eq!(kinase.canonical(), "motif:R-x(2)-[ST]*");
    assert_eq!(kinase.reach(), (3, 0));
    assert!(SiteMotif::parse("N-{P}-[ST]").is_err());
    assert!(SiteMotif::parse("N*-[ST]*").is_err());
    assert!(SiteMotif::parse("x(2)*").is_err());
    assert!(SiteMotif::parse("N*-[ZB]").is_err());
    assert!(SiteMotif::parse("N*-x(64)").is_err());
    assert_eq!(
        SiteMotif::intern("N*-{P}-[ST]").unwrap() as *const _,
        SiteMotif::intern("N*-{P}-[TS]").unwrap() as *const _
    );
}

#[test]
fn sequon_matches_peptide_and_flanks() {
    let motif = SiteMotif::parse("N*-{P}-[ST]").unwrap();
    // NGT sequon, NPS is excluded by the proline rule.
    assert_eq!(motif.sites(b"ANGTANPSK", ctx(b"", b"")), vec![1]);
    // N second-to-last needs one flanking residue after the peptide.
    assert_eq!(motif.sites(b"AANK", ctx(b"", b"")), Vec::<u32>::new());
    assert_eq!(motif.sites(b"AANK", ctx(b"", b"S")), vec![2]);
    assert_eq!(motif.sites(b"AANK", ctx(b"", b"P")), Vec::<u32>::new());
}

#[test]
fn kinase_motif_uses_left_flank() {
    let motif = SiteMotif::parse("R-x-x-[ST]*").unwrap();
    assert_eq!(motif.sites(b"GASPK", ctx(b"R", b"")), vec![2]);
    assert_eq!(motif.sites(b"GASPK", ctx(b"", b"")), Vec::<u32>::new());
}

#[test]
fn protein_terminal_anchor() {
    let motif = SiteMotif::parse("C*-x-x-x>").unwrap();
    let terminal = MotifContext {
        right_boundary: true,
        ..Default::default()
    };
    assert_eq!(motif.sites(b"GKCVLS", terminal), vec![2]);
    assert_eq!(motif.sites(b"GKCVLS", ctx(b"", b"")), Vec::<u32>::new());
    assert_eq!(motif.sites(b"GKCVLSA", terminal), Vec::<u32>::new());
}

#[test]
fn canonical_classes_reparse_identically() {
    for pattern in [
        "N*-{P}-[ST]",
        "<M-x-[ST]*",
        "x-{ACDEFGHIKLMPQRSTVWYUO}*",
        "C*-x(3)>",
    ] {
        let motif = SiteMotif::parse(pattern).unwrap();
        let spelled = motif.canonical().strip_prefix("motif:").unwrap();
        let reparsed = SiteMotif::parse(spelled).unwrap();
        assert_eq!(reparsed.canonical(), motif.canonical());
        assert_eq!(reparsed.elements, motif.elements);
    }
    // `x` and exclusions match unknown FASTA residues; positive classes do not.
    let motif = SiteMotif::parse("N*-x-[ST]").unwrap();
    assert_eq!(motif.sites(b"NXS", ctx(b"", b"")), vec![0]);
    assert!(SiteMotif::parse("{ACDEFGHIKLMNPQRSTVWYUO}*").is_err());
}

mod database {
    use crate::database::{Builder, Parameters};
    use crate::enzyme::Position;
    use crate::fasta::Fasta;
    use crate::modification::ModificationSpecificity;
    use crate::peptide::{Peptide, Site};
    use std::collections::HashSet;

    const HEXNAC: &str = "motif:N*-{P}-[ST]";

    fn parameters(extra: serde_json::Value) -> Parameters {
        let mut config = serde_json::json!({
            "enzyme": {"min_len": 4},
            "peptide_min_mass": 100.0,
            "static_mods": {},
            "variable_mods": {
                "HexNAc": {"mass": 203.079373, "sites": [HEXNAC], "max_count": 1},
                "Phospho": {"mass": 79.966331, "sites": ["motif:R-x(2)-[ST]*"], "max_count": 1}
            }
        });
        for (key, value) in extra.as_object().unwrap() {
            config[key] = value.clone();
        }
        serde_json::from_value::<Builder>(config)
            .unwrap()
            .make_parameters()
    }

    fn expand(parameters: &Parameters, fasta: &str, decoys: bool) -> Vec<Peptide> {
        let fasta = Fasta::parse(fasta.into(), "rev_", decoys).unwrap();
        parameters.modify_digests(parameters.digest_unmodified(&fasta))
    }

    fn rendered(peptides: &[Peptide]) -> HashSet<(bool, String)> {
        peptides
            .iter()
            .map(|peptide| (peptide.decoy, peptide.to_string()))
            .collect()
    }

    fn motif() -> ModificationSpecificity {
        HEXNAC.parse().unwrap()
    }

    #[test]
    fn motifs_use_the_full_source_protein() {
        let parameters = parameters(serde_json::json!({}));
        // MR|AK|SLLGGK|ANPSK|AANK|STTR: R-x-x-S needs three residues before
        // SLLGGK; N-K-S in AANK needs the residue after it.
        let peptides = expand(&parameters, ">P1\nMRAKSLLGGKANPSKAANKSTTR\n", true);
        let shown = rendered(&peptides);
        assert!(
            shown.contains(&(false, "S[Phospho]LLGGK".into())),
            "{shown:?}"
        );
        assert!(shown.contains(&(false, "AAN[HexNAc]K".into())), "{shown:?}");
        assert!(!shown
            .iter()
            .any(|(decoy, p)| !decoy && p.starts_with("AN[")));

        // An Asn-C style digest leaves N last; the sequon spans two residues
        // beyond the peptide.
        let parameters = parameters_with_enzyme("N");
        let peptides = expand(&parameters, ">P1\nGGGGNKSGGGGNPS\n", true);
        let shown = rendered(&peptides);
        assert!(
            shown.contains(&(false, "GGGGN[HexNAc]".into())),
            "{shown:?}"
        );
        assert!(
            !shown.contains(&(false, "KSGGGGN[HexNAc]".into())),
            "{shown:?}"
        );
    }

    fn parameters_with_enzyme(cleave: &str) -> Parameters {
        parameters(serde_json::json!({
            "enzyme": {"min_len": 4, "cleave_at": cleave, "restrict": null, "missed_cleavages": 0}
        }))
    }

    #[test]
    fn protein_terminal_anchor_follows_protein_coordinates() {
        let parameters = parameters(serde_json::json!({
            "variable_mods": {
                "Farnesyl": {"mass": 204.187801, "sites": ["motif:C*-x(3)>"], "max_count": 1}
            }
        }));
        let peptides = expand(&parameters, ">P1\nMAGCKGGCVLSK\n>P2\nGGCVLSKAAAAK\n", false);
        let shown = rendered(&peptides);
        assert!(
            !shown.iter().any(|(_, p)| p.contains("C[Farnesyl]VLSK")),
            "{shown:?}"
        );
        let peptides = expand(&parameters, ">P1\nMAGCKGGCVLS\n", false);
        let shown = rendered(&peptides);
        assert!(
            shown.contains(&(false, "GGC[Farnesyl]VLS".into())),
            "{shown:?}"
        );
        assert!(!shown.iter().any(|(_, p)| p.contains("AGC[")), "{shown:?}");
    }

    #[test]
    fn shared_peptides_are_eligible_if_any_occurrence_matches() {
        let parameters = parameters(serde_json::json!({}));
        // GAANK is followed by S in P1 but by P in P2.
        let fasta = ">P1\nMRGAANKSLLR\n>P2\nMRGAANKPLLR\n";
        let peptides = expand(&parameters, fasta, true);
        let modified = peptides
            .iter()
            .find(|peptide| !peptide.decoy && peptide.to_string() == "GAAN[HexNAc]K")
            .expect("shared peptide carries the motif site");
        assert_eq!(modified.protein_sites.len(), 2);
        assert_eq!(modified.rule_sites(motif()), vec![Site::Sequence(3)]);
        let rules = parameters.localization_rules(modified);
        assert!(rules
            .iter()
            .any(|rule| rule.specificity == motif() && rule.sites == vec![Site::Sequence(3)]));
        // Retention and mobility features count the placed motif modification.
        assert_eq!(
            crate::ml::retention_model::variable_mod_count(modified, motif(), 203.079373),
            1.0
        );
        // Each occurrence alone decides differently.
        let only = |index: usize| {
            motif().sites_for_occurrences(
                &modified.sequence,
                modified.position,
                false,
                &modified.protein_sites[index..=index],
            )
        };
        assert_ne!(only(0).is_empty(), only(1).is_empty());

        // The memory estimate sees the same union.
        let groups =
            parameters.digest_unmodified(&Fasta::parse(fasta.into(), "rev_", true).unwrap());
        let estimate = parameters.estimate_modified_memory(&groups);
        let targets = peptides.iter().filter(|peptide| !peptide.decoy).count() as u64;
        assert!(estimate.modified_peptides >= 2 * targets, "{estimate:?}");
    }

    #[test]
    fn generated_decoys_mirror_target_motif_sites() {
        let parameters = parameters(serde_json::json!({}));
        let peptides = expand(&parameters, ">P1\nMRGAANKSLLRGANGTAAK\n", true);
        let shown = rendered(&peptides);
        // GAANK -> GNAAK: the mirrored N is not literally a sequon.
        let decoy = peptides
            .iter()
            .find(|peptide| peptide.decoy && peptide.to_string() == "GN[HexNAc]AAK")
            .unwrap_or_else(|| panic!("{shown:?}"));
        assert_eq!(decoy.rule_sites(motif()), vec![Site::Sequence(1)]);
        let literal = motif().sites_for_occurrences(&decoy.sequence, decoy.position, false, &[]);
        assert!(literal.is_empty());
        let rules = parameters.localization_rules(decoy);
        assert!(rules
            .iter()
            .any(|rule| rule.specificity == motif() && rule.sites == vec![Site::Sequence(1)]));
        let count = |decoy: bool| {
            peptides
                .iter()
                .filter(|peptide| peptide.decoy == decoy)
                .map(|peptide| peptide.modifications.len())
                .sum::<usize>()
        };
        assert_eq!(count(false), count(true));
    }

    #[test]
    fn fasta_decoys_are_matched_literally() {
        let parameters = parameters(serde_json::json!({"generate_decoys": false}));
        // The decoy protein is a literal reversal: SKNAAGRM contains no sequon
        // at GRNAAK... but rev_P1 below contains its own literal N-G-T.
        let fasta = ">P1\nMRGAANKSLLR\n>rev_P1\nMRNGTAAKR\n";
        let peptides = expand(&parameters, fasta, false);
        let shown = rendered(&peptides);
        assert!(
            shown.contains(&(false, "GAAN[HexNAc]K".into())),
            "{shown:?}"
        );
        let decoy = peptides
            .iter()
            .find(|peptide| peptide.decoy && peptide.to_string() == "N[HexNAc]GTAAK")
            .unwrap_or_else(|| panic!("{shown:?}"));
        assert_eq!(decoy.rule_sites(motif()), vec![Site::Sequence(0)]);
        assert!(!shown.iter().any(|(decoy, p)| *decoy && p.contains("GNAAK")));
    }

    #[test]
    fn mass_offset_motifs_use_the_same_sites() {
        let parameters = parameters(serde_json::json!({
            "variable_mods": {
                "HexNAc": {"mass": 203.079373, "sites": [HEXNAC], "search_mode": "mass_offset"}
            }
        }));
        let fasta = Fasta::parse(">P1\nMRGAANKSLLRGANPTK\n".into(), "rev_", true).unwrap();
        let database = parameters.clone().build(fasta);
        let offsets = parameters.mass_offset_modifications();
        assert_eq!(offsets.len(), 1);
        let sites = |sequence: &str, decoy: bool| {
            let peptide = database
                .peptides
                .iter()
                .find(|peptide| peptide.decoy == decoy && peptide.sequence == sequence)
                .unwrap();
            database.mass_offset_sites(peptide, &offsets[0])
        };
        assert_eq!(sites("GAANK", false), vec![Site::Sequence(3)]);
        assert_eq!(sites("GNAAK", true), vec![Site::Sequence(1)]);
        assert!(sites("GANPTK", false).is_empty());
        assert!(sites("GTPNAK", true).is_empty());
    }

    #[test]
    fn peptide_lists_without_coordinates_use_the_peptide_only() {
        let parameters = parameters(serde_json::json!({}));
        let peptides = parameters.peptides_from_tsv("sequence\nGANGTK\nGAANK\n");
        let shown = rendered(&peptides);
        assert!(
            shown.contains(&(false, "GAN[HexNAc]GTK".into())),
            "{shown:?}"
        );
        assert!(
            shown.contains(&(true, "GTGN[HexNAc]AK".into())),
            "{shown:?}"
        );
        assert!(!shown.iter().any(|(_, p)| p.contains("AAN[")));
        let _ = Position::Internal;
    }

    #[test]
    fn motif_sites_require_named_definitions() {
        let legacy = serde_json::from_value::<Builder>(serde_json::json!({
            "variable_mods": {"motif:N*-{P}-[ST]": [203.079373]}
        }));
        let error = legacy.err().unwrap().to_string();
        assert!(error.contains("requires a named definition"), "{error}");
        for sites in [
            ["motif:N*-{P}-[TS]"],
            ["motif:N-{P}-[ST]"],
            ["motif:N*-[ST]*"],
        ] {
            assert!(serde_json::from_value::<Builder>(serde_json::json!({
                "variable_mods": {"HexNAc": {"mass": 203.0, "sites": sites}}
            }))
            .is_err());
        }
        let error = serde_json::from_value::<Builder>(serde_json::json!({
            "variable_mods": {"HexNAc": {"mass": 203.0, "sites": ["motif:N-{P}-[ST]"]}}
        }))
        .err()
        .unwrap()
        .to_string();
        assert!(error.contains("must mark the modified residue"), "{error}");
    }
}
