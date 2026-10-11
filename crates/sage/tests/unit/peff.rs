use super::*;

const COMPLEX_PEFF: &str = "# PEFF 1.0\n\
        # GeneralComment=test\n\
        # //\n\
        # DbName=cx\n\
        # //\n\
        >cx:ENTRY1 \\PName=Test \\Length=20 \\ModResUnimod=(2,5|UNIMOD:21|Phospho) \\ModResPsi=(7|MOD:00048|whatever)\n\
        ACDEFGHIKLMNPQRSTVWY\n";

#[test]
fn detects_peff_header() {
    assert!(looks_like_peff(COMPLEX_PEFF));
    assert!(!looks_like_peff(">sp:P12345\nMKTL\n"));
}

#[test]
fn parses_unimod_mods() {
    let (targets, mods) = parse(COMPLEX_PEFF, "rev_", true);
    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].0.as_ref(), "cx:ENTRY1");
    assert_eq!(targets[0].1, "ACDEFGHIKLMNPQRSTVWY");
    let entry_mods = mods.get(&targets[0].0).expect("mods missing");
    assert_eq!(entry_mods.len(), 2);
    // 1-based -> 0-based: 2 -> 1, 5 -> 4
    let positions: Vec<u32> = entry_mods.iter().map(|m| m.protein_pos).collect();
    assert!(positions.contains(&1));
    assert!(positions.contains(&4));
    let phospho = 79.966_33_f32;
    for m in entry_mods {
        assert!((m.mass - phospho).abs() < 1e-3, "got {}", m.mass);
    }
}

#[test]
fn skips_non_integer_positions_and_unknown_accessions() {
    let peff = "# PEFF 1.0\n\
            # //\n\
            >sp:P0 \\ModResUnimod=(N,3|UNIMOD:35|Oxidation)(7|UNIMOD:99999999|Bogus)\n\
            ACDEFGHIK\n";
    let (_, mods) = parse(peff, "rev_", true);
    let entry = mods.values().next().expect("entry mods missing");
    assert_eq!(entry.len(), 1, "{:?}", entry);
    assert_eq!(entry[0].protein_pos, 2); // 3 - 1
}

#[test]
fn entry_without_mods_yields_no_map_entry() {
    let peff = "# PEFF 1.0\n\
            # //\n\
            >sp:P0 \\PName=Plain\n\
            ACDEFGHIK\n";
    let (targets, mods) = parse(peff, "rev_", true);
    assert_eq!(targets.len(), 1);
    assert!(mods.is_empty());
}

#[test]
fn split_paren_items_handles_nested() {
    let v = split_paren_items("(380||N-linked (GlcNAc...))(5|UNIMOD:21|Phospho)");
    assert_eq!(v.len(), 2);
    assert_eq!(v[0], "380||N-linked (GlcNAc...)");
    assert_eq!(v[1], "5|UNIMOD:21|Phospho");
}

#[test]
fn looks_like_peff_skips_leading_blank_lines() {
    assert!(looks_like_peff("\n   \n# PEFF 1.0\n>x\nAA\n"));
    // First non-empty line decides; a later PEFF header does not count.
    assert!(!looks_like_peff("\n>x\n# PEFF 1.0\n"));
    assert!(!looks_like_peff(""));
    assert!(!looks_like_peff("   \n\n"));
}

#[test]
fn multi_entry_multi_line_sequences_keep_mods_per_protein() {
    let peff = "# PEFF 1.0\n\
            >db:A \\ModResUnimod=(1|UNIMOD:35|Oxidation)\n\
            MKT\n\
            LLK\n\
            \n\
            >db:B \\PName=NoMods\n\
            PEP\n\
            >db:C \\ModResUnimod=(4|UNIMOD:21|Phospho)\n\
            GGGS\n";
    let (targets, mods) = parse(peff, "rev_", true);
    let accs: Vec<&str> = targets.iter().map(|(a, _)| a.as_ref()).collect();
    assert_eq!(accs, vec!["db:A", "db:B", "db:C"]);
    assert_eq!(targets[0].1, "MKTLLK");
    assert_eq!(targets[1].1, "PEP");
    assert_eq!(targets[2].1, "GGGS");

    assert_eq!(mods.len(), 2);
    let a = &mods[&targets[0].0];
    assert_eq!(a.len(), 1);
    assert_eq!(a[0].protein_pos, 0);
    assert_eq!(a[0].name, "Oxidation");
    assert!((a[0].mass - 15.994_915).abs() < 1e-3, "got {}", a[0].mass);
    assert!(!mods.contains_key(&targets[1].0));
    let c = &mods[&targets[2].0];
    assert_eq!(c[0].protein_pos, 3);
    assert_eq!(c[0].name, "Phospho");
}

#[test]
fn decoy_entries_dropped_only_when_generating_decoys() {
    let peff = "# PEFF 1.0\n\
            >db:T \\ModResUnimod=(1|UNIMOD:35|Oxidation)\n\
            MK\n\
            >rev_db:T \\ModResUnimod=(2|UNIMOD:35|Oxidation)\n\
            KM\n";
    let (targets, mods) = parse(peff, "rev_", true);
    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].0.as_ref(), "db:T");
    assert_eq!(mods.len(), 1);

    // Without decoy generation, decoys present in the file are kept (with mods).
    let (targets, mods) = parse(peff, "rev_", false);
    assert_eq!(targets.len(), 2);
    assert_eq!(targets[1].0.as_ref(), "rev_db:T");
    assert_eq!(targets[1].1, "KM");
    assert_eq!(mods[&targets[1].0][0].protein_pos, 1);
}

#[test]
fn sequence_lines_before_first_header_are_discarded() {
    let peff = "# PEFF 1.0\nORPHAN\n>db:X\nACD\n";
    let (targets, _) = parse(peff, "rev_", true);
    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].1, "ACD");
}

#[test]
fn canonicalizes_names_and_accepts_lowercase_prefix_and_bare_accession() {
    let peff = "# PEFF 1.0\n\
            >db:X \\ModResUnimod=(1|unimod:21|phospho)(2|35|MyOwnLabel)(3|UNIMOD: 4 |carbamidomethyl)\n\
            SMC\n";
    let (_, mods) = parse(peff, "rev_", true);
    let m = mods.values().next().unwrap();
    assert_eq!(m.len(), 3, "{m:?}");
    // A lowercased name maps back to the canonical Unimod spelling.
    assert_eq!(m[0].name, "Phospho");
    assert!((m[0].mass - 79.966_33).abs() < 1e-3);
    // An unknown label is kept verbatim; the mass still comes from the accession.
    assert_eq!(m[1].name, "MyOwnLabel");
    assert!((m[1].mass - 15.994_915).abs() < 1e-3);
    assert_eq!(m[2].name, "Carbamidomethyl");
    assert!((m[2].mass - 57.021_465).abs() < 1e-3);
    assert_eq!(
        m.iter().map(|m| m.protein_pos).collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
}

#[test]
fn skips_zero_position_short_items_and_non_numeric_accessions() {
    let peff = "# PEFF 1.0\n\
            >db:X \\ModResUnimod=(0|UNIMOD:21|Phospho)(1|UNIMOD:21)(2|UNIMOD:abc|Phospho)(3, 4 |UNIMOD:21|Phospho)\n\
            SSSSS\n";
    let (_, mods) = parse(peff, "rev_", true);
    let m = mods.values().next().unwrap();
    let pos: Vec<u32> = m.iter().map(|m| m.protein_pos).collect();
    assert_eq!(pos, vec![2, 3]);
}

#[test]
fn strips_annotation_id_prefix_from_positions() {
    let v = parse_mod_res_unimod("x \\ModResUnimod=(ModId1:7|UNIMOD:21|Phospho)");
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].protein_pos, 6);
}

#[test]
fn mod_res_unimod_parsing_stops_at_next_key_or_garbage() {
    // A following annotation must not be consumed as more items.
    let v = parse_mod_res_unimod(
        "\\ModResUnimod=(1|UNIMOD:21|Phospho) \\ModResPsi=(2|UNIMOD:21|Phospho)",
    );
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].protein_pos, 0);
    // Bare text between items ends the list.
    let v = parse_mod_res_unimod("\\ModResUnimod=(1|UNIMOD:21|Phospho)junk(2|UNIMOD:21|Phospho)");
    assert_eq!(v.len(), 1);
    // No key at all.
    assert!(parse_mod_res_unimod("\\ModRes=(1||Phospho)").is_empty());
}

#[test]
fn split_paren_items_drops_unterminated_item() {
    assert_eq!(split_paren_items("(a)(b"), vec!["a".to_string()]);
    assert_eq!(split_paren_items("  (a)  (b)"), vec!["a", "b"]);
    assert!(split_paren_items("").is_empty());
}

#[test]
fn split_pipe_fields_ignores_pipes_inside_parens() {
    assert_eq!(split_pipe_fields("1|X (a|b)|c"), vec!["1", "X (a|b)", "c"]);
    assert_eq!(split_pipe_fields("a||"), vec!["a", "", ""]);
    // An unbalanced close paren does not underflow the depth counter.
    assert_eq!(split_pipe_fields(")|a"), vec![")", "a"]);
}
