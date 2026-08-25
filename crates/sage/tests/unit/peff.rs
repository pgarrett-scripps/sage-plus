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
