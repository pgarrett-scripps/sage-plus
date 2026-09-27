use super::*;

#[test]
fn counts_variants_per_ambiguous_residue() {
    assert_eq!(variant_count(b"PEPTIDE"), 1);
    assert_eq!(variant_count(b"PEPJIDE"), 1);
    assert_eq!(variant_count(b"PEPBIDE"), 2);
    assert_eq!(variant_count(b"PEPZIDE"), 2);
    assert_eq!(variant_count(b"PEPXIDE"), 20);
    assert_eq!(variant_count(b"BXZ"), 80);
    assert_eq!(variant_count(&[b'X'; 64]), usize::MAX);
}

#[test]
fn expands_every_combination_in_order() {
    assert_eq!(expand(b"PEPTIDE"), vec![b"PEPTIDE".to_vec()]);
    assert_eq!(
        expand(b"ABZ"),
        vec![
            b"ADE".to_vec(),
            b"ADQ".to_vec(),
            b"ANE".to_vec(),
            b"ANQ".to_vec()
        ]
    );
    let x = expand(b"PEPXIDE");
    assert_eq!(x.len(), 20);
    assert!(x.contains(&b"PEPTIDE".to_vec()));
    assert!(x.iter().all(|variant| !is_ambiguous(variant)));
}

#[test]
fn recognizes_expanded_sequences() {
    assert!(is_expansion_of(b"PEPXIDE", b"PEPTIDE"));
    assert!(is_expansion_of(b"PEPTIDE", b"PEPTIDE"));
    assert!(is_expansion_of(b"BZ", b"NQ"));
    assert!(!is_expansion_of(b"BZ", b"EQ"));
    assert!(!is_expansion_of(b"PEPXIDE", b"PEPTIDES"));
    assert!(!is_expansion_of(b"PEPTIDE", b"PEPXIDE"));
}

#[test]
fn formats_substitutions_in_position_order() {
    assert_eq!(format_substitutions(b"PEPTIDE", b"PEPTIDE"), "");
    assert_eq!(format_substitutions(b"PEPXIDE", b"PEPTIDE"), "X4T");
    assert_eq!(
        format_substitutions(b"GBPXIDZK", b"GDPKIDQK"),
        "B2D;X4K;Z7Q"
    );
    // Two-digit positions, and an X kept as a residue it may stand for.
    assert_eq!(
        format_substitutions(b"AAAAAAAAAXB", b"AAAAAAAAAAN"),
        "X10A;B11N"
    );
    // J is scored as I/L, not substituted.
    assert_eq!(format_substitutions(b"PEPJXDE", b"PEPLTDE"), "X5T");
}

#[test]
fn orients_spans_like_the_peptide() {
    assert_eq!(
        oriented_span(b"PEPXIDEK", b"PEPTIDEK").as_deref(),
        Some(&b"PEPXIDEK"[..])
    );
    // A generated decoy reverses the interior.
    assert_eq!(
        oriented_span(b"PEPXIDEK", b"PEDITPEK").as_deref(),
        Some(&b"PEDIXPEK"[..])
    );
    // Merged I/L twins may show another twin's residues.
    assert_eq!(
        oriented_span(b"PEPXJDEK", b"PEPTLDEK").as_deref(),
        Some(&b"PEPXJDEK"[..])
    );
    assert_eq!(
        oriented_span(b"PEPXIDEK", b"PEPTLDEK").as_deref(),
        Some(&b"PEPXIDEK"[..])
    );
    assert_eq!(oriented_span(b"PEPBIDEK", b"PEPTIDEK"), None);
    assert_eq!(oriented_span(b"PEPXIDEK", b"PEPTIDE"), None);
}
