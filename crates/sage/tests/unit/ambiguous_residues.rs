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
