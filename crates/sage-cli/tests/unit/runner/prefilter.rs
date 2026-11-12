use super::*;

#[test]
fn index_budget_follows_max_memory() {
    const GIB: u64 = 1024 * 1024 * 1024;
    // A quarter of 16 GiB indexes spectra; 8 GiB without a limit.
    assert_eq!(
        PrefilterBudgets::from_max_memory(Some(16.0)).index_bytes,
        4 * GIB
    );
    assert_eq!(PrefilterBudgets::from_max_memory(None).index_bytes, 8 * GIB);
    assert_eq!(
        PrefilterBudgets::from_max_memory(Some(0.0)).index_bytes,
        8 * GIB
    );
}
#[test]
fn isoleucine_leucine_twins_share_a_prefilter_key() {
    // APEPIDEK and APEPJDEK are twins; DFNEEGIK is a twin of the decoy
    // of DLGEENFK. Each digest is unique as written.
    let fasta = Fasta::parse(
        ">a\nGGRAPEPIDEKR\n>b\nWWKAPEPJDEKR\n>c\nDLGEENFKR\n>d\nMSSRDFNEEGIKR\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let enzyme = sage_core::database::Builder::default()
        .make_parameters()
        .enzyme_parameters();
    for merge in [false, true] {
        let shared = SharedSequences::scan(&fasta, &enzyme, None, true, merge);
        let key = |protein: usize, sequence: &str| {
            let digest = fasta
                .digest_protein(protein, &enzyme, None)
                .into_iter()
                .find(|digest| digest.sequence.as_str() == sequence)
                .unwrap();
            shared.shared_key(&digest)
        };
        let (i, j) = (key(0, "APEPIDEK"), key(1, "APEPJDEK"));
        let (target, collision) = (key(2, "DLGEENFK"), key(3, "DFNEEGIK"));
        match merge {
            true => {
                assert!(i.is_some() && i == j);
                assert!(target.is_some() && target == collision);
            }
            false => assert_eq!([i, j, target, collision], [None; 4]),
        }
    }
}
