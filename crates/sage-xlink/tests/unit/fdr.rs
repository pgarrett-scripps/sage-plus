use super::*;

#[test]
fn fdr_subtracts_dd_from_td() {
    use Class::*;
    let items = [
        (10.0, TT),
        (9.0, TT),
        (8.0, TD),
        (7.0, TT),
        (6.0, TT),
        (5.0, DD),
        (4.0, TD),
        (3.0, TT),
    ];
    let q = q_values(&items);
    assert_eq!(q[0], 0.0);
    assert_eq!(q[1], 0.0);
    // At score 6: 4 TT, 1 TD, 0 DD -> 0.25; at score 5 still 1 TD - 1 DD.
    assert_eq!(q[4], 0.0);
    assert_eq!(q[2], 0.0);
    // At score 3: 5 TT, 2 TD, 1 DD -> 0.2.
    assert!((q[7] - 0.2).abs() < 1e-6);
}

#[test]
fn q_values_are_monotone_and_tie_aware() {
    use Class::*;
    let items = [(5.0, TT), (5.0, TD), (4.0, TT), (3.0, TT)];
    let q = q_values(&items);
    assert_eq!(q[0], q[1]);
    assert!((q[3] - 1.0 / 3.0).abs() < 1e-6);
    assert!((q[0] - 1.0 / 3.0).abs() < 1e-6);
}

fn csm(score: f64, q: f32, class: Class) -> Csm {
    let chain = crate::search::ChainMatch {
        peptide: sage_core::database::PeptideIx(0),
        site: sage_core::peptide::Site::Sequence(0),
        hyperscore: 0.0,
        matched_peaks: 0,
        doublet: false,
        length: 8,
    };
    Csm {
        file_id: 0,
        spectrum_id: String::new(),
        rt: 0.0,
        charge: 3,
        expmass: 0.0,
        calcmass: 0.0,
        isotope_error: 0,
        precursor_ppm: 0.0,
        alpha: chain.clone(),
        beta: chain,
        hyperscore: 0.0,
        delta_next: 0.0,
        matched_peaks: 0,
        matched_intensity_pct: 0.0,
        doublets: 0,
        class,
        intra: true,
        discriminant_score: score,
        csm_q: q,
        residue_pair_q: 1.0,
    }
}

#[test]
fn pair_evidence_counts_only_supporting_csms() {
    let mut pair = PairEvidence::new(&csm(2.0, 0.2, Class::TT));
    // A single CSM that does not pass support is still scored by its best score.
    assert_eq!(pair.features(), [2.0, 0.0]);
    pair.add(&csm(5.0, 0.01, Class::TT));
    pair.add(&csm(1.0, 0.04, Class::TT));
    assert_eq!(pair.best, 5.0);
    assert_eq!(pair.supporting, 2);
    assert!((pair.features()[1] - 3f64.ln()).abs() < 1e-12);
}

#[test]
fn weak_chain_coverage_uses_the_less_covered_chain() {
    let mut c = csm(0.0, 1.0, Class::TT);
    c.alpha.matched_peaks = 8;
    c.alpha.length = 8;
    c.beta.matched_peaks = 3;
    c.beta.length = 12;
    assert_eq!(features(&c)[FEATURES - 1], 0.25);
    c.beta.length = 0;
    assert_eq!(per_residue(&c.beta), 3.0);
}

#[test]
fn best_csm_sets_pair_class() {
    let mut pair = PairEvidence::new(&csm(1.0, 0.5, Class::TD));
    pair.add(&csm(3.0, 0.5, Class::TT));
    assert_eq!(pair.class, Class::TT);
}

#[test]
fn groups_are_independent() {
    use Class::*;
    let items = [(10.0, TT, true), (9.0, TD, false), (8.0, TT, false)];
    let q = grouped_q_values(&items);
    assert_eq!(q[0], 0.0);
    assert_eq!(q[2], 1.0);
}
