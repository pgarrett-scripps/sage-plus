use super::*;

#[test]
fn q_values_follow_target_decoy_competition() {
    let scores = [5.0, 4.0, 3.0, 2.0, 1.0];
    let decoy = [false, false, true, false, true];
    let q = q_values(&scores, &decoy);
    // (1 + decoys) / targets, then cumulative minimum from the bottom.
    assert_eq!(q, vec![0.5, 0.5, 0.6666667, 0.6666667, 1.0]);
}

#[test]
fn core_only_means_no_y_ion_beyond_the_core() {
    let mut counts = ClassCounts::default();
    counts.generated[CORE] = 2;
    // Oxonium classes do not count as Y ions.
    counts.generated[NEUAC_ABSENT] = 2;
    assert!(is_core_only(&counts));
    counts.generated[crate::evidence::EXTENDED] = 1;
    assert!(!is_core_only(&counts));
}

#[test]
fn tied_scores_share_a_q_value() {
    let q = q_values(&[2.0, 2.0, 1.0], &[false, true, false]);
    assert_eq!(q[0], q[1]);
}

fn explanation(matched: u16, decoy_matched: u16) -> Explanation {
    let mut target = ClassCounts::default();
    target.generated[2] = 10;
    target.matched[2] = matched;
    let mut decoy = target;
    decoy.matched[2] = decoy_matched;
    Explanation {
        composition: 0,
        isotope: 0,
        adducts: 0,
        error_ppm: 1.0,
        target,
        decoy,
        y_intensity: 0.0,
    }
}

#[test]
fn more_matched_ions_win_and_ties_are_a_fair_coin() {
    let model = GlycanModel::new((0, 1), 1, 20.0);
    let score = |e: &Explanation, c: &ClassCounts| model.score(e, c, 0.1, 0.1);
    let e = explanation(6, 1);
    assert!(score(&e, &e.target) > score(&e, &e.decoy));

    let candidate = |explanations| GlycoCandidate {
        feature: sage_core::scoring::Feature {
            spec_id: "scan=1".into(),
            ..Default::default()
        },
        peptide_mass: 1000.0,
        oxonium_ions: 3,
        oxonium_fraction: 0.1,
        anchored: true,
        y1: true,
        hex_ratio: 0.0,
        random_y: 0.1,
        random_oxonium: 0.1,
        explanations,
        site: Default::default(),
        twin_hyperscore: None,
        ladder: Default::default(),
        siblings: 0,
    };
    let a = model.assign(&candidate(vec![explanation(1, 6)]));
    assert!(a.decoy);
    // Ties: deterministic, and roughly half go to the decoy.
    let tie = candidate(vec![explanation(3, 3)]);
    assert_eq!(model.assign(&tie), model.assign(&tie));
    let decoys = (0..1000)
        .filter(|scan| {
            let mut tie = tie.clone();
            tie.feature.spec_id = format!("scan={scan}");
            model.assign(&tie).decoy
        })
        .count();
    assert!((400..600).contains(&decoys), "{decoys}");
    let a = model.assign(&candidate(vec![explanation(1, 1), explanation(7, 2)]));
    assert_eq!((a.explanation, a.decoy), (1, false));
    // The competition statistic adds the lead over the best decoy.
    let (low, high) = (explanation(1, 1), explanation(7, 2));
    let best_decoy = score(&low, &low.decoy).max(score(&high, &high.decoy));
    assert_eq!(a.opposite, best_decoy);
    assert_eq!(a.competition(), 2.0 * a.score - a.opposite);
}

#[test]
fn selection_keeps_the_best_supported_candidate_per_spectrum() {
    let candidate = |scan: &str| GlycoCandidate {
        feature: sage_core::scoring::Feature {
            spec_id: scan.into(),
            ..Default::default()
        },
        peptide_mass: 1000.0,
        oxonium_ions: 3,
        oxonium_fraction: 0.1,
        anchored: true,
        y1: true,
        hex_ratio: 0.0,
        random_y: 0.1,
        random_oxonium: 0.1,
        explanations: vec![explanation(1, 0)],
        site: Default::default(),
        twin_hyperscore: None,
        ladder: Default::default(),
        siblings: 0,
    };
    let assignment = |score: f64| Assignment {
        explanation: 0,
        decoy: false,
        score,
        runner_up: f64::NEG_INFINITY,
        opposite: f64::NEG_INFINITY,
    };
    let candidates = [
        candidate("a"),
        candidate("a"),
        candidate("a"),
        candidate("b"),
        candidate("c"),
        candidate("c"),
    ];
    let assignments: Vec<_> = [1.0, 5.0, 3.0, 2.0, 4.0, 4.0]
        .into_iter()
        .map(assignment)
        .collect();
    let selected = select(&candidates, &assignments);
    assert_eq!(
        selected,
        vec![
            Selection {
                index: 1,
                lead: 2.0
            },
            // A lone candidate leads over zero.
            Selection {
                index: 3,
                lead: 2.0
            },
            // Ties keep the better-ranked (earlier) peptide.
            Selection {
                index: 4,
                lead: 0.0
            },
        ]
    );
    // Every candidate gets the lead its selection would have.
    let leads = glycan_leads(&candidates, &assignments);
    for s in &selected {
        assert_eq!(leads[s.index], s.lead);
    }
    let groups: Vec<_> = spectrum_groups(&candidates)
        .into_iter()
        .map(|g| (g.start, g.end))
        .collect();
    assert_eq!(groups, vec![(0, 3), (3, 4), (4, 6)]);

    // Siblings: other spectra whose pick is the same peptide.
    let mut candidates = candidates.to_vec();
    for (i, (candidate, peptide)) in candidates.iter_mut().zip([7, 8, 9, 8, 8, 7]).enumerate() {
        candidate.feature.peptide_idx = sage_core::database::PeptideIx(peptide);
        candidate.feature.expmass = 1000.0 + 100.0 * i as f32;
        candidate.feature.rt = i as f32 * 5.0;
    }
    let chosen: Vec<usize> = selected.iter().map(|s| s.index).collect();
    count_siblings(&mut candidates, &chosen, 0.0);
    let siblings: Vec<u16> = candidates.iter().map(|c| c.siblings).collect();
    // Picks 1, 3 and 4 are all peptide 8, at three different precursors.
    assert_eq!(siblings, vec![0, 2, 0, 2, 2, 0]);
    // Within 6 min, pick 1 (5 min) is alone; picks 3 and 4 (15, 20) pair.
    count_siblings(&mut candidates, &chosen, 6.0);
    let siblings: Vec<u16> = candidates.iter().map(|c| c.siblings).collect();
    assert_eq!(siblings, vec![0, 0, 0, 1, 1, 0]);
    // A repeat scan of pick 3's precursor counts once.
    candidates[4].feature.expmass = candidates[3].feature.expmass + 0.005;
    candidates[4].feature.rt = candidates[3].feature.rt + 0.5;
    count_siblings(&mut candidates, &chosen, 0.0);
    let siblings: Vec<u16> = candidates.iter().map(|c| c.siblings).collect();
    assert_eq!(siblings, vec![0, 1, 0, 1, 1, 0]);
}
