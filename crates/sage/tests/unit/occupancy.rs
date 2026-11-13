use super::*;
use crate::database::Builder;
use crate::enzyme::Digest;
use crate::modification::ModificationDefinition;
use crate::peptide::{AppliedModification, CompactModifications, ModificationKind, Site};

fn parameters(database: serde_json::Value) -> Parameters {
    serde_json::from_value::<Builder>(database)
        .unwrap()
        .make_parameters()
        .expect("valid test configuration")
}

fn peptide(sequence: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: sequence.into(),
        ..Default::default()
    })
    .unwrap()
}

/// `peptide` carrying `mass` at each site.
fn modified(sequence: &str, mass: f32, sites: &[Site]) -> Peptide {
    let mut peptide = peptide(sequence);
    let definition = Arc::new(ModificationDefinition::bare(mass));
    peptide.modifications =
        CompactModifications::from_applied(sites.iter().map(|&site| AppliedModification {
            site,
            modification: definition.clone(),
            kind: ModificationKind::Ordinary,
        }))
        .unwrap();
    for site in sites {
        match site {
            Site::Nterm => peptide.nterm = Some(mass),
            Site::Cterm => peptide.cterm = Some(mass),
            Site::Sequence(_) => (),
        }
    }
    peptide
}

const OXIDATION: f32 = 15.9949;
const ACETYL: f32 = 42.010565;

fn oxidation_and_n_term_acetyl() -> Parameters {
    parameters(serde_json::json!({
        "variable_mods": {"M": [OXIDATION], "^": [ACETYL]}
    }))
}

#[test]
fn every_searched_definition_and_site_rule_gets_a_rule() {
    let rules = rules(&oxidation_and_n_term_acetyl());
    let described = rules
        .iter()
        .map(|rule| {
            (
                rule.site_rule_name(),
                rule.mass(),
                rule.search_mode_name(),
                rule.name().map(str::to_string),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        described,
        vec![
            ("peptide_n_term".to_string(), ACETYL, "database", None),
            ("M".to_string(), OXIDATION, "database", None),
        ]
    );
}

#[test]
fn a_definition_configured_twice_at_one_site_rule_is_one_rule() {
    let parameters = parameters(serde_json::json!({
        "variable_mods": {"M": [OXIDATION, OXIDATION]}
    }));
    assert_eq!(rules(&parameters).len(), 1);
}

#[test]
fn occupancy_counts_eligible_and_modified_sites_of_each_rule() {
    let parameters = oxidation_and_n_term_acetyl();
    let rules = rules(&parameters);
    // Two eligible M sites in each peptidoform, four in all; the second
    // peptidoform carries the oxidation on one of them. Each peptidoform has
    // one eligible peptide N-terminus, neither of them acetylated.
    let peptidoforms = [
        peptide("MKMPEPTIDEK"),
        modified("MKMPEPTIDEK", OXIDATION, &[Site::Sequence(0)]),
    ];
    let counts = count(&parameters, &rules, &peptidoforms);

    let n_term = counts[0];
    assert_eq!(rules[0].site_rule_name(), "peptide_n_term");
    assert_eq!(n_term.peptidoforms, 2);
    assert_eq!(n_term.eligible_sites, 2);
    assert_eq!(n_term.modified_sites, 0);
    assert_eq!(n_term.modified_peptidoforms, 0);
    assert_eq!(n_term.occupancy_pct(), 0.0);

    let oxidation = counts[1];
    assert_eq!(rules[1].site_rule_name(), "M");
    assert_eq!(oxidation.peptidoforms, 2);
    assert_eq!(oxidation.eligible_sites, 4);
    assert_eq!(oxidation.modified_sites, 1);
    assert_eq!(oxidation.modified_peptidoforms, 1);
    assert_eq!(oxidation.occupancy_pct(), 25.0);
}

#[test]
fn a_peptidoform_without_an_eligible_site_does_not_count() {
    let parameters = oxidation_and_n_term_acetyl();
    let rules = rules(&parameters);
    let counts = count(&parameters, &rules, &[peptide("PEPTIDEK")]);
    assert_eq!(counts[1].peptidoforms, 0);
    assert_eq!(counts[1].eligible_sites, 0);
    assert_eq!(counts[1].occupancy_pct(), 0.0);
    // The peptide N-terminus is eligible in every peptide.
    assert_eq!(counts[0].peptidoforms, 1);
    assert_eq!(counts[0].eligible_sites, 1);
}

#[test]
fn a_terminal_rule_counts_its_own_site_only() {
    let parameters = oxidation_and_n_term_acetyl();
    let rules = rules(&parameters);
    let counts = count(
        &parameters,
        &rules,
        &[modified("MPEPTIDEK", ACETYL, &[Site::Nterm])],
    );
    assert_eq!(counts[0].modified_sites, 1);
    assert_eq!(counts[0].eligible_sites, 1);
    assert_eq!(counts[0].occupancy_pct(), 100.0);
    // The acetyl sits on the N-terminus, not on the M.
    assert_eq!(counts[1].modified_sites, 0);
    assert_eq!(counts[1].eligible_sites, 1);
}

#[test]
fn a_different_modification_at_an_eligible_site_is_not_this_rule() {
    let parameters = oxidation_and_n_term_acetyl();
    let rules = rules(&parameters);
    let counts = count(
        &parameters,
        &rules,
        &[modified("MPEPTIDEK", 79.96633, &[Site::Sequence(0)])],
    );
    assert_eq!(counts[1].eligible_sites, 1);
    assert_eq!(counts[1].modified_sites, 0);
}

#[test]
fn a_mass_offset_definition_is_counted_and_named_as_one() {
    let parameters = parameters(serde_json::json!({
        "variable_mods": {
            "Oxidation": {"mass": OXIDATION, "sites": ["M"], "search_mode": "mass_offset"}
        }
    }));
    let rules = rules(&parameters);
    assert_eq!(rules[0].search_mode_name(), "mass_offset");
    assert_eq!(rules[0].name(), Some("Oxidation"));

    // A mass-offset PSM is materialized into a peptidoform carrying the
    // definition, so the same counting applies.
    let mut peptide = peptide("MPEPTIDEM");
    let definition = Arc::new(rules[0].definition.as_ref().clone());
    peptide.modifications = CompactModifications::from_applied([AppliedModification {
        site: Site::Sequence(0),
        modification: definition,
        kind: ModificationKind::Ordinary,
    }])
    .unwrap();
    let counts = count(&parameters, &rules, &[peptide]);
    assert_eq!(counts[0].eligible_sites, 2);
    assert_eq!(counts[0].modified_sites, 1);
    assert_eq!(counts[0].occupancy_pct(), 50.0);
}

#[test]
fn an_internal_residue_rule_excludes_the_terminal_residues() {
    let parameters = parameters(serde_json::json!({
        "variable_mods": {"Acetyl": {"mass": ACETYL, "sites": ["internal_residue:K"]}}
    }));
    let rules = rules(&parameters);
    assert_eq!(rules[0].site_rule_name(), "internal_residue:K");
    let counts = count(&parameters, &rules, &[peptide("KPEPKTIDEK")]);
    assert_eq!(counts[0].eligible_sites, 1);
}

#[test]
fn no_variable_modification_gives_no_rules() {
    let parameters = parameters(serde_json::json!({}));
    assert!(rules(&parameters).is_empty());
    assert!(count(&parameters, &[], &[peptide("PEPTIDEK")]).is_empty());
}
