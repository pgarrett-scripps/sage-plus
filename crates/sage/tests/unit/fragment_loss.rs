use super::*;
use crate::enzyme::Digest;
use crate::ion_series::{IonGroupSeries, Kind};

const WATER: f32 = 18.010565;
const AMMONIA: f32 = 17.026549;

fn peptide(s: &str) -> Peptide {
    Peptide::try_from(Digest {
        sequence: s.into(),
        ..Default::default()
    })
    .unwrap()
}

fn entry(mass: f32, sites: &[&str], kinds: &[Kind]) -> FragmentLossEntry {
    FragmentLossEntry {
        mass,
        sites: sites.iter().map(|s| s.to_string()).collect(),
        ion_kinds: kinds.to_vec(),
        allow_modified: false,
    }
}

fn approved() -> BTreeMap<String, FragmentLossEntry> {
    serde_json::from_value(serde_json::json!({
        "Water":   { "mass": 18.010565, "sites": ["S","T","E","D"], "ion_kinds": ["b","y"], "allow_modified": false },
        "Ammonia": { "mass": 17.026549, "sites": ["R","K","N","Q"], "ion_kinds": ["y"] }
    }))
    .unwrap()
}

fn losses(entries: BTreeMap<String, FragmentLossEntry>, max: Option<usize>) -> FragmentLosses {
    resolve(Some(&entries), max, &[Kind::B, Kind::Y])
        .unwrap()
        .unwrap()
}

/// (series index, sorted generic-loss totals) for every group of `kind`.
fn loss_forms(peptide: &Peptide, kind: Kind, losses: &FragmentLosses) -> Vec<(usize, Vec<f32>)> {
    IonGroupSeries::with_fragment_losses(peptide, kind, Some(losses))
        .map(|group| {
            let mut totals = group
                .variants
                .iter()
                .filter(|variant| variant.fragment_loss)
                .map(|variant| variant.neutral_loss.unwrap())
                .collect::<Vec<_>>();
            totals.sort_by(f32::total_cmp);
            (group.series_index, totals)
        })
        .collect()
}

#[test]
fn approved_config_resolves() {
    let resolved = losses(approved(), None);
    assert_eq!(resolved.max_losses, 1);
    let ammonia = &resolved.losses[0];
    assert_eq!(&*ammonia.name, "Ammonia");
    assert_eq!(ammonia.ion_kinds, vec![Kind::Y]);
    assert!(!ammonia.allow_modified);
    assert_eq!(resolved.losses[1].sites.len(), 4);
}

#[test]
fn absent_key_is_off() {
    assert_eq!(resolve(None, None, &[Kind::B, Kind::Y]), Ok(None));
}

#[test]
fn invalid_configurations_are_rejected() {
    let kinds = [Kind::B, Kind::Y];
    let one = |e: FragmentLossEntry| BTreeMap::from([("Loss".to_string(), e)]);
    let error = |map: BTreeMap<String, FragmentLossEntry>, max: Option<usize>| {
        resolve(Some(&map), max, &kinds).unwrap_err()
    };

    assert!(resolve(None, Some(1), &kinds)
        .unwrap_err()
        .contains("requires `database.fragment_losses`"));
    assert!(error(BTreeMap::new(), None).contains("at least one loss"));
    assert!(error(one(entry(0.0, &["S"], &[Kind::B])), None).contains("positive, finite"));
    assert!(error(one(entry(-18.0, &["S"], &[Kind::B])), None).contains("positive, finite"));
    assert!(error(one(entry(f32::NAN, &["S"], &[Kind::B])), None).contains("positive, finite"));
    assert!(error(one(entry(f32::INFINITY, &["S"], &[Kind::B])), None).contains("positive"));
    assert!(error(one(entry(18.0, &[], &[Kind::B])), None).contains("no sites"));
    assert!(error(one(entry(18.0, &["B"], &[Kind::B])), None).contains("invalid site `B`"));
    assert!(error(one(entry(18.0, &["H2O"], &[Kind::B])), None).contains("invalid site"));
    assert!(error(one(entry(18.0, &["^E"], &[Kind::B])), None)
        .contains("use explicit site `first_residue:E`"));
    assert!(error(one(entry(18.0, &["motif:N[^P][ST]"], &[Kind::B])), None).contains("motif"));
    assert!(error(one(entry(18.0, &["S"], &[])), None).contains("no `ion_kinds`"));
    assert!(error(one(entry(18.0, &["S"], &[Kind::C])), None)
        .contains("ion kind `c`, which is not in `database.ion_kinds`"));
    assert!(error(one(entry(18.0, &["S"], &[Kind::B])), Some(0)).contains("at least 1"));
    let named = |name: &str| BTreeMap::from([(name.to_string(), entry(18.0, &["S"], &[Kind::B]))]);
    assert!(error(named(" Water"), None).contains("surrounding whitespace"));
    assert!(error(named(""), None).contains("nonempty"));
}

#[test]
fn unknown_fields_are_rejected() {
    let result = serde_json::from_value::<BTreeMap<String, FragmentLossEntry>>(
        serde_json::json!({"Water": {"mass": 18.0, "sites": ["S"], "ion_kinds": ["b"], "formula": "H2O"}}),
    );
    assert!(result.is_err());
}

#[test]
fn fragments_carry_losses_only_when_they_contain_a_site() {
    let losses = losses(approved(), None);
    // P E P T I D E K: E at 1, T at 3, D at 5, E at 6, K at 7.
    let peptide = peptide("PEPTIDEK");
    let b = loss_forms(&peptide, Kind::B, &losses);
    // b1 = P has no water site; b2 onwards contains E.
    assert_eq!(b[0], (0, vec![]));
    for (_, totals) in &b[1..] {
        assert_eq!(totals, &vec![WATER]);
    }
    let y = loss_forms(&peptide, Kind::Y, &losses);
    // y1 = K: ammonia only. y2 = EK onwards: ammonia or water.
    assert_eq!(y.last().unwrap(), &(6, vec![AMMONIA]));
    for (_, totals) in &y[..6] {
        assert_eq!(totals, &vec![AMMONIA, WATER]);
    }
}

#[test]
fn ion_kinds_restrict_losses() {
    let losses = losses(approved(), None);
    // No water sites: only ammonia, and only on y ions.
    let peptide = peptide("PAGLAK");
    assert!(loss_forms(&peptide, Kind::B, &losses)
        .iter()
        .all(|(_, totals)| totals.is_empty()));
    assert!(loss_forms(&peptide, Kind::Y, &losses)
        .iter()
        .all(|(_, totals)| totals == &vec![AMMONIA]));
}

#[test]
fn loss_masses_are_subtracted_from_the_retained_ion() {
    let losses = losses(approved(), None);
    let peptide = peptide("PEPTIDEK");
    for group in IonGroupSeries::with_fragment_losses(&peptide, Kind::Y, Some(&losses)) {
        let retained = group.variants[0].monoisotopic_mass;
        assert!(!group.variants[0].fragment_loss);
        for variant in group.variants.iter().filter(|v| v.fragment_loss) {
            let loss = variant.neutral_loss.unwrap();
            assert!((retained - loss - variant.monoisotopic_mass).abs() < 1e-3);
        }
    }
}

#[test]
fn stacking_is_capped_and_limited_by_site_count() {
    let water = BTreeMap::from([("Water".to_string(), entry(WATER, &["S"], &[Kind::B]))]);
    let peptide = peptide("ASGSGK");
    let one = losses(water.clone(), None);
    let two = losses(water.clone(), Some(2));
    let three = losses(water, Some(3));
    // b2 = AS has one S: one water at most, whatever the cap.
    assert_eq!(loss_forms(&peptide, Kind::B, &two)[1].1, vec![WATER]);
    // b4 = ASGS has two S.
    assert_eq!(loss_forms(&peptide, Kind::B, &one)[3].1, vec![WATER]);
    assert_eq!(
        loss_forms(&peptide, Kind::B, &two)[3].1,
        vec![WATER, 2.0 * WATER]
    );
    assert_eq!(
        loss_forms(&peptide, Kind::B, &three)[3].1,
        vec![WATER, 2.0 * WATER]
    );
}

#[test]
fn stacking_combines_different_losses() {
    let map = BTreeMap::from([
        ("Water".to_string(), entry(WATER, &["E"], &[Kind::Y])),
        ("Ammonia".to_string(), entry(AMMONIA, &["K"], &[Kind::Y])),
    ]);
    let losses = losses(map, Some(2));
    let peptide = peptide("GGEK");
    let y = loss_forms(&peptide, Kind::Y, &losses);
    // y2 = EK: water, ammonia, or both.
    let mut expected = vec![WATER, AMMONIA, WATER + AMMONIA];
    expected.sort_by(f32::total_cmp);
    assert_eq!(y[1].1, expected);
}

#[test]
fn modified_residues_are_sites_only_when_allowed() {
    let static_mods = [(ModificationSpecificity::Residue(b'S'), 79.966_33)].into();
    let modified = peptide("ASGK").apply(&[], &static_mods, 1, None).remove(0);
    let mut water = entry(WATER, &["S"], &[Kind::B]);
    let blocked = losses(BTreeMap::from([("Water".to_string(), water.clone())]), None);
    assert!(loss_forms(&modified, Kind::B, &blocked)
        .iter()
        .all(|(_, totals)| totals.is_empty()));
    water.allow_modified = true;
    let allowed = losses(BTreeMap::from([("Water".to_string(), water)]), None);
    assert_eq!(loss_forms(&modified, Kind::B, &allowed)[1].1, vec![WATER]);
}

#[test]
fn terminal_sites_select_the_terminal_fragments() {
    let map = BTreeMap::from([(
        "Water".to_string(),
        entry(WATER, &["peptide_c_term"], &[Kind::B, Kind::Y]),
    )]);
    let losses = losses(map, None);
    let peptide = peptide("GAGK");
    assert!(loss_forms(&peptide, Kind::B, &losses)
        .iter()
        .all(|(_, totals)| totals.is_empty()));
    assert!(loss_forms(&peptide, Kind::Y, &losses)
        .iter()
        .all(|(_, totals)| totals == &vec![WATER]));
}

/// Generic losses stack on each modification neutral-loss form of the group,
/// and a required modification loss removes the generic-only form too.
#[test]
fn generic_losses_combine_with_modification_losses() {
    use crate::modification::{ModificationDefinition, NeutralLossMode};
    use std::collections::HashMap;
    const MOD_LOSS: f32 = 10.0;
    let with_mod = |mode| {
        let modification = Arc::new(ModificationDefinition {
            mass: 20.0,
            name: Some(Arc::from("TestMod")),
            neutral_losses: Arc::from([MOD_LOSS]),
            site_losses: None,
            neutral_loss_mode: mode,
            channel_offsets: Arc::default(),
        });
        peptide("AMEK")
            .apply(
                &[(
                    ModificationSpecificity::Residue(b'M'),
                    modification,
                    Some(1),
                )],
                &HashMap::default(),
                1,
                None,
            )
            .into_iter()
            .find(|peptide| peptide.to_string().contains("TestMod"))
            .unwrap()
    };
    let water = losses(
        BTreeMap::from([("Water".to_string(), entry(WATER, &["E"], &[Kind::B]))]),
        None,
    );
    let totals = |groups: Vec<crate::ion_series::IonGroup>, generic: bool| {
        groups
            .iter()
            .map(|group| {
                group
                    .variants
                    .iter()
                    .filter(|variant| variant.fragment_loss == generic)
                    .map(|variant| variant.neutral_loss.unwrap_or(0.0))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };

    // b1 = A, b2 = AM (modified, no water site), b3 = AME.
    let optional = with_mod(NeutralLossMode::Optional);
    let forms = loss_forms(&optional, Kind::B, &water);
    assert_eq!(forms[0].1, Vec::<f32>::new());
    assert_eq!(forms[1].1, Vec::<f32>::new());
    assert_eq!(forms[2].1, vec![WATER, WATER + MOD_LOSS]);
    // The modification-loss forms are the same as without generic losses.
    assert_eq!(
        totals(
            IonGroupSeries::with_fragment_losses(&optional, Kind::B, Some(&water)).collect(),
            false
        ),
        totals(IonGroupSeries::new(&optional, Kind::B).collect(), false)
    );
    assert_eq!(
        totals(IonGroupSeries::new(&optional, Kind::B).collect(), false)[2],
        vec![0.0, MOD_LOSS]
    );

    // Required: only the forms that carry the modification loss.
    let required = with_mod(NeutralLossMode::Required);
    assert_eq!(
        loss_forms(&required, Kind::B, &water)[2].1,
        vec![WATER + MOD_LOSS]
    );
}

#[test]
fn without_losses_the_groups_are_unchanged() {
    let peptide = peptide("PEPTIDEK");
    for kind in [Kind::B, Kind::Y] {
        let plain = IonGroupSeries::new(&peptide, kind).collect::<Vec<_>>();
        let none = IonGroupSeries::with_fragment_losses(&peptide, kind, None).collect::<Vec<_>>();
        assert_eq!(plain.len(), none.len());
        for (a, b) in plain.iter().zip(&none) {
            assert_eq!(a.variants.len(), b.variants.len());
            for (x, y) in a.variants.iter().zip(&b.variants) {
                assert_eq!(x.monoisotopic_mass.to_bits(), y.monoisotopic_mass.to_bits());
                assert!(!y.fragment_loss);
            }
        }
    }
}
