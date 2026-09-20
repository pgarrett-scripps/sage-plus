//! Synthetic regression benchmark for positional modification search and localization.
use sage_core::database::Builder;
use sage_core::enzyme::{group_digests, Digest, Position};
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::{Tolerance, PROTON};
use sage_core::modification::ModificationDefinition;
use sage_core::peptide::{Peptide, Site};
use sage_core::scoring::{ScoreType, Scorer};
use sage_core::spectrum::{Precursor, ProcessedSpectrum};
use serde_json::json;
use std::sync::Arc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut results = Vec::new();
    for (name, sequence, position, site, keys) in [
        (
            "H3K9_first",
            "KSTGGKAPR",
            Position::Internal,
            0,
            vec!["^K", "~K"],
        ),
        (
            "H3_internal",
            "KSTGGKAPR",
            Position::Internal,
            5,
            vec!["~K"],
        ),
        (
            "protein_C",
            "ASQKSTGGK",
            Position::Cterm,
            8,
            vec!["^K", "~K", "]K"],
        ),
        (
            "peptide_last",
            "ASQKSTGGK",
            Position::Internal,
            8,
            vec!["~K", "$K"],
        ),
    ] {
        let digest = Digest {
            sequence: sequence.into(),
            position,
            protein: name.into(),
            ..Default::default()
        };
        let definition = Arc::new(ModificationDefinition {
            name: Some("Acetyl".into()),
            ..ModificationDefinition::bare(42.0106)
        });
        let truth =
            Peptide::try_from(digest.clone())?.with_mass_offset(Site::Sequence(site), &definition);
        let mut masses = [Kind::B, Kind::Y]
            .into_iter()
            .flat_map(|kind| IonSeries::new(&truth, kind))
            .map(|ion| ion.monoisotopic_mass)
            .collect::<Vec<_>>();
        masses.sort_by(f32::total_cmp);
        let spectrum = ProcessedSpectrum {
            level: 2,
            id: name.into(),
            precursors: vec![Precursor {
                mz: truth.monoisotopic / 2.0 + PROTON,
                charge: Some(2),
                ..Default::default()
            }],
            intensities: vec![100.0; masses.len()],
            charges: vec![1; masses.len()],
            total_ion_current: 100.0 * masses.len() as f32,
            masses,
            ..Default::default()
        };
        for mode in ["database", "mass_offset"] {
            let mods = keys
                .iter()
                .map(|key| {
                    (
                        (*key).into(),
                        json!([{
                            "mass":42.0106, "name":"Acetyl", "max_count":1, "search_mode":mode
                        }]),
                    )
                })
                .collect::<serde_json::Map<_, _>>();
            let parameters = serde_json::from_value::<Builder>(json!({
                "variable_mods":mods, "generate_decoys":false, "max_variable_mods":1,
                "peptide_min_mass":0
            }))?
            .make_parameters();
            let peptides = parameters.modify_digests(group_digests(vec![digest.clone()]));
            let database = parameters.build_from_peptides(peptides);
            let scorer = Scorer {
                db: &database,
                precursor_tol: Tolerance::Ppm(-10.0, 10.0),
                fragment_tol: Tolerance::Ppm(-10.0, 10.0),
                min_matched_peaks: 4,
                min_isotope_err: 0,
                max_isotope_err: 0,
                min_precursor_charge: 2,
                max_precursor_charge: 2,
                override_precursor_charge: false,
                max_fragment_charge: Some(1),
                chimera: false,
                report_psms: 1,
                wide_window: false,
                annotate_matches: false,
                mass_shift_ppm: 20.0,
                score_type: ScoreType::SageHyperScore,
            };
            let hits = scorer.score(&spectrum);
            let hit = hits.first().expect("synthetic truth must be identified");
            let observed = database.resolve_peptide(hit);
            assert_eq!(observed.to_string(), truth.to_string(), "{name} {mode}");
            let localization = sage_core::ptm::localize(
                &observed,
                &spectrum,
                &[Kind::B, Kind::Y],
                &database.localization_mods,
                Tolerance::Ppm(-10.0, 10.0),
                Some(1),
                2,
            );
            let localized = &localization.mods[0];
            assert_eq!(
                localized.best_sites[0].position, site as usize,
                "{name} {mode}"
            );
            results.push(json!({"case":name,"mode":mode,"peptide":observed.to_string(),
                "site":site+1,"matched_peaks":hit.matched_peaks,"candidate_sites":localized.candidate_sites,
                "indexed_peptides":database.peptides.len(),"passed":true}));
        }
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "version":env!("CARGO_PKG_VERSION"),"scope":"synthetic regression, not empirical calibration", "results":results
        }))?
    );
    Ok(())
}
