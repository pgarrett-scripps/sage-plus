//! Scoping probe for the rustyms successor crates (`mzcore`, `mzannotate`).
//! Not part of any release; exists so the size and link cost can be measured.

use mzcore::{
    chemistry::Molecule,
    ontology::{Ontology, STATIC_ONTOLOGIES},
    prelude::*,
    system::dalton,
};

/// Parse a ProForma string and return its monoisotopic neutral masses (Da).
pub fn proforma_masses(proforma: &str) -> Result<Vec<f64>, String> {
    let (peptidoform, _warnings) = PeptidoformIonSet::pro_forma(proforma, &STATIC_ONTOLOGIES)
        .map_err(|errors| {
            errors
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<_>>()
                .join("; ")
        })?;
    Ok(peptidoform
        .monoisotopic_masses()
        .iter()
        .map(|m| m.mass().get::<dalton>())
        .collect())
}

/// Monoisotopic delta mass of a Unimod modification looked up by name.
pub fn unimod_mass(name: &str) -> Option<f64> {
    STATIC_ONTOLOGIES
        .get_by_name(&[Ontology::Unimod], name)
        .map(|m| m.monoisotopic_mass().mass().get::<dalton>())
}

/// Probe entry point: parse, look up, print.
pub fn probe(proforma: &str, mod_name: &str) -> String {
    format!(
        "proforma {proforma}: {:?}\nunimod {mod_name}: {:?}",
        proforma_masses(proforma),
        unimod_mass(mod_name)
    )
}

/// Number of singly charged b/y fragments mzannotate generates for `proforma`.
#[cfg(feature = "rustyms-fragments")]
pub fn by_fragment_count(proforma: &str) -> Result<usize, String> {
    use mzannotate::{
        annotation::model::{FragmentationModel, PrimaryIonSeries},
        prelude::*,
    };
    use mzcore::{
        chemistry::OutputMonoIsotopic,
        system::{e, isize::Charge},
    };
    let (peptidoform, _) = PeptidoformIonSet::pro_forma(proforma, &STATIC_ONTOLOGIES)
        .map_err(|errors| errors[0].to_string())?;
    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default());
    Ok(peptidoform
        .generate_theoretical_fragments::<OutputMonoIsotopic>(Charge::new::<e>(1), &model)
        .len())
}
