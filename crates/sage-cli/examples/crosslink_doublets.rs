//! Exploratory harness for the MS-cleavable crosslink prototype
//! (`sage_core::crosslink`). Not part of the search.
//!
//! Usage: crosslink_doublets <file.raw> <db.fasta> <out.jsonl> [linker]
//!
//! For every MS2 spectrum with a precursor charge it writes one JSON line:
//! doublets found at the real linker spacing, the same count at decoy
//! spacings (to estimate chance doublets), pair hypotheses, and the top five
//! chain candidates for every observed or inferred chain mass.

use rayon::prelude::*;
use sage_cloudpath::thermoraw::ThermoRawReader;
use sage_core::crosslink::{
    find_doublets, pair_hypotheses, rank_chain_candidates, CleavableLinker, DSBU, DSSO,
};
use sage_core::database::Builder;
use sage_core::fasta::Fasta;
use sage_core::mass::{Tolerance, H2O, NEUTRON, PROTON};
use sage_core::spectrum::SpectrumProcessor;
use serde_json::json;
use std::io::Write;

/// Offsets added to the real doublet spacing to count chance doublets. They
/// avoid common residue and modification mass differences near 32 Da.
const DECOY_SPACING_DELTAS: [f32; 4] = [-1.9, -0.7, 0.55, 1.3];
const MIN_CHAIN_MASS: f32 = 400.0;
const TOP_K: usize = 5;
const MAX_PAIRS: usize = 10;

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let args: Vec<String> = std::env::args().collect();
    anyhow::ensure!(
        args.len() >= 4,
        "usage: crosslink_doublets <raw> <fasta> <out.jsonl> [DSSO|DSBU]"
    );
    let linker = match args.get(4).map(String::as_str) {
        Some("DSBU") => DSBU,
        _ => DSSO,
    };

    let builder: Builder = serde_json::from_value(json!({
        "bucket_size": 8192,
        "enzyme": { "missed_cleavages": 3, "min_len": 4, "max_len": 40 },
        "peptide_min_mass": 350.0,
        "peptide_max_mass": 6000.0,
        "static_mods": { "C": 57.021464 },
        "variable_mods": { "M": [15.994915] },
        "max_variable_mods": 2,
        "decoy_tag": "rev_",
        "generate_decoys": true,
    }))?;
    let parameters = builder.make_parameters();
    let fasta = Fasta::parse(std::fs::read_to_string(&args[2])?, "rev_", true)?;
    let db = parameters.build(fasta);
    eprintln!("indexed peptides: {}", db.peptides.len());

    let raw = ThermoRawReader::with_file_id(0).parse(&args[1])?;
    let processor = SpectrumProcessor::new(150, true, 0.0);
    let spectra: Vec<_> = raw
        .into_iter()
        .filter(|s| s.ms_level == 2)
        .map(|s| processor.process(s))
        .collect();
    eprintln!("MS2 spectra: {}", spectra.len());

    let precursor_tol = Tolerance::Ppm(-10.0, 10.0);
    let fragment_tol = Tolerance::Ppm(-20.0, 20.0);
    let linkable = |p: &sage_core::peptide::Peptide| p.sequence.contains(&b'K');

    let decoy_linkers: Vec<CleavableLinker> = DECOY_SPACING_DELTAS
        .iter()
        .map(|delta| CleavableLinker {
            heavy_stub: linker.heavy_stub + delta,
            ..linker
        })
        .collect();

    let lines: Vec<String> = spectra
        .par_iter()
        .filter_map(|spectrum| {
            let precursor = spectrum.precursors.first()?;
            let charge = precursor.charge?;
            let mass = (precursor.mz - PROTON) * charge as f32;
            let max_charge = charge.max(2);
            let doublets = find_doublets(spectrum, &linker, max_charge, fragment_tol);
            let decoy_counts: Vec<usize> = decoy_linkers
                .iter()
                .map(|l| find_doublets(spectrum, l, max_charge, fragment_tol).len())
                .collect();

            let mut pairs = Vec::new();
            for isotope in 0..=2 {
                let precursor_mass = mass - isotope as f32 * NEUTRON;
                for pair in pair_hypotheses(
                    &doublets,
                    precursor_mass,
                    &linker,
                    precursor_tol,
                    MIN_CHAIN_MASS,
                ) {
                    pairs.push((isotope, pair));
                }
            }
            pairs.sort_by(|a, b| {
                b.1.both_observed
                    .cmp(&a.1.both_observed)
                    .then(b.1.intensity.total_cmp(&a.1.intensity))
            });
            pairs.truncate(MAX_PAIRS);

            let lookup = |chain: f32, partner: f32| {
                let shifts = [
                    linker.light_stub,
                    linker.heavy_stub,
                    linker.heavy_stub + H2O,
                    partner + linker.crosslink_mass,
                ];
                rank_chain_candidates(
                    &db,
                    spectrum,
                    chain,
                    &shifts,
                    precursor_tol,
                    fragment_tol,
                    max_charge,
                    linkable,
                    TOP_K,
                )
                .into_iter()
                .map(|c| {
                    let p = &db[c.peptide];
                    json!({ "peptide": p.to_string(), "decoy": p.decoy, "matched": c.matched })
                })
                .collect::<Vec<_>>()
            };

            let pairs_json: Vec<_> = pairs
                .iter()
                .map(|(isotope, pair)| {
                    json!({
                        "isotope": isotope,
                        "alpha_mass": pair.alpha_mass,
                        "beta_mass": pair.beta_mass,
                        "both_observed": pair.both_observed,
                        "intensity": pair.intensity,
                        "alpha_top": lookup(pair.alpha_mass, pair.beta_mass),
                        "beta_top": lookup(pair.beta_mass, pair.alpha_mass),
                    })
                })
                .collect();

            let line = json!({
                "id": spectrum.id,
                "charge": charge,
                "precursor_mass": mass,
                "doublets": doublets.len(),
                "doublet_chain_masses": doublets.iter().map(|d| d.chain_mass).collect::<Vec<_>>(),
                "decoy_doublets": decoy_counts,
                "pairs": pairs_json,
            });
            Some(line.to_string())
        })
        .collect();

    let mut out = std::io::BufWriter::new(std::fs::File::create(&args[3])?);
    for line in lines {
        writeln!(out, "{line}")?;
    }
    Ok(())
}
