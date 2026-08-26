use sage_core::database::{Builder, Parameters};
use sage_core::enzyme::Position;
use sage_core::peptide::Peptide;
use std::sync::Arc;
use std::time::Instant;

fn sequence(mut value: usize, len: usize) -> Arc<[u8]> {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut sequence = vec![b'A'; len];
    for residue in &mut sequence {
        *residue = RESIDUES[value % RESIDUES.len()];
        value /= RESIDUES.len();
    }
    Arc::from(sequence.into_boxed_slice())
}

fn main() {
    let count = std::env::args()
        .nth(1)
        .and_then(|value| value.parse().ok())
        .unwrap_or(250_000);
    let peptide_len = std::env::args()
        .nth(2)
        .and_then(|value| value.parse().ok())
        .unwrap_or(20);

    let peptides = (0..count)
        .map(|idx| Peptide {
            sequence: sequence(idx, peptide_len),
            modifications: Vec::new(),
            monoisotopic: 1_000.0 + idx as f32 / count as f32 * 1_000.0,
            position: Position::Internal,
            proteins: vec![Arc::from("benchmark")].into(),
            ..Peptide::default()
        })
        .collect::<Vec<_>>();

    let parameters: Parameters = Builder::default().make_parameters();
    let started = Instant::now();
    let database = parameters.build_from_peptides(peptides);
    println!(
        "peptides={} fragments={} elapsed_ms={}",
        database.peptides.len(),
        database.fragments.len(),
        started.elapsed().as_millis()
    );
}
