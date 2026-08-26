use sage_core::database::Builder;
use sage_core::fasta::Fasta;
use sage_core::modification::StaticModEntry;
use std::collections::HashMap;
use std::fmt::Write;
use std::time::Instant;

fn sequence(mut value: usize, len: usize) -> String {
    const RESIDUES: &[u8] = b"DEFGHILMNPQSTVWY";
    let mut sequence = vec![b'A'; len];
    sequence[1] = b'C';
    for residue in &mut sequence[2..len - 1] {
        *residue = RESIDUES[value % RESIDUES.len()];
        value /= RESIDUES.len();
    }
    sequence[len - 1] = b'K';
    String::from_utf8(sequence).unwrap()
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
    let modification_mode = std::env::args().nth(3).unwrap_or_else(|| "normal".into());

    let mut fasta_text = String::with_capacity(count * (peptide_len + 16));
    for index in 0..count {
        writeln!(fasta_text, ">benchmark_{index}").unwrap();
        writeln!(fasta_text, "{}", sequence(index, peptide_len)).unwrap();
    }
    let fasta = Fasta::parse(fasta_text, "rev_", true).unwrap();

    let residues = match modification_mode.as_str() {
        "normal" => "AC",
        "heavy" => "ACDEFGHIKLMNPQRSTVWY",
        mode => panic!("unknown modification mode {mode}"),
    };
    let static_mods = residues
        .chars()
        .map(|residue| (residue.to_string(), StaticModEntry::Mass(1.0)))
        .collect::<HashMap<_, _>>();
    let parameters = Builder {
        static_mods: Some(static_mods),
        generate_decoys: Some(false),
        ..Builder::default()
    }
    .make_parameters();

    let digest_started = Instant::now();
    let peptides = parameters.digest(&fasta);
    let digest_ms = digest_started.elapsed().as_millis();
    let modification_entries = peptides
        .iter()
        .map(|peptide| peptide.modifications.len())
        .sum::<usize>();
    let spilled_collections = peptides
        .iter()
        .filter(|peptide| peptide.modifications.spilled())
        .count();

    let index_started = Instant::now();
    let database = parameters.build_from_peptides(peptides);
    println!(
        "mode={modification_mode} peptides={} modifications={modification_entries} spilled={spilled_collections} fragments={} digest_ms={digest_ms} index_ms={}",
        database.peptides.len(),
        database.fragments.len(),
        index_started.elapsed().as_millis()
    );
}
