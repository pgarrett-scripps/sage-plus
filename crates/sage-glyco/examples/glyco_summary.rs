//! Summarize a `glyco.sage.parquet` file.
//!
//! `cargo run --release -p sage-glyco --example glyco_summary -- glyco.sage.parquet [rows.tsv]`
//!
//! Prints counts at the file's FDR thresholds and the share of passing
//! glycoPSMs that are not high-mannose/paucimannose (HexNAc(2)Hex(n)), which
//! is the error proxy for yeast samples, whose N-glycans are all
//! HexNAc(2)Hex(n). With a second argument, every row is also written
//! there as tab-separated text for ad hoc analysis.

use std::collections::{BTreeMap, HashSet};

use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::Field;
use sage_glyco::composition::{GlycanComposition, Monosaccharide};

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .ok_or_else(|| anyhow::anyhow!("usage: glyco_summary FILE"))?;
    let reader = SerializedFileReader::new(std::fs::File::open(&path)?)?;
    let mut tsv = match std::env::args().nth(2) {
        Some(out) => Some(std::io::BufWriter::new(std::fs::File::create(out)?)),
        None => None,
    };
    if let Some(tsv) = tsv.as_mut() {
        use std::io::Write;
        let names: Vec<_> = reader
            .metadata()
            .file_metadata()
            .schema_descr()
            .columns()
            .iter()
            .map(|c| c.name().to_string())
            .collect();
        writeln!(tsv, "{}", names.join("\t"))?;
    }
    let mut rows = 0usize;
    let mut peptide_pass = 0usize;
    let mut glycan_decoys = 0usize;
    let mut passing = 0usize;
    let mut non_hm = 0usize;
    let mut isotope1 = 0usize;
    let mut adduct = 0usize;
    let mut ambiguous = 0usize;
    let mut unique = HashSet::new();
    let mut peptides = HashSet::new();
    let mut compositions: BTreeMap<String, usize> = BTreeMap::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        rows += 1;
        if let Some(tsv) = tsv.as_mut() {
            use std::io::Write;
            let values: Vec<_> = row
                .get_column_iter()
                .map(|(_, field)| match field {
                    Field::Str(v) => v.clone(),
                    other => other.to_string(),
                })
                .collect();
            writeln!(tsv, "{}", values.join("\t"))?;
        }
        let mut get = std::collections::HashMap::new();
        for (name, field) in row.get_column_iter() {
            get.insert(name.as_str(), field);
        }
        let int = |name: &str| match get[name] {
            Field::Int(v) => *v as i64,
            Field::Long(v) => *v,
            _ => 0,
        };
        let float = |name: &str| match get[name] {
            Field::Float(v) => *v as f64,
            Field::Double(v) => *v,
            _ => f64::NAN,
        };
        let text = |name: &str| match get[name] {
            Field::Str(v) => v.clone(),
            _ => String::new(),
        };
        let boolean = |name: &str| matches!(get[name], Field::Bool(true));
        if int("label") == 1 && float("peptide_q") <= 0.01 {
            peptide_pass += 1;
            glycan_decoys += usize::from(boolean("glycan_decoy"));
        }
        if !boolean("passes_fdr") {
            continue;
        }
        passing += 1;
        let glycan = text("glycan");
        let composition = GlycanComposition::parse(&glycan).map_err(anyhow::Error::msg)?;
        let high_mannose = composition.count(Monosaccharide::HexNAc) == 2
            && composition.count(Monosaccharide::Fuc) == 0
            && composition.count(Monosaccharide::NeuAc) == 0
            && composition.count(Monosaccharide::NeuGc) == 0;
        non_hm += usize::from(!high_mannose);
        isotope1 += usize::from(int("isotope_error") != 0);
        adduct += usize::from(int("ammonium_adducts") != 0);
        ambiguous += usize::from(int("explanations") > 1);
        unique.insert((text("stripped_peptide"), glycan.clone()));
        peptides.insert(text("stripped_peptide"));
        *compositions.entry(glycan).or_default() += 1;
    }
    let pct = |n: usize| 100.0 * n as f64 / passing.max(1) as f64;
    println!("rows\t{rows}");
    println!("target_peptide_q_1pct\t{peptide_pass}");
    println!("glycan_decoy_winners\t{glycan_decoys}");
    println!("glyco_psms_passing\t{passing}");
    println!("unique_peptide_glycan\t{}", unique.len());
    println!("unique_peptides\t{}", peptides.len());
    println!("non_high_mannose\t{non_hm} ({:.1}%)", pct(non_hm));
    println!("isotope_error_nonzero\t{isotope1} ({:.1}%)", pct(isotope1));
    println!("ammonium_adduct\t{adduct} ({:.1}%)", pct(adduct));
    println!("ambiguous_mass\t{ambiguous} ({:.1}%)", pct(ambiguous));
    let mut top: Vec<_> = compositions.into_iter().collect();
    top.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    println!("top_compositions");
    for (glycan, n) in top.iter().take(15) {
        println!("  {glycan}\t{n}");
    }
    Ok(())
}
