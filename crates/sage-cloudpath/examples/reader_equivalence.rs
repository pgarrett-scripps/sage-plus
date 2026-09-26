use std::io::{BufWriter, Write};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let path = args.next().ok_or("input path required")?;
    let output = args.next().ok_or("output path required")?;
    let url = sage_cloudpath::to_url(&path)?;
    let spectra = sage_cloudpath::util::read_spectra(&url, 0, None, Default::default(), true)?;
    let mut out = BufWriter::new(std::fs::File::create(output)?);
    let mut peaks = 0usize;
    let mut levels = std::collections::BTreeMap::new();
    for spectrum in &spectra {
        *levels.entry(spectrum.ms_level).or_insert(0usize) += 1;
        peaks += spectrum.mz.len();
        writeln!(out, "{spectrum:?}")?;
    }
    out.flush()?;
    println!(
        "{}",
        serde_json::json!({"spectra":spectra.len(),"peaks":peaks,"levels":levels})
    );
    Ok(())
}
