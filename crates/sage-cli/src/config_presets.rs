//! Starting configurations printed by `sage --write-config <NAME>`.

/// Name, one-line description and JSON text of each preset.
pub const PRESETS: [(&str, &str, &str); 7] = [
    (
        "minimal",
        "Only the required settings: trypsin, carbamidomethyl C, b/y ions, ±10/±20 ppm",
        include_str!("../presets/minimal.json"),
    ),
    (
        "full",
        "Every option, set to the value Sage uses when it is left out",
        include_str!("../presets/full.json"),
    ),
    (
        "trypsin-hcd",
        "Tryptic HCD search with variable Met oxidation",
        include_str!("../presets/trypsin-hcd.json"),
    ),
    (
        "trypsin-hcd-tmt",
        "Tryptic HCD search with TMTpro on K and peptide N-termini and MS2 TMTpro 18-plex quant",
        include_str!("../presets/trypsin-hcd-tmt.json"),
    ),
    (
        "phospho",
        "Tryptic HCD search with variable STY phosphorylation and site localization",
        include_str!("../presets/phospho.json"),
    ),
    (
        "etd",
        "Tryptic ETD search with c and z-dot ions (add b and y for EThcD)",
        include_str!("../presets/etd.json"),
    ),
    (
        "nonspecific",
        "Non-specific digest of 8-15 residues, no alkylation, with the prefilter",
        include_str!("../presets/nonspecific.json"),
    ),
];

/// The JSON text of a preset.
pub fn preset(name: &str) -> Option<&'static str> {
    PRESETS
        .iter()
        .find(|(preset, _, _)| *preset == name)
        .map(|(_, _, json)| *json)
}

/// The preset names with their descriptions, one per line.
pub fn list() -> String {
    let width = PRESETS
        .iter()
        .map(|(name, _, _)| name.len())
        .max()
        .unwrap_or(0);
    PRESETS
        .iter()
        .map(|(name, description, _)| format!("{name:width$}  {description}\n"))
        .collect()
}

/// Write a preset to `path`, or to stdout when `path` is `None`. An existing
/// file is replaced only with `overwrite`.
pub fn write(name: &str, path: Option<&str>, overwrite: bool) -> anyhow::Result<()> {
    if name == "list" {
        print!("{}", list());
        return Ok(());
    }
    let Some(json) = preset(name) else {
        anyhow::bail!(
            "unknown configuration preset `{name}`; choose one of:\n{}",
            list()
        );
    };
    match path {
        None => print!("{json}"),
        Some(path) => {
            let mut file = std::fs::OpenOptions::new();
            file.write(true);
            if overwrite {
                file.create(true).truncate(true);
            } else {
                file.create_new(true);
            }
            let mut file = file.open(path).map_err(|error| match error.kind() {
                std::io::ErrorKind::AlreadyExists => {
                    anyhow::anyhow!("{path} already exists; pass --overwrite to replace it")
                }
                _ => anyhow::Error::new(error).context(format!("cannot write {path}")),
            })?;
            std::io::Write::write_all(&mut file, json.as_bytes())?;
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "../tests/unit/config_presets.rs"]
mod tests;
