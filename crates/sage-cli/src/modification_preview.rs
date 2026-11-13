use anyhow::{ensure, Context};
use sage_core::database::{Builder, Parameters};
use sage_core::enzyme::{group_digests, Digest, Position};
use sage_core::modification::{SearchMode, SiteMode};
use sage_core::peptide::{Peptide, Site};
use serde_json::{json, Value};

fn site_value(site: Site) -> Value {
    match site {
        Site::Nterm => json!({"terminal_group": "N"}),
        Site::Cterm => json!({"terminal_group": "C"}),
        Site::Sequence(index) => json!({"position": index + 1}),
    }
}

/// Preview one peptide without loading a FASTA or spectra or writing search outputs.
pub fn preview(
    config: &str,
    sequence: &str,
    position: &str,
    limit: usize,
) -> anyhow::Result<Value> {
    preview_with_context(config, sequence, position, limit, None)
}

pub fn preview_with_context(
    config: &str,
    sequence: &str,
    position: &str,
    limit: usize,
    context: Option<(&str, u32)>,
) -> anyhow::Result<Value> {
    preview_with_flanks(config, sequence, position, limit, context, ("", ""))
}

/// Preview with known protein residues before and after the peptide. Motif
/// sites see these residues; residues beyond them are unknown and never match.
pub fn preview_with_flanks(
    config: &str,
    sequence: &str,
    position: &str,
    limit: usize,
    context: Option<(&str, u32)>,
    (before, after): (&str, &str),
) -> anyhow::Result<Value> {
    ensure!(
        before
            .bytes()
            .chain(after.bytes())
            .all(|residue| residue.is_ascii_uppercase()),
        "preview flanking residues must be uppercase one-letter codes"
    );
    ensure!(
        (1..=10000).contains(&limit),
        "preview limit must be between 1 and 10000"
    );
    let config: Value = serde_json::from_str(config).context("invalid JSON configuration")?;
    let builder: Builder = serde_json::from_value(
        config
            .get("database")
            .cloned()
            .context("`database` must be configured")?,
    )?;
    builder
        .validate_modification_keys()
        .map_err(anyhow::Error::msg)?;
    builder.validate_enzyme().map_err(anyhow::Error::msg)?;
    let mut parameters = builder.make_parameters().map_err(anyhow::Error::msg)?;
    parameters.validate_channels().map_err(anyhow::Error::msg)?;
    parameters
        .validate_compact_modifications()
        .map_err(anyhow::Error::msg)?;
    if let Some(settings) = &parameters.ptm_library {
        ensure!(
            context.is_some(),
            "library preview requires --preview-protein and --preview-start"
        );
        let library = if sage_core::ptm_library::is_tsv_path(&settings.path) {
            sage_core::ptm_library::PtmLibrary::from_tsv(&sage_cloudpath::util::read_text(
                &settings.path,
            )?)
            .map_err(anyhow::Error::msg)?
        } else {
            sage_cloudpath::parquet::deserialize_ptm_library(sage_cloudpath::util::read_bytes(
                &settings.path,
            )?)?
        };
        parameters
            .validate_ptm_library(&library)
            .map_err(anyhow::Error::msg)?;
        parameters.loaded_ptm_library = Some(std::sync::Arc::new(library));
    } else {
        parameters
            .validate_ptm_library(&Default::default())
            .map_err(anyhow::Error::msg)?;
        ensure!(
            parameters
                .variable_mods
                .values()
                .flatten()
                .all(|entry| entry.site_mode() == SiteMode::Exhaustive),
            "preview requires a PTM library for non-exhaustive site_mode"
        );
    }
    let protein_position = match position {
        "internal" => Position::Internal,
        "nterm" => Position::Nterm,
        "cterm" => Position::Cterm,
        "full" => Position::Full,
        _ => anyhow::bail!("unknown peptide position `{position}`"),
    };
    let protein_start = context
        .map(|(_, start)| start.checked_sub(1).context("preview-start is one-based"))
        .transpose()?;
    let (sequence_span, protein_start) = if before.is_empty() && after.is_empty() {
        (sequence.into(), protein_start)
    } else {
        // Build a protein-backed span so motif rules see the flanks. Unknown
        // residues before the supplied ones are padded with a non-residue byte.
        let start = protein_start.unwrap_or(before.len() as u32) as usize;
        ensure!(
            before.len() <= start,
            "--preview-before is longer than the residues before --preview-start"
        );
        let protein = format!(
            "{}{before}{sequence}{after}",
            "-".repeat(start - before.len())
        );
        let span = sage_core::sequence::ProteinSequence::from(protein)
            .peptide(start..start + sequence.len())
            .context("preview peptide span")?;
        (span, Some(start as u32))
    };
    let digest = Digest {
        sequence: sequence_span,
        position: protein_position,
        protein: context.map_or("", |(protein, _)| protein).into(),
        protein_start,
        prev_aa: before.bytes().last(),
        next_aa: after.bytes().next(),
        ..Default::default()
    };
    let peptide = Peptide::try_from(digest.clone()).map_err(anyhow::Error::msg)?;
    let rules = rules(&parameters, &peptide);
    let configured_cap = parameters.max_combinations;
    parameters.max_combinations = Some(configured_cap.unwrap_or(usize::MAX).min(limit + 1));
    parameters.generate_decoys = false;
    parameters.peptide_min_mass = 0.0;
    parameters.peptide_max_mass = f32::MAX;
    let indexed = parameters.modify_digests(group_digests(vec![digest]));
    let database = parameters.clone().build_from_peptides(indexed);
    let mut variants = database.peptides.clone();
    'offsets: for peptide in &database.peptides {
        for offset in &database.mass_offsets {
            for site in database.mass_offset_sites(peptide, offset) {
                variants.push(peptide.with_mass_offset(site, &offset.definition));
                if variants.len() > limit {
                    break 'offsets;
                }
            }
        }
    }
    let truncated = variants.len() > limit;
    variants.truncate(limit);
    Ok(json!({
        "sequence": sequence,
        "protein_position": position,
        "rules": rules,
        "limits": {
            "max_variable_mods": parameters.max_variable_mods,
            "max_total_variable_mods": parameters.max_total_variable_mods,
            "max_combinations": configured_cap,
            "preview_limit": limit
        },
        "truncated": truncated,
        "variants": variants.iter().map(|peptide| json!({
            "peptide": peptide.to_string(),
            "monoisotopic_mass": peptide.monoisotopic,
            "modifications": peptide.applied_modifications().map(|applied| json!({
                "site": site_value(applied.site),
                "mass": applied.modification.mass,
                "name": applied.modification.name.as_deref()
            })).collect::<Vec<_>>()
        })).collect::<Vec<_>>()
    }))
}

fn rules(parameters: &Parameters, peptide: &Peptide) -> Vec<Value> {
    let mut output = Vec::new();
    for (specificity, entry) in &parameters.static_mods {
        let mut sites = Vec::new();
        peptide.compatible_sites(*specificity, &mut sites);
        let definition = entry.definition();
        output.push(json!({"kind": "static", "key": specificity.explicit_name(),
            "mass": definition.mass, "name": definition.name.as_deref(),
            "eligible_sites": sites.into_iter().map(site_value).collect::<Vec<_>>() }));
    }
    let resolved = parameters.localization_rules(peptide);
    for (specificity, entries) in &parameters.variable_mods {
        for entry in entries {
            let mut sites = Vec::new();
            let definition = entry.definition();
            if let Some(rule) = resolved.iter().find(|rule| {
                rule.specificity == *specificity && rule.definition.as_ref() == &definition
            }) {
                sites.extend_from_slice(&rule.sites);
            }
            output.push(json!({"kind": "variable", "key": specificity.explicit_name(),
                "mass": definition.mass, "name": definition.name.as_deref(),
                "max_count": entry.max_count(),
                "max_total_count": entry.max_total_count(),
                "search_mode": entry.search_mode(),
                "effective_max_count": if entry.search_mode() == SearchMode::MassOffset { Some(1) } else { entry.max_count() },
                "eligible_sites": sites.into_iter().map(site_value).collect::<Vec<_>>() }));
        }
    }
    output.sort_by_key(|rule| rule.to_string());
    output
}

#[cfg(test)]
#[path = "../tests/unit/modification_preview.rs"]
mod test;
