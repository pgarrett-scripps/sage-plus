use super::*;

#[test]
fn index_budget_follows_max_memory() {
    const GIB: u64 = 1024 * 1024 * 1024;
    // A quarter of 16 GiB indexes spectra; 8 GiB without a limit.
    assert_eq!(
        PrefilterBudgets::from_max_memory(Some(16.0)).index_bytes,
        4 * GIB
    );
    assert_eq!(PrefilterBudgets::from_max_memory(None).index_bytes, 8 * GIB);
    assert_eq!(
        PrefilterBudgets::from_max_memory(Some(0.0)).index_bytes,
        8 * GIB
    );
}
#[test]
fn isoleucine_leucine_twins_share_a_prefilter_key() {
    // APEPIDEK and APEPJDEK are twins; DFNEEGIK is a twin of the decoy
    // of DLGEENFK. Each digest is unique as written.
    let fasta = Fasta::parse(
        ">a\nGGRAPEPIDEKR\n>b\nWWKAPEPJDEKR\n>c\nDLGEENFKR\n>d\nMSSRDFNEEGIKR\n".into(),
        "rev_",
        true,
    )
    .unwrap();
    let enzyme = sage_core::database::Builder::default()
        .make_parameters()
        .expect("valid test configuration")
        .enzyme_parameters();
    for merge in [false, true] {
        let shared = SharedSequences::scan(&fasta, &enzyme, None, true, merge);
        let key = |protein: usize, sequence: &str| {
            let digest = fasta
                .digest_protein(protein, &enzyme, None)
                .into_iter()
                .find(|digest| digest.sequence.as_str() == sequence)
                .unwrap();
            shared.shared_slot(&digest).map(|(key, _)| key)
        };
        let (i, j) = (key(0, "APEPIDEK"), key(1, "APEPJDEK"));
        let (target, collision) = (key(2, "DLGEENFK"), key(3, "DFNEEGIK"));
        match merge {
            true => {
                assert!(i.is_some() && i == j);
                assert!(target.is_some() && target == collision);
            }
            false => assert_eq!([i, j, target, collision], [None; 4]),
        }
    }
}

/// The property the early test rests on: a unit keeps a peptide under the
/// exact filter only if some copy of its key, expanded on its own, matches
/// the spectrum index. Checked over every shared key of a database built to
/// share digests in each way, at the loosest and the default match count.
#[test]
fn a_unit_none_of_whose_copies_matches_keeps_nothing() {
    let workspace = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let p = "LQSRPAAPPAPGPGQLTLR";
    let text = [
        format!(">sp|A|INTERNAL\nMSGGWGGK{p}AAGEDK\n"),
        format!(">sp|B|NTERM\n{p}WSPEDK\n"),
        format!(">sp|C|CLIPPED\nM{p}GGHEK\n"),
        format!(">sp|D|CTERM\nGGWGGK{p}\n"),
        ">sp|E|ITWIN\nAAKLQSRPAAPPAPGPGQITLRAAK\n".to_string(),
        ">sp|P02768|ALBU_HUMAN\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK\n".to_string(),
        ">sp|REV|REVERSED_DIGEST\nMSSRDFNEEGLKAGSEYR\n".to_string(),
        ">sp|H|PALINDROMES\nGAPAKGAPAKWPEDGKSPAPSR\n".to_string(),
        std::fs::read_to_string(workspace.join("tests/Q99536.fasta")).unwrap(),
    ]
    .concat();
    let fasta = Fasta::parse(text.clone(), "rev_", true).unwrap();
    let mut live_units = 0;
    for min_matched in [1u16, 4] {
        let mut config: serde_json::Value =
            serde_json::from_slice(&std::fs::read(workspace.join("tests/config.json")).unwrap())
                .unwrap();
        config["database"]["fasta"] = workspace
            .join("tests/Q99536.fasta")
            .display()
            .to_string()
            .into();
        config["database"]["prefilter"] = true.into();
        config["database"]["prefilter_min_matched_peaks"] = min_matched.into();
        config["database"]["variable_mods"] = serde_json::json!({
            "Acetyl": {"mass": 42.010565, "sites": ["protein_n_term"], "max_count": 1},
            "Oxidation": {"mass": 15.994915, "sites": ["M"], "max_count": 1}
        });
        config["mzml_paths"] = serde_json::json!([workspace
            .join("tests/LQSRPAAPPAPGPGQLTLR.mzML")
            .display()
            .to_string()]);
        let input: crate::input::Input = serde_json::from_value(config).unwrap();
        let runner = Runner::new(input.build().unwrap(), 1).unwrap();
        let db_params = runner.database_parameters.clone();
        let enzyme = db_params.enzyme_parameters();
        let mut builder = runner.spectrum_index_builder(&db_params);
        let urls = runner.parameters.mzml_paths.clone();
        let spectra = runner
            .read_processed_spectra_with_ms1(&urls, 0, 1, false, false)
            .unwrap();
        builder.add(&spectra.1);
        let index = Runner::finish_spectrum_index(builder);
        let shared = SharedSequences::scan(
            &fasta,
            &enzyme,
            None,
            db_params.generate_decoys,
            db_params.merge_isoleucine_leucine,
        );
        let mut deferred = Vec::new();
        for protein in 0..fasta.targets.len() {
            for digest in fasta.digest_protein(protein, &enzyme, None) {
                if let Some((key, _)) = shared.shared_slot(&digest) {
                    deferred.push((key, digest));
                }
            }
        }
        // Every copy of each key, as streamed.
        let mut copies: HashMap<u64, Vec<Digest>> = HashMap::new();
        for (key, digest) in &deferred {
            copies.entry(*key).or_default().push(digest.clone());
        }
        let units = shared_units(deferred);
        let unit_key = |unit: &SharedUnit| {
            shared
                .shared_slot(&unit.groups[0].reference)
                .map(|(key, _)| key)
                .unwrap()
        };
        assert!(
            units.len() >= 3,
            "the database shares digests ({} units)",
            units.len()
        );
        db_params.with_digest_expander(|expander| {
            let mut scratch = index.scratch();
            for mut unit in units {
                // The early test of every copy, as the stream runs it.
                let early = copies
                    .remove(&unit_key(&unit))
                    .unwrap_or_default()
                    .into_iter()
                    .any(|digest| copy_could_match(expander, &index, &mut scratch, digest));
                // The exact filter of the unit, as the stream runs it.
                let targets = unit
                    .groups
                    .iter()
                    .filter(|group| !group.reference.decoy)
                    .map(|group| {
                        expander
                            .parameters()
                            .decoy_collision_key(&group.reference.sequence)
                    })
                    .collect::<HashSet<_>>();
                let mut peptides = std::mem::take(&mut unit.groups)
                    .into_iter()
                    .flat_map(|group| expander.expand(group, Some(&targets)))
                    .collect::<Vec<_>>();
                filter_unit(
                    expander,
                    &index,
                    &mut scratch,
                    &mut peptides,
                    &mut unit.keep,
                    true,
                    true,
                );
                if !peptides.is_empty() {
                    live_units += 1;
                    assert!(
                        early,
                        "a unit kept peptides that no copy's early test matched"
                    );
                }
            }
        });
    }
    assert!(
        live_units > 0,
        "some unit keeps peptides, so the property is tested"
    );
}
