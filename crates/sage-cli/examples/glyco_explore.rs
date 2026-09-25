//! Exploration harness for intact N-glycopeptide search (docs/explore/GLYCO_SEARCH.md).
//!
//! Milestone-1 path, built from existing pieces and no core changes:
//! 1. Gate MS2 spectra on oxonium ions and measure a peptide-free Y-ion ladder rate.
//! 2. Search gated spectra against sequon-bearing peptides only, with an open
//!    precursor window spanning the glycan library and a labile HexNAc mass
//!    offset (fragments carry either nothing or the innermost HexNAc).
//! 3. Explain each candidate's precursor delta with a glycan composition,
//!    choose among explanations with Y ions and sialic-acid oxonium ions,
//!    and estimate PSM FDR with LDA and target-decoy competition.
//!
//! ```text
//! cargo run --release -p sage-cli --example glyco_explore -- \
//!     --raw FILE.raw --fasta FILE.fasta --glycans LIST.glyc [--glycans ...] \
//!     [--high-mannose] --out DIR [--max-spectra N]
//! ```

use rayon::prelude::*;
use sage_core::database::Parameters;
use sage_core::glycan::{
    oxonium_evidence, y_ion_evidence, GlycanComposition, GlycanLibrary, Monosaccharide,
    YIonEvidence,
};
use sage_core::mass::{Tolerance, NEUTRON};
use sage_core::scoring::{Feature, ScoreType, Scorer};
use sage_core::spectrum::{ProcessedSpectrum, SpectrumProcessor};
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

const HEXNAC: f64 = 203.079_373;
/// NH3, the neutral mass of a noncovalent ammonium adduct on the precursor.
const AMMONIA: f64 = 17.026_549;
const TOP_PEAKS: usize = 300;
const REPORT_PSMS: usize = 5;

struct Args {
    raw: String,
    fasta: String,
    glycans: Vec<String>,
    high_mannose: bool,
    out: String,
    max_spectra: Option<usize>,
    bucket_size: usize,
    /// Search every spectrum, not only oxonium-gated ones (timing tests).
    no_gate: bool,
    /// Also explain precursor deltas with one noncovalent ammonium adduct.
    ammonium: bool,
}

fn parse_args() -> Args {
    let mut args = Args {
        raw: String::new(),
        fasta: String::new(),
        glycans: Vec::new(),
        high_mannose: false,
        out: String::new(),
        max_spectra: None,
        bucket_size: 65536,
        no_gate: false,
        ammonium: false,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().expect("missing flag value");
        match flag.as_str() {
            "--raw" => args.raw = value(),
            "--fasta" => args.fasta = value(),
            "--glycans" => args.glycans.push(value()),
            "--high-mannose" => args.high_mannose = true,
            "--no-gate" => args.no_gate = true,
            "--ammonium" => args.ammonium = true,
            "--out" => args.out = value(),
            "--max-spectra" => args.max_spectra = Some(value().parse().unwrap()),
            "--bucket-size" => args.bucket_size = value().parse().unwrap(),
            other => panic!("unknown flag {other}"),
        }
    }
    assert!(!args.raw.is_empty() && !args.fasta.is_empty() && !args.out.is_empty());
    args
}

fn load_library(args: &Args) -> GlycanLibrary {
    let mut compositions = Vec::new();
    for path in &args.glycans {
        let text = std::fs::read_to_string(path).expect("read glycan list");
        for line in text.lines().map(str::trim).filter(|l| !l.is_empty()) {
            match GlycanComposition::parse(line) {
                Ok(c) => compositions.push(c),
                Err(e) => eprintln!("skipping glycan '{line}': {e}"),
            }
        }
    }
    if args.high_mannose {
        for hex in 3..=20 {
            compositions.push(GlycanComposition([2, hex, 0, 0, 0]));
        }
    }
    compositions.sort_by_key(|c| c.0);
    compositions.dedup();
    GlycanLibrary::new(compositions)
}

/// Longest chain of peaks (at one charge state) separated by single
/// monosaccharide masses, above 500 Da. Peptide-free Y-ion ladder evidence.
fn ladder_steps(query: &ProcessedSpectrum, max_charge: u8) -> usize {
    const STEPS: [f32; 3] = [203.079_4, 162.052_8, 146.057_9];
    let mut best = 0;
    for z in 1..=max_charge.max(1) {
        let mut masses: Vec<f32> = (0..query.masses.len())
            .filter_map(|i| {
                if query.has_known_charge(i) {
                    (query.charges[i] == z).then_some(query.masses[i])
                } else {
                    Some(query.masses[i] * z as f32)
                }
            })
            .filter(|m| *m > 500.0)
            .collect();
        masses.sort_by(f32::total_cmp);
        let mut chain = vec![0usize; masses.len()];
        for j in 0..masses.len() {
            for i in 0..j {
                let diff = masses[j] - masses[i];
                let tol = masses[j] * 20e-6 * 2.0;
                if STEPS.iter().any(|s| (diff - s).abs() <= tol) {
                    chain[j] = chain[j].max(chain[i] + 1);
                }
            }
            best = best.max(chain[j]);
        }
    }
    best
}

fn is_yeast_plausible(c: &GlycanComposition) -> bool {
    c.count(Monosaccharide::HexNAc) == 2
        && c.count(Monosaccharide::Fuc) == 0
        && c.count(Monosaccharide::NeuAc) == 0
        && c.count(Monosaccharide::NeuGc) == 0
}

struct Candidate {
    feature: Feature,
    peptide: String,
    bare_mass: f64,
    composition: GlycanComposition,
    isotope: i8,
    adducts: u8,
    error_ppm: f64,
    explanations: usize,
    y: YIonEvidence,
    oxonium: usize,
}

fn main() -> anyhow::Result<()> {
    let args = parse_args();
    std::fs::create_dir_all(&args.out)?;
    let started = Instant::now();
    let library = load_library(&args);
    let (gmin, gmax) = library.mass_range().expect("empty glycan library");
    eprintln!(
        "glycan library: {} compositions, {gmin:.2}-{gmax:.2} Da",
        library.len()
    );

    // Configuration through the normal input path, so the HexNAc offset and
    // its motif are parsed exactly as a user config would be.
    let config = serde_json::json!({
        "database": {
            "bucket_size": args.bucket_size,
            "fasta": args.fasta,
            "decoy_tag": "rev_",
            "generate_decoys": true,
            "enzyme": {"missed_cleavages": 2, "min_len": 5, "max_len": 50},
            "peptide_min_mass": 500.0,
            "peptide_max_mass": 5000.0,
            "static_mods": {"C": 57.021464},
            "variable_mods": {
                "HexNAc": {
                    "mass": HEXNAC,
                    "sites": ["motif:N*-{P}-[ST]"],
                    "max_count": 1,
                    "search_mode": "mass_offset",
                    "neutral_losses": [HEXNAC]
                }
            },
            "max_variable_mods": 1
        },
        "precursor_tol": {"ppm": [-20.0, 20.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "mzml_paths": [args.raw]
    });
    let config_path = format!("{}/config.json", args.out);
    std::fs::write(&config_path, serde_json::to_string_pretty(&config)?)?;
    let search = sage_cli::input::Input::load(&config_path)?.build()?;
    let parameters: Parameters = search.database;

    let fasta = sage_cloudpath::util::read_fasta(
        &sage_cloudpath::to_url(&args.fasta)?,
        &parameters.decoy_tag,
        parameters.generate_decoys,
    )?;
    let offsets = parameters.mass_offset_modifications();
    let all = parameters.digest(&fasta);
    let total_peptides = all.len();
    let mut sites = Vec::new();
    let sequon: Vec<_> = all
        .into_iter()
        .filter(|peptide| {
            sites.clear();
            for offset in &offsets {
                for specificity in &offset.specificities {
                    peptide.compatible_sites(*specificity, &mut sites);
                }
            }
            !sites.is_empty()
        })
        .collect();
    eprintln!(
        "sequon peptides: {} of {} ({:.1}%), {:.0}s",
        sequon.len(),
        total_peptides,
        100.0 * sequon.len() as f64 / total_peptides as f64,
        started.elapsed().as_secs_f64()
    );
    let n_sequon = sequon.len();
    let db = parameters.build_from_peptides(sequon);
    eprintln!("index built, {:.0}s", started.elapsed().as_secs_f64());

    let raw = sage_cloudpath::util::read_thermoraw(&sage_cloudpath::to_url(&args.raw)?, 0, None)
        .map_err(|e| anyhow::anyhow!("{e:?}"))?;
    let processor = SpectrumProcessor::new(TOP_PEAKS, true, 0.0);
    let mut spectra: Vec<ProcessedSpectrum> = raw
        .into_par_iter()
        .filter(|s| s.ms_level == 2 && !s.precursors.is_empty())
        .map(|s| processor.process(s))
        .collect();
    if let Some(n) = args.max_spectra {
        spectra.truncate(n);
    }
    eprintln!(
        "{} MS2 spectra read, {:.0}s",
        spectra.len(),
        started.elapsed().as_secs_f64()
    );

    let ppm20 = Tolerance::Ppm(-20.0, 20.0);
    // Offset-0 hypothesis covers bare peptide + [gmin - HexNAc, gmax]; the HexNAc
    // offset hypothesis covers bare peptide + HexNAc + the same range.
    let precursor_tol = Tolerance::Da(-(gmax as f32) - 1.5, -((gmin - HEXNAC) as f32) + 0.1);
    let scorer = Scorer {
        db: &db,
        precursor_tol,
        fragment_tol: ppm20,
        min_matched_peaks: 4,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 2,
        max_precursor_charge: 6,
        override_precursor_charge: false,
        max_fragment_charge: Some(3),
        chimera: false,
        report_psms: REPORT_PSMS,
        wide_window: false,
        annotate_matches: false,
        mass_shift_ppm: 20.0,
        score_type: ScoreType::SageHyperScore,
        mass_recalibration: None,
    };

    // Per-spectrum: gate, ladder, and (for gated spectra) the best explained candidate.
    let with_features = std::sync::atomic::AtomicUsize::new(0);
    let per_spectrum: Vec<(bool, usize, Option<Candidate>)> = spectra
        .par_iter()
        .map(|query| {
            let oxonium = oxonium_evidence(query, ppm20);
            let gated = oxonium.is_glyco(1) || args.no_gate;
            let precursor_charge = query.precursors.first().and_then(|p| p.charge).unwrap_or(3);
            let ladder = ladder_steps(query, precursor_charge.min(4));
            if !gated {
                return (false, ladder, None);
            }
            let sialic = oxonium.matched[6] || oxonium.matched[7];
            let neugc = oxonium.matched[8];
            let mut best = None;
            let features = scorer.score(query);
            if !features.is_empty() {
                with_features.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            for feature in features {
                let base = &db.peptides[feature.peptide_idx.0 as usize];
                let bare_mass = base.monoisotopic as f64;
                let delta = feature.expmass as f64 - bare_mass;
                let tol_da = feature.expmass as f64 * 20e-6;
                // Each explanation: (library index, isotope error, ammonium adducts).
                let explanations: Vec<(usize, i8, u8)> = (0..=u8::from(args.ammonium))
                    .flat_map(|adducts| {
                        library
                            .explain(delta - adducts as f64 * AMMONIA, tol_da, 0..=1)
                            .into_iter()
                            .map(move |(index, isotope)| (index, isotope, adducts))
                    })
                    .collect();
                if explanations.is_empty() {
                    continue;
                }
                // Rank: sialic-acid oxonium consistency, Y-ion matches, then the
                // simpler explanation (no adduct, isotope 0), then mass error.
                let scored = explanations
                    .iter()
                    .map(|&(index, isotope, adducts)| {
                        let composition = *library.get(index).unwrap();
                        let y = y_ion_evidence(
                            query,
                            bare_mass as f32,
                            &composition,
                            feature.charge,
                            ppm20,
                        );
                        let consistent =
                            u8::from((composition.count(Monosaccharide::NeuAc) > 0) == sialic)
                                + u8::from((composition.count(Monosaccharide::NeuGc) > 0) == neugc);
                        let error = delta
                            - isotope as f64 * NEUTRON as f64
                            - adducts as f64 * AMMONIA
                            - composition.mass();
                        (composition, isotope, error, y, consistent, adducts)
                    })
                    .max_by(|a, b| {
                        (a.4, a.3.matched, -(a.5 as i32), -(a.1 as i32))
                            .cmp(&(b.4, b.3.matched, -(b.5 as i32), -(b.1 as i32)))
                            .then(b.2.abs().total_cmp(&a.2.abs()))
                    })
                    .unwrap();
                let (composition, isotope, error, y, _, adducts) = scored;
                best = Some(Candidate {
                    peptide: db.resolve_peptide(&feature).to_string(),
                    feature,
                    bare_mass,
                    composition,
                    isotope,
                    adducts,
                    error_ppm: error / (bare_mass + delta) * 1e6,
                    explanations: explanations.len(),
                    y,
                    oxonium: oxonium.count(),
                });
                // Features arrive best-first; keep the top explained one.
                break;
            }
            (true, ladder, best)
        })
        .collect();
    eprintln!("scored, {:.0}s", started.elapsed().as_secs_f64());

    let n = per_spectrum.len();
    let gated = per_spectrum.iter().filter(|s| s.0).count();
    let ladder = |gate: bool, steps: usize| {
        per_spectrum
            .iter()
            .filter(|s| s.0 == gate && s.1 >= steps)
            .count()
    };
    let ladder_rows: Vec<(usize, usize, usize)> = (2..=5)
        .map(|steps| (steps, ladder(true, steps), ladder(false, steps)))
        .collect();
    let mut candidates: Vec<Candidate> = per_spectrum.into_iter().filter_map(|s| s.2).collect();

    // FDR: LDA over the standard feature set, then target-decoy q-values.
    let mut features: Vec<Feature> = candidates.iter().map(|c| c.feature.clone()).collect();
    let lda = sage_core::ml::linear_discriminant::score_psms(&mut features, precursor_tol);
    if lda.is_none() {
        eprintln!("LDA failed; ranking by hyperscore");
        for f in &mut features {
            f.discriminant_score = f.hyperscore as f32;
        }
    }
    for (c, f) in candidates.iter_mut().zip(features) {
        c.feature = f;
    }
    candidates.sort_by(|a, b| {
        b.feature
            .discriminant_score
            .total_cmp(&a.feature.discriminant_score)
    });
    let mut features: Vec<Feature> = candidates.iter().map(|c| c.feature.clone()).collect();
    sage_core::ml::qvalue::spectrum_q_value(&mut features);
    for (c, f) in candidates.iter_mut().zip(features) {
        c.feature.spectrum_q = f.spectrum_q;
    }

    let mut tsv = std::io::BufWriter::new(std::fs::File::create(format!(
        "{}/glyco_psms.tsv",
        args.out
    ))?);
    writeln!(
        tsv,
        "spec_id\tlabel\tq\tdiscriminant\thyperscore\tcharge\tpeptide\tbare_mass\texpmass\tglycan\tglycan_mass\tisotope\tadducts\terror_ppm\texplanations\ty_matched\ty_generated\ty_anchored\toxonium"
    )?;
    for c in &candidates {
        let f = &c.feature;
        writeln!(
            tsv,
            "{}\t{}\t{:.5}\t{:.4}\t{:.3}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{:.4}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}",
            f.spec_id,
            f.label,
            f.spectrum_q,
            f.discriminant_score,
            f.hyperscore,
            f.charge,
            c.peptide,
            c.bare_mass,
            f.expmass,
            c.composition,
            c.composition.mass(),
            c.isotope,
            c.adducts,
            c.error_ppm,
            c.explanations,
            c.y.matched,
            c.y.generated,
            c.y.anchored,
            c.oxonium
        )?;
    }
    tsv.flush()?;

    let passing: Vec<&Candidate> = candidates
        .iter()
        .filter(|c| c.feature.spectrum_q <= 0.01 && c.feature.label == 1)
        .collect();
    let decoys_passing = candidates
        .iter()
        .filter(|c| c.feature.spectrum_q <= 0.01 && c.feature.label == -1)
        .count();
    let mut by_composition: HashMap<String, usize> = HashMap::new();
    let mut unique: HashMap<(String, String), ()> = HashMap::new();
    let mut peptides: HashMap<String, ()> = HashMap::new();
    for c in &passing {
        *by_composition.entry(c.composition.to_string()).or_default() += 1;
        unique.insert((c.peptide.clone(), c.composition.to_string()), ());
        peptides.insert(c.peptide.clone(), ());
    }
    let mut top: Vec<_> = by_composition.into_iter().collect();
    top.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let anchored = passing.iter().filter(|c| c.y.anchored).count();
    let y_any = passing.iter().filter(|c| c.y.matched >= 2).count();
    let ambiguous = passing.iter().filter(|c| c.explanations > 1).count();
    let iso1 = passing.iter().filter(|c| c.isotope != 0).count();
    let adducted = passing.iter().filter(|c| c.adducts != 0).count();
    let implausible_yeast = passing
        .iter()
        .filter(|c| !is_yeast_plausible(&c.composition))
        .count();
    let pct = |a: usize, b: usize| 100.0 * a as f64 / b.max(1) as f64;

    let mut summary = String::new();
    use std::fmt::Write as _;
    writeln!(summary, "raw\t{}", args.raw)?;
    writeln!(
        summary,
        "glycan_library\t{} compositions ({gmin:.1}-{gmax:.1} Da)",
        library.len()
    )?;
    writeln!(summary, "sequon_peptides\t{n_sequon} of {total_peptides}")?;
    writeln!(summary, "ms2_spectra\t{n}")?;
    writeln!(summary, "oxonium_gated\t{gated} ({:.1}%)", pct(gated, n))?;
    for (steps, on, off) in &ladder_rows {
        writeln!(
            summary,
            "ladder_ge{steps}_steps\tgated {on} ({:.1}%)\tungated {off} ({:.1}%)",
            pct(*on, gated),
            pct(*off, n - gated)
        )?;
    }
    writeln!(
        summary,
        "searched_with_any_candidate\t{}",
        with_features.load(std::sync::atomic::Ordering::Relaxed)
    )?;
    writeln!(summary, "explained_candidates\t{}", candidates.len())?;
    writeln!(summary, "lda_fit\t{}", lda.is_some())?;
    writeln!(
        summary,
        "glyco_psms_1pct\t{} (decoys {decoys_passing})",
        passing.len()
    )?;
    writeln!(summary, "unique_peptide_glycan_1pct\t{}", unique.len())?;
    writeln!(summary, "unique_peptides_1pct\t{}", peptides.len())?;
    writeln!(
        summary,
        "y_anchored_1pct\t{anchored} ({:.1}%)",
        pct(anchored, passing.len())
    )?;
    writeln!(
        summary,
        "y_ge2_1pct\t{y_any} ({:.1}%)",
        pct(y_any, passing.len())
    )?;
    writeln!(
        summary,
        "ambiguous_composition_1pct\t{ambiguous} ({:.1}%)",
        pct(ambiguous, passing.len())
    )?;
    writeln!(
        summary,
        "isotope_1_assigned_1pct\t{iso1} ({:.1}%)",
        pct(iso1, passing.len())
    )?;
    writeln!(
        summary,
        "ammonium_adduct_1pct\t{adducted} ({:.1}%)",
        pct(adducted, passing.len())
    )?;
    writeln!(
        summary,
        "non_high_mannose_1pct\t{implausible_yeast} ({:.1}%)",
        pct(implausible_yeast, passing.len())
    )?;
    writeln!(summary, "elapsed_s\t{:.0}", started.elapsed().as_secs_f64())?;
    writeln!(summary, "top_compositions_1pct")?;
    for (name, count) in top.iter().take(15) {
        writeln!(summary, "  {name}\t{count}")?;
    }
    print!("{summary}");
    std::fs::write(format!("{}/summary.tsv", args.out), summary)?;
    Ok(())
}
