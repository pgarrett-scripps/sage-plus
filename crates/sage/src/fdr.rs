//! False discovery rate control using double-competition (picked-peptide &
//! picked-protein) approaches
//!
//! Lin et al., <https://pubmed.ncbi.nlm.nih.gov/36166314/>
//! Savitski et al., <https://pubmed.ncbi.nlm.nih.gov/25987413/>

use crate::database::{IndexedDatabase, PeptideIx};
use crate::lfq::{PrecursorId, QuantifiedPeak};
use crate::ml::kde::Estimator;
use crate::peptide::Peptide;
use crate::scoring::Feature;
use fnv::FnvHashMap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::BuildHasher;

#[derive(Copy, Clone, Debug)]
pub struct Competition<Ix> {
    pub forward: f32,
    pub foward_ix: Option<Ix>,
    pub reverse: f32,
    pub reverse_ix: Option<Ix>,
}

struct Row<Ix> {
    ix: Ix,
    decoy: bool,
    score: f32,
    q: f32,
}

impl<Ix: Default + Send> Default for Competition<Ix> {
    fn default() -> Self {
        Self {
            forward: f32::MIN,
            reverse: f32::MIN,
            foward_ix: None,
            reverse_ix: None,
        }
    }
}

impl<Ix: Default + Send> Competition<Ix> {
    fn score(&self) -> f32 {
        self.forward.max(self.reverse)
    }

    fn is_decoy(&self) -> bool {
        self.reverse >= self.forward
    }

    fn fit_kde<K, B>(scores: &HashMap<K, Self, B>) -> Option<Estimator> {
        let (scores, decoys): (Vec<f64>, Vec<bool>) = scores
            .values()
            .map(|score| (score.score() as f64, score.is_decoy()))
            .unzip();
        for class in [false, true] {
            let sample = scores
                .iter()
                .zip(&decoys)
                .filter_map(|(score, decoy)| (*decoy == class).then_some(*score))
                .collect::<Vec<_>>();
            if sample.len() < 2
                || sample.iter().any(|score| !score.is_finite())
                || !sample.iter().any(|score| *score != sample[0])
            {
                return None;
            }
        }
        let estimator = crate::ml::kde::Builder::default().build(&scores, &decoys);
        scores
            .iter()
            .all(|score| estimator.posterior_error(*score).is_finite())
            .then_some(estimator)
    }

    fn assign_q_value<K, B>(
        scores: HashMap<K, Self, B>,
        threshold: f32,
    ) -> (HashMap<Ix, f32, B>, usize)
    where
        K: Eq + std::hash::Hash + Send,
        Ix: Eq + std::hash::Hash,
        B: BuildHasher + Default + Send,
    {
        let estimator = Self::fit_kde(&scores);
        let mut scores = scores
            .into_par_iter()
            .flat_map(|(_, comp)| {
                [
                    (comp.foward_ix, false, comp.forward),
                    (comp.reverse_ix, true, comp.reverse),
                ]
            })
            .filter_map(|(ix, decoy, score)| {
                ix.map(|ix| Row {
                    ix,
                    decoy,
                    score,
                    q: 1.0,
                })
            })
            .collect::<Vec<Row<Ix>>>();

        scores.par_sort_by(|a, b| b.score.total_cmp(&a.score));

        if estimator.is_none() {
            log::warn!("peptide/protein confidence model is underdetermined, using target-decoy counts with a +1 correction");
            let passing = assign_count_q_values(&mut scores, threshold);
            return (
                scores
                    .into_iter()
                    .map(|score| (score.ix, score.q))
                    .collect(),
                passing,
            );
        }
        let estimator = estimator.unwrap();

        let mut decoy = 1.0;
        let mut target = 0.0;
        for score in scores.iter_mut() {
            let pep = estimator.posterior_error(score.score as f64) as f32;

            // Cumulative sum of PEP ~ # of decoys
            decoy += pep;
            if !score.decoy {
                target += 1.0;
            }
            score.q = decoy / target;
        }
        // Q-value is the minimum q-value at any given score threshold
        // `q = q[::-1].cummin()[::-1] in python`
        let mut q_min = 1.0f32;
        let mut passing = 0;
        for score in scores.iter_mut().rev() {
            q_min = q_min.min(score.q);
            score.q = q_min;
            if q_min <= threshold && !score.decoy {
                passing += 1;
            }
        }

        (
            scores
                .into_par_iter()
                .map(|score| (score.ix, score.q))
                .collect(),
            passing,
        )
    }
}

/// Assign count-based q-values at complete score thresholds, keeping ties together.
fn assign_count_q_values<Ix>(scores: &mut [Row<Ix>], threshold: f32) -> usize {
    scores.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut decoys = 1.0_f64;
    let mut targets = 0.0_f64;
    let mut start = 0;
    while start < scores.len() {
        let mut end = start + 1;
        while end < scores.len() && scores[end].score == scores[start].score {
            end += 1;
        }
        for row in &scores[start..end] {
            if row.decoy {
                decoys += 1.0;
            } else {
                targets += 1.0;
            }
        }
        let q = (decoys / targets).min(1.0) as f32;
        for row in &mut scores[start..end] {
            row.q = q;
        }
        start = end;
    }
    let mut minimum = 1.0_f32;
    let mut passing = 0;
    for row in scores.iter_mut().rev() {
        minimum = minimum.min(row.q);
        row.q = minimum;
        passing += usize::from(!row.decoy && row.q <= threshold);
    }
    passing
}

pub fn picked_peptide(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    let mut map: FnvHashMap<String, Competition<PeptideIx>> = FnvHashMap::default();
    for feat in features.iter() {
        let peptide = &db[feat.peptide_idx];
        let key = peptide_competition_key(db, feat.peptide_idx);

        let entry = map.entry(key).or_default();
        match peptide.decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some(feat.peptide_idx);
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some(feat.peptide_idx);
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);
    let label_scores = scores
        .iter()
        .filter_map(|(peptide_idx, q)| {
            let peptide = &db[*peptide_idx];
            peptide.label_channel.as_ref().map(|_| {
                (
                    (peptide_competition_key(db, *peptide_idx), peptide.decoy),
                    *q,
                )
            })
        })
        .collect::<FnvHashMap<_, _>>();

    features.par_iter_mut().for_each(|feat| {
        let peptide = &db[feat.peptide_idx];
        feat.peptide_q = peptide
            .label_channel
            .as_ref()
            .and_then(|_| {
                label_scores
                    .get(&(peptide_competition_key(db, feat.peptide_idx), peptide.decoy))
                    .copied()
            })
            .or_else(|| scores.get(&feat.peptide_idx).copied())
            .unwrap_or(1.0);
    });

    passing
}

fn peptide_competition_key(db: &IndexedDatabase, peptide_idx: PeptideIx) -> String {
    let peptide = &db[peptide_idx];
    if peptide.decoy {
        db.decoy_pairing
            .get(peptide_idx.0 as usize)
            .and_then(|paired| db.peptides.get(paired.0 as usize))
            .map(Peptide::label_group)
            .unwrap_or_else(|| {
                if db.generate_decoys {
                    peptide.reverse().label_group()
                } else {
                    peptide.label_group()
                }
            })
    } else {
        peptide.label_group()
    }
}

pub fn picked_protein(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    // Critical: All non-proteotypic, non-unique, or shared peptides are discarded
    // else the assumptions of picked protein FDR are invalid. Shared peptides are
    // still reported, albeit with protein FDR = 1.0
    let mut map: FnvHashMap<_, Competition<String>> = FnvHashMap::default();
    for feat in features
        .iter()
        .filter(|x| db[x.peptide_idx].proteins.len() == 1)
    {
        let decoy = db[feat.peptide_idx].decoy;
        let entry = map.entry(&db[feat.peptide_idx].proteins).or_default();
        let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        match decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some(proteins);
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some(proteins);
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features
        .par_iter_mut()
        .filter(|x| db[x.peptide_idx].proteins.len() == 1)
        .for_each(|feat| {
            let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
            feat.protein_q = scores[&proteins];
        });

    passing
}

pub fn picked_protein_group(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    // Critical: All non-proteotypic, non-unique, or shared peptides are discarded
    // else the assumptions of picked group FDR are invalid. Shared peptides are
    // still reported, albeit with protein group FDR = 1.0
    let target_groups = target_protein_groups(db, features);
    let competition_group = |feat: &Feature| match db[feat.peptide_idx].decoy {
        true => decoy_competition_group(db, feat, &target_groups),
        false => feat.protein_groups.clone(),
    };
    // Several decoy accessions can compete under one target group, so q-values
    // are keyed by the competition group and side rather than by each
    // feature's reported group string.
    let mut map: FnvHashMap<_, Competition<(String, bool)>> = FnvHashMap::default();
    for feat in features
        .iter()
        .filter(|x| x.num_protein_groups == 1 && x.protein_groups.is_some())
    {
        let decoy = db[feat.peptide_idx].decoy;
        let key = competition_group(feat).unwrap_or_default();
        let entry = map.entry(key.clone()).or_default();
        match decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some((key, true));
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some((key, false));
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    let groups = features
        .iter()
        .map(|feat| {
            (feat.num_protein_groups == 1 && feat.protein_groups.is_some())
                .then(|| competition_group(feat).unwrap_or_default())
        })
        .collect::<Vec<_>>();
    features
        .par_iter_mut()
        .zip(groups)
        .for_each(|(feat, group)| {
            if let Some(group) = group {
                let decoy = db[feat.peptide_idx].decoy;
                feat.protein_group_q = scores[&(group, decoy)];
            }
        });

    passing
}

/// Target protein accession to the single protein group reported for it.
/// Only groups that list the accession as a member are considered. When an
/// accession belongs to several reported groups, the group with the best
/// target score wins (ties go to the smallest group string), so a fallback
/// group from a low-confidence peptide does not take the pairing away from
/// the group that carries the protein's real evidence.
fn target_protein_groups<'a>(
    db: &'a IndexedDatabase,
    features: &'a [Feature],
) -> FnvHashMap<&'a str, &'a str> {
    let mut best: FnvHashMap<&str, f32> = FnvHashMap::default();
    for feat in features.iter().filter(|x| x.num_protein_groups == 1) {
        let Some(group) = feat.protein_groups.as_deref() else {
            continue;
        };
        if db[feat.peptide_idx].decoy {
            continue;
        }
        let score = best.entry(group).or_insert(f32::NEG_INFINITY);
        *score = score.max(feat.discriminant_score);
    }

    let mut groups: FnvHashMap<&str, (&str, f32)> = FnvHashMap::default();
    for (&group, &score) in &best {
        for protein in group.split('/') {
            groups
                .entry(protein)
                .and_modify(|current| {
                    let (current_group, current_score) = *current;
                    if score > current_score || (score == current_score && group < current_group) {
                        *current = (group, score);
                    }
                })
                .or_insert((group, score));
        }
    }
    groups
        .into_iter()
        .map(|(protein, (group, _))| (protein, group))
        .collect()
}

/// Decoy features are not grouped, so a decoy competes under the group of the
/// target protein it was reversed from. Decoys without a reported target group
/// compete under their own accession and remain unpaired.
fn decoy_competition_group(
    db: &IndexedDatabase,
    feat: &Feature,
    target_groups: &FnvHashMap<&str, &str>,
) -> Option<String> {
    let peptide = &db[feat.peptide_idx];
    let target = match peptide.proteins.as_slice() {
        [protein] if db.generate_decoys => Some(protein.to_string()),
        [protein] => protein
            .find(db.decoy_tag.as_str())
            .map(|at| format!("{}{}", &protein[..at], &protein[at + db.decoy_tag.len()..])),
        _ => None,
    };
    target
        .and_then(|protein| target_groups.get(protein.as_str()))
        .map(|group| group.to_string())
        .or_else(|| feat.protein_groups.clone())
}

pub fn picked_precursor(peaks: &mut FnvHashMap<(PrecursorId, bool), QuantifiedPeak>) -> usize {
    let mut scores = peaks
        .par_iter()
        .map(|(&(ix, decoy), quantified)| Row {
            ix,
            decoy,
            score: quantified.peak.score as f32,
            q: 1.0,
        })
        .collect::<Vec<_>>();

    let passing = assign_count_q_values(&mut scores, 0.05);

    let scores = scores
        .into_par_iter()
        .map(|score| ((score.ix, score.decoy), score.q))
        .collect::<FnvHashMap<_, _>>();

    peaks.par_iter_mut().for_each(|(ix, quantified)| {
        quantified.peak.q_value = scores[ix];
    });
    passing
}

#[cfg(test)]
#[path = "../tests/unit/fdr.rs"]
mod tests;
