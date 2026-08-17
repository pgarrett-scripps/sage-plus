//! False discovery rate control using double-competition (picked-peptide &
//! picked-protein) approaches
//!
//! Lin et al., https://pubmed.ncbi.nlm.nih.gov/36166314/
//! Savitski et al., https://pubmed.ncbi.nlm.nih.gov/25987413/

use crate::database::{IndexedDatabase, PeptideIx};
use crate::lfq::{PrecursorId, QuantifiedPeak};
use crate::ml::kde::Estimator;
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

    fn fit_kde<K, B>(scores: &HashMap<K, Self, B>) -> Estimator {
        let (scores, decoys): (Vec<f64>, Vec<bool>) = scores
            .values()
            .map(|score| (score.score() as f64, score.is_decoy()))
            .unzip();
        crate::ml::kde::Builder::default().build(&scores, &decoys)
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

pub fn picked_peptide(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    let mut map: FnvHashMap<String, Competition<PeptideIx>> = FnvHashMap::default();
    for feat in features.iter() {
        let peptide = &db[feat.peptide_idx];
        // Only reverse the peptide sequence if we generated decoys ourselves
        let key = match db.generate_decoys && peptide.decoy {
            true => peptide.reverse().to_string(),
            false => peptide.to_string(),
        };

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

    features.par_iter_mut().for_each(|feat| {
        feat.peptide_q = scores.get(&feat.peptide_idx).copied().unwrap_or(1.0);
    });

    passing
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
    let mut map: FnvHashMap<_, Competition<String>> = FnvHashMap::default();
    for feat in features
        .iter()
        .filter(|x| x.num_protein_groups == 1 && x.protein_groups.is_some())
    {
        let decoy = db[feat.peptide_idx].decoy;
        let entry = map.entry(feat.protein_groups.clone()).or_default();
        match decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = feat.protein_groups.clone();
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = feat.protein_groups.clone();
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features
        .par_iter_mut()
        .filter(|x| x.num_protein_groups == 1 && x.protein_groups.is_some())
        .for_each(|feat| {
            let protein_groups = feat.protein_groups.as_deref().unwrap().to_string();
            feat.protein_group_q = scores[&protein_groups];
        });

    passing
}

pub fn picked_precursor(peaks: &mut FnvHashMap<(PrecursorId, bool), QuantifiedPeak>) -> usize {
    // let mut map: FnvHashMap<PeptideIx, Competition<(PeptideIx, bool)>> = FnvHashMap::default();
    // for (key, (peak, _)) in peaks.iter() {
    //     let entry = map.entry(key.0).or_default();
    //     match key.1 {
    //         true => {
    //             entry.reverse = entry.reverse.max(peak.score as f32);
    //             entry.reverse_ix = Some(*key);
    //         }
    //         false => {
    //             entry.forward = entry.forward.max(peak.score as f32);
    //             entry.foward_ix = Some(*key);
    //         }
    //     }
    // }
    let mut scores = peaks
        .par_iter()
        .map(|(&(ix, decoy), quantified)| Row {
            ix,
            decoy,
            score: quantified.peak.score as f32,
            q: 1.0,
        })
        .collect::<Vec<_>>();

    scores.par_sort_by(|a, b| b.score.total_cmp(&a.score));

    let mut decoy = 1.0;
    let mut target = 0.0;
    for score in scores.iter_mut() {
        match score.decoy {
            true => decoy += 1.0,
            false => target += 1.0,
        };
        score.q = decoy / target;
    }
    // Q-value is the minimum q-value at any given score threshold
    // `q = q[::-1].cummin()[::-1] in python`
    let mut q_min = 1.0f32;
    let mut passing = 0;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.q);
        score.q = q_min;
        if q_min <= 0.05 && !score.decoy {
            passing += 1;
        }
    }

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
mod tests {
    use super::*;
    use crate::enzyme::Digest;
    use crate::modification::{ModificationDefinition, NeutralLossMode};
    use crate::peptide::{AppliedModification, Peptide, Site};
    use std::sync::Arc;

    fn peptide(sequence: &str, decoy: bool) -> Peptide {
        let mut peptide = Peptide::try_from(Digest {
            sequence: sequence.into(),
            ..Default::default()
        })
        .unwrap();
        peptide.decoy = decoy;
        peptide
    }

    #[test]
    fn picked_peptide_assigns_one_to_orphaned_competition_twins() {
        let mut twin_a = peptide("PEPTIDE", false);
        twin_a.modifications = vec![0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        twin_a.applied_modifications = Arc::new(vec![AppliedModification {
            site: Site::Sequence(1),
            modification: Arc::new(ModificationDefinition {
                mass: 10.0,
                name: None,
                neutral_losses: Arc::from([5.0]),
                neutral_loss_mode: NeutralLossMode::Optional,
            }),
        }]);
        let mut twin_b = twin_a.clone();
        twin_b.applied_modifications = Arc::new(vec![AppliedModification {
            site: Site::Sequence(1),
            modification: Arc::new(ModificationDefinition {
                mass: 10.0,
                name: None,
                neutral_losses: Arc::from([6.0]),
                neutral_loss_mode: NeutralLossMode::Optional,
            }),
        }]);

        assert_eq!(twin_a.to_string(), twin_b.to_string());
        let mut twins = vec![twin_a.clone(), twin_b.clone()];
        crate::database::Parameters::reorder_peptides(&mut twins);
        assert_eq!(twins.len(), 2);

        let db = IndexedDatabase {
            peptides: vec![
                twin_a,
                twin_b,
                peptide("AAAAA", false),
                peptide("CCCCC", true),
                peptide("GGGGG", true),
            ],
            generate_decoys: false,
            ..Default::default()
        };
        let mut features = [
            Feature {
                peptide_idx: PeptideIx(0),
                discriminant_score: 10.0,
                ..Default::default()
            },
            Feature {
                peptide_idx: PeptideIx(1),
                discriminant_score: 9.0,
                ..Default::default()
            },
            Feature {
                peptide_idx: PeptideIx(2),
                discriminant_score: 8.0,
                ..Default::default()
            },
            Feature {
                peptide_idx: PeptideIx(3),
                discriminant_score: 7.0,
                ..Default::default()
            },
            Feature {
                peptide_idx: PeptideIx(4),
                discriminant_score: 2.0,
                ..Default::default()
            },
        ];

        picked_peptide(&db, &mut features);

        assert_eq!(features[0].peptide_q, 1.0);
    }
}
