//! The fast chimeric search: repeated scoring with matched peaks removed.
//!
//! Moved out of `scoring.rs` unchanged; every item keeps its `scoring::` path.

use crate::database::MassOffsetAssignment;
use crate::spectrum::ProcessedSpectrum;
use std::borrow::Cow;

use super::*;

impl<'db> Scorer<'db> {
    /// Return multiple PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed, etc
    pub fn score_chimera_fast(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        // As in `score_standard`: no selected precursor means no searchable
        // precursor mass, so report no match.
        let Some(precursor) = query.precursors.first() else {
            return Vec::new();
        };

        let mut query = query.clone();
        let hits = self.initial_hits(&query, precursor);

        let mut candidates: Vec<Feature> = Vec::with_capacity(self.report_psms);

        let mut prev = 0;
        while candidates.len() < self.report_psms {
            self.build_features(&query, precursor, &hits, 1, &mut candidates);
            if candidates.len() > prev {
                if let Some(feat) = candidates.get_mut(prev) {
                    self.remove_matched_peaks(&mut query, feat);
                    feat.rank = prev as u32 + 1;
                }
                prev = candidates.len()
            } else {
                break;
            }
        }
        candidates
    }

    /// Calculate full hyperscore for an indexed candidate.
    #[cfg(test)]
    pub(super) fn score_candidate(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
        collect_fragments: bool,
    ) -> (Score, Option<Fragments>, Coverage) {
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, pre_score.precursor_charge);
        let fragment_index = FragmentMatchIndex::new(query, max_fragment_charge);
        let peptide = &self.db[pre_score.peptide];
        self.score_peptide(
            query,
            peptide,
            pre_score,
            collect_fragments,
            &fragment_index,
        )
    }

    /// Fully score a preliminary hypothesis. Every compatible placement of
    /// an offset is scored and competes as its own candidate, exactly as the
    /// placements of an expanded database would; site confidence is assessed
    /// later by PTM localization.
    pub(super) fn score_hypothesis(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
        collect_fragments: bool,
        fragment_index: &FragmentMatchIndex,
        candidates: &mut Vec<Candidate<'db>>,
    ) {
        let base = &self.db[pre_score.peptide];
        if pre_score.offset == 0 {
            let (score, fragments, coverage) =
                self.score_peptide(query, base, pre_score, collect_fragments, fragment_index);
            candidates.push(Candidate {
                score,
                fragments,
                coverage,
                peptide: Cow::Borrowed(base),
                mass_offset: None,
            });
            return;
        }

        let offset_index = pre_score.offset as usize - 1;
        let offset = &self.db.mass_offsets[offset_index];
        for site in self.db.mass_offset_sites(base, offset) {
            let peptide = base.with_mass_offset(site, &offset.definition);
            let (score, fragments, coverage) = self.score_peptide(
                query,
                &peptide,
                pre_score,
                collect_fragments,
                fragment_index,
            );
            candidates.push(Candidate {
                score,
                fragments,
                coverage,
                peptide: Cow::Owned(peptide),
                mass_offset: Some(MassOffsetAssignment {
                    offset: offset_index as u16,
                    site,
                }),
            });
        }
    }
}
