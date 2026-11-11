//! The standard, non-chimeric search, including open-search mass offsets.
//!
//! Moved out of `scoring.rs` unchanged; every item keeps its `scoring::` path.

use crate::database::PeptideIx;
use crate::heap::bounded_min_heapify;
use crate::spectrum::ProcessedSpectrum;

use super::*;

impl<'db> Scorer<'db> {
    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        // An MS2 scan with no selected precursor cannot be searched: there is no
        // precursor mass to query the index with. Report no match instead of
        // panicking, so a library caller that skips `is_searchable` sees an empty
        // result rather than an abort.
        let Some(precursor) = query.precursors.first() else {
            return Vec::new();
        };

        let hits = self.initial_hits(query, precursor);
        let mut features = Vec::with_capacity(self.report_psms);
        self.build_features(query, precursor, &hits, self.report_psms, &mut features);
        features
    }

    /// Score spectrum-specific mass-offset hypotheses.
    ///
    /// Each hypothesis searches indexed peptides at `precursor_mass -
    /// offset.mass()` (within [`Scorer::precursor_tol`], evaluated at
    /// `precursor_mass`), counting preliminary matches unshifted and at every
    /// `preliminary_shifts` entry, and then fully scores the best
    /// `preliminary_k` candidates with the offset placed at each compatible
    /// site. Configured offsets ([`IndexedDatabase::mass_offsets`]) are the
    /// special case of one fixed hypothesis per offset shared by every
    /// spectrum; this entry point lets a caller derive offsets from the
    /// spectrum itself, for example a crosslinked partner chain.
    ///
    /// Returns, for every hypothesis, its matches sorted by hyperscore and
    /// truncated to `report_k`. The spectrum is scored as given: search-time
    /// recalibration is applied by the caller through [`Scorer::recalibrated`].
    pub fn score_offset_hypotheses(
        &self,
        query: &ProcessedSpectrum,
        hypotheses: &[OffsetHypothesis],
        preliminary_k: usize,
        report_k: usize,
    ) -> Vec<Vec<OffsetMatch>> {
        let mut indexes: Vec<(u8, FragmentMatchIndex)> = Vec::new();
        hypotheses
            .iter()
            .map(|hypothesis| {
                let charge =
                    max_fragment_charge(self.max_fragment_charge, hypothesis.precursor_charge);
                if !indexes.iter().any(|(c, _)| *c == charge) {
                    indexes.push((charge, FragmentMatchIndex::new(query, charge)));
                }
                let fragment_index = &indexes
                    .iter()
                    .find(|(c, _)| *c == charge)
                    .expect("index was inserted")
                    .1;
                self.score_offset_hypothesis(
                    query,
                    hypothesis,
                    fragment_index,
                    preliminary_k,
                    report_k,
                )
            })
            .collect()
    }

    fn score_offset_hypothesis(
        &self,
        query: &ProcessedSpectrum,
        hypothesis: &OffsetHypothesis,
        fragment_index: &FragmentMatchIndex,
        preliminary_k: usize,
        report_k: usize,
    ) -> Vec<OffsetMatch> {
        let (mass, tolerance) = offset_query(
            hypothesis.precursor_mass,
            self.precursor_tol,
            hypothesis.offset.mass(),
        );
        let mut hits = self.preliminary_hits(
            fragment_index,
            mass,
            tolerance,
            &hypothesis.preliminary_shifts,
            PreScore {
                precursor_charge: hypothesis.precursor_charge,
                ..Default::default()
            },
        );
        let k = preliminary_k.min(hits.preliminary.len());
        bounded_min_heapify(&mut hits.preliminary, k);
        hits.preliminary.truncate(k);

        let mut matches = Vec::new();
        for pre in hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
        {
            let base = &self.db[pre.peptide];
            for site in self.db.mass_offset_sites(base, &hypothesis.offset) {
                let peptide = base.with_mass_offset(site, &hypothesis.offset.definition);
                let (score, fragments, _) =
                    self.score_peptide(query, &peptide, pre, true, fragment_index);
                if score.matched_b + score.matched_y < self.min_matched_peaks {
                    continue;
                }
                matches.push(OffsetMatch {
                    peptide_idx: pre.peptide,
                    peptide,
                    site,
                    preliminary_matched: pre.matched,
                    hyperscore: score.hyperscore,
                    matched_b: score.matched_b,
                    matched_y: score.matched_y,
                    summed_b: score.summed_b,
                    summed_y: score.summed_y,
                    longest_b: score.longest_b as u32,
                    longest_y: score.longest_y as u32,
                    fragment_ppm: score.ppm_difference,
                    fragments: fragments.expect("fragments were requested"),
                });
            }
        }
        matches.sort_by(|a, b| {
            b.hyperscore
                .total_cmp(&a.hyperscore)
                .then_with(|| a.peptide_idx.cmp(&b.peptide_idx))
                .then_with(|| a.site.cmp(&b.site))
        });
        matches.truncate(report_k);
        matches
    }
}
