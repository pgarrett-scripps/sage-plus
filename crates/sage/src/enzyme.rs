use fnv::FnvHashSet;
use regex::Regex;
use std::sync::Arc;

use crate::ambiguous_residues;
use crate::mass::VALID_AA;
use crate::sequence::{PeptideSequence, ProteinSequence};

#[derive(Clone, PartialOrd, Ord, Debug, Default)]
/// An enzymatic digest
///
/// # Important invariant about [`Digest`]:
/// * two digests are equal if and only if their sequences and position are equal
///   i.e., decoy status is ignored for equality and hashing
pub struct Digest {
    /// Is this a decoy peptide?
    pub decoy: bool,
    /// Semi-enzymatic?
    pub semi_enzymatic: bool,
    /// Cleaved peptide sequence
    pub sequence: PeptideSequence,
    /// Protein accession
    pub protein: Arc<str>,
    /// Zero-based start offset of this peptide in the source protein.
    pub protein_start: Option<u32>,
    /// Amino acid immediately before the peptide in the source protein.
    pub prev_aa: Option<u8>,
    /// Amino acid immediately after the peptide in the source protein.
    pub next_aa: Option<u8>,
    /// Missed cleavages
    pub missed_cleavages: u8,
    /// Is this an N-terminal peptide of the protein?
    pub position: Position,
    /// Source protein of a digest expanded from ambiguous FASTA residues
    /// (B, X or Z). Its `sequence` is then a standalone allocation holding
    /// one expansion, and `protein_start` locates the ambiguous span as
    /// written in this protein. `None` for every other digest.
    pub expanded_from: Option<ProteinSequence>,
}

#[derive(Clone)]
pub struct DigestGroup {
    pub reference: Digest,
    pub origins: Vec<ProteinOccurrence>,
}

#[derive(Clone)]
pub struct ProteinOccurrence {
    pub protein: Arc<str>,
    pub start: Option<u32>,
    pub prev_aa: Option<u8>,
    pub next_aa: Option<u8>,
    /// Shared full source protein, when known, for motif site rules. This is a
    /// reference to the FASTA allocation, not a copy, and it does not take part
    /// in equality, ordering, or reporting.
    pub source: Option<ProteinSequence>,
    /// The occurrence was digested from the protein with its initiator
    /// methionine clipped (see [`metap_clips`]), so it starts the mature
    /// protein at offset 1. Set only when clipping was enabled for an
    /// enzymatic digest; it does not take part in equality or ordering.
    pub met_clipped: bool,
}

impl ProteinOccurrence {
    /// Occurrence of a digest, keeping its source protein when the digest
    /// sequence is a span of that protein at `protein_start`.
    ///
    /// A sequence that spans its whole allocation is ambiguous: it may be a
    /// whole-protein digest or a standalone sequence, such as a reversed
    /// [`Digest::reverse`] decoy or a preview peptide, whose allocation is not
    /// a protein. Such digests get no source here; use
    /// [`Self::of_protein_digest`] for digests cut from a protein.
    pub fn of(digest: &Digest) -> Self {
        Self::with_source(digest, false)
    }

    /// Occurrence of a digest cut from its protein by [`EnzymeParameters`],
    /// so the sequence always views the protein allocation. This keeps the
    /// source also when the peptide is the whole protein, which matters for a
    /// single-peptide FASTA decoy that must be matched literally.
    pub fn of_protein_digest(digest: &Digest) -> Self {
        Self::with_source(digest, true)
    }

    fn with_source(digest: &Digest, protein_backed: bool) -> Self {
        let source = match &digest.expanded_from {
            // An expansion views its own allocation; its source is the
            // protein the ambiguous span was cut from.
            Some(protein) => digest.protein_start.map(|_| protein.clone()),
            None => {
                let (storage, offset) = digest.sequence.source();
                (digest.protein_start == Some(offset)
                    && (protein_backed || storage.as_bytes().len() > digest.sequence.len()))
                .then_some(storage)
            }
        };
        Self {
            protein: digest.protein.clone(),
            start: digest.protein_start,
            prev_aa: digest.prev_aa,
            next_aa: digest.next_aa,
            source,
            met_clipped: digest.protein_start == Some(1)
                && matches!(digest.position, Position::Nterm | Position::Full),
        }
    }

    fn key(&self) -> (&str, Option<u32>, Option<u8>, Option<u8>) {
        (&self.protein, self.start, self.prev_aa, self.next_aa)
    }
}

impl std::fmt::Debug for ProteinOccurrence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ProteinOccurrence")
            .field("protein", &self.protein)
            .field("start", &self.start)
            .field("prev_aa", &self.prev_aa)
            .field("next_aa", &self.next_aa)
            .field("source", &self.source.is_some())
            .field("met_clipped", &self.met_clipped)
            .finish()
    }
}

impl PartialEq for ProteinOccurrence {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}

impl Eq for ProteinOccurrence {}

impl PartialOrd for ProteinOccurrence {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ProteinOccurrence {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.key().cmp(&other.key())
    }
}

fn normalize_origins(group: &mut DigestGroup) {
    group.origins.sort_unstable();
    group.origins.dedup();
}

/// Group digests of any origin. See [`ProteinOccurrence::of`].
pub fn group_digests(digests: Vec<Digest>) -> Vec<DigestGroup> {
    group_digests_by(digests, ProteinOccurrence::of)
}

/// Group digests cut from FASTA proteins, keeping every occurrence's source
/// protein. See [`ProteinOccurrence::of_protein_digest`].
pub fn group_protein_digests(digests: Vec<Digest>) -> Vec<DigestGroup> {
    group_digests_by(digests, ProteinOccurrence::of_protein_digest)
}

fn group_digests_by(
    mut digests: Vec<Digest>,
    occurrence: fn(&Digest) -> ProteinOccurrence,
) -> Vec<DigestGroup> {
    if digests.is_empty() {
        return Vec::new();
    }
    let mut groups = Vec::new();
    // A total order, so the group reference (which supplies semi-enzymatic
    // and missed-cleavage state) does not depend on which other digests are
    // grouped alongside, e.g. in prefilter sequence buckets. Within a group
    // the most enzymatic occurrence comes first and becomes the reference.
    digests.sort_unstable_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then(a.decoy.cmp(&b.decoy))
            .then(a.sequence.cmp(&b.sequence))
            .then(a.semi_enzymatic.cmp(&b.semi_enzymatic))
            .then(a.missed_cleavages.cmp(&b.missed_cleavages))
            .then_with(|| a.protein.cmp(&b.protein))
            .then(a.protein_start.cmp(&b.protein_start))
    });
    let mut digests = digests.into_iter();
    let first = digests.next().expect("checked non-empty above");
    let first_origin = occurrence(&first);
    let mut curr_group = DigestGroup {
        reference: first,
        origins: vec![first_origin],
    };
    for digest in digests {
        if digest.decoy == curr_group.reference.decoy
            && digest.position == curr_group.reference.position
            && digest.sequence == curr_group.reference.sequence
        {
            curr_group.origins.push(occurrence(&digest));
        } else {
            normalize_origins(&mut curr_group);
            groups.push(curr_group);
            let origin = occurrence(&digest);
            curr_group = DigestGroup {
                reference: digest,
                origins: vec![origin],
            };
        }
    }
    normalize_origins(&mut curr_group);
    groups.push(curr_group);
    groups
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Default)]
pub enum Position {
    Nterm,
    Cterm,
    Full,
    #[default]
    Internal,
}

impl Digest {
    /// Generate an internal decoy sequence by reversing the sequence
    /// inside the first and last amino acids
    pub fn reverse(&self) -> Self {
        if self.decoy {
            return self.clone();
        }

        Digest {
            decoy: true,
            semi_enzymatic: self.semi_enzymatic,
            protein: self.protein.clone(),
            protein_start: self.protein_start,
            prev_aa: self.prev_aa,
            next_aa: self.next_aa,
            sequence: self.sequence.reversed_internal(),
            missed_cleavages: self.missed_cleavages,
            position: self.position,
            expanded_from: self.expanded_from.clone(),
        }
    }
}

impl PartialEq for Digest {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.position == other.position
    }
}

impl Eq for Digest {}

impl std::hash::Hash for Digest {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sequence.hash(state);
        self.position.hash(state);
    }
}

/// Longest peptide Sage can search: modification sites are stored as one
/// byte per residue position.
pub const MAX_PEPTIDE_LEN: usize = u8::MAX as usize;

#[derive(Clone)]
pub struct EnzymeParameters {
    /// Number of missed cleavages to produce; `None` for unlimited, bounded
    /// only by `max_len`
    pub missed_cleavages: Option<u8>,
    /// Inclusive
    pub min_len: usize,
    /// Inclusive
    pub max_len: usize,
    pub enzyme: Option<Enzyme>,
    /// Also digest each protein as if methionine aminopeptidase removed its
    /// initiator methionine. See [`metap_clips`]. Ignored by non-specific
    /// digests.
    pub clip_n_term_met: bool,
    /// When set, each digest with ambiguous residues (B, X or Z) is replaced
    /// by one digest per residue combination, and digests with more
    /// combinations than this are dropped. When `None`, such digests are
    /// returned as written and dropped when converted to peptides.
    pub ambiguous_variants: Option<usize>,
}

/// Residues that let methionine aminopeptidase (MetAP) remove an initiator
/// methionine when they sit at the second position of a protein.
pub const METAP_SECOND_RESIDUES: &[u8] = b"GASTCPV";

/// Does MetAP remove the initiator methionine of `protein`? True when the
/// protein starts with M followed by a small residue in
/// [`METAP_SECOND_RESIDUES`]. The clipped protein starts at offset 1.
pub fn metap_clips(protein: &[u8]) -> bool {
    matches!(protein, [b'M', second, ..] if METAP_SECOND_RESIDUES.contains(second))
}

#[derive(Clone)]
pub struct Enzyme {
    // Skip cleaving if the site is followed by one of these AAs
    pub skip_suffix: [bool; 26],
    // Regex for matching cleavage sites
    regex: Regex,
    // Cleave at c-terminal?
    pub c_terminal: bool,
    // Semi-enzymatic cleavage?
    pub semi_enzymatic: bool,
}

#[derive(Clone)]
pub struct DigestSite {
    // Range defining cleavage position
    pub site: std::ops::Range<usize>,
    // Number of missed cleavages
    pub missed_cleavages: u8,

    pub semi_enzymatic: bool,
}

impl Enzyme {
    /// Check that `cleave` (`cleave_at`) and `skip_suffix` (`restrict`) name
    /// only residues Sage can digest. Ambiguity codes such as B, Z, J and X,
    /// lowercase letters and other symbols are rejected; `cleave` may also be
    /// empty (no digestion) or `$` (cleave at the C-terminus only).
    pub fn validate_residues(cleave: &str, skip_suffix: &str) -> Result<(), String> {
        let invalid = |residues: &str| {
            residues
                .chars()
                .filter(|x| !x.is_ascii() || !VALID_AA.contains(&(*x as u8)))
                .collect::<String>()
        };
        let bad = invalid(cleave);
        if !bad.is_empty() && cleave != "$" {
            return Err(format!(
                "`database.enzyme.cleave_at` contains unsupported residues `{bad}` in `{cleave}`; \
                 use one-letter codes from {}, `$`, or an empty string",
                std::str::from_utf8(&VALID_AA).unwrap_or_default()
            ));
        }
        let bad = invalid(skip_suffix);
        if !bad.is_empty() {
            return Err(format!(
                "`database.enzyme.restrict` contains unsupported residues `{bad}` in `{skip_suffix}`; \
                 use one-letter codes from {}",
                std::str::from_utf8(&VALID_AA).unwrap_or_default()
            ));
        }
        Ok(())
    }

    /// Build an enzyme, returning an error for unsupported residues instead of
    /// panicking. `Ok(None)` means no digestion (empty `cleave`).
    pub fn try_new(
        cleave: &str,
        skip_suffix: &str,
        c_terminal: bool,
        semi_enzymatic: bool,
    ) -> Result<Option<Self>, String> {
        Self::validate_residues(cleave, skip_suffix)?;
        Ok(Self::build(cleave, skip_suffix, c_terminal, semi_enzymatic))
    }

    /// Build an enzyme from validated residues.
    ///
    /// # Panics
    ///
    /// Panics when `cleave` or `skip_suffix` contains unsupported residues;
    /// use [`Enzyme::try_new`] for untrusted input.
    pub fn new(
        cleave: &str,
        skip_suffix: &str,
        c_terminal: bool,
        semi_enzymatic: bool,
    ) -> Option<Self> {
        Self::try_new(cleave, skip_suffix, c_terminal, semi_enzymatic)
            .unwrap_or_else(|error| panic!("{error}"))
    }

    fn build(
        cleave: &str,
        skip_suffix: &str,
        c_terminal: bool,
        semi_enzymatic: bool,
    ) -> Option<Self> {
        // At this point, cleave can be three things: empty, "$", or a string of valid AA's
        match cleave {
            "" => None,
            "$" => Some(Enzyme {
                regex: Regex::new("$").unwrap(),
                skip_suffix: [false; 26],
                // Allowing this to be set to false could cause unexpected behavior
                c_terminal: true,
                // Do not allow strange behavior
                semi_enzymatic: false,
            }),
            _ => Some(Enzyme {
                regex: Regex::new(&format!("[{}]", cleave.replace('?', ""))).unwrap(),
                skip_suffix: {
                    let mut arr = [false; 26];
                    for b in skip_suffix.bytes() {
                        arr[(b - b'A') as usize] = true;
                    }
                    arr
                },
                c_terminal,
                semi_enzymatic,
            }),
        }
    }

    /// Whether this enzyme cuts the bond between residues `left` and `right`,
    /// by the rule [`Enzyme::cleavage_sites`] applies inside a protein: the
    /// residue on the cleaved side matches the cleavage set and `right` is not
    /// a restricted residue. The `$` (no-cleavage) enzyme never cuts a bond.
    pub fn cleaves_between(&self, left: u8, right: u8) -> bool {
        let site = if self.c_terminal { left } else { right };
        if !site.is_ascii() {
            return false;
        }
        let mut buffer = [0u8; 4];
        let site = (site as char).encode_utf8(&mut buffer);
        let matched = self.regex.find(site).is_some_and(|found| !found.is_empty());
        let restricted = right.is_ascii_uppercase() && self.skip_suffix[(right - b'A') as usize];
        matched && !restricted
    }

    pub fn cleavage_sites(&self, sequence: &str) -> Vec<DigestSite> {
        let mut sites = Vec::new();
        let mut left = 0;
        for mat in self.regex.find_iter(sequence) {
            let right = match self.c_terminal {
                true => mat.end(),
                false => mat.start(),
            };
            if sequence
                .as_bytes()
                .get(right)
                .is_some_and(|b| self.skip_suffix[(b - b'A') as usize])
            {
                continue;
            }
            sites.push(DigestSite {
                site: left..right,
                missed_cleavages: 0,
                semi_enzymatic: false,
            });
            left = right;
        }
        sites.push(DigestSite {
            site: left..sequence.len(),
            missed_cleavages: 0,
            semi_enzymatic: false,
        });
        sites
    }
}

impl EnzymeParameters {
    pub fn cleavage_sites(&self, sequence: &str) -> Vec<DigestSite> {
        match &self.enzyme {
            Some(enzyme) => enzyme.cleavage_sites(sequence),
            None => {
                // Perform a non-specific digest
                let mut v = Vec::new();
                for len in self.min_len..=self.max_len.min(sequence.len()) {
                    for i in 0..=sequence.len().saturating_sub(len) {
                        v.push(DigestSite {
                            site: i..i + len,
                            missed_cleavages: 0,
                            semi_enzymatic: false,
                        });
                    }
                }
                v
            }
        }
    }

    /// Longest span a digest site may have and still produce a peptide of at
    /// most `max_len` residues: initiator methionine clipping shortens
    /// N-terminal spans by one residue.
    fn max_span(&self) -> usize {
        self.max_len.saturating_add(self.clip_n_term_met as usize)
    }

    /// Join runs of adjacent cleavage sites into missed-cleavage sites, in
    /// order of increasing run length and then start. Each run grows one site
    /// at a time and stops once it can no longer yield a peptide within
    /// [`Self::max_span`]; a semi-enzymatic run can, while either its newly
    /// added end or its newly added start half still fits.
    fn missed_cleavage_sites(&self, sites: &mut Vec<DigestSite>, missed_cleavages: Option<u8>) {
        let max_span = self.max_span();
        let semi_enzymatic = self.is_semi_enzymatic();
        let limit = missed_cleavages.map_or(usize::MAX, usize::from);
        let mut missed_cleavage_sites = Vec::new();
        let mut open = (0..sites.len()).collect::<Vec<_>>();
        let mut missed = 0;
        while !open.is_empty() && missed <= limit {
            open.retain(|&first| {
                let Some(last) = sites.get(first + missed) else {
                    return false;
                };
                let (start, end) = (sites[first].site.start, last.site.end);
                let fits = match semi_enzymatic {
                    true => {
                        last.site.start - start <= max_span
                            || end - sites[first].site.end <= max_span
                    }
                    false => end - start <= max_span,
                };
                if fits {
                    missed_cleavage_sites.push(DigestSite {
                        site: start..end,
                        missed_cleavages: u8::try_from(missed).unwrap_or(u8::MAX),
                        semi_enzymatic: false,
                    });
                }
                fits
            });
            missed += 1;
        }
        sites.append(&mut missed_cleavage_sites);
    }

    fn is_semi_enzymatic(&self) -> bool {
        match &self.enzyme {
            Some(enzyme) => enzyme.semi_enzymatic,
            None => false,
        }
    }

    fn semi_enzymatic_sites(&self, sites: &mut Vec<DigestSite>) {
        let mut semi_enzymatic_sites = Vec::new();
        for site in sites.iter() {
            let start = site.site.start;
            let end = site.site.end;
            for cut_site in start..end {
                // Missed cleavages should actually be split across the sites,
                // but we rely on the fact that we generate missed cleavages in ascending
                // order, and eliminate previously seen digests
                semi_enzymatic_sites.push(DigestSite {
                    site: start..cut_site,
                    missed_cleavages: site.missed_cleavages,
                    semi_enzymatic: true,
                });
                semi_enzymatic_sites.push(DigestSite {
                    site: cut_site..end,
                    missed_cleavages: site.missed_cleavages,
                    semi_enzymatic: true,
                });
            }
        }
        sites.append(&mut semi_enzymatic_sites);
    }

    pub fn digest(&self, sequence: &str, protein: Arc<str>) -> Vec<Digest> {
        self.digest_with_custom_cleavages(sequence, protein, &[])
    }

    /// Perform the configured digest and add peptides anchored at the supplied
    /// protein-specific cleavage boundaries. Boundaries are zero-based offsets
    /// into `sequence` (one greater than the user-facing residue index).
    pub fn digest_with_custom_cleavages(
        &self,
        sequence: &str,
        protein: Arc<str>,
        custom_boundaries: &[usize],
    ) -> Vec<Digest> {
        let sequence: ProteinSequence = sequence.into();
        self.digest_protein_with_custom_cleavages(&sequence, protein, custom_boundaries)
    }

    /// Digest a sequence already owned by the FASTA database. Every emitted
    /// peptide retains only a shared storage pointer and a start/end range.
    pub fn digest_protein_with_custom_cleavages(
        &self,
        protein_sequence: &ProteinSequence,
        protein: Arc<str>,
        custom_boundaries: &[usize],
    ) -> Vec<Digest> {
        let sequence = protein_sequence.as_str();
        let n = sequence.len();
        let mut digests = Vec::new();
        let mut sites = self.cleavage_sites(sequence);
        let mut enzyme_boundaries = Vec::new();
        if self.enzyme.is_some() {
            enzyme_boundaries.push(0);
            enzyme_boundaries.extend(sites.iter().map(|site| site.site.end));
            enzyme_boundaries.sort_unstable();
            enzyme_boundaries.dedup();
        }
        // Allowing missed_cleavages with non-specific digest causes OOB panics
        // in the below indexing code
        let missed_cleavages = match self.enzyme {
            None => Some(0),
            _ => self.missed_cleavages,
        };

        if missed_cleavages != Some(0) {
            self.missed_cleavage_sites(&mut sites, missed_cleavages);
        }

        if self.is_semi_enzymatic() {
            self.semi_enzymatic_sites(&mut sites);
        }

        // A non-specific digest already contains every possible peptide. For
        // enzymatic and no-digest modes, add only peptides that terminate at a
        // nominated custom boundary and an ordinary enzyme/protein boundary.
        if !enzyme_boundaries.is_empty() {
            let max_span = self.max_span();
            let runs = self
                .missed_cleavages
                .map_or(usize::MAX, |missed| usize::from(missed) + 1);
            for &boundary in custom_boundaries {
                if boundary == 0 || boundary >= n {
                    continue;
                }
                let semi_enzymatic = enzyme_boundaries.binary_search(&boundary).is_err();
                for (missed, &start) in enzyme_boundaries
                    .iter()
                    .rev()
                    .filter(|&&candidate| candidate < boundary)
                    .take_while(|&&candidate| boundary - candidate <= max_span)
                    .take(runs)
                    .enumerate()
                {
                    sites.push(DigestSite {
                        site: start..boundary,
                        missed_cleavages: u8::try_from(missed).unwrap_or(u8::MAX),
                        semi_enzymatic,
                    });
                }
                for (missed, &end) in enzyme_boundaries
                    .iter()
                    .filter(|&&candidate| candidate > boundary)
                    .take_while(|&&candidate| candidate - boundary <= max_span)
                    .take(runs)
                    .enumerate()
                {
                    sites.push(DigestSite {
                        site: boundary..end,
                        missed_cleavages: u8::try_from(missed).unwrap_or(u8::MAX),
                        semi_enzymatic,
                    });
                }
            }
        }

        // Initiator methionine clipping: every site anchored at the protein
        // N-terminus is repeated from offset 1, where the clipped protein
        // starts. These peptides are protein N-terminal, keep the enzymatic
        // state of their unclipped counterpart, and stay in real protein
        // coordinates. A span that the digest already produced from offset 1
        // (a semi-enzymatic or custom cleavage after the Met) is emitted once,
        // as the clipped N-terminal peptide.
        let mut clipped_sites = Vec::new();
        if self.clip_n_term_met && self.enzyme.is_some() && metap_clips(sequence.as_bytes()) {
            // An enzyme boundary after the Met is not a missed cleavage of the
            // clipped protein.
            let cut_after_met = enzyme_boundaries.binary_search(&1).is_ok();
            clipped_sites.extend(
                sites
                    .iter()
                    .filter(|site| site.site.start == 0 && site.site.end > 1)
                    .map(|site| DigestSite {
                        site: 1..site.site.end,
                        missed_cleavages: site.missed_cleavages.saturating_sub(cut_after_met as u8),
                        semi_enzymatic: site.semi_enzymatic,
                    }),
            );
        }

        // Keep ranges unique while preserving repeated peptide sequences at
        // different protein positions for protein-coordinate PTM libraries.
        let mut seen = FnvHashSet::default();
        seen.extend(
            clipped_sites
                .iter()
                .map(|site| (site.site.start, site.site.end)),
        );
        let mut seen_clipped = FnvHashSet::default();

        let regular = sites.iter().map(|site| (site, false));
        let clipped = clipped_sites.iter().map(|site| (site, true));
        for (site, is_clipped) in regular.chain(clipped) {
            let start = site.site.start;
            let end = site.site.end;

            let peptide_sequence = match protein_sequence.peptide(start..end) {
                Some(peptide_sequence) => peptide_sequence,
                None => continue,
            };

            let len = peptide_sequence.len();

            let position = match (start == 0 || is_clipped, end == n) {
                (true, true) => Position::Full,
                (true, false) => Position::Nterm,
                (false, true) => Position::Cterm,
                (false, false) => Position::Internal,
            };

            let unique = match is_clipped {
                true => seen_clipped.insert((start, end)),
                false => seen.insert((start, end)),
            };
            if len >= self.min_len && len <= self.max_len && len > 0 && unique {
                digests.push(Digest {
                    sequence: peptide_sequence,
                    missed_cleavages: site.missed_cleavages,
                    decoy: false,
                    semi_enzymatic: site.semi_enzymatic,
                    position,
                    protein: protein.clone(),
                    protein_start: Some(start as u32),
                    prev_aa: start.checked_sub(1).map(|index| sequence.as_bytes()[index]),
                    next_aa: sequence.as_bytes().get(end).copied(),
                    expanded_from: None,
                });
            }
        }
        match self.ambiguous_variants {
            Some(max_variants) => expand_ambiguous(digests, protein_sequence, max_variants),
            None => digests,
        }
    }
}

/// Replace every digest containing B, X or Z by its expansions (see
/// [`ambiguous_residues::expand`]), dropping digests with more than
/// `max_variants` of them. Cleavage sites were already chosen from the
/// residues as written, so an X is never a trypsin site.
fn expand_ambiguous(
    digests: Vec<Digest>,
    protein: &ProteinSequence,
    max_variants: usize,
) -> Vec<Digest> {
    if !ambiguous_residues::is_ambiguous(protein.as_bytes()) {
        return digests;
    }
    let mut expanded = Vec::with_capacity(digests.len());
    for digest in digests {
        let count = ambiguous_residues::variant_count(digest.sequence.as_bytes());
        if count == 1 {
            expanded.push(digest);
            continue;
        }
        if count > max_variants {
            continue;
        }
        for variant in ambiguous_residues::expand(digest.sequence.as_bytes()) {
            expanded.push(Digest {
                sequence: variant.into(),
                expanded_from: Some(protein.clone()),
                ..digest.clone()
            });
        }
    }
    expanded
}

#[cfg(test)]
#[path = "../tests/unit/enzyme.rs"]
mod test;
