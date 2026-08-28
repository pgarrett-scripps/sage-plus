use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt::{Debug, Formatter};
use std::hash::{Hash, Hasher};
use std::ops::{Deref, Range};
use std::sync::{Arc, OnceLock};

/// One protein-sized allocation shared by every peptide span derived from it.
#[derive(Debug, PartialEq, Eq)]
struct SequenceStorage {
    bytes: Box<[u8]>,
}

/// An immutable protein sequence that can cheaply produce peptide spans.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProteinSequence {
    storage: Arc<SequenceStorage>,
}

impl ProteinSequence {
    pub fn peptide(&self, range: Range<usize>) -> Option<PeptideSequence> {
        if range.start > range.end || range.end > self.storage.bytes.len() {
            return None;
        }
        Some(PeptideSequence {
            storage: self.storage.clone(),
            start: u32::try_from(range.start).ok()?,
            end: u32::try_from(range.end).ok()?,
        })
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.storage.bytes
    }

    pub fn as_str(&self) -> &str {
        // FASTA parsing and the string constructors only admit UTF-8 input.
        std::str::from_utf8(self.as_bytes()).expect("protein sequence is not valid UTF-8")
    }
}

impl From<String> for ProteinSequence {
    fn from(sequence: String) -> Self {
        Self {
            storage: Arc::new(SequenceStorage {
                bytes: sequence.into_bytes().into_boxed_slice(),
            }),
        }
    }
}

impl From<&str> for ProteinSequence {
    fn from(sequence: &str) -> Self {
        sequence.to_owned().into()
    }
}

impl Deref for ProteinSequence {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        self.as_str()
    }
}

impl AsRef<str> for ProteinSequence {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

/// A content-addressable view into shared immutable sequence storage.
///
/// Equality, ordering, and hashing use the viewed residues. The backing
/// allocation and offsets are representation details, so identical peptides
/// from different proteins remain identical database keys.
#[derive(Clone)]
pub struct PeptideSequence {
    storage: Arc<SequenceStorage>,
    start: u32,
    end: u32,
}

impl PeptideSequence {
    pub fn as_bytes(&self) -> &[u8] {
        &self.storage.bytes[self.start as usize..self.end as usize]
    }

    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).expect("peptide sequence is not valid UTF-8")
    }

    pub fn starts_with(&self, prefix: &str) -> bool {
        self.as_bytes().starts_with(prefix.as_bytes())
    }

    pub fn reversed_internal(&self) -> Self {
        let mut sequence = self.as_bytes().to_vec();
        let last = sequence.len().saturating_sub(1);
        if last > 1 {
            sequence[1..last].reverse();
        }
        sequence.into()
    }

    #[cfg(test)]
    pub(crate) fn shares_storage_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.storage, &other.storage)
    }

    #[cfg(test)]
    pub(crate) fn storage_len(&self) -> usize {
        self.storage.bytes.len()
    }
}

impl Default for PeptideSequence {
    fn default() -> Self {
        static EMPTY: OnceLock<Arc<SequenceStorage>> = OnceLock::new();
        Self {
            storage: EMPTY
                .get_or_init(|| {
                    Arc::new(SequenceStorage {
                        bytes: Box::default(),
                    })
                })
                .clone(),
            start: 0,
            end: 0,
        }
    }
}

impl From<Vec<u8>> for PeptideSequence {
    fn from(sequence: Vec<u8>) -> Self {
        let end = u32::try_from(sequence.len()).expect("peptide sequence exceeds u32::MAX bytes");
        Self {
            storage: Arc::new(SequenceStorage {
                bytes: sequence.into_boxed_slice(),
            }),
            start: 0,
            end,
        }
    }
}

impl From<Box<[u8]>> for PeptideSequence {
    fn from(sequence: Box<[u8]>) -> Self {
        Vec::from(sequence).into()
    }
}

impl From<Arc<[u8]>> for PeptideSequence {
    fn from(sequence: Arc<[u8]>) -> Self {
        sequence.as_ref().to_vec().into()
    }
}

impl From<String> for PeptideSequence {
    fn from(sequence: String) -> Self {
        sequence.into_bytes().into()
    }
}

impl From<&str> for PeptideSequence {
    fn from(sequence: &str) -> Self {
        sequence.as_bytes().to_vec().into()
    }
}

impl From<&[u8]> for PeptideSequence {
    fn from(sequence: &[u8]) -> Self {
        sequence.to_vec().into()
    }
}

impl<const N: usize> From<&[u8; N]> for PeptideSequence {
    fn from(sequence: &[u8; N]) -> Self {
        sequence.as_slice().into()
    }
}

impl Deref for PeptideSequence {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.as_bytes()
    }
}

impl AsRef<[u8]> for PeptideSequence {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl Borrow<[u8]> for PeptideSequence {
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl Debug for PeptideSequence {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match std::str::from_utf8(self.as_bytes()) {
            Ok(sequence) => formatter
                .debug_tuple("PeptideSequence")
                .field(&sequence)
                .finish(),
            Err(_) => formatter
                .debug_tuple("PeptideSequence")
                .field(&self.as_bytes())
                .finish(),
        }
    }
}

impl std::fmt::Display for PeptideSequence {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&String::from_utf8_lossy(self.as_bytes()))
    }
}

impl PartialEq for PeptideSequence {
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl Eq for PeptideSequence {}

impl PartialEq<str> for PeptideSequence {
    fn eq(&self, other: &str) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<&str> for PeptideSequence {
    fn eq(&self, other: &&str) -> bool {
        self == *other
    }
}

impl PartialEq<PeptideSequence> for str {
    fn eq(&self, other: &PeptideSequence) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<PeptideSequence> for &str {
    fn eq(&self, other: &PeptideSequence) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<String> for PeptideSequence {
    fn eq(&self, other: &String) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialEq<PeptideSequence> for String {
    fn eq(&self, other: &PeptideSequence) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl PartialOrd for PeptideSequence {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PeptideSequence {
    fn cmp(&self, other: &Self) -> Ordering {
        self.as_bytes().cmp(other.as_bytes())
    }
}

impl Hash for PeptideSequence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn spans_share_protein_storage_but_compare_by_content() {
        let protein: ProteinSequence = "PEPTIDEPEPTIDE".into();
        let first = protein.peptide(0..7).unwrap();
        let second = protein.peptide(7..14).unwrap();
        let owned: PeptideSequence = "PEPTIDE".into();

        assert!(first.shares_storage_with(&second));
        assert_eq!(first.storage_len(), 14);
        assert_eq!(first, second);
        assert_eq!(first, owned);

        let sequences = HashSet::from([first, second, owned]);
        assert_eq!(sequences.len(), 1);
    }

    #[test]
    fn peptide_sequence_remains_two_machine_words() {
        assert_eq!(
            std::mem::size_of::<PeptideSequence>(),
            2 * std::mem::size_of::<usize>()
        );
    }
}
