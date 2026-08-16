use crate::database::PeptideIx;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::peptide::Peptide;
use rayon::prelude::*;

/// Bitmap-based preliminary search index.
///
/// Encodes theoretical fragment ions as packed bitsets (one bit per mass bin)
/// per peptide, sorted by precursor mass. Experimental spectra are likewise
/// encoded as a bitset at search time. Scoring is then a bitwise AND + popcount,
/// which is extremely fast and SIMD-friendly.
///
/// **Requires deisotoped spectra.** After deisotoping each peak carries a
/// resolved charge state, so `peak.mass` already stores the neutral
/// monoisotopic mass: M = (mz − H) × z.  Non-deisotoped spectra are not
/// supported.
///
/// Configuration:
/// - `bitmap_size`: number of `u64` words per peptide (default 30 → 1920 bins).
/// - Mass range: `[min_mass, max_mass)` spans `bitmap_size * 64` equal-width bins.
pub struct BitmapIndex {
    /// Precursor monoisotopic masses, one entry per peptide, sorted ascending.
    pub precursor_masses: Vec<f32>,
    /// Peptide indices corresponding to each entry in `precursor_masses`.
    pub peptide_indices: Vec<PeptideIx>,
    /// Packed forward (A/B/C) bitmaps; stride = `bitmap_size`.
    pub forward_bitmaps: Vec<u64>,
    /// Packed reverse (X/Y/Z) bitmaps; stride = `bitmap_size`.
    pub reverse_bitmaps: Vec<u64>,
    /// Number of `u64` words per peptide bitmap.
    pub bitmap_size: usize,
    /// Lower bound of the mass range mapped to the bitmap.
    pub min_mass: f32,
    /// Upper bound of the mass range mapped to the bitmap.
    pub max_mass: f32,
}

impl Default for BitmapIndex {
    fn default() -> Self {
        Self {
            precursor_masses: Vec::new(),
            peptide_indices: Vec::new(),
            forward_bitmaps: Vec::new(),
            reverse_bitmaps: Vec::new(),
            bitmap_size: 30,
            min_mass: 500.0,
            max_mass: 5000.0,
        }
    }
}

impl BitmapIndex {
    /// Build a `BitmapIndex` from a slice of peptides (already sorted by
    /// monoisotopic mass, as produced by `Parameters::build_from_peptides`).
    pub fn build(
        peptides: &[Peptide],
        ion_kinds: &[Kind],
        bitmap_size: usize,
        min_mass: f32,
        max_mass: f32,
    ) -> Self {
        let n = peptides.len();
        let total_bins = bitmap_size * 64;

        // Build bitmaps in parallel — each peptide gets its own (forward, reverse) pair.
        let bitmaps: Vec<(Vec<u64>, Vec<u64>)> = peptides
            .par_iter()
            .map(|peptide| {
                let mut fwd = vec![0u64; bitmap_size];
                let mut rev = vec![0u64; bitmap_size];
                for &kind in ion_kinds {
                    let target = match kind {
                        Kind::A | Kind::B | Kind::C => &mut fwd,
                        Kind::X | Kind::Y | Kind::Z => &mut rev,
                    };
                    for ion in IonSeries::new(peptide, kind) {
                        let mass = ion.monoisotopic_mass;
                        if let Some(bin) = mass_to_bin(mass, min_mass, max_mass, total_bins) {
                            set_bit(target, bin);
                        }
                    }
                }
                (fwd, rev)
            })
            .collect();

        let mut precursor_masses = Vec::with_capacity(n);
        let mut peptide_indices = Vec::with_capacity(n);
        let mut forward_bitmaps = Vec::with_capacity(n * bitmap_size);
        let mut reverse_bitmaps = Vec::with_capacity(n * bitmap_size);

        for (idx, (fwd, rev)) in bitmaps.into_iter().enumerate() {
            precursor_masses.push(peptides[idx].monoisotopic);
            peptide_indices.push(PeptideIx(idx as u32));
            forward_bitmaps.extend_from_slice(&fwd);
            reverse_bitmaps.extend_from_slice(&rev);
        }

        BitmapIndex {
            precursor_masses,
            peptide_indices,
            forward_bitmaps,
            reverse_bitmaps,
            bitmap_size,
            min_mass,
            max_mass,
        }
    }

    /// Total number of bins in the bitmap.
    #[inline]
    pub fn total_bins(&self) -> usize {
        self.bitmap_size * 64
    }

    /// Build the experimental bitmap from a deisotoped spectrum.
    ///
    /// Each peak's `mass` field is treated as a neutral monoisotopic mass
    /// (M = (mz − H) × z, already resolved by deisotoping).  The tolerance
    /// window `[lo, hi]` is converted to a bin range; all bins in that range
    /// are set to 1, handling the bin-edge case automatically.
    pub fn experimental_bitmap(&self, masses: &[f32], tol: Tolerance) -> Vec<u64> {
        let mut bitmap = vec![0u64; self.bitmap_size];
        let total_bins = self.total_bins();

        for &mass in masses {
            let (lo, hi) = tol.bounds(mass);
            let bin_lo = mass_to_bin_clamped(lo, self.min_mass, self.max_mass, total_bins);
            let bin_hi = mass_to_bin_clamped(hi, self.min_mass, self.max_mass, total_bins);
            for bin in bin_lo..=bin_hi {
                set_bit(&mut bitmap, bin);
            }
        }
        bitmap
    }

    /// Score `exp_bitmap` against the theoretical bitmap for peptide `i`
    /// (both forward and reverse), returning `(matched_forward, matched_reverse)`.
    #[inline]
    pub fn score_peptide(&self, exp_bitmap: &[u64], i: usize) -> (u16, u16) {
        let offset = i * self.bitmap_size;
        let fwd = bitmap_score(
            exp_bitmap,
            &self.forward_bitmaps[offset..offset + self.bitmap_size],
        );
        let rev = bitmap_score(
            exp_bitmap,
            &self.reverse_bitmaps[offset..offset + self.bitmap_size],
        );
        (fwd as u16, rev as u16)
    }

    /// Find all candidate peptide indices within the given precursor mass
    /// tolerance (adjusted for isotope error) and score them against
    /// `exp_bitmap`.
    ///
    /// Returns an iterator of `(matched_forward, matched_reverse, PeptideIx)`.
    pub fn search<'a>(
        &'a self,
        exp_bitmap: &'a [u64],
        precursor_mass: f32,
        precursor_tol: Tolerance,
        isotope_error: i8,
    ) -> impl Iterator<Item = (u16, u16, PeptideIx)> + 'a {
        use crate::mass::NEUTRON;
        let search_mass = precursor_mass - isotope_error as f32 * NEUTRON;
        let (lo, hi) = precursor_tol.bounds(search_mass);

        let start = self.precursor_masses.partition_point(|&m| m < lo);
        let end = self.precursor_masses.partition_point(|&m| m <= hi);

        (start..end).map(move |i| {
            let (fwd, rev) = self.score_peptide(exp_bitmap, i);
            (fwd, rev, self.peptide_indices[i])
        })
    }
}

/// Map a mass value to a bin index, returning `None` if out of range.
#[inline]
fn mass_to_bin(mass: f32, min_mass: f32, max_mass: f32, total_bins: usize) -> Option<usize> {
    if mass < min_mass || mass >= max_mass {
        return None;
    }
    let frac = (mass - min_mass) / (max_mass - min_mass);
    Some((frac * total_bins as f32) as usize)
}

/// Map a mass to a bin, clamping to `[0, total_bins - 1]`.
#[inline]
fn mass_to_bin_clamped(mass: f32, min_mass: f32, max_mass: f32, total_bins: usize) -> usize {
    if mass <= min_mass {
        return 0;
    }
    if mass >= max_mass {
        return total_bins - 1;
    }
    let frac = (mass - min_mass) / (max_mass - min_mass);
    ((frac * total_bins as f32) as usize).min(total_bins - 1)
}

/// Set a single bit in a packed `u64` slice.
#[inline]
fn set_bit(bitmap: &mut [u64], bin: usize) {
    bitmap[bin / 64] |= 1u64 << (bin % 64);
}

/// Count the number of bits set in the bitwise AND of two equal-length `u64` slices.
#[inline]
pub fn bitmap_score(exp: &[u64], theo: &[u64]) -> u32 {
    exp.iter()
        .zip(theo)
        .map(|(e, t)| (e & t).count_ones())
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mass_to_bin() {
        // min=500, max=5000, total_bins=1920
        let total_bins = 1920usize;
        assert!(mass_to_bin(499.9, 500.0, 5000.0, total_bins).is_none());
        assert!(mass_to_bin(5000.0, 500.0, 5000.0, total_bins).is_none());
        assert_eq!(mass_to_bin(500.0, 500.0, 5000.0, total_bins), Some(0));
        assert_eq!(
            mass_to_bin(4999.9, 500.0, 5000.0, total_bins),
            Some(total_bins - 1)
        );
    }

    #[test]
    fn test_set_bit_and_score() {
        let mut bitmap = vec![0u64; 2]; // 128 bins
        set_bit(&mut bitmap, 0);
        set_bit(&mut bitmap, 63);
        set_bit(&mut bitmap, 64);
        set_bit(&mut bitmap, 127);

        let mut other = vec![0u64; 2];
        set_bit(&mut other, 0);
        set_bit(&mut other, 127);

        assert_eq!(bitmap_score(&bitmap, &other), 2);
    }

    #[test]
    fn test_experimental_bitmap_tolerance() {
        let index = BitmapIndex {
            bitmap_size: 2,
            min_mass: 0.0,
            max_mass: 128.0, // 1 Da per bin (128 bins)
            ..BitmapIndex::default()
        };

        // Peak at mass=10.0; tolerance ±0.6 Da → bins 9, 10
        let masses = vec![10.0];
        let tol = Tolerance::Da(-0.6, 0.6);
        let bm = index.experimental_bitmap(&masses, tol);

        // bin 9 = word 0, bit 9; bin 10 = word 0, bit 10
        assert!(bm[0] & (1u64 << 9) != 0);
        assert!(bm[0] & (1u64 << 10) != 0);
        // bin 8 should NOT be set
        assert!(bm[0] & (1u64 << 8) == 0);
    }
}
