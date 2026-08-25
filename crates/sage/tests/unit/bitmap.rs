use super::*;
use crate::enzyme::Digest;
use crate::ion_series::IonSeries;
use crate::modification::{ModificationDefinition, ModificationSpecificity, NeutralLossMode};
use std::{collections::HashMap, sync::Arc};

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

#[test]
fn required_neutral_loss_is_used_by_bitmap_index() {
    let modification = Arc::new(ModificationDefinition {
        mass: 20.0,
        name: None,
        neutral_losses: Arc::from([10.0]),
        neutral_loss_mode: NeutralLossMode::Required,
        channel_offsets: Arc::default(),
    });
    let peptide = Peptide::try_from(Digest {
        sequence: "AMK".into(),
        ..Digest::default()
    })
    .unwrap()
    .apply(
        &[(
            ModificationSpecificity::Residue(b'M'),
            modification,
            Some(1),
        )],
        &HashMap::default(),
        1,
        None,
    )
    .into_iter()
    .find(|peptide| peptide.modification_at(1) != 0.0)
    .unwrap();

    let retained = IonSeries::new(&peptide, Kind::B)
        .nth(1)
        .unwrap()
        .monoisotopic_mass;
    let loss = IonGroupSeries::new(&peptide, Kind::B)
        .nth(1)
        .unwrap()
        .variants[0]
        .monoisotopic_mass;
    let index = BitmapIndex::build(&[peptide], &[Kind::B], 1024, 0.0, 1000.0);

    let loss_bitmap = index.experimental_bitmap(&[loss], Tolerance::Da(0.0, 0.0));
    assert_eq!(index.score_peptide(&loss_bitmap, 0).0, 1);
    let retained_bitmap = index.experimental_bitmap(&[retained], Tolerance::Da(0.0, 0.0));
    assert_eq!(index.score_peptide(&retained_bitmap, 0).0, 0);
}
