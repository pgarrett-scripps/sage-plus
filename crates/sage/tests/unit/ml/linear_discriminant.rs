use super::*;
use crate::ml::*;

#[test]
fn linear_discriminant() {
    let a = Matrix::new([1., 2., 3., 4.], 2, 2);
    let eigenvector = [0.4159736, 0.90937671];
    assert!(all_close(
        &a.power_method(&[0.54, 0.34]),
        &eigenvector,
        1E-5
    ));

    #[rustfmt::skip]
        let feats: [[f64; 4]; 8] = [
            [5., 4., 3., 2.],
            [4., 5., 4., 3.],
            [6., 3., 4., 5.],
            [1., 0., 2., 9.],
            [5., 4., 4., 3.],
            [2., 1., 1., 9.5],
            [1., 0., 2., 8.],
            [3., 2., -2., 10.],
        ];

    let lda = LinearDiscriminantAnalysis::train::<_, 4>(
        &feats,
        &[false, false, false, true, false, true, true, true],
        |row| *row,
    )
    .expect("error training LDA");

    let mut scores: Vec<f64> = feats.iter().map(|row| lda.score(row)).collect();
    let norm = norm(&scores);
    scores = scores.into_iter().map(|s| s / norm).collect();

    let expected = [
        0.49706043,
        0.48920177,
        0.48920177,
        -0.07209359,
        0.51204672,
        -0.02849527,
        -0.04924864,
        -0.06055943,
    ];

    assert!(
        all_close(&scores, &expected, 1E-8),
        "{:?} {:?}",
        scores,
        expected
    );
}
