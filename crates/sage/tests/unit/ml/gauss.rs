use super::*;

fn assert_matrix_close(actual: &Matrix, expected: &[f64]) {
    assert_eq!(actual.data.len(), expected.len());
    assert!(
        actual
            .data
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (actual - expected).abs() < 1e-5),
        "actual matrix was {:?}",
        actual.data
    );
}

#[test]
fn solves_a_known_linear_system() {
    let left = Matrix::new([2.0, 1.0, 5.0, 7.0], 2, 2);
    let right = Matrix::col_vector(vec![11.0, 13.0]);

    let solution = Gauss::solve(left, right).unwrap();

    assert_matrix_close(&solution, &[64.0 / 9.0, -29.0 / 9.0]);
}

#[test]
fn pivots_when_the_leading_diagonal_is_zero() {
    let left = Matrix::new([0.0, 2.0, 1.0, 3.0], 2, 2);
    let right = Matrix::col_vector(vec![4.0, 5.0]);

    let solution = Gauss::solve(left, right).unwrap();

    assert_matrix_close(&solution, &[-1.0, 2.0]);
}

#[test]
fn solves_multiple_right_hand_sides() {
    let left = Matrix::diagonal(2, 2.0);
    let right = Matrix::identity(2);

    let solution = Gauss::solve(left, right).unwrap();

    assert_eq!(solution.shape(), (2, 2));
    assert_matrix_close(&solution, &[0.5, 0.0, 0.0, 0.5]);
}
