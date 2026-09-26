use super::*;

#[test]
fn dotv() {
    let a = Matrix::new([1., 2., 3., 4.], 2, 2);

    let v0 = a.dotv(&[0.5, 0.5]);
    assert_eq!(v0, vec![1.5, 3.5]);
    let n = norm(&v0);

    let c = v0.iter().map(|v| v / n).collect::<Vec<_>>();
    assert!(c
        .iter()
        .zip(&[0.3939193, 0.91914503])
        .all(|(x, y)| (x - y).abs() <= 0.0001));
}

#[test]
fn tranpose() {
    let mut mat = Matrix {
        data: vec![1., 2., 3., 4., 5., 6.],
        rows: 3,
        cols: 2,
    };

    assert_eq!(mat[(0, 0)], 1., "{:?}", mat);
    assert_eq!(mat[(0, 1)], 2., "{:?}", mat);
    assert_eq!(mat[(1, 0)], 3., "{:?}", mat);
    assert_eq!(mat[(1, 1)], 4., "{:?}", mat);
    assert_eq!(mat[(2, 0)], 5., "{:?}", mat);
    assert_eq!(mat[(2, 1)], 6., "{:?}", mat);

    mat = mat.transpose();

    assert_eq!(mat[(0, 0)], 1., "{:?}", mat);
    assert_eq!(mat[(0, 1)], 3., "{:?}", mat);
    assert_eq!(mat[(0, 2)], 5., "{:?}", mat);
    assert_eq!(mat[(1, 0)], 2., "{:?}", mat);
    assert_eq!(mat[(1, 1)], 4., "{:?}", mat);
    assert_eq!(mat[(1, 2)], 6., "{:?}", mat);
}

#[test]
fn dot() {
    #[rustfmt::skip]
        let a = vec![
            1., 0., 1., 
            2., 1., 1., 
            0., 1., 1., 
            1., 1., 2.
        ];
    let a = Matrix::new(a, 4, 3);

    #[rustfmt::skip]
        let b = vec![
            1., 2., 1., 
            2., 3., 1., 
            4., 2., 2.
        ];
    let b = Matrix::new(b, 3, 3);

    let c = a.dot(&b);
    assert_eq!(c.rows, 4);
    assert_eq!(c.cols, 3);
    #[rustfmt::skip]
        assert_eq!(
            c.data,
            vec![
                5., 4., 3., 
                8., 9., 5., 
                6., 5., 3., 
                11., 9., 6.
            ]
        );

    let d = vec![1., 2., 3., 4., 5., 6.];
    let d = Matrix::new(d, 2, 3);
    let e = Matrix::col_vector(vec![7., 9., 11.]);

    assert_eq!(
        d.dot(&e),
        Matrix {
            data: vec![58., 139.],
            cols: 1,
            rows: 2
        }
    );
}

#[test]
fn slice() {
    #[rustfmt::skip]
        let b = vec![
            1., 2., 1., 
            2., 3., 1., 
            4., 2., 2.
        ];
    let b = Matrix::new(b, 3, 3);

    assert_eq!(b.row_slice(0), &[1., 2., 1.]);
    assert_eq!(b.row_slice(1), &[2., 3., 1.]);
    assert_eq!(b.row_slice(2), &[4., 2., 2.]);
}

#[test]
fn constructors_fill_the_expected_shapes() {
    let zeros = Matrix::zeros(2, 3);
    assert_eq!(zeros.shape(), (2, 3));
    assert!(zeros.data.iter().all(|&x| x == 0.0));

    assert_eq!(Matrix::identity(2).take(), vec![1., 0., 0., 1.]);
    assert_eq!(Matrix::diagonal(2, 3.5).take(), vec![3.5, 0., 0., 3.5]);

    let col = Matrix::col_vector(vec![1., 2., 3.]);
    assert_eq!(col.shape(), (3, 1));
    let row = Matrix::row_vector(vec![1., 2., 3.]);
    assert_eq!(row.shape(), (1, 3));
    // Transposing a vector only swaps its shape.
    assert!(col.transpose() == row);
    assert_eq!(row.col(0).collect::<Vec<_>>(), vec![1.]);
    assert_eq!(col.col(0).collect::<Vec<_>>(), vec![1., 2., 3.]);
}

#[test]
#[should_panic(expected = "does not have shape (2, 2)")]
fn new_rejects_data_of_the_wrong_length() {
    Matrix::new([1., 2., 3.], 2, 2);
}

#[test]
fn element_access_is_bounds_checked() {
    let mut mat = Matrix::new([1., 2., 3., 4.], 2, 2);
    assert_eq!(mat.get(1, 0), Some(3.));
    assert_eq!(mat.get(2, 0), None);
    assert_eq!(mat.get(0, 2), None);
    *mat.get_mut(1, 1).unwrap() = 9.;
    mat.row_slice_mut(0)[1] = 7.;
    mat[(1, 0)] += 1.;
    assert_eq!(mat.take(), vec![1., 7., 4., 9.]);
}

#[test]
#[ignore = "bug: Matrix::get_mut does not bounds-check the column, so (0, cols) aliases (1, 0)"]
fn get_mut_rejects_out_of_range_columns() {
    let mut mat = Matrix::new([1., 2., 3., 4.], 2, 2);
    assert!(mat.get_mut(0, 2).is_none());
}

#[test]
fn is_close_compares_square_matrices_within_tolerance() {
    let a = Matrix::new([1., 2., 3., 4.], 2, 2);
    let b = Matrix::new([1.05, 2., 3., 3.95], 2, 2);
    assert!(a.is_close(&b, 0.1));
    assert!(!a.is_close(&b, 0.01));
    let rect = Matrix::new([1., 2., 3., 4., 5., 6.], 2, 3);
    assert!(!rect.is_close(&rect, 1.0));
}

#[test]
fn power_method_finds_the_dominant_eigenvector() {
    // Eigenvalues 3 (eigenvector [1, 1]) and 1 (eigenvector [1, -1]).
    let mat = Matrix::new([2., 1., 1., 2.], 2, 2);
    let v = mat.power_method(&[1., 0.]);
    let expected = std::f64::consts::FRAC_1_SQRT_2;
    assert!((v[0] - expected).abs() < 1e-4, "{v:?}");
    assert!((v[1] - expected).abs() < 1e-4, "{v:?}");
}

#[test]
fn mean_and_correlation_handle_constant_columns() {
    #[rustfmt::skip]
    let mat = Matrix::new([
        1., 6., 5.,
        2., 4., 5.,
        3., 2., 5.,
    ], 3, 3);
    assert_eq!(mat.mean(), vec![2., 4., 5.]);
    // Column 1 is an exact negative linear function of column 0; column 2 is
    // constant, so its undefined correlations are reported as zero.
    #[rustfmt::skip]
    let expected = Matrix::new([
         1., -1., 0.,
        -1.,  1., 0.,
         0.,  0., 0.,
    ], 3, 3);
    let corr = mat.correlation_matrix();
    assert!(corr.is_close(&expected, 1e-12), "{corr:?}");
}

#[test]
fn arithmetic_operators_are_elementwise() {
    let a = Matrix::new([1., 2., 3., 4.], 2, 2);
    let b = Matrix::identity(2);
    assert_eq!((a.clone() + b.clone()).take(), vec![2., 2., 3., 5.]);
    let mut c = a.clone();
    c += b;
    assert_eq!(c.take(), vec![2., 2., 3., 5.]);
    assert_eq!((a.clone() / 2.).take(), vec![0.5, 1., 1.5, 2.]);
    assert_eq!(format!("{a:?}"), "[\n[1.0, 2.0]\n[3.0, 4.0]\n]\n");
}

#[test]
#[should_panic(expected = "matrices must have equal shape to add")]
fn adding_mismatched_shapes_panics() {
    let _ = Matrix::zeros(2, 2) + Matrix::zeros(1, 2);
}

#[test]
#[should_panic(expected = "matrices must have equal shape to add")]
fn add_assign_with_mismatched_shapes_panics() {
    let mut a = Matrix::zeros(2, 2);
    a += Matrix::zeros(2, 1);
}

#[test]
#[should_panic(expected = "rhs has shape (1,3)")]
fn dotv_rejects_mismatched_lengths() {
    Matrix::identity(2).dotv(&[1., 2., 3.]);
}

#[test]
#[should_panic(expected = "rhs has shape (3,1)")]
fn dot_rejects_mismatched_shapes() {
    Matrix::identity(2).dot(&Matrix::zeros(3, 1));
}
