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
