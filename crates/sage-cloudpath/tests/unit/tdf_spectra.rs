use super::*;

#[test]
fn group_and_sum_merges_equal_tof_indices() {
    let (tofs, ints) = group_and_sum(vec![5, 3, 5, 1], vec![1, 2, 3, 4]);
    assert_eq!(tofs, vec![1, 3, 5]);
    assert_eq!(ints, vec![4, 2, 4]);
}

#[test]
fn smooth_adds_neighbors_within_window() {
    let smoothed = smooth(&[10, 11, 20], &[1, 2, 4], 1);
    assert_eq!(smoothed, vec![3, 3, 4]);
}

#[test]
fn centroid_keeps_local_maxima() {
    let (tofs, ints) = centroid(&[10, 11, 12, 30], &[1, 5, 2, 7], 1);
    assert_eq!(tofs, vec![11, 30]);
    assert_eq!(ints, vec![5, 7]);
}

#[test]
fn centroid_drops_the_later_point_on_ties() {
    let (tofs, _) = centroid(&[10, 11], &[3, 3], 1);
    assert_eq!(tofs, vec![10]);
}

#[test]
fn isolation_window_round_trips_center_and_width() {
    let window = IsolationWindow::from_center(500.0, 25.0, 30.0);
    assert_eq!(window.center(), 500.0);
    assert_eq!(window.width(), 25.0);
    let bounds = IsolationWindow::from_bounds(400.0, 410.0, 0.0);
    assert_eq!(bounds.center(), 405.0);
}
