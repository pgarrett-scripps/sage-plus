use std::fmt::Debug;

use quickcheck_macros::quickcheck;

use super::bounded_min_heapify;
use super::check_heap;

fn check<T: Ord + Clone + Debug>(mut data: Vec<T>, k: usize) {
    let k = k.min(data.len());
    let mut cloned = data.clone();
    // Stable sort the data
    cloned.sort_by(|a, b| b.cmp(a));

    bounded_min_heapify(&mut data, k);

    // Take the heap part, and sort it
    let top_k = &mut data[..k];

    // Check that heap property is maintained, or that k == length of the data
    assert!(check_heap(top_k) || k == cloned.len());

    top_k.sort_by(|a, b| b.cmp(a));
    assert_eq!(top_k, &mut cloned[..k]);
}

#[quickcheck]
fn run_quickcheck(data: Vec<i32>, k: usize) {
    check(data, k);
}

#[test]
fn smoke() {
    let asc = (0..500).collect::<Vec<_>>();
    let desc = (0..500).rev().collect::<Vec<_>>();
    check(asc, 50);
    check(desc, 50);
}
