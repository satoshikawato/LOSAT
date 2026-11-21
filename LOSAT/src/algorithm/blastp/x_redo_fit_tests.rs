use super::*;

#[test]
fn a_fitted_vector_keeps_its_elements_in_order() {
    let mut hsps = Vec::with_capacity(4);
    hsps.extend([3u32, 1, 2]);
    hsps.remove(1);
    x_fit_redone_hsps(&mut hsps);
    assert_eq!(hsps, [3, 2]);
    assert_eq!(hsps.capacity(), hsps.len());
}
