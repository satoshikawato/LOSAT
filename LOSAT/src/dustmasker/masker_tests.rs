use super::*;

#[test]
fn runs_of_n_at_the_ends_always_count_and_inner_ones_only_when_longer_than_max() {
    assert_eq!(
        long_n_runs(2, b"NNACNNTNNNG"),
        vec![(0, 1), (7, 9)],
        "the run of 2 inside is not longer than 2"
    );
    assert_eq!(long_n_runs(64, b"ACGTN"), vec![(4, 4)]);
    assert_eq!(long_n_runs(u32::MAX, b"NNNN"), vec![(0, 3)]);
    assert!(long_n_runs(0, b"ACGT").is_empty());
}

#[test]
fn merges_take_the_linker_as_an_unsigned_sum() {
    let mut list = vec![(0, 9)];
    insert_merge(&mut list, (10, 20), 1);
    assert_eq!(list, vec![(0, 20)]);
    // -1 as Uint4: 20 + 4294967295 wraps to 19.
    insert_merge(&mut list, (19, 30), u32::MAX);
    assert_eq!(list, vec![(0, 30)]);
    insert_merge(&mut list, (40, 50), 1);
    assert_eq!(list, vec![(0, 30), (40, 50)]);
}
