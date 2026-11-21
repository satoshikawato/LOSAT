use super::*;

#[test]
fn xmerge_keeps_descending_order() {
    let mut a = vec![3, 9, 5, 2];
    midl_xmerge(&mut a, &[3, 8, 4, 1]);
    assert_eq!(&a[..7], &[6, 9, 8, 5, 4, 2, 1]);
    let mut e = midl_new();
    midl_xmerge(&mut e, &[2, 7, 3]);
    assert_eq!(&e[..3], &[2, 7, 3]);
}

#[test]
fn append_range_is_descending() {
    let mut a = vec![1, 50];
    midl_append_range(&mut a, 10, 3);
    assert_eq!(&a[..5], &[4, 50, 12, 11, 10]);
}

#[test]
fn search_finds_insert_position_in_descending_list() {
    let a = [4, 40, 30, 20, 10];
    assert_eq!(midl_search(&a, 30), 2);
    assert_eq!(midl_search(&a, 25), 3);
    assert_eq!(midl_search(&a, 50), 1);
    assert_eq!(midl_search(&a, 5), 5);
}
