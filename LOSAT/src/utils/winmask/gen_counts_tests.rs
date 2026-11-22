use super::*;

#[test]
fn the_unit_size_follows_the_genome_length() {
    assert_eq!(unit_size_for(1), 1);
    assert_eq!(unit_size_for(4), 1);
    assert_eq!(unit_size_for(20), 2);
    assert_eq!(unit_size_for(21), 2);
    assert_eq!(unit_size_for(200_000), 8);
    assert_eq!(unit_size_for(400_000), 9);
    assert_eq!(unit_size_for(2_000_000_000), 15);
    assert_eq!(unit_size_for((-5i64) as u64), 16);
}

#[test]
fn the_suffix_fits_the_memory() {
    assert_eq!(suffix_size(12, 1), 9);
    assert_eq!(suffix_size(12, 64), 12);
    assert_eq!(suffix_size(16, 1536), 14);
    assert_eq!(suffix_size(16, 16384), 16);
}
