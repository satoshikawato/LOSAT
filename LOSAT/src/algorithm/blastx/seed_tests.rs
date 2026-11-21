use super::*;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
#[test]
fn diagonal_overflow_threshold_clears_flags_and_negative_window() {
    let mut d = Diagonals::new(500, 40).unwrap();
    d.offset = i32::MAX / 4 - 1;
    let mut entry = DiagStruct::default();
    entry.set_last_hit(123);
    entry.set_flag(1);
    d.entries[0] = entry;
    d.finish_subject(83).unwrap();
    assert_eq!(d.offset, i32::MAX / 4 - 1 + 83 + 40);
    assert_eq!(d.entries[0], entry);
    d.finish_subject(83).unwrap();
    assert_eq!(d.offset, 40);
    assert!(d
        .entries
        .iter()
        .all(|x| x.last_hit() == -40 && x.flag() == 0));
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:49-59
#[test]
fn diagonal_size_power_and_int4_overflow_boundary() {
    assert_eq!(Diagonals::new(512, 0).unwrap().entries.len(), 512);
    assert_eq!(Diagonals::new(512, 1).unwrap().entries.len(), 1024);
    assert!(Diagonals::new(i32::MAX as usize, 1).is_err());
    assert!(Diagonals::new(usize::MAX, 0).is_err());
}
