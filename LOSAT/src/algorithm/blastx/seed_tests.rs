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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:49-59,62-63
// A negative window (a search strategy's) makes the table shorter than the query, down to one entry
// (the doubling starts at 1 and never runs for a sum of 1 or less), and starts the offset at the
// window; the entries start as `calloc` leaves them.
#[test]
fn negative_window_shortens_the_table_and_starts_the_offset_at_the_window() {
    for (query, window, length) in [
        (500, -40, 512),
        (500, -244, 256),
        (500, -499, 1),
        (500, -498, 2),
        (500, -500, 1),
        (500, -100_000, 1),
        (3, i32::MIN, 1),
    ] {
        let d = Diagonals::new(query, window).unwrap();
        assert_eq!(d.entries.len(), length, "query={query} window={window}");
        assert_eq!(
            (d.offset, d.window, d.wide_offset),
            (window, window, i64::from(window))
        );
        assert!(d.entries.iter().all(|x| x.last_hit() == 0 && x.flag() == 0));
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:101-103,171-176
// With a negative window the offset moves by the subject's length plus the window (down over
// subjects shorter than the window's size), and at the threshold it is the window again with every
// entry at `last_hit` -window and no flag.
#[test]
fn negative_window_offset_moves_by_length_plus_window_and_clears_at_the_threshold() {
    let mut d = Diagonals::new(500, -40).unwrap();
    d.finish_subject_wide(100);
    assert_eq!((d.wide_offset, d.offset), (20, 20));
    d.finish_subject_wide(10);
    assert_eq!((d.wide_offset, d.offset), (-10, -10));
    d.finish_subject_wide(0);
    assert_eq!((d.wide_offset, d.offset), (-50, -50));
    d.wide_offset = i64::from(i32::MAX / 4) - 1;
    d.finish_subject_wide(100);
    assert_eq!(d.wide_offset, i64::from(i32::MAX / 4) - 1 + 100 - 40);
    let mut entry = DiagStruct::default();
    entry.set_last_hit(7);
    entry.set_flag(1);
    d.entries[3] = entry;
    d.finish_subject_wide(100);
    assert_eq!((d.wide_offset, d.offset), (-40, -40));
    assert!(d
        .entries
        .iter()
        .all(|x| x.last_hit() == 40 && x.flag() == 0));
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107,171-176
// Where NCBI's `Int4` sum overflows (undefined: a window of `INT_MIN`, or a long run of short
// subjects) the offset is a 64-bit one; `offset` is its low 32 bits, the entries' `last_hit` of
// the window `INT_MIN` is 0 (its negation truncated to 31 bits), and nothing panics.
#[test]
fn negative_window_offset_is_widened_where_the_int4_sum_overflows() {
    let mut d = Diagonals::new(10, i32::MIN).unwrap();
    d.finish_subject_wide(0);
    assert_eq!(d.wide_offset, 2 * i64::from(i32::MIN));
    assert_eq!(d.offset, 0);
    d.finish_subject_wide(7);
    assert_eq!(d.wide_offset, 3 * i64::from(i32::MIN) + 7);
    d.wide_offset = i64::from(i32::MAX / 4);
    d.entries[0].set_flag(1);
    d.finish_subject_wide(1);
    assert_eq!((d.wide_offset, d.offset), (i64::from(i32::MIN), i32::MIN));
    assert!(d.entries.iter().all(|x| x.last_hit() == 0 && x.flag() == 0));
}
