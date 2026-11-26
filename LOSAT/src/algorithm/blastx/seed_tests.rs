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
    d.finish_chunk_wide(100);
    assert_eq!((d.wide_offset, d.offset), (20, 20));
    d.finish_chunk_wide(10);
    assert_eq!((d.wide_offset, d.offset), (-10, -10));
    d.finish_chunk_wide(0);
    assert_eq!((d.wide_offset, d.offset), (-50, -50));
    d.wide_offset = i64::from(i32::MAX / 4) - 1;
    d.finish_chunk_wide(100);
    assert_eq!(d.wide_offset, i64::from(i32::MAX / 4) - 1 + 100 - 40);
    let mut entry = DiagStruct::default();
    entry.set_last_hit(7);
    entry.set_flag(1);
    d.entries[3] = entry;
    d.finish_chunk_wide(100);
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
    d.finish_chunk_wide(0);
    assert_eq!(d.wide_offset, 2 * i64::from(i32::MIN));
    assert_eq!(d.offset, 0);
    d.finish_chunk_wide(7);
    assert_eq!(d.wide_offset, 3 * i64::from(i32::MIN) + 7);
    d.wide_offset = i64::from(i32::MAX / 4);
    d.entries[0].set_flag(1);
    d.finish_chunk_wide(1);
    assert_eq!((d.wide_offset, d.offset), (i64::from(i32::MIN), i32::MIN));
    assert!(d.entries.iter().all(|x| x.last_hit() == 0 && x.flag() == 0));
}

/// BLASTX's table as the `-db_hard_mask` search moves it between sequences
/// (`search_db_masked.rs`: `TableStart::apply`, `TableEnd::apply`): NCBI's table at its one-thread
/// offset, the clear before the move, the moves after a sequence, and the move past Int4 an error.
#[test]
fn hard_mask_starts_and_ends_move_the_table() {
    use crate::algorithm::common::hard_chunk_walk::{HardTable, NcbiStart, TableEnd, TableStart};
    let window = 40;
    let mut d = Diagonals::new(500, window).unwrap();
    let mut written = DiagStruct::default();
    written.set_last_hit(5000);
    written.set_flag(1);
    d.entries.fill(written);
    d.offset = 5000;
    TableStart::Keep.apply(&mut d, window).unwrap();
    assert_eq!((d.offset, d.entries[3]), (5000, written));
    TableStart::Move {
        clear: false,
        by: 70,
    }
    .apply(&mut d, window)
    .unwrap();
    assert_eq!((d.offset, d.entries[3]), (5070, written));
    TableStart::Move {
        clear: true,
        by: 70,
    }
    .apply(&mut d, window)
    .unwrap();
    assert_eq!(d.offset, window + 70);
    assert!(d
        .entries
        .iter()
        .all(|&entry| entry == DiagStruct::clear(window)));
    TableEnd::Move(98).apply(&mut d, window).unwrap();
    assert_eq!(d.offset, window + 168);
    d.entries.fill(written);
    TableEnd::Clear.apply(&mut d, window).unwrap();
    assert_eq!(d.offset, window);
    assert!(d
        .entries
        .iter()
        .all(|&entry| entry == DiagStruct::clear(window)));
    for (start, entry) in [
        (
            NcbiStart {
                offset: -98,
                cleared: false,
            },
            DiagStruct::default(),
        ),
        (
            NcbiStart {
                offset: 5,
                cleared: true,
            },
            DiagStruct::clear(window),
        ),
    ] {
        d.entries.fill(written);
        TableStart::New(start).apply(&mut d, window).unwrap();
        assert_eq!(d.offset, start.offset);
        assert!(d.entries.iter().all(|&e| e == entry), "{start:?}");
        d.entries.fill(written);
        d.ncbi_table(start, window);
        assert_eq!(d.offset, start.offset);
        assert!(d.entries.iter().all(|&e| e == entry), "{start:?}");
    }
    d.offset = i32::MAX - 10;
    let error = TableEnd::Move(98).apply(&mut d, window).unwrap_err();
    assert_eq!(error.to_string(), "BLASTX diagonal offset exceeds Int4");
}

// The refusal of a batch whose words fill NCBI's overflow banks comes after the reports of the
// batches before and is LOSAT's refusal, not an engine error (line RG item X-short-bat-2).
#[test]
fn filled_overflow_banks_are_refused_after_the_batches_before() {
    let error = x_banks_filled(1.0, 5);
    assert!(error.is::<crate::blastinput::app::ErrorAfterBatchReports>());
    assert!(super::super::native::search_error_as_is(&error));
    assert_eq!(
        error.to_string(),
        "-threshold 1 with the compressed lookup table of -word_size 5 is not supported by LOSAT's BLASTX: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap"
    );
}
