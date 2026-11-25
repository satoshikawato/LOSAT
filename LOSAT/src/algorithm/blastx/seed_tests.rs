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
