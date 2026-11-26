use super::*;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60,62-71
#[test]
fn s31_is_the_31_bit_signed_field() {
    assert_eq!(s31(0), 0);
    assert_eq!(s31(5), 5);
    assert_eq!(s31(-5), -5);
    assert_eq!(s31((1 << 30) - 1), (1 << 30) - 1);
    assert_eq!(s31(1 << 30), -(1 << 30));
    assert_eq!(s31(-(1 << 30)), -(1 << 30));
    assert_eq!(s31(-(1 << 30) - 1), (1 << 30) - 1);
    // `-INT_MIN` widened (NCBI's overflow) is 2^31: its low 31 bits are 0.
    assert_eq!(s31(1 << 31), 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,109-157
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:670-671
#[test]
fn fresh_array_with_window_minus_40_skips_words_below_40() {
    let table = NegTable::new(-40, Some(64), 0);
    assert_eq!(table.offset(), -40);
    for s_off in 0..40 {
        assert_eq!(table.seed_pos(s_off, 3), None, "s_off {s_off}");
    }
    assert_eq!(table.seed_pos(40, 3), Some(0));
    assert_eq!(table.seed_pos(41, 3), Some(1));
}

// The same for the hash, whose missing record reads 0.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:830-836
#[test]
fn fresh_hash_with_window_minus_40_skips_words_below_40() {
    let table = NegTable::new(-40, None, -90);
    for s_off in 0..40 {
        assert_eq!(table.seed_pos(s_off, -7), None, "s_off {s_off}");
    }
    assert_eq!(table.seed_pos(40, -7), Some(0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-185
#[test]
fn state_carries_from_one_call_to_the_next() {
    let mut table = NegTable::new(-40, Some(64), 0);
    // A hit on diagonal 5 ends at subject offset 100 (position 60 with the offset -40).
    table.store(5, 100 - 40, 60);
    // The call ends: the offset passes on by the subject length plus the window.
    table.exit(1000);
    assert_eq!(table.offset(), -40 + 1000 - 40);
    // The next call's words are at positions s_off + 920: all of them past the stored 60.
    assert_eq!(table.seed_pos(0, 5), Some(920));
    // A second short call: the offset moves by 10 - 40 and goes back below the first one's end.
    table.exit(10);
    assert_eq!(table.offset(), 920 + 10 - 40);
    // A diagonal whose last hit was stored beyond the offset blocks the next call's first words
    // (positions 890 to 899 are before the stored 900).
    table.store(9, 900, 0);
    for s_off in 0..10 {
        assert_eq!(table.seed_pos(s_off, 9), None, "s_off {s_off}");
    }
    assert_eq!(table.seed_pos(10, 9), Some(900));
    assert_eq!(table.seed_pos(10, 5), Some(900));
    // After one more call the offset is lower again and blocks 10 more words of the diagonal.
    table.exit(20);
    assert_eq!(table.offset(), 890 + 20 - 40);
    assert_eq!(table.seed_pos(29, 9), None);
    assert_eq!(table.seed_pos(30, 9), Some(900));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-62
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:663-665
#[test]
fn array_shorter_than_the_query_aliases_diagonals() {
    // A query block of 30 letters and a window of -40: the length stays 1, the mask is 0.
    let mut table = NegTable::new(-40, Some(1), 0);
    table.store(5, 25, 0);
    assert_eq!(table.last_hit(5), 25);
    assert_eq!(table.last_hit(-13), 25);
    assert_eq!(table.last_hit(1000), 25);
    // A table of 4 cells: diagonals 3 and 7 and -1 share a cell, 2 does not.
    let mut table = NegTable::new(-40, Some(4), 0);
    table.store(3, 11, 0);
    assert_eq!(table.last_hit(7), 11);
    assert_eq!(table.last_hit(-1), 11);
    assert_eq!(table.last_hit(2), 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:90-106,162-185
#[test]
fn exit_clears_the_table_once_the_offset_passes_a_quarter_of_int_max() {
    let mut table = NegTable::new(-40, Some(8), 0);
    table.store(2, 77, 0);
    table.offset = i64::from(i32::MAX / 4) - 1;
    table.exit(100);
    assert_eq!(table.offset(), i64::from(i32::MAX / 4) - 1 + 100 - 40);
    assert_eq!(table.last_hit(2), 77);
    table.exit(100);
    assert_eq!(table.offset(), -40);
    // The array's cells hold -window after the clear: words below 80 are blocked.
    assert_eq!(table.last_hit(2), 40);
    assert_eq!(table.seed_pos(79, 2), None);
    assert_eq!(table.seed_pos(80, 2), Some(40));
    let mut hash = NegTable::new(-40, None, -90);
    hash.store(2, 77, 0);
    hash.offset = i64::from(i32::MAX / 4);
    hash.exit(100);
    assert_eq!(hash.offset(), -40);
    assert_eq!(hash.last_hit(2), 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
#[test]
fn hash_replaces_a_stale_cell_of_its_bucket_only() {
    // Diagonals 3 and 515 share a bucket. Word length 11, scan range 0: Delta is -51 and the
    // stale window -40 - 51 + 1 = -90.
    assert_eq!(bucket(3), bucket(515));
    let mut table = NegTable::new(-40, None, -90);
    table.store(3, 100, 90);
    // 20 - 100 = -80 is above -90: the cell is stale and is taken over.
    table.store(515, 200, 20);
    assert_eq!(table.last_hit(515), 200);
    assert_eq!(table.last_hit(3), 0);
    let mut table = NegTable::new(-40, None, -90);
    table.store(3, 100, 90);
    // 5 - 100 = -95 is not above -90: a new cell, both kept.
    table.store(515, 200, 5);
    assert_eq!(table.last_hit(515), 200);
    assert_eq!(table.last_hit(3), 100);
    // The same diagonal is updated in place.
    table.store(3, 300, 5);
    assert_eq!(table.last_hit(3), 300);
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:62-71
#[test]
fn stored_levels_are_31_bit() {
    let mut table = NegTable::new(-40, Some(2), 0);
    table.store(0, 1 << 30, 0);
    assert_eq!(table.last_hit(0), -(1 << 30));
    let mut hash = NegTable::new(-40, None, -90);
    hash.store(4, (1 << 30) + 3, 0);
    assert_eq!(hash.last_hit(4), -(1 << 30) + 3);
}

// A window at the end of the `Int4` range never panics and keeps NCBI's values where they are
// defined (the widened offset is not wrapped).
#[test]
fn extreme_windows_do_not_overflow() {
    let mut table = NegTable::new(i32::MIN, Some(1), 0);
    assert_eq!(table.offset(), i64::from(i32::MIN));
    assert_eq!(table.seed_pos(5, 0), None);
    table.store(0, i64::from(i32::MIN) + 5, 0);
    table.exit(i32::MAX as usize);
    assert!(table.offset() < i64::from(i32::MIN));
    let mut hash = NegTable::new(i32::MIN, None, i64::from(i32::MIN) * 2);
    hash.store(i64::from(i32::MIN), 7, i64::from(i32::MIN));
    assert_eq!(hash.last_hit(i64::from(i32::MIN)), 7);
}
