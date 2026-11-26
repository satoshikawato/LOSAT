use super::*;

fn hsp(query: (i32, i32), score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: query.0,
        query_end: query.1,
        subject_start: query.0,
        subject_end: query.1,
        gapped_query_start: query.0,
        gapped_subject_start: query.0,
        raw_score: score,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-999 (cited at
// `delete_failed_hsp_tests`): `Blast_HSPTest` of a BLASTP HSP fails for a negative length only.
#[test]
fn hsp_test_deletes_a_negative_length_and_keeps_zero_in_order() {
    let mut hsps = vec![
        hsp((10, 20), 30),
        hsp((12, 11), 25),
        hsp((5, 5), 20),
        hsp((40, 39), 15),
        hsp((7, 9), 10),
    ];
    delete_failed_hsp_tests(&mut hsps);
    let kept: Vec<(i32, i32, i32)> = hsps
        .iter()
        .map(|hsp| (hsp.query_start, hsp.query_end, hsp.raw_score))
        .collect();
    assert_eq!(kept, [(10, 20, 30), (5, 5, 20), (7, 9, 10)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:776-782
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:36-66
// A table made for a negative window starts at offset `window` with every entry 0, so an entry
// reads `-window` letters ahead: a word of the first `-window` letters of a subject is not
// extended on a diagonal that nothing has written (the offset moves on after each subject).
#[test]
fn a_fresh_table_with_a_negative_window_skips_the_words_of_the_first_letters() {
    for window in [-1, -40] {
        let offset = i64::from(window);
        assert_eq!(
            neg_window_diff(-window - 1, 0, offset),
            -1,
            "window {window}"
        );
        assert_eq!(neg_window_diff(-window, 0, offset), 0, "window {window}");
        assert!(neg_window_diff(0, 0, offset) < 0, "window {window}");
    }
    // A window of 0 is the command line's one-hit search: nothing is skipped.
    assert_eq!(neg_window_diff(0, 0, 0), 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-172
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:799-800
// The table's state passes from one subject to the next: after a subject of 300 letters the offset
// is 220 (from -40, less the window twice), an extension that ended at letter 290 left
// `290 - 2 + (-40)` in its entry, and the next subject's words on that diagonal before letter 28
// are not extended, as are none of the entries a subject has not written.
#[test]
fn a_negative_window_carries_the_table_from_a_subject_to_the_next() {
    let window = -40;
    let (offset, cleared) = neg_window_exit(i64::from(window), window, 300);
    assert_eq!((offset, cleared), (220, false));
    let stored = neg_window_last_hit(290, 3, i64::from(window));
    assert_eq!(stored, 248);
    assert_eq!(neg_window_diff(27, stored, offset), -1);
    assert_eq!(neg_window_diff(28, stored, offset), 0);
    assert!(neg_window_diff(0, 0, offset) > 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-106,162-172
// From an offset of a quarter of `INT4_MAX` the table is cleared to minus the window and the offset goes
// back to the window.
#[test]
fn a_negative_window_clears_the_table_at_a_quarter_of_the_int4_range() {
    let limit = i64::from(i32::MAX / 4);
    assert_eq!(neg_window_exit(limit, -40, 100), (-40, true));
    assert_eq!(
        neg_window_exit(limit - 1, -40, 100),
        (limit - 1 + 100 - 40, false)
    );
    let entry = DiagStruct::clear(-40);
    assert_eq!((entry.last_hit(), entry.flag()), (40, 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:776-782,799-800
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
// The stored value is the sum's low 31 bits, a difference that NCBI's `Uint4` subtraction wraps
// wraps (where the `Int4` subtraction before it is defined), and an `Int4` that would overflow is
// the 64-bit value: a window of `INT_MIN` runs without a panic.
#[test]
fn a_negative_window_stores_31_bits_wraps_where_ncbi_does_and_widens_where_it_overflows() {
    let mut entry = DiagStruct::default();
    entry.set_last_hit(neg_window_last_hit(10, 3, 1 << 30));
    assert_eq!(i64::from(entry.last_hit()), (1i64 << 30) + 8 - (1i64 << 31));
    // `last_hit - offset` is `INT_MIN`, defined; the subject offset less it wraps in 32 bits.
    assert_eq!(neg_window_diff(i32::MAX, -(1 << 30), 1 << 30), -1);
    // `last_hit - offset` would overflow `Int4`: the 64-bit difference.
    let least = i64::from(i32::MIN);
    assert_eq!(neg_window_diff(5, 0, least), 5 - (1i64 << 31));
    assert_eq!(
        neg_window_exit(least, i32::MIN, 300),
        (2 * least + 300, false)
    );
    let entry = DiagStruct::clear(i32::MIN);
    assert_eq!(entry.last_hit(), 0);
    assert_eq!(neg_window_last_hit(7, 3, least), (5 + least) as i32);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-3035 (`Blast_HSPListsMerge`)
#[test]
fn merge_prelim_subject_chunks_joins_overlapping_hsp_across_split() {
    let mut first = hsp((10, 80), 120);
    first.subject_start = 4_999_880;
    first.subject_end = 4_999_950;
    first.gapped_subject_start = 4_999_880;
    let mut second = hsp((40, 120), 140);
    second.subject_start = 4_999_910;
    second.subject_end = 4_999_990;
    second.gapped_subject_start = 4_999_910;
    let mut combined = vec![first];
    let mut incoming = vec![second];
    merge_prelim_subject_chunks(&mut combined, &mut incoming, 4_999_900, 100, true);
    assert_eq!(combined.len(), 1);
    assert_eq!(
        (
            combined[0].query_start,
            combined[0].query_end,
            combined[0].subject_start,
            combined[0].subject_end
        ),
        (10, 120, 4_999_880, 4_999_990)
    );
}
