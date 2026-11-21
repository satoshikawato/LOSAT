use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213 (cited at
// `scan_ranges`): the left end is the last letter of the run before, included.
#[test]
fn scan_ranges_run_from_the_last_lower_case_letter() {
    assert!(scan_ranges(b"MKV").is_empty());
    assert_eq!(scan_ranges(b"MKvvA"), [(0, 2), (3, 5)]);
    assert_eq!(scan_ranges(b"mKVvA"), [(0, 0), (0, 3), (3, 5)]);
    assert_eq!(scan_ranges(b"MKv"), [(0, 2), (2, 3)]);
    assert_eq!(scan_ranges(b"mkv"), [(0, 0), (2, 3)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-310 (cited at
// `chunk_scans`): an 11,200,000-letter subject is cut at 4,999,900 and 9,999,800. A chunk
// starts with the last range of the chunk before; when that range reaches into it by fewer
// letters than a word, its one-hit scan stops at once (aa_ungapped.c:762-766) and its two-hit
// scan takes one word (`:499-500`).
#[test]
fn each_chunk_starts_its_own_scan() {
    let length = 11_200_000;
    // A run of 4,999,890..=5,000,050: the second chunk's first range starts 150 letters in,
    // and the scan is that of the whole subject.
    let wide = [(0, 4_999_890), (5_000_050, length)];
    for raise in [false, true] {
        assert_eq!(
            chunk_scans(&wide, length, 3, raise),
            (vec![(0, 4_999_887), (5_000_050, length - 3)], vec![])
        );
    }
    // A run of 4,999,902..=5,000,001: the second chunk's first range is its letters 0..2.
    let short = [(0, 4_999_902), (5_000_001, length)];
    assert_eq!(
        chunk_scans(&short, length, 3, false),
        (vec![(0, 4_999_899), (9_999_800, length - 3)], vec![])
    );
    assert_eq!(
        chunk_scans(&short, length, 3, true),
        (
            vec![(0, 4_999_900), (5_000_001, length - 3)],
            vec![4_999_900]
        )
    );
    // A run of 4,999,800..=9,999,999 covers the second chunk (SUBJECT_SPLIT_NO_RANGE).
    let covered = [(0, 4_999_800), (9_999_999, length)];
    for raise in [false, true] {
        assert_eq!(
            chunk_scans(&covered, length, 3, raise),
            (vec![(0, 4_999_797), (9_999_999, length - 3)], vec![])
        );
    }
    // 5,000,001 letters: a second chunk of 101 letters, whose first range is 0..2.
    let tail = [(0, 4_999_902), (5_000_000, 5_000_001)];
    assert_eq!(
        chunk_scans(&tail, 5_000_001, 3, false),
        (vec![(0, 4_999_899)], vec![])
    );
    assert_eq!(
        chunk_scans(&tail, 5_000_001, 3, true),
        (vec![(0, 4_999_900)], vec![4_999_900])
    );
}
