use super::{RpsPastRows, RpsRows};
use crate::algorithm::pssm_search::rows::PssmRows;
use crate::algorithm::pssm_search::ungapped::pssm_extend_one_hit;
use crate::core::blast_util::BLASTAA_SIZE;
use std::sync::atomic::{AtomicBool, Ordering};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:938-975, 1036-1083 (quoted
// at `pssm_search::ungapped`): the extension to the right reads up to `query->length`, the
// concatenated length `Lc`, which `RpsRows::length` holds within the rows, and past them the
// row after the rows (`RpsRows::past`). Nine rows that score letter 1 (A) 5,000 and one
// sentinel row of -32,767 (`-BLAST_SCORE_MAX`); a word of three letters at 0 against a query
// of A.
#[test]
fn the_extension_to_the_right_stops_at_the_concatenated_length() {
    let mut row = [-4i32; BLASTAA_SIZE];
    row[1] = 5_000;
    let mut rows = vec![row; 9];
    rows.push([-32_767; BLASTAA_SIZE]);
    let query = [1u8; 20];
    let zeros = [Some(0i32); BLASTAA_SIZE];
    let faults = [None; BLASTAA_SIZE];
    let extend = |length: usize,
                  past: Option<(&[Option<i32>; BLASTAA_SIZE], i32)>,
                  subject: &[u8],
                  dropoff: i32| {
        let past = past.map(|(row, concat_length)| RpsPastRows {
            row,
            concat_length,
            overran: AtomicBool::new(false),
        });
        let rows = RpsRows {
            rows: &rows,
            length,
            past: past.as_ref(),
        };
        let hit = pssm_extend_one_hit(&rows, subject, 0, 0, dropoff, 3);
        let flagged = past.is_some_and(|past| past.overran.load(Ordering::Relaxed));
        (hit, flagged)
    };
    // `Lc` 6 below the 10 rows: the walk from 3 stops at 6.
    let (short, _) = extend(6, None, &query, 40_000);
    assert_eq!(
        (short.ungapped_data.length, short.ungapped_data.score),
        (6, 30_000)
    );
    assert_eq!(short.s_last_off, 6);
    // `Lc` 11: the sentinel row stops an X-drop of 20.
    let (stopped, flagged) = extend(10, Some((&zeros, 11)), &query, 20);
    assert_eq!(
        (stopped.ungapped_data.length, stopped.ungapped_data.score),
        (9, 45_000)
    );
    assert_eq!(stopped.s_last_off, 9);
    assert!(!flagged);
    // An X-drop above the sentinel's scores: the walk scores the sentinel row, then NCBI's row
    // after the rows (zeros of the `.rps` file's last page) and ends at `Lc` 11.
    let (through, flagged) = extend(10, Some((&zeros, 11)), &query, 40_000);
    assert_eq!(through.s_last_off, 11);
    assert_eq!(through.ungapped_data.length, 9);
    assert!(!flagged);
    // `Lc` 10 (the rows): NCBI's bound is the rows' end.
    let (at_rows, _) = extend(10, None, &query, 40_000);
    assert_eq!(at_rows.s_last_off, 10);
    // The row after the rows past the file's last page: NCBI's read faults.
    assert!(extend(10, Some((&faults, 11)), &query, 40_000).1);
    // `Lc` 12: NCBI reads past the pointer array after that row.
    assert!(extend(10, Some((&zeros, 12)), &query, 40_000).1);
    // The query ends with the rows: NCBI's bound is the query's end too.
    assert!(!extend(10, Some((&faults, 12)), &query[..10], 40_000).1);
    // A row after the rows of other bytes (a `.rps` file with bytes after its rows) scores
    // them: 40,000 for A, above the sentinel's -32,767, takes the best score to that row.
    let mut high = zeros;
    high[1] = Some(40_000);
    let (scored, flagged) = extend(10, Some((&high, 11)), &query, 1_000_000);
    assert_eq!(scored.ungapped_data.length, 11);
    assert!(!flagged);
}

// `RpsRows::whole` gives the rows as one slice only where the slice's own length is every
// bound of the scan: `Lc` is the rows' count and there is no row after the rows. An `Lc`
// below the rows (an alias file's `LENGTH`) and an `Lc` past them keep `RpsRows`.
#[test]
fn the_rows_are_one_slice_only_where_the_concatenated_length_is_their_count() {
    let rows = vec![[0i32; BLASTAA_SIZE]; 10];
    let zeros = [Some(0i32); BLASTAA_SIZE];
    let past = RpsPastRows {
        row: &zeros,
        concat_length: 11,
        overran: AtomicBool::new(false),
    };
    let whole = |length: usize, past: Option<&RpsPastRows<'_>>| {
        let rows = RpsRows {
            rows: &rows,
            length,
            past,
        };
        rows.whole().map(|whole| PssmRows::len(&whole))
    };
    assert_eq!(whole(10, None), Some(10));
    assert_eq!(whole(9, None), None);
    assert_eq!(whole(0, None), None);
    assert_eq!(whole(10, Some(&past)), None);
}
