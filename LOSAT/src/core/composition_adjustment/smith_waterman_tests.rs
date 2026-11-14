//! Tests of the Smith-Waterman kernels on hand-made rows (the brute-force reference below fills
//! NCBI's three-score recurrence on whole matrices instead of one rolling row).

use super::*;

/// A PSSM row source: `rows[query_pos]`.
fn pssm<'a>(rows: &'a [[i32; BLASTAA_SIZE]]) -> impl Fn(usize) -> &'a [i32; BLASTAA_SIZE] + 'a {
    move |pos| &rows[pos]
}

/// A row that scores `hit` for the letters of `letters` and `miss` for every other.
fn row(letters: &[u8], hit: i32, miss: i32) -> [i32; BLASTAA_SIZE] {
    let mut row = [miss; BLASTAA_SIZE];
    for &letter in letters {
        row[letter as usize] = hit;
    }
    row
}

/// The rows of an identity scoring of `query` (a square matrix read by the query's letters, NCBI's
/// `matrix[query[queryPos]]`).
fn identity_matrix(hit: i32, miss: i32) -> Vec<[i32; BLASTAA_SIZE]> {
    (0..BLASTAA_SIZE)
        .map(|letter| row(&[letter as u8], hit, miss))
        .collect()
}

/// NCBI's recurrence on whole matrices: `h` the cell's score, `e` the gap in the subject (along
/// the row), `f` the gap in the query (down the column), the forbidden cells' diagonal at
/// `COMPO_SCORE_MIN`; the best cell is the first strict maximum in the scan order.
fn reference_score_only(
    subject: &[u8],
    query_rows: &[[i32; BLASTAA_SIZE]],
    gap_open: i32,
    gap_extend: i32,
    forbidden: &dyn Fn(usize, usize) -> bool,
) -> (i32, usize, usize) {
    let m = query_rows.len();
    let n = subject.len();
    let cost = gap_open + gap_extend;
    let mut h = vec![vec![0i32; n + 1]; m + 1];
    let mut f_prev = vec![-gap_open; n + 1];
    let mut best = (0, 0, 0);
    for i in 1..=m {
        let mut e = -gap_open;
        let mut f_row = vec![-gap_open; n + 1];
        for j in 1..=n {
            e = (h[i][j - 1] - cost).max(e - gap_extend);
            let f = (h[i - 1][j] - cost).max(f_prev[j] - gap_extend);
            let diagonal = if forbidden(i - 1, j - 1) {
                COMPO_SCORE_MIN
            } else {
                h[i - 1][j - 1] + query_rows[i - 1][subject[j - 1] as usize]
            };
            h[i][j] = diagonal.max(0).max(e).max(f);
            f_row[j] = f;
            if h[i][j] > best.0 {
                best = (h[i][j], j - 1, i - 1);
            }
        }
        f_prev = f_row;
    }
    best
}

fn score_only(
    subject: &[u8],
    query_rows: &[[i32; BLASTAA_SIZE]],
    forbidden: &ForbiddenRanges,
) -> SwEnd {
    smith_waterman_score_only(
        subject,
        query_rows.len(),
        pssm(query_rows),
        11,
        1,
        forbidden,
        &mut SwScratch::new(),
    )
}

fn find_start(
    subject: &[u8],
    query_rows: &[[i32; BLASTAA_SIZE]],
    end: SwEnd,
    forbidden: &ForbiddenRanges,
) -> SwStart {
    smith_waterman_find_start(
        subject,
        pssm(query_rows),
        11,
        1,
        end.match_end,
        end.query_end,
        end.score,
        forbidden,
        &mut SwScratch::new(),
    )
}

// A position-specific row source: the profile A C D (one row each) on the subject W A C D W:
// score 3 x 5 ending at (2, 3), starting at (0, 1).
#[test]
fn a_profile_aligns_where_its_rows_score() {
    let rows = [row(&[1], 5, -3), row(&[3], 5, -3), row(&[4], 5, -3)];
    let subject = [17, 1, 3, 4, 17];
    let none = ForbiddenRanges::new();
    let end = score_only(&subject, &rows, &none);
    assert_eq!(
        end,
        SwEnd {
            score: 15,
            match_end: 3,
            query_end: 2
        }
    );
    assert_eq!(
        find_start(&subject, &rows, end, &none),
        SwStart {
            score: 15,
            match_start: 1,
            query_start: 0
        }
    );
}

// The square matrix through the query's letters (NCBI's `matrix[query[queryPos]]`) is a row
// source of the same kernel.
#[test]
fn a_square_matrix_is_a_row_source_too() {
    let matrix = identity_matrix(4, -2);
    let query = [1u8, 3, 4, 5];
    let subject = [1u8, 3, 4, 5];
    let end = smith_waterman_score_only(
        &subject,
        query.len(),
        |pos| &matrix[query[pos] as usize],
        11,
        1,
        &ForbiddenRanges::new(),
        &mut SwScratch::new(),
    );
    assert_eq!((end.score, end.match_end, end.query_end), (16, 3, 3));
}

// Two cells of the best score: the score-only scan keeps the first (the lower query position,
// then the lower subject position); the find-start scan from the second best end runs backwards.
#[test]
fn ties_keep_the_first_cell_of_each_scan() {
    let rows = [row(&[1], 6, -9), row(&[1], 6, -9)];
    // Query rows 0 and 1 both score A; the subject A x A: the best 6 first at (0, 0).
    let subject = [1, 20, 1];
    let none = ForbiddenRanges::new();
    let end = score_only(&subject, &rows, &none);
    assert_eq!((end.score, end.query_end, end.match_end), (6, 0, 0));
    // From the end (1, 2) with 6: the reverse scan reaches 6 at once.
    let start = find_start(
        &subject,
        &rows,
        SwEnd {
            score: 6,
            match_end: 2,
            query_end: 1,
        },
        &none,
    );
    assert_eq!(
        (start.score, start.query_start, start.match_start),
        (6, 1, 2)
    );
}

// A forbidden range over the best cell moves the alignment to the next best; a range pushed for
// query positions [start, end) holds subject positions [match_start, match_end] (the end
// included).
#[test]
fn forbidden_ranges_cover_the_best_cell() {
    let rows = [row(&[1], 5, -3), row(&[3], 5, -3), row(&[4], 5, -3)];
    // A C D twice; the first copy scores 15 and ends first.
    let subject = [1, 3, 4, 17, 17, 1, 3, 4];
    let mut forbidden = ForbiddenRanges::new();
    forbidden.initialize(rows.len());
    let first = score_only(&subject, &rows, &forbidden);
    assert_eq!((first.score, first.query_end, first.match_end), (15, 2, 2));
    forbidden.push(0, 3, 0, 2);
    let second = score_only(&subject, &rows, &forbidden);
    assert_eq!(
        (second.score, second.query_end, second.match_end),
        (15, 2, 7)
    );
    let start = find_start(&subject, &rows, second, &forbidden);
    assert_eq!((start.query_start, start.match_start), (0, 5));
    // The range's end is included: [0, 4] leaves the second copy only from subject position 5.
    forbidden.clear();
    forbidden.push(0, 3, 0, 4);
    assert_eq!(score_only(&subject, &rows, &forbidden).match_end, 7);
    // The query end is not included: position 2 keeps no range, so D at subject position 2 still
    // scores 5 alone (and the second copy wins).
    forbidden.clear();
    forbidden.push(0, 2, 0, 7);
    let left = score_only(&subject, &rows, &forbidden);
    assert_eq!((left.score, left.query_end, left.match_end), (5, 2, 2));
    // Cleared ranges forbid nothing.
    forbidden.clear();
    assert_eq!(score_only(&subject, &rows, &forbidden), first);
}

// Ranges at the ends of the subject: position 0 and the last position.
#[test]
fn ranges_at_the_ends() {
    let rows = [row(&[1], 7, -4)];
    let subject = [1, 2, 2, 1];
    let mut forbidden = ForbiddenRanges::new();
    forbidden.initialize(1);
    forbidden.push(0, 1, 0, 0);
    let end = score_only(&subject, &rows, &forbidden);
    assert_eq!((end.score, end.match_end), (7, 3));
    forbidden.push(0, 1, 3, 3);
    let end = score_only(&subject, &rows, &forbidden);
    assert_eq!((end.score, end.query_end, end.match_end), (0, 0, 0));
}

// Empty windows: no rows or no subject letters give score 0 at (0, 0).
#[test]
fn empty_windows_score_zero() {
    let none = ForbiddenRanges::new();
    let rows = [row(&[1], 5, -3)];
    assert_eq!(
        score_only(&[], &rows, &none),
        SwEnd {
            score: 0,
            match_end: 0,
            query_end: 0
        }
    );
    assert_eq!(
        score_only(&[1, 2], &[], &none),
        SwEnd {
            score: 0,
            match_end: 0,
            query_end: 0
        }
    );
    // A score of 0 to reach (the score-only scan found no positive cell): the scan stops after its
    // first cell, whose score 0 leaves the start at (0, 0).
    let start = find_start(
        &[1, 2],
        &rows,
        SwEnd {
            score: 0,
            match_end: 1,
            query_end: 0,
        },
        &none,
    );
    assert_eq!(
        start,
        SwStart {
            score: 0,
            match_start: 0,
            query_start: 0
        }
    );
}

// The reverse scan stops where the best score first reaches the score to find, even inside a row,
// and gaps count as NCBI's (open 11 + extend 1 per letter).
#[test]
fn the_find_start_scan_stops_at_the_score() {
    // Query rows: A A A ... (each 9 for A), subject A A A A: the best 36 on the diagonal.
    let rows = [row(&[1], 9, -9); 4];
    let subject = [1, 1, 1, 1];
    let none = ForbiddenRanges::new();
    let end = score_only(&subject, &rows, &none);
    assert_eq!((end.score, end.query_end, end.match_end), (36, 3, 3));
    let start = find_start(&subject, &rows, end, &none);
    assert_eq!(
        (start.score, start.query_start, start.match_start),
        (36, 0, 0)
    );
    // With 18 to reach, the scan stops at the second cell of the diagonal.
    let start = find_start(&subject, &rows, SwEnd { score: 18, ..end }, &none);
    assert_eq!(
        (start.score, start.query_start, start.match_start),
        (18, 2, 2)
    );
}

/// A small xorshift for the random rows of the reference test.
struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

// The rolling row of the kernel gives the whole-matrix recurrence's score and end, with and
// without forbidden ranges, on random rows with gaps cheap enough to be taken.
#[test]
fn the_score_only_scan_equals_the_whole_matrix_recurrence() {
    let mut random = XorShift(0x5eed_1234_abcd_0001);
    for case in 0..400 {
        let m = random.below(12) as usize;
        let n = random.below(14) as usize;
        let rows: Vec<[i32; BLASTAA_SIZE]> = (0..m)
            .map(|_| {
                let mut row = [0; BLASTAA_SIZE];
                for score in row.iter_mut() {
                    *score = random.below(15) as i32 - 9;
                }
                row
            })
            .collect();
        let subject: Vec<u8> = (0..n).map(|_| random.below(6) as u8 + 1).collect();
        let (gap_open, gap_extend) = (random.below(4) as i32, random.below(3) as i32 + 1);
        let mut forbidden = ForbiddenRanges::new();
        forbidden.initialize(m);
        for _ in 0..random.below(3) {
            if m == 0 || n == 0 {
                break;
            }
            let q0 = random.below(m as u64) as i32;
            let q1 = q0 + random.below((m as i32 - q0) as u64 + 1) as i32;
            let s0 = random.below(n as u64) as i32;
            let s1 = s0 + random.below(3) as i32;
            forbidden.push(q0, q1, s0, s1);
        }
        let is_forbidden = |q: usize, s: usize| {
            forbidden
                .at(q)
                .iter()
                .any(|&(a, b)| s as i32 >= a && s as i32 <= b)
        };
        let expected = reference_score_only(&subject, &rows, gap_open, gap_extend, &is_forbidden);
        let end = smith_waterman_score_only(
            &subject,
            m,
            pssm(&rows),
            gap_open,
            gap_extend,
            &forbidden,
            &mut SwScratch::new(),
        );
        assert_eq!(
            (end.score, end.match_end, end.query_end),
            expected,
            "case {case}"
        );
        // The find-start scan returns the score it was asked for, at a cell no later than the end.
        if end.score > 0 {
            let start = smith_waterman_find_start(
                &subject,
                pssm(&rows),
                gap_open,
                gap_extend,
                end.match_end,
                end.query_end,
                end.score,
                &forbidden,
                &mut SwScratch::new(),
            );
            assert_eq!(start.score, end.score, "case {case}");
            assert!(start.query_start <= end.query_end && start.match_start <= end.match_end);
            assert!(
                !is_forbidden(start.query_start, start.match_start),
                "case {case}"
            );
        }
    }
}
