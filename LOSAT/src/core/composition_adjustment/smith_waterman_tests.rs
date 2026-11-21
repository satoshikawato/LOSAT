//! Tests of the Smith-Waterman kernels on hand-made rows (the brute-force reference below fills
//! NCBI's three-score recurrence on whole matrices instead of one rolling row).

use super::super::compo_scale::COMPO_SCORE_MIN;
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
// without forbidden ranges, on random rows with gaps cheap enough to be taken; one scratch serves
// every case (the rows' marks of their ranges are cleared after each row).
#[test]
fn the_score_only_scan_equals_the_whole_matrix_recurrence() {
    let mut random = XorShift(0x5eed_1234_abcd_0001);
    let mut scratch = SwScratch::new();
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
            &mut scratch,
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
                &mut scratch,
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

/// The scans as they were before the gap along the row was carried alone (dev e6e9a61ad): each
/// cell takes the gap along the row from the score of the cell before it, in NCBI's order of the
/// cell's terms. The reference of the kernel's cell (`cell`) in the tests below; its maxima are
/// `i32::max`, not the kernel's `cmov_max`, so that the reference shares no operation of the cell.
mod before_k1 {
    use super::*;

    fn padded_row(row: &[i32; BLASTAA_SIZE]) -> [i32; ROW_SLOTS_BEFORE] {
        let mut padded = [0; ROW_SLOTS_BEFORE];
        padded[..BLASTAA_SIZE].copy_from_slice(row);
        padded
    }

    const ROW_SLOTS_BEFORE: usize = 32;

    /// The scans' memory as it was: the subject's letters, the columns and the marks of the
    /// forbidden ranges.
    fn scan(subject: &[u8], gap_open: i32) -> (Vec<u8>, Vec<SwGapInfo>, Vec<bool>) {
        let letters = subject.iter().map(|&l| row_column(l) as u8).collect();
        let column = SwGapInfo {
            no_gap: 0,
            gap_exists: -gap_open,
        };
        (
            letters,
            vec![column; subject.len()],
            vec![false; subject.len()],
        )
    }

    fn blocked_row(blocked: &mut [bool], ranges: &[(i32, i32)], value: bool) {
        let last = blocked.len() as i64 - 1;
        for &(start, end) in ranges {
            let (start, end) = (i64::from(start).max(0), i64::from(end).min(last));
            if start <= end {
                blocked[start as usize..=end as usize].fill(value);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn cell(
        new_score: &mut i32,
        prev_score_gap_match_seq: &mut i32,
        prev_score_no_gap_match_seq: &mut i32,
        column: &mut SwGapInfo,
        substitution: Option<i32>,
        new_gap_cost: i32,
        gap_extend: i32,
    ) {
        *prev_score_gap_match_seq =
            (*prev_score_gap_match_seq - gap_extend).max(*new_score - new_gap_cost);
        let continue_gap_score = (column.gap_exists - gap_extend).max(column.no_gap - new_gap_cost);
        let diagonal = match substitution {
            None => COMPO_SCORE_MIN,
            Some(substitution) => *prev_score_no_gap_match_seq + substitution,
        };
        let score = diagonal
            .max(0)
            .max(continue_gap_score)
            .max(*prev_score_gap_match_seq);
        *prev_score_no_gap_match_seq = column.no_gap;
        column.no_gap = score;
        column.gap_exists = continue_gap_score;
        *new_score = score;
    }

    fn substitution(row: &[i32; ROW_SLOTS_BEFORE], blocked: bool, letter: u8) -> Option<i32> {
        if blocked {
            None
        } else {
            Some(row[usize::from(letter) & (ROW_SLOTS_BEFORE - 1)])
        }
    }

    fn score_only_row(
        subject: &[u8],
        row: &[i32; ROW_SLOTS_BEFORE],
        blocked: &[bool],
        score_vector: &mut [SwGapInfo],
        (gap_open, gap_extend): (i32, i32),
        mut best: i32,
    ) -> Option<(i32, usize)> {
        let new_gap_cost = gap_open + gap_extend;
        let mut found = None;
        let mut new_score = 0;
        let mut prev_score_no_gap_match_seq = 0;
        let mut prev_score_gap_match_seq = -gap_open;
        for (match_seq_pos, ((&letter, column), &blocked)) in subject
            .iter()
            .zip(score_vector.iter_mut())
            .zip(blocked)
            .enumerate()
        {
            cell(
                &mut new_score,
                &mut prev_score_gap_match_seq,
                &mut prev_score_no_gap_match_seq,
                column,
                substitution(row, blocked, letter),
                new_gap_cost,
                gap_extend,
            );
            if new_score > best {
                best = new_score;
                found = Some((best, match_seq_pos));
            }
        }
        found
    }

    #[allow(clippy::too_many_arguments)]
    fn find_start_row(
        subject: &[u8],
        row: &[i32; ROW_SLOTS_BEFORE],
        blocked: &[bool],
        score_vector: &mut [SwGapInfo],
        match_end: usize,
        (gap_open, gap_extend): (i32, i32),
        mut best: i32,
        score_in: i32,
    ) -> (Option<(i32, usize)>, bool) {
        let new_gap_cost = gap_open + gap_extend;
        let mut found = None;
        let mut new_score = 0;
        let mut prev_score_no_gap_match_seq = 0;
        let mut prev_score_gap_match_seq = -gap_open;
        for (match_seq_pos, ((&letter, column), &blocked)) in subject[..=match_end]
            .iter()
            .zip(score_vector[..=match_end].iter_mut())
            .zip(&blocked[..=match_end])
            .enumerate()
            .rev()
        {
            cell(
                &mut new_score,
                &mut prev_score_gap_match_seq,
                &mut prev_score_no_gap_match_seq,
                column,
                substitution(row, blocked, letter),
                new_gap_cost,
                gap_extend,
            );
            if new_score > best {
                best = new_score;
                found = Some((best, match_seq_pos));
            }
            if best >= score_in {
                return (found, true);
            }
        }
        (found, false)
    }

    pub(super) fn score_only<'r, F>(
        subject: &[u8],
        query_length: usize,
        rows: F,
        gap_open: i32,
        gap_extend: i32,
        forbidden: &ForbiddenRanges,
    ) -> SwEnd
    where
        F: Fn(usize) -> &'r [i32; BLASTAA_SIZE],
    {
        let (subject, mut score_vector, mut blocked) = scan(subject, gap_open);
        let (mut best_score, mut best_query_pos, mut best_match_seq_pos) = (0, 0, 0);
        for query_pos in 0..query_length {
            let row = &padded_row(rows(query_pos));
            let ranges = forbidden.at(query_pos);
            blocked_row(&mut blocked, ranges, true);
            let found = score_only_row(
                &subject,
                row,
                &blocked,
                &mut score_vector,
                (gap_open, gap_extend),
                best_score,
            );
            blocked_row(&mut blocked, ranges, false);
            if let Some((score, match_seq_pos)) = found {
                best_score = score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
        }
        SwEnd {
            score: best_score.max(0),
            match_end: best_match_seq_pos,
            query_end: best_query_pos,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn find_start<'r, F>(
        subject: &[u8],
        rows: F,
        gap_open: i32,
        gap_extend: i32,
        match_end: usize,
        query_end: usize,
        score_in: i32,
        forbidden: &ForbiddenRanges,
    ) -> SwStart
    where
        F: Fn(usize) -> &'r [i32; BLASTAA_SIZE],
    {
        if subject.is_empty() {
            return SwStart {
                score: 0,
                match_start: 0,
                query_start: 0,
            };
        }
        let (subject, mut score_vector, mut blocked) = scan(subject, gap_open);
        let (mut best_score, mut best_query_pos, mut best_match_seq_pos) = (0, 0, 0);
        for query_pos in (0..=query_end).rev() {
            let row = &padded_row(rows(query_pos));
            let ranges = forbidden.at(query_pos);
            blocked_row(&mut blocked, ranges, true);
            let (found, reached) = find_start_row(
                &subject,
                row,
                &blocked,
                &mut score_vector,
                match_end,
                (gap_open, gap_extend),
                best_score,
                score_in,
            );
            blocked_row(&mut blocked, ranges, false);
            if let Some((score, match_seq_pos)) = found {
                best_score = score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if reached {
                break;
            }
        }
        SwStart {
            score: best_score.max(0),
            match_start: best_match_seq_pos,
            query_start: best_query_pos,
        }
    }
}

/// Gap costs `(open, extend)` of the comparison with `before_k1`: 0 and negative openings (the
/// identity of `cell` holds for every integer opening), NCBI's defaults, the redo's costs scaled by
/// 32 (BLOSUM62's 11/1), and the largest pair that the matrix tables accept (32767/32767, the
/// ungapped row) scaled by 32 and by 100.
const GAP_COSTS: [(i32, i32); 10] = [
    (0, 1),
    (0, 7),
    (-1, 1),
    (-3, 4),
    (11, 1),
    (9, 2),
    (352, 32),
    (32767 * 32, 32767 * 32),
    (32767 * 100, 32767 * 100),
    (1, 3),
];

/// `count` rows of a square matrix read through random query letters (NCBI's
/// `matrix[query[queryPos]]`) or of a random PSSM (`matrix[queryPos]`), each score `scale` times a
/// value in `low..=high`.
fn random_rows(
    random: &mut XorShift,
    count: usize,
    square: bool,
    (low, high, scale): (i32, i32, i32),
) -> Vec<[i32; BLASTAA_SIZE]> {
    let score =
        |random: &mut XorShift| (low + random.below((high - low + 1) as u64) as i32) * scale;
    if square {
        let matrix: Vec<[i32; BLASTAA_SIZE]> = (0..BLASTAA_SIZE)
            .map(|_| std::array::from_fn(|_| score(random)))
            .collect();
        (0..count)
            .map(|_| matrix[random.below(BLASTAA_SIZE as u64) as usize])
            .collect()
    } else {
        (0..count)
            .map(|_| std::array::from_fn(|_| score(random)))
            .collect()
    }
}

/// Random forbidden ranges as the round loop pushes them (query positions `[q0, q1)`, subject
/// positions `[s0, s1]`), some past the subject's ends.
fn random_ranges(random: &mut XorShift, m: usize, n: usize) -> ForbiddenRanges {
    let mut forbidden = ForbiddenRanges::new();
    forbidden.initialize(m);
    if m > 0 && n > 0 {
        for _ in 0..random.below(4) {
            let q0 = random.below(m as u64) as i32;
            let q1 = q0 + random.below((m as i32 - q0) as u64 + 1) as i32;
            let s0 = random.below(n as u64 + 2) as i32 - 1;
            let s1 = s0 + random.below(n as u64 / 2 + 2) as i32;
            forbidden.push(q0, q1, s0, s1);
        }
    }
    forbidden
}

/// Compares both scans with `before_k1` on one window: the score-only scan's
/// `(score, query_end, match_end)`, then the find-start scan from that end for its score, for half
/// of it and for one more (a score not reached: the scan runs to its last row). Each scan runs as
/// the search runs it (on a CPU with AVX2, compiled for it with its cells in `Lane0`) and with its
/// cells in `i32`, the scans of every other CPU.
fn assert_same_as_before_k1(
    subject: &[u8],
    rows: &[[i32; BLASTAA_SIZE]],
    (gap_open, gap_extend): (i32, i32),
    forbidden: &ForbiddenRanges,
    scratch: &mut SwScratch,
    case: &str,
) {
    let m = rows.len();
    let before = before_k1::score_only(subject, m, pssm(rows), gap_open, gap_extend, forbidden);
    let end = smith_waterman_score_only(
        subject,
        m,
        pssm(rows),
        gap_open,
        gap_extend,
        forbidden,
        scratch,
    );
    assert_eq!(end, before, "score-only scan, {case}");
    let plain = score_only_scan::<i32, _>(
        subject,
        m,
        pssm(rows),
        gap_open,
        gap_extend,
        forbidden,
        scratch,
    );
    assert_eq!(plain, before, "score-only scan in i32, {case}");
    if m == 0 {
        return;
    }
    for score_in in [end.score, end.score / 2, end.score + 1] {
        let before = before_k1::find_start(
            subject,
            pssm(rows),
            gap_open,
            gap_extend,
            end.match_end,
            end.query_end,
            score_in,
            forbidden,
        );
        let start = smith_waterman_find_start(
            subject,
            pssm(rows),
            gap_open,
            gap_extend,
            end.match_end,
            end.query_end,
            score_in,
            forbidden,
            scratch,
        );
        assert_eq!(start, before, "find-start scan to {score_in}, {case}");
        let plain = find_start_scan::<i32, _>(
            subject,
            pssm(rows),
            gap_open,
            gap_extend,
            end.match_end,
            end.query_end,
            score_in,
            forbidden,
            scratch,
        );
        assert_eq!(
            plain, before,
            "find-start scan in i32 to {score_in}, {case}"
        );
    }
}

// The cell that carries the gap along the row alone (`cell`) and the two cells of each pass give
// every value of the kernel before it (NCBI's order of the cell's terms, `before_k1`): short
// windows of every length up to 17 (odd and even, so the last single cell of both scans runs) with
// random forbidden ranges, square-matrix and PSSM rows, every gap cost of `GAP_COSTS` and random
// small ones; rows of few distinct scores, whose best cells tie; long windows with scores of a
// scaled matrix's range. One scratch serves every case, as in the round loop.
#[test]
fn the_k1_cell_gives_every_value_of_the_cell_before() {
    let mut random = XorShift(0x6b31_c0de_2026_1120);
    let mut scratch = SwScratch::new();
    for case in 0..6000 {
        let m = random.below(14) as usize;
        let n = random.below(18) as usize;
        let square = case % 2 == 0;
        // Few distinct scores (ties) in every third case.
        let range = if case % 3 == 0 {
            (-1, 2, 3)
        } else {
            (-9, 7, 1)
        };
        let rows = random_rows(&mut random, m, square, range);
        let subject: Vec<u8> = (0..n)
            .map(|_| random.below(BLASTAA_SIZE as u64) as u8)
            .collect();
        let gaps = if case % 4 == 0 {
            GAP_COSTS[(case / 4) % GAP_COSTS.len()]
        } else {
            (random.below(6) as i32 - 2, 1 + random.below(4) as i32)
        };
        let forbidden = if case % 5 < 2 {
            ForbiddenRanges::new()
        } else {
            random_ranges(&mut random, m, n)
        };
        let label = format!("case {case} gaps {gaps:?} m {m} n {n}");
        assert_same_as_before_k1(&subject, &rows, gaps, &forbidden, &mut scratch, &label);
    }
    for case in 0..24 {
        let m = 150 + random.below(400) as usize;
        let n = 150 + random.below(400) as usize;
        let rows = random_rows(&mut random, m, case % 2 == 0, (-4, 11, 32));
        let mut subject: Vec<u8> = (0..n)
            .map(|_| random.below(BLASTAA_SIZE as u64) as u8)
            .collect();
        // A stretch that scores well against the rows, so that the best cell is far from the ends.
        let stretch = m.min(n) / 3;
        for (offset, letter) in subject.iter_mut().skip(n / 4).take(stretch).enumerate() {
            let row = &rows[m / 4 + offset];
            let best = (0..BLASTAA_SIZE as u8)
                .max_by_key(|&l| row[usize::from(l)])
                .unwrap_or(1);
            if random.below(10) < 7 {
                *letter = best;
            }
        }
        let gaps = GAP_COSTS[case % GAP_COSTS.len()];
        let forbidden = if case % 3 == 0 {
            random_ranges(&mut random, m, n)
        } else {
            ForbiddenRanges::new()
        };
        let label = format!("long case {case} gaps {gaps:?} m {m} n {n}");
        assert_same_as_before_k1(&subject, &rows, gaps, &forbidden, &mut scratch, &label);
    }
}

// A row whose best score ties in many cells: both scans keep the cell that the kernel before them
// kept (the first of the score-only scan, the last new best of the find-start scan), at every
// parity of the tied cells' positions.
#[test]
fn tied_best_cells_stay_where_they_were() {
    let mut scratch = SwScratch::new();
    let rows = [row(&[1], 4, -4), row(&[1], 4, -4), row(&[2], 4, -4)];
    for n in 1..12 {
        for shift in 0..3 {
            // A at every third position from `shift`: many cells of the best score 4 (or 8, 12).
            let subject: Vec<u8> = (0..n)
                .map(|j| if (j + shift) % 3 == 0 { 1 } else { 2 })
                .collect();
            for gaps in [(0, 1), (11, 1), (-1, 1)] {
                let label = format!("n {n} shift {shift} gaps {gaps:?}");
                assert_same_as_before_k1(
                    &subject,
                    &rows,
                    gaps,
                    &ForbiddenRanges::new(),
                    &mut scratch,
                    &label,
                );
            }
        }
    }
}

// The negated gap costs are the costs' negations (`opaque` changes no value), and a row starts
// with NCBI's first gap along the row: `-max(open, 0)` plus the row's step is `-(open + extend)`.
#[test]
fn the_cells_costs_are_the_gap_costs() {
    for (open, extend) in GAP_COSTS {
        let costs = CellCosts::<i32>::new(open, extend);
        assert_eq!(costs.open_extend, -(open + extend));
        assert_eq!(costs.extend, -extend);
        assert_eq!(costs.row_extend, -(extend + open.min(0)));
        assert_eq!(costs.row_start + costs.row_extend, -(open + extend));
        assert_eq!(opaque(open), open);
    }
}

// `cmov_max` is `max` (a conditional move on x86-64), at the ends of the range too.
#[test]
fn the_conditional_move_is_the_maximum() {
    let values = [
        i32::MIN,
        i32::MIN + 1,
        -32768,
        -385,
        -1,
        0,
        1,
        384,
        i32::MAX - 1,
        i32::MAX,
    ];
    for &a in &values {
        for &b in &values {
            assert_eq!(cmov_max(a, b), a.max(b), "{a} {b}");
        }
    }
}

// A forbidden cell reads `i32::MIN` (`FORBIDDEN_LETTER`'s slot), whose diagonal floors at 0 as
// NCBI's `COMPO_SCORE_MIN` does for every column score, and a scan with forbidden ranges leaves
// the scratch's letters as the subject's, for the next scan of the round loop.
#[test]
fn forbidden_cells_floor_at_zero_and_leave_the_letters() {
    let slots = ScanRow::new();
    assert_eq!(slots.0[usize::from(FORBIDDEN_LETTER)], i32::MIN);
    for prev in [0, 1, 32767, 1 << 30, i32::MAX] {
        assert_eq!((prev + i32::MIN).max(0), COMPO_SCORE_MIN.max(0), "{prev}");
    }
    let rows: Vec<[i32; BLASTAA_SIZE]> = [1u8, 2, 3, 4, 5, 6]
        .iter()
        .map(|&letter| row(&[letter], 5, -3))
        .collect();
    let subject = [1u8, 2, 3, 4, 5, 6, 7, 1, 2, 3];
    let mut forbidden = ForbiddenRanges::new();
    forbidden.initialize(rows.len());
    forbidden.push(0, 3, 2, 4);
    forbidden.push(2, 6, -5, 30);
    let mut scratch = SwScratch::new();
    let end = smith_waterman_score_only(&subject, 6, pssm(&rows), 11, 1, &forbidden, &mut scratch);
    let expected: Vec<u8> = subject.iter().map(|&l| row_column(l) as u8).collect();
    assert_eq!(scratch.letters, expected);
    smith_waterman_find_start(
        &subject,
        pssm(&rows),
        11,
        1,
        end.match_end,
        end.query_end,
        end.score,
        &forbidden,
        &mut scratch,
    );
    assert_eq!(scratch.letters, expected);
}

// `Lane0`'s operations are `i32`'s (an addition that wraps, the signed maximum), at the ends of the
// range too, where the CPU has AVX2 (the scans that use `Lane0` run only there).
#[cfg(target_arch = "x86_64")]
#[test]
fn the_lane_operations_are_the_integer_operations() {
    if !std::is_x86_feature_detected!("avx2") {
        return;
    }
    let values = [
        i32::MIN,
        i32::MIN + 1,
        -32768,
        -385,
        -1,
        0,
        1,
        384,
        i32::MAX - 1,
        i32::MAX,
    ];
    for &a in &values {
        assert_eq!(Lane0::of(a).get(), a);
        for &b in &values {
            assert_eq!(Lane0::of(a).larger(Lane0::of(b)).get(), a.max(b), "{a} {b}");
            assert_eq!(
                Lane0::of(a).plus(Lane0::of(b)).get(),
                a.wrapping_add(b),
                "{a} {b}"
            );
        }
    }
}

// The scans read a letter's own column; a letter past NCBIstdaa fails a debug build (the tests),
// so an encoding change cannot read another column than NCBI's in silence.
#[test]
fn letters_read_their_own_column() {
    for letter in 0..BLASTAA_SIZE as u8 {
        assert_eq!(row_column(letter), usize::from(letter));
    }
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "is not NCBIstdaa")]
fn a_letter_past_ncbistdaa_fails_a_debug_build() {
    row_column(BLASTAA_SIZE as u8);
}
