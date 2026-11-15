//! NCBI's Smith-Waterman kernels of the composition-based redo (`-use_sw_tback`): the score and
//! right-hand ends of the locally optimal alignment (`Blast_SmithWatermanScoreOnly`), its
//! left-hand ends (`Blast_SmithWatermanFindStart`) and the ranges that later alignments of the
//! same window must not use (`Blast_ForbiddenRanges`), ported line by line (series PR, item
//! PS11c-3; DECISIONS 26 D7 of line PS).
//!
//! The kernels read the score rows through the caller's row source: a closure that gives the row
//! of a query position. NCBI's `positionSpecific` reads `matrix[queryPos]` (a PSSM: RPS-BLAST's
//! profiles now) and otherwise `matrix[query[queryPos]]` (a square matrix: blastp's and psiblast's
//! redo later); both are row sources of one kernel.
//!
//! NCBI's version without forbidden ranges (`BLbasicSmithWatermanScoreOnly`,
//! `BLSmithWatermanFindStart`) and the version with them (`BLspecialSmithWatermanScoreOnly`,
//! `BLspecialSmithWatermanFindStart`) differ only in the forbidden test of each cell, and a
//! position without forbidden ranges passes that test, so one scan of each kind serves both.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:544-575
//! ```c
//! int
//! Blast_SmithWatermanScoreOnly(int *score,
//!                              int *matchSeqEnd, int *queryEnd,
//!                              const Uint1 * subject_data, int subject_length,
//!                              const Uint1 * query_data, int query_length,
//!                              int **matrix,
//!                              int gapOpen,
//!                              int gapExtend,
//!                              int positionSpecific,
//!                              const Blast_ForbiddenRanges * forbiddenRanges )
//! {
//!     if (forbiddenRanges->isEmpty) {
//!         return BLbasicSmithWatermanScoreOnly(score, matchSeqEnd,
//!                                              queryEnd, subject_data,
//!                                              subject_length,
//!                                              query_data, query_length,
//!                                              matrix, gapOpen,
//!                                              gapExtend,
//!                                              positionSpecific);
//!     } else {
//!         return BLspecialSmithWatermanScoreOnly(score, matchSeqEnd,
//!                                                queryEnd, subject_data,
//!                                                subject_length,
//!                                                query_data,
//!                                                query_length, matrix,
//!                                                gapOpen, gapExtend,
//!                                                forbiddenRanges->numForbidden,
//!                                                forbiddenRanges->ranges,
//!                                                positionSpecific);
//!     }
//! }
//! ```

use super::compo_scale::COMPO_SCORE_MIN;
use crate::core::blast_util::BLASTAA_SIZE;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:38-43
// ```c
// typedef struct SwGapInfo {
//     int noGap;         /**< score if opening a gap */
//     int gapExists;     /**< score if continuing a gap */
// } SwGapInfo;
// ```
/// `SwGapInfo`: one column of the row that the scans keep.
#[derive(Debug, Clone, Copy)]
struct SwGapInfo {
    no_gap: i32,
    gap_exists: i32,
}

/// The memory of the scans (`scoreVector`, which NCBI allocates for each call).
#[derive(Debug, Default)]
pub(crate) struct SwScratch {
    score_vector: Vec<SwGapInfo>,
    /// For each subject position, whether a forbidden range of the row being scanned holds it
    /// (`blocked_row`); false between rows.
    blocked: Vec<bool>,
    /// The subject's letters, each at most `BLASTAA_SIZE - 1` (`padded_row`).
    letters: Vec<u8>,
}

/// The slots of a row as the scans read it (`padded_row`): a power of two above `BLASTAA_SIZE`.
const ROW_SLOTS: usize = 32;

/// A score row as the scans read it: the row's `BLASTAA_SIZE` scores, then zeros. The scans read
/// it at `letter & (ROW_SLOTS - 1)` for letters that `SwScratch::scan` brought to at most
/// `BLASTAA_SIZE - 1` (NCBIstdaa letters are below 28, and a letter above reads the last score,
/// as before), so the mask changes no letter and lets the compiler drop the bound check; the row
/// is a local array, which frees the register of its address in the cell loop.
#[inline(always)]
fn padded_row(row: &[i32; BLASTAA_SIZE]) -> [i32; ROW_SLOTS] {
    let mut padded = [0; ROW_SLOTS];
    padded[..BLASTAA_SIZE].copy_from_slice(row);
    padded
}

impl SwScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// `scoreVector` of `length` columns, each `{0, -gapOpen}`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:81-85
    /// ```c
    ///     newGapCost = gapOpen + gapExtend;
    ///     for (matchSeqPos = 0;  matchSeqPos < matchSeqLength;  matchSeqPos++) {
    ///         scoreVector[matchSeqPos].noGap = 0;
    ///         scoreVector[matchSeqPos].gapExists = -gapOpen;
    ///     }
    /// ```
    ///
    /// With it the subject's letters as the scans read them (`padded_row`) and the cleared marks
    /// of the forbidden ranges (`blocked_row`).
    fn scan(&mut self, subject: &[u8], gap_open: i32) -> SwScan<'_> {
        let length = subject.len();
        self.score_vector.clear();
        self.score_vector.resize(
            length,
            SwGapInfo {
                no_gap: 0,
                gap_exists: -gap_open,
            },
        );
        self.blocked.clear();
        self.blocked.resize(length, false);
        self.letters.clear();
        self.letters.extend(
            subject
                .iter()
                .map(|&letter| letter.min(BLASTAA_SIZE as u8 - 1)),
        );
        (&self.letters, &mut self.score_vector, &mut self.blocked)
    }
}

/// The subject's letters, the columns and the marks of one scan (`SwScratch::scan`).
type SwScan<'a> = (&'a [u8], &'a mut [SwGapInfo], &'a mut [bool]);

/// Marks (`value` true) or clears (false) the subject positions that the forbidden ranges `ranges`
/// of one query position hold, so that the row's cells test one flag each where NCBI's test every
/// range (`forbidden = TRUE` where `ranges[2f] <= matchSeqPos <= ranges[2f + 1]`, quoted at
/// `cell`): the same cells are forbidden.
fn blocked_row(blocked: &mut [bool], ranges: &[(i32, i32)], value: bool) {
    let last = blocked.len() as i64 - 1;
    for &(start, end) in ranges {
        let (start, end) = (i64::from(start).max(0), i64::from(end).min(last));
        if start <= end {
            blocked[start as usize..=end as usize].fill(value);
        }
    }
}

/// `Blast_ForbiddenRanges`: for each query position, the subject ranges `[start, end]` (both
/// ends included) that the scans do not align it to, in the order they were pushed.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/smith_waterman.h:44-51
/// ```c
/// typedef struct Blast_ForbiddenRanges {
///     int   isEmpty;             /**< True if there are no forbidden ranges */
///     int  *numForbidden;        /**< how many forbidden ranges at each
///                                     database position */
///     int **ranges;              /**< forbidden ranges for each database
///                                    position */
///     int   capacity;         /**< length of the query sequence */
/// } Blast_ForbiddenRanges;
/// ```
#[derive(Debug, Default)]
pub(crate) struct ForbiddenRanges {
    is_empty: bool,
    /// `ranges[f]` with its `numForbidden[f]` pairs; `capacity` is the vector's length.
    ranges: Vec<Vec<(i32, i32)>>,
}

impl ForbiddenRanges {
    pub(crate) fn new() -> Self {
        Self {
            is_empty: true,
            ranges: Vec::new(),
        }
    }

    /// `Blast_ForbiddenRangesInitialize`: `capacity` query positions without ranges (NCBI
    /// allocates them once for each call of the redo, with the query's length).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:473-500
    /// ```c
    ///     self->capacity  = capacity;
    ///     self->numForbidden = NULL;
    ///     self->ranges       = NULL;
    ///     self->isEmpty      = TRUE;
    ///     ...
    ///     for (f = 0;  f < capacity;  f++) {
    ///         self->numForbidden[f] = 0;
    ///         self->ranges[f]       = (int *) malloc(2 * sizeof(int));
    ///         ...
    ///         self->ranges[f][0]    = 0;
    ///         self->ranges[f][1]    = 0;
    ///     }
    /// ```
    pub(crate) fn initialize(&mut self, capacity: usize) {
        self.ranges.truncate(capacity);
        self.ranges.resize_with(capacity, Vec::new);
        self.clear();
    }

    /// `Blast_ForbiddenRangesClear`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:505-512
    /// ```c
    /// void Blast_ForbiddenRangesClear(Blast_ForbiddenRanges * self)
    /// {
    ///     int f;
    ///     for (f = 0;  f < self->capacity;  f++) {
    ///         self->numForbidden[f] = 0;
    ///     }
    ///     self->isEmpty = TRUE;
    /// }
    /// ```
    pub(crate) fn clear(&mut self) {
        for ranges in &mut self.ranges {
            ranges.clear();
        }
        self.is_empty = true;
    }

    /// `Blast_ForbiddenRangesPush`: the subject range `[match_start, match_end]` for each query
    /// position from `query_start` up to, not including, `query_end`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:516-540
    /// ```c
    ///     int f;
    ///     for (f = queryStart;  f < queryEnd;  f++) {
    ///         int last = 2 * self->numForbidden[f];
    ///         if (0 != last) {    /* we must resize the array */
    ///             int * new_ranges =
    ///                 realloc(self->ranges[f], (last + 2) * sizeof(int));
    ///             if (new_ranges == NULL)
    ///                 return -1;
    ///             self->ranges[f] = new_ranges;
    ///         }
    ///         self->ranges[f][last]     = matchStart;
    ///         self->ranges[f][last + 1] = matchEnd;
    ///
    ///         self->numForbidden[f]++;
    ///     }
    ///     self->isEmpty = FALSE;
    /// ```
    pub(crate) fn push(
        &mut self,
        query_start: i32,
        query_end: i32,
        match_start: i32,
        match_end: i32,
    ) {
        for f in query_start..query_end {
            self.ranges[f as usize].push((match_start, match_end));
        }
        self.is_empty = false;
    }

    /// The ranges of query position `query_pos` (none while `isEmpty`).
    #[inline(always)]
    fn at(&self, query_pos: usize) -> &[(i32, i32)] {
        if self.is_empty {
            &[]
        } else {
            &self.ranges[query_pos]
        }
    }
}

/// `Blast_SmithWatermanScoreOnly`'s results: the score and the right-hand ends (`*score`,
/// `*matchSeqEnd`, `*queryEnd`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwEnd {
    pub(crate) score: i32,
    pub(crate) match_end: usize,
    pub(crate) query_end: usize,
}

/// `Blast_SmithWatermanFindStart`'s results: the score and the left-hand ends (`*score_out`,
/// `*matchSeqStart`, `*queryStart`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SwStart {
    pub(crate) score: i32,
    pub(crate) match_start: usize,
    pub(crate) query_start: usize,
}

/// One cell of the scans: on entry `new_score` is the score of the cell before it in the scan's
/// row, on exit this cell's score; `substitution` is the row's score of the subject letter, or
/// `None` for a forbidden cell.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:291-331
/// ```c
///             /* testing scores with a gap in matchSeq, either starting
///              * a new gap or extending an existing gap */
///             if ((newScore = newScore - newGapCost) >
///                 (prevScoreGapMatchSeq = prevScoreGapMatchSeq - gapExtend))
///                 prevScoreGapMatchSeq = newScore;
///             /* testing scores with a gap in query, either starting a
///              * new gap or extending an existing gap */
///             if ((newScore = scoreVector[matchSeqPos].noGap - newGapCost) >
///                 (continueGapScore =
///                  scoreVector[matchSeqPos].gapExists - gapExtend))
///                 continueGapScore = newScore;
///             /* compute new score extending one position in matchSeq
///              * and query */
///             forbidden = FALSE;
///             for (f = 0;  f < numForbidden[queryPos];  f++) {
///                 if ((matchSeqPos >= forbiddenRanges[queryPos][2 * f]) &&
///                     (matchSeqPos <= forbiddenRanges[queryPos][2*f + 1])) {
///                     forbidden = TRUE;
///                     break;
///                 }
///             }
///             if (forbidden)
///                 newScore = COMPO_SCORE_MIN;
///             else
///                 newScore =
///                     prevScoreNoGapMatchSeq + matrixRow[matchSeq[matchSeqPos]];
///             if (newScore < 0)
///                 newScore = 0; /* Smith-Waterman locality condition */
///             /* test two alternatives */
///             if (newScore < prevScoreGapMatchSeq)
///                 newScore = prevScoreGapMatchSeq;
///             if (newScore < continueGapScore)
///                 newScore = continueGapScore;
///             prevScoreNoGapMatchSeq = scoreVector[matchSeqPos].noGap;
///             scoreVector[matchSeqPos].noGap = newScore;
///             scoreVector[matchSeqPos].gapExists = continueGapScore;
/// ```
/// (The basic scans, smith_waterman.c:96-119 and 190-213, are these lines without the forbidden
/// test.)
#[inline(always)]
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
    // NCBI's `if (a > b) b = a` and `if (s < t) s = t` are maxima of integers (`cmov_max`), so
    // the order of the cell's alternatives does not change its score; the gap along the row,
    // which the cell before computes, comes last (the other terms are ready before it).
    *prev_score_gap_match_seq = cmov_max(
        *prev_score_gap_match_seq - gap_extend,
        *new_score - new_gap_cost,
    );
    let continue_gap_score = cmov_max(column.gap_exists - gap_extend, column.no_gap - new_gap_cost);
    let diagonal = match substitution {
        None => COMPO_SCORE_MIN,
        Some(substitution) => *prev_score_no_gap_match_seq + substitution,
    };
    let score = cmov_max(
        cmov_max(cmov_max(diagonal, 0), continue_gap_score),
        *prev_score_gap_match_seq,
    );
    *prev_score_no_gap_match_seq = column.no_gap;
    column.no_gap = score;
    column.gap_exists = continue_gap_score;
    *new_score = score;
}

/// The larger of `a` and `b` by a conditional move on x86-64. The cell's maxima depend on the
/// scores, so a branch for them is mispredicted about every other cell, and LLVM's x86 back end
/// turns some of them into branches in this loop whatever the source's form: `max` and `if`
/// (the review's build took 1.1 to 1.9 times NCBI's time, PS11c-3 review S5), and also
/// `core::hint::select_unpredictable(a < b, b, a)`, built from which the gap along the row and
/// the maximum with 0 still became `jg` and `jle`, and the scans took 2.6 to 3.6 times as long
/// as with this `asm!` (PS11c-3 fix round 2, ws, 3 rounds: q40 x rps_mini `-evalue 1000` 18.7 s
/// against 7.0 s, NCBI 23.2 s; rpstest_le 67.7 s against 18.9 s, NCBI 24.7 s; WSSV x prof_rps
/// 21.5 s against 6.3 s, NCBI 6.5 s). NCBI's build keeps conditional moves and `pmaxsd`. Other
/// targets take `max`. The value is `max(a, b)` on every target.
#[inline(always)]
fn cmov_max(a: i32, b: i32) -> i32 {
    #[cfg(target_arch = "x86_64")]
    {
        let mut larger = a;
        // SAFETY: two instructions on the two registers given (`cmp`, then `cmovl`: `larger`
        // takes `b` where `larger < b`, signed); no memory or stack is touched, and the flags
        // that `cmp` sets are clobbered as `asm!` assumes by default.
        unsafe {
            core::arch::asm!(
                "cmp {larger:e}, {b:e}",
                "cmovl {larger:e}, {b:e}",
                larger = inout(reg) larger,
                b = in(reg) b,
                options(pure, nomem, nostack),
            );
        }
        larger
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        a.max(b)
    }
}

/// The row's score of subject letter `letter`, or `None` (`COMPO_SCORE_MIN` in the cell) where
/// `CHECKED` and the position is `blocked` (`blocked_row`). A row without ranges runs without the
/// test (NCBI's basic scans); `row` and `letter` are as `padded_row` says.
#[inline(always)]
fn substitution<const CHECKED: bool>(
    row: &[i32; ROW_SLOTS],
    blocked: bool,
    letter: u8,
) -> Option<i32> {
    if CHECKED && blocked {
        None
    } else {
        Some(row[usize::from(letter) & (ROW_SLOTS - 1)])
    }
}

/// One row of the score-only scan (subject positions up) after the cells of the rows before it,
/// whose best score is `best`: the row's first cell above it and its score, if any.
#[inline(always)]
fn score_only_row<const CHECKED: bool>(
    subject: &[u8],
    row: &[i32; ROW_SLOTS],
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
            substitution::<CHECKED>(row, blocked, letter),
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

/// One row of the find-start scan (subject positions from `match_end` down) after the rows before
/// it, whose best score is `best`: the row's last new best cell and its score, if any, and whether
/// the best score reached `score_in` (the scan stops there).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn find_start_row<const CHECKED: bool>(
    subject: &[u8],
    row: &[i32; ROW_SLOTS],
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
            substitution::<CHECKED>(row, blocked, letter),
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

/// `Blast_SmithWatermanScoreOnly`: the score of the locally optimal alignment of the
/// `query_length` query positions (rows `rows(0)`, `rows(1)`, ...) with `subject`, and its last
/// cell, the first of the best cells in the scan (query positions up, then subject positions
/// up), avoiding `forbidden`. `gap_open` and `gap_extend` are NCBI's (a gap of `k` letters costs
/// `gap_open + k * gap_extend`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:276-341
/// ```c
///     bestMatchSeqPos = 0;
///     bestQueryPos = 0;
///     bestScore = 0;
///     newGapCost = gapOpen + gapExtend;
///     ...
///     for (queryPos = 0;  queryPos < queryLength;  queryPos++) {
///         if (positionSpecific)
///             matrixRow = matrix[queryPos];
///         else
///             matrixRow = matrix[query[queryPos]];
///         newScore = 0;
///         prevScoreNoGapMatchSeq = 0;
///         prevScoreGapMatchSeq = -(gapOpen);
///         for (matchSeqPos = 0;  matchSeqPos < matchSeqLength;  matchSeqPos++) {
///             ...
///             if (newScore > bestScore) {
///                 bestScore = newScore;
///                 bestQueryPos = queryPos;
///                 bestMatchSeqPos = matchSeqPos;
///             }
///         }
///     }
///     free(scoreVector);
///     if (bestScore < 0)
///         bestScore = 0;
///     *matchSeqEnd = bestMatchSeqPos;
///     *queryEnd = bestQueryPos;
///     *score = bestScore;
/// ```
// One call per scan, out of line: the row loop keeps its values in registers.
#[inline(never)]
pub(crate) fn smith_waterman_score_only<'r, F>(
    subject: &[u8],
    query_length: usize,
    rows: F,
    gap_open: i32,
    gap_extend: i32,
    forbidden: &ForbiddenRanges,
    scratch: &mut SwScratch,
) -> SwEnd
where
    F: Fn(usize) -> &'r [i32; BLASTAA_SIZE],
{
    let (subject, score_vector, blocked) = scratch.scan(subject, gap_open);
    let mut best_score = 0;
    let mut best_query_pos = 0;
    let mut best_match_seq_pos = 0;
    for query_pos in 0..query_length {
        let row = &padded_row(rows(query_pos));
        let ranges = forbidden.at(query_pos);
        let gaps = (gap_open, gap_extend);
        let found = if ranges.is_empty() {
            score_only_row::<false>(subject, row, blocked, score_vector, gaps, best_score)
        } else {
            blocked_row(blocked, ranges, true);
            let found =
                score_only_row::<true>(subject, row, blocked, score_vector, gaps, best_score);
            blocked_row(blocked, ranges, false);
            found
        };
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

/// `Blast_SmithWatermanFindStart`: the first cell of the alignment that ends at `(query_end,
/// match_end)` with the score `score_in` (`Blast_SmithWatermanScoreOnly`'s), from the scan of the
/// reversed sequences that starts there (query positions down, then subject positions down) and
/// stops at the first cell where the best score reaches `score_in`, avoiding `forbidden`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:380-454
/// ```c
///     bestMatchSeqPos = 0;
///     bestQueryPos = 0;
///     bestScore = 0;
///     newGapCost = gapOpen + gapExtend;
///     ...
///     for (queryPos = queryEnd;  queryPos >= 0;  queryPos--) {
///         if (positionSpecific)
///             matrixRow = matrix[queryPos];
///         else
///             matrixRow = matrix[query[queryPos]];
///         newScore = 0;
///         prevScoreNoGapMatchSeq = 0;
///         prevScoreGapMatchSeq = -(gapOpen);
///         for (matchSeqPos = matchSeqEnd;  matchSeqPos >= 0;  matchSeqPos--) {
///             ...
///             if (newScore > bestScore) {
///                 bestScore = newScore;
///                 bestQueryPos = queryPos;
///                 bestMatchSeqPos = matchSeqPos;
///             }
///             if (bestScore >= score_in)
///                 break;
///         }
///         if (bestScore >= score_in)
///             break;
///     }
///     free(scoreVector);
///     if (bestScore < 0)
///         bestScore = 0;
///     *matchSeqStart = bestMatchSeqPos;
///     *queryStart = bestQueryPos;
///     *score_out = bestScore;
/// ```
// One call per scan, out of line: the row loop keeps its values in registers.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn smith_waterman_find_start<'r, F>(
    subject: &[u8],
    rows: F,
    gap_open: i32,
    gap_extend: i32,
    match_end: usize,
    query_end: usize,
    score_in: i32,
    forbidden: &ForbiddenRanges,
    scratch: &mut SwScratch,
) -> SwStart
where
    F: Fn(usize) -> &'r [i32; BLASTAA_SIZE],
{
    // The scan reads subject positions `match_end` down to 0; a subject without letters has no
    // alignment to start (no caller asks: it has no HSP).
    if subject.is_empty() {
        return SwStart {
            score: 0,
            match_start: 0,
            query_start: 0,
        };
    }
    let (subject, score_vector, blocked) = scratch.scan(subject, gap_open);
    let mut best_score = 0;
    let mut best_query_pos = 0;
    let mut best_match_seq_pos = 0;
    for query_pos in (0..=query_end).rev() {
        let row = &padded_row(rows(query_pos));
        let ranges = forbidden.at(query_pos);
        let gaps = (gap_open, gap_extend);
        let (found, reached) = if ranges.is_empty() {
            find_start_row::<false>(
                subject,
                row,
                blocked,
                score_vector,
                match_end,
                gaps,
                best_score,
                score_in,
            )
        } else {
            blocked_row(blocked, ranges, true);
            let outcome = find_start_row::<true>(
                subject,
                row,
                blocked,
                score_vector,
                match_end,
                gaps,
                best_score,
                score_in,
            );
            blocked_row(blocked, ranges, false);
            outcome
        };
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

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;
