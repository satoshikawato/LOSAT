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
    fn score_vector(&mut self, length: usize, gap_open: i32) -> &mut [SwGapInfo] {
        self.score_vector.clear();
        self.score_vector.resize(
            length,
            SwGapInfo {
                no_gap: 0,
                gap_exists: -gap_open,
            },
        );
        &mut self.score_vector
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
    let opened = *new_score - new_gap_cost;
    *prev_score_gap_match_seq -= gap_extend;
    if opened > *prev_score_gap_match_seq {
        *prev_score_gap_match_seq = opened;
    }
    let opened = column.no_gap - new_gap_cost;
    let mut continue_gap_score = column.gap_exists - gap_extend;
    if opened > continue_gap_score {
        continue_gap_score = opened;
    }
    let mut score = match substitution {
        None => COMPO_SCORE_MIN,
        Some(substitution) => *prev_score_no_gap_match_seq + substitution,
    };
    if score < 0 {
        score = 0;
    }
    if score < *prev_score_gap_match_seq {
        score = *prev_score_gap_match_seq;
    }
    if score < continue_gap_score {
        score = continue_gap_score;
    }
    *prev_score_no_gap_match_seq = column.no_gap;
    column.no_gap = score;
    column.gap_exists = continue_gap_score;
    *new_score = score;
}

/// The row's score of subject position `match_pos` (letter `letter`), or `None` where a range of
/// `ranges` holds the position.
#[inline(always)]
fn substitution(
    row: &[i32; BLASTAA_SIZE],
    ranges: &[(i32, i32)],
    match_pos: usize,
    letter: u8,
) -> Option<i32> {
    let position = match_pos as i64;
    let forbidden = ranges
        .iter()
        .any(|&(start, end)| position >= i64::from(start) && position <= i64::from(end));
    (!forbidden).then(|| row[letter as usize])
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
    let score_vector = scratch.score_vector(subject.len(), gap_open);
    let new_gap_cost = gap_open + gap_extend;
    let mut best_score = 0;
    let mut best_query_pos = 0;
    let mut best_match_seq_pos = 0;
    for query_pos in 0..query_length {
        let row = rows(query_pos);
        let ranges = forbidden.at(query_pos);
        let mut new_score = 0;
        let mut prev_score_no_gap_match_seq = 0;
        let mut prev_score_gap_match_seq = -gap_open;
        for (match_seq_pos, (&letter, column)) in
            subject.iter().zip(score_vector.iter_mut()).enumerate()
        {
            cell(
                &mut new_score,
                &mut prev_score_gap_match_seq,
                &mut prev_score_no_gap_match_seq,
                column,
                substitution(row, ranges, match_seq_pos, letter),
                new_gap_cost,
                gap_extend,
            );
            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
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
    let score_vector = scratch.score_vector(subject.len(), gap_open);
    let new_gap_cost = gap_open + gap_extend;
    let mut best_score = 0;
    let mut best_query_pos = 0;
    let mut best_match_seq_pos = 0;
    'rows: for query_pos in (0..=query_end).rev() {
        let row = rows(query_pos);
        let ranges = forbidden.at(query_pos);
        let mut new_score = 0;
        let mut prev_score_no_gap_match_seq = 0;
        let mut prev_score_gap_match_seq = -gap_open;
        for match_seq_pos in (0..=match_end).rev() {
            cell(
                &mut new_score,
                &mut prev_score_gap_match_seq,
                &mut prev_score_no_gap_match_seq,
                &mut score_vector[match_seq_pos],
                substitution(row, ranges, match_seq_pos, subject[match_seq_pos]),
                new_gap_cost,
                gap_extend,
            );
            if new_score > best_score {
                best_score = new_score;
                best_query_pos = query_pos;
                best_match_seq_pos = match_seq_pos;
            }
            if best_score >= score_in {
                break 'rows;
            }
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
