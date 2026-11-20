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

use super::compo_scale::COMPO_SCORE_MIN;
use crate::core::blast_util::BLASTAA_SIZE;
use crate::utils::glibc_fma_libm::GlibcLibm;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:38-43
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
    /// The subject's letters as the scans read them (`row_column`, `padded_row`).
    letters: Vec<u8>,
}

/// The column of a score row that subject letter `letter` reads, in the scans (`SwScratch::scan`)
/// and in the redo's test of the first cell (`redo1::smith_waterman`): the letter. The letters
/// are NCBIstdaa, below `BLASTAA_SIZE` (the query's encoding, masked letters as X), which a debug
/// build checks, so that an encoding change fails the tests; a release build keeps the read
/// inside the row (at its last column) where NCBI would read past its row.
#[inline(always)]
pub(crate) fn row_column(letter: u8) -> usize {
    debug_assert!(
        usize::from(letter) < BLASTAA_SIZE,
        "subject letter {letter} is not NCBIstdaa"
    );
    usize::from(letter).min(BLASTAA_SIZE - 1)
}

/// The slots of a row as the scans read it (`padded_row`): a power of two above `BLASTAA_SIZE`.
const ROW_SLOTS: usize = 32;

/// A score row as the scans read it: the row's `BLASTAA_SIZE` scores, then zeros. The scans read
/// it at `letter & (ROW_SLOTS - 1)` for letters that `SwScratch::scan` took through `row_column`
/// (below `BLASTAA_SIZE`), so the mask changes no letter and lets the compiler drop the bound
/// check; the row is a local array, which frees the register of its address in the cell loop.
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
        // `row_column` is below `BLASTAA_SIZE`, so the cast keeps it.
        self.letters
            .extend(subject.iter().map(|&letter| row_column(letter) as u8));
        (&self.letters, &mut self.score_vector, &mut self.blocked)
    }
}

/// The subject's letters, the columns and the marks of one scan (`SwScratch::scan`).
type SwScan<'a> = (&'a [u8], &'a mut [SwGapInfo], &'a mut [bool]);

/// Marks (`value` true) or clears (false) the subject positions that the forbidden ranges `ranges`
/// of one query position hold, so that the row's cells test one flag each where NCBI's test every
/// range (`forbidden = TRUE` where `ranges[2f] <= matchSeqPos <= ranges[2f + 1]`, cited at
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
/// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/smith_waterman.h:44-53
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
    pub(crate) fn initialize(&mut self, capacity: usize) {
        self.ranges.truncate(capacity);
        self.ranges.resize_with(capacity, Vec::new);
        self.clear();
    }

    /// `Blast_ForbiddenRangesClear`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:505-512
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:290-331
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

/// The memory of NCBI's Smith-Waterman round loop that the redo keeps from one window to the next
/// (`smith_waterman_redo_window`): the forbidden ranges (NCBI's thread-local `forbidden_tld`) and
/// the scans' memory (`scoreVector`).
#[derive(Debug)]
pub(crate) struct SwRounds {
    forbidden: ForbiddenRanges,
    scan: SwScratch,
}

impl Default for SwRounds {
    // No ranges (`ForbiddenRanges::new`: `isEmpty` true, as `clear` leaves them at each window).
    fn default() -> Self {
        Self {
            forbidden: ForbiddenRanges::new(),
            scan: SwScratch::new(),
        }
    }
}

/// What the Smith-Waterman redo keeps from one match to the next: the round loop's memory, one
/// owner for BLASTP's and RPS-BLAST's redo (`rounds`), and the caller's own memory of its windows
/// (`window`: RPS-BLAST's copy of the subject between NCBI's two sentinels, BLASTP's rows of the
/// adjusted square matrix and its X-drop traceback's memory).
#[derive(Debug, Default)]
pub(crate) struct SwRedoScratch<W> {
    pub(crate) rounds: SwRounds,
    pub(crate) window: W,
}

impl<W: Default> SwRedoScratch<W> {
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

/// `s_EvalueFromScore`: the e-value of a Smith-Waterman score in the redo's significance test.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:943-947
pub(crate) fn evalue_from_score(score: i32, lambda: f64, log_k: f64, searchsp: f64) -> f64 {
    searchsp * (-(lambda * f64::from(score)) + log_k).glibc_exp()
}

/// The significance test of a round (redo_alignment.c:1456-1474, cited at
/// `smith_waterman_redo_window`).
#[derive(Debug, Clone, Copy)]
pub(crate) enum SwSignificance {
    /// With linked HSPs (`do_link_hsps`): `aSwScore >= params->cutoff_s`.
    CutoffScore(i32),
    /// Otherwise the e-value of the score (`s_EvalueFromScore` with the redo's scaled `Lambda`,
    /// `logK` and the query's `eff_search_space`) below `params->cutoff_e` (`-evalue`), and before
    /// the query's first alignment the heap test (`BlastCompo_HeapWouldInsert`).
    Evalue {
        lambda: f64,
        log_k: f64,
        searchsp: f64,
        cutoff_e: f64,
    },
}

impl SwSignificance {
    /// The e-value of `score` (`None` with linked HSPs) and whether it passes the threshold.
    fn test(self, score: i32) -> (Option<f64>, bool) {
        match self {
            Self::CutoffScore(cutoff_s) => (None, score >= cutoff_s),
            Self::Evalue {
                lambda,
                log_k,
                searchsp,
                cutoff_e,
            } => {
                let evalue = evalue_from_score(score, lambda, log_k, searchsp);
                (Some(evalue), evalue < cutoff_e)
            }
        }
    }
}

/// One window of the Smith-Waterman redo as the round loop reads it: NCBI's `subject.data - 1`
/// (the sentinel before the window's letters, then the `subject_length` letters and a sentinel;
/// `ALIGN_EX` reads it so), the window's query length and the score rows of its positions (a PSSM's
/// `pssm[queryPos]`, or a square matrix's `matrix[query[queryPos]]`: one row source of the kernels),
/// the gap costs and the window's count of preliminary HSPs (`window->hspcnt`).
pub(crate) struct SwWindow<'a, F> {
    pub(crate) subject_from_sentinel: &'a [u8],
    pub(crate) subject_length: usize,
    pub(crate) query_length: usize,
    pub(crate) rows: F,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) hspcnt: usize,
}

/// The traceback of a Smith-Waterman alignment (`s_SWFindFinalEndsUsingXdrop`): its extents from
/// the start in the query and the subject (`*queryAlignmentExtent`, `*matchAlignmentExtent`), the
/// X-drop alignment's score and its edit script (`pssm_search::gapalign::pssm_sw_xdrop_traceback`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SwXdropTraceback {
    pub(crate) query_extent: i32,
    pub(crate) match_extent: i32,
    pub(crate) score: i32,
    pub(crate) edit_script: Vec<crate::common::GapEditOp>,
}

/// An alignment of the round loop, in the window's coordinates (`queryStart`, the updated
/// `queryEnd`, `matchStart`, the updated `matchEnd` of `s_NewAlignmentUsingXdrop`), with the X-drop
/// score and edit script; the caller adds the window's offsets, its frame and matrix rule
/// (`BlastCompo_AlignmentNew`, blast_kappa.c:1856-1873).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SwAlignment {
    pub(crate) score: i32,
    pub(crate) query_start: i32,
    pub(crate) query_end: i32,
    pub(crate) match_start: i32,
    pub(crate) match_end: i32,
    pub(crate) edit_script: Vec<crate::common::GapEditOp>,
}

/// What the round loop of one window gives (`smith_waterman_redo_window`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum SwWindowRedo {
    /// The window's alignments in the order the rounds made them (NCBI prepends each one to
    /// `alignments[query_index]`, so its list starts with the last).
    Alignments(Vec<SwAlignment>),
    /// NCBI's failed redo: a significant alignment whose X-drop traceback has no operation
    /// (`GapEditScriptNew(0)` is NULL, `s_NewAlignmentUsingXdrop` returns -1); NCBI frees the
    /// match's alignments and returns -1 (redo_alignment.c:1542-1548).
    Failed,
    /// NCBI's redo that does not end: a significant score of 0 in a window of two or more HSPs
    /// whose first query position scores the subject's first letter above 0 (see
    /// `smith_waterman_redo_window`); NCBI makes that one-cell alignment until an allocation fails.
    Endless,
}

/// The round loop of one window of `Blast_RedoOneMatchSmithWaterman`, after its composition
/// adjustment, for BLASTP's square matrix and RPS-BLAST's PSSM: the score-only scan, the
/// significance test (`significance`; before the query's first alignment also
/// `would_insert(evalue, score)`, NCBI's `BlastCompo_HeapWouldInsert` of the query's heap, which
/// the window's first alignment is), the start scan, the X-drop traceback from the start
/// (`traceback((query_start, query_end), (match_start, match_end), score)`, NCBI's
/// `new_xdrop_align`), and for a window of two or more HSPs the range of each alignment made
/// forbidden and the next round. The forbidden ranges are set to the window's query length
/// (NCBI's `Blast_ForbiddenRangesInitialize` with the queries' longest length, every range of the
/// window inside of it) and cleared.
///
/// The regime of a significant score of 0 is decided at the window's first significant round. When
/// the score of 0 passes the threshold (`s_EvalueFromScore(0, ..) = searchsp * K < cutoff_e`, or
/// `0 >= cutoff_s`) and the window has two or more HSPs, every round is significant (a later score
/// is at least 0, its e-value at most that of 0, and the heap is asked at the first round only),
/// and every round of a positive score forbids its start cell, which was free and whose diagonal
/// step is positive (the reverse scan's first best comes from the diagonal), so the loop makes
/// every positive alignment and comes to a round of score 0. That round's score-only scan ends at
/// (0, 0) (`bestQueryPos = bestMatchSeqPos = 0`, a best only for `>`), its find-start stops at
/// (0, 0) (`bestScore >= score_in` with 0), and `s_SWFindFinalEndsUsingXdrop` aligns the 1 x 1
/// rectangle there: the first query position's row against the subject's first letter (after
/// SEG). A positive score makes one substitution, prepends it and pushes `[0, 1]` onto row 0 again,
/// which changes no score the scan sees (every cell scores 0 already): the next round is the same,
/// and NCBI makes it until an allocation fails (`Endless`). Otherwise `ALIGN_EX` makes no
/// operation, `GapEditScriptNew(0)` is NULL, the X-drop callback returns -1 and NCBI frees the
/// match's alignments (`Failed`). Both outcomes follow from the cell, so they are answered at the
/// first round, without the rounds. With one HSP the loop makes one round, whose score of 0 (all
/// cells 0 or less) also fails, at its traceback.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1438-1529
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:437-446
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:428-432
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/gapinfo.c:61-62
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1856-1876
// Generic over the caller's rows, heap and traceback: one instance per redo (BLASTP's, RPS-BLAST's),
// each on the Smith-Waterman path only.
pub(crate) fn smith_waterman_redo_window<'r, F, H, T>(
    window: &SwWindow<'_, F>,
    significance: SwSignificance,
    mut would_insert: H,
    mut traceback: T,
    rounds: &mut SwRounds,
) -> SwWindowRedo
where
    F: Fn(usize) -> &'r [i32; BLASTAA_SIZE] + Copy,
    H: FnMut(f64, i32) -> bool,
    T: FnMut((usize, usize), (usize, usize), i32) -> SwXdropTraceback,
{
    let letters = &window.subject_from_sentinel[1..=window.subject_length];
    let rows = window.rows;
    rounds.forbidden.initialize(window.query_length);
    rounds.forbidden.clear();
    let mut alignments: Vec<SwAlignment> = Vec::new();
    loop {
        let end = smith_waterman_score_only(
            letters,
            window.query_length,
            rows,
            window.gap_open,
            window.gap_extend,
            &rounds.forbidden,
            &mut rounds.scan,
        );
        let (evalue, mut alignment_is_significant) = significance.test(end.score);
        if let (Some(evalue), true) = (evalue, alignments.is_empty()) {
            alignment_is_significant = alignment_is_significant && would_insert(evalue, end.score);
        }
        if !alignment_is_significant {
            break;
        }
        // The regime of a significant score of 0 (see above), at the first round.
        if alignments.is_empty() && window.hspcnt > 1 && significance.test(0).1 {
            // `ALIGN_EX`'s cell (1, 1): row 0 (`query_offset = queryStart - 1`) against `B[1]`, the
            // subject's first letter (the sentinel when it has none).
            let first_letter = window.subject_from_sentinel.get(1).copied().unwrap_or(0);
            let first_cell = if window.query_length == 0 {
                0
            } else {
                rows(0)[row_column(first_letter)]
            };
            return if first_cell > 0 {
                SwWindowRedo::Endless
            } else {
                SwWindowRedo::Failed
            };
        }
        let start = smith_waterman_find_start(
            letters,
            rows,
            window.gap_open,
            window.gap_extend,
            end.match_end,
            end.query_end,
            end.score,
            &rounds.forbidden,
            &mut rounds.scan,
        );
        let traced = traceback(
            (start.query_start, end.query_end),
            (start.match_start, end.match_end),
            end.score,
        );
        // An X-drop alignment without operations has no edit script (`GapEditScriptNew(0)` is
        // NULL): the redo of the match fails and frees the match's alignments.
        if traced.edit_script.is_empty() {
            return SwWindowRedo::Failed;
        }
        let query_start = start.query_start as i32;
        let match_start = start.match_start as i32;
        let query_end = query_start + traced.query_extent;
        let match_end = match_start + traced.match_extent;
        alignments.push(SwAlignment {
            score: traced.score,
            query_start,
            query_end,
            match_start,
            match_end,
            edit_script: traced.edit_script,
        });
        if window.hspcnt <= 1 {
            break;
        }
        // A significant score of 0 is decided above, so this round's score is positive and its
        // range forbids its free start cell: the loop ends.
        debug_assert!(end.score > 0);
        rounds
            .forbidden
            .push(query_start, query_end, match_start, match_end);
    }
    SwWindowRedo::Alignments(alignments)
}

#[cfg(test)]
#[path = "smith_waterman_tests.rs"]
mod tests;
