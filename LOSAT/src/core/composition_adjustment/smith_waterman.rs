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
    /// The subject's letters as the scans read them (`row_column`, `ScanRow`); while a row with
    /// forbidden ranges is scanned, `FORBIDDEN_LETTER` at the positions they hold
    /// (`forbid_cells`).
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

/// The slots of a row as the scans read it (`ScanRow`): one for each value of a byte.
const ROW_SLOTS: usize = 256;

/// The letter of a forbidden cell (`forbid_cells`): its slot of `ScanRow` holds `i32::MIN`.
const FORBIDDEN_LETTER: u8 = u8::MAX;

// No letter of the alphabet reads the forbidden slot: the build fails if the alphabet reaches it.
const _: () = assert!(FORBIDDEN_LETTER as usize >= BLASTAA_SIZE);

/// The score row that a scan reads: the row's `BLASTAA_SIZE` scores, then zeros, one slot for each
/// value of a byte, so that a letter of `SwScratch::scan` (below `BLASTAA_SIZE`, `row_column`)
/// reads its own column with no bound check and no mask. One local array for each scan, whose
/// first `BLASTAA_SIZE` slots each row overwrites (`load`): the cell loop reads it from the stack,
/// which leaves the register of its address free.
///
/// The slot of `FORBIDDEN_LETTER` holds `i32::MIN`: a forbidden cell's diagonal is then
/// `prev + i32::MIN` for the score `prev` of the column in the row before, which is never
/// negative, so the diagonal is negative without overflow and `max(diagonal, 0)` is 0. NCBI gives
/// a forbidden cell the diagonal `COMPO_SCORE_MIN` before the same floor at 0, so the cell's
/// value is the same. Every row of a scan, with forbidden ranges or without, runs the same cell
/// with no test.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:303-317
struct ScanRow([i32; ROW_SLOTS]);

impl ScanRow {
    #[inline(always)]
    fn new() -> Self {
        let mut slots = [0; ROW_SLOTS];
        slots[usize::from(FORBIDDEN_LETTER)] = i32::MIN;
        Self(slots)
    }

    /// The scores of `row` in the first slots; the other slots keep theirs.
    #[inline(always)]
    fn load(&mut self, row: &[i32; BLASTAA_SIZE]) -> &[i32; ROW_SLOTS] {
        self.0[..BLASTAA_SIZE].copy_from_slice(row);
        &self.0
    }
}

impl SwScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// `scoreVector` of `length` columns, each `{0, -gapOpen}`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:81-85
    ///
    /// With it the subject's letters as the scans read them (`ScanRow`).
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
        self.letters.clear();
        // `row_column` is below `BLASTAA_SIZE`, so the cast keeps it.
        self.letters
            .extend(subject.iter().map(|&letter| row_column(letter) as u8));
        (&mut self.letters, &mut self.score_vector)
    }
}

/// The subject's letters and the columns of one scan (`SwScratch::scan`).
type SwScan<'a> = (&'a mut [u8], &'a mut [SwGapInfo]);

/// Writes `FORBIDDEN_LETTER` (`forbid` true) at the subject positions that the forbidden ranges
/// `ranges` of one query position hold, or the subject's own letters back (false), so that the
/// row's cells read `i32::MIN` there (`ScanRow`). NCBI's cell instead tests its subject position
/// against each range of the query position, both ends included; a position in any range is
/// forbidden. The ranges here are the same, so the same cells are forbidden.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:303-317
fn forbid_cells(letters: &mut [u8], subject: &[u8], ranges: &[(i32, i32)], forbid: bool) {
    let last = letters.len() as i64 - 1;
    for &(start, end) in ranges {
        let (start, end) = (i64::from(start).max(0), i64::from(end).min(last));
        if start <= end {
            let (start, end) = (start as usize, end as usize);
            if forbid {
                letters[start..=end].fill(FORBIDDEN_LETTER);
            } else {
                for (letter, &original) in
                    letters[start..=end].iter_mut().zip(&subject[start..=end])
                {
                    // `row_column` is below `BLASTAA_SIZE`, so the cast keeps it.
                    *letter = row_column(original) as u8;
                }
            }
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

/// The values of a scan's cells (`cell`): `i32`, whose maxima are `cmov_max`, or, in the scans
/// compiled for AVX2 (`score_only_avx2`, `find_start_avx2`), `Lane0`, whose maxima are `pmaxsd`.
/// Both compute the same integers: an addition wraps as `i32`'s does in the release build, and
/// `larger` is the signed maximum.
trait CellValue: Copy {
    fn of(value: i32) -> Self;
    fn get(self) -> i32;
    fn plus(self, other: Self) -> Self;
    fn larger(self, other: Self) -> Self;
}

impl CellValue for i32 {
    #[inline(always)]
    fn of(value: i32) -> Self {
        value
    }

    #[inline(always)]
    fn get(self) -> i32 {
        self
    }

    #[inline(always)]
    fn plus(self, other: Self) -> Self {
        self + other
    }

    #[inline(always)]
    fn larger(self, other: Self) -> Self {
        cmov_max(self, other)
    }
}

/// A cell's value in lane 0 of an SSE register, in the scans compiled for AVX2: each maximum is one
/// `pmaxsd` (ports 0 and 1) where `cmov_max` takes a `cmp` and a `cmovl` (`cmovl` on ports 0 and 6
/// only), so the cells' operations share the ALU ports less. The scalar cell is bound by its ALU
/// operations (ws, X4-mt kernel benchmark, 230.7 million cells: about 2.85 of them a cycle in every
/// form of the cell tried), and this one takes 1.090 G cycles where `i32`'s takes 1.387 G.
///
/// Invariant: a `Lane0` is made only where the CPU has AVX2, in `score_only_avx2` and
/// `find_start_avx2` (which run only where `is_x86_feature_detected!("avx2")` holds) or in a test
/// under the same detection. Its cells need SSE4.1 (`pmaxsd`) at least; AVX2 is the gate because
/// it is the one x86-64 fast level of the project's SIMD policy (2026-10-14), so each scan has one
/// compiled form for it, and only CPUs with AVX2 were measured. A CPU with SSE4.1 or AVX but not
/// AVX2 runs the `i32` cells, with the same values.
#[cfg(target_arch = "x86_64")]
#[derive(Debug, Clone, Copy)]
struct Lane0(core::arch::x86_64::__m128i);

#[cfg(target_arch = "x86_64")]
impl CellValue for Lane0 {
    #[inline(always)]
    fn of(value: i32) -> Self {
        // SAFETY: SSE2, which every x86-64 CPU has (the target's baseline).
        Self(unsafe { core::arch::x86_64::_mm_cvtsi32_si128(value) })
    }

    #[inline(always)]
    fn get(self) -> i32 {
        // SAFETY: SSE2, which every x86-64 CPU has (the target's baseline).
        unsafe { core::arch::x86_64::_mm_cvtsi128_si32(self.0) }
    }

    #[inline(always)]
    fn plus(self, other: Self) -> Self {
        // SAFETY: SSE2, which every x86-64 CPU has (the target's baseline).
        Self(unsafe { core::arch::x86_64::_mm_add_epi32(self.0, other.0) })
    }

    #[inline(always)]
    fn larger(self, other: Self) -> Self {
        // SAFETY: `pmaxsd` is SSE4.1, which AVX2 includes, and a `Lane0` exists only where the
        // CPU has AVX2 (the type's invariant).
        Self(unsafe { core::arch::x86_64::_mm_max_epi32(self.0, other.0) })
    }
}

/// The gap costs as the scans' cells use them (`cell`), from NCBI's `gapOpen` (`open` below) and
/// `gapExtend` (`extend`), negated, so that a cell adds them (`opaque`).
#[derive(Debug, Clone, Copy)]
struct CellCosts<V> {
    /// `-(open + extend)`: NCBI's `newGapCost`, negated.
    open_extend: V,
    /// `-extend`.
    extend: V,
    /// `-(extend + min(open, 0))`: what the gap along the row takes from one cell to the next in
    /// the form that `cell` computes it.
    row_extend: V,
    /// `-max(open, 0)`: the gap along the row before a row's first cell in that form.
    row_start: V,
}

impl<V: CellValue> CellCosts<V> {
    #[inline(always)]
    fn new(gap_open: i32, gap_extend: i32) -> Self {
        Self {
            open_extend: V::of(opaque(-(gap_open + gap_extend))),
            extend: V::of(opaque(-gap_extend)),
            row_extend: V::of(opaque(-(gap_extend + gap_open.min(0)))),
            row_start: V::of(-gap_open.max(0)),
        }
    }
}

/// `value`, with how it was computed hidden from the compiler: an empty `asm!` on x86-64, the value
/// itself on other targets. A cell adds the negated gap costs (`CellCosts`), and x86-64 adds into a
/// new register (`lea`) where a value must stay for later: the score of the row before for the
/// next cell's diagonal, and the cell's score without the gap along the row for the next cell's
/// gap. LLVM turns `x + (-c)` back into `x - c` when it sees the negation, and copies `x` before
/// such a subtraction (ws, X4-mt kernel benchmark, 230.7 million cells: 6.02 G instructions
/// against 6.25 G).
#[inline(always)]
fn opaque(value: i32) -> i32 {
    #[cfg(target_arch = "x86_64")]
    {
        let mut value = value;
        // SAFETY: no instruction (an assembler comment that names the register): the register
        // keeps the value; no memory, stack or flags are touched.
        unsafe {
            core::arch::asm!(
                "/* {value:e} */",
                value = inout(reg) value,
                options(pure, nomem, nostack, preserves_flags),
            );
        }
        value
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        value
    }
}

/// What a cell of the scans passes on to the next cell of its row (`cell`).
#[derive(Debug, Clone, Copy)]
struct RowCarry<V> {
    /// The cell's gap along the row (`G` at `cell`; NCBI's `prevScoreGapMatchSeq`).
    gap_row: V,
    /// The cell's score without the gap along the row (`D` at `cell`).
    no_row_gap: V,
    /// The cell's column in the row before (NCBI's `prevScoreNoGapMatchSeq`).
    prev_no_gap: V,
}

impl<V: CellValue> RowCarry<V> {
    /// Before the row's first cell (`D(-1) = 0` and `G(-1) = -max(open, 0)` at `cell`).
    #[inline(always)]
    fn start(costs: CellCosts<V>) -> Self {
        Self {
            gap_row: costs.row_start,
            no_row_gap: V::of(0),
            prev_no_gap: V::of(0),
        }
    }
}

/// One cell of the scans: NCBI's cell with its values computed in another order. `substitution`
/// is the row's score of the subject letter (`i32::MIN` for a forbidden cell, `ScanRow`); the
/// result is the cell's score.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:291-331
/// (the cells of the basic scans, :96-119 and :190-213, are those lines without the forbidden
/// test).
///
/// NCBI's cell `j` of a row takes its gap along the row `G(j) = max(S(j-1) - (open + extend),
/// G(j-1) - extend)` from the score `S(j-1)` of the cell before it, its gap down the column
/// `C(j)` from the column, and scores `S(j) = max(diag(j), 0, G(j), C(j))` (the diagonal is
/// `COMPO_SCORE_MIN` in a forbidden cell); its `if` tests are maxima of integers, whose order
/// changes no value. With `D(j) = max(diag(j), 0, C(j))`, `S(j-1) = max(D(j-1), G(j-1))`, so
/// `G(j) = max(D(j-1) - (open + extend), G(j-1) - (open + extend), G(j-1) - extend)`, which is
/// `max(D(j-1) - (open + extend), G(j-1) - (extend + min(open, 0)))` for every integer `open`.
/// Computed so, the gap along the row goes from one cell to the next through one addition and
/// one maximum, and the cell's score, which only the column, the best cell and the next cell's
/// `D` read, leaves that chain (in NCBI's order it runs through two maxima). Before a row's first
/// cell NCBI has `S(-1) = 0` and `G(-1) = -open`, so `G(0) = -(open + extend)`; the start
/// `D(-1) = 0`, `G(-1) = -max(open, 0)` (`RowCarry::start`) gives the same `G(0)`.
///
/// No value overflows where NCBI's order does not: `D(j-1) - (open + extend)` lies between
/// `0 - (open + extend)` and `S(j-1) - (open + extend)`, both of which NCBI computes;
/// `G(j-1) - (extend + min(open, 0))` is NCBI's `G(j-1) - extend` where `open >= 0`, and where
/// `open < 0` lies between `G(j-1) - extend` and `S(j-1) - (open + extend)`; at a row's start the
/// sum is NCBI's `-(open + extend)` and `-max(open, 0)` is `-open` or 0. The negated costs are
/// representable (each cost is above `i32::MIN`), and `x + (-c)` overflows where `x - c` does.
/// A forbidden cell's diagonal is `prev + i32::MIN` (`ScanRow`), negative for every column score
/// `prev >= 0`, so `max(diagonal, 0)` is NCBI's `max(COMPO_SCORE_MIN, 0)`, 0.
#[inline(always)]
fn cell<V: CellValue>(
    carry: &mut RowCarry<V>,
    column: &mut SwGapInfo,
    substitution: i32,
    costs: CellCosts<V>,
) -> i32 {
    carry.gap_row = carry
        .gap_row
        .plus(costs.row_extend)
        .larger(carry.no_row_gap.plus(costs.open_extend));
    let no_gap = V::of(column.no_gap);
    let continue_gap_score = V::of(column.gap_exists)
        .plus(costs.extend)
        .larger(no_gap.plus(costs.open_extend));
    let diagonal = carry.prev_no_gap.plus(V::of(substitution));
    let no_row_gap = diagonal.larger(V::of(0)).larger(continue_gap_score);
    let score = no_row_gap.larger(carry.gap_row).get();
    carry.no_row_gap = no_row_gap;
    carry.prev_no_gap = no_gap;
    column.no_gap = score;
    column.gap_exists = continue_gap_score.get();
    score
}

/// The larger of `a` and `b` by a conditional move on x86-64. The cell's maxima depend on the
/// scores, so a branch for them is mispredicted about every other cell, and LLVM's x86 back end
/// turns some of them into branches in this loop whatever the source's form: `max` and `if`
/// (the review's build took 1.1 to 1.9 times NCBI's time, PS11c-3 review S5), and also
/// `core::hint::select_unpredictable(a < b, b, a)`, built from which the gap along the row and
/// the maximum with 0 still became `jg` and `jle`, and the scans took 2.6 to 3.6 times as long
/// as with this `asm!` (PS11c-3 fix round 2, ws, 3 rounds: q40 x rps_mini `-evalue 1000` 18.7 s
/// against 7.0 s, NCBI 23.2 s; rpstest_le 67.7 s against 18.9 s, NCBI 24.7 s; WSSV x prof_rps
/// 21.5 s against 6.3 s, NCBI 6.5 s). With the cell of X4-mt (`cell`), `max` still made more
/// instructions and cycles (ws, X4-mt kernel benchmark: 7.86 G instructions and 1.67 G cycles
/// against 7.40 G and 1.58 G). NCBI's build keeps conditional moves and `pmaxsd`; on a CPU with
/// AVX2 the scans take `pmaxsd` too (`Lane0`), and this serves the other x86-64 CPUs. Other
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

/// One row of the score-only scan (subject positions up) after the cells of the rows before it,
/// whose best score is `best`: the row's first cell above it and its score, if any. `letters` and
/// `row` are as `ScanRow` says. Two cells for each pass of the loop: the loop's own count, test and
/// jump, and the moves of the carried values between registers, come once for two cells (ws,
/// X4-mt kernel benchmark: 6.25 G instructions and 1.42 G cycles against 7.17 G and 1.52 G for one
/// cell a pass).
#[inline(always)]
fn score_only_row<V: CellValue>(
    letters: &[u8],
    row: &[i32; ROW_SLOTS],
    score_vector: &mut [SwGapInfo],
    costs: CellCosts<V>,
    mut best: i32,
) -> Option<(i32, usize)> {
    let mut found = None;
    let mut carry = RowCarry::start(costs);
    let length = letters.len();
    let (pairs, last_letter) = letters.as_chunks::<2>();
    let (columns, last_column) = score_vector[..length].as_chunks_mut::<2>();
    for (pair, (&[letter0, letter1], [column0, column1])) in pairs.iter().zip(columns).enumerate() {
        let score = cell(&mut carry, column0, row[usize::from(letter0)], costs);
        if score > best {
            best = score;
            found = Some((best, 2 * pair));
        }
        let score = cell(&mut carry, column1, row[usize::from(letter1)], costs);
        if score > best {
            best = score;
            found = Some((best, 2 * pair + 1));
        }
    }
    if let (Some(&letter), Some(column)) = (last_letter.first(), last_column.first_mut()) {
        let score = cell(&mut carry, column, row[usize::from(letter)], costs);
        if score > best {
            best = score;
            found = Some((best, length - 1));
        }
    }
    found
}

/// One row of the find-start scan (subject positions from `match_end` down) after the rows before
/// it, whose best score is `best`: the row's last new best cell and its score, if any, and whether
/// the best score reached `score_in` (the scan stops there). Two cells for each pass of the loop,
/// as `score_only_row`.
#[inline(always)]
fn find_start_row<V: CellValue>(
    letters: &[u8],
    row: &[i32; ROW_SLOTS],
    score_vector: &mut [SwGapInfo],
    match_end: usize,
    costs: CellCosts<V>,
    mut best: i32,
    score_in: i32,
) -> (Option<(i32, usize)>, bool) {
    let mut found = None;
    let mut carry = RowCarry::start(costs);
    // The positions `0..=match_end` in pairs from the top; with an odd count, position 0 is the
    // last cell, alone.
    let (first_letter, pairs) = letters[..=match_end].as_rchunks::<2>();
    let (first_column, columns) = score_vector[..=match_end].as_rchunks_mut::<2>();
    let base = first_letter.len();
    for (pair, (&[letter0, letter1], [column0, column1])) in
        pairs.iter().zip(columns).enumerate().rev()
    {
        let score = cell(&mut carry, column1, row[usize::from(letter1)], costs);
        if score > best {
            best = score;
            found = Some((best, base + 2 * pair + 1));
        }
        if best >= score_in {
            return (found, true);
        }
        let score = cell(&mut carry, column0, row[usize::from(letter0)], costs);
        if score > best {
            best = score;
            found = Some((best, base + 2 * pair));
        }
        if best >= score_in {
            return (found, true);
        }
    }
    if let (Some(&letter), Some(column)) = (first_letter.first(), first_column.first_mut()) {
        let score = cell(&mut carry, column, row[usize::from(letter)], costs);
        if score > best {
            best = score;
            found = Some((best, 0));
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
// One call per scan, out of line: the row loop keeps its values in registers. On a CPU with AVX2
// the scan compiled for it runs (`score_only_avx2`); both compute the same values.
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
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has AVX2, which `score_only_avx2`'s `#[target_feature]` needs.
        return unsafe {
            score_only_avx2(
                subject,
                query_length,
                rows,
                gap_open,
                gap_extend,
                forbidden,
                scratch,
            )
        };
    }
    score_only_scan::<i32, F>(
        subject,
        query_length,
        rows,
        gap_open,
        gap_extend,
        forbidden,
        scratch,
    )
}

/// `smith_waterman_score_only` compiled for AVX2, its cells in `Lane0`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn score_only_avx2<'r, F>(
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
    score_only_scan::<Lane0, F>(
        subject,
        query_length,
        rows,
        gap_open,
        gap_extend,
        forbidden,
        scratch,
    )
}

/// The score-only scan (`smith_waterman_score_only`) with its cells in `V`.
#[inline(always)]
fn score_only_scan<'r, V: CellValue, F>(
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
    let (letters, score_vector) = scratch.scan(subject, gap_open);
    let costs = CellCosts::<V>::new(gap_open, gap_extend);
    let mut scan_row = ScanRow::new();
    let mut best_score = 0;
    let mut best_query_pos = 0;
    let mut best_match_seq_pos = 0;
    for query_pos in 0..query_length {
        let row = scan_row.load(rows(query_pos));
        let ranges = forbidden.at(query_pos);
        forbid_cells(letters, subject, ranges, true);
        let found = score_only_row::<V>(letters, row, score_vector, costs, best_score);
        forbid_cells(letters, subject, ranges, false);
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
// One call per scan, out of line: the row loop keeps its values in registers. On a CPU with AVX2
// the scan compiled for it runs (`find_start_avx2`); both compute the same values.
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
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has AVX2, which `find_start_avx2`'s `#[target_feature]` needs.
        return unsafe {
            find_start_avx2(
                subject, rows, gap_open, gap_extend, match_end, query_end, score_in, forbidden,
                scratch,
            )
        };
    }
    find_start_scan::<i32, F>(
        subject, rows, gap_open, gap_extend, match_end, query_end, score_in, forbidden, scratch,
    )
}

/// `smith_waterman_find_start` compiled for AVX2, its cells in `Lane0`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn find_start_avx2<'r, F>(
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
    find_start_scan::<Lane0, F>(
        subject, rows, gap_open, gap_extend, match_end, query_end, score_in, forbidden, scratch,
    )
}

/// The find-start scan (`smith_waterman_find_start`) with its cells in `V`.
#[inline(always)]
fn find_start_scan<'r, V: CellValue, F>(
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
    let (letters, score_vector) = scratch.scan(subject, gap_open);
    let costs = CellCosts::<V>::new(gap_open, gap_extend);
    let mut scan_row = ScanRow::new();
    let mut best_score = 0;
    let mut best_query_pos = 0;
    let mut best_match_seq_pos = 0;
    for query_pos in (0..=query_end).rev() {
        let row = scan_row.load(rows(query_pos));
        let ranges = forbidden.at(query_pos);
        forbid_cells(letters, subject, ranges, true);
        let (found, reached) = find_start_row::<V>(
            letters,
            row,
            score_vector,
            match_end,
            costs,
            best_score,
            score_in,
        );
        forbid_cells(letters, subject, ranges, false);
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
    /// With linked HSPs (`do_link_hsps`): `aSwScore` is at least `params->cutoff_s`.
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
/// (`BlastCompo_AlignmentNew`, blast_kappa.c:1858-1876).
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
    /// (`GapEditScriptNew` of 0 operations gives NULL, `s_NewAlignmentUsingXdrop` returns -1); NCBI
    /// frees the match's alignments and returns -1 (redo_alignment.c:1542-1548).
    Failed,
    /// NCBI's redo that does not end: a significant score of 0 in a window of two or more HSPs
    /// whose first query position scores the subject's first letter above 0 (see
    /// `smith_waterman_redo_window`); NCBI makes that one-cell alignment until an allocation fails.
    Endless,
}

/// The round loop of one window of `Blast_RedoOneMatchSmithWaterman`, after its composition
/// adjustment, for the square matrix of BLASTP's and TBLASTN's windows and RPS-BLAST's PSSM: the
/// score-only scan, the significance test (`significance`; before the query's first alignment also
/// `would_insert(evalue, score)`, NCBI's `BlastCompo_HeapWouldInsert` of the query's heap), the start
/// scan, the X-drop traceback from the start
/// (`traceback((query_start, query_end), (match_start, match_end), score)`, NCBI's
/// `new_xdrop_align`), and for a window of two or more HSPs the range of each alignment made
/// forbidden and the next round. The forbidden ranges are set to the window's query length
/// (NCBI's `Blast_ForbiddenRangesInitialize` with the queries' longest length, every range of the
/// window inside of it) and cleared.
///
/// The query's first alignment is the first alignment of the match's call for that query (NCBI
/// asks while the query's list of alignments is empty), not of the window: `query_had_none` says
/// that the query had no alignment when the window started (true for one window per query, as
/// BLASTP's and RPS-BLAST's protein windows are; TBLASTN's translated windows of one query follow
/// each other in one call), and the heap is asked only while that holds and the window has made
/// none.
///
/// The regime of a significant score of 0 is decided at the window's first significant round. When
/// the score of 0 passes the threshold (the e-value of 0 from `s_EvalueFromScore`, `searchsp` times
/// `K`, is below `cutoff_e`, or 0 is at least `cutoff_s`) and the window has two or more HSPs,
/// every round is significant (a later score is at least 0, its e-value at most that of 0, and the
/// heap is asked at the first round only), and every round of a positive score forbids its start
/// cell, which was free and whose diagonal step is positive (the reverse scan's first best comes
/// from the diagonal), so the loop makes every positive alignment and comes to a round of score 0.
/// That round's score-only scan ends at (0, 0) (`bestQueryPos` and `bestMatchSeqPos` are 0, a best
/// only for `>`), its find-start stops at (0, 0) (its test that `bestScore` is at least `score_in`,
/// here 0), and `s_SWFindFinalEndsUsingXdrop` aligns the 1 x 1 rectangle there: the first query
/// position's row against the subject's first letter (after SEG). A positive score makes one
/// substitution, prepends it and pushes `[0, 1]` onto row 0 again, which changes no score the scan
/// sees (every cell scores 0 already): the next round is the same, and NCBI makes it until an
/// allocation fails (`Endless`). Otherwise `ALIGN_EX` makes no operation, `GapEditScriptNew` of 0
/// operations gives NULL, the X-drop callback returns -1 and NCBI frees the match's alignments
/// (`Failed`). Both outcomes follow from the cell, so they are answered at the first round, without
/// the rounds. With one HSP the loop makes one round, whose score of 0 (all cells 0 or less) also
/// fails, at its traceback.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1438-1529
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/smith_waterman.c:437-446
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:428-432
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/gapinfo.c:61-62
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1856-1876
// Generic over the caller's rows, heap and traceback: one instance per redo (the square matrix's
// match body, `smith_waterman_match`, and RPS-BLAST's), each on the Smith-Waterman path only.
pub(crate) fn smith_waterman_redo_window<'r, F, H, T>(
    window: &SwWindow<'_, F>,
    significance: SwSignificance,
    query_had_none: bool,
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
        if let (Some(evalue), true) = (evalue, query_had_none && alignments.is_empty()) {
            alignment_is_significant = alignment_is_significant && would_insert(evalue, end.score);
        }
        if !alignment_is_significant {
            break;
        }
        // The regime of a significant score of 0 (see above), at the first round.
        if alignments.is_empty() && window.hspcnt > 1 && significance.test(0).1 {
            // `ALIGN_EX`'s cell (1, 1): row 0 (`query_offset` is `queryStart` less 1) against
            // `B[1]`, the subject's first letter (the sentinel when it has none).
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
        // An X-drop alignment without operations has no edit script (`GapEditScriptNew` of 0
        // operations gives NULL): the redo of the match fails and frees the match's alignments.
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
