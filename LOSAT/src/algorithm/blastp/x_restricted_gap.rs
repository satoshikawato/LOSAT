//! EXPERIMENT (LOSAT_X_BPGAP): the restricted gapped alignment of the blastp preliminary stage
//! (NCBI `s_RestrictedGappedAlign`) with its score rows in separate arrays, eight or sixteen cells
//! of a row at a time (x86_64 with AVX2; with LOSAT_X_BPGAP128 the sixteen-bit rows eight cells at
//! a time in 128 bits on wasm32 simd128 and x86_64 without AVX2; elsewhere the reference function
//! runs).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1107-1122,1143-1147,1171
//!
//! What the reference computes. Row `a` of the alignment holds the cells `first_b_index ..
//! b_size`; cell `b` takes the diagonal `D(b)` = the previous row's `best` of cell `b - 1` plus
//! the matrix score of query letter `a` against subject letter `b` (`MININT` for the first cell of
//! the row), the column gap `best_gap` of cell `b` when `b` is a multiple of `RESTRICT_SIZE`
//! (`b_gap` starts at the first multiple at or after `first_b_index` and steps by
//! `RESTRICT_SIZE`), and the row gap `score_gap_row` only in rows `a` that are multiples of
//! `RESTRICT_SIZE`. A cell is dropped when `best_score - score > x_dropoff`.
//!
//! Why this module returns the same `(a_offset, b_offset, best_score)` for every input.
//!
//! - Two rows instead of one array updated in place. The reference reads `score_array[b].best`
//!   (for `next_score`) before it writes it in the same row, so it reads the previous row; here
//!   the previous row is its own array (`h[(a - 1) % 2]`). Every cell the next row reads (cells
//!   `first_b_index - 1 .. b_size - 1` of this row) is written by this row's loop, by the band
//!   extension after it or by the `MININT` cell after that, as in the reference; the diagonal of
//!   the first cell is `MININT` (the reference starts each row with `score = MININT`), so the
//!   value stored left of the band is never used. `best_gap` is read only in the columns that
//!   are multiples of `RESTRICT_SIZE` (`b_index == b_gap`), so only those are kept (`e10[b / 10]`).
//! - Rows where no gap opens in A (nine of ten). A cell's score depends on the previous row and
//!   on its own column gap only; along the row only the X-drop test links the cells, through
//!   `best_score`. When cell `b` is tested, `best_score` is the larger of its value at the start
//!   of the row and the scores of the cells left of `b` in this row: a dropped cell has
//!   `score < best_score - x_dropoff <= best_score` (`x_dropoff >= gap_open + gap_extend >= 0`),
//!   so it could not have raised `best_score` either. That is a prefix maximum, computed for eight
//!   cells at a time; when no cell of the eight exceeds `best_score`, it is `best_score` for all
//!   of them and is not computed. `best_score` changes at the first cell of the row that holds the row's
//!   largest score, if that score exceeds `best_score` (`score > best_score` is strict), which
//!   gives `a_offset` and `b_offset`. `first_b_index` passes the dropped cells at the start of the
//!   row, `last_b_index` is the last cell kept, and the score after the loop
//!   (`score_gap_row = score`) is the diagonal into cell `b_size`.
//! - Rows where a gap may open in A. The row gap is the prefix maximum of
//!   `H0(j) - gap_open_extend - (b - 1 - j) * gap_extend` over the cells `j < b` of the row, where
//!   `H0` is the cell's score without the row gap (diagonal, and column gap in the columns that
//!   are multiples of `RESTRICT_SIZE`). Let `thr = best_score - x_dropoff` at the moment a cell is
//!   tested; it never decreases. The reference takes the row gap from the kept cells only and does
//!   not lower it at a dropped cell; everything this module computes differently from that
//!   (contributions of dropped cells, and of cells left of a dropped cell) is below `thr`, because
//!   a dropped cell's score, which is at least its row gap, is below `thr`. A value below `thr`
//!   cannot be the score of a kept cell, cannot change which cells are dropped, and cannot start
//!   the band extension after the row (`score_gap_row >= best_score - x_dropoff`), so the kept
//!   cells, their scores, the column gaps written from them and the extension are the reference's.
//!   A row gap that continues into a cell and opens again from it is never larger than the row gap
//!   carried on (`gap_open >= 0`), which is why `H0` may stand for the score.
//! - The values of these rows are 32-bit integers as in the reference: no rebasing, no saturation
//!   (the sixteen-bit rows below are a separate kernel). The matrix row is read from a table of
//!   16-bit entries built from the reference's own score function (`blastp_score_from_row`), so
//!   every score is the reference's (BLOSUM62's `BLAST_SCORE_MIN` = -32768 for the gap letter
//!   fits). A matrix with an entry outside 16 bits, a letter code of 28 or more, a negative gap
//!   cost, a gap cost above 65,536 or an X-drop above 2^28 are not handled here; the reference
//!   function runs instead.
//! - Sixteen-bit rows (`LOSAT_X_BPGAP16`, `restricted_rows16_avx2`). The same rows with the cells,
//!   the column gaps (`e16`, one entry per column; the columns that are not multiples of
//!   `RESTRICT_SIZE` hold `i16::MIN` and are never written) and the row gap in 16-bit lanes,
//!   sixteen cells per step, with saturating adds and subtracts; `i16::MIN` stands where the
//!   32-bit rows store `MININT`. They run where `fits16` holds (X-drop and `gap_open_extend` at
//!   most 4,096, `gap_extend` at most 256, largest matrix entry `max` at most 4,096) and hold the
//!   32-bit rows' values:
//!   1. Upward: let `best` be `best_score` at the start of a row. A cell of the previous row is
//!      at most `max(best, max)`: a kept cell is at most `best`; a row gap, and so a cell of the
//!      band extension of a row where a gap opens, is at most a cell of its row minus
//!      `gap_open_extend`; the band extension of a row without one starts from the cell before
//!      `b_size` of the row before it, which is the `MININT` cell that row wrote (or there is no
//!      extension, `b_size > N`), except after row 0, which writes none: row 1's extension is at
//!      most row 0's last cell (at most 0) plus an entry. So the row computes nothing above
//!      `max(best, max) + max` (a diagonal adds an entry to a cell, a column gap is a cell minus
//!      `gap_open_extend`), and the row gap's lanes add at most `16 * gap_extend` before they
//!      subtract it again. A row starts only while `best` is at most `limit16` =
//!      `i16::MAX - max - 16 * gap_extend - 1` (and `2 * max + 16 * gap_extend <= 12,288` under
//!      `fits16`), so no lane saturates upward. When
//!      `best` is above it, the sixteen-bit rows stop and the 32-bit rows compute the alignment
//!      again from row 0 (the sixteen-bit rows wrote only their own arrays and `col`, which the
//!      32-bit rows fill again from cell 0).
//!   2. Downward: every value the 32-bit rows store is `MININT` or at least
//!      `-(x_dropoff + gap_open_extend) >= -8,192` (a kept cell is at least `best - x_dropoff`,
//!      a cell of row 0 or of the band extension at least `-x_dropoff`, a column gap is written as
//!      at least a kept or extension cell minus `gap_open_extend`), so it is stored as is. A sum
//!      whose true value is below `i16::MIN` saturates there; a lane derived from `i16::MIN` by
//!      adding an entry and a gap ramp stays at most `i16::MIN + 8,192`, below
//!      `-x_dropoff <= best_score - x_dropoff` = `thr`, as its true value is. A value below `thr`
//!      changes nothing (the argument above): it loses every `max` against a kept cell's value,
//!      its cell is dropped and stored as `i16::MIN`, a column gap is written only for kept cells,
//!      whose `score - gap_open_extend >= -8,192` wins the `max`, and a row gap whose true value
//!      is below `i16::MIN` leaves the row at most `i16::MIN + 8,192` (`MININT` for `i16::MIN`),
//!      below `thr`, so it starts no band extension, as the 32-bit value does not.
//!
//!   Neither item depends on the number of cells per step: the sixteen-bit rows in 128 bits
//!   (`LOSAT_X_BPGAP128`, `restricted_rows16_w128`) compute the same rows eight cells per step
//!   (their row gap's lanes add at most `8 * gap_extend`, within the `16 * gap_extend` that
//!   `limit16` leaves), so they store the same values and return the same result; `x_w128_16`
//!   lists what their steps do differently from the AVX2 rows'. They are the only rows of that
//!   build: where they do not take the problem (`fits16` fails, or `best` passes `limit16`) the
//!   reference function runs, as on a CPU without any vector rows (no supported matrix and gap
//!   cost pair fails `fits16`, and few restricted alignments score above 32,000).
//! - The scalar rows (`row_nogap_scalar`, `row_gap_scalar`) are the reference's loops on these
//!   arrays, cell by cell. They were slower than the reference function, so a CPU without vector
//!   rows runs the reference function; the tests run them to check the arrays apart from the
//!   vector code.
//!
//! The variables are read once. `LOSAT_X_BPGAP` (a strict-set switch, `utils/xswitch.rs`: on
//! unless it is `0`) runs this module, and
//! `LOSAT_X_BPGAP16` (another strict-set switch) its sixteen-bit rows first where they fit;
//! `LOSAT_X_BPGAP128` (another, of class platform) runs the 128-bit rows where the AVX2 rows do not
//! run (`2`: in place of them, a diagnostic).
//! With `LOSAT_X_STATS`, `main` prints `[X_STATS] BPGAP fast= not_taken=`: the restricted
//! alignments this module computed and those it passed on.

#![allow(clippy::too_many_arguments)]

use super::{
    blastp_query_residue, blastp_score_from_row, blastp_score_row_for_query,
    blastp_subject_residue, restricted_gapped_align_protein, BlastpScoreMatrix, GapAlignScratch,
    GAP_MININT, RESTRICT_SIZE,
};
use crate::config::ScoringMatrix;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1096-1097
// The reference's minus infinity (`GAP_MININT`, as in `restricted_gapped_align_protein_impl`).
const MININT: i32 = GAP_MININT;
/// No NCBI counterpart: cells kept before cell 0 of the two rows (the eight-cell block of the first
/// cell reads the cell left of it) and after the last one (a block may run up to seven cells past
/// `b_size`); they hold no value NCBI computes.
const PAD: usize = 8;
const TAIL: usize = 24;
/// No NCBI counterpart: the cells past `b_size` that a row touches in arrays holding `cells`
/// cells (LOSAT_X_BPGAPBAND): the column letters are copied up to `need + 56` with `need =
/// b_size + 17` (the sixteen-bit rows; `+ 9` in the 32-bit rows), and `col` keeps `2 + TAIL`
/// cells past its last cell.
const ROW_AHEAD: usize = 17 + 56 - (2 + TAIL);
/// No NCBI counterpart: the largest gap cost and X-drop this module takes, so that no sum below
/// overflows 32 bits (otherwise the reference function runs).
const MAX_GAP_COST: i32 = 1 << 16;
const MAX_X_DROP: i32 = 1 << 28;

// No NCBI counterpart: reads the strict-set switch LOSAT_X_BPGAP once (on: this module runs, the
// reference function when it does not take the problem); it does not change any value NCBI
// computes.
pub(super) fn x_bpgap() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPGAP"))
}

// No NCBI counterpart: reads the strict-set switch LOSAT_X_BPGAP16 once (the 16-bit rows of this
// module); it does not change any value NCBI computes.
pub(super) fn x_bpgap16() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPGAP16"))
}

// No NCBI counterpart: reads the strict-set switch LOSAT_X_BPGAPBAND once (on: the arrays of
// `XRestrictedScratch` hold the cells the band reaches, as NCBI's `dp_mem`; off: the subject
// letters `len2`); it does not change any value NCBI computes.
pub(super) fn x_bpgapband() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPGAPBAND"))
}

// No NCBI counterpart: counters of LOSAT_X_STATS; they do not change any value NCBI computes.
/// Restricted alignments computed by this module.
static FAST_CALLS: AtomicU64 = AtomicU64::new(0);
/// Restricted alignments this module did not take (the reference ran alone).
static FALLBACK_CALLS: AtomicU64 = AtomicU64::new(0);

/// The summary line `main` prints at exit with LOSAT_X_STATS: the number of alignments this module
/// computed and passed on.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    let fast = FAST_CALLS.load(Ordering::Relaxed);
    let fallback = FALLBACK_CALLS.load(Ordering::Relaxed);
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!("[X_STATS] BPGAP fast={fast} not_taken={fallback}");
    }
}

/// Rows and column letters of one restricted alignment (no NCBI counterpart: NCBI keeps them in
/// `gap_align->dp_mem`; these arrays hold the same values, see the module comment).
pub(super) struct XRestrictedScratch {
    /// `best` of the previous and the current row (`PAD` cells before cell 0).
    h: [Vec<i32>; 2],
    /// `best_gap` of the cells `10 k`.
    e10: Vec<i32>,
    /// `col[c]`: the subject letter whose score cell `c` adds (`c >= 1`).
    col: Vec<u8>,
    /// LOSAT_X_BPGAP16: `best` of the previous and the current row in 16 bits (`PAD` cells before
    /// cell 0; `i16::MIN` stands for `MININT`).
    #[cfg_attr(
        all(
            not(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            )),
            not(test)
        ),
        expect(
            dead_code,
            reason = "read only by the 16-bit rows of x86_64 and wasm32 simd128 (and a unit test)"
        )
    )]
    h16: [Vec<i16>; 2],
    /// LOSAT_X_BPGAP16: `best_gap` per column in 16 bits; `i16::MIN` in every column that is not
    /// a multiple of `RESTRICT_SIZE` (never written there).
    #[cfg_attr(
        not(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        )),
        expect(
            dead_code,
            reason = "read only by the 16-bit rows of x86_64 and wasm32 simd128"
        )
    )]
    e16: Vec<i16>,
}

impl XRestrictedScratch {
    // No NCBI counterpart: empty until the first alignment that runs here.
    pub(super) fn new() -> Self {
        Self {
            h: [Vec::new(), Vec::new()],
            e10: Vec::new(),
            col: Vec::new(),
            h16: [Vec::new(), Vec::new()],
            e16: Vec::new(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1078-1083,1115
/// `matrix[q][s]` for the letters `q, s < 28`, read through the reference's own functions; for
/// the vector kernel each row also as four 16-byte tables (low and high bytes of entries 0-15,
/// then of entries 16-31).
struct XScoreTables {
    rows: [[i32; 32]; 28],
    lut: [[u8; 64]; 28],
    /// The largest entry (LOSAT_X_BPGAP16's overflow bound).
    max: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1078-1083
// Built once per matrix. None: an entry does not fit 16 bits (the reference runs).
fn x_score_tables(score_matrix: BlastpScoreMatrix<'_>) -> Option<&'static XScoreTables> {
    fn build<const BLOSUM62: bool>(
        score_matrix: BlastpScoreMatrix<'_>,
    ) -> Option<Box<XScoreTables>> {
        let mut t = Box::new(XScoreTables {
            rows: [[0; 32]; 28],
            lut: [[0; 64]; 28],
            max: i32::MIN,
        });
        for q in 0..28u8 {
            let matrix_row = blastp_score_row_for_query::<BLOSUM62>(score_matrix, q);
            for s in 0..28u8 {
                let v = blastp_score_from_row::<BLOSUM62>(score_matrix, matrix_row, q, s);
                let v16 = i16::try_from(v).ok()?;
                t.rows[q as usize][s as usize] = v;
                t.max = t.max.max(v);
                let [lo, hi] = v16.to_le_bytes();
                let half = (s as usize / 16) * 32;
                t.lut[q as usize][half + s as usize % 16] = lo;
                t.lut[q as usize][half + 16 + s as usize % 16] = hi;
            }
        }
        Some(t)
    }
    static TABLES: [OnceLock<Option<Box<XScoreTables>>>; 10] = [const { OnceLock::new() }; 10];
    let slot = match score_matrix {
        BlastpScoreMatrix::Blosum62 => 0,
        BlastpScoreMatrix::Standard(m) => match m {
            ScoringMatrix::Blosum45 => 1,
            ScoringMatrix::Blosum50 => 2,
            ScoringMatrix::Blosum62 => 3,
            ScoringMatrix::Blosum80 => 4,
            ScoringMatrix::Blosum90 => 5,
            ScoringMatrix::Pam30 => 6,
            ScoringMatrix::Pam70 => 7,
            ScoringMatrix::Pam250 => 8,
            ScoringMatrix::Identity => 9,
        },
        // A composition-adjusted matrix belongs to one subject; the preliminary stage does not
        // use one (`BlastpScoreMatrix::standard`).
        BlastpScoreMatrix::Adjusted(_) => return None,
    };
    TABLES[slot]
        .get_or_init(|| {
            if score_matrix.is_blosum62() {
                build::<true>(score_matrix)
            } else {
                build::<false>(score_matrix)
            }
        })
        .as_deref()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4273-4278,4307-4316
/// Dispatch point of the two `s_RestrictedGappedAlign` calls of `s_BlastProtGappedAlignment`: the
/// reference function or this module. The arguments are those of
/// `restricted_gapped_align_protein`.
pub(super) fn restricted_gapped_align_protein_x(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    scratch: &mut GapAlignScratch,
) -> (usize, usize, i32) {
    if x_bpgap() {
        // Without vector rows (`x_rows`) the reference function runs: the scalar rows on the same
        // arrays were slower than it (they serve the tests as the cell-by-cell version of the
        // arrays).
        let fast = match x_rows() {
            // SAFETY: `x_rows` chooses only rows the CPU runs.
            Some(rows) => unsafe {
                restricted_fast(
                    query,
                    query_base,
                    subject,
                    subject_base,
                    len1,
                    len2,
                    score_matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    reverse,
                    &mut scratch.x_restricted,
                    rows,
                    x_bpgap16(),
                    x_bpgapband(),
                )
            },
            None => None,
        };
        crate::utils::xreach::x_reach!("LOSAT_X_BPGAP", fast.is_some());
        if let Some(fast) = fast {
            FAST_CALLS.fetch_add(1, Ordering::Relaxed);
            return fast;
        }
        FALLBACK_CALLS.fetch_add(1, Ordering::Relaxed);
    }
    restricted_gapped_align_protein(
        query,
        query_base,
        subject,
        subject_base,
        len1,
        len2,
        score_matrix,
        gap_open,
        gap_extend,
        x_drop,
        reverse,
        &mut scratch.dp_mem,
        &mut scratch.dp_mem_alloc,
    )
}

// No NCBI counterpart: the CPU check of the AVX2 rows; it does not change any value NCBI
// computes.
#[cfg(target_arch = "x86_64")]
fn cpu_has_avx2() -> bool {
    static AVX2: OnceLock<bool> = OnceLock::new();
    *AVX2.get_or_init(|| {
        is_x86_feature_detected!("avx2")
            && is_x86_feature_detected!("bmi1")
            && is_x86_feature_detected!("lzcnt")
    })
}

// No NCBI counterpart: the CPU check of the 128-bit rows (LOSAT_X_BPGAP128): SSSE3 and SSE4.1 on
// x86_64, always on wasm32 simd128; it does not change any value NCBI computes.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn cpu_has_w128() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        static W128: OnceLock<bool> = OnceLock::new();
        *W128
            .get_or_init(|| is_x86_feature_detected!("ssse3") && is_x86_feature_detected!("sse4.1"))
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        true
    }
}

// No NCBI counterpart: reads LOSAT_X_BPGAP128 once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPGAP128): 0 = off, 1 = the 128-bit sixteen-bit rows where the AVX2 rows do
/// not run (wasm32 simd128; x86_64 with SSSE3 and SSE4.1 but not AVX2, BMI1 and LZCNT), 2 = the
/// 128-bit rows on any x86_64 CPU with SSSE3 and SSE4.1, AVX2 or not (diagnostic: the tests run
/// them end to end). A switch of the strict set, on unless it is `0`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_bpgap128() -> usize {
    static ON: OnceLock<usize> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BPGAP128", 1..=2, 1).unwrap_or(0))
}

/// The rows `restricted_fast` computes with (no NCBI counterpart: every choice gives the result
/// of `s_RestrictedGappedAlign`, module comment).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum XRows {
    /// The scalar rows (`row_nogap_scalar`, `row_gap_scalar`): the tests' cell-by-cell form of the
    /// arrays.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "slower than the reference function, the scalar rows run in the tests only"
        )
    )]
    Scalar,
    /// The AVX2 rows (x86_64 with AVX2, BMI1 and LZCNT): the sixteen-bit rows first where they fit
    /// (LOSAT_X_BPGAP16), then the 32-bit rows.
    #[cfg(target_arch = "x86_64")]
    Avx2,
    /// LOSAT_X_BPGAP128: the sixteen-bit rows in 128 bits (wasm32 simd128; x86_64 with SSSE3 and
    /// SSE4.1), the only rows of this kind (module comment, "Sixteen-bit rows").
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    W128,
}

// No NCBI counterpart: chooses LOSAT's rows once; every choice gives the reference's result.
/// The rows this module runs here: the AVX2 rows where the CPU has AVX2, BMI1 and LZCNT (unless
/// LOSAT_X_BPGAP128 is `2`), else the 128-bit rows where LOSAT_X_BPGAP128 is on and the CPU has
/// SSSE3 and SSE4.1 (always on wasm32 simd128). None: the reference function runs.
fn x_rows() -> Option<XRows> {
    static ROWS: OnceLock<Option<XRows>> = OnceLock::new();
    *ROWS.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        {
            let w128 = x_bpgap128();
            if w128 != 2 && cpu_has_avx2() {
                return Some(XRows::Avx2);
            }
            if w128 != 0 && cpu_has_w128() {
                return Some(XRows::W128);
            }
            None
        }
        #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
        {
            (x_bpgap128() != 0 && cpu_has_w128()).then_some(XRows::W128)
        }
        #[cfg(not(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        )))]
        {
            None
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1096-1098,1250-1251
/// The band and the best cell, the variables of `s_RestrictedGappedAlign` that live across rows.
struct XBand {
    first: usize,
    b_size: usize,
    last: usize,
    best: i32,
    a_offset: usize,
    b_offset: usize,
}

// No NCBI counterpart: why LOSAT's rows stopped before the last row; it holds no value NCBI
// computes.
/// Why the rows of `restricted_fast` stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum XStop {
    /// LOSAT_X_BPGAPBAND: the rows reached the letters `col` holds (`fill_col`); on doubled
    /// arrays they run again from row 0.
    Band,
    /// The rows do not take the problem: a letter code of 28 or more, or (the sixteen-bit rows)
    /// `best` passed `limit16`. Larger arrays change neither: the rows run the same up to there.
    Decline,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1028-1098
/// `s_RestrictedGappedAlign` on the arrays of `sc`. `rows` selects the rows (the AVX2 rows only
/// when the CPU has AVX2, BMI1 and LZCNT, the 128-bit rows only when it has SSSE3 and SSE4.1 or
/// on wasm32 simd128); with the AVX2 rows, `narrow` (LOSAT_X_BPGAP16) runs the sixteen-bit rows
/// first where `fits16` holds, and the 32-bit rows when they stop. The 128-bit rows
/// (LOSAT_X_BPGAP128) are sixteen-bit rows only: without `narrow`, where `fits16` fails or where
/// they stop, the problem is not taken. `band` (LOSAT_X_BPGAPBAND) starts the
/// arrays at the cells the first rows can reach; where the rows reach the letters `col` holds
/// (`band_col`) they stop, the arrays double and the rows run again from row 0. Without it they
/// hold `len2` cells. None: the problem is not taken (module comment); nothing outside `sc` has
/// been changed.
///
/// # Safety
/// The CPU runs `rows`: `XRows::Avx2` only where `cpu_has_avx2` holds, `XRows::W128` only where
/// `cpu_has_w128` holds (`x_rows`, and in the tests `kernels`, choose them so).
unsafe fn restricted_fast(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    sc: &mut XRestrictedScratch,
    rows: XRows,
    narrow: bool,
    band: bool,
) -> Option<(usize, usize, i32)> {
    if !(0..=MAX_GAP_COST).contains(&gap_open)
        || !(0..=MAX_GAP_COST).contains(&gap_extend)
        || x_drop > MAX_X_DROP
    {
        return None;
    }
    let gap_open_extend = gap_open + gap_extend;
    let x_dropoff = x_drop.max(gap_open_extend);
    if len1 == 0 || len2 == 0 {
        return Some((0, 0, 0));
    }
    let tables = x_score_tables(score_matrix)?;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1040-1051
    // LOSAT_X_BPGAPBAND: `ahead` bounds the cells past a row's `b_size` that the row and its band
    // extension touch: the row reads `ROW_AHEAD` cells past `b_size`; the extension starts at
    // `b_size` and adds at most `(x_dropoff + 2 * max) / gap_extend + 1` cells and the `MININT`
    // cell (it starts from a score at most `best_score + 2 * max`, module comment: a cell of the
    // previous row is at most `max(best, max)`, plus an entry). The arrays start at `ahead + 100`
    // cells (NCBI's `num_extra_cells + 100`); `col` lets the rows copy `ahead` letters fewer than
    // the cells held (`band_col`), and where they stop the arrays double, up to `full` cells (the
    // rows then copy every letter, `len2 + 2 + TAIL`), and the rows run again from row 0. Without
    // a gap extension cost the band can run to `len2`, as NCBI's `N + 3`: the arrays hold `len2`.
    let banded = band && gap_extend > 0;
    let ahead = if banded {
        ((x_dropoff + 2 * tables.max.max(0)) / gap_extend) as usize + 3 + ROW_AHEAD
    } else {
        0
    };
    let full = len2 + ahead + 2 + TAIL;
    let first = if banded { full.min(ahead + 100) } else { len2 };
    let cells = PAD + first + 2 + TAIL;
    let p = XParams {
        x_dropoff,
        gap_open_extend,
        gap_extend,
    };
    // LOSAT_X_BPGAP128: the sixteen-bit rows in 128 bits, as the AVX2 sixteen-bit rows below
    // (the same arrays, sizes and band growth); they use no 32-bit array. Where they do not take
    // the problem, None: the reference function runs (no 32-bit rows in 128 bits, module comment).
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if rows == XRows::W128 {
        if !narrow {
            return None;
        }
        let narrow_result = if fits16(&p, tables) {
            for h in sc.h16.iter_mut() {
                if h.len() < cells {
                    h.resize(cells, i16::MIN);
                }
            }
            if sc.e16.len() < cells {
                sc.e16.resize(cells, i16::MIN);
            }
            loop {
                let [row0, row1] = &sc.h16;
                let held = band_cells16(row0, row1, &sc.e16);
                band_col(&mut sc.col, banded, held, ahead, len2);
                // SAFETY: `rows` is `W128` only on wasm32 simd128 or when the CPU has SSSE3 and
                // SSE4.1 (`cpu_has_w128`, which the tests also call); the 16-bit arrays hold
                // `held` cells and `col` the letters for them (`band_col`); `fits16` holds.
                // A stop: as for the AVX2 sixteen-bit rows below.
                let narrow = unsafe {
                    restricted_rows16_w128(
                        query,
                        query_base,
                        subject,
                        subject_base,
                        len1,
                        len2,
                        reverse,
                        tables,
                        &p,
                        sc,
                    )
                };
                match narrow {
                    Err(XStop::Band) if banded && held < full => {}
                    narrow => break narrow.ok(),
                }
                let [row0, row1] = &mut sc.h16;
                band_grow16(row0, row1, &mut sc.e16, held, full);
            }
        } else {
            None
        };
        crate::utils::xreach::x_reach!("LOSAT_X_BPGAP16", narrow_result.is_some());
        return narrow_result;
    }
    for h in sc.h.iter_mut() {
        if h.len() < cells {
            h.resize(cells, MININT);
        }
    }
    if sc.e10.len() < first / RESTRICT_SIZE + 2 {
        sc.e10.resize(first / RESTRICT_SIZE + 2, MININT);
    }
    #[cfg(target_arch = "x86_64")]
    if rows == XRows::Avx2 && narrow {
        let narrow_result = if fits16(&p, tables) {
            for h in sc.h16.iter_mut() {
                if h.len() < cells {
                    h.resize(cells, i16::MIN);
                }
            }
            if sc.e16.len() < cells {
                sc.e16.resize(cells, i16::MIN);
            }
            loop {
                let [row0, row1] = &sc.h16;
                let held = band_cells16(row0, row1, &sc.e16);
                band_col(&mut sc.col, banded, held, ahead, len2);
                // SAFETY: as for the 32-bit rows below; the 16-bit arrays hold `held` cells and
                // `col` the letters for them (`band_col`).
                // `XStop::Decline`: a score came near the 16-bit bound, and the 32-bit rows start
                // again from row 0; `XStop::Band` (LOSAT_X_BPGAPBAND, arrays short of `len2`): the
                // rows reached the letters `col` holds, and they run again on doubled arrays.
                let narrow = unsafe {
                    restricted_rows16_avx2(
                        query,
                        query_base,
                        subject,
                        subject_base,
                        len1,
                        len2,
                        reverse,
                        tables,
                        &p,
                        sc,
                    )
                };
                match narrow {
                    Err(XStop::Band) if banded && held < full => {}
                    narrow => break narrow.ok(),
                }
                let [row0, row1] = &mut sc.h16;
                band_grow16(row0, row1, &mut sc.e16, held, full);
            }
        } else {
            None
        };
        crate::utils::xreach::x_reach!("LOSAT_X_BPGAP16", narrow_result.is_some());
        if narrow_result.is_some() {
            return narrow_result;
        }
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = narrow;
    // LOSAT_X_BPGAPBAND: the 32-bit rows; where they reach the letters `col` holds (`XStop::Band`,
    // arrays short of `len2`) they run again on doubled arrays.
    loop {
        let [row0, row1] = &sc.h;
        let held = band_cells32(row0, row1, &sc.e10);
        band_col(&mut sc.col, banded, held, ahead, len2);
        #[cfg(target_arch = "x86_64")]
        let wide = if rows == XRows::Avx2 {
            // SAFETY: `rows` is `Avx2` only when the CPU has AVX2, BMI1 and LZCNT
            // (`cpu_has_avx2`, which the tests also call). The arrays of `sc` hold `held` cells
            // and `col` the letters for them (`band_col`).
            unsafe {
                restricted_rows_avx2(
                    query,
                    query_base,
                    subject,
                    subject_base,
                    len1,
                    len2,
                    reverse,
                    tables,
                    &p,
                    sc,
                )
            }
        } else {
            restricted_rows::<false, COL32>(
                query,
                query_base,
                subject,
                subject_base,
                len1,
                len2,
                reverse,
                tables,
                &p,
                sc,
            )
        };
        #[cfg(not(target_arch = "x86_64"))]
        let wide = {
            let _ = rows;
            restricted_rows::<false, COL32>(
                query,
                query_base,
                subject,
                subject_base,
                len1,
                len2,
                reverse,
                tables,
                &p,
                sc,
            )
        };
        match wide {
            Err(XStop::Band) if banded && held < full => {}
            wide => return wide.ok(),
        }
        let [row0, row1] = &mut sc.h;
        band_grow32(row0, row1, &mut sc.e10, held, full);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1076
/// The rows with the vector kernel, compiled for AVX2 as a whole (no call per row).
///
/// # Safety
/// AVX2, BMI1 and LZCNT; `sc` sized by `restricted_fast`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,bmi1,lzcnt")]
unsafe fn restricted_rows_avx2(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    reverse: bool,
    tables: &XScoreTables,
    p: &XParams,
    sc: &mut XRestrictedScratch,
) -> Result<(usize, usize, i32), XStop> {
    restricted_rows::<true, COL32_AVX2>(
        query,
        query_base,
        subject,
        subject_base,
        len1,
        len2,
        reverse,
        tables,
        p,
        sc,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1053-1098
/// Row 0 and the row loop of `s_RestrictedGappedAlign` (`AVX2`: the vector rows; the caller
/// guarantees the CPU features and has sized the arrays of `sc` for the cells the band reaches:
/// `len2`, or with LOSAT_X_BPGAPBAND the letters `col` holds, `fill_col`). `COL`: the instance of
/// `fill_col` of these rows (`COL32_AVX2` for the vector rows, `COL32` for the scalar ones).
/// `XStop::Decline`: a letter code of 28 or more; `XStop::Band`: the rows reached the letters
/// `col` holds. Nothing outside `sc` has been changed.
#[inline(always)]
fn restricted_rows<const AVX2: bool, const COL: u8>(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    reverse: bool,
    tables: &XScoreTables,
    p: &XParams,
    sc: &mut XRestrictedScratch,
) -> Result<(usize, usize, i32), XStop> {
    let x_dropoff = p.x_dropoff;
    let gap_open_extend = p.gap_open_extend;
    let gap_extend = p.gap_extend;
    let col_cap = len2 + 2 + TAIL;
    let mut col_ready = 0usize;

    // Cell `b` of a row is `h[PAD + b]`. Row 0 goes to `h[0]`.
    let [row0, row1] = &mut sc.h;
    let e10 = &mut sc.e10[..];
    let col = &mut sc.col[..];
    let mut score = -gap_open_extend;
    row0[PAD] = 0;
    e10[0] = -gap_open_extend;
    let mut b_size = 1usize;
    for i in 1..=len2 {
        if score < -x_dropoff {
            break;
        }
        row0[PAD + i] = score;
        if i % RESTRICT_SIZE == 0 {
            e10[i / RESTRICT_SIZE] = score - gap_open_extend;
        }
        score -= gap_extend;
        b_size = i + 1;
    }

    let mut band = XBand {
        first: 0,
        b_size,
        last: 0,
        best: 0,
        a_offset: 0,
        b_offset: 0,
    };

    for a_index in 1..=len1 {
        let qc = if reverse {
            blastp_query_residue::<true>(query, query_base, len1, a_index)
        } else {
            blastp_query_residue::<false>(query, query_base, len1, a_index)
        };
        if qc >= 28 {
            return Err(XStop::Decline);
        }
        // Letters for the cells up to `b_size` and the block that may run past it. A letter is
        // below 28 (checked when it is copied), so `& 31` in the scalar rows changes no index.
        let need = (band.b_size + 9).min(col_cap);
        if col_ready < need {
            let upto = (need + 56).min(col_cap);
            fill_col::<COL>(col, col_ready, upto, subject, subject_base, len2, reverse)?;
            col_ready = upto;
        }
        let (prev, cur): (&[i32], &mut [i32]) = if a_index % 2 == 1 {
            (&row0[..], &mut row1[..])
        } else {
            (&row1[..], &mut row0[..])
        };
        let gap_row = a_index % RESTRICT_SIZE == 0;
        let mut score_gap_row = run_row::<AVX2>(
            gap_row,
            &mut band,
            p,
            a_index,
            &tables.rows[qc as usize],
            &tables.lut[qc as usize],
            prev,
            cur,
            col,
            e10,
        );

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1250-1251,1271-1296
        // The next `b_gap` (lines 1256-1259) is the first multiple of RESTRICT_SIZE at or after
        // `first_b_index`, which the rows find from `first_b_index` itself.
        if band.first == band.b_size {
            break;
        }
        if band.last < band.b_size.saturating_sub(1) {
            band.b_size = band.last + 1;
        } else {
            while score_gap_row >= band.best - x_dropoff && band.b_size <= len2 {
                cur[PAD + band.b_size] = score_gap_row;
                if band.b_size.is_multiple_of(RESTRICT_SIZE) {
                    e10[band.b_size / RESTRICT_SIZE] = score_gap_row - gap_open_extend;
                }
                score_gap_row -= gap_extend;
                band.b_size += 1;
            }
        }
        if band.b_size <= len2 {
            cur[PAD + band.b_size] = MININT;
            if band.b_size.is_multiple_of(RESTRICT_SIZE) {
                e10[band.b_size / RESTRICT_SIZE] = MININT;
            }
            band.b_size += 1;
        }
    }

    Ok((band.a_offset, band.b_offset, band.best))
}

/// `gap_open_extend`, `gap_extend` and `x_dropoff` of `s_RestrictedGappedAlign`.
struct XParams {
    x_dropoff: i32,
    gap_open_extend: i32,
    gap_extend: i32,
}

// No NCBI counterpart (LOSAT_X_BPGAPBAND): the cells the 32-bit arrays hold, by the sizes
// `restricted_fast` gives them for `cells`: `PAD + cells + 2 + TAIL` rows and `cells /
// RESTRICT_SIZE + 2` column gaps.
fn band_cells32(row0: &[i32], row1: &[i32], e10: &[i32]) -> usize {
    let rows = row0.len().min(row1.len()).saturating_sub(PAD + 2 + TAIL);
    let gaps = match e10.len() {
        0 | 1 => 0,
        n => (n - 2) * RESTRICT_SIZE + RESTRICT_SIZE - 1,
    };
    rows.min(gaps)
}

// No NCBI counterpart (LOSAT_X_BPGAPBAND): as `band_cells32` for the sixteen-bit arrays
// (`PAD + cells + 2 + TAIL` rows and column gaps).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn band_cells16(row0: &[i16], row1: &[i16], e16: &[i16]) -> usize {
    row0.len()
        .min(row1.len())
        .min(e16.len())
        .saturating_sub(PAD + 2 + TAIL)
}

// No NCBI counterpart (LOSAT_X_BPGAPBAND): `col` given the length that bounds the letters the
// rows copy (`fill_col` stops past `col.len() - 2 - TAIL`). `banded`: arrays holding `held` cells
// get `held - ahead` letters, since a row and its band extension touch at most `ahead` cells past
// the row's `b_size` and the rows copy letters past `b_size` before they use them; the length
// depends on the arrays only, so it changes only when they grow. Otherwise (arrays of `len2`
// cells) `col` holds at least `2 + TAIL` more than the rows' `len2 + 2 + TAIL` letters and the
// rows never stop. Only the length changes; the vector keeps its capacity.
fn band_col(col: &mut Vec<u8>, banded: bool, held: usize, ahead: usize, len2: usize) {
    if banded {
        let n = held.saturating_sub(ahead) + 2 + TAIL;
        if col.len() != n {
            col.resize(n, 0);
        }
    } else if col.len() < len2 + 2 * (2 + TAIL) {
        col.resize(len2 + 2 * (2 + TAIL), 0);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1261-1269
// LOSAT_X_BPGAPBAND: the 32-bit arrays grown to twice the `held` cells (at most `full`, where the
// rows copy every letter), as NCBI's `2 * dp_mem_alloc` (`held` is at least `ahead + 100`, so it
// is also past NCBI's first term).
// NCBI grows the block within the row loop; here the rows stopped where the letters ran out and
// run again from row 0. NCBI's block is written only where the band reaches; the vectors write
// the cells they hold.
#[cold]
#[inline(never)]
fn band_grow32(
    row0: &mut Vec<i32>,
    row1: &mut Vec<i32>,
    e10: &mut Vec<i32>,
    held: usize,
    full: usize,
) {
    let cells = full.min(2 * held.max(1));
    for h in [&mut *row0, &mut *row1] {
        if h.len() < PAD + cells + 2 + TAIL {
            h.resize(PAD + cells + 2 + TAIL, MININT);
        }
    }
    if e10.len() < cells / RESTRICT_SIZE + 2 {
        e10.resize(cells / RESTRICT_SIZE + 2, MININT);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1261-1269 (as `band_grow32`)
// LOSAT_X_BPGAPBAND: `band_grow32` for the sixteen-bit arrays.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cold]
#[inline(never)]
fn band_grow16(
    row0: &mut Vec<i16>,
    row1: &mut Vec<i16>,
    e16: &mut Vec<i16>,
    held: usize,
    full: usize,
) {
    let cells = full.min(2 * held.max(1));
    for h in [&mut *row0, &mut *row1, &mut *e16] {
        if h.len() < PAD + cells + 2 + TAIL {
            h.resize(PAD + cells + 2 + TAIL, i16::MIN);
        }
    }
}

// No NCBI counterpart: the instances of `fill_col`, one per row kernel; they hold no value NCBI
// computes.
const COL32: u8 = 0;
#[cfg(target_arch = "x86_64")]
const COL32_AVX2: u8 = 3;
#[cfg(target_arch = "x86_64")]
const COL16: u8 = 1;
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const COL16_W128: u8 = 2;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1091-1094,1114-1115
/// `col[c]` for `c` in `from..to`: the letter the reference reads for `next_score` at
/// `b_index = c - 1` (`blastp_subject_residue`, sentinels included); 0 for cell 0 and past cell
/// `len2 + 1`, which no cell of the reference reads. The cells whose letter lies inside `subject`
/// are one slice copy (forward: `subject[subject_base + c]`; reverse: `subject[subject_base + len2 -
/// c]`, the indices that function computes); the others go through that function.
/// `XStop::Decline`: a letter code of 28 or more (not taken); `XStop::Band`: `to` is past the
/// letters `col` holds.
///
/// `ROWS`: one instance per row kernel (`COL32_AVX2`: the AVX2 32-bit rows, `COL32`: the scalar
/// 32-bit rows, `COL16`: the AVX2 sixteen-bit rows, `COL16_W128`: the 128-bit sixteen-bit rows),
/// each with one caller, so that fat LTO inlines each into its rows. Shared by the 32-bit and the
/// sixteen-bit rows, it was kept out of line, which changed the machine code of the 32-bit rows
/// that `LOSAT_X_BPGAP16=0` runs; shared by the AVX2 and the scalar 32-bit rows, it was kept out of
/// line once `restricted_fast` took the 128-bit rows. With one instance each, the AVX2 rows keep
/// the machine code of the build before LOSAT_X_BPGAP128.
fn fill_col<const ROWS: u8>(
    col: &mut [u8],
    from: usize,
    to: usize,
    subject: &[u8],
    subject_base: usize,
    len2: usize,
    reverse: bool,
) -> Result<(), XStop> {
    // LOSAT_X_BPGAPBAND: `col`'s length bounds the letters copied (`band_col`): past it the rows
    // stop (`XStop::Band`), before a row's band passes the cells the arrays hold.
    if to + 2 + TAIL > col.len() {
        return Err(XStop::Band);
    }
    let reference = |c: usize| -> u8 {
        if c == 0 || c > len2 + 1 {
            0
        } else if reverse {
            blastp_subject_residue::<true>(subject, subject_base, len2, c - 1)
        } else {
            blastp_subject_residue::<false>(subject, subject_base, len2, c - 1)
        }
    };
    // Cells `fa..fb` read inside `subject`.
    let (fa, fb) = if reverse {
        (
            (subject_base + len2 + 1)
                .saturating_sub(subject.len())
                .max(1),
            len2 + 1,
        )
    } else {
        (
            1,
            (len2 + 2).min(subject.len().saturating_sub(subject_base)),
        )
    };
    let a = fa.max(from);
    let b = fb.min(to);
    if a < b {
        if reverse {
            let src = &subject[subject_base + len2 + 1 - b..=subject_base + len2 - a];
            for (slot, &v) in col[a..b].iter_mut().zip(src.iter().rev()) {
                *slot = v;
            }
        } else {
            col[a..b].copy_from_slice(&subject[subject_base + a..subject_base + b]);
        }
        for c in (from..a).chain(b..to) {
            col[c] = reference(c);
        }
    } else {
        for c in from..to {
            col[c] = reference(c);
        }
    }
    if col[from..to].iter().all(|&v| v < 28) {
        Ok(())
    } else {
        Err(XStop::Decline)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1107,1173
/// One row: the vector kernel when `AVX2`, else the scalar loops. Returns `score_gap_row` after
/// the row.
#[inline(always)]
fn run_row<const AVX2: bool>(
    gap_row: bool,
    band: &mut XBand,
    p: &XParams,
    a_index: usize,
    row: &[i32; 32],
    lut: &[u8; 64],
    prev: &[i32],
    cur: &mut [i32],
    col: &[u8],
    e10: &mut [i32],
) -> i32 {
    #[cfg(target_arch = "x86_64")]
    if AVX2 {
        // SAFETY: `AVX2` is true only inside `restricted_rows_avx2`, whose caller checked AVX2,
        // BMI1 and LZCNT. `prev` and `cur` hold `PAD` cells before cell 0 and `TAIL` cells after
        // cell `cells + 1`, `col` and `e10` cover every index the row reads (sizes set in
        // `restricted_fast` for `cells` = `len2`, or with LOSAT_X_BPGAPBAND for `cells` where the
        // letters `col` holds keep `b_size + ahead` within `cells`, `band_col`; `b_size <= cells +
        // 1`).
        // SAFETY: `x_avx2::row`'s contract: AVX2, BMI1 and LZCNT, and `prev`, `cur`, `col` and
        // `e10` long enough for `band.b_size` (paragraph above; `col`: `fill_col` succeeded up to
        // `b_size + 9`, which it does only when `to + 2 + TAIL <= col.len()`).
        unsafe {
            return if gap_row {
                x_avx2::row::<true>(band, p, a_index, row, lut, prev, cur, col, e10)
            } else {
                x_avx2::row::<false>(band, p, a_index, row, lut, prev, cur, col, e10)
            };
        }
    }
    let _ = lut;
    if gap_row {
        row_gap_scalar(band, p, a_index, row, prev, cur, col, e10)
    } else {
        row_nogap_scalar(band, p, a_index, row, prev, cur, col, e10)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1112-1171
/// A row where no gap opens in A, cell by cell (`b_gap` steps over the multiples of 10 from the
/// first one at or after `first_b_index`, as in the reference; `e10[gk]` is its column gap).
fn row_nogap_scalar(
    band: &mut XBand,
    p: &XParams,
    a_index: usize,
    row: &[i32; 32],
    prev: &[i32],
    cur: &mut [i32],
    col: &[u8],
    e10: &mut [i32],
) -> i32 {
    // The band in locals while the row runs.
    let b_size = band.b_size;
    let mut first = band.first;
    let mut last = first;
    let mut best = band.best;
    let mut a_offset = band.a_offset;
    let mut b_offset = band.b_offset;
    let mut score = MININT;
    let mut gk = first.div_ceil(RESTRICT_SIZE);
    let mut b_gap = gk * RESTRICT_SIZE;
    #[expect(
        clippy::mut_range_bound,
        reason = "NCBI's loop reads its start once and moves first_b_index past a first column that fails X-drop (blast_gapalign.c:1112,1124-1125)"
    )]
    for b_index in first..b_size {
        let next_score = prev[PAD + b_index] + row[(col[b_index + 1] & 31) as usize];
        if b_index != b_gap {
            if best - score > p.x_dropoff {
                cur[PAD + b_index] = MININT;
                if b_index == first {
                    first += 1;
                }
            } else {
                last = b_index;
                if score > best {
                    best = score;
                    a_offset = a_index;
                    b_offset = b_index;
                }
                cur[PAD + b_index] = score;
            }
        } else {
            b_gap += RESTRICT_SIZE;
            let mut score_gap_col = e10[gk];
            if score < score_gap_col {
                score = score_gap_col;
            }
            if best - score > p.x_dropoff {
                cur[PAD + b_index] = MININT;
                if b_index == first {
                    first += 1;
                }
            } else {
                last = b_index;
                if score > best {
                    best = score;
                    a_offset = a_index;
                    b_offset = b_index;
                }
                score_gap_col -= p.gap_extend;
                e10[gk] = (score - p.gap_open_extend).max(score_gap_col);
                cur[PAD + b_index] = score;
            }
            gk += 1;
        }
        score = next_score;
    }
    band.first = first;
    band.last = last;
    band.best = best;
    band.a_offset = a_offset;
    band.b_offset = b_offset;
    score
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1177-1247
/// A row where a gap may open in A, cell by cell.
fn row_gap_scalar(
    band: &mut XBand,
    p: &XParams,
    a_index: usize,
    row: &[i32; 32],
    prev: &[i32],
    cur: &mut [i32],
    col: &[u8],
    e10: &mut [i32],
) -> i32 {
    // The band in locals while the row runs.
    let b_size = band.b_size;
    let mut first = band.first;
    let mut last = first;
    let mut best = band.best;
    let mut a_offset = band.a_offset;
    let mut b_offset = band.b_offset;
    let mut score = MININT;
    let mut score_gap_row = MININT;
    let mut gk = first.div_ceil(RESTRICT_SIZE);
    let mut b_gap = gk * RESTRICT_SIZE;
    #[expect(
        clippy::mut_range_bound,
        reason = "NCBI's loop reads its start once and moves first_b_index past a first column that fails X-drop (blast_gapalign.c:1177,1192-1193)"
    )]
    for b_index in first..b_size {
        let next_score = prev[PAD + b_index] + row[(col[b_index + 1] & 31) as usize];
        if b_index != b_gap {
            if score < score_gap_row {
                score = score_gap_row;
            }
            if best - score > p.x_dropoff {
                cur[PAD + b_index] = MININT;
                if b_index == first {
                    first += 1;
                }
            } else {
                last = b_index;
                if score > best {
                    best = score;
                    a_offset = a_index;
                    b_offset = b_index;
                }
                score_gap_row -= p.gap_extend;
                score_gap_row = (score - p.gap_open_extend).max(score_gap_row);
                cur[PAD + b_index] = score;
            }
        } else {
            b_gap += RESTRICT_SIZE;
            let mut score_gap_col = e10[gk];
            if score < score_gap_col {
                score = score_gap_col;
            }
            if score < score_gap_row {
                score = score_gap_row;
            }
            if best - score > p.x_dropoff {
                cur[PAD + b_index] = MININT;
                if b_index == first {
                    first += 1;
                }
            } else {
                last = b_index;
                if score > best {
                    best = score;
                    a_offset = a_index;
                    b_offset = b_index;
                }
                score_gap_row -= p.gap_extend;
                score_gap_col -= p.gap_extend;
                e10[gk] = (score - p.gap_open_extend).max(score_gap_col);
                score_gap_row = (score - p.gap_open_extend).max(score_gap_row);
                cur[PAD + b_index] = score;
            }
            gk += 1;
        }
        score = next_score;
    }
    band.first = first;
    band.last = last;
    band.best = best;
    band.a_offset = a_offset;
    band.b_offset = b_offset;
    score_gap_row
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1112-1247
// The two row loops eight cells at a time (AVX2); the module comment gives the argument.
#[cfg(target_arch = "x86_64")]
mod x_avx2 {
    use super::{XBand, XParams, MININT, PAD, RESTRICT_SIZE};
    use core::arch::x86_64::*;

    // No NCBI counterpart: lane masks; they hold no value NCBI computes.
    /// `VALID[8 - n..]`: the first `n` lanes set.
    static VALID: [i32; 16] = [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0];
    /// `ONE[7 - j..]`: lane `j` set.
    static ONE: [i32; 15] = [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0];

    /// Inclusive prefix maximum over the eight lanes (lanes left of lane 0 count as `MININT`).
    #[inline(always)]
    unsafe fn prefix_max(v: __m256i, minv: __m256i) -> __m256i {
        let sh1 = _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6);
        let sh2 = _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5);
        let t = _mm256_max_epi32(
            v,
            _mm256_blend_epi32::<0x01>(_mm256_permutevar8x32_epi32(v, sh1), minv),
        );
        let t = _mm256_max_epi32(
            t,
            _mm256_blend_epi32::<0x03>(_mm256_permutevar8x32_epi32(t, sh2), minv),
        );
        _mm256_max_epi32(t, _mm256_permute2x128_si256::<0x02>(t, minv))
    }

    /// Lane `i` = `max(carry, p[i - 1])`, lane 0 = `carry` (exclusive prefix from an inclusive one).
    #[inline(always)]
    unsafe fn exclusive(p: __m256i, carry: __m256i) -> __m256i {
        let sh1 = _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6);
        _mm256_max_epi32(
            _mm256_blend_epi32::<0x01>(_mm256_permutevar8x32_epi32(p, sh1), carry),
            carry,
        )
    }

    /// One row of `s_RestrictedGappedAlign`, eight cells per step. `GAP`: a row where a gap may
    /// open in A. Cell `b` of `prev` and `cur` is at `PAD + b`. Reads the cells `first - 1 ..
    /// b_size + 7` of `prev` and `col[first .. b_size + 8]`, writes the cells `first .. b_size + 8`
    /// of `cur` (cells from `b_size` on get `MININT`, which nothing reads before the band
    /// extension writes them) and `e10` at the multiples of 10 inside the band. Returns
    /// `score_gap_row` after the row (the diagonal into cell `b_size`, or the row gap).
    ///
    /// # Safety
    /// Inlined into code compiled for AVX2, BMI1 and LZCNT (`restricted_rows_avx2`); `prev` and
    /// `cur` at least `PAD + b_size + 8` long, `col` at least
    /// `b_size + 8`, `e10` at least `b_size / 10 + 1`.
    #[inline(always)]
    pub(super) unsafe fn row<const GAP: bool>(
        band: &mut XBand,
        p: &XParams,
        a_index: usize,
        row: &[i32; 32],
        lut: &[u8; 64],
        prev: &[i32],
        cur: &mut [i32],
        col: &[u8],
        e10: &mut [i32],
    ) -> i32 {
        let first0 = band.first;
        let end = band.b_size;
        debug_assert!(first0 < end);
        debug_assert!(prev.len() >= PAD + end + 8 && cur.len() >= PAD + end + 8);
        debug_assert!(col.len() >= end + 8);
        // Cell `b` is at `PAD + b`; the step of cell 0 reads `PAD - 1`.
        let prev_p = prev.as_ptr().add(PAD);
        let cur_p = cur.as_mut_ptr().add(PAD);
        let col_p = col.as_ptr();
        let minv = _mm256_set1_epi32(MININT);
        let xv = _mm256_set1_epi32(p.x_dropoff);
        let goev = _mm256_set1_epi32(p.gap_open_extend);
        let t0l = _mm_loadu_si128(lut.as_ptr() as *const __m128i);
        let t0h = _mm_loadu_si128(lut.as_ptr().add(16) as *const __m128i);
        let t1l = _mm_loadu_si128(lut.as_ptr().add(32) as *const __m128i);
        let t1h = _mm_loadu_si128(lut.as_ptr().add(48) as *const __m128i);
        let k15 = _mm_set1_epi8(15);
        let ge = p.gap_extend;
        // (j + 1) * gap_extend and j * gap_extend for lane j.
        let ramp1 = _mm256_setr_epi32(ge, 2 * ge, 3 * ge, 4 * ge, 5 * ge, 6 * ge, 7 * ge, 8 * ge);
        let ramp0 = _mm256_setr_epi32(0, ge, 2 * ge, 3 * ge, 4 * ge, 5 * ge, 6 * ge, 7 * ge);

        let mut best = band.best;
        let mut rv = _mm256_set1_epi32(best);
        // best_score - x_dropoff, the drop threshold while `best` holds.
        let mut thr = _mm256_sub_epi32(rv, xv);
        let mut leading = true;
        let mut first = first0;
        let mut last = first0;
        // The next column where a gap may start in B (`b_gap`).
        let mut gk = first0.div_ceil(RESTRICT_SIZE);
        let mut gcol = gk * RESTRICT_SIZE;
        // Row gap entering the block (rows where a gap opens in A).
        let mut f_in = MININT;
        let mut score_gap_row = MININT;

        let mut b0 = first0;
        while b0 < end {
            let n = (end - b0).min(8);
            // Diagonal: prev[b - 1] + matrix_row[col[b]].
            let idx = _mm_loadl_epi64(col_p.add(b0) as *const __m128i);
            let sel = _mm_cmpgt_epi8(idx, k15);
            let lo = _mm_blendv_epi8(_mm_shuffle_epi8(t0l, idx), _mm_shuffle_epi8(t1l, idx), sel);
            let hi = _mm_blendv_epi8(_mm_shuffle_epi8(t0h, idx), _mm_shuffle_epi8(t1h, idx), sel);
            let m = _mm256_cvtepi16_epi32(_mm_unpacklo_epi8(lo, hi));
            let mut h0 = _mm256_add_epi32(
                _mm256_loadu_si256(prev_p.add(b0).sub(1) as *const __m256i),
                m,
            );
            if b0 == first0 {
                // Initial score is `MININT` at the start of the row.
                h0 = _mm256_blend_epi32::<0x01>(h0, minv);
            }
            // Column gap in the column that is a multiple of 10 (at most one per block).
            let has_gap_col = gcol < b0 + n;
            if has_gap_col {
                let one = _mm256_loadu_si256(ONE.as_ptr().add(7 - (gcol - b0)) as *const __m256i);
                let ev = _mm256_set1_epi32(*e10.get_unchecked(gk));
                h0 = _mm256_max_epi32(h0, _mm256_blendv_epi8(minv, ev, one));
            }
            // Lanes past `b_size` (the last step of the row only) hold no cell.
            let valid = if n < 8 {
                Some(_mm256_loadu_si256(
                    VALID.as_ptr().add(8 - n) as *const __m256i
                ))
            } else {
                None
            };
            let mut s = h0;
            if GAP {
                // Row gap: F(i) = max(f_in, max_{j<i} (h0[j] - goe + (j+1) ge)) - i ge.
                let mut w = _mm256_add_epi32(_mm256_sub_epi32(h0, goev), ramp1);
                if let Some(valid) = valid {
                    w = _mm256_blendv_epi8(minv, w, valid);
                }
                let pw = prefix_max(w, minv);
                let fv = _mm256_set1_epi32(f_in);
                let f = _mm256_sub_epi32(exclusive(pw, fv), ramp0);
                s = _mm256_max_epi32(s, f);
                // Row gap after the last cell of the block (cell b0 + n).
                let pl = _mm256_cvtsi256_si32(_mm256_permutevar8x32_epi32(
                    pw,
                    _mm256_set1_epi32(n as i32 - 1),
                ));
                let f_out = pl.max(f_in) - n as i32 * ge;
                f_in = f_out.max(MININT);
                score_gap_row = f_in;
            }
            if let Some(valid) = valid {
                s = _mm256_blendv_epi8(minv, s, valid);
            }
            // X-drop: best_score when a cell is tested = max(best at the block, cells left of it).
            // When no cell of the block exceeds `best`, that is `best` for every cell.
            let above = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(s, rv)));
            let dropped = if above == 0 {
                _mm256_cmpgt_epi32(thr, s)
            } else {
                let ps = prefix_max(s, minv);
                let ex = exclusive(ps, rv);
                // The block's largest score exceeds `best`: `score > best_score` holds last at the
                // first cell that holds it.
                let block_max = _mm256_extract_epi32::<7>(ps);
                let eq = _mm256_cmpeq_epi32(s, _mm256_set1_epi32(block_max));
                let lane = (_mm256_movemask_ps(_mm256_castsi256_ps(eq)) as u32).trailing_zeros();
                best = block_max;
                band.a_offset = a_index;
                band.b_offset = b0 + lane as usize;
                rv = _mm256_set1_epi32(best);
                thr = _mm256_sub_epi32(rv, xv);
                _mm256_cmpgt_epi32(_mm256_sub_epi32(ex, xv), s)
            };
            _mm256_storeu_si256(
                cur_p.add(b0) as *mut __m256i,
                _mm256_blendv_epi8(s, minv, dropped),
            );
            let kept =
                !(_mm256_movemask_ps(_mm256_castsi256_ps(dropped)) as u32) & ((1u32 << n) - 1);
            if kept != 0 {
                last = b0 + 31 - kept.leading_zeros() as usize;
                if leading {
                    first = b0 + kept.trailing_zeros() as usize;
                    leading = false;
                }
            } else if leading {
                first = b0 + n;
            }
            if has_gap_col {
                // A kept cell is >= best - x_dropoff >= -x_dropoff > MININT.
                let kept_score = *cur_p.add(gcol);
                if kept_score != MININT {
                    let e = e10.get_unchecked_mut(gk);
                    *e = (kept_score - p.gap_open_extend).max(*e - ge);
                }
                gk += 1;
                gcol += RESTRICT_SIZE;
            }
            b0 += 8;
        }
        band.best = best;
        band.first = first;
        band.last = last;
        if !GAP {
            // score after the loop: next_score of cell b_size - 1.
            score_gap_row = prev[PAD + end - 1] + row[col[end] as usize];
        }
        score_gap_row
    }
}

// No NCBI counterpart: the bound under which the 16-bit rows (LOSAT_X_BPGAP16) hold the values of
// the 32-bit rows. A row computes nothing above `max(best, max) + max` (`max` the largest matrix
// entry; module comment, "Sixteen-bit rows", item 1), and the row gap's intermediate adds at most
// `16 * gap_extend` (the 128-bit rows', LOSAT_X_BPGAP128, at most `8 * gap_extend`); both stay below
// `i16::MAX` while `best` at the start of the row is at most this bound (`fits16` keeps `2 * max +
// 16 * gap_extend` at most 12,288).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn limit16(p: &XParams, tables: &XScoreTables) -> i32 {
    i16::MAX as i32 - tables.max.max(0) - 16 * p.gap_extend - 1
}

// No NCBI counterpart: the costs the 16-bit rows take. With them, a value that saturates at
// `i16::MIN` (or `i16::MIN` plus a gap ramp or an entry) stays below `best - x_dropoff >=
// -x_dropoff`, the threshold below which no value changes what the reference keeps (module
// comment). They also keep `limit16` above 24,000.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn fits16(p: &XParams, tables: &XScoreTables) -> bool {
    p.x_dropoff <= 4096 && p.gap_open_extend <= 4096 && p.gap_extend <= 256 && tables.max <= 4096
}

// No NCBI counterpart: the kernel of `restricted_rows16`; it holds no value NCBI computes.
/// A kernel of the sixteen-bit rows (`x_avx2_16::Rows`: sixteen cells per step, AVX2;
/// `x_w128_16::Rows`: eight, LOSAT_X_BPGAP128): what `restricted_rows16` runs for the letters and
/// the rows. Each kernel's function (`restricted_rows16_avx2`, `restricted_rows16_w128`) is
/// compiled for its instructions and inlines the driver with its kernel.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
trait XRows16 {
    /// The vectors the rows of one alignment use, made once (`()` where the rows make none).
    type Consts;

    /// # Safety
    /// Inlined into code compiled for the kernel's instructions; `fits16`.
    unsafe fn consts(p: &XParams) -> Self::Consts;

    /// The kernel's own instance of `fill_col` (`fill_col`, `ROWS`).
    fn fill_col(
        col: &mut [u8],
        from: usize,
        to: usize,
        subject: &[u8],
        subject_base: usize,
        len2: usize,
        reverse: bool,
    ) -> Result<(), XStop>;

    /// One row (`x_avx2_16::row`, `x_w128_16::row`). Returns `score_gap_row` after the row.
    ///
    /// # Safety
    /// As the kernel's `row`.
    unsafe fn row<const GAP: bool>(
        band: &mut XBand,
        p: &XParams,
        a_index: usize,
        row: &[i32; 32],
        lut: &[u8; 64],
        prev: &[i16],
        cur: &mut [i16],
        col: &[u8],
        e16: &mut [i16],
        c: &Self::Consts,
    ) -> i32;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1053-1067,1076-1098,1250-1287
/// LOSAT_X_BPGAP16: `restricted_rows` with the rows, the column gaps and the row gap
/// in 16 bits, one step of the kernel `K` at a time. `i16::MIN` stands for `MININT`; every other
/// stored value is the 32-bit rows' value (module comment, "Sixteen-bit rows").
/// `XStop::Decline`: `best` passed `limit16` at the start of a row, or a letter code of 28 or more
/// (the AVX2 caller then runs the 32-bit rows from row 0); `XStop::Band`: the rows reached the
/// letters `col` holds (they run again on doubled arrays). Nothing outside `sc` has been changed.
///
/// # Safety
/// Inlined into code compiled for `K`'s instructions (`restricted_rows16_avx2`,
/// `restricted_rows16_w128`); the arrays of `sc` (16-bit ones included) sized by
/// `restricted_fast`; `fits16`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
unsafe fn restricted_rows16<K: XRows16>(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    reverse: bool,
    tables: &XScoreTables,
    p: &XParams,
    sc: &mut XRestrictedScratch,
) -> Result<(usize, usize, i32), XStop> {
    let x_dropoff = p.x_dropoff;
    let gap_open_extend = p.gap_open_extend;
    let gap_extend = p.gap_extend;
    let limit = limit16(p, tables);
    let col_cap = len2 + 2 + TAIL;
    let mut col_ready = 0usize;
    let consts = K::consts(p);

    let [row0, row1] = &mut sc.h16;
    let e16 = &mut sc.e16[..];
    let col = &mut sc.col[..];
    let mut score = -gap_open_extend;
    row0[PAD] = 0;
    e16[0] = -gap_open_extend as i16;
    let mut b_size = 1usize;
    for i in 1..=len2 {
        if score < -x_dropoff {
            break;
        }
        row0[PAD + i] = score as i16;
        if i % RESTRICT_SIZE == 0 {
            e16[i] = (score - gap_open_extend) as i16;
        }
        score -= gap_extend;
        b_size = i + 1;
    }

    let mut band = XBand {
        first: 0,
        b_size,
        last: 0,
        best: 0,
        a_offset: 0,
        b_offset: 0,
    };

    for a_index in 1..=len1 {
        if band.best > limit {
            return Err(XStop::Decline);
        }
        let qc = if reverse {
            blastp_query_residue::<true>(query, query_base, len1, a_index)
        } else {
            blastp_query_residue::<false>(query, query_base, len1, a_index)
        };
        if qc >= 28 {
            return Err(XStop::Decline);
        }
        let need = (band.b_size + 17).min(col_cap);
        if col_ready < need {
            let upto = (need + 56).min(col_cap);
            K::fill_col(col, col_ready, upto, subject, subject_base, len2, reverse)?;
            col_ready = upto;
        }
        let (prev, cur): (&[i16], &mut [i16]) = if a_index % 2 == 1 {
            (&row0[..], &mut row1[..])
        } else {
            (&row1[..], &mut row0[..])
        };
        let mut score_gap_row = if a_index % RESTRICT_SIZE == 0 {
            K::row::<true>(
                &mut band,
                p,
                a_index,
                &tables.rows[qc as usize],
                &tables.lut[qc as usize],
                prev,
                cur,
                col,
                e16,
                &consts,
            )
        } else {
            K::row::<false>(
                &mut band,
                p,
                a_index,
                &tables.rows[qc as usize],
                &tables.lut[qc as usize],
                prev,
                cur,
                col,
                e16,
                &consts,
            )
        };

        if band.first == band.b_size {
            break;
        }
        if band.last < band.b_size.saturating_sub(1) {
            band.b_size = band.last + 1;
        } else {
            // Written values are at least `best - x_dropoff` and at most `best + max`: 16 bits.
            while score_gap_row >= band.best - x_dropoff && band.b_size <= len2 {
                cur[PAD + band.b_size] = score_gap_row as i16;
                if band.b_size.is_multiple_of(RESTRICT_SIZE) {
                    e16[band.b_size] = (score_gap_row - gap_open_extend) as i16;
                }
                score_gap_row -= gap_extend;
                band.b_size += 1;
            }
        }
        if band.b_size <= len2 {
            cur[PAD + band.b_size] = i16::MIN;
            if band.b_size.is_multiple_of(RESTRICT_SIZE) {
                e16[band.b_size] = i16::MIN;
            }
            band.b_size += 1;
        }
    }

    Ok((band.a_offset, band.b_offset, band.best))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1076
/// LOSAT_X_BPGAP16: `restricted_rows16` with the AVX2 rows (`x_avx2_16`), compiled for AVX2 as a
/// whole.
///
/// # Safety
/// AVX2, BMI1 and LZCNT; the arrays of `sc` (16-bit ones included) sized by `restricted_fast`;
/// `fits16`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,bmi1,lzcnt")]
unsafe fn restricted_rows16_avx2(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    reverse: bool,
    tables: &XScoreTables,
    p: &XParams,
    sc: &mut XRestrictedScratch,
) -> Result<(usize, usize, i32), XStop> {
    restricted_rows16::<x_avx2_16::Rows>(
        query,
        query_base,
        subject,
        subject_base,
        len1,
        len2,
        reverse,
        tables,
        p,
        sc,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1076
/// LOSAT_X_BPGAP128: `restricted_rows16` with the rows of `x_w128_16`, eight cells per step,
/// written once over the 128-bit layer of `utils/xdrop_simd.rs` for x86_64 (SSSE3 and SSE4.1) and
/// wasm32 (simd128), compiled for those instructions as a whole.
///
/// # Safety
/// On x86_64 the CPU has SSSE3 and SSE4.1; the arrays of `sc` (16-bit ones included) sized by
/// `restricted_fast`; `fits16`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1"))]
unsafe fn restricted_rows16_w128(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    reverse: bool,
    tables: &XScoreTables,
    p: &XParams,
    sc: &mut XRestrictedScratch,
) -> Result<(usize, usize, i32), XStop> {
    restricted_rows16::<x_w128_16::Rows>(
        query,
        query_base,
        subject,
        subject_base,
        len1,
        len2,
        reverse,
        tables,
        p,
        sc,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1112-1247
// The two row loops sixteen 16-bit cells at a time (AVX2); the module comment ("Sixteen-bit rows")
// gives the argument.
#[cfg(target_arch = "x86_64")]
mod x_avx2_16 {
    use super::{XBand, XParams, XStop, MININT, PAD, RESTRICT_SIZE};
    use core::arch::x86_64::*;

    // `GAPCOL16` and the step of `r10` in `row` are written for columns of ten.
    const _: () = assert!(RESTRICT_SIZE == 10);
    // No NCBI counterpart: the stand-in for `MININT` and lane masks; they hold no value NCBI
    // computes.
    const M16: i16 = i16::MIN;
    /// `VALID16[16 - n..]`: the first `n` lanes set.
    static VALID16: [i16; 32] = [
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
    ];
    /// `GAPCOL16[r..]`: the lanes `i` with `(r + i) % RESTRICT_SIZE == 0` (`r < 10`).
    static GAPCOL16: [i16; 26] = [
        -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0,
    ];

    /// Inclusive prefix maximum over the sixteen lanes (lanes left of lane 0 count as `M16`).
    #[inline(always)]
    unsafe fn prefix_max(v: __m256i) -> __m256i {
        let f1 = _mm256_setr_epi16(M16, 0, 0, 0, 0, 0, 0, 0, M16, 0, 0, 0, 0, 0, 0, 0);
        let f2 = _mm256_setr_epi16(M16, M16, 0, 0, 0, 0, 0, 0, M16, M16, 0, 0, 0, 0, 0, 0);
        let f4 = _mm256_setr_epi16(
            M16, M16, M16, M16, 0, 0, 0, 0, M16, M16, M16, M16, 0, 0, 0, 0,
        );
        let mut t = _mm256_max_epi16(v, _mm256_or_si256(_mm256_slli_si256::<2>(v), f1));
        t = _mm256_max_epi16(t, _mm256_or_si256(_mm256_slli_si256::<4>(t), f2));
        t = _mm256_max_epi16(t, _mm256_or_si256(_mm256_slli_si256::<8>(t), f4));
        // Lane 7 (the low half's maximum) into every lane of the high half; M16 in the low half.
        let b7 = _mm256_setr_epi8(
            14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15,
            14, 15, 14, 15, 14, 15, 14, 15, 14, 15,
        );
        let low_min = _mm256_setr_epi16(
            M16, M16, M16, M16, M16, M16, M16, M16, 0, 0, 0, 0, 0, 0, 0, 0,
        );
        let carry = _mm256_or_si256(
            _mm256_shuffle_epi8(_mm256_permute2x128_si256::<0x08>(t, t), b7),
            low_min,
        );
        _mm256_max_epi16(t, carry)
    }

    /// Lane `i` = `max(carry, p[i - 1])`, lane 0 = `carry` (exclusive prefix from an inclusive one).
    #[inline(always)]
    unsafe fn exclusive(p: __m256i, carry: __m256i) -> __m256i {
        // [0 | low half of p], then each half shifted up one lane, the lane below it filling lane 0.
        let x = _mm256_permute2x128_si256::<0x08>(p, p);
        let up = _mm256_alignr_epi8::<14>(p, x);
        let lane0 = _mm256_setr_epi16(M16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        _mm256_max_epi16(_mm256_or_si256(up, lane0), carry)
    }

    /// One row of `s_RestrictedGappedAlign` (as `x_avx2::row`), sixteen 16-bit cells per step, with
    /// saturating sums. The column gaps are a vector of `e16` (`i16::MIN` outside the columns that
    /// are multiples of `RESTRICT_SIZE`), updated for the kept cells of those columns by a lane mask.
    /// Reads the cells `first - 1 .. b_size + 15` of `prev`, `col[first .. b_size + 16]` and
    /// `e16[first .. b_size + 16]`; writes the cells `first .. b_size + 16` of `cur` (cells from
    /// `b_size` on get `M16`) and `e16` (lanes outside the kept gap columns keep their value).
    /// Returns `score_gap_row` after the row in 32 bits (`MININT` for `M16`).
    ///
    /// # Safety
    /// Inlined into code compiled for AVX2, BMI1 and LZCNT (`restricted_rows16_avx2`); `prev`,
    /// `cur` at least `PAD + b_size + 16` long, `col` and `e16` at least `b_size + 16`; `fits16`
    /// and `best <= limit16`.
    #[inline(always)]
    pub(super) unsafe fn row<const GAP: bool>(
        band: &mut XBand,
        p: &XParams,
        a_index: usize,
        row: &[i32; 32],
        lut: &[u8; 64],
        prev: &[i16],
        cur: &mut [i16],
        col: &[u8],
        e16: &mut [i16],
    ) -> i32 {
        let first0 = band.first;
        let end = band.b_size;
        debug_assert!(first0 < end);
        debug_assert!(prev.len() >= PAD + end + 16 && cur.len() >= PAD + end + 16);
        debug_assert!(col.len() >= end + 16 && e16.len() >= end + 16);
        let prev_p = prev.as_ptr().add(PAD);
        let cur_p = cur.as_mut_ptr().add(PAD);
        let col_p = col.as_ptr();
        let e_p = e16.as_mut_ptr();
        let minv = _mm256_set1_epi16(M16);
        let xv = _mm256_set1_epi16(p.x_dropoff as i16);
        let goev = _mm256_set1_epi16(p.gap_open_extend as i16);
        let gev = _mm256_set1_epi16(p.gap_extend as i16);
        let t0l = _mm_loadu_si128(lut.as_ptr() as *const __m128i);
        let t0h = _mm_loadu_si128(lut.as_ptr().add(16) as *const __m128i);
        let t1l = _mm_loadu_si128(lut.as_ptr().add(32) as *const __m128i);
        let t1h = _mm_loadu_si128(lut.as_ptr().add(48) as *const __m128i);
        let k15 = _mm_set1_epi8(15);
        let ramp0 = _mm256_mullo_epi16(
            _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
            gev,
        );
        let ramp1 = _mm256_add_epi16(ramp0, gev);
        let lane0 = _mm256_setr_epi16(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        let mut best = band.best;
        let mut rv = _mm256_set1_epi16(best as i16);
        let mut thr = _mm256_subs_epi16(rv, xv);
        let mut leading = true;
        let mut first = first0;
        let mut last = first0;
        let mut f_in = M16;
        let mut r10 = first0 % RESTRICT_SIZE;

        let mut b0 = first0;
        while b0 < end {
            let n = (end - b0).min(16);
            let idx = _mm_loadu_si128(col_p.add(b0) as *const __m128i);
            let sel = _mm_cmpgt_epi8(idx, k15);
            let lo = _mm_blendv_epi8(_mm_shuffle_epi8(t0l, idx), _mm_shuffle_epi8(t1l, idx), sel);
            let hi = _mm_blendv_epi8(_mm_shuffle_epi8(t0h, idx), _mm_shuffle_epi8(t1h, idx), sel);
            let m = _mm256_set_m128i(_mm_unpackhi_epi8(lo, hi), _mm_unpacklo_epi8(lo, hi));
            let mut h0 = _mm256_adds_epi16(
                _mm256_loadu_si256(prev_p.add(b0).sub(1) as *const __m256i),
                m,
            );
            if b0 == first0 {
                h0 = _mm256_blendv_epi8(h0, minv, lane0);
            }
            let e = _mm256_loadu_si256(e_p.add(b0) as *const __m256i);
            h0 = _mm256_max_epi16(h0, e);
            let valid = _mm256_loadu_si256(VALID16.as_ptr().add(16 - n) as *const __m256i);
            let mut s = h0;
            if GAP {
                let w = _mm256_blendv_epi8(
                    minv,
                    _mm256_adds_epi16(_mm256_subs_epi16(h0, goev), ramp1),
                    valid,
                );
                let pw = prefix_max(w);
                let f = _mm256_subs_epi16(exclusive(pw, _mm256_set1_epi16(f_in)), ramp0);
                s = _mm256_max_epi16(s, f);
                // Lanes from n on are M16 in w, so lane 15 of the prefix is lane n - 1.
                let pl = _mm256_extract_epi16::<15>(pw) as i16;
                f_in = (pl.max(f_in) as i32 - n as i32 * p.gap_extend).max(M16 as i32) as i16;
            }
            s = _mm256_blendv_epi8(minv, s, valid);
            let above = _mm256_movemask_epi8(_mm256_cmpgt_epi16(s, rv));
            let dropped = if above == 0 {
                _mm256_cmpgt_epi16(thr, s)
            } else {
                let ps = prefix_max(s);
                let ex = exclusive(ps, rv);
                let block_max = _mm256_extract_epi16::<15>(ps) as i16;
                let eq = _mm256_cmpeq_epi16(s, _mm256_set1_epi16(block_max));
                let lane = (_mm256_movemask_epi8(eq) as u32).trailing_zeros() as usize / 2;
                best = block_max as i32;
                band.a_offset = a_index;
                band.b_offset = b0 + lane;
                rv = _mm256_set1_epi16(block_max);
                thr = _mm256_subs_epi16(rv, xv);
                _mm256_cmpgt_epi16(_mm256_subs_epi16(ex, xv), s)
            };
            _mm256_storeu_si256(
                cur_p.add(b0) as *mut __m256i,
                _mm256_blendv_epi8(s, minv, dropped),
            );
            let gap_cols = _mm256_loadu_si256(GAPCOL16.as_ptr().add(r10) as *const __m256i);
            let e_new = _mm256_max_epi16(_mm256_subs_epi16(s, goev), _mm256_subs_epi16(e, gev));
            _mm256_storeu_si256(
                e_p.add(b0) as *mut __m256i,
                _mm256_blendv_epi8(e, e_new, _mm256_andnot_si256(dropped, gap_cols)),
            );
            // Two mask bits per lane; lanes from n on are M16, below the threshold, so dropped.
            let kept = !(_mm256_movemask_epi8(dropped) as u32);
            let any = kept != 0;
            let last_c = b0 + (31 - (kept | 1).leading_zeros() as usize) / 2;
            last = if any { last_c } else { last };
            let first_c = if any {
                b0 + kept.trailing_zeros() as usize / 2
            } else {
                b0 + n
            };
            first = if leading { first_c } else { first };
            leading &= !any;
            r10 = if r10 >= 4 { r10 - 4 } else { r10 + 6 };
            b0 += 16;
        }
        band.best = best;
        band.first = first;
        band.last = last;
        if GAP {
            if f_in == M16 {
                MININT
            } else {
                f_in as i32
            }
        } else {
            // score after the loop: next_score of cell b_size - 1.
            let v = prev[PAD + end - 1];
            (if v == M16 { MININT } else { v as i32 }) + row[col[end] as usize]
        }
    }

    // No NCBI counterpart: the AVX2 kernel of `restricted_rows16`; it holds no value NCBI computes.
    /// The AVX2 sixteen-bit rows (sixteen cells per step) for `super::restricted_rows16` (`super::restricted_rows16_avx2`).
    pub(super) struct Rows;

    impl super::XRows16 for Rows {
        type Consts = ();

        #[inline(always)]
        unsafe fn consts(_: &XParams) -> Self::Consts {}

        #[inline(always)]
        fn fill_col(
            col: &mut [u8],
            from: usize,
            to: usize,
            subject: &[u8],
            subject_base: usize,
            len2: usize,
            reverse: bool,
        ) -> Result<(), XStop> {
            super::fill_col::<{ super::COL16 }>(col, from, to, subject, subject_base, len2, reverse)
        }

        #[inline(always)]
        unsafe fn row<const GAP: bool>(
            band: &mut XBand,
            p: &XParams,
            a_index: usize,
            row: &[i32; 32],
            lut: &[u8; 64],
            prev: &[i16],
            cur: &mut [i16],
            col: &[u8],
            e16: &mut [i16],
            _: &Self::Consts,
        ) -> i32 {
            self::row::<GAP>(band, p, a_index, row, lut, prev, cur, col, e16)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1112-1247
// LOSAT_X_BPGAP128: the two row loops of `x_avx2_16` with eight 16-bit cells per step, over the
// 128-bit layer of `utils/xdrop_simd.rs` (`arch`: x86_64 SSSE3 and SSE4.1, wasm32 simd128). A step
// does for its eight cells what a block of `x_avx2_16::row` does for its sixteen: the same lane
// operations, saturating as there, so the module comment ("Sixteen-bit rows") holds for these
// rows as stated. What is written differently gives the same values:
// 1. Letters are looked up sixteen at a time for two steps (`look16`; `lut32k`: wasm32's swizzle
//    gives 0 for an index of 16 or more, which selects the table half, so each table is the OR of
//    two swizzles).
// 2. Steps whose eight lanes all hold cells (`FULL`) apply no lane mask, and only the first step
//    of a row (`FIRST`) clears lane 0's diagonal; both are const parameters, so no test of them
//    runs per step.
// 3. Lanes moved up across the vector take the lanes moved in from a second vector (`up16_*`,
//    `palignr`): `M16` for the prefix maxima, the carry for the exclusive ones.
// 4. The row gap carried into the next step is a vector of one value (`fv`):
//    `max(carry, lane 7 of the prefix) - n * gap_extend`, saturated at `i16::MIN`, the AVX2 rows'
//    scalar `f_in`.
// 5. A step that raises `best_score` takes `b_offset` from the last lane whose score is above the
//    best before it (above `best_score` and every lane left of it), which is the first lane that
//    holds the step's maximum, the lane the AVX2 rows find with an equality test: the lanes left
//    of it are below that maximum and the lanes right of it at most that maximum.
// 6. `first_b_index` and `last_b_index` are made after the row from the first and the last step
//    with a kept cell and their drop masks, the lanes the AVX2 rows read from the same masks.
// 7. The vector constants are made once per alignment (`Consts`), from values the compiler cannot
//    see (`black_box`) where they are constants: wasm32 engines make a constant again at each use
//    inside a loop.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
mod x_w128_16 {
    use super::{XBand, XParams, XStop, MININT, PAD, RESTRICT_SIZE};
    use crate::utils::xdrop_simd::arch::*;

    // `GAPCOL8` and the step of `r10` in `step` are written for columns of ten.
    const _: () = assert!(RESTRICT_SIZE == 10);
    // No NCBI counterpart: the stand-in for `MININT` and lane masks; they hold no value NCBI
    // computes.
    const M16: i16 = i16::MIN;
    /// `VALID8[8 - n..]`: the first `n` lanes set.
    static VALID8: [i16; 16] = [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0];
    /// `GAPCOL8[r..]`: the lanes `i` with `(r + i) % RESTRICT_SIZE == 0` (`r < 10`).
    static GAPCOL8: [i16; 18] = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0];
    /// Lane 0 set.
    static LANE0: [i16; 8] = [-1, 0, 0, 0, 0, 0, 0, 0];
    /// Byte shuffle: lane 7 to every lane (x86_64's `bcast_last16`; wasm32's does not read it).
    static KB7: [u8; 16] = [
        14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15,
    ];

    /// No NCBI counterpart: the vectors the rows of one alignment use, made once (item 7).
    pub(super) struct Consts {
        /// `M16` in every lane.
        minv: V,
        /// `x_dropoff`, `gap_open_extend`, `gap_extend` and `8 * gap_extend` in every lane.
        xv: V,
        goev: V,
        gev: V,
        ge8: V,
        /// `j * gap_extend` and `(j + 1) * gap_extend` in lane `j`.
        ramp0: V,
        ramp1: V,
        lane0: V,
        kb7: V,
        /// -16 in every byte (`lut32k`).
        minus16: V,
    }

    impl Consts {
        /// # Safety
        /// Inlined into code compiled for the build's 128-bit instructions (x86_64: SSSE3 and
        /// SSE4.1); `fits16` (every value fits 16 bits).
        #[inline(always)]
        pub(super) unsafe fn new(p: &XParams) -> Self {
            let ge = p.gap_extend;
            let mut ramp0 = [0i16; 8];
            let mut ramp1 = [0i16; 8];
            for j in 0..8 {
                ramp0[j] = (j as i32 * ge) as i16;
                ramp1[j] = ((j as i32 + 1) * ge) as i16;
            }
            Consts {
                minv: splat16(core::hint::black_box(M16)),
                xv: splat16(p.x_dropoff as i16),
                goev: splat16(p.gap_open_extend as i16),
                gev: splat16(ge as i16),
                ge8: splat16((8 * ge) as i16),
                ramp0: load(ramp0.as_ptr()),
                ramp1: load(ramp1.as_ptr()),
                lane0: load(core::hint::black_box(&LANE0).as_ptr()),
                kb7: load(core::hint::black_box(&KB7).as_ptr()),
                minus16: splat8(core::hint::black_box(-16)),
            }
        }
    }

    /// No NCBI counterpart: what a row carries from one step to the next.
    struct Steps {
        /// `best_score` in every lane, and `best_score - x_dropoff` (saturated).
        rv: V,
        thr: V,
        /// Rows where a gap opens in A: the row gap entering the step, in every lane (item 4).
        fv: V,
        /// `b0 % RESTRICT_SIZE` of the step.
        r10: usize,
        /// No cell of the row kept yet; else the first step with a kept cell and its kept lanes
        /// (two mask bits per lane).
        leading: bool,
        first_b0: usize,
        first_kept: u32,
        /// The last step with a kept cell and its kept lanes (`first_b_index` and lane 0 before
        /// any).
        last_b0: usize,
        last_kept: u32,
        /// A step raised `best_score`; `b_offset` of the last one.
        raised: bool,
        b_offset: usize,
    }

    /// Inclusive prefix maximum over the eight lanes (lanes left of lane 0 count as `M16`, the
    /// lanes of `minv`).
    #[inline(always)]
    unsafe fn prefix_max(v: V, minv: V) -> V {
        let t = max16(v, up16_1(v, minv));
        let t = max16(t, up16_2(t, minv));
        max16(t, up16_4(t, minv))
    }

    /// The matrix entries of the letters `col[b0 .. b0 + 16]` (all below 28): those of the cells
    /// `b0 .. b0 + 8`, then those of `b0 + 8 .. b0 + 16` (item 1). `t` holds the row's table as
    /// `XScoreTables::lut` does: low bytes of entries 0-15, their high bytes, low bytes of
    /// entries 16-31, their high bytes.
    #[inline(always)]
    unsafe fn look16(t: &[V; 4], col_p: *const u8, b0: usize, minus16: V) -> (V, V) {
        let idx = load(col_p.add(b0));
        let lo = lut32k(t[0], t[2], idx, minus16);
        let hi = lut32k(t[1], t[3], idx, minus16);
        (unpacklo8(lo, hi), unpackhi8(lo, hi))
    }

    /// One step: the cells `b0 .. b0 + n` of a row (`n` = 8 when `FULL`, else 1 to 7: the last
    /// step of the row), `m` the entries of their letters; the operations of a block of
    /// `x_avx2_16::row` on eight lanes (items 2-6).
    ///
    /// # Safety
    /// As `row`; `b0 < b_size`.
    #[inline(always)]
    unsafe fn step<const GAP: bool, const FIRST: bool, const FULL: bool>(
        c: &Consts,
        st: &mut Steps,
        m: V,
        b0: usize,
        n: usize,
        gap_extend: i32,
        prev_p: *const i16,
        cur_p: *mut i16,
        e_p: *mut i16,
    ) {
        // Diagonal: prev[b - 1] + matrix_row[col[b]].
        let mut h0 = adds16(load(prev_p.add(b0).sub(1)), m);
        if FIRST {
            // Initial score is `MININT` at the start of the row.
            h0 = blend8(h0, c.minv, c.lane0);
        }
        // Column gap (`i16::MIN` outside the columns that are multiples of 10).
        let e = load(e_p.add(b0));
        h0 = max16(h0, e);
        // Lanes past `b_size` (the last step of the row only) hold no cell.
        let valid = if FULL {
            None
        } else {
            Some(load(VALID8.as_ptr().add(8 - n)))
        };
        let mut s = h0;
        if GAP {
            // Row gap: F(i) = max(f_in, max_{j<i} (h0[j] - goe + (j+1) ge)) - i ge.
            let mut w = adds16(subs16(h0, c.goev), c.ramp1);
            if let Some(valid) = valid {
                w = blend8(c.minv, w, valid);
            }
            let pw = prefix_max(w, c.minv);
            let f = subs16(max16(up16_1(pw, st.fv), st.fv), c.ramp0);
            s = max16(s, f);
            // Row gap after the step's last cell (cell b0 + n): lanes from n on are M16 in w, so
            // lane 7 of the prefix is lane n - 1.
            let gn = if FULL {
                c.ge8
            } else {
                splat16((n as i32 * gap_extend) as i16)
            };
            st.fv = subs16(max16(bcast_last16(pw, c.kb7), st.fv), gn);
        }
        if let Some(valid) = valid {
            s = blend8(c.minv, s, valid);
        }
        // X-drop: best_score when a cell is tested = max(best at the step, cells left of it).
        // When no cell of the step exceeds `best`, that is `best` for every cell.
        let dropped = if mask8(gt16(s, st.rv)) == 0 {
            gt16(st.thr, s)
        } else {
            let ps = prefix_max(s, c.minv);
            let ex = max16(up16_1(ps, st.rv), st.rv);
            // Item 5: two mask bits per lane; the last lane above `ex` holds the maximum first.
            let above = mask8(gt16(s, ex));
            st.raised = true;
            st.b_offset = b0 + (31 - above.leading_zeros() as usize) / 2;
            st.rv = bcast_last16(ps, c.kb7);
            st.thr = subs16(st.rv, c.xv);
            gt16(subs16(ex, c.xv), s)
        };
        store(cur_p.add(b0), blend8(s, c.minv, dropped));
        let gap_cols = load(GAPCOL8.as_ptr().add(st.r10));
        let e_new = max16(subs16(s, c.goev), subs16(e, c.gev));
        store(e_p.add(b0), blend8(e, e_new, andnot(dropped, gap_cols)));
        // Two mask bits per lane; lanes from n on are M16, below the threshold, so dropped.
        let kept = !mask8(dropped) & 0xFFFF;
        if kept != 0 {
            if st.leading {
                st.leading = false;
                st.first_b0 = b0;
                st.first_kept = kept;
            }
            st.last_b0 = b0;
            st.last_kept = kept;
        }
        st.r10 = if st.r10 >= 2 { st.r10 - 2 } else { st.r10 + 8 };
    }

    /// The steps of the sixteen cells from `b0` whose entries are `lo` and `hi`, of which the first
    /// `n` (1 to 16) are cells of the row (`FIRST`: `b0` is the row's first cell).
    ///
    /// # Safety
    /// As `row`.
    #[inline(always)]
    unsafe fn pair<const GAP: bool, const FIRST: bool>(
        c: &Consts,
        st: &mut Steps,
        (lo, hi): (V, V),
        b0: usize,
        n: usize,
        gap_extend: i32,
        prev_p: *const i16,
        cur_p: *mut i16,
        e_p: *mut i16,
    ) {
        if n >= 8 {
            step::<GAP, FIRST, true>(c, st, lo, b0, 8, gap_extend, prev_p, cur_p, e_p);
            if n == 16 {
                step::<GAP, false, true>(c, st, hi, b0 + 8, 8, gap_extend, prev_p, cur_p, e_p);
            } else if n > 8 {
                step::<GAP, false, false>(c, st, hi, b0 + 8, n - 8, gap_extend, prev_p, cur_p, e_p);
            }
        } else {
            step::<GAP, FIRST, false>(c, st, lo, b0, n, gap_extend, prev_p, cur_p, e_p);
        }
    }

    /// One row of `s_RestrictedGappedAlign` (as `x_avx2_16::row`), eight 16-bit cells per step,
    /// with saturating sums. The column gaps are a vector of `e16` (`i16::MIN` outside the columns
    /// that are multiples of `RESTRICT_SIZE`), updated for the kept cells of those columns by a lane
    /// mask. Reads the cells `first - 1 .. b_size + 7` of `prev`, `col[first .. b_size + 16]` and
    /// `e16[first .. b_size + 8]`; writes the cells `first .. b_size + 8` of `cur` (cells from
    /// `b_size` on get `M16`) and `e16` (lanes outside the kept gap columns keep their value).
    /// Returns `score_gap_row` after the row in 32 bits (`MININT` for `M16`).
    ///
    /// # Safety
    /// Inlined into code compiled for the build's 128-bit instructions (`restricted_rows16_w128`;
    /// x86_64: SSSE3 and SSE4.1); `prev`, `cur` at least `PAD + b_size + 16` long, `col` and `e16`
    /// at least `b_size + 16`; `fits16` and `best <= limit16`.
    #[inline(always)]
    pub(super) unsafe fn row<const GAP: bool>(
        band: &mut XBand,
        p: &XParams,
        a_index: usize,
        row: &[i32; 32],
        lut: &[u8; 64],
        prev: &[i16],
        cur: &mut [i16],
        col: &[u8],
        e16: &mut [i16],
        c: &Consts,
    ) -> i32 {
        let first0 = band.first;
        let end = band.b_size;
        debug_assert!(first0 < end);
        debug_assert!(prev.len() >= PAD + end + 16 && cur.len() >= PAD + end + 16);
        debug_assert!(col.len() >= end + 16 && e16.len() >= end + 16);
        let prev_p = prev.as_ptr().add(PAD);
        let cur_p = cur.as_mut_ptr().add(PAD);
        let col_p = col.as_ptr();
        let e_p = e16.as_mut_ptr();
        let lut_p = lut.as_ptr();
        let t = [
            load(lut_p),
            load(lut_p.add(16)),
            load(lut_p.add(32)),
            load(lut_p.add(48)),
        ];
        let ge = p.gap_extend;
        let rv = splat16(band.best as i16);
        let mut st = Steps {
            rv,
            thr: subs16(rv, c.xv),
            fv: c.minv,
            r10: first0 % RESTRICT_SIZE,
            leading: true,
            first_b0: 0,
            first_kept: 0,
            last_b0: first0,
            last_kept: 1,
            raised: false,
            b_offset: 0,
        };

        // The first sixteen cells (the first step clears lane 0's diagonal), the full pairs, and
        // the cells left (item 2).
        let entries = look16(&t, col_p, first0, c.minus16);
        let n = (end - first0).min(16);
        pair::<GAP, true>(c, &mut st, entries, first0, n, ge, prev_p, cur_p, e_p);
        let mut b0 = first0 + 16;
        while b0 + 16 <= end {
            let (lo, hi) = look16(&t, col_p, b0, c.minus16);
            step::<GAP, false, true>(c, &mut st, lo, b0, 8, ge, prev_p, cur_p, e_p);
            step::<GAP, false, true>(c, &mut st, hi, b0 + 8, 8, ge, prev_p, cur_p, e_p);
            b0 += 16;
        }
        if b0 < end {
            let entries = look16(&t, col_p, b0, c.minus16);
            pair::<GAP, false>(c, &mut st, entries, b0, end - b0, ge, prev_p, cur_p, e_p);
        }

        // Every lane of `rv` holds `best_score`.
        band.best = lane0(st.rv);
        if st.raised {
            band.a_offset = a_index;
            band.b_offset = st.b_offset;
        }
        // Item 6: no kept cell leaves `first_b_index` at `b_size`.
        band.first = if st.leading {
            end
        } else {
            st.first_b0 + st.first_kept.trailing_zeros() as usize / 2
        };
        band.last = st.last_b0 + (31 - st.last_kept.leading_zeros() as usize) / 2;
        if GAP {
            let f = lane0(st.fv);
            if f == M16 as i32 {
                MININT
            } else {
                f
            }
        } else {
            // score after the loop: next_score of cell b_size - 1.
            let v = prev[PAD + end - 1];
            (if v == M16 { MININT } else { v as i32 }) + row[col[end] as usize]
        }
    }

    // No NCBI counterpart: the 128-bit kernel of `restricted_rows16`; it holds no value NCBI computes.
    /// The 128-bit sixteen-bit rows (eight cells per step) for `super::restricted_rows16` (`super::restricted_rows16_w128`).
    pub(super) struct Rows;

    impl super::XRows16 for Rows {
        type Consts = Consts;

        #[inline(always)]
        unsafe fn consts(p: &XParams) -> Self::Consts {
            Consts::new(p)
        }

        #[inline(always)]
        fn fill_col(
            col: &mut [u8],
            from: usize,
            to: usize,
            subject: &[u8],
            subject_base: usize,
            len2: usize,
            reverse: bool,
        ) -> Result<(), XStop> {
            super::fill_col::<{ super::COL16_W128 }>(
                col,
                from,
                to,
                subject,
                subject_base,
                len2,
                reverse,
            )
        }

        #[inline(always)]
        unsafe fn row<const GAP: bool>(
            band: &mut XBand,
            p: &XParams,
            a_index: usize,
            row: &[i32; 32],
            lut: &[u8; 64],
            prev: &[i16],
            cur: &mut [i16],
            col: &[u8],
            e16: &mut [i16],
            c: &Self::Consts,
        ) -> i32 {
            self::row::<GAP>(band, p, a_index, row, lut, prev, cur, col, e16, c)
        }
    }
}

#[cfg(test)]
#[path = "x_restricted_gap_tests.rs"]
mod tests;
