//! EXPERIMENT (LOSAT_X_BPGAP / LOSAT_X_BPGAPSHADOW): the restricted gapped alignment of the blastp
//! preliminary stage (NCBI `s_RestrictedGappedAlign`) with its score rows in separate arrays,
//! eight cells of a row at a time (x86_64 with AVX2; elsewhere the reference function runs).
//!
//! NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1107-1122,1143-1147,1171
//! ```c
//!         if (a_index % RESTRICT_SIZE != 0) {
//! ...
//!             for (b_index = first_b_index; b_index < b_size; b_index++) {
//!
//!                 b_ptr += b_increment;
//!                 next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
//!
//!                 if (b_index != b_gap) {
//! ...
//!                     if (best_score - score > x_dropoff) {
//! ...
//!                     b_gap += RESTRICT_SIZE;
//!                     score_gap_col = score_array[b_index].best_gap;
//!
//!                     if (score < score_gap_col)
//!                         score = score_gap_col;
//! ...
//!             score_gap_row = score;
//! ```
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
//! - The values are 32-bit integers as in the reference: no lanes of 8 or 16 bits, no rebasing,
//!   no saturation. The matrix row is read from a table of 16-bit entries built from the
//!   reference's own score function (`blastp_score_from_row`), so every score is the reference's
//!   (BLOSUM62's `BLAST_SCORE_MIN` = -32768 for the gap letter fits). A matrix with an entry
//!   outside 16 bits, a letter code of 28 or more, a negative gap cost, a gap cost above 65,536 or
//!   an X-drop above 2^28 are not handled here; the reference function runs instead.
//! - The scalar rows (`row_nogap_scalar`, `row_gap_scalar`) are the reference's loops on these
//!   arrays, cell by cell. They were slower than the reference function, so a CPU without AVX2
//!   runs the reference function; the tests run them to check the arrays apart from the vector
//!   code.
//!
//! The variables are read once. `LOSAT_X_BPGAP` (a strict-set switch, `utils/xswitch.rs`: on
//! unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set) runs this module;
//! `LOSAT_X_BPGAPSHADOW` runs both for every restricted alignment and panics if they differ, and
//! `main` prints one `[X_SHADOW] BPGAP` line.

#![allow(clippy::too_many_arguments)]

use super::{
    blastp_query_residue, blastp_score_from_row, blastp_score_row_for_query,
    blastp_subject_residue, restricted_gapped_align_protein, BlastpScoreMatrix, GapAlignScratch,
    GAP_MININT, RESTRICT_SIZE,
};
use crate::config::ScoringMatrix;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1096-1097
// ```c
//         score = MININT;
//         score_gap_row = MININT;
// ```
// The reference's minus infinity (`GAP_MININT`, as in `restricted_gapped_align_protein_impl`).
const MININT: i32 = GAP_MININT;
/// No NCBI counterpart: cells kept before cell 0 of the two rows (the eight-cell block of the first
/// cell reads the cell left of it) and after the last one (a block may run up to seven cells past
/// `b_size`); they hold no value NCBI computes.
const PAD: usize = 8;
const TAIL: usize = 24;
/// No NCBI counterpart: the largest gap cost and X-drop this module takes, so that no sum below
/// overflows 32 bits (otherwise the reference function runs).
const MAX_GAP_COST: i32 = 1 << 16;
const MAX_X_DROP: i32 = 1 << 28;

/// EXPERIMENT (LOSAT_X_BPGAP / LOSAT_X_BPGAPSHADOW): which restricted alignment runs.
// No NCBI counterpart: names the LOSAT code path; it does not change any value NCBI computes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XBpGap {
    /// `restricted_gapped_align_protein` (the port of NCBI's loop).
    Reference,
    /// This module (the reference when it does not take the problem).
    Fast,
    /// Both; panics if they differ.
    Shadow,
}

// No NCBI counterpart: reads LOSAT_X_BPGAP / LOSAT_X_BPGAPSHADOW once (the shadow variable
// first, set to any value; LOSAT_X_BPGAP is a strict-set switch); it does not change any value NCBI
// computes.
pub(super) fn x_bpgap() -> XBpGap {
    static MODE: OnceLock<XBpGap> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPGAPSHADOW").is_some() {
            XBpGap::Shadow
        } else if crate::utils::xswitch::on("LOSAT_X_BPGAP") {
            XBpGap::Fast
        } else {
            XBpGap::Reference
        }
    })
}

// No NCBI counterpart: counters of the shadow summary and of LOSAT_X_STATS; they do not change any
// value NCBI computes.
/// Restricted alignments computed by this module (in shadow mode, those compared).
static FAST_CALLS: AtomicU64 = AtomicU64::new(0);
/// Restricted alignments this module did not take (the reference ran alone).
static FALLBACK_CALLS: AtomicU64 = AtomicU64::new(0);

/// The summary lines `main` prints at exit: the shadow summary, and with LOSAT_X_STATS the
/// number of alignments this module computed and passed on.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    let fast = FAST_CALLS.load(Ordering::Relaxed);
    let fallback = FALLBACK_CALLS.load(Ordering::Relaxed);
    if x_bpgap() == XBpGap::Shadow {
        eprintln!(
            "[X_SHADOW] BPGAP compared={fast} not_taken={fallback} ({})",
            if fast == 0 {
                "nothing compared: no restricted gapped alignment ran"
            } else {
                "every restricted gapped alignment identical"
            }
        );
    }
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
}

impl XRestrictedScratch {
    // No NCBI counterpart: empty until the first alignment that runs here.
    pub(super) fn new() -> Self {
        Self {
            h: [Vec::new(), Vec::new()],
            e10: Vec::new(),
            col: Vec::new(),
        }
    }
}

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1078-1083,1115
// ```c
//         if (!(gap_align->positionBased)) {
//             if(reverse_sequence)
//                 matrix_row = matrix[ A[ M - a_index ] ];
//             else
//                 matrix_row = matrix[ A[ a_index ] ];
//         }
// ...
//                 next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
// ```
/// `matrix[q][s]` for the letters `q, s < 28`, read through the reference's own functions; for
/// the vector kernel each row also as four 16-byte tables (low and high bytes of entries 0-15,
/// then of entries 16-31).
struct XScoreTables {
    rows: [[i32; 32]; 28],
    lut: [[u8; 64]; 28],
}

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1078-1083
// ```c
//             if(reverse_sequence)
//                 matrix_row = matrix[ A[ M - a_index ] ];
//             else
//                 matrix_row = matrix[ A[ a_index ] ];
// ```
// Built once per matrix. None: an entry does not fit 16 bits (the reference runs).
fn x_score_tables(score_matrix: BlastpScoreMatrix<'_>) -> Option<&'static XScoreTables> {
    fn build<const BLOSUM62: bool>(
        score_matrix: BlastpScoreMatrix<'_>,
    ) -> Option<Box<XScoreTables>> {
        let mut t = Box::new(XScoreTables {
            rows: [[0; 32]; 28],
            lut: [[0; 64]; 28],
        });
        for q in 0..28u8 {
            let matrix_row = blastp_score_row_for_query::<BLOSUM62>(score_matrix, q);
            for s in 0..28u8 {
                let v = blastp_score_from_row::<BLOSUM62>(score_matrix, matrix_row, q, s);
                let v16 = i16::try_from(v).ok()?;
                t.rows[q as usize][s as usize] = v;
                let [lo, hi] = v16.to_le_bytes();
                let half = (s as usize / 16) * 32;
                t.lut[q as usize][half + s as usize % 16] = lo;
                t.lut[q as usize][half + 16 + s as usize % 16] = hi;
            }
        }
        Some(t)
    }
    static TABLES: [OnceLock<Option<Box<XScoreTables>>>; 9] = [const { OnceLock::new() }; 9];
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

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:4273-4278,4307-4316
// ```c
//          if (restricted_alignment) {
//             score_left = s_RestrictedGappedAlign(query, subject+subject_shift,
//                                    q_length, s_length,
//                                    &private_q_start, &private_s_start,
//                                    gap_align, score_params,
//                                    init_hsp->offsets.qs_offsets.q_off, TRUE);
// ...
//          if (restricted_alignment) {
//             score_right = s_RestrictedGappedAlign(
//                                    query+init_hsp->offsets.qs_offsets.q_off,
//                                    subject+init_hsp->offsets.qs_offsets.s_off,
//                                    query_length-q_length,
//                                    subject_length-s_length,
// ```
/// Dispatch point of the two `s_RestrictedGappedAlign` calls of `s_BlastProtGappedAlignment`: the
/// reference function, this module, or both (shadow). The arguments are those of
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
    let mode = x_bpgap();
    if mode != XBpGap::Reference {
        // Without AVX2 the reference function runs: the scalar rows on the same arrays were
        // slower than it (they serve the tests as the cell-by-cell version of the arrays).
        let fast = if cpu_has_avx2() {
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
                true,
            )
        } else {
            None
        };
        if let Some(fast) = fast {
            FAST_CALLS.fetch_add(1, Ordering::Relaxed);
            if mode == XBpGap::Fast {
                return fast;
            }
            let reference = restricted_gapped_align_protein(
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
            );
            assert_eq!(
                fast, reference,
                "LOSAT_X_BPGAPSHADOW: the restricted gapped alignment differs from the reference \
                 (len1={len1} len2={len2} gap_open={gap_open} gap_extend={gap_extend} \
                 x_drop={x_drop} reverse={reverse} query_base={query_base} subject_base={subject_base})"
            );
            return reference;
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

// No NCBI counterpart: the CPU check of the vector kernel; it does not change any value NCBI
// computes.
fn cpu_has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        static AVX2: OnceLock<bool> = OnceLock::new();
        *AVX2.get_or_init(|| {
            is_x86_feature_detected!("avx2")
                && is_x86_feature_detected!("bmi1")
                && is_x86_feature_detected!("lzcnt")
        })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1096-1098,1250-1251
// ```c
//         score = MININT;
//         score_gap_row = MININT;
//         last_b_index = first_b_index;
// ...
//         if (first_b_index == b_size)
//             break;
// ```
/// The band and the best cell, the variables of `s_RestrictedGappedAlign` that live across rows.
struct XBand {
    first: usize,
    b_size: usize,
    last: usize,
    best: i32,
    a_offset: usize,
    b_offset: usize,
}

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1028-1038,1053-1067,1076-1098
// ```c
//     gap_open = score_params->gap_open;
//     gap_extend = score_params->gap_extend;
//     gap_open_extend = gap_open + gap_extend;
//     x_dropoff = gap_align->gap_x_dropoff;
//
//     if (x_dropoff < gap_open_extend)
//         x_dropoff = gap_open_extend;
//
//     if(N <= 0 || M <= 0)
//         return 0;
// ...
//     score = -gap_open_extend;
//     score_array[0].best = 0;
//     score_array[0].best_gap = -gap_open_extend;
//
//     for (i = 1; i <= N; i++) {
//         if (score < -x_dropoff)
//             break;
//
//         score_array[i].best = score;
//         score_array[i].best_gap = score - gap_open_extend;
//         score -= gap_extend;
//     }
//
//     b_size = i;
//     best_score = 0;
//     first_b_index = 0;
//     b_gap = 0;
// ...
//     for (a_index = 1; a_index <= M; a_index++) {
// ```
/// `s_RestrictedGappedAlign` on the arrays of `sc`. `use_avx2` selects the vector rows (only
/// when the CPU has AVX2). None: the problem is not taken (module comment); nothing outside
/// `sc` has been changed.
fn restricted_fast(
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
    use_avx2: bool,
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

    let cells = PAD + len2 + 2 + TAIL;
    for h in sc.h.iter_mut() {
        if h.len() < cells {
            h.resize(cells, MININT);
        }
    }
    if sc.e10.len() < len2 / RESTRICT_SIZE + 2 {
        sc.e10.resize(len2 / RESTRICT_SIZE + 2, MININT);
    }
    if sc.col.len() < len2 + 2 + TAIL {
        sc.col.resize(len2 + 2 + TAIL, 0);
    }
    let p = XParams {
        x_dropoff,
        gap_open_extend,
        gap_extend,
    };
    #[cfg(target_arch = "x86_64")]
    if use_avx2 {
        // SAFETY: `use_avx2` is true only when the CPU has AVX2, BMI1 and LZCNT (`cpu_has_avx2`,
        // which the tests also call). The arrays of `sc` were sized for `len2` above.
        return unsafe {
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
        };
    }
    let _ = use_avx2;
    restricted_rows::<false>(
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

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1076
// ```c
//     for (a_index = 1; a_index <= M; a_index++) {
// ```
/// The rows with the vector kernel, compiled for AVX2 as a whole (no call per row).
///
/// # Safety
/// AVX2, BMI1 and LZCNT; `sc` sized for `len2` by `restricted_fast`.
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
) -> Option<(usize, usize, i32)> {
    restricted_rows::<true>(
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

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1053-1067,1076-1098
// ```c
//     score = -gap_open_extend;
//     score_array[0].best = 0;
//     score_array[0].best_gap = -gap_open_extend;
//
//     for (i = 1; i <= N; i++) {
//         if (score < -x_dropoff)
//             break;
//
//         score_array[i].best = score;
//         score_array[i].best_gap = score - gap_open_extend;
//         score -= gap_extend;
//     }
//
//     b_size = i;
//     best_score = 0;
//     first_b_index = 0;
//     b_gap = 0;
// ...
//     for (a_index = 1; a_index <= M; a_index++) {
// ```
/// Row 0 and the row loop of `s_RestrictedGappedAlign` (`AVX2`: the vector rows; the caller
/// guarantees the CPU features and has sized the arrays of `sc` for `len2`).
#[inline(always)]
fn restricted_rows<const AVX2: bool>(
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
) -> Option<(usize, usize, i32)> {
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
            return None;
        }
        // Letters for the cells up to `b_size` and the block that may run past it. A letter is
        // below 28 (checked when it is copied), so `& 31` in the scalar rows changes no index.
        let need = (band.b_size + 9).min(col_cap);
        if col_ready < need {
            let upto = (need + 56).min(col_cap);
            if !fill_col(col, col_ready, upto, subject, subject_base, len2, reverse) {
                return None;
            }
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

        // NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1250-1251,1271-1296
        // ```c
        //         if (first_b_index == b_size)
        //             break;
        // ...
        //         if (last_b_index < b_size - 1) {
        //             b_size = last_b_index + 1;
        //         }
        //         else {
        // ...
        //             while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
        //                 score_array[b_size].best = score_gap_row;
        //                 score_array[b_size].best_gap = score_gap_row - gap_open_extend;
        //                 score_gap_row -= gap_extend;
        //                 b_size++;
        //             }
        //         }
        //
        //         if (b_size <= N) {
        //             score_array[b_size].best = MININT;
        //             score_array[b_size].best_gap = MININT;
        //             b_size++;
        //         }
        // ```
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
                if band.b_size % RESTRICT_SIZE == 0 {
                    e10[band.b_size / RESTRICT_SIZE] = score_gap_row - gap_open_extend;
                }
                score_gap_row -= gap_extend;
                band.b_size += 1;
            }
        }
        if band.b_size <= len2 {
            cur[PAD + band.b_size] = MININT;
            if band.b_size % RESTRICT_SIZE == 0 {
                e10[band.b_size / RESTRICT_SIZE] = MININT;
            }
            band.b_size += 1;
        }
    }

    Some((band.a_offset, band.b_offset, band.best))
}

/// `gap_open_extend`, `gap_extend` and `x_dropoff` of `s_RestrictedGappedAlign`.
struct XParams {
    x_dropoff: i32,
    gap_open_extend: i32,
    gap_extend: i32,
}

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1091-1094,1114-1115
// ```c
//         if(reverse_sequence)
//             b_ptr = &B[N - first_b_index];
//         else
//             b_ptr = &B[first_b_index];
// ...
//                 b_ptr += b_increment;
//                 next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
// ```
/// `col[c]` for `c` in `from..to`: the letter the reference reads for `next_score` at
/// `b_index = c - 1` (`blastp_subject_residue`, sentinels included); 0 for cell 0 and past cell
/// `len2 + 1`, which no cell of the reference reads. The cells whose letter lies inside `subject`
/// are one slice copy (forward: `subject[subject_base + c]`; reverse: `subject[subject_base + len2
/// - c]`, the indices that function computes); the others go through that function. False: a
/// letter code of 28 or more (not taken).
fn fill_col(
    col: &mut [u8],
    from: usize,
    to: usize,
    subject: &[u8],
    subject_base: usize,
    len2: usize,
    reverse: bool,
) -> bool {
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
    col[from..to].iter().all(|&v| v < 28)
}

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1107,1173
// ```c
//         if (a_index % RESTRICT_SIZE != 0) {
// ...
//         else {
// ```
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
        // cell `len2 + 1`, `col` and `e10` cover every index the row reads (sizes set in
        // `restricted_fast`; `b_size <= len2 + 1`).
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

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1112-1171
// ```c
//             for (b_index = first_b_index; b_index < b_size; b_index++) {
//
//                 b_ptr += b_increment;
//                 next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
//
//                 if (b_index != b_gap) {
//                     if (best_score - score > x_dropoff) {
//                         score_array[b_index].best = MININT;
//                         if (b_index == first_b_index)
//                             first_b_index++;
//                     }
//                     else {
//                         last_b_index = b_index;
//                         if (score > best_score) {
//                             best_score = score;
//                             *a_offset = a_index;
//                             *b_offset = b_index;
//                         }
//                         score_array[b_index].best = score;
//                     }
//                 }
//                 else {
//                     b_gap += RESTRICT_SIZE;
//                     score_gap_col = score_array[b_index].best_gap;
//
//                     if (score < score_gap_col)
//                         score = score_gap_col;
//
//                     if (best_score - score > x_dropoff) {
// ...
//                         score_gap_col -= gap_extend;
//                         score_array[b_index].best_gap =
//                                    MAX(score - gap_open_extend, score_gap_col);
//                         score_array[b_index].best = score;
//                     }
//                 }
//                 score = next_score;
//             }
//
//             score_gap_row = score;
// ```
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

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1177-1247
// ```c
//             for (b_index = first_b_index; b_index < b_size; b_index++) {
//
//                 b_ptr += b_increment;
//                 next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
//
//                 if (b_index != b_gap) {
//                     if (score < score_gap_row)
//                         score = score_gap_row;
//
//                     if (best_score - score > x_dropoff) {
// ...
//                         score_gap_row -= gap_extend;
//                         score_gap_row = MAX(score - gap_open_extend,
//                                             score_gap_row);
//                         score_array[b_index].best = score;
//                     }
//                 }
//                 else {
//                     b_gap += RESTRICT_SIZE;
//                     score_gap_col = score_array[b_index].best_gap;
//
//                     if (score < score_gap_col)
//                         score = score_gap_col;
//                     if (score < score_gap_row)
//                         score = score_gap_row;
// ...
//                         score_gap_row -= gap_extend;
//                         score_gap_col -= gap_extend;
//                         score_array[b_index].best_gap =
//                                  MAX(score - gap_open_extend, score_gap_col);
//                         score_gap_row = MAX(score - gap_open_extend,
//                                             score_gap_row);
//                         score_array[b_index].best = score;
// ```
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

// NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:1112-1247
// ```c
//             for (b_index = first_b_index; b_index < b_size; b_index++) {
//
//                 b_ptr += b_increment;
//                 next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
// ...
//                     if (best_score - score > x_dropoff) {
//                         score_array[b_index].best = MININT;
//                         if (b_index == first_b_index)
//                             first_b_index++;
//                     }
//                     else {
//                         last_b_index = b_index;
//                         if (score > best_score) {
//                             best_score = score;
// ...
//                         score_gap_row -= gap_extend;
//                         score_gap_row = MAX(score - gap_open_extend,
//                                             score_gap_row);
// ```
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
                // score = MININT at the start of the row.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::gapalign::GapAlignScratch;

    /// xorshift64*: a fixed generator, so every failure can be replayed.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn reference(
        q: &[u8],
        qb: usize,
        s: &[u8],
        sb: usize,
        len1: usize,
        len2: usize,
        m: BlastpScoreMatrix<'_>,
        go: i32,
        ge: i32,
        x: i32,
        rev: bool,
        scratch: &mut GapAlignScratch,
    ) -> (usize, usize, i32) {
        restricted_gapped_align_protein(
            q,
            qb,
            s,
            sb,
            len1,
            len2,
            m,
            go,
            ge,
            x,
            rev,
            &mut scratch.dp_mem,
            &mut scratch.dp_mem_alloc,
        )
    }

    /// Every kernel this machine runs (scalar, and AVX2 when present) against the reference.
    fn check(
        q: &[u8],
        qb: usize,
        s: &[u8],
        sb: usize,
        len1: usize,
        len2: usize,
        m: BlastpScoreMatrix<'_>,
        go: i32,
        ge: i32,
        x: i32,
        rev: bool,
        scratch: &mut GapAlignScratch,
    ) -> (usize, usize, i32) {
        let want = reference(q, qb, s, sb, len1, len2, m, go, ge, x, rev, scratch);
        let kernels: &[bool] = if cpu_has_avx2() {
            &[false, true]
        } else {
            &[false]
        };
        for &avx2 in kernels {
            let got = restricted_fast(
                q,
                qb,
                s,
                sb,
                len1,
                len2,
                m,
                go,
                ge,
                x,
                rev,
                &mut scratch.x_restricted,
                avx2,
            )
            .expect("problem taken");
            assert_eq!(
                got, want,
                "avx2={avx2} len1={len1} len2={len2} qb={qb} sb={sb} go={go} ge={ge} x={x} rev={rev}"
            );
        }
        want
    }

    /// A random protein over `alphabet` letters (codes 1..=alphabet), with a few gap letters (0)
    /// and stop letters (25) when `odd` is set.
    fn random_seq(rng: &mut Rng, len: usize, alphabet: u64, odd: bool) -> Vec<u8> {
        (0..len)
            .map(|_| {
                if odd && rng.below(50) == 0 {
                    [0u8, 25, 21, 23, 24, 26, 27][rng.below(7) as usize]
                } else {
                    1 + rng.below(alphabet) as u8
                }
            })
            .collect()
    }

    /// `src` with substitutions and insertions / deletions, so that the alignment runs long and
    /// uses gaps.
    fn mutate(rng: &mut Rng, src: &[u8], rate: u64) -> Vec<u8> {
        let mut out = Vec::with_capacity(src.len() + 16);
        for &c in src {
            match rng.below(100) {
                r if r < rate => out.push(1 + rng.below(24) as u8),
                r if r < rate + 2 => {}
                r if r < rate + 4 => {
                    out.push(c);
                    for _ in 0..1 + rng.below(4) {
                        out.push(1 + rng.below(24) as u8);
                    }
                }
                _ => out.push(c),
            }
        }
        out
    }

    /// The two calls of `s_BlastProtGappedAlignment` (left, reverse; right, forward) at a seed
    /// inside `q` and `s`, with the arguments `blastp_score_only_gapped_alignment_with_scratch`
    /// gives them (no subject shift).
    fn check_seed(
        q: &[u8],
        s: &[u8],
        q_start: usize,
        s_start: usize,
        m: BlastpScoreMatrix<'_>,
        go: i32,
        ge: i32,
        x: i32,
        scratch: &mut GapAlignScratch,
    ) {
        let q_length = q_start + 1;
        let s_length = s_start + 1;
        check(q, 0, s, 0, q_length, s_length, m, go, ge, x, true, scratch);
        if q_length < q.len() && s_length < s.len() {
            check(
                q,
                q_start,
                s,
                s_start,
                q.len() - q_length,
                s.len() - s_length,
                m,
                go,
                ge,
                x,
                false,
                scratch,
            );
        }
    }

    const MATRICES: [BlastpScoreMatrix<'static>; 4] = [
        BlastpScoreMatrix::Blosum62,
        BlastpScoreMatrix::Standard(ScoringMatrix::Pam30),
        BlastpScoreMatrix::Standard(ScoringMatrix::Blosum45),
        BlastpScoreMatrix::Standard(ScoringMatrix::Blosum80),
    ];

    #[test]
    fn x_bpgap_matches_reference_on_random_problems() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut scratch = GapAlignScratch::new();
        let gaps = [
            (11, 1),
            (10, 1),
            (9, 2),
            (5, 2),
            (0, 1),
            (1, 0),
            (0, 0),
            (15, 3),
        ];
        for round in 0..3000 {
            let alphabet = [2u64, 4, 20, 24][rng.below(4) as usize];
            let ql = 1 + rng.below(200) as usize;
            let sl = 1 + rng.below(200) as usize;
            let q = random_seq(&mut rng, ql, alphabet, round % 3 == 0);
            let s = random_seq(&mut rng, sl, alphabet, round % 5 == 0);
            let m = MATRICES[rng.below(4) as usize];
            let (go, ge) = gaps[rng.below(gaps.len() as u64) as usize];
            let x = [0, 5, 12, 16, 38, 60, 150][rng.below(7) as usize];
            let qs = rng.below(ql as u64) as usize;
            let ss = rng.below(sl as u64) as usize;
            check_seed(&q, &s, qs, ss, m, go, ge, x, &mut scratch);
        }
    }

    #[test]
    fn x_bpgap_matches_reference_on_related_sequences() {
        // Long alignments: wide bands, many rows where a gap opens, scores past 255.
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        let mut scratch = GapAlignScratch::new();
        for round in 0..400 {
            let len = 20 + rng.below(if round % 10 == 0 { 3000 } else { 400 }) as usize;
            let q = random_seq(&mut rng, len, 20, round % 4 == 0);
            let s = mutate(&mut rng, &q, [5, 15, 30, 50][round % 4]);
            if s.is_empty() {
                continue;
            }
            let m = MATRICES[round % 4];
            let (go, ge) = [(11, 1), (10, 1), (8, 2)][round % 3];
            let x = [38, 16, 100, 300][rng.below(4) as usize];
            let qs = rng.below(q.len() as u64) as usize;
            let ss = qs.min(s.len() - 1);
            check_seed(&q, &s, qs, ss, m, go, ge, x, &mut scratch);
        }
    }

    #[test]
    fn x_bpgap_matches_reference_on_ties_and_drop_boundaries() {
        // Two letters with scores that tie often; every X-drop from 0 to 40 meets cells exactly
        // at `best_score - x_dropoff` and best scores reached again later in a row.
        let mut rng = Rng(0x1234_5678_9ABC_DEF1);
        let mut scratch = GapAlignScratch::new();
        for x in 0..=40 {
            for _ in 0..40 {
                let ql = 1 + rng.below(80) as usize;
                let sl = 1 + rng.below(80) as usize;
                let q = random_seq(&mut rng, ql, 2, false);
                let s = random_seq(&mut rng, sl, 2, false);
                let qs = rng.below(ql as u64) as usize;
                let ss = rng.below(sl as u64) as usize;
                for (go, ge) in [(11, 1), (3, 1), (2, 2)] {
                    check_seed(&q, &s, qs, ss, MATRICES[0], go, ge, x, &mut scratch);
                }
            }
        }
        // Repeats: the same score on many diagonals.
        let unit = [1u8, 7, 9, 3, 1, 7];
        let q: Vec<u8> = unit.iter().cycle().take(300).copied().collect();
        let s: Vec<u8> = unit.iter().cycle().skip(2).take(280).copied().collect();
        for x in [0, 1, 7, 12, 38, 90] {
            check_seed(&q, &s, 150, 140, MATRICES[0], 11, 1, x, &mut scratch);
            check_seed(&q, &s, 0, 0, MATRICES[0], 11, 1, x, &mut scratch);
            check_seed(&q, &s, 299, 279, MATRICES[0], 11, 1, x, &mut scratch);
        }
    }

    #[test]
    fn x_bpgap_matches_reference_at_the_ends() {
        // Lengths 1-3, seeds at both ends (sentinels read past the window), subject windows that
        // start inside the subject (the reverse call's `subject_base`), and the gap letter.
        let mut rng = Rng(0x0F0F_1E1E_2D2D_3C3C);
        let mut scratch = GapAlignScratch::new();
        for ql in 1..=4usize {
            for sl in 1..=4usize {
                for _ in 0..30 {
                    let q = random_seq(&mut rng, ql, 20, true);
                    let s = random_seq(&mut rng, sl, 20, true);
                    for qs in 0..ql {
                        for ss in 0..sl {
                            check_seed(&q, &s, qs, ss, MATRICES[0], 11, 1, 38, &mut scratch);
                        }
                    }
                }
            }
        }
        for _ in 0..300 {
            let q = random_seq(&mut rng, 60, 20, true);
            let s = random_seq(&mut rng, 90, 20, true);
            let sb = rng.below(30) as usize;
            let len2 = 1 + rng.below((90 - sb - 1) as u64) as usize;
            let len1 = 1 + rng.below(59) as usize;
            for rev in [false, true] {
                let qb = if rev {
                    0
                } else {
                    rng.below((60 - len1) as u64) as usize
                };
                let qb = qb.min(59 - len1.min(59));
                check(
                    &q,
                    qb,
                    &s,
                    sb,
                    len1,
                    len2,
                    MATRICES[0],
                    11,
                    1,
                    38,
                    rev,
                    &mut scratch,
                );
            }
        }
        // Empty problems.
        let q = [1u8, 2, 3];
        assert_eq!(
            check(
                &q,
                0,
                &q,
                0,
                0,
                2,
                MATRICES[0],
                11,
                1,
                38,
                false,
                &mut scratch
            ),
            (0, 0, 0)
        );
        assert_eq!(
            check(
                &q,
                0,
                &q,
                0,
                2,
                0,
                MATRICES[0],
                11,
                1,
                38,
                true,
                &mut scratch
            ),
            (0, 0, 0)
        );
    }

    #[test]
    fn x_bpgap_reused_scratch_matches_fresh_scratch() {
        // One scratch over problems of very different sizes (stale cells from a longer problem
        // must not leak into a shorter one) against a fresh scratch per problem.
        let mut rng = Rng(0xAAAA_5555_3333_CCCC);
        let mut shared = GapAlignScratch::new();
        for round in 0..300 {
            let len = [5usize, 2000, 30, 700, 1][round % 5] + rng.below(20) as usize;
            let q = random_seq(&mut rng, len, 20, false);
            let s = mutate(&mut rng, &q, 10);
            if s.is_empty() {
                continue;
            }
            let qs = rng.below(q.len() as u64) as usize;
            let ss = qs.min(s.len() - 1);
            let x = [16, 38, 200][round % 3];
            check_seed(&q, &s, qs, ss, MATRICES[0], 11, 1, x, &mut shared);
            let mut fresh = GapAlignScratch::new();
            check_seed(&q, &s, qs, ss, MATRICES[0], 11, 1, x, &mut fresh);
        }
    }

    #[test]
    fn x_bpgap_column_letters_are_the_reference_reads() {
        // Every window position, including bases and lengths that reach past the subject.
        let mut rng = Rng(0x7777_1111_2222_9999);
        let mut col = vec![0u8; 200];
        for _ in 0..4000 {
            let sl = rng.below(40) as usize;
            let subject = random_seq(&mut rng, sl, 27, false);
            let subject_base = rng.below(50) as usize;
            let len2 = rng.below(60) as usize;
            let from = rng.below(len2 as u64 + 30) as usize;
            let to = from + rng.below(100) as usize;
            for reverse in [false, true] {
                let ok = fill_col(&mut col, from, to, &subject, subject_base, len2, reverse);
                for (c, &v) in col.iter().enumerate().take(to).skip(from) {
                    let want = if c == 0 || c > len2 + 1 {
                        0
                    } else if reverse {
                        blastp_subject_residue::<true>(&subject, subject_base, len2, c - 1)
                    } else {
                        blastp_subject_residue::<false>(&subject, subject_base, len2, c - 1)
                    };
                    assert_eq!(
                        v, want,
                        "c={c} base={subject_base} len2={len2} rev={reverse}"
                    );
                }
                assert!(ok);
            }
        }
    }

    #[test]
    fn x_bpgap_leaves_unusual_problems_to_the_reference() {
        let mut scratch = GapAlignScratch::new();
        let q = [1u8, 2, 3, 4];
        let s = [1u8, 2, 30, 4];
        // A letter code of 28 or more.
        assert!(restricted_fast(
            &q,
            0,
            &s,
            0,
            3,
            3,
            MATRICES[0],
            11,
            1,
            38,
            false,
            &mut scratch.x_restricted,
            false
        )
        .is_none());
        // Gap costs and X-drop outside the checked ranges.
        for (go, ge, x) in [
            (-1, 1, 38),
            (11, -1, 38),
            (MAX_GAP_COST + 1, 1, 38),
            (11, 1, MAX_X_DROP + 1),
        ] {
            assert!(restricted_fast(
                &q,
                0,
                &q,
                0,
                3,
                3,
                MATRICES[0],
                go,
                ge,
                x,
                false,
                &mut scratch.x_restricted,
                false
            )
            .is_none());
        }
    }

    /// Timing aid (ignored; run with `cargo test --release -- --ignored x_bpgap_time_kernels
    /// --nocapture`): the reference, the scalar rows and the AVX2 rows on the same unrelated
    /// sequence pairs (the restricted alignments of the preliminary stage are low-scoring ones),
    /// one after the other, three times. Prints ns per call and per row.
    #[test]
    #[ignore]
    fn x_bpgap_time_kernels() {
        let mut rng = Rng(0x5151_7272_9393_B4B4);
        let mut problems = Vec::new();
        for _ in 0..4000 {
            let (ql, sl) = (100 + rng.below(500) as usize, 100 + rng.below(500) as usize);
            let q = random_seq(&mut rng, ql, 20, false);
            let s = random_seq(&mut rng, sl, 20, false);
            let qs = rng.below(q.len() as u64) as usize;
            let ss = rng.below(s.len() as u64) as usize;
            problems.push((q, s, qs, ss));
        }
        let mut scratch = GapAlignScratch::new();
        let m = MATRICES[0];
        let calls = |scratch: &mut GapAlignScratch, kernel: u8| -> (u64, i64) {
            let mut sum = 0i64;
            let mut n = 0u64;
            for (q, s, qs, ss) in &problems {
                let (qs, ss) = (*qs, *ss);
                for (qb, sb, len1, len2, rev) in [
                    (0, 0, qs + 1, ss + 1, true),
                    (qs, ss, q.len() - qs - 1, s.len() - ss - 1, false),
                ] {
                    let r = if kernel == 0 {
                        reference(q, qb, s, sb, len1, len2, m, 11, 1, 38, rev, scratch)
                    } else {
                        restricted_fast(
                            q,
                            qb,
                            s,
                            sb,
                            len1,
                            len2,
                            m,
                            11,
                            1,
                            38,
                            rev,
                            &mut scratch.x_restricted,
                            kernel == 2,
                        )
                        .unwrap()
                    };
                    sum += r.2 as i64 + r.0 as i64;
                    n += 1;
                }
            }
            (n, sum)
        };
        let kernels: &[u8] = if cpu_has_avx2() { &[0, 1, 2] } else { &[0, 1] };
        for round in 0..3 {
            for &k in kernels {
                let t = std::time::Instant::now();
                let (n, sum) = calls(&mut scratch, k);
                let ns = t.elapsed().as_nanos() as f64;
                eprintln!(
                    "round {round} kernel {} calls {n} {:.1} ns/call checksum {sum}",
                    ["reference", "scalar", "avx2"][k as usize],
                    ns / n as f64
                );
            }
        }
    }
}
