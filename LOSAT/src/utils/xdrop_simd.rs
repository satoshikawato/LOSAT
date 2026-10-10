//! EXPERIMENT (LOSAT_X_DPFAST / LOSAT_X_DPSHADOW): row-wise SIMD X-drop dynamic programming
//! that is intended to return exactly what the scalar transliterations of
//! NCBI `Blast_SemiGappedAlign` / `ALIGN_EX` return.
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:563-635
//! ```c
//! for (b_index = first_b_index; b_index < b_size; b_index++) {
//!     score_gap_col = score_array[b_index].best_gap;
//!     next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
//!     if (score < score_gap_col) { script = SCRIPT_GAP_IN_B; score = score_gap_col; }
//!     if (score < score_gap_row) { script = SCRIPT_GAP_IN_A; score = score_gap_row; }
//!     if (best_score - score > x_dropoff) {
//!         if (first_b_index == b_index) first_b_index++;
//!         else score_array[b_index].best = MININT;
//!     } else { ... }
//! }
//! ```
//!
//! Why the result is the same. Call `thr = best_score - x_dropoff` at the
//! moment a cell is evaluated; `thr` never decreases. A cell that survives has
//! `score >= thr`, so its score is decided by candidates that are themselves
//! `>= thr`. Everything this kernel computes differently from the scalar loop
//! (the gap states behind a dropped cell, saturated values, `-inf + score`)
//! is `< thr` in both computations, so it can neither win a `max` of a
//! surviving cell nor change which cells are dropped, and the traceback only
//! walks surviving cells.
//!
//! Representation: stored value = `score - base + bias`, a non-negative
//! 16-bit number, 0 standing for minus infinity. `bias` is chosen so that
//! every value `>= thr` stays exact after one gap cost is subtracted and so
//! that `0 + max_score < thr`.
//!
//! Per row, blocks of 8 lanes, left to right:
//!   D  = H_prev[b-1] + score            (diagonal)
//!   H0 = max(D, E)
//!   F  = decayed prefix maximum of H0 - (go+ge)   (row gap)
//!   H  = max(H0, F)
//!   dropped = H < running best - X      (dropped cells store 0)
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770,862-874,903-908
//! ```c
//!     if (!score_only) {
//!         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
//! ...
//!         for (b_index = first_b_index; b_index < b_size; b_index++) {
//! ...
//!             if (score < score_gap_col)
//!                 score = score_gap_col;
//!             if (score < score_gap_row)
//!                 score = score_gap_row;
//!             if (best_score - score > x_dropoff) {
//! ...
//!                 score_gap_row -= gap_extend;
//!                 score_gap_col -= gap_extend;
//!                 score_array[b_index].best_gap = MAX(score - gap_open_extend,
//!                                                     score_gap_col);
//!                 score_gap_row = MAX(score - gap_open_extend, score_gap_row);
//!                 score_array[b_index].best = score;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3125-3134,3146-3158,3187-3193
//! ```c
//!         if(reverse_sequence) {
//!             a_base_pair = NCBI2NA_UNPACK_BASE(A[(M-a_index)/4],
//!                                                ((a_index-1)%4));
//!             matrix_row = matrix[a_base_pair];
//! ...
//!         for (b_index = first_b_index; b_index < b_size; b_index++) {
//! ...
//!             if (best_score - score > x_dropoff) {
//! ...
//!                 score_gap_row = MAX(score - gap_open_extend, score_gap_row);
//!                 score_array[b_index].best = score;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:476-500,563-576,601-636,682-727
//! ```c
//!     score_array[0].best = 0;
//!     score_array[0].best_gap = -gap_open_extend;
//! ...
//!     for (a_index = 1; a_index <= M; a_index++) {
//! ...
//!         for (b_index = first_b_index; b_index < b_size; b_index++) {
//! ...
//!             if (matrix_index == FENCE_SENTRY) {
//!                 if (fence_hit) {
//!                     *fence_hit = 1;
//!                 }
//!                 break;
//!             }
//! ...
//!             if (best_score - score > x_dropoff) {
//! ...
//!             edit_script_row[b_index] = script;
//! ...
//!     while (a_index > 0 || b_index > 0) {
//! ...
//!         GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
//! ```
//!
//! Which NCBI code this module serves. The protein programs (blastp, tblastn, blastx)
//! use `Blast_SemiGappedAlign` with `score_only` TRUE (preliminary gapped score) and with
//! `score_only` FALSE, which calls `ALIGN_EX` (traceback). blastn uses
//! `s_BlastAlignPackedNucl` (preliminary score on the packed subject) and `ALIGN_EX` on
//! one byte per base (traceback). LOSAT_X_DPFAST runs the vector kernels below for all
//! of these.
//!
//! How the module relates to that code. Each row of the DP is the same recurrence as the
//! cited loops. The diagonal term (`next_score`), the column gap (`best_gap`) and the
//! X-drop test are evaluated lane by lane, in blocks of 8 (16-bit) or 16 (8-bit) cells.
//! The row gap (`score_gap_row`), which the C loop carries from cell to cell, is a decayed
//! prefix maximum over the lanes. The script bits have the values of the NCBI `SCRIPT_*`
//! constants, and the traceback walk is the C walk. The argument above is why the cells that
//! are kept have the C values. LOSAT_X_DPSHADOW and the random test
//! (`simd_kernels_match_scalar_loops_on_random_problems`) compare every result with the
//! scalar loop. The module is integer arithmetic only. A problem that does not fit the lanes
//! is not handled here and the scalar port runs.

#![allow(clippy::too_many_arguments)]

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:362-371
// ```c
// /** Values for the editing script operations in traceback */
// enum {
//     SCRIPT_SUB           = eGapAlignSub,     /**< Substitution */
//     SCRIPT_GAP_IN_A      = eGapAlignDel,     /**< Deletion */
//     SCRIPT_GAP_IN_B      = eGapAlignIns,     /**< Insertion */
//     SCRIPT_OP_MASK       = 0x07, /**< Mask for edit script operations */
// ...
//     SCRIPT_EXTEND_GAP_A  = 0x10, /**< continue a gap in A */
//     SCRIPT_EXTEND_GAP_B  = 0x40  /**< continue a gap in B */
// };
// ```
// The same script bits as NCBI. The three operation values are NCBI's `eGapAlignSub` (3),
// `eGapAlignDel` (0) and `eGapAlignIns` (6) from `gapinfo.h`.
pub(crate) const SCRIPT_SUB: u8 = 3;
pub(crate) const SCRIPT_GAP_IN_A: u8 = 0;
pub(crate) const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:364-364
// ```c
// #define FENCE_SENTRY 201
// ```
// The same sentinel value as NCBI. The kernels stop at it as `ALIGN_EX` does (lines 571-575 of
// `blast_gapalign.c`).
const FENCE_SENTRY: u8 = 201;

// No NCBI counterpart: lane padding, 16-bit value range and the largest score the vector
// path accepts. A problem outside these limits is not handled and the scalar port runs; the
// limits do not change any value NCBI computes.
const PAD: usize = 32;
const VMAX: i64 = 32767;
/// Largest substitution score the 16-bit path accepts.
const MAX_M: i32 = 2000;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:541-545,3125-3134
/// ```c
///         if (!(gap_align->positionBased)) {
///             if(reverse_sequence)
///                 matrix_row = matrix[ A[ M - a_index ] ];
///             else
///                 matrix_row = matrix[ A[ a_index ] ];
/// ...
///             a_base_pair = NCBI2NA_UNPACK_BASE(A[1+((a_index-1)/4)],
///                                                (3-((a_index-1)%4)));
///             matrix_row = matrix[a_base_pair];
/// ```
/// The residue of `A` (the row sequence) that selects the matrix row of row `a_index`, one byte
/// per residue or unpacked from NCBI2NA bytes, forward or reversed as in the C code.
///
/// Row residues: `data[base + step * a]` for `a` in `1..=len1`.
#[derive(Clone, Copy)]
pub(crate) enum RowSeq<'a> {
    Bytes {
        data: &'a [u8],
        base: isize,
        step: isize,
    },
    /// NCBI2NA packed rows, forward (`s_BlastAlignPackedNucl`).
    PackedFwd { data: &'a [u8], byte_offset: isize },
    /// NCBI2NA packed rows, reverse.
    PackedRev { data: &'a [u8], len: usize },
}

impl RowSeq<'_> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3125-3134
    /// ```c
    ///         if(reverse_sequence) {
    ///             a_base_pair = NCBI2NA_UNPACK_BASE(A[(M-a_index)/4],
    ///                                                ((a_index-1)%4));
    ///             matrix_row = matrix[a_base_pair];
    ///         }
    ///         else {
    ///             a_base_pair = NCBI2NA_UNPACK_BASE(A[1+((a_index-1)/4)],
    ///                                                (3-((a_index-1)%4)));
    /// ```
    /// The same index arithmetic for the packed forms.
    #[inline(always)]
    fn get(&self, a: usize) -> Option<u8> {
        match *self {
            RowSeq::Bytes { data, base, step } => {
                let i = base + step * a as isize;
                if i < 0 {
                    return None;
                }
                data.get(i as usize).copied()
            }
            RowSeq::PackedFwd { data, byte_offset } => {
                let byte_index = byte_offset + 1 + ((a - 1) / 4) as isize;
                let shift = 3 - ((a - 1) % 4);
                let byte = if byte_index >= 0 {
                    data.get(byte_index as usize).copied().unwrap_or(0)
                } else {
                    0
                };
                Some((byte >> (2 * shift)) & 0x03)
            }
            RowSeq::PackedRev { data, len } => {
                let byte_index = (len - a) / 4;
                let shift = (a - 1) % 4;
                let byte = data.get(byte_index).copied().unwrap_or(0);
                Some((byte >> (2 * shift)) & 0x03)
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:554-557,566-576
/// ```c
///         if(reverse_sequence)
///             b_ptr = &B[N - first_b_index];
///         else
///             b_ptr = &B[first_b_index];
/// ...
///             b_ptr += b_increment;
/// ...
///             matrix_index = *b_ptr;
///
///             if (matrix_index == FENCE_SENTRY) {
/// ```
/// The residue of `B` (the column sequence) at column `k`, read forward or backward as `b_ptr`
/// is in the C loop. Out of range reads 0; the fence sentinel is recognised by the kernel.
///
/// Column residues: column `k` reads `data[base + step * k]`; out of range
/// or `k >= zero_from` reads 0.
#[derive(Clone, Copy)]
pub(crate) struct ColSeq<'a> {
    pub data: &'a [u8],
    pub base: isize,
    pub step: isize,
    pub zero_from: usize,
}

impl ColSeq<'_> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:566-576
    /// ```c
    ///             b_ptr += b_increment;
    /// ...
    ///             matrix_index = *b_ptr;
    ///
    ///             if (matrix_index == FENCE_SENTRY) {
    ///                 if (fence_hit) {
    ///                     *fence_hit = 1;
    ///                 }
    ///                 break;
    /// ```
    /// One column residue (`*b_ptr`).
    #[inline(always)]
    fn get(&self, k: usize) -> u8 {
        if k >= self.zero_from {
            return 0;
        }
        let i = self.base + self.step * k as isize;
        if i < 0 {
            return 0;
        }
        self.data.get(i as usize).copied().unwrap_or(0)
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,569,578
/// ```c
///                 matrix_row = matrix[ A[ M - a_index ] ];
/// ...
///             matrix_index = *b_ptr;
/// ...
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ```
/// `matrix[row residue][column residue]`: the score matrix NCBI reads through `matrix_row`.
///
/// Substitution scores by (row residue, column residue).
#[derive(Clone, Copy)]
pub(crate) enum Scores<'a> {
    /// 28x28 (NCBIstdaa) rows.
    Rows28(&'a [[i32; 28]; 28]),
    /// Flat `n x n` matrix (BLASTNA: n = 16).
    Flat { data: &'a [i32], n: usize },
    /// Arbitrary function; must be defined for residues `< n`.
    Func {
        f: &'a dyn Fn(u8, u8) -> i32,
        n: usize,
    },
    /// Prebuilt lookup tables for an `n x n` matrix, e.g. BLOSUM62.
    Tables { t: &'a StaticTables, n: usize },
}

// No NCBI counterpart: the layout of the score tables (bytes per row residue) used by the
// vector lookups; the tables hold the same matrix values.
/// Bytes per row residue in a lookup table: low and high bytes of the 16-bit
/// score, then the score clamped to 8 bits.
const TAB: usize = 96;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
// ```c
//                 matrix_row = matrix[ A[ M - a_index ] ];
// ...
//             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
// ```
// The matrix rows of a fixed matrix (for example BLOSUM62) copied into lookup tables for the
// vector kernels, built once. The values are the matrix entries.
pub(crate) struct StaticTables {
    data: [u8; 32 * TAB],
    /// Largest score in the matrix.
    max: i32,
}

impl Scores<'_> {
    /// Residues `>= ncols()` are not handled by the SIMD kernel.
    #[inline]
    fn ncols(&self) -> usize {
        match *self {
            Scores::Rows28(_) => 28,
            Scores::Flat { n, .. } | Scores::Func { n, .. } | Scores::Tables { n, .. } => n.min(32),
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
/// ```c
///                 matrix_row = matrix[ A[ M - a_index ] ];
/// ...
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ```
/// Copies one row of the score matrix (the `matrix_row` of one row residue) into the lookup
/// planes. It returns false, and the scalar port runs, if a score is outside the range of the
/// 16-bit path.
///
/// Build the per-row-residue lookup planes (32 low bytes, 32 high bytes, 32
/// 8-bit scores). Returns false when a score is outside the 16-bit path's range.
fn build_row_table(scores: &Scores<'_>, q: usize, out: &mut [u8]) -> bool {
    let mut tmp = [0i32; 32];
    match *scores {
        Scores::Rows28(m) => {
            if q >= 28 {
                return false;
            }
            tmp[..28].copy_from_slice(&m[q]);
        }
        Scores::Flat { data, n } => {
            if q >= n || n > 32 {
                return false;
            }
            tmp[..n].copy_from_slice(&data[q * n..q * n + n]);
        }
        Scores::Func { f, n } => {
            if q >= n || n > 32 {
                return false;
            }
            for (r, slot) in tmp.iter_mut().enumerate().take(n) {
                *slot = f(q as u8, r as u8);
            }
        }
        Scores::Tables { .. } => return false,
    }
    for r in 0..32 {
        let v = tmp[r];
        if v > MAX_M {
            return false;
        }
        // Any score <= -32768 puts the diagonal below every threshold of this
        // path (X < 32768) in both the exact and the saturated computation,
        // so "minus infinity" entries (BLASTNA gap column, COMPO_SCORE_MIN)
        // can be clamped.
        let v = v.max(i16::MIN as i32);
        let u = v as i16 as u16;
        out[r] = (u & 0xFF) as u8;
        out[32 + r] = (u >> 8) as u8;
        out[64 + r] = v.clamp(-128, 127) as i8 as u8;
    }
    true
}

/// No NCBI counterpart: the largest matrix score, used to decide whether the 8-bit lanes can
/// hold the problem; it does not change any value NCBI computes.
///
/// Largest score of the matrix, when it is cheap to know.
fn max_score(scores: &Scores<'_>) -> Option<i32> {
    match *scores {
        Scores::Tables { t, .. } => Some(t.max),
        Scores::Flat { data, n } => data[..n * n].iter().copied().max(),
        Scores::Rows28(m) => m.iter().flat_map(|row| row.iter().copied()).max(),
        Scores::Func { .. } => None,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:543,578
/// ```c
///                 matrix_row = matrix[ A[ M - a_index ] ];
/// ...
///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
/// ```
/// Builds the lookup planes of every matrix row; see `build_row_table`.
///
/// Build all row tables for a fixed matrix (used once per matrix).
pub(crate) fn build_tables(f: &dyn Fn(u8, u8) -> i32, n: usize) -> Option<Box<StaticTables>> {
    let mut t = Box::new(StaticTables {
        data: [0u8; 32 * TAB],
        max: i32::MIN,
    });
    let scores = Scores::Func { f, n };
    for q in 0..n.min(32) {
        if !build_row_table(&scores, q, &mut t.data[q * TAB..(q + 1) * TAB]) {
            return None;
        }
        for r in 0..n.min(32) {
            t.max = t.max.max(f(q as u8, r as u8));
        }
    }
    Some(t)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:462-468,476-479
// ```c
//     if (num_extra_cells > gap_align->dp_mem_alloc) {
//         gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
//                                       2 * gap_align->dp_mem_alloc);
// ...
//     score_array = gap_align->dp_mem;
//     score_array[0].best = 0;
//     score_array[0].best_gap = -gap_open_extend;
// ```
// The vector kernels' rows (`best`, `best_gap` of `score_array`, here `H` and `E`), the traceback
// rows (`edit_script`) and the operation list. Kept between calls like `gap_align->dp_mem`.
pub(crate) struct XdropScratch {
    ha: Vec<i16>,
    hb: Vec<i16>,
    e: Vec<i16>,
    // 8-bit path
    h8a: Vec<u8>,
    h8b: Vec<u8>,
    e8: Vec<u8>,
    dirty8: usize,
    sb: Vec<u8>,
    sb_filled: usize,
    fence_k: usize,
    bad_k: usize,
    tabs: Vec<u8>,
    tab_ok: u32,
    dirty: usize,
    trace: Vec<u8>,
    rows: Vec<(u32, u32, u32)>,
    surv: Vec<u64>,
    /// Traceback operations in walk order (alignment end -> start), run-length encoded.
    pub ops: Vec<(u8, u32)>,
}

impl XdropScratch {
    /// No NCBI counterpart: an empty scratch when the kernel is off, so that a search without the
    /// switch allocates nothing for it; it does not change any value NCBI computes.
    ///
    /// The scratch of a search: empty unless the kernel is switched on, so
    /// that a search without the switch allocates nothing for it.
    pub(crate) fn for_search() -> Self {
        if mode() != 0 {
            return Self::new();
        }
        Self {
            ha: Vec::new(),
            hb: Vec::new(),
            e: Vec::new(),
            h8a: Vec::new(),
            h8b: Vec::new(),
            e8: Vec::new(),
            dirty8: 0,
            sb: Vec::new(),
            sb_filled: 0,
            fence_k: usize::MAX,
            bad_k: usize::MAX,
            tabs: Vec::new(),
            tab_ok: 0,
            dirty: 0,
            trace: Vec::new(),
            rows: Vec::new(),
            surv: Vec::new(),
            ops: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:462-468
    // ```c
    //     if (num_extra_cells > gap_align->dp_mem_alloc) {
    //         gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
    //                                       2 * gap_align->dp_mem_alloc);
    //         sfree(gap_align->dp_mem);
    // ```
    // The initial allocation, like the first `dp_mem` allocation.
    pub(crate) fn new() -> Self {
        Self {
            ha: vec![0; 1024],
            hb: vec![0; 1024],
            e: vec![0; 1024],
            h8a: vec![0; 1024],
            h8b: vec![0; 1024],
            e8: vec![0; 1024],
            dirty8: 0,
            sb: vec![0; 1024],
            sb_filled: 0,
            fence_k: usize::MAX,
            bad_k: usize::MAX,
            tabs: vec![0; 32 * TAB],
            tab_ok: 0,
            dirty: 0,
            trace: Vec::new(),
            rows: Vec::new(),
            surv: vec![0; 64],
            ops: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:641-649
    // ```c
    //         if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {
    //
    //             gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
    //                                           2 * gap_align->dp_mem_alloc);
    //             score_array = (BlastGapDP *)realloc(score_array,
    //                                                gap_align->dp_mem_alloc *
    //                                                sizeof(BlastGapDP));
    // ```
    // Grows the rows when the band needs more room, as the C code enlarges `score_array`.
    #[inline(always)]
    fn ensure(&mut self, lanes: usize) {
        let need = lanes + 2 * PAD + 64;
        if self.ha.len() < need {
            self.grow(need);
        }
    }

    #[cold]
    #[inline(never)]
    fn grow(&mut self, need: usize) {
        {
            let n = need.max(self.ha.len() * 2);
            self.ha.resize(n, 0);
            self.hb.resize(n, 0);
            self.e.resize(n, 0);
            self.h8a.resize(n, 0);
            self.h8b.resize(n, 0);
            self.e8.resize(n, 0);
            self.sb.resize(n, 0);
            self.surv.resize(n / 32 + 8, 0);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:438,476-479
    // ```c
    //     s_GapPurgeState(gap_align->state_struct);
    // ...
    //     score = -gap_open_extend;
    //     score_array = gap_align->dp_mem;
    //     score_array[0].best = 0;
    //     score_array[0].best_gap = -gap_open_extend;
    // ```
    // Clears what the previous call wrote, so that the next call starts from empty rows.
    fn reset(&mut self) {
        let hi = (self.dirty + 2 * PAD + 32).min(self.ha.len());
        self.ha[..hi].fill(0);
        self.hb[..hi].fill(0);
        self.e[..hi].fill(0);
        self.dirty = 0;
        let hi8 = (self.dirty8 + 2 * PAD + 32).min(self.h8a.len());
        self.h8a[..hi8].fill(0);
        self.h8b[..hi8].fill(0);
        self.e8[..hi8].fill(0);
        self.dirty8 = 0;
        self.sb_filled = 0;
        self.fence_k = usize::MAX;
        self.bad_k = usize::MAX;
        self.tab_ok = 0;
        self.trace.clear();
        self.rows.clear();
        self.ops.clear();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:612-613,732
// ```c
//                     *a_offset = a_index;
//                     *b_offset = b_index;
// ...
//     return best_score;
// ```
// The output of `ALIGN_EX` / `Blast_SemiGappedAlign`: `a_offset`, `b_offset`, the best score, and
// whether the fence sentinel was hit. `cells` only counts work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct XdropResult {
    pub a_offset: usize,
    pub b_offset: usize,
    pub score: i32,
    pub fence_hit: bool,
    /// Band cells the scalar loop visits (rows that run to completion).
    pub cells: u64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:726
// ```c
//         GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
// ```
// Appends one operation to the edit script, merging equal neighbours as `GapPrelimEditBlockAdd`
// does.
#[inline(always)]
fn push_op(ops: &mut Vec<(u8, u32)>, op: u8) {
    if let Some(last) = ops.last_mut() {
        if last.0 == op {
            last.1 += 1;
            return;
        }
    }
    ops.push((op, 1));
}

// ---------------------------------------------------------------------------
// 128-bit vector layer
// ---------------------------------------------------------------------------

// No NCBI counterpart: vector instruction wrappers (one operation on all lanes); they do not
// change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
mod arch {
    use core::arch::x86_64::*;
    pub type V = __m128i;
    #[inline(always)]
    pub unsafe fn zero() -> V {
        _mm_setzero_si128()
    }
    #[inline(always)]
    pub unsafe fn splat16(x: i16) -> V {
        _mm_set1_epi16(x)
    }
    #[inline(always)]
    pub unsafe fn splat8(x: i8) -> V {
        _mm_set1_epi8(x)
    }
    #[inline(always)]
    pub unsafe fn load<T>(p: *const T) -> V {
        _mm_loadu_si128(p as *const V)
    }
    #[inline(always)]
    pub unsafe fn store<T>(p: *mut T, v: V) {
        _mm_storeu_si128(p as *mut V, v)
    }
    #[inline(always)]
    pub unsafe fn adds16(a: V, b: V) -> V {
        _mm_adds_epi16(a, b)
    }
    #[inline(always)]
    pub unsafe fn add16(a: V, b: V) -> V {
        _mm_add_epi16(a, b)
    }
    #[inline(always)]
    pub unsafe fn subsu16(a: V, b: V) -> V {
        _mm_subs_epu16(a, b)
    }
    #[inline(always)]
    pub unsafe fn max16(a: V, b: V) -> V {
        _mm_max_epi16(a, b)
    }
    #[inline(always)]
    pub unsafe fn gt16(a: V, b: V) -> V {
        _mm_cmpgt_epi16(a, b)
    }
    #[inline(always)]
    pub unsafe fn and(a: V, b: V) -> V {
        _mm_and_si128(a, b)
    }
    #[inline(always)]
    pub unsafe fn or(a: V, b: V) -> V {
        _mm_or_si128(a, b)
    }
    /// `!mask & v`
    #[inline(always)]
    pub unsafe fn andnot(mask: V, v: V) -> V {
        _mm_andnot_si128(mask, v)
    }
    #[inline(always)]
    pub unsafe fn shlq16(v: V) -> V {
        _mm_slli_epi64(v, 16)
    }
    #[inline(always)]
    pub unsafe fn shlq32(v: V) -> V {
        _mm_slli_epi64(v, 32)
    }
    #[inline(always)]
    pub unsafe fn bsl2(v: V) -> V {
        _mm_slli_si128(v, 2)
    }
    #[inline(always)]
    pub unsafe fn bsl4(v: V) -> V {
        _mm_slli_si128(v, 4)
    }
    #[inline(always)]
    pub unsafe fn bsl8(v: V) -> V {
        _mm_slli_si128(v, 8)
    }
    /// Byte shuffle; an index byte with the top bit set gives 0.
    #[inline(always)]
    pub unsafe fn pick(v: V, idx: V) -> V {
        _mm_shuffle_epi8(v, idx)
    }
    /// 32-entry byte table lookup (indices >= 32 give an unspecified table byte or 0).
    #[inline(always)]
    pub unsafe fn lut32(t0: V, t1: V, idx: V) -> V {
        let i0 = _mm_add_epi8(idx, _mm_set1_epi8(0x70));
        let i1 = _mm_sub_epi8(idx, _mm_set1_epi8(0x10));
        _mm_or_si128(_mm_shuffle_epi8(t0, i0), _mm_shuffle_epi8(t1, i1))
    }
    /// Lane 7 (16-bit) to every lane.
    #[inline(always)]
    pub unsafe fn bcast_last16(v: V, mask: V) -> V {
        _mm_shuffle_epi8(v, mask)
    }
    /// Lane 3 (16-bit) to lanes 4..8, zero elsewhere.
    #[inline(always)]
    pub unsafe fn cross16(v: V, mask: V) -> V {
        _mm_shuffle_epi8(v, mask)
    }
    /// Lane 15 (8-bit) to every lane.
    #[inline(always)]
    pub unsafe fn bcast_last8(v: V, mask: V) -> V {
        _mm_shuffle_epi8(v, mask)
    }
    /// Lane 7 (8-bit) to lanes 8..16, zero elsewhere.
    #[inline(always)]
    pub unsafe fn cross8(v: V, mask: V) -> V {
        _mm_shuffle_epi8(v, mask)
    }
    #[inline(always)]
    pub unsafe fn unpacklo8(a: V, b: V) -> V {
        _mm_unpacklo_epi8(a, b)
    }
    #[inline(always)]
    pub unsafe fn unpackhi8(a: V, b: V) -> V {
        _mm_unpackhi_epi8(a, b)
    }
    #[inline(always)]
    pub unsafe fn packus16(a: V, b: V) -> V {
        _mm_packus_epi16(a, b)
    }
    #[inline(always)]
    pub unsafe fn mask8(v: V) -> u32 {
        _mm_movemask_epi8(v) as u32
    }
    #[inline(always)]
    pub unsafe fn lane0(v: V) -> i32 {
        _mm_cvtsi128_si32(v) as i16 as i32
    }
    #[inline(always)]
    pub unsafe fn adds8(a: V, b: V) -> V {
        _mm_adds_epi8(a, b)
    }
    #[inline(always)]
    pub unsafe fn add8(a: V, b: V) -> V {
        _mm_add_epi8(a, b)
    }
    #[inline(always)]
    pub unsafe fn subsu8(a: V, b: V) -> V {
        _mm_subs_epu8(a, b)
    }
    #[inline(always)]
    pub unsafe fn max8(a: V, b: V) -> V {
        _mm_max_epi8(a, b)
    }
    #[inline(always)]
    pub unsafe fn gt8(a: V, b: V) -> V {
        _mm_cmpgt_epi8(a, b)
    }
    #[inline(always)]
    pub unsafe fn shlq8(v: V) -> V {
        _mm_slli_epi64(v, 8)
    }
    #[inline(always)]
    pub unsafe fn bsl1(v: V) -> V {
        _mm_slli_si128(v, 1)
    }
    #[inline(always)]
    pub unsafe fn lane0_u8(v: V) -> i32 {
        _mm_cvtsi128_si32(v) & 0xFF
    }
    #[inline(always)]
    pub unsafe fn lane15_u8(v: V) -> i32 {
        _mm_extract_epi8(v, 15) & 0xFF
    }
    #[inline(always)]
    pub unsafe fn lane7(v: V) -> i32 {
        _mm_extract_epi16(v, 7) as i16 as i32
    }
}

// No NCBI counterpart: vector instruction wrappers (one operation on all lanes); they do not
// change any value NCBI computes.
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
mod arch {
    use core::arch::wasm32::*;
    pub type V = v128;
    #[inline(always)]
    pub unsafe fn zero() -> V {
        i16x8_splat(0)
    }
    #[inline(always)]
    pub unsafe fn splat16(x: i16) -> V {
        i16x8_splat(x)
    }
    #[inline(always)]
    pub unsafe fn splat8(x: i8) -> V {
        i8x16_splat(x)
    }
    #[inline(always)]
    pub unsafe fn load<T>(p: *const T) -> V {
        v128_load(p as *const v128)
    }
    #[inline(always)]
    pub unsafe fn store<T>(p: *mut T, v: V) {
        v128_store(p as *mut v128, v)
    }
    #[inline(always)]
    pub unsafe fn adds16(a: V, b: V) -> V {
        i16x8_add_sat(a, b)
    }
    #[inline(always)]
    pub unsafe fn add16(a: V, b: V) -> V {
        i16x8_add(a, b)
    }
    #[inline(always)]
    pub unsafe fn subsu16(a: V, b: V) -> V {
        u16x8_sub_sat(a, b)
    }
    #[inline(always)]
    pub unsafe fn max16(a: V, b: V) -> V {
        i16x8_max(a, b)
    }
    #[inline(always)]
    pub unsafe fn gt16(a: V, b: V) -> V {
        i16x8_gt(a, b)
    }
    #[inline(always)]
    pub unsafe fn and(a: V, b: V) -> V {
        v128_and(a, b)
    }
    #[inline(always)]
    pub unsafe fn or(a: V, b: V) -> V {
        v128_or(a, b)
    }
    /// `!mask & v`
    #[inline(always)]
    pub unsafe fn andnot(mask: V, v: V) -> V {
        v128_andnot(v, mask)
    }
    #[inline(always)]
    pub unsafe fn shlq16(v: V) -> V {
        i64x2_shl(v, 16)
    }
    #[inline(always)]
    pub unsafe fn shlq32(v: V) -> V {
        i64x2_shl(v, 32)
    }
    #[inline(always)]
    pub unsafe fn bsl2(v: V) -> V {
        i8x16_shuffle::<16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13>(v, i16x8_splat(0))
    }
    #[inline(always)]
    pub unsafe fn bsl4(v: V) -> V {
        i8x16_shuffle::<16, 17, 18, 19, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11>(v, i16x8_splat(0))
    }
    #[inline(always)]
    pub unsafe fn bsl8(v: V) -> V {
        i8x16_shuffle::<16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7>(v, i16x8_splat(0))
    }
    /// Byte shuffle; an index byte >= 16 gives 0.
    #[inline(always)]
    pub unsafe fn pick(v: V, idx: V) -> V {
        i8x16_swizzle(v, idx)
    }
    #[inline(always)]
    pub unsafe fn lut32(t0: V, t1: V, idx: V) -> V {
        v128_or(
            i8x16_swizzle(t0, idx),
            i8x16_swizzle(t1, i8x16_sub(idx, i8x16_splat(16))),
        )
    }
    /// Lane 7 (16-bit) to every lane (constant shuffle; `mask` unused).
    #[inline(always)]
    pub unsafe fn bcast_last16(v: V, _mask: V) -> V {
        i8x16_shuffle::<14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15>(v, v)
    }
    /// Lane 3 (16-bit) to lanes 4..8, zero elsewhere.
    #[inline(always)]
    pub unsafe fn cross16(v: V, _mask: V) -> V {
        i8x16_shuffle::<16, 16, 16, 16, 16, 16, 16, 16, 6, 7, 6, 7, 6, 7, 6, 7>(v, i16x8_splat(0))
    }
    /// Lane 15 (8-bit) to every lane.
    #[inline(always)]
    pub unsafe fn bcast_last8(v: V, _mask: V) -> V {
        i8x16_shuffle::<15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15>(v, v)
    }
    /// Lane 7 (8-bit) to lanes 8..16, zero elsewhere.
    #[inline(always)]
    pub unsafe fn cross8(v: V, _mask: V) -> V {
        i8x16_shuffle::<16, 16, 16, 16, 16, 16, 16, 16, 7, 7, 7, 7, 7, 7, 7, 7>(v, i16x8_splat(0))
    }
    #[inline(always)]
    pub unsafe fn unpacklo8(a: V, b: V) -> V {
        i8x16_shuffle::<0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23>(a, b)
    }
    #[inline(always)]
    pub unsafe fn unpackhi8(a: V, b: V) -> V {
        i8x16_shuffle::<8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31>(a, b)
    }
    #[inline(always)]
    pub unsafe fn packus16(a: V, b: V) -> V {
        u8x16_narrow_i16x8(a, b)
    }
    #[inline(always)]
    pub unsafe fn mask8(v: V) -> u32 {
        i8x16_bitmask(v) as u32
    }
    #[inline(always)]
    pub unsafe fn lane0(v: V) -> i32 {
        i16x8_extract_lane::<0>(v) as i32
    }
    #[inline(always)]
    pub unsafe fn adds8(a: V, b: V) -> V {
        i8x16_add_sat(a, b)
    }
    #[inline(always)]
    pub unsafe fn add8(a: V, b: V) -> V {
        i8x16_add(a, b)
    }
    #[inline(always)]
    pub unsafe fn subsu8(a: V, b: V) -> V {
        u8x16_sub_sat(a, b)
    }
    #[inline(always)]
    pub unsafe fn max8(a: V, b: V) -> V {
        i8x16_max(a, b)
    }
    #[inline(always)]
    pub unsafe fn gt8(a: V, b: V) -> V {
        i8x16_gt(a, b)
    }
    #[inline(always)]
    pub unsafe fn shlq8(v: V) -> V {
        i64x2_shl(v, 8)
    }
    #[inline(always)]
    pub unsafe fn bsl1(v: V) -> V {
        i8x16_shuffle::<16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14>(v, i16x8_splat(0))
    }
    #[inline(always)]
    pub unsafe fn lane0_u8(v: V) -> i32 {
        u8x16_extract_lane::<0>(v) as i32
    }
    #[inline(always)]
    pub unsafe fn lane15_u8(v: V) -> i32 {
        u8x16_extract_lane::<15>(v) as i32
    }
    #[inline(always)]
    pub unsafe fn lane7(v: V) -> i32 {
        i16x8_extract_lane::<7>(v) as i32
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-619,631
// ```c
//             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
// ...
//             if (score < score_gap_col) {
//                 script = SCRIPT_GAP_IN_B;
//                 score = score_gap_col;
//             }
//             if (score < score_gap_row) {
//                 script = SCRIPT_GAP_IN_A;
//                 score = score_gap_row;
//             }
//             if (best_score - score > x_dropoff) {
// ...
//                 score_gap_row -= gap_extend;
//                 score_gap_col -= gap_extend;
//                 if (score_gap_col < (score - gap_open_extend)) {
//                     score_array[b_index].best_gap = score - gap_open_extend;
// ...
//                 score_array[b_index].best = score;
// ```
// The vector kernels: ALIGN_EX / Blast_SemiGappedAlign (and `s_BlastAlignPackedNucl`) as
// described at the head of this file, with 8-bit lanes (`align_body8`) or 16-bit lanes
// (`align_body`) and with (`TB`) or without a traceback.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
mod kernel {
    use super::arch::*;
    use super::*;

    // No NCBI counterpart: shuffle and mask constants of the 16-bit lanes; they do not change any
    // value NCBI computes.
    static KCROSS: [u8; 16] = [
        128, 128, 128, 128, 128, 128, 128, 128, 6, 7, 6, 7, 6, 7, 6, 7,
    ];
    static KB7: [u8; 16] = [
        14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15, 14, 15,
    ];
    static LANE_MASK: [[i16; 8]; 9] = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1],
    ];

    // No NCBI counterpart: the gap costs and constants of one call, splatted over the lanes
    // (`gap_open_extend`, `gap_extend`, `x_dropoff`, the script bits).
    struct Consts {
        goe: V,
        ge1: V,
        ge2: V,
        ramp8: V,
        ramp_hi: V,
        kcross: V,
        kb7: V,
        xv: V,
        c3: V,
        c10: V,
        c40: V,
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:601,610-614
    // ```c
    //             if (best_score - score > x_dropoff) {
    // ...
    //                 if (score > best_score) {
    //                     best_score = score;
    //                     *a_offset = a_index;
    //                     *b_offset = b_index;
    // ```
    // The running `best_score` of the row, the threshold `best_score - x_dropoff` and the cell where
    // the best score was last raised, kept in vectors across the blocks of one row.
    struct RowState {
        tprev: V,
        bestv: V,
        thrv: V,
        rec_b: usize,
        slow: u32,
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-619,631
    /// ```c
    ///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
    /// ...
    ///             if (score < score_gap_col) {
    ///                 script = SCRIPT_GAP_IN_B;
    ///                 score = score_gap_col;
    ///             }
    ///             if (score < score_gap_row) {
    ///                 script = SCRIPT_GAP_IN_A;
    ///                 score = score_gap_row;
    ///             }
    ///             if (best_score - score > x_dropoff) {
    /// ...
    ///                 score_gap_row -= gap_extend;
    ///                 score_gap_col -= gap_extend;
    ///                 if (score_gap_col < (score - gap_open_extend)) {
    ///                     score_array[b_index].best_gap = score - gap_open_extend;
    /// ...
    ///                 score_array[b_index].best = score;
    /// ```
    /// One block of 8 cells of one row. `D` is `next_score` of the diagonal, `E` is `score_gap_col`
    /// (`best_gap`), `F` is `score_gap_row`, `H` is `score`. The comparisons set the same script bits
    /// as the `if`s of the C loop. The cells that fail the X-drop test are returned as a mask and
    /// stored as 0 (minus infinity). The argument at the head of this file says why the kept cells
    /// have the C values.
    ///
    /// One block of 8 lanes at `b`. Returns (drop mask, script).
    #[inline(always)]
    unsafe fn block<const TB: bool>(
        c: &Consts,
        rs: &mut RowState,
        hr: *const i16,
        hw: *mut i16,
        e: *mut i16,
        b: usize,
        s: V,
        lane_mask: Option<V>,
    ) -> (V, V) {
        let hp = load(hr.add(b).sub(1));
        let ep = load(e.add(b));
        let d = adds16(hp, s);
        let h0 = max16(d, ep);
        let v = subsu16(h0, c.goe);
        // decayed prefix inside each 4-lane half
        let i1 = max16(v, subsu16(shlq16(v), c.ge1));
        let i2 = max16(i1, subsu16(shlq32(i1), c.ge2));
        let fx = shlq16(i2);
        let cross = subsu16(cross16(i2, c.kcross), c.ramp_hi);
        let fc = subsu16(bcast_last16(rs.tprev, c.kb7), c.ramp8);
        let f = max16(max16(fx, cross), fc);
        let mut hh = max16(h0, f);
        if let Some(m) = lane_mask {
            hh = and(hh, m);
        }
        let hg = subsu16(hh, c.goe);
        let fs = subsu16(f, c.ge1);
        rs.tprev = max16(hg, fs);
        let es = subsu16(ep, c.ge1);
        let en = max16(hg, es);
        let gt = gt16(hh, rs.bestv);
        let thr;
        if mask8(gt) != 0 {
            // A new best inside this block: exact running best per lane.
            let y0 = bsl2(hh);
            let y1 = max16(y0, bsl2(y0));
            let y2 = max16(y1, bsl4(y1));
            let y3 = max16(y2, bsl8(y2));
            let pfx = max16(y3, rs.bestv);
            thr = subsu16(pfx, c.xv);
            let rec = gt16(hh, pfx);
            let rm = mask8(rec);
            rs.rec_b = b + ((31 - rm.leading_zeros()) >> 1) as usize;
            rs.bestv = bcast_last16(max16(pfx, hh), c.kb7);
            rs.thrv = subsu16(rs.bestv, c.xv);
            rs.slow += 1;
        } else {
            thr = rs.thrv;
        }
        let drop = gt16(thr, hh);
        store(hw.add(b), andnot(drop, hh));
        store(e.add(b), en);
        let scr = if TB {
            let me = gt16(ep, d);
            let mf = gt16(f, h0);
            let op = andnot(mf, add16(c.c3, and(me, c.c3)));
            let fb = andnot(gt16(hg, es), c.c40);
            let fa = andnot(gt16(hg, fs), c.c10);
            or(op, or(fb, fa))
        } else {
            zero()
        };
        (drop, scr)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:476-500,641-675
    // ```c
    //     score_array[0].best = 0;
    //     score_array[0].best_gap = -gap_open_extend;
    // ...
    //     for (a_index = 1; a_index <= M; a_index++) {
    // ...
    //         if (last_b_index < b_size - 1) {
    //             b_size = last_b_index + 1;
    //         }
    //         else {
    //             while (score_gap_row >= (best_score - x_dropoff) && b_size <= N) {
    // ...
    //         if (b_size <= N) {
    //             score_array[b_size].best = MININT;
    //             score_array[b_size].best_gap = MININT;
    //             b_size++;
    // ```
    // The whole DP on 16-bit lanes: row 0, then one row per `a_index`, with the band
    // `first_b_index .. b_size` narrowed and widened as in the C loop. It returns None for a
    // problem that does not fit (the caller then runs the scalar port).
    #[inline(always)]
    pub(super) unsafe fn align_body<const TB: bool>(
        q: &RowSeq<'_>,
        s: &ColSeq<'_>,
        scores: &Scores<'_>,
        len1: usize,
        len2: usize,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        check_fence: bool,
        sc: &mut XdropScratch,
    ) -> Option<XdropResult> {
        sc.ops.clear();
        if len1 == 0 || len2 == 0 {
            return Some(XdropResult {
                a_offset: 0,
                b_offset: 0,
                score: 0,
                fence_hit: false,
                cells: 0,
            });
        }
        let goe = gap_open.checked_add(gap_extend)?;
        let mut x = x_drop;
        if x < goe {
            x = goe;
        }
        if gap_extend < 1 || gap_open < 0 {
            return None;
        }
        // Stored value = score - base + bias. Every value that can still
        // matter (>= best - X) stays exact after subtracting a gap cost, and
        // anything derived from -inf stays below every threshold.
        let bias = x as i64 + goe as i64 + MAX_M as i64 + 1;
        let limit_rebase = VMAX - bias - MAX_M as i64;
        if limit_rebase < MAX_M as i64 || gap_extend as i64 * 8 > VMAX {
            return None;
        }
        let bias = bias as i32;
        let limit_rebase = limit_rebase as i32;
        let num_extra = (x / gap_extend + 3) as usize;
        sc.reset();

        // Row 0 goes to buffer A.
        sc.ensure(num_extra.min(len2 + 2) + 16);
        let mut b_size = 1usize;
        {
            let h = sc.ha.as_mut_ptr().add(PAD);
            let e = sc.e.as_mut_ptr().add(PAD);
            *h = bias as i16;
            *e = (bias - goe) as i16;
            let mut score = -goe;
            for i in 1..=len2 {
                if score < -x {
                    break;
                }
                *h.add(i) = (bias + score) as i16;
                *e.add(i) = (bias + score - goe) as i16;
                score -= gap_extend;
                b_size = i + 1;
            }
            // In row 1 the scalar loop never forms the diagonal into column
            // b_size (that column is only created by the row-gap extension).
            *h.add(b_size - 1) = 0;
        }
        sc.dirty = b_size;
        if TB {
            sc.trace.resize(b_size, SCRIPT_GAP_IN_A);
            sc.rows.push((0, 0, b_size as u32));
        }

        let g = gap_extend as i16;
        let ramp8: [i16; 8] = [0, g, 2 * g, 3 * g, 4 * g, 5 * g, 6 * g, 7 * g];
        let ramp_hi: [i16; 8] = [0, 0, 0, 0, 0, g, 2 * g, 3 * g];
        let c = Consts {
            goe: splat16(goe as i16),
            ge1: splat16(g),
            ge2: splat16(2 * g),
            ramp8: load(ramp8.as_ptr()),
            ramp_hi: load(ramp_hi.as_ptr()),
            kcross: load(core::hint::black_box(&KCROSS).as_ptr()),
            kb7: load(core::hint::black_box(&KB7).as_ptr()),
            xv: splat16(x as i16),
            c3: splat16(SCRIPT_SUB as i16),
            c10: splat16(SCRIPT_EXTEND_GAP_A as i16),
            c40: splat16(SCRIPT_EXTEND_GAP_B as i16),
        };
        let zero_v = zero();
        let static_tabs: Option<*const u8> = match scores {
            Scores::Tables { t, .. } => Some(t.data.as_ptr()),
            _ => None,
        };
        let ncols = scores.ncols();
        let mut cells: u64 = 0;

        let mut base: i32 = 0;
        let mut best_rel: i32 = 0;
        let mut a_off = 0usize;
        let mut b_off = 0usize;
        let mut first = 0usize;
        let mut fence_hit = false;
        // exclusive end of the lanes last written in each H buffer
        let mut wend = [b_size, 0usize];
        let mut cur = 0usize; // buffer holding the previous row

        for a in 1..=len1 {
            sc.ensure(b_size + num_extra + 64);
            if sc.sb_filled < b_size {
                let to = (b_size + 64).min(len2 + 1);
                for k in sc.sb_filled..to {
                    let ch = s.get(k);
                    if ch as usize >= ncols {
                        if check_fence && ch == FENCE_SENTRY {
                            if sc.fence_k == usize::MAX {
                                sc.fence_k = k;
                            }
                        } else if sc.bad_k == usize::MAX {
                            sc.bad_k = k;
                        }
                    }
                    *sc.sb.get_unchecked_mut(PAD + k + 1) = ch;
                }
                sc.sb_filled = to;
            }
            if sc.bad_k < b_size {
                return None;
            }
            let mut limit = len2 + 1;
            let mut bend = b_size;
            let mut fence_row = false;
            if sc.fence_k < b_size {
                limit = sc.fence_k;
                bend = sc.fence_k;
                fence_row = true;
                if bend <= first {
                    fence_hit = true;
                    break;
                }
            }
            let qa = q.get(a)? as usize;
            if qa >= ncols {
                return None;
            }
            cells += (b_size - first) as u64;
            let t: *const u8 = if let Some(st) = static_tabs {
                st.add(qa * TAB)
            } else {
                if sc.tab_ok & (1u32 << qa) == 0 {
                    if !build_row_table(scores, qa, &mut sc.tabs[qa * TAB..(qa + 1) * TAB]) {
                        return None;
                    }
                    sc.tab_ok |= 1u32 << qa;
                }
                sc.tabs.as_ptr().add(qa * TAB)
            };
            let (hr, hw) = if cur == 0 {
                (sc.ha.as_mut_ptr().add(PAD), sc.hb.as_mut_ptr().add(PAD))
            } else {
                (sc.hb.as_mut_ptr().add(PAD), sc.ha.as_mut_ptr().add(PAD))
            };
            let e = sc.e.as_mut_ptr().add(PAD);
            let sb = sc.sb.as_ptr().add(PAD);
            let t0lo = load(t);
            let t1lo = load(t.add(16));
            let t0hi = load(t.add(32));
            let t1hi = load(t.add(48));

            if best_rel > limit_rebase {
                let dv = splat16(best_rel as i16);
                let mut b = first;
                while b < b_size {
                    store(hr.add(b), subsu16(load(hr.add(b)), dv));
                    store(e.add(b), subsu16(load(e.add(b)), dv));
                    b += 8;
                }
                base += best_rel;
                best_rel = 0;
            }
            // no diagonal into the first lane of the band
            *hr.add(first).sub(1) = 0;

            let row_off = sc.trace.len();
            if TB {
                sc.trace.reserve(bend - first + num_extra + 96);
            }
            let tr = sc.trace.as_mut_ptr().add(row_off);
            let survw = sc.surv.as_mut_ptr();

            let bestv0 = splat16((bias + best_rel) as i16);
            let mut rs = RowState {
                tprev: zero_v,
                bestv: bestv0,
                thrv: subsu16(bestv0, c.xv),
                rec_b: 0,
                slow: 0,
            };
            let mut b0 = first;
            let mut it = 0usize;
            loop {
                // scores for the 16 lanes b0..b0+16 (lane b reads column b-1)
                let idx = load(sb.add(b0));
                let lo = lut32(t0lo, t1lo, idx);
                let hi = lut32(t0hi, t1hi, idx);
                let s_lo = unpacklo8(lo, hi);
                let s_hi = unpackhi8(lo, hi);
                let (d0, sc0, d1, sc1);
                if b0 + 16 <= limit {
                    (d0, sc0) = block::<TB>(&c, &mut rs, hr, hw, e, b0, s_lo, None);
                    (d1, sc1) = block::<TB>(&c, &mut rs, hr, hw, e, b0 + 8, s_hi, None);
                } else {
                    let n = limit - b0;
                    let m0 = load(LANE_MASK[n.min(8)].as_ptr());
                    let m1 = load(LANE_MASK[n.saturating_sub(8)].as_ptr());
                    (d0, sc0) = block::<TB>(&c, &mut rs, hr, hw, e, b0, s_lo, Some(m0));
                    (d1, sc1) = block::<TB>(&c, &mut rs, hr, hw, e, b0 + 8, s_hi, Some(m1));
                }
                if TB {
                    store(tr.add(b0 - first), packus16(sc0, sc1));
                }
                // 2 bits per lane, 32 bits per iteration
                let dm = mask8(d0) | (mask8(d1) << 16);
                let w = it >> 1;
                if it & 1 == 0 {
                    *survw.add(w) = (!dm) as u64;
                } else {
                    *survw.add(w) |= ((!dm) as u64) << 32;
                }
                it += 1;
                b0 += 16;
                if b0 >= bend {
                    if b0 >= limit {
                        break;
                    }
                    if lane7(rs.tprev) < lane0(rs.thrv) {
                        break;
                    }
                }
            }
            let proc_end = b0;
            if proc_end > sc.dirty {
                sc.dirty = proc_end;
            }
            // lanes this buffer held beyond what the row wrote are stale
            let wi = 1 - cur;
            if wend[wi] > proc_end {
                let mut b = proc_end;
                while b < wend[wi] {
                    store(hw.add(b), zero_v);
                    b += 8;
                }
            }
            wend[wi] = proc_end;
            cur = wi;
            if rs.slow != 0 {
                b_off = rs.rec_b;
                a_off = a;
                best_rel = lane0(rs.bestv) - bias;
            }
            if TB {
                sc.trace.set_len(row_off + proc_end - first);
                sc.rows
                    .push((row_off as u32, first as u32, (proc_end - first) as u32));
            }
            if fence_row {
                fence_hit = true;
                break;
            }
            // survivors
            let nwords = (it + 1) >> 1;
            let mut wl = 0usize;
            while wl < nwords && *survw.add(wl) == 0 {
                wl += 1;
            }
            if wl == nwords {
                break;
            }
            let first_new = first + wl * 32 + ((*survw.add(wl)).trailing_zeros() as usize >> 1);
            let mut wr = nwords - 1;
            while *survw.add(wr) == 0 {
                wr -= 1;
            }
            let last = first + wr * 32 + ((63 - (*survw.add(wr)).leading_zeros() as usize) >> 1);
            // columns right of the last survivor leave the band: E = -inf
            store(e.add(last + 1), zero_v);
            store(e.add(last + 9), zero_v);
            if proc_end > last + 17 {
                let mut b = last + 17;
                while b < proc_end {
                    *e.add(b) = 0;
                    b += 1;
                }
            }
            first = first_new;
            b_size = last + 1;
            if b_size <= len2 {
                b_size += 1;
            }
        }

        let best = base + best_rel;
        finish::<TB>(sc, a_off, b_off, best, fence_hit, cells)
    }

    // No NCBI counterpart: shuffle and mask constants of the 8-bit lanes; they do not change any
    // value NCBI computes.
    // -----------------------------------------------------------------
    // 8-bit lanes (16 per vector): same algorithm, used when
    // x_dropoff + gap costs + 2 * max_score fit in 7 bits.
    // -----------------------------------------------------------------

    static KCROSS8: [u8; 16] = [
        128, 128, 128, 128, 128, 128, 128, 128, 7, 7, 7, 7, 7, 7, 7, 7,
    ];
    static KB15: [u8; 16] = [15; 16];
    /// `LANE_MASK8[16 - n ..]` has the first `n` bytes set.
    static LANE_MASK8: [u8; 32] = [
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    // No NCBI counterpart: the gap costs and constants of one call on 8-bit lanes.
    struct Consts8 {
        goe: V,
        ge1: V,
        ge2: V,
        ge4: V,
        ramp16: V,
        ramp_hi: V,
        kcross: V,
        kb15: V,
        xv: V,
        c3: V,
        c10: V,
        c40: V,
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-601,616-619,631
    /// ```c
    ///             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
    /// ...
    ///             if (score < score_gap_col) {
    ///                 script = SCRIPT_GAP_IN_B;
    ///                 score = score_gap_col;
    ///             }
    ///             if (score < score_gap_row) {
    ///                 script = SCRIPT_GAP_IN_A;
    ///                 score = score_gap_row;
    ///             }
    ///             if (best_score - score > x_dropoff) {
    /// ...
    ///                 score_gap_row -= gap_extend;
    ///                 score_gap_col -= gap_extend;
    ///                 if (score_gap_col < (score - gap_open_extend)) {
    ///                     score_array[b_index].best_gap = score - gap_open_extend;
    /// ...
    ///                 score_array[b_index].best = score;
    /// ```
    /// The same block as `block`, on 16 cells of 8 bits. Scores are kept relative to a base so that
    /// every value that can still matter fits in 7 bits; `fits_8bit` checks this before the
    /// call.
    ///
    /// One block of 16 lanes at `b`. Returns (drop mask, script).
    #[inline(always)]
    unsafe fn block8<const TB: bool>(
        c: &Consts8,
        rs: &mut RowState,
        hr: *const u8,
        hw: *mut u8,
        e: *mut u8,
        b: usize,
        s: V,
        lane_mask: Option<V>,
    ) -> (V, V) {
        let hp = load(hr.add(b).sub(1));
        let ep = load(e.add(b));
        let d = adds8(hp, s);
        let h0 = max8(d, ep);
        let v = subsu8(h0, c.goe);
        // decayed prefix inside each 8-lane half
        let i1 = max8(v, subsu8(shlq8(v), c.ge1));
        let i2 = max8(i1, subsu8(shlq16(i1), c.ge2));
        let i3 = max8(i2, subsu8(shlq32(i2), c.ge4));
        let fx = shlq8(i3);
        let cross = subsu8(cross8(i3, c.kcross), c.ramp_hi);
        let fc = subsu8(bcast_last8(rs.tprev, c.kb15), c.ramp16);
        let f = max8(max8(fx, cross), fc);
        let mut hh = max8(h0, f);
        if let Some(m) = lane_mask {
            hh = and(hh, m);
        }
        let hg = subsu8(hh, c.goe);
        let fs = subsu8(f, c.ge1);
        rs.tprev = max8(hg, fs);
        let es = subsu8(ep, c.ge1);
        let en = max8(hg, es);
        let gt = gt8(hh, rs.bestv);
        let thr;
        if mask8(gt) != 0 {
            // A new best inside this block: exact running best per lane.
            let y0 = bsl1(hh);
            let y1 = max8(y0, bsl1(y0));
            let y2 = max8(y1, bsl2(y1));
            let y3 = max8(y2, bsl4(y2));
            let y4 = max8(y3, bsl8(y3));
            let pfx = max8(y4, rs.bestv);
            thr = subsu8(pfx, c.xv);
            let rec = gt8(hh, pfx);
            let rm = mask8(rec);
            rs.rec_b = b + (31 - rm.leading_zeros()) as usize;
            rs.bestv = bcast_last8(max8(pfx, hh), c.kb15);
            rs.thrv = subsu8(rs.bestv, c.xv);
            rs.slow += 1;
        } else {
            thr = rs.thrv;
        }
        let drop = gt8(thr, hh);
        store(hw.add(b), andnot(drop, hh));
        store(e.add(b), en);
        let scr = if TB {
            let me = gt8(ep, d);
            let mf = gt8(f, h0);
            let op = andnot(mf, add8(c.c3, and(me, c.c3)));
            let fb = andnot(gt8(hg, es), c.c40);
            let fa = andnot(gt8(hg, fs), c.c10);
            or(op, or(fb, fa))
        } else {
            zero()
        };
        (drop, scr)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:476-500,641-675
    // ```c
    //     score_array[0].best = 0;
    //     score_array[0].best_gap = -gap_open_extend;
    // ...
    //     for (a_index = 1; a_index <= M; a_index++) {
    // ...
    //         if (last_b_index < b_size - 1) {
    //             b_size = last_b_index + 1;
    // ```
    // The same DP as `align_body` on 8-bit lanes (16 cells per vector). It returns None for a
    // problem that does not fit.
    /// `max_m` is the largest score of the matrix (>= 0).
    #[inline(always)]
    pub(super) unsafe fn align_body8<const TB: bool>(
        q: &RowSeq<'_>,
        s: &ColSeq<'_>,
        scores: &Scores<'_>,
        max_m: i32,
        len1: usize,
        len2: usize,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        check_fence: bool,
        sc: &mut XdropScratch,
    ) -> Option<XdropResult> {
        sc.ops.clear();
        if len1 == 0 || len2 == 0 {
            return Some(XdropResult {
                a_offset: 0,
                b_offset: 0,
                score: 0,
                fence_hit: false,
                cells: 0,
            });
        }
        let goe = gap_open.checked_add(gap_extend)?;
        let mut x = x_drop;
        if x < goe {
            x = goe;
        }
        if gap_extend < 1 || gap_open < 0 || max_m < 0 {
            return None;
        }
        let bias = x as i64 + goe as i64 + max_m as i64 + 1;
        let limit_rebase = 127 - bias - max_m as i64;
        if limit_rebase < 0 {
            return None;
        }
        let bias = bias as i32;
        let limit_rebase = limit_rebase as i32;
        let num_extra = (x / gap_extend + 3) as usize;
        let ncols = scores.ncols();
        sc.reset();

        // Row 0 goes to buffer A.
        sc.ensure(num_extra.min(len2 + 2) + 16);
        let mut b_size = 1usize;
        {
            let h = sc.h8a.as_mut_ptr().add(PAD);
            let e = sc.e8.as_mut_ptr().add(PAD);
            *h = bias as u8;
            *e = (bias - goe) as u8;
            let mut score = -goe;
            for i in 1..=len2 {
                if score < -x {
                    break;
                }
                *h.add(i) = (bias + score) as u8;
                *e.add(i) = (bias + score - goe) as u8;
                score -= gap_extend;
                b_size = i + 1;
            }
            *h.add(b_size - 1) = 0;
        }
        sc.dirty8 = b_size;
        if TB {
            sc.trace.resize(b_size, SCRIPT_GAP_IN_A);
            sc.rows.push((0, 0, b_size as u32));
        }

        let g = gap_extend;
        let mut ramp16 = [0u8; 16];
        let mut ramp_hi = [0u8; 16];
        for k in 0..16 {
            ramp16[k] = (k as i32 * g).min(255) as u8;
            if k >= 8 {
                ramp_hi[k] = ((k as i32 - 8) * g).min(255) as u8;
            }
        }
        let sat = |v: i32| v.min(255) as u8 as i8;
        let c = Consts8 {
            goe: splat8(sat(goe)),
            ge1: splat8(sat(g)),
            ge2: splat8(sat(2 * g)),
            ge4: splat8(sat(4 * g)),
            ramp16: load(ramp16.as_ptr()),
            ramp_hi: load(ramp_hi.as_ptr()),
            kcross: load(core::hint::black_box(&KCROSS8).as_ptr()),
            kb15: load(core::hint::black_box(&KB15).as_ptr()),
            xv: splat8(sat(x)),
            c3: splat8(SCRIPT_SUB as i8),
            c10: splat8(SCRIPT_EXTEND_GAP_A as i8),
            c40: splat8(SCRIPT_EXTEND_GAP_B as i8),
        };
        let zero_v = zero();
        let static_tabs: Option<*const u8> = match scores {
            Scores::Tables { t, .. } => Some(t.data.as_ptr()),
            _ => None,
        };
        let small_alphabet = ncols <= 16;
        let mut cells: u64 = 0;

        let mut base: i32 = 0;
        let mut best_rel: i32 = 0;
        let mut a_off = 0usize;
        let mut b_off = 0usize;
        let mut first = 0usize;
        let mut fence_hit = false;
        let mut wend = [b_size, 0usize];
        let mut cur = 0usize;

        for a in 1..=len1 {
            sc.ensure(b_size + num_extra + 64);
            if sc.sb_filled < b_size {
                let to = (b_size + 64).min(len2 + 1);
                for k in sc.sb_filled..to {
                    let ch = s.get(k);
                    if ch as usize >= ncols {
                        if check_fence && ch == FENCE_SENTRY {
                            if sc.fence_k == usize::MAX {
                                sc.fence_k = k;
                            }
                        } else if sc.bad_k == usize::MAX {
                            sc.bad_k = k;
                        }
                    }
                    *sc.sb.get_unchecked_mut(PAD + k + 1) = ch;
                }
                sc.sb_filled = to;
            }
            if sc.bad_k < b_size {
                return None;
            }
            let mut limit = len2 + 1;
            let mut bend = b_size;
            let mut fence_row = false;
            if sc.fence_k < b_size {
                limit = sc.fence_k;
                bend = sc.fence_k;
                fence_row = true;
                if bend <= first {
                    fence_hit = true;
                    break;
                }
            }
            let qa = q.get(a)? as usize;
            if qa >= ncols {
                return None;
            }
            cells += (b_size - first) as u64;
            let t: *const u8 = if let Some(st) = static_tabs {
                st.add(qa * TAB)
            } else {
                if sc.tab_ok & (1u32 << qa) == 0 {
                    if !build_row_table(scores, qa, &mut sc.tabs[qa * TAB..(qa + 1) * TAB]) {
                        return None;
                    }
                    sc.tab_ok |= 1u32 << qa;
                }
                sc.tabs.as_ptr().add(qa * TAB)
            };
            let (hr, hw) = if cur == 0 {
                (sc.h8a.as_mut_ptr().add(PAD), sc.h8b.as_mut_ptr().add(PAD))
            } else {
                (sc.h8b.as_mut_ptr().add(PAD), sc.h8a.as_mut_ptr().add(PAD))
            };
            let e = sc.e8.as_mut_ptr().add(PAD);
            let sb = sc.sb.as_ptr().add(PAD);
            let t0 = load(t.add(64));
            let t1 = load(t.add(80));

            if best_rel > limit_rebase {
                let dv = splat8(best_rel as u8 as i8);
                let mut b = first;
                while b < b_size {
                    store(hr.add(b), subsu8(load(hr.add(b)), dv));
                    store(e.add(b), subsu8(load(e.add(b)), dv));
                    b += 16;
                }
                base += best_rel;
                best_rel = 0;
            }
            *hr.add(first).sub(1) = 0;

            let row_off = sc.trace.len();
            if TB {
                sc.trace.reserve(bend - first + num_extra + 96);
            }
            let tr = sc.trace.as_mut_ptr().add(row_off);
            let survw = sc.surv.as_mut_ptr();

            let bestv0 = splat8((bias + best_rel) as u8 as i8);
            let mut rs = RowState {
                tprev: zero_v,
                bestv: bestv0,
                thrv: subsu8(bestv0, c.xv),
                rec_b: 0,
                slow: 0,
            };
            let mut b0 = first;
            let mut it = 0usize;
            loop {
                // scores for the 16 lanes b0..b0+16 (lane b reads column b-1)
                let idx = load(sb.add(b0));
                let sv = if small_alphabet {
                    pick(t0, idx)
                } else {
                    lut32(t0, t1, idx)
                };
                let (d0, sc0);
                if b0 + 16 <= limit {
                    (d0, sc0) = block8::<TB>(&c, &mut rs, hr, hw, e, b0, sv, None);
                } else {
                    let n = limit - b0;
                    let m0 = load(LANE_MASK8.as_ptr().add(16 - n));
                    (d0, sc0) = block8::<TB>(&c, &mut rs, hr, hw, e, b0, sv, Some(m0));
                }
                if TB {
                    store(tr.add(b0 - first), sc0);
                }
                // 1 bit per lane, 16 bits per iteration
                let live = (!mask8(d0)) & 0xFFFF;
                let w = it >> 2;
                let sh = (it & 3) * 16;
                if sh == 0 {
                    *survw.add(w) = live as u64;
                } else {
                    *survw.add(w) |= (live as u64) << sh;
                }
                it += 1;
                b0 += 16;
                if b0 >= bend {
                    if b0 >= limit {
                        break;
                    }
                    if lane15_u8(rs.tprev) < lane0_u8(rs.thrv) {
                        break;
                    }
                }
            }
            let proc_end = b0;
            if proc_end > sc.dirty8 {
                sc.dirty8 = proc_end;
            }
            let wi = 1 - cur;
            if wend[wi] > proc_end {
                let mut b = proc_end;
                while b < wend[wi] {
                    store(hw.add(b), zero_v);
                    b += 16;
                }
            }
            wend[wi] = proc_end;
            cur = wi;
            if rs.slow != 0 {
                b_off = rs.rec_b;
                a_off = a;
                best_rel = lane0_u8(rs.bestv) - bias;
            }
            if TB {
                sc.trace.set_len(row_off + proc_end - first);
                sc.rows
                    .push((row_off as u32, first as u32, (proc_end - first) as u32));
            }
            if fence_row {
                fence_hit = true;
                break;
            }
            let nwords = (it + 3) >> 2;
            let mut wl = 0usize;
            while wl < nwords && *survw.add(wl) == 0 {
                wl += 1;
            }
            if wl == nwords {
                break;
            }
            let first_new = first + wl * 64 + (*survw.add(wl)).trailing_zeros() as usize;
            let mut wr = nwords - 1;
            while *survw.add(wr) == 0 {
                wr -= 1;
            }
            let last = first + wr * 64 + (63 - (*survw.add(wr)).leading_zeros() as usize);
            store(e.add(last + 1), zero_v);
            if proc_end > last + 17 {
                let mut b = last + 17;
                while b < proc_end {
                    *e.add(b) = 0;
                    b += 1;
                }
            }
            first = first_new;
            b_size = last + 1;
            if b_size <= len2 {
                b_size += 1;
            }
        }

        let best = base + best_rel;
        finish::<TB>(sc, a_off, b_off, best, fence_hit, cells)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:682-727
    /// ```c
    ///     a_index = *a_offset;
    ///     b_index = *b_offset;
    ///     script = SCRIPT_SUB;
    /// ...
    ///     while (a_index > 0 || b_index > 0) {
    /// ...
    ///         GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
    /// ```
    /// The end of `ALIGN_EX`: the walk from the best cell back to the start, with the same state
    /// machine on `script`. Score-only calls return the offsets and the score, as
    /// `Blast_SemiGappedAlign` does.
    ///
    /// Result assembly and traceback walk shared by both lane widths.
    #[inline(always)]
    unsafe fn finish<const TB: bool>(
        sc: &mut XdropScratch,
        a_off: usize,
        b_off: usize,
        best: i32,
        fence_hit: bool,
        cells: u64,
    ) -> Option<XdropResult> {
        if !TB {
            return Some(XdropResult {
                a_offset: a_off,
                b_offset: b_off,
                score: best,
                fence_hit: false,
                cells,
            });
        }
        if best <= 0 {
            return Some(XdropResult {
                a_offset: 0,
                b_offset: 0,
                score: 0,
                fence_hit,
                cells,
            });
        }
        if fence_hit {
            return Some(XdropResult {
                a_offset: a_off,
                b_offset: b_off,
                score: best,
                fence_hit: true,
                cells,
            });
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:689-726
        // ```c
        // while (a_index > 0 || b_index > 0) {
        //     next_script = edit_script[a_index][b_index - edit_start_offset[a_index]];
        //     switch(script) { ... }
        //     GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
        // }
        // ```
        let mut a_index = a_off;
        let mut b_index = b_off;
        let mut script = SCRIPT_SUB;
        while a_index > 0 || b_index > 0 {
            if a_index >= sc.rows.len() {
                break;
            }
            let (off, start, len) = sc.rows[a_index];
            if b_index < start as usize {
                break;
            }
            let b_rel = b_index - start as usize;
            if b_rel >= len as usize {
                break;
            }
            let next = *sc.trace.get_unchecked(off as usize + b_rel);
            match script & SCRIPT_OP_MASK {
                x if x == SCRIPT_GAP_IN_A => {
                    script = next & SCRIPT_OP_MASK;
                    if next & SCRIPT_EXTEND_GAP_A != 0 {
                        script = SCRIPT_GAP_IN_A;
                    }
                }
                x if x == SCRIPT_GAP_IN_B => {
                    script = next & SCRIPT_OP_MASK;
                    if next & SCRIPT_EXTEND_GAP_B != 0 {
                        script = SCRIPT_GAP_IN_B;
                    }
                }
                _ => {
                    script = next & SCRIPT_OP_MASK;
                }
            }
            let op = script & SCRIPT_OP_MASK;
            if op == SCRIPT_GAP_IN_A {
                b_index = b_index.saturating_sub(1);
            } else if op == SCRIPT_GAP_IN_B {
                a_index = a_index.saturating_sub(1);
            } else {
                a_index = a_index.saturating_sub(1);
                b_index = b_index.saturating_sub(1);
            }
            push_op(&mut sc.ops, op);
        }
        Some(XdropResult {
            a_offset: a_off,
            b_offset: b_off,
            score: best,
            fence_hit: false,
            cells,
        })
    }

    // No NCBI counterpart: entry points compiled with CPU features enabled (AVX or SSE) around
    // `align_body` and `align_body8`. The same code runs; only the instructions the compiler may use
    // change.
    macro_rules! x86_entry {
        ($name:ident, $name8:ident, $feat:literal) => {
            #[cfg(target_arch = "x86_64")]
            #[target_feature(enable = $feat)]
            pub(super) unsafe fn $name<const TB: bool>(
                q: &RowSeq<'_>,
                s: &ColSeq<'_>,
                scores: &Scores<'_>,
                len1: usize,
                len2: usize,
                gap_open: i32,
                gap_extend: i32,
                x_drop: i32,
                check_fence: bool,
                sc: &mut XdropScratch,
            ) -> Option<XdropResult> {
                align_body::<TB>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                )
            }

            #[cfg(target_arch = "x86_64")]
            #[target_feature(enable = $feat)]
            pub(super) unsafe fn $name8<const TB: bool>(
                q: &RowSeq<'_>,
                s: &ColSeq<'_>,
                scores: &Scores<'_>,
                max_m: i32,
                len1: usize,
                len2: usize,
                gap_open: i32,
                gap_extend: i32,
                x_drop: i32,
                check_fence: bool,
                sc: &mut XdropScratch,
            ) -> Option<XdropResult> {
                align_body8::<TB>(
                    q,
                    s,
                    scores,
                    max_m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                )
            }
        };
    }
    x86_entry!(align_avx, align8_avx, "avx,ssse3,sse4.1,bmi1,lzcnt,popcnt");
    x86_entry!(align_sse, align8_sse, "ssse3,sse4.1");
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
/// ```c
///     if (!score_only) {
///         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
/// ```
/// Reads the LOSAT_X_DPFAST / LOSAT_X_DPSHADOW switches once. The mode only chooses whether the
/// scalar port of these functions (0), the vector kernels (1), or both with a comparison (2) run.
///
/// 0 = scalar only, 1 = SIMD, 2 = SIMD and compare every call with the scalar kernel.
pub(crate) fn mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_DPSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_DPFAST").is_some() {
            1
        } else {
            0
        }
    })
}

/// No NCBI counterpart: checks the CPU features at run time (LOSAT_X_DPNOAVX turns the AVX
/// level off for timing); it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
fn cpu_level() -> u8 {
    use std::sync::OnceLock;
    static LEVEL: OnceLock<u8> = OnceLock::new();
    *LEVEL.get_or_init(|| {
        let sse = is_x86_feature_detected!("ssse3") && is_x86_feature_detected!("sse4.1");
        if sse
            && std::env::var_os("LOSAT_X_DPNOAVX").is_none()
            && is_x86_feature_detected!("avx")
            && is_x86_feature_detected!("bmi1")
            && is_x86_feature_detected!("lzcnt")
            && is_x86_feature_detected!("popcnt")
        {
            2
        } else if sse {
            1
        } else {
            0
        }
    })
}

/// No NCBI counterpart: switch reader for timing; it does not change any value NCBI computes.
///
/// LOSAT_X_DPNO8 keeps every call on the 16-bit lanes (for A/B timing).
fn use_8bit() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_DPNO8").is_none())
}

/// No NCBI counterpart: the test of whether the problem fits 8-bit lanes; a problem that does not
/// fit runs on 16-bit lanes or on the scalar port.
///
/// Largest matrix score if the 8-bit lanes can hold this problem.
fn fits_8bit(scores: &Scores<'_>, gap_open: i32, gap_extend: i32, x_drop: i32) -> Option<i32> {
    if !use_8bit() {
        return None;
    }
    let goe = gap_open.checked_add(gap_extend)?;
    let x = x_drop.max(goe);
    // cheap rejection before looking at the matrix
    if x as i64 + goe as i64 + 1 > 127 {
        return None;
    }
    let max_m = max_score(scores)?.max(0);
    let bias = x as i64 + goe as i64 + max_m as i64 + 1;
    if 127 - bias - max_m as i64 >= 0 {
        Some(max_m)
    } else {
        None
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:736-741,766-770
/// ```c
/// Blast_SemiGappedAlign(const Uint1* A, const Uint1* B, Int4 M, Int4 N,
///    Int4* a_offset, Int4* b_offset, Boolean score_only,
///    GapPrelimEditBlock *edit_block, BlastGapAlignStruct* gap_align,
///    const BlastScoringParameters* score_params,
///    Int4 query_offset, Boolean reversed, Boolean reverse_sequence,
///    Boolean * fence_hit)
/// ...
///     if (!score_only) {
///         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
/// ```
/// Dispatch point of LOSAT_X_DPFAST / LOSAT_X_DPSHADOW for all callers. `tb == false` stands for
/// `score_only` (`Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl`), `tb == true` for `ALIGN_EX`.
/// It picks the AVX or SSE build and the 8-bit or 16-bit lanes. `None` means that the caller must
/// run the scalar port of the C function; callers do so, and compare when the mode is 2.
///
/// Score-only (`tb == false`) or traceback X-drop alignment.
/// `None` means "not handled here": the caller must run the scalar kernel.
pub(crate) fn xdrop_align(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    tb: bool,
    check_fence: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    #[cfg(target_arch = "x86_64")]
    {
        let level = cpu_level();
        if level == 0 {
            return None;
        }
        let m8 = fits_8bit(scores, gap_open, gap_extend, x_drop);
        // SAFETY: the required CPU features were detected at run time.
        unsafe {
            match (level, tb, m8) {
                (2, true, Some(m)) => kernel::align8_avx::<true>(
                    q,
                    s,
                    scores,
                    m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (2, false, Some(m)) => kernel::align8_avx::<false>(
                    q,
                    s,
                    scores,
                    m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (_, true, Some(m)) => kernel::align8_sse::<true>(
                    q,
                    s,
                    scores,
                    m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (_, false, Some(m)) => kernel::align8_sse::<false>(
                    q,
                    s,
                    scores,
                    m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (2, true, None) => kernel::align_avx::<true>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (2, false, None) => kernel::align_avx::<false>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (_, true, None) => kernel::align_sse::<true>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (_, false, None) => kernel::align_sse::<false>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
            }
        }
    }
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        let m8 = fits_8bit(scores, gap_open, gap_extend, x_drop);
        // SAFETY: simd128 is a compile-time feature of this build.
        unsafe {
            match (tb, m8) {
                (true, Some(m)) => kernel::align_body8::<true>(
                    q,
                    s,
                    scores,
                    m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (false, Some(m)) => kernel::align_body8::<false>(
                    q,
                    s,
                    scores,
                    m,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (true, None) => kernel::align_body::<true>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
                (false, None) => kernel::align_body::<false>(
                    q,
                    s,
                    scores,
                    len1,
                    len2,
                    gap_open,
                    gap_extend,
                    x_drop,
                    check_fence,
                    sc,
                ),
            }
        }
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    {
        let _ = (
            q,
            s,
            scores,
            len1,
            len2,
            gap_open,
            gap_extend,
            x_drop,
            tb,
            check_fence,
            sc,
        );
        None
    }
}

#[cfg(test)]
mod tests {
    //! Differential test: the SIMD kernels against a plain scalar
    //! transliteration of the NCBI loops, on random related sequences.
    use super::*;

    const MININT: i32 = i32::MIN / 2;

    struct Problem {
        q: Vec<u8>,
        s: Vec<u8>,
        n: usize,
        m: Vec<i32>,
    }

    impl Problem {
        fn q(&self, a: usize) -> u8 {
            self.q[a - 1]
        }
        fn s(&self, k: usize, len2: usize, sentinel: bool) -> u8 {
            if k >= len2 && !sentinel {
                return 0;
            }
            self.s.get(k).copied().unwrap_or(0)
        }
        fn score(&self, q: u8, s: u8) -> i32 {
            self.m[q as usize * self.n + s as usize]
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-912,563-636
    // ```c
    //         for (b_index = first_b_index; b_index < b_size; b_index++) {
    //             b_ptr += b_increment;
    //             score_gap_col = score_array[b_index].best_gap;
    //             next_score = score_array[b_index].best + matrix_row[ *b_ptr ];
    //
    //             if (score < score_gap_col)
    //                 score = score_gap_col;
    // ...
    //             if (matrix_index == FENCE_SENTRY) {
    // ...
    //             script = SCRIPT_SUB;
    //             script_col = SCRIPT_EXTEND_GAP_B;
    //             script_row = SCRIPT_EXTEND_GAP_A;
    // ```
    // The plain scalar transliteration of the C loops that the differential test compares the vector
    // kernels with.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:835-959
    // (Blast_SemiGappedAlign, score-only) and :500-727 (ALIGN_EX).
    fn reference(
        p: &Problem,
        len1: usize,
        len2: usize,
        gap_open: i32,
        gap_extend: i32,
        x_drop: i32,
        tb: bool,
        sentinel: bool,
    ) -> (usize, usize, i32, bool, Vec<(u8, u32)>, u64) {
        let goe = gap_open + gap_extend;
        let x = x_drop.max(goe);
        let mut best = vec![MININT; len2 + 4];
        let mut best_gap = vec![MININT; len2 + 4];
        let mut rows: Vec<(usize, Vec<u8>)> = Vec::new();
        let mut score = -goe;
        best[0] = 0;
        best_gap[0] = -goe;
        let mut b_size = 1usize;
        for i in 1..=len2 {
            if score < -x {
                break;
            }
            best[i] = score;
            best_gap[i] = score - goe;
            score -= gap_extend;
            b_size = i + 1;
        }
        rows.push((0, vec![SCRIPT_GAP_IN_A; b_size]));
        let (mut best_score, mut a_off, mut b_off, mut first) = (0i32, 0usize, 0usize, 0usize);
        let mut fence_hit = false;
        let mut cells = 0u64;
        for a in 1..=len1 {
            let orig = first;
            let mut row = Vec::new();
            let qc = p.q(a);
            let mut score_val = MININT;
            let mut gap_row = MININT;
            let mut last = first;
            cells += (b_size - first) as u64;
            for b in orig..b_size {
                let sc = p.s(b, len2, sentinel);
                if tb && sc == FENCE_SENTRY {
                    fence_hit = true;
                    break;
                }
                let gap_col = best_gap[b];
                let next = best[b].wrapping_add(p.score(qc, sc));
                let mut script = SCRIPT_SUB;
                if score_val < gap_col {
                    script = SCRIPT_GAP_IN_B;
                    score_val = gap_col;
                }
                if score_val < gap_row {
                    script = SCRIPT_GAP_IN_A;
                    score_val = gap_row;
                }
                if best_score.wrapping_sub(score_val) > x {
                    if b == first {
                        first += 1;
                    } else {
                        best[b] = MININT;
                    }
                } else {
                    last = b;
                    if score_val > best_score {
                        best_score = score_val;
                        a_off = a;
                        b_off = b;
                    }
                    gap_row -= gap_extend;
                    let ext = gap_col - gap_extend;
                    let open = score_val - goe;
                    if ext < open {
                        best_gap[b] = open;
                    } else {
                        best_gap[b] = ext;
                        script += SCRIPT_EXTEND_GAP_B;
                    }
                    if gap_row < open {
                        gap_row = open;
                    } else {
                        script += SCRIPT_EXTEND_GAP_A;
                    }
                    best[b] = score_val;
                }
                score_val = next;
                row.push(script);
            }
            if first == b_size || fence_hit {
                rows.push((orig, row));
                break;
            }
            if last < b_size - 1 {
                b_size = last + 1;
            } else {
                while gap_row >= best_score - x && b_size <= len2 {
                    best[b_size] = gap_row;
                    best_gap[b_size] = gap_row - goe;
                    gap_row -= gap_extend;
                    let idx = b_size - orig;
                    if idx < row.len() {
                        row[idx] = SCRIPT_GAP_IN_A;
                    } else {
                        row.resize(idx, SCRIPT_GAP_IN_A);
                        row.push(SCRIPT_GAP_IN_A);
                    }
                    b_size += 1;
                }
            }
            rows.push((orig, row));
            if b_size <= len2 {
                best[b_size] = MININT;
                best_gap[b_size] = MININT;
                b_size += 1;
            }
        }
        let mut ops = Vec::new();
        if !tb {
            return (a_off, b_off, best_score, false, ops, cells);
        }
        if best_score <= 0 {
            return (0, 0, 0, fence_hit, ops, cells);
        }
        if fence_hit {
            return (a_off, b_off, best_score, true, ops, cells);
        }
        let (mut a, mut b, mut script) = (a_off, b_off, SCRIPT_SUB);
        while a > 0 || b > 0 {
            let (start, row) = &rows[a];
            let next = row[b - start];
            match script & SCRIPT_OP_MASK {
                x if x == SCRIPT_GAP_IN_A => {
                    script = next & SCRIPT_OP_MASK;
                    if next & SCRIPT_EXTEND_GAP_A != 0 {
                        script = SCRIPT_GAP_IN_A;
                    }
                }
                x if x == SCRIPT_GAP_IN_B => {
                    script = next & SCRIPT_OP_MASK;
                    if next & SCRIPT_EXTEND_GAP_B != 0 {
                        script = SCRIPT_GAP_IN_B;
                    }
                }
                _ => script = next & SCRIPT_OP_MASK,
            }
            let op = script & SCRIPT_OP_MASK;
            if op == SCRIPT_GAP_IN_A {
                b = b.saturating_sub(1);
            } else if op == SCRIPT_GAP_IN_B {
                a = a.saturating_sub(1);
            } else {
                a = a.saturating_sub(1);
                b = b.saturating_sub(1);
            }
            push_op(&mut ops, op);
        }
        (a_off, b_off, best_score, false, ops, cells)
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn random_problem(rng: &mut Rng, n: usize, scale: i32, len: usize) -> Problem {
        // random symmetric-ish matrix with a positive diagonal
        let mut m = vec![0i32; n * n];
        for i in 0..n {
            for j in 0..n {
                let base = if i == j {
                    2 + rng.below(8) as i32
                } else {
                    -(rng.below(5) as i32)
                };
                let noise = if scale > 1 {
                    rng.below(scale as u64 / 2 + 1) as i32 - scale / 4
                } else {
                    0
                };
                m[i * n + j] = base * scale + noise;
            }
        }
        // one "minus infinity" row/column, like the gap letter
        for i in 0..n {
            m[i] = i16::MIN as i32;
            m[i * n] = i32::MIN / 2;
        }
        let letters = (n - 1) as u64;
        let q: Vec<u8> = (0..len).map(|_| 1 + rng.below(letters) as u8).collect();
        let (ps, pi) = (rng.unit() * 0.8, rng.unit() * 0.05);
        let mut s = Vec::new();
        let mut i = 0;
        while i < q.len() {
            let u = rng.unit();
            if u < pi {
                for _ in 0..=rng.below(5) {
                    s.push(1 + rng.below(letters) as u8);
                }
            } else if u < 2.0 * pi {
                i += 1 + rng.below(5) as usize;
                continue;
            }
            s.push(if rng.unit() < ps {
                1 + rng.below(letters) as u8
            } else {
                q[i]
            });
            i += 1;
        }
        if rng.below(6) == 0 && !s.is_empty() {
            let k = rng.below(s.len() as u64) as usize;
            s[k] = 0;
        }
        Problem { q, s, n, m }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:736-741,766-770
    // ```c
    // Blast_SemiGappedAlign(const Uint1* A, const Uint1* B, Int4 M, Int4 N,
    //    Int4* a_offset, Int4* b_offset, Boolean score_only,
    // ...
    //     if (!score_only) {
    //         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
    // ```
    // Random problems (matrices of 16 and 28 letters, several gap costs and X-drop values,
    // some with a fence sentinel) are run through `reference` and `xdrop_align`. The
    // offsets, score, fence flag, edit script and cell count are compared.
    #[test]
    fn simd_kernels_match_scalar_loops_on_random_problems() {
        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1500);
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut sc = XdropScratch::new();
        let (mut handled, mut fences, mut ops_total) = (0usize, 0usize, 0usize);
        for case in 0..cases {
            let n = if rng.below(2) == 0 { 16 } else { 28 };
            let scale = if rng.below(3) == 0 { 32 } else { 1 };
            let max_len = if rng.below(5) == 0 { 2500 } else { 300 };
            let len = 1 + rng.below(max_len) as usize;
            let mut p = random_problem(&mut rng, n, scale, len);
            let (go, ge, x) = match rng.below(6) {
                0 => (11 * scale, scale, 38 * scale),
                1 => (11 * scale, scale, 64 * scale + rng.below(5) as i32),
                2 => (5 * scale, 2 * scale, 33 * scale),
                3 => (5 * scale, 2 * scale, 110 * scale),
                4 => (
                    rng.below(12) as i32 * scale,
                    (1 + rng.below(3) as i32) * scale,
                    (3 + rng.below(150) as i32) * scale,
                ),
                _ => (9 * scale, scale, 3 * scale),
            };
            let len1 = p.q.len();
            let mut len2 = p.s.len();
            if len2 == 0 {
                continue;
            }
            for tb in [false, true] {
                let mut sentinel = false;
                if tb {
                    match rng.below(6) {
                        0 => {
                            p.s.push(FENCE_SENTRY);
                            sentinel = true;
                        }
                        1 if len2 > 3 => {
                            let k = 1 + rng.below(len2 as u64 - 1) as usize;
                            p.s[k] = FENCE_SENTRY;
                        }
                        2 if len2 > 3 => {
                            len2 = 1 + rng.below(len2 as u64 - 1) as usize;
                            sentinel = rng.below(2) == 0;
                        }
                        _ => {}
                    }
                }
                let expected = reference(&p, len1, len2, go, ge, x, tb, sentinel);
                let rows = RowSeq::Bytes {
                    data: &p.q,
                    base: -1,
                    step: 1,
                };
                let cols = ColSeq {
                    data: &p.s,
                    base: 0,
                    step: 1,
                    zero_from: if sentinel { usize::MAX } else { len2 },
                };
                let scores = Scores::Flat { data: &p.m, n };
                let Some(r) = xdrop_align(
                    &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut sc,
                ) else {
                    continue;
                };
                handled += 1;
                fences += usize::from(r.fence_hit);
                ops_total += sc.ops.len();
                assert_eq!(
                    (r.a_offset, r.b_offset, r.score, r.fence_hit),
                    (expected.0, expected.1, expected.2, expected.3),
                    "case {case} tb={tb} n={n} scale={scale} go={go} ge={ge} x={x} len1={len1} len2={len2}"
                );
                assert_eq!(sc.ops, expected.4, "case {case} traceback");
                if !r.fence_hit {
                    assert_eq!(r.cells, expected.5, "case {case} cell count");
                }
            }
        }
        // On targets without the SIMD kernel nothing is handled.
        if cfg!(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        )) {
            assert!(
                handled > cases,
                "the SIMD kernel handled {handled} of {cases} cases"
            );
            assert!(fences > 0 && ops_total > 0);
        }
    }
}
