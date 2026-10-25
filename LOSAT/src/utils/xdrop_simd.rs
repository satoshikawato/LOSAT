//! Row-wise SIMD X-drop dynamic programming
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
//! one byte per base (traceback). The vector kernels below run for all of these.
//! (The scalar ports stay for the calls the kernels decline.)
//!
//! How the module relates to that code. Each row of the DP is the same recurrence as the
//! cited loops. The diagonal term (`next_score`), the column gap (`best_gap`) and the
//! X-drop test are evaluated lane by lane, in blocks of 8 (16-bit) or 16 (8-bit) cells.
//! The row gap (`score_gap_row`), which the C loop carries from cell to cell, is a decayed
//! prefix maximum over the lanes. The script bits have the values of the NCBI `SCRIPT_*`
//! constants (or are the four-bit codes of them, `LOSAT_X_DPNIBBLE`), and the traceback walk is
//! the C walk. The argument above is why the cells that
//! are kept have the C values. The random test
//! (`simd_kernels_match_scalar_loops_on_random_problems`) compares results with the
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

// No NCBI counterpart: the traceback rows at four bits per cell (`LOSAT_X_DPNIBBLE`). NCBI keeps one
// byte of `SCRIPT_*` bits per cell (`edit_script`). With the switch on, the kernels store a code of
// four bits instead, two cells per byte (cell 2k in the low half of byte k of the row, cell 2k + 1
// in the high half), and the walk follows the codes (`kernel::finish`), giving each operation its
// NCBI value. Bits 0-1 hold
// the operation (0 = `SCRIPT_GAP_IN_A`, 1 = `SCRIPT_SUB`, 2 = `SCRIPT_GAP_IN_B`: the NCBI value / 3),
// bit 2 `SCRIPT_EXTEND_GAP_A`, bit 3 `SCRIPT_EXTEND_GAP_B`. The kernels compute the code with the
// same comparisons as the byte, so each cell keeps the script NCBI gives it.
const NIB_SUB: u8 = 1;
const NIB_EXTEND_GAP_A: u8 = 4;
const NIB_EXTEND_GAP_B: u8 = 8;

/// No NCBI counterpart: the values a kernel adds up for a cell's script (`SCRIPT_SUB`,
/// `SCRIPT_EXTEND_GAP_A`, `SCRIPT_EXTEND_GAP_B`), as NCBI bytes or as four-bit codes (`NIB`).
/// The kernels read them in `const` blocks: a branch on `NIB` that the MIR inliner copies into a
/// kernel changes its register allocation (one byte per cell: +0.85 % instructions at one
/// thread).
#[inline(always)]
pub(crate) const fn script_codes(nib: bool) -> (u8, u8, u8) {
    if nib {
        (NIB_SUB, NIB_EXTEND_GAP_A, NIB_EXTEND_GAP_B)
    } else {
        (SCRIPT_SUB, SCRIPT_EXTEND_GAP_A, SCRIPT_EXTEND_GAP_B)
    }
}

/// No NCBI counterpart: the bytes a traceback row of `cells` cells takes (two cells per byte with
/// `NIB`, rounded up), without a branch (`script_codes`).
#[inline(always)]
pub(crate) const fn trace_bytes(nib: bool, cells: usize) -> usize {
    (cells + nib as usize) >> nib as usize
}

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

// No NCBI counterpart: `Rows` and the only code that sets its fields.
mod rows {
    /// No NCBI counterpart (memory placement only): a row of the vector kernels whose first
    /// element lies on a 64-byte boundary. The 256-bit kernels store 32 bytes at `PAD + first`
    /// plus multiples of 32 bytes and load the same places in the next row; from a buffer at
    /// glibc's 16-byte alignment half of those stores could cross a cache line, and a split store
    /// does not forward to the load after it (an offset set by the heap's history: one more
    /// allocation elsewhere made blastn's 16-thread search 2 % slower with three times the blocked
    /// store forwards). The values and the lengths are those of the `Vec` before. The fields are
    /// private to `rows`, so that only `new` and `resize` set `ptr` and `len`.
    pub(super) struct Rows<T> {
        v: Vec<T>,
        /// `v`'s element on the boundary (`v` is never pushed to, so its buffer does not move).
        ptr: *mut T,
        len: usize,
    }

    // SAFETY: `ptr` points into `v` (or dangles with `len` 0) and is sent and shared only together
    // with `v`, so `Rows<T>` is `Send` / `Sync` when `Vec<T>` is.
    unsafe impl<T: Send> Send for Rows<T> {}
    unsafe impl<T: Sync> Sync for Rows<T> {}

    impl<T: Copy + Default> Rows<T> {
        pub(super) fn new(len: usize) -> Self {
            let mut r = Rows {
                v: Vec::new(),
                ptr: std::ptr::NonNull::dangling().as_ptr(),
                len: 0,
            };
            r.resize(len, T::default());
            r
        }

        /// `Vec::resize`: the rows only grow (`grow`).
        pub(super) fn resize(&mut self, len: usize, value: T) {
            if len <= self.len {
                self.len = len;
                return;
            }
            let slack = 64 / std::mem::size_of::<T>();
            let mut v = vec![value; len + slack];
            // `align_offset` gives `usize::MAX` only for a pointer it cannot align, which `v`'s
            // buffer (aligned for `T`) never is; if it ever did, the rows would start at `v`'s first
            // element, which changes only the speed.
            let off = v.as_ptr().align_offset(64);
            let off = if off <= slack { off } else { 0 };
            v[off..off + self.len].copy_from_slice(self);
            // SAFETY: `off + len <= v.len()`. `Vec::as_mut_ptr` makes no reference to `v`'s
            // elements, so the pointer may reach all of them.
            let ptr = unsafe { v.as_mut_ptr().add(off) };
            *self = Rows { v, ptr, len };
        }
    }

    impl<T> Rows<T> {
        /// `Vec::as_ptr`: the row's first element without a slice made first, so that a pointer a
        /// kernel takes stays valid when the row is then read or written through `Deref` /
        /// `DerefMut` or another `as_ptr` / `as_mut_ptr` call, as Miri accepts for a `Vec` (a
        /// pointer from the slice's method comes from a slice that the next `deref` or `deref_mut`
        /// of the row can invalidate).
        #[inline(always)]
        pub(super) fn as_ptr(&self) -> *const T {
            self.ptr
        }

        /// `Vec::as_mut_ptr`; see `as_ptr`.
        #[inline(always)]
        pub(super) fn as_mut_ptr(&mut self) -> *mut T {
            self.ptr
        }
    }

    impl<T> std::ops::Deref for Rows<T> {
        type Target = [T];
        #[inline(always)]
        fn deref(&self) -> &[T] {
            // SAFETY: `ptr .. ptr + len` lies in `v`, or `len` is 0 and `ptr` is dangling (`new`);
            // only `new` and `resize` set them.
            unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
        }
    }

    impl<T> std::ops::DerefMut for Rows<T> {
        #[inline(always)]
        fn deref_mut(&mut self) -> &mut [T] {
            // SAFETY: as `deref`, and `&mut self` borrows `v` uniquely.
            unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
        }
    }

    impl<T: PartialEq> PartialEq for Rows<T> {
        fn eq(&self, other: &Self) -> bool {
            **self == **other
        }
    }

    impl<T: std::fmt::Debug> std::fmt::Debug for Rows<T> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            (**self).fmt(f)
        }
    }
}
use rows::Rows;

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
    ha: Rows<i16>,
    hb: Rows<i16>,
    e: Rows<i16>,
    // 8-bit path
    h8a: Rows<u8>,
    h8b: Rows<u8>,
    e8: Rows<u8>,
    dirty8: usize,
    sb: Rows<u8>,
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
    /// The column of lane 0 of the rows (`ha` .. `e8`, `sb`): a kernel moves its lanes down when
    /// the band moves right (`x_shift`), so the rows hold the band and not every column before it.
    x_sh: usize,
    /// The band start, in lanes, at which a kernel moves its rows (`LOSAT_X_DPWINDOW`); `usize::MAX`
    /// when they hold every column.
    x_win_at: usize,
    /// Whether the kernels write `trace` as four-bit codes, two cells per byte (`LOSAT_X_DPNIBBLE`;
    /// `NIB_SUB`; a search's scratch follows `x_dp_nibble_here`); the walk reads the rows the same
    /// way.
    x_nib: bool,
}

impl XdropScratch {
    /// No NCBI counterpart: an empty scratch when the target has no vector kernel, so that a
    /// search there allocates nothing for it; it does not change any value NCBI computes.
    ///
    /// The scratch of a search: empty unless the vector kernels are compiled for the target.
    pub(crate) fn for_search() -> Self {
        if KERNEL_COMPILED {
            let mut sc = Self::new();
            sc.x_nib = x_dp_nibble_here();
            return sc;
        }
        Self {
            ha: Rows::new(0),
            hb: Rows::new(0),
            e: Rows::new(0),
            h8a: Rows::new(0),
            h8b: Rows::new(0),
            e8: Rows::new(0),
            dirty8: 0,
            sb: Rows::new(0),
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
            x_sh: 0,
            x_win_at: x_dp_window_at(),
            x_nib: x_dp_nibble(),
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
            ha: Rows::new(1024),
            hb: Rows::new(1024),
            e: Rows::new(1024),
            h8a: Rows::new(1024),
            h8b: Rows::new(1024),
            e8: Rows::new(1024),
            dirty8: 0,
            sb: Rows::new(1024),
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
            x_sh: 0,
            x_win_at: x_dp_window_at(),
            x_nib: x_dp_nibble(),
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
        self.x_sh = 0;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:641-649
    // ```c
    //         if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {
    //
    //             gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
    //                                           2 * gap_align->dp_mem_alloc);
    // ```
    // NCBI's `score_array` is indexed by column and grows to the last column of the band, like
    // these rows before `LOSAT_X_DPWINDOW`. Here a kernel whose band start `first` reached
    // `x_win_at` moves the lanes from `d` on down to lane 0 (`d` a multiple of 64 with at least 64
    // lanes left below `first`) and lowers every column it keeps by `d`; `x_sh` adds `d` back
    // where a column leaves the kernel (the best cell, the start of a traceback row). The lanes
    // keep their values, the `d` lanes freed at the top are zero as lanes never written are, and
    // the lanes below `first` are not read, so the DP is the same; only the rows' length changes.
    // `wide`: the 16-bit rows (`ha`, `hb`, `e`), else the 8-bit ones (`h8a`, `h8b`, `e8`).
    #[cold]
    #[inline(never)]
    fn x_shift(&mut self, d: usize, wide: bool) {
        crate::utils::xreach::x_reach!("LOSAT_X_DPWINDOW", true);
        if wide {
            let n = self.dirty + 2 * PAD + 32;
            x_shift_lanes(&mut self.ha, d, n);
            x_shift_lanes(&mut self.hb, d, n);
            x_shift_lanes(&mut self.e, d, n);
            self.dirty = self.dirty.saturating_sub(d);
        } else {
            let n = self.dirty8 + 2 * PAD + 32;
            x_shift_lanes(&mut self.h8a, d, n);
            x_shift_lanes(&mut self.h8b, d, n);
            x_shift_lanes(&mut self.e8, d, n);
            self.dirty8 = self.dirty8.saturating_sub(d);
        }
        // `sb[PAD + k + 1]` holds column `k + x_sh`; lanes past `sb_filled` are filled before use
        let n = (PAD + self.sb_filled + 1 + 32).min(self.sb.len());
        if PAD + d < n {
            self.sb.copy_within(PAD + d..n, PAD);
        }
        self.sb_filled = self.sb_filled.saturating_sub(d);
        if self.fence_k != usize::MAX {
            self.fence_k = self.fence_k.saturating_sub(d);
        }
        if self.bad_k != usize::MAX {
            self.bad_k = self.bad_k.saturating_sub(d);
        }
        self.x_sh += d;
    }
}

// No NCBI counterpart: lanes `PAD + d .. n` of one row move to `PAD ..`, and the `d` lanes at the
// top become zero (`x_shift`).
fn x_shift_lanes<T: Copy + Default>(v: &mut [T], d: usize, n: usize) {
    let n = n.min(v.len());
    if PAD + d < n {
        v.copy_within(PAD + d..n, PAD);
        v[n - d..n].fill(T::default());
    } else if PAD < n {
        v[PAD..n].fill(T::default());
    }
}

// PROBE (probe/b16e-pool, line B s18; never a PR as is): the trace and rows of a long traceback
// come from a pool the scratches share. A scratch whose trace passes `x_pool_from()` bytes while
// it grows takes the largest idle pair (its rows so far are copied over); after the traceback
// call (`xdrop_align`) a pair above that size goes back, at most `x_pool_idle()` idle (the
// smallest freed). Placement of memory only. PROBE_B16E3_FROM=bytes (0 = off; default 1 MiB),
// PROBE_B16E3_IDLE=n (default 4).
type XTracePair = (Vec<u8>, Vec<(u32, u32, u32)>);
static X_TRACE_POOL: std::sync::Mutex<Vec<XTracePair>> = std::sync::Mutex::new(Vec::new());

fn x_probe_env(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn x_pool_from() -> usize {
    use std::sync::OnceLock;
    static FROM: OnceLock<usize> = OnceLock::new();
    *FROM.get_or_init(|| match x_probe_env("PROBE_B16E3_FROM", 1 << 20) {
        0 => usize::MAX,
        n => n,
    })
}

fn x_pool_idle() -> usize {
    use std::sync::OnceLock;
    static IDLE: OnceLock<usize> = OnceLock::new();
    *IDLE.get_or_init(|| x_probe_env("PROBE_B16E3_IDLE", 4))
}

fn x_pool_lock() -> std::sync::MutexGuard<'static, Vec<XTracePair>> {
    X_TRACE_POOL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

thread_local! {
    // The calling thread's own small trace + rows while it holds a pooled pair (one xdrop_align
    // call: taken in `x_trace_grow`, put back in `XdropScratch::x_trace_release`).
    static X_OWN: std::cell::RefCell<Option<XTracePair>> = const { std::cell::RefCell::new(None) };
}

/// Room for `n` more trace bytes (the kernels' row reserve); only the two buffers escape, as
/// with `Vec::reserve`.
#[inline(always)]
pub(crate) fn x_trace_reserve(trace: &mut Vec<u8>, rows: &mut Vec<(u32, u32, u32)>, n: usize) {
    if trace.capacity() - trace.len() < n {
        x_trace_grow(trace, rows, n);
    }
}

#[cold]
#[inline(never)]
fn x_trace_grow(trace: &mut Vec<u8>, rows: &mut Vec<(u32, u32, u32)>, n: usize) {
    let from = x_pool_from();
    if trace.len() + n > from && trace.capacity() <= from {
        let taken = {
            let mut pool = x_pool_lock();
            let best = pool
                .iter()
                .enumerate()
                .max_by_key(|(_, p)| p.0.capacity())
                .map(|(i, _)| i);
            best.map(|i| pool.swap_remove(i))
        };
        if let Some((mut t, mut r)) = taken {
            t.clear();
            t.extend_from_slice(trace);
            r.clear();
            r.extend_from_slice(rows);
            std::mem::swap(trace, &mut t);
            std::mem::swap(rows, &mut r);
            X_OWN.with(|own| *own.borrow_mut() = Some((t, r)));
        }
    }
    trace.reserve(n);
}

impl XdropScratch {
    /// After a traceback call: a pooled pair goes back (the thread's own buffers return), and a
    /// trace of a pool thread that grew past the pool size by itself goes to the pool. A
    /// one-thread search (no pool) keeps its trace: it has nothing to share it with.
    #[inline]
    pub(crate) fn x_trace_release(&mut self) {
        if self.trace.capacity() > x_pool_from() {
            self.x_trace_give_back();
        }
    }

    #[cold]
    #[inline(never)]
    fn x_trace_give_back(&mut self) {
        let own = X_OWN.with(|own| own.borrow_mut().take());
        let pair = match own {
            Some((t, r)) => (
                std::mem::replace(&mut self.trace, t),
                std::mem::replace(&mut self.rows, r),
            ),
            None if x_in_pool() => (
                std::mem::take(&mut self.trace),
                std::mem::take(&mut self.rows),
            ),
            None => return,
        };
        let dropped = {
            let mut pool = x_pool_lock();
            pool.push(pair);
            if pool.len() > x_pool_idle() {
                let small = pool
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, p)| p.0.capacity())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                Some(pool.swap_remove(small))
            } else {
                None
            }
        };
        drop(dropped);
    }
}

/// No NCBI counterpart: `LOSAT_X_DPWINDOW` (strict set, a value switch): the band start, in
/// lanes, at which the vector kernels move their rows down (`XdropScratch::x_shift`); built-in
/// 1024, at least 128. Off: the rows hold every column up to the band's end, as before.
pub(crate) fn x_dp_window_at() -> usize {
    use std::sync::OnceLock;
    static AT: OnceLock<usize> = OnceLock::new();
    *AT.get_or_init(|| {
        crate::utils::xswitch::value("LOSAT_X_DPWINDOW", 128..=usize::MAX, 1024)
            .unwrap_or(usize::MAX)
    })
}

/// No NCBI counterpart: `LOSAT_X_DPNIBBLE` (strict set): the vector kernels keep the traceback rows
/// at four bits per cell, two cells per byte (`NIB_SUB`). Off: one byte per cell, as before.
pub(crate) fn x_dp_nibble() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DPNIBBLE"))
}

/// No NCBI counterpart: the rule of `LOSAT_X_DPNIBBLE` for the scratch of a search
/// (`XdropScratch::for_search`): four-bit codes when the scratch is made on a thread of a Rayon
/// pool, that is in a search on two threads or more (its pool, or the one-thread pools of the
/// blastn query-batch slots); one byte per cell on the calling thread of a one-thread search,
/// which has no pool. One thread keeps the bytes because there the codes' packing cost time
/// (Sakai x EDL933 x1.009) and the memory they save is one trace. Memory only: both layouts give
/// every script NCBI gives.
fn x_dp_nibble_here() -> bool {
    if !x_dp_nibble() {
        return false;
    }
    let pooled = x_in_pool();
    crate::utils::xreach::x_reach!("LOSAT_X_DPNIBBLE", pooled);
    pooled
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_in_pool() -> bool {
    rayon::current_thread_index().is_some()
}

#[cfg(any(
    not(feature = "parallel"),
    all(target_arch = "wasm32", not(feature = "wasm-threads"))
))]
fn x_in_pool() -> bool {
    false
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
    /// 16 four-bit codes (bytes <= 15) to 8 bytes in the low half, byte k = code 2k | code 2k+1 << 4.
    #[inline(always)]
    pub unsafe fn nibs(v: V) -> V {
        let m = _mm_maddubs_epi16(v, _mm_set1_epi16(0x1001));
        _mm_packus_epi16(m, m)
    }
    /// The low 8 bytes of `v`.
    #[inline(always)]
    pub unsafe fn store64<T>(p: *mut T, v: V) {
        _mm_storel_epi64(p as *mut V, v)
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
    /// 16 four-bit codes (bytes <= 15) to 8 bytes in the low half, byte k = code 2k | code 2k+1 << 4.
    #[inline(always)]
    pub unsafe fn nibs(v: V) -> V {
        let w = v128_or(v, u16x8_shr(v, 4));
        i8x16_shuffle::<0, 2, 4, 6, 8, 10, 12, 14, 0, 2, 4, 6, 8, 10, 12, 14>(w, w)
    }
    /// The low 8 bytes of `v`.
    #[inline(always)]
    pub unsafe fn store64<T>(p: *mut T, v: V) {
        v128_store64_lane::<0>(v, p as *mut u64)
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

// No NCBI counterpart: the same vector instruction wrappers for aarch64 NEON (Advanced SIMD, part
// of every aarch64 target Rust builds for, so no run-time check). One 128-bit register type for
// all lane widths; the reinterpretations between them emit no instruction. Where NEON has no
// single instruction for an x86 one, the comment gives the sequence; the result is the wasm result
// for every input (`pick` and `lut32` included) and the x86 result where x86 and wasm agree.
#[cfg(target_arch = "aarch64")]
mod arch {
    use core::arch::aarch64::*;
    pub type V = uint8x16_t;
    #[inline(always)]
    unsafe fn s16(v: V) -> int16x8_t {
        vreinterpretq_s16_u8(v)
    }
    #[inline(always)]
    unsafe fn u16(v: V) -> uint16x8_t {
        vreinterpretq_u16_u8(v)
    }
    #[inline(always)]
    unsafe fn s8(v: V) -> int8x16_t {
        vreinterpretq_s8_u8(v)
    }
    #[inline(always)]
    unsafe fn u64(v: V) -> uint64x2_t {
        vreinterpretq_u64_u8(v)
    }
    #[inline(always)]
    pub unsafe fn zero() -> V {
        vdupq_n_u8(0)
    }
    #[inline(always)]
    pub unsafe fn splat16(x: i16) -> V {
        vreinterpretq_u8_s16(vdupq_n_s16(x))
    }
    #[inline(always)]
    pub unsafe fn splat8(x: i8) -> V {
        vreinterpretq_u8_s8(vdupq_n_s8(x))
    }
    #[inline(always)]
    pub unsafe fn load<T>(p: *const T) -> V {
        vld1q_u8(p as *const u8)
    }
    #[inline(always)]
    pub unsafe fn store<T>(p: *mut T, v: V) {
        vst1q_u8(p as *mut u8, v)
    }
    /// 16 four-bit codes (bytes <= 15) to 8 bytes in the low half, byte k = code 2k | code 2k+1 << 4.
    #[inline(always)]
    pub unsafe fn nibs(v: V) -> V {
        let w = vreinterpretq_u16_u8(v);
        let n = vmovn_u16(vsraq_n_u16::<4>(w, w));
        vcombine_u8(n, n)
    }
    /// The low 8 bytes of `v`.
    #[inline(always)]
    pub unsafe fn store64<T>(p: *mut T, v: V) {
        vst1_u8(p as *mut u8, vget_low_u8(v))
    }
    #[inline(always)]
    pub unsafe fn adds16(a: V, b: V) -> V {
        vreinterpretq_u8_s16(vqaddq_s16(s16(a), s16(b)))
    }
    #[inline(always)]
    pub unsafe fn add16(a: V, b: V) -> V {
        vreinterpretq_u8_s16(vaddq_s16(s16(a), s16(b)))
    }
    #[inline(always)]
    pub unsafe fn subsu16(a: V, b: V) -> V {
        vreinterpretq_u8_u16(vqsubq_u16(u16(a), u16(b)))
    }
    #[inline(always)]
    pub unsafe fn max16(a: V, b: V) -> V {
        vreinterpretq_u8_s16(vmaxq_s16(s16(a), s16(b)))
    }
    #[inline(always)]
    pub unsafe fn gt16(a: V, b: V) -> V {
        vreinterpretq_u8_u16(vcgtq_s16(s16(a), s16(b)))
    }
    #[inline(always)]
    pub unsafe fn and(a: V, b: V) -> V {
        vandq_u8(a, b)
    }
    #[inline(always)]
    pub unsafe fn or(a: V, b: V) -> V {
        vorrq_u8(a, b)
    }
    /// `!mask & v` (`bic` takes the operands the other way round).
    #[inline(always)]
    pub unsafe fn andnot(mask: V, v: V) -> V {
        vbicq_u8(v, mask)
    }
    #[inline(always)]
    pub unsafe fn shlq16(v: V) -> V {
        vreinterpretq_u8_u64(vshlq_n_u64::<16>(u64(v)))
    }
    #[inline(always)]
    pub unsafe fn shlq32(v: V) -> V {
        vreinterpretq_u8_u64(vshlq_n_u64::<32>(u64(v)))
    }
    /// Bytes up by 2 / 4 / 8 lanes, zeros shifted in (x86 `pslldq`): the window of `0 ++ v` that
    /// starts 16 - n bytes in.
    #[inline(always)]
    pub unsafe fn bsl2(v: V) -> V {
        vextq_u8::<14>(vdupq_n_u8(0), v)
    }
    #[inline(always)]
    pub unsafe fn bsl4(v: V) -> V {
        vextq_u8::<12>(vdupq_n_u8(0), v)
    }
    #[inline(always)]
    pub unsafe fn bsl8(v: V) -> V {
        vextq_u8::<8>(vdupq_n_u8(0), v)
    }
    /// Byte shuffle; an index byte >= 16 gives 0. `tbl` with one table register is wasm
    /// `swizzle` for every index byte. x86 `pshufb` differs for 16..128 only.
    #[inline(always)]
    pub unsafe fn pick(v: V, idx: V) -> V {
        vqtbl1q_u8(v, idx)
    }
    /// 32-entry byte table lookup; an index byte >= 32 gives 0. `tbl` with two table registers is
    /// the two `swizzle`s of wasm for every index byte. x86 can differ from both for 32 and up:
    /// those bytes (the fence sentinel 201, a letter the kernel declines, columns past the band)
    /// reach only lanes the kernels mask or do not read, so all three give the same results.
    #[inline(always)]
    pub unsafe fn lut32(t0: V, t1: V, idx: V) -> V {
        vqtbl2q_u8(uint8x16x2_t(t0, t1), idx)
    }
    /// Lane 7 (16-bit) to every lane (`mask` unused).
    #[inline(always)]
    pub unsafe fn bcast_last16(v: V, _mask: V) -> V {
        vreinterpretq_u8_u16(vdupq_laneq_u16::<7>(u16(v)))
    }
    /// Lane 3 (16-bit) to lanes 4..8, zero elsewhere: the x86 shuffle mask (bytes 128 or 6, 7)
    /// through `tbl`, which gives 0 for the out-of-range bytes as `pshufb` does for 128.
    #[inline(always)]
    pub unsafe fn cross16(v: V, mask: V) -> V {
        vqtbl1q_u8(v, mask)
    }
    /// Lane 15 (8-bit) to every lane (`mask` unused).
    #[inline(always)]
    pub unsafe fn bcast_last8(v: V, _mask: V) -> V {
        vdupq_laneq_u8::<15>(v)
    }
    /// Lane 7 (8-bit) to lanes 8..16, zero elsewhere (the x86 mask through `tbl`, as `cross16`).
    #[inline(always)]
    pub unsafe fn cross8(v: V, mask: V) -> V {
        vqtbl1q_u8(v, mask)
    }
    #[inline(always)]
    pub unsafe fn unpacklo8(a: V, b: V) -> V {
        vzip1q_u8(a, b)
    }
    #[inline(always)]
    pub unsafe fn unpackhi8(a: V, b: V) -> V {
        vzip2q_u8(a, b)
    }
    /// i16 to u8 with unsigned saturation: the lanes of `a`, then those of `b` (`packuswb`).
    #[inline(always)]
    pub unsafe fn packus16(a: V, b: V) -> V {
        vqmovun_high_s16(vqmovun_s16(s16(a)), s16(b))
    }
    /// The top bit of each byte, byte i at bit i (`pmovmskb`): the top bits spread over each
    /// byte and weighted 1..128 within each half; byte i of the low half and byte i of the high
    /// half are interleaved into 16-bit lane i, and one `addv` over the eight lanes sums the
    /// disjoint bits. (Of three forms timed on the arm runners, this one was fastest: 1.5-14 %
    /// less time than one `addv` per half on ubuntu-24.04-arm, 69 of 70 paired rounds.)
    #[inline(always)]
    pub unsafe fn mask8(v: V) -> u32 {
        const W: [u8; 16] = [1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128];
        let m = vandq_u8(
            vreinterpretq_u8_s8(vshrq_n_s8::<7>(s8(v))),
            vld1q_u8(W.as_ptr()),
        );
        let z = vzip1q_u8(m, vextq_u8::<8>(m, m));
        vaddvq_u16(vreinterpretq_u16_u8(z)) as u32
    }
    #[inline(always)]
    pub unsafe fn lane0(v: V) -> i32 {
        vgetq_lane_s16::<0>(s16(v)) as i32
    }
    #[inline(always)]
    pub unsafe fn adds8(a: V, b: V) -> V {
        vreinterpretq_u8_s8(vqaddq_s8(s8(a), s8(b)))
    }
    #[inline(always)]
    pub unsafe fn add8(a: V, b: V) -> V {
        vaddq_u8(a, b)
    }
    #[inline(always)]
    pub unsafe fn subsu8(a: V, b: V) -> V {
        vqsubq_u8(a, b)
    }
    #[inline(always)]
    pub unsafe fn max8(a: V, b: V) -> V {
        vreinterpretq_u8_s8(vmaxq_s8(s8(a), s8(b)))
    }
    #[inline(always)]
    pub unsafe fn gt8(a: V, b: V) -> V {
        vcgtq_s8(s8(a), s8(b))
    }
    #[inline(always)]
    pub unsafe fn shlq8(v: V) -> V {
        vreinterpretq_u8_u64(vshlq_n_u64::<8>(u64(v)))
    }
    #[inline(always)]
    pub unsafe fn bsl1(v: V) -> V {
        vextq_u8::<15>(vdupq_n_u8(0), v)
    }
    #[inline(always)]
    pub unsafe fn lane0_u8(v: V) -> i32 {
        vgetq_lane_u8::<0>(v) as i32
    }
    #[inline(always)]
    pub unsafe fn lane15_u8(v: V) -> i32 {
        vgetq_lane_u8::<15>(v) as i32
    }
    #[inline(always)]
    pub unsafe fn lane7(v: V) -> i32 {
        vgetq_lane_s16::<7>(s16(v)) as i32
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
    all(target_arch = "wasm32", target_feature = "simd128"),
    target_arch = "aarch64"
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
        if !TB || sc.x_nib {
            align_body_rows::<TB, TB>(
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
        } else {
            align_body_bytes(
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
    }

    // No NCBI counterpart: the traceback of `align_body` with one byte per cell (`LOSAT_X_DPNIBBLE=0`). Off
    // x86_64 it is a function of its own, so that the default entry keeps one traceback body.
    #[cfg_attr(target_arch = "x86_64", inline(always))]
    #[cfg_attr(not(target_arch = "x86_64"), inline(never))]
    unsafe fn align_body_bytes(
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
        align_body_rows::<true, false>(
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

    // No NCBI counterpart: the kernel of `align_body` for one trace layout (`NIB`: four-bit codes,
    // `LOSAT_X_DPNIBBLE`).
    #[inline(always)]
    unsafe fn align_body_rows<const TB: bool, const NIB: bool>(
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
            sc.trace.resize(trace_bytes(NIB, b_size), SCRIPT_GAP_IN_A);
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
            c3: splat16(const { script_codes(NIB).0 } as i16),
            c10: splat16(const { script_codes(NIB).1 } as i16),
            c40: splat16(const { script_codes(NIB).2 } as i16),
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

        // band window (LOSAT_X_DPWINDOW, `XdropScratch::x_shift`): columns below are relative to
        // `sc.x_sh`
        let win_at = sc.x_win_at;
        let mut len2 = len2;
        for a in 1..=len1 {
            if first >= win_at {
                let d = (first - 64) & !63;
                sc.x_shift(d, true);
                first -= d;
                b_size -= d;
                len2 -= d;
                wend[0] = wend[0].saturating_sub(d);
                wend[1] = wend[1].saturating_sub(d);
            }
            sc.ensure(b_size + num_extra + 64);
            if sc.sb_filled < b_size {
                let to = (b_size + 64).min(len2 + 1);
                let sh = sc.x_sh;
                for k in sc.sb_filled..to {
                    let ch = s.get(k + sh);
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
                crate::utils::xdrop_simd::x_trace_reserve(
                    &mut sc.trace,
                    &mut sc.rows,
                    bend - first + num_extra + 96,
                );
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
                if TB && NIB {
                    store64(tr.add((b0 - first) >> 1), nibs(packus16(sc0, sc1)));
                } else if TB {
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
                b_off = rs.rec_b + sc.x_sh;
                a_off = a;
                best_rel = lane0(rs.bestv) - bias;
            }
            if TB {
                sc.trace
                    .set_len(row_off + trace_bytes(NIB, proc_end - first));
                sc.rows.push((
                    row_off as u32,
                    (first + sc.x_sh) as u32,
                    (proc_end - first) as u32,
                ));
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
        finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
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
        if !TB || sc.x_nib {
            align_body8_rows::<TB, TB>(
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
        } else {
            align_body8_bytes(
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
    }

    // No NCBI counterpart: the traceback of `align_body8` with one byte per cell (`LOSAT_X_DPNIBBLE=0`). Off
    // x86_64 it is a function of its own, so that the default entry keeps one traceback body.
    #[cfg_attr(target_arch = "x86_64", inline(always))]
    #[cfg_attr(not(target_arch = "x86_64"), inline(never))]
    unsafe fn align_body8_bytes(
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
        align_body8_rows::<true, false>(
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

    // No NCBI counterpart: the kernel of `align_body8` for one trace layout (`NIB`: four-bit codes,
    // `LOSAT_X_DPNIBBLE`).
    #[inline(always)]
    unsafe fn align_body8_rows<const TB: bool, const NIB: bool>(
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
            sc.trace.resize(trace_bytes(NIB, b_size), SCRIPT_GAP_IN_A);
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
            c3: splat8(const { script_codes(NIB).0 } as i8),
            c10: splat8(const { script_codes(NIB).1 } as i8),
            c40: splat8(const { script_codes(NIB).2 } as i8),
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

        // band window (LOSAT_X_DPWINDOW, `XdropScratch::x_shift`): columns below are relative to
        // `sc.x_sh`
        let win_at = sc.x_win_at;
        let mut len2 = len2;
        for a in 1..=len1 {
            if first >= win_at {
                let d = (first - 64) & !63;
                sc.x_shift(d, false);
                first -= d;
                b_size -= d;
                len2 -= d;
                wend[0] = wend[0].saturating_sub(d);
                wend[1] = wend[1].saturating_sub(d);
            }
            sc.ensure(b_size + num_extra + 64);
            if sc.sb_filled < b_size {
                let to = (b_size + 64).min(len2 + 1);
                let sh = sc.x_sh;
                for k in sc.sb_filled..to {
                    let ch = s.get(k + sh);
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
                crate::utils::xdrop_simd::x_trace_reserve(
                    &mut sc.trace,
                    &mut sc.rows,
                    bend - first + num_extra + 96,
                );
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
                if TB && NIB {
                    store64(tr.add((b0 - first) >> 1), nibs(sc0));
                } else if TB {
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
                b_off = rs.rec_b + sc.x_sh;
                a_off = a;
                best_rel = lane0_u8(rs.bestv) - bias;
            }
            if TB {
                sc.trace
                    .set_len(row_off + trace_bytes(NIB, proc_end - first));
                sc.rows.push((
                    row_off as u32,
                    (first + sc.x_sh) as u32,
                    (proc_end - first) as u32,
                ));
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
        finish::<TB, NIB>(sc, a_off, b_off, best, fence_hit, cells)
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
    pub(super) unsafe fn finish<const TB: bool, const NIB: bool>(
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
        // the script bits of the rows: NCBI bytes, or four-bit codes (`NIB`, `LOSAT_X_DPNIBBLE`); the
        // walk below runs on the codes themselves, so a step reads its cell with no table between
        let (op_mask, ext_a, ext_b, gap_a, gap_b) = if NIB {
            (3, NIB_EXTEND_GAP_A, NIB_EXTEND_GAP_B, 0, 2)
        } else {
            (
                SCRIPT_OP_MASK,
                SCRIPT_EXTEND_GAP_A,
                SCRIPT_EXTEND_GAP_B,
                SCRIPT_GAP_IN_A,
                SCRIPT_GAP_IN_B,
            )
        };
        let mut a_index = a_off;
        let mut b_index = b_off;
        let mut script = const { script_codes(NIB).0 };
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
            let next = if NIB {
                let pair = *sc.trace.get_unchecked(off as usize + (b_rel >> 1));
                (pair >> ((b_rel & 1) << 2)) & 15
            } else {
                *sc.trace.get_unchecked(off as usize + b_rel)
            };
            match script & op_mask {
                x if x == gap_a => {
                    script = next & op_mask;
                    if next & ext_a != 0 {
                        script = gap_a;
                    }
                }
                x if x == gap_b => {
                    script = next & op_mask;
                    if next & ext_b != 0 {
                        script = gap_b;
                    }
                }
                _ => {
                    script = next & op_mask;
                }
            }
            let op = script & op_mask;
            if op == gap_a {
                b_index = b_index.saturating_sub(1);
            } else if op == gap_b {
                a_index = a_index.saturating_sub(1);
            } else {
                a_index = a_index.saturating_sub(1);
                b_index = b_index.saturating_sub(1);
            }
            // the NCBI operation: a code's operation times SCRIPT_SUB (0, 3, 6)
            push_op(&mut sc.ops, if NIB { op * SCRIPT_SUB } else { op });
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

/// No NCBI counterpart: whether the vector kernels of this file (`mod kernel`) are compiled for
/// the target: x86_64, wasm32 built with simd128, and aarch64 (NEON). The call sites run the
/// scalar port of the X-drop where it is false, because `xdrop_align` returns `None` for every
/// call there and would only build the kernel's tables before the scalar port runs.
pub(crate) const KERNEL_COMPILED: bool = cfg!(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128"),
    target_arch = "aarch64"
));

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
/// ```c
///     if (!score_only) {
///         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
/// ```
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

// EXPERIMENT: 256-bit AVX2 kernels (16 x i16, 32 x i8), x86_64 only. A child
// module, so that it uses this module's private items (lookup tables, scratch, traceback walk)
// without changing them.
#[cfg(target_arch = "x86_64")]
#[path = "x_xdrop_avx2.rs"]
mod x_avx2;
// EXPERIMENT: the 128-bit kernels with the per-row savings of the 256-bit ones,
// over the same `arch` layer (x86_64 SSE4.1 / AVX-128 builds, wasm32 simd128 and aarch64 NEON). A
// child module like `x_avx2`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128"),
    target_arch = "aarch64"
))]
#[path = "x_xdrop_w128.rs"]
mod x_w128;
// LOSAT_X_DPNEON2 (on unless `=0`): the 256-bit kernels of `x_xdrop_avx2.rs` on aarch64 NEON, each
// vector in two registers. A child module like `x_avx2`.
#[cfg(target_arch = "aarch64")]
#[path = "x_xdrop_neon2.rs"]
mod x_neon2;
// EXPERIMENT (LOSAT_X_DPCAPTURE, LOSAT_X_DPCAPTURE_EVERY): diagnostic capture of DP problems.
#[path = "x_xdrop_capture.rs"]
mod x_capture;

/// No NCBI counterpart: reads the CPU tests and the LOSAT_X_DPCAPTURE switch once; it does not
/// change any value NCBI computes. Bit 0: try the 256-bit kernels first (set only when the CPU
/// has AVX2); bit 1: record a sample of the calls (LOSAT_X_DPCAPTURE).
/// Bit 2 (EXPERIMENT): try the 128-bit kernels of `x_xdrop_w128.rs` after the
/// 256-bit ones (set only when the CPU runs the 128-bit kernels; always on wasm32 simd128 and
/// aarch64). Bit 3 (LOSAT_X_DPNEON2, on unless `=0`, aarch64 only): try the two-register NEON
/// kernels of `x_xdrop_neon2.rs` first, where bit 0 stands on x86_64. Bit 4 (LOSAT_X_DPEMASK, with bit 0):
/// the 8-bit 256-bit kernel's `EM` instances (no E clearing at the row end, see `x_xdrop_avx2.rs`).
fn x_extra_switches() -> u8 {
    use std::sync::OnceLock;
    static ON: OnceLock<u8> = OnceLock::new();
    *ON.get_or_init(|| {
        let mut on = 0u8;
        #[cfg(target_arch = "x86_64")]
        {
            if x_avx2::available() {
                on |= 1;
                if x_avx2::emask_on() {
                    on |= 16;
                }
            }
        }
        #[cfg(target_arch = "aarch64")]
        {
            if x_neon2::switch_on() {
                on |= 8;
            }
        }
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128"),
            target_arch = "aarch64"
        ))]
        {
            if x_w128::available() {
                on |= 4;
            }
        }
        if x_capture::switch_on() {
            on |= 2;
        }
        on
    })
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
/// EXPERIMENT dispatch point of the 256-bit kernels and LOSAT_X_DPCAPTURE, ahead of the 128-bit
/// dispatch in `xdrop_align`. On a CPU with AVX2 the 256-bit kernel (`x_xdrop_avx2.rs`) runs
/// the same recurrence (`tb == false`: `Blast_SemiGappedAlign` score-only and
/// `s_BlastAlignPackedNucl`; `tb == true`: `ALIGN_EX`) and returns exactly what the 128-bit
/// kernel returns; when it returns None the 128-bit dispatch runs (and may itself leave the call
/// to the scalar port). The argument for the equality is at the head of `x_xdrop_avx2.rs`.
/// LOSAT_X_DPCAPTURE records every Nth call of a thread
/// (inputs and result); the call itself runs through `xdrop_align` unchanged. The 128-bit
/// kernels of `x_xdrop_w128.rs` are tried next (see the comment at that step).
///
/// `Some(r)`: the call is answered with `r`; `None`: run the 128-bit dispatch.
#[inline]
fn x_extra_dispatch(
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
    on: u8,
) -> Option<Option<XdropResult>> {
    if on & 2 != 0 && x_capture::take_sample() {
        let r = x_capture::without_capture(|| {
            xdrop_align(
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
            )
        });
        x_capture::record(
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
            r.as_ref(),
            &sc.ops,
        );
        return Some(r);
    }
    #[cfg(target_arch = "x86_64")]
    {
        if on & 1 != 0 {
            // SAFETY: bit 0 is set only when `x_avx2::available()` held (AVX2 detected).
            let r = unsafe {
                x_avx2::align(
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
                    on & 16 != 0,
                    sc,
                )
            };
            if r.is_some() {
                return Some(r);
            }
        }
    }
    // LOSAT_X_DPNEON2 (on unless `=0`): on aarch64 the two-register NEON kernels of
    // `x_xdrop_neon2.rs` take the place of the 256-bit kernels: the same recurrence, the same lanes per block and
    // band steps as `x_xdrop_avx2.rs`, so exactly what the 128-bit kernel returns (the argument is
    // at the head of that file). When they return None (a row-gap precondition, or a problem the
    // 128-bit kernel declines too) the next dispatch runs.
    #[cfg(target_arch = "aarch64")]
    {
        if on & 8 != 0 {
            // SAFETY: NEON is a baseline feature of the aarch64 targets.
            let r = unsafe {
                x_neon2::align(
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
                )
            };
            if r.is_some() {
                return Some(r);
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:766-770
    // ```c
    //     if (!score_only) {
    //         return ALIGN_EX(A, B, M, N, a_offset, b_offset, edit_block, gap_align,
    // ```
    // EXPERIMENT: where the 128-bit dispatch below would run (after the 256-bit
    // kernels, which keep priority), the 128-bit kernels of `x_xdrop_w128.rs` run the same
    // recurrence (`tb == false`: `Blast_SemiGappedAlign` score-only and `s_BlastAlignPackedNucl`;
    // `tb == true`: `ALIGN_EX`) with the per-row savings of the 256-bit kernels and return exactly
    // what the 128-bit kernel returns (the argument, item by item, is at the head of that file).
    // When they return None (the row-gap precondition, or a problem the 128-bit kernel declines
    // too) the 128-bit dispatch runs.
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    ))]
    {
        if on & 4 != 0 {
            // SAFETY: bit 2 is set only when `x_w128::available()` held (the 128-bit kernels' CPU
            // features on x86_64; simd128 on wasm32 and NEON on aarch64 are compile-time features).
            let r = unsafe {
                x_w128::align(
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
                )
            };
            if r.is_some() {
                return Some(r);
            }
        }
    }
    None
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
/// Dispatch point of the vector kernels for all callers. `tb == false` stands for
/// `score_only` (`Blast_SemiGappedAlign`, `s_BlastAlignPackedNucl`), `tb == true` for `ALIGN_EX`.
/// It picks the AVX or SSE build and the 8-bit or 16-bit lanes. `None` means that the caller must
/// run the scalar port of the C function; callers do so.
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
    // EXPERIMENT (256-bit and 128-bit kernels, LOSAT_X_DPCAPTURE): see `x_extra_dispatch`. This
    // is one read of a value set once; when the dispatch returns None the code below runs.
    let on = x_extra_switches();
    if on != 0 {
        if let Some(r) = x_extra_dispatch(
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
            on,
        ) {
            return r;
        }
    }
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
    #[cfg(any(
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    ))]
    {
        let m8 = fits_8bit(scores, gap_open, gap_extend, x_drop);
        // SAFETY: simd128 (wasm32) and NEON (aarch64) are compile-time features of this build.
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
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
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

    pub(super) struct Problem {
        pub(super) q: Vec<u8>,
        pub(super) s: Vec<u8>,
        pub(super) n: usize,
        pub(super) m: Vec<i32>,
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
    pub(super) fn reference(
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

    pub(super) struct Rng(pub(super) u64);
    impl Rng {
        pub(super) fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        pub(super) fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        pub(super) fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    pub(super) fn random_problem(rng: &mut Rng, n: usize, scale: i32, len: usize) -> Problem {
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
    // `Rows` (memory placement only): the first element on a 64-byte boundary at every length
    // the scratch takes, the contents kept when the rows grow, as `Vec::resize`.
    #[test]
    fn rows_start_on_64_byte_boundaries_and_keep_their_values_when_they_grow() {
        fn check<T: Copy + Default + PartialEq + std::fmt::Debug + From<u8>>() {
            let mut r = Rows::<T>::new(1024);
            assert_eq!(r.len(), 1024);
            assert_eq!(r.as_ptr() as usize % 64, 0);
            for (i, x) in r.iter_mut().enumerate() {
                *x = T::from((i % 251) as u8);
            }
            for len in [1030usize, 2100, 4096, 9001] {
                r.resize(len, T::default());
                assert_eq!(r.len(), len);
                assert_eq!(r.as_ptr() as usize % 64, 0, "len {len}");
                assert!(r[..1024]
                    .iter()
                    .enumerate()
                    .all(|(i, &x)| x == T::from((i % 251) as u8)));
                assert!(r[1024..].iter().all(|&x| x == T::default()));
            }
            assert_eq!(Rows::<T>::new(0).len(), 0);
        }
        check::<u8>();
        check::<i16>();
    }

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
            // band window (LOSAT_X_DPWINDOW): rows moved from band start 128 or 320 on, or never
            sc.x_win_at = [128, 320, usize::MAX][case % 3];
            // trace layout (LOSAT_X_DPNIBBLE): four-bit codes or one byte per cell
            sc.x_nib = case % 2 == 0;
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
            all(target_arch = "wasm32", target_feature = "simd128"),
            target_arch = "aarch64"
        )) {
            assert!(
                handled > cases,
                "the SIMD kernel handled {handled} of {cases} cases"
            );
            assert!(fences > 0 && ops_total > 0);
        }
    }

    // No NCBI counterpart: the band window of the vector kernels (`XdropScratch::x_shift`,
    // `LOSAT_X_DPWINDOW`). Long problems whose band moves right (both scores of the blastn task
    // and a protein matrix, score-only and traceback, forward columns) run with the rows moved
    // from band start 128 on and with every column held: the offsets, score, fence flag, edit
    // script and cell count are equal, equal to the scalar loops, and the moved rows stay short.
    #[test]
    fn band_window_keeps_every_result() {
        let mut rng = Rng(0x5851_F42D_4C95_7F2D);
        let (mut on, mut off) = (XdropScratch::new(), XdropScratch::new());
        (on.x_win_at, off.x_win_at) = (128, usize::MAX);
        let (mut shifted, mut handled) = (0usize, 0usize);
        for case in 0..40 {
            let n = if case % 2 == 0 { 16 } else { 28 };
            let len = 1500 + rng.below(2500) as usize;
            let p = random_problem(&mut rng, n, 1, len);
            let (go, ge, x) = if case % 4 < 2 {
                (5, 2, 110)
            } else {
                (11, 1, 38)
            };
            let (len1, len2) = (p.q.len(), p.s.len());
            if len2 == 0 {
                continue;
            }
            let rows = RowSeq::Bytes {
                data: &p.q,
                base: -1,
                step: 1,
            };
            let cols = ColSeq {
                data: &p.s,
                base: 0,
                step: 1,
                zero_from: len2,
            };
            let scores = Scores::Flat { data: &p.m, n };
            for tb in [false, true] {
                let expected = reference(&p, len1, len2, go, ge, x, tb, false);
                let a = xdrop_align(
                    &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut on,
                );
                let moved = on.x_sh;
                let rows_len = on.h8a.len().max(on.ha.len());
                let b = xdrop_align(
                    &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut off,
                );
                let (Some(a), Some(b)) = (a, b) else {
                    continue;
                };
                handled += 1;
                assert_eq!(a, b, "case {case} tb={tb} len1={len1} len2={len2}");
                assert_eq!(on.ops, off.ops, "case {case} tb={tb} traceback");
                assert_eq!(
                    (a.a_offset, a.b_offset, a.score, a.fence_hit),
                    (expected.0, expected.1, expected.2, expected.3),
                    "case {case} tb={tb} against the scalar loops"
                );
                assert_eq!(on.ops, expected.4, "case {case} tb={tb} scalar traceback");
                if moved > 0 {
                    shifted += 1;
                    assert!(
                        rows_len <= 2048,
                        "case {case}: rows of {rows_len} lanes after moving {moved} columns"
                    );
                }
            }
        }
        // every vector kernel moves its rows, the NEON2 ones (LOSAT_X_DPNEON2) too
        if cfg!(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128"),
            target_arch = "aarch64"
        )) {
            assert!(
                handled > 40 && shifted > 20,
                "handled {handled}, moved {shifted}"
            );
        }
    }

    // No NCBI counterpart: the traceback rows at four bits per cell (`LOSAT_X_DPNIBBLE`). Long
    // tracebacks (both scores of the blastn task and a protein matrix, 8-bit and 16-bit lanes) run
    // with four-bit codes and with one byte per cell: the offsets, score, fence flag, edit script
    // and cell count are equal, equal to the scalar loops; the rows have the same bands and each
    // takes half its cells' bytes, rounded up.
    #[test]
    fn nibble_trace_keeps_every_result() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        let (mut on, mut off) = (XdropScratch::new(), XdropScratch::new());
        (on.x_nib, off.x_nib) = (true, false);
        let mut handled = 0usize;
        for case in 0..40 {
            let n = if case % 2 == 0 { 16 } else { 28 };
            // scale 32 leaves the 8-bit lanes (16-bit kernels)
            let scale = if case % 3 == 2 { 32 } else { 1 };
            let len = 1500 + rng.below(2500) as usize;
            let p = random_problem(&mut rng, n, scale, len);
            let (go, ge, x) = if case % 4 < 2 {
                (5 * scale, 2 * scale, 110 * scale)
            } else {
                (11 * scale, scale, 38 * scale)
            };
            let (len1, len2) = (p.q.len(), p.s.len());
            if len2 == 0 {
                continue;
            }
            let rows = RowSeq::Bytes {
                data: &p.q,
                base: -1,
                step: 1,
            };
            let cols = ColSeq {
                data: &p.s,
                base: 0,
                step: 1,
                zero_from: len2,
            };
            let scores = Scores::Flat { data: &p.m, n };
            let expected = reference(&p, len1, len2, go, ge, x, true, false);
            let a = xdrop_align(
                &rows, &cols, &scores, len1, len2, go, ge, x, true, true, &mut on,
            );
            let b = xdrop_align(
                &rows, &cols, &scores, len1, len2, go, ge, x, true, true, &mut off,
            );
            assert_eq!(a.is_some(), b.is_some(), "case {case}: one layout declined");
            let (Some(a), Some(b)) = (a, b) else {
                continue;
            };
            handled += 1;
            assert_eq!(a, b, "case {case} len1={len1} len2={len2}");
            assert_eq!(on.ops, off.ops, "case {case} traceback");
            assert_eq!(
                (a.a_offset, a.b_offset, a.score, a.fence_hit),
                (expected.0, expected.1, expected.2, expected.3),
                "case {case} against the scalar loops"
            );
            assert_eq!(on.ops, expected.4, "case {case} scalar traceback");
            let bands = |r: &[(u32, u32, u32)]| r.iter().map(|x| (x.1, x.2)).collect::<Vec<_>>();
            assert_eq!(bands(&on.rows), bands(&off.rows), "case {case} rows");
            for w in on.rows.windows(2) {
                assert_eq!(w[1].0 - w[0].0, (w[0].2 + 1) / 2, "case {case} row bytes");
            }
            assert_eq!(
                off.trace.len(),
                off.rows.iter().map(|r| r.2 as usize).sum::<usize>()
            );
            assert_eq!(
                on.trace.len(),
                on.rows
                    .iter()
                    .map(|r| (r.2 as usize + 1) / 2)
                    .sum::<usize>(),
                "case {case} trace bytes"
            );
        }
        if cfg!(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128"),
            target_arch = "aarch64"
        )) {
            assert!(handled > 30, "handled {handled}");
        }
    }

    // No NCBI counterpart: `nibs` of the `arch` layer (`LOSAT_X_DPNIBBLE`) against its definition
    // on random codes: byte k of the low half = code 2k | code 2k+1 << 4.
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    ))]
    #[test]
    fn nibs_pairs_adjacent_codes() {
        // the x86_64 `arch` layer runs inside the kernels' SSSE3 / SSE4.1 entries (`x86_entry!`)
        #[cfg(target_arch = "x86_64")]
        if !is_x86_feature_detected!("ssse3") {
            return;
        }
        let mut rng = Rng(0x0123_4567_89AB_CDEF);
        for _ in 0..1000 {
            let mut v = [0u8; 16];
            for b in v.iter_mut() {
                *b = rng.below(16) as u8;
            }
            let mut out = [0xAAu8; 16];
            // SAFETY: SSSE3 was detected above (x86_64); simd128 and NEON are compile-time features
            // of the wasm32 and aarch64 builds; the pointers cover 16 bytes.
            unsafe {
                super::arch::store64(
                    out.as_mut_ptr(),
                    super::arch::nibs(super::arch::load(v.as_ptr())),
                );
            }
            for k in 0..8 {
                assert_eq!(out[k], v[2 * k] | v[2 * k + 1] << 4, "byte {k} of {v:?}");
            }
            assert!(
                out[8..].iter().all(|&b| b == 0xAA),
                "store64 wrote past 8 bytes"
            );
        }
    }
}

// No NCBI counterpart: the aarch64 `arch` layer against a scalar model of each operation (the x86
// and wasm definition of the wrapper; for `pick` and `lut32`, wasm's for every index byte), on
// random vectors. The kernels themselves are compared with the scalar port by `tests` above.
#[cfg(all(test, target_arch = "aarch64"))]
mod neon_arch_tests {
    use super::arch::*;

    struct Rng(u64);
    impl Rng {
        fn bytes(&mut self) -> [u8; 16] {
            let mut b = [0u8; 16];
            for x in b.iter_mut() {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                // mostly small values and values near the i8 / i16 limits
                *x = match self.0 >> 61 {
                    0 => 0,
                    1 => 0xFF,
                    2 => 0x7F,
                    3 => 0x80,
                    _ => (self.0 >> 8) as u8,
                };
            }
            b
        }
    }

    unsafe fn v(b: [u8; 16]) -> V {
        load(b.as_ptr())
    }
    unsafe fn b(x: V) -> [u8; 16] {
        let mut o = [0u8; 16];
        store(o.as_mut_ptr(), x);
        o
    }
    fn w(b: [u8; 16]) -> [i16; 8] {
        core::array::from_fn(|i| i16::from_le_bytes([b[2 * i], b[2 * i + 1]]))
    }
    fn from_w(x: [i16; 8]) -> [u8; 16] {
        let mut o = [0u8; 16];
        for (i, l) in x.iter().enumerate() {
            o[2 * i..2 * i + 2].copy_from_slice(&l.to_le_bytes());
        }
        o
    }
    fn mask(t: bool) -> i16 {
        if t {
            -1
        } else {
            0
        }
    }

    #[test]
    fn neon_arch_layer_matches_the_scalar_definitions() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        let kcross: [u8; 16] = [
            128, 128, 128, 128, 128, 128, 128, 128, 6, 7, 6, 7, 6, 7, 6, 7,
        ];
        let kcross8: [u8; 16] = [
            128, 128, 128, 128, 128, 128, 128, 128, 7, 7, 7, 7, 7, 7, 7, 7,
        ];
        for _ in 0..20_000 {
            let (x, y, z) = (rng.bytes(), rng.bytes(), rng.bytes());
            let (xw, yw) = (w(x), w(y));
            unsafe {
                let (vx, vy, vz) = (v(x), v(y), v(z));
                let m = x
                    .iter()
                    .enumerate()
                    .map(|(i, &c)| u32::from(c >> 7) << i)
                    .sum::<u32>();
                assert_eq!(mask8(vx), m, "mask8 {x:?}");
                let p: [u8; 16] = core::array::from_fn(|i| *x.get(z[i] as usize).unwrap_or(&0));
                assert_eq!(b(pick(vx, vz)), p, "pick");
                let t: Vec<u8> = x.iter().chain(y.iter()).copied().collect();
                let l: [u8; 16] = core::array::from_fn(|i| *t.get(z[i] as usize).unwrap_or(&0));
                assert_eq!(b(lut32(vx, vy, vz)), l, "lut32");
                for (n, got) in [(1, bsl1(vx)), (2, bsl2(vx)), (4, bsl4(vx)), (8, bsl8(vx))] {
                    let e: [u8; 16] = core::array::from_fn(|i| if i >= n { x[i - n] } else { 0 });
                    assert_eq!(b(got), e, "bsl{n}");
                }
                for (k, got) in [(8, shlq8(vx)), (16, shlq16(vx)), (32, shlq32(vx))] {
                    let lo = u64::from_le_bytes(x[..8].try_into().unwrap()) << k;
                    let hi = u64::from_le_bytes(x[8..].try_into().unwrap()) << k;
                    let mut e = [0u8; 16];
                    e[..8].copy_from_slice(&lo.to_le_bytes());
                    e[8..].copy_from_slice(&hi.to_le_bytes());
                    assert_eq!(b(got), e, "shlq{k}");
                }
                let e: [i16; 8] = core::array::from_fn(|i| xw[i].saturating_add(yw[i]));
                assert_eq!(b(adds16(vx, vy)), from_w(e), "adds16");
                let e: [i16; 8] = core::array::from_fn(|i| xw[i].wrapping_add(yw[i]));
                assert_eq!(b(add16(vx, vy)), from_w(e), "add16");
                let e: [i16; 8] =
                    core::array::from_fn(|i| (xw[i] as u16).saturating_sub(yw[i] as u16) as i16);
                assert_eq!(b(subsu16(vx, vy)), from_w(e), "subsu16");
                let e: [i16; 8] = core::array::from_fn(|i| xw[i].max(yw[i]));
                assert_eq!(b(max16(vx, vy)), from_w(e), "max16");
                let e: [i16; 8] = core::array::from_fn(|i| mask(xw[i] > yw[i]));
                assert_eq!(b(gt16(vx, vy)), from_w(e), "gt16");
                let e: [u8; 16] =
                    core::array::from_fn(|i| (x[i] as i8).saturating_add(y[i] as i8) as u8);
                assert_eq!(b(adds8(vx, vy)), e, "adds8");
                let e: [u8; 16] = core::array::from_fn(|i| x[i].wrapping_add(y[i]));
                assert_eq!(b(add8(vx, vy)), e, "add8");
                let e: [u8; 16] = core::array::from_fn(|i| x[i].saturating_sub(y[i]));
                assert_eq!(b(subsu8(vx, vy)), e, "subsu8");
                let e: [u8; 16] = core::array::from_fn(|i| (x[i] as i8).max(y[i] as i8) as u8);
                assert_eq!(b(max8(vx, vy)), e, "max8");
                let e: [u8; 16] =
                    core::array::from_fn(|i| if (x[i] as i8) > (y[i] as i8) { 0xFF } else { 0 });
                assert_eq!(b(gt8(vx, vy)), e, "gt8");
                let e: [u8; 16] = core::array::from_fn(|i| x[i] & y[i]);
                assert_eq!(b(and(vx, vy)), e, "and");
                let e: [u8; 16] = core::array::from_fn(|i| x[i] | y[i]);
                assert_eq!(b(or(vx, vy)), e, "or");
                let e: [u8; 16] = core::array::from_fn(|i| !x[i] & y[i]);
                assert_eq!(b(andnot(vx, vy)), e, "andnot");
                let e: [u8; 16] =
                    core::array::from_fn(|i| if i % 2 == 0 { x[i / 2] } else { y[i / 2] });
                assert_eq!(b(unpacklo8(vx, vy)), e, "unpacklo8");
                let e: [u8; 16] = core::array::from_fn(|i| {
                    if i % 2 == 0 {
                        x[8 + i / 2]
                    } else {
                        y[8 + i / 2]
                    }
                });
                assert_eq!(b(unpackhi8(vx, vy)), e, "unpackhi8");
                let e: [u8; 16] = core::array::from_fn(|i| {
                    let l = if i < 8 { xw[i] } else { yw[i - 8] };
                    l.clamp(0, 255) as u8
                });
                assert_eq!(b(packus16(vx, vy)), e, "packus16");
                assert_eq!(b(bcast_last16(vx, vz)), from_w([xw[7]; 8]), "bcast_last16");
                assert_eq!(b(bcast_last8(vx, vz)), [x[15]; 16], "bcast_last8");
                let e: [i16; 8] = core::array::from_fn(|i| if i >= 4 { xw[3] } else { 0 });
                assert_eq!(b(cross16(vx, v(kcross))), from_w(e), "cross16");
                let e: [u8; 16] = core::array::from_fn(|i| if i >= 8 { x[7] } else { 0 });
                assert_eq!(b(cross8(vx, v(kcross8))), e, "cross8");
                assert_eq!(lane0(vx), xw[0] as i32, "lane0");
                assert_eq!(lane7(vx), xw[7] as i32, "lane7");
                assert_eq!(lane0_u8(vx), x[0] as i32, "lane0_u8");
                assert_eq!(lane15_u8(vx), x[15] as i32, "lane15_u8");
                assert_eq!(b(splat16(xw[2])), from_w([xw[2]; 8]), "splat16");
                assert_eq!(b(splat8(x[3] as i8)), [x[3]; 16], "splat8");
                assert_eq!(b(zero()), [0u8; 16], "zero");
            }
        }
    }
}
