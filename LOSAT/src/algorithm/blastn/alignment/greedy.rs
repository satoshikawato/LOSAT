//! Greedy alignment algorithms for high-identity nucleotide sequences
//!
//! This module implements NCBI BLAST's greedy alignment algorithm (Zhang et al., 2000)
//! for fast alignment of high-identity sequences.
//!
//! ## NCBI BLAST Compatibility Notes
//!
//! This implementation follows NCBI BLAST's greedy_align.c and blast_gapalign.c:
//!
//! - **Non-affine greedy** (`greedy_align_one_direction_with_max_dist`):
//!   Implements BLAST_GreedyAlign with distance-based tracking and X-drop
//!
//! - **Affine greedy** (`affine_greedy_align_one_direction_with_max_dist`):
//!   Implements BLAST_AffineGreedyAlign with three-state tracking (insert/match/delete)
//!
//! Key NCBI BLAST behaviors preserved:
//! - Score normalization when reward is odd (doubles all scores)
//! - Dynamic max_dist doubling for long alignments
//! - X-drop scaled to match score units
//! - Convergence detection via diagonal bounds
//!
//! Statistics estimation:
//! - Since we don't perform full traceback, statistics (matches, mismatches,
//!   gap_opens, gap_letters) are estimated from the alignment endpoints
//! - NCBI BLAST's edit script traceback could be added for exact statistics if needed

use std::cell::RefCell;

use super::super::constants::{
    GREEDY_MAX_COST, GREEDY_MAX_COST_FRACTION, INVALID_DIAG, INVALID_OFFSET,
};
use super::super::sequence_compare::{find_first_mismatch, find_first_mismatch_ex};
use super::utilities::gdb3;
use crate::common::GapEditOp;
use crate::core::blast_encoding::COMPRESSION_RATIO;

/// Thread-local memory pool for non-affine greedy alignment.
/// This avoids per-call allocation overhead by reusing memory across calls.
/// Similar to NCBI BLAST's SGreedyAlignMem structure.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/greedy_align.h:88-99
/// ```c
/// typedef struct SGreedyAlignMem {
///    Int4 max_dist;
///    Int4 xdrop;
///    Int4** last_seq2_off;
///    Int4* max_score;
///    SGreedyOffset** last_seq2_off_affine;
///    Int4* diag_bounds;
///    SMBSpace* space;
/// } SGreedyAlignMem;
/// ```
pub struct GreedyAlignMem {
    /// Two rows for last_seq2_off (we swap between them)
    last_seq2_off_a: Vec<i32>,
    last_seq2_off_b: Vec<i32>,
    /// Array of maximum scores at each distance
    max_score: Vec<i32>,
    /// Current allocated size
    allocated_size: usize,
}

impl GreedyAlignMem {
    fn new() -> Self {
        Self {
            last_seq2_off_a: Vec::new(),
            last_seq2_off_b: Vec::new(),
            max_score: Vec::new(),
            allocated_size: 0,
        }
    }

    /// Ensure the memory pool has enough capacity for the given array_size and max_score_size
    fn ensure_capacity(&mut self, array_size: usize, max_score_size: usize) {
        if self.allocated_size < array_size {
            self.last_seq2_off_a.resize(array_size, -2);
            self.last_seq2_off_b.resize(array_size, -2);
            self.allocated_size = array_size;
        }
        if self.max_score.len() < max_score_size {
            self.max_score.resize(max_score_size, 0);
        }
    }

    /// Reset arrays to initial state (fill with sentinel values)
    /// Uses slice::fill for better performance than element-by-element loops
    fn reset(&mut self, array_size: usize, max_score_size: usize) {
        // Fill with sentinel values using slice::fill (more efficient than loops)
        let a_len = array_size.min(self.last_seq2_off_a.len());
        let b_len = array_size.min(self.last_seq2_off_b.len());
        let score_len = max_score_size.min(self.max_score.len());

        self.last_seq2_off_a[..a_len].fill(-2);
        self.last_seq2_off_b[..b_len].fill(-2);
        self.max_score[..score_len].fill(0);
    }
}

thread_local! {
    /// Thread-local memory pool for non-affine greedy alignment
    static GREEDY_MEM: RefCell<GreedyAlignMem> = RefCell::new(GreedyAlignMem::new());
}

/// Scratch memory for affine greedy alignment.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/greedy_align.h:88-99
/// ```c
/// typedef struct SGreedyAlignMem {
///    SGreedyOffset** last_seq2_off_affine;
///    Int4* diag_bounds;
///    Int4* max_score;
/// } SGreedyAlignMem;
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:240-251
/// ```c
///       gamp->diag_bounds = (Int4*) calloc(2*(max_d+1+max_cost), sizeof(Int4));
///       gamp->last_seq2_off_affine = (SGreedyOffset**)
/// 	 malloc((MAX(max_d, max_cost) + 2) * sizeof(SGreedyOffset*));
///     ...
///       gamp->last_seq2_off_affine[0] = (SGreedyOffset*)
/// 	 calloc((2*max_d_1 + 6) , sizeof(SGreedyOffset) * (max_cost+1));
///       for (i = 1; i <= max_cost; i++)
/// 	 gamp->last_seq2_off_affine[i] =
/// 	    gamp->last_seq2_off_affine[i-1] + 2*max_d_1 + 6;
/// ```
/// NCBI allocates these arrays once with `calloc`, so only the parts that an extension
/// reaches take memory, and it initializes only the entries that it reads before writing
/// (`blast_affine_greedy_align`). The rows of `last_seq2_off` are `AffineRows`.
struct GreedyAffineMem {
    rows: Vec<AffineRow>,
    diag_lower: Vec<i32>,
    diag_upper: Vec<i32>,
    max_score: Vec<i32>,
}

impl GreedyAffineMem {
    fn new() -> Self {
        Self {
            rows: Vec::new(),
            diag_lower: Vec::new(),
            diag_upper: Vec::new(),
            max_score: Vec::new(),
        }
    }

    /// Makes the bound and score arrays at least as long as given. They are zeroed, so
    /// their untouched parts take no memory, as with NCBI's `calloc`.
    fn ensure_capacity(&mut self, diag_len: usize, max_score_len: usize) {
        if self.diag_lower.len() < diag_len {
            self.diag_lower = vec![0; diag_len];
            self.diag_upper = vec![0; diag_len];
        }
        if self.max_score.len() < max_score_len {
            self.max_score = vec![0; max_score_len];
        }
    }
}

const INVALID_GREEDY_OFFSET: GreedyOffset = GreedyOffset {
    insert_off: INVALID_OFFSET,
    match_off: INVALID_OFFSET,
    delete_off: INVALID_OFFSET,
};

/// The offsets of one distance for the diagonals `lower..lower + offsets.len()`.
#[derive(Default)]
struct AffineRow {
    lower: i32,
    offsets: Vec<GreedyOffset>,
}

/// The rows of `last_seq2_off` of one affine greedy extension.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:1181-1197
/// ```c
///         if (d > max_penalty) {
///             if (edit_block == NULL) {
///
///                 /* if no traceback is required, the next row of
///                    last_seq2_off can reuse previously allocated memory */
///
///                 last_seq2_off[d] = last_seq2_off[d - max_penalty - 1];
///             }
///             else {
///
///                 /* traceback requires all rows of last_seq2_off to be saved,
///                    so a new row must be allocated */
///
///                 last_seq2_off[d] = s_GetMBSpace(mem_pool,
///                                    curr_diag_upper - curr_diag_lower + 1) -
///                                    curr_diag_lower;
///             }
///         }
/// ```
/// NCBI's rows for the distances up to `max_penalty` span every diagonal, and a later
/// row reuses the row of `d - max_penalty - 1` or spans the diagonals that distance `d`
/// tests. Every row is read only inside the diagonal bounds of its distance, which lie
/// in the diagonals that the distance tests, so here every row holds only those
/// diagonals: a ring of `max_penalty + 1` rows without traceback (`ring`), one row per
/// distance with traceback.
struct AffineRows<'m> {
    rows: &'m mut Vec<AffineRow>,
    ring: Option<usize>,
}

impl AffineRows<'_> {
    fn slot(&self, d: i32) -> usize {
        let d = d as usize;
        self.ring.map_or(d, |period| d % period)
    }

    /// Starts the row of distance `d` for the diagonals `lower..=upper`.
    fn start(&mut self, d: i32, lower: i32, upper: i32) {
        let slot = self.slot(d);
        if self.rows.len() <= slot {
            self.rows.resize_with(slot + 1, AffineRow::default);
        }
        let row = &mut self.rows[slot];
        row.lower = lower;
        row.offsets.clear();
        row.offsets
            .resize((upper - lower + 1).max(0) as usize, INVALID_GREEDY_OFFSET);
    }

    fn get(&self, d: i32, k: i32) -> GreedyOffset {
        let row = &self.rows[self.slot(d)];
        usize::try_from(k - row.lower)
            .ok()
            .and_then(|index| row.offsets.get(index))
            .copied()
            .unwrap_or(INVALID_GREEDY_OFFSET)
    }

    fn get_mut(&mut self, d: i32, k: i32) -> &mut GreedyOffset {
        let slot = self.slot(d);
        let row = &mut self.rows[slot];
        &mut row.offsets[(k - row.lower) as usize]
    }
}

/// Bookkeeping structure for affine greedy alignment (NCBI BLAST's SGreedyOffset).
/// When aligning two sequences, stores the largest offset into the second sequence
/// that leads to a high-scoring alignment for a given start point, tracking
/// different path endings separately for affine gap penalties.
#[derive(Clone, Copy, Default)]
pub struct GreedyOffset {
    insert_off: i32, // Best offset for a path ending in an insertion (gap in seq2)
    match_off: i32,  // Best offset for a path ending in a match/mismatch
    delete_off: i32, // Best offset for a path ending in a deletion (gap in seq1)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_util.h:358-364 (FENCE_SENTRY)
const FENCE_SENTRY: u8 = 201;

/// Greedy seed descriptor for locating the best start point.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/greedy_align.h:102-106 (SGreedySeed)
#[derive(Clone, Copy, Default)]
struct GreedySeed {
    start_q: i32,
    start_s: i32,
    match_length: i32,
}

/// Edit script operation types used by greedy traceback.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:44-54 (EGapAlignOpType)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GapAlignOpType {
    Del,     // eGapAlignDel
    Sub,     // eGapAlignSub
    Ins,     // eGapAlignIns
    Invalid, // eGapAlignInvalid
}

impl GapAlignOpType {
    #[inline]
    fn to_gap_edit_op(self, num: i32) -> GapEditOp {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:44-54
        match self {
            GapAlignOpType::Sub => GapEditOp::Sub(num as u32),
            GapAlignOpType::Del => GapEditOp::Del(num as u32),
            GapAlignOpType::Ins => GapEditOp::Ins(num as u32),
            GapAlignOpType::Invalid => GapEditOp::Sub(num as u32),
        }
    }
}

/// Preliminary edit operation for greedy traceback.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:63-68 (GapPrelimEditScript)
#[derive(Clone, Copy, Debug)]
struct GapPrelimEditOp {
    op_type: GapAlignOpType,
    num: i32,
}

/// Offset-indexed row storage for non-affine greedy traceback.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:674-678
/// ```c
/// last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
///                          (diag_upper - diag_lower + 7) / 3);
/// last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
/// ```
#[derive(Clone, Copy)]
enum NonAffineGreedyRowStorage {
    Base(usize),
    Pool { start: usize },
}

// NCBI reference: c++/src/algo/blast/core/greedy_align.c:673-678
// last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
//                          (diag_upper - diag_lower + 7) / 3);
// last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
// Store only the biased row descriptor. The persistent scratch owns all cells.
#[derive(Clone, Copy)]
struct NonAffineGreedyRow {
    origin: i32,
    len: usize,
    storage: NonAffineGreedyRowStorage,
}

// NCBI reference: c++/src/algo/blast/core/greedy_align.c:523-526,548-550,589
// last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
// seq2_index = MAX(last_seq2_off[d - 1][k + 1], last_seq2_off[d - 1][k]) + 1;
// last_seq2_off[d][k] = seq2_index;
impl NonAffineGreedyRow {
    fn from_base(origin: i32, len: usize, base_index: usize) -> Self {
        Self {
            origin,
            len,
            storage: NonAffineGreedyRowStorage::Base(base_index),
        }
    }

    fn from_pool(origin: i32, len: usize, start: usize) -> Self {
        Self {
            origin,
            len,
            storage: NonAffineGreedyRowStorage::Pool { start },
        }
    }

    #[inline]
    fn values<'a>(&self, mem: &'a GreedyNonAffineMem) -> &'a [i32] {
        match self.storage {
            NonAffineGreedyRowStorage::Base(index) => &mem.last_seq2_off[index][..self.len],
            NonAffineGreedyRowStorage::Pool { start } => {
                &mem.traceback_pool[start..start + self.len]
            }
        }
    }

    #[inline]
    fn get(&self, mem: &GreedyNonAffineMem, diag: i32) -> i32 {
        let idx = diag - self.origin;
        if idx < 0 {
            return INVALID_OFFSET;
        }
        self.values(mem)
            .get(idx as usize)
            .copied()
            .unwrap_or(INVALID_OFFSET)
    }

    #[inline]
    fn set(&self, mem: &mut GreedyNonAffineMem, diag: i32, value: i32) {
        let idx = diag - self.origin;
        debug_assert!(idx >= 0 && (idx as usize) < self.len);
        match self.storage {
            NonAffineGreedyRowStorage::Base(index) => {
                mem.last_seq2_off[index][idx as usize] = value
            }
            NonAffineGreedyRowStorage::Pool { start } => {
                mem.traceback_pool[start + idx as usize] = value
            }
        }
    }
}

/// Persistent non-affine greedy scratch memory.
///
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/greedy_align.h:88-99
/// ```c
/// typedef struct SGreedyAlignMem {
///    Int4 max_dist;
///    Int4 xdrop;
///    Int4** last_seq2_off;
///    Int4* max_score;
///    SGreedyOffset** last_seq2_off_affine;
///    Int4* diag_bounds;
///    SMBSpace* space;
/// } SGreedyAlignMem;
/// ```
struct GreedyNonAffineMem {
    last_seq2_off: [Vec<i32>; 2],
    max_score: Vec<i32>,
    traceback_pool: Vec<i32>,
    traceback_used: usize,
}

impl GreedyNonAffineMem {
    fn new() -> Self {
        Self {
            last_seq2_off: [Vec::new(), Vec::new()],
            max_score: Vec::new(),
            traceback_pool: Vec::new(),
            traceback_used: 0,
        }
    }

    fn ensure_base_capacity(&mut self, row_len: usize) {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:211-228
        // ```c
        // gamp->last_seq2_off = (Int4**) malloc((max_d + 2) * sizeof(Int4*));
        // gamp->last_seq2_off[0] =
        //    (Int4*) malloc((max_d + max_d + 6) * sizeof(Int4) * 2);
        // gamp->last_seq2_off[1] = gamp->last_seq2_off[0] + max_d + max_d + 6;
        // ```
        // The two rolling non-affine rows live for the lifetime of
        // SGreedyAlignMem; BLAST_GreedyAlign overwrites only the cells it uses.
        for row in &mut self.last_seq2_off {
            if row.len() < row_len {
                row.resize(row_len, 0);
            }
        }
    }

    fn ensure_max_score_capacity(&mut self, len: usize) {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:257
        // ```c
        // gamp->max_score = (Int4*) malloc(sizeof(Int4) * (max_d + 1 + d_diff));
        // ```
        // BLAST_GreedyAlign clears only the leading xdrop-offset cells and
        // writes the per-distance cells it reaches.
        if self.max_score.len() < len {
            self.max_score.resize(len, 0);
        }
    }

    fn refresh_traceback_pool(&mut self) {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:67-77
        // ```c
        // static void s_RefreshMBSpace(SMBSpace* space)
        // {
        //     while (space != NULL) {
        //         space->space_used = 0;
        //         space = space->next;
        //     }
        // }
        // ```
        // Refreshing the pool rewinds allocation state but does not clear cells.
        self.traceback_used = 0;
    }

    // NCBI reference: c++/src/algo/blast/core/greedy_align.c:548-550,666-678
    // seq2_index = MAX(last_seq2_off[d - 1][k + 1], last_seq2_off[d - 1][k]) + 1;
    // if (edit_block == NULL) last_seq2_off[d + 1] = last_seq2_off[d - 1];
    // else last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool, ...);
    // Consecutive rows are disjoint. Borrow them only within one distance;
    // no reference survives the next pool allocation (which may move its Vec).
    fn row_pair_mut(
        &mut self,
        previous: NonAffineGreedyRow,
        current: NonAffineGreedyRow,
    ) -> (&mut [i32], &mut [i32]) {
        use NonAffineGreedyRowStorage::{Base, Pool};
        match (previous.storage, current.storage) {
            (Base(i), Base(j)) => {
                let (left, right) = self.last_seq2_off.split_at_mut(1);
                match (i, j) {
                    (0, 1) => (&mut left[0][..previous.len], &mut right[0][..current.len]),
                    (1, 0) => (&mut right[0][..previous.len], &mut left[0][..current.len]),
                    _ => unreachable!("greedy rolling rows must be distinct"),
                }
            }
            (Base(i), Pool { start }) => (
                &mut self.last_seq2_off[i][..previous.len],
                &mut self.traceback_pool[start..start + current.len],
            ),
            (Pool { start: p }, Pool { start: c }) => {
                let (left, right) = self.traceback_pool.split_at_mut(c);
                (&mut left[p..p + previous.len], &mut right[..current.len])
            }
            (Pool { .. }, Base(_)) => unreachable!("traceback rows never return to base storage"),
        }
    }

    fn alloc_traceback_row(&mut self, int_len: usize) -> usize {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:99-125
        // ```c
        // out_ptr = pool->space_array + pool->space_used;
        // pool->space_used += num_alloc;
        // ```
        // The C pool is allocated as SGreedyOffset cells. The caller converts
        // that allocation to Int4 storage, so Rust stores the same data as a
        // contiguous Int4 pool and preserves stale contents across refreshes.
        let start = self.traceback_used;
        let end = start + int_len;
        if self.traceback_pool.len() < end {
            self.traceback_pool.resize(end, 0);
        }
        self.traceback_used = end;
        start
    }
}

/// Preliminary edit block for greedy traceback.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:70-78 (GapPrelimEditBlock)
struct GapPrelimEditBlock {
    edit_ops: Vec<GapPrelimEditOp>,
    num_ops: usize,
    last_op: GapAlignOpType,
}

impl GapPrelimEditBlock {
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/gapinfo.c:187-198 (GapPrelimEditBlockNew)
    fn new() -> Self {
        Self {
            edit_ops: Vec::with_capacity(100),
            num_ops: 0,
            last_op: GapAlignOpType::Invalid,
        }
    }

    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/gapinfo.c:212-218 (GapPrelimEditBlockReset)
    fn reset(&mut self) {
        self.num_ops = 0;
        self.last_op = GapAlignOpType::Invalid;
        self.edit_ops.clear();
    }

    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/gapinfo.c:174-185 (GapPrelimEditBlockAdd)
    fn add(&mut self, op_type: GapAlignOpType, num_ops: i32) {
        if num_ops == 0 {
            return;
        }

        if self.last_op == op_type && self.num_ops > 0 {
            let idx = self.num_ops - 1;
            self.edit_ops[idx].num += num_ops;
            return;
        }

        self.last_op = op_type;
        self.edit_ops.push(GapPrelimEditOp {
            op_type,
            num: num_ops,
        });
        self.num_ops += 1;
    }

    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/gapinfo.c:221-231 (GapPrelimEditBlockAppend)
    #[allow(dead_code)]
    fn append(&mut self, other: &GapPrelimEditBlock) {
        for op in &other.edit_ops {
            self.add(op.op_type, op.num);
        }
    }
}

/// Scratch space for greedy gapped alignment (traceback + affine DP arrays).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:356-357
/// ```c
/// gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
/// gap_align->rev_prelim_tback = GapPrelimEditBlockNew();
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2799-2802
/// ```c
/// fwd_prelim_tback = gap_align->fwd_prelim_tback;
/// rev_prelim_tback = gap_align->rev_prelim_tback;
/// GapPrelimEditBlockReset(fwd_prelim_tback);
/// GapPrelimEditBlockReset(rev_prelim_tback);
/// ```
pub struct GreedyAlignScratch {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/greedy_align.h:88-99
    // ```c
    // typedef struct SGreedyAlignMem {
    //    Int4 max_dist;
    //    Int4 xdrop;
    //    ...
    // } SGreedyAlignMem;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2818-2829
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2851-2862
    // ```c
    // new_dist = gap_align->greedy_align_mem->max_dist * 2;
    // ...
    // gap_align->greedy_align_mem =
    //    s_BlastGreedyAlignMemAlloc(score_params, NULL, new_dist, xdrop);
    // ```
    // NCBI keeps greedy scratch max_dist in BlastGapAlignStruct and only grows it.
    max_dist: i32,
    non_affine_mem: GreedyNonAffineMem,
    affine_mem: GreedyAffineMem,
    fwd_prelim_tback: GapPrelimEditBlock,
    rev_prelim_tback: GapPrelimEditBlock,
}

impl GreedyAlignScratch {
    pub fn new() -> Self {
        Self {
            max_dist: 0,
            non_affine_mem: GreedyNonAffineMem::new(),
            affine_mem: GreedyAffineMem::new(),
            fwd_prelim_tback: GapPrelimEditBlock::new(),
            rev_prelim_tback: GapPrelimEditBlock::new(),
        }
    }
}

/// Greedy edit script container.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:56-61 (GapEditScript)
struct GapEditScript {
    op_type: Vec<GapAlignOpType>,
    num: Vec<i32>,
    size: usize,
}

impl GapEditScript {
    /// NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:88-90 (GapEditScriptNew)
    fn new(size: usize) -> Self {
        Self {
            op_type: vec![GapAlignOpType::Invalid; size],
            num: vec![0; size],
            size,
        }
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/gapinfo.h:56-61
// ```c
// typedef struct GapEditScript {
//    EGapAlignOpType* op_type;
//    Int4* num;
//    Int4 size;
// } GapEditScript;
// ```
fn format_gap_edit_script_for_trace(esp: &GapEditScript) -> String {
    let mut out = String::from("[");
    for i in 0..esp.size {
        if i > 0 {
            out.push(',');
        }
        let tag = match esp.op_type[i] {
            GapAlignOpType::Sub => 'S',
            GapAlignOpType::Del => 'D',
            GapAlignOpType::Ins => 'I',
            GapAlignOpType::Invalid => 'X',
        };
        out.push(tag);
        out.push_str(&esp.num[i].to_string());
    }
    out.push(']');
    out
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/gapinfo.h:70-78
// ```c
// typedef struct GapPrelimEditBlock {
//    GapPrelimEditScript* edit_ops;
//    Int4 num_ops_allocated;
//    Int4 num_ops;
//    EGapAlignOpType last_op;
// } GapPrelimEditBlock;
// ```
fn format_gap_prelim_edit_block_for_trace(block: &GapPrelimEditBlock) -> String {
    let mut out = String::from("[");
    for i in 0..block.num_ops {
        if i > 0 {
            out.push(',');
        }
        let op = block.edit_ops[i];
        let tag = match op.op_type {
            GapAlignOpType::Sub => 'S',
            GapAlignOpType::Del => 'D',
            GapAlignOpType::Ins => 'I',
            GapAlignOpType::Invalid => 'X',
        };
        out.push(tag);
        out.push_str(&op.num.to_string());
    }
    out.push(']');
    out
}

fn debug_greedy_traceback_enabled(q_off: usize, s_off: usize) -> bool {
    // No NCBI counterpart: debug print switch (LOSAT_DEBUG_COORDS, read once with
    // LOSAT_X_ENVCACHE); it does not change any value NCBI computes.
    let debug_all = crate::utils::xenv::debug_coords_is_some();
    let Some(filter) = crate::utils::xenv::debug_coords_start() else {
        return debug_all;
    };
    let filter = filter.to_string_lossy();
    let Some((q_filter, s_filter)) = filter.split_once(',') else {
        return false;
    };

    q_filter.trim().parse::<usize>().ok() == Some(q_off)
        && s_filter.trim().parse::<usize>().ok() == Some(s_off)
}

/// Signal that a diagonal/offset is invalid
// INVALID_OFFSET and INVALID_DIAG are imported from constants

/// Greedy alignment for high-identity nucleotide sequences.
///
/// This implements NCBI BLAST's greedy alignment algorithm (Zhang et al., 2000):
/// - Uses distance-based tracking instead of score-based
/// - Quickly scans for exact matches before doing any DP
/// - Only explores diagonals that can achieve a given distance
/// - Much faster than full DP for similar sequences
///
/// For sequences with >90% identity, this is typically 10-100x faster than full DP.
///
/// This is a wrapper that implements NCBI BLAST's dynamic max_dist doubling strategy:
/// - Start with initial max_dist based on sequence length
/// - If alignment doesn't converge, double max_dist and retry
/// - This allows finding arbitrarily long alignments
///
/// For affine gap penalties (gap_open != 0 || gap_extend != 0), uses the affine
/// greedy algorithm (NCBI BLAST's BLAST_AffineGreedyAlign).
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters)
/// Greedy alignment for high-identity nucleotide sequences.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// This avoids the need to copy and reverse sequences for left extension.
pub fn greedy_align_one_direction_ex(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
) -> (usize, usize, i32, usize, usize, usize, usize) {
    // Calculate initial max_dist like NCBI BLAST:
    // max_dist = MIN(GREEDY_MAX_COST, max(len1, len2) / GREEDY_MAX_COST_FRACTION + 1)
    let max_len = len1.max(len2);
    let mut max_dist = GREEDY_MAX_COST.min(max_len / GREEDY_MAX_COST_FRACTION + 1);

    // Retry with doubled max_dist until convergence (NCBI BLAST approach)
    loop {
        // Choose between affine and non-affine greedy based on gap penalties
        // This matches NCBI BLAST's BLAST_AffineGreedyAlign which falls back to
        // BLAST_GreedyAlign when gap_open == 0 && gap_extend == 0
        let (result, converged) = if gap_open != 0 || gap_extend != 0 {
            affine_greedy_align_one_direction_with_max_dist(
                q_seq, s_seq, reward, penalty, gap_open, gap_extend, x_drop, max_dist,
            )
        } else {
            greedy_align_one_direction_with_max_dist(
                q_seq, s_seq, len1, len2, reward, penalty, gap_open, gap_extend, x_drop, max_dist,
                reverse,
            )
        };

        if converged {
            return result;
        }

        // Double max_dist and retry (NCBI BLAST's approach)
        // NCBI reference: blast_gapalign.c:2820-2825
        // /* double the max distance */
        // new_dist = gap_align->greedy_align_mem->max_dist * 2;
        // NCBI BLAST has NO upper limit - it continues until convergence or memory failure
        max_dist *= 2;

        // NCBI BLAST has no explicit upper limit on max_dist
        // The algorithm will naturally converge for any finite alignment
        // Memory allocation failure would be the only practical limit (handled elsewhere)
        // For safety against pathological cases, we allow a very high limit that should
        // never be reached for any reasonable biological sequence comparison
        if max_dist > 100_000_000 {
            // Return best result found so far (this should never happen in practice)
            return result;
        }
    }
}

/// Wrapper for backward compatibility - uses slice lengths and forward direction
pub fn greedy_align_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize) {
    greedy_align_one_direction_ex(
        q_seq,
        s_seq,
        q_seq.len(),
        s_seq.len(),
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
        false,
    )
}

/// Internal greedy alignment function with explicit max_dist parameter.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// Returns: ((q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters), converged)
pub fn greedy_align_one_direction_with_max_dist(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    x_drop: i32,
    max_dist: usize,
    reverse: bool,
) -> ((usize, usize, i32, usize, usize, usize, usize), bool) {
    if len1 == 0 || len2 == 0 {
        return ((0, 0, 0, 0, 0, 0, 0), true);
    }

    // Calculate match and mismatch costs for distance-based tracking
    // NCBI reference: greedy_align.c:380-453
    // NCBI BLAST doubles scores if reward is odd to avoid fractions
    // IMPORTANT: x_drop must also be doubled to maintain consistent units
    // NCBI: if (match_score % 2 == 1) { match_score *= 2; mismatch_score *= 2; xdrop_threshold *= 2; }
    let (match_cost, mismatch_cost, scaled_xdrop) = if reward % 2 == 1 {
        (reward * 2, (-penalty) * 2, x_drop * 2)
    } else {
        (reward, -penalty, x_drop)
    };

    let op_cost = match_cost + mismatch_cost; // Cost of a mismatch in distance terms

    // X-drop offset for score comparison (using scaled_xdrop for consistent units)
    // NCBI reference: greedy_align.c:452-453
    // xdrop_offset = (xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1;
    let xdrop_offset = ((scaled_xdrop + match_cost / 2) / op_cost.max(1) + 1) as usize;

    // Find initial run of matches
    let initial_matches = find_first_mismatch_ex(q_seq, s_seq, len1, len2, 0, 0, reverse);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (
                initial_matches,
                initial_matches,
                score,
                initial_matches,
                0,
                0,
                0,
            ),
            true,
        );
    }

    // Diagonal origin (center of diagonal space)
    // NCBI BLAST uses max_dist (not scaled_max_dist) for diag_origin
    let diag_origin = max_dist + 2;
    let array_size = 2 * diag_origin + 4;
    let max_score_size = max_dist + xdrop_offset + 2;

    // Use thread-local memory pool to avoid per-call allocation overhead
    // This is critical for performance - without pooling, non-affine greedy is 20x slower
    GREEDY_MEM.with(|mem_cell| {
        let mut mem = mem_cell.borrow_mut();
        mem.ensure_capacity(array_size, max_score_size);
        mem.reset(array_size, max_score_size);

        // Initialize distance 0
        mem.last_seq2_off_a[diag_origin] = initial_matches as i32;
        mem.max_score[xdrop_offset] = (initial_matches as i32) * match_cost;

        let mut best_dist = 0usize;
        let mut best_seq1_len = initial_matches;
        let mut best_seq2_len = initial_matches;

        let mut diag_lower = diag_origin as i32 - 1;
        let mut diag_upper = diag_origin as i32 + 1;
        let mut end1_reached = initial_matches == len1;
        let mut end2_reached = initial_matches == len2;

        // Track convergence
        let mut converged = false;

        // Use flag to track which array is "prev" (true = a is prev, false = b is prev)
        let mut use_a_as_prev = true;

        // For each distance (use max_dist for non-affine greedy)
        for d in 1..=max_dist {
            if diag_lower > diag_upper {
                converged = true;
                break; // Converged
            }

            // Set sentinel values on the "prev" array
            if use_a_as_prev {
                if diag_lower >= 1 {
                    mem.last_seq2_off_a[(diag_lower - 1) as usize] = -2;
                    mem.last_seq2_off_a[diag_lower as usize] = -2;
                }
                if (diag_upper as usize) < array_size - 1 {
                    mem.last_seq2_off_a[diag_upper as usize] = -2;
                    mem.last_seq2_off_a[(diag_upper + 1) as usize] = -2;
                }
            } else {
                if diag_lower >= 1 {
                    mem.last_seq2_off_b[(diag_lower - 1) as usize] = -2;
                    mem.last_seq2_off_b[diag_lower as usize] = -2;
                }
                if (diag_upper as usize) < array_size - 1 {
                    mem.last_seq2_off_b[diag_upper as usize] = -2;
                    mem.last_seq2_off_b[(diag_upper + 1) as usize] = -2;
                }
            }

            // X-drop score threshold (using scaled_xdrop for consistent units)
            // NCBI reference: greedy_align.c:496, 531-533
            // max_score = aux_data->max_score + xdrop_offset;  (pointer arithmetic)
            // xdrop_score = max_score[d - xdrop_offset] + (match_cost + mismatch_cost) * d - xdrop_threshold;
            // xdrop_score = (Int4)ceil((double)xdrop_score / (match_cost / 2));
            // Note: In NCBI, max_score[d - xdrop_offset] accesses aux_data->max_score[d] when d >= xdrop_offset
            // When d < xdrop_offset, d - xdrop_offset is negative, accessing before the array start
            // NCBI initializes max_score[0..xdrop_offset] to 0, so negative index access returns 0
            // LOSAT uses: mem.max_score[xdrop_idx + xdrop_offset] where xdrop_idx = max(0, d - xdrop_offset)
            // This is equivalent: when d < xdrop_offset, we use max_score[xdrop_offset] = 0 (initialized)
            let xdrop_idx = if d >= xdrop_offset {
                d - xdrop_offset
            } else {
                0 // When d < xdrop_offset, use index 0 which maps to max_score[xdrop_offset] = 0
            };
            let xdrop_score =
                mem.max_score[xdrop_idx + xdrop_offset] + (op_cost * d as i32) - scaled_xdrop;
            // NCBI uses ceil division: (Int4)ceil((double)xdrop_score / (match_cost / 2))
            // Rust equivalent: (xdrop_score + match_cost / 2 - 1) / (match_cost / 2)
            let xdrop_score = (xdrop_score + match_cost / 2 - 1) / (match_cost / 2);

            let mut curr_extent = 0i32;
            let mut curr_seq2_index = 0i32;
            let mut curr_diag = diag_origin as i32;
            let tmp_diag_lower = diag_lower;
            let tmp_diag_upper = diag_upper;

            // For each diagonal
            for k in tmp_diag_lower..=tmp_diag_upper {
                let ku = k as usize;
                if ku >= array_size || ku == 0 {
                    continue;
                }

                // Find largest seq2 offset that increases distance from d-1 to d
                // NCBI reference: greedy_align.c:548-550
                // seq2_index = MAX(last_seq2_off[d - 1][k + 1], last_seq2_off[d - 1][k]) + 1;
                // seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
                // Access the "prev" array based on the flag
                let (prev_k_plus, prev_k, prev_k_minus) = if use_a_as_prev {
                    let p_plus = if ku + 1 < array_size {
                        mem.last_seq2_off_a[ku + 1]
                    } else {
                        -2
                    };
                    let p = mem.last_seq2_off_a[ku];
                    let p_minus = if ku >= 1 {
                        mem.last_seq2_off_a[ku - 1]
                    } else {
                        -2
                    };
                    (p_plus, p, p_minus)
                } else {
                    let p_plus = if ku + 1 < array_size {
                        mem.last_seq2_off_b[ku + 1]
                    } else {
                        -2
                    };
                    let p = mem.last_seq2_off_b[ku];
                    let p_minus = if ku >= 1 {
                        mem.last_seq2_off_b[ku - 1]
                    } else {
                        -2
                    };
                    (p_plus, p, p_minus)
                };

                // NCBI: seq2_index = MAX(prev_k_plus, prev_k) + 1; then MAX(seq2_index, prev_k_minus)
                let mut seq2_index = prev_k_plus.max(prev_k) + 1;
                seq2_index = seq2_index.max(prev_k_minus);

                let seq1_index = seq2_index + k - diag_origin as i32;

                if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                    // X-drop test failed or invalid diagonal
                    if k == diag_lower {
                        diag_lower += 1;
                    } else {
                        // Write to "curr" array
                        if use_a_as_prev {
                            mem.last_seq2_off_b[ku] = -2;
                        } else {
                            mem.last_seq2_off_a[ku] = -2;
                        }
                    }
                    continue;
                }

                diag_upper = k;

                // Slide down diagonal until mismatch
                let seq1_idx = seq1_index as usize;
                let seq2_idx = seq2_index as usize;

                if seq1_idx < len1 && seq2_idx < len2 {
                    let matches = find_first_mismatch_ex(
                        q_seq, s_seq, len1, len2, seq1_idx, seq2_idx, reverse,
                    );
                    let new_seq1_index = seq1_index + matches as i32;
                    let new_seq2_index = seq2_index + matches as i32;

                    // Write to "curr" array
                    if use_a_as_prev {
                        mem.last_seq2_off_b[ku] = new_seq2_index;
                    } else {
                        mem.last_seq2_off_a[ku] = new_seq2_index;
                    }

                    // Track best extent
                    let extent = new_seq1_index + new_seq2_index;
                    if extent > curr_extent {
                        curr_extent = extent;
                        curr_seq2_index = new_seq2_index;
                        curr_diag = k;
                    }

                    // Clamp bounds
                    // NCBI reference: greedy_align.c:607-614
                    // if (seq2_index == len2) {
                    //     diag_lower = k + 1;
                    //     end2_reached = TRUE;
                    // }
                    // if (seq1_index == len1) {
                    //     diag_upper = k - 1;
                    //     end1_reached = TRUE;
                    // }
                    // Note: new_seq2_index = seq2_index + matches, so we check == len2 (exact match)
                    if new_seq2_index as usize == len2 {
                        diag_lower = k + 1;
                        end2_reached = true;
                    }
                    if new_seq1_index as usize == len1 {
                        diag_upper = k - 1;
                        end1_reached = true;
                    }
                } else {
                    // Write to "curr" array
                    if use_a_as_prev {
                        mem.last_seq2_off_b[ku] = seq2_index;
                    } else {
                        mem.last_seq2_off_a[ku] = seq2_index;
                    }
                }
            }

            // Compute max score for this distance
            // NCBI reference: greedy_align.c:619-634
            // curr_score = curr_extent * (match_cost / 2) - d * (match_cost + mismatch_cost);
            // if (curr_score >= max_score[d - 1]) {
            //     max_score[d] = curr_score;
            //     best_dist = d;
            //     best_diag = curr_diag;
            //     *seq2_align_len = curr_seq2_index;
            //     *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
            // } else {
            //     max_score[d] = max_score[d - 1];
            // }
            let curr_score = (curr_extent * match_cost) / 2 - (d as i32) * op_cost;

            if curr_score >= mem.max_score[d - 1 + xdrop_offset] {
                mem.max_score[d + xdrop_offset] = curr_score;
                best_dist = d;
                best_seq2_len = curr_seq2_index as usize;
                best_seq1_len = (curr_seq2_index + curr_diag - diag_origin as i32) as usize;
            } else {
                mem.max_score[d + xdrop_offset] = mem.max_score[d - 1 + xdrop_offset];
            }

            // Check convergence
            if diag_lower > diag_upper {
                converged = true;
                break;
            }

            // Expand bounds for next distance
            if !end2_reached {
                diag_lower -= 1;
            }
            if !end1_reached {
                diag_upper += 1;
            }

            // Swap arrays by toggling the flag
            use_a_as_prev = !use_a_as_prev;
        }

        // Calculate final statistics
        // For greedy alignment, distance = mismatches + gaps
        // gap_letters is the difference in consumed lengths (indicates indels)
        let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
        // If there are gap letters, there's at least one gap open
        // Without traceback, we estimate gap_opens = 1 if any gaps exist
        let gap_opens = if gap_letters > 0 { 1 } else { 0 };
        // Mismatches = distance - gap_letters (distance includes both mismatches and gaps)
        let mismatches = best_dist.saturating_sub(gap_letters);
        let matches = best_seq1_len.min(best_seq2_len).saturating_sub(mismatches);

        // Calculate score
        let score = (matches as i32) * reward + (mismatches as i32) * penalty;

        (
            (
                best_seq1_len,
                best_seq2_len,
                score,
                matches,
                mismatches,
                gap_opens,
                gap_letters,
            ),
            converged,
        )
    })
}

/// Affine greedy alignment function (NCBI BLAST's BLAST_AffineGreedyAlign).
/// This handles affine gap penalties properly by tracking three separate offsets
/// for each diagonal: paths ending in insertion, match, or deletion.
///
/// Returns: ((q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters), converged)
pub fn affine_greedy_align_one_direction_with_max_dist(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    in_gap_open: i32,
    in_gap_extend: i32,
    x_drop: i32,
    max_dist: usize,
) -> ((usize, usize, i32, usize, usize, usize, usize), bool) {
    let len1 = q_seq.len();
    let len2 = s_seq.len();

    if len1 == 0 || len2 == 0 {
        return ((0, 0, 0, 0, 0, 0, 0), true);
    }

    // Score normalization (NCBI BLAST approach):
    // NCBI reference: greedy_align.c:799-808
    // if (match_score % 2 == 1) {
    //     match_score *= 2;
    //     mismatch_score *= 2;
    //     xdrop_threshold *= 2;
    //     in_gap_open *= 2;
    //     in_gap_extend *= 2;
    // }
    // Make sure bits of match_score don't disappear if divided by 2
    // IMPORTANT: In NCBI BLAST, mismatch_score is passed as -score_params->penalty,
    // where score_params->penalty is stored as a NEGATIVE value (e.g., -3).
    // So -(-3) = 3, meaning mismatch_score is the POSITIVE penalty magnitude.
    // In LOSAT, penalty is stored as a POSITIVE value (e.g., 3), so we use it directly.
    let (match_score, mismatch_penalty, xdrop_threshold, gap_open_in, gap_extend_in) =
        if reward % 2 == 1 {
            (
                reward * 2,
                penalty * 2,
                x_drop * 2,
                in_gap_open * 2,
                in_gap_extend * 2,
            )
        } else {
            (reward, penalty, x_drop, in_gap_open, in_gap_extend)
        };

    // Fill in derived scores and penalties (NCBI BLAST approach)
    // NCBI reference: greedy_align.c:835-843
    // match_score_half = match_score / 2;
    // op_cost = match_score + mismatch_score;
    // gap_open = in_gap_open;
    // gap_extend = in_gap_extend + match_score_half;
    // score_common_factor = BLAST_Gdb3(&op_cost, &gap_open, &gap_extend);
    // gap_open_extend = gap_open + gap_extend;
    // max_penalty = MAX(op_cost, gap_open_extend);
    // op_cost = match_score + mismatch_penalty (both positive, so op_cost is positive)
    // This represents the "cost" of a mismatch in distance units
    let match_score_half = match_score / 2;
    let mut op_cost = match_score + mismatch_penalty;
    let mut gap_open = gap_open_in;
    let mut gap_extend = gap_extend_in + match_score_half;
    let score_common_factor = gdb3(&mut op_cost, &mut gap_open, &mut gap_extend);
    let gap_open_extend = gap_open + gap_extend;
    let max_penalty = op_cost.max(gap_open_extend);

    // Scaled max_dist for affine alignment
    // With the correct sign convention (op_cost positive), gap_extend should always be positive
    let scaled_max_dist = (max_dist as i32) * gap_extend;

    // Diagonal origin (center of diagonal space)
    let diag_origin = max_dist + 2;
    let array_size = 2 * diag_origin + 4;

    // For affine greedy, we need to track diag_lower and diag_upper for ALL distances
    // (not just current), because contributions can come from distances < d-1
    // With correct sign convention, max_penalty should always be positive

    // Debug assertions to catch sign issues
    debug_assert!(op_cost > 0, "op_cost must be positive, got {}", op_cost);
    debug_assert!(
        gap_extend > 0,
        "gap_extend must be positive, got {}",
        gap_extend
    );
    debug_assert!(
        max_penalty > 0,
        "max_penalty must be positive, got {}",
        max_penalty
    );
    debug_assert!(
        scaled_max_dist >= 0,
        "scaled_max_dist must be non-negative, got {}",
        scaled_max_dist
    );

    // Safe conversion from i32 to usize with bounds checking
    if max_penalty < 0 || scaled_max_dist < 0 {
        // Sign error - return early with non-convergence
        return ((0, 0, 0, 0, 0, 0, 0), false);
    }

    let max_penalty_usize = max_penalty as usize;
    let scaled_max_dist_usize = scaled_max_dist as usize;
    let bounds_size = scaled_max_dist_usize + 1 + max_penalty_usize + 1;

    let mut diag_lower_arr: Vec<i32> = vec![INVALID_DIAG; bounds_size];
    let mut diag_upper_arr: Vec<i32> = vec![-INVALID_DIAG; bounds_size];

    // Initialize negative distance bounds with empty ranges
    // (for distances < max_penalty, which map to indices 0..max_penalty_usize)
    for i in 0..max_penalty_usize {
        diag_lower_arr[i] = INVALID_DIAG;
        diag_upper_arr[i] = -INVALID_DIAG;
    }

    // last_seq2_off[d][k] stores GreedyOffset for diagonal k at distance d
    // We need to keep all rows for affine alignment (contributions from d - max_penalty)
    // Allocate enough rows: scaled_max_dist + max_penalty + 2
    let num_rows = scaled_max_dist_usize + max_penalty_usize + 2;

    // Check allocation sizes to prevent memory exhaustion
    let total_elements = num_rows.checked_mul(array_size);
    if total_elements.is_none() || total_elements.unwrap() > 100_000_000 {
        // Allocation would be too large, return early with non-convergence
        return ((0, 0, 0, 0, 0, 0, 0), false);
    }

    let mut last_seq2_off: Vec<Vec<GreedyOffset>> = vec![
        vec![
            GreedyOffset {
                insert_off: INVALID_OFFSET,
                match_off: INVALID_OFFSET,
                delete_off: INVALID_OFFSET,
            };
            array_size
        ];
        num_rows
    ];

    // Max score at each distance for X-drop
    // NCBI reference: greedy_align.c:872-873
    // xdrop_offset = (xdrop_threshold + match_score_half) / score_common_factor + 1;
    let xdrop_offset =
        ((xdrop_threshold + match_score_half) / score_common_factor.max(1) + 1) as usize;
    let max_score_size = scaled_max_dist_usize + xdrop_offset + 2;
    let mut max_score: Vec<i32> = vec![0; max_score_size];

    // Find initial run of matches
    let initial_matches = find_first_mismatch(q_seq, s_seq, 0, 0);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (
                initial_matches,
                initial_matches,
                score,
                initial_matches,
                0,
                0,
                0,
            ),
            true,
        );
    }

    // Initialize distance 0
    last_seq2_off[0][diag_origin].match_off = initial_matches as i32;
    last_seq2_off[0][diag_origin].insert_off = INVALID_OFFSET;
    last_seq2_off[0][diag_origin].delete_off = INVALID_OFFSET;
    max_score[xdrop_offset] = (initial_matches as i32) * match_score;
    diag_lower_arr[max_penalty_usize] = diag_origin as i32;
    diag_upper_arr[max_penalty_usize] = diag_origin as i32;

    let mut best_dist = 0i32;
    let mut best_diag = diag_origin as i32;
    let mut best_seq1_len = initial_matches;
    let mut best_seq2_len = initial_matches;
    let _ = best_diag; // Suppress unused warning for initial value

    // Set up for distance 1
    let mut curr_diag_lower = diag_origin as i32 - 1;
    let mut curr_diag_upper = diag_origin as i32 + 1;
    let mut end1_diag = 0i32;
    let mut end2_diag = 0i32;
    let mut num_nonempty_dist = 1i32;
    let mut d = 1i32;

    let mut converged = false;

    // Helper function to safely access diag bounds
    let get_diag_lower = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() {
            arr[idx]
        } else {
            INVALID_DIAG
        }
    };
    let get_diag_upper = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() {
            arr[idx]
        } else {
            -INVALID_DIAG
        }
    };

    // For each distance
    while d <= scaled_max_dist {
        // Compute X-dropoff score threshold
        // NCBI reference: greedy_align.c:920, 979-983
        // max_score = aux_data->max_score + xdrop_offset;  (pointer arithmetic)
        // xdrop_score = max_score[d - xdrop_offset] + score_common_factor * d - xdrop_threshold;
        // xdrop_score = (Int4)ceil((double)xdrop_score / match_score_half);
        // if (xdrop_score < 0) xdrop_score = 0;
        // Note: In NCBI, max_score[d - xdrop_offset] accesses aux_data->max_score[d] when d >= xdrop_offset
        // When d < xdrop_offset, d - xdrop_offset is negative, accessing before the array start
        // NCBI initializes max_score[0..xdrop_offset] to 0, so negative index access returns 0
        // LOSAT uses: max_score[xdrop_idx + xdrop_offset] where xdrop_idx = max(0, d - xdrop_offset)
        // This is equivalent: when d < xdrop_offset, we use max_score[xdrop_offset] = 0 (initialized)
        let xdrop_idx = if d as usize >= xdrop_offset {
            (d as usize) - xdrop_offset
        } else {
            0 // When d < xdrop_offset, use index 0 which maps to max_score[xdrop_offset] = 0
        };
        let xdrop_score_raw =
            max_score[xdrop_idx + xdrop_offset] + score_common_factor * d - xdrop_threshold;
        // NCBI: (Int4)ceil((double)xdrop_score / match_score_half)
        let xdrop_score = ((xdrop_score_raw as f64) / (match_score_half as f64)).ceil() as i32;
        // NCBI: if (xdrop_score < 0) xdrop_score = 0;
        let xdrop_score = xdrop_score.max(0);

        let mut curr_extent = 0i32;
        let mut curr_seq2_index = 0i32;
        let mut curr_diag = 0i32;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;

        // For each valid diagonal
        for k in tmp_diag_lower..=tmp_diag_upper {
            let ku = k as usize;
            if ku >= array_size || ku == 0 {
                continue;
            }

            // Find best offset for DELETE (gap in seq1) - look at k+1 diagonal
            // NCBI reference: greedy_align.c:1000-1021
            // seq2_index = kInvalidOffset;
            // if (k + 1 <= diag_upper[d - gap_open_extend] && k + 1 >= diag_lower[d - gap_open_extend]) {
            //     seq2_index = last_seq2_off[d - gap_open_extend][k+1].match_off;
            // }
            // if (k + 1 <= diag_upper[d - gap_extend] && k + 1 >= diag_lower[d - gap_extend] &&
            //     seq2_index < last_seq2_off[d - gap_extend][k+1].delete_off) {
            //     seq2_index = last_seq2_off[d - gap_extend][k+1].delete_off;
            // }
            // if (seq2_index == kInvalidOffset)
            //     last_seq2_off[d][k].delete_off = kInvalidOffset;
            // else
            //     last_seq2_off[d][k].delete_off = seq2_index + 1;
            let mut seq2_index_del = INVALID_OFFSET;

            // From gap opening (match -> delete)
            let d_open = d - gap_open_extend;
            if d_open >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_open, max_penalty_usize);
                let du = get_diag_upper(&diag_upper_arr, d_open, max_penalty_usize);
                // NCBI: k + 1 <= diag_upper && k + 1 >= diag_lower
                if k + 1 >= dl && k + 1 <= du && (d_open as usize) < last_seq2_off.len() {
                    let ku1 = (k + 1) as usize;
                    if ku1 < array_size {
                        seq2_index_del = last_seq2_off[d_open as usize][ku1].match_off;
                    }
                }
            }

            // From gap extension (delete -> delete)
            // NCBI: if (k + 1 <= diag_upper[d - gap_extend] && k + 1 >= diag_lower[d - gap_extend] &&
            //      seq2_index < last_seq2_off[d - gap_extend][k+1].delete_off)
            let d_ext = d - gap_extend;
            if d_ext >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_ext, max_penalty_usize);
                let du = get_diag_upper(&diag_upper_arr, d_ext, max_penalty_usize);
                if k + 1 >= dl && k + 1 <= du && (d_ext as usize) < last_seq2_off.len() {
                    let ku1 = (k + 1) as usize;
                    if ku1 < array_size {
                        let ext_off = last_seq2_off[d_ext as usize][ku1].delete_off;
                        // NCBI: seq2_index < last_seq2_off[d - gap_extend][k+1].delete_off
                        if ext_off > seq2_index_del {
                            seq2_index_del = ext_off;
                        }
                    }
                }
            }

            // Save delete offset (deletion means seq2 offset slips by one)
            // NCBI: if (seq2_index == kInvalidOffset) ... else ... seq2_index + 1
            let du = d as usize;
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].delete_off = if seq2_index_del == INVALID_OFFSET {
                    INVALID_OFFSET
                } else {
                    seq2_index_del + 1
                };
            }

            // Find best offset for INSERT (gap in seq2) - look at k-1 diagonal
            // NCBI reference: greedy_align.c:1026-1036
            // seq2_index = kInvalidOffset;
            // if (k - 1 <= diag_upper[d - gap_open_extend] && k - 1 >= diag_lower[d - gap_open_extend]) {
            //     seq2_index = last_seq2_off[d - gap_open_extend][k-1].match_off;
            // }
            // if (k - 1 <= diag_upper[d - gap_extend] && k - 1 >= diag_lower[d - gap_extend] &&
            //     seq2_index < last_seq2_off[d - gap_extend][k-1].insert_off) {
            //     seq2_index = last_seq2_off[d - gap_extend][k-1].insert_off;
            // }
            // last_seq2_off[d][k].insert_off = seq2_index;
            let mut seq2_index_ins = INVALID_OFFSET;

            // From gap opening (match -> insert)
            if d_open >= 0 && k >= 1 {
                let dl = get_diag_lower(&diag_lower_arr, d_open, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_open, max_penalty_usize);
                // NCBI: k - 1 <= diag_upper && k - 1 >= diag_lower
                if k - 1 >= dl && k - 1 <= du_bound && (d_open as usize) < last_seq2_off.len() {
                    let km1 = (k - 1) as usize;
                    if km1 < array_size {
                        seq2_index_ins = last_seq2_off[d_open as usize][km1].match_off;
                    }
                }
            }

            // From gap extension (insert -> insert)
            // NCBI: if (k - 1 <= diag_upper[d - gap_extend] && k - 1 >= diag_lower[d - gap_extend] &&
            //      seq2_index < last_seq2_off[d - gap_extend][k-1].insert_off)
            if d_ext >= 0 && k >= 1 {
                let dl = get_diag_lower(&diag_lower_arr, d_ext, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_ext, max_penalty_usize);
                if k - 1 >= dl && k - 1 <= du_bound && (d_ext as usize) < last_seq2_off.len() {
                    let km1 = (k - 1) as usize;
                    if km1 < array_size {
                        let ext_off = last_seq2_off[d_ext as usize][km1].insert_off;
                        // NCBI: seq2_index < last_seq2_off[d - gap_extend][k-1].insert_off
                        if ext_off > seq2_index_ins {
                            seq2_index_ins = ext_off;
                        }
                    }
                }
            }

            // Save insert offset (insertion doesn't change seq2 offset)
            // NCBI: last_seq2_off[d][k].insert_off = seq2_index;
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].insert_off = seq2_index_ins;
            }

            // Compare with mismatch path (from diagonal k at d - op_cost)
            // NCBI reference: greedy_align.c:1041-1047
            // seq2_index = MAX(last_seq2_off[d][k].insert_off, last_seq2_off[d][k].delete_off);
            // if (k <= diag_upper[d - op_cost] && k >= diag_lower[d - op_cost]) {
            //     seq2_index = MAX(seq2_index, last_seq2_off[d - op_cost][k].match_off + 1);
            // }
            let mut seq2_index = last_seq2_off[du][ku]
                .insert_off
                .max(last_seq2_off[du][ku].delete_off);

            let d_mismatch = d - op_cost;
            if d_mismatch >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_mismatch, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_mismatch, max_penalty_usize);
                // NCBI: k <= diag_upper[d - op_cost] && k >= diag_lower[d - op_cost]
                if k >= dl && k <= du_bound && (d_mismatch as usize) < last_seq2_off.len() {
                    let match_off = last_seq2_off[d_mismatch as usize][ku].match_off;
                    if match_off != INVALID_OFFSET {
                        // NCBI: MAX(seq2_index, last_seq2_off[d - op_cost][k].match_off + 1)
                        seq2_index = seq2_index.max(match_off + 1);
                    }
                }
            }

            // Choose seq1 offset to remain on diagonal k
            let seq1_index = seq2_index + k - diag_origin as i32;

            // X-dropoff test
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else if du < last_seq2_off.len() {
                    last_seq2_off[du][ku].match_off = INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;

            // Slide down diagonal until mismatch
            let seq1_idx = seq1_index as usize;
            let seq2_idx = seq2_index as usize;

            let matches = if seq1_idx < len1 && seq2_idx < len2 {
                find_first_mismatch(q_seq, s_seq, seq1_idx, seq2_idx)
            } else {
                0
            };

            let new_seq1_index = seq1_index + matches as i32;
            let new_seq2_index = seq2_index + matches as i32;

            // Save match offset
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].match_off = new_seq2_index;
            }

            // Track best extent
            let extent = new_seq1_index + new_seq2_index;
            if extent > curr_extent {
                curr_extent = extent;
                curr_seq2_index = new_seq2_index;
                curr_diag = k;
            }

            // Clamp bounds to avoid walking off sequences
            // NCBI reference: greedy_align.c:1103-1110
            // if (seq1_index == len1) {
            //     curr_diag_upper = k;
            //     end1_diag = k - 1;
            // }
            // if (seq2_index == len2) {
            //     curr_diag_lower = k;
            //     end2_diag = k + 1;
            // }
            // Note: new_seq1_index = seq1_index + matches, so we check == len1 (exact match)
            if new_seq1_index as usize == len1 {
                curr_diag_upper = k;
                end1_diag = k - 1;
            }
            if new_seq2_index as usize == len2 {
                curr_diag_lower = k;
                end2_diag = k + 1;
            }
        }

        // Compute maximum score for distance d
        // NCBI reference: greedy_align.c:1115-1129
        // curr_score = curr_extent * match_score_half - d * score_common_factor;
        // if (curr_score > max_score[d - 1]) {
        //     max_score[d] = curr_score;
        //     best_dist = d;
        //     best_diag = curr_diag;
        //     *seq2_align_len = curr_seq2_index;
        //     *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        // } else {
        //     max_score[d] = max_score[d - 1];
        // }
        let curr_score = curr_extent * match_score_half - d * score_common_factor;

        // Update best if this is better
        // NCBI: if (curr_score > max_score[d - 1]) (strict greater than)
        let prev_max = if d >= 1 && (d - 1) as usize + xdrop_offset < max_score.len() {
            max_score[(d - 1) as usize + xdrop_offset]
        } else {
            0
        };

        if curr_score > prev_max {
            if (d as usize) + xdrop_offset < max_score.len() {
                max_score[(d as usize) + xdrop_offset] = curr_score;
            }
            best_dist = d;
            best_diag = curr_diag;
            best_seq2_len = curr_seq2_index as usize;
            best_seq1_len = (curr_seq2_index + best_diag - diag_origin as i32) as usize;
        } else if (d as usize) + xdrop_offset < max_score.len() {
            max_score[(d as usize) + xdrop_offset] = prev_max;
        }

        // Save diagonal bounds for this distance
        // NCBI reference: greedy_align.c:1139-1147
        // if (curr_diag_lower <= curr_diag_upper) {
        //     num_nonempty_dist++;
        //     diag_lower[d] = curr_diag_lower;
        //     diag_upper[d] = curr_diag_upper;
        // } else {
        //     diag_lower[d] = kInvalidDiag;
        //     diag_upper[d] = -kInvalidDiag;
        // }
        // Note: d maps to index (d + max_penalty) in the array
        let bounds_idx = (d + max_penalty) as usize;
        if bounds_idx < diag_lower_arr.len() {
            if curr_diag_lower <= curr_diag_upper {
                num_nonempty_dist += 1;
                diag_lower_arr[bounds_idx] = curr_diag_lower;
                diag_upper_arr[bounds_idx] = curr_diag_upper;
            } else {
                diag_lower_arr[bounds_idx] = INVALID_DIAG;
                diag_upper_arr[bounds_idx] = -INVALID_DIAG;
            }
        }

        // Check if we should decrement num_nonempty_dist
        // NCBI reference: greedy_align.c:1149-1150
        // if (diag_lower[d - max_penalty] <= diag_upper[d - max_penalty])
        //     num_nonempty_dist--;
        // Note: d - max_penalty maps to index (d - max_penalty + max_penalty) = d in the array
        let old_bounds_idx = (d - max_penalty + max_penalty) as usize;
        if old_bounds_idx < diag_lower_arr.len() && d >= max_penalty {
            let old_lower = diag_lower_arr[old_bounds_idx];
            let old_upper = diag_upper_arr[old_bounds_idx];
            // NCBI: if (diag_lower[d - max_penalty] <= diag_upper[d - max_penalty])
            if old_lower <= old_upper {
                num_nonempty_dist -= 1;
            }
        }

        // Convergence check: max_penalty consecutive empty ranges
        if num_nonempty_dist == 0 {
            converged = true;
            break;
        }

        // Compute diagonal range for next distance
        d += 1;

        let d_goe = d - gap_open_extend;
        let d_ge = d - gap_extend;
        let d_op = d - op_cost;

        let lower_goe = get_diag_lower(&diag_lower_arr, d_goe, max_penalty_usize);
        let lower_ge = get_diag_lower(&diag_lower_arr, d_ge, max_penalty_usize);
        let lower_op = get_diag_lower(&diag_lower_arr, d_op, max_penalty_usize);

        curr_diag_lower = lower_goe.min(lower_ge) - 1;
        curr_diag_lower = curr_diag_lower.min(lower_op);

        if end2_diag > 0 {
            curr_diag_lower = curr_diag_lower.max(end2_diag);
        }

        let upper_goe = get_diag_upper(&diag_upper_arr, d_goe, max_penalty_usize);
        let upper_ge = get_diag_upper(&diag_upper_arr, d_ge, max_penalty_usize);
        let upper_op = get_diag_upper(&diag_upper_arr, d_op, max_penalty_usize);

        curr_diag_upper = upper_goe.max(upper_ge) + 1;
        curr_diag_upper = curr_diag_upper.max(upper_op);

        if end1_diag > 0 {
            curr_diag_upper = curr_diag_upper.min(end1_diag);
        }
    }

    if !converged {
        // Did not converge - return best result found so far
        // Calculate statistics from best alignment
        let alignment_len = best_seq1_len.max(best_seq2_len);
        let gap_letters = best_seq1_len.abs_diff(best_seq2_len);

        // Estimate statistics (without full traceback)
        // For affine greedy, we estimate based on distance and gap penalties
        let estimated_gaps = if gap_letters > 0 { 1 } else { 0 };
        let estimated_mismatches = (best_dist as usize).saturating_sub(estimated_gaps);
        let matches = alignment_len
            .saturating_sub(estimated_mismatches)
            .saturating_sub(gap_letters);

        let score = (matches as i32) * reward + (estimated_mismatches as i32) * penalty
            - (estimated_gaps as i32) * in_gap_open
            - (gap_letters as i32) * in_gap_extend;

        return (
            (
                best_seq1_len,
                best_seq2_len,
                score,
                matches,
                estimated_mismatches,
                estimated_gaps,
                gap_letters,
            ),
            false,
        );
    }

    // Calculate final statistics
    // For affine greedy, we need to estimate statistics from the alignment
    // Since we don't have full traceback, we estimate based on the distance and gap penalties
    let alignment_len = best_seq1_len.max(best_seq2_len);
    let gap_letters = best_seq1_len.abs_diff(best_seq2_len);

    // Estimate gap opens and mismatches from distance
    // In affine greedy: distance = sum of (gap_open_extend for each gap) + (gap_extend for each additional gap letter) + (op_cost for each mismatch)
    // We approximate: if there are gaps, assume one gap open
    let estimated_gap_opens = if gap_letters > 0 { 1 } else { 0 };

    // Remaining distance after accounting for gaps
    let gap_cost = if gap_letters > 0 {
        gap_open_extend + (gap_letters.saturating_sub(1) as i32) * gap_extend
    } else {
        0
    };
    let remaining_dist = (best_dist - gap_cost).max(0);
    let estimated_mismatches = if op_cost > 0 {
        (remaining_dist / op_cost) as usize
    } else {
        0
    };

    let matches = alignment_len
        .saturating_sub(estimated_mismatches)
        .saturating_sub(gap_letters);

    // Calculate final score
    let score = (matches as i32) * reward + (estimated_mismatches as i32) * penalty
        - (estimated_gap_opens as i32) * in_gap_open
        - (gap_letters as i32) * in_gap_extend;

    (
        (
            best_seq1_len,
            best_seq2_len,
            score,
            matches,
            estimated_mismatches,
            estimated_gap_opens,
            gap_letters,
        ),
        converged,
    )
}

// =============================================================================
// NCBI Greedy Gapped Alignment with Traceback (BLAST_GreedyGappedAlignment)
// =============================================================================

/// Convert preliminary greedy edit blocks to a gap edit script.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2471-2536
fn prelim_edit_block_to_gap_edit_script(
    rev_prelim_tback: &GapPrelimEditBlock,
    fwd_prelim_tback: &GapPrelimEditBlock,
) -> Option<GapEditScript> {
    let mut merge_ops = false;

    if rev_prelim_tback.num_ops > 0 && fwd_prelim_tback.num_ops > 0 {
        let rev_last = rev_prelim_tback.edit_ops[rev_prelim_tback.num_ops - 1].op_type;
        let fwd_last = fwd_prelim_tback.edit_ops[fwd_prelim_tback.num_ops - 1].op_type;
        if rev_last == fwd_last {
            merge_ops = true;
        }
    }

    let mut size = rev_prelim_tback.num_ops + fwd_prelim_tback.num_ops;
    if merge_ops {
        size = size.saturating_sub(1);
    }

    let mut esp = GapEditScript::new(size);
    let mut index = 0usize;

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2509-2515
    // ```c
    // index = 0;
    // for (i=0; i < rev_prelim_tback->num_ops; i++) {
    //    op = rev_prelim_tback->edit_ops + i;
    //    esp->op_type[index] = op->op_type;
    //    esp->num[index] = op->num;
    //    index++;
    // }
    // ```
    for op in rev_prelim_tback
        .edit_ops
        .iter()
        .take(rev_prelim_tback.num_ops)
    {
        esp.op_type[index] = op.op_type;
        esp.num[index] = op.num;
        index += 1;
    }

    if fwd_prelim_tback.num_ops == 0 {
        return Some(esp);
    }

    if merge_ops && index > 0 {
        esp.num[index - 1] += fwd_prelim_tback.edit_ops[fwd_prelim_tback.num_ops - 1].num;
    }

    let mut i: isize = if merge_ops {
        fwd_prelim_tback.num_ops as isize - 2
    } else {
        fwd_prelim_tback.num_ops as isize - 1
    };

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2523-2534
    // ```c
    // if (merge_ops)
    //    i = fwd_prelim_tback->num_ops - 2;
    // else
    //    i = fwd_prelim_tback->num_ops - 1;
    //
    // for (; i >= 0; i--) {
    //    op = fwd_prelim_tback->edit_ops + i;
    //    esp->op_type[index] = op->op_type;
    //    esp->num[index] = op->num;
    //    index++;
    // }
    // ```
    while i >= 0 {
        let op = fwd_prelim_tback.edit_ops[i as usize];
        esp.op_type[index] = op.op_type;
        esp.num[index] = op.num;
        index += 1;
        i -= 1;
    }

    Some(esp)
}

/// Update edit script around a substitution run.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2572-2632
fn update_edit_script(esp: &mut GapEditScript, pos: isize, bf: i32, af: i32) {
    if bf > 0 {
        let mut op = pos;
        let mut qd = bf;
        let mut sd = bf;
        loop {
            op -= 1;
            if op < 0 {
                return;
            }
            let idx = op as usize;
            match esp.op_type[idx] {
                GapAlignOpType::Sub => {
                    qd -= esp.num[idx];
                    sd -= esp.num[idx];
                }
                GapAlignOpType::Ins => {
                    qd -= esp.num[idx];
                }
                GapAlignOpType::Del => {
                    sd -= esp.num[idx];
                }
                GapAlignOpType::Invalid => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }

        let idx = op as usize;
        esp.num[idx] = -qd.max(sd);
        esp.op_type[idx] = GapAlignOpType::Sub;

        let mut op_idx = op + 1;
        while op_idx < pos - 1 {
            esp.num[op_idx as usize] = 0;
            op_idx += 1;
        }

        let pos_idx = pos as usize;
        if pos_idx < esp.num.len() {
            esp.num[pos_idx] += bf;
        }

        qd -= sd;
        let before_idx = (pos - 1) as usize;
        if before_idx < esp.op_type.len() {
            esp.op_type[before_idx] = if qd > 0 {
                GapAlignOpType::Del
            } else {
                GapAlignOpType::Ins
            };
            esp.num[before_idx] = if qd > 0 { qd } else { -qd };
        }
    }

    if af > 0 {
        let mut op = pos;
        let mut qd = af;
        let mut sd = af;
        loop {
            op += 1;
            if op >= esp.size as isize {
                return;
            }
            let idx = op as usize;
            match esp.op_type[idx] {
                GapAlignOpType::Sub => {
                    qd -= esp.num[idx];
                    sd -= esp.num[idx];
                }
                GapAlignOpType::Ins => {
                    qd -= esp.num[idx];
                }
                GapAlignOpType::Del => {
                    sd -= esp.num[idx];
                }
                GapAlignOpType::Invalid => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }

        let idx = op as usize;
        esp.num[idx] = -qd.max(sd);
        esp.op_type[idx] = GapAlignOpType::Sub;

        let mut op_idx = op - 1;
        while op_idx > pos + 1 {
            esp.num[op_idx as usize] = 0;
            op_idx -= 1;
        }

        let pos_idx = pos as usize;
        if pos_idx < esp.num.len() {
            esp.num[pos_idx] += af;
        }

        qd -= sd;
        let after_idx = (pos + 1) as usize;
        if after_idx < esp.op_type.len() {
            esp.op_type[after_idx] = if qd > 0 {
                GapAlignOpType::Del
            } else {
                GapAlignOpType::Ins
            };
            esp.num[after_idx] = if qd > 0 { qd } else { -qd };
        }
    }
}

/// Rebuild edit script after updates.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2634-2667
fn rebuild_edit_script(esp: &mut GapEditScript) {
    let mut j: isize = -1;
    for i in 0..esp.size {
        if esp.num[i] == 0 {
            continue;
        }
        if j >= 0 && esp.op_type[i] == esp.op_type[j as usize] {
            esp.num[j as usize] += esp.num[i];
        } else if j == -1
            || esp.op_type[i] == GapAlignOpType::Sub
            || esp.op_type[j as usize] == GapAlignOpType::Sub
        {
            j += 1;
            let j_idx = j as usize;
            esp.op_type[j_idx] = esp.op_type[i];
            esp.num[j_idx] = esp.num[i];
        } else {
            let j_idx = j as usize;
            let d = esp.num[j_idx] - esp.num[i];
            if d > 0 {
                esp.num[j_idx - 1] += esp.num[i];
                esp.num[j_idx] = d;
            } else if d < 0 {
                if j == 0 && i as isize - j > 0 {
                    esp.op_type[j_idx] = GapAlignOpType::Sub;
                    j += 1;
                } else {
                    esp.num[j_idx - 1] += esp.num[j_idx];
                }
                let j_idx = j as usize;
                esp.num[j_idx] = -d;
                esp.op_type[j_idx] = esp.op_type[i];
            } else {
                esp.num[j_idx - 1] += esp.num[j_idx];
                j -= 1;
            }
        }
    }
    esp.size = (j + 1) as usize;
}

/// Reduce small gaps in greedy edit script.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2669-2758
fn reduce_gaps(
    esp: &mut GapEditScript,
    q: &[u8],
    s: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
) {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2669-2758
    // ```c
    // for (q1=q, s1=s, i=0; i<esp->size; i++) {
    //     if (esp->num[i] == 0) continue;
    //     if (esp->op_type[i] == eGapAlignSub) {
    //         if(esp->num[i] >= 12) {
    //             nm1 = 1;
    //             if (i > 0) {
    //                 while (q1-nm1>=q && (*(q1-nm1) == *(s1-nm1))) ++nm1;
    //             }
    //             q1 += esp->num[i];
    //             s1 += esp->num[i];
    //             nm2 = 0;
    //             if (i < esp->size -1) {
    //                 while ((q1+1<qf) && (s1+1<sf) && (*(q1++) == *(s1++))) ++nm2;
    //             }
    //             if (nm1>1 || nm2>0) s_UpdateEditScript(esp, i, nm1-1, nm2);
    //             q1--; s1--;
    //         } else {
    //             q1 += esp->num[i];
    //             s1 += esp->num[i];
    //         }
    //     } else if (esp->op_type[i] == eGapAlignIns) {
    //         q1 += esp->num[i];
    //     } else {
    //         s1 += esp->num[i];
    //     }
    // }
    // ```
    let mut q_idx: isize = 0;
    let mut s_idx: isize = 0;
    let q_start = q_start as isize;
    let s_start = s_start as isize;
    let q_len = q_end as isize - q_start;
    let s_len = s_end as isize - s_start;

    for i in 0..esp.size {
        if esp.num[i] == 0 {
            continue;
        }

        match esp.op_type[i] {
            GapAlignOpType::Sub => {
                if esp.num[i] >= 12 {
                    let mut nm1: i32 = 1;
                    if i > 0 {
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2679-2681
                        // ```c
                        //                if (i > 0) {
                        //                    while (q1-nm1>=q && (*(q1-nm1) == *(s1-nm1))) ++nm1;
                        //                }
                        // ```
                        // NCBI bounds only the query; the subject is read before the
                        // alignment start, down to its leading sentinel, which matches no
                        // query letter, so the scan stops at the subject's first letter.
                        while (q_idx - nm1 as isize) >= 0
                            && (s_start + s_idx - nm1 as isize) >= 0
                            && q[(q_start + q_idx - nm1 as isize) as usize]
                                == s[(s_start + s_idx - nm1 as isize) as usize]
                        {
                            nm1 += 1;
                        }
                    }

                    q_idx += esp.num[i] as isize;
                    s_idx += esp.num[i] as isize;

                    let mut nm2: i32 = 0;
                    if i < esp.size - 1 {
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2687-2689
                        // ```c
                        // while ((q1+1<qf) && (s1+1<sf) && (*(q1++) == *(s1++))) ++nm2;
                        // ...
                        // q1--; s1--;
                        // ```
                        // The probe advances `q1`/`s1` even on the first failing
                        // comparison, then backs up by one after the loop.
                        while q_idx + 1 < q_len && s_idx + 1 < s_len {
                            let q_match =
                                q[(q_start + q_idx) as usize] == s[(s_start + s_idx) as usize];
                            q_idx += 1;
                            s_idx += 1;
                            if q_match {
                                nm2 += 1;
                            } else {
                                break;
                            }
                        }
                    }

                    if nm1 > 1 || nm2 > 0 {
                        update_edit_script(esp, i as isize, nm1 - 1, nm2);
                    }

                    q_idx -= 1;
                    s_idx -= 1;
                } else {
                    q_idx += esp.num[i] as isize;
                    s_idx += esp.num[i] as isize;
                }
            }
            GapAlignOpType::Ins => {
                q_idx += esp.num[i] as isize;
            }
            GapAlignOpType::Del => {
                s_idx += esp.num[i] as isize;
            }
            GapAlignOpType::Invalid => {}
        }
    }

    rebuild_edit_script(esp);

    let mut q_pos: isize = 0;
    let mut s_pos: isize = 0;
    let q_len = q_end as isize - q_start;
    let s_len = s_end as isize - s_start;

    for i in 0..esp.size {
        if esp.op_type[i] == GapAlignOpType::Sub {
            q_pos += esp.num[i] as isize;
            s_pos += esp.num[i] as isize;
            continue;
        }

        if i > 1 && esp.op_type[i] != esp.op_type[i - 2] && esp.num[i - 2] > 0 {
            let mut d = esp.num[i] + esp.num[i - 1] + esp.num[i - 2];
            if d == 3 {
                esp.num[i - 2] = 0;
                esp.num[i - 1] = 2;
                esp.num[i] = 0;
                if esp.op_type[i] == GapAlignOpType::Ins {
                    q_pos += 1;
                } else {
                    s_pos += 1;
                }
            } else if d < 12 {
                let mut nm1 = 0;
                let mut nm2 = 0;
                d = esp.num[i].min(esp.num[i - 2]);

                q_pos -= esp.num[i - 1] as isize;
                s_pos -= esp.num[i - 1] as isize;
                let mut q1 = q_pos;
                let mut s1 = s_pos;

                if esp.op_type[i] == GapAlignOpType::Ins {
                    s_pos -= d as isize;
                } else {
                    q_pos -= d as isize;
                }

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2729-2736
                // ```c
                // for (j=0; j<esp->num[i-1]; ++j, ++q1, ++s1, ++q, ++s) {
                //    if (*q1 == *s1) nm1++;
                //    if (*q == *s) nm2++;
                // }
                // ```
                for _ in 0..esp.num[i - 1] {
                    if q[(q_start + q1) as usize] == s[(s_start + s1) as usize] {
                        nm1 += 1;
                    }
                    if q[(q_start + q_pos) as usize] == s[(s_start + s_pos) as usize] {
                        nm2 += 1;
                    }
                    q1 += 1;
                    s1 += 1;
                    q_pos += 1;
                    s_pos += 1;
                }

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2737-2739
                // ```c
                // for (j=0; j<d; ++j, ++q, ++s) {
                //    if (*q == *s) nm2++;
                // }
                // ```
                for _ in 0..d {
                    if q[(q_start + q_pos) as usize] == s[(s_start + s_pos) as usize] {
                        nm2 += 1;
                    }
                    q_pos += 1;
                    s_pos += 1;
                }

                if nm2 >= nm1 - d {
                    esp.num[i - 2] -= d;
                    esp.num[i - 1] += d;
                    esp.num[i] -= d;
                } else {
                    q_pos = q1;
                    s_pos = s1;
                }
            }
        }

        if esp.op_type[i] == GapAlignOpType::Ins {
            q_pos += esp.num[i] as isize;
        } else {
            s_pos += esp.num[i] as isize;
        }
    }

    rebuild_edit_script(esp);
}

/// Unpack a base from NCBI2NA packed byte.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:365-369 (NCBI2NA_UNPACK_BASE usage)
#[inline]
fn ncbi2na_unpack_base(byte: u8, pos: u8) -> u8 {
    (byte >> (pos * 2)) & 0x03
}

/// Find first mismatch with fence detection.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:297-375 (s_FindFirstMismatch)
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:313-317
// ```c
// static NCBI_INLINE Int4 s_FindFirstMismatch(const Uint1 *seq1, const Uint1 *seq2,
//                                 Int4 len1, Int4 len2,
//                                 Int4 seq1_index, Int4 seq2_index,
//                                 Boolean *fence_hit,
//                                 Boolean reverse, Uint1 rem)
// ```
// Keep the short mismatch scans in the recurrence to avoid a call boundary
// for each surviving diagonal, as intended by NCBI_INLINE.
#[inline(always)]
fn find_first_mismatch_greedy(
    seq1: &[u8],
    seq2: &[u8],
    len1: i32,
    len2: i32,
    seq1_index: i32,
    seq2_index: i32,
    reverse: bool,
    rem: u8,
    fence_hit: &mut bool,
) -> i32 {
    let mut seq1_index = seq1_index;
    let mut seq2_index = seq2_index;
    let tmp = seq1_index;

    if reverse {
        if rem == 4 {
            while seq1_index < len1
                && seq2_index < len2
                && seq1[(len1 - 1 - seq1_index) as usize] < 4
                && seq1[(len1 - 1 - seq1_index) as usize] == seq2[(len2 - 1 - seq2_index) as usize]
            {
                seq1_index += 1;
                seq2_index += 1;
            }
            if seq2_index < len2 && seq2[(len2 - 1 - seq2_index) as usize] == FENCE_SENTRY {
                *fence_hit = true;
            }
        } else {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:338-344
            // ```c
            // while (seq1_index < len1 && seq2_index < len2 &&
            //     seq1[len1-1 - seq1_index] ==
            //     NCBI2NA_UNPACK_BASE(seq2[(len2-1-seq2_index) / 4],
            //                         3 - (len2-1-seq2_index) % 4)) {
            // ```
            // Reverse traversal over compressed subject does not add `rem`.
            while seq1_index < len1
                && seq2_index < len2
                && seq1[(len1 - 1 - seq1_index) as usize]
                    == ncbi2na_unpack_base(
                        seq2[((len2 - 1 - seq2_index) / 4) as usize],
                        (3 - (len2 - 1 - seq2_index) % 4) as u8,
                    )
            {
                seq1_index += 1;
                seq2_index += 1;
            }
        }
    } else {
        if rem == 4 {
            while seq1_index < len1
                && seq2_index < len2
                && seq1[seq1_index as usize] < 4
                && seq1[seq1_index as usize] == seq2[seq2_index as usize]
            {
                seq1_index += 1;
                seq2_index += 1;
            }
            if seq2_index < len2 && seq2[seq2_index as usize] == FENCE_SENTRY {
                *fence_hit = true;
            }
        } else {
            while seq1_index < len1
                && seq2_index < len2
                && seq1[seq1_index as usize]
                    == ncbi2na_unpack_base(
                        seq2[((seq2_index + rem as i32) / 4) as usize],
                        (3 - (seq2_index + rem as i32) % 4) as u8,
                    )
            {
                seq1_index += 1;
                seq2_index += 1;
            }
        }
    }

    seq1_index - tmp
}

/// Affine traceback helper from match state.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:148-178
fn get_next_affine_tback_from_match(
    last_seq2_off: &AffineRows<'_>,
    diag_lower: &[i32],
    diag_upper: &[i32],
    diag_offset: i32,
    d: &mut i32,
    diag: i32,
    op_cost: i32,
    seq2_index: &mut i32,
) -> GapAlignOpType {
    let mut new_seq2_index;

    let idx = (*d) - op_cost;
    let dl_idx = idx + diag_offset;
    if idx >= 0 && dl_idx >= 0 && (dl_idx as usize) < diag_lower.len() {
        if diag >= diag_lower[dl_idx as usize] && diag <= diag_upper[dl_idx as usize] {
            new_seq2_index = last_seq2_off.get(idx, diag).match_off;
            if new_seq2_index
                >= last_seq2_off
                    .get(*d, diag)
                    .insert_off
                    .max(last_seq2_off.get(*d, diag).delete_off)
            {
                *d -= op_cost;
                *seq2_index = new_seq2_index;
                return GapAlignOpType::Sub;
            }
        }
    }

    if last_seq2_off.get(*d, diag).insert_off > last_seq2_off.get(*d, diag).delete_off {
        *seq2_index = last_seq2_off.get(*d, diag).insert_off;
        GapAlignOpType::Ins
    } else {
        *seq2_index = last_seq2_off.get(*d, diag).delete_off;
        GapAlignOpType::Del
    }
}

/// Affine traceback helper from insertion/deletion state.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:199-258
fn get_next_affine_tback_from_indel(
    last_seq2_off: &AffineRows<'_>,
    diag_lower: &[i32],
    diag_upper: &[i32],
    diag_offset: i32,
    d: &mut i32,
    diag: i32,
    gap_open: i32,
    gap_extend: i32,
    state: GapAlignOpType,
) -> GapAlignOpType {
    let gap_open_extend = gap_open + gap_extend;
    let new_diag = if state == GapAlignOpType::Ins {
        diag - 1
    } else {
        diag + 1
    };

    let last_d = (*d) - gap_extend;
    let mut new_seq2_index = INVALID_OFFSET;
    let dl_idx = last_d + diag_offset;
    if last_d >= 0 && dl_idx >= 0 && (dl_idx as usize) < diag_lower.len() {
        if new_diag >= diag_lower[dl_idx as usize] && new_diag <= diag_upper[dl_idx as usize] {
            new_seq2_index = if state == GapAlignOpType::Ins {
                last_seq2_off.get(last_d, new_diag).insert_off
            } else {
                last_seq2_off.get(last_d, new_diag).delete_off
            };
        }
    }

    let last_d = (*d) - gap_open_extend;
    let dl_idx = last_d + diag_offset;
    if last_d >= 0 && dl_idx >= 0 && (dl_idx as usize) < diag_lower.len() {
        if new_diag >= diag_lower[dl_idx as usize]
            && new_diag <= diag_upper[dl_idx as usize]
            && new_seq2_index < last_seq2_off.get(last_d, new_diag).match_off
        {
            *d -= gap_open_extend;
            return GapAlignOpType::Sub;
        }
    }

    *d -= gap_extend;
    state
}

/// Non-affine traceback helper.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:260-294
fn get_next_non_affine_tback(
    last_seq2_off: &[NonAffineGreedyRow],
    non_affine_mem: &GreedyNonAffineMem,
    d: i32,
    diag: i32,
    seq2_index: &mut i32,
) -> i32 {
    let prev_row = &last_seq2_off[(d - 1) as usize];
    if prev_row.get(non_affine_mem, diag - 1)
        > prev_row
            .get(non_affine_mem, diag)
            .max(prev_row.get(non_affine_mem, diag + 1))
    {
        *seq2_index = prev_row.get(non_affine_mem, diag - 1);
        return diag - 1;
    }
    if prev_row.get(non_affine_mem, diag) > prev_row.get(non_affine_mem, diag + 1) {
        *seq2_index = prev_row.get(non_affine_mem, diag);
        return diag;
    }
    *seq2_index = prev_row.get(non_affine_mem, diag + 1);
    diag + 1
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-614,351-362,375
// ```c
// for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
//     seq2_index = MAX(last_seq2_off[d - 1][k + 1],
//                      last_seq2_off[d - 1][k    ]) + 1;
//     seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
// ...
//     index = s_FindFirstMismatch(seq1, seq2, len1, len2,
// ...
//     last_seq2_off[d][k] = seq2_index;
// ...
//     if (seq1_index + seq2_index > curr_extent) {
// ...
// while (seq1_index < len1 && seq2_index < len2 &&
//        seq1[seq1_index] < 4 &&
//        seq1[seq1_index] == seq2[seq2_index]) {
// ```
// The code below is a faster form of the loop over diagonals `k` of `BLAST_GreedyAlign`
// and of its `s_FindFirstMismatch` (`rem == 4`). Pass 1 reads only the row of distance
// d-1 (the `seq2_index` choice and the X-drop test). Pass 2 reads only the sequences
// (the slide length). Pass 3 applies the writes and the state updates (`diag_lower`,
// `diag_upper`, `curr_extent`, seed, sequence ends, fence) diagonal by diagonal in the
// order of the C loop. By construction this gives the cells, bounds and extent of the C
// loop; LOSAT_X_GREEDYSHADOW compares them row by row. All of it is integer arithmetic.
// ---------------------------------------------------------------------------
// EXPERIMENT (LOSAT_X_GREEDYFAST / LOSAT_X_GREEDYSHADOW): one distance of
// `BLAST_GreedyAlign` over an uncompressed subject, with the same cell order,
// the same writes and the same bookkeeping as the loop below, but as a small
// leaf routine that compares eight bases per step.
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:537-611
// ```c
// for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
//     seq2_index = MAX(last_seq2_off[d - 1][k + 1],
//                      last_seq2_off[d - 1][k    ]) + 1;
//     seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
//     seq1_index = seq2_index + k - diag_origin;
//     if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {
//         if (k == diag_lower) diag_lower++;
//         else last_seq2_off[d][k] = kInvalidOffset;
//         continue;
//     }
//     diag_upper = k;
//     index = s_FindFirstMismatch(seq1, seq2, len1, len2, seq1_index,
//                                 seq2_index, fence_hit, reverse, rem);
//     if (fence_hit && *fence_hit) return 0;
//     ...
// }
// ```
// ---------------------------------------------------------------------------

#[cfg(test)]
thread_local! {
    /// Lets a test choose the mode per call.
    static X_GREEDY_TEST_MODE: std::cell::Cell<Option<u8>> = const { std::cell::Cell::new(None) };
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537,571-573
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ...
///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
/// ```
/// Reads the LOSAT_X_GREEDYFAST / LOSAT_X_GREEDYSHADOW switches once. The mode only
/// chooses whether this loop of `BLAST_GreedyAlign` runs as ported (0), as the fast row
/// (1), or as both with a comparison (2).
/// 0 = reference loop, 1 = fast row, 2 = both and compare every row.
fn x_greedy_mode() -> u8 {
    #[cfg(test)]
    if let Some(mode) = X_GREEDY_TEST_MODE.with(|mode| mode.get()) {
        return mode;
    }
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_GREEDYSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_GREEDYFAST").is_some() {
            1
        } else {
            0
        }
    })
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:565,578-582,596-614
/// ```c
///             diag_upper = k;
/// ...
///             if (index > longest_match_run) {
///                 seed->start_q = seq1_index;
/// ...
///             if (seq1_index + seq2_index > curr_extent) {
///                 curr_extent = seq1_index + seq2_index;
///                 curr_seq2_index = seq2_index;
///                 curr_diag = k;
/// ...
///                 diag_lower = k + 1;
///                 end2_reached = TRUE;
/// ```
/// These are the variables the C loop updates across diagonals, copied in and out of the
/// fast row so that it can leave exactly the state the C loop would leave.
/// Everything a distance reads and updates besides the two rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct XGreedyRowState {
    diag_lower: i32,
    diag_upper: i32,
    end1_reached: bool,
    end2_reached: bool,
    curr_extent: i32,
    curr_seq2_index: i32,
    curr_diag: i32,
    longest_match_run: i32,
    seed_start_q: i32,
    seed_start_s: i32,
    fence_hit: bool,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:327-334,351-358,375
/// ```c
///     if (reverse) {
///         if (rem == 4) {
///             while (seq1_index < len1 && seq2_index < len2 &&
///                    seq1[len1-1 - seq1_index] < 4 &&
///                    seq1[len1-1 - seq1_index] == seq2[len2-1 - seq2_index]) {
/// ...
///             while (seq1_index < len1 && seq2_index < len2 &&
///                    seq1[seq1_index] < 4 &&
///                    seq1[seq1_index] == seq2[seq2_index]) {
/// ...
///     return seq1_index - tmp;
/// ```
/// Same result as the byte loop: the number of positions before the first one where
/// the bytes differ or `seq1` is not a base. Eight positions are tested per step.
/// `s_FindFirstMismatch` for `rem == 4`, without the fence test.
///
/// A position matches when the two bytes are equal and the `seq1` byte is
/// below 4, i.e. when `(b1 ^ b2) | (b1 & 0xFC)` is zero; eight positions are
/// tested per step and the first non-zero byte is the first mismatch.
///
/// SAFETY: `seq1` and `seq2` must be readable for `len1` and `len2` bytes.
#[inline(always)]
unsafe fn x_first_mismatch<const REVERSE: bool>(
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    seq1_index: i32,
    seq2_index: i32,
) -> i32 {
    const NOT_A_BASE: u64 = 0xFCFC_FCFC_FCFC_FCFC;
    // Number of positions both sequences still have (may be <= 0).
    let n = (len1 - seq1_index).min(len2 - seq2_index);
    let mut i: i32 = 0;
    if !REVERSE {
        while n - i >= 8 {
            let w1 =
                u64::from_le((seq1.add((seq1_index + i) as usize) as *const u64).read_unaligned());
            let w2 =
                u64::from_le((seq2.add((seq2_index + i) as usize) as *const u64).read_unaligned());
            let diff = (w1 ^ w2) | (w1 & NOT_A_BASE);
            if diff != 0 {
                return i + (diff.trailing_zeros() >> 3) as i32;
            }
            i += 8;
        }
        while i < n {
            let b1 = *seq1.add((seq1_index + i) as usize);
            if b1 >= 4 || b1 != *seq2.add((seq2_index + i) as usize) {
                break;
            }
            i += 1;
        }
    } else {
        // Position `i` reads seq1[len1 - 1 - seq1_index - i]: the eight
        // positions i..i+8 are the eight bytes ending there, last byte first.
        while n - i >= 8 {
            let w1 = u64::from_le(
                (seq1.add((len1 - seq1_index - i - 8) as usize) as *const u64).read_unaligned(),
            );
            let w2 = u64::from_le(
                (seq2.add((len2 - seq2_index - i - 8) as usize) as *const u64).read_unaligned(),
            );
            let diff = (w1 ^ w2) | (w1 & NOT_A_BASE);
            if diff != 0 {
                return i + (diff.leading_zeros() >> 3) as i32;
            }
            i += 8;
        }
        while i < n {
            let b1 = *seq1.add((len1 - 1 - seq1_index - i) as usize);
            if b1 >= 4 || b1 != *seq2.add((len2 - 1 - seq2_index - i) as usize) {
                break;
            }
            i += 1;
        }
    }
    i
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:578-582,596-614
/// ```c
///             if (index > longest_match_run) {
///                 seed->start_q = seq1_index;
///                 seed->start_s = seq2_index;
///                 seed->match_length = longest_match_run = index;
/// ...
///             if (seq1_index == len1) {
///                 diag_upper = k - 1;
///                 end1_reached = TRUE;
/// ```
/// The variables that the C loop changes only on some diagonals (seed, ends, fence). The
/// fast row keeps them here and remembers which cell last changed them.
/// What a distance updates only on a few of its diagonals. The row loop
/// keeps it in memory and leaves those diagonals to `x_greedy_cell_slow`.
struct XGreedyRowCold {
    diag_lower: i32,
    end1_reached: bool,
    end2_reached: bool,
    curr_extent: i32,
    /// Cell that set `curr_extent` in this distance.
    best_cell: usize,
    /// Last cell whose slide ended at `len1`.
    end1_cell: usize,
    longest_match_run: i32,
    seed_start_q: i32,
    seed_start_s: i32,
    fence_hit: bool,
    /// The reference loop returns from inside this distance.
    stopped: bool,
}

const X_NO_CELL: usize = usize::MAX;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:571-614
/// ```c
///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
///                                         seq1_index, seq2_index,
///                                         fence_hit, reverse, rem);
///             if(fence_hit && *fence_hit){
///                 return 0;
///             }
/// ...
///             last_seq2_off[d][k] = seq2_index;
/// ...
///             if (seq2_index == len2) {
/// ```
/// The C code after the X-drop test, for one diagonal, with the same order of tests and
/// writes (slide, fence, seed, cell, extent, ends).
/// One surviving diagonal, exactly as the reference loop handles it after
/// its X-drop test: any slide length, the fence test, the seed, the extent
/// and both sequence ends.
///
/// SAFETY: `seq1` and `seq2` must be readable for `len1` and `len2` bytes and
/// `current` must be writable at `cell`.
#[inline(never)]
unsafe fn x_greedy_cell_slow<const REVERSE: bool>(
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    seq1_index: i32,
    seq2_index: i32,
    cell: usize,
    k: i32,
    current: *mut i32,
    cold: &mut XGreedyRowCold,
) {
    // The reference loop indexes the slices with these offsets.
    assert!(seq1_index >= 0 && seq2_index >= 0);
    let index = x_first_mismatch::<REVERSE>(seq1, seq2, len1, len2, seq1_index, seq2_index);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:331-334,356-359
    // ```c
    // if (seq2_index < len2 && seq2[len2-1 - seq2_index] == FENCE_SENTRY) {
    //     ASSERT(fence_hit);
    //     *fence_hit = TRUE;
    // }
    // ```
    let stop1 = seq1_index + index;
    let stop2 = seq2_index + index;
    if stop2 < len2 {
        let at = if REVERSE { len2 - 1 - stop2 } else { stop2 };
        if *seq2.add(at as usize) == FENCE_SENTRY {
            cold.fence_hit = true;
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:571-576
    // ```c
    // if (fence_hit && *fence_hit) {
    //     return 0;
    // }
    // ```
    if cold.fence_hit {
        cold.stopped = true;
        return;
    }
    if index > cold.longest_match_run {
        cold.seed_start_q = seq1_index;
        cold.seed_start_s = seq2_index;
        cold.longest_match_run = index;
    }
    *current.add(cell) = stop2;
    if stop1 + stop2 > cold.curr_extent {
        cold.curr_extent = stop1 + stop2;
        cold.best_cell = cell;
    }
    if stop2 == len2 {
        cold.diag_lower = k + 1;
        cold.end2_reached = true;
    }
    if stop1 == len1 {
        cold.end1_cell = cell;
        cold.end1_reached = true;
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ```
/// The band of diagonals is cut into blocks of this many cells for the three passes.
/// Pass 3 still visits the diagonals in increasing `k`, as the C loop does.
/// Cells per pass over the band.
const X_GREEDY_BLOCK: usize = 256;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
///             seq2_index = MAX(last_seq2_off[d - 1][k + 1],
///                              last_seq2_off[d - 1][k    ]) + 1;
///             seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
///             seq1_index = seq2_index + k - diag_origin;
///             if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {
/// ...
///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
/// ...
///             last_seq2_off[d][k] = seq2_index;
/// ```
/// One distance `d` of the C loop above. Pass 1 is the first five statements, pass 2
/// is the slide, pass 3 is the rest in the order of the diagonals. The comments below
/// say which cases pass 2 may settle by itself and which go to `x_greedy_cell_slow`.
/// One distance. `previous` holds diagonals `tmp_diag_lower - 1 ..=
/// tmp_diag_upper + 1` of distance `d - 1`, `current` diagonals
/// `tmp_diag_lower ..= tmp_diag_upper` of distance `d`. Returns true when the
/// reference loop returns from inside the distance (fence).
///
/// Three passes per block of the band.
///
/// 1. For every diagonal, the offset the reference loop starts its slide
///    from, or -1 when the diagonal fails the X-drop test. This reads only
///    distance `d - 1`, which this distance never writes.
/// 2. For every diagonal, the slide length when it can be read off one
///    eight-byte comparison, else -1. It is read off when all of this holds:
///      * the diagonal survived, and both sequences have at least eight
///        positions left,
///      * the slide stops within those eight (so it stops inside both
///        sequences: neither end test of the reference loop can fire),
///      * the byte the slide stops at in `seq2` is not the fence sentry,
///      * the slide is not longer than the longest run when the block
///        started (no seed update: the longest run only grows).
///    This pass reads the sequences only and keeps no state between
///    diagonals.
/// 3. The diagonals in the reference order: dropped diagonals as in the
///    reference loop; diagonals with a slide length from pass 2 store their
///    cell and compare their extent, which is all the reference loop does for
///    them; every other surviving diagonal goes through `x_greedy_cell_slow`.
///
/// What the reference loop assigns on every diagonal is recovered from the
/// cell where it was last assigned:
///   * `diag_upper` is `k` of the last surviving diagonal that was visited,
///     minus one when that diagonal ended at `len1`;
///   * `curr_seq2_index` / `curr_diag` are the stored value and the diagonal
///     of the cell that last raised `curr_extent`.
///
/// SAFETY: `seq1` and `seq2` must be readable for `len1` and `len2` bytes.
#[inline(always)]
unsafe fn x_greedy_row_impl<const REVERSE: bool>(
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    previous: &[i32],
    current: &mut [i32],
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: &mut XGreedyRowState,
) -> bool {
    const NOT_A_BASE: u64 = 0xFCFC_FCFC_FCFC_FCFC;
    let band_len = current.len();
    debug_assert!(previous.len() >= band_len + 2);
    let previous = previous.as_ptr();
    let current = current.as_mut_ptr();
    // seq1_index = seq2_index + k - diag_origin = seq2_index + first + cell
    let first = tmp_diag_lower - diag_origin;
    // Pass 2 reads eight bytes at offset 0 for the diagonals it does not
    // handle, so both sequences must have eight bytes. With the fence already
    // hit, no diagonal may skip `x_greedy_cell_slow`: the reference loop
    // returns at its first surviving diagonal.
    let inline_ok = len1 >= 8 && len2 >= 8 && !st.fence_hit;
    let mut cold = XGreedyRowCold {
        diag_lower: st.diag_lower,
        end1_reached: st.end1_reached,
        end2_reached: st.end2_reached,
        curr_extent: st.curr_extent,
        best_cell: X_NO_CELL,
        end1_cell: X_NO_CELL,
        longest_match_run: st.longest_match_run,
        seed_start_q: st.seed_start_q,
        seed_start_s: st.seed_start_s,
        fence_hit: st.fence_hit,
        stopped: false,
    };
    let mut last_valid = X_NO_CELL;

    // Both arrays are written for `0..block_len` before they are read.
    let mut start = [std::mem::MaybeUninit::<i32>::uninit(); X_GREEDY_BLOCK];
    let mut slide = [std::mem::MaybeUninit::<i32>::uninit(); X_GREEDY_BLOCK];
    let start = start.as_mut_ptr() as *mut i32;
    let slide = slide.as_mut_ptr() as *mut i32;
    let mut block = 0usize;
    'row: while block < band_len {
        let block_len = (band_len - block).min(X_GREEDY_BLOCK);
        let block_first = first + block as i32;

        // pass 1
        for i in 0..block_len {
            let p = previous.add(block + i);
            let mut seq2_index = (*p.add(2)).max(*p.add(1)) + 1;
            seq2_index = seq2_index.max(*p);
            let extent = seq2_index + seq2_index + block_first + i as i32;
            *start.add(i) = if seq2_index < 0 || extent < xdrop_score {
                -1
            } else {
                seq2_index
            };
        }

        // pass 2
        if !inline_ok {
            for i in 0..block_len {
                *slide.add(i) = -1;
            }
        } else {
            let longest_match_run = cold.longest_match_run;
            for i in 0..block_len {
                let seq2_index = *start.add(i);
                let seq1_index = seq2_index + block_first + i as i32;
                // Sign bit set when the diagonal was dropped, when either
                // sequence has fewer than eight positions left, or when an
                // offset is negative.
                let unfit =
                    (len1 - 8 - seq1_index) | (len2 - 8 - seq2_index) | seq1_index | seq2_index;
                let at1 = if unfit < 0 { 0 } else { seq1_index };
                let at2 = if unfit < 0 { 0 } else { seq2_index };
                let (w1, w2) = if REVERSE {
                    (
                        u64::from_le(
                            (seq1.add((len1 - 8 - at1) as usize) as *const u64).read_unaligned(),
                        ),
                        u64::from_le(
                            (seq2.add((len2 - 8 - at2) as usize) as *const u64).read_unaligned(),
                        ),
                    )
                } else {
                    (
                        u64::from_le((seq1.add(at1 as usize) as *const u64).read_unaligned()),
                        u64::from_le((seq2.add(at2 as usize) as *const u64).read_unaligned()),
                    )
                };
                let diff = (w1 ^ w2) | (w1 & NOT_A_BASE);
                // 64 when all eight positions match
                let bits = if REVERSE {
                    diff.leading_zeros()
                } else {
                    diff.trailing_zeros()
                };
                let shift = bits & 0x38;
                let stop_byte = if REVERSE {
                    (w2 >> (56 - shift)) as u8
                } else {
                    (w2 >> shift) as u8
                };
                let index = (bits >> 3) as i32;
                let other = (unfit < 0)
                    | (diff == 0)
                    | (stop_byte == FENCE_SENTRY)
                    | (index > longest_match_run);
                *slide.add(i) = if other { -1 } else { index };
            }
        }

        // pass 3
        let mut diag_lower = cold.diag_lower;
        let mut curr_extent = cold.curr_extent;
        let mut best_cell = cold.best_cell;
        for i in 0..block_len {
            let cell = block + i;
            let seq2_index = *start.add(i);
            if seq2_index < 0 {
                if tmp_diag_lower + cell as i32 == diag_lower {
                    diag_lower += 1;
                } else {
                    *current.add(cell) = INVALID_OFFSET;
                }
                continue;
            }
            let index = *slide.add(i);
            if index < 0 {
                cold.diag_lower = diag_lower;
                cold.curr_extent = curr_extent;
                cold.best_cell = best_cell;
                x_greedy_cell_slow::<REVERSE>(
                    seq1,
                    seq2,
                    len1,
                    len2,
                    seq2_index + block_first + i as i32,
                    seq2_index,
                    cell,
                    tmp_diag_lower + cell as i32,
                    current,
                    &mut cold,
                );
                if cold.stopped {
                    last_valid = cell;
                    break 'row;
                }
                diag_lower = cold.diag_lower;
                curr_extent = cold.curr_extent;
                best_cell = cold.best_cell;
                continue;
            }
            let stop2 = seq2_index + index;
            *current.add(cell) = stop2;
            // seq1_index + index + seq2_index + index
            let extent = stop2 + stop2 + block_first + i as i32;
            if extent > curr_extent {
                curr_extent = extent;
                best_cell = cell;
            }
        }
        cold.diag_lower = diag_lower;
        cold.curr_extent = curr_extent;
        cold.best_cell = best_cell;

        // last surviving diagonal of this block
        let mut j = block_len;
        while j > 0 {
            j -= 1;
            if *start.add(j) >= 0 {
                last_valid = block + j;
                break;
            }
        }
        block += block_len;
    }

    st.diag_lower = cold.diag_lower;
    if last_valid != X_NO_CELL {
        let k = tmp_diag_lower + last_valid as i32;
        st.diag_upper = if cold.end1_cell == last_valid {
            k - 1
        } else {
            k
        };
    }
    st.end1_reached = cold.end1_reached;
    st.end2_reached = cold.end2_reached;
    st.curr_extent = cold.curr_extent;
    if cold.best_cell != X_NO_CELL {
        st.curr_seq2_index = *current.add(cold.best_cell);
        st.curr_diag = tmp_diag_lower + cold.best_cell as i32;
    }
    st.longest_match_run = cold.longest_match_run;
    st.seed_start_q = cold.seed_start_q;
    st.seed_start_s = cold.seed_start_s;
    st.fence_hit = cold.fence_hit;
    cold.stopped
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537,571-573
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ...
///             index = s_FindFirstMismatch(seq1, seq2, len1, len2,
/// ```
/// Chooses the `REVERSE` instance of `x_greedy_row_impl` (the C `reverse` argument).
///
/// SAFETY: `seq1` and `seq2` must be readable for `len1` and `len2` bytes.
#[inline(always)]
unsafe fn x_greedy_row_any(
    reverse: bool,
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    previous: &[i32],
    current: &mut [i32],
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: &mut XGreedyRowState,
) -> bool {
    if reverse {
        x_greedy_row_impl::<true>(
            seq1,
            seq2,
            len1,
            len2,
            previous,
            current,
            tmp_diag_lower,
            diag_origin,
            xdrop_score,
            st,
        )
    } else {
        x_greedy_row_impl::<false>(
            seq1,
            seq2,
            len1,
            len2,
            previous,
            current,
            tmp_diag_lower,
            diag_origin,
            xdrop_score,
            st,
        )
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ```
/// The same integer code as `x_greedy_row_any`, so the same result as the C loop; the
/// compiler may only use wider registers and `tzcnt` / `lzcnt`.
///
/// The same code compiled for AVX2/BMI (wider first pass, `tzcnt`/`lzcnt`).
///
/// SAFETY: as `x_greedy_row_any`, and the CPU must support the features.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,bmi1,bmi2,lzcnt")]
unsafe fn x_greedy_row_avx2(
    reverse: bool,
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    previous: &[i32],
    current: &mut [i32],
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: &mut XGreedyRowState,
) -> bool {
    x_greedy_row_any(
        reverse,
        seq1,
        seq2,
        len1,
        len2,
        previous,
        current,
        tmp_diag_lower,
        diag_origin,
        xdrop_score,
        st,
    )
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537
/// ```c
///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
/// ```
/// Picks the AVX2 or the plain build of the same row code (the CPU is checked once at
/// run time). Both give the same result as the C loop; this is a choice of machine code.
///
/// SAFETY: `seq1` and `seq2` must be readable for `len1` and `len2` bytes.
#[inline(never)]
unsafe fn x_greedy_row(
    reverse: bool,
    seq1: *const u8,
    seq2: *const u8,
    len1: i32,
    len2: i32,
    previous: &[i32],
    current: &mut [i32],
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: &mut XGreedyRowState,
) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        use std::sync::OnceLock;
        static AVX2: OnceLock<bool> = OnceLock::new();
        if *AVX2.get_or_init(|| {
            std::is_x86_feature_detected!("avx2")
                && std::is_x86_feature_detected!("bmi1")
                && std::is_x86_feature_detected!("bmi2")
                && std::is_x86_feature_detected!("lzcnt")
        }) {
            return x_greedy_row_avx2(
                reverse,
                seq1,
                seq2,
                len1,
                len2,
                previous,
                current,
                tmp_diag_lower,
                diag_origin,
                xdrop_score,
                st,
            );
        }
    }
    x_greedy_row_any(
        reverse,
        seq1,
        seq2,
        len1,
        len2,
        previous,
        current,
        tmp_diag_lower,
        diag_origin,
        xdrop_score,
        st,
    )
}

/// Non-affine greedy alignment with optional traceback.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:379-751 (BLAST_GreedyAlign)
fn blast_greedy_align(
    seq1: &[u8],
    len1: i32,
    seq2: &[u8],
    len2: i32,
    reverse: bool,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    seq1_align_len: &mut i32,
    seq2_align_len: &mut i32,
    max_dist: i32,
    non_affine_mem: &mut GreedyNonAffineMem,
    mut edit_block: Option<&mut GapPrelimEditBlock>,
    rem: u8,
    fence_hit: &mut bool,
    seed: &mut GreedySeed,
) -> i32 {
    let mut seq1_index: i32;
    let mut seq2_index: i32;
    let mut index: i32;
    let mut d: i32;
    let mut k: i32;
    let mut diag_lower: i32;
    let mut diag_upper: i32;
    let diag_origin: i32 = max_dist + 2;
    let mut best_dist: i32 = 0;
    let mut best_diag: i32 = 0;
    let mut longest_match_run: i32;
    let mut end1_reached = false;
    let mut end2_reached = false;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:211-225
    // ```c
    // gamp->last_seq2_off[0] =
    //    (Int4*) malloc((max_d + max_d + 6) * sizeof(Int4) * 2);
    // gamp->last_seq2_off[1] = gamp->last_seq2_off[0] + max_d + max_d + 6;
    // ```
    // NCBI preallocates two full-width rows of `2 * max_dist + 6` cells.
    let base_row_len = (2 * max_dist + 6) as usize;
    let store_traceback = edit_block.is_some();
    let rows = if store_traceback {
        (max_dist + 2) as usize
    } else {
        2
    };
    non_affine_mem.ensure_base_capacity(base_row_len);
    non_affine_mem.ensure_max_score_capacity(
        (max_dist as usize)
            + (((xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1) as usize)
            + 1,
    );
    if store_traceback {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:479-490
        // ```c
        // mem_pool = aux_data->space;
        // ...
        // else {
        //     s_RefreshMBSpace(mem_pool);
        // }
        // ```
        // Traceback row allocations reuse the same pool after rewinding
        // space_used, preserving NCBI's allocation timing and stale cells.
        non_affine_mem.refresh_traceback_pool();
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:216-228
    // ```c
    // gamp->last_seq2_off[0] =
    //    (Int4*) malloc((max_d + max_d + 6) * sizeof(Int4) * 2);
    // gamp->last_seq2_off[1] = gamp->last_seq2_off[0] + max_d + max_d + 6;
    // ```
    // Only the two base rows exist initially. Reserve the descriptor capacity
    // once, then retain each traceback row when NCBI allocates it. Unused
    // distance slots have no cells to read, persist, or destroy.
    let mut last_seq2_off = Vec::with_capacity(rows);
    last_seq2_off.push(NonAffineGreedyRow::from_base(0, base_row_len, 0));
    last_seq2_off.push(NonAffineGreedyRow::from_base(0, base_row_len, 1));

    let xdrop_offset = (xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1;
    let max_score_len = (max_dist as usize) + (xdrop_offset as usize) + 1;
    non_affine_mem.ensure_max_score_capacity(max_score_len);
    let max_score_offset = xdrop_offset as usize;
    for index in 0..max_score_offset {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:496-498
        // ```c
        // max_score = aux_data->max_score + xdrop_offset;
        // for (index = 0; index < xdrop_offset; index++)
        //     aux_data->max_score[index] = 0;
        // ```
        non_affine_mem.max_score[index] = 0;
    }

    index = find_first_mismatch_greedy(seq1, seq2, len1, len2, 0, 0, reverse, rem, fence_hit);

    *seq1_align_len = index;
    *seq2_align_len = index;
    seq1_index = index;

    seed.start_q = 0;
    seed.start_s = 0;
    longest_match_run = index;
    seed.match_length = longest_match_run;

    if index == len1 || index == len2 {
        if let Some(block) = edit_block.as_mut() {
            block.add(GapAlignOpType::Sub, index);
        }
        return 0;
    }

    last_seq2_off[0].set(non_affine_mem, diag_origin, seq1_index);
    non_affine_mem.max_score[max_score_offset] = seq1_index * match_cost;
    diag_lower = diag_origin - 1;
    diag_upper = diag_origin + 1;

    let mut converged = false;
    for d_val in 1..=max_dist {
        d = d_val;
        let mut curr_score: i32;
        let mut curr_extent: i32 = 0;
        let mut curr_seq2_index: i32 = 0;
        let mut curr_diag: i32 = 0;
        let tmp_diag_lower = diag_lower;
        let tmp_diag_upper = diag_upper;

        let prev_row_idx = if store_traceback {
            (d - 1) as usize
        } else {
            ((d - 1) & 1) as usize
        };
        let row_idx = if store_traceback {
            d as usize
        } else {
            (d & 1) as usize
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:548-550,657-678
        // ```c
        // seq2_index = MAX(last_seq2_off[d - 1][k + 1],
        //                  last_seq2_off[d - 1][k]) + 1;
        // seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
        // if (edit_block == NULL) last_seq2_off[d + 1] = last_seq2_off[d - 1];
        // ```
        // Consecutive distances use distinct rows, including the two rolling
        // rows without traceback. Borrow their descriptors once per distance.
        let previous_row = last_seq2_off[prev_row_idx];
        let current_row = last_seq2_off[row_idx];

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:523-526
        // ```c
        // last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
        // last_seq2_off[d - 1][diag_lower] = kInvalidOffset;
        // last_seq2_off[d - 1][diag_upper] = kInvalidOffset;
        // last_seq2_off[d - 1][diag_upper+1] = kInvalidOffset;
        // ```
        previous_row.set(non_affine_mem, diag_lower - 1, INVALID_OFFSET);
        previous_row.set(non_affine_mem, diag_lower, INVALID_OFFSET);
        previous_row.set(non_affine_mem, diag_upper, INVALID_OFFSET);
        previous_row.set(non_affine_mem, diag_upper + 1, INVALID_OFFSET);

        let xdrop_score = {
            let raw = non_affine_mem.max_score[d as usize] + (match_cost + mismatch_cost) * d
                - xdrop_threshold;
            ((raw as f64) / (match_cost as f64 / 2.0)).ceil() as i32
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:537-550,562,589,673-678
        // ```c
        // for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
        //     seq2_index = MAX(last_seq2_off[d - 1][k + 1],
        //                      last_seq2_off[d - 1][k]) + 1;
        //     seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
        //     /* ... */ last_seq2_off[d][k] = kInvalidOffset;
        //     /* ... */ last_seq2_off[d][k] = seq2_index;
        // }
        // last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
        //                          (diag_upper - diag_lower + 7) / 3);
        // last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
        // ```
        // Each allocated row covers its band plus two cells on each side;
        // the next band expands by at most one. The previous slice therefore
        // includes every k-1/k/k+1, and the current slice every writable k.
        // Check these ranges once without changing lengths, cells or storage.
        let band_len = (tmp_diag_upper - tmp_diag_lower + 1) as usize;
        let previous_start = (tmp_diag_lower - 1 - previous_row.origin) as usize;
        let current_start = (tmp_diag_lower - current_row.origin) as usize;
        let (previous_values, current_values) =
            non_affine_mem.row_pair_mut(previous_row, current_row);
        let previous = &previous_values[previous_start..previous_start + band_len + 2];
        let current = &mut current_values[current_start..current_start + band_len];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-615
        // ```c
        // for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
        //     seq2_index = MAX(last_seq2_off[d - 1][k + 1],
        //                      last_seq2_off[d - 1][k    ]) + 1;
        // ...
        //     index = s_FindFirstMismatch(seq1, seq2, len1, len2,
        //                                 seq1_index, seq2_index,
        //                                 fence_hit, reverse, rem);
        //     if(fence_hit && *fence_hit){
        //         return 0;
        //     }
        // ...
        // }   /* end loop over diagonals */
        // ```
        // Dispatch point of LOSAT_X_GREEDYFAST / LOSAT_X_GREEDYSHADOW. The `for` loop that follows
        // the `if !x_done` below is the port of this C loop and runs when no switch is set, when
        // `rem != 4`, or when the lengths do not fit the slices. With the switch, `x_greedy_row`
        // does the same distance; a fence stop returns 0 as the C code does. The loop over
        // distances `d` and everything after the diagonals (scores, convergence) is shared.
        // EXPERIMENT (LOSAT_X_GREEDYFAST / LOSAT_X_GREEDYSHADOW)
        let x_mode = x_greedy_mode();
        let mut x_done = false;
        let mut x_shadow: Option<(XGreedyRowState, Vec<i32>, bool)> = None;
        if x_mode != 0
            && rem == 4
            && len1 >= 0
            && len2 >= 0
            && len1 as usize <= seq1.len()
            && len2 as usize <= seq2.len()
        {
            let mut st = XGreedyRowState {
                diag_lower,
                diag_upper,
                end1_reached,
                end2_reached,
                curr_extent,
                curr_seq2_index,
                curr_diag,
                longest_match_run,
                seed_start_q: seed.start_q,
                seed_start_s: seed.start_s,
                fence_hit: *fence_hit,
            };
            if x_mode == 1 {
                // SAFETY: both lengths were checked against the slices above.
                let stopped = unsafe {
                    x_greedy_row(
                        reverse,
                        seq1.as_ptr(),
                        seq2.as_ptr(),
                        len1,
                        len2,
                        previous,
                        current,
                        tmp_diag_lower,
                        diag_origin,
                        xdrop_score,
                        &mut st,
                    )
                };
                diag_lower = st.diag_lower;
                diag_upper = st.diag_upper;
                end1_reached = st.end1_reached;
                end2_reached = st.end2_reached;
                curr_extent = st.curr_extent;
                curr_seq2_index = st.curr_seq2_index;
                curr_diag = st.curr_diag;
                longest_match_run = st.longest_match_run;
                seed.start_q = st.seed_start_q;
                seed.start_s = st.seed_start_s;
                seed.match_length = longest_match_run;
                *fence_hit = st.fence_hit;
                if stopped {
                    return 0;
                }
                x_done = true;
            } else {
                let mut copy = current.to_vec();
                // SAFETY: both lengths were checked against the slices above.
                let stopped = unsafe {
                    x_greedy_row(
                        reverse,
                        seq1.as_ptr(),
                        seq2.as_ptr(),
                        len1,
                        len2,
                        previous,
                        &mut copy,
                        tmp_diag_lower,
                        diag_origin,
                        xdrop_score,
                        &mut st,
                    )
                };
                x_shadow = Some((st, copy, stopped));
            }
        }
        if !x_done {
            for (cell_index, (previous, cell)) in
                previous.windows(3).zip(current.iter_mut()).enumerate()
            {
                k = tmp_diag_lower + cell_index as i32;
                seq2_index = previous[2].max(previous[1]) + 1;
                seq2_index = seq2_index.max(previous[0]);
                seq1_index = seq2_index + k - diag_origin;

                if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                    if k == diag_lower {
                        diag_lower += 1;
                    } else {
                        *cell = INVALID_OFFSET;
                    }
                    continue;
                }

                diag_upper = k;

                index = find_first_mismatch_greedy(
                    seq1, seq2, len1, len2, seq1_index, seq2_index, reverse, rem, fence_hit,
                );
                if *fence_hit {
                    // NCBI reference: c++/src/algo/blast/core/greedy_align.c:523-526,571-576
                    // last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
                    // if (fence_hit && *fence_hit) { return 0; }
                    // NCBI's row writes already reside in persistent scratch.
                    return 0;
                }

                if index > longest_match_run {
                    seed.start_q = seq1_index;
                    seed.start_s = seq2_index;
                    longest_match_run = index;
                    seed.match_length = longest_match_run;
                }
                seq1_index += index;
                seq2_index += index;

                *cell = seq2_index;

                if seq1_index + seq2_index > curr_extent {
                    curr_extent = seq1_index + seq2_index;
                    curr_seq2_index = seq2_index;
                    curr_diag = k;
                }

                if seq2_index == len2 {
                    diag_lower = k + 1;
                    end2_reached = true;
                }
                if seq1_index == len1 {
                    diag_upper = k - 1;
                    end1_reached = true;
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:574-576,589,596-614
        // ```c
        //             if(fence_hit && *fence_hit){
        //                 return 0;
        //             }
        // ...
        //             last_seq2_off[d][k] = seq2_index;
        // ...
        //             if (seq1_index + seq2_index > curr_extent) {
        // ```
        // LOSAT_X_GREEDYSHADOW: after the ported loop has run, the state and the row it left are
        // compared with what the fast row left on its copy of the row. The reference loop did
        // not return inside this distance, so the fast row must not have stopped either.
        if let Some((st, copy, stopped)) = x_shadow {
            // LOSAT_X_GREEDYSHADOW: the fast row must leave exactly what the
            // reference loop left (which did not return inside this row).
            assert!(
                !stopped,
                "LOSAT_X_GREEDYSHADOW: only the fast row stopped at d={d}"
            );
            let reference = XGreedyRowState {
                diag_lower,
                diag_upper,
                end1_reached,
                end2_reached,
                curr_extent,
                curr_seq2_index,
                curr_diag,
                longest_match_run,
                seed_start_q: seed.start_q,
                seed_start_s: seed.start_s,
                fence_hit: *fence_hit,
            };
            assert!(
                st == reference,
                "LOSAT_X_GREEDYSHADOW: state differs at d={d}: fast {st:?} reference {reference:?}"
            );
            assert!(
                copy[..] == current[..],
                "LOSAT_X_GREEDYSHADOW: row differs at d={d}"
            );
        }

        curr_score = curr_extent * (match_cost / 2) - d * (match_cost + mismatch_cost);
        let prev_max = non_affine_mem.max_score[(d as usize - 1) + max_score_offset];
        if curr_score >= prev_max {
            non_affine_mem.max_score[(d as usize) + max_score_offset] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        } else {
            non_affine_mem.max_score[(d as usize) + max_score_offset] = prev_max;
        }

        if diag_lower > diag_upper {
            converged = true;
            break;
        }

        if !end2_reached {
            diag_lower -= 1;
        }
        if !end1_reached {
            diag_upper += 1;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:666-678
        // ```c
        // if (edit_block == NULL) {
        //     last_seq2_off[d + 1] = last_seq2_off[d - 1];
        // } else {
        //     last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
        //                              (diag_upper - diag_lower + 7) / 3);
        //     last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
        // }
        // ```
        if store_traceback {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:666-678
            // ```c
            // last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
            //                          (diag_upper - diag_lower + 7) / 3);
            // last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
            // ```
            // s_GetMBSpace allocates SGreedyOffset cells; each cell stores
            // three Int4 values, so the reachable Int4 row length is rounded
            // down to that 3-Int4 allocation unit exactly as in NCBI.
            let row_len = (((diag_upper - diag_lower + 7) / 3) * 3).max(0) as usize;
            let pool_start = non_affine_mem.alloc_traceback_row(row_len);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:673-678
            // ```c
            // last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
            //                          (diag_upper - diag_lower + 7) / 3);
            // last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
            // ```
            // At distance d, rows 0..=d already exist. Append d+1, including
            // the last nonconverged iteration, preserving pool allocation and
            // cell contents on retries and fence returns.
            debug_assert_eq!(last_seq2_off.len(), (d + 1) as usize);
            last_seq2_off.push(NonAffineGreedyRow::from_pool(
                diag_lower - 2,
                row_len,
                pool_start,
            ));
        }
    }

    if !converged {
        return -1;
    }

    if edit_block.is_none() {
        return best_dist;
    }

    d = best_dist;
    seq2_index = *seq2_align_len;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:697-751
    // ```c
    // d = best_dist;
    // seq2_index = *seq2_align_len;
    //
    // if (fence_hit && *fence_hit)
    //     goto done;
    // ...
    // done:
    // GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
    //                       last_seq2_off[0][diag_origin]);
    // ```
    // Even when a fence was hit, NCBI still writes the final substitution
    // block before returning from traceback.
    if *fence_hit {
        edit_block.as_mut().unwrap().add(
            GapAlignOpType::Sub,
            last_seq2_off[0].get(non_affine_mem, diag_origin),
        );
        return best_dist;
    }

    while d > 0 {
        let mut new_seq2_index: i32 = 0;
        let new_diag = get_next_non_affine_tback(
            &last_seq2_off,
            non_affine_mem,
            d,
            best_diag,
            &mut new_seq2_index,
        );

        if new_diag == best_diag {
            if seq2_index - new_seq2_index > 0 {
                edit_block
                    .as_mut()
                    .unwrap()
                    .add(GapAlignOpType::Sub, seq2_index - new_seq2_index);
            }
        } else if new_diag < best_diag {
            if seq2_index - new_seq2_index > 0 {
                edit_block
                    .as_mut()
                    .unwrap()
                    .add(GapAlignOpType::Sub, seq2_index - new_seq2_index);
            }
            edit_block.as_mut().unwrap().add(GapAlignOpType::Ins, 1);
        } else {
            if seq2_index - new_seq2_index - 1 > 0 {
                edit_block
                    .as_mut()
                    .unwrap()
                    .add(GapAlignOpType::Sub, seq2_index - new_seq2_index - 1);
            }
            edit_block.as_mut().unwrap().add(GapAlignOpType::Del, 1);
        }

        d -= 1;
        best_diag = new_diag;
        seq2_index = new_seq2_index;
    }

    edit_block.as_mut().unwrap().add(
        GapAlignOpType::Sub,
        last_seq2_off[0].get(non_affine_mem, diag_origin),
    );

    best_dist
}

/// Affine greedy alignment with optional traceback.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:755-1249 (BLAST_AffineGreedyAlign)
fn blast_affine_greedy_align(
    seq1: &[u8],
    len1: i32,
    seq2: &[u8],
    len2: i32,
    reverse: bool,
    mut xdrop_threshold: i32,
    mut match_score: i32,
    mut mismatch_score: i32,
    mut in_gap_open: i32,
    mut in_gap_extend: i32,
    seq1_align_len: &mut i32,
    seq2_align_len: &mut i32,
    max_dist: i32,
    affine_mem: &mut GreedyAffineMem,
    non_affine_mem: &mut GreedyNonAffineMem,
    mut edit_block: Option<&mut GapPrelimEditBlock>,
    rem: u8,
    fence_hit: &mut bool,
    seed: &mut GreedySeed,
) -> i32 {
    let mut seq1_index: i32;
    let mut seq2_index: i32;
    let mut index: i32;
    let mut d: i32;
    let mut k: i32;
    let mut longest_match_run: i32;

    if match_score % 2 == 1 {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold *= 2;
        in_gap_open *= 2;
        in_gap_extend *= 2;
    }

    if in_gap_open == 0 && in_gap_extend == 0 {
        return blast_greedy_align(
            seq1,
            len1,
            seq2,
            len2,
            reverse,
            xdrop_threshold,
            match_score,
            mismatch_score,
            seq1_align_len,
            seq2_align_len,
            max_dist,
            non_affine_mem,
            edit_block,
            rem,
            fence_hit,
            seed,
        );
    }

    let match_score_half = match_score / 2;
    let mut op_cost = match_score + mismatch_score;
    let mut gap_open = in_gap_open;
    let mut gap_extend = in_gap_extend + match_score_half;
    let score_common_factor = gdb3(&mut op_cost, &mut gap_open, &mut gap_extend);
    let gap_open_extend = gap_open + gap_extend;
    let max_penalty = op_cost.max(gap_open_extend);

    let scaled_max_dist = max_dist * gap_extend;
    let diag_origin = max_dist + 2;

    let xdrop_offset = (xdrop_threshold + match_score_half) / score_common_factor + 1;
    let max_score_len = (scaled_max_dist as usize) + (xdrop_offset as usize) + 2;
    let max_score_offset = xdrop_offset as usize;

    let diag_len = (scaled_max_dist + max_penalty + 2) as usize;
    let diag_offset = max_penalty;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:848-856
    // ```c
    // max_dist = aux_data->max_dist;
    // scaled_max_dist = max_dist * gap_extend;
    // diag_origin = max_dist + 2;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:920-935
    // ```c
    // max_score = aux_data->max_score + xdrop_offset;
    // diag_lower = aux_data->diag_bounds;
    // diag_upper = aux_data->diag_bounds + scaled_max_dist + 1 + max_penalty;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:920-942
    // ```c
    //     max_score = aux_data->max_score + xdrop_offset;
    //     for (index = 0; index < xdrop_offset; index++)
    //         aux_data->max_score[index] = 0;
    //     ...
    //     for (index = 0; index < max_penalty; index++) {
    //         diag_lower[index] = kInvalidDiag;
    //         diag_upper[index] = -kInvalidDiag;
    //     }
    // ```
    // The other entries are written for a distance before they are read.
    affine_mem.ensure_capacity(diag_len, max_score_len);
    let GreedyAffineMem {
        rows,
        diag_lower,
        diag_upper,
        max_score,
    } = affine_mem;
    let max_score_base = &mut max_score[..max_score_len];
    max_score_base[..max_score_offset].fill(0);
    let diag_lower = &mut diag_lower[..diag_len];
    let diag_upper = &mut diag_upper[..diag_len];
    diag_lower[..diag_offset as usize].fill(INVALID_DIAG);
    diag_upper[..diag_offset as usize].fill(-INVALID_DIAG);
    let mut last_seq2_off = AffineRows {
        rows,
        ring: edit_block.is_none().then_some((max_penalty + 1) as usize),
    };

    index = find_first_mismatch_greedy(seq1, seq2, len1, len2, 0, 0, reverse, rem, fence_hit);
    if *fence_hit {
        return -1;
    }

    *seq1_align_len = index;
    *seq2_align_len = index;
    seq1_index = index;

    seed.start_q = 0;
    seed.start_s = 0;
    longest_match_run = index;
    seed.match_length = longest_match_run;

    if index == len1 || index == len2 {
        if let Some(block) = edit_block.as_mut() {
            block.add(GapAlignOpType::Sub, index);
        }
        return index * match_score;
    }

    let diag0_idx = (0 + diag_offset) as usize;
    diag_lower[diag0_idx] = diag_origin;
    diag_upper[diag0_idx] = diag_origin;

    last_seq2_off.start(0, diag_origin, diag_origin);
    *last_seq2_off.get_mut(0, diag_origin) = GreedyOffset {
        insert_off: INVALID_OFFSET,
        match_off: seq1_index,
        delete_off: INVALID_OFFSET,
    };
    max_score_base[max_score_offset] = seq1_index * match_score;

    let mut best_dist: i32 = 0;
    let mut best_diag: i32 = diag_origin;
    let mut curr_diag_lower: i32 = diag_origin - 1;
    let mut curr_diag_upper: i32 = diag_origin + 1;
    let mut end1_diag: i32 = 0;
    let mut end2_diag: i32 = 0;
    let mut num_nonempty_dist: i32 = 1;
    d = 1;
    let mut converged = false;

    while d <= scaled_max_dist {
        let xdrop_raw = max_score_base[d as usize] + score_common_factor * d - xdrop_threshold;
        let mut xdrop_score = ((xdrop_raw as f64) / (match_score_half as f64)).ceil() as i32;
        if xdrop_score < 0 {
            xdrop_score = 0;
        }

        let mut curr_extent: i32 = 0;
        let mut curr_seq2_index: i32 = 0;
        let mut curr_diag: i32 = 0;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;
        last_seq2_off.start(d, tmp_diag_lower, tmp_diag_upper);

        for k_val in tmp_diag_lower..=tmp_diag_upper {
            k = k_val;

            let mut seq2_index_del = INVALID_OFFSET;
            let d_open = d - gap_open_extend;
            let idx_open = d_open + diag_offset;
            if d_open >= 0 && idx_open >= 0 && (idx_open as usize) < diag_lower.len() {
                if k + 1 >= diag_lower[idx_open as usize] && k + 1 <= diag_upper[idx_open as usize]
                {
                    seq2_index_del = last_seq2_off.get(d_open, k + 1).match_off;
                }
            }

            let d_ext = d - gap_extend;
            let idx_ext = d_ext + diag_offset;
            if d_ext >= 0 && idx_ext >= 0 && (idx_ext as usize) < diag_lower.len() {
                if k + 1 >= diag_lower[idx_ext as usize] && k + 1 <= diag_upper[idx_ext as usize] {
                    let ext_off = last_seq2_off.get(d_ext, k + 1).delete_off;
                    if ext_off > seq2_index_del {
                        seq2_index_del = ext_off;
                    }
                }
            }

            if seq2_index_del == INVALID_OFFSET {
                last_seq2_off.get_mut(d, k).delete_off = INVALID_OFFSET;
            } else {
                last_seq2_off.get_mut(d, k).delete_off = seq2_index_del + 1;
            }

            let mut seq2_index_ins = INVALID_OFFSET;
            if d_open >= 0 && idx_open >= 0 && (idx_open as usize) < diag_lower.len() {
                if k - 1 >= diag_lower[idx_open as usize] && k - 1 <= diag_upper[idx_open as usize]
                {
                    seq2_index_ins = last_seq2_off.get(d_open, k - 1).match_off;
                }
            }
            if d_ext >= 0 && idx_ext >= 0 && (idx_ext as usize) < diag_lower.len() {
                if k - 1 >= diag_lower[idx_ext as usize] && k - 1 <= diag_upper[idx_ext as usize] {
                    let ext_off = last_seq2_off.get(d_ext, k - 1).insert_off;
                    if ext_off > seq2_index_ins {
                        seq2_index_ins = ext_off;
                    }
                }
            }
            if seq2_index_ins == INVALID_OFFSET {
                last_seq2_off.get_mut(d, k).insert_off = INVALID_OFFSET;
            } else {
                last_seq2_off.get_mut(d, k).insert_off = seq2_index_ins;
            }

            seq2_index = last_seq2_off
                .get(d, k)
                .insert_off
                .max(last_seq2_off.get(d, k).delete_off);

            let d_match = d - op_cost;
            let idx_match = d_match + diag_offset;
            if d_match >= 0 && idx_match >= 0 && (idx_match as usize) < diag_lower.len() {
                if k >= diag_lower[idx_match as usize] && k <= diag_upper[idx_match as usize] {
                    seq2_index = seq2_index.max(last_seq2_off.get(d_match, k).match_off + 1);
                }
            }

            seq1_index = seq2_index + k - diag_origin;

            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else {
                    last_seq2_off.get_mut(d, k).match_off = INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;

            index = find_first_mismatch_greedy(
                seq1, seq2, len1, len2, seq1_index, seq2_index, reverse, rem, fence_hit,
            );
            if *fence_hit {
                return -1;
            }

            if index > longest_match_run {
                seed.start_q = seq1_index;
                seed.start_s = seq2_index;
                longest_match_run = index;
                seed.match_length = longest_match_run;
            }
            seq1_index += index;
            seq2_index += index;

            last_seq2_off.get_mut(d, k).match_off = seq2_index;

            if seq1_index + seq2_index > curr_extent {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }

            if seq1_index == len1 {
                curr_diag_upper = k;
                end1_diag = k - 1;
            }
            if seq2_index == len2 {
                curr_diag_lower = k;
                end2_diag = k + 1;
            }
        }

        let curr_score = curr_extent * match_score_half - d * score_common_factor;
        let prev_max = max_score_base[(d as usize - 1) + max_score_offset];
        if curr_score > prev_max {
            max_score_base[(d as usize) + max_score_offset] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        } else {
            max_score_base[(d as usize) + max_score_offset] = prev_max;
        }

        let diag_idx = (d + diag_offset) as usize;
        if curr_diag_lower <= curr_diag_upper {
            num_nonempty_dist += 1;
            diag_lower[diag_idx] = curr_diag_lower;
            diag_upper[diag_idx] = curr_diag_upper;
        } else {
            diag_lower[diag_idx] = INVALID_DIAG;
            diag_upper[diag_idx] = -INVALID_DIAG;
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:1145-1150
        // ```c
        // if (curr_diag_lower <= curr_diag_upper) {
        //     num_nonempty_dist++;
        //     diag_lower[d] = curr_diag_lower;
        //     diag_upper[d] = curr_diag_upper;
        // } else {
        //     diag_lower[d] = kInvalidDiag;
        //     diag_upper[d] = -kInvalidDiag;
        // }
        //
        // if (diag_lower[d - max_penalty] <= diag_upper[d - max_penalty])
        //     num_nonempty_dist--;
        // ```
        // `diag_lower`/`diag_upper` are already shifted by `max_penalty`,
        // so distance `d - max_penalty` lives at `d - max_penalty + diag_offset`.
        if d >= max_penalty {
            let old_idx = (d - max_penalty + diag_offset) as usize;
            if diag_lower[old_idx] <= diag_upper[old_idx] {
                num_nonempty_dist -= 1;
            }
        }

        if num_nonempty_dist == 0 {
            converged = true;
            break;
        }

        d += 1;

        let idx_goe = d - gap_open_extend + diag_offset;
        let idx_ge = d - gap_extend + diag_offset;
        let idx_op = d - op_cost + diag_offset;

        let mut lower = diag_lower[idx_goe as usize].min(diag_lower[idx_ge as usize]) - 1;
        lower = lower.min(diag_lower[idx_op as usize]);
        if end2_diag > 0 {
            lower = lower.max(end2_diag);
        }
        curr_diag_lower = lower;

        let mut upper = diag_upper[idx_goe as usize].max(diag_upper[idx_ge as usize]) + 1;
        upper = upper.max(diag_upper[idx_op as usize]);
        if end1_diag > 0 {
            upper = upper.min(end1_diag);
        }
        curr_diag_upper = upper;
    }

    if !converged {
        return -1;
    }

    if let Some(block) = edit_block.as_mut() {
        d = best_dist;
        seq2_index = *seq2_align_len;
        let mut state = GapAlignOpType::Sub;

        while d > 0 {
            if state == GapAlignOpType::Sub {
                let mut new_seq2_index = 0;
                state = get_next_affine_tback_from_match(
                    &last_seq2_off,
                    &diag_lower,
                    &diag_upper,
                    diag_offset,
                    &mut d,
                    best_diag,
                    op_cost,
                    &mut new_seq2_index,
                );
                block.add(GapAlignOpType::Sub, seq2_index - new_seq2_index);
                seq2_index = new_seq2_index;
            } else if state == GapAlignOpType::Ins {
                block.add(GapAlignOpType::Ins, 1);
                state = get_next_affine_tback_from_indel(
                    &last_seq2_off,
                    &diag_lower,
                    &diag_upper,
                    diag_offset,
                    &mut d,
                    best_diag,
                    gap_open,
                    gap_extend,
                    GapAlignOpType::Ins,
                );
                best_diag -= 1;
            } else {
                block.add(GapAlignOpType::Del, 1);
                state = get_next_affine_tback_from_indel(
                    &last_seq2_off,
                    &diag_lower,
                    &diag_upper,
                    diag_offset,
                    &mut d,
                    best_diag,
                    gap_open,
                    gap_extend,
                    GapAlignOpType::Del,
                );
                best_diag += 1;
                seq2_index -= 1;
            }
        }

        block.add(
            GapAlignOpType::Sub,
            last_seq2_off.get(0, diag_origin).match_off,
        );
    }

    max_score_base[(best_dist as usize) + max_score_offset]
}

/// Compute alignment stats from a gap edit script.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:745-818 (identity counting)
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1055-1074 (gap counts)
fn stats_from_edit_ops(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    edit_ops: &[GapEditOp],
) -> (usize, usize, usize, usize) {
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_letters = 0usize;

    let mut qi = q_start;
    let mut si = s_start;

    for op in edit_ops {
        match *op {
            GapEditOp::Sub(n) => {
                for _ in 0..n {
                    if qi < q_seq.len() && si < s_seq.len() {
                        if q_seq[qi] == s_seq[si] {
                            matches += 1;
                        } else {
                            mismatches += 1;
                        }
                    }
                    qi += 1;
                    si += 1;
                }
            }
            GapEditOp::Del(n) => {
                gap_opens += 1;
                gap_letters += n as usize;
                si += n as usize;
            }
            GapEditOp::Ins(n) => {
                gap_opens += 1;
                gap_letters += n as usize;
                qi += n as usize;
            }
        }
    }

    (matches, mismatches, gap_opens, gap_letters)
}

/// Core greedy gapped alignment result (internal).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2539-2569 (s_BlastGreedyGapAlignStructFill)
struct GreedyGappedCore {
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
    q_seed_start: i32,
    s_seed_start: i32,
    score: i32,
    edit_script: Option<GapEditScript>,
}

/// Greedy gapped alignment core (with or without traceback).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
fn greedy_gapped_alignment_internal(
    query: &[u8],
    subject: &[u8],
    subject_len: usize,
    q_off: usize,
    s_off: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    compressed_subject: bool,
    scratch: &mut GreedyAlignScratch,
    do_traceback: bool,
) -> Option<GreedyGappedCore> {
    let q_avail = query.len().saturating_sub(q_off) as i32;
    let s_avail = subject_len.saturating_sub(s_off) as i32;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2793
    // ```c
    // if (!compressed_subject) {
    //    s = subject + s_off;
    //    rem = 4;
    // } else {
    //    s = subject + s_off/4;
    //    rem = s_off % 4;
    // }
    // ```
    let rem_forward: u8 = if compressed_subject {
        (s_off % COMPRESSION_RATIO) as u8
    } else {
        4
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2833-2834
    // ```c
    // if (compressed_subject)
    //    rem = 0;
    // ```
    let rem_reverse: u8 = if compressed_subject { 0 } else { 4 };
    let subject_offset = if compressed_subject {
        s_off / COMPRESSION_RATIO
    } else {
        s_off
    };
    if subject_offset >= subject.len() {
        return None;
    }
    let subject_full = subject;
    let subject_forward = &subject_full[subject_offset..];

    if q_avail <= 0 || s_avail <= 0 {
        return None;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:323-330
    // ```c
    // max_subject_length = MIN(max_subject_length, MAX_DBSEQ_LEN);
    // max_subject_length = MIN(GREEDY_MAX_COST,
    //              max_subject_length / GREEDY_MAX_COST_FRACTION + 1);
    // gap_align->greedy_align_mem =
    //     s_BlastGreedyAlignMemAlloc(score_params, ext_params,
    //                                max_subject_length, 0);
    // ```
    // NCBI reuses `gap_align->greedy_align_mem` across HSPs and does not shrink
    // `max_dist` between calls. Keep the same grow-only behavior in scratch.
    let initial_max_dist =
        (GREEDY_MAX_COST as i32).min((subject_len as i32) / GREEDY_MAX_COST_FRACTION as i32 + 1);
    let mut max_dist = scratch.max_dist.max(initial_max_dist);

    let mut q_ext_r = 0;
    let mut s_ext_r = 0;
    let mut q_ext_l = 0;
    let mut s_ext_l = 0;

    let affine_mem = &mut scratch.affine_mem;
    let non_affine_mem = &mut scratch.non_affine_mem;
    let fwd_prelim_tback = &mut scratch.fwd_prelim_tback;
    let rev_prelim_tback = &mut scratch.rev_prelim_tback;
    if do_traceback {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2799-2802
        // ```c
        // GapPrelimEditBlockReset(fwd_prelim_tback);
        // GapPrelimEditBlockReset(rev_prelim_tback);
        // ```
        fwd_prelim_tback.reset();
        rev_prelim_tback.reset();
    }
    let mut fwd_prelim_tback = if do_traceback {
        Some(fwd_prelim_tback)
    } else {
        None
    };
    let mut rev_prelim_tback = if do_traceback {
        Some(rev_prelim_tback)
    } else {
        None
    };

    let mut fwd_start_point = GreedySeed::default();
    let mut rev_start_point = GreedySeed::default();

    let mut fence_hit = false;
    let mut score: i32;

    loop {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2808-2817
        // ```c
        // score = BLAST_AffineGreedyAlign(q, q_avail, s, s_avail, FALSE, X,
        //        score_params->reward, -score_params->penalty,
        //        score_params->gap_open, score_params->gap_extend,
        //        &q_ext_r, &s_ext_r, gap_align->greedy_align_mem,
        //        fwd_prelim_tback, rem, fence_hit, &fwd_start_point);
        // ```
        score = blast_affine_greedy_align(
            &query[q_off..],
            q_avail,
            subject_forward,
            s_avail,
            false,
            x_drop,
            reward,
            -penalty,
            gap_open,
            gap_extend,
            &mut q_ext_r,
            &mut s_ext_r,
            max_dist,
            affine_mem,
            non_affine_mem,
            fwd_prelim_tback.as_deref_mut(),
            rem_forward,
            &mut fence_hit,
            &mut fwd_start_point,
        );
        if fence_hit {
            return None;
        }
        if debug_greedy_traceback_enabled(q_off, s_off) {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2805-2831
            // ```c
            // score = BLAST_AffineGreedyAlign(..., &q_ext_r, &s_ext_r, ...);
            // if (score >= 0) break;
            // ```
            eprintln!(
                "[GREEDY_DIRECTION] right q_ext={} s_ext={} score={}",
                q_ext_r, s_ext_r, score
            );
        }
        if score >= 0 {
            scratch.max_dist = max_dist;
            break;
        }
        max_dist *= 2;
        scratch.max_dist = max_dist;
    }

    loop {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2839-2847
        // ```c
        // score1 = BLAST_AffineGreedyAlign(query, q_off,
        //         subject, s_off, TRUE, X,
        //         score_params->reward, -score_params->penalty,
        //         score_params->gap_open, score_params->gap_extend,
        //         &q_ext_l, &s_ext_l, gap_align->greedy_align_mem,
        //         rev_prelim_tback, rem, fence_hit, &rev_start_point);
        // ```
        let score_left = blast_affine_greedy_align(
            query,
            q_off as i32,
            subject_full,
            s_off as i32,
            true,
            x_drop,
            reward,
            -penalty,
            gap_open,
            gap_extend,
            &mut q_ext_l,
            &mut s_ext_l,
            max_dist,
            affine_mem,
            non_affine_mem,
            rev_prelim_tback.as_deref_mut(),
            rem_reverse,
            &mut fence_hit,
            &mut rev_start_point,
        );
        if fence_hit {
            return None;
        }
        if debug_greedy_traceback_enabled(q_off, s_off) {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2836-2864
            // ```c
            // score1 = BLAST_AffineGreedyAlign(..., &q_ext_l, &s_ext_l, ...);
            // if (score1 >= 0) { score += score1; break; }
            // ```
            eprintln!(
                "[GREEDY_DIRECTION] left q_ext={} s_ext={} score={}",
                q_ext_l, s_ext_l, score_left
            );
        }
        if score_left >= 0 {
            score += score_left;
            scratch.max_dist = max_dist;
            break;
        }
        max_dist *= 2;
        scratch.max_dist = max_dist;
    }

    if gap_open == 0 && gap_extend == 0 {
        score = (q_ext_r + s_ext_r + q_ext_l + s_ext_l) * reward / 2 - score * (reward - penalty);
    } else if reward % 2 == 1 {
        score /= 2;
    }

    let mut q_seed_start = q_off as i32;
    let mut s_seed_start = s_off as i32;
    let mut edit_script = None;

    if do_traceback {
        let rev = rev_prelim_tback.as_ref().unwrap();
        let fwd = fwd_prelim_tback.as_ref().unwrap();
        if let Some(mut esp) = prelim_edit_block_to_gap_edit_script(rev, fwd) {
            let q_start = q_off as i32 - q_ext_l;
            let s_start = s_off as i32 - s_ext_l;
            let q_end = q_off as i32 + q_ext_r;
            let s_end = s_off as i32 + s_ext_r;
            let debug_traceback = debug_greedy_traceback_enabled(q_off, s_off);
            if debug_traceback {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2876-2882
                // ```c
                // esp = Blast_PrelimEditBlockToGapEditScript(rev_prelim_tback,
                //                                      fwd_prelim_tback);
                // if (esp) s_ReduceGaps(esp, query+q_off-q_ext_l,
                //                          subject+s_off-s_ext_l, ...);
                // ```
                eprintln!(
                    "[GREEDY_TRACEBACK] start=({}, {}) q_ext_l={} q_ext_r={} s_ext_l={} s_ext_r={} rev_block={} fwd_block={} pre_reduce={}",
                    q_off,
                    s_off,
                    q_ext_l,
                    q_ext_r,
                    s_ext_l,
                    s_ext_r,
                    format_gap_prelim_edit_block_for_trace(rev),
                    format_gap_prelim_edit_block_for_trace(fwd),
                    format_gap_edit_script_for_trace(&esp)
                );
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2881 (s_ReduceGaps query+q_off-q_ext_l, subject+s_off-s_ext_l)
            reduce_gaps(
                &mut esp,
                query,
                subject,
                q_start as usize,
                q_end as usize,
                s_start as usize,
                s_end as usize,
            );
            if debug_traceback {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2669-2758
                // ```c
                // static void s_ReduceGaps(GapEditScript* esp, const Uint1 *q,
                //                         const Uint1 *s, const Uint1 *qf,
                //                         const Uint1 *sf) { ... }
                // ```
                eprintln!(
                    "[GREEDY_TRACEBACK] start=({}, {}) post_reduce={}",
                    q_off,
                    s_off,
                    format_gap_edit_script_for_trace(&esp)
                );
            }
            edit_script = Some(esp);
        }
    } else {
        let q_box_l = q_off as i32 - q_ext_l;
        let s_box_l = s_off as i32 - s_ext_l;
        let q_box_r = q_off as i32 + q_ext_r;
        let s_box_r = s_off as i32 + s_ext_r;
        let mut q_seed_start_l = q_off as i32 - rev_start_point.start_q;
        let mut s_seed_start_l = s_off as i32 - rev_start_point.start_s;
        let mut q_seed_start_r = q_off as i32 + fwd_start_point.start_q;
        let mut s_seed_start_r = s_off as i32 + fwd_start_point.start_s;
        let mut valid_seed_len_l = 0;
        let mut valid_seed_len_r = 0;

        if q_seed_start_r < q_box_r && s_seed_start_r < s_box_r {
            valid_seed_len_r = (q_box_r - q_seed_start_r)
                .min(s_box_r - s_seed_start_r)
                .min(fwd_start_point.match_length)
                / 2;
        } else {
            q_seed_start_r = q_off as i32;
            s_seed_start_r = s_off as i32;
        }

        if q_seed_start_l > q_box_l && s_seed_start_l > s_box_l {
            valid_seed_len_l = (q_seed_start_l - q_box_l)
                .min(s_seed_start_l - s_box_l)
                .min(rev_start_point.match_length)
                / 2;
        } else {
            q_seed_start_l = q_off as i32;
            s_seed_start_l = s_off as i32;
        }

        if valid_seed_len_r > valid_seed_len_l {
            q_seed_start = q_seed_start_r + valid_seed_len_r;
            s_seed_start = s_seed_start_r + valid_seed_len_r;
        } else {
            q_seed_start = q_seed_start_l - valid_seed_len_l;
            s_seed_start = s_seed_start_l - valid_seed_len_l;
        }
    }

    Some(GreedyGappedCore {
        q_start: q_off as i32 - q_ext_l,
        q_end: q_off as i32 + q_ext_r,
        s_start: s_off as i32 - s_ext_l,
        s_end: s_off as i32 + s_ext_r,
        q_seed_start,
        s_seed_start,
        score,
        edit_script,
    })
}

/// Greedy gapped alignment score-only (preliminary).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936
pub fn greedy_gapped_alignment_score_only(
    query: &[u8],
    subject: &[u8],
    subject_len: usize,
    q_off: usize,
    s_off: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GreedyAlignScratch,
) -> Option<(usize, usize, usize, usize, i32, usize, usize)> {
    let core = greedy_gapped_alignment_internal(
        query,
        subject,
        subject_len,
        q_off,
        s_off,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2793
        // ```c
        // if (!compressed_subject) {
        //    s = subject + s_off;
        //    rem = 4;
        // } else {
        //    s = subject + s_off/4;
        //    rem = s_off % 4;
        // }
        // ```
        true,
        scratch,
        false,
    )?;

    if core.q_start < 0 || core.s_start < 0 {
        return None;
    }

    Some((
        core.q_start as usize,
        core.q_end as usize,
        core.s_start as usize,
        core.s_end as usize,
        core.score,
        core.q_seed_start as usize,
        core.s_seed_start as usize,
    ))
}

/// Greedy gapped alignment with traceback.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936
pub fn greedy_gapped_alignment_with_traceback(
    query: &[u8],
    subject: &[u8],
    subject_len: usize,
    q_off: usize,
    s_off: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GreedyAlignScratch,
) -> Option<(usize, usize, usize, usize, i32, usize, Vec<GapEditOp>)> {
    let core = greedy_gapped_alignment_internal(
        query,
        subject,
        subject_len,
        q_off,
        s_off,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop,
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2798-2802
        // ```c
        // if (do_traceback) {
        //    fwd_prelim_tback = gap_align->fwd_prelim_tback;
        //    rev_prelim_tback = gap_align->rev_prelim_tback;
        //    GapPrelimEditBlockReset(fwd_prelim_tback);
        //    GapPrelimEditBlockReset(rev_prelim_tback);
        // }
        // ```
        false,
        scratch,
        true,
    )?;

    let edit_script = core.edit_script?;
    let mut edit_ops: Vec<GapEditOp> = Vec::with_capacity(edit_script.size);
    for i in 0..edit_script.size {
        let op = edit_script.op_type[i];
        let num = edit_script.num[i];
        edit_ops.push(op.to_gap_edit_op(num));
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:583-597
    // ```c
    // Blast_HSPUpdateWithTraceback(gap_align, hsp);
    //
    // if (!delete_hsp && !kGreedyTraceback) {
    //     Int4 align_length = 0;
    //     Blast_HSPGetNumIdentitiesAndPositives(..., &align_length, ...);
    //     delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
    // }
    // ```
    // For greedy traceback NCBI does not walk the edit script here to compute
    // identities/mismatches. Keep only the traceback script and alignment
    // length; identity statistics are computed later in post-traceback passes.
    let align_length = edit_ops.iter().map(|op| op.num() as usize).sum();

    Some((
        core.q_start as usize,
        core.q_end as usize,
        core.s_start as usize,
        core.s_end as usize,
        core.score,
        align_length,
        edit_ops,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:380-387,537
    /// ```c
    /// Int4 BLAST_GreedyAlign(const Uint1* seq1, Int4 len1,
    ///                        const Uint1* seq2, Int4 len2,
    ///                        Boolean reverse, Int4 xdrop_threshold,
    /// ...
    ///         for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {
    /// ```
    /// The test runs the port of `BLAST_GreedyAlign` (and `BLAST_AffineGreedyAlign`) in
    /// reference, fast and shadow mode on random inputs and compares the alignments.
    /// EXPERIMENT (LOSAT_X_GREEDYFAST / LOSAT_X_GREEDYSPEC).
    ///
    /// On random pairs of related sequences, with and without traceback:
    ///   * the fast row gives the alignment of the reference loop (compared
    ///     as a whole, and row by row in the shadow mode);
    ///   * the reference loop gives the same alignment on a scratch that
    ///     earlier, unrelated alignments have used as on a new one, which is
    ///     what running tracebacks ahead of time on other scratches relies on.
    #[test]
    fn x_fast_greedy_rows_match_reference_on_random_alignments() {
        struct Rng(u64);
        impl Rng {
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
        fn run(
            mode: u8,
            query: &[u8],
            subject: &[u8],
            q_off: usize,
            s_off: usize,
            scoring: (i32, i32, i32, i32),
            x_drop: i32,
            scratch: &mut GreedyAlignScratch,
            traceback: bool,
        ) -> Option<(
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
            i32,
            Option<(Vec<GapAlignOpType>, Vec<i32>)>,
        )> {
            X_GREEDY_TEST_MODE.with(|m| m.set(Some(mode)));
            let core = greedy_gapped_alignment_internal(
                query,
                subject,
                subject.len(),
                q_off,
                s_off,
                scoring.0,
                scoring.1,
                scoring.2,
                scoring.3,
                x_drop,
                false,
                scratch,
                traceback,
            );
            X_GREEDY_TEST_MODE.with(|m| m.set(None));
            core.map(|c| {
                (
                    c.q_start,
                    c.q_end,
                    c.s_start,
                    c.s_end,
                    c.q_seed_start,
                    c.s_seed_start,
                    c.score,
                    c.edit_script
                        .map(|e| (e.op_type[..e.size].to_vec(), e.num[..e.size].to_vec())),
                )
            })
        }

        let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1500);
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        let mut used_fast = GreedyAlignScratch::new();
        let mut used_reference = GreedyAlignScratch::new();
        let mut aligned = 0usize;
        let mut with_gaps = 0usize;
        for case in 0..cases {
            // a core, a mutated copy of it, and unrelated flanks on both
            let core_len = 20 + rng.below(if case % 10 == 0 { 6000 } else { 600 }) as usize;
            let substitution = [0u64, 2, 5, 20, 60, 150][rng.below(6) as usize]; // per mille
            let indel = [0u64, 0, 2, 10, 30][rng.below(5) as usize]; // per mille
            let ambiguity = if case % 7 == 0 { 8 } else { 0 }; // per mille
            let random_base = |rng: &mut Rng| rng.below(4) as u8;
            let core: Vec<u8> = (0..core_len).map(|_| random_base(&mut rng)).collect();
            let mut query = Vec::new();
            let mut subject = Vec::new();
            for _ in 0..rng.below(40) {
                query.push(random_base(&mut rng));
            }
            for _ in 0..rng.below(40) {
                subject.push(random_base(&mut rng));
            }
            // positions of an unchanged core base in both sequences
            let mut anchors: Vec<(usize, usize)> = Vec::new();
            for &base in &core {
                let roll = rng.below(1000);
                if roll < indel {
                    // drop from the subject
                    query.push(base);
                } else if roll < 2 * indel {
                    // extra base in the subject
                    subject.push(random_base(&mut rng));
                    anchors.push((query.len(), subject.len()));
                    query.push(base);
                    subject.push(base);
                } else if roll < 2 * indel + substitution {
                    query.push(base);
                    subject.push((base + 1 + rng.below(3) as u8) & 3);
                } else if roll < 2 * indel + substitution + ambiguity {
                    // an ambiguity code in the query (never matches)
                    query.push(4 + rng.below(12) as u8);
                    subject.push(base);
                } else {
                    anchors.push((query.len(), subject.len()));
                    query.push(base);
                    subject.push(base);
                }
            }
            for _ in 0..rng.below(40) {
                query.push(random_base(&mut rng));
            }
            for _ in 0..rng.below(40) {
                subject.push(random_base(&mut rng));
            }
            if anchors.is_empty() {
                continue;
            }
            let (q_off, s_off) = anchors[rng.below(anchors.len() as u64) as usize];
            // megablast, and two scorings that take the affine route
            let scoring = match case % 5 {
                0 => (2, -3, 5, 2),
                1 => (1, -3, 0, 0),
                _ => (1, -2, 0, 0),
            };
            let x_drop = 5 + rng.below(60) as i32;
            for traceback in [false, true] {
                let reference = run(
                    0,
                    &query,
                    &subject,
                    q_off,
                    s_off,
                    scoring,
                    x_drop,
                    &mut GreedyAlignScratch::new(),
                    traceback,
                );
                let fast = run(
                    1,
                    &query,
                    &subject,
                    q_off,
                    s_off,
                    scoring,
                    x_drop,
                    &mut used_fast,
                    traceback,
                );
                assert_eq!(
                    fast, reference,
                    "case {case}: fast rows, traceback={traceback}"
                );
                let shadow = run(
                    2,
                    &query,
                    &subject,
                    q_off,
                    s_off,
                    scoring,
                    x_drop,
                    &mut GreedyAlignScratch::new(),
                    traceback,
                );
                assert_eq!(
                    shadow, reference,
                    "case {case}: shadow, traceback={traceback}"
                );
                let reused = run(
                    0,
                    &query,
                    &subject,
                    q_off,
                    s_off,
                    scoring,
                    x_drop,
                    &mut used_reference,
                    traceback,
                );
                assert_eq!(
                    reused, reference,
                    "case {case}: used scratch, traceback={traceback}"
                );
                if let Some(result) = &reference {
                    aligned += 1;
                    if let Some((ops, _)) = &result.7 {
                        if ops.iter().any(|op| *op != GapAlignOpType::Sub) {
                            with_gaps += 1;
                        }
                    }
                }
            }
        }
        // the generator must produce real alignments, some of them gapped
        assert!(
            aligned > cases,
            "only {aligned} alignments in {cases} cases"
        );
        assert!(
            with_gaps * 20 > cases,
            "only {with_gaps} gapped alignments in {cases} cases"
        );
    }

    // NCBI reference: c++/src/algo/blast/core/greedy_align.c:500-501,523-526,571-576
    // last_seq2_off[0][diag_origin] = seq1_index;
    // last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
    // last_seq2_off[d - 1][diag_lower] = kInvalidOffset;
    // last_seq2_off[d - 1][diag_upper] = kInvalidOffset;
    // last_seq2_off[d - 1][diag_upper+1] = kInvalidOffset;
    // if (fence_hit && *fence_hit) { return 0; }
    // C writes into persistent scratch before the fence return.
    #[test]
    fn non_affine_fence_retains_written_scratch() {
        for traceback in [false, true] {
            for reverse in [false, true] {
                let mut mem = GreedyNonAffineMem::new();
                mem.ensure_base_capacity(14);
                mem.last_seq2_off[0].fill(123);
                let mut edit = GapPrelimEditBlock::new();
                let mut seed = GreedySeed::default();
                let mut fence = false;
                let (mut qlen, mut slen) = (0, 0);
                let score = blast_greedy_align(
                    &[0, 1, 0],
                    3,
                    &[0, FENCE_SENTRY, 0],
                    3,
                    reverse,
                    20,
                    2,
                    4,
                    &mut qlen,
                    &mut slen,
                    4,
                    &mut mem,
                    traceback.then_some(&mut edit),
                    4,
                    &mut fence,
                    &mut seed,
                );
                assert!(fence);
                assert_eq!((score, qlen, slen), (0, 1, 1));
                assert_eq!(mem.last_seq2_off[0][6], 1);
                for index in [4, 5, 7, 8] {
                    assert_eq!(mem.last_seq2_off[0][index], INVALID_OFFSET);
                }
                assert_eq!(mem.last_seq2_off[0][3], 123);
            }
        }
    }

    #[test]
    fn test_greedy_traceback_ncbi_regression_cases() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936
        // ```c
        // BLAST_AffineGreedyAlign(..., fwd_prelim_tback, ...);
        // BLAST_AffineGreedyAlign(..., rev_prelim_tback, ...);
        // esp = Blast_PrelimEditBlockToGapEditScript(rev_prelim_tback,
        //                                             fwd_prelim_tback);
        // if (esp) s_ReduceGaps(esp, query+q_off-q_ext_l,
        //                          subject+s_off-s_ext_l, ...);
        // ```
        let cases: &[(
            &[u8],
            &[u8],
            usize,
            usize,
            i32,
            i32,
            i32,
            i32,
            i32,
            (usize, usize, usize, usize, i32, usize),
            &[GapEditOp],
        )] = &[
            (
                &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                5,
                5,
                1,
                -2,
                0,
                0,
                54,
                (0, 11, 0, 11, 11, 11),
                &[GapEditOp::Sub(11)],
            ),
            (
                &[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                &[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                6,
                5,
                1,
                -2,
                0,
                0,
                54,
                (1, 12, 0, 11, 11, 11),
                &[GapEditOp::Sub(11)],
            ),
            (
                &[0, 1, 2, 3, 0, 1, 2],
                &[0, 1, 2, 0, 1, 2],
                3,
                3,
                1,
                -2,
                0,
                0,
                54,
                (0, 7, 0, 6, 3, 7),
                &[GapEditOp::Sub(3), GapEditOp::Ins(1), GapEditOp::Sub(3)],
            ),
            (
                &[0, 1, 2, 3, 0, 1],
                &[0, 1, 2, 3, 0, 1],
                3,
                3,
                1,
                -2,
                5,
                2,
                54,
                (0, 6, 0, 6, 6, 6),
                &[GapEditOp::Sub(6)],
            ),
            (
                &[0],
                &[0],
                0,
                0,
                1,
                -2,
                5,
                2,
                54,
                (0, 1, 0, 1, 1, 1),
                &[GapEditOp::Sub(1)],
            ),
        ];

        for (
            query,
            subject,
            q_off,
            s_off,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
            expected,
            expected_ops,
        ) in cases
        {
            let mut scratch = GreedyAlignScratch::new();
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:479-490
            // ```c
            // else { s_RefreshMBSpace(mem_pool); }
            // ```
            // Reuse the same scratch for the existing frozen cases and edits.
            for _ in 0..2 {
                let result = greedy_gapped_alignment_with_traceback(
                    query,
                    subject,
                    subject.len(),
                    *q_off,
                    *s_off,
                    *reward,
                    *penalty,
                    *gap_open,
                    *gap_extend,
                    *x_drop,
                    &mut scratch,
                )
                .expect("NCBI greedy traceback converges for the bounded fixture");
                assert_eq!(
                    (result.0, result.1, result.2, result.3, result.4, result.5),
                    *expected
                );
                assert_eq!(result.6, *expected_ops);
            }
        }
    }

    #[test]
    fn test_greedy_reverse_and_empty_boundary_cases() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:380-453
        // ```c
        // Int4 BLAST_GreedyAlign(..., Boolean reverse, ...)
        // {
        //     index = s_FindFirstMismatch(..., reverse, ...);
        //     ...
        // }
        // ```
        let reverse = greedy_align_one_direction_ex(
            &[0, 1, 2, 3],
            &[0, 1, 2, 3],
            4,
            4,
            1,
            -2,
            0,
            0,
            54,
            true,
        );
        assert_eq!(reverse, (4, 4, 4, 4, 0, 0, 0));

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2793
        // ```c
        // q_avail = query_length - q_off;
        // s_avail = subject_length - s_off;
        // if (q_avail <= 0 || s_avail <= 0) return -1;
        // ```
        let mut scratch = GreedyAlignScratch::new();
        assert!(greedy_gapped_alignment_with_traceback(
            &[],
            &[],
            0,
            0,
            0,
            1,
            -2,
            5,
            2,
            54,
            &mut scratch,
        )
        .is_none());
    }
}
