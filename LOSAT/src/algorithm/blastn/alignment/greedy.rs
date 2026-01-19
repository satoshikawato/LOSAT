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

use crate::common::GapEditOp;
use crate::core::blast_encoding::COMPRESSION_RATIO;
use super::super::constants::{GREEDY_MAX_COST, GREEDY_MAX_COST_FRACTION, INVALID_OFFSET, INVALID_DIAG};
use super::super::sequence_compare::{find_first_mismatch_ex, find_first_mismatch};
use super::utilities::gdb3;

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

/// Scratch memory for affine greedy alignment (preallocated arrays).
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/greedy_align.h:88-99
/// ```c
/// typedef struct SGreedyAlignMem {
///    SGreedyOffset** last_seq2_off_affine;
///    Int4* diag_bounds;
///    Int4* max_score;
/// } SGreedyAlignMem;
/// ```
struct GreedyAffineMem {
    last_seq2_off: Vec<Vec<GreedyOffset>>,
    diag_lower: Vec<i32>,
    diag_upper: Vec<i32>,
    max_score: Vec<i32>,
}

impl GreedyAffineMem {
    fn new() -> Self {
        Self {
            last_seq2_off: Vec::new(),
            diag_lower: Vec::new(),
            diag_upper: Vec::new(),
            max_score: Vec::new(),
        }
    }

    fn ensure_capacity(
        &mut self,
        num_rows: usize,
        array_size: usize,
        diag_len: usize,
        max_score_len: usize,
    ) {
        if self.last_seq2_off.len() < num_rows {
            self.last_seq2_off.resize_with(num_rows, Vec::new);
        }
        for row in self.last_seq2_off.iter_mut().take(num_rows) {
            if row.len() < array_size {
                row.resize(
                    array_size,
                    GreedyOffset {
                        insert_off: INVALID_OFFSET,
                        match_off: INVALID_OFFSET,
                        delete_off: INVALID_OFFSET,
                    },
                );
            }
        }
        if self.diag_lower.len() < diag_len {
            self.diag_lower.resize(diag_len, INVALID_DIAG);
        }
        if self.diag_upper.len() < diag_len {
            self.diag_upper.resize(diag_len, -INVALID_DIAG);
        }
        if self.max_score.len() < max_score_len {
            self.max_score.resize(max_score_len, 0);
        }
    }

    fn reset(
        &mut self,
        num_rows: usize,
        array_size: usize,
        diag_len: usize,
        max_score_len: usize,
    ) {
        let invalid = GreedyOffset {
            insert_off: INVALID_OFFSET,
            match_off: INVALID_OFFSET,
            delete_off: INVALID_OFFSET,
        };
        for row in self.last_seq2_off.iter_mut().take(num_rows) {
            row[..array_size].fill(invalid);
        }
        self.diag_lower[..diag_len].fill(INVALID_DIAG);
        self.diag_upper[..diag_len].fill(-INVALID_DIAG);
        self.max_score[..max_score_len].fill(0);
    }
}
/// Bookkeeping structure for affine greedy alignment (NCBI BLAST's SGreedyOffset).
/// When aligning two sequences, stores the largest offset into the second sequence
/// that leads to a high-scoring alignment for a given start point, tracking
/// different path endings separately for affine gap penalties.
#[derive(Clone, Copy, Default)]
pub struct GreedyOffset {
    insert_off: i32,  // Best offset for a path ending in an insertion (gap in seq2)
    match_off: i32,   // Best offset for a path ending in a match/mismatch
    delete_off: i32,  // Best offset for a path ending in a deletion (gap in seq1)
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
    affine_mem: GreedyAffineMem,
    fwd_prelim_tback: GapPrelimEditBlock,
    rev_prelim_tback: GapPrelimEditBlock,
}

impl GreedyAlignScratch {
    pub fn new() -> Self {
        Self {
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
                q_seq, s_seq, len1, len2, reward, penalty, gap_open, gap_extend, x_drop, max_dist, reverse,
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
        q_seq, s_seq, q_seq.len(), s_seq.len(),
        reward, penalty, gap_open, gap_extend, x_drop, false,
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
                0  // When d < xdrop_offset, use index 0 which maps to max_score[xdrop_offset] = 0
            };
            let xdrop_score = mem.max_score[xdrop_idx + xdrop_offset] + (op_cost * d as i32) - scaled_xdrop;
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
                let p_plus = if ku + 1 < array_size { mem.last_seq2_off_a[ku + 1] } else { -2 };
                let p = mem.last_seq2_off_a[ku];
                let p_minus = if ku >= 1 { mem.last_seq2_off_a[ku - 1] } else { -2 };
                (p_plus, p, p_minus)
            } else {
                let p_plus = if ku + 1 < array_size { mem.last_seq2_off_b[ku + 1] } else { -2 };
                let p = mem.last_seq2_off_b[ku];
                let p_minus = if ku >= 1 { mem.last_seq2_off_b[ku - 1] } else { -2 };
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
                    let matches = find_first_mismatch_ex(q_seq, s_seq, len1, len2, seq1_idx, seq2_idx, reverse);
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
            (reward * 2, penalty * 2, x_drop * 2, in_gap_open * 2, in_gap_extend * 2)
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
    debug_assert!(gap_extend > 0, "gap_extend must be positive, got {}", gap_extend);
    debug_assert!(max_penalty > 0, "max_penalty must be positive, got {}", max_penalty);
    debug_assert!(scaled_max_dist >= 0, "scaled_max_dist must be non-negative, got {}", scaled_max_dist);
    
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
    
    let mut last_seq2_off: Vec<Vec<GreedyOffset>> = vec![vec![GreedyOffset {
        insert_off: INVALID_OFFSET,
        match_off: INVALID_OFFSET,
        delete_off: INVALID_OFFSET,
    }; array_size]; num_rows];

    // Max score at each distance for X-drop
    // NCBI reference: greedy_align.c:872-873
    // xdrop_offset = (xdrop_threshold + match_score_half) / score_common_factor + 1;
    let xdrop_offset = ((xdrop_threshold + match_score_half) / score_common_factor.max(1) + 1) as usize;
    let max_score_size = scaled_max_dist_usize + xdrop_offset + 2;
    let mut max_score: Vec<i32> = vec![0; max_score_size];

    // Find initial run of matches
    let initial_matches = find_first_mismatch(q_seq, s_seq, 0, 0);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (initial_matches, initial_matches, score, initial_matches, 0, 0, 0),
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
        if idx < arr.len() { arr[idx] } else { INVALID_DIAG }
    };
    let get_diag_upper = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() { arr[idx] } else { -INVALID_DIAG }
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
            0  // When d < xdrop_offset, use index 0 which maps to max_score[xdrop_offset] = 0
        };
        let xdrop_score_raw = max_score[xdrop_idx + xdrop_offset] 
            + score_common_factor * d - xdrop_threshold;
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
            let mut seq2_index = last_seq2_off[du][ku].insert_off.max(last_seq2_off[du][ku].delete_off);
            
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
        let matches = alignment_len.saturating_sub(estimated_mismatches).saturating_sub(gap_letters);
        
        let score = (matches as i32) * reward 
            + (estimated_mismatches as i32) * penalty
            - (estimated_gaps as i32) * in_gap_open
            - (gap_letters as i32) * in_gap_extend;
        
        return (
            (best_seq1_len, best_seq2_len, score, matches, estimated_mismatches, estimated_gaps, gap_letters),
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
    
    let matches = alignment_len.saturating_sub(estimated_mismatches).saturating_sub(gap_letters);
    
    // Calculate final score
    let score = (matches as i32) * reward 
        + (estimated_mismatches as i32) * penalty
        - (estimated_gap_opens as i32) * in_gap_open
        - (gap_letters as i32) * in_gap_extend;

    (
        (best_seq1_len, best_seq2_len, score, matches, estimated_mismatches, estimated_gap_opens, gap_letters),
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

    for op in &rev_prelim_tback.edit_ops {
        if index >= esp.size {
            break;
        }
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

    while i >= 0 {
        let op = fwd_prelim_tback.edit_ops[i as usize];
        if index >= esp.size {
            break;
        }
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
            esp.op_type[before_idx] = if qd > 0 { GapAlignOpType::Del } else { GapAlignOpType::Ins };
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
            esp.op_type[after_idx] = if qd > 0 { GapAlignOpType::Del } else { GapAlignOpType::Ins };
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2669-2758 (q/s pointers with qf/sf end pointers)
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
                        while (q_idx - nm1 as isize) >= 0
                            && (s_idx - nm1 as isize) >= 0
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
                        while q_idx + 1 < q_len
                            && s_idx + 1 < s_len
                            && q[(q_start + q_idx) as usize] == s[(s_start + s_idx) as usize]
                        {
                            nm2 += 1;
                            q_idx += 1;
                            s_idx += 1;
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

                for _ in 0..esp.num[i - 1] {
                    if q1 >= 0 && s1 >= 0 && q1 < q_len && s1 < s_len {
                        if q[(q_start + q1) as usize] == s[(s_start + s1) as usize] {
                            nm1 += 1;
                        }
                    }
                    if q_pos >= 0 && s_pos >= 0 && q_pos < q_len && s_pos < s_len {
                        if q[(q_start + q_pos) as usize] == s[(s_start + s_pos) as usize] {
                            nm2 += 1;
                        }
                    }
                    q1 += 1;
                    s1 += 1;
                    q_pos += 1;
                    s_pos += 1;
                }

                for _ in 0..d {
                    if q_pos >= 0 && s_pos >= 0 && q_pos < q_len && s_pos < s_len {
                        if q[(q_start + q_pos) as usize] == s[(s_start + s_pos) as usize] {
                            nm2 += 1;
                        }
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
                && seq1[(len1 - 1 - seq1_index) as usize]
                    == seq2[(len2 - 1 - seq2_index) as usize]
            {
                seq1_index += 1;
                seq2_index += 1;
            }
            if seq2_index < len2 && seq2[(len2 - 1 - seq2_index) as usize] == FENCE_SENTRY {
                *fence_hit = true;
            }
        } else {
            while seq1_index < len1
                && seq2_index < len2
                && seq1[(len1 - 1 - seq1_index) as usize] == ncbi2na_unpack_base(
                    seq2[((len2 - 1 - seq2_index + rem as i32) / 4) as usize],
                    (3 - (len2 - 1 - seq2_index + rem as i32) % 4) as u8,
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
    last_seq2_off: &[Vec<GreedyOffset>],
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
            new_seq2_index = last_seq2_off[idx as usize][diag as usize].match_off;
            if new_seq2_index
                >= last_seq2_off[*d as usize][diag as usize]
                    .insert_off
                    .max(last_seq2_off[*d as usize][diag as usize].delete_off)
            {
                *d -= op_cost;
                *seq2_index = new_seq2_index;
                return GapAlignOpType::Sub;
            }
        }
    }

    if last_seq2_off[*d as usize][diag as usize].insert_off
        > last_seq2_off[*d as usize][diag as usize].delete_off
    {
        *seq2_index = last_seq2_off[*d as usize][diag as usize].insert_off;
        GapAlignOpType::Ins
    } else {
        *seq2_index = last_seq2_off[*d as usize][diag as usize].delete_off;
        GapAlignOpType::Del
    }
}

/// Affine traceback helper from insertion/deletion state.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c:199-258
fn get_next_affine_tback_from_indel(
    last_seq2_off: &[Vec<GreedyOffset>],
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
    let new_diag = if state == GapAlignOpType::Ins { diag - 1 } else { diag + 1 };

    let last_d = (*d) - gap_extend;
    let mut new_seq2_index = INVALID_OFFSET;
    let dl_idx = last_d + diag_offset;
    if last_d >= 0 && dl_idx >= 0 && (dl_idx as usize) < diag_lower.len() {
        if new_diag >= diag_lower[dl_idx as usize] && new_diag <= diag_upper[dl_idx as usize] {
            new_seq2_index = if state == GapAlignOpType::Ins {
                last_seq2_off[last_d as usize][new_diag as usize].insert_off
            } else {
                last_seq2_off[last_d as usize][new_diag as usize].delete_off
            };
        }
    }

    let last_d = (*d) - gap_open_extend;
    let dl_idx = last_d + diag_offset;
    if last_d >= 0 && dl_idx >= 0 && (dl_idx as usize) < diag_lower.len() {
        if new_diag >= diag_lower[dl_idx as usize]
            && new_diag <= diag_upper[dl_idx as usize]
            && new_seq2_index < last_seq2_off[last_d as usize][new_diag as usize].match_off
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
    last_seq2_off: &[Vec<i32>],
    d: i32,
    diag: i32,
    seq2_index: &mut i32,
) -> i32 {
    if last_seq2_off[(d - 1) as usize][(diag - 1) as usize]
        > last_seq2_off[(d - 1) as usize][diag as usize]
            .max(last_seq2_off[(d - 1) as usize][(diag + 1) as usize])
    {
        *seq2_index = last_seq2_off[(d - 1) as usize][(diag - 1) as usize];
        return diag - 1;
    }
    if last_seq2_off[(d - 1) as usize][diag as usize]
        > last_seq2_off[(d - 1) as usize][(diag + 1) as usize]
    {
        *seq2_index = last_seq2_off[(d - 1) as usize][diag as usize];
        return diag;
    }
    *seq2_index = last_seq2_off[(d - 1) as usize][(diag + 1) as usize];
    diag + 1
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

    let array_size = (2 * diag_origin + 4) as usize;
    let store_traceback = edit_block.is_some();
    let rows = if store_traceback { (max_dist + 2) as usize } else { 2 };
    let mut last_seq2_off: Vec<Vec<i32>> = vec![vec![INVALID_OFFSET; array_size]; rows];

    let xdrop_offset =
        (xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1;
    let max_score_len = (max_dist as usize) + (xdrop_offset as usize) + 2;
    let mut max_score_base = vec![0; max_score_len];
    let max_score_offset = xdrop_offset as usize;

    index = find_first_mismatch_greedy(
        seq1,
        seq2,
        len1,
        len2,
        0,
        0,
        reverse,
        rem,
        fence_hit,
    );

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

    last_seq2_off[0][diag_origin as usize] = seq1_index;
    max_score_base[max_score_offset] = seq1_index * match_cost;
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

        if (diag_lower - 1) >= 0 {
            last_seq2_off[prev_row_idx][(diag_lower - 1) as usize] = INVALID_OFFSET;
        }
        if diag_lower >= 0 {
            last_seq2_off[prev_row_idx][diag_lower as usize] = INVALID_OFFSET;
        }
        if diag_upper >= 0 && (diag_upper as usize) < array_size {
            last_seq2_off[prev_row_idx][diag_upper as usize] = INVALID_OFFSET;
        }
        if (diag_upper + 1) >= 0 && ((diag_upper + 1) as usize) < array_size {
            last_seq2_off[prev_row_idx][(diag_upper + 1) as usize] = INVALID_OFFSET;
        }

        let xdrop_score = {
            let raw = max_score_base[d as usize]
                + (match_cost + mismatch_cost) * d
                - xdrop_threshold;
            ((raw as f64) / (match_cost as f64 / 2.0)).ceil() as i32
        };

        for k_val in tmp_diag_lower..=tmp_diag_upper {
            k = k_val;
            let prev_row = &last_seq2_off[prev_row_idx];
            seq2_index = prev_row[(k + 1) as usize].max(prev_row[k as usize]) + 1;
            seq2_index = seq2_index.max(prev_row[(k - 1) as usize]);
            seq1_index = seq2_index + k - diag_origin;

            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == diag_lower {
                    diag_lower += 1;
                } else if (k as usize) < array_size {
                    last_seq2_off[row_idx][k as usize] = INVALID_OFFSET;
                }
                continue;
            }

            diag_upper = k;

            index = find_first_mismatch_greedy(
                seq1,
                seq2,
                len1,
                len2,
                seq1_index,
                seq2_index,
                reverse,
                rem,
                fence_hit,
            );
            if *fence_hit {
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

            last_seq2_off[row_idx][k as usize] = seq2_index;

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

        curr_score =
            curr_extent * (match_cost / 2) - d * (match_cost + mismatch_cost);
        let prev_max = max_score_base[(d as usize - 1) + max_score_offset];
        if curr_score >= prev_max {
            max_score_base[(d as usize) + max_score_offset] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        } else {
            max_score_base[(d as usize) + max_score_offset] = prev_max;
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
    }

    if !converged {
        return -1;
    }

    if edit_block.is_none() {
        return best_dist;
    }

    d = best_dist;
    seq2_index = *seq2_align_len;
    if *fence_hit {
        return best_dist;
    }

    while d > 0 {
        let mut new_seq2_index: i32 = 0;
        let new_diag =
            get_next_non_affine_tback(&last_seq2_off, d, best_diag, &mut new_seq2_index);

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

    edit_block
        .as_mut()
        .unwrap()
        .add(GapAlignOpType::Sub, last_seq2_off[0][diag_origin as usize]);

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
    let array_size = (2 * diag_origin + 4) as usize;

    let xdrop_offset = (xdrop_threshold + match_score_half) / score_common_factor + 1;
    let max_score_len = (scaled_max_dist as usize) + (xdrop_offset as usize) + 2;
    let max_score_offset = xdrop_offset as usize;

    let diag_len = (scaled_max_dist + max_penalty + 2) as usize;
    let diag_offset = max_penalty;
    let num_rows = (scaled_max_dist + max_penalty + 2) as usize;

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
    affine_mem.ensure_capacity(num_rows, array_size, diag_len, max_score_len);
    affine_mem.reset(num_rows, array_size, diag_len, max_score_len);
    let max_score_base = &mut affine_mem.max_score[..max_score_len];
    let diag_lower = &mut affine_mem.diag_lower[..diag_len];
    let diag_upper = &mut affine_mem.diag_upper[..diag_len];

    index = find_first_mismatch_greedy(
        seq1,
        seq2,
        len1,
        len2,
        0,
        0,
        reverse,
        rem,
        fence_hit,
    );
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

    let last_seq2_off = &mut affine_mem.last_seq2_off;

    let diag0_idx = (0 + diag_offset) as usize;
    diag_lower[diag0_idx] = diag_origin;
    diag_upper[diag0_idx] = diag_origin;

    last_seq2_off[0][diag_origin as usize].match_off = seq1_index;
    last_seq2_off[0][diag_origin as usize].insert_off = INVALID_OFFSET;
    last_seq2_off[0][diag_origin as usize].delete_off = INVALID_OFFSET;
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
        let xdrop_raw = max_score_base[d as usize]
            + score_common_factor * d
            - xdrop_threshold;
        let mut xdrop_score =
            ((xdrop_raw as f64) / (match_score_half as f64)).ceil() as i32;
        if xdrop_score < 0 {
            xdrop_score = 0;
        }

        let mut curr_extent: i32 = 0;
        let mut curr_seq2_index: i32 = 0;
        let mut curr_diag: i32 = 0;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;

        for k_val in tmp_diag_lower..=tmp_diag_upper {
            k = k_val;

            let mut seq2_index_del = INVALID_OFFSET;
            let d_open = d - gap_open_extend;
            let idx_open = d_open + diag_offset;
            if d_open >= 0 && idx_open >= 0 && (idx_open as usize) < diag_lower.len() {
                if k + 1 >= diag_lower[idx_open as usize]
                    && k + 1 <= diag_upper[idx_open as usize]
                {
                    seq2_index_del =
                        last_seq2_off[d_open as usize][(k + 1) as usize].match_off;
                }
            }

            let d_ext = d - gap_extend;
            let idx_ext = d_ext + diag_offset;
            if d_ext >= 0 && idx_ext >= 0 && (idx_ext as usize) < diag_lower.len() {
                if k + 1 >= diag_lower[idx_ext as usize]
                    && k + 1 <= diag_upper[idx_ext as usize]
                {
                    let ext_off =
                        last_seq2_off[d_ext as usize][(k + 1) as usize].delete_off;
                    if ext_off > seq2_index_del {
                        seq2_index_del = ext_off;
                    }
                }
            }

            if seq2_index_del == INVALID_OFFSET {
                last_seq2_off[d as usize][k as usize].delete_off = INVALID_OFFSET;
            } else {
                last_seq2_off[d as usize][k as usize].delete_off = seq2_index_del + 1;
            }

            let mut seq2_index_ins = INVALID_OFFSET;
            if d_open >= 0 && idx_open >= 0 && (idx_open as usize) < diag_lower.len() {
                if k - 1 >= diag_lower[idx_open as usize]
                    && k - 1 <= diag_upper[idx_open as usize]
                {
                    seq2_index_ins =
                        last_seq2_off[d_open as usize][(k - 1) as usize].match_off;
                }
            }
            if d_ext >= 0 && idx_ext >= 0 && (idx_ext as usize) < diag_lower.len() {
                if k - 1 >= diag_lower[idx_ext as usize]
                    && k - 1 <= diag_upper[idx_ext as usize]
                {
                    let ext_off =
                        last_seq2_off[d_ext as usize][(k - 1) as usize].insert_off;
                    if ext_off > seq2_index_ins {
                        seq2_index_ins = ext_off;
                    }
                }
            }
            if seq2_index_ins == INVALID_OFFSET {
                last_seq2_off[d as usize][k as usize].insert_off = INVALID_OFFSET;
            } else {
                last_seq2_off[d as usize][k as usize].insert_off = seq2_index_ins;
            }

            seq2_index = last_seq2_off[d as usize][k as usize]
                .insert_off
                .max(last_seq2_off[d as usize][k as usize].delete_off);

            let d_match = d - op_cost;
            let idx_match = d_match + diag_offset;
            if d_match >= 0 && idx_match >= 0 && (idx_match as usize) < diag_lower.len() {
                if k >= diag_lower[idx_match as usize]
                    && k <= diag_upper[idx_match as usize]
                {
                    seq2_index =
                        seq2_index.max(last_seq2_off[d_match as usize][k as usize].match_off + 1);
                }
            }

            seq1_index = seq2_index + k - diag_origin;

            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else {
                    last_seq2_off[d as usize][k as usize].match_off = INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;

            index = find_first_mismatch_greedy(
                seq1,
                seq2,
                len1,
                len2,
                seq1_index,
                seq2_index,
                reverse,
                rem,
                fence_hit,
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

            last_seq2_off[d as usize][k as usize].match_off = seq2_index;

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

        if d > max_penalty {
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
            last_seq2_off[0][diag_origin as usize].match_off,
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

    let mut max_len = q_avail.max(s_avail);
    if max_len <= 0 {
        return None;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:325-330
    let mut max_dist =
        (GREEDY_MAX_COST as i32).min(max_len / GREEDY_MAX_COST_FRACTION as i32 + 1);

    let mut q_ext_r = 0;
    let mut s_ext_r = 0;
    let mut q_ext_l = 0;
    let mut s_ext_l = 0;

    let GreedyAlignScratch {
        affine_mem,
        fwd_prelim_tback,
        rev_prelim_tback,
    } = scratch;
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
            fwd_prelim_tback.as_deref_mut(),
            rem_forward,
            &mut fence_hit,
            &mut fwd_start_point,
        );
        if fence_hit {
            return None;
        }
        if score >= 0 {
            break;
        }
        max_dist *= 2;
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
            rev_prelim_tback.as_deref_mut(),
            rem_reverse,
            &mut fence_hit,
            &mut rev_start_point,
        );
        if fence_hit {
            return None;
        }
        if score_left >= 0 {
            score += score_left;
            break;
        }
        max_dist *= 2;
    }

    if gap_open == 0 && gap_extend == 0 {
        score = (q_ext_r + s_ext_r + q_ext_l + s_ext_l) * reward / 2
            - score * (reward - penalty);
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
) -> Option<(usize, usize, usize, usize, i32, usize, usize, usize, usize, Vec<GapEditOp>)> {
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

    let (matches, mismatches, gap_opens, gap_letters) = stats_from_edit_ops(
        query,
        subject,
        core.q_start as usize,
        core.s_start as usize,
        &edit_ops,
    );

    Some((
        core.q_start as usize,
        core.q_end as usize,
        core.s_start as usize,
        core.s_end as usize,
        core.score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
        edit_ops,
    ))
}
