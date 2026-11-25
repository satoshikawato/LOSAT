//! Sparse 64-bit affine greedy alignment for large gap costs (`gap_open > 32767` or
//! `gap_extend > 32767`, or any dense setup whose 32-bit distance arithmetic would overflow).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:148-258,755-1249
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:171-267
//!
//! NCBI's dense affine greedy alignment stores arrays indexed by every integer distance up to
//! `scaled_max_dist + max_penalty`, which overflows signed `Int4` or exhausts memory when gap
//! costs are large (see `PD-LOSAT-NCBI-DEFECTS`). Storing only non-empty distances in strictly
//! increasing order and stepping with three monotone cursors reproduces NCBI's exact DP
//! recurrence and tie-breaking without allocating dense distance tables or changing the hot
//! dense path used for ordinary gap costs.

use super::{
    affine_first_mismatch, debug_greedy_traceback_enabled, finish_greedy_gapped_core, gdb3,
    gdb3_i64, prior_at, prior_row_base, AffineRow, ContextQuery, GapAlignOpType,
    GapPrelimEditBlock, GreedyAlignScratch, GreedyGappedCore, GreedyOffset, GreedySeed, SeqPointer,
    INVALID_DIAG, INVALID_GREEDY_OFFSET, INVALID_OFFSET,
};

/// Sentinel for an absent parent node in `AffineNode`.
const NO_AFFINE_NODE: u32 = u32::MAX;

/// Cap on cached row offset buffers retained across extensions so traceback calls do not
/// accumulate unbounded heap memory in the scratch pool.
const MAX_POOLED_ROWS: usize = 64;

/// Metadata for one non-empty distance in `blast_affine_greedy_align_sparse`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:920-955,1139-1180
/// Every empty distance leaves the best score unchanged (its extent is 0, so its score is
/// strictly negative), so storing the running maximum at non-empty distances reproduces NCBI's
/// step function exactly.
#[derive(Clone, Copy)]
struct AffineNode {
    dist: i64,
    lower: i32,
    upper: i32,
    max_score: i64,
    parent_op: u32,
    parent_ext: u32,
    parent_open: u32,
}

/// Reusable scratch space for sparse 64-bit affine greedy extensions.
pub(super) struct SparseAffineMem {
    nodes: Vec<AffineNode>,
    rows: Vec<AffineRow>,
    free_offsets: Vec<Vec<GreedyOffset>>,
}

impl SparseAffineMem {
    pub(super) fn new() -> Self {
        Self {
            nodes: Vec::new(),
            rows: Vec::new(),
            free_offsets: Vec::new(),
        }
    }

    fn reset(&mut self) {
        while let Some(row) = self.rows.pop() {
            self.recycle_row(row);
        }
        self.nodes.clear();
    }

    fn reclaim_before(&mut self, dead_cursor: &mut usize, min_live: usize) {
        while *dead_cursor < min_live {
            let old = std::mem::take(&mut self.rows[*dead_cursor]);
            self.recycle_row(old);
            *dead_cursor += 1;
        }
    }

    fn take_row(&mut self, lower: i32, upper: i32) -> AffineRow {
        let mut offsets = self.free_offsets.pop().unwrap_or_default();
        offsets.clear();
        offsets.resize((upper - lower + 1).max(0) as usize, INVALID_GREEDY_OFFSET);
        AffineRow { lower, offsets }
    }

    fn recycle_row(&mut self, row: AffineRow) {
        if row.offsets.capacity() > 0 && self.free_offsets.len() < MAX_POOLED_ROWS {
            self.free_offsets.push(row.offsets);
        }
    }
}

/// Sparse affine traceback helper from match state.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:148-178
fn get_next_sparse_affine_tback_from_match(
    nodes: &[AffineNode],
    rows: &[AffineRow],
    cur_idx: &mut usize,
    diag: i32,
    seq2_index: &mut i32,
) -> GapAlignOpType {
    let cur_cell = rows[*cur_idx].at(diag);
    let parent_op = nodes[*cur_idx].parent_op;
    if parent_op != NO_AFFINE_NODE {
        let p_idx = parent_op as usize;
        let p = &nodes[p_idx];
        if diag >= p.lower && diag <= p.upper {
            let new_seq2_index = rows[p_idx].at(diag).match_off;
            if new_seq2_index >= cur_cell.insert_off.max(cur_cell.delete_off) {
                *cur_idx = p_idx;
                *seq2_index = new_seq2_index;
                return GapAlignOpType::Sub;
            }
        }
    }

    if cur_cell.insert_off > cur_cell.delete_off {
        *seq2_index = cur_cell.insert_off;
        GapAlignOpType::Ins
    } else {
        *seq2_index = cur_cell.delete_off;
        GapAlignOpType::Del
    }
}

/// Sparse affine traceback helper from insertion/deletion state.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:199-258
fn get_next_sparse_affine_tback_from_indel(
    nodes: &[AffineNode],
    rows: &[AffineRow],
    cur_idx: &mut usize,
    diag: i32,
    state: GapAlignOpType,
) -> GapAlignOpType {
    let new_diag = if state == GapAlignOpType::Ins {
        diag - 1
    } else {
        diag + 1
    };

    let mut new_seq2_index = INVALID_OFFSET;
    let parent_ext = nodes[*cur_idx].parent_ext;
    if parent_ext != NO_AFFINE_NODE {
        let p_idx = parent_ext as usize;
        let p = &nodes[p_idx];
        if new_diag >= p.lower && new_diag <= p.upper {
            let cell = rows[p_idx].at(new_diag);
            new_seq2_index = if state == GapAlignOpType::Ins {
                cell.insert_off
            } else {
                cell.delete_off
            };
        }
    }

    let parent_open = nodes[*cur_idx].parent_open;
    if parent_open != NO_AFFINE_NODE {
        let p_idx = parent_open as usize;
        let p = &nodes[p_idx];
        if new_diag >= p.lower
            && new_diag <= p.upper
            && new_seq2_index < rows[p_idx].at(new_diag).match_off
        {
            *cur_idx = p_idx;
            return GapAlignOpType::Sub;
        }
    }

    debug_assert!(parent_ext != NO_AFFINE_NODE);
    *cur_idx = parent_ext as usize;
    state
}

/// Returns true when the dense 32-bit affine greedy table sizes would overflow `i32`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:802-873
#[cold]
#[inline(never)]
pub(super) fn dense_affine_overflows(
    max_dist: i32,
    mut xdrop_threshold: i32,
    mut match_score: i32,
    mut mismatch_score: i32,
    mut in_gap_open: i32,
    mut in_gap_extend: i32,
) -> bool {
    if match_score % 2 == 1 {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold = xdrop_threshold.wrapping_mul(2);
        in_gap_open *= 2;
        in_gap_extend *= 2;
    }
    let match_score_half = match_score / 2;
    let mut op_cost = match_score + mismatch_score;
    let mut gap_open = in_gap_open;
    let mut gap_extend = in_gap_extend + match_score_half;
    let score_common_factor = gdb3(&mut op_cost, &mut gap_open, &mut gap_extend);
    let max_penalty = op_cost.max(gap_open + gap_extend);
    let xdrop_offset = (xdrop_threshold.wrapping_add(match_score_half)) / score_common_factor + 1;
    max_dist
        .checked_mul(gap_extend)
        .and_then(|smd| {
            smd.checked_add(xdrop_offset)?.checked_add(2)?;
            smd.checked_add(max_penalty)?.checked_add(2)?;
            Some(())
        })
        .is_none()
}

/// Sparse 64-bit two-direction affine greedy gapped alignment, invoked out-of-line from
/// `greedy_gapped_alignment_internal` when gap costs exceed `MAX_DENSE_GREEDY_GAP_COST` or when
/// the dense 32-bit table sizing would overflow.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2762-2885
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn greedy_gapped_alignment_sparse(
    query: ContextQuery<'_>,
    subject: SeqPointer<'_>,
    q_off: i32,
    s_off: i32,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut GreedyAlignScratch,
    do_traceback: bool,
) -> Option<GreedyGappedCore> {
    let (query_right, q_len_right) = query.right_of(q_off);
    let (query_left, q_len_left) = query.left_of(q_off);
    let (subject_right, s_len_right) = subject.right_of(s_off);
    let (subject_left, s_len_left) = subject.left_of(s_off);

    let mut q_ext_r = 0;
    let mut s_ext_r = 0;
    let mut q_ext_l = 0;
    let mut s_ext_l = 0;

    if do_traceback {
        scratch.fwd_prelim_tback.reset();
        scratch.rev_prelim_tback.reset();
    }

    let mut fwd_start_point = GreedySeed::default();
    let mut rev_start_point = GreedySeed::default();
    let mut fence_hit = false;
    let trace_dir = debug_greedy_traceback_enabled(q_off, s_off);

    let score_right = blast_affine_greedy_align_sparse::<false>(
        query_right,
        q_len_right,
        subject_right,
        s_len_right,
        x_drop,
        reward,
        -penalty,
        gap_open,
        gap_extend,
        &mut q_ext_r,
        &mut s_ext_r,
        &mut scratch.sparse_mem,
        do_traceback.then_some(&mut scratch.fwd_prelim_tback),
        &mut fence_hit,
        &mut fwd_start_point,
    );
    if fence_hit {
        super::super::x_fence::x_fence_met();
        return None;
    }
    if trace_dir {
        eprintln!(
            "[GREEDY_DIRECTION] right q_ext={} s_ext={} score={}",
            q_ext_r, s_ext_r, score_right
        );
    }

    let score_left = blast_affine_greedy_align_sparse::<true>(
        query_left,
        q_len_left,
        subject_left,
        s_len_left,
        x_drop,
        reward,
        -penalty,
        gap_open,
        gap_extend,
        &mut q_ext_l,
        &mut s_ext_l,
        &mut scratch.sparse_mem,
        do_traceback.then_some(&mut scratch.rev_prelim_tback),
        &mut fence_hit,
        &mut rev_start_point,
    );
    if fence_hit {
        super::super::x_fence::x_fence_met();
        return None;
    }
    if trace_dir {
        eprintln!(
            "[GREEDY_DIRECTION] left q_ext={} s_ext={} score={}",
            q_ext_l, s_ext_l, score_left
        );
    }

    let score = if reward % 2 == 1 {
        ((i64::from(score_right) + i64::from(score_left)) / 2) as i32
    } else {
        score_right + score_left
    };

    Some(finish_greedy_gapped_core(
        query,
        subject,
        q_off,
        s_off,
        q_ext_l,
        s_ext_l,
        q_ext_r,
        s_ext_r,
        score,
        fwd_start_point,
        rev_start_point,
        scratch,
        do_traceback,
    ))
}

/// Sparse 64-bit affine greedy alignment for one direction.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:755-1249
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn blast_affine_greedy_align_sparse<const REVERSE: bool>(
    seq1: &[u8],
    len1: i32,
    seq2: &[u8],
    len2: i32,
    mut xdrop_threshold: i32,
    mut match_score: i32,
    mut mismatch_score: i32,
    in_gap_open: i32,
    in_gap_extend: i32,
    seq1_align_len: &mut i32,
    seq2_align_len: &mut i32,
    sparse_mem: &mut Option<Box<SparseAffineMem>>,
    mut edit_block: Option<&mut GapPrelimEditBlock>,
    fence_hit: &mut bool,
    seed: &mut GreedySeed,
) -> i32 {
    let mut seq1_index: i32;
    let mut seq2_index: i32;
    let mut index: i32;
    let mut k: i32;
    let mut longest_match_run: i32;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:801-812
    // Gap costs are widened to 64 bits so values above 32767 and up to `i32::MAX` do not hit
    // NCBI's signed 32-bit overflow.
    let mut gap_open_in = i64::from(in_gap_open);
    let mut gap_extend_in = i64::from(in_gap_extend);
    if match_score % 2 == 1 {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold = xdrop_threshold.wrapping_mul(2);
        gap_open_in *= 2;
        gap_extend_in *= 2;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:837-873
    let match_score_half = i64::from(match_score / 2);
    let mut op_cost = i64::from(match_score) + i64::from(mismatch_score);
    let mut gap_open = gap_open_in;
    let mut gap_extend = gap_extend_in + match_score_half;
    let score_common_factor = gdb3_i64(&mut op_cost, &mut gap_open, &mut gap_extend);
    let gap_open_extend = gap_open + gap_extend;
    debug_assert!(gap_extend >= 1 && gap_open_extend >= 1 && op_cost >= 1);

    let xdrop_offset =
        i64::from(xdrop_threshold.wrapping_add(match_score / 2)) / score_common_factor + 1;

    index = affine_first_mismatch::<REVERSE>(seq1, seq2, len1, len2, 0, 0, fence_hit);
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:920-1201
    let sparse = sparse_mem.get_or_insert_with(|| Box::new(SparseAffineMem::new()));
    sparse.reset();
    let keep_all_rows = edit_block.is_some();

    let mut best_score: i64 = i64::from(seq1_index) * i64::from(match_score);
    let mut best_node: usize = 0;
    let mut best_diag: i32 = 0;

    let mut row0 = sparse.take_row(0, 0);
    row0.offsets[0] = GreedyOffset {
        insert_off: INVALID_OFFSET,
        match_off: seq1_index,
        delete_off: INVALID_OFFSET,
    };
    sparse.nodes.push(AffineNode {
        dist: 0,
        lower: 0,
        upper: 0,
        max_score: best_score,
        parent_op: NO_AFFINE_NODE,
        parent_ext: NO_AFFINE_NODE,
        parent_open: NO_AFFINE_NODE,
    });
    sparse.rows.push(row0);

    let mut end1_diag: Option<i32> = None;
    let mut end2_diag: Option<i32> = None;
    let mut p_op: usize = 0;
    let mut p_ext: usize = 0;
    let mut p_open: usize = 0;
    let mut dead_cursor: usize = 0;
    let mut xdrop_ptr: usize = 0;
    let empty_row = AffineRow::default();

    loop {
        let cand_op = sparse.nodes.get(p_op).map(|n| n.dist + op_cost);
        let cand_ext = sparse.nodes.get(p_ext).map(|n| n.dist + gap_extend);
        let cand_open = sparse.nodes.get(p_open).map(|n| n.dist + gap_open_extend);
        let Some(d) = [cand_op, cand_ext, cand_open].into_iter().flatten().min() else {
            break;
        };

        let idx_op = if cand_op == Some(d) {
            let idx = p_op;
            p_op += 1;
            Some(idx)
        } else {
            None
        };
        let idx_ext = if cand_ext == Some(d) {
            let idx = p_ext;
            p_ext += 1;
            Some(idx)
        } else {
            None
        };
        let idx_open = if cand_open == Some(d) {
            let idx = p_open;
            p_open += 1;
            Some(idx)
        } else {
            None
        };

        let (open_lo, open_hi) = idx_open.map_or((INVALID_DIAG, -INVALID_DIAG), |i| {
            let n = &sparse.nodes[i];
            (n.lower, n.upper)
        });
        let (ext_lo, ext_hi) = idx_ext.map_or((INVALID_DIAG, -INVALID_DIAG), |i| {
            let n = &sparse.nodes[i];
            (n.lower, n.upper)
        });
        let (match_lo, match_hi) = idx_op.map_or((INVALID_DIAG, -INVALID_DIAG), |i| {
            let n = &sparse.nodes[i];
            (n.lower, n.upper)
        });

        let mut curr_diag_lower = open_lo.min(ext_lo) - 1;
        curr_diag_lower = curr_diag_lower.min(match_lo);
        if let Some(e2) = end2_diag {
            curr_diag_lower = curr_diag_lower.max(e2);
        }

        let mut curr_diag_upper = open_hi.max(ext_hi) + 1;
        curr_diag_upper = curr_diag_upper.max(match_hi);
        if let Some(e1) = end1_diag {
            curr_diag_upper = curr_diag_upper.min(e1);
        }

        if curr_diag_lower > curr_diag_upper {
            if !keep_all_rows {
                sparse.reclaim_before(&mut dead_cursor, p_op.min(p_ext).min(p_open));
            }
            continue;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:979-983
        let prev_xdrop_max = if d < xdrop_offset {
            0i64
        } else {
            let target = d - xdrop_offset;
            while xdrop_ptr + 1 < sparse.nodes.len() && sparse.nodes[xdrop_ptr + 1].dist <= target {
                xdrop_ptr += 1;
            }
            sparse.nodes[xdrop_ptr].max_score
        };
        let xdrop_raw = prev_xdrop_max
            .saturating_add(score_common_factor.saturating_mul(d))
            .saturating_sub(i64::from(xdrop_threshold));
        let xdrop_score = (((xdrop_raw as f64) / (match_score_half as f64)).ceil() as i64).max(0);

        let mut curr_extent: i32 = 0;
        let mut curr_seq2_index: i32 = 0;
        let mut curr_diag: i32 = 0;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;
        let mut cur = sparse.take_row(tmp_diag_lower, tmp_diag_upper);

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:1000-1047
        let row_open = idx_open.map_or(&empty_row, |i| &sparse.rows[i]);
        let row_ext = idx_ext.map_or(&empty_row, |i| &sparse.rows[i]);
        let row_match = idx_op.map_or(&empty_row, |i| &sparse.rows[i]);
        let open_base = prior_row_base(row_open, open_lo, open_hi);
        let ext_base = prior_row_base(row_ext, ext_lo, ext_hi);
        let match_base = prior_row_base(row_match, match_lo, match_hi);
        assert!(
            cur.lower == tmp_diag_lower
                && cur.offsets.len() == (tmp_diag_upper - tmp_diag_lower + 1).max(0) as usize
        );
        let mut fence = false;

        for (k_val, cell) in (tmp_diag_lower..=tmp_diag_upper).zip(cur.offsets.iter_mut()) {
            k = k_val;

            let mut seq2_index_del = INVALID_OFFSET;
            if k + 1 >= open_lo && k < open_hi {
                // SAFETY: `k + 1` lies in `open_lo..=open_hi`, checked for `open_base`.
                seq2_index_del = unsafe { prior_at(open_base, k + 1) }.match_off;
            }
            if k + 1 >= ext_lo && k < ext_hi {
                // SAFETY: `k + 1` lies in `ext_lo..=ext_hi`, checked for `ext_base`.
                let ext_off = unsafe { prior_at(ext_base, k + 1) }.delete_off;
                if ext_off > seq2_index_del {
                    seq2_index_del = ext_off;
                }
            }
            cell.delete_off = if seq2_index_del == INVALID_OFFSET {
                INVALID_OFFSET
            } else {
                seq2_index_del + 1
            };

            let mut seq2_index_ins = INVALID_OFFSET;
            if k > open_lo && k - 1 <= open_hi {
                // SAFETY: `k - 1` lies in `open_lo..=open_hi`, checked for `open_base`.
                seq2_index_ins = unsafe { prior_at(open_base, k - 1) }.match_off;
            }
            if k > ext_lo && k - 1 <= ext_hi {
                // SAFETY: `k - 1` lies in `ext_lo..=ext_hi`, checked for `ext_base`.
                let ext_off = unsafe { prior_at(ext_base, k - 1) }.insert_off;
                if ext_off > seq2_index_ins {
                    seq2_index_ins = ext_off;
                }
            }
            cell.insert_off = seq2_index_ins;

            seq2_index = cell.insert_off.max(cell.delete_off);
            if k >= match_lo && k <= match_hi {
                // SAFETY: `k` lies in `match_lo..=match_hi`, checked for `match_base`.
                seq2_index = seq2_index.max(unsafe { prior_at(match_base, k) }.match_off + 1);
            }

            seq1_index = seq2_index + k;

            if seq2_index < 0 || i64::from(seq1_index) + i64::from(seq2_index) < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else {
                    cell.match_off = INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;

            index = affine_first_mismatch::<REVERSE>(
                seq1, seq2, len1, len2, seq1_index, seq2_index, fence_hit,
            );
            if *fence_hit {
                fence = true;
                break;
            }

            if index > longest_match_run {
                seed.start_q = seq1_index;
                seed.start_s = seq2_index;
                longest_match_run = index;
                seed.match_length = longest_match_run;
            }
            seq1_index += index;
            seq2_index += index;

            cell.match_off = seq2_index;

            if seq1_index + seq2_index > curr_extent {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }

            if seq1_index == len1 {
                curr_diag_upper = k;
                end1_diag = Some(k - 1);
            }
            if seq2_index == len2 {
                curr_diag_lower = k;
                end2_diag = Some(k + 1);
            }
        }
        if fence {
            sparse.recycle_row(cur);
            return -1;
        }

        if curr_diag_lower <= curr_diag_upper {
            let curr_score = i64::from(curr_extent)
                .saturating_mul(match_score_half)
                .saturating_sub(d.saturating_mul(score_common_factor));
            if curr_score > best_score {
                best_score = curr_score;
                best_node = sparse.nodes.len();
                best_diag = curr_diag;
                *seq2_align_len = curr_seq2_index;
                *seq1_align_len = curr_seq2_index + best_diag;
            }
            sparse.nodes.push(AffineNode {
                dist: d,
                lower: curr_diag_lower,
                upper: curr_diag_upper,
                max_score: best_score,
                parent_op: idx_op.map_or(NO_AFFINE_NODE, |i| i as u32),
                parent_ext: idx_ext.map_or(NO_AFFINE_NODE, |i| i as u32),
                parent_open: idx_open.map_or(NO_AFFINE_NODE, |i| i as u32),
            });
            sparse.rows.push(cur);
        } else {
            sparse.recycle_row(cur);
        }

        if !keep_all_rows {
            sparse.reclaim_before(&mut dead_cursor, p_op.min(p_ext).min(p_open));
        }
    }

    if let Some(block) = edit_block.as_mut() {
        let mut cur_idx = best_node;
        seq2_index = *seq2_align_len;
        let mut state = GapAlignOpType::Sub;

        while cur_idx > 0 {
            if state == GapAlignOpType::Sub {
                let mut new_seq2_index = 0;
                state = get_next_sparse_affine_tback_from_match(
                    &sparse.nodes,
                    &sparse.rows,
                    &mut cur_idx,
                    best_diag,
                    &mut new_seq2_index,
                );
                block.add(GapAlignOpType::Sub, seq2_index - new_seq2_index);
                seq2_index = new_seq2_index;
            } else if state == GapAlignOpType::Ins {
                block.add(GapAlignOpType::Ins, 1);
                state = get_next_sparse_affine_tback_from_indel(
                    &sparse.nodes,
                    &sparse.rows,
                    &mut cur_idx,
                    best_diag,
                    GapAlignOpType::Ins,
                );
                best_diag -= 1;
            } else {
                block.add(GapAlignOpType::Del, 1);
                state = get_next_sparse_affine_tback_from_indel(
                    &sparse.nodes,
                    &sparse.rows,
                    &mut cur_idx,
                    best_diag,
                    GapAlignOpType::Del,
                );
                best_diag += 1;
                seq2_index -= 1;
            }
        }

        block.add(GapAlignOpType::Sub, sparse.rows[0].at(0).match_off);
    }

    best_score as i32
}
