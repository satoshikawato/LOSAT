//! TBLASTN gapped uneven-gap linking before preliminary E-value reap.

use crate::utils::glibc_fma_libm::GlibcLibm;
use std::cmp::Ordering;

use anyhow::{ensure, Result};

use super::search_gapped::GappedHsp;
use super::stage_d_stats::{LocalContextLength, LocalLinkParameters};
use crate::stats::spouge::{blast_spouge_stoe, BlastGumbelBlk};
use crate::stats::sum_statistics::{gap_decay_divisor, uneven_gap_sum_e};
use crate::stats::tables::KarlinParams;

// The candidates of a head chosen by subject position.
mod x_link_window;
// Each query context linked on its own (LOSAT_X_TBNLINKCTX).
mod x_link_ctx;
pub use x_link_ctx::print_summary as x_link_ctx_print_summary;
pub use x_link_window::print_summary as x_link_window_print_summary;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-143;
// c++/src/algo/blast/core/link_hsps.c:1765-1810:
// Initializes each HSP's linked count to 1, computes individual E-values via
// `Blast_HSPListGetEvalues`, and links uneven-gap sets via `s_BlastUnevenGapLinkHSPs`.
#[derive(Clone, Copy, Debug)]
pub(super) struct LinkedHsp {
    pub context: usize,
    pub hsp: GappedHsp,
    pub num: i32,
    pub evalue: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1098-1103,1765-1810:
    // The Rust input position identifies the retained HSP value.
    pub source_index: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1810:
#[derive(Clone, Debug)]
pub(super) struct LinkedHspList {
    pub hsps: Vec<LinkedHsp>,
    pub best_evalue: f64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1090-1100,1469-1492:
// Doubly-linked chain node links (`prev` and `next` of `BlastLinkedHSPSet`).
#[derive(Clone, Copy)]
struct WorkHsp {
    value: LinkedHsp,
    prev: Option<usize>,
    next: Option<usize>,
    sum_score: f64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1356:
pub(super) fn score_compare(a: &GappedHsp, b: &GappedHsp) -> Ordering {
    b.score
        .cmp(&a.score)
        .then(a.s_start.cmp(&b.s_start))
        .then(b.s_end.cmp(&a.s_end))
        .then(a.q_start.cmp(&b.q_start))
        .then(b.q_end.cmp(&a.q_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1567-1594:
// Tracks the running prefix index and maximum query end across each context.
fn query_end_prefix_index(work: &[WorkHsp], offset_order: &[usize]) -> Vec<usize> {
    let mut index = vec![0; offset_order.len()];
    let mut current = 0;
    let mut current_end = work[offset_order[0]].value.hsp.q_end;
    for i in 1..offset_order.len() {
        let candidate = offset_order[i];
        let current_id = offset_order[current];
        if work[candidate].value.context > work[current_id].value.context
            || work[candidate].value.hsp.q_end > current_end
        {
            current = i;
            current_end = work[candidate].value.hsp.q_end;
        }
        index[i] = current;
    }
    index
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1234-1263:
fn first_start_at_or_after(
    work: &[WorkHsp],
    order: &[usize],
    context: usize,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, order.len());
    while begin < end {
        let middle = (begin + end) / 2;
        let h = &work[order[middle]].value;
        if h.context < context || (h.context == context && h.hsp.q_start < offset) {
            begin = middle + 1;
        } else {
            end = middle;
        }
    }
    end
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1282-1306:
fn first_end_at_or_after(
    work: &[WorkHsp],
    order: &[usize],
    end_index: &[usize],
    context: usize,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, order.len());
    #[expect(
        clippy::if_same_then_else,
        reason = "NCBI's binary search sets end = left_index in two branches (link_hsps.c:1291-1299, s_HSPOffsetEndBinarySearch)"
    )]
    while begin < end {
        let right = (begin + end) / 2;
        let left = end_index[right];
        let right_context = work[order[right]].value.context;
        if right_context < context {
            begin = right + 1;
        } else if right_context > context {
            end = left;
        } else if work[order[left]].value.hsp.q_end >= offset {
            end = left;
        } else {
            begin = right + 1;
        }
    }
    end
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1310-1354:
fn chain_head(work: &[WorkHsp], mut id: usize) -> usize {
    while let Some(previous) = work[id].prev {
        id = previous;
    }
    id
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1310-1354:
fn merged_chain(work: &[WorkHsp], first: usize, second: usize) -> Vec<usize> {
    let (mut a, mut b) = (
        Some(chain_head(work, first)),
        Some(chain_head(work, second)),
    );
    let mut merged = Vec::new();
    while a.is_some() || b.is_some() {
        if b.is_none()
            || a.is_some_and(|id| work[id].value.hsp.q_start < work[b.unwrap()].value.hsp.q_start)
        {
            let id = a.unwrap();
            merged.push(id);
            a = work[id].next;
        } else {
            let id = b.unwrap();
            merged.push(id);
            b = work[id].next;
        }
    }
    merged
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1414-1491:
fn sets_admissible(
    work: &[WorkHsp],
    first: usize,
    second: usize,
    params: &LocalLinkParameters,
) -> bool {
    if work[first].prev.is_some() {
        return false;
    }
    let second = chain_head(work, second);
    if first == second {
        return false;
    }
    if work[first].value.context != work[second].value.context
        || work[first].value.hsp.frame.signum() != work[second].value.hsp.frame.signum()
    {
        return false;
    }
    let merged = merged_chain(work, first, second);
    for pair in merged.windows(2) {
        let left = work[pair[0]].value.hsp;
        let right = work[pair[1]].value.hsp;
        if left.q_end < right.q_start - params.gap_size
            || left.q_start >= right.q_start
            || left.q_end > right.q_start + params.overlap_size
            || left.q_end >= right.q_end
            || left.s_end > right.s_start + params.overlap_size
            || left.s_end < right.s_start - params.longest_intron
            || left.s_start >= right.s_start
            || left.s_end >= right.s_end
        {
            return false;
        }
    }
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1117-1155:
fn linked_sum_evalue(
    work: &[WorkHsp],
    head: usize,
    candidate: usize,
    query_lengths: &[i32],
    lengths: &[LocalContextLength],
    subject_nt_length: i32,
    params: &LocalLinkParameters,
) -> (f64, f64) {
    let context = work[head].value.context;
    let adjustment = lengths[context].length_adjustment as i32;
    let query_eff_length = (query_lengths[context] - adjustment).max(1);
    let subject_eff_length = (subject_nt_length / 3 - adjustment).max(1);
    let num = (work[head].value.num + work[candidate].value.num) as i16;
    let sum_score = work[candidate].sum_score + work[head].sum_score;
    let evalue = uneven_gap_sum_e(
        params.overlap_size + params.gap_size + 1,
        params.overlap_size + params.longest_intron + 1,
        num,
        sum_score,
        query_eff_length,
        subject_eff_length,
        lengths[context].eff_searchsp,
        gap_decay_divisor(params.gap_decay_rate, num as usize),
    );
    (evalue, sum_score)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1360-1404:
fn combine_sets(
    work: &mut [WorkHsp],
    first: usize,
    second: usize,
    sum_score: f64,
    evalue: f64,
) -> usize {
    let merged = merged_chain(work, first, second);
    let num = merged.len() as i32;
    for (position, &id) in merged.iter().enumerate() {
        work[id].prev = position.checked_sub(1).map(|i| merged[i]);
        work[id].next = merged.get(position + 1).copied();
        work[id].sum_score = sum_score;
        work[id].value.num = num;
        work[id].value.evalue = evalue;
    }
    merged[0]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1602-1757,1765-1810;
// c++/src/algo/blast/core/blast_hits.c:1811-1926:
#[allow(dead_code)] // The public TBLASTN path stays gated until C/D/E pass.
pub(super) fn link_preliminary_hsps(
    input: &[(usize, GappedHsp)],
    query_lengths: &[i32],
    lengths: &[LocalContextLength],
    subject_nt_length: i32,
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    link: &LocalLinkParameters,
) -> Result<LinkedHspList> {
    link_preliminary_hsps_on(
        None,
        input,
        query_lengths,
        lengths,
        subject_nt_length,
        gapped_params,
        gumbel,
        link,
    )
}

/// `link_preliminary_hsps` for a subject whose units ran on `pool` (if any): with
/// `LOSAT_X_TBNLINKCTX` the query contexts are linked on that pool side by side.
#[allow(clippy::too_many_arguments)]
pub(super) fn link_preliminary_hsps_on(
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    input: &[(usize, GappedHsp)],
    query_lengths: &[i32],
    lengths: &[LocalContextLength],
    subject_nt_length: i32,
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    link: &LocalLinkParameters,
) -> Result<LinkedHspList> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1774-1777:
    // Preserve the Stage C-to-D call on an allocated empty subject list.
    if input.is_empty() {
        return Ok(LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        });
    }
    ensure!(
        link.longest_intron > 0,
        "TBLASTN uneven-gap link parameters required"
    );
    ensure!(
        query_lengths.len() == lengths.len() && lengths.len() == gapped_params.len(),
        "one statistical context per query"
    );
    let link_values = |input: &[(usize, GappedHsp)]| {
        x_link_values(
            input,
            query_lengths,
            lengths,
            subject_nt_length,
            gapped_params,
            gumbel,
            link,
        )
    };
    // No NCBI counterpart: LOSAT_X_TBNLINKCTX links each query context's HSPs on their own
    // (`stage_d_linking/x_link_ctx.rs`); a chain never leaves its context, so the values are the
    // same.
    let hsps = if x_link_ctx::x_link_ctx_on() {
        x_link_ctx::x_link_by_context(pool, input, link_values)
    } else {
        link_values(input)
    };
    Ok(x_link_list(hsps))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1779-1782
/// The linking of an ungapped search with a `longest_intron` of 0 (`-max_intron_length` 4 or
/// less): `s_BlastEvenGapLinkHSPs` (`core::link_hsps::even_gap_link_hsps`) of the HSPs of `input`
/// on their subject frames, with the ungapped Karlin blocks `kbp` (`sbp->kbp`), the contexts'
/// lengths and the subject's length in nucleotides (`EvenGapLinkParams::subject_translated`), then
/// sorted by score with the best e-value (link_hsps.c:1801-1810).
pub(super) fn link_even_gap_translated(
    input: &[(usize, GappedHsp)],
    query_lengths: &[i32],
    lengths: &[LocalContextLength],
    subject_nt_length: i32,
    kbp: &[KarlinParams],
    params: &crate::core::link_hsps::EvenGapLinkParams,
) -> Result<LinkedHspList> {
    use crate::core::link_hsps::{even_gap_link_hsps, LinkContext, LinkHsp};
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1774-1777
    // An empty list returns at once (as in `link_preliminary_hsps_on`).
    if input.is_empty() {
        return Ok(LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        });
    }
    let contexts: Vec<LinkContext> = (0..query_lengths.len())
        .map(|context| {
            Ok(LinkContext {
                lambda: kbp[context].lambda,
                log_k: kbp[context].k.glibc_ln(),
                query_length: query_lengths[context],
                length_adjustment: i32::try_from(lengths[context].length_adjustment)?,
                eff_searchsp: lengths[context].eff_searchsp,
            })
        })
        .collect::<Result<_>>()?;
    let hsps: Vec<LinkHsp> = input
        .iter()
        .map(|&(context, hsp)| LinkHsp {
            context,
            subject_frame: hsp.frame,
            query_offset: hsp.q_start,
            query_end: hsp.q_end,
            subject_offset: hsp.s_start,
            subject_end: hsp.s_end,
            score: hsp.score,
        })
        .collect();
    let mut linked = Vec::with_capacity(input.len());
    even_gap_link_hsps(
        &hsps,
        &|context| contexts[context],
        subject_nt_length,
        params,
        &mut |l| {
            let (context, hsp) = input[l.index];
            linked.push(LinkedHsp {
                context,
                hsp,
                num: l.num,
                evalue: l.evalue,
                source_index: l.index,
            })
        },
    );
    Ok(x_link_list(linked))
}

/// The HSPs of `input` linked (`s_BlastUnevenGapLinkHSPs`), with their sum e-values and chain
/// sizes, in the order of `input` (`source_index` = the position in `input`).
fn x_link_values(
    input: &[(usize, GappedHsp)],
    query_lengths: &[i32],
    lengths: &[LocalContextLength],
    subject_nt_length: i32,
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    link: &LocalLinkParameters,
) -> Vec<LinkedHsp> {
    let mut work: Vec<WorkHsp> = input
        .iter()
        .enumerate()
        .map(|(source_index, &(context, hsp))| {
            let kbp = &gapped_params[context];
            let evalue = match gumbel {
                Some(gumbel) => blast_spouge_stoe(
                    hsp.score,
                    kbp,
                    gumbel,
                    query_lengths[context],
                    subject_nt_length / 3,
                ),
                None => {
                    lengths[context].eff_searchsp as f64
                        * (-kbp.lambda * hsp.score as f64 + kbp.k.glibc_ln()).glibc_exp()
                }
            } / gap_decay_divisor(link.gap_decay_rate, 1);
            WorkHsp {
                value: LinkedHsp {
                    context,
                    hsp,
                    num: 1,
                    evalue,
                    source_index,
                },
                prev: None,
                next: None,
                sum_score: kbp.lambda * hsp.score as f64 - kbp.k.glibc_ln(),
            }
        })
        .collect();
    if work.len() > 1 {
        let mut score_order: Vec<usize> = (0..work.len()).collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1217-1243:
        score_order.sort_by(|&a, &b| {
            if work[a].sum_score < work[b].sum_score {
                Ordering::Greater
            } else if work[a].sum_score > work[b].sum_score {
                Ordering::Less
            } else {
                score_compare(&work[a].value.hsp, &work[b].value.hsp)
            }
        });
        let mut offset_order: Vec<usize> = (0..work.len()).collect();
        offset_order.sort_by(|&a, &b| {
            work[a]
                .value
                .context
                .cmp(&work[b].value.context)
                .then(work[a].value.hsp.q_start.cmp(&work[b].value.hsp.q_start))
                .then(work[a].value.hsp.s_start.cmp(&work[b].value.hsp.s_start))
        });
        let end_index = query_end_prefix_index(&work, &offset_order);
        // The lookup by subject position.
        let mut x_window = Some(x_link_window::XLinkWindow::new(&work, &offset_order, link));
        #[cfg(test)]
        if x_link_window::x_whole_window_forced() {
            x_window = None;
        }
        let mut head = None;
        let mut index = 0;
        while index < score_order.len() {
            if head.is_none() {
                while index < score_order.len() {
                    let id = score_order[index];
                    if work[id].prev.is_none() && work[id].next.is_none() {
                        break;
                    }
                    index += 1;
                }
                if index == score_order.len() {
                    break;
                }
                head = Some(score_order[index]);
            }
            let current = head.unwrap();
            let mut tail = current;
            while let Some(next) = work[tail].next {
                tail = next;
            }
            let left_offset = work[current].value.hsp.q_start - link.gap_size;
            let left = first_end_at_or_after(
                &work,
                &offset_order,
                &end_index,
                work[current].value.context,
                left_offset,
            );
            let right = first_start_at_or_after(
                &work,
                &offset_order,
                work[tail].value.context,
                work[tail].value.hsp.q_end + link.gap_size,
            );
            let mut best_evalue = work[current].value.evalue;
            let mut best = None;
            let mut best_sum_score = 0.0;
            // The candidates of the window that can be
            // admissible, in the same order.
            let candidates = match x_window.as_mut() {
                Some(window) => window.candidates(&work, &offset_order, current, left, right),
                // NCBI's whole window: the reference of the unit test only.
                None => &offset_order[left..right],
            };
            for &candidate in candidates {
                if work[candidate]
                    .prev
                    .is_some_and(|prev| work[prev].value.hsp.q_end >= left_offset)
                {
                    continue;
                }
                if sets_admissible(&work, current, candidate, link) {
                    let (evalue, sum_score) = linked_sum_evalue(
                        &work,
                        current,
                        candidate,
                        query_lengths,
                        lengths,
                        subject_nt_length,
                        link,
                    );
                    if evalue < best_evalue.min(work[candidate].value.evalue) {
                        best = Some(candidate);
                        best_evalue = evalue;
                        best_sum_score = sum_score;
                    }
                }
            }
            if let Some(candidate) = best {
                let combined =
                    combine_sets(&mut work, current, candidate, best_sum_score, best_evalue);
                if let Some(window) = x_window.as_mut() {
                    window.linked(&work, combined);
                }
                head = Some(combined);
            } else {
                head = None;
                index += 1;
            }
        }
    }
    work.into_iter().map(|item| item.value).collect()
}

/// The list of linked HSPs `hsps` (in the order of the input): sorted by score
/// (`Blast_HSPListSortByScore`, stable), with the best e-value.
fn x_link_list(mut hsps: Vec<LinkedHsp>) -> LinkedHspList {
    hsps.sort_by(|a, b| score_compare(&a.hsp, &b.hsp));
    let best_evalue = hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .reduce(f64::min)
        .unwrap_or(f64::MAX);
    LinkedHspList { hsps, best_evalue }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1385-1404,3071-3115,3404-3417;
// c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319
// The reverse-sort array is consumed from its end, yielding this comparator's
// ascending best-E-value, descending score, descending OID order.
#[allow(dead_code)] // Used by the internal TBLASTN composition result pipeline.
pub(super) fn compare_preliminary_lists_for_kappa(
    oid_a: i32,
    a: &LinkedHspList,
    oid_b: i32,
    b: &LinkedHspList,
) -> Ordering {
    match (a.hsps.is_empty(), b.hsps.is_empty()) {
        (true, true) => return Ordering::Equal,
        (true, false) => return Ordering::Greater,
        (false, true) => return Ordering::Less,
        (false, false) => {}
    }
    let evalue = if a.best_evalue < 1.0e-180 && b.best_evalue < 1.0e-180 {
        Ordering::Equal
    } else if a.best_evalue < b.best_evalue {
        Ordering::Less
    } else if a.best_evalue > b.best_evalue {
        Ordering::Greater
    } else {
        Ordering::Equal
    };
    evalue
        .then_with(|| b.hsps[0].hsp.score.cmp(&a.hsps[0].hsp.score))
        .then_with(|| oid_b.cmp(&oid_a))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-676 and
// c++/src/algo/blast/core/blast_hits.c:1976-2010:
// Reaps HSPs exceeding `prelim_evalue` in the preliminary pass or `expect_value`
// in the post-Kappa pass, compacting surviving entries in place. A NaN E-value (a PSSM whose
// gapped block is not valid) is not above the limit, and its HSP stays. The preliminary stage of
// each subject reaps with it (`stage_d_pipeline`).
pub(super) fn reap_by_evalue(list: &mut LinkedHspList, prelim_evalue: f64) {
    list.hsps
        .retain(|hsp| hsp.evalue <= prelim_evalue || hsp.evalue.is_nan());
}

#[cfg(test)]
#[path = "stage_d_linking_tests.rs"]
mod tests;
