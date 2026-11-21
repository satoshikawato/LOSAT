//! BLASTX's linking with sum statistics in internal context coordinates: the uneven-gap linking
//! (`link_uneven`) and the even-gap linking of the ungapped search (`link_even`, `core::link_hsps`).
use super::{
    args::ResolvedOptions,
    parameters::ContextParameters,
    preliminary::{compare_score, PreliminaryHsp},
    query_setup::PreparedQueryBatch,
};
use crate::core::link_hsps::{
    calculate_link_hsp_cutoffs, even_gap_link_hsps, EvenGapLinkParams, LinkContext, LinkCutoffs,
    LinkHsp,
};
use crate::stats::{
    spouge::{blast_spouge_stoe, BlastGumbelBlk},
    sum_statistics::{gap_decay_divisor, uneven_gap_sum_e},
    tables::KarlinParams,
};
use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::{ensure, Result};
use std::cmp::Ordering;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:598-617
#[derive(Clone, Debug)]
pub struct LinkParameters {
    pub gap_size: i32,
    pub overlap_size: i32,
    pub longest_intron: i32,
    pub gap_decay_rate: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:781-815
pub fn link_parameters(options: &ResolvedOptions) -> Option<LinkParameters> {
    if !options.sum_stats {
        return None;
    }
    // The conversion of `-max_intron_length` above has one owner, shared with TBLASTN.
    let intron = crate::core::link_hsps::translated_longest_intron(
        options.max_intron_length,
        options.gapped,
    )?;
    Some(LinkParameters {
        gap_size: 40,
        overlap_size: 9,
        longest_intron: intron,
        gap_decay_rate: if options.gapped {
            crate::core::link_hsps::BLAST_GAP_DECAY_RATE_GAPPED
        } else {
            crate::core::link_hsps::BLAST_GAP_DECAY_RATE
        },
    })
}

/// Whether the search computes the linking cutoffs of each subject (`CalculateLinkHSPCutoffs`,
/// `link_cutoffs`): an ungapped search that links. With `-ungapped`, a search links
/// exactly when `-sum_stats` is on (`link_parameters`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
pub(crate) fn computes_link_cutoffs(options: &ResolvedOptions) -> bool {
    !options.gapped && link_parameters(options).is_some()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1789
#[derive(Clone, Debug)]
pub struct LinkedHsp {
    pub hsp: PreliminaryHsp,
    pub num: i32,
    pub evalue: f64,
    pub source_index: usize,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1802-1810
#[derive(Clone, Debug)]
pub struct LinkedHspList {
    pub hsps: Vec<LinkedHsp>,
    pub best_evalue: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1093-1101
#[derive(Clone)]
struct WorkHsp {
    value: LinkedHsp,
    prev: Option<usize>,
    next: Option<usize>,
    sum_score: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1567-1594
fn query_end_prefix_index(work: &[WorkHsp], offset_order: &[usize]) -> Vec<usize> {
    let mut index = vec![0; offset_order.len()];
    let mut current = 0;
    let mut current_end = work[offset_order[0]].value.hsp.q_end;
    for i in 1..offset_order.len() {
        let candidate = offset_order[i];
        let current_id = offset_order[current];
        if work[candidate].value.hsp.context / 3 > work[current_id].value.hsp.context / 3
            || work[candidate].value.hsp.q_end > current_end
        {
            current = i;
            current_end = work[candidate].value.hsp.q_end;
        }
        index[i] = current;
    }
    index
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1234-1263
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
        if h.hsp.context / 3 < context || (h.hsp.context / 3 == context && h.hsp.q_start < offset) {
            begin = middle + 1;
        } else {
            end = middle;
        }
    }
    end
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1282-1306
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
        let right_context = work[order[right]].value.hsp.context / 3;
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1310-1325
fn chain_head(work: &[WorkHsp], mut id: usize) -> usize {
    while let Some(previous) = work[id].prev {
        id = previous;
    }
    id
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1310-1354
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1416-1492
fn sets_admissible(work: &[WorkHsp], first: usize, second: usize, params: &LinkParameters) -> bool {
    if work[first].prev.is_some() {
        return false;
    }
    let second = chain_head(work, second);
    if first == second {
        return false;
    }
    if work[first].value.hsp.context / 3 != work[second].value.hsp.context / 3 {
        return false;
    }
    let merged = merged_chain(work, first, second);
    for pair in merged.windows(2) {
        let left = &work[pair[0]].value.hsp;
        let right = &work[pair[1]].value.hsp;
        if left.q_end < right.q_start - params.longest_intron
            || left.q_start >= right.q_start
            || left.q_end > right.q_start + params.overlap_size
            || left.q_end >= right.q_end
            || left.s_end > right.s_start + params.overlap_size
            || left.s_end < right.s_start - params.gap_size
            || left.s_start >= right.s_start
            || left.s_end >= right.s_end
        {
            return false;
        }
    }
    true
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1114-1155
fn linked_sum_evalue(
    work: &[WorkHsp],
    head: usize,
    candidate: usize,
    query_lengths: &[i32],
    lengths: &[ContextParameters],
    subject_length: i32,
    params: &LinkParameters,
) -> (f64, f64) {
    let context = work[head].value.hsp.context;
    let adjustment = lengths[context].length_adjustment;
    let query_eff_length = (query_lengths[context] - adjustment).max(1);
    let subject_eff_length = (subject_length - adjustment).max(1);
    let num = (work[head].value.num + work[candidate].value.num) as i16;
    let sum_score = work[candidate].sum_score + work[head].sum_score;
    let evalue = uneven_gap_sum_e(
        params.overlap_size + params.gap_size + 1,
        params.overlap_size + params.longest_intron + 1,
        num,
        sum_score,
        query_eff_length,
        subject_eff_length,
        lengths[context].search_space,
        gap_decay_divisor(params.gap_decay_rate, num as usize),
    );
    (evalue, sum_score)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1361-1404
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
// The public BLASTX path stays gated until C/D/E pass.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1627-1666
pub fn link_uneven(
    input: &[PreliminaryHsp],
    query_lengths: &[i32],
    lengths: &[ContextParameters],
    subject_length: i32,
    gapped_params: &[KarlinParams],
    gumbel: Option<&BlastGumbelBlk>,
    link: &LinkParameters,
) -> Result<LinkedHspList> {
    if input.is_empty() {
        return Ok(LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        });
    }
    ensure!(
        link.longest_intron > 0,
        "BLASTX uneven-gap link parameters required"
    );
    ensure!(
        query_lengths.len() == lengths.len() && lengths.len() == gapped_params.len(),
        "one statistical context per query"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1524-1540
    let mut work: Vec<WorkHsp> = input
        .iter()
        .enumerate()
        .map(|(source_index, hsp)| {
            let context = hsp.context;
            let hsp = hsp.clone();
            let kbp = &gapped_params[context];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1851-1896
            let evalue = if let Some(gumbel) = gumbel {
                blast_spouge_stoe(
                    hsp.score,
                    kbp,
                    gumbel,
                    query_lengths[context],
                    subject_length,
                )
            } else {
                lengths[context].search_space as f64
                    * (-kbp.lambda * hsp.score as f64 + kbp.k.glibc_ln()).glibc_exp()
            } / gap_decay_divisor(link.gap_decay_rate, 1);
            let sum_score = kbp.lambda * hsp.score as f64 - kbp.k.glibc_ln();
            WorkHsp {
                value: LinkedHsp {
                    hsp,
                    num: 1,
                    evalue,
                    source_index,
                },
                prev: None,
                next: None,
                sum_score,
            }
        })
        .collect();
    if work.len() > 1 {
        let mut score_order: Vec<usize> = (0..work.len()).collect();
        score_order.sort_by(|&a, &b| {
            if work[a].sum_score < work[b].sum_score {
                Ordering::Greater
            } else if work[a].sum_score > work[b].sum_score {
                Ordering::Less
            } else {
                compare_score(&work[a].value.hsp, &work[b].value.hsp)
            }
        });
        let mut offset_order: Vec<usize> = (0..work.len()).collect();
        offset_order.sort_by(|&a, &b| {
            work[a]
                .value
                .hsp
                .context
                .div_euclid(3)
                .cmp(&work[b].value.hsp.context.div_euclid(3))
                .then(work[a].value.hsp.q_start.cmp(&work[b].value.hsp.q_start))
                .then(work[a].value.hsp.s_start.cmp(&work[b].value.hsp.s_start))
        });
        let end_index = query_end_prefix_index(&work, &offset_order);
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
            let left_offset = work[current].value.hsp.q_start - link.longest_intron;
            let left = first_end_at_or_after(
                &work,
                &offset_order,
                &end_index,
                work[current].value.hsp.context / 3,
                left_offset,
            );
            let right = first_start_at_or_after(
                &work,
                &offset_order,
                work[tail].value.hsp.context / 3,
                work[tail].value.hsp.q_end + link.longest_intron,
            );
            let mut best_evalue = work[current].value.evalue;
            let mut best = None;
            let mut best_sum_score = 0.0;
            for &candidate in &offset_order[left..right] {
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
                        subject_length,
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
                head = Some(combine_sets(
                    &mut work,
                    current,
                    candidate,
                    best_sum_score,
                    best_evalue,
                ));
            } else {
                head = None;
                index += 1;
            }
        }
    }
    work.sort_by(|a, b| compare_score(&a.value.hsp, &b.value.hsp));
    let hsps: Vec<_> = work.into_iter().map(|item| item.value).collect();
    let best_evalue = hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .reduce(f64::min)
        .unwrap_or(f64::MAX);
    Ok(LinkedHspList { hsps, best_evalue })
}

/// The even-gap linking of BLASTX's ungapped search (`BLAST_LinkHsps` with `longest_intron <= 0`,
/// `core::link_hsps`): the query is translated, so the three frames of a strand link together,
/// and every context reads its ungapped Karlin block (`sbp->kbp`, `gapped_calculation` false).
/// The list is then sorted by score (`Blast_HSPListSortByScore`) and its best e-value is the
/// smallest of its e-values (link_hsps.c:1801-1811).
pub(crate) fn link_even(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    subject_length: i32,
    params: &EvenGapLinkParams,
) -> LinkedHspList {
    let hsps: Vec<LinkHsp> = input
        .iter()
        .map(|h| LinkHsp {
            context: h.context,
            query_offset: h.q_start,
            query_end: h.q_end,
            subject_offset: h.s_start,
            subject_end: h.s_end,
            score: h.score,
        })
        .collect();
    let mut linked: Vec<LinkedHsp> = Vec::with_capacity(input.len());
    even_gap_link_hsps(
        &hsps,
        &|context| {
            let p = &parameters[context];
            LinkContext {
                lambda: p.ungapped.lambda,
                log_k: p.ungapped.k.glibc_ln(),
                query_length: batch.contexts[context].length as i32,
                length_adjustment: p.length_adjustment,
                eff_searchsp: p.search_space,
            }
        },
        subject_length,
        params,
        &mut |l| {
            linked.push(LinkedHsp {
                hsp: input[l.index].clone(),
                num: l.num,
                evalue: l.evalue,
                source_index: l.index,
            })
        },
    );
    let mut hsps = linked;
    hsps.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
    // `best_evalue` from the first HSP, lowered by each smaller one (the empty list keeps 0).
    let best_evalue = hsps.split_first().map_or(0.0, |(first, rest)| {
        rest.iter().fold(
            first.evalue,
            |best, h| {
                if h.evalue < best {
                    h.evalue
                } else {
                    best
                }
            },
        )
    });
    LinkedHspList { hsps, best_evalue }
}

/// The linking cutoffs of BLASTX's ungapped search with one subject
/// (`CalculateLinkHSPCutoffs`, `core::link_hsps::calculate_link_hsp_cutoffs`): the smallest Lambda, the average
/// query length and `cutoff_score_min` are those of the query batch, the reader batch's when the
/// unit is a part of one (`reader_batch.rs`). Where NCBI computes the statistics of each subject (a
/// database length of 0), `subject_cutoff_min` is the reader batch's minimum of the subject's word
/// cutoffs (`parameters::subject_word_cutoff_min`; line RG D1a-3); a batch searched as one unit
/// passes its subject's parameters as `params`. The subject is a protein, so its lengths are not
/// divided by 3.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:101-112
pub(crate) fn link_cutoffs(
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    subject_length: i32,
    db_length: i64,
    subject_cutoff_min: Option<i32>,
) -> LinkCutoffs {
    let shape = batch.reader_shape.as_deref();
    let kbp = match shape.and_then(|shape| shape.smallest_lambda) {
        Some(smallest) => smallest,
        None => {
            params
                .iter()
                .filter(|p| p.valid && p.ungapped.lambda > 0.0)
                .min_by(|a, b| a.ungapped.lambda.partial_cmp(&b.ungapped.lambda).unwrap())
                .expect("valid context")
                .ungapped
        }
    };
    let avg_query_length =
        (batch.batch_concatenated_length() as i32 - 1) / batch.batch_context_count() as i32;
    let cutoff_score_min = || match subject_cutoff_min.or(shape.and_then(|s| s.word_cutoff_min)) {
        Some(batch_min) => batch_min,
        None => params
            .iter()
            .filter(|p| p.valid)
            .map(|p| p.word_cutoff)
            .min()
            .expect("valid context"),
    };
    calculate_link_hsp_cutoffs(
        &kbp,
        avg_query_length,
        subject_length as i64,
        db_length,
        false,
        cutoff_score_min,
    )
}

#[cfg(test)]
#[path = "linking_link_cutoffs_tests.rs"]
mod link_cutoffs_tests;
