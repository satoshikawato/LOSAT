//! TBLASTN gapped uneven-gap linking before preliminary E-value reap.

use std::cmp::Ordering;

use anyhow::{ensure, Result};

use super::search_gapped::GappedHsp;
use super::stage_d_stats::{LocalContextLength, LocalLinkParameters};
use crate::stats::spouge::{blast_spouge_stoe, BlastGumbelBlk};
use crate::stats::sum_statistics::{gap_decay_divisor, uneven_gap_sum_e};
use crate::stats::tables::KarlinParams;

// NCBI c++/include/algo/blast/core/blast_hits.h:126-143;
// c++/src/algo/blast/core/link_hsps.c:1765-1810:
// hsp->num = 1; Blast_HSPListGetEvalues(...); s_BlastUnevenGapLinkHSPs(...);
#[derive(Clone, Copy, Debug)]
pub(super) struct LinkedHsp {
    pub context: usize,
    pub hsp: GappedHsp,
    pub num: i32,
    pub evalue: f64,
}

// NCBI c++/src/algo/blast/core/link_hsps.c:1802-1810:
// Blast_HSPListSortByScore(hsp_list);
// hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
// if (hsp_list->hsp_array[index]->evalue < hsp_list->best_evalue)
//     hsp_list->best_evalue = hsp_list->hsp_array[index]->evalue;
#[derive(Clone, Debug)]
pub(super) struct LinkedHspList {
    pub hsps: Vec<LinkedHsp>,
    pub best_evalue: f64,
}

// NCBI c++/src/algo/blast/core/link_hsps.c:1090-1100,1469-1492:
// BlastLinkedHSPSet* prev; BlastLinkedHSPSet* next;
// double sum_score;  /* Sum bit score for the linked set. */
#[derive(Clone, Copy)]
struct WorkHsp {
    value: LinkedHsp,
    prev: Option<usize>,
    next: Option<usize>,
    sum_score: f64,
}

// NCBI c++/src/algo/blast/core/blast_hits.c:1330-1356:
// if (0 == (result = BLAST_CMP(hsp2->score, hsp1->score)) &&
//     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//     0 == (result = BLAST_CMP(hsp2->subject.end, hsp1->subject.end)) &&
//     0 == (result = BLAST_CMP(hsp1->query.offset, hsp2->query.offset))) {
//     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
// }
pub(super) fn score_compare(a: &GappedHsp, b: &GappedHsp) -> Ordering {
    b.score
        .cmp(&a.score)
        .then(a.s_start.cmp(&b.s_start))
        .then(b.s_end.cmp(&a.s_end))
        .then(a.q_start.cmp(&b.q_start))
        .then(b.q_end.cmp(&a.q_end))
}

// NCBI c++/src/algo/blast/core/link_hsps.c:1567-1594:
// if (link->queryId > hsp_array[current_index]->queryId ||
//     link->hsp->query.end > current_end) {
//     current_index = index; current_end = link->hsp->query.end;
// }
// qend_index_array[index] = current_index;
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1234-1263:
// if (hsp_array[index]->queryId < queryId) begin = index + 1;
// else if (hsp_array[index]->queryId > queryId) end = index;
// else if (hsp_array[index]->hsp->query.offset >= offset) end = index;
// else begin = index + 1;
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1282-1306:
// Int4 right_index = (begin + end) / 2;
// Int4 left_index = qend_index_array[right_index];
// if (hsp_array[left_index]->hsp->query.end >= offset)
//     end = left_index;
fn first_end_at_or_after(
    work: &[WorkHsp],
    order: &[usize],
    end_index: &[usize],
    context: usize,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, order.len());
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1310-1354:
// while (hsp_set1->prev) hsp_set1 = hsp_set1->prev;
// while (hsp_set2->prev) hsp_set2 = hsp_set2->prev;
fn chain_head(work: &[WorkHsp], mut id: usize) -> usize {
    while let Some(previous) = work[id].prev {
        id = previous;
    }
    id
}

// NCBI c++/src/algo/blast/core/link_hsps.c:1310-1354:
// if (!hsp_set2 || (hsp_set1 &&
//     hsp_set1->hsp->query.offset < hsp_set2->hsp->query.offset)) {
//     merged_hsps[index] = hsp_set1; hsp_set1 = hsp_set1->next;
// } else { merged_hsps[index] = hsp_set2; hsp_set2 = hsp_set2->next; }
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1414-1491:
// if (hsp_set1->queryId != hsp_set2->queryId) return FALSE;
// if (SIGN(hsp_set1->hsp->subject.frame) !=
//     SIGN(hsp_set2->hsp->subject.frame)) return FALSE;
// if (left_hsp->hsp->query.end < right_hsp->hsp->query.offset - gap_q) break;
// if (left_hsp->hsp->subject.end < right_hsp->hsp->subject.offset - gap_s) break;
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1117-1155:
// subject_eff_length = (Blast_SubjectIsTranslated(program_number)) ?
//     subject_length/3 : subject_length;
// subject_eff_length = MAX(subject_eff_length - len_adj, 1);
// *sum_score = new_hsp->sum_score + head_hsp->sum_score;
// BLAST_UnevenGapSumE(..., num, *sum_score, query_eff_length,
//                     subject_eff_length, ...);
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1360-1404:
// link->sum_score = sum_score;
// link->hsp->evalue = evalue;
// link->hsp->num = new_num;
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

// NCBI c++/src/algo/blast/core/link_hsps.c:1602-1757,1765-1810;
// c++/src/algo/blast/core/blast_hits.c:1811-1926:
// hsp_list->hsp_array[index]->num = 1;
// Blast_HSPListGetEvalues(..., subject_length / CODON_LENGTH, ...);
// s_BlastUnevenGapLinkHSPs(...);
// Blast_HSPListSortByScore(hsp_list);
// hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
#[allow(dead_code)] // The public TBLASTN path stays gated until C/D/E pass.
pub(super) fn link_preliminary_hsps(
    input: &[(usize, GappedHsp)],
    query_lengths: &[i32],
    lengths: &[LocalContextLength],
    subject_nt_length: i32,
    gapped_params: &[KarlinParams],
    gumbel: &BlastGumbelBlk,
    link: &LocalLinkParameters,
) -> Result<LinkedHspList> {
    // NCBI c++/src/algo/blast/core/link_hsps.c:1774-1777:
    // if (!hsp_list || hsp_list->hspcnt == 0) return 0;
    // This internal helper is entered only for the nonempty branch.
    ensure!(
        !input.is_empty(),
        "TBLASTN link helper requires nonempty HSP list"
    );
    ensure!(
        link.longest_intron > 0,
        "TBLASTN uneven-gap link parameters required"
    );
    ensure!(
        query_lengths.len() == lengths.len() && lengths.len() == gapped_params.len(),
        "one statistical context per query"
    );
    let mut work: Vec<WorkHsp> = input
        .iter()
        .map(|&(context, hsp)| {
            let kbp = &gapped_params[context];
            let evalue = blast_spouge_stoe(
                hsp.score,
                kbp,
                gumbel,
                query_lengths[context],
                subject_nt_length / 3,
            ) / gap_decay_divisor(link.gap_decay_rate, 1);
            WorkHsp {
                value: LinkedHsp {
                    context,
                    hsp,
                    num: 1,
                    evalue,
                },
                prev: None,
                next: None,
                sum_score: kbp.lambda * hsp.score as f64 - kbp.k.ln(),
            }
        })
        .collect();
    if work.len() > 1 {
        let mut score_order: Vec<usize> = (0..work.len()).collect();
        // NCBI c++/src/algo/blast/core/link_hsps.c:1217-1243:
        // if (h1->sum_score < h2->sum_score) return 1;
        // if (h1->sum_score > h2->sum_score) return -1;
        // return ScoreCompareHSPs(&h1->hsp, &h2->hsp);
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
    work.sort_by(|a, b| score_compare(&a.value.hsp, &b.value.hsp));
    let hsps: Vec<_> = work.into_iter().map(|item| item.value).collect();
    let best_evalue = hsps
        .iter()
        .map(|hsp| hsp.evalue)
        .reduce(f64::min)
        .unwrap_or(f64::MAX);
    Ok(LinkedHspList { hsps, best_evalue })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1385-1404,3071-3115,3404-3417;
// c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319
// ```c
// if (h1->hspcnt == 0 && h2->hspcnt == 0) return 0;
// else if (h1->hspcnt == 0) return 1;
// else if (h2->hspcnt == 0) return -1;
// if (evalue1 < 1.0e-180 && evalue2 < 1.0e-180) return 0;
// if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
//     return retval;
// if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
// if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
// return BLAST_CMP(h2->oid, h1->oid);
// Blast_HSPResultsReverseSort(results);
// *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
// ```
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

// NCBI c++/src/algo/blast/core/blast_engine.c:643-676 and
// c++/src/algo/blast/core/blast_hits.c:1976-2010:
// cutoff = hit_params->prelim_evalue;   /* preliminary pass */
// cutoff = hit_options->expect_value;   /* post-Kappa pass */
// if (hsp->evalue > cutoff) hsp_array[index] = Blast_HSPFree(hsp_array[index]);
// else { if (index > hsp_cnt) hsp_array[hsp_cnt] = hsp_array[index]; hsp_cnt++; }
// hsp_list->hspcnt = hsp_cnt;
#[allow(dead_code)]
pub(super) fn reap_by_evalue(list: &mut LinkedHspList, prelim_evalue: f64) {
    list.hsps.retain(|hsp| hsp.evalue <= prelim_evalue);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::stage_d_stats::{
        local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
    };
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    use crate::stats::spouge::lookup_protein_gumbel_params;
    use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};
    use std::fs;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1385-1404,3071-3115;
    // c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319
    // ```c
    // if (evalue1 < 1.0e-180 && evalue2 < 1.0e-180) return 0;
    // if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
    // return BLAST_CMP(h2->oid, h1->oid);
    // ```
    #[test]
    fn preliminary_kappa_stream_comparator_uses_evalue_score_and_oid() {
        let list = |evalue: f64, score: i32| LinkedHspList {
            hsps: vec![LinkedHsp {
                context: 0,
                hsp: GappedHsp {
                    frame: 1,
                    score,
                    q_start: 0,
                    q_end: 1,
                    q_gapped_start: 0,
                    s_start: 0,
                    s_end: 1,
                    s_gapped_start: 0,
                },
                num: 1,
                evalue,
            }],
            best_evalue: evalue,
        };
        assert_eq!(
            compare_preliminary_lists_for_kappa(0, &list(1e-90, 10), 1, &list(2e-90, 100)),
            Ordering::Less
        );
        assert_eq!(
            compare_preliminary_lists_for_kappa(0, &list(1e-90, 100), 1, &list(1e-90, 90)),
            Ordering::Less
        );
        assert_eq!(
            compare_preliminary_lists_for_kappa(2, &list(1e-90, 100), 1, &list(1e-90, 100)),
            Ordering::Less
        );
        assert_eq!(
            compare_preliminary_lists_for_kappa(0, &list(2e-181, 100), 1, &list(1e-181, 90)),
            Ordering::Less
        );
        let empty = LinkedHspList {
            hsps: Vec::new(),
            best_evalue: f64::MAX,
        };
        assert_eq!(
            compare_preliminary_lists_for_kappa(0, &empty, 1, &list(1e-90, 10)),
            Ordering::Greater
        );
        assert_eq!(
            compare_preliminary_lists_for_kappa(0, &empty, 1, &empty),
            Ordering::Equal
        );
    }

    // NCBI c++/src/algo/blast/core/link_hsps.c:1602-1810;
    // c++/src/algo/blast/core/blast_engine.c:870-906:
    // compare the complete ordered link output including num and exact E-value,
    // with both a positive uneven-gap chain and a 21-HSP no-chain list.
    #[test]
    fn saved_code1_link_lists_match_ncbi_score_order_num_and_double_evalue() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/"
        );
        for (path, subject_nt_length, query_contexts) in [
            (
                "run_20260924/multi_query_20260924_default.trace",
                6377,
                vec![(120, true), (70, true), (120, false)],
            ),
            (
                "uneven_gap_run_20260924/uneven_gap_20260924_default.trace",
                452,
                vec![(120, true)],
            ),
        ] {
            let trace = fs::read_to_string(format!("{root}{path}")).unwrap();
            let before: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("D_HSP\t0\tlink_before\t"))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[4].parse::<usize>().unwrap(),
                        GappedHsp {
                            frame: f[5].parse().unwrap(),
                            q_start: f[6].parse().unwrap(),
                            q_end: f[7].parse().unwrap(),
                            q_gapped_start: 0,
                            s_start: f[8].parse().unwrap(),
                            s_end: f[9].parse().unwrap(),
                            s_gapped_start: f[16].parse().unwrap(),
                            score: f[10].parse().unwrap(),
                        },
                    )
                })
                .collect();
            let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
            let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
            let gumbel = lookup_protein_gumbel_params(
                &ProteinScoringSpec {
                    matrix: ScoringMatrix::Blosum62,
                    gap_open: 11,
                    gap_extend: 1,
                },
                (subject_nt_length / 3) as i64,
            )
            .unwrap();
            let parameters = local_parameters_for_call(
                &query_contexts,
                subject_nt_length,
                &vec![gapped; query_contexts.len()],
                &vec![ungapped; query_contexts.len()],
                LocalParameterOptions {
                    expect_value: 10.0,
                    do_sum_stats: true,
                    max_intron_length: 0,
                    gap_trigger_bits: 22.0,
                    word_xdrop_bits: 7.0,
                    scale_factor: 1.0,
                    gumbel: Some(&gumbel),
                },
                LocalParameterCall::Initial {
                    min_subject_length: (subject_nt_length / 3) as i32,
                    composition_based_stats: 2,
                },
            );
            let query_lengths: Vec<_> = query_contexts
                .iter()
                .map(|&(length, _)| length as i32)
                .collect();
            let actual = link_preliminary_hsps(
                &before,
                &query_lengths,
                &parameters.lengths,
                subject_nt_length as i32,
                &vec![gapped; query_contexts.len()],
                &gumbel,
                &parameters.link.unwrap(),
            )
            .unwrap();
            let after: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("D_HSP\t0\tlink_after\t"))
                .collect();
            assert_eq!(actual.hsps.len(), after.len(), "{path}");
            let list_line = trace
                .lines()
                .find(|line| line.starts_with("D_LIST\t0\tlink_after\t"))
                .unwrap();
            let list_fields: Vec<_> = list_line.split('\t').collect();
            assert_eq!(
                actual.best_evalue.to_bits(),
                list_fields[6].parse::<f64>().unwrap().to_bits(),
                "{path} list best E-value"
            );
            for (index, (observed, line)) in actual.hsps.iter().zip(after).enumerate() {
                let f: Vec<_> = line.split('\t').collect();
                assert_eq!(
                    observed.context,
                    f[4].parse::<usize>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.hsp.frame,
                    f[5].parse::<i8>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.hsp.q_start,
                    f[6].parse::<i32>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.hsp.q_end,
                    f[7].parse::<i32>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.hsp.s_start,
                    f[8].parse::<i32>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.hsp.s_end,
                    f[9].parse::<i32>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.hsp.score,
                    f[10].parse::<i32>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.num,
                    f[12].parse::<i32>().unwrap(),
                    "{path} HSP {index}"
                );
                assert_eq!(
                    observed.evalue.to_bits(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    "{path} HSP {index}: actual {:.17e} expected {}",
                    observed.evalue,
                    f[13]
                );
            }
            if subject_nt_length == 6377 {
                // NCBI blast_engine.c:643-676: one of 21 preliminary HSPs has
                // E-value 50.221... > prelim_evalue 50 and is removed in place.
                let mut reaped = actual;
                reap_by_evalue(&mut reaped, parameters.prelim_evalue);
                assert_eq!(parameters.prelim_evalue, 50.0);
                assert_eq!(reaped.hsps.len(), 20);
                assert!(reaped.hsps.iter().all(|hsp| hsp.evalue <= 50.0));
            } else {
                assert_eq!(
                    actual.hsps.iter().map(|hsp| hsp.num).collect::<Vec<_>>(),
                    [2, 2]
                );
            }
        }
    }

    // NCBI c++/src/algo/blast/core/blast_kappa.c:102-114,2117-2134,
    // 3100-3109,395-445:
    // #define SCALING_FACTOR 32
    // kbp->Lambda /= scale_factor;
    // s_RescaleSearch(sbp, scoringParams, numContexts, localScalingFactor);
    // BLAST_LinkHsps(..., subject_length, ...);
    // Blast_HSPListReapByEvalue(hsp_list, hitParams->options);
    // Pinned D trace contains the post-redo input and ordered output for
    // query contexts 0 and 1. The unported redo itself is not invoked here.
    #[test]
    fn scaled_postredo_relink_and_reap_match_ncbi_trace() {
        let trace = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.trace"
        ));
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let scaled = KarlinParams {
            lambda: gapped.lambda / 32.0,
            ..gapped
        };
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            2_125,
        )
        .unwrap();
        let query_contexts = [(120, true), (70, true), (120, false)];
        let parameters = local_parameters_for_call(
            &query_contexts,
            6_377,
            &[gapped; 3],
            &[ungapped; 3],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: 2_125,
                composition_based_stats: 2,
            },
        );
        let query_lengths = [120, 70, 120];
        for (link_event, reap_event) in [(3, 5), (6, 8)] {
            let before_prefix = format!("D_HSP\t{link_event}\tlink_before\t");
            let input: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with(&before_prefix))
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[4].parse::<usize>().unwrap(),
                        GappedHsp {
                            frame: f[5].parse().unwrap(),
                            q_start: f[6].parse().unwrap(),
                            q_end: f[7].parse().unwrap(),
                            q_gapped_start: 0,
                            s_start: f[8].parse().unwrap(),
                            s_end: f[9].parse().unwrap(),
                            s_gapped_start: f[16].parse().unwrap(),
                            score: f[10].parse().unwrap(),
                        },
                    )
                })
                .collect();
            let mut actual = link_preliminary_hsps(
                &input,
                &query_lengths,
                &parameters.lengths,
                6_377,
                &[scaled; 3],
                &gumbel,
                &parameters.link.unwrap(),
            )
            .unwrap();
            let after_prefix = format!("D_HSP\t{link_event}\tlink_after\t");
            let expected: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with(&after_prefix))
                .collect();
            assert_eq!(actual.hsps.len(), expected.len(), "link event {link_event}");
            for (i, (observed, row)) in actual.hsps.iter().zip(expected).enumerate() {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    observed.context,
                    f[4].parse::<usize>().unwrap(),
                    "event {link_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.score,
                    f[10].parse::<i32>().unwrap(),
                    "event {link_event} HSP {i}"
                );
                assert_eq!(
                    observed.num,
                    f[12].parse::<i32>().unwrap(),
                    "event {link_event} HSP {i}"
                );
                assert_eq!(
                    observed.evalue.to_bits(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    "event {link_event} HSP {i}"
                );
            }
            let list_prefix = format!("D_LIST\t{link_event}\tlink_after\t");
            let list_row = trace
                .lines()
                .find(|line| line.starts_with(&list_prefix))
                .unwrap();
            assert_eq!(
                actual.best_evalue.to_bits(),
                list_row
                    .split('\t')
                    .nth(6)
                    .unwrap()
                    .parse::<f64>()
                    .unwrap()
                    .to_bits()
            );
            let reap_before_best = actual.best_evalue.to_bits();
            reap_by_evalue(&mut actual, 10.0);
            let reap_prefix = format!("D_HSP\t{reap_event}\treap_after\t");
            let kept: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with(&reap_prefix))
                .collect();
            assert_eq!(actual.hsps.len(), kept.len(), "reap event {reap_event}");
            assert_eq!(actual.best_evalue.to_bits(), reap_before_best);
            for (i, (observed, row)) in actual.hsps.iter().zip(kept).enumerate() {
                let f: Vec<_> = row.split('\t').collect();
                assert_eq!(
                    observed.context,
                    f[4].parse::<usize>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.frame,
                    f[5].parse::<i8>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.q_start,
                    f[6].parse::<i32>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.q_end,
                    f[7].parse::<i32>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.s_start,
                    f[8].parse::<i32>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.s_end,
                    f[9].parse::<i32>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.hsp.score,
                    f[10].parse::<i32>().unwrap(),
                    "reap event {reap_event} HSP {i}"
                );
                assert_eq!(
                    observed.evalue.to_bits(),
                    f[13].parse::<f64>().unwrap().to_bits(),
                    "reap event {reap_event} HSP {i}"
                );
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-906;
    // c++/src/algo/blast/core/link_hsps.c:1765-1810;
    // c++/src/algo/blast/core/blast_hspstream.c:144-152,289-319:
    // BLAST_LinkHsps(..., hsp_list_out, ...);
    // s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
    // Blast_HSPResultsReverseSort(hsp_stream->results);
    // *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
    #[test]
    fn natural_heap_replacement_fixture_matches_ncbi_preliminary_stream() {
        use crate::algorithm::tblastn::kappa_heap::{
            compo_early_termination, CompoHeap, CompoHeapRecord,
        };
        use crate::algorithm::tblastn::search_gapped::{
            preliminary_protein_hsps_in_ncbi_order, PreliminaryProfile,
        };
        use crate::utils::seg::SegParams;
        use std::collections::{HashMap, HashSet};

        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/"
        );
        let read_fasta = |name: &str| -> Vec<Vec<u8>> {
            let mut records: Vec<Vec<u8>> = Vec::new();
            for line in fs::read_to_string(format!("{root}{name}")).unwrap().lines() {
                if line.starts_with('>') {
                    records.push(Vec::new());
                } else {
                    records.last_mut().unwrap().extend(line.bytes());
                }
            }
            records
        };
        let query = read_fasta("query.faa").remove(0);
        let subjects = read_fasta("subjects.fna");
        let total_nt: usize = subjects.iter().map(Vec::len).sum();
        assert_eq!((query.len(), subjects.len(), total_nt), (120, 112, 40320));

        let trace =
            fs::read_to_string(format!("{root}natural_c_d_early_20260925/calls.tsv")).unwrap();
        let params_trace =
            fs::read_to_string(format!("{root}natural_c_d_early_20260925/parameters.tsv")).unwrap();
        let early_trace =
            fs::read_to_string(format!("{root}natural_c_d_early_20260925/early.tsv")).unwrap();
        let kappa_trace =
            fs::read_to_string(format!("{root}result_order_20260925/ncbi.trace")).unwrap();
        let prelim_events: HashSet<usize> = trace
            .lines()
            .filter(|line| line.starts_with("D_CALL\t") && line.split('\t').nth(2) == Some("link"))
            .take(subjects.len())
            .map(|line| line.split('\t').nth(1).unwrap().parse().unwrap())
            .collect();
        assert_eq!(prelim_events.len(), 112);
        let mut event_oid = HashMap::new();
        for line in trace.lines().filter(|line| line.starts_with("D_LIST\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let event: usize = f[1].parse().unwrap();
            if prelim_events.contains(&event) && f[2] == "link_before" && f[3] != "NULL" {
                event_oid.insert(event, f[3].parse::<usize>().unwrap());
            }
        }
        let mut before: HashMap<usize, Vec<Vec<&str>>> = HashMap::new();
        let mut after: HashMap<usize, Vec<Vec<&str>>> = HashMap::new();
        for line in trace.lines().filter(|line| line.starts_with("D_HSP\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let event: usize = f[1].parse().unwrap();
            let Some(&oid) = event_oid.get(&event) else {
                continue;
            };
            match f[2] {
                "link_before" => before.entry(oid).or_default().push(f),
                "link_after" => after.entry(oid).or_default().push(f),
                _ => {}
            }
        }
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (total_nt / 3) as i64,
        )
        .unwrap();
        let parameters = local_parameters_for_call(
            &[(query.len(), true)],
            total_nt,
            &[gapped],
            &[ungapped],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: 120,
                composition_based_stats: 2,
            },
        );
        assert!(params_trace.contains("D_PARAM_CONTEXT\t0\t25\t25\n"));
        assert_eq!(
            (
                parameters.cutoffs[0].word_cutoff,
                parameters.cutoffs[0].hit_cutoff
            ),
            (25, 25)
        );
        let seg = SegParams::default();
        let word_xdrop = [parameters.cutoffs[0].word_xdrop];
        let word_cutoff = [parameters.cutoffs[0].word_cutoff];
        let hit_cutoff = [parameters.cutoffs[0].hit_cutoff];
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &hit_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        let mut retained: Vec<(usize, LinkedHspList)> = Vec::new();
        for (oid, subject) in subjects.iter().enumerate() {
            let (preliminary, _) =
                preliminary_protein_hsps_in_ncbi_order(&[&query], subject, 1, profile).unwrap();
            let expected = before.get(&oid).map(Vec::as_slice).unwrap_or(&[]);
            assert_eq!(preliminary.len(), expected.len(), "OID {oid} prelink count");
            for ((context, hsp), row) in preliminary.iter().zip(expected) {
                assert_eq!(*context, row[4].parse().unwrap(), "OID {oid} context");
                assert_eq!(hsp.frame, row[5].parse().unwrap(), "OID {oid} frame");
                assert_eq!(hsp.q_start, row[6].parse().unwrap(), "OID {oid} qstart");
                assert_eq!(hsp.q_end, row[7].parse().unwrap(), "OID {oid} qend");
                assert_eq!(hsp.s_start, row[8].parse().unwrap(), "OID {oid} sstart");
                assert_eq!(hsp.s_end, row[9].parse().unwrap(), "OID {oid} send");
                assert_eq!(hsp.score, row[10].parse().unwrap(), "OID {oid} score");
                assert_eq!(
                    hsp.s_gapped_start,
                    row[16].parse().unwrap(),
                    "OID {oid} gapped start"
                );
            }
            if preliminary.is_empty() {
                continue;
            }
            let mut linked = link_preliminary_hsps(
                &preliminary,
                &[query.len() as i32],
                &parameters.lengths,
                subject.len() as i32,
                &[gapped],
                &gumbel,
                parameters.link.as_ref().unwrap(),
            )
            .unwrap();
            let expected_linked = after.get(&oid).map(Vec::as_slice).unwrap_or(&[]);
            assert_eq!(
                linked.hsps.len(),
                expected_linked.len(),
                "OID {oid} linked count"
            );
            for (hsp, row) in linked.hsps.iter().zip(expected_linked) {
                assert_eq!(
                    hsp.hsp.score,
                    row[10].parse().unwrap(),
                    "OID {oid} linked score"
                );
                assert_eq!(hsp.num, row[12].parse().unwrap(), "OID {oid} linked num");
                assert_eq!(
                    hsp.evalue.to_bits(),
                    row[13].parse::<f64>().unwrap().to_bits(),
                    "OID {oid} linked E-value"
                );
            }
            reap_by_evalue(&mut linked, parameters.prelim_evalue);
            if !linked.hsps.is_empty() {
                retained.push((oid, linked));
            }
        }
        retained.sort_by(|(oid_a, a), (oid_b, b)| {
            compare_preliminary_lists_for_kappa(*oid_a as i32, a, *oid_b as i32, b)
        });
        let expected_stream: Vec<_> = early_trace
            .lines()
            .filter(|line| line.starts_with("K_EARLY_STREAM\t0\t"))
            .collect();
        assert_eq!(retained.len(), 27);
        assert_eq!(retained.len(), expected_stream.len());
        for ((oid, linked), row) in retained.iter().zip(expected_stream) {
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(*oid, f[2].parse().unwrap());
            assert_eq!(linked.hsps.len(), f[5].parse().unwrap());
            assert_eq!(
                linked.best_evalue.to_bits(),
                f[4].parse::<f64>().unwrap().to_bits(),
                "OID {oid} stream best E-value"
            );
        }
        let early_decisions: Vec<_> = early_trace
            .lines()
            .filter(|line| line.starts_with("K_EARLY_EVAL\t"))
            .collect();
        assert_eq!(early_decisions.len(), retained.len());
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3383-3427,
        // 3525-3535; composition_adjustment/redo_alignment.c:64,1560-1582;
        // composition_adjustment/compo_heap.c:252-275,330-391,405-410:
        // if (BlastCompo_EarlyTermination(localMatch->best_evalue,
        //     redoneMatches, numQueries)) { Blast_HSPListFree(localMatch); continue; }
        // if (BlastCompo_HeapWouldInsert(...)) BlastCompo_HeapInsert(...);
        let would_rows: Vec<_> = kappa_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_WOULD\t"))
            .collect();
        let insert_rows: Vec<_> = kappa_trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_INSERT\t"))
            .collect();
        let mut heap = CompoHeap::new(2, 0.002).unwrap();
        let mut redo_index = 0;
        let mut insert_index = 0;
        for ((oid, linked), row) in retained.iter().zip(early_decisions) {
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(
                linked.best_evalue.to_bits(),
                f[2].parse::<f64>().unwrap().to_bits(),
                "OID {oid} early-termination input"
            );
            assert_eq!(heap.len(), f[3].parse().unwrap(), "OID {oid} heap count");
            assert_eq!(
                heap.worst_evalue().to_bits(),
                f[5].parse::<f64>().unwrap().to_bits(),
                "OID {oid} heap worst E-value"
            );
            let terminated =
                compo_early_termination(linked.best_evalue, std::slice::from_ref(&heap));
            assert_eq!(terminated, f[7] == "1", "OID {oid} early decision");
            if terminated {
                continue;
            }
            let would: Vec<_> = would_rows[redo_index].split('\t').collect();
            assert_eq!(*oid, would[1].parse::<usize>().unwrap());
            let expected_prelim: Vec<_> = kappa_trace
                .lines()
                .filter(|line| line.starts_with(&format!("K_TRACE_PRELIM\t{redo_index}\t")))
                .collect();
            assert_eq!(linked.hsps.len(), expected_prelim.len());
            for (linked, row) in linked.hsps.iter().zip(expected_prelim) {
                let p: Vec<_> = row.split('\t').collect();
                assert_eq!(linked.hsp.score, p[3].parse().unwrap());
                assert_eq!(linked.hsp.frame, p[5].parse().unwrap());
                assert_eq!(linked.hsp.q_start, p[6].parse().unwrap());
                assert_eq!(linked.hsp.q_end, p[7].parse().unwrap());
                assert_eq!(linked.hsp.q_gapped_start, p[8].parse().unwrap());
                assert_eq!(linked.hsp.s_start, p[9].parse().unwrap());
                assert_eq!(linked.hsp.s_end, p[10].parse().unwrap());
                assert_eq!(linked.hsp.s_gapped_start, p[11].parse().unwrap());
            }
            let candidate = CompoHeapRecord {
                subject_index: *oid as i32,
                best_evalue: would[2].parse().unwrap(),
                best_score: would[3].parse().unwrap(),
            };
            assert_eq!(heap.would_insert(candidate), would[9] == "1");
            if would[9] == "1" {
                let insert: Vec<_> = insert_rows[insert_index].split('\t').collect();
                assert_eq!(*oid, insert[1].parse::<usize>().unwrap());
                assert_eq!(
                    candidate.best_evalue.to_bits(),
                    insert[2].parse::<f64>().unwrap().to_bits()
                );
                assert_eq!(candidate.best_score, insert[3].parse().unwrap());
                assert!(heap.insert(candidate).is_none());
                insert_index += 1;
            }
            redo_index += 1;
        }
        assert_eq!((redo_index, insert_index, heap.len()), (11, 10, 10));
    }
}
