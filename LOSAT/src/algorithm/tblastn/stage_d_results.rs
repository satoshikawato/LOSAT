//! Internal TBLASTN Kappa result hitlist order and replacement boundary.

use std::cmp::Ordering;

use anyhow::{ensure, Result};

use crate::common::GapEditOp;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;

use super::stage_d_linking::{
    compare_preliminary_lists_for_kappa, score_compare, LinkedHsp, LinkedHspList,
};

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.h:152-181;
// core/blast_kappa.c:2494-2515:
// while (NULL != (hsp_list = BlastCompo_HeapPop(heap)))
//     Blast_HitListUpdate(hitlist, hsp_list);
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:305-358,
// 3690-3706; core/blast_hits.h:126-143:
// GapEditScript *editScript = align->context; align->context = NULL;
// Blast_HSPInit(..., &editScript, &new_hsp);
// new_hsp->num_ident = 0; /* filled after postredo evaluation */
// s_HSPListNormalizeScores(...); s_ComputeNumIdentities(...);
#[derive(Clone, Debug, PartialEq)]
pub(super) struct KappaHspPayload {
    pub bit_score: f64,
    pub num_ident: i32,
    // NCBI c++/src/algo/blast/core/blast_kappa.c:515-526;
    // c++/src/algo/blast/core/blast_hits.c:767-811,966-989:
    // Blast_HSPGetNumIdentitiesAndPositives(query, target_sequence, hsp,
    //                                       scoring_options, 0, sbp);
    // *num_pos_ptr = num_pos + num_ident;
    // Alignment length and gap counts follow the same gap_info operation lengths.
    pub num_positives: usize,
    pub align_length: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub gap_letters: usize,
    pub edit_script: Vec<GapEditOp>,
    pub matrix_adjust_rule: EMatrixAdjustRule,
}

#[derive(Clone, Debug)]
pub(super) struct KappaResultList {
    pub oid: i32,
    pub hsps: LinkedHspList,
    // One payload per HSP, in the same order as `hsps` after each sort.
    // Empty only in the older numeric-only comparison fixtures.
    pub payloads: Vec<KappaHspPayload>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3125-3135,
// 3243-3297,3420-3437:
// new_hitlist->hsplist_max = hitlist_size;
// new_hitlist->low_score = INT4_MAX;
// new_hitlist->hsplist_count = 0;
#[allow(dead_code)] // The public TBLASTN result path remains gated until C/D/E pass.
pub(super) struct KappaResultHitList {
    lists: Vec<KappaResultList>,
    max: usize,
    worst_evalue: f64,
    low_score: i32,
    heapified: bool,
}

#[allow(dead_code)] // The public TBLASTN result path remains gated until C/D/E pass.
impl KappaResultHitList {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3125-3135:
    // new_hitlist->hsplist_max = hitlist_size;
    // new_hitlist->low_score = INT4_MAX;
    // new_hitlist->hsplist_count = 0;
    pub(super) fn new(hitlist_size: usize) -> Result<Self> {
        ensure!(hitlist_size > 0, "TBLASTN hitlist size must be positive");
        Ok(Self {
            lists: Vec::new(),
            max: hitlist_size,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
        })
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3243-3297:
    // hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    // if (hit_list->hsplist_count < hit_list->hsplist_max)
    //     hit_list->hsplist_array[hit_list->hsplist_count++] = hsp_list;
    // else {
    //     if (!hit_list->heapified) s_CreateHeap(..., s_EvalueCompareHSPLists);
    //     evalue_order = s_EvalueCompareHSPLists(&hit_list->hsplist_array[0], &hsp_list);
    //     if (evalue_order < 0) Blast_HSPListFree(hsp_list);
    //     else s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
    // }
    pub(super) fn update(
        &mut self,
        mut incoming: KappaResultList,
    ) -> Result<Option<KappaResultList>> {
        ensure!(
            !incoming.hsps.hsps.is_empty(),
            "Kappa result list must contain an HSP"
        );
        ensure!(
            incoming.payloads.is_empty() || incoming.payloads.len() == incoming.hsps.hsps.len(),
            "Kappa result payload count must match HSP count"
        );
        refresh_best_evalue(&mut incoming.hsps);
        if self.lists.len() < self.max {
            self.worst_evalue = self.worst_evalue.max(incoming.hsps.best_evalue);
            self.low_score = self.low_score.min(incoming.hsps.hsps[0].hsp.score);
            self.lists.push(incoming);
            return Ok(None);
        }
        if !self.heapified {
            for list in &mut self.lists {
                sort_hsps_by_evalue(list);
                refresh_best_evalue(&mut list.hsps);
            }
            for start in (0..self.lists.len() / 2).rev() {
                heapify_down(&mut self.lists, start);
            }
            self.heapified = true;
        }
        sort_hsps_by_evalue(&mut incoming);
        refresh_best_evalue(&mut incoming.hsps);
        let discarded = if compare_result_lists(&self.lists[0], &incoming) == Ordering::Less {
            incoming
        } else {
            let old = std::mem::replace(&mut self.lists[0], incoming);
            heapify_down(&mut self.lists, 0);
            old
        };
        self.worst_evalue = self.lists[0].hsps.best_evalue;
        self.low_score = self.lists[0].hsps.hsps[0].hsp.score;
        Ok(Some(discarded))
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3420-3437:
    // for (index1 = 0; index1 < hit_list->hsplist_count/2; ++index1) {
    //     hsplist = hit_list->hsplist_array[index1];
    //     hit_list->hsplist_array[index1] =
    //         hit_list->hsplist_array[hit_list->hsplist_count-index1-1];
    //     hit_list->hsplist_array[hit_list->hsplist_count-index1-1] = hsplist;
    // }
    pub(super) fn reverse_order(&mut self) {
        self.lists.reverse();
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.h:169-181:
    // Int4 hsplist_count, hsplist_max;
    // double worst_evalue; Int4 low_score; Boolean heapified;
    // BlastHSPList** hsplist_array;
    pub(super) fn lists(&self) -> &[KappaResultList] {
        &self.lists
    }

    pub(super) fn worst_evalue(&self) -> f64 {
        self.worst_evalue
    }

    pub(super) fn low_score(&self) -> i32 {
        self.low_score
    }

    pub(super) fn heapified(&self) -> bool {
        self.heapified
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1734-1750:
// double best_evalue = (double) INT4_MAX;
// for (index=0; index<hsp_list->hspcnt; index++)
//     best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
fn refresh_best_evalue(list: &mut LinkedHspList) {
    list.best_evalue = list
        .hsps
        .iter()
        .fold(i32::MAX as f64, |best, hsp| best.min(hsp.evalue));
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1385-1455:
// if (evalue1 < 1.0e-180 && evalue2 < 1.0e-180) return 0;
// if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0) return retval;
// return ScoreCompareHSPs(v1, v2);
fn sort_hsps_by_evalue(list: &mut KappaResultList) {
    let compare = |a: &LinkedHsp, b: &LinkedHsp| {
        let evalue = if a.evalue < 1.0e-180 && b.evalue < 1.0e-180 {
            Ordering::Equal
        } else if a.evalue < b.evalue {
            Ordering::Less
        } else if a.evalue > b.evalue {
            Ordering::Greater
        } else {
            Ordering::Equal
        };
        evalue.then_with(|| score_compare(&a.hsp, &b.hsp))
    };
    if list.payloads.is_empty() {
        list.hsps.hsps.sort_by(compare);
    } else {
        // NCBI c++/src/algo/blast/core/blast_hits.c:1385-1455:
        // return ScoreCompareHSPs(v1, v2); /* HSP pointer sort */
        // Thus each
        // bit score, identity count and edit script moves with its HSP.
        let mut paired: Vec<_> = std::mem::take(&mut list.hsps.hsps)
            .into_iter()
            .zip(std::mem::take(&mut list.payloads))
            .collect();
        paired.sort_by(|a, b| compare(&a.0, &b.0));
        let (hsps, payloads) = paired.into_iter().unzip();
        list.hsps.hsps = hsps;
        list.payloads = payloads;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3071-3111:
// if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0) return retval;
// if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
// if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
// return BLAST_CMP(h2->oid, h1->oid);
fn compare_result_lists(a: &KappaResultList, b: &KappaResultList) -> Ordering {
    compare_preliminary_lists_for_kappa(a.oid, &a.hsps, b.oid, &b.hsps)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1627-1676,
// 3194-3210:
// large_son = (*compar)(left_son, left_son+width) >= 0 ? left_son : left_son+width;
// if ((*compar)(base, large_son) < 0) SWAP(base, large_son);
fn heapify_down(lists: &mut [KappaResultList], start: usize) {
    let mut root = start;
    loop {
        let left = root * 2 + 1;
        if left >= lists.len() {
            break;
        }
        let right = left + 1;
        let large = if right < lists.len()
            && compare_result_lists(&lists[left], &lists[right]) == Ordering::Less
        {
            right
        } else {
            left
        };
        if compare_result_lists(&lists[root], &lists[large]) != Ordering::Less {
            break;
        }
        lists.swap(root, large);
        root = large;
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use crate::algorithm::tblastn::search_gapped::GappedHsp;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2494-2515;
    // core/blast_hits.c:3420-3437:
    // for (query_index = 0; query_index < num_queries; query_index++) {
    //     while ((hsp_list = BlastCompo_HeapPop(heap)) != NULL)
    //         Blast_HitListUpdate(hitlist, hsp_list);
    // }
    // Blast_HSPResultsReverseOrder(results);
    #[test]
    fn saved_multi_query_result_updates_follow_ncbi_query_order() {
        let trace = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/multi_query_20260924_default.tsv"
        ))
        .unwrap();
        let mut results = [
            KappaResultHitList::new(500).unwrap(),
            KappaResultHitList::new(500).unwrap(),
            KappaResultHitList::new(500).unwrap(),
        ];
        let mut query_index = 0;
        let mut updates = [0usize; 3];
        for row in trace.lines() {
            let f: Vec<_> = row.split('\t').collect();
            match f[0] {
                "K_TRACE_HEAP_POP" if f[1] == "-1" && query_index < 3 => query_index += 1,
                "K_TRACE_RESULT_UPDATE_IN" => {
                    assert!(query_index < 3);
                    let result = &mut results[query_index];
                    assert_eq!(result.lists().len(), f[2].parse().unwrap());
                    assert_eq!(
                        result.worst_evalue().to_bits(),
                        f[4].parse::<f64>().unwrap().to_bits()
                    );
                    assert_eq!(result.low_score(), f[5].parse().unwrap());
                    let score = f[8].parse().unwrap();
                    let evalue = f[7].parse().unwrap();
                    result
                        .update(KappaResultList {
                            oid: f[1].parse().unwrap(),
                            hsps: LinkedHspList {
                                hsps: vec![LinkedHsp {
                                    context: query_index,
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
                                    source_index: 0,
                                }],
                                best_evalue: evalue,
                            },
                            payloads: Vec::new(),
                        })
                        .unwrap();
                    updates[query_index] += 1;
                }
                "K_TRACE_RESULT_UPDATE_OUT" => {
                    let result = &results[query_index];
                    assert_eq!(result.lists().len(), f[3].parse().unwrap());
                    assert_eq!(
                        result.worst_evalue().to_bits(),
                        f[4].parse::<f64>().unwrap().to_bits()
                    );
                    assert_eq!(result.low_score(), f[5].parse().unwrap());
                }
                _ => {}
            }
        }
        assert_eq!(updates, [1, 1, 0]);
        for query_index in 0..3 {
            results[query_index].reverse_order();
            let prefix = format!("K_TRACE_RESULT_REVERSE_OUT\t{query_index}\t");
            let row = trace
                .lines()
                .find(|line| line.starts_with(&prefix))
                .unwrap();
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(results[query_index].lists().len(), f[2].parse().unwrap());
            assert_eq!(
                results[query_index]
                    .lists()
                    .iter()
                    .map(|list| list.oid)
                    .collect::<Vec<_>>(),
                f[3..]
                    .iter()
                    .map(|value| value.parse::<i32>().unwrap())
                    .collect::<Vec<_>>()
            );
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2494-2515;
    // core/blast_hits.c:3243-3297,3420-3437:
    // while (NULL != (hsp_list = BlastCompo_HeapPop(heap)))
    //     Blast_HitListUpdate(hitlist, hsp_list);
    // Blast_HSPResultsReverseOrder(results);
    #[test]
    fn saved_local_kappa_pop_hitlist_update_and_reverse_match_ncbi() {
        for (trace, hitlist_size) in [
            (
                concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/run_20260923_default.tsv"),
                500,
            ),
            (
                concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/seg_hard_query_20260924_default.tsv"),
                500,
            ),
            (
                concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/result_order_20260925/ncbi.trace"),
                2,
            ),
        ] {
            let trace_text = std::fs::read_to_string(trace).unwrap();
            let mut hsp_map: HashMap<i32, Vec<LinkedHsp>> = HashMap::new();
            for row in trace_text.lines().filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t")) {
                let f: Vec<_> = row.split('\t').collect();
                hsp_map.entry(f[1].parse().unwrap()).or_default().push(LinkedHsp {
                    context: f[7].parse().unwrap(),
                    hsp: GappedHsp {
                        frame: f[8].parse().unwrap(),
                        score: f[3].parse().unwrap(),
                        q_start: f[9].parse().unwrap(),
                        q_end: f[10].parse().unwrap(),
                        q_gapped_start: 0,
                        s_start: f[11].parse().unwrap(),
                        s_end: f[12].parse().unwrap(),
                        s_gapped_start: 0,
                    },
                    num: f[13].parse().unwrap(),
                    evalue: f[5].parse().unwrap(),
                    source_index: 0,
                });
            }
            let mut results = KappaResultHitList::new(hitlist_size).unwrap();
            let mut updates = 0;
            let mut reversed = false;
            for row in trace_text.lines().filter(|row| row.starts_with("K_TRACE_RESULT_")) {
                let f: Vec<_> = row.split('\t').collect();
                let expected_oids = |fields: &[&str]| -> Vec<i32> {
                    fields.iter().map(|value| value.parse().unwrap()).collect()
                };
                let actual_oids = |results: &KappaResultHitList| -> Vec<i32> {
                    results.lists().iter().map(|list| list.oid).collect()
                };
                match f[0] {
                    "K_TRACE_RESULT_UPDATE_IN" => {
                        assert_eq!(results.lists().len(), f[2].parse().unwrap(), "{trace}");
                        assert_eq!(hitlist_size, f[3].parse().unwrap(), "{trace}");
                        assert_eq!(results.worst_evalue().to_bits(), f[4].parse::<f64>().unwrap().to_bits(), "{trace}");
                        assert_eq!(results.low_score(), f[5].parse().unwrap(), "{trace}");
                        assert_eq!(results.heapified(), f[6] == "1", "{trace}");
                        let oid = f[1].parse().unwrap();
                        let hsps = hsp_map.remove(&oid).unwrap();
                        assert_eq!(hsps[0].hsp.score, f[8].parse().unwrap(), "{trace}");
                        results.update(KappaResultList {
                            oid,
                            hsps: LinkedHspList { hsps, best_evalue: f[7].parse().unwrap() },
                            payloads: Vec::new(),
                        }).unwrap();
                        updates += 1;
                    }
                    "K_TRACE_RESULT_UPDATE_OUT" => {
                        assert_eq!(f[2], "0", "{trace}");
                        assert_eq!(results.lists().len(), f[3].parse().unwrap(), "{trace}");
                        assert_eq!(results.worst_evalue().to_bits(), f[4].parse::<f64>().unwrap().to_bits(), "{trace}");
                        assert_eq!(results.low_score(), f[5].parse().unwrap(), "{trace}");
                        assert_eq!(results.heapified(), f[6] == "1", "{trace}");
                        assert_eq!(actual_oids(&results), expected_oids(&f[7..]), "{trace}");
                    }
                    "K_TRACE_RESULT_HSP" => {
                        let list = &results.lists()[f[1].parse::<usize>().unwrap()];
                        assert_eq!(list.oid, f[2].parse().unwrap(), "{trace}");
                        let hsp = &list.hsps.hsps[f[3].parse::<usize>().unwrap()];
                        assert_eq!(hsp.hsp.score, f[4].parse().unwrap(), "{trace}");
                        assert_eq!(hsp.evalue.to_bits(), f[5].parse::<f64>().unwrap().to_bits(), "{trace}");
                        assert_eq!(hsp.hsp.q_start, f[6].parse().unwrap(), "{trace}");
                        assert_eq!(hsp.hsp.q_end, f[7].parse().unwrap(), "{trace}");
                        assert_eq!(hsp.hsp.s_start, f[8].parse().unwrap(), "{trace}");
                        assert_eq!(hsp.hsp.s_end, f[9].parse().unwrap(), "{trace}");
                        assert_eq!(hsp.hsp.frame, f[10].parse().unwrap(), "{trace}");
                        assert_eq!(hsp.num, f[11].parse().unwrap(), "{trace}");
                    }
                    "K_TRACE_RESULT_REVERSE_IN" => {
                        assert_eq!(f[1], "0", "{trace}");
                        assert_eq!(actual_oids(&results), expected_oids(&f[3..]), "{trace}");
                        results.reverse_order();
                        reversed = true;
                    }
                    "K_TRACE_RESULT_REVERSE_OUT" => {
                        assert_eq!(f[1], "0", "{trace}");
                        assert_eq!(actual_oids(&results), expected_oids(&f[3..]), "{trace}");
                    }
                    _ => panic!("unexpected result trace row: {row}"),
                }
            }
            assert!(updates > 0 && reversed, "{trace}");
        }
    }
}
