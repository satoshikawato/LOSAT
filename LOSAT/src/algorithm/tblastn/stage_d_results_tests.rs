use std::collections::HashMap;

use super::*;
use crate::algorithm::tblastn::search_gapped::GappedHsp;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2494-2515;
// core/blast_hits.c:3420-3437:
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2494-2515;
// core/blast_hits.c:3243-3297,3420-3437:
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
