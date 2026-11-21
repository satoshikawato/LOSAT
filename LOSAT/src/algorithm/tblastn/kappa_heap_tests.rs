use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:64,
// 1560-1580; composition_adjustment/compo_heap.c:405-410:
#[test]
fn natural_112_subject_early_termination_matches_ncbi_call_state() {
    let trace = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/natural_c_d_early_20260925/early.tsv"
    ))
    .unwrap();
    let mut redone = 0;
    let mut terminated = 0;
    for row in trace
        .lines()
        .filter(|line| line.starts_with("K_EARLY_EVAL\t"))
    {
        let f: Vec<_> = row.split('\t').collect();
        assert_eq!(f[1], "0");
        let evalue = f[2].parse().unwrap();
        let n = f[3].parse().unwrap();
        let threshold = f[4].parse().unwrap();
        let worst_evalue = f[5].parse().unwrap();
        let ecutoff = f[6].parse().unwrap();
        let expected = f[7] == "1";
        let mut heap = CompoHeap::new(threshold, ecutoff).unwrap();
        heap.records = vec![
            CompoHeapRecord {
                best_evalue: worst_evalue,
                best_score: 0,
                subject_index: 0,
            };
            n
        ];
        heap.worst_evalue = worst_evalue;
        assert_eq!(compo_early_termination(evalue, &[heap]), expected);
        if expected {
            terminated += 1;
        } else {
            redone += 1;
        }
    }
    assert_eq!((redone, terminated), (11, 16));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3705-3736,2500-2515;
// composition_adjustment/compo_heap.c:252-275,330-391,439-466
#[test]
fn saved_ncbi_local_subject_heap_rank_and_pop_order() {
    for (path, threshold, expected) in [
        (
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../docs/evidence/tlosan_stage_d/kappa_traceback_20260924/run_20260923_default.tsv"
            ),
            500,
            (11, 11, 11, 0),
        ),
        (
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/run_20260925/ncbi.trace"
            ),
            2,
            (11, 10, 10, 1),
        ),
    ] {
        let trace = std::fs::read_to_string(path).unwrap();
        let mut heap = CompoHeap::new(threshold, 0.002).unwrap();
        let mut decisions = 0usize;
        let mut inserts = 0usize;
        let mut pops = 0usize;
        let mut rejected = 0usize;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "K_TRACE_HEAP_WOULD" => {
                    assert_eq!(heap.len(), f[4].parse().unwrap(), "{path}");
                    assert_eq!(heap.capacity(), f[6].parse().unwrap(), "{path}");
                    assert_eq!(
                        heap.worst_evalue().to_bits(),
                        f[8].parse::<f64>().unwrap().to_bits(),
                        "{path}"
                    );
                    let candidate = CompoHeapRecord {
                        best_evalue: f[2].parse().unwrap(),
                        best_score: f[3].parse().unwrap(),
                        subject_index: f[1].parse().unwrap(),
                    };
                    let accepted = heap.would_insert(candidate);
                    assert_eq!(accepted, f[9] == "1", "{path}");
                    if !accepted {
                        rejected += 1;
                    }
                    decisions += 1;
                }
                "K_TRACE_HEAP_INSERT" => {
                    let candidate = CompoHeapRecord {
                        best_evalue: f[2].parse().unwrap(),
                        best_score: f[3].parse().unwrap(),
                        subject_index: f[1].parse().unwrap(),
                    };
                    assert_eq!(heap.insert(candidate).map(|item| item.subject_index), None);
                    inserts += 1;
                }
                "K_TRACE_HEAP_INSERT_RETURN" => {
                    assert_eq!(heap.len(), f[3].parse().unwrap(), "{path}");
                    assert_eq!(
                        heap.worst_evalue().to_bits(),
                        f[4].parse::<f64>().unwrap().to_bits(),
                        "{path}"
                    );
                    assert_eq!(f[5], "-1");
                }
                "K_TRACE_HEAP_POP" => {
                    let actual = heap.pop();
                    assert_eq!(
                        actual.map_or(-1, |item| item.subject_index),
                        f[1].parse().unwrap(),
                        "{path}"
                    );
                    assert_eq!(heap.len(), f[2].parse().unwrap(), "{path}");
                    if let Some(item) = actual {
                        assert_eq!(
                            item.best_evalue.to_bits(),
                            f[3].parse::<f64>().unwrap().to_bits(),
                            "{path}"
                        );
                        pops += 1;
                    }
                }
                _ => {}
            }
        }
        assert_eq!((decisions, inserts, pops, rejected), expected, "{path}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101,252-275,330-391
#[test]
fn heap_threshold_ties_and_discarded_record_follow_ncbi_comparator() {
    let mut heap = CompoHeap::new(2, 0.001).unwrap();
    let rec = |oid, e, score| CompoHeapRecord {
        subject_index: oid,
        best_evalue: e,
        best_score: score,
    };
    heap.insert(rec(4, 0.1, 30));
    heap.insert(rec(3, 0.1, 30));
    assert!(!heap.would_insert(rec(2, 0.1, 30)));
    assert!(heap.would_insert(rec(5, 0.1, 30)));
    assert_eq!(heap.insert(rec(5, 0.1, 30)).unwrap().subject_index, 3);
    assert_eq!(heap.pop().unwrap().subject_index, 4);
    assert_eq!(heap.pop().unwrap().subject_index, 5);
}
