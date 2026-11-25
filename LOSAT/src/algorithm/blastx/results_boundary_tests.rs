use super::*;
use std::fmt::Write;
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:75-85
fn dump(out: &mut String, stage: &str, oid: i32, result: i32, h: &CompoHeap) {
    write!(
        out,
        "H\t{stage}\t{oid}\t{result}\t{}\t{:016x}\t{}",
        h.entries.len(),
        h.worst.to_bits(),
        i32::from(h.heapified)
    )
    .unwrap();
    for e in &h.entries {
        write!(
            out,
            "\t{}:{}:{:016x}",
            e.list.oid,
            e.best_score,
            e.list.best_evalue.to_bits()
        )
        .unwrap();
    }
    out.push('\n');
}
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:439-461
#[test]
fn pinned_c_heap_threshold_cutoff_tie_and_empty_pop() {
    let expected = include_str!("../../../tests/unit/blastx_d_boundary_expected.tsv")
        .lines()
        .filter(|l| !l.starts_with("C\t"))
        .map(|l| format!("{l}\n"))
        .collect::<String>();
    let mut actual = String::new();
    for threshold in 1..=2 {
        let mut heap = CompoHeap::new(threshold);
        writeln!(actual, "T\t{threshold}").unwrap();
        dump(&mut actual, "INIT", -1, 0, &heap);
        let es = [0.1, 0.1, 0.2, 0.09, 0.0021, 0.002, 0.001, 0.002, 0.003];
        let scores = [100, 100, 101, 99, 99, 98, 97, 98, 1000];
        for (i, (&e, &score)) in es.iter().zip(&scores).enumerate() {
            let yes = heap.would_insert(e, score, i);
            dump(&mut actual, "WOULD", i as i32, i32::from(yes), &heap);
            if yes {
                let discarded = heap.insert(HeapEntry {
                    list: HspList {
                        query_index: 0,
                        oid: i,
                        hsps: Vec::new(),
                        best_evalue: e,
                    },
                    best_score: score,
                });
                dump(
                    &mut actual,
                    "INSERT",
                    i as i32,
                    discarded.map_or(-1, |e| e.list.oid as i32),
                    &heap,
                );
            }
        }
        loop {
            let entry = heap.pop();
            dump(
                &mut actual,
                "POP",
                entry.as_ref().map_or(-1, |e| e.list.oid as i32),
                0,
                &heap,
            );
            if entry.is_none() {
                break;
            }
        }
    }
    assert_eq!(actual, expected);
}

// `-use_sw_tback` in e-value mode (line RG item X4b-2): a window's query context indexes the batch's
// heaps, one per query. Inside the array the heap's own answer: a full heap turns away a window that
// is no better than its worst entry. Past the array an empty heap (approved exception 23 of
// PD-LOSAT-NCBI-DEFECTS), which takes every window, whatever the heaps of the batch hold.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1465-1474
#[test]
fn a_context_past_the_batch_tests_an_empty_heap() {
    let full = || {
        let mut heap = CompoHeap::new(1);
        heap.insert(HeapEntry {
            list: HspList {
                query_index: 0,
                oid: 0,
                hsps: Vec::new(),
                best_evalue: 0.1,
            },
            best_score: 50,
        });
        heap
    };
    let mut heaps = vec![full(), full()];
    // Context 1 reads query 1's full heap: a worse window is turned away, a better one taken.
    assert!(!CompoHeap::would_insert_at(&mut heaps, 1, 0.5, 40, 1));
    assert!(CompoHeap::would_insert_at(&mut heaps, 1, 0.05, 60, 1));
    // Contexts 2 to 11 (frames of the two queries past the array): an empty heap's answer.
    for context in 2..12 {
        assert!(
            CompoHeap::would_insert_at(&mut heaps, context, 0.5, 40, 1),
            "{context}"
        );
    }
    assert!(heaps.iter().all(|heap| heap.entries.len() == 1));
}

// No NCBI counterpart (scheduling only): a heap filled to its cutoff (`CompoHeap::early`) is one with
// at least `threshold` entries at or below the cutoff, and each insert at or below it adds one, so
// `x_early_ahead` (the guess of the composition redo's matches redone ahead) with `ahead[q]` matches
// counted ahead is `early` on the heaps after that many more such inserts, and `early` itself with
// none. Random inserts into two heaps of random thresholds, E-values about the cutoff and five times
// it, ties included.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:404-409
#[test]
fn filled_to_cutoff_counts_the_entries_at_or_below_it() {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut random = move |below: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % below as u64) as usize
    };
    let evalues = [
        0.0, 1e-30, 0.001, 0.002, 0.0021, 0.005, 0.01, 0.0101, 0.1, 1.0,
    ];
    // A heap that `insert` made from empty, as the redo makes it (an insert where `would_insert`).
    let build = |threshold: usize, inserts: &[(f64, i32)]| {
        let mut heap = CompoHeap::new(threshold);
        for (oid, &(e, score)) in inserts.iter().enumerate() {
            if heap.would_insert(e, score, oid) {
                heap.insert(HeapEntry {
                    list: HspList {
                        query_index: 0,
                        oid,
                        hsps: Vec::new(),
                        best_evalue: e,
                    },
                    best_score: score,
                });
            }
        }
        heap
    };
    for round in 0..500 {
        let shapes: Vec<(usize, Vec<(f64, i32)>)> = (0..2)
            .map(|_| {
                let threshold = 1 + random(5);
                let inserts = (0..random(20))
                    .map(|_| (evalues[random(evalues.len())], random(40) as i32))
                    .collect();
                (threshold, inserts)
            })
            .collect();
        let heaps: Vec<CompoHeap> = shapes.iter().map(|(t, i)| build(*t, i)).collect();
        for heap in &heaps {
            let at_or_below = heap
                .entries
                .iter()
                .filter(|entry| entry.list.best_evalue <= heap.cutoff)
                .count();
            assert_eq!(
                heap.entries.len() >= heap.threshold && heap.worst <= heap.cutoff,
                at_or_below >= heap.threshold,
                "round {round}"
            );
        }
        assert!(heaps[0].x_counts_ahead(0.002) && !heaps[0].x_counts_ahead(0.0021));
        for &e in &evalues {
            assert_eq!(
                CompoHeap::x_early_ahead(e, &heaps, &[0, 0]),
                CompoHeap::early(e, &heaps),
                "round {round} e {e}"
            );
            for ahead in [[1, 0], [0, 2], [1, 1], [3, 5], [6, 6]] {
                // The heaps after `ahead[q]` more inserts at or below the cutoff.
                let filled: Vec<CompoHeap> = shapes
                    .iter()
                    .zip(ahead)
                    .map(|((threshold, inserts), more)| {
                        let mut inserts = inserts.clone();
                        inserts.extend((0..more).map(|k| (evalues[random(4)], k as i32)));
                        build(*threshold, &inserts)
                    })
                    .collect();
                assert_eq!(
                    CompoHeap::x_early_ahead(e, &heaps, &ahead),
                    CompoHeap::early(e, &filled),
                    "round {round} e {e} ahead {ahead:?}"
                );
            }
        }
    }
}
