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
