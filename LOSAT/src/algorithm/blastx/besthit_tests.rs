use super::*;
use crate::algorithm::blastx::preliminary::PreliminaryHsp;

/// The query information of queries of `dna` letters each, both strands searched or one
/// (`plus`, `minus`), as `query_setup::prepare_queries` lays the contexts out.
fn frames(dna: &[usize], strand: &str) -> QueryFrames {
    let mut lengths = Vec::new();
    let mut minus = Vec::new();
    for &letters in dna {
        for i in 0..6 {
            let (length, enabled) =
                crate::algorithm::blastx::query_setup::context_frame(letters, i, strand);
            lengths.push(if enabled { length as i32 } else { 0 });
            minus.push(i >= 3);
        }
    }
    QueryFrames { lengths, minus }
}

fn hsp(context: usize, q_start: i32, q_end: i32, score: i32, evalue: f64) -> Hsp {
    Hsp {
        hsp: PreliminaryHsp {
            context,
            frame: [1, 2, 3, -1, -2, -3][context % 6],
            score,
            q_start,
            q_end,
            q_gapped_start: q_start,
            s_start: 5,
            s_end: 5 + (q_end - q_start),
            s_gapped_start: 5,
        },
        num: 1,
        evalue,
        bit_score: 0.0,
        identity: 0,
        positive: 0,
        edit_script: Vec::new(),
        composition_method: 0,
    }
}

fn params() -> BestHitParams {
    BestHitParams {
        overhang: 0.1,
        score_edge: 0.1,
    }
}

fn list(oid: usize, hsps: Vec<Hsp>) -> HspList {
    let mut list = HspList {
        query_index: hsps[0].hsp.context / 6,
        oid,
        hsps,
        best_evalue: 0.0,
    };
    list.refresh();
    list
}

/// Each list of a hit list: its OID and its HSPs' scores, in order.
fn scores(hit_list: &HitList) -> Vec<(usize, Vec<i32>)> {
    hit_list
        .lists
        .iter()
        .map(|l| (l.oid, l.hsps.iter().map(|h| h.hsp.score).collect()))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266: a plus frame's
// HSP starts at its own offset in the frame.
#[test]
fn a_plus_frame_key_is_the_offset_in_the_frame() {
    let (query, key) = frames(&[100], "both").key(&hsp(1, 10, 40, 77, 1e-9));
    assert_eq!(query, 0);
    assert_eq!(
        key,
        BestHitKey {
            begin: 10,
            len: 30,
            score: 77,
            evalue: 1e-9,
        }
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266: a minus frame's
// HSP starts at the query's length (`BlastQueryInfoGetQueryLength`: 2 plus the three plus frames'
// lengths, 33 + 33 + 32 for 100 letters) less its end in the frame.
#[test]
fn a_minus_frame_key_folds_on_the_dna_length() {
    let (query, key) = frames(&[100], "both").key(&hsp(4, 10, 40, 77, 1e-9));
    assert_eq!(query, 0);
    assert_eq!((key.begin, key.len), (100 - 40, 30));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:141-161: with only the
// minus strand searched the plus contexts have length 0 and the length is the minus frames'
// (33 + 33 + 33 + 2 for 101 letters); with only the plus strand it is the plus frames'.
#[test]
fn one_strand_takes_the_length_of_its_frames() {
    let (_, key) = frames(&[101], "minus").key(&hsp(3, 0, 20, 50, 1e-3));
    assert_eq!(key.begin, 101 - 20);
    let (_, key) = frames(&[101], "plus").key(&hsp(2, 0, 20, 50, 1e-3));
    assert_eq!(key.begin, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:40-49: the query of an
// HSP is its context over the six frames, and its length is that query's own.
#[test]
fn the_second_query_reads_its_own_frames() {
    let (query, key) = frames(&[100, 61], "both").key(&hsp(9, 5, 15, 40, 1e-2));
    assert_eq!(query, 1);
    assert_eq!(key.begin, 61 - 15);
}

/// Three HSPs of one query frame where the order of the saves decides (a chain): `a` holds `b`
/// (within its overhang, with a better e-value and a score per letter better by the edge), `b`
/// holds `c`, and `a` does not hold `c`. Saved `a`, `b`, `c`, `b` is dropped and `c` kept;
/// saved `c`, `b`, `a`, `b` removes `c` and `a` removes `b`.
fn chain() -> [Hsp; 3] {
    [
        hsp(0, 100, 200, 500, 1e-30),
        hsp(0, 95, 195, 400, 1e-20),
        hsp(0, 86, 186, 300, 1e-10),
    ]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539: the pipe saves
// the HSPs of a list in the list's order (here the traceback's, not sorted by e-value); it does
// not sort them.
#[test]
fn the_pipe_keeps_the_order_of_a_lists_hsps() {
    let piped = |order: [usize; 3]| {
        let [a, b, c] = chain();
        let mut all = [Some(a), Some(b), Some(c)];
        let hsps: Vec<Hsp> = order.iter().map(|&i| all[i].take().unwrap()).collect();
        let mut queries = vec![HitList::new(500)];
        queries[0].update(list(0, hsps));
        traceback_pipe(&mut queries, &frames(&[900], "both"), params(), 500);
        scores(&queries[0])
    };
    assert_eq!(piped([2, 1, 0]), vec![(0, vec![500])]);
    assert_eq!(piped([0, 1, 2]), vec![(0, vec![500, 300])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400: the pipe sorts a
// query's lists by their best e-values before it saves their HSPs (`a`'s list first, so `c` is
// kept), and the final hit list is in e-value order.
#[test]
fn the_pipe_saves_the_lists_in_evalue_order() {
    let [a, b, c] = chain();
    let mut queries = vec![HitList::new(500)];
    queries[0].update(list(4, vec![c]));
    queries[0].update(list(5, vec![b]));
    queries[0].update(list(7, vec![a]));
    traceback_pipe(&mut queries, &frames(&[900], "both"), params(), 500);
    assert_eq!(scores(&queries[0]), vec![(7, vec![500]), (4, vec![300])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:194-231: the writer's
// final step makes each query's hit list apart: one list per subject sorted by score, the lists
// by e-value, a query without HSPs empty.
#[test]
fn the_writer_makes_each_querys_hit_list() {
    let mut writer = BestHitWriter::new(params(), frames(&[300, 300, 300], "both"), 500);
    writer.write(
        2,
        vec![hsp(0, 100, 130, 60, 1e-4), hsp(12, 0, 30, 90, 1e-8)],
    );
    writer.write(3, vec![hsp(0, 0, 30, 50, 1e-3), hsp(0, 200, 240, 80, 1e-6)]);
    let queries = writer.finalize(3, 500);
    assert_eq!(queries.len(), 3);
    assert_eq!(scores(&queries[0]), vec![(3, vec![80, 50]), (2, vec![60])]);
    assert!(queries[1].lists.is_empty());
    assert_eq!(scores(&queries[2]), vec![(2, vec![90])]);
    assert!(queries[0].lists.iter().all(|l| l.query_index == 0));
    assert_eq!(queries[0].lists[0].best_evalue, 1e-6);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298: the final hit list
// holds `prelim_hitlist_size` lists, those of the best e-values.
#[test]
fn the_final_hit_list_keeps_the_best_lists() {
    let mut writer = BestHitWriter::new(params(), frames(&[900], "both"), 2);
    for (oid, start) in [(0usize, 0), (1, 100), (2, 200)] {
        let evalue = [1e-3, 1e-9, 1e-6][oid];
        writer.write(oid, vec![hsp(0, start, start + 40, 70, evalue)]);
    }
    let queries = writer.finalize(1, 2);
    let mut oids: Vec<usize> = queries[0].lists.iter().map(|l| l.oid).collect();
    oids.sort_unstable();
    assert_eq!(oids, vec![1, 2]);
}
