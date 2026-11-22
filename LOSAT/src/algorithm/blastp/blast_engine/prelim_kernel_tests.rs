use super::*;

fn hsp(query: (i32, i32), score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: query.0,
        query_end: query.1,
        subject_start: query.0,
        subject_end: query.1,
        gapped_query_start: query.0,
        gapped_subject_start: query.0,
        raw_score: score,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-999 (cited at
// `delete_failed_hsp_tests`): `Blast_HSPTest` of a BLASTP HSP fails for a negative length only.
#[test]
fn hsp_test_deletes_a_negative_length_and_keeps_zero_in_order() {
    let mut hsps = vec![
        hsp((10, 20), 30),
        hsp((12, 11), 25),
        hsp((5, 5), 20),
        hsp((40, 39), 15),
        hsp((7, 9), 10),
    ];
    delete_failed_hsp_tests(&mut hsps);
    let kept: Vec<(i32, i32, i32)> = hsps
        .iter()
        .map(|hsp| (hsp.query_start, hsp.query_end, hsp.raw_score))
        .collect();
    assert_eq!(kept, [(10, 20, 30), (5, 5, 20), (7, 9, 10)]);
}
