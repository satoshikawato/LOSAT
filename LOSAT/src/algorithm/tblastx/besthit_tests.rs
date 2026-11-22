use super::*;
use crate::algorithm::tblastx::hit_collector::HitCollector;
use crate::common::Hit;

fn hsp(
    s_idx: u32,
    frame: i32,
    raw_score: i32,
    e_value: f64,
    offset: usize,
    len: usize,
) -> TblastxHsp {
    TblastxHsp {
        hit: Hit {
            identity: 0.0,
            length: len,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 30,
            s_start: 1,
            s_end: 30,
            e_value,
            bit_score: raw_score as f64 / 2.0,
            num_ident: 0,
            query_frame: frame,
            query_length: 0,
            q_idx: 0,
            s_idx,
            raw_score,
            sort_query_offset: offset,
            sort_query_end: offset + len,
            sort_subject_offset: offset,
            sort_subject_end: offset + len,
            has_sort_offsets: true,
            gap_info: None,
            num_positives: 0,
        },
        subject_frame: 1,
        num: 1,
    }
}

fn params(overhang: f64, score_edge: f64) -> BestHitParams {
    BestHitParams {
        overhang,
        score_edge,
    }
}

/// The kept HSPs as (subject, score) in report order.
fn run(
    hitlist_size: usize,
    query_length: usize,
    subjects: Vec<Vec<TblastxHsp>>,
) -> Vec<(u32, i32)> {
    let mut collector =
        HitCollector::with_options(hitlist_size, 0, 0, 0.0, vec![query_length], 10.0)
            .with_best_hit(Some(params(0.1, 0.1)));
    for (oid, hits) in subjects.into_iter().enumerate() {
        collector.push(oid, hits);
    }
    collector
        .finish()
        .iter()
        .map(|hsp| (hsp.hit.s_idx, hsp.hit.raw_score))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266,
// c++/src/algo/blast/core/blast_query_info.c:141-161: a minus frame folds its frame range over
// the query's length, 2 plus its three frames (a DNA length; the offsets are residues).
#[test]
fn a_minus_frame_key_folds_over_the_translated_query_length() {
    // 100 nucleotides: frames of 33, 33 and 32 residues, so the length is 100.
    let qlen = translated_query_length(&tblastx_context_lengths(100));
    assert_eq!(qlen, 100);
    let key = best_hit_key(&hsp(0, -2, 40, 1e-5, 5, 10), qlen);
    assert_eq!((key.begin, key.len, key.score), (85, 10, 40));
    let key = best_hit_key(&hsp(0, 3, 40, 1e-5, 5, 10), qlen);
    assert_eq!((key.begin, key.len), (5, 10));
}

// A subject whose HSP lies within a better HSP of an earlier subject is dropped, a subject
// elsewhere on the query is kept.
#[test]
fn a_contained_worse_subject_is_dropped() {
    let kept = run(
        500,
        300,
        vec![
            vec![hsp(0, 1, 100, 1e-20, 10, 50)],
            vec![hsp(1, 1, 30, 1e-3, 20, 20)],
            vec![hsp(2, 1, 40, 1e-5, 70, 20)],
        ],
    );
    assert_eq!(kept, vec![(0, 100), (2, 40)]);
}

// The same range on the other strand is another place: a minus-frame HSP is not held by a plus
// one whose frame range it shares unless the fold puts it there.
#[test]
fn a_minus_frame_hsp_is_compared_where_it_folds() {
    // qlen 300: frame range 10..60 on the minus strand folds to 240..290.
    let kept = run(
        500,
        300,
        vec![
            vec![hsp(0, 1, 100, 1e-20, 10, 50)],
            vec![hsp(1, -1, 30, 1e-3, 10, 50)],
            vec![hsp(2, -1, 10, 1e-2, 15, 30)],
        ],
    );
    assert_eq!(kept, vec![(0, 100), (1, 30)]);
}

// `-max_target_seqs` keeps the subjects of the best e-values (the writer's hit list and the
// pipe's are of that size).
#[test]
fn the_hit_list_keeps_the_best_subjects() {
    let kept = run(
        2,
        3000,
        vec![
            vec![hsp(0, 1, 40, 1e-5, 0, 20)],
            vec![hsp(1, 1, 60, 1e-9, 100, 20)],
            vec![hsp(2, 1, 50, 1e-7, 200, 20)],
        ],
    );
    assert_eq!(kept, vec![(1, 60), (2, 50)]);
}
