use super::*;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;

fn prelim(query: u32, qlen: usize, qs: i32, qe: i32, ss: i32, score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: qs,
        query_end: qe,
        subject_start: ss,
        subject_end: ss + (qe - qs),
        gapped_query_start: qs,
        gapped_subject_start: ss,
        raw_score: score,
        query_context: query as i32,
        query_frame: 0,
        subject_frame: 0,
        query_length: qlen,
        q_idx: query,
        s_idx: 0,
    }
}

fn full_hsp(hsp: &BlastpPreliminaryHsp, oid: u32, e_value: f64) -> BlastpHsp {
    let mut full = super::super::blastp_hsp_from_preliminary(hsp, e_value);
    full.s_idx = oid;
    full
}

/// The subject result of `oid` with the HSPs of one query.
fn subject(oid: usize, query: usize, hsps: &[(BlastpPreliminaryHsp, f64)]) -> BlastpSubjectResult {
    let mut result = BlastpSubjectResult::new(oid);
    result.x_compact_hsp_lists.push((
        query,
        BlastpHspList {
            oid: oid as u32,
            query_index: query as u32,
            hsps: hsps
                .iter()
                .map(|(h, e)| XPrelimHsp {
                    hsp: BlastpPreliminaryHsp {
                        s_idx: oid as u32,
                        ..*h
                    },
                    e_value: *e,
                })
                .collect(),
            best_evalue: 0.0,
        },
    ));
    result
}

fn finished(culling: PrelimCulling, queries: usize) -> Vec<Option<BlastpHitList<XPrelimHsp>>> {
    let mut compact: Vec<Option<BlastpHitList<XPrelimHsp>>> = (0..queries).map(|_| None).collect();
    culling.finish(&mut compact);
    compact
}

fn scores(lists: &[Option<BlastpHitList<XPrelimHsp>>], query: usize) -> Vec<(u32, Vec<i32>)> {
    lists[query]
        .as_ref()
        .map(|hit_list| {
            hit_list
                .hsplist_array
                .iter()
                .map(|l| (l.oid, l.hsps.iter().map(|h| h.hsp.raw_score).collect()))
                .collect()
        })
        .unwrap_or_default()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:612-633: the key of a protein HSP is its
// own query range (0-based start, exclusive end), its score and its 0-based subject offset;
// `qlen` is the HSP's query length.
#[test]
fn the_key_is_the_context_relative_range_score_and_subject_offset() {
    let hsp = full_hsp(&prelim(0, 321, 10, 60, 7, 55), 3, 1e-5);
    let (qlen, key) = cull_key(&hsp);
    assert_eq!(qlen, 321);
    assert_eq!(
        key,
        CullKey {
            begin: 10,
            end: 60,
            score: 55,
            subject_offset: 7
        }
    );
    let compact = XPrelimHsp {
        hsp: prelim(0, 321, 10, 60, 7, 55),
        e_value: 1e-5,
    };
    assert_eq!(cull_key(&compact), (qlen, key));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-581 and c++/src/algo/blast/core/blast_hits.c:3125-3133: the hit list is made with
// `prelim_hitlist_size` and holds more lists than that, in the order of the first HSP of
// each subject in the rip, with the best e-value, worst e-value and low score set, and not
// heapified.
#[test]
fn the_hit_list_is_not_cut_and_keeps_the_rip_order() {
    let mut culling = PrelimCulling::new(2, 5);
    for oid in 0..4usize {
        let a = prelim(
            0,
            1000,
            100 + 50 * oid as i32,
            120 + 50 * oid as i32,
            1,
            40 - oid as i32,
        );
        let b = prelim(0, 1000, 100 + 50 * oid as i32, 110 + 50 * oid as i32, 2, 20);
        let c = prelim(1, 500, 10, 30, 1, 30);
        let mut first = subject(oid, 0, &[(a, 1e-5), (b, 1e-2)]);
        let second = subject(oid, 1, &[(c, 1e-3)]);
        first.x_compact_hsp_lists.extend(second.x_compact_hsp_lists);
        culling.add(&mut first);
        assert!(first.x_compact_hsp_lists.is_empty());
    }
    let lists = finished(culling, 2);
    let q0 = lists[0].as_ref().unwrap();
    assert_eq!(q0.hsplist_max, 2);
    assert_eq!(q0.hsplist_count, 4);
    assert!(!q0.heapified);
    assert_eq!(q0.low_score, 37);
    assert!((q0.worst_evalue - 1e-5).abs() < 1e-30);
    // Limit 5 (+3 = 8) keeps everything. The root list holds the newest HSP first, so the
    // rip meets the last subject first; every list is score sorted.
    let expected = vec![
        (3, vec![37, 20]),
        (2, vec![38, 20]),
        (1, vec![39, 20]),
        (0, vec![40, 20]),
    ];
    assert_eq!(scores(&lists, 0), expected);
    assert_eq!(lists[1].as_ref().unwrap().hsplist_count, 4);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:100-112: the lower OID wins a tie, and the subjects must come in
// OID order for the later one to lose.
#[test]
fn identical_hsps_of_two_subjects_keep_the_lower_oid() {
    let mut culling = PrelimCulling::new(10, 1);
    let h = prelim(0, 1000, 100, 200, 5, 80);
    culling.add(&mut subject(0, 0, &[(h, 1e-10)]));
    culling.add(&mut subject(1, 0, &[(h, 1e-10)]));
    let lists = finished(culling, 1);
    assert_eq!(scores(&lists, 0), vec![(0, vec![80])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:602-644: the tree of a query is made with the length of its HSPs,
// and a chunk call's parts are queries of their own: a better HSP of query 1 does not
// dominate an HSP of query 0 over the same range.
#[test]
fn each_query_has_its_own_tree() {
    let mut culling = PrelimCulling::new(10, 1);
    culling.add(&mut subject(0, 0, &[(prelim(0, 400, 10, 90, 0, 50), 1e-3)]));
    culling.add(&mut subject(1, 1, &[(prelim(1, 300, 10, 90, 0, 90), 1e-9)]));
    culling.add(&mut subject(2, 0, &[(prelim(0, 400, 10, 90, 0, 60), 1e-5)]));
    let lists = finished(culling, 2);
    // Limit 1: score 60 dominates score 50 in query 0.
    assert_eq!(scores(&lists, 0), vec![(2, vec![60])]);
    assert_eq!(scores(&lists, 1), vec![(1, vec![90])]);
}

fn pipe_hsp(oid: u32, qs: usize, qe: usize, score: i32, e_value: f64) -> BlastpHsp {
    let hsp = prelim(0, 1000, qs as i32 - 1, qe as i32, 0, score);
    full_hsp(&hsp, oid, e_value)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:699-733: the lists are sorted by e-value (best first, then the
// first score, the higher OID first), the HSPs are saved in that order into a tree with the
// limit itself (no +3), and the kept HSPs are made into lists again, each sorted by score.
#[test]
fn the_pipe_saves_in_e_value_order_with_the_plain_limit() {
    let mut hit_list = BlastpHitList::new(10);
    // OID 1 has the better e-value, so it is saved first; OID 0's weaker HSP over the same
    // range is dominated by it (limit 1).
    for (oid, score, e_value) in [(0u32, 60, 1e-6), (1, 90, 1e-12)] {
        hit_list.hsplist_array.push(BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![pipe_hsp(oid, 11, 90, score, e_value)],
            best_evalue: e_value,
        });
    }
    hit_list.hsplist_count = 2;
    cull_traceback(&mut hit_list, 1);
    assert_eq!(hit_list.hsplist_count, 1);
    assert_eq!(hit_list.hsplist_array[0].oid, 1);
    assert_eq!(hit_list.hsplist_max, 10);
    assert!(!hit_list.heapified);
    assert_eq!(hit_list.low_score, 90);
}

// With the limit 2 the second dominated HSP stays; the lists keep the rip order and each is
// sorted by score.
#[test]
fn the_pipe_keeps_the_hsps_that_fewer_than_the_limit_dominate() {
    let mut hit_list = BlastpHitList::new(10);
    for (oid, score, e_value) in [(0u32, 60, 1e-6), (1, 90, 1e-12), (2, 70, 1e-8)] {
        hit_list.hsplist_array.push(BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![pipe_hsp(oid, 11, 90, score, e_value)],
            best_evalue: e_value,
        });
    }
    hit_list.hsplist_count = 3;
    cull_traceback(&mut hit_list, 2);
    let mut kept: Vec<u32> = hit_list.hsplist_array.iter().map(|l| l.oid).collect();
    kept.sort_unstable();
    // Saved in the order 1 (1e-12), 2 (1e-8), 0 (1e-6): OID 0 is dominated by both.
    assert_eq!(kept, vec![1, 2]);
    assert_eq!(hit_list.hsplist_count, 2);
}
