use super::*;

fn prelim(
    query: u32,
    minus: bool,
    qs: i32,
    qe: i32,
    ss: usize,
    score: i32,
    evalue: f64,
) -> PrelimHit {
    PrelimHit {
        context_idx: 2 * query + u32::from(minus),
        query_idx: query,
        query_frame: if minus { -1 } else { 1 },
        query_context_offset: 0,
        prelim_qs: qs,
        prelim_qe: qe,
        prelim_ss: ss,
        prelim_se: ss + (qe - qs) as usize,
        prelim_score: score,
        seed_qs: qs,
        seed_ss: ss,
        prelim_evalue: evalue,
    }
}

fn scores(lists: &PrelimHitLists, query: usize) -> Vec<(u32, Vec<i32>)> {
    lists[query]
        .as_ref()
        .map(|hit_list| {
            hit_list
                .hsplist_array
                .iter()
                .map(|list| (list.oid, list.hsps.iter().map(|h| h.prelim_score).collect()))
                .collect()
        })
        .unwrap_or_default()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:617-633: a minus HSP is folded to
// the plus strand, so a minus HSP over the same plus range as a better plus HSP is
// dominated by it, and one elsewhere is not.
#[test]
fn a_minus_hsp_is_folded_to_the_plus_strand() {
    let mut culling = PrelimCulling::new(&[1000], 10, 1);
    // Minus (700, 800) is plus (200, 300) of a query of 1000.
    culling.add(
        0,
        vec![
            prelim(0, false, 200, 300, 5, 100, 1e-20),
            prelim(0, true, 700, 800, 5, 60, 1e-10),
            prelim(0, true, 200, 300, 9, 50, 1e-9),
        ],
    );
    let lists = culling.finish();
    assert_eq!(scores(&lists, 0), vec![(0, vec![100, 50])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-581, c++/src/algo/blast/core/blast_hits.c:3125-3133: the
// hit list is made with `prelim_hitlist_size` and holds more lists than that, in the
// order of the first HSP of each subject, with the best e-value, worst e-value and low
// score set.
#[test]
fn the_hit_list_is_not_cut_and_keeps_the_rip_order() {
    let mut culling = PrelimCulling::new(&[1000, 500], 2, 5);
    for oid in 0..4usize {
        culling.add(
            oid,
            vec![
                prelim(
                    0,
                    false,
                    100 + 50 * oid as i32,
                    120 + 50 * oid as i32,
                    1,
                    40 - oid as i32,
                    1e-5,
                ),
                prelim(1, false, 10, 30, 1, 30, 1e-3),
                prelim(
                    0,
                    false,
                    100 + 50 * oid as i32,
                    110 + 50 * oid as i32,
                    2,
                    20,
                    1e-2,
                ),
            ],
        );
    }
    let lists = culling.finish();
    let q0 = lists[0].as_ref().unwrap();
    assert_eq!(q0.hsplist_max, 2);
    assert_eq!(q0.hsplist_count, 4);
    assert!(!q0.heapified);
    assert_eq!(q0.low_score, 37);
    assert!((q0.worst_evalue - 1e-5).abs() < 1e-30);
    // Limit 5 (+3 = 8) keeps everything. The root list holds the newest HSP first, so
    // the rip meets the last subject first; every list is score sorted.
    assert_eq!(
        scores(&lists, 0),
        vec![
            (3, vec![37, 20]),
            (2, vec![38, 20]),
            (1, vec![39, 20]),
            (0, vec![40, 20])
        ]
    );
    assert_eq!(lists[1].as_ref().unwrap().hsplist_count, 4);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:100-112: the lower OID wins a tie,
// and the subject must come in OID order for the later one to lose.
#[test]
fn identical_hsps_of_two_subjects_keep_the_lower_oid() {
    let mut culling = PrelimCulling::new(&[1000], 10, 1);
    culling.add(0, vec![prelim(0, false, 100, 200, 5, 80, 1e-10)]);
    culling.add(1, vec![prelim(0, false, 100, 200, 5, 80, 1e-10)]);
    let lists = culling.finish();
    assert_eq!(scores(&lists, 0), vec![(0, vec![80])]);
}
