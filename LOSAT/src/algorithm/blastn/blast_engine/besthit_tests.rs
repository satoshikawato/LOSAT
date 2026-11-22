use super::*;

fn prelim(query: u32, minus: bool, qs: i32, qe: i32, score: i32, evalue: f64) -> PrelimHit {
    PrelimHit {
        context_idx: 2 * query + u32::from(minus),
        query_idx: query,
        query_frame: if minus { -1 } else { 1 },
        query_context_offset: 0,
        prelim_qs: qs,
        prelim_qe: qe,
        prelim_ss: 5,
        prelim_se: 5 + (qe - qs) as usize,
        prelim_score: score,
        seed_qs: qs,
        seed_ss: 5,
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

fn params() -> BestHitParams {
    BestHitParams {
        overhang: 0.1,
        score_edge: 0.1,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266: a minus
// HSP is folded to the plus strand, so a minus HSP over the plus range of a better plus HSP
// is dropped, and one elsewhere is not.
#[test]
fn a_minus_hsp_is_folded_to_the_plus_strand() {
    let mut writer = PrelimBestHit::new(&[1000], 10, params());
    // Minus (700, 800) is plus (200, 300) of a query of 1000.
    writer.add(
        0,
        vec![
            prelim(0, false, 200, 300, 200, 1e-40),
            prelim(0, true, 700, 800, 60, 1e-10),
            prelim(0, true, 200, 300, 50, 1e-9),
        ],
    );
    let lists = writer.finish();
    assert_eq!(scores(&lists, 0), vec![(0, vec![200, 50])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:207-221: the lists
// are sorted by e-value and saved into a hit list of `prelim_hitlist_size`, which keeps the
// best ones.
#[test]
fn the_hit_list_keeps_the_best_lists() {
    let mut writer = PrelimBestHit::new(&[1000], 2, params());
    for oid in 0..4usize {
        let start = 100 * oid as i32;
        writer.add(
            oid,
            vec![prelim(
                0,
                false,
                start,
                start + 50,
                40 + oid as i32,
                1e-5 / (oid + 1) as f64,
            )],
        );
    }
    let lists = writer.finish();
    let hit_list = lists[0].as_ref().unwrap();
    assert_eq!(hit_list.hsplist_count, 2);
    let mut kept = scores(&lists, 0);
    kept.sort();
    assert_eq!(kept, vec![(2, vec![42]), (3, vec![43])]);
}
