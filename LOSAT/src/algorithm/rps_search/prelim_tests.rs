use super::*;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;

fn prelim_list(oid: u32, score: i32, e_value: f64) -> BlastpHspList<XPrelimHsp> {
    BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![XPrelimHsp {
            hsp: BlastpPreliminaryHsp {
                query_start: 0,
                query_end: 10,
                subject_start: 0,
                subject_end: 10,
                gapped_query_start: 0,
                gapped_subject_start: 0,
                raw_score: score,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 100,
                q_idx: 0,
                s_idx: oid,
            },
            e_value,
        }],
        best_evalue: f64::MAX,
    }
}

// blast_engine.c:874-899, blast_hits.c (`Blast_HitListUpdate`, `s_EvalueCompareHSPLists`):
// without the Gumbel block every preliminary HSP has the e-value 0, so the lists are ordered
// by the first HSP's score (down), then the ordinal id (down), which is the order of
// rpsblast's composition-based redo's reads: the last list read is the lowest, and its
// profile's K is the batch's (`stream_lists`, `RpsTraceback::trace_query_redo`).
#[test]
fn lists_of_e_value_zero_are_ordered_by_score_then_ordinal_id() {
    let mut hit_list = BlastpHitList::new(10);
    for (oid, score) in [(0, 50), (1, 80), (2, 80), (3, 65)] {
        hit_list.update(prelim_list(oid, score, 0.0));
    }
    hit_list.sort_by_evalue();
    let oids: Vec<u32> = hit_list.hsplist_array[..hit_list.hsplist_count]
        .iter()
        .map(|list| list.oid)
        .collect();
    assert_eq!(oids, vec![2, 1, 3, 0]);
    // With e-values the best e-value decides first.
    let mut hit_list = BlastpHitList::new(10);
    for (oid, score, e_value) in [(0, 50, 1e-5), (1, 80, 1e-2), (2, 80, 1e-9)] {
        hit_list.update(prelim_list(oid, score, e_value));
    }
    hit_list.sort_by_evalue();
    let oids: Vec<u32> = hit_list.hsplist_array[..hit_list.hsplist_count]
        .iter()
        .map(|list| list.oid)
        .collect();
    assert_eq!(oids, vec![2, 0, 1]);
}

fn kept_oids(lists: &[(u32, i32, f64)]) -> Vec<u32> {
    let size = crate::algorithm::blastn::hsp::get_prelim_hitlist_size(1, false, true);
    assert_eq!(size, 10);
    let mut hit_list = BlastpHitList::new(size);
    for &(oid, score, e_value) in lists {
        hit_list.update(prelim_list(oid, score, e_value));
    }
    hit_list.sort_by_evalue();
    hit_list.hsplist_array[..hit_list.hsplist_count]
        .iter()
        .map(|list| list.oid)
        .collect()
}

// The prelim hit list holds `min(max(2*size, 10), size + 50)` lists, 10 for `-max_target_seqs`
// 1 to 5 (blast_hits.c:68), and a list past that capacity goes through
// `Blast_HitListUpdate` (blast_hits.c:3243-3297) into the heap whose root is the worst list:
// with `s_EvalueCompareHSPLists` (blast_hits.c:3078-3106) ordering by best E-value (up), then
// the first HSP's score (down), then `BLAST_CMP(h2->oid, h1->oid)`: of two lists equal in
// E-value and score the one of the higher ordinal id is the better, so the root is the lowest
// id of the worst tie and a new list of a higher id replaces it ("newer hits with score and
// e-value both equal to the current worst will be saved"), one of a lower id is discarded.
// The kept lists are therefore the best 10 of the total order whatever the insertion order,
// and sorted they read best first.
#[test]
fn a_full_prelim_list_of_equal_e_values_keeps_the_best_ten_of_the_total_order() {
    let ascending: Vec<u32> = (0..18).collect();
    let descending: Vec<u32> = (0..18).rev().collect();
    let shuffled = [
        9u32, 3, 17, 0, 12, 5, 16, 1, 8, 14, 2, 11, 6, 15, 4, 13, 7, 10,
    ];
    for order in [ascending.as_slice(), descending.as_slice(), &shuffled] {
        // Equal e-value 0 and equal score: the ten highest ordinal ids, highest first.
        let equal: Vec<(u32, i32, f64)> = order.iter().map(|&oid| (oid, 100, 0.0)).collect();
        assert_eq!(
            kept_oids(&equal),
            vec![17, 16, 15, 14, 13, 12, 11, 10, 9, 8],
            "order {order:?}"
        );
        // Even ids score 90 and odd ids 70: the nine of score 90 (high id first), then
        // the highest id of score 70.
        let mixed: Vec<(u32, i32, f64)> = order
            .iter()
            .map(|&oid| (oid, if oid % 2 == 0 { 90 } else { 70 }, 0.0))
            .collect();
        assert_eq!(
            kept_oids(&mixed),
            vec![16, 14, 12, 10, 8, 6, 4, 2, 0, 17],
            "order {order:?}"
        );
        // The E-value decides before the score and the id: the ten lists of E-value 0 and
        // score 10 (ids 8 to 17) beat the eight of E-value 1e-30 and score 100 (ids 0 to 7).
        let with_e: Vec<(u32, i32, f64)> = order
            .iter()
            .map(|&oid| {
                if oid < 8 {
                    (oid, 100, 1e-30)
                } else {
                    (oid, 10, 0.0)
                }
            })
            .collect();
        let expected = vec![17, 16, 15, 14, 13, 12, 11, 10, 9, 8];
        assert_eq!(kept_oids(&with_e), expected, "order {order:?}");
    }
}

#[test]
fn the_context_of_an_offset_is_the_profile_that_holds_it() {
    let offsets = [0, 399, 1130, 3482, 3772, 4621, 5035];
    assert_eq!(context_of(&offsets, 0), 0);
    assert_eq!(context_of(&offsets, 398), 0);
    assert_eq!(context_of(&offsets, 399), 1);
    assert_eq!(context_of(&offsets, 5033), 5);
    assert_eq!(context_of(&offsets, 5034), 5);
}
