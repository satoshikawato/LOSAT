use super::*;
use crate::common::Hit;

fn hsp_of(
    q_idx: u32,
    s_idx: u32,
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
            query_frame: 1,
            query_length: 0,
            q_idx,
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

fn hsp(q_idx: u32, s_idx: u32, raw_score: i32, e_value: f64, offset: usize) -> TblastxHsp {
    hsp_of(q_idx, s_idx, raw_score, e_value, offset, 10)
}

fn key(hsps: &[TblastxHsp]) -> Vec<(u32, u32, i32, f64, usize)> {
    hsps.iter()
        .map(|hsp| {
            (
                hsp.hit.q_idx,
                hsp.hit.s_idx,
                hsp.hit.raw_score,
                hsp.hit.e_value,
                hsp.hit.sort_query_offset,
            )
        })
        .collect()
}

const THRESHOLD: f64 = 1e-5;

// The final HSPs of `subjects` subjects for two queries, as the subject loop sends them: one
// list per subject, the two queries' HSPs interleaved, not in score order. Every seventh
// subject has none. Most best e-values tie (0.0, as 18S rRNA against its own database), and
// in every list the HSP with the highest score has a worse e-value (9e-6, under the
// threshold) than the best one, so the comparison of two lists reads a different first HSP
// before and after the hit list fills (score order, then e-value order). Query 0's lists
// have three or four HSPs (`-max_hsps` cuts them), two with equal score and e-value (a tie
// the stable sort keeps in list order), and one of 30 residues (`-qcov_hsp_perc` keeps it
// alone). Every eleventh subject has an HSP above the threshold.
fn subject_lists(subjects: usize) -> Vec<(usize, Vec<TblastxHsp>)> {
    (0..subjects)
        .map(|s| {
            let oid = s as u32;
            if s % 7 == 6 {
                return (s, Vec::new());
            }
            let best = if s % 3 == 0 {
                0.0
            } else {
                1e-6 * (1 + s % 5) as f64
            };
            let mut hsps = vec![
                hsp(1, oid, 130 + (s % 4) as i32, 9e-6, 40),
                hsp(0, oid, 120 + (s % 4) as i32, 9e-6, 0),
                hsp_of(0, oid, 100 + (s % 6) as i32, best, 20, 30),
                hsp(1, oid, 110 + (s % 6) as i32, best * 0.5, 60),
                hsp(0, oid, 90, best * 2.0, 80),
            ];
            if s % 2 == 1 {
                hsps.push(hsp(0, oid, 100 + (s % 6) as i32, best, 10));
            }
            if s % 5 == 0 {
                hsps.push(hsp(1, oid, 105, best * 0.5, 70));
            }
            if s % 11 == 10 {
                hsps.push(hsp(0, oid, 50, 1e-2, 90));
            }
            (s, hsps)
        })
        .collect()
}

fn non_empty(lists: &[(usize, Vec<TblastxHsp>)]) -> usize {
    lists.iter().filter(|(_, hsps)| !hsps.is_empty()).count()
}

// `final_hit_order` on the HSPs of every subject at once (the lists of a run), after the
// `-evalue` reap.
fn reference(
    lists: &[(usize, Vec<TblastxHsp>)],
    hitlist_size: usize,
    filters: &HitFilters<'_>,
) -> Vec<TblastxHsp> {
    #[expect(
        clippy::neg_cmp_op_on_partial_ord,
        reason = "the same filter as NCBI's evalue > cutoff removal (blast_hits.c:1996); `<=` would drop NaN values"
    )]
    let hits: Vec<TblastxHsp> = lists
        .iter()
        .flat_map(|(_, hsps)| hsps.iter().cloned())
        .filter(|hsp| !(hsp.hit.e_value > THRESHOLD))
        .collect();
    report::final_hit_order(hits, hitlist_size, filters)
}

// A permutation of 0..n (7919 is a prime above every n used here).
fn shuffled(n: usize) -> Vec<usize> {
    (0..n).map(|i| i * 7919 % n).collect()
}

fn collect(
    lists: &[(usize, Vec<TblastxHsp>)],
    order: &[usize],
    mut collector: HitCollector,
) -> Vec<TblastxHsp> {
    for &s in order {
        collector.push(s, lists[s].1.clone());
    }
    collector.finish()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:143-150 and c++/src/algo/blast/core/blast_hits.c:3243-3297
// (Blast_HitListUpdate): the subjects enter the hit lists in OID order whatever order their
// searches end in, and the lists are those of the whole batch's HSPs ordered at once
// (`final_hit_order`): with hit lists that never fill (sizes at and above the number of
// subjects with HSPs) and with ones that fill and heapify (sizes 1, 2 and one below).
#[test]
fn subjects_in_any_order_give_final_hit_orders_lists() {
    for subjects in [12, 520] {
        let lists = subject_lists(subjects);
        let full = non_empty(&lists);
        assert!(full >= 10 && full < subjects);
        for size in [1, 2, full - 1, full, full + 1] {
            let want = reference(&lists, size, &HitFilters::default());
            assert!(
                want.len() > 2 * size.min(full),
                "{subjects} subjects, size {size}"
            );
            for order in [(0..subjects).collect::<Vec<_>>(), shuffled(subjects)] {
                let collector = HitCollector::with_options(size, 0, 0, 0.0, Vec::new(), THRESHOLD);
                let got = collect(&lists, &order, collector);
                assert_eq!(key(&got), key(&want), "{subjects} subjects, size {size}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3078-3105 (s_EvalueCompareHSPLists) and 3266-3278: pinned by hand. Subject
// 0's best e-value is on its lower-scoring HSP. Below the size the lists keep their score
// order and subject 0 compares with 120; when the third subject fills a hit list of two, the
// lists are sorted by e-value, subject 0 compares with 100 and is the worst, so subject 2
// replaces it. The subjects come in the order 1, 2, 0.
#[test]
fn hit_list_keeps_ncbis_subjects_on_ties_in_any_order() {
    let lists = vec![
        (0, vec![hsp(0, 0, 120, 1e-3, 0), hsp(0, 0, 100, 1e-5, 20)]),
        (1, vec![hsp(0, 1, 110, 1e-5, 0)]),
        (2, vec![hsp(0, 2, 115, 1e-5, 0)]),
    ];
    let order = |size, subjects: &[usize]| -> Vec<(u32, i32)> {
        let collector = HitCollector::with_options(size, 0, 0, 0.0, Vec::new(), 1.0);
        collect(&lists, subjects, collector)
            .iter()
            .map(|hsp| (hsp.hit.s_idx, hsp.hit.raw_score))
            .collect()
    };
    assert_eq!(order(2, &[1, 0]), vec![(0, 100), (0, 120), (1, 110)]);
    assert_eq!(order(1, &[1, 0]), vec![(1, 110)]);
    assert_eq!(
        order(3, &[1, 2, 0]),
        vec![(0, 100), (0, 120), (2, 115), (1, 110)]
    );
    assert_eq!(order(2, &[1, 2, 0]), vec![(2, 115), (1, 110)]);
    assert_eq!(order(2, &[0, 1, 2]), vec![(2, 115), (1, 110)]);
}

// `-max_hsps` cuts each kept list (query 0's lists have three or four HSPs), `-qcov_hsp_perc`
// reaps a subject's HSPs before its list enters the hit list (query 0 keeps its HSP of 30
// residues alone, query 1 keeps every HSP): as `final_hit_order` applies them, with hit lists
// that fill.
#[test]
fn filters_apply_as_in_final_hit_order() {
    let subjects = 503;
    let lists = subject_lists(subjects);
    let full = non_empty(&lists);
    let size = full - 20;
    let query_lengths = vec![300usize, 90];
    let filters = HitFilters {
        max_hsps: 2,
        query_cov_hsp_perc: 25.0,
        query_lengths: &query_lengths,
        subject_besthit: false,
    };
    let want = reference(&lists, size, &filters);
    let q0: Vec<&TblastxHsp> = want.iter().filter(|hsp| hsp.hit.q_idx == 0).collect();
    assert_eq!(q0.len(), size, "one HSP of 30 residues per kept subject");
    assert!(q0.iter().all(|hsp| hsp.hit.length == 30));
    let q1 = want.iter().filter(|hsp| hsp.hit.q_idx == 1).count();
    assert!(
        q1 > size && q1 <= 2 * size,
        "{q1} HSPs of query 1 after -max_hsps 2"
    );
    let collector = HitCollector::with_options(size, 0, 2, 25.0, query_lengths.clone(), THRESHOLD);
    assert_eq!(
        key(&collect(&lists, &shuffled(subjects), collector)),
        key(&want)
    );
}

// `-culling_limit`: the batch's HSPs are culled at its end (`culled_hit_order`).
#[test]
fn culling_limit_culls_the_batch_at_its_end() {
    let subjects = 60;
    let lists = subject_lists(subjects);
    let query_lengths = vec![300usize, 90];
    #[expect(
        clippy::neg_cmp_op_on_partial_ord,
        reason = "the same filter as NCBI's evalue > cutoff removal (blast_hits.c:1996); `<=` would drop NaN values"
    )]
    let hits: Vec<TblastxHsp> = lists
        .iter()
        .flat_map(|(_, hsps)| hsps.iter().cloned())
        .filter(|hsp| !(hsp.hit.e_value > THRESHOLD))
        .collect();
    let want = report::culled_hit_order(
        hits,
        500,
        1,
        &query_lengths,
        &HitFilters {
            max_hsps: 0,
            query_cov_hsp_perc: 0.0,
            query_lengths: &query_lengths,
            subject_besthit: false,
        },
        true,
    );
    assert!(!want.is_empty());
    let collector = HitCollector::with_options(500, 1, 0, 0.0, query_lengths.clone(), THRESHOLD);
    assert_eq!(
        key(&collect(&lists, &shuffled(subjects), collector)),
        key(&want)
    );
}

// A subject that never sends (no subject of the loop today) does not lose the ones past it:
// they enter at the end, in OID order.
#[test]
fn subjects_past_a_missing_one_enter_at_the_end() {
    let lists: Vec<(usize, Vec<TblastxHsp>)> = subject_lists(30)
        .into_iter()
        .map(|(s, hsps)| {
            if s == 4 || s == 17 {
                (s, Vec::new())
            } else {
                (s, hsps)
            }
        })
        .collect();
    let want = reference(&lists, 500, &HitFilters::default());
    let order: Vec<usize> = shuffled(30)
        .into_iter()
        .filter(|&s| s != 4 && s != 17)
        .collect();
    let collector = HitCollector::with_options(500, 0, 0, 0.0, Vec::new(), THRESHOLD);
    assert_eq!(key(&collect(&lists, &order, collector)), key(&want));
}
