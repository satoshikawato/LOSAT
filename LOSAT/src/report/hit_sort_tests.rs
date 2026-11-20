//! Tests of `hit_sort`: NCBI's comparators, the orders of `sort_query_hits` (ties, the
//! comparisons' side effects on the hits' HSP orders, HSPs without identities) and the values that
//! the description table's columns read.

use super::*;

fn hsp(
    bits: f64,
    evalue: f64,
    num_ident: i64,
    length: i64,
    query: (u64, u64),
    subject: u64,
) -> HspScores {
    HspScores {
        bits,
        evalue,
        num_ident,
        length,
        query_from: query.0,
        query_to: query.1,
        subject_from: subject,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:828-830 (quoted
// at `PairwiseHit::seqalign_num_ident`): `GetAlnScores` reads the Seq-align's `num_ident`, the
// engine's identities (here 270, as TBLASTN's engine counts SEG-masked letters as mismatches),
// not the 284 identities of the displayed rows; the coordinates are 0-based, from <= to.
#[test]
fn hsp_scores_read_the_seqalign_identities() {
    let mut hit = PairwiseHit::from(crate::common::Hit {
        identity: 99.3,
        length: 286,
        mismatch: 2,
        gapopen: 0,
        q_start: 1,
        q_end: 286,
        s_start: 1858,
        s_end: 1001,
        e_value: 0.0,
        bit_score: 557.0,
        num_ident: 284,
        query_frame: 0,
        query_length: 286,
        q_idx: 0,
        s_idx: 0,
        raw_score: 1435,
        sort_query_offset: 0,
        sort_query_end: 286,
        sort_subject_offset: 0,
        sort_subject_end: 286,
        has_sort_offsets: true,
        gap_info: None,
        num_positives: 285,
    });
    // `From<Hit>` takes the hit's identities (BLASTN's rows are not counted again).
    assert_eq!(HspScores::of(&hit).num_ident, 284);
    hit.seqalign_num_ident = 270;
    assert_eq!(
        HspScores::of(&hit),
        hsp(557.0, 0.0, 270, 286, (0, 285), 1000)
    );
}

#[test]
fn hsp_comparators_break_ties_by_e_value() {
    let a = hsp(50.0, 1e-5, 40, 50, (10, 59), 3);
    let b = hsp(50.0, 1e-6, 40, 50, (10, 59), 3);
    // Equal scores: `>` is false both ways, equal.
    assert_eq!(hsp_by_score_descending(&a, &b), Ordering::Equal);
    let higher = hsp(51.0, 1e-5, 40, 50, (10, 59), 3);
    assert_eq!(hsp_by_score_descending(&higher, &a), Ordering::Less);
    assert_eq!(hsp_by_score_descending(&a, &higher), Ordering::Greater);
    // Equal starts: the lower E-value first.
    assert_eq!(hsp_by_query_start_ascending(&b, &a), Ordering::Less);
    assert_eq!(hsp_by_query_start_ascending(&a, &b), Ordering::Greater);
    assert_eq!(hsp_by_subject_start_ascending(&b, &a), Ordering::Less);
    assert_eq!(hsp_by_subject_start_ascending(&a, &b), Ordering::Greater);
    assert_eq!(hsp_by_query_start_ascending(&a, &a), Ordering::Equal);
    let c = hsp(10.0, 1.0, 1, 5, (9, 9), 2);
    assert_eq!(hsp_by_query_start_ascending(&c, &a), Ordering::Less);
    assert_eq!(hsp_by_subject_start_ascending(&c, &a), Ordering::Less);
    // -0.0 and 0.0 are equal, as C++ compares them.
    let zero = hsp(0.0, 0.0, 1, 5, (9, 9), 2);
    let minus_zero = hsp(-0.0, -0.0, 1, 5, (9, 9), 2);
    assert_eq!(hsp_by_score_descending(&zero, &minus_zero), Ordering::Equal);
    assert_eq!(
        hsp_by_query_start_ascending(&zero, &minus_zero),
        Ordering::Equal
    );
    // Equal identity ratios (40/50 and 4/5): the E-values.
    let d = hsp(20.0, 1e-7, 4, 5, (0, 4), 0);
    assert!(hsp_by_percent_identity_descending(&d, &a, false, false));
    assert!(!hsp_by_percent_identity_descending(&a, &d, false, false));
    // Equal ratios and E-values: false both ways.
    assert!(!hsp_by_percent_identity_descending(&a, &a, false, false));
    // Different ratios: the higher first (`>=` on unequal ratios).
    let e = hsp(20.0, 1e-9, 3, 5, (0, 4), 0);
    assert!(hsp_by_percent_identity_descending(&a, &e, false, false));
    assert!(!hsp_by_percent_identity_descending(&e, &a, false, false));
    // No identities on one side: the E-values decide, whatever the ratios.
    let none = hsp(20.0, 1e-3, 0, 5, (0, 4), 0);
    assert!(!hsp_by_percent_identity_descending(&none, &a, false, false));
    assert!(hsp_by_percent_identity_descending(&a, &none, false, false));
    let early = hsp(20.0, 1e-30, 0, 5, (0, 4), 0);
    assert!(hsp_by_percent_identity_descending(&early, &a, false, false));
}

#[test]
fn tblastx_lengths_are_nucleotides_unless_translated() {
    let a = hsp(50.0, 1e-5, 10, 30, (0, 89), 0);
    assert_eq!(a.alignment_length(false, false), 30);
    assert_eq!(a.alignment_length(false, true), 30);
    assert_eq!(a.alignment_length(true, false), 90);
    assert_eq!(a.alignment_length(true, true), 30);
    // The identity of the description table (no translation) is a third of the shown one.
    assert_eq!(
        percent_identity(30, a.alignment_length(true, false)) as i64,
        33
    );
    assert_eq!(percent_identity(30, a.alignment_length(true, true)), 100.0);
    assert_eq!(percent_identity(7, 7), 100.0);
    assert_eq!(percent_identity(0, 1), 0.0);
}

#[test]
fn query_coverage_is_the_union_of_the_query_ranges() {
    let ranges = [(10, 19), (15, 30), (40, 40), (31, 35), (0, 2)];
    let hsps: Vec<HspScores> = ranges.iter().map(|&q| hsp(1.0, 1.0, 1, 1, q, 0)).collect();
    // 0-2 (3), 10-35 (26, adjacent ranges count once each), 40 (1).
    assert_eq!(query_coverage(&hsps), 30);
    assert_eq!(query_coverage(&hsps[..1]), 10);
    assert_eq!(query_coverage(&[]), 0);
    let contained = [
        hsp(1.0, 1.0, 1, 1, (0, 99), 0),
        hsp(1.0, 1.0, 1, 1, (10, 20), 0),
    ];
    assert_eq!(query_coverage(&contained), 100);
}

fn sorted(sort: ReportSort, tblastx: bool, hits: &[Vec<HspScores>]) -> SortedHits {
    sort_query_hits(sort, tblastx, hits).expect("a sort")
}

fn by(hits: Option<u8>, hsps: Option<u8>) -> ReportSort {
    ReportSort { hits, hsps }
}

#[test]
fn nothing_is_sorted_without_a_positive_option() {
    let hits = vec![vec![hsp(1.0, 1.0, 1, 1, (0, 0), 0)]; 3];
    for sort in [
        by(None, None),
        by(Some(0), None),
        by(None, Some(0)),
        by(Some(0), Some(0)),
    ] {
        assert!(sort_query_hits(sort, false, &hits).is_none(), "{sort:?}");
        assert!(!sort.sorts());
    }
    assert!(by(Some(0), None).shows_columns() && !by(None, Some(3)).shows_columns());
    assert!(by(Some(0), Some(1)).sorts() && by(Some(4), None).sorts());
}

#[test]
fn equal_total_scores_come_in_reverse() {
    // Three copies of one subject (dup9, dup1, msA in the search's order) and a lower one.
    let copy = vec![
        hsp(30.0, 1e-3, 20, 25, (0, 24), 0),
        hsp(10.0, 1.0, 5, 8, (30, 37), 0),
    ];
    let low = vec![hsp(35.0, 1e-4, 20, 25, (0, 24), 0)];
    let hits = vec![copy.clone(), copy.clone(), low, copy];
    let order = sorted(by(Some(2), None), false, &hits);
    // Totals 40, 40, 35, 40: `>=` puts each later 40 before the earlier ones.
    assert_eq!(order.hits, vec![3, 1, 0, 2]);
    // The HSPs keep their order.
    assert_eq!(
        order.hsps,
        vec![vec![0, 1], vec![0, 1], vec![0], vec![0, 1]]
    );
}

#[test]
fn the_score_sort_orders_the_hsps_of_compared_hits_only() {
    // The second HSP has the higher score.
    let two = vec![
        hsp(20.0, 1e-3, 10, 12, (0, 11), 0),
        hsp(25.0, 1e-2, 10, 12, (20, 31), 5),
    ];
    // One hit: no comparison, the HSPs stay in their order.
    let order = sorted(by(Some(1), None), false, std::slice::from_ref(&two));
    assert_eq!(
        order,
        SortedHits {
            hits: vec![0],
            hsps: vec![vec![0, 1]]
        }
    );
    // Two hits: the comparison sorts both hits' HSPs by score.
    let other = vec![hsp(22.0, 1e-4, 10, 12, (0, 11), 0)];
    let order = sorted(by(Some(1), None), false, &[other.clone(), two.clone()]);
    assert_eq!(order.hits, vec![1, 0]);
    assert_eq!(order.hsps, vec![vec![0], vec![1, 0]]);
    // Then -sorthsps 2 puts them back by query start.
    let order = sorted(by(Some(1), Some(2)), false, &[other, two]);
    assert_eq!(order.hsps, vec![vec![0], vec![0, 1]]);
}

#[test]
fn the_percent_identity_sort_reads_the_best_scoring_hsp() {
    // Hit 0: best HSP 18/20 (90 %); hit 1: best HSP 30/40 (75 %) despite an HSP of 100 %.
    let hit0 = vec![hsp(30.0, 1e-5, 18, 20, (0, 19), 0)];
    let hit1 = vec![
        hsp(40.0, 1e-8, 30, 40, (10, 49), 0),
        hsp(5.0, 1.0, 4, 4, (0, 3), 0),
    ];
    let order = sorted(by(Some(3), None), false, &[hit1.clone(), hit0.clone()]);
    assert_eq!(order.hits, vec![1, 0]);
    // The comparison sorted hit 1's HSPs by percent identity (100 % first).
    assert_eq!(order.hsps[0], vec![1, 0]);
    // Equal percentages compare E-values of -1: false, so the order stays.
    let same = vec![hsp(30.0, 1e-9, 18, 20, (0, 19), 0)];
    let order = sorted(by(Some(3), None), false, &[hit0.clone(), same.clone()]);
    assert_eq!(order.hits, vec![0, 1]);
    let order = sorted(by(Some(3), None), false, &[same, hit0]);
    assert_eq!(order.hits, vec![0, 1]);
    // No identities: a percentage of 0 is not above 0, the order stays.
    let zero = vec![hsp(50.0, 1e-20, 0, 20, (0, 19), 0)];
    let best = vec![hsp(10.0, 1.0, 20, 20, (0, 19), 0)];
    let order = sorted(by(Some(3), None), false, &[zero.clone(), best.clone()]);
    assert_eq!(order.hits, vec![0, 1]);
}

#[test]
fn tblastx_s_identity_sort_takes_amino_acid_lengths() {
    // 9/10 against 25/30 amino acids: by amino acids 90 % against 83 %; in nucleotides (no
    // translation) 30 % against 27.8 %, the same order. The HSP sort with -sorthits 3 reads
    // the same lengths (kTranslation), without it nucleotides: here an HSP of 9 identities in 10
    // letters and one of 28 in 30 order the same either way, but the hit sort's percentages
    // differ by the factor three (GetPercentIdentity of 9 and 30 is 30 %).
    let a = vec![hsp(30.0, 1e-5, 9, 10, (0, 29), 0)];
    let b = vec![hsp(40.0, 1e-7, 25, 30, (0, 89), 0)];
    assert_eq!(
        sorted(by(Some(3), None), true, &[b.clone(), a.clone()]).hits,
        vec![1, 0]
    );
    assert_eq!(set_percent_identity(&a, &[0], true, true), 90.0);
    assert_eq!(set_percent_identity(&a, &[0], true, false), 30.0);
    // An HSP whose identities equal its amino-acid length is 100 % only with the translation.
    let full = vec![hsp(30.0, 1e-5, 10, 10, (0, 29), 0)];
    assert_eq!(set_percent_identity(&full, &[0], true, true), 100.0);
    assert!((set_percent_identity(&full, &[0], true, false) - 100.0 / 3.0).abs() < 1e-12);
}

#[test]
fn the_coverage_sort_breaks_ties_by_the_first_hsp_s_e_value() {
    let wide = vec![hsp(20.0, 1e-2, 10, 50, (0, 49), 0)];
    let narrow = vec![hsp(60.0, 1e-20, 10, 20, (0, 19), 0)];
    let tie_low = vec![hsp(20.0, 1e-3, 10, 50, (100, 149), 0)];
    let order = sorted(by(Some(4), None), false, &[narrow, wide, tie_low]);
    // 50, 50 (the lower E-value first), then 20.
    assert_eq!(order.hits, vec![2, 1, 0]);
}

/// A BLASTP HSP of the query's letters 1-400 with E-value `e_value`.
fn pairwise_hit(e_value: f64) -> PairwiseHit {
    PairwiseHit::from(crate::common::Hit {
        identity: 100.0,
        length: 400,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 400,
        s_start: 1,
        s_end: 400,
        e_value,
        bit_score: 800.0,
        num_ident: 400,
        query_frame: 0,
        query_length: 400,
        q_idx: 0,
        s_idx: 0,
        raw_score: 2000,
        sort_query_offset: 0,
        sort_query_end: 400,
        sort_subject_offset: 0,
        sort_subject_end: 400,
        has_sort_offsets: true,
        gap_info: None,
        num_positives: 400,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1186
// ```c++
//     double evalue = (hsp->evalue < SMALLEST_EVALUE) ? 0.0 : hsp->evalue;
// ```
// The sorts compare the Seq-align's E-values: below 1e-180 both are 0.0, so `0.0 < 0.0` is false
// and the order stays (NCBI 2.17.0 `blastp -comp_based_stats 0 -sorthits 4` on a query and two
// subjects, the query itself and the query between 800 letters on each side, keeps the second
// subject, of the larger E-value, first: case d1.blastp.h4 of line RG T1).
#[test]
fn the_sorts_compare_the_seqalign_s_e_values() {
    assert_eq!(HspScores::of(&pairwise_hit(1e-181)).evalue, 0.0);
    assert_eq!(HspScores::of(&pairwise_hit(1e-180)).evalue, 1e-180);
    assert_eq!(HspScores::of(&pairwise_hit(2e-3)).evalue, 2e-3);
    // Equal coverage and scores, the larger E-value first (the search's order).
    let worse = HspScores::of(&pairwise_hit(1e-200));
    let better = HspScores::of(&pairwise_hit(1e-250));
    let hits = vec![vec![worse], vec![better]];
    assert_eq!(sorted(by(Some(4), None), false, &hits).hits, vec![0, 1]);
    // An E-value at 1e-180 or above still decides.
    let kept = HspScores::of(&pairwise_hit(1e-180));
    let hits = vec![vec![kept], vec![better]];
    assert_eq!(sorted(by(Some(4), None), false, &hits).hits, vec![1, 0]);
    // The HSP sorts' ties on equal starts and identities alike.
    let hsps = vec![vec![worse, better]];
    for hsp_sort in [2, 3, 4] {
        assert_eq!(
            sorted(by(None, Some(hsp_sort)), false, &hsps).hsps,
            vec![vec![0, 1]],
            "-sorthsps {hsp_sort}"
        );
    }
}

#[test]
fn hsp_sorts_apply_to_every_hit() {
    let hit = vec![
        hsp(10.0, 1e-1, 5, 10, (50, 59), 7),
        hsp(30.0, 1e-5, 9, 10, (5, 14), 90),
        hsp(20.0, 1e-3, 8, 10, (5, 14), 3),
    ];
    let hits = vec![hit.clone(), hit];
    assert_eq!(
        sorted(by(None, Some(1)), false, &hits).hsps,
        vec![vec![1, 2, 0]; 2]
    );
    // Equal query starts: the lower E-value.
    assert_eq!(
        sorted(by(None, Some(2)), false, &hits).hsps,
        vec![vec![1, 2, 0]; 2]
    );
    assert_eq!(
        sorted(by(None, Some(3)), false, &hits).hsps,
        vec![vec![1, 2, 0]; 2]
    );
    assert_eq!(
        sorted(by(Some(0), Some(4)), false, &hits).hsps,
        vec![vec![2, 0, 1]; 2]
    );
    // The hits keep their order.
    assert_eq!(sorted(by(None, Some(4)), false, &hits).hits, vec![0, 1]);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1282-1293
// ```c++
//     if(length1 > 0 && length2 > 0 && num_ident1 > 0 &&num_ident2 > 0 ) {
//         if (((double)num_ident1)/length1 == ((double)num_ident2)/length2) {
//
//             retval = evalue1 < evalue2;
//
//         } else {
//             retval = ((double)num_ident1)/length1 >= ((double)num_ident2)/length2;
//
//         }
//     } else {
//         retval = evalue1 < evalue2;
//     }
// ```
// An HSP without identities compares by E-value with every HSP. Where that orders the list, the
// order is the comparator's; where it orders three HSPs in a cycle, the HSPs with identities come
// first.
#[test]
fn hsps_without_identities_sort_by_e_value_where_the_comparator_orders_them() {
    let sort3 = |hsps: Vec<HspScores>| sorted(by(None, Some(3)), false, &[hsps]).hsps[0].clone();
    // 90 %, E 1e-3; no identities, E 1e-5: the E-values decide.
    let high = hsp(30.0, 1e-3, 9, 10, (0, 9), 0);
    let none = hsp(35.0, 1e-5, 0, 10, (20, 29), 20);
    assert_eq!(sort3(vec![high, none]), vec![1, 0]);
    assert_eq!(sort3(vec![none, high]), vec![0, 1]);
    // Two without identities and one with: by E-value (an order).
    let none2 = hsp(20.0, 1e-1, 0, 10, (40, 49), 40);
    assert_eq!(sort3(vec![none2, high, none]), vec![2, 1, 0]);
    // 50 % with E 1e-10 and 90 % with E 1e-3, the one without identities at E 1e-20: the
    // comparator orders the three (E 1e-20 first, then 90 % before 50 %).
    let low = hsp(60.0, 1e-10, 25, 50, (60, 109), 60);
    let best = hsp(70.0, 1e-20, 0, 10, (120, 129), 120);
    assert_eq!(sort3(vec![low, high, best]), vec![2, 1, 0]);
    // The one without identities at E 1e-5, between the two: a cycle (1e-10 before 1e-5 before
    // 1e-3 by E-value, 90 % before 50 %). The HSPs with identities first, by ratio, then the other.
    assert_eq!(sort3(vec![low, none, high]), vec![2, 0, 1]);
    // At an end of the interval (E 1e-3, equal to the 90 % HSP's): no order either.
    let tie = hsp(25.0, 1e-3, 0, 10, (140, 149), 140);
    assert_eq!(sort3(vec![tie, low, high]), vec![2, 1, 0]);
}

#[test]
fn percent_identity_order_holds_only_without_a_cycle() {
    let hsps = vec![
        hsp(30.0, 1e-3, 9, 10, (0, 9), 0),
        hsp(60.0, 1e-10, 25, 50, (0, 49), 0),
        hsp(35.0, 1e-5, 0, 10, (0, 9), 0),
        hsp(70.0, 1e-20, 0, 10, (0, 9), 0),
        hsp(71.0, f64::NAN, 5, 10, (0, 9), 0),
    ];
    let identities = |k: usize| {
        let h = &hsps[k];
        (h.num_ident > 0).then(|| h.num_ident as f64 / h.length as f64)
    };
    let holds = |order: &[usize]| percent_identity_order_holds(&hsps, order, identities);
    assert!(holds(&[0, 1]) && holds(&[0, 2]) && holds(&[2, 3]) && holds(&[0, 1, 3]));
    assert!(!holds(&[0, 1, 2]));
    // A NaN E-value: not an order.
    assert!(!holds(&[0, 4]));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1584-1593
// ```c++
//     if(percentIdent1 > 0 &&percentIdent2 > 0) {
//         if (percentIdent1 == percentIdent2) {
//             retval = evalue1 < evalue2;
//
//         } else {
//             retval = percentIdent1 >= percentIdent2;
//         }
//     } else {
//         retval = evalue1 < evalue2;
//     }
// ```
// A hit whose highest-scoring HSP has no identities (0 %) ties with every hit.
#[test]
fn a_hit_of_zero_percent_ties_with_every_hit() {
    let zero = vec![hsp(50.0, 1e-20, 0, 20, (0, 19), 0)];
    let half = vec![hsp(30.0, 1e-5, 10, 20, (0, 19), 0)];
    let most = vec![hsp(40.0, 1e-7, 18, 20, (0, 19), 0)];
    let order = |hits: &[Vec<HspScores>]| sorted(by(Some(3), None), false, hits).hits;
    // The comparator orders these lists: the order stays.
    assert_eq!(order(&[zero.clone(), half.clone()]), vec![0, 1]);
    assert_eq!(
        order(&[half.clone(), zero.clone(), half.clone()]),
        vec![0, 1, 2]
    );
    assert_eq!(order(&[half.clone(), most.clone()]), vec![1, 0]);
    // 50 % and 90 % with a 0 % hit between: not an order (90 % before 50 %, both tie with 0 %);
    // the hits above 0 % first, by percentage, then the others.
    assert_eq!(
        order(&[half.clone(), zero.clone(), most.clone()]),
        vec![2, 0, 1]
    );
    assert_eq!(order(&[zero.clone(), most, half]), vec![1, 2, 0]);
}

#[test]
fn nan_values_sort_without_a_panic() {
    let nan = hsp(f64::NAN, f64::NAN, 1, 2, (0, 1), 0);
    let one = hsp(10.0, 1e-3, 1, 2, (0, 1), 0);
    let two = hsp(20.0, 1e-4, 2, 2, (0, 1), 0);
    let hits = vec![vec![nan, one, two], vec![one], vec![nan], vec![two, nan]];
    for hit_sort in 0..=4 {
        for hsp_sort in 0..=4 {
            for tblastx in [false, true] {
                let sort = by(Some(hit_sort), Some(hsp_sort));
                if let Some(order) = sort_query_hits(sort, tblastx, &hits) {
                    let mut seen = order.hits.clone();
                    seen.sort_unstable();
                    assert_eq!(seen, vec![0, 1, 2, 3], "{sort:?}");
                }
            }
        }
    }
    // The bit scores: a NaN counts as above every number.
    assert_eq!(
        sorted(by(None, Some(1)), false, &hits).hsps[0],
        vec![0, 2, 1]
    );
}
