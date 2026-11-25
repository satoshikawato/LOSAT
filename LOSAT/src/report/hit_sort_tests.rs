//! Tests of `hit_sort`: the replay of the calls that the C++ library's list sort was seen to make
//! (`LOSAT/tests/unit/list_sort_observed/README.md`), NCBI's comparators, the orders of
//! `sort_query_hits` (ties, the comparisons' side effects on the hits' HSP orders, HSPs without
//! identities) and the values that the description table's columns read.

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
        query_sort_start: query.0.min(query.1) as i64,
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
    // Equal scores: `>` is false both ways.
    assert!(!hsp_by_score_descending(&a, &b) && !hsp_by_score_descending(&b, &a));
    // Equal starts: the lower E-value first.
    assert!(hsp_by_query_start_ascending(&b, &a) && !hsp_by_query_start_ascending(&a, &b));
    assert!(hsp_by_subject_start_ascending(&b, &a) && !hsp_by_subject_start_ascending(&a, &b));
    let c = hsp(10.0, 1.0, 1, 5, (9, 9), 2);
    assert!(hsp_by_query_start_ascending(&c, &a) && hsp_by_subject_start_ascending(&c, &a));
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

/// The fields of a line of `list_sort_calls.txt`: the list's length, the answers, the calls (the
/// input positions of the comparator's first and second elements) and the order after the sort.
type ObservedSort = (usize, Vec<bool>, Vec<(usize, usize)>, Vec<usize>);

fn observed_sort(line: &str) -> ObservedSort {
    let fields: Vec<&str> = line.split(' ').collect();
    let (n, answers, calls, order) = match fields[0] {
        "T" => (fields[1], fields[2], fields[3], fields[4]),
        "R" => (fields[2], fields[4], fields[5], fields[6]),
        kind => panic!("a line of kind {kind}: {line}"),
    };
    let split = |field: &str| -> Vec<String> {
        if field == "-" {
            Vec::new()
        } else {
            field.split(',').map(str::to_string).collect()
        }
    };
    let answers = if answers == "-" {
        Vec::new()
    } else {
        answers.bytes().map(|answer| answer == b'1').collect()
    };
    let calls = split(calls)
        .iter()
        .map(|call| {
            let (first, second) = call.split_once(':').unwrap();
            (first.parse().unwrap(), second.parse().unwrap())
        })
        .collect();
    let order = split(order).iter().map(|k| k.parse().unwrap()).collect();
    (n.parse().unwrap(), answers, calls, order)
}

// For each recorded sequence of answers (every one for lists of 0 to 5 elements, and random ones
// for lists of 6 to 257 elements), `list_sort` hands the comparator the elements that g++ 13's
// `std::list::sort` handed it, in the same argument order and sequence, and leaves the same order.
#[test]
fn list_sort_calls_as_observed() {
    let records = include_str!("../../tests/unit/list_sort_observed/list_sort_calls.txt");
    let (mut lists, mut all_calls) = (0, 0);
    for line in records.lines() {
        let (n, answers, calls, order) = observed_sort(line);
        assert_eq!(answers.len(), calls.len(), "{line}");
        let mut list: Vec<usize> = (0..n).collect();
        let mut asked = Vec::new();
        list_sort(&mut list, |&first, &second| {
            asked.push((first, second));
            answers[asked.len() - 1]
        });
        assert_eq!(asked, calls, "{line}");
        assert_eq!(list, order, "{line}");
        lists += 1;
        all_calls += calls.len();
    }
    assert_eq!((lists, all_calls), (182, 7563));
}

// NCBI BLAST+ 2.17.0 (blastp -threshold 9 -comp_based_stats 0 -seg no, outfmt 0 -sorthsps 3)
// printed one subject's HSPs of 100 % (E-value 2e-14), 51 % (2e-13), 0 % (2e-06) and 80 % (3e-04)
// in that order, the order without -sorthsps. SortHspByPercentIdentityDescending compares the
// HSP without identities with every HSP by E-value, so it orders these four in a cycle (80 %
// before 51 % by ratio, 51 % before 0 % and 0 % before 80 % by E-value), and the order is the one
// that the calls of the list sort give.
#[test]
fn hsps_in_a_cycle_take_the_order_of_ncbi_s_calls() {
    let hsps = vec![
        hsp(55.5, 2e-14, 26, 26, (0, 25), 0),
        hsp(53.1, 2e-13, 20, 39, (30, 68), 30),
        hsp(32.3, 2e-6, 0, 24, (70, 93), 70),
        hsp(25.8, 3e-4, 8, 10, (100, 109), 100),
    ];
    let hits = vec![hsps];
    assert_eq!(
        sorted(by(None, Some(3)), false, &hits).hsps,
        vec![vec![0, 1, 2, 3]]
    );
    assert_eq!(
        sorted(by(Some(3), Some(3)), false, &hits).hsps,
        vec![vec![0, 1, 2, 3]]
    );
}

// NCBI BLAST+ 2.17.0 (the same options, -sorthits 3) printed the subjects of 100 % (E-value
// 3e-16), 0 % (6e-08) and 80 % (1e-05) in that order: a hit of 0 % ties with every hit
// (SortHitByPercentIdentityDescendingEx compares its E-values of -1), while 80 % goes after 100 %,
// so the order is the one that the calls of the list sort give.
#[test]
fn a_hit_of_zero_percent_takes_the_order_of_ncbi_s_calls() {
    let hits = vec![
        vec![hsp(55.5, 3e-16, 26, 26, (0, 25), 0)],
        vec![hsp(32.3, 6e-8, 0, 24, (70, 93), 0)],
        vec![hsp(25.8, 1e-5, 8, 10, (100, 109), 0)],
    ];
    assert_eq!(sorted(by(Some(3), None), false, &hits).hits, vec![0, 1, 2]);
}

// NCBI BLAST+ 2.17.0 (the same options) printed one subject's HSPs of 51 % (E-value 5e-14), 0 %
// (1e-06) and 80 % (2e-04), in that E-value order, as 80 %, 51 %, 0 % with -sorthsps 3, and three
// subjects of one such HSP each (E-values 3e-15, 8e-08, 2e-05) as 80 %, 51 %, 0 % with -sorthits
// 3. Both lists are cycles of the comparator that the list sort's calls move: a stable sort by the
// comparator, or a merge in another order, leaves another order.
#[test]
fn cycles_move_as_ncbi_s_calls_move_them() {
    let hsps = vec![
        hsp(53.1, 5e-14, 20, 39, (30, 68), 0),
        hsp(32.3, 1e-6, 0, 24, (70, 93), 77),
        hsp(25.8, 2e-4, 8, 10, (100, 109), 132),
    ];
    assert_eq!(
        sorted(by(None, Some(3)), false, &[hsps]).hsps,
        vec![vec![2, 0, 1]]
    );
    let hits = vec![
        vec![hsp(53.1, 3e-15, 20, 39, (30, 68), 0)],
        vec![hsp(32.3, 8e-8, 0, 24, (70, 93), 0)],
        vec![hsp(25.8, 2e-5, 8, 10, (100, 109), 0)],
    ];
    assert_eq!(sorted(by(Some(3), None), false, &hits).hits, vec![2, 0, 1]);
}

// Two more cycles of SortHspByPercentIdentityDescending (HSPs in E-value order, one without
// identities), in the orders that g++ 13's std::list::sort leaves them.
#[test]
fn hsp_cycles_take_the_order_of_the_list_sort_s_calls() {
    let first = vec![
        hsp(100.0, 1e-20, 7, 30, (0, 29), 0),
        hsp(90.0, 1e-14, 0, 30, (40, 69), 40),
        hsp(80.0, 1e-10, 4, 10, (80, 89), 80),
    ];
    let second = vec![
        hsp(80.0, 1e-5, 0, 20, (0, 19), 0),
        hsp(90.0, 1e-20, 7, 10, (30, 39), 30),
        hsp(70.0, 1e-3, 29, 30, (50, 79), 50),
    ];
    assert_eq!(
        sorted(by(None, Some(3)), false, &[first, second]).hsps,
        vec![vec![2, 0, 1], vec![2, 1, 0]]
    );
}
