use super::*;
use crate::algorithm::tblastn::search_gapped::GappedHsp;

fn hsp(context: usize, frame: i8, qs: i32, qe: i32, ss: i32, score: i32) -> LinkedHsp {
    LinkedHsp {
        context,
        hsp: GappedHsp {
            frame,
            score,
            q_start: qs,
            q_end: qe,
            q_gapped_start: qs,
            s_start: ss,
            s_end: ss + (qe - qs),
            s_gapped_start: ss,
        },
        num: 1,
        evalue: 1e-5 / f64::from(score),
        source_index: 0,
    }
}

fn summary(list: &KappaResultHitList) -> Vec<(i32, Vec<i32>)> {
    list.lists()
        .iter()
        .map(|l| (l.oid, l.hsps.hsps.iter().map(|h| h.hsp.score).collect()))
        .collect()
}

#[test]
fn the_key_is_the_query_range_and_the_frame_local_subject_offset() {
    let key = cull_key(&hsp(0, -2, 10, 50, 7, 33));
    assert_eq!(
        key,
        CullKey {
            begin: 10,
            end: 50,
            score: 33,
            subject_offset: 7
        }
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120, 612-633: the frame is no key,
// so an HSP on another frame over the same query range is dominated by a better one,
// and the tree is the query's (`cid` = context), with the length of that query.
#[test]
fn the_frames_of_a_subject_share_the_tree_of_the_query() {
    let mut culling = PrelimCulling::new(&[1000, 500], 10, 1);
    culling
        .add(
            0,
            &[
                hsp(0, 1, 100, 200, 50, 90),
                hsp(0, -1, 100, 200, 60, 40),
                hsp(0, 2, 400, 500, 10, 30),
                // Another query: its own tree.
                hsp(1, -3, 100, 200, 50, 20),
            ],
        )
        .unwrap();
    let lists = culling.finish().unwrap();
    assert_eq!(summary(lists[0].as_ref().unwrap()), vec![(0, vec![90, 30])]);
    assert_eq!(summary(lists[1].as_ref().unwrap()), vec![(0, vec![20])]);
}

// Two frames tie on all five keys when their frame-local offsets are equal; the one
// written later loses (`s_DominateTest`, the later insert is the dominated one).
#[test]
fn frames_that_tie_on_the_five_keys_keep_the_first_written() {
    let mut culling = PrelimCulling::new(&[1000], 10, 1);
    let mut first = hsp(0, 1, 100, 200, 50, 90);
    first.num = 11;
    let mut second = hsp(0, -1, 100, 200, 50, 90);
    second.num = 22;
    culling.add(0, &[first, second]).unwrap();
    let lists = culling.finish().unwrap();
    let kept = lists[0].as_ref().unwrap();
    assert_eq!(kept.lists().len(), 1);
    assert_eq!(kept.lists()[0].hsps.hsps.len(), 1);
    assert_eq!(kept.lists()[0].hsps.hsps[0].num, 11);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-581: the hit list is made with
// `prelim_hitlist_size` and holds more lists than that, in the rip order, with the best
// e-value, worst e-value and low score set, and is not heapified.
#[test]
fn the_hit_list_is_not_cut_and_keeps_the_rip_order() {
    let mut culling = PrelimCulling::new(&[1000], 2, 5);
    for oid in 0..4usize {
        culling
            .add(
                oid,
                &[
                    hsp(
                        0,
                        1,
                        100 + 50 * oid as i32,
                        120 + 50 * oid as i32,
                        1,
                        40 - oid as i32,
                    ),
                    hsp(0, 1, 100 + 50 * oid as i32, 110 + 50 * oid as i32, 2, 20),
                ],
            )
            .unwrap();
    }
    let lists = culling.finish().unwrap();
    let q0 = lists[0].as_ref().unwrap();
    assert_eq!(q0.max, 2);
    assert_eq!(q0.lists().len(), 4);
    assert!(!q0.heapified());
    assert_eq!(q0.low_score(), 37);
    // The best e-value of subject 3's list: the score 37 HSP.
    assert!((q0.worst_evalue() - 1e-5 / 37.0).abs() < 1e-18);
    // Limit 5 (+3 = 8) keeps everything: the newest subject first, scores descending.
    assert_eq!(
        summary(q0),
        vec![
            (3, vec![37, 20]),
            (2, vec![38, 20]),
            (1, vec![39, 20]),
            (0, vec![40, 20])
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:100-112: the lower OID wins a tie, and
// the subjects must come in OID order for the later one to lose.
#[test]
fn identical_hsps_of_two_subjects_keep_the_lower_oid() {
    let mut culling = PrelimCulling::new(&[1000], 10, 1);
    culling.add(0, &[hsp(0, 1, 100, 200, 5, 80)]).unwrap();
    culling.add(1, &[hsp(0, 1, 100, 200, 5, 80)]).unwrap();
    let lists = culling.finish().unwrap();
    assert_eq!(summary(lists[0].as_ref().unwrap()), vec![(0, vec![80])]);
}

fn payload(bit_score: f64) -> KappaHspPayload {
    use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
    KappaHspPayload {
        bit_score,
        num_ident: 0,
        num_positives: 0,
        report_num_ident: 0,
        report_mismatches: 0,
        align_length: 0,
        mismatches: 0,
        gap_opens: 0,
        gap_letters: 0,
        edit_script: Vec::new(),
        matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
    }
}

// The pipe: lists sorted by e-value (payloads stay with their HSPs), the culled HSPs
// dropped with their payloads, the survivors score sorted in the new lists.
#[test]
fn the_pipe_culls_in_evalue_order_and_keeps_the_payloads_paired() {
    let mut hit_list = KappaResultHitList::new(500).unwrap();
    // Subject 0: a weak HSP (e-value order puts it after the strong one of subject 1).
    let weak = hsp(0, 1, 100, 200, 5, 40);
    let strong = hsp(0, 1, 100, 200, 5, 90);
    let other = hsp(0, 1, 400, 480, 5, 30);
    let list0 = KappaResultList {
        oid: 0,
        hsps: LinkedHspList {
            hsps: vec![weak, other],
            best_evalue: weak.evalue,
        },
        payloads: vec![payload(40.0), payload(30.0)],
    };
    let list1 = KappaResultList {
        oid: 1,
        hsps: LinkedHspList {
            hsps: vec![strong],
            best_evalue: strong.evalue,
        },
        payloads: vec![payload(90.0)],
    };
    hit_list.update(list0).unwrap();
    hit_list.update(list1).unwrap();
    hit_list.cull(1, 1000).unwrap();
    // The weak HSP is dominated by the strong one (overlap 100 %, better score). The
    // lists leave in the rip order of the tree: the node of the 400-480 HSP holds it
    // before the left subtree that has the strong HSP.
    assert_eq!(summary(&hit_list), vec![(0, vec![30]), (1, vec![90])]);
    let bits: Vec<Vec<f64>> = hit_list
        .lists()
        .iter()
        .map(|l| l.payloads.iter().map(|p| p.bit_score).collect())
        .collect();
    assert_eq!(bits, vec![vec![30.0], vec![90.0]]);
    assert_eq!(hit_list.max, 500);
    assert_eq!(hit_list.low_score(), 30);
    assert!(!hit_list.heapified());
}
