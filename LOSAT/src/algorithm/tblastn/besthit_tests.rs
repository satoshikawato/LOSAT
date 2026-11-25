use super::*;
use crate::algorithm::tblastn::search_gapped::GappedHsp;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;

fn hsp(context: usize, frame: i8, qs: i32, qe: i32, score: i32, evalue: f64) -> LinkedHsp {
    LinkedHsp {
        context,
        hsp: GappedHsp {
            frame,
            score,
            q_start: qs,
            q_end: qe,
            q_gapped_start: qs,
            s_start: qs + 7,
            s_end: qe + 7,
            s_gapped_start: qs + 7,
        },
        num: 1,
        evalue,
        source_index: 0,
    }
}

fn payload(num_ident: i32) -> KappaHspPayload {
    KappaHspPayload {
        bit_score: 0.0,
        num_ident,
        num_positives: 0,
        report_num_ident: 0,
        report_mismatches: 0,
        align_length: 0,
        mismatches: 0,
        gap_opens: 0,
        gap_letters: 0,
        edit_script: Vec::new(),
        matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
        query_left_out: 0,
    }
}

fn summary(list: &KappaResultHitList) -> Vec<(i32, Vec<i32>)> {
    list.lists()
        .iter()
        .map(|l| (l.oid, l.hsps.hsps.iter().map(|h| h.hsp.score).collect()))
        .collect()
}

fn params() -> BestHitParams {
    BestHitParams {
        overhang: 0.1,
        score_edge: 0.1,
    }
}

/// A final hit list of `lists` (OID, HSPs with the payload's `num_ident` as a tag), saved as the
/// traceback saves them.
fn final_list(lists: Vec<(i32, Vec<(LinkedHsp, i32)>)>, max: usize) -> KappaResultHitList {
    let mut hit_list = KappaResultHitList::new(max).unwrap();
    for (oid, hsps) in lists {
        let (hsps, tags): (Vec<_>, Vec<_>) = hsps.into_iter().unzip();
        hit_list
            .update(KappaResultList {
                oid,
                hsps: LinkedHspList {
                    hsps,
                    best_evalue: 0.0,
                },
                payloads: tags.into_iter().map(payload).collect(),
            })
            .unwrap();
    }
    hit_list
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266: the key of a
// tblastn HSP (a protein query, frame 0) is its own 0-based query offset, its query extent, its
// score and e-value; the subject frame is not read.
#[test]
fn the_key_is_the_query_offset_and_extent() {
    assert_eq!(
        best_hit_key(&hsp(0, -2, 10, 60, 55, 1e-5)),
        BestHitKey {
            begin: 10,
            len: 50,
            score: 55,
            evalue: 1e-5
        }
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:271-288: a later
// subject's HSP inside a better HSP of an earlier subject is dropped, whatever their frames (the
// frames of a subject are not contexts); one of another query over the same range is kept, and a
// query without HSPs has no hit list.
#[test]
fn the_writer_drops_an_hsp_inside_a_better_one_of_its_query() {
    let mut writer = PrelimBestHit::new(3, 10, params());
    writer.add(0, &[hsp(0, 1, 100, 300, 400, 1e-60)]).unwrap();
    writer
        .add(
            1,
            &[
                hsp(0, -3, 150, 250, 100, 1e-10),
                hsp(0, 2, 500, 600, 90, 1e-9),
                hsp(2, -3, 150, 250, 100, 1e-10),
            ],
        )
        .unwrap();
    let lists = writer.finish().unwrap();
    assert_eq!(lists.len(), 3);
    assert_eq!(
        summary(lists[0].as_ref().unwrap()),
        vec![(0, vec![400]), (1, vec![90])]
    );
    assert!(lists[1].is_none());
    assert_eq!(summary(lists[2].as_ref().unwrap()), vec![(1, vec![100])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:207-221: the lists of
// a query are sorted by e-value and saved into a hit list of `prelim_hitlist_size`, which keeps
// the best ones (blast_hits.c:3243-3297: the full hit list sorts the HSPs of its lists by
// e-value).
#[test]
fn the_writer_hit_list_keeps_the_best_lists() {
    let mut writer = PrelimBestHit::new(1, 2, params());
    for oid in 0..4usize {
        let start = 1000 * oid as i32;
        writer
            .add(
                oid,
                &[
                    hsp(
                        0,
                        1,
                        start,
                        start + 50,
                        40 + oid as i32,
                        1e-5 / (oid + 1) as f64,
                    ),
                    hsp(0, -1, start + 500, start + 600, 60, 1e-3),
                ],
            )
            .unwrap();
    }
    let lists = writer.finish().unwrap();
    let mut kept = summary(lists[0].as_ref().unwrap());
    kept.sort();
    assert_eq!(kept, vec![(2, vec![42, 60]), (3, vec![43, 60])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539: the pipe sorts
// the lists by e-value first, so the HSP of a worse list saved first in the hit list is dropped
// by the better list's HSP that holds it with its overhang; the payloads follow their HSPs and
// the kept lists leave sorted by e-value.
#[test]
fn the_pipe_reads_the_lists_by_evalue() {
    let mut hit_list = final_list(
        vec![
            (7, vec![(hsp(0, 1, 100, 200, 80, 1e-12), 71)]),
            (
                4,
                vec![
                    (hsp(0, -2, 95, 205, 300, 1e-50), 41),
                    (hsp(0, 3, 400, 500, 70, 1e-11), 42),
                ],
            ),
            (9, vec![(hsp(0, 1, 800, 900, 50, 1e-3), 91)]),
        ],
        10,
    );
    hit_list.best_hit(params(), 10).unwrap();
    assert_eq!(summary(&hit_list), vec![(4, vec![300, 70]), (9, vec![50])]);
    let tags: Vec<Vec<i32>> = hit_list
        .lists()
        .iter()
        .map(|l| l.payloads.iter().map(|p| p.num_ident).collect())
        .collect();
    assert_eq!(tags, vec![vec![41, 42], vec![91]]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:88-186,520-539: the pipe
// saves the HSPs of a list in the list's own order (by score here), not by e-value. With a hit
// list size of 2 the fifth HSP rebuilds the query's list, which drops the overhang of the HSPs
// saved so far: the HSP of score 300 saved first then no longer holds the HSP of subject 2 that
// starts 5 letters before it, which is kept; saved after the rebuild (the e-value order of the
// list) it would still hold it with its overhang and drop it.
#[test]
fn the_pipe_keeps_the_order_of_the_hsps_of_a_list() {
    let far = |i: i32| {
        (
            hsp(0, 1, 1000 * i, 1000 * i + 100, 100, 10f64.powi(-40 - i)),
            i,
        )
    };
    let mut hit_list = final_list(
        vec![
            (
                1,
                vec![
                    (hsp(0, 1, 100, 200, 300, 1e-10), 0),
                    far(1),
                    far(2),
                    far(3),
                    far(4),
                    far(5),
                ],
            ),
            (2, vec![(hsp(0, 2, 95, 150, 55, 1e-5), 20)]),
        ],
        10,
    );
    hit_list.best_hit(params(), 2).unwrap();
    assert_eq!(
        summary(&hit_list),
        vec![(1, vec![300, 100, 100, 100, 100, 100]), (2, vec![55])]
    );
}
