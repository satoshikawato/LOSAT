use super::*;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;

fn prelim(query: u32, qs: i32, qe: i32, score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: qs,
        query_end: qe,
        subject_start: 5,
        subject_end: 5 + (qe - qs),
        gapped_query_start: qs,
        gapped_subject_start: 5,
        raw_score: score,
        query_context: query as i32,
        query_frame: 0,
        subject_frame: 0,
        query_length: 1000,
        q_idx: query,
        s_idx: 0,
    }
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

fn finished(writer: PrelimBestHit, queries: usize) -> Vec<Option<BlastpHitList<XPrelimHsp>>> {
    let mut compact: Vec<Option<BlastpHitList<XPrelimHsp>>> = (0..queries).map(|_| None).collect();
    writer.finish(&mut compact);
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

fn params() -> BestHitParams {
    BestHitParams {
        overhang: 0.1,
        score_edge: 0.1,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:256-266: the key of a
// blastp HSP (frame 0) is its own 0-based query offset, its query extent, its score and e-value.
#[test]
fn the_key_is_the_query_offset_and_extent() {
    let hsp = XPrelimHsp {
        hsp: prelim(0, 10, 60, 55),
        e_value: 1e-5,
    };
    assert_eq!(
        best_hit_key(&hsp),
        BestHitKey {
            begin: 10,
            len: 50,
            score: 55,
            evalue: 1e-5
        }
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:271-288: a later
// subject's HSP inside a better HSP of an earlier subject is dropped; one of the same query
// elsewhere, and one of another query over the same range, are kept.
#[test]
fn the_writer_drops_an_hsp_inside_a_better_one_of_its_query() {
    let mut writer = PrelimBestHit::new(10, params());
    writer.add(&mut subject(0, 0, &[(prelim(0, 100, 300, 400), 1e-60)]));
    writer.add(&mut subject(
        1,
        0,
        &[
            (prelim(0, 150, 250, 100), 1e-10),
            (prelim(0, 500, 600, 90), 1e-9),
        ],
    ));
    writer.add(&mut subject(2, 1, &[(prelim(1, 150, 250, 100), 1e-10)]));
    let lists = finished(writer, 2);
    assert_eq!(scores(&lists, 0), vec![(0, vec![400]), (1, vec![90])]);
    assert_eq!(scores(&lists, 1), vec![(2, vec![100])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:207-221: the lists
// are sorted by e-value and saved into a hit list of `prelim_hitlist_size`, which keeps the
// best ones.
#[test]
fn the_hit_list_keeps_the_best_lists() {
    let mut writer = PrelimBestHit::new(2, params());
    for oid in 0..4usize {
        let start = 100 * oid as i32;
        writer.add(&mut subject(
            oid,
            0,
            &[(
                prelim(0, start, start + 50, 40 + oid as i32),
                1e-5 / (oid + 1) as f64,
            )],
        ));
    }
    let lists = finished(writer, 1);
    assert_eq!(lists[0].as_ref().unwrap().hsplist_count, 2);
    let mut kept = scores(&lists, 0);
    kept.sort();
    assert_eq!(kept, vec![(2, vec![42]), (3, vec![43])]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c:520-539: the pipe
// sorts the lists by e-value first, so a better list met later removes the HSPs inside it of
// a worse list met earlier, with the overhang of the options; the new hit list keeps the
// query index.
#[test]
fn the_pipe_reads_the_lists_by_evalue() {
    let full = |qs: i32, qe: i32, score: i32, e_value: f64, oid: u32| {
        let mut hsp = super::super::blastp_hsp_from_preliminary(&prelim(3, qs, qe, score), e_value);
        hsp.s_idx = oid;
        hsp
    };
    let mut hit_list = BlastpHitList::new(10);
    hit_list.update(BlastpHspList {
        oid: 7,
        query_index: 3,
        hsps: vec![full(100, 200, 80, 1e-12, 7)],
        best_evalue: 0.0,
    });
    hit_list.update(BlastpHspList {
        oid: 4,
        query_index: 3,
        hsps: vec![full(95, 205, 300, 1e-50, 4), full(400, 500, 70, 1e-11, 4)],
        best_evalue: 0.0,
    });
    best_hit_traceback(&mut hit_list, params(), 10);
    let kept: Vec<(u32, u32, Vec<i32>)> = hit_list
        .hsplist_array
        .iter()
        .map(|l| {
            (
                l.oid,
                l.query_index,
                l.hsps.iter().map(|h| h.raw_score).collect(),
            )
        })
        .collect();
    assert_eq!(kept, vec![(4, 3, vec![300, 70])]);
    assert_eq!(hit_list.hsplist_count, 1);
}
