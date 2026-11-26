use super::super::empty_rows::SubjectLoc;
use super::super::input::FastaRecord;
use super::super::runtime::BatchResults;
use super::*;
use crate::common::GapEditOp::{Del, Ins, Sub};

fn hsp(evalue: f64, edit_script: Vec<GapEditOp>) -> Hsp {
    let [q, s] = edit_script.iter().fold([0, 0], |[q, s], op| match *op {
        Sub(n) => [q + n as i32, s + n as i32],
        Ins(n) => [q + n as i32, s],
        Del(n) => [q, s + n as i32],
    });
    Hsp {
        hsp: PreliminaryHsp {
            context: 0,
            frame: 1,
            score: 40,
            q_start: 0,
            q_end: q,
            q_gapped_start: 0,
            s_start: 5,
            s_end: 5 + s,
            s_gapped_start: 5,
        },
        num: 1,
        evalue,
        bit_score: 0.0,
        identity: 0,
        positive: 0,
        edit_script,
        composition_method: 0,
        query_left_out: 0,
    }
}

fn batch(hsps: Vec<Hsp>) -> BatchResults {
    let record = FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence: vec![b'A'; 60],
        lowercase_masks: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    BatchResults {
        query_ordinal: 0,
        in_report_order: false,
        queries: vec![vec![HspList {
            query_index: 0,
            oid: 0,
            hsps,
            best_evalue: 0.0,
        }]],
        prepared: super::super::query_setup::prepare_queries(
            &[record],
            &super::super::stage_c_tests::options(),
        )
        .unwrap(),
        parameters: Vec::new(),
        search_skipped: false,
        request_space: None,
        unlisted: Vec::new(),
        unlisted_subject_row_throws: false,
        shown: Vec::new(),
        gapped_k: None,
        fetched_utf8: None,
    }
}

/// (start, end) of the query and the subject, the left-out query letters and the script of each
/// HSP of the batch, in the list's order.
fn hsps(batch: &BatchResults) -> Vec<((i32, i32), (i32, i32), u32, Vec<GapEditOp>)> {
    batch.queries[0][0]
        .hsps
        .iter()
        .map(|h| {
            (
                (h.hsp.q_start, h.hsp.q_end),
                (h.hsp.s_start, h.hsp.s_end),
                h.query_left_out,
                h.edit_script.clone(),
            )
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1632-1638
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1568-1580
/// `-subject`'s records are sorted, then mapped once: an end segment with a gap row goes (its
/// letters move the starts; the query's are kept as `query_left_out` for the coverage), and a
/// second call, which a report of several formats makes, neither sorts nor maps again (the script
/// `[Del, Ins, ..]` would lose its `Ins` too). `-db`'s records are sorted only.
#[test]
fn put_in_report_order_maps_subject_records_once() {
    let made = || {
        vec![
            hsp(1e-3, vec![Del(3), Ins(2), Sub(10)]),
            hsp(1e-5, vec![Ins(4), Sub(10), Del(1)]),
        ]
    };
    let mapped = vec![
        ((4, 14), (5, 15), 4, vec![Sub(10)]),
        ((0, 12), (8, 18), 0, vec![Ins(2), Sub(10)]),
    ];
    let mut subject = batch(made());
    for _ in 0..2 {
        subject.put_in_report_order(SubjectLoc::Interval);
        assert_eq!(hsps(&subject), mapped);
    }

    let mut db = batch(made());
    db.put_in_report_order(SubjectLoc::Whole);
    assert_eq!(
        hsps(&db),
        vec![
            ((0, 14), (5, 16), 0, vec![Ins(4), Sub(10), Del(1)]),
            ((0, 12), (5, 18), 0, vec![Del(3), Ins(2), Sub(10)]),
        ]
    );
}
