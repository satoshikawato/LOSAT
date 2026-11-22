use super::*;
use crate::algorithm::blastp::hsp::{BlastpHsp, BlastpHspList};
use crate::blastinput::seq_range::RecordPlacement;
use crate::common::Hit;

/// An HSP of subject `s_idx` whose query range is `q` and subject range `s` (1-based, closed;
/// `(n + 1, n)` is the range of length 0 at offset `n`).
fn hsp(q: (usize, usize), s: (usize, usize), s_idx: u32) -> BlastpHsp {
    BlastpHsp::from_hit(Hit {
        identity: 0.0,
        length: s.1 + 1 - s.0,
        mismatch: 0,
        gapopen: 0,
        q_start: q.0,
        q_end: q.1,
        s_start: s.0,
        s_end: s.1,
        e_value: 1.0,
        bit_score: 10.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx,
        raw_score: 20,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 0,
    })
}

fn lists(hsps: Vec<BlastpHsp>) -> Vec<Option<BlastpHitList>> {
    let mut hit_list = BlastpHitList::new(10);
    let oid = hsps[0].s_idx;
    hit_list.update(BlastpHspList {
        oid,
        query_index: 0,
        hsps,
        best_evalue: 1.0,
    });
    vec![None, Some(hit_list)]
}

// NCBI's exception, from the official 2.17.0 binary (`-subject` with an HSP of length 0 at
// offset 0 of the subject): exit 255 and this text after `Error: `.
#[test]
fn a_length_zero_hsp_at_the_start_of_a_subject_throws_the_mappers_exception() {
    let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
    let exception = subject_start_exception(
        &lists(vec![hsp((5, 4), (1, 0), 0)]),
        &Placements::default(),
        true,
    )
    .expect("an HSP of length 0 at the first letter");
    // LOSAT's words of NCBI's text (`crate::wording`, as `what` gives it).
    let ncbi = format!(
        "NCBI C++ Exception:\n    T0 \"{root}/src/objects/seq/seq_align_mapper_base.cpp\", line 1292: Error: (CAnnotMapperException::eBadAlignment) Mapped alignment contains gaps and can not be converted to dense-diag.\n"
    );
    assert_eq!(
        exception.what(),
        crate::wording::losat_bytes(ncbi.as_bytes()).into_owned()
    );
}

#[test]
fn only_a_shown_length_zero_hsp_at_the_first_letter_of_a_searched_subject_throws() {
    let none = Placements::default();
    let at = |hsps: Vec<BlastpHsp>, placements: &Placements, subject_set: bool| {
        subject_start_exception(&lists(hsps), placements, subject_set).is_some()
    };
    // Not at the second letter, not a letter of an ordinary HSP, not with a database.
    assert!(!at(vec![hsp((5, 4), (2, 1), 0)], &none, true));
    assert!(!at(vec![hsp((1, 10), (1, 10), 0)], &none, true));
    assert!(!at(vec![hsp((5, 4), (1, 0), 0)], &none, false));
    // The query's first letter does not matter by itself.
    assert!(!at(vec![hsp((1, 0), (2, 1), 0)], &none, true));
    assert!(at(vec![hsp((1, 0), (1, 0), 0)], &none, true));
    // `-subject_loc 24-189` of subject 0: record offset 23 is relative offset 0.
    let ranged = Placements(vec![RecordPlacement {
        offset: 23,
        length: 200,
    }]);
    assert!(at(vec![hsp((5, 4), (24, 23), 0)], &ranged, true));
    assert!(!at(vec![hsp((5, 4), (25, 24), 0)], &ranged, true));
    assert!(!at(vec![hsp((5, 4), (1, 0), 0)], &ranged, true));
}
