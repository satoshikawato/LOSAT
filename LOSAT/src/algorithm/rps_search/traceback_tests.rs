use super::*;

fn prelim_list(oid: u32, score: i32) -> BlastpHspList<XPrelimHsp> {
    let hsp = crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp {
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
    };
    BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![XPrelimHsp { hsp, e_value: 0.0 }],
        best_evalue: 0.0,
    }
}

// blast_hspstream.c:676-686 and 196-203 (cited at `stream_lists`): only rpsblast's redo
// reads a query's lists by score; rpstblastn's redo and every traceback without
// composition-based statistics read them by ordinal id, so the last list, whose K the
// batch's first footer may print, is the highest profile.
#[test]
fn the_stream_reads_by_score_only_for_a_protein_querys_redo() {
    let mut prelim = BlastpHitList::new(10);
    for (oid, score) in [(0, 50), (3, 80), (1, 80), (2, 65)] {
        prelim.update(prelim_list(oid, score));
    }
    let oids = |by_score| -> Vec<u32> {
        stream_lists(prelim.clone(), by_score)
            .iter()
            .map(|list| list.oid)
            .collect()
    };
    assert_eq!(oids(true), [3, 1, 2, 0]);
    assert_eq!(oids(false), [0, 1, 2, 3]);
}

// blast_traceback.c:346-354 (cited at `evalue_lengths`, DECISIONS 46): with `isRPS`, `m` is
// the traceback's `stat_length` and `n` the HSP's context of the user query, which the
// shallow copy of the query info turns into the profile's length for context 0 (frame +1).
#[test]
fn a_translated_lists_evalues_read_the_profile_length_for_frame_plus_1_only() {
    // rpsblast: the query's length and the profile's.
    assert_eq!(evalue_lengths(false, 0, 200, 180, 999), (200, 180));
    // rpstblastn: m = stat_length (the frame length of the list's last HSP aligned).
    assert_eq!(evalue_lengths(true, 1, 190, 180, 189), (189, 180));
    for (frame, length) in [(2, 189), (3, 188), (-1, 190), (-2, 189), (-3, 188)] {
        assert_eq!(
            evalue_lengths(true, frame, length, 180, 189),
            (189, length),
            "{frame}"
        );
    }
}
