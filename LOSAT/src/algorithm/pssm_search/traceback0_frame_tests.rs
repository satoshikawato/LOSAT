use super::*;

fn hsp(frame: i32, query: (i32, i32), subject: (i32, i32), score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: query.0,
        query_end: query.1,
        subject_start: subject.0,
        subject_end: subject.1,
        gapped_query_start: query.0,
        gapped_subject_start: subject.0,
        raw_score: score,
        query_context: 0,
        query_frame: 0,
        subject_frame: frame,
        query_length: 100,
        q_idx: 0,
        s_idx: 0,
    }
}

// blast_traceback.c:402-405 (cited at `trace_hsps`): before each alignment of a translated
// subject, an HSP inside one traced before is dropped when their subject frames have the
// same sign (`BlastIntervalTreeContainsHSP` compares the signs only, blast_itree.c), so a
// frame +2 HSP inside a frame +1 one is contained and a frame -1 one is not.
#[test]
fn containment_reads_the_sign_of_the_subject_frame_only() {
    let mut tree = BlastIntervalTree::new(0, 101, 0, 101);
    tree.add_hsp(
        tree_hsp(&hsp(1, (10, 60), (10, 60), 500)),
        0,
        IndexMethod::QueryAndSubject,
    );
    assert!(tree.contains_hsp(&tree_hsp(&hsp(2, (20, 40), (20, 40), 100)), 0, 0));
    assert!(!tree.contains_hsp(&tree_hsp(&hsp(-1, (20, 40), (20, 40), 100)), 0, 0));
    assert!(!tree.contains_hsp(&tree_hsp(&hsp(1, (5, 40), (20, 40), 100)), 0, 0));
}

// blast_traceback.c:416-421 (cited at `TraceSubject`): each HSP aligns on its own frame.
#[test]
fn a_translated_subject_gives_each_hsp_its_frames_letters() {
    let frames: [&[u8]; 6] = [b"a", b"b", b"c", b"d", b"e", b"f"];
    let subject = TraceSubject::Frames {
        frames: &frames,
        tree_bound: 3,
    };
    let letters: Vec<&[u8]> = [1, 2, 3, -1, -2, -3]
        .map(|frame| subject.of(frame))
        .to_vec();
    assert_eq!(letters, frames.to_vec());
    assert_eq!(TraceSubject::One(b"xyz").of(0), b"xyz");
}
