use super::*;

// NCBIstdaa: A 1, C 3, W 20 (BLOSUM62: W/W 11, A/A 4, W/C -2, A/W -3).
const A: u8 = 1;
const C: u8 = 3;
const W: u8 = 20;

fn hsp(q: (i32, i32), s: (i32, i32)) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: q.0,
        query_end: q.1,
        subject_start: s.0,
        subject_end: s.1,
        gapped_query_start: 0,
        gapped_subject_start: 0,
        raw_score: 0,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
    }
}

#[test]
fn short_hsp_starts_at_its_middle() {
    let seq = [W; 20];
    assert_eq!(
        get_offsets_for_gapped_alignment(
            &seq,
            &seq,
            &hsp((3, 14), (5, 16)),
            ScoringMatrix::Blosum62
        ),
        Some((8, 10))
    );
}

#[test]
fn long_hsp_starts_at_the_end_of_its_first_best_window() {
    // Every window scores 121: the first one's last letter.
    let seq = [W; 30];
    assert_eq!(
        get_offsets_for_gapped_alignment(
            &seq,
            &seq,
            &hsp((2, 22), (4, 24)),
            ScoringMatrix::Blosum62
        ),
        Some((12, 14))
    );
    // A run of A then W: the windows grow until the first all-W one (ending at 21); the
    // later ones only tie.
    let mut seq = [A; 30];
    seq[11..].fill(W);
    assert_eq!(
        get_offsets_for_gapped_alignment(
            &seq,
            &seq,
            &hsp((0, 25), (0, 25)),
            ScoringMatrix::Blosum62
        ),
        Some((21, 21))
    );
}

#[test]
fn hsp_without_a_positive_window_has_no_start() {
    let query = [W; 30];
    let subject = [C; 30];
    assert_eq!(
        get_offsets_for_gapped_alignment(
            &query,
            &subject,
            &hsp((0, 20), (0, 20)),
            ScoringMatrix::Blosum62
        ),
        None
    );
}
