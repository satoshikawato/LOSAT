use super::*;
use crate::algorithm::blastp::encoding::EncodedProtein;

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:806-833 (cited at
// `final_hsp`) and blast_gapalign.c:4719-4775 (`BLAST_GetUngappedHSPList`: an ungapped HSP is
// one run of aligned letters from `q_start`, its score the ungapped score). The final record
// of an ungapped HSP counts the identities and the positives of its one run on the unmasked
// query, with no gap.
#[test]
fn ungapped_hsp_counts_identities_and_positives_of_one_run() {
    // NCBIstdaa S is 17 (BLOSUM62: A/S 1, W/C -2).
    const S: u8 = 17;
    let query = [W, A, C, A, W];
    let subject = [C, W, A, W, S, A];
    let mut prelim = hsp((0, 4), (1, 5));
    prelim.raw_score = 7;
    let record = final_hsp(
        &prelim,
        vec![GapEditOp::Sub(4)],
        &query,
        &subject,
        ScoringMatrix::Blosum62,
    );
    // W/W A/A: identities; C/W: none; A/S: positive.
    assert_eq!((record.q_start, record.q_end), (1, 4));
    assert_eq!((record.s_start, record.s_end), (2, 5));
    assert_eq!((record.length, record.gapopen), (4, 0));
    assert_eq!((record.num_ident, record.num_positives), (2, 3));
    assert_eq!(record.mismatch, 2);
    assert_eq!(record.raw_score, 7);
    assert!((record.identity - 50.0).abs() < 1e-12);
}

/// A context of the query `letters` (NCBIstdaa, without sentinels) with BLOSUM62's ungapped
/// block, as an ungapped search's contexts hold it.
fn context_of(letters: &[u8]) -> QueryContext {
    let mut aa_seq = vec![0u8];
    aa_seq.extend_from_slice(letters);
    aa_seq.push(0);
    QueryContext {
        q_idx: 0,
        frame: 0,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: letters.len(),
        orig_len: letters.len(),
        frame_base: 0,
        is_valid: true,
        karlin_params: crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62),
    }
}

fn encoded(letters: &[u8]) -> EncodedProtein {
    let mut aa_seq = vec![0u8];
    aa_seq.extend_from_slice(letters);
    aa_seq.push(0);
    EncodedProtein {
        aa_seq,
        aa_len: letters.len(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1685-1704 and
// c++/src/algo/blast/core/blast_hits.c:1918-1928 (cited at `UngappedFinish`): an ungapped
// search's final HSP is the preliminary one with its bit score from the context's own block
// and the E-value that the preliminary stage gave it; an HSP of length 0 is kept.
#[test]
fn ungapped_finish_keeps_the_evalue_and_scores_with_the_contexts_block() {
    let contexts = [context_of(&[W, A, C, A, W])];
    let subjects = [encoded(&[C, W, A, C, A, W])];
    let finish = UngappedFinish {
        contexts: &contexts,
        subjects: XRedoSubjects::Encoded(&subjects),
        matrix: ScoringMatrix::Blosum62,
        hitlist_size: 10,
        sorts_results: true,
    };
    let mut long = hsp((0, 5), (1, 6));
    long.raw_score = 28;
    let mut empty = hsp((2, 2), (3, 3));
    empty.raw_score = 4;
    let list = finish
        .finish_list(0, 0, 3, &[(long, 0.125), (empty, 0.5)])
        .expect("the list is finished")
        .expect("a list with HSPs");
    assert_eq!((list.oid, list.query_index), (0, 3));
    assert_eq!(list.hsps.len(), 2);
    let karlin = contexts[0].karlin_params;
    let (first, second) = (&list.hsps[0], &list.hsps[1]);
    assert_eq!(first.e_value, 0.125);
    assert_eq!(second.e_value, 0.5);
    let bits = |score: i32| {
        (f64::from(score) * karlin.lambda - karlin.k.glibc_ln()) / std::f64::consts::LN_2
    };
    assert!(
        (first.bit_score - bits(28)).abs() < 1e-9,
        "{}",
        first.bit_score
    );
    assert!(
        (second.bit_score - bits(4)).abs() < 1e-9,
        "{}",
        second.bit_score
    );
    assert_eq!((first.q_start, first.q_end, first.length), (1, 5, 5));
    assert_eq!((first.num_ident, first.gapopen), (5, 0));
    // Length 0: no letter, the query range [2, 2) written as q_start 3, q_end 2.
    assert_eq!((second.q_start, second.q_end, second.length), (3, 2, 0));
    assert_eq!(second.num_ident, 0);
    assert_eq!(second.identity, 0.0);
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
