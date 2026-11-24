use super::*;
use crate::algorithm::blastn::interval_tree::IndexMethod;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819-834
#[test]
fn containment_groups_three_frames_by_query_strand() {
    let (_, batch, _) = super::super::stage_c_tests::fixture("S06_seed.fna");
    let h = PreliminaryHsp {
        context: 0,
        frame: 1,
        score: 100,
        q_start: 0,
        q_end: 70,
        q_gapped_start: 2,
        s_start: 0,
        s_end: 70,
        s_gapped_start: 2,
    };
    let mut tree = BlastIntervalTree::new(0, 600, 0, 100);
    let node = tree_hsp(&h, &batch);
    tree.add_hsp(
        node,
        node.query_context_offset,
        IndexMethod::QueryAndSubject,
    );
    for c in 0..6 {
        let mut candidate = h.clone();
        candidate.context = c;
        candidate.frame = batch.contexts[c].frame;
        candidate.score = 90;
        candidate.q_start = 5;
        candidate.q_end = 60;
        candidate.s_start = 5;
        candidate.s_end = 60;
        let node = tree_hsp(&candidate, &batch);
        assert_eq!(
            tree.contains_hsp(&node, node.query_context_offset, 0),
            c < 3
        );
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1358
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1687-1706
#[test]
fn capped_save_heap_ties_preserve_ncbi_comparator() {
    let h = |score, s_start| PreliminaryHsp {
        context: 0,
        frame: 1,
        score,
        q_start: 0,
        q_end: 5,
        q_gapped_start: 1,
        s_start,
        s_end: s_start + 5,
        s_gapped_start: s_start + 1,
    };
    let mut hits = Vec::new();
    let mut heap = false;
    for x in [h(30, 5), h(20, 0), h(30, 1), h(10, 0)] {
        save_hsp(&mut hits, x, 2, &mut heap);
    }
    sort_by_score(&mut hits);
    assert_eq!(hits, vec![h(30, 1), h(30, 5)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3405-3437 (cited at
// `gapped_start_negative`)
/// The gapped start of a one-hit HSP of negative length: `q_start + 5` for a window inside both
/// sequences that sums above 0, `q_start` for one that sums to 0 or less or holds a sequence's
/// NULLB, whatever follows it, and so also for a window that reads past NCBI's buffers, whose
/// letters past the sequence read as the sentinel (PD-LOSAT-NCBI-DEFECTS exception 3).
#[test]
fn a_one_hit_gapped_start_of_negative_length_scores_its_first_window() {
    use crate::utils::matrix::aa_char_to_ncbistdaa;
    let start = |query: &[u8], subject: &[u8], q_start: usize, s_start: usize| {
        gapped_start_negative(query, subject, q_start, s_start, ScoringMatrix::Blosum62)
    };
    let code =
        |letters: &[u8]| -> Vec<u8> { letters.iter().map(|&l| aa_char_to_ncbistdaa(l)).collect() };
    // BLOSUM62: W/W 11, W/G -2.
    let query = code(b"MKWWWWWWWWWWWK");
    let subject = code(b"WWWWWWWWWWWWW");
    assert_eq!(start(&query, &subject, 2, 0), 2 + 5);
    // The subject ends one letter before the window's end: its NULLB is in the window.
    assert_eq!(start(&query, &subject[..10], 2, 0), 2);
    // The window passes the query's end (its context's NULLB).
    assert_eq!(start(&query, &subject, 4, 0), 4);
    // Inside both, sums of -22, -9 (one W/W) and 4 (two).
    assert_eq!(start(&query, &code(b"GGGGGGGGGGGG"), 2, 0), 2);
    assert_eq!(start(&query, &code(b"WGGGGGGGGGGG"), 2, 0), 2);
    assert_eq!(start(&query, &code(b"WWGGGGGGGGGG"), 2, 0), 2 + 5);
    // From s_start on the subject side.
    assert_eq!(start(&query, &code(b"GGWWGGGGGGGGGG"), 2, 2), 2 + 5);
    // Windows whose letters lie past NCBI's buffers (the query batch's last context near its
    // end, a `-subject` protein near its end, every matrix): the letters past the sequence read
    // as NULLBs, as the first one past it does, and the start is the HSP's.
    for matrix in ScoringMatrix::ALL {
        let wide = code(b"WWWWWWWWWWWWWWWWWWWW");
        for (q_start, s_start) in [(19, 0), (0, 19), (13, 9), (19, 19)] {
            assert_eq!(
                gapped_start_negative(&wide, &wide, q_start, s_start, matrix),
                q_start,
                "{matrix:?} {q_start} {s_start}"
            );
        }
        assert_eq!(gapped_start_negative(&wide, &wide, 9, 9, matrix), 9 + 5);
    }
}
