//! Unit tests for tblastx/lookup.rs

use crate::algorithm::tblastx::lookup::{build_ncbi_lookup, pv_test};
use crate::algorithm::tblastx::translation::generate_frames;
use crate::stats::KarlinParams;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::ncbistdaa;

#[test]
fn test_build_ncbi_lookup_indexes_low_self_score_exact_word() {
    // NCBI behavior: even if a word's self-score is below the neighborhood threshold,
    // the exact word occurrences must still be indexed (blast_aalookup.c s_AddWordHits()).
    //
    // Example: AAA has self-score 4+4+4 = 12 which is < 13 (tblastx default threshold),
    // so it must still appear in the lookup table as an exact match.
    let code = GeneticCode::from_id(1);

    // 5 alanines (A) => 3 occurrences of AAA as a 3-mer.
    let mut dna = Vec::new();
    for _ in 0..5 {
        dna.extend_from_slice(b"GCT"); // Ala in standard code
    }

    let mut frames = generate_frames(&dna, &code);
    frames.retain(|f| f.frame == 1);
    assert_eq!(frames.len(), 1);

    let queries = vec![frames];
    let karlin = KarlinParams::default();
    let (lookup, _ctx) = build_ncbi_lookup(&queries, 13, &karlin, true);

    // NCBISTDAA: A=1, so AAA index = (1<<10)|(1<<5)|1 = 1057.
    let aaa_idx: usize =
        ((ncbistdaa::A as usize) << 10) | ((ncbistdaa::A as usize) << 5) | (ncbistdaa::A as usize);
    assert!(pv_test(&lookup.pv, aaa_idx));
    let hits = lookup.get_hits(aaa_idx);
    assert_eq!(hits.len(), 3);
}

#[test]
fn test_build_ncbi_lookup_longest_chain_tracks_max_cell() {
    // longest_chain should equal the max num_used over all cells (NCBI Finalize semantics).
    let code = GeneticCode::from_id(1);

    // 100 alanines => 98 AAA 3-mers in frame 1.
    let mut dna = Vec::new();
    for _ in 0..100 {
        dna.extend_from_slice(b"GCT");
    }

    let mut frames = generate_frames(&dna, &code);
    frames.retain(|f| f.frame == 1);
    assert_eq!(frames.len(), 1);

    let queries = vec![frames];
    let karlin = KarlinParams::default();
    let (lookup, _ctx) = build_ncbi_lookup(&queries, 13, &karlin, true);

    // AAA index in NCBISTDAA (A=1) is 1057.
    let aaa_idx: usize =
        ((ncbistdaa::A as usize) << 10) | ((ncbistdaa::A as usize) << 5) | (ncbistdaa::A as usize);
    let hits = lookup.get_hits(aaa_idx);
    assert_eq!(lookup.longest_chain as usize, hits.len());
    assert_eq!(lookup.longest_chain, 98);
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/queryinfo_unit_test.cpp:174-180 BlastnSearchContextInfo
#[test]
fn test_get_context_idx_matches_ncbi() {
    // Test get_context_idx against NCBI BSearchContextInfo behavior
    // Reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/queryinfo_unit_test.cpp:174-180
    //
    // NCBI test cases:
    // - BSearchContextInfo(0, query_info) → 0 (first context)
    // - BSearchContextInfo(length_1 + 1, query_info) → 1 (second context)
    // - BSearchContextInfo(2 * (length_1 + 1), query_info) → 2 (third context)
    //
    // For tblastx, we have 6 frames per query, so contexts are:
    // - Context 0: query 0, frame 0 (base = 0)
    // - Context 1: query 0, frame 1 (base = frame_0_len - 1)
    // - Context 2: query 0, frame 2 (base = frame_0_len - 1 + frame_1_len - 1)
    // - ...

    let code = GeneticCode::from_id(1);
    let karlin = KarlinParams::default();

    // Create two queries with different lengths
    let dna_seq1 = b"ATGCGATCGATCGATCGATGCGATCGATCGATCG"; // 36 bases = 12 codons
    let dna_seq2 = b"ATGCGATCGATCG"; // 14 bases = 4 codons (shorter)

    let frames1 = generate_frames(dna_seq1, &code);
    let frames2 = generate_frames(dna_seq2, &code);

    let queries = vec![frames1, frames2];
    let (lookup, contexts) = build_ncbi_lookup(&queries, 13, &karlin, true);

    // Verify frame_bases structure
    assert_eq!(lookup.frame_bases.len(), contexts.len());
    assert_eq!(lookup.num_contexts, contexts.len());

    // Test 1: Offset 0 should return context 0 (first frame of first query)
    let ctx_idx_0 = lookup.get_context_idx(0);
    assert_eq!(ctx_idx_0, 0, "Offset 0 should return context 0");
    assert_eq!(contexts[ctx_idx_0].q_idx, 0);

    // Test 2: Offset at first frame base should return context 0
    let frame_0_base = lookup.frame_bases[0];
    let ctx_idx_frame0 = lookup.get_context_idx(frame_0_base);
    assert_eq!(
        ctx_idx_frame0, 0,
        "Offset at first frame base should return context 0"
    );

    // Test 3: Offset at second frame base should return context 1
    if lookup.frame_bases.len() > 1 {
        let frame_1_base = lookup.frame_bases[1];
        let ctx_idx_frame1 = lookup.get_context_idx(frame_1_base);
        assert_eq!(
            ctx_idx_frame1, 1,
            "Offset at second frame base should return context 1"
        );
        assert_eq!(contexts[ctx_idx_frame1].q_idx, 0);
    }

    // Test 4: Offset just before second frame base should return context 0
    if lookup.frame_bases.len() > 1 {
        let frame_1_base = lookup.frame_bases[1];
        let ctx_idx_before_frame1 = lookup.get_context_idx(frame_1_base - 1);
        assert_eq!(
            ctx_idx_before_frame1, 0,
            "Offset just before second frame base should return context 0"
        );
    }

    // Test 5: Offset at last frame base should return last context
    let last_base = lookup.frame_bases[lookup.frame_bases.len() - 1];
    let ctx_idx_last = lookup.get_context_idx(last_base);
    assert_eq!(
        ctx_idx_last,
        lookup.frame_bases.len() - 1,
        "Offset at last frame base should return last context"
    );

    // Test 6: Offset beyond last frame base should return last context (saturating)
    let beyond_last = last_base + 1000;
    let ctx_idx_beyond = lookup.get_context_idx(beyond_last);
    assert_eq!(
        ctx_idx_beyond,
        lookup.frame_bases.len() - 1,
        "Offset beyond last frame base should return last context"
    );

    // Test 7: Verify all frame bases are in ascending order (required for binary search)
    for i in 1..lookup.frame_bases.len() {
        assert!(
            lookup.frame_bases[i] >= lookup.frame_bases[i - 1],
            "Frame bases must be in ascending order for binary search"
        );
    }

    // Test 8: For each context, verify get_context_idx returns correct index
    for (idx, &base) in lookup.frame_bases.iter().enumerate() {
        let ctx_idx = lookup.get_context_idx(base);
        assert_eq!(
            ctx_idx, idx,
            "get_context_idx should return correct index for frame base {}",
            base
        );
    }
}

#[test]
fn test_get_context_idx_edge_cases() {
    // Test edge cases for get_context_idx
    let code = GeneticCode::from_id(1);
    let karlin = KarlinParams::default();

    // Single query with single frame (minimal case)
    let dna_seq = b"ATGCGATCG"; // 9 bases = 3 codons
    let mut frames = generate_frames(dna_seq, &code);
    frames.retain(|f| f.frame == 1); // Keep only frame 1
    assert_eq!(frames.len(), 1);

    let queries = vec![frames];
    let (lookup, contexts) = build_ncbi_lookup(&queries, 13, &karlin, true);

    // With single context, all offsets should return 0
    assert_eq!(lookup.get_context_idx(0), 0);
    assert_eq!(lookup.get_context_idx(100), 0);
    assert_eq!(lookup.get_context_idx(-100), 0); // saturating_sub(1) handles negative

    // Verify context structure
    assert_eq!(contexts.len(), 1);
    assert_eq!(contexts[0].q_idx, 0);
}

#[test]
fn test_get_context_idx_multiple_queries() {
    // Test get_context_idx with multiple queries (similar to NCBI's multi-query test)
    let code = GeneticCode::from_id(1);
    let karlin = KarlinParams::default();

    let dna_seq1 = b"ATGCGATCGATCGATCG"; // 18 bases
    let dna_seq2 = b"ATGCGATCGATCGATCGATGCGATCGATCGATCG"; // 36 bases (longer)

    let frames1 = generate_frames(dna_seq1, &code);
    let frames2 = generate_frames(dna_seq2, &code);

    let queries = vec![frames1, frames2];
    let (lookup, contexts) = build_ncbi_lookup(&queries, 13, &karlin, true);

    // Should have 12 contexts (6 frames × 2 queries)
    assert_eq!(contexts.len(), 12);
    assert_eq!(lookup.frame_bases.len(), 12);

    // Verify query indices are correct
    for (idx, ctx) in contexts.iter().enumerate() {
        let expected_q_idx = (idx / 6) as u32;
        assert_eq!(
            ctx.q_idx, expected_q_idx,
            "Context {} should belong to query {}",
            idx, expected_q_idx
        );
    }

    // Test: Offset at start of second query (context 6) should return context 6
    if lookup.frame_bases.len() > 6 {
        let query_1_base = lookup.frame_bases[6];
        let ctx_idx = lookup.get_context_idx(query_1_base);
        assert_eq!(
            ctx_idx, 6,
            "Offset at start of second query should return context 6"
        );
        assert_eq!(contexts[ctx_idx].q_idx, 1);
    }
}
