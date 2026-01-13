//! Unit tests for tblastx/lookup.rs

use LOSAT::algorithm::tblastx::lookup::{encode_aa_kmer, build_direct_lookup};
use LOSAT::algorithm::tblastx::lookup::{build_ncbi_lookup, pv_test, BlastAaLookupTable};
use LOSAT::algorithm::tblastx::translation::generate_frames;
use LOSAT::utils::genetic_code::GeneticCode;
use LOSAT::utils::matrix::ncbistdaa;
use LOSAT::stats::KarlinParams;

#[test]
fn test_encode_aa_kmer_basic() {
    // Sequences are encoded in NCBISTDAA (0-27). Example: A=1, R=16, N=13.
    // Lookup index uses bit-shift encoding (charsize=5 for alphabet_size=28).
    let seq = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N];
    let code = encode_aa_kmer(&seq, 0);
    let expected = ((ncbistdaa::A as usize) << 10)
        | ((ncbistdaa::R as usize) << 5)
        | (ncbistdaa::N as usize);
    assert_eq!(code, Some(expected));
}

#[test]
fn test_encode_aa_kmer_different_positions() {
    // A R N D C Q in NCBISTDAA encoding.
    let seq = [
        ncbistdaa::A,
        ncbistdaa::R,
        ncbistdaa::N,
        ncbistdaa::D,
        ncbistdaa::C,
        ncbistdaa::Q,
    ];
    
    // ARN
    let code1 = encode_aa_kmer(&seq, 0);
    let expected1 = ((ncbistdaa::A as usize) << 10)
        | ((ncbistdaa::R as usize) << 5)
        | (ncbistdaa::N as usize);
    assert_eq!(code1, Some(expected1));
    
    // RND
    let code2 = encode_aa_kmer(&seq, 1);
    let expected2 = ((ncbistdaa::R as usize) << 10)
        | ((ncbistdaa::N as usize) << 5)
        | (ncbistdaa::D as usize);
    assert_eq!(code2, Some(expected2));
    
    // NDC
    let code3 = encode_aa_kmer(&seq, 2);
    let expected3 = ((ncbistdaa::N as usize) << 10)
        | ((ncbistdaa::D as usize) << 5)
        | (ncbistdaa::C as usize);
    assert_eq!(code3, Some(expected3));
}

#[test]
fn test_encode_aa_kmer_with_stop_codon() {
    // NCBI tblastx uses BLASTAA_SIZE=28 lookup alphabet and includes STOP ('*') in encoding.
    // Stop codon is 25 in NCBISTDAA.
    let seq_with_stop = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::STOP];
    let code = encode_aa_kmer(&seq_with_stop, 0);
    let expected = ((ncbistdaa::A as usize) << 10)
        | ((ncbistdaa::R as usize) << 5)
        | (ncbistdaa::STOP as usize);
    assert_eq!(code, Some(expected));
    
    // Stop codon in middle
    let seq_with_stop_middle = [ncbistdaa::A, ncbistdaa::STOP, ncbistdaa::N];
    let code = encode_aa_kmer(&seq_with_stop_middle, 0);
    let expected = ((ncbistdaa::A as usize) << 10)
        | ((ncbistdaa::STOP as usize) << 5)
        | (ncbistdaa::N as usize);
    assert_eq!(code, Some(expected));
    
    // Stop codon at start
    let seq_with_stop_start = [ncbistdaa::STOP, ncbistdaa::R, ncbistdaa::N];
    let code = encode_aa_kmer(&seq_with_stop_start, 0);
    let expected = ((ncbistdaa::STOP as usize) << 10)
        | ((ncbistdaa::R as usize) << 5)
        | (ncbistdaa::N as usize);
    assert_eq!(code, Some(expected));
}

#[test]
fn test_encode_aa_kmer_out_of_bounds() {
    let seq = [0u8, 1u8, 2u8];
    
    // Try to encode beyond sequence length
    let code = encode_aa_kmer(&seq, 1);
    assert_eq!(code, None);
    
    // Try to encode at position that would exceed bounds
    let code = encode_aa_kmer(&seq, 2);
    assert_eq!(code, None);
}

#[test]
fn test_encode_aa_kmer_edge_cases() {
    // Empty sequence
    let empty_seq: [u8; 0] = [];
    let code = encode_aa_kmer(&empty_seq, 0);
    assert_eq!(code, None);
    
    // Sequence too short
    let short_seq = [0u8, 1u8];
    let code = encode_aa_kmer(&short_seq, 0);
    assert_eq!(code, None);
    
    // Exactly 3 amino acids
    let exact_seq = [0u8, 1u8, 2u8];
    let code = encode_aa_kmer(&exact_seq, 0);
    assert!(code.is_some());
}

#[test]
fn test_encode_aa_kmer_max_values() {
    // Test with high-value amino acids in NCBISTDAA encoding
    // NCBISTDAA encoding (blast_encoding.c:115-118):
    // '-','A','B','C','D','E','F','G','H','I','K','L','M',
    // 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'
    // X=21 (Unknown), Y=22, Z=23 (Glu or Gln), U=24 (Selenocysteine)
    // Note: Lookup table uses BLASTAA_SIZE=28, so all residues 0-27 are valid
    // but stop codon (*=25) and other special residues may be filtered
    let seq = [21u8, 22u8, 23u8]; // X, Y, Z
    let code = encode_aa_kmer(&seq, 0);
    // (21<<10) | (22<<5) | 23 = 21,591
    assert_eq!(code, Some((21 << 10) | (22 << 5) | 23));
    
    // Verify it's within backbone size (max index + 1 for 0..=27 with charsize=5)
    assert!(code.unwrap() < 28_540);
}

#[test]
fn test_build_direct_lookup_basic() {
    let code = GeneticCode::from_id(1);
    
    // Create a simple DNA sequence and generate frames
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases = 6 codons
    let frames = generate_frames(dna_seq, &code);
    
    // Build lookup table
    let queries = vec![frames];
    let lookup = build_direct_lookup(&queries);
    
    // Verify table size (NCBI-style backbone size for 0..=23 with charsize=5)
    assert_eq!(lookup.len(), 28_540);
    
    // Check that some k-mers were found
    let total_hits: usize = lookup.iter().map(|bucket| bucket.len()).sum();
    assert!(total_hits > 0);
}

#[test]
fn test_build_direct_lookup_multiple_frames() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG";
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let lookup = build_direct_lookup(&queries);
    
    // Each frame should contribute k-mers
    // With 6 frames and ~6 codons per frame, we should have multiple hits
    let total_hits: usize = lookup.iter().map(|bucket| bucket.len()).sum();
    assert!(total_hits >= 6); // At least one k-mer per frame
}

#[test]
fn test_build_direct_lookup_multiple_queries() {
    let code = GeneticCode::from_id(1);
    
    let dna_seq1 = b"ATGCGATCGATCGATCG";
    let dna_seq2 = b"ATGCGATCGATCGATCG";
    
    let frames1 = generate_frames(dna_seq1, &code);
    let frames2 = generate_frames(dna_seq2, &code);
    
    let queries = vec![frames1, frames2];
    let lookup = build_direct_lookup(&queries);
    
    // Both queries should contribute hits
    let total_hits: usize = lookup.iter().map(|bucket| bucket.len()).sum();
    assert!(total_hits > 0);
    
    // Verify query indices are correct
    for bucket in &lookup {
        for (query_idx, _frame_idx, _aa_pos) in bucket {
            assert!(*query_idx < 2); // Should be 0 or 1
        }
    }
}

#[test]
fn test_build_direct_lookup_with_seg_masking() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases
    let mut frames_nomask = generate_frames(dna_seq, &code);
    // Keep just one frame to make the hit count comparison stable.
    frames_nomask.retain(|f| f.frame == 1);
    assert_eq!(frames_nomask.len(), 1);

    let mut frames_masked = frames_nomask.clone();
    // Mask first few logical AA positions; seeds should not be generated from them.
    frames_masked[0].seg_masks.push((0, 3));

    let lookup_nomask = build_direct_lookup(&vec![frames_nomask]);
    let lookup_masked = build_direct_lookup(&vec![frames_masked]);

    let total_hits_nomask: usize = lookup_nomask.iter().map(|bucket| bucket.len()).sum();
    let total_hits_masked: usize = lookup_masked.iter().map(|bucket| bucket.len()).sum();

    assert!(total_hits_masked <= total_hits_nomask);
}

#[test]
fn test_build_direct_lookup_short_sequence() {
    let code = GeneticCode::from_id(1);
    
    // Very short sequence (only 3 bases = 1 codon)
    let dna_seq = b"ATG";
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let lookup = build_direct_lookup(&queries);
    
    // With only 1 codon, we can't form any 3-mer k-mers
    // So lookup should be mostly empty
    let total_hits: usize = lookup.iter().map(|bucket| bucket.len()).sum();
    assert_eq!(total_hits, 0);
}

#[test]
fn test_build_direct_lookup_frame_indices() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG";
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let lookup = build_direct_lookup(&queries);
    
    // Verify frame indices are valid
    let num_frames = queries[0].len();
    for bucket in &lookup {
        for (_query_idx, frame_idx, _aa_pos) in bucket {
            // Frame index should be less than number of frames
            assert!(*frame_idx < num_frames as u8);
        }
    }
}

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
    let (lookup, _ctx) = build_ncbi_lookup(&queries, 13, true, true, &karlin);

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
    let (lookup, _ctx) = build_ncbi_lookup(&queries, 13, true, true, &karlin);

    // AAA index in NCBISTDAA (A=1) is 1057.
    let aaa_idx: usize =
        ((ncbistdaa::A as usize) << 10) | ((ncbistdaa::A as usize) << 5) | (ncbistdaa::A as usize);
    let hits = lookup.get_hits(aaa_idx);
    assert_eq!(lookup.longest_chain as usize, hits.len());
    assert_eq!(lookup.longest_chain, 98);
}

#[test]
fn test_build_direct_lookup_amino_acid_positions() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases = 6 codons
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let lookup = build_direct_lookup(&queries);
    
    // Verify amino acid positions are LOGICAL (0-indexed, not counting sentinels)
    for bucket in &lookup {
        for (_query_idx, frame_idx, aa_pos) in bucket {
            let frame = &queries[0][*frame_idx as usize];
            // AA position should be within actual amino acid count (not including sentinels)
            assert!((*aa_pos as usize) < frame.aa_len);
            // AA position should allow for 3-mer (pos + 3 <= aa_len)
            assert!((*aa_pos as usize) + 3 <= frame.aa_len);
        }
    }
}

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
    let (lookup, contexts) = build_ncbi_lookup(&queries, 13, true, true, &karlin);
    
    // Verify frame_bases structure
    assert_eq!(lookup.frame_bases.len(), contexts.len());
    assert_eq!(lookup.num_contexts, contexts.len());
    
    // Test 1: Offset 0 should return context 0 (first frame of first query)
    let ctx_idx_0 = lookup.get_context_idx(0);
    assert_eq!(ctx_idx_0, 0, "Offset 0 should return context 0");
    assert_eq!(contexts[ctx_idx_0].q_idx, 0);
    assert_eq!(contexts[ctx_idx_0].f_idx, 0);
    
    // Test 2: Offset at first frame base should return context 0
    let frame_0_base = lookup.frame_bases[0];
    let ctx_idx_frame0 = lookup.get_context_idx(frame_0_base);
    assert_eq!(ctx_idx_frame0, 0, "Offset at first frame base should return context 0");
    
    // Test 3: Offset at second frame base should return context 1
    if lookup.frame_bases.len() > 1 {
        let frame_1_base = lookup.frame_bases[1];
        let ctx_idx_frame1 = lookup.get_context_idx(frame_1_base);
        assert_eq!(ctx_idx_frame1, 1, "Offset at second frame base should return context 1");
        assert_eq!(contexts[ctx_idx_frame1].q_idx, 0);
        assert_eq!(contexts[ctx_idx_frame1].f_idx, 1);
    }
    
    // Test 4: Offset just before second frame base should return context 0
    if lookup.frame_bases.len() > 1 {
        let frame_1_base = lookup.frame_bases[1];
        let ctx_idx_before_frame1 = lookup.get_context_idx(frame_1_base - 1);
        assert_eq!(ctx_idx_before_frame1, 0, "Offset just before second frame base should return context 0");
    }
    
    // Test 5: Offset at last frame base should return last context
    let last_base = lookup.frame_bases[lookup.frame_bases.len() - 1];
    let ctx_idx_last = lookup.get_context_idx(last_base);
    assert_eq!(ctx_idx_last, lookup.frame_bases.len() - 1, 
               "Offset at last frame base should return last context");
    
    // Test 6: Offset beyond last frame base should return last context (saturating)
    let beyond_last = last_base + 1000;
    let ctx_idx_beyond = lookup.get_context_idx(beyond_last);
    assert_eq!(ctx_idx_beyond, lookup.frame_bases.len() - 1,
               "Offset beyond last frame base should return last context");
    
    // Test 7: Verify all frame bases are in ascending order (required for binary search)
    for i in 1..lookup.frame_bases.len() {
        assert!(lookup.frame_bases[i] >= lookup.frame_bases[i-1],
                "Frame bases must be in ascending order for binary search");
    }
    
    // Test 8: For each context, verify get_context_idx returns correct index
    for (idx, &base) in lookup.frame_bases.iter().enumerate() {
        let ctx_idx = lookup.get_context_idx(base);
        assert_eq!(ctx_idx, idx, 
                   "get_context_idx should return correct index for frame base {}", base);
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
    let (lookup, contexts) = build_ncbi_lookup(&queries, 13, true, true, &karlin);
    
    // With single context, all offsets should return 0
    assert_eq!(lookup.get_context_idx(0), 0);
    assert_eq!(lookup.get_context_idx(100), 0);
    assert_eq!(lookup.get_context_idx(-100), 0); // saturating_sub(1) handles negative
    
    // Verify context structure
    assert_eq!(contexts.len(), 1);
    assert_eq!(contexts[0].q_idx, 0);
    assert_eq!(contexts[0].f_idx, 0);
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
    let (lookup, contexts) = build_ncbi_lookup(&queries, 13, true, true, &karlin);
    
    // Should have 12 contexts (6 frames × 2 queries)
    assert_eq!(contexts.len(), 12);
    assert_eq!(lookup.frame_bases.len(), 12);
    
    // Verify query indices are correct
    for (idx, ctx) in contexts.iter().enumerate() {
        let expected_q_idx = (idx / 6) as u32;
        assert_eq!(ctx.q_idx, expected_q_idx, 
                   "Context {} should belong to query {}", idx, expected_q_idx);
    }
    
    // Test: Offset at start of second query (context 6) should return context 6
    if lookup.frame_bases.len() > 6 {
        let query_1_base = lookup.frame_bases[6];
        let ctx_idx = lookup.get_context_idx(query_1_base);
        assert_eq!(ctx_idx, 6, "Offset at start of second query should return context 6");
        assert_eq!(contexts[ctx_idx].q_idx, 1);
        assert_eq!(contexts[ctx_idx].f_idx, 0);
    }
}
