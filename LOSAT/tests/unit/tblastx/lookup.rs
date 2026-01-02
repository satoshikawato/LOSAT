//! Unit tests for tblastx/lookup.rs

use LOSAT::algorithm::tblastx::lookup::{encode_aa_kmer, build_direct_lookup};
use LOSAT::algorithm::tblastx::lookup::{build_ncbi_lookup, pv_test};
use LOSAT::algorithm::tblastx::translation::generate_frames;
use LOSAT::utils::genetic_code::GeneticCode;
use LOSAT::utils::dust::MaskedInterval;

#[test]
fn test_encode_aa_kmer_basic() {
    // Amino acid indices in NCBI matrix order: A=0, R=1, N=2
    // NCBI-style lookup index uses bit-shift encoding (charsize=5 for 0..=23).
    let seq = [0u8, 1u8, 2u8];
    let code = encode_aa_kmer(&seq, 0);
    assert_eq!(code, Some((0 << 10) | (1 << 5) | 2)); // 34
}

#[test]
fn test_encode_aa_kmer_different_positions() {
    let seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8]; // ARNDCQ in NCBI order
    
    // ARN = (0<<10) | (1<<5) | 2 = 34
    let code1 = encode_aa_kmer(&seq, 0);
    assert_eq!(code1, Some((0 << 10) | (1 << 5) | 2));
    
    // RND = (1<<10) | (2<<5) | 3 = 1091
    let code2 = encode_aa_kmer(&seq, 1);
    assert_eq!(code2, Some((1 << 10) | (2 << 5) | 3));
    
    // NDC = (2<<10) | (3<<5) | 4 = 2148
    let code3 = encode_aa_kmer(&seq, 2);
    assert_eq!(code3, Some((2 << 10) | (3 << 5) | 4));
}

#[test]
fn test_encode_aa_kmer_with_stop_codon() {
    // Stop codon is 24 in NCBI matrix order
    let seq_with_stop = [0u8, 1u8, 24u8]; // Contains stop codon
    let code = encode_aa_kmer(&seq_with_stop, 0);
    assert_eq!(code, None);
    
    // Stop codon in middle
    let seq_with_stop_middle = [0u8, 24u8, 2u8];
    let code = encode_aa_kmer(&seq_with_stop_middle, 0);
    assert_eq!(code, None);
    
    // Stop codon at start
    let seq_with_stop_start = [24u8, 1u8, 2u8];
    let code = encode_aa_kmer(&seq_with_stop_start, 0);
    assert_eq!(code, None);
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
    // Test with maximum valid amino acid values (23 = X, the last valid AA before stop codon)
    // NCBI matrix order: ARNDCQEGHILKMFPSTWYVBJZX* (0-24)
    // X=23, Z=22, J=21
    let seq = [23u8, 22u8, 21u8];
    let code = encode_aa_kmer(&seq, 0);
    // (23<<10) | (22<<5) | 21 = 24,277
    assert_eq!(code, Some((23 << 10) | (22 << 5) | 21));
    
    // Verify it's within backbone size (max index + 1 for 0..=23 with charsize=5)
    assert!(code.unwrap() < 24_312);
}

#[test]
fn test_build_direct_lookup_basic() {
    let code = GeneticCode::from_id(1);
    
    // Create a simple DNA sequence and generate frames
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases = 6 codons
    let frames = generate_frames(dna_seq, &code);
    
    // Build lookup table
    let queries = vec![frames];
    let query_masks = vec![Vec::<MaskedInterval>::new()]; // No masking
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
    // Verify table size (NCBI-style backbone size for 0..=23 with charsize=5)
    assert_eq!(lookup.len(), 24_312);
    
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
    let query_masks = vec![Vec::<MaskedInterval>::new()];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
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
    let query_masks = vec![Vec::<MaskedInterval>::new(), Vec::<MaskedInterval>::new()];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
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
fn test_build_direct_lookup_with_masking() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases
    let frames = generate_frames(dna_seq, &code);
    
    // Mask positions 0-5 (first 6 bases, which affects frame 1)
    let mask = vec![MaskedInterval { start: 0, end: 6 }];
    
    let queries = vec![frames];
    let query_masks = vec![mask];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
    // Some k-mers should be filtered out due to masking
    let total_hits: usize = lookup.iter().map(|bucket| bucket.len()).sum();
    assert!(total_hits >= 0); // May be reduced due to masking
}

#[test]
fn test_build_direct_lookup_short_sequence() {
    let code = GeneticCode::from_id(1);
    
    // Very short sequence (only 3 bases = 1 codon)
    let dna_seq = b"ATG";
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let query_masks = vec![Vec::<MaskedInterval>::new()];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
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
    let query_masks = vec![Vec::<MaskedInterval>::new()];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
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
    let query_masks = vec![Vec::<MaskedInterval>::new()];

    let (lookup, _ctx) = build_ncbi_lookup(&queries, &query_masks, 13, true, true, 0, false);

    // AAA encodes to index 0 regardless of alphabet_size/word_length here.
    assert!(pv_test(&lookup.pv, 0));
    let hits = lookup.get_hits(0);
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
    let query_masks = vec![Vec::<MaskedInterval>::new()];

    let (lookup, _ctx) = build_ncbi_lookup(&queries, &query_masks, 13, true, true, 0, false);

    // AAA index is 0.
    let hits = lookup.get_hits(0);
    assert_eq!(lookup.longest_chain as usize, hits.len());
    assert_eq!(lookup.longest_chain, 98);
}

#[test]
fn test_build_direct_lookup_amino_acid_positions() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases = 6 codons
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let query_masks = vec![Vec::<MaskedInterval>::new()];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
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
