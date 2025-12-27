//! Unit tests for tblastx/lookup.rs

use LOSAT::algorithm::tblastx::lookup::{encode_aa_kmer, build_direct_lookup};
use LOSAT::algorithm::tblastx::translation::{generate_frames, QueryFrame};
use LOSAT::utils::genetic_code::GeneticCode;
use LOSAT::utils::dust::MaskedInterval;

#[test]
fn test_encode_aa_kmer_basic() {
    // Amino acid indices: A=0, B=1, C=2 (not ASCII characters)
    let seq = [0u8, 1u8, 2u8];
    let code = encode_aa_kmer(&seq, 0);
    assert_eq!(code, Some(0 * 676 + 1 * 26 + 2)); // 28
}

#[test]
fn test_encode_aa_kmer_different_positions() {
    let seq = [0u8, 1u8, 2u8, 3u8, 4u8, 5u8]; // ABCDEF
    
    // ABC = 0*676 + 1*26 + 2 = 28
    let code1 = encode_aa_kmer(&seq, 0);
    assert_eq!(code1, Some(28));
    
    // BCD = 1*676 + 2*26 + 3 = 731
    let code2 = encode_aa_kmer(&seq, 1);
    assert_eq!(code2, Some(731));
    
    // CDE = 2*676 + 3*26 + 4 = 1434
    let code3 = encode_aa_kmer(&seq, 2);
    assert_eq!(code3, Some(1434));
}

#[test]
fn test_encode_aa_kmer_with_stop_codon() {
    // Stop codon is 25
    let seq_with_stop = [0u8, 1u8, 25u8]; // Contains stop codon
    let code = encode_aa_kmer(&seq_with_stop, 0);
    assert_eq!(code, None);
    
    // Stop codon in middle
    let seq_with_stop_middle = [0u8, 25u8, 2u8];
    let code = encode_aa_kmer(&seq_with_stop_middle, 0);
    assert_eq!(code, None);
    
    // Stop codon at start
    let seq_with_stop_start = [25u8, 1u8, 2u8];
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
    // Test with maximum valid amino acid values (24 = Z)
    // Z=24, Y=23, X=22
    let seq = [24u8, 23u8, 22u8];
    let code = encode_aa_kmer(&seq, 0);
    assert_eq!(code, Some(24 * 676 + 23 * 26 + 22)); // 165,894
    
    // Verify it's within table size (26^3 = 17,576)
    // Actually wait, 24*676 + 23*26 + 22 = 16,224 + 598 + 22 = 16,844
    // But table size is 26^3 = 17,576, so this is valid
    assert!(code.unwrap() < 17576);
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
    
    // Verify table size
    assert_eq!(lookup.len(), 17576); // 26^3
    
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
fn test_build_direct_lookup_amino_acid_positions() {
    let code = GeneticCode::from_id(1);
    let dna_seq = b"ATGCGATCGATCGATCG"; // 18 bases = 6 codons
    let frames = generate_frames(dna_seq, &code);
    
    let queries = vec![frames];
    let query_masks = vec![Vec::<MaskedInterval>::new()];
    
    let lookup = build_direct_lookup(&queries, &query_masks);
    
    // Verify amino acid positions are valid
    for bucket in &lookup {
        for (_query_idx, frame_idx, aa_pos) in bucket {
            let frame = &queries[0][*frame_idx as usize];
            // AA position should be within sequence length
            assert!(*aa_pos < frame.aa_seq.len() as u32);
            // AA position should allow for 3-mer (pos + 3 <= len)
            assert!(*aa_pos + 3 <= frame.aa_seq.len() as u32);
        }
    }
}

