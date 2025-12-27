//! Unit tests for tblastx/translation.rs

use LOSAT::algorithm::tblastx::translation::generate_frames;
use LOSAT::utils::genetic_code::GeneticCode;

#[test]
fn test_generate_frames_basic() {
    let code = GeneticCode::from_id(1); // Standard genetic code
    let seq = b"ATGCGATCGATCGATCG"; // 18 bases, divisible by 3
    
    let frames = generate_frames(seq, &code);
    
    // Should generate up to 6 frames (3 forward + 3 reverse)
    assert!(frames.len() <= 6);
    assert!(frames.len() >= 3); // At least forward frames should be generated
    
    // Check frame numbers
    let frame_numbers: Vec<i8> = frames.iter().map(|f| f.frame).collect();
    for frame_num in frame_numbers {
        assert!(frame_num >= -3 && frame_num <= 3 && frame_num != 0);
    }
}

#[test]
fn test_generate_frames_forward_frames() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATGCGATCGATCGATCG";
    
    let frames = generate_frames(seq, &code);
    
    // Check that forward frames (1, 2, 3) are present
    let forward_frames: Vec<_> = frames.iter().filter(|f| f.frame > 0).collect();
    assert!(!forward_frames.is_empty());
    
    // All forward frames should have frame numbers 1, 2, or 3
    for frame in forward_frames {
        assert!(frame.frame >= 1 && frame.frame <= 3);
    }
}

#[test]
fn test_generate_frames_reverse_frames() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATGCGATCGATCGATCG";
    
    let frames = generate_frames(seq, &code);
    
    // Check that reverse frames (-1, -2, -3) are present
    let reverse_frames: Vec<_> = frames.iter().filter(|f| f.frame < 0).collect();
    assert!(!reverse_frames.is_empty());
    
    // All reverse frames should have frame numbers -1, -2, or -3
    for frame in reverse_frames {
        assert!(frame.frame >= -3 && frame.frame <= -1);
    }
}

#[test]
fn test_generate_frames_original_length() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATGCGATCGATCGATCG";
    let seq_len = seq.len();
    
    let frames = generate_frames(seq, &code);
    
    // All frames should preserve original sequence length
    for frame in &frames {
        assert_eq!(frame.orig_len, seq_len);
    }
}

#[test]
fn test_generate_frames_amino_acid_sequence() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATGCGATCGATCGATCG"; // 18 bases = 6 codons
    
    let frames = generate_frames(seq, &code);
    
    // Each frame should have translated amino acid sequence
    for frame in &frames {
        // Length should be approximately seq_len / 3 (may vary by frame offset)
        assert!(!frame.aa_seq.is_empty());
        assert!(frame.aa_seq.len() <= seq.len() / 3);
    }
}

#[test]
fn test_generate_frames_short_sequence() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATG"; // Only 3 bases, minimum for one codon
    
    let frames = generate_frames(seq, &code);
    
    // Should still generate some frames if sequence is long enough
    // For 3 bases, only frame 1 (forward, offset 0) can be translated
    assert!(frames.len() <= 6);
}

#[test]
fn test_generate_frames_very_short_sequence() {
    let code = GeneticCode::from_id(1);
    let seq = b"AT"; // Only 2 bases, too short for translation
    
    let frames = generate_frames(seq, &code);
    
    // Should handle gracefully, may return empty or minimal frames
    // The exact behavior depends on implementation
    assert!(frames.len() <= 6);
}

#[test]
fn test_generate_frames_different_genetic_codes() {
    // Test with different genetic code tables
    let code1 = GeneticCode::from_id(1); // Standard
    let code11 = GeneticCode::from_id(11); // Bacterial/Archaeal
    
    let seq = b"ATGCGATCGATCGATCG";
    
    let frames1 = generate_frames(seq, &code1);
    let frames11 = generate_frames(seq, &code11);
    
    // Should generate same number of frames regardless of genetic code
    assert_eq!(frames1.len(), frames11.len());
    
    // Frame numbers should be the same
    let frame_nums1: Vec<i8> = frames1.iter().map(|f| f.frame).collect();
    let frame_nums11: Vec<i8> = frames11.iter().map(|f| f.frame).collect();
    assert_eq!(frame_nums1, frame_nums11);
}

