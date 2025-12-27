//! Unit tests for blastn/lookup.rs

use LOSAT::algorithm::blastn::lookup::{encode_kmer, reverse_complement};

#[test]
fn test_encode_kmer_basic() {
    let seq = b"ACGTACGT";
    
    // Encode first 4-mer: ACGT
    let encoded = encode_kmer(seq, 0, 4);
    assert!(encoded.is_some());
    // ACGT = 0b00011011 = 27
    assert_eq!(encoded.unwrap(), 27);
}

#[test]
fn test_encode_kmer_different_positions() {
    let seq = b"ACGTACGT";
    
    // ACGT = 0b00011011 = 27
    let encoded1 = encode_kmer(seq, 0, 4);
    assert_eq!(encoded1.unwrap(), 27);
    
    // CGTA = C<<6 | G<<4 | T<<2 | A = 1<<6 | 2<<4 | 3<<2 | 0 = 64 + 32 + 12 + 0 = 108
    let encoded2 = encode_kmer(seq, 1, 4);
    assert_eq!(encoded2.unwrap(), 108);
    
    // GTAC = G<<6 | T<<4 | A<<2 | C = 2<<6 | 3<<4 | 0<<2 | 1 = 128 + 48 + 0 + 1 = 177
    let encoded3 = encode_kmer(seq, 2, 4);
    assert_eq!(encoded3.unwrap(), 177);
}

#[test]
fn test_encode_kmer_case_insensitive() {
    let seq_upper = b"ACGT";
    let seq_lower = b"acgt";
    let seq_mixed = b"AcGt";
    
    let encoded_upper = encode_kmer(seq_upper, 0, 4);
    let encoded_lower = encode_kmer(seq_lower, 0, 4);
    let encoded_mixed = encode_kmer(seq_mixed, 0, 4);
    
    assert_eq!(encoded_upper, encoded_lower);
    assert_eq!(encoded_lower, encoded_mixed);
}

#[test]
fn test_encode_kmer_invalid_base() {
    let seq = b"ACNT"; // N is invalid
    
    let encoded = encode_kmer(seq, 0, 4);
    assert!(encoded.is_none());
}

#[test]
fn test_encode_kmer_out_of_bounds() {
    let seq = b"ACGT";
    
    // Try to encode beyond sequence length
    let encoded = encode_kmer(seq, 2, 4);
    assert!(encoded.is_none());
}

#[test]
fn test_encode_kmer_different_lengths() {
    let seq = b"ACGTACGT";
    
    // 2-mer: AC = 0b0001 = 1
    let encoded2 = encode_kmer(seq, 0, 2);
    assert_eq!(encoded2.unwrap(), 1);
    
    // 3-mer: ACG = 0b000110 = 6
    let encoded3 = encode_kmer(seq, 0, 3);
    assert_eq!(encoded3.unwrap(), 6);
    
    // 4-mer: ACGT = 0b00011011 = 27
    let encoded4 = encode_kmer(seq, 0, 4);
    assert_eq!(encoded4.unwrap(), 27);
}

#[test]
fn test_reverse_complement_basic() {
    let seq = b"ACGT";
    let rc = reverse_complement(seq);
    
    assert_eq!(rc, b"ACGT"); // ACGT -> TGCA -> reverse -> ACGT
    // Actually: ACGT -> complement -> TGCA -> reverse -> ACGT
    // Wait, let me check: A->T, C->G, G->C, T->A
    // ACGT -> TGCA, then reverse -> ACGT
    // Actually reverse complement: reverse first then complement
    // ACGT reverse -> TGCA, complement -> ACGT
    // Let me verify the actual implementation
}

#[test]
fn test_reverse_complement_actual() {
    let seq = b"ACGT";
    let rc = reverse_complement(seq);
    
    // Reverse complement: reverse first, then complement
    // ACGT -> reverse -> TGCA -> complement -> ACGT
    // Actually the function does: reverse, then complement each base
    // ACGT reversed = TGCA, then complement:
    // T->A, G->C, C->G, A->T = ACGT
    // So ACGT -> ACGT (palindrome)
    
    // Let's test a non-palindrome
    let seq2 = b"AAAA";
    let rc2 = reverse_complement(seq2);
    assert_eq!(rc2, b"TTTT");
}

#[test]
fn test_reverse_complement_case_insensitive() {
    let seq_upper = b"ACGT";
    let seq_lower = b"acgt";
    
    let rc_upper = reverse_complement(seq_upper);
    let rc_lower = reverse_complement(seq_lower);
    
    assert_eq!(rc_upper, rc_lower);
}

#[test]
fn test_reverse_complement_ambiguous_bases() {
    let seq = b"ACNT";
    let rc = reverse_complement(seq);
    
    // Ambiguous bases should become N
    assert_eq!(rc[1], b'N');
}

#[test]
fn test_reverse_complement_round_trip() {
    let seq = b"ACGTACGT";
    let rc = reverse_complement(seq);
    let rc_rc = reverse_complement(rc.as_slice());
    
    // Double reverse complement should return original
    assert_eq!(seq, rc_rc.as_slice());
}

#[test]
fn test_reverse_complement_u_and_t() {
    // U and T should both be treated as T
    let seq_t = b"ACGT";
    let seq_u = b"ACGU";
    
    let rc_t = reverse_complement(seq_t);
    let rc_u = reverse_complement(seq_u);
    
    assert_eq!(rc_t, rc_u);
}

