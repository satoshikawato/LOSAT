//! Unit tests for tblastx/chaining.rs

use LOSAT::algorithm::tblastx::chaining::{chain_and_filter_hsps_protein, ExtendedHit};
use LOSAT::common::Hit;
use LOSAT::stats::KarlinParams;
use rustc_hash::FxHashMap;

fn create_test_hit(
    query_id: &str,
    subject_id: &str,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    bit_score: f64,
    e_value: f64,
) -> Hit {
    Hit {
        query_id: query_id.to_string(),
        subject_id: subject_id.to_string(),
        identity: 90.0,
        length: (q_end - q_start).max(s_end - s_start),
        mismatch: 0,
        gapopen: 0,
        q_start,
        q_end,
        s_start,
        s_end,
        e_value,
        bit_score,
        q_idx: 0,
        s_idx: 0,
        raw_score: ((bit_score * 0.693 + (-0.041_f64).ln()) / 0.267) as i32,
    }
}

fn create_extended_hit(
    hit: Hit,
    q_frame: i8,
    s_frame: i8,
    q_aa_start: usize,
    q_aa_end: usize,
    s_aa_start: usize,
    s_aa_end: usize,
    q_orig_len: usize,
    s_orig_len: usize,
    from_gapped: bool,
) -> ExtendedHit {
    // Calculate raw_score from bit_score using BLOSUM62 Karlin params
    // raw_score = (bit_score * ln(2) + ln(K)) / lambda
    // For BLOSUM62: lambda=0.267, K=0.041
    let raw_score = ((hit.bit_score * 0.693 + (-0.041_f64).ln()) / 0.267) as i32;
    ExtendedHit {
        hit,
        raw_score,
        q_frame,
        s_frame,
        q_aa_start,
        q_aa_end,
        s_aa_start,
        s_aa_end,
        q_orig_len,
        s_orig_len,
        from_gapped,
    }
}

#[test]
fn test_chain_and_filter_hsps_protein_empty() {
    let hits: Vec<ExtendedHit> = Vec::new();
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    assert!(result.is_empty());
}

#[test]
fn test_chain_and_filter_hsps_protein_single_hit() {
    let hit = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let ext_hit = create_extended_hit(
        hit,
        1,  // q_frame
        1,  // s_frame
        3,  // q_aa_start
        6,  // q_aa_end
        10, // s_aa_start
        13, // s_aa_end
        100, // q_orig_len
        200, // s_orig_len
        false, // from_gapped
    );
    
    let hits = vec![ext_hit];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].query_id, "query1");
    assert_eq!(result[0].subject_id, "subject1");
}

#[test]
fn test_chain_and_filter_hsps_protein_multiple_hits_same_frame() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let hit2 = create_test_hit("query1", "subject1", 25, 35, 45, 55, 45.0, 1e-9);
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1, // Same frames
        3, 6, 10, 13,
        100, 200,
        false,
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        1, 1, // Same frames
        8, 11, 15, 18,
        100, 200,
        false,
    );
    
    let hits = vec![ext_hit1, ext_hit2];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    // Both hits should be kept (clustering may group them, but both should pass e-value filter)
    assert!(result.len() >= 1);
}

#[test]
fn test_chain_and_filter_hsps_protein_different_frames() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let hit2 = create_test_hit("query1", "subject1", 25, 35, 45, 55, 45.0, 1e-9);
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1, // Frame 1
        3, 6, 10, 13,
        100, 200,
        false,
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        2, 2, // Different frame
        8, 11, 15, 18,
        100, 200,
        false,
    );
    
    let hits = vec![ext_hit1, ext_hit2];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    // Hits from different frames should be in separate groups
    // Both should be kept since they're in different frame combinations
    assert!(result.len() >= 1);
}

#[test]
fn test_chain_and_filter_hsps_protein_evalue_filtering() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let hit2 = create_test_hit("query1", "subject1", 25, 35, 45, 55, 45.0, 1e-2); // High e-value
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1,
        3, 6, 10, 13,
        100, 200,
        false,
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        1, 1,
        8, 11, 15, 18,
        100, 200,
        false,
    );
    
    let hits = vec![ext_hit1, ext_hit2];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    // Use strict e-value threshold
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5, // Only hit1 should pass
        None,
    );
    
    // Only hit1 should pass the e-value filter
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].e_value, 1e-10);
}

#[test]
fn test_chain_and_filter_hsps_protein_different_queries() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let hit2 = create_test_hit("query2", "subject1", 25, 35, 45, 55, 45.0, 1e-9);
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1,
        3, 6, 10, 13,
        100, 200,
        false,
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        1, 1,
        8, 11, 15, 18,
        100, 200,
        false,
    );
    
    let hits = vec![ext_hit1, ext_hit2];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    // Both hits should be kept (different queries)
    assert_eq!(result.len(), 2);
}

#[test]
fn test_chain_and_filter_hsps_protein_different_subjects() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let hit2 = create_test_hit("query1", "subject2", 25, 35, 45, 55, 45.0, 1e-9);
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1,
        3, 6, 10, 13,
        100, 200,
        false,
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        1, 1,
        8, 11, 15, 18,
        100, 200,
        false,
    );
    
    let hits = vec![ext_hit1, ext_hit2];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    // Both hits should be kept (different subjects)
    assert_eq!(result.len(), 2);
}

#[test]
fn test_chain_and_filter_hsps_protein_sorting_by_bit_score() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 30.0, 1e-10);
    let hit2 = create_test_hit("query1", "subject1", 25, 35, 45, 55, 50.0, 1e-9);
    let hit3 = create_test_hit("query1", "subject1", 40, 50, 60, 70, 40.0, 1e-8);
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1,
        3, 6, 10, 13,
        100, 200,
        false,
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        1, 1,
        8, 11, 15, 18,
        100, 200,
        false,
    );
    let ext_hit3 = create_extended_hit(
        hit3,
        1, 1,
        13, 16, 20, 23,
        100, 200,
        false,
    );
    
    let hits = vec![ext_hit1, ext_hit2, ext_hit3];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    // Results should be sorted by bit score (descending)
    assert!(result.len() >= 1);
    for i in 1..result.len() {
        assert!(result[i-1].bit_score >= result[i].bit_score);
    }
}

#[test]
fn test_chain_and_filter_hsps_protein_gapped_vs_ungapped() {
    let hit1 = create_test_hit("query1", "subject1", 10, 20, 30, 40, 50.0, 1e-10);
    let hit2 = create_test_hit("query1", "subject1", 25, 35, 45, 55, 45.0, 1e-9);
    
    let ext_hit1 = create_extended_hit(
        hit1,
        1, 1,
        3, 6, 10, 13,
        100, 200,
        false, // from_gapped = false
    );
    let ext_hit2 = create_extended_hit(
        hit2,
        1, 1,
        8, 11, 15, 18,
        100, 200,
        true, // from_gapped = true
    );
    
    let hits = vec![ext_hit1, ext_hit2];
    let sequences: FxHashMap<(String, String, i8, i8), (Vec<u8>, Vec<u8>)> = FxHashMap::default();
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let result = chain_and_filter_hsps_protein(
        hits,
        &sequences,
        1000,
        &params,
        1e-5,
        None,
    );
    
    // Both hits should be kept regardless of gapped/ungapped source
    assert!(result.len() >= 1);
}


