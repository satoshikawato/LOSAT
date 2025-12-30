//! Manual test to verify LOSAT's effective search space calculation
//!
//! This test manually calculates effective search space using LOSAT's implementation
//! and compares it with expected values. Since NCBI BLAST doesn't always output
//! effective search space in the header, we calculate it directly using LOSAT.

use LOSAT::stats::search_space::SearchSpace;
use LOSAT::stats::tables::KarlinParams;

/// Test effective search space calculation for BLASTN
/// 
/// Parameters: reward=1, penalty=-2, gapopen=0, gapextend=0
/// These correspond to BLASTN default parameters
#[test]
fn test_effective_search_space_blastn_manual() {
    // BLASTN parameters for reward=1, penalty=-2, gapopen=0, gapextend=0
    // These values are from NCBI BLAST's precomputed tables
    let params = KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    
    // Test case: query length 6188, subject length 6188, 1 sequence
    let query_len = 6188;
    let db_len = 6188;
    let db_num_seqs = 1;
    
    // Calculate using LOSAT
    let search_space = SearchSpace::for_database_search(
        query_len,
        db_len,
        db_num_seqs,
        &params,
        true, // use_length_adjustment
    );
    
    println!("Query length: {}", query_len);
    println!("Database length: {}", db_len);
    println!("Number of sequences: {}", db_num_seqs);
    println!("Length adjustment: {}", search_space.length_adjustment);
    println!("Effective query length: {:.2}", search_space.effective_query_len);
    println!("Effective database length: {:.2}", search_space.effective_db_len);
    println!("Effective search space: {:.2}", search_space.effective_space);
    
    // Verify that length adjustment is reasonable
    assert!(search_space.length_adjustment >= 0);
    assert!(search_space.length_adjustment < query_len as i64);
    
    // Verify that effective lengths are positive
    assert!(search_space.effective_query_len > 0.0);
    assert!(search_space.effective_db_len > 0.0);
    assert!(search_space.effective_space > 0.0);
    
    // Verify that effective lengths are less than original lengths
    assert!(search_space.effective_query_len < query_len as f64);
    assert!(search_space.effective_db_len < db_len as f64);
    
    // For single sequence comparison, effective lengths should be:
    // effective_query_len = query_len - length_adjustment
    // effective_db_len = db_len - length_adjustment
    let expected_effective_query = (query_len as f64) - (search_space.length_adjustment as f64);
    let expected_effective_db = (db_len as f64) - (search_space.length_adjustment as f64);
    
    // Allow small floating-point differences
    assert!((search_space.effective_query_len - expected_effective_query).abs() < 0.01);
    assert!((search_space.effective_db_len - expected_effective_db).abs() < 0.01);
    
    // Effective search space should be the product
    let expected_space = expected_effective_query * expected_effective_db;
    assert!((search_space.effective_space - expected_space).abs() < 0.01);
}

/// Test with different sequence lengths
#[test]
fn test_effective_search_space_various_lengths() {
    let params = KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    
    let test_cases = vec![
        (100, 1000, 1),
        (1000, 10000, 1),
        (1000, 10000, 10),
        (10000, 100000, 100),
    ];
    
    for (q_len, db_len, db_num_seqs) in test_cases {
        let search_space = SearchSpace::for_database_search(
            q_len,
            db_len,
            db_num_seqs,
            &params,
            true,
        );
        
        println!("\nTest case: q_len={}, db_len={}, db_num_seqs={}", q_len, db_len, db_num_seqs);
        println!("  Length adjustment: {}", search_space.length_adjustment);
        println!("  Effective space: {:.2}", search_space.effective_space);
        
        // Basic sanity checks
        assert!(search_space.length_adjustment >= 0);
        assert!(search_space.effective_query_len > 0.0);
        assert!(search_space.effective_db_len > 0.0);
        assert!(search_space.effective_space > 0.0);
        
        // Effective space should be less than simple product
        let simple_space = (q_len as f64) * (db_len as f64);
        assert!(search_space.effective_space < simple_space);
    }
}


