//! Direct comparison test with NCBI BLAST's effective search space
//!
//! This test compares LOSAT's effective search space calculation
//! with values extracted from actual NCBI BLAST output.

use LOSAT::stats::search_space::SearchSpace;
use LOSAT::stats::tables::KarlinParams;

/// Test case from actual NCBI BLAST run
/// 
/// Query: test_sequences/query.fasta (6188 bp)
/// Subject: test_sequences/subject.fasta (6188 bp)
/// Parameters: reward=1, penalty=-2, gapopen=0, gapextend=0
/// NCBI BLAST reported: "Effective search space used: 38081241"
#[test]
fn test_effective_search_space_ncbi_comparison() {
    // BLASTN parameters for reward=1, penalty=-2, gapopen=0, gapextend=0
    let params = KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    
    // Test case from actual NCBI BLAST run
    let query_len = 6188;
    let db_len = 6188;
    let db_num_seqs = 1;
    let ncbi_effective_space = 38081241.0;
    
    // Calculate using LOSAT
    let search_space = SearchSpace::for_database_search(
        query_len,
        db_len,
        db_num_seqs,
        &params,
        true, // use_length_adjustment
    );
    
    println!("Test case: NCBI BLAST comparison");
    println!("  Query length: {}", query_len);
    println!("  Database length: {}", db_len);
    println!("  Number of sequences: {}", db_num_seqs);
    println!("  Length adjustment: {}", search_space.length_adjustment);
    println!("  LOSAT effective space: {:.2}", search_space.effective_space);
    println!("  NCBI effective space: {:.2}", ncbi_effective_space);
    
    let diff = (search_space.effective_space - ncbi_effective_space).abs();
    let relative_diff = diff / ncbi_effective_space;
    
    println!("  Absolute difference: {:.2}", diff);
    println!("  Relative difference: {:.4}%", relative_diff * 100.0);
    
    // Verify with 0.1% tolerance (very strict, since implementation should match exactly)
    assert!(
        relative_diff <= 0.001,
        "Effective search space mismatch: LOSAT={:.2}, NCBI={:.2}, diff={:.4}% (tolerance: 0.1%)",
        search_space.effective_space,
        ncbi_effective_space,
        relative_diff * 100.0
    );
    
    // Also verify length adjustment is reasonable
    assert!(search_space.length_adjustment >= 0);
    assert!(search_space.length_adjustment < query_len as i64);
}


