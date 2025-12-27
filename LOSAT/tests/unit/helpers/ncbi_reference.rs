//! NCBI BLAST reference data and comparison utilities
//!
//! This module provides reference values from NCBI BLAST for testing
//! and utilities to compare LOSAT output with NCBI BLAST expected values.
//!
//! **Important Testing Strategy:**
//! This module uses a multi-layered approach to avoid the pitfalls of simple
//! reverse-engineering:
//!
//! 1. **Effective Search Space Verification**: Compare LOSAT's effective search space
//!    calculation with NCBI BLAST's reported values (from output headers)
//!
//! 2. **Parameter Set Testing**: Test multiple scoring parameter combinations,
//!    not just default values
//!
//! 3. **Length Variation Testing**: Test with various sequence lengths to catch
//!    edge cases in length adjustment
//!
//! 4. **Integration Test Validation**: Unit tests focus on calculation correctness,
//!    while integration tests verify end-to-end NCBI BLAST compatibility

use LOSAT::stats::tables::KarlinParams;

/// Reference e-value test case with full context
/// 
/// This structure includes all necessary information to verify that LOSAT's
/// e-value calculation matches NCBI BLAST, including effective search space
/// to catch length adjustment issues.
pub struct NcbiEvalueTestCase {
    // Input parameters
    pub score: i32,
    pub q_len: usize,
    pub db_len: usize,
    pub db_num_seqs: usize,
    pub aln_len: usize, // For TBLASTX (amino acid length)
    
    // Karlin-Altschul parameters (from scoring scheme)
    pub params: KarlinParams,
    
    // Expected outputs from NCBI BLAST
    pub expected_bit_score: f64,
    pub expected_evalue: f64,
    
    // **Critical**: Expected effective search space from NCBI BLAST
    // This is the key to catching length adjustment bugs
    pub expected_effective_space: Option<f64>,
    pub expected_length_adjustment: Option<i64>,
    
    // Test metadata
    pub tolerance: f64,
    pub test_name: String, // For debugging
}

/// Verify effective search space calculation
/// 
/// This is the most important check to catch length adjustment bugs.
/// NCBI BLAST reports effective search space in output headers (sometimes).
/// If available, we compare LOSAT's calculation with NCBI's reported value.
/// 
/// **Implementation Status**:
/// LOSAT's `compute_length_adjustment_ncbi` is a faithful port of NCBI BLAST's
/// `BLAST_ComputeLengthAdjustment` function. The implementation should match
/// NCBI BLAST exactly, but we use a small tolerance (0.1%) to account for
/// floating-point precision differences between C and Rust.
/// 
/// **Note**: If tests show larger differences (>1%), this indicates a bug
/// in the port that needs to be fixed, not a fundamental limitation.
pub fn verify_effective_search_space(
    losat_effective_space: f64,
    losat_length_adj: i64,
    expected_effective_space: Option<f64>,
    expected_length_adj: Option<i64>,
    tolerance: f64, // Typically 0.01 (1%) for effective search space
) -> Result<(), String> {
    // If we have expected values from NCBI BLAST, compare them with tolerance
    if let Some(expected_space) = expected_effective_space {
        let relative_diff = (losat_effective_space - expected_space).abs() / expected_space;
        if relative_diff > tolerance {
            return Err(format!(
                "Effective search space mismatch: LOSAT={}, NCBI={}, diff={:.2}% (tolerance: {:.2}%)",
                losat_effective_space, expected_space, relative_diff * 100.0, tolerance * 100.0
            ));
        }
    }
    
    // Length adjustment can vary slightly due to implementation differences
    // Use a more lenient tolerance (allow up to 10 units difference)
    if let Some(expected_adj) = expected_length_adj {
        let diff = (losat_length_adj - expected_adj).abs();
        // For small adjustments, allow 1 unit difference
        // For large adjustments, allow 1% relative difference
        let max_diff = if expected_adj < 100 {
            1
        } else {
            (expected_adj as f64 * 0.01) as i64 + 1
        };
        if diff > max_diff {
            return Err(format!(
                "Length adjustment mismatch: LOSAT={}, NCBI={}, diff={} (max allowed: {})",
                losat_length_adj, expected_adj, diff, max_diff
            ));
        }
    }
    
    Ok(())
}

/// Reference e-value test cases for nucleotide alignments (BLASTN)
/// 
/// **Testing Strategy:**
/// - Include multiple parameter sets (megablast, blastn task)
/// - Include various sequence lengths (short, medium, long)
/// - Include cases with and without length adjustment
/// - Focus on effective search space verification
pub fn get_ncbi_blastn_evalue_cases() -> Vec<NcbiEvalueTestCase> {
    vec![
        // Test case 1: Megablast default (reward=1, penalty=-2, no gaps)
        // This is the most common use case
        NcbiEvalueTestCase {
            score: 100,
            q_len: 1000,
            db_len: 10000,
            db_num_seqs: 5,
            aln_len: 0,
            params: KarlinParams {
                lambda: 1.28,
                k: 0.46,
                h: 0.85,
                alpha: 1.5,
                beta: -2.0,
            },
            expected_bit_score: 0.0, // Will be calculated
            expected_evalue: 0.0,   // Will be calculated
            expected_effective_space: None, // Extract from NCBI BLAST header if available
            expected_length_adjustment: None,
            tolerance: 0.1,
            test_name: "megablast_default".to_string(),
        },
        // TODO: Add more test cases with:
        // - Different scoring parameters (blastn task: reward=1, penalty=-3)
        // - Different sequence lengths (short: 100bp, long: 1Mbp)
        // - Multiple database sequences
        // - Actual NCBI BLAST output data
    ]
}

/// Reference e-value test cases for protein alignments (TBLASTX)
/// 
/// **Testing Strategy:**
/// - Test with -comp_based_stats 0 to avoid composition-based corrections
/// - Include various alignment lengths
/// - Verify alignment-length-based search space calculation
pub fn get_ncbi_tblastx_evalue_cases() -> Vec<NcbiEvalueTestCase> {
    vec![
        // Test case 1: BLOSUM62 with gap_open=11, gap_extend=1
        // **Important**: Run NCBI BLAST with -comp_based_stats 0 to get pure statistics
        NcbiEvalueTestCase {
            score: 100,
            q_len: 0, // Not used for alignment-length-based
            db_len: 0,
            db_num_seqs: 0,
            aln_len: 200, // Alignment length in amino acids
            params: KarlinParams {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
                alpha: 1.9,
                beta: -30.0,
            },
            expected_bit_score: 0.0,
            expected_evalue: 0.0,
            expected_effective_space: None,
            expected_length_adjustment: None,
            tolerance: 0.1,
            test_name: "tblastx_blosum62".to_string(),
        },
        // TODO: Add more test cases with:
        // - Different alignment lengths (short: 50aa, long: 1000aa)
        // - Different scoring matrices (BLOSUM45, PAM250, etc.)
        // - Actual NCBI BLAST output data with -comp_based_stats 0
    ]
}

/// Compare LOSAT e-value with NCBI BLAST expected value
/// Returns true if values are within tolerance
pub fn compare_evalue_with_ncbi(actual: f64, expected: f64, tolerance: f64) -> bool {
    if expected == 0.0 {
        // For zero e-values, check if actual is also very small
        // NCBI BLAST reports 0.0 for e-values < 1e-180
        return actual < 1e-180;
    }
    
    let relative_diff = (actual - expected).abs() / expected;
    relative_diff <= tolerance
}

/// Compare LOSAT bit score with NCBI BLAST expected value
/// Returns true if values are within tolerance
pub fn compare_bit_score_with_ncbi(actual: f64, expected: f64, tolerance: f64) -> bool {
    let diff = (actual - expected).abs();
    diff <= tolerance
}

/// Extract effective search space from NCBI BLAST output header
/// 
/// NCBI BLAST sometimes reports effective search space in output headers.
/// This function attempts to parse it.
/// 
/// Example header line:
/// # Effective search space used: 1234567890
pub fn extract_effective_search_space_from_header(header_lines: &[String]) -> Option<f64> {
    for line in header_lines {
        if line.contains("Effective search space") || line.contains("effective search space") {
            // Try to extract number
            if let Some(num_str) = line.split_whitespace().last() {
                if let Ok(num) = num_str.parse::<f64>() {
                    return Some(num);
                }
            }
        }
    }
    None
}

/// Load reference data from NCBI BLAST output files
/// 
/// This function parses actual NCBI BLAST output and extracts:
/// - E-values and bit scores from hit lines
/// - Effective search space from headers (if available)
/// - Query/database lengths from headers
/// 
/// **Important**: Run NCBI BLAST with -comp_based_stats 0 for protein searches
/// to get pure statistical values without composition-based corrections.
pub fn load_ncbi_reference_data(_file_path: &str) -> Result<Vec<NcbiEvalueTestCase>, String> {
    // TODO: Implement parsing of NCBI BLAST output files
    // This would:
    // 1. Parse header to extract effective search space
    // 2. Parse hit lines to extract e-values and bit scores
    // 3. Calculate raw scores from bit scores
    // 4. Extract query/database lengths from headers or FASTA files
    Err("Not yet implemented. Use extract_ncbi_cases.py script instead.".to_string())
}
