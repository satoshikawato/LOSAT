//! Unit tests for common/evalue.rs with NCBI BLAST reference data
//!
//! **Testing Strategy:**
//!
//! These tests use a multi-layered approach to verify NCBI BLAST compatibility:
//!
//! 1. **Effective Search Space Verification**: The most critical check to catch
//!    length adjustment bugs. We compare LOSAT's effective search space calculation
//!    with NCBI BLAST's reported values (when available).
//!
//! 2. **Parameter Set Testing**: Test multiple scoring parameter combinations,
//!    not just default values, to catch parameter-dependent bugs.
//!
//! 3. **Length Variation Testing**: Test with various sequence lengths to catch
//!    edge cases in length adjustment algorithms.
//!
//! 4. **Integration Test Validation**: These unit tests focus on calculation
//!    correctness. End-to-end NCBI BLAST compatibility is verified in integration
//!    tests (tests/run_comparison.sh).

use LOSAT::algorithm::common::evalue::{
    calculate_evalue_alignment_length, calculate_evalue_database_search,
};
use LOSAT::stats::search_space::SearchSpace;
use super::super::helpers::ncbi_reference::{
    compare_bit_score_with_ncbi, compare_evalue_with_ncbi, get_ncbi_blastn_evalue_cases,
    get_ncbi_tblastx_evalue_cases, verify_effective_search_space,
};

#[test]
#[ignore] // Ignore until reference data is populated
fn test_evalue_against_ncbi_blastn_with_effective_space() {
    // This test verifies both e-value calculation AND effective search space
    // to catch length adjustment bugs
    
    let test_cases = get_ncbi_blastn_evalue_cases();
    
    for case in test_cases {
        // Calculate e-value using LOSAT
        let (bit_score, e_value) = calculate_evalue_database_search(
            case.score,
            case.q_len,
            case.db_len,
            case.db_num_seqs,
            &case.params,
        );
        
        // **Critical**: Verify effective search space calculation
        // This catches the most common bug: incorrect length adjustment
        let search_space = SearchSpace::for_database_search(
            case.q_len,
            case.db_len,
            case.db_num_seqs,
            &case.params,
            true, // use_length_adjustment
        );
        
        // Verify effective search space with small tolerance (0.1%)
        // LOSAT's implementation is a direct port of NCBI BLAST's code,
        // so we expect very close agreement (within floating-point precision)
        // Larger differences (>1%) indicate a bug that needs to be fixed
        if let Err(msg) = verify_effective_search_space(
            search_space.effective_space,
            search_space.length_adjustment,
            case.expected_effective_space,
            case.expected_length_adjustment,
            0.001, // 0.1% tolerance (very strict, since implementation should match exactly)
        ) {
            panic!("Effective search space mismatch for {}: {}", case.test_name, msg);
        }
        
        // Compare bit scores
        if case.expected_bit_score > 0.0 {
            assert!(
                compare_bit_score_with_ncbi(bit_score, case.expected_bit_score, 0.01),
                "Bit score mismatch for {}: {} vs {}",
                case.test_name,
                bit_score,
                case.expected_bit_score
            );
        }
        
        // Compare e-values
        if case.expected_evalue > 0.0 {
            assert!(
                compare_evalue_with_ncbi(e_value, case.expected_evalue, case.tolerance),
                "E-value mismatch for {}: {} vs {} (tolerance: {})",
                case.test_name,
                e_value,
                case.expected_evalue,
                case.tolerance
            );
        }
    }
}

#[test]
#[ignore] // Ignore until reference data is populated
fn test_evalue_against_ncbi_tblastx_with_effective_space() {
    // This test verifies TBLASTX e-value calculations
    // **Important**: Reference data should be from NCBI BLAST runs with
    // -comp_based_stats 0 to avoid composition-based corrections
    
    let test_cases = get_ncbi_tblastx_evalue_cases();
    
    for case in test_cases {
        let (bit_score, e_value) = calculate_evalue_alignment_length(
            case.score,
            case.aln_len, // Used as alignment length
            &case.params,
        );
        
        // Compare with NCBI BLAST expected values
        if case.expected_bit_score > 0.0 {
            assert!(
                compare_bit_score_with_ncbi(bit_score, case.expected_bit_score, 0.01),
                "Bit score mismatch for {}: {} vs {}",
                case.test_name,
                bit_score,
                case.expected_bit_score
            );
        }
        
        if case.expected_evalue > 0.0 {
            assert!(
                compare_evalue_with_ncbi(e_value, case.expected_evalue, case.tolerance),
                "E-value mismatch for {}: {} vs {} (tolerance: {})",
                case.test_name,
                e_value,
                case.expected_evalue,
                case.tolerance
            );
        }
    }
}

#[test]
fn test_effective_search_space_variation() {
    // Test that effective search space changes correctly with sequence length
    // This catches bugs where length adjustment is not applied correctly
    //
    // **Important**: We test the *direction* of change and *sign* of values,
    // not exact values, because NCBI BLAST's effective length calculation
    // is extremely complex and may differ slightly in implementation.
    
    let params = LOSAT::stats::tables::KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    
    // Test with different sequence lengths
    let test_lengths = vec![
        (100, 1000, 1),   // Short sequences
        (1000, 10000, 5), // Medium sequences
        (10000, 100000, 10), // Long sequences
    ];
    
    for (q_len, db_len, db_num_seqs) in test_lengths {
        let ss1 = SearchSpace::for_database_search(q_len, db_len, db_num_seqs, &params, false);
        let ss2 = SearchSpace::for_database_search(q_len, db_len, db_num_seqs, &params, true);
        
        // With length adjustment, effective space should be smaller
        assert!(
            ss2.effective_space < ss1.effective_space,
            "Length adjustment should reduce effective search space"
        );
        
        // But still positive
        assert!(ss2.effective_space > 0.0);
        assert!(ss2.effective_query_len > 0.0);
        assert!(ss2.effective_db_len > 0.0);
        
        // Length adjustment should be reasonable (not larger than sequence length)
        assert!(ss2.length_adjustment >= 0);
        assert!(ss2.length_adjustment < q_len.min(db_len) as i64);
    }
}

#[test]
fn test_evalue_with_fixed_alignment() {
    // Test E-value calculation with a FIXED alignment (not dynamically computed)
    // This separates "alignment extension logic" from "score/E-value calculation logic"
    //
    // **Important**: This test uses hardcoded alignment coordinates to avoid
    // X-drop butterfly effects where slight differences in alignment termination
    // cause score differences.
    
    let params = LOSAT::stats::tables::KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    
    // Fixed alignment: q_start=100, q_end=200, s_start=100, s_end=200
    // This represents a perfect 100bp match
    let fixed_alignment_length = 100;
    let fixed_raw_score = 100; // reward=1, perfect match
    
    let q_len = 1000;
    let db_len = 10000;
    let db_num_seqs = 5;
    
    // Calculate e-value for this fixed alignment
    let (bit_score, e_value) = calculate_evalue_database_search(
        fixed_raw_score,
        q_len,
        db_len,
        db_num_seqs,
        &params,
    );
    
    // Verify that calculations are consistent (same inputs = same outputs)
    let (bit_score2, e_value2) = calculate_evalue_database_search(
        fixed_raw_score,
        q_len,
        db_len,
        db_num_seqs,
        &params,
    );
    
    assert_eq!(bit_score, bit_score2);
    assert_eq!(e_value, e_value2);
    
    // Verify that bit score is positive for positive raw score
    assert!(bit_score > 0.0);
    
    // Verify that e-value is positive
    assert!(e_value > 0.0);
}

#[test]
fn test_parameter_set_independence() {
    // Test that different parameter sets produce different but consistent results
    // This catches bugs where parameters are hardcoded or incorrectly looked up
    
    let params1 = LOSAT::stats::tables::KarlinParams {
        lambda: 1.28, // Megablast
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    
    let params2 = LOSAT::stats::tables::KarlinParams {
        lambda: 0.267, // BLOSUM62
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    
    let score = 100;
    let q_len = 1000;
    let db_len = 10000;
    let db_num_seqs = 5;
    
    let (bs1, ev1) = calculate_evalue_database_search(score, q_len, db_len, db_num_seqs, &params1);
    let (bs2, ev2) = calculate_evalue_database_search(score, q_len, db_len, db_num_seqs, &params2);
    
    // Different parameters should produce different results
    // (This is a sanity check that parameters are actually being used)
    assert_ne!(bs1, bs2, "Different parameters should produce different bit scores");
    assert_ne!(ev1, ev2, "Different parameters should produce different e-values");
}
