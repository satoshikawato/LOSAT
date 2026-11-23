use super::*;

fn default_blosum62_params() -> KarlinParams {
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

fn default_blastn_params() -> KarlinParams {
    KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    }
}

#[test]
fn test_length_adjustment_basic() {
    let params = default_blosum62_params();
    let result = compute_length_adjustment_ncbi(100, 10000, 10, &params);

    // Length adjustment should be positive
    assert!(result.length_adjustment >= 0);
    // Length adjustment should be less than query length
    assert!(result.length_adjustment < 100);
}

#[test]
fn test_length_adjustment_nucleotide() {
    let params = default_blastn_params();
    let result = compute_length_adjustment_ncbi(1000, 100000, 100, &params);

    assert!(result.length_adjustment >= 0);
}

#[test]
fn test_length_adjustment_edge_cases() {
    let params = default_blosum62_params();

    // Very short query
    let result = compute_length_adjustment_ncbi(10, 10000, 10, &params);
    assert!(result.length_adjustment >= 0);
    assert!(result.length_adjustment <= 10);

    // Single sequence database
    let result = compute_length_adjustment_ncbi(100, 100, 1, &params);
    assert!(result.length_adjustment >= 0);
}

#[test]
fn test_length_adjustment_zero_inputs() {
    let params = default_blosum62_params();

    // Zero query length
    let result = compute_length_adjustment_ncbi(0, 10000, 10, &params);
    assert_eq!(result.length_adjustment, 0);

    // Zero database length
    let result = compute_length_adjustment_ncbi(100, 0, 10, &params);
    assert_eq!(result.length_adjustment, 0);
}

#[test]
fn test_length_adjustment_includes_beta() {
    // Test that beta is properly included in the calculation
    let params_with_beta = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0, // Negative beta
    };

    let params_zero_beta = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: 0.0, // Zero beta
    };

    let result_with_beta = compute_length_adjustment_ncbi(100, 10000, 10, &params_with_beta);
    let result_zero_beta = compute_length_adjustment_ncbi(100, 10000, 10, &params_zero_beta);

    // With negative beta, the length adjustment should be smaller
    assert!(result_with_beta.length_adjustment <= result_zero_beta.length_adjustment);
}

#[test]
fn test_length_adjustment_db_num_seqs_effect() {
    let params = default_blosum62_params();

    // Same total database length, different number of sequences
    let result_few_seqs = compute_length_adjustment_ncbi(100, 10000, 10, &params);
    let result_many_seqs = compute_length_adjustment_ncbi(100, 10000, 100, &params);

    // More sequences should generally result in different (often smaller) adjustment
    // because the constraint `K*(m-A)*(n-N*A) > max(m, n)` is tighter
    assert!(result_many_seqs.length_adjustment <= result_few_seqs.length_adjustment);
}
