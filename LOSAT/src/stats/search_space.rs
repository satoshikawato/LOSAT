use super::length_adjustment::{compute_length_adjustment_ncbi, compute_length_adjustment_simple};
use super::tables::KarlinParams;

/// Effective search space with length-adjusted values
#[derive(Debug, Clone, Copy)]
pub struct SearchSpace {
    /// Effective query length after adjustment
    pub effective_query_len: f64,
    /// Effective database length after adjustment
    pub effective_db_len: f64,
    /// Effective search space (product of effective lengths)
    pub effective_space: f64,
    /// The length adjustment value used (for debugging/verification)
    pub length_adjustment: i64,
}

impl SearchSpace {
    /// Create a simple search space without length adjustment (BLEMIR default)
    ///
    /// This is the original BLEMIR behavior for backward compatibility.
    /// For NCBI-compatible E-values, use `with_length_adjustment` or `for_database_search`.
    pub fn simple(query_len: usize, db_len: usize) -> Self {
        let q = query_len as f64;
        let d = db_len as f64;
        Self {
            effective_query_len: q,
            effective_db_len: d,
            effective_space: q * d,
            length_adjustment: 0,
        }
    }

    /// Create search space with NCBI-style length adjustment for single-sequence comparison.
    ///
    /// This uses the NCBI BLAST algorithm for computing length adjustment,
    /// which accounts for the fact that alignments cannot extend beyond
    /// sequence boundaries.
    ///
    /// The algorithm solves for the fixed point:
    ///   ell = alpha/lambda * (logK + log(ss)) + beta
    /// where ss = (m - ell) * (n - ell)
    ///
    /// For database searches with multiple sequences, use `for_database_search` instead.
    pub fn with_length_adjustment(query_len: usize, db_len: usize, params: &KarlinParams) -> Self {
        let m = query_len as f64;
        let n = db_len as f64;

        // Use NCBI-compatible length adjustment (single sequence = db_num_seqs = 1)
        let result = compute_length_adjustment_simple(query_len as i64, db_len as i64, params);
        let length_adj = result.length_adjustment as f64;

        let effective_m = (m - length_adj).max(1.0);
        let effective_n = (n - length_adj).max(1.0);

        Self {
            effective_query_len: effective_m,
            effective_db_len: effective_n,
            effective_space: effective_m * effective_n,
            length_adjustment: result.length_adjustment,
        }
    }

    /// Create search space for database search (multiple subjects) using NCBI BLAST algorithm.
    ///
    /// This is the recommended method for computing E-values that match NCBI BLAST.
    ///
    /// The algorithm properly accounts for:
    /// - The beta parameter in the length adjustment formula
    /// - The number of database sequences (N) in the search space calculation
    /// - Bisection-like iteration for robust convergence
    ///
    /// Formula: ell = alpha/lambda * (logK + log(ss)) + beta
    /// where ss = (m - ell) * (n - N * ell)
    pub fn for_database_search(
        query_len: usize,
        total_db_len: usize,
        num_sequences: usize,
        params: &KarlinParams,
        use_length_adjustment: bool,
    ) -> Self {
        if !use_length_adjustment {
            return Self::simple(query_len, total_db_len);
        }

        let m = query_len as f64;
        let n = total_db_len as f64;
        let n_seqs = num_sequences as f64;

        // Use NCBI-compatible length adjustment with proper N handling
        let result = compute_length_adjustment_ncbi(
            query_len as i64,
            total_db_len as i64,
            num_sequences as i64,
            params,
        );
        let length_adj = result.length_adjustment as f64;

        // For database search, subtract length adjustment from query,
        // and N * length_adjustment from total database length
        let effective_m = (m - length_adj).max(1.0);
        let effective_n = (n - (length_adj * n_seqs)).max(1.0);

        Self {
            effective_query_len: effective_m,
            effective_db_len: effective_n,
            effective_space: effective_m * effective_n,
            length_adjustment: result.length_adjustment,
        }
    }
}

/// Compute effective search space for translated searches (TBLASTX)
///
/// For translated searches, the effective lengths are based on the
/// amino acid sequence lengths (nucleotide length / 3)
pub fn compute_tblastx_search_space(
    query_nucl_len: usize,
    subject_nucl_len: usize,
    params: &KarlinParams,
    use_length_adjustment: bool,
) -> SearchSpace {
    // Convert to amino acid lengths (6 frames each)
    let query_aa_len = query_nucl_len / 3;
    let subject_aa_len = subject_nucl_len / 3;

    // Total search space is 6 * 6 = 36 frame combinations
    // But we use the single-frame length for statistics
    if use_length_adjustment {
        SearchSpace::with_length_adjustment(query_aa_len, subject_aa_len, params)
    } else {
        SearchSpace::simple(query_aa_len, subject_aa_len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_search_space() {
        let ss = SearchSpace::simple(100, 1000);
        assert_eq!(ss.effective_query_len, 100.0);
        assert_eq!(ss.effective_db_len, 1000.0);
        assert_eq!(ss.effective_space, 100000.0);
    }

    #[test]
    fn test_length_adjustment() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        let ss = SearchSpace::with_length_adjustment(100, 1000, &params);

        // Effective lengths should be smaller than original
        assert!(ss.effective_query_len < 100.0);
        assert!(ss.effective_db_len < 1000.0);
        assert!(ss.effective_space < 100000.0);

        // But still positive
        assert!(ss.effective_query_len > 0.0);
        assert!(ss.effective_db_len > 0.0);
        assert!(ss.effective_space > 0.0);
    }

    #[test]
    fn test_length_adjustment_convergence() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        // Test with various sequence lengths using the new NCBI-compatible function
        for (m, n) in [(100, 1000), (1000, 10000), (50, 500)] {
            let result = compute_length_adjustment_simple(m, n, &params);

            // Adjustment should be non-negative
            assert!(result.length_adjustment >= 0);

            // Adjustment should be less than the shorter sequence
            assert!(result.length_adjustment < m.min(n));
        }
    }

    #[test]
    fn test_database_search_space() {
        let params = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        // Test database search with multiple sequences
        let ss = SearchSpace::for_database_search(100, 10000, 10, &params, true);

        // Effective lengths should be smaller than original
        assert!(ss.effective_query_len < 100.0);
        assert!(ss.effective_db_len < 10000.0);
        assert!(ss.effective_space < 1000000.0);

        // Length adjustment should be recorded
        assert!(ss.length_adjustment >= 0);
    }
}
