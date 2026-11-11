use super::length_adjustment::compute_length_adjustment_ncbi;
use super::tables::KarlinParams;

/// Effective search space with length-adjusted values
#[derive(Debug, Clone, Copy)]
pub struct SearchSpace {
    /// Effective search space (product of effective lengths)
    pub effective_space: f64,
    /// The length adjustment value used (for debugging/verification)
    pub length_adjustment: i64,
}

impl SearchSpace {
    /// Create a simple search space without length adjustment (BLEMIR default)
    ///
    /// This is the original BLEMIR behavior for backward compatibility.
    /// For NCBI-compatible E-values, use `for_database_search`.
    pub fn simple(query_len: usize, db_len: usize) -> Self {
        let q = query_len as f64;
        let d = db_len as f64;
        Self {
            effective_space: q * d,
            length_adjustment: 0,
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
            effective_space: effective_m * effective_n,
            length_adjustment: result.length_adjustment,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_search_space() {
        let ss = SearchSpace::simple(100, 1000);
        assert_eq!(ss.effective_space, 100.0 * 1000.0);
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
        assert!(ss.length_adjustment > 0);
        assert!(ss.effective_space < 1000000.0);

        // The space is (m - adj) * (n - N * adj)
        let adj = ss.length_adjustment as f64;
        let expected_m = (100.0_f64 - adj).max(1.0);
        let expected_n = (10000.0_f64 - adj * 10.0).max(1.0);
        assert_eq!(ss.effective_space, expected_m * expected_n);
    }

    #[test]
    fn test_effective_search_space_variation() {
        // Test that effective search space changes correctly with sequence length
        // This catches bugs where length adjustment is not applied correctly
        //
        // **Important**: We test the *direction* of change and *sign* of values,
        // not exact values, because NCBI BLAST's effective length calculation
        // is extremely complex and may differ slightly in implementation.

        let params = crate::stats::tables::KarlinParams {
            lambda: 1.28,
            k: 0.46,
            h: 0.85,
            alpha: 1.5,
            beta: -2.0,
        };

        // Test with different sequence lengths
        let test_lengths = vec![
            (100, 1000, 1),      // Short sequences
            (1000, 10000, 5),    // Medium sequences
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

            // The space is (m - adj) * (n - N * adj), as the constructor computes it
            let adj = ss2.length_adjustment as f64;
            let expected_m = (q_len as f64 - adj).max(1.0);
            let expected_n = (db_len as f64 - adj * db_num_seqs as f64).max(1.0);
            assert_eq!(ss2.effective_space, expected_m * expected_n);

            // Length adjustment should be reasonable (not larger than sequence length)
            assert!(ss2.length_adjustment >= 0);
            assert!(ss2.length_adjustment < q_len.min(db_len) as i64);
        }
    }
}
