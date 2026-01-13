//! Search space definitions for BLAST.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c

use super::karlin_params::KarlinParams;
use super::length_adjustment::{compute_length_adjustment_ncbi, compute_length_adjustment_simple};

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
    pub fn with_length_adjustment(query_len: usize, db_len: usize, params: &KarlinParams) -> Self {
        let m = query_len as f64;
        let n = db_len as f64;

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

        let result = compute_length_adjustment_ncbi(
            query_len as i64,
            total_db_len as i64,
            num_sequences as i64,
            params,
        );
        let length_adj = result.length_adjustment as f64;

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
pub fn compute_tblastx_search_space(
    query_nucl_len: usize,
    subject_nucl_len: usize,
    params: &KarlinParams,
    use_length_adjustment: bool,
) -> SearchSpace {
    let query_aa_len = query_nucl_len / 3;
    let subject_aa_len = subject_nucl_len / 3;

    if use_length_adjustment {
        SearchSpace::with_length_adjustment(query_aa_len, subject_aa_len, params)
    } else {
        SearchSpace::simple(query_aa_len, subject_aa_len)
    }
}
