/// Effective search space with length-adjusted values
#[derive(Debug, Clone, Copy)]
pub struct SearchSpace {
    /// Effective search space (product of effective lengths)
    pub effective_space: f64,
    /// The length adjustment value used (for debugging/verification)
    pub length_adjustment: i64,
}

#[cfg(test)]
#[path = "search_space_ncbi_tests.rs"]
mod search_space_ncbi_tests;
