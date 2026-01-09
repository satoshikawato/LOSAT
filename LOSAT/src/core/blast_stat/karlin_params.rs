//! Karlin-Altschul statistical parameters.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c

/// Karlin-Altschul statistical parameters
#[derive(Debug, Clone, Copy)]
pub struct KarlinParams {
    /// Lambda parameter for bit score calculation
    pub lambda: f64,
    /// K parameter for E-value calculation
    pub k: f64,
    /// H parameter (entropy) for length adjustment
    pub h: f64,
    /// Alpha parameter for length correction mean
    pub alpha: f64,
    /// Beta parameter for length correction
    pub beta: f64,
}

impl Default for KarlinParams {
    fn default() -> Self {
        Self {
            lambda: 0.625,
            k: 0.041,
            h: 0.85,
            alpha: 1.5,
            beta: -2.0,
        }
    }
}
