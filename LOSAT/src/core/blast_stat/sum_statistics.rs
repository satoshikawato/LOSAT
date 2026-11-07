//! P-value and E-value conversion for blastp's composition-based statistics.
//!
//! `p_to_e` and `e_to_p` are the only items of the old sum-statistics layer that
//! are still called (blastp `adjust_evalues_for_composition`). They use libm
//! `ln_1p` / `exp_m1`, where NCBI uses `BLAST_Log1p` / `BLAST_Expm1` (the
//! `stats::p_to_e` / `stats::e_to_p` versions). Switching changes result bits,
//! so it is a separate change.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:174-176 converts with
//! `BLAST_KarlinEtoP` and `BLAST_KarlinPtoE` (c++/src/algo/blast/core/blast_stat.c:4174-4192).

/// Convert P-value to E-value.
pub fn p_to_e(p: f64) -> f64 {
    #[expect(
        clippy::manual_range_contains,
        reason = "NaN: `!(a..=b).contains(&x)` is true for NaN, the comparisons are false"
    )]
    if p < 0.0 || p > 1.0 {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return i32::MAX as f64;
    }
    -(-p).ln_1p()
}

/// Convert E-value to P-value.
pub fn e_to_p(e: f64) -> f64 {
    if e < 0.0 {
        return 0.0;
    }
    -(-e).exp_m1()
}
