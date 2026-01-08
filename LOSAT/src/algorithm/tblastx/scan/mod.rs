//! SIMD-accelerated scanning submodule for TBLASTX
//!
//! Contains offset pair handling, PV test functions, and k-mer index computation.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c

mod offset_pairs;
mod simd_helpers;

pub use offset_pairs::OffsetPair;

#[cfg(target_arch = "x86_64")]
pub use offset_pairs::{copy_offset_pairs_overflow_avx2, copy_offset_pairs_overflow_sse2};

#[cfg(target_arch = "x86_64")]
pub use simd_helpers::{
    compute_3mer_indices_16_avx2, compute_3mer_indices_4_scalar, compute_3mer_indices_8_sse2,
    pv_test_mask16_avx2, pv_test_mask4_avx2, pv_test_mask8_avx2,
};
