//! BLAST Utility Functions
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c
//!
//! This module provides:
//! - `BLASTAA_SIZE`, the size of the amino acid alphabet of NCBI's lookup tables
//! - `BLAST_SCORE_MIN`, the score of the `'-'` row and column of a matrix
//! - `ncbi_int4_from_double`, the `(Int4)` cast of a `double` that `BLAST_Cutoffs` uses
//!
//! The BLOSUM62 matrix and its NCBIstdaa score lookup live in `utils::matrix` (`blosum62_score`).

/// NCBI BLASTAA_SIZE - size of amino acid alphabet for lookup tables
pub const BLASTAA_SIZE: usize = 28;

/// NCBI blast core keeps the `'-'`/NULLB row and column at BLAST_SCORE_MIN.
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:121-122
///   #define BLAST_SCORE_MIN INT2_MIN
pub const BLAST_SCORE_MIN: i32 = i16::MIN as i32;

/// NCBI's `(Int4)` conversion of a double, as the x86-64 build of NCBI BLAST+ 2.17.0 (the
/// oracle) performs it with `cvttsd2si`: truncation toward zero, and `INT_MIN` for NaN and
/// for values outside `Int4` (undefined in C). Rust's `as` saturates instead.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4047-4049,4059-4061
/// ```c
/// /* Smallest float that might not cause a floating point exception in
///    S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda )); below.  */
///    const double kSmallFloat = 1.0e-297;
///    ...
///    E = MAX(E, kSmallFloat);
///
///    S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ));
/// ```
/// With `E` at `kSmallFloat`, `K * searchsp / E` overflows to infinity once `K * searchsp`
/// exceeds about 1.8e11; `S` is then `INT_MIN`, and `BLAST_Cutoffs` keeps its initial
/// cutoff of 1 (the S07+ fifteenth audit round).
pub fn ncbi_int4_from_double(value: f64) -> i32 {
    let truncated = value.trunc();
    if truncated.is_nan() || truncated < f64::from(i32::MIN) || truncated > f64::from(i32::MAX) {
        i32::MIN
    } else {
        truncated as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn int4_conversion_is_the_x86_64_one() {
        assert_eq!(ncbi_int4_from_double(41.9), 41);
        assert_eq!(ncbi_int4_from_double(-41.9), -41);
        assert_eq!(ncbi_int4_from_double(f64::from(i32::MAX)), i32::MAX);
        for value in [
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            2_147_483_648.0,
            -2_147_483_649.0,
        ] {
            assert_eq!(ncbi_int4_from_double(value), i32::MIN, "{value}");
        }
    }
}
