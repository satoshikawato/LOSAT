use super::*;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/ncbi_std.h:112
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1739-1749
#[test]
fn the_best_evalue_is_ncbis_min_from_int4_max_so_a_nan_takes_its_place() {
    assert_eq!(best_evalue([]), 2147483647.0);
    assert_eq!(best_evalue([3.0, 1e-5, 2.0]), 1e-5);
    assert_eq!(best_evalue([0.5, 0.5]), 0.5);
    // `MIN(evalue, best)` keeps `best` only when the e-value is above it: a NaN e-value becomes
    // the best, and the e-value after a NaN best replaces it (`<` from INT4_MAX, or `f64::min`,
    // kept INT4_MAX and 1e-5).
    assert!(best_evalue([f64::NAN]).is_nan());
    assert!(best_evalue([1e-5, f64::NAN]).is_nan());
    assert_eq!(best_evalue([f64::NAN, 1e-5]), 1e-5);
    assert_eq!(best_evalue([f64::NAN, 1e12]), 1e12);
    for (best, evalue) in [(1.0, 2.0), (2.0, 1.0), (1.0, 1.0), (0.0, -0.0), (-0.0, 0.0)] {
        assert_eq!(
            ncbi_min_evalue(best, evalue).to_bits(),
            (if evalue > best { best } else { evalue }).to_bits()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1996
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1567
#[test]
fn the_reap_keeps_a_nan_and_the_early_termination_ends_one() {
    assert!(evalue_not_above(1.0, 1.0));
    assert!(!evalue_not_above(1.5, 1.0));
    assert!(evalue_not_above(f64::NAN, 1.0));
    assert!(evalue_not_above(1e300, f64::NAN));
    // The early termination goes on (does not end the match) only for `evalue <= limit`: a NaN
    // best e-value ends it, as INT4_MAX did only under a limit below it (an `-evalue` below
    // about 4.29e8, five times which is INT4_MAX).
    assert!(!evalue_not_at_most(5.0, 5.0));
    assert!(evalue_not_at_most(5.5, 5.0));
    assert!(evalue_not_at_most(f64::NAN, 5.0));
    assert!(evalue_not_at_most(f64::NAN, 5.0 * 1e10));
    assert!(!evalue_not_at_most(best_evalue([]), 5.0 * 1e9));
}
