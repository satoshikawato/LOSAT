//! Utility functions for alignment operations

/// GCD for i32 values
pub(crate) fn gcd_i32(a: i32, b: i32) -> i32 {
    let mut a = a.abs();
    let mut b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

/// Signal that a diagonal/offset is invalid
/// Reduces three integers by their GCD if GCD > 1
pub fn gdb3(a: &mut i32, b: &mut i32, c: &mut i32) -> i32 {
    let g = if *b == 0 {
        gcd_i32(*a, *c)
    } else {
        gcd_i32(*a, gcd_i32(*b, *c))
    };
    if g > 1 {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    g
}

/// Greatest common divisor for 64-bit integers.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:383-402
pub(crate) fn gcd_i64(a: i64, b: i64) -> i64 {
    let mut a = a.unsigned_abs();
    let mut b = b.unsigned_abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a as i64
}

/// Reduces three 64-bit integers by their greatest common divisor when it is above 1, widening
/// NCBI's 32-bit routine so large megablast gap costs cannot overflow.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:447-468
pub(crate) fn gdb3_i64(a: &mut i64, b: &mut i64, c: &mut i64) -> i64 {
    let g = if *b == 0 {
        gcd_i64(*a, *c)
    } else {
        gcd_i64(*a, gcd_i64(*b, *c))
    };
    if g > 1 {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    g
}

#[cfg(test)]
#[path = "utilities_tests.rs"]
mod utilities_tests;
