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


