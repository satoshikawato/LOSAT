//! Unit tests for blastn/alignment modules

use crate::algorithm::blastn::alignment::utilities::{gcd_i64, gdb3, gdb3_i64};

#[test]
fn test_gdb3_no_reduction() {
    // Numbers with GCD = 1 should not be reduced
    let mut a = 3;
    let mut b = 5;
    let mut c = 7;

    let g = gdb3(&mut a, &mut b, &mut c);

    assert_eq!(g, 1);
    assert_eq!(a, 3);
    assert_eq!(b, 5);
    assert_eq!(c, 7);
}

#[test]
fn test_gdb3_with_reduction() {
    // Numbers with GCD > 1 should be reduced
    let mut a = 6;
    let mut b = 9;
    let mut c = 12;

    let g = gdb3(&mut a, &mut b, &mut c);

    assert_eq!(g, 3);
    assert_eq!(a, 2);
    assert_eq!(b, 3);
    assert_eq!(c, 4);
}

#[test]
fn test_gdb3_zero_b() {
    // When b = 0, should compute GCD of a and c
    let mut a = 8;
    let mut b = 0;
    let mut c = 12;

    let g = gdb3(&mut a, &mut b, &mut c);

    assert_eq!(g, 4);
    assert_eq!(a, 2);
    assert_eq!(b, 0);
    assert_eq!(c, 3);
}

#[test]
fn test_gdb3_all_zero() {
    // All zeros should return GCD = 0 (or handle gracefully)
    let mut a = 0;
    let mut b = 0;
    let mut c = 0;

    let g = gdb3(&mut a, &mut b, &mut c);

    // GCD of zeros is typically 0 or undefined
    // The function should handle this without panicking
    assert!(g >= 0);
}

#[test]
fn test_gdb3_negative_numbers() {
    // Should handle negative numbers (GCD uses absolute values)
    let mut a = -6;
    let mut b = 9;
    let mut c = -12;

    let g = gdb3(&mut a, &mut b, &mut c);

    assert_eq!(g, 3);
    // After reduction, signs may be preserved or normalized
    // The exact behavior depends on implementation
}

#[test]
fn test_gdb3_large_numbers() {
    // Test with larger numbers
    let mut a = 100;
    let mut b = 150;
    let mut c = 200;

    let g = gdb3(&mut a, &mut b, &mut c);

    assert_eq!(g, 50);
    assert_eq!(a, 2);
    assert_eq!(b, 3);
    assert_eq!(c, 4);
}

#[test]
fn test_gdb3_prime_numbers() {
    // Prime numbers should have GCD = 1
    let mut a = 17;
    let mut b = 19;
    let mut c = 23;

    let g = gdb3(&mut a, &mut b, &mut c);

    assert_eq!(g, 1);
    assert_eq!(a, 17);
    assert_eq!(b, 19);
    assert_eq!(c, 23);
}

#[test]
fn test_gdb3_i64_large_gap_costs() {
    // Doubled INT32_MAX gap costs with odd reward (1/-2): no 32-bit overflow.
    let mut op = 6i64;
    let mut go = 2 * i64::from(i32::MAX);
    let mut ge = 2 * i64::from(i32::MAX) + 1;
    let g = gdb3_i64(&mut op, &mut go, &mut ge);
    assert_eq!(g, 1);
    assert_eq!((op, go, ge), (6, 4_294_967_294, 4_294_967_295));

    // Scaled nucleotide scoring (reward 1000, penalty -2000, gapopen 40000, gapextend 2000):
    // op_cost = 3000, gap_open = 40000, gap_extend = 2500 -> gcd = 500.
    let (mut op, mut go, mut ge) = (3000i64, 40_000i64, 2500i64);
    let g = gdb3_i64(&mut op, &mut go, &mut ge);
    assert_eq!(g, 500);
    assert_eq!((op, go, ge), (6, 80, 5));

    // Zero gap_open with huge gap_extend reduces cleanly.
    let (mut op, mut go, mut ge) = (6i64, 0i64, 4_000_000_002i64);
    let g = gdb3_i64(&mut op, &mut go, &mut ge);
    assert_eq!(g, 6);
    assert_eq!((op, go, ge), (1, 0, 666_666_667));

    assert_eq!(gcd_i64(4_000_000_000, 6), 2);
}
