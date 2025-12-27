//! Unit tests for blastn/alignment modules

use LOSAT::algorithm::blastn::alignment::utilities::gdb3;

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

