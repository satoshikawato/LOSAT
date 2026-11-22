// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
// The expected values are the results of the x87 unit itself (gcc long double on x86-64
// Linux, control word 0x037F): `fmul`, `fdiv`, `fadd` of two doubles and `fstpl`.
use super::X87;

fn ld(exponent: i32, significand: u64) -> X87 {
    X87::from_parts(significand, exponent)
}

#[test]
fn mul_rounds_ties_to_even_and_carries() {
    // 1.5 * (1 + 2^-63): the product lies halfway; the kept significand is odd, so up.
    assert_eq!(
        ld(0, 0xC000_0000_0000_0000)
            .mul(ld(0, 0x8000_0000_0000_0001))
            .parts(),
        (false, 0, 0xC000_0000_0000_0002)
    );
    // 1.5 * (1 + 3 * 2^-63): halfway with an even kept significand, so it stays.
    assert_eq!(
        ld(0, 0xC000_0000_0000_0000)
            .mul(ld(0, 0x8000_0000_0000_0003))
            .parts(),
        (false, 0, 0xC000_0000_0000_0004)
    );
    // The product rounds up to 2^64, which carries into the exponent.
    assert_eq!(
        ld(0, 0xFFFF_FFFF_FFFF_FFFE)
            .mul(ld(0, 0x8000_0000_0000_0001))
            .parts(),
        (false, 1, 0x8000_0000_0000_0000)
    );
}

// gcc long double on x86-64 (swebok-improve-20261105/sw8/nfix/probe/x87add/add.c): `fadd`
// of nonnegative values with ties to even, a carry, bits far below the last place, and the
// infinity of an overflow (`None`); `StringToDoublePosix`'s `third` times 10 plus `c`.
#[test]
fn checked_add_and_mul_round_as_fadd_and_fmul_and_overflow() {
    let add = |a: X87, b: X87| a.checked_add(b).map(X87::parts);
    let one = ld(0, 0x8000_0000_0000_0000);
    assert_eq!(
        add(X87::from_u64(u64::MAX), X87::from_u32(1)),
        Some((false, 64, 0x8000_0000_0000_0000))
    );
    assert_eq!(
        add(one, ld(-64, 0x8000_0000_0000_0000)),
        Some((false, 0, 0x8000_0000_0000_0000))
    );
    assert_eq!(
        add(ld(0, 0x8000_0000_0000_0001), ld(-64, 0x8000_0000_0000_0000)),
        Some((false, 0, 0x8000_0000_0000_0002))
    );
    assert_eq!(
        add(one, ld(-64, 0x8000_0000_0800_0000)),
        Some((false, 0, 0x8000_0000_0000_0001))
    );
    assert_eq!(
        add(one, ld(-200, 0x8000_0000_0000_0000)),
        Some((false, 0, 0x8000_0000_0000_0000))
    );
    assert_eq!(add(X87::ZERO, one), Some(one.parts()));
    let max = ld(16383, u64::MAX);
    assert_eq!(add(max, max), None);
    assert!(max.checked_mul(X87::from_u32(2)).is_none());
    assert_eq!(
        max.checked_mul(one).map(X87::parts),
        Some((false, 16383, u64::MAX))
    );
    let ten = X87::from_u32(10);
    let step = |t: X87, c: u32| t.checked_mul(ten).unwrap().checked_add(X87::from_u32(c));
    let t = step(X87::from_u64(1_234_567_890_123_456_789), 7).unwrap();
    assert_eq!(t.parts(), (false, 63, 0xab54_a98c_eb1f_0ad9));
    assert_eq!(
        step(t, 3).map(X87::parts),
        Some((false, 66, 0xd629_d3f0_25e6_cd8f))
    );
}

#[test]
fn div_rounds_to_nearest() {
    // 1 / 1.5 rounds up, 1 / 7 down, 3 / 1.75 up (quotient below 2^64 of the significands).
    assert_eq!(
        ld(0, 0x8000_0000_0000_0000)
            .div(ld(0, 0xC000_0000_0000_0000))
            .parts(),
        (false, -1, 0xAAAA_AAAA_AAAA_AAAB)
    );
    assert_eq!(
        ld(0, 0x8000_0000_0000_0000)
            .div(ld(2, 0xE000_0000_0000_0000))
            .parts(),
        (false, -3, 0x9249_2492_4924_9249)
    );
    assert_eq!(
        ld(1, 0xC000_0000_0000_0000)
            .div(ld(0, 0xE000_0000_0000_0000))
            .parts(),
        (false, 0, 0xDB6D_B6DB_6DB6_DB6E)
    );
    // Quotients of 2^64 and above of the significands drop one more bit.
    assert_eq!(
        ld(0, 0xFFFF_FFFF_FFFF_FFFF)
            .div(ld(0, 0xAAAA_AAAA_AAAA_AAAB))
            .parts(),
        (false, 0, 0xBFFF_FFFF_FFFF_FFFF)
    );
    assert_eq!(
        ld(0, 0xFFFF_FFFF_FFFF_FFFF)
            .div(ld(0, 0xC000_0000_0000_0001))
            .parts(),
        (false, 0, 0xAAAA_AAAA_AAAA_AAA9)
    );
    // -36 / 3 is exact.
    assert_eq!(
        X87::from_f64(-36.0).div(X87::from_f64(3.0)).parts(),
        (true, 3, 0xC000_0000_0000_0000)
    );
    assert_eq!(X87::ZERO.div(X87::from_f64(3.0)).parts(), (false, 0, 0));
}

#[test]
fn nearest_to_sum_rounds_ties_to_even_with_lower_bits_and_borrows() {
    let sum = |hi: f64, lo: f64| X87::nearest_to_sum(hi, lo).parts();
    let p = |e: i32| 2f64.powi(e);
    // 1 + 2^-64: halfway, even kept significand.
    assert_eq!(sum(1.0, p(-64)), (false, 0, 0x8000_0000_0000_0000));
    // 1 + 2^-63 + 2^-64 and (1 + 2^-52) + 2^-63 + 2^-64: halfway, odd, so up.
    assert_eq!(sum(1.0, 1.5 * p(-63)), (false, 0, 0x8000_0000_0000_0002));
    assert_eq!(
        sum(1.0 + p(-52), 1.5 * p(-63)),
        (false, 0, 0x8000_0000_0000_0802)
    );
    // Bits of lo far below the last place (above 2^-126, exact) break the tie.
    assert_eq!(
        sum(1.0, p(-64) + p(-100)),
        (false, 0, 0x8000_0000_0000_0001)
    );
    // A negative lo below a power of two: halfway (to even = 1), exact, below halfway.
    assert_eq!(sum(1.0, -p(-65)), (false, 0, 0x8000_0000_0000_0000));
    assert_eq!(sum(1.0, -p(-64)), (false, -1, 0xFFFF_FFFF_FFFF_FFFF));
    assert_eq!(
        sum(1.0, -(p(-65) + p(-117))),
        (false, -1, 0xFFFF_FFFF_FFFF_FFFF)
    );
    assert_eq!(
        sum(1.0 + p(-52), -(p(-64) + p(-110))),
        (false, 0, 0x8000_0000_0000_07FF)
    );
    // A subnormal hi, and zero (the exponential below the double range).
    assert_eq!(
        sum(f64::from_bits(0x10), 0.0),
        (false, -1070, 0x8000_0000_0000_0000)
    );
    assert_eq!(sum(0.0, 0.0), (false, 0, 0));
}

#[test]
fn to_f64_rounds_ties_to_even_carries_and_gives_subnormals() {
    let store = |exponent: i32, significand: u64| ld(exponent, significand).to_f64().to_bits();
    assert_eq!(store(0, 0x8000_0000_0000_0400), 0x3FF0_0000_0000_0000);
    assert_eq!(store(0, 0x8000_0000_0000_0C00), 0x3FF0_0000_0000_0002);
    assert_eq!(store(0, 0x8000_0000_0000_0401), 0x3FF0_0000_0000_0001);
    assert_eq!(store(0, 0xFFFF_FFFF_FFFF_FC00), 0x4000_0000_0000_0000);
    assert_eq!(store(1023, 0xFFFF_FFFF_FFFF_FC00), 0x7FF0_0000_0000_0000);
    assert_eq!(store(-1074, 0x8000_0000_0000_0000), 0x0000_0000_0000_0001);
    assert_eq!(store(-1075, 0x8000_0000_0000_0000), 0x0000_0000_0000_0000);
    assert_eq!(store(-1074, 0xC000_0000_0000_0000), 0x0000_0000_0000_0002);
    assert_eq!(store(-1075, 0x8000_0000_0000_0001), 0x0000_0000_0000_0001);
    assert_eq!(store(-1023, 0xFFFF_FFFF_FFFF_FFFF), 0x0010_0000_0000_0000);
    assert_eq!(store(-1023, 0xFFFF_FFFF_FFFF_E800), 0x000F_FFFF_FFFF_FFFE);
    assert_eq!(store(-1100, 0x8000_0000_0000_0000), 0x0000_0000_0000_0000);
    assert_eq!(
        X87::from_f64(-1.0)
            .mul(ld(0, 0xC000_0000_0000_0400))
            .to_f64()
            .to_bits(),
        0xBFF8_0000_0000_0000
    );
    assert_eq!(X87::from_f64(-0.0).to_f64().to_bits(), (-0.0f64).to_bits());
}

#[test]
fn f64_round_trips_exactly() {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    for _ in 0..100_000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        // Every finite double, subnormals included (the top exponent field, 0x7FF, is skipped).
        let bits = (state & !(0x7FF << 52)) | (((state >> 52) % 0x7FF) << 52);
        let x = f64::from_bits(bits);
        assert_eq!(X87::from_f64(x).to_f64().to_bits(), bits, "{x:e}");
        if x > 0.0 {
            assert_eq!(
                X87::nearest_to_sum(x, 0.0).parts(),
                X87::from_f64(x).parts()
            );
        }
    }
}
