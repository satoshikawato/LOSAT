use super::*;

#[test]
fn every_carriage_return_and_line_feed_ends_a_line() {
    let got: Vec<&[u8]> = lines(b"A 1\r\nB 2\rC 3\nD")
        .filter(|line| !line.is_empty())
        .collect();
    assert_eq!(got, [&b"A 1"[..], b"B 2", b"C 3", b"D"]);
}

/// Row a_ov1: FIRST_OID 65535 and LAST_OID 65535 multiply to 2^32, 0 in `int`: the exact and the
/// wrapped products decide the pairing check apart, which LOSAT refuses until Q-D2b-OV.
#[test]
fn an_overflow_that_decides_a_check_is_refused() {
    let product = Both::int(65536).mul(Both::int(65536), false);
    assert_eq!(product.wrapped, 0);
    assert!(decide(product, |v| v == 0, b"ov1.pal", "FIRST_OID and LAST_OID").is_err());
}

/// Row a_ov3: FIRST_OID 2147483647 wraps to INT_MIN after `+ 1`, but neither product is 0: the
/// check gives NCBI's result.
#[test]
fn an_overflow_that_decides_nothing_gives_the_result() {
    let product = Both::int(i32::MAX)
        .add(Both::int(1), false)
        .mul(Both::int(1), false);
    assert_eq!(product.wrapped, i64::from(i32::MIN));
    assert_eq!(
        decide(product, |v| v == 0, b"ov3.pal", "FIRST_OID and LAST_OID").ok(),
        Some(false)
    );
}

/// The OIDLIST product is taken in `Int8` after 1 is added to MAXOID in `int` (row a_ov2:
/// 2^16 * 2^32 * 2^16 wraps to 0).
#[test]
fn the_oidlist_product_wraps_in_int8() {
    let product = Both::int(65535)
        .add(Both::int(1), false)
        .mul(Both::int8(4_294_967_295).add(Both::int8(1), true), true)
        .mul(Both::int(65535).add(Both::int(1), false), true);
    assert_eq!(product.wrapped, 0);
    assert_eq!(product.exact, 1 << 64);
}
