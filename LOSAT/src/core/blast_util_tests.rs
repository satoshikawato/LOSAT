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
