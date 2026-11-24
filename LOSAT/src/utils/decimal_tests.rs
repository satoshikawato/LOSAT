use super::push_decimal;

#[test]
fn push_decimal_writes_the_digits_display_writes() {
    let mut values = vec![0, 1, -1, 9, 10, -10, 99, 100, 12_345, -987_654_321];
    values.extend([
        i64::from(i32::MIN),
        i64::from(i32::MAX),
        i64::from(u32::MAX),
    ]);
    values.extend([i64::MIN, i64::MIN + 1, i64::MAX]);
    let mut power = 1i64;
    while let Some(next) = power.checked_mul(10) {
        values.extend([power - 1, power, -power, -(power - 1)]);
        power = next;
    }
    for value in values {
        let mut out = b"x".to_vec();
        push_decimal(&mut out, value);
        assert_eq!(out, format!("x{value}").into_bytes(), "{value}");
    }
}
