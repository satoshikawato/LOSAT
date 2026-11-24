//! The decimal digits of an integer appended to bytes, without `core::fmt`.

/// Appends `value` in decimal: its digits, after `-` where it is negative (the `%d` format of
/// `sprintf`, `NStr::NumericToString`, `operator<<`; `Display` writes the same bytes). The one
/// form of the Seq-id labels' numbers (`seqid::label`), `losat blastdbcmd`'s `%` fields and the
/// ISAM keys of makeblastdb (`writedb::isam`): `write!` into a `Vec` took about 450 instructions
/// a number through `core::fmt` and the `io::Write` adapter (line DT, DT1-14).
pub(crate) fn push_decimal(out: &mut Vec<u8>, value: i64) {
    let mut digits = [0u8; 20];
    let mut n = value.unsigned_abs();
    let mut at = digits.len();
    loop {
        at -= 1;
        digits[at] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if value < 0 {
        out.push(b'-');
    }
    out.extend_from_slice(&digits[at..]);
}

#[cfg(test)]
#[path = "decimal_tests.rs"]
mod tests;
