use super::*;

/// The lines of `bytes` as the test reads them (`alias_file_test`'s loop).
fn read_lines(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut reader = LineReader::new(FastaStream::from_bytes(bytes));
    let mut lines = Vec::new();
    loop {
        lines.push(reader.next_line().to_vec());
        if reader.at_eof() {
            return lines;
        }
    }
}

/// NCBI's `CStreamLineReader` learns the end of line from the first line; the other byte inside a
/// later line cuts it, and the rest after that byte runs into the next line (review F1 of row
/// RG-D2b, rows mix1 and mix2: `unknown keyword encountered: NSEQMAXLEN`).
#[test]
fn the_rest_after_the_other_end_of_line_runs_into_the_next_line() {
    let got = read_lines(b"DBLIST p5\nTITLE a\rNSEQ\nMAXLEN 5\n");
    assert!(got.contains(&b"NSEQMAXLEN 5".to_vec()), "{got:?}");
    let got = read_lines(b"A 1\r\nB 2\r\n");
    let lines: Vec<&[u8]> = got
        .iter()
        .map(Vec::as_slice)
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!(lines, [&b"A 1"[..], b"B 2"]);
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
