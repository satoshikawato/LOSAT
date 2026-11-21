use super::*;

#[test]
fn the_grammar_is_ncbi_s() {
    let range = |text: &str| parse_sequence_range_open_end(text);
    assert_eq!(range("10-30"), Ok(SeqRange { from: 9, to: 29 }));
    assert_eq!(range("5"), Ok(SeqRange { from: 4, to: 4 }));
    assert_eq!(range("5-5"), Ok(SeqRange { from: 4, to: 4 }));
    assert_eq!(range("1-2-3"), Ok(SeqRange { from: 0, to: 2 }));
    assert_eq!(range("5--7"), Ok(SeqRange { from: 4, to: 6 }));
    assert_eq!(range("+5-7"), Ok(SeqRange { from: 4, to: 6 }));
    assert_eq!(
        range("130-"),
        Ok(SeqRange {
            from: 129,
            to: u32::MAX
        })
    );
    // `0-`: the start is stored as the largest TSeqPos; the clipped range is empty.
    let zero_open = range("0-").unwrap();
    assert_eq!(zero_open.from, u32::MAX);
    assert_eq!(zero_open.clipped(137), None);
    assert_eq!(
        range("0-10").unwrap_err(),
        "Failed to parse sequence range (range elements cannot be less than or equal to 0)"
    );
    assert_eq!(
        range("5-0").unwrap_err(),
        "Failed to parse sequence range (range elements cannot be less than or equal to 0)"
    );
    assert_eq!(
        range("20-10").unwrap_err(),
        "Failed to parse sequence range (start cannot be larger than stop)"
    );
    assert_eq!(
        range("-5").unwrap_err(),
        "Failed to parse sequence range (start cannot be empty)"
    );
    assert_eq!(
        range("abc").unwrap_err(),
        "Cannot convert string 'abc' to Int8"
    );
    assert_eq!(
        range("5-x").unwrap_err(),
        "Cannot convert string 'x' to Int8"
    );
    assert_eq!(
        range("3000000000-5").unwrap_err(),
        "Cannot convert string '3000000000' to int, overflow"
    );
    assert_eq!(
        range("99999999999999999999-5").unwrap_err(),
        "Cannot convert string '99999999999999999999' to Int8, overflow"
    );
}

#[test]
fn clipping_follows_the_formatters() {
    assert_eq!(NO_RANGE.clipped(100), None);
    assert_eq!(SeqRange { from: 4, to: 19 }.clipped(488), Some((4, 19)));
    // `130-` of a 142-letter sequence; `200-` of one of 142.
    assert_eq!(
        SeqRange {
            from: 129,
            to: u32::MAX
        }
        .clipped(142),
        Some((129, 141))
    );
    assert_eq!(
        SeqRange {
            from: 199,
            to: u32::MAX
        }
        .clipped(142),
        None
    );
    // `488-600` of a 488-letter sequence: the last letter.
    assert_eq!(
        SeqRange { from: 487, to: 599 }.clipped(488),
        Some((487, 487))
    );
}
