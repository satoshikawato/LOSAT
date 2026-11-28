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
    // `0-`: the start is stored as the largest TSeqPos; the range is empty (NCBI's `Empty`), so
    // `x_GetSeq` reads no length, and the clipped range is empty.
    let zero_open = range("0-").unwrap();
    assert_eq!(zero_open.from, u32::MAX);
    assert!(zero_open.is_empty());
    assert_eq!(zero_open.clipped(137), None);
    assert!(!range("1-").unwrap().is_empty());
    assert!(!range("10-30").unwrap().is_empty());
    assert!(NO_RANGE.is_empty());
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
    // Zero-length sequence: NCBI's `SetTo` sets the open end to zero, making `Empty` true.
    assert_eq!(SeqRange { from: 0, to: 9 }.clipped(0), None);
    assert_eq!(
        SeqRange {
            from: 0,
            to: u32::MAX
        }
        .clipped(0),
        None
    );
    assert_eq!(
        SeqRange {
            from: u32::MAX,
            to: u32::MAX
        }
        .clipped(0),
        None
    );
    // `2^31 + 10` protein sequence: `0-` clips to `None`; open tail clips to `length - 1`.
    let over_2gib = (1u32 << 31) + 10;
    assert_eq!(
        SeqRange {
            from: u32::MAX,
            to: u32::MAX
        }
        .clipped(over_2gib),
        None
    );
    assert_eq!(
        SeqRange {
            from: 2_147_483_646,
            to: u32::MAX
        }
        .clipped(over_2gib),
        Some((2_147_483_646, over_2gib - 1))
    );
}
