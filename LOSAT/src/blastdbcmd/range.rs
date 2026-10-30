//! `-range` and the range of a batch line: NCBI's `ParseSequenceRangeOpenEnd` and the
//! clipping of the range to a sequence (`CBlastDB_FastaFormatter::Write`,
//! `CBlastDB_SeqFormatter::x_GetSeq`).

use super::args;

/// `TSeqRange` as blastdbcmd holds it: 0-based, both ends inclusive, `to` unset (`u32::MAX`)
/// for an open end; `from` above `to` is the empty range (`TSeqRange::GetEmpty()`, the range
/// of no `-range`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct SeqRange {
    pub from: u32,
    pub to: u32,
}

impl SeqRange {
    /// `TSeqRange::GetEmpty()`.
    pub const EMPTY: SeqRange = SeqRange {
        from: u32::MAX,
        to: u32::MAX,
    };

    /// `TSeqRange::Empty()`: `from` above `to`. The empty range has `from == to`
    /// (`GetEmpty` is `[max, max)` with the open end as `max - 1`: `TSeqRange::GetToOpen`
    /// unset); NCBI's `NotEmpty` is what `-range` tests, and the range of no `-range` holds
    /// `to` below `from` there (`CRange<TSeqPos>::GetEmpty()` is `from = max, to = max - 1`).
    pub fn is_empty(self) -> bool {
        self.from > self.to
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:477-483
    // ```c++
    //     CRef<CSeq_loc> range;
    //     if(config.m_SeqRange.NotEmpty()) {
    //     	TSeqPos length = m_BlastDb.GetSeqLength(oid);
    //     	TSeqRange r = config.m_SeqRange;
    //         if((TSeqPos)length <= config.m_SeqRange.GetTo()) {
    //         	r.SetTo(length-1);
    //     	}
    //         if(!r.Empty()) {
    // ```
    /// The range clipped to a sequence of `length` letters: `None` where NCBI then takes the
    /// whole sequence (no range, or a range whose start is at or past the end).
    pub fn clipped(self, length: u32) -> Option<(u32, u32)> {
        if self.is_empty() {
            return None;
        }
        let to = if length <= self.to {
            length.wrapping_sub(1)
        } else {
            self.to
        };
        (self.from <= to).then_some((self.from, to))
    }
}

/// `CRange<TSeqPos>::GetEmpty()`: the range of no `-range` (`m_SeqRange = TSeqRange::GetEmpty()`).
pub(super) const NO_RANGE: SeqRange = SeqRange {
    from: u32::MAX,
    to: u32::MAX - 1,
};

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1246-1270
// ```c++
// int NStr::StringToInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Int8 value = StringToInt8(str, flags, base);
//     if ( value < kMin_Int  ||  value > kMax_Int ) {
//         S2N_CONVERT_ERROR_OVERFLOW(int);
//     }
//     return (int) value;
// }
// ```
/// `NStr::StringToInt(text)`: the number, or the message of NCBI's exception (`Int8` for the
/// text that is not a number or passes a 64-bit integer, `int` for one that fits 64 bits but
/// not 32).
fn string_to_int(text: &str) -> Result<i32, String> {
    let number = args::string_to_int8_message(text)?;
    i32::try_from(number).map_err(|_| {
        format!(
            "Cannot convert string '{}' to int, overflow",
            args::printable_string(text)
        )
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:181-219
// ```c++
// TSeqRange
// ParseSequenceRangeOpenEnd(const string& range_str,
//                    	   	  const char* error_prefix /* = NULL */)
// {
//     static const char* kDfltErrorPrefix = "Failed to parse sequence range";
//     static const string kDelimiters("-");
//     string error_msg(error_prefix ? error_prefix : kDfltErrorPrefix);
//
//     vector<string> tokens;
//     NStr::Split(range_str, kDelimiters, tokens);
//     if (tokens.front().empty()) {
//         error_msg += " (start cannot be empty)";
//         NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
//     }
//
//     TSeqRange retval;
//     int from = NStr::StringToInt(tokens.front());
//     int to = 0;
//
//     if(!tokens.back().empty()) {
//     	to = NStr::StringToInt(tokens.back());
//
//     	if (from <= 0 || to <= 0) {
//     		error_msg += " (range elements cannot be less than or equal to 0)";
//     		NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
//     	}
//     	if (from > to) {
//     	    error_msg += " (start cannot be larger than stop)";
//     	    NCBI_THROW(CBlastException, eInvalidArgument, error_msg);
//     	}
//     	to--;
//     	retval.SetTo(to);
//     }
//     //Note: TSeqRange is defaulted to max value.
//
//     from--;
//     retval.SetFrom(from);
//     return retval;
// }
// ```
/// `ParseSequenceRangeOpenEnd(text)`: the range, or the message of the exception NCBI
/// throws (a `CBlastException` of its own, or the `CStringException` of `StringToInt`).
/// The first and the last `-`-separated pieces are the start and the end (a middle piece is
/// ignored); an empty last piece is an open end. The `TSeqRange` is built from `int`s:
/// `0-` passes the check of the end and stores the start as `-1`, which the range holds as
/// the largest `TSeqPos` (so the clipped range is empty and the whole sequence is written).
///
/// `-range ""` is not reproduced: NCBI reads `tokens.front()` of an empty vector (it
/// crashed, exit status 139, probe D61).
pub(super) fn parse_sequence_range_open_end(text: &str) -> Result<SeqRange, String> {
    const PREFIX: &str = "Failed to parse sequence range";
    let tokens: Vec<&str> = text.split('-').collect();
    let (Some(&front), Some(&back)) = (tokens.first(), tokens.last()) else {
        unreachable!("split gives at least one piece");
    };
    if front.is_empty() {
        return Err(format!("{PREFIX} (start cannot be empty)"));
    }
    let from = string_to_int(front)?;
    let mut range = SeqRange {
        from: 0,
        to: u32::MAX,
    };
    if !back.is_empty() {
        let to = string_to_int(back)?;
        if from <= 0 || to <= 0 {
            return Err(format!(
                "{PREFIX} (range elements cannot be less than or equal to 0)"
            ));
        }
        if from > to {
            return Err(format!("{PREFIX} (start cannot be larger than stop)"));
        }
        range.to = (to - 1) as u32;
    }
    // `from--` on an `int`, stored as a `TSeqPos` (unsigned 32-bit).
    range.from = (from.wrapping_sub(1)) as u32;
    Ok(range)
}

#[cfg(test)]
mod tests {
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
}
