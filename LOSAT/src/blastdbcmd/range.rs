//! `-range` and the range of a batch line: NCBI's `ParseSequenceRangeOpenEnd` and the
//! clipping of the range to a sequence (`CBlastDB_FastaFormatter::Write`,
//! `CBlastDB_SeqFormatter::x_GetSeq`).

/// `TSeqRange` as blastdbcmd holds it: 0-based, both ends inclusive, `to` unset (`u32::MAX`)
/// for an open end; `from` above `to` (`TSeqRange::GetEmpty()`, the range of no `-range`) and
/// the open range from `u32::MAX` (`0-`) are empty ([`SeqRange::is_empty`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct SeqRange {
    pub from: u32,
    pub to: u32,
}

impl SeqRange {
    // NCBI reference (598d8ae6): c++/include/util/range.hpp:148-155
    /// `TSeqRange::Empty()`: the open end at or below the start. A closed range is empty where
    /// `from` is above `to` (the range of no `-range`, `CRange<TSeqPos>::GetEmpty()`, has `from`
    /// at the maximum and `to` one below). An open one (`to` unset, NCBI's open end at the
    /// maximum) is empty only from the maximum: `0-`, whose start `ParseSequenceRangeOpenEnd`
    /// stores as -1, so `x_GetSeq` reads no length for it.
    pub fn is_empty(self) -> bool {
        self.from > self.to || self.from == u32::MAX
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:285-291
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:477-486
    // NCBI reference (598d8ae6): c++/include/util/range.hpp:148-155
    // NCBI reference (598d8ae6): c++/include/util/range.hpp:180-183
    /// The range clipped to a sequence of `length` letters: `None` where NCBI then takes the
    /// whole sequence (no range, a zero-length sequence where `SetTo` sets the open end to zero,
    /// or a range whose start is at or past the end).
    pub fn clipped(self, length: u32) -> Option<(u32, u32)> {
        if length == 0 || self.is_empty() {
            return None;
        }
        let to = if length <= self.to {
            length - 1
        } else {
            self.to
        };
        (self.from <= to).then_some((self.from, to))
    }
}

/// `CRange<TSeqPos>::GetEmpty()`: the range of no `-range` (`m_SeqRange` is set from
/// `TSeqRange::GetEmpty`).
pub(super) const NO_RANGE: SeqRange = SeqRange {
    from: u32::MAX,
    to: u32::MAX - 1,
};

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643
/// `NStr::StringToInt` of `text` (`ncbi_str::string_to_int`): the number, or the message of NCBI's
/// exception as `GetMsg` gives it (`Int8` for the text that is not a number or passes a 64-bit
/// integer, `int` for one that fits 64 bits but not 32).
fn string_to_int(text: &[u8]) -> Result<i32, String> {
    crate::ncbi_str::string_to_int(text, 0)
        .map_err(|error| String::from_utf8_lossy(error.get_msg()).into_owned())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:181-219
/// `ParseSequenceRangeOpenEnd` of `text`: the range, or the message of the exception NCBI
/// throws (a `CBlastException` of its own, or the `CStringException` of `StringToInt`).
/// The first and the last `-`-separated pieces are the start and the end (a middle piece is
/// ignored); an empty last piece is an open end. The `TSeqRange` is built from `int`s:
/// `0-` passes the check of the end and stores the start as `-1`, which the range holds as
/// the largest `TSeqPos` (so the clipped range is empty and the whole sequence is written).
///
/// `-range ""` is not reproduced: NCBI reads `tokens.front()` of an empty vector (it
/// crashed, exit status 139, probe D61).
pub(super) fn parse_sequence_range_open_end(text: impl AsRef<[u8]>) -> Result<SeqRange, String> {
    const PREFIX: &str = "Failed to parse sequence range";
    let tokens: Vec<&[u8]> = text.as_ref().split(|&byte| byte == b'-').collect();
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
    // A decrement of `from` on an `int`, stored as a `TSeqPos` (unsigned 32-bit).
    range.from = (from.wrapping_sub(1)) as u32;
    Ok(range)
}

#[cfg(test)]
#[path = "range_tests.rs"]
mod tests;
