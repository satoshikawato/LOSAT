//! The numeric id list files of the alias files' `GILIST`, `TILIST` and `TAXIDLIST` keys, read as
//! NCBI's `SeqDB_ReadMemoryGiList`, `SeqDB_ReadMemoryTiList` and `SeqDB_ReadMemoryTaxIdList` read
//! them: a binary list (its first byte 0xFF: a marker word, the big-endian count, then one
//! big-endian id per word) or a text list (decimal numbers, each ended by a space, a line feed or
//! a carriage return). The one owner of their layout.

use std::collections::BTreeSet;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:854-887
// A first byte that is a digit or '#' makes a text list (where '#' is then an invalid byte), a
// first byte 0xFF in a file of 8 bytes or more a binary one; its fourth byte 0xFE or 0xFC names
// 8-byte ids.
/// Whether the bytes of a numeric list are a binary list (`Some`), and of 8-byte ids (its fourth
/// byte); NCBI's `CSeqDBException` text for an empty or unknown file.
fn binary_kind(bytes: &[u8]) -> Result<Option<bool>, String> {
    let Some(&first) = bytes.first() else {
        return Err("Specified file is empty.".to_string());
    };
    if first.is_ascii_digit() || first == b'#' {
        return Ok(None);
    }
    if bytes.len() >= 8 && first == 0xFF {
        return Ok(Some(matches!(bytes[3], 0xFE | 0xFC)));
    }
    Err("Specified file is not a valid GI/TI list.".to_string())
}

/// The big-endian 32-bit word `at` of `bytes`.
fn word(bytes: &[u8], at: usize) -> u32 {
    u32::from_be_bytes(bytes[4 * at..4 * at + 4].try_into().expect("4 bytes"))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:889-922
/// The digit of a text list's byte, `None` for a separator (a space, a line feed, a carriage
/// return); NCBI's `CSeqDBException` text for any other byte.
fn read_digit(byte: u8, list_type: &str) -> Result<Option<u8>, String> {
    match byte {
        b'0'..=b'9' => Ok(Some(byte - b'0')),
        b' ' | b'\n' | b'\r' => Ok(None),
        _ => Err(format!("Invalid byte in text{list_type} list [{byte}].")),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:983-1001
// The loop keeps a number at the separator after it and none after the last byte: a list whose
// last number has no separator after it loses that number. A 0 is no number.
/// The numbers of a text list in file order, `step` adding a digit to the number read so far.
fn read_text<T: Copy + Default + PartialEq>(
    bytes: &[u8],
    list_type: &str,
    step: impl Fn(T, u8) -> T,
) -> Result<Vec<T>, String> {
    let mut out = Vec::new();
    let mut elem = T::default();
    for &byte in bytes {
        match read_digit(byte, list_type)? {
            Some(digit) => elem = step(elem, digit),
            None => {
                if elem != T::default() {
                    out.push(elem);
                }
                elem = T::default();
            }
        }
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:924-1002
// The binary list's words are counted with the file's length rounded down to whole words; its
// first word is 0xFFFFFFFF and its second the count of the words after the two. A text list's
// number is an unsigned 32-bit value that wraps past 2^32 - 1, as NCBI's `Uint4` does.
/// The GIs of a GI list file (`SeqDB_ReadMemoryGiList`), in file order; NCBI's
/// `CSeqDBException` text for a file that is not one.
pub(crate) fn read_gi_list(bytes: &[u8]) -> Result<Vec<u32>, String> {
    if binary_kind(bytes)?.is_some() {
        let words = bytes.len() / 4;
        let num_gis = words - 2;
        if word(bytes, 0) != 0xFFFF_FFFF || word(bytes, 1) as usize != num_gis {
            return Err("Specified file is not a valid binary GI file.".to_string());
        }
        return Ok((2..words).map(|at| word(bytes, at)).collect());
    }
    read_text(bytes, "GI", |elem: u32, digit| {
        elem.wrapping_mul(10).wrapping_add(u32::from(digit))
    })
}

/// A number of a TI or taxid list that NCBI's signed accumulator cannot hold: NCBI's overflow is
/// undefined, and LOSAT keeps the number it reads past the limit (widened), which no TI or taxid
/// of a database is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ListNumber {
    /// A number in the accumulator's range.
    Value(i64),
    /// A number past it.
    Past,
}

impl Default for ListNumber {
    fn default() -> Self {
        ListNumber::Value(0)
    }
}

/// `elem * 10 + digit` of NCBI's signed accumulator of `max` (`Int4` or `Int8`), widened.
fn widened_step(elem: ListNumber, digit: u8, max: i64) -> ListNumber {
    match elem {
        ListNumber::Value(value) => value
            .checked_mul(10)
            .and_then(|value| value.checked_add(i64::from(digit)))
            .filter(|&value| value <= max)
            .map_or(ListNumber::Past, ListNumber::Value),
        ListNumber::Past => ListNumber::Past,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1149-1282
// A binary TI list's first word is -3 (4-byte TIs) or -4 (8-byte TIs, whose count is that of the
// 8-byte TIs and whose words must be even); its count is compared as an unsigned value. A text
// list's number is an `Int8`.
/// The TIs of a TI list file (`SeqDB_ReadMemoryTiList`), in file order; NCBI's
/// `CSeqDBException` text for a file that is not one.
pub(crate) fn read_ti_list(bytes: &[u8]) -> Result<Vec<ListNumber>, String> {
    if let Some(long_ids) = binary_kind(bytes)? {
        let words = bytes.len() / 4;
        let data = words - 2;
        let num_tis = if long_ids { data / 2 } else { data };
        let marker = word(bytes, 0) as i32;
        let bad = (marker != -3 && marker != -4)
            || word(bytes, 1) as usize != num_tis
            || (data % 2 == 1 && long_ids);
        if bad {
            return Err("Specified file is not a valid binary GI or TI file.".to_string());
        }
        if long_ids {
            return Ok((0..num_tis)
                .map(|at| {
                    let high = u64::from(word(bytes, 2 + 2 * at));
                    let low = u64::from(word(bytes, 3 + 2 * at));
                    ListNumber::Value(((high << 32) | low) as i64)
                })
                .collect());
        }
        return Ok((2..words)
            .map(|at| ListNumber::Value(i64::from(word(bytes, at) as i32)))
            .collect());
    }
    read_text(bytes, "TI", |elem, digit| {
        widened_step(elem, digit, i64::MAX)
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1080-1121
// A binary list's ids are big-endian `Int4`s; a text list's number is an `Int4`.
/// The taxids of a taxid list file (`SeqDB_ReadMemoryTaxIdList`), as the set NCBI keeps;
/// NCBI's `CSeqDBException` text for a file that is not one.
pub(crate) fn read_tax_id_list(bytes: &[u8]) -> Result<BTreeSet<ListNumber>, String> {
    if binary_kind(bytes)?.is_some() {
        let words = bytes.len() / 4;
        let num_taxids = words - 2;
        if word(bytes, 0) != 0xFFFF_FFFF || word(bytes, 1) as i32 != num_taxids as i32 {
            return Err("Specified file is not a valid binary Tax Id List file.".to_string());
        }
        return Ok((2..words)
            .map(|at| ListNumber::Value(i64::from(word(bytes, at) as i32)))
            .collect());
    }
    Ok(read_text(bytes, "TAXID", |elem, digit| {
        widened_step(elem, digit, i64::from(i32::MAX))
    })?
    .into_iter()
    .collect())
}

#[cfg(test)]
#[path = "id_list_tests.rs"]
mod tests;
