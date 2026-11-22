//! The numeric id list files of the alias files' `GILIST` and `TILIST` keys, read as NCBI's
//! `SeqDB_ReadMemoryGiList` and `SeqDB_ReadMemoryTiList` read them: a binary list (its first byte 0xFF: a marker word, the big-endian count, then one
//! big-endian id per word) or a text list (decimal numbers, each ended by a space, a line feed or
//! a carriage return). The one owner of their layout, and of the lists the alias files attach to
//! the volumes they reach for NCBI's filtered header read (`AliasIdList`).

use crate::seqid::model::ObjectId;
use crate::seqid::ParsedSeqId;

/// The four keys of an alias file that name an id list (`CSeqDB_AliasMask`'s `eGiList`,
/// `eTiList`, `eSiList`, `eTaxIdList`), in the order NCBI's `ComputeMasks` reads them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum IdListKind {
    Gi,
    Ti,
    Si,
    TaxId,
}

impl IdListKind {
    /// Every kind, in NCBI's order.
    pub(crate) const ALL: [IdListKind; 4] = [
        IdListKind::Gi,
        IdListKind::Ti,
        IdListKind::Si,
        IdListKind::TaxId,
    ];

    /// The alias file's key.
    pub(crate) fn key(self) -> &'static str {
        match self {
            IdListKind::Gi => "GILIST",
            IdListKind::Ti => "TILIST",
            IdListKind::Si => "SEQIDLIST",
            IdListKind::TaxId => "TAXIDLIST",
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:2016-2090
    /// The words of NCBI's error for a value with a space: `has multiple <words> (<value>).`.
    pub(crate) fn multiple_words(self) -> &'static str {
        match self {
            IdListKind::Gi => "GI lists",
            IdListKind::Ti => "TI lists",
            IdListKind::Si => "SEQID lists",
            IdListKind::TaxId => "Tax ids lists",
        }
    }
}

/// An id list of the alias files as NCBI attaches it to a volume it reaches
/// (`CSeqDBVol::AttachVolumeGiList`): the ids its filtered header read finds a definition line's
/// Seq-ids in (`CSeqDBGiList::FindId`). Each kind sorted as NCBI's `InsureOrder` leaves it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct AliasIdList {
    /// The GIs of a `GILIST`.
    pub gis: Vec<i64>,
    /// The TIs of a `TILIST`, as NCBI's unsigned `TTi` (a number past `Int8` is none).
    pub tis: Vec<u64>,
    /// The Seq-id strings of a `SEQIDLIST`, as the lookups left them.
    pub sis: Vec<Vec<u8>>,
}

impl AliasIdList {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1633-1665
    // A `ti` general id's string tag is converted by `NStr::StringToInt8`, which throws for a tag
    // that is not a number: LOSAT refuses the line (`Err`). The Seq-id strings are looked up in
    // the forms of `filter::find_si_forms`.
    /// `FindId`: whether the list holds `id`; `Err` for a `ti` id whose tag NCBI's conversion
    /// throws on.
    pub(crate) fn find_id(&self, id: &ParsedSeqId) -> Result<bool, Vec<u8>> {
        match id {
            ParsedSeqId::Gi(gi) => Ok(self.gis.binary_search(gi).is_ok()),
            ParsedSeqId::General(dbtag) if dbtag.db == b"ti" => {
                let ti = match &dbtag.tag {
                    ObjectId::Id(value) => i64::from(*value),
                    ObjectId::Str(text) => {
                        crate::ncbi_str::string_to_int8(text, 0).map_err(|_| text.clone())?
                    }
                };
                Ok(self.tis.binary_search(&(ti as u64)).is_ok())
            }
            _ => Ok(super::filter::find_si_forms(&self.sis, id)),
        }
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:957-968
/// `have_vol` of NCBI's id test of a definition line with the Seq-ids `ids`: one of `lists` holds
/// one of them, asked in the order of the ids; `Err` with the tag of a `ti` id whose conversion
/// throws before one is found.
pub(crate) fn lists_hold(
    lists: &[std::sync::Arc<AliasIdList>],
    ids: &[ParsedSeqId],
) -> Result<bool, Vec<u8>> {
    for id in ids {
        for list in lists {
            if list.find_id(id)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

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

#[cfg(test)]
#[path = "id_list_tests.rs"]
mod tests;
