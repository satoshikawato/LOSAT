//! A seqid list file (line DI-8b: `-seqidlist`, `-negative_seqidlist`), read as NCBI's
//! `CSeqDBFileGiList(fname, eSiList)` reads it (`SeqDB_ReadSiList`): a binary list made by
//! `blastdb_aliastool -seqid_file_in` (its first byte is NUL; `CSeqidlistRead`), else text
//! (`SeqDB_ReadMemorySiList`: white-space separated tokens, a leading `>` skipped, `#` to
//! the end of the line a comment). The ids of a text list are then parsed and rewritten
//! for a version 5 database as `s_ProcessSeqIDsForV5` rewrites them: `gi|N` dropped, `pir`
//! and `prf` ids in FASTA form, every other id as `CSeq_id::GetSeqIdString(true)` (the
//! upper-cased accession with its version when given), sorted and unique.

use crate::seqid::accguide::parse_flags::{ANY_LOCAL, PARTIAL_OK, RAW_TEXT};
use crate::seqid::accguide::SeqIdWarning;
use crate::seqid::label;
use crate::seqid::parse::{seq_id_from_text, SeqIdError};
use crate::seqid::ParsedSeqId;

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:134-145
/// What the list file says of itself (`SBlastSeqIdListInfo`): a text list is `is_v4`, a binary
/// list carries its header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqidListInfo {
    /// A text list (`is_v4 = true`, the constructor's default); `false` for a binary list.
    pub is_v4: bool,
    pub file_size: u64,
    pub num_ids: u64,
    pub title: Vec<u8>,
    pub create_date: Vec<u8>,
    /// The total length of the database the binary list was made for (0: none recorded).
    pub db_vol_length: u64,
    pub db_create_date: Vec<u8>,
    pub db_vol_names: Vec<u8>,
}

impl Default for SeqidListInfo {
    fn default() -> Self {
        SeqidListInfo {
            is_v4: true,
            file_size: 0,
            num_ids: 0,
            title: Vec::new(),
            create_date: Vec::new(),
            db_vol_length: 0,
            db_create_date: Vec::new(),
            db_vol_names: Vec::new(),
        }
    }
}

/// A seqid list file as read: its ids as the file gives them and its header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqidListFile {
    /// The ids: the tokens of a text list (as written), the ids of a binary list.
    pub ids: Vec<Vec<u8>>,
    pub info: SeqidListInfo,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:58-85
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:87-108
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:111-120
// A 0-byte file has no memory map (`CMemoryFile` leaves the pointer null), hence the first
// error. NCBI reads the header and the ids from the map without bounds; a header or an id
// cut off by the end of the file reads memory past it. LOSAT gives the second error for such
// a file instead.
/// The header of a binary list and its ids: `Ok(None)` for a file that is not binary (its
/// first byte is not NUL), `Err` with NCBI's message for a damaged one.
fn read_binary(bytes: &[u8]) -> Result<Option<(SeqidListInfo, Vec<Vec<u8>>)>, String> {
    const INVALID: &str = "Invalid seqidlist file";
    if bytes.first() != Some(&0) {
        return Ok(None);
    }
    let mut cursor = Cursor { bytes, at: 1 };
    let mut info = SeqidListInfo {
        is_v4: false,
        ..SeqidListInfo::default()
    };
    info.file_size = cursor.u64().ok_or(INVALID)?;
    if info.file_size != bytes.len() as u64 {
        return Err(INVALID.to_string());
    }
    info.num_ids = cursor.u64().ok_or(INVALID)?;
    let title_length = cursor.u32().ok_or(INVALID)? as usize;
    info.title = cursor.take(title_length).ok_or(INVALID)?.to_vec();
    let create_date_length = usize::from(cursor.u8().ok_or(INVALID)?);
    info.create_date = cursor.take(create_date_length).ok_or(INVALID)?.to_vec();
    info.db_vol_length = cursor.u64().ok_or(INVALID)?;
    if info.db_vol_length != 0 {
        let db_create_date_length = usize::from(cursor.u8().ok_or(INVALID)?);
        info.db_create_date = cursor.take(db_create_date_length).ok_or(INVALID)?.to_vec();
        let vol_names_length = cursor.u32().ok_or(INVALID)? as usize;
        info.db_vol_names = cursor.take(vol_names_length).ok_or(INVALID)?.to_vec();
    }
    let mut ids = Vec::with_capacity(usize::try_from(info.num_ids).unwrap_or(0).min(1 << 20));
    while cursor.at < bytes.len() && (ids.len() as u64) < info.num_ids {
        let id_len = cursor.u8().ok_or(INVALID)?;
        let len = if id_len == 0xFF {
            cursor.u32().ok_or(INVALID)? as usize
        } else {
            usize::from(id_len)
        };
        ids.push(cursor.take(len).ok_or(INVALID)?.to_vec());
    }
    if ids.len() as u64 != info.num_ids {
        return Err("Invalid total num of ids in seqidlist file".to_string());
    }
    Ok(Some((info, ids)))
}

/// A read position in the binary list (`m_Ptr`); every read is bounded by the file's end.
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, len: usize) -> Option<&'a [u8]> {
        let end = self
            .at
            .checked_add(len)
            .filter(|&end| end <= self.bytes.len())?;
        let piece = &self.bytes[self.at..end];
        self.at = end;
        Some(piece)
    }

    fn u8(&mut self) -> Option<u8> {
        self.take(1).map(|b| b[0])
    }

    fn u32(&mut self) -> Option<u32> {
        self.take(4)
            .map(|b| u32::from_le_bytes(b.try_into().expect("4 bytes")))
    }

    fn u64(&mut self) -> Option<u64> {
        self.take(8)
            .map(|b| u64::from_le_bytes(b.try_into().expect("8 bytes")))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1284-1322
// A token holds no space, tab, line feed or carriage return, so `TruncateSpaces` (which
// trims `isspace`: also form feed and vertical tab) can empty it only when it is made of
// those two characters; such a token is reported and dropped.
/// The tokens of a text list.
pub(super) fn read_memory_si_list(bytes: &[u8], warnings: &mut Vec<String>) -> Vec<Vec<u8>> {
    let mut sis = Vec::with_capacity(bytes.len() / 7);
    let mut p = 0;
    while p < bytes.len() {
        while p < bytes.len() && matches!(bytes[p], b'>' | b' ' | b'\t' | b'\n' | b'\r') {
            p += 1;
        }
        if p < bytes.len() && bytes[p] == b'#' {
            while p < bytes.len() && bytes[p] != b'\n' {
                p += 1;
            }
            continue;
        }
        let head = p;
        while p < bytes.len() && !matches!(bytes[p], b' ' | b'\t' | b'\n' | b'\r') {
            p += 1;
        }
        if p > head {
            let acc = &bytes[head..p];
            let str_id = crate::ncbi_str::truncate_spaces(acc);
            if !str_id.is_empty() {
                sis.push(str_id.to_vec());
            } else {
                warnings.push(format!(
                    "WARNING:  {} is not a valid seqid string.",
                    String::from_utf8_lossy(acc)
                ));
            }
        }
    }
    sis
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1476-1491
// A binary list of no ids (`GetSeqidlist` returns 0) is read again as text: its header
// bytes become tokens. The header of a text file is the default (`is_v4`).
/// Reads a seqid list file's bytes: a binary list, else a text list. `Err`: NCBI's
/// `CSeqDBException` text (`BLAST Database error: ...`, exit status 2). The lines NCBI
/// writes to `cerr` for tokens it drops are pushed to `warnings`.
pub fn read_si_list(bytes: &[u8], warnings: &mut Vec<String>) -> Result<SeqidListFile, String> {
    if bytes.is_empty() {
        return Err("Failed to map seqidlist file ".to_string());
    }
    if let Some((info, ids)) = read_binary(bytes)? {
        if info.num_ids > 0 {
            return Ok(SeqidListFile { ids, info });
        }
    }
    Ok(SeqidListFile {
        ids: read_memory_si_list(bytes, warnings),
        info: SeqidListInfo::default(),
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:148-181
/// The parse of a token failed: NCBI posts the exception's message, then throws
/// `Invalid seq id: <token>`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InvalidSeqId {
    pub token: Vec<u8>,
    pub error: SeqIdError,
}

impl InvalidSeqId {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:165-167
    /// The `CSeqDBException` text (`BLAST Database error: Invalid seq id: <token>`), with the
    /// token's bytes (probe R3 of line LIB-B: a word that is not UTF-8).
    pub fn message(&self) -> Vec<u8> {
        [b"Invalid seq id: ".as_slice(), &self.token].concat()
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2532-2541
/// The ids of a text list as a version 5 database looks them up (`s_ProcessSeqIDsForV5`):
/// sorted and unique. An empty result is NCBI's `Empty seqid list` warning (the caller
/// writes it). `warn` gets the warnings the parse of each token posts, in order (a FASTA-style
/// id of several parts: one for each non-empty part after the first id), the gi ids' too.
pub fn process_for_v5(
    ids: &[Vec<u8>],
    warn: &mut dyn FnMut(SeqIdWarning),
) -> Result<Vec<Vec<u8>>, InvalidSeqId> {
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(ids.len());
    for token in ids {
        let id = seq_id_from_text(token, LIST_ID_FLAGS, warn).map_err(|error| InvalidSeqId {
            token: token.clone(),
            error,
        })?;
        let mut text = Vec::new();
        match &id {
            ParsedSeqId::Gi(_) => continue,
            ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => label::fasta(&id, &mut text),
            _ => label::content(&id, true, &mut text),
        }
        out.push(text);
    }
    out.sort();
    out.dedup();
    Ok(out)
}

/// `fParse_RawText | fParse_AnyLocal | fParse_PartialOK`, the flags of a list id's `CSeq_id`.
const LIST_ID_FLAGS: u32 = RAW_TEXT | ANY_LOCAL | PARTIAL_OK;

/// `CSeq_id(text, fParse_RawText | fParse_AnyLocal | fParse_PartialOK)`, the parse of a
/// list id where nothing is posted: the lookup files' ids and the processed ids, one part
/// each ([`process_for_v5`] posts the warnings of the list's tokens).
pub fn parse_list_id(text: &[u8]) -> Result<ParsedSeqId, SeqIdError> {
    seq_id_from_text(text, LIST_ID_FLAGS, &mut |_| {})
}

#[cfg(test)]
#[path = "seqidlist_tests.rs"]
mod tests;
