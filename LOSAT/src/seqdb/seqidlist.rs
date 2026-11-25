//! A seqid list file (line DI-8b: `-seqidlist`, `-negative_seqidlist`), read as NCBI's
//! `CSeqDBFileGiList` with `fname` and `eSiList` reads it (`SeqDB_ReadSiList`): a binary list made
//! by `blastdb_aliastool -seqid_file_in` (its first byte is NUL; `CSeqidlistRead`), else text
//! (`SeqDB_ReadMemorySiList`: white-space separated tokens, a leading `>` skipped, `#` to
//! the end of the line a comment). The ids of a text list are then parsed and rewritten
//! for a version 5 database as `s_ProcessSeqIDsForV5` rewrites them: `gi|N` dropped, `pir`
//! and `prf` ids in FASTA form, every other id as `CSeq_id::GetSeqIdString` with true (the
//! upper-cased accession with its version when given), sorted and unique.

use super::{SearchPath, SeqDbError};
use crate::ncbi_exception::{NcbiException, ThrowSite};
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
    /// A text list (`is_v4` is true, the constructor's default); `false` for a binary list.
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
    let Some((info, mut cursor)) = read_header(bytes)? else {
        return Ok(None);
    };
    let mut ids = Vec::with_capacity(usize::try_from(info.num_ids).unwrap_or(0).min(1 << 20));
    while cursor.at < bytes.len() && (ids.len() as u64) < info.num_ids {
        let id_len = cursor.u8().ok_or(INVALID)?;
        let len = if id_len == LONG_ID {
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

/// NCBI's message of a binary list whose header is damaged.
const INVALID: &str = "Invalid seqidlist file";

/// The length byte of an id of 255 bytes or more, which a 4-byte length follows (`byte_max`).
const LONG_ID: u8 = 0xFF;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:58-85
/// The header of a binary list (`CSeqidlistRead`'s constructor) and the position of its first id:
/// `Ok(None)` for a file that is not binary, `Err` with NCBI's message for a damaged header.
fn read_header(bytes: &[u8]) -> Result<Option<(SeqidListInfo, Cursor<'_>)>, String> {
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
    Ok(Some((info, cursor)))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5748-5755
// NCBI 2.17.0 (the conda build) writes `Error: NCBI C++ Exception:\n    T0
// "/opt/conda/conda-bld/blast_1787399212531/work/c++/src/corelib/ncbifile.cpp", line 5754:
// Error: (CFileException::eMemoryMap) ncbi::CMemoryFileMap::CMemoryFileMap() - To be memory
// mapped the file must exist: ''\n\n` (dropin-survey-20261029/DI8 probe S22, and the official
// Linux x64 build in probe I04 of row RG-D2a); LOSAT names the build tree of the official Linux
// x64 release (`ncbi_exception::NCBI_BUILD_ROOT`), as the Seq-id parser's exceptions do.
const MEMORY_MAP_MISSING: ThrowSite = ThrowSite {
    file: "src/corelib/ncbifile.cpp",
    line: 5754,
    kind: "CFileException::eMemoryMap",
    function: Some("ncbi::CMemoryFileMap::CMemoryFileMap"),
};

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5748-5755
// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5886-5947
// Constructing `CMemoryFile` throws `CFileException::eMemoryMap` at line 5754 when `path` is not
// an existing regular file (a missing file or a directory), and at line 5945 / 5650 when opening
// or mapping fails; an empty regular file leaves the mapped pointer null and length 0.
/// Reads `path` as NCBI's `CMemoryFile` maps it: empty regular file as `Vec::new()`, missing path
/// or directory as `CFileException::eMemoryMap` (`MEMORY_MAP_MISSING`), and open/map failure as
/// `memory_file_exception` ([`SeqDbError::Ncbi`]), else the reader's own I/O error.
pub(super) fn read_memory_file(path: &[u8]) -> Result<Vec<u8>, SeqDbError> {
    let file_path = crate::utils::os_bytes::path_of_bytes(path);
    if !super::alias::exists(&file_path) {
        return Err(SeqDbError::Ncbi(NcbiException::new(
            MEMORY_MAP_MISSING,
            super::oid_list::not_a_file(path),
        )));
    }
    let map_err = |error: SeqDbError| match error.memory_file_exception() {
        Some(exception) => SeqDbError::Ncbi(exception),
        None => error,
    };
    let file = super::file::DbFile::open(&file_path).map_err(map_err)?;
    if file.is_empty() {
        return Ok(Vec::new());
    }
    file.read_all().map_err(map_err)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:453-462
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:122-130
// The file is looked for on the database search path (`SeqDB_ResolveDbPath`). A file that is not
// there resolves to the empty name, whose memory map throws a `CFileException` (not a database
// error): a search reports it as `Error: <what>` with exit status 255 (`BLAST_UNKNOWN_ERROR`),
// `blastdb_aliastool -seqid_file_info` too. A file that is there but does not open (no permission)
// throws `CMemoryFile`'s exception with its path and the system's text, reported the same way.
/// The bytes of the seqid list file `name` as NCBI's readers map them (`CMemoryFile` of the path
/// that `SeqDB_ResolveDbPath` gives for `name`): the first file of that name on the search path
/// `path`. `Err`: NCBI's exception of the memory map ([`SeqDbError::Ncbi`]: of the empty name when
/// there is no such file, or of the file that does not open), else the reader's error of a file
/// that does not read without a system text for NCBI's message.
pub(crate) fn read_resolved(name: &[u8], path: &SearchPath) -> Result<Vec<u8>, SeqDbError> {
    let Some(file) = path.resolve_file(name) else {
        return Err(SeqDbError::Ncbi(NcbiException::new(
            MEMORY_MAP_MISSING,
            super::oid_list::not_a_file(b""),
        )));
    };
    read_memory_file(&file)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/seqidlist_writer.cpp:145-183
// The lengths of the title, the dates and the volume names are NCBI's `Uint4` and `char` of the
// strings' sizes (a date is never 256 bytes long); the database part is written when the list
// names volumes, as `CSeqidlistRead` reads it when the total length is not 0.
/// The bytes of a binary list (`WriteBlastSeqidlistFile`, the format that [`read_si_list`]
/// reads) of `info` and `ids`: the NUL byte, the file's size, the number of ids, the title, the
/// creation date, the database's total length and, with `info.db_vol_names`, its date and
/// volume names, then each id after its length (a byte, or 0xFF and 4 bytes from 255 on). The
/// size and the number of ids are those of the bytes made (`info.file_size` and `info.num_ids`
/// are not read).
pub fn write_binary(info: &SeqidListInfo, ids: &[Vec<u8>]) -> Vec<u8> {
    let mut out = vec![0u8];
    out.extend_from_slice(&0u64.to_le_bytes());
    out.extend_from_slice(&(ids.len() as u64).to_le_bytes());
    out.extend_from_slice(&(info.title.len() as u32).to_le_bytes());
    out.extend_from_slice(&info.title);
    out.push(info.create_date.len() as u8);
    out.extend_from_slice(&info.create_date);
    out.extend_from_slice(&info.db_vol_length.to_le_bytes());
    if !info.db_vol_names.is_empty() {
        out.push(info.db_create_date.len() as u8);
        out.extend_from_slice(&info.db_create_date);
        out.extend_from_slice(&(info.db_vol_names.len() as u32).to_le_bytes());
        out.extend_from_slice(&info.db_vol_names);
    }
    for id in ids {
        match u8::try_from(id.len()).ok().filter(|&len| len != LONG_ID) {
            Some(len) => out.push(len),
            None => {
                out.push(LONG_ID);
                out.extend_from_slice(&(id.len() as u32).to_le_bytes());
            }
        }
        out.extend_from_slice(id);
    }
    let size = out.len() as u64;
    out[1..9].copy_from_slice(&size.to_le_bytes());
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:122-155
// The file is mapped whole (`CMemoryFile`); only the header is read, so a list whose ids are
// cut short still prints.
/// What `blastdb_aliastool -seqid_file_info` prints of a list's bytes `bytes`
/// (`PrintSeqidlistInfo`): the header of a binary list of at least one id, else NCBI's line
/// that the file is not of version 5. `Err`: NCBI's `CSeqDBException` text for an empty file or
/// a damaged header.
pub fn info_text(bytes: &[u8]) -> Result<Vec<u8>, String> {
    if bytes.is_empty() {
        return Err("Failed to map seqidlist file ".to_string());
    }
    let mut out = Vec::new();
    // NCBI's test is of the count as an `int` (`GetSeqidlistInfo`'s `static_cast<int>`): a count
    // of 2^31 or more whose low 32 bits are not a positive `int` prints the line of version 4.
    match read_header(bytes)? {
        Some((info, _)) if (info.num_ids as i32) > 0 => {
            out.extend_from_slice(format!("Num of Ids: {}\nTitle: ", info.num_ids).as_bytes());
            out.extend_from_slice(&info.title);
            out.extend_from_slice(b"\nCreate Date: ");
            out.extend_from_slice(&info.create_date);
            out.push(b'\n');
            if info.db_vol_length > 0 {
                out.extend_from_slice(
                    format!(
                        "DB Info: \n\tTotal Vol Length: {}\n\tDB Create Date: ",
                        info.db_vol_length
                    )
                    .as_bytes(),
                );
                out.extend_from_slice(&info.db_create_date);
                out.extend_from_slice(b"\n\tDB Vols: ");
                // `NStr::Split` of the names at spaces, with its default flags: empty pieces
                // between two spaces are kept, and empty names are no piece.
                let vols = (!info.db_vol_names.is_empty())
                    .then(|| info.db_vol_names.split(|&byte| byte == b' '))
                    .into_iter()
                    .flatten();
                for vol in vols {
                    out.extend_from_slice(b"\n\t\t");
                    out.extend_from_slice(vol);
                }
            }
        }
        _ => out.extend_from_slice(b"Seqidlist file is not in blast db version 5 format"),
    }
    out.push(b'\n');
    Ok(out)
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

/// `fParse_RawText`, `fParse_AnyLocal` and `fParse_PartialOK`, the flags of a list id's `CSeq_id`.
const LIST_ID_FLAGS: u32 = RAW_TEXT | ANY_LOCAL | PARTIAL_OK;

/// `CSeq_id` of `text` with the flags `fParse_RawText`, `fParse_AnyLocal` and `fParse_PartialOK`,
/// the parse of a list id where nothing is posted: the lookup files' ids and the processed ids, one
/// part each ([`process_for_v5`] posts the warnings of the list's tokens).
pub fn parse_list_id(text: &[u8]) -> Result<ParsedSeqId, SeqIdError> {
    seq_id_from_text(text, LIST_ID_FLAGS, &mut |_| {})
}

#[cfg(test)]
#[path = "seqidlist_tests.rs"]
mod tests;
