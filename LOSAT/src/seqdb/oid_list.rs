//! The OID list file of an alias file's `OIDLIST` key (`.msk`, and makeblastdb's `.pxm` / `.nxm`
//! of `-oid_masks`, `seqdb::oid_mask`), the one owner of its layout: the big-endian `Uint4` index
//! of its last OID (makeblastdb writes it as an `Int4`, `writedb::oid_list`), then one bit per OID,
//! the first OID of each byte in its most significant bit, set for an OID the list keeps. The file
//! is read relative to each volume it is applied to: bit `i` is OID `i` of that volume, whichever
//! volume it is (`alias::alias_oids`).

use std::path::{Path, PathBuf};

use super::file::{mapped_length, DbFile};
use super::SeqDbError;

// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5748-5756
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:699-722
// `CMemoryFileMap` throws this text for a path that is not an existing file (a directory has the
// length -1, ncbifile.cpp:5886-5890), and `CSeqDBFileMemMap::Init` makes it the message of a
// `CSeqDBException` of the code `eMemoryMappingFailure` (`SeqDbError::MemoryMap`).
/// The message of NCBI's error for an OID list file that is not a file.
fn not_a_file(path: &str) -> Vec<u8> {
    [
        b"To be memory mapped the file must exist: '".as_slice(),
        path.as_bytes(),
        b"'",
    ]
    .concat()
}

/// Why NCBI's OID list of the alias files cannot be built: the error of an OID list file, or
/// LOSAT's refusal where NCBI crashes.
/// It is kept with the database where NCBI meets it later (`SeqDb::oid_list_error`), so it is
/// `Clone` and makes its `SeqDbError` on demand (`OidListError::to_error`).
#[derive(Debug, Clone)]
pub(super) enum OidListError {
    /// The path is not a regular file (missing, a directory, the directory of an empty value):
    /// NCBI's `eMemoryMappingFailure` (`not_a_file`).
    NotAFile(String),
    /// The file holds no byte: NCBI maps nothing and reads the header through a null pointer
    /// (SIGSEGV, probe G5 of line RG D1b).
    Empty(PathBuf),
    /// The file holds 1 to 4 bytes, no byte of bits: NCBI's assertion that the file is longer
    /// than its 4-byte header fails (probes G3, G4).
    Short(PathBuf, u64),
    /// The bits NCBI copies for a volume reach past the last memory page of the file's mapping,
    /// where NCBI reads memory that is not the file's.
    PastPage {
        path: PathBuf,
        /// The bytes NCBI reads, from offset 4.
        read: u64,
        len: u64,
    },
    /// An alias file's `FIRST_OID`-`LAST_OID` range lies outside the bit array of the OIDs of the
    /// alias files below it, where NCBI clears bits outside that array (`alias::OidSet`): a
    /// crash, or a write into memory that is not the array's. The path is the alias file's.
    RangeOutsideBits(PathBuf),
    /// The file cannot be opened or read (`DbFile`): `SeqDbError::Open` with the system's code,
    /// else `SeqDbError::Io`.
    Read {
        path: PathBuf,
        open: bool,
        code: Option<i32>,
        text: String,
    },
    /// The error of an id list of the alias files (`alias::ListSet`): its file's
    /// `CSeqDBException` text, the error of a volume's ISAM that translates it, LOSAT's refusal
    /// (copied each time it is met, `SeqDbError::duplicate`).
    Reader(std::sync::Arc<SeqDbError>),
}

impl OidListError {
    /// `error` of an id list of the alias files, kept ([`OidListError::Reader`]).
    pub(super) fn reader(error: SeqDbError) -> Self {
        OidListError::Reader(std::sync::Arc::new(error))
    }

    /// The reader error of this failure: NCBI's text where NCBI has one (`SeqDbError::MemoryMap`,
    /// `SeqDbError::Open`), else LOSAT's (`SeqDbError::Integrity` for NCBI's failed assertion,
    /// `SeqDbError::NcbiCrashes` where NCBI crashes).
    pub(super) fn to_error(&self) -> SeqDbError {
        match self {
            OidListError::NotAFile(path) => SeqDbError::MemoryMap(not_a_file(path)),
            OidListError::Empty(path) => SeqDbError::NcbiCrashes {
                path: path.clone(),
                message: "the OID list file of an alias file's OIDLIST is empty".to_string(),
            },
            OidListError::Short(path, len) => SeqDbError::Integrity {
                path: path.clone(),
                message: format!(
                    "the OID list file of an alias file's OIDLIST holds {len} bytes, no bit \
                     after its 4-byte header"
                ),
            },
            OidListError::RangeOutsideBits(path) => SeqDbError::NcbiCrashes {
                path: path.clone(),
                message: "its FIRST_OID-LAST_OID range lies outside the OIDs that the alias files \
                          below it keep, where NCBI clears bits outside its array of them, over \
                          memory that is not the array's"
                    .to_string(),
            },
            OidListError::PastPage { path, read, len } => SeqDbError::NcbiCrashes {
                path: path.clone(),
                message: format!(
                    "NCBI reads {read} bytes of bits of the OID list file of an alias file's \
                     OIDLIST, past the memory page after its {len} bytes"
                ),
            },
            OidListError::Reader(kept) => kept.duplicate(),
            OidListError::Read {
                path,
                open,
                code,
                text,
            } => {
                let source = match code {
                    Some(code) => std::io::Error::from_raw_os_error(*code),
                    None => std::io::Error::other(text.clone()),
                };
                if *open {
                    SeqDbError::Open {
                        path: path.clone(),
                        source,
                    }
                } else {
                    SeqDbError::Io {
                        path: path.clone(),
                        source,
                    }
                }
            }
        }
    }

    /// The failure of a `DbFile` open or read.
    fn read(error: SeqDbError, path: &Path) -> Self {
        let (open, source) = match error {
            SeqDbError::Open { source, .. } => (true, Some(source)),
            SeqDbError::Io { source, .. } => (false, Some(source)),
            _ => (false, None),
        };
        OidListError::Read {
            path: path.to_path_buf(),
            open,
            code: source.as_ref().and_then(std::io::Error::raw_os_error),
            text: source.map_or_else(|| "read failed".to_string(), |s| s.to_string()),
        }
    }
}

/// One OID list file, read whole (a few bytes a thousand OIDs).
#[derive(Debug)]
pub(super) struct OidListFile {
    path: PathBuf,
    /// The bytes of the file.
    bytes: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:464-505
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:71-79
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:123-133
// `x_GetOidMask` maps the file (`not_a_file`), reads the header with `ReadSwapped` (big-endian;
// a file of 1 to 3 bytes is read past its end, inside the mapped page, before the check), then
// asserts that the file holds a byte after the 4-byte header.
/// Reads the OID list file `path` (the path of the alias file's `OIDLIST`, as
/// `alias::compute_masks` resolved it).
pub(super) fn read(path: &str) -> Result<OidListFile, OidListError> {
    let at = Path::new(path);
    if !super::alias::exists(at) {
        return Err(OidListError::NotAFile(path.to_string()));
    }
    let file = DbFile::open(at).map_err(|error| OidListError::read(error, at))?;
    if file.is_empty() {
        return Err(OidListError::Empty(at.to_path_buf()));
    }
    if file.len() <= 4 {
        return Err(OidListError::Short(at.to_path_buf(), file.len()));
    }
    let bytes = file
        .read_all()
        .map_err(|error| OidListError::read(error, at))?;
    Ok(OidListFile {
        path: at.to_path_buf(),
        bytes,
    })
}

impl OidListFile {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:37-58
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.hpp:248-251
    // The number of OIDs is the header plus one, and the end of the bits lies that number rounded
    // up to whole 32-bit words, in bytes, after their start, in `Uint4` arithmetic: a header of
    // 0xFFFFFFFF leaves no byte (probe H3), and the bits past the header's last OID count up to
    // their 32-bit word (H1). `CSeqDB_BitSet` copies
    // the bytes of the volume's OIDs, `ceil(oids / 8)`, or fewer when `bitend` comes first. The bytes
    // past the end of the file are the zeros of its last mapped page (H2, H7, `file::mapped_length`);
    // a read past that page reads memory that is not the file's, which LOSAT refuses
    // (`OidListError::PastPage`). The bits of the last byte past the volume's OIDs are outside the
    // bit set (`m_End`), so no OID is set there.
    /// The OIDs of a volume of `oids` OIDs that the file keeps: bit `oid % 64` of word `oid / 64`
    /// is set for a kept OID `oid` of the volume.
    pub(super) fn volume_words(&self, oids: u32) -> Result<Vec<u64>, OidListError> {
        let header = u32::from_be_bytes(self.bytes[..4].try_into().expect("4 header bytes"));
        let num_oids = header.wrapping_add(1);
        let available = u64::from(num_oids.wrapping_add(31) / 32) * 4;
        let count = u64::from(oids).div_ceil(8).min(available);
        let len = self.bytes.len() as u64;
        if 4 + count > mapped_length(len) {
            return Err(OidListError::PastPage {
                path: self.path.clone(),
                read: count,
                len,
            });
        }
        let mut words = vec![0u64; (oids as usize).div_ceil(64)];
        let stored = &self.bytes[4..];
        let count = (count as usize).min(stored.len());
        // The first OID of a byte is its most significant bit, of a word its least significant.
        for (at, &byte) in stored[..count].iter().enumerate() {
            words[at / 8] |= u64::from(byte.reverse_bits()) << ((at % 8) * 8);
        }
        let tail = oids % 64;
        if tail > 0 {
            if let Some(last) = words.last_mut() {
                *last &= (1u64 << tail) - 1;
            }
        }
        Ok(words)
    }
}

#[cfg(test)]
#[path = "oid_list_tests.rs"]
mod tests;
