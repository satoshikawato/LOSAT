//! The OID list file of an alias file's `OIDLIST` key (`.msk`, and makeblastdb's `.pxm` / `.nxm`
//! of `-oid_masks`, `seqdb::oid_mask`): the big-endian `Uint4` index of its last OID, then one bit
//! per OID, the first OID of each byte in its most significant bit, set for an OID the list keeps
//! (`writedb::oid_list` writes it). The file is read relative to each volume it is applied to: bit
//! `i` is OID `i` of that volume, whichever volume it is (`alias::alias_oid_mask`).

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

/// Why an OID list file cannot be applied: NCBI's error, or LOSAT's refusal where NCBI crashes.
/// It is kept with the database where NCBI meets it later (`SeqDb::oid_list_error`), so it is
/// `Clone` and makes its `SeqDbError` on demand (`OidListError::to_error`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum OidListError {
    /// The path is not a regular file (missing, a directory, the directory of an empty value):
    /// NCBI's `eMemoryMappingFailure` (`not_a_file`).
    NotAFile(String),
    /// The file holds no byte: NCBI maps nothing and reads the header through a null pointer
    /// (SIGSEGV, probe G5 of line RG D1b).
    Empty(PathBuf),
    /// The file holds 1 to 4 bytes, no byte of bits: NCBI's `SEQDB_FILE_ASSERT(start < end)` of
    /// `GetFileDataPtr(lease, 4, file_length)` fails (probes G3, G4).
    Short(PathBuf, u64),
    /// The bits NCBI copies for a volume reach past the last memory page of the file's mapping,
    /// where NCBI reads memory that is not the file's.
    PastPage {
        path: PathBuf,
        /// The bytes NCBI reads, from offset 4.
        read: u64,
        len: u64,
    },
    /// The file cannot be opened or read (`DbFile`): `SeqDbError::Open` with the system's code,
    /// else `SeqDbError::Io`.
    Read {
        path: PathBuf,
        open: bool,
        code: Option<i32>,
        text: String,
    },
}

impl OidListError {
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
            OidListError::PastPage { path, read, len } => SeqDbError::NcbiCrashes {
                path: path.clone(),
                message: format!(
                    "NCBI reads {read} bytes of bits of the OID list file of an alias file's \
                     OIDLIST, past the memory page after its {len} bytes"
                ),
            },
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
// `GetFileDataPtr(lease, 4, file_length)` asserts `4 < file_length`.
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
    // `num_oids` is the header plus one and `bitend` lies `((num_oids + 31) / 32) * 4` bytes after
    // the bits, in `Uint4` arithmetic: a header of 0xFFFFFFFF leaves no byte (probe H3), and the
    // bits past the header's last OID count up to their 32-bit word (H1). `CSeqDB_BitSet` copies
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
        for at in 0..count as usize {
            let byte = stored.get(at).copied().unwrap_or(0);
            if byte == 0 {
                continue;
            }
            for bit in 0..8 {
                if byte & (0x80 >> bit) != 0 {
                    let oid = at * 8 + bit;
                    if oid < oids as usize {
                        words[oid / 64] |= 1u64 << (oid % 64);
                    }
                }
            }
        }
        Ok(words)
    }
}

#[cfg(test)]
#[path = "oid_list_tests.rs"]
mod tests;
