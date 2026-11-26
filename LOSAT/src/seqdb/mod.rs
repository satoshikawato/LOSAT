//! Reader of NCBI BLAST databases (the volume layout that makeblastdb writes,
//! format versions 4 and 5).
//!
//! NCBI's reader is `c++/src/objtools/blast/seqdb_reader/` (CSeqDB). This
//! module reads the same files with positioned reads instead of memory
//! mapping, and hands out sequences in blocks: contiguous ranges of ordinal
//! ids whose bytes are read with one call.
//!
//! Read: the index file (`.pin` / `.nin`), the sequence file (`.psq`:
//! NCBIstdaa; `.nsq`: 2-bit bases followed by ambiguity data), the header
//! file (`.phr` / `.nhr`: title and `gnl|BL_ORD_ID|n` with `deflines`, every
//! member and Seq-id with `defline_set`), and alias files (`.pal` / `.nal`:
//! `TITLE` and `DBLIST`).
//!
//! Not read: the LMDB file and the id and taxonomy files of version 5, the
//! ISAM files of `-parse_seqids`, column files (mask data and user columns: the
//! mask reader of `losat blastdbcmd` opens them on first use, see [`mask`]; LOSAT's
//! searches never do, as NCBI's do not without a database mask, whose options
//! LOSAT refuses), group alias files, and alias keys that filter
//! sequences or replace the totals. A database
//! that needs one of these to give NCBI's result is rejected with
//! [`SeqDbError::Unsupported`]; the version 5 files that only serve lookups
//! by id or taxonomy are left unread. The LMDB file of a version 5 database
//! must exist, as for NCBI, but its contents are not checked.
//!
//! EXPERIMENT (LOSAT_X_DBIDXREAD, in the strict set): an index file stays open
//! after it is checked, and the ranges of its offset arrays that a block, a
//! header or a planned block needs are read by position, so the arrays take
//! no memory (NCBI maps the file). Without the switch the arrays are held
//! (4 bytes per sequence each).
//!
//! A search program opens the database named by its `-db` with
//! [`SeqDb::open_in`] and a [`SearchPath`]: the name is looked for as NCBI
//! BLAST+ looks for it (current directory, `BLASTDB`, the registry value
//! `[BLAST] BLASTDB`). [`SeqDbError::ncbi_message`] tells NCBI's own errors,
//! which a search program reports as `BLAST Database error: <text>` with exit
//! status 2, from LOSAT's refusals.

pub mod alias;
pub mod column;
pub mod date;
pub mod file;
pub mod filter;
pub(crate) mod find;
pub mod gi_mask;
pub mod header;
pub(crate) mod id_list;
pub mod index;
pub mod isam;
pub mod isam_num;
pub mod lmdb_read;
pub mod lookup;
pub mod mask;
pub(crate) mod oid_list;
pub(crate) mod oid_mask;
pub(crate) mod rps;
pub mod seqidlist;
pub mod tax_files;
pub mod tax_sqlite;
pub mod taxdb;
pub mod volume;

pub use filter::{DbList, GiList, ListWarning, PigList, SeqidList, TaxidList};
pub(crate) use filter::{DbSequence, Oid, OidMap};

#[cfg(test)]
mod damaged_files_tests;
#[cfg(test)]
mod mask_tests;
#[cfg(test)]
mod seqid_taxid_damaged_tests;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod x_dbidxread_tests;

use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::sync::OnceLock;

use crate::ncbi_exception::ThrowSite;
use crate::seqid::model::ObjectId;
use crate::seqid::{DeflineSet, ParsedSeqId};
use crate::utils::os_bytes::{path_from_bytes, path_of_bytes};

pub use alias::SearchPath;
pub use column::ColumnBlock;
pub use header::Defline;
pub use mask::{BlockMasks, MaskAlgorithm, MaskAlgorithms, MaskBlock, MaskRanges, NameMaskError};
pub(crate) use volume::ncbi_protein_length;
pub use volume::{AmbiguityRun, NucleotideAmbiguity, SequenceBlock, Volume, VolumeHeader};
#[cfg(test)]
pub use volume::{BLASTNA_SENTINEL, FENCE_SENTRY};

/// EXPERIMENT (LOSAT_X_BNDBSMALL): offsets per read of `SeqDb::x_next_blocks` (64 KiB).
#[cfg(any(
    test,
    all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )
))]
const X_LIST_ENTRIES: usize = 1 << 14;

/// EXPERIMENT (LOSAT_X_DBIDXREAD): how the volumes read their index files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum XDbIdxRead {
    /// The whole file is read and the offset arrays are held (`IndexFile::parse`).
    Off,
    /// The arrays are read from the open file where they are needed (`IndexFile::open`).
    On,
}

// No NCBI counterpart: reads the LOSAT_X_DBIDXREAD switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_DBIDXREAD): a strict-set switch, on unless it is `0`.
pub(crate) fn x_dbidxread() -> XDbIdxRead {
    static MODE: OnceLock<XDbIdxRead> = OnceLock::new();
    *MODE.get_or_init(|| {
        if crate::utils::xswitch::on("LOSAT_X_DBIDXREAD") {
            XDbIdxRead::On
        } else {
            XDbIdxRead::Off
        }
    })
}

// No NCBI counterpart: reads the LOSAT_X_BPDBKEEP switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): a strict-set switch, on unless it is `0`
/// (`algorithm/blastp/blast_engine/database/x_db_keep.rs`).
pub(crate) fn x_bp_db_keep() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDBKEEP"))
}

/// EXPERIMENT (LOSAT_X_DBIDXREAD): the most sequences of a range of
/// `SeqDb::header_ranges` that the checks of every title read at once.
pub const HEADER_RANGE_SEQS: u32 = 1024;

/// EXPERIMENT (LOSAT_X_DBIDXREAD): the most bytes of the header file that
/// `SeqDb::visit_deflines` reads at once (one longer record alone).
const HEADER_PIECE_BYTES: u64 = 1 << 20;

/// EXPERIMENT (LOSAT_X_BPDBKEEP): `HEADER_PIECE_BYTES` with the switch (64 KiB): the reads stay
/// below glibc's first mmap threshold (128 KiB), so that freeing them does not raise it (a range
/// of 1,024 UniRef50 headers is about 190 KiB). The size changes nothing but the number of reads.
const X_KEEP_HEADER_PIECE_BYTES: u64 = 1 << 16;

/// EXPERIMENT (LOSAT_X_BPDBKEEP): the header piece size (`X_KEEP_HEADER_PIECE_BYTES` with the
/// switch).
fn header_piece_bytes() -> u64 {
    if x_bp_db_keep() {
        X_KEEP_HEADER_PIECE_BYTES
    } else {
        HEADER_PIECE_BYTES
    }
}

/// Sequence type of a database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    Protein,
    Nucleotide,
}

impl SeqType {
    /// First letter of the file extensions: `p` or `n`.
    pub fn letter(self) -> &'static str {
        match self {
            SeqType::Protein => "p",
            SeqType::Nucleotide => "n",
        }
    }

    /// The word NCBI's messages use.
    pub fn name(self) -> &'static str {
        match self {
            SeqType::Protein => "protein",
            SeqType::Nucleotide => "nucleotide",
        }
    }
}

/// Error of the database reader.
///
/// [`SeqDbError::ncbi_message`] says whether an error is one that NCBI BLAST+
/// raises too, and with which text.
#[derive(Debug)]
pub enum SeqDbError {
    /// A database, volume or file was not found, or no database was named.
    /// The message is the text of NCBI's exception, as bytes: a name it holds keeps the bytes
    /// of NCBI's `std::string` (a seqid list's word, an LMDB file's name that is not UTF-8:
    /// probes R3 and R8 of line LIB-B).
    NotFound(Vec<u8>),
    /// A file that exists could not be opened.
    Open {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The operating system failed to tell the length of a file that it opened.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The operating system failed to read the first byte of a file of `len` bytes where it was
    /// opened (`file::DbFile::open`): what NCBI meets where it maps the file
    /// ([`SeqDbError::ncbi_message`]).
    Map {
        path: PathBuf,
        len: u64,
        source: std::io::Error,
    },
    /// The operating system failed to read bytes of a file after its first byte was read where
    /// it was opened (a file that a disk or the network stopped giving, a file that became
    /// shorter): where NCBI reads mapped bytes that the system no longer gives and ends by
    /// SIGBUS. LOSAT's own error ("cannot be read"; line WQA D18, 2026-11-21).
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
    /// A file is not a valid database file.
    Format { path: PathBuf, message: String },
    /// A list of volumes of both version 4 and version 5 (`CSeqDBLMDBSet`): NCBI's text
    /// [`MIXED_VERSIONS`], from one of two throws, which the volume `path`, the first whose
    /// version is not the first volume's, chooses: a version 4 volume after version 5 ones
    /// (`version4_after_version5`), or a version 5 volume after version 4 ones.
    MixedVersions {
        path: PathBuf,
        version4_after_version5: bool,
    },
    /// A file fails one of the integrity checks of NCBI's reader (`SEQDB_FILE_ASSERT`, a
    /// `CSeqDBException` whose text names a line of NCBI's source): LOSAT's words, and NCBI's
    /// text where NCBI's search replaces it ([`SeqDbError::ncbi_search_message`]).
    Integrity { path: PathBuf, message: String },
    /// A valid database that needs something LOSAT does not read.
    Unsupported { path: PathBuf, message: String },
    /// A file that NCBI BLAST+ crashes on where it reads it (an ISAM whose header does not fit its
    /// files, as an empty `.psi` or `.psd`: NCBI's constructor reads its index bounds whatever
    /// `x_InitSearch` returned, `isam_num::NumericIsam::open_volume`; probe I of line BDC, review 4
    /// of line DBID loadord), which LOSAT refuses as an `Unsupported` database: the text of one
    /// whose message is `message` followed by [`NCBI_CRASHES_ON_IT`]
    /// ([`SeqDbError::ncbi_crashes`]).
    NcbiCrashes { path: PathBuf, message: String },
    /// An exception that NCBI's reader throws and whose report LOSAT does not give (the
    /// `CSeqIdException` of NCBI's simplification of a seqid list's id,
    /// `SeqDB_SimplifyAccession`): LOSAT's refusal wherever a program reports it, in the words of an
    /// `Unsupported` database; a caller that catches every exception of NCBI's (the alias file
    /// writer's catch of any exception, `writedb::alias_file`) takes it for NCBI's failure
    /// ([`SeqDbError::ncbi_throws`]).
    NcbiThrows { path: PathBuf, message: String },
    /// A `CSeqDBException` of the code `eMemoryMappingFailure` whose message is NCBI's (an OID
    /// list file of an alias file's `OIDLIST` that is not a file, `oid_list::OidListError`; the
    /// GI-mask data file of a volume that is not found, `gi_mask::GiMaskFiles::open`): a
    /// search program reports `Database memory map file error - <message>`
    /// ([`SeqDbError::ncbi_message`]), blastdbcmd the message alone
    /// ([`SeqDbError::ncbi_request_message`]).
    MemoryMap(Vec<u8>),
    /// A request outside the database (ordinal id, block range).
    Range(String),
    /// A LOSAT internal error: a caller asked for what it did not read (a sequence outside the
    /// block of masks it read, `mask::MaskBlock`). LOSAT's words, never a database error of
    /// NCBI's ([`SeqDbError::ncbi_message`] is `None`).
    Internal(String),
    /// An exception of NCBI's C++ Toolkit that is not a `CSeqDBException`: the
    /// `CStringException` of an alias file value that does not convert (`alias::walk_totals`,
    /// `alias::compute_masks`). NCBI's programs report it as `Error: <what()>` with exit
    /// status 255 ([`SeqDbError::ncbi_exception`]).
    Ncbi(crate::ncbi_exception::NcbiException),
    /// A file of the taxonomy database (`taxdb.bti`, `taxdb.btd`) that does not open or whose
    /// first byte cannot be read (`taxdb::TaxDb::open`), with the system's text `reason`: NCBI's
    /// taxonomy file object throws `CMemoryFile`'s exception `exception`
    /// ([`SeqDbError::memory_file_exception`]) at every name lookup of the run, and so does every
    /// load of a database sequence with a positive taxid (line WQA D18b, 2026-11-21). Where
    /// nothing catches it, NCBI's programs report it as `Error: <what()>` with exit status 255
    /// ([`SeqDbError::ncbi_exception`]). The exception is boxed: the error stays small.
    TaxdbThrows {
        path: PathBuf,
        reason: String,
        exception: Box<crate::ncbi_exception::NcbiException>,
    },
    /// A header record that fails a check of NCBI's reader that the run's settings turn off,
    /// after which NCBI reads on (`SeqDb::refuse_unchecked_reads`): what LOSAT refuses, a noun
    /// phrase (`ReadVerification::refused`), which the program's refusal names
    /// (`algorithm::common::database::losat_reader_error`).
    Setting(String),
}

impl std::fmt::Display for SeqDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqDbError::NotFound(message) | SeqDbError::MemoryMap(message) => {
                f.write_str(&String::from_utf8_lossy(message))
            }
            SeqDbError::Range(message) | SeqDbError::Internal(message) => f.write_str(message),
            SeqDbError::Setting(refused) => write!(f, "{refused} is not supported by LOSAT"),
            SeqDbError::Open { path, source } | SeqDbError::Io { path, source } => {
                write!(f, "{}: {}", path.display(), source)
            }
            SeqDbError::Map { path, source, .. } | SeqDbError::Read { path, source } => {
                write!(f, "{} cannot be read: {}", path.display(), source)
            }
            SeqDbError::Format { path, message } | SeqDbError::Integrity { path, message } => {
                write!(f, "{}: {}", path.display(), message)
            }
            SeqDbError::MixedVersions { path, .. } => {
                write!(f, "{}: {MIXED_VERSIONS}", path.display())
            }
            SeqDbError::Unsupported { path, message }
            | SeqDbError::NcbiThrows { path, message } => {
                if path.as_os_str().is_empty() {
                    write!(f, "unsupported database: {message}")
                } else {
                    write!(f, "{}: unsupported database: {}", path.display(), message)
                }
            }
            SeqDbError::NcbiCrashes { path, message } => write!(
                f,
                "{}: unsupported database: {} {NCBI_CRASHES_ON_IT}",
                path.display(),
                message
            ),
            SeqDbError::Ncbi(exception) => {
                f.write_str(String::from_utf8_lossy(&exception.what()).trim_end())
            }
            SeqDbError::TaxdbThrows { path, reason, .. } => {
                write!(f, "{} cannot be opened: {reason}", path.display())
            }
        }
    }
}

/// The end of the text of LOSAT's refusal of a database file that NCBI BLAST+ crashes on where it
/// reads it ([`SeqDbError::NcbiCrashes`]).
pub const NCBI_CRASHES_ON_IT: &str = "(NCBI crashes on it)";

impl SeqDbError {
    /// A copy of the error, met again later (`oid_list::OidListError::Reader`): the system's
    /// error of an I/O failure is copied as its code, else as its kind and text.
    pub(super) fn duplicate(&self) -> SeqDbError {
        let io = |source: &std::io::Error| match source.raw_os_error() {
            Some(code) => std::io::Error::from_raw_os_error(code),
            None => std::io::Error::new(source.kind(), source.to_string()),
        };
        match self {
            SeqDbError::NotFound(message) => SeqDbError::NotFound(message.clone()),
            SeqDbError::MemoryMap(message) => SeqDbError::MemoryMap(message.clone()),
            SeqDbError::Open { path, source } => SeqDbError::Open {
                path: path.clone(),
                source: io(source),
            },
            SeqDbError::Io { path, source } => SeqDbError::Io {
                path: path.clone(),
                source: io(source),
            },
            SeqDbError::Map { path, len, source } => SeqDbError::Map {
                path: path.clone(),
                len: *len,
                source: io(source),
            },
            SeqDbError::Read { path, source } => SeqDbError::Read {
                path: path.clone(),
                source: io(source),
            },
            SeqDbError::Format { path, message } => SeqDbError::Format {
                path: path.clone(),
                message: message.clone(),
            },
            SeqDbError::MixedVersions {
                path,
                version4_after_version5,
            } => SeqDbError::MixedVersions {
                path: path.clone(),
                version4_after_version5: *version4_after_version5,
            },
            SeqDbError::Integrity { path, message } => SeqDbError::Integrity {
                path: path.clone(),
                message: message.clone(),
            },
            SeqDbError::Unsupported { path, message } => SeqDbError::Unsupported {
                path: path.clone(),
                message: message.clone(),
            },
            SeqDbError::NcbiCrashes { path, message } => SeqDbError::NcbiCrashes {
                path: path.clone(),
                message: message.clone(),
            },
            SeqDbError::NcbiThrows { path, message } => SeqDbError::NcbiThrows {
                path: path.clone(),
                message: message.clone(),
            },
            SeqDbError::Range(message) => SeqDbError::Range(message.clone()),
            SeqDbError::Internal(message) => SeqDbError::Internal(message.clone()),
            SeqDbError::Setting(refused) => SeqDbError::Setting(refused.clone()),
            SeqDbError::Ncbi(exception) => SeqDbError::Ncbi(exception.clone()),
            SeqDbError::TaxdbThrows {
                path,
                reason,
                exception,
            } => SeqDbError::TaxdbThrows {
                path: path.clone(),
                reason: reason.clone(),
                exception: exception.clone(),
            },
        }
    }

    /// Whether this is the error of the ordinal-id decoders for a record with a Seq-id other
    /// than the `gnl|BL_ORD_ID|n` of the sequence's own number in its volume (a database built
    /// with `-parse_seqids`; `SeqDb::read_header`).
    pub fn is_parsed_seqids(&self) -> bool {
        matches!(self, SeqDbError::Unsupported { message, .. } if message == header::PARSED_SEQIDS)
    }

    /// Whether this is LOSAT's refusal of a file that NCBI BLAST+ crashes on where it reads it
    /// ([`SeqDbError::NcbiCrashes`]; review 5 of line DBID loadord, L4: the variant, not the text,
    /// says it).
    pub fn ncbi_crashes(&self) -> bool {
        matches!(self, SeqDbError::NcbiCrashes { .. })
    }

    /// Whether this is LOSAT's refusal of an exception that NCBI's reader throws
    /// ([`SeqDbError::NcbiThrows`]): the variant, not the text, says it.
    pub fn ncbi_throws(&self) -> bool {
        matches!(self, SeqDbError::NcbiThrows { .. })
    }

    /// Whether this is a file whose bytes the system does not give, in LOSAT's words: a read
    /// that fails after the file's first byte was read ([`SeqDbError::Read`], where NCBI ends by
    /// SIGBUS), or a first byte that cannot be read without NCBI's text for it
    /// ([`SeqDbError::Map`] whose [`SeqDbError::ncbi_message`] is `None`). A program reports it
    /// as LOSAT's error, not as a database that LOSAT does not read (line WQA D18 review M1,
    /// 2026-11-21).
    pub fn is_unreadable(&self) -> bool {
        match self {
            SeqDbError::Read { .. } => true,
            SeqDbError::Map { .. } => self.ncbi_message().is_none(),
            _ => false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
    // `BLAST_UNKNOWN_ERROR` is 255 (blast_app_util.hpp:161); blastdbcmd catches the exception the
    // same way (`CATCH_ALL`). NCBI's guess of a database type moves on to the nucleotide
    // database only after a `CSeqDBException` (seqdb.cpp:121-123), so not after this one.
    /// The NCBI exception of this error when it is not a database error of NCBI's: a search
    /// program or blastdbcmd writes `Error: <what()>` and ends with exit status 255.
    pub fn ncbi_exception(&self) -> Option<&crate::ncbi_exception::NcbiException> {
        match self {
            SeqDbError::Ncbi(exception) => Some(exception),
            SeqDbError::TaxdbThrows { exception, .. } => Some(exception.as_ref()),
            _ => None,
        }
    }
}

impl std::error::Error for SeqDbError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            SeqDbError::Open { source, .. }
            | SeqDbError::Io { source, .. }
            | SeqDbError::Map { source, .. }
            | SeqDbError::Read { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// The texts of `SeqDbError::Format` that are NCBI's own: the errors NCBI
/// raises for the same file with the same words. Every other format error is
/// LOSAT's (a damaged file that NCBI reads past the end of its mapping).
const NCBI_FORMAT_MESSAGES: [&str; 15] = [
    // seqdbfile.cpp:181-185
    "Error: Not a valid version 4 or 5 database.",
    // seqdbfile.cpp:222-226
    "Error: requested sequence type does not match DB.",
    // seqdbalias.cpp:810-814
    "Illegal configuration: DB alias files are mutually recursive.",
    // seqdbcol.cpp:220-238
    "Column file uses unknown format_version.",
    "Column file uses unknown data type.",
    "Column file uses unsupported offset size.",
    // seqdbcol.cpp:267-269, 288-290, 316-318
    "CSeqDBColumn: File format error.",
    // seqdbcol.cpp:300-302
    "CSeqDBColumn: Error; duplicate metadata key.",
    // seqdbblob.cpp:112-114, 193-195, 220-222
    "CBlastDbBlob::ReadVarInt: eof while reading integer.",
    "CBlastDbBlob::ReadString: Unterminated string.",
    "CBlastDbBlob::x_ReadRaw: hit end of data",
    // seqdbimpl.cpp:2375-2376
    "Error in stored mask algorithm description data.",
    // seqdb_lmdb.cpp:418-420, 678-680: an empty `.pos`/`.nos` or `.pot`/`.not`
    // (`tax_files::OidIndex::open`; NCBI BLAST+ 2.17.0 with -negative_seqidlist
    // and -negative_taxids: line WQA D18 review M4, probe 10 of round 2).
    "Failed to open oid-to-seqid lookup file",
    "Failed to open oid-to-taxids lookup file",
    // seqdbgimask.cpp:213-217
    "Gi-mask file uses unknown format_version.",
];

/// seqdblmdbset.cpp:224, 231: NCBI's text of a list of volumes of both versions
/// ([`SeqDbError::MixedVersions`]).
pub(crate) const MIXED_VERSIONS: &str = "DB list contains both Version 4 and Version 5 dbs";

/// seqdbvol.cpp:3365-3382 (the title or the counts follow). The duplicate of seqdbimpl.cpp:2456-2459,
/// whose text holds the volume's path, is a [`SeqDbError::NotFound`] of its bytes (`mask.rs`).
const NCBI_COLUMN_MESSAGE_PREFIXES: [&str; 2] = [
    "Error: duplicate column title [",
    "Error: column has wrong #oids [",
];

/// ncbifile.cpp:5893-5943: `CMemoryFile`'s text for a file it cannot open; the path, `':` and the
/// `strerror` text follow (`SeqDbError::memory_map_message`).
const MEMORY_FILE_PREFIX: &str = "CMemoryFile: Cannot memory map file '";

/// The text of `strerror` for EMFILE and ENFILE, which makes NCBI's exception of a file it cannot
/// map `eTooManyOpenFiles` (seqdbatlas.hpp:713-715).
const TOO_MANY_OPEN_FILES: &[u8] = b"Too many open files";

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:712-718
/// Whether NCBI's exception of a file it cannot map, whose message is `message`, is
/// `eTooManyOpenFiles`: whatever the message names, also a path (a missing file's message holds
/// it), NCBI tests the whole text.
fn too_many_open_files(message: &[u8]) -> bool {
    message
        .windows(TOO_MANY_OPEN_FILES.len())
        .any(|window| window == TOO_MANY_OPEN_FILES)
}

impl SeqDbError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:709-720
    // `BLAST_DATABASE_ERROR` is 2 (blast_app_util.hpp:153). A database file
    // that cannot be opened or mapped (`file::DbFile::open`) is the reader's
    // rethrow of `CMemoryFile`'s exception (`memory_file_exception`): its code
    // is `eTooManyOpenFiles` when that text holds "Too many open files" (the
    // `strerror` text of EMFILE and of ENFILE), the text then NCBI's own, and
    // `eMemoryMappingFailure` otherwise, `Database memory map file error - `
    // and the exception's text. NCBI maps a file whole where it opens it:
    // LOSAT reads the file's first byte there, and a first byte that the
    // system does not give has the text of that mapping (NCBI BLAST+ 2.17.0
    // with `mmap()` failing with EIO, and a 3 GiB header file:
    // `length=3221225472`; line WQA D18 probe 3, 2026-11-21). A read that
    // fails later ([`SeqDbError::Read`]: the system stopped giving the bytes,
    // the file became shorter) is where NCBI reads mapped bytes that are gone
    // and ends by SIGBUS: LOSAT's own error (review M1 of D18, 2026-11-21). On
    // Windows NCBI's text comes from another call, which LOSAT does not
    // reproduce: the error is LOSAT's there.
    /// The text NCBI BLAST+ gives for this error, when it is one of NCBI's
    /// own: a search program prints `BLAST Database error: <text>` and ends
    /// with exit status 2. `None`: the error is LOSAT's (a database that
    /// LOSAT refuses, a damaged file, a read that fails after the file's
    /// first byte: [`SeqDbError::is_unreadable`]), to be reported
    /// in LOSAT's words. The text is bytes ([`SeqDbError::NotFound`]).
    pub fn ncbi_message(&self) -> Option<Vec<u8>> {
        match self {
            SeqDbError::NotFound(message) => Some(message.clone()),
            SeqDbError::Format { message, .. } => (NCBI_FORMAT_MESSAGES
                .contains(&message.as_str())
                || NCBI_COLUMN_MESSAGE_PREFIXES
                    .iter()
                    .any(|prefix| message.starts_with(prefix)))
            .then(|| message.clone().into_bytes()),
            SeqDbError::MixedVersions { .. } => Some(MIXED_VERSIONS.as_bytes().to_vec()),
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-197
            SeqDbError::MemoryMap(_) | SeqDbError::Open { .. } | SeqDbError::Map { .. } => {
                let message = self.memory_map_message()?;
                if too_many_open_files(&message) {
                    Some(b"Too many open files, please raise the open file limit".to_vec())
                } else {
                    Some([b"Database memory map file error - ".as_slice(), &message].concat())
                }
            }
            SeqDbError::Io { .. }
            | SeqDbError::Read { .. }
            | SeqDbError::Integrity { .. }
            | SeqDbError::Unsupported { .. }
            | SeqDbError::NcbiCrashes { .. }
            | SeqDbError::NcbiThrows { .. }
            | SeqDbError::Range(_)
            | SeqDbError::Internal(_)
            | SeqDbError::Ncbi(_)
            | SeqDbError::TaxdbThrows { .. }
            | SeqDbError::Setting(_) => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5893-5947
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5600-5656
    /// The `CFileException` (`eMemoryMap`) of NCBI's `CMemoryFile` for this error, where it is
    /// one: a file that does not open ([`SeqDbError::Open`], `CMemoryFileMap::x_Open`) or whose
    /// first byte cannot be had where it was opened ([`SeqDbError::Map`], the mapping of the
    /// whole file, `CMemoryFileSegment`), with the system's text (`utils::os_error_text::ncbi_strerror`;
    /// `None` without one). The database reader's atlas rethrows it as its own error
    /// ([`SeqDbError::ncbi_message`]); the lookup files that NCBI maps without the atlas give it
    /// to the program as it is (`filter::lookup_file_error`).
    pub(crate) fn memory_file_exception(&self) -> Option<crate::ncbi_exception::NcbiException> {
        let (path, strerror) = match self {
            SeqDbError::Open { path, source } | SeqDbError::Map { path, source, .. } => {
                (path, crate::utils::os_error_text::ncbi_strerror(source)?)
            }
            _ => return None,
        };
        let path = crate::utils::os_bytes::path_bytes(path);
        let (line, function, message) = match self {
            SeqDbError::Map { len, .. } => (
                5650,
                "ncbi::CMemoryFileSegment::CMemoryFileSegment",
                [
                    b"Cannot map file '".as_slice(),
                    &path,
                    format!("' to memory (offset=0, length={len}): {strerror}").as_bytes(),
                ]
                .concat(),
            ),
            _ => (
                5945,
                "ncbi::CMemoryFileMap::x_Open",
                [
                    MEMORY_FILE_PREFIX.as_bytes(),
                    &path,
                    b"':",
                    strerror.as_bytes(),
                ]
                .concat(),
            ),
        };
        Some(crate::ncbi_exception::NcbiException::new(
            ThrowSite {
                file: "src/corelib/ncbifile.cpp",
                line,
                kind: "CFileException::eMemoryMap",
                function: Some(function),
            },
            message,
        ))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
    /// Whether NCBI's constructor fails here with a `CSeqDBException`, which NCBI's guess of a
    /// database's type catches to look for the name as a nucleotide database: an error in NCBI's
    /// words ([`Self::ncbi_message`]) or a failed integrity check ([`SeqDbError::Integrity`], as an
    /// OID list file of its header alone at the constructor's scan: probes RB_G3 and RB_G4 of
    /// review B of line RG D1b-1). LOSAT's refusals are not: NCBI may open that database.
    pub fn is_seqdb_exception(&self) -> bool {
        self.ncbi_message().is_some() || matches!(self, SeqDbError::Integrity { .. })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:229-252
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:227-236
    // A reader error inside NCBI's preliminary search loses its own text: the
    // search stage replaces it by one of two messages, and the program
    // reports a `CBlastException` (`BLAST engine error: <text>`, exit status
    // `BLAST_ENGINE_ERROR`, 3; blast_app_util.hpp:155, 216-219). A volume
    // whose sequence file is missing is reported this way, because NCBI opens
    // that file when the search reads its first sequence, and so is a failed
    // integrity check of a file the search reads (`SeqDbError::Integrity`: the
    // mask column's offsets and blobs, seqsrc_seqdb.cpp:318-321).
    /// The text NCBI BLAST+ gives when this error is raised while the search
    /// reads the sequences (or their masks) of the database: a search program
    /// prints `BLAST engine error: <text>` and ends with exit status 3. `None`:
    /// the error is LOSAT's ([`Self::ncbi_message`] is `None` too).
    pub fn ncbi_search_message(&self) -> Option<&'static str> {
        const TOO_MANY: &str = "Too many open files, please raise the open file limit";
        if let SeqDbError::Integrity { .. } = self {
            return Some("Database memory map file error");
        }
        match self.ncbi_message()? {
            message if message == TOO_MANY.as_bytes() => Some(TOO_MANY),
            _ => Some("Database memory map file error"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:241-256
    // `CLocalBlast::Run` catches `CSeqDBException` into `Database memory map file error` (or
    // `Too many open files`), and catches any other `std::exception` (such as `CMemoryFile`'s
    // `CFileException` on a lazy alias `SEQIDLIST`/`TAXIDLIST` or `CSeqIdException`) into a
    // `CBlastException` with code `eSystem` carrying its `what` text: engine error, exit status 3.
    /// The bytes after `BLAST engine error: ` when this error is raised during `CLocalBlast`'s
    /// preliminary stage: [`Self::ncbi_search_message`], or `what()` of a known NCBI exception
    /// ([`Self::ncbi_exception`]). `None`: the error is LOSAT's.
    pub fn ncbi_preliminary_message(&self) -> Option<std::borrow::Cow<'static, [u8]>> {
        if let Some(text) = self.ncbi_search_message() {
            return Some(std::borrow::Cow::Borrowed(text.as_bytes()));
        }
        self.ncbi_exception()
            .filter(|exception| exception.places_known())
            .map(|exception| std::borrow::Cow::Owned(exception.what()))
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:709-719
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:5893-5943
    // A file that NCBI cannot map is one exception of `CSeqDBFileMemMap::Init`, whose message is
    // the text of `CMemoryFile` (`eMemoryMappingFailure`, or `eTooManyOpenFiles` when that text
    // holds "Too many open files"): `To be memory mapped the file must exist: '<path>'` for a path
    // that is not a file (`SeqDbError::MemoryMap`), `CMemoryFile: Cannot memory map file
    // '<path>':<strerror>` for one that cannot be opened (`SeqDbError::Open`, cited at
    // `file::DbFile::open`; probe S2 of review A of line RG D1b-1: an OID list file of mode 000),
    // `Cannot map file '<path>' to memory (offset=0, length=<bytes>): <strerror>` for one whose
    // first byte cannot be read where it was opened (`SeqDbError::Map`; line WQA D18). Both of
    // these are the message of `memory_file_exception`, with the system's text of
    // `utils::os_error_text::ncbi_strerror` (Linux's on WASI too); on Windows NCBI's text comes
    // from another call (`WIN_LAST_ERROR_STR`), which LOSAT does not reproduce: the error is then
    // LOSAT's own.
    /// The message of NCBI's exception for a file it cannot map, without the text that the
    /// programs' `catch` puts before it; `None` for another error.
    fn memory_map_message(&self) -> Option<Vec<u8>> {
        match self {
            SeqDbError::MemoryMap(message) => Some(message.clone()),
            SeqDbError::Open { .. } | SeqDbError::Map { .. } => self
                .memory_file_exception()
                .map(|exception| exception.get_msg().to_vec()),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:888-895
    /// The text blastdbcmd gives for this error where a request meets it (its `catch` posts
    /// `Error: [blastdbcmd] <GetMsg()>`, exit status 1), which is also the message of NCBI's
    /// exception that dustmasker and segmasker report (`throw_site`): the message of a file NCBI
    /// cannot map without the search programs' text before it (`memory_map_message`), the message
    /// of an exception of NCBI's that is not a `CSeqDBException` (`SeqDbError::Ncbi`: a lookup
    /// file's `CFileException`, `filter::lookup_file_error`; NCBI BLAST+ 2.17.0's blastdbcmd with
    /// a `.pto` of mode 000: line WQA D18 probe 8, 2026-11-21) or of a taxonomy file
    /// (`SeqDbError::TaxdbThrows`, made by `memory_file_exception`: blastdbcmd's FASTA of a
    /// sequence with a taxid, line WQA D18b probe 12), else [`Self::ncbi_message`].
    pub fn ncbi_request_message(&self) -> Option<Vec<u8>> {
        match self {
            SeqDbError::Ncbi(exception) => Some(exception.get_msg().to_vec()),
            SeqDbError::TaxdbThrows { exception, .. } => Some(exception.get_msg().to_vec()),
            _ => self.memory_map_message().or_else(|| self.ncbi_message()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:444-450
// (`CSeqDBException`, seqdbcommon.hpp:73-100, does not name its codes, so every one is reported
// as `eInvalid`: `(CSeqDBException::eInvalid) No alias or index file found ...` for the
// `eFileErr` of seqdbalias.cpp:311-321, fixture db_missing.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:311-322
// (`NCBI_THROW` is on line 320; its function is `BLASTDB::ncbi::CSeqDBAliasNode::x_ResolveNames`,
// the module `BLASTDB` before it: segmasker 2.17.0's `e.what()`, fixture db_missing.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:740-752, 921-925
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:155-164 (the
// constructor of `CSeqDB` from `dbname` and `seqtype`)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:123-128
// (Each of these texts has this one place in NCBI's database reader; segmasker 2.17.0 names the
// lines 750, 925, 162 and 128 and these functions, fixtures db_alias_*, db_empty_name,
// db_space_name; review of line RG item K2, L2.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:124-142
// (The three texts of a version 5 database without the table, `seqdb::tax_files::NO_VOLUME_INFO`
// and `NO_ACCESSION_INFO`, come from the one throw after the `switch`; NCBI's official build
// names its line 142, not the 141 of the source of 598d8ae6, and the module and class before the
// function: blastp 2.17.0's `Error: [blastp] DB contains no accession info.NCBI C++ Exception:`
// of outfmt 15 at a lookup of a list of version 5 databases with and without the accession table,
// `blastdb-ids-20261021/R/mixdb2/review1/out/ncbi/bp.lim1.U2.15.err`, DECISIONS 62 of line DBID.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:56-66
// (A reader opens each LMDB file read-only through `InitDbi`, which leaves only the accession
// table unset: a file without `volname` or `volinfo` throws lmdb's `mdb_dbi_open: MDB_NOTFOUND:
// No matching key/data pair found` there, before `GetDbi` can throw `DB contains no vol info.`
// (segmasker and dustmasker 2.17.0: `Error: (106.7) Application's execution failed: mdb_dbi_open:
// ...`, exit 3; blastp and blastn: `Error: mdb_dbi_open: ...`, exit 255; review 1 of line DBID
// pfoff-1, L4). So no text of the reader's table below is GetDbi's: its one throw that NCBI's
// programs reach is the lookup's (`lookup_throw_site`).)
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:699-722
// (The OID list file of an alias file's `OIDLIST` that is not a file, `SeqDbError::MemoryMap`:
// segmasker 2.17.0 names the line 718 of the official build and the function without a module,
// probe mk of line RG D1b.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:178-182
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:218-226
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:861-864
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:78-96
// (An index file of another format version, an index of the other molecule type, an empty binary
// list file, and the LMDB file of a version 5 volume that is missing or not an LMDB file:
// blastdbcheck 2.17.0 names these lines and functions, rows pp_p4_ver7, d_p4_pin_nucl_prot,
// a_gizero, d_n5_no_ndb and pp_p5_pdb_cut of row RG-D2b. The missing
// LMDB file's text starts with its path: `LMDB_MISSING_TAIL`.)
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:879-885
// (`SetMemBit` throws in the inline header without a function name in NCBI's report; blastdbcheck
// 2.17.0 at verbosity 3/4 reports `include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp`, line
// 881, with no function prefix.)
/// The places of NCBI's throws of a database error of the reader, by the start of its text.
const DATABASE_ERROR_SITES: [(&str, ThrowSite); 12] = [
    (
        "MEMB_BIT error: conflicting bit found.",
        ThrowSite {
            file: "include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp",
            line: 881,
            kind: "CSeqDBException::eInvalid",
            function: None,
        },
    ),
    (
        "Error: requested sequence type does not match DB.",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdbfile.cpp",
            223,
            "BLASTDB::ncbi::CSeqDBIdxFile::CSeqDBIdxFile",
        ),
    ),
    (
        "Error: Not a valid version 4 or 5 database.",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdbfile.cpp",
            182,
            "BLASTDB::ncbi::CSeqDBIdxFile::CSeqDBIdxFile",
        ),
    ),
    (
        "Specified file is empty.",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdbcommon.cpp",
            864,
            "BLASTDB::ncbi::s_SeqDB_IsBinaryNumericList",
        ),
    ),
    (
        "LMDB runtime error: ",
        seqdb_site(SEQDB_LMDB, 95, LMDB_ENV_OPEN),
    ),
    (
        "No alias or index file found for ",
        seqdb_site(
            SEQDBALIAS,
            320,
            "BLASTDB::ncbi::CSeqDBAliasNode::x_ResolveNames",
        ),
    ),
    (
        "Could not find volume or alias file (",
        seqdb_site(
            SEQDBALIAS,
            925,
            "BLASTDB::ncbi::CSeqDBAliasNode::x_ExpandAliases",
        ),
    ),
    (
        "No database names were ",
        seqdb_site(
            SEQDBALIAS,
            750,
            "BLASTDB::ncbi::CSeqDBAliasNode::x_ExpandAliases",
        ),
    ),
    (
        "Database name is required.",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdb.cpp",
            162,
            "BLASTDB::ncbi::CSeqDB::CSeqDB",
        ),
    ),
    (
        "Error: File (",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdbfile.cpp",
            128,
            "BLASTDB::ncbi::CSeqDBExtFile::CSeqDBExtFile",
        ),
    ),
    (
        "To be memory mapped the file must exist: '",
        seqdb_site(SEQDBATLAS, 718, "ncbi::CSeqDBFileMemMap::Init"),
    ),
    (
        MEMORY_FILE_PREFIX,
        seqdb_site(SEQDBATLAS, 718, "ncbi::CSeqDBFileMemMap::Init"),
    ),
];

/// NCBI's reader of the LMDB files of version 5 volumes.
const SEQDB_LMDB: &str = "src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp";

/// The function that opens a version 5 volume's LMDB file for reading.
const LMDB_ENV_OPEN: &str = "BLASTDB::ncbi::CBlastLMDBManager::CBlastEnv::CBlastEnv";

/// The end of the text of a version 5 volume's missing LMDB file, after `File <path>`.
const LMDB_MISSING_TAIL: &str = " not found. If you renamed any BLAST database files, please use \
original file names, and makeblastdb to rename the database. If you deleted any BLAST database \
files, you need to recreate the database.";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:84-89
/// NCBI's error of a version 5 volume's LMDB file `path` that is not there (`CBlastEnv`): the one
/// producer of its text, which `throw_site` places.
pub(crate) fn missing_lmdb_file(path: &Path) -> SeqDbError {
    SeqDbError::NotFound(
        [
            b"File ".as_slice(),
            &crate::utils::os_bytes::path_bytes(path),
            LMDB_MISSING_TAIL.as_bytes(),
        ]
        .concat(),
    )
}

/// NCBI's header of the memory-mapped files of the reader.
const SEQDBATLAS: &str = "include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp";

/// The place of `CBlastLMDBManager::CBlastEnv::GetDbi`'s throw of a version 5 database's missing
/// table (`DB contains no ...`).
const GET_DBI_SITE: ThrowSite = seqdb_site(
    "src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp",
    142,
    "BLASTDB::ncbi::CBlastLMDBManager::CBlastEnv::GetDbi",
);

/// NCBI's alias reader.
const SEQDBALIAS: &str = "src/objtools/blast/seqdb_reader/seqdbalias.cpp";

/// The place of a `CSeqDBException` (its codes are all reported as `eInvalid`).
const fn seqdb_site(file: &'static str, line: u32, function: &'static str) -> ThrowSite {
    ThrowSite {
        file,
        line,
        kind: "CSeqDBException::eInvalid",
        function: Some(function),
    }
}

/// The place in NCBI's database reader that throws the `CSeqDBException` whose message is
/// `message` (a text of [`SeqDbError::ncbi_message`]), where LOSAT knows it: the one owner of
/// those places (dustmasker's and segmasker's reports of the reader's errors,
/// `seqmasks_io::bdb_reader`; a lookup's own throw is [`lookup_throw_site`]'s). Each place
/// reports the code `eInvalid` (`CSeqDBException` names none). `None` for every other text, of
/// which LOSAT does not know the place: also the texts that `ncbi_message` gives of the codes
/// `eTooManyOpenFiles` and `eMemoryMappingFailure` (but the exception's own message of a file
/// NCBI cannot map, `ncbi_request_message`), which are the texts of the programs' `CATCH_ALL`
/// and not the exception's own message, the message of `eTooManyOpenFiles` (another line of
/// `CSeqDBFileMemMap::Init`), and `DB contains no ...`, which NCBI's reader never throws (it
/// throws lmdb's error first, `InitDbi` cited above).
pub(crate) fn throw_site(message: &[u8]) -> Option<ThrowSite> {
    if too_many_open_files(message) {
        return None;
    }
    if message.starts_with(b"File ") && message.ends_with(LMDB_MISSING_TAIL.as_bytes()) {
        return Some(seqdb_site(SEQDB_LMDB, 88, LMDB_ENV_OPEN));
    }
    DATABASE_ERROR_SITES
        .iter()
        .find(|(start, _)| message.starts_with(start.as_bytes()))
        .map(|&(_, site)| site)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:213-233
// (A version 4 volume after version 5 ones throws on line 225, a version 5 volume after version 4
// ones on line 232: blastdbcheck 2.17.0 names these lines, rows mix45_v4 and mix45_v5 of row
// RG-D2b.)
impl SeqDbError {
    /// The place in NCBI's database reader that throws this error's `CSeqDBException`, where
    /// LOSAT knows it: the throw that [`SeqDbError::MixedVersions`] names, else the place of its
    /// text ([`throw_site`]).
    pub(crate) fn throw_site(&self) -> Option<ThrowSite> {
        if let SeqDbError::MixedVersions {
            version4_after_version5,
            ..
        } = self
        {
            let line = if *version4_after_version5 { 225 } else { 232 };
            return Some(seqdb_site(
                "src/objtools/blast/seqdb_reader/seqdblmdbset.cpp",
                line,
                "BLASTDB::ncbi::CSeqDBLMDBSet::CSeqDBLMDBSet",
            ));
        }
        throw_site(&self.ncbi_request_message()?)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:167-172
/// The place of NCBI's throw where its lookup of a Seq-align id fails with the message `message`
/// (a text of `report::db_subject::LookupFailure::ncbi_message`), for the texts that NCBI's
/// lookup was seen to throw: the accession table of a version 5 database (`GetReadEnvAcc`'s
/// `GetDbi`, "DB contains no accession info.": blastp 2.17.0's outfmt 15 at a lookup of a list of
/// version 5 databases with and without the table, DECISIONS 62 (1) and (2) of line DBID). `None`
/// for every other text, also those of [`throw_site`], which no lookup was seen to throw.
pub(crate) fn lookup_throw_site(message: &[u8]) -> Option<ThrowSite> {
    (message == tax_files::NO_ACCESSION_INFO.as_bytes()).then_some(GET_DBI_SITE)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
// `NStr::TruncateSpaces` (`eTrunc_Both`) is `ncbi_str::truncate_spaces_str`: the `isspace` of
// the "C" locale (space, tab, line feed, vertical tab, form feed and carriage return). The
// reader itself (`CSeqDB`, `SeqDb::open_in`) takes the name as it is given: a name of spaces
// alone is "No database names were passed to CSeqDB::CSeqDB()." there, while the programs,
// which trim first, report "Database name is required.".
/// The database name that a BLAST program hands to the reader for the value
/// of its `-db` argument: the value without the white space at its ends.
/// A search program passes the result to [`SeqDb::open_in`].
pub fn db_argument_name(value: &str) -> &str {
    crate::ncbi_str::truncate_spaces_str(value)
}

/// A contiguous range of sequences inside one volume, as planned by
/// [`SeqDb::next_block`] and read by [`SeqDb::read_block`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockRange {
    /// Index of the volume in the database.
    pub volume: usize,
    /// First sequence of the range, relative to the volume.
    pub first: u32,
    /// One past the last sequence of the range, relative to the volume.
    pub end: u32,
    /// Database ordinal id of the first sequence.
    pub first_ordinal: u32,
}

impl BlockRange {
    /// Number of sequences in the range.
    pub fn len(&self) -> u32 {
        self.end - self.first
    }

    /// Whether the range holds no sequence (never true for a planned range).
    pub fn is_empty(&self) -> bool {
        self.end == self.first
    }

    /// Database ordinal id one past the last sequence of the range.
    pub fn end_ordinal(&self) -> u32 {
        self.first_ordinal + self.len()
    }
}

/// Where the reads of one sequence find it (NCBI's `CSeqDBVolSet::FindVol`): its volume, its
/// OID in the volume, and the id that the readers' messages and blocks name (the ordinal id or
/// the OID that the caller asked for). [`SeqDb::locate_oid`] finds it by NCBI's OID, as every
/// read of NCBI's `CSeqDB` does; `SeqDb::find_volume` by ordinal id.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeqLocation {
    /// Index of the volume in the database.
    volume: usize,
    /// The sequence's OID in the volume.
    vol_oid: u32,
    /// The id the caller gave (an ordinal id or an OID): the first ordinal id of the block that
    /// `read_sequence_at` reads. The messages name the OID (`SeqDb::location_oid`).
    id: u32,
}

/// A database resolved with only the headers of its volumes read
/// ([`SeqDb::locate_one_name_metadata`]).
#[derive(Debug)]
pub(crate) struct LocatedDatabase {
    /// `GetTitle()`.
    pub title: Vec<u8>,
    /// `GetDate()`, or the exception of the `CTime` that NCBI's `GetDBMetaData` builds from it.
    pub date: Result<String, date::DateError>,
    /// Whether every volume names an LMDB file: a version 5 database.
    pub version5: bool,
    /// The volume paths without extension, in NCBI's order.
    pub volumes: Vec<Vec<u8>>,
    /// The number of sequences and the total length that the header of each volume holds,
    /// in the order of `volumes`.
    pub volume_totals: Vec<(u32, u64)>,
    /// The alias files (paths with extension), sorted.
    pub alias_files: Vec<Vec<u8>>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:375-391
/// What NCBI's reports read of one name of a `-db` list (`CAlignFormatUtil::SDbInfo`): the
/// database of that name opened alone and without a list. The `Database:` lines of the
/// pairwise report, the database counts of XML and JSON and the title of the Seq-annot
/// formats are made of these, one per name.
#[derive(Debug)]
pub struct NameRecord {
    /// The name as the reports take it (`alias::report_names`): `GetDBNameList`.
    pub name: String,
    /// The title of the database (`GetTitle`), empty when it has none.
    pub title: Vec<u8>,
    /// `GetDate`.
    pub date: String,
    /// `GetNumSeqs`: NCBI's `int` (the walk of `NSEQ` converted, so negative for a negative
    /// `NSEQ` or one of 2^31 or more).
    pub num_seqs: i32,
    /// `GetTotalLength`: NCBI's `Uint8`.
    pub total_length: u64,
    /// The filtering algorithms of the name's database, which the formatter of a search with a
    /// database mask reads (`GetMaskAlgorithmDetails`, `SeqDb::name_mask_details`), or the error
    /// of reading them, met only where the formatter asks for them; read with the record when
    /// the search has a database mask (`read_name_records`), `None` otherwise.
    pub mask_algorithms: Option<Result<MaskAlgorithms, SeqDbError>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:163-209
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:905-938 (cited at
// `SeqDb::statistics_total_length`) and c++/src/algo/blast/core/blast_engine.c:1407, 1434
/// The numbers of a database that NCBI's statistics read, with NCBI's types (line RG D1a-3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DbStatistics {
    /// The setup's total length (`BLAST_GapAlignSetUp`, `Int8`): `STATS_TOTLEN`, else
    /// `GetTotalLength`, when positive; else the length of OID 0 (0 without an OID:
    /// `SeqDb::first_oid_read_fails`).
    pub setup_total: i64,
    /// The setup's number of sequences (`Int4`): `STATS_NSEQ`, else `GetNumSeqs`, with a
    /// positive total; else 1 (0 without an OID). Negative for a negative `NSEQ` or one of 2^31
    /// or more (the `int` of NCBI's walk).
    pub setup_count: i32,
    /// The Gumbel block's database length before `-dbsize` (`sbp->gbp->db_length`, before the
    /// division of a translated subject): the setup's total before its fallback to OID 0, so 0
    /// for a `LENGTH` of 0 and negative for one of 2^63 or more.
    pub gumbel_db_length: i64,
    /// `BlastSeqSrcGetTotLen` (`Int8` of `GetTotalLength`), which the preliminary search and the
    /// traceback read: 0 computes the statistics of each subject on its own, and only a
    /// positive one sorts the results by e-value at the end of the traceback.
    pub engine_db_length: i64,
    /// The number of sequences of the requests' effective search space (`STATS_NSEQ`, else
    /// `GetNumSeqs`, without the setup's fallback to 1; split_query_aux_priv.cpp:164-173, quoted
    /// at `OptionsHandle::set_database_filter`).
    pub request_count: i32,
    /// `GetMinLength` (`SeqDb::min_length`).
    pub min_length: i32,
    /// `BlastSeqSrcGetAvgSeqLen` (`SeqDb::average_length`).
    pub average: i32,
}

impl DbStatistics {
    /// Whether NCBI's preliminary search computes the statistics of each subject on its own
    /// (`db_length == 0`, blast_engine.c:1434).
    pub fn per_subject(&self) -> bool {
        self.engine_db_length == 0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1772
    /// Whether NCBI's traceback sorts each query's hit list by e-value at its end: not for a
    /// total length of 0 or of 2^63 or more, where the hit lists keep the traceback's order.
    pub fn sorts_results(&self) -> bool {
        self.engine_db_length > 0
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
/// What NCBI's `CSeqDB` answers for the shortest length when the alias walk gives none
/// (`SeqDb::min_length`).
pub(crate) const BLAST_SEQSRC_MINLENGTH: i32 = 10;

/// An alias tree after NCBI's `ComputeMasks` and the paths of its volumes in the order of
/// `SeqDb::volumes` (`SeqDb::listed_tree`).
#[derive(Debug)]
struct ListedTree {
    root: alias::AliasNode,
    names: Vec<Vec<u8>>,
}

/// An open database: one volume, or the volumes of an alias file.
#[derive(Debug)]
pub struct SeqDb {
    seqtype: SeqType,
    volumes: Vec<Volume>,
    /// Ordinal id of the first sequence of each volume, then the total.
    volume_starts: Vec<u32>,
    title: Vec<u8>,
    date: String,
    /// The number of included OIDs: ordinal ids run from 0 to it (the OIDs of the alias files'
    /// ranges and of a list, `filter`). Iterations and schedules count these; NCBI's counts are
    /// `counts`.
    num_seqs: u32,
    /// What NCBI's `CSeqDB` of the search answers (`GetNumSeqs`, `GetTotalLength`,
    /// `GetMaxLength`, `GetMinLength`): the walks of the alias tree without a list
    /// (`alias::walk_totals`; the longest length is the volumes'), the scan of the included OIDs
    /// with one (`scan_totals`).
    counts: Counts,
    /// `GetNumSeqsStats` and `GetTotalLengthStats`: the walks of `STATS_NSEQ` and
    /// `STATS_TOTLEN`, 0 without the keys; a list does not change them.
    stats_num_seqs: i32,
    stats_total_length: u64,
    /// `GetNumSeqs` and `GetTotalLength` of the database opened again without a list, which the
    /// reports print (`printed_num_seqs`, `printed_total_length`): the walks.
    printed_num_seqs: i32,
    printed_total_length: u64,
    /// The sum of the volumes' lengths (`CSeqDBVolSet::GetVolumeSetLength`).
    volume_set_length: u64,
    /// `SeqDb::included_length` of a `filter`, made at its first use (none before). Only the
    /// schedules of the parallel searches read it (`SeqDb::set_included_length`).
    #[cfg(any(test, feature = "parallel"))]
    included_length: OnceLock<u64>,
    /// The length of OID 0 of the volumes, which the statistics read when their total length
    /// is not positive (`statistics_numbers`); `None` when it is.
    first_oid_length: Option<i64>,
    /// The walks of the alias tree, with the exception of each one whose value does not
    /// convert (`check_alias_values`).
    alias_totals: alias::AliasTotals,
    /// Whether NCBI's OID bit set is the special all-set one after the alias files' OID ranges
    /// (`alias::AliasOids`): a positive list left with no id then keeps every OID
    /// (`SeqDb::restrict_seqids`).
    all_bits_special_set: bool,
    /// NCBI's `m_OidListSetup`: whether its OID list was built, which a database opened
    /// without a list does at its first use (`SeqDb::check_or_find_oid`).
    oid_list_set_up: AtomicBool,
    /// Whether a list was applied (`SeqDb::restrict`): NCBI's `CSeqDB` then scans the totals.
    list_applied: bool,
    /// Whether the alias files make NCBI's `CSeqDB` scan the totals without a list (an
    /// `OIDLIST` not covered by `NSEQ` and `LENGTH`, `alias::needs_totals_scan`): the counts and
    /// the printed counts are then the scan of the included OIDs.
    alias_scanned: bool,
    /// The iteration range that NCBI's first build of the alias files' OID list sets
    /// (`alias::top_iteration_range`, every OID without it), inside which a list's scan starts
    /// (`SeqDb::restrict`).
    alias_range: (u32, u32),
    /// The failure of an OID list file of the alias files, which every build of NCBI's OID list
    /// meets (`SeqDb::build_oid_list`; a search's reads, `SeqDb::fail_search_reads_with_oid_list`).
    oid_list_error: Option<oid_list::OidListError>,
    /// Whether NCBI's constructor builds the list to scan the totals
    /// (`alias::needs_totals_scan`) and so fails with `oid_list_error`: after the warnings and
    /// errors of a list (`SeqDb::restrict`), else before the walks of the alias files' counts
    /// (`SeqDb::check_alias_values`).
    oid_list_fails_scan: bool,
    /// The alias tree, kept where it names id lists that NCBI's OID list translates
    /// (`alias::reads_id_lists`): a positive list translates them again (`SeqDb::restrict`).
    listed_tree: Option<Box<ListedTree>>,
    /// The failure of the alias files' OID list without a list, kept where a positive list then
    /// translated the alias files' id lists without failing and NCBI's constructor scans the
    /// totals without a list (`oid_list_fails_scan`): the reports' database, opened again without
    /// the list, meets it (`SeqDb::read_name_records`).
    unlisted_oid_list_error: Option<oid_list::OidListError>,
    /// NCBI's OID list (the alias files' ranges, a list) and its iteration range, when there
    /// is one: ordinal id `i` is the `i`-th OID of the list inside the range.
    filter: Option<Box<filter::OidFilter>>,
    /// The search path the database was found on (`open_in`), for files NCBI looks for
    /// beside it (`taxdb.bti`).
    search_path: Option<SearchPath>,
    /// `-show_gis` of the search that reads the headers for its reports.
    reports_show_gis: bool,
    /// Whether the database was opened with a positive list (`SeqDb::restrict`).
    positive_list: bool,
    /// `SeqDb::list_taxids`.
    list_taxids: Vec<i32>,
    /// `SeqDb::list_gis`.
    list_gis: Vec<i64>,
    /// `SeqDb::name_records`.
    name_records: Vec<NameRecord>,
    /// The LMDB files of a version 5 database for `SeqDb::seqid_lookup`, opened at its first use.
    lookup_files: OnceLock<lookup::LookupFiles>,
    /// The ISAMs of the volumes for `SeqDb::seqid_lookup`, each opened at its first use.
    lookup_isams: OnceLock<lookup::VolumeIsams>,
    /// The filtering algorithms of the mask columns, built on first use (`mask.rs`).
    mask_algorithms: OnceLock<MaskAlgorithms>,
    /// Whether the run posted NCBI's warning of a UTF8String read for a VisibleString
    /// (`SeqDb::take_utf8_post`).
    utf8_posted: AtomicBool,
    /// Whether a reader met a reported sequence whose Seq-align id the lookup gives first to
    /// another sequence with such a header (`SeqDb::note_scope_utf8`).
    scope_utf8: AtomicBool,
    /// Whether the model of the run's loads has reached the batch whose own reads post that
    /// warning (`SeqDb::end_scope_utf8_loads`).
    scope_utf8_ended: AtomicBool,
    /// Whether a load of that model posts the warning (`SeqDb::claim_utf8_post`).
    utf8_claimed: AtomicBool,
    /// Whether the loads modelled now are not NCBI's (`SeqDb::without_scope_utf8_loads`).
    scope_utf8_held: AtomicBool,
}

/// NCBI's counts of a database (`CSeqDBImpl::m_NumSeqs`, `m_TotalLength`, `m_MaxLength`,
/// `m_MinLength`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Counts {
    num_seqs: i32,
    total_length: u64,
    max_length: u32,
    min_length: i32,
}

/// The walks of the alias files' counts that NCBI's `CSeqDB` constructor makes without a list,
/// in its order (seqdbimpl.cpp:117-139, cited at [`SeqDb::open_every_oid_in`]).
const UNLISTED_WALK_KEYS: &[alias::WalkKey] = &[
    alias::WalkKey::Nseq,
    alias::WalkKey::Length,
    alias::WalkKey::MinSeqLength,
    alias::WalkKey::StatsNseq,
    alias::WalkKey::StatsTotlen,
];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1093-1102
fn fix_string(s: &[u8]) -> &[u8] {
    match s.iter().position(|&byte| byte == 0) {
        Some(at) => &s[..at],
        None => s,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:217-233
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:85-89
// NCBI opens the LMDB file of a version 5 database when it opens the
// database (seqdblmdbset.cpp:235-262, 35-74) and compares the file's
// volume list with the volumes: `lmdb_entries`.
/// The checks of `CSeqDBLMDBSet`'s constructor on the volumes' `bases` (paths without
/// extension) and the names of the LMDB files their indexes hold: a version 4 volume among
/// version 5 ones (or the other way round), and an LMDB file that is not there. A volume is
/// version 5 when its index names an LMDB file.
fn check_lmdb_names(volumes: &[(PathBuf, Option<&[u8]>)]) -> Result<(), SeqDbError> {
    let names_lmdb = |name: &Option<&[u8]>| name.is_some_and(|name| !name.is_empty());
    if let Some((_, first)) = volumes.first() {
        let version5 = names_lmdb(first);
        if let Some((other, _)) = volumes
            .iter()
            .find(|(_, name)| names_lmdb(name) != version5)
        {
            return Err(SeqDbError::MixedVersions {
                path: other.clone(),
                version4_after_version5: version5,
            });
        }
    }
    for (base, name) in volumes {
        if !names_lmdb(name) {
            continue;
        }
        // The name's bytes, as the volume's header holds them (probe R8 of line LIB-B: a
        // name that is not UTF-8).
        let lmdb_name = name.unwrap_or_default();
        let lmdb_path = base.with_file_name(crate::utils::os_bytes::path_from_bytes(lmdb_name));
        if !std::fs::metadata(&lmdb_path).is_ok_and(|meta| meta.is_file()) {
            return Err(missing_lmdb_file(&lmdb_path));
        }
    }
    Ok(())
}

/// Resolved volume bases and file paths that reading a BLAST database can open
/// ([`SeqDb::read_paths`]).
#[derive(Debug, Default)]
pub(crate) struct DbReadPaths {
    /// Resolved volume paths without extension, in NCBI's `FindVolumePaths` order.
    pub(crate) volumes: Vec<Vec<u8>>,
    /// All files that reading the database can open: alias files, alias filter list files
    /// (`OIDLIST`, `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`), `MASKLIST` GI-mask files,
    /// volume component files (`VOLUME_EXTENSIONS`, `LMDB_EXTENSIONS`), and for protein databases,
    /// RPS files (`rps::RPS_EXTENSIONS`).
    pub(crate) files: Vec<Vec<u8>>,
}

fn push_mol_extensions(
    out: &mut std::collections::BTreeSet<Vec<u8>>,
    base: &[u8],
    mol: &[u8],
    extensions: &[&str],
) {
    for ext in extensions {
        out.insert([base, b".", mol, ext.as_bytes()].concat());
    }
}

fn strip_ascii_suffix<'a>(bytes: &'a [u8], suffix: &[u8]) -> Option<&'a [u8]> {
    let len = bytes.len().checked_sub(suffix.len())?;
    let (head, tail) = bytes.split_at(len);
    tail.eq_ignore_ascii_case(suffix).then_some(head)
}

impl SeqDb {
    /// Opens the database `name`: a path without extension, absolute or
    /// relative to the current directory. `name.{p|n}al` is read when it
    /// exists, else the single volume `name.{p|n}in`. This is
    /// [`SeqDb::open_in`] with a search path that holds the current directory
    /// only: `BLASTDB` and the registry are not consulted.
    #[cfg(test)]
    pub fn open(name: &str, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let cwd = std::env::current_dir().unwrap_or_default();
        Self::open_in(name.as_bytes(), seqtype, &SearchPath::new(&cwd, None, None))
    }

    /// Opens the database that a search program's `-db NAME` names, as NCBI
    /// BLAST+ finds it: `name` is looked for in the directories of `path`
    /// (`alias::resolve_in`). The name is a list in NCBI's syntax (spaces
    /// separate names, double quotes keep one name together): the databases of
    /// several names are searched as one. A name carries no file extension.
    ///
    /// A database that stores parsed Seq-ids opens; `has_parsed_seqids` tells
    /// it from one that holds only the ordinal ids.
    ///
    /// Column files (the mask data of `makeblastdb -mask_data`, user columns) are not read
    /// when the database is opened: [`SeqDb::mask_algorithms`] and [`SeqDb::mask_ranges_at`]
    /// (`losat blastdbcmd`) open them on first use, as NCBI does, and LOSAT's searches never
    /// ask for them (`Volume::open_as`). The name is the bytes of NCBI's `std::string`
    /// (`utils::os_bytes`), UTF-8 or not.
    pub fn open_in(name: &[u8], seqtype: SeqType, path: &SearchPath) -> Result<Self, SeqDbError> {
        let db = Self::open_search_in(name, seqtype, path)?;
        db.check_alias_values()?;
        Ok(db)
    }

    /// [`SeqDb::open_in`] without the errors of the alias files' counts, which depend on a
    /// list that the caller applies next ([`SeqDb::restrict`]); the caller then calls
    /// [`SeqDb::check_alias_values`] before it reads the database. The errors of the OID
    /// ranges (`FIRST_OID`, `LAST_OID`) come here, as NCBI's `ComputeMasks` runs first; those of
    /// the alias files' OID lists where NCBI's constructor scans the totals come after the list
    /// (`restrict`) or in `check_alias_values`.
    pub fn open_search_in(
        name: &[u8],
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<Self, SeqDbError> {
        let mut tree = alias::resolve_in(name, seqtype, path)?;
        let gi_masks = tree.gi_masks.take();
        let mut db = Self::open_tree(tree, name, seqtype)?;
        db.search_path = Some(path.clone());
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:81-91
        // (`m_UseGiMask` from the alias file's mask list): a database whose alias file names GI
        // masks has them as its filtering algorithms, and its mask columns are not read. The mask
        // files are looked for along the search path; nothing is read here.
        if let Some(names) = gi_masks {
            let masks = gi_mask::GiMasks::new(names, path.clone(), seqtype, db.volumes.len());
            db.mask_algorithms = OnceLock::from(MaskAlgorithms::gi(std::sync::Arc::new(masks)));
        }
        Ok(db)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:150-153
    // In NCBI's formatter (`CBlastFormat`), made by a search program after it opened the
    // database with its lists and after `Query is Empty!`; `blastdbcmd` makes none.
    /// Reads what the reports show of each name of the `-db` value `name` (`NameRecord`,
    /// `name_records`): each name of the list is opened again, alone, along the search path
    /// the database was opened with; a database of one name is this one, before any list. With
    /// `masks` (a search with `-db_soft_mask` or `-db_hard_mask`) each record also keeps its
    /// database's filtering algorithms (`NameRecord::mask_algorithms`): the one opening of each
    /// name of the formatter (`GetBlastDbInfo`).
    pub fn read_name_records(&mut self, name: &str, masks: bool) -> Result<(), SeqDbError> {
        let names = alias::report_names(name);
        self.name_records = match names.as_slice() {
            [one] if one == name => {
                // The reports' database is this one opened again without the list: after a
                // list, its walks of the counts meet their errors here (`check_alias_values`),
                // unless its alias files make it scan the totals too, where the build of its OID
                // list without the list fails first.
                if self.list_applied {
                    if let Some(failure) = &self.unlisted_oid_list_error {
                        return Err(failure.to_error());
                    }
                    if !self.alias_scanned {
                        use alias::WalkKey;
                        self.walk_errors(&[WalkKey::Nseq, WalkKey::Length, WalkKey::MinSeqLength])?;
                    }
                }
                vec![self.name_record(one, masks)]
            }
            _ => {
                let path = match &self.search_path {
                    Some(path) => path.clone(),
                    None => SearchPath::current_dir(),
                };
                let mut records = Vec::with_capacity(names.len());
                for one in &names {
                    let mut alone = Self::open_in(one.as_bytes(), self.seqtype, &path)?;
                    alone.read_name_records(one, masks)?;
                    records.append(&mut alone.name_records);
                }
                records
            }
        };
        Ok(())
    }

    /// What the reports read of this database under the name `name`, before any list; with
    /// `masks`, its filtering algorithms too ([`SeqDb::mask_algorithms`], a copy of the table).
    fn name_record(&self, name: &str, masks: bool) -> NameRecord {
        NameRecord {
            name: name.to_string(),
            title: self.title.clone(),
            date: self.date.clone(),
            num_seqs: self.printed_num_seqs(),
            total_length: self.printed_total_length(),
            mask_algorithms: masks.then(|| self.mask_algorithms().cloned()),
        }
    }

    /// One record per name of the `-db` list in the order of the list, a name given twice
    /// twice (`NameRecord`; empty before `read_name_records`).
    pub fn name_records(&self) -> &[NameRecord] {
        &self.name_records
    }

    /// Whether `taxdb.bti` is found on the search path the database was opened with
    /// (`SeqDB_ResolveDbPath` of `"taxdb.bti"`, data4xml2format.cpp:x_InitCommon): NCBI's XML2 and
    /// JSON formats then write the scientific name of a taxid.
    pub fn taxdb_found(&self) -> bool {
        self.search_path
            .as_ref()
            .is_some_and(|path| path.resolve_file(b"taxdb.bti").is_some())
    }

    /// Whether the reports of the search show the headers it reads with `-show_gis`.
    pub fn reports_show_gis(&self) -> bool {
        self.reports_show_gis
    }

    /// Sets `reports_show_gis`.
    pub fn set_reports_show_gis(&mut self, show_gis: bool) {
        self.reports_show_gis = show_gis;
    }

    /// Whether the database was opened with a positive list (`CSeqDB::GetGiList` is set), whatever
    /// the list then selects: a list of GIs alone also searches the whole database of version 5.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:87-90
    pub fn has_positive_list(&self) -> bool {
        self.positive_list
    }

    /// The taxids of the search's list object after the volumes opened, ascending, which the
    /// BLAST archive writes (`DbFilterArgs::archive_parameters`): of a positive list those the
    /// database holds (`CSeqDBLMDBSet::TaxIdsToOids`, `tax_ids.swap` with `rv_tax_ids`), of a
    /// negative list all (`CSearchDatabase::x_InitializeDb` opens the database with a copy of
    /// it). Empty without a taxid list.
    pub fn list_taxids(&self) -> &[i32] {
        &self.list_taxids
    }

    /// The GIs of the search's list object after the volumes opened, which the BLAST archive and
    /// search strategy write (`DbFilterArgs::archive_parameters`): sorted ascending (with
    /// duplicates) for a positive `-gilist` (`InsureOrder` ran when the database opened), in file
    /// order (with duplicates) for `-negative_gilist` (`CSearchDatabase::x_InitializeDb` opened
    /// the database with a `CSeqDBIdSet` copy). Empty without a GI list or for a 0-GI file.
    pub fn list_gis(&self) -> &[i64] {
        &self.list_gis
    }

    /// The search path the database was opened with (`open_in`), none for a database opened
    /// by another way.
    pub fn search_path(&self) -> Option<&SearchPath> {
        self.search_path.as_ref()
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:566-572
    /// What `blastdbcmd -info` lists under "Volumes:" for the database that
    /// [`SeqDb::open_in`] opens with the same arguments: the volumes of the first two levels
    /// of its alias tree and the alias files of the third, sorted by base name
    /// (`alias::info_volumes`). The alias files are read again; no volume file is opened.
    pub fn info_volumes(
        name: &[u8],
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<Vec<Vec<u8>>, SeqDbError> {
        alias::info_volumes(name, seqtype, path)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
    /// NCBI's guess of a database's type (`CSeqDB` of `eUnknown`, the programs' `guess`): `open`
    /// of each type of `types` in order. An exception of NCBI's reader
    /// ([`SeqDbError::is_seqdb_exception`]) moves on to the next type, and the last type's error is
    /// the one reported; a refusal of LOSAT's is reported at once, since NCBI may open that
    /// database as the type it tried.
    pub(crate) fn open_guessing<T>(
        types: &[SeqType],
        mut open: impl FnMut(SeqType) -> Result<T, SeqDbError>,
    ) -> Result<T, SeqDbError> {
        let mut types = types.iter().peekable();
        while let Some(&seqtype) = types.next() {
            match open(seqtype) {
                Err(error) if types.peek().is_some() && error.is_seqdb_exception() => {}
                opened => return opened,
            }
        }
        Err(SeqDbError::Internal(
            "no sequence type to open the database as".to_string(),
        ))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:929-936
    /// `FindVolumePaths` (with `paths`, `alias` and true) of the database that [`SeqDb::open_in`]
    /// opens with the same arguments: every volume path (without extension) and every alias file
    /// of its alias tree, each once, in NCBI's order (`SeqDB_CompareVolume`:
    /// [`SeqDb::compare_volume_paths`]). No volume file is opened.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:952-966
    pub fn volume_and_alias_paths(
        name: &[u8],
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>), SeqDbError> {
        let tree = alias::resolve_in(name, seqtype, path)?;
        let volumes: std::collections::BTreeSet<Vec<u8>> = tree.volumes.into_iter().collect();
        let mut alias_files = std::collections::BTreeSet::new();
        alias::collect_alias_files(&tree.root, &mut alias_files);
        Ok((
            Self::sort_volume_set(volumes),
            Self::sort_volume_set(alias_files),
        ))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:929-966
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2629-2668
    /// Single owner of the files that reading the BLAST database `name` (`seqtype`, `path`) can
    /// open: calls `alias::resolve_in` once and collects every alias file, every alias filter list
    /// (`OIDLIST`, `GILIST`, `TILIST`, `SEQIDLIST`, `TAXIDLIST`), every `MASKLIST` GI-mask file,
    /// and every volume component file (`VOLUME_EXTENSIONS`, `LMDB_EXTENSIONS`, and for protein
    /// databases, `rps::RPS_EXTENSIONS`), along with the resolved volume bases in NCBI order.
    pub(crate) fn read_paths(
        name: &[u8],
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<DbReadPaths, SeqDbError> {
        use crate::writedb::build_db::{LMDB_EXTENSIONS, VOLUME_EXTENSIONS};
        let tree = alias::resolve_in(name, seqtype, path)?;
        let volumes = Self::sort_volume_set(tree.volumes.into_iter().collect());
        let mut alias_files = std::collections::BTreeSet::new();
        alias::collect_alias_files(&tree.root, &mut alias_files);

        let mut files = std::collections::BTreeSet::new();
        let mol = seqtype.letter().as_bytes();
        let alias_suffix = match seqtype {
            SeqType::Nucleotide => b".nal".as_slice(),
            SeqType::Protein => b".pal".as_slice(),
        };
        for alias in &alias_files {
            files.insert(alias.clone());
            if let Some(stem) = strip_ascii_suffix(alias, alias_suffix) {
                push_mol_extensions(&mut files, stem, mol, &LMDB_EXTENSIONS);
                files.insert([stem, b".", mol, b"js"].concat());
            }
        }
        alias::collect_alias_list_files(&tree.root, &mut files);
        if let Some(gi_masks) = &tree.gi_masks {
            for mask_name in gi_masks {
                gi_mask::read_paths(mask_name, path, &mut files);
            }
        }
        for volume in &volumes {
            push_mol_extensions(&mut files, volume, mol, &VOLUME_EXTENSIONS);
            push_mol_extensions(&mut files, volume, mol, &LMDB_EXTENSIONS);
            if seqtype == SeqType::Protein {
                for ext in rps::RPS_EXTENSIONS {
                    files.insert([volume.as_slice(), b".", ext.as_bytes()].concat());
                }
            }
        }
        Ok(DbReadPaths {
            volumes,
            files: files.into_iter().collect(),
        })
    }

    fn sort_volume_set(set: std::collections::BTreeSet<Vec<u8>>) -> Vec<Vec<u8>> {
        // `std::sort` of the set's (distinct) members: no two compare equal, so the order is
        // that of the comparison alone.
        let mut paths: Vec<Vec<u8>> = set.into_iter().collect();
        paths.sort_by(|a, b| Self::compare_volume_names(a, b));
        paths
    }

    /// The order in which NCBI searches and numbers the volumes of a database
    /// (`SeqDB_CompareVolume`, the sort at the end of `FindVolumePaths`): by base name (the
    /// text after the last path separator, without an extension `.nal`, `.pal`, `.nin`,
    /// `.pin`, `.ndb` or `.pdb` after at least one other byte), then by the whole path. For
    /// volumes in one directory it is the byte order of the paths.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:190-197
    /// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbgeneral.hpp:975-979
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:76-108
    pub fn compare_volume_paths(a: &Path, b: &Path) -> std::cmp::Ordering {
        Self::compare_volume_names(
            a.as_os_str().as_encoded_bytes(),
            b.as_os_str().as_encoded_bytes(),
        )
    }

    /// [`SeqDb::compare_volume_paths`] of two paths given as their bytes.
    pub(crate) fn compare_volume_names(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
        fn base_name(path: &[u8]) -> &[u8] {
            let separator = |byte: &u8| *byte == b'/' || (cfg!(windows) && *byte == b'\\');
            let name = match path.iter().rposition(separator) {
                Some(at) => &path[at + 1..],
                None => path,
            };
            if name.len() > 4 {
                let (stem, extension) = name.split_at(name.len() - 4);
                if extension[0] == b'.'
                    && matches!(extension[1], b'n' | b'p')
                    && matches!(&extension[2..], b"al" | b"in" | b"db")
                {
                    return stem;
                }
            }
            name
        }
        base_name(a).cmp(base_name(b)).then_with(|| a.cmp(b))
    }

    /// The number of the volume that holds `ordinal` (`CSeqDBVolSet::FindVol`).
    pub fn volume_index(&self, ordinal: u32) -> Result<usize, SeqDbError> {
        self.find_volume(ordinal).map(|at| at.volume)
    }

    /// The number of the volume that holds NCBI's OID `oid` ([`SeqDb::locate_oid`]).
    pub fn oid_volume(&self, oid: u32) -> Result<usize, SeqDbError> {
        self.locate_oid(oid).map(|at| at.volume)
    }

    /// Opens the one database at the path `name` (without extension), taken
    /// as it is: a space or a double quote in it is a character of the path
    /// (`alias::resolve_one_name`). Everything else is as in [`SeqDb::open`],
    /// including NCBI's rule that a volume named in an alias file is looked
    /// for in the current directory first.
    #[cfg(test)]
    pub(crate) fn open_one_name(name: &str, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let tree = alias::resolve_one_name(
            name.as_bytes(),
            seqtype,
            alias::LocalCopy::CurrentDir,
            &SearchPath::current_dir(),
        )?;
        Self::open_tree(tree, name.as_bytes(), seqtype)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1612-1616
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1726
    /// Resolves the one database at the path `name` as `SeqDb::open_one_name` does and reads
    /// the header of each volume's index file as NCBI's constructor of `CSeqDB` reads it
    /// ([`index::IndexFile::read_mapped_header`]: a file that ends early reads as zeros to the end of
    /// its page), without opening the other files of the volumes and without reading the
    /// offset arrays: the cost does not grow with the number of sequences. Returns what
    /// `CSeqDB::GetDBMetaData` reads of the opened database: the title (`GetTitle`: the alias
    /// files' `TITLE`, else the volumes' titles), the date, the volume paths with their totals
    /// and the paths of the alias files (`FindVolumePaths` with `paths`, `alias` and true).
    /// `search` is the search path that the message of a database which is not there names; only
    /// `name` itself is looked for.
    pub(crate) fn locate_one_name_metadata(
        name: &[u8],
        seqtype: SeqType,
        search: &SearchPath,
    ) -> Result<LocatedDatabase, SeqDbError> {
        let tree = alias::resolve_one_name(name, seqtype, alias::LocalCopy::CurrentDir, search)?;
        let mut headers = Vec::with_capacity(tree.volumes.len());
        for volume_name in &tree.volumes {
            let index_path = volume::with_extension(&path_of_bytes(volume_name), seqtype, "in");
            let file = file::DbFile::open(&index_path)?;
            headers.push(index::IndexFile::read_mapped_header(&file, seqtype)?);
        }
        let bases: Vec<(PathBuf, Option<&[u8]>)> = tree
            .volumes
            .iter()
            .zip(&headers)
            .map(|(name, header)| {
                (
                    path_from_bytes(name.as_slice()),
                    header.lmdb_file.as_deref(),
                )
            })
            .collect();
        check_lmdb_names(&bases)?;
        let mut title = Vec::new();
        alias::walk_title(
            &tree.root,
            &|volume_name| {
                tree.volumes
                    .iter()
                    .position(|name| name == volume_name)
                    .map(|at| headers[at].title.clone())
            },
            &mut title,
        );
        let title = fix_string(&title).to_vec();
        let dates: Vec<String> = headers
            .iter()
            .map(|header| String::from_utf8_lossy(fix_string(&header.date)).into_owned())
            .collect();
        let date_refs: Vec<&str> = dates.iter().map(String::as_str).collect();
        let mut alias_files = std::collections::BTreeSet::new();
        alias::collect_alias_files(&tree.root, &mut alias_files);
        Ok(LocatedDatabase {
            title,
            date: date::database_date(&date_refs),
            version5: headers.first().is_some_and(|header| {
                header
                    .lmdb_file
                    .as_ref()
                    .is_some_and(|name| !name.is_empty())
            }),
            volume_totals: headers
                .iter()
                .map(|header| (header.num_oids, header.volume_length))
                .collect(),
            volumes: tree.volumes,
            alias_files: alias_files.into_iter().collect(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.cpp:101-116
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.hpp:503-523
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:984-988
    // The number of sequences and the total length that NCBI takes from its
    // alias walk (seqdbalias.cpp:1189-1192, 1290-1293, 1575-1596) are the sums
    // over the volumes when no alias file has the keys NSEQ and LENGTH, or
    // when their values equal the sums below their nodes (`alias::check_sums`,
    // which rejects every other value).
    /// Opens the volumes of a resolved database (`name` is for messages).
    fn open_tree(
        tree: alias::AliasTree,
        name: &[u8],
        seqtype: SeqType,
    ) -> Result<Self, SeqDbError> {
        let mut volumes = Vec::with_capacity(tree.volumes.len());
        let mut volume_starts = Vec::with_capacity(tree.volumes.len() + 1);
        let mut num_oids: u64 = 0;
        let mut total_length: u64 = 0;
        let mut max_length: u32 = 0;
        let mut protein_max_length: i32 = 0;
        for volume_name in &tree.volumes {
            let volume = Volume::open(&path_of_bytes(volume_name), seqtype)?;
            volume_starts.push(num_oids as u32);
            num_oids += u64::from(volume.num_oids());
            // NCBI holds the number of sequences in a signed 32-bit `int`
            // (seqdbimpl.cpp:948-949).
            if num_oids > i32::MAX as u64 {
                return Err(SeqDbError::Unsupported {
                    path: volume.base().to_path_buf(),
                    message: "the database holds more sequences than a signed 32-bit count"
                        .to_string(),
                });
            }
            total_length = total_length
                .checked_add(volume.index().volume_length)
                .ok_or_else(|| SeqDbError::Format {
                    path: volume.base().to_path_buf(),
                    message: "the total length of the volumes overflows".to_string(),
                })?;
            // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbfile.hpp:557-560
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.hpp:514-523
            // Each volume's `GetMaxLength()` returns `m_MaxLen` as a signed `int`, and the volume
            // set starts its maximum at `0`, so a protein volume whose longest sequence wraps
            // negative (`>= 2^31`) contributes `0` rather than its unsigned length. Nucleotide
            // volumes keep the widened length.
            match seqtype {
                SeqType::Protein => {
                    protein_max_length =
                        protein_max_length.max(ncbi_protein_length(volume.index().max_length));
                    max_length = protein_max_length as u32;
                }
                SeqType::Nucleotide => {
                    max_length = max_length.max(volume.index().max_length);
                }
            }
            volumes.push(volume);
        }
        volume_starts.push(num_oids as u32);

        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:217-233
        let bases: Vec<(PathBuf, Option<&[u8]>)> = volumes
            .iter()
            .map(|volume| {
                (
                    volume.base().to_path_buf(),
                    volume.index().lmdb_file.as_deref(),
                )
            })
            .collect();
        check_lmdb_names(&bases)?;

        let mut title = Vec::new();
        alias::walk_title(
            &tree.root,
            &|volume_name| {
                tree.volumes
                    .iter()
                    .position(|name| name == volume_name)
                    .map(|at| volumes[at].index().title.clone())
            },
            &mut title,
        );
        let title = fix_string(&title).to_vec();

        let dates: Vec<String> = volumes
            .iter()
            .map(|volume| String::from_utf8_lossy(fix_string(&volume.index().date)).into_owned())
            .collect();
        let date_refs: Vec<&str> = dates.iter().map(String::as_str).collect();
        let date = date::database_date(&date_refs).map_err(|error| SeqDbError::Unsupported {
            path: path_from_bytes(name),
            message: format!("the volumes carry different dates and {error}"),
        })?;

        let counts = Counts {
            num_seqs: num_oids as i32,
            total_length,
            max_length,
            min_length: BLAST_SEQSRC_MINLENGTH,
        };
        let mut db = Self {
            seqtype,
            volumes,
            volume_starts,
            title,
            date,
            num_seqs: num_oids as u32,
            counts,
            stats_num_seqs: 0,
            stats_total_length: 0,
            printed_num_seqs: counts.num_seqs,
            printed_total_length: total_length,
            volume_set_length: total_length,
            #[cfg(any(test, feature = "parallel"))]
            included_length: OnceLock::new(),
            first_oid_length: None,
            alias_totals: alias::AliasTotals::default(),
            all_bits_special_set: false,
            oid_list_set_up: AtomicBool::new(false),
            list_applied: false,
            alias_scanned: false,
            alias_range: (0, num_oids as u32),
            oid_list_error: None,
            oid_list_fails_scan: false,
            listed_tree: None,
            unlisted_oid_list_error: None,
            filter: None,
            search_path: None,
            reports_show_gis: false,
            positive_list: false,
            list_taxids: Vec::new(),
            list_gis: Vec::new(),
            name_records: Vec::new(),
            lookup_files: OnceLock::new(),
            lookup_isams: OnceLock::new(),
            mask_algorithms: OnceLock::new(),
            utf8_posted: AtomicBool::new(false),
            scope_utf8: AtomicBool::new(false),
            scope_utf8_ended: AtomicBool::new(false),
            utf8_claimed: AtomicBool::new(false),
            scope_utf8_held: AtomicBool::new(false),
        };
        // A volume list that NCBI opens and whose OIDs it then gets wrong
        // (`SeqDbError::Unsupported`) is refused by the lookups that need them, not here.
        if db.is_version5() {
            match db.lmdb_entries() {
                Ok(_) | Err(SeqDbError::Unsupported { .. } | SeqDbError::NcbiCrashes { .. }) => {}
                Err(error) => return Err(error),
            }
        }
        db.read_alias_tree(tree.root, tree.volumes)?;
        Ok(db)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:94-139
    // `HasFilters` runs `ComputeMasks` (the OID ranges and the paths of the OID lists,
    // `alias::compute_masks`) after the volumes and the LMDB files are open. The walks of a value
    // that does not convert fail where NCBI's constructor reaches them, which depends on the list
    // and on the scan (`check_alias_values`). `x_GetNumSeqs` and `x_GetNumSeqsStats` are the `int`
    // of the walks' `Int8` (seqdbimpl.cpp:948-949); the filters make the OID list
    // (`alias::alias_oids`), which a list later intersects (`restrict`), and the iteration
    // range that its first build sets (`alias::top_iteration_range`), inside which a search reads
    // the list's OIDs. An `OIDLIST` that no node with both `NSEQ` and `LENGTH` covers makes the
    // constructor build the list and scan the totals (`alias::needs_totals_scan`, `scan_alias`):
    // an OID list file that fails is then the constructor's error (probes G1, G3, G6, G7 of line
    // RG D1b, exit status 2), which comes after a list's warnings and errors (`CSeqDBGiListSet`
    // resolves the list before the alias files' OID lists are read; probe S1 of review A of line
    // RG D1b-1) and before the walks of the counts (`oid_list_fails_scan`); else NCBI builds the
    // list at its first use, where that error comes (`build_oid_list`: `-info` never meets it,
    // probe G2a).
    /// Reads the counts, the OID ranges and lists and their errors from the alias tree `root` of
    /// the volumes `names` (the volume paths in the order of `self.volumes`).
    fn read_alias_tree(
        &mut self,
        mut root: alias::AliasNode,
        names: Vec<Vec<u8>>,
    ) -> Result<(), SeqDbError> {
        alias::compute_masks(&mut root)?;
        let tree = ListedTree { root, names };
        let read = self.read_listed_tree(&tree);
        if alias::reads_id_lists(&tree.root, self.is_version5()) {
            self.listed_tree = Some(Box::new(tree));
        }
        read
    }

    /// [`SeqDb::read_alias_tree`] of the alias tree `tree`, after NCBI's `ComputeMasks`.
    fn read_listed_tree(&mut self, tree: &ListedTree) -> Result<(), SeqDbError> {
        let (root, names) = (&tree.root, tree.names.as_slice());
        let at = |path: &[u8]| {
            names
                .iter()
                .position(|name| name == path)
                .expect("every volume of the alias tree is open")
        };
        let totals = alias::walk_totals(root, &|path| {
            let volume = &self.volumes[at(path)];
            (volume.num_oids(), volume.index().volume_length)
        });
        let num_oids = self.num_oids();
        let scan = alias::needs_totals_scan(root);
        let oids = match self.tree_oids(tree, alias::UserList::None, None) {
            Ok(oids) => oids,
            Err(failure) => {
                self.oid_list_error = Some(failure);
                self.oid_list_fails_scan = scan;
                alias::AliasOids {
                    words: None,
                    all_bits_special_set: false,
                    mask_types: Vec::new(),
                    memb_bits: Vec::new(),
                    lists: Vec::new(),
                }
            }
        };
        let oids = self.set_tree_oids(oids);
        // `(int)` of the `Int8` walks.
        let nseq = totals.nseq.as_ref().map_or(0, |&value| value as i32);
        let length = totals.length.as_ref().map_or(0, |&value| value);
        let min = totals.min_seq_length.as_ref().map_or(0, |&value| value);
        self.counts = Counts {
            num_seqs: nseq,
            total_length: length,
            max_length: self.counts.max_length,
            min_length: if min <= 0 {
                BLAST_SEQSRC_MINLENGTH
            } else {
                min
            },
        };
        self.stats_num_seqs = totals.stats_nseq.as_ref().map_or(0, |&value| value as i32);
        self.stats_total_length = totals.stats_totlen.as_ref().map_or(0, |&value| value);
        self.printed_num_seqs = nseq;
        self.printed_total_length = length;
        self.alias_totals = totals;
        if self.oid_list_fails_scan {
            return Ok(());
        }
        if let Some(words) = oids.words {
            let range = alias::top_iteration_range(root, num_oids).unwrap_or((0, num_oids));
            self.alias_range = range;
            let filter = filter::OidFilter::new(words, num_oids).with_range(range);
            if scan {
                return self.scan_alias(filter);
            }
            self.num_seqs = filter.count();
            self.filter = Some(Box::new(filter));
            self.set_included_length(None);
        }
        self.update_first_oid_length()
    }

    /// NCBI's OID list of the alias tree `tree` ([`alias::alias_oids`]) under the search's
    /// Pass 1 positive list `user` and Pass 2 rule `user_rule`.
    fn tree_oids(
        &self,
        tree: &ListedTree,
        user: alias::UserList<'_>,
        user_rule: Option<&filter::DeflineRule>,
    ) -> Result<alias::AliasOids, oid_list::OidListError> {
        let volumes: Vec<alias::TreeVolume<'_>> = tree
            .names
            .iter()
            .enumerate()
            .map(|(at, name)| alias::TreeVolume {
                name: name.as_slice(),
                start: self.volume_starts[at],
                end: self.volume_starts[at + 1],
            })
            .collect();
        alias::alias_oids(
            &tree.root,
            &volumes,
            self.num_oids(),
            self.seqtype,
            self.is_version5(),
            user,
            user_rule,
            &|| self.lmdb_entries(),
        )
    }

    /// Gives the volumes their `OID_MASK_TYPE`, `MEMB_BIT` and the id lists of `oids`, and the
    /// database the shape of NCBI's bit set; returns `oids` without them.
    fn set_tree_oids(&mut self, mut oids: alias::AliasOids) -> alias::AliasOids {
        for (volume, &mask_type) in self.volumes.iter_mut().zip(&oids.mask_types) {
            volume.set_oid_mask_type(mask_type);
        }
        for (volume, &memb_bit) in self.volumes.iter_mut().zip(&oids.memb_bits) {
            volume.set_memb_bit(memb_bit);
        }
        for (at, lists) in std::mem::take(&mut oids.lists).into_iter().enumerate() {
            let start = self.volume_starts[at];
            self.volumes[at].set_alias_lists(lists, start);
        }
        self.all_bits_special_set = oids.all_bits_special_set;
        oids
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:394-466
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:578-716
    // NCBI's `CSeqDB` opened with a list translates the alias files' id lists under it
    // (positive lists in Pass 1; both positive and negative seqid/taxid lists in version 5 Pass 2),
    // when its constructor builds the OID list to scan the totals, after the
    // list's own resolution: the alias files' OID list is made again under the list. Its failure
    // is the constructor's (`build_oid_list`); the failure without the list stays the reports'
    // where their database scans the totals (`unlisted_oid_list_error`).
    /// Makes the alias files' OID list again under the search's Pass 1 list `user` and Pass 2 rule
    /// `user_rule` of a database whose alias tree names id lists (`listed_tree`): the alias files'
    /// bits and iteration range (`filter`), the volumes' id lists and the shape of NCBI's bit set,
    /// or the failure (`oid_list_error`). Nothing changes for any other database.
    fn translate_alias_lists(
        &mut self,
        user: alias::UserList<'_>,
        user_rule: Option<&filter::DeflineRule>,
    ) {
        let Some(tree) = self.listed_tree.take() else {
            return;
        };
        let unlisted = self.oid_list_error.take();
        self.unlisted_oid_list_error = unlisted.filter(|_| self.oid_list_fails_scan);
        match self.tree_oids(&tree, user, user_rule) {
            Ok(oids) => {
                let oids = self.set_tree_oids(oids);
                if let Some(words) = oids.words {
                    let num_oids = self.num_oids();
                    let range =
                        alias::top_iteration_range(&tree.root, num_oids).unwrap_or((0, num_oids));
                    self.alias_range = range;
                    self.filter = Some(Box::new(
                        filter::OidFilter::new(words, num_oids).with_range(range),
                    ));
                }
            }
            Err(failure) => self.oid_list_error = Some(failure),
        }
        self.listed_tree = Some(tree);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-938 (at
    // `statistics_numbers`)
    /// Reads the length of OID 0 when the statistics read it (`first_oid_length`).
    fn update_first_oid_length(&mut self) -> Result<(), SeqDbError> {
        self.first_oid_length = None;
        if self.statistics_total_length() > 0 {
            return Ok(());
        }
        // OID 0 is the first OID of the first volume that has any.
        if let Some(volume) = self.volumes.iter().find(|volume| volume.num_oids() > 0) {
            self.first_oid_length = Some(match self.seqtype {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:932-938
                // `BlastSeqSrcGetSeqLen` returns NCBI's signed 32-bit `GetSeqLengthProt`; any
                // negative length becomes `-1`.
                SeqType::Protein => {
                    let len = volume.protein_length_i32(0)?;
                    if len < 0 {
                        -1
                    } else {
                        i64::from(len)
                    }
                }
                SeqType::Nucleotide => volume.nucleotide_length(0, 0)? as i64,
            });
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:126-133
    // NCBI opens a header file at the first read of a header of its volume
    // (`x_GetHdrAsn1Binary`, seqdbvol.cpp:2302-2310): a search program reads
    // the first one after its search, so a header file that is missing or
    // cannot be opened fails there, not here (`Volume::lazy_file`).
    /// Whether the database holds parsed Seq-ids, which a search program asks right after
    /// opening it: the first header of every volume is decoded with the ordinal-id decoder
    /// (`deflines`), whose error for any other Seq-id is `SeqDbError::Unsupported`
    /// (`SeqDbError::is_parsed_seqids`). `Ok(true)` when some volume's first header holds
    /// parsed Seq-ids; any other error of a first header is returned. A volume whose header
    /// file cannot be opened, mapped or read (missing, an `open` that fails, a first byte that
    /// the system does not give: [`SeqDbError::Map`], or a later read: [`SeqDbError::Read`]) is
    /// not checked: NCBI opens and maps a header file at its first read, and its error, or its
    /// SIGBUS, comes there (after the outfmt 0 prolog; line WQA D18). NCBI reads no header when
    /// it opens a database: a header that its read fails on (`SeqDbError::Ncbi`) fails where a
    /// report reads it, so the volume is classified by its next header (`first_deflines`).
    pub fn has_parsed_seqids(&self) -> Result<bool, SeqDbError> {
        let mut parsed = false;
        for (at, volume) in self.volumes.iter().enumerate() {
            if volume.num_oids() > 0 {
                match self.first_deflines(at) {
                    Ok(_)
                    | Err(
                        SeqDbError::NotFound(_)
                        | SeqDbError::Open { .. }
                        | SeqDbError::Map { .. }
                        | SeqDbError::Read { .. },
                    ) => {}
                    Err(error) if error.is_parsed_seqids() => parsed = true,
                    Err(error) => return Err(error),
                }
            }
        }
        Ok(parsed)
    }

    /// `deflines` of the first sequence of volume `at` (non-empty) whose header NCBI's read does
    /// not fail on, whether a list includes it or not, and with every definition line: the check
    /// reads the header as stored, before a taxid list's `DeflineRule`
    /// (`Volume::stored_deflines`), which would refuse a first sequence that the list leaves out.
    /// Its messages name the sequence's OID. No definition line when NCBI's read fails on every
    /// header of the volume (each report of them fails).
    fn first_deflines(&self, at: usize) -> Result<Vec<OrdinalDefline>, SeqDbError> {
        let volume = &self.volumes[at];
        for vol_oid in 0..volume.num_oids() {
            let oid = self.volume_starts[at] + vol_oid;
            match self.ordinal_deflines(at, oid, volume.stored_deflines(vol_oid)) {
                Err(SeqDbError::Ncbi(_)) => {}
                other => return other,
            }
        }
        Ok(Vec::new())
    }

    /// Whether the database is of format version 5 (its volumes' index files name an LMDB
    /// file; `open` refuses a mix of versions).
    pub fn is_version5(&self) -> bool {
        self.volumes.first().is_some_and(|volume| {
            volume
                .index()
                .lmdb_file
                .as_ref()
                .is_some_and(|name| !name.is_empty())
        })
    }

    /// Sequence type of the database.
    pub fn seqtype(&self) -> SeqType {
        self.seqtype
    }

    /// The number of included OIDs, the sequences a search reads: ordinal ids run from 0 to
    /// this number minus one. They are the OIDs of the volumes that the alias files' ranges
    /// (`FIRST_OID`, `LAST_OID`) and a list ([`SeqDb::restrict`]) include, inside the iteration
    /// range NCBI sets for the only alias file of a `-db` (without a list); iterations and
    /// schedules count them. NCBI's `GetNumSeqs` is another number, which the statistics
    /// ([`SeqDb::statistics_numbers`]) and the warning of [`SeqDb::filtered_empty`] read.
    pub fn num_seqs(&self) -> u32 {
        self.num_seqs
    }

    /// The OIDs of the volumes (`CSeqDB::GetNumOIDs`), included or not: the volume set holds
    /// no sequence when it is 0.
    pub fn num_oids(&self) -> u32 {
        self.volume_starts[self.volumes.len()]
    }

    // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp:44-49
    /// Every OID of the volumes as an ordinal id ([`SeqDb::num_oids`] of them): the alias
    /// files' OID ranges (`FIRST_OID`, `LAST_OID`) no longer leave any out. For a reader that
    /// loops `GetNumOIDs` and reads each OID with `GetBioseq` of `oid`, which reads an OID outside
    /// the ranges too (dustmasker's `CMaskBDBReader`). Only before a list
    /// ([`SeqDb::restrict`]): the alias keys of lists are refused when the database opens.
    fn include_every_oid(&mut self) {
        debug_assert!(!self.list_applied, "include_every_oid after a list");
        self.filter = None;
        self.set_included_length(None);
        self.num_seqs = self.num_oids();
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:117-139
    /// Opens the database `name` for a reader that loops `GetNumOIDs` and reads each OID with
    /// `GetBioseq` of `oid` (dustmasker's `CMaskBDBReader`, which makes `CSeqDB` from `name` and
    /// `type`): [`SeqDb::open_search_in`], the errors of the walks of the alias files' counts that
    /// NCBI's constructor makes (those of [`SeqDb::check_alias_values`] without a list), and
    /// `include_every_oid`. LOSAT's refusals of the values that only a search reads
    /// differently are not made, as such a reader reads none of them: a negative `NSEQ`, a
    /// `LENGTH` of 0 or of 2^63 or more (NCBI 2.17.0 dustmasker reads them as it reads the
    /// database without the alias file: review 2 of line RG item K1, F5).
    pub fn open_every_oid_in(
        name: &[u8],
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<Self, SeqDbError> {
        let mut db = Self::open_search_in(name, seqtype, path)?;
        db.check_alias_values()?;
        db.include_every_oid();
        Ok(db)
    }

    /// `CSeqDB::GetNumSeqs` of the search: the walk of `NSEQ` (the volumes' OIDs below the
    /// alias files without the key) without a list, the included OIDs with one. It is not the
    /// number of included OIDs ([`SeqDb::num_seqs`]) when an alias file has `NSEQ` or a range.
    pub fn ncbi_num_seqs(&self) -> i32 {
        self.counts.num_seqs
    }

    /// `CSeqDB::GetTotalLength` of the search: the walk of `LENGTH` (the volumes' lengths below
    /// the alias files without the key) without a list, NCBI's scanned total of the included
    /// OIDs with one (approximate for nucleotides: `restrict_taxids`).
    pub fn total_length(&self) -> u64 {
        self.counts.total_length
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:872-876
    /// `CSeqDB::GetVolumeLength`: the total length of the volumes, whatever the alias files and a
    /// list include (blastdbcheck's metadata test).
    pub fn volume_length(&self) -> u64 {
        self.volume_set_length
    }

    /// `CSeqDB::GetMaxLength`: the longest length of the volumes (an alias file's
    /// `MAX_SEQ_LENGTH` is never read), with a list of NCBI's scan as `total_length`.
    pub fn max_length(&self) -> u32 {
        self.counts.max_length
    }

    /// `CSeqDB::GetMinLength`: without a list the walk of `MIN_SEQ_LENGTH` (0 for a branch of
    /// volumes alone: NCBI's index files store no shortest length), `BLAST_SEQSRC_MINLENGTH`
    /// when that is not positive; with a list the scanned shortest length of the included OIDs
    /// (`INT4_MAX` without any), which the preliminary cutoff of the protein programs reads
    /// (`BLAST_GapAlignSetUp`, `min_subject_length`).
    pub fn min_length(&self) -> i32 {
        self.counts.min_length
    }

    /// `CSeqDB::GetNumSeqsStats`: the walk of `STATS_NSEQ`, 0 without the key.
    pub fn num_seqs_stats(&self) -> i32 {
        self.stats_num_seqs
    }

    /// `CSeqDB::GetTotalLengthStats`: the walk of `STATS_TOTLEN`, 0 without the key.
    pub fn total_length_stats(&self) -> u64 {
        self.stats_total_length
    }

    // NCBI's reports open the database again without the list (`CBlastFormatUtil::
    // GetBlastDbInfo` reads a `CSeqDB` of the name alone): the `Database:` lines of the
    // pairwise report and the database counts of XML and JSON are the walks of the alias tree.
    /// The number of sequences the reports print: `GetNumSeqs` of the database without a list,
    /// NCBI's `int` (negative for a negative `NSEQ` or one of 2^31 or more).
    pub fn printed_num_seqs(&self) -> i32 {
        self.printed_num_seqs
    }

    /// The total length the reports print: `GetTotalLength` of the database without a list.
    pub fn printed_total_length(&self) -> u64 {
        self.printed_total_length
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:206-215
    /// Whether NCBI warns that the database is empty (`GetNumSeqs` of the search is 0:
    /// `report::query_warnings::EMPTY_DATABASE_MESSAGE` for each query): a list that kept no
    /// sequence of a database that has some, an alias file's `NSEQ 0`, or its `LENGTH` or
    /// `STATS_TOTLEN` over volumes without sequences. The search goes on over the included OIDs,
    /// which may be none. Where the first query batch fails instead
    /// ([`SeqDb::first_oid_read_fails`]), NCBI writes no warning.
    pub fn filtered_empty(&self) -> bool {
        self.counts.num_seqs == 0 && !self.first_oid_read_fails()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:905-938
    // `BlastSeqSrcGetTotLenStats` and `BlastSeqSrcGetTotLen` are the `Int8` of NCBI's `Uint8`
    // (seqsrc_seqdb.cpp:183-197). `CSeqDB::GetSeqLength` of OID 0 does not ask the list or the
    // ranges.
    /// The total length the statistics read before the OID 0 fallback: `STATS_TOTLEN`, else
    /// `GetTotalLength`, as NCBI's `Int8`.
    fn statistics_total_length(&self) -> i64 {
        let stats = self.stats_total_length as i64;
        if stats > 0 {
            stats
        } else {
            self.counts.total_length as i64
        }
    }

    /// The numbers of the database that NCBI's statistics read, with NCBI's types
    /// (`DbStatistics`): the setup's total and count (`BLAST_GapAlignSetUp`), the Gumbel block's
    /// length, the engine's `BlastSeqSrcGetTotLen`, the shortest and the average length.
    pub fn statistics(&self) -> DbStatistics {
        let total = self.statistics_total_length();
        let request_count = if self.stats_num_seqs > 0 {
            self.stats_num_seqs
        } else {
            self.counts.num_seqs
        };
        let (setup_total, setup_count) = if total > 0 {
            (total, request_count)
        } else {
            match self.first_oid_length {
                Some(length) => (length, 1),
                None => (0, 0),
            }
        };
        DbStatistics {
            setup_total,
            setup_count,
            gumbel_db_length: total,
            engine_db_length: self.counts.total_length as i64,
            request_count,
            min_length: self.counts.min_length,
            average: self.average_length(),
        }
    }

    /// The total length and the number of sequences that the statistics of a search read
    /// (`BLAST_CalcEffLengths`, `BlastSetup_GetScoreBlock`): `DbStatistics::setup_total` and
    /// `setup_count`, NCBI's `Int8` and `Int4`. The number is 0 for an alias file's `NSEQ 0`,
    /// which NCBI's average subject length divides by (`common::database::setup_count_refused`), and
    /// negative for a negative `NSEQ` or one of 2^31 or more.
    pub fn statistics_numbers(&self) -> (i64, i32) {
        let statistics = self.statistics();
        (statistics.setup_total, statistics.setup_count)
    }

    /// Whether every OID of the volumes is included and is its own ordinal id: no alias file's OID
    /// range (`FIRST_OID`, `LAST_OID`) and no list (`SeqDb::included_ordinal`).
    pub fn ordinal_ids_are_oids(&self) -> bool {
        self.filter.is_none()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:164-173
    /// Whether the total length of the statistics is not positive before their fallback to the
    /// length of OID 0 (`STATS_TOTLEN` and `GetTotalLength` are not: a list that kept no
    /// sequence, without `STATS_TOTLEN`, or of 2^63 or more): with no OID either, the first query
    /// batch's setup fails (`first_oid_read_fails`).
    pub fn statistics_length_is_empty(&self) -> bool {
        self.statistics().gumbel_db_length <= 0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:932-936 (at
    // `statistics_numbers`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251 (at
    // `common::database::PreliminaryReadError`): the `CSeqDBException` of that read is
    // `BLAST engine error: Database memory map file error` (exit status 3; psiblast's
    // `BLAST Database error: OID not found`, exit status 2).
    /// Whether the setup of the first query batch with a valid context fails: the statistics'
    /// total is not positive ([`SeqDb::statistics_length_is_empty`]), so NCBI reads the length of
    /// OID 0 (`GetSeqLength` of 0, which asks the volumes, not the list), and the volumes have no
    /// OID. An alias file's `LENGTH` or `STATS_TOTLEN` over such volumes makes the total positive:
    /// NCBI then reads no sequence and reports no hit (probes p4c_dp, p4c_dn of line RG D1a-2).
    pub fn first_oid_read_fails(&self) -> bool {
        self.statistics_length_is_empty() && self.num_oids() == 0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:199-209
    /// `BlastSeqSrcGetAvgSeqLen`: `GetTotalLength` divided by the larger of 1 and `GetNumSeqs` (the
    /// alias files' `STATS_*` keys are not read), as NCBI's `Int4`.
    pub fn average_length(&self) -> i32 {
        let total = self.counts.total_length as i64;
        (total / i64::from(self.counts.num_seqs.max(1))) as i32
    }

    /// The errors and LOSAT's refusals of the alias files' values that NCBI's `CSeqDB` of the
    /// search meets (`read_alias_tree`): without a list the walks of `NSEQ`, `LENGTH`,
    /// `MIN_SEQ_LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN` in that order; after a list
    /// ([`SeqDb::restrict`]) or with an `OIDLIST` that makes NCBI scan the totals
    /// (`alias::needs_totals_scan`), those of `STATS_NSEQ` and `STATS_TOTLEN` (the reports'
    /// database, opened again without the list, reads the others unless it scans too:
    /// [`SeqDb::read_name_records`]). NCBI's error is `SeqDbError::Ncbi`. Before them, an alias
    /// file's OID list file that fails where NCBI's constructor scans the totals (its build of
    /// the OID list, [`SeqDb::build_oid_list`]). [`SeqDb::open_in`] calls this; a search program
    /// opens with [`SeqDb::open_search_in`] and calls it after its list.
    pub fn check_alias_values(&self) -> Result<(), SeqDbError> {
        use alias::WalkKey;
        if self.oid_list_fails_scan {
            self.build_oid_list()?;
        }
        let keys: &[WalkKey] = if self.list_applied || self.alias_scanned {
            &[WalkKey::StatsNseq, WalkKey::StatsTotlen]
        } else {
            UNLISTED_WALK_KEYS
        };
        self.walk_errors(keys)
    }

    /// The first exception of the walks of `keys`, in that order.
    fn walk_errors(&self, keys: &[alias::WalkKey]) -> Result<(), SeqDbError> {
        match keys.iter().find_map(|&key| self.alias_totals.error(key)) {
            Some(exception) => Err(SeqDbError::Ncbi(exception.clone())),
            None => Ok(()),
        }
    }

    /// Title of the database.
    pub fn title(&self) -> &[u8] {
        &self.title
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1010-1013
    /// Date of the database: the latest date of its volumes (`date`).
    pub fn date(&self) -> &str {
        &self.date
    }

    /// The volumes, in the order of their ordinal ids.
    pub fn volumes(&self) -> &[Volume] {
        &self.volumes
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.hpp:245-275
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:416-425
    /// The place of NCBI's OID `oid` (`FindVol`), whether the alias files' ranges and a list
    /// include it or not: NCBI's reads by OID do not ask them. `OID not found` (`kOidNotFound`,
    /// seqdb.cpp:53) for an OID past the volumes.
    pub fn locate_oid(&self, oid: u32) -> Result<SeqLocation, SeqDbError> {
        if oid >= self.num_oids() {
            return Err(SeqDbError::Range("OID not found".to_string()));
        }
        // The last volume whose first OID is not above `oid`; empty volumes share their
        // start with the next volume.
        let volume = self.volume_starts[..self.volumes.len()].partition_point(|&s| s <= oid) - 1;
        Ok(SeqLocation {
            volume,
            vol_oid: oid - self.volume_starts[volume],
            id: oid,
        })
    }

    /// The place of the ordinal id `ordinal`: the place of its OID ([`SeqDb::locate_oid`]),
    /// named by `ordinal`.
    fn find_volume(&self, ordinal: u32) -> Result<SeqLocation, SeqDbError> {
        if ordinal >= self.num_seqs {
            return Err(SeqDbError::Range(format!(
                "ordinal id {ordinal} is outside the database ({} sequences)",
                self.num_seqs
            )));
        }
        let at = self.locate_oid(self.oid(ordinal).get())?;
        Ok(SeqLocation { id: ordinal, ..at })
    }

    /// The ordinal id of the OID `oid` (`SeqDb::locate_oid`), `None` where the alias files' ranges
    /// or a list leave it out (NCBI's `x_CheckOrFindOID` of a lookup by Seq-id, which then finds
    /// nothing: seqdbimpl.cpp:1414-1418) or past the volumes. RPS-BLAST searches every profile of
    /// its volumes by OID and its reports look them up by ordinal id
    /// (`rps_search::volumes::report_ordinal_ids`).
    pub fn included_ordinal(&self, oid: u32) -> Option<u32> {
        match &self.filter {
            None => (oid < self.num_oids()).then_some(oid),
            Some(filter) => filter.contains(oid).then(|| filter.rank(oid)),
        }
    }

    /// The OID in the database of the sequence at `at`, which the readers' messages name (the id
    /// that `at` was found by can be an ordinal id).
    fn location_oid(&self, at: &SeqLocation) -> u32 {
        self.volume_starts[at.volume] + at.vol_oid
    }

    /// `end` cut at the first OID of `volume` from `first` on that a taxid list excludes:
    /// a block of a restricted database holds consecutive OIDs.
    fn included_end(&self, volume: usize, first: u32, end: u32) -> u32 {
        match &self.filter {
            Some(filter) => {
                let start = self.volume_starts[volume];
                filter.run_end(start + first, start + end) - start
            }
            None => end,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:328-350
    // NCBI hands its search threads chunks of consecutive ordinal ids of a
    // fixed count. A LOSAT block is such a chunk cut by bytes instead: it
    // ends before the sequence that would pass `max_bytes`, and at the end of
    // its volume, because one block is one read of one file. The cut does not
    // change which sequences are searched or their order.
    /// The block that starts at `first_ordinal`: as many following sequences
    /// of the same volume as fit `max_bytes` of the sequence file, at least
    /// one. `None` when `first_ordinal` is past the last sequence.
    ///
    /// When the index file cannot be read here (`LOSAT_X_DBIDXREAD`), the
    /// block is the one sequence `first_ordinal`: reading it
    /// ([`SeqDb::read_block`]) reads the same index entries and returns the
    /// error where a block's read errors are reported.
    pub fn next_block(&self, first_ordinal: u32, max_bytes: u64) -> Option<BlockRange> {
        let SeqLocation {
            volume,
            vol_oid: first,
            ..
        } = self.find_volume(first_ordinal).ok()?;
        let end = self.volumes[volume]
            .block_end(first, max_bytes)
            .unwrap_or(first + 1);
        let end = self.included_end(volume, first, end);
        Some(BlockRange {
            volume,
            first,
            end,
            first_ordinal,
        })
    }

    /// EXPERIMENT (LOSAT_X_BNDBSMALL): the blocks that `next_block` gives from `first_ordinal`
    /// on, one after the other, up to `max` of them. The sequence offsets of a volume are read
    /// `X_LIST_ENTRIES` at a time (64 KiB) and kept for the next blocks, where `next_block` reads
    /// two entries and 1,024 more from the file for each block (`LOSAT_X_DBIDXREAD`; with the
    /// arrays held, both cut in memory): listing the 81,218 blocks of 256 KiB of
    /// ref_prok_rep_genomes took 183 ms, on one thread while the others waited, at every query
    /// batch (`v4-D12/STATE.md`, s30). The blocks are those of `next_block` (the same ends,
    /// `Volume::x_listed_block_end`, and the same cut at excluded OIDs).
    #[cfg(any(
        test,
        all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )
    ))]
    pub fn x_next_blocks(&self, first_ordinal: u32, max_bytes: u64, max: usize) -> Vec<BlockRange> {
        self.x_next_blocks_in(first_ordinal, max_bytes, max, X_LIST_ENTRIES)
    }

    /// `x_next_blocks` with windows of `entries` offsets.
    #[cfg(any(
        test,
        all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )
    ))]
    fn x_next_blocks_in(
        &self,
        first_ordinal: u32,
        max_bytes: u64,
        max: usize,
        entries: usize,
    ) -> Vec<BlockRange> {
        let mut blocks = Vec::new();
        let mut window = volume::XOffsetWindow::default();
        let mut window_volume = usize::MAX;
        let mut next_ordinal = first_ordinal;
        while blocks.len() < max {
            let Ok(SeqLocation {
                volume,
                vol_oid: first,
                ..
            }) = self.find_volume(next_ordinal)
            else {
                break;
            };
            if volume != window_volume {
                window = volume::XOffsetWindow::default();
                window_volume = volume;
            }
            let end = self.volumes[volume]
                .x_listed_block_end(first, max_bytes, &mut window, entries)
                .unwrap_or(first + 1);
            let end = self.included_end(volume, first, end);
            let range = BlockRange {
                volume,
                first,
                end,
                first_ordinal: next_ordinal,
            };
            blocks.push(range);
            next_ordinal = range.end_ordinal();
        }
        blocks
    }

    /// The block of the sequences `first_ordinal..end_ordinal`
    /// that the volume of `first_ordinal` holds (to `end_ordinal` or to the end of the volume),
    /// found without reading the index (`next_block` reads it to cut by bytes). `None` when the
    /// range is empty or `first_ordinal` is past the last sequence.
    pub fn x_block_range(&self, first_ordinal: u32, end_ordinal: u32) -> Option<BlockRange> {
        if first_ordinal >= end_ordinal {
            return None;
        }
        let SeqLocation {
            volume,
            vol_oid: first,
            ..
        } = self.find_volume(first_ordinal).ok()?;
        let end = self.volumes[volume]
            .num_oids()
            .min(first.saturating_add(end_ordinal - first_ordinal));
        let end = self.included_end(volume, first, end);
        Some(BlockRange {
            volume,
            first,
            end,
            first_ordinal,
        })
    }

    /// The numbers of residues of the proteins of
    /// `x_block_range(first_ordinal, end_ordinal)`, from the index alone
    /// (`Volume::x_protein_lengths`).
    pub fn x_protein_lengths(
        &self,
        first_ordinal: u32,
        end_ordinal: u32,
    ) -> Result<Vec<u32>, SeqDbError> {
        let Some(range) = self.x_block_range(first_ordinal, end_ordinal) else {
            return Err(SeqDbError::Range(format!(
                "ordinal ids {first_ordinal}..{end_ordinal} are not a non-empty range of the database ({} sequences)",
                self.num_seqs()
            )));
        };
        self.volumes[range.volume].x_protein_lengths(range.first, range.end)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
    /// [`Self::x_protein_lengths`] as NCBI's signed 32-bit `GetSeqLengthProt` lengths
    /// (`Volume::x_protein_lengths_i32`).
    pub fn x_protein_lengths_i32(
        &self,
        first_ordinal: u32,
        end_ordinal: u32,
    ) -> Result<Vec<i32>, SeqDbError> {
        let Some(range) = self.x_block_range(first_ordinal, end_ordinal) else {
            return Err(SeqDbError::Range(format!(
                "ordinal ids {first_ordinal}..{end_ordinal} are not a non-empty range of the database ({} sequences)",
                self.num_seqs()
            )));
        };
        self.volumes[range.volume].x_protein_lengths_i32(range.first, range.end)
    }

    /// Reads the sequences of `range` with one positioned read.
    pub fn read_block(&self, range: BlockRange) -> Result<SequenceBlock, SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        volume.read_block(
            range.first,
            range.end,
            range.first_ordinal,
            self.volume_starts[range.volume] + range.first,
        )
    }

    fn volume_and_first_oid(&self, range: BlockRange) -> Result<(&Volume, u32), SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        Ok((volume, self.volume_starts[range.volume] + range.first))
    }

    /// Reads the protein sequences of `range` for a search (`Volume::read_protein_search_block`),
    /// validating `.psq` ranges without allocating bytes for proteins whose NCBI 32-bit length is
    /// not positive (`<= 0`).
    pub fn read_protein_search_block(
        &self,
        range: BlockRange,
    ) -> Result<SequenceBlock, SeqDbError> {
        let (volume, first_oid) = self.volume_and_first_oid(range)?;
        volume.read_protein_search_block(range.first, range.end, range.first_ordinal, first_oid)
    }

    /// Reads the one protein sequence `ordinal` for a search, skipping byte reads when its NCBI
    /// 32-bit length is not positive (`Volume::read_protein_search_block`).
    pub fn read_protein_search_sequence(&self, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].read_protein_search_block(
            at.vol_oid,
            at.vol_oid + 1,
            at.id,
            self.location_oid(&at),
        )
    }

    /// The length and the ambiguity data that
    /// `read_block(range)` gives the one nucleotide of `range`, read without its 2-bit data, with
    /// the errors of `read_block` (`Volume::x_nucleotide_check`).
    pub fn x_nucleotide_check(
        &self,
        range: BlockRange,
    ) -> Result<(usize, SequenceBlock), SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        volume.x_nucleotide_check(
            range.first,
            range.end,
            range.first_ordinal,
            self.volume_starts[range.volume] + range.first,
        )
    }

    /// Reads the one sequence `ordinal` (for traceback and reports).
    pub fn read_sequence(&self, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
        self.read_sequence_at(self.find_volume(ordinal)?)
    }

    /// Reads the one sequence at `at` (a block of one sequence whose first ordinal id is the id
    /// `at` was found by).
    pub fn read_sequence_at(&self, at: SeqLocation) -> Result<SequenceBlock, SeqDbError> {
        self.volumes[at.volume].read_block(
            at.vol_oid,
            at.vol_oid + 1,
            at.id,
            self.location_oid(&at),
        )
    }

    /// The number of residues of protein `ordinal`, read from the index
    /// alone (`read_sequence(ordinal)`'s `protein(0).len()`).
    pub fn protein_length(&self, ordinal: u32) -> Result<u32, SeqDbError> {
        self.protein_length_at(self.find_volume(ordinal)?)
    }

    /// [`SeqDb::protein_length`] of the protein at `at`.
    pub fn protein_length_at(&self, at: SeqLocation) -> Result<u32, SeqDbError> {
        self.volumes[at.volume].protein_length(at.vol_oid)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
    /// The signed 32-bit length of protein `ordinal` as NCBI's `GetSeqLengthProt` returns it
    /// (`Volume::protein_length_i32`).
    pub fn protein_length_i32(&self, ordinal: u32) -> Result<i32, SeqDbError> {
        self.protein_length_i32_at(self.find_volume(ordinal)?)
    }

    /// [`SeqDb::protein_length_i32`] of the protein at `at`.
    pub fn protein_length_i32_at(&self, at: SeqLocation) -> Result<i32, SeqDbError> {
        self.volumes[at.volume].protein_length_i32(at.vol_oid)
    }

    /// Validates that the sequence file opens and covers protein `ordinal` without reading its
    /// bytes (`Volume::check_protein_range`).
    pub fn check_protein_range(&self, ordinal: u32) -> Result<(), SeqDbError> {
        self.check_protein_range_at(self.find_volume(ordinal)?)
    }

    /// [`SeqDb::check_protein_range`] of the protein at `at`.
    pub fn check_protein_range_at(&self, at: SeqLocation) -> Result<(), SeqDbError> {
        self.volumes[at.volume].check_protein_range(at.vol_oid)
    }

    /// Reads `len` NCBIstdaa bytes from residue `offset` of the protein at `at` without reading
    /// the whole sequence (`Volume::read_protein_slice`).
    pub fn read_protein_slice_at(
        &self,
        at: SeqLocation,
        offset: u64,
        len: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        self.volumes[at.volume].read_protein_slice(at.vol_oid, offset, len)
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): `read_block` with the bytes in the
    /// buffer that a worker keeps from one block to the next, `kept` (`Volume::x_read_block_into`;
    /// a block of 1 MiB or more that `kept` does not hold, in a buffer of `utils::x_mapped::X_BYTES`);
    /// `SequenceBlock::x_into_kept` gives the buffer back.
    pub fn x_read_block_into(
        &self,
        range: BlockRange,
        kept: &mut Vec<u8>,
    ) -> Result<SequenceBlock, SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        volume.x_read_block_into(
            range.first,
            range.end,
            range.first_ordinal,
            self.volume_starts[range.volume] + range.first,
            kept,
        )
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP, LOSAT_X_BNMAPPED): [`Self::read_protein_search_block`] into
    /// the worker's buffer `kept` (`Volume::x_read_protein_search_block_into`).
    pub fn x_read_protein_search_block_into(
        &self,
        range: BlockRange,
        kept: &mut Vec<u8>,
    ) -> Result<SequenceBlock, SeqDbError> {
        let (volume, first_oid) = self.volume_and_first_oid(range)?;
        volume.x_read_protein_search_block_into(
            range.first,
            range.end,
            range.first_ordinal,
            first_oid,
            kept,
        )
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `read_sequence` with the bytes in a buffer that glibc maps
    /// from 1 MiB on.
    pub fn x_read_sequence_mapped(&self, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].x_read_block_mapped(
            at.vol_oid,
            at.vol_oid + 1,
            ordinal,
            self.location_oid(&at),
        )
    }

    /// The length in bases of nucleotide `ordinal`, read from the index and the last byte of
    /// its packed bases (no letters).
    pub fn nucleotide_length(&self, ordinal: u32) -> Result<usize, SeqDbError> {
        self.nucleotide_length_at(self.find_volume(ordinal)?)
    }

    /// [`SeqDb::nucleotide_length`] of the nucleotide at `at`.
    pub fn nucleotide_length_at(&self, at: SeqLocation) -> Result<usize, SeqDbError> {
        self.volumes[at.volume].nucleotide_length(at.vol_oid, self.location_oid(&at))
    }

    /// [`SeqDb::nucleotide_iupac_range_at`] of nucleotide `ordinal`.
    #[cfg(test)]
    pub fn nucleotide_iupac_range(
        &self,
        ordinal: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        self.nucleotide_iupac_range_at(self.find_volume(ordinal)?, from, to)
    }

    /// The IUPAC letters `[from, to)` of the nucleotide at `at`: the bases of the range and the
    /// ambiguity runs that overlap it, read without the rest of the sequence.
    pub fn nucleotide_iupac_range_at(
        &self,
        at: SeqLocation,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        self.volumes[at.volume].nucleotide_iupac_range(at.vol_oid, self.location_oid(&at), from, to)
    }

    /// [`SeqDb::nucleotide_iupac_at`] of nucleotide `ordinal`.
    pub fn nucleotide_iupac(&self, ordinal: u32) -> Result<Vec<u8>, SeqDbError> {
        self.nucleotide_iupac_at(self.find_volume(ordinal)?)
    }

    /// The IUPAC letters of the whole nucleotide at `at`: [`SeqDb::nucleotide_iupac_range_at`]
    /// of `[0, length)`, the length read once, with the sequence's ambiguity data
    /// (`Volume::nucleotide_iupac`).
    pub fn nucleotide_iupac_at(&self, at: SeqLocation) -> Result<Vec<u8>, SeqDbError> {
        self.volumes[at.volume].nucleotide_iupac(at.vol_oid, self.location_oid(&at))
    }

    /// The length and the ambiguity runs of nucleotide `ordinal`, read once for any number of
    /// [`SeqDb::nucleotide_iupac_range_with`] calls on it.
    pub fn nucleotide_ambiguity(&self, ordinal: u32) -> Result<NucleotideAmbiguity, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].nucleotide_ambiguity(at.vol_oid, self.location_oid(&at))
    }

    /// [`SeqDb::nucleotide_iupac_range_at`] of nucleotide `ordinal` with the length and the runs
    /// of [`SeqDb::nucleotide_ambiguity`] for the same `ordinal`: one read of the range's bytes.
    pub fn nucleotide_iupac_range_with(
        &self,
        ordinal: u32,
        from: usize,
        to: usize,
        ambiguity: &NucleotideAmbiguity,
    ) -> Result<Vec<u8>, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].nucleotide_iupac_range_with(
            at.vol_oid,
            self.location_oid(&at),
            from,
            to,
            ambiguity,
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2222-2241
    /// The definition lines of `ordinal`: the title, `n` of `gnl|BL_ORD_ID|n`, where `n` is the
    /// stored number plus the first OID of the volume, as NCBI reports it, and the taxid (a
    /// record whose taxid does not decode is an `Err`, as in `defline_taxids`).
    pub fn deflines(&self, ordinal: u32) -> Result<Vec<OrdinalDefline>, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.ordinal_deflines(
            at.volume,
            self.location_oid(&at),
            self.volumes[at.volume].deflines(at.vol_oid),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2221-2240
    /// The header of `ordinal` as a report reads it: the ordinal-id definition lines
    /// (`deflines`, with their taxids), or every member and Seq-id (`defline_set`) of a record
    /// with parsed Seq-ids.
    pub fn read_header(&self, ordinal: u32) -> Result<SequenceHeader, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.sequence_header(
            at.volume,
            self.location_oid(&at),
            self.volumes[at.volume].header(at.vol_oid),
        )
    }

    /// The header of the OID `oid` (`locate_oid`), whether the alias files' ranges and a list
    /// include it or not, as NCBI's `GetSeqIDs` reads it (seqdbimpl.cpp:793-813, cited at
    /// `rps_search::volumes::report_ordinal_ids`): RPS-BLAST's Seq-aligns name every profile of a
    /// volume's results by it.
    pub fn read_header_of_oid(&self, oid: u32) -> Result<SequenceHeader, SeqDbError> {
        let at = self.locate_oid(oid)?;
        self.sequence_header(at.volume, oid, self.volumes[at.volume].header(at.vol_oid))
    }

    /// [`SeqDb::read_header_of_oid`] with every definition line the header stores
    /// (`Volume::stored_header`): the read of a volume that NCBI's RPS search of a database of
    /// several volumes searches bare, where neither a list nor the alias files' `OID_MASK_TYPE`
    /// reaches the Seq-aligns' ids.
    pub fn read_stored_header_of_oid(&self, oid: u32) -> Result<SequenceHeader, SeqDbError> {
        let at = self.locate_oid(oid)?;
        self.sequence_header(
            at.volume,
            oid,
            self.volumes[at.volume].stored_header(at.vol_oid),
        )
    }

    /// The `SequenceHeader` of the header `header` of the sequence `oid` (its OID, which the
    /// messages name) in `volume`.
    fn sequence_header(
        &self,
        volume: usize,
        oid: u32,
        header: Result<VolumeHeader, SeqDbError>,
    ) -> Result<SequenceHeader, SeqDbError> {
        match header? {
            VolumeHeader::Ordinal(deflines) => self
                .ordinal_deflines(volume, oid, Ok(deflines))
                .map(SequenceHeader::Ordinal),
            VolumeHeader::Parsed(set) => self
                .ordinal_set(volume, oid, set)
                .map(SequenceHeader::Parsed),
        }
    }

    /// The definition lines of `ordinal` with every member and every Seq-id
    /// (`seqid::ber::decode_full_set`), as NCBI's `GetHdr` gives them: the
    /// tag of a `general` id `BL_ORD_ID` is its OID in the database (the
    /// stored number plus the first OID of the volume).
    pub fn defline_set(&self, ordinal: u32) -> Result<DeflineSet, SeqDbError> {
        self.defline_set_at(self.find_volume(ordinal)?)
    }

    /// [`SeqDb::defline_set`] written over `set`, whose buffers are kept from one sequence to the
    /// next (`seqid::ber::decode_full_set_in`: sequential readers such as `seqmasks_io`).
    pub fn defline_set_in(&self, ordinal: u32, set: &mut DeflineSet) -> Result<(), SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.header_set_in(at, set)?;
        Ok(())
    }

    /// The member that NCBI's warning of a UTF8String read where a VisibleString is expected
    /// names (`seqid::ber::utf8_string_path`) for the header of `ordinal`: where NCBI's read of
    /// the record meets one before any failure. `None` also where the record cannot be read.
    /// The record is read again: for a read that noted such a string (`OrdinalDefline::utf8`,
    /// `DeflineSet::utf8`) or failed.
    #[cold]
    pub fn utf8_string_path(&self, ordinal: u32) -> Option<String> {
        self.utf8_string_path_at(self.find_volume(ordinal).ok()?)
    }

    // NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:617-622
    /// Whether NCBI's warning of a UTF8String read for a VisibleString is to be posted now:
    /// true once, at the first call (`ERR_POST_X_ONCE` posts once per process; a search opens
    /// its database once, so once per run of a program).
    pub fn take_utf8_post(&self) -> bool {
        !self.utf8_claimed.load(std::sync::atomic::Ordering::Relaxed)
            && !self
                .utf8_posted
                .swap(true, std::sync::atomic::Ordering::Relaxed)
    }

    /// Whether the run has posted that warning (`take_utf8_post`), or a load of the model of its
    /// loads posts it (`claim_utf8_post`), so that no read need name its member again.
    pub fn utf8_posted(&self) -> bool {
        self.utf8_posted.load(std::sync::atomic::Ordering::Relaxed)
            || self.utf8_claimed.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// `take_utf8_post` for a load of the model of the run's loads (`scope_utf8_loads`), which
    /// writes the warning with the load's text where the program writes it: true once, and none
    /// of the run's other posts is taken after (`take_utf8_post`). A search of the batches again
    /// from the start forgets it (`restart_utf8_loads`).
    pub fn claim_utf8_post(&self) -> bool {
        !self.utf8_posted.load(std::sync::atomic::Ordering::Relaxed)
            && !self
                .utf8_claimed
                .swap(true, std::sync::atomic::Ordering::Relaxed)
    }

    /// The batches are searched and their loads modelled again from the start, nothing of the
    /// first model written (BLASTP's search of the batches before a failed one): no load posts
    /// the warning yet (`claim_utf8_post`, `end_scope_utf8_loads`).
    pub fn restart_utf8_loads(&self) {
        self.utf8_claimed
            .store(false, std::sync::atomic::Ordering::Relaxed);
        self.scope_utf8_ended
            .store(false, std::sync::atomic::Ordering::Relaxed);
    }

    /// Notes that a reader met a reported sequence whose Seq-align id the lookup gives first to
    /// another sequence whose header holds a UTF8String where a VisibleString is expected
    /// (`algorithm::common::database::HeaderReader`'s scope read, `ParsedSubject::scope_utf8`).
    pub fn note_scope_utf8(&self) {
        self.scope_utf8
            .store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// The model of the run's loads reaches a batch whose own reads (its fetches, Seq-aligns and
    /// the reads of its reported sequences, `algorithm::common::database::Utf8Read`) post NCBI's
    /// warning of a UTF8String header, before any load of the batch: no load posts it after
    /// (`scope_utf8_loads`), whether the run has written it yet or not (a program that models
    /// every batch before it writes a report).
    pub fn end_scope_utf8_loads(&self) {
        self.scope_utf8_ended
            .store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// Whether a load whose Seq-align id's lookup reads such a header may post NCBI's warning: a
    /// reader met such a subject (`scope_utf8_seen`), the run has not posted it
    /// (`take_utf8_post`), no batch before the load posts it with its own reads
    /// (`end_scope_utf8_loads`), and NCBI makes the loads modelled now
    /// (`without_scope_utf8_loads`).
    pub fn scope_utf8_loads(&self) -> bool {
        self.scope_utf8_seen()
            && !self.utf8_posted()
            && !self
                .scope_utf8_ended
                .load(std::sync::atomic::Ordering::Relaxed)
            && !self
                .scope_utf8_held
                .load(std::sync::atomic::Ordering::Relaxed)
    }

    /// `model` (a model of loads, `algorithm::common::database::batch_prefetch`) where NCBI does
    /// not make the loads whose Seq-align id's lookup reads such a header: none of them posts the
    /// warning (`scope_utf8_loads`), which a later load or read may still post.
    #[cold]
    #[inline(never)]
    pub fn without_scope_utf8_loads<T>(&self, model: impl FnOnce() -> T) -> T {
        let held = self
            .scope_utf8_held
            .swap(true, std::sync::atomic::Ordering::Relaxed);
        let modelled = model();
        self.scope_utf8_held
            .store(held, std::sync::atomic::Ordering::Relaxed);
        modelled
    }

    /// The sequence at `at` (by its ordinal id).
    pub(crate) fn sequence_at(&self, at: SeqLocation) -> DbSequence {
        self.sequence(self.ordinal_of(self.location_oid(&at)))
    }

    /// Whether a reader noted such a sequence (`note_scope_utf8`): only then do the batches'
    /// prefetches look for it among their subjects (`prefetch_messages`).
    pub fn scope_utf8_seen(&self) -> bool {
        self.scope_utf8.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// The member of NCBI's warning of a read of the header at `at` that took a UTF8String for a
    /// VisibleString (`Utf8Reads`), where the run posts it now (`take_utf8_post`, once per
    /// process): `None` once posted, or where the record holds no such string before a failure.
    #[cold]
    #[inline(never)]
    pub fn take_utf8_post_at(&self, at: SeqLocation) -> Option<String> {
        if self.utf8_posted() {
            return None;
        }
        let path = self.utf8_string_path_at(at)?;
        self.take_utf8_post().then_some(path)
    }

    /// [`SeqDb::utf8_string_path`] of the sequence at `at`.
    #[cold]
    pub fn utf8_string_path_at(&self, at: SeqLocation) -> Option<String> {
        let record = self.volumes[at.volume].header_bytes(at.vol_oid).ok()?;
        crate::seqid::ber::utf8_string_path(&record)
    }

    /// The place of the sequence `ordinal` (`None` outside the database).
    pub fn location_of(&self, ordinal: u32) -> Option<SeqLocation> {
        self.find_volume(ordinal).ok()
    }

    /// [`SeqDb::defline_set`] of the sequence at `at`.
    pub fn defline_set_at(&self, at: SeqLocation) -> Result<DeflineSet, SeqDbError> {
        let set = self.volumes[at.volume].defline_set(at.vol_oid)?;
        self.ordinal_set(at.volume, self.location_oid(&at), set)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1995-1999
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:907-922
    /// Whether NCBI's header read of the sequence at `at` filters the set it read
    /// (`x_GetFilteredHeader`): the database is restricted by a list or the alias files attached
    /// an id list to the volume (`x_HaveIdFilter`, a list attached even where it kept nothing),
    /// the volume has a non-zero `MEMB_BIT`, or the volume has an `OID_MASK_TYPE`. The filter
    /// needs a set: a stored header of no bytes, no set, makes NCBI's read throw there
    /// (`Source::null_set` of blastdbcmd).
    pub fn filters_header_sets(&self, at: SeqLocation) -> bool {
        let volume = &self.volumes[at.volume];
        self.is_restricted()
            || volume.has_alias_lists()
            || volume.has_memb_bit()
            || volume.has_oid_mask_type()
    }

    /// [`SeqDb::defline_set_at`] as NCBI's `GetHdr` gives it, written over `set`, whose buffers
    /// are kept from one sequence to the next (`seqid::ber::decode_full_set_in`: blastdbcmd
    /// reads the header of one sequence after another): false where the stored header has no
    /// bytes, which `defline_set_at` gives as a set of no lines (`Volume::header_set_in`).
    pub fn header_set_in(&self, at: SeqLocation, set: &mut DeflineSet) -> Result<bool, SeqDbError> {
        let volume = at.volume;
        if self.volumes[volume]
            .header_set_in(at.vol_oid, set, None)?
            .is_none()
        {
            return Ok(false);
        }
        self.adjust_ordinals(volume, self.location_oid(&at), set)?;
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2319-2346
    /// The definition lines of the sequence at `at` (`defline_set_at`) and the header that
    /// NCBI's `GetBioseq` keeps in the Bioseq's `ASN1_BlastDefLine` user object
    /// (`x_GetFilteredBinaryHeader`): the stored record, or the set serialized again
    /// (`writedb::header::encode_defline_set`, NCBI's `CObjectOStreamAsnBinary`) where NCBI's read
    /// changed it: a `general` id `BL_ORD_ID` of a volume after the first, whose tag the read made
    /// the OID in the database (`x_GetHdrAsn1`), or a line that the definition-line rule dropped
    /// (`x_GetFilteredHeader`). `None` where the stored header has no bytes: NCBI's `GetHdr` gives
    /// no set there (`Volume::header_set_in`). The set is written over `set`, as
    /// [`SeqDb::header_set_in`] writes it.
    pub fn defline_set_with_header_in(
        &self,
        at: SeqLocation,
        set: &mut DeflineSet,
    ) -> Result<Option<Vec<u8>>, SeqDbError> {
        let volume = at.volume;
        let mut record = Vec::new();
        let Some(dropped) =
            self.volumes[volume].header_set_in(at.vol_oid, set, Some(&mut record))?
        else {
            return Ok(None);
        };
        let adjusted = self.adjust_ordinals(volume, self.location_oid(&at), set)?;
        if dropped || adjusted {
            record = crate::writedb::header::encode_defline_set(set, &mut Vec::new());
        }
        Ok(Some(record))
    }

    /// `set`, the definition lines of the sequence `oid` (its OID, which the messages name) in
    /// `volume`, with the tag of a `general` id `BL_ORD_ID` made its OID in the database
    /// (`defline_set`).
    fn ordinal_set(
        &self,
        volume: usize,
        oid: u32,
        mut set: DeflineSet,
    ) -> Result<DeflineSet, SeqDbError> {
        self.adjust_ordinals(volume, oid, &mut set)?;
        Ok(set)
    }

    /// Makes the tag of each `general` id `BL_ORD_ID` of `set`, the definition lines of the
    /// sequence `oid` in `volume`, its OID in the database (`adjust_ordinal_tags`, with the
    /// reader's error naming `oid`); whether a tag changed.
    fn adjust_ordinals(
        &self,
        volume: usize,
        oid: u32,
        set: &mut DeflineSet,
    ) -> Result<bool, SeqDbError> {
        adjust_ordinal_tags(set, self.volume_starts[volume])
            .map_err(|error| error.error(self.volumes[volume].base(), oid))
    }

    /// The Blast-def-line `taxid` of each definition line of `ordinal`, 0 for none: read
    /// apart from `deflines`, for the XML2 and JSON formats, the prefetch's taxonomy lookups and
    /// an empty stored title made from its source (`HeaderReader::stored_source`).
    pub fn defline_taxids(&self, ordinal: u32) -> Result<Vec<i32>, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].defline_taxids(at.vol_oid)
    }

    /// The definition lines of a volume's sequence as `deflines` gives them for the sequence
    /// `oid` (its OID, which the messages name).
    fn ordinal_deflines(
        &self,
        volume: usize,
        oid: u32,
        deflines: Result<Vec<Defline>, SeqDbError>,
    ) -> Result<Vec<OrdinalDefline>, SeqDbError> {
        let volume_start = self.volume_starts[volume];
        let vol = &self.volumes[volume];
        deflines?
            .into_iter()
            .map(|defline| {
                let bl_ord_id = u64::from(volume_start) + u64::from(defline.volume_ordinal);
                if bl_ord_id > i32::MAX as u64 {
                    return Err(SeqDbError::Format {
                        path: vol.base().to_path_buf(),
                        message: format!("BL_ORD_ID of sequence {oid} passes a signed 32-bit id"),
                    });
                }
                Ok(OrdinalDefline {
                    title: defline.title,
                    bl_ord_id: bl_ord_id as u32,
                    replaced_title_bytes: defline.replaced_title_bytes,
                    taxid: defline.taxid,
                    utf8: defline.utf8,
                })
            })
            .collect()
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): the ordinal ids of the database cut
    /// into ranges of at most `max_seqs` sequences (at least one) inside one
    /// volume each, in order: the ranges that a check of every sequence's
    /// header reads one at a time ([`SeqDb::visit_deflines`]).
    ///
    /// With a taxid list (`restrict_taxids`) a range is still `max_seqs` OIDs of a volume:
    /// the visits skip the OIDs that the list excludes, `first_ordinal` is the ordinal id of
    /// the first included OID from `first` on, and a range without one is left out.
    pub fn header_ranges(&self, max_seqs: u32) -> Vec<BlockRange> {
        let max_seqs = max_seqs.max(1);
        let mut ranges = Vec::new();
        for (volume, vol) in self.volumes.iter().enumerate() {
            let mut first = 0;
            while first < vol.num_oids() {
                let end = first + max_seqs.min(vol.num_oids() - first);
                let range = self.oid_range(volume, first, end);
                if !self.excludes_all(range) {
                    ranges.push(range);
                }
                first = end;
            }
        }
        ranges
    }

    /// The range of the OIDs `first..end` of `volume` (`header_ranges`).
    fn oid_range(&self, volume: usize, first: u32, end: u32) -> BlockRange {
        let oid = self.volume_starts[volume] + first;
        BlockRange {
            volume,
            first,
            end,
            first_ordinal: match &self.filter {
                Some(filter) => filter.rank(oid),
                None => oid,
            },
        }
    }

    /// Whether a taxid list excludes every OID of `range` (a range of `header_ranges`).
    fn excludes_all(&self, range: BlockRange) -> bool {
        let Some(filter) = &self.filter else {
            return false;
        };
        let Some(&start) = self.volume_starts.get(range.volume) else {
            return false;
        };
        filter.rank(start + range.end) == range.first_ordinal
    }

    /// The ordinal id of the sequence `oid` (volume-relative) of `range` (a range of
    /// `header_ranges`); `None` for an OID that a taxid list excludes.
    fn range_ordinal(&self, range: BlockRange, oid: u32) -> Option<u32> {
        match &self.filter {
            None => Some(range.first_ordinal + (oid - range.first)),
            Some(filter) => {
                let oid = self.volume_starts[range.volume] + oid;
                filter.contains(oid).then(|| filter.rank(oid))
            }
        }
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP): the number of ranges of `header_ranges(max_seqs)`.
    pub fn x_header_range_count(&self, max_seqs: u32) -> usize {
        let max_seqs = max_seqs.max(1);
        self.volumes
            .iter()
            .map(|vol| vol.num_oids().div_ceil(max_seqs) as usize)
            .sum()
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP): `header_ranges(max_seqs)[index]`, without the list (a list
    /// of UniRef50's 65,000 ranges is 1.5 MiB, whose free raises glibc's mmap threshold).
    pub fn x_header_range(&self, max_seqs: u32, index: usize) -> BlockRange {
        let max_seqs = max_seqs.max(1);
        let mut index = index;
        for (volume, vol) in self.volumes.iter().enumerate() {
            let count = vol.num_oids().div_ceil(max_seqs) as usize;
            if index < count {
                let first = index as u32 * max_seqs;
                let end = first + max_seqs.min(vol.num_oids() - first);
                return self.oid_range(volume, first, end);
            }
            index -= count;
        }
        panic!("header range {index} past the database's ranges");
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): calls `visit(ordinal, self.deflines(ordinal))` for
    /// the sequences of `range` in order until `visit` returns false; returns whether every
    /// call returned true. The headers of the range are read with one read of the index file
    /// and one read of the header file per MiB (`Volume::visit_deflines`) instead of two reads
    /// per sequence.
    /// The titles and the taxids of a range come in one pass over its headers (a record whose
    /// taxid does not decode is an `Err`, as in `defline_taxids`).
    pub fn visit_deflines(
        &self,
        range: BlockRange,
        mut visit: impl FnMut(u32, Result<Vec<OrdinalDefline>, SeqDbError>) -> bool,
    ) -> bool {
        let Some(vol) = self.volumes.get(range.volume) else {
            return if range.is_empty() {
                true
            } else {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            };
        };
        if self.excludes_all(range) {
            return true;
        }
        vol.visit_deflines_as(
            range.first,
            range.end,
            header_piece_bytes(),
            true,
            &mut |oid, deflines| {
                let Some(ordinal) = self.range_ordinal(range, oid) else {
                    return true;
                };
                let deflines = self.ordinal_deflines(
                    range.volume,
                    self.volume_starts[range.volume] + oid,
                    deflines,
                );
                visit(ordinal, deflines)
            },
        )
    }

    /// `visit(ordinal, self.defline_taxids(ordinal))` for the sequences of `range` in order until
    /// `visit` returns false; returns whether every call returned true. The headers of the range
    /// are read as [`SeqDb::visit_deflines`] reads them (one read of the index file and one read of
    /// the header file per MiB), decoded as `defline_taxids` decodes them: the taxdb pre-pass of
    /// BLASTX 15/16 (`algorithm::common::database::taxdb_may_refuse`) walks the database with it.
    pub fn visit_defline_taxids(
        &self,
        range: BlockRange,
        mut visit: impl FnMut(u32, Result<Vec<i32>, SeqDbError>) -> bool,
    ) -> bool {
        let Some(vol) = self.volumes.get(range.volume) else {
            return if range.is_empty() {
                true
            } else {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            };
        };
        if self.excludes_all(range) {
            return true;
        }
        vol.visit_defline_taxids(
            range.first,
            range.end,
            header_piece_bytes(),
            &mut |oid, taxids| match self.range_ordinal(range, oid) {
                Some(ordinal) => visit(ordinal, taxids),
                None => true,
            },
        )
    }

    /// [`SeqDb::visit_defline_taxids`] for a database of any Seq-ids: the taxid of every
    /// definition line of each sequence of `range` (0 for none), decoded as [`SeqDb::defline_set`]
    /// decodes the definition lines (parsed Seq-ids too), read as `visit_deflines` reads them.
    pub fn visit_defline_set_taxids(
        &self,
        range: BlockRange,
        mut visit: impl FnMut(u32, Result<Vec<i32>, SeqDbError>) -> bool,
    ) -> bool {
        let Some(vol) = self.volumes.get(range.volume) else {
            return if range.is_empty() {
                true
            } else {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            };
        };
        if self.excludes_all(range) {
            return true;
        }
        vol.visit_defline_set_taxids(
            range.first,
            range.end,
            header_piece_bytes(),
            &mut |oid, taxids| match self.range_ordinal(range, oid) {
                Some(ordinal) => visit(ordinal, taxids),
                None => true,
            },
        )
    }

    /// `visit(ordinal, self.read_header(ordinal))` for the sequences of `range` in order until
    /// `visit` returns false; returns whether every call returned true. The headers of the range
    /// are read as [`SeqDb::visit_deflines`] reads them (one read of the index file and one read
    /// of the header file per MiB).
    pub fn visit_sequence_headers(
        &self,
        range: BlockRange,
        mut visit: impl FnMut(u32, Result<SequenceHeader, SeqDbError>) -> bool,
    ) -> bool {
        let Some(vol) = self.volumes.get(range.volume) else {
            return if range.is_empty() {
                true
            } else {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            };
        };
        if self.excludes_all(range) {
            return true;
        }
        vol.visit_sequence_headers(
            range.first,
            range.end,
            header_piece_bytes(),
            &mut |oid, header| match self.range_ordinal(range, oid) {
                Some(ordinal) => {
                    let oid = self.volume_starts[range.volume] + oid;
                    visit(ordinal, self.sequence_header(range.volume, oid, header))
                }
                None => true,
            },
        )
    }
}

fn volume_range_error(range: BlockRange, db: &SeqDb) -> SeqDbError {
    SeqDbError::Range(format!(
        "volume {} is outside the database ({} volumes)",
        range.volume,
        db.volumes.len()
    ))
}

/// Number of leading values of the ascending `range` for which `fits` holds.
fn partition_point_u32(range: std::ops::Range<u32>, fits: impl Fn(u32) -> bool) -> usize {
    let (mut lo, mut hi) = (range.start, range.end);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if fits(mid) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    (lo - range.start) as usize
}

/// The header of one sequence (`SeqDb::read_header`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SequenceHeader {
    /// A record whose only Seq-ids are `gnl|BL_ORD_ID|n` of the sequence's number n in its volume
    /// (`SeqDb::deflines`).
    Ordinal(Vec<OrdinalDefline>),
    /// A record with parsed Seq-ids (`SeqDb::defline_set`).
    Parsed(DeflineSet),
}

/// Whether a definition line of `set` holds a `general` id `BL_ORD_ID`, whose tag NCBI's
/// header reads of a volume after the first change (`adjust_ordinal_tags`).
pub(super) fn has_ordinal_tag(set: &DeflineSet) -> bool {
    set.deflines
        .iter()
        .flat_map(|defline| &defline.ids)
        .any(|id| matches!(id, ParsedSeqId::General(dbtag) if dbtag.db == b"BL_ORD_ID"))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2222-2236
// The overload with a stream (seqdbvol.cpp:2275-2289, for `GetSeqIDs`) does the same. The filtered
// header read (`x_GetFilteredHeader`, seqdbvol.cpp:1992-1993 for `GetBioseq`, 2108-2109 for
// `GetSeqIDs`) asks for this change before a list's rule reads the Seq-ids (`Volume::line_mask`),
// as `GetHdr` does for the reports.
/// `set`, the definition lines of a sequence of a volume whose first OID is `volume_start`, with
/// the tag of each `general` id `BL_ORD_ID` made the stored tag plus `volume_start` (`x_GetHdrAsn1`
/// with `adjust_oids`): the OID in the database where makeblastdb stored the number in the volume.
/// Fails where NCBI's `GetId` throws (a string tag) and where the sum passes a signed 32-bit id
/// (`OrdinalTagError`), with `set` then partly changed. Ok: whether a tag was changed (NCBI's
/// `changed`, which makes `x_GetFilteredBinaryHeader` serialize the set again).
pub(super) fn adjust_ordinal_tags(
    set: &mut DeflineSet,
    volume_start: u32,
) -> Result<bool, OrdinalTagError> {
    let mut changed = false;
    if volume_start == 0 {
        return Ok(changed);
    }
    for id in set
        .deflines
        .iter_mut()
        .flat_map(|defline| defline.ids.iter_mut())
    {
        let ParsedSeqId::General(dbtag) = id else {
            continue;
        };
        if dbtag.db != b"BL_ORD_ID" {
            continue;
        }
        // `GetId` of a string tag throws in NCBI (an invalid CHOICE selection). makeblastdb
        // -parse_seqids stores a FASTA tag with a leading 0 (`gnl|BL_ORD_ID|0`) as a string.
        let ObjectId::Id(tag) = dbtag.tag else {
            return Err(OrdinalTagError::StringTag);
        };
        let Ok(adjusted) = i32::try_from(i64::from(volume_start) + i64::from(tag)) else {
            return Err(OrdinalTagError::PastSigned32);
        };
        dbtag.tag = ObjectId::Id(adjusted);
        changed = true;
    }
    Ok(changed)
}

/// A `BL_ORD_ID` tag that NCBI's header read of a volume after the first cannot make an OID in
/// the database (`adjust_ordinal_tags`). NCBI makes that read only for a sequence it reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum OrdinalTagError {
    /// A string tag: NCBI's `GetId` throws.
    StringTag,
    /// The stored tag plus the volume's first OID passes `i32::MAX`.
    PastSigned32,
}

impl OrdinalTagError {
    /// The reader's error for the sequence `oid` (its OID in the database, which the message
    /// names) of the volume `path`.
    pub(super) fn error(self, path: &Path, oid: u32) -> SeqDbError {
        let message = match self {
            OrdinalTagError::StringTag => format!("BL_ORD_ID of sequence {oid} has a string tag"),
            OrdinalTagError::PastSigned32 => {
                format!("BL_ORD_ID of sequence {oid} passes a signed 32-bit id")
            }
        };
        SeqDbError::Format {
            path: path.to_path_buf(),
            message,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:304-309 (cited at
// `algorithm::common::database::Utf8Read`)
/// The first, in NCBI's order, of a reader's reads of header records that hold a UTF8String
/// where a VisibleString is expected (`DeflineSet::utf8`, `OrdinalDefline::utf8`). NCBI's reader
/// posts its warning at the first such read of the process (`ERR_POST_X_ONCE`,
/// [`SeqDb::take_utf8_post`]), naming the member the record holds
/// ([`SeqDb::utf8_string_path`]), wherever the read is: this is the one note of every LOSAT read
/// that stands for one of NCBI's reads of a header. Its readers: a search's reads of its
/// reported sequences and of the other sequences its scope reads
/// (`algorithm::common::database::HeaderReader`), the fetches with `check_oid_exclusion` and the
/// Seq-aligns under a seqid list (`LineCheck`), the lookups of an accession whose version NCBI
/// checks against the candidates' Seq-ids (`IdLookup`, `x_CheckVersions`), and blastdbcmd's
/// reads (`Source`). A reader gives each read its place in NCBI's order (`order`, the least
/// first; of equal places the one noted first), so that reads on several threads note the same
/// first; a record is noted only where its read found the string (cold).
#[derive(Debug, Default)]
pub struct Utf8Reads {
    first: std::sync::Mutex<Option<Utf8Note>>,
    /// Whether `first` may hold a read (a reader asks this before it locks).
    noted: AtomicBool,
}

/// A read that `Utf8Reads` noted: its place in the reader's order and the sequence read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Utf8Note {
    pub order: u128,
    pub at: SeqLocation,
}

impl Utf8Reads {
    /// Notes the read of the sequence at `at` at `order`, where it comes before the reads noted
    /// so far.
    #[cold]
    #[inline(never)]
    pub fn note(&self, order: u128, at: SeqLocation) {
        let mut first = self
            .first
            .lock()
            .expect("the UTF8String note is not poisoned");
        if first.is_none_or(|first| order < first.order) {
            *first = Some(Utf8Note { order, at });
        }
        self.noted.store(true, std::sync::atomic::Ordering::Release);
    }

    /// [`Utf8Reads::note`] of the sequence `ordinal` of `db` (none outside the database).
    #[cold]
    #[inline(never)]
    pub fn note_ordinal(&self, db: &SeqDb, order: u128, ordinal: u32) {
        if let Ok(at) = db.find_volume(ordinal) {
            self.note(order, at);
        }
    }

    /// The first read noted.
    pub fn first(&self) -> Option<Utf8Note> {
        if !self.noted.load(std::sync::atomic::Ordering::Acquire) {
            return None;
        }
        *self
            .first
            .lock()
            .expect("the UTF8String note is not poisoned")
    }

    /// `first`, which the reader then forgets (each batch of a search takes its own).
    pub fn take(&self) -> Option<Utf8Note> {
        if !self.noted.load(std::sync::atomic::Ordering::Acquire) {
            return None;
        }
        self.take_noted()
    }

    #[cold]
    #[inline(never)]
    fn take_noted(&self) -> Option<Utf8Note> {
        let mut first = self
            .first
            .lock()
            .expect("the UTF8String note is not poisoned");
        self.noted
            .store(false, std::sync::atomic::Ordering::Release);
        first.take()
    }
}

/// A definition line with its database-wide `BL_ORD_ID`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrdinalDefline {
    /// The title: the whole FASTA definition line after `>`.
    pub title: Vec<u8>,
    /// `n` of `gnl|BL_ORD_ID|n`.
    pub bl_ord_id: u32,
    /// Number of stored title bytes outside ' ' ..= '~' that were read as '#'
    /// (`Defline::replaced_title_bytes`).
    pub replaced_title_bytes: usize,
    /// `Defline::taxid` (0 for none), read with the title by [`SeqDb::deflines`],
    /// [`SeqDb::visit_deflines`] and [`SeqDb::read_header`].
    pub taxid: i32,
    /// `Defline::utf8`: the record holds a UTF8String where a VisibleString is expected
    /// ([`SeqDb::utf8_string_path`] names the member of NCBI's warning).
    pub utf8: bool,
}
