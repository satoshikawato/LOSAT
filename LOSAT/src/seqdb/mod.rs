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
pub mod header;
pub mod index;
pub mod isam;
pub mod isam_num;
pub mod lmdb_read;
pub mod lookup;
pub mod mask;
pub(crate) mod rps;
pub mod seqidlist;
pub mod tax_files;
pub mod tax_sqlite;
pub mod taxdb;
pub mod volume;

pub use filter::{DbList, ListWarning, SeqidList, TaxidList};
pub(crate) use filter::{Oid, OidMap};

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

use crate::seqid::model::ObjectId;
use crate::seqid::{DeflineSet, ParsedSeqId};

pub use alias::SearchPath;
pub use header::Defline;
pub use mask::{BlockMasks, MaskAlgorithm, MaskAlgorithms, MaskBlock, MaskRanges, NameMaskError};
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
    /// The message is the text of NCBI's exception.
    NotFound(String),
    /// A file that exists could not be opened.
    Open {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The operating system failed to read a file.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// A file is not a valid database file.
    Format { path: PathBuf, message: String },
    /// A file fails one of the integrity checks of NCBI's reader (`SEQDB_FILE_ASSERT`, a
    /// `CSeqDBException` whose text names a line of NCBI's source): LOSAT's words, and NCBI's
    /// text where NCBI's search replaces it ([`SeqDbError::ncbi_search_message`]).
    Integrity { path: PathBuf, message: String },
    /// A valid database that needs something LOSAT does not read.
    Unsupported { path: PathBuf, message: String },
    /// A request outside the database (ordinal id, block range).
    Range(String),
    /// An exception of NCBI's C++ Toolkit that is not a `CSeqDBException`: the
    /// `CStringException` of an alias file value that does not convert (`alias::walk_totals`,
    /// `alias::compute_masks`). NCBI's programs report it as `Error: <what()>` with exit
    /// status 255 ([`SeqDbError::ncbi_exception`]).
    Ncbi(crate::ncbi_exception::NcbiException),
}

impl std::fmt::Display for SeqDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqDbError::NotFound(message) | SeqDbError::Range(message) => f.write_str(message),
            SeqDbError::Open { path, source } | SeqDbError::Io { path, source } => {
                write!(f, "{}: {}", path.display(), source)
            }
            SeqDbError::Format { path, message } | SeqDbError::Integrity { path, message } => {
                write!(f, "{}: {}", path.display(), message)
            }
            SeqDbError::Unsupported { path, message } => {
                write!(f, "{}: unsupported database: {}", path.display(), message)
            }
            SeqDbError::Ncbi(exception) => f.write_str(exception.what().trim_end()),
        }
    }
}

impl SeqDbError {
    /// Whether this is the error of the ordinal-id decoders for a record with a Seq-id other
    /// than the `gnl|BL_ORD_ID|n` of the sequence's own number in its volume (a database built
    /// with `-parse_seqids`; `SeqDb::read_header`).
    pub fn is_parsed_seqids(&self) -> bool {
        matches!(self, SeqDbError::Unsupported { message, .. } if message == header::PARSED_SEQIDS)
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    // `BLAST_UNKNOWN_ERROR` is 255 (blast_app_util.hpp:161); blastdbcmd catches the exception the
    // same way (`CATCH_ALL`). NCBI's guess of a database type moves on to the nucleotide
    // database only after a `CSeqDBException` (seqdb.cpp:121-123), so not after this one.
    /// The NCBI exception of this error when it is not a database error of NCBI's: a search
    /// program or blastdbcmd writes `Error: <what()>` and ends with exit status 255.
    pub fn ncbi_exception(&self) -> Option<&crate::ncbi_exception::NcbiException> {
        match self {
            SeqDbError::Ncbi(exception) => Some(exception),
            _ => None,
        }
    }
}

impl std::error::Error for SeqDbError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            SeqDbError::Open { source, .. } | SeqDbError::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// The texts of `SeqDbError::Format` that are NCBI's own: the errors NCBI
/// raises for the same file with the same words. Every other format error is
/// LOSAT's (a damaged file that NCBI reads past the end of its mapping).
const NCBI_FORMAT_MESSAGES: [&str; 13] = [
    // seqdbfile.cpp:181-185
    "Error: Not a valid version 4 or 5 database.",
    // seqdbfile.cpp:222-226
    "Error: requested sequence type does not match DB.",
    // seqdblmdbset.cpp:224, 231
    "DB list contains both Version 4 and Version 5 dbs",
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
];

/// seqdbvol.cpp:3365-3382 (the title or the counts follow), seqdbimpl.cpp:2456-2459 (the volume
/// and the value follow).
const NCBI_COLUMN_MESSAGE_PREFIXES: [&str; 3] = [
    "Error: duplicate column title [",
    "Error: column has wrong #oids [",
    "Error: volume (",
];

/// seqdbalias.cpp:740-754; the name of the alias file follows.
const NCBI_EMPTY_DBLIST_PREFIX: &str = "No database names were found in alias file [";

impl SeqDbError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
    // ```c++
    //     catch (const CSeqDBException& e) {                                      \
    // 		if (e.GetErrCode() == CSeqDBException::eTooManyOpenFiles) {         \
    // 			string err_msg =                                                \
    //                 "Too many open files, please raise the open file limit";    \
    //         	LOG_POST(Error << "BLAST Database error: " << err_msg);         \
    //             exit_code = BLAST_DATABASE_ERROR;                               \
    // 		}                                                                   \
    // 		else if (e.GetErrCode() == CSeqDBException::eMemoryMappingFailure){ \
    // 			string err_msg = "Database memory map file error";              \
    //             if (!e.GetMsg().empty())                                        \
    //                 err_msg += " - " + e.GetMsg();                              \
    //         	LOG_POST(Error << "BLAST Database error: " << err_msg);         \
    //             exit_code = BLAST_DATABASE_ERROR;                               \
    // 		}                                                                   \
    //         else {                                                              \
    //         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
    //         	exit_code = BLAST_DATABASE_ERROR;                               \
    //         }                                                                   \
    //     }                                                                       \
    // ```
    // `BLAST_DATABASE_ERROR` is 2 (blast_app_util.hpp:153). The exception of
    // a file that cannot be opened is quoted at `file::DbFile::open`: its
    // message is `CMemoryFile: Cannot memory map file '<path>':<strerror>`,
    // and it has the code `eTooManyOpenFiles` when that text holds "Too many
    // open files" (the `strerror` text of EMFILE and of ENFILE).
    //
    // `std::io::Error` prints an operating system error as the `strerror`
    // text followed by ` (os error N)`; the text before that suffix is what
    // NCBI appends. On Windows NCBI's text comes from another call
    // (`WIN_LAST_ERROR_STR`), which LOSAT does not reproduce: the error is
    // then reported as LOSAT's own.
    /// The text NCBI BLAST+ gives for this error, when it is one of NCBI's
    /// own: a search program prints `BLAST Database error: <text>` and ends
    /// with exit status 2. `None`: the error is LOSAT's (a database that
    /// LOSAT refuses, a damaged file, a failed read), to be reported in
    /// LOSAT's words.
    pub fn ncbi_message(&self) -> Option<String> {
        match self {
            SeqDbError::NotFound(message) => Some(message.clone()),
            SeqDbError::Format { message, .. } => (NCBI_FORMAT_MESSAGES
                .contains(&message.as_str())
                || message.starts_with(NCBI_EMPTY_DBLIST_PREFIX)
                || NCBI_COLUMN_MESSAGE_PREFIXES
                    .iter()
                    .any(|prefix| message.starts_with(prefix)))
            .then(|| message.clone()),
            SeqDbError::Open { path, source } => {
                let code = source.raw_os_error()?;
                let text = source.to_string();
                let strerror = text.strip_suffix(&format!(" (os error {code})"))?;
                if strerror.contains("Too many open files") {
                    Some("Too many open files, please raise the open file limit".to_string())
                } else if cfg!(unix) {
                    Some(format!(
                        "Database memory map file error - CMemoryFile: Cannot memory map file \
                         '{}':{strerror}",
                        path.display()
                    ))
                } else {
                    None
                }
            }
            SeqDbError::Io { .. }
            | SeqDbError::Integrity { .. }
            | SeqDbError::Unsupported { .. }
            | SeqDbError::Range(_)
            | SeqDbError::Ncbi(_) => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:229-252
    // ```c++
    // 	try {
    // 	    m_PrelimSearch->SetNumberOfThreads(GetNumberOfThreads());
    // 	    m_InternalData = m_PrelimSearch->Run();
    // 	} catch( CIndexedDbException & e) {
    //         _TRACE(e.GetMsg());
    // 	    throw;
    // 	} catch (CBlastException & e) {
    //         _TRACE(e.GetMsg());
    // 		if(e.GetErrCode() == CBlastException::eCoreBlastError) {
    // 			throw;
    // 		}
    // 	}
    // 	catch (CSeqDBException & e) {
    // 		if (e.GetErrCode() == CSeqDBException::eTooManyOpenFiles) {
    //             _TRACE(e.GetMsg());
    // 			string err_msg = BlastErrorCode2String(BLASTERR_DB_TOO_MANY_OPEN_FILES);
    // 			NCBI_THROW(CBlastException, eCoreBlastError, err_msg);
    // 		}
    // 		else {
    //             _TRACE(e.GetMsg());
    // 			string err_msg = BlastErrorCode2String(BLASTERR_DB_MEMORY_MAP);
    // 			NCBI_THROW(CBlastException, eCoreBlastError, err_msg);
    // 		}
    // 	}
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:227-236
    // ```c++
    //     case BLASTERR_DB_MEMORY_MAP :
    //     	new_msg->message = strdup("Database memory map file error");
    //         new_msg->severity = eBlastSevFatal;
    //         new_msg->context = context;
    //         break;
    //     case BLASTERR_DB_TOO_MANY_OPEN_FILES :
    //     	new_msg->message = strdup("Too many open files, please raise the open file limit");
    //         new_msg->severity = eBlastSevFatal;
    //         new_msg->context = context;
    //         break;
    // ```
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
            message if message == TOO_MANY => Some(TOO_MANY),
            _ => Some("Database memory map file error"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
// ```c++
//     if (args.Exist(kArgDb) && args[kArgDb]) {
// 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
//
//         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
//                                              mol_type));
// ```
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
    /// The id the caller gave (an ordinal id or an OID), for messages.
    id: u32,
}

/// A database resolved without opening its volumes
/// ([`SeqDb::locate_one_name_metadata`]).
#[derive(Debug)]
pub(crate) struct LocatedDatabase {
    /// `GetTitle()`.
    pub title: Vec<u8>,
    /// The volume paths without extension, in NCBI's order.
    pub volumes: Vec<String>,
    /// The alias files (paths with extension), sorted.
    pub alias_files: Vec<String>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:377-391
// ```c++
// static bool
// s_FillDbInfoLocally(const string& dbname,
//                     CAlignFormatUtil::SDbInfo& info,
//                     int dbfilt_algorithm)
// {
//     CRef<CSeqDB> seqdb(new CSeqDB(dbname, info.is_protein
//                           ? CSeqDB::eProtein : CSeqDB::eNucleotide));
//     if ( !seqdb ) {
//         return false;
//     }
//     info.name = seqdb->GetDBNameList();
//     info.definition = seqdb->GetTitle();
//     if (info.definition.empty())
//         info.definition = info.name;
//     info.date = seqdb->GetDate();
//     info.total_length = seqdb->GetTotalLength();
//     info.number_seqs = seqdb->GetNumSeqs();
// ```
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
    /// `GetNumSeqs`.
    pub num_seqs: u32,
    /// `GetTotalLength`.
    pub total_length: u64,
    /// The filtering algorithms of the name's database, which the formatter of a search with a
    /// database mask reads (`GetMaskAlgorithmDetails`, `SeqDb::name_mask_details`), or the error
    /// of reading them, met only where the formatter asks for them; read with the record when
    /// the search has a database mask (`read_name_records`), `None` otherwise.
    pub mask_algorithms: Option<Result<MaskAlgorithms, SeqDbError>>,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
// ```c
// #define BLAST_SEQSRC_MINLENGTH  10    /**< Default minimal sequence length */
// ```
/// What NCBI's `CSeqDB` answers for the shortest length when the alias walk gives none
/// (`SeqDb::min_length`).
pub(crate) const BLAST_SEQSRC_MINLENGTH: i32 = 10;

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
    first_oid_length: Option<u64>,
    /// The walks of the alias tree, with the exception of each one whose value does not
    /// convert (`check_alias_values`).
    alias_totals: alias::AliasTotals,
    /// Whether the alias files give an OID range (`alias::range_mask`).
    alias_ranges: bool,
    /// NCBI's `m_OidListSetup`: whether its OID list was built, which a database opened
    /// without a list does at its first use (`SeqDb::check_or_find_oid`).
    oid_list_set_up: AtomicBool,
    /// Whether a list was applied (`SeqDb::restrict`): NCBI's `CSeqDB` then scans the totals.
    list_applied: bool,
    /// The name the database was opened with (`-db`), for LOSAT's refusals.
    name: String,
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
    /// `SeqDb::name_records`.
    name_records: Vec<NameRecord>,
    /// The LMDB files of a version 5 database for `SeqDb::seqid_lookup`, opened at its first use.
    lookup_files: OnceLock<lookup::LookupFiles>,
    /// The ISAMs of the volumes for `SeqDb::seqid_lookup`, each opened at its first use.
    lookup_isams: OnceLock<lookup::VolumeIsams>,
    /// The filtering algorithms of the mask columns, built on first use (`mask.rs`).
    mask_algorithms: OnceLock<MaskAlgorithms>,
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
/// in its order (seqdbimpl.cpp:117-139, quoted at [`SeqDb::open_every_oid_in`]).
const UNLISTED_WALK_KEYS: &[alias::WalkKey] = &[
    alias::WalkKey::Nseq,
    alias::WalkKey::Length,
    alias::WalkKey::MinSeqLength,
    alias::WalkKey::StatsNseq,
    alias::WalkKey::StatsTotlen,
];

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1093-1102
// ```c++
// string CSeqDBImpl::x_FixString(const string & s) const
// {
//     CHECK_MARKER();
//     for(int i = 0; i < (int) s.size(); i++) {
//         if (s[i] == char(0)) {
//             return string(s,0,i);
//         }
//     }
//     return s;
// }
// ```
fn fix_string(s: &[u8]) -> &[u8] {
    match s.iter().position(|&byte| byte == 0) {
        Some(at) => &s[..at],
        None => s,
    }
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
        Self::open_in(name, seqtype, &SearchPath::new(&cwd, None, None)?)
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
    /// ask for them (`Volume::open_as`).
    pub fn open_in(name: &str, seqtype: SeqType, path: &SearchPath) -> Result<Self, SeqDbError> {
        let db = Self::open_search_in(name, seqtype, path)?;
        db.check_alias_values()?;
        Ok(db)
    }

    /// [`SeqDb::open_in`] without the errors of the alias files' counts, which depend on a
    /// list that the caller applies next ([`SeqDb::restrict`]); the caller then calls
    /// [`SeqDb::check_alias_values`] before it reads the database. The errors of the OID
    /// ranges (`FIRST_OID`, `LAST_OID`) come here, as NCBI's `ComputeMasks` runs first.
    pub fn open_search_in(
        name: &str,
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<Self, SeqDbError> {
        let mut db = Self::open_tree(alias::resolve_in(name, seqtype, path)?, name, seqtype)?;
        db.search_path = Some(path.clone());
        Ok(db)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:150-153
    // ```c++
    //         CBlastFormatUtil::GetBlastDbInfo(m_DbInfo, m_DbName, m_DbIsAA,
    //                                    dbfilt_algorithm, is_remote_search);
    // ```
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
                // list, its walks of the counts meet their errors here (`check_alias_values`).
                if self.list_applied {
                    use alias::WalkKey;
                    self.walk_errors(&[WalkKey::Nseq, WalkKey::Length, WalkKey::MinSeqLength])?;
                    self.check_printed_counts()?;
                }
                vec![self.name_record(one, masks)]
            }
            _ => {
                let path = match &self.search_path {
                    Some(path) => path.clone(),
                    None => {
                        SearchPath::new(&std::env::current_dir().unwrap_or_default(), None, None)?
                    }
                };
                let mut records = Vec::with_capacity(names.len());
                for one in &names {
                    let mut alone = Self::open_in(one, self.seqtype, &path)?;
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
    /// (`SeqDB_ResolveDbPath("taxdb.bti")`, data4xml2format.cpp:x_InitCommon): NCBI's XML2 and
    /// JSON formats then write the scientific name of a taxid.
    pub fn taxdb_found(&self) -> bool {
        self.search_path
            .as_ref()
            .is_some_and(|path| path.resolve_file("taxdb.bti").is_some())
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
    /// ```c++
    /// bool CSeqDbSeqInfoSrc::HasGiList() const
    /// {
    ///     return !! m_iSeqDb->GetGiList();
    /// }
    /// ```
    pub fn has_positive_list(&self) -> bool {
        self.positive_list
    }

    /// The taxids of the search's list object after the volumes opened, ascending, which the
    /// BLAST archive writes (`DbFilterArgs::archive_parameter`): of a positive list those the
    /// database holds (`CSeqDBLMDBSet::TaxIdsToOids`, `tax_ids.swap(rv_tax_ids)`), of a negative
    /// list all (`CSearchDatabase::x_InitializeDb` opens the database with a copy of it). Empty
    /// without a taxid list.
    pub fn list_taxids(&self) -> &[i32] {
        &self.list_taxids
    }

    /// The search path the database was opened with (`open_in`), none for a database opened
    /// by another way.
    pub fn search_path(&self) -> Option<&SearchPath> {
        self.search_path.as_ref()
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:566-572
    // ```c++
    //     // Print volume names
    //     vector<string> volumes;
    //     m_BlastDb->FindVolumePaths(volumes,false);
    //     out << endl << "Volumes:" << endl;
    //     ITERATE(vector<string>, file_name, volumes) {
    //         out << "\t" << *file_name << endl;
    //     }
    // ```
    /// What `blastdbcmd -info` lists under "Volumes:" for the database that
    /// [`SeqDb::open_in`] opens with the same arguments: the volumes of the first two levels
    /// of its alias tree and the alias files of the third, sorted by base name
    /// (`alias::info_volumes`). The alias files are read again; no volume file is opened.
    pub fn info_volumes(
        name: &str,
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<Vec<String>, SeqDbError> {
        alias::info_volumes(name, seqtype, path)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:929-934
    // ```c++
    // void CSeqDBAliasNode::FindVolumePaths(vector<string> & vols, vector<string> * alias, bool recursive) const
    // {
    //
    //     set<string> volset;
    //     set<string> aliset;
    //
    //     if (recursive) {
    //         x_FindVolumePaths(volset, aliset);
    // ```
    /// `FindVolumePaths(paths, alias, true)` of the database that [`SeqDb::open_in`] opens
    /// with the same arguments: every volume path (without extension) and every alias file
    /// of its alias tree, each once, sorted. No volume file is opened.
    pub fn volume_and_alias_paths(
        name: &str,
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<(Vec<String>, Vec<String>), SeqDbError> {
        let tree = alias::resolve_in(name, seqtype, path)?;
        let volumes: std::collections::BTreeSet<String> = tree.volumes.into_iter().collect();
        let mut alias_files = std::collections::BTreeSet::new();
        alias::collect_alias_files(&tree.root, &mut alias_files);
        Ok((
            volumes.into_iter().collect(),
            alias_files.into_iter().collect(),
        ))
    }

    /// The order in which NCBI searches and numbers the volumes of a database
    /// (`SeqDB_CompareVolume`, the sort at the end of `FindVolumePaths`): by file name (the
    /// text after the last path separator), then by the whole path. For volumes in one
    /// directory it is the byte order of the paths.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:190-197
    /// ```c++
    /// bool SeqDB_CompareVolume(const string & s1, const string & s2)
    /// {
    ///     string x1, x2;
    ///     CSeqDB_Path(s1).FindBaseName().GetString(x1);
    ///     CSeqDB_Path(s2).FindBaseName().GetString(x2);
    ///     if (x1 != x2) return (x1 < x2);
    ///     else return (s1 < s2);
    /// }
    /// ```
    pub fn compare_volume_paths(a: &Path, b: &Path) -> std::cmp::Ordering {
        fn bytes(path: &Path) -> &[u8] {
            path.as_os_str().as_encoded_bytes()
        }
        fn base_name(path: &[u8]) -> &[u8] {
            let separator = |byte: &u8| *byte == b'/' || (cfg!(windows) && *byte == b'\\');
            match path.iter().rposition(separator) {
                Some(at) => &path[at + 1..],
                None => path,
            }
        }
        let (a, b) = (bytes(a), bytes(b));
        base_name(a).cmp(base_name(b)).then_with(|| a.cmp(b))
    }

    /// The number of the volume that holds `ordinal` (`CSeqDBVolSet::FindVol`).
    pub fn volume_index(&self, ordinal: u32) -> Result<usize, SeqDbError> {
        self.find_volume(ordinal).map(|at| at.volume)
    }

    /// Opens the one database at the path `name` (without extension), taken
    /// as it is: a space or a double quote in it is a character of the path
    /// (`alias::resolve_one_name`). Everything else is as in [`SeqDb::open`],
    /// including NCBI's rule that a volume named in an alias file is looked
    /// for in the current directory first.
    #[cfg(test)]
    pub(crate) fn open_one_name(name: &str, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let tree = alias::resolve_one_name(
            name,
            seqtype,
            alias::LocalCopy::CurrentDir,
            &SearchPath::current_dir()?,
        )?;
        Self::open_tree(tree, name, seqtype)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1613-1616
    // ```c++
    //     vector<string> paths;
    //     vector<string> alias;
    //     m_Impl->FindVolumePaths(paths, alias, true);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1726
    // ```c++
    // 	m->SetDescription(GetTitle());
    // ```
    /// Resolves the one database at the path `name` as `SeqDb::open_one_name` does, without
    /// opening its volumes (the alias files are read and the index file of every volume they name
    /// is located), and returns what `CSeqDB::GetDBMetaData` reads of its alias tree: the title
    /// (`GetTitle`: the alias files' `TITLE`, else the volumes' titles, which `volume_title` gives
    /// for a volume path), the volume paths and the paths of the alias files
    /// (`FindVolumePaths(paths, alias, true)`). No volume file is read, so the cost does not grow
    /// with the number of sequences. `search` is the search path that the message of a database
    /// which is not there names; only `name` itself is looked for.
    pub(crate) fn locate_one_name_metadata(
        name: &str,
        seqtype: SeqType,
        search: &SearchPath,
        volume_title: &dyn Fn(&str) -> Option<Vec<u8>>,
    ) -> Result<LocatedDatabase, SeqDbError> {
        let tree = alias::resolve_one_name(name, seqtype, alias::LocalCopy::CurrentDir, search)?;
        let mut title = Vec::new();
        alias::walk_title(&tree.root, volume_title, &mut title);
        let mut alias_files = std::collections::BTreeSet::new();
        alias::collect_alias_files(&tree.root, &mut alias_files);
        Ok(LocatedDatabase {
            title,
            volumes: tree.volumes,
            alias_files: alias_files.into_iter().collect(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.cpp:101-116
    // ```c++
    // void CSeqDBVolSet::x_AddVolume(CSeqDBAtlas        & atlas,
    //                                const string       & nm,
    //                                char                 pn,
    //                                CSeqDBGiList       * user_list,
    //                                CSeqDBNegativeList * neg_list,
    //                                CSeqDBLockHold     & locked)
    // {
    //     int num_oids = x_GetNumOIDs();
    //
    //     CSeqDBVol * new_volp =
    //         new CSeqDBVol(atlas, nm, pn, user_list, neg_list, num_oids, locked);
    //
    //     CSeqDBVolEntry new_vol( new_volp );
    //     new_vol.SetStartAndEnd( num_oids );
    //     m_VolList.push_back( new_vol );
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.hpp:503-523
    // ```c++
    //     Uint8 GetVolumeSetLength() const
    //     {
    //         Uint8 vol_total = 0;
    //
    //         for(int index = 0; index < (int) m_VolList.size(); index++) {
    //             vol_total += m_VolList[index].Vol()->GetVolumeLength();
    //         }
    //
    //         return vol_total;
    //     }
    //
    //     int GetMaxLength() const
    //     {
    //         int max_len = 0;
    //
    //         for(int index = 0; index < (int) m_VolList.size(); index++) {
    //             max_len = max( max_len, m_VolList[index].Vol()->GetMaxLength());
    //         }
    //
    //         return max_len;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:984-988
    // ```c++
    // string CSeqDBImpl::GetTitle() const
    // {
    //     CHECK_MARKER();
    //     return x_FixString( m_Aliases.GetTitle(m_VolSet) );
    // }
    // ```
    // The number of sequences and the total length that NCBI takes from its
    // alias walk (seqdbalias.cpp:1189-1192, 1290-1293, 1575-1596) are the sums
    // over the volumes when no alias file has the keys NSEQ and LENGTH, or
    // when their values equal the sums below their nodes (`alias::check_sums`,
    // which rejects every other value).
    /// Opens the volumes of a resolved database (`name` is for messages).
    fn open_tree(tree: alias::AliasTree, name: &str, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let mut volumes = Vec::with_capacity(tree.volumes.len());
        let mut volume_starts = Vec::with_capacity(tree.volumes.len() + 1);
        let mut num_oids: u64 = 0;
        let mut total_length: u64 = 0;
        let mut max_length: u32 = 0;
        for volume_name in &tree.volumes {
            let volume = Volume::open(Path::new(volume_name), seqtype)?;
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
            max_length = max_length.max(volume.index().max_length);
            volumes.push(volume);
        }
        volume_starts.push(num_oids as u32);

        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:217-233
        // ```c++
        // 	bool version5 = true;
        // 	for(int i=0; i < volSet.GetNumVols(); i++) {
        // 		const CSeqDBVol * vol = volSet.GetVol(i);
        // 		string vol_lmdb_name = vol->GetLMDBFileName();
        // 		if(vol_lmdb_name == kEmptyStr) {
        // 			version5 = false;
        // 			if ((!m_LMDBEntrySet.empty())|| (!lmdb_vol_names.empty())) {
        // 				m_LMDBEntrySet.clear();
        // 				NCBI_THROW(CSeqDBException, eVersionErr, "DB list contains both Version 4 and Version 5 dbs");
        // 			}
        // 			continue;
        // 		}
        //
        // 		if(version5 == false) {
        // 			m_LMDBEntrySet.clear();
        // 			NCBI_THROW(CSeqDBException, eVersionErr, "DB list contains both Version 4 and Version 5 dbs");
        // 		}
        // ```
        // A volume counts as version 5 when its index names an LMDB file.
        let names_lmdb = |volume: &Volume| {
            volume
                .index()
                .lmdb_file
                .as_ref()
                .is_some_and(|name| !name.is_empty())
        };
        if let Some(first) = volumes.first() {
            let version5 = names_lmdb(first);
            if let Some(other) = volumes.iter().find(|volume| names_lmdb(volume) != version5) {
                return Err(SeqDbError::Format {
                    path: other.base().to_path_buf(),
                    message: "DB list contains both Version 4 and Version 5 dbs".to_string(),
                });
            }
        }

        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:85-89
        // ```c++
        //         CFile tf(fname);
        // 		Uint8 readMapSize = (tf.GetLength()/10000 + 1) *10000;
        //         if (!tf.Exists()) {
        //             NCBI_THROW( CSeqDBException, eFileErr, "File " + fname + " not found. If you renamed any BLAST database files, please use original file names, and makeblastdb to rename the database. If you deleted any BLAST database files, you need to recreate the database.");
        //         }
        // ```
        // NCBI opens the LMDB file of a version 5 database when it opens the
        // database (seqdblmdbset.cpp:235-262, 35-74) and compares the file's
        // volume list with the volumes: below, `lmdb_entries`.
        for volume in &volumes {
            if !names_lmdb(volume) {
                continue;
            }
            let lmdb_name = volume.index().lmdb_file.as_deref().unwrap_or_default();
            let lmdb_path = volume
                .base()
                .with_file_name(String::from_utf8_lossy(lmdb_name).as_ref());
            if !std::fs::metadata(&lmdb_path).is_ok_and(|meta| meta.is_file()) {
                return Err(SeqDbError::NotFound(format!(
                    "File {} not found. If you renamed any BLAST database files, please use \
                     original file names, and makeblastdb to rename the database. If you deleted \
                     any BLAST database files, you need to recreate the database.",
                    lmdb_path.display()
                )));
            }
        }

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
            path: PathBuf::from(name),
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
            alias_ranges: false,
            oid_list_set_up: AtomicBool::new(false),
            list_applied: false,
            name: name.to_string(),
            filter: None,
            search_path: None,
            reports_show_gis: false,
            positive_list: false,
            list_taxids: Vec::new(),
            name_records: Vec::new(),
            lookup_files: OnceLock::new(),
            lookup_isams: OnceLock::new(),
            mask_algorithms: OnceLock::new(),
        };
        // A volume list that NCBI opens and whose OIDs it then gets wrong
        // (`SeqDbError::Unsupported`) is refused by the lookups that need them, not here.
        if db.is_version5() {
            match db.lmdb_entries() {
                Ok(_) | Err(SeqDbError::Unsupported { .. }) => {}
                Err(error) => return Err(error),
            }
        }
        let mut root = tree.root;
        db.read_alias_tree(&mut root, &tree.volumes)?;
        Ok(db)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:94-139
    // ```c++
    //     m_OidListSetup = ! (m_Aliases.HasFilters() || gi_list || neg_list);
    //     ...
    //         if (m_NeedTotalsScan) {
    //             ...
    //             x_ScanTotals(true, & m_NumSeqs, & m_TotalLength,
    //                                & m_MaxLength, & m_MinLength, locked);
    //
    //         } else {
    //             m_NumSeqs     = x_GetNumSeqs();
    //             m_TotalLength = x_GetTotalLength();
    //             m_MaxLength   = x_GetMaxLength();
    //             m_MinLength   = x_GetMinLength();
    //
    //             // Do not bother scanning the db... it would be slow
    //             // FIXME: future implementation should probably have the
    //             // shortest length encoded in index file...
    //             if ( m_MinLength <= 0)  m_MinLength = BLAST_SEQSRC_MINLENGTH;
    //         }
    //         m_NumSeqsStats     = x_GetNumSeqsStats();
    //         m_TotalLengthStats = x_GetTotalLengthStats();
    // ```
    // `HasFilters` runs `ComputeMasks` (the OID ranges, `alias::compute_masks`) after the
    // volumes and the LMDB files are open. The walks of a value that does not convert fail
    // where NCBI's constructor reaches them, which depends on the list (`check_alias_values`).
    // `x_GetNumSeqs` and `x_GetNumSeqsStats` are the `int` of the walks' `Int8`
    // (seqdbimpl.cpp:948-949); the ranges make the OID list (`alias::range_mask`), which a
    // list later intersects (`restrict`), and the iteration range that its first build sets
    // (`alias::top_iteration_range`), inside which a search reads the list's OIDs.
    /// Reads the counts, the OID ranges and their errors from the alias tree `root` of the
    /// volumes `names` (the volume paths in the order of `self.volumes`).
    fn read_alias_tree(
        &mut self,
        root: &mut alias::AliasNode,
        names: &[String],
    ) -> Result<(), SeqDbError> {
        alias::compute_masks(root).map_err(SeqDbError::Ncbi)?;
        let at = |path: &str| {
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
        let mask = alias::range_mask(
            root,
            &|path| {
                let at = at(path);
                (self.volume_starts[at], self.volume_starts[at + 1])
            },
            num_oids,
        );
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
        self.alias_ranges = mask.is_some();
        if let Some(words) = mask {
            let range = alias::top_iteration_range(root, num_oids).unwrap_or((0, num_oids));
            let filter = filter::OidFilter::new(words, num_oids).with_range(range);
            self.num_seqs = filter.count();
            self.filter = Some(Box::new(filter));
            self.set_included_length(None);
        }
        self.alias_totals = totals;
        self.update_first_oid_length()?;
        Ok(())
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
                SeqType::Protein => u64::from(volume.protein_length(0)?),
                SeqType::Nucleotide => volume.nucleotide_length(0, 0)? as u64,
            });
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:126-133
    // ```c++
    // void
    // CSeqDBVol::x_OpenHdrFile(void) const {
    //     CFastMutexGuard mtx_gurad(m_MtxHdr);
    //     if (!m_HdrFileOpened && m_Idx->GetNumOIDs() != 0) {
    //         m_Hdr.Reset(new CSeqDBHdrFile(m_Atlas, m_VolName, (m_IsAA?'p':'n')));
    //     }
    //     m_HdrFileOpened = true;
    // }
    // ```
    // NCBI opens a header file at the first read of a header of its volume
    // (`x_GetHdrAsn1Binary`, seqdbvol.cpp:2302-2310): a search program reads
    // the first one after its search, so a header file that is missing or
    // cannot be opened fails there, not here (`Volume::lazy_file`).
    /// Whether the database holds parsed Seq-ids, which a search program asks right after
    /// opening it: the first header of every volume is decoded with the ordinal-id decoder
    /// (`deflines`), whose error for any other Seq-id is `SeqDbError::Unsupported`
    /// (`SeqDbError::is_parsed_seqids`). `Ok(true)` when some volume's first header holds
    /// parsed Seq-ids; any other error of a first header is returned. A volume whose header
    /// file cannot be opened (missing, or an `open` that fails) is not checked: NCBI opens a
    /// header file at its first read, and the error comes there.
    pub fn has_parsed_seqids(&self) -> Result<bool, SeqDbError> {
        let mut parsed = false;
        for (at, volume) in self.volumes.iter().enumerate() {
            if volume.num_oids() > 0 {
                match self.first_deflines(at) {
                    Ok(_) | Err(SeqDbError::NotFound(_) | SeqDbError::Open { .. }) => {}
                    Err(error) if error.is_parsed_seqids() => parsed = true,
                    Err(error) => return Err(error),
                }
            }
        }
        Ok(parsed)
    }

    /// `deflines` of the first sequence of volume `at` (non-empty), whether a list includes it
    /// or not, and with every definition line: the check reads the header as stored, before a
    /// taxid list's `DeflineRule` (`Volume::stored_deflines`), which would refuse a first sequence
    /// that the list leaves out. Its messages name the sequence's OID.
    fn first_deflines(&self, at: usize) -> Result<Vec<OrdinalDefline>, SeqDbError> {
        let oid = self.volume_starts[at];
        self.ordinal_deflines(at, oid, self.volumes[at].stored_deflines(0))
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
    // ```c++
    // CRef< CSeq_entry > CMaskBDBReader::GetNextSequence()
    // {
    //     CRef< CSeq_entry > entry( null );
    //     if( oid_ < seqdb_->GetNumOIDs() ) {
    //         CRef< CBioseq > seq = seqdb_->GetBioseq( oid_++ );
    // ```
    /// Every OID of the volumes as an ordinal id ([`SeqDb::num_oids`] of them): the alias
    /// files' OID ranges (`FIRST_OID`, `LAST_OID`) no longer leave any out. For a reader that
    /// loops `GetNumOIDs` and reads each OID with `GetBioseq(oid)`, which reads an OID outside
    /// the ranges too (dustmasker's `CMaskBDBReader`). Only before a list
    /// ([`SeqDb::restrict`]): the alias keys of lists are refused when the database opens.
    fn include_every_oid(&mut self) {
        debug_assert!(!self.list_applied, "include_every_oid after a list");
        self.filter = None;
        self.set_included_length(None);
        self.num_seqs = self.num_oids();
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:117-139
    // ```c++
    //         if (m_NeedTotalsScan) {
    //             CSeqDBLockHold locked(m_Atlas);
    //
    //             // This is a whole-database scan; it's always done in
    //             // approximate length mode.
    //             x_ScanTotals(true, & m_NumSeqs, & m_TotalLength,
    //                                & m_MaxLength, & m_MinLength, locked);
    //
    //         } else {
    //             m_NumSeqs     = x_GetNumSeqs();
    //             m_TotalLength = x_GetTotalLength();
    //             m_MaxLength   = x_GetMaxLength();
    //             m_MinLength   = x_GetMinLength();
    //             ...
    //         }
    //         m_NumSeqsStats     = x_GetNumSeqsStats();
    //         m_TotalLengthStats = x_GetTotalLengthStats();
    // ```
    /// Opens the database `name` for a reader that loops `GetNumOIDs` and reads each OID with
    /// `GetBioseq(oid)` (dustmasker's `CMaskBDBReader`, which makes `CSeqDB(name, type)`):
    /// [`SeqDb::open_search_in`], the errors of the walks of the alias files' counts that NCBI's
    /// constructor makes (those of [`SeqDb::check_alias_values`] without a list), and
    /// `include_every_oid`. LOSAT's refusals of the values that only a search reads
    /// differently are not made, as such a reader reads none of them: a negative `NSEQ`, a
    /// `LENGTH` of 0 or of 2^63 or more (NCBI 2.17.0 dustmasker reads them as it reads the
    /// database without the alias file: review 2 of line RG item K1, F5).
    pub fn open_every_oid_in(
        name: &str,
        seqtype: SeqType,
        path: &SearchPath,
    ) -> Result<Self, SeqDbError> {
        let mut db = Self::open_search_in(name, seqtype, path)?;
        db.walk_errors(UNLISTED_WALK_KEYS)?;
        db.include_every_oid();
        Ok(db)
    }

    /// `CSeqDB::GetNumSeqs` of the search: the walk of `NSEQ` (the volumes' OIDs below the
    /// alias files without the key) without a list, the included OIDs with one. It is not the
    /// number of included OIDs ([`SeqDb::num_seqs`]) when an alias file has `NSEQ` or a range.
    #[cfg(test)]
    pub fn ncbi_num_seqs(&self) -> i32 {
        self.counts.num_seqs
    }

    /// `CSeqDB::GetTotalLength` of the search: the walk of `LENGTH` (the volumes' lengths below
    /// the alias files without the key) without a list, NCBI's scanned total of the included
    /// OIDs with one (approximate for nucleotides: `restrict_taxids`).
    pub fn total_length(&self) -> u64 {
        self.counts.total_length
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
    #[cfg(test)]
    pub fn num_seqs_stats(&self) -> i32 {
        self.stats_num_seqs
    }

    /// `CSeqDB::GetTotalLengthStats`: the walk of `STATS_TOTLEN`, 0 without the key.
    #[cfg(test)]
    pub fn total_length_stats(&self) -> u64 {
        self.stats_total_length
    }

    // NCBI's reports open the database again without the list (`CBlastFormatUtil::
    // GetBlastDbInfo` reads a `CSeqDB` of the name alone): the `Database:` lines of the
    // pairwise report and the database counts of XML and JSON are the walks of the alias tree.
    /// The number of sequences the reports print: `GetNumSeqs` of the database without a list
    /// (`check_alias_values` refuses a negative one).
    pub fn printed_num_seqs(&self) -> u32 {
        self.printed_num_seqs.max(0) as u32
    }

    /// The total length the reports print: `GetTotalLength` of the database without a list.
    pub fn printed_total_length(&self) -> u64 {
        self.printed_total_length
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:206-215
    // ```c++
    //     if (! BlastSeqSrcGetNumSeqs(m_InternalData->m_SeqSrc->GetPointer())) {
    //         string msg =
    //             "Filtering resulted in an empty database.";
    //
    //         m_Messages.AddMessageAllQueries(eBlastSevWarning,
    //                                         kBlastMessageNoContext,
    //                                         msg);
    //     }
    // ```
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
    // ```c
    //    Int8 total_length = -1;
    //    Int4 num_seqs = -1;
    //
    //    if (seq_src) {
    //       total_length = BlastSeqSrcGetTotLenStats(seq_src);
    //       if (total_length <= 0)
    //           total_length = BlastSeqSrcGetTotLen(seq_src);
    //       ...
    //       if (total_length > 0) {
    //           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
    //           if (num_seqs <= 0)
    //               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
    //       } else {
    //           /* Not a database search; each subject sequence is considered
    //              individually */
    //           Int4 oid = 0;  /* Get length of first sequence. */
    //           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
    //               total_length = -1;
    //               num_seqs = -1;
    //           }
    //           num_seqs = 1;
    //       }
    //    }
    // ```
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

    /// The total length and the number of sequences that the statistics of a search read
    /// (`BLAST_CalcEffLengths`, `BlastSetup_GetScoreBlock`): `STATS_TOTLEN` or `GetTotalLength`,
    /// and `STATS_NSEQ` or `GetNumSeqs`; the length of OID 0 and one sequence where the total
    /// is not positive (a list that kept no sequence). The number is 0 for an alias file's
    /// `NSEQ 0`, which NCBI's average subject length divides by (`SeqDb::average_length`).
    pub fn statistics_numbers(&self) -> (u64, u32) {
        let total = self.statistics_total_length();
        if total > 0 {
            let stats = self.stats_num_seqs;
            let num_seqs = if stats > 0 {
                stats
            } else {
                self.counts.num_seqs
            };
            return (total as u64, num_seqs.max(0) as u32);
        }
        match self.first_oid_length {
            Some(length) => (length, 1),
            None => (0, 0),
        }
    }

    /// Whether the alias files give an OID range (`FIRST_OID`, `LAST_OID`) or a key that
    /// changes the statistics or the shortest length (`STATS_NSEQ`, `STATS_TOTLEN`,
    /// `MIN_SEQ_LENGTH`): RPS-BLAST, which searches the volumes' profiles by OID, refuses them
    /// (`rps_hooks::check_alias_keys`).
    pub fn alias_ranges_or_statistics(&self) -> bool {
        self.alias_ranges
            || self.stats_num_seqs != 0
            || self.stats_total_length != 0
            || !matches!(self.alias_totals.min_seq_length, Ok(0))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:164-173
    // ```c++
    //     Int8 total_length = BlastSeqSrcGetTotLenStats(seqsrc);
    //     if (total_length <= 0)
    //         total_length = BlastSeqSrcGetTotLen(seqsrc);
    // ```
    /// Whether the total length of the statistics is not positive before their fallback to the
    /// length of OID 0 (`STATS_TOTLEN` and `GetTotalLength` are not: a list that kept no
    /// sequence, without `STATS_TOTLEN`): the effective search space that NCBI's BLAST archive
    /// records is then computed with a database length of 0 (`note_archive_database_filter`).
    pub fn statistics_length_is_empty(&self) -> bool {
        self.statistics_total_length() <= 0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:932-936 (at
    // `statistics_numbers`)
    // ```c
    //           Int4 oid = 0;  /* Get length of first sequence. */
    //           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251 (at
    // `common::database::PreliminaryReadError`): the `CSeqDBException` of that read is
    // `BLAST engine error: Database memory map file error` (exit status 3; psiblast's
    // `BLAST Database error: OID not found`, exit status 2).
    /// Whether the setup of the first query batch with a valid context fails: the statistics'
    /// total is not positive ([`SeqDb::statistics_length_is_empty`]), so NCBI reads the length of
    /// OID 0 (`GetSeqLength(0)`, which asks the volumes, not the list), and the volumes have no
    /// OID. An alias file's `LENGTH` or `STATS_TOTLEN` over such volumes makes the total positive:
    /// NCBI then reads no sequence and reports no hit (probes p4c_dp, p4c_dn of line RG D1a-2).
    pub fn first_oid_read_fails(&self) -> bool {
        self.statistics_length_is_empty() && self.num_oids() == 0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-975
    // ```c
    //    if (sbp->gbp) {
    //    ...
    //    } else {
    //        min_subject_length = (Int4) (total_length/num_seqs);
    //    }
    // ```
    /// Whether the statistics read 0 sequences of a positive total length
    /// ([`SeqDb::statistics_numbers`]: an alias file's `NSEQ 0` without `STATS_NSEQ`, its `LENGTH`
    /// or `STATS_TOTLEN` over volumes without sequences, or a list that kept no sequence under a
    /// `STATS_TOTLEN`). The programs without the Gumbel block
    /// (BLASTN, TBLASTX, and the protein programs with `-ungapped`, `-searchsp` or `OLD_FSC`)
    /// divide the total by it, and NCBI BLAST+ 2.17.0 crashes (SIGFPE).
    pub fn statistics_count_is_zero(&self) -> bool {
        let (total, count) = self.statistics_numbers();
        total > 0 && count == 0
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:199-209
    // ```c++
    // static Int4
    // s_SeqDbGetAvgLength(void* seqdb_handle, void* ignoreme)
    // {
    //    Int8 total_length = s_SeqDbGetTotLen(seqdb_handle, ignoreme);
    //    Int4 num_seqs = MAX(1, s_SeqDbGetNumSeqs(seqdb_handle, ignoreme));
    //
    //    return (Int4) (total_length/num_seqs);
    // }
    // ```
    /// `BlastSeqSrcGetAvgSeqLen`: `GetTotalLength / MAX(1, GetNumSeqs)` (the alias files'
    /// `STATS_*` keys are not read), as NCBI's `Int4`.
    pub fn average_length(&self) -> i32 {
        let total = self.counts.total_length as i64;
        (total / i64::from(self.counts.num_seqs.max(1))) as i32
    }

    /// The errors and LOSAT's refusals of the alias files' values that NCBI's `CSeqDB` of the
    /// search meets (`read_alias_tree`): without a list the walks of `NSEQ`, `LENGTH`,
    /// `MIN_SEQ_LENGTH`, `STATS_NSEQ`, `STATS_TOTLEN` in that order; after a list
    /// ([`SeqDb::restrict`]), which makes NCBI scan the totals, those of `STATS_NSEQ` and
    /// `STATS_TOTLEN` (the reports' database, opened again without the list, reads the others:
    /// [`SeqDb::read_name_records`]). NCBI's error is `SeqDbError::Ncbi`. [`SeqDb::open_in`]
    /// calls this; a search program opens with [`SeqDb::open_search_in`] and calls it after its
    /// list.
    pub fn check_alias_values(&self) -> Result<(), SeqDbError> {
        use alias::WalkKey;
        let keys: &[WalkKey] = if self.list_applied {
            &[WalkKey::StatsNseq, WalkKey::StatsTotlen]
        } else {
            UNLISTED_WALK_KEYS
        };
        self.walk_errors(keys)?;
        // Values NCBI reads in ways LOSAT does not reproduce yet (line RG, item D1a-2).
        let refuse = |message: &str| {
            Err(SeqDbError::Unsupported {
                path: PathBuf::from(self.alias_name()),
                message: format!("{message}; LOSAT does not read such alias files yet"),
            })
        };
        if self.counts.num_seqs < 0 {
            return refuse(&format!(
                "the number of sequences of the alias files (NSEQ) is {} as NCBI's int",
                self.counts.num_seqs
            ));
        }
        // A total of 2^63 or more is negative as NCBI's Int8: NCBI's search reads that negative
        // `db_length` the usual way (blast_setup.c:925-935, the cutoffs of blast_engine.c:1454),
        // not per subject; LOSAT's statistics take no negative total (line RG D1a-3).
        let negative_total = (self.counts.total_length as i64) < 0 && self.num_seqs > 0;
        if !self.list_applied && (self.statistics_per_subject() || negative_total) {
            return refuse(&format!(
                "the total length of the alias files (LENGTH) is {} as NCBI's Int8 while the \
                 database includes sequences",
                self.counts.total_length as i64
            ));
        }
        if !self.list_applied {
            self.check_printed_counts()?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443
    // ```c
    //        if (db_length == 0) {
    //            /* This is not a database search, hence need to recalculate and save
    //             the effective search spaces and length adjustments for all
    //             queries based on the length of the current single subject
    //             sequence. */
    //            if ((status = BLAST_OneSubjectUpdateParameters(program_number,
    //                           seq_arg.seq->length, score_options, query_info,
    //                           sbp, hit_params, word_params,
    //                           eff_len_options)) != 0)
    // ```
    // `db_length` is `BlastSeqSrcGetTotLen`, NCBI's `GetTotalLength` (`counts`): an alias file's
    // `LENGTH 0`, or a list whose scan of the totals visits no sequence (the walk of
    // `OidFilter::walk` stops before the listed OIDs that the search reads: line RG D1a-2,
    // survey A, the iteration range of the only alias file of a `-db`).
    /// Whether NCBI's search computes the statistics of each subject on its own: `GetTotalLength`
    /// is 0 while the search reads sequences. LOSAT's database searches do not (an alias file's
    /// `LENGTH 0` is refused by [`SeqDb::check_alias_values`], a list's scan by the search's
    /// opener). Every included sequence counts as read here, while NCBI's TBLASTN and TBLASTX
    /// skip the subjects shorter than `min_subj_seq_length` before this check
    /// (blast_engine.c:1429): a database whose searched subjects are all that short is refused
    /// though NCBI searches it the usual way, until the statistics of each subject are ported
    /// (line RG D1a-3, which removes this refusal).
    pub fn statistics_per_subject(&self) -> bool {
        self.counts.total_length == 0 && self.num_seqs > 0
    }

    /// The first exception of the walks of `keys`, in that order.
    fn walk_errors(&self, keys: &[alias::WalkKey]) -> Result<(), SeqDbError> {
        match keys.iter().find_map(|&key| self.alias_totals.error(key)) {
            Some(exception) => Err(SeqDbError::Ncbi(exception.clone())),
            None => Ok(()),
        }
    }

    /// LOSAT's refusals of the counts the reports print (line RG, item D1a-2): a negative
    /// `GetNumSeqs`, and a total length of 10^18 or more, which NCBI's commas print with the
    /// upper digits lost.
    fn check_printed_counts(&self) -> Result<(), SeqDbError> {
        let refuse = |message: String| {
            Err(SeqDbError::Unsupported {
                path: PathBuf::from(self.alias_name()),
                message: format!("{message}; LOSAT does not read such alias files yet"),
            })
        };
        if self.printed_num_seqs < 0 {
            return refuse(format!(
                "the number of sequences of the alias files (NSEQ) is {} as NCBI's int",
                self.printed_num_seqs
            ));
        }
        if self.printed_total_length >= 1_000_000_000_000_000_000 {
            return refuse(format!(
                "the total length of the alias files (LENGTH) is {}, which NCBI prints with \
                 digits lost",
                self.printed_total_length
            ));
        }
        Ok(())
    }

    /// The name of the database for LOSAT's refusals of its alias files.
    fn alias_name(&self) -> &str {
        &self.name
    }

    /// Title of the database.
    pub fn title(&self) -> &[u8] {
        &self.title
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1010-1013
    // ```c++
    //     string date;
    //
    //     for(int i = 0; i < m_VolSet.GetNumVols(); i++) {
    //         string d = x_FixString( m_VolSet.GetVol(i)->GetDate() );
    // ```
    /// Date of the database: the latest date of its volumes (`date`).
    pub fn date(&self) -> &str {
        &self.date
    }

    /// The volumes, in the order of their ordinal ids.
    pub fn volumes(&self) -> &[Volume] {
        &self.volumes
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.hpp:245-275
    // ```c++
    //     const CSeqDBVol * FindVol(int oid, int & vol_oid, int & vol_idx) const
    //     {
    //         ...
    //         for(int index = 0; index < (int) m_VolList.size(); index++) {
    //             if ((m_VolList[index].OIDStart() <= oid) &&
    //                 (m_VolList[index].OIDEnd()   >  oid)) {
    //                 ...
    //                 vol_oid = oid - m_VolList[index].OIDStart();
    //                 vol_idx = index;
    //
    //                 return m_VolList[index].Vol();
    //             }
    //         }
    //
    //         return NULL;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:416-425
    // ```c++
    //         if (const CSeqDBVol * vol = m_VolSet.FindVol(oid, vol_oid)) {
    //             return vol->GetSeqLengthExact(vol_oid);
    //         }
    //     }
    //
    //     NCBI_THROW(CSeqDBException, eArgErr, CSeqDB::kOidNotFound);
    // ```
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
    // ```c++
    //     // Case 1: Iteration's End.
    //
    //     if (*state_obj >= m_RestrictEnd) {
    //         begin_chunk = 0;
    //         end_chunk   = 0;
    //         return CSeqDB::eOidRange;
    //     }
    //
    //     begin_chunk = * state_obj;
    //     ...
    //         end_chunk = begin_chunk + oid_size;
    //     }
    //
    //     if (end_chunk > m_RestrictEnd) {
    //         end_chunk = m_RestrictEnd;
    //     }
    //     *state_obj = end_chunk;
    // ```
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

    /// Reads the sequences of `range` with one positioned read.
    pub fn read_block(&self, range: BlockRange) -> Result<SequenceBlock, SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        volume.read_block(range.first, range.end, range.first_ordinal)
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
        volume.x_nucleotide_check(range.first, range.end, range.first_ordinal)
    }

    /// Reads the one sequence `ordinal` (for traceback and reports).
    pub fn read_sequence(&self, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
        self.read_sequence_at(self.find_volume(ordinal)?)
    }

    /// Reads the one sequence at `at` (a block of one sequence whose first ordinal id is the id
    /// `at` was found by).
    pub fn read_sequence_at(&self, at: SeqLocation) -> Result<SequenceBlock, SeqDbError> {
        self.volumes[at.volume].read_block(at.vol_oid, at.vol_oid + 1, at.id)
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

    /// The bytes of the volume's sequence file after the sentinel byte of protein `ordinal`
    /// (`Volume::protein_bytes_after`).
    pub fn protein_bytes_after(&self, ordinal: u32) -> Result<u64, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].protein_bytes_after(at.vol_oid)
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `read_block` with the bytes in a buffer that glibc maps
    /// from 1 MiB on (`Volume::x_read_block_mapped`).
    pub fn x_read_block_mapped(&self, range: BlockRange) -> Result<SequenceBlock, SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        volume.x_read_block_mapped(range.first, range.end, range.first_ordinal)
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP): `read_block` with the bytes in `buf`
    /// (`Volume::x_read_block_into`).
    #[cfg(any(
        test,
        all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )
    ))]
    pub fn x_read_block_into(
        &self,
        range: BlockRange,
        buf: Vec<u8>,
    ) -> Result<SequenceBlock, SeqDbError> {
        let Some(volume) = self.volumes.get(range.volume) else {
            return Err(SeqDbError::Range(format!(
                "volume {} is outside the database ({} volumes)",
                range.volume,
                self.volumes.len()
            )));
        };
        volume.x_read_block_into(range.first, range.end, range.first_ordinal, buf)
    }

    /// EXPERIMENT (LOSAT_X_BNMAPPED): `read_sequence` with the bytes in a buffer that glibc maps
    /// from 1 MiB on.
    pub fn x_read_sequence_mapped(&self, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].x_read_block_mapped(at.vol_oid, at.vol_oid + 1, ordinal)
    }

    /// The length in bases of nucleotide `ordinal`, read from the index and the last byte of
    /// its packed bases (no letters).
    pub fn nucleotide_length(&self, ordinal: u32) -> Result<usize, SeqDbError> {
        self.nucleotide_length_at(self.find_volume(ordinal)?)
    }

    /// [`SeqDb::nucleotide_length`] of the nucleotide at `at`.
    pub fn nucleotide_length_at(&self, at: SeqLocation) -> Result<usize, SeqDbError> {
        self.volumes[at.volume].nucleotide_length(at.vol_oid, at.id)
    }

    /// The IUPAC letters `[from, to)` of nucleotide `ordinal`: the bases of the range and the
    /// ambiguity runs that overlap it, read without the rest of the sequence.
    pub fn nucleotide_iupac_range(
        &self,
        ordinal: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        self.nucleotide_iupac_range_at(self.find_volume(ordinal)?, from, to)
    }

    /// [`SeqDb::nucleotide_iupac_range`] of the nucleotide at `at`.
    pub fn nucleotide_iupac_range_at(
        &self,
        at: SeqLocation,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        self.volumes[at.volume].nucleotide_iupac_range(at.vol_oid, at.id, from, to)
    }

    /// The length and the ambiguity runs of nucleotide `ordinal`, read once for any number of
    /// [`SeqDb::nucleotide_iupac_range_with`] calls on it.
    pub fn nucleotide_ambiguity(&self, ordinal: u32) -> Result<NucleotideAmbiguity, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].nucleotide_ambiguity(at.vol_oid, ordinal)
    }

    /// [`SeqDb::nucleotide_iupac_range`] with the length and the runs of
    /// [`SeqDb::nucleotide_ambiguity`] for the same `ordinal`: one read of the range's bytes.
    pub fn nucleotide_iupac_range_with(
        &self,
        ordinal: u32,
        from: usize,
        to: usize,
        ambiguity: &NucleotideAmbiguity,
    ) -> Result<Vec<u8>, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume]
            .nucleotide_iupac_range_with(at.vol_oid, ordinal, from, to, ambiguity)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2222-2241
    // ```c++
    //     if (adjust_oids && bdls.NotEmpty() && m_VolStart) {
    //         NON_CONST_ITERATE(list< CRef<CBlast_def_line> >, dl, bdls->Set()) {
    //             if (! (**dl).CanGetSeqid()) {
    //                 continue;
    //             }
    //
    //             NON_CONST_ITERATE(list< CRef<CSeq_id> >, id, (*dl)->SetSeqid()) {
    //                 CSeq_id & seqid = **id;
    //
    //                 if (seqid.Which() == CSeq_id::e_General) {
    //                     CDbtag & dbt = seqid.SetGeneral();
    //
    //                     if (dbt.GetDb() == "BL_ORD_ID") {
    //                         int vol_oid = dbt.GetTag().GetId();
    //                         dbt.SetTag().SetId(m_VolStart + vol_oid);
    //                         if (changed) {
    //                             *changed = true;
    //                         }
    //                     }
    //                 }
    // ```
    /// The definition lines of `ordinal`: the title and `n` of
    /// `gnl|BL_ORD_ID|n`, where `n` is the stored number plus the first
    /// OID of the volume, as NCBI reports it.
    pub fn deflines(&self, ordinal: u32) -> Result<Vec<OrdinalDefline>, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.ordinal_deflines(
            at.volume,
            ordinal,
            self.volumes[at.volume].deflines(at.vol_oid),
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2221-2240
    // ```c++
    //     if (adjust_oids && bdls.NotEmpty() && m_VolStart) {
    //         NON_CONST_ITERATE(list< CRef<CBlast_def_line> >, dl, bdls->Set()) {
    //             if (! (**dl).CanGetSeqid()) {
    //                 continue;
    //             }
    //
    //             NON_CONST_ITERATE(list< CRef<CSeq_id> >, id, (*dl)->SetSeqid()) {
    //                 CSeq_id & seqid = **id;
    //
    //                 if (seqid.Which() == CSeq_id::e_General) {
    //                     CDbtag & dbt = seqid.SetGeneral();
    //
    //                     if (dbt.GetDb() == "BL_ORD_ID") {
    //                         int vol_oid = dbt.GetTag().GetId();
    //                         dbt.SetTag().SetId(m_VolStart + vol_oid);
    //                         if (changed) {
    //                             *changed = true;
    //                         }
    //                     }
    // ```
    /// The header of `ordinal` as a report reads it: the ordinal-id definition lines
    /// (`deflines`), or every member and Seq-id (`defline_set`) of a record with parsed
    /// Seq-ids. A record without parsed Seq-ids is read as before, with the same errors.
    pub fn read_header(&self, ordinal: u32) -> Result<SequenceHeader, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.sequence_header(
            at.volume,
            ordinal,
            self.volumes[at.volume].header(at.vol_oid),
        )
    }

    /// The `SequenceHeader` of the header `header` of the sequence `ordinal` in `volume`.
    fn sequence_header(
        &self,
        volume: usize,
        ordinal: u32,
        header: Result<VolumeHeader, SeqDbError>,
    ) -> Result<SequenceHeader, SeqDbError> {
        match header? {
            VolumeHeader::Ordinal(deflines) => self
                .ordinal_deflines(volume, ordinal, Ok(deflines))
                .map(SequenceHeader::Ordinal),
            VolumeHeader::Parsed(set) => self
                .ordinal_set(volume, ordinal, set)
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

    /// [`SeqDb::defline_set`] of the sequence at `at`.
    pub fn defline_set_at(&self, at: SeqLocation) -> Result<DeflineSet, SeqDbError> {
        let set = self.volumes[at.volume].defline_set(at.vol_oid)?;
        self.ordinal_set(at.volume, at.id, set)
    }

    /// `set`, the definition lines of the sequence `ordinal` (the id its messages name) in
    /// `volume`, with the tag of a `general` id `BL_ORD_ID` made its OID in the database
    /// (`defline_set`).
    fn ordinal_set(
        &self,
        volume: usize,
        ordinal: u32,
        mut set: DeflineSet,
    ) -> Result<DeflineSet, SeqDbError> {
        let volume_start = self.volume_starts[volume];
        if volume_start == 0 {
            return Ok(set);
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
            let format = |message: String| SeqDbError::Format {
                path: self.volumes[volume].base().to_path_buf(),
                message,
            };
            // `GetId` of a string tag throws in NCBI (an invalid CHOICE selection).
            let ObjectId::Id(tag) = dbtag.tag else {
                return Err(format(format!(
                    "BL_ORD_ID of sequence {ordinal} has a string tag"
                )));
            };
            let Ok(adjusted) = i32::try_from(i64::from(volume_start) + i64::from(tag)) else {
                return Err(format(format!(
                    "BL_ORD_ID of sequence {ordinal} passes a signed 32-bit id"
                )));
            };
            dbtag.tag = ObjectId::Id(adjusted);
        }
        Ok(set)
    }

    /// The Blast-def-line `taxid` of each definition line of `ordinal`, 0 for none: read
    /// apart from `deflines`, for the XML2 and JSON formats, the prefetch's taxonomy lookups and
    /// an empty stored title (`HeaderReader::stored_taxonomy_title`).
    pub fn defline_taxids(&self, ordinal: u32) -> Result<Vec<i32>, SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].defline_taxids(at.vol_oid)
    }

    /// The definition lines of a volume's sequence as `deflines` gives them for `ordinal`.
    fn ordinal_deflines(
        &self,
        volume: usize,
        ordinal: u32,
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
                        message: format!(
                            "BL_ORD_ID of sequence {ordinal} passes a signed 32-bit id"
                        ),
                    });
                }
                Ok(OrdinalDefline {
                    title: defline.title,
                    bl_ord_id: bl_ord_id as u32,
                    replaced_title_bytes: defline.replaced_title_bytes,
                    taxid: defline.taxid,
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
    pub fn visit_deflines(
        &self,
        range: BlockRange,
        visit: impl FnMut(u32, Result<Vec<OrdinalDefline>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_deflines_as(range, false, visit)
    }

    /// `visit_deflines`, with the `taxid` of each definition line decoded (the decode of
    /// `defline_taxids`, which also gives the titles): the titles and the taxids of a range in one
    /// pass over its headers. A record whose taxid does not decode is an `Err`, as in
    /// `defline_taxids`.
    pub fn visit_deflines_taxids(
        &self,
        range: BlockRange,
        visit: impl FnMut(u32, Result<Vec<OrdinalDefline>, SeqDbError>) -> bool,
    ) -> bool {
        self.visit_deflines_as(range, true, visit)
    }

    fn visit_deflines_as(
        &self,
        range: BlockRange,
        want_taxid: bool,
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
            want_taxid,
            &mut |oid, deflines| {
                let Some(ordinal) = self.range_ordinal(range, oid) else {
                    return true;
                };
                let deflines = self.ordinal_deflines(range.volume, ordinal, deflines);
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
                    visit(ordinal, self.sequence_header(range.volume, ordinal, header))
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
    /// `Defline::taxid`: 0 unless the definition lines were read with their taxids
    /// ([`SeqDb::visit_deflines_taxids`]).
    pub taxid: i32,
}
