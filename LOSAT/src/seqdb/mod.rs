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
//! ISAM files of `-parse_seqids`, column files (mask data), group alias files,
//! and alias keys that filter sequences or replace the totals. A database
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
pub mod date;
pub mod file;
pub mod header;
pub mod index;
pub mod taxdb;
pub mod volume;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod x_dbidxread_tests;

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use crate::seqid::model::ObjectId;
use crate::seqid::{DeflineSet, ParsedSeqId};

pub use alias::SearchPath;
pub use header::Defline;
pub use volume::{
    AmbiguityRun, NucleotideAmbiguity, SequenceBlock, Volume, BLASTNA_SENTINEL, FENCE_SENTRY,
};

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
    /// A valid database that needs something LOSAT does not read.
    Unsupported { path: PathBuf, message: String },
    /// A request outside the database (ordinal id, block range).
    Range(String),
}

impl std::fmt::Display for SeqDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqDbError::NotFound(message) | SeqDbError::Range(message) => f.write_str(message),
            SeqDbError::Open { path, source } | SeqDbError::Io { path, source } => {
                write!(f, "{}: {}", path.display(), source)
            }
            SeqDbError::Format { path, message } => write!(f, "{}: {}", path.display(), message),
            SeqDbError::Unsupported { path, message } => {
                write!(f, "{}: unsupported database: {}", path.display(), message)
            }
        }
    }
}

impl SeqDbError {
    /// Whether this is the error of the ordinal-id decoders for a record with a Seq-id other
    /// than `gnl|BL_ORD_ID|n` (a database built with `-parse_seqids`; `SeqDb::read_header`).
    pub fn is_parsed_seqids(&self) -> bool {
        matches!(self, SeqDbError::Unsupported { message, .. } if message == header::PARSED_SEQIDS)
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
const NCBI_FORMAT_MESSAGES: [&str; 4] = [
    // seqdbfile.cpp:181-185
    "Error: Not a valid version 4 or 5 database.",
    // seqdbfile.cpp:222-226
    "Error: requested sequence type does not match DB.",
    // seqdblmdbset.cpp:224, 231
    "DB list contains both Version 4 and Version 5 dbs",
    // seqdbalias.cpp:810-814
    "Illegal configuration: DB alias files are mutually recursive.",
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
                || message.starts_with(NCBI_EMPTY_DBLIST_PREFIX))
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
            SeqDbError::Io { .. } | SeqDbError::Unsupported { .. } | SeqDbError::Range(_) => None,
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
    // that file when the search reads its first sequence.
    /// The text NCBI BLAST+ gives when this error is raised while the search
    /// reads the sequences of the database: a search program prints
    /// `BLAST engine error: <text>` and ends with exit status 3. `None`: the
    /// error is LOSAT's ([`Self::ncbi_message`] is `None` too).
    pub fn ncbi_search_message(&self) -> Option<&'static str> {
        const TOO_MANY: &str = "Too many open files, please raise the open file limit";
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
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3153-3172
// ```c++
//     SIZE_TYPE beg = 0;
//     if (where == NStr::eTrunc_Begin  ||  where == NStr::eTrunc_Both) {
//         _ASSERT(beg < length);
//         while ( isspace((unsigned char) str[beg]) ) {
//             if (++beg == length) {
//                 return empty_str;
//             }
//         }
//     }
//     SIZE_TYPE end = length;
//     if ( where == NStr::eTrunc_End  ||  where == NStr::eTrunc_Both ) {
//         _ASSERT(beg < end);
//         while (isspace((unsigned char) str[--end])) {
//             if (beg == end) {
//                 return empty_str;
//             }
//         }
//         _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
//         ++end;
//     }
// ```
// `isspace` in the "C" locale: space, tab, line feed, vertical tab, form feed
// and carriage return. The reader itself (`CSeqDB`, `SeqDb::open_in`) takes
// the name as it is given: a name of spaces alone is "No database names were
// passed to CSeqDB::CSeqDB()." there, while the programs, which trim first,
// report "Database name is required.".
/// The database name that a BLAST program hands to the reader for the value
/// of its `-db` argument: the value without the white space at its ends.
/// A search program passes the result to [`SeqDb::open_in`].
pub fn db_argument_name(value: &str) -> &str {
    value.trim_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r'])
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

/// An open database: one volume, or the volumes of an alias file.
#[derive(Debug)]
pub struct SeqDb {
    seqtype: SeqType,
    volumes: Vec<Volume>,
    /// Ordinal id of the first sequence of each volume, then the total.
    volume_starts: Vec<u32>,
    title: Vec<u8>,
    date: String,
    total_length: u64,
    max_length: u32,
    /// The search path the database was found on (`open_in`), for files NCBI looks for
    /// beside it (`taxdb.bti`).
    search_path: Option<SearchPath>,
}

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
    pub fn open(name: &str, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let cwd = std::env::current_dir().unwrap_or_default();
        Self::open_in(name, seqtype, &SearchPath::new(&cwd, None, None)?)
    }

    /// Opens the database that a search program's `-db NAME` names, as NCBI
    /// BLAST+ finds it: `name` is looked for in the directories of `path`
    /// (`alias::resolve_in`). The name is a list in NCBI's syntax (spaces
    /// separate names, double quotes keep one name together), of which LOSAT
    /// accepts one name; it carries no file extension.
    ///
    /// A database that stores parsed Seq-ids opens; `check_ordinal_ids` tells
    /// it from one that holds only the ordinal ids.
    pub fn open_in(name: &str, seqtype: SeqType, path: &SearchPath) -> Result<Self, SeqDbError> {
        let mut db = Self::open_tree(alias::resolve_in(name, seqtype, path)?, name, seqtype)?;
        db.search_path = Some(path.clone());
        Ok(db)
    }

    /// Whether `taxdb.bti` is found on the search path the database was opened with
    /// (`SeqDB_ResolveDbPath("taxdb.bti")`, data4xml2format.cpp:x_InitCommon): NCBI's XML2 and
    /// JSON formats then write the scientific name of a taxid.
    pub fn taxdb_found(&self) -> bool {
        self.search_path
            .as_ref()
            .is_some_and(|path| path.resolve_file("taxdb.bti").is_some())
    }

    /// The search path the database was opened with (`open_in`), none for a database opened
    /// by another way.
    pub fn search_path(&self) -> Option<&SearchPath> {
        self.search_path.as_ref()
    }

    /// The number of the volume that holds `ordinal` (`CSeqDBVolSet::FindVol`).
    pub fn volume_index(&self, ordinal: u32) -> Result<usize, SeqDbError> {
        self.find_volume(ordinal).map(|(volume, _)| volume)
    }

    /// Opens the one database at the path `name` (without extension), taken
    /// as it is: a space or a double quote in it is a character of the path
    /// (`alias::resolve_one_name`). Everything else is as in [`SeqDb::open`],
    /// including NCBI's rule that a volume named in an alias file is looked
    /// for in the current directory first.
    #[cfg(test)]
    pub(crate) fn open_one_name(name: &str, seqtype: SeqType) -> Result<Self, SeqDbError> {
        let tree = alias::resolve_one_name(name, seqtype, alias::LocalCopy::CurrentDir)?;
        Self::open_tree(tree, name, seqtype)
    }

    /// Resolves the one database at the path `name` as [`SeqDb::open_one_name`]
    /// does, without opening its volumes: the alias files are read and the
    /// index file of every volume they name is located. No volume file is
    /// read, so the cost does not grow with the number of sequences.
    pub(crate) fn locate_one_name(name: &str, seqtype: SeqType) -> Result<(), SeqDbError> {
        alias::resolve_one_name(name, seqtype, alias::LocalCopy::CurrentDir).map(|_| ())
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

        alias::check_sums(&tree.root, &|volume_name| {
            tree.volumes
                .iter()
                .position(|name| name == volume_name)
                .map(|at| {
                    (
                        u64::from(volumes[at].num_oids()),
                        volumes[at].index().volume_length,
                    )
                })
        })?;

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
        // volume list with the volumes. LOSAT requires the file, in the
        // volume's directory, and does not read it.
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

        let db = Self {
            seqtype,
            volumes,
            volume_starts,
            title,
            date,
            total_length,
            max_length,
            search_path: None,
        };
        Ok(db)
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
    /// Checks that the database holds the ordinal ids of a database built without
    /// `-parse_seqids`: the first header of every volume is decoded with the
    /// ordinal-id decoder (`deflines`), whose error for any other Seq-id is
    /// `SeqDbError::Unsupported`; the first error is returned. A volume whose
    /// header file cannot be opened (missing, or an `open` that fails) is not
    /// checked: NCBI opens a header file at its first read, and the error comes
    /// there. A search program that reads only ordinal ids calls this right
    /// after opening the database.
    pub fn check_ordinal_ids(&self) -> Result<(), SeqDbError> {
        for (at, volume) in self.volumes.iter().enumerate() {
            if volume.num_oids() > 0 {
                match self.deflines(self.volume_starts[at]) {
                    Ok(_) | Err(SeqDbError::NotFound(_) | SeqDbError::Open { .. }) => {}
                    Err(error) => return Err(error),
                }
            }
        }
        Ok(())
    }

    /// `check_ordinal_ids` of a program that reads parsed Seq-ids too: `Ok(true)` when the first
    /// header of some volume holds parsed Seq-ids, any other error as `check_ordinal_ids` gives
    /// it.
    pub fn has_parsed_seqids(&self) -> Result<bool, SeqDbError> {
        let mut parsed = false;
        for (at, volume) in self.volumes.iter().enumerate() {
            if volume.num_oids() > 0 {
                match self.deflines(self.volume_starts[at]) {
                    Ok(_) | Err(SeqDbError::NotFound(_) | SeqDbError::Open { .. }) => {}
                    Err(error) if error.is_parsed_seqids() => parsed = true,
                    Err(error) => return Err(error),
                }
            }
        }
        Ok(parsed)
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

    /// Number of sequences; ordinal ids run from 0 to this number minus one.
    pub fn num_seqs(&self) -> u32 {
        self.volume_starts[self.volumes.len()]
    }

    /// Total number of residues or bases.
    pub fn total_length(&self) -> u64 {
        self.total_length
    }

    /// Length of the longest sequence.
    pub fn max_length(&self) -> u32 {
        self.max_length
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvolset.hpp:210-226
    // ```c++
    //     CSeqDBVol * FindVol(int oid, int & vol_oid) const
    //     {
    //         ...
    //         int vol_idx(0);
    //         return const_cast<CSeqDBVol*>(FindVol(oid, vol_oid, vol_idx));
    //     }
    // ```
    /// The volume that holds `ordinal` and the sequence's number in it.
    fn find_volume(&self, ordinal: u32) -> Result<(usize, u32), SeqDbError> {
        if ordinal >= self.num_seqs() {
            return Err(SeqDbError::Range(format!(
                "ordinal id {ordinal} is outside the database ({} sequences)",
                self.num_seqs()
            )));
        }
        // The last volume whose first ordinal id is not above `ordinal`;
        // empty volumes share their start with the next volume.
        let volume =
            self.volume_starts[..self.volumes.len()].partition_point(|&s| s <= ordinal) - 1;
        Ok((volume, ordinal - self.volume_starts[volume]))
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
        let (volume, first) = self.find_volume(first_ordinal).ok()?;
        let end = self.volumes[volume]
            .block_end(first, max_bytes)
            .unwrap_or(first + 1);
        Some(BlockRange {
            volume,
            first,
            end,
            first_ordinal,
        })
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the block of the sequences `first_ordinal..end_ordinal`
    /// that the volume of `first_ordinal` holds (to `end_ordinal` or to the end of the volume),
    /// found without reading the index (`next_block` reads it to cut by bytes). `None` when the
    /// range is empty or `first_ordinal` is past the last sequence.
    pub fn x_block_range(&self, first_ordinal: u32, end_ordinal: u32) -> Option<BlockRange> {
        if first_ordinal >= end_ordinal {
            return None;
        }
        let (volume, first) = self.find_volume(first_ordinal).ok()?;
        let end = self.volumes[volume]
            .num_oids()
            .min(first.saturating_add(end_ordinal - first_ordinal));
        Some(BlockRange {
            volume,
            first,
            end,
            first_ordinal,
        })
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the numbers of residues of the proteins of
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

    /// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the length and the ambiguity data that
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
        let (volume, first) = self.find_volume(ordinal)?;
        self.volumes[volume].read_block(first, first + 1, ordinal)
    }

    /// EXPERIMENT (LOSAT_X_BPREDODBSUBJ): the number of residues of protein `ordinal`, read from the index
    /// alone (`read_sequence(ordinal)`'s `protein(0).len()`).
    pub fn protein_length(&self, ordinal: u32) -> Result<u32, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].protein_length(oid)
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
        let (volume, first) = self.find_volume(ordinal)?;
        self.volumes[volume].x_read_block_mapped(first, first + 1, ordinal)
    }

    /// Whether the ambiguity data of nucleotide `ordinal` decode
    /// (`Volume::check_ambiguity`: `SequenceBlock::ambiguity_runs`' error,
    /// read without the letters).
    pub fn check_ambiguity(&self, ordinal: u32) -> Result<(), SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].check_ambiguity(oid, ordinal)
    }

    /// The length in bases of nucleotide `ordinal`, read from the index and the last byte of
    /// its packed bases (no letters).
    pub fn nucleotide_length(&self, ordinal: u32) -> Result<usize, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].nucleotide_length(oid, ordinal)
    }

    /// The IUPAC letters `[from, to)` of nucleotide `ordinal`: the bases of the range and the
    /// ambiguity runs that overlap it, read without the rest of the sequence.
    pub fn nucleotide_iupac_range(
        &self,
        ordinal: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].nucleotide_iupac_range(oid, ordinal, from, to)
    }

    /// The length and the ambiguity runs of nucleotide `ordinal`, read once for any number of
    /// [`SeqDb::nucleotide_iupac_range_with`] calls on it.
    pub fn nucleotide_ambiguity(&self, ordinal: u32) -> Result<NucleotideAmbiguity, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].nucleotide_ambiguity(oid, ordinal)
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
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].nucleotide_iupac_range_with(oid, ordinal, from, to, ambiguity)
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
    /// ordinal id of the volume, as NCBI reports it.
    pub fn deflines(&self, ordinal: u32) -> Result<Vec<OrdinalDefline>, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.ordinal_deflines(volume, ordinal, self.volumes[volume].deflines(oid))
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
        match self.deflines(ordinal) {
            Err(error) if error.is_parsed_seqids() => {
                self.defline_set(ordinal).map(SequenceHeader::Parsed)
            }
            read => read.map(SequenceHeader::Ordinal),
        }
    }

    /// The definition lines of `ordinal` with every member and every Seq-id
    /// (`seqid::ber::decode_full_set`), as NCBI's `GetHdr` gives them: the
    /// tag of a `general` id `BL_ORD_ID` is its database ordinal id (the
    /// stored number plus the first ordinal id of the volume).
    pub fn defline_set(&self, ordinal: u32) -> Result<DeflineSet, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        let mut set = self.volumes[volume].defline_set(oid)?;
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
    /// apart from `deflines`, only for the XML2 and JSON formats.
    pub fn defline_taxids(&self, ordinal: u32) -> Result<Vec<i32>, SeqDbError> {
        let (volume, oid) = self.find_volume(ordinal)?;
        self.volumes[volume].defline_taxids(oid)
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
    /// header or ambiguity data reads one at a time
    /// ([`SeqDb::visit_deflines`], [`SeqDb::visit_ambiguity`]).
    pub fn header_ranges(&self, max_seqs: u32) -> Vec<BlockRange> {
        let max_seqs = max_seqs.max(1);
        let mut ranges = Vec::new();
        for (volume, vol) in self.volumes.iter().enumerate() {
            let mut first = 0;
            while first < vol.num_oids() {
                let end = first + max_seqs.min(vol.num_oids() - first);
                ranges.push(BlockRange {
                    volume,
                    first,
                    end,
                    first_ordinal: self.volume_starts[volume] + first,
                });
                first = end;
            }
        }
        ranges
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
                return BlockRange {
                    volume,
                    first,
                    end,
                    first_ordinal: self.volume_starts[volume] + first,
                };
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
            return range.is_empty().then_some(true).unwrap_or_else(|| {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            });
        };
        vol.visit_deflines_as(
            range.first,
            range.end,
            header_piece_bytes(),
            want_taxid,
            &mut |oid, deflines| {
                let ordinal = range.first_ordinal + (oid - range.first);
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
            return range.is_empty().then_some(true).unwrap_or_else(|| {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            });
        };
        vol.visit_defline_taxids(
            range.first,
            range.end,
            header_piece_bytes(),
            &mut |oid, taxids| visit(range.first_ordinal + (oid - range.first), taxids),
        )
    }

    /// EXPERIMENT (LOSAT_X_DBIDXREAD): calls `visit(ordinal, self.check_ambiguity(ordinal))`
    /// for the sequences of `range` in order until `visit` returns false; returns whether every
    /// call returned true. The offsets of the range are read at once
    /// (`Volume::visit_ambiguity`).
    pub fn visit_ambiguity(
        &self,
        range: BlockRange,
        mut visit: impl FnMut(u32, Result<(), SeqDbError>) -> bool,
    ) -> bool {
        let Some(vol) = self.volumes.get(range.volume) else {
            return range.is_empty().then_some(true).unwrap_or_else(|| {
                visit(range.first_ordinal, Err(volume_range_error(range, self)))
            });
        };
        vol.visit_ambiguity(
            range.first,
            range.end,
            range.first_ordinal,
            &mut |oid, result| {
                let ordinal = range.first_ordinal + (oid - range.first);
                visit(ordinal, result)
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
    /// A record whose only Seq-ids are `gnl|BL_ORD_ID|n` (`SeqDb::deflines`).
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
