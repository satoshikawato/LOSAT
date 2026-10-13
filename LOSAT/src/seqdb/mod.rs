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
//! file (`.phr` / `.nhr`: title and `gnl|BL_ORD_ID|n`), and alias files
//! (`.pal` / `.nal`: `TITLE` and `DBLIST`).
//!
//! Not read: the LMDB file and the id and taxonomy files of version 5, the
//! ISAM files of `-parse_seqids`, column files (mask data), group alias files,
//! and alias keys that filter sequences or replace the totals. A database
//! that needs one of these to give NCBI's result is rejected with
//! [`SeqDbError::Unsupported`]; the version 5 files that only serve lookups
//! by id or taxonomy are left unread. The LMDB file of a version 5 database
//! must exist, as for NCBI, but its contents are not checked.
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
pub mod volume;

#[cfg(test)]
mod tests;

use std::path::{Path, PathBuf};

pub use alias::SearchPath;
pub use header::Defline;
pub use volume::{AmbiguityRun, SequenceBlock, Volume, BLASTNA_SENTINEL};

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
    /// A database that stores parsed Seq-ids is rejected here: the first
    /// header of every volume is decoded once. A volume whose header file
    /// cannot be opened (missing, or an `open` that fails) is not checked:
    /// NCBI opens a header file at its first read, and the error comes there.
    pub fn open_in(name: &str, seqtype: SeqType, path: &SearchPath) -> Result<Self, SeqDbError> {
        Self::open_tree(alias::resolve_in(name, seqtype, path)?, name, seqtype)
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
        };
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
        for (at, volume) in db.volumes.iter().enumerate() {
            if volume.num_oids() > 0 {
                match db.deflines(db.volume_starts[at]) {
                    Ok(_) | Err(SeqDbError::NotFound(_) | SeqDbError::Open { .. }) => {}
                    Err(error) => return Err(error),
                }
            }
        }
        Ok(db)
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
    pub fn next_block(&self, first_ordinal: u32, max_bytes: u64) -> Option<BlockRange> {
        let (volume, first) = self.find_volume(first_ordinal).ok()?;
        let vol = &self.volumes[volume];
        // The ends `first + 1 ..= num_oids`; the offsets ascend, so the ends
        // whose bytes fit form a prefix.
        let candidates = first + 1..vol.num_oids() + 1;
        let fitting =
            partition_point_u32(candidates, |end| vol.block_bytes(first, end) <= max_bytes);
        let end = first + (fitting.max(1) as u32);
        Some(BlockRange {
            volume,
            first,
            end,
            first_ordinal,
        })
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

    /// Reads the one sequence `ordinal` (for traceback and reports).
    pub fn read_sequence(&self, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
        let (volume, first) = self.find_volume(ordinal)?;
        self.volumes[volume].read_block(first, first + 1, ordinal)
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
        let volume_start = self.volume_starts[volume];
        let vol = &self.volumes[volume];
        vol.deflines(oid)?
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
                })
            })
            .collect()
    }
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
}
