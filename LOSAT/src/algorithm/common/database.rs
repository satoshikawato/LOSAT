//! What the search programs share in a search of an NCBI BLAST database (`-db NAME`, read
//! through `crate::seqdb`): the opening of the database by its name, the reader's errors
//! as the command line reports them, the `Database:` lines of the pairwise report, and the
//! rules for the stored titles of the reported sequences. BLASTP (stage DB-3) and BLASTN
//! (stage DB-4) use them with their own name in LOSAT's refusals.
//!
//! The reader (`crate::seqdb`) knows nothing of the command line; this module maps its
//! errors to NCBI's application errors (`crate::cli::NativeError`) and to LOSAT's
//! refusals.

use anyhow::{bail, Result};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
// ```c++
//     static const string GenerateSearchPath() {
//         string splitter;
//         string path;
// #if defined(NCBI_OS_UNIX)
//         splitter = ":";
// #else
//         splitter = ";";
// #endif
//         // Local directory first;
//         path  = CDirEntry::NormalizePath(CDir::GetCwd(),eFollowLinks);
//         path += splitter;
//         // Then, BLASTDB;
//         CNcbiEnvironment env;
//         path += CDirEntry::NormalizePath(env.Get("BLASTDB"),eFollowLinks);
//         path += splitter;
//         // Finally, the config file.
//         CNcbiApplicationAPI* app = CNcbiApplicationAPI::Instance();
//         if (app) {
//             const CNcbiRegistry& registry = app->GetConfig();
//             if (registry.HasEntry("BLAST", "BLASTDB")) {
//                 path += CDirEntry::NormalizePath(registry.Get("BLAST", "BLASTDB"),eFollowLinks);
//                 path += splitter;
//             }
//         }
//         return path;
//     }
// ```
/// Opens the database `name` (the `-db` value as NCBI keeps it, `seqdb::db_argument_name`)
/// of the sequence type `seqtype` for a search: the one function through which a search
/// program opens a database. The name is looked for as NCBI looks for it: in the current
/// directory, then in the directories of the `BLASTDB` variable, then in those of the
/// registry entry `[BLAST] BLASTDB` (`settings.blastdb`, read from
/// `NCBI_CONFIG__BLAST__BLASTDB`, the program's `.ini` file and `.ncbirc`).
pub(crate) fn open_search_database(
    name: &str,
    seqtype: SeqType,
    settings: &ApplicationSettings,
) -> Result<SeqDb, DatabaseOpenError> {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => return Err(DatabaseOpenError::Registry(reason.clone())),
    };
    let path = SearchPath::from_environment(registry).map_err(DatabaseOpenError::Reader)?;
    SeqDb::open_in(name, seqtype, &path).map_err(DatabaseOpenError::Reader)
}

/// Why `open_search_database` could not open a database.
#[derive(Debug)]
pub(crate) enum DatabaseOpenError {
    /// The registry entry `[BLAST] BLASTDB` has a value that LOSAT does not read as NCBI
    /// does (`ApplicationSettings::blastdb`).
    Registry(String),
    /// An error of the reader (NCBI's own, or LOSAT's refusal: `database_error`).
    Reader(SeqDbError),
}

/// Opens the database of `-db name` (`open_search_database`) with the command line's
/// errors (`database_error`; LOSAT's refusals name `program`).
pub(crate) fn open_database(
    name: &str,
    seqtype: SeqType,
    settings: &ApplicationSettings,
    program: &str,
) -> Result<SeqDb> {
    open_search_database(name, seqtype, settings).map_err(|error| match error {
        DatabaseOpenError::Registry(reason) => refusal(&reason, program),
        DatabaseOpenError::Reader(error) => database_error(error, program),
    })
}

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
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:153
// ```c++
// #define BLAST_DATABASE_ERROR                    2
// ```
/// A reader error outside the preliminary search (when the database is opened, and when
/// the sequences and titles of the hit lists are read) as the command line reports it:
/// NCBI's own errors (`SeqDbError::ncbi_message`: a database that is not found, a file
/// that is no database file or cannot be opened) as `BLAST Database error: <text>` with
/// exit status 2; everything else is LOSAT's refusal, in LOSAT's wording (`program` is
/// the program that refuses).
pub(crate) fn database_error(error: SeqDbError, program: &str) -> anyhow::Error {
    match error.ncbi_message() {
        Some(text) => NativeError {
            exit: 2,
            message: format!("BLAST Database error: {text}\n"),
        }
        .into(),
        None => refusal(&error.to_string(), program),
    }
}

/// LOSAT's refusal of a database that `program` does not search as NCBI does.
pub(crate) fn refusal(reason: &str, program: &str) -> anyhow::Error {
    anyhow::anyhow!("{reason}; this database is not supported by LOSAT's {program}")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251
// ```c++
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
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
/// A reader error while the preliminary search reads the sequences: NCBI's search stage
/// replaces the reader's text (`SeqDbError::ncbi_search_message`), and the application
/// reports `BLAST engine error: <text>` with exit status 3. The error is kept apart from
/// the others (`PreliminaryReadError`), because NCBI has written the outfmt 0 prolog and
/// read the query batch before it.
pub(crate) fn search_read_error(error: SeqDbError, program: &str) -> anyhow::Error {
    match error.ncbi_search_message() {
        Some(text) => PreliminaryReadError(text).into(),
        None => refusal(&error.to_string(), program),
    }
}

/// NCBI's error of a database file that cannot be read during the preliminary search
/// (`search_read_error`), with the text after `BLAST engine error: `.
#[derive(Debug)]
pub(crate) struct PreliminaryReadError(pub(crate) &'static str);

impl std::fmt::Display for PreliminaryReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BLAST engine error: {}", self.0)
    }
}

impl std::error::Error for PreliminaryReadError {}

/// A reader error of NCBI's own while the definition lines of the final hit lists are
/// read (`read_title`): `query` is the query whose hit list reached it. NCBI fails that
/// query's batch, after the reports of the batches before.
#[derive(Debug)]
pub(crate) struct HeaderReadError {
    pub(crate) query: usize,
    pub(crate) error: anyhow::Error,
}

impl std::fmt::Display for HeaderReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for HeaderReadError {}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
// ```c++
//     info.name = seqdb->GetDBNameList();
//     info.definition = seqdb->GetTitle();
//     if (info.definition.empty())
//         info.definition = info.name;
//     info.date = seqdb->GetDate();
//     info.total_length = seqdb->GetTotalLength();
//     info.number_seqs = seqdb->GetNumSeqs();
// ```
/// What the `Database:` lines of the pairwise report show of the database `db` opened by
/// the name `name`: the title (the name as given for an empty title), the date and the
/// real totals, also with `-dbsize`.
pub(crate) fn pairwise_database(db: &SeqDb, name: &str) -> PairwiseDatabase {
    let title = db.title();
    PairwiseDatabase {
        definition: if title.is_empty() {
            name.as_bytes().to_vec()
        } else {
            title.to_vec()
        },
        date: db.date().to_string(),
        num_sequences: u64::from(db.num_seqs()),
        total_letters: db.total_length(),
    }
}

/// The stored title of a database sequence, as the reports read it (`read_title`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct StoredTitle {
    /// The title of the sequence's definition line, with a byte outside ' ' ..= '~' read
    /// as '#'.
    pub(crate) title: Vec<u8>,
    /// Whether the stored title held a byte outside ' ' ..= '~' (`check_reported_title`).
    pub(crate) replaced_bytes: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310
// ```c++
// CTempString
// CSeqDBVol::x_GetHdrAsn1Binary(int oid) const
// {
//     TIndx hdr_start = 0;
//     TIndx hdr_end   = 0;
//
//     //m_Atlas.Lock(locked);
//
//     if (!m_HdrFileOpened) x_OpenHdrFile();
// ```
// NCBI reference for the titles:
// c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:970-1017 (`s_GetBioseqTitle`: the
// title of the first definition line, the others appended as ` >id title`).
/// Reads the definition line of the database sequence `ordinal` for a report, where the
/// hit list of query `query` reaches it. A reader error of NCBI's own is a
/// `HeaderReadError` that names the query (NCBI fails that query's batch); LOSAT's
/// refusals (`program` names the program) are errors of their own.
pub(crate) fn read_title(
    db: &SeqDb,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<StoredTitle> {
    let mut deflines = db.deflines(ordinal).map_err(|error| {
        if error.ncbi_message().is_some() {
            anyhow::Error::new(HeaderReadError {
                query,
                error: database_error(error, program),
            })
        } else {
            database_error(error, program)
        }
    })?;
    // makeblastdb without -parse_seqids writes one definition line per sequence; NCBI
    // shows the others of a sequence as ` >id title` after the first.
    if deflines.len() != 1 {
        bail!(
            "database sequence {ordinal} has {} definition lines; a sequence with more or fewer than one is not supported by LOSAT's {program}",
            deflines.len()
        );
    }
    let defline = deflines.remove(0);
    Ok(StoredTitle {
        title: defline.title,
        replaced_bytes: defline.replaced_title_bytes > 0,
    })
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1865-1886
// ```c++
//         string message;
//         if (io != NULL) {
//             message += io->GetStackTrace() + "\n";
//         }
//         message += "Bad char [0x" +
//                    NStr::NumericToString((unsigned char)c,0,16)+
//                    "] in string";
//         if (io != NULL) {
//             message += " at " + io->GetPosition();
//         }
//         if (!str.empty()) {
//             message += "\n" + str;
//         }
//         switch (fix_method) {
//         case eFNP_ReplaceAndWarn:
// #if 0
//             CNcbiDiag(eDiag_Error, eDPF_Default)
//                 << ErrCode(NCBI_ERRCODE_X, 7) << message << Endm;
// #else
//             ERR_POST_X(7, message);
// #endif
//             break;
// ```
// NCBI BLAST+ 2.17.0 `blastp -db` on a version 4 database whose header file holds
// the bytes 0xE9 and 0x7F in two titles (campaign dbsearch-20261012,
// `$BUILD_ROOT/dbsearch/db3/exp_badbyte`): the report shows `#` for the byte, and
// standard error gets, for each reported subject and bad byte,
//   Error: (803.7) [blastp] Blast-def-line-set.E.title
//   Bad char [0xE9] in string at byte 21
//   <the stored title>
// as often as the formatter reads the definition line: twice per subject for the
// default outfmt 6, three times for `6 sseqid stitle`, five times for outfmt 0. The
// number of readings is the formatter's, not ported; makeblastdb and `LOSAT makedb`
// store `#` for such a byte, so only a header file written otherwise holds one.
/// LOSAT's rejection of the reported database sequence `ordinal` (one with rows in the
/// final hit lists) whose stored title held a byte outside ' ' ..= '~' (`replaced_bytes`,
/// `StoredTitle`): NCBI writes an error message for the byte each time its formatter
/// reads the definition line.
pub(crate) fn check_reported_title(
    replaced_bytes: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if replaced_bytes {
        bail!(
            "database sequence {ordinal} has a stored title with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line); this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of the database
/// sequence `ordinal` with the stored title `title` (`report::defline::database_subject_id`;
/// `protein` for a protein database). A title whose id LOSAT cannot write as NCBI is
/// LOSAT's refusal (`program` names the program).
pub(crate) fn tabular_subject_id(
    title: &[u8],
    ordinal: u32,
    protein: bool,
    program: &str,
) -> Result<Vec<u8>> {
    crate::report::defline::database_subject_id(title, protein).map_err(|problem| {
        anyhow::anyhow!(
            "database sequence {ordinal} has a title {problem}, which is not supported by LOSAT's {program}"
        )
    })
}

/// LOSAT's rejection of the outfmt 0 title `title` of the protein database sequence
/// `ordinal` that it does not write as NCBI: a title of punctuation that NCBI's
/// `x_CleanAndCompress` reads past the end of when it makes the title of the description
/// table (NCBI BLAST+ 2.17.0 crashes there), as for the records of a subject set
/// (`check_shown_subject_titles`). Approved exception 2 of PD-LOSAT-NCBI-DEFECTS covers
/// the nucleotide subjects of BLASTN, so BLASTN does not use it.
pub(crate) fn check_shown_protein_title(title: &[u8], ordinal: u32, program: &str) -> Result<()> {
    let reads_past_end = crate::report::defline::generate_defline(title, true, true)
        .is_ok_and(|defline| defline.reads_past_end);
    if reads_past_end {
        bail!(
            "database sequence {ordinal} has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0); this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}
