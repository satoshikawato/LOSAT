//! What the search programs share in a search of an NCBI BLAST database (`-db NAME`, read
//! through `crate::seqdb`): the opening of the database by its name, the reader's errors
//! as the command line reports them, the `Database:` lines of the pairwise report, and the
//! rules for the stored titles of the reported sequences. BLASTP (stage DB-3), BLASTN
//! (stage DB-4) and TBLASTN (stage DB-5) use them with their own name in LOSAT's refusals;
//! TBLASTN also reads a nucleotide database whole into packed subjects
//! (`read_nucleotide_stores`).
//!
//! The reader (`crate::seqdb`) knows nothing of the command line; this module maps its
//! errors to NCBI's application errors (`crate::cli::NativeError`) and to LOSAT's
//! refusals.

use anyhow::{bail, Result};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

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

/// LOSAT's refusal of a database with parsed Seq-ids (`SeqDb::check_ordinal_ids`), which
/// `program` does not search: psiblast calls it right after it opens its database
/// (`open_database`), so the refusal comes where the reader gave it before; the five search
/// programs read such a database (`check_seqids`, line DBID ID2a and ID2b).
pub(crate) fn refuse_parsed_seqids(db: &SeqDb, program: &str) -> Result<()> {
    db.check_ordinal_ids()
        .map_err(|error| database_error(error, program))
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

/// Bytes of the sequence file that one read of `read_nucleotide_stores` holds at most
/// (`SeqDb::next_block`; a longer sequence is a read of its own). No NCBI counterpart (NCBI
/// maps the volume files); the size does not change what is read.
const STORE_BLOCK_BYTES: u64 = 1 << 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1405
// ```c
//     memset((void*) &seq_arg, 0, sizeof(seq_arg));
//
//     /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
//        sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
//     seq_arg.encoding = eBlastEncodingProtein;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:771-775
// ```c
//         } else {
//             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                                      backup.full_range.right,
//                                      subject->gen_code_string, &translation_buffer,
//                                      &frame_offsets, NULL);
// ```
/// Every sequence of the nucleotide database `db`, in ordinal order, as a packed subject
/// (`PackedNucleotide::from_stored`): the stored 2-bit data, which a translated search's
/// preliminary stage translates as NCBI does (makeblastdb's bases for the ambiguous ones),
/// and the stored ambiguity runs, from which the traceback and the reports read the real
/// codes. The file is read in blocks of at most `STORE_BLOCK_BYTES`, each dropped once its
/// sequences are packed, so the database is held packed only. NCBI reads every sequence's
/// 2-bit data in the first query batch's preliminary search, so a block that cannot be read
/// is that search's error (`search_read_error`).
pub(crate) fn read_nucleotide_stores(db: &SeqDb, program: &str) -> Result<Vec<PackedNucleotide>> {
    let mut stores = Vec::with_capacity(db.num_seqs() as usize);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let (length, runs) = stored_parts(&block, i, range.first_ordinal + i as u32, program)?;
            stores.push(PackedNucleotide::from_stored(
                block.nucleotide_packed(i),
                length,
                &runs,
            ));
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(stores)
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the checks of `read_nucleotide_stores` without the copy. The
/// database is read in the same blocks, and every sequence's length and ambiguity runs are read
/// with the same errors in the same order; only the lengths are kept, by ordinal id. A search
/// that reads its subjects where it searches them (`x_read_nucleotide_store`) fails the query
/// batch that the copy would fail, with the same text. A block of one sequence (one longer than
/// a block, or the last of its volume) is checked without its 2-bit data
/// (`SeqDb::x_nucleotide_check`): a human chromosome is not read here only to be dropped.
pub(crate) fn x_check_nucleotide_stores(db: &SeqDb, program: &str) -> Result<Vec<usize>> {
    let mut lengths = Vec::with_capacity(db.num_seqs() as usize);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, STORE_BLOCK_BYTES) {
        if range.len() == 1 {
            let (length, ambiguity) = db
                .x_nucleotide_check(range)
                .map_err(|error| search_read_error(error, program))?;
            // The error of `stored_parts`.
            ambiguity.ambiguity_runs(0).map_err(|error| {
                refusal(
                    &format!(
                        "the ambiguity data of database sequence {} cannot be read ({error})",
                        range.first_ordinal
                    ),
                    program,
                )
            })?;
            lengths.push(length);
            next_ordinal = range.end_ordinal();
            continue;
        }
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let (length, _) = stored_parts(&block, i, range.first_ordinal + i as u32, program)?;
            lengths.push(length);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(lengths)
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): sequence `ordinal` of `db` as `read_nucleotide_stores` packs
/// it (`PackedNucleotide::from_stored` of the stored 2-bit data, length and ambiguity runs), read
/// with one read of the sequence file. The errors are those of `read_nucleotide_stores` for the
/// same sequence. With `mapped`, the bytes are read into a buffer that glibc maps from 1 MiB on,
/// taken from `utils::x_mapped::X_BYTES` when it holds one (`SeqDb::x_read_sequence_mapped`), and
/// the store is that buffer cut to the 2-bit data (`PackedNucleotide::x_from_stored_buffer`): a
/// long sequence read for a while and dropped gives its pages back (a store is not given to the
/// pool, which forgets the buffer). Without it they are read as `SeqDb::read_sequence` reads
/// them (the heap of the thread), where a dropped store's pages serve the next read.
pub(crate) fn x_read_nucleotide_store(
    db: &SeqDb,
    ordinal: u32,
    program: &str,
    mapped: bool,
) -> Result<PackedNucleotide> {
    let block = if mapped {
        db.x_read_sequence_mapped(ordinal)
    } else {
        db.read_sequence(ordinal)
    }
    .map_err(|error| search_read_error(error, program))?;
    let (length, runs) = stored_parts(&block, 0, ordinal, program)?;
    let packed = block.x_into_nucleotide_packed();
    if crate::utils::x_mapped::x_is_mapped(packed.capacity()) {
        crate::utils::x_mapped::X_BYTES.x_disown(&packed);
    }
    Ok(PackedNucleotide::x_from_stored_buffer(
        packed, length, &runs,
    ))
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): sequence `i` of `block` (ordinal id `ordinal`), a block of
/// `SeqDb::read_block`, as `read_nucleotide_stores` packs it, with its errors.
pub(crate) fn x_nucleotide_store_in(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<PackedNucleotide> {
    let (length, runs) = stored_parts(block, i, ordinal, program)?;
    Ok(PackedNucleotide::from_stored(
        block.nucleotide_packed(i),
        length,
        &runs,
    ))
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the length of sequence `i` of `block` (ordinal id
/// `ordinal`), with the errors of `x_nucleotide_store_in`, without the store.
pub(crate) fn x_nucleotide_length_in(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<usize> {
    stored_parts(block, i, ordinal, program).map(|(length, _)| length)
}

/// The length and the ambiguity runs of sequence `i` of `block` (ordinal id `ordinal`), with the
/// errors of `read_nucleotide_stores` (which reads every sequence through it). NCBI's preliminary
/// search reads the 2-bit data only; it reads the ambiguity data in the traceback, of the
/// reported sequences. LOSAT reads them for every sequence, so data it cannot read are LOSAT's
/// refusal, not an NCBI error.
fn stored_parts(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<(usize, Vec<crate::seqdb::AmbiguityRun>)> {
    let length = block
        .nucleotide_len(i)
        .map_err(|error| search_read_error(error, program))?;
    let runs = block.ambiguity_runs(i).map_err(|error| {
        refusal(
            &format!("the ambiguity data of database sequence {ordinal} cannot be read ({error})"),
            program,
        )
    })?;
    Ok((length, runs))
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
    stored_title(db.deflines(ordinal), ordinal, query, program)
}

/// `read_title` of the sequence `ordinal` whose definition lines `deflines` were read
/// (`SeqDb::deflines`, or `SeqDb::visit_deflines` with the same results).
pub(crate) fn stored_title(
    deflines: std::result::Result<Vec<crate::seqdb::OrdinalDefline>, crate::seqdb::SeqDbError>,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<StoredTitle> {
    let mut deflines = deflines.map_err(|error| header_read_error(error, query, program))?;
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

/// A reader error of the header of a reported sequence where the hit list of query `query`
/// reaches it: NCBI's own is a `HeaderReadError` that names the query; LOSAT's refusals are
/// errors of their own.
fn header_read_error(error: SeqDbError, query: usize, program: &str) -> anyhow::Error {
    if error.ncbi_message().is_some() {
        anyhow::Error::new(HeaderReadError {
            query,
            error: database_error(error, program),
        })
    } else {
        database_error(error, program)
    }
}

/// The header of a reported database sequence (`read_header`).
#[derive(Debug, Clone)]
pub(crate) enum ReportedHeader {
    /// A sequence of a database without parsed Seq-ids: its one title (`read_title`).
    Stored(StoredTitle),
    /// A sequence with parsed Seq-ids (`report::db_subject`).
    Parsed(std::sync::Arc<crate::report::db_subject::ParsedSubject>),
}

#[cfg(test)]
impl ReportedHeader {
    /// The stored title of a sequence without parsed Seq-ids.
    pub(crate) fn stored(&self) -> &StoredTitle {
        match self {
            ReportedHeader::Stored(stored) => stored,
            ReportedHeader::Parsed(_) => panic!("a sequence with parsed Seq-ids"),
        }
    }
}

/// `read_title` of a program that also reads databases with parsed Seq-ids: the header of the
/// sequence `ordinal`, read as `read_title` reads it, or with every Seq-id where the record has
/// parsed ones (`SeqDb::read_header`). A parsed subject that LOSAT does not write as NCBI
/// (`ParsedSubject::new`) is LOSAT's refusal.
pub(crate) fn read_header(
    db: &SeqDb,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<ReportedHeader> {
    match db.read_header(ordinal) {
        Ok(crate::seqdb::SequenceHeader::Parsed(set)) => {
            crate::report::db_subject::ParsedSubject::new(set, db.is_version5())
                .map(|subject| ReportedHeader::Parsed(std::sync::Arc::new(subject)))
                .map_err(|reason| parsed_refusal(ordinal, &reason, program))
        }
        Ok(crate::seqdb::SequenceHeader::Ordinal(deflines)) => {
            stored_title(Ok(deflines), ordinal, query, program).map(ReportedHeader::Stored)
        }
        Err(error) => Err(header_read_error(error, query, program)),
    }
}

/// LOSAT's refusal of a reported database sequence with parsed Seq-ids (`parsed_refusal`,
/// `refuse_joined_hits`, `refuse_parsed_taxonomy`). BLASTX finds one only when it reports a
/// batch, after the start of its document: a run that holds its XML, XML2 or JSON document
/// drops it on this error (`ReportedSubjects::for_batch`).
#[derive(Debug)]
pub(crate) struct ParsedRefusal(String);

impl std::fmt::Display for ParsedRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ParsedRefusal {}

/// LOSAT's refusal of the database sequence `ordinal` with parsed Seq-ids for `reason` (what the
/// sequence has, `ParsedSubject`).
pub(crate) fn parsed_refusal(ordinal: u32, reason: &str, program: &str) -> anyhow::Error {
    anyhow::Error::new(ParsedRefusal(format!(
        "database sequence {ordinal} {reason}; this is not supported by LOSAT's {program}"
    )))
}

/// `refuse_parsed_seqids` of a program that searches databases with parsed Seq-ids too: the
/// first header of each volume is read right after the database is opened, as before, and only
/// a header with parsed Seq-ids is let through (`SeqDb::has_parsed_seqids`). Returns whether the
/// database has parsed Seq-ids; such a database is LOSAT's refusal where a variable changes
/// NCBI's choice of its Seq-ids (`ApplicationSettings::seq_id_rank_variable`).
pub(crate) fn check_seqids(
    db: &SeqDb,
    settings: &ApplicationSettings,
    program: &str,
) -> Result<bool> {
    let parsed = db
        .has_parsed_seqids()
        .map_err(|error| database_error(error, program))?;
    if let (true, Some(name)) = (parsed, settings.seq_id_rank_variable) {
        bail!(
            "the environment variable {name}, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's {program}"
        );
    }
    Ok(parsed)
}

/// LOSAT's refusal of two parsed subjects that follow each other in the hit list of query
/// `query` and whose Seq-align ids match (`report::db_subject::joined_in_reports`): NCBI's XML,
/// XML2, JSON and outfmt 0 reports show their alignments as one hit, which LOSAT does not.
/// `subjects` gives the ordinal id and the parsed subject of each HSP list of the hit list, in
/// its order (`None` for a sequence without parsed Seq-ids).
pub(crate) fn refuse_joined_hits<'a>(
    query: usize,
    subjects: impl IntoIterator<Item = Option<(u32, &'a crate::report::db_subject::ParsedSubject)>>,
    program: &str,
) -> Result<()> {
    let mut previous: Option<(u32, &crate::report::db_subject::ParsedSubject)> = None;
    for subject in subjects {
        if let (Some((before, a)), Some((ordinal, b))) = (previous, subject) {
            if crate::report::db_subject::joined_in_reports(a, b) {
                return Err(anyhow::Error::new(ParsedRefusal(format!(
                    "database sequences {before} and {ordinal} follow each other in the hits of query {} and have Seq-align ids that NCBI BLAST+ takes for one sequence (its XML and outfmt 0 reports show their alignments as one hit); this is not supported by LOSAT's {program}",
                    query + 1
                ))));
            }
        }
        previous = subject;
    }
    Ok(())
}

/// LOSAT's refusal of a parsed subject whose outfmt 0 title NCBI makes from the source
/// (`ParsedSubject::title_from_taxonomy`): with `taxdb.bti` on the search path NCBI names the
/// organism of the taxids in it, which LOSAT does not do.
pub(crate) fn refuse_parsed_taxonomy(ordinal: u32, program: &str) -> anyhow::Error {
    anyhow::Error::new(ParsedRefusal(format!(
        "database sequence {ordinal} has an empty title that NCBI BLAST+'s outfmt 0 makes from its taxids and the taxonomy database (taxdb.bti is on the search path); this is not supported by LOSAT's {program}"
    )))
}

/// LOSAT's refusals of the outfmt 0 description and alignment heading of the subject with
/// parsed Seq-ids `subject` (the database sequence `ordinal`; `protein` for a protein database)
/// that it does not write as NCBI: a protein title of punctuation that NCBI's `x_CleanAndCompress`
/// reads past the end of (NCBI BLAST+ 2.17.0 crashes or does not end; approved exception 2 of
/// PD-LOSAT-NCBI-DEFECTS has BLASTN, TBLASTN and TBLASTX stop such a nucleotide title at the end
/// of the string, as for a database without parsed Seq-ids), a description id longer
/// than NCBI's cut of the line leaves room for (NCBI then pads the line without end), a
/// generated heading title whose end NCBI rewrites from the taxids
/// (`ParsedSubject::title_suffix_adjusted`), and a title that NCBI makes from the sequence's
/// taxids and the taxonomy database (`taxdb_found`, `ParsedSubject::title_from_taxonomy`).
pub(crate) fn check_shown_parsed_subject(
    subject: &crate::report::db_subject::ParsedSubject,
    protein: bool,
    taxdb_found: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    let (start, title) = subject.description(protein);
    let generated = subject.heading_generates_title();
    let heading_reads_past_end = generated
        && subject
            .bioseq_title(protein, false)
            .is_ok_and(|defline| defline.reads_past_end);
    if protein && (title.reads_past_end || heading_reads_past_end) {
        return Err(parsed_refusal(
            ordinal,
            "has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0)",
            program,
        ));
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:913-916
    // ```c++
    //         line_component = (line_component.empty() ? "" : "  ") + sdl->defline;
    //         string actual_line_component;
    //         if(line_component.size()+line_length > m_LineLen){
    //             actual_line_component = line_component.substr(0, m_LineLen -
    //                                                           line_length - 3);
    // ```
    // The unsigned length wraps for an id longer than 65 bytes: NCBI 2.17.0 then pads the
    // line with spaces without end (probe of the ID2a review, `tests/db_search/ids/refuse.faa`).
    if start.len() > 65 {
        return Err(parsed_refusal(
            ordinal,
            "has an outfmt 0 description id longer than 65 bytes (NCBI BLAST+ 2.17.0 pads the description line without end)",
            program,
        ));
    }
    if generated && subject.title_suffix_adjusted(protein) {
        return Err(parsed_refusal(
            ordinal,
            "has a definition line without a title and a title ending in ']' that NCBI BLAST+'s outfmt 0 heading rewrites from the taxids (x_AdjustProteinTitleSuffix)",
            program,
        ));
    }
    if taxdb_found && subject.title_from_taxonomy() {
        return Err(refuse_parsed_taxonomy(ordinal, program));
    }
    Ok(())
}

/// What the reports of a run show of their reported database sequences, for LOSAT's refusals of
/// those with parsed Seq-ids (`check_reported_parsed`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ParsedShown {
    /// The tabular formats (`sseqid`, `ParsedSubject::tabular_ids`).
    pub(crate) tabular: bool,
    /// Outfmt 0 (`check_shown_parsed_subject`).
    pub(crate) pairwise: bool,
    /// The formats that show the alignments of following hits whose Seq-align ids match as one
    /// hit: outfmt 0, 5 and 13 to 16 (`refuse_joined_hits`).
    pub(crate) joined: bool,
    /// SAM (outfmt 17), which names a subject by an id of its Bioseq (`ParsedSubject::sam_id`).
    pub(crate) sam: bool,
}

impl ParsedShown {
    /// Whether a report reads the subject's Bioseq, which NCBI's formatters find by its
    /// Seq-align id (`ParsedSubject::fetched`); the Seq-annots and the archive write the
    /// Seq-align id alone.
    pub(crate) fn reads_bioseq(self) -> bool {
        self.tabular || self.pairwise || self.joined || self.sam
    }
}

/// LOSAT's refusals of the reported database sequence `ordinal` with parsed Seq-ids `subject` in
/// the reports `shown` (`protein` for a protein database; `taxdb_found`: `SeqDb::taxdb_found`),
/// as BLASTX's (`blastx::database::ReportedSubjects::for_batch`): a sequence that NCBI's
/// formatters do not find by its Seq-align id, an outfmt 0 description or heading that LOSAT
/// does not write as NCBI (`check_shown_parsed_subject`), and a tabular id that it does not
/// write as NCBI (`ParsedSubject::tabular_ids`).
pub(crate) fn check_reported_parsed(
    subject: &crate::report::db_subject::ParsedSubject,
    shown: ParsedShown,
    protein: bool,
    taxdb_found: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if shown.reads_bioseq() {
        subject
            .fetched()
            .map_err(|reason| parsed_refusal(ordinal, &reason, program))?;
    }
    if shown.pairwise {
        check_shown_parsed_subject(subject, protein, taxdb_found, ordinal, program)?;
    }
    if shown.tabular {
        subject
            .tabular_ids(protein)
            .map_err(|problem| parsed_refusal(ordinal, &format!("has {problem}"), program))?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40 (the fields of
// `-outfmt 6` without a list)
// ```c++
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// ```
/// The subject id of the default tabular fields (`saccver`) of each reported database sequence by
/// rank, for the programs that write no custom fields (BLASTN, TBLASTN, TBLASTX): the first word
/// of the stored title (`tabular_subject_id`, the same for `sseqid`, `sacc` and `saccver`), or for
/// a sequence with parsed Seq-ids (`parsed[rank]`) its `saccver` (`ParsedSubject::tabular_ids`).
/// `records[rank]` holds the stored title, `ordinals[rank]` the ordinal id; `protein` for a
/// protein database.
pub(crate) fn tabular_subject_ids(
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    ordinals: &[u32],
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    protein: bool,
    program: &str,
) -> Result<Vec<std::sync::Arc<[u8]>>> {
    records
        .iter()
        .zip(ordinals)
        .enumerate()
        .map(
            |(rank, (record, &ordinal))| match parsed.get(rank).and_then(Option::as_deref) {
                Some(subject) => {
                    let [_, _, accession_version] =
                        subject.tabular_ids(protein).map_err(|problem| {
                            parsed_refusal(ordinal, &format!("has {problem}"), program)
                        })?;
                    Ok(std::sync::Arc::from(accession_version))
                }
                None => tabular_subject_id(&record.title, ordinal, protein, program)
                    .map(std::sync::Arc::from),
            },
        )
        .collect()
}

/// Why a report of some sequence of a database could be refused (`reports_may_refuse`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReportRefusal {
    /// A header that does not read, or a sequence whose stored title, tabular id or parsed
    /// Seq-ids the reports refuse.
    Headers,
    /// Two sequences whose Seq-align ids NCBI takes for one (`refuse_joined_hits`).
    Joined,
    /// A taxid that does not decode, or a lookup of the taxonomy database that could be refused,
    /// where the batches' prefetch reads them (`prefetch_messages`).
    Taxonomy,
}

// No NCBI counterpart: LOSAT's check whether one of its own refusals can come where a report of
// the search names a sequence of the database; it changes no value NCBI computes.
/// Whether a report that names any sequence of `db` could be refused by LOSAT, with the reports
/// `shown`: every sequence's header is read as the reports read it (`read_header`) and checked as
/// they check it (`check_reported_title`, `tabular_subject_id` for the tabular formats,
/// `check_reported_parsed` with `taxdb_found` where outfmt 0 is shown, as `reported_subjects`
/// asks). Where the reports join following hits (`shown.joined`), no two sequences with parsed
/// Seq-ids may be joined (`refuse_joined_hits`): the sequences that share a key of their Seq-align
/// ids (`db_subject::join_keys`; joined ones always share one) are read again and compared
/// (`joined_in_reports`). With `prefetch` (the batches' prefetch of a search with `taxdb.bti` on
/// the search path) the taxids that the prefetch looks up must read (a stored title's first one,
/// every one of a sequence with parsed Seq-ids, as `prefetch_messages` reads them) and no lookup
/// of them may be refused (`taxids_may_refuse`). A reader error is a refusal too. The headers are
/// read in `pool`; nothing is written and no state of the run changes. `None`: no report of any
/// sequence can be refused.
pub(crate) fn reports_may_refuse(
    db: &SeqDb,
    shown: ParsedShown,
    protein: bool,
    taxdb_found: bool,
    prefetch: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
    program: &str,
) -> Option<ReportRefusal> {
    /// What the pass keeps: the join keys with their sequence, and the prefetch's taxids.
    #[derive(Default)]
    struct Kept {
        keys: Vec<(u64, u32)>,
        taxids: Vec<i32>,
    }
    impl Kept {
        fn join(mut self, mut other: Kept) -> std::result::Result<Kept, ReportRefusal> {
            self.keys.append(&mut other.keys);
            self.taxids.append(&mut other.taxids);
            Ok(self)
        }
    }
    let check = |mut kept: Kept, ordinal: u32| -> std::result::Result<Kept, ReportRefusal> {
        match read_header(db, ordinal, 0, program) {
            Ok(ReportedHeader::Stored(stored)) => {
                check_reported_title(stored.replaced_bytes, ordinal, program)
                    .map_err(|_| ReportRefusal::Headers)?;
                if shown.tabular {
                    tabular_subject_id(&stored.title, ordinal, protein, program)
                        .map_err(|_| ReportRefusal::Headers)?;
                }
                if prefetch {
                    let taxids = db
                        .defline_taxids(ordinal)
                        .map_err(|_| ReportRefusal::Taxonomy)?;
                    kept.taxids
                        .extend(taxids.first().copied().filter(|&taxid| taxid > 0));
                }
            }
            Ok(ReportedHeader::Parsed(subject)) => {
                check_reported_parsed(
                    &subject,
                    shown,
                    protein,
                    shown.pairwise && taxdb_found,
                    ordinal,
                    program,
                )
                .map_err(|_| ReportRefusal::Headers)?;
                if shown.joined {
                    let mut keys = Vec::new();
                    crate::report::db_subject::join_keys(&subject, &mut keys);
                    kept.keys.extend(keys.into_iter().map(|key| (key, ordinal)));
                }
                if prefetch {
                    kept.taxids.extend(subject.volume_taxids());
                }
            }
            Err(_) => return Err(ReportRefusal::Headers),
        }
        Ok(kept)
    };
    let all = || -> std::result::Result<Kept, ReportRefusal> {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            use rayon::prelude::*;
            return pool.install(|| {
                (0..db.num_seqs())
                    .into_par_iter()
                    .with_min_len(4096)
                    .try_fold(Kept::default, check)
                    .try_reduce(Kept::default, Kept::join)
            });
        }
        let _ = pool;
        (0..db.num_seqs()).try_fold(Kept::default(), check)
    };
    let mut kept = match all() {
        Ok(kept) => kept,
        Err(refusal) => return Some(refusal),
    };
    kept.keys.sort_unstable();
    for shared_key in kept.keys.chunk_by(|a, b| a.0 == b.0) {
        let mut ordinals: Vec<u32> = shared_key.iter().map(|&(_, ordinal)| ordinal).collect();
        ordinals.dedup();
        if ordinals.len() < 2 {
            continue;
        }
        let mut subjects = Vec::with_capacity(ordinals.len());
        for ordinal in ordinals {
            match read_header(db, ordinal, 0, program) {
                Ok(ReportedHeader::Parsed(subject)) => subjects.push(subject),
                _ => return Some(ReportRefusal::Headers),
            }
        }
        for (index, a) in subjects.iter().enumerate() {
            if subjects[index + 1..]
                .iter()
                .any(|b| crate::report::db_subject::joined_in_reports(a, b))
            {
                return Some(ReportRefusal::Joined);
            }
        }
    }
    kept.taxids.sort_unstable();
    kept.taxids.dedup();
    (prefetch && taxids_may_refuse(db, &kept.taxids)).then_some(ReportRefusal::Taxonomy)
}

/// Whether the lookup of one of `taxids` in the taxonomy database on the search path of `db`
/// could be refused, as the batches' prefetch looks them up (`prefetch_messages`); `false`
/// without `taxdb.bti` there. Each one is looked up once with a taxonomy handle of its own
/// (as `taxdb_may_refuse`'s database side does), so the run's handle and messages are as they
/// would be without this call. The prefetch's refusal of NCBI's cache
/// (`TaxdbRun::cache_collision`) comes at the run's first prefetch, before any report.
fn taxids_may_refuse(db: &SeqDb, taxids: &[i32]) -> bool {
    if !db.taxdb_found() {
        return false;
    }
    let Some(search_path) = db.search_path() else {
        return false;
    };
    let taxdb = crate::seqdb::taxdb::LazyTaxDb::new(search_path.clone());
    let mut messages = Vec::new();
    taxids
        .iter()
        .any(|&taxid| taxdb.names(taxid, &mut messages).is_err())
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
// number of readings is the formatter's, not ported; makeblastdb and `losat makeblastdb`
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

/// Byte 0 of NCBIstdaa: the gap character (`-`), and the sentinel byte between the
/// sequences of a protein volume's sequence file.
const NCBISTDAA_GAP: u8 = 0;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593
// ```c
//             /* skip special characters */
//             if (i == AMINOACID_TO_NCBISTDAA['U'] ||
//                 i == AMINOACID_TO_NCBISTDAA['O'] ||
//                 i == AMINOACID_TO_NCBISTDAA['-'] ||
// ```
// A stored gap character (0, which makeblastdb writes for `-` in a protein FASTA)
// never aligns in NCBI (row and column 0 stay `BLAST_SCORE_MIN`) and NCBI's formatter
// then counts it as a gap of the alignment; LOSAT's kernels read byte 0 as the end of a
// sequence, so a sequence that holds one is rejected. A byte that is no NCBIstdaa
// residue (28 or more) indexes NCBI's matrix past its rows and is rejected too.
/// LOSAT's refusal of the protein database sequence `ordinal` (`with_sentinels`: its
/// residues with the byte before and the byte after them in the sequence file) when
/// `program` cannot search it as NCBI does: a sequence not between two 0 bytes, one
/// that holds the gap character, or one with a byte that is no NCBIstdaa residue.
pub(crate) fn check_protein_sequence(
    with_sentinels: &[u8],
    ordinal: u64,
    program: &str,
) -> Result<()> {
    let (first, last) = (with_sentinels[0], with_sentinels[with_sentinels.len() - 1]);
    if first != NCBISTDAA_GAP || last != NCBISTDAA_GAP {
        bail!(
            "database sequence {ordinal} is not between two NUL bytes in its sequence file (the volume is damaged), which is not supported by LOSAT's {program}"
        );
    }
    let residues = &with_sentinels[1..with_sentinels.len() - 1];
    // Every residue must be a byte other than the gap and below BLASTAA_SIZE: with the gap byte
    // 0 (asserted below), the residues 1..BLASTAA_SIZE, the bytes whose `residue - 1` (wrapping,
    // so the gap becomes 255) is below BLASTAA_SIZE - 1. The largest such difference, taken over
    // all bytes without an early exit, is a loop the compiler vectorizes; the early exit of
    // `all` made it a loop over single bytes (LOSAT-only cost; NCBI checks nothing here).
    const _: () = assert!(
        NCBISTDAA_GAP == 0
            && crate::utils::matrix::BLASTAA_SIZE > 1
            && crate::utils::matrix::BLASTAA_SIZE <= 256
    );
    const FIRST: u8 = NCBISTDAA_GAP + 1;
    const LETTERS: usize = crate::utils::matrix::BLASTAA_SIZE - FIRST as usize;
    let largest = residues
        .iter()
        .map(|&residue| residue.wrapping_sub(FIRST))
        .max()
        .unwrap_or(0);
    if usize::from(largest) < LETTERS {
        return Ok(());
    }
    if let Some(position) = residues
        .iter()
        .position(|&residue| residue == NCBISTDAA_GAP)
    {
        bail!(
            "database sequence {ordinal} holds a gap character ('-') at position {}; a protein database sequence with gap characters is not supported by LOSAT's {program} (NCBI BLAST+ never aligns the character and then counts it as a gap of the alignment)",
            position + 1
        );
    }
    if let Some(position) = residues
        .iter()
        .position(|&residue| usize::from(residue) >= crate::utils::matrix::BLASTAA_SIZE)
    {
        bail!(
            "database sequence {ordinal} holds a byte that is no NCBIstdaa residue at position {}, which is not supported by LOSAT's {program}",
            position + 1
        );
    }
    Ok(())
}

/// Bytes of the sequence file that one read of `read_protein_stores` holds at most
/// (`SeqDb::next_block`). No NCBI counterpart: NCBI maps the volume files. The size
/// changes nothing but the number of reads.
const PROTEIN_STORE_BLOCK_BYTES: u64 = 1 << 20;

/// The sequences of a protein database as its sequence files store them: NCBIstdaa
/// bytes, each sequence followed by the sentinel byte 0 (and the first preceded by
/// one), in ordinal order, with the totals of the index files. NCBI searches these
/// bytes in place (blast_engine.c:1404, `eBlastEncodingProtein`), and its traceback and
/// composition redo read them again; a program that holds them all searches, redoes and
/// shows the bytes NCBI does.
#[derive(Debug)]
pub(crate) struct StoredProteins {
    /// The sentinel, then each sequence and the sentinel after it.
    bytes: Vec<u8>,
    /// Where each sequence's residues start in `bytes`, and one more entry: sequence
    /// `i` is `bytes[starts[i]..starts[i + 1] - 1]`.
    starts: Vec<usize>,
    /// The total length of the index files (`CSeqDB::GetTotalLength`).
    total_length: u64,
}

impl StoredProteins {
    /// The number of sequences.
    pub(crate) fn len(&self) -> usize {
        self.starts.len() - 1
    }

    /// Whether the database holds no sequence.
    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The length of sequence `ordinal`.
    pub(crate) fn length(&self, ordinal: usize) -> usize {
        self.starts[ordinal + 1] - self.starts[ordinal] - 1
    }

    /// The residues of sequence `ordinal` followed by the sentinel byte 0.
    pub(crate) fn residues_with_sentinel(&self, ordinal: usize) -> &[u8] {
        &self.bytes[self.starts[ordinal]..self.starts[ordinal + 1]]
    }

    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): `residues_with_sentinel` with the sentinel before the
    /// residues (the store's first byte, or the sentinel after the sequence before).
    pub(crate) fn x_residues_with_sentinels(&self, ordinal: usize) -> &[u8] {
        &self.bytes[self.starts[ordinal] - 1..self.starts[ordinal + 1]]
    }

    /// The total length of the index files.
    pub(crate) fn total_length(&self) -> u64 {
        self.total_length
    }

    /// The length of the longest sequence (0 for an empty database).
    pub(crate) fn max_length(&self) -> usize {
        (0..self.len()).map(|i| self.length(i)).max().unwrap_or(0)
    }
}

/// Reads every sequence of the protein database `db` once, in blocks of
/// `PROTEIN_STORE_BLOCK_BYTES` in ordinal order, as its sequence files store them. A
/// reader error is NCBI's error of the preliminary search (`search_read_error`: NCBI
/// reads every sequence there); a sequence `program` cannot search as NCBI does is
/// LOSAT's refusal (`check_protein_sequence`).
pub(crate) fn read_protein_stores(db: &SeqDb, program: &str) -> Result<StoredProteins> {
    let sequences = db.num_seqs() as usize;
    let mut bytes = Vec::with_capacity(
        usize::try_from(db.total_length())
            .unwrap_or(0)
            .saturating_add(sequences + 1),
    );
    let mut starts = Vec::with_capacity(sequences + 1);
    bytes.push(NCBISTDAA_GAP);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, PROTEIN_STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let with_sentinels = block.protein_with_sentinels(i);
            check_protein_sequence(
                with_sentinels,
                u64::from(range.first_ordinal) + i as u64,
                program,
            )?;
            starts.push(bytes.len());
            bytes.extend_from_slice(&with_sentinels[1..]);
        }
        next_ordinal = range.end_ordinal();
    }
    starts.push(bytes.len());
    Ok(StoredProteins {
        bytes,
        starts,
        total_length: db.total_length(),
    })
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAM): `read_protein_stores` without the store: the same reads in
/// the same blocks and the same checks in the same order, with the same errors, the number of
/// residues of each sequence handed to `each` (ordinal id, residues) and its bytes dropped with
/// its block. A search that reads the sequences again later (BLASTX's passes,
/// `blastx/x_db_stream.rs`) knows from it that they read, and what the store's
/// `max_length` would be.
pub(crate) fn x_check_protein_stores(
    db: &SeqDb,
    program: &str,
    mut each: impl FnMut(u32, usize),
) -> Result<()> {
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, PROTEIN_STORE_BLOCK_BYTES) {
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let with_sentinels = block.protein_with_sentinels(i);
            let ordinal = range.first_ordinal + i as u32;
            check_protein_sequence(with_sentinels, u64::from(ordinal), program)?;
            each(ordinal, with_sentinels.len() - 2);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::matrix::BLASTAA_SIZE;

    /// The residues between the two sentinel bytes.
    fn stored(residues: &[u8]) -> Vec<u8> {
        let mut bytes = vec![NCBISTDAA_GAP];
        bytes.extend_from_slice(residues);
        bytes.push(NCBISTDAA_GAP);
        bytes
    }

    fn check(residues: &[u8]) -> std::result::Result<(), String> {
        check_protein_sequence(&stored(residues), 7, "blastx").map_err(|e| e.to_string())
    }

    /// The check before the single pass (an `all` over the residues).
    fn accepted_before(residues: &[u8]) -> bool {
        residues
            .iter()
            .all(|&residue| residue != NCBISTDAA_GAP && usize::from(residue) < BLASTAA_SIZE)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593 (quoted at
    // `check_protein_sequence`): the gap byte and bytes past the matrix are what the check
    // refuses. Not a port of an NCBI unit case.
    #[test]
    fn the_check_accepts_the_letters_and_names_the_first_bad_byte() {
        let first = NCBISTDAA_GAP + 1;
        let last = (BLASTAA_SIZE - 1) as u8;
        assert_eq!(check(&[]), Ok(()));
        assert_eq!(check(&[first]), Ok(()));
        assert_eq!(check(&[last]), Ok(()));
        let letters: Vec<u8> = (0..1_000).map(|i| [first, last][i % 2]).collect();
        assert_eq!(check(&letters), Ok(()));
        // (byte, words of the refusal, the character after its position)
        for (bad, says, after) in [
            (NCBISTDAA_GAP, "gap character", ';'),
            (BLASTAA_SIZE as u8, "no NCBIstdaa residue", ','),
            (255, "no NCBIstdaa residue", ','),
        ] {
            for length in [1, 2, 31, 32, 33, 100, 1_000] {
                for at in [0, length / 2, length - 1] {
                    let mut residues = letters[..length].to_vec();
                    residues[at] = bad;
                    let error = check(&residues).unwrap_err();
                    assert!(
                        error.contains(says)
                            && error.contains(&format!("position {}{after}", at + 1)),
                        "byte {bad} at {at} of {length}: {error}"
                    );
                }
            }
        }
        // A gap before a byte past the matrix is reported as the gap, as before.
        let mut residues = letters[..10].to_vec();
        residues[7] = NCBISTDAA_GAP;
        residues[3] = 200;
        assert!(check(&residues).unwrap_err().contains("gap character"));
        // The sentinels.
        assert!(check_protein_sequence(&[1, first, 0], 7, "blastx").is_err());
        assert!(check_protein_sequence(&[0, first, 1], 7, "blastx").is_err());
    }

    // Random residues around the letters' bounds: the single pass accepts exactly what the
    // `all` accepted.
    #[test]
    fn the_check_accepts_what_the_check_before_it_accepted() {
        let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for case in 0..2_000 {
            let length = (next() % 200) as usize;
            let mut residues: Vec<u8> = (0..length).map(|_| 1 + (next() % 27) as u8).collect();
            if case % 2 == 0 && length > 0 {
                let at = (next() as usize) % length;
                residues[at] = [0, 27, 28, 29, 254, 255][(next() % 6) as usize];
            }
            assert_eq!(
                check(&residues).is_ok(),
                accepted_before(&residues),
                "{residues:?}"
            );
        }
    }
}

/// The taxonomy database of a run (`CTaxDBFileInfo`, a static of NCBI's `GetTaxNames`) and the
/// taxids the volumes have looked up (`CSeqDBVol::x_GetTaxonomy` caches them per volume).
///
/// The volume lookups are a set: a (volume, taxid) is looked up once per run. NCBI's cache of a
/// volume is not a set but 256 slots, one per `taxid & 255` (`CSeqDBVol::m_TaxCache`, below), so
/// a volume looks a taxid up again after another taxid of its slot. The only lookup whose
/// repetition writes a text is that of the taxid of the offset error
/// (`TaxDb::offset_error_taxid`), so the set gives NCBI's bytes unless a volume holds that taxid
/// and another of its slot; such a run is refused at its first prefetch
/// ([`TaxdbRun::cache_collision`]).
///
/// It lives on the database of the search, so a `run_local` call that writes outfmt 15 and 16
/// of one `-db` search shares it: the open and the volume-lookup messages are written for the
/// first of the two only (the CLI writes one format). NCBI, which runs each format alone, writes
/// them for both; the ABI callers that ask for both do not see them twice.
#[derive(Default)]
pub(crate) struct TaxdbRun {
    lazy: std::sync::Mutex<Option<std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>>>,
    volume_lookups: std::sync::Mutex<std::collections::HashSet<(usize, i32)>>,
    found: std::sync::OnceLock<bool>,
    /// What the lookups of the run can write ([`TaxdbRun::prefetch_texts`]).
    texts: std::sync::OnceLock<PrefetchTexts>,
    /// The refusal of a run whose prefetch NCBI's cache makes write the offset error more than
    /// once, as a reader's refusal ([`TaxdbRun::cache_collision`]).
    collision: std::sync::OnceLock<Option<(std::path::PathBuf, String)>>,
}

/// What the lookups of a run's taxonomy database can write ([`TaxdbRun::prefetch_texts`]).
#[derive(Clone, Copy)]
struct PrefetchTexts {
    /// The positive taxid whose lookup writes the offset error (`TaxDb::offset_error_taxid`).
    offset_error_taxid: Option<i32>,
}

impl TaxdbRun {
    /// [`SeqDb::taxdb_found`] of the run's database, asked once: the search path is not
    /// looked through again for each batch.
    pub(crate) fn found(&self, db: &SeqDb) -> bool {
        *self.found.get_or_init(|| db.taxdb_found())
    }

    /// [`PrefetchTexts`] of the run, read once.
    fn prefetch_texts(&self, db: &SeqDb) -> PrefetchTexts {
        *self.texts.get_or_init(|| {
            let refused = PrefetchTexts {
                offset_error_taxid: None,
            };
            if !self.found(db) {
                return PrefetchTexts {
                    offset_error_taxid: None,
                };
            }
            let Some(search_path) = db.search_path() else {
                return PrefetchTexts {
                    offset_error_taxid: None,
                };
            };
            let Ok(taxdb) = crate::seqdb::taxdb::TaxDb::open(search_path) else {
                return refused;
            };
            let Ok(taxid) = taxdb.offset_error_taxid() else {
                return refused;
            };
            // `x_GetTaxonomy` looks a positive taxid up only.
            let taxid = taxid.filter(|&taxid| taxid > 0);
            PrefetchTexts {
                offset_error_taxid: taxid,
            }
        })
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbgeneral.hpp:1184-1196
    // ```c++
    //     TValue & Lookup(int key)
    //     {
    //         _ASSERT(IS_POWER_OF_TWO(m_Slots.size()));
    //
    //         TSlot & slot = m_Slots[key & (m_Slots.size()-1)];
    //
    //         if (slot.first != key) {
    //             slot.first = key;
    //             slot.second = TValue();
    //         }
    //
    //         return slot.second;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:84 (the
    // constructor's `m_TaxCache     (256),`) and 1142-1182
    // ```c++
    //         if (use_taxinfo_cache && m_TaxCache.Lookup(TAX_ID_TO(int, taxid)).NotEmpty()) {
    //             have_org_desc = true;
    //         }
    //
    //         SSeqDBTaxInfo tnames(taxid);
    //         bool found_taxid_in_taxonomy_blastdb = true;
    //
    //         if ((! have_org_desc) && provide_new_taxonomy_info) {
    //             try {
    //                 found_taxid_in_taxonomy_blastdb = CSeqDBTaxInfo::GetTaxNames(taxid, tnames);
    //     ...
    //                 if (use_taxinfo_cache) {
    //                     m_TaxCache.Lookup(TAX_ID_TO(int, taxid)) = desc;
    //                 }
    // ```
    // A volume's cache (`CSeqDBIntCache`) keeps one taxid per slot `taxid & 255`; a fetch of a
    // sequence with another taxid of the slot clears it, and the next fetch of the first taxid
    // asks `GetTaxNames` again. Of the texts of the lookups only the offset error repeats with it
    // (the texts of opening come once per process). How often it repeats follows the order in
    // which the object manager fetches the batch's subjects, which is neither the order of the
    // reports nor that of the ordinal ids (review 1 of item T: one query on a volume of the taxids
    // 99990, 150, 99990 writes the message twice, two queries once; fix 1: the same query once
    // with a version 4 database or an absolute `BLASTDB`), so LOSAT does not model it.
    /// LOSAT's refusal of a run whose prefetch NCBI would make write the offset error more than
    /// once in a volume: the taxonomy database has the offset error of a positive taxid `L`
    /// ([`TaxdbRun::prefetch_texts`]), and a volume of the database holds `L` and another
    /// positive taxid of the same cache slot (`T & 255 == L & 255`), in any definition line of any
    /// sequence, whether or not a search reports them. One pass over the headers
    /// ([`SeqDb::visit_defline_set_taxids`] in ranges of [`crate::seqdb::HEADER_RANGE_SEQS`]), made
    /// only for such a damaged taxonomy database and once for the run; a header whose taxids do
    /// not decode is left to the lookups of the batches, which fail on it where it is reported.
    /// `(path of taxdb.bti, reason)` of the refusal, which [`prefetch_messages`] makes as a
    /// reader's refusal ([`TaxonomyRefusal`]) at the run's first prefetch.
    fn cache_collision(&self, db: &SeqDb) -> Option<&(std::path::PathBuf, String)> {
        self.collision
            .get_or_init(|| {
                let last = self.prefetch_texts(db).offset_error_taxid?;
                let (volume, other) = taxonomy_cache_collision(db, last)?;
                let path = db.search_path()?.resolve_file("taxdb.bti")?;
                let volume_name = db
                    .volumes()
                    .get(volume)
                    .map_or_else(String::new, |vol| vol.base().display().to_string());
                Some((
                    std::path::PathBuf::from(path),
                    format!(
                        "taxid {last}, the last entry of the taxonomy database, begins past the \
                         end of taxdb.btd, and the volume {volume_name} holds it and taxid \
                         {other}, which takes the same slot of NCBI BLAST+'s 256-slot taxonomy \
                         cache of a volume (CSeqDBVol::m_TaxCache, CSeqDBIntCache): NCBI writes \
                         \"{}\" again for each fetch of taxid {last} after one of {other}, in an \
                         order of fetches that LOSAT does not reproduce",
                        crate::seqdb::taxdb::MSG_OFFSET_ERROR
                    ),
                ))
            })
            .as_ref()
    }

    /// Whether a prefetch of the run is refused because of NCBI's taxonomy cache
    /// ([`TaxdbRun::cache_collision`]): BLASTX holds its outfmt 5 document then, which the refusal
    /// drops.
    pub(crate) fn prefetch_refused(&self, db: &SeqDb) -> bool {
        self.cache_collision(db).is_some()
    }

    /// The run's taxonomy database, a handle that opens it at its first name asked for.
    fn taxdb(
        &self,
        search_path: &crate::seqdb::SearchPath,
    ) -> std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb> {
        let mut lazy = self.lazy.lock().expect("the taxonomy database");
        std::sync::Arc::clone(lazy.get_or_insert_with(|| {
            std::sync::Arc::new(crate::seqdb::taxdb::LazyTaxDb::new(search_path.clone()))
        }))
    }
}

/// LOSAT's refusal of a taxonomy database name or record (`read_reported_taxids`): a scientific
/// name NCBI's XML stream cannot write, or a record NCBI reads from unmapped memory. A run that
/// holds its XML2 or JSON document drops it on this error and no other (BLASTX,
/// `XmlRefusal::of_taxonomy`); the other errors of the lookups (a read error of a file) are not
/// refusals, and the held document is written as NCBI writes what it has.
#[derive(Debug)]
pub(crate) struct TaxonomyRefusal(String);

impl std::fmt::Display for TaxonomyRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TaxonomyRefusal {}

/// [`taxdb_error`] for the tests of the callers' error handling.
#[cfg(test)]
pub(crate) fn taxdb_error_for_tests(error: SeqDbError) -> anyhow::Error {
    taxdb_error(error, "BLASTX")
}

/// The error of a taxonomy database lookup or open: [`TaxonomyRefusal`] for a refusal of the
/// reader (`SeqDbError::Unsupported`), [`database_error`] for the others.
fn taxdb_error(error: SeqDbError, program: &str) -> anyhow::Error {
    let refused = matches!(error, SeqDbError::Unsupported { .. });
    let error = database_error(error, program);
    if refused {
        anyhow::Error::new(TaxonomyRefusal(format!("{error:#}")))
    } else {
        error
    }
}

/// The Blast-def-line taxid of each reported subject by rank (0 for none and for a subject
/// that is not in `reported`), which the XML2 and JSON formats write in the `Hit` description,
/// and what the taxonomy database says about it. They are read apart from the titles, for those
/// formats only. `ordinals[rank]` is the ordinal id of the subject of rank `rank`.
///
/// NCBI writes the scientific name of a taxid when it finds `taxdb.bti` on the search path
/// (`data4xml2format.cpp x_InitCommon`, `blastxml2_format.cpp:380-386`; validity is not
/// checked): `taxonomy[rank]` is then `Named`, or `NotFound` when the taxonomy database does
/// not know the taxid (the description is the id alone). Without `taxdb.bti`, or without a
/// taxid, `taxonomy` is empty. A scientific name that NCBI's XML stream cannot write, and a
/// taxonomy database that LOSAT cannot read as NCBI does, are refusals made here, before any
/// byte is written.
///
/// NCBI's messages about the taxonomy database go to standard error, and not here: those of
/// opening and of the lookup that the volume makes for each taxid of a fetched subject
/// (`CSeqDBVol::x_GetTaxonomy`, once per volume and taxid) are in the returned [`TaxdbNotes`],
/// which writes them at the first `Hit` that needs them, and those of the lookup of each `Hit`
/// are in `taxonomy` and written at the hit ([`TaxdbPrinter::at_hit`]). `program` is the
/// program that refuses, as LOSAT's refusals name it (`BLASTP`); its messages name it in lower
/// case.
pub(crate) fn read_reported_taxids(
    db: &SeqDb,
    run: &TaxdbRun,
    ordinals: &[u32],
    reported: &std::collections::BTreeSet<usize>,
    program: &str,
) -> Result<(
    Vec<i32>,
    Vec<crate::report::blast_xml::HitTaxonomy>,
    TaxdbNotes,
)> {
    let (taxids, taxonomy, _, notes) =
        read_reported_taxonomy(db, run, ordinals, reported, &[], program)?;
    Ok((taxids, taxonomy, notes))
}

/// [`read_reported_taxids`] of a database whose reported sequences can have parsed Seq-ids
/// (`parsed[rank]`, empty where none has): the third vector holds, by rank, what the taxonomy
/// database said about the taxid of each definition line of such a sequence (empty where nothing
/// was asked), as NCBI's `HitDescr` loop asks it (`SubjectNames::of_parsed`): `Unnamed` for a
/// definition line without a taxid, and no answer after the first `NotFound` (NCBI's loop ends
/// there). The volume's lookups are those of `ParsedSubject::volume_taxids`.
pub(crate) fn read_reported_taxonomy(
    db: &SeqDb,
    run: &TaxdbRun,
    ordinals: &[u32],
    reported: &std::collections::BTreeSet<usize>,
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    program: &str,
) -> Result<(
    Vec<i32>,
    Vec<crate::report::blast_xml::HitTaxonomy>,
    Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
    TaxdbNotes,
)> {
    use crate::report::blast_xml::HitTaxonomy;
    let parsed_subject = |rank: usize| parsed.get(rank).and_then(Option::as_deref);
    let defline_taxids = |subject: &crate::report::db_subject::ParsedSubject| -> Vec<i32> {
        subject
            .deflines()
            .iter()
            .map(|defline| defline.taxid.unwrap_or(0))
            .collect()
    };
    let mut taxids = vec![0; ordinals.len()];
    let mut any_taxid = false;
    for &rank in reported {
        if let Some(subject) = parsed_subject(rank) {
            any_taxid |= defline_taxids(subject).iter().any(|&taxid| taxid != 0);
            continue;
        }
        let defline_taxids = db
            .defline_taxids(ordinals[rank])
            .map_err(|error| database_error(error, program))?;
        // One definition line per reported sequence (`stored_title`).
        taxids[rank] = defline_taxids.first().copied().unwrap_or(0);
        any_taxid |= taxids[rank] != 0;
    }
    if !any_taxid || !db.taxdb_found() {
        return Ok((taxids, Vec::new(), Vec::new(), TaxdbNotes::default()));
    }
    let Some(search_path) = db.search_path() else {
        return Ok((taxids, Vec::new(), Vec::new(), TaxdbNotes::default()));
    };
    let taxdb = run.taxdb(search_path);
    let mut notes = TaxdbNotes {
        program: program.to_ascii_lowercase(),
        ..TaxdbNotes::default()
    };
    // NCBI opens the database at its first lookup; the texts are written at the first `Hit`.
    notes.open = taxdb.open().map_err(|error| taxdb_error(error, program))?;
    let mut asked: std::collections::BTreeMap<i32, HitTaxonomy> = Default::default();
    // The lookup of a `Hit` description (`CSeqDB::GetTaxInfo`).
    let mut hit_lookup = |taxid: i32| -> Result<HitTaxonomy> {
        if !asked.contains_key(&taxid) {
            let mut messages = Vec::new();
            let found = taxdb
                .names(taxid, &mut messages)
                .map_err(|error| taxdb_error(error, program))?;
            let answer = match found {
                Some(names) => {
                    crate::report::blast_xml::names::check_text(&names.scientific_name).map_err(
                        |refusal| {
                            anyhow::Error::new(TaxonomyRefusal(
                                refusal
                                    .error(program, "a taxonomy database scientific name")
                                    .to_string(),
                            ))
                        },
                    )?;
                    HitTaxonomy::Named(names.scientific_name)
                }
                None => HitTaxonomy::NotFound(messages),
            };
            asked.insert(taxid, answer);
        }
        Ok(asked[&taxid].clone())
    };
    let mut taxonomy = vec![HitTaxonomy::Unnamed; taxids.len()];
    let mut parsed_taxonomy = vec![Vec::new(); if parsed.is_empty() { 0 } else { taxids.len() }];
    notes.volumes = vec![0; taxids.len()];
    for &rank in reported {
        let subject = parsed_subject(rank);
        // The volume's lookups when it fetches the sequence (once per volume and taxid;
        // `x_GetTaxonomy` asks for a positive taxid only, seqdbvol.cpp:1136).
        let volume_taxids = match subject {
            Some(subject) if defline_taxids(subject).iter().all(|&taxid| taxid == 0) => continue,
            Some(subject) => subject.volume_taxids(),
            None if taxids[rank] == 0 => continue,
            None => Some(taxids[rank])
                .filter(|&taxid| taxid > 0)
                .into_iter()
                .collect(),
        };
        let volume = db
            .volume_index(ordinals[rank])
            .map_err(|error| database_error(error, program))?;
        notes.volumes[rank] = volume;
        for taxid in volume_taxids {
            let first = run
                .volume_lookups
                .lock()
                .expect("the volume lookups")
                .insert((volume, taxid));
            if first {
                let mut messages = Vec::new();
                taxdb
                    .names(taxid, &mut messages)
                    .map_err(|error| taxdb_error(error, program))?;
                notes.lookups.insert((volume, taxid), messages);
            }
        }
        let Some(subject) = subject else {
            taxonomy[rank] = hit_lookup(taxids[rank])?;
            continue;
        };
        let mut answers = Vec::new();
        for taxid in defline_taxids(subject) {
            let answer = if taxid == 0 {
                HitTaxonomy::Unnamed
            } else {
                hit_lookup(taxid)?
            };
            let stop = matches!(answer, HitTaxonomy::NotFound(_));
            answers.push(answer);
            if stop {
                break;
            }
        }
        parsed_taxonomy[rank] = answers;
    }
    Ok((taxids, taxonomy, parsed_taxonomy, notes))
}

/// NCBI's messages about the taxonomy database when a query batch prefetches the subjects it
/// reports, for the formats that read the subject Bioseqs but show no taxonomy (outfmt 0, 5, 6,
/// 7 and 17): the texts of opening it and of the lookup that a volume makes the first time it
/// fetches a subject with a taxid (`CSeqDBVol::x_GetTaxonomy`, once per volume and positive
/// taxid in the run, `TaxdbRun`), on standard error. `subjects` are the batch's reported
/// sequences (ordinal id, and the definition lines of one with parsed Seq-ids, whose volume
/// lookups are `ParsedSubject::volume_taxids`; the others have one definition line, whose
/// taxid is read here), in any order: the only text of a lookup is the same offset error, and
/// a volume looks the taxid of that error up once in the run unless it holds another taxid of
/// that taxid's slot of NCBI's cache, which is refused (`TaxdbRun::cache_collision`, at the
/// run's first call, before its subjects are looked at, so also in a first batch that reports
/// nothing), so the order of the subjects does not change the bytes. The caller writes the
/// text after the reader's messages of the batch and before the warnings of its first query,
/// where NCBI's prefetch comes. Empty, and nothing read, without `taxdb.bti` on the search
/// path. A taxonomy database that LOSAT cannot read as NCBI does is refused as in
/// [`read_reported_taxonomy`].
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:716-772
/// ```c++
/// s_IsPrefetchFormat(blast::CFormattingArgs::EOutputFormat format_type)
/// {
/// 	if ((format_type == CFormattingArgs::eAsnText) ||
/// 	    (format_type == CFormattingArgs::eAsnBinary) ||
/// 	    (format_type == CFormattingArgs::eArchiveFormat)||
/// 	    (format_type == CFormattingArgs::eJsonSeqalign)) {
/// 		return false;
/// 	}
/// 	return true;
/// }
/// ...
///        s_ExtractSeqidsAndRanges(results, ids, ranges);
///        _TRACE("Prefetching " << ids.size() << " sequence lengths");
///        LoadSequencesToScope(ids, ranges, scope);
/// ```
/// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
/// ```c++
///                 BlastFormatter_PreFetchSequenceData(*results, scope,
///                 		                            fmt_args->GetFormattedOutputChoice());
///                 ITERATE(CSearchResultSet, result, *results) {
///                     formatter.PrintOneResultSet(**result, query_batch);
///                 }
/// ```
pub(crate) fn prefetch_messages<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: impl IntoIterator<Item = (u32, Option<&'s crate::report::db_subject::ParsedSubject>)>,
    program: &str,
) -> Result<Vec<u8>> {
    if !run.found(db) {
        return Ok(Vec::new());
    }
    let Some(search_path) = db.search_path() else {
        return Ok(Vec::new());
    };
    if let Some((path, message)) = run.cache_collision(db) {
        return Err(taxdb_error(
            SeqDbError::Unsupported {
                path: path.clone(),
                message: message.clone(),
            },
            program,
        ));
    }
    let taxdb = run.taxdb(search_path);
    let lower = program.to_ascii_lowercase();
    let mut text = String::new();
    let mut fetched: std::collections::HashSet<u32> = Default::default();
    for (ordinal, parsed) in subjects {
        if !fetched.insert(ordinal) {
            continue;
        }
        // `x_GetTaxonomy` asks for a positive taxid only (seqdbvol.cpp:1136); one definition
        // line per sequence without parsed Seq-ids (`stored_title`).
        let volume_taxids: Vec<i32> = match parsed {
            Some(subject) => subject.volume_taxids(),
            None => db
                .defline_taxids(ordinal)
                .map_err(|error| database_error(error, program))?
                .first()
                .copied()
                .filter(|&taxid| taxid > 0)
                .into_iter()
                .collect(),
        };
        if volume_taxids.is_empty() {
            continue;
        }
        let volume = db
            .volume_index(ordinal)
            .map_err(|error| database_error(error, program))?;
        for taxid in volume_taxids {
            let first = run
                .volume_lookups
                .lock()
                .expect("the volume lookups")
                .insert((volume, taxid));
            if !first {
                continue;
            }
            // The texts of opening come with the run's first lookup (`LazyTaxDb::names`).
            let mut messages = Vec::new();
            taxdb
                .names(taxid, &mut messages)
                .map_err(|error| taxdb_error(error, program))?;
            for line in messages {
                text.push_str(&crate::seqdb::taxdb::err_post_line(&lower, &line));
                text.push('\n');
            }
        }
    }
    Ok(text.into_bytes())
}

/// The volume of `db` that holds the taxid `last` and another positive taxid of the same slot of
/// NCBI's taxonomy cache of a volume (`taxid & 255`, `TaxdbRun::cache_collision`), with that
/// taxid; `None` when no volume does. Every definition line's taxid of every sequence, parsed
/// Seq-ids or not, in one pass over the headers; a header whose taxids do not decode is skipped.
fn taxonomy_cache_collision(db: &SeqDb, last: i32) -> Option<(usize, i32)> {
    /// `CSeqDBVol::m_TaxCache(256)`.
    const SLOTS: i32 = 256;
    let slot = |taxid: i32| taxid & (SLOTS - 1);
    let volumes = db.volumes().len();
    let mut holds_last = vec![false; volumes];
    let mut same_slot: Vec<Option<i32>> = vec![None; volumes];
    let mut found = None;
    for range in db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS) {
        let volume = range.volume;
        if volume >= volumes {
            continue;
        }
        db.visit_defline_set_taxids(range, |_ordinal, taxids| {
            let Ok(taxids) = taxids else {
                return true;
            };
            for taxid in taxids {
                if taxid == last {
                    holds_last[volume] = true;
                } else if taxid > 0 && slot(taxid) == slot(last) {
                    same_slot[volume].get_or_insert(taxid);
                }
            }
            match same_slot[volume] {
                Some(other) if holds_last[volume] => {
                    found = Some((volume, other));
                    false
                }
                _ => true,
            }
        });
        if found.is_some() {
            return found;
        }
    }
    None
}

/// Whether the first format of a call is one whose batches prefetch their reported subjects
/// and that shows no taxonomy of them (`prefetch_messages`): outfmt 0, 5, 6, 7 and 17. The
/// query warnings, which the batch's messages go with, are written with the first format; the
/// XML2 and JSON formats write theirs at each `Hit` (`TaxdbPrinter`), and the Seq-annots and the
/// archive prefetch nothing.
pub(crate) fn prefetches_subjects(outfmt: Option<i32>) -> bool {
    matches!(outfmt, Some(0 | 5 | 6 | 7 | 17))
}

/// How many taxonomy database entries cost as much to check as one database sequence: the
/// taxonomy side is chosen when `entries <= TAXDB_ENTRIES_PER_SEQUENCE * sequences`. Measured
/// (OF5b fix, 1 thread, hot page cache): the pass over NCBI's taxdb takes 0.15 s for 2,314,695
/// entries (about 0.065 us each), the database side 0.25 s for 2,000,000 sequences (about
/// 0.125 us each), so an entry costs about half a sequence.
const TAXDB_ENTRIES_PER_SEQUENCE: u64 = 2;

/// The most database sequences the per-sequence side of [`taxdb_may_refuse`] reads after the
/// taxonomy side found a record that could be refused (about 0.3 us each, so a few seconds);
/// a larger database is held without reading it.
const PREPASS_DB_SEQUENCES_MAX: u64 = 1 << 22;

/// Whether [`read_reported_taxids`] (`names`) or [`prefetch_messages`] (not `names`) could
/// refuse for some sequence of the database, which a run that writes its batches one by one
/// cannot take back (OF5b BLASTX, DECISIONS 32): with `names` a scientific name NCBI's XML
/// stream cannot write (the XML2 and JSON formats write the names; the prefetch of outfmt 0, 5,
/// 6, 7 and 17 asks for none), and in both a record past the last mapped page of `taxdb.btd`,
/// or a header or the taxonomy database that cannot be read. A superset of what a search
/// reports: limits (`-gilist`, `-taxidlist`) and the hits themselves are ignored. `false`
/// without `taxdb.bti` on the search path and for a database without taxids. The prefetch's
/// refusal of NCBI's cache (`TaxdbRun::cache_collision`) comes at the run's first prefetch,
/// before the first batch's reports, and is not asked here.
///
/// The answer comes from whichever side is cheaper: the taxonomy side when the entries of
/// `taxdb.bti` (from its header) number at most [`TAXDB_ENTRIES_PER_SEQUENCE`] times the
/// sequences of the database (`SeqDb::num_seqs`), else the database side:
/// * [`taxdb_side`]: one sequential pass over the taxonomy database asks whether any taxid's
///   lookup could be refused; none -> `false`, without reading the database. This is the
///   normal case for NCBI's taxdb, so the cost does not grow with the database.
/// * [`database_side`]: the first definition line's taxid of every sequence, each distinct
///   non-zero one looked up. It decides the cases the taxonomy side leaves open. When the
///   taxonomy side ran first (the taxonomy database has at most twice as many entries as the
///   database has sequences) and found something, a database of more than
///   [`PREPASS_DB_SEQUENCES_MAX`] sequences is held without this pass.
///
/// "No hold" means that no lookup of the run can refuse after a byte: the taxonomy side proves it
/// for every taxid, the database side for every taxid the database holds.
///
/// It has its own taxonomy handle (`LazyTaxDb`), so the run's [`TaxdbRun`], its volume lookups
/// and the notes that [`TaxdbPrinter`] writes are as they would be without this call. That
/// handle and the run's are two opens of the same files, which differ only if `taxdb.bti` or
/// `taxdb.btd` is replaced while the run goes on (the run opens its handle at its first lookup,
/// possibly many batches after this call): the run can then refuse after bytes were written.
/// NCBI reads the two files at its first lookup too, and a replaced file gives it a mixed
/// result as well; this is not guarded (DECISIONS 32).
pub(crate) fn taxdb_may_refuse(db: &SeqDb, names: bool) -> bool {
    use crate::seqdb::taxdb::LazyTaxDb;
    if !db.taxdb_found() {
        return false;
    }
    let Some(search_path) = db.search_path() else {
        return false;
    };
    let taxdb = LazyTaxDb::new(search_path.clone());
    let sequences = u64::from(db.num_seqs());
    // An index that cannot be opened is left to the database side, which asks for a name only
    // when a sequence has a taxid.
    let taxdb_first = taxdb.entry_count().map_or(false, |entries| {
        entries <= TAXDB_ENTRIES_PER_SEQUENCE * sequences
    });
    if taxdb_first {
        if taxdb_side(&taxdb, names) == Some(false) {
            return false;
        }
        if sequences > PREPASS_DB_SEQUENCES_MAX {
            return true;
        }
    }
    database_side(db, &taxdb, names)
}

/// The taxonomy side of [`taxdb_may_refuse`]: `Some(false)` when no lookup of any taxid can be
/// refused or, with `names`, give a scientific name that `names::check_text` refuses
/// ([`crate::seqdb::taxdb::TaxDb::any_lookup_may_refuse`]); `Some(true)` when one can; `None` when
/// the taxonomy database cannot be opened or read, which only the lookups of the database's taxids
/// can answer.
fn taxdb_side(taxdb: &crate::seqdb::taxdb::LazyTaxDb, names: bool) -> Option<bool> {
    taxdb
        .any_lookup_may_refuse(&|name| {
            names && crate::report::blast_xml::names::check_text(name).is_err()
        })
        .ok()
}

/// The database side of [`taxdb_may_refuse`]: whether the lookup of the first definition line's
/// taxid of some sequence would be refused, looking each distinct non-zero taxid up with the
/// `names` + `names::check_text` code `read_reported_taxids` uses (`names::check_text` only
/// with `names`). The headers are read in
/// ranges of [`crate::seqdb::HEADER_RANGE_SEQS`] sequences, one read per MiB
/// ([`SeqDb::visit_defline_taxids`]). A header that cannot be read answers `true`. Memory: the
/// distinct taxids, nothing per sequence.
fn database_side(db: &SeqDb, taxdb: &crate::seqdb::taxdb::LazyTaxDb, names: bool) -> bool {
    let mut seen: std::collections::HashSet<i32> = Default::default();
    let mut refused = false;
    for range in db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS) {
        let whole = db.visit_defline_taxids(range, |_ordinal, defline_taxids| {
            let Ok(defline_taxids) = defline_taxids else {
                refused = true;
                return false;
            };
            // One definition line per sequence (`read_reported_taxids`).
            let taxid = defline_taxids.first().copied().unwrap_or(0);
            if taxid == 0 || !seen.insert(taxid) {
                return true;
            }
            let mut messages = Vec::new();
            refused = match taxdb.names(taxid, &mut messages) {
                Err(_) => true,
                Ok(Some(found)) => {
                    names
                        && crate::report::blast_xml::names::check_text(&found.scientific_name)
                            .is_err()
                }
                Ok(None) => false,
            };
            !refused
        });
        if !whole {
            return refused;
        }
    }
    false
}

/// NCBI's messages about the taxonomy database that are not the lookup of a `Hit`: the texts of
/// opening it (at the first lookup of the run) and the lookup that a volume makes the first
/// time it fetches a subject with a taxid (`CSeqDBVol::x_GetTaxonomy`). NCBI writes them when
/// the first `Hit` that needs them is formatted, after the reader's messages of its batch, so
/// they are kept here and written by [`TaxdbPrinter::at_hit`].
#[derive(Default)]
pub(crate) struct TaxdbNotes {
    /// The program's name in lower case, as the messages show it.
    program: String,
    /// The texts of opening the database, when this batch opened it.
    open: Vec<String>,
    /// The volume of each rank with a taxid.
    volumes: Vec<usize>,
    /// The texts of the lookup of each (volume, taxid) that no earlier batch had looked up.
    lookups: std::collections::BTreeMap<(usize, i32), Vec<String>>,
}

impl TaxdbNotes {
    /// A printer that writes these notes as the `Hit`s of a report come.
    pub(crate) fn printer(&self) -> TaxdbPrinter<'_> {
        TaxdbPrinter {
            notes: self,
            opened: false,
            written: Default::default(),
        }
    }
}

/// Writes the taxonomy database's messages of a report in NCBI's order: the notes of the
/// taxid of a `Hit` before the messages of the `Hit`'s own lookup.
pub(crate) struct TaxdbPrinter<'a> {
    notes: &'a TaxdbNotes,
    opened: bool,
    written: std::collections::HashSet<(usize, i32)>,
}

impl TaxdbPrinter<'_> {
    /// NCBI's messages at the `Hit` of subject `rank` (taxid `taxid`, `taxonomy` what its lookup
    /// said): the texts of opening and of the volume's lookup, where they first apply, then the
    /// messages of the `Hit`'s lookup (`CSeqDB::GetTaxInfo`, every time), to `sink`.
    pub(crate) fn at_hit(
        &mut self,
        rank: usize,
        taxid: i32,
        taxonomy: &crate::report::blast_xml::HitTaxonomy,
        sink: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        self.write(rank, taxid != 0, &[taxid], std::iter::once(taxonomy), sink)
    }

    /// `at_hit` of the `Hit` of subject `rank` named `names`: a subject with parsed Seq-ids
    /// has the volume's lookup of each of its taxids (`SubjectNames::volume_taxids`) and the
    /// lookup of each `HitDescr`.
    pub(crate) fn at_names(
        &mut self,
        rank: usize,
        names: &crate::report::blast_xml::subjects::SubjectNames,
        sink: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        use crate::report::blast_xml::HitTaxonomy;
        if names.descriptions.is_empty() {
            return self.at_hit(rank, names.taxid, &names.taxonomy, sink);
        }
        let asks = !names.volume_taxids.is_empty()
            || names.descriptions.iter().any(|description| {
                description.taxid != 0 || matches!(description.taxonomy, HitTaxonomy::NotFound(_))
            });
        let lookups = names
            .descriptions
            .iter()
            .map(|description| &description.taxonomy);
        self.write(rank, asks, &names.volume_taxids, lookups, sink)
    }

    /// The messages of a `Hit` that asks the taxonomy database (`asks`): opening and the
    /// volume's lookup of each of `volume_taxids`, where they first apply, then the messages of
    /// the `Hit`'s own lookups (`lookups`, in order).
    fn write<'t>(
        &mut self,
        rank: usize,
        asks: bool,
        volume_taxids: &[i32],
        lookups: impl Iterator<Item = &'t crate::report::blast_xml::HitTaxonomy>,
        sink: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        use crate::seqdb::taxdb::err_post_line;
        let mut text = String::new();
        let program = &self.notes.program;
        let mut push = |line: &str| {
            text.push_str(&err_post_line(program, line));
            text.push('\n');
        };
        if asks {
            if !std::mem::replace(&mut self.opened, true) {
                self.notes.open.iter().for_each(|line| push(line));
            }
            if let Some(volume) = self.notes.volumes.get(rank).copied() {
                for &taxid in volume_taxids {
                    if let Some(lines) = self.notes.lookups.get(&(volume, taxid)) {
                        if self.written.insert((volume, taxid)) {
                            lines.iter().for_each(|line| push(line));
                        }
                    }
                }
            }
        }
        for taxonomy in lookups {
            if let crate::report::blast_xml::HitTaxonomy::NotFound(messages) = taxonomy {
                messages.iter().for_each(|line| push(line));
            }
        }
        if text.is_empty() {
            return Ok(());
        }
        sink.write_all(text.as_bytes())?;
        sink.flush()
    }
}

#[cfg(test)]
mod taxdb_notes_tests {
    use super::*;
    use crate::report::blast_xml::HitTaxonomy;

    fn lines(printer: &mut TaxdbPrinter<'_>, rank: usize, taxid: i32, hit: &HitTaxonomy) -> String {
        let mut sink = Vec::new();
        printer.at_hit(rank, taxid, hit, &mut sink).unwrap();
        String::from_utf8(sink).unwrap()
    }

    /// NCBI writes the texts of opening at the first `Hit` with a taxid, the volume's lookup
    /// at the first `Hit` of each (volume, taxid), and the `Hit`'s own lookup every time, each
    /// after the lines before it (and so after the reader's messages of the batch).
    #[test]
    fn the_notes_come_at_the_first_hit_that_needs_them() {
        let notes = TaxdbNotes {
            program: "blastp".to_string(),
            open: vec!["Error: open".to_string()],
            volumes: vec![0, 0, 1, 0],
            lookups: [
                ((0, 9606), vec!["Error: lookup 9606".to_string()]),
                ((1, 9606), vec!["Error: lookup 9606".to_string()]),
            ]
            .into_iter()
            .collect(),
        };
        let mut printer = notes.printer();
        let missing = HitTaxonomy::NotFound(vec!["Error: hit".to_string()]);
        // a subject without a taxid writes nothing and does not open the database
        assert_eq!(lines(&mut printer, 3, 0, &HitTaxonomy::Unnamed), "");
        assert_eq!(
            lines(&mut printer, 0, 9606, &missing),
            "Error: [blastp] Error: open\nError: [blastp] Error: lookup 9606\nError: [blastp] Error: hit\n"
        );
        // the same (volume, taxid): the volume does not look it up again, the `Hit` does
        assert_eq!(
            lines(&mut printer, 1, 9606, &missing),
            "Error: [blastp] Error: hit\n"
        );
        assert_eq!(
            lines(&mut printer, 1, 9606, &HitTaxonomy::Named(Vec::new())),
            ""
        );
        // another volume: its own lookup
        assert_eq!(
            lines(&mut printer, 2, 9606, &HitTaxonomy::Named(Vec::new())),
            "Error: [blastp] Error: lookup 9606\n"
        );
        // no notes (a batch after the first, or a database without a taxonomy database)
        let none = TaxdbNotes::default();
        assert_eq!(
            lines(&mut none.printer(), 0, 9606, &HitTaxonomy::Unnamed),
            ""
        );
    }
}

#[cfg(test)]
mod taxdb_prepass_tests {
    //! `taxdb_may_refuse`: the answer of the lookup of every sequence's taxid, with the
    //! databases of `tests/fixtures/outfmt16/db` (`sprot_tax`: two sequences, taxid 9606;
    //! `sprot_mixtax`: those two with taxid 10090 and a third with taxid 9606) and the taxonomy
    //! databases of `tests/seqdb/taxdb`.
    use super::*;
    use crate::seqdb::alias::SearchPath;
    use crate::seqdb::SeqType;
    use std::path::{Path, PathBuf};

    fn root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
    }

    /// The database `name` of `tests/fixtures/outfmt16/db` opened with `BLASTDB=taxdb`
    /// (`None`: not set), from an empty current directory.
    fn database(name: &str, taxdb: Option<&Path>) -> SeqDb {
        let cwd = std::env::temp_dir().join(format!("losat_prepass_cwd_{}", std::process::id()));
        std::fs::create_dir_all(&cwd).unwrap();
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
        let path = root().join("fixtures/outfmt16/db").join(name);
        SeqDb::open_in(path.to_str().unwrap(), SeqType::Protein, &search).unwrap()
    }

    fn variant(name: &str) -> PathBuf {
        root().join("seqdb/taxdb").join(name)
    }

    /// An empty current directory of the test `tag`, which the test removes at its end.
    fn empty_cwd(tag: &str) -> PathBuf {
        let cwd = std::env::temp_dir().join(format!("losat_{tag}_cwd_{}", std::process::id()));
        std::fs::create_dir_all(&cwd).unwrap();
        cwd
    }

    /// The subject of each sequence of a database with parsed Seq-ids, by ordinal id.
    fn parsed_subjects(db: &SeqDb) -> Vec<crate::report::db_subject::ParsedSubject> {
        (0..db.num_seqs())
            .map(|ordinal| {
                let set = db.defline_set(ordinal).unwrap();
                crate::report::db_subject::ParsedSubject::new(set, db.is_version5()).unwrap()
            })
            .collect()
    }

    /// The database `name` of `tests/db_search` opened with `BLASTDB=taxdb` (`None`: not set),
    /// from the empty current directory `cwd`.
    fn db_search(name: &str, taxdb: Option<&Path>, cwd: &Path) -> SeqDb {
        let search = SearchPath::new(cwd, taxdb.map(Path::as_os_str), None).unwrap();
        let path = root().join("db_search").join(name);
        SeqDb::open_in(path.to_str().unwrap(), SeqType::Protein, &search).unwrap()
    }

    #[test]
    fn names_that_can_be_written_do_not_hold() {
        for name in ["sprot_tax", "sprot_mixtax"] {
            let db = database(name, Some(&variant("good")));
            assert!(db.taxdb_found());
            assert!(!taxdb_may_refuse(&db, true), "{name}");
        }
        // 10090 has a common name and a kingdom with bytes >= 0x80 in this taxonomy database,
        // which XML2 does not write: only the scientific name is checked.
        let db = database("sprot_mixtax", Some(&variant("high_bytes_latin1")));
        assert!(!taxdb_may_refuse(&db, true));
    }

    #[test]
    fn a_name_that_cannot_be_written_holds_whether_or_not_the_taxid_is_reported() {
        // sprot_tax: 9606 is the taxid of both sequences. sprot_mixtax: 9606 belongs to the
        // third sequence only, which a search for the first two never reports.
        for name in ["sprot_tax", "sprot_mixtax"] {
            let db = database(name, Some(&variant("high_bytes_utf8")));
            assert!(taxdb_may_refuse(&db, true), "{name}");
            // The prefetch asks for no name: nothing it looks up can be refused.
            assert!(!taxdb_may_refuse(&db, false), "{name}");
        }
    }

    #[test]
    fn a_record_past_the_last_page_holds() {
        // The record of 10090 (entry 3) ends at the offset of entry 4: 5000 bytes, past the
        // 4096-byte page that ends the mapping of the small `taxdb.btd`.
        let dir = std::env::temp_dir().join(format!("losat_prepass_btd_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::copy(variant("good").join("taxdb.btd"), dir.join("taxdb.btd")).unwrap();
        let mut index = std::fs::read(variant("good").join("taxdb.bti")).unwrap();
        let at = 24 + 8 * 4 + 4;
        index[at..at + 4].copy_from_slice(&5000u32.to_be_bytes());
        std::fs::write(dir.join("taxdb.bti"), &index).unwrap();
        for names in [true, false] {
            assert!(taxdb_may_refuse(
                &database("sprot_mixtax", Some(&dir)),
                names
            ));
            // 9606 (entry 2) is fine: the database with that taxid alone does not hold.
            assert!(!taxdb_may_refuse(&database("sprot_tax", Some(&dir)), names));
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_missing_or_damaged_taxdb_that_ncbi_survives_does_not_hold() {
        // No `taxdb.bti` on the search path: the formats write no name.
        let db = database("sprot_mixtax", None);
        assert!(!db.taxdb_found());
        assert!(!taxdb_may_refuse(&db, true));
        // NCBI turns a wrong magic number into a message and N/A (`expected.tsv`).
        let db = database("sprot_mixtax", Some(&variant("bad_magic")));
        assert!(!taxdb_may_refuse(&db, true));
    }

    #[test]
    fn the_two_sides_agree_and_the_taxonomy_side_alone_clears_a_good_taxdb() {
        use crate::seqdb::taxdb::LazyTaxDb;
        for (taxonomy, expected) in [
            ("good", Some(false)),
            ("high_bytes_latin1", Some(false)),
            ("high_bytes_utf8", Some(true)),
            ("offset_next_past_end", Some(true)),
            ("bad_magic", Some(false)),
            // Not checked: whatever the pass says, the answer is the database side's.
            ("unsorted_adjacent", None),
        ] {
            for name in ["sprot", "sprot_tax", "sprot_mixtax"] {
                let db = database(name, Some(&variant(taxonomy)));
                let taxdb = LazyTaxDb::new(db.search_path().unwrap().clone());
                for names in [true, false] {
                    let side = taxdb_side(&taxdb, names);
                    if let Some(side) = side {
                        // The taxonomy side is a superset of the database side.
                        assert!(
                            side || !database_side(&db, &taxdb, names),
                            "{taxonomy} {name}"
                        );
                        if !side {
                            assert!(!taxdb_may_refuse(&db, names), "{taxonomy} {name}");
                        }
                    }
                    // Without the names only the record past the end refuses.
                    if name == "sprot_mixtax" && expected.is_some() {
                        let expected = expected
                            .map(|refused| refused && (names || taxonomy != "high_bytes_utf8"));
                        assert_eq!(side, expected, "{taxonomy} {names}");
                    }
                    // The database side alone and the dispatch give one answer.
                    let alone = LazyTaxDb::new(db.search_path().unwrap().clone());
                    assert_eq!(
                        taxdb_may_refuse(&db, names),
                        database_side(&db, &alone, names),
                        "{taxonomy} {name} {names}"
                    );
                }
            }
        }
    }

    /// Times the taxonomy side on a real NCBI `taxdb` (`LOSAT_TEST_NCBI_TAXDB=<dir with taxdb.bti
    /// and taxdb.btd>`): `cargo test --release`-like timing of one pass, and whether any record
    /// of it could be refused. Ignored: it needs the 190 MB download of OF5b's fix.
    #[test]
    #[ignore]
    fn real_ncbi_taxdb_walk() {
        use crate::seqdb::taxdb::LazyTaxDb;
        let Some(dir) = std::env::var_os("LOSAT_TEST_NCBI_TAXDB") else {
            return;
        };
        let cwd = std::env::temp_dir();
        let search = SearchPath::new(&cwd, Some(&dir), None).unwrap();
        let taxdb = LazyTaxDb::new(search);
        let entries = taxdb.entry_count().unwrap();
        for _ in 0..3 {
            let start = std::time::Instant::now();
            let answer = taxdb_side(&taxdb, true);
            eprintln!(
                "real taxdb: {entries} entries, taxonomy side {answer:?} in {:?}",
                start.elapsed()
            );
            assert_eq!(answer, Some(false));
        }
    }

    /// `prefetch_messages`: the open texts at the run's first lookup, a volume's lookup once per
    /// taxid over the batches of a run where no volume holds another taxid of the slot of NCBI's
    /// cache that the offset error's taxid takes, nothing without a taxonomy database or a taxid.
    #[test]
    fn the_prefetch_writes_each_text_once_without_a_cache_collision() {
        let text = |db: &SeqDb, run: &TaxdbRun, ordinals: &[u32]| -> String {
            let subjects = ordinals.iter().map(|&ordinal| (ordinal, None));
            String::from_utf8(prefetch_messages(db, run, subjects, "BLASTP").unwrap()).unwrap()
        };
        let db = database("sprot_mixtax", Some(&variant("bad_magic")));
        let run = TaxdbRun::default();
        assert_eq!(
            text(&db, &run, &[0, 1, 2, 0]),
            "Error: [blastp] Error: Tax database file has wrong magic number.\n"
        );
        assert_eq!(text(&db, &run, &[2, 1]), "");
        // No taxonomy database on the search path, or one that writes nothing.
        for taxdb in [None, Some(variant("good"))] {
            let db = database("sprot_mixtax", taxdb.as_deref());
            assert_eq!(text(&db, &TaxdbRun::default(), &[0, 1, 2]), "");
        }
        // A database without taxids does not open the taxonomy database.
        let db = database("sprot", Some(&variant("bad_magic")));
        assert_eq!(text(&db, &TaxdbRun::default(), &[0, 1]), "");
        let cwd = empty_cwd("prefetch");
        // The last entry's offset error (taxid 99990 of x16_taxid_99990), once per volume; the
        // two volumes of tm_2vol (copies of x16_taxid_99990) write it once each.
        let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
        for (name, volumes) in [("x16_taxid_99990", 1), ("tm_2vol", 2)] {
            let db = db_search(name, Some(&variant("offset_last_past_end")), &cwd);
            let run = TaxdbRun::default();
            let all: Vec<u32> = (0..db.num_seqs()).collect();
            assert_eq!(text(&db, &run, &all), offset.repeat(volumes), "{name}");
            assert_eq!(text(&db, &run, &all), "", "{name}");
        }
        // tm_cache151 (parsed Seq-ids): 99990 (st1, st3) and 151 (st2), which takes another slot
        // than 99990's (99990 & 255 = 150): once, whatever the order.
        let db = db_search("tm_cache151", Some(&variant("offset_last_past_end")), &cwd);
        let parsed = parsed_subjects(&db);
        for order in [[0, 1, 2], [2, 1, 0], [1, 0, 2]] {
            let run = TaxdbRun::default();
            for expected in [offset, ""] {
                let subjects = order
                    .iter()
                    .map(|&ordinal| (ordinal, Some(&parsed[ordinal as usize])));
                let text = prefetch_messages(&db, &run, subjects, "BLASTP").unwrap();
                assert_eq!(String::from_utf8(text).unwrap(), expected, "{order:?}");
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    /// `prefetch_messages`: a run whose volume holds the offset error's taxid (99990, slot 150 of
    /// NCBI's 256) and another taxid of that slot (150) is refused at its first prefetch, before
    /// its subjects are looked at, and at every later one; the taxonomy databases that write no
    /// offset error, and the same database with 151, are not.
    #[test]
    fn a_cache_collision_with_the_offset_error_is_refused_at_the_first_prefetch() {
        // The databases have parsed Seq-ids: the subjects with their definition lines.
        let prefetch = |db: &SeqDb, run: &TaxdbRun, ordinals: &[u32]| {
            let parsed = parsed_subjects(db);
            let subjects = ordinals
                .iter()
                .map(|&ordinal| (ordinal, Some(&parsed[ordinal as usize])))
                .collect::<Vec<_>>();
            prefetch_messages(db, run, subjects, "BLASTP")
        };
        let cwd = empty_cwd("cache_collision");
        let db = db_search("tm_cache150", Some(&variant("offset_last_past_end")), &cwd);
        let run = TaxdbRun::default();
        assert!(run.prefetch_refused(&db));
        for ordinals in [&[][..], &[0], &[1, 2]] {
            let error = prefetch(&db, &run, ordinals).unwrap_err();
            assert!(error.is::<TaxonomyRefusal>(), "{error:#}");
            let text = format!("{error:#}");
            assert!(text.contains("taxid 99990"), "{text}");
            assert!(text.contains("taxid 150,"), "{text}");
            assert!(text.contains("256-slot taxonomy cache"), "{text}");
            assert!(text.contains("LOSAT's BLASTP"), "{text}");
        }
        // Nothing was looked up: the run's own handle is not open.
        assert!(run.lazy.lock().unwrap().is_none());
        for taxdb in ["good", "bad_magic", "high_bytes_utf8"] {
            let db = db_search("tm_cache150", Some(&variant(taxdb)), &cwd);
            let run = TaxdbRun::default();
            assert!(!run.prefetch_refused(&db), "{taxdb}");
            assert!(prefetch(&db, &run, &[0, 1, 2]).is_ok(), "{taxdb}");
        }
        let db = db_search("tm_cache150", None, &cwd);
        assert!(!TaxdbRun::default().prefetch_refused(&db));
        let db = db_search("tm_cache151", Some(&variant("offset_last_past_end")), &cwd);
        assert!(!TaxdbRun::default().prefetch_refused(&db));
        // The two volumes of tm_2vol hold 99990 alone.
        let db = db_search("tm_2vol", Some(&variant("offset_last_past_end")), &cwd);
        assert!(!TaxdbRun::default().prefetch_refused(&db));
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn a_database_without_taxids_does_not_hold() {
        let db = database("sprot", Some(&variant("high_bytes_utf8")));
        assert!(db.taxdb_found());
        assert!(!taxdb_may_refuse(&db, true));
    }
}
