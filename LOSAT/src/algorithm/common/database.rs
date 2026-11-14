//! What the search programs share in a search of an NCBI BLAST database (`-db NAME`, read
//! through `crate::seqdb`): the opening of the database by its name, the reader's errors
//! as the command line reports them, the `Database:` lines of the pairwise report, and the
//! rules for the stored titles of the reported sequences. BLASTP (stage DB-3), BLASTN
//! (stage DB-4) and TBLASTN (stage DB-5) use them with their own name in LOSAT's refusals;
//! TBLASTN and TBLASTX check a nucleotide database as a copy of it into packed subjects would
//! be read (`x_check_nucleotide_stores`; the copy, `read_nucleotide_stores`, is the reference of
//! the unit tests) and read each sequence where they search it (`x_read_nucleotide_store`,
//! `common::x_db_stream`).
//!
//! The reader (`crate::seqdb`) knows nothing of the command line; this module maps its
//! errors to NCBI's application errors (`crate::cli::NativeError`) and to LOSAT's
//! refusals.

use anyhow::{bail, Result};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::pairwise::{DatabaseInfo, PairwiseDatabase};
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
    // The errors of the alias files' counts depend on the list, and those of the headers and the
    // accession tables come after them (`restrict_database`).
    SeqDb::open_search_in(name, seqtype, &path).map_err(DatabaseOpenError::Reader)
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:109-139
// ```c++
//         if (gi_list || neg_list || m_Aliases.NeedTotalsScan(m_VolSet)) {
//             m_NeedTotalsScan = true;
//             x_InitIdSet();
//         }
//         ...
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
//             ...
//         }
//         m_NumSeqsStats     = x_GetNumSeqsStats();
//         m_TotalLengthStats = x_GetTotalLengthStats();
// ```
// The list is resolved in the totals scan (`x_ScanTotals`, `CSeqDBGiListSet`: its errors and
// warnings), before the alias files' statistics keys are read.
/// Restricts the search's database `db` (opened as `name` by `open_database`) to the taxid or
/// seqid list of the command line (line DI-8, `DbFilterArgs::resolve`; none: the whole
/// database), as NCBI's `CSeqDB` opened with the list (`SeqDb::restrict`), then reads what its
/// alias files' counts and the reports need: NCBI's errors are `BLAST Database error: <text>`
/// with exit status 2 (`database_error`), and the warnings NCBI writes go to `diagnostics`
/// (`LOG_POST(Warning << ...)`: the text alone; `ERR_POST`: `Warning: [<program>] <text>`,
/// the program in lower case as `SetDiagPostPrefix` names it; an `ERR_POST` of an error
/// before a throw: `Error: [<program>] <text>`), before an error too.
///
/// The one owner of the order of a search database's errors, which is NCBI's (probes p4c_e*,
/// p4c_f2 of line RG D1a-2): the alias files' tree and ranges (`open_database`), the list, the
/// alias files' counts (`SeqDb::check_alias_values`), the formatter's reopening of each name
/// (`SeqDb::read_name_records`), then what NCBI meets only when it reads the sequences' headers
/// and reports: LOSAT's refusal of mixed accession tables (`SeqDb::check_accession_tables`) and
/// the first headers (`check_seqids`, whose answer it returns: whether the database has parsed
/// Seq-ids).
pub(crate) fn restrict_database(
    db: &mut SeqDb,
    name: &str,
    list: Option<&crate::seqdb::DbList>,
    diagnostics: &mut dyn std::io::Write,
    settings: &ApplicationSettings,
    program: &str,
) -> Result<bool> {
    if let Some(list) = list {
        apply_list(db, list, diagnostics, program)?;
    }
    // NCBI's `CSeqDB` of the search meets the errors of the alias files' counts after the
    // list (`SeqDb::check_alias_values`); its formatter then opens each name of the list again
    // without it, before the search (`pairwise_database`), and meets the others there.
    db.check_alias_values()
        .map_err(|error| database_error(error, program))?;
    if db.is_restricted() && db.statistics_per_subject() {
        return Err(database_error(
            SeqDbError::Unsupported {
                path: std::path::PathBuf::from(name),
                message: format!(
                    "the total length of NCBI's scan of the list is 0 while the search reads {} \
                     sequences (the scan stops at the FIRST_OID-LAST_OID range of the only alias \
                     file of the -db), so NCBI computes the statistics of each subject on its \
                     own; LOSAT does not search such a database with this list yet",
                    db.num_seqs()
                ),
            },
            program,
        ));
    }
    db.read_name_records(name)
        .map_err(|error| database_error(error, program))?;
    db.check_accession_tables()
        .map_err(|error| database_error(error, program))?;
    // A database with parsed Seq-ids is read (line DBID), unless a variable changes NCBI's choice
    // of its Seq-ids.
    check_seqids(db, settings, program)
}

/// `restrict_database` with a list.
fn apply_list(
    db: &mut SeqDb,
    list: &crate::seqdb::DbList,
    diagnostics: &mut dyn std::io::Write,
    program: &str,
) -> Result<()> {
    let mut warnings = Vec::new();
    let restricted = db.restrict(list, &mut warnings);
    let prefix = program.to_ascii_lowercase();
    for warning in warnings {
        let line = match warning {
            crate::seqdb::ListWarning::Log(text) => text,
            crate::seqdb::ListWarning::Posted(text) => format!("Warning: [{prefix}] {text}"),
            crate::seqdb::ListWarning::Error(text) => format!("Error: [{prefix}] {text}"),
        };
        diagnostics.write_all(line.as_bytes())?;
        diagnostics.write_all(b"\n")?;
    }
    restricted.map_err(|error| database_error(error, program))
}

/// Tells the BLAST archive of the run (`archive`, `None` without one) and its search strategy
/// export (`strategy`) what the database's list makes of their requests: the run's options
/// handle, made once per run and shared by both (`report::archive::OptionsHandle::set_database_filter`).
/// `empty_length` is whether the
/// total length of the statistics is not positive (`SeqDb::statistics_length_is_empty`: a list
/// that kept no sequence, without `STATS_TOTLEN`). `dbsize` is whether `-dbsize` is given: the
/// `EffectiveSearchSpace` option then counts the length adjustment with no sequence, while the
/// searches read the length of OID 0 as one sequence (probes E14 and E17 of line DI-8): LOSAT
/// refuses it (the strategy's requests after the first carry that option too).
pub(crate) fn note_archive_database_filter(
    archive: Option<&crate::report::archive::ArchiveHandle>,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    filter: &crate::blastinput::db_filter::DbFilterArgs,
    list_taxids: &[i32],
    empty_length: bool,
    searchsp: bool,
    dbsize: bool,
    program: &str,
) -> Result<()> {
    if archive.is_none() && strategy.is_none() {
        return Ok(());
    }
    if empty_length && dbsize {
        let option = if archive.is_some() {
            "-outfmt 11"
        } else {
            "-export_search_strategy"
        };
        anyhow::bail!(
            "-dbsize with {option} and a list that keeps no sequence of the database is not supported by LOSAT's {program}"
        );
    }
    let options = match (archive, strategy) {
        (Some(handle), strategy) => {
            anyhow::ensure!(
                strategy.is_none_or(|strategy| strategy.options().same(handle.options())),
                "LOSAT internal error: the archive and the search strategy of one run have two options handles"
            );
            handle.options()
        }
        (None, Some(strategy)) => strategy.options(),
        (None, None) => return Ok(()),
    };
    options.set_database_filter(
        filter.archive_parameter(list_taxids),
        empty_length && !searchsp,
    )?;
    Ok(())
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
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    if let Some(exception) = error.ncbi_exception() {
        return NativeError {
            exit: 255,
            message: format!("Error: {}\n", exception.what()),
        }
        .into();
    }
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

/// What LOSAT refuses where the statistics read 0 sequences of a positive total length
/// (`check_statistics_count`, and BLASTX's own check of its statistics).
pub(crate) const STATISTICS_COUNT_ZERO: &str = "the statistics count 0 sequences of the database (an alias file's NSEQ 0, or its LENGTH or STATS_TOTLEN over volumes without sequences), by which NCBI BLAST+ divides its length";

/// LOSAT's refusal of statistics that read 0 sequences of a positive total length
/// (`SeqDb::statistics_count_is_zero`), for a search of `program` without the Gumbel block:
/// NCBI BLAST+ 2.17.0 divides the total by 0 there and crashes (DECISIONS 14 of line RG; an
/// alias file's `LENGTH` 100 alone over a volume of no OID, probe p4c_dn_nB of D1a-2: SIGFPE).
pub(crate) fn check_statistics_count(db: &SeqDb, program: &str) -> Result<()> {
    if db.statistics_count_is_zero() {
        return Err(refusal(STATISTICS_COUNT_ZERO, program));
    }
    Ok(())
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
/// is that search's error (`search_read_error`). The searches read their sequences where they
/// search them (`common::x_db_stream`); this copy is the reference of the unit tests.
#[cfg(test)]
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

/// What `x_check_nucleotide_stores` keeps of the lengths it reads: the number of sequences and
/// the lengths of the first, of the longest and of all of them, not a length per sequence (a
/// `u32` each was 529 MB on core_nt, held through the search).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct XChecked {
    pub(crate) count: usize,
    pub(crate) first: u32,
    pub(crate) longest: u32,
    pub(crate) total: u64,
}

impl XChecked {
    fn add(&mut self, length: u32) {
        if self.count == 0 {
            self.first = length;
        }
        self.count += 1;
        self.longest = self.longest.max(length);
        self.total += u64::from(length);
    }
}

/// The checks of `read_nucleotide_stores` without the copy (the `-db` searches of TBLASTN and
/// TBLASTX, which read their subjects where they search them, `x_read_nucleotide_store`). The
/// database is read in the same blocks, and every sequence's length and ambiguity runs are read
/// with the same errors in the same order (a length that does not fit a signed 32-bit integer is
/// an error of the read); the lengths are summed up (`XChecked`), not kept. The search
/// fails the query batch that the copy would fail, with the same text. A block of one sequence (one longer than
/// a block, or the last of its volume) is checked without its 2-bit data
/// (`SeqDb::x_nucleotide_check`): a human chromosome is not read here only to be dropped.
pub(crate) fn x_check_nucleotide_stores(db: &SeqDb, program: &str) -> Result<XChecked> {
    let mut checked = XChecked::default();
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
            checked.add(length as u32);
            next_ordinal = range.end_ordinal();
            continue;
        }
        let block = db
            .read_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let (length, _) = stored_parts(&block, i, range.first_ordinal + i as u32, program)?;
            checked.add(length as u32);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(checked)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:312-322
// ```c++
//     int whole_bytes = int(end_offset - start_offset - 1);
//
//     // The last byte is partially full; the last two bits of
//     // the last byte store the number of nucleotides in the
//     // last byte (0 to 3).
//
//     char amb_char = 0;
//     m_Seq->ReadBytes(& amb_char, end_offset - 1, end_offset);
//
//     int remainder = amb_char & 3;
//     return (whole_bytes * 4) + remainder;
// ```
/// The length of nucleotide `ordinal` of `db`, read from the index and the last byte of its 2-bit
/// data (`SeqDb::nucleotide_length`), as NCBI reads the length of a sequence it does not fetch:
/// the length of the sequence's store, with the error of its read (`search_read_error`).
pub(crate) fn x_nucleotide_length(db: &SeqDb, ordinal: u32, program: &str) -> Result<usize> {
    db.nucleotide_length(ordinal)
        .map_err(|error| search_read_error(error, program))
}

/// Sequence `ordinal` of `db` as `read_nucleotide_stores` packs
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

/// The length and the ambiguity runs (as stored) of sequence `i` of `block` (ordinal id
/// `ordinal`), with the errors of the sequence's store (`x_read_nucleotide_store`), without the
/// store.
pub(crate) fn x_nucleotide_parts_in(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    ordinal: u32,
    program: &str,
) -> Result<(usize, Vec<crate::seqdb::AmbiguityRun>)> {
    stored_parts(block, i, ordinal, program)
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
/// read (`stored_title`): `query` is the query whose hit list reached it. NCBI fails that
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
/// What the `Database:` lines of the pairwise report show of the database `db`: for each
/// name of its `-db` list (`SeqDb::name_records`), that database's title (its name for an
/// empty title), date and real totals, also with `-dbsize` or a list.
pub(crate) fn pairwise_database(db: &SeqDb) -> PairwiseDatabase {
    PairwiseDatabase::new(
        db.name_records()
            .iter()
            .map(|record| DatabaseInfo {
                definition: if record.title.is_empty() {
                    record.name.as_bytes().to_vec()
                } else {
                    record.title.clone()
                },
                date: record.date.clone(),
                num_sequences: u64::from(record.num_seqs),
                total_letters: record.total_length,
            })
            .collect(),
    )
}

/// The stored title of a database sequence, as the reports read it (`stored_title`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct StoredTitle {
    /// The title of the sequence's definition line, with a byte outside ' ' ..= '~' read
    /// as '#'.
    pub(crate) title: Vec<u8>,
    /// Whether the stored title held a byte outside ' ' ..= '~' (`check_reported_title`).
    pub(crate) replaced_bytes: bool,
    /// Whether NCBI makes the title that outfmt 0 shows from the source of the sequence's
    /// Bioseq, where the taxonomy database names the taxid of its definition line
    /// (`HeaderReader::stored_taxonomy_title`, `check_shown_stored_title`).
    pub(crate) taxonomy_title: bool,
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
/// hit list of query `query` reaches it; `deflines` are its definition lines as read
/// (`SeqDb::deflines`, or `SeqDb::visit_deflines` with the same results). A reader error of
/// NCBI's own is a `HeaderReadError` that names the query (NCBI fails that query's batch);
/// LOSAT's refusals (`program` names the program) are errors of their own.
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
        taxonomy_title: false,
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
    /// A sequence of a database without parsed Seq-ids: its one title (`stored_title`).
    Stored(StoredTitle),
    /// A sequence with parsed Seq-ids (`report::db_subject`).
    Parsed(std::sync::Arc<crate::report::db_subject::ParsedSubject>),
}

#[cfg(all(test, unix))]
impl ReportedHeader {
    /// The stored title of a sequence without parsed Seq-ids.
    pub(crate) fn stored(&self) -> &StoredTitle {
        match self {
            ReportedHeader::Stored(stored) => stored,
            ReportedHeader::Parsed(_) => panic!("a sequence with parsed Seq-ids"),
        }
    }
}

/// The reader of the headers of a search's reported database sequences (`read_header`): the
/// database, with NCBI's lookup of the Seq-align ids of those with parsed Seq-ids
/// (`SeqDb::seqid_lookup`), made at the first of them and kept for the reader's other reads.
pub(crate) struct HeaderReader<'a> {
    db: &'a SeqDb,
    lookup: std::cell::OnceCell<crate::seqdb::lookup::SeqidLookup<'a>>,
    /// Whether NCBI's formatters may read some reported sequences one by one
    /// (`prefetch_may_be_off`), where the weak matches of its scope are looked for
    /// (`weakly_held`), which the reports refuse where a batch that reports the sequence has
    /// its prefetch off (`ParsedShown::one_by_one`).
    prefetch_off: bool,
    /// The taxonomy database of the search path where `taxdb.bti` is found and it opens, opened
    /// at the first record whose title NCBI can make from its taxids (`source_named`).
    taxdb: std::cell::OnceCell<Option<crate::seqdb::taxdb::TaxDb>>,
    /// `source_named` of the taxids of empty stored titles (`stored_taxonomy_title`), by taxid.
    stored_named: std::cell::RefCell<std::collections::HashMap<i32, bool>>,
    /// `LookupTexts::of`, asked at the first record found with other taxids or in another
    /// volume (`copy_of`).
    lookup_texts: std::cell::OnceCell<LookupTexts>,
}

/// How a sequence of a database differs from another as the reports show them
/// (`HeaderReader::copy_of`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RecordCopy {
    /// NCBI's formatters show the one as the other.
    Same,
    /// The same but where a report writes the taxids or their names.
    OtherTaxids,
    /// The same definition lines but for their taxids, and letters, that a report can still
    /// show apart (`FoundInstead::CopyApart`).
    Apart(crate::report::db_subject::CopyApart),
    /// A report can show the one apart.
    Differs,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:216-220
// ```c++
// bool CSeqDBTaxInfo::GetTaxNames(TTaxId           tax_id,
//                                 SSeqDBTaxInfo  & info )
// {
// 	static CTaxDBFileInfo t;
//     if (t.IsMissingTaxInfo()) return false;
// ```
/// What the lookups of the taxonomy database of the search path of a database write
/// (`TaxdbRun`): the texts of opening it, once at the run's first lookup (the static
/// `CTaxDBFileInfo`), and the offset error of its last entry at each volume's lookup of that
/// taxid. Nothing without `taxdb.bti` on the path.
#[derive(Clone, Copy, Debug, Default)]
struct LookupTexts {
    /// The database writes texts of opening.
    open: bool,
    /// The taxid whose lookup writes the offset error (`TaxDb::offset_error_taxid`).
    offset_error_taxid: Option<i32>,
    /// The database does not open or read: taken to write a text at any lookup.
    unknown: bool,
}

impl LookupTexts {
    /// The texts of the lookups of the taxonomy database of `db`.
    fn of(db: &SeqDb) -> Self {
        if !db.taxdb_found() {
            return Self::default();
        }
        let Some(search_path) = db.search_path() else {
            return Self::default();
        };
        let unknown = Self {
            unknown: true,
            ..Self::default()
        };
        let Ok(taxdb) = crate::seqdb::taxdb::TaxDb::open(search_path) else {
            return unknown;
        };
        let Ok(offset_error_taxid) = taxdb.offset_error_taxid() else {
            return unknown;
        };
        Self {
            open: !taxdb.open_messages().is_empty(),
            offset_error_taxid,
            unknown: false,
        }
    }

    /// Whether the texts depend on which of two copies NCBI fetches, the one with the taxids
    /// `taxids` looked up by its volume or the other with `other` by another volume (each
    /// `ParsedSubject::volume_taxids`): the offset error where either holds its taxid, the texts
    /// of opening where only one of them looks a taxid up (the run's first lookup can move).
    fn differ(self, taxids: &[i32], other: &[i32]) -> bool {
        (self.open && taxids.is_empty() != other.is_empty())
            || self.per_volume(taxids)
            || self.per_volume(other)
    }

    /// Whether a volume's lookups of the taxids `taxids` can write a text that another volume's
    /// lookups of them write again: the offset error.
    fn per_volume(self, taxids: &[i32]) -> bool {
        !taxids.is_empty()
            && (self.unknown
                || self
                    .offset_error_taxid
                    .is_some_and(|taxid| taxids.contains(&taxid)))
    }
}

/// Whether the lookups of the taxonomy database of `db` can write other texts for two of the
/// taxid lists `choices`, each what a prefetch can read for one subject (a sequence's positive
/// taxids, `ParsedSubject::volume_taxids`; empty for reading nothing): `LookupTexts::differ`.
/// `false` without `taxdb.bti` on the search path and for a taxonomy database whose lookups
/// write nothing.
pub(crate) fn lookup_texts_differ(db: &SeqDb, choices: &[&[i32]]) -> bool {
    let texts = LookupTexts::of(db);
    choices.iter().enumerate().any(|(index, taxids)| {
        choices[index + 1..]
            .iter()
            .any(|other| texts.differ(taxids, other))
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:728-753
// ```c++
// static bool
// s_PreFetchSeqs(const blast::CSearchResultSet& results,
// 		       blast::CFormattingArgs::EOutputFormat format_type)
// {
// 	{
// 		char * pre_fetch_limit_str = getenv("PRE_FETCH_SEQS_LIMIT");
// 		if (pre_fetch_limit_str) {
// 			int pre_fetch_limit = NStr::StringToInt(pre_fetch_limit_str);
// 			if(pre_fetch_limit == 0) {
// 				return false;
// 			}
// 			if(pre_fetch_limit == INT_MAX){
// 				return true;
// 			}
// 			int num_of_seqs = 0;
// 			for(unsigned int i=0; i < results.GetNumResults(); i++) {
// 				if(results[i].HasAlignments()) {
// 					num_of_seqs += results[i].GetSeqAlign()->Size();
// 				}
// 			}
// 			if(num_of_seqs > pre_fetch_limit) {
// 				return false;
// 			}
// 		}
// 	}
//
// 	return s_IsPrefetchFormat(format_type);
// }
// ```
/// Whether NCBI's formatters may read a batch's reported sequences one by one, not all at once
/// before its report: `PRE_FETCH_SEQS_LIMIT` is an integer other than `INT_MAX` (0, or one below
/// the batch's number of Seq-aligns turns NCBI's prefetch off, [`prefetch_off`]; the readers of
/// the headers do not count them, the reports of each batch do, [`BatchPrefetch`]). A value
/// that is no integer stops NCBI before its report; LOSAT refuses it in
/// `blastinput::app::check_search_environment` (not for outfmt 11, whose archive NCBI writes
/// without the prefetch).
pub(crate) fn prefetch_may_be_off() -> bool {
    prefetch_seqs_limit().is_some_and(|limit| limit != i32::MAX)
}

/// Whether NCBI's prefetch is off for a query batch whose results hold `seq_aligns` Seq-aligns
/// (one for each HSP of a gapped search, the count of outfmt 7's `# N hits found`, one for each
/// subject of an ungapped one: [`subject_seq_aligns`]):
/// `PRE_FETCH_SEQS_LIMIT` is an integer other than `INT_MAX` that is 0 or below `seq_aligns`.
/// Its formatters then fetch the sequences that each query's report shows, after the query's
/// own warnings, where the report reads them (`TabularSpec::subject_fetch`).
pub(crate) fn prefetch_off(seq_aligns: usize) -> bool {
    prefetch_off_at(prefetch_seqs_limit(), seq_aligns)
}

/// [`prefetch_off`] with the limit `limit` (`prefetch_seqs_limit`).
fn prefetch_off_at(limit: Option<i32>, seq_aligns: usize) -> bool {
    limit.is_some_and(|limit| {
        limit != i32::MAX
            && (limit == 0 || i64::try_from(seq_aligns).unwrap_or(i64::MAX) > i64::from(limit))
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1417-1446
// ```c++
//     /* All HSPs are put in one seqalign, containing a list of
//      * DenseDiag for same molecule search, or StdSeg for translated searches.
//      */
//     ...
//     sa_vector.push_back(seqalign);
// ```
/// The number of Seq-aligns NCBI makes of a subject's `hsps` HSPs, which [`prefetch_off`]
/// counts: one for each HSP of a gapped search (`BLASTHspListToSeqAlign`), one for all of them
/// in an ungapped search (`BLASTUngappedHspListToSeqAlign`; TBLASTX's is always ungapped).
pub(crate) fn subject_seq_aligns(hsps: usize, gapped: bool) -> usize {
    if gapped {
        hsps
    } else {
        usize::from(hsps > 0)
    }
}

/// Whether NCBI's prefetch is off for every query batch of the run: `PRE_FETCH_SEQS_LIMIT` is 0 or
/// negative ([`prefetch_off`] of a batch without Seq-aligns). A run whose reports then read no
/// subject looks no taxid up, so it has no refusal of NCBI's taxonomy cache to ask for.
pub(crate) fn prefetch_always_off() -> bool {
    prefetch_off(0)
}

/// Whether the prefetch of each query's batch reads the batch's reported subjects
/// (`BatchPrefetch::reads`): NCBI counts each batch's Seq-aligns against a positive
/// `PRE_FETCH_SEQS_LIMIT` ([`prefetch_off`]), so a batch above it reads nothing of its subjects
/// before its reports, whose formatters then read them one by one. That matters for a call whose
/// reports read nothing of the subjects ([`ParsedShown::prefetched`], refused by
/// [`check_prefetched_parsed`]) and for the weak matches of NCBI's scope, which only such
/// reads make ([`ParsedShown::one_by_one`]).
#[derive(Debug)]
pub(crate) struct BatchPrefetch(Batches);

/// The prefetch of the batches of a call ([`BatchPrefetch`]).
#[derive(Debug)]
enum Batches {
    /// Every batch prefetches (no limit, `INT_MAX`, a subject set, a test).
    Every,
    /// No batch does (0 or below).
    None,
    /// Whether the batch of each query of the call prefetches (a positive limit); `None` for a
    /// query of no batch range, which is taken both ways ([`BatchPrefetch::reads`]).
    Each(Vec<Option<bool>>),
}

impl BatchPrefetch {
    /// Every batch prefetches (no limit that counts, a subject set, a test).
    pub(crate) fn every() -> Self {
        Self(Batches::Every)
    }

    /// No batch prefetches (a limit of 0 or below).
    fn none() -> Self {
        Self(Batches::None)
    }

    /// The batches of a call of `queries` queries (its hit lists, numbered from 0): `batches`,
    /// the ranges of the queries of each batch numbered from `first` (the call's first query;
    /// ranges before it are the batches written before), or the whole call as one batch where
    /// `batches` is empty; `seq_aligns` counts the Seq-aligns of the hit lists of a range of
    /// the call ([`subject_seq_aligns`]). A query of no range is taken both as prefetched and as
    /// read one by one.
    pub(crate) fn of(
        queries: usize,
        first: usize,
        batches: &[std::ops::Range<usize>],
        seq_aligns: impl Fn(std::ops::Range<usize>) -> usize,
    ) -> Self {
        Self::at(prefetch_seqs_limit(), queries, first, batches, seq_aligns)
    }

    /// [`BatchPrefetch::of`] with the limit `limit` (`prefetch_seqs_limit`).
    fn at(
        limit: Option<i32>,
        queries: usize,
        first: usize,
        batches: &[std::ops::Range<usize>],
        seq_aligns: impl Fn(std::ops::Range<usize>) -> usize,
    ) -> Self {
        if limit.is_none_or(|limit| limit == i32::MAX) {
            return Self::every();
        }
        if prefetch_off_at(limit, 0) {
            return Self::none();
        }
        let mut on = vec![None; queries];
        let whole = [first..first + queries];
        let batches = if batches.is_empty() {
            &whole[..]
        } else {
            batches
        };
        for batch in batches {
            let span = batch.start.saturating_sub(first).min(queries)
                ..batch.end.saturating_sub(first).min(queries);
            if !span.is_empty() {
                let prefetches = !prefetch_off_at(limit, seq_aligns(span.clone()));
                on[span].fill(Some(prefetches));
            }
        }
        Self(Batches::Each(on))
    }

    /// Whether the batch of the call's query `query` prefetches; `None` for a query of no batch.
    fn batch_of(&self, query: usize) -> Option<bool> {
        match &self.0 {
            Batches::Every => Some(true),
            Batches::None => Some(false),
            Batches::Each(on) => on.get(query).copied().flatten(),
        }
    }

    /// Whether the batch of the call's query `query` may read the query's subjects one by one:
    /// its prefetch is off, or the query is of no batch (taken both ways, as in
    /// [`BatchPrefetch::reads`]).
    pub(crate) fn reads_one_by_one(&self, query: usize) -> bool {
        self.batch_of(query) != Some(true)
    }

    /// What the batches read of the subjects of `reported`, the (query of the call, ordinal id)
    /// of each reported subject of each query: the subjects that a batch whose prefetch is on
    /// reports, and those that a batch whose prefetch is off reports (`BatchReads::shown`). The
    /// subjects are listed only where a positive limit counts the batches' Seq-aligns.
    pub(crate) fn reads(&self, reported: impl IntoIterator<Item = (usize, u32)>) -> BatchReads {
        let each = match &self.0 {
            Batches::Every => return BatchReads(BatchesRead::Every),
            Batches::None => return BatchReads(BatchesRead::None),
            Batches::Each(_) => {
                let mut on = std::collections::HashSet::new();
                let mut off = std::collections::HashSet::new();
                for (query, ordinal) in reported {
                    let prefetches = self.batch_of(query);
                    if prefetches != Some(false) {
                        on.insert(ordinal);
                    }
                    if prefetches != Some(true) {
                        off.insert(ordinal);
                    }
                }
                BatchesRead::Each { on, off }
            }
        };
        BatchReads(each)
    }
}

/// What the batches of a call read of each subject they report ([`BatchPrefetch::reads`]).
#[derive(Debug)]
pub(crate) struct BatchReads(BatchesRead);

#[derive(Debug)]
enum BatchesRead {
    Every,
    None,
    /// The ordinal ids that a batch whose prefetch is on reports, and those that one whose
    /// prefetch is off reports.
    Each {
        on: std::collections::HashSet<u32>,
        off: std::collections::HashSet<u32>,
    },
}

impl BatchReads {
    /// Whether a batch whose prefetch is on reports the subject `ordinal`.
    pub(crate) fn prefetched(&self, ordinal: u32) -> bool {
        match &self.0 {
            BatchesRead::Every => true,
            BatchesRead::None => false,
            BatchesRead::Each { on, .. } => on.contains(&ordinal),
        }
    }

    /// Whether a batch whose prefetch is off reports the subject `ordinal`, whose reports then
    /// read it one by one ([`ParsedShown::one_by_one`]).
    pub(crate) fn one_by_one(&self, ordinal: u32) -> bool {
        match &self.0 {
            BatchesRead::Every => false,
            BatchesRead::None => true,
            BatchesRead::Each { off, .. } => off.contains(&ordinal),
        }
    }

    /// What the reports `shown` show of the reported subject `ordinal` ([`ParsedShown::of_batches`]).
    pub(crate) fn shown(&self, shown: ParsedShown, ordinal: u32) -> ParsedShown {
        shown.of_batches(self.prefetched(ordinal), self.one_by_one(ordinal))
    }
}

/// `PRE_FETCH_SEQS_LIMIT` as `NStr::StringToInt` reads it; `None` when it is not set or is no
/// integer.
fn prefetch_seqs_limit() -> Option<i32> {
    std::env::var_os("PRE_FETCH_SEQS_LIMIT")
        .and_then(|value| crate::blastinput::query_batch::ncbi_string_to_int(&value))
}

impl<'a> HeaderReader<'a> {
    /// The reader of the headers of `db`.
    pub(crate) fn new(db: &'a SeqDb) -> Self {
        Self {
            db,
            lookup: std::cell::OnceCell::new(),
            prefetch_off: prefetch_may_be_off(),
            taxdb: std::cell::OnceCell::new(),
            stored_named: Default::default(),
            lookup_texts: std::cell::OnceCell::new(),
        }
    }

    /// The taxonomy database of the search path (`taxdb`), opened at the first call.
    fn taxdb(&self) -> Option<&crate::seqdb::taxdb::TaxDb> {
        self.taxdb
            .get_or_init(|| {
                let path = self.db.search_path().filter(|_| self.db.taxdb_found())?;
                crate::seqdb::taxdb::TaxDb::open(path).ok()
            })
            .as_ref()
    }

    /// Whether the taxonomy database gives a scientific name for the first of the taxids
    /// `taxids` that a volume looks up for a Bioseq (`ParsedSubject::volume_taxids`), which NCBI
    /// then puts in the source its titles are made from (`ParsedSubject::with_taxonomy_title`);
    /// also where the lookup fails, which LOSAT cannot tell from a name.
    fn source_named(&self, taxids: &[i32]) -> bool {
        let Some(&taxid) = taxids.first() else {
            return false;
        };
        let Some(taxdb) = self.taxdb() else {
            return false;
        };
        let mut messages = Vec::new();
        taxdb.names(taxid, &mut messages).map_or(true, |names| {
            names.is_some_and(|names| !names.scientific_name.is_empty())
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // ```c++
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
    // The Bioseq of a sequence without parsed Seq-ids has its definition line's title and a
    // source for its taxid as one with them has (`CSeqDBVol::x_GetTaxonomy`), so NCBI's outfmt 0
    // makes an empty title from the source too: `Homo sapiens` for a nucleotide and `unnamed
    // protein product [Homo sapiens]` for a protein with taxid 9606, in the description and the
    // heading, where the taxonomy database names the taxid (review 1 of line DBID N1b L1,
    // `blastdb-ids-20261021/R/titlesrc/review1/r2.tsv`, `E.tax.nne.blastn.t.f0`).
    /// Whether NCBI makes the outfmt 0 title of the sequence `ordinal` without parsed Seq-ids,
    /// whose stored title is `title`, from its Bioseq's source and the taxonomy database's name of
    /// its definition line's taxid (`report::defline::bioseq_title_from_source`, `source_named`).
    /// The taxid is read only for a title that the source replaces, where a taxonomy database
    /// opens, and each taxid is looked up once; one that does not read counts as named, as a
    /// failed lookup does.
    fn stored_taxonomy_title(&self, title: &[u8], ordinal: u32) -> bool {
        if !crate::report::defline::bioseq_title_from_source(title, &Default::default())
            || self.taxdb().is_none()
        {
            return false;
        }
        let Ok(taxids) = self.db.defline_taxids(ordinal) else {
            return true;
        };
        let Some(taxid) = taxids.into_iter().find(|&taxid| taxid > 0) else {
            return false;
        };
        *self
            .stored_named
            .borrow_mut()
            .entry(taxid)
            .or_insert_with(|| self.source_named(&[taxid]))
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:213-221
    // ```c++
    //             // Try to open the BLAST database even for remote searches, as if
    //             // it is available locally, it will be better to fetch the
    //             // sequence data for formatting from this (local) source
    //             CRef<CSeqDB> seqdb = search_db->GetSeqDb();
    //             db_adapter.Reset(new CLocalDbAdapter(*search_db));
    //             dl = RegisterOMDataLoader(seqdb);
    //             scope->AddDataLoader(dl);
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:541-551
    // ```c++
    //     for ( TTSE_BySeqId::const_iterator it = m_TSE_BySeqId.lower_bound(idh);
    //           it != m_TSE_BySeqId.end() && it->first == idh; ++it ) {
    //         if ( !tse || x_IsBetter(idh, *it->second, *tse) ) {
    //             tse = it->second;
    //         }
    //     }
    // ```
    // NCBI's scope finds the Seq-align id first among the Bioseqs it already holds, the reported
    // sequences of the batch (`CDataSource_ScopeInfo::BestResolve` and `x_GetSeqMatch`,
    // scope_info.cpp:832-850), the first one read where several hold it (`x_IsBetter`), without
    // case, and asks the data loader only where none holds it. So where another reported
    // sequence also holds the Seq-align id, NCBI 2.17.0 can write for the subject the names and
    // identities of that sequence, as its prefetch reads them: `gb|X.1|` with `lcl|555` in its
    // second definition line for `lcl|555`, also for the hit of another query, and with `lcl|ABC`
    // for `lcl|abc` (reviews 1 and 2 of line DBID N2b, `blastdb-ids-20261021/R/unfetch_review1`,
    // `unfetch_review2`), and `gb|X.1|` for a later `emb|X.1|` (line DBID N1). Where no sequence it
    // holds has the Seq-align id, it also takes one that holds an id the Seq-align id matches
    // weakly (`x_GetSeqMatch`, scope_info.cpp:865-882; `seqid::id_scope::weak_keys`), which its
    // prefetch never leaves to it but a scope that reads the reported sequences one by one
    // (`PRE_FETCH_SEQS_LIMIT`) does: `sp|P00001.1|NAME_B` read before for a later `sp||NAME_B`
    // (probe of line DBID N2c, `blastdb-ids-20261021/R/weak/p1.out`).
    /// What NCBI's lookup of the Seq-align id of the sequence `ordinal`, `subject`, finds in the
    /// database where that is not the sequence (`ParsedSubject::found_instead`): NCBI's
    /// formatters look it up in the search's own database, list included
    /// (`SeqidLookup::seqid_to_oids`), and its scope can show another sequence that holds an id it
    /// takes for the Seq-align id (`defline_set::scope_keys`) where both are reported. Another
    /// sequence that holds the same definition lines and letters is shown as the subject itself
    /// (the same record in two volumes or two databases of an alias), one that differs only in
    /// their taxids as the subject but by the reports that write taxids
    /// (`FoundInstead::OtherTaxids`, `copy_of`) and, where a title can be made from them, those
    /// of a generated title (`CopyApart::Title`); these come after any other finding, of which
    /// the first is returned. Where NCBI may read the reported sequences one by one
    /// (`prefetch_off`), also, apart, another sequence that holds an id the Seq-align id matches
    /// weakly (`weakly_held`, `ParsedSubject::weak_match`), which the reports refuse where a
    /// batch that reports the subject has its prefetch off (`ParsedShown::one_by_one`), after
    /// what the lookup finds. A lookup that fails is kept for the formats that read the Bioseq;
    /// the others never look the id up. Where the lookup gives another sequence first or none
    /// (`FoundInstead::Sequence`, `Nothing`), the scan for other holders and the weak match go on
    /// (the third value: another sequence that holds a scope key of the Seq-align id, which
    /// psiblast's PSSM read refuses, `ParsedSubject::pssm_read`). A subject whose Bioseq NCBI's
    /// formatters do not find (`unfetchable`) is looked up too, for that read (its reports refuse
    /// it whatever the lookup gives); a giim id never (NCBI's handles never take it again).
    fn found_instead(
        &self,
        subject: &crate::report::db_subject::ParsedSubject,
        ordinal: u32,
    ) -> (
        Option<crate::report::db_subject::FoundInstead>,
        Option<crate::report::db_subject::WeakMatch>,
        Option<u32>,
    ) {
        use crate::report::db_subject::FoundInstead;
        if matches!(subject.alignment_id(), crate::seqid::ParsedSeqId::Giim(_)) {
            return (None, None, None);
        }
        let failed = |error: SeqDbError| {
            (
                Some(FoundInstead::Failed(
                    error.ncbi_message().unwrap_or_else(|| error.to_string()),
                )),
                None,
                None,
            )
        };
        let lookup = match self.lookup.get() {
            Some(lookup) => lookup,
            None => match self.db.seqid_lookup() {
                Ok(lookup) => {
                    let _ = self.lookup.set(lookup);
                    self.lookup.get().expect("set above")
                }
                Err(error) => return failed(error),
            },
        };
        let id = subject.alignment_id();
        let oids = match lookup.seqid_to_oids(id, true) {
            Ok(oids) => oids,
            Err(error) => {
                return match lookup.throws_on_id(id) {
                    Some(message) => (Some(FoundInstead::ThrowsOnId(message)), None, None),
                    None => failed(error),
                }
            }
        };
        let own = self.db.oid(ordinal);
        let first = oids.first().copied();
        // What the lookup gives instead, where the scan of the other holders goes on.
        let mut given = first.is_none().then_some(FoundInstead::Nothing);
        // The first other sequence found that only the reports of taxids show apart, and the
        // first that those of a generated title also show apart, which those of taxids show
        // apart too (`ParsedSubject::fetched`).
        let mut other_taxids = None;
        let mut title_apart = None;
        let mut copy = |other: u32| -> Option<RecordCopy> {
            let copy = self.copy_of(subject, ordinal, other);
            match copy {
                RecordCopy::Same => None,
                RecordCopy::OtherTaxids => {
                    other_taxids.get_or_insert(other);
                    None
                }
                RecordCopy::Apart(crate::report::db_subject::CopyApart::Title) => {
                    title_apart.get_or_insert(other);
                    None
                }
                RecordCopy::Apart(_) | RecordCopy::Differs => Some(copy),
            }
        };
        if let Some(first) = first.filter(|&first| first != own) {
            let other = self.db.ordinal_of(first);
            match copy(other) {
                Some(RecordCopy::Apart(cause)) => {
                    return (Some(FoundInstead::CopyApart(other, cause)), None, None);
                }
                Some(_) => given = Some(FoundInstead::Sequence(other)),
                None => {}
            }
        }
        let folded = match lookup.case_folded_oids(id) {
            Ok(oids) => oids,
            Err(error) => return failed(error),
        };
        let mut keys = Vec::new();
        crate::report::db_subject::scope_keys(subject, &mut keys);
        // The first was compared above. Besides a sequence given instead, any other holder
        // (also a copy of the subject's record) is one the scope can read for it.
        let mut also_held = None;
        for other in oids
            .iter()
            .chain(&folded)
            .filter(|&&oid| oid != own && Some(oid) != first)
            .map(|&oid| self.db.ordinal_of(oid))
        {
            if !self.holds_key(other, &keys, crate::seqid::defline_set::scope_keys) {
                continue;
            }
            if given.is_some() {
                also_held = Some(other);
                break;
            }
            match copy(other) {
                Some(RecordCopy::Apart(cause)) => {
                    return (Some(FoundInstead::CopyApart(other, cause)), None, None);
                }
                Some(_) => return (Some(FoundInstead::AlsoHeld(other)), None, None),
                None => {}
            }
        }
        // A failed lookup of the weak candidates (LOSAT's own, of the name alone) leaves them
        // unknown where the batch reads one by one (review 1 of line DBID weakbatch, L1).
        let weak = if self.prefetch_off {
            self.weakly_held(lookup, subject, ordinal, &oids)
                .unwrap_or_else(|error| {
                    Some(crate::report::db_subject::WeakMatch::Failed(
                        error.ncbi_message().unwrap_or_else(|| error.to_string()),
                    ))
                })
        } else {
            None
        };
        let found = given.or_else(|| {
            title_apart
                .map(|other| {
                    FoundInstead::CopyApart(other, crate::report::db_subject::CopyApart::Title)
                })
                .or(other_taxids.map(FoundInstead::OtherTaxids))
        });
        (found, weak, also_held)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:128-135
    // ```c++
    //     if(seqid->IsPir() || seqid->IsPrf()) {
    //     	SKeyValuePair kv_pir;
    //     	kv_pir.id = seqid->AsFastaString();
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:160-163
    // ```c++
    //     if(seqid->GetTextseq_Id() != NULL) {
    //     	SKeyValuePair kv;
    //     	kv.id = seqid->GetSeqIdString(false);
    // ```
    // A version 5 database keys a text Seq-id as written, in its letter case, with no key for the
    // name of an id with an accession, and a `pir` or `prf` id by its FASTA text; the scope matches
    // names without case and across the alternatives (`CSeq_id_Mapper::GetMatchingHandles`), so the
    // lookup cannot list the holders of a name (review 1 of line DBID N2c,
    // `blastdb-ids-20261021/R/weak/REVIEW_weak_1.md` H1). A version 4 database's ISAM keys are in
    // lower case without the alternative.
    /// For `found_instead` where NCBI's formatters may read the reported sequences one by one:
    /// another sequence of the database, by ordinal id, that holds an id the Seq-align id `id` of
    /// the sequence `ordinal` matches weakly (`seqid::id_scope::weak_keys`), found among the
    /// sequences that NCBI's lookup gives for `id` (`oids`: in both versions an accession without a
    /// version gives its versions, in version 4 a name its holders) and for its name alone. In
    /// version 5 a Seq-align id with a name, or a `pir` or `prf` one, with weak keys is refused
    /// whatever the database holds (`WeakMatch::Unlisted`).
    fn weakly_held(
        &self,
        lookup: &crate::seqdb::lookup::SeqidLookup<'_>,
        subject: &crate::report::db_subject::ParsedSubject,
        ordinal: u32,
        oids: &[u32],
    ) -> Result<Option<crate::report::db_subject::WeakMatch>, SeqDbError> {
        use crate::report::db_subject::WeakMatch;
        let id = subject.alignment_id();
        let Some(text) = id.textseq() else {
            return Ok(None);
        };
        let mut keys = Vec::new();
        crate::seqid::id_scope::weak_keys(id, &mut keys);
        if keys.is_empty() {
            return Ok(None);
        }
        let pir_or_prf = matches!(
            id,
            crate::seqid::model::ParsedSeqId::Pir(_) | crate::seqid::model::ParsedSeqId::Prf(_)
        );
        if self.db.is_version5() && (text.name.is_some() || pir_or_prf) {
            return Ok(Some(WeakMatch::Unlisted));
        }
        let mut found = oids.to_vec();
        let name = match (&text.name, pir_or_prf) {
            (Some(name), _) => Some(name),
            (None, true) => text.accession.as_ref(),
            (None, false) => None,
        };
        if let Some(name) = name {
            let mut name_alone = id.clone();
            *name_alone.textseq_mut().expect("a text Seq-id") = crate::seqid::model::TextseqId {
                name: Some(name.clone()),
                accession: None,
                release: None,
                version: None,
            };
            found.extend(lookup.seqid_to_oids(&name_alone, true)?);
        }
        let own = self.db.oid(ordinal);
        Ok(found
            .iter()
            .filter(|&&oid| oid != own)
            .map(|&oid| self.db.ordinal_of(oid))
            .find(|&other| {
                self.holds_key(other, &keys, crate::seqid::id_scope::kind_keys)
                    && self.copy_of(subject, ordinal, other) == RecordCopy::Differs
            })
            .map(WeakMatch::Held))
    }

    /// Whether a Seq-id of the sequence `ordinal` has one of the keys `keys` that `key_of` gives
    /// (`defline_set::scope_keys`: NCBI's scope can show it for a Seq-align id of those keys;
    /// `id_scope::kind_keys`: for a Seq-align id whose `weak_keys` they are).
    fn holds_key(
        &self,
        ordinal: u32,
        keys: &[u64],
        key_of: fn(&crate::seqid::model::ParsedSeqId, &mut Vec<u64>),
    ) -> bool {
        let Ok(crate::seqdb::SequenceHeader::Parsed(set)) = self.db.read_header(ordinal) else {
            return false;
        };
        let mut held = Vec::new();
        set.deflines
            .iter()
            .flat_map(|defline| &defline.ids)
            .any(|id| {
                held.clear();
                key_of(id, &mut held);
                held.iter().any(|key| keys.contains(key))
            })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1130-1139
    // ```c++
    //     for(TBDLLConstIter iter = dl.begin(); iter != dl.end(); iter ++) {
    //         TTaxId taxid = ZERO_TAX_ID;
    //
    //         if ((*iter)->CanGetTaxid()) {
    //             taxid = (*iter)->GetTaxid();
    //         }
    //         if (taxid <= ZERO_TAX_ID) {
    //             continue;
    //         }
    // ```
    // NCBI's formatters show for the subject the sequence that its lookup or its scope gives
    // (`found_instead`): with the same definition lines and letters, they write the subject's
    // names, and the taxids of that sequence, which its volume looks up in the taxonomy database
    // when it fetches the sequence (`CSeqDBVol::x_GetTaxonomy`, the sources of its Bioseq). Of the
    // same record in two databases with other taxids, NCBI 2.17.0 writes the first database's
    // taxids in both rows (outfmt 6 `staxids`, `sscinames`, 15, 16) and the same bytes as for one
    // taxid in outfmt 0, 5, 6, 7 and 10, but the outfmt 0 title of an empty one made from the
    // source with the taxonomy database (probe of line DBID N1b, `blastdb-ids-20261021/R/copytaxid/p1`).
    /// How the sequence `b` differs from the sequence `a`, `subject`, as NCBI's formatters show
    /// one for the other: the same where both hold the same definition lines, as the reports
    /// read them, and letters; also where only the taxids of the definition lines differ
    /// (`RecordCopy::OtherTaxids`, which the reports of taxids and their names alone show apart).
    /// Shown apart (`RecordCopy::Apart`): by every report that reads the Bioseq where the
    /// taxonomy database's texts depend on which copy NCBI fetches, whose volume looks its taxids
    /// up (`LookupTexts::differ`; for the same taxids in another volume, the offset error,
    /// `LookupTexts::per_volume`), and by the reports that show a generated title where it can
    /// depend on the Bioseq's source (`ParsedSubject::title_reads_source`) and the sources
    /// differ.
    fn copy_of(
        &self,
        subject: &crate::report::db_subject::ParsedSubject,
        a: u32,
        b: u32,
    ) -> RecordCopy {
        let (
            Ok(crate::seqdb::SequenceHeader::Parsed(set_a)),
            Ok(crate::seqdb::SequenceHeader::Parsed(set_b)),
        ) = (self.db.read_header(a), self.db.read_header(b))
        else {
            return RecordCopy::Differs;
        };
        let without_taxids = |set: &crate::seqid::DeflineSet| -> Vec<crate::seqid::BlastDefline> {
            set.deflines
                .iter()
                .map(|defline| crate::seqid::BlastDefline {
                    taxid: None,
                    ..defline.clone()
                })
                .collect()
        };
        let copy = if set_a == set_b {
            RecordCopy::Same
        } else if without_taxids(&set_a) == without_taxids(&set_b) {
            RecordCopy::OtherTaxids
        } else {
            return RecordCopy::Differs;
        };
        let (Ok(block_a), Ok(block_b)) = (self.db.read_sequence(a), self.db.read_sequence(b))
        else {
            return RecordCopy::Differs;
        };
        let protein = self.db.seqtype() == SeqType::Protein;
        let same_letters = if protein {
            block_a.protein(0) == block_b.protein(0)
        } else {
            matches!(
                (block_a.nucleotide_iupac(0), block_b.nucleotide_iupac(0)),
                (Ok(letters_a), Ok(letters_b)) if letters_a == letters_b
            )
        };
        if !same_letters {
            return RecordCopy::Differs;
        }
        use crate::report::db_subject::CopyApart;
        let texts = || *self.lookup_texts.get_or_init(|| LookupTexts::of(self.db));
        let taxids = subject.volume_taxids();
        match copy {
            RecordCopy::Same => {
                let other_volume = !matches!(
                    (self.db.volume_index(a), self.db.volume_index(b)),
                    (Ok(volume_a), Ok(volume_b)) if volume_a == volume_b
                );
                // The same taxids look the same up first: of the texts only the offset error
                // comes per volume.
                if other_volume && texts().per_volume(&taxids) {
                    return RecordCopy::Apart(CopyApart::TaxonomyTexts);
                }
                RecordCopy::Same
            }
            _ => {
                let Ok(copy) = crate::report::db_subject::ParsedSubject::new(set_b) else {
                    return RecordCopy::Differs;
                };
                let other = copy.volume_taxids();
                // The texts first: every report that reads the Bioseq shows them, a generated
                // title only some (`ParsedSubject::generated_title_fetched`).
                if texts().differ(&taxids, &other) {
                    return RecordCopy::Apart(CopyApart::TaxonomyTexts);
                }
                let named = self.source_named(&taxids) || self.source_named(&other);
                if other != taxids && subject.title_reads_source(protein, named) {
                    return RecordCopy::Apart(CopyApart::Title);
                }
                RecordCopy::OtherTaxids
            }
        }
    }
}

/// `stored_title` of a program that also reads databases with parsed Seq-ids: the header of the
/// sequence `ordinal`, read as `stored_title` reads it, or with every Seq-id where the record has
/// parsed ones (`SeqDb::read_header`), with what NCBI's lookup of its Seq-align id finds
/// (`HeaderReader`). A parsed subject that LOSAT does not write as NCBI (`ParsedSubject::new`)
/// is LOSAT's refusal.
pub(crate) fn read_header(
    headers: &HeaderReader<'_>,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<ReportedHeader> {
    reported_header(
        headers,
        headers.db.read_header(ordinal),
        ordinal,
        query,
        program,
    )
}

/// The `read_header` of the header `header` of the sequence `ordinal` of the database of
/// `headers`.
fn reported_header(
    headers: &HeaderReader<'_>,
    header: std::result::Result<crate::seqdb::SequenceHeader, SeqDbError>,
    ordinal: u32,
    query: usize,
    program: &str,
) -> Result<ReportedHeader> {
    let db = headers.db;
    match header {
        Ok(crate::seqdb::SequenceHeader::Parsed(set)) => {
            crate::report::db_subject::ParsedSubject::new(set)
                .map_err(|reason| parsed_refusal(ordinal, &reason, program))
                .map(|subject| {
                    let (found, weak, also_held) = headers.found_instead(&subject, ordinal);
                    ReportedHeader::Parsed(std::sync::Arc::new(
                        subject
                            .with_taxonomy_title(|taxids| headers.source_named(taxids))
                            .showing_gis(db.reports_show_gis())
                            .under_positive_list(db.has_positive_list())
                            .found_instead(found)
                            .weak_match(weak)
                            .also_held(also_held),
                    ))
                })
        }
        Ok(crate::seqdb::SequenceHeader::Ordinal(deflines)) => {
            let mut stored = stored_title(Ok(deflines), ordinal, query, program)?;
            stored.taxonomy_title = headers.stored_taxonomy_title(&stored.title, ordinal);
            Ok(ReportedHeader::Stored(stored))
        }
        Err(error) => Err(header_read_error(error, query, program)),
    }
}

/// With `-parse_deflines`, LOSAT's refusal of a reported database sequence whose Seq-align id is
/// the id of a query of the scope `queries` (`query_batch::QueryScopes::of_query` of a query that
/// reports it; a query whose ids the reports show: `-parse_deflines`, and the PSSM query of
/// `-in_pssm` whose Seq-id is not one local id): NCBI's scope finds the query instead, so its
/// reports show the query's title and letters for the database sequence (`dropin-survey-20261029/DI7/DECISIONS.md` 8). Every report
/// of a query batch asks for its own scope, also of a sequence whose header an earlier batch read.
pub(crate) fn refuse_id_of_a_query(
    subject: &crate::report::db_subject::ParsedSubject,
    ordinal: u32,
    queries: &crate::blastinput::query_batch::QueryScope,
    program: &str,
) -> Result<()> {
    let keys = subject.query_scope_keys();
    if queries.has_key(&keys.own) {
        return Err(parsed_refusal(
            ordinal,
            "has the Seq-id of a query read with -parse_deflines (or of the PSSM query of \
             -in_pssm) in the query batch that reports it (in any with outfmt 5; NCBI BLAST+ shows \
             the query for it)",
            program,
        ));
    }
    if queries.has_key(&keys.weak) {
        return Err(parsed_refusal(
            ordinal,
            &format!(
                "has the Seq-align id {}, which NCBI BLAST+'s scope can match weakly to the id \
                 of a query read with -parse_deflines (or of the PSSM query of -in_pssm) in the \
                 query batch that reports it (in any with outfmt 5): an accession with a version for one without, a name with or \
                 without an accession, a pir or prf accession as a name; NCBI BLAST+ then shows \
                 the query for it",
                String::from_utf8_lossy(&subject.alignment_fasta())
            ),
            program,
        ));
    }
    Ok(())
}

/// `refuse_id_of_a_query` of each reported sequence of `reported` (a searched query, by its index
/// in `scopes`, and a sequence that its report names, by ordinal id) whose header `titles` holds
/// with parsed Seq-ids, with the scope of that query. A program calls it for each query batch
/// after it reads the batch's titles: its cache of the titles keeps those of the batches before.
pub(crate) fn refuse_reported_ids_of_queries(
    titles: &std::collections::BTreeMap<u32, ReportedHeader>,
    reported: impl IntoIterator<Item = (usize, u32)>,
    scopes: &crate::blastinput::query_batch::QueryScopes,
    program: &str,
) -> Result<()> {
    for (query, ordinal) in reported {
        if let Some(ReportedHeader::Parsed(subject)) = titles.get(&ordinal) {
            refuse_id_of_a_query(subject, ordinal, scopes.of_query(query), program)?;
        }
    }
    Ok(())
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

/// The check of a search program's database at the end of its opening (`restrict_database`):
/// the first header of each volume is read, as before line DBID, and only a header with parsed
/// Seq-ids is let through (`SeqDb::has_parsed_seqids`; every search program reads such a
/// database). Returns whether the database has parsed Seq-ids; such a database is LOSAT's refusal
/// where a variable changes NCBI's choice of its Seq-ids
/// (`ApplicationSettings::seq_id_rank_variable`).
fn check_seqids(db: &SeqDb, settings: &ApplicationSettings, program: &str) -> Result<bool> {
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
/// `query`, whose Seq-align ids match (`report::db_subject::joined_in_reports`) and that are not
/// copies of one record (`joins_as_copy`): NCBI's XML, XML2, JSON and outfmt 0 reports show
/// their alignments as one hit, and its tabular `qcovs` and `qcovus` the coverage of both on the
/// first row (`N/A` on the other), which LOSAT does only for copies. `subjects` gives the ordinal
/// id and the parsed subject of each HSP list of the hit list, in its order (`None` for a
/// sequence without parsed Seq-ids).
pub(crate) fn refuse_joined_hits<'a>(
    query: usize,
    subjects: impl IntoIterator<Item = Option<(u32, &'a crate::report::db_subject::ParsedSubject)>>,
    program: &str,
) -> Result<()> {
    let mut previous: Option<(u32, &crate::report::db_subject::ParsedSubject)> = None;
    for subject in subjects {
        if let Some((before, ordinal)) = joined_apart(previous, subject) {
            return Err(anyhow::Error::new(ParsedRefusal(format!(
                "database sequences {before} and {ordinal} follow each other in the hits of query {} and have Seq-align ids that NCBI BLAST+ takes for one sequence (its XML and outfmt 0 reports show their alignments as one hit, its tabular qcovs and qcovus count them as one) but are not copies of one record; this is not supported by LOSAT's {program}",
                query + 1
            ))));
        }
        previous = subject;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1025-1033
// ```c++
//         } else {
//             subid = &((*iter)->GetSeq_id(1));
//             if(is_first_aln || (!is_first_aln && !subid->Match(*previous_id))){
//                 ++num_align;
//             }
//
//             if(num_align > number) {
//                  break;
//             }
// ```
/// The number of HSP lists of a hit list (query `query`; `subjects`, as `refuse_joined_hits`
/// takes them) that its cut to `max_target_seqs` subjects keeps, as NCBI's formatter cuts it
/// (`PruneSeqalign`): it counts following Seq-aligns whose Seq-align ids match as those of one
/// subject. A copy of one record that follows another (`joins_as_copy`) is kept with it; two
/// other parsed subjects that the count reaches with matching Seq-align ids are LOSAT's refusal.
/// psiblast cuts its traceback hit lists so, which are longer than `-max_target_seqs`.
pub(crate) fn joined_cut<'a>(
    query: usize,
    max_target_seqs: usize,
    subjects: impl IntoIterator<Item = Option<(u32, &'a crate::report::db_subject::ParsedSubject)>>,
    program: &str,
) -> Result<usize> {
    let mut previous: Option<(u32, &crate::report::db_subject::ParsedSubject)> = None;
    let mut counted = 0usize;
    let mut kept = 0usize;
    for subject in subjects {
        if let Some((before, ordinal)) = joined_apart(previous, subject) {
            return Err(anyhow::Error::new(ParsedRefusal(format!(
                "database sequences {before} and {ordinal} follow each other in the hits of query {} that the cut to -max_target_seqs {max_target_seqs} counts and have Seq-align ids that NCBI BLAST+ takes for one sequence but are not copies of one record (it counts them as one sequence and keeps one more); this is not supported by LOSAT's {program}",
                query + 1
            ))));
        }
        let joined = matches!(
            (previous, subject),
            (Some((_, a)), Some((_, b))) if crate::report::db_subject::joins_as_copy(a, b)
        );
        if !joined {
            counted += 1;
            if counted > max_target_seqs {
                break;
            }
        }
        kept += 1;
        previous = subject;
    }
    Ok(kept)
}

/// The ordinal ids of `previous` and `subject`, parsed subjects that follow each other in a hit
/// list, where their Seq-align ids match (`report::db_subject::joined_in_reports`) but they are
/// not copies of one record (`joins_as_copy`).
fn joined_apart(
    previous: Option<(u32, &crate::report::db_subject::ParsedSubject)>,
    subject: Option<(u32, &crate::report::db_subject::ParsedSubject)>,
) -> Option<(u32, u32)> {
    use crate::report::db_subject::{joined_in_reports, joins_as_copy};
    match (previous, subject) {
        (Some((before, a)), Some((ordinal, b)))
            if joined_in_reports(a, b) && !joins_as_copy(a, b) =>
        {
            Some((before, ordinal))
        }
        _ => None,
    }
}

/// LOSAT's refusal of a parsed subject whose outfmt 0 title NCBI makes from the source
/// (`ParsedSubject::taxonomy_title`): where the taxonomy database names the taxid NCBI names the
/// organism in it, which LOSAT does not do.
pub(crate) fn refuse_parsed_taxonomy(ordinal: u32, program: &str) -> anyhow::Error {
    anyhow::Error::new(ParsedRefusal(format!(
        "database sequence {ordinal} has an empty title that NCBI BLAST+'s outfmt 0 makes from its taxids and the taxonomy database (which names its taxid); this is not supported by LOSAT's {program}"
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
/// taxids and the taxonomy database (`taxdb_found`, `ParsedSubject::title_from_taxonomy`) or
/// from those of the same record that it can show for the sequence
/// (`ParsedSubject::generated_title_fetched`).
pub(crate) fn check_shown_parsed_subject(
    subject: &crate::report::db_subject::ParsedSubject,
    protein: bool,
    taxdb_found: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    let generated = subject.heading_generates_title();
    if protein {
        let (_, title) = subject.description(protein);
        let heading_reads_past_end = generated
            && subject
                .bioseq_title(protein, false)
                .is_ok_and(|defline| defline.reads_past_end);
        if title.reads_past_end || heading_reads_past_end {
            return Err(parsed_refusal(
                ordinal,
                "has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0)",
                program,
            ));
        }
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
    if subject.description_start().len() > 65 {
        return Err(parsed_refusal(
            ordinal,
            "has an outfmt 0 description id longer than 65 bytes (NCBI BLAST+ 2.17.0 pads the description line without end)",
            program,
        ));
    }
    if subject.pairwise_title_reads_source(protein, taxdb_found) {
        subject
            .generated_title_fetched()
            .map_err(|reason| parsed_refusal(ordinal, &reason, program))?;
    }
    if generated && subject.title_suffix_adjusted(protein) {
        return Err(parsed_refusal(
            ordinal,
            "has a definition line without a title and a title ending in ']' that NCBI BLAST+'s outfmt 0 heading rewrites from the taxids (x_AdjustProteinTitleSuffix)",
            program,
        ));
    }
    if subject.taxonomy_title() {
        return Err(refuse_parsed_taxonomy(ordinal, program));
    }
    Ok(())
}

/// What the reports of a run show of their reported database sequences, for LOSAT's refusals of
/// those with parsed Seq-ids (`check_reported_parsed`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ParsedShown {
    /// The tabular formats whose fields read the sequence found by the subject's Seq-align id
    /// (`report::tabular::reads_subject`: the default fields, the subject's ids, titles,
    /// taxonomy, strand and length, the alignment's letters; `sseqid`,
    /// `ParsedSubject::tabular_ids`). A table of other fields (`6 qseqid evalue`) writes the same
    /// whatever that lookup finds.
    pub(crate) tabular: bool,
    /// Outfmt 0 (`check_shown_parsed_subject`).
    pub(crate) pairwise: bool,
    /// The structured documents, outfmt 5 and 13 to 16, which read the subject's Bioseq for its
    /// ids and title.
    pub(crate) structured: bool,
    /// The formats that show the alignments of following hits whose Seq-align ids match as one
    /// hit, outfmt 0, 5 and 13 to 16, and the tabular formats that count their coverage as one,
    /// with `qcovs` or `qcovus` (`TabularNeeds::joins_subjects`; `refuse_joined_hits`). The
    /// coverage groups the Seq-aligns by their ids without reading a Bioseq
    /// (`CBlastFormatUtil::InsertSubjectScores`).
    pub(crate) joined: bool,
    /// SAM (outfmt 17), which names a subject by an id of its Bioseq (`ParsedSubject::sam_id`).
    pub(crate) sam: bool,
    /// The formats that write the subjects' taxids or their names: outfmt 13 to 16 and the
    /// tabular `staxid`, `staxids` and the names (`ParsedSubject::fetched`).
    pub(crate) taxonomy: bool,
    /// The subject ids that the tabular formats print (`ParsedSubject::check_printed_word`).
    pub(crate) words: crate::report::tabular::PrintedIds,
    /// Whether the batches' prefetch reads the subjects by their Seq-align ids (the first
    /// format's: outfmt 0, 5, 6, 7, 10, 17 and 20, with `PRE_FETCH_SEQS_LIMIT=2147483647` also the
    /// Seq-annots 8, 9 and 12, [`prefetches_subjects`]), unless `PRE_FETCH_SEQS_LIMIT` turns it
    /// off for every batch ([`prefetch_always_off`]). Where no report reads their Bioseqs (the
    /// Seq-annots, a table of fields that read nothing of the subject), what the prefetch
    /// cannot read as NCBI does is refused (`check_prefetched_parsed`) for a subject that a
    /// batch whose prefetch is on reports ([`BatchReads::shown`], [`BatchPrefetch`]:
    /// a batch whose Seq-aligns outnumber a positive limit reads nothing).
    pub(crate) prefetched: bool,
    /// Whether NCBI's formatters may read the subject one by one where the reports read it: a
    /// batch that reports it has its prefetch off ([`BatchPrefetch`]; before the batches are
    /// counted, [`prefetch_may_be_off`]). Only there can NCBI's scope show another sequence read
    /// before for a Seq-align id that matches an id of it weakly
    /// (`ParsedSubject::fetched`, `WeakMatch`): the prefetch looks for the ids
    /// it has resolved, then loads the others by their own ids.
    pub(crate) one_by_one: bool,
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1773-1783
// ```c++
//     CInitGuard init(info.second.m_Bioseq_Info, m_MutexPool, CInitGuard::force);
//     if ( init || info.second.m_Bioseq_Info->NeedsReResolve(m_BioseqChangeCounter) ) {
//         if ( get_flag != CScope::eGetBioseq_Resolved ) {
//             // Resolve only if the flag allows
//             x_ResolveSeq_id(info, get_flag, match);
//         }
//         else {
//             // outdated
//             return null;
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3289-3300
// ```c++
//     for ( size_t i = from; i < from + count; ++i ) {
//         ret[i] = GetBioseqHandle(ids[i], CScope::eGetBioseq_Resolved);
//         if ( !ret[i] ) {
//             match_map[ids[i]];
//         }
//     }
//     if ( match_map.empty() ) {
//         return;
//     }
//     for (CPriority_I it(m_setDataSrc); it; ++it) {
//         it->GetBlobs(match_map);
//     }
// ```
// The prefetch (`LoadSequencesToScope`'s `GetBioseqHandles`) asks the scope only for the ids it
// has resolved, which it gets back as they were resolved, and the data loader for the others by
// their own ids; the data source's weak match (data_source.cpp:1439-1471) looks only among the
// blobs that the loader gives for that id. So the weak matches of
// `CDataSource_ScopeInfo::x_GetSeqMatch` (scope_info.cpp:865-882) come from the formatters' (and
// psiblast's PSSM's) own reads where the batch's prefetch is off. A later batch can get them back
// only where the scope is kept: psiblast's rounds of one query, and outfmt 5. Both read one by one
// only subjects that LOSAT refuses there (a round's reports and PSSM, an off batch's reported
// subjects), so no weak match of a run that LOSAT writes reaches a later batch. Probe w0 of line
// DBID weakbatch (`blastdb-ids-20261021/R/weakbatch/w0/run_w0.sh`): a batch of 2 Seq-aligns under
// limit 2 shows each subject, the next batch of 4 shows `sp|P00001.1|NAME_B` for `sp||NAME_B`;
// outfmt 5 with the batches before kept, and a batch of 252 subjects, which the prefetch resolves
// 100 at a time, show each subject.
impl ParsedShown {
    /// What the reports show of a subject that the batches' prefetch reads only where `read` (a
    /// batch that reports it prefetches) and that the formatters may read one by one where
    /// `one_by_one` (a batch that reports it does not; [`BatchReads::shown`]).
    pub(crate) fn of_batches(self, read: bool, one_by_one: bool) -> Self {
        Self {
            prefetched: self.prefetched && read,
            one_by_one,
            ..self
        }
    }

    /// Whether a report reads the subject's Bioseq, which NCBI's formatters find by its
    /// Seq-align id (`ParsedSubject::fetched`); the Seq-annots, the archive and the tables of
    /// fields that read nothing of the subject write the Seq-align id or nothing of it.
    pub(crate) fn reads_bioseq(self) -> bool {
        self.tabular || self.pairwise || self.structured || self.sam || self.taxonomy
    }
}

/// LOSAT's refusals of the reported database sequence `ordinal` with parsed Seq-ids `subject` in
/// the reports `shown` (`protein` for a protein database; `taxdb_found`: `SeqDb::taxdb_found`),
/// as BLASTX's (`blastx::database::ReportedSubjects::for_batch`): a sequence that NCBI's
/// formatters do not find by its Seq-align id, an outfmt 0 description or heading that LOSAT
/// does not write as NCBI (`check_shown_parsed_subject`), a tabular id that it does not
/// write as NCBI (`ParsedSubject::tabular_ids`), and where only the prefetch reads the subject
/// (`shown.prefetched`), what it cannot read as NCBI does (`check_prefetched_parsed`). `db` is
/// the subject's database.
pub(crate) fn check_reported_parsed(
    db: &SeqDb,
    subject: &crate::report::db_subject::ParsedSubject,
    shown: ParsedShown,
    protein: bool,
    taxdb_found: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if shown.reads_bioseq() {
        subject
            .fetched(shown.taxonomy, shown.one_by_one)
            .map_err(|reason| parsed_refusal(ordinal, &reason, program))?;
    } else if shown.prefetched {
        check_prefetched_parsed(db, subject, ordinal, program)?;
    }
    if shown.pairwise {
        check_shown_parsed_subject(subject, protein, taxdb_found, ordinal, program)?;
    }
    if shown.tabular {
        subject
            .check_tabular_ids(protein)
            .and_then(|()| subject.check_printed_word(protein, shown.words))
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
/// protein database. An id that LOSAT cannot write as NCBI is its refusal where a tabular report
/// `reads` the subject (`report::tabular::reads_subject`), else empty: a table of fields that
/// read nothing of the subject prints none of its ids.
pub(crate) fn tabular_subject_ids(
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    ordinals: &[u32],
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    protein: bool,
    reads: bool,
    program: &str,
) -> Result<Vec<std::sync::Arc<[u8]>>> {
    records
        .iter()
        .zip(ordinals)
        .enumerate()
        .map(|(rank, (record, &ordinal))| {
            let id = match parsed.get(rank).and_then(Option::as_deref) {
                Some(subject) => subject
                    .tabular_saccver(protein)
                    .map_err(|problem| parsed_refusal(ordinal, &format!("has {problem}"), program)),
                None => tabular_subject_id(&record.title, ordinal, protein, program)
                    .map(std::sync::Arc::from),
            };
            match id {
                Err(_) if !reads => Ok(std::sync::Arc::from(&b""[..])),
                id => id,
            }
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
// ```c++
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
// ```c++
//        		const CBioseq_Handle& subject_bh =
//                 scope.GetBioseqHandle(align.GetSeq_id(1));
//             if(setSubjectId) {
//             	SetSubjectId(subject_bh);
//             }
//             m_SubjectLength = subject_bh.GetBioseqLength();
// ```
/// A reported sequence as a caller of the library reads it (`FormatObserver::database_subjects`):
/// the id of the default tabular report (`tabular_subject_ids`) of the sequence `ordinal` with the
/// stored title `title` (with parsed Seq-ids: `parsed`), and `length` letters. A sequence whose
/// tabular id LOSAT refuses shows the first word of its title: only a tabular report refuses
/// that id, not the search.
pub(crate) fn library_subject(
    db: &SeqDb,
    title: &[u8],
    ordinal: u32,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    length: usize,
    protein: bool,
    program: &str,
) -> crate::api::local_blast::ReportedSubject {
    let id = match parsed {
        Some(subject) => subject.tabular_saccver(protein).ok().map(|id| id.to_vec()),
        None => tabular_subject_id(title, ordinal, protein, program).ok(),
    }
    .unwrap_or_else(|| {
        title
            .split(|&byte| byte == b' ')
            .next()
            .unwrap_or_default()
            .to_vec()
    });
    crate::api::local_blast::ReportedSubject {
        oid: db.oid(ordinal),
        id: std::sync::Arc::from(String::from_utf8_lossy(&id)),
        length,
    }
}

/// `library_subject` of each reported sequence of a nucleotide database, by rank: `records` hold
/// the stored titles, `ordinals` and `parsed` as in `tabular_subject_ids`.
pub(crate) fn library_subjects(
    db: &SeqDb,
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    ordinals: &[u32],
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    length: impl Fn(usize) -> usize,
    program: &str,
) -> Vec<crate::api::local_blast::ReportedSubject> {
    records
        .iter()
        .zip(ordinals)
        .enumerate()
        .map(|(rank, (record, &ordinal))| {
            let subject = parsed.get(rank).and_then(Option::as_deref);
            library_subject(
                db,
                &record.title,
                ordinal,
                subject,
                length(rank),
                false,
                program,
            )
        })
        .collect()
}

// No NCBI counterpart: bounds LOSAT's own check before the search (`reports_may_refuse`).
/// The most sequences of a database with parsed Seq-ids whose headers a program's check before
/// the search reads with `reports_may_refuse` (tblastx's `report_stream::check`; blastn's
/// report stream read none before the search since 2026-11-03). The check decodes each header (2.3-5.6 us at one thread with
/// tblastx's 18S, line DBID ID4t fix 2b) and keeps 16 bytes per join key with outfmt 0, so this
/// bounds the pass to a few seconds and some tens of MiB; a larger database (NCBI's `nt` and
/// `core_nt`) keeps its reports after the batch loop, as before the stream took databases with
/// parsed Seq-ids. The volumes' sequences count (`SeqDb::printed_num_seqs`), also under a taxid
/// list: the range walks decode the excluded sequences of every range they read.
pub(crate) const PARSED_CHECK_SEQUENCES_MAX: u32 = 1 << 20;

/// Why a report of some sequence of a database could be refused (`reports_may_refuse`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReportRefusal {
    /// A header that does not read, or a sequence whose stored title, tabular id or parsed
    /// Seq-ids the reports refuse.
    Headers,
    /// Two sequences whose Seq-align ids NCBI takes for one and that are not copies of one record
    /// (`refuse_joined_hits`).
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
/// `check_shown_stored_title` and `check_reported_parsed` with `taxdb_found` where outfmt 0 is
/// shown, as `reported_subjects` asks). Where the reports join following hits (`shown.joined`), no two sequences with parsed
/// Seq-ids may be joined but copies of one record (`refuse_joined_hits`): the sequences that share
/// a key of their Seq-align ids (`db_subject::join_keys`; joined ones always share one) are read
/// again and compared (`joined_in_reports`, `joins_as_copy`). With `prefetch` (the batches' prefetch of a search with `taxdb.bti` on
/// the search path) the taxids that the prefetch looks up must read (a stored title's first one,
/// every one of a sequence with parsed Seq-ids, as `prefetch_messages` reads them) and no lookup
/// of them may be refused (`taxids_may_refuse`). Under `-parse_deflines` no sequence may have the
/// id of a query of `queries` (`refuse_id_of_a_query`; the scope of every query of the run, a
/// superset of each batch's). A reader error is a refusal too. The headers are read in `pool`;
/// nothing is written and no state of the run changes. `None`: no report of any sequence can be
/// refused.
pub(crate) fn reports_may_refuse(
    db: &SeqDb,
    shown: ParsedShown,
    protein: bool,
    taxdb_found: bool,
    prefetch: bool,
    queries: Option<&crate::blastinput::query_batch::QueryScope>,
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
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        fn join(mut self, mut other: Kept) -> std::result::Result<Kept, ReportRefusal> {
            self.keys.append(&mut other.keys);
            self.taxids.append(&mut other.taxids);
            Ok(self)
        }
    }
    let check = |kept: &mut Kept,
                 headers: &HeaderReader<'_>,
                 ordinal: u32,
                 header: std::result::Result<crate::seqdb::SequenceHeader, SeqDbError>|
     -> std::result::Result<(), ReportRefusal> {
        match reported_header(headers, header, ordinal, 0, program) {
            Ok(ReportedHeader::Stored(stored)) => {
                check_reported_title(stored.replaced_bytes, ordinal, program)
                    .map_err(|_| ReportRefusal::Headers)?;
                if shown.tabular {
                    tabular_subject_id(&stored.title, ordinal, protein, program)
                        .map_err(|_| ReportRefusal::Headers)?;
                }
                if shown.pairwise {
                    check_shown_stored_title(stored.taxonomy_title, ordinal, program)
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
                    db,
                    &subject,
                    shown,
                    protein,
                    shown.pairwise && taxdb_found,
                    ordinal,
                    program,
                )
                .map_err(|_| ReportRefusal::Headers)?;
                if let Some(queries) = queries {
                    refuse_id_of_a_query(&subject, ordinal, queries, program)
                        .map_err(|_| ReportRefusal::Headers)?;
                }
                if shown.joined {
                    let mut keys = Vec::new();
                    crate::report::db_subject::join_keys(&subject, &mut keys);
                    kept.keys.extend(keys.into_iter().map(|key| (key, ordinal)));
                }
                if prefetch {
                    match prefetch_reads(db, ordinal, &subject, true, program) {
                        Ok(Ok(Some(blob))) if blob == ordinal => {
                            kept.taxids.extend(subject.volume_taxids())
                        }
                        Ok(Ok(Some(other))) => match db.defline_set(other) {
                            Ok(set) => kept.taxids.extend(defline_set_taxids(&set)),
                            Err(_) => return Err(ReportRefusal::Taxonomy),
                        },
                        Ok(Ok(None)) => {}
                        Ok(Err(_)) | Err(_) => return Err(ReportRefusal::Taxonomy),
                    }
                }
            }
            Err(_) => return Err(ReportRefusal::Headers),
        }
        Ok(())
    };
    let check_range = |mut kept: Kept,
                       range: &crate::seqdb::BlockRange|
     -> std::result::Result<Kept, ReportRefusal> {
        let headers = HeaderReader::new(db);
        let mut refusal = None;
        db.visit_sequence_headers(*range, |ordinal, header| {
            refusal = check(&mut kept, &headers, ordinal, header).err();
            refusal.is_none()
        });
        refusal.map_or(Ok(kept), Err)
    };
    // The database's lookup files, opened once before the ranges look the ids up in parallel
    // (`HeaderReader`); an error comes again at the first lookup.
    let _ = db.seqid_lookup();
    let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
    let all = || -> std::result::Result<Kept, ReportRefusal> {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            use rayon::prelude::*;
            return pool.install(|| {
                ranges
                    .par_iter()
                    .try_fold(Kept::default, check_range)
                    .try_reduce(Kept::default, Kept::join)
            });
        }
        let _ = pool;
        ranges.iter().try_fold(Kept::default(), check_range)
    };
    let mut kept = match all() {
        Ok(kept) => kept,
        Err(refusal) => return Some(refusal),
    };
    kept.keys.sort_unstable();
    let headers = HeaderReader::new(db);
    for shared_key in kept.keys.chunk_by(|a, b| a.0 == b.0) {
        let mut ordinals: Vec<u32> = shared_key.iter().map(|&(_, ordinal)| ordinal).collect();
        ordinals.dedup();
        if ordinals.len() < 2 {
            continue;
        }
        let mut subjects = Vec::with_capacity(ordinals.len());
        for ordinal in ordinals {
            match read_header(&headers, ordinal, 0, program) {
                Ok(ReportedHeader::Parsed(subject)) => subjects.push(subject),
                _ => return Some(ReportRefusal::Headers),
            }
        }
        for (index, a) in subjects.iter().enumerate() {
            if subjects[index + 1..].iter().any(|b| {
                crate::report::db_subject::joined_in_reports(a, b)
                    && !crate::report::db_subject::joins_as_copy(a, b)
            }) {
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
/// without `taxdb.bti` there. One pass over the index usually finds that no lookup of any taxid
/// can fail (`TaxDb::any_lookup_may_fail`); else each one is looked up once. Both use a taxonomy
/// handle of their own (as `taxdb_may_refuse`'s database side does), so the run's handle and
/// messages are as they would be without this call. The prefetch's refusal of NCBI's cache
/// (`TaxdbRun::cache_collision`) comes at the run's first prefetch, before any report.
fn taxids_may_refuse(db: &SeqDb, taxids: &[i32]) -> bool {
    if !db.taxdb_found() || taxids.is_empty() {
        return false;
    }
    let Some(search_path) = db.search_path() else {
        return false;
    };
    let taxdb = crate::seqdb::taxdb::LazyTaxDb::new(search_path.clone());
    if matches!(taxdb.any_lookup_may_fail(), Ok(false)) {
        return false;
    }
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

/// LOSAT's refusal of the outfmt 0 title of the shown database sequence `ordinal` without parsed
/// Seq-ids that NCBI makes from its taxid and the taxonomy database (`taxonomy_title`,
/// `StoredTitle::taxonomy_title`), as for one with them (`refuse_parsed_taxonomy`).
pub(crate) fn check_shown_stored_title(
    taxonomy_title: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if taxonomy_title {
        bail!(
            "database sequence {ordinal} has an empty title that NCBI BLAST+'s outfmt 0 makes from its taxid and the taxonomy database (which names that taxid, or whose lookup of it fails); this is not supported by LOSAT's {program}"
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
// A stored gap character (0, which makeblastdb writes for `-` in a protein FASTA) is a
// residue that scores `BLAST_SCORE_MIN` against every letter (row and column 0), as LOSAT's
// kernels score it: an ungapped extension stops at it and a gapped alignment goes round it
// with a gap. A byte that is no NCBIstdaa residue (28 or more) indexes NCBI's matrix past
// its rows and is rejected.
/// LOSAT's refusal of the protein database sequence `ordinal` (`with_sentinels`: its
/// residues with the byte before and the byte after them in the sequence file) when
/// `program` cannot search it as NCBI does: a sequence not between two 0 bytes, or one
/// with a byte that is no NCBIstdaa residue.
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
    // Every residue must be a byte below BLASTAA_SIZE. The largest byte, taken over all bytes
    // without an early exit, is a loop the compiler vectorizes; the early exit of `all` made it a
    // loop over single bytes (LOSAT-only cost; NCBI checks nothing here).
    const _: () = assert!(crate::utils::matrix::BLASTAA_SIZE <= 256);
    let largest = residues.iter().copied().max().unwrap_or(0);
    if usize::from(largest) < crate::utils::matrix::BLASTAA_SIZE {
        return Ok(());
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

/// Bytes of the sequence file that one read of `x_check_protein_stores` holds at most
/// (`SeqDb::next_block`). No NCBI counterpart: NCBI maps the volume files. The size
/// changes nothing but the number of reads.
const PROTEIN_STORE_BLOCK_BYTES: u64 = 1 << 20;

/// Reads every sequence of the protein database `db` once, in blocks of
/// `PROTEIN_STORE_BLOCK_BYTES` in ordinal order, as its sequence files store them, and checks
/// them: a reader error is NCBI's error of the preliminary search (`search_read_error`: NCBI
/// reads every sequence there); a sequence `program` cannot search as NCBI does is LOSAT's
/// refusal (`check_protein_sequence`). The number of residues of each sequence is handed to
/// `each` (ordinal id, residues) and its bytes are dropped with its block. BLASTX's passes
/// (`blastx/x_db_stream.rs`) check the blocks they read the same way and take this order's first
/// error where one fails (`BlastxDatabase::x_check_stores`). NCBI searches the stored bytes in place
/// (blast_engine.c:1404, `eBlastEncodingProtein`) and its traceback and composition redo read
/// them again; the passes read the bytes NCBI does.
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

    // Line DBID pfbatch: each batch's prefetch counts its own Seq-aligns against a positive
    // `PRE_FETCH_SEQS_LIMIT` (`s_PreFetchSeqs`, quoted at `prefetch_may_be_off`).
    #[test]
    fn each_batch_prefetches_by_its_own_seq_aligns() {
        // Seq-aligns of the queries 0 to 4 of a call.
        let counts = [2, 1, 0, 3, 1];
        let of = |limit, first, batches: &[std::ops::Range<usize>]| {
            let prefetch = BatchPrefetch::at(limit, counts.len(), first, batches, |span| {
                counts[span].iter().sum()
            });
            // A query of no batch (`None`, taken both ways) shows as `true`.
            (0..counts.len() + 1)
                .map(|query| prefetch.batch_of(query).unwrap_or(true))
                .collect::<Vec<_>>()
        };
        let every = vec![true; counts.len() + 1];
        // No limit and `INT_MAX`: every batch prefetches; 0 or below: none does.
        for limit in [None, Some(i32::MAX)] {
            assert_eq!(of(limit, 0, &[0..2, 2..5]), every, "{limit:?}");
        }
        for limit in [Some(0), Some(-1)] {
            assert_eq!(of(limit, 0, &[0..2, 2..5]), [false; 6], "{limit:?}");
        }
        // Batches of 3 and 4 Seq-aligns: limit 3 turns the second off; a query beyond the call
        // is of no batch.
        assert_eq!(
            of(Some(3), 0, &[0..2, 2..5]),
            [true, true, false, false, false, true]
        );
        assert_eq!(of(Some(4), 0, &[0..2, 2..5]), every);
        // No batches: the call is one batch of 7 Seq-aligns.
        assert_eq!(
            of(Some(6), 0, &[]),
            [false, false, false, false, false, true]
        );
        assert_eq!(of(Some(7), 0, &[]), every);
        // Ranges numbered from the call's first query 10; a batch before it is not this call's.
        assert_eq!(
            of(Some(1), 10, &[3..10, 10..11, 11..13, 13..15]),
            [false, true, true, false, false, true]
        );
        // Line DBID weakbatch: what the batches read of each subject they report, by ordinal id:
        // 7 in both batches, 8 in the first (3 Seq-aligns, on under limit 3), 9 in the second
        // (4, off), whose reports read it one by one.
        let sum = |span: std::ops::Range<usize>| counts[span].iter().sum();
        let reads = BatchPrefetch::at(Some(3), 5, 0, &[0..2, 2..5], sum).reads([
            (0, 7),
            (4, 7),
            (1, 8),
            (3, 9),
        ]);
        assert!(reads.prefetched(7) && reads.one_by_one(7));
        assert!(reads.prefetched(8) && !reads.one_by_one(8));
        assert!(!reads.prefetched(9) && reads.one_by_one(9));
        let every = BatchPrefetch::at(None, 5, 0, &[0..2, 2..5], sum).reads([(3, 9)]);
        assert!(every.prefetched(9) && !every.one_by_one(9));
        let none = BatchPrefetch::at(Some(0), 5, 0, &[0..2, 2..5], sum).reads([(0, 7)]);
        assert!(!none.prefetched(7) && none.one_by_one(7));
        // Review 1 M1: a query of no batch (5, beyond the call's 5) is taken both ways.
        let beyond = BatchPrefetch::at(Some(3), 5, 0, &[0..2, 2..5], sum).reads([(5, 6)]);
        assert!(beyond.prefetched(6) && beyond.one_by_one(6));
        let shown = ParsedShown {
            prefetched: true,
            ..ParsedShown::default()
        };
        let of_9 = reads.shown(shown, 9);
        assert!(!of_9.prefetched && of_9.one_by_one);
        assert!(prefetch_off_at(Some(1), 2) && !prefetch_off_at(Some(2), 2));
        assert!(prefetch_off_at(Some(0), 0) && !prefetch_off_at(Some(i32::MAX), usize::MAX));
        assert!(!prefetch_off_at(None, usize::MAX));
    }

    /// The subjects of the database `name` of `tests/db_search/ids` by ordinal id.
    fn parsed_subjects(name: &str) -> Vec<crate::report::db_subject::ParsedSubject> {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids")
            .join(name);
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
        (0..db.num_seqs())
            .map(|ordinal| {
                crate::report::db_subject::ParsedSubject::new(db.defline_set(ordinal).unwrap())
                    .unwrap()
            })
            .collect()
    }

    // Line DBID N3 review 1 M1: the joined copies' letters are the same because `copy_of` reads
    // them: the definition lines of r1 and r3 with other letters (pjoin_x: four substitutions, ten
    // more letters) are no copies, those with other taxids (pjoin_abt) are.
    #[test]
    fn a_copy_has_the_same_letters_and_definition_lines_but_for_taxids() {
        let open = |name: &str| {
            let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/db_search/ids")
                .join(name);
            SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap()
        };
        let copy = |db: &SeqDb, a: u32, b: u32| {
            let subject = match db.read_header(a).unwrap() {
                crate::seqdb::SequenceHeader::Parsed(set) => {
                    crate::report::db_subject::ParsedSubject::new(set).unwrap()
                }
                _ => panic!("parsed Seq-ids"),
            };
            HeaderReader::new(db).copy_of(&subject, a, b)
        };
        // pjoin_ax: a (r1 0, r2 1, r3 2), x (r1 3, r3 4).
        let ax = open("pjoin_ax");
        assert_eq!(copy(&ax, 0, 3), RecordCopy::Differs);
        assert_eq!(copy(&ax, 2, 4), RecordCopy::Differs);
        // pjoin_abc / pjoin_abt: a (r1 0, r2 1, r3 2), b (r1 3, r3 4, r4 5, r6 6), c (r1 7).
        let abc = open("pjoin_abc");
        assert_eq!(copy(&abc, 0, 3), RecordCopy::Same);
        assert_eq!(copy(&abc, 2, 4), RecordCopy::Same);
        assert_eq!(copy(&abc, 0, 6), RecordCopy::Differs);
        let abt = open("pjoin_abt");
        assert_eq!(copy(&abt, 0, 3), RecordCopy::OtherTaxids);
        assert_eq!(copy(&abt, 0, 7), RecordCopy::Same);
    }

    // Line DBID N3: in NCBI's order of pjoin_abc's hits (r1 of c, r6, r1 of b and a, r2, r4,
    // r3 of b and a; `make_join_fasta.py`) the cut to N subjects keeps the copies that follow a
    // kept subject (psiblast, `PruneSeqalign`), and the hits show them in one hit; pfound_v5's
    // pair 4 (gb|DDD00004| then gb|DDD00004.1|) is refused where the count reaches it.
    #[test]
    fn copies_count_as_one_subject_in_the_cut_and_other_matching_ids_are_refused() {
        let s = parsed_subjects("pjoin_abc");
        let order = [7u32, 6, 3, 0, 1, 5, 4, 2];
        let hits = || {
            order
                .iter()
                .map(|&ordinal| Some((ordinal, &s[ordinal as usize])))
        };
        for (max, kept) in [(1, 1), (2, 2), (3, 4), (4, 5), (5, 6), (6, 8), (7, 8)] {
            assert_eq!(
                joined_cut(0, max, hits(), "PSIBLAST").unwrap(),
                kept,
                "{max}"
            );
        }
        assert_eq!(joined_cut(0, 0, hits(), "PSIBLAST").unwrap(), 0);
        refuse_joined_hits(0, hits(), "BLASTP").unwrap();
        // A subject without parsed Seq-ids between copies parts them.
        let parted = [Some((3, &s[3])), None, Some((0, &s[0]))];
        assert_eq!(joined_cut(0, 1, parted, "PSIBLAST").unwrap(), 1);
        let f = parsed_subjects("pfound_v5");
        let pair = f
            .iter()
            .enumerate()
            .filter(|(_, subject)| subject.alignment_fasta().starts_with(b"gb|DDD00004"))
            .map(|(ordinal, subject)| Some((ordinal as u32, subject)))
            .collect::<Vec<_>>();
        assert_eq!(pair.len(), 2);
        let error = refuse_joined_hits(1, pair.clone(), "BLASTP")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("but are not copies of one record"),
            "{error}"
        );
        // The cut to one subject reaches the second of the pair: NCBI would keep both.
        let error = joined_cut(1, 1, pair.clone(), "PSIBLAST")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("the cut to -max_target_seqs 1 counts"),
            "{error}"
        );
        // A cut that ends before the pair does not reach it.
        let before = [Some((7, &s[7]))].into_iter().chain(pair);
        assert_eq!(joined_cut(1, 1, before, "PSIBLAST").unwrap(), 1);
    }

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

    /// The check as an `all` over the residues.
    fn accepted_by_all(residues: &[u8]) -> bool {
        residues
            .iter()
            .all(|&residue| usize::from(residue) < BLASTAA_SIZE)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593 (quoted at
    // `check_protein_sequence`): the bytes past the matrix are what the check refuses; the gap
    // byte is a residue of the matrix. Not a port of an NCBI unit case.
    #[test]
    fn the_check_accepts_the_letters_and_names_the_first_bad_byte() {
        let first = NCBISTDAA_GAP + 1;
        let last = (BLASTAA_SIZE - 1) as u8;
        assert_eq!(check(&[]), Ok(()));
        assert_eq!(check(&[first]), Ok(()));
        assert_eq!(check(&[last]), Ok(()));
        assert_eq!(check(&[NCBISTDAA_GAP]), Ok(()));
        let letters: Vec<u8> = (0..1_000)
            .map(|i| [first, last, NCBISTDAA_GAP][i % 3])
            .collect();
        assert_eq!(check(&letters), Ok(()));
        for bad in [BLASTAA_SIZE as u8, 255] {
            for length in [1, 2, 31, 32, 33, 100, 1_000] {
                for at in [0, length / 2, length - 1] {
                    let mut residues = letters[..length].to_vec();
                    residues[at] = bad;
                    let error = check(&residues).unwrap_err();
                    assert!(
                        error.contains("no NCBIstdaa residue")
                            && error.contains(&format!("position {},", at + 1)),
                        "byte {bad} at {at} of {length}: {error}"
                    );
                }
            }
        }
        // The sentinels.
        assert!(check_protein_sequence(&[1, first, 0], 7, "blastx").is_err());
        assert!(check_protein_sequence(&[0, first, 1], 7, "blastx").is_err());
    }

    // Random residues around the letters' bounds: the single pass accepts exactly what an `all`
    // accepts.
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
                accepted_by_all(&residues),
                "{residues:?}"
            );
        }
    }
}

/// The taxonomy database of a run (`CTaxDBFileInfo`, a static of NCBI's `GetTaxNames`) and the
/// lookups the volumes have made (`CSeqDBVol::x_GetTaxonomy` caches them per volume), in the
/// order in which LOSAT records them ([`TaxdbRun::load`]).
///
/// NCBI's cache of a volume has 256 slots, one per `taxid & 255` (`CSeqDBVol::m_TaxCache`,
/// below): a lookup asks the taxonomy database where its slot holds nothing or another taxid,
/// and the taxid then takes the slot. Of the texts a lookup can write, those of opening come once per
/// process and only the offset error of the last entry (`TaxDb::offset_error_taxid`, `L`) comes
/// again at each lookup that misses, so the record keeps what decides the bytes
/// ([`VolumeLookups`]): each (volume, taxid) looked up (the taxonomy database is asked at the
/// first), and per volume what `L`'s slot holds, which only the lookups of `L` and of another
/// positive taxid of its slot change. A run whose volumes hold no such other taxid looks `L` up
/// once per volume, in any order.
///
/// Where the order of the fetches is not known, a run with a volume that holds `L` and another
/// taxid of its slot is refused where NCBI first fetches a subject: a prefetch
/// ([`prefetch_messages`]) or, in the XML2 and JSON formats, a report with a subject that has a
/// positive taxid ([`read_reported_taxonomy`]) ([`TaxdbRun::cache_collision`]). psiblast's PSSM
/// reads its subjects one by one in a known order after each query's report
/// ([`pssm_load_texts`]), so its loads follow the slot; there only a read of a sequence that the
/// run loaded before, which NCBI's scope may keep or load again, is refused where the two give
/// other bytes ([`VolumeLookups::offset_slot`]).
///
/// It lives on the database of the search, so a `run_local` call that writes outfmt 15 and 16
/// of one `-db` search shares it: the open and the volume-lookup messages are written for the
/// first of the two only (the CLI writes one format). NCBI, which runs each format alone, writes
/// them for both; the ABI callers that ask for both do not see them twice.
#[derive(Default)]
pub(crate) struct TaxdbRun {
    lazy: std::sync::Mutex<Option<std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>>>,
    volume_lookups: std::sync::Mutex<VolumeLookups>,
    found: std::sync::OnceLock<bool>,
    /// What the lookups of the run can write ([`TaxdbRun::prefetch_texts`]).
    texts: std::sync::OnceLock<PrefetchTexts>,
    /// The refusal of a run whose fetches NCBI's cache makes write the offset error more than
    /// once, as a reader's refusal ([`TaxdbRun::cache_collision`]).
    collision: std::sync::OnceLock<Option<(std::path::PathBuf, String)>>,
    /// The texts of opening the database when the rows of a tabular report opened it, kept for
    /// those rows ([`TaxdbRun::open_for_rows`], [`tabular_names`]).
    rows_open: std::sync::Arc<std::sync::Mutex<Vec<String>>>,
}

/// What the lookups of a run's taxonomy database can write ([`TaxdbRun::prefetch_texts`]).
#[derive(Clone, Copy)]
struct PrefetchTexts {
    /// The positive taxid whose lookup writes the offset error (`TaxDb::offset_error_taxid`).
    offset_error_taxid: Option<i32>,
}

/// `CSeqDBVol::m_TaxCache(256)`: the slots of a volume's taxonomy cache.
const TAX_CACHE_SLOTS: i32 = 256;

/// The slot of NCBI's taxonomy cache of a volume that `taxid` takes (`CSeqDBIntCache::Lookup`,
/// `key & (m_Slots.size()-1)`).
fn tax_cache_slot(taxid: i32) -> i32 {
    taxid & (TAX_CACHE_SLOTS - 1)
}

/// The volume lookups of a run ([`TaxdbRun`]), in the order in which LOSAT records them
/// ([`TaxdbRun::load`]).
#[derive(Default)]
struct VolumeLookups {
    /// Each (volume, taxid) looked up: NCBI asks the taxonomy database at the first.
    asked: std::collections::HashSet<(usize, i32)>,
    /// By volume, what the slot of the offset error's taxid holds ([`OffsetSlot`]); only where
    /// the taxonomy database has that error.
    offset_slots: std::collections::HashMap<usize, OffsetSlot>,
    /// The sequences, by ordinal id, whose loads looked up a taxid of that slot: a later read of
    /// one may find it in NCBI's scope and load nothing, or load it again.
    slot_loads: std::collections::HashSet<u32>,
}

/// What the slot of the offset error's taxid `L` holds in a volume's taxonomy cache.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum SlotHolds {
    /// Nothing or another taxid of the slot: a lookup of `L` asks the taxonomy database.
    #[default]
    Other,
    /// `L`: a lookup of `L` finds it.
    Last,
    /// `L` or another taxid, as NCBI's scope kept a sequence read before or loaded it again.
    Either,
}

/// A load that put a taxid of the slot of the offset error's taxid there, other than that taxid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SlotFetch {
    /// The sequence loaded, by ordinal id.
    ordinal: u32,
    /// The other taxid.
    taxid: i32,
    /// Whether it was a read again of a sequence loaded before (NCBI may have loaded nothing).
    again: bool,
}

/// The slot of the offset error's taxid in a volume's taxonomy cache.
#[derive(Clone, Copy, Debug, Default)]
struct OffsetSlot {
    holds: SlotHolds,
    /// The latest load that put another taxid of the slot there, for LOSAT's refusal.
    other: Option<SlotFetch>,
}

/// A load of the sequence `ordinal` whose texts of the offset error NCBI's scope decides
/// ([`VolumeLookups::offset_slot`]): `again`, a read of a sequence loaded before, where keeping
/// it and loading it again give other texts; else a load after such a read. `after`: the load of
/// the other taxid of the slot before it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SlotUncertain {
    ordinal: u32,
    again: bool,
    after: Option<SlotFetch>,
}

/// One load's lookups of the slot of the offset error's taxid ([`SlotHolds::walk`]).
struct SlotWalk {
    /// What the slot holds after them.
    holds: SlotHolds,
    /// Which lookups of the taxid of the offset error miss the slot (NCBI asks the taxonomy
    /// database again).
    misses: Vec<bool>,
    /// The last other taxid of the slot that they look up.
    other: Option<i32>,
}

impl SlotHolds {
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
    // A volume's cache (`CSeqDBIntCache`) keeps one taxid per slot `taxid & 255`, which a lookup
    // of another taxid of the slot clears; `x_GetTaxonomy` asks `GetTaxNames` where it misses.
    // `use_taxinfo_cache` is true on the main thread, where psiblast makes its PSSMs.
    /// The lookups `lookups` of one load (in `x_GetTaxonomy`'s order), from this state of the slot
    /// of the offset error's taxid `last`: a lookup of `last` misses where the slot holds another
    /// taxid, and takes the slot; a lookup of another taxid of the slot takes it. `Err` where a
    /// lookup of `last` meets `Either`, with the other taxid of the slot that this load looked up
    /// before it.
    fn walk(self, lookups: &[i32], last: i32) -> std::result::Result<SlotWalk, Option<i32>> {
        let mut holds = self;
        let mut misses = vec![false; lookups.len()];
        let mut other = None;
        for (index, &taxid) in lookups.iter().enumerate() {
            if taxid == last {
                match holds {
                    SlotHolds::Last => {}
                    SlotHolds::Other => {
                        misses[index] = true;
                        holds = SlotHolds::Last;
                    }
                    SlotHolds::Either => return Err(other),
                }
            } else if tax_cache_slot(taxid) == tax_cache_slot(last) {
                holds = SlotHolds::Other;
                other = Some(taxid);
            }
        }
        Ok(SlotWalk {
            holds,
            misses,
            other,
        })
    }
}

impl VolumeLookups {
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:453-462
    // ```c++
    // CBlastDbDataLoader::GetBlobById(const TBlobId& blob_id)
    // {
    //     CTSE_LoadLock lock = GetDataSource()->GetTSE_LoadLock(blob_id);
    //     if ( !lock.IsLoaded() ) {
    //         const TBlastDbId& id =
    //             dynamic_cast<const CBlobIdBlastDb&>(*blob_id).GetValue();
    //         x_LoadData(id.second, id.first, lock, kSequenceSliceSize);
    //     }
    //     return lock;
    // }
    // ```
    // A sequence that no fetch of the run has loaded is loaded (`x_LoadData`, whose
    // `CCachedSequence` asks `GetBioseqNoData` and so `x_GetTaxonomy`). One loaded before may
    // still be in NCBI's scope or data source, which then gives it without a load, or have been
    // released: review 4 of line DBID pssmread found a second query's PSSM loading a sequence again
    // after 40 other included subjects and not after none (`cc_big40`, `cc_ac`), which LOSAT does
    // not model.
    /// Which of the lookups `lookups` of a load of the sequence `ordinal` in the volume `volume`
    /// ask the taxonomy database again for the taxid of the offset error `last` (a lookup that
    /// misses the slot after an earlier lookup of `last` in the volume), as NCBI's cache of the
    /// volume holds the slot of `last` after the lookups recorded before. The first load of a
    /// sequence in the run follows the slot. A read again of a sequence whose load looks up a
    /// taxid of the slot is a load or none, as NCBI's scope decides: where loading it again gives
    /// a text, LOSAT's refusal; else no text, and the slot holds what either leaves (`Either` where
    /// they differ). `Err` names the load that is not certain ([`SlotUncertain`]).
    fn offset_slot(
        &mut self,
        volume: usize,
        ordinal: u32,
        lookups: &[i32],
        last: i32,
    ) -> std::result::Result<Vec<bool>, SlotUncertain> {
        let in_slot = |taxid: i32| tax_cache_slot(taxid) == tax_cache_slot(last);
        if !lookups.iter().any(|&taxid| in_slot(taxid)) {
            return Ok(vec![false; lookups.len()]);
        }
        let again = !self.slot_loads.insert(ordinal);
        let slot = self.offset_slots.entry(volume).or_default();
        let before = *slot;
        let uncertain = |own: Option<i32>| SlotUncertain {
            ordinal,
            again,
            after: own
                .map(|taxid| SlotFetch {
                    ordinal,
                    taxid,
                    again,
                })
                .or(before.other),
        };
        let walk = before.holds.walk(lookups, last).map_err(uncertain)?;
        if !again {
            slot.holds = walk.holds;
            if let Some(taxid) = walk.other {
                slot.other = Some(SlotFetch {
                    ordinal,
                    taxid,
                    again: false,
                });
            }
            return Ok(walk.misses);
        }
        if let Some(index) = walk.misses.iter().position(|&miss| miss) {
            let own = lookups[..index]
                .iter()
                .rev()
                .find(|&&taxid| taxid != last && in_slot(taxid))
                .copied();
            return Err(uncertain(own));
        }
        if walk.holds != before.holds {
            slot.holds = SlotHolds::Either;
            if let Some(taxid) = walk.other {
                slot.other = Some(SlotFetch {
                    ordinal,
                    taxid,
                    again: true,
                });
            }
        }
        Ok(vec![false; lookups.len()])
    }
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

    // NCBI's 256-slot cache of a volume (quoted at `SlotHolds::walk`): a fetch of a sequence with
    // another taxid of the slot of `L` clears it, and the next fetch of `L` asks `GetTaxNames`
    // again, which writes the offset error again. How often it repeats in a prefetch follows the
    // order in which the object manager loads the batch's subjects, which is neither the order of
    // the reports nor that of the ordinal ids (review 1 of item T: one query on a volume of the
    // taxids 99990, 150, 99990 writes the message twice, two queries once; fix 1: the same query
    // once with a version 4 database or an absolute `BLASTDB`), so LOSAT refuses such a run there.
    /// LOSAT's refusal of a run whose fetches (a prefetch, or the reports of outfmt 13 to 16)
    /// NCBI would make write the offset error more than once in a volume: the taxonomy database
    /// has the offset error of a positive taxid `L` ([`TaxdbRun::prefetch_texts`]), and a volume
    /// of the database holds `L` and another
    /// positive taxid of the same cache slot (`T & 255 == L & 255`), in any definition line of any
    /// sequence, whether or not a search reports them. One pass over the headers
    /// ([`SeqDb::visit_defline_set_taxids`] in ranges of [`crate::seqdb::HEADER_RANGE_SEQS`]), made
    /// only for such a damaged taxonomy database and once for the run; a header whose taxids do
    /// not decode is left to the lookups of the batches, which fail on it where it is reported.
    /// `(path of taxdb.bti, reason)` of the refusal, which [`prefetch_messages`] and
    /// [`read_reported_taxonomy`] make as a reader's refusal ([`TaxonomyRefusal`]) where NCBI
    /// fetches ([`TaxdbRun::refuse_cache_collision`]).
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

    /// Whether a fetch of the run is refused because of NCBI's taxonomy cache
    /// ([`TaxdbRun::cache_collision`]): BLASTX holds its XML document (outfmt 5, and 13 to 16)
    /// then, which the refusal drops.
    pub(crate) fn prefetch_refused(&self, db: &SeqDb) -> bool {
        self.cache_collision(db).is_some()
    }

    /// LOSAT's refusal of a fetch of the run where NCBI's taxonomy cache would repeat the offset
    /// error ([`TaxdbRun::cache_collision`]), as a reader's refusal ([`TaxonomyRefusal`]) that names
    /// `program`; `Ok` for any other run. The callers ask it only where NCBI fetches a subject.
    fn refuse_cache_collision(&self, db: &SeqDb, program: &str) -> Result<()> {
        match self.cache_collision(db) {
            Some((path, message)) => Err(taxdb_error(
                SeqDbError::Unsupported {
                    path: path.clone(),
                    message: message.clone(),
                },
                program,
            )),
            None => Ok(()),
        }
    }

    /// The open of the database by the rows of a query batch whose report looks names up, made
    /// after the batch's prefetch ([`prefetch_messages`]) and before the next batch's: NCBI writes
    /// a batch's rows before it fetches the next batch, and a plural name looks taxid 0 up, which
    /// the prefetch does not, so the run's first lookup can be such a row's. The texts of opening,
    /// when this call opens the database, are kept for the rows ([`tabular_names`]). Nothing
    /// without `taxdb.bti` on the search path.
    pub(crate) fn open_for_rows(&self, db: &SeqDb, program: &str) -> Result<()> {
        if !self.found(db) {
            return Ok(());
        }
        let Some(search_path) = db.search_path() else {
            return Ok(());
        };
        let texts = self
            .taxdb(search_path)
            .open()
            .map_err(|error| taxdb_error(error, program))?;
        self.rows_open
            .lock()
            .expect("the texts of opening")
            .extend(texts);
        Ok(())
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

    /// NCBI's lookups of the volume's taxonomy cache where it loads the sequence `ordinal` of `db`:
    /// the positive taxids `lookups` in the order of `CSeqDBVol::x_GetTaxonomy`, a taxid of several
    /// definition lines once for each ([`blob_taxids`]). The taxonomy database is asked at the
    /// volume's first lookup of a taxid in the run (the texts of opening come with the run's first,
    /// `LazyTaxDb::names`) and again where a lookup of the offset error's taxid misses its slot of
    /// the cache ([`VolumeLookups::offset_slot`]). Returns what each of those asks wrote, by taxid,
    /// in the order of the lookups. A load whose texts NCBI's scope decides (a read again of a
    /// sequence loaded before, where its scope keeping it and loading it again write other texts,
    /// or a load after such a read) is LOSAT's refusal ([`TaxonomyRefusal`]), naming `program`.
    fn load(
        &self,
        db: &SeqDb,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        ordinal: u32,
        lookups: &[i32],
        program: &str,
    ) -> Result<Vec<(i32, Vec<String>)>> {
        if lookups.is_empty() {
            return Ok(Vec::new());
        }
        let volume = db
            .volume_index(ordinal)
            .map_err(|error| database_error(error, program))?;
        let last = self.prefetch_texts(db).offset_error_taxid;
        let mut record = self.volume_lookups.lock().expect("the volume lookups");
        let misses = match last {
            Some(last) => record
                .offset_slot(volume, ordinal, lookups, last)
                .map_err(|cause| slot_refusal(db, volume, last, cause, program))?,
            None => vec![false; lookups.len()],
        };
        let mut answers = Vec::new();
        for (&taxid, miss) in lookups.iter().zip(misses) {
            if record.asked.insert((volume, taxid)) || miss {
                let mut messages = Vec::new();
                taxdb
                    .names(taxid, &mut messages)
                    .map_err(|error| taxdb_error(error, program))?;
                answers.push((taxid, messages));
            }
        }
        Ok(answers)
    }

    /// [`TaxdbRun::load`] of the sequence `ordinal`, its texts appended to `text` as `ERR_POST`
    /// lines of the program `lower`.
    #[allow(clippy::too_many_arguments)]
    fn volume_lookups_of(
        &self,
        db: &SeqDb,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        ordinal: u32,
        lookups: &[i32],
        lower: &str,
        text: &mut String,
        program: &str,
    ) -> Result<()> {
        for (_, messages) in self.load(db, taxdb, ordinal, lookups, program)? {
            for line in messages {
                text.push_str(&crate::seqdb::taxdb::err_post_line(lower, &line));
                text.push('\n');
            }
        }
        Ok(())
    }
}

/// LOSAT's refusal of a load in the volume `volume` of `db` whose texts of the offset error of
/// taxid `last` NCBI's scope decides ([`SlotUncertain`], [`VolumeLookups::offset_slot`]), as a
/// reader's refusal ([`TaxonomyRefusal`]) of `taxdb.bti` that names the read and the load of the
/// other taxid of the slot before it, and `program`.
fn slot_refusal(
    db: &SeqDb,
    volume: usize,
    last: i32,
    cause: SlotUncertain,
    program: &str,
) -> anyhow::Error {
    let volume_name = db
        .volumes()
        .get(volume)
        .map_or_else(String::new, |vol| vol.base().display().to_string());
    let path = db
        .search_path()
        .and_then(|search_path| search_path.resolve_file("taxdb.bti"))
        .map_or_else(
            || std::path::PathBuf::from("taxdb.bti"),
            std::path::PathBuf::from,
        );
    let read = if cause.again {
        format!(
            "reads database sequence {} of the volume {volume_name}, which holds it, again",
            cause.ordinal
        )
    } else {
        format!(
            "fetches database sequence {} of the volume {volume_name}, which holds it,",
            cause.ordinal
        )
    };
    let after = match cause.after {
        Some(fetch) if fetch.ordinal == cause.ordinal => {
            format!("after its own taxid {}", fetch.taxid)
        }
        Some(SlotFetch {
            ordinal,
            taxid,
            again: false,
        }) => format!("after the run's fetch of database sequence {ordinal} with taxid {taxid}"),
        Some(SlotFetch {
            ordinal,
            taxid,
            again: true,
        }) => format!(
            "after the run read database sequence {ordinal} again, whose load looks up taxid {taxid}"
        ),
        None => "after another taxid".to_string(),
    };
    let scope = match cause.after {
        _ if cause.again => {
            "there again where its scope loads that sequence again and not where it keeps it"
                .to_string()
        }
        Some(fetch) => format!(
            "there or not as its scope loaded database sequence {} again or kept it",
            fetch.ordinal
        ),
        None => "there or not as its scope kept the sequences read before".to_string(),
    };
    taxdb_error(
        SeqDbError::Unsupported {
            path,
            message: format!(
                "taxid {last}, the last entry of the taxonomy database, begins past the end of \
                 taxdb.btd, and the run {read} {after}, which takes the same slot of NCBI \
                 BLAST+'s 256-slot taxonomy cache of a volume (CSeqDBVol::m_TaxCache, \
                 CSeqDBIntCache): NCBI writes \"{}\" {scope}, which LOSAT does not reproduce",
                crate::seqdb::taxdb::MSG_OFFSET_ERROR
            ),
        },
        program,
    )
}

/// LOSAT's refusal of a taxonomy database name or record (`read_reported_taxonomy`): a scientific
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
/// byte is written; so is a run whose fetches NCBI's taxonomy cache would make repeat the
/// offset error ([`TaxdbRun::cache_collision`]), when a reported subject has a positive taxid.
///
/// NCBI's messages about the taxonomy database go to standard error, and not here: those of
/// opening and of the lookup that the volume makes for each taxid of a fetched subject
/// (`CSeqDBVol::x_GetTaxonomy`, once per volume and taxid) are in the returned [`TaxdbNotes`],
/// which writes them at the first `Hit` that needs them, and those of the lookup of each `Hit`
/// are in `taxonomy` and written at the hit ([`TaxdbPrinter::at_hit`]). `program` is the
/// program that refuses, as LOSAT's refusals name it (`BLASTP`); its messages name it in lower
/// case.
///
/// The database's reported sequences can have parsed Seq-ids (`parsed[rank]`, empty where none
/// has): the third vector holds, by rank, what the taxonomy database said about the taxid of each
/// definition line of such a sequence (empty where nothing was asked), as NCBI's `HitDescr` loop
/// asks it (`SubjectNames::of_parsed`): `Unnamed` for a definition line without a taxid, and no
/// answer after the first `NotFound` (NCBI's loop ends there). The volume's lookups are those of
/// `ParsedSubject::volume_taxids`.
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
    // Whether a volume looks a taxid up when NCBI fetches the reported subjects (a positive one,
    // `x_GetTaxonomy`; `ParsedSubject::volume_taxids`).
    let mut volume_lookup = false;
    for &rank in reported {
        if let Some(subject) = parsed_subject(rank) {
            any_taxid |= defline_taxids(subject).iter().any(|&taxid| taxid != 0);
            volume_lookup |= !subject.volume_taxids().is_empty();
            continue;
        }
        let defline_taxids = db
            .defline_taxids(ordinals[rank])
            .map_err(|error| database_error(error, program))?;
        // One definition line per reported sequence (`stored_title`).
        taxids[rank] = defline_taxids.first().copied().unwrap_or(0);
        any_taxid |= taxids[rank] != 0;
        volume_lookup |= taxids[rank] > 0;
    }
    if !any_taxid || !db.taxdb_found() {
        return Ok((taxids, Vec::new(), Vec::new(), TaxdbNotes::default()));
    }
    let Some(search_path) = db.search_path() else {
        return Ok((taxids, Vec::new(), Vec::new(), TaxdbNotes::default()));
    };
    // NCBI fetches every reported subject of these formats (`s_IsPrefetchFormat` is true for
    // them; with the prefetch off the report fetches each one), so a volume's lookups can repeat
    // the offset error as in the prefetch of the other formats: the same refusal, asked where a
    // volume looks a taxid up.
    if volume_lookup {
        run.refuse_cache_collision(db, program)?;
    }
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
        if let std::collections::btree_map::Entry::Vacant(slot) = asked.entry(taxid) {
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
            slot.insert(answer);
        }
        Ok(asked[&taxid].clone())
    };
    let mut taxonomy = vec![HitTaxonomy::Unnamed; taxids.len()];
    let mut parsed_taxonomy = vec![Vec::new(); if parsed.is_empty() { 0 } else { taxids.len() }];
    notes.volumes = vec![0; taxids.len()];
    for &rank in reported {
        let subject = parsed_subject(rank);
        // The volume's lookups when it fetches the sequence (`TaxdbRun::load`; `x_GetTaxonomy`
        // asks for a positive taxid only, seqdbvol.cpp:1136).
        let volume_lookups = match subject {
            Some(subject) if defline_taxids(subject).iter().all(|&taxid| taxid == 0) => continue,
            Some(subject) => subject.volume_lookups(),
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
        for (taxid, messages) in run.load(db, &taxdb, ordinals[rank], &volume_lookups, program)? {
            notes
                .lookups
                .entry((volume, taxid))
                .or_default()
                .extend(messages);
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

/// Whether `taxdb.bti` is on NCBI's search path for a `-subject` search (`x_CheckTaxDB`'s
/// `SeqDB_ResolveDbPath("taxdb.bti")`: the current directory, `BLASTDB` and the registry's
/// `[BLAST] BLASTDB`, as a database is looked for, `seqdb::SearchPath`). A registry value that
/// LOSAT does not read as NCBI is its refusal, as for `-db`; `program` names it.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
/// ```c++
/// void CBlastTabularInfo::x_CheckTaxDB()
/// {
///    	if( x_IsFieldRequested(eSubjectSciNames) ||
///     ...
///    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
///    		if(resolved.empty())
/// ```
pub(crate) fn subject_taxdb_found(settings: &ApplicationSettings, program: &str) -> Result<bool> {
    Ok(subject_search_path(settings, program)?
        .resolve_file("taxdb.bti")
        .is_some())
}

/// The search path of NCBI's databases of a `-subject` search (`subject_taxdb_found`).
fn subject_search_path(settings: &ApplicationSettings, program: &str) -> Result<SearchPath> {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => bail!("{reason}; this is not supported by LOSAT's {program}"),
    };
    SearchPath::from_environment(registry).map_err(|error| database_error(error, program))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-220
// ```c++
// bool CSeqDBTaxInfo::GetTaxNames(TTaxId           tax_id,
//                                 SSeqDBTaxInfo  & info  )
// {
// 	static CTaxDBFileInfo t;
// ```
/// The taxonomy database of a `-subject` search's organism report (outfmt 18): NCBI's
/// `GetTaxNames` opens its static `CTaxDBFileInfo` through the search path of NCBI's databases at
/// the report's first lookup. One per run, as a `-db` search's (`TaxdbRun`), made by the command
/// line with `subject_taxdb_found` (`run_local`'s callers give no search path: none, DECISIONS 7
/// of line CLI C2, every lookup fails).
#[derive(Clone, Debug)]
pub struct SubjectTaxdb(std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>);

impl SubjectTaxdb {
    /// The run's database where `taxdb.bti` is on the search path of `settings` (else `None`),
    /// not opened before the report's first lookup. A registry value that LOSAT does not read as
    /// NCBI is its refusal (`subject_taxdb_found`); `program` names it.
    pub(crate) fn of(settings: &ApplicationSettings, program: &str) -> Result<Option<Self>> {
        let path = subject_search_path(settings, program)?;
        Ok(path.resolve_file("taxdb.bti").is_some().then(|| {
            SubjectTaxdb(std::sync::Arc::new(crate::seqdb::taxdb::LazyTaxDb::new(
                path,
            )))
        }))
    }

    /// The names that the organism report of a query batch asks for in it
    /// (`TabularNames::for_report`); a database that LOSAT cannot read as NCBI does is refused as
    /// in [`tabular_names`].
    pub(crate) fn report_names(&self, program: &str) -> crate::report::tabular::TabularNames {
        crate::report::tabular::TabularNames::for_report(
            crate::report::tax_report::TaxNamesDb::new(
                std::sync::Arc::clone(&self.0),
                None,
                taxdb_error,
            ),
            program,
        )
    }
}

/// The names of the taxonomy database that the rows of a tabular report with a name field look
/// up (`report::tabular::TabularNames`) for its subjects `taxonomies`, through the run's
/// database; `None` without a name field or without `taxdb.bti` on the search path (every lookup
/// fails, silently). A taxonomy database that LOSAT cannot read as NCBI does is refused as in
/// [`read_reported_taxonomy`], before any byte of the report.
pub(crate) fn tabular_names<'t>(
    db: &SeqDb,
    run: &TaxdbRun,
    spec: &crate::report::tabular::TabularSpec,
    taxonomies: impl IntoIterator<Item = &'t crate::report::tabular::TabularTaxonomy>,
    program: &str,
) -> Result<Option<crate::report::tabular::TabularNames>> {
    if !spec.needs().taxonomy_names || !run.found(db) {
        return Ok(None);
    }
    let Some(search_path) = db.search_path() else {
        return Ok(None);
    };
    // The organism report (outfmt 18) looks each taxid up itself where NCBI's report does
    // (`report::tax_report::TaxLookups`), the texts of an open by earlier rows first.
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        return Ok(Some(crate::report::tabular::TabularNames::for_report(
            crate::report::tax_report::TaxNamesDb::new(
                run.taxdb(search_path),
                Some(std::sync::Arc::clone(&run.rows_open)),
                taxdb_error,
            ),
            program,
        )));
    }
    let names = crate::report::tabular::TabularNames::look_up(
        spec,
        &run.taxdb(search_path),
        taxonomies,
        program,
    )
    .map_err(|error| taxdb_error(error, program))?;
    // No NCBI counterpart: NCBI writes these texts where its first lookup opens the taxonomy
    // database; LOSAT's report calls are per query batch. The texts of opening of an open by the
    // rows of an earlier batch ([`TaxdbRun::open_for_rows`]), which the first row that looks a
    // name up writes: a report without such a row leaves them for the next report of the run.
    if !names.looks_up() {
        return Ok(Some(names));
    }
    let opened = std::mem::take(&mut *run.rows_open.lock().expect("the texts of opening"));
    Ok(Some(names.opened_before(opened)))
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
/// nothing), so the order of the subjects does not change the bytes. The callers decide where
/// they ask: BLASTX, which writes each batch's reports before it searches the next, calls it
/// from its first batch on, unless no batch can fetch a subject (`prefetch_always_off` with
/// rows that read none); BLASTN, TBLASTN and TBLASTX call it only where NCBI fetches a subject,
/// so a run that fetches none is not refused: TBLASTN, and TBLASTX, whose report stream is off
/// for a run with this refusal, make every batch's lookups before any report, and BLASTN's
/// report stream writes a later batch's refusal after the reports of the batches before, as
/// its other refusals of a batch's reported subjects (`x_rpt_stream`). The caller writes the
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
    let reads = subject_read_texts(
        db,
        run,
        subjects
            .into_iter()
            .map(|(ordinal, parsed)| (0, ordinal, parsed)),
        false,
        program,
    )?;
    Ok(reads.into_iter().flat_map(|(_, text)| text).collect())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:819-838
// ```c++
//     bool setSubjectTaxInfoAll = (x_IsFieldRequested(eSubjectTaxIds) ||
//     						     x_IsFieldRequested(eSubjectSciNames) ||
//     ...
//     if(setSubjectIds || setSubjectTaxInfo || setSubjectTaxInfoAll || setSubjectTitle ||
//        x_IsFieldRequested(eSubjectStrand) || setSubjectId)
//     {
//         try {
//        		const CBioseq_Handle& subject_bh =
//                 scope.GetBioseqHandle(align.GetSeq_id(1));
// ```
/// [`prefetch_messages`] of `subjects`, each `(key, ordinal, parsed Seq-ids)`, by subject: the
/// text that the first read of each subject writes, with the subject's `key` (the name a report
/// writer gives it), in the order of `subjects`; a subject whose read writes nothing has no
/// entry. A query's report whose batch's prefetch is off writes each where it reads the subject
/// (`report::query_warnings::SubjectReads`). `rows_open`: the report's rows look the plural
/// names up (taxid 0 too, `TaxdbRun::open_for_rows`), each row after it loads its subject, so
/// the first row opens the taxonomy database after the first subject's read, where that read
/// did not.
pub(crate) fn subject_read_texts<'s, K>(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: impl IntoIterator<Item = (K, u32, Option<&'s crate::report::db_subject::ParsedSubject>)>,
    rows_open: bool,
    program: &str,
) -> Result<Vec<(K, Vec<u8>)>> {
    if !run.found(db) {
        return Ok(Vec::new());
    }
    let Some(search_path) = db.search_path() else {
        return Ok(Vec::new());
    };
    run.refuse_cache_collision(db, program)?;
    let taxdb = run.taxdb(search_path);
    let lower = program.to_ascii_lowercase();
    let mut reads = Vec::new();
    let mut fetched: std::collections::HashSet<u32> = Default::default();
    for (key, ordinal, parsed) in subjects {
        if !fetched.insert(ordinal) {
            continue;
        }
        let first = fetched.len() == 1;
        // A subject with parsed Seq-ids is read by its Seq-align id (`prefetch_reads`: the blob of
        // the sequence the scope gets, whose taxids its volume looks up). In
        // the formats that read the Bioseq, `ParsedSubject::fetched` refuses every subject that
        // is not itself, with the same reason (BLASTN and TBLASTX ask it after the batch loop, so
        // such a subject can be refused here first); psiblast's kept beyond its cut are itself
        // (`prefetched_ranks`). One without them is read by its ordinal id, its own blob.
        let blob = match parsed {
            Some(subject) => match prefetch_reads(db, ordinal, subject, true, program)? {
                Ok(Some(blob)) => Some(blob),
                Ok(None) => None,
                Err(reason) => return Err(parsed_refusal(ordinal, &reason, program)),
            },
            None => Some(ordinal),
        };
        if let Some(blob) = blob {
            let lookups = blob_taxids(db, ordinal, parsed, blob, program)?;
            let mut text = String::new();
            run.volume_lookups_of(db, &taxdb, blob, &lookups, &lower, &mut text, program)?;
            if !text.is_empty() {
                reads.push((key, text.into_bytes()));
            }
        }
        if first && rows_open {
            run.open_for_rows(db, program)?;
        }
    }
    Ok(reads)
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:513
// ```c++
//             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
// ```
/// The texts of the taxonomy database where psiblast's PSSM of a query whose batch reads one by
/// one reads the sequence `blob` by a Seq-align id, the lookups of its load `lookups` (the
/// positive taxids in `x_GetTaxonomy`'s order, [`blob_taxids`]; `ListedTitles::pssm_load_texts`),
/// which NCBI writes there, after the query's report: the run's record of the volume lookups
/// ([`TaxdbRun::load`], as [`prefetch_messages`] for the reports' reads) asks the taxonomy
/// database at the volume's first lookup of a taxid (the texts of opening with the run's first)
/// and where a lookup of the offset error's taxid misses its slot of the volume's 256-slot cache,
/// which keeps the taxid across the run's queries and rounds. The PSSM's reads come in the order of
/// its Seq-aligns, so the first load of a sequence in the run writes NCBI's texts, a taxid of the
/// slot looked up between two lookups of that taxid included. NCBI's scope does not keep what a
/// query loaded for its later rounds or queries (`ResetScopeHistory` after each query; review 4 of
/// line DBID pssmread: a second query or round loads a sequence again after 40 other included
/// subjects and not after none): a sequence loaded before (the query's report, an earlier query or
/// round) is loaded again or not, which LOSAT does not model, so such a read is refused where the
/// two give other texts (the slot then holds another taxid) and else writes nothing. The lookups are the run's from then on. Empty without `taxdb.bti` on the search
/// path.
pub(crate) fn pssm_load_texts(
    db: &SeqDb,
    run: &TaxdbRun,
    blob: u32,
    lookups: &[i32],
    program: &str,
) -> Result<Vec<u8>> {
    if lookups.is_empty() || !run.found(db) {
        return Ok(Vec::new());
    }
    let Some(search_path) = db.search_path() else {
        return Ok(Vec::new());
    };
    let taxdb = run.taxdb(search_path);
    let mut text = String::new();
    run.volume_lookups_of(
        db,
        &taxdb,
        blob,
        lookups,
        &program.to_ascii_lowercase(),
        &mut text,
        program,
    )?;
    Ok(text.into_bytes())
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3304-3320
// ```c++
//     for ( size_t i = from; i < from + count; ++i ) {
//         ret[i] = GetBioseqHandle(ids[i], CScope::eGetBioseq_Resolved);
//         if ( !ret[i] ) {
//             match_map[ids[i]];
//         }
//     }
//     ...
//     for (CPriority_I it(m_setDataSrc); it; ++it) {
//         it->GetBlobs(match_map);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
// ```c++
// void CDataLoader::GetBlobs(TTSE_LockSets& tse_sets)
// {
//     NON_CONST_ITERATE(TTSE_LockSets, tse_set, tse_sets) {
//         tse_set->second = GetRecords(tse_set->first, eBlob);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:201-208
// ```c++
//     bool operator<  (const CSeq_id_Handle& handle) const
//         {
//             // Packed (m_Packed != 0) first:
//             // zeroes are converted to a highest unsigned value by decrement.
//             TUintId p1 = INT_ID_TO(TUintId, m_Packed-1);
//             TUintId p2 = INT_ID_TO(TUintId, handle.m_Packed-1);
//             return p1 < p2 || (p1 == p2 && m_Info < handle.m_Info);
//         }
// ```
// `LoadSequencesToScope` (quoted at `check_prefetched_parsed`) first loads the blob of every
// subject of the batch (`GetBioseqHandles`), in the order of a map of Seq-id handles: the ids
// that NCBI's Seq-id mapper packs (GIs, accessions of a prefix and digits, the text tag's
// `gnl|BL_ORD_ID|0`) by their packed value, then the others by the address of each id's
// `CSeq_id_Info` (review 1 of line DBID pfwarn, M1: `gb|FFF00093.1|` loads before the text tag,
// `lcl|zzz` never), which LOSAT does not model. A lookup that throws on its id there ends the
// loads, so which subjects were loaded before it is not known here; the throw is caught with
// `Error pre-fetching sequence data` (`ERR_POST`, on standard error, after the texts of the
// loads before). A lookup that fails on a file is a `CSeqDBException`, which the prefetch
// rethrows (`FoundInstead::Failed`, refused where it is read). A minus-strand Seq-align of a subject
// without a Bioseq throws after every blob was loaded (`prefetch_gets_no_bioseq`): the texts of
// the whole batch, then the warning. Either way the prefetch reads no sequence data, which no
// text comes with (probe of line DBID pfwarn, `R/pfwarn/w0`).
/// One reported sequence of a query batch whose prefetch is on ([`batch_prefetch`]).
#[derive(Clone, Copy)]
pub(crate) struct PrefetchSubject<'s> {
    /// The ordinal id.
    pub(crate) ordinal: u32,
    /// The definition lines of a sequence with parsed Seq-ids.
    pub(crate) parsed: Option<&'s crate::report::db_subject::ParsedSubject>,
    /// Whether the batch has a Seq-align of it on the query's minus strand and the subject's plus
    /// strand, whose length NCBI's prefetch asks for.
    pub(crate) flipped: bool,
}

impl<'s> PrefetchSubject<'s> {
    /// A sequence without a minus-strand Seq-align (a protein subject or query).
    pub(crate) fn unflipped(
        ordinal: u32,
        parsed: Option<&'s crate::report::db_subject::ParsedSubject>,
    ) -> Self {
        Self {
            ordinal,
            parsed,
            flipped: false,
        }
    }
}

/// NCBI's prefetch of one query batch's reported subjects on standard error, for a first format
/// that prefetches ([`prefetches_subjects`]) and a batch whose prefetch is on: the taxonomy
/// database's texts of the loads ([`prefetch_messages`], empty without `taxdb.bti` on the search
/// path), then `Warning: [<program>] Error pre-fetching sequence data` where the prefetch fails.
/// `subjects` are the batch's reported sequences ([`PrefetchSubject`]), in any order. The
/// prefetch fails at a subject whose lookup throws on its id (`ParsedSubject::prefetch_fails`),
/// which loads nothing, and at a flipped Seq-align of a subject without a Bioseq
/// ([`prefetch_gets_no_bioseq`]). A batch with a failed lookup whose other subjects' loads write
/// a text is refused: LOSAT does not model which of them NCBI loads before the failure. The caller
/// writes the text where the batch's prefetch messages go, after the reader's messages of the
/// batch and before the warnings of its first query.
pub(crate) fn batch_prefetch<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: impl IntoIterator<Item = PrefetchSubject<'s>>,
    program: &str,
) -> Result<Vec<u8>> {
    let mut failed_lookup: Option<(u32, String)> = None;
    let mut fails = false;
    let mut uncertain: Option<(u32, String)> = None;
    let mut asked: std::collections::HashSet<u32> = Default::default();
    let mut loaded = Vec::new();
    for PrefetchSubject {
        ordinal,
        parsed,
        flipped,
    } in subjects
    {
        if let Some(subject) = parsed {
            if subject.prefetch_fails() {
                failed_lookup.get_or_insert_with(|| {
                    (
                        ordinal,
                        String::from_utf8_lossy(&subject.alignment_fasta()).into_owned(),
                    )
                });
                continue;
            }
            if flipped && !fails && asked.insert(ordinal) {
                match prefetch_gets_no_bioseq(db, ordinal, subject, program)? {
                    Ok(no_bioseq) => fails = no_bioseq,
                    Err(reason) => {
                        uncertain.get_or_insert((ordinal, reason));
                    }
                }
            }
        }
        loaded.push((ordinal, parsed));
    }
    // A handle that depends on the scope's earlier reads matters only where no other flipped
    // Seq-align fails the prefetch.
    if let Some((ordinal, reason)) = uncertain.filter(|_| !fails) {
        return Err(parsed_refusal(
            ordinal,
            &format!("{reason}, and has a Seq-align on the query's minus strand and its plus strand, whose length NCBI BLAST+ 2.17.0's prefetch asks for"),
            program,
        ));
    }
    let mut text = prefetch_messages(db, run, loaded, program)?;
    if let Some((ordinal, id)) = &failed_lookup {
        if !text.is_empty() {
            return Err(parsed_refusal(
                *ordinal,
                &format!(
                    "has the Seq-align id {id}, whose lookup in the database fails where NCBI BLAST+ 2.17.0's prefetch loads the batch's subjects, in an order that LOSAT does not reproduce: the taxonomy database's texts of the other subjects of the batch, which the prefetch may load before or not at all, are not known"
                ),
                program,
            ));
        }
    }
    if failed_lookup.is_some() || fails {
        text.extend_from_slice(
            format!(
                "Warning: [{}] Error pre-fetching sequence data\n",
                program.to_ascii_lowercase()
            )
            .as_bytes(),
        );
    }
    Ok(text)
}

/// The volume of `db` that holds the taxid `last` and another positive taxid of the same slot of
/// NCBI's taxonomy cache of a volume (`taxid & 255`, `TaxdbRun::cache_collision`), with that
/// taxid; `None` when no volume does. Every definition line's taxid of every sequence, parsed
/// Seq-ids or not, in one pass over the headers; a header whose taxids do not decode is skipped.
fn taxonomy_cache_collision(db: &SeqDb, last: i32) -> Option<(usize, i32)> {
    let slot = tax_cache_slot;
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

/// The end of the reason of a refusal by [`prefetch_reads`] where no report reads the
/// subject's Bioseq.
const PREFETCH_READS_IT: &str = "and NCBI BLAST+ 2.17.0's prefetch reads it by that id though the report reads nothing of it (outfmt 8, 9 and 12 with PRE_FETCH_SEQS_LIMIT=2147483647, a tabular field list that reads nothing of the subject)";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1340-1355
// ```c++
//         ITERATE(TDeflines, iter, orig_deflines->Get()) {
//             if ((seqid.NotEmpty() && s_SeqDB_SeqIdIn((**iter).GetSeqid(), *seqid)) ||
//                 (seqid_lcl.NotEmpty() && s_SeqDB_SeqIdIn((**iter).GetSeqid(), *seqid_lcl))) {
//                 filt_dl = *iter;
//                 break;
//             }
//         }
//
//         if (filt_dl.Empty()) {
//             NCBI_THROW(CSeqDBException, eArgErr,
//                        "Error: oid headers do not contain target gi/seq_id.");
// ```
// The prefetch rethrows a `CSeqDBException` (`BlastFormatter_PreFetchSequenceData`, quoted at
// `ParsedSubject::prefetched`): NCBI stops with `BLAST Database error` (exit 2; probe of line DBID
// P3, `pdb|1ABC|A` and `pdb|1abc|A` in the two databases of `ids/pshare_cd`).
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1465-1466
// ```c++
//     list< CRef<CSeqdesc> > tax =
//         x_GetTaxonomy(oid, target_gi, target_seq_id);
// ```
// The blob that the loader reads for the Seq-align id is that of the sequence the scope gets
// (`scope_bioseq`), whose taxids its volume looks up: another sequence's where the lookup gives
// that one first and the loader keeps it (review 2 of line DBID pssmread, H3: `gi|777` of
// `lcl|777 ^A gi|777` after `gi|777|gb|YYY00001.1|` of taxid 99990 and a taxonomy database with a
// wrong magic number writes its text at the prefetch).
/// What NCBI's prefetch reads for the database sequence `ordinal` with parsed Seq-ids `subject`,
/// which it loads by its Seq-align id where no report may read its Bioseq (psiblast's traceback
/// hit lists beyond `-max_target_seqs`, the reports of `ParsedShown::prefetched` that read
/// nothing of it; `ParsedSubject::prefetched`): `Ok(Some(blob))` the sequence, by ordinal id, whose
/// taxids it looks up (`ordinal` itself, `ParsedSubject::volume_taxids`, or another sequence that
/// the lookup gives first and the data loader keeps, [`scope_bioseq`]), `Ok(None)` nothing whose
/// lookups write a text, and `Err` with the reason of LOSAT's refusal: a lookup that fails (NCBI's
/// prefetch writes `Error pre-fetching sequence data`, [`batch_prefetch`]; the reports that read
/// it stop), another sequence that NCBI's data loader can take for it whose definition lines hold
/// no id that `s_SeqDB_SeqIdIn` takes for the Seq-align id (`DeflineSet::bioseq_defline`:
/// `GetBioseq` throws and NCBI stops), and with `texts` (`taxdb.bti` on the search path) where the
/// taxonomy database's texts can differ between what the prefetch can read (`lookup_texts_differ`:
/// the subject or another sequence that also holds the id, as the blobs read before decide; the
/// sequence the lookup gives or none where another holder's Bioseq can stand for it; the same
/// record in another volume). Without `texts` nothing is looked up (`Ok(None)` but for the subject
/// itself). The other sequence's definition lines are read from `db`.
pub(crate) fn prefetch_reads(
    db: &SeqDb,
    ordinal: u32,
    subject: &crate::report::db_subject::ParsedSubject,
    texts: bool,
    program: &str,
) -> Result<std::result::Result<Option<u32>, String>> {
    use crate::report::db_subject::Prefetched;
    let (other, or_itself, reason) = match subject.prefetched() {
        Prefetched::Itself => return Ok(Ok(Some(ordinal))),
        Prefetched::Nothing => return Ok(Ok(None)),
        Prefetched::Failed(reason) => return Ok(Err(reason)),
        Prefetched::Another {
            other: None,
            reason,
            ..
        } => return Ok(if texts { Err(reason) } else { Ok(None) }),
        Prefetched::Another {
            other: Some(other),
            or_itself,
            reason,
        } => (other, or_itself, reason),
    };
    let set = db
        .defline_set(other)
        .map_err(|error| database_error(error, program))?;
    // `x_GetOid` keeps the OID only where one of its ids compares equal to the Seq-align id
    // (`scope_bioseq`), else the prefetch reads nothing; the subject that another sequence also
    // holds can be read as that sequence through the ids of a blob read before (`or_itself`).
    let read = other_bioseq(&set, other, subject);
    let throws = match read {
        Some(read) if !or_itself => read.throws,
        _ => or_itself && set.bioseq_defline(subject.alignment_id()).is_none(),
    };
    if throws {
        return Ok(Err(format!(
            "{reason}; NCBI BLAST+ 2.17.0 can read that sequence for it, whose definition lines do not hold the Seq-align id, and then stops (oid headers do not contain target gi/seq_id)"
        )));
    }
    if !texts {
        return Ok(Ok(None));
    }
    let other_taxids = defline_set_taxids(&set);
    if or_itself {
        let own = subject.volume_taxids();
        if lookup_texts_differ(db, &[&own, &other_taxids]) {
            return Ok(Err(reason));
        }
        return Ok(Ok(Some(ordinal)));
    }
    match read {
        // Another holder's Bioseq that the scope read before stands for it without a load.
        Some(_) if subject.other_holder().is_some() => {
            if lookup_texts_differ(db, &[&other_taxids, &[]]) {
                return Ok(Err(reason));
            }
            Ok(Ok(None))
        }
        Some(_) => Ok(Ok(Some(other))),
        None => Ok(Ok(None)),
    }
}

/// LOSAT's refusal of the reported sequence `ordinal` with parsed Seq-ids `subject` where the
/// batches' prefetch reads it though no report reads its Bioseq (`ParsedShown::prefetched`), with
/// or without `taxdb.bti` on the search path: another sequence whose definition lines do not
/// hold its Seq-align id, which NCBI reads for it and stops ([`prefetch_reads`] without
/// texts; the texts are `prefetch_messages`'s). A lookup that fails is not refused: the batch's
/// prefetch writes its warning ([`batch_prefetch`]).
pub(crate) fn check_prefetched_parsed(
    db: &SeqDb,
    subject: &crate::report::db_subject::ParsedSubject,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if subject.prefetch_fails() {
        return Ok(());
    }
    match prefetch_reads(db, ordinal, subject, false, program)? {
        Ok(_) => Ok(()),
        Err(reason) => Err(parsed_refusal(
            ordinal,
            &format!("{reason}, {PREFETCH_READS_IT}"),
            program,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1113-1134
// ```c++
// 	CScope::TBioseqHandles bhs = scope->GetBioseqHandles(ids);
//     ...
//     ITERATE(CScope::TBioseqHandles, it, bhs) {
//     	CRef<CDelta_seq> seq(new CDelta_seq);
// 	    CSeq_interval& interval = seq->SetLoc().SetInt();
// 	    interval.SetId(*SerialClone(*it->GetAccessSeq_id_Handle().GetSeqId()));
// 	    if (ranges[i].GetFrom() > ranges[i].GetToOpen()) {
// 	    	TSeqPos length = it->GetBioseqLength();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:679-687
// ```c++
//             CConstRef<CSeq_id> subj(&(*aln)->GetSeq_id(kSubjRow));
//             TSeqRange subj_range((*aln)->GetSeqRange(kSubjRow));
//             if ((*aln)->GetSeqStrand(kQueryRow) == eNa_strand_minus &&
//                 (*aln)->GetSeqStrand(kSubjRow) == eNa_strand_plus) {
//                 TSeqRange r(subj_range);
//                 // flag the range as needed to be flipped once the sequence
//                 // length is known
// ```
// A Seq-align on the query's minus strand and the subject's plus strand (BLASTN's minus-strand
// HSPs, TBLASTX's of a minus query frame and a plus subject frame: the strand of the first
// Std-seg) asks the handle for the length; a null handle throws, which the prefetch catches with
// `Error pre-fetching sequence data` (probe of line DBID P3: `lcl|777` before `gi|777` in
// `ids/nrefuse_v5`, TBLASTX). The protein side of TBLASTN and BLASTX has no strand.
// Where the data loader reads a blob for the Seq-align id, the scope's handle is that blob's
// Bioseq or none, whatever else holds the id: `gi|777` of `lcl|777 ^A gi|777` warns in every run of
// review 2 of line DBID pssmread (`nunf_ah_a nunf_ah_c8`, `nx_unf_0y_4 nu_a_4`, `nz_unf_0y
// nz_unf_a`: BLASTN's three tasks and TBLASTX, versions 4 and 5) though another sequence holds it.
// Where it reads none, a Bioseq that the scope holds from an earlier read can stand for it: for
// `gb|AAA00001.1|` after `lcl|AAA00001.1` with `gb|ZZZ00001.1| ^A emb|AAA00001.1|` also held
// (`nah4`), NCBI warns in some runs and not in others of the same database.
/// Whether NCBI's prefetch gets no Bioseq of the database sequence `ordinal` with parsed Seq-ids
/// `subject` by its Seq-align id (a null handle of `GetBioseqHandles`, [`scope_bioseq`]): `Err`
/// with the reason of LOSAT's refusal where that depends on what the scope read before (the data
/// loader reads no blob for the id and another sequence holds it, [`ScopeBioseq::instead`]) or
/// where `GetBioseq` throws (NCBI's prefetch stops with the `CSeqDBException`). Where the batch has
/// a Seq-align of it on the query's minus strand and the subject's plus strand, a null handle
/// fails NCBI's prefetch ([`batch_prefetch`]; review 2 of line DBID pssmread, H5).
fn prefetch_gets_no_bioseq(
    db: &SeqDb,
    ordinal: u32,
    subject: &crate::report::db_subject::ParsedSubject,
    program: &str,
) -> Result<std::result::Result<bool, String>> {
    let read = scope_bioseq(db, ordinal, subject, program)?;
    let id = || String::from_utf8_lossy(&subject.alignment_fasta()).into_owned();
    if let Some(other) = read.blob.filter(|_| read.throws) {
        return Ok(Err(format!(
            "has the Seq-align id {}, by which NCBI BLAST+'s lookup gives database sequence {other}, whose definition lines do not hold that id: NCBI BLAST+ 2.17.0's prefetch stops there (oid headers do not contain target gi/seq_id)",
            id()
        )));
    }
    match read.instead {
        Some(holder) if read.blob.is_none() => Ok(Err(format!(
            "has the Seq-align id {}, by which NCBI BLAST+'s scope gets no sequence unless it read database sequence {holder}, which also holds that id, before: NCBI BLAST+ 2.17.0's prefetch of a minus-strand alignment then writes Error pre-fetching sequence data or not, as its reads decide",
            id()
        ))),
        _ => Ok(Ok(!read.bioseq)),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1129-1140
// ```c++
//     for(TBDLLConstIter iter = dl.begin(); iter != dl.end(); iter ++) {
//         TTaxId taxid = ZERO_TAX_ID;
//
//         if ((*iter)->CanGetTaxid()) {
//             taxid = (*iter)->GetTaxid();
//         }
//         if (taxid <= ZERO_TAX_ID) {
//             continue;
//         }
// ```
/// The lookups of the volume's taxonomy cache where NCBI loads the blob of the sequence `blob` for a
/// read of the subject `ordinal` by its Seq-align id (`CSeqDBVol::x_GetTaxonomy`, of the definition
/// lines of `x_GetTaxDefline`): their positive taxids in that order, a taxid of several lines once
/// for each ([`TaxdbRun::load`]). The subject's own (`ParsedSubject::volume_lookups`, the line
/// that holds the Seq-align id first) where it has parsed Seq-ids (`parsed`) and the read loads its
/// own blob, another sequence's lines where the read loads that one (`DeflineSet::taxonomy_lookups`
/// of the subject's Seq-align id, which the data loader asks that blob for), and for a subject
/// without parsed Seq-ids, which has one definition line (`stored_title`) and is read by its
/// ordinal id, the first positive taxid of that line. The one producer of these for the prefetch
/// ([`prefetch_messages`]) and psiblast's PSSM ([`pssm_load_texts`]).
pub(crate) fn blob_taxids(
    db: &SeqDb,
    ordinal: u32,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    blob: u32,
    program: &str,
) -> Result<Vec<i32>> {
    Ok(match parsed {
        Some(subject) if blob == ordinal => subject.volume_lookups(),
        Some(subject) => db
            .defline_set(blob)
            .map_err(|error| database_error(error, program))?
            .taxonomy_lookups(subject.alignment_id()),
        None => db
            .defline_taxids(blob)
            .map_err(|error| database_error(error, program))?
            .first()
            .copied()
            .filter(|&taxid| taxid > 0)
            .into_iter()
            .collect(),
    })
}

/// The positive taxids of the definition lines `set`, each once, in their order: the taxids that
/// NCBI looks up where it loads that sequence's blob (`DeflineSet::taxonomy_lookups` without the
/// repeats and the place of the line that holds the id), for the checks that ask which taxids a
/// load can look up.
pub(crate) fn defline_set_taxids(set: &crate::seqid::DeflineSet) -> Vec<i32> {
    let mut taxids: Vec<i32> = Vec::new();
    for defline in &set.deflines {
        if let Some(taxid) = defline.taxid.filter(|&taxid| taxid > 0) {
            if !taxids.contains(&taxid) {
                taxids.push(taxid);
            }
        }
    }
    taxids
}

/// What NCBI's scope gets by the Seq-align id of a database sequence with parsed Seq-ids
/// ([`scope_bioseq`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ScopeBioseq {
    /// The sequence, by ordinal id, whose blob NCBI's data loader reads for the id: the subject's
    /// own where the lookup gives it first (or the same record), another sequence that the lookup
    /// gives first where `x_GetOid` keeps it (one of its ids compares equal to the Seq-align id,
    /// [`loader_keeps`]); `None` where it reads none (the lookup gives nothing, fails, or gives a
    /// sequence that `x_GetOid` does not keep).
    pub(crate) blob: Option<u32>,
    /// Whether the scope gets the Bioseq of `blob` (a handle): the definition line that
    /// `GetBioseqNoData`'s filter keeps holds an id that compares equal to the Seq-align id (the
    /// Bioseq's ids are that line's). Else a null handle.
    pub(crate) bioseq: bool,
    /// Whether `GetBioseqNoData`'s filter keeps no definition line of `blob` (no line holds an id
    /// that `s_SeqDB_SeqIdIn` takes for the Seq-align id or its local form): it throws a
    /// `CSeqDBException` (`oid headers do not contain target gi/seq_id`) before the taxonomy of
    /// the blob is looked up.
    pub(crate) throws: bool,
    /// Another sequence, by ordinal id, that holds an id the scope takes for the Seq-align id
    /// besides `blob` (`ParsedSubject::other_holder`): where the scope holds its Bioseq from an
    /// earlier read, it can get that one instead (psiblast's PSSM refuses every such read,
    /// `ParsedSubject::pssm_read`; the prefetch's handle where `blob` is none,
    /// `prefetch_gets_no_bioseq`).
    pub(crate) instead: Option<u32>,
}

impl ScopeBioseq {
    /// The sequence whose taxids `GetBioseqNoData` looks up (`x_GetTaxonomy`, after the filter):
    /// `blob` but where the filter throws (review 2 of line DBID pssmread, M1).
    pub(crate) fn taxonomy(&self) -> Option<u32> {
        self.blob.filter(|_| !self.throws)
    }

    /// The sequence whose letters the scope gives for the id (`blob` where it gets its Bioseq).
    pub(crate) fn letters(&self) -> Option<u32> {
        self.blob.filter(|_| self.bioseq)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1340-1356 and
// 1437-1439
// ```c++
//         ITERATE(TDeflines, iter, orig_deflines->Get()) {
//             if ((seqid.NotEmpty() && s_SeqDB_SeqIdIn((**iter).GetSeqid(), *seqid)) ||
//                 (seqid_lcl.NotEmpty() && s_SeqDB_SeqIdIn((**iter).GetSeqid(), *seqid_lcl))) {
//                 filt_dl = *iter;
//                 break;
//             }
//         }
//
//         if (filt_dl.Empty()) {
//             NCBI_THROW(CSeqDBException, eArgErr,
//                        "Error: oid headers do not contain target gi/seq_id.");
//     ...
//     seqids = defline->GetSeqid();
//     ...
//     bioseq->SetId().swap(seqids);
// ```
// NCBI's scope gets a Seq-align id's Bioseq from the data loader (`x_GetOid`, quoted at
// `ParsedSubject::found_instead` and [`loader_keeps`]): the blob of the first sequence the
// lookup gives where one of its ids compares equal to the id, whose Bioseq `GetBioseqNoData`
// builds with the ids of the first definition line that holds the id or its local form; the scope
// finds the Bioseq where those ids hold the id. So `gnl|DB|SEQA` reads a sequence whose second
// line holds it (`gnl|DB|OTHER ^A gnl|DB|SEQA`), not the formatters' shown line (H2), and
// `gi|777` reads none of a sequence `lcl|777 ^A gi|777` (the local form's line, review 1 H3 (b));
// probes of review 2 of line DBID pssmread (`R/pssmread/review2`: `xg_0y`, `nz_unf_0y`,
// `nu_0y`).
/// What NCBI's scope gets by the Seq-align id of the database sequence `ordinal` with parsed
/// Seq-ids `subject` ([`ScopeBioseq`]), from what the lookup finds (`ParsedSubject::found`; a
/// giim id, which NCBI's handles never take again, reads its own blob and gets no Bioseq): the one
/// producer of psiblast's PSSM reads (`ListedTitles::pssm_reads`) and of the prefetch's handles
/// ([`prefetch_gets_no_bioseq`], [`prefetch_reads`]). A subject whose definition line the
/// formatters do not find (`unfetchable`) is read by the scope as any other.
pub(crate) fn scope_bioseq(
    db: &SeqDb,
    ordinal: u32,
    subject: &crate::report::db_subject::ParsedSubject,
    program: &str,
) -> Result<ScopeBioseq> {
    use crate::report::db_subject::FoundInstead;
    let instead = subject.other_holder();
    let none = ScopeBioseq {
        blob: None,
        bioseq: false,
        throws: false,
        instead,
    };
    match subject.found() {
        Some(FoundInstead::Sequence(other)) => {
            let set = db
                .defline_set(*other)
                .map_err(|error| database_error(error, program))?;
            Ok(other_bioseq(&set, *other, subject).unwrap_or(none))
        }
        Some(FoundInstead::Nothing | FoundInstead::ThrowsOnId(_) | FoundInstead::Failed(_)) => {
            Ok(none)
        }
        None
        | Some(
            FoundInstead::AlsoHeld(_) | FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(..),
        ) => Ok(ScopeBioseq {
            blob: Some(ordinal),
            bioseq: subject.bioseq_holds_alignment_id(),
            throws: false,
            instead,
        }),
    }
}

/// [`scope_bioseq`] of the sequence `other`, of the definition lines `set`, that the lookup of
/// the Seq-align id of `subject` gives first: `None` where `x_GetOid` does not keep it.
fn other_bioseq(
    set: &crate::seqid::DeflineSet,
    other: u32,
    subject: &crate::report::db_subject::ParsedSubject,
) -> Option<ScopeBioseq> {
    if !loader_keeps(set, subject) {
        return None;
    }
    let id = subject.alignment_id();
    let kept = set.bioseq_defline(id);
    Some(ScopeBioseq {
        blob: Some(other),
        bioseq: kept.is_some_and(|line| {
            set.deflines[line].ids.iter().any(|held| {
                crate::seqid::defline_set::compare(held, id) == crate::seqid::defline_set::Sic::Yes
            })
        }),
        throws: kept.is_none(),
        instead: subject.other_holder(),
    })
}

/// Whether NCBI's data loader keeps the sequence of the definition lines `set` that the lookup of
/// the Seq-align id of `subject` gives (`CBlastDbDataLoader::x_GetOid`, bdbloader.cpp): one of its
/// ids compares equal to that id.
fn loader_keeps(
    set: &crate::seqid::DeflineSet,
    subject: &crate::report::db_subject::ParsedSubject,
) -> bool {
    set.deflines
        .iter()
        .flat_map(|defline| &defline.ids)
        .any(|id| {
            crate::seqid::defline_set::compare(id, subject.alignment_id())
                == crate::seqid::defline_set::Sic::Yes
        })
}

/// Whether the batches of a `-db` search prefetch their reported subjects though the first format
/// `outfmt` writes the Seq-align ids alone: outfmt 8, 9 and 12 with
/// `PRE_FETCH_SEQS_LIMIT=2147483647` ([`prefetches_subjects`]; BLASTP's block search declines
/// such a run, as it does a table of fields that read nothing of the subject on a database with
/// parsed Seq-ids, whose definition lines it never takes).
pub(crate) fn seq_annots_prefetch(outfmt: Option<i32>) -> bool {
    matches!(outfmt, Some(8 | 9 | 12)) && prefetches_subjects(outfmt)
}

/// Whether the first format of a call is one whose batches prefetch their reported subjects
/// and that shows no taxonomy of them (`prefetch_messages`): outfmt 0, 5, 6, 7 and 17, and with
/// `PRE_FETCH_SEQS_LIMIT=2147483647` the Seq-annots 8, 9 and 12 too (`s_PreFetchSeqs`, quoted at
/// [`prefetch_may_be_off`], returns before its test of the format). The query warnings, which
/// the batch's messages go with, are written with the first format (the Seq-annots write only
/// the reader's); the XML2 and JSON formats write theirs at each `Hit` (`TaxdbPrinter`), and the
/// archive, which NCBI's apps but psiblast write without calling the prefetch, prefetches nothing
/// (psiblast_app.cpp:262-264 prefetches before its test of the archive; LOSAT's PSIBLAST refuses
/// outfmt 8, 9, 11 and 12).
pub(crate) fn prefetches_subjects(outfmt: Option<i32>) -> bool {
    match outfmt {
        Some(0 | 5 | 6 | 7 | 17) => true,
        Some(8 | 9 | 12) => prefetch_seqs_limit() == Some(i32::MAX),
        _ => false,
    }
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

/// Whether [`read_reported_taxonomy`] (`names`) or [`prefetch_messages`] (not `names`) could
/// refuse for some sequence of the database, which a run that writes its batches one by one
/// cannot take back (OF5b BLASTX, DECISIONS 32): with `names` a scientific name NCBI's XML
/// stream cannot write (the XML2 and JSON formats write the names; the prefetch of outfmt 0, 5,
/// 6, 7 and 17 asks for none), and in both a record past the last mapped page of `taxdb.btd`,
/// or a header or the taxonomy database that cannot be read. A superset of what a search
/// reports: limits (`-gilist`, `-taxidlist`) and the hits themselves are ignored. `false`
/// without `taxdb.bti` on the search path and for a database without taxids. The refusal of
/// NCBI's cache (`TaxdbRun::cache_collision`) is not asked here: it is a property of the run
/// (`TaxdbRun::prefetch_refused`), which BLASTX and BLASTP's block search ask apart.
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
    let taxdb_first = taxdb
        .entry_count()
        .is_ok_and(|entries| entries <= TAXDB_ENTRIES_PER_SEQUENCE * sequences);
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
/// `names` + `names::check_text` code `read_reported_taxonomy` uses (`names::check_text` only
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
            // One definition line per sequence (`read_reported_taxonomy`).
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
                crate::report::db_subject::ParsedSubject::new(set).unwrap()
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
    fn a_reported_sequence_whose_seq_align_id_finds_another_record_first() {
        // tests/db_search/ids/make_ids.py (line DBID N1): the second database of an alias holds
        // the record of the first, which NCBI shows for it (`pdup_ab`), or the record with another
        // title, which it shows in its place (`pdup_ac`).
        let cwd = empty_cwd("found_instead");
        for (name, want) in [
            ("ids/pdup_ab", None),
            (
                "ids/pdup_ac",
                Some("lookup gives database sequence 0 first"),
            ),
            (
                "ids/pfound_v5",
                Some("lookup gives database sequence 0 first"),
            ),
        ] {
            let db = db_search(name, None, &cwd);
            let headers = HeaderReader::new(&db);
            let ReportedHeader::Parsed(subject) = read_header(&headers, 1, 0, "BLASTP").unwrap()
            else {
                panic!("{name}: a sequence with parsed Seq-ids");
            };
            match (subject.fetched(false, true), want) {
                (Ok(()), None) => {}
                (Err(error), Some(text)) => assert!(error.contains(text), "{name}: {error}"),
                (got, _) => panic!("{name}: {got:?}"),
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn the_same_record_with_other_taxids_is_shown_apart_by_the_reports_of_taxids_only() {
        // tests/db_search/ids/make_ids.py (line DBID N1b): `pdup_at` holds `gb|EEE00005.1|`
        // without a taxid then with taxid 9606, `pdup_tu` the other way, `pdup_no` twice with
        // taxid 99990 (the last entry of every variant of tests/seqdb/taxdb, whose lookup writes
        // the offset error of `offset_last_past_end`), `pdupe_at` an empty title without a taxid
        // then with 9606; `taxdb`, the variant on BLASTDB. `without`, `with`: the refusal without
        // the reports of taxids (`fetched(false)`) and with them.
        let cwd = empty_cwd("other_taxids");
        let first = |n: u32| format!("lookup gives database sequence {n} first");
        let taxids = |n: u32| {
            format!("database sequence {n} also holds with the same definition lines and letters but other taxids (")
        };
        let title = |n: u32| {
            format!("database sequence {n} also holds with the same definition lines and letters but other taxids, from which")
        };
        let texts = |n: u32| {
            format!("database sequence {n} also holds with the same definition lines and letters: where")
        };
        let both = |text: String| (Some(text.clone()), Some(text));
        for (name, taxdb, ordinal, (without, with)) in [
            ("ids/pdup_at", None, 0, (None, Some(taxids(1)))),
            ("ids/pdup_at", None, 1, (None, Some(taxids(0)))),
            ("ids/pdup_tu", None, 0, (None, Some(taxids(1)))),
            ("ids/pdup_tu", None, 1, (None, Some(taxids(0)))),
            ("ids/pdup_at", Some("good"), 1, (None, Some(taxids(0)))),
            // Texts of opening at the run's first lookup, which only the copy with a taxid makes.
            ("ids/pdup_at", Some("bad_magic"), 0, both(texts(1))),
            ("ids/pdup_at", Some("bad_magic"), 1, both(texts(0))),
            ("ids/pdup_tu", Some("count_larger"), 1, both(texts(0))),
            // The offset error of a taxid that neither copy holds.
            (
                "ids/pdup_at",
                Some("offset_last_past_end"),
                1,
                (None, Some(taxids(0))),
            ),
            // The same taxids: the texts of opening come at the same lookup, the offset error
            // once per volume.
            ("ids/pdup_no", None, 1, (None, None)),
            ("ids/pdup_no", Some("good"), 1, (None, None)),
            ("ids/pdup_no", Some("bad_magic"), 1, (None, None)),
            ("ids/pdup_no", Some("count_larger"), 0, (None, None)),
            (
                "ids/pdup_no",
                Some("offset_last_past_end"),
                1,
                both(texts(0)),
            ),
            (
                "ids/pdup_no",
                Some("offset_last_past_end"),
                0,
                both(texts(1)),
            ),
            // The same record without taxids looks nothing up.
            ("ids/pdup_ab", Some("offset_last_past_end"), 1, (None, None)),
            // An empty title made from the taxids where the taxonomy database is found, which the
            // reports of generated titles refuse (below).
            ("ids/pdupe_at", None, 0, (None, Some(taxids(1)))),
            ("ids/pdupe_at", Some("good"), 0, (None, Some(taxids(1)))),
            ("ids/pdupe_at", Some("good"), 1, (None, Some(taxids(0)))),
            // A taxonomy database that names no taxid: NCBI makes no title from them.
            ("ids/pdupe_at", Some("no_btd"), 0, (None, Some(taxids(1)))),
            (
                "ids/pdupe_at",
                Some("empty_fields"),
                0,
                (None, Some(taxids(1))),
            ),
            // One that writes a text at the run's first lookup, which only the copy with a taxid
            // makes.
            ("ids/pdupe_at", Some("bad_magic"), 0, both(texts(1))),
            // Another title is another record whatever the taxids.
            ("ids/pdup_ac", None, 1, both(first(0))),
        ] {
            let taxdb = taxdb.map(|variant| root().join("seqdb/taxdb").join(variant));
            let db = db_search(name, taxdb.as_deref(), &cwd);
            let headers = HeaderReader::new(&db);
            let ReportedHeader::Parsed(subject) =
                read_header(&headers, ordinal, 0, "BLASTP").unwrap()
            else {
                panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
            };
            for (taxonomy, want) in [(false, &without), (true, &with)] {
                match (subject.fetched(taxonomy, true), want) {
                    (Ok(()), None) => {}
                    (Err(error), Some(text)) => {
                        assert!(error.contains(text.as_str()), "{name} {ordinal}: {error}")
                    }
                    (got, _) => panic!("{name} {ordinal} {taxdb:?} {taxonomy}: {got:?}"),
                }
            }
            // Line DBID N1b L1: outfmt 0 and the tabular word of a replaced id
            // (`ParsedSubject::generated_title_fetched`).
            let shown =
                name == "ids/pdupe_at" && taxdb.as_deref().is_some_and(|dir| dir.ends_with("good"));
            match subject.generated_title_fetched() {
                Ok(()) => assert!(!shown, "{name} {ordinal} {taxdb:?}"),
                Err(error) => assert!(
                    shown && error.contains(&title(1 - ordinal)),
                    "{name} {ordinal} {taxdb:?}: {error}"
                ),
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn a_seq_align_id_of_any_kind_is_left_to_the_lookup() {
        // tests/db_search/ids/make_look_fasta.py (line DBID N2b): `None`, NCBI's formatters find the
        // sequence itself; else the text of the refusal.
        let cwd = empty_cwd("any_kind");
        let other = |n: u32| format!("lookup gives database sequence {n} first");
        let none = "formatters do not find the sequence".to_string();
        let held = |n: u32| format!("which database sequence {n} also holds");
        for (name, ordinal, want) in [
            ("ids/plook_v5", 1, None),
            ("ids/plook_v5", 2, None),
            ("ids/plook_v5", 3, None),
            ("ids/plook_v5", 4, Some(none.clone())),
            ("ids/plook_v5", 5, None),
            ("ids/plook_v5", 6, Some(other(5))),
            ("ids/plook_v5", 8, Some(other(7))),
            ("ids/plook_v5", 9, None),
            ("ids/plook_v5", 10, None),
            ("ids/plook_v5", 11, Some(other(4))),
            ("ids/plook_v5", 12, Some(none.clone())),
            ("ids/plook_v5", 13, None),
            ("ids/plook_v4", 1, Some(none.clone())),
            ("ids/plook_v4", 2, Some(none.clone())),
            ("ids/plook_v4", 3, None),
            ("ids/plook_v4", 4, Some(none.clone())),
            ("ids/plook_v4", 5, Some(other(6))),
            ("ids/plook_v4", 9, None),
            ("ids/plook_v4", 11, Some(other(4))),
            ("ids/plook_v4", 12, Some(none.clone())),
            ("ids/plook_vz", 16, None),
            // Another sequence holds the Seq-align id (review 1 of N2b): NCBI's scope shows it
            // where both are reported.
            ("ids/pshare_v5", 0, Some(held(1))),
            ("ids/pshare_v5", 1, None),
            ("ids/pshare_v4", 0, Some(held(1))),
            ("ids/pshare_ab", 0, Some(held(1))),
            // Ids that differ in case alone (review 2 of N2b), and a later copy of a record whose
            // id a third sequence also holds.
            ("ids/pshare_v5", 4, Some(held(5))),
            ("ids/pshare_v4", 4, Some(held(5))),
            ("ids/pshare_cd", 0, Some(held(2))),
            ("ids/pshare_cd", 1, Some(held(3))),
            ("ids/pdup_abc", 0, Some(held(2))),
            ("ids/pdup_abc", 1, Some(held(2))),
            (
                "ids/pshare_v5",
                2,
                Some("whose lookup in the database fails".to_string()),
            ),
        ] {
            let db = db_search(name, None, &cwd);
            let headers = HeaderReader::new(&db);
            let ReportedHeader::Parsed(subject) =
                read_header(&headers, ordinal, 0, "BLASTP").unwrap()
            else {
                panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
            };
            assert!(!subject.unfetchable(), "{name} {ordinal}");
            match (subject.fetched(false, true), want) {
                (Ok(()), None) => {}
                (Err(error), Some(text)) => {
                    assert!(error.contains(&text), "{name} {ordinal}: {error}")
                }
                (got, _) => panic!("{name} {ordinal}: {got:?}"),
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn the_prefetch_reads_a_subject_by_its_seq_align_id_as_ncbis_loader_does() {
        // Line DBID P3 (probes of `blastdb-ids-20261021/R/intmaxpf`): where no report reads the
        // Bioseq, `prefetch_reads` without texts (no `taxdb.bti` on this search path).
        let cwd = empty_cwd("prefetch_reads");
        let stops = "and then stops (oid headers do not contain target gi/seq_id)";
        for (name, ordinal, want) in [
            // `pdb|1ABC|A` of pshare_c and `pdb|1abc|A` of pshare_d: the scope can read the
            // other, whose definition lines do not hold the id (NCBI exit 2).
            ("ids/pshare_cd", 0, Err(stops)),
            // `lcl|555`, which sequence 1 also holds in a second definition line: read as it or
            // as itself, nothing to refuse without texts.
            ("ids/pshare_v5", 0, Ok(None)),
            ("ids/pshare_v5", 1, Ok(Some(1))),
            // `gnl|BL_ORD_ID|0` stored with a text tag.
            (
                "ids/pshare_v5",
                2,
                Err("whose lookup in the database fails"),
            ),
            // An id that the lookup finds nowhere, and one that it finds first on sequence 5,
            // none of whose ids compares equal: nothing is read.
            ("ids/plook_v5", 4, Ok(None)),
            ("ids/plook_v5", 6, Ok(None)),
        ] {
            let db = db_search(name, None, &cwd);
            let headers = HeaderReader::new(&db);
            let ReportedHeader::Parsed(subject) =
                read_header(&headers, ordinal, 0, "BLASTP").unwrap()
            else {
                panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
            };
            match (
                prefetch_reads(&db, ordinal, &subject, false, "BLASTP").unwrap(),
                want,
            ) {
                (Ok(got), Ok(want)) => assert_eq!(got, want, "{name} {ordinal}"),
                (Err(reason), Err(text)) => {
                    assert!(reason.contains(text), "{name} {ordinal}: {reason}")
                }
                (got, _) => panic!("{name} {ordinal}: {got:?}"),
            }
            // A lookup that fails is the batch's warning (`batch_prefetch`), not a refusal.
            let check = check_prefetched_parsed(&db, &subject, ordinal, "BLASTP");
            assert_eq!(
                check.is_err(),
                want.is_err() && !subject.prefetch_fails(),
                "{name} {ordinal}"
            );
            assert_eq!(
                subject.prefetch_fails(),
                want == Err("whose lookup in the database fails"),
                "{name} {ordinal}"
            );
        }
        // `lcl|777` before `gi|777`: the formatters do not find the Bioseq, which the prefetch
        // reads as itself; it fails where it asks for the length (a minus-strand alignment).
        let search = SearchPath::new(&cwd, None, None).unwrap();
        let path = root().join("db_search/ids/nrefuse_v5");
        let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
        let headers = HeaderReader::new(&db);
        let ReportedHeader::Parsed(subject) = read_header(&headers, 0, 0, "TBLASTX").unwrap()
        else {
            panic!("nrefuse_v5 0: a sequence with parsed Seq-ids");
        };
        assert!(subject.unfetchable());
        assert_eq!(
            prefetch_reads(&db, 0, &subject, true, "TBLASTX").unwrap(),
            Ok(Some(0))
        );
        let run = TaxdbRun::default();
        let warning = b"Warning: [tblastx] Error pre-fetching sequence data\n".to_vec();
        let prefetch = |flipped: bool| {
            let subject = PrefetchSubject {
                ordinal: 0,
                parsed: Some(&*subject),
                flipped,
            };
            batch_prefetch(&db, &run, [subject], "TBLASTX").unwrap()
        };
        assert_eq!(prefetch(false), b"");
        assert_eq!(prefetch(true), warning);
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn a_batch_whose_prefetch_fails_writes_ncbis_warning_after_its_texts() {
        // tests/db_search/ids/make_ids.py (line DBID pfwarn): `ntag_t5` 0 `gnl|BL_ORD_ID|0`
        // stored with a text tag (its lookup fails), 1 `gb|FFF00092.1|`; `nlook_t4` 1 `bbs|123`,
        // which a version 4 database does not look up (no Bioseq), 0 `gb|AF000070.1|`; every
        // sequence with taxid 99990, whose lookup in `offset_last_past_end` writes the offset
        // error.
        let cwd = std::env::temp_dir().join(format!(
            "losat-batch-prefetch-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&cwd).unwrap();
        let taxdb = root().join("seqdb/taxdb/offset_last_past_end");
        let open = |name: &str, taxdb: Option<&Path>| {
            let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
            let path = root().join("db_search/ids").join(name);
            SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap()
        };
        let offset = "Error: [blastn] Error: Offset error at end of taxdb file.\n";
        let warning = "Warning: [blastn] Error pre-fetching sequence data\n";
        let text = |bytes: Vec<u8>| String::from_utf8(bytes).unwrap();
        for taxdb in [None, Some(taxdb.as_path())] {
            let db = open("ntag_t5", taxdb);
            let headers = HeaderReader::new(&db);
            let parsed = |ordinal: u32| match read_header(&headers, ordinal, 0, "BLASTN").unwrap() {
                ReportedHeader::Parsed(subject) => subject,
                _ => panic!("ntag_t5 {ordinal}: a sequence with parsed Seq-ids"),
            };
            let (failed, other) = (parsed(0), parsed(1));
            assert!(failed.prefetch_fails() && !other.prefetch_fails());
            let run = TaxdbRun::default();
            let alone = batch_prefetch(
                &db,
                &run,
                [PrefetchSubject {
                    ordinal: 0,
                    parsed: Some(&*failed),
                    flipped: false,
                }],
                "BLASTN",
            );
            assert_eq!(text(alone.unwrap()), warning);
            // With the other subject, whose load writes the offset error before the failed
            // lookup or not at all: refused where the taxonomy database writes it.
            let both = batch_prefetch(
                &db,
                &run,
                [
                    PrefetchSubject {
                        ordinal: 0,
                        parsed: Some(&*failed),
                        flipped: false,
                    },
                    PrefetchSubject {
                        ordinal: 1,
                        parsed: Some(&*other),
                        flipped: false,
                    },
                ],
                "BLASTN",
            );
            match taxdb {
                None => assert_eq!(text(both.unwrap()), warning),
                Some(_) => {
                    let refused = both.unwrap_err().to_string();
                    assert!(
                        refused.contains("in an order that LOSAT does not reproduce"),
                        "{refused}"
                    );
                }
            }
            // A later batch of the run, whose lookup of the other subject wrote its text before:
            // nothing to place.
            if taxdb.is_some() {
                let run = TaxdbRun::default();
                let first = batch_prefetch(
                    &db,
                    &run,
                    [PrefetchSubject {
                        ordinal: 1,
                        parsed: Some(&*other),
                        flipped: false,
                    }],
                    "BLASTN",
                );
                assert_eq!(text(first.unwrap()), offset);
                let later = batch_prefetch(
                    &db,
                    &run,
                    [
                        PrefetchSubject {
                            ordinal: 1,
                            parsed: Some(&*other),
                            flipped: false,
                        },
                        PrefetchSubject {
                            ordinal: 0,
                            parsed: Some(&*failed),
                            flipped: false,
                        },
                    ],
                    "BLASTN",
                );
                assert_eq!(text(later.unwrap()), warning);
            }
        }
        // A minus-strand Seq-align of a subject without a Bioseq fails after every load: the
        // texts of the batch, then the warning.
        let db = open("nlook_t4", Some(&taxdb));
        let headers = HeaderReader::new(&db);
        let parsed = |ordinal: u32| match read_header(&headers, ordinal, 0, "BLASTN").unwrap() {
            ReportedHeader::Parsed(subject) => subject,
            _ => panic!("nlook_t4 {ordinal}: a sequence with parsed Seq-ids"),
        };
        let (control, bbs) = (parsed(0), parsed(1));
        let run = TaxdbRun::default();
        let prefetch = batch_prefetch(
            &db,
            &run,
            [
                PrefetchSubject {
                    ordinal: 1,
                    parsed: Some(&*bbs),
                    flipped: true,
                },
                PrefetchSubject {
                    ordinal: 0,
                    parsed: Some(&*control),
                    flipped: false,
                },
            ],
            "BLASTN",
        );
        assert_eq!(text(prefetch.unwrap()), format!("{offset}{warning}"));
        let run = TaxdbRun::default();
        let plus = batch_prefetch(
            &db,
            &run,
            [PrefetchSubject {
                ordinal: 1,
                parsed: Some(&*bbs),
                flipped: false,
            }],
            "BLASTN",
        );
        assert_eq!(text(plus.unwrap()), "");
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn a_weak_match_is_refused_where_ncbi_may_read_the_reported_sequences_one_by_one() {
        // tests/db_search/ids/make_ids.py (line DBID N2c): `pweak_v5` and `pweak_v4` 0
        // `sp|P00041.1|NAMEW`, 1 `gb|FFF00040|`, 2 `pir|PIRW1|`, 3 `sp||NAMEX`, 4 `gb|FFF00041.1|`,
        // 5 `sp||NAMEY`, 6 `sp|P00042.1|NAMEY`, 7 `gb|FFF00043|`, 8 `gb|FFF00043.1|`, 9
        // `sp|P00044.1|NAMEZ`, 10 `sp||NAMEZ`. With the prefetch on, NCBI writes each sequence's
        // own report (in version 4 the lookup of 10 gives 9 first: refused either way, not listed).
        let cwd = empty_cwd("weak_match");
        // `None`: fetched; `Some(Some(n))`: another sequence n matches weakly; `Some(None)`: weak
        // keys that version 5 does not look up.
        let unlisted = Some(None);
        for (name, ordinal, want) in [
            ("ids/pweak_v5", 0, unlisted),
            ("ids/pweak_v5", 1, None),
            ("ids/pweak_v5", 2, unlisted),
            ("ids/pweak_v5", 3, unlisted),
            ("ids/pweak_v5", 4, None),
            ("ids/pweak_v5", 6, unlisted),
            ("ids/pweak_v5", 7, Some(Some(8))),
            ("ids/pweak_v5", 8, None),
            ("ids/pweak_v4", 0, None),
            ("ids/pweak_v4", 1, None),
            ("ids/pweak_v4", 2, None),
            ("ids/pweak_v4", 3, None),
            ("ids/pweak_v4", 5, Some(Some(6))),
            ("ids/pweak_v4", 6, Some(Some(5))),
            ("ids/pweak_v4", 7, Some(Some(8))),
            ("ids/pweak_v4", 9, Some(Some(10))),
        ] {
            let db = db_search(name, None, &cwd);
            for prefetch_off in [false, true] {
                let headers = HeaderReader {
                    prefetch_off,
                    ..HeaderReader::new(&db)
                };
                let ReportedHeader::Parsed(subject) =
                    read_header(&headers, ordinal, 0, "BLASTP").unwrap()
                else {
                    panic!("{name} {ordinal}: no parsed Seq-ids");
                };
                // A batch whose prefetch is on loads the subject by its own id.
                assert!(subject.fetched(false, false).is_ok(), "{name} {ordinal}");
                match (subject.fetched(false, true), want.filter(|_| prefetch_off)) {
                    (Ok(()), None) => {}
                    (Err(error), Some(Some(other))) => assert!(
                        error.contains(&format!("match an id of database sequence {other} weakly")),
                        "{name} {ordinal}: {error}"
                    ),
                    (Err(error), Some(None)) => assert!(
                        error.contains("a version 5 database does not look up"),
                        "{name} {ordinal}: {error}"
                    ),
                    (got, want) => {
                        panic!("{name} {ordinal} {prefetch_off}: {got:?}, want {want:?}")
                    }
                }
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn a_stored_ordinal_id_of_another_number_is_not_fetched() {
        // tests/db_search/ids/make_ids.py (line DBID N2a): sequences whose only Seq-id is
        // `gnl|BL_ORD_ID|n`. With n the sequence's number in its volume it reads as a sequence
        // without parsed Seq-ids (in `pord_vx`, sequence 6 is the second of `pord_x`); another n is
        // a parsed Seq-id, n plus the first ordinal id of the volume, by which NCBI's formatters
        // find no sequence or another one.
        let cwd = empty_cwd("stored_ordinal_id");
        // `Ok(true)`: read as a sequence without parsed Seq-ids; `Ok(false)`: a parsed subject
        // that the reports fetch; `Err(id)`: refused, with its Seq-align id.
        for (name, ordinal, want) in [
            ("ids/pord_v5", 0, Ok(false)),
            ("ids/pord_v5", 2, Ok(true)),
            ("ids/pord_v5", 3, Err("gnl|BL_ORD_ID|6")),
            ("ids/pord_v5", 4, Err("gnl|BL_ORD_ID|1")),
            ("ids/pord_v4", 2, Ok(true)),
            ("ids/pord_v4", 3, Err("gnl|BL_ORD_ID|6")),
            ("ids/pord_vx", 3, Err("gnl|BL_ORD_ID|6")),
            ("ids/pord_vx", 5, Ok(false)),
            ("ids/pord_vx", 6, Ok(true)),
            ("ids/pord_vx", 7, Err("gnl|BL_ORD_ID|12")),
        ] {
            let db = db_search(name, None, &cwd);
            let headers = HeaderReader::new(&db);
            match (read_header(&headers, ordinal, 0, "BLASTP").unwrap(), want) {
                (ReportedHeader::Stored(_), Ok(true)) => {}
                (ReportedHeader::Parsed(subject), Ok(false)) => {
                    assert!(subject.fetched(false, true).is_ok(), "{name} {ordinal}")
                }
                (ReportedHeader::Parsed(subject), Err(id)) => {
                    let error = subject.fetched(false, true).unwrap_err();
                    assert!(
                        error.contains(&format!("the Seq-align id {id},")),
                        "{name} {ordinal}: {error}"
                    );
                }
                (got, want) => panic!("{name} {ordinal}: {got:?}, want {want:?}"),
            }
        }
        // The first sequence of `nord_v5` holds `gnl|BL_ORD_ID|3`: the database has parsed Seq-ids.
        let search = SearchPath::new(&cwd, None, None).unwrap();
        for name in ["nord_v5", "nord_v4"] {
            let path = root().join("db_search/ids").join(name);
            let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
            assert!(db.has_parsed_seqids().unwrap(), "{name}");
            assert!(matches!(
                db.read_header(1).unwrap(),
                crate::seqdb::SequenceHeader::Ordinal(_)
            ));
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn an_empty_stored_title_is_made_from_a_named_taxid_in_outfmt_0() {
        // tests/db_search/ids/make_ids.py (line DBID H3): `pbare_t` and `nbare_t` hold one
        // sequence without parsed Seq-ids, without a title, with taxid 9606, which the good
        // taxonomy database names and `no_btd` and `bad_magic` do not;
        // tests/fixtures/outfmt16/db/other_tax9606 has titles with that taxid.
        let cwd = empty_cwd("stored_taxonomy");
        for (name, seq_type, taxdb, want) in [
            (
                "db_search/ids/pbare_t",
                SeqType::Protein,
                Some("good"),
                true,
            ),
            (
                "db_search/ids/pbare_t",
                SeqType::Protein,
                Some("no_btd"),
                false,
            ),
            ("db_search/ids/pbare_t", SeqType::Protein, None, false),
            (
                "db_search/ids/nbare_t",
                SeqType::Nucleotide,
                Some("good"),
                true,
            ),
            (
                "db_search/ids/nbare_t",
                SeqType::Nucleotide,
                Some("bad_magic"),
                false,
            ),
            (
                "fixtures/outfmt16/db/other_tax9606",
                SeqType::Protein,
                Some("good"),
                false,
            ),
        ] {
            let taxdb = taxdb.map(variant);
            let search =
                SearchPath::new(&cwd, taxdb.as_deref().map(Path::as_os_str), None).unwrap();
            let path = root().join(name);
            let db = SeqDb::open_in(path.to_str().unwrap(), seq_type, &search).unwrap();
            let headers = HeaderReader::new(&db);
            let ReportedHeader::Stored(stored) = read_header(&headers, 0, 0, "BLASTP").unwrap()
            else {
                panic!("{name}: a sequence without parsed Seq-ids");
            };
            assert_eq!(stored.taxonomy_title, want, "{name} {taxdb:?}");
            let refused = check_shown_stored_title(stored.taxonomy_title, 0, "BLASTP");
            assert_eq!(refused.is_err(), want, "{name} {taxdb:?}");
            // TBLASTX's check before the search (`report_stream`) for outfmt 0 and for a report
            // that shows no title (outfmt 5 of a nucleotide database).
            if seq_type == SeqType::Nucleotide {
                let pool = crate::utils::threading::SearchPool::x_serial();
                let shown = |pairwise| ParsedShown {
                    pairwise,
                    ..Default::default()
                };
                let found = db.taxdb_found();
                let refusal = |pairwise| {
                    reports_may_refuse(
                        &db,
                        shown(pairwise),
                        false,
                        found,
                        false,
                        None,
                        &pool,
                        "TBLASTX",
                    )
                };
                assert_eq!(
                    refusal(true),
                    want.then_some(ReportRefusal::Headers),
                    "{taxdb:?}"
                );
                assert_eq!(refusal(false), None, "{taxdb:?}");
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
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

    /// `read_reported_taxonomy` (XML2 and JSON, line DBID L1d): NCBI fetches every reported
    /// subject, so a run whose volume has the cache collision is refused where a reported subject
    /// has a positive taxid (`nz_rep`: a1 99990, b1 150, c1 99990, d1 without a taxid, e1 -5;
    /// `nz_repo`: a1 to d1 without parsed Seq-ids), and not for d1 or e1 alone or for nothing
    /// reported; the taxonomy databases without the offset error read the same subjects.
    #[test]
    fn the_xml2_lookups_ask_for_the_cache_refusal_where_a_subject_has_a_taxid() {
        let cwd = empty_cwd("xml2_cache");
        let open = |taxdb: &str| {
            let search = SearchPath::new(&cwd, Some(variant(taxdb).as_os_str()), None).unwrap();
            let path = root().join("db_search/nz_rep");
            SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap()
        };
        let read = |db: &SeqDb, run: &TaxdbRun, ranks: &[usize]| {
            let parsed: Vec<_> = parsed_subjects(db)
                .into_iter()
                .map(|subject| Some(std::sync::Arc::new(subject)))
                .collect();
            let ordinals: Vec<u32> = (0..db.num_seqs()).collect();
            let reported: std::collections::BTreeSet<usize> = ranks.iter().copied().collect();
            read_reported_taxonomy(db, run, &ordinals, &reported, &parsed, "BLASTN").map(|_| ())
        };
        let db = open("offset_last_past_end");
        for ranks in [&[0][..], &[1], &[3, 2], &[0, 1, 2, 3]] {
            let run = TaxdbRun::default();
            let error = read(&db, &run, ranks).unwrap_err();
            assert!(error.is::<TaxonomyRefusal>(), "{ranks:?} {error:#}");
            let text = format!("{error:#}");
            assert!(text.contains("256-slot taxonomy cache"), "{text}");
            assert!(text.contains("LOSAT's BLASTN"), "{text}");
            // Refused before any lookup: the run's own handle is not open.
            assert!(run.lazy.lock().unwrap().is_none());
        }
        // Nothing reported, d1 (no taxid), e1 (taxid -5: a volume looks a positive taxid up only).
        for ranks in [&[][..], &[3], &[4], &[3, 4]] {
            assert!(read(&db, &TaxdbRun::default(), ranks).is_ok(), "{ranks:?}");
        }
        for taxdb in ["good", "bad_magic"] {
            assert!(
                read(&open(taxdb), &TaxdbRun::default(), &[0, 1, 2, 3, 4]).is_ok(),
                "{taxdb}"
            );
        }
        // nz_repo: a1 to d1 without parsed Seq-ids (the taxid of the definition line).
        let taxdb = variant("offset_last_past_end");
        let search = SearchPath::new(&cwd, Some(taxdb.as_os_str()), None).unwrap();
        let path = root().join("db_search/nz_repo");
        let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
        let ordinals: Vec<u32> = (0..db.num_seqs()).collect();
        for (ranks, refused) in [(&[0][..], true), (&[1, 3], true), (&[3], false)] {
            let reported: std::collections::BTreeSet<usize> = ranks.iter().copied().collect();
            let run = TaxdbRun::default();
            let result = read_reported_taxonomy(&db, &run, &ordinals, &reported, &[], "BLASTN");
            assert_eq!(result.is_err(), refused, "{ranks:?}");
            if let Err(error) = result {
                assert!(error.is::<TaxonomyRefusal>(), "{ranks:?} {error:#}");
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    /// `VolumeLookups::offset_slot`: the slot of the offset error's taxid L = 99990 (slot 150) of
    /// NCBI's cache of each volume, in the order of the loads (review 4 of line DBID pssmread):
    /// the first load of a sequence writes L's text at each lookup of L that misses the slot (twice
    /// for L, 150, L in one load), a read again of a sequence whose load looks L up writes nothing
    /// where the slot holds L whether NCBI's scope keeps it or loads it again, and is refused where
    /// loading it again would write the text; a load after a read again that may have put 150
    /// there is refused; other slots and other volumes are apart.
    #[test]
    fn the_slot_of_the_offset_errors_taxid_follows_the_loads() {
        const L: i32 = 99990;
        let fetch = |ordinal, taxid, again| SlotFetch {
            ordinal,
            taxid,
            again,
        };
        // H1 (cc_m): one load of L, 150, L; again: its own 150 is between.
        let mut record = VolumeLookups::default();
        assert_eq!(
            record.offset_slot(0, 0, &[L, 150, L], L),
            Ok(vec![true, false, true])
        );
        assert_eq!(
            record.offset_slot(0, 0, &[L, 150, L], L),
            Err(SlotUncertain {
                ordinal: 0,
                again: true,
                after: Some(fetch(0, 150, true)),
            })
        );
        // M1 (b) (cc_m2): A (L, 150) and B (L), first loads: two texts. Read again, A may put 150
        // there (Either), and B's read again is not certain.
        let mut record = VolumeLookups::default();
        assert_eq!(
            record.offset_slot(0, 0, &[L, 150], L),
            Ok(vec![true, false])
        );
        assert_eq!(record.offset_slot(0, 1, &[L], L), Ok(vec![true]));
        assert_eq!(
            record.offset_slot(0, 0, &[L, 150], L),
            Ok(vec![false, false])
        );
        assert_eq!(
            record.offset_slot(0, 1, &[L], L),
            Err(SlotUncertain {
                ordinal: 1,
                again: true,
                after: Some(fetch(0, 150, true)),
            })
        );
        // A first load after such a read again (C of 99990 + 256 * 2 is not L).
        let mut record = VolumeLookups::default();
        assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![true]));
        assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
        assert_eq!(record.offset_slot(0, 2, &[L], L), Ok(vec![true]));
        assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
        assert_eq!(
            record.offset_slot(0, 3, &[9606, L + 512, L], L),
            Ok(vec![false, false, true])
        );
        assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
        assert_eq!(
            record.offset_slot(0, 4, &[L], L),
            Err(SlotUncertain {
                ordinal: 4,
                again: false,
                after: Some(fetch(1, 150, true)),
            })
        );
        // M1 (a) (cc_ax): 150 never looked up: L's reads again write nothing.
        let mut record = VolumeLookups::default();
        assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![true]));
        for _ in 0..3 {
            assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![false]));
        }
        // cc_ac: A (L), C (150), then A again: refused, after C's fetch.
        let mut record = VolumeLookups::default();
        assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![true]));
        assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
        assert_eq!(
            record.offset_slot(0, 0, &[L], L),
            Err(SlotUncertain {
                ordinal: 0,
                again: true,
                after: Some(fetch(1, 150, false)),
            })
        );
        // C, then A in another volume (cc_c1 cc_a1), then both again; taxids of other slots.
        let mut record = VolumeLookups::default();
        assert_eq!(record.offset_slot(1, 0, &[150], L), Ok(vec![false]));
        assert_eq!(
            record.offset_slot(0, 1, &[L, 9606], L),
            Ok(vec![true, false])
        );
        assert_eq!(record.offset_slot(1, 0, &[150], L), Ok(vec![false]));
        assert_eq!(
            record.offset_slot(0, 1, &[L, 9606], L),
            Ok(vec![false, false])
        );
        assert_eq!(
            record.offset_slot(0, 5, &[9606, 151], L),
            Ok(vec![false, false])
        );
        assert!(!record.slot_loads.contains(&5));
    }

    /// `pssm_load_texts` on `tm_cache150` (st1 99990, st2 150, st3 99990, parsed Seq-ids, one
    /// volume) with the offset error of `offset_last_past_end`: the PSSM's loads in their order
    /// write NCBI's texts, and a read again that NCBI's scope decides is refused, naming the fetch
    /// of 150; with 151 (`tm_cache151`) the read again writes nothing.
    #[test]
    fn the_pssms_loads_write_the_offset_error_as_the_slot_holds_it() {
        let cwd = empty_cwd("pssm_slot");
        let offset = "Error: [psiblast] Error: Offset error at end of taxdb file.\n";
        let loads = |name: &str, ordinals: &[u32]| -> Result<String> {
            let db = db_search(name, Some(&variant("offset_last_past_end")), &cwd);
            let parsed = parsed_subjects(&db);
            let run = TaxdbRun::default();
            let mut text = Vec::new();
            for &ordinal in ordinals {
                let subject = &parsed[ordinal as usize];
                let lookups = blob_taxids(&db, ordinal, Some(subject), ordinal, "PSIBLAST")?;
                text.extend(pssm_load_texts(&db, &run, ordinal, &lookups, "PSIBLAST")?);
            }
            Ok(String::from_utf8(text).unwrap())
        };
        assert_eq!(loads("tm_cache150", &[0, 1, 2]).unwrap(), offset.repeat(2));
        assert_eq!(loads("tm_cache150", &[0, 2, 1]).unwrap(), offset);
        assert_eq!(loads("tm_cache150", &[0, 0, 2]).unwrap(), offset);
        let error = loads("tm_cache150", &[0, 1, 0]).unwrap_err();
        assert!(error.is::<TaxonomyRefusal>(), "{error:#}");
        let text = format!("{error:#}");
        for part in [
            "taxid 99990, the last entry of the taxonomy database",
            "reads database sequence 0 of the volume",
            "again after the run's fetch of database sequence 1 with taxid 150",
            "256-slot taxonomy cache",
            "LOSAT's PSIBLAST",
        ] {
            assert!(text.contains(part), "{part}: {text}");
        }
        assert_eq!(loads("tm_cache151", &[0, 1, 0, 2]).unwrap(), offset);
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    #[test]
    fn a_database_without_taxids_does_not_hold() {
        let db = database("sprot", Some(&variant("high_bytes_utf8")));
        assert!(db.taxdb_found());
        assert!(!taxdb_may_refuse(&db, true));
    }
}
