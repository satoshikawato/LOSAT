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

use crate::wording::{DATABASE_ERROR, ENGINE_ERROR};
use anyhow::{bail, Result};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::defline::SubjectWord;
use crate::report::pairwise::{DatabaseInfo, PairwiseDatabase};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
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
    // The errors of the alias files' counts depend on the list, and those of the headers come
    // after them (`restrict_database`).
    let mut db = SeqDb::open_search_in(name, seqtype, &path).map_err(DatabaseOpenError::Reader)?;
    db.fail_search_reads_with_oid_list();
    Ok(db)
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
// The list is resolved in the totals scan (`x_ScanTotals`, `CSeqDBGiListSet`: its errors and
// warnings), before the alias files' statistics keys are read.
/// Restricts the search's database `db` (opened as `name` by `open_database`) to the taxid or
/// seqid list of the command line (line DI-8, `DbFilterArgs::resolve`; none: the whole
/// database), as NCBI's `CSeqDB` opened with the list (`SeqDb::restrict`), then reads what its
/// alias files' counts and the reports need: NCBI's errors are `BLAST Database error: <text>`
/// with exit status 2 (`database_error`), and the warnings NCBI writes go to `diagnostics`
/// (`LOG_POST` at the warning level: the text alone; `ERR_POST`: `Warning: [<program>] <text>`,
/// the program in lower case as `SetDiagPostPrefix` names it; an `ERR_POST` of an error
/// before a throw: `Error: [<program>] <text>`), before an error too.
///
/// The one owner of the order of a search database's errors, which is NCBI's (probes p4c_e*,
/// p4c_f2 of line RG D1a-2): the alias files' tree and ranges (`open_database`), the list, the
/// alias files' counts (`SeqDb::check_alias_values`), the id of a database mask (`mask`, BLASTN's
/// `-db_soft_mask` / `-db_hard_mask`: `DbMaskOption::validate_at_open`), the formatter's
/// reopening of each name (`SeqDb::read_name_records`), then the first headers (`check_seqids`,
/// whose answer it returns: whether the database has parsed Seq-ids). A version 5 list that mixes
/// files with and without the accession table opens as in NCBI: the lookup of a parsed Seq-id
/// fails only where a report or the prefetch looks it up (`FoundInstead::Failed`).
pub(crate) fn restrict_database(
    db: &mut SeqDb,
    name: &str,
    list: Option<&crate::seqdb::DbList>,
    mask: Option<&crate::blastinput::db_mask::DbMaskOption>,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/uniform_search.cpp:290-298
    if let Some(mask) = mask {
        mask.validate_at_open(db, name, db.seqtype() == SeqType::Nucleotide, program)?;
    }
    db.read_name_records(name, mask.is_some())
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
    // The lines keep the bytes of a seqid list's words (probe R3 of line LIB-B).
    for warning in warnings {
        let line = match warning {
            crate::seqdb::ListWarning::Log(text) => text,
            crate::seqdb::ListWarning::Posted(text) => {
                format!("Warning: [{prefix}] {text}").into_bytes()
            }
            crate::seqdb::ListWarning::Error(text) => {
                [format!("Error: [{prefix}] ").as_bytes(), &text].concat()
            }
            crate::seqdb::ListWarning::SeqId(warning) => {
                let mut line = Vec::new();
                crate::writedb::fasta_source::defline_ids::post_seq_id_warning_of(
                    &prefix,
                    warning,
                    &mut |text| line.extend_from_slice(text),
                );
                line
            }
        };
        diagnostics.write_all(&line)?;
        diagnostics.write_all(b"\n")?;
    }
    restricted.map_err(|error| database_error(error, program))
}

/// Tells the BLAST archive of the run (`archive`, `None` without one) and its search strategy
/// export (`strategy`) what the database's list and statistics make of their requests: the run's
/// options handle, made once per run and shared by both
/// (`report::archive::OptionsHandle::set_database_filter`). `statistics` are the database's
/// (`SeqDb::statistics`),
/// `translated` whether its subjects are translated (TBLASTN, TBLASTX). The requests'
/// `EffectiveSearchSpace` is the search's own value where the statistics' length before the
/// fallback to OID 0 (`DbStatistics::gumbel_db_length`) is positive or `searchsp` (`-searchsp`
/// is given, or an RPS-BLAST search hands the space of its last search unit,
/// `rps_search::volumes::RpsUnits::request_space`), else
/// `BLAST_CalcEffLengths` with that length (`archive::RequestLength`: 0 for a list that kept no
/// sequence without `STATS_TOTLEN`; the context's query length for an alias file's `LENGTH` of
/// 2^63 or more). `dbsize` is whether `-dbsize` is given: the
/// `EffectiveSearchSpace` option then counts the length adjustment with no sequence, while the
/// searches read the length of OID 0 as one sequence (probes E14 and E17 of line DI-8): LOSAT
/// refuses it (the strategy's requests after the first carry that option too).
#[allow(clippy::too_many_arguments)]
pub(crate) fn note_archive_database_filter(
    archive: Option<&crate::report::archive::ArchiveHandle>,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    filter: &crate::blastinput::db_filter::DbFilterArgs,
    list_taxids: &[i32],
    statistics: &crate::seqdb::DbStatistics,
    translated: bool,
    searchsp: bool,
    dbsize: bool,
    program: &str,
) -> Result<()> {
    if archive.is_none() && strategy.is_none() {
        return Ok(());
    }
    let statistics_length = statistics.gumbel_db_length;
    if statistics_length <= 0 && dbsize {
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
    let request_length = if statistics_length > 0 || searchsp {
        crate::report::archive::RequestLength::Search
    } else {
        crate::report::archive::RequestLength::Lengths(crate::stats::setup::eff_lengths_database(
            statistics_length,
            statistics.request_count,
            0,
            false,
            translated,
        ))
    };
    options.set_database_filter(filter.archive_parameter(list_taxids), request_length)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:499-516 (cited at
// `blastinput::db_mask::ResolvedMask::archive_params`)
/// The database masking's program options in the run's options handle, which the archive
/// (`archive`) and the search strategy (`strategy`) share (`note_archive_database_filter`): the
/// archive's (`ResolvedMask::archive_params`, none after the formatter's warning) and the
/// strategy's (`DbMaskOption::strategy_params`), set where the formatter resolves the masking.
pub(crate) fn note_archive_database_masking(
    archive: Option<&crate::report::archive::ArchiveHandle>,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    archive_params: Vec<crate::report::archive::options::Param>,
    strategy_params: Vec<crate::report::archive::options::Param>,
) -> Result<()> {
    let options = match (archive, strategy) {
        (Some(handle), _) => handle.options(),
        (None, Some(strategy)) => strategy.options(),
        (None, None) => return Ok(()),
    };
    options.set_database_masking(archive_params, strategy_params)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153 (cited at
// `blastinput::db_mask::DbMaskOption::resolve`)
/// The formatter's step of the database masking `option` of a search of `db` (the `-db` value
/// `db_name`) by `program` (NCBI's upper-case name; series NP, M4b: BLASTP and BLASTX, and the
/// translated programs that take it): the resolution (`DbMaskOption::resolve`), whose warning goes
/// to `warning` (one diagnostic, which an archive saves as one message), or whose database error
/// ends the run before the formatter's other work; then the program options of the archive
/// `archive` and of the search strategy `strategy` (`note_archive_database_masking`). The result
/// is the search's masking (`None`: not masked).
pub(crate) fn resolve_database_masking(
    db: &SeqDb,
    db_name: &str,
    option: &crate::blastinput::db_mask::DbMaskOption,
    program: &str,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    warning: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<Option<crate::blastinput::db_mask::ResolvedMask>> {
    let resolved = option.resolve(db, db_name, program, warning)?;
    note_archive_database_masking(
        archive,
        strategy,
        resolved
            .as_ref()
            .map(crate::blastinput::db_mask::ResolvedMask::archive_params)
            .unwrap_or_default(),
        option.strategy_params(resolved.as_ref()),
    )?;
    Ok(resolved)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:153
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:123-129
/// The error of the first volume of `db` (in the order NCBI reads them) whose sequence file
/// cannot be opened, as NCBI's first read of a sequence meets it. A PSSM search (psiblast's and
/// tblastn's `-in_pssm`, `CPsiBlast`) runs its preliminary stage without `CLocalBlast`, whose catch
/// makes such an error `Database memory map file error` for a FASTA query (local_blast.cpp:229-252,
/// cited at `seqdb::SeqDbError`): the reader's own text reaches the application
/// (`database_error`: `BLAST Database error: Error: File (<path>) not found.`, exit status 2; NCBI
/// 2.17.0, probe `y.g.ps.tk.nopsq.0` of line PS PS5c, and review 4 of line RG F1b-2 for tblastn).
fn unopenable_sequence_file(db: &SeqDb) -> Option<SeqDbError> {
    db.volumes()
        .iter()
        .find_map(|volume| volume.open_sequence_file().err())
}

/// The error of the first read of the sequences of `db` in the preliminary search of a query batch
/// with a valid context, through `CLocalBlast` (a FASTA query of blastp, blastx, blastn,
/// tblastn and tblastx), after the batch's set-up checks: a sequence file that cannot be opened
/// (`unopenable_sequence_file`, the open of a volume's sequence file that the search's first read
/// of the volume makes, `seqdb::Volume::open_sequence_file`), as the search's read reports it
/// (`search_read_error`: `PreliminaryReadError`, NCBI's `Database memory map file error`; LOSAT's
/// refusal of an error that is not NCBI's). The set-up of a query batch that LOSAT refuses asks it
/// where the search would read (`query_batch::refused_batch_error`; line WQA D10 round 6,
/// 2026-11-20, review 4 M2: `-db` of a database without its `.psq` or `.nsq` file, NCBI BLAST+
/// 2.17.0 exits 3 with that error, and at outfmt 11 writes its archive first).
pub(crate) fn sequence_file_read_error(db: &SeqDb, program: &str) -> Option<anyhow::Error> {
    unopenable_sequence_file(db).map(|error| search_read_error(error, program))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:425
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:241-251 (the catch that a PSSM
// search's preliminary stage has not)
/// What the first read of the database of a search without `CLocalBlast` meets (psiblast, whose
/// `CPsiBlast` searches FASTA queries too, and PSI-TBLASTN), in NCBI's order: no OID to read
/// (`SeqDb::first_oid_read_fails`: the setup's `GetSeqLength` of OID 0, SeqDB's `OID not found`),
/// else a sequence file that cannot be opened (`unopenable_sequence_file`). SeqDB's own exception
/// reaches the application (`error`), where `CLocalBlast` gives a FASTA query of the other programs
/// `Database memory map file error` (`PreliminaryReadError`).
#[derive(Debug)]
pub(crate) enum PsiSearchReadFailure {
    EmptyDatabase,
    Unopenable(SeqDbError),
}

impl PsiSearchReadFailure {
    /// The failure of the first read of `db`, if any.
    pub(crate) fn of(db: &SeqDb) -> Option<Self> {
        if db.first_oid_read_fails() {
            Some(PsiSearchReadFailure::EmptyDatabase)
        } else {
            unopenable_sequence_file(db).map(PsiSearchReadFailure::Unopenable)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:930-961
    /// Whether NCBI meets this failure before the query batch's setup checks (`BLAST_GapAlignSetUp`:
    /// `BlastExtensionParametersNew`'s valid Karlin-Altschul block, which a PSSM query whose file's
    /// ungapped values are invalid fails, then the average subject length): the empty database in
    /// the setup's read of OID 0, before them; an unopenable sequence file in the search's first
    /// read, after them. The one order of psiblast (`psiblast_stop_at_invalid_query`) and PSI-TBLASTN
    /// (`tblastn_app`; review 4 R7-1 of line RG F1b-2).
    pub(crate) fn precedes_setup_checks(&self) -> bool {
        matches!(self, PsiSearchReadFailure::EmptyDatabase)
    }

    /// NCBI's error (`BLAST Database error: ...`, exit status 2; `program` names LOSAT's program in
    /// a refusal of an error that is not NCBI's, `database_error`).
    pub(crate) fn error(self, program: &str) -> anyhow::Error {
        match self {
            PsiSearchReadFailure::EmptyDatabase => NativeError {
                exit: 2,
                message: format!("{DATABASE_ERROR} OID not found\n").into(),
            }
            .into(),
            PsiSearchReadFailure::Unopenable(error) => database_error(error, program),
        }
    }
}

/// A reader error outside the preliminary search (when the database is opened, and when
/// the sequences and titles of the hit lists are read) as the command line reports it:
/// NCBI's own errors (`SeqDbError::ncbi_message`: a database that is not found, a file
/// that is no database file or cannot be opened) as `BLAST Database error: <text>` with
/// exit status 2; everything else is LOSAT's refusal, in LOSAT's wording (`program` is
/// the program that refuses).
pub(crate) fn database_error(error: SeqDbError, program: &str) -> anyhow::Error {
    database_error_of(&error, program)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
// `CATCH_ALL` writes an exception that is not a `CSeqDBException` with `what()`, which names the
// place of each throw: one thrown by the code NCBI's `datatool` generates (`ThrowSite::generated`,
// the `CObject_id::GetId` of a `BL_ORD_ID` or `PIG` id with a string tag under an id list, line
// RG D1c-1) names a place LOSAT does not have, and LOSAT refuses it with its text, as segmasker
// does (`segmasker::work::database_failure`).
/// The exception of NCBI's that a search program reports as NCBI's does
/// (`SeqDbError::ncbi_exception` whose places LOSAT knows), the one test of [`database_error_of`]
/// and `header_read_error`.
fn reported_exception(error: &SeqDbError) -> Option<&crate::ncbi_exception::NcbiException> {
    error
        .ncbi_exception()
        .filter(|exception| exception.places_known())
}

/// [`database_error`] of an error that its owner keeps (`seqdb::NameMaskError::Table`).
pub(crate) fn database_error_of(error: &SeqDbError, program: &str) -> anyhow::Error {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    // A header record that ends before NCBI's read of it enters the set's elements throws the
    // `CEofException` that `HandleEOF` rethrows (`NCBI_RETHROW_SAME`, objistr.cpp:604): a
    // `CIOException` whose `GetErrCode()` is `eInvalid` (a derived class, ncbiexpt.hpp), which
    // this branch takes before the `CException` one and posts nothing: the run ends there with exit
    // status 0 (NCBI 2.17.0, probe RG-T5d fix 1 le3, le7, le8: after the outfmt 0 prolog and the
    // reports of the batches before).
    if error
        .ncbi_exception()
        .is_some_and(|exception| exception.is_class("CEofException"))
    {
        return NativeError::silent_end().into();
    }
    if let Some(exception) = reported_exception(error) {
        return NativeError::exception(exception).into();
    }
    if let Some(exception) = error.ncbi_exception() {
        return refusal(
            &format!(
                "the report of the exception '{}', which names the place in NCBI's source that \
                 throws it",
                String::from_utf8_lossy(exception.get_msg())
            ),
            program,
        );
    }
    match error.ncbi_message() {
        Some(text) => NativeError {
            exit: 2,
            message: [DATABASE_ERROR.as_bytes(), b" ", &text, b"\n"].concat(),
        }
        .into(),
        None => anyhow::Error::msg(losat_reader_error(error, program)),
    }
}

/// LOSAT's refusal of a database that `program` does not search as NCBI does.
pub(crate) fn refusal(reason: &str, program: &str) -> anyhow::Error {
    anyhow::anyhow!("{reason}; this database is not supported by LOSAT's {program}")
}

/// The text of a reader error that is not NCBI's (`SeqDbError::ncbi_message` `None`), as LOSAT
/// reports it: a file whose bytes the system does not give is LOSAT's error as it is
/// (`SeqDbError::is_unreadable`: "cannot be read", where NCBI ends by SIGBUS; review M1 of line
/// WQA D18, 2026-11-21), any other LOSAT's refusal of the database in `program`'s name
/// ([`refusal`]). The one place of that choice for the search programs, blastdbcmd and the
/// maskers' database reader; each ends with LOSAT's exit status of its own errors.
pub(crate) fn losat_reader_error(error: &SeqDbError, program: &str) -> String {
    if error.is_unreadable() {
        error.to_string()
    } else {
        refusal(&error.to_string(), program).to_string()
    }
}

/// What LOSAT refuses where the statistics read 0 sequences of a positive total length
/// (`check_setup_count`).
pub(crate) const STATISTICS_COUNT_ZERO: &str = "the statistics count 0 sequences of the database (an alias file's NSEQ 0, or its LENGTH or STATS_TOTLEN over volumes without sequences), by which NCBI BLAST+ divides its length";

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:322-327
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:853-877
/// Whether the setup of a search divides its total length `setup_total` by a count `setup_count`
/// of 0 (`SeqDb::statistics_numbers`, `DbStatistics::setup_total` and `setup_count`: an alias
/// file's `NSEQ 0` without `STATS_NSEQ`, its `LENGTH` or `STATS_TOTLEN` over volumes without
/// sequences, a list that kept no sequence under a `STATS_TOTLEN`; an RPS unit's numbers): without
/// the Gumbel block (`gumbel` false: BLASTN, TBLASTX, and the protein programs with `-searchsp`,
/// `OLD_FSC`, BLASTX also with `-ungapped`) `BLAST_GapAlignSetUp`'s shortest subject length is the
/// total over the count, and NCBI BLAST+ 2.17.0 crashes there (SIGFPE, DECISIONS 14 of line RG; an
/// alias file's `LENGTH` 100 alone over a volume of no OID, probe p4c_dn_nB of D1a-2). NCBI makes
/// that division only in the setup of a search with a valid context (`CLocalBlast::Run` stops a
/// query batch whose `CheckInternalData`, `BlastScoreBlkCheck`, finds none; PSI-BLAST's rounds:
/// `blastp::blast_engine::call_setup::setup_fails`), before it reads a sequence: a total that is
/// not positive counts 1, so the count is 0 only of a positive total. The one test of every program
/// (`check_setup_count`).
pub(crate) fn setup_count_refused(setup_total: i64, setup_count: i32, gumbel: bool) -> bool {
    !gumbel && setup_total > 0 && setup_count == 0
}

/// LOSAT's refusal of `setup_count_refused` for `program`, which each program makes where NCBI's
/// setup meets the division: in the setup of the first query batch (or query chunk, RPS unit or
/// PSI-BLAST round) with a valid context. BLASTN, TBLASTX, TBLASTN and RPSTBLASTN write the
/// reports of the batches before it (NCBI's output before its crash in the tabular formats; in the
/// others the crash cuts NCBI's buffered output short), BLASTP, PSI-BLAST and RPSBLAST refuse
/// before any output. A run whose queries are all invalid is reported as NCBI reports it (exit 0;
/// D1a-5 of line RG, DECISIONS D1a4-12).
pub(crate) fn check_setup_count(
    setup_total: i64,
    setup_count: i32,
    gumbel: bool,
    program: &str,
) -> Result<()> {
    if setup_count_refused(setup_total, setup_count, gumbel) {
        return Err(refusal(STATISTICS_COUNT_ZERO, program));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
/// A reader error while the preliminary search reads the sequences: NCBI's search stage
/// replaces the reader's text (`SeqDbError::ncbi_search_message`), and the application
/// reports `BLAST engine error: <text>` with exit status 3. The error is kept apart from
/// the others (`PreliminaryReadError`), because NCBI has written the outfmt 0 prolog and
/// read the query batch before it.
pub(crate) fn search_read_error(error: SeqDbError, program: &str) -> anyhow::Error {
    search_read_error_of(&error, program)
}

/// [`search_read_error`] of an error that its owner keeps (`seqdb::MaskBlock`: the masks that the
/// preliminary search reads with each sequence, seqsrc_seqdb.cpp:318-321).
pub(crate) fn search_read_error_of(error: &SeqDbError, program: &str) -> anyhow::Error {
    match error.ncbi_search_message() {
        Some(text) => PreliminaryReadError(text).into(),
        None => anyhow::Error::msg(losat_reader_error(error, program)),
    }
}

/// Bytes of the sequence file that one read of `read_nucleotide_stores` holds at most
/// (`SeqDb::next_block`; a longer sequence is a read of its own). No NCBI counterpart (NCBI
/// maps the volume files); the size does not change what is read.
const STORE_BLOCK_BYTES: u64 = 1 << 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1405
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:771-775
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
            let ordinal = range.first_ordinal + i as u32;
            let (length, runs) = stored_parts(&block, i, || db.oid(ordinal), program)?;
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
                        db.oid(range.first_ordinal)
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
            let ordinal = range.first_ordinal + i as u32;
            let (length, _) = stored_parts(&block, i, || db.oid(ordinal), program)?;
            checked.add(length as u32);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(checked)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:312-322
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
    let (length, runs) = stored_parts(&block, 0, || db.oid(ordinal), program)?;
    let packed = block.x_into_nucleotide_packed();
    if crate::utils::x_mapped::x_is_mapped(packed.capacity()) {
        crate::utils::x_mapped::X_BYTES.x_disown(&packed);
    }
    Ok(PackedNucleotide::x_from_stored_buffer(
        packed, length, &runs,
    ))
}

/// The length and the ambiguity runs (as stored) of sequence `i` of `block` (`oid` gives its
/// OID, which a refusal names), with the errors of the sequence's store
/// (`x_read_nucleotide_store`), without the store.
pub(crate) fn x_nucleotide_parts_in(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<(usize, Vec<crate::seqdb::AmbiguityRun>)> {
    stored_parts(block, i, oid, program)
}

/// The length and the ambiguity runs of sequence `i` of `block` (`oid` gives its OID, which a
/// refusal names), with the errors of `read_nucleotide_stores` (which reads every sequence
/// through it). NCBI's preliminary
/// search reads the 2-bit data only; it reads the ambiguity data in the traceback, of the
/// reported sequences. LOSAT reads them for every sequence, so data it cannot read are LOSAT's
/// refusal, not an NCBI error.
fn stored_parts(
    block: &crate::seqdb::SequenceBlock,
    i: usize,
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<(usize, Vec<crate::seqdb::AmbiguityRun>)> {
    let length = block
        .nucleotide_len(i)
        .map_err(|error| search_read_error(error, program))?;
    let runs = block.ambiguity_runs(i).map_err(|error| {
        refusal(
            &format!(
                "the ambiguity data of database sequence {} cannot be read ({error})",
                oid()
            ),
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
        write!(f, "{ENGINE_ERROR} {}", self.0)
    }
}

impl std::error::Error for PreliminaryReadError {}

/// A reader error of NCBI's own while the definition lines of the final hit lists are
/// read (`stored_title`), or while rpsblast's traceback of a batch reads a profile volume's
/// sequence file (`rps_hooks::failed_batch_error`), or where a batch's prefetch looks up a
/// Seq-align id whose lookup fails on a file (`batch_prefetch`): `query` is the query whose hit
/// list reached it (rpsblast and the prefetch: the batch's first query). NCBI fails that query's
/// batch, after the reports of the batches before. BLASTP's `-use_sw_tback` redo whose batch's
/// last redone match fails (`Traceback failed`, `blastp::blast_engine::redo_run`) is such an error
/// of NCBI's search too, and BLASTP's refusals of a step where NCBI's search of a batch stops,
/// crashes or never ends (a Smith-Waterman redo without end, a failed batch whose first query
/// holds a redone match, a query chunk split again or without a query) take its form (the batch's
/// first query): the type of an error of a batch's search after the reports of the batches before
/// in the drivers that search all batches together, which name the failing batch. The batch loops
/// of TBLASTX and TBLASTN, which search a batch at a time, hold the same error without the query
/// (`blastinput::app::ErrorAfterBatchReports`).
#[derive(Debug)]
pub(crate) struct HeaderReadError {
    pub(crate) query: usize,
    pub(crate) error: anyhow::Error,
    /// The first UTF8String that the reads before it (in NCBI's order), or the failed read
    /// itself before its failure, took for a VisibleString (`HeaderReader::first_utf8`): NCBI
    /// posted its warning before the error (`Utf8Read::post`).
    pub(crate) utf8: Option<Utf8Read>,
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:304-309
/// NCBI's first read, in a run of reported sequences' headers (`HeaderReader`), of a record that
/// holds a UTF8String where a VisibleString is expected: the query whose hit list reached it
/// (numbered as `HeaderReadError::query`) and the member NCBI's warning names
/// (`seqid::ber::utf8_string_path`). NCBI reads the headers of a batch's final hit lists when
/// its traceback's results become Seq-aligns (`read_titles` of each program), so it posts the
/// warning there, once per process: after the batch's reader messages and the posts of its
/// search strategy request, before anything the batch writes after its search (`post`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Utf8Read {
    pub(crate) query: usize,
    pub(crate) path: String,
}

/// The places of a batch's reads of header records in NCBI's order (`HeaderReader::first_utf8`):
/// the fetches with `check_oid_exclusion` of its traceback or redo (`FetchExcluded`), then its
/// Seq-aligns of the final hit lists (`SeqAligned` under a seqid list, which reads them all),
/// the reads of the reported sequences (`read_header`, which the Seq-aligns of the other
/// databases read). The reads of NCBI's scope that come after them (the data loader's, where it
/// resolves a Seq-align id: a lookup's version check, the sequence the lookup gives first) are
/// the loads' (`ParsedSubject::scope_utf8`, posted where the model of the loads makes the load,
/// `load_read_texts`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Utf8Phase {
    Fetch,
    Aligned,
    Reported,
}

impl Utf8Phase {
    /// The place in `Utf8Reads`' order of a read of this phase for the query `query` of the
    /// batch that starts at `batch` (of equal places, the first noted is the first): a batch's
    /// reads come before the next batch's.
    fn order(self, batch: usize, query: usize) -> u128 {
        (batch as u128) << 96 | (self as u128) << 88 | ((query as u128) & ((1 << 56) - 1)) << 32
    }

    /// The query of the place `order`.
    fn query(order: u128) -> usize {
        ((order >> 32) & ((1 << 56) - 1)) as usize
    }
}

/// The reads of a search call's batches that come before those of their reported sequences in
/// NCBI's order and that took a UTF8String for a VisibleString (`HeaderReader::note_earlier`),
/// each with a query of its batch: the first of its fetches with `check_oid_exclusion`
/// (`FetchExcluded`: those of the traceback, `FetchedUtf8::traceback_first`, or of the redo,
/// `redo_fetches`) and of its Seq-aligns under a seqid list (`SeqAligned`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct EarlierUtf8 {
    pub(crate) fetched: Option<(usize, crate::seqdb::SeqLocation)>,
    pub(crate) aligned: Option<(usize, crate::seqdb::SeqLocation)>,
}

impl Utf8Read {
    /// The read `note` of a `HeaderReader` (its query from the place, `Utf8Phase::order`), with
    /// the member of NCBI's warning; `None` where the record holds no such string before a
    /// failure.
    #[cold]
    #[inline(never)]
    fn of_note(db: &SeqDb, note: crate::seqdb::Utf8Note) -> Option<Self> {
        Some(Utf8Read {
            query: Utf8Phase::query(note.order),
            path: db.utf8_string_path_at(note.at)?,
        })
    }

    /// NCBI's warning of this read, with the post prefix `program` (the lower-case program
    /// name, `[blastp]`), where the run of `db` has not posted it (`SeqDb::take_utf8_post`):
    /// the posts that a program writes there (`Posts::write` makes it a message of a BLAST
    /// archive), `None` once posted.
    pub(crate) fn post(&self, db: &SeqDb, program: &str) -> Option<crate::report::diag_log::Posts> {
        db.take_utf8_post().then(|| {
            crate::report::diag_log::Posts::of_lines(crate::seqid::ber::utf8_string_warning(
                Some(program),
                &self.path,
            ))
        })
    }

    /// Before the model of the loads of the batch of the queries `batch` (a program that models
    /// every batch before it writes a report: BLASTP's and TBLASTN's call): where the read `read`
    /// (the call's first, `HeaderReader::first_utf8`) is of that batch or one before, its batch
    /// posts NCBI's warning before its loads, so no load of the model posts it from here on
    /// (`SeqDb::end_scope_utf8_loads`).
    pub(crate) fn reach(read: Option<&Utf8Read>, batch: &std::ops::Range<usize>, db: &SeqDb) {
        if read.is_some_and(|read| read.query < batch.end) {
            db.end_scope_utf8_loads();
        }
    }

    /// `post` of the read `read`, where there is one.
    pub(crate) fn post_of(
        read: Option<&Utf8Read>,
        db: &SeqDb,
        program: &str,
    ) -> Option<crate::report::diag_log::Posts> {
        read?.post(db, program)
    }

    /// The first read of a batch's read of its reported sequences' titles where NCBI makes its
    /// Seq-aligns (a program's `read_titles`, `read`) that met a UTF8String: `utf8` of what it
    /// read, or the one that a failed read made before it failed (`HeaderReadError::utf8`). NCBI's
    /// warning of it comes after the batch's reader messages and request posts, before its
    /// prefetch and before its error.
    pub(crate) fn of_read<'r, T>(
        read: &'r Result<T>,
        utf8: impl FnOnce(&'r T) -> Option<&'r Utf8Read>,
    ) -> Option<&'r Utf8Read> {
        match read {
            Ok(read) => utf8(read),
            Err(error) => error
                .downcast_ref::<HeaderReadError>()
                .and_then(|failure| failure.utf8.as_ref()),
        }
    }

    /// `post` of the read `read`, where its query is one of `queries` (the batch that reaches
    /// it).
    pub(crate) fn post_in(
        read: Option<&Utf8Read>,
        queries: std::ops::Range<usize>,
        db: &SeqDb,
        program: &str,
    ) -> Option<crate::report::diag_log::Posts> {
        read.filter(|read| queries.contains(&read.query))
            .and_then(|read| read.post(db, program))
    }
}

impl EarlierUtf8 {
    /// The earlier reads of the batch whose first query is `query`: the note of its first fetch
    /// with `check_oid_exclusion` that met a UTF8String (`fetched`) and that of its Seq-aligns
    /// under a seqid list (`aligned`, a program's `drop_unaligned`).
    pub(crate) fn of_batch(
        query: usize,
        fetched: Option<crate::seqdb::SeqLocation>,
        aligned: Option<crate::seqdb::SeqLocation>,
    ) -> Self {
        Self {
            fetched: fetched.map(|at| (query, at)),
            aligned: aligned.map(|at| (query, at)),
        }
    }
}

/// The first read of a search call's batches, among their own reads (their fetches, Seq-aligns
/// and the reads of their reported sequences, `HeaderReader::first_utf8`), of a header with a
/// UTF8String where a VisibleString is expected, with the first query of its batch, for a program
/// that models every batch before it writes a report (BLASTP's and TBLASTN's call): NCBI's
/// warning comes in that query, after the batch's reader messages and request posts (`post`,
/// `query_warnings::query_posts`), or before the error of the batch where it fails after the read
/// (`failed_post`), where no load of the model of the loads of a batch before posted it
/// (`load_read_texts`; no load of the read's batch or a later one does, `reach`).
#[derive(Debug, Default)]
pub(crate) struct CallUtf8(std::sync::OnceLock<(usize, Utf8Read)>);

impl CallUtf8 {
    /// Notes `read` of the batch whose first query is `batch`, where no read of the call came
    /// before.
    pub(crate) fn note(&self, batch: usize, read: Option<&Utf8Read>) {
        if let Some(read) = read {
            let _ = self.0.set((batch, read.clone()));
        }
    }

    /// The read, where there is one.
    pub(crate) fn read(&self) -> Option<&Utf8Read> {
        self.0.get().map(|(_, read)| read)
    }

    /// [`Utf8Read::reach`] of the read, before the model of the loads of the batch of the
    /// queries `batch`.
    pub(crate) fn reach(&self, batch: &std::ops::Range<usize>, db: &SeqDb) {
        Utf8Read::reach(self.read(), batch, db);
    }

    /// NCBI's warning of the read in the first query of its batch, where that batch is reported
    /// (its first query below `reported`) and the run has not posted it (`Utf8Read::post`), with
    /// the post prefix `program`.
    pub(crate) fn post(
        &self,
        db: Option<&SeqDb>,
        program: &str,
        reported: usize,
    ) -> Option<(usize, crate::report::diag_log::Posts)> {
        let (query, read) = self.0.get().filter(|(query, _)| *query < reported)?;
        Some((*query, read.post(db?, program)?))
    }

    /// NCBI's warning of the read where the failed batch whose first query is `batch` made it
    /// before its failure: written after the batch's reader messages and request posts.
    pub(crate) fn failed_post(
        &self,
        db: Option<&SeqDb>,
        program: &str,
        batch: usize,
    ) -> Option<crate::report::diag_log::Posts> {
        let (_, read) = self.0.get().filter(|(query, _)| *query == batch)?;
        read.post(db?, program)
    }
}

impl HeaderReadError {
    /// `error` of a batch that fails after NCBI's reads of its Seq-aligns' headers (its
    /// prefetch), with `read`, the first of those reads that met a UTF8String (`first_utf8`),
    /// where it is a `HeaderReadError` without one: NCBI posted its warning before the failure.
    pub(crate) fn after_reads(error: anyhow::Error, read: Option<&Utf8Read>) -> anyhow::Error {
        match error.downcast::<HeaderReadError>() {
            Ok(mut failure) => {
                if failure.utf8.is_none() {
                    failure.utf8 = read.cloned();
                }
                anyhow::Error::new(failure)
            }
            Err(error) => error,
        }
    }
}

impl std::fmt::Display for HeaderReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for HeaderReadError {}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
/// What the `Database:` lines of the pairwise report show of the database `db`: for each
/// name of its `-db` list (`SeqDb::name_records`), that database's title (its name for an
/// empty title), date and real totals, also with `-dbsize` or a list.
pub(crate) fn pairwise_database(db: &SeqDb) -> PairwiseDatabase {
    pairwise_database_masked(db, None)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:393-406
/// [`pairwise_database`] of a search masked with an algorithm whose details in each name's
/// database are `shown` (`SeqDb::name_mask_details`, one per name record; `None` without database
/// masking: `dbfilt_algorithm` is -1): the `Masked using:` lines of each name.
pub(crate) fn pairwise_database_masked(
    db: &SeqDb,
    shown: Option<&[crate::seqdb::MaskAlgorithm]>,
) -> PairwiseDatabase {
    let records = db.name_records();
    if let Some(shown) = shown {
        // One algorithm per name (`name_mask_details` reads them from these records): a list of
        // another length would drop or shift the names' `Masked using:` lines.
        assert_eq!(
            shown.len(),
            records.len(),
            "one masking algorithm per name of the -db list"
        );
    }
    PairwiseDatabase::new(
        records
            .iter()
            .enumerate()
            .map(|(index, record)| {
                let algorithm = shown.map(|shown| &shown[index]);
                DatabaseInfo {
                    definition: if record.title.is_empty() {
                        record.name.as_bytes().to_vec()
                    } else {
                        record.title.clone()
                    },
                    date: record.date.clone(),
                    num_sequences: record.num_seqs,
                    // `SDbInfo::total_length` is an `Int8`: a `LENGTH` of 2^63 or more is
                    // negative.
                    total_letters: record.total_length as i64,
                    filt_algorithm_name: algorithm
                        .map(|algorithm| algorithm.program_name.clone())
                        .unwrap_or_default(),
                    filt_algorithm_options: algorithm
                        .map(|algorithm| algorithm.options.clone())
                        .unwrap_or_default(),
                }
            })
            .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:513-522
// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:116-119 (`m_NumSequences`
// is an `int`, `m_NumBases` an `Int8`) and c++/src/algo/blast/format/data4xml2format.cpp:166-169
// (`m_NumSequences` and `m_NumBases` are both `Int8`)
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1002-1010
// The sums are signed overflow where they pass their type (NCBI 2.17.0 prints them wrapped,
// probes p5_pr_sum2 and p5b_sum_int of line RG D1a-3): LOSAT refuses these (DECISIONS D1a3-2 of
// line RG: Q2 = A, rule 3). The serial XML writer's `-n` of INT64_MIN is signed overflow too (it
// writes `-0`, probes p5_neg_2p63 and p5b_neg_fmt): LOSAT writes the digits of a summed length of
// exactly -2^63 (`big_int`; PD-LOSAT-NCBI-DEFECTS approved exception 17 (Int8 -2^63 written as
// digits), the Owner's Q1 = A of 2026-11-17, DECISIONS OQ-1117 of line RG).
/// What LOSAT refuses of the database counts that the report of `outfmt` (NCBI's format number)
/// would print where NCBI's arithmetic overflows (`num_sequences` is the `PairwiseDatabase` sum
/// over the names of the `-db` list, `total_letters_overflow` whether its lengths' sum passes
/// `Int8`): a sum of the total lengths past `Int8` (the pairwise prolog of
/// outfmt 0 to 4, BLAST XML 5, the XML2 and JSON documents 13 to 16) and a sum of the numbers of
/// sequences past NCBI's `int` (BLAST XML).
pub(crate) fn report_sum_refusal(
    num_sequences: i64,
    total_letters_overflow: bool,
    outfmt: i32,
) -> Option<&'static str> {
    if matches!(outfmt, 0..=5 | 13..=16) && total_letters_overflow {
        return Some(
            "the total lengths of the databases of the -db list (their alias files' LENGTH) sum to more than NCBI BLAST+'s 64-bit report count holds",
        );
    }
    if outfmt == 5 && i32::try_from(num_sequences).is_err() {
        return Some(
            "the numbers of sequences of the databases of the -db list (their alias files' NSEQ) sum to more than BLAST XML's 32-bit count holds",
        );
    }
    None
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-747
// NCBI sums the counts of the pairwise prolog (outfmt 0 to 4) before it searches, those of the
// XML and JSON documents when it formats a query's report: a search that fails before (an engine
// error of the setup, "The average subject length is too short", probe p5_pr_sum2#6 of line RG
// D1a-3) never reaches the documents' sums. The documents' refusals are LOSAT's report writers'
// (`StructuredXml::start` of `report::blast_xml::writer` and `report::blast_xml2::writer`), made
// before the first byte of the document (the XML2 and JSON documents' start included, which NCBI
// writes before the search: BLASTX holds it, `blastx::blast_xml::HeldOutput`).
/// LOSAT's refusal of the database counts of the report of `outfmt` that NCBI's sums overflow
/// (`report_sum_refusal`): its own type, so that a run that holds its document drops it
/// (`blastx::blast_xml::drops_held_document`) and the refusal follows no byte.
#[derive(Debug)]
pub(crate) struct ReportSumRefusal(String);

impl std::fmt::Display for ReportSumRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ReportSumRefusal {}

/// LOSAT's refusal (`ReportSumRefusal`) of the database counts that the report of `outfmt` would
/// print where NCBI's sums overflow (`report_sum_refusal`), for `program` (its name in the
/// message): `num_sequences` and `total_letters_overflow` are the `PairwiseDatabase` sums of the
/// `-db` list.
pub(crate) fn check_sums_of(
    num_sequences: i64,
    total_letters_overflow: bool,
    outfmt: i32,
    program: &str,
) -> Result<()> {
    match report_sum_refusal(num_sequences, total_letters_overflow, outfmt) {
        Some(reason) => Err(ReportSumRefusal(refusal(reason, program).to_string()).into()),
        None => Ok(()),
    }
}

/// LOSAT's refusals of the counts of the pairwise prolog of `outfmt` (`report_sum_refusal` for
/// outfmt 0 to 4, which NCBI writes before the search) for the `-db` list's databases `database`.
pub(crate) fn check_report_sums(
    database: &PairwiseDatabase,
    outfmt: i32,
    program: &str,
) -> Result<()> {
    if !prints_report_sums(outfmt) {
        return Ok(());
    }
    check_sums_of(
        database.num_sequences,
        database.total_letters_overflow,
        outfmt,
        program,
    )
}

/// Whether the report format `outfmt` prints the database's sums (`check_report_sums`).
fn prints_report_sums(outfmt: i32) -> bool {
    (0..=4).contains(&outfmt)
}

/// [`check_report_sums`] for each report format of `formats` (its `-outfmt` value and its
/// refusal), with the database's report made by `shown` once, only where a format prints its
/// sums: each format's refusal goes to that format's `FormatRefusal::check`, which fails the call
/// (the command line's format) or that format alone, as a program's other checks of a format
/// before its rows.
pub(crate) fn check_formats_report_sums<'f, 'r>(
    shown: impl FnOnce() -> PairwiseDatabase,
    formats: impl IntoIterator<Item = (&'f str, &'r mut crate::api::local_blast::FormatRefusal)>,
    program: &str,
) -> Result<()> {
    let mut shown = Some(shown);
    let mut made = None;
    for (outfmt, refusal) in formats {
        let number = crate::blastinput::app::parse_formatting_string(outfmt)?.number;
        if !prints_report_sums(number) || refusal.refused() {
            continue;
        }
        let database = made.get_or_insert_with(|| (shown.take().expect("made once"))());
        refusal.check(check_report_sums(database, number, program))?;
    }
    Ok(())
}

/// The stored title of a database sequence, as the reports read it (`stored_title`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct StoredTitle {
    /// The title of the sequence's definition line, with a byte outside ' ' ..= '~' read
    /// as '#'.
    pub(crate) title: Vec<u8>,
    /// Whether the stored title held a byte outside ' ' ..= '~' (`check_reported_title`).
    pub(crate) replaced_bytes: bool,
    /// The name of the organism of the source of the sequence's Bioseq, from which NCBI makes an
    /// empty title (`HeaderReader::stored_source`; `report::defline::BioseqFacts::of_source`).
    pub(crate) source: crate::report::defline::SourceTaxname,
    /// The `taxid` of the definition line (0 for none), read with the title: the Bioseq's source
    /// where it is positive (`HeaderReader::stored_source`).
    pub(crate) taxid: i32,
}

impl StoredTitle {
    /// LOSAT's refusal of the sequence (`oid` gives its OID, which the refusal names; `program`
    /// names the program) whose title NCBI makes from a source that LOSAT refuses
    /// (`SourceTaxname::Refused`), for a report that shows that title where it cannot refuse it
    /// itself: outfmt 0 (the tabular subject ids refuse it in `tabular_subject_id`, outfmt 18 in
    /// `report::tax_report`).
    pub(crate) fn check_shown_source(
        &self,
        oid: impl Fn() -> crate::seqdb::Oid,
        program: &str,
    ) -> Result<()> {
        source_name(&self.source, oid, program).map(drop)
    }
}

/// The name of the source `source` of a database sequence without parsed Seq-ids (`oid` gives its
/// OID, which a refusal names), from which NCBI makes its empty title; `Err` with LOSAT's refusal
/// of a source that it refuses (`SourceTaxname::Refused`; `program` names the program).
pub(crate) fn source_name(
    source: &crate::report::defline::SourceTaxname,
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<Option<std::sync::Arc<[u8]>>> {
    source.name().map_err(|reason| {
        anyhow::anyhow!(
            "database sequence {} has {reason}; this is not supported by LOSAT's {program}",
            oid()
        )
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310
// NCBI reference for the titles:
// c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:970-1017 (`s_GetBioseqTitle`: the
// title of the first definition line, the others appended as ` >id title`).
/// Reads the definition line of a database sequence (`oid` gives its OID, which a refusal names)
/// for a report, where the hit list of query `query` reaches it; `deflines` are its definition
/// lines as read
/// (`SeqDb::deflines`, or `SeqDb::visit_deflines` with the same results). A reader error of
/// NCBI's own is a `HeaderReadError` that names the query (NCBI fails that query's batch);
/// LOSAT's refusals (`program` names the program) are errors of their own.
pub(crate) fn stored_title(
    deflines: std::result::Result<Vec<crate::seqdb::OrdinalDefline>, crate::seqdb::SeqDbError>,
    oid: impl Fn() -> crate::seqdb::Oid,
    query: usize,
    program: &str,
) -> Result<StoredTitle> {
    let mut deflines = deflines.map_err(|error| header_read_error(error, query, program, None))?;
    // makeblastdb without -parse_seqids writes one definition line per sequence; NCBI
    // shows the others of a sequence as ` >id title` after the first.
    if deflines.len() != 1 {
        bail!(
            "database sequence {} has {} definition lines; a sequence with more or fewer than one is not supported by LOSAT's {program}",
            oid(),
            deflines.len()
        );
    }
    let defline = deflines.remove(0);
    Ok(StoredTitle {
        title: defline.title,
        replaced_bytes: defline.replaced_title_bytes > 0,
        source: Default::default(),
        taxid: defline.taxid,
    })
}

/// Keeps the source `source` of the reported sequence of rank `rank` (of `ranks`) in `sources`,
/// the side data of the reports by rank (`StoredTitle::source`), which stays empty, and is not
/// allocated, until a sequence has a source that names its organism or that LOSAT refuses.
/// Inlined: the reports of every search ask it for each reported sequence, almost all unnamed.
#[inline]
pub(crate) fn keep_source(
    sources: &mut Vec<crate::report::defline::SourceTaxname>,
    rank: usize,
    ranks: usize,
    source: crate::report::defline::SourceTaxname,
) {
    if matches!(source, crate::report::defline::SourceTaxname::Unnamed) {
        return;
    }
    if sources.is_empty() {
        sources.resize(ranks, Default::default());
    }
    sources[rank] = source;
}

/// A reader error of the header of a reported sequence where the hit list of query `query`
/// reaches it: NCBI's own (a database error, `SeqDbError::ncbi_message`, or another exception
/// of NCBI's that LOSAT reports as NCBI's, `reported_exception`: the `CSerialException` of its
/// read of the record, or the conversion of a `ti` tag under an id list, line RG D1c-1) is a
/// `HeaderReadError` that names the query, so the batch fails after the outfmt 0 prolog as NCBI's
/// does; LOSAT's refusals are errors of their own.
fn header_read_error(
    error: SeqDbError,
    query: usize,
    program: &str,
    utf8: Option<&Utf8Read>,
) -> anyhow::Error {
    if error.ncbi_message().is_some() || reported_exception(&error).is_some() {
        anyhow::Error::new(HeaderReadError {
            query,
            error: database_error(error, program),
            utf8: utf8.cloned(),
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

/// Whether the sequences of a database search keep a definition line under its seqid list
/// (`SeqDb::line_read`), each read once: the one check of NCBI's two places that drop a
/// sequence without a line, the fetches with `check_oid_exclusion` (`FetchExcluded`) and the
/// Seq-aligns (`SeqAligned`), with whether NCBI's read of the record there takes a UTF8String
/// for a VisibleString (its warning, `seqdb::Utf8Reads`). A header that cannot be read as
/// NCBI's header read adjusts it is not kept in the cache: each place decides what its error
/// does.
struct LineCheck<'a> {
    db: &'a SeqDb,
    kept: std::sync::Mutex<std::collections::HashMap<u32, LineRead>>,
}

/// What a `LineCheck` read of a header found: whether the sequence keeps a line, and whether
/// the record holds a UTF8String where a VisibleString is expected.
#[derive(Clone, Copy, Debug)]
struct LineRead {
    kept: bool,
    utf8: bool,
}

impl<'a> LineCheck<'a> {
    fn new(db: &'a SeqDb) -> Self {
        LineCheck {
            db,
            kept: std::sync::Mutex::new(std::collections::HashMap::new()),
        }
    }

    /// What the read of `ordinal`'s header finds (read at its first check; the lock is not held
    /// while the header is read and decoded), or the reader's error.
    fn read(&self, ordinal: u32) -> Result<LineRead, SeqDbError> {
        let cached = self
            .kept
            .lock()
            .expect("the line cache is not poisoned")
            .get(&ordinal)
            .copied();
        if let Some(read) = cached {
            return Ok(read);
        }
        let (kept, utf8) = self.db.line_read(ordinal)?;
        let read = LineRead { kept, utf8 };
        self.kept
            .lock()
            .expect("the line cache is not poisoned")
            .insert(ordinal, read);
        Ok(read)
    }
}

/// A subject of the items of `FetchExcluded::retain` whose header holds a UTF8String where a
/// VisibleString is expected: NCBI's fetch of it reads the header (`GetSeqIDs`) and posts its
/// warning, where the fetch is made (the ordinary traceback fetches every subject,
/// `FetchedUtf8::traceback_first`; the composition-based redo only those its early termination
/// test leaves, `redo_fetches`). Its item's index before the retain, whether the fetch keeps it,
/// the number of items kept before it (its index after the retain where it is kept) and its
/// ordinal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct FetchedUtf8 {
    pub(crate) index: usize,
    pub(crate) kept: bool,
    pub(crate) before: usize,
    pub(crate) ordinal: u32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1580
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:354-359
// A match inserts into its own query's heap only. A heap is filled to its cutoff once the
// records within the cutoff that the redo inserted reach `heapThreshold`: an insert pops a record
// only where the heap holds one beyond the cutoff (its worst), so every record within the cutoff
// that the redo inserted stays to the end, and the heap at a match's turn is filled where enough
// of them came before the match (the redo's order: query by query for a protein query,
// blast_hspstream.c:289-319; by subject for blastx, whose stream does not sort by score,
// blast_hspstream.c:676-686).
/// One heap of a batch's composition-based redo at a match's turn (`redo_fetches`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct RedoHeap {
    pub(crate) ecutoff: f64,
    pub(crate) heap_threshold: usize,
    /// The records within `ecutoff` that the redo had inserted before the match: those of the
    /// heap as the redo left it whose matches come before the match's.
    pub(crate) inserted_before: usize,
    /// No NCBI counterpart: `false` for a heap of a redo that ends no match early (BLASTP's
    /// unified P-values, `blastp::hsp::BlastCompoHeap::without_early_end`), which the test then
    /// reads as not filled to its cutoff.
    pub(crate) ends_early: bool,
}

/// Whether NCBI's composition-based redo fetches (and reads the header of) a match whose
/// preliminary best E-value is `evalue`, `heaps` those of its batch at its turn (the ones
/// `BlastCompo_EarlyTermination` reads): not where the early termination test ends it. A match
/// without HSPs is not fetched (blast_kappa.c:3513-3523); the caller leaves it out.
pub(crate) fn redo_fetches(evalue: f64, heaps: &[RedoHeap]) -> bool {
    // `EVALUE_STRETCH` (redo_alignment.c:64).
    const EVALUE_STRETCH: f64 = 5.0;
    !heaps.iter().all(|heap| {
        heap.ends_early
            && heap.inserted_before >= heap.heap_threshold
            && evalue > EVALUE_STRETCH * heap.ecutoff
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3575 (cited at
// `blastp::blast_engine::database::exclude_fetched`)
/// The first of `fetched` (`FetchMarks::of_noting` of a redo's stream `stream`, in its order,
/// with the query of each list, its best E-value and whether it has HSPs) that NCBI's redo
/// (composition-based or Smith-Waterman, one loop over the stream) fetches, of the lists up to
/// the stream's place `until` (the match whose redo failed the batch, or every list): the
/// subjects whose turn its early termination test leaves (`redo_fetches`, with the heaps that
/// `heaps_at` makes from the lists of `stream` before the turn), kept or excluded. The one walk
/// of TBLASTN's redo (`stage_d_pipeline::note_redo_fetched`) and BLASTX's Smith-Waterman redo
/// (`kappa_batch::sw_fetched`).
pub(crate) fn first_redo_fetch<K: Copy + Eq + std::hash::Hash>(
    fetched: &[(FetchedUtf8, (usize, f64, bool))],
    stream: &[K],
    until: usize,
    heaps_at: impl Fn(&std::collections::HashSet<K>) -> Vec<RedoHeap>,
) -> Option<u32> {
    let mut before = std::collections::HashSet::new();
    let mut at = 0;
    for (read, (_, evalue, has_hsps)) in fetched {
        if read.index > until {
            break;
        }
        while at < read.index {
            before.insert(stream[at]);
            at += 1;
        }
        if *has_hsps && redo_fetches(*evalue, &heaps_at(&before)) {
            return Some(read.ordinal);
        }
    }
    None
}

impl FetchedUtf8 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream_mt_utils.c:145-170
    // and blast_hspstream.c:196-203 (`s_SortHSPListByOid`: decreasing OID, read from the end)
    // The ordinary traceback (blast_traceback.c:1505-1577) fetches the subjects of the batch's
    // preliminary hit lists one OID at a time, in increasing OID.
    /// The first of `fetched` in the ordinary traceback's order (the least ordinal), of every
    /// fetch with `check_oid_exclusion` that reads all of its subjects (also the ungapped
    /// re-evaluation of a nucleotide subject, in the order of the subject loop).
    pub(crate) fn traceback_first(fetched: &[FetchedUtf8]) -> Option<u32> {
        fetched.iter().map(|read| read.ordinal).min()
    }
}

/// NCBI's fetches with `check_oid_exclusion` of a database search (line DBID item lineless2):
/// under a negative seqid list with ids (`SeqDb::excludes_lineless_at_fetch`), a sequence that
/// the list keeps without any of its definition lines is excluded where the gapped traceback,
/// the composition-based redo or the ungapped re-evaluation of a nucleotide subject fetches it
/// (seqsrc_seqdb.cpp:302-313). It was searched and counts in the statistics; a gapped search's
/// preliminary hit lists hold it (it took its place there), and its hits leave the search
/// before the final hit lists are made, so they take no place in them. A program removes
/// such subjects from the preliminary hit lists right before its traceback (`retain`), after
/// their cut and trim (`ADAPTIVE_CBS`); for the composition-based redo the places of their
/// lists in its stream count (`marks`, `FetchMarks`), since they change NCBI's status of the
/// batch; an ungapped search of a nucleotide database drops a subject's HSPs before they
/// reach the preliminary hit lists. Shared by the threads of a search (`Sync`); each subject's
/// header is read once.
///
/// NCBI's fetch reads the header of every subject it fetches (`GetSeqIDs`), so a header that
/// cannot be read ends NCBI's search there (exit 255 in the traceback, `BLAST engine error`
/// in the ungapped re-evaluation) even where the subject would not be reported; LOSAT refuses
/// such a search there (`refusal`). The composition-based redo fetches only the matches that
/// its early termination test leaves (blast_kappa.c:3525-3540), so LOSAT, which checks every
/// subject of the preliminary hit lists, may refuse a search with such a header that NCBI ends
/// without reading it.
pub(crate) struct FetchExcluded<'a> {
    check: LineCheck<'a>,
    /// The first of the reads of `excludes` that took a UTF8String for a VisibleString, in the
    /// ordinary traceback's order (`FetchedUtf8::traceback_first`).
    utf8: crate::seqdb::Utf8Reads,
}

impl<'a> FetchExcluded<'a> {
    /// The check of `db`'s fetches; `None` where they exclude nothing (no negative seqid list
    /// with ids).
    pub(crate) fn of(db: &'a SeqDb) -> Option<Self> {
        db.excludes_lineless_at_fetch().then(|| FetchExcluded {
            check: LineCheck::new(db),
            utf8: Default::default(),
        })
    }

    /// Whether the fetch of `ordinal` excludes it; the reader's error where its header cannot
    /// be read (`refusal`). A header that holds a UTF8String where a VisibleString is expected
    /// is noted, in increasing ordinal (`take_utf8`): of a fetch that reads every subject with
    /// hits, as the ungapped re-evaluation of a nucleotide subject does.
    pub(crate) fn excludes(&self, ordinal: u32) -> Result<bool, SeqDbError> {
        self.excludes_noting(ordinal, &self.utf8)
    }

    /// `excludes`, with the read of a UTF8String header noted in `reads` (in increasing ordinal):
    /// a check that the batches of a search share and that are searched side by side, each with
    /// its own note (TBLASTN's `stage_d_pipeline::BatchFetches::excludes`).
    pub(crate) fn excludes_noting(
        &self,
        ordinal: u32,
        reads: &crate::seqdb::Utf8Reads,
    ) -> Result<bool, SeqDbError> {
        let read = self
            .check
            .read(ordinal)
            .map_err(|error| self.read_error(ordinal, error))?;
        if read.utf8 {
            reads.note_ordinal(self.check.db, u128::from(ordinal), ordinal);
        }
        Ok(!read.kept)
    }

    /// Notes in `reads` NCBI's fetch of the sequence `ordinal`, whose header holds a UTF8String
    /// where a VisibleString is expected (the batch's first such fetch: `FetchedUtf8`).
    pub(crate) fn note_fetched(&self, reads: &crate::seqdb::Utf8Reads, ordinal: u32) {
        reads.note_ordinal(self.check.db, 0, ordinal);
    }

    /// The first of the reads of `excludes` since the last call that took a UTF8String for a
    /// VisibleString (a batch takes its own after its search).
    pub(crate) fn take_utf8(&self) -> Option<crate::seqdb::SeqLocation> {
        self.utf8.take().map(|note| note.at)
    }

    /// `items` (preliminary hit lists, or HSP lists) without those on a subject the fetch
    /// excludes, in their order; `ordinal` names an item's subject (`None`: kept). The first
    /// header that cannot be read, in the order of `items`, is the error (with its ordinal),
    /// and `items` is left as it was. The items whose header holds a UTF8String where a
    /// VisibleString is expected are returned (`FetchedUtf8`: those that NCBI's fetches of the
    /// item read post its warning).
    pub(crate) fn retain<T>(
        &self,
        items: &mut Vec<T>,
        ordinal: impl Fn(&T) -> Option<u32>,
    ) -> Result<Vec<FetchedUtf8>, (u32, SeqDbError)> {
        self.retain_noting(items, ordinal, |_| ())
            .map(|fetched| fetched.into_iter().map(|(read, ())| read).collect())
    }

    /// `retain`, each `FetchedUtf8` with `noted(item)` of its item (what its redo's order and
    /// early termination test read of it, `redo_fetches`).
    pub(crate) fn retain_noting<T, N>(
        &self,
        items: &mut Vec<T>,
        ordinal: impl Fn(&T) -> Option<u32>,
        noted: impl Fn(&T) -> N,
    ) -> Result<Vec<(FetchedUtf8, N)>, (u32, SeqDbError)> {
        let (marks, utf8) = self.marks_noting(items, ordinal, noted)?;
        let mut index = 0;
        items.retain(|_| {
            index += 1;
            !marks.excluded(index - 1)
        });
        Ok(utf8)
    }

    /// Which of `items` (the preliminary hit lists of a redo's stream, in the stream's order)
    /// are on a subject the fetch excludes (`FetchMarks`); `ordinal` names an item's subject
    /// (`None`: fetched). The first header that cannot be read, in the order of `items`, is the
    /// error (with its ordinal), as `retain`'s. The items whose header holds a UTF8String where a
    /// VisibleString is expected come with them, as `retain_noting` returns them: the one read
    /// of the items' headers for both.
    pub(crate) fn marks_noting<T, N>(
        &self,
        items: &[T],
        ordinal: impl Fn(&T) -> Option<u32>,
        noted: impl Fn(&T) -> N,
    ) -> Result<(FetchMarks, Vec<(FetchedUtf8, N)>), (u32, SeqDbError)> {
        let mut excluded = Vec::with_capacity(items.len());
        let mut utf8 = Vec::new();
        let mut kept = 0;
        for (index, item) in items.iter().enumerate() {
            let keeps = match ordinal(item) {
                Some(ordinal) => {
                    let read = self
                        .check
                        .read(ordinal)
                        .map_err(|error| (ordinal, self.read_error(ordinal, error)))?;
                    if read.utf8 {
                        utf8.push((
                            FetchedUtf8 {
                                index,
                                kept: read.kept,
                                before: kept,
                                ordinal,
                            },
                            noted(item),
                        ));
                    }
                    read.kept
                }
                None => true,
            };
            kept += usize::from(keeps);
            excluded.push(!keeps);
        }
        Ok((FetchMarks(excluded), utf8))
    }

    /// The error of the read of `ordinal`'s header for `refusal`: a record that NCBI's read
    /// fails on (`SeqDbError::Ncbi`, the text of NCBI's exception, which names no file) is named
    /// by its header file, as the reader's other errors are.
    #[cold]
    #[inline(never)]
    fn read_error(&self, ordinal: u32, error: SeqDbError) -> SeqDbError {
        if !matches!(error, SeqDbError::Ncbi(_)) {
            return error;
        }
        match self.check.db.header_file(ordinal) {
            Some(path) => SeqDbError::Format {
                path,
                message: error.to_string(),
            },
            None => error,
        }
    }

    /// LOSAT's refusal (in `program`'s name) of a search whose fetch of `ordinal` cannot read
    /// its header (`error`), where NCBI's search fails (`FetchRefusal`); a header file whose
    /// bytes the system does not give is LOSAT's error as it is (`losat_reader_error`).
    pub(crate) fn refusal((ordinal, error): (u32, SeqDbError), program: &str) -> anyhow::Error {
        if error.is_unreadable() {
            return FetchRefusal(losat_reader_error(&error, program)).into();
        }
        // A header that LOSAT refuses where it is read (`seqdb::Volume::filtered_mask`) keeps
        // that refusal's words.
        if matches!(error, SeqDbError::Unsupported { .. }) {
            return FetchRefusal(database_error(error, program).to_string()).into();
        }
        FetchRefusal(
            refusal(
                &format!(
                    "under the negative seqid list NCBI reads the definition lines of each \
                     sequence that its search aligns again, and those of sequence {ordinal} \
                     cannot be read ({error})"
                ),
                program,
            )
            .to_string(),
        )
        .into()
    }
}

/// Which lists of a composition-based redo's stream, by their place in it, are on a subject
/// that NCBI's fetch excludes (`FetchExcluded::marks_noting`; none without a negative seqid
/// list with ids, `FetchMarks::default`). NCBI's redo loop tests such a list's early termination first,
/// like any other's; where the test does not stop the list, the fetch excludes the subject and
/// NCBI sets the status of the thread back to 0 and goes on with the next list, with no redo and
/// no heap insert, so a failed Smith-Waterman redo before it no longer fails the batch; where
/// the test stops the list, nothing changes (core's `SwBatchStatus::reaches_redo`, the order of
/// the two steps for every program's redo loop). One owner of the question "is this list's
/// subject excluded at the fetch", with two readers: TBLASTN's `stage_d_pipeline` keeps the marked
/// lists in its stream and asks `SwBatchStatus::reaches_redo` at each (the X-drop redo has no
/// status, so a marked list is a list skipped there); BLASTP's `exclude_fetched` removes them and
/// keeps their places (`blastp::blast_engine::database::FetchExcludedMatch`), which its redo loop
/// reads after each query's loop through the same rule. The other programs remove such lists with
/// `FetchExcluded::retain`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:306-314
#[derive(Debug, Default)]
pub(crate) struct FetchMarks(Vec<bool>);

impl FetchMarks {
    /// The marks of `items` under `excluded` (`FetchExcluded::marks_noting`), with the items
    /// whose header holds a UTF8String where a VisibleString is expected; none without it.
    pub(crate) fn of_noting<T, N>(
        excluded: Option<&FetchExcluded<'_>>,
        items: &[T],
        ordinal: impl Fn(&T) -> Option<u32>,
        noted: impl Fn(&T) -> N,
    ) -> Result<(Self, Vec<(FetchedUtf8, N)>), (u32, SeqDbError)> {
        match excluded {
            Some(excluded) => excluded.marks_noting(items, ordinal, noted),
            None => Ok((Self::default(), Vec::new())),
        }
    }

    /// Whether the fetch excludes the subject of the stream's item `index`.
    #[inline(always)]
    pub(crate) fn excluded(&self, index: usize) -> bool {
        self.0.get(index).copied().unwrap_or(false)
    }
}

/// LOSAT's refusal of `FetchExcluded::refusal`, a type of its own so that a program that
/// reports the other errors of its search as NCBI's engine errors passes it on as it is.
#[derive(Debug)]
pub(crate) struct FetchRefusal(String);

impl std::fmt::Display for FetchRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for FetchRefusal {}

/// NCBI's Seq-aligns of the final hit lists of a database search (line DBID item lineless):
/// a hit on a sequence that the database's seqid list keeps without any of its definition
/// lines gets none (`SeqDb::keeps_line`, blast_seqalign.cpp:1582-1586), so it leaves every
/// report, the prefetch, the taxonomy lookups and psiblast's PSSM, while it was searched, took
/// its place in the hit list before the cut to `-max_target_seqs` and counts in the
/// statistics. Each program drops such hit lists right after that cut, before the reads of
/// the reported sequences; a check reads each sequence's header once (one check per query
/// batch or search call; the reads of the reports decode the kept ones again). Under a
/// negative seqid list with ids the search's fetches have excluded them before
/// (`FetchExcluded`), except a protein subject of an ungapped search, which no fetch checks.
pub(crate) struct SeqAligned<'a> {
    check: LineCheck<'a>,
    /// The first read that took a UTF8String for a VisibleString, before any read that failed
    /// (`first_utf8`).
    utf8: Option<crate::seqdb::SeqLocation>,
    /// Whether a read failed: NCBI's Seq-aligns stop there, so no later read posts.
    failed: bool,
}

impl<'a> SeqAligned<'a> {
    /// The check of `db`'s hits; `None` where no hit is dropped here: no seqid list removes
    /// lines, or the search's fetches have excluded the sequences without a line before
    /// (`fetched`, the caller's: the search fetches its subjects with `check_oid_exclusion`, a
    /// gapped search or an ungapped one of a nucleotide database; and the database's list is
    /// negative with ids, `SeqDb::excludes_lineless_at_fetch`).
    pub(crate) fn of(db: &'a SeqDb, fetched: bool) -> Option<Self> {
        (db.drops_lineless_hits() && !(fetched && db.excludes_lineless_at_fetch())).then(|| {
            SeqAligned {
                check: LineCheck::new(db),
                utf8: None,
                failed: false,
            }
        })
    }

    /// `items` (hit lists, or the HSPs of a batch) without those on a subject that gets no
    /// Seq-align, in their order; `ordinal` names an item's subject (`None`: no database
    /// ordinal, kept). The calls of a batch give its items in NCBI's order of the Seq-aligns
    /// (query by query, each in hit-list order), so that the first read that took a
    /// UTF8String for a VisibleString is NCBI's (`first_utf8`).
    pub(crate) fn retain<T>(&mut self, items: &mut Vec<T>, ordinal: impl Fn(&T) -> Option<u32>) {
        // A header that cannot be read counts as kept: the reads of the reports then give its
        // error, as before (NCBI's `GetSeqIDs` fails at the same Seq-align).
        let (check, utf8, failed) = (&self.check, &mut self.utf8, &mut self.failed);
        items.retain(|item| {
            ordinal(item).is_none_or(|ordinal| match check.read(ordinal) {
                Ok(read) => {
                    if read.utf8 && utf8.is_none() && !*failed {
                        *utf8 = check.db.location_of(ordinal);
                    }
                    read.kept
                }
                Err(_) => {
                    *failed = true;
                    true
                }
            })
        });
    }

    /// The batch's first read of a header record (in NCBI's order of the Seq-aligns, before a
    /// read that fails) that took a UTF8String for a VisibleString: NCBI posts its warning
    /// there, before the reads of the reports (`HeaderReader::note_earlier`), also where the
    /// sequence gets no Seq-align.
    pub(crate) fn first_utf8(&self) -> Option<crate::seqdb::SeqLocation> {
        self.utf8
    }
}

/// The reader of the headers of a search's reported database sequences (`read_header`): the
/// database, with NCBI's lookup of the Seq-align ids of those with parsed Seq-ids
/// (`SeqDb::seqid_lookup`), made at the first of them and kept for the reader's other reads.
pub(crate) struct HeaderReader<'a> {
    db: &'a SeqDb,
    /// The first of the batch's reads, in NCBI's order, that met a UTF8String where a
    /// VisibleString is expected (`first_utf8`): the reader's own and those before it
    /// (`note_earlier`), in the order of `Utf8Phase`.
    utf8: crate::seqdb::Utf8Reads,
    /// The first queries of NCBI's query batches among the reader's queries, ascending
    /// (`batch_starts`; none: one batch), which order its reads before their phases.
    batches: std::cell::OnceCell<Vec<usize>>,
    lookup: std::cell::OnceCell<crate::seqdb::lookup::SeqidLookup<'a>>,
    /// Whether NCBI's formatters may read some reported sequences one by one
    /// (`prefetch_may_be_off`), where the weak matches of its scope are looked for
    /// (`weakly_held`), which the reports refuse where a batch that reports the sequence has
    /// its prefetch off (`ParsedShown::one_by_one`).
    prefetch_off: bool,
    /// The taxonomy database of the search path where `taxdb.bti` is found and it opens, opened
    /// at the first record whose title NCBI can make from its taxids (`source_taxname`).
    taxdb: std::cell::OnceCell<Option<crate::seqdb::taxdb::TaxDb>>,
    /// What `source_taxname` looked up, as NCBI's volume keeps it (`seqdb::taxdb::TaxCache`).
    source_names: crate::seqdb::taxdb::TaxCache,
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:216-221
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

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:727-753
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

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:764-772 (cited at
// `prefetch_stop`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-904 (cited at
// `report::tabular::TabularSpec::subject_fetch`): where nothing was loaded before, a row reads
// the subject's Bioseq only for the fields that need it.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-813
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1216-1217 (the pairwise
// report's alignments; its description lines likewise, 1550-1551)
/// What NCBI's BLAST database data loader loaded for a query batch (`RegisterIds`; the run's
/// record, [`TaxdbRun::note_loaded`]), as [`batch_loads`] tells it: the one answer to whether
/// the batch's prefetch ran.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BatchLoads {
    /// The prefetch ran: every subject of the batch's Seq-aligns, those that `s_PreFetchSeqs`
    /// counts (psiblast's traceback hit lists, before the formatter's cut to `-max_target_seqs`).
    Prefetched,
    /// `PRE_FETCH_SEQS_LIMIT` turned the prefetch off ([`prefetch_off`]): the subjects that the
    /// reports NCBI writes read, after the formatter's cut to the hitlist size (`PruneSeqalign`),
    /// by the first format's rows (`TabularSpec::subject_fetch`; `Every` for a format that is not
    /// a field list), and those that psiblast's PSSM reads; a report that NCBI's XML stream cut
    /// before its hits loads nothing ([`XmlCuts`], [`LoadedBefore::subjects`]).
    Read(crate::report::tabular::SubjectFetch),
}

impl BatchLoads {
    /// The rows that these loads read of the hit lists that a batch's reports show, in a program
    /// whose prefetched set is its reported set (BLASTN, TBLASTN, BLASTX, TBLASTX): every subject
    /// where the prefetch ran, else the first format's rows (`Read`). Not psiblast's, whose
    /// prefetch reads its traceback beyond the reports' cut and whose PSSM reads more
    /// (`blastp::blast_engine::report_prep`). The one owner of that rule (review 7 of line DBID
    /// loadord, L2).
    pub(crate) fn reported_rows(self) -> crate::report::tabular::SubjectFetch {
        match self {
            BatchLoads::Prefetched => crate::report::tabular::SubjectFetch::Every,
            BatchLoads::Read(rows) => rows,
        }
    }

    /// The queries of the batch of the queries `span` (by query of the call, as `cuts` counts
    /// them) that these loads are made for: the batch at once where the prefetch ran (it loads the
    /// subjects of the batch's Seq-aligns before the reports); where it is off, each query whose
    /// report NCBI writes, alone: not one whose outfmt 5 report NCBI's XML stream cut before its
    /// hits (`cuts`), which loads nothing. The one walk of a batch's reports for the run's record
    /// ([`LoadedBefore::subjects`], BLASTX's `ReportedSubjects::note_loaded`; review 7 of line DBID
    /// loadord, L2).
    pub(crate) fn report_spans<'c>(
        self,
        cuts: XmlCuts<'c>,
        span: std::ops::Range<usize>,
    ) -> impl Iterator<Item = std::ops::Range<usize>> + 'c {
        let (batch, each) = match self {
            BatchLoads::Prefetched => (Some(span), None),
            BatchLoads::Read(_) => (None, Some(span)),
        };
        batch.into_iter().chain(
            each.into_iter()
                .flatten()
                .filter(move |&query| !cuts.cut(query))
                .map(|query| query..query + 1),
        )
    }
}

/// What NCBI's BLAST database data loader loads for a query batch with `seq_aligns` Seq-aligns
/// ([`BatchLoads`]): `Prefetched` where its prefetch is on; where `PRE_FETCH_SEQS_LIMIT` turns it
/// off ([`prefetch_off`]), what its reports read, `rows` (the first format's), which each program
/// applies to the hit lists that its reports show (review 3 of line DBID loadord, H2; review 6,
/// H1: psiblast's prefetched set holds its traceback beyond the reports' cut).
pub(crate) fn batch_loads(
    seq_aligns: usize,
    rows: crate::report::tabular::SubjectFetch,
) -> BatchLoads {
    batch_loads_at(prefetch_seqs_limit(), seq_aligns, rows)
}

/// [`batch_loads`] with the limit `limit` (`prefetch_seqs_limit`).
fn batch_loads_at(
    limit: Option<i32>,
    seq_aligns: usize,
    rows: crate::report::tabular::SubjectFetch,
) -> BatchLoads {
    if prefetch_off_at(limit, seq_aligns) {
        BatchLoads::Read(rows)
    } else {
        BatchLoads::Prefetched
    }
}

/// [`prefetch_off`] with the limit `limit` (`prefetch_seqs_limit`).
fn prefetch_off_at(limit: Option<i32>, seq_aligns: usize) -> bool {
    limit.is_some_and(|limit| {
        limit != i32::MAX
            && (limit == 0 || i64::try_from(seq_aligns).unwrap_or(i64::MAX) > i64::from(limit))
    })
}

/// The query batches of a call before a batch's prefetch that the call holds, for the run's record
/// of what NCBI's BLAST database data loader loaded for them ([`batch_prefetch`]'s `before`,
/// [`LoadedBefore::subjects`]): the hits of the call (`hits`, in each program's form: BLASTP's,
/// BLASTN's and TBLASTN's final hit lists by query, TBLASTX's reported hits in report order), the
/// queries of each batch (`spans`), which rows of the first format read a subject where a
/// batch's prefetch is off (`rows`, `TabularSpec::subject_fetch`), and the outfmt 5 reports that
/// NCBI's XML stream cuts before their hits (`cuts`, the call's [`XmlCuts`], which the texts that
/// they leave ask too: [`XmlReportsBefore`]).
pub(crate) struct LoadedBefore<'a, T> {
    pub(crate) hits: &'a [T],
    pub(crate) spans: &'a [std::ops::Range<usize>],
    pub(crate) rows: crate::report::tabular::SubjectFetch,
    pub(crate) cuts: XmlCuts<'a>,
}

// By hand: a derive would ask `T: Copy` of the hits it only borrows.
impl<T> Clone for LoadedBefore<'_, T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for LoadedBefore<'_, T> {}

impl<'a, T> LoadedBefore<'a, T> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:764-766 (cited above
    // `batch_loads`: no prefetch, the reports load what they read)
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375 (cited at
    // `report::query_warnings::ReadCursor::cut`: an `Iteration` that the stream cut makes no
    // `Hit`, so it loads nothing)
    /// What NCBI's BLAST database data loader loaded for each held batch, in the order of the
    /// batches: `loaded(span, hits, loads)` of the batch of the queries `span`, whose hits
    /// `batch` gives (`None`: the call does not hold them) and whose loads are [`batch_loads`] of
    /// their Seq-aligns (`seq_aligns(span, hits)`, [`subject_seq_aligns`]); where they are the
    /// reports' reads ([`BatchLoads::Read`]), of each query of the batch whose report NCBI writes
    /// (`span` its query alone; not one whose outfmt 5 report the stream cut before its hits,
    /// `cuts`, which loads nothing: review 6 of line DBID loadord, H2). The walk over the
    /// batches and the one question whether a batch's prefetch is off are here, which of its
    /// reports load is [`BatchLoads::report_spans`]; how each program counts its Seq-aligns and
    /// which of its subjects a batch's loads read stay with its hit lists (review 4 of line DBID
    /// loadord, L4; [`BatchLoads::reported_rows`] where its prefetched set is its reported set).
    pub(crate) fn subjects<I>(
        self,
        batch: impl Fn(&'a [T], std::ops::Range<usize>) -> Option<&'a [T]> + 'a,
        seq_aligns: impl Fn(std::ops::Range<usize>, &'a [T]) -> usize + 'a,
        loaded: impl Fn(std::ops::Range<usize>, &'a [T], BatchLoads) -> I + 'a,
    ) -> impl Iterator<Item = I::Item> + 'a
    where
        I: IntoIterator + 'a,
        I::IntoIter: 'a,
    {
        self.subjects_at(prefetch_seqs_limit(), batch, seq_aligns, loaded)
    }

    /// [`LoadedBefore::subjects`] with the limit `limit` (`prefetch_seqs_limit`).
    fn subjects_at<I>(
        self,
        limit: Option<i32>,
        batch: impl Fn(&'a [T], std::ops::Range<usize>) -> Option<&'a [T]> + 'a,
        seq_aligns: impl Fn(std::ops::Range<usize>, &'a [T]) -> usize + 'a,
        loaded: impl Fn(std::ops::Range<usize>, &'a [T], BatchLoads) -> I + 'a,
    ) -> impl Iterator<Item = I::Item> + 'a
    where
        I: IntoIterator + 'a,
        I::IntoIter: 'a,
    {
        let Self {
            hits,
            spans,
            rows,
            cuts,
        } = self;
        spans
            .iter()
            .flat_map(move |span| {
                let Some(held) = batch(hits, span.clone()) else {
                    return Vec::new();
                };
                let loads = batch_loads_at(limit, seq_aligns(span.clone(), held), rows);
                loads
                    .report_spans(cuts, span.clone())
                    .filter_map(|queries| Some((queries.clone(), batch(hits, queries)?, loads)))
                    .collect()
            })
            .flat_map(move |(span, hits, loads)| loaded(span, hits, loads))
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1417-1446
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

/// The limit of `PRE_FETCH_SEQS_LIMIT` (`app::prefetch_seqs_limit`); `None` when it is not set. A
/// value that NCBI cannot convert stops the first batch before its prefetch
/// (`app::check_prefetch_environment`), so no batch's prefetch is decided by it: the plans made
/// before the search count it as unset.
fn prefetch_seqs_limit() -> Option<i32> {
    crate::blastinput::app::prefetch_seqs_limit().ok().flatten()
}

/// Whether a Seq-id of the sequence `ordinal` of `db` is one `CSeq_id_Handle` with the key
/// `handle` (`seqid::mask_handle::mask_handle_key`).
fn holds_handle(db: &SeqDb, ordinal: u32, handle: &[u8]) -> bool {
    let Ok(crate::seqdb::SequenceHeader::Parsed(set)) = db.read_header(ordinal) else {
        return false;
    };
    set.deflines
        .iter()
        .flat_map(|defline| &defline.ids)
        .any(|id| {
            crate::seqid::mask_handle::mask_handle_key(id)
                .is_ok_and(|key| key.as_deref() == Some(handle))
        })
}

impl<'a> HeaderReader<'a> {
    /// The reader of the headers of `db`.
    pub(crate) fn new(db: &'a SeqDb) -> Self {
        Self {
            db,
            utf8: Default::default(),
            batches: std::cell::OnceCell::new(),
            lookup: std::cell::OnceCell::new(),
            prefetch_off: prefetch_may_be_off(),
            taxdb: std::cell::OnceCell::new(),
            source_names: Default::default(),
            lookup_texts: std::cell::OnceCell::new(),
        }
    }

    /// The first of the batch's reads, in NCBI's order, of a record that holds a UTF8String where
    /// a VisibleString is expected, with the query it is posted with (`Utf8Read`): the reads of
    /// the batch's fetches and Seq-aligns (`note_earlier`), then the reader's own (`read_header`,
    /// in their order: the programs' `read_titles` read in NCBI's). The scope's reads come after
    /// them, at the loads (`ParsedSubject::scope_utf8`). NCBI's warning names it; `None` where
    /// the run has posted it.
    pub(crate) fn first_utf8(&self) -> Option<Utf8Read> {
        let note = self.utf8.first()?;
        if self.db.utf8_posted() {
            return None;
        }
        Utf8Read::of_note(self.db, note)
    }

    /// `first_utf8` of the reads that NCBI makes up to the read of `phase` by `query` (a read
    /// that fails ends NCBI's Seq-aligns, before the scope reads anything).
    fn first_utf8_before(&self, phase: Utf8Phase, query: usize) -> Option<Utf8Read> {
        let note = self.utf8.first()?;
        if note.order > self.order(phase, query) || self.db.utf8_posted() {
            return None;
        }
        Utf8Read::of_note(self.db, note)
    }

    /// Sets the first queries of NCBI's query batches among the reader's queries (`batches`):
    /// a search call that holds several batches reads each batch's fetches, then its
    /// Seq-aligns, before the next batch's (`Utf8Phase`).
    pub(crate) fn batch_starts(&self, starts: Vec<usize>) {
        let _ = self.batches.set(starts);
    }

    /// The place in NCBI's order of a read of `phase` by `query`: its batch's (`batches`), then
    /// its phase's, then its query's.
    fn order(&self, phase: Utf8Phase, query: usize) -> u128 {
        let batch = self.batches.get().map_or(0, |starts| {
            starts[..starts.partition_point(|&start| start <= query)]
                .last()
                .copied()
                .unwrap_or(0)
        });
        phase.order(batch, query)
    }

    /// Notes the reads of the batches' fetches and Seq-aligns that come before those of the
    /// reader in NCBI's order and that took a UTF8String for a VisibleString (`EarlierUtf8`).
    pub(crate) fn note_earlier(&self, earlier: EarlierUtf8) {
        for (phase, read) in [
            (Utf8Phase::Fetch, earlier.fetched),
            (Utf8Phase::Aligned, earlier.aligned),
        ] {
            if let Some((query, at)) = read {
                self.utf8.note(self.order(phase, query), at);
            }
        }
    }

    /// Notes the read of the sequence `ordinal` by query `query`, whose record holds a
    /// UTF8String or whose read failed (`first_utf8`), where it is the first.
    #[cold]
    #[inline(never)]
    fn note_utf8(&self, ordinal: u32, query: usize) {
        if !self.db.utf8_posted() {
            self.utf8
                .note_ordinal(self.db, self.order(Utf8Phase::Reported, query), ordinal);
        }
    }

    /// The first read of the lookup `lookup` that took a UTF8String for a VisibleString
    /// (`SeqidLookup::take_utf8`: the version check of an accession), where the run has not
    /// posted NCBI's warning: the sequence read, which the load of the looked-up Seq-align id
    /// reads with the warning (`ParsedSubject::scope_utf8`).
    fn lookup_utf8(
        &self,
        lookup: &crate::seqdb::lookup::SeqidLookup<'_>,
    ) -> Option<crate::seqdb::DbSequence> {
        let note = lookup.take_utf8()?;
        (!self.db.utf8_posted()).then(|| self.scoped_utf8(self.db.sequence_at(note.at)))
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1169
    // The Bioseq gets one source for each positive taxid, in the order of the definition lines
    // that `x_GetTaxDefline` gives (the one that holds the Seq-align id first), and
    // `CDeflineGenerator::x_SetFlags` takes the first (create_defline.cpp:1004-1011). An empty
    // scientific name is no name (`m_Taxname` empty). Without `taxdb.bti` on the search path
    // `GetTaxNames` finds nothing; a database that does not open names nothing either (its
    // texts are `LookupTexts`'). A lookup whose scientific name NCBI reads past its mapping of
    // `taxdb.btd` (`Err`: an empty `taxdb.btd`, where NCBI faults, or a record without its third
    // tab inside the mapping, where NCBI's result depends on its memory map) is refused where the
    // title is shown.
    /// The name that the taxonomy database gives the taxid `taxid` of a Bioseq's first source
    /// (`DeflineSet::first_source_taxid`: of a stored title its definition line's, of a sequence
    /// with parsed Seq-ids `ParsedSubject::first_source_taxid`; none: `Unnamed`), its titles made
    /// from (`report::defline::BioseqFacts::taxname`). The lookups are kept as NCBI's volume keeps
    /// them (`source_names`); their texts are `TaxdbRun`'s and `LookupTexts`'.
    fn source_taxname(&self, taxid: Option<i32>) -> crate::report::defline::SourceTaxname {
        use crate::report::defline::SourceTaxname;
        let Some(taxid) = taxid else {
            return SourceTaxname::Unnamed;
        };
        let Some(taxdb) = self.taxdb() else {
            return SourceTaxname::Unnamed;
        };
        use crate::seqdb::taxdb::TaxLookup;
        let lookup = |messages: &mut Vec<String>| {
            match taxdb.lookup(taxid, messages) {
            Ok(TaxLookup::Named(names)) => Ok(Some(names)),
            Ok(TaxLookup::NotInIndex | TaxLookup::Unnamed) => Ok(None),
            Ok(TaxLookup::PastMapping(_)) => Err(
                "a title that NCBI BLAST+ makes from its taxid with the taxonomy database (an empty one, or the end of a protein title ending in ']'), whose lookup of that taxid reads the record past NCBI BLAST+'s mapping of taxdb.btd (NCBI BLAST+ 2.17.0 reads memory past its mapping there: its bytes, or whether it faults, depend on its memory map)",
            ),
            Err(_) => Err(
                "a title that NCBI BLAST+ makes from its taxid with the taxonomy database (an empty one, or the end of a protein title ending in ']'), whose lookup of that taxid reads a record of an empty taxdb.btd (NCBI BLAST+ 2.17.0 maps the file as a null pointer and reads unmapped memory)",
            ),
        }
        };
        match self.source_names.names(taxid, &mut Vec::new(), lookup) {
            Ok(names) => SourceTaxname::of_names(&names),
            Err(reason) => SourceTaxname::Refused(reason),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // The Bioseq of a sequence without parsed Seq-ids has its definition line's title and a
    // source for its taxid as one with them has (`CSeqDBVol::x_GetTaxonomy`), so NCBI makes an
    // empty title from the source too: `Homo sapiens` for a nucleotide and `unnamed protein
    // product [Homo sapiens]` for a protein with taxid 9606 (outfmt 0's description and heading,
    // outfmt 18, the first word of the tabular subject ids), where the taxonomy database names
    // the taxid (review 1 of line DBID N1b L1, `blastdb-ids-20261021/R/titlesrc/review1/r2.tsv`;
    // probe of line RG T5, `remaining-gaps-20261112/T/t5`).
    // A protein title ending in `]` loses its end to `x_AdjustProteinTitleSuffix` where the
    // Bioseq has a source, named or not (`report::defline::adjust_protein_title_suffix`): the
    // first word that replaces `gnl|BL_ORD_ID|n` in the tabular ids is the adjusted title's
    // (`alpha` of `alpha, partial [Bacillus sp.]` with taxid 100002; `chi` of `chi, partial [ ]`
    // also without a name: probe of line RG T5c, `T/t5c/out/ncbi/psall.blastp.f6.*`). The other
    // reports show the stored title, or the description's without the adjustment. So the source
    // of such a title is asked only where the first word can depend on it: most titles ending in
    // `]` (`hypothetical protein [Bacillus sp.]`) keep their first word whatever the source
    // (`report::defline::protein_title_word_needs_source`) and ask nothing, as for any other title.
    // NCBI looks the taxid up where it loads the Bioseq, whatever the title: that lookup, and its
    // refusal of a record of an empty `taxdb.btd` (NCBI faults), are the prefetch's and the XML2
    // reports' (`TaxdbRun::load`).
    /// The source of the sequence `ordinal` without parsed Seq-ids, whose stored title is `title`
    /// and whose definition line's taxid is `taxid` (`StoredTitle::taxid`, read with the title;
    /// one definition line, so its positive taxid is the first source's,
    /// `DeflineSet::first_source_taxid` without a target: NCBI reads the sequence by its ordinal
    /// id), where a title that NCBI makes depends on it: the name from which it makes an empty
    /// title (`report::defline::bioseq_title_from_source`, where a taxonomy database opens), and
    /// for a protein title ending in `]` whose first word its rewritten end can change
    /// (`report::defline::protein_title_word_needs_source`) whether there is a source
    /// (`SourceTaxname::Nameless` without a name) and its name (`source_taxname`). Other titles
    /// ask nothing (`Unnamed`).
    pub(crate) fn stored_source(
        &self,
        title: &[u8],
        taxid: i32,
        ordinal: u32,
    ) -> crate::report::defline::SourceTaxname {
        use crate::report::defline::SourceTaxname;
        let taxid = (taxid > 0).then_some(taxid);
        if let Some(taxid) = taxid.filter(|_| {
            self.db.seqtype() == SeqType::Protein
                && crate::report::defline::protein_title_word_needs_source(title)
        }) {
            return match self.source_taxname(Some(taxid)) {
                SourceTaxname::Unnamed => SourceTaxname::Nameless,
                source => source,
            };
        }
        if !crate::report::defline::bioseq_title_from_source(title, Default::default)
            || self.taxdb().is_none()
        {
            return SourceTaxname::Unnamed;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1976-1981
        // `CSeq_id::GetSeqIdString()` (Seq_id.cpp:2145-2154, `seqid::label::content`) of the
        // Seq-align's `gnl|BL_ORD_ID|n`, the OID (with a taxid list, of the ordinal id):
        // `seqdb::Oid::content_label`.
        self.source_taxname(taxid)
            .of_stored(|| self.db.oid(ordinal).content_label().into_bytes())
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:213-221
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:541-551
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
    /// (`SeqidLookup::seqid_to_oid`, which returns at the first OID found), and its scope can show
    /// another sequence that holds an id it takes for the Seq-align id (`defline_set::scope_keys`)
    /// where both are reported, which LOSAT lists by every OID of the lookup
    /// (`SeqidLookup::seqid_to_oids` with `multi`) and the keys in another letter case: where its
    /// list fails on a file that NCBI's lookup does not read (an ISAM of a volume after the first
    /// one that holds the id: a version 4 string ISAM that does not open; in version 5 the pages
    /// of the other keys), the other holders are unknown (`FoundInstead::HoldersUnknown`, with
    /// the sequence the lookup gives where it is another one). Another
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
    /// it whatever the lookup gives). A giim id of a version 5 database is looked up only for a
    /// lookup that fails: NCBI's handles never take another sequence for it, but its prefetch looks
    /// it up as a string id in every LMDB file of the list, which stops at a file without the
    /// accession table (review 1 of line DBID mixdb, H1); of a version 4 database never.
    /// `scope_utf8` gets the sequence that the lookup gives first where that is another one
    /// whose header holds a UTF8String where a VisibleString is expected (`scope_read`).
    fn found_instead(
        &self,
        subject: &crate::report::db_subject::ParsedSubject,
        ordinal: u32,
        scope_utf8: &mut Option<crate::seqdb::DbSequence>,
    ) -> (
        Option<crate::report::db_subject::FoundInstead>,
        Option<crate::report::db_subject::WeakMatch>,
        Option<crate::seqdb::DbSequence>,
    ) {
        use crate::report::db_subject::FoundInstead;
        let giim = matches!(subject.alignment_id(), crate::seqid::ParsedSeqId::Giim(_));
        if giim && !self.db.is_version5() {
            return (None, None, None);
        }
        let failed = |error: SeqDbError| {
            (
                Some(FoundInstead::Failed(
                    crate::report::db_subject::LookupFailure::of_error(&error),
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
        // NCBI's lookup (`SeqidToOid`), whose version check reads its candidates' headers, also
        // where it fails after them (`lookup_utf8`): the load of the id reads them first.
        let found = lookup.seqid_to_oid(id);
        *scope_utf8 = self.lookup_utf8(lookup);
        let first = match found {
            Ok(_) if giim => return (None, None, None),
            Ok(first) => first,
            Err(error) => {
                return match lookup.throws_on_id(id) {
                    Some(message) => (Some(FoundInstead::ThrowsOnId(message)), None, None),
                    None => failed(error),
                }
            }
        };
        let own = self.db.oid(ordinal).get();
        // The data loader reads the header of the sequence that this lookup gives first
        // (`scope_read`), whatever the scan below finds.
        if let Some(first) = first.filter(|&first| first != own && scope_utf8.is_none()) {
            *scope_utf8 = self.scope_read(self.db.ordinal_of(first));
        }
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
                    return (
                        Some(FoundInstead::CopyApart(self.db.sequence(other), cause)),
                        None,
                        None,
                    );
                }
                Some(_) => given = Some(FoundInstead::Sequence(self.db.sequence(other))),
                None => {}
            }
        }
        // Where LOSAT's lists of the other holders below fail on a file that NCBI's lookup did not
        // read, the other holders are unknown, and what the lookup gives is kept.
        let instead = match &given {
            Some(FoundInstead::Sequence(other)) => Some(*other),
            _ => None,
        };
        let unknown = |why: String, error: SeqDbError| {
            (
                Some(FoundInstead::HoldersUnknown {
                    instead,
                    // `why` is only written in LOSAT's own refusal text (`db_subject`).
                    why: format!(
                        "{why} ({})",
                        error.ncbi_message().map_or_else(
                            || error.to_string(),
                            |text| String::from_utf8_lossy(&text).into_owned()
                        )
                    ),
                }),
                None,
                None,
            )
        };
        // Every holder: the lookup goes on to the ISAMs of the volumes after the first one that
        // holds the id, which NCBI's never opens (review 3 of line DBID loadord, H3: a version 4
        // string ISAM that does not open). Its version checks there are LOSAT's reads, not
        // NCBI's: they note no warning (`SeqidLookup::take_utf8`, dropped).
        let holders = lookup.seqid_to_oids(id, true);
        let _ = lookup.take_utf8();
        let oids = match holders {
            Ok(oids) => oids,
            Err(error) => {
                return unknown(
                    "the other sequences that hold it, as a file of a volume after the first one that holds it, which NCBI BLAST+'s lookup does not read, fails".to_owned(),
                    error,
                )
            }
        };
        // LOSAT's own searches of the keys in another letter case (version 5: the upper-case key
        // of an id without a text accession, the other keys of a name or a `pir` or `prf` id):
        // where one fails, NCBI's exact lookup, which read the files before, still succeeds, so the
        // holders are unknown (review 1 of line DBID pircase, L1; review 1 of line DBID v4late, L2).
        let folded = match lookup.case_folded_oids(id) {
            Ok(oids) => oids,
            Err(error) => {
                return unknown(
                    "the sequences that hold it in another letter case".to_owned(),
                    error,
                )
            }
        };
        let variants = match lookup.case_variant_oids(id) {
            Ok(oids) => oids,
            Err(error) => {
                return unknown(
                    "the sequences that hold it in another letter case".to_owned(),
                    error,
                )
            }
        };
        let mut keys = Vec::new();
        crate::report::db_subject::scope_keys(subject, &mut keys);
        // A sequence found under a key in another letter case holds such an id where it holds one
        // of the Seq-align id's own handle: the key of a name is the same in every alternative and
        // release, which the handles keep apart (`gb||locusy` and `emb||LOCUSY`, `sp||NAMEY` and
        // `tr||namey`: review 1 of line DBID pircase, H1).
        let handle = (!variants.is_empty())
            .then(|| {
                crate::seqid::mask_handle::mask_handle_key(id)
                    .ok()
                    .flatten()
            })
            .flatten();
        // The first was compared above. Besides a sequence given instead, any other holder
        // (also a copy of the subject's record) is one the scope can read for it.
        let mut also_held = None;
        for (other, variant) in oids
            .iter()
            .chain(&folded)
            .map(|oid| (oid, false))
            .chain(variants.iter().map(|oid| (oid, true)))
            .filter(|&(&oid, _)| oid != own && Some(oid) != first)
            .map(|(&oid, variant)| (self.db.ordinal_of(oid), variant))
        {
            let held = if variant {
                handle
                    .as_deref()
                    .is_some_and(|handle| self.holds_handle(other, handle))
            } else {
                self.holds_key(other, &keys, crate::seqid::defline_set::scope_keys)
            };
            if !held {
                continue;
            }
            if given.is_some() {
                also_held = Some(self.db.sequence(other));
                break;
            }
            match copy(other) {
                Some(RecordCopy::Apart(cause)) => {
                    return (
                        Some(FoundInstead::CopyApart(self.db.sequence(other), cause)),
                        None,
                        None,
                    );
                }
                Some(_) => {
                    return (
                        Some(FoundInstead::AlsoHeld(self.db.sequence(other))),
                        None,
                        None,
                    )
                }
                None => {}
            }
        }
        // A failed lookup of the weak candidates (LOSAT's own, of the name alone) leaves them
        // unknown where the batch reads one by one (review 1 of line DBID weakbatch, L1).
        let weak = if self.prefetch_off {
            self.weakly_held(lookup, subject, ordinal, &oids)
                .unwrap_or_else(|error| {
                    Some(crate::report::db_subject::WeakMatch::Failed(
                        // Written only in LOSAT's own refusal text (`db_subject`).
                        error.ncbi_message().map_or_else(
                            || error.to_string(),
                            |text| String::from_utf8_lossy(&text).into_owned(),
                        ),
                    ))
                })
        } else {
            None
        };
        let found = given.or_else(|| {
            title_apart
                .map(|other| {
                    FoundInstead::CopyApart(
                        self.db.sequence(other),
                        crate::report::db_subject::CopyApart::Title,
                    )
                })
                .or(other_taxids.map(|other| FoundInstead::OtherTaxids(self.db.sequence(other))))
        });
        (found, weak, also_held)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:128-135
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:160-163
    // A version 5 database keys a text Seq-id as written, in its letter case, with no key for the
    // name of an id with an accession, and a `pir` or `prf` id by its FASTA text; the scope matches
    // names without case and across the alternatives (`CSeq_id_Mapper::GetMatchingHandles`), so the
    // lookup cannot list the holders of a name (review 1 of line DBID N2c,
    // `blastdb-ids-20261021/R/weak/REVIEW_weak_1.md` H1). A version 4 database's ISAM keys are in
    // lower case without the alternative.
    /// For `found_instead` where NCBI's formatters may read the reported sequences one by one:
    /// another sequence of the database (`DbSequence`) that holds an id the Seq-align id `id` of
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
            // LOSAT's own scan, not a read of NCBI's lookup: its version checks note nothing.
            let alone = lookup.seqid_to_oids(&name_alone, true);
            let _ = lookup.take_utf8();
            found.extend(alone?);
        }
        let own = self.db.oid(ordinal).get();
        Ok(found
            .iter()
            .filter(|&&oid| oid != own)
            .map(|&oid| self.db.ordinal_of(oid))
            .find(|&other| {
                self.holds_key(other, &keys, crate::seqid::id_scope::kind_keys)
                    && self.copy_of(subject, ordinal, other) == RecordCopy::Differs
            })
            .map(|other| WeakMatch::Held(self.db.sequence(other))))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393-406
    // The data loader resolves a Seq-align id to the first OID of `SeqidToOid` (seqdb.cpp:903-916,
    // `oids[0]`) and reads that sequence's header whatever it then keeps; it reads no other holder
    // of the id. So the header of another sequence that the lookup gives first is read wherever
    // the scope resolves the id: the prefetch and every format that reads the Bioseq, not the
    // archive (outfmt 11) nor the Seq-annots without the prefetch, psiblast's PSSM also (probe of
    // line RG item T5d fix 3, `T/t5d/fix3`: 580 NCBI runs; the other sequence after the subject
    // is never read).
    /// The read by NCBI's scope of the header of the sequence `other`, which the lookup of a
    /// reported sequence's Seq-align id gives first: `other` where its record holds a UTF8String
    /// where a VisibleString is expected and the run has not posted NCBI's warning. The load of
    /// the id posts it (`load_read_texts`), where NCBI's scope resolves the id: a batch's
    /// prefetch, or a report, or psiblast's PSSM, that reads the subject.
    #[cold]
    #[inline(never)]
    fn scope_read(&self, other: u32) -> Option<crate::seqdb::DbSequence> {
        if self.db.utf8_posted() {
            return None;
        }
        let at = self.db.location_of(other)?;
        self.db.utf8_string_path_at(at)?;
        Some(self.scoped_utf8(self.db.sequence(other)))
    }

    /// `other`, a sequence whose header the load of a reported sequence's Seq-align id reads with
    /// NCBI's warning, noted for the run (`SeqDb::note_scope_utf8`: its loads then look for it).
    #[cold]
    #[inline(never)]
    fn scoped_utf8(&self, other: crate::seqdb::DbSequence) -> crate::seqdb::DbSequence {
        self.db.note_scope_utf8();
        other
    }

    /// Whether a Seq-id of the sequence `ordinal` is one `CSeq_id_Handle` with the key `handle`
    /// (`seqid::mask_handle::mask_handle_key`).
    fn holds_handle(&self, ordinal: u32, handle: &[u8]) -> bool {
        holds_handle(self.db, ordinal, handle)
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
        // A UTF8String where a VisibleString is expected reads as that string (`DeflineSet::utf8`):
        // a record that differs in that alone is the same to every report (probe of line RG item
        // T5d fix 3: NCBI writes the same bytes, with its warning).
        let copy = if set_a.deflines == set_b.deflines {
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
                // The titles made from the sources are the same where the names are.
                let named = self.source_taxname(subject.first_source_taxid())
                    != self.source_taxname(copy.first_source_taxid());
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
    match &header {
        Ok(crate::seqdb::SequenceHeader::Parsed(set)) if set.utf8 => {
            headers.note_utf8(ordinal, query)
        }
        Ok(crate::seqdb::SequenceHeader::Ordinal(deflines))
            if deflines.first().is_some_and(|defline| defline.utf8) =>
        {
            headers.note_utf8(ordinal, query)
        }
        _ => {}
    }
    match header {
        Ok(crate::seqdb::SequenceHeader::Parsed(set)) => {
            crate::report::db_subject::ParsedSubject::new(set)
                .map_err(|reason| parsed_refusal(db.oid(ordinal), &reason, program))
                .map(|subject| {
                    let mut scope_utf8 = None;
                    let (found, weak, also_held) =
                        headers.found_instead(&subject, ordinal, &mut scope_utf8);
                    ReportedHeader::Parsed(std::sync::Arc::new(
                        subject
                            .with_source(db.seqtype() == SeqType::Protein, |taxid| {
                                headers.source_taxname(Some(taxid))
                            })
                            .showing_gis(db.reports_show_gis())
                            .under_positive_list(db.has_positive_list())
                            .found_instead(found)
                            .weak_match(weak)
                            .also_held(also_held)
                            .scope_utf8(scope_utf8),
                    ))
                })
        }
        Ok(crate::seqdb::SequenceHeader::Ordinal(deflines)) => {
            let mut stored = stored_title(Ok(deflines), || db.oid(ordinal), query, program)?;
            stored.source = headers.stored_source(&stored.title, stored.taxid, ordinal);
            Ok(ReportedHeader::Stored(stored))
        }
        Err(error) => {
            // NCBI's read posted the warning of a UTF8String before it failed (where the record
            // holds one before the failure), after the reads before it; the scope reads nothing
            // after a Seq-align fails.
            if headers.db.utf8_string_path(ordinal).is_some() {
                headers.note_utf8(ordinal, query);
            }
            Err(header_read_error(
                error,
                query,
                program,
                headers
                    .first_utf8_before(Utf8Phase::Reported, query)
                    .as_ref(),
            ))
        }
    }
}

/// With `-parse_deflines`, LOSAT's refusal of a reported database sequence whose Seq-align id is
/// the id of a query of the scope `queries` (`query_batch::QueryScopes::of_query` of a query that
/// reports it; a query whose ids the reports show: `-parse_deflines`, and the PSSM query of
/// `-in_pssm` whose Seq-id is not one local id): NCBI's scope finds the query instead, so its
/// reports show the query's title and letters for the database sequence (`dropin-survey-20261029/DI7/DECISIONS.md` 8). Every report
/// of a query batch asks for its own scope, also of a sequence whose header an earlier batch read.
/// `oid` gives the sequence's OID, which the refusal names.
pub(crate) fn refuse_id_of_a_query(
    subject: &crate::report::db_subject::ParsedSubject,
    oid: impl Fn() -> crate::seqdb::Oid,
    queries: &crate::blastinput::query_batch::QueryScope,
    program: &str,
) -> Result<()> {
    let keys = subject.query_scope_keys();
    let query = query_with_shown_ids(program);
    if queries.has_key(&keys.own) {
        return Err(parsed_refusal(
            oid(),
            &format!(
                "has the Seq-id of {query} in the query batch that reports it (in any with \
                 outfmt 5; NCBI BLAST+ shows the query for it)"
            ),
            program,
        ));
    }
    if queries.has_key(&keys.weak) {
        return Err(parsed_refusal(
            oid(),
            &format!(
                "has the Seq-align id {}, which NCBI BLAST+'s scope can match weakly to the id \
                 of {query} in the \
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

/// How the refusals of `refuse_id_of_a_query` and `refuse_display_id_of_a_query` name a query
/// whose ids the reports show, for `program` (its name in the refusals): one read with
/// `-parse_deflines`, or the PSSM query, which the user gives with `-in_pssm` (PSIBLAST, TBLASTN)
/// or, to PSIBLAST alone, as an alignment with `-in_msa` (line PS item PS6; the PSSM made of it
/// keeps the master sequence's id). The text of every other program is the one before PS6.
fn query_with_shown_ids(program: &str) -> &'static str {
    if program == "PSIBLAST" {
        "a query read with -parse_deflines (or of the PSSM query of -in_pssm or -in_msa)"
    } else {
        "a query read with -parse_deflines (or of the PSSM query of -in_pssm)"
    }
}

/// What the refusals of the query-anchored formats (`refuse_generated_query_id`,
/// `refuse_display_id_of_a_query`) call the id of a query that NCBI's scope finds: the local id
/// that LOSAT generates for a FASTA query read without `-parse_deflines` (`Query_N`), or the id of
/// the PSSM query of `-in_pssm` that a report shows as a local id (`lcl|myprot`).
const GENERATED_QUERY_ID: &str = "the id of a query (generated, or the PSSM's)";

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1588-1589, 432-433, 497-498
/// LOSAT's refusal of a database sequence that a description table shows (outfmt 0 to 4, the
/// first `-num_descriptions` subjects) whose display id (`ParsedSubject::display_id`, the id of its
/// description line) is, or is matched weakly by NCBI's scope to, the id of a query of the scope
/// `queries` (`QueryScopeKeys::display`): the description line looks its title up again by that
/// id, a second lookup after the one by the Seq-align id (`refuse_id_of_a_query`, every format),
/// and NCBI's scope finds the query, so NCBI BLAST+ shows the query's title there (`gi|779|lcl|abc`
/// and a query `lcl|abc`; review 3 of line RG item T3a). `queries` holds the ids of queries read
/// with `-parse_deflines` (or of the PSSM query of `-in_pssm`), or, with `generated`, the local
/// ids `Query_N` of queries read without it, which BLASTP refuses in outfmt 1 to 4 only (line RG
/// item X7 owns NCBI's behaviour for them in outfmt 0).
pub(crate) fn refuse_display_id_of_a_query(
    subject: &crate::report::db_subject::ParsedSubject,
    oid: crate::seqdb::Oid,
    queries: &crate::blastinput::query_batch::QueryScope,
    generated: bool,
    program: &str,
) -> Result<()> {
    if !queries.has_key(&subject.query_scope_keys().display) {
        return Ok(());
    }
    let mut display = Vec::new();
    if let Some(id) = subject.display_id() {
        crate::seqid::label::fasta(id, &mut display);
    }
    let shown = format!("the id of {}", query_with_shown_ids(program));
    let (queries, weak, formats) = if generated {
        (
            GENERATED_QUERY_ID,
            "",
            "the query-anchored formats, outfmt 1 to 4",
        )
    } else {
        (
            shown.as_str(),
            ", or an id that NCBI BLAST+'s scope matches weakly to it",
            "outfmt 0 to 4",
        )
    };
    Err(parsed_refusal(
        oid,
        &format!(
            "has the display id {} (the id of its description line), {queries} in the query \
             batch that reports it{weak}: NCBI BLAST+ shows that query's title in its \
             description line ({formats})",
            String::from_utf8_lossy(&display)
        ),
        program,
    ))
}

/// Outfmt 1 to 4 of queries read without `-parse_deflines` (line RG T3a, review 2 D2): LOSAT's
/// refusal of a reported database sequence `subject` whose Seq-align id is the generated id
/// (`lcl|Query_N`) of a query of the batch that reports it (`scope`, `QueryScope::of_generated_ids`,
/// as NCBI's scope finds an id: the keys of `refuse_id_of_a_query`). NCBI's formatter adds the
/// batch's queries to its scope before their reports and empties it after them and finds such a
/// subject there: its reports show the query's title and letters in the subject's place, and the
/// query-anchored merge of the query's own Bioseq can abort (`double free or corruption`). NCBI's
/// display of such a sequence in outfmt 0, 6 and 7 is line RG item X7's. `effect` says what NCBI
/// does with it: `QUERY_ANCHORED_EFFECT`, `XML_EFFECT` for the RPS programs' BLAST XML, XML2
/// and JSON (`-outfmt 5`, 13 to 16: RPSBLAST gets the query's title for the entry, RPSTBLASTN an
/// exception from its Seq-align conversion and half a document or none; PS15a-1, PS15a-2), or
/// `TABLE_EFFECT` for their outfmt 18 and 20 (PS15a-2).
pub(crate) fn refuse_generated_query_id(
    subject: &crate::report::db_subject::ParsedSubject,
    oid: crate::seqdb::Oid,
    scope: &crate::blastinput::query_batch::QueryScope,
    program: &str,
    effect: &str,
) -> Result<()> {
    if !scope.has_key(&subject.query_scope_keys().own) {
        return Ok(());
    }
    Err(parsed_refusal(
        oid,
        &format!(
            "has the Seq-id {}, {GENERATED_QUERY_ID} in the query batch that reports it: NCBI \
             BLAST+ {effect}",
            String::from_utf8_lossy(&subject.alignment_fasta())
        ),
        program,
    ))
}

/// What NCBI does with such a sequence in the query-anchored formats (`refuse_generated_query_id`).
pub(crate) const QUERY_ANCHORED_EFFECT: &str =
    "shows that query in its place in the query-anchored formats (outfmt 1 to 4)";

/// What NCBI does with it in the BLAST XML, XML2 and JSON documents of the RPS programs
/// (`refuse_generated_query_id`).
pub(crate) const XML_EFFECT: &str =
    "takes that query for it in BLAST XML, XML2 and JSON (outfmt 5, 13 to 16)";

/// What NCBI does with it in the organism report and in the comma-separated values with a header
/// line of the RPS programs where they differ for it (`report::tabular::TabularSpec::
/// refuses_generated_id`, `refuse_generated_query_id`).
pub(crate) const TABLE_EFFECT: &str = "takes that query for it in the organism report (outfmt \
     18) of a profile with a taxid and in comma-separated values with a header line (outfmt 20) \
     whose fields show what it finds by that id (RPSTBLASTN's outfmt 20 can stop with an \
     exception)";

/// The query-anchored formats (outfmt 1 to 4) of a program whose reports name a database's
/// sequences by rank (BLASTN, TBLASTN: line RG T3b; BLASTP's `ListedTitles` asks the same
/// refusals): for each query `0..num_queries`, the sequences of its final hit list (`reported`,
/// by rank) with parsed Seq-ids (`parsed`, by rank; `oid_of` a rank's OID), with the scope of the
/// query's batch (`scope_of`): queries read without `-parse_deflines` (`generated`, their
/// `lcl|Query_N`): `refuse_generated_query_id` of each; then, for a query that reports a sequence
/// whose display id its scope has (`QueryScopeKeys::display`), `refuse_display_id_of_a_query` of
/// each sequence whose description line its report shows (`described`, by rank: the description
/// line looks the title up again by that id and NCBI shows the query's title there). A query read
/// with `-parse_deflines` whose id a reported sequence has as its Seq-align id is refused in every
/// format by `refuse_reported_ids_of_queries` (line RG X7).
#[allow(clippy::too_many_arguments)]
pub(crate) fn refuse_query_ids_of_anchored_reports<'s>(
    num_queries: usize,
    mut reported: impl FnMut(usize) -> Vec<u32>,
    mut described: impl FnMut(usize) -> Result<Vec<u32>>,
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    oid_of: impl Fn(u32) -> crate::seqdb::Oid,
    scope_of: impl Fn(usize) -> &'s crate::blastinput::query_batch::QueryScope,
    generated: bool,
    program: &str,
) -> Result<()> {
    if parsed.is_empty() {
        return Ok(());
    }
    let parsed_of = |rank: u32| parsed.get(rank as usize).and_then(Option::as_deref);
    for query in 0..num_queries {
        let scope = scope_of(query);
        let ranks = reported(query);
        if generated {
            for &rank in &ranks {
                if let Some(subject) = parsed_of(rank) {
                    refuse_generated_query_id(
                        subject,
                        oid_of(rank),
                        scope,
                        program,
                        QUERY_ANCHORED_EFFECT,
                    )?;
                }
            }
        }
        if !ranks.iter().any(|&rank| {
            parsed_of(rank)
                .is_some_and(|subject| scope.has_key(&subject.query_scope_keys().display))
        }) {
            continue;
        }
        for rank in described(query)? {
            if let Some(subject) = parsed_of(rank) {
                refuse_display_id_of_a_query(subject, oid_of(rank), scope, generated, program)?;
            }
        }
    }
    Ok(())
}

/// `refuse_id_of_a_query` of each reported sequence of `reported` (a searched query, by its index
/// in `scopes`, and a sequence that its report names, by ordinal id; `oids` the database's OIDs)
/// whose header `titles` holds with parsed Seq-ids, with the scope of that query. A program calls it for each query batch
/// after it reads the batch's titles: its cache of the titles keeps those of the batches before.
pub(crate) fn refuse_reported_ids_of_queries(
    titles: &std::collections::BTreeMap<u32, ReportedHeader>,
    reported: impl IntoIterator<Item = (usize, u32)>,
    oids: crate::seqdb::OidMap<'_>,
    scopes: &crate::blastinput::query_batch::QueryScopes,
    program: &str,
) -> Result<()> {
    for (query, ordinal) in reported {
        if let Some(ReportedHeader::Parsed(subject)) = titles.get(&ordinal) {
            refuse_id_of_a_query(
                subject,
                || oids.oid(ordinal),
                scopes.of_query(query),
                program,
            )?;
        }
    }
    Ok(())
}

/// LOSAT's refusal of a reported database sequence with parsed Seq-ids (`parsed_refusal`,
/// `refuse_joined_hits`). BLASTX finds one only when it reports a
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

/// LOSAT's refusal of the database sequence `oid` with parsed Seq-ids for `reason` (what the
/// sequence has, `ParsedSubject`).
pub(crate) fn parsed_refusal(oid: crate::seqdb::Oid, reason: &str, program: &str) -> anyhow::Error {
    anyhow::Error::new(ParsedRefusal(format!(
        "database sequence {oid} {reason}; this is not supported by LOSAT's {program}"
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
/// sequence without parsed Seq-ids); `oids` the OIDs, by which the refusal names them.
pub(crate) fn refuse_joined_hits<'a>(
    query: usize,
    subjects: impl IntoIterator<Item = Option<(u32, &'a crate::report::db_subject::ParsedSubject)>>,
    oids: crate::seqdb::OidMap<'_>,
    program: &str,
) -> Result<()> {
    let mut previous: Option<(u32, &crate::report::db_subject::ParsedSubject)> = None;
    for subject in subjects {
        if let Some((before, ordinal)) = joined_apart(previous, subject) {
            return Err(anyhow::Error::new(ParsedRefusal(format!(
                "database sequences {} and {} follow each other in the hits of query {} and have Seq-align ids that NCBI BLAST+ takes for one sequence (its XML and outfmt 0 reports show their alignments as one hit, its tabular qcovs and qcovus count them as one) but are not copies of one record; this is not supported by LOSAT's {program}",
                oids.oid(before),
                oids.oid(ordinal),
                query + 1
            ))));
        }
        previous = subject;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1025-1033
/// NCBI's count of a cut to `cut` subjects (`PruneSeqalign`; `TrimSeqAlign` counts so too): the
/// number of the entries at the start of `entries` (one for each HSP list of a hit list, in its
/// order) that the cut keeps, where an entry that follows one of the same subject
/// (`joined(previous, entry)`: their Seq-align ids match) is counted with it. An error of
/// `entries` or `joined` (a caller's read, or its refusal of the pair) ends the count. The one
/// count of that cut, which the reports (`joined_cut`) and RPS-BLAST's check of the profiles that
/// the reports show (`rps_search::volumes::formatter_kept`) take.
pub(crate) fn pruned_kept<T, E>(
    cut: usize,
    entries: impl IntoIterator<Item = std::result::Result<T, E>>,
    mut joined: impl FnMut(&T, &T) -> std::result::Result<bool, E>,
) -> std::result::Result<usize, E> {
    let mut previous: Option<T> = None;
    let mut counted = 0usize;
    let mut kept = 0usize;
    for entry in entries {
        let entry = entry?;
        let joins = match &previous {
            Some(previous) => joined(previous, &entry)?,
            None => false,
        };
        if !joins {
            counted += 1;
            if counted > cut {
                break;
            }
        }
        kept += 1;
        previous = Some(entry);
    }
    Ok(kept)
}

/// The number of HSP lists of a hit list (query `query`; `subjects` and `oids`, as
/// `refuse_joined_hits` takes them) that its cut to `max_target_seqs` subjects keeps, as NCBI's
/// formatter cuts it (`pruned_kept`): following Seq-aligns whose Seq-align ids match count as
/// those of one subject. A copy of one record that follows another (`joins_as_copy`) is kept with
/// it; two other parsed subjects that the count reaches with matching Seq-align ids are LOSAT's
/// refusal. psiblast cuts its traceback hit lists so, which are longer than `-max_target_seqs`.
pub(crate) fn joined_cut<'a>(
    query: usize,
    max_target_seqs: usize,
    subjects: impl IntoIterator<Item = Option<(u32, &'a crate::report::db_subject::ParsedSubject)>>,
    oids: crate::seqdb::OidMap<'_>,
    program: &str,
) -> Result<usize> {
    pruned_kept(
        max_target_seqs,
        subjects.into_iter().map(Ok),
        |&previous, &subject| {
            if let Some((before, ordinal)) = joined_apart(previous, subject) {
                return Err(anyhow::Error::new(ParsedRefusal(format!(
                    "database sequences {} and {} follow each other in the hits of query {} that the cut to -max_target_seqs {max_target_seqs} counts and have Seq-align ids that NCBI BLAST+ takes for one sequence but are not copies of one record (it counts them as one sequence and keeps one more); this is not supported by LOSAT's {program}",
                    oids.oid(before),
                    oids.oid(ordinal),
                    query + 1
                ))));
            }
            Ok(matches!(
                (previous, subject),
                (Some((_, a)), Some((_, b))) if crate::report::db_subject::joins_as_copy(a, b)
            ))
        },
    )
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

/// The refusal of a description id that NCBI's cut of the description line leaves no room for
/// (`check_shown_parsed_subject`).
const LONG_DESCRIPTION_ID: &str = "has an outfmt 0 description id longer than 65 bytes (NCBI BLAST+ 2.17.0 pads the description line without end)";

/// LOSAT's refusals of the outfmt 0 description and alignment heading of the subject with
/// parsed Seq-ids `subject` (`oid` gives its OID, which a refusal names; `protein` for a protein
/// database)
/// that it does not write as NCBI: a protein title of punctuation that NCBI's `x_CleanAndCompress`
/// reads past the end of (NCBI BLAST+ 2.17.0 crashes or does not end; approved exception 2 of
/// PD-LOSAT-NCBI-DEFECTS has BLASTN, TBLASTN and TBLASTX stop such a nucleotide title at the end
/// of the string, as for a database without parsed Seq-ids), a description id longer
/// than NCBI's cut of the line leaves room for (NCBI then pads the line without end), a title
/// that NCBI makes from a source that LOSAT refuses (`ParsedSubject::source`,
/// `SourceTaxname::Refused`), and one that it makes from the
/// taxids of the same record that it can show for the sequence where `taxdb.bti` is found
/// (`taxdb_found`, `ParsedSubject::generated_title_fetched`). A title made from the sequence's
/// own source names its organism, and a generated heading title ending in `]` gets its end from
/// it (`ParsedSubject::with_source`, `report::defline::adjust_protein_title_suffix`).
pub(crate) fn check_shown_parsed_subject(
    subject: &crate::report::db_subject::ParsedSubject,
    protein: bool,
    taxdb_found: bool,
    oid: impl Fn() -> crate::seqdb::Oid,
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
                oid(),
                "has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0)",
                program,
            ));
        }
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:913-917
    // The unsigned length wraps for an id longer than 65 bytes: NCBI 2.17.0 then pads the
    // line with spaces without end (probe of the ID2a review, `tests/db_search/ids/refuse.faa`).
    if subject.description_start().len() > 65 {
        return Err(parsed_refusal(oid(), LONG_DESCRIPTION_ID, program));
    }
    if subject.pairwise_title_reads_source(protein, taxdb_found) {
        subject
            .generated_title_fetched()
            .map_err(|reason| parsed_refusal(oid(), &reason, program))?;
    }
    if let Err(reason) = subject.source().name() {
        return Err(parsed_refusal(oid(), &format!("has {reason}"), program));
    }
    Ok(())
}

/// [`check_shown_parsed_subject`] of a subject whose load in the pairwise reports is `load`
/// (`HitLoads::load`): the description line of one whose load throws shows its Seq-align id and
/// `Unknown` (`ParsedSubject::thrown_description`), with no heading, so that only that id's
/// length is refused there.
pub(crate) fn check_pairwise_parsed(
    subject: &crate::report::db_subject::ParsedSubject,
    load: crate::report::db_subject::HitLoad,
    protein: bool,
    taxdb_found: bool,
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<()> {
    if load != crate::report::db_subject::HitLoad::Throws {
        return check_shown_parsed_subject(subject, protein, taxdb_found, oid, program);
    }
    if subject.thrown_description().0.len() > 65 {
        return Err(parsed_refusal(oid(), LONG_DESCRIPTION_ID, program));
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
    /// Whether a tabular format can stop at a row whose subject lookups fail
    /// (`report::tabular::may_stop_at_failed_lookup`: `sacc`, `saccver`, the default fields, a
    /// first row's `staxid`), so that such a sequence may be refused (`reports_may_refuse`).
    pub(crate) stops_at_failed_lookup: bool,
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
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3289-3300
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
/// the subject's database; `load`, what the reports do at the subject's `Hit`
/// (`report::db_subject::HitLoads`, the one answer: its `Hit` throws, or no report loads it,
/// `ParsedSubject::fetched`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn check_reported_parsed(
    db: &SeqDb,
    subject: &crate::report::db_subject::ParsedSubject,
    shown: ParsedShown,
    load: crate::report::db_subject::HitLoad,
    protein: bool,
    taxdb_found: bool,
    ordinal: u32,
    program: &str,
) -> Result<()> {
    if shown.reads_bioseq() {
        subject
            .fetched(shown.taxonomy, shown.one_by_one, load)
            .map_err(|reason| parsed_refusal(db.oid(ordinal), &reason, program))?;
    } else if shown.prefetched {
        check_prefetched_parsed(db, subject, ordinal, program)?;
    }
    if shown.pairwise {
        check_pairwise_parsed(
            subject,
            load,
            protein,
            taxdb_found,
            || db.oid(ordinal),
            program,
        )?;
    }
    if shown.tabular {
        subject
            .check_tabular_ids(protein)
            .and_then(|()| subject.check_printed_word(protein, shown.words))
            .map_err(|problem| {
                parsed_refusal(db.oid(ordinal), &format!("has {problem}"), program)
            })?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40 (the fields of
// `-outfmt 6` without a list)
/// The subject id of the default tabular fields (`saccver`) of each reported database sequence by
/// rank, for BLASTN, TBLASTN and TBLASTX: the `saccver` of a sequence without parsed Seq-ids
/// (`TabularSubjectId::accession`: the first word of the title that NCBI makes of the stored one,
/// an empty one made from the name of its source, `sources[rank]`, or `BL_ORD_ID:n`), or of a
/// sequence with parsed Seq-ids (`parsed[rank]`, `ParsedSubject::tabular_ids`). The custom fields
/// ask `sseqid` of it (`tabular_seqid`). `records[rank]` holds the stored title, `ordinals[rank]`
/// the ordinal id, `oids` the database's OIDs; `protein` for a protein database. `None` where the
/// lookup fails (`common::tabular::SubjectId`). An id that LOSAT cannot write as NCBI is its
/// refusal where a tabular report `reads` the subject (`report::tabular::reads_subject`), else
/// empty: a table of fields that read nothing of the subject prints none of its ids.
pub(crate) fn tabular_subject_ids(
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    ordinals: &[u32],
    oids: crate::seqdb::OidMap<'_>,
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    sources: &[crate::report::defline::SourceTaxname],
    protein: bool,
    reads: bool,
    program: &str,
) -> Result<Vec<crate::algorithm::common::tabular::SubjectId>> {
    records
        .iter()
        .zip(ordinals)
        .enumerate()
        .map(|(rank, (record, &ordinal))| {
            let id = match parsed.get(rank).and_then(Option::as_deref) {
                Some(subject) => subject.tabular_saccver(protein).map_err(|problem| {
                    parsed_refusal(oids.oid(ordinal), &format!("has {problem}"), program)
                }),
                None => tabular_subject_id(
                    &record.title,
                    sources
                        .get(rank)
                        .unwrap_or(&crate::report::defline::UNNAMED_SOURCE),
                    oids.oid(ordinal),
                    protein,
                    program,
                )
                .map(|id| match id {
                    TabularSubjectId::Unknown(_) => None,
                    id => Some(std::sync::Arc::from(id.accession())),
                }),
            };
            match id {
                Err(_) if !reads => Ok(Some(std::sync::Arc::from(&b""[..]))),
                id => id,
            }
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
/// A reported sequence as a caller of the library reads it (`FormatObserver::database_subjects`):
/// the id of the default tabular report (`tabular_subject_ids`) of the sequence `ordinal` with the
/// stored title `title` and the source `source` of its Bioseq (with parsed Seq-ids: `parsed`), and
/// `length` letters. A sequence whose tabular id LOSAT refuses shows the first word of its title:
/// only a tabular report refuses that id, not the search.
#[allow(clippy::too_many_arguments)]
pub(crate) fn library_subject(
    db: &SeqDb,
    title: &[u8],
    source: &crate::report::defline::SourceTaxname,
    ordinal: u32,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    length: usize,
    protein: bool,
    program: &str,
) -> crate::api::local_blast::ReportedSubject {
    let id = match parsed {
        // No `saccver`: the lookup throws; the library names it as one without parsed ids.
        Some(subject) => subject
            .tabular_saccver(protein)
            .ok()
            .flatten()
            .map(|id| id.to_vec()),
        None => tabular_subject_id(title, source, db.oid(ordinal), protein, program)
            .ok()
            .filter(|id| !matches!(id, TabularSubjectId::Unknown(_)))
            .map(|id| id.accession().into_owned()),
    }
    .unwrap_or_else(|| {
        title
            .split(|&byte| byte == b' ')
            .next()
            .unwrap_or_default()
            .to_vec()
    });
    crate::api::local_blast::ReportedSubject {
        oid: db.oid(ordinal).get(),
        id: crate::api::Id::new(&id),
        length,
    }
}

/// `library_subject` of each reported sequence of a nucleotide database, by rank: `records` hold
/// the stored titles, `ordinals`, `parsed` and `sources` as in `tabular_subject_ids`.
pub(crate) fn library_subjects(
    db: &SeqDb,
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    ordinals: &[u32],
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    sources: &[crate::report::defline::SourceTaxname],
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
                sources
                    .get(rank)
                    .unwrap_or(&crate::report::defline::UNNAMED_SOURCE),
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
/// `StoredTitle::check_shown_source` and `check_reported_parsed` with `taxdb_found` where outfmt 0
/// is shown, as `reported_subjects` asks). Where the reports join following hits (`shown.joined`), no two sequences with parsed
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
    let lineless_dropped = db.drops_lineless_hits();
    let check = |kept: &mut Kept,
                 headers: &HeaderReader<'_>,
                 ordinal: u32,
                 header: std::result::Result<crate::seqdb::SequenceHeader, SeqDbError>|
     -> std::result::Result<(), ReportRefusal> {
        // A sequence left without a line gets no Seq-align (`SeqAligned`): none of its hits is
        // reported, so it cannot make a report refuse (an error of its header comes below).
        if lineless_dropped && matches!(db.keeps_line(ordinal), Ok(false)) {
            return Ok(());
        }
        match reported_header(headers, header, ordinal, 0, program) {
            Ok(ReportedHeader::Stored(stored)) => {
                check_reported_title(stored.replaced_bytes, || db.oid(ordinal), program)
                    .map_err(|_| ReportRefusal::Headers)?;
                // A failed lookup (`TabularSubjectId::Unknown`) is refused only by the field lists
                // that can stop at it (`report::tabular::failed_lookup_stop`).
                if shown.tabular {
                    match tabular_subject_id(
                        &stored.title,
                        &stored.source,
                        db.oid(ordinal),
                        protein,
                        program,
                    ) {
                        Ok(TabularSubjectId::Word(_) | TabularSubjectId::General(_)) => {}
                        Ok(TabularSubjectId::Unknown(_)) if !shown.stops_at_failed_lookup => {}
                        _ => return Err(ReportRefusal::Headers),
                    }
                }
                if shown.pairwise {
                    stored
                        .check_shown_source(|| db.oid(ordinal), program)
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
                // Before the search no `Hit` is known: every subject is taken as loaded.
                check_reported_parsed(
                    db,
                    &subject,
                    shown,
                    crate::report::db_subject::HitLoad::Loaded,
                    protein,
                    shown.pairwise && taxdb_found,
                    ordinal,
                    program,
                )
                .map_err(|_| ReportRefusal::Headers)?;
                // A sequence whose lookups fail (`ParsedSubject::lookup_failure`), where a field
                // list can stop at it.
                if shown.tabular
                    && shown.stops_at_failed_lookup
                    && !matches!(subject.lookup_failure(protein), Ok(None))
                {
                    return Err(ReportRefusal::Headers);
                }
                if let Some(queries) = queries {
                    refuse_id_of_a_query(&subject, || db.oid(ordinal), queries, program)
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
/// could be refused, as the batches' prefetch looks them up (`prefetch_messages`; its names reach
/// no output, so only a fault is refused: a record of an empty `taxdb.btd`, `TaxDb::lookup`);
/// `false` without `taxdb.bti` there. One pass over the index usually finds that no lookup of any
/// taxid can fault (`TaxDb::any_lookup_faults`); else each one is looked up once. Both use a taxonomy
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
    if matches!(taxdb.any_lookup_faults(), Ok(false)) {
        return false;
    }
    let mut messages = Vec::new();
    taxids
        .iter()
        .any(|&taxid| taxdb.lookup(taxid, &mut messages).is_err())
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1865-1886
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
/// LOSAT's rejection of a reported database sequence (one with rows in the final hit lists;
/// `oid` gives its OID, which the rejection names, where it rejects) whose stored title held a byte outside ' ' ..= '~' (`replaced_bytes`,
/// `StoredTitle`): NCBI writes an error message for the byte each time its formatter
/// reads the definition line.
pub(crate) fn check_reported_title(
    replaced_bytes: bool,
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<()> {
    if replaced_bytes {
        bail!(
            "database sequence {} has a stored title with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line); this is not supported by LOSAT's {program}",
            oid()
        );
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:169-182 (the ids of a
// subject's list; for a sequence without parsed Seq-ids, the local id of the title's first word
// that `GetSeqIdList` makes, or its general id, showdefline.cpp:218-242 at
// `report::defline::database_subject_id`)
/// The subject id of the tabular formats of a database sequence without parsed Seq-ids
/// (`tabular_subject_id`): the first word of the title that NCBI makes of its stored title (an
/// empty one from the name of its source), which `GetSeqIdList` puts in the place of its id, or
/// where that word is empty its general id `gnl|BL_ORD_ID|n`, n its OID (`seqdb::Oid`; NCBI
/// BLAST+ 2.17.0, line DBID ordname, `blastdb-ids-20261021/R/ordname`).
/// `Unknown` where `GetSeqIdList` makes no id: the generated title throws
/// (`report::defline::SubjectWord::Unknown`, `report::tabular::SubjectLookup::Failed`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum TabularSubjectId {
    Word(Vec<u8>),
    General(crate::seqdb::Oid),
    Unknown(crate::seqdb::Oid),
}

impl TabularSubjectId {
    /// `sseqid`, and `sallseqid` of the one definition line: the word (a local id prints its
    /// label), or `gnl|BL_ORD_ID|n`; empty where the lookup fails (`Unknown`).
    pub(crate) fn id(&self) -> std::borrow::Cow<'_, [u8]> {
        match self {
            TabularSubjectId::Word(word) => std::borrow::Cow::Borrowed(word),
            TabularSubjectId::General(oid) => std::borrow::Cow::Owned(oid.fasta_id().into_bytes()),
            TabularSubjectId::Unknown(_) => std::borrow::Cow::Borrowed(&[]),
        }
    }

    /// `sacc` and `saccver`, and `sallacc`: the word, or the general id's content label
    /// `BL_ORD_ID:n`; empty where the lookup fails (`Unknown`, which the programs' id lists hold
    /// as `None`, `common::tabular::SubjectId`).
    pub(crate) fn accession(&self) -> std::borrow::Cow<'_, [u8]> {
        match self {
            TabularSubjectId::Word(word) => std::borrow::Cow::Borrowed(word),
            TabularSubjectId::General(oid) => {
                std::borrow::Cow::Owned(oid.content_label().into_bytes())
            }
            TabularSubjectId::Unknown(_) => std::borrow::Cow::Borrowed(&[]),
        }
    }
}

/// The failed lookup (`report::tabular::FailedLookup`) of the database sequence `oid` without
/// parsed Seq-ids: its Seq-align id `gnl|BL_ORD_ID|n` (`sallseqid`) and `BL_ORD_ID:n` (`sallacc`),
/// whose definition line also fails in `x_SetSubjectIds`.
pub(crate) fn failed_general_lookup(oid: crate::seqdb::Oid) -> [Vec<u8>; 2] {
    [
        oid.fasta_id().into_bytes(),
        oid.content_label().into_bytes(),
    ]
}

/// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of the database
/// sequence `oid` with the stored title `title` and the source `source` of its Bioseq, from
/// which NCBI makes an empty title and rewrites the end of a protein title ending in `]`
/// (`report::defline::database_subject_id` with `stored_facts`; `protein` for a protein
/// database; `HeaderReader::stored_source`). A title whose id LOSAT cannot write as NCBI, and a
/// source that it refuses (`SourceTaxname::Refused`), are LOSAT's refusals (`program` names the
/// program).
pub(crate) fn tabular_subject_id(
    title: &[u8],
    source: &crate::report::defline::SourceTaxname,
    oid: crate::seqdb::Oid,
    protein: bool,
    program: &str,
) -> Result<TabularSubjectId> {
    let facts = stored_facts(source, source_name(source, || oid, program)?);
    match crate::report::defline::database_subject_id(title, &facts, protein) {
        Ok(SubjectWord::Word(word)) => Ok(TabularSubjectId::Word(word)),
        Ok(SubjectWord::General) => Ok(TabularSubjectId::General(oid)),
        Ok(SubjectWord::Unknown) => Ok(TabularSubjectId::Unknown(oid)),
        Err(problem) => bail!(
            "database sequence {oid} has a title {problem}, which is not supported by LOSAT's {program}"
        ),
    }
}

/// The facts of the Bioseq of a database sequence without parsed Seq-ids whose source is
/// `source`, named `name` (`SourceTaxname::name`): `BioseqFacts::of_source`, with a source
/// without a name too (`SourceTaxname::Nameless`), for which `x_AdjustProteinTitleSuffix` runs as
/// well (`report::defline::adjust_protein_title_suffix`). One form for the tabular ids
/// (`tabular_subject_id`, `tabular_seqid`).
fn stored_facts(
    source: &crate::report::defline::SourceTaxname,
    name: Option<std::sync::Arc<[u8]>>,
) -> crate::report::defline::BioseqFacts {
    crate::report::defline::BioseqFacts {
        has_source: source.has_source(),
        ..crate::report::defline::BioseqFacts::of_source(name)
    }
}

/// `sseqid` of the database sequence `oid` without parsed Seq-ids with the stored title `title`
/// and the source `source` of its Bioseq, whose `saccver` is `accession` (`tabular_subject_ids`):
/// the accession, unless it is the general id's (`BL_ORD_ID:n`, `TabularSubjectId::General`),
/// where the title decides whether it is (an id of that form asks it again; a title whose first
/// word is `BL_ORD_ID:n` prints that word).
pub(crate) fn tabular_seqid<'a>(
    accession: &'a [u8],
    oid: crate::seqdb::Oid,
    title: &[u8],
    source: &crate::report::defline::SourceTaxname,
    protein: bool,
) -> std::borrow::Cow<'a, [u8]> {
    if accession.starts_with(b"BL_ORD_ID:")
        && accession == oid.content_label().as_bytes()
        && source.name().is_ok_and(|name| {
            crate::report::defline::database_subject_id(title, &stored_facts(source, name), protein)
                == Ok(SubjectWord::General)
        })
    {
        return std::borrow::Cow::Owned(oid.fasta_id().into_bytes());
    }
    std::borrow::Cow::Borrowed(accession)
}

/// LOSAT's rejection of the outfmt 0 title `title` of a protein database sequence without parsed
/// Seq-ids (`oid` gives its OID, which the rejection names) that it does not write as NCBI: a
/// title of punctuation that NCBI's `x_CleanAndCompress` reads past the end of when it makes the
/// title of the description table (NCBI BLAST+ 2.17.0 crashes there), as for the records of a
/// subject set (`check_shown_subject_titles`). Approved exception 2 of PD-LOSAT-NCBI-DEFECTS
/// covers the nucleotide subjects of BLASTN, so BLASTN does not use it.
///
/// The name of the Bioseq's source does not change the answer, so none is asked for: the title
/// is made from the source only where the stored one is empty, and the Bioseq of a sequence
/// without parsed Seq-ids (its one id `gnl|BL_ORD_ID|n` sets no flag of `BioseqFacts`) then gets
/// `unnamed protein product` and ` [name]`, which ends in `]` or a letter whatever the name is
/// (`report::defline::source_title`). A title made of the name alone (a RefSeq record's, which
/// can read past its end) has parsed Seq-ids: `check_shown_parsed_subject`, which every reader of
/// such a sequence asks (the block search's check of the titles, `x_titles_allow_blocks`, takes
/// none: its definition lines are not `SeqDb::deflines`). Pinned by the unit test
/// `a_stored_protein_title_from_any_source_ends_in_its_own_text`.
pub(crate) fn check_shown_protein_title(
    title: &[u8],
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<()> {
    let reads_past_end = crate::report::defline::generate_bioseq_defline(
        title,
        true,
        true,
        &crate::report::defline::BioseqFacts::default(),
    )
    .is_ok_and(|defline| defline.reads_past_end);
    if reads_past_end {
        bail!(
            "database sequence {} has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0); this is not supported by LOSAT's {program}",
            oid()
        );
    }
    Ok(())
}

/// Byte 0 of NCBIstdaa: the gap character (`-`), and the sentinel byte between the
/// sequences of a protein volume's sequence file.
const NCBISTDAA_GAP: u8 = 0;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593
// A stored gap character (0, which makeblastdb writes for `-` in a protein FASTA) is a
// residue that scores `BLAST_SCORE_MIN` against every letter (row and column 0), as LOSAT's
// kernels score it: an ungapped extension stops at it and a gapped alignment goes round it
// with a gap. A byte that is no NCBIstdaa residue (28 or more) indexes NCBI's matrix past
// its rows and is rejected.
/// LOSAT's refusal of a protein database sequence (`oid` gives its OID, which the refusal names;
/// `with_sentinels`: its residues with the byte before and the byte after them in the sequence
/// file) when
/// `program` cannot search it as NCBI does: a sequence not between two 0 bytes, or one
/// with a byte that is no NCBIstdaa residue.
pub(crate) fn check_protein_sequence(
    with_sentinels: &[u8],
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<()> {
    let (first, last) = (with_sentinels[0], with_sentinels[with_sentinels.len() - 1]);
    if first != NCBISTDAA_GAP || last != NCBISTDAA_GAP {
        bail!(
            "database sequence {} is not between two NUL bytes in its sequence file (the volume is damaged), which is not supported by LOSAT's {program}",
            oid()
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
            "database sequence {} holds a byte that is no NCBIstdaa residue at position {}, which is not supported by LOSAT's {program}",
            oid(),
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
            check_protein_sequence(with_sentinels, || db.oid(ordinal), program)?;
            each(ordinal, with_sentinels.len() - 2);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

#[cfg(test)]
#[path = "database_tests.rs"]
mod tests;

/// The taxonomy database of a run (`CTaxDBFileInfo`, a static of NCBI's `GetTaxNames`) and the
/// lookups the volumes have made (`CSeqDBVol::x_GetTaxonomy` caches them per volume), in the
/// order in which LOSAT records them ([`TaxdbRun::load`]).
///
/// NCBI's cache of a volume has 256 slots, one per `taxid & 255` (`CSeqDBVol::m_TaxCache`,
/// below): a lookup asks the taxonomy database where its slot holds nothing or another taxid,
/// and the taxid then takes the slot. Of the texts a lookup can write, those of opening come once per
/// process and only the offset error of the last entry (`TaxDb::offset_error_taxid`, `L`) comes
/// again at each lookup that misses, so the record keeps what decides the bytes
/// ([`VolumeLookups`]): per volume what `L`'s slot holds, which only the lookups of `L` and of
/// another positive taxid of its slot change; the taxonomy database is asked where `L` misses, and
/// at every other lookup only where a lookup can fault ([`TaxdbRun::volume_load`]). A run whose
/// volumes hold no such other taxid looks `L` up once per volume, in any order.
///
/// On a database with a volume that holds `L` and another taxid of its slot
/// ([`TaxdbRun::cache_collision`]), a batch's prefetch follows the slot in NCBI's order of its
/// loads, by the keys of their Seq-id handles; loads of one key, which NCBI orders by the
/// addresses of the handles' infos, are refused only where their order changes the texts or
/// leaves the slot to a later load's texts ([`prefetch_slot_plan`]). The reports' reads where a
/// prefetch is off ([`subject_read_texts`], [`prefetch_messages`]) and, in the XML2 and JSON
/// formats, a report with a subject that has a positive taxid ([`read_reported_taxonomy`]) are
/// refused on such a database where NCBI first fetches a subject. psiblast's PSSM
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
    lazy: std::sync::OnceLock<std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>>,
    volume_lookups: std::sync::Mutex<VolumeLookups>,
    found: std::sync::OnceLock<bool>,
    /// What the lookups of the run can write ([`TaxdbRun::prefetch_texts`]).
    texts: std::sync::OnceLock<PrefetchTexts>,
    /// The refusal of the reads of a run with a volume where NCBI's cache can make the fetches
    /// write the offset error more than once, as a reader's refusal ([`TaxdbRun::cache_collision`]).
    collision: std::sync::OnceLock<Option<(std::path::PathBuf, String)>>,
    /// The texts of opening the database when the rows of a tabular report opened it, kept for
    /// those rows ([`TaxdbRun::open_for_rows`], [`tabular_names`]).
    rows_open: std::sync::OnceLock<std::sync::Arc<std::sync::Mutex<Vec<String>>>>,
    /// The sequences with parsed Seq-ids, by ordinal id, that NCBI's BLAST database data loader
    /// has loaded in the run before a batch whose prefetch stops ([`TaxdbRun::note_loaded`]).
    loaded: std::sync::Mutex<std::collections::HashSet<u32>>,
    /// Whether a prefetch of the run can stop, so that the run keeps `loaded`
    /// ([`TaxdbRun::prefetch_can_stop`]), decided once.
    can_stop: std::sync::OnceLock<bool>,
    /// Whether every load of a sequence with a positive taxid throws
    /// ([`TaxdbRun::load_throws`]), decided once.
    load_throws: std::sync::OnceLock<bool>,
    /// The texts of the loads that the model of the reads placed where NCBI makes no load
    /// ([`TaxdbRun::load_texts`], `report::query_warnings::LoadTexts`).
    load_texts: std::sync::Arc<crate::report::query_warnings::LoadTexts>,
}

/// What the lookups of a run's taxonomy database can write ([`TaxdbRun::prefetch_texts`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct PrefetchTexts {
    /// The positive taxid whose lookup writes the offset error (`TaxDb::offset_error_taxid`).
    offset_error_taxid: Option<i32>,
    /// Whether a lookup can fault: NCBI maps an empty `taxdb.btd` as a null pointer
    /// (`TaxDb::maps_null`).
    maps_null: bool,
    /// Whether opening the database writes texts (`TaxDb::open_messages`).
    opening: bool,
}

impl PrefetchTexts {
    /// What the lookups of the opened taxonomy database `taxdb` can write: the one producer for
    /// the search's loads ([`TaxdbRun::prefetch_texts`]) and the readers' ([`TaxdbRun::bioseq_loads`]).
    /// A lookup of the last entry's taxid that fails writes no offset error, and the null mapping
    /// and the texts of opening are the file's whatever that lookup gives: with an empty
    /// `taxdb.btd` its failure is NCBI's fault reading a record (an index whose search for that
    /// taxid ends on another entry), which the run's loads meet at their first lookup that reads
    /// one ([`TaxdbRun::volume_load`]).
    fn of(taxdb: &crate::seqdb::taxdb::TaxDb) -> Self {
        PrefetchTexts {
            // `x_GetTaxonomy` looks a positive taxid up only.
            offset_error_taxid: taxdb
                .offset_error_taxid()
                .ok()
                .flatten()
                .filter(|&taxid| taxid > 0),
            maps_null: taxdb.maps_null(),
            opening: !taxdb.open_messages().is_empty(),
        }
    }
}

/// The size 256 of `CSeqDBVol::m_TaxCache`: the slots of a volume's taxonomy cache.
const TAX_CACHE_SLOTS: i32 = 256;

/// The slot of NCBI's taxonomy cache of a volume that `taxid` takes (`CSeqDBIntCache::Lookup`,
/// `key` ANDed with one less than the size of `m_Slots`).
fn tax_cache_slot(taxid: i32) -> i32 {
    taxid & (TAX_CACHE_SLOTS - 1)
}

/// The volume lookups of a run ([`TaxdbRun`]), in the order in which LOSAT records them
/// ([`TaxdbRun::load`]): only what decides the texts of the offset error, which only a taxonomy
/// database with that error records. Nothing grows with the taxids looked up.
#[derive(Default)]
struct VolumeLookups {
    /// By volume, what the slot of the offset error's taxid holds ([`OffsetSlot`]); only where
    /// the taxonomy database has that error.
    offset_slots: rustc_hash::FxHashMap<usize, OffsetSlot>,
    /// The sequences, by ordinal id, whose loads looked up a taxid of that slot: a later read of
    /// one may find it in NCBI's scope and load nothing, or load it again.
    slot_loads: rustc_hash::FxHashSet<u32>,
    /// A run of a program without NCBI's object manager ([`TaxdbRun::sequential`]): every read
    /// of a sequence is a load (`CSeqDB::GetBioseq`), a read again too.
    sequential: bool,
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
    /// Where `holds` is `Either` because of the order of a prefetch's loads
    /// ([`VolumeLookups::unordered_slot`]): two of them, for LOSAT's refusal.
    unordered: Option<SlotOrder>,
}

/// Two loads of one batch's prefetch in the same volume whose order NCBI takes from the addresses
/// of their Seq-id handles ([`load_before_stop`]), which LOSAT does not reproduce, named by its
/// refusal: `first` a load that looks up the offset error's taxid where the group has one (else
/// the group's first load), or the load read again; `second` another load, one that looks up
/// the other taxid `taxid` of the slot where the group has one (else `taxid` is the offset
/// error's).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SlotOrder {
    first: u32,
    second: u32,
    taxid: i32,
}

/// A load of the sequence `ordinal` whose texts of the offset error NCBI's scope decides
/// ([`VolumeLookups::offset_slot`]): `again`, a read of a sequence loaded before, where keeping
/// it and loading it again give other texts; else a load after such a read. `after`: the load of
/// the other taxid of the slot before it.
///
/// `order`: the cause is the order of a prefetch's loads instead ([`SlotOrder`]): the load
/// `ordinal` is one of them where their order decides how often its offset error is written
/// (`ordinal == order.first`), else a later load that meets what that order left in the slot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SlotUncertain {
    ordinal: u32,
    again: bool,
    after: Option<SlotFetch>,
    order: Option<SlotOrder>,
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1142-1182 (the
    // constructor's initializer that gives `m_TaxCache` 256 slots) and 1142-1182
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
    /// they differ). `Err` names the load that is not certain ([`SlotUncertain`]). In a
    /// sequential run every read is a first load: the slot is followed, never `Err`.
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
        let again = !self.sequential && !self.slot_loads.insert(ordinal);
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
            order: before.unordered,
        };
        let walk = before.holds.walk(lookups, last).map_err(uncertain)?;
        if !again {
            slot.holds = walk.holds;
            // A walk with a lookup of the slot leaves `L` or another taxid there.
            slot.unordered = None;
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
            slot.unordered = None;
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

    /// The misses of the loads `loads` of one batch's prefetch that look a taxid of the slot of
    /// the offset error's taxid `last` up ([`SlotLoad`]), by the ordinal id of the reported
    /// subject, in NCBI's order of the prefetch's loads (cited above [`LoadBeforeStop`]): by the
    /// keys of their Seq-id handles, loads of one key (`order`) in the order of the addresses of
    /// the handles' infos, which LOSAT does not reproduce (all of a volume's loads where
    /// `unordered`: the batch has more than [`PREFETCH_ONE_MAP_IDS`] distinct Seq-align ids, or
    /// such a load's id is another subject's too). A load of its own key follows the slot
    /// ([`VolumeLookups::offset_slot`]); loads of one key are one step ([`Self::unordered_slot`]).
    /// `Err`: the volume and the load whose texts NCBI's order or scope decides.
    fn prefetch_slots(
        &mut self,
        loads: &[SlotLoad],
        unordered: impl Fn(usize) -> bool,
        last: i32,
    ) -> std::result::Result<rustc_hash::FxHashMap<u32, Vec<bool>>, (usize, SlotUncertain)> {
        let mut misses = rustc_hash::FxHashMap::default();
        let mut volumes: Vec<usize> = loads.iter().map(|load| load.volume).collect();
        volumes.sort_unstable();
        volumes.dedup();
        for volume in volumes {
            let all = unordered(volume);
            let mut members: Vec<&SlotLoad> =
                loads.iter().filter(|load| load.volume == volume).collect();
            // A sort that keeps the batch's order among loads of one key.
            members.sort_by_key(|load| if all { 0 } else { load.order });
            for group in members.chunk_by(|a, b| all || a.order == b.order) {
                if let [load] = group {
                    let own = self
                        .offset_slot(volume, load.blob, &load.lookups, last)
                        .map_err(|cause| (volume, cause))?;
                    misses.insert(load.ordinal, own);
                } else {
                    for (ordinal, own) in self
                        .unordered_slot(volume, group, last)
                        .map_err(|cause| (volume, cause))?
                    {
                        misses.insert(ordinal, own);
                    }
                }
            }
        }
        Ok(misses)
    }

    // NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:201-208 (cited at
    // `load_before_stop`: handles of one packed number by their infos' addresses)
    /// The misses of the loads `group` of one prefetch in the volume `volume`, all of one key of
    /// NCBI's map of the batch's Seq-id handles, which NCBI makes in an order LOSAT does not
    /// reproduce: each load looks the slot of `last` up as its lookups go (`L` for `last`, `O`
    /// for another taxid of the slot), and a lookup of `L` misses after an `O` or where the slot
    /// held another taxid before the group; only the boundaries between the loads depend on their
    /// order ([`unordered_boundary_misses`]). Where every order gives the same number of misses
    /// (from each state the slot may hold before), the texts are the same in every order: each load
    /// gets the misses of the batch's order (the prefetch writes the texts of its loads as one, in
    /// any order of the loads: the texts of opening first, then offset errors of one text), and the
    /// slot holds what the last load leaves, or `Either` where the loads leave other taxids.
    /// `Err`: a load loaded before in the run (NCBI's scope may keep it) or an order that changes
    /// the number of misses.
    fn unordered_slot(
        &mut self,
        volume: usize,
        group: &[&SlotLoad],
        last: i32,
    ) -> std::result::Result<Vec<(u32, Vec<bool>)>, SlotUncertain> {
        let in_slot = |taxid: i32| tax_cache_slot(taxid) == tax_cache_slot(last);
        let before = self.offset_slots.get(&volume).copied().unwrap_or_default();
        // A load with `L` and one with another taxid of the slot (or a second with `L`), named by
        // LOSAT's refusal.
        let first = group
            .iter()
            .find(|load| load.lookups.contains(&last))
            .unwrap_or(&group[0]);
        let (second, taxid) = group
            .iter()
            .filter(|load| load.ordinal != first.ordinal)
            .find_map(|load| {
                load.lookups
                    .iter()
                    .find(|&&taxid| taxid != last && in_slot(taxid))
                    .map(|&taxid| (load, taxid))
            })
            .unwrap_or((
                group
                    .iter()
                    .find(|load| load.ordinal != first.ordinal)
                    .expect("a group of two loads or more"),
                last,
            ));
        let order = SlotOrder {
            first: first.blob,
            second: second.blob,
            taxid,
        };
        // A load of a sequence loaded before may be kept by NCBI's data loader (`offset_slot`), at
        // a place LOSAT does not know either: the refusal names that load.
        let mut blobs = std::collections::HashSet::new();
        for load in group {
            if !self.sequential
                && (self.slot_loads.contains(&load.blob) || !blobs.insert(load.blob))
            {
                let other = group
                    .iter()
                    .find(|other| other.blob != load.blob)
                    .map_or(load.blob, |other| other.blob);
                return Err(SlotUncertain {
                    ordinal: load.blob,
                    again: true,
                    after: before.other,
                    order: Some(SlotOrder {
                        first: load.blob,
                        second: other,
                        taxid,
                    }),
                });
            }
        }
        // Each load's lookups of the slot, a run of one taxid as one: its first and last, and the
        // misses inside it.
        let mut counts = UnorderedLoads::default();
        let mut ends = Vec::with_capacity(group.len());
        for load in group {
            let marks: Vec<bool> = load
                .lookups
                .iter()
                .filter(|&&taxid| in_slot(taxid))
                .map(|&taxid| taxid == last)
                .collect();
            let (starts_last, ends_last) = (marks[0], marks[marks.len() - 1]);
            counts.add(starts_last, ends_last);
            ends.push(ends_last);
        }
        let starts: &[bool] = match before.holds {
            SlotHolds::Other => &[false],
            SlotHolds::Last => &[true],
            SlotHolds::Either => &[false, true],
        };
        let ranges: Vec<(usize, usize)> = starts
            .iter()
            .map(|&held_last| unordered_boundary_misses(held_last, counts))
            .collect();
        if ranges.iter().any(|&(fewest, most)| fewest != most) {
            // The group's own order decides.
            return Err(SlotUncertain {
                ordinal: first.blob,
                again: false,
                after: before.other,
                order: Some(order),
            });
        }
        if ranges.iter().any(|&range| range != ranges[0]) {
            // Any order gives the same misses from each state, which an earlier load left open:
            // that load is the cause (`offset_slot`'s refusal of a load after it).
            return Err(SlotUncertain {
                ordinal: first.blob,
                again: false,
                after: before.other,
                order: before.unordered,
            });
        }
        // The batch's order gives the same number of misses as NCBI's (from `Other` where the slot
        // may hold either: the same number from each).
        let mut holds = match before.holds {
            SlotHolds::Last => SlotHolds::Last,
            _ => SlotHolds::Other,
        };
        let mut other = before.other;
        let mut misses = Vec::with_capacity(group.len());
        for load in group {
            let walk = holds
                .walk(&load.lookups, last)
                .expect("a walk from a definite state");
            holds = walk.holds;
            if let Some(taxid) = walk.other {
                other = Some(SlotFetch {
                    ordinal: load.blob,
                    taxid,
                    again: false,
                });
            }
            self.slot_loads.insert(load.blob);
            misses.push((load.ordinal, walk.misses));
        }
        let slot = self.offset_slots.entry(volume).or_default();
        slot.other = other;
        if ends.iter().all(|&end| end == ends[0]) {
            slot.holds = holds;
            slot.unordered = None;
        } else {
            slot.holds = SlotHolds::Either;
            slot.unordered = Some(order);
        }
        Ok(misses)
    }
}

/// One load of a batch's prefetch that looks a taxid of the slot of the offset error's taxid up
/// ([`VolumeLookups::prefetch_slots`]).
struct SlotLoad {
    /// The reported subject, by ordinal id.
    ordinal: u32,
    /// The sequence loaded for it, by ordinal id ([`prefetch_reads`]).
    blob: u32,
    /// The blob's volume.
    volume: usize,
    /// The blob's lookups ([`blob_taxids`]).
    lookups: Vec<i32>,
    /// The key of the subject's Seq-id handle in NCBI's map of the batch's loads
    /// ([`prefetch_load_key`]).
    order: u64,
}

/// The loads of one step of [`VolumeLookups::unordered_slot`], by the first and the last of their
/// lookups of the slot (`L`: the offset error's taxid; `O`: another taxid of the slot).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct UnorderedLoads {
    /// First and last `L`.
    ll: usize,
    /// First `L`, last `O`.
    lo: usize,
    /// First `O`, last `L`.
    ol: usize,
    /// First and last `O`.
    oo: usize,
}

impl UnorderedLoads {
    fn add(&mut self, starts_last: bool, ends_last: bool) {
        match (starts_last, ends_last) {
            (true, true) => self.ll += 1,
            (true, false) => self.lo += 1,
            (false, true) => self.ol += 1,
            (false, false) => self.oo += 1,
        }
    }
}

/// The fewest and the most misses at the boundaries of the loads `loads` over every order of
/// them, after a slot that holds `L` (`held_last`) or another taxid: a load whose first lookup is
/// `L` misses where the load before it (or the slot before the first) left another taxid. The
/// misses inside each load are the same in every order. Most: each `O` before an `L` gives one,
/// a chain of `LO` loads after one `O` (an `OO` load or the slot) all, and each `O` that heads a
/// chain one more for an `LL` load at its end; without such an `O` the first `LO` heads one
/// unmissed. Fewest: each `L` before (an `OL` load or the slot) protects a chain of `LL` loads and
/// one `LO` load at its end; without such an `L` an `LL` load takes one miss and protects as
/// much, and `LO` loads without a chain miss each (`tests::unordered_misses_are_every_order`
/// counts every order of up to seven loads).
fn unordered_boundary_misses(held_last: bool, loads: UnorderedLoads) -> (usize, usize) {
    let heads_o = loads.oo + usize::from(!held_last);
    let most = if heads_o > 0 {
        loads.lo + loads.ll.min(heads_o)
    } else if loads.lo > 0 {
        loads.lo - 1 + loads.ll.min(1)
    } else {
        0
    };
    let heads_l = loads.ol + usize::from(held_last);
    let fewest = if heads_l > 0 {
        loads.lo.saturating_sub(heads_l)
    } else if loads.ll > 0 {
        loads.lo.max(1)
    } else {
        loads.lo
    };
    (fewest, most)
}

impl TaxdbRun {
    /// Back to the state of a run that has looked nothing up, for a search of the batches before
    /// again (BLASTP's `HeaderReadError`: a later batch failed, and the batches before it are
    /// searched again alone): their prefetch and reports write the texts of opening the database
    /// and of each volume's first lookups again, as the first pass made them, before the error.
    /// What does not depend on the lookups (`found`, the texts the database can write, the
    /// refusal of its cache) is kept.
    /// The run's handle stays and writes the texts of opening again at its next lookup
    /// (`LazyTaxDb::restart`), as a new handle would.
    pub(crate) fn restart(&self) {
        if let Some(lazy) = self.lazy.get() {
            lazy.restart();
        }
        let mut record = self.volume_lookups.lock().expect("the volume lookups");
        *record = VolumeLookups {
            sequential: record.sequential,
            ..VolumeLookups::default()
        };
        drop(record);
        if let Some(rows_open) = self.rows_open.get() {
            rows_open.lock().expect("the rows' texts").clear();
        }
        self.load_texts.clear();
    }

    /// Whether a prefetch of the run can stop at a lookup that fails on a file with NCBI's
    /// exception ([`batch_prefetch`]; every producer of `FoundInstead::Failed` with NCBI's message,
    /// `HeaderReader::found_instead`), the one gate of what the run keeps for such a stop (its
    /// record of loads, [`TaxdbRun::note_loaded`]; BLASTX's scopes of the batches before in outfmt
    /// 5): the lookup of `db` cannot be opened, or one of its LMDB files has no accession table
    /// (`SeqidLookup::lacks_an_accession_table`), or an ISAM of a version 4 volume may fail to open
    /// (`SeqDb::isams_may_fail_to_open`: the lookup's own open of a string or PIG ISAM, which
    /// fails with NCBI's `Error: Could not open input file` where both files are regular files
    /// and one does not open; review 3 of line DBID loadord, M1, and review 4, H1). Asked once, at
    /// the first record of the run's loads (`note_loaded`, before it makes its subjects, review
    /// 4, M1) or of psiblast's PSSM's reads: the files are opened, never read, and a database
    /// whose lookups cannot stop keeps nothing for the rest of the run.
    pub(crate) fn prefetch_can_stop(&self, db: &SeqDb) -> bool {
        *self.can_stop.get_or_init(|| {
            db.seqid_lookup()
                .map_or(true, |lookup| lookup.lacks_an_accession_table())
                || db.isams_may_fail_to_open()
        })
    }

    /// Records that NCBI's BLAST database data loader has loaded `subjects` in the run (a query
    /// batch's, or a psiblast round's; [`PrefetchSubject`], in any order), for the stop of a later
    /// batch's prefetch ([`stop_lookup_unsure`]): the loader keeps the Seq-ids of every sequence
    /// that it loaded (`RegisterIds`). A batch's loads are every subject of its Seq-aligns where
    /// its prefetch is on, else the subjects that the reports NCBI writes (after their cut to
    /// `-max_target_seqs`; not an outfmt 5 report that its XML stream cut before its hits) and
    /// psiblast's PSSM read ([`batch_loads`], [`LoadedBefore::subjects`]; review 3 of line DBID
    /// loadord, H2; review 6, H1 and H2). A caller that searches the queries of
    /// a run in several calls, or writes each batch's reports before it searches the next
    /// (BLASTX), records each batch that does not stop; the callers that hold a call's batches give
    /// the loads of the batches before to the stop itself ([`batch_prefetch`]'s `before`). Only the
    /// sequences with parsed Seq-ids are kept (one without them holds no id that a lookup can fail
    /// on), and only in a run that can stop ([`TaxdbRun::prefetch_can_stop`], asked first): the
    /// subjects are made (`subjects`, called once) only there, so a run that cannot stop builds
    /// and walks nothing (review 4 of line DBID loadord, M1). The callers ask it only of a
    /// database with parsed Seq-ids: nothing is opened for any other.
    pub(crate) fn note_loaded<'s, I>(&self, db: &SeqDb, subjects: impl FnOnce() -> I)
    where
        I: IntoIterator<Item = PrefetchSubject<'s>>,
    {
        if self.prefetch_can_stop(db) {
            self.keep_loaded(subjects());
        }
    }

    /// The record of [`TaxdbRun::note_loaded`] without its gate: at a stop
    /// ([`batch_prefetch`], [`held_batch_prefetch_stop`]), which only a run that can stop reaches,
    /// the subjects that the call's batches before loaded.
    fn keep_loaded<'s>(&self, subjects: impl IntoIterator<Item = PrefetchSubject<'s>>) {
        let mut parsed = subjects
            .into_iter()
            .filter(|subject| subject.parsed.is_some())
            .map(|subject| subject.ordinal)
            .peekable();
        if parsed.peek().is_none() {
            return;
        }
        self.loaded
            .lock()
            .expect("the loaded sequences")
            .extend(parsed);
    }

    /// The sequences that [`TaxdbRun::note_loaded`] recorded, in the order of their ordinal ids.
    fn loaded(&self) -> Vec<u32> {
        let mut loaded: Vec<u32> = self
            .loaded
            .lock()
            .expect("the loaded sequences")
            .iter()
            .copied()
            .collect();
        loaded.sort_unstable();
        loaded
    }

    /// [`SeqDb::taxdb_found`] of the run's database, asked once: the search path is not
    /// looked through again for each batch.
    pub(crate) fn found(&self, db: &SeqDb) -> bool {
        *self.found.get_or_init(|| db.taxdb_found())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1155 (cited
    // at `taxdb_error`)
    /// Whether a file of the run's taxonomy database does not open where NCBI's `CMemoryFile`
    /// throws (`SeqDbError::TaxdbThrows`, `seqdb::taxdb`'s `throws_at_every_lookup`): every load
    /// of a sequence with a positive taxid throws then, so psiblast's PSSM reads no letters of
    /// such a subject (`Failed to retrieve sequence`, `ListedTitles::pssm_reads` of BLASTP), its
    /// loads write no text ([`pssm_load_texts`]) and the organism report leaves its HSPs out
    /// (`report::tax_report::SubjectTable::with_loads_that_throw`). Decided once, without
    /// `taxdb.bti` on the search path as false; the files are opened apart from the run's handle,
    /// whose texts of opening stay for its first lookup.
    pub(crate) fn load_throws(&self, db: &SeqDb) -> bool {
        *self.load_throws.get_or_init(|| {
            self.found(db)
                && db.search_path().is_some_and(|search_path| {
                    matches!(
                        crate::seqdb::taxdb::TaxDb::open(search_path),
                        Err(SeqDbError::TaxdbThrows { .. })
                    )
                })
        })
    }

    /// [`PrefetchTexts`] of the run, read once ([`PrefetchTexts::of`]): none without `taxdb.bti`
    /// on the search path.
    fn prefetch_texts(&self, db: &SeqDb) -> PrefetchTexts {
        *self.texts.get_or_init(|| {
            if !self.found(db) {
                return PrefetchTexts::default();
            }
            let Some(search_path) = db.search_path() else {
                return PrefetchTexts::default();
            };
            // A database that does not open fails at the run's first lookup, through its handle.
            let Ok(taxdb) = crate::seqdb::taxdb::TaxDb::open(search_path) else {
                return PrefetchTexts::default();
            };
            PrefetchTexts::of(&taxdb)
        })
    }

    // NCBI's 256-slot cache of a volume (cited at `SlotHolds::walk`): a fetch of a sequence with
    // another taxid of the slot of `L` clears it, and the next fetch of `L` asks `GetTaxNames`
    // again, which writes the offset error again. How often it repeats in a prefetch follows the
    // order in which the object manager loads the batch's subjects, the keys of their Seq-id
    // handles and, for one key (every unpacked id), the addresses of their infos (review 1 of item
    // T: one query on a volume of the taxids 99990, 150, 99990 writes the message twice, two
    // queries once; fix 1: the same query once with a version 4 database or an absolute
    // `BLASTDB`): `prefetch_slot_plan` follows the keys and refuses where the addresses decide.
    /// Whether a volume of the database holds the taxid `L` of the offset error and another
    /// positive taxid of its slot, where NCBI's fetches can write the offset error more than once
    /// in a volume: the taxonomy database has the offset error of a positive taxid `L`
    /// ([`TaxdbRun::prefetch_texts`]), and a volume of the database holds `L` and another
    /// positive taxid of the same cache slot (`T & 255 == L & 255`), in any definition line of any
    /// sequence, whether or not a search reports them. One pass over the headers
    /// ([`SeqDb::visit_defline_set_taxids`] in ranges of [`crate::seqdb::HEADER_RANGE_SEQS`]), made
    /// only for such a damaged taxonomy database and once for the run; a header whose taxids do
    /// not decode is left to the lookups of the batches, which fail on it where it is reported.
    /// `(path of taxdb.bti, reason)` of the refusal of the fetches whose order LOSAT does not
    /// follow there (the reports' reads, [`load_read_texts`], and [`read_reported_taxonomy`]),
    /// made as a reader's refusal ([`TaxonomyRefusal`]) where NCBI fetches
    /// ([`TaxdbRun::refuse_cache_collision`]); a batch's prefetch follows NCBI's order instead
    /// ([`prefetch_slot_plan`]).
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

    /// Whether a fetch of the run may be refused because of NCBI's taxonomy cache: a volume holds
    /// the offset error's taxid and another taxid of its slot ([`TaxdbRun::cache_collision`]).
    /// The reports' reads and the XML2 and JSON lookups of such a run are refused, and a batch's
    /// prefetch where NCBI's order of its loads decides the texts ([`prefetch_slot_plan`]), which
    /// can be a later batch: BLASTX holds its XML document (outfmt 5, and 13 to 16) then, which
    /// the refusal drops; TBLASTX turns its report stream off and BLASTP takes the search of all
    /// queries, which make every batch's lookups before the first byte.
    pub(crate) fn cache_may_refuse(&self, db: &SeqDb) -> bool {
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
        let texts = match self.taxdb(search_path).open() {
            Ok(texts) => texts,
            // A file that does not open throws at each lookup and writes nothing: the rows decide
            // (the organism report leaves the HSP out, `report::tax_report::TaxLookups::look_up`;
            // a name field is refused, `tabular_names`).
            Err(SeqDbError::TaxdbThrows { .. }) => return Ok(()),
            Err(error) => return Err(taxdb_error(error, program)),
        };
        self.rows_open()
            .lock()
            .expect("the texts of opening")
            .extend(texts);
        Ok(())
    }

    /// The texts of opening kept for the rows ([`TaxdbRun::open_for_rows`]), made at the first use.
    fn rows_open(&self) -> &std::sync::Arc<std::sync::Mutex<Vec<String>>> {
        self.rows_open.get_or_init(Default::default)
    }

    /// The run's taxonomy database, a handle that opens it at its first name asked for.
    fn taxdb(
        &self,
        search_path: &crate::seqdb::SearchPath,
    ) -> std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb> {
        std::sync::Arc::clone(self.lazy_taxdb(search_path))
    }

    /// [`TaxdbRun::taxdb`], borrowed.
    fn lazy_taxdb(
        &self,
        search_path: &crate::seqdb::SearchPath,
    ) -> &std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb> {
        self.lazy.get_or_init(|| {
            std::sync::Arc::new(crate::seqdb::taxdb::LazyTaxDb::new(search_path.clone()))
        })
    }

    /// NCBI's lookups of the volume's taxonomy cache where it loads the sequence `ordinal` of `db`:
    /// the positive taxids `lookups` in the order of `CSeqDBVol::x_GetTaxonomy`, a taxid of several
    /// definition lines once for each ([`blob_taxids`]). The texts of opening come with the run's
    /// first lookup; the taxonomy database is asked where a lookup of the offset error's taxid
    /// misses its slot of the cache ([`VolumeLookups::offset_slot`]) and, where `taxdb.btd` is mapped
    /// as a null pointer, at every lookup of another taxid ([`TaxdbRun::volume_load`]). Returns what
    /// each lookup that wrote a text wrote, by taxid, in the order of the lookups. A load whose texts NCBI's scope decides (a read again of a
    /// sequence loaded before, where its scope keeping it and loading it again write other texts,
    /// or a load after such a read) is LOSAT's refusal ([`TaxonomyRefusal`]), naming `program`.
    /// `planned`: the misses of a load of a batch's prefetch that NCBI's order of the loads gave
    /// ([`prefetch_slot_plan`]), recorded there.
    fn load(
        &self,
        db: &SeqDb,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        ordinal: u32,
        lookups: &[i32],
        planned: Option<&[bool]>,
        program: &str,
    ) -> Result<Vec<(i32, Vec<String>)>> {
        if lookups.is_empty() {
            return Ok(Vec::new());
        }
        let volume = db
            .volume_index(ordinal)
            .map_err(|error| database_error(error, program))?;
        let texts = self.prefetch_texts(db);
        let misses = match (texts.offset_error_taxid, planned) {
            (Some(_), Some(planned)) => Some(planned.to_vec()),
            (Some(last), None) => Some(self.slot_misses(volume, ordinal, lookups, last).map_err(
                |cause| taxdb_error(slot_uncertain_error(db, volume, last, cause), program),
            )?),
            (None, _) => None,
        };
        let open = taxdb.open().map_err(|error| taxdb_error(error, program))?;
        let mut answers = Vec::new();
        self.volume_load(
            taxdb,
            texts,
            misses.as_deref(),
            lookups.iter().copied(),
            open,
            &mut |taxid, messages| answers.push((taxid, messages)),
        )
        .map_err(|error| taxdb_error(error, program))?;
        Ok(answers)
    }

    /// Which of the lookups `lookups` of a load of the sequence `ordinal` in the volume `volume`
    /// ask the taxonomy database again for the offset error's taxid `last`
    /// ([`VolumeLookups::offset_slot`]).
    fn slot_misses(
        &self,
        volume: usize,
        ordinal: u32,
        lookups: &[i32],
        last: i32,
    ) -> std::result::Result<Vec<bool>, SlotUncertain> {
        self.volume_lookups
            .lock()
            .expect("the volume lookups")
            .offset_slot(volume, ordinal, lookups, last)
    }

    // NCBI asks `GetTaxNames` where the volume's 256-slot cache misses (cited at
    // `SlotHolds::walk`). Of what such an ask does, only three things reach the output of a load,
    // whose names nothing writes (`GetBioseq`'s sources; a title made from one asks the names
    // itself): the texts of opening the database at the process's first ask, the offset error at
    // each ask of its taxid `L`, and a fault where `taxdb.btd` is mapped as a null pointer. A
    // record read past a mapping that is not null is not a fault for certain and writes nothing
    // (`TaxDb::lookup`). So the taxonomy database is asked where `L` misses its slot
    // (`VolumeLookups::offset_slot`) and, with a null mapping, at every lookup of another taxid (the
    // first one that reads a record is NCBI's fault); the cache's other slots decide nothing.
    /// The asks of a load ([`TaxdbRun::load`]) whose lookups are `lookups`, by what the run's
    /// taxonomy database `taxdb` can write (`texts`) and which lookups of the offset error's taxid
    /// miss its slot (`misses`, [`TaxdbRun::slot_misses`]; `None` without that error): each ask
    /// given to `answer` with its texts, `open` (the texts of opening the database, where the
    /// load opened it) with the first lookup. Nothing is recorded per taxid, and a load of a run
    /// whose taxonomy database writes no offset error and maps a `taxdb.btd` that is not empty
    /// asks nothing.
    #[inline]
    fn volume_load(
        &self,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        texts: PrefetchTexts,
        misses: Option<&[bool]>,
        lookups: impl Iterator<Item = i32>,
        mut open: Vec<String>,
        answer: &mut impl FnMut(i32, Vec<String>),
    ) -> std::result::Result<(), SeqDbError> {
        for (index, taxid) in lookups.enumerate() {
            let ask = if Some(taxid) == texts.offset_error_taxid {
                misses.is_some_and(|misses| misses[index])
            } else {
                texts.maps_null
            };
            if ask || !open.is_empty() {
                let mut messages = std::mem::take(&mut open);
                if ask {
                    taxdb.lookup(taxid, &mut messages)?;
                }
                answer(taxid, messages);
            }
        }
        Ok(())
    }

    /// The run of a program that reads sequences without NCBI's object manager: dustmasker's and
    /// segmasker's `CMaskBDBReader`, blastdbcmd's formatters. Each read of a sequence is a load
    /// (`CSeqDB::GetBioseq`), whose lookups follow the volume's cache in the order of the reads
    /// ([`TaxdbRun::bioseq_texts`]); the names, `GetTaxInfo` and `-tax_info` ask the same
    /// taxonomy database ([`TaxdbRun::lookup`]), NCBI's one per process.
    pub(crate) fn sequential() -> Self {
        TaxdbRun {
            volume_lookups: std::sync::Mutex::new(VolumeLookups {
                sequential: true,
                ..VolumeLookups::default()
            }),
            ..TaxdbRun::default()
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1462-1466
    // `GetBioseq` makes the sources of a Bioseq it gives (one with definition lines and, with its
    // data, letters) through `x_GetTaxonomy` (cited at `SlotHolds::walk`), whose `GetTaxNames`
    // writes its texts through `ERR_POST` where it reads them.
    /// The `ERR_POST` texts of `CSeqDB::GetBioseq` of NCBI's OID `oid` of `db` in a sequential run
    /// ([`TaxdbRun::sequential`]), whose definition lines are `set` and target `target`: the
    /// lookups of its positive taxids in the order of `x_GetTaxonomy`
    /// (`DeflineSet::taxonomy_order`, the target's line first) write the texts of opening the
    /// taxonomy database at the run's first lookup, and the offset error where a lookup of its
    /// taxid misses the volume's cache ([`VolumeLookups::offset_slot`]). Err: a lookup that reads a
    /// record of an empty `taxdb.btd` (NCBI maps it as a null pointer and faults) or a taxonomy
    /// database that does not open or read: LOSAT's refusal. A record read past NCBI's mapping of
    /// a `taxdb.btd` that is not empty is not refused here: nothing of its names is written (NCBI
    /// reads memory past its mapping there; whether it faults depends on its memory map). One path
    /// for every load, the search's ([`TaxdbRun::volume_load`]): the taxonomy database is asked
    /// only where a lookup can write a text or fault; a load of a run without the offset error
    /// asks nothing and allocates nothing. Inlined up to its first test: a sequence without a
    /// positive taxid, or a run without `taxdb.bti` on the search path once a load found none,
    /// asks nothing (the readers ask it of every sequence). After the run's first load, a load of
    /// a run whose taxonomy database has no offset error makes its asks here, in the reader's
    /// code: there is no slot to follow.
    #[inline]
    pub(crate) fn bioseq_texts(
        &self,
        db: &SeqDb,
        oid: u32,
        set: &crate::seqid::DeflineSet,
        target: Option<&crate::seqid::ParsedSeqId>,
    ) -> std::result::Result<Vec<String>, SeqDbError> {
        if self.found.get() == Some(&false) {
            return Ok(Vec::new());
        }
        let taxids = set.taxonomy_order(target);
        if let (Some(taxdb), Some(&prefetch)) = (self.lazy.get(), self.texts.get()) {
            if prefetch.offset_error_taxid.is_none() {
                let mut texts = Vec::new();
                self.volume_load(
                    taxdb,
                    prefetch,
                    None,
                    taxids,
                    Vec::new(),
                    &mut |_, messages| texts.extend(messages),
                )?;
                return Ok(texts);
            }
        }
        self.bioseq_loads(db, oid, taxids)
    }

    /// [`TaxdbRun::bioseq_texts`] of the run's first load of a sequence with a positive taxid
    /// (it resolves `taxdb.bti` and opens the database) and of each load where the slot of the
    /// offset error's taxid is followed (out of line: the readers' own code holds the rest).
    #[inline(never)]
    fn bioseq_loads(
        &self,
        db: &SeqDb,
        oid: u32,
        taxids: crate::seqid::defline_set::TaxonomyOrder<'_>,
    ) -> std::result::Result<Vec<String>, SeqDbError> {
        // NCBI's first lookup opens the database: at the first sequence with a positive taxid.
        if taxids.clone().next().is_none() {
            return Ok(Vec::new());
        }
        let mut open = Vec::new();
        let (taxdb, prefetch) = match (self.lazy.get(), self.texts.get()) {
            (Some(taxdb), Some(prefetch)) => (taxdb, *prefetch),
            _ => {
                let Some(search_path) = db.search_path().filter(|_| self.found(db)) else {
                    return Ok(Vec::new());
                };
                // The run's first lookup opens the database (its texts of opening come first);
                // the offset error's taxid and the null mapping are the same database's (no
                // second open, as the search's prefetch makes before its first lookup).
                let taxdb = self.lazy_taxdb(search_path);
                open = taxdb.open()?;
                let prefetch = PrefetchTexts::of(taxdb.database()?);
                let _ = self.texts.set(prefetch);
                (taxdb, prefetch)
            }
        };
        // Only the slot of the offset error's taxid is followed: its volume, where the run's
        // taxonomy database has that error.
        let misses = match prefetch.offset_error_taxid {
            Some(last) => {
                let volume = db.oid_volume(oid)?;
                let lookups: Vec<i32> = taxids.clone().collect();
                Some(
                    self.slot_misses(volume, oid, &lookups, last)
                        .map_err(|cause| slot_uncertain_error(db, volume, last, cause))?,
                )
            }
            None => None,
        };
        let mut texts = Vec::new();
        self.volume_load(
            taxdb,
            prefetch,
            misses.as_deref(),
            taxids,
            open,
            &mut |_, messages| texts.extend(messages),
        )?;
        Ok(texts)
    }

    /// `CSeqDBTaxInfo::GetTaxNames` of `taxid` through the run's taxonomy database (`GetTaxInfo`, the
    /// names of a title, `-tax_info`): the texts of opening it, when this call opens it, and of the
    /// lookup are appended to `messages`. [`TaxLookup::NotInIndex`] without a search path.
    ///
    /// [`TaxLookup::NotInIndex`]: crate::seqdb::taxdb::TaxLookup::NotInIndex
    pub(crate) fn lookup(
        &self,
        db: &SeqDb,
        taxid: i32,
        messages: &mut Vec<String>,
    ) -> std::result::Result<crate::seqdb::taxdb::TaxLookup, SeqDbError> {
        match db.search_path() {
            Some(search_path) => self.lazy_taxdb(search_path).lookup(taxid, messages),
            None => Ok(crate::seqdb::taxdb::TaxLookup::NotInIndex),
        }
    }

    /// The run's texts of the loads that NCBI does not make (`report::query_warnings::LoadTexts`),
    /// where its lookups can write a text (the offset error, or texts of opening the database):
    /// the model then records what each read's load looks up (`subject_read_texts`).
    fn load_texts(
        &self,
        db: &SeqDb,
    ) -> Option<std::sync::Arc<crate::report::query_warnings::LoadTexts>> {
        let texts = self.prefetch_texts(db);
        (texts.offset_error_taxid.is_some() || texts.opening).then(|| self.ledger())
    }

    /// The run's ledger of the texts of the loads that NCBI does not make (`load_texts`), also
    /// where the lookups write none: a load whose Seq-align id's lookup reads a UTF8String header
    /// leaves NCBI's warning there (`ReadLoad::of_scope_utf8`).
    fn ledger(&self) -> std::sync::Arc<crate::report::query_warnings::LoadTexts> {
        std::sync::Arc::clone(&self.load_texts)
    }

    /// [`TaxdbRun::volume_lookups_of`] of a load, with the length of the start of its text that
    /// is the texts of opening the database (the run's first load writes them before the texts of
    /// its asks, `TaxdbRun::load`; the only text of a lookup is the offset error).
    #[allow(clippy::too_many_arguments)]
    fn load_text(
        &self,
        db: &SeqDb,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        ordinal: u32,
        lookups: &[i32],
        planned: Option<&[bool]>,
        lower: &str,
        program: &str,
    ) -> Result<(Vec<u8>, usize)> {
        let mut text = String::new();
        let mut opening = 0;
        for (_, messages) in self.load(db, taxdb, ordinal, lookups, planned, program)? {
            for line in messages {
                let opens = line != crate::seqdb::taxdb::MSG_OFFSET_ERROR;
                text.push_str(&crate::seqdb::taxdb::err_post_line(lower, &line));
                text.push('\n');
                if opens {
                    opening = text.len();
                }
            }
        }
        Ok((text.into_bytes(), opening))
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
        for (_, messages) in self.load(db, taxdb, ordinal, lookups, None, program)? {
            for line in messages {
                text.push_str(&crate::seqdb::taxdb::err_post_line(lower, &line));
                text.push('\n');
            }
        }
        Ok(())
    }
}

/// LOSAT's refusal of a load in the volume `volume` of `db` whose texts of the offset error of
/// taxid `last` NCBI's scope decides ([`SlotUncertain`], [`VolumeLookups::offset_slot`]): the
/// reader's refusal of `taxdb.bti` that names the read and the load of the other taxid of the
/// slot before it ([`taxdb_error`] makes it a [`TaxonomyRefusal`] that names the program).
fn slot_uncertain_error(db: &SeqDb, volume: usize, last: i32, cause: SlotUncertain) -> SeqDbError {
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
    if let Some(order) = cause.order {
        return SeqDbError::Unsupported {
            path,
            message: slot_order_message(db, &volume_name, last, cause, order),
        };
    }
    let read = if cause.again {
        format!(
            "reads database sequence {} of the volume {volume_name}, which holds it, again",
            db.oid(cause.ordinal)
        )
    } else {
        format!(
            "fetches database sequence {} of the volume {volume_name}, which holds it,",
            db.oid(cause.ordinal)
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
        }) => format!(
            "after the run's fetch of database sequence {} with taxid {taxid}",
            db.oid(ordinal)
        ),
        Some(SlotFetch {
            ordinal,
            taxid,
            again: true,
        }) => format!(
            "after the run read database sequence {} again, whose load looks up taxid {taxid}",
            db.oid(ordinal)
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
            db.oid(fetch.ordinal)
        ),
        None => "there or not as its scope kept the sequences read before".to_string(),
    };
    SeqDbError::Unsupported {
        path,
        message: format!(
            "taxid {last}, the last entry of the taxonomy database, begins past the end of \
             taxdb.btd, and the run {read} {after}, which takes the same slot of NCBI \
             BLAST+'s 256-slot taxonomy cache of a volume (CSeqDBVol::m_TaxCache, \
             CSeqDBIntCache): NCBI writes \"{}\" {scope}, which LOSAT does not reproduce",
            crate::seqdb::taxdb::MSG_OFFSET_ERROR
        ),
    }
}

/// The reason of [`slot_uncertain_error`] where the order of a batch's prefetch is the cause
/// ([`SlotOrder`]): of the loads of one key whose order changes the number of offset errors, or of
/// a later load after them.
fn slot_order_message(
    db: &SeqDb,
    volume_name: &str,
    last: i32,
    cause: SlotUncertain,
    order: SlotOrder,
) -> String {
    let pair = if order.first == order.second {
        format!(
            "database sequence {} of the volume {volume_name} twice",
            db.oid(order.first)
        )
    } else {
        format!(
            "database sequences {} and {} of the volume {volume_name}",
            db.oid(order.first),
            db.oid(order.second)
        )
    };
    let taxa = if order.taxid == last {
        format!("taxid {last}")
    } else {
        format!("taxids {last} and {}", order.taxid)
    };
    let what = if cause.again && (cause.ordinal == order.first || cause.ordinal == order.second) {
        format!(
            "the prefetch of a query batch reads {pair}, database sequence {} again (NCBI's data \
             loader may keep it from an earlier load)",
            db.oid(cause.ordinal)
        )
    } else if cause.again {
        format!(
            "the run reads database sequence {} again (NCBI's data loader may keep it from an \
             earlier load), after the prefetch of a query batch read {pair}",
            db.oid(cause.ordinal)
        )
    } else if cause.ordinal == order.first {
        format!("the prefetch of a query batch reads {pair}")
    } else {
        format!(
            "the run fetches database sequence {}, which holds it, after the prefetch of a query \
             batch read {pair}",
            db.oid(cause.ordinal)
        )
    };
    format!(
        "taxid {last}, the last entry of the taxonomy database, begins past the end of taxdb.btd, \
         and {what}: their lookups of {taxa} take the same slot of NCBI BLAST+'s 256-slot taxonomy \
         cache of a volume (CSeqDBVol::m_TaxCache, CSeqDBIntCache), and NCBI loads them in the \
         order of the addresses of their Seq-id handles, which decides how often it writes \
         \"{}\"; LOSAT does not reproduce that order",
        crate::seqdb::taxdb::MSG_OFFSET_ERROR
    )
}

/// LOSAT's refusal of a taxonomy database record (`read_reported_taxonomy`): a name NCBI reads
/// past its mapping of `taxdb.btd`, or a record of an empty one (NCBI faults). A run that
/// holds its XML2 or JSON document drops it on this error and no other (BLASTX,
/// `XmlRefusal::of_taxonomy`); the other errors of the lookups (a read error of a file) are not
/// refusals, and the held document is written as NCBI writes what it has.
#[derive(Debug)]
pub(crate) struct TaxonomyRefusal {
    text: String,
    /// The refusal is of a load of a subject with a taxid where a file of the taxonomy database
    /// does not open (`SeqDbError::TaxdbThrows`): NCBI's load throws there, which its prefetch
    /// catches ([`batch_prefetch`]) and its reports do not read as LOSAT does.
    load_throws: bool,
}

impl TaxonomyRefusal {
    /// Whether `error` is the refusal of a load that throws where a file of the taxonomy database
    /// does not open ([`TaxonomyRefusal::load_throws`]).
    fn of_load_that_throws(error: &anyhow::Error) -> bool {
        error
            .downcast_ref::<TaxonomyRefusal>()
            .is_some_and(|refusal| refusal.load_throws)
    }
}

impl std::fmt::Display for TaxonomyRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.text)
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
    if let SeqDbError::TaxdbThrows { .. } = &error {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1155
        // The Bioseq's taxonomy catches only a `CSeqDBException`: the load throws `CMemoryFile`'s
        // exception, and each report then reads the subject as its own `catch` decides (line WQA
        // D18b probe 12: `Unknown` ids and titles, `N/A` names, alignments without sequences, or
        // the exception, exit status 255), which LOSAT does not write. The formats that NCBI
        // writes the same are not refused (`report::db_subject::unchanged_by_load_throws`); the
        // words say what the user can do.
        return anyhow::Error::new(TaxonomyRefusal {
            text: format!(
                "{error}; NCBI BLAST+ 2.17.0 then fails at each taxonomy lookup and writes this \
                 output differently or stops, which LOSAT's {program} does not support; make the \
                 file readable, or move it off the BLASTDB path"
            ),
            load_throws: true,
        });
    }
    let refused = matches!(error, SeqDbError::Unsupported { .. });
    let error = database_error(error, program);
    if refused {
        anyhow::Error::new(TaxonomyRefusal {
            text: format!("{error:#}"),
            load_throws: false,
        })
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
/// taxid, `taxonomy` is empty. A taxonomy database that LOSAT cannot read as NCBI does is a
/// refusal made here, before any byte is written (a scientific name of any bytes is written as
/// NCBI's streams write it); so is a run whose fetches NCBI's taxonomy cache would make repeat the
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
pub(crate) fn read_reported_taxonomy<P>(
    db: &SeqDb,
    run: &TaxdbRun,
    ordinals: &[u32],
    reported: &std::collections::BTreeSet<usize>,
    parsed: &[Option<P>],
    program: &str,
) -> Result<(
    Vec<i32>,
    Vec<crate::report::blast_xml::HitTaxonomy>,
    Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
    TaxdbNotes,
)>
where
    P: std::ops::Deref<Target = crate::report::db_subject::ParsedSubject>,
{
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
            // The scientific name reaches NCBI's stream as it is (raw bytes of `taxdb.btd`): the
            // XML stream writes or throws on them (`report::serial::xml`), JSON writes them.
            let answer = match found {
                Some(names) => HitTaxonomy::Named(names.scientific_name),
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
        for (taxid, messages) in
            run.load(db, &taxdb, ordinals[rank], &volume_lookups, None, program)?
        {
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

/// The `[BLAST] BLASTDB` value of a run's application settings (`ApplicationSettings::blastdb`):
/// the registry's part of the search path through which a `-subject` search's reports look for the
/// taxonomy database (`subject_taxonomy`). The command line takes it from its settings; the library
/// API gives it to a search of records whose reports look every taxid up (`api::run`, with the
/// call's start); LOSAT Web gives none (it has no files: the database is not found).
#[derive(Clone, Debug)]
pub struct TaxdbSetting(Result<Option<String>, String>);

impl TaxdbSetting {
    /// The value of `settings`.
    pub(crate) fn of(settings: &ApplicationSettings) -> Self {
        TaxdbSetting(settings.blastdb.clone())
    }
}

/// What the tabular reports of a `-subject` search know of the taxonomy database
/// (`subject_taxonomy`).
#[derive(Clone, Debug, Default)]
pub(crate) struct SubjectTaxonomy {
    /// The organism report's database (outfmt 18, `SubjectTaxdb`), when `taxdb.bti` is on the
    /// search path.
    pub(crate) taxdb: Option<SubjectTaxdb>,
    /// Whether `taxdb.bti` is on the search path; `None` when no report looks every taxid up.
    pub(crate) found: Option<bool>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:582-606
/// What the tabular reports `specs` of a `-subject` search of `program` (upper case) ask of the
/// taxonomy database where NCBI makes their formatters: whether `taxdb.bti` is on the search path
/// of `setting` (`x_CheckTaxDB`'s `SeqDB_ResolveDbPath` of `taxdb.bti`: the current directory,
/// `BLASTDB` and the registry's `[BLAST] BLASTDB`, as a database is looked for, `seqdb::SearchPath`)
/// for the plural names, which warn without it, and the organism report, whose lookups open it
/// (`SubjectTaxdb`). A report of singular names reads nothing: NCBI looks a singular name up by
/// the subject's taxid, and the records of `-subject` have none. A registry value that LOSAT does
/// not read as NCBI is its refusal, as for `-db`. `setting` is `None` in LOSAT Web: not found.
///
/// The command line and the library's searches of records call it at the same point, after the
/// options are processed (each program's `open_search_database`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-134
pub(crate) fn subject_taxonomy<'s>(
    setting: Option<&TaxdbSetting>,
    program: &str,
    specs: impl IntoIterator<Item = &'s crate::report::tabular::TabularSpec>,
) -> Result<SubjectTaxonomy> {
    let (mut looks_up, mut tax_report) = (false, false);
    for spec in specs {
        looks_up |= spec.looks_up_every_taxid();
        tax_report |= spec.kind == crate::report::tabular::TabularKind::TaxReport;
    }
    if !looks_up {
        return Ok(SubjectTaxonomy::default());
    }
    let Some(setting) = setting else {
        return Ok(SubjectTaxonomy {
            taxdb: None,
            found: Some(false),
        });
    };
    let path = subject_search_path(setting, program)?;
    let found = path.resolve_file("taxdb.bti").is_some();
    // The organism report's lookups open the run's own database.
    let taxdb = (tax_report && found).then(|| {
        SubjectTaxdb(std::sync::Arc::new(crate::seqdb::taxdb::LazyTaxDb::new(
            path,
        )))
    });
    Ok(SubjectTaxonomy {
        taxdb,
        found: Some(found),
    })
}

/// The search path of NCBI's databases of a `-subject` search (`subject_taxonomy`).
fn subject_search_path(setting: &TaxdbSetting, program: &str) -> Result<SearchPath> {
    let registry = match &setting.0 {
        Ok(value) => value.as_deref(),
        Err(reason) => bail!("{reason}; this is not supported by LOSAT's {program}"),
    };
    SearchPath::from_environment(registry).map_err(|error| database_error(error, program))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-220
/// The taxonomy database of a `-subject` search's organism report (outfmt 18): NCBI's
/// `GetTaxNames` opens its static `CTaxDBFileInfo` through the search path of NCBI's databases at
/// the report's first lookup. One per run, as a `-db` search's (`TaxdbRun`), made by
/// `subject_taxonomy` (LOSAT Web gives no search path: none, every lookup fails).
#[derive(Clone, Debug)]
pub struct SubjectTaxdb(std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>);

impl SubjectTaxdb {
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
                Some(std::sync::Arc::clone(run.rows_open())),
                taxdb_error,
            )
            .where_loads_throw(run.load_throws(db)),
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
    let opened = std::mem::take(&mut *run.rows_open().lock().expect("the texts of opening"));
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
/// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
///
/// Where a subject's Seq-align id has a lookup that reads a header with a UTF8String where a
/// VisibleString is expected (`ParsedSubject::scope_utf8`), NCBI's warning (803.10) comes at that
/// subject's load, in the order of `subjects` (a query's report, whose batch's prefetch is off,
/// writes the texts of its reads before the report: `load_read_texts`).
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
    Ok(reads.into_iter().flat_map(|read| read.text).collect())
}

/// The loads of a batch's prefetch of `subjects` (each `(ordinal, parsed Seq-ids)`;
/// [`subject_read_texts`] with no key): [`batch_prefetch`]'s. Where the load of a subject's
/// Seq-align id posts NCBI's warning of a UTF8String header (`ParsedSubject::scope_utf8`,
/// `load_read_texts`), the loads come in NCBI's order ([`scope_utf8_prefetch_loads`]); else in
/// the order of `subjects` (the bytes of a set of loads do not depend on their order:
/// [`prefetch_messages`]). `at_posts`: the loads' texts are made again where the batch's posts are
/// written (`diag_log::Posts::loads`), where a load that asks the taxonomy database takes the texts
/// that a load NCBI does not make left (`ReadLoad::placed`).
fn prefetch_loads<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: impl IntoIterator<Item = (u32, Option<&'s crate::report::db_subject::ParsedSubject>)>,
    at_posts: bool,
    order: LoadsOrder<'_>,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    let subjects: Vec<(u32, Option<&'s crate::report::db_subject::ParsedSubject>)> =
        subjects.into_iter().collect();
    if db.scope_utf8_loads() && scope_utf8_among(&subjects) {
        return scope_utf8_prefetch_loads(db, run, &subjects, at_posts, order, program);
    }
    load_read_texts(
        db,
        run,
        subjects
            .into_iter()
            .map(|(ordinal, parsed)| (0, ordinal, SubjectLoad::Read(parsed))),
        false,
        order,
        program,
    )
}

/// Whether one of `subjects` has a Seq-align id whose lookup reads a header with a UTF8String
/// where a VisibleString is expected (`ParsedSubject::scope_utf8`).
fn scope_utf8_among(subjects: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)]) -> bool {
    subjects
        .iter()
        .any(|(_, parsed)| parsed.is_some_and(|subject| subject.scope_utf8_read().is_some()))
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3346-3366 (cited at `load_before_stop`)
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-406 (cited
// at `report::db_subject::ParsedSubject::scope_utf8` and above `StopUnsure`)
// The prefetch loads the batch's Seq-align ids in the order of a map of Seq-id handles (by packed
// number, then by the address of each handle's info, `load_before_stop`); each load looks its id
// up (`x_GetOid`: the version check and the header of the first OID, where the warning comes)
// before it loads the blob, whose taxonomy writes the load's texts (`TaxdbRun::load`).
/// [`prefetch_loads`] where a load posts NCBI's warning of a UTF8String header: the loads in
/// NCBI's order of their Seq-id handles ([`prefetch_load_key`]; a batch of more than
/// [`PREFETCH_ONE_MAP_IDS`] ids, or two subjects of one Seq-align id, in an order that LOSAT does
/// not know), so that the warning comes after the texts of the loads before it. The batch is
/// refused where a load whose order to the warning's comes from the addresses of NCBI's Seq-id
/// infos (a handle of the same packed number, both unpacked included; the probe of line DBID
/// loadord: `gnl|BL_ORD_ID|3` before `ref|NP_000003.1|` with one command line and after it with
/// another) or is not known could change the bytes: another load whose lookup posts a warning that
/// names another member (one of the same member writes the same bytes, fix 5), or one that asks
/// the taxonomy database (in a run whose lookups can write a text, `ReadLoad`)
/// where the batch's loads write a text of it or may take one that a report before left
/// (`at_posts`).
#[cold]
#[inline(never)]
fn scope_utf8_prefetch_loads(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    at_posts: bool,
    loads_order: LoadsOrder<'_>,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    let places = prefetch_places(db, subjects);
    let mut order: Vec<usize> = (0..subjects.len()).collect();
    order.sort_by_key(|&index| places[index].order);
    let mut reads = load_read_texts(
        db,
        run,
        order.iter().map(|&index| {
            (
                index as u32,
                subjects[index].0,
                SubjectLoad::Read(subjects[index].1),
            )
        }),
        false,
        loads_order,
        program,
    )?;
    let lead = |read: &crate::report::query_warnings::SubjectRead| {
        read.load.as_ref().map_or(0, |load| load.lead)
    };
    let Some(posting) = reads.iter().position(|read| lead(read) > 0) else {
        return Ok(reads);
    };
    let at = places[reads[posting].subject as usize];
    // `at_posts`: the texts that the loads NCBI did not make before left (a BLAST XML report that
    // the stream cut, `ReadLoad::left_out`) go to the first load that asks for them where the
    // posts are written, which the model of the batch does not know yet.
    let taxonomy_texts = at_posts || reads.iter().any(|read| read.text.len() > lead(read));
    for (index, read) in reads.iter().enumerate() {
        let place = places[read.subject as usize];
        let Some(why) = place.order_unknown(at).filter(|_| index != posting) else {
            continue;
        };
        // Another load whose lookup reads such a header changes what the prefetch writes only
        // where its warning names another member (the run writes the first load's).
        let scoped = subjects[read.subject as usize]
            .1
            .is_some_and(|subject| subject.scope_utf8_read().is_some())
            && scope_utf8_path(db, subjects[read.subject as usize].1)
                != scope_utf8_path(db, subjects[reads[posting].subject as usize].1);
        let asks = read
            .load
            .as_ref()
            .is_some_and(|load| load.looks_up || load.offset);
        if scoped || (asks && taxonomy_texts) {
            let (ordinal, parsed) = subjects[reads[posting].subject as usize];
            let other = parsed
                .and_then(|subject| subject.scope_utf8_read())
                .expect("the load that posts the warning reads such a header");
            return Err(parsed_refusal(
                db.oid(ordinal),
                &format!(
                    "has a Seq-align id by which NCBI BLAST+'s lookup gives database sequence {other} first, whose header holds a UTF8String where a VisibleString is expected: NCBI BLAST+ 2.17.0's prefetch writes its warning (803.10) at the load of that id, and the load of database sequence {} ({}) may come before it or after it ({why}), which changes what the prefetch writes before the warning",
                    db.oid(subjects[read.subject as usize].0),
                    if scoped {
                        "whose lookup reads such a header with another member"
                    } else {
                        "which asks the taxonomy database"
                    },
                ),
                program,
            ));
        }
    }
    for read in &mut reads {
        read.subject = 0;
    }
    Ok(reads)
}

/// Why the order of two loads of NCBI's prefetch whose Seq-id handles have the same packed number
/// (or are both unpacked) is not known ([`scope_utf8_prefetch_loads`]).
const ADDRESS_ORDER: &str = "their Seq-id handles have the same packed number, or neither is packed, and NCBI's map orders such handles by the address of their info";

/// The place of each of `subjects` (a batch's prefetch's loads) in NCBI's map of Seq-id handles
/// ([`LoadPlace`] without the stops, as [`load_before_stop`] makes it).
fn prefetch_places(
    db: &SeqDb,
    subjects: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
) -> Vec<LoadPlace> {
    let prefetched: Vec<PrefetchSubject<'_>> = subjects
        .iter()
        .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed))
        .collect();
    load_places(db, &prefetched).0
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:819-838
/// [`prefetch_messages`] of `subjects`, each `(key, ordinal, parsed Seq-ids)`, by subject: the
/// text that the first read of each subject writes, with the subject's `key` (the name a report
/// writer gives it), in the order of `subjects`; a subject whose read writes nothing has no
/// entry, unless the run's lookups can write a text (`TaxdbRun::load_texts`): then each read
/// that loads a blob has one, with what its load looks up (`ReadLoad`), and its text holds the
/// texts that loads NCBI did not make left to it (a BLAST XML report that the stream cut,
/// `report::query_warnings::ReadCursor::cut`, before this batch). A query's report whose batch's
/// prefetch is off writes each where it reads the subject
/// (`report::query_warnings::SubjectReads`). `rows_open`: the report's rows look the plural
/// names up (taxid 0 too, `TaxdbRun::open_for_rows`), each row after it loads its subject, so
/// the first row opens the taxonomy database after the first subject's read, where that read
/// did not.
pub(crate) fn subject_read_texts<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: impl IntoIterator<
        Item = (
            u32,
            u32,
            Option<&'s crate::report::db_subject::ParsedSubject>,
        ),
    >,
    rows_open: bool,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    load_read_texts(
        db,
        run,
        subjects
            .into_iter()
            .map(|(key, ordinal, parsed)| (key, ordinal, SubjectLoad::Read(parsed))),
        rows_open,
        LoadsOrder::Reads,
        program,
    )
}

/// How NCBI orders the loads of one call of [`load_read_texts`], as far as it decides the texts of
/// the offset error on a database whose volume holds its taxid and another taxid of its slot
/// ([`TaxdbRun::cache_collision`]).
#[derive(Clone, Copy)]
enum LoadsOrder<'p> {
    /// A report's reads (also those of a prefetch that stops, [`prefetch_stop`]): on such a
    /// database LOSAT's refusal, before the first load ([`TaxdbRun::refuse_cache_collision`]).
    /// A BLAST XML report that the stream cuts leaves its later reads out, which moves their
    /// texts and not what they put in the slot (`report::query_warnings::ReadLoad::left_out`),
    /// and some loads before a stop may not be made.
    Reads,
    /// A batch's prefetch ([`batch_prefetch`]): the misses of each load in the slot of the offset
    /// error's taxid that NCBI's order of the loads gives (`planned`, by the loaded blob's ordinal
    /// id, [`prefetch_slot_plan`]; `None` where no volume holds that taxid and another of its slot,
    /// whose loads give the same texts in any order and follow the slot one by one).
    Prefetch(Option<&'p PrefetchSlots>),
}

/// What one load of [`load_read_texts`] reads.
#[derive(Clone, Copy)]
enum SubjectLoad<'s> {
    /// The read of a subject by its Seq-align id (the definition lines of one with parsed
    /// Seq-ids), or by its ordinal id (`None`).
    Read(Option<&'s crate::report::db_subject::ParsedSubject>),
    /// The load of a sequence that NCBI's BLAST database data loader gives for the Seq-align id
    /// of the stopping subject `stop` (`by`) without its lookup, as an earlier load registered
    /// that id (`RegisterIds`; [`StopUnsure::instead`]): the sequence's own blob, whose taxids
    /// the loader asks for by that id ([`blob_taxids`]).
    Instead {
        stop: u32,
        by: &'s crate::report::db_subject::ParsedSubject,
    },
}

/// [`subject_read_texts`] of the loads `loads`, each `(key, ordinal, what it reads)`: also the
/// loads of the sequences that a skipped stopping lookup is given instead ([`prefetch_stop`]).
/// The first load of a subject whose Seq-align id's lookup reads a header with a UTF8String where a
/// VisibleString is expected (`ParsedSubject::scope_utf8`) posts NCBI's warning (803.10) where the
/// run has not: the start of its text (`ReadLoad::lead`), before the texts of its blob's taxonomy,
/// also without `taxdb.bti` on the search path. The loads come in NCBI's order (the callers': a
/// report's reads, a prefetch's Seq-id handles, [`prefetch_loads`], [`prefetch_stop`]); `order` says
/// how that order decides the misses in the slot of the offset error's taxid ([`LoadsOrder`]).
fn load_read_texts<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    loads: impl IntoIterator<Item = (u32, u32, SubjectLoad<'s>)>,
    rows_open: bool,
    order: LoadsOrder<'_>,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    use crate::report::query_warnings::{ReadLoad, SubjectRead};
    let lower = program.to_ascii_lowercase();
    let taxonomy = run.found(db).then(|| db.search_path()).flatten();
    // NCBI's warning of a UTF8String header that a load's lookup reads (`ParsedSubject::scope_utf8`):
    // only in a run whose reader met such a subject and that has not posted it.
    let scope_utf8 = db.scope_utf8_loads();
    let Some(search_path) = taxonomy else {
        if !scope_utf8 {
            return Ok(Vec::new());
        }
        return Ok(scope_utf8_reads(db, run, loads, &lower));
    };
    let planned = match order {
        LoadsOrder::Reads => {
            run.refuse_cache_collision(db, program)?;
            None
        }
        LoadsOrder::Prefetch(planned) => planned,
    };
    let taxdb = run.taxdb(search_path);
    let texts = run.load_texts(db);
    let offset_taxid = run.prefetch_texts(db).offset_error_taxid;
    let mut reads = Vec::new();
    let mut fetched: std::collections::HashSet<u32> = Default::default();
    let mut thrown: std::collections::HashSet<u32> = Default::default();
    for (key, ordinal, load) in loads {
        // A subject whose lookup throws (`ParsedSubject::lookup_throws`, a report that reads one by
        // one reaches it: `HitLoads::reached`) has read the candidates' headers first, so its load
        // posts NCBI's warning of a UTF8String among them and nothing else: no blob is read, no
        // taxonomy text comes, and a row of it opens nothing (`rows_open` waits for the first
        // subject that loads).
        if let SubjectLoad::Read(Some(subject)) = load {
            if subject.lookup_throws().is_some() {
                if thrown.insert(ordinal) && scope_utf8 {
                    reads.extend(scope_utf8_read(db, run, key, subject, &lower));
                }
                continue;
            }
        }
        if !fetched.insert(ordinal) {
            continue;
        }
        let first = fetched.len() == 1;
        let (parsed, lookups_of) = match load {
            SubjectLoad::Read(parsed) => (parsed, None),
            SubjectLoad::Instead { stop, by } => (None, Some((stop, by))),
        };
        // The lookup of the Seq-align id comes before the blob's load (`x_GetOid`).
        let lead = if scope_utf8 {
            scope_utf8_lead(db, parsed, &lower)
        } else {
            Vec::new()
        };
        let scoped = parsed.is_some_and(|subject| subject.scope_utf8_read().is_some());
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
                Err(reason) => return Err(parsed_refusal(db.oid(ordinal), &reason, program)),
            },
            None => Some(ordinal),
        };
        let (text, load) = match blob {
            Some(blob) => {
                let lookups = match lookups_of {
                    Some((stop, by)) => blob_taxids(db, stop, Some(by), blob, program)?,
                    None => blob_taxids(db, ordinal, parsed, blob, program)?,
                };
                let misses = planned.map(|planned| planned.misses(ordinal, lookups.len()));
                let (text, opening) = run.load_text(
                    db,
                    &taxdb,
                    blob,
                    &lookups,
                    misses.as_deref(),
                    &lower,
                    program,
                )?;
                let load = texts.as_ref().map(|texts| {
                    Ok::<_, anyhow::Error>(ReadLoad {
                        texts: std::sync::Arc::clone(texts),
                        volume: db
                            .volume_index(blob)
                            .map_err(|error| database_error(error, program))?
                            as u32,
                        looks_up: !lookups.is_empty(),
                        offset: offset_taxid.is_some_and(|taxid| lookups.contains(&taxid)),
                        lead: 0,
                        opening,
                        scope_utf8: scoped,
                    })
                });
                (text, load.transpose()?)
            }
            None => (Vec::new(), None),
        };
        // A load that NCBI does not make leaves its warning to the next that reads such a header
        // (`ReadLoad::left_out`): every load of such a subject has the facts of its load.
        let load = match load {
            None if scoped => Some(ReadLoad::of_scope_utf8(run.ledger())),
            load => load,
        };
        match load {
            Some(mut load) => {
                let text = load.place_after(lead, text);
                reads.push(SubjectRead {
                    subject: key,
                    text,
                    load: Some(load),
                });
            }
            None if !text.is_empty() || !lead.is_empty() => reads.push(SubjectRead {
                subject: key,
                text: [lead, text].concat(),
                load: None,
            }),
            None => {}
        }
        if first && rows_open {
            run.open_for_rows(db, program)?;
        }
    }
    Ok(reads)
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393-406 (cited
// at `report::db_subject::ParsedSubject::scope_utf8`)
/// The warning (803.10) that the load of `parsed`'s Seq-align id posts, where its lookup reads a
/// header with a UTF8String where a VisibleString is expected (`ParsedSubject::scope_utf8`) and
/// the run has not posted it (`SeqDb::claim_utf8_post`, once): `Warning: (803.10) [<program>] ...`
/// naming the member of that header (`seqid::ber::utf8_string_warning`); else empty.
#[cold]
#[inline(never)]
fn scope_utf8_lead(
    db: &SeqDb,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    lower: &str,
) -> Vec<u8> {
    let Some(path) = scope_utf8_path(db, parsed) else {
        return Vec::new();
    };
    if !db.claim_utf8_post() {
        return Vec::new();
    }
    crate::seqid::ber::utf8_string_warning(Some(lower), &path)
}

/// The member of the UTF8String that the lookup of `parsed`'s Seq-align id reads in the header of
/// the sequence it gives first (`ParsedSubject::scope_utf8`, `SeqDb::utf8_string_path_at`): what
/// NCBI's warning of that load names. `None` where it reads none.
#[cold]
#[inline(never)]
fn scope_utf8_path(
    db: &SeqDb,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
) -> Option<String> {
    let other = parsed.and_then(|subject| subject.scope_utf8_read())?;
    db.location_of(other.ordinal())
        .and_then(|at| db.utf8_string_path_at(at))
}

/// The reads of a query's report whose batch's prefetch is off, for the XML2 and JSON formats
/// (whose writers write the taxonomy texts of each `Hit`'s load themselves, `TaxdbPrinter`): the
/// `Hit`s' loads in their order (`subjects`, each `(key, ordinal, parsed Seq-ids)`), of which the
/// first whose Seq-align id's lookup reads a header with a UTF8String where a VisibleString is
/// expected posts NCBI's warning there (`report::query_warnings::ReadCursor`, at the `Hit`); empty
/// where none does (`SeqDb::scope_utf8_loads`). The `Hit`s that the report reaches
/// (`HitLoads::of` of this report alone, `after`: what the run's formats do after a subject whose
/// lookup throws, `FormatThrow::of_run`): those up to the first `Hit` whose lookup throws, that
/// one included (its lookup read the candidates' headers first, `HitLoads::reached`), not those
/// after it where the report or the run ends there (line DBID pfoff, DECISIONS 62); the one owner
/// of that answer is `HitLoads`.
pub(crate) fn scope_utf8_hit_reads<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    after: Option<crate::report::db_subject::RunThrow>,
    subjects: impl IntoIterator<
        Item = (
            u32,
            u32,
            Option<&'s crate::report::db_subject::ParsedSubject>,
        ),
    >,
    program: &str,
) -> crate::report::query_warnings::QueryReads {
    if !db.scope_utf8_loads() {
        return Vec::new();
    }
    let subjects: Vec<_> = subjects.into_iter().collect();
    let loads = crate::report::db_subject::HitLoads::of(
        after,
        [(
            true,
            subjects
                .iter()
                .map(|&(key, ..)| Some(key as usize))
                .collect::<Vec<_>>(),
        )],
        |key| {
            subjects
                .iter()
                .find(|&&(subject, ..)| subject as usize == key)
                .and_then(|&(.., parsed)| parsed)
        },
        |_, _| 0,
    );
    scope_utf8_reads(
        db,
        run,
        subjects
            .iter()
            .filter(|&&(key, ..)| loads.reached(key as usize))
            .map(|&(key, ordinal, parsed)| (key, ordinal, SubjectLoad::Read(parsed))),
        &program.to_ascii_lowercase(),
    )
}

/// [`load_read_texts`] without `taxdb.bti` on the search path, in a run whose reader met a subject
/// whose Seq-align id's lookup reads a UTF8String header (`SeqDb::scope_utf8_seen`): the loads
/// write no taxonomy text, and each load of such a subject has an entry, the first posting NCBI's
/// warning (`scope_utf8_lead`).
#[cold]
#[inline(never)]
fn scope_utf8_reads<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    loads: impl IntoIterator<Item = (u32, u32, SubjectLoad<'s>)>,
    lower: &str,
) -> crate::report::query_warnings::QueryReads {
    let mut fetched: std::collections::HashSet<u32> = Default::default();
    let mut reads = Vec::new();
    for (key, ordinal, load) in loads {
        if !fetched.insert(ordinal) {
            continue;
        }
        let SubjectLoad::Read(Some(parsed)) = load else {
            continue;
        };
        reads.extend(scope_utf8_read(db, run, key, parsed, lower));
    }
    reads
}

/// The read of `parsed` (key `key`) that only NCBI's lookup of its Seq-align id makes: where that
/// lookup reads a header with a UTF8String where a VisibleString is expected
/// (`ParsedSubject::scope_utf8_read`), the warning (`scope_utf8_lead`) in a load of the run's
/// ledger; none otherwise. The one producer of the read of a subject whose lookup then throws
/// (`ParsedSubject::lookup_throws`: it read the candidates' headers first, and loads nothing of
/// the subject, `load_read_texts`) and of every load without `taxdb.bti` (`scope_utf8_reads`).
fn scope_utf8_read(
    db: &SeqDb,
    run: &TaxdbRun,
    key: u32,
    parsed: &crate::report::db_subject::ParsedSubject,
    lower: &str,
) -> Option<crate::report::query_warnings::SubjectRead> {
    use crate::report::query_warnings::{ReadLoad, SubjectRead};
    parsed.scope_utf8_read()?;
    let mut load = ReadLoad::of_scope_utf8(run.ledger());
    let text = load.place_after(scope_utf8_lead(db, Some(parsed), lower), Vec::new());
    Some(SubjectRead {
        subject: key,
        text,
        load: Some(load),
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:513
/// The texts of the taxonomy database where psiblast's PSSM of a query whose batch reads one by
/// one reads the sequence `blob` by a Seq-align id, the lookups of its load `lookups` (the
/// positive taxids in `x_GetTaxonomy`'s order, [`blob_taxids`]; `ListedTitles::pssm_load_texts`),
/// which NCBI writes there, after the query's report: the run's loads ([`TaxdbRun::load`], as
/// [`prefetch_messages`] for the reports' reads) write the texts of opening with the run's first
/// lookup and ask the taxonomy database where a lookup of the offset error's taxid misses its slot
/// of the volume's 256-slot cache, which keeps the taxid across the run's queries and rounds, and,
/// where `taxdb.btd` is mapped as a null pointer, at every lookup of another taxid
/// ([`TaxdbRun::volume_load`]). The PSSM's reads come in the order of its Seq-aligns, so the first
/// load of a sequence in the run writes NCBI's texts, a taxid of the slot looked up between two
/// lookups of that taxid included. NCBI's scope does not keep what a query loaded for its later
/// rounds or queries (`ResetScopeHistory` after each query; review 4 of line DBID pssmread: a second
/// query or round loads a sequence again after 40 other included subjects and not after none): a
/// sequence loaded before (the query's report, an earlier query or round) is loaded again or not,
/// which LOSAT does not model, so such a read is refused where the two give other texts (the slot
/// then holds another taxid) and else writes nothing. The lookups are the run's from then on. Empty
/// without `taxdb.bti` on the search path, and where every such load throws
/// ([`TaxdbRun::load_throws`]). Where the lookup of the read's Seq-align id (of the subject
/// `parsed`) reads a header with a UTF8String where a VisibleString is expected and the run has not
/// posted NCBI's warning (an off batch whose report read nothing of the subject), the read posts it
/// before those texts (`load_read_texts`); a load that throws posts it too, since the lookup reads
/// the header before the taxonomy database throws.
pub(crate) fn pssm_load_texts(
    db: &SeqDb,
    run: &TaxdbRun,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    blob: u32,
    lookups: &[i32],
    program: &str,
) -> Result<Vec<u8>> {
    let lead = if db.scope_utf8_loads() {
        scope_utf8_lead(db, parsed, &program.to_ascii_lowercase())
    } else {
        Vec::new()
    };
    // A load that throws at its first lookup writes no text of the taxonomy database
    // ([`TaxdbRun::load_throws`]).
    if lookups.is_empty() || !run.found(db) || run.load_throws(db) {
        return Ok(lead);
    }
    let Some(search_path) = db.search_path() else {
        return Ok(lead);
    };
    let taxdb = run.taxdb(search_path);
    let mut text = String::from_utf8(lead).expect("NCBI's warning is ASCII");
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

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3289-3299
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:201-208
// `LoadSequencesToScope` (cited at `check_prefetched_parsed`) first loads the blob of every
// subject of the batch (`GetBioseqHandles`), in the order of a map of Seq-id handles: the ids
// that NCBI's Seq-id mapper packs (GIs, accessions of a prefix and digits, the text tag's
// `gnl|BL_ORD_ID|0`) by their packed value, then the others by the address of each id's
// `CSeq_id_Info` (review 1 of line DBID pfwarn, M1: `gb|FFF00093.1|` loads before the text tag,
// `lcl|zzz` never), which LOSAT models only where a lookup stops the run ([`load_before_stop`]).
// A lookup that throws on its id there ends the loads, and the subjects after it are read by the
// reports, so which subjects were loaded before it is not known here; the throw is caught with
// `Error pre-fetching sequence data` (`ERR_POST`, on standard error, after the texts of the
// loads before). A lookup that fails on a file is a `CSeqDBException`, which the prefetch
// rethrows (`FoundInstead::Failed`, refused where it is read). A minus-strand Seq-align of a subject
// without a Bioseq throws after every blob was loaded (`prefetch_gets_no_bioseq`): the texts of
// the whole batch, then the warning. Either way the prefetch reads no sequence data, which no
// text comes with (probe of line DBID pfwarn, `R/pfwarn/w0`).
/// What the stop of a query batch's prefetch ([`prefetch_stop`]) asks of the batch's queries, made
/// by each program where it prefetches a batch ([`batch_prefetch`]).
#[derive(Clone, Copy)]
pub(crate) struct PrefetchQueries<'a> {
    /// The query that the stop's [`HeaderReadError`] names (the batch's first query of the
    /// searched queries of the call; BLASTX, which writes each batch's reports before it searches
    /// the next, names 0).
    pub(crate) error_query: usize,
    /// The query records that NCBI's scope holds at the batch's prefetch, asked only where the
    /// batch stops ([`stop_lookup_unsure`]): the one owner of the ids of the queries of a scope is
    /// `query_batch::QueryScope` (`QueryScope::of_shown_records`). They are the input records that
    /// NCBI's reader has added for the batch (a record that `-query_loc` drops included, and the
    /// query of a PSSM), and in outfmt 5, whose scope the formatter does not reset
    /// (`CBlastFormat::ResetScopeHistory`, cited at [`stop_lookup_unsure`]), those of the batches
    /// before too; never those of the batches after.
    scope: &'a dyn Fn() -> crate::blastinput::query_batch::QueryScope,
    /// Whether calls before this one reported sequences of the run that the run's record does not
    /// hold ([`PrefetchQueries::after_unrecorded_calls`]).
    unrecorded_calls: bool,
    /// The outfmt 5 reports of the queries before the batch, where outfmt 5 is the first format
    /// of a program that writes them after the stop ([`PrefetchQueries::after_xml_reports`]).
    xml_before: Option<&'a XmlReportsBefore<'a>>,
    /// Whether the reports of the run write the same where every load of a subject with a
    /// positive taxid throws ([`PrefetchQueries::where_reports_unchanged_by_load_throws`]).
    reports_unchanged_by_load_throws: bool,
}

impl<'a> PrefetchQueries<'a> {
    /// The batch whose stop names the query `error_query` and whose scope of queries at its
    /// prefetch `scope` makes.
    pub(crate) fn new(
        error_query: usize,
        scope: &'a dyn Fn() -> crate::blastinput::query_batch::QueryScope,
    ) -> Self {
        Self {
            error_query,
            scope,
            unrecorded_calls: false,
            xml_before: None,
            reports_unchanged_by_load_throws: false,
        }
    }

    /// The same batch of a run whose formats write the same where every load of a subject with a
    /// positive taxid throws (`unchanged`: `report::db_subject::unchanged_by_load_throws` of each
    /// format). Where a file of the taxonomy database does not open, such a load throws at the
    /// batch's prefetch, which catches it with the warning ([`batch_prefetch`]), and then at the
    /// loads of the reports, which such a run writes as LOSAT does; another run is refused there
    /// ([`TaxonomyRefusal`]).
    pub(crate) fn where_reports_unchanged_by_load_throws(self, unchanged: bool) -> Self {
        Self {
            reports_unchanged_by_load_throws: self.reports_unchanged_by_load_throws || unchanged,
            ..self
        }
    }

    /// The same batch of a call that comes after other calls of the run (`unrecorded`: the
    /// library's batches, each searched as a call of its own with its own database handle) whose
    /// loads neither the call's batches before ([`batch_prefetch`]'s `before`) nor the run's
    /// record ([`TaxdbRun::note_loaded`]) hold: every sequence of the database may have been
    /// loaded before the batch.
    pub(crate) fn after_unrecorded_calls(self, unrecorded: bool) -> Self {
        Self {
            unrecorded_calls: self.unrecorded_calls || unrecorded,
            ..self
        }
    }

    /// The same batch after the outfmt 5 reports `xml_before` of the queries before it
    /// ([`XmlReportsBefore`]; `None` for other first formats, and for BLASTX, whose model of a
    /// batch follows the reports of the batches before), asked only where the batch stops
    /// ([`prefetch_stop`], line RG K4h's E2).
    pub(crate) fn after_xml_reports(self, xml_before: Option<&'a XmlReportsBefore<'a>>) -> Self {
        Self { xml_before, ..self }
    }
}

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

/// The BLAST XML (outfmt 5) reports of the queries before a batch whose prefetch stops
/// ([`batch_prefetch`], [`PrefetchQueries::after_xml_reports`]), where outfmt 5 is the first
/// format of a program that models every batch before it writes a report (BLASTP, BLASTN,
/// TBLASTN, TBLASTX: the run's ledger of the texts left out, `report::query_warnings::LoadTexts`,
/// holds none of these reports' yet): by query, the model's reads of its report (`reads`: those
/// of the reports of a batch whose prefetch is off, which load their subjects), whether NCBI's XML
/// stream throws in its `Iteration` before the first `Hit` (`cut`,
/// `report::blast_xml::names::iteration_cut_before_hits`), and the posts written before its
/// report (`prefetched`: those of the prefetch of the batch it begins, whose loads' texts are
/// made where they are written, `diag_log::Posts::loads`). Asked only where the batch's prefetch
/// stops. Made by [`XmlCuts::reports_before`], so that its cuts are those of the run's record of
/// the batches' loads ([`LoadedBefore`]'s `cuts`).
pub(crate) struct XmlReportsBefore<'a> {
    pub reads: &'a [crate::report::query_warnings::QueryReads],
    cuts: XmlCuts<'a>,
    pub prefetched: &'a dyn Fn(usize) -> Option<&'a crate::report::diag_log::Posts>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375 (cited at
// `report::query_warnings::ReadCursor::cut`)
/// The outfmt 5 reports of a call's queries that NCBI's XML stream cuts before their first `Hit`
/// (`report::blast_xml::names::iteration_cut_before_hits` of the query's record, at its
/// `Iteration_query-ID` or `Iteration_query-def`), where outfmt 5 is the first format; none for
/// any other first format. Such a report makes none of its loads: the one answer to whether a
/// report before a stopping batch loaded its subjects, for the run's record of the batches'
/// loads ([`LoadedBefore`]'s `cuts`) and for the texts that the loads not made leave
/// ([`XmlReportsBefore`], `report::query_warnings::LoadTexts::left_before`), so that one stop
/// never holds two answers (review 6 of line DBID loadord, H2).
#[derive(Clone, Copy)]
pub(crate) struct XmlCuts<'a> {
    cut: Option<&'a dyn Fn(usize) -> bool>,
}

impl<'a> XmlCuts<'a> {
    /// No report is cut: the first format is not outfmt 5, or the record is of batches whose
    /// reports no outfmt 5 writes (the report streams of BLASTN and TBLASTX, which never write
    /// it, `x_rpt_stream`, `report_stream`).
    pub(crate) const NONE: XmlCuts<'static> = XmlCuts { cut: None };

    /// The cuts where outfmt 5 is the first format (`xml_first`): `cut(query)` by query of the
    /// call (the predicate of the query's record); else none.
    pub(crate) fn of(xml_first: bool, cut: &'a dyn Fn(usize) -> bool) -> Self {
        Self {
            cut: xml_first.then_some(cut),
        }
    }

    /// Whether NCBI's XML stream cuts the report of query `query` before its first `Hit`.
    pub(crate) fn cut(&self, query: usize) -> bool {
        self.cut.is_some_and(|cut| cut(query))
    }

    /// The outfmt 5 reports before a batch (`reads`, `prefetched`: [`XmlReportsBefore`]) with
    /// these cuts, where outfmt 5 is the first format (`None` for any other).
    pub(crate) fn reports_before(
        self,
        reads: &'a [crate::report::query_warnings::QueryReads],
        prefetched: &'a dyn Fn(usize) -> Option<&'a crate::report::diag_log::Posts>,
    ) -> Option<XmlReportsBefore<'a>> {
        self.cut.map(|_| XmlReportsBefore {
            reads,
            cuts: self,
            prefetched,
        })
    }
}

impl XmlReportsBefore<'_> {
    /// The texts that the loads NCBI did not make (a report that the stream cut before its hits,
    /// `report::query_warnings::ReadLoad::left_out`) left and that no load of the reports and
    /// prefetches of the queries before `query` took (`ReadLoad::placed`), where they are written
    /// in NCBI's order (`report::query_warnings::LoadTexts::left_before`, the writers' model);
    /// `None` where none is left.
    fn left_texts(
        &self,
        query: usize,
    ) -> Option<std::sync::Arc<crate::report::query_warnings::LoadTexts>> {
        let left = crate::report::query_warnings::LoadTexts::left_before(
            self.reads,
            &|before| self.cuts.cut(before),
            self.prefetched,
            query,
        );
        (!left.is_empty()).then_some(left)
    }
}

/// NCBI's prefetch of one query batch's reported subjects on standard error, for a first format
/// that prefetches ([`prefetches_subjects`]) and a batch whose prefetch is on: the posts of the
/// taxonomy database's texts of the loads (those of [`prefetch_messages`], none without
/// `taxdb.bti` on the search path), then `Warning: [<program>] Error pre-fetching sequence data`
/// where the prefetch fails. The loads' texts are made where the posts are written
/// (`diag_log::Posts::loads`): a BLAST XML report of a batch before that the stream cut left the
/// texts of loads NCBI did not make to them (`report::query_warnings::ReadCursor::cut`), which
/// a program that models every batch before it writes a report cannot know here.
/// `subjects` are the batch's reported sequences ([`PrefetchSubject`]), in any order. The
/// prefetch fails at a subject whose lookup throws on its id (`ParsedSubject::prefetch_fails`),
/// which loads nothing, and at a flipped Seq-align of a subject without a Bioseq
/// ([`prefetch_gets_no_bioseq`]). A batch with a failed lookup whose other subjects' loads write
/// a text is refused: LOSAT does not model which of them NCBI loads before the failure. The caller
/// writes the text where the batch's prefetch messages go, after the reader's messages of the
/// batch and before the warnings of its first query.
///
/// A subject whose lookup fails on a file of the database with NCBI's own exception
/// (`ParsedSubject::prefetch_stops`: a list of version 5 databases with and without the accession
/// table) stops the run there: `Err` with a [`HeaderReadError`] of the query that `queries` names
/// (`PrefetchQueries::error_query`) whose error is NCBI's `BLAST Database error: <message>`
/// (exit status 2), after the texts of the loads that NCBI's prefetch makes before it in the order
/// of the batch's Seq-id handles ([`load_before_stop`]). Each program writes then what NCBI wrote
/// before: the prologs, the reports of the batches before (without the epilog) and the batch's
/// reader messages. That batch is refused where a load whose order to the stop LOSAT does not know
/// writes a text, and where another subject's lookup throws on its id and may come before it (the
/// prefetch catches that one and ends its loads); where the batch has several stopping ids, its
/// loads end at the first lookup that NCBI makes ([`prefetch_stop`]). `before`: the subjects that
/// NCBI's BLAST database data loader loaded for the batches of the call before this one
/// ([`batch_loads`]), read only where the batch stops (the run's record of what its loader may
/// hold, [`TaxdbRun::note_loaded`]); a caller that records each batch itself gives none. The
/// outfmt 5 reports of the queries before the batch, for that stop, are the queries'
/// ([`PrefetchQueries::after_xml_reports`], [`XmlReportsBefore`]).
pub(crate) fn batch_prefetch<'s, 'b>(
    db: &SeqDb,
    run: &TaxdbRun,
    before: impl IntoIterator<Item = PrefetchSubject<'b>>,
    subjects: impl IntoIterator<Item = PrefetchSubject<'s>>,
    program: &str,
    queries: &PrefetchQueries<'_>,
) -> Result<crate::report::diag_log::Posts> {
    let subjects: Vec<PrefetchSubject<'s>> = subjects.into_iter().collect();
    if stops(&subjects) {
        run.keep_loaded(before);
        return Err(prefetch_stop(db, run, run, &subjects, program, queries));
    }
    let mut failed_lookup: Option<(u32, String)> = None;
    let mut failing = Vec::new();
    let mut fails = false;
    let mut uncertain: Option<(u32, String)> = None;
    let mut asked: std::collections::HashSet<u32> = Default::default();
    let mut loaded = Vec::new();
    for prefetched @ PrefetchSubject {
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
                failing.push(prefetched);
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
            db.oid(ordinal),
            &format!("{reason}, and has a Seq-align on the query's minus strand and its plus strand, whose length NCBI BLAST+ 2.17.0's prefetch asks for"),
            program,
        ));
    }
    // The failed lookup ends the loads: NCBI's warning of a UTF8String header that a load's
    // lookup reads comes before it where that load does, else not in the prefetch.
    let utf8_loads = db.scope_utf8_loads() && scope_utf8_among(&loaded);
    let posted = !(failed_lookup.is_some() && utf8_loads)
        || scope_utf8_before_failure(db, &loaded, &failing, program)?;
    let first_utf8 = loaded
        .iter()
        .find(|(_, parsed)| parsed.is_some_and(|subject| subject.scope_utf8_read().is_some()))
        .map(|&(ordinal, _)| ordinal);
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-790
    // A load that throws where a file of the taxonomy database does not open ends the prefetch's
    // loads, which catches it as a `CException` with the warning (line WQA D18b probe 12). The
    // loads write no text then: the taxonomy database never opens.
    let planned = prefetch_slot_plan(db, run, &loaded, failed_lookup.as_ref(), program);
    let reads = match planned.and_then(|planned| {
        let order = LoadsOrder::Prefetch(planned.as_ref());
        if posted {
            prefetch_loads(db, run, loaded, true, order, program)
        } else {
            db.without_scope_utf8_loads(|| prefetch_loads(db, run, loaded, true, order, program))
        }
    }) {
        Err(error)
            if queries.reports_unchanged_by_load_throws
                && TaxonomyRefusal::of_load_that_throws(&error) =>
        {
            // Which loads NCBI makes before the first that throws (a subject with a taxid, in
            // the order of the Seq-id handles) decides whether a lookup that reads a UTF8String
            // header posts NCBI's warning here (merge 6 of line RG item T5d).
            if let Some(ordinal) = first_utf8.filter(|_| utf8_loads) {
                return Err(parsed_refusal(
                    db.oid(ordinal),
                    "has a Seq-align id whose lookup reads a header with a UTF8String where a VisibleString is expected, in a prefetch whose loads throw where a file of the taxonomy database does not open: NCBI BLAST+ 2.17.0 writes its warning (803.10) there or not by the order of its loads, which LOSAT does not reproduce; make the file readable, or move it off the BLASTDB path",
                    program,
                ));
            }
            fails = true;
            Vec::new()
        }
        reads => reads?,
    };
    if let Some(stop) = &failed_lookup {
        // NCBI's warning of a UTF8String header before the failure is not a text of the
        // taxonomy database (`scope_utf8_before_failure`).
        if reads
            .iter()
            .any(|read| read.text.len() > read.load.as_ref().map_or(0, |load| load.lead))
        {
            return Err(failed_lookup_refusal(db, stop, program));
        }
    }
    let mut posts = crate::report::diag_log::Posts::new();
    posts.loads(reads);
    if failed_lookup.is_some() || fails {
        posts.lines(&prefetch_warning(program));
    }
    Ok(posts)
}

/// LOSAT's refusal of a batch's prefetch with a Seq-align id whose lookup fails (`stop`: the
/// subject's ordinal id and the id, [`batch_prefetch`]) where the other subjects' loads matter:
/// NCBI's loads end at the failed lookup, and which of them come before it LOSAT does not follow.
fn failed_lookup_refusal(db: &SeqDb, stop: &(u32, String), program: &str) -> anyhow::Error {
    let (ordinal, id) = stop;
    parsed_refusal(
        db.oid(*ordinal),
        &format!(
            "has the Seq-align id {id}, whose lookup in the database fails where NCBI BLAST+ 2.17.0's prefetch loads the batch's subjects, in an order that LOSAT does not reproduce: the taxonomy database's texts of the other subjects of the batch, which the prefetch may load before or not at all, are not known"
        ),
        program,
    )
}

/// The misses of a batch's prefetch in the slot of the offset error's taxid, by the ordinal id of
/// the reported subject ([`prefetch_slot_plan`]).
struct PrefetchSlots {
    misses: rustc_hash::FxHashMap<u32, Vec<bool>>,
}

impl PrefetchSlots {
    /// The misses of the load for the subject `ordinal`, whose lookups are `lookups` long: none
    /// for a load that looks no taxid of the slot up.
    fn misses(&self, ordinal: u32, lookups: usize) -> Vec<bool> {
        self.misses
            .get(&ordinal)
            .cloned()
            .unwrap_or_else(|| vec![false; lookups])
    }
}

/// The misses in the slot of the offset error's taxid of the loads of a batch's prefetch
/// ([`batch_prefetch`]), `loaded` its subjects in the batch's order (each ordinal once, as
/// [`load_read_texts`] reads them), on a database whose volume holds that taxid and another of
/// its slot ([`TaxdbRun::cache_collision`]); `None` on any other database, whose loads give the
/// same texts in any order. NCBI loads the batch's Seq-align ids in the order of its map of their
/// Seq-id handles (cited above [`LoadBeforeStop`]; [`load_places`]), ties by the addresses of
/// the handles' infos; so LOSAT follows the slot through the loads by their handles' keys and
/// takes the loads of one key as one step whose order only the boundaries between them decide
/// ([`VolumeLookups::prefetch_slots`]). LOSAT's refusal where that order changes the texts, or
/// leaves a slot that a later load's texts depend on (`SlotHolds::Either`). The reads of each
/// subject are [`load_read_texts`]'s, whose refusals come first. `stop`: a Seq-align id of the
/// batch whose lookup fails ([`batch_prefetch`]), where NCBI's loads end, so a load of the slot
/// that may come after it (one that only looks up the other taxid writes no text, but changes
/// what the slot holds) is refused ([`failed_lookup_refusal`]).
fn prefetch_slot_plan(
    db: &SeqDb,
    run: &TaxdbRun,
    loaded: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    stop: Option<&(u32, String)>,
    program: &str,
) -> Result<Option<PrefetchSlots>> {
    if !run.found(db) || db.search_path().is_none() || run.cache_collision(db).is_none() {
        return Ok(None);
    }
    let Some(last) = run.prefetch_texts(db).offset_error_taxid else {
        return Ok(None);
    };
    // NCBI's map of the batch's Seq-id handles ([`load_places`]): each load's key, and whether its
    // place among the loads of its key is not known (a batch of more than `PREFETCH_ONE_MAP_IDS`
    // distinct ids, or a Seq-align id of two subjects), so that LOSAT takes all of its volume's
    // loads as one key.
    let places = prefetch_places(db, loaded);
    let mut seen = std::collections::HashSet::new();
    let mut loads = Vec::new();
    let mut unordered = std::collections::HashSet::new();
    for (&(ordinal, parsed), place) in loaded.iter().zip(&places) {
        if !seen.insert(ordinal) {
            continue;
        }
        let blob = match parsed {
            Some(subject) => match prefetch_reads(db, ordinal, subject, true, program)? {
                Ok(Some(blob)) => blob,
                Ok(None) => continue,
                Err(reason) => return Err(parsed_refusal(db.oid(ordinal), &reason, program)),
            },
            None => ordinal,
        };
        let lookups = blob_taxids(db, ordinal, parsed, blob, program)?;
        if !lookups
            .iter()
            .any(|&taxid| tax_cache_slot(taxid) == tax_cache_slot(last))
        {
            continue;
        }
        let volume = db
            .volume_index(blob)
            .map_err(|error| database_error(error, program))?;
        if place.unknown.is_some() {
            unordered.insert(volume);
        }
        loads.push(SlotLoad {
            ordinal,
            blob,
            volume,
            lookups,
            order: place.order,
        });
    }
    if let (Some(stop), false) = (stop, loads.is_empty()) {
        return Err(failed_lookup_refusal(db, stop, program));
    }
    let misses = run
        .volume_lookups
        .lock()
        .expect("the volume lookups")
        .prefetch_slots(&loads, |volume| unordered.contains(&volume), last)
        .map_err(|(volume, cause)| {
            taxdb_error(slot_uncertain_error(db, volume, last, cause), program)
        })?;
    Ok(Some(PrefetchSlots { misses }))
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3279-3298 (`x_GetBioseqHandlesSorted`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
// The prefetch's loads go in the order of the map of their Seq-id handles; the lookup that throws
// on its id (`GetRecords` -> `x_GetOid`) ends them, and the prefetch catches it with its warning
// (blast_app_util.cpp:778-790): a load after it is not made, so its lookup reads no header
// there. Probe `T/t5d/fix5/probe` (NCBI 2.17.0 n1, 5 programs): `gnl|BL_ORD_ID|x002` (packed 2)
// before `gb|DUP00005.1|` (5) writes `Error pre-fetching sequence data` and no warning in a batch
// whose reports read nothing of the subjects; the warning comes at the next batch's prefetch.
/// Where a load of a batch's prefetch whose lookup throws on its id (`failing`, which ends the
/// loads, [`batch_prefetch`]) and one whose lookup reads a header with a UTF8String where a
/// VisibleString is expected (`ParsedSubject::scope_utf8`) are among its loads `loaded`: whether
/// the first such load in NCBI's order of the Seq-id handles ([`load_places`]) comes before every
/// failing lookup (`true`: the prefetch writes NCBI's warning there, before its own,
/// [`prefetch_loads`]) or after one (`false`: NCBI makes none of these loads, so the prefetch
/// writes no warning, which a later load or read of the run may write; the reports of such a
/// batch that read a subject read the one whose lookup throws, which LOSAT refuses,
/// `ParsedSubject::fetched`). LOSAT's refusal where their order comes from the addresses of
/// NCBI's Seq-id infos or is not known.
#[cold]
#[inline(never)]
fn scope_utf8_before_failure(
    db: &SeqDb,
    loaded: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    failing: &[PrefetchSubject<'_>],
    program: &str,
) -> Result<bool> {
    let subjects: Vec<PrefetchSubject<'_>> = loaded
        .iter()
        .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed))
        .chain(failing.iter().copied())
        .collect();
    let (places, _) = load_places(db, &subjects);
    let (loads, fails) = places.split_at(loaded.len());
    let Some((index, at)) = loads
        .iter()
        .enumerate()
        .filter(|(index, _)| {
            loaded[*index]
                .1
                .is_some_and(|subject| subject.scope_utf8_read().is_some())
        })
        .min_by_key(|(_, place)| place.order)
    else {
        return Ok(true);
    };
    let (ordinal, parsed) = loaded[index];
    let other = parsed
        .and_then(|subject| subject.scope_utf8_read())
        .expect("a subject whose lookup reads such a header");
    let failure = fails
        .iter()
        .zip(failing)
        .min_by_key(|(place, _)| place.order)
        .expect("a failing lookup");
    match at.order_unknown(*failure.0) {
        None => Ok(at.order < failure.0.order),
        Some(why) => Err(parsed_refusal(
            db.oid(ordinal),
            &format!(
                "has a Seq-align id by which NCBI BLAST+'s lookup gives database sequence {other} first, whose header holds a UTF8String where a VisibleString is expected: NCBI BLAST+ 2.17.0's prefetch writes its warning (803.10) at the load of that id where it loads that id before the lookup of {}, which throws and ends its loads, and no warning where it loads it after ({why})",
                String::from_utf8_lossy(
                    &failure
                        .1
                        .parsed
                        .expect("a failing lookup has parsed Seq-ids")
                        .alignment_fasta()
                ),
            ),
            program,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:778-786
/// NCBI's warning where the prefetch of a query batch fails ([`batch_prefetch`], and BLASTX's
/// Seq-align of empty subject rows, `blastx::empty_rows::prefetch_reads_empty_row`): the one
/// producer of its text.
pub(crate) fn prefetch_warning(program: &str) -> Vec<u8> {
    // A post of NCBI's main thread (`ncbi_exception::post_line`: the header alone after the run's diagnostic
    // stream failed).
    crate::ncbi_exception::post_line(
        crate::ncbi_exception::Severity::Warning,
        None,
        Some(&program.to_ascii_lowercase()),
        b"Error pre-fetching sequence data",
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:728-753
/// Whether NCBI's prefetch of a query batch reads the batch's Seq-aligns (their subject ids and
/// ranges, `s_ExtractSeqidsAndRanges`) for the first format `outfmt`: every format but the
/// archive (outfmt 11, which the apps write without calling the prefetch) and, unless
/// `PRE_FETCH_SEQS_LIMIT=2147483647`, the Seq-annots 8, 9 and 12 (`s_IsPrefetchFormat`, cited at
/// `blastx::native::reported_batch`); a batch whose prefetch is off ([`prefetch_off`]) reads
/// nothing. Wider than [`prefetches_subjects`], which tells where the loads write the taxonomy
/// database's texts.
pub(crate) fn prefetch_reads_seq_aligns(outfmt: i32) -> bool {
    match outfmt {
        11 => false,
        8 | 9 | 12 => prefetch_seqs_limit() == Some(i32::MAX),
        _ => true,
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:781-784 (cited at
// `ParsedSubject::prefetched`)
/// The stop of NCBI's prefetch of one query batch's reported subjects `subjects` for a first
/// format that writes the taxonomy database's texts at each `Hit` (`TaxdbPrinter`: the XML2 and
/// JSON formats, outfmt 13 to 16, which NCBI's `s_IsPrefetchFormat`, cited at
/// [`prefetch_messages`], prefetches for and [`prefetches_subjects`] leaves out) and a batch whose
/// prefetch is on (`PRE_FETCH_SEQS_LIMIT`, [`prefetch_off`]): `Err` as in [`batch_prefetch`] where
/// a subject's lookup fails on a file with NCBI's exception, else `Ok` (the batch's other loads
/// and their texts are the formats' own, at each `Hit`). `run` records the lookups of the reports
/// of the batches before ([`read_reported_taxonomy`]: the loads of NCBI's prefetch and reports of
/// those batches, and the open of the taxonomy database at their first taxid) and nothing of this
/// batch: BLASTX's own run, which writes each batch's reports before it searches the next; the
/// programs that write those reports after the stop ask [`held_batch_prefetch_stop`]. As in
/// [`batch_prefetch`], the texts of the loads before the stop after those lookups come before the
/// error, and the batch is refused where a load whose order LOSAT does not know writes a text.
/// `run` also holds the subjects that NCBI's data loader loaded for the batches before
/// ([`TaxdbRun::note_loaded`], which BLASTX records as it writes each batch). A batch that does not
/// stop: `Ok` with NCBI's warning of a UTF8String header that a load's lookup reads
/// ([`scope_utf8_prefetch_post`]; empty where none), which the caller writes where the batch's
/// prefetch comes.
pub(crate) fn batch_prefetch_stop<'s>(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: impl IntoIterator<Item = PrefetchSubject<'s>>,
    program: &str,
    queries: &PrefetchQueries<'_>,
) -> Result<Vec<u8>> {
    let subjects: Vec<PrefetchSubject<'s>> = subjects.into_iter().collect();
    if stops(&subjects) {
        Err(prefetch_stop(db, run, run, &subjects, program, queries))
    } else {
        scope_utf8_prefetch_post(db, &subjects, program)
    }
}

/// NCBI's warning (803.10) where a load of the prefetch of a batch that does not stop, of the XML2
/// and JSON formats (whose loads' taxonomy texts are the formats' own, at each `Hit`,
/// [`batch_prefetch_stop`]), has a lookup that reads a header with a UTF8String where a
/// VisibleString is expected (`ParsedSubject::scope_utf8`), and the run has not posted it: the
/// first such load in NCBI's order of the Seq-id handles ([`load_places`]) posts it
/// (`load_read_texts`), the only text of the prefetch. Refused where two such loads whose
/// warnings differ have an order that comes from the addresses of NCBI's Seq-id infos or is not
/// known. Empty where no load posts it.
pub(crate) fn scope_utf8_prefetch_post(
    db: &SeqDb,
    subjects: &[PrefetchSubject<'_>],
    program: &str,
) -> Result<Vec<u8>> {
    if !db.scope_utf8_loads() {
        return Ok(Vec::new());
    }
    let scoped: Vec<(u32, Option<&crate::report::db_subject::ParsedSubject>)> = subjects
        .iter()
        .filter(|subject| !subject.parsed.is_some_and(|parsed| parsed.prefetch_fails()))
        .map(|subject| (subject.ordinal, subject.parsed))
        .collect();
    if !scope_utf8_among(&scoped) {
        return Ok(Vec::new());
    }
    let places = prefetch_places(db, &scoped);
    let path = |index: usize| scope_utf8_path(db, scoped[index].1);
    let mut candidates: Vec<usize> = (0..scoped.len())
        .filter(|&index| path(index).is_some())
        .collect();
    candidates.sort_by_key(|&index| places[index].order);
    let Some(&first) = candidates.first() else {
        return Ok(Vec::new());
    };
    for &other in &candidates[1..] {
        let unknown = places[other].order_unknown(places[first]);
        if let Some(why) = unknown.filter(|_| path(other) != path(first)) {
            return Err(parsed_refusal(
                db.oid(scoped[first].0),
                &format!(
                    "has a Seq-align id whose lookup reads a header with a UTF8String where a VisibleString is expected, as the Seq-align id of database sequence {} does another one, and NCBI BLAST+ 2.17.0's prefetch writes the warning (803.10) of the first of their loads, whose order LOSAT does not know ({why})",
                    db.oid(scoped[other].0)
                ),
                program,
            ));
        }
    }
    Ok(scope_utf8_lead(
        db,
        scoped[first].1,
        &program.to_ascii_lowercase(),
    ))
}

/// [`batch_prefetch_stop`] where the XML2 and JSON reports of the batches before are written after
/// the stop (BLASTP, BLASTN, TBLASTN and TBLASTX hold those documents to the end of the run, and
/// write the batches before a failed batch there): their lookups are recorded where the batch
/// stops, in a run of their own ([`reports_record`]) of the subjects `before` that those batches
/// report ([`PrefetchSubject`], in any order). Nothing is read where the batch does not stop.
/// `run`, the run's own record, takes those subjects as sequences that the run loaded before the
/// batch ([`TaxdbRun::note_loaded`]: the reports of these formats read every subject they show),
/// beside those of the calls before, and the other loads of the batches before (`loaded`, read
/// only where the batch stops: psiblast's traceback beyond the cut, which no report reads), all
/// without the record's gate, as every stop records the call's loads before it (`keep_loaded`).
pub(crate) fn held_batch_prefetch_stop<'s, 'b, 'l>(
    db: &SeqDb,
    run: &TaxdbRun,
    before: impl IntoIterator<Item = PrefetchSubject<'b>>,
    loaded: impl IntoIterator<Item = PrefetchSubject<'l>>,
    subjects: impl IntoIterator<Item = PrefetchSubject<'s>>,
    program: &str,
    queries: &PrefetchQueries<'_>,
) -> Result<Vec<u8>> {
    let subjects: Vec<PrefetchSubject<'s>> = subjects.into_iter().collect();
    if !stops(&subjects) {
        return scope_utf8_prefetch_post(db, &subjects, program);
    }
    let before: Vec<PrefetchSubject<'b>> = before.into_iter().collect();
    run.keep_loaded(before.iter().copied());
    run.keep_loaded(loaded);
    let lookups = reports_record(db, before, program)?;
    Err(prefetch_stop(
        db, &lookups, run, &subjects, program, queries,
    ))
}

/// The record of the taxonomy lookups of the XML2 and JSON reports of the subjects `reported`
/// (those of the batches before a stopping batch, [`held_batch_prefetch_stop`]):
/// [`read_reported_taxonomy`] of each once, in a run of its own that has looked nothing up before;
/// the texts are dropped (the reports write their own). Its refusals and errors are those that the
/// reports of those subjects make.
fn reports_record<'s>(
    db: &SeqDb,
    reported: impl IntoIterator<Item = PrefetchSubject<'s>>,
    program: &str,
) -> Result<TaxdbRun> {
    let mut seen = std::collections::HashSet::new();
    let (ordinals, parsed): (
        Vec<u32>,
        Vec<Option<&crate::report::db_subject::ParsedSubject>>,
    ) = reported
        .into_iter()
        .filter(|subject| seen.insert(subject.ordinal))
        .map(|subject| (subject.ordinal, subject.parsed))
        .unzip();
    let ranks = (0..ordinals.len()).collect();
    let run = TaxdbRun::default();
    read_reported_taxonomy(db, &run, &ordinals, &ranks, &parsed, program)?;
    Ok(run)
}

/// Whether one of `subjects` has a lookup that fails on a file with NCBI's exception
/// (`ParsedSubject::prefetch_stops`): the batch's prefetch stops ([`prefetch_stop`]).
fn stops(subjects: &[PrefetchSubject<'_>]) -> bool {
    subjects.iter().any(|subject| {
        subject
            .parsed
            .is_some_and(|parsed| parsed.prefetch_stops().is_some())
    })
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3346-3366
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2014-2026
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:389-396
// The prefetch's `GetBioseqHandles` (`LoadSequencesToScope`, cited at `check_prefetched_parsed`)
// loads the distinct Seq-align ids of a batch of at most 200 of them in one map of Seq-id handles
// (`x_GetBioseqHandlesSorted`, scope_impl.cpp:3280-3299; `CDataLoader::GetBlobs`, cited above
// `PrefetchSubject`): by `m_Packed` less 1 as an unsigned number (`seqid::mask_handle`), then by the
// address of the handle's info. Each load asks the BLAST database for the id's OID first, so the
// first id in that order whose lookup fails with NCBI's exception ends the loads, and the loads
// before it have written their texts (probe of line DBID loadord, `R/loadord/logs/probe.tsv`:
// `gnl|BL_ORD_ID|1` and `|2` load before `ref|NP_000003.1|`, `|4` to `|7` never, `|3` before it
// with one command line and after it with another (`-outfmt 6` and `0`); every packed id before
// `pdb|1ABC|A`, `sp|P00007.1|NAME_X` and `lcl|abc`). More ids are first sorted by
// `CSortableSeq_id` (seq_id_sort.cpp: a `general` id by its database and tag parts, the others by
// `CSeq_id_Handle::PLessOrdered`) and loaded in pieces of 100, an order that LOSAT does not
// model: of the probe's 201 ids, the first 100 `gnl|BL_ORD_ID|n` load before `dbj|AB000001.1|`
// (a `ddbj` id sorts after every `general` one), where 200 ids load it first.
/// Whether NCBI's prefetch of a batch that stops ([`prefetch_stop`]) loads one of the batch's
/// other subjects before the first lookup that stops it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LoadBeforeStop {
    Yes,
    No,
    /// LOSAT does not know, with the reason.
    Unknown(&'static str),
}

/// The most distinct Seq-align ids that NCBI's prefetch loads in one map of Seq-id handles
/// (`CScope_Impl::GetBioseqHandles`, cited at [`load_before_stop`]).
const PREFETCH_ONE_MAP_IDS: usize = 200;

/// The key of the Seq-id handle of the Seq-align id of `subject` in NCBI's map of a batch's loads:
/// its position among handles of different `m_Packed` (`seqid::mask_handle::mask_handle_order`)
/// and its identity (`mask_handle_key`; `None` where no pure function knows it: a new handle for
/// each Seq-id, or a giim id). A sequence without parsed Seq-ids is `gnl|BL_ORD_ID|<OID>`.
fn prefetch_load_key(db: &SeqDb, subject: &PrefetchSubject<'_>) -> (u64, Option<Vec<u8>>) {
    use crate::seqid::mask_handle::{mask_handle_key, mask_handle_order, mask_handle_packed};
    let ordinal_id;
    let id = match subject.parsed {
        Some(parsed) => parsed.alignment_id(),
        None => {
            ordinal_id = crate::seqid::ParsedSeqId::General(crate::seqid::model::Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: crate::seqid::model::ObjectId::Id(db.oid(subject.ordinal).get() as i32),
            });
            &ordinal_id
        }
    };
    (
        mask_handle_order(mask_handle_packed(id)),
        mask_handle_key(id).ok().flatten(),
    )
}

/// [`LoadBeforeStop`] of each of `subjects` (a stopping batch's reported sequences, in any order,
/// a sequence possibly more than once) in NCBI's order of the loads (cited above
/// [`LoadBeforeStop`]), and the order of the lookups that may stop it ([`LoadOrder`]): the
/// subjects whose lookups stop the batch (`ParsedSubject::prefetch_stops`) load nothing; the
/// others load by their handle's place to the stops ([`LoadOrder::before`]). Unknown besides: a
/// Seq-align id that another subject has too, and every subject of a batch of more than
/// [`PREFETCH_ONE_MAP_IDS`] distinct Seq-align ids (whose stops' order LOSAT does not know either).
/// The one owner of NCBI's order of a stopping batch's loads and lookups.
fn load_before_stop(db: &SeqDb, subjects: &[PrefetchSubject<'_>]) -> LoadOrder {
    let (places, many) = load_places(db, subjects);
    // The stopping subjects, each ordinal once, by their index in `subjects`.
    let mut ordinals = std::collections::HashSet::new();
    let mut stopping: Vec<usize> = (0..subjects.len())
        .filter(|&index| ordinals.insert(subjects[index].ordinal) && places[index].stops)
        .collect();
    assert!(
        !stopping.is_empty(),
        "a stopping batch has a stopping subject"
    );
    // A sort that keeps the order of the batch among the handles of one key (ties).
    stopping.sort_by_key(|&index| if many { 0 } else { places[index].order });
    LoadOrder {
        places,
        stops: stopping,
        many,
    }
}

/// The place of each of a batch's reported sequences `subjects` (in any order, a sequence possibly
/// more than once) in NCBI's map of the prefetch's Seq-id handles (cited above
/// [`LoadBeforeStop`]): its handle's key ([`prefetch_load_key`]), whether its lookup stops the
/// batch (`ParsedSubject::prefetch_stops`), and why its place is not known: a Seq-align id that
/// another subject has too, and every subject of a batch of more than [`PREFETCH_ONE_MAP_IDS`]
/// distinct Seq-align ids (the second value). The one owner of that map for the stops
/// ([`load_before_stop`]) and for the place of NCBI's warning of a UTF8String header among the
/// loads ([`scope_utf8_prefetch_loads`]).
fn load_places(db: &SeqDb, subjects: &[PrefetchSubject<'_>]) -> (Vec<LoadPlace>, bool) {
    let stops = |subject: &PrefetchSubject<'_>| {
        subject
            .parsed
            .is_some_and(|parsed| parsed.prefetch_stops().is_some())
    };
    let keys: Vec<(u64, Option<Vec<u8>>)> = subjects
        .iter()
        .map(|subject| prefetch_load_key(db, subject))
        .collect();
    let mut ordinals = std::collections::HashSet::new();
    let mut handles: std::collections::HashMap<&[u8], u32> = Default::default();
    let mut shared = std::collections::HashSet::new();
    let mut ids = 0;
    for (subject, (_, handle)) in subjects.iter().zip(&keys) {
        if !ordinals.insert(subject.ordinal) {
            continue;
        }
        match handle {
            Some(handle) => match handles.entry(handle) {
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(subject.ordinal);
                    ids += 1;
                }
                std::collections::hash_map::Entry::Occupied(entry) => {
                    shared.insert(*entry.get());
                    shared.insert(subject.ordinal);
                }
            },
            None => ids += 1,
        }
    }
    let many = ids > PREFETCH_ONE_MAP_IDS;
    let places = subjects
        .iter()
        .zip(&keys)
        .map(|(subject, (order, handle))| LoadPlace {
            order: *order,
            handle: handle.as_ref().map(|handle| {
                use std::hash::{Hash, Hasher};
                let mut hasher = std::collections::hash_map::DefaultHasher::new();
                handle.hash(&mut hasher);
                hasher.finish()
            }),
            many,
            stops: stops(subject),
            unknown: if many {
                Some("the batch has more than 200 Seq-align ids, which NCBI's prefetch sorts by CSortableSeq_id and loads in pieces of 100")
            } else if shared.contains(&subject.ordinal) {
                Some("another subject of the batch has the same Seq-align id")
            } else {
                None
            },
        })
        .collect::<Vec<_>>();
    (places, many)
}

/// NCBI's order of the loads and lookups of a stopping batch's subjects ([`load_before_stop`]).
struct LoadOrder {
    /// The place of each subject's load, in the order of the batch's subjects.
    places: Vec<LoadPlace>,
    /// The stopping subjects (their indexes in the batch's subjects, each sequence once) in the
    /// order of NCBI's lookups of their Seq-align ids: by their handles' keys
    /// ([`LoadOrder::stop_key`]), a group of one key in the batch's order (NCBI orders a group by
    /// the addresses of the handles' infos).
    stops: Vec<usize>,
    /// The batch has more than [`PREFETCH_ONE_MAP_IDS`] distinct Seq-align ids.
    many: bool,
}

/// The place of one subject's load in NCBI's map of a stopping batch's Seq-id handles
/// ([`LoadOrder`]).
#[derive(Clone, Copy, Debug)]
struct LoadPlace {
    /// Its handle's key among handles of different `m_Packed` ([`prefetch_load_key`]).
    order: u64,
    /// Its lookup stops the batch (`ParsedSubject::prefetch_stops`): it loads nothing.
    stops: bool,
    /// Why LOSAT does not know its place, whatever the stops' keys.
    unknown: Option<&'static str>,
    /// Its handle's identity (a hash of [`prefetch_load_key`]'s; `None` where no pure function
    /// knows it): two subjects of one handle are one load of NCBI's map.
    handle: Option<u64>,
    /// The batch has more than [`PREFETCH_ONE_MAP_IDS`] distinct Seq-align ids.
    many: bool,
}

impl LoadPlace {
    /// Why NCBI's order of this load and the load `other` is not known, for the place of NCBI's
    /// warning of a UTF8String header among a prefetch's loads ([`scope_utf8_prefetch_loads`]):
    /// a batch of more than [`PREFETCH_ONE_MAP_IDS`] ids, or two handles of one key (one packed
    /// number, or neither packed), which NCBI's map orders by the address of their info
    /// ([`ADDRESS_ORDER`]); `None` where the keys order them or they are one handle (one load).
    fn order_unknown(self, other: LoadPlace) -> Option<&'static str> {
        if self.many {
            return Some("the batch has more than 200 Seq-align ids, which NCBI's prefetch sorts by CSortableSeq_id and loads in pieces of 100");
        }
        if self.order != other.order || (self.handle.is_some() && self.handle == other.handle) {
            return None;
        }
        Some(ADDRESS_ORDER)
    }
}

impl LoadOrder {
    /// The key of the stopping subject `index` (of the batch's subjects) in the order of the
    /// lookups: its handle's, or one key for all where the batch has more than
    /// [`PREFETCH_ONE_MAP_IDS`] ids, whose order LOSAT does not know.
    fn stop_key(&self, index: usize) -> u64 {
        if self.many {
            0
        } else {
            self.places[index].order
        }
    }

    /// The key of the first stopping lookup in NCBI's order.
    fn first_stop_key(&self) -> u64 {
        self.stop_key(self.stops[0])
    }

    /// [`LoadOrder::before`] where the loads end at the first stopping lookup.
    #[cfg(test)]
    fn at_first_stop(&self) -> Vec<LoadBeforeStop> {
        let first = self.first_stop_key();
        self.before(first, first)
    }

    /// [`LoadBeforeStop`] of each subject (in the order of the batch's subjects) where NCBI's loads
    /// end at a stopping lookup whose key lies between `first` (the first stopping key) and `last`
    /// (where they end for sure, [`prefetch_stop`]): `No` for the stopping subjects and the
    /// handles after `last`, `Yes` for those before `first`; `Unknown` for a handle of the key of
    /// `first` or of `last` (both unpacked included), whose order NCBI takes from its info's
    /// address, for one between them (loaded where the lookups before it are not made), and for
    /// the subjects whose place is not known ([`LoadPlace::unknown`]).
    fn before(&self, first: u64, last: u64) -> Vec<LoadBeforeStop> {
        self.places
            .iter()
            .map(|place| {
                if place.stops {
                    LoadBeforeStop::No
                } else if let Some(why) = place.unknown {
                    LoadBeforeStop::Unknown(why)
                } else if place.order < first {
                    LoadBeforeStop::Yes
                } else if place.order > last {
                    LoadBeforeStop::No
                } else if place.order == first || place.order == last {
                    LoadBeforeStop::Unknown("its Seq-id handle has the same packed number as a stopping one's (or neither is packed), which NCBI orders by their address")
                } else {
                    LoadBeforeStop::Unknown("its Seq-id handle comes after the first stopping one's and before the one whose lookup NCBI surely makes: NCBI loads it only where it makes none of the stopping lookups before it, and LOSAT does not know whether it makes them")
                }
            })
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-387
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:98-110
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3288-3297 (cited at
// `load_before_stop`): `GetBioseqHandle` of `ids[i]` with `CScope::eGetBioseq_Resolved` first, then the
// data sources by priority, the queries' entries (`AddTopLevelSeqEntry`) before the BLAST
// database's loader.
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2032-2068 and 1449-1468
// A data source without a loader (the queries' entries) finds a Bioseq of the id itself or, where
// it holds none, of an id that the id matches weakly (`x_GetSeqMatch`), as the reports' reads
// find a query (`refuse_id_of_a_query`'s keys, `QueryScopeKeys::own` and `weak`).
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2310-2313
/// Why NCBI's prefetch may never make the lookup of the Seq-align id of the stopping subject
/// `ordinal` (`parsed`) of the batch of `queries` that fails with its exception (reviews 1 to 3 of
/// line DBID loadord, line DBID DECISIONS 60 (5) to (10)), and what NCBI then loads instead:
/// NCBI's scope finds a query by its id before it asks the BLAST database's data loader, and that
/// loader keeps every Seq-id of each sequence that the run has loaded (`RegisterIds`) and gives
/// that sequence for the id without a lookup. So the lookup may not be made where
/// - a query that the scope holds at the batch's prefetch ([`PrefetchQueries::new`]'s scope, made
///   here only) has an id that the scope takes for the stopping id, or that the stopping id
///   matches weakly (the keys of [`refuse_id_of_a_query`]): nothing is loaded instead;
/// - a load of the batch that may come before the stop (`before` of `subjects`, from
///   [`load_before_stop`]: `Yes` or `Unknown`) gives a sequence (the first OID of its Seq-align
///   id's lookup, `x_GetOid`; an unparsed subject's own) that holds an id of the stopping id's
///   Seq-id handle;
/// - another sequence of the database that holds such an id may have been loaded before the batch:
///   the run's record of its loads holds it (`record`, [`TaxdbRun::note_loaded`]: the prefetches,
///   the reads of the reports where a prefetch is off and the PSSMs' reads of the batches,
///   queries and rounds before) and its own Seq-align id's lookup does not fail (a sequence whose
///   own lookup fails loads nothing, and a copy with the same Seq-align id is the same lookup,
///   which would have stopped the run there), or it is the first sequence that the lookup of an id
///   of a sequence loaded before gives (NCBI's reads by another id of that sequence, a superset);
///   and, after calls whose loads nothing recorded ([`PrefetchQueries::after_unrecorded_calls`]),
///   any such sequence whose own lookup does not fail (found in the LMDB files that have the
///   accession table, [`crate::seqdb::lookup::SeqidLookup::in_indexed_files`]).
///
/// The sequences that may then be given instead are each of those that hold such an id
/// ([`StopUnsure::instead`]), loaded by the stopping id where the batch's loads go on to a later
/// stop ([`prefetch_stop`]). `None`: NCBI looks it up.
fn stop_lookup_unsure(
    db: &SeqDb,
    record: &TaxdbRun,
    subjects: &[PrefetchSubject<'_>],
    before: &[LoadBeforeStop],
    queries: &PrefetchQueries<'_>,
    ordinal: u32,
    parsed: &crate::report::db_subject::ParsedSubject,
) -> Result<Option<StopUnsure>> {
    let id = parsed.alignment_id();
    let keys = parsed.query_scope_keys();
    let scope = (queries.scope)();
    if scope.has_key(&keys.own) || scope.has_key(&keys.weak) {
        return Ok(Some(StopUnsure {
            why: "NCBI's scope holds a query (of the batch, or in outfmt 5 of a batch before) with an id that it takes for this one, or that this one matches weakly, at the batch's prefetch, where it finds that query before it asks the database".to_owned(),
            instead: Some(Vec::new()),
        }));
    }
    let handle = crate::seqid::mask_handle::mask_handle_key(id)
        .ok()
        .flatten();
    let holds = |oid: u32| {
        handle
            .as_deref()
            .is_none_or(|handle| holds_handle(db, db.ordinal_of(oid), handle))
    };
    let lookup = db.seqid_lookup()?;
    // The first reason, and every sequence that may be given instead (by OID).
    let mut why: Option<String> = None;
    let mut instead: Option<Vec<u32>> = Some(Vec::new());
    let mut give = |oid: u32, reason: &dyn Fn() -> String| {
        why.get_or_insert_with(reason);
        if let Some(instead) = instead.as_mut() {
            instead.push(oid);
        }
    };
    for (subject, before) in subjects.iter().zip(before) {
        if *before == LoadBeforeStop::No {
            continue;
        }
        let given = match subject.parsed {
            Some(parsed) => lookup.seqid_to_oid(parsed.alignment_id()).ok().flatten(),
            None => Some(db.oid(subject.ordinal).get()),
        };
        if let Some(oid) = given.filter(|&oid| holds(oid)) {
            give(oid, &|| {
                format!(
                "a load of the batch's prefetch that may come before it gives a sequence of the database (OID {oid}) that holds an id that NCBI's Seq-id handles take for it, which registers the id in NCBI's BLAST database data loader, which then gives that sequence for it without the lookup"
            )
            });
        }
    }
    let loaded = record.loaded();
    if loaded.is_empty() && !queries.unrecorded_calls {
        return Ok(why.map(|why| StopUnsure {
            why,
            instead: instead.map(|oids| ordinals_of(db, oids)),
        }));
    }
    let own = db.oid(ordinal).get();
    // The other sequences of the database that hold such an id, found in the LMDB files that have
    // the accession table and kept where they hold an id of the handle; the lookup's error where
    // it fails (a version 4 database's ISAM that does not open): each sequence is then asked.
    let indexed = db.seqid_lookup()?.in_indexed_files();
    let holders = indexed
        .seqid_to_oids(id, true)
        .and_then(|mut holders| {
            holders.extend(indexed.case_folded_oids(id)?);
            holders.extend(indexed.case_variant_oids(id)?);
            Ok(holders)
        })
        .map(|mut holders| {
            holders.retain(|&oid| oid != own && holds(oid));
            holders.sort_unstable();
            holders.dedup();
            holders
        });
    let holder = |oid: u32| {
        oid != own
            && match &holders {
                Ok(holders) => holders.binary_search(&oid).is_ok(),
                Err(_) => holds(oid),
            }
    };
    let ids_of = |oid_ordinal: u32| match db.read_header(oid_ordinal) {
        Ok(crate::seqdb::SequenceHeader::Parsed(set)) => Some(set),
        _ => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:389-396
    // (`SeqidToOid`, cited at `SeqidLookup::seqid_to_oid`, gives the first OID of the lookup,
    // which every lookup of NCBI's loads here is.)
    // A sequence of a call before whose reports LOSAT has no record of: its own Seq-align id's
    // lookup does not fail (a load of it by that id may be made).
    let own_lookup_ok = |oid: u32| {
        ids_of(db.ordinal_of(oid)).is_none_or(|set| {
            set.alignment_id()
                .is_none_or(|id| lookup.seqid_to_oid(id).is_ok())
        })
    };
    // A load of the recorded sequence `oid` by its own Seq-align id loads `oid` itself: the
    // loader loads the first OID that the id's lookup gives, which may be another sequence's
    // (review 4 of line DBID loadord, L3: two holders of one gi); that one is asked below, as the
    // first sequence that a lookup of an id of a sequence loaded before gives.
    let loads_itself = |oid: u32| {
        ids_of(db.ordinal_of(oid)).is_none_or(|set| {
            set.alignment_id().is_none_or(|id| {
                lookup
                    .seqid_to_oid(id)
                    .is_ok_and(|first| first == Some(oid))
            })
        })
    };
    if !holders.as_ref().is_ok_and(Vec::is_empty) {
        if queries.unrecorded_calls {
            match &holders {
                Ok(holders) => {
                    for &oid in holders.iter().filter(|&&oid| own_lookup_ok(oid)) {
                        give(oid, &|| {
                            format!("another sequence of the database (OID {oid}) holds an id that NCBI's Seq-id handles take for it, and a load of that sequence by its own Seq-align id in a call of the run before this one, whose reports LOSAT has no record of here, registers the id in NCBI's BLAST database data loader, which then gives that sequence for it without the lookup")
                        });
                    }
                }
                Err(error) => {
                    why.get_or_insert_with(|| format!(
                        "a call of the run before this one, whose reports LOSAT has no record of here, may have loaded another sequence of the database that holds an id that NCBI's Seq-id handles take for it, and LOSAT cannot look for one ({error})"
                    ));
                    instead = None;
                }
            }
        } else {
            let registered = |oid: u32| {
                format!("another sequence of the database (OID {oid}) holds an id that NCBI's Seq-id handles take for it, and the run loaded that sequence, or one that shares an id with it, before the batch: its load there registers the id in NCBI's BLAST database data loader, which then gives that sequence for it without the lookup")
            };
            for &before in &loaded {
                let oid = db.oid(before).get();
                if holder(oid) && loads_itself(oid) {
                    give(oid, &|| registered(oid));
                }
            }
            // The first sequence that a lookup of any id of a sequence loaded before gives.
            for &before in &loaded {
                let Some(set) = ids_of(before) else {
                    continue;
                };
                for id in set.deflines.iter().flat_map(|defline| &defline.ids) {
                    let first = lookup.seqid_to_oid(id).ok().flatten();
                    if let Some(oid) = first.filter(|&oid| holder(oid)) {
                        give(oid, &|| registered(oid));
                    }
                }
            }
        }
    }
    Ok(why.map(|why| StopUnsure {
        why,
        instead: instead.map(|oids| ordinals_of(db, oids)),
    }))
}

/// The ordinal ids of the OIDs `oids`, each once, in their order.
fn ordinals_of(db: &SeqDb, mut oids: Vec<u32>) -> Vec<u32> {
    oids.sort_unstable();
    oids.dedup();
    oids.into_iter().map(|oid| db.ordinal_of(oid)).collect()
}

/// Why NCBI's prefetch may not make a stopping lookup ([`stop_lookup_unsure`]).
#[derive(Debug, PartialEq, Eq)]
struct StopUnsure {
    /// The first reason.
    why: String,
    /// The sequences (ordinal ids) that NCBI's BLAST database data loader may give for the
    /// stopping id without its lookup (empty where only the scope finds it, which loads nothing);
    /// `None` where LOSAT cannot tell which.
    instead: Option<Vec<u32>>,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:767-777
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-396 (quoted
// in parts at `stop_lookup_unsure` and above `LoadBeforeStop`): `x_GetOid` of each handle of the
// map in its order, `m_Ids.find` first, then `SeqidToOid`, whose exception ends the loads; the
// handles after it are never looked up.
/// The end of a batch's prefetch whose subjects `subjects` hold one or more whose lookups fail on
/// a file with NCBI's exception ([`batch_prefetch`], [`batch_prefetch_stop`]), after the lookups
/// that `run` records and the loads that `record` records ([`TaxdbRun::note_loaded`]; the same
/// run, or the run's own where `run` is the record of the held reports,
/// [`held_batch_prefetch_stop`]): the [`HeaderReadError`] of `queries` (its `error_query`) with
/// the texts of the loads before the stop in NCBI's order ([`load_before_stop`]; with the texts
/// that loads NCBI did not make left to them, `report::query_warnings::ReadLoad::place`: those
/// of an outfmt 5 report before the batch that NCBI's XML stream cut, also where the program
/// writes that report after the stop, [`PrefetchQueries::after_xml_reports`], line RG K4h) and
/// NCBI's `BLAST Database error`, or LOSAT's refusal where what NCBI writes depends on what LOSAT
/// does not know.
///
/// NCBI's loads end at the first stopping lookup that it makes, in the order of the handles
/// (`LoadSequencesToScope` rethrows the `CSeqDBException`, cited at `check_prefetched_parsed`):
/// the stops are asked in that order ([`LoadOrder::stops`]) until one that NCBI surely makes
/// ([`stop_lookup_unsure`]: `None`), with the others of its key. A stop before it that NCBI may not
/// make (review 3 of line DBID loadord, H1) does not refuse the batch by itself: the loads from
/// the first stop's key to the sure one's and the sequences that such a skipped lookup is given
/// instead ([`StopUnsure::instead`]) are loads whose place LOSAT does not know, which refuse the
/// batch where one of them writes a text (as a load of the first stop's packed number), and the
/// stops that NCBI may end at must have one message. Refused: no stop that NCBI surely makes (it
/// may then go on), a sequence given instead that LOSAT cannot name.
fn prefetch_stop(
    db: &SeqDb,
    run: &TaxdbRun,
    record: &TaxdbRun,
    subjects: &[PrefetchSubject<'_>],
    program: &str,
    queries: &PrefetchQueries<'_>,
) -> anyhow::Error {
    let stopping = |index: usize| {
        let subject = &subjects[index];
        let parsed = subject
            .parsed
            .expect("a stopping subject has parsed Seq-ids");
        let message = parsed
            .prefetch_stops()
            .expect("the stopping subject's lookup fails with NCBI's exception");
        (subject.ordinal, parsed, message)
    };
    let may_not_look_up = |index: usize, why: &str| {
        let (ordinal, parsed, message) = stopping(index);
        parsed_refusal(
            db.oid(ordinal),
            &format!(
                "has the Seq-align id {}, whose lookup in the database fails ({}; NCBI BLAST+ 2.17.0's prefetch stops with a BLAST Database error), but which NCBI's prefetch may not look up: {why}",
                String::from_utf8_lossy(&parsed.alignment_fasta()),
                // LOSAT's own refusal text.
                String::from_utf8_lossy(message),
            ),
            program,
        )
    };
    let order = load_before_stop(db, subjects);
    let first = order.first_stop_key();
    // The stops that NCBI's loads may end at, in the order of the lookups, up to the key of the
    // first that NCBI surely makes; the sequences given instead for those it may not make.
    let mut ends = Vec::new();
    let mut sure: Option<(u64, usize)> = None;
    let mut unsure: Option<(usize, String)> = None;
    let mut unnamed: Option<(usize, String)> = None;
    let mut instead: Vec<(usize, Vec<u32>)> = Vec::new();
    for &index in &order.stops {
        let key = order.stop_key(index);
        if sure.is_some_and(|(last, _)| key > last) {
            break;
        }
        let (ordinal, parsed, _) = stopping(index);
        match stop_lookup_unsure(
            db,
            record,
            subjects,
            &order.before(first, key),
            queries,
            ordinal,
            parsed,
        ) {
            Ok(None) => {
                sure.get_or_insert((key, index));
            }
            Ok(Some(StopUnsure {
                why,
                instead: given,
            })) => {
                match given {
                    Some(given) => instead.push((index, given)),
                    None => {
                        unnamed.get_or_insert((index, why.clone()));
                    }
                }
                unsure.get_or_insert((index, why));
            }
            Err(error) => return error,
        }
        ends.push(index);
    }
    let Some((last, deciding)) = sure else {
        let (index, why) = unsure.expect("a stop that NCBI may not look up");
        return may_not_look_up(index, &why);
    };
    if let Some((index, why)) = unnamed {
        return may_not_look_up(
            index,
            &format!("{why}; where it does not, NCBI's loads go on to a later stop after a load that LOSAT cannot name"),
        );
    }
    let (_, parsed, message) = stopping(deciding);
    let refuse = |what: &str| {
        parsed_refusal(
            db.oid(subjects[deciding].ordinal),
            &format!(
                "has the Seq-align id {}, whose lookup in the database fails ({}; NCBI BLAST+ 2.17.0's prefetch stops with a BLAST Database error), in a batch of which LOSAT does not know what that prefetch writes before it: {what}",
                String::from_utf8_lossy(&parsed.alignment_fasta()),
                // LOSAT's own refusal text.
                String::from_utf8_lossy(message)
            ),
            program,
        )
    };
    if let Some(&other) = ends.iter().find(|&&index| stopping(index).2 != message) {
        let (_, other_parsed, other_message) = stopping(other);
        return refuse(&format!(
            "the prefetch may also end at the lookup of {}, which fails with another message ({}), as NCBI makes it or not",
            String::from_utf8_lossy(&other_parsed.alignment_fasta()),
            // LOSAT's own refusal text.
            String::from_utf8_lossy(other_message)
        ));
    }
    let before = order.before(first, last);
    if subjects.iter().zip(&before).any(|(subject, before)| {
        *before != LoadBeforeStop::No
            && subject.parsed.is_some_and(|parsed| parsed.prefetch_fails())
    }) {
        return refuse("another of its subjects has an id whose lookup throws, which the prefetch catches and which ends its loads, and which comes before it or may");
    }
    // The loads before the stop, then those whose order LOSAT does not know, then the sequences
    // that a skipped stopping lookup is given instead (loaded by that stopping id). The bytes of a
    // set of loads do not depend on their order ([`prefetch_messages`]: a text comes with the
    // run's first open of the taxonomy database or where a lookup asks it, a volume's first
    // lookup of a taxid, and the lookup of the offset error's taxid that misses its slot needs
    // another taxid of the slot, which `TaxdbRun::cache_collision` refuses; a lookup that reads a
    // record of an empty `taxdb.btd`, where NCBI faults, is refused whichever load makes it,
    // `TaxdbRun::volume_load`), so where those of unknown order write nothing after the others, no
    // part of them writes anything before the stop.
    const INSTEAD: &str = "the sequence that NCBI's BLAST database data loader gives for a stopping id before it whose lookup NCBI may not make, which NCBI loads where it does not";
    let mut unknown = None;
    // Where a load's lookup reads a header with a UTF8String where a VisibleString is expected
    // (`ParsedSubject::scope_utf8`), NCBI's warning comes at that load (`load_read_texts`), after
    // the texts of the loads before it: the loads before the stop come in the order of their
    // Seq-id handles (by key; of one key, in an order that LOSAT does not know, below).
    let scope_utf8 = db.scope_utf8_loads();
    let mut yes: Vec<usize> = (0..subjects.len())
        .filter(|&index| before[index] == LoadBeforeStop::Yes)
        .collect();
    if scope_utf8 {
        yes.sort_by_key(|&index| order.places[index].order);
    }
    let reads = yes
        .iter()
        .map(|&index| (None, index))
        .chain((0..subjects.len()).filter_map(|index| match before[index] {
            LoadBeforeStop::Unknown(why) => Some((Some(why), index)),
            _ => None,
        }))
        .map(|(why, index)| {
            let subject = &subjects[index];
            (why, subject.ordinal, SubjectLoad::Read(subject.parsed))
        });
    // The subject of each load (its index in `subjects`), by its place in `whys`.
    let indexes: Vec<usize> = yes
        .iter()
        .copied()
        .chain(
            (0..subjects.len())
                .filter(|&index| matches!(before[index], LoadBeforeStop::Unknown(_))),
        )
        .collect();
    let given = instead.iter().flat_map(|(index, holders)| {
        let (ordinal, parsed, _) = stopping(*index);
        holders.iter().map(move |&holder| {
            (
                Some(INSTEAD),
                holder,
                SubjectLoad::Instead {
                    stop: ordinal,
                    by: parsed,
                },
            )
        })
    });
    // Nothing loads before the stop: no read, so no refusal of the run's taxonomy cache either
    // (review 1 of line DBID loadord, L5: NCBI writes its error alone). Each load is named by its
    // place in `whys` (`load_read_texts` reads each sequence once).
    let (whys, loads): (Vec<Option<&'static str>>, Vec<_>) = reads
        .chain(given)
        .enumerate()
        .map(|(index, (why, ordinal, load))| (why, (index as u32, ordinal, load)))
        .unzip();
    let mut reads = if loads.is_empty() {
        Vec::new()
    } else {
        match load_read_texts(db, run, loads, false, LoadsOrder::Reads, program) {
            Ok(reads) => reads,
            Err(error) => return error,
        }
    };
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1440-1466
    // (cited at `TaxdbRun::load`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375 (cited at
    // `report::query_warnings::ReadCursor::cut`)
    // A BLAST XML report before the batch that NCBI's stream cuts before its hits makes none of
    // its loads, whose texts (in the model, those of its reads: the texts of opening the taxonomy
    // database at the run's first ask, a volume's offset error at its first lookup of that taxid)
    // go to NCBI's next loads that make the same asks (`report::query_warnings::ReadLoad::left_out`;
    // line RG K4h, review 3 H1 = E2): the loads of the reports and prefetches after it, then those
    // of this prefetch before the stop. BLASTX, which writes each batch's reports before it models
    // the next, has placed them in the run's ledger already (`ReadLoad::place` in
    // `load_read_texts`). The four programs that model every batch before they write a report
    // (`XmlReportsBefore`, [`PrefetchQueries::after_xml_reports`]) place them here by the same
    // model, run over the reports and prefetches of the queries before the batch in the order the
    // writers write them (`LoadTexts::left_before`; review 5 of line DBID loadord, M1). The texts
    // go to the first load that makes the same ask in the order of `reads`: the loads before the
    // stop (whose order among themselves does not change the bytes, as the texts left are the
    // texts of opening, which no load's own text holds after the run's first ask, and offset
    // errors of one text), then those whose place LOSAT does not know, which are refused below
    // where they take one. The loads after the stop ([`load_before_stop`]'s `No`) are not made: a
    // text that no load before the stop asks for is never written (NCBI writes its error alone).
    let mut left_to_unknown = None;
    if let Some(left) = queries
        .xml_before
        .and_then(|before| before.left_texts(queries.error_query))
    {
        for read in &mut reads {
            let own = read.text.len();
            read.place_left(&left);
            if read.text.len() > own {
                if let Some(why) = whys[read.subject as usize] {
                    left_to_unknown.get_or_insert(why);
                }
            }
        }
    }
    if let Some(why) = left_to_unknown {
        return refuse(&format!(
            "the load of another of its subjects, which LOSAT does not know NCBI's prefetch to make before the stop ({why}), asks the taxonomy database for what a BLAST XML report before it did not (NCBI's XML stream cut its Iteration at a query title or id outside US-ASCII, so it made none of its loads), and NCBI writes there the texts of those asks where it makes it"
        ));
    }
    let lead = |read: &crate::report::query_warnings::SubjectRead| {
        read.load.as_ref().map_or(0, |load| load.lead)
    };
    // NCBI's warning of a UTF8String header at a load before the stop: another load of the same
    // key (whose order to it comes from the addresses of NCBI's Seq-id infos) may come before it,
    // which changes the bytes where that load's lookup reads such a header too, or where it asks
    // the taxonomy database and the loads write a text of it.
    if let Some(posting) = reads
        .iter()
        .find(|read| lead(read) > 0 && whys[read.subject as usize].is_none())
    {
        let at = order.places[indexes[posting.subject as usize]];
        let taxonomy_texts = reads.iter().any(|read| read.text.len() > lead(read));
        for read in &reads {
            if std::ptr::eq(read, posting) || whys[read.subject as usize].is_some() {
                continue;
            }
            let index = indexes[read.subject as usize];
            if order.places[index].order_unknown(at).is_none() {
                continue;
            }
            let scoped = subjects[index]
                .parsed
                .is_some_and(|parsed| parsed.scope_utf8_read().is_some());
            let asks = read
                .load
                .as_ref()
                .is_some_and(|load| load.looks_up || load.offset);
            if scoped || (asks && taxonomy_texts) {
                return refuse(&format!(
                    "NCBI's prefetch writes the warning (803.10) of a UTF8String header at the load of database sequence {} before the stop, and the load of database sequence {} may come before it or after it ({ADDRESS_ORDER})",
                    db.oid(subjects[indexes[posting.subject as usize]].ordinal),
                    db.oid(subjects[index].ordinal),
                ));
            }
        }
    }
    let mut text = Vec::new();
    let mut unknown_utf8 = None;
    for read in &reads {
        match whys[read.subject as usize] {
            Some(why) if lead(read) > 0 => {
                unknown_utf8.get_or_insert((indexes[read.subject as usize], why));
            }
            Some(why) if !read.text.is_empty() => {
                unknown.get_or_insert(why);
            }
            Some(_) => {}
            None => text.extend_from_slice(&read.text),
        }
    }
    if let Some((index, why)) = unknown_utf8 {
        let other = subjects[index]
            .parsed
            .and_then(|parsed| parsed.scope_utf8_read())
            .expect("a load whose lookup reads such a header");
        return refuse(&format!(
            "the lookup of another of its subjects' Seq-align id gives database sequence {other} first, whose header holds a UTF8String where a VisibleString is expected, and LOSAT does not know NCBI's prefetch to make that load before the stop ({why}), where it writes its warning (803.10)"
        ));
    }
    if let Some(why) = unknown {
        return refuse(&format!(
            "the load of another of its subjects, which LOSAT does not know NCBI's prefetch to make before the stop ({why}), makes the taxonomy database write a text"
        ));
    }
    anyhow::Error::new(HeaderReadError {
        query: queries.error_query,
        error: NativeError {
            exit: 2,
            message: [
                text.as_slice(),
                DATABASE_ERROR.as_bytes(),
                b" ",
                message,
                b"\n",
            ]
            .concat(),
        }
        .into(),
        utf8: None,
    })
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
// The prefetch rethrows a `CSeqDBException` (`BlastFormatter_PreFetchSequenceData`, cited at
// `ParsedSubject::prefetched`): NCBI stops with `BLAST Database error` (exit 2; probe of line DBID
// P3, `pdb|1ABC|A` and `pdb|1abc|A` in the two databases of `ids/pshare_cd`).
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1465-1466
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
    let (other, or_itself, or_none, reason) = match subject.prefetched() {
        Prefetched::Itself => return Ok(Ok(Some(ordinal))),
        Prefetched::Nothing => return Ok(Ok(None)),
        Prefetched::Failed(reason) => return Ok(Err(reason)),
        // The subject's own blob, or none where another holder's Bioseq read before stands for it
        // (review 1 of line DBID v4late, H2).
        Prefetched::ItselfOrNothing { reason } => {
            if texts && lookup_texts_differ(db, &[&subject.volume_taxids(), &[]]) {
                return Ok(Err(reason));
            }
            return Ok(Ok(Some(ordinal)));
        }
        Prefetched::Another {
            other: None,
            reason,
            ..
        } => return Ok(if texts { Err(reason) } else { Ok(None) }),
        Prefetched::Another {
            other: Some(other),
            or_itself,
            or_none,
            reason,
        } => (other, or_itself, or_none, reason),
    };
    let other = other.ordinal();
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
        Some(_) if or_none => {
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
    // A lookup that throws on the id is the prefetch's warning, one that fails on a file with
    // NCBI's exception its stop (`batch_prefetch`).
    if subject.prefetch_fails() || subject.prefetch_stops().is_some() {
        return Ok(());
    }
    match prefetch_reads(db, ordinal, subject, false, program)? {
        Ok(_) => Ok(()),
        Err(reason) => Err(parsed_refusal(
            db.oid(ordinal),
            &format!("{reason}, {PREFETCH_READS_IT}"),
            program,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1113-1134
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:678-687
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
            "has the Seq-align id {}, by which NCBI BLAST+'s lookup gives database sequence {}, whose definition lines do not hold that id: NCBI BLAST+ 2.17.0's prefetch stops there (oid headers do not contain target gi/seq_id)",
            id(),
            db.oid(other)
        )));
    }
    match read.instead {
        Some(holder) if read.blob.is_none() => Ok(Err(format!(
            "has the Seq-align id {}, by which NCBI BLAST+'s scope gets no sequence unless it read database sequence {holder}, which also holds that id, before: NCBI BLAST+ 2.17.0's prefetch of a minus-strand alignment then writes Error pre-fetching sequence data or not, as its reads decide",
            id()
        ))),
        // Another holder that LOSAT cannot list (review 1 of line DBID v4late, L1).
        None if read.blob.is_none()
            && matches!(
                subject.found(),
                Some(crate::report::db_subject::FoundInstead::HoldersUnknown { .. })
            ) =>
        {
            Ok(Err(format!(
                "has the Seq-align id {}, by which NCBI BLAST+'s scope gets no sequence unless it read before another sequence that holds that id, which LOSAT cannot list: NCBI BLAST+ 2.17.0's prefetch of a minus-strand alignment then writes Error pre-fetching sequence data or not, as its reads decide",
                id()
            )))
        }
        _ => Ok(Ok(!read.bioseq)),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1129-1140
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
            .taxonomy_lookups(Some(subject.alignment_id())),
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
    /// Another sequence that holds an id the scope takes for the Seq-align id besides `blob`
    /// (`ParsedSubject::other_holder`): where the scope holds its Bioseq from an earlier read, it
    /// can get that one instead (psiblast's PSSM refuses every such read,
    /// `ParsedSubject::pssm_read`; the prefetch's handle where `blob` is none,
    /// `prefetch_gets_no_bioseq`).
    pub(crate) instead: Option<crate::seqdb::DbSequence>,
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1340-1372 and
// 1437-1439
// NCBI's scope gets a Seq-align id's Bioseq from the data loader (`x_GetOid`, cited at
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
        Some(
            FoundInstead::Sequence(other)
            | FoundInstead::HoldersUnknown {
                instead: Some(other),
                ..
            },
        ) => {
            let set = db
                .defline_set(other.ordinal())
                .map_err(|error| database_error(error, program))?;
            Ok(other_bioseq(&set, other.ordinal(), subject).unwrap_or(none))
        }
        Some(FoundInstead::Nothing | FoundInstead::ThrowsOnId(_) | FoundInstead::Failed(_)) => {
            Ok(none)
        }
        None
        | Some(
            FoundInstead::AlsoHeld(_)
            | FoundInstead::HoldersUnknown { instead: None, .. }
            | FoundInstead::OtherTaxids(_)
            | FoundInstead::CopyApart(..),
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
/// and that shows no taxonomy of them (`prefetch_messages`): outfmt 0 to 7 and 17, and with
/// `PRE_FETCH_SEQS_LIMIT=2147483647` the Seq-annots 8, 9 and 12 too (`s_PreFetchSeqs`, cited at
/// [`prefetch_may_be_off`], returns before its test of the format). The query warnings, which
/// the batch's messages go with, are written with the first format (the Seq-annots write only
/// the reader's); the XML2 and JSON formats write theirs at each `Hit` (`TaxdbPrinter`), and the
/// archive, which NCBI's apps but psiblast write without calling the prefetch, prefetches nothing
/// (psiblast_app.cpp:262-264 prefetches before its test of the archive; LOSAT's PSIBLAST refuses
/// outfmt 8, 9, 11 and 12).
pub(crate) fn prefetches_subjects(outfmt: Option<i32>) -> bool {
    match outfmt {
        // `s_IsPrefetchFormat` takes the query-anchored formats 1 to 4 as 0 (DECISIONS T3-4); the
        // limit's rule of the Seq-annots is that of the prefetch's reads.
        Some(format @ (0..=9 | 12 | 17)) => prefetch_reads_seq_aligns(format),
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

/// Whether [`read_reported_taxonomy`] or [`prefetch_messages`] could refuse for some sequence of
/// the database, which a run that writes its batches one by one cannot take back (OF5b BLASTX,
/// DECISIONS 32): a name read past NCBI's mapping of `taxdb.btd` or a record of an empty one
/// (the prefetch's loads refuse only the second, `TaxdbRun::volume_load`), or a header or the
/// taxonomy database that cannot be read. A superset of what a search
/// reports: limits (`-gilist`, `-taxidlist`) and the hits themselves are ignored. `false`
/// without `taxdb.bti` on the search path and for a database without taxids. The refusal of
/// NCBI's cache (`TaxdbRun::cache_collision`) is not asked here: it is a property of the run
/// (`TaxdbRun::cache_may_refuse`), which BLASTX and BLASTP's block search ask apart.
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
pub(crate) fn taxdb_may_refuse(db: &SeqDb) -> bool {
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
        if taxdb_side(&taxdb) == Some(false) {
            return false;
        }
        if sequences > PREPASS_DB_SEQUENCES_MAX {
            return true;
        }
    }
    database_side(db, &taxdb)
}

/// The taxonomy side of [`taxdb_may_refuse`]: `Some(false)` when no lookup of any taxid can be
/// refused, which is only where a lookup reads a record past NCBI's mapping of `taxdb.btd`
/// ([`crate::seqdb::taxdb::TaxDb::any_lookup_may_fail`], whose pass holds for any order of the
/// index; no scientific name is refused for its bytes: NCBI's streams write them); `Some(true)`
/// when one can; `None` when the taxonomy database cannot be opened or read, which only the
/// lookups of the database's taxids can answer.
fn taxdb_side(taxdb: &crate::seqdb::taxdb::LazyTaxDb) -> Option<bool> {
    taxdb.any_lookup_may_fail().ok()
}

/// The database side of [`taxdb_may_refuse`]: whether the lookup of the first definition line's
/// taxid of some sequence would be refused, looking each distinct non-zero taxid up with the
/// `names` code `read_reported_taxonomy` uses (refused: a record of an empty `taxdb.btd`, or one
/// past NCBI's mapping without its third tab inside it, `TaxLookup::into_names`; the prefetch's
/// loads, which write no name, refuse only the first). The headers are read in
/// ranges of [`crate::seqdb::HEADER_RANGE_SEQS`] sequences, one read per MiB
/// ([`SeqDb::visit_defline_taxids`]). A header that cannot be read answers `true`. Memory: the
/// distinct taxids, nothing per sequence.
fn database_side(db: &SeqDb, taxdb: &crate::seqdb::taxdb::LazyTaxDb) -> bool {
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
            refused = taxdb.names(taxid, &mut messages).is_err();
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
#[path = "database_taxdb_notes_tests.rs"]
mod taxdb_notes_tests;

#[cfg(test)]
#[path = "database_taxdb_prepass_tests.rs"]
mod taxdb_prepass_tests;
