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

use crate::blastinput::app::ErrorAfterBatchReports;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::defline::SubjectWord;
use crate::report::pairwise::{DatabaseInfo, PairwiseDatabase};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

// NCBI reference (598d8ae6): c++/src/corelib/ncbi_system.cpp:620-621
/// Parses Linux `/sys/devices/system/cpu/online` CPU lists (for example `0`, `0-31`, `0,2-5`).
#[cfg(any(target_os = "linux", test))]
fn parse_linux_cpu_list(text: &str) -> Option<usize> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return None;
    }
    let mut total = 0usize;
    for part in trimmed.split(',') {
        let part = part.trim();
        let count = if let Some((lo, hi)) = part.split_once('-') {
            let lo = lo.trim().parse::<usize>().ok()?;
            let hi = hi.trim().parse::<usize>().ok()?;
            hi.checked_sub(lo)?.checked_add(1)?
        } else {
            let _ = part.parse::<usize>().ok()?;
            1usize
        };
        total = total.checked_add(count)?;
    }
    (total > 0).then_some(total)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbi_system.cpp:616-630
/// Counts `cpu<N>` entries in Linux `/proc/stat`, matching glibc's `_SC_NPROCESSORS_ONLN`
/// fallback when `/sys/devices/system/cpu/online` is unavailable.
#[cfg(any(target_os = "linux", test))]
fn parse_linux_proc_stat_cpus(text: &str) -> Option<usize> {
    let count = text
        .lines()
        .filter(|line| {
            line.strip_prefix("cpu")
                .is_some_and(|rest| rest.as_bytes().first().is_some_and(u8::is_ascii_digit))
        })
        .count();
    (count > 0).then_some(count)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbi_system.cpp:598-639
/// Online CPU count matching NCBI's `CSystemInfo::GetCpuCount` (process-cached as in NCBI).
///
/// On Linux, `std::thread::available_parallelism` caps by `sched_getaffinity` (`taskset`) and
/// cgroup quotas, whereas NCBI calls `sysconf` with `_SC_NPROCESSORS_ONLN`. Match glibc by reading
/// `/sys/devices/system/cpu/online`, falling back to `/proc/stat` and then `available_parallelism`.
#[cfg(not(target_arch = "wasm32"))]
fn ncbi_online_cpu_count() -> usize {
    static ONLINE_CPUS: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *ONLINE_CPUS.get_or_init(|| {
        #[cfg(target_os = "linux")]
        {
            if let Some(n) = std::fs::read_to_string("/sys/devices/system/cpu/online")
                .ok()
                .and_then(|s| parse_linux_cpu_list(&s))
            {
                return n;
            }
            if let Some(n) = std::fs::read_to_string("/proc/stat")
                .ok()
                .and_then(|s| parse_linux_proc_stat_cpus(&s))
            {
                return n;
            }
        }
        std::thread::available_parallelism().map_or(1, usize::from)
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// NCBI reference (598d8ae6): c++/src/corelib/ncbi_system.cpp:598-638
/// Requested `-num_threads` capped at the host's online CPU count, matching NCBI's
/// `CMTArgs::x_ExtractAlgorithmOptions` (`CSystemInfo::GetCpuCount`).
///
/// On `wasm32` (WASI Preview 1), `wasi-libc` hardcodes `sysconf` for `_SC_NPROCESSORS_ONLN` to 1
/// while `validate_threads` and the host cap already bound `num_threads`, so keep `num_threads`
/// uncapped there to match native CLI DB-split behavior.
pub(crate) fn ncbi_effective_threads(num_threads: usize) -> usize {
    #[cfg(target_arch = "wasm32")]
    {
        num_threads
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        num_threads.min(ncbi_online_cpu_count())
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:96-171
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:99-174
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:96-173
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:153-154
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:339-344,699-707
/// Whether NCBI's application runs this search in multithreaded DB-split mode (positive
/// `m_NumThreads` on `CSeqDB`, so `CSeqDBImpl::x_FillSeqBuffer` stops scanning at the first OID
/// whose `GetSequence` length is negative).
///
/// With effective threads `<= 1` (after capping `-num_threads` at the host's CPU count) or
/// `-mt_mode 1` (split by queries), each search node runs with 1 thread on `CSeqDB`
/// (`m_NumThreads` stays zero). With `-mt_mode 0` (auto), NCBI splits by queries only when
/// `total_length < max_db_size` and the query file size exceeds `min_q_size`.
pub(crate) fn ncbi_mt_splits_db(
    num_threads: usize,
    mt_mode: Option<u8>,
    total_length: u64,
    query_file: Option<&std::path::Path>,
    max_db_size: u64,
    min_q_size: u64,
) -> bool {
    if ncbi_effective_threads(num_threads) <= 1 {
        return false;
    }
    match mt_mode {
        Some(1) => false,
        Some(2) => true,
        _ => {
            if total_length >= max_db_size {
                return true;
            }
            let query_size = query_file
                .filter(|path| path.as_os_str() != "-")
                .and_then(|path| std::fs::metadata(path).ok())
                .map(|meta| meta.len())
                .unwrap_or(0);
            query_size <= min_q_size
        }
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
/// Opens the database `name` (the `-db` value as NCBI keeps it, `seqdb::db_argument_name`)
/// of the sequence type `seqtype` for a search: the one function through which a search
/// program opens a database. The name is looked for as NCBI looks for it: in the current
/// directory, then in the directories of the `BLASTDB` variable, then in those of the
/// registry entry `[BLAST] BLASTDB` (`settings.blastdb`, read from
/// `NCBI_CONFIG__BLAST__BLASTDB`, the program's `.ini` file and `.ncbirc`).
pub(crate) fn open_search_database(
    name: &[u8],
    seqtype: SeqType,
    settings: &ApplicationSettings,
) -> Result<SeqDb, SeqDbError> {
    let path = SearchPath::from_environment(settings.blastdb.as_deref());
    // The errors of the alias files' counts depend on the list, and those of the headers come
    // after them (`restrict_database`).
    let mut db = SeqDb::open_search_in(name, seqtype, &path)?;
    db.fail_search_reads_with_oid_list();
    // The checks of NCBI's reads of the definition lines (`SERIAL_VERIFY_DATA_READ`).
    db.refuse_unchecked_reads(&settings.read_verification);
    Ok(db)
}

/// Opens the database of `-db name` (`open_search_database`) with the command line's
/// errors (`database_error`; LOSAT's refusals name `program`).
pub(crate) fn open_database(
    name: &[u8],
    seqtype: SeqType,
    settings: &ApplicationSettings,
    program: &str,
) -> Result<SeqDb> {
    open_search_database(name, seqtype, settings).map_err(|error| database_error(error, program))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:109-139
// The list is resolved in the totals scan (`x_ScanTotals`, `CSeqDBGiListSet`: its errors and
// warnings), before the alias files' statistics keys are read.
/// Restricts the search's database `db` (opened as `name` by `open_database`) to the taxid or
/// seqid list of the command line (line DI-8, `DbFilterArgs::resolve`; none: the whole
/// database), as NCBI's `CSeqDB` opened with the list (`SeqDb::restrict`), then reads what its
/// alias files' counts and the reports need: NCBI's errors are `BLAST Database error: <text>`
/// with exit status 2 (`database_error`), and the warnings NCBI writes go to `diagnostics`
/// (`LOG_POST` at the warning level: the text alone; `ERR_POST`: `Warning: [<app>] <text>`,
/// the running application `app` in lower case as `SetDiagPostPrefix` names it (the program's
/// own, or that of a command line that runs another program's strategy); an `ERR_POST` of an
/// error before a throw: `Error: [<app>] <text>`), before an error too.
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
    name: &[u8],
    list: Option<&crate::seqdb::DbList>,
    mask: Option<&crate::blastinput::db_mask::DbMaskOption>,
    diagnostics: &mut dyn std::io::Write,
    settings: &ApplicationSettings,
    program: &str,
    app: &str,
) -> Result<bool> {
    if let Some(list) = list {
        apply_list(db, list, diagnostics, program, app)?;
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
    app: &str,
) -> Result<()> {
    let mut warnings = Vec::new();
    let restricted = db.restrict(list, &mut warnings);
    let prefix = app.to_ascii_lowercase();
    for warning in warnings {
        diagnostics.write_all(&list_warning_line(warning, Some(&prefix)))?;
    }
    restricted.map_err(|error| database_error(error, program))
}

/// A warning of `SeqDb::restrict` as NCBI posts it, with its line feed: after the diagnostic
/// prefix `[prefix] ` of the application, or none where it sets none (blast_formatter's `Run`).
/// The lines keep the bytes of a seqid list's words (probe R3 of line LIB-B).
pub(crate) fn list_warning_line(
    warning: crate::seqdb::ListWarning,
    prefix: Option<&str>,
) -> Vec<u8> {
    let tag = prefix.map(|p| format!("[{p}] ")).unwrap_or_default();
    let mut line = match warning {
        crate::seqdb::ListWarning::Log(text) => text,
        crate::seqdb::ListWarning::Posted(text) => format!("Warning: {tag}{text}").into_bytes(),
        crate::seqdb::ListWarning::Error(text) => {
            [format!("Error: {tag}").as_bytes(), &text].concat()
        }
        crate::seqdb::ListWarning::SeqId(warning) => {
            crate::writedb::fasta_source::defline_ids::seq_id_warning_line(prefix, warning)
        }
    };
    line.push(b'\n');
    line
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
    list_gis: &[i64],
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
    options.set_database_filter(
        filter.archive_parameters(list_gis, list_taxids),
        request_length,
    )?;
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
/// translated programs that take it) run by the application `app` (NCBI's lower-case name, the
/// prefix of the warning): the resolution (`DbMaskOption::resolve`), whose warning goes
/// to `warning` (one diagnostic, which an archive saves as one message), or whose database error
/// ends the run before the formatter's other work; then the program options of the archive
/// `archive` and of the search strategy `strategy` (`note_archive_database_masking`). The result
/// is the search's masking (`None`: not masked).
pub(crate) fn resolve_database_masking(
    db: &SeqDb,
    db_name: &[u8],
    option: &crate::blastinput::db_mask::DbMaskOption,
    program: &str,
    app: &str,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    warning: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<Option<crate::blastinput::db_mask::ResolvedMask>> {
    let resolved = option.resolve(db, db_name, program, app, warning)?;
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
// place of each throw: one thrown by the code NCBI's `datatool` generates (`ThrowSite::generated`)
// names a place LOSAT does not have, and LOSAT refuses it with its text, as segmasker does
// (`segmasker::work::database_failure`); the `CObject_id::GetId` of a `BL_ORD_ID` or `PIG` id with
// a string tag (under an id list, line RG D1c-1; a Seq-align id, line DBID pfoff-7) has the place
// that NCBI's report names (`seqdb::isam::OBJECT_ID_GET_ID`).
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
/// WQA D18, 2026-11-21), and so is a LOSAT internal error ([`SeqDbError::Internal`]); a setting
/// of the run that LOSAT refuses ([`SeqDbError::Setting`]) is refused in `program`'s name; any
/// other is LOSAT's refusal of the database in `program`'s name ([`refusal`]). The one place of that choice
/// for the search programs, blastdbcmd and the maskers' database reader; each ends with LOSAT's
/// exit status of its own errors.
pub(crate) fn losat_reader_error(error: &SeqDbError, program: &str) -> String {
    if error.is_unreadable() || matches!(error, SeqDbError::Internal(_)) {
        error.to_string()
    } else if let SeqDbError::Setting(refused) = error {
        // A setting of the run, not the database, is refused (`SeqDb::refuse_unchecked_reads`).
        format!("{refused} is not supported by LOSAT's {program}")
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-256
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
/// A reader error while the preliminary search reads the sequences: NCBI's search stage
/// replaces a `CSeqDBException` (`SeqDbError::ncbi_search_message`) and wraps any other known
/// `std::exception` (`SeqDbError::ncbi_preliminary_message`), and the application
/// reports `BLAST engine error: <text>` with exit status 3. The error is kept apart from
/// the others (`PreliminaryReadError`), because NCBI has written the outfmt 0 prolog and
/// read the query batch before it.
pub(crate) fn search_read_error(error: SeqDbError, program: &str) -> anyhow::Error {
    search_read_error_of(&error, program)
}

/// [`search_read_error`] of an error that its owner keeps (`seqdb::MaskBlock`: the masks that the
/// preliminary search reads with each sequence, seqsrc_seqdb.cpp:318-321).
pub(crate) fn search_read_error_of(error: &SeqDbError, program: &str) -> anyhow::Error {
    match error.ncbi_preliminary_message() {
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
/// (`search_read_error`), with the bytes after `BLAST engine error: `.
#[derive(Debug)]
pub(crate) struct PreliminaryReadError(pub(crate) std::borrow::Cow<'static, [u8]>);

impl PreliminaryReadError {
    /// A static ASCII text of `CLocalBlast`'s preliminary search error.
    pub(crate) const fn from_static(text: &'static str) -> Self {
        Self(std::borrow::Cow::Borrowed(text.as_bytes()))
    }

    /// The bytes after `BLAST engine error: `.
    pub(crate) fn as_bytes(&self) -> &[u8] {
        self.0.as_ref()
    }
}

impl std::fmt::Display for PreliminaryReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{ENGINE_ERROR} {}", String::from_utf8_lossy(&self.0))
    }
}

impl std::error::Error for PreliminaryReadError {}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:304-309
/// NCBI's first read, in a run of reported sequences' headers (`HeaderReader`), of a record that
/// holds a UTF8String where a VisibleString is expected: the query whose hit list reached it
/// (numbered as `ErrorAfterBatchReports::query`) and the member NCBI's warning names
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
    /// read, or the one that a failed read made before it failed (`ErrorAfterBatchReports::utf8`). NCBI's
    /// warning of it comes after the batch's reader messages and request posts, before its
    /// prefetch and before its error.
    pub(crate) fn of_read<'r, T>(
        read: &'r Result<T>,
        utf8: impl FnOnce(&'r T) -> Option<&'r Utf8Read>,
    ) -> Option<&'r Utf8Read> {
        match read {
            Ok(read) => utf8(read),
            Err(error) => error
                .downcast_ref::<ErrorAfterBatchReports>()
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

    /// Whether the read is of the batch whose first query is `batch`, which posts its warning.
    pub(crate) fn of_batch(&self, batch: usize) -> bool {
        self.0.get().is_some_and(|(first, _)| *first == batch)
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1568-1582
/// NCBI's posts of the titles that a failed batch's Seq-aligns read before the error
/// (`ErrorAfterBatchReports::title_reads`), which posted NCBI's warning of a UTF8String header or
/// not (`utf8_posted`), with the post prefix of `program`, as [`batch_title_posts`] makes a
/// batch's: the failed read posts what it read of its titles before it failed
/// (`seqid::ber::title_replaced_chars`). No prefetch or formatter of the batch reads anything
/// after them, so they are NCBI's in every format. LOSAT's refusal where a posting read's volume
/// filters the header records it reads ([`check_reported_title`]: NCBI's fetches read them too),
/// and `batch_title_posts`' with that warning. The caller writes them after the warning, before
/// the error (and in the run's BLAST archive, `Posts::write`).
#[cold]
#[inline(never)]
pub(crate) fn failed_batch_title_posts(
    title_reads: &[u32],
    db: &SeqDb,
    utf8_posted: bool,
    program: &str,
) -> Result<crate::report::diag_log::Posts> {
    if title_reads.is_empty() {
        return Ok(crate::report::diag_log::Posts::default());
    }
    let prefix = program.to_ascii_lowercase();
    let posts_at = |at: crate::seqdb::SeqLocation| !db.title_posts_at(at, &prefix).is_empty();
    // A run whose reader met no such title before (none of its records in the cache posts)
    // models the cache from these reads on, from an empty one.
    if !db.replaced_titles()
        && title_reads
            .iter()
            .any(|&ordinal| db.location_of(ordinal).is_some_and(posts_at))
    {
        db.note_replaced_title();
    }
    if db.replaced_titles() {
        for &ordinal in title_reads {
            if db
                .location_of(ordinal)
                .is_some_and(|at| db.filters_header_sets(at) && posts_at(at))
            {
                check_reported_title(db, ordinal, true, TitleReads::SeqAligns, program)?;
            }
        }
    }
    batch_title_posts(db, title_reads.iter().copied(), utf8_posted, None, program)
}

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
                        record.name.clone()
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
pub(crate) fn check_formats_report_sums<'r, O: AsRef<[u8]>>(
    shown: impl FnOnce() -> PairwiseDatabase,
    formats: impl IntoIterator<Item = (O, &'r mut crate::api::local_blast::FormatRefusal)>,
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
/// NCBI's own is an `ErrorAfterBatchReports` that names the query (NCBI fails that query's batch);
/// LOSAT's refusals (`program` names the program) are errors of their own.
pub(crate) fn stored_title(
    deflines: std::result::Result<Vec<crate::seqdb::OrdinalDefline>, crate::seqdb::SeqDbError>,
    oid: impl Fn() -> crate::seqdb::Oid,
    query: usize,
    program: &str,
) -> Result<StoredTitle> {
    let mut deflines =
        deflines.map_err(|error| header_read_error(error, query, program, (None, None)))?;
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
/// `ErrorAfterBatchReports` that names the query, so the batch fails after the outfmt 0 prolog as NCBI's
/// does; so is LOSAT's refusal of a setting of the run that NCBI's read of that record would
/// meet ([`SeqDbError::Setting`], `SeqDb::refuse_unchecked_reads`), which stops where NCBI's
/// default error does, after the same output; LOSAT's other refusals are errors of their own.
/// `read`: the sequence whose read failed where NCBI makes the batch's Seq-aligns
/// (`ErrorAfterBatchReports::title_reads`).
fn header_read_error(
    error: SeqDbError,
    query: usize,
    program: &str,
    (utf8, read): (Option<&Utf8Read>, Option<u32>),
) -> anyhow::Error {
    if error.ncbi_message().is_some()
        || reported_exception(&error).is_some()
        || matches!(error, SeqDbError::Setting(_))
    {
        anyhow::Error::new(ErrorAfterBatchReports {
            query: Some(query),
            error: database_error(error, program),
            utf8: utf8.cloned(),
            title_reads: read.into_iter().collect(),
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

/// Whether the sequences of a database search keep a definition line with a Seq-id
/// (`SeqDb::line_read`): the one check of NCBI's two places that drop a sequence without one,
/// the fetches with `check_oid_exclusion` (`FetchExcluded`) and the Seq-aligns (`SeqAligned`),
/// with whether NCBI's read of the record there takes a UTF8String for a VisibleString (its
/// warning, `seqdb::Utf8Reads`). Under a list (`SeqDb::drops_lineless_hits`) each header is
/// decoded once and kept in the cache; without one the check walks the record to its first
/// Seq-id at each ask (`seqdb::header::holds_a_seqid`) and keeps nothing. A header that cannot
/// be read as NCBI's header read adjusts it is not kept in the cache: each place decides what
/// its error does.
struct LineCheck<'a> {
    db: &'a SeqDb,
    /// `SeqDb::drops_lineless_hits`: whether the reads decode the header (and the cache keeps it).
    filters: bool,
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
            filters: db.drops_lineless_hits(),
            kept: std::sync::Mutex::new(std::collections::HashMap::new()),
        }
    }

    /// What the read of `ordinal`'s header finds (under a list read at its first check; the lock
    /// is not held while the header is read and decoded), or the reader's error.
    fn read(&self, ordinal: u32) -> Result<LineRead, SeqDbError> {
        if !self.filters {
            let (kept, utf8) = self.db.line_read(ordinal)?;
            return Ok(LineRead { kept, utf8 });
        }
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
/// under a negative seqid list (`SeqDb::excludes_lineless_at_fetch`), a sequence that
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
/// such a search there (`refusal`). The composition-based redo (and the Smith-Waterman redo)
/// fetches only the matches that its early termination test leaves (blast_kappa.c:3525-3540):
/// a stream's list whose header cannot be read stays in it (`marks_noting`,
/// `FetchedReads::unreadable`), and the redo's program refuses the search only where NCBI's redo
/// fetches it (`refusal_at`, found as the reads of a UTF8String header are: `first_redo_fetch`).
pub(crate) struct FetchExcluded<'a> {
    check: LineCheck<'a>,
    /// The first of the reads of `excludes` that took a UTF8String for a VisibleString, in the
    /// ordinary traceback's order (`FetchedUtf8::traceback_first`).
    utf8: crate::seqdb::Utf8Reads,
}

impl<'a> FetchExcluded<'a> {
    /// The check of `db`'s fetches; `None` where they exclude nothing (no negative seqid
    /// list).
    pub(crate) fn of(db: &'a SeqDb) -> Option<Self> {
        db.excludes_lineless_at_fetch().then(|| FetchExcluded {
            check: LineCheck::new(db),
            utf8: Default::default(),
        })
    }

    /// The database whose fetches this checks.
    pub(crate) fn db(&self) -> &'a SeqDb {
        self.check.db
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
            .map_err(|error| read_error(self.check.db, ordinal, error))?;
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
        // The traceback fetches every subject: the first header that cannot be read ends it.
        let (marks, reads) = self.read_marks(items, ordinal, |_| (), false)?;
        marks.retain(items);
        Ok(reads.utf8.into_iter().map(|(read, ())| read).collect())
    }

    /// `retain` of the lists of a composition-based redo's stream (which its fetches read in
    /// the stream's order), each `FetchedUtf8` with `noted(item)` of its item (what its redo's
    /// order and early termination test read of it, `redo_fetches`): a list whose header cannot
    /// be read stays (`FetchedReads::unreadable`), and its program refuses the search where the
    /// redo fetches it (`refusal_at`).
    pub(crate) fn retain_noting<T, N>(
        &self,
        items: &mut Vec<T>,
        ordinal: impl Fn(&T) -> Option<u32>,
        noted: impl Fn(&T) -> N,
    ) -> FetchedReads<N> {
        let (marks, reads) = self.marks_noting(items, ordinal, noted);
        marks.retain(items);
        reads
    }

    /// Which of `items` (the preliminary hit lists of a redo's stream, in the stream's order)
    /// are on a subject the fetch excludes (`FetchMarks`); `ordinal` names an item's subject
    /// (`None`: fetched). The items whose header holds a UTF8String where a VisibleString is
    /// expected and those whose header cannot be read come with them (`FetchedReads`, each with
    /// `noted(item)`), as `retain_noting` returns them: the one read of the items' headers for
    /// all. A list whose header cannot be read is neither kept nor excluded by NCBI's fetch,
    /// which ends the search where its redo reaches it; it stays in the stream (marked as kept),
    /// so that the redo tests its early termination as NCBI's does.
    pub(crate) fn marks_noting<T, N>(
        &self,
        items: &[T],
        ordinal: impl Fn(&T) -> Option<u32>,
        noted: impl Fn(&T) -> N,
    ) -> (FetchMarks, FetchedReads<N>) {
        match self.read_marks(items, ordinal, noted, true) {
            Ok(marks) => marks,
            Err(_) => unreachable!("a deferred read error is noted, not returned"),
        }
    }

    /// `marks_noting`; with `defer` false, the first header that cannot be read, in the order
    /// of `items`, is the error (with its ordinal) instead (a fetch of every item: `retain`).
    fn read_marks<T, N>(
        &self,
        items: &[T],
        ordinal: impl Fn(&T) -> Option<u32>,
        noted: impl Fn(&T) -> N,
        defer: bool,
    ) -> Result<(FetchMarks, FetchedReads<N>), (u32, SeqDbError)> {
        let mut excluded = Vec::with_capacity(items.len());
        let mut reads = FetchedReads {
            utf8: Vec::new(),
            unreadable: Vec::new(),
        };
        let mut kept = 0;
        for (index, item) in items.iter().enumerate() {
            let keeps = match ordinal(item) {
                Some(ordinal) => match self.check.read(ordinal) {
                    Ok(read) => {
                        if read.utf8 {
                            reads.utf8.push((
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
                    Err(error) if !defer => {
                        return Err((ordinal, read_error(self.check.db, ordinal, error)));
                    }
                    Err(_) => {
                        reads.unreadable.push((
                            FetchedUtf8 {
                                index,
                                kept: true,
                                before: kept,
                                ordinal,
                            },
                            noted(item),
                        ));
                        true
                    }
                },
                None => true,
            };
            kept += usize::from(keeps);
            excluded.push(!keeps);
        }
        Ok((FetchMarks(excluded), reads))
    }

    /// LOSAT's refusal (in `program`'s name) of a search of `db` whose redo fetches the sequence
    /// `ordinal`, whose header cannot be read (`FetchedReads::unreadable`): `refusal` of the
    /// error of a new read of it.
    #[cold]
    #[inline(never)]
    pub(crate) fn refusal_at(db: &SeqDb, ordinal: u32, program: &str) -> anyhow::Error {
        match db.line_read(ordinal) {
            Err(error) => Self::refusal((ordinal, read_error(db, ordinal, error)), program),
            // Not met: the reader keeps no state between two reads of a header.
            Ok(_) => FetchRefusal(
                refusal(
                    &format!(
                        "under the negative seqid list NCBI reads the definition lines of each \
                         sequence that its search aligns again, and those of sequence {ordinal} \
                         cannot be read"
                    ),
                    program,
                )
                .to_string(),
            )
            .into(),
        }
    }

    /// LOSAT's refusal (in `program`'s name) of a search whose fetch of `ordinal` cannot read
    /// its header (`error`), where NCBI's search fails (`FetchRefusal`); a header file whose
    /// bytes the system does not give, and a setting of the run that LOSAT refuses, keep their
    /// own words (`losat_reader_error`).
    pub(crate) fn refusal((ordinal, error): (u32, SeqDbError), program: &str) -> anyhow::Error {
        if error.is_unreadable() || matches!(error, SeqDbError::Setting(_)) {
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
/// list, `FetchMarks::default`). NCBI's redo loop tests such a list's early termination first,
/// like any other's; where the test does not stop the list, the fetch excludes the subject and
/// NCBI sets the status of the thread back to 0 and goes on with the next list, with no redo and
/// no heap insert, so a failed Smith-Waterman redo before it no longer fails the batch; where
/// the test stops the list, nothing changes (core's `SwBatchStatus::reaches_redo`, the order of
/// the two steps for every program's redo loop). One owner of the question "is this list's
/// subject excluded at the fetch", with two readers: TBLASTN's `stage_d_pipeline` keeps the marked
/// lists in its stream and asks `SwBatchStatus::reaches_redo` at each (the X-drop redo has no
/// status, so a marked list is a list skipped there); BLASTP's `exclude_fetched` removes them and
/// keeps their places (`blastp::blast_engine::database::FetchExcludedMatch`), which both of its
/// redo loops read through the same rule with one function
/// (`blastp::blast_engine::redo_run::query_fetch_excluded_reached`): the one-thread loop after
/// each query's loop, and the replay of the parallel redo's windows when it passes the query
/// (`blastp::x_redo_par`, `XSwBatch`). The other programs remove such lists with
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
    ) -> (Self, FetchedReads<N>) {
        match excluded {
            Some(excluded) => excluded.marks_noting(items, ordinal, noted),
            None => (Self::default(), FetchedReads::default()),
        }
    }

    /// `items` (those the marks were made of) without the excluded ones, in their order.
    fn retain<T>(&self, items: &mut Vec<T>) {
        let mut index = 0;
        items.retain(|_| {
            index += 1;
            !self.excluded(index - 1)
        });
    }

    /// Whether the fetch excludes the subject of the stream's item `index`.
    #[inline(always)]
    pub(crate) fn excluded(&self, index: usize) -> bool {
        self.0.get(index).copied().unwrap_or(false)
    }
}

/// What the one read of the headers of a stream's lists (`FetchExcluded::marks_noting`) found
/// beyond the marks, each subject with its item's note (`FetchedUtf8`, and what its caller's
/// redo reads of the item).
#[derive(Debug)]
pub(crate) struct FetchedReads<N> {
    /// The subjects whose header holds a UTF8String where a VisibleString is expected: the first
    /// that NCBI's fetches read posts its warning.
    pub(crate) utf8: Vec<(FetchedUtf8, N)>,
    /// The subjects whose header cannot be read as NCBI's read adjusts it (a string `BL_ORD_ID`
    /// tag in a volume after the first, a record NCBI's read fails on) or that LOSAT refuses
    /// where it reads it (`seqdb::Volume::filtered_mask`, whose warning or error NCBI gives at
    /// that read): their lists stay in the stream, and where NCBI's redo fetches one (found as the first fetched UTF8String read is,
    /// `first_redo_fetch`) its search ends; LOSAT refuses it there
    /// (`FetchExcluded::refusal_at`). Their `kept` is true.
    pub(crate) unreadable: Vec<(FetchedUtf8, N)>,
}

impl<N> Default for FetchedReads<N> {
    fn default() -> Self {
        FetchedReads {
            utf8: Vec::new(),
            unreadable: Vec::new(),
        }
    }
}

/// The error of the read of `ordinal`'s header of `db` for `FetchExcluded::refusal`: a record
/// that NCBI's read fails on (`SeqDbError::Ncbi`, the text of NCBI's exception, which names no
/// file) is named by its header file, as the reader's other errors are.
#[cold]
#[inline(never)]
fn read_error(db: &SeqDb, ordinal: u32, error: SeqDbError) -> SeqDbError {
    if !matches!(error, SeqDbError::Ncbi(_)) {
        return error;
    }
    match db.header_file(ordinal) {
        Some(path) => SeqDbError::Format {
            path,
            message: error.to_string(),
        },
        None => error,
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
/// a hit on a sequence that the database's list keeps without any of its definition
/// lines, or whose lines hold no Seq-id (an edited header, without a list too; line RG
/// D1c-1c-b), gets none (`SeqDb::keeps_line`, blast_seqalign.cpp:1582-1586), so it leaves every
/// report, the prefetch, the taxonomy lookups and psiblast's PSSM, while it was searched, took
/// its place in the hit list before the cut to `-max_target_seqs` and counts in the
/// statistics. Each program drops such hit lists right after that cut, before the reads of
/// the reported sequences; a check reads each sequence's header once (one check per query
/// batch or search call; the reads of the reports decode the kept ones again). Under a
/// negative seqid list the search's fetches have excluded them before
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
    /// The check of `db`'s hits; `None` where no hit is dropped here: the search's fetches have
    /// excluded the sequences without a line before (`fetched`, the caller's: the search fetches
    /// its subjects with `check_oid_exclusion`, a gapped search or an ungapped one of a nucleotide
    /// database; and the database's list is negative, `SeqDb::excludes_lineless_at_fetch`, which
    /// erases every line without a Seq-id).
    pub(crate) fn of(db: &'a SeqDb, fetched: bool) -> Option<Self> {
        (!(fetched && db.excludes_lineless_at_fetch())).then(|| SeqAligned {
            check: LineCheck::new(db),
            utf8: None,
            failed: false,
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
    /// The lookup in the LMDB files that have the accession table
    /// (`SeqidLookup::in_indexed_files`), made at the first lookup that fails on a file
    /// (`failed_holders`), or why it is not made.
    indexed_lookup:
        std::cell::OnceCell<std::result::Result<crate::seqdb::lookup::SeqidLookup<'a>, String>>,
    /// Whether calls of the run before this one loaded sequences that no record holds
    /// (`TaxdbRun::unrecorded_calls`), so that the other holders of each Seq-align id are listed
    /// from the database (`unrecorded_holders`).
    unrecorded: bool,
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

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:715-753
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:279-288
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262-266
// NCBI's prefetch of a query batch's reported subjects (`BlastFormatter_PreFetchSequenceData`)
// runs where `s_PreFetchSeqs` says so: never at `PRE_FETCH_SEQS_LIMIT` 0 or below one of the
// batch's Seq-aligns, always at `INT_MAX`, else for every format that `s_IsPrefetchFormat` takes:
// not the Seq-annots 8, 9 and 12, nor the archive 11. The applications but psiblast write the
// archive without calling the prefetch (blastp's, and so tblastn's, tblastx's, rpsblast's,
// blastn's and blastx's, cited at `blastinput::app::check_prefetch_environment`); psiblast's
// rounds call it before their test of the archive.
/// The first format of a run as NCBI's prefetch of each query batch asks it (`format_type`), with
/// whether the application calls the prefetch before it writes the archive (psiblast): the one
/// owner of whether a batch's prefetch runs ([`PrefetchFormat::batch_prefetches`]), which every
/// question of it takes ([`BatchPrefetch`], [`batch_loads`], [`prefetches_subjects`],
/// [`prefetch_off`], [`prefetch_may_be_off`]; review 5 of line DBID pfoff-4, M2: the batches'
/// answer from the limit and the count alone took psiblast's outfmt 8, 9 and 11 for prefetched).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PrefetchFormat {
    /// `-outfmt`'s number, as the programs give it to the formatter (the query-anchored 1 to 4
    /// as 0, DECISIONS T3-4; the tables 10 and 20 as 6).
    outfmt: i32,
    /// Whether the application calls the prefetch before it writes the archive (psiblast).
    archive_after_prefetch: bool,
}

impl PrefetchFormat {
    /// The first format `outfmt` of a run of the application `app`.
    pub(crate) fn of(outfmt: i32, app: crate::report::archive::options::RequestKind) -> Self {
        Self {
            outfmt,
            archive_after_prefetch: app == crate::report::archive::options::RequestKind::Psiblast,
        }
    }

    /// The format's number (`PrefetchFormat::of`).
    pub(crate) fn outfmt(self) -> i32 {
        self.outfmt
    }

    /// Whether a batch's prefetch reads its Seq-aligns (their subject ids and ranges,
    /// `s_ExtractSeqidsAndRanges`) where the limit `limit` allows its count: the application
    /// calls it for the format, which `s_IsPrefetchFormat` takes or which the limit `INT_MAX`
    /// takes whatever it is.
    fn reads_at(self, limit: Option<i32>) -> bool {
        use crate::blastinput::app::output_format::{
            ARCHIVE_FORMAT, ASN_BINARY, ASN_TEXT, JSON_SEQALIGN,
        };
        let called = self.outfmt != ARCHIVE_FORMAT || self.archive_after_prefetch;
        called
            && (limit == Some(i32::MAX)
                || !matches!(
                    self.outfmt,
                    ASN_TEXT | ASN_BINARY | ARCHIVE_FORMAT | JSON_SEQALIGN
                ))
    }

    /// Whether the prefetch of a batch of `seq_aligns` Seq-aligns runs at the limit `limit`
    /// (`s_PreFetchSeqs`): the one rule.
    fn runs_at(self, limit: Option<i32>, seq_aligns: usize) -> bool {
        self.reads_at(limit) && !prefetch_off_at(limit, seq_aligns)
    }

    /// Whether a batch's prefetch reads its Seq-aligns where `PRE_FETCH_SEQS_LIMIT` allows its
    /// count (`PrefetchFormat::reads_at`): the format's part of the rule.
    pub(crate) fn reads_seq_aligns(self) -> bool {
        self.reads_at(prefetch_seqs_limit())
    }

    /// Whether the prefetch of a query batch whose results hold `seq_aligns` Seq-aligns (one for
    /// each HSP of a gapped search, the count of outfmt 7's `# N hits found`, one for each subject
    /// of an ungapped one: [`subject_seq_aligns`]) runs: the one answer of every program.
    pub(crate) fn batch_prefetches(self, seq_aligns: usize) -> bool {
        self.runs_at(prefetch_seqs_limit(), seq_aligns)
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:727-753 (cited at
// `PrefetchFormat`)
/// Whether NCBI's formatters, or psiblast's PSSM, may read a batch's reported sequences one by
/// one, not all at once before its report: a batch's prefetch can be off for the first format
/// `format` ([`PrefetchFormat::batch_prefetches`]: `PRE_FETCH_SEQS_LIMIT` is an integer other than
/// `INT_MAX`, or the format is one that NCBI does not prefetch for; the readers of the headers do
/// not count the Seq-aligns, the reports of each batch do, [`BatchPrefetch`]). A value that is no
/// integer stops NCBI before its report; LOSAT refuses it in
/// `blastinput::app::check_prefetch_environment` (not for an archive that NCBI writes without the
/// prefetch).
pub(crate) fn prefetch_may_be_off(format: PrefetchFormat) -> bool {
    let limit = prefetch_seqs_limit();
    !format.reads_at(limit) || limit.is_some_and(|limit| limit != i32::MAX)
}

/// Whether NCBI's prefetch is off for a query batch whose results hold `seq_aligns` Seq-aligns
/// in a run of the first format `format` ([`PrefetchFormat::batch_prefetches`]). Its formatters
/// then fetch the sequences that each query's report shows, after the query's own warnings, where
/// the report reads them (`TabularSpec::subject_fetch`).
pub(crate) fn prefetch_off(format: PrefetchFormat, seq_aligns: usize) -> bool {
    !format.batch_prefetches(seq_aligns)
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
    /// The prefetch did not run ([`prefetch_off`]: `PRE_FETCH_SEQS_LIMIT` turned it off, or the
    /// format is one that NCBI does not prefetch for): the subjects that the reports NCBI writes
    /// read, after the formatter's cut to the hitlist size (`PruneSeqalign`), by the first
    /// format's rows (`TabularSpec::subject_fetch`; `Every` for a format that reads the subjects
    /// and is not a field list, `Never` for the Seq-annots and the archive, which read nothing of
    /// them), and those that psiblast's PSSM reads; a report that NCBI's XML stream cut before
    /// its hits loads nothing ([`XmlCuts`], [`LoadedBefore::subjects`]).
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
/// in a run of the first format `format` ([`BatchLoads`]): `Prefetched` where its prefetch runs
/// ([`PrefetchFormat::batch_prefetches`]); where it does not ([`prefetch_off`]), what its reports
/// read, `rows` (the first format's), which each program applies to the hit lists that its
/// reports show (review 3 of line DBID loadord, H2; review 6, H1: psiblast's prefetched set holds
/// its traceback beyond the reports' cut; review 5 of line DBID pfoff-4, H2: psiblast's outfmt 8,
/// 9 and 11 prefetch nothing).
pub(crate) fn batch_loads(
    format: PrefetchFormat,
    seq_aligns: usize,
    rows: crate::report::tabular::SubjectFetch,
) -> BatchLoads {
    batch_loads_at(prefetch_seqs_limit(), format, seq_aligns, rows)
}

/// [`batch_loads`] with the limit `limit` (`prefetch_seqs_limit`).
fn batch_loads_at(
    limit: Option<i32>,
    format: PrefetchFormat,
    seq_aligns: usize,
    rows: crate::report::tabular::SubjectFetch,
) -> BatchLoads {
    if format.runs_at(limit, seq_aligns) {
        BatchLoads::Prefetched
    } else {
        BatchLoads::Read(rows)
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
/// queries of each batch (`spans`), the first format, whose prefetch the batches ask
/// (`format`, [`PrefetchFormat`]), which rows of it read a subject where a batch's prefetch is
/// off (`rows`, `TabularSpec::subject_fetch`), and the outfmt 5 reports that NCBI's XML stream
/// cuts before their hits (`cuts`, the call's [`XmlCuts`], which the texts that they leave ask
/// too: [`XmlReportsBefore`]).
pub(crate) struct LoadedBefore<'a, T> {
    pub(crate) hits: &'a [T],
    pub(crate) spans: &'a [std::ops::Range<usize>],
    pub(crate) format: PrefetchFormat,
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
            format,
            rows,
            cuts,
        } = self;
        spans
            .iter()
            .flat_map(move |span| {
                let Some(held) = batch(hits, span.clone()) else {
                    return Vec::new();
                };
                let loads = batch_loads_at(limit, format, seq_aligns(span.clone(), held), rows);
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

/// Whether NCBI's prefetch is off for every query batch of a run of the first format `format`:
/// `PRE_FETCH_SEQS_LIMIT` is 0 or negative, or the format is one that NCBI does not prefetch for
/// ([`prefetch_off`] of a batch without Seq-aligns). A run whose reports then read no subject
/// looks no taxid up, so it has no refusal of NCBI's taxonomy cache to ask for.
pub(crate) fn prefetch_always_off(format: PrefetchFormat) -> bool {
    prefetch_off(format, 0)
}

/// Whether the prefetch of each query's batch reads the batch's reported subjects
/// (`BatchPrefetch::reads`), by the one rule ([`PrefetchFormat::batch_prefetches`]): NCBI counts
/// each batch's Seq-aligns against a positive `PRE_FETCH_SEQS_LIMIT` ([`prefetch_off`]), so a
/// batch above it reads nothing of its subjects before its reports, whose formatters then read
/// them one by one, as every batch of a format that NCBI does not prefetch for (psiblast's PSSM
/// reads them so in outfmt 8, 9 and 11: review 5 of line DBID pfoff-4, H2 and H3). That matters for a call whose
/// reports read nothing of the subjects ([`ParsedShown::prefetched`], refused by
/// [`check_prefetched_parsed`]) and for the weak matches of NCBI's scope, which only such
/// reads make ([`BatchPrefetch::loads_one_by_one`], walked by
/// [`crate::report::db_subject::ScopeWalk`]). A batch whose prefetch runs and fails at a subject
/// whose lookup throws on the id (`ParsedSubject::prefetch_fails`, [`batch_prefetch`]'s warning)
/// leaves its later loads to its reports, which load its subjects one by one too
/// ([`BatchPrefetch::loads_one_by_one`]; line DBID pfoff-7).
#[derive(Debug)]
pub(crate) struct BatchPrefetch(Batches, BatchStarts, BatchFails);

/// Whether the prefetch of the batch of each query of a call runs and fails at a subject whose
/// lookup throws on the id ([`BatchPrefetch`]).
#[derive(Debug)]
struct BatchFails(Vec<bool>);

/// The query batch of each query of a call ([`BatchPrefetch::report_read`]): the first query of
/// its batch, as the call's batch ranges number it; `None` for a query of no batch range.
#[derive(Debug)]
struct BatchStarts(Vec<Option<usize>>);

/// The prefetch of the batches of a call ([`BatchPrefetch`]).
#[derive(Debug)]
enum Batches {
    /// Every batch prefetches (no limit, `INT_MAX`, a subject set, a test).
    Every,
    /// No batch does (0 or below, a format that NCBI does not prefetch for).
    None,
    /// Whether the batch of each query of the call prefetches (a positive limit); `None` for a
    /// query of no batch range, which is taken both ways ([`BatchPrefetch::reads`]).
    Each(Vec<Option<bool>>),
}

impl BatchPrefetch {
    /// The batches of a call of `queries` queries (its hit lists, numbered from 0) of a run of the
    /// first format `format`: `batches`, the ranges of the queries of each batch numbered from
    /// `first` (the call's first query; ranges before it are the batches written before), or the
    /// whole call as one batch where `batches` is empty; `seq_aligns` counts the Seq-aligns of the
    /// hit lists of a range of the call ([`subject_seq_aligns`]), and `fails` tells whether the
    /// prefetch of a range would fail at one of its reported subjects whose lookup throws on the
    /// id ([`BatchPrefetch::prefetch_fails_among`]), asked of the ranges that prefetch. A query of no range is
    /// taken both as prefetched and as read one by one.
    pub(crate) fn of(
        format: PrefetchFormat,
        queries: usize,
        first: usize,
        batches: &[std::ops::Range<usize>],
        seq_aligns: impl Fn(std::ops::Range<usize>) -> usize,
        fails: impl Fn(std::ops::Range<usize>) -> bool,
    ) -> Self {
        Self::at(
            prefetch_seqs_limit(),
            format,
            queries,
            first,
            batches,
            (seq_aligns, fails),
        )
    }

    /// [`BatchPrefetch::of`] with the limit `limit` (`prefetch_seqs_limit`).
    fn at(
        limit: Option<i32>,
        format: PrefetchFormat,
        queries: usize,
        first: usize,
        batches: &[std::ops::Range<usize>],
        (seq_aligns, fails): (
            impl Fn(std::ops::Range<usize>) -> usize,
            impl Fn(std::ops::Range<usize>) -> bool,
        ),
    ) -> Self {
        let whole = [first..first + queries];
        let batches = if batches.is_empty() {
            &whole[..]
        } else {
            batches
        };
        let span = |batch: &std::ops::Range<usize>| {
            batch.start.saturating_sub(first).min(queries)
                ..batch.end.saturating_sub(first).min(queries)
        };
        let mut starts = vec![None; queries];
        for batch in batches {
            starts[span(batch)].fill(Some(batch.start));
        }
        let starts = BatchStarts(starts);
        // No batch's count decides where every count gives the same answer.
        let every = if !format.runs_at(limit, 0) {
            Batches::None
        } else if format.runs_at(limit, usize::MAX) {
            Batches::Every
        } else {
            let mut on = vec![None; queries];
            for batch in batches {
                let span = span(batch);
                if !span.is_empty() {
                    let prefetches = format.runs_at(limit, seq_aligns(span.clone()));
                    on[span].fill(Some(prefetches));
                }
            }
            Batches::Each(on)
        };
        // The batches that prefetch and whose prefetch fails.
        let mut failing = vec![false; queries];
        if !matches!(every, Batches::None) {
            for batch in batches {
                let span = span(batch);
                let prefetches = match &every {
                    Batches::Each(on) => on.get(span.start).copied().flatten() == Some(true),
                    Batches::Every | Batches::None => true,
                };
                if !span.is_empty() && prefetches && fails(span.clone()) {
                    failing[span].fill(true);
                }
            }
        }
        Self(every, starts, BatchFails(failing))
    }

    /// Whether the batch of the call's query `query` prefetches; `None` for a query of no batch.
    fn batch_of(&self, query: usize) -> Option<bool> {
        match &self.0 {
            Batches::Every => Some(true),
            Batches::None => Some(false),
            Batches::Each(on) => on.get(query).copied().flatten(),
        }
    }

    /// Whether the prefetch of the batch of the call's query `query` runs (it may fail: a lookup
    /// that throws on the id ends its loads, one that fails on a file stops the run there,
    /// [`batch_prefetch`]); not for a query of no batch.
    pub(crate) fn prefetches(&self, query: usize) -> bool {
        self.batch_of(query) == Some(true)
    }

    /// Whether the batch of the call's query `query` prefetches and its prefetch fails at a
    /// subject whose lookup throws on the id.
    fn prefetch_fails(&self, query: usize) -> bool {
        self.2 .0.get(query).copied().unwrap_or(false)
    }

    /// Whether the prefetch of a batch of the call fails at a subject whose lookup throws on the
    /// id: [`batch_prefetch`] writes its warning there, also without `taxdb.bti`.
    pub(crate) fn fails_somewhere(&self) -> bool {
        self.2 .0.contains(&true)
    }

    /// Whether the reports of the batch of the call's query `query` may load the query's subjects
    /// one by one: its prefetch is off, or the query is of no batch (taken both ways, as in
    /// [`BatchPrefetch::reads`]), or its prefetch fails at a subject whose lookup throws on the id
    /// (NCBI's prefetch catches it and ends its loads; the reports load the subjects it did not).
    /// The one answer of every reader of the reports' loads (review 1 of line DBID pf7, H1 and
    /// M2): where a subject's load throws (`ParsedSubject::lookup_throws`), where NCBI's scope can
    /// match an id weakly ([`BatchPrefetch::report_read`], [`BatchPrefetch::reads`]), the texts of
    /// the reports' reads (each program's report model, after the batch's prefetch where it runs:
    /// the loads that a failed prefetch never made write theirs there) and psiblast's PSSM reads.
    pub(crate) fn loads_one_by_one(&self, query: usize) -> bool {
        !self.prefetches(query) || self.prefetch_fails(query)
    }

    /// How the report of the call's query `query` reads its subjects, for the walk of the
    /// reports' loads (`report::db_subject::HitLoads::of`): one by one or not
    /// ([`BatchPrefetch::loads_one_by_one`]), after a prefetch that fails or not
    /// (`ReportRead::prefetch_fails`), in the scope of its query batch, named by the batch's
    /// first query as the call's batch ranges number it plus `base` (the number of the call's
    /// first batch range in the run's record, `TaxdbRun::loads_before`; `None` for a query of no
    /// batch range).
    pub(crate) fn report_read(
        &self,
        query: usize,
        base: usize,
    ) -> crate::report::db_subject::ReportRead {
        crate::report::db_subject::ReportRead {
            one_by_one: self.loads_one_by_one(query),
            batch: self
                .1
                 .0
                .get(query)
                .copied()
                .flatten()
                .map(|start| start + base),
            prefetch_fails: self.prefetch_fails(query),
        }
    }

    /// What the batches read of the subjects of `reported`, the (query of the call, ordinal id)
    /// of each reported subject of each query: the subjects that a batch whose prefetch is on
    /// reports, and those that a batch whose reports load them one by one reports (its prefetch
    /// is off or fails, [`BatchPrefetch::loads_one_by_one`]; `BatchReads::shown`). The subjects are
    /// listed only where a positive limit counts the batches' Seq-aligns or a prefetch fails.
    pub(crate) fn reads(&self, reported: impl IntoIterator<Item = (usize, u32)>) -> BatchReads {
        let each = match &self.0 {
            Batches::Every if !self.fails_somewhere() => return BatchReads(BatchesRead::Every),
            Batches::None => return BatchReads(BatchesRead::None),
            Batches::Every | Batches::Each(_) => {
                let mut on = std::collections::HashSet::new();
                for (query, ordinal) in reported {
                    if self.batch_of(query) != Some(false) {
                        on.insert(ordinal);
                    }
                }
                BatchesRead::Each { on }
            }
        };
        BatchReads(each)
    }

    /// Whether the prefetch of a batch whose reported subjects are `subjects` fails at one whose
    /// lookup throws on the id ([`batch_prefetch`]'s warning; `ParsedSubject::prefetch_fails`):
    /// the `fails` of [`BatchPrefetch::of`] of each program.
    pub(crate) fn prefetch_fails_among<'s>(
        subjects: impl IntoIterator<Item = Option<&'s crate::report::db_subject::ParsedSubject>>,
    ) -> bool {
        subjects
            .into_iter()
            .flatten()
            .any(crate::report::db_subject::ParsedSubject::prefetch_fails)
    }
}

/// What the batches of a call read of each subject they report ([`BatchPrefetch::reads`]).
#[derive(Debug)]
pub(crate) struct BatchReads(BatchesRead);

#[derive(Debug)]
enum BatchesRead {
    Every,
    None,
    /// The ordinal ids that a batch whose prefetch is on reports.
    Each {
        on: std::collections::HashSet<u32>,
    },
}

impl BatchReads {
    /// Whether a batch whose prefetch is on reports the subject `ordinal`.
    pub(crate) fn prefetched(&self, ordinal: u32) -> bool {
        match &self.0 {
            BatchesRead::Every => true,
            BatchesRead::None => false,
            BatchesRead::Each { on } => on.contains(&ordinal),
        }
    }

    /// What the reports `shown` show of the reported subject `ordinal` ([`ParsedShown::of_batches`]).
    pub(crate) fn shown(&self, shown: ParsedShown, ordinal: u32) -> ParsedShown {
        shown.of_batches(self.prefetched(ordinal))
    }
}

/// The limit of `PRE_FETCH_SEQS_LIMIT` (`app::prefetch_seqs_limit`); `None` when it is not set. A
/// value that NCBI cannot convert stops the first batch before its prefetch
/// (`app::check_prefetch_environment`), so no batch's prefetch is decided by it: the plans made
/// before the search count it as unset.
fn prefetch_seqs_limit() -> Option<i32> {
    crate::blastinput::app::prefetch_seqs_limit().ok().flatten()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1328-1372, 1439
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:58, 75-110
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:257-290, 441-462
// The BLAST database data loader loads a sequence for the id it is asked for (the blob is the OID
// with that id's handle) and asks its volume for the Bioseq with that id as the target: the
// Bioseq's ids are those of the first definition line that holds the target, so NCBI's scope
// indexes that line's ids alone; the loader registers in `m_Ids` the handle of every Seq-id of
// every definition line of the set (`RegisterIds` of the whole ASN.1 set that the Bioseq carries).
// Every load of a reported subject in NCBI's searches is a load by its Seq-align id (the reports,
// the batches' prefetch, psiblast's PSSM reads; `GetBioseq` of an OID without a target, whose
// Bioseq takes the first line's ids, puts nothing into a scope there), so the scope knows a
// sequence by the line of its own Seq-align id (review 3 of line DBID pfoff-4, M1). A load of a
// holder that the loader gives for a failing subject's id targets the line of that id; LOSAT
// refuses every such hold (`ParsedSubject::fetched`), so nothing after it is written.
/// The definition lines of a sequence in which one of NCBI's two places that answer a load before
/// its lookup finds its ids (`report::db_subject::ScopeHolders`): the one owner of that choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Lines {
    /// Every line: the BLAST database data loader's `m_Ids`, by exact handle.
    Every,
    /// The line by which a load by the sequence's own Seq-align id makes its Bioseq
    /// (`DeflineSet::bioseq_defline` of `DeflineSet::alignment_id`; none where no line holds it,
    /// where that load throws): NCBI's scope, by every handle that matches the id it is asked
    /// for.
    Loaded,
}

impl Lines {
    /// The Seq-ids of these lines of the set `set`, in order.
    fn ids(
        self,
        set: &crate::seqid::DeflineSet,
    ) -> impl Iterator<Item = &crate::seqid::ParsedSeqId> {
        let loaded = match self {
            Lines::Every => None,
            Lines::Loaded => Some(set.alignment_id().and_then(|id| set.bioseq_defline(id))),
        };
        set.deflines
            .iter()
            .enumerate()
            .filter(move |&(at, _)| loaded.is_none_or(|line| line == Some(at)))
            .flat_map(|(_, defline)| &defline.ids)
    }
}

/// Whether a Seq-id of the lines `lines` of the sequence `ordinal` of `db` is one `CSeq_id_Handle`
/// with the key `handle` (`seqid::mask_handle::mask_handle_key`).
fn holds_handle(db: &SeqDb, ordinal: u32, handle: &[u8], lines: Lines) -> bool {
    let Ok(crate::seqdb::SequenceHeader::Parsed(set)) = db.read_header(ordinal) else {
        return false;
    };
    let held = lines.ids(&set).any(|id| {
        crate::seqid::mask_handle::mask_handle_key(id)
            .is_ok_and(|key| key.as_deref() == Some(handle))
    });
    held
}

/// The words of a lookup's error `error` in LOSAT's own refusal texts: NCBI's message where it
/// has one.
fn error_words(error: &SeqDbError) -> String {
    error.ncbi_message().map_or_else(
        || error.to_string(),
        |text| String::from_utf8_lossy(&text).into_owned(),
    )
}

/// The sequences that the searches of [`HeaderReader::case_holders`] find for a Seq-align id, and
/// what tells which hold an id that NCBI's scope takes for it.
struct CaseHolders<'r, 'd> {
    reader: &'r HeaderReader<'d>,
    folded: Vec<u32>,
    variants: Vec<u32>,
    /// The scope keys of the Seq-align id (`defline_set::scope_keys`).
    keys: Vec<u64>,
    /// The key of the Seq-align id's own handle, where a key in another letter case found one.
    handle: Option<Vec<u8>>,
}

impl<'r, 'd> CaseHolders<'r, 'd> {
    /// The OIDs, of `oids` (every OID of the lookup of the Seq-align id) and those found in
    /// another letter case, of the sequences whose lines `lines` hold an id that NCBI's scope
    /// takes for it: a scope key of it, or for one found under a key in another letter case, an
    /// id of its own handle. In their order, each as often as found.
    fn holders<'s>(
        &'s self,
        oids: &'s [u32],
        lines: Lines,
    ) -> impl Iterator<Item = u32> + use<'s, 'r, 'd> {
        let reader = self.reader;
        oids.iter()
            .chain(&self.folded)
            .map(|&oid| (oid, false))
            .chain(self.variants.iter().map(|&oid| (oid, true)))
            .filter(move |&(oid, variant)| {
                let other = reader.db.ordinal_of(oid);
                if variant {
                    self.handle
                        .as_deref()
                        .is_some_and(|handle| holds_handle(reader.db, other, handle, lines))
                } else {
                    reader.holds_key(
                        other,
                        &self.keys,
                        crate::seqid::defline_set::scope_keys,
                        lines,
                    )
                }
            })
            .map(|(oid, _)| oid)
    }
}

/// What [`HeaderReader::weak_matches`] finds for a Seq-align id.
enum WeakMatches {
    /// It matches nothing weakly (no weak keys).
    None,
    /// A version 5 database cannot list what it matches weakly (a name, a `pir` or `prf` id).
    Unlisted,
    /// Its weak keys (`id_scope::weak_keys`) and the OIDs of the sequences that may hold one.
    Found { keys: Vec<u64>, oids: Vec<u32> },
}

impl WeakMatches {
    /// The OIDs found whose lines `lines` hold an id of a weak key (`id_scope::kind_keys`), in
    /// their order, each as often as found.
    fn holders<'s, 'a>(
        &'s self,
        reader: &'s HeaderReader<'a>,
        lines: Lines,
    ) -> impl Iterator<Item = u32> + use<'s, 'a> {
        let (keys, oids): (&[u64], &[u32]) = match self {
            WeakMatches::Found { keys, oids } => (keys, oids),
            WeakMatches::None | WeakMatches::Unlisted => (&[], &[]),
        };
        oids.iter().copied().filter(move |&oid| {
            reader.holds_key(
                reader.db.ordinal_of(oid),
                keys,
                crate::seqid::id_scope::kind_keys,
                lines,
            )
        })
    }
}

impl<'a> HeaderReader<'a> {
    /// The reader of the headers of `db` for a call after calls of the run whose loads no record
    /// holds where `unrecorded` (`TaxdbRun::unrecorded_calls`).
    pub(crate) fn new(db: &'a SeqDb, unrecorded: bool) -> Self {
        Self {
            db,
            utf8: Default::default(),
            batches: std::cell::OnceCell::new(),
            lookup: std::cell::OnceCell::new(),
            indexed_lookup: std::cell::OnceCell::new(),
            unrecorded,
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:903-916
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:389-406
    // The data loader resolves a Seq-align id that no place before it answers by NCBI's lookup in
    // the search's own database, list included (`SeqidToOid`, the first OID found), and loads that
    // sequence's blob (step 4 of `report::db_subject::ScopeWalk::load`). Which place answers a load
    // is the walk's; the lookup's answer is the subject's own, made once here. Another sequence
    // that holds the same definition lines and letters is shown as the subject itself (the same
    // record in two volumes or two databases of an alias), one that differs only in their taxids
    // as the subject but by the reports that write taxids (`FoundInstead::OtherTaxids`, `copy_of`)
    // and, where a title can be made from them, those of a generated title (`CopyApart::Title`).
    /// What NCBI's lookup of the Seq-align id of the sequence `ordinal`, `subject`, gives
    /// ([`crate::report::db_subject::ScopeResolution`], the one producer): its answer where that
    /// is not the sequence (`found`), the sequence it gives first where that is another one
    /// (`gives`), the header with a UTF8String where a VisibleString is expected that it reads
    /// (`utf8`: its version check of an accession, `lookup_utf8`, or the sequence it gives first,
    /// `scope_read`), and in an API call after the first, whose calls before LOSAT has no record
    /// of (`unrecorded`), the other sequences that NCBI's places may give for it
    /// (`unrecorded_holders`). A lookup that fails is kept for the formats that read the Bioseq;
    /// the others never look the id up. A subject whose Bioseq NCBI's formatters do not find
    /// (`unfetchable`) is looked up too, for psiblast's PSSM read (its reports refuse it whatever
    /// the lookup gives). A giim id of a version 5 database is looked up only for a lookup that
    /// fails: NCBI's handles never take another sequence for it, but its prefetch looks it up as a
    /// string id in every LMDB file of the list, which stops at a file without the accession table
    /// (review 1 of line DBID mixdb, H1); of a version 4 database never.
    fn resolution(
        &self,
        subject: &crate::report::db_subject::ParsedSubject,
        ordinal: u32,
    ) -> crate::report::db_subject::ScopeResolution {
        use crate::report::db_subject::{FoundInstead, ScopeResolution};
        let giim = matches!(subject.alignment_id(), crate::seqid::ParsedSeqId::Giim(_));
        let mut resolution = ScopeResolution::default();
        if giim && !self.db.is_version5() {
            return resolution;
        }
        let failed = |resolution: ScopeResolution, error: SeqDbError| ScopeResolution {
            found: Some(FoundInstead::Failed(
                crate::report::db_subject::LookupFailure::of_error(&error),
            )),
            unrecorded: self.unrecorded_holders(subject, ordinal, None),
            ..resolution
        };
        let lookup = match self.lookup.get() {
            Some(lookup) => lookup,
            None => match self.db.seqid_lookup() {
                Ok(lookup) => {
                    let _ = self.lookup.set(lookup);
                    self.lookup.get().expect("set above")
                }
                Err(error) => return failed(resolution, error),
            },
        };
        let id = subject.alignment_id();
        // NCBI's lookup (`SeqidToOid`), whose version check reads its candidates' headers, also
        // where it fails after them (`lookup_utf8`): the load of the id reads them first.
        let found = lookup.seqid_to_oid(id);
        resolution.utf8 = self.lookup_utf8(lookup);
        let first = match found {
            Ok(_) if giim => return resolution,
            Ok(first) => first,
            // A lookup that throws on the id (`SeqidLookup::throws_on_id`) loads nothing where the
            // formatters load the subject one by one, as one that fails on a file: the holders
            // that NCBI's scope may give for it without the lookup are listed alike.
            Err(error) => {
                return match lookup.throws_on_id(id) {
                    Some(exception) => ScopeResolution {
                        found: Some(FoundInstead::ThrowsOnId(exception)),
                        unrecorded: self.unrecorded_holders(subject, ordinal, None),
                        ..resolution
                    },
                    None => failed(resolution, error),
                }
            }
        };
        let own = self.db.oid(ordinal).get();
        let other = first
            .filter(|&first| first != own)
            .map(|first| self.db.ordinal_of(first));
        // The data loader reads the header of the sequence that this lookup gives first.
        if let Some(other) = other.filter(|_| resolution.utf8.is_none()) {
            resolution.utf8 = self.scope_read(other);
        }
        resolution.gives = other.map(|other| self.db.sequence(other));
        resolution.found = match other {
            None => first.is_none().then_some(FoundInstead::Nothing),
            Some(other) => {
                let sequence = self.db.sequence(other);
                match self.copy_of(subject, ordinal, other) {
                    RecordCopy::Same => None,
                    RecordCopy::OtherTaxids => Some(FoundInstead::OtherTaxids(sequence)),
                    RecordCopy::Apart(cause) => Some(FoundInstead::CopyApart(sequence, cause)),
                    RecordCopy::Differs => {
                        // Its blob registers its lines where `x_GetOid` keeps it.
                        if let Ok(crate::seqdb::SequenceHeader::Parsed(set)) =
                            self.db.read_header(other)
                        {
                            resolution.gives_keys = loader_keeps(&set, subject)
                                .then(|| subject.shown_keys(&set, sequence));
                        }
                        Some(FoundInstead::Sequence(sequence))
                    }
                }
            }
        };
        if self.unrecorded {
            resolution.unrecorded =
                self.unrecorded_holders(subject, ordinal, Some((lookup, first)));
        }
        resolution
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-396
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:864-880
    // NCBI's scope finds a Bioseq it holds for the Seq-align id before the data loader looks it
    // up, by the id's handle, an id that it takes for it, or a weak match, among the ids of the
    // definition line that each Bioseq was loaded by (`Lines::Loaded`), and the loader gives a
    // sequence it loaded for every Seq-id of every line of it (`RegisterIds`), by the handle
    // (reviews 1 to 3 of line DBID pfoff-4, H1 and H2, H1 and M1, M1). In an API call after the
    // first, LOSAT has no record of what the calls before loaded, so it lists the sequences that
    // either place may give from the database.
    /// For `resolution` in an API call after the first ([`HeaderReader::new`]'s `unrecorded`):
    /// the other sequences of the database that NCBI may give for the Seq-align id of the sequence
    /// `ordinal`, `subject`, without the lookup, where the calls before loaded them
    /// ([`crate::report::db_subject::ScopeHolders`]). `lookup`: the lookup that does not fail,
    /// with the OID that it gives first (`None`: the lookup fails on a file, and the holders are
    /// found in the LMDB files that have the accession table, `SeqidLookup::in_indexed_files`: a
    /// file without it holds no text id, and NCBI's lookup fails there, as the stop of a prefetch
    /// finds them, [`stop_lookup_unsure`]). The searches are a lookup's (every OID of the lookup and
    /// the keys in another letter case, [`HeaderReader::case_holders`], with the scope keys): the
    /// scope's, those whose definition line of their own Seq-align id (`Lines::Loaded`) holds a
    /// scope key of it, with those of its weak matches ([`HeaderReader::weak_matches`]) by that
    /// line; the data loader's, those whose lines (any, `Lines::Every`) hold an id of the Seq-align
    /// id's handle (`seqid::mask_handle::mask_handle_key`; every one where the id has no handle
    /// that a function can know). Each list keeps the sequences whose own Seq-align id's lookup
    /// does not fail (a load of one by that id would throw and load nothing), but the subject and
    /// the sequence that its lookup gives first (a load of either is the subject's own).
    /// `Holders::Unknown` where a list fails on a file, and for the scope's where the weak matches
    /// cannot be listed.
    fn unrecorded_holders(
        &self,
        subject: &crate::report::db_subject::ParsedSubject,
        ordinal: u32,
        lookup: Option<(&crate::seqdb::lookup::SeqidLookup<'a>, Option<u32>)>,
    ) -> crate::report::db_subject::ScopeHolders {
        use crate::report::db_subject::{Holders, ScopeHolders};
        if !self.unrecorded {
            return ScopeHolders::default();
        }
        let unknown = |error: &str| {
            let why = format!("the other sequences that hold it, which a call of the run before this one may have loaded ({error})");
            ScopeHolders {
                loader: Holders::Unknown(why.clone()),
                scope: Holders::Unknown(why),
            }
        };
        let (searched, first) = match lookup {
            Some((lookup, first)) => (lookup, first),
            None => {
                let indexed = self.indexed_lookup.get_or_init(|| {
                    self.db
                        .seqid_lookup()
                        .map(crate::seqdb::lookup::SeqidLookup::in_indexed_files)
                        .map_err(|error| error_words(&error))
                });
                match indexed {
                    Ok(indexed) => (indexed, None),
                    Err(error) => return unknown(error),
                }
            }
        };
        let id = subject.alignment_id();
        // LOSAT's own scans, not reads of NCBI's lookup: their version checks note nothing. A
        // Seq-align id whose lookup throws on the id: the sequences under the keys that
        // makeblastdb writes for it (`SeqidLookup::thrown_id_oids`), which `list` below keeps to
        // those whose lines hold it.
        let oids = if searched.throws_on_id(id).is_some() {
            searched.thrown_id_oids(id)
        } else {
            searched.seqid_to_oids(id, true)
        };
        let listed = oids.and_then(|oids| {
            let weak = self.weak_matches(searched, subject, &oids)?;
            Ok((self.case_holders(searched, subject)?, oids, weak))
        });
        let _ = searched.take_utf8();
        let (in_case, oids, weak) = match listed {
            Ok(listed) => listed,
            Err(error) => return unknown(&error_words(&error)),
        };
        let own = self.db.oid(ordinal).get();
        let lookup = self.lookup.get();
        let list = |found: &mut dyn Iterator<Item = u32>| {
            let mut holders: Vec<u32> = found
                .filter(|&oid| oid != own && Some(oid) != first)
                .map(|oid| self.db.ordinal_of(oid))
                .filter(|&other| {
                    let Ok(crate::seqdb::SequenceHeader::Parsed(set)) = self.db.read_header(other)
                    else {
                        return true;
                    };
                    set.alignment_id().is_none_or(|id| {
                        lookup.is_some_and(|lookup| lookup.seqid_to_oid(id).is_ok())
                    })
                })
                .collect();
            if let Some(lookup) = lookup {
                let _ = lookup.take_utf8();
            }
            holders.sort_unstable();
            holders.dedup();
            holders
                .into_iter()
                .map(|other| self.db.sequence(other))
                .collect::<Vec<_>>()
        };
        let handle = crate::seqid::mask_handle::mask_handle_key(id)
            .ok()
            .flatten();
        // The loader knows every line of a sequence it loaded, the scope the line of its load
        // (`Lines`; review 3 of line DBID pfoff-4, M1).
        let loader = list(&mut in_case.holders(&oids, Lines::Every).filter(|&oid| {
            handle.as_deref().is_none_or(|handle| {
                holds_handle(self.db, self.db.ordinal_of(oid), handle, Lines::Every)
            })
        }));
        let scope = match weak {
            WeakMatches::Unlisted => Holders::Unknown(
                "the sequences that it matches weakly by its name, which a version 5 database cannot list"
                    .to_owned(),
            ),
            weak => Holders::Listed(list(
                &mut in_case
                    .holders(&oids, Lines::Loaded)
                    .chain(weak.holders(self, Lines::Loaded)),
            )),
        };
        ScopeHolders {
            loader: Holders::Listed(loader),
            scope,
        }
    }

    /// The sequences that `lookup` finds for the Seq-align id of `subject` under its keys in
    /// another letter case, LOSAT's own searches (version 5: the upper-case key of an id without
    /// a text accession, the other keys of a name or a `pir` or `prf` id), with what tells which
    /// of them and of the OIDs of its lookup hold an id that NCBI's scope takes for it
    /// ([`CaseHolders::holders`]): where one of these searches fails, NCBI's exact lookup, which
    /// read the files before, still succeeds, so the holders are unknown (review 1 of line DBID
    /// pircase, L1; review 1 of line DBID v4late, L2).
    fn case_holders<'r>(
        &'r self,
        lookup: &crate::seqdb::lookup::SeqidLookup<'_>,
        subject: &crate::report::db_subject::ParsedSubject,
    ) -> std::result::Result<CaseHolders<'r, 'a>, SeqDbError> {
        let id = subject.alignment_id();
        let folded = lookup.case_folded_oids(id)?;
        let variants = lookup.case_variant_oids(id)?;
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
        Ok(CaseHolders {
            reader: self,
            folded,
            variants,
            keys,
            handle,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:128-135
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:160-163
    // A version 5 database keys a text Seq-id as written, in its letter case, with no key for the
    // name of an id with an accession, and a `pir` or `prf` id by its FASTA text; the scope matches
    // names without case and across the alternatives (`CSeq_id_Mapper::GetMatchingHandles`), so the
    // lookup cannot list the holders of a name (review 1 of line DBID N2c,
    // `blastdb-ids-20261021/R/weak/REVIEW_weak_1.md` H1). A version 4 database's ISAM keys are in
    // lower case without the alternative.
    /// The sequences of the database that hold an id that the Seq-align id of `subject` matches
    /// weakly (`seqid::id_scope::weak_keys`), as `lookup` finds them: among the sequences that its
    /// lookup gives (`oids`: in both versions an accession without a version gives its versions,
    /// in version 4 a name its holders) and that the lookup of its name alone gives. In version 5
    /// a Seq-align id with a name, or a `pir` or `prf` one, with weak keys cannot be listed
    /// (`WeakMatches::Unlisted`). The one search of the weak matches of NCBI's scope that LOSAT
    /// makes in the database (`unrecorded_holders`).
    fn weak_matches(
        &self,
        lookup: &crate::seqdb::lookup::SeqidLookup<'_>,
        subject: &crate::report::db_subject::ParsedSubject,
        oids: &[u32],
    ) -> Result<WeakMatches, SeqDbError> {
        let id = subject.alignment_id();
        let Some(text) = id.textseq() else {
            return Ok(WeakMatches::None);
        };
        let mut keys = Vec::new();
        crate::seqid::id_scope::weak_keys(id, &mut keys);
        if keys.is_empty() {
            return Ok(WeakMatches::None);
        }
        let pir_or_prf = matches!(
            id,
            crate::seqid::model::ParsedSeqId::Pir(_) | crate::seqid::model::ParsedSeqId::Prf(_)
        );
        if self.db.is_version5() && (text.name.is_some() || pir_or_prf) {
            return Ok(WeakMatches::Unlisted);
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
        Ok(WeakMatches::Found { keys, oids: found })
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

    /// Whether a Seq-id of the lines `lines` of the sequence `ordinal` has one of the keys `keys`
    /// that `key_of` gives (`defline_set::scope_keys`: NCBI's scope can show it for a Seq-align id
    /// of those keys; `id_scope::kind_keys`: for a Seq-align id whose `weak_keys` they are).
    fn holds_key(
        &self,
        ordinal: u32,
        keys: &[u64],
        key_of: fn(&crate::seqid::model::ParsedSeqId, &mut Vec<u64>),
        lines: Lines,
    ) -> bool {
        let Ok(crate::seqdb::SequenceHeader::Parsed(set)) = self.db.read_header(ordinal) else {
            return false;
        };
        let mut held = Vec::new();
        let holds = lines.ids(&set).any(|id| {
            held.clear();
            key_of(id, &mut held);
            held.iter().any(|key| keys.contains(key))
        });
        holds
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
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1380-1386
        // No NCBI counterpart: LOSAT's copy-equivalence check (`copy_of`). `GetBioseq` returns a
        // null Bioseq when the signed 32-bit `GetSeqLengthProt` length is less than 1, so treat
        // such proteins as distinct without loading their sequence bytes into memory.
        let protein = self.db.seqtype() == SeqType::Protein;
        if protein {
            let (Ok(len_a), Ok(len_b)) =
                (self.db.protein_length_i32(a), self.db.protein_length_i32(b))
            else {
                return RecordCopy::Differs;
            };
            if len_a < 1 || len_b < 1 {
                return RecordCopy::Differs;
            }
        }
        let (Ok(block_a), Ok(block_b)) = (self.db.read_sequence(a), self.db.read_sequence(b))
        else {
            return RecordCopy::Differs;
        };
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
                // A copy whose title holds a byte read as '#' is taken to differ: what the reports
                // read of it posts (`check_reported_title`).
                if copy.replaced_title_bytes() {
                    return RecordCopy::Differs;
                }
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
                    if subject.replaced_title_bytes() {
                        db.note_replaced_title();
                    }
                    let resolution = headers.resolution(&subject, ordinal);
                    ReportedHeader::Parsed(std::sync::Arc::new(
                        subject
                            .with_source(db.seqtype() == SeqType::Protein, |taxid| {
                                headers.source_taxname(Some(taxid))
                            })
                            .showing_gis(db.reports_show_gis())
                            .under_positive_list(db.has_positive_list())
                            .resolved(resolution)
                            .in_database(db.sequence(ordinal)),
                    ))
                })
        }
        Ok(crate::seqdb::SequenceHeader::Ordinal(deflines)) => {
            let mut stored = stored_title(Ok(deflines), || db.oid(ordinal), query, program)?;
            if stored.replaced_bytes {
                db.note_replaced_title();
            }
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
                (
                    headers
                        .first_utf8_before(Utf8Phase::Reported, query)
                        .as_ref(),
                    Some(ordinal),
                ),
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
/// `oid` gives the sequence's OID, which the refusal names. The callers ask only where the
/// reports resolve the sequence's Seq-align id (`ParsedShown::resolves_subject`: a report reads
/// its Bioseq, or the batches' prefetch loads it); the Seq-annots, the archive and a table of
/// fields that read nothing of the subject, without the prefetch, write NCBI's bytes (review 2 of
/// line RG-KB item 12, N1). A PDB id's key is its handle's (`defline_set::scope_keys`), so the
/// query `pdb|1AB1|A` takes the place of a sequence whose Seq-align id is `pdb` 1AB1 with the chain
/// 321 or 65 and no chain-id.
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
    anchored: bool,
    program: &str,
) -> Result<()> {
    if parsed.is_empty() || (!anchored && generated) {
        return Ok(());
    }
    let parsed_of = |rank: u32| parsed.get(rank as usize).and_then(Option::as_deref);
    for query in 0..num_queries {
        let scope = scope_of(query);
        let ranks = reported(query);
        if anchored && generated {
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
/// whose header `titles` holds with parsed Seq-ids, with the scope of that query, where the
/// reports `shown` resolve the sequences' Seq-align ids (`ParsedShown::resolves_subject`), unless
/// the scope finds a query of its batch for the sequence and the reports show NCBI's bytes for it:
/// that query is noted in `seen` (`report::query_bioseq::note_or_refuse` with `run`). A program
/// calls it for each query batch after it reads the batch's titles: its cache of the titles keeps
/// those of the batches before.
#[allow(clippy::too_many_arguments)]
pub(crate) fn refuse_reported_ids_of_queries(
    titles: &std::collections::BTreeMap<u32, ReportedHeader>,
    reported: impl IntoIterator<Item = (usize, u32)>,
    oids: crate::seqdb::OidMap<'_>,
    scopes: &crate::blastinput::query_batch::QueryScopes,
    shown: ParsedShown,
    seen: &mut crate::report::query_bioseq::QueryBioseqs,
    run: &crate::report::query_bioseq::SeenRun<'_, '_>,
    program: &str,
) -> Result<()> {
    if !shown.resolves_subject() {
        return Ok(());
    }
    for (query, ordinal) in reported {
        if let Some(ReportedHeader::Parsed(subject)) = titles.get(&ordinal) {
            crate::report::query_bioseq::note_or_refuse(
                seen,
                subject,
                || oids.oid(ordinal),
                query,
                ordinal,
                scopes,
                run,
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
/// (`description_line_refusal`).
const LONG_DESCRIPTION_ID: &str = "has an outfmt 0 description id longer than 65 bytes (NCBI BLAST+ 2.17.0 pads the description line without end)";

/// LOSAT's refusals of the outfmt 0 description and alignment heading of the subject with
/// parsed Seq-ids `subject` (`oid` gives its OID, which a refusal names; `protein` for a protein
/// database)
/// that it does not write as NCBI: a protein title of punctuation that NCBI's `x_CleanAndCompress`
/// reads past the end of (NCBI BLAST+ 2.17.0 crashes or does not end; approved exception 2 of
/// PD-LOSAT-NCBI-DEFECTS has BLASTN, TBLASTN and TBLASTX stop such a nucleotide title at the end
/// of the string, as for a database without parsed Seq-ids), a title
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
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1586-1600
    // Check the protein description title for punctuation-only titles that NCBI's
    // `x_CleanAndCompress` reads past the end of; the `> 65`-byte description-id check lives in
    // `description_line_refusal` so it runs only where a report prints the description line.
    if protein {
        let (_, title, _) = subject.description_line(protein);
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

/// The refusal of a description line whose `-show_gis` GI NCBI never sets
/// (`ParsedSubject::description_gi_unset`, `description_line_refusal`).
const UNSET_DESCRIPTION_GI: &str = "has a Seq-align id by which NCBI BLAST+'s scope gets no Bioseq, whose outfmt 0 to 4 description line with -show_gis prints a GI that NCBI BLAST+ 2.17.0 never sets (a different number in each run)";

/// [`check_shown_parsed_subject`] of a subject whose load in the pairwise reports is `load`
/// (`HitLoads::load`): one whose load throws, or whose Bioseq the scope does not get
/// (`ParsedSubject::scope_finds_none`), shows only its Seq-align id and `Unknown` in its
/// description line (`ParsedSubject::thrown_description`), with no heading, so that nothing is
/// refused here. What its description line shows is refused only where a report prints that line
/// (`refuse_described_lines`).
pub(crate) fn check_pairwise_parsed(
    subject: &crate::report::db_subject::ParsedSubject,
    load: crate::report::db_subject::HitLoad,
    protein: bool,
    taxdb_found: bool,
    oid: impl Fn() -> crate::seqdb::Oid,
    program: &str,
) -> Result<()> {
    if load != crate::report::db_subject::HitLoad::Throws && !subject.scope_finds_none() {
        return check_shown_parsed_subject(subject, protein, taxdb_found, oid, program);
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:913-917
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:886-892
/// The reason of LOSAT's refusal of the outfmt 0 to 4 description line of the subject `subject`
/// whose load in the pairwise reports is `load` (`HitLoads::load`), `None` where LOSAT writes it
/// as NCBI: a description id longer than 65 bytes (`description_start`, or of a subject whose
/// load throws or whose Bioseq the scope does not get, `thrown_description`: NCBI's unsigned
/// length wraps and NCBI BLAST+ 2.17.0 pads the line with spaces without end; probe of the ID2a
/// review, `tests/db_search/ids/refuse.faa`), and the `-show_gis` GI of a subject whose Bioseq the
/// scope does not get, which NCBI never sets (`ParsedSubject::description_gi_unset`). The title of
/// a Bioseq that NCBI's data loader makes of another definition line than the Seq-align id's is
/// that line's, as LOSAT writes it (`ParsedSubject::target`, line RG-KB item 12 review 2).
pub(crate) fn description_line_refusal(
    subject: &crate::report::db_subject::ParsedSubject,
    load: crate::report::db_subject::HitLoad,
    protein: bool,
) -> Option<String> {
    let unloaded = load == crate::report::db_subject::HitLoad::Throws || subject.scope_finds_none();
    if unloaded && subject.description_gi_unset() {
        return Some(UNSET_DESCRIPTION_GI.to_owned());
    }
    let id_length = if unloaded {
        subject.thrown_description().0.len()
    } else if subject.description_may_throw() {
        subject.description_line(protein).0.len()
    } else {
        subject.description_start().len()
    };
    (id_length > 65).then(|| LONG_DESCRIPTION_ID.to_owned())
}

/// LOSAT's refusals of the outfmt 0 to 4 description lines that the reports of a call print
/// (`description_line_refusal`), before any byte of the report that is not NCBI's own prefix (a
/// program that writes its reports batch by batch asks at each batch, after the prolog and the
/// earlier batches, which NCBI writes as LOSAT does before it pads the line): for each query
/// `0..num_queries` whose final hit list (`reported`, by rank) holds a sequence with parsed Seq-ids
/// (`parsed`, by rank; `load_of` its load in the pairwise reports, `oid_of` its OID) whose line
/// LOSAT does not write as NCBI, the sequences whose description lines its report shows
/// (`described`, by rank: `pairwise::blastp_described_subjects`, the first `-num_descriptions`); a
/// sequence below `-num_descriptions` prints no line and is not refused (review 1 of line RG-KB
/// item 12, R1). Nothing is asked where no sequence has such a line. The `-subject` records go
/// through the same refusals (`refuse_described_records`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn refuse_described_lines(
    num_queries: usize,
    reported: impl FnMut(usize) -> Vec<u32>,
    described: impl FnMut(usize) -> Result<Vec<u32>>,
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    load_of: impl Fn(u32) -> crate::report::db_subject::HitLoad,
    oid_of: impl Fn(u32) -> crate::seqdb::Oid,
    program: &str,
) -> Result<()> {
    let protein = matches!(program, "BLASTP" | "PSIBLAST" | "BLASTX" | "RPSBLAST");
    match described_refusal(num_queries, reported, described, parsed, load_of, protein)? {
        Some((rank, reason)) => Err(parsed_refusal(oid_of(rank), &reason, program)),
        None => Ok(()),
    }
}

/// `refuse_described_lines` of the `-subject` records of a run under `-parse_deflines`
/// (`records`, `db_subject::parsed_records`, by subject index: each is its own Bioseq), whose
/// refusal names the record (`ParsedSubject::record_refusal`); `reported` and `described` give
/// subject indexes. One owner with the database's sequences: a record whose description line no
/// report prints is not refused (review 2 of line RG-KB item 12, N2).
pub(crate) fn refuse_described_records(
    num_queries: usize,
    reported: impl FnMut(usize) -> Vec<u32>,
    described: impl FnMut(usize) -> Result<Vec<u32>>,
    records: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    program: &str,
) -> Result<()> {
    // A record that the reports show as one read without `-parse_deflines` has its local id there
    // (`ParsedSubject::shown_as_record`).
    let shown: Vec<_> = records
        .iter()
        .map(|record| record.clone().filter(|record| !record.shown_as_record()))
        .collect();
    let loaded = |_| crate::report::db_subject::HitLoad::Loaded;
    let protein = matches!(program, "BLASTP" | "PSIBLAST" | "BLASTX" | "RPSBLAST");
    match described_refusal(num_queries, reported, described, &shown, loaded, protein)? {
        Some((index, reason)) => Err(shown[index as usize]
            .as_deref()
            .expect("a refused record has parsed ids")
            .record_refusal(&reason, program)),
        None => Ok(()),
    }
}

/// The first refused description line of `refuse_described_lines`: its rank and the reason.
fn described_refusal(
    num_queries: usize,
    mut reported: impl FnMut(usize) -> Vec<u32>,
    mut described: impl FnMut(usize) -> Result<Vec<u32>>,
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    load_of: impl Fn(u32) -> crate::report::db_subject::HitLoad,
    protein: bool,
) -> Result<Option<(u32, String)>> {
    let refusals: std::collections::BTreeMap<u32, String> = parsed
        .iter()
        .enumerate()
        .filter_map(|(rank, subject)| {
            let rank = u32::try_from(rank).ok()?;
            description_line_refusal(subject.as_deref()?, load_of(rank), protein)
                .map(|reason| (rank, reason))
        })
        .collect();
    if refusals.is_empty() {
        return Ok(None);
    }
    for query in 0..num_queries {
        if !reported(query)
            .iter()
            .any(|rank| refusals.contains_key(rank))
        {
            continue;
        }
        for rank in described(query)? {
            if let Some(reason) = refusals.get(&rank) {
                return Ok(Some((rank, reason.clone())));
            }
        }
    }
    Ok(None)
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
    /// Whether a tabular format can stop at a row whose subject lookup throws
    /// (`report::tabular::stops_at_lookup_failure_for`: `sacc`, `saccver`, the default fields, a
    /// first row's `staxid`), so that such a sequence may be refused (`reports_may_refuse`).
    pub(crate) stops_at_failed_lookup: bool,
    /// Whether a tabular format stops at a subject that has a definition line without a Seq-id
    /// (`report::tabular::stops_at_lineless_accession`: `sallacc`), so that such a sequence may be
    /// refused (`reports_may_refuse`).
    pub(crate) stops_at_lineless: bool,
    /// Whether the batches' prefetch reads the subjects by their Seq-align ids (the first
    /// format's: outfmt 0, 5, 6, 7, 10, 17 and 20, with `PRE_FETCH_SEQS_LIMIT=2147483647` also the
    /// Seq-annots 8, 9 and 12, [`prefetches_subjects`]), unless `PRE_FETCH_SEQS_LIMIT` turns it
    /// off for every batch ([`prefetch_always_off`]). Where no report reads their Bioseqs (the
    /// Seq-annots, a table of fields that read nothing of the subject), what the prefetch
    /// cannot read as NCBI does is refused (`check_prefetched_parsed`) for a subject that a
    /// batch whose prefetch is on reports ([`BatchReads::shown`], [`BatchPrefetch`]:
    /// a batch whose Seq-aligns outnumber a positive limit reads nothing).
    pub(crate) prefetched: bool,
    /// What the reports read of a definition line besides the Seq-aligns ([`TitleReads`],
    /// [`check_reported_title`]; by default the batches' prefetch or a formatter).
    pub(crate) titles: TitleReads,
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1773-1783
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3289-3300
// The prefetch (`LoadSequencesToScope`'s `GetBioseqHandles`) asks the scope only for the ids it
// has resolved, which it gets back as they were resolved, and the data loader for the others by
// their own ids; the data source's weak match (data_source.cpp:1439-1471) looks only among the
// blobs that the loader gives for that id. So the weak matches of
// `CDataSource_ScopeInfo::x_GetSeqMatch` (scope_info.cpp:865-882) come from the formatters' (and
// psiblast's PSSM's) own reads where the batch's prefetch is off. A later batch can get them back
// only where the scope is kept: psiblast's rounds of one query, and outfmt 5. Probe w0 of line
// DBID weakbatch (`blastdb-ids-20261021/R/weakbatch/w0/run_w0.sh`): a batch of 2 Seq-aligns under
// limit 2 shows each subject, the next batch of 4 shows `sp|P00001.1|NAME_B` for `sp||NAME_B`;
// outfmt 5 with the batches before kept, and a batch of 252 subjects, which the prefetch resolves
// 100 at a time, show each subject. `report::db_subject::ScopeWalk::bulk` resolves the prefetch's
// ids without the scope's index, and its `load` the formatters' reads with it.
impl ParsedShown {
    /// What the reports show of a subject that the batches' prefetch reads only where `read` (a
    /// batch that reports it prefetches; [`BatchReads::shown`]).
    pub(crate) fn of_batches(self, read: bool) -> Self {
        Self {
            prefetched: self.prefetched && read,
            ..self
        }
    }

    /// Whether a report reads the subject's Bioseq, which NCBI's formatters find by its
    /// Seq-align id (`ParsedSubject::fetched`); the Seq-annots, the archive and the tables of
    /// fields that read nothing of the subject write the Seq-align id or nothing of it.
    pub(crate) fn reads_bioseq(self) -> bool {
        self.tabular || self.pairwise || self.structured || self.sam || self.taxonomy
    }

    /// Whether NCBI's scope resolves the subject's Seq-align id: a report reads its Bioseq
    /// (`reads_bioseq`), or the batches' prefetch loads it (`prefetched`), whose texts (the
    /// taxonomy database's lookups, the warnings and errors of the load) are those of the record
    /// that the scope finds. Where a query of the batch holds the id, the scope finds the query
    /// (`refuse_id_of_a_query`); the Seq-annots, the archive and a table of fields that read
    /// nothing of the subject, without the prefetch, write the Seq-align id alone.
    pub(crate) fn resolves_subject(self) -> bool {
        self.reads_bioseq() || self.prefetched
    }
}

/// LOSAT's refusals of the reported database sequence `ordinal` with parsed Seq-ids `subject` in
/// the reports `shown` (`protein` for a protein database; `taxdb_found`: `SeqDb::taxdb_found`),
/// as BLASTX's (`blastx::database::ReportedSubjects::for_batch`): a title with a byte outside
/// ' ' ..= '~' that the reports read besides the Seq-aligns (`check_reported_title`), a sequence that NCBI's
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
    seen_by_scope: bool,
    program: &str,
) -> Result<()> {
    check_reported_title(
        db,
        ordinal,
        subject.replaced_title_bytes(),
        shown.titles,
        program,
    )?;
    if shown.reads_bioseq() {
        subject
            .fetched_with_scope(shown.taxonomy, load, seen_by_scope)
            .map_err(|reason| parsed_refusal(db.oid(ordinal), &reason, program))?;
    } else if shown.prefetched {
        check_prefetched_parsed(db, subject, load, ordinal, seen_by_scope, program)?;
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
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
                Some(subject) if subject.scope_finds_none() => Ok(None),
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
/// `StoredTitle::check_shown_source`, `check_reported_parsed` with `taxdb_found` and
/// `description_line_refusal` where outfmt 0 is shown, as `reported_subjects` asks). Where the reports join following hits (`shown.joined`), no two sequences with parsed
/// Seq-ids may be joined but copies of one record (`refuse_joined_hits`): the sequences that share
/// a key of their Seq-align ids (`db_subject::join_keys`; joined ones always share one) are read
/// again and compared (`joined_in_reports`, `joins_as_copy`). With `prefetch` (the batches' prefetch of a search with `taxdb.bti` on
/// the search path) the taxids that the prefetch looks up must read (a stored title's first one,
/// every one of a sequence with parsed Seq-ids, as `prefetch_messages` reads them) and no lookup
/// of them may be refused (`taxids_may_refuse`). Under `-parse_deflines` no sequence may have the
/// id of a query of `queries` (`refuse_id_of_a_query`; the scope of every query of the run, a
/// superset of each batch's). No load may be given another sequence by NCBI's places
/// (`ParsedSubject::place_keys`): the data loader's `m_Ids` always, the scope's index only where a
/// batch may read its subjects one by one (`one_by_one`: [`prefetch_may_be_off`] of the first
/// format) or after a prefetch that fails (`ParsedSubject::prefetch_fails`), since the prefetch
/// resolves its ids without that index and the reports' loads then find them in the scope's
/// cache (`report::db_subject::ScopeWalk::bulk`). A reader error is a refusal too. The headers are
/// read in `pool`; nothing is written, and the one state of the run that changes is the note of a
/// title with a byte outside ' ' ..= '~' (`SeqDb::note_replaced_title`): the model of NCBI's cache
/// of header records then starts at the run's first Seq-align read, from an empty cache, which
/// gives the same posts. `None`: no report of any sequence can be refused.
pub(crate) fn reports_may_refuse(
    db: &SeqDb,
    shown: ParsedShown,
    protein: bool,
    taxdb_found: bool,
    (prefetch, one_by_one): (bool, bool),
    queries: Option<&crate::blastinput::query_batch::QueryScope>,
    pool: &crate::utils::threading::SearchPool<'_>,
    program: &str,
) -> Option<ReportRefusal> {
    /// What the pass keeps: the join keys with their sequence, the prefetch's taxids, the keys
    /// of NCBI's places (`ParsedSubject::place_keys`: those that each loaded blob holds, by its
    /// ordinal id, and those that each subject's load looks for in `m_Ids` and in the scope's
    /// index, with its own ordinal id and the one its lookup gives), and whether a subject's
    /// prefetch fails.
    #[derive(Default)]
    struct Kept {
        keys: Vec<(u64, u32)>,
        taxids: Vec<i32>,
        held: Vec<(u64, u32)>,
        sought: Vec<(u64, u32, u32)>,
        scope_sought: Vec<(u64, u32, u32)>,
        prefetch_fails: bool,
    }
    impl Kept {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        fn join(mut self, mut other: Kept) -> std::result::Result<Kept, ReportRefusal> {
            self.keys.append(&mut other.keys);
            self.taxids.append(&mut other.taxids);
            self.held.append(&mut other.held);
            self.sought.append(&mut other.sought);
            self.scope_sought.append(&mut other.scope_sought);
            self.prefetch_fails |= other.prefetch_fails;
            Ok(self)
        }
    }
    let check = |kept: &mut Kept,
                 headers: &HeaderReader<'_>,
                 ordinal: u32,
                 header: std::result::Result<crate::seqdb::SequenceHeader, SeqDbError>|
     -> std::result::Result<(), ReportRefusal> {
        // A sequence left without a line with a Seq-id gets no Seq-align (`SeqAligned`): none of
        // its hits is reported, so it cannot make a report refuse (an error of its header comes
        // below). The header read here is the filtered one: a sequence with a line that holds a
        // Seq-id is kept without asking again (review 1 of line RG D1c-1b, F2; an ordinal line
        // holds its `BL_ORD_ID`). A set of no such line and a read that failed ask
        // `SeqDb::keeps_line` (an empty record is kept, and a seqid list's read of a sequence it
        // leaves no line fails; review 2 of D1c-1b, G3).
        if match &header {
            Ok(crate::seqdb::SequenceHeader::Ordinal(lines)) if !lines.is_empty() => false,
            Ok(crate::seqdb::SequenceHeader::Parsed(set))
                if set.deflines.iter().any(|defline| !defline.ids.is_empty()) =>
            {
                false
            }
            _ => matches!(db.keeps_line(ordinal), Ok(false)),
        } {
            return Ok(());
        }
        match reported_header(headers, header, ordinal, 0, program) {
            Ok(ReportedHeader::Stored(stored)) => {
                check_reported_title(db, ordinal, stored.replaced_bytes, shown.titles, program)
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
                    false,
                    program,
                )
                .map_err(|_| ReportRefusal::Headers)?;
                // A description line that LOSAT refuses where a report prints it
                // (`refuse_described_lines`, after the search). The pre-pass reads no `-show_gis`:
                // the GI that NCBI never sets (`ParsedSubject::description_gi_unset`) is refused in
                // the batch that reports it, after the earlier batches' reports, which are NCBI's
                // own prefix (review 2 of line RG-KB item 12, D1).
                if shown.pairwise
                    && description_line_refusal(
                        &subject,
                        crate::report::db_subject::HitLoad::Loaded,
                        protein,
                    )
                    .is_some()
                {
                    return Err(ReportRefusal::Headers);
                }
                // A sequence whose lookups fail (`ParsedSubject::lookup_failure`), where a field
                // list can stop at it (`failed_lookup_stop`: `sacc`/`saccver`/default/`staxid` at
                // a thrown lookup, `sallacc` at a definition line without a Seq-id).
                if shown.tabular
                    && (shown.stops_at_failed_lookup || shown.stops_at_lineless)
                    && match subject.lookup_failure(protein) {
                        Ok(None) => false,
                        Ok(Some(failure)) => {
                            (shown.stops_at_failed_lookup
                                && (failure.subject_id || failure.all_ids))
                                || (shown.stops_at_lineless && failure.lineless)
                        }
                        Err(_) => true,
                    }
                {
                    return Err(ReportRefusal::Headers);
                }
                if let Some(queries) = queries.filter(|_| shown.resolves_subject()) {
                    refuse_id_of_a_query(&subject, || db.oid(ordinal), queries, program)
                        .map_err(|_| ReportRefusal::Headers)?;
                    if shown.pairwise && queries.has_key(&subject.query_scope_keys().display) {
                        return Err(ReportRefusal::Headers);
                    }
                }
                if shown.joined {
                    let mut keys = Vec::new();
                    crate::report::db_subject::join_keys(&subject, &mut keys);
                    kept.keys.extend(keys.into_iter().map(|key| (key, ordinal)));
                }
                // NCBI's places may give another loaded sequence for a subject (a hold, which the
                // reports refuse; `report::db_subject::ScopeWalk`): kept for the check below.
                let (mut sought, mut scope_sought) = (Vec::new(), Vec::new());
                subject.place_keys(&mut kept.held, (&mut sought, &mut scope_sought));
                let gives = subject.lookup_gives_ordinal().unwrap_or(ordinal);
                kept.sought
                    .extend(sought.into_iter().map(|key| (key, ordinal, gives)));
                kept.scope_sought
                    .extend(scope_sought.into_iter().map(|key| (key, ordinal, gives)));
                kept.prefetch_fails |= subject.prefetch_fails();
                if prefetch {
                    match prefetch_reads(db, ordinal, &subject, true, false, program) {
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
        let headers = HeaderReader::new(db, false);
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
    // A key that a load looks for held by another blob than the subject's or its lookup's: a
    // load of some run can be given that sequence (a hold, `ParsedSubject::fetched`).
    kept.held.sort_unstable();
    kept.held.dedup();
    let scope_index = one_by_one || kept.prefetch_fails;
    let scope_sought = kept.scope_sought.iter().filter(|_| scope_index);
    let held_by_other = kept
        .sought
        .iter()
        .chain(scope_sought)
        .any(|&(key, own, gives)| {
            let start = kept.held.partition_point(|&(held, _)| held < key);
            kept.held[start..]
                .iter()
                .take_while(|&&(held, _)| held == key)
                .any(|&(_, holder)| holder != own && holder != gives)
        });
    if held_by_other {
        return Some(ReportRefusal::Headers);
    }
    drop((
        std::mem::take(&mut kept.held),
        std::mem::take(&mut kept.sought),
        std::mem::take(&mut kept.scope_sought),
    ));
    kept.keys.sort_unstable();
    let headers = HeaderReader::new(db, false);
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
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:67-70
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:306-313
// NCBI BLAST+ 2.17.0 shows `#` for a byte outside ' ' ..= '~' of a stored title and posts, for
// each such byte, each time it reads the definition line (`seqid::ber::title_replaced_chars`):
// where it makes a batch's Seq-aligns (`GetSeqIDs` of each subject of each query's hit list, through
// the volume's cache of header records, `SeqDb::read_header_record`; `seq_align_title_posts`), at
// each load of the subject's Bioseq (the prefetch, or the first read of a formatter: the header
// through that cache, and the registration of its ids, `prefetch_title_posts`), and at each read of
// the formatters ([`FormatterTitleReads`]). Probes of lines RG-KB 9(b) and 9b-2,
// `K/kb/i9b/agent1/probe` (1,368 cases, 5 programs) and `K/kb/i9b2/probe`.
/// What the reports of a run read of a reported database sequence's definition line besides
/// NCBI's Seq-aligns of its batch ([`check_reported_title`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum TitleReads {
    /// Nothing: every format writes the Seq-aligns as they are (outfmt 8, 9, 11 and 12,
    /// `ReportFormat::writes_seqaligns`), the first one does not prefetch the subjects
    /// ([`prefetches_subjects`]: a `PRE_FETCH_SEQS_LIMIT` of `INT_MAX` makes 8, 9 and 12 prefetch),
    /// and the application is one of the five search programs (`blastn`, `blastp`, `blastx`,
    /// `tblastn`, `tblastx`; not `psiblast`, `rpsblast` or `rpstblastn`).
    SeqAligns,
    /// Each batch's prefetch loads the subjects (`PRE_FETCH_SEQS_LIMIT` turns it off for no
    /// batch: [`prefetch_may_be_off`]; `prefetch_title_posts`), and the run's one format reads the
    /// definition line where LOSAT writes NCBI's posts of it ([`FormatterTitleReads::written`]).
    Prefetched,
    /// A formatter that LOSAT does not write the posts of, a prefetch that may be off, several
    /// formats, `psiblast`, `rpsblast`, or `rpstblastn` reads it too.
    #[default]
    Reports,
}

impl TitleReads {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-322
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-447
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-707
    // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:195-215
    // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:191-209
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:250-262
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:374-400
    /// What the reports of the formats `outfmts` (NCBI's numbers, the first format first) of a
    /// run of the application `app` read.
    ///
    /// `RequestKind::Psiblast`, `RequestKind::Rpsblast`, and `RequestKind::Rpstblastn` always
    /// return [`TitleReads::Reports`]:
    /// - `psiblast` reads subject definition lines when `CPsiBlastInputData` loads each included
    ///   subject through `CScope` to compute a PSSM (`-num_iterations` != 1, including 0, or
    ///   `-save_pssm_after_last_round` with `-out_pssm` or `-out_ascii_pssm`), and its traceback
    ///   makes Seq-aligns for the longer hit list before `CBlastFormatUtil::PruneSeqalign` cuts it
    ///   to `-max_target_seqs`.
    /// - `rpsblast` and `rpstblastn` construct a fresh `CLocalRPSBlast` per query batch and open a
    ///   fresh `CSearchDatabase` (`CSeqDB`, empty header cache) per volume in `s_RunLocalRpsSearch`,
    ///   making Seq-aligns in volume order before `s_CombineSearchSets` merges them and before
    ///   `PruneSeqalign` cuts them to `-max_target_seqs`.
    ///
    /// None of those extra reads or per-batch/per-volume caches is modelled.
    pub(crate) fn of_outfmts(
        app: crate::report::archive::options::RequestKind,
        outfmts: impl IntoIterator<Item = i32>,
    ) -> Self {
        use crate::blastinput::app::ReportFormat;
        if matches!(
            app,
            crate::report::archive::options::RequestKind::Psiblast
                | crate::report::archive::options::RequestKind::Rpsblast
                | crate::report::archive::options::RequestKind::Rpstblastn
        ) {
            return TitleReads::Reports;
        }
        let mut outfmts = outfmts.into_iter();
        let Some(first) = outfmts.next() else {
            return TitleReads::Reports;
        };
        let second = outfmts.next();
        // `s_IsPrefetchFormat` takes the tables 10 and 20 as 6 (cited at `PrefetchFormat`).
        let prefetch = PrefetchFormat::of(
            match ReportFormat::from_number(first) {
                Some(
                    ReportFormat::CommaSeparatedValues
                    | ReportFormat::CommaSeparatedValuesWithHeader,
                ) => ReportFormat::Tabular as i32,
                _ => first,
            },
            app,
        );
        // The Seq-annots and the archive are formats that the prefetch takes as they are.
        let seq_aligns = |number: i32| {
            ReportFormat::from_number(number).is_some_and(ReportFormat::writes_seqaligns)
        };
        if seq_aligns(first)
            && second.into_iter().chain(outfmts).all(seq_aligns)
            && !prefetches_subjects(Some(prefetch))
        {
            return TitleReads::SeqAligns;
        }
        let written = ReportFormat::from_number(first)
            .is_some_and(|format| FormatterTitleReads::of(format).written());
        if second.is_none()
            && written
            && prefetches_subjects(Some(prefetch))
            && !prefetch_may_be_off(prefetch)
        {
            TitleReads::Prefetched
        } else {
            TitleReads::Reports
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:809-842
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1231-1268
// Every read of a formatter decodes the header that the subject's Bioseq holds
// (`CSeqDB::ExtractBlastDefline`), not through the volume's cache: one that holds the stored bytes
// posts each time (`SeqDb::rewrites_header_at`). Counted by the probes of line RG-KB 9(b) (the
// 05:31 entry of `K/kb/i9b/agent1/logs/progress.log`) and 9b-2, the same in the five programs,
// with and without `-html`.
/// The reads of a reported subject's definition line that NCBI's formatter of one format makes:
/// the one table of them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FormatterTitleReads {
    /// None: SAM (17), and the formats that write the Seq-aligns as they are (8, 9, 11, 12).
    None,
    /// One at each row whose fields read the definition line (`TabularSpec::reads_defline`): the
    /// tables 6, 7, 10 and 20.
    Rows,
    /// Two at each description line, and one at each alignment heading where `headings`: outfmt
    /// 0 (`headings`), 1 to 4.
    Pairwise { headings: bool },
    /// One at each `Hit`: BLAST XML (5), and the XML2 and JSON formats (13 to 16), which read
    /// all of a query's before its report.
    Hits,
    /// Two at each Seq-align of a subject: the organism report (18).
    Hsps,
}

impl FormatterTitleReads {
    /// The reads of the formatter of `format`.
    pub(crate) fn of(format: crate::blastinput::app::ReportFormat) -> Self {
        use crate::blastinput::app::ReportFormat::*;
        match format {
            Pairwise => FormatterTitleReads::Pairwise { headings: true },
            QueryAnchoredIdentities
            | QueryAnchoredNoIdentities
            | FlatQueryAnchoredIdentities
            | FlatQueryAnchoredNoIdentities => FormatterTitleReads::Pairwise { headings: false },
            Xml | Json | Xml2 | JsonS | Xml2S => FormatterTitleReads::Hits,
            Tabular
            | TabularWithComments
            | CommaSeparatedValues
            | CommaSeparatedValuesWithHeader => FormatterTitleReads::Rows,
            TaxFormat => FormatterTitleReads::Hsps,
            AsnText | AsnBinary | ArchiveFormat | JsonSeqalign | Sam => FormatterTitleReads::None,
        }
    }

    /// Whether LOSAT writes NCBI's posts of these reads ([`TitleReads::Prefetched`]): the tables'
    /// rows, and the formats that read nothing.
    pub(crate) fn written(self) -> bool {
        matches!(self, FormatterTitleReads::None | FormatterTitleReads::Rows)
    }
}

/// LOSAT's refusal of a reported database sequence `ordinal` of `db` whose stored title held a
/// byte outside ' ' ..= '~' (`replaced_bytes`, `StoredTitle`) where the run's reports read its
/// definition line where LOSAT does not write NCBI's posts of it (`reads`, [`TitleReads`]; with
/// the prefetch, in a database with parsed Seq-ids, or where an environment variable of NCBI's
/// object manager is set, which changes what it keeps between the batches:
/// `default_object_manager`), or where its volume filters the header records it reads (a seqid or
/// taxid list, of the search or of an alias, or an OID mask type; `SeqDb::filters_header_sets`):
/// NCBI's traceback then reads the headers of its fetches too (`check_oid_exclusion`). NCBI's
/// bytes are determined in all of these; LOSAT does not model them yet (line RG-KB item 9b-3: the
/// prefetch's order of parsed Seq-ids, known for packed ids, with the address order of the
/// handles of one key; the object manager's sizes; the traceback's header reads).
pub(crate) fn check_reported_title(
    db: &SeqDb,
    ordinal: u32,
    replaced_bytes: bool,
    reads: TitleReads,
    program: &str,
) -> Result<()> {
    if replaced_bytes
        && (reads == TitleReads::Reports
            || db
                .location_of(ordinal)
                .is_none_or(|at| db.filters_header_sets(at))
            || (reads == TitleReads::Prefetched && !prefetch_title_posts_known(db)))
    {
        bail!(
            "database sequence {} has a stored title with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line); this is not supported by LOSAT's {program}",
            db.oid(ordinal)
        );
    }
    Ok(())
}

/// Whether LOSAT models what the loads of NCBI's prefetch of the subjects of `db` post
/// ([`prefetch_title_posts`]): a database without parsed Seq-ids (whose `gnl|BL_ORD_ID|n` handles
/// the prefetch orders by their numbers alone, [`ordinal_prefetch_places`]), with NCBI's object
/// manager at its defaults. The others are determined too but not modelled yet ([`check_reported_title`]).
#[cold]
fn prefetch_title_posts_known(db: &SeqDb) -> bool {
    matches!(db.has_parsed_seqids(), Ok(false)) && default_object_manager()
}

/// NCBI's posts where it makes the Seq-aligns of a batch's results (`GetSeqIDs` of each subject,
/// [`check_reported_title`]): `reads` are the ordinal ids of the subjects of the batch's final hit
/// lists that get a Seq-align, query by query in the order of each list, every one a read of its
/// header record through its volume's cache (`SeqDb::read_header_record`); one that reads the
/// record posts for each byte outside ' ' ..= '~' of its titles (`SeqDb::title_posts_at`), with
/// the post prefix of `program`. Empty, and nothing read, where no reader met such a title
/// (`SeqDb::replaced_titles`). The caller writes them after the batch's reader messages, its
/// search strategy request's posts and NCBI's warning of its first read of a UTF8String header,
/// before its prefetch's texts.
#[inline]
pub(crate) fn seq_align_title_posts(
    db: &SeqDb,
    reads: impl IntoIterator<Item = u32>,
    program: &str,
) -> crate::report::diag_log::Posts {
    if !db.replaced_titles() {
        return crate::report::diag_log::Posts::default();
    }
    seq_align_title_reads(db, &mut reads.into_iter(), program)
}

/// [`seq_align_title_posts`] of a batch whose reads posted NCBI's warning of its first read of a
/// UTF8String header or not (`utf8_posted`), followed where the batch's prefetch loads its subjects
/// (`prefetch`, the run's taxonomy database: the run's reads are [`TitleReads::Prefetched`]) by
/// the posts of those loads ([`prefetch_title_posts`]): with that warning, LOSAT's refusal, which
/// names the first of `reads` whose title posts (the order of that warning among the Seq-aligns'
/// posts is not modelled).
pub(crate) fn batch_title_posts(
    db: &SeqDb,
    reads: impl IntoIterator<Item = u32> + Clone,
    utf8_posted: bool,
    prefetch: Option<&TaxdbRun>,
    program: &str,
) -> Result<crate::report::diag_log::Posts> {
    let mut posts = seq_align_title_posts(db, reads.clone(), program);
    if let Some(taxonomy) = prefetch {
        if db.replaced_titles() {
            prefetch_title_posts(
                db,
                taxonomy,
                &mut reads.clone().into_iter(),
                &mut posts,
                program,
            )?;
        } else {
            db.note_unmodelled_prefetch();
            taxonomy.note_title_loads(db, reads.clone());
        }
    }
    if utf8_posted && !posts.is_empty() {
        let prefix = program.to_ascii_lowercase();
        let posting = reads.into_iter().find(|&ordinal| {
            db.location_of(ordinal)
                .is_some_and(|at| !db.title_posts_at(at, &prefix).is_empty())
        });
        bail!(
            "database sequence {} has a stored title with a byte outside ' ' to '~', in a query batch whose reads also post NCBI BLAST+'s warning of a UTF8String header (803.10), whose order among their error messages LOSAT does not reproduce; this is not supported by LOSAT's {program}",
            posting.map_or_else(|| "?".to_string(), |ordinal| db.oid(ordinal).to_string())
        );
    }
    Ok(posts)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:755-772
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1111-1149
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-406
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:259-271
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:58
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:74-110
// The prefetch loads the batch's Seq-align ids in the order of their Seq-id handles (by packed
// number: `gnl|BL_ORD_ID|1` to the last OID, then `gnl|BL_ORD_ID|0`; [`prefetch_places`]), each
// subject once, unless NCBI's data source still holds its blob from a batch before
// (`SeqDb::batch_blobs`). A load reads the subject's header through its volume's cache (`x_GetOid`'s
// `GetSeqIDs`, `GetBioseqNoData`: the first of them may read the record, the others find it), then
// decodes the header that the new Bioseq holds to register its ids (`RegisterIds`), which posts
// unless the Bioseq holds a header written again (`SeqDb::rewrites_header_at`).
/// The posts of the loads of a batch's prefetch of the subjects of its Seq-aligns `reads` (as
/// [`seq_align_title_posts`] takes them), after its Seq-aligns' `posts`. LOSAT's refusal where
/// what a load that posts does is not known: NCBI's data source may hold its blob or not (a release
/// of more than it keeps, in an order from the addresses of NCBI's objects; rule 3), or a subject
/// whose blob it may hold or not takes the same slot of the header cache as one whose title posts
/// (whether its load reads the record decides whether the other's finds its own). And LOSAT's
/// refusal of what NCBI determines but LOSAT does not model yet: a batch whose loads of a subject
/// with a positive taxid write texts of the run's taxonomy database (`taxonomy`,
/// [`TaxdbRun::title_loads_write_texts`]: NCBI writes them in the load that looks the taxid up,
/// after its header read and before the decode, `CSeqDBVol::GetBioseq`'s `x_GetTaxonomy`, among
/// the loads' posts, where LOSAT
/// writes the prefetch's texts after its posts), and a batch whose map LOSAT does not know
/// ([`UNNAMED_MAP`]).
#[cold]
#[inline(never)]
fn prefetch_title_posts(
    db: &SeqDb,
    taxonomy: &TaxdbRun,
    reads: &mut dyn Iterator<Item = u32>,
    posts: &mut crate::report::diag_log::Posts,
    program: &str,
) -> Result<()> {
    use crate::seqdb::BlobHeld;
    let prefix = program.to_ascii_lowercase();
    let mut seen = std::collections::HashSet::new();
    // Each subject once, with its place and what a read of its title posts.
    let subjects: Vec<(u32, Option<crate::seqdb::SeqLocation>, Vec<Vec<u8>>)> = reads
        .filter(|&ordinal| seen.insert(ordinal))
        .map(|ordinal| {
            let at = db.location_of(ordinal);
            let texts = at.map_or_else(Vec::new, |at| db.title_posts_at(at, &prefix));
            (ordinal, at, texts)
        })
        .collect();
    let held = db.batch_blobs(
        &subjects
            .iter()
            .map(|(ordinal, _, texts)| (*ordinal, !texts.is_empty()))
            .collect::<Vec<_>>(),
    );
    let places = ordinal_prefetch_places(
        db,
        &subjects
            .iter()
            .map(|&(ordinal, _, _)| ordinal)
            .collect::<Vec<_>>(),
    );
    let refuse = |ordinal: u32, why: &str| -> anyhow::Error {
        anyhow::anyhow!(
            "database sequence {} has a stored title with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line), in a query batch whose prefetch LOSAT does not model: {why}; this is not supported by LOSAT's {program}",
            db.oid(ordinal)
        )
    };
    let slot = |at: Option<crate::seqdb::SeqLocation>| at.map(|at| db.header_slot(at));
    let load_texts = std::cell::OnceCell::new();
    for (index, (ordinal, at, texts)) in subjects.iter().enumerate() {
        if texts.is_empty() {
            continue;
        }
        if held[index] == BlobHeld::Unsure {
            return Err(refuse(*ordinal, "NCBI's data source may still hold its sequence from a batch before, after a release of more sequences than it keeps (10), in an order that comes from the addresses of NCBI's objects"));
        }
        if let (BlobHeld::Out, Some(why)) = (held[index], places[index].many) {
            return Err(refuse(*ordinal, why));
        }
        if held[index] == BlobHeld::Out
            && *load_texts.get_or_init(|| {
                taxonomy.title_loads_write_texts(
                    db,
                    subjects
                        .iter()
                        .zip(&held)
                        .filter(|(_, state)| **state != BlobHeld::Held)
                        .map(|((ordinal, _, _), _)| *ordinal),
                )
            })
        {
            return Err(refuse(*ordinal, "a load of the batch's prefetch writes a text of the run's taxonomy database (NCBI writes it in the load that looks a taxid up, after that load's header read and before its decode: the texts of opening the database at the run's first such load, the offset error at a lookup of its taxid that the volume's cache does not hold), among the loads' error messages, which LOSAT does not place among those messages yet"));
        }
        // Whether such a load reads its record decides what the slot holds for this subject's
        // load, and for its reads in the batches after, held or not.
        if let Some((other, _)) =
            subjects
                .iter()
                .zip(&held)
                .find(|((other, at_other, _), state)| {
                    **state == BlobHeld::Unsure
                        && *other != *ordinal
                        && slot(*at_other) == slot(*at)
                })
        {
            return Err(refuse(*ordinal, &format!("NCBI's data source may still hold the sequence {} from a batch before, whose load would read the header record that takes the same slot of the volume's cache", db.oid(other.0))));
        }
    }
    let mut order: Vec<usize> = (0..subjects.len()).collect();
    order.sort_by_key(|&index| places[index].order);
    for index in order {
        let (_, Some(at), texts) = &subjects[index] else {
            continue;
        };
        if held[index] != BlobHeld::Out {
            continue;
        }
        let read = db.read_header_record(*at);
        if texts.is_empty() {
            continue;
        }
        // An error (`eDiag_Error`), which a BLAST archive saves as one message.
        if read {
            for text in texts {
                posts.post(crate::report::diag_log::ERROR, text.clone());
            }
        }
        if !db.rewrites_header_at(*at) {
            for text in texts {
                posts.post(crate::report::diag_log::ERROR, text.clone());
            }
        }
    }
    taxonomy.note_title_loads(db, subjects.iter().map(|(ordinal, _, _)| *ordinal));
    Ok(())
}

/// NCBI's posts of each read of the definition lines of the reported database subjects that a
/// formatter makes ([`FormatterTitleReads`]), as [`crate::report::query_warnings::TitleDecodes`]
/// by rank, the key of every program's report writers: `ordinals` are the ordinal ids in `db` of
/// the reported subjects, by rank (the ranks that the final hit lists name them by). Only those
/// whose read posts (a title with a byte outside ' ' ..= '~', in a Bioseq that holds the stored
/// header: `SeqDb::rewrites_header_at`); none where no reader met such a title
/// (`SeqDb::replaced_titles`).
pub(crate) fn title_decodes(
    db: &SeqDb,
    ordinals: &[u32],
    program: &str,
) -> crate::report::query_warnings::TitleDecodes {
    if !db.replaced_titles() {
        return Default::default();
    }
    title_decode_reads(db, ordinals, program)
}

/// [`title_decodes`] of a run whose reader met such a title.
#[cold]
#[inline(never)]
fn title_decode_reads(
    db: &SeqDb,
    ordinals: &[u32],
    program: &str,
) -> crate::report::query_warnings::TitleDecodes {
    let prefix = program.to_ascii_lowercase();
    let mut decodes = crate::report::query_warnings::TitleDecodes::default();
    for (key, &ordinal) in (0u32..).zip(ordinals) {
        let Some(at) = db.location_of(ordinal) else {
            continue;
        };
        let texts = db.title_posts_at(at, &prefix);
        if !texts.is_empty() && !db.rewrites_header_at(at) {
            decodes.insert(key, texts.concat());
        }
    }
    decodes
}

/// [`seq_align_title_posts`] of a run whose reader met such a title.
#[cold]
#[inline(never)]
fn seq_align_title_reads(
    db: &SeqDb,
    reads: &mut dyn Iterator<Item = u32>,
    program: &str,
) -> crate::report::diag_log::Posts {
    let prefix = program.to_ascii_lowercase();
    let mut posts = crate::report::diag_log::Posts::default();
    for ordinal in reads {
        let Some(at) = db.location_of(ordinal) else {
            continue;
        };
        if db.read_header_record(at) {
            for text in db.title_posts_at(at, &prefix) {
                // An error (`eDiag_Error`), which a BLAST archive saves as one message.
                posts.post(crate::report::diag_log::ERROR, text);
            }
        }
    }
    posts
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
            .read_protein_search_block(range)
            .map_err(|error| search_read_error(error, program))?;
        for i in 0..block.len() {
            let with_sentinels = block.protein_with_sentinels(i);
            let ordinal = range.first_ordinal + i as u32;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1425-1432
            // No NCBI counterpart: LOSAT's stored-byte check (`x_check_protein_stores`). Proteins
            // with non-positive NCBI 32-bit length are skipped by NCBI and returned as `[0, 0]` by
            // `read_protein_search_block`, so skip byte validation and report length 0.
            if with_sentinels.len() == 2 {
                each(ordinal, 0);
                continue;
            }
            check_protein_sequence(with_sentinels, || db.oid(ordinal), program)?;
            each(ordinal, with_sentinels.len() - 2);
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:339-344,699-707
        if block.stopped_at_negative() {
            break;
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
    /// NCBI's places after the loads of the run's batches, rounds and PSSM reads that the program
    /// recorded ([`TaxdbRun::note_loaded`], `report::db_subject::ScopeRecord`).
    loaded: crate::report::db_subject::SharedRecord,
    /// Whether a prefetch of the run can stop, so that the run keeps `loaded`
    /// ([`TaxdbRun::prefetch_can_stop`]), decided once.
    can_stop: std::sync::OnceLock<bool>,
    /// Whether every load of a sequence with a positive taxid throws
    /// ([`TaxdbRun::load_throws`]), decided once.
    load_throws: std::sync::OnceLock<bool>,
    /// The texts of the loads that the model of the reads placed where NCBI makes no load
    /// ([`TaxdbRun::load_texts`], `report::query_warnings::LoadTexts`).
    load_texts: std::sync::Arc<crate::report::query_warnings::LoadTexts>,
    /// What the prefetch loads of the run's batches before looked up, as the posts of their titles
    /// meet them ([`TaxdbRun::title_loads_write_texts`]).
    title_lookups: std::sync::Mutex<TitleLookups>,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:126
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:128
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:100
    /// NCBI's name of the running application, in lower case, which every text of the run's
    /// lookups names as its prefix (`SetDiagPostPrefix`, set once per process by the application)
    /// where it is not the program's own: the application of a command line that runs another
    /// program's search strategy ([`TaxdbRun::running`]). `None`: the program's name in lower
    /// case ([`TaxdbRun::prefix`]).
    app: Option<&'static str>,
    /// Whether calls of the run before this one loaded sequences that no record here holds
    /// ([`TaxdbRun::note_calls_before`]).
    calls_before: bool,
    /// Whether XML2/JSON taxonomy reads (`read_reported_taxonomy`) belong to a secondary format
    /// behind a non-XML2 first format (e.g. a library call requesting `["5", "16"]` or `["6", "16"]`).
    /// On a cache-colliding volume (`TaxdbRun::cache_collision`), secondary XML2/JSON is refused
    /// (`TaxdbRun::refuse_cache_collision`); non-collision runs still share the run's taxonomy
    /// handle, open state, volume lookups, and `calls_before`.
    separate_xml2: bool,
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

/// What the prefetch loads of a run's batches looked up for certain in its taxonomy database, as
/// the posts of their titles meet the batches in order (`batch_title_posts`,
/// [`TaxdbRun::note_title_loads`]): what decides whether a later batch's loads write a text
/// ([`TaxdbRun::title_loads_write_texts`]).
#[derive(Debug, Default)]
struct TitleLookups {
    /// A load looked a positive taxid up: the database is open, its texts of opening written.
    opened: bool,
    /// The volumes where a load looked the offset error's taxid up, whose slot of the volume's
    /// cache holds it from then on where no other taxid of the volume takes that slot.
    offset_volumes: std::collections::HashSet<usize>,
}

/// The size 256 of `CSeqDBVol::m_TaxCache`: the slots of a volume's taxonomy cache.
const TAX_CACHE_SLOTS: i32 = 256;

/// The slot of NCBI's taxonomy cache of a volume that `taxid` takes (`CSeqDBIntCache::Lookup`,
/// `key` ANDed with one less than the size of `m_Slots`).
fn tax_cache_slot(taxid: i32) -> i32 {
    taxid & (TAX_CACHE_SLOTS - 1)
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:77-79
// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:555-575
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:93-95
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2326-2369
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:984-987
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:434-437
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1111-1149
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
// A sequence that NCBI's BLAST database data loader loaded stays loaded while a handle holds it.
// A batch's prefetch links every sequence it loads to the Bioseq it adds to the scope, which holds
// them to the scope's reset after the batch; a report that reads with the prefetch off drops each
// handle after its read. Released, a sequence goes to the scope's queue of 10
// (`OBJMGR_SCOPE_AUTORELEASE_SIZE`; the loader lets its sequences be unloaded); one pushed out of
// that queue, and every one at the reset of the formatter's scope after each query batch (not in
// BLAST XML), goes to the data source's cache of 10 (`OBJMGR_BLOB_CACHE`, the smaller of it and the
// loader's own limit, which has none), whose oldest is dropped when it holds more. A read of a
// sequence that is in neither loads it again (its taxonomy lookups run again); one in either loads
// nothing.
/// The sequences that NCBI's scope keeps in its queue of released sequences, and the data source
/// in its cache, each (cited above): a released sequence is dropped only after this many others
/// were released after it twice over.
const RELEASED_KEPT: usize = 10;

// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:93-100
// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:66-97
/// Whether NCBI's object manager keeps the sizes of its defaults ([`RELEASED_KEPT`]): no
/// environment variable of its parameters is set (`OBJMGR_BLOB_CACHE`,
/// `OBJMGR_SCOPE_AUTORELEASE`, `OBJMGR_SCOPE_AUTORELEASE_SIZE`, `OBJMGR_SCOPE_POSTPONE_DELETE`,
/// any `OBJMGR_` one), whose values LOSAT does not read; their registry entries are refused as
/// every entry LOSAT does not know is (`blastinput::ncbi_environment`).
fn default_object_manager() -> bool {
    object_manager_default_in(std::env::vars_os())
}

/// [`default_object_manager`] of the environment `vars`.
fn object_manager_default_in(
    vars: impl IntoIterator<Item = (std::ffi::OsString, std::ffi::OsString)>,
) -> bool {
    !vars.into_iter().any(|(name, _)| {
        name.to_str()
            .is_some_and(|name| name.to_ascii_uppercase().starts_with("OBJMGR_"))
    })
}

/// Whether NCBI keeps every sequence of `db` that a run loads, so that a read again of one by the
/// same Seq-align id loads nothing: the data source's cache is dropped from only where it holds
/// more than [`RELEASED_KEPT`] sequences, so a database of no more sequences than that never loses
/// one, where the object manager keeps its sizes ([`default_object_manager`]).
fn keeps_every_load(db: &SeqDb) -> bool {
    db.num_oids() as usize <= RELEASED_KEPT && default_object_manager()
}

/// Whether a report's read of a sequence of `db` with the prefetch off may come after NCBI dropped
/// a sequence that the run read before, so that a later read loads it again where LOSAT does not
/// model it (a report's second read of a subject): each read releases what it read, and a
/// released sequence is dropped only after more than the scope's queue and the data source's cache
/// hold together ([`RELEASED_KEPT`]), so a database of no more sequences than that loses none
/// where the object manager keeps its sizes ([`default_object_manager`]).
fn reads_may_load_again(db: &SeqDb) -> bool {
    db.num_oids() as usize > 2 * RELEASED_KEPT || !default_object_manager()
}

/// The volume lookups of a run ([`TaxdbRun`]), in the order in which LOSAT records them
/// ([`TaxdbRun::load`]): only what decides the texts of the offset error, which only a taxonomy
/// database with that error records. Nothing grows with the taxids looked up.
#[derive(Default)]
struct VolumeLookups {
    /// By volume, what the slot of the offset error's taxid holds ([`OffsetSlot`]); only where
    /// the taxonomy database has that error.
    offset_slots: rustc_hash::FxHashMap<usize, OffsetSlot>,
    /// The blobs ([`BlobName`]), by the ordinal id of their sequence, whose loads looked up a
    /// taxid of that slot: a later read by the same Seq-id handle may find the blob in NCBI's
    /// scope or data source and load nothing, or load it again; one by another id is a load of
    /// its own, unless NCBI's scope still holds a Bioseq of those blobs for that id
    /// ([`VolumeLookups::blob_read`]).
    slot_blobs: rustc_hash::FxHashMap<u32, Vec<SlotBlob>>,
    /// A run of a program without NCBI's object manager ([`TaxdbRun::sequential`]): every read
    /// of a sequence is a load (`CSeqDB::GetBioseq`), a read again too.
    sequential: bool,
    /// NCBI keeps every sequence that the run loads across scope resets ([`keeps_every_load`] of
    /// the run's database, set at each use, [`TaxdbRun::volume_lookups`]): a read again loads
    /// nothing.
    kept: bool,
    /// NCBI keeps every sequence released in the current live scope until `ResetScopeHistory`
    /// (`!reads_may_load_again` of the run's database: at most `2 * RELEASED_KEPT` sequences with
    /// default `OBJMGR_*`, set at each use, [`TaxdbRun::volume_lookups`]).
    in_scope_kept: bool,
    /// Whether NCBI's formatter keeps its scope across query batches (`RunThrow::keeps_scope`:
    /// BLAST XML, outfmt 5), so [`Self::reset_scope`] leaves the scope intact.
    keeps_scope: bool,
    /// While the loads of one batch's prefetch are walked ([`Self::prefetch_slots`]): the number
    /// of blobs that each sequence had before it, the only ones that [`Self::blob_read`] meets
    /// (NCBI's prefetch asks its scope only for the ids that it resolved before, then its loader
    /// for every other handle: no blob of one prefetch stands in for another of it).
    prefetch_from: Option<rustc_hash::FxHashMap<u32, usize>>,
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
    // A prefetch resolves in its scope the ids that the scope resolved before, which holds what it
    // read until `ResetScopeHistory` (never in BLAST XML), and asks the data source for the others.
    /// Whether NCBI's live scope may hold a blob that `slot_blobs` does not name: a prefetch or a
    /// read of the run since the scope's last reset ([`Self::reset_scope`], which never clears it
    /// where the run keeps its scope), whose sequences may look no taxid of the slot up
    /// ([`TaxdbRun::note_scope_read`]). Without one, a load of a failed prefetch asks NCBI's data
    /// source for its blob whatever it loads ([`Self::release_failed`]).
    unrecorded_reads: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:209-209
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:441-462
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2224-2246
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:986-1006
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1489-1491
// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:832-882
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1785
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3288-3297
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1328-1371
// NCBI's data source answers an id by its BLAST database data loader first, which names a blob by
// the OID with the id's Seq-id handle and loads a blob that the data source does not hold. A read
// of the reports finds a Bioseq that the scope holds without the loader, by the ids of that
// Bioseq (those of the definition line that its load targeted) or a handle that matches the id; a
// prefetch only by an id that the scope resolved before, and asks the loader for the batch's
// other handles. So a read of a sequence by another Seq-align id than its loads before is a load
// of its own where no Bioseq of those loads holds that id, and in one prefetch always (probe of
// line NP item O4b2-7, `tm_hd` and `tm_hdr`, four heap layouts: `gb|XYZ00002.1|` of the second
// volume, which the loader gives as the first sequence that holds it in its second line, looks
// that sequence's slot up again after its load by `gb|AAA00001.1|`: twice; with the handles'
// order the other way round, once).
/// NCBI's name of the blob that a load of a sequence reads for a Seq-align id
/// ([`VolumeLookups::offset_slot`]): the id's Seq-id handle, and the keys by which NCBI's scope
/// matches that id with the ids of a Bioseq that it holds ([`blob_name`]).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct BlobName {
    /// The key of the id's `CSeq_id_Handle` (`seqid::mask_handle::mask_handle_key`); `None`
    /// where NCBI makes a new handle at each call (a Textseq-id with neither accession nor name)
    /// or where no pure function knows it (a giim id).
    handle: Option<Vec<u8>>,
    /// The id's scope keys (`defline_set::scope_keys`) and its weak keys
    /// (`id_scope::weak_keys`).
    finds: (Vec<u64>, Vec<u64>),
    /// The keys of the ids of the blob's Bioseq (the definition line that the load targets).
    holds: BioseqKeys,
}

/// The keys of the ids of a Bioseq by which NCBI's scope finds it for an id ([`BlobName`]): their
/// handles, their scope keys, and the keys that weak matches find (`id_scope::kind_keys`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct BioseqKeys {
    handles: Vec<Vec<u8>>,
    scope: Vec<u64>,
    kinds: Vec<u64>,
}

impl BlobName {
    /// Whether NCBI's scope can give the Bioseq whose ids have the keys `bioseq` for this id
    /// without a load: one of its handle, an id that the scope takes for it, or a weak match.
    fn found_in(&self, bioseq: &BioseqKeys) -> bool {
        self.handle
            .as_ref()
            .is_some_and(|handle| bioseq.handles.contains(handle))
            || self.finds.0.iter().any(|key| bioseq.scope.contains(key))
            || self.finds.1.iter().any(|key| bioseq.kinds.contains(key))
    }
}

/// A blob of [`VolumeLookups::slot_blobs`]: its handle, its Bioseq's keys ([`BlobName`]), and
/// its scope lifetime state.
#[derive(Clone, Debug)]
struct SlotBlob {
    handle: Vec<u8>,
    holds: BioseqKeys,
    /// Whether NCBI loaded it for certain (else its scope may have given another blob's Bioseq
    /// for the id, so a later read by the handle may be its first load).
    sure: bool,
    /// Whether the blob's TSE is still indexed in NCBI's live scope since the last
    /// `ResetScopeHistory` ([`VolumeLookups::reset_scope`]).
    in_scope: bool,
    /// Whether a succeeding prefetch in the current live scope linked the blob's TSE via
    /// `SetLinkUsedTSE` (`blast_aux.cpp:1144-1153`), pinning it until the scope is reset.
    pinned: bool,
    /// Whether NCBI's data source holds the blob in its cache of released blobs without the
    /// scope indexing it: a load of a prefetch whose lookup failed, which the scope never took
    /// and whose cache kept it ([`VolumeLookups::release_failed`]), until the end of the batch.
    cached: bool,
}

/// How a load of a blob meets the run's blobs before ([`VolumeLookups::blob_read`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BlobRead {
    /// NCBI loads it.
    Fresh,
    /// NCBI holds it from a load before ([`VolumeLookups::kept`]): no load.
    Kept,
    /// A load or none, as NCBI's scope and data source decide: the blob was read before, or
    /// NCBI's scope can give another blob's Bioseq of the sequence for its id.
    Again,
}

/// NCBI's name of the blob of the sequence `blob` of `db` that a read by the Seq-align id of the
/// subject `by` loads ([`BlobName`]): the id's handle and keys, and the keys of the ids of the
/// definition line of `blob` that the load targets (`DeflineSet::bioseq_defline`; a sequence
/// without parsed Seq-ids is `gnl|BL_ORD_ID|<OID>`, its Bioseq's one id).
fn blob_name(db: &SeqDb, blob: u32, by: u32, program: &str) -> Result<BlobName> {
    use crate::seqid::ParsedSeqId;
    let ordinal_id = |ordinal: u32| {
        ParsedSeqId::General(crate::seqid::model::Dbtag {
            db: b"BL_ORD_ID".to_vec(),
            tag: crate::seqid::model::ObjectId::Id(db.oid(ordinal).get() as i32),
        })
    };
    let header = |ordinal: u32| {
        db.read_header(ordinal)
            .map_err(|error| database_error(error, program))
    };
    let id = match header(by)? {
        crate::seqdb::SequenceHeader::Parsed(set) => set.alignment_id().cloned(),
        crate::seqdb::SequenceHeader::Ordinal(_) => None,
    }
    .unwrap_or_else(|| ordinal_id(by));
    let line: Vec<ParsedSeqId> = match header(blob)? {
        crate::seqdb::SequenceHeader::Parsed(set) => set
            .bioseq_defline(&id)
            .map(|line| set.deflines[line].ids.clone())
            .unwrap_or_default(),
        crate::seqdb::SequenceHeader::Ordinal(_) => vec![ordinal_id(blob)],
    };
    let mut holds = BioseqKeys::default();
    for held in &line {
        holds.handles.extend(
            crate::seqid::mask_handle::mask_handle_key(held)
                .ok()
                .flatten(),
        );
        crate::seqid::defline_set::scope_keys(held, &mut holds.scope);
        crate::seqid::id_scope::kind_keys(held, &mut holds.kinds);
    }
    let mut finds = (Vec::new(), Vec::new());
    crate::seqid::defline_set::scope_keys(&id, &mut finds.0);
    crate::seqid::id_scope::weak_keys(&id, &mut finds.1);
    Ok(BlobName {
        handle: crate::seqid::mask_handle::mask_handle_key(&id)
            .ok()
            .flatten(),
        finds,
        holds,
    })
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
    /// Whether it was a load of a prefetch whose lookup of a Seq-align id failed, which NCBI may
    /// make before that failure or not at all ([`SlotLoad::maybe`]).
    maybe: bool,
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
/// `unnamed`: the cause is the load's Seq-align id, whose Seq-id handle LOSAT does not know
/// ([`BlobName::handle`]): NCBI makes a new one at each read, each a load of its own.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SlotUncertain {
    ordinal: u32,
    again: bool,
    after: Option<SlotFetch>,
    order: Option<SlotOrder>,
    unnamed: bool,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1144-1153
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:555-575
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2326-2371
    /// Whether NCBI is guaranteed to keep `blob` loaded without asking the data loader again:
    /// it was loaded for certain (`sure`) and is either pinned in the current scope by a
    /// succeeding prefetch (`pinned`), kept in the data source's cache across resets (`kept`),
    /// or kept in the current live scope's queue + data source cache (`in_scope && in_scope_kept`)
    /// or in that cache alone from a failed prefetch (`cached`, [`Self::release_failed`]).
    fn blob_kept(&self, blob: &SlotBlob) -> bool {
        blob.sure
            && (blob.pinned || self.kept || ((blob.in_scope || blob.cached) && self.in_scope_kept))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2638-2653
    /// Mirrors `CBlastFormat::ResetScopeHistory` at the end of a query batch: clears `in_scope`
    /// and `pinned` on every recorded blob unless `keeps_scope` is set (BLAST XML, outfmt 5), and
    /// takes a blob that a failed prefetch left in the data source's cache (`cached`) for one
    /// that NCBI may have dropped since (a read again is a load or none unless NCBI keeps every
    /// load, `kept`).
    fn reset_scope(&mut self) {
        for blobs in self.slot_blobs.values_mut() {
            for blob in blobs {
                // What a failed prefetch left in the data source's cache is followed to the end of
                // its batch only ([`Self::release_failed`]).
                blob.cached = false;
                if !self.keeps_scope {
                    blob.in_scope = false;
                    blob.pinned = false;
                }
            }
        }
        // The scope that the reset empties holds no read of the run.
        if !self.keeps_scope {
            self.unrecorded_reads = false;
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:453-462
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1785
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2638-2653
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:540-551
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:885-953
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1489-1517
    // A sequence that no fetch of the run has loaded is loaded (`x_LoadData`, whose
    // `CCachedSequence` asks `GetBioseqNoData` and so `x_GetTaxonomy`). One loaded before may
    // still be in NCBI's scope or data source, which then gives it without a load, or have been
    // released: review 4 of line DBID pssmread found a second query's PSSM loading a sequence again
    // after 40 other included subjects and not after none (`cc_big40`, `cc_ac`), which LOSAT does
    // not model. After `ResetScopeHistory` (`m_Scope->ResetDataAndHistory`), only blobs looked up
    // by the same `(oid, handle)` can hit `CDataSource::m_Blob_Map`; another handle cannot match a
    // cleared scope (`in_scope == false`) and is loaded fresh.
    /// How a load of the sequence `ordinal` by the blob name `name` (with a handle) meets the
    /// blobs of the run before whose loads looked up a taxid of the slot (`slot_blobs`): the same
    /// handle's blob is held where [`Self::blob_kept`] holds (`Kept`), else a load or none
    /// (`Again`); another handle's is `Kept` during a report read (`!in_prefetch`) when at least
    /// one `in_scope` blob holds that exact handle in its Bioseq (`holds.handles`) and every such
    /// `in_scope` blob satisfies [`Self::blob_kept`] (`x_FindBestTSEInIndex` picks the earliest
    /// loaded matching TSE), `Again` where an `in_scope` blob's Bioseq can match it
    /// ([`BlobName::found_in`]), and `Fresh` otherwise (including whenever all earlier blobs were
    /// cleared by `ResetScopeHistory`). In a prefetch, only the blobs before it count
    /// ([`Self::prefetch_from`]). A sequential run loads at every read.
    fn blob_read(&self, ordinal: u32, name: &BlobName) -> BlobRead {
        if self.sequential {
            return BlobRead::Fresh;
        }
        let (Some(blobs), Some(handle)) = (self.slot_blobs.get(&ordinal), &name.handle) else {
            return BlobRead::Fresh;
        };
        let in_prefetch = self.prefetch_from.is_some();
        let blobs = match &self.prefetch_from {
            Some(from) => &blobs[..from.get(&ordinal).copied().unwrap_or(0)],
            None => &blobs[..],
        };
        if let Some(same) = blobs.iter().find(|blob| &blob.handle == handle) {
            return if self.blob_kept(same) {
                BlobRead::Kept
            } else {
                BlobRead::Again
            };
        }
        if !in_prefetch {
            let mut exact = blobs
                .iter()
                .filter(|blob| blob.in_scope && blob.holds.handles.contains(handle))
                .peekable();
            if exact.peek().is_some() && exact.all(|blob| self.blob_kept(blob)) {
                return BlobRead::Kept;
            }
        }
        if blobs
            .iter()
            .any(|blob| blob.in_scope && name.found_in(&blob.holds))
        {
            BlobRead::Again
        } else {
            BlobRead::Fresh
        }
    }

    /// Marks the existing blob of `ordinal` matching `name.handle` as attached to the current live
    /// scope (and pinned when `pins` is true) after a [`BlobRead::Kept`] hit. Does nothing when
    /// `Kept` came from another blob's `holds.handles` (which is already in scope).
    fn touch_same_handle(&mut self, ordinal: u32, name: &BlobName, pins: bool) {
        let Some(handle) = &name.handle else {
            return;
        };
        if let Some(same) = self
            .slot_blobs
            .get_mut(&ordinal)
            .and_then(|blobs| blobs.iter_mut().find(|blob| &blob.handle == handle))
        {
            same.in_scope = true;
            same.cached = false;
            if pins {
                same.pinned = true;
            }
        }
    }

    /// Records the blob `name` of the sequence `ordinal` as loaded (`sure`) or maybe loaded
    /// ([`SlotBlob`]) in the current live scope (and pinned when `pins` is true), or updates its
    /// scope flags when its handle already has a blob of the sequence.
    fn record_blob(&mut self, ordinal: u32, name: &BlobName, sure: bool, pins: bool) {
        let Some(handle) = &name.handle else {
            return;
        };
        let blobs = self.slot_blobs.entry(ordinal).or_default();
        if let Some(same) = blobs.iter_mut().find(|blob| &blob.handle == handle) {
            same.in_scope = true;
            same.cached = false;
            if pins {
                same.pinned = true;
            }
        } else {
            blobs.push(SlotBlob {
                handle: handle.clone(),
                holds: name.holds.clone(),
                sure,
                in_scope: true,
                pinned: pins,
                cached: false,
            });
        }
    }

    /// Which of the lookups `lookups` of a load of the sequence `ordinal` in the volume `volume`
    /// by the blob name `name` ([`BlobName`]: NCBI's blob of the Seq-align id that reads it) ask
    /// the taxonomy database again for the taxid of the offset error `last` (a lookup that
    /// misses the slot after an earlier lookup of `last` in the volume), as NCBI's cache of the
    /// volume holds the slot of `last` after the lookups recorded before. A first load of a blob
    /// follows the slot ([`Self::blob_read`]: also another Seq-align id's read of a sequence
    /// loaded before, where no Bioseq of its blobs holds that id). A read again of a blob whose
    /// load looks up a taxid of the slot, or of one whose Bioseq NCBI's scope can give instead, is
    /// a load or none, as NCBI's scope decides: where loading it gives a text, LOSAT's refusal;
    /// else no text, and the slot holds what either leaves (`Either` where they differ). Where
    /// NCBI keeps every sequence of the run (`kept`), a read again of a blob it loaded for certain
    /// loads nothing: no lookup, the slot as it was. `Err` names the load that is not certain
    /// ([`SlotUncertain`]), also one by an id whose handle LOSAT does not know (`unnamed`: NCBI
    /// loads a blob at each read). In a sequential run every read is a first load: the slot is
    /// followed, never `Err`.
    fn offset_slot(
        &mut self,
        volume: usize,
        ordinal: u32,
        name: &BlobName,
        lookups: &[i32],
        last: i32,
    ) -> std::result::Result<Vec<bool>, SlotUncertain> {
        self.offset_slot_with_pins(volume, ordinal, name, lookups, last, false)
    }

    /// [`Self::offset_slot`], also pinning the blob in the current live scope when `pins` is true
    /// (for a prefetch load that `CSeqMap::CanResolveRange` links via `AddUsedTSE`).
    fn offset_slot_with_pins(
        &mut self,
        volume: usize,
        ordinal: u32,
        name: &BlobName,
        lookups: &[i32],
        last: i32,
        pins: bool,
    ) -> std::result::Result<Vec<bool>, SlotUncertain> {
        // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1785
        // The read resolves its id in the scope, where a later prefetch finds it whatever its taxids
        // (`unrecorded_reads`).
        self.unrecorded_reads = true;
        let in_slot = |taxid: i32| tax_cache_slot(taxid) == tax_cache_slot(last);
        if !lookups.iter().any(|&taxid| in_slot(taxid)) {
            return Ok(vec![false; lookups.len()]);
        }
        let before = self.offset_slots.get(&volume).copied().unwrap_or_default();
        if !self.sequential && name.handle.is_none() {
            return Err(SlotUncertain {
                ordinal,
                again: false,
                after: before.other,
                order: None,
                unnamed: true,
            });
        }
        let read = self.blob_read(ordinal, name);
        if read == BlobRead::Kept {
            self.touch_same_handle(ordinal, name, pins);
            return Ok(vec![false; lookups.len()]);
        }
        let again = read == BlobRead::Again;
        self.record_blob(ordinal, name, !again, pins);
        let slot = self.offset_slots.entry(volume).or_default();
        let uncertain = |own: Option<i32>| SlotUncertain {
            ordinal,
            again,
            after: own
                .map(|taxid| SlotFetch {
                    ordinal,
                    taxid,
                    again,
                    maybe: false,
                })
                .or(before.other),
            order: before.unordered,
            unnamed: false,
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
                    maybe: false,
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
                    maybe: false,
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
    /// ([`VolumeLookups::offset_slot_with_pins`]); loads of one key are one step
    /// ([`Self::unordered_slot`]). A blob of the prefetch is loaded whatever the prefetch loaded
    /// before it ([`Self::prefetch_from`]). Each load's [`SlotLoad::pins`] records whether
    /// `CSeqMap::CanResolveRange` pins its TSE in the live scope. `failed`: where a lookup of the
    /// batch's Seq-align ids failed, every load that LOSAT's model makes ([`FailedLoad`]), whose
    /// blobs the scope never takes ([`Self::release_failed`]). `Err`: the volume and the load
    /// whose texts NCBI's order or scope decides.
    fn prefetch_slots(
        &mut self,
        loads: &[SlotLoad],
        unordered: impl Fn(usize) -> bool,
        last: i32,
        failed: Option<&[FailedLoad]>,
    ) -> std::result::Result<rustc_hash::FxHashMap<u32, Vec<bool>>, (usize, SlotUncertain)> {
        self.prefetch_from = Some(
            self.slot_blobs
                .iter()
                .map(|(&ordinal, blobs)| (ordinal, blobs.len()))
                .collect(),
        );
        // What the scope and the data source held before the walk decides which loads go
        // through the data source (the walk records the prefetch's own blobs).
        let releases: Vec<FailedRelease> = match failed {
            Some(made) => made
                .iter()
                .map(|made| self.failed_release(loads, made))
                .collect(),
            None => Vec::new(),
        };
        let misses = self.prefetch_walks(loads, unordered, last);
        self.prefetch_from = None;
        if let (Ok(_), Some(made)) = (&misses, failed) {
            self.release_failed(loads, made, &releases);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1110-1153
        // A prefetch that does not fail leaves every load in the scope (`unrecorded_reads`).
        self.unrecorded_reads = true;
        misses
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2000-2026
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2374-2396
    // The prefetch asks its scope for the ids that it resolved before, then the data source for
    // every other handle of the piece, which locks the blob that its loader loads or finds (taking
    // it out of the cache of released blobs, `x_SetLock`), each handle a blob of its own.
    /// How the load `made` of a prefetch whose lookup of a Seq-align id failed meets NCBI's data
    /// source ([`FailedRelease`]), by the blobs before the prefetch ([`Self::prefetch_from`]). A
    /// load of the slot (one of `loads`): a blob that the scope indexes for its handle is resolved
    /// there; one that only the data source holds, or a first load, goes through the data source
    /// where NCBI surely makes the load, and a first load that NCBI may make or not through it if
    /// made; any other read may. Another load: through the data source where NCBI surely makes it
    /// ([`FailedLoad::surely`]) and the scope has read nothing that the records do not name since
    /// its reset ([`Self::unrecorded_reads`]).
    fn failed_release(&self, loads: &[SlotLoad], made: &FailedLoad) -> FailedRelease {
        let Some(load) = loads.iter().find(|load| load.ordinal == made.ordinal) else {
            return if made.surely && !self.unrecorded_reads {
                FailedRelease::Released
            } else {
                FailedRelease::Maybe
            };
        };
        let Some(handle) = &load.name.handle else {
            return FailedRelease::Maybe;
        };
        let held = match self.blob_read(load.blob, &load.name) {
            BlobRead::Fresh => false,
            BlobRead::Kept => {
                let from = self
                    .prefetch_from
                    .as_ref()
                    .and_then(|from| from.get(&load.blob).copied())
                    .unwrap_or(0);
                let in_scope = self.slot_blobs.get(&load.blob).is_some_and(|blobs| {
                    blobs[..from]
                        .iter()
                        .any(|blob| &blob.handle == handle && (blob.in_scope || blob.pinned))
                });
                if in_scope {
                    return FailedRelease::Scope;
                }
                true
            }
            BlobRead::Again => return FailedRelease::Maybe,
        };
        match (load.maybe, held) {
            (false, _) => FailedRelease::Released,
            (true, false) => FailedRelease::Unheld,
            (true, true) => FailedRelease::Maybe,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2000-2026
    // NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:974-999
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2326-2369
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2374-2396
    // `CDataSource::GetBlobs` holds the locks of the blobs that the prefetch's loads lock in its
    // local map of the batch's Seq-id handles until it returns them to the scope
    // (`CDataSource_ScopeInfo::GetBlobs`); a lookup that throws ends `CDataLoader::GetBlobs`
    // before that, so the map's destruction releases every blob that the loads before it locked,
    // whatever their taxids (one that the data loader loaded, or found in the data source's cache,
    // which takes it out of the cache: `x_SetLock`), into that cache of released blobs
    // (`x_ReleaseLastTSELock`), which takes each at its end and drops its oldest where it holds
    // more than `OBJMGR_BLOB_CACHE` (10, cited above `RELEASED_KEPT`): a blob is dropped where 10
    // or more are released after it, whatever the cache held before. The scope never indexes them,
    // and no `SetLinkUsedTSE` pins them. The map releases its entries from the largest key down:
    // probe `O4b2/impl10/r2` (NCBI 2.17.0, four heap layouts each, outfmt 0 and 5): of 11 loads of
    // taxid 99990 (`gb|FFF00011.1|` packed 11 to `gb|FFF00001.1|` packed 1) before the failed
    // lookup of `gnl|BL_ORD_ID|x050` (packed 50), the reports' read of `gb|FFF00011.1|` after a
    // load of taxid 150 of the same slot writes a second offset error (`pslotm_t5`: NCBI loads it
    // again), that of `gb|FFF00001.1|` none (`pslotr_t5`: the cache kept it); with `gb|FFF00001.1|`
    // of taxid 9606 (another slot), the same (`pslotmix_t5`), and of the failure's key, made or not
    // before 10 such loads, `gb|FFF00050.1|` is loaded again (`pslottie_t5`). Within the batch
    // nothing more is dropped where the scope keeps what it releases (`in_scope_kept`: the scope's
    // queue and the cache hold every blob).
    /// What NCBI's data source holds, after a prefetch whose lookup of a Seq-align id failed, of the
    /// blobs of its loads `loads` that look the slot up, among the loads that LOSAT's model makes,
    /// `made` ([`FailedLoad`], each met as `releases` says before the walk, [`FailedRelease`]): a
    /// blob that NCBI holds nowhere after [`RELEASED_KEPT`] others or more are released after it for
    /// certain (one released, or not made of a blob that NCBI did not hold) loses its record (a read
    /// again loads it); one released after fewer whatever NCBI made stays in the cache to the end of
    /// the batch ([`SlotBlob::cached`]); any other, and every load that NCBI may not make or whose
    /// id its scope may resolve, is a load or none at a read again (not `sure`). None of them is in
    /// the scope. A blob that a failed prefetch before left in the cache may be dropped by this one
    /// (no longer `cached`). A run with sizes of the object manager that LOSAT does not read
    /// ([`default_object_manager`]) keeps none of them for certain.
    fn release_failed(
        &mut self,
        loads: &[SlotLoad],
        made: &[FailedLoad],
        releases: &[FailedRelease],
    ) {
        if self.sequential || made.is_empty() {
            return;
        }
        for blobs in self.slot_blobs.values_mut() {
            for blob in blobs {
                blob.cached = false;
            }
        }
        let sizes = default_object_manager();
        for load in loads {
            let Some(at) = made.iter().position(|made| made.ordinal == load.ordinal) else {
                continue;
            };
            let release = releases[at];
            if release == FailedRelease::Scope {
                continue;
            }
            let (surely, most) = released_after(made, releases, at);
            let Some(handle) = &load.name.handle else {
                continue;
            };
            let Some(blobs) = self.slot_blobs.get_mut(&load.blob) else {
                continue;
            };
            let Some(index) = blobs.iter().position(|blob| &blob.handle == handle) else {
                continue;
            };
            // Released, or not made of a blob that NCBI held nowhere: dropped either way.
            let unheld = matches!(release, FailedRelease::Released | FailedRelease::Unheld);
            if sizes && unheld && surely >= RELEASED_KEPT {
                blobs.remove(index);
                continue;
            }
            let blob = &mut blobs[index];
            blob.in_scope = false;
            blob.pinned = false;
            blob.cached = sizes && release == FailedRelease::Released && most < RELEASED_KEPT;
            if !blob.cached {
                blob.sure = false;
            }
        }
    }

    /// [`Self::prefetch_slots`] once [`Self::prefetch_from`] is set.
    fn prefetch_walks(
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
                match group {
                    [load] if !load.maybe => {
                        let own = self
                            .offset_slot_with_pins(
                                volume,
                                load.blob,
                                &load.name,
                                &load.lookups,
                                last,
                                load.pins,
                            )
                            .map_err(|cause| (volume, cause))?;
                        misses.insert(load.ordinal, own);
                    }
                    // Loads of one key, or a load that NCBI makes or not: one step whose cases
                    // LOSAT counts.
                    _ => {
                        for (ordinal, own) in self
                            .unordered_slot(volume, group, last)
                            .map_err(|cause| (volume, cause))?
                        {
                            misses.insert(ordinal, own);
                        }
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
    /// order ([`unordered_boundary_misses`]). A load of a blob that NCBI holds from before loads
    /// nothing ([`BlobRead::Kept`]), as does a second load of one handle in the group (NCBI's map
    /// holds a handle once); another handle of a sequence of the group is a load of its own (the
    /// prefetch asks the loader for every handle that its scope had not resolved before it); one
    /// that NCBI makes or not ([`BlobRead::Again`]: also a load of a prefetch whose lookup of a
    /// Seq-align id failed that NCBI may make before the failure or not at all, [`SlotLoad::maybe`],
    /// review of line NP item O4b2-10, F2) is counted both ways (review 1 of line NP item O4b2-3,
    /// F2), so a group may be one such load. Where every
    /// order and every such choice gives the same number of misses (from each state the slot may
    /// hold before), the texts are the same in every case: each load gets the misses of the
    /// batch's order with every such load made (the prefetch writes the texts of its loads as one,
    /// in any order of the loads: the texts of opening first, then offset errors of one text), and
    /// the slot holds what the last load leaves, or `Either` where the cases leave other taxids.
    /// `Err`: a Seq-align id whose handle LOSAT does not know (`SlotUncertain::unnamed`), more
    /// than [`UNORDERED_CHOICES`] loads made or not, or an order or a choice that changes the
    /// number of misses.
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
                    .unwrap_or(first),
                last,
            ));
        let order = SlotOrder {
            first: first.blob,
            second: second.blob,
            taxid,
        };
        let refusal = |load: &SlotLoad, again: bool, order: Option<SlotOrder>| SlotUncertain {
            ordinal: load.blob,
            again,
            after: before.other,
            order,
            unnamed: false,
        };
        if let Some(load) = group
            .iter()
            .find(|load| !self.sequential && load.name.handle.is_none())
        {
            return Err(SlotUncertain {
                unnamed: true,
                ..refusal(load, false, None)
            });
        }
        let reads: Vec<BlobRead> = group
            .iter()
            .enumerate()
            .map(|(at, load)| {
                let of_sequence = |other: &&&SlotLoad| other.blob == load.blob;
                if self.sequential {
                    BlobRead::Fresh
                } else if group[..at]
                    .iter()
                    .filter(of_sequence)
                    .any(|other| other.name.handle == load.name.handle)
                {
                    BlobRead::Kept
                } else {
                    match self.blob_read(load.blob, &load.name) {
                        BlobRead::Fresh if load.maybe => BlobRead::Again,
                        read => read,
                    }
                }
            })
            .collect();
        if reads.iter().all(|&read| read == BlobRead::Kept) {
            for load in group {
                self.touch_same_handle(load.blob, &load.name, load.pins);
            }
            return Ok(group
                .iter()
                .map(|load| (load.ordinal, vec![false; load.lookups.len()]))
                .collect());
        }
        let maybe: Vec<usize> = (0..group.len())
            .filter(|&at| reads[at] == BlobRead::Again)
            .collect();
        if maybe.len() > UNORDERED_CHOICES {
            // More cases than LOSAT counts: the first of those loads, which NCBI makes or not.
            return Err(refusal(group[maybe[0]], true, None));
        }
        // Each load's lookups of the slot (`true` for `L`).
        let marks: Vec<Vec<bool>> = group
            .iter()
            .map(|load| {
                load.lookups
                    .iter()
                    .filter(|&&taxid| in_slot(taxid))
                    .map(|&taxid| taxid == last)
                    .collect()
            })
            .collect();
        let starts: &[bool] = match before.holds {
            SlotHolds::Other => &[false],
            SlotHolds::Last => &[true],
            SlotHolds::Either => &[false, true],
        };
        let mut total = None;
        let mut leaves = None;
        let mut own_order = false;
        for choice in 0..1usize << maybe.len() {
            let made = |at: usize| match reads[at] {
                BlobRead::Fresh => true,
                BlobRead::Kept => false,
                BlobRead::Again => {
                    let bit = maybe
                        .iter()
                        .position(|&index| index == at)
                        .expect("a choice");
                    choice >> bit & 1 == 1
                }
            };
            // The loads made, each by its first and last lookup of the slot, a run of one taxid as
            // one, and the misses inside each (an `L` after an `O` of the same load).
            let mut counts = UnorderedLoads::default();
            let mut inside = 0;
            let mut ends = Vec::with_capacity(group.len());
            for marks in (0..group.len()).filter(|&at| made(at)).map(|at| &marks[at]) {
                let (starts_last, ends_last) = (marks[0], marks[marks.len() - 1]);
                counts.add(starts_last, ends_last);
                inside += marks.windows(2).filter(|pair| !pair[0] && pair[1]).count();
                ends.push(ends_last);
            }
            let ranges: Vec<(usize, usize)> = starts
                .iter()
                .map(|&held_last| unordered_boundary_misses(held_last, counts))
                .collect();
            if ranges.iter().any(|&(fewest, most)| fewest != most) {
                // The group's own order decides.
                return Err(refusal(first, false, Some(order)));
            }
            if ranges.iter().any(|&range| range != ranges[0]) {
                // Any order gives the same misses from each state, which an earlier load left
                // open: that load is the cause (`offset_slot`'s refusal of a load after it).
                return Err(refusal(first, false, before.unordered));
            }
            let misses = inside + ranges[0].0;
            if total.is_some_and(|total| total != misses) {
                // Whether NCBI makes a load decides.
                return Err(refusal(group[maybe[0]], true, None));
            }
            total = Some(misses);
            let end = match ends.first() {
                None => before.holds,
                Some(&end) if ends.iter().all(|&other| other == end) => {
                    if end {
                        SlotHolds::Last
                    } else {
                        SlotHolds::Other
                    }
                }
                Some(_) => {
                    own_order = true;
                    SlotHolds::Either
                }
            };
            leaves = Some(match leaves {
                Some(leaves) if leaves != end => SlotHolds::Either,
                _ => end,
            });
        }
        // The batch's order with every load made or not made gives the number of misses of every
        // case (from `Other` where the slot may hold either: the same number from each).
        let mut holds = match before.holds {
            SlotHolds::Last => SlotHolds::Last,
            _ => SlotHolds::Other,
        };
        let mut other = before.other;
        let mut misses = Vec::with_capacity(group.len());
        for (load, &read) in group.iter().zip(&reads) {
            if read == BlobRead::Kept {
                self.touch_same_handle(load.blob, &load.name, load.pins);
                misses.push((load.ordinal, vec![false; load.lookups.len()]));
                continue;
            }
            let walk = holds
                .walk(&load.lookups, last)
                .expect("a walk from a definite state");
            holds = walk.holds;
            if let Some(taxid) = walk.other {
                other = Some(SlotFetch {
                    ordinal: load.blob,
                    taxid,
                    again: read == BlobRead::Again && !load.maybe,
                    maybe: load.maybe,
                });
            }
            self.record_blob(load.blob, &load.name, read == BlobRead::Fresh, load.pins);
            misses.push((load.ordinal, walk.misses));
        }
        let slot = self.offset_slots.entry(volume).or_default();
        slot.other = other;
        slot.holds = leaves.expect("one case at least");
        slot.unordered = (slot.holds == SlotHolds::Either && own_order).then_some(order);
        Ok(misses)
    }
}

/// The most loads of one step of [`VolumeLookups::unordered_slot`] that NCBI makes or not
/// ([`BlobRead::Again`]), each of whose cases LOSAT counts: more are LOSAT's refusal.
const UNORDERED_CHOICES: usize = 8;

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
    /// NCBI's name of the blob ([`blob_name`] of the subject's Seq-align id).
    name: BlobName,
    /// Whether `CSeqMap::CanResolveRange` links this load's TSE via `AddUsedTSE` (`SetLinkUsedTSE`,
    /// pinning it in the live scope until `ResetScopeHistory`): true when the prefetch reaches
    /// `AddBioseq` without throwing and this subject precedes every empty `CBioseq_Handle` of the
    /// batch in `CSeq_id::CompareOrdered` order.
    pins: bool,
    /// Whether NCBI may not make the load: a prefetch whose lookup of a Seq-align id failed, and a
    /// load whose place to that failure LOSAT does not know ([`made_before_failure`]), which
    /// LOSAT takes as made or not ([`VolumeLookups::unordered_slot`]).
    maybe: bool,
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2000-2026
// The loads of the handles before the failing one lock their blobs in the data source's map of
// the batch's handles, which the failure releases from the largest key down.
/// A load that LOSAT's model makes in a prefetch whose lookup of a Seq-align id failed
/// ([`made_before_failure`]), for what NCBI's data source holds after the failure
/// ([`VolumeLookups::release_failed`]), whether or not it looks the slot up.
#[derive(Clone, Copy, Debug)]
struct FailedLoad {
    /// The reported subject, by ordinal id.
    ordinal: u32,
    /// The key of its Seq-id handle in NCBI's map of the batch ([`prefetch_load_key`]), by which
    /// that map releases its blobs (the largest first); `None` where LOSAT does not know the map
    /// ([`LoadPlace::many`]).
    order: Option<u64>,
    /// Whether NCBI surely makes it and loads a blob for it: a place before the failure that LOSAT
    /// knows ([`LoadBeforeStop::Yes`]), a sequence that its data loader gives for the Seq-align id
    /// ([`prefetch_reads`]), in a call with no calls of the run before it, whose scope LOSAT has
    /// not followed ([`TaxdbRun::unrecorded_calls`]).
    surely: bool,
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2326-2369
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2374-2396
// A handle that the scope resolved before asks no data source; every other one locks a blob
// there (taking it out of the cache of released blobs), which the failure releases into that
// cache unless another lock holds it.
/// How a load of a prefetch whose lookup of a Seq-align id failed meets NCBI's data source
/// ([`VolumeLookups::failed_release`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FailedRelease {
    /// The scope resolves its id with a blob that it holds: the load asks no data source, and
    /// nothing is released.
    Scope,
    /// NCBI makes it through the data source, which releases the blob at the failure.
    Released,
    /// NCBI makes it through the data source or does not make it at all, and holds no blob of it
    /// before (a first load).
    Unheld,
    /// Any of these.
    Maybe,
}

// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2000-2026
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:2326-2369
// The failure releases the locks of the data source's map from the largest key down (cited at
// `VolumeLookups::release_failed`), each into the cache of released blobs at its end.
/// Of the loads `made` of a prefetch whose lookup of a Seq-align id failed (each met as `releases`
/// says, [`FailedRelease`]), how many NCBI's data source releases after that of `made[at]`: for
/// certain (released, of a smaller key of a map that LOSAT knows) and at most (any that the scope
/// does not resolve, of a key not larger or not known).
fn released_after(made: &[FailedLoad], releases: &[FailedRelease], at: usize) -> (usize, usize) {
    let own = made[at].order;
    let mut surely = 0;
    let mut most = 0;
    for (index, (other, &release)) in made.iter().zip(releases).enumerate() {
        if index == at || release == FailedRelease::Scope {
            continue;
        }
        match (own, other.order) {
            (Some(own), Some(key)) => {
                surely += usize::from(key < own && release == FailedRelease::Released);
                most += usize::from(key <= own);
            }
            _ => most += 1,
        }
    }
    (surely, most)
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
    /// The run of the application `app` (NCBI's lower-case name), whose posts all carry its
    /// prefix, also where the search it runs is another program's.
    pub(crate) fn running(app: &'static str) -> Self {
        TaxdbRun {
            app: Some(app),
            ..TaxdbRun::default()
        }
    }

    /// The prefix of the run's posts: the running application's name, else `program`'s in lower
    /// case.
    pub(crate) fn prefix(&self, program: &str) -> String {
        self.app
            .map_or_else(|| program.to_ascii_lowercase(), str::to_owned)
    }

    /// Back to the state of a run that has looked nothing up, for a search of the batches before
    /// again (BLASTP's `ErrorAfterBatchReports`: a later batch failed, and the batches before it are
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
            keeps_scope: record.keeps_scope,
            ..VolumeLookups::default()
        };
        drop(record);
        if let Some(rows_open) = self.rows_open.get() {
            rows_open.lock().expect("the rows' texts").clear();
        }
        *self.title_lookups.lock().expect("the title lookups") = TitleLookups::default();
        self.load_texts.clear();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2310-2313
    /// Records whether NCBI's formatter keeps its scope across query batches (`RunThrow::
    /// keeps_scope`: BLAST XML, outfmt 5, where `CBlastFormat::ResetScopeHistory` leaves the scope
    /// intact). Set once where the program opens its database.
    pub(crate) fn set_keeps_scope(&mut self, keeps: bool) {
        self.volume_lookups
            .get_mut()
            .expect("the volume lookups")
            .keeps_scope = keeps;
    }

    /// Records whether the run's [`crate::report::db_subject::ScopeRecord`] keeps earlier batches' scopes across batches
    /// (`FormatThrow::record_keeps_scope`: `Some(false)` drops earlier batches on every new batch;
    /// `Some(true)` or `None` keeps them for outfmt 5 or unknown/multi-format calls).
    pub(crate) fn set_record_keeps_scope(&mut self, keeps: Option<bool>) {
        self.loaded
            .write()
            .expect("the run's record")
            .set_keeps_scope(keeps);
    }

    /// Records that the run's first format is not XML2/JSON (`!FormatThrow::is_xml2`), so any
    /// secondary XML2/JSON lookup on a cache-colliding volume (`read_reported_taxonomy`) is
    /// refused (`TaxdbRun::refuse_cache_collision`), while non-collision runs still share the
    /// run's taxonomy handle, open state, volume lookups, and `calls_before`.
    pub(crate) fn set_separate_xml2(&mut self, separate: bool) {
        self.separate_xml2 = separate;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2638-2653
    /// Mirrors `CBlastFormat::ResetScopeHistory` at the end of a query batch (or of a query's
    /// rounds in psiblast, `psiblast_app.cpp:513`): clears live-scope state on recorded blobs
    /// unless the run keeps its scope (`set_keeps_scope`). Returns without locking when the run
    /// has not recorded an offset-error taxid.
    pub(crate) fn reset_scope(&self) {
        if self
            .texts
            .get()
            .is_none_or(|texts| texts.offset_error_taxid.is_none())
        {
            return;
        }
        self.volume_lookups
            .lock()
            .expect("the volume lookups")
            .reset_scope();
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1785
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
    // A read resolves its id in the scope, where a later prefetch finds it without the data source.
    /// Records that NCBI's live scope may hold a blob that the run's records of the slot do not
    /// name (`VolumeLookups::unrecorded_reads`): a read of the reports, of a prefetch or of
    /// psiblast's PSSM, of any sequence, until the scope's reset ([`TaxdbRun::reset_scope`]).
    fn note_scope_read(&self) {
        self.volume_lookups
            .lock()
            .expect("the volume lookups")
            .unrecorded_reads = true;
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
    /// batch's, or a psiblast round's; [`PrefetchSubject`], in any order), for the walks of the
    /// later reports' loads ([`TaxdbRun::loads_before`]) and the stop of a later batch's prefetch
    /// ([`stop_lookup_unsure`]): the loader keeps the Seq-ids of every sequence that it loaded
    /// (`RegisterIds`), the scope the Bioseqs of the batch (of the run in outfmt 5). A batch's
    /// loads are every subject of its Seq-aligns where its prefetch is on, else the subjects that
    /// the reports NCBI writes (after their cut to `-max_target_seqs`; not an outfmt 5 report that
    /// its XML stream cut before its hits) and psiblast's PSSM read ([`batch_loads`],
    /// [`LoadedBefore::subjects`]; review 3 of line DBID loadord, H2; review 6, H1 and H2). A
    /// caller that searches the queries of a run in several calls, or writes each batch's reports
    /// before it searches the next (BLASTX), records each batch that does not stop; the callers
    /// that hold a call's batches give the loads of the batches before to the stop itself
    /// ([`batch_prefetch`]'s `before`). Only the sequences with parsed Seq-ids are kept (one
    /// without them holds no id that a load can find), each as `ScopeRecord::note` keeps it; the
    /// callers ask it only of a database with parsed Seq-ids. `batch`: the query batch of these
    /// loads, as the program numbers its batches for its walks of the reports' loads
    /// (`BatchPrefetch::report_read`; psiblast's rounds of a query, its batch, share the number),
    /// whose scope holds them until a record of another batch ([`TaxdbRun::loads_before`]).
    pub(crate) fn note_loaded<'s, I>(&self, batch: usize, subjects: impl FnOnce() -> I)
    where
        I: IntoIterator<Item = PrefetchSubject<'s>>,
    {
        self.keep_loaded(Some(batch), subjects());
    }

    /// The record of [`TaxdbRun::note_loaded`], also at a stop ([`batch_prefetch`]): the subjects
    /// that the call's batches before loaded (`batch` `None`: the loader's alone; the run ends at
    /// the stop).
    fn keep_loaded<'s>(
        &self,
        batch: Option<usize>,
        subjects: impl IntoIterator<Item = PrefetchSubject<'s>>,
    ) {
        let parsed: Vec<&crate::report::db_subject::ParsedSubject> = subjects
            .into_iter()
            .filter_map(|subject| subject.parsed)
            .collect();
        if parsed.is_empty() {
            return;
        }
        self.loaded
            .write()
            .expect("the run's record")
            .note(batch, parsed);
    }

    /// Records that calls of the run before this one loaded sequences that no record here holds
    /// where `unrecorded` (`api::run_batches`' calls after the first, each a query batch searched
    /// with a database handle of its own: the call's first query is not the run's first,
    /// `first_query_ordinal` above 0), so that every sequence of the database may have been loaded
    /// before the call: the one owner of that fact for the walks of the reports' loads
    /// ([`TaxdbRun::loads_before`], [`TaxdbRun::calls_before`]) and the stops of the prefetch
    /// ([`stop_lookup_unsure`]). Each program's database sets it where it opens.
    pub(crate) fn note_calls_before(&mut self, unrecorded: bool) {
        self.calls_before = unrecorded;
    }

    /// Whether calls of the run before this one loaded what no record here holds
    /// ([`TaxdbRun::note_calls_before`]).
    pub(crate) fn unrecorded_calls(&self) -> bool {
        self.calls_before
    }

    /// The loads before a walk of the reports' loads that takes no record of this call's own
    /// (TBLASTN and TBLASTX, which keep none where a lookup can fail, and the walks of the batches
    /// that the call still holds at a stop): those of the calls before
    /// ([`TaxdbRun::note_calls_before`]).
    pub(crate) fn calls_before(&self) -> crate::report::db_subject::LoadsBefore {
        crate::report::db_subject::LoadsBefore::after_calls(self.calls_before)
    }

    /// The sequences whose blobs the run's record registered, by ordinal id, in order (tests).
    #[cfg(test)]
    pub(crate) fn loaded(&self) -> Vec<u32> {
        self.loaded
            .read()
            .expect("the run's record")
            .registered_sequences()
    }

    /// The loads before a walk of the reports' loads whose first report is of the query batch
    /// `batch` (`db_subject::LoadsBefore`, the one producer): the run's record as it is now, whose
    /// last batch's scope the reports share where it is `batch` (psiblast's rounds before of the
    /// query; `None`: the record's last batch), after the calls of the run before
    /// ([`TaxdbRun::note_calls_before`]).
    pub(crate) fn loads_before(
        &self,
        batch: Option<usize>,
    ) -> crate::report::db_subject::LoadsBefore {
        crate::report::db_subject::LoadsBefore::new(&self.loaded, batch, self.calls_before)
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

    /// Whether a load of a sequence with a positive taxid can write a text or throw in this run
    /// ([`PrefetchTexts`] of a damaged taxonomy database, [`TaxdbRun::load_throws`]): where NCBI's
    /// scope finds a query for a reported sequence it loads nothing of that sequence, which LOSAT's
    /// loads do not follow, so such a run keeps LOSAT's refusal of it
    /// (`report::query_bioseq::note_or_refuse`).
    pub(crate) fn loads_write_texts(&self, db: &SeqDb) -> bool {
        self.prefetch_texts(db) != PrefetchTexts::default() || self.load_throws(db)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1306-1324, 1455-1466
    // A load looks each positive taxid of its definition lines up after it reads the header,
    // before it decodes it (`x_GetTaxonomy`), through the volume's 256-slot cache (cited at
    // `SlotHolds::walk`): the run's first lookup opens the database (its texts of opening), a
    // lookup of the offset error's taxid that misses its slot writes that error, and every lookup
    // faults or throws where `taxdb.btd` is mapped as a null pointer or a file throws.
    /// Whether a load of the sequences `loads` (ordinal ids in `db`: the loads of a batch's
    /// prefetch, those that NCBI's data source may hold included) may write a text of the run's
    /// taxonomy database or throw, after the loads of the batches before
    /// ([`TaxdbRun::note_title_loads`]): a load of a sequence with a positive taxid (or taxids
    /// LOSAT cannot read) where every lookup can (a null mapping, a file that throws), or before
    /// the run's first lookup where opening the database writes texts, and a load with the offset
    /// error's taxid in a volume where no load before looked it up, or where another taxid of the
    /// volume can take its slot ([`TaxdbRun::cache_collision`]).
    fn title_loads_write_texts(&self, db: &SeqDb, loads: impl IntoIterator<Item = u32>) -> bool {
        if !self.loads_write_texts(db) {
            return false;
        }
        let texts = self.prefetch_texts(db);
        let every = texts.maps_null || self.load_throws(db);
        let collision = self.cache_collision(db).is_some();
        let before = self.title_lookups.lock().expect("the title lookups");
        loads.into_iter().any(|ordinal| {
            let taxids = db.defline_taxids(ordinal).ok();
            if !taxids
                .as_ref()
                .is_none_or(|taxids| taxids.iter().any(|&taxid| taxid > 0))
            {
                return false;
            }
            every
                || (texts.opening && !before.opened)
                || texts.offset_error_taxid.is_some_and(|last| {
                    taxids.as_ref().is_none_or(|taxids| taxids.contains(&last))
                        && (collision
                            || db
                                .volume_index(ordinal)
                                .map_or(true, |volume| !before.offset_volumes.contains(&volume)))
                })
        })
    }

    /// Notes the subjects of a batch's prefetch `ordinals` (every one of its Seq-aligns: each was
    /// loaded by this batch or a batch before) for the batches after it
    /// ([`TaxdbRun::title_loads_write_texts`]): a positive taxid opens the database, the offset
    /// error's taxid fills its slot of the volume. Nothing where no load writes a text.
    pub(crate) fn note_title_loads(&self, db: &SeqDb, ordinals: impl IntoIterator<Item = u32>) {
        if !self.loads_write_texts(db) {
            return;
        }
        let last = self.prefetch_texts(db).offset_error_taxid;
        let mut seen = self.title_lookups.lock().expect("the title lookups");
        for ordinal in ordinals {
            let Ok(taxids) = db.defline_taxids(ordinal) else {
                continue;
            };
            if taxids.iter().any(|&taxid| taxid > 0) {
                seen.opened = true;
            }
            if let (Some(last), Ok(volume)) = (last, db.volume_index(ordinal)) {
                if taxids.contains(&last) {
                    seen.offset_volumes.insert(volume);
                }
            }
        }
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
    /// `(path of taxdb.bti, reason)` for the taxonomy cache collision, where a batch's prefetch
    /// ([`prefetch_slot_plan`]) and the reports' reads ([`load_read_texts`]) follow NCBI's order,
    /// the XML2 and JSON reports' too (they prefetch as the other formats,
    /// [`prefetches_subjects`]), and refuse what that order does not decide.
    fn cache_collision(&self, db: &SeqDb) -> Option<&(std::path::PathBuf, String)> {
        self.collision
            .get_or_init(|| {
                let last = self.prefetch_texts(db).offset_error_taxid?;
                let (volume, other) = taxonomy_cache_collision(db, last)?;
                let path = db.search_path()?.resolve_file(b"taxdb.bti")?;
                let volume_name = db
                    .volumes()
                    .get(volume)
                    .map_or_else(String::new, |vol| vol.base().display().to_string());
                Some((
                    crate::utils::os_bytes::path_from_bytes(path),
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
    /// Such a run refuses a batch's prefetch or a report's read (in every format) where NCBI's
    /// order of its loads, or its keeping of a sequence read before, decides the texts ([`prefetch_slot_plan`], [`load_read_texts`]), which can be a later batch: BLASTX holds its XML document (outfmt 5, and 13 to 16) then, which
    /// the refusal drops; TBLASTX turns its report stream off and BLASTP takes the search of all
    /// queries, which make every batch's lookups before the first byte.
    pub(crate) fn cache_may_refuse(&self, db: &SeqDb) -> bool {
        self.cache_collision(db).is_some()
    }

    /// LOSAT's refusal of a fetch of the run where NCBI's taxonomy cache would repeat the offset
    /// error ([`TaxdbRun::cache_collision`]), as a reader's refusal ([`TaxonomyRefusal`]) that names
    /// `program`; `Ok` for any other run. Used for secondary XML2/JSON taxonomy reads
    /// (`separate_xml2`, [`read_reported_taxonomy`]).
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
    /// ([`prefetch_slot_plan`]), recorded there. `by`: the subject whose Seq-align id asks for
    /// the sequence (NCBI's blob of that id, [`blob_name`]).
    #[allow(clippy::too_many_arguments)]
    fn load(
        &self,
        db: &SeqDb,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        ordinal: u32,
        by: u32,
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
            (Some(last), None) => {
                let name = if lookups
                    .iter()
                    .any(|&taxid| tax_cache_slot(taxid) == tax_cache_slot(last))
                {
                    blob_name(db, ordinal, by, program)?
                } else {
                    BlobName::default()
                };
                Some(
                    self.slot_misses(db, volume, ordinal, &name, lookups, last)
                        .map_err(|cause| {
                            taxdb_error(slot_uncertain_error(db, volume, last, cause), program)
                        })?,
                )
            }
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

    /// Which of the lookups `lookups` of a load of the sequence `ordinal` by the blob name `name`
    /// in the volume `volume` ask the taxonomy database again for the offset error's taxid `last`
    /// ([`VolumeLookups::offset_slot`]).
    fn slot_misses(
        &self,
        db: &SeqDb,
        volume: usize,
        ordinal: u32,
        name: &BlobName,
        lookups: &[i32],
        last: i32,
    ) -> std::result::Result<Vec<bool>, SlotUncertain> {
        self.volume_lookups(db)
            .offset_slot(volume, ordinal, name, lookups, last)
    }

    /// The run's [`VolumeLookups`] on `db`, whose `kept` and `in_scope_kept` are NCBI's for `db`
    /// ([`keeps_every_load`], `!reads_may_load_again`; never in a sequential run, whose reads are
    /// all loads).
    fn volume_lookups(&self, db: &SeqDb) -> std::sync::MutexGuard<'_, VolumeLookups> {
        let mut record = self.volume_lookups.lock().expect("the volume lookups");
        record.kept = !record.sequential && keeps_every_load(db);
        record.in_scope_kept = !record.sequential && !reads_may_load_again(db);
        record
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

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:218-222 (one
    // `CTaxDBFileInfo` per process, made at the first `GetTaxNames`)
    /// A sequential run ([`TaxdbRun::sequential`]) of another `CSeqDB` of the process of this run:
    /// the taxonomy database is the one this run opened (NCBI's one per process), so its texts of
    /// opening are not written again; the volumes' caches are the new database's, whose lookups
    /// write their texts again (NCBI 2.17.0, windowmasker `-mk_counts -infmt blastdb`, line RG
    /// T8 survey W1: the open error once in three passes, the offset error in each).
    pub(crate) fn sequential_after(&self) -> Self {
        let run = TaxdbRun::sequential();
        if let Some(lazy) = self.lazy.get() {
            let _ = run.lazy.set(std::sync::Arc::clone(lazy));
        }
        if let Some(&texts) = self.texts.get() {
            let _ = run.texts.set(texts);
        }
        if let Some(&found) = self.found.get() {
            let _ = run.found.set(found);
        }
        run
    }

    /// [`TaxdbRun::sequential`] of one `CSeqDB` of a program that opens several, each with its
    /// volumes' caches, while the taxonomy database `taxdb` is the process's (NCBI's static
    /// `CTaxDBFileInfo`): blastdbcheck opens each volume as a database of its own. The texts of
    /// opening `taxdb` come with the first lookup of any of them.
    pub(crate) fn sequential_with(taxdb: std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>) -> Self {
        let run = Self::sequential();
        let _ = run.lazy.set(taxdb);
        run
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
                // A read by the OID (`CSeqDB::GetBioseq`), not by a Seq-id handle: one name for
                // each sequence.
                let name = BlobName {
                    handle: Some([b"oid".as_slice(), &oid.to_le_bytes()].concat()),
                    ..BlobName::default()
                };
                Some(
                    self.slot_misses(db, volume, oid, &name, &lookups, last)
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
        (ordinal, by): (u32, u32),
        lookups: &[i32],
        planned: Option<&[bool]>,
        lower: &str,
        program: &str,
    ) -> Result<(Vec<u8>, usize)> {
        let mut text = String::new();
        let mut opening = 0;
        for (_, messages) in self.load(db, taxdb, ordinal, by, lookups, planned, program)? {
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

    /// [`TaxdbRun::load`] of the sequence `ordinal` by the Seq-align id of the subject `by`, its
    /// texts appended to `text` as `ERR_POST` lines of the program `lower`.
    #[allow(clippy::too_many_arguments)]
    fn volume_lookups_of(
        &self,
        db: &SeqDb,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        (ordinal, by): (u32, u32),
        lookups: &[i32],
        lower: &str,
        text: &mut String,
        program: &str,
    ) -> Result<()> {
        for (_, messages) in self.load(db, taxdb, ordinal, by, lookups, None, program)? {
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
    let volume_name = volume_name(db, volume);
    let path = taxdb_bti_path(db);
    if let Some(order) = cause.order {
        return SeqDbError::Unsupported {
            path,
            message: slot_order_message(db, &volume_name, last, cause, order),
        };
    }
    if cause.unnamed {
        return SeqDbError::Unsupported {
            path,
            message: format!(
                "taxid {last}, the last entry of the taxonomy database, begins past the end of \
                 taxdb.btd, and the run reads database sequence {} of the volume {volume_name}, \
                 which holds it or another taxid of its slot of NCBI BLAST+'s 256-slot taxonomy \
                 cache of a volume (CSeqDBVol::m_TaxCache, CSeqDBIntCache), by a Seq-align id \
                 whose Seq-id handle LOSAT does not know (a text id with neither accession nor \
                 name, for which NCBI makes a new handle and so a new load at each read, or a giim \
                 id): those loads decide how often NCBI writes \"{}\", which LOSAT does not \
                 reproduce",
                db.oid(cause.ordinal),
                crate::seqdb::taxdb::MSG_OFFSET_ERROR
            ),
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
            maybe: true,
            ..
        }) => format!(
            "after the prefetch of a query batch, whose lookup of a Seq-align id failed, loaded \
             database sequence {} with taxid {taxid} before that failure or not at all",
            db.oid(ordinal)
        ),
        Some(SlotFetch {
            ordinal,
            taxid,
            again: false,
            ..
        }) => format!(
            "after the run's fetch of database sequence {} with taxid {taxid}",
            db.oid(ordinal)
        ),
        Some(SlotFetch {
            ordinal,
            taxid,
            again: true,
            ..
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
        Some(fetch) if fetch.maybe => format!(
            "there or not as that prefetch loaded database sequence {} or not",
            db.oid(fetch.ordinal)
        ),
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
/// NCBI's streams write it).
///
/// NCBI's messages about the taxonomy database go to standard error, and not here: those of the
/// loads of the reported subjects (the lookup that the volume makes for each taxid of a fetched
/// subject, `CSeqDBVol::x_GetTaxonomy`, once per volume and taxid) are the batch's prefetch's
/// ([`batch_prefetch`]) or, where it is off, the reads of the `Hit`s ([`subject_read_texts`]),
/// which the caller models before this call; the texts of opening the database where a `Hit`'s
/// own lookup opens it ([`TaxdbRun::open_for_rows`], [`RowLookups::Taxids`]) are in the returned
/// [`TaxdbNotes`], which writes them at the first `Hit` that asks, and those of the lookup of each
/// `Hit` are in `taxonomy` and written at the hit ([`TaxdbPrinter::at_hit`]). `program` is the
/// program that refuses, as LOSAT's refusals name it (`BLASTP`); its messages name it in lower
/// case.
///
/// The database's reported sequences can have parsed Seq-ids (`parsed[rank]`, empty where none
/// has): the third vector holds, by rank, what the taxonomy database said about the taxid of each
/// definition line of such a sequence (empty where nothing was asked), as NCBI's `HitDescr` loop
/// asks it (`SubjectNames::of_parsed`): `Unnamed` for a definition line without a taxid, and no
/// answer after the first `NotFound` nor from the first definition line without a Seq-id on
/// (NCBI's loop ends there). The volume's lookups are those of
/// `ParsedSubject::volume_taxids`.
///
/// On a database whose volume holds the offset error's taxid and another taxid of its slot
/// ([`TaxdbRun::cache_collision`]), a secondary XML2/JSON format after a non-XML2 first format
/// (`TaxdbRun::set_separate_xml2`) is refused (`TaxdbRun::refuse_cache_collision`); otherwise those
/// loads are the caller's model too: a batch's prefetch in the order of its Seq-id handles
/// ([`prefetch_slot_plan`]), or with the prefetch off each `Hit`'s read in the report's order
/// ([`load_read_texts`]), as for the other formats, which refuse what that order does not decide.
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
    if run.separate_xml2 && run.cache_collision(db).is_some() {
        run.refuse_cache_collision(db, program)?;
    }
    let taxdb = run.taxdb(search_path);
    // NCBI opens the database at its first lookup: a load's (the texts are the load's, modelled
    // before this call) or a `Hit`'s, whose texts the caller's model of the batches kept
    // ([`TaxdbRun::open_for_rows`]); the texts are written at the first `Hit` that asks. The open
    // here makes texts only where that model opened nothing (no `Hit` asks).
    let mut open = std::mem::take(&mut *run.rows_open().lock().expect("the texts of opening"));
    open.extend(taxdb.open().map_err(|error| taxdb_error(error, program))?);
    let notes = TaxdbNotes {
        program: run.prefix(program),
        open,
    };
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
    for &rank in reported {
        let subject = parsed_subject(rank);
        match subject {
            Some(subject) if defline_taxids(subject).iter().all(|&taxid| taxid == 0) => continue,
            None if taxids[rank] == 0 => continue,
            _ => {}
        }
        let Some(subject) = subject else {
            taxonomy[rank] = hit_lookup(taxids[rank])?;
            continue;
        };
        let mut answers = Vec::new();
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:370-371
        // The loop throws at a line without a Seq-id before it looks its taxid up.
        let ends = subject
            .deflines()
            .iter()
            .position(|defline| defline.ids.is_empty());
        for taxid in defline_taxids(subject)
            .into_iter()
            .take(ends.unwrap_or(usize::MAX))
        {
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
pub struct TaxdbSetting(Option<Vec<u8>>);

impl TaxdbSetting {
    /// The value of `settings`.
    pub(crate) fn of(settings: &ApplicationSettings) -> Self {
        TaxdbSetting(settings.blastdb.clone())
    }
}

/// How a `-outfmt` value consults NCBI's taxonomy database (`taxdb.bti` / `taxdb.btd`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OutfmtTaxdbUse {
    /// The format neither checks `taxdb.bti` presence nor reads taxonomy names for its output
    /// (outfmt 5, 8, 9, 11, 12, 17, and 6/7/10/20 without taxonomy name fields — including
    /// `staxid` and `staxids`).
    None,
    /// The format reads taxonomy names only when `taxdb.bti` already exists on the search path
    /// (tabular singular name fields `ssciname`/`scomname`/`sblastname`/`sskingdom`, organism
    /// report 18, and pairwise/query-anchored 0..=4 for titleless database sequences with a
    /// taxid): creating a missing `taxdb.bti` when none is on the search path leaves CLI output
    /// unchanged, whereas shadowing an existing `taxdb.bti` (or creating its paired
    /// `taxdb.btd`) changes it.
    WhenPresent,
    /// The format checks `taxdb.bti` presence even when no `taxdb.bti` is on the search path,
    /// so creating a missing `taxdb.bti` changes CLI output or stderr:
    /// - XML2/JSON formats 13, 14, 15, 16 on `-db` (`Xml2Kind::of`, `x_InitCommon`);
    /// - Tabular/CSV formats 6, 7, 10, 20 with plural taxonomy name fields
    ///   (`TabularSpec::checks_taxdb`, `x_CheckTaxDB`), for both `-db` and `-subject`.
    Always,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-135
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:582-606
// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xml2format.cpp:95-110
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:380-386
/// Parsed `-outfmt` report format and optional tabular specification for `app`, shared by
/// [`outfmt_taxdb_use`] and [`outfmts_need_subject_taxdb`].
fn parse_outfmt_for_taxdb(
    outfmt: impl AsRef<[u8]>,
    app: crate::report::archive::options::RequestKind,
) -> Option<(
    crate::blastinput::app::ReportFormat,
    Option<crate::report::tabular::TabularSpec>,
)> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    let report = crate::blastinput::app::ReportFormat::from_number(choice.number)?;
    let tabular = crate::report::tabular::TabularSpec::of_app(report, &choice, app);
    Some((report, tabular))
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:522-532
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
/// How an `-outfmt` string consults NCBI's taxonomy database for `app` (and whether the search
/// uses `-db`, `is_db`), sharing the engine's own `parse_formatting_string`, `Xml2Kind::of`,
/// `ReportFormat::writes_html_links`, `TabularSpec::of_app` and `PrintedIds::of` rules
/// (line WQA D19b).
///
/// On `-db` searches, tabular formats (6, 7, 10, 20) that print a subject id column
/// (`PrintedIds`: `sseqid`, `sacc`, `saccver`, `sallseqid`, `sallacc`, including the default
/// `std` fields of outfmt 6 and 7) derive the replacement word for `gnl|BL_ORD_ID|n` on
/// databases without `-parse_seqids` from `GenerateDefline`, which uses `taxdb` when a sequence's
/// stored title is empty, so they return [`OutfmtTaxdbUse::WhenPresent`].
pub(crate) fn outfmt_taxdb_use(
    outfmt: impl AsRef<[u8]>,
    app_kind: crate::report::archive::options::RequestKind,
    is_db: bool,
) -> OutfmtTaxdbUse {
    use crate::blastinput::app::ReportFormat;
    use crate::report::tabular::PrintedIds;
    let Some((report, tabular)) = parse_outfmt_for_taxdb(outfmt, app_kind) else {
        return OutfmtTaxdbUse::None;
    };
    if is_db && crate::report::blast_xml2::Xml2Kind::of(report).is_some() {
        return OutfmtTaxdbUse::Always;
    }
    if let Some(spec) = tabular {
        if spec.checks_taxdb() {
            return OutfmtTaxdbUse::Always;
        }
        if is_db && (spec.needs().taxonomy_names || PrintedIds::of(Some(&spec)).any()) {
            return OutfmtTaxdbUse::WhenPresent;
        }
        return OutfmtTaxdbUse::None;
    }
    if is_db
        && (report.writes_html_links()
            || (matches!(
                report,
                ReportFormat::Tabular | ReportFormat::TabularWithComments
            ) && PrintedIds::of(None).any()))
    {
        return OutfmtTaxdbUse::WhenPresent;
    }
    OutfmtTaxdbUse::None
}

/// Whether any format in `outfmts` looks every taxid up for a `-subject` search (`x_CheckTaxDB`
/// plural names or outfmt 18, `TabularSpec::looks_up_every_taxid`), shared by
/// `api::run::taxdb_setting`.
pub(crate) fn outfmts_need_subject_taxdb<O: AsRef<[u8]>>(
    outfmts: impl IntoIterator<Item = O>,
    kind: crate::report::archive::options::RequestKind,
) -> bool {
    outfmts.into_iter().any(|outfmt| {
        parse_outfmt_for_taxdb(outfmt, kind)
            .and_then(|(_, spec)| spec)
            .is_some_and(|spec| spec.looks_up_every_taxid())
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:756-772
/// Whether a `-db` search with `-outfmt` `outfmt` prefetches reported subjects (`prefetch_messages`,
/// outfmt 0..=7 and 17): when `taxdb.btd` already exists on disk in a search directory, creating a
/// `taxdb.bti` log file beside it causes `TaxDb::open` to read that file during prefetch.
pub(crate) fn outfmt_prefetches_taxdb(outfmt: impl AsRef<[u8]>) -> bool {
    let Ok(choice) = crate::blastinput::app::parse_formatting_string(outfmt) else {
        return false;
    };
    matches!(choice.number, 0..=7 | 17)
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
/// What the tabular reports `specs` of a `-subject` search ask of the
/// taxonomy database where NCBI makes their formatters: whether `taxdb.bti` is on the search path
/// of `setting` (`x_CheckTaxDB`'s `SeqDB_ResolveDbPath` of `taxdb.bti`: the current directory,
/// `BLASTDB` and the registry's `[BLAST] BLASTDB`, as a database is looked for, `seqdb::SearchPath`)
/// for the plural names, which warn without it, and the organism report, whose lookups open it
/// (`SubjectTaxdb`). A report of singular names reads nothing: NCBI looks a singular name up by
/// the subject's taxid, and the records of `-subject` have none. `setting` is `None` in LOSAT Web:
/// not found.
///
/// The command line and the library's searches of records call it at the same point, after the
/// options are processed (each program's `open_search_database`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-134
pub(crate) fn subject_taxonomy<'s>(
    setting: Option<&TaxdbSetting>,
    specs: impl IntoIterator<Item = &'s crate::report::tabular::TabularSpec>,
) -> SubjectTaxonomy {
    let (mut looks_up, mut tax_report) = (false, false);
    for spec in specs {
        looks_up |= spec.looks_up_every_taxid();
        tax_report |= spec.kind == crate::report::tabular::TabularKind::TaxReport;
    }
    if !looks_up {
        return SubjectTaxonomy::default();
    }
    let Some(setting) = setting else {
        return SubjectTaxonomy {
            taxdb: None,
            found: Some(false),
        };
    };
    let path = SearchPath::from_environment(setting.0.as_deref());
    let found = path.resolve_file(b"taxdb.bti").is_some();
    // The organism report's lookups open the run's own database.
    let taxdb = (tax_report && found).then(|| {
        SubjectTaxdb(std::sync::Arc::new(crate::seqdb::taxdb::LazyTaxDb::new(
            path,
        )))
    });
    SubjectTaxonomy {
        taxdb,
        found: Some(found),
    }
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
            &run.prefix(program),
        )));
    }
    let names = crate::report::tabular::TabularNames::look_up(
        spec,
        &run.taxdb(search_path),
        taxonomies,
        &run.prefix(program),
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
/// reports, for the formats that prefetch ([`prefetches_subjects`]): the texts of opening it and of the lookup that a volume makes the first time it
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
        RowLookups::None,
        program,
    )?;
    Ok(reads.into_iter().flat_map(|read| read.text).collect())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1110-1113
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (cited at
// `loads_before_failure`)
// `LoadSequencesToScope` loads the batch's Seq-align ids with one `GetBioseqHandles`, in the order
// of the map of their Seq-id handles; the texts of a set of loads do not depend on their order
// but where a load posts NCBI's warning of a UTF8String header (`scope_utf8_prefetch_loads`).
/// The loads of a batch's prefetch of `subjects` (each `(ordinal, parsed Seq-ids)`;
/// [`subject_read_texts`] with no key): [`batch_prefetch`]'s. Where the load of a subject's
/// Seq-align id posts NCBI's warning of a UTF8String header (`ParsedSubject::scope_utf8`,
/// `load_read_texts`), the loads come in NCBI's order ([`scope_utf8_prefetch_loads`]); else in
/// the order of `subjects` (the bytes of a set of loads do not depend on their order:
/// [`prefetch_messages`]). `at_posts`: the loads' texts are made again where the batch's posts are
/// written (`diag_log::Posts::loads`), where a load that asks the taxonomy database takes the texts
/// that a load NCBI does not make left (`ReadLoad::placed`). Each read's `subject` is the index of
/// its load in `subjects` (the caller's checks of a load whose place LOSAT does not know). `map`:
/// why LOSAT does not know NCBI's map of the batch ([`batch_map`]).
fn prefetch_loads(
    db: &SeqDb,
    run: &TaxdbRun,
    subjects: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    at_posts: bool,
    order: LoadsOrder<'_>,
    map: Option<&'static str>,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    if db.scope_utf8_loads() && scope_utf8_among(subjects) {
        return scope_utf8_prefetch_loads(db, run, subjects, at_posts, order, map, program);
    }
    load_read_texts(
        db,
        run,
        subjects
            .iter()
            .enumerate()
            .map(|(index, &(ordinal, parsed))| (index as u32, ordinal, SubjectLoad::Read(parsed))),
        RowLookups::None,
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
    map: Option<&'static str>,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    let places = prefetch_places(db, subjects, map);
    let mut order: Vec<usize> = (0..subjects.len()).collect();
    order.sort_by_key(|&index| places[index].order);
    let reads = load_read_texts(
        db,
        run,
        order.iter().map(|&index| {
            (
                index as u32,
                subjects[index].0,
                SubjectLoad::Read(subjects[index].1),
            )
        }),
        RowLookups::None,
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
    Ok(reads)
}

/// Why the order of two loads of NCBI's prefetch whose Seq-id handles have the same packed number
/// (or are both unpacked) is not known ([`scope_utf8_prefetch_loads`]).
const ADDRESS_ORDER: &str = "their Seq-id handles have the same packed number, or neither is packed, and NCBI's map orders such handles by the address of their info";

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3279-3300, 3345-3368
// NCBI reference (598d8ae6): c++/src/objmgr/seq_id_sort.cpp:47-121, 144-163
// NCBI reference (598d8ae6): c++/include/objmgr/impl/seq_id_sort.hpp:100
// The prefetch sorts the batch's Seq-id handles by `CSortableSeq_id` (a general id by its
// database, then its numeric tag: `gnl|BL_ORD_ID|n` by n, OID 0 first) without repeats; of more
// than 200 it takes pieces of 100 while more than 150 are left, then the rest, and each piece
// loads its subjects in the order of its map of Seq-id handles (by packed number: OID 0 last).
/// The place of each of `ordinals` (a batch's prefetch's loads, a sequence possibly more than
/// once) of a database without parsed Seq-ids, whose Seq-align ids are `gnl|BL_ORD_ID|n`, in
/// NCBI's order of the loads ([`load_places`] in the map of [`batch_map`]); of a batch of more
/// than [`PREFETCH_ONE_MAP_IDS`] ids, the place of its piece first, then its place in the piece's
/// map.
fn ordinal_prefetch_places(db: &SeqDb, ordinals: &[u32]) -> Vec<LoadPlace> {
    let subjects: Vec<PrefetchSubject<'_>> = ordinals
        .iter()
        .map(|&ordinal| PrefetchSubject::unflipped(ordinal, None))
        .collect();
    let map = batch_map(db, &subjects);
    let mut places = load_places(db, &subjects, map);
    if map != Some(MANY_IDS) {
        return places;
    }
    let mut oids: Vec<u32> = ordinals
        .iter()
        .map(|&ordinal| db.oid(ordinal).get())
        .collect();
    oids.sort_unstable();
    oids.dedup();
    // The first sorted index of each piece.
    let mut starts = Vec::new();
    let mut pos = 0;
    while pos < oids.len() {
        starts.push(pos);
        let left = oids.len() - pos;
        pos += if left > 150 { 100 } else { left };
    }
    for (place, &ordinal) in places.iter_mut().zip(ordinals) {
        let index = oids
            .binary_search(&db.oid(ordinal).get())
            .expect("every OID is sorted");
        let piece = starts.partition_point(|&start| start <= index) - 1;
        // The handle's place in its piece's map: `gnl|BL_ORD_ID|n` at n - 1, OID 0 last.
        let in_piece = place.order.min(u64::from(u32::MAX));
        place.order = ((piece as u64) << 32) | in_piece;
        place.many = None;
        place.unknown = None;
    }
    places
}

/// The place of each of `subjects` (a batch's prefetch's loads) in NCBI's map of Seq-id handles
/// ([`LoadPlace`] without the stops, as [`load_before_stop`] makes it), of a batch whose map LOSAT
/// does not know for the reason `map` ([`batch_map`]; `None`: it knows it).
fn prefetch_places(
    db: &SeqDb,
    subjects: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    map: Option<&'static str>,
) -> Vec<LoadPlace> {
    let prefetched: Vec<PrefetchSubject<'_>> = subjects
        .iter()
        .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed))
        .collect();
    load_places(db, &prefetched, map)
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
/// (`report::query_warnings::SubjectReads`). `rows`: which of the report's rows look their
/// subject's taxid up themselves ([`RowLookups`]), each after it loads its subject, so the first
/// such row opens the taxonomy database after its subject's read, where no read before did.
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
    rows: RowLookups,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    load_read_texts(
        db,
        run,
        subjects
            .into_iter()
            .map(|(key, ordinal, parsed)| (key, ordinal, SubjectLoad::Read(parsed))),
        rows,
        LoadsOrder::Reads,
        program,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:476-509
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1444
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:750-763
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:877-890
/// Filters a single query's one-by-one report subject keys (`ranks`, in hit-list order) to match
/// the `Hit` reads NCBI actually makes:
/// - When `joins_copies` is true (BLAST XML, XML2, JSON), drops adjacent keys that join the
///   preceding `Hit` as the same record (`joins_as_copy`, matching `blast_xml::hit_ranges`).
/// - Drops keys never reached by any report (`loads.reached(key)`).
/// - When `ends_on_throw` is true (`AfterThrow::EndsReport` / `EndsRun`), truncates immediately
///   after the first key whose lookup throws (`loads.thrown(key).is_some()`), so subjects after a
///   throwing `Hit` in that query do not perform phantom taxonomy/slot loads even if a later query
///   reaches them.
pub(crate) fn one_by_one_query_reads<P>(
    ranks: impl IntoIterator<Item = usize>,
    parsed: impl Fn(usize) -> Option<P>,
    loads: &crate::report::db_subject::HitLoads,
    joins_copies: bool,
    ends_on_throw: bool,
) -> Vec<usize>
where
    P: std::ops::Deref<Target = crate::report::db_subject::ParsedSubject>,
{
    let mut out = Vec::new();
    let mut previous: Option<usize> = None;
    for rank in ranks {
        if joins_copies {
            let joins = previous.is_some_and(|before| {
                before == rank
                    || matches!(
                        (parsed(before), parsed(rank)),
                        (Some(earlier), Some(later))
                            if crate::report::db_subject::joins_as_copy(&earlier, &later)
                    )
            });
            previous = Some(rank);
            if joins {
                continue;
            }
        }
        if !loads.reached(rank) {
            continue;
        }
        out.push(rank);
        if ends_on_throw && loads.thrown(rank).is_some() {
            break;
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:819-838 (cited above)
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:380-386 (the `Hit`
// asks `CSeqDB::GetTaxInfo` for each taxid that is not 0)
/// Which rows of a report look their subject's taxid up in the taxonomy database themselves,
/// after the row's read of its subject: the first such row opens the database where no load
/// before it did ([`TaxdbRun::open_for_rows`]), so its texts of opening come there.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum RowLookups {
    /// None: the rows show no name of the taxonomy database.
    #[default]
    None,
    /// Every row, taxid 0 too: a tabular report with a plural name, the organism report.
    Every,
    /// A row whose subject has a taxid that is not 0: the `Hit` of the XML2 and JSON formats
    /// (outfmt 13 to 16; a negative taxid too, which no load asks for).
    Taxids,
}

impl RowLookups {
    /// The rows of the first format `first` that look a taxid up, of a tabular report whose
    /// field list has a plural name (`plural_names`, `TabularSpec::looks_up_every_taxid`, the
    /// organism report included): the plural names, or the XML2 and JSON `Hit`s.
    pub(crate) fn of_first_format(
        first: Option<crate::blastinput::app::ReportFormat>,
        plural_names: bool,
    ) -> Self {
        use crate::blastinput::app::ReportFormat;
        if plural_names {
            RowLookups::Every
        } else if matches!(
            first,
            Some(
                ReportFormat::Json | ReportFormat::Xml2 | ReportFormat::JsonS | ReportFormat::Xml2S
            )
        ) {
            RowLookups::Taxids
        } else {
            RowLookups::None
        }
    }

    /// Whether the row of the subject `ordinal` (`parsed`: its definition lines, where it has
    /// parsed Seq-ids) looks a taxid up.
    pub(crate) fn asks(
        self,
        db: &SeqDb,
        ordinal: u32,
        parsed: Option<&crate::report::db_subject::ParsedSubject>,
        program: &str,
    ) -> Result<bool> {
        Ok(match self {
            RowLookups::None => false,
            RowLookups::Every => true,
            RowLookups::Taxids => match parsed {
                Some(subject) => subject
                    .deflines()
                    .iter()
                    .any(|defline| defline.taxid.is_some_and(|taxid| taxid != 0)),
                // One definition line per reported sequence (`read_reported_taxonomy`).
                None => db
                    .defline_taxids(ordinal)
                    .map_err(|error| database_error(error, program))?
                    .first()
                    .is_some_and(|&taxid| taxid != 0),
            },
        })
    }

    /// The open of the taxonomy database by the rows of the subjects `subjects` (each its ordinal
    /// id and parsed Seq-ids), after their batch's prefetch and before the next batch's
    /// ([`TaxdbRun::open_for_rows`]), where one of the rows asks.
    pub(crate) fn open_after_prefetch<'s>(
        self,
        db: &SeqDb,
        run: &TaxdbRun,
        subjects: impl IntoIterator<Item = (u32, Option<&'s crate::report::db_subject::ParsedSubject>)>,
        program: &str,
    ) -> Result<()> {
        if self == RowLookups::None {
            return Ok(());
        }
        for (ordinal, parsed) in subjects {
            if self.asks(db, ordinal, parsed, program)? {
                return run.open_for_rows(db, program);
            }
        }
        Ok(())
    }
}

/// How NCBI orders the loads of one call of [`load_read_texts`], as far as it decides the texts of
/// the offset error on a database whose volume holds its taxid and another taxid of its slot
/// ([`TaxdbRun::cache_collision`]).
#[derive(Clone, Copy)]
enum LoadsOrder<'p> {
    /// A report's reads, in NCBI's order of the reads: each load follows the slot
    /// ([`VolumeLookups::offset_slot`]; a read of a sequence that NCBI keeps loads nothing). On
    /// such a database, LOSAT's refusal of a load of that slot where a query batch may read more
    /// sequences than NCBI keeps ([`slot_load_checks`]), and where a BLAST XML report that the
    /// stream cuts leaves such a load out (`report::query_warnings::ReadLoad::cut_refusal`).
    Reads,
    /// A batch's prefetch ([`batch_prefetch`], and the loads before a stop, [`prefetch_stop`]):
    /// the misses of each load in the slot of the offset error's taxid that NCBI's order of the
    /// loads gives (`planned`, by the reported subject's ordinal id, [`prefetch_slot_plan`]; `None`
    /// where no volume holds that taxid and another of its slot, whose loads give the same texts
    /// in any order and follow the slot one by one). A load of that slot that `planned` does not
    /// hold (one whose place before a stop LOSAT does not know) is refused.
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
    rows: RowLookups,
    order: LoadsOrder<'_>,
    program: &str,
) -> Result<crate::report::query_warnings::QueryReads> {
    use crate::report::query_warnings::{ReadLoad, SubjectRead};
    let lower = run.prefix(program);
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
    let (planned, report_reads) = match order {
        LoadsOrder::Reads => (None, true),
        LoadsOrder::Prefetch(planned) => (planned, false),
    };
    let taxdb = run.taxdb(search_path);
    let texts = run.load_texts(db);
    let offset_taxid = run.prefetch_texts(db).offset_error_taxid;
    // On a database whose volume holds the offset error's taxid and another of its slot, the
    // taxid whose slot a load's lookups take, if they take it ([`slot_load_checks`]).
    let collision = offset_taxid.filter(|_| run.cache_collision(db).is_some());
    let mut reads = Vec::new();
    // NCBI's blobs read (the sequence with the subject whose Seq-align id reads it, `blob_name`):
    // a load that a stopping id is given instead is another blob than the sequence's own read.
    let mut fetched: std::collections::HashSet<(u32, u32)> = Default::default();
    let mut thrown: std::collections::HashSet<u32> = Default::default();
    let mut row_opened = false;
    let mut scope_read = false;
    for (key, ordinal, load) in loads {
        // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1785
        // The scope keeps the ids that its reads resolve until its reset
        // ([`TaxdbRun::note_scope_read`]).
        if !std::mem::replace(&mut scope_read, true) {
            run.note_scope_read();
        }
        // A subject whose lookup throws (`ParsedSubject::lookup_throws`, a report that reads one by
        // one reaches it: `HitLoads::reached`) has read the candidates' headers first, so its load
        // posts NCBI's warning of a UTF8String among them and nothing else: no blob is read, no
        // taxonomy text comes, and a row of it opens nothing (`rows` waits for the first
        // subject that loads).
        if let SubjectLoad::Read(Some(subject)) = load {
            if subject.lookup_throws().is_some() {
                if thrown.insert(ordinal) && scope_utf8 {
                    reads.extend(scope_utf8_read(db, run, key, subject, &lower));
                }
                continue;
            }
        }
        // A load that a stopping id is given instead is NCBI's blob of that id's handle, not the
        // one of the sequence's own read, at a place that LOSAT does not know: on such a database
        // it is refused where it takes the slot, also after the sequence's own load.
        if let (SubjectLoad::Instead { stop, by }, Some(last)) = (load, collision) {
            instead_slot_check(db, ordinal, stop, by, last, program)?;
        }
        let by = match load {
            SubjectLoad::Read(_) => ordinal,
            SubjectLoad::Instead { stop, .. } => stop,
        };
        if !fetched.insert((ordinal, by)) {
            continue;
        }
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
            Some(subject) => match prefetch_reads(db, ordinal, subject, true, false, program)? {
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
                let slot = collision.filter(|&last| {
                    lookups
                        .iter()
                        .any(|&taxid| tax_cache_slot(taxid) == tax_cache_slot(last))
                });
                let cut_refusal = match slot {
                    Some(last) => {
                        slot_load_checks(db, ordinal, blob, last, report_reads, planned, program)?
                    }
                    None => None,
                };
                let misses = planned.map(|planned| planned.misses(ordinal, lookups.len()));
                // The subject whose Seq-align id asks for the blob (`blob_name`).
                let (text, opening) = run.load_text(
                    db,
                    &taxdb,
                    (blob, by),
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
                        cut_refusal: cut_refusal.clone(),
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
        // The first row that looks a taxid up opens the database after its read (a later call
        // of the run finds it open: `open_for_rows` makes the texts once).
        if !row_opened && rows.asks(db, ordinal, parsed, program)? {
            row_opened = true;
            run.open_for_rows(db, program)?;
        }
    }
    Ok(reads)
}

/// The checks of a load of [`load_read_texts`] of the sequence `blob` for the subject `ordinal`
/// whose lookups take the slot of the offset error's taxid `last`, on a database whose volume
/// holds that taxid and another of its slot ([`TaxdbRun::cache_collision`]):
/// - a report's read (`reads`) where NCBI may have dropped a sequence that the run read before
///   ([`reads_may_load_again`]): NCBI may load it again at a read that LOSAT does not model, which
///   looks the slot up again, so LOSAT's refusal;
/// - a load of a batch's prefetch whose misses NCBI's order does not give (`planned` without it:
///   a load before a stop whose place LOSAT does not know, [`prefetch_stop`]): LOSAT's refusal.
///
/// Else the refusal that the load's read carries where NCBI does not make it after all (a BLAST
/// XML report that the stream cuts, `report::query_warnings::ReadLoad::left_out`): the slot that
/// LOSAT followed through it is then not NCBI's.
fn slot_load_checks(
    db: &SeqDb,
    ordinal: u32,
    blob: u32,
    last: i32,
    reads: bool,
    planned: Option<&PrefetchSlots>,
    program: &str,
) -> Result<Option<std::sync::Arc<str>>> {
    let volume = db
        .volume_index(blob)
        .map_err(|error| database_error(error, program))?;
    let fetch = |what: &str| {
        slot_refusal(
            db,
            format!(
                "taxid {last}, the last entry of the taxonomy database, begins past the end of \
                 taxdb.btd, and {what} database sequence {} of the volume {}, which holds it and \
                 another taxid of its slot of NCBI BLAST+'s 256-slot taxonomy cache of a volume \
                 (CSeqDBVol::m_TaxCache, CSeqDBIntCache)",
                db.oid(blob),
                volume_name(db, volume)
            ),
        )
    };
    let message = crate::seqdb::taxdb::MSG_OFFSET_ERROR;
    if reads && reads_may_load_again(db) {
        return Err(taxdb_error(
            fetch(&format!(
                "a report reads, of a database of more than {} sequences or with an environment \
                 variable of NCBI BLAST+'s object manager set (OBJMGR_*): NCBI BLAST+ may release \
                 a sequence that the run read and load it again at a later read, which looks the \
                 slot up again and can write \"{message}\" again, as the reads come; LOSAT does \
                 not reproduce that",
                2 * RELEASED_KEPT
            )),
            program,
        ));
    }
    if planned.is_some_and(|planned| !planned.misses.contains_key(&ordinal)) {
        return Err(taxdb_error(
            fetch(&format!(
                "the prefetch of a query batch that stops at a lookup that fails loads, at a \
                 place among its loads that LOSAT does not know (before the stop or not): that \
                 place decides how often NCBI writes \"{message}\"; LOSAT does not reproduce that"
            )),
            program,
        ));
    }
    let refusal = taxdb_error(
        fetch(&format!(
            "a report that NCBI's XML stream cut reads, which LOSAT's model of the loads had \
             made: NCBI does not make that load, so the later loads look the slot up from another \
             state and may write \"{message}\" more or less often; LOSAT does not reproduce that"
        )),
        program,
    );
    Ok(Some(std::sync::Arc::from(format!("{refusal:#}"))))
}

/// LOSAT's refusal of the load of the sequence `holder` that the stopping subject `stop` (its
/// parsed Seq-ids `by`) is given instead ([`SubjectLoad::Instead`]) where its lookups take the slot
/// of the offset error's taxid `last`, on a database whose volume holds that taxid and another of
/// its slot ([`TaxdbRun::cache_collision`]): NCBI makes that load or not, at a place among the
/// prefetch's loads that LOSAT does not know, and as the blob of the stopping id's handle, which
/// looks the slot up again after the sequence's own load (`blob_name`).
fn instead_slot_check(
    db: &SeqDb,
    holder: u32,
    stop: u32,
    by: &crate::report::db_subject::ParsedSubject,
    last: i32,
    program: &str,
) -> Result<()> {
    let lookups = blob_taxids(db, stop, Some(by), holder, program)?;
    if !lookups
        .iter()
        .any(|&taxid| tax_cache_slot(taxid) == tax_cache_slot(last))
    {
        return Ok(());
    }
    let volume = db
        .volume_index(holder)
        .map_err(|error| database_error(error, program))?;
    Err(taxdb_error(
        slot_refusal(
            db,
            format!(
                "taxid {last}, the last entry of the taxonomy database, begins past the end of \
                 taxdb.btd, and the prefetch of a query batch that stops at a lookup that fails \
                 may load, for the Seq-align id {} of database sequence {} whose lookup NCBI \
                 BLAST+ may not make, database sequence {} of the volume {}, which holds it and \
                 another taxid of its slot of NCBI BLAST+'s 256-slot taxonomy cache of a volume \
                 (CSeqDBVol::m_TaxCache, CSeqDBIntCache): whether and where NCBI makes that load \
                 decides how often it writes \"{}\"; LOSAT does not reproduce that",
                String::from_utf8_lossy(&by.alignment_fasta()),
                db.oid(stop),
                db.oid(holder),
                volume_name(db, volume),
                crate::seqdb::taxdb::MSG_OFFSET_ERROR
            ),
        ),
        program,
    ))
}

/// The reader's refusal of `taxdb.bti` of the run's search path with `message` (as
/// [`slot_uncertain_error`]'s).
fn slot_refusal(db: &SeqDb, message: String) -> SeqDbError {
    SeqDbError::Unsupported {
        path: taxdb_bti_path(db),
        message,
    }
}

/// The path of `taxdb.bti` on the search path of `db` (`taxdb.bti` where none resolves), which
/// LOSAT's refusals of the taxonomy cache name.
fn taxdb_bti_path(db: &SeqDb) -> std::path::PathBuf {
    db.search_path()
        .and_then(|search_path| search_path.resolve_file(b"taxdb.bti"))
        .map_or_else(
            || std::path::PathBuf::from("taxdb.bti"),
            crate::utils::os_bytes::path_from_bytes,
        )
}

/// The name of the volume `volume` of `db` (empty where it has none), as LOSAT's refusals of the
/// taxonomy cache name it.
fn volume_name(db: &SeqDb, volume: usize) -> String {
    db.volumes()
        .get(volume)
        .map_or_else(String::new, |vol| vol.base().display().to_string())
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
/// one reads the sequence `blob` by the Seq-align id of the subject `by`, the lookups of its load
/// `lookups` (the positive taxids in `x_GetTaxonomy`'s order, [`blob_taxids`];
/// `ListedTitles::pssm_load_texts`),
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
    (blob, by): (u32, u32),
    lookups: &[i32],
    program: &str,
) -> Result<Vec<u8>> {
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1785
    // The scope keeps the id that the PSSM's read resolves until its reset
    // ([`TaxdbRun::note_scope_read`]).
    run.note_scope_read();
    let lead = if db.scope_utf8_loads() {
        scope_utf8_lead(db, parsed, &run.prefix(program))
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
        (blob, by),
        lookups,
        &run.prefix(program),
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
    /// The query that the stop's [`ErrorAfterBatchReports`] names (the batch's first query of the
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
    /// What the batch's prefetch gives for its Seq-align id in the run's one walk
    /// (`report::db_subject::HitLoads::prefetched`, `RunLoads::prefetched`), which the caller
    /// sets where it is flipped (`None`: the subject's own blob, as for a subject without parsed
    /// Seq-ids).
    pub(crate) load: Option<crate::report::db_subject::HitLoad>,
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
            load: None,
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
    /// `Err`: LOSAT's refusal of a load that such a report left out and whose lookups took the
    /// slot of the offset error's taxid on a database whose volume holds that taxid and another of
    /// its slot (`report::query_warnings::ReadLoad::cut_refusal`).
    fn left_texts(
        &self,
        query: usize,
    ) -> Result<Option<std::sync::Arc<crate::report::query_warnings::LoadTexts>>> {
        let left = crate::report::query_warnings::LoadTexts::left_before(
            self.reads,
            &|before| self.cuts.cut(before),
            self.prefetched,
            query,
        )
        .map_err(|refusal| {
            anyhow::Error::new(TaxonomyRefusal {
                text: refusal.to_string(),
                load_throws: false,
            })
        })?;
        Ok((!left.is_empty()).then_some(left))
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
/// ([`prefetch_gets_no_bioseq`]). A failed lookup ends the loads in NCBI's order of the batch's
/// Seq-id handles (of a map that counts every Seq-align id of the batch, [`batch_map`]): the loads
/// before it write their texts, those after it are not made ([`made_before_failure`]), and the
/// blobs of those made go to NCBI's data source's cache, not to the scope
/// ([`VolumeLookups::release_failed`]); a load whose place to it LOSAT does not know (a handle of
/// its key, a Seq-align id of two subjects, a map that LOSAT does not know) is made or not,
/// refused where it writes a text, where it reads a header with a UTF8String before a load that
/// NCBI surely makes does ([`unknown_utf8_load`]), and where it changes what a later lookup of the
/// slot of the offset error's taxid writes on a database whose volume holds that taxid and
/// another of its slot ([`prefetch_slot_plan`]). Such a batch asks no flipped Seq-align for its
/// length (NCBI's lookup throws before). The caller writes the text
/// where the batch's prefetch messages go, after the reader's messages of the batch and before the
/// warnings of its first query.
///
/// A subject whose lookup fails on a file of the database with NCBI's own exception
/// (`ParsedSubject::prefetch_stops`: a list of version 5 databases with and without the accession
/// table) stops the run there: `Err` with an [`ErrorAfterBatchReports`] of the query that `queries` names
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
        run.keep_loaded(None, before);
        return Err(prefetch_stop(db, run, &subjects, program, queries));
    }
    // NCBI's map of the batch counts every Seq-align id of it, those that a query resolves and
    // those whose lookup fails included.
    let map = batch_map(db, &subjects);
    let mut failed_lookup: Option<(u32, String)> = None;
    let mut failing = Vec::new();
    let mut fails = false;
    let mut uncertain: Option<(u32, String)> = None;
    let mut asked: std::collections::HashSet<u32> = Default::default();
    let mut loaded = Vec::new();
    let mut flips = Vec::new();
    let mut query_scope: Option<crate::blastinput::query_batch::QueryScope> = None;
    for prefetched @ PrefetchSubject {
        ordinal,
        parsed,
        flipped,
        load,
    } in subjects
    {
        if let Some(subject) = parsed {
            // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3288-3297
            // A query Bioseq in the batch's scope resolves the Seq-align id before the BLAST
            // database data loader is called, so a subject whose loader would throw or return an
            // empty handle on its normalized PDB handle is resolved to the query Bioseq and never
            // fails `PrefetchSeqData` (line RG-KB item 12).
            if subject.load_throws() || subject.scope_finds_none() {
                let scope = query_scope.get_or_insert_with(|| (queries.scope)());
                if scope.has_key(&subject.query_scope_keys().own) {
                    continue;
                }
            }
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
            if flipped {
                flips.push((ordinal, subject, load));
            }
        }
        loaded.push((ordinal, parsed));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1110-1141
    // `LoadSequencesToScope` asks the handle of a flipped Seq-align for its length only after
    // `GetBioseqHandles` returned, which a failed lookup ends with its exception: a batch with
    // one asks no length. What NCBI's places give at the prefetch for the id is the run's one
    // walk's (`PrefetchSubject::load`).
    if failed_lookup.is_none() {
        for (ordinal, subject, load) in flips {
            if fails {
                break;
            }
            if asked.insert(ordinal) {
                match prefetch_gets_no_bioseq(db, ordinal, subject, load, program)? {
                    Ok(no_bioseq) => fails = no_bioseq,
                    Err(reason) => {
                        uncertain.get_or_insert((ordinal, reason));
                    }
                }
            }
        }
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
        || scope_utf8_before_failure(db, &loaded, &failing, map, program)?;
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (cited at
    // `loads_before_failure`)
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
    // The failed lookup ends the loads in the order of the batch's Seq-id handles: the loads
    // before it are made, those after it are not (their subjects are loaded where a report reads
    // them, one by one), and those whose place LOSAT does not know are made in the model after
    // the others and refused where they matter ([`made_before_failure`]).
    let unknown = match &failed_lookup {
        Some(stop) => {
            let whys;
            (loaded, whys) = made_before_failure(db, loaded, &failing, map);
            unknown_utf8_load(db, stop, &loaded, &whys, map, program)?;
            whys
        }
        None => Vec::new(),
    };
    let failure = failed_lookup.as_ref().map(|stop| FailedPrefetch {
        stop,
        unknown: &unknown,
    });
    let first_utf8 = loaded
        .iter()
        .find(|(_, parsed)| parsed.is_some_and(|subject| subject.scope_utf8_read().is_some()))
        .map(|&(ordinal, _)| ordinal);
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
    // A load that throws where a file of the taxonomy database does not open ends the prefetch's
    // loads, which catches it as a `CException` with the warning (line WQA D18b probe 12). The
    // loads write no text then: the taxonomy database never opens.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1131-1153
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:1085-1097
    // A prefetch pins its loaded TSEs via `SetLinkUsedTSE` only when it reaches `AddBioseq`
    // without throwing on a failed lookup or a flipped no-bioseq alignment; within `AddBioseq`,
    // `CanResolveRange` pins segments in `CSeq_id::CompareOrdered` order up to the first empty
    // `CBioseq_Handle` (computed per load in `prefetch_slot_plan`).
    let pins = !fails && failed_lookup.is_none();
    let planned = prefetch_slot_plan(db, run, &loaded, failure.as_ref(), pins, map, program);
    let mut reads = match planned.and_then(|planned| {
        let order = LoadsOrder::Prefetch(planned.as_ref());
        if posted {
            prefetch_loads(db, run, &loaded, true, order, map, program)
        } else {
            db.without_scope_utf8_loads(|| {
                prefetch_loads(db, run, &loaded, true, order, map, program)
            })
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
    if let Some(failure) = &failure {
        // NCBI's warning of a UTF8String header before the failure is not a text of the
        // taxonomy database (`scope_utf8_before_failure`).
        if let Some((read, why)) = reads
            .iter()
            .find_map(|read| Some((read, failure.unknown_text(read)?)))
        {
            return Err(failed_lookup_refusal(
                db,
                failure.stop,
                (loaded[read.subject as usize].0, why),
                "makes the taxonomy database write a text",
                program,
            ));
        }
    }
    // The reads' keys are the loads' indexes (`prefetch_loads`), which nothing reads after
    // these checks.
    for read in &mut reads {
        read.subject = 0;
    }
    let mut posts = crate::report::diag_log::Posts::new();
    posts.loads(reads);
    if failed_lookup.is_some() || fails {
        posts.lines(&prefetch_warning(&run.prefix(program)));
    }
    Ok(posts)
}

/// The subjects `subjects` of a query batch's prefetch (BLASTN's and TBLASTX's, whose
/// Seq-aligns can flip, [`PrefetchSubject::flipped`]) with what the prefetch gives for each
/// Seq-align id in the run's one walk (`run`, `RunLoads::prefetched` after the run's record on
/// `taxdb`, [`TaxdbRun::loads_before`]), for [`batch_prefetch`].
pub(crate) fn walk_prefetch<'s>(
    taxdb: &TaxdbRun,
    run: &crate::report::db_subject::RunLoads,
    subjects: &[PrefetchSubject<'s>],
) -> Vec<PrefetchSubject<'s>> {
    // Each sequence once (TBLASTX gives one for each HSP).
    let set: std::collections::BTreeMap<usize, &crate::report::db_subject::ParsedSubject> =
        subjects
            .iter()
            .filter_map(|subject| Some((subject.ordinal as usize, subject.parsed?)))
            .collect();
    let answers = run.prefetched(taxdb.loads_before(None), set);
    subjects
        .iter()
        .map(|&subject| PrefetchSubject {
            load: answers.get(&(subject.ordinal as usize)).copied(),
            ..subject
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (cited at
// `loads_before_failure`)
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
// The failed lookup ends the prefetch's loads in the order of the batch's Seq-id handles, and the
// prefetch writes its warning: a load whose place to it LOSAT does not know is made or not.
/// LOSAT's refusal of a batch's prefetch with a Seq-align id whose lookup fails (`stop`: the
/// subject's ordinal id and the id, [`batch_prefetch`]) where a load whose place to the failure
/// LOSAT does not know matters (`load`: the loaded subject's ordinal id and why its place is not
/// known, [`made_before_failure`]; `what`: what that load does that NCBI does only where it makes
/// the load before the failure).
fn failed_lookup_refusal(
    db: &SeqDb,
    stop: &(u32, String),
    (load, why): (u32, &str),
    what: &str,
    program: &str,
) -> anyhow::Error {
    let (ordinal, id) = stop;
    parsed_refusal(
        db.oid(*ordinal),
        &format!(
            "has the Seq-align id {id}, whose lookup in the database fails where NCBI BLAST+ 2.17.0's prefetch loads the batch's subjects, in an order that LOSAT does not reproduce: database sequence {}, which the prefetch may load before or not at all ({why}), {what}",
            db.oid(load)
        ),
        program,
    )
}

// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769 (cited at
// `loads_before_failure`)
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
// The lookup that throws ends `CDataLoader::GetBlobs`, and the prefetch catches the exception
// with its warning; the loads before it have written their texts.
/// A batch's prefetch with a Seq-align id whose lookup fails ([`batch_prefetch`]): the id
/// (`stop`: the failing subject's ordinal id and the id) and, for each load that LOSAT's model
/// makes (in the order of [`made_before_failure`]'s loads), why NCBI may not make it before the
/// failure (`None`: it surely does).
struct FailedPrefetch<'a> {
    stop: &'a (u32, String),
    unknown: &'a [Option<&'static str>],
}

// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769 (cited at
// `loads_before_failure`)
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
impl FailedPrefetch<'_> {
    // NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769 (cited at
    // `loads_before_failure`)
    // A load that NCBI may not make before the failed lookup writes its texts only where it is
    // made: a text of the taxonomy database (opening, an offset error) then depends on its place.
    /// Why NCBI may not make the load of the prefetch's read `read` (keyed by its load's index,
    /// [`prefetch_loads`]) whose text holds more than NCBI's warning of a UTF8String header, a text
    /// of the taxonomy database (`None`: the load writes none, or NCBI surely makes it).
    fn unknown_text(
        &self,
        read: &crate::report::query_warnings::SubjectRead,
    ) -> Option<&'static str> {
        let lead = read.load.as_ref().map_or(0, |load| load.lead);
        self.unknown[read.subject as usize].filter(|_| read.text.len() > lead)
    }
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

/// The Seq-align id `gnl|BL_ORD_ID|<OID>` of a database sequence without parsed Seq-ids.
fn ordinal_seqid(db: &SeqDb, ordinal: u32) -> crate::seqid::ParsedSeqId {
    crate::seqid::ParsedSeqId::General(crate::seqid::model::Dbtag {
        db: b"BL_ORD_ID".to_vec(),
        tag: crate::seqid::model::ObjectId::Id(db.oid(ordinal).get() as i32),
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:667-697
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_aux.hpp:186-200
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1113-1153
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3346-3369
// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:1036-1097
// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_handle.cpp:941-944
/// Computes the empty `CBioseq_Handle` Seq-align ids that stop `CSeqMap::CanResolveRange` from
/// calling `AddUsedTSE` on later segments of a prefetch's delta Bioseq. `s_ExtractSeqidsAndRanges`
/// orders the batch's ids by `CConstRefCSeqId_LessThan` (`CSeq_id::CompareOrdered`),
/// `CScope_Impl::GetBioseqHandles` restores that order (`RestoreOrder`) even beyond 200 ids, and
/// `CanResolveRange` walks the delta segments in that order until the first empty handle.
/// Returns `Ok(None)` when the batch cannot be assumed to pin anything (`!pins`, a `GetBioseq`
/// throw before `AddBioseq`, or `ParsedSubject::holders_unknown`), and `Ok(Some(empty_ids))`
/// otherwise (a load pins iff its id strictly precedes every id in `empty_ids` under
/// `compare_ordered`).
fn prefetch_pin_stops(
    db: &SeqDb,
    loaded: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    pins: bool,
    program: &str,
) -> Result<Option<Vec<crate::seqid::ParsedSeqId>>> {
    if !pins {
        return Ok(None);
    }
    let mut empty_ids = Vec::new();
    for &(ordinal, parsed) in loaded {
        if let Some(subject) = parsed {
            let bioseq = scope_bioseq(db, ordinal, subject, program)?;
            if bioseq.throws || subject.holders_unknown() {
                return Ok(None);
            }
            if !bioseq.bioseq {
                empty_ids.push(subject.alignment_id().clone());
            }
        }
    }
    Ok(Some(empty_ids))
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
/// batch whose lookup fails ([`batch_prefetch`]), where NCBI's loads end: `loaded` are then the
/// loads that LOSAT's model makes ([`made_before_failure`]); a load of the slot whose place to
/// the failure LOSAT does not know is made or not ([`SlotLoad::maybe`]), refused where that
/// changes the texts or leaves a slot that a later lookup depends on ([`failed_lookup_refusal`]),
/// and NCBI's data source holds what the loads leave as it releases them at the failure
/// ([`VolumeLookups::release_failed`]). `pins`: whether the batch's
/// prefetch reaches `AddBioseq` without throwing on a failed lookup or a flipped no-bioseq
/// alignment, so that loads strictly preceding every empty `CBioseq_Handle` in
/// `CSeq_id::CompareOrdered` order ([`prefetch_pin_stops`]) are marked pinned via `SetLinkUsedTSE`
/// (`SlotLoad::pins`). `map`: why LOSAT does not know NCBI's map of the batch ([`batch_map`]).
fn prefetch_slot_plan(
    db: &SeqDb,
    run: &TaxdbRun,
    loaded: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    stop: Option<&FailedPrefetch<'_>>,
    pins: bool,
    map: Option<&'static str>,
    program: &str,
) -> Result<Option<PrefetchSlots>> {
    if !run.found(db) || db.search_path().is_none() || run.cache_collision(db).is_none() {
        return Ok(None);
    }
    let Some(last) = run.prefetch_texts(db).offset_error_taxid else {
        return Ok(None);
    };
    let pin_stops = prefetch_pin_stops(db, loaded, pins, program)?;
    // NCBI's map of the batch's Seq-id handles ([`load_places`]): each load's key, and whether its
    // place among the loads of its key is not known (a batch of more than `PREFETCH_ONE_MAP_IDS`
    // distinct ids, or a Seq-align id of two subjects), so that LOSAT takes all of its volume's
    // loads as one key.
    let places = prefetch_places(db, loaded, map);
    let mut seen = std::collections::HashSet::new();
    let mut loads = Vec::new();
    let mut unordered = std::collections::HashSet::new();
    // Every load that the model makes, for what NCBI's data source holds after a failure, and
    // each load of the slot that NCBI may not make (its blob, its subject and why).
    let mut made = Vec::new();
    let mut whys = Vec::new();
    for (index, (&(ordinal, parsed), place)) in loaded.iter().zip(&places).enumerate() {
        if !seen.insert(ordinal) {
            continue;
        }
        let maybe = stop.and_then(|failure| failure.unknown[index]);
        let blob = match parsed {
            Some(subject) => match prefetch_reads(db, ordinal, subject, true, false, program)? {
                Ok(blob) => blob,
                Err(reason) => return Err(parsed_refusal(db.oid(ordinal), &reason, program)),
            },
            None => Some(ordinal),
        };
        made.push(FailedLoad {
            ordinal,
            order: place.many.is_none().then_some(place.order),
            surely: maybe.is_none() && blob.is_some() && !run.unrecorded_calls(),
        });
        let Some(blob) = blob else {
            continue;
        };
        let lookups = blob_taxids(db, ordinal, parsed, blob, program)?;
        if !lookups
            .iter()
            .any(|&taxid| tax_cache_slot(taxid) == tax_cache_slot(last))
        {
            continue;
        }
        if let Some(why) = maybe {
            whys.push((blob, ordinal, why));
        }
        let volume = db
            .volume_index(blob)
            .map_err(|error| database_error(error, program))?;
        if place.unknown.is_some() {
            unordered.insert(volume);
        }
        // A load pins iff `CanResolveRange` reaches it before any empty handle in
        // `CSeq_id::CompareOrdered` order (uncertain patent ties return `None` -> unpinned).
        let load_pins = pin_stops.as_ref().is_some_and(|empty| {
            if empty.is_empty() {
                return true;
            }
            let ordinal_id;
            let id = match parsed {
                Some(subject) => subject.alignment_id(),
                None => {
                    ordinal_id = ordinal_seqid(db, ordinal);
                    &ordinal_id
                }
            };
            empty.iter().all(|other| {
                crate::seqid::defline_set::compare_ordered(id, other)
                    == Some(std::cmp::Ordering::Less)
            })
        });
        loads.push(SlotLoad {
            ordinal,
            blob,
            volume,
            lookups,
            order: place.order,
            name: blob_name(db, blob, ordinal, program)?,
            pins: load_pins,
            maybe: maybe.is_some(),
        });
    }
    // A prefetch that does not fail holds what it loads to the end of the batch (cited above
    // `RELEASED_KEPT`): the batch's reports read the sequences again without a load. One whose
    // lookup failed leaves them to the data source's cache.
    let misses = run
        .volume_lookups(db)
        .prefetch_slots(
            &loads,
            |volume| unordered.contains(&volume),
            last,
            stop.map(|_| &made[..]),
        )
        .map_err(|(volume, cause)| {
            // Whether NCBI makes a load that it may make before the failure or not at all decides.
            let maybe = whys
                .iter()
                .find(|&&(blob, _, _)| blob == cause.ordinal)
                .filter(|_| cause.again && cause.order.is_none());
            match stop.zip(maybe) {
                Some((failure, &(_, ordinal, why))) => failed_lookup_refusal(
                    db,
                    failure.stop,
                    (ordinal, why),
                    &format!("looks up a taxid of the slot of taxid {last} in NCBI BLAST+'s 256-slot taxonomy cache of a volume that holds that taxid and another of its slot, which decides what the prefetch writes"),
                    program,
                ),
                None => taxdb_error(slot_uncertain_error(db, volume, last, cause), program),
            }
        })?;
    Ok(Some(PrefetchSlots { misses }))
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
// The prefetch's loads go in the order of the map of their Seq-id handles; the lookup that throws
// on its id (`GetRecords` -> `x_GetOid`) ends them, and the prefetch catches it with its warning
// (blast_app_util.cpp:778-788): a load after it is not made, so its lookup reads no header
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
    map: Option<&'static str>,
    program: &str,
) -> Result<bool> {
    let subjects: Vec<PrefetchSubject<'_>> = loaded
        .iter()
        .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed))
        .chain(failing.iter().copied())
        .collect();
    let places = load_places(db, &subjects, map);
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

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (`x_GetBioseqHandlesSorted`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
// The prefetch asks the data loader for the blobs of the batch's Seq-align ids in the order of its
// map of their Seq-id handles; the first lookup in that order that throws on its id ends the asks,
// and the prefetch catches the exception with its warning. So the loads of the handles before
// that lookup are made (their texts, their steps through the taxonomy cache's slot), and those
// after it are not, here (a report that reads such a subject later loads it then). Probe
// `O4b2/p10/probe` (NCBI 2.17.0, six programs, four heap layouts each, BLASTDB with the offset
// error of taxid 99990): `gb|FFF00092.1|` (packed 92) before `gnl|BL_ORD_ID|0` (a text tag, not
// packed) writes the offset error, then the warning; `gb|FFF00094.1|` (94) after
// `gnl|BL_ORD_ID|x002` (packed 2) writes the warning alone; of `gb|FFF00010.1|` (10),
// `gnl|BL_ORD_ID|x050` (50, fails) and `gb|FFF00090.1|` (90, taxid 150 of the same slot), one offset
// error, then the warning.
/// [`LoadBeforeStop`] of each of the loads `loaded` (each `(ordinal, parsed Seq-ids)`, in any order)
/// of a batch's prefetch whose lookups that throw on their ids are `failing` ([`batch_prefetch`]):
/// `Yes` for a load whose Seq-id handle comes before the first failing one's in NCBI's map
/// ([`load_places`]), `No` for one after it, `Unknown` for a handle of the same key (the same packed
/// number, or neither packed: [`ADDRESS_ORDER`]; one handle is one load) and for every load whose
/// place is not known ([`LoadPlace::unknown`]: a Seq-align id that another subject has too, a batch
/// whose map LOSAT does not know).
fn loads_before_failure(
    db: &SeqDb,
    loaded: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    failing: &[PrefetchSubject<'_>],
    map: Option<&'static str>,
) -> Vec<LoadBeforeStop> {
    let subjects: Vec<PrefetchSubject<'_>> = loaded
        .iter()
        .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed))
        .chain(failing.iter().copied())
        .collect();
    let places = load_places(db, &subjects, map);
    let (loads, fails) = places.split_at(loaded.len());
    let failure = *fails
        .iter()
        .min_by_key(|place| place.order)
        .expect("a failing lookup");
    loads
        .iter()
        .map(
            |place| match place.unknown.or_else(|| place.order_unknown(failure)) {
                Some(why) => LoadBeforeStop::Unknown(why),
                None if place.order < failure.order => LoadBeforeStop::Yes,
                None => LoadBeforeStop::No,
            },
        )
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 (cited at
// `loads_before_failure`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_loader.cpp:764-769 (cited at
// `loads_before_failure`)
// The loads before the failed lookup in the order of the map are made, those after it are not;
// one of the failure's key goes before it or after it by the address of its handle's info.
/// The loads of a batch's prefetch with a Seq-align id whose lookup fails ([`batch_prefetch`])
/// that LOSAT's model makes, of the loads `loaded` (each `(ordinal, parsed Seq-ids)`, in the
/// batch's order) and the failing lookups `failing` ([`loads_before_failure`]): those that NCBI
/// surely makes before the failure (`Yes`), in the batch's order, then those whose place to it
/// LOSAT does not know (`Unknown`), which NCBI makes, if at all, after every `Yes` one, each with
/// why (`None` for a `Yes` one). The loads after the failure (`No`) are not made: they write no
/// text there and take no step through the taxonomy cache's slot. One rule for every program
/// (they share [`batch_prefetch`]).
fn made_before_failure<'s>(
    db: &SeqDb,
    loaded: Vec<(u32, Option<&'s crate::report::db_subject::ParsedSubject>)>,
    failing: &[PrefetchSubject<'_>],
    map: Option<&'static str>,
) -> (
    Vec<(u32, Option<&'s crate::report::db_subject::ParsedSubject>)>,
    Vec<Option<&'static str>>,
) {
    let before = loads_before_failure(db, &loaded, failing, map);
    let mut made: Vec<_> = loaded
        .iter()
        .zip(&before)
        .filter(|(_, before)| **before == LoadBeforeStop::Yes)
        .map(|(&load, _)| (load, None))
        .collect();
    made.extend(
        loaded
            .iter()
            .zip(&before)
            .filter_map(|(&load, before)| match before {
                LoadBeforeStop::Unknown(why) => Some((load, Some(*why))),
                _ => None,
            }),
    );
    made.into_iter().unzip()
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:299-310
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:621-622
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393-406 (cited
// at `report::db_subject::ParsedSubject::scope_utf8`)
// The reader posts the warning with `ERR_POST_X_ONCE`: the first read of such a header in the
// process posts it, and no later one. So where a load that NCBI surely makes before the failure
// reads such a header and comes first in the map's order (a smaller key), the prefetch posts it
// there, whether a later load is made or not.
/// LOSAT's refusal of a batch's prefetch with a Seq-align id whose lookup fails (`stop`) where a
/// load whose place to the failure LOSAT does not know (`unknown`, by the loads `loaded` that the
/// model makes, [`made_before_failure`]) reads by its lookup a header with a UTF8String where a
/// VisibleString is expected (`ParsedSubject::scope_utf8`), in a run that posts NCBI's warning
/// (803.10) of such a read (`SeqDb::scope_utf8_loads`): NCBI posts it at that load where it makes
/// it first, so where the prefetch writes it depends on that place. Not where a load that NCBI
/// surely makes reads such a header with a smaller key of the batch's map (which LOSAT knows,
/// `map`, [`batch_map`]): that load posts it, and no later one does.
#[cold]
#[inline(never)]
fn unknown_utf8_load(
    db: &SeqDb,
    stop: &(u32, String),
    loaded: &[(u32, Option<&crate::report::db_subject::ParsedSubject>)],
    unknown: &[Option<&'static str>],
    map: Option<&'static str>,
    program: &str,
) -> Result<()> {
    if !db.scope_utf8_loads() {
        return Ok(());
    }
    let utf8 = |index: usize| {
        loaded[index]
            .1
            .is_some_and(|subject| subject.scope_utf8_read().is_some())
    };
    let places = prefetch_places(db, loaded, map);
    // The smallest key of a load that NCBI surely makes and whose lookup reads such a header.
    let first = (0..loaded.len())
        .filter(|&index| unknown[index].is_none() && utf8(index))
        .map(|index| places[index].order)
        .min();
    let posted_before = |index: usize| {
        places[index].many.is_none() && first.is_some_and(|first| first < places[index].order)
    };
    match (0..loaded.len()).find_map(|index| {
        let why = unknown[index]?;
        (utf8(index) && !posted_before(index)).then_some((loaded[index].0, why))
    }) {
        Some(load) => Err(failed_lookup_refusal(
            db,
            stop,
            load,
            "reads by its lookup a header with a UTF8String where a VisibleString is expected, whose warning (803.10) NCBI writes at that load",
            program,
        )),
        None => Ok(()),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:778-786
/// NCBI's warning where the prefetch of a query batch fails ([`batch_prefetch`], and BLASTX's and
/// RPSTBLASTN's Seq-align of a subject row whose range it cannot take,
/// `blastx::empty_rows::prefetch_reads_unranged_row`): the one producer of its text.
pub(crate) fn prefetch_warning(app: &str) -> Vec<u8> {
    // A post of NCBI's main thread (`ncbi_exception::post_line`: the header alone after the run's diagnostic
    // stream failed).
    crate::ncbi_exception::post_line(
        crate::ncbi_exception::Severity::Warning,
        None,
        Some(&app.to_ascii_lowercase()),
        b"Error pre-fetching sequence data",
    )
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
/// other subjects before the first lookup that stops it, or that of a batch whose lookup of a
/// Seq-align id throws on the id, before the first such lookup ([`loads_before_failure`]).
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
/// [`PREFETCH_ONE_MAP_IDS`] distinct Seq-align ids or with one whose handle LOSAT does not know
/// (whose stops' order LOSAT does not know either).
/// The one owner of NCBI's order of a stopping batch's loads and lookups.
fn load_before_stop(db: &SeqDb, subjects: &[PrefetchSubject<'_>]) -> LoadOrder {
    let map = batch_map(db, subjects);
    let places = load_places(db, subjects, map);
    let many = map.is_some();
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
/// another subject has too, and every subject of a batch whose map LOSAT does not know (`many`,
/// [`batch_map`] of the whole batch). The one owner of that map for the stops
/// ([`load_before_stop`]), for the loads before a lookup that throws on its id
/// ([`loads_before_failure`]) and for the place of NCBI's warning of a UTF8String header among
/// the loads ([`scope_utf8_prefetch_loads`]).
fn load_places(
    db: &SeqDb,
    subjects: &[PrefetchSubject<'_>],
    many: Option<&'static str>,
) -> Vec<LoadPlace> {
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
    for (subject, (_, handle)) in subjects.iter().zip(&keys) {
        if !ordinals.insert(subject.ordinal) {
            continue;
        }
        let Some(handle) = handle else {
            continue;
        };
        match handles.entry(handle) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(subject.ordinal);
            }
            std::collections::hash_map::Entry::Occupied(entry) => {
                shared.insert(*entry.get());
                shared.insert(subject.ordinal);
            }
        }
    }
    subjects
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
            unknown: many.or_else(|| {
                shared
                    .contains(&subject.ordinal)
                    .then_some("another subject of the batch has the same Seq-align id")
            }),
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3346-3366
// NCBI reference (598d8ae6): c++/src/objmgr/seq_id_sort.cpp:154-163
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:658-699
// `GetBioseqHandles` counts the distinct Seq-id handles of every id that the prefetch gives it
// (`CSortedSeq_ids`, whose `GetSortedIds` keeps one of equal handles), one for each range of the
// batch's Seq-aligns (`s_ExtractSeqidsAndRanges`), before it asks the scope or a data source for
// any: the ids that a query of the scope resolves (`x_GetBioseqHandlesSorted` takes them from the
// scope) and those whose lookup fails count too.
/// Why LOSAT does not know NCBI's map of the Seq-id handles of the batch whose reported sequences
/// are `subjects` (every one of the batch, in any order, a sequence possibly more than once):
/// [`map_unknown`] of their Seq-align ids, counted as NCBI counts them. `None` where it knows it.
fn batch_map(db: &SeqDb, subjects: &[PrefetchSubject<'_>]) -> Option<&'static str> {
    let keys: Vec<(u64, Option<Vec<u8>>)> = subjects
        .iter()
        .map(|subject| prefetch_load_key(db, subject))
        .collect();
    let ids: std::collections::HashSet<&[u8]> = keys
        .iter()
        .filter_map(|(_, handle)| handle.as_deref())
        .collect();
    map_unknown(&keys, ids.len())
}

/// Why LOSAT does not know NCBI's map of the Seq-id handles of a batch whose Seq-align ids have
/// the keys `keys` ([`prefetch_load_key`]), `ids` of them distinct handles that LOSAT knows: more
/// than [`PREFETCH_ONE_MAP_IDS`] ids, which NCBI's prefetch sorts otherwise and loads in pieces;
/// or an id without such a handle, as NCBI makes a new handle for each Seq-align of a Textseq-id
/// with neither accession nor name (one for each HSP of a gapped search) and LOSAT does not know a
/// giim id's, so the number of handles in the map is not known (review 1 of line NP item O4b2-3,
/// F6; re-review of item O4b2-7, R2-H1).
fn map_unknown(keys: &[(u64, Option<Vec<u8>>)], ids: usize) -> Option<&'static str> {
    if keys.iter().any(|(_, handle)| handle.is_none()) {
        Some(UNNAMED_MAP)
    } else {
        (ids > PREFETCH_ONE_MAP_IDS).then_some(MANY_IDS)
    }
}

/// Why NCBI's order of a batch's loads is not known where it has more than
/// [`PREFETCH_ONE_MAP_IDS`] Seq-align ids ([`load_places`]).
const MANY_IDS: &str = "the batch has more than 200 Seq-align ids, which NCBI's prefetch sorts by CSortableSeq_id and loads in pieces of 100";

/// Why NCBI's order of a batch's loads is not known where a Seq-align id has no Seq-id handle
/// that LOSAT knows ([`load_places`]).
const UNNAMED_MAP: &str = "a Seq-align id of the batch has no Seq-id handle that LOSAT knows (NCBI makes one for each Seq-align of a text id with neither accession nor name; a giim id), so the number of ids in NCBI's map, which decides whether its prefetch sorts them by CSortableSeq_id in pieces of 100, is not known";

/// NCBI's order of the loads and lookups of a stopping batch's subjects ([`load_before_stop`]).
struct LoadOrder {
    /// The place of each subject's load, in the order of the batch's subjects.
    places: Vec<LoadPlace>,
    /// The stopping subjects (their indexes in the batch's subjects, each sequence once) in the
    /// order of NCBI's lookups of their Seq-align ids: by their handles' keys
    /// ([`LoadOrder::stop_key`]), a group of one key in the batch's order (NCBI orders a group by
    /// the addresses of the handles' infos).
    stops: Vec<usize>,
    /// LOSAT does not know NCBI's map of the batch ([`LoadPlace::many`]).
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
    /// Why LOSAT knows the order of none of the batch's loads ([`MANY_IDS`], [`UNNAMED_MAP`]).
    many: Option<&'static str>,
}

impl LoadPlace {
    /// Why NCBI's order of this load and the load `other` is not known, for the place of NCBI's
    /// warning of a UTF8String header among a prefetch's loads ([`scope_utf8_prefetch_loads`]) and
    /// for the place of a load to a lookup that throws on its id ([`loads_before_failure`]):
    /// a batch whose map LOSAT does not know (`many`), or two handles of one key (one packed
    /// number, or neither packed), which NCBI's map orders by the address of their info
    /// ([`ADDRESS_ORDER`]); `None` where the keys order them or they are one handle (one load).
    fn order_unknown(self, other: LoadPlace) -> Option<&'static str> {
        if let Some(why) = self.many {
            return Some(why);
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
///   and, after calls whose loads nothing recorded ([`TaxdbRun::note_calls_before`]),
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
    let keys = parsed.query_scope_keys();
    let scope = (queries.scope)();
    if scope.has_key(&keys.own) || scope.has_key(&keys.weak) {
        return Ok(Some(StopUnsure {
            why: "NCBI's scope holds a query (of the batch, or in outfmt 5 of a batch before) with an id that it takes for this one, or that this one matches weakly, at the batch's prefetch, where it finds that query before it asks the database".to_owned(),
            instead: Some(Vec::new()),
        }));
    }
    // The first reason, and every sequence that may be given instead (by OID).
    let mut why: Option<String> = None;
    let mut instead: Option<Vec<u32>> = Some(Vec::new());
    let mut give = |oid: u32, reason: &dyn Fn() -> String| {
        why.get_or_insert_with(reason);
        if let Some(instead) = instead.as_mut() {
            instead.push(oid);
        }
    };
    // A load of the batch that may come before the stop registers the blob that its lookup
    // gives (`ParsedSubject::registered_by`; an id without a handle that a function can know
    // may be any of them).
    let unknown_handle = parsed.id_handle_unknown();
    for (subject, before) in subjects.iter().zip(before) {
        if *before == LoadBeforeStop::No {
            continue;
        }
        let given = match subject.parsed {
            Some(other) => parsed
                .registered_by(other)
                .map(|given| db.oid(given.ordinal()).get()),
            None => unknown_handle.then(|| db.oid(subject.ordinal).get()),
        };
        if let Some(oid) = given {
            give(oid, &|| {
                format!(
                "a load of the batch's prefetch that may come before it gives a sequence of the database (OID {oid}) that holds an id that NCBI's Seq-id handles take for it, which registers the id in NCBI's BLAST database data loader, which then gives that sequence for it without the lookup"
            )
            });
        }
    }
    let own = db.oid(ordinal).get();
    // The data loader's record of the run's loads before the batch (`ScopeRecord`): the
    // sequences that registered the id's handle there, the last one (each where a note of
    // several registered it).
    for sequence in record.loads_before(None).registrants_of(parsed) {
        let oid = db.oid(sequence.ordinal()).get();
        if oid != own {
            give(oid, &|| {
                format!("another sequence of the database (OID {oid}) holds an id that NCBI's Seq-id handles take for it, and the run loaded that sequence before the batch: its load there registers the id in NCBI's BLAST database data loader, which then gives that sequence for it without the lookup")
            });
        }
    }
    // After calls whose loads nothing recorded, every other holder whose own lookup does not
    // fail (`HeaderReader::unrecorded_holders`, found in the LMDB files that have the accession
    // table).
    if record.unrecorded_calls() {
        match parsed.unrecorded_loader() {
            Ok(holders) => {
                for oid in holders.iter().map(|holder| db.oid(holder.ordinal()).get()) {
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
/// a file with NCBI's exception ([`batch_prefetch`]), after the lookups and the loads that `run`
/// records ([`TaxdbRun::note_loaded`]): the [`ErrorAfterBatchReports`] of `queries` (its `error_query`) with
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
            run,
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
    // another taxid of the slot, whose loads follow NCBI's order below; a lookup that reads a
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
    // On a database whose volume holds the offset error's taxid and another of its slot, the loads
    // that surely come before the stop follow NCBI's order of a prefetch's loads
    // ([`prefetch_slot_plan`]); one whose place LOSAT does not know and that takes that slot is
    // refused (`slot_load_checks`).
    let mut reads = if loads.is_empty() {
        Vec::new()
    } else {
        let sure: Vec<(u32, Option<&crate::report::db_subject::ParsedSubject>)> = yes
            .iter()
            .map(|&index| (subjects[index].ordinal, subjects[index].parsed))
            .collect();
        let map = batch_map(db, subjects);
        let read =
            prefetch_slot_plan(db, run, &sure, None, false, map, program).and_then(|planned| {
                load_read_texts(
                    db,
                    run,
                    loads,
                    RowLookups::None,
                    LoadsOrder::Prefetch(planned.as_ref()),
                    program,
                )
            });
        match read {
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
    let left = match queries
        .xml_before
        .map(|before| before.left_texts(queries.error_query))
        .transpose()
    {
        Ok(left) => left.flatten(),
        Err(error) => return error,
    };
    if let Some(left) = left {
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
    anyhow::Error::new(ErrorAfterBatchReports {
        query: Some(queries.error_query),
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
        title_reads: Vec::new(),
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
    seen_by_scope: bool,
    program: &str,
) -> Result<std::result::Result<Option<u32>, String>> {
    use crate::report::db_subject::Prefetched;
    let (other, reason) = match subject.prefetched_with_scope(seen_by_scope) {
        Prefetched::Itself => return Ok(Ok(Some(ordinal))),
        Prefetched::Nothing => return Ok(Ok(None)),
        Prefetched::Failed(reason) => return Ok(Err(reason)),
        Prefetched::Another {
            other: None,
            reason,
        } => return Ok(if texts { Err(reason) } else { Ok(None) }),
        Prefetched::Another {
            other: Some(other),
            reason,
        } => (other, reason),
    };
    let other = other.ordinal();
    let set = db
        .defline_set(other)
        .map_err(|error| database_error(error, program))?;
    // `x_GetOid` keeps the OID only where one of its ids compares equal to the Seq-align id
    // (`scope_bioseq`), else the prefetch reads nothing.
    let read = other_bioseq(&set, other, subject);
    if read.is_some_and(|read| read.throws) {
        return Ok(Err(format!(
            "{reason}; NCBI BLAST+ 2.17.0 can read that sequence for it, whose definition lines do not hold the Seq-align id, and then stops (oid headers do not contain target gi/seq_id)"
        )));
    }
    if !texts {
        return Ok(Ok(None));
    }
    Ok(Ok(read.map(|_| other)))
}

/// LOSAT's refusal of the reported sequence `ordinal` with parsed Seq-ids `subject` where the
/// batches' prefetch reads it though no report reads its Bioseq (`ParsedShown::prefetched`), with
/// or without `taxdb.bti` on the search path: another sequence whose definition lines do not
/// hold its Seq-align id, which NCBI reads for it and stops ([`prefetch_reads`] without
/// texts; the texts are `prefetch_messages`'s), and where NCBI's places give another sequence for
/// its Seq-align id at the prefetch (`load`, `HitLoad::Shown`: `report::db_subject::ScopeWalk`)
/// with `taxdb.bti` on the search path, a sequence whose taxonomy texts can differ from the
/// subject's (`lookup_texts_differ`; any where LOSAT does not name it). A lookup that fails is not
/// refused: the batch's prefetch writes its warning ([`batch_prefetch`]).
pub(crate) fn check_prefetched_parsed(
    db: &SeqDb,
    subject: &crate::report::db_subject::ParsedSubject,
    load: crate::report::db_subject::HitLoad,
    ordinal: u32,
    seen_by_scope: bool,
    program: &str,
) -> Result<()> {
    if let crate::report::db_subject::HitLoad::Shown(shown) = load {
        // The holder's blob is loaded by the id: where no line of it holds the id as
        // `GetBioseqNoData`'s filter takes it, NCBI's prefetch stops (`other_bioseq`); else the
        // prefetch looks its taxids up instead.
        let other = match shown.holder {
            Some(holder) => Some(
                db.defline_set(holder.ordinal())
                    .map_err(|error| database_error(error, program))?,
            ),
            None => None,
        };
        let stops = other
            .as_ref()
            .is_some_and(|set| set.bioseq_defline(subject.alignment_id()).is_none());
        let differs = stops
            || (db.taxdb_found()
                && other.as_ref().is_none_or(|set| {
                    lookup_texts_differ(db, &[&subject.volume_taxids(), &defline_set_taxids(set)])
                }));
        if differs {
            let reason = subject.shown_reason(shown);
            let reason = if stops {
                format!(
                    "{reason}; NCBI BLAST+ 2.17.0 can read that sequence for it, whose definition lines do not hold the Seq-align id, and then stops (oid headers do not contain target gi/seq_id)"
                )
            } else {
                reason
            };
            return Err(parsed_refusal(
                db.oid(ordinal),
                &format!("{reason}, {PREFETCH_READS_IT}"),
                program,
            ));
        }
        return Ok(());
    }
    // A lookup that throws on the id is the prefetch's warning, one that fails on a file with
    // NCBI's exception its stop (`batch_prefetch`).
    if subject.prefetch_fails() || subject.prefetch_stops().is_some() {
        return Ok(());
    }
    match prefetch_reads(db, ordinal, subject, false, seen_by_scope, program)? {
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
/// with the reason of LOSAT's refusal where that depends on what LOSAT does not know (NCBI's
/// places give another sequence for the id at the prefetch, `load`, but LOSAT does not know
/// which: `HitLoad::Shown` without a holder, or in an order LOSAT does not know) or where
/// `GetBioseq` throws (NCBI's prefetch stops with the `CSeqDBException`). A holder that the
/// places give is loaded by the line that holds the id's handle, whose Bioseq holds the id: a
/// handle. Where the batch has a Seq-align of it on the query's minus strand and the subject's
/// plus strand, a null handle fails NCBI's prefetch ([`batch_prefetch`]; review 2 of line DBID
/// pssmread, H5).
fn prefetch_gets_no_bioseq(
    db: &SeqDb,
    ordinal: u32,
    subject: &crate::report::db_subject::ParsedSubject,
    load: Option<crate::report::db_subject::HitLoad>,
    program: &str,
) -> Result<std::result::Result<bool, String>> {
    use crate::report::db_subject::{HitLoad, ShownBy};
    if let Some(HitLoad::Shown(shown)) = load {
        return Ok(match (shown.holder, shown.by) {
            (Some(_), ShownBy::Holder | ShownBy::Weak) => Ok(false),
            _ => Err(format!(
                "{}: NCBI BLAST+ 2.17.0's prefetch of a minus-strand alignment then writes Error pre-fetching sequence data or not, as its reads decide",
                subject.shown_reason(shown)
            )),
        });
    }
    let read = scope_bioseq(db, ordinal, subject, program)?;
    if let Some(other) = read.blob.filter(|_| read.throws) {
        return Ok(Err(format!(
            "has the Seq-align id {}, by which NCBI BLAST+'s lookup gives database sequence {}, whose definition lines do not hold that id: NCBI BLAST+ 2.17.0's prefetch stops there (oid headers do not contain target gi/seq_id)",
            String::from_utf8_lossy(&subject.alignment_fasta()),
            db.oid(other)
        )));
    }
    Ok(Ok(!read.bioseq))
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
    let none = ScopeBioseq {
        blob: None,
        bioseq: false,
        throws: false,
    };
    match subject.found() {
        Some(FoundInstead::Sequence(other)) => {
            let set = db
                .defline_set(other.ordinal())
                .map_err(|error| database_error(error, program))?;
            Ok(other_bioseq(&set, other.ordinal(), subject).unwrap_or(none))
        }
        Some(FoundInstead::Nothing | FoundInstead::ThrowsOnId(_) | FoundInstead::Failed(_)) => {
            Ok(none)
        }
        // The loader does not keep the sequence that the lookup gives for the id's handle
        // (`ParsedSubject::scope_finds_none`): no blob is loaded.
        _ if subject.scope_finds_none() => Ok(none),
        None | Some(FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(..)) => {
            Ok(ScopeBioseq {
                blob: Some(ordinal),
                bioseq: subject.bioseq_holds_alignment_id(),
                throws: false,
            })
        }
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
    let id = subject.loader_id();
    let id = &*id;
    let kept = set.bioseq_defline(id);
    Some(ScopeBioseq {
        blob: Some(other),
        bioseq: kept.is_some_and(|line| {
            set.deflines[line].ids.iter().any(|held| {
                crate::seqid::defline_set::compare(held, id) == crate::seqid::defline_set::Sic::Yes
            })
        }),
        throws: kept.is_none(),
    })
}

/// Whether NCBI's data loader keeps the sequence of the definition lines `set` that the lookup of
/// the Seq-align id of `subject` gives (`CBlastDbDataLoader::x_GetOid`, bdbloader.cpp): one of its
/// ids compares equal to the Seq-id of that id's handle (`ParsedSubject::loader_id`).
fn loader_keeps(
    set: &crate::seqid::DeflineSet,
    subject: &crate::report::db_subject::ParsedSubject,
) -> bool {
    let loader_id = subject.loader_id();
    set.deflines
        .iter()
        .flat_map(|defline| &defline.ids)
        .any(|id| {
            crate::seqid::defline_set::compare(id, &loader_id)
                == crate::seqid::defline_set::Sic::Yes
        })
}

/// Whether the batches of a `-db` search prefetch their reported subjects though the first format
/// `format` writes the Seq-align ids alone: outfmt 8, 9 and 12 with
/// `PRE_FETCH_SEQS_LIMIT=2147483647` ([`prefetches_subjects`]; BLASTP's block search declines
/// such a run, as it does a table of fields that read nothing of the subject on a database with
/// parsed Seq-ids, whose definition lines it never takes).
pub(crate) fn seq_annots_prefetch(format: Option<PrefetchFormat>) -> bool {
    format.is_some_and(|format| matches!(format.outfmt(), 8 | 9 | 12))
        && prefetches_subjects(format)
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262-266 (cited at
// `PrefetchFormat`)
/// Whether the first format `format` of a call is one whose batches prefetch their reported
/// subjects where the limit allows ([`PrefetchFormat::reads_seq_aligns`]): outfmt 0 to 7 and 13
/// to 18 (`s_IsPrefetchFormat` leaves out 8, 9, 11 and 12 alone), with
/// `PRE_FETCH_SEQS_LIMIT=2147483647` the Seq-annots 8, 9 and 12 too, and there psiblast's archive
/// 11, whose rounds prefetch before their test of the archive. The posts of a batch's prefetch
/// (`prefetch_messages`) go with the query warnings of its first query, written with the first
/// format (the Seq-annots write only the reader's, the archive keeps them as its messages): after
/// the prolog that NCBI wrote before the batch loop, the XML2 and JSON documents' start too
/// (`CBlastFormat::PrintProlog`, before the reports of the batch; DECISIONS 48 (1) said these
/// formats prefetch nothing, which `s_IsPrefetchFormat` does not say: line DBID taxpf-1). The
/// archive of the other applications, which write it without calling the prefetch, prefetches
/// nothing.
pub(crate) fn prefetches_subjects(format: Option<PrefetchFormat>) -> bool {
    // `s_IsPrefetchFormat` takes the query-anchored formats 1 to 4 as 0 (DECISIONS T3-4); the
    // limit's rule of the Seq-annots and of psiblast's archive is that of the prefetch's reads.
    format.is_some_and(|format| {
        matches!(format.outfmt(), 0..=9 | 11..=17) && format.reads_seq_aligns()
    })
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

/// NCBI's texts of opening the taxonomy database where the lookup of a `Hit` opens it (the run's
/// first lookup, where no load of the reported subjects asked the database before,
/// [`TaxdbRun::open_for_rows`]): NCBI's XML2 and JSON writers make a query's whole document, and
/// its lookups, before they write it, so the texts come before the document, where the first
/// `Hit` that asks is made; they are kept here and written by [`TaxdbPrinter::at_hit`]. The
/// texts of the loads are the batch's prefetch's, or of the reads of the `Hit`s where it is off
/// ([`read_reported_taxonomy`]).
#[derive(Default)]
pub(crate) struct TaxdbNotes {
    /// The program's name in lower case, as the messages show it.
    program: String,
    /// The texts of opening the database, when a `Hit` of these reports opened it.
    open: Vec<String>,
}

impl TaxdbNotes {
    /// A printer that writes these notes as the `Hit`s of a report come.
    pub(crate) fn printer(&self) -> TaxdbPrinter<'_> {
        TaxdbPrinter {
            notes: self,
            opened: false,
        }
    }
}

/// Writes the taxonomy database's messages of a report in NCBI's order: the texts of opening at
/// the first `Hit` that asks, before the messages of the `Hit`'s own lookup.
pub(crate) struct TaxdbPrinter<'a> {
    notes: &'a TaxdbNotes,
    opened: bool,
}

impl TaxdbPrinter<'_> {
    /// NCBI's messages at the `Hit` of a subject of taxid `taxid` (`taxonomy` what its lookup
    /// said): the texts of opening, where a `Hit` first asks, then the messages of the `Hit`'s
    /// lookup (`CSeqDB::GetTaxInfo`, every time), to `sink`.
    pub(crate) fn at_hit(
        &mut self,
        taxid: i32,
        taxonomy: &crate::report::blast_xml::HitTaxonomy,
        sink: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        self.write(taxid != 0, std::iter::once(taxonomy), sink)
    }

    /// `at_hit` of the `Hit` of a subject named `names`: a subject with parsed Seq-ids has the
    /// lookup of each `HitDescr`.
    pub(crate) fn at_names(
        &mut self,
        names: &crate::report::blast_xml::subjects::SubjectNames,
        sink: &mut dyn std::io::Write,
    ) -> std::io::Result<()> {
        use crate::report::blast_xml::HitTaxonomy;
        if names.descriptions.is_empty() {
            return self.at_hit(names.taxid, &names.taxonomy, sink);
        }
        let asks = !names.volume_taxids.is_empty()
            || names.descriptions.iter().any(|description| {
                description.taxid != 0 || matches!(description.taxonomy, HitTaxonomy::NotFound(_))
            });
        let lookups = names
            .descriptions
            .iter()
            .map(|description| &description.taxonomy);
        self.write(asks, lookups, sink)
    }

    /// The messages of a `Hit` that asks the taxonomy database (`asks`): opening, where a `Hit`
    /// first asks, then the messages of the `Hit`'s own lookups (`lookups`, in order).
    fn write<'t>(
        &mut self,
        asks: bool,
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
        if asks && !std::mem::replace(&mut self.opened, true) {
            self.notes.open.iter().for_each(|line| push(line));
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
