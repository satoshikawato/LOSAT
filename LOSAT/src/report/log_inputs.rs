//! Single owner of the input files a run reads and the check that `-logfile` does not point to
//! one of them (line WQA D23, Owner decision `Q-D23=A` / WQA DECISIONS 49, 2026-11-24).
//!
//! No NCBI counterpart: NCBI BLAST+ 2.17.0 opens `-logfile` in append mode before loading its
//! configuration or reading its inputs, so naming an input file (such as `-query q.fa -logfile
//! q.fa` in `tblastx`) appends diagnostics into the input and can enter an endless read loop (5
//! min -> 10 GiB on NCBI 2.17.0 in probe `s2-query`), or corrupt `.ncbirc` on its second read
//! (`s2-ncbirc-self`). LOSAT refuses any `-logfile` that would open or create a file read by the
//! run before writing or creating it, using the same rule on the command line, in the library API,
//! and in LOSAT Web. The files a run reads include those that the strategy of
//! `-import_search_strategy` names (its database, the repeat filter's database, the window
//! masker's counts of its taxid and its MegaBLAST index: `import::StrategyReads`, line RG-F item
//! F2f), which NCBI reads after it has opened the file of `-logfile`, and the library API's inputs
//! that its argv does not name (`Supplied`).

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

use crate::blastinput::ncbi_environment::{
    log_file_name, log_path, ncbi_cwd, preparse, special_log_target, ApplicationStart, PathRules,
    Preparsed,
};
use crate::dbindex::format::{volume_name as mb_index_volume_name, SuperHeader};
use crate::ncbi_app::args::{scan_argv_kinds, Kind};
use crate::report::archive::read::import::{
    imports_strategies, kept_input, strategy_reads, KeptDatabase, MbIndexRequest, StrategyReads,
};
use crate::seqdb::{SearchPath, SeqDb, SeqType};
use crate::utils::same_file::same_file_target;

/// Data files that `rpsbproc` reads from its `-d` directory (default `data`).
const RPSBPROC_DATA_FILES: &[&str] = &[
    "cddid.tbl",
    "cdtrack.txt",
    "family_superfamily_links",
    "cddannot.dat",
    "cddannot_generic.dat",
    "bitscore_specific.txt",
];

/// The input files and database names that a run of `program` with `words` reads.
#[derive(Debug, Default)]
pub(crate) struct RunInputs {
    /// Direct file paths read by the run (command-line input files, list files, WindowMasker
    /// counts, registry files, and `rpsbproc` data files).
    pub(crate) files: Vec<PathBuf>,
    /// BLAST database arguments read by the run (`-db`, `-filtering_db`, `-dblist`, `-seqid_db`,
    /// or `-in`/`-input` with `blastdb` format), passed verbatim to `SeqDb::read_paths`.
    pub(crate) databases: Vec<OsString>,
    /// Old-style MegaBLAST index names (`blastn -index_name`, or `-db` when `-index_name` is not
    /// given), resolved via `SearchPath`.
    pub(crate) mb_indexes: Vec<OsString>,
    /// Explicit `-dbtype` (`prot` or `nucl`) for tools whose database type is not fixed by program.
    pub(crate) dbtype: Option<SeqType>,
    /// BLAST databases whose type an imported strategy fixes (`add_strategy`), which may not be
    /// the running program's (blastx running a TBLASTX strategy searches a nucleotide database).
    pub(crate) typed_databases: Vec<TypedDatabase>,
}

/// A BLAST database of [`RunInputs::typed_databases`].
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct TypedDatabase {
    /// The name, passed verbatim to `SeqDb::read_paths`.
    pub(crate) name: OsString,
    /// The type it is opened as.
    pub(crate) seqtype: SeqType,
    /// Whether the files of its volumes' MegaBLAST indexes count (a BLASTN task's database, as
    /// BLASTN's `-db` and `-filtering_db` on the command line).
    pub(crate) mb_volumes: bool,
}

/// The inputs of a library call that its argv does not name, which the call gives as `api::Inputs`
/// (its argv holds none of `-db`, `-in_pssm`, `-in_msa` and `-import_search_strategy`): the
/// command line's words for them. The command line gives none (`Supplied::default`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Supplied<'a> {
    /// The database (`Inputs::database`), the command line's `-db`.
    pub(crate) database: Option<&'a OsStr>,
    /// The query's PSSM or alignment file, with its option (`in_pssm`, `in_msa`).
    pub(crate) query_file: Option<(&'static str, &'a Path)>,
    /// What the strategy of `Inputs::strategy` makes the run read (its import's
    /// `ImportRun::reads`), the command line's `-import_search_strategy FILE`.
    pub(crate) strategy: Option<&'a StrategyReads>,
}

/// Checks whether `-logfile` `log_name` on `Disk` would open or create an input file read by the
/// run of `program` with `words` and `start`, returning `Err(reason)` (without `is not supported`
/// suffix, matching `RunFiles::log_target`) if refused.
///
/// Fast path: when `log_name` is `None` or a special target (`-`, `""`, `/dev/null`), returns
/// `Ok(())` immediately without scanning `words` or touching the filesystem.
pub(crate) fn check_disk_logfile<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    start: &ApplicationStart,
    log_name: Option<&OsStr>,
) -> Result<(), String> {
    check_disk_logfile_with(program, words, Supplied::default(), start, log_name)
}

/// [`check_disk_logfile`] of a run that also reads the inputs `supplied` (a library call's).
pub(crate) fn check_disk_logfile_with<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    supplied: Supplied<'_>,
    start: &ApplicationStart,
    log_name: Option<&OsStr>,
) -> Result<(), String> {
    if special_log_target(log_name).is_some() {
        return Ok(());
    }
    let Some(raw_name) = log_name else {
        return Ok(());
    };
    let candidates = disk_log_candidates(program, raw_name);
    if candidates.is_empty() {
        return Ok(());
    }
    let inputs = collect_run_inputs(program, words, supplied, start);
    if matches_any_input(program, &candidates, &inputs, start) {
        return Err(format!(
            "a -logfile file '{}' that is also an input file read by this run (appending \
             diagnostics to an input file corrupts it and can cause an endless read loop)",
            raw_name.display()
        ));
    }
    Ok(())
}

/// The filesystem paths that opening `-logfile raw_name` on `Disk` would open or create for
/// writing:
/// - NCBI's `log_file_name` path (when its checked directory exists and the path is not a
///   directory, matching `ncbi_environment::open_log_file`);
/// - for `makeblastdb`, also the raw `raw_name` path opened by `makeblastdb::log_of` (which opens
///   the raw OS path directly in addition to `diag_target`).
fn disk_log_candidates(program: &str, raw_name: &OsStr) -> Vec<PathBuf> {
    let mut candidates = Vec::new();
    let cwd = ncbi_cwd();
    if let Some(file) = log_file_name(
        PathRules::HOST,
        raw_name.as_encoded_bytes(),
        cwd.as_deref().map(OsStr::as_encoded_bytes),
    ) {
        let dir_ok = file
            .directory
            .as_ref()
            .is_none_or(|dir| !dir.is_empty() && log_path(dir).is_dir());
        if dir_ok {
            let path = log_path(&file.path).to_path_buf();
            if !path.is_dir() {
                candidates.push(path);
            }
        }
    }
    if program == "makeblastdb" && !raw_name.is_empty() && raw_name != "-" {
        let raw = PathBuf::from(raw_name);
        let parent_ok = raw
            .parent()
            .is_none_or(|dir| dir.as_os_str().is_empty() || dir.is_dir());
        if parent_ok && !raw.is_dir() && !candidates.contains(&raw) {
            candidates.push(raw);
        }
    }
    candidates
}

/// Collects all direct input file paths and BLAST database names that `program` with `words`,
/// `supplied` and `start` reads. Version requests (`-version`, `-version-full`,
/// `-version-full-xml`, `-version-full-json`) exit before reading command-line inputs, so only the
/// registry inputs recorded by `start` apply.
pub(crate) fn collect_run_inputs<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    supplied: Supplied<'_>,
    start: &ApplicationStart,
) -> RunInputs {
    let mut inputs = RunInputs::default();
    for path in start.registry_inputs() {
        push_file(&mut inputs.files, path.clone());
    }
    let word_refs: Vec<&OsStr> = words.iter().map(AsRef::as_ref).collect();
    if matches!(preparse(&word_refs), Preparsed::Version(_)) {
        return inputs;
    }
    match program {
        "windowmasker" => collect_windowmasker_inputs(words, &mut inputs),
        "makeblastdb" => collect_makeblastdb_inputs(words, &mut inputs),
        "rpsbproc" => collect_rpsbproc_inputs(words, &mut inputs),
        "makembindex" => collect_makembindex_inputs(words, &mut inputs),
        "dustmasker" | "segmasker" => collect_masker_inputs(program, words, &mut inputs),
        _ => collect_general_inputs(program, words, supplied, start, &mut inputs),
    }
    inputs
}

/// Collects input files for `windowmasker` (single owner shared by CLI D23 and LOSAT Web D19/D19b):
/// returns `(for_logfile, for_out, databases, databases_for_out)` where `for_out` is the subset
/// opened and read after opening `-out` (`-in` of `-mk_counts`, including entries read from
/// `-fa_list true`, and `-ustat` of masking) and `databases_for_out` is true when `-mk_counts`
/// without `-convert` reads its `-infmt blastdb` entries after opening `-out`.
pub(crate) fn windowmasker_read_parts<S: AsRef<OsStr>>(
    words: &[S],
) -> (Vec<PathBuf>, Vec<PathBuf>, Vec<OsString>, bool) {
    let mut mk_counts = false;
    let mut convert = false;
    let mut fa_list = false;
    let mut infmt_blastdb = false;
    let mut input: Option<OsString> = None;
    let mut ustat: Option<OsString> = None;
    let mut ids: Option<OsString> = None;
    let mut exclude_ids: Option<OsString> = None;
    scan_argv_kinds("windowmasker", words, |name, _, value| {
        match (name, value) {
            ("mk_counts", None) => mk_counts = true,
            ("convert", None) => convert = true,
            ("fa_list", Some(v)) => fa_list = parse_ncbi_bool(v),
            ("infmt", Some(v)) => infmt_blastdb = v == "blastdb",
            ("in", Some(v)) => input = Some(v.to_os_string()),
            ("ustat", Some(v)) => ustat = Some(v.to_os_string()),
            ("ids", Some(v)) => ids = Some(v.to_os_string()),
            ("exclude_ids", Some(v)) => exclude_ids = Some(v.to_os_string()),
            _ => {}
        }
    });
    let mut for_log = Vec::new();
    let mut for_out = Vec::new();
    let mut databases = Vec::new();
    let databases_for_out = mk_counts && !convert;
    if let Some(ref raw_in) = input {
        if mk_counts && !convert {
            if fa_list {
                push_non_stdin_file(&mut for_log, raw_in);
            }
            for entry in crate::windowmasker::counts::file_list(raw_in, fa_list) {
                let entry_os = crate::utils::os_bytes::os_string_from_bytes(entry);
                if infmt_blastdb {
                    push_db_list(&mut databases, &entry_os);
                } else {
                    push_non_stdin_file(&mut for_log, &entry_os);
                    push_non_stdin_file(&mut for_out, &entry_os);
                }
            }
        } else if infmt_blastdb && !convert {
            push_db_list(&mut databases, raw_in);
        } else {
            push_non_stdin_file(&mut for_log, raw_in);
        }
    }
    if !mk_counts && !convert {
        if let Some(ref u) = ustat {
            push_non_empty_file(&mut for_log, u);
            push_non_empty_file(&mut for_out, u);
        }
    }
    if !convert {
        for list in [ids.as_ref(), exclude_ids.as_ref()].into_iter().flatten() {
            push_non_empty_file(&mut for_log, list);
        }
    }
    (for_log, for_out, databases, databases_for_out)
}

fn collect_windowmasker_inputs<S: AsRef<OsStr>>(words: &[S], inputs: &mut RunInputs) {
    let (for_log, _for_out, databases, _databases_for_out) = windowmasker_read_parts(words);
    for path in for_log {
        push_file(&mut inputs.files, path);
    }
    for db in databases {
        if !inputs.databases.contains(&db) {
            inputs.databases.push(db);
        }
    }
}

/// Collects input files for `makeblastdb`: `-in` (quoted-space split or BLAST DB when
/// `-input_type blastdb`), `-taxid_map`, and comma-separated `-mask_data`.
fn collect_makeblastdb_inputs<S: AsRef<OsStr>>(words: &[S], inputs: &mut RunInputs) {
    let mut input_type_blastdb = false;
    let mut raw_in: Vec<OsString> = Vec::new();
    scan_argv_kinds("makeblastdb", words, |name, _, value| {
        let Some(value) = value else {
            return;
        };
        match name {
            "input_type" => input_type_blastdb = value == "blastdb",
            "dbtype" => inputs.set_dbtype(value),
            "in" => raw_in.push(value.to_os_string()),
            "taxid_map" => push_non_stdin_file(&mut inputs.files, value),
            "mask_data" => {
                for part in split_os_byte(value, b',') {
                    let trimmed = trim_ascii_spaces(&part);
                    push_non_empty_file(&mut inputs.files, &trimmed);
                }
            }
            _ => {}
        }
    });
    for val in raw_in {
        if input_type_blastdb {
            inputs.add_db(&val);
        } else {
            for part in split_quoted_os(&val) {
                inputs.add_non_stdin(&part);
            }
        }
    }
}

/// Collects input files for `rpsbproc`: `-i` (input archive/XML) and the CDD data files in `-d`
/// (default `data`).
fn collect_rpsbproc_inputs<S: AsRef<OsStr>>(words: &[S], inputs: &mut RunInputs) {
    let mut data_dir: Option<OsString> = None;
    scan_argv_kinds("rpsbproc", words, |name, _, value| {
        let Some(value) = value else {
            return;
        };
        match name {
            "i" => push_non_stdin_file(&mut inputs.files, value),
            "d" => data_dir = Some(value.to_os_string()),
            _ => {}
        }
    });
    let base = data_dir
        .filter(|d| !d.is_empty())
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("data"));
    for file in RPSBPROC_DATA_FILES {
        push_file(&mut inputs.files, base.join(file));
    }
}

/// Collects input files for `makembindex`: `-input` (file when `-iformat fasta` or
/// `-old_style_index`, else BLAST DB).
fn collect_makembindex_inputs<S: AsRef<OsStr>>(words: &[S], inputs: &mut RunInputs) {
    let mut iformat: Option<OsString> = None;
    let mut old_style = false;
    let mut raw_input: Vec<OsString> = Vec::new();
    scan_argv_kinds("makembindex", words, |name, _, value| match (name, value) {
        ("iformat", Some(v)) => iformat = Some(v.to_os_string()),
        ("old_style_index", Some(v)) => old_style = parse_ncbi_bool(v),
        ("old_style_index", None) => old_style = true,
        ("input", Some(v)) => raw_input.push(v.to_os_string()),
        _ => {}
    });
    let is_blastdb = match iformat.as_deref() {
        Some(v) => v == "blastdb",
        None => !old_style,
    };
    for val in raw_input {
        if is_blastdb {
            inputs.add_db(&val);
        } else {
            inputs.add_non_stdin(&val);
        }
    }
}

/// Collects input files for `dustmasker` and `segmasker`: `-in` (file or BLAST DB when
/// `-infmt blastdb`).
fn collect_masker_inputs<S: AsRef<OsStr>>(program: &str, words: &[S], inputs: &mut RunInputs) {
    let mut infmt_blastdb = false;
    let mut masker_inputs = Vec::new();
    scan_argv_kinds(program, words, |name, _, value| {
        let Some(value) = value else {
            return;
        };
        match name {
            "infmt" => infmt_blastdb = value == "blastdb",
            "in" => masker_inputs.push(value.to_os_string()),
            _ => {}
        }
    });
    for val in masker_inputs {
        if infmt_blastdb {
            inputs.add_db(&val);
        } else {
            inputs.add_non_stdin(&val);
        }
    }
}

/// Collects input files for search programs (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`,
/// `psiblast`, `rpsblast`, `rpstblastn`) and all remaining tools (`blastdbcmd`,
/// `blastdb_aliastool`, `blastdbcheck`, `convert2blastmask`, `makeprofiledb`, and any future
/// programs whose input files use `Kind::InputFile` or `-archive`), with the library call's inputs
/// `supplied` read as the command line's words for them. A run that imports a strategy
/// ([`StrategyRun`]) reads of the command line only what the import keeps, and also what the
/// strategy makes it read (`strategy_names`, the library call's `Supplied::strategy`).
fn collect_general_inputs<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    supplied: Supplied<'_>,
    start: &ApplicationStart,
    inputs: &mut RunInputs,
) {
    let is_blastn = program == "blastn";
    let strategy_run = StrategyRun::of(program, words, supplied);
    let mut reads_formatter_rc = false;
    let mut blastn_dbs: Vec<OsString> = Vec::new();
    let mut kept_dbs: Vec<OsString> = Vec::new();
    let mut visit = |name: &'static str, kind: Kind, value: Option<&OsStr>| {
        if name == "html" {
            reads_formatter_rc = true;
            return;
        }
        let Some(value) = value else {
            return;
        };
        if name == "conffile" || name == "logfile" {
            return;
        }
        if name == "outfmt" {
            if let Some(text) = value.to_str() {
                if text.trim_start().starts_with("18") {
                    reads_formatter_rc = true;
                }
            }
            return;
        }
        if let Some(run) = &strategy_run {
            if !kept_input(name, run.database_given) {
                return;
            }
            if name == "db" {
                // Opened as the type of the strategy's database (`strategy_database`).
                push_db_list(&mut kept_dbs, value);
                if is_blastn {
                    push_db_list(&mut blastn_dbs, value);
                }
                return;
            }
        }
        if kind == Kind::InputFile {
            push_non_stdin_file(&mut inputs.files, value);
            return;
        }
        match name {
            "gilist" | "negative_gilist" | "ipglist" | "negative_ipglist" | "seqidlist"
            | "negative_seqidlist" | "taxidlist" | "negative_taxidlist" | "seqid_file_info"
            | "archive" => {
                inputs.add_file(value);
            }
            "window_masker_db" if is_blastn => {
                inputs.add_file(value);
            }
            "window_masker_taxid" if is_blastn => {
                if let Some(taxid) = value.to_str().and_then(|s| s.parse::<i32>().ok()) {
                    inputs.add_window_masker_taxid(taxid, start);
                }
            }
            "dbtype" => {
                inputs.set_dbtype(value);
            }
            "db" | "dblist" | "seqid_db" => {
                inputs.add_db(value);
                if is_blastn && name == "db" {
                    push_db_list(&mut blastn_dbs, value);
                }
            }
            "filtering_db" if is_blastn => {
                inputs.add_db(value);
            }
            "index_name" if is_blastn => {
                inputs.add_mb_index(value);
            }
            _ => {}
        }
    };
    scan_argv_kinds(program, words, &mut visit);
    if let Some(database) = supplied.database {
        visit("db", Kind::Text, Some(database));
    }
    if let Some((option, path)) = supplied.query_file {
        visit(option, Kind::InputFile, Some(path.as_os_str()));
    }
    let mut run_index = None;
    if let Some(run) = strategy_run {
        let reads = supplied.strategy.cloned().or_else(|| {
            run.path
                .and_then(|path| strategy_names(program, &path, run.database_given))
        });
        for db in &kept_dbs {
            match &reads {
                Some(reads) => inputs.add_typed_db(db, reads.database_type, is_blastn),
                None => inputs.add_db(db),
            }
        }
        if let Some(reads) = &reads {
            inputs.add_strategy(reads, start);
        }
        let forced = reads.and_then(|reads| reads.mb_index);
        run_index = run.mb_index.run_index(forced).ok().flatten();
    }
    // BLASTN's `-db` counts its old-style index files unless `-index_name` names the index (an
    // imported strategy's run: unless its index is the old-style one of `-index_name`).
    if is_blastn
        && inputs.mb_indexes.is_empty()
        && !run_index.as_ref().is_some_and(|&(_, old_style)| old_style)
    {
        for db in blastn_dbs {
            inputs.add_mb_index(&db);
        }
    }
    if let Some((name, old_style)) = run_index {
        inputs.add_run_index(&name, old_style);
    }
    if reads_formatter_rc {
        push_file(&mut inputs.files, PathBuf::from(".ncbirc"));
        if let Some(fmtcfg) = std::env::var_os("FMTCFG").filter(|v| !v.is_empty()) {
            push_file(&mut inputs.files, PathBuf::from(&fmtcfg));
            if let Some(ncbi) = std::env::var_os("NCBI").filter(|v| !v.is_empty()) {
                push_file(&mut inputs.files, PathBuf::from(ncbi).join(&fmtcfg));
            }
        }
    }
}

/// The command line's arguments of a run that imports a strategy: the strategy of
/// `-import_search_strategy` (or the library call's `Supplied::strategy`) of a program whose
/// command line imports strategies (`import::imports_strategies`, which `strategy_reads` checks),
/// whether or not the strategy reads. The import keeps of the command line only the inputs of
/// `import::kept_input` and the MegaBLAST index of `MbIndexRequest::of` (`import::Kept`).
struct StrategyRun {
    /// The file of `-import_search_strategy` (`None` for the library call's strategy).
    path: Option<OsString>,
    /// `KeptDatabase::has_been_set`: the command line gives `-db` (or the library call
    /// `Inputs::database`) or `-subject`.
    database_given: bool,
    /// BLASTN's `CMbIndexArgs` of the command line (`MbIndexRequest::None` for the other
    /// programs).
    mb_index: MbIndexRequest,
}

impl StrategyRun {
    /// The strategy run of `program` with `words` and the library call's inputs `supplied`, `None`
    /// for a run that imports no strategy.
    fn of<S: AsRef<OsStr>>(program: &str, words: &[S], supplied: Supplied<'_>) -> Option<Self> {
        let mut path = None;
        let mut db = supplied.database.map(OsStr::to_os_string);
        let mut subject_given = false;
        let mut use_index = None;
        let mut task = None;
        let mut index_name = None;
        scan_argv_kinds(program, words, |name, _, value| {
            let Some(value) = value else {
                return;
            };
            match name {
                "import_search_strategy" => path = Some(value.to_os_string()),
                "db" => db = Some(value.to_os_string()),
                "subject" => subject_given = true,
                "use_index" => use_index = value.to_str().and_then(crate::ncbi_str::string_to_bool),
                "task" => task = Some(value.to_os_string()),
                "index_name" => index_name = Some(value.to_os_string()),
                _ => {}
            }
        });
        let imports =
            (path.is_some() && imports_strategies(program)) || supplied.strategy.is_some();
        if !imports {
            return None;
        }
        let mb_index = if program == "blastn" {
            MbIndexRequest::of(
                use_index,
                task.as_deref()
                    .and_then(OsStr::to_str)
                    .unwrap_or("megablast"),
                subject_given,
                index_name.as_deref(),
                db.as_deref(),
            )
        } else {
            MbIndexRequest::None
        };
        Some(Self {
            path: path.filter(|_| supplied.strategy.is_none()),
            database_given: KeptDatabase::has_been_set(db.is_some(), subject_given),
            mb_index,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:843-879
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:606-628
/// What the strategy file `path` of `program`'s `-import_search_strategy` makes the run read
/// (`import::strategy_reads`, `None` for a program whose command line does not import strategies;
/// `override_subject`: the command line gives `-db` or `-subject`), read before the file of
/// `-logfile` is opened, as NCBI's application opens it before its run reads the strategy. Only a
/// regular file is read, which a second reading finds as the import does: the standard input
/// (`-`), a pipe or a device is read once, by the import, and adds nothing here. A file that does
/// not read as a strategy adds nothing either (the import reports it).
fn strategy_names(program: &str, path: &OsStr, override_subject: bool) -> Option<StrategyReads> {
    if path.is_empty() || path == "-" {
        return None;
    }
    let path = Path::new(path);
    if !std::fs::metadata(path).is_ok_and(|metadata| metadata.is_file()) {
        return None;
    }
    let file = std::fs::File::open(path).ok()?;
    strategy_reads(program, std::io::BufReader::new(file), override_subject)
}

/// Whether any path in `candidates` refers to the same file (or missing target) as any input in
/// `inputs`.
fn matches_any_input(
    program: &str,
    candidates: &[PathBuf],
    inputs: &RunInputs,
    start: &ApplicationStart,
) -> bool {
    for candidate in candidates {
        if inputs
            .files
            .iter()
            .any(|input| same_file_target(candidate, input))
        {
            return true;
        }
    }
    if !inputs.databases.is_empty()
        || !inputs.mb_indexes.is_empty()
        || !inputs.typed_databases.is_empty()
    {
        let registry_blastdb = start.settings().ok().and_then(|s| s.blastdb);
        let search_path = SearchPath::from_environment(registry_blastdb.as_deref());
        let seqtypes = db_seqtypes(program, inputs.dbtype);
        let is_blastn = program == "blastn";
        for db in &inputs.databases {
            if matches_database(candidates, db, seqtypes, &search_path, is_blastn) {
                return true;
            }
        }
        for db in &inputs.typed_databases {
            let seqtypes = std::slice::from_ref(&db.seqtype);
            if matches_database(candidates, &db.name, seqtypes, &search_path, db.mb_volumes) {
                return true;
            }
        }
        for index in &inputs.mb_indexes {
            if matches_mb_index(candidates, index, &search_path) {
                return true;
            }
        }
    }
    false
}

/// The sequence types tried when resolving a BLAST database for `program` (`SeqDb::open_in` /
/// `SeqDb::open_guessing` order).
fn db_seqtypes(program: &str, dbtype: Option<SeqType>) -> &'static [SeqType] {
    match program {
        "blastn" | "tblastn" | "tblastx" | "dustmasker" | "windowmasker" | "makembindex" => {
            &[SeqType::Nucleotide]
        }
        "blastp" | "blastx" | "psiblast" | "rpsblast" | "rpstblastn" | "segmasker" => {
            &[SeqType::Protein]
        }
        _ => match dbtype {
            Some(SeqType::Protein) => &[SeqType::Protein],
            Some(SeqType::Nucleotide) => &[SeqType::Nucleotide],
            None => &[SeqType::Protein, SeqType::Nucleotide],
        },
    }
}

/// Whether any candidate `-logfile` path matches an alias file, alias filter list, GI-mask file,
/// volume component file, or (for `blastn`) new-style MegaBLAST index file of BLAST database
/// `db_name` resolved through `SeqDb::read_paths`. Databases that cannot be resolved fail to open
/// in the run itself and are not refused here.
fn matches_database(
    candidates: &[PathBuf],
    db_name: &OsStr,
    seqtypes: &[SeqType],
    search_path: &SearchPath,
    is_blastn: bool,
) -> bool {
    let db_bytes = crate::utils::os_bytes::os_str_bytes(db_name);
    let Some(read_paths) = seqtypes
        .iter()
        .find_map(|&seqtype| SeqDb::read_paths(db_bytes.as_ref(), seqtype, search_path).ok())
    else {
        return false;
    };
    for file in &read_paths.files {
        if matches_bytes_path(candidates, file) {
            return true;
        }
    }
    if is_blastn && matches_volume_mb_indexes(candidates, &read_paths.volumes, search_path) {
        return true;
    }
    false
}

/// Whether any candidate `-logfile` path matches a new-style MegaBLAST index super header
/// (`<volume>.shd`) or index volume (`<volume>.NN.idx`) of any resolved volume in `volumes`.
fn matches_volume_mb_indexes(
    candidates: &[PathBuf],
    volumes: &[Vec<u8>],
    search_path: &SearchPath,
) -> bool {
    for volume in volumes {
        let shd = [volume.as_slice(), b".shd"].concat();
        let shd_path = crate::utils::os_bytes::path_of_bytes(&shd);
        if matches_candidate(candidates, &shd_path) {
            return true;
        }
        let vol_count = SuperHeader::load(&shd_path).map_or(100u32, |h| h.num_vol.min(10_000));
        for vol in 0..vol_count {
            let idx = mb_index_volume_name(volume, vol);
            if matches_bytes_path(candidates, &idx)
                || search_path
                    .resolve_file(&idx)
                    .is_some_and(|resolved| matches_bytes_path(candidates, &resolved))
            {
                return true;
            }
        }
    }
    false
}

/// Whether any candidate `-logfile` path matches an old-style MegaBLAST index volume
/// (`<base>.NN.idx`, 00..=99) resolved through `search_path`.
fn matches_mb_index(candidates: &[PathBuf], index_name: &OsStr, search_path: &SearchPath) -> bool {
    let bytes = crate::utils::os_bytes::os_str_bytes(index_name);
    for token in bytes.split(|&b| b == b' ').filter(|s| !s.is_empty()) {
        let base = token.split(|&b| b == b',').next().unwrap_or(token);
        if base.is_empty() {
            continue;
        }
        for vol in 0..100u32 {
            let file = mb_index_volume_name(base, vol);
            let target = search_path.resolve_file(&file).unwrap_or(file);
            if matches_bytes_path(candidates, &target) {
                return true;
            }
        }
    }
    false
}

fn matches_bytes_path(candidates: &[PathBuf], path: &[u8]) -> bool {
    let path = crate::utils::os_bytes::path_of_bytes(path);
    matches_candidate(candidates, &path)
}

fn matches_candidate(candidates: &[PathBuf], target: &Path) -> bool {
    candidates
        .iter()
        .any(|candidate| same_file_target(candidate, target))
}

impl RunInputs {
    fn add_file(&mut self, value: &OsStr) {
        push_non_empty_file(&mut self.files, value);
    }

    fn add_non_stdin(&mut self, value: &OsStr) {
        push_non_stdin_file(&mut self.files, value);
    }

    fn add_db(&mut self, value: &OsStr) {
        push_db_list(&mut self.databases, value);
    }

    fn add_mb_index(&mut self, value: &OsStr) {
        push_db_list(&mut self.mb_indexes, value);
    }

    fn add_typed_db(&mut self, name: &OsStr, seqtype: SeqType, mb_volumes: bool) {
        let entry = TypedDatabase {
            name: name.to_os_string(),
            seqtype,
            mb_volumes,
        };
        if !name.is_empty() && name != "-" && !self.typed_databases.contains(&entry) {
            self.typed_databases.push(entry);
        }
    }

    /// The counts file of BLASTN's window masker for `taxid` (`-window_masker_taxid`, a strategy's
    /// `WindowMaskerTaxId`): `<path>/<taxid>/wmasker.obinary` under the registry's
    /// `[WINDOW_MASKER] WINDOW_MASKER_PATH` (`.` without one); none for 0.
    fn add_window_masker_taxid(&mut self, taxid: i32, start: &ApplicationStart) {
        if taxid == 0 {
            return;
        }
        let base = start
            .settings()
            .ok()
            .and_then(|s| s.window_masker_registry.ok())
            .and_then(|r| r.path)
            .filter(|p| !p.is_empty())
            .map(crate::utils::os_bytes::path_from_bytes)
            .unwrap_or_else(|| PathBuf::from("."));
        push_file(
            &mut self.files,
            base.join(taxid.to_string()).join("wmasker.obinary"),
        );
    }

    /// The names that an imported strategy makes the run read (`reads`): its database as the type
    /// of its program's subjects, the repeat filter's database and the window masker's counts of
    /// its taxid. Its MegaBLAST index is the run's unless the command line names one
    /// (`add_run_index`).
    fn add_strategy(&mut self, reads: &StrategyReads, start: &ApplicationStart) {
        if let Some(database) = &reads.database {
            self.add_typed_db(OsStr::new(database), reads.database_type, reads.blastn_task);
        }
        if let Some(database) = &reads.filtering_db {
            self.add_typed_db(OsStr::new(database), SeqType::Nucleotide, reads.blastn_task);
        }
        if let Some(taxid) = reads.window_masker_taxid {
            self.add_window_masker_taxid(taxid, start);
        }
    }

    /// The MegaBLAST index `name` that a run importing a strategy searches with
    /// (`MbIndexRequest::run_index`): the old-style index files of the name, and for a new-style
    /// index (`old_style` false) also the database of the name, whose volumes' index files the
    /// loader reads before it goes on to the old-style files.
    fn add_run_index(&mut self, name: &[u8], old_style: bool) {
        let name = crate::utils::os_bytes::os_string_from_bytes(name.to_vec());
        if !old_style {
            self.add_typed_db(&name, SeqType::Nucleotide, true);
        }
        self.add_mb_index(&name);
    }

    fn set_dbtype(&mut self, value: &OsStr) {
        self.dbtype = match value.to_str() {
            Some("prot") => Some(SeqType::Protein),
            Some("nucl") => Some(SeqType::Nucleotide),
            _ => None,
        };
    }
}

fn parse_ncbi_bool(value: &OsStr) -> bool {
    value
        .to_str()
        .and_then(crate::ncbi_str::string_to_bool)
        .unwrap_or(false)
}

pub(crate) fn push_file(out: &mut Vec<PathBuf>, path: PathBuf) {
    if !path.as_os_str().is_empty() && !out.contains(&path) {
        out.push(path);
    }
}

fn push_non_empty_file(out: &mut Vec<PathBuf>, value: &OsStr) {
    if !value.is_empty() {
        push_file(out, PathBuf::from(value));
    }
}

fn push_non_stdin_file(out: &mut Vec<PathBuf>, value: &OsStr) {
    if !value.is_empty() && value != "-" {
        push_file(out, PathBuf::from(value));
    }
}

fn push_db_list(out: &mut Vec<OsString>, value: &OsStr) {
    if !value.is_empty() && value != "-" {
        let os = value.to_os_string();
        if !out.contains(&os) {
            out.push(os);
        }
    }
}

/// Splits an `OsStr` on a single ASCII delimiter byte without losing non-UTF-8 bytes.
fn split_os_byte(value: &OsStr, delim: u8) -> Vec<OsString> {
    let bytes = value.as_encoded_bytes();
    bytes
        .split(|&b| b == delim)
        .map(|slice| {
            // SAFETY: `slice` is delimited by an ASCII byte in `value`, so both ends are valid
            // encoded boundaries.
            unsafe { OsStr::from_encoded_bytes_unchecked(slice) }.to_os_string()
        })
        .collect()
}

/// Trims leading and trailing ASCII spaces from an `OsStr`.
fn trim_ascii_spaces(value: &OsStr) -> OsString {
    let bytes = crate::ncbi_str::truncate_spaces(value.as_encoded_bytes());
    // SAFETY: `truncate_spaces` trims only ASCII whitespace bytes from the ends of `value`.
    unsafe { OsStr::from_encoded_bytes_unchecked(bytes) }.to_os_string()
}

/// Splits space-separated tokens with optional double quotes (as NCBI splits `-db` and
/// `makeblastdb -in`), preserving non-UTF-8 bytes.
fn split_quoted_os(value: &OsStr) -> Vec<OsString> {
    let bytes = value.as_encoded_bytes();
    let mut parts = Vec::new();
    let mut current: Vec<u8> = Vec::new();
    let mut in_quotes = false;
    let mut at = 0;
    while at < bytes.len() {
        match bytes[at] {
            b'"' => in_quotes = !in_quotes,
            b' ' | b'\t' if !in_quotes => {
                if !current.is_empty() {
                    // SAFETY: `current` contains verbatim slices from `bytes` split only at ASCII
                    // `"` / `' '` / `'\t'`.
                    let os = unsafe {
                        OsString::from_encoded_bytes_unchecked(std::mem::take(&mut current))
                    };
                    parts.push(os);
                }
            }
            b => current.push(b),
        }
        at += 1;
    }
    if !current.is_empty() {
        // SAFETY: as above.
        parts.push(unsafe { OsString::from_encoded_bytes_unchecked(current) });
    }
    parts
}

#[cfg(test)]
#[path = "log_inputs_tests.rs"]
mod tests;
