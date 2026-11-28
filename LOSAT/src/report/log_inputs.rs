//! Single owner of the input files a run reads and the check that `-logfile` does not point to
//! one of them (line WQA D23, Owner decision `Q-D23=A` / WQA DECISIONS 49, 2026-11-24).
//!
//! No NCBI counterpart: NCBI BLAST+ 2.17.0 opens `-logfile` in append mode before loading its
//! configuration or reading its inputs, so naming an input file (such as `-query q.fa -logfile
//! q.fa` in `tblastx`) appends diagnostics into the input and can enter an endless read loop (5
//! min -> 10 GiB on NCBI 2.17.0 in probe `s2-query`), or corrupt `.ncbirc` on its second read
//! (`s2-ncbirc-self`). LOSAT refuses any `-logfile` that would open or create a file read by the
//! run before writing or creating it, using the same rule on the command line, in the library API,
//! and in LOSAT Web.

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

use crate::algorithm::blastn::blast_engine::{
    index_volume_names, old_style_names, volume_super_header,
};
use crate::blastinput::ncbi_environment::{
    log_file_name, log_path, ncbi_cwd, preparse, special_log_target, ApplicationStart, PathRules,
    Preparsed,
};
use crate::ncbi_app::args::{scan_argv_kinds, Kind};
use crate::seqdb::{DbUse, SearchPath, SeqDb, SeqType};
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
#[derive(Debug)]
pub(crate) struct RunInputs {
    /// Direct file paths read by the run (command-line input files, list files, WindowMasker
    /// counts, registry files, and `rpsbproc` data files).
    pub(crate) files: Vec<PathBuf>,
    /// BLAST database arguments read by the run with the program's sequence types and
    /// [`RunInputs::db_use`] (`-db`, `-dblist`, `-seqid_db`, or `-in`/`-input` with `blastdb`
    /// format), passed verbatim to `SeqDb::read_paths`.
    pub(crate) databases: Vec<OsString>,
    /// The databases whose sequence type and use come from where the run finds their names
    /// ([`NamedDatabase`]).
    pub(crate) named_databases: Vec<NamedDatabase>,
    /// The MegaBLAST index that a BLASTN run's loader looks for (`None`: the run reads no index).
    pub(crate) mb_index: Option<MbIndexNames>,
    /// Explicit `-dbtype` (`prot` or `nucl`) for tools whose database type is not fixed by program.
    pub(crate) dbtype: Option<SeqType>,
    /// What the run can open of the files of its databases ([`db_use`]).
    pub(crate) db_use: DbUse,
}

/// A database that the run reads with a sequence type and a use of its own: BLASTN's
/// `-filtering_db` (nucleotide, the program's use without the options, which filter `-db` alone),
/// the databases of an imported search strategy (the strategy's program's molecule) and those of
/// blast_formatter's archive (each request's).
#[derive(Debug)]
pub(crate) struct NamedDatabase {
    pub(crate) name: OsString,
    /// The sequence types tried, in order (one: the molecule is known).
    pub(crate) seqtypes: &'static [SeqType],
    pub(crate) db_use: DbUse,
}

impl RunInputs {
    fn new(db_use: DbUse) -> Self {
        Self {
            files: Vec::new(),
            databases: Vec::new(),
            named_databases: Vec::new(),
            mb_index: None,
            dbtype: None,
            db_use,
        }
    }
}

/// What the run of `program` with `words` can open of the files of its databases: what the
/// program opens ([`DbUse::of_program`]: a search, a tool that reads the sequences in OID order,
/// makeblastdb's copy, any other tool) and what its options and flags ask for
/// ([`DbUse::note_option`]). The one reader of the options, for the command line and for LOSAT
/// Web.
pub(crate) fn db_use<S: AsRef<OsStr>>(program: &str, words: &[S]) -> DbUse {
    let mut db_use = DbUse::of_program(program);
    scan_argv_kinds(program, words, |name, kind, value| {
        if value.is_some() || kind == Kind::Flag {
            db_use.note_option(name);
        }
    });
    db_use
}

/// The names of the MegaBLAST index that NCBI's handler asks the loader for
/// (`BlastnArgs::megablast_index_name`, cited at `blast_engine::mb_index::megablast_index`).
#[derive(Debug)]
pub(crate) enum MbIndexNames {
    /// The `-index_name` values: the old-style index volumes of each (the loader looks for no
    /// new-style index then).
    OldStyle(Vec<OsString>),
    /// The `-db` values: the new-style index of each database volume (its super header, and the
    /// index volumes it names), then the old-style index volumes of the value, which the loader
    /// looks for when the new style fails.
    NewStyle(Vec<OsString>),
}

/// What BLASTN's arguments say of the MegaBLAST index (`collect_general_inputs`): `-use_index`,
/// `-task`, a discontiguous template, an imported search strategy and the `-index_name` values.
#[derive(Default)]
struct IndexArgs {
    use_index: bool,
    task: Option<OsString>,
    template: bool,
    import: bool,
    index_names: Vec<OsString>,
}

impl IndexArgs {
    fn note(&mut self, name: &str, value: Option<&OsStr>) {
        match (name, value) {
            ("use_index", Some(value)) => self.use_index = parse_ncbi_bool(value),
            ("task", Some(value)) => self.task = Some(value.to_os_string()),
            ("template_type" | "template_length", Some(_)) => self.template = true,
            ("import_search_strategy", Some(_)) => self.import = true,
            ("index_name", Some(value)) => push_db_list(&mut self.index_names, value),
            _ => {}
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3381-3414
    /// The index that the loader looks for with the `-db` values `dbs`: with `-use_index true`,
    /// task megablast (the default) and a database, without a template (`megablast_index`
    /// looks for no index then), the `-index_name` values, else the `-db` values. An imported
    /// strategy may ask for an index of its own (`ImportRun::mb_index`), which the strategy file
    /// holds: the names of the command line are compared then whatever its other options say.
    fn names(self, dbs: Vec<OsString>) -> Option<MbIndexNames> {
        let asked = self.import
            || (!dbs.is_empty()
                && self.use_index
                && self.task.as_deref().is_none_or(|task| task == "megablast")
                && !self.template);
        if !asked {
            None
        } else if !self.index_names.is_empty() {
            Some(MbIndexNames::OldStyle(self.index_names))
        } else if !dbs.is_empty() {
            Some(MbIndexNames::NewStyle(dbs))
        } else {
            None
        }
    }
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
    let inputs = collect_run_inputs(program, words, start);
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
/// - for `makeblastdb`, also the raw `raw_name` path of its log (`OpenFiles::log_output`, which opens
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

/// Collects all direct input file paths and BLAST database names that `program` with `words` and
/// `start` reads. Version requests (`-version`, `-version-full`, `-version-full-xml`,
/// `-version-full-json`) exit before reading command-line inputs, so only the registry inputs
/// recorded by `start` apply.
pub(crate) fn collect_run_inputs<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    start: &ApplicationStart,
) -> RunInputs {
    let mut inputs = RunInputs::new(db_use(program, words));
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
        _ => collect_general_inputs(program, words, start, &mut inputs),
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
    let (for_log, _for_out, databases, _databases_for_out) = masker_read_parts(program, words);
    for path in for_log {
        push_file(&mut inputs.files, path);
    }
    for db in databases {
        if !inputs.databases.contains(&db) {
            inputs.databases.push(db);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:276-317 (the output stream,
// then `x_GetReader`)
// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:206-250 (`x_GetReader`, then the
// output stream)
/// Collects the input files of `dustmasker` and `segmasker` (one owner for the command line's
/// `-logfile` check, D23, and LOSAT Web's tool runs, line NP WW4c, 2026-11-28): `(for_logfile,
/// for_out, databases, databases_for_out)`. `-in` is a FASTA file (not the standard input `-`) or,
/// with `-infmt blastdb`, a database. dustmasker opens `-out` before it opens `-in`, so its input
/// is also read after `-out` opens (`for_out`, `databases_for_out`: a missing `-in` named as `-out`
/// is the file that the run made); segmasker opens `-in` first, where a missing input ends the run
/// before `-out` opens.
pub(crate) fn masker_read_parts<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
) -> (Vec<PathBuf>, Vec<PathBuf>, Vec<OsString>, bool) {
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
    let out_first = program == "dustmasker";
    let (mut for_log, mut for_out, mut databases) = (Vec::new(), Vec::new(), Vec::new());
    for val in masker_inputs {
        if infmt_blastdb {
            push_db_list(&mut databases, &val);
        } else {
            push_non_stdin_file(&mut for_log, &val);
            if out_first {
                push_non_stdin_file(&mut for_out, &val);
            }
        }
    }
    (for_log, for_out, databases, out_first)
}

/// The input file names of a tool's run in LOSAT Web (`FileReceiver::read_files` and
/// `output_read_files`), from the parts of `windowmasker_read_parts` or `masker_read_parts` of
/// `program` with `words`: `(for_logfile, for_out)`, each database given as the files that
/// reading it checks or opens for the run's use ([`db_use`], `seqdb::alias::database_files`, in
/// the search path of the run's `BLASTDB` setting), the command line's files of `-logfile`
/// ([`check_disk_logfile`]); in `for_out` too when `databases_for_out`, only the files whose
/// presence the database's resolution checks (`DatabaseFiles::looked_up`): a missing component,
/// GI-mask or alias LMDB file that `-out` makes is never read by the run's reading of the
/// sequences (line WQA D26).
pub(crate) fn tool_read_files<S: AsRef<OsStr> + Clone>(
    program: &str,
    seqtype: SeqType,
    words: &[S],
    (log_paths, out_paths, databases, databases_for_out): (
        Vec<PathBuf>,
        Vec<PathBuf>,
        Vec<OsString>,
        bool,
    ),
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let encode = |paths: Vec<PathBuf>| -> Vec<Vec<u8>> {
        paths
            .into_iter()
            .map(|p| p.into_os_string().into_encoded_bytes())
            .collect()
    };
    let mut for_log = encode(log_paths);
    let mut for_out = encode(out_paths);
    if !databases.is_empty() {
        let blastdb = run_blastdb(program, words);
        let use_of_run = db_use(program, words);
        for db in databases {
            let files = crate::seqdb::alias::database_files(
                &crate::utils::os_bytes::os_str_bytes(&db),
                seqtype,
                blastdb.as_deref(),
                use_of_run,
            );
            for path in files.looked_up.into_iter().filter(|path| !path.is_empty()) {
                if databases_for_out {
                    for_out.push(path.clone());
                }
                for_log.push(path);
            }
            for_log.extend(files.found.into_iter().filter(|path| !path.is_empty()));
        }
    }
    (for_log, for_out)
}

/// The `BLASTDB` setting of the run of `program` with `words` (its registry), along which its
/// databases are looked for.
fn run_blastdb<S: AsRef<OsStr> + Clone>(program: &str, words: &[S]) -> Option<Vec<u8>> {
    ApplicationStart::of_process(program, words)
        .settings()
        .ok()
        .and_then(|s| s.blastdb)
}

/// The files that reading each of `databases` (a name and its molecule) checks or opens for a
/// tool's run of `program` with `words` in LOSAT Web, for the run's use ([`db_use`],
/// `seqdb::alias::database_lookup_files`, in the search path of the run's `BLASTDB` setting): the
/// command line's files of `-logfile` ([`check_disk_logfile`], `RunInputs::named_databases`).
pub(crate) fn tool_database_files<S: AsRef<OsStr> + Clone>(
    program: &str,
    words: &[S],
    databases: &[(Vec<u8>, SeqType)],
) -> Vec<Vec<u8>> {
    if databases.is_empty() {
        return Vec::new();
    }
    let blastdb = run_blastdb(program, words);
    let use_of_run = db_use(program, words);
    databases
        .iter()
        .flat_map(|(name, seqtype)| {
            crate::seqdb::alias::database_lookup_files(
                name,
                *seqtype,
                blastdb.as_deref(),
                use_of_run,
            )
        })
        .filter(|path| !path.is_empty())
        .collect()
}

/// Collects input files for search programs (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`,
/// `psiblast`, `rpsblast`, `rpstblastn`) and all remaining tools (`blastdbcmd`,
/// `blastdb_aliastool`, `blastdbcheck`, `convert2blastmask`, `makeprofiledb`, and any future
/// programs whose input files use `Kind::InputFile` or `-archive`).
fn collect_general_inputs<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    start: &ApplicationStart,
    inputs: &mut RunInputs,
) {
    let is_blastn = program == "blastn";
    let mut reads_formatter_rc = false;
    let mut blastn_dbs: Vec<OsString> = Vec::new();
    let mut index_args = IndexArgs::default();
    let mut strategy: Option<OsString> = None;
    scan_argv_kinds(program, words, |name, kind, value| {
        if is_blastn {
            index_args.note(name, value);
        }
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
        if name == "import_search_strategy" {
            strategy = Some(value.to_os_string());
        }
        if kind == Kind::InputFile {
            push_non_stdin_file(&mut inputs.files, value);
            return;
        }
        match name {
            "gilist" | "negative_gilist" | "ipglist" | "negative_ipglist" | "seqidlist"
            | "negative_seqidlist" | "taxidlist" | "negative_taxidlist" | "seqid_file_info" => {
                inputs.add_file(value);
            }
            "archive" => {
                inputs.add_file(value);
                if program == "blast_formatter" && !value.is_empty() && value != "-" {
                    let archive_use = inputs.db_use;
                    for (name, seqtype) in crate::blast_formatter::archive_databases(value) {
                        inputs.named_databases.push(NamedDatabase {
                            name: crate::utils::os_bytes::os_string_from_bytes(name),
                            seqtypes: seqtype_slice(seqtype),
                            db_use: archive_use,
                        });
                    }
                }
            }
            "window_masker_db" if is_blastn => {
                inputs.add_file(value);
            }
            "window_masker_taxid" if is_blastn => {
                if let Some(taxid) = value
                    .to_str()
                    .and_then(|s| s.parse::<i32>().ok())
                    .filter(|&id| id != 0)
                {
                    let base = start
                        .settings()
                        .ok()
                        .and_then(|s| s.window_masker_registry.ok())
                        .and_then(|r| r.path)
                        .filter(|p| !p.is_empty())
                        .map(crate::utils::os_bytes::path_from_bytes)
                        .unwrap_or_else(|| PathBuf::from("."));
                    push_file(
                        &mut inputs.files,
                        base.join(taxid.to_string()).join("wmasker.obinary"),
                    );
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
                push_named_database(
                    &mut inputs.named_databases,
                    value,
                    &[SeqType::Nucleotide],
                    DbUse::of_program(program),
                );
            }
            _ => {}
        }
    });
    if is_blastn {
        inputs.mb_index = index_args.names(blastn_dbs);
    }
    if let Some(path) = strategy.filter(|path| !path.is_empty() && path.as_os_str() != "-") {
        collect_strategy_databases(program, words, &path, inputs);
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

/// The databases of a search that imports the strategy file `path` (`-import_search_strategy`,
/// not the standard input, which the check cannot read before the run): the run's own reading of
/// the arguments and its import of the file (`report::archive::read::import::apply`), whose
/// arguments name the database that the search opens, the strategy's or the command line's,
/// with the strategy's program's molecule, BLASTN's `-filtering_db` and MegaBLAST index
/// (`api::strategy::search_databases`). They take the place of the command line's `-db` and
/// `-filtering_db`; where the arguments or the import fail, the run ends before it opens a
/// database.
fn collect_strategy_databases<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    path: &OsStr,
    inputs: &mut RunInputs,
) {
    inputs.databases.clear();
    inputs.named_databases.clear();
    inputs.mb_index = None;
    let argv = [OsStr::new("losat"), OsStr::new(program)]
        .into_iter()
        .chain(words.iter().map(AsRef::as_ref))
        .map(OsStr::to_os_string);
    let Ok(cli) = crate::cli::try_parse_with_own_start(
        argv,
        &crate::blastinput::ncbi_environment::DiagTarget::Discard,
    ) else {
        return;
    };
    let Ok(bytes) = std::fs::read(path) else {
        return;
    };
    let Ok(applied) = crate::report::archive::read::import::apply(
        cli.command,
        Some(&bytes),
        &mut std::io::sink(),
    ) else {
        return;
    };
    let found = crate::api::strategy::search_databases(&applied);
    if let Some((db, searched_by)) = found.db {
        push_named_database(
            &mut inputs.named_databases,
            &db,
            db_seqtypes(searched_by.name(), None),
            inputs.db_use,
        );
    }
    if let Some(filtering_db) = found.filtering_db {
        push_named_database(
            &mut inputs.named_databases,
            &filtering_db,
            &[SeqType::Nucleotide],
            DbUse::of_program("blastn"),
        );
    }
    inputs.mb_index = found.mb_index.map(|(name, old_style)| {
        let name = vec![crate::utils::os_bytes::os_string_from_bytes(name)];
        if old_style {
            MbIndexNames::OldStyle(name)
        } else {
            MbIndexNames::NewStyle(name)
        }
    });
}

/// The sequence types of a database of `seqtype`.
fn seqtype_slice(seqtype: SeqType) -> &'static [SeqType] {
    match seqtype {
        SeqType::Protein => &[SeqType::Protein],
        SeqType::Nucleotide => &[SeqType::Nucleotide],
    }
}

/// Adds the database `name` (not empty, not `-`) of `seqtypes` and `db_use` to `out`.
fn push_named_database(
    out: &mut Vec<NamedDatabase>,
    name: &OsStr,
    seqtypes: &'static [SeqType],
    db_use: DbUse,
) {
    if !name.is_empty() && name != "-" {
        out.push(NamedDatabase {
            name: name.to_os_string(),
            seqtypes,
            db_use,
        });
    }
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
    if inputs.databases.is_empty() && inputs.named_databases.is_empty() && inputs.mb_index.is_none()
    {
        return false;
    }
    let registry_blastdb = start.settings().ok().and_then(|s| s.blastdb);
    let search_path = SearchPath::from_environment(registry_blastdb.as_deref());
    let seqtypes = db_seqtypes(program, inputs.dbtype);
    let databases = inputs
        .databases
        .iter()
        .map(|db| (db.as_os_str(), seqtypes, inputs.db_use))
        .chain(
            inputs
                .named_databases
                .iter()
                .map(|db| (db.name.as_os_str(), db.seqtypes, db.db_use)),
        );
    let mut volumes_of: Vec<(&OsStr, Vec<Vec<u8>>)> = Vec::new();
    for (db, seqtypes, db_use) in databases {
        let (files, volumes) = database_read_paths(db, seqtypes, &search_path, db_use);
        if files
            .iter()
            .any(|file| matches_bytes_path(candidates, file))
        {
            return true;
        }
        volumes_of.push((db, volumes));
    }
    inputs
        .mb_index
        .as_ref()
        .is_some_and(|index| matches_mb_index(candidates, index, &volumes_of, &search_path))
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

/// The files that reading BLAST database `db_name` checks or opens (`SeqDb::read_paths`) for each
/// of `seqtypes` that the run tries, in the order and with the rule of `SeqDb::open_guessing`,
/// and the volumes of the type that opens: `(files, volumes)`. A database that does not resolve
/// gives the files that its resolution checked, which `-logfile` would create and the run then
/// read.
fn database_read_paths(
    db_name: &OsStr,
    seqtypes: &[SeqType],
    search_path: &SearchPath,
    db_use: DbUse,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let db_bytes = crate::utils::os_bytes::os_str_bytes(db_name);
    let mut files = Vec::new();
    let volumes = SeqDb::open_guessing(seqtypes, |seqtype| {
        let read = SeqDb::read_paths(db_bytes.as_ref(), seqtype, search_path, db_use);
        files.extend(read.files);
        match read.error {
            None => Ok(read.volumes),
            Some(error) => Err(error),
        }
    })
    .unwrap_or_default();
    (files, volumes)
}

/// Whether any candidate `-logfile` path matches a file of the MegaBLAST index that the loader
/// looks for (`index`): the old-style index volumes of each name, or the new-style index of the
/// volumes of each `-db` value (`volumes_of`, from `database_read_paths`) and then its
/// old-style index volumes.
fn matches_mb_index(
    candidates: &[PathBuf],
    index: &MbIndexNames,
    volumes_of: &[(&OsStr, Vec<Vec<u8>>)],
    search_path: &SearchPath,
) -> bool {
    match index {
        MbIndexNames::OldStyle(names) => names
            .iter()
            .any(|name| matches_old_style_index(candidates, name, search_path)),
        MbIndexNames::NewStyle(dbs) => dbs.iter().any(|db| {
            let volumes = volumes_of
                .iter()
                .find(|(name, _)| *name == db.as_os_str())
                .map_or(&[][..], |(_, volumes)| volumes.as_slice());
            matches_new_style_index(candidates, volumes, search_path)
                || matches_old_style_index(candidates, db, search_path)
        }),
    }
}

/// Whether any candidate `-logfile` path matches a file of the new-style index of the database
/// `volumes` that the loader reads (`AddIndexInfo`): each volume's super header (`<volume>.shd`),
/// and where it loads, the index volumes it names along `search_path` up to the first one not
/// found (`volume_super_header`, `index_volume_names`). A volume without a super header has no
/// index volume that the loader reads. The names are a superset where the loader stops before it
/// reads them: it compares the volume's number of sequences with the super header's first, which
/// needs the volume's own index file, and the check does not read it; and the old-style names of
/// the `-db` value that [`matches_mb_index`] adds are read only after the new style fails.
fn matches_new_style_index(
    candidates: &[PathBuf],
    volumes: &[Vec<u8>],
    search_path: &SearchPath,
) -> bool {
    for volume in volumes {
        let (shd, header) = volume_super_header(volume);
        if matches_bytes_path(candidates, &shd) {
            return true;
        }
        let Some(header) = header else {
            continue;
        };
        for name in index_volume_names(volume, &header) {
            let resolved = search_path.resolve_file(&name);
            if matches_bytes_path(candidates, resolved.as_deref().unwrap_or(&name)) {
                return true;
            }
            if resolved.is_none() {
                break;
            }
        }
    }
    false
}

/// Whether any candidate `-logfile` path matches an old-style MegaBLAST index volume that the
/// loader looks for along `search_path` for the value `index_name` (`old_style_names`: every
/// volume of each name's range, `<base>.NN.idx`, 00 to 99 by default); the loader refuses a range
/// that it never finishes before it looks for that name's volumes.
fn matches_old_style_index(
    candidates: &[PathBuf],
    index_name: &OsStr,
    search_path: &SearchPath,
) -> bool {
    let bytes = crate::utils::os_bytes::os_str_bytes(index_name);
    for name in old_style_names(&bytes) {
        let Some(volumes) = name.volume_files() else {
            return false;
        };
        for (_, file) in volumes {
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
