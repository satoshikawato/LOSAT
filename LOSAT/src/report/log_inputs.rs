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

use crate::blastinput::ncbi_environment::{
    log_file_name, log_path, ncbi_cwd, preparse, special_log_target, ApplicationStart, PathRules,
    Preparsed,
};
use crate::ncbi_app::args::{scan_argv_kinds, Kind};
use crate::seqdb::{SearchPath, SeqDb, SeqType};
use crate::utils::same_file::same_file_target;

/// Nucleotide BLAST database volume, mask, and index file extensions checked for resolved volumes
/// of `-db` (and other database arguments) when `-logfile` is given.
const NUCL_VOLUME_EXTENSIONS: &[&str] = &[
    "nhr", "nin", "nsq", "ndb", "not", "nos", "ntf", "nto", "njs", "nhd", "nhi", "nnd", "nni",
    "nog", "nsd", "nsi", "nma", "nmb", "nmc", "shd", "idx", "msk",
];

/// Protein BLAST database volume, RPS, and mask file extensions checked for resolved volumes of
/// `-db` (and other database arguments) when `-logfile` is given.
const PROT_VOLUME_EXTENSIONS: &[&str] = &[
    "phr", "pin", "psq", "pdb", "pot", "pos", "ptf", "pto", "pjs", "phd", "phi", "pnd", "pni",
    "pog", "psd", "psi", "pma", "pmb", "pmc", "rps", "loo", "aux", "freq", "blocks", "msk",
];

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
    /// or `-in`/`-input` with `blastdb` format), passed verbatim to `SeqDb::volume_and_alias_paths`.
    pub(crate) databases: Vec<OsString>,
    /// Old-style MegaBLAST index names (`blastn -index_name`), resolved via `SearchPath`.
    pub(crate) mb_indexes: Vec<OsString>,
    /// Explicit `-dbtype` (`prot` or `nucl`) for tools whose database type is not fixed by program.
    pub(crate) dbtype: Option<SeqType>,
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

/// Collects all direct input file paths and BLAST database names that `program` with `words` and
/// `start` reads. Version requests (`-version`, `-version-full`, `-version-full-xml`,
/// `-version-full-json`) exit before reading command-line inputs, so only the registry inputs
/// recorded by `start` apply.
pub(crate) fn collect_run_inputs<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
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
        _ => collect_general_inputs(program, words, start, &mut inputs),
    }
    inputs
}

/// Collects input files for `windowmasker` (single owner shared by CLI D23 and LOSAT Web D19):
/// returns `(for_logfile, for_out, databases)` where `for_out` is the subset opened and read after
/// opening `-out` (`-in` of `-mk_counts` without `-fa_list` or `-infmt blastdb`, and `-ustat` of
/// masking).
pub(crate) fn windowmasker_read_parts<S: AsRef<OsStr>>(
    words: &[S],
) -> (Vec<PathBuf>, Vec<PathBuf>, Vec<OsString>) {
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
    if let Some(ref raw_in) = input {
        if infmt_blastdb {
            push_db_list(&mut databases, raw_in);
        } else if mk_counts && !convert && !fa_list {
            for part in split_os_byte(raw_in, b',') {
                push_non_stdin_file(&mut for_log, &part);
                push_non_stdin_file(&mut for_out, &part);
            }
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
    (for_log, for_out, databases)
}

fn collect_windowmasker_inputs<S: AsRef<OsStr>>(words: &[S], inputs: &mut RunInputs) {
    let (for_log, _for_out, databases) = windowmasker_read_parts(words);
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
/// programs whose input files use `Kind::InputFile` or `-archive`).
fn collect_general_inputs<S: AsRef<OsStr>>(
    program: &str,
    words: &[S],
    start: &ApplicationStart,
    inputs: &mut RunInputs,
) {
    let is_blastn = program == "blastn";
    let mut reads_formatter_rc = false;
    scan_argv_kinds(program, words, |name, kind, value| {
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
                        .map(PathBuf::from)
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
            }
            "filtering_db" if is_blastn => {
                inputs.add_db(value);
            }
            "index_name" if is_blastn => {
                inputs.add_mb_index(value);
            }
            _ => {}
        }
    });
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
    if !inputs.databases.is_empty() || !inputs.mb_indexes.is_empty() {
        let registry_blastdb = start.settings().ok().and_then(|s| s.blastdb);
        let search_path = SearchPath::from_environment(registry_blastdb.as_deref());
        let seqtypes = db_seqtypes(program, inputs.dbtype);
        for db in &inputs.databases {
            if matches_database(candidates, db, seqtypes, &search_path) {
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

/// Whether any candidate `-logfile` path matches an alias file or volume component file of BLAST
/// database `db_name` resolved through `SeqDb::volume_and_alias_paths` (NCBI `FindVolumePaths`).
/// Databases that cannot be resolved fail to open in the run itself and are not refused here.
fn matches_database(
    candidates: &[PathBuf],
    db_name: &OsStr,
    seqtypes: &[SeqType],
    search_path: &SearchPath,
) -> bool {
    let db_bytes = crate::utils::os_bytes::os_str_bytes(db_name);
    let resolved = seqtypes.iter().find_map(|&seqtype| {
        SeqDb::volume_and_alias_paths(db_bytes.as_ref(), seqtype, search_path)
            .ok()
            .map(|paths| (seqtype, paths))
    });
    let Some((seqtype, (volumes, alias_files))) = resolved else {
        return false;
    };
    let (extensions, alias_suffix) = match seqtype {
        SeqType::Nucleotide => (NUCL_VOLUME_EXTENSIONS, b".nal".as_slice()),
        SeqType::Protein => (PROT_VOLUME_EXTENSIONS, b".pal".as_slice()),
    };
    for alias in &alias_files {
        let alias_path = crate::utils::os_bytes::path_of_bytes(alias);
        if matches_candidate(candidates, &alias_path) {
            return true;
        }
        if let Some(stem) = strip_ascii_suffix(alias, alias_suffix) {
            if matches_volume_components(candidates, stem, extensions) {
                return true;
            }
        }
    }
    for volume in &volumes {
        if matches_volume_components(candidates, volume, extensions) {
            return true;
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
        for vol in 0..100u8 {
            let file = [base, format!(".{vol:02}.idx").as_bytes()].concat();
            if let Some(resolved) = search_path.resolve_file(&file) {
                let idx_path = crate::utils::os_bytes::path_of_bytes(&resolved);
                if matches_candidate(candidates, &idx_path) {
                    return true;
                }
            }
        }
    }
    false
}

fn matches_candidate(candidates: &[PathBuf], target: &Path) -> bool {
    candidates
        .iter()
        .any(|candidate| same_file_target(candidate, target))
}

fn matches_volume_components(candidates: &[PathBuf], base: &[u8], extensions: &[&str]) -> bool {
    let mut buf = Vec::with_capacity(base.len() + 8);
    buf.extend_from_slice(base);
    buf.push(b'.');
    let prefix_len = buf.len();
    for ext in extensions {
        buf.truncate(prefix_len);
        buf.extend_from_slice(ext.as_bytes());
        let path = crate::utils::os_bytes::path_of_bytes(&buf);
        if matches_candidate(candidates, &path) {
            return true;
        }
    }
    false
}

fn strip_ascii_suffix<'a>(bytes: &'a [u8], suffix: &[u8]) -> Option<&'a [u8]> {
    let len = bytes.len().checked_sub(suffix.len())?;
    let (head, tail) = bytes.split_at(len);
    tail.eq_ignore_ascii_case(suffix).then_some(head)
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
