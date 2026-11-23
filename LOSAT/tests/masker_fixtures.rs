//! `losat dustmasker`, `losat segmasker`, `losat windowmasker`, `losat convert2blastmask` and
//! `losat blast_formatter` against NCBI BLAST+ 2.17.0 dustmasker 1.0.0, segmasker 1.0.0,
//! windowmasker 1.0.0, convert2blastmask and blast_formatter: each row of
//! `tests/<program>/cases.tsv` (its arguments as shell words, run in `tests/<program>`; `< file` is
//! the standard input) gives NCBI's standard output, standard error and exit status, frozen in
//! `tests/<program>/ncbi/<name>.{out,err,rc}` (probes of line RG, items K1, K2 and K3, 2026-11-12
//! and 11-13, on ws and the local official release; windowmasker's of line NP item M5a, 2026-11-20
//! on ws, `tests/windowmasker/make_fixtures.py`; blast_formatter's of line RG item F3, 2026-11-20
//! on ws and 2026-11-23 on the local official release). A name the table repeats is compared on its
//! last row, whose output the probe kept. NCBI is not run here. Every row runs through the library
//! API too (`losat::api::run_dustmasker`, `run_segmasker`, `run_windowmasker`,
//! `run_convert2blastmask`, `run_blast_formatter`, by `examples/library_command_line`), with the
//! same bytes.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::library_examples::example_binary;
use crate::ncbi_windows_paths::{alias_list_path_on_windows, canonical, ncbi_paths_on_windows};

/// The fixtures of `program`.
fn fixtures(program: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(program)
}

/// One row of the table: its name, the words after the program's name, the file of `< file`,
/// and the variables of its first `NAME=value` words.
struct Case {
    name: String,
    words: Vec<String>,
    stdin: Option<String>,
    variables: Vec<(String, String)>,
}

/// The variables of the first words of `words` of the form `NAME=value` (`NAME` in upper case,
/// digits and `_`, not starting with a digit, such as `BLASTDB`, where the taxonomy database is,
/// or `BL2SEQ_LEGACY`), taken out of `words`: the shell's assignments before the command.
pub(crate) fn take_variables(words: &mut Vec<String>) -> Vec<(String, String)> {
    let mut variables = Vec::new();
    while let Some((name, value)) = words.first().and_then(|word| word.split_once('=')) {
        if name.is_empty()
            || name.as_bytes()[0].is_ascii_digit()
            || !name
                .bytes()
                .all(|byte| byte.is_ascii_uppercase() || byte.is_ascii_digit() || byte == b'_')
        {
            break;
        }
        variables.push((name.to_string(), value.to_string()));
        words.remove(0);
    }
    variables
}

/// The words of `line` as the probe's shell split them: blanks between words, double quotes
/// around a word with blanks or none, and `<` before the file of the standard input (no other
/// shell syntax is in the table).
pub(crate) fn shell_words(line: &str) -> (Vec<String>, Option<String>) {
    let mut words = Vec::new();
    let mut word: Option<String> = None;
    // The open quote: inside `"..."` a `'` is a letter, inside `'...'` a `"`.
    let mut quote: Option<char> = None;
    for c in line.chars() {
        match (quote, c) {
            (None, '"' | '\'') => {
                quote = Some(c);
                word.get_or_insert_with(String::new);
            }
            (Some(open), _) if c == open => quote = None,
            (None, ' ' | '\t') => {
                if let Some(done) = word.take() {
                    words.push(done);
                }
            }
            _ => word.get_or_insert_with(String::new).push(c),
        }
    }
    assert!(quote.is_none(), "unbalanced quotes: {line}");
    words.extend(word);
    let mut stdin = None;
    if let Some(at) = words.iter().position(|word| word == "<") {
        stdin = Some(words[at + 1].clone());
        words.drain(at..at + 2);
    }
    (words, stdin)
}

/// The rows of `program`'s table, the last of each name.
fn cases(program: &str) -> Vec<Case> {
    let table = fs::read_to_string(fixtures(program).join("cases.tsv")).unwrap();
    let mut cases: Vec<Case> = Vec::new();
    for line in table.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let (name, arguments) = line.split_once('\t').unwrap();
        let (mut words, stdin) = shell_words(arguments);
        let variables = take_variables(&mut words);
        cases.retain(|case| case.name != name);
        cases.push(Case {
            name: name.to_string(),
            words,
            stdin,
            variables,
        });
    }
    cases
}

/// Stdout, stderr and exit status of a run.
type Streams = (Vec<u8>, Vec<u8>, i32);

/// One run of `losat <program>` for `case` in `tests/<program>`: stdout, stderr, exit status.
/// No LOSAT, NCBI or BLAST variable, no registry file (`NCBI_CONFIG_PATH` names the fixture
/// directory, which has none).
fn run(program: &str, case: &Case) -> Streams {
    run_binary(Path::new(env!("CARGO_BIN_EXE_losat")), program, case)
}

/// `run` of `binary <program>` (the command line, or the library's example). A row whose words
/// name `@OUT@` (windowmasker's counts files) runs with a directory of its own there; the files
/// the run leaves in it follow the standard output (`with_files`), and its path in the standard
/// error is written `@OUT@`.
fn run_binary(binary: &Path, program: &str, case: &Case) -> Streams {
    let Some(dir) = out_directory(program, case) else {
        let output = command(binary, program, case).output().unwrap();
        return (output.stdout, output.stderr, output.status.code().unwrap());
    };
    let path = dir.to_string_lossy().into_owned();
    let words = case
        .words
        .iter()
        .map(|word| word.replace(OUT, &path))
        .collect();
    let named = Case {
        name: case.name.clone(),
        words,
        stdin: case.stdin.clone(),
        variables: case.variables.clone(),
    };
    let output = command(binary, program, &named).output().unwrap();
    let stdout = with_files(output.stdout, &dir);
    fs::remove_dir_all(&dir).unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr)
        .replace(&path, OUT)
        .into_bytes();
    (stdout, stderr, output.status.code().unwrap())
}

/// The word of a row that names the directory of its output files.
const OUT: &str = "@OUT@";

/// A new empty directory for a run of `case`, when its words name `OUT`.
fn out_directory(program: &str, case: &Case) -> Option<PathBuf> {
    static RUNS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    if !case.words.iter().any(|word| word.contains(OUT)) {
        return None;
    }
    let run = RUNS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "masker-out-{program}-{}-{}-{run}",
        case.name,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    Some(dir)
}

/// `stdout` followed by each file of `dir` in the order of the names: a line `@OUT@/<name>
/// <length>` after a newline, then the bytes of the file (`make_fixtures.py` writes NCBI's runs
/// the same way).
fn with_files(mut stdout: Vec<u8>, dir: &Path) -> Vec<u8> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    for name in names {
        let bytes = fs::read(dir.join(&name)).unwrap();
        stdout.extend_from_slice(format!("\n{OUT}/{name} {}\n", bytes.len()).as_bytes());
        stdout.extend_from_slice(&bytes);
    }
    stdout
}

/// The command of `run_binary`.
fn command(binary: &Path, program: &str, case: &Case) -> Command {
    let dir = fixtures(program);
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    let stdin = match &case.stdin {
        Some(file) => Stdio::from(fs::File::open(dir.join(file)).unwrap()),
        None => Stdio::null(),
    };
    command
        .env("NCBI_CONFIG_PATH", &dir)
        .envs(case.variables.iter().map(|(name, value)| (name, value)))
        .current_dir(&dir)
        .arg(program)
        .args(&case.words)
        .stdin(stdin);
    command
}

/// NCBI's stdout, stderr and exit status for `case` of `program`. NCBI's bytes come from Linux.
/// Where NCBI's code does otherwise on Windows, the row expects what it does there (the rest of
/// the text is the Linux build's, which LOSAT writes on every platform):
/// - a directory does not open as a file (`CArg_InputFile::x_Open`: the `ifstream::open` of
///   MSVC's library fails on a directory, `is_open()` is false, the exception
///   `CArgException::eNoFile`), so `-in` naming one is "File is not accessible" as a file that
///   does not exist: the row `MISSING_INPUT` of the program with the name of `case`
///   (segmasker and convert2blastmask open `-in` first; dustmasker opens `-out` first, and the
///   interval writer of its directory row writes nothing before a record; windowmasker opens
///   `-in` itself, `CWinMaskConfig`'s `ifstream`, whose failure names the file after ` - `; the test
///   `the_rows_that_read_a_directory_are_the_ones_ncbi_read_as_one` checks that the standard
///   output is the same); windowmasker's `-mk_counts` opens each file of `-in` itself
///   (`CWinMaskUtil::CInputBioseq_CI`), whose failed stream gives the read error of
///   `MK_COUNTS_MISSING_INPUT`; on Linux the directory opens and its first read fails;
/// - the paths of the database messages and the search path (`ncbi_paths_on_windows`).
fn expected(program: &str, case: &Case) -> Streams {
    if let Some((_, near, posts)) = empty_id_row(program, &case.name) {
        let (out, err, status) = expected_row(program, near);
        let mut posted: Vec<u8> = posts
            .iter()
            .flat_map(|row| expected_row(program, row).1)
            .collect();
        posted.extend(err);
        return (out, posted, status);
    }
    if cfg!(windows) {
        if let Some(directory) = input_directory(program, case) {
            let (missing, missing_name) = missing_input(program, case);
            let named = |bytes: Vec<u8>| -> Vec<u8> {
                String::from_utf8(bytes)
                    .unwrap()
                    .replace(&format!("`{missing_name}'"), &format!("`{directory}'"))
                    .replace(
                        &format!(" - {missing_name}\n"),
                        &format!(" - {directory}\n"),
                    )
                    .replace(&format!(": {missing_name}\n"), &format!(": {directory}\n"))
                    .into_bytes()
            };
            let (out, err, status) = expected_row(program, missing);
            return (named(out), named(err), status);
        }
    }
    expected_row(program, &case.name)
}

/// windowmasker's rows of `PD-LOSAT-NCBI-DEFECTS` approved exception 22: an `-ids` or
/// `-exclude_ids` list matched as text with lines that give an empty id (`in/ids_e_*.txt`,
/// make_inputs.py). NCBI BLAST+ 2.17.0 posts its `bad id` error for the first such line and then
/// crashes (SIGSEGV, exit status 139: the row's own freeze); LOSAT posts the error for each and goes
/// on. Each entry: the row, the row of its near input (the list without those lines, which NCBI
/// completes), and the rows of one such line each whose errors LOSAT posts, in the order of the
/// row's list. LOSAT's bytes are the near row's, its standard error after those errors as NCBI
/// posts them (`expected`; `the_empty_id_rows_are_ncbis_crash_and_near_input`).
const EMPTY_ID_ROWS: [EmptyIdRow; 14] = [
    ("e_ids_gt", "i_ids_a", &["e_ids_gt"]),
    ("e_ids_gtsp", "i_ids_a", &["e_ids_gtsp"]),
    ("e_ids_tab", "i_ids_a", &["e_ids_tab"]),
    ("e_ids_sp", "i_ids_a", &["e_ids_sp"]),
    ("e_ids_pipe", "i_ids_a", &["e_ids_pipe"]),
    ("e_ids_gtpipe", "i_ids_a", &["e_ids_gtpipe"]),
    ("e_ids_gtgt", "i_ids_a", &["e_ids_gtgt"]),
    ("e_ids_gtgtpipe", "i_ids_a", &["e_ids_gtgtpipe"]),
    ("e_ids_all", "i_ids_a", EVERY_EMPTY_ID),
    ("e_ids_only", "i_ids_emptyfile", EVERY_EMPTY_ID),
    ("e_excl_gt", "i_excl_a", &["e_excl_gt"]),
    ("e_excl_all", "i_excl_a", EVERY_EMPTY_ID),
    ("e_k_gt", "k_ids", &["e_k_gt"]),
    ("e_k_all", "k_ids", EVERY_EMPTY_ID),
];

/// An entry of `EMPTY_ID_ROWS`: the row, its near row, the rows of its errors.
type EmptyIdRow = (&'static str, &'static str, &'static [&'static str]);

/// The rows of one line each of `in/ids_e_all.txt` and `in/ids_e_only.txt`, in their order.
const EVERY_EMPTY_ID: &[&str] = &[
    "e_ids_gt",
    "e_ids_gtsp",
    "e_ids_tab",
    "e_ids_sp",
    "e_ids_pipe",
    "e_ids_gtpipe",
    "e_ids_gtgt",
    "e_ids_gtgtpipe",
];

/// The entry of `EMPTY_ID_ROWS` for the row `name` of `program`.
fn empty_id_row(program: &str, name: &str) -> Option<&'static EmptyIdRow> {
    EMPTY_ID_ROWS
        .iter()
        .find(|(row, _, _)| program == "windowmasker" && *row == name)
}

/// The row of each program's table whose `-in` names a file that does not exist.
const MISSING_INPUT: [(&str, &str); 4] = [
    ("dustmasker", "noin_file"),
    ("segmasker", "noin_file"),
    ("windowmasker", "m_in_missing"),
    ("convert2blastmask", "in_missing"),
];

/// windowmasker's `-mk_counts` row whose `-in` names a file that does not exist, with the
/// settings of its directory row.
const MK_COUNTS_MISSING_INPUT: &str = "k_in_missing";

/// Whether `case` is a row of windowmasker's `-mk_counts`.
fn makes_counts(program: &str, case: &Case) -> bool {
    program == "windowmasker" && case.words.iter().any(|word| word == "-mk_counts")
}

/// The name of the row of `program` whose `-in` names a missing file in the mode of `case`
/// (`MISSING_INPUT`, `MK_COUNTS_MISSING_INPUT`) and the file its `-in` names.
fn missing_input(program: &str, case: &Case) -> (&'static str, String) {
    let (_, row) = MISSING_INPUT
        .into_iter()
        .find(|(name, _)| *name == program)
        .expect("a masking program of MISSING_INPUT");
    let row = if makes_counts(program, case) {
        MK_COUNTS_MISSING_INPUT
    } else {
        row
    };
    let case = cases(program)
        .into_iter()
        .find(|case| case.name == row)
        .expect("the row of MISSING_INPUT in the table");
    let name = input(&case).expect("the -in of the row").to_string();
    (row, name)
}

/// The value of `-in` in `case`.
fn input(case: &Case) -> Option<&str> {
    let at = case.words.iter().position(|word| word == "-in")?;
    case.words.get(at + 1).map(String::as_str)
}

/// The value of `-in` in `case` when it names a directory that the program opens as a file: the
/// FASTA input (`-infmt blastdb` opens no file of that name, windowmasker's `-convert` reads
/// `-in` as a counts file: `counts_directory`), a name that is not blank (an empty name joins
/// to the fixture directory itself, and so does `" "` on Windows).
fn input_directory<'a>(program: &str, case: &'a Case) -> Option<&'a str> {
    let database = case
        .words
        .windows(2)
        .any(|pair| pair[0] == "-infmt" && pair[1] == "blastdb");
    input(case).filter(|name| {
        !database
            && !converts_counts(program, case)
            && !name.trim().is_empty()
            && fixtures(program).join(name).is_dir()
    })
}

/// Whether `case` of windowmasker names a list that NCBI reads as the lines of a text stream
/// (`-ids`, `-exclude_ids`, and `-in` with `-fa_list`) in a file that has a `\r\n`: the text
/// mode of the stream on Windows reads each `\r\n` as `\n` (`ncbi_app::files::LineEnds`), where
/// on Linux the `\r` stays in the line, so NCBI's bytes from Linux are not those of its build for
/// Windows there.
fn crlf_text_list(program: &str, case: &Case) -> bool {
    let fa_list = case.words.windows(2).any(|pair| {
        pair[0] == "-fa_list" && ["t", "true", "1"].contains(&pair[1].to_ascii_lowercase().as_str())
    });
    program == "windowmasker"
        && case.words.windows(2).any(|pair| {
            let list = matches!(pair[0].as_str(), "-ids" | "-exclude_ids")
                || (fa_list && pair[0] == "-in");
            list && fs::read(fixtures(program).join(&pair[1]))
                .is_ok_and(|bytes| bytes.windows(2).any(|two| two == b"\r\n"))
        })
}

/// Whether `case` is a row of windowmasker's `-convert`.
fn converts_counts(program: &str, case: &Case) -> bool {
    program == "windowmasker" && case.words.iter().any(|word| word == "-convert")
}

/// Whether `case` is a `-convert` row whose `-in` names a directory: NCBI's counts reader
/// (`CSeqMaskerIstatFactory`) opens it on Linux and finds no format in it; what its build for
/// Windows, whose library opens no directory, says was not probed, so the row runs on Unix.
fn counts_directory(program: &str, case: &Case) -> bool {
    converts_counts(program, case)
        && input(case)
            .is_some_and(|name| !name.trim().is_empty() && fixtures(program).join(name).is_dir())
}

/// The rows whose error quotes an `OIDLIST` value that starts with `/`, which NCBI joins to a
/// directory on Windows (`alias_list_path_on_windows`): the program, the row, that directory
/// relative to the program's fixtures and the value. The directory is the alias file's, or, where
/// the first name of its `DBLIST` is a volume, that volume's (`seqdb::alias::compute_masks`):
/// `n_miss` and `n_missl` name the alias `tfn`, `p_miss` and `p_missl` the volume `tfp`.
const OIDLIST_ROWS: &[(&str, &str, &str, &str)] = &[
    (
        "dustmasker",
        "db_oidlist_miss_interval",
        "../db_search/oidlist",
        MISSING_OIDLIST,
    ),
    (
        "dustmasker",
        "db_oidlist_missl_interval",
        "../db_search/oidlist",
        MISSING_OIDLIST,
    ),
    (
        "segmasker",
        "db_oidlist_miss_interval",
        "../db_search/taxfilter",
        MISSING_OIDLIST,
    ),
    (
        "segmasker",
        "db_oidlist_missl_interval",
        "../db_search/taxfilter",
        MISSING_OIDLIST,
    ),
];

/// The `OIDLIST` value of `tests/db_search/oidlist/{n,p}_miss{,l}` (make_oidlist.py `MISSING`).
const MISSING_OIDLIST: &str = "/losat-oidlist-fixture-missing/none.msk";

/// NCBI's stdout, stderr and exit status for the row `name` of `program`, as frozen; the standard
/// output and error in LOSAT's words (NCBI's bytes through `losat::wording`, AGENTS.md rule 8:
/// `-xmlhelp`'s version of the program, the prefixes). A row that LOSAT refuses by name has its
/// message in `refused/<name>.err` (no standard output, exit status 1), and a row of
/// `REFUSED_ROWS` that says what NCBI's record of it must show; NCBI's bytes of that row stay in
/// `ncbi/` as that record.
fn expected_row(program: &str, name: &str) -> Streams {
    if let Ok(err) = fs::read(
        fixtures(program)
            .join("refused")
            .join(format!("{name}.err")),
    ) {
        let refused = refused_row(program, name);
        assert!(
            refused.fits(&fixtures(program), &err),
            "{program} {name}: NCBI's record or LOSAT's message is not the one of a refusal of {}",
            refused.kind
        );
        return (Vec::new(), err, 1);
    }
    let file = |extension: &str| {
        fixtures(program)
            .join("ncbi")
            .join(format!("{name}.{extension}"))
    };
    let status = fs::read_to_string(file("rc")).unwrap();
    // The database rows name the probe's directory in NCBI's search path; here it is this one.
    let here = |bytes: Vec<u8>| -> Vec<u8> {
        let text = String::from_utf8_lossy(&bytes);
        let text = if cfg!(windows) {
            ncbi_paths_on_windows(&text, &["@CWD@"])
        } else {
            text.into_owned()
        };
        let dir = fixtures(program);
        text.replace("@CWD@", &dir.to_string_lossy()).into_bytes()
    };
    let err = fs::read(file("err")).unwrap();
    let err = if err.windows(5).any(|w| w == b"@CWD@") {
        here(err)
    } else {
        err
    };
    let err = match OIDLIST_ROWS
        .iter()
        .find(|(each, row, ..)| *each == program && *row == name)
    {
        Some(&(_, _, dir, value)) if cfg!(windows) => {
            let dir = canonical(&fixtures(program).join(dir));
            let text = String::from_utf8(err).unwrap();
            alias_list_path_on_windows(&text, &dir.to_string_lossy(), value).into_bytes()
        }
        _ => err,
    };
    (
        losat::wording::losat_output_bytes(&fs::read(file("out")).unwrap()).into_owned(),
        losat::wording::losat_bytes(&err).into_owned(),
        status.trim().parse().unwrap(),
    )
}

/// A row that LOSAT refuses by name (`refused/<name>.err`): the kind of its refusal, the exit
/// status and a text of NCBI's frozen record that show what NCBI's run depends on there, and a
/// text of LOSAT's message. A row cannot move to `refused/` over a result of NCBI's that LOSAT
/// should reproduce without a row here that says why.
struct RefusedRow {
    program: &'static str,
    name: &'static str,
    kind: &'static str,
    ncbi_status: i32,
    ncbi_err: &'static str,
    losat_err: &'static str,
}

/// The rows that LOSAT refuses by name.
const REFUSED_ROWS: [RefusedRow; 5] = [
    // NCBI's database reader does not open it: NCBI then asks its servers (its remote loader's
    // failure under `unshare -rn`).
    RefusedRow {
        program: "blast_formatter",
        name: "refused_database",
        kind: "Database",
        ncbi_status: 255,
        ncbi_err: "s_InitializeSubject",
        losat_err: "asks NCBI's servers",
    },
    // A stored Fatal message: NCBI's post carries a stack trace of its own binary.
    RefusedRow {
        program: "blast_formatter",
        name: "refused_fatal_message",
        kind: "FatalMessage",
        ncbi_status: 255,
        ncbi_err: "Note[F]",
        losat_err: "severity Fatal",
    },
    // An output format outside what LOSAT's blast_formatter writes (8, 9 and 12): NCBI writes it.
    RefusedRow {
        program: "blast_formatter",
        name: "refused_of6",
        kind: "OutputFormat",
        ncbi_status: 0,
        ncbi_err: "",
        losat_err: "output format 6 is not supported",
    },
    // A NaN or an infinite real number in the results: NCBI writes it in ASN.1 binary (its one
    // special octet) and stops at it in JSON with an error of its stream.
    RefusedRow {
        program: "blast_formatter",
        name: "refused_real_9",
        kind: "Unwritable",
        ncbi_status: 0,
        ncbi_err: "",
        losat_err: "output format 9 of a BLAST archive with a NaN or an infinite real number",
    },
    RefusedRow {
        program: "blast_formatter",
        name: "refused_real_12",
        kind: "Unwritable",
        ncbi_status: 1,
        ncbi_err: "invalid double: not a number",
        losat_err: "output format 12 of a BLAST archive with a NaN or an infinite real number",
    },
];

/// The row of `REFUSED_ROWS` for `name` of `program`.
fn refused_row(program: &str, name: &str) -> &'static RefusedRow {
    REFUSED_ROWS
        .iter()
        .find(|row| row.program == program && row.name == name)
        .unwrap_or_else(|| panic!("{program} {name}: refused/ without a row of REFUSED_ROWS"))
}

impl RefusedRow {
    /// Whether NCBI's frozen record of the row in `dir` and LOSAT's message `err` are those of
    /// the refusal.
    fn fits(&self, dir: &Path, err: &[u8]) -> bool {
        let ncbi = |extension: &str| {
            fs::read(dir.join("ncbi").join(format!("{}.{extension}", self.name))).unwrap()
        };
        let status: i32 = String::from_utf8(ncbi("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        let ncbi_err = String::from_utf8_lossy(&ncbi("err")).into_owned();
        let err = String::from_utf8_lossy(err);
        status == self.ncbi_status
            && ncbi_err.contains(self.ncbi_err)
            && err.starts_with("Error: ")
            && err.contains(self.losat_err)
    }
}

#[test]
fn every_refused_row_names_its_refusal_and_ncbi_s_record_fits_it() {
    for program in [
        "dustmasker",
        "segmasker",
        "windowmasker",
        "convert2blastmask",
        "blast_formatter",
    ] {
        let Ok(entries) = fs::read_dir(fixtures(program).join("refused")) else {
            continue;
        };
        for entry in entries {
            let file = entry.unwrap().file_name().to_string_lossy().into_owned();
            let name = file.strip_suffix(".err").expect("refused/<name>.err");
            let err = fs::read(fixtures(program).join("refused").join(&file)).unwrap();
            let row = refused_row(program, name);
            assert!(row.fits(&fixtures(program), &err), "{program} {name}");
        }
    }
    for row in &REFUSED_ROWS {
        assert!(
            fixtures(row.program)
                .join("refused")
                .join(format!("{}.err", row.name))
                .is_file(),
            "{} {}: a row of REFUSED_ROWS without its refused/ file",
            row.program,
            row.name
        );
    }
}

/// Whether `case` of `program` runs on this platform: a row that writes its output to a device
/// file runs where the file is, `/dev/full` on Linux, `/dev/null` on Unix (`-logfile /dev/null`
/// means no log file on every platform: NCBI's `SetLogFile` compares the name,
/// corelib/ncbidiag.cpp:7594); a `counts_directory` row runs on Unix; a `crlf_text_list` row runs
/// off Windows.
fn runs_here(program: &str, case: &Case) -> bool {
    if counts_directory(program, case) && !cfg!(unix) {
        return false;
    }
    if crlf_text_list(program, case) && cfg!(windows) {
        return false;
    }
    let out = case
        .words
        .iter()
        .position(|word| word == "-out")
        .and_then(|at| case.words.get(at + 1))
        .map(String::as_str);
    // A table, not a `match` or an `if` chain: the `cfg!` are constants, which clippy would read
    // as `matches!` (Windows) or as an `if` that is always true (Unix).
    let devices = [
        ("/dev/full", cfg!(target_os = "linux")),
        ("/dev/null", cfg!(unix)),
    ];
    devices
        .iter()
        .find(|(device, _)| out == Some(*device))
        .is_none_or(|&(_, here)| here)
}

/// The rows of `program` that run here (`runs_here`).
fn cases_here(program: &str) -> Vec<Case> {
    cases(program)
        .into_iter()
        .filter(|case| runs_here(program, case))
        .collect()
}

/// `rows`, the rows that differ, each with what LOSAT and NCBI wrote (the start of each stream),
/// for the message of a failed test.
fn describe(rows: &[(String, Streams, Streams)]) -> String {
    let start = |bytes: &[u8]| {
        let text = String::from_utf8_lossy(bytes);
        text.chars().take(600).collect::<String>()
    };
    let mut text = String::new();
    for (name, got, want) in rows {
        text.push_str(&format!(
            "\n--- {name}: LOSAT status {}, NCBI {}\nLOSAT stdout: {:?}\nNCBI stdout:  {:?}\n\
             LOSAT stderr: {:?}\nNCBI stderr:  {:?}\n",
            got.2,
            want.2,
            start(&got.0),
            start(&want.0),
            start(&got.1),
            start(&want.1),
        ));
    }
    text
}

/// The rows of `program`'s table whose run differs from NCBI's.
fn rows_that_differ(program: &str) -> Vec<String> {
    let mut differ = Vec::new();
    for case in cases_here(program) {
        let (got, want) = (run(program, &case), expected(program, &case));
        if got != want {
            differ.push((case.name.clone(), got, want));
        }
    }
    if !differ.is_empty() {
        eprintln!("{}", describe(&differ));
    }
    differ.into_iter().map(|(name, _, _)| name).collect()
}

#[test]
fn dustmasker_gives_ncbis_bytes() {
    let differ = rows_that_differ("dustmasker");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn segmasker_gives_ncbis_bytes() {
    let differ = rows_that_differ("segmasker");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn windowmasker_gives_ncbis_bytes() {
    let differ = rows_that_differ("windowmasker");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

/// The rows of `EMPTY_ID_ROWS` are what the table says: each row is NCBI's crash (exit status 139,
/// no output) after the error of its first such line, the one `bad id` error of a row of one line;
/// and each row is its near row with only its list changed, the list of the near row with the lines
/// of the rows whose errors it takes, in order, one line each.
#[test]
fn the_empty_id_rows_are_ncbis_crash_and_near_input() {
    let rows = cases("windowmasker");
    let row = |name: &str| rows.iter().find(|case| case.name == name).expect(name);
    let list = |case: &Case| {
        let at = case
            .words
            .iter()
            .position(|word| word == "-ids" || word == "-exclude_ids");
        let name = &case.words[at.expect("an id list") + 1];
        let text = fs::read_to_string(fixtures("windowmasker").join(name)).unwrap();
        text.lines().map(str::to_string).collect::<Vec<_>>()
    };
    // The lines of `case`'s list that are not those of `near`'s, which keeps them in order.
    let added = |case: &Case, near: &Case| {
        let (lines, kept) = (list(case), list(near));
        let mut next = kept.iter().peekable();
        let added: Vec<String> = lines
            .into_iter()
            .filter(|line| next.next_if(|kept| *kept == line).is_none())
            .collect();
        assert!(
            next.peek().is_none(),
            "{}: lines of {} missing",
            case.name,
            near.name
        );
        added
    };
    for (name, near, posts) in EMPTY_ID_ROWS {
        let (case, near) = (row(name), row(near));
        let changed: Vec<_> = case
            .words
            .iter()
            .zip(&near.words)
            .filter(|(a, b)| a != b)
            .collect();
        assert!(
            case.words.len() == near.words.len() && changed.len() == 1,
            "{name}"
        );
        let lines: Vec<String> = posts
            .iter()
            .flat_map(|&post| {
                let post = row(post);
                let &(_, near_of_post, _) = empty_id_row("windowmasker", &post.name).unwrap();
                added(post, row(near_of_post))
            })
            .collect();
        assert_eq!(added(case, near), lines, "{name}");
        assert_eq!(lines.len(), posts.len(), "{name}");
        // NCBI's own run: its crash after the error of the first such line.
        let (out, err, status) = expected_row("windowmasker", name);
        assert_eq!((out.len(), status), (0, 139), "{name}");
        assert_eq!(err, expected_row("windowmasker", posts[0]).1, "{name}");
        if *posts == [name] {
            let err = String::from_utf8(err).unwrap();
            assert_eq!(err.lines().count(), 1, "{name}");
            assert!(
                err.starts_with("Error: CWinMaskConfig::CIdSet_TextMatch::insert(): bad id: ")
                    && err.ends_with(": ignoring\n"),
                "{name}: {err}"
            );
        }
    }
}

#[test]
fn convert2blastmask_gives_ncbis_bytes() {
    let differ = rows_that_differ("convert2blastmask");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn blast_formatter_gives_ncbis_bytes() {
    let differ = rows_that_differ("blast_formatter");
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

/// blast_formatter's standard output on a full device and on a closed pipe: NCBI's runtime
/// aborts (`std::terminate` inside the serial stream's destructor, its post's flush of `cout`
/// through `cerr`'s tie: exit status 134 by SIGABRT, the runtime's report on the standard error,
/// no post; probes of line RG item F3, agent C6, 2026-11-23, NCBI BLAST+ 2.17.0, the official
/// Linux x64 release: every archive of the table with outfmt 8, 9 and 12).
#[cfg(target_os = "linux")]
#[test]
fn blast_formatter_standard_output_failure_aborts() {
    use std::os::unix::process::ExitStatusExt;
    let want: &[u8] = b"terminate called after throwing an instance of 'std::__ios_failure'\n  \
        what():  basic_ios::clear: iostream error\n";
    let dir = fixtures("blast_formatter");
    for (archive, outfmt) in [("db", "8"), ("bn", "12")] {
        let words = format!("{archive}.asn -outfmt {outfmt}");
        let case = Case {
            name: "standard_output_failure".to_string(),
            words: vec![
                "-archive".to_string(),
                dir.join("in")
                    .join(format!("{archive}.asn"))
                    .display()
                    .to_string(),
                "-outfmt".to_string(),
                outfmt.to_string(),
            ],
            stdin: None,
            // The run is out of the fixtures: the archive's database `in/db/s3` from `BLASTDB`.
            variables: vec![("BLASTDB".to_string(), dir.display().to_string())],
        };
        let binary = Path::new(env!("CARGO_BIN_EXE_losat"));
        let full = fs::OpenOptions::new()
            .write(true)
            .open("/dev/full")
            .unwrap();
        let (reader, writer) = std::io::pipe().unwrap();
        drop(reader);
        for (what, stdout) in [
            ("/dev/full", Stdio::from(full)),
            ("a closed pipe", Stdio::from(writer)),
        ] {
            // Out of the fixtures, where an abort may leave a core file.
            let output = command(binary, "blast_formatter", &case)
                .current_dir(env!("CARGO_TARGET_TMPDIR"))
                .stdout(stdout)
                .stderr(Stdio::piped())
                .output()
                .unwrap();
            assert_eq!(
                (output.status.signal(), output.stderr.as_slice()),
                (Some(6), want),
                "{words} on {what}: {:?}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }
}

// `run_dustmasker`, `run_segmasker`, `run_windowmasker`, `run_convert2blastmask` and
// `run_blast_formatter` run NCBI's application start as the command line does
// (`ncbi_app::run_with`, line RG item E1w-2): the rows of `-logfile`,
// `-conffile` and `-dryrun` give the command line's bytes through the library too. `-xmlhelp`,
// whose XML only the command line writes, is refused (line RG DECISIONS E2x-1). The example ends
// every failure with exit status 1. The rows of `program`'s table whose run through the library
// differs.
fn rows_that_differ_through_the_library(program: &str, refusal_name: &str) -> Vec<String> {
    let example = example_binary("library_command_line");
    let mut differ = Vec::new();
    let mut started = 0;
    let xmlhelp_refused: Streams = (
        Vec::new(),
        format!(
            "Error: the option -xmlhelp runs on the command line of \
             LOSAT's {refusal_name} only: the library API and the Web engine do not write \
             NCBI's XML description of the arguments\n"
        )
        .into_bytes(),
        1,
    );
    // `-xmlhelp` alone, also where the table has no row of it (blast_formatter's XML is not
    // NCBI's bytes: it lists no `-rid`).
    let xmlhelp = Case {
        name: "xmlhelp_through_the_library".to_string(),
        words: vec!["-xmlhelp".to_string()],
        stdin: None,
        variables: Vec::new(),
    };
    let got = run_binary(&example, program, &xmlhelp);
    if got != xmlhelp_refused {
        differ.push((xmlhelp.name.clone(), got, xmlhelp_refused.clone()));
    }
    for case in cases_here(program) {
        let (out, err, status) = run_binary(&example, program, &case);
        let help = case
            .words
            .iter()
            .find(|given| matches!(given.as_str(), "-h" | "-help" | "-help-full" | "-xmlhelp"));
        let expected = if help.is_some_and(|given| given == "-xmlhelp") {
            xmlhelp_refused.clone()
        } else {
            if ["-conffile", "-dryrun"].into_iter().any(|word| {
                case.words
                    .iter()
                    .any(|given| given == word || given.starts_with(&format!("{word}=")))
            }) {
                started += 1;
            }
            let (out, err, status) = expected(program, &case);
            (out, err, i32::from(status != 0))
        };
        let got = (out, err, status);
        if got != expected {
            differ.push((case.name.clone(), got, expected));
        }
    }
    assert!(started > 0, "{program}: no row of -conffile or -dryrun");
    if !differ.is_empty() {
        eprintln!("{}", describe(&differ));
    }
    differ.into_iter().map(|(name, _, _)| name).collect()
}

#[test]
fn dustmasker_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("dustmasker", "DUSTMASKER");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn segmasker_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("segmasker", "SEGMASKER");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn windowmasker_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("windowmasker", "WINDOWMASKER");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn convert2blastmask_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("convert2blastmask", "CONVERT2BLASTMASK");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

#[test]
fn blast_formatter_through_the_library_gives_the_command_line_s_bytes() {
    let differ = rows_that_differ_through_the_library("blast_formatter", "BLAST_FORMATTER");
    assert!(
        differ.is_empty(),
        "rows that differ through the library: {differ:?}"
    );
}

/// NCBI posts a byte that its serial stream replaces when it writes it, its `cerr` tied to
/// `cout`: with the standard output and error in one file (`2>&1`) the post follows the bytes of
/// the document that NCBI's serial buffer (`COStreamBuffer`) had written by then. The rows of
/// `program` with `ncbi/<name>.merged` (NCBI's `2>&1` stream of a run of its own) whose run
/// through `binary` into one file differs, and the number of such rows.
fn merged_rows_that_differ(binary: &Path, program: &str) -> (usize, Vec<String>) {
    let ncbi = fixtures(program).join("ncbi");
    let mut ran = 0;
    let mut differ = Vec::new();
    for case in cases_here(program) {
        let Ok(want) = fs::read(ncbi.join(format!("{}.merged", case.name))) else {
            continue;
        };
        // NCBI's stream in LOSAT's words (`losat::wording`, AGENTS.md rule 8).
        let want = losat::wording::losat_output_bytes(&want).into_owned();
        ran += 1;
        let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
            "masker-merged-{program}-{}-{}",
            case.name,
            std::process::id()
        ));
        let one = fs::File::create(&path).unwrap();
        let status = command(binary, program, &case)
            .stdout(one.try_clone().unwrap())
            .stderr(one)
            .status()
            .unwrap();
        let got = fs::read(&path).unwrap();
        let _ = fs::remove_file(&path);
        if got != want || status.code() != Some(0) {
            eprintln!(
                "--- {}: status {status}\nLOSAT: {:?}\nNCBI:  {:?}",
                case.name,
                String::from_utf8_lossy(&got),
                String::from_utf8_lossy(&want)
            );
            differ.push(case.name.clone());
        }
    }
    (ran, differ)
}

#[test]
fn convert2blastmask_posts_where_ncbi_does_in_one_file() {
    for binary in [
        PathBuf::from(env!("CARGO_BIN_EXE_losat")),
        example_binary("library_command_line"),
    ] {
        let (ran, differ) = merged_rows_that_differ(&binary, "convert2blastmask");
        assert!(ran >= 5, "{ran} rows with a merged stream");
        assert!(
            differ.is_empty(),
            "{}: rows whose merged stream differs from NCBI's: {differ:?}",
            binary.display()
        );
    }
}

/// dustmasker's `seqloc_*` and `maskinfo_*` documents of a sequence whose id NCBI's stream writes
/// with `#` in each of its 300 Seq-locs (`order_*`): each post among the document's bytes.
#[test]
fn dustmasker_posts_where_ncbi_does_in_one_file() {
    for binary in [
        PathBuf::from(env!("CARGO_BIN_EXE_losat")),
        example_binary("library_command_line"),
    ] {
        let (ran, differ) = merged_rows_that_differ(&binary, "dustmasker");
        assert!(ran >= 4, "{ran} rows with a merged stream");
        assert!(
            differ.is_empty(),
            "{}: rows whose merged stream differs from NCBI's: {differ:?}",
            binary.display()
        );
    }
}

/// windowmasker's `-mk_counts` writes the counts files of `tests/winmask/counts` that NCBI BLAST+
/// 2.17.0 wrote (`tests/winmask/README.md`, 2026-10-19 and 11-12): the units of 9 letters of
/// `tests/fasta/MG1655.fna` (4.6 Mb) in the optimized binary, ascii and optimized ascii formats,
/// and the units of 2 letters of `dm2L_60k.fa` (optimized binary, a cache bit array of 0 words).
/// These formats write the same bytes to the standard output and to a file.
#[test]
fn windowmasker_writes_ncbis_counts_files_of_the_winmask_fixtures() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    for (input, unit, sformat, counts) in [
        ("fasta/MG1655.fna", "9", "obinary", "mg1655_u9.obinary"),
        ("fasta/MG1655.fna", "9", "ascii", "mg1655_u9.ascii"),
        ("fasta/MG1655.fna", "9", "oascii", "mg1655_u9.oascii"),
        ("winmask/fasta/dm2L_60k.fa", "2", "obinary", "dm_u2.obinary"),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .current_dir(&tests)
            .args(["windowmasker", "-mk_counts", "-in", input, "-unit", unit])
            .args(["-sformat", sformat])
            .env("BLAST_USAGE_REPORT", "false")
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let want = fs::read(tests.join("winmask").join("counts").join(counts)).unwrap();
        assert_eq!(output.status.code(), Some(0), "{counts}");
        assert!(output.stdout == want, "{counts}: the bytes differ");
    }
}

/// windowmasker writes the FASTA reader's warnings between records (post level Warning) while
/// earlier records' output may still be in NCBI's stdout buffer: in one file (`2>&1`) each
/// warning falls where NCBI's does (`ncbi/<name>.merged`, make_fixtures.py's `MERGED`).
#[test]
fn windowmasker_posts_where_ncbi_does_in_one_file() {
    for binary in [
        PathBuf::from(env!("CARGO_BIN_EXE_losat")),
        example_binary("library_command_line"),
    ] {
        let (ran, differ) = merged_rows_that_differ(&binary, "windowmasker");
        assert!(ran >= 6, "{ran} rows with a merged stream");
        assert!(
            differ.is_empty(),
            "{}: rows whose merged stream differs from NCBI's: {differ:?}",
            binary.display()
        );
    }
}

// The rows that Windows answers as a missing file (`expected`) are the ones whose `-in` NCBI
// read as a directory on Linux, and their standard output on Linux is that of the missing file
// (a writer that wrote something before the failed read, or as the program ended, would write it
// on Windows too, where the row would need its own standard output).
#[test]
fn the_rows_that_read_a_directory_are_the_ones_ncbi_read_as_one() {
    for program in [
        "dustmasker",
        "segmasker",
        "windowmasker",
        "convert2blastmask",
    ] {
        for case in cases(program) {
            let (out, err, _) = expected_row(program, &case.name);
            let read_a_directory = String::from_utf8_lossy(&err)
                .contains("basic_filebuf::xsgetn error reading the file: Is a directory");
            assert_eq!(
                input_directory(program, &case).is_some(),
                read_a_directory,
                "{program} {}",
                case.name
            );
            if read_a_directory {
                let (missing, _) = missing_input(program, &case);
                let (missing_out, _, _) = expected_row(program, missing);
                assert_eq!(out, missing_out, "{program} {}: standard output", case.name);
            }
        }
    }
}

#[test]
fn the_shell_words_of_the_table_are_read_as_the_shell_reads_them() {
    let words = |line: &str| shell_words(line);
    assert_eq!(
        words("-in in/basic.fa -window \" 5\""),
        (
            vec![
                "-in".to_string(),
                "in/basic.fa".into(),
                "-window".into(),
                " 5".into()
            ],
            None
        )
    );
    assert_eq!(
        words("-outfmt \"\""),
        (vec!["-outfmt".to_string(), String::new()], None)
    );
    assert_eq!(
        words("-masking_options 'a\"b' \"c'd\" ''"),
        (
            vec![
                "-masking_options".to_string(),
                "a\"b".into(),
                "c'd".into(),
                String::new()
            ],
            None
        )
    );
    assert_eq!(
        words("-in - < in/basic.fa"),
        (
            vec!["-in".to_string(), "-".into()],
            Some("in/basic.fa".into())
        )
    );
}
