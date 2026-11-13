//! `-export_search_strategy FILE` of blastn, blastp, blastx, tblastn and tblastx (line RG, item F1),
//! and of psiblast with a FASTA query, rpsblast and rpstblastn (item F1b-1), against the bytes
//! that NCBI BLAST+ 2.17.0 `-num_threads 1` writes: the strategy file (one
//! `Blast4-request` per query batch, written before the batch's search), the standard output, the
//! standard error and the exit status, with `-subject` and `-db`, several batches (`BATCH_SIZE`;
//! from the second request on, `EffectiveSearchSpace` is the batch before's), `-searchsp`,
//! `-lcase_masking` (blastn's `LCaseMask` holds the lower-case runs only, not the DUST or
//! WindowMasker intervals), `-query_loc`, the database lists, `-num_descriptions`, the archive
//! (outfmt 11) of the same batches, the request of a batch whose search fails (a database without
//! its sequence or header file, also after a batch with no valid query; blastx's batch of a record
//! without letters; tblastx `-matrix BLOSUM99`, whose gap costs are 0 and 0), and the order of the
//! errors (the file is made after the argument parser and before the other options are read: an
//! error after it leaves the file empty, the parser's errors leave an existing file alone).
//! The first batch's request alone when NCBI fails to format that batch in BLAST XML (a matrix
//! name its XML data does not know). psiblast's request of each query (before its rounds, with the
//! search space of the query before; `Web_StepNumber`), rpsblast's and rpstblastn's (the profile
//! database's length and `.aux` costs in the requests after the first; `-mt_mode 1` at 4 and 8
//! threads writes NCBI's `-num_threads 1` file, DECISIONS F1b-2). The rows that LOSAT refuses
//! (DECISIONS 15, 21 and 26 of line RG; `-in_pssm` of tblastn and psiblast, item F1b-2) check the
//! refusal. blastp's query blocks and the `0` side of
//! the default-on switches that the strategy's paths cross write the same files.
//!
//! The cases are the rows of `tests/fixtures/strategy/cases.tsv` (its header says what the
//! columns mean); the frozen bytes are `tests/fixtures/strategy/<name>.ss`, `.out` and `.err`,
//! written by NCBI's program run from a copy of `LOSAT/` (the RG-F1 freeze script
//! `remaining-gaps-20261112/F/f1/freeze_strategy.py` of the line's records, NCBI BLAST+ 2.17.0).
//! Every row runs at 1, 4 and 8 threads. An env value `{rps:NAMES[-EXT]}` is a folder of the
//! run's own with the profile databases NAMES (`+`-separated) of `tests/db_search/rps` inflated,
//! without their files of extension EXT (`rps_databases`); NCBI's messages name it `{RPS}`.

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

mod reach_support;

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixtures() -> PathBuf {
    engine().join("tests/fixtures/strategy")
}

struct Case {
    name: String,
    program: String,
    env: Vec<(String, String)>,
    args: Vec<String>,
    check: String,
    rc: i32,
    pre: Option<String>,
    expect: String,
}

fn cases() -> Vec<Case> {
    let table = std::fs::read_to_string(fixtures().join("cases.tsv")).expect("cases.tsv");
    let mut cases = Vec::new();
    for line in table.lines() {
        if line.starts_with('#') || line.starts_with("name\t") || line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), 8, "cases.tsv row: {line}");
        let env = if cells[2] == "-" {
            Vec::new()
        } else {
            cells[2]
                .split(',')
                .map(|item| {
                    let (name, value) = item.split_once('=').expect("NAME=VALUE");
                    (name.to_string(), value.to_string())
                })
                .collect()
        };
        cases.push(Case {
            name: cells[0].to_string(),
            program: cells[1].to_string(),
            env,
            args: cells[3]
                .split(' ')
                .map(|arg| arg.replace("%20", " "))
                .collect(),
            check: cells[4].to_string(),
            rc: cells[5].parse().expect("rc"),
            pre: (cells[6] != "-").then(|| cells[6].to_string()),
            expect: cells[7].to_string(),
        });
    }
    cases
}

fn frozen(name: &str, ext: &str) -> Option<Vec<u8>> {
    std::fs::read(fixtures().join(format!("{name}.{ext}"))).ok()
}

/// One run of the case from `LOSAT/` in a clean environment plus the case's and `extra_env`;
/// `{ss}` is `ss`. Standard input is empty.
/// NCBI's messages were frozen on Linux. On Windows LOSAT writes the platform's separators, as
/// NCBI's Windows build does: `\` in paths and `;` between the directories of a search path
/// (`in search path [{ENGINE};;]` for Linux's `[{ENGINE}::]`).
fn unix_separators(text: &str) -> String {
    let text = text.replace('\\', "/");
    let mut out = String::with_capacity(text.len());
    let mut rest = text.as_str();
    const MARK: &str = "in search path [";
    while let Some(at) = rest.find(MARK) {
        let (head, tail) = rest.split_at(at + MARK.len());
        out.push_str(head);
        let end = tail.find(']').unwrap_or(tail.len());
        out.push_str(&tail[..end].replace(';', ":"));
        rest = &tail[end..];
    }
    out.push_str(rest);
    out
}

/// `{rps:NAMES[-EXT]}`: the profile databases NAMES (separated by `+`) of `tests/db_search/rps`,
/// their `NAME.*.gz` files inflated into `dir/rpsdb` without those of extension EXT (`-aux`: no
/// `.aux` file). The folder's path.
fn rps_databases(spec: &str, dir: &Path) -> PathBuf {
    use std::io::Read;
    let (names, drop) = match spec.split_once('-') {
        Some((names, extension)) => (names, Some(format!(".{extension}"))),
        None => (spec, None),
    };
    let target = dir.join("rpsdb");
    std::fs::create_dir_all(&target).unwrap();
    let source = engine().join("tests/db_search/rps");
    for entry in std::fs::read_dir(&source).unwrap() {
        let file = entry.unwrap().file_name().to_string_lossy().into_owned();
        let Some(stem) = file.strip_suffix(".gz") else {
            continue;
        };
        if !names
            .split('+')
            .any(|name| stem.starts_with(&format!("{name}.")))
            || drop
                .as_deref()
                .is_some_and(|extension| stem.ends_with(extension))
        {
            continue;
        }
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(std::fs::File::open(source.join(&file)).unwrap())
            .read_to_end(&mut bytes)
            .unwrap();
        std::fs::write(target.join(stem), bytes).unwrap();
    }
    target
}

/// The value of an env cell in the run's folder `dir` (`rps_databases`).
fn env_value(value: &str, dir: &Path) -> String {
    match value
        .strip_prefix("{rps:")
        .and_then(|spec| spec.strip_suffix('}'))
    {
        Some(spec) => rps_databases(spec, dir)
            .to_str()
            .expect("UTF-8 path")
            .to_string(),
        None => value.to_string(),
    }
}

fn run(case: &Case, ss: &Path, extra_env: &[(&str, &str)], threads: &str) -> Output {
    let dir = ss.parent().expect("the run's folder");
    let ss = ss.to_str().expect("UTF-8 path");
    let mut command = Command::new(BIN);
    command
        .current_dir(engine())
        .arg(&case.program)
        .args(case.args.iter().map(|arg| arg.replace("{ss}", ss)))
        .args(["-num_threads", threads])
        .stdin(Stdio::null());
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE" | "BLASTDB"
            )
        {
            command.env_remove(name);
        }
    }
    for (name, value) in &case.env {
        command.env(name, env_value(value, dir));
    }
    for (name, value) in extra_env {
        command.env(name, value);
    }
    command.output().expect("run losat")
}

/// A folder of this run's own (the tests of one binary run at the same time).
fn run_dir(case: &Case, label: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat_strategy_{}_{}_{label}",
        case.name,
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn check(case: &Case, extra_env: &[(&str, &str)], threads: &str, label: &str) {
    let dir = run_dir(case, label);
    let ss = dir.join("s.ss");
    if let Some(pre) = &case.pre {
        std::fs::write(&ss, pre).unwrap();
    }
    let output = run(case, &ss, extra_env, threads);
    let label = format!("{} threads {threads} {extra_env:?}", case.name);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(
        output.status.code(),
        Some(case.rc),
        "{label}: exit status; stderr: {stderr}"
    );
    match case.check.as_str() {
        "ncbi" | "parser" => {
            let written = std::fs::read(&ss).ok();
            let want = frozen(&case.name, "ss");
            assert_eq!(
                written.is_some(),
                want.is_some(),
                "{label}: the strategy file is there"
            );
            if let (Some(written), Some(want)) = (&written, &want) {
                assert!(
                    written == want,
                    "{label}: strategy file differs ({} bytes, NCBI wrote {}): {}",
                    written.len(),
                    want.len(),
                    crate::outfmt11_blastx::first_difference(written, want)
                );
            }
            if case.check == "parser" {
                assert!(output.stdout.is_empty(), "{label}: standard output");
                assert!(stderr.starts_with("error: "), "{label}: {stderr}");
            } else {
                let want_out = frozen(&case.name, "out").unwrap_or_default();
                assert!(
                    output.stdout == want_out,
                    "{label}: stdout differs ({} bytes, NCBI wrote {}): {}",
                    output.stdout.len(),
                    want_out.len(),
                    crate::outfmt11_blastx::first_difference(&output.stdout, &want_out)
                );
                let engine = engine();
                let rps = dir.join("rpsdb");
                let mut stderr = stderr;
                for form in [std::fs::canonicalize(&rps).unwrap_or_default(), rps] {
                    if !form.as_os_str().is_empty() {
                        stderr = stderr.replace(form.to_str().unwrap(), "{RPS}");
                    }
                }
                let stderr = stderr.replace(engine.to_str().unwrap(), "{ENGINE}");
                let stderr = if cfg!(windows) {
                    unix_separators(&stderr)
                } else {
                    stderr
                };
                let want_err = frozen(&case.name, "err").unwrap_or_default();
                assert_eq!(
                    stderr,
                    String::from_utf8_lossy(&want_err),
                    "{label}: stderr"
                );
            }
        }
        "refuse" => {
            assert!(output.stdout.is_empty(), "{label}: standard output");
            assert!(
                stderr.contains(&case.expect),
                "{label}: the refusal says {:?}: {stderr}",
                case.expect
            );
        }
        other => panic!("{label}: check {other}"),
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn strategy_files_equal_ncbis_at_1_4_and_8_threads() {
    let cases = cases();
    assert!(cases.len() >= 100, "cases.tsv rows");
    for case in &cases {
        for threads in ["1", "4", "8"] {
            check(case, &[], threads, threads);
        }
    }
}

/// blastp's query-block search (`LOSAT_X_BPQBLOCK`, a block of one batch) writes the same
/// strategy file: each batch's request is the same document in or out of a block.
#[test]
fn blastp_query_blocks_do_not_change_the_strategy_file() {
    let cases = cases();
    let multi = cases
        .iter()
        .filter(|case| case.name == "blastp.sub.batches" || case.name == "blastp.db.batches");
    for case in multi {
        for threads in ["1", "4"] {
            let reach = reach_support::Reach::new("strategy_bp");
            check(
                case,
                &[("LOSAT_X_BPQBLOCK", "1"), reach.var()],
                threads,
                "qblock1",
            );
            let sides = reach.sides("LOSAT_X_BPQBLOCK");
            assert!(
                sides.0,
                "{} threads {threads}: LOSAT_X_BPQBLOCK did not take the block search: {sides:?}",
                case.name
            );
            check(case, &[("LOSAT_X_BPQBLOCK", "0")], threads, "qblock0");
        }
    }
    // NCBI fails when it formats the first batch in BLAST XML: the file holds the first batch's
    // request alone, in or out of the blocks.
    let xml = cases
        .iter()
        .find(|case| case.name == "blastp.sub.xmlmatrix.batches")
        .expect("row blastp.sub.xmlmatrix.batches");
    for block in ["1", "0"] {
        check(
            xml,
            &[("LOSAT_X_BPQBLOCK", block)],
            "1",
            &format!("xmlqblock{block}"),
        );
    }
}

/// The default-on switches that the strategy's paths cross, at their `0` side (the paths they
/// replaced, which give the same bytes): blastx's flow of batches (`LOSAT_X_BXFLOW`), blastn's
/// streamed report and query reading (`LOSAT_X_RPTSTREAM`, `LOSAT_X_QSTREAM`, `LOSAT_X_BNQPACK`),
/// and the packed `-subject` stores of tblastn and tblastx (`LOSAT_X_TBNSTORE`,
/// `LOSAT_X_TBXSTORE`).
#[test]
fn switches_off_write_the_same_strategy_files() {
    let off: &[(&str, &str)] = &[
        ("blastx", "LOSAT_X_BXFLOW"),
        ("blastn", "LOSAT_X_RPTSTREAM"),
        ("blastn", "LOSAT_X_QSTREAM"),
        ("blastn", "LOSAT_X_BNQPACK"),
        ("tblastn", "LOSAT_X_TBNSTORE"),
        ("tblastx", "LOSAT_X_TBXSTORE"),
    ];
    let cases = cases();
    for &(program, switch) in off {
        let rows = cases
            .iter()
            .filter(|case| case.program == program && case.check == "ncbi");
        for case in rows {
            for threads in ["1", "4"] {
                check(
                    case,
                    &[(switch, "0")],
                    threads,
                    &format!("{switch}_{threads}"),
                );
            }
        }
    }
}

/// blastx refuses a query that is not a regular file (a pipe) with the strategy, as with outfmt 11
/// (DECISIONS 26 of line RG): the titles of the whole input are read before the first batch.
#[test]
fn blastx_refuses_a_piped_query() {
    use std::io::Write;
    let dir = std::env::temp_dir().join(format!("losat_strategy_pipe_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let ss = dir.join("s.ss");
    let query = std::fs::read(engine().join("tests/fasta/outfmt11/arch_bx_q.fna")).unwrap();
    let mut child = Command::new(BIN)
        .current_dir(engine())
        .args([
            "blastx",
            "-query",
            "-",
            "-subject",
            "tests/fasta/outfmt11/arch_bx_s.faa",
        ])
        .args([
            "-outfmt",
            "6",
            "-export_search_strategy",
            ss.to_str().unwrap(),
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run losat");
    // The refusal may come before the whole input is read.
    let _ = child.stdin.take().unwrap().write_all(&query);
    let output = child.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(output.stdout.is_empty());
    assert!(
        stderr.contains("unsupported BLASTX query input for -export_search_strategy"),
        "{stderr}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
