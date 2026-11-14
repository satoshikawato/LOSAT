//! Where the organism report (outfmt 18) reaches its stream (line RG T2, review 3 D1): NCBI's
//! buffer of the stream writes at the end of an insertion of NCBI's (`out << value`), the
//! `-out` file's `basic_filebuf` at the insertion that reaches its 8,191st byte, `cout`'s glibc
//! buffer in whole blocks, so the insertions decide where the taxonomy database's texts on
//! standard error come between the report's bytes, and which write fails on a full file. The
//! `Query=` block's insertions are `"Query"`, `"= "`, each wrapped line and `"\n"`, `"\nLength="`,
//! the length and `"\n"`.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:729
//! ```c++
//!         out << label << "= ";
//! ```
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:743-744
//! ```c++
//!             out << "\nLength=";
//!             out << cbs.GetInst().GetLength() <<"\n";
//! ```
//!
//! The cases (NCBI BLAST+ 2.17.0, review 3's probe B, `tests/fixtures/tax_report_streams`): blastp
//! of `ids/txr_p` with the taxonomy database `offset_last_past_end`, whose lookup of B's taxid
//! writes the offset error after each `Length=` (two queries: `txr_q.faa`'s sequence twice, one
//! with a title of "abcdefg" words). `length`: q2's title puts the `"\nLength="` of q2 at byte
//! 8,190 of the report, so the `-out` file's buffer writes 8,198 bytes there, the error text
//! after `Length=`; `query`: q1's title puts q2's `"Query"` at byte 8,186, so the buffer writes
//! at its `"= "`. Under `prlimit --fsize` (SIGXFSZ ignored) the `-out` file takes 8,195 bytes
//! (`length`: the write of 8,198 fails, "BLAST failed to write output", exit 6) or 8,192
//! (`query`: the write at `"= "` is the end of a later block, exit 0). And a query without
//! warnings writes nothing to `cerr` (NCBI posts none), so it does not flush `cout` (fix 3's
//! `T/t2/fix3/probe_tie.sh`, `blastx_tie.stdout.txt`), and a blastx run that holds its output for
//! `-export_search_strategy` writes the same bytes (review 4 D1, `blastx_strategy.*`), also when it
//! ends in NCBI's error or a failed write (review 5 D1, `blastx_hold_error.*`). A query batch that
//! NCBI does not search stops NCBI's report after its `Length=` line, where LOSAT refuses (review 5
//! D2, `unsearched.*`).

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// The taxonomy database's text of B's lookup (`offset_last_past_end`).
#[cfg(target_os = "linux")]
const OFFSET_ERROR: &str = "Error: [blastp] Error: Offset error at end of taxdb file.\n";

/// A case: its name, whether the second query (else the first) has the title, the title's
/// length, and (Linux, `prlimit`) its `-out` file under `prlimit --fsize`.
struct Case {
    name: &'static str,
    second_titled: bool,
    title: usize,
    #[cfg(target_os = "linux")]
    fsize: Fsize,
}

/// The `-out` file of a case under `prlimit --fsize`: its room, and the run's exit status and
/// standard error.
#[cfg(target_os = "linux")]
struct Fsize {
    room: u64,
    exit: i32,
    stderr: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        name: "length",
        second_titled: true,
        title: 7441,
        #[cfg(target_os = "linux")]
        fsize: Fsize {
            room: 8195,
            exit: 6,
            stderr: &[OFFSET_ERROR, OFFSET_ERROR, "BLAST failed to write output\n"],
        },
    },
    Case {
        name: "query",
        second_titled: false,
        title: 7447,
        #[cfg(target_os = "linux")]
        fsize: Fsize {
            room: 8192,
            exit: 0,
            stderr: &[OFFSET_ERROR, OFFSET_ERROR, OFFSET_ERROR],
        },
    },
];

fn manifest() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

/// A fresh directory of the target's temporary space.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-tax-report-streams-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A query file in `dir` of two queries, q1 and q2, each the sequence of the record `record`
/// (under `tests/db_search/ids`), the first (`second_titled`: the second) with a title of
/// `title` bytes of "abcdefg" words.
fn titled_query(record: &str, title: usize, second_titled: bool, dir: &Path) -> PathBuf {
    let record =
        std::fs::read_to_string(manifest().join("tests/db_search/ids").join(record)).unwrap();
    let sequence = record.split_once('\n').unwrap().1.trim_end();
    let words = "abcdefg ".repeat(title / 8 + 1);
    let title = words[..title].trim_end();
    let (first, second) = match second_titled {
        true => (String::new(), format!(" {title}")),
        false => (format!(" {title}"), String::new()),
    };
    let path = dir.join("query.fa");
    std::fs::write(
        &path,
        format!(">q1{first}\n{sequence}\n>q2{second}\n{sequence}\n"),
    )
    .unwrap();
    path
}

/// The command of `program` at `threads` threads, in `tests/db_search` (which holds the
/// databases), with the home `dir` and none of the caller's variables that change a run.
fn losat(program: &str, dir: &Path, threads: &str) -> Command {
    let mut command = Command::new(BIN);
    command
        .current_dir(manifest().join("tests/db_search"))
        .args([program, "-num_threads", threads]);
    for (variable, _) in std::env::vars() {
        // `LOSAT_BLASTX_WORKER_LOG` stays (set below to a file of the case).
        if variable.starts_with("LOSAT_X_")
            || variable.starts_with("NCBI")
            || matches!(
                variable.as_str(),
                "BATCH_SIZE" | "BLASTDB" | "CHUNK_SIZE" | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(variable);
        }
    }
    // The probe of an `--all-features` binary goes to a file of the case: `prlimit --fsize`
    // limits every file the run writes, and the caller's log may already be past the room.
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_some() {
        command.env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"));
    }
    command.env("HOME", dir).env("USERPROFILE", dir);
    command
}

/// The command of `program`'s organism report of the database `ids/<db>` and the query file
/// `query`, with the taxonomy database `taxdb` (under `tests/seqdb/taxdb`), in `dir`, at
/// `threads` threads, with `options`.
fn search(
    program: &str,
    db: &str,
    query: &Path,
    taxdb: &str,
    dir: &Path,
    threads: &str,
    options: &[&str],
) -> Command {
    let mut command = losat(program, dir, threads);
    command
        .args(["-db", &format!("ids/{db}"), "-outfmt", "18"])
        .arg("-query")
        .arg(query)
        .args(options)
        .env("BLASTDB", format!("../seqdb/taxdb/{taxdb}"));
    command
}

/// The command of a case's run in `dir` at `threads` threads, with `options`.
fn command(case: &Case, dir: &Path, threads: &str, options: &[&str]) -> Command {
    let query = titled_query("txr_q.faa", case.title, case.second_titled, dir);
    search(
        "blastp",
        "txr_p",
        &query,
        "offset_last_past_end",
        dir,
        threads,
        options,
    )
}

/// `command` run where a file takes `room` bytes (`prlimit --fsize`, SIGXFSZ ignored, so the
/// write that passes it fails), with standard output to `stdout`; its standard error and exit
/// status.
#[cfg(target_os = "linux")]
fn under_fsize(command: &Command, room: u64, stdout: Stdio) -> std::process::Output {
    let mut limited = Command::new("sh");
    limited
        .args(["-c", "trap '' XFSZ; exec prlimit --fsize=\"$0\" \"$@\""])
        .arg(room.to_string())
        .arg(command.get_program())
        .args(command.get_args())
        .current_dir(command.get_current_dir().unwrap())
        .stdout(stdout);
    for (variable, value) in command.get_envs() {
        match value {
            Some(value) => limited.env(variable, value),
            None => limited.env_remove(variable),
        };
    }
    limited.output().expect("run LOSAT CLI under prlimit")
}

/// The bytes of standard output and standard error of `command`, in one pipe (as `2>&1 | cat`),
/// and its exit status.
fn merged(mut command: Command) -> (Vec<u8>, Option<i32>) {
    use std::io::Read;
    let (mut reader, writer) = std::io::pipe().unwrap();
    command
        .stdout(Stdio::from(writer.try_clone().unwrap()))
        .stderr(Stdio::from(writer));
    let mut child = command.spawn().expect("run LOSAT CLI");
    // The parent's ends of the pipe go with the command.
    drop(command);
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes).unwrap();
    (bytes, child.wait().unwrap().code())
}

fn fixture(name: &str) -> Vec<u8> {
    std::fs::read(
        manifest()
            .join("tests/fixtures/tax_report_streams")
            .join(name),
    )
    .unwrap()
}

/// `-out /dev/stdout` with standard error in the same pipe: the `-out` file's buffer is not
/// written before a write to standard error (`cerr` is tied to `cout` only).
#[cfg(unix)]
#[test]
fn the_out_files_buffer_writes_at_the_end_of_ncbis_insertions() {
    for case in CASES {
        for threads in ["1", "4"] {
            let dir = scratch(&format!("{}-out-{threads}", case.name));
            let (bytes, exit) = merged(command(case, &dir, threads, &["-out", "/dev/stdout"]));
            std::fs::remove_dir_all(&dir).unwrap();
            assert_eq!(exit, Some(0), "{}", case.name);
            assert!(
                bytes == fixture(&format!("{}.out.txt", case.name)),
                "{} at {threads} threads:\n{}",
                case.name,
                String::from_utf8_lossy(&bytes)
            );
        }
    }
}

/// Standard output with standard error in the same pipe: the texts come after every byte inserted
/// before them (the tie).
#[test]
fn the_standard_outputs_texts_come_after_the_report_before_them() {
    for case in CASES {
        for threads in ["1", "4"] {
            let dir = scratch(&format!("{}-stdout-{threads}", case.name));
            let (bytes, exit) = merged(command(case, &dir, threads, &[]));
            std::fs::remove_dir_all(&dir).unwrap();
            assert_eq!(exit, Some(0), "{}", case.name);
            assert!(
                bytes == fixture(&format!("{}.stdout.txt", case.name)),
                "{} at {threads} threads:\n{}",
                case.name,
                String::from_utf8_lossy(&bytes)
            );
        }
    }
}

/// `-out F` that takes `room` bytes (`prlimit --fsize`, SIGXFSZ ignored, so the write that
/// passes it fails): the exit status says whether NCBI's buffer wrote past it during the run (6)
/// or only at the end (0, the failure unreported).
#[cfg(target_os = "linux")]
#[test]
fn a_full_out_file_fails_where_ncbis_buffer_writes() {
    for case in CASES {
        let dir = scratch(&format!("{}-fsize", case.name));
        let out = dir.join("report.out");
        let report = command(case, &dir, "1", &["-out", out.to_str().unwrap()]);
        let output = under_fsize(&report, case.fsize.room, Stdio::null());
        let written = std::fs::read(&out).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            (
                output.status.code(),
                String::from_utf8_lossy(&output.stderr)
            ),
            (Some(case.fsize.exit), case.fsize.stderr.concat().into()),
            "{}",
            case.name
        );
        // The file holds the report's first `room` bytes (the report: NCBI's `-out` stream without
        // the texts on standard error).
        let report = String::from_utf8(fixture(&format!("{}.out.txt", case.name)))
            .unwrap()
            .replace(OFFSET_ERROR, "");
        assert!(
            written == report.as_bytes()[..case.fsize.room as usize],
            "{}",
            case.name
        );
    }
}

/// A query without warnings writes nothing to standard error (NCBI posts none), so standard
/// output keeps what glibc's buffer holds: NCBI BLAST+ 2.17.0 (fix 3's `probe_tie.sh`), blastx of
/// two queries (q1's report 4,682 bytes) with the taxonomy database `good`, standard output a file
/// that takes 4,200 bytes: glibc writes the first block of 4,096 bytes while it writes q1's report,
/// the rest at the exit, whose failure is unreported (exit 0).
#[cfg(target_os = "linux")]
#[test]
fn a_query_without_warnings_leaves_standard_output_held() {
    let dir = scratch("blastx-tie");
    let query = titled_query("txr_q.fna", 3800, false, &dir);
    let out = dir.join("stdout");
    let run = search("blastx", "txr_p", &query, "good", &dir, "1", &[]);
    let output = under_fsize(
        &run,
        4200,
        Stdio::from(std::fs::File::create(&out).unwrap()),
    );
    let written = std::fs::read(&out).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(
        (
            output.status.code(),
            String::from_utf8_lossy(&output.stderr)
        ),
        (Some(0), "".into())
    );
    assert!(written == fixture("blastx_tie.stdout.txt")[..4200]);
}

/// blastx with `-export_search_strategy`, whose run holds its output while a later batch's title
/// could still be refused (`archive::ArchiveHold`): NCBI's buffer comes before the hold and the
/// texts go through it, so the merged stream is the run's without the option. NCBI BLAST+ 2.17.0
/// (line RG T2 review 4, `T/t2/review4/probe_strategy.sh`; fix 4, `T/t2/fix4/
/// probe_strategy_stdout.sh`), the `length` case of blastx (q2's title of 7,441 bytes): with
/// `-out /dev/stdout` the first three texts come before the report, which the file's buffer still
/// holds, and the fourth inside q2's `Length=`; on standard output each comes after the report
/// written before it (the tie). 9,002 bytes, exit 0 both.
#[cfg(unix)]
#[test]
fn a_held_blastx_report_is_written_where_ncbis_buffer_writes_it() {
    for (name, options) in [
        ("blastx_strategy.out.txt", &["-out", "/dev/stdout"][..]),
        ("blastx_strategy.stdout.txt", &[][..]),
    ] {
        for threads in ["1", "4"] {
            let dir = scratch(&format!("blastx-strategy-{}-{threads}", options.len()));
            let query = titled_query("txr_q.fna", 7441, true, &dir);
            let strategy = dir.join("search.strategy");
            let mut run = search(
                "blastx",
                "txr_p",
                &query,
                "offset_last_past_end",
                &dir,
                threads,
                options,
            );
            run.arg("-export_search_strategy").arg(&strategy);
            let (bytes, exit) = merged(run);
            std::fs::remove_dir_all(&dir).unwrap();
            assert_eq!(exit, Some(0), "{name}");
            assert!(
                bytes == fixture(name),
                "{name} at {threads} threads:\n{}",
                String::from_utf8_lossy(&bytes)
            );
        }
    }
}

/// The query file of review 5's probe (`T/t2/review5/probe_hold2.sh`, `in2/q_hold.fna`): q1 encodes
/// A (the first volume of `vol2_nophr`), q2 encodes B (the second volume, whose header file is
/// missing) and has a hyphen, q3's title is `q3 café` (`held`), which the search strategy cannot
/// write, so the run holds its output from the start, or `q3 cafe`; written into `dir`.
#[cfg(unix)]
fn hold_error_query(dir: &Path, held: bool) -> PathBuf {
    let query = String::from_utf8(fixture("blastx_hold_error.fna")).unwrap();
    let path = dir.join("query.fna");
    std::fs::write(
        &path,
        match held {
            true => query,
            false => query.replace("q3 café", "q3 cafe"),
        },
    )
    .unwrap();
    path
}

/// Review 5's probe: blastx `-outfmt 18 -max_target_seqs 1 -export_search_strategy` of
/// `hold_error_query` against `vol2_nophr` with `BATCH_SIZE=1`, in `dir`, at `threads` threads.
#[cfg(unix)]
fn hold_error_run(dir: &Path, threads: &str, held: bool) -> Command {
    let query = hold_error_query(dir, held);
    let mut run = losat("blastx", dir, threads);
    run.args([
        "-db",
        "vol2_nophr",
        "-outfmt",
        "18",
        "-max_target_seqs",
        "1",
    ])
    .arg("-query")
    .arg(query)
    .arg("-export_search_strategy")
    .arg(dir.join("search.strategy"))
    .env("BATCH_SIZE", "1");
    run
}

/// A blastx run that holds its output for `-export_search_strategy` (a later title the strategy
/// cannot write) and ends in one of NCBI's errors: NCBI BLAST+ 2.17.0 (line RG T2 review 5, D1,
/// `T/t2/review5/probe_hold2.sh`, `out2/ncbi.txt`, 603 bytes, exit status 2) writes q1's report,
/// the reader's text of batch 2 after it (`cerr`'s tie) and the database error. The hold releases
/// the report and the texts with the run's flushes, in their order (`ArchiveHold::release`), so
/// the held run (`q3 café`) writes the bytes of the run without a hold (`q3 cafe`) and NCBI's.
#[cfg(unix)]
#[test]
fn a_held_blastx_run_that_ends_in_ncbis_error_writes_ncbis_order() {
    let db = std::fs::canonicalize(manifest().join("tests/db_search")).unwrap();
    let expected = String::from_utf8(fixture("blastx_hold_error.txt"))
        .unwrap()
        .replace("{DB}", db.to_str().unwrap());
    for held in [true, false] {
        for threads in ["1", "4"] {
            let dir = scratch(&format!("blastx-hold-error-{held}-{threads}"));
            let (bytes, exit) = merged(hold_error_run(&dir, threads, held));
            std::fs::remove_dir_all(&dir).unwrap();
            assert_eq!(exit, Some(2), "held {held} at {threads} threads");
            assert!(
                bytes == expected.as_bytes(),
                "held {held} at {threads} threads:\n{}",
                String::from_utf8_lossy(&bytes)
            );
        }
    }
}

/// The run of `a_held_blastx_run_that_ends_in_ncbis_error_writes_ncbis_order` with standard output
/// a full device: NCBI BLAST+ 2.17.0 aborts where `cerr`'s tie writes q1's report (exit status 134,
/// `terminate called after throwing an instance of 'std::__ios_failure'`, review 5's `out2/
/// ncbi.full.err`), after the warning of `-max_target_seqs 1`. LOSAT reports the failed write with a
/// non-zero exit status where NCBI aborts (as PD-LOSAT-CLI-NONSEARCH-DIFFERENCES's outfmt 6 and 7;
/// `cli::NcbiReportBuffer`) and writes nothing after it: the run without a hold stops at the tie
/// (`Error: No space left on device`, exit status 1), the held run at the release's write of the
/// same bytes ("BLAST failed to write output", exit status 6; that the two texts differ is item T6
/// of line RG).
#[cfg(target_os = "linux")]
#[test]
fn a_held_blastx_run_on_a_full_standard_output_reports_the_failed_write() {
    const WARNING: &str = "Warning: [blastx] Examining 5 or more matches is recommended\n";
    for (held, exit, failure) in [
        (true, 6, "BLAST failed to write output\n"),
        (false, 1, "Error: No space left on device (os error 28)\n"),
    ] {
        let dir = scratch(&format!("blastx-hold-full-{held}"));
        let output = hold_error_run(&dir, "1", held)
            .stdout(
                std::fs::OpenOptions::new()
                    .write(true)
                    .open("/dev/full")
                    .unwrap(),
            )
            .output()
            .expect("run LOSAT CLI");
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            (
                output.status.code(),
                String::from_utf8_lossy(&output.stderr)
            ),
            (Some(exit), format!("{WARNING}{failure}").into()),
            "held {held}"
        );
    }
}

/// LOSAT's refusal of a query whose batch is not searched (`TaxReport::write_query`) in `program`.
#[cfg(unix)]
fn unsearched_refusal(program: &str) -> String {
    format!(
        "Error: the organism report (outfmt 18) of a query that is not searched (no query of its batch has a valid context: \"Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options\"), where NCBI BLAST+ stops with a null pointer exception; this is not supported by LOSAT's {program}\n"
    )
}

/// A query batch that NCBI does not search (no query of it has a valid context: 60 N, or 60 X for
/// a protein query) has a null Seq-align-set, which NCBI's organism report dereferences after the
/// query's `Length=` line (`x_PrintTaxReport`): NCBI BLAST+ 2.17.0 then writes its null pointer
/// exception with a stack trace of its build and exits with status 255 (line RG T2 fix 5,
/// `T/t2/fix5/probe_d2.sh`, n1; the fixtures `unsearched.*` hold NCBI's bytes before the
/// exception). LOSAT writes those bytes and refuses (exit status 1), in all five programs: with
/// `-subject`, and blastx with `-db` after a searched batch (`BATCH_SIZE=1`). A query without a
/// valid context in a batch that NCBI searches (`unsearched.*mix.txt`, NCBI's whole output, exit
/// status 0) gets the headers of a query without hits.
#[cfg(unix)]
#[test]
fn a_query_batch_that_ncbi_does_not_search_is_refused_after_its_length_line() {
    let ids = manifest().join("tests/db_search/ids");
    let invalid = |dir: &Path, letter: char| {
        let path = dir.join(format!("{letter}.fa"));
        let name = match letter {
            'N' => "xn all N",
            _ => "xp all X",
        };
        std::fs::write(
            &path,
            format!(">{name}\n{}\n", letter.to_string().repeat(60)),
        )
        .unwrap();
        path
    };
    // q1 (`ids/txr_q.fna`, which blastx's search finds in `ids/txr_p`), then 60 N.
    let searched_then_invalid = |dir: &Path| {
        let record = std::fs::read_to_string(ids.join("txr_q.fna")).unwrap();
        let path = dir.join("q_bx.fna");
        std::fs::write(
            &path,
            format!(
                "{}>xn all N\n{}\n",
                record.replacen(">qn\n", ">q1 searched\n", 1),
                "N".repeat(60)
            ),
        )
        .unwrap();
        path
    };
    for (name, program, refused) in [
        ("bn", "blastn", true),
        ("bp", "blastp", true),
        ("tbn", "tblastn", true),
        ("tbx", "tblastx", true),
        ("bx", "blastx", true),
        ("bnmix", "blastn", false),
        ("bxmix", "blastx", false),
    ] {
        for threads in ["1", "4"] {
            let dir = scratch(&format!("unsearched-{name}-{threads}"));
            let (query, subject) = match name {
                "bn" | "tbx" => (invalid(&dir, 'N'), Some("txr_q.fna")),
                "bp" => (invalid(&dir, 'X'), Some("txr_q.faa")),
                "tbn" => (invalid(&dir, 'X'), Some("txr_q.fna")),
                "bnmix" => (searched_then_invalid(&dir), Some("txr_q.fna")),
                _ => (searched_then_invalid(&dir), None),
            };
            let mut run = losat(program, &dir, threads);
            run.args(["-outfmt", "18", "-query"])
                .arg(&query)
                .env("BLASTDB", "../seqdb/taxdb/good");
            match subject {
                Some(subject) => run.arg("-subject").arg(ids.join(subject)),
                None => run.args(["-db", "ids/txr_p"]),
            };
            if name == "bx" {
                run.env("BATCH_SIZE", "1");
            }
            let (bytes, exit) = merged(run);
            std::fs::remove_dir_all(&dir).unwrap();
            let mut expected = fixture(&format!("unsearched.{name}.txt"));
            if refused {
                expected.extend_from_slice(unsearched_refusal(&program.to_uppercase()).as_bytes());
            }
            assert_eq!(exit, Some(if refused { 1 } else { 0 }), "{name}");
            assert!(
                bytes == expected,
                "{name} at {threads} threads:\n{}",
                String::from_utf8_lossy(&bytes)
            );
        }
    }
}
