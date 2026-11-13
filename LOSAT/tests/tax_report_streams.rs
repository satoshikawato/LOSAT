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
//! `T/t2/fix3/probe_tie.sh`, `blastx_tie.stdout.txt`).

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// The taxonomy database's text of B's lookup (`offset_last_past_end`).
const OFFSET_ERROR: &str = "Error: [blastp] Error: Offset error at end of taxdb file.\n";

/// A case: its name, whether the second query (else the first) has the title, the title's
/// length, and the `-out` file's room, exit status and standard error under `prlimit --fsize`.
struct Case {
    name: &'static str,
    second_titled: bool,
    title: usize,
    room: u64,
    fsize_exit: i32,
    fsize_stderr: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        name: "length",
        second_titled: true,
        title: 7441,
        room: 8195,
        fsize_exit: 6,
        fsize_stderr: &[OFFSET_ERROR, OFFSET_ERROR, "BLAST failed to write output\n"],
    },
    Case {
        name: "query",
        second_titled: false,
        title: 7447,
        room: 8192,
        fsize_exit: 0,
        fsize_stderr: &[OFFSET_ERROR, OFFSET_ERROR, OFFSET_ERROR],
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
    let mut command = Command::new(BIN);
    command
        .current_dir(manifest().join("tests/db_search"))
        .args([program, "-db", &format!("ids/{db}"), "-outfmt", "18"])
        .arg("-query")
        .arg(query)
        .args(["-num_threads", threads])
        .args(options);
    for (variable, _) in std::env::vars() {
        // `LOSAT_BLASTX_WORKER_LOG` stays: an `--all-features` binary writes its probe there.
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
    command
        .env("BLASTDB", format!("../seqdb/taxdb/{taxdb}"))
        .env("HOME", dir)
        .env("USERPROFILE", dir);
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
        let output = under_fsize(&report, case.room, Stdio::null());
        let written = std::fs::read(&out).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(
            (
                output.status.code(),
                String::from_utf8_lossy(&output.stderr)
            ),
            (Some(case.fsize_exit), case.fsize_stderr.concat().into()),
            "{}",
            case.name
        );
        // The file holds the report's first `room` bytes (the report: NCBI's `-out` stream without
        // the texts on standard error).
        let report = String::from_utf8(fixture(&format!("{}.out.txt", case.name)))
            .unwrap()
            .replace(OFFSET_ERROR, "");
        assert!(
            written == report.as_bytes()[..case.room as usize],
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
