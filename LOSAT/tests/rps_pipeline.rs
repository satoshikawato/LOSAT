//! `losat rpsblast` and `losat rpstblastn` into `losat rpsbproc` (series PR, item PS15a-1): the
//! search writes NCBI BLAST+ 2.17.0's BLAST archive (`-outfmt 11`) and BLAST XML (`-outfmt 5`)
//! byte for byte (the frozen NCBI bytes of the rows `db.rps.of11.xojd`, `db.rps.of5.xojd`,
//! `db.rpstbn.of11.d` and `db.rpstbn.of5.d` of `tests/db_search_manifest/`, which the
//! manifest's runner also checks), and `losat rpsbproc` reads them, from a file and through a
//! pipe from the search, with the bytes of NCBI's rpsbproc 0.5.1 (`tests/fixtures/rps_pipeline/`: its
//! standard output and standard error for each input; the data files are `tests/fixtures/rpsbproc/data`, the profile
//! database `rps_cdd` of `tests/db_search/rps/`, whose profiles carry the PSSM ids of three rows of
//! CDD v3.22). NCBI is not run here.

use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::{Command, Output, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

mod temp_support;

fn engine() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

/// A directory of the test's own holding the inflated `rps_cdd` files and a copy of the CDD data
/// files, removed at the end (`temp_support::TempPath`).
struct RunDir(temp_support::TempPath);

impl RunDir {
    fn new(stem: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "losat_rps_pipeline_{stem}_{}_{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(dir.join("data")).unwrap();
        let rps = engine().join("tests/db_search/rps");
        for entry in fs::read_dir(&rps).unwrap() {
            let name = entry.unwrap().file_name().into_string().unwrap();
            if let Some(plain) = name
                .strip_prefix("rps_cdd.")
                .and_then(|rest| rest.strip_suffix(".gz"))
            {
                let mut bytes = Vec::new();
                flate2::read::GzDecoder::new(fs::File::open(rps.join(&name)).unwrap())
                    .read_to_end(&mut bytes)
                    .unwrap();
                fs::write(dir.join(format!("rps_cdd.{plain}")), bytes).unwrap();
            }
        }
        let data = engine().join("tests/fixtures/rpsbproc/data");
        for entry in fs::read_dir(&data).unwrap() {
            let entry = entry.unwrap();
            fs::copy(entry.path(), dir.join("data").join(entry.file_name())).unwrap();
        }
        RunDir(temp_support::TempPath::new(dir))
    }
}

/// `losat words...` in `dir`, without LOSAT or NCBI variables and registry files.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI coverage) needs that file name as soon as a search
/// has more than one thread, and it selects nothing.
fn losat(dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || text.starts_with("NCBI_")
            || matches!(text.as_ref(), "BLASTDB" | "BATCH_SIZE" | "ARCHIVE_FORMAT")
        {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .env("BLASTDB", dir)
        .current_dir(dir);
    command
}

fn assert_bytes(what: &str, got: &[u8], want: &[u8]) {
    if got != want {
        let at = got
            .iter()
            .zip(want)
            .position(|(a, b)| a != b)
            .unwrap_or(got.len().min(want.len()));
        panic!(
            "{what}: {} bytes, NCBI's {}; first difference at byte {at}: {:?} against {:?}",
            got.len(),
            want.len(),
            String::from_utf8_lossy(&got[at.saturating_sub(40)..(at + 40).min(got.len())]),
            String::from_utf8_lossy(&want[at.saturating_sub(40)..(at + 40).min(want.len())]),
        );
    }
}

/// NCBI's frozen bytes in LOSAT's words (`losat::wording`, AGENTS.md rule 8: the archive's
/// `ident`, the SESSION line, rpsbproc's banner).
fn frozen(path: &str) -> Vec<u8> {
    let ncbi = fs::read(engine().join("tests/fixtures").join(path)).unwrap();
    losat::wording::losat_output_bytes(&ncbi).into_owned()
}

/// The search of a case: the program, its query and options, the frozen archive and XML of
/// NCBI's run, and the stems of NCBI's rpsbproc outputs.
struct Case {
    program: &'static str,
    query: &'static str,
    options: &'static [&'static str],
    archive: &'static str,
    xml: &'static str,
    stem: &'static str,
}

const CASES: [Case; 2] = [
    Case {
        program: "rpsblast",
        query: "tests/fasta/rps/xoj.faa",
        options: &[],
        archive: "db_search/db.rps.of11.xojd.out",
        xml: "db_search/db.rps.of5.xojd.out",
        stem: "prot",
    },
    Case {
        program: "rpstblastn",
        query: "tests/fasta/rps/ps15n.fna",
        options: &["-evalue", "1e-3"],
        archive: "db_search/db.rpstbn.of11.d.out",
        xml: "db_search/db.rpstbn.of5.d.out",
        stem: "nucl",
    },
];

const REPORT: [&str; 7] = ["-d", "data", "-t", "both", "-m", "std", "-f"];

fn search(dir: &Path, case: &Case, outfmt: &str, threads: &str) -> Command {
    let mut command = losat(dir);
    command
        .arg(case.program)
        .args(["-db", "rps_cdd", "-query"])
        .arg(engine().join(case.query))
        .args(case.options)
        .args(["-outfmt", outfmt, "-num_threads", threads]);
    command
}

fn run(mut command: Command) -> Output {
    let output = command.stdin(Stdio::null()).output().unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{command:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn check_report(output: &Output, stem: &str, input: &str) {
    assert_bytes(
        &format!("rpsbproc {input}: standard output"),
        &output.stdout,
        &frozen(&format!("rps_pipeline/{stem}.out")),
    );
    assert_bytes(
        &format!("rpsbproc {input}: standard error"),
        &output.stderr,
        &frozen(&format!("rps_pipeline/{stem}.err")),
    );
}

/// The archive and the XML are NCBI's bytes at one, four and eight threads, and rpsbproc's
/// reports of them (read from a file) are NCBI's.
#[test]
fn the_archive_and_the_xml_are_ncbis_and_rpsbproc_reads_them() {
    let dir = RunDir::new("files");
    for case in &CASES {
        for (outfmt, file, frozen_path) in [
            ("11", format!("{}.asn", case.stem), case.archive),
            ("5", format!("{}.xml", case.stem), case.xml),
        ] {
            for threads in ["1", "4", "8"] {
                let mut command = search(&dir.0, case, outfmt, threads);
                command.args(["-out", &file]);
                run(command);
                assert_bytes(
                    &format!("{} -outfmt {outfmt} -num_threads {threads}", case.program),
                    &fs::read(dir.0.join(&file)).unwrap(),
                    &frozen(frozen_path),
                );
            }
            let mut command = losat(&dir.0);
            command.arg("rpsbproc").args(["-i", &file]).args(REPORT);
            check_report(
                &run(command),
                &format!("{}.{}", case.stem, &file[5..]),
                &file,
            );
        }
    }
}

/// `rpsblast -outfmt 11 | rpsbproc` and `rpstblastn -outfmt 11 | rpsbproc`: the search's standard
/// output is the report's standard input, as NCBI's `rpsbproc` reads an archive from a pipe.
#[test]
fn the_search_pipes_into_rpsbproc() {
    let dir = RunDir::new("pipe");
    for case in &CASES {
        let mut search = search(&dir.0, case, "11", "2")
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let mut report = losat(&dir.0);
        report.arg("rpsbproc").args(REPORT);
        let report = report
            .stdin(Stdio::from(search.stdout.take().unwrap()))
            .output()
            .unwrap();
        assert_eq!(search.wait().unwrap().code(), Some(0));
        assert_eq!(report.status.code(), Some(0));
        check_report(
            &report,
            &format!("{}.stdin", case.stem),
            "the standard input",
        );
    }
}
