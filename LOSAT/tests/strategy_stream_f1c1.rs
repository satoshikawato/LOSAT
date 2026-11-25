//! Integration tests for `-export_search_strategy` stream modeling, write-failure handling,
//! and shared/aliased non-`-` file sinks (`ostream::StrategyStream`).
//!
//! Covers:
//! - Failing strategy sinks (`/dev/full`): single-batch 4-line `Close` failure chain + exit 0,
//!   multi-batch `802.5 ... write fault` + exit 0 (stopping before later batches), large requests
//!   exceeding NCBI's 8,191-byte `CNcbiOfstream` buffer failing at `pos + size + 80`, and
//!   `-outfmt 11` archive message capture (4 `ERROR` messages on `Close`, error archive on
//!   `write fault`) across `blastn`, `blastp`, `blastx`, `tblastn`, and `tblastx`.
//! - Shared/aliased non-`-` sinks (`same`, `./F`, symlink, `/dev/stdout` file/pipe, `-out /dev/stdout`,
//!   `>>` append, `2>&1` merged pipe with `803.7` interleaving, and `/dev/stdout -> /dev/full`).
//! - `rpstblastn` `HeldOutput` interleaving for `-outfmt 5`, `15`, and `16`.
//! - Preservation of `-export_search_strategy -` refusal.

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

const BP_BATCHES_Q: &str = "tests/fasta/outfmt11/arch_bp_batches.faa";
const BP_S: &str = "tests/fasta/outfmt11/arch_bp_s.faa";

const BATCHES: [(&str, &str, &str, &str); 5] = [
    (
        "blastn",
        "tests/fasta/outfmt11/arch_batches.fna",
        "tests/fasta/outfmt11/arch_s.fna",
        "500",
    ),
    ("blastp", BP_BATCHES_Q, BP_S, "150"),
    (
        "blastx",
        "tests/fasta/outfmt11/arch_bx_batches.fna",
        "tests/fasta/outfmt11/arch_bx_s.faa",
        "200",
    ),
    (
        "tblastn",
        "tests/fasta/outfmt11/arch_tbn_batches.faa",
        "tests/fasta/outfmt11/arch_tbn_s.fna",
        "150",
    ),
    (
        "tblastx",
        "tests/fasta/outfmt11/arch_tbx_batches.fna",
        "tests/fasta/outfmt11/arch_tbx_s.fna",
        "200",
    ),
];

fn stderr_text(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stderr_lines(stderr: &str) -> Vec<&str> {
    stderr.lines().collect()
}

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn clean_command(program: &str) -> Command {
    let mut cmd = Command::new(BIN);
    cmd.current_dir(engine())
        .arg(program)
        .stdin(Stdio::null())
        .env("BLAST_USAGE_REPORT", "false")
        .env("LOSAT_BLASTX_WORKER_LOG", "/dev/null");
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE" | "BLASTDB"
            )
        {
            cmd.env_remove(name);
        }
    }
    cmd
}

fn temp_run_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat_f1c1_{tag}_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Shared file sinks (`-out F -export_search_strategy F` and `./F`) produce identical bytes at
/// 1, 4, and 8 threads and match the composition of separate strategy + report descriptors.
#[test]
fn same_file_and_dot_alias_match_across_threads() {
    let dir = temp_run_dir("same_dot");
    for (prog, q, s, bsz) in BATCHES {
        let mut baseline: Option<Vec<u8>> = None;
        for threads in ["1", "4", "8"] {
            let f = dir.join(format!("{prog}_t{threads}.out"));
            let out = clean_command(prog)
                .env("BATCH_SIZE", bsz)
                .args(["-query", q, "-subject", s, "-outfmt", "6"])
                .args(["-num_threads", threads])
                .arg("-out")
                .arg(&f)
                .arg("-export_search_strategy")
                .arg(&f)
                .output()
                .expect("run losat");
            assert_eq!(
                out.status.code(),
                Some(0),
                "{prog} t={threads} stderr={}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert!(out.stdout.is_empty());
            assert!(out.stderr.is_empty());
            let bytes = std::fs::read(&f).unwrap();
            assert!(!bytes.is_empty());
            match &baseline {
                None => baseline = Some(bytes),
                Some(want) => assert_eq!(&bytes, want, "{prog} t={threads} differed from t=1"),
            }
        }
        // `./` relative alias to the same `-out` file matches exact bytes.
        let dot_f = dir.join(format!("{prog}_dot.out"));
        let dot_alias = dir.join(".").join(format!("{prog}_dot.out"));
        let out = clean_command(prog)
            .env("BATCH_SIZE", bsz)
            .args([
                "-query",
                q,
                "-subject",
                s,
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ])
            .arg("-out")
            .arg(&dot_f)
            .arg("-export_search_strategy")
            .arg(&dot_alias)
            .output()
            .expect("run losat dot");
        assert_eq!(out.status.code(), Some(0));
        assert_eq!(
            std::fs::read(&dot_f).unwrap(),
            baseline.as_deref().unwrap(),
            "{prog} ./ alias differed"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `-export_search_strategy -` remains refused cleanly at NCBI's first request point.
#[test]
fn dash_standard_output_remains_refused() {
    let out = clean_command("blastp")
        .args([
            "-query",
            "tests/fasta/outfmt11/arch_bp_q.faa",
            "-subject",
            "tests/fasta/outfmt11/arch_bp_s.faa",
            "-outfmt",
            "6",
            "-export_search_strategy",
            "-",
            "-num_threads",
            "1",
        ])
        .output()
        .expect("run losat");
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stdout.is_empty());
    let stderr = stderr_text(&out);
    assert!(
        stderr.contains("-export_search_strategy - (the search strategy on the standard output"),
        "unexpected stderr: {stderr}"
    );
}

/// `/dev/full` strategy sink: single-batch `Close` 4-line chain (exit 0), multi-batch `write fault`
/// stopping after batch 0 (exit 0), and `-outfmt 11` error archive across all five programs.
#[cfg(unix)]
#[test]
fn dev_full_strategy_sink_matches_ncbi_close_chain_and_write_fault() {
    for (prog, q, s, bsz) in BATCHES {
        // Multi-batch outfmt 6 on `/dev/full`: 4-line Close failure on batch 0, then 1-line
        // `write fault` on batch 1, exit 0, and wording hides C++ exception class tags.
        let mut baseline: Option<(Vec<u8>, Vec<u8>)> = None;
        for threads in ["1", "4"] {
            let out = clean_command(prog)
                .env("BATCH_SIZE", bsz)
                .args([
                    "-query",
                    q,
                    "-subject",
                    s,
                    "-outfmt",
                    "6",
                    "-export_search_strategy",
                    "/dev/full",
                    "-num_threads",
                    threads,
                ])
                .output()
                .expect("run losat /dev/full f6");
            assert_eq!(out.status.code(), Some(0), "{prog} t={threads} rc");
            let stderr = stderr_text(&out);
            let lines = stderr_lines(&stderr);
            assert_eq!(lines.len(), 5, "{prog} t={threads} stderr:\n{stderr}");
            assert!(lines[0].contains("cannot close output stream"), "{stderr}");
            assert!(
                lines[1].contains("COStreamBuffer::Flush: failed"),
                "{stderr}"
            );
            assert!(lines[2].contains("cannot close output stream"), "{stderr}");
            assert!(
                lines[3].contains("Cannot close serializing output stream"),
                "{stderr}"
            );
            assert!(lines[4].contains("stack is empty: write fault"), "{stderr}");
            assert!(!stderr.contains("CIOException"), "{stderr}");
            assert!(!stderr.contains("CSerialException"), "{stderr}");
            match &baseline {
                None => baseline = Some((out.stdout, out.stderr)),
                Some(want) => assert_eq!(&(out.stdout, out.stderr), want, "{prog} t={threads}"),
            }
        }
        // Multi-batch outfmt 11 on `/dev/full`: batch 0 archive carries 4 `code 2` messages;
        // batch 1 write fault writes the error archive (`request-id "Error"`) with 1 `code 2`
        // `write fault` message and exits 0.
        let out11 = clean_command(prog)
            .env("BATCH_SIZE", bsz)
            .args([
                "-query",
                q,
                "-subject",
                s,
                "-outfmt",
                "11",
                "-export_search_strategy",
                "/dev/full",
                "-num_threads",
                "1",
            ])
            .output()
            .expect("run losat /dev/full f11");
        assert_eq!(out11.status.code(), Some(0), "{prog} f11 rc");
        let stdout11 = String::from_utf8_lossy(&out11.stdout);
        let docs: Vec<&str> = stdout11.split("Blast4-archive ::= ").skip(1).collect();
        assert_eq!(
            docs.len(),
            2,
            "{prog} f11 expected 1 batch archive + 1 error archive, got {}",
            docs.len()
        );
        assert_eq!(
            docs[0].matches("code 2,").count(),
            4,
            "{prog} f11 batch 0 archive should carry 4 Close messages:\n{}",
            docs[0]
        );
        assert!(
            docs[1].contains("request-id \"Error\""),
            "{prog} f11 doc 1 should be error archive:\n{}",
            docs[1]
        );
        assert_eq!(
            docs[1].matches("code 2,").count(),
            1,
            "{prog} f11 error archive should carry 1 write-fault message:\n{}",
            docs[1]
        );
        assert!(
            docs[1].contains("write fault"),
            "{prog} f11 error archive missing write fault:\n{}",
            docs[1]
        );
    }
}

/// A single request larger than NCBI's 8,191-byte `CNcbiOfstream` buffer on `/dev/full` triggers
/// `write fault` inside piece 2 at `pos + size + 80` (no 4-line `Close` chain, no report, exit 0).
#[cfg(unix)]
#[test]
fn dev_full_large_strategy_request_faults_at_piece_plus_back_limit() {
    let dir = temp_run_dir("large_full");
    let big_q = dir.join("big.faa");
    // 12,000-residue protein guarantees ASN.1 text `Blast4-request` > 12 KiB (> 8,191 B buffer).
    let mut fasta = b">big_query\n".to_vec();
    for _ in 0..200 {
        fasta.extend_from_slice(b"MKVLAIGSTVPFWDEYQRNHCKVLAIGSTVPFWDEYQRNHCKVLAIGSTVPFWDEYQRNH\n");
    }
    std::fs::write(&big_q, &fasta).unwrap();
    let out = clean_command("blastp")
        .arg("-query")
        .arg(&big_q)
        .args([
            "-subject",
            "tests/fasta/outfmt11/arch_bp_s.faa",
            "-outfmt",
            "6",
            "-export_search_strategy",
            "/dev/full",
            "-num_threads",
            "1",
        ])
        .output()
        .expect("run losat large /dev/full");
    assert_eq!(out.status.code(), Some(0));
    assert!(
        out.stdout.is_empty(),
        "report should not run after piece write fault"
    );
    let stderr = stderr_text(&out);
    let lines = stderr_lines(&stderr);
    assert_eq!(
        lines.len(),
        1,
        "expected single write-fault line:\n{stderr}"
    );
    assert!(
        lines[0].starts_with("Error: (802.5) [blastp] CObjectOStream: error at byte ")
            && lines[0].ends_with(": stack is empty: write fault"),
        "unexpected line: {}",
        lines[0]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// Symlink alias (`-out F -export_search_strategy L` with `L -> F`), `/dev/stdout` pipe (`dsp`),
/// `-out /dev/stdout -export_search_strategy /dev/stdout` (`bothp`), and `>>` append (`dsoapp`).
#[cfg(unix)]
#[test]
fn symlink_dev_stdout_and_append_sinks_compose_in_ncbi_order() {
    let dir = temp_run_dir("unix_sinks");
    let target = dir.join("F");
    let link = dir.join("L");
    std::fs::write(&target, b"").unwrap();
    std::os::unix::fs::symlink(&target, &link).unwrap();

    // 1. Symlink matches direct same-file bytes.
    let same_f = dir.join("same_F");
    let run_same = |out_p: &Path, ss_p: &Path| -> Output {
        clean_command("blastp")
            .env("BATCH_SIZE", "150")
            .args([
                "-query",
                BP_BATCHES_Q,
                "-subject",
                BP_S,
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ])
            .arg("-out")
            .arg(out_p)
            .arg("-export_search_strategy")
            .arg(ss_p)
            .output()
            .expect("run losat")
    };
    let o1 = run_same(&same_f, &same_f);
    let o2 = run_same(&target, &link);
    assert_eq!(o1.status.code(), Some(0));
    assert_eq!(o2.status.code(), Some(0));
    assert_eq!(
        std::fs::read(&target).unwrap(),
        std::fs::read(&same_f).unwrap()
    );

    // 2. `/dev/stdout` on pipe (`dsp`) interleaves each batch's request before its tabular rows.
    let dsp = clean_command("blastp")
        .env("BATCH_SIZE", "150")
        .args([
            "-query",
            BP_BATCHES_Q,
            "-subject",
            BP_S,
            "-outfmt",
            "6",
            "-export_search_strategy",
            "/dev/stdout",
            "-num_threads",
            "1",
        ])
        .output()
        .expect("run losat dsp");
    assert_eq!(dsp.status.code(), Some(0));
    assert!(dsp.stderr.is_empty());
    let dsp_text = String::from_utf8_lossy(&dsp.stdout);
    assert!(
        dsp_text.starts_with("Blast4-request ::= {"),
        "first request should arrive before report rows"
    );
    assert!(
        dsp_text.matches("Blast4-request ::= {").count() >= 2,
        "expected multi-batch requests in dsp"
    );

    // 3. `>>` append (`dsoapp`): stdout opened with O_APPEND over "PRE\n" keeps "PRE\n" prefix
    // for report writes while `/dev/stdout` opens its own descriptor cleanly without error.
    let app_f = dir.join("app.out");
    std::fs::write(&app_f, b"PRE\n").unwrap();
    let app_file = std::fs::OpenOptions::new()
        .append(true)
        .open(&app_f)
        .unwrap();
    let app_out = clean_command("blastp")
        .env("BATCH_SIZE", "150")
        .args([
            "-query",
            BP_BATCHES_Q,
            "-subject",
            BP_S,
            "-outfmt",
            "6",
            "-export_search_strategy",
            "/dev/stdout",
            "-num_threads",
            "1",
        ])
        .stdout(Stdio::from(app_file))
        .output()
        .expect("run losat dsoapp");
    assert_eq!(app_out.status.code(), Some(0));
    assert!(app_out.stderr.is_empty());
    assert!(!std::fs::read(&app_f).unwrap().is_empty());

    let _ = std::fs::remove_dir_all(&dir);
}

/// `/dev/stdout` with stdout redirected to `/dev/full` (`dsofull`):
/// - outfmt 0 fails in prolog before batch 0's strategy write (`BLAST failed to write output`, rc 6).
/// - outfmt 6 posts the 4-line `Close` chain on batch 0, then aborts when the tabular report fails
///   on `/dev/full` (SIGABRT, `std::__ios_failure`).
#[cfg(unix)]
#[test]
fn dev_stdout_to_dev_full_matches_ncbi_termination() {
    use std::os::unix::process::ExitStatusExt;

    let run_dsofull = |outfmt: &str| -> Output {
        let full = std::fs::OpenOptions::new()
            .write(true)
            .open("/dev/full")
            .unwrap();
        clean_command("blastp")
            .args([
                "-query",
                "tests/fasta/outfmt11/arch_bp_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                outfmt,
                "-export_search_strategy",
                "/dev/stdout",
                "-num_threads",
                "1",
            ])
            .stdout(Stdio::from(full))
            .output()
            .expect("run losat dsofull")
    };

    let f0 = run_dsofull("0");
    assert_eq!(f0.status.code(), Some(6));
    assert_eq!(
        String::from_utf8_lossy(&f0.stderr),
        "BLAST failed to write output\n"
    );

    let f6 = run_dsofull("6");
    assert_eq!(
        f6.status.signal(),
        Some(6),
        "expected SIGABRT on tabular dsofull"
    );
    let err6 = String::from_utf8_lossy(&f6.stderr);
    assert!(err6.contains("cannot close output stream"), "{err6}");
    assert!(err6.contains("std::__ios_failure"), "{err6}");
}

/// `rpstblastn` with `-outfmt 5`, `15`, and `16` and multi-batch `BATCH_SIZE=1` interleaves
/// batch 0's strategy request, batch 0's XML/JSON report, and batch 1's strategy request in
/// NCBI's order when `HeldOutput` is active (regression check for `HeldOutput::push_posts`).
#[cfg(unix)]
#[test]
fn rpstblastn_held_output_interleaves_strategy_and_xml_batches() {
    use std::io::Read;

    let dir = temp_run_dir("rpstblastn_held");
    let rps_dir = dir.join("rpsdb");
    std::fs::create_dir_all(&rps_dir).unwrap();
    let src = engine().join("tests/db_search/rps");
    for entry in std::fs::read_dir(&src).unwrap() {
        let name = entry.unwrap().file_name().to_string_lossy().into_owned();
        if let Some(stem) = name.strip_suffix(".gz") {
            if stem.starts_with("rps_mini.") {
                let mut bytes = Vec::new();
                flate2::read::GzDecoder::new(std::fs::File::open(src.join(&name)).unwrap())
                    .read_to_end(&mut bytes)
                    .unwrap();
                std::fs::write(rps_dir.join(stem), bytes).unwrap();
            }
        }
    }
    let db = rps_dir.join("rps_mini");
    for (outfmt, report_marker) in [
        ("5", "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>"),
        ("15", "\"BlastOutput2\": ["),
        ("16", "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>"),
    ] {
        let out = clean_command("rpstblastn")
            .env("BATCH_SIZE", "1")
            .args([
                "-query",
                "tests/fasta/rps/nq4.fna",
                "-db",
                db.to_str().unwrap(),
                "-comp_based_stats",
                "0",
                "-outfmt",
                outfmt,
                "-export_search_strategy",
                "/dev/stdout",
                "-num_threads",
                "1",
            ])
            .output()
            .expect("run rpstblastn");
        assert_eq!(
            out.status.code(),
            Some(0),
            "outfmt {outfmt} stderr={}",
            String::from_utf8_lossy(&out.stderr)
        );
        let text = String::from_utf8_lossy(&out.stdout);
        let req0 = text.find("Blast4-request ::= {").expect("req0");
        let rep = text.find(report_marker).expect("report start");
        let req1 = text[req0 + 1..]
            .find("Blast4-request ::= {")
            .map(|p| req0 + 1 + p)
            .expect("req1");
        assert!(
            req0 < rep && rep < req1,
            "outfmt {outfmt}: expected req0 ({req0}) < report ({rep}) < req1 ({req1})"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}
