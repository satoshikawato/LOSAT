//! Integration tests for `-export_search_strategy` stream modeling, write-failure handling,
//! shared/aliased non-`-` file sinks, and `2>&1` merged diagnostics (`ostream::StrategyStream`).
//!
//! - Cross-platform: shared `-out`/strategy files (`same` and `./F` aliases) across 1, 4, and 8
//!   threads, and preservation of the `-export_search_strategy -` refusal.
//! - Linux (`tests/fixtures/strategy/f1c1_cases.tsv` and `tests/fixtures/strategy/f1c1/`):
//!   byte-for-byte comparison against frozen official NCBI BLAST+ 2.17.0 outputs (mapped via
//!   `losat::wording::losat_output_bytes`) at 1 and 4 threads across all 8 programs, covering:
//!   - `2>&1` (`merged`) strategy write failures (`/dev/full` and `RLIMIT_FSIZE=3000`) for
//!     `blastp`, `tblastn`, `psiblast`, `rpsblast` outfmt `0..=4` and all programs across formats;
//!   - `/dev/full` and `RLIMIT_FSIZE` (`100`, `3000`, `5000`, `10000`) strategy failures, including
//!     4-line `Close` failure chains, subsequent-batch `write fault`, `-outfmt 11` error archives,
//!     and large single-batch requests (> 8,191 B) faulting at `byte 12128` (`ncbieaa`) or
//!     `line N: <stack>: write fault` (`ncbi2na` `~ByteBlock()` unwind);
//!   - pipe shared sinks (`/dev/stdout`, `-out /dev/stdout -export_search_strategy /dev/stdout`,
//!     `/dev/stdout -out F`, and `blastx`/`rpstblastn` held `outfmt 5/15/16` interleaving);
//!   - regular-file, `./F`, symlink (`L -> F`), `/dev/stdout` file, and `>>` append (`dsoapp` /
//!     `bothapp`, where `/dev/stdout` opens `O_TRUNC` and erases pre-existing `PRE\n`, matching
//!     NCBI's output starting with `Blast4-request ::= {`);
//!   - `803.7` replaced-byte title interleaving on `2>&1` (both before and after the 4,096 B
//!     `COStreamBuffer` piece boundary);
//!   - stdout redirected to `/dev/full` with `-export_search_strategy /dev/stdout` (`dsofull`:
//!     exit 6 for outfmt 0, SIGABRT with Close failure chain + `std::__ios_failure` for tabular /
//!     XML / archive / organism reports).

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

const BATCHES: [(&str, &str, &str, &str); 5] = [
    (
        "blastn",
        "tests/fasta/outfmt11/arch_batches.fna",
        "tests/fasta/outfmt11/arch_s.fna",
        "500",
    ),
    (
        "blastp",
        "tests/fasta/outfmt11/arch_bp_batches.faa",
        "tests/fasta/outfmt11/arch_bp_s.faa",
        "150",
    ),
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

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn clean_env(cmd: &mut Command) {
    cmd.env("BLAST_USAGE_REPORT", "false")
        .env("LOSAT_BLASTX_WORKER_LOG", "/dev/null");
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || name.starts_with("NCBI")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE" | "BLASTDB" | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            cmd.env_remove(name);
        }
    }
}

fn clean_command(program: &str) -> Command {
    let mut cmd = Command::new(BIN);
    cmd.current_dir(engine()).arg(program).stdin(Stdio::null());
    clean_env(&mut cmd);
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
/// 1, 4, and 8 threads across platforms.
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
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("-export_search_strategy - (the search strategy on the standard output"),
        "unexpected stderr: {stderr}"
    );
}

#[cfg(target_os = "linux")]
mod ncbi_frozen {
    use super::*;
    use std::io::Read;
    use std::os::unix::fs::MetadataExt;
    use std::os::unix::process::ExitStatusExt;

    struct Case {
        name: String,
        program: String,
        env: Vec<(String, String)>,
        stdout: String,
        fsize: Option<String>,
        setup: String,
        args: Vec<String>,
        rc: String,
    }

    fn fixtures_dir() -> PathBuf {
        engine().join("tests/fixtures/strategy")
    }

    fn parse_cases() -> Vec<Case> {
        let text =
            std::fs::read_to_string(fixtures_dir().join("f1c1_cases.tsv")).expect("f1c1_cases.tsv");
        let mut out = Vec::new();
        for line in text.lines() {
            if line.is_empty() || line.starts_with('#') || line.starts_with("name\t") {
                continue;
            }
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 8, "bad f1c1_cases.tsv row: {line}");
            let env = if cells[2] == "-" {
                Vec::new()
            } else {
                cells[2]
                    .split(',')
                    .map(|item| {
                        let (k, v) = item.split_once('=').expect("KEY=VAL");
                        (k.to_string(), v.to_string())
                    })
                    .collect()
            };
            out.push(Case {
                name: cells[0].to_string(),
                program: cells[1].to_string(),
                env,
                stdout: cells[3].to_string(),
                fsize: (cells[4] != "-").then(|| cells[4].to_string()),
                setup: cells[5].to_string(),
                args: cells[6].split(' ').map(str::to_string).collect(),
                rc: cells[7].to_string(),
            });
        }
        out
    }

    fn inflate_rps(spec: &str, dir: &Path) -> PathBuf {
        let target = dir.join("rpsdb");
        std::fs::create_dir_all(&target).unwrap();
        let src = engine().join("tests/db_search/rps");
        for entry in std::fs::read_dir(&src).unwrap() {
            let file = entry.unwrap().file_name().to_string_lossy().into_owned();
            let Some(stem) = file.strip_suffix(".gz") else {
                continue;
            };
            if !spec
                .split('+')
                .any(|name| stem.starts_with(&format!("{name}.")))
            {
                continue;
            }
            let mut bytes = Vec::new();
            flate2::read::GzDecoder::new(std::fs::File::open(src.join(&file)).unwrap())
                .read_to_end(&mut bytes)
                .unwrap();
            std::fs::write(target.join(stem), bytes).unwrap();
        }
        target
    }

    fn frozen(name: &str, ext: &str) -> Option<Vec<u8>> {
        let path = fixtures_dir().join("f1c1").join(format!("{name}.{ext}.gz"));
        let file = std::fs::File::open(path).ok()?;
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(file)
            .read_to_end(&mut bytes)
            .expect("inflate f1c1 fixture");
        Some(bytes)
    }

    fn map_ncbi(bytes: Vec<u8>) -> Vec<u8> {
        crate::outfmt11_blastx::expected_losat_bytes(bytes)
    }

    fn run_case(case: &Case, threads: &str, root: &Path, rps_shared: &Path) {
        let dir = root.join(format!("{}_t{threads}", case.name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let f_path = dir.join("F");
        let dot_path = format!("{}/./F", dir.display());
        let lnk_path = dir.join("L");
        if case.setup == "symlink" {
            std::fs::write(&f_path, b"").unwrap();
            std::os::unix::fs::symlink("F", &lnk_path).unwrap();
        }

        let uses_f = case
            .args
            .iter()
            .any(|a| a.contains("{f}") || a.contains("{dot}") || a.contains("{lnk}"));
        let mut argv: Vec<String> = Vec::with_capacity(case.args.len() + 3);
        argv.push(case.program.clone());
        for a in &case.args {
            argv.push(
                a.replace("{f}", f_path.to_str().unwrap())
                    .replace("{dot}", &dot_path)
                    .replace("{lnk}", lnk_path.to_str().unwrap()),
            );
        }
        argv.push("-num_threads".to_string());
        argv.push(threads.to_string());

        let mut cmd = match (case.stdout.as_str(), case.fsize.as_deref()) {
            ("merged", Some(limit)) => {
                let mut c = Command::new("sh");
                c.arg("-c")
                    .arg(r#"trap '' XFSZ; l=$1; shift; exec prlimit --fsize="$l" -- "$@" 2>&1"#)
                    .arg("sh")
                    .arg(limit)
                    .arg(BIN)
                    .args(&argv);
                c
            }
            ("merged", None) => {
                let mut c = Command::new("sh");
                c.arg("-c")
                    .arg(r#"exec "$@" 2>&1"#)
                    .arg("sh")
                    .arg(BIN)
                    .args(&argv);
                c
            }
            ("pipe", Some(limit)) => {
                let mut c = Command::new("sh");
                c.arg("-c")
                    .arg(r#"trap '' XFSZ; l=$1; shift; exec prlimit --fsize="$l" -- "$@""#)
                    .arg("sh")
                    .arg(limit)
                    .arg(BIN)
                    .args(&argv);
                c
            }
            (_, None) => {
                let mut c = Command::new(BIN);
                c.args(&argv);
                c
            }
            (other, Some(limit)) => panic!("unsupported stdout={other} with fsize={limit}"),
        };
        cmd.current_dir(engine()).stdin(Stdio::null());
        clean_env(&mut cmd);
        cmd.env("HOME", &dir);
        for (k, v) in &case.env {
            let val = if let Some(spec) = v.strip_prefix("{rps:").and_then(|s| s.strip_suffix('}'))
            {
                if spec == "rps_mini" {
                    rps_shared.to_str().unwrap().to_string()
                } else {
                    inflate_rps(spec, &dir).to_str().unwrap().to_string()
                }
            } else {
                v.clone()
            };
            cmd.env(k, val);
        }

        let so_path = dir.join("stdout.out");
        let out_file = match case.stdout.as_str() {
            "pipe" | "merged" => None,
            "file" => Some(std::fs::File::create(&so_path).unwrap()),
            "append" => {
                std::fs::write(&so_path, b"PRE\n").unwrap();
                Some(
                    std::fs::OpenOptions::new()
                        .append(true)
                        .open(&so_path)
                        .unwrap(),
                )
            }
            "full" => Some(
                std::fs::OpenOptions::new()
                    .write(true)
                    .open("/dev/full")
                    .unwrap(),
            ),
            other => panic!("unknown stdout mode: {other}"),
        };
        if let Some(f) = out_file {
            cmd.stdout(Stdio::from(f));
        }
        let output = cmd.output().expect("run losat");

        let label = format!("{} (threads={threads})", case.name);
        let got_rc = match output.status.code() {
            Some(code) => code.to_string(),
            None => format!("sig:{}", output.status.signal().unwrap_or(0)),
        };
        let mut got_out = match case.stdout.as_str() {
            "file" | "append" => std::fs::read(&so_path).unwrap(),
            "full" => Vec::new(),
            _ => output.stdout,
        };
        let mut got_err = output.stderr;
        if case.fsize.is_some() {
            if case.stdout == "merged" {
                got_out = crate::llvm_profile::without_profile_errors(&got_out);
            } else {
                got_err = crate::llvm_profile::without_profile_errors(&got_err);
            }
        }

        assert_eq!(
            got_rc,
            case.rc,
            "{label}: exit status; stderr={}",
            String::from_utf8_lossy(&got_err)
        );

        let want_out = map_ncbi(frozen(&case.name, "out").unwrap_or_default());
        assert!(
            got_out == want_out,
            "{label}: stdout differs ({} B vs NCBI {} B): {}",
            got_out.len(),
            want_out.len(),
            crate::outfmt11_blastx::first_difference(&got_out, &want_out)
        );

        let want_err = map_ncbi(frozen(&case.name, "err").unwrap_or_default());
        assert!(
            got_err == want_err,
            "{label}: stderr differs ({} B vs NCBI {} B): {}",
            got_err.len(),
            want_err.len(),
            crate::outfmt11_blastx::first_difference(&got_err, &want_err)
        );

        if uses_f {
            let got_file = std::fs::read(&f_path).ok();
            let want_file = frozen(&case.name, "file").map(map_ncbi);
            assert_eq!(
                got_file.is_some(),
                want_file.is_some(),
                "{label}: file presence mismatch"
            );
            if let (Some(got), Some(want)) = (&got_file, &want_file) {
                assert!(
                    got == want,
                    "{label}: file F differs ({} B vs NCBI {} B): {}",
                    got.len(),
                    want.len(),
                    crate::outfmt11_blastx::first_difference(got, want)
                );
            }
        }

        // Explicitly confirm `dsoapp`/`bothapp` matches NCBI's `/dev/stdout` O_TRUNC behavior
        // (erasing pre-existing `PRE\n`; for `dsoapp`, `cout` appends after the strategy so output
        // begins with `Blast4-request ::= {`, whereas `bothapp` also opens `-out /dev/stdout` at 0).
        if case.stdout == "append" {
            assert!(
                !got_out.starts_with(b"PRE\n"),
                "{label}: expected /dev/stdout O_TRUNC to erase PRE\\n"
            );
            if case.name.starts_with("dsoapp.") {
                assert!(
                    got_out.starts_with(b"Blast4-request ::= {"),
                    "{label}: expected dsoapp output to begin with Blast4-request ::= {{"
                );
            }
        }

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn f1c1_frozen_cases_match_ncbi_at_1_and_4_threads() {
        let cases = parse_cases();
        assert!(cases.len() >= 250, "expected full F1c-1 frozen matrix");
        let root = temp_run_dir("frozen_ncbi");
        assert_eq!(
            std::fs::metadata(&root).unwrap().blksize(),
            4096,
            "F1c-1 file/append rows require a temporary directory with 4096-byte blocks"
        );
        let rps_shared = inflate_rps("rps_mini", &root);
        let work: Vec<(&Case, &str)> = cases
            .iter()
            .flat_map(|case| [(case, "1"), (case, "4")])
            .collect();
        let next = std::sync::atomic::AtomicUsize::new(0);
        let workers = std::thread::available_parallelism().map_or(4, |n| n.get().clamp(1, 4));
        std::thread::scope(|scope| {
            for _ in 0..workers {
                scope.spawn(|| loop {
                    let idx = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(&(case, threads)) = work.get(idx) else {
                        break;
                    };
                    run_case(case, threads, &root, &rps_shared);
                });
            }
        });
        let _ = std::fs::remove_dir_all(&root);
    }
}
