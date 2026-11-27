//! Integration tests for `-export_search_strategy` stream modeling, write-failure handling,
//! shared/aliased non-`-` file sinks, and `2>&1` merged diagnostics (`ostream::StrategyStream`).
//!
//! - Cross-platform: shared `-out`/strategy files (`same` and `./F` aliases) and
//!   `-export_search_strategy -` (NCBI's `cout`, line RG item F1c-2) across 1, 4, and 8 threads.
//! - Linux, item F1c-2 (`tests/fixtures/strategy/f1c2_cases.tsv` and `tests/fixtures/strategy/f1c2/`,
//!   frozen with the official NCBI BLAST+ 2.17.0 at one thread): `-export_search_strategy -` with the
//!   report on the same `cout` (`dash`, every program and format), in `-out F` (`dasho`), in
//!   `-out /dev/stdout` (`revp` on a pipe, `rev` on a file), with `2>&1` and the `803.7` posts before
//!   and after the first piece (4,016 B: the 4,096 B `COStreamBuffer` less the 80 B it keeps back;
//!   `dashm`, `dashom`), standard output on `/dev/full` (`dashfull`:
//!   exit 6 for outfmt 0, SIGABRT for the formats whose request write throws), and `>>` (`dashapp`,
//!   `PRE\n` kept); psiblast's `-in_pssm` and `-in_msa` with rounds (the prolog before the request,
//!   also with `-html`, `-out_pssm`, `2>&1` and `/dev/stdout`); a report that cannot be written
//!   (`outfull`: `-out /dev/full`; `outfsz`: `-out F` at its size limit), whose run ends with the
//!   requests up to the batch whose report failed (`-`, `/dev/stdout` and a file); SAM's `@PG CL:`
//!   tag (the program's path) is left out on both sides.
//! - Linux, item F1c-1 (`tests/fixtures/strategy/f1c1_cases.tsv` and `tests/fixtures/strategy/f1c1/`):
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
//!   - `803.7` replaced-byte title interleaving on `2>&1` (both before and after the end of the
//!     first `COStreamBuffer` piece, its 4,096 B buffer less the 80 B it keeps back);
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

/// `-export_search_strategy -` writes the requests into the report's own standard output (NCBI's
/// `cout`), with the same bytes at 1, 4, and 8 threads across platforms.
#[test]
fn dash_standard_output_matches_across_threads() {
    for (prog, q, s, bsz) in BATCHES {
        let mut baseline: Option<Vec<u8>> = None;
        for threads in ["1", "4", "8"] {
            let out = clean_command(prog)
                .env("BATCH_SIZE", bsz)
                .args(["-query", q, "-subject", s, "-outfmt", "6"])
                .args(["-num_threads", threads])
                .args(["-export_search_strategy", "-"])
                .output()
                .expect("run losat");
            assert_eq!(
                out.status.code(),
                Some(0),
                "{prog} t={threads} stderr={}",
                String::from_utf8_lossy(&out.stderr)
            );
            assert!(out.stderr.is_empty());
            assert!(
                out.stdout.starts_with(b"Blast4-request ::= {"),
                "{prog} t={threads}: the first request comes first"
            );
            match &baseline {
                None => baseline = Some(out.stdout),
                Some(want) => assert_eq!(&out.stdout, want, "{prog} t={threads} differed from t=1"),
            }
        }
    }
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

    /// A table of frozen cases: `<table>_cases.tsv` and the gzipped bytes in `<table>/`.
    #[derive(Clone, Copy)]
    struct Table(&'static str);

    fn fixtures_dir() -> PathBuf {
        engine().join("tests/fixtures/strategy")
    }

    fn parse_cases(table: Table) -> Vec<Case> {
        let file = format!("{}_cases.tsv", table.0);
        let text = std::fs::read_to_string(fixtures_dir().join(&file)).expect(&file);
        let mut out = Vec::new();
        for line in text.lines() {
            if line.is_empty() || line.starts_with('#') || line.starts_with("name\t") {
                continue;
            }
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 8, "bad {file} row: {line}");
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

    fn frozen(table: Table, name: &str, ext: &str) -> Option<Vec<u8>> {
        let path = fixtures_dir()
            .join(table.0)
            .join(format!("{name}.{ext}.gz"));
        let file = std::fs::File::open(path).ok()?;
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(file)
            .read_to_end(&mut bytes)
            .expect("inflate frozen fixture");
        Some(bytes)
    }

    /// The bytes without SAM's `@PG` CL tag (the program's path and arguments): from the `\tCL:`
    /// after `\n@PG\tID:0\tVN:` up to the line's `\tPN:` (`import_search_strategy.rs`'s rule).
    fn without_sam_cl(document: &[u8]) -> Vec<u8> {
        let find = |needle: &[u8], from: usize| {
            document
                .get(from..)?
                .windows(needle.len())
                .position(|window| window == needle)
                .map(|at| at + from)
        };
        let Some(pg) = find(b"\n@PG\tID:0\tVN:", 0) else {
            return document.to_vec();
        };
        let line_end = find(b"\n", pg + 1).unwrap_or(document.len());
        let Some(cl) = find(b"\tCL:", pg + 1).filter(|&at| at < line_end) else {
            return document.to_vec();
        };
        let Some(pn) = find(b"\tPN:", cl).filter(|&at| at < line_end) else {
            return document.to_vec();
        };
        [&document[..cl], &document[pn..]].concat()
    }

    fn map_ncbi(bytes: Vec<u8>) -> Vec<u8> {
        crate::outfmt11_blastx::expected_losat_bytes(bytes)
    }

    fn run_case(table: Table, case: &Case, threads: &str, root: &Path, rps_shared: &Path) {
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

        let sam = case
            .args
            .windows(2)
            .any(|pair| pair[0] == "-outfmt" && pair[1] == "17");
        let sam_form = |bytes: Vec<u8>| if sam { without_sam_cl(&bytes) } else { bytes };
        let got_out = sam_form(got_out);
        let want_out = sam_form(map_ncbi(
            frozen(table, &case.name, "out").unwrap_or_default(),
        ));
        assert!(
            got_out == want_out,
            "{label}: stdout differs ({} B vs NCBI {} B): {}",
            got_out.len(),
            want_out.len(),
            crate::outfmt11_blastx::first_difference(&got_out, &want_out)
        );

        let want_err = map_ncbi(frozen(table, &case.name, "err").unwrap_or_default());
        assert!(
            got_err == want_err,
            "{label}: stderr differs ({} B vs NCBI {} B): {}",
            got_err.len(),
            want_err.len(),
            crate::outfmt11_blastx::first_difference(&got_err, &want_err)
        );

        if uses_f {
            let got_file = std::fs::read(&f_path).ok().map(sam_form);
            let want_file = frozen(table, &case.name, "file")
                .map(map_ncbi)
                .map(sam_form);
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
        // begins with `Blast4-request ::= {`, whereas `bothapp` also opens `-out /dev/stdout` at 0);
        // `dashapp` (`-`, `cout` alone) keeps `PRE\n` and appends.
        if case.stdout == "append" && case.name.starts_with("dashapp.") {
            assert!(
                got_out.starts_with(b"PRE\nBlast4-request ::= {"),
                "{label}: expected cout to append after PRE\\n"
            );
        } else if case.stdout == "append" {
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

    /// Every case of `table` at 1 and 4 threads (at least `min` cases).
    fn frozen_cases_match_ncbi(table: Table, min: usize) {
        let cases = parse_cases(table);
        assert!(
            cases.len() >= min,
            "expected the full {} frozen matrix",
            table.0
        );
        let root = temp_run_dir(&format!("frozen_ncbi_{}", table.0));
        assert_eq!(
            std::fs::metadata(&root).unwrap().blksize(),
            4096,
            "file/append rows require a temporary directory with 4096-byte blocks"
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
                    run_case(table, case, threads, &root, &rps_shared);
                });
            }
        });
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn f1c1_frozen_cases_match_ncbi_at_1_and_4_threads() {
        frozen_cases_match_ncbi(Table("f1c1"), 250);
    }

    #[test]
    fn f1c2_frozen_cases_match_ncbi_at_1_and_4_threads() {
        frozen_cases_match_ncbi(Table("f1c2"), 250);
    }
}
