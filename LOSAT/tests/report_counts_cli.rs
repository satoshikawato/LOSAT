//! `-num_descriptions` and `-num_alignments` on the command line of every program: the hit
//! list size and the warnings and the error that NCBI's `ExtractAlgorithmOptions` and
//! `BlastHitSavingOptionsValidate` raise (blast_args.cpp:2894-2977, blast_options.c:1511-1517).
//! The output bytes against NCBI are the fixtures of `outfmt0_manifest.tsv` and
//! `db_search_manifest.tsv`.

use std::path::PathBuf;
use std::process::{Command, Output};

const F: &str = "tests/fasta/outfmt0";
const D: &str = "tests/db_search";

fn losat_binary() -> PathBuf {
    if let Some(path) = option_env!("CARGO_BIN_EXE_LOSAT") {
        return PathBuf::from(path);
    }
    let mut path = std::env::current_exe().expect("current test executable");
    path.pop();
    if path.file_name().is_some_and(|name| name == "deps") {
        path.pop();
    }
    path.push(format!("LOSAT{}", std::env::consts::EXE_SUFFIX));
    path
}

/// The inputs of each program (queries with hits in the subjects).
fn inputs(program: &str) -> Vec<String> {
    let (query, subject) = match program {
        "blastn" => (
            format!("{F}/multi_query.fasta"),
            format!("{F}/multi_subject.fasta"),
        ),
        "blastp" => (
            format!("{F}/e2e_many_query.faa"),
            format!("{F}/e2e_many_subject.faa"),
        ),
        "tblastn" => (
            format!("{F}/e2e_many_query.faa"),
            format!("{F}/e2e_many_subject.fna"),
        ),
        "blastx" => (format!("{D}/bslice_q.fna"), format!("{D}/bslice.faa")),
        "tblastx" => (
            format!("{F}/tblastx_multi_query.fasta"),
            format!("{F}/tblastx_multi_subject.fasta"),
        ),
        _ => unreachable!("{program}"),
    };
    ["-query".to_string(), query, "-subject".to_string(), subject].to_vec()
}

/// A temporary directory of one test (removed at its end), so that tests running in
/// parallel share no file.
struct WorkDir(PathBuf);

impl WorkDir {
    fn new(test: &str) -> Self {
        let path =
            std::env::temp_dir().join(format!("report_counts_cli_{}_{test}", std::process::id()));
        std::fs::create_dir_all(&path).expect("create the test directory");
        Self(path)
    }
}

impl Drop for WorkDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn run(work: &WorkDir, program: &str, extra: &[&str]) -> Output {
    let mut command = Command::new(losat_binary());
    // Like the sibling tests (`blastn_qgroup_par.rs`): no `LOSAT_*` switch and none of the
    // batch variables reaches the child. A build with the worker probe (`--all-features`)
    // panics when a pool starts without `LOSAT_BLASTX_WORKER_LOG`, so every spawn gets a
    // file of the test's own directory.
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_X_DEFAULT")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", work.0.join("worker.log"))
        .arg(program)
        .args(inputs(program))
        .args(extra)
        .output()
        .expect("run LOSAT")
}

const PROGRAMS: [&str; 5] = ["blastn", "blastp", "tblastn", "blastx", "tblastx"];

// NCBI blast_options.c:1511-1517: a hit list size of 0 is an option error, after the
// warning for a hit list size below 5 (blast_args.cpp:2975-2977); nothing is searched.
#[test]
fn a_hit_list_size_of_zero_is_ncbis_error_after_the_few_matches_warning() {
    let work = WorkDir::new("a_hit_list_size_of_zero_is_ncb");
    for program in PROGRAMS {
        for extra in [
            &[
                "-outfmt",
                "0",
                "-num_descriptions",
                "0",
                "-num_alignments",
                "0",
            ][..],
            &["-outfmt", "6", "-num_alignments", "0"],
            &[
                "-outfmt",
                "7",
                "-num_alignments",
                "0",
                "-num_descriptions",
                "9",
            ],
        ] {
            let output = run(&work, program, extra);
            let ignored = if extra.contains(&"-num_descriptions") && extra[1] != "0" {
                format!("Warning: [{program}] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n")
            } else {
                String::new()
            };
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                format!(
                    "{ignored}Warning: [{program}] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"
                ),
                "{program} {extra:?}"
            );
            assert_eq!(output.status.code(), Some(1), "{program} {extra:?}");
            assert!(output.stdout.is_empty(), "{program} {extra:?}");
        }
    }
}

// NCBI blast_args.cpp:2949-2955,2975-2977: -num_descriptions is ignored (with a warning)
// for the formats above 4, and the hit list size below 5 is warned about; the hit list size
// of outfmt 0 is the larger of the two counts, so -num_descriptions 3 does not warn.
#[test]
fn the_warnings_of_the_formatting_handler() {
    let work = WorkDir::new("the_warnings_of_the_formatting");
    for program in PROGRAMS {
        let output = run(
            &work,
            program,
            &[
                "-outfmt",
                "6",
                "-num_descriptions",
                "3",
                "-num_alignments",
                "3",
            ],
        );
        assert_eq!(output.status.code(), Some(0), "{program}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.starts_with(&format!(
                "Warning: [{program}] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\nWarning: [{program}] Examining 5 or more matches is recommended\n"
            )),
            "{program}: {stderr}"
        );
        let output = run(&work, program, &["-outfmt", "0", "-num_descriptions", "3"]);
        assert_eq!(output.status.code(), Some(0), "{program}");
        assert!(output.stderr.is_empty(), "{program}");
    }
}

// NCBI blast_args.cpp:2960-2973: the hit list size of a tabular format is -num_alignments,
// so a subject limit of N shows at most N subjects for each query.
#[test]
fn a_tabular_format_shows_num_alignments_subjects() {
    let work = WorkDir::new("a_tabular_format_shows_num_ali");
    for program in PROGRAMS {
        let count = |extra: &[&str]| {
            let output = run(&work, program, extra);
            assert_eq!(output.status.code(), Some(0), "{program} {extra:?}");
            let mut subjects: Vec<(String, String)> = Vec::new();
            for line in String::from_utf8(output.stdout).unwrap().lines() {
                let mut fields = line.split('\t');
                let pair = (
                    fields.next().unwrap().to_string(),
                    fields.next().unwrap().to_string(),
                );
                if !subjects.contains(&pair) {
                    subjects.push(pair);
                }
            }
            let mut per_query = std::collections::BTreeMap::<String, usize>::new();
            for (query, _) in &subjects {
                *per_query.entry(query.clone()).or_default() += 1;
            }
            per_query.values().copied().max().unwrap_or(0)
        };
        let all = count(&["-outfmt", "6"]);
        assert!(
            all > 1,
            "{program}: the input has hits in two subjects of a query"
        );
        assert_eq!(
            count(&["-outfmt", "6", "-num_alignments", "1"]),
            1,
            "{program}"
        );
        assert_eq!(
            count(&["-outfmt", "6", "-max_target_seqs", "1"]),
            1,
            "{program}"
        );
    }
}
