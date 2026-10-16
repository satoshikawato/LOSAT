//! BLASTN hands each subject's HSPs to the hit lists in OID order as it is searched (one
//! thread) or from a window of subjects searched ahead (`LOSAT_X_BNSTREAM`, two or more
//! threads, on by default; `LOSAT_X_BNSTREAM=0` or `LOSAT_X_DEFAULT=0` hands them over after
//! all subjects are searched). The hit lists, and so the report, must not depend on any of
//! these: the searches have more subjects than the window of two or four threads holds (64
//! subjects per thread), and a short hit list that drops subjects.

use std::path::{Path, PathBuf};
use std::process::Command;

/// One FASTA file of 300 subjects cut from the NZ_CP006932 genome (2,000 bases each, the
/// first starting at 0, then every 2,100 bases, wrapping), and a query of three of them.
fn inputs() -> (PathBuf, PathBuf) {
    let genome = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/NZ_CP006932.fasta"),
    )
    .expect("read the genome");
    let letters: Vec<u8> = genome
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect();
    let piece = |index: usize| -> String {
        let start = index * 2_100 % (letters.len() - 2_000);
        String::from_utf8(letters[start..start + 2_000].to_vec()).expect("ASCII")
    };
    let dir = std::env::temp_dir().join(format!("losat-bnstream-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("make the directory");
    let subjects = dir.join("subjects.fna");
    let query = dir.join("query.fna");
    let mut text = String::new();
    for index in 0..300 {
        text.push_str(&format!(">s{index}\n{}\n", piece(index)));
    }
    std::fs::write(&subjects, text).expect("write the subjects");
    let mut text = String::new();
    for index in [7, 150, 299] {
        text.push_str(&format!(">q{index}\n{}\n", piece(index)));
    }
    std::fs::write(&query, text).expect("write the query");
    (query, subjects)
}

fn blastn(
    query: &Path,
    subjects: &Path,
    task: &str,
    threads: u32,
    extra: &[&str],
    switches: &[(&str, &str)],
) -> Vec<u8> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bnstream-worker.log"),
        );
    }
    let output = command
        .arg("blastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects)
        .args(["-task", task, "-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .args(extra)
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{task} {extra:?} {switches:?} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

#[test]
fn the_report_does_not_depend_on_when_subjects_reach_the_hit_lists() {
    let (query, subjects) = inputs();
    // `-task blastn` is left out: on these subjects a debug build stops at the
    // `ungapped.q_start >= q_context_start` assertion of the ungapped stage at one thread on main as
    // well (a release build gives NCBI's bytes); the hit-list path under test is shared by all tasks.
    for task in ["megablast", "dc-megablast"] {
        for extra in [&[][..], &["-max_target_seqs", "5"][..]] {
            let expected = blastn(&query, &subjects, task, 1, extra, &[]);
            assert!(!expected.is_empty(), "{task} {extra:?}: no hits");
            for threads in [2, 4] {
                for switches in [
                    &[][..],
                    &[("LOSAT_X_BNSTREAM", "0")][..],
                    &[("LOSAT_X_DEFAULT", "0")][..],
                    &[("LOSAT_X_DEFAULT", "0"), ("LOSAT_X_BNSTREAM", "1")][..],
                ] {
                    let report = blastn(&query, &subjects, task, threads, extra, switches);
                    assert_eq!(
                        report, expected,
                        "{task} {extra:?} {switches:?} threads {threads}"
                    );
                }
            }
        }
    }
    let _ = std::fs::remove_dir_all(query.parent().expect("a directory"));
}
