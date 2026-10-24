//! tblastn's units of a long packed subject translate the residues they read
//! (`LOSAT_X_TBNUNITXLATE`, in the strict set; `src/algorithm/tblastn/x_subject_translate.rs`),
//! on whole searches.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
//! ```c
//!             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//!                                      backup.full_range.right,
//!                                      subject->gen_code_string, &translation_buffer,
//!                                      &frame_offsets, NULL);
//! ...
//!             subject->sequence = translation_buffer + frame_offsets[context] + 1;
//!             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
//! ```
//! NCBI translates the six frames of a subject and each unit reads its chunk of a frame. With the
//! switch, the plan of a packed subject whose first frame has 1 MiB or more (3,145,722 bases)
//! made for its search alone (the search keeps no subject plans for several query batches) has
//! no residues, and each unit translates those of its chunk (and
//! a margin) into a buffer of its thread, so the report does not depend on the switch. It is on
//! by default and off with `LOSAT_X_TBNUNITXLATE=0`. This test compares the reports and reads
//! the `[X_STATS] TBNUNITXLATE plans= units= residues=` line.
//!
//! The query is the first protein of `tests/fasta/PemoMJNVA.faa`; the subjects are eleven copies
//! of `tests/fasta/MelaMJNV.fasta`, forward and reverse complement in turn (one long subject),
//! and the genome once more (a short one, whose plan keeps its frames).

use std::path::{Path, PathBuf};
use std::process::Command;

/// The query file and the subject file (written once per test process), and the length of the
/// long subject.
fn inputs() -> &'static (PathBuf, PathBuf, usize) {
    static FILES: std::sync::OnceLock<(PathBuf, PathBuf, usize)> = std::sync::OnceLock::new();
    FILES.get_or_init(|| {
        let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_unit_xlate");
        std::fs::create_dir_all(&dir).expect("work directory");
        let pid = std::process::id();
        let text = std::fs::read_to_string(fasta.join("PemoMJNVA.faa")).expect("query file");
        let second = text[1..].find('>').map_or(text.len(), |at| at + 1);
        let query = dir.join(format!("one.{pid}.faa"));
        std::fs::write(&query, &text[..second]).expect("query file");
        let genome: Vec<u8> = std::fs::read_to_string(fasta.join("MelaMJNV.fasta"))
            .expect("subject file")
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.trim().bytes())
            .collect();
        let reverse: Vec<u8> = genome
            .iter()
            .rev()
            .map(|&base| match base.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect();
        let long: Vec<u8> = (0..11)
            .flat_map(|copy| {
                if copy % 2 == 0 { &genome } else { &reverse }
                    .iter()
                    .copied()
            })
            .collect();
        assert!(long.len() / 3 + 2 >= 1 << 20, "{} bases", long.len());
        let mut file = Vec::new();
        for (name, bases) in [("long", &long), ("genome", &genome)] {
            file.extend_from_slice(format!(">{name}\n").as_bytes());
            for line in bases.chunks(80) {
                file.extend_from_slice(line);
                file.push(b'\n');
            }
        }
        let subject = dir.join(format!("subjects.{pid}.fna"));
        std::fs::write(&subject, file).expect("subject file");
        (query, subject, long.len())
    })
}

/// One TBLASTN run (`-outfmt 6` on standard output) with `variables` and `LOSAT_X_STATS` as its
/// only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build with the worker probe
/// (`--all-features`) writes to that file as soon as a pool starts, so the caller's value is
/// kept, or a file of the tests is named. Returns the report and the standard error.
fn tblastn(threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let (query, subject, _) = inputs();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .env("LOSAT_X_STATS", "1")
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// The fields of the one line that starts with `prefix`; `None` without the line.
fn fields<const N: usize>(stderr: &str, prefix: &str, names: [&str; N]) -> Option<[u64; N]> {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| line.starts_with(prefix))
        .collect();
    assert!(lines.len() <= 1, "one line per process:\n{stderr}");
    let line = lines.first()?;
    Some(names.map(|field| {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(field)?.strip_prefix('='))
            .unwrap_or_else(|| panic!("no {field} in {line}"))
            .parse()
            .expect("a number")
    }))
}

fn stats(stderr: &str) -> Option<[u64; 3]> {
    fields(
        stderr,
        "[X_STATS] TBNUNITXLATE ",
        ["plans", "units", "residues"],
    )
}

/// The residues the six units of a subject of `length` bases with one chunk per frame read: each
/// frame's residues and its end sentinel.
fn unit_residues(length: usize) -> u64 {
    (0..3).map(|i| 2 * ((length - i) / 3 + 1) as u64).sum()
}

const OFF: (&str, &str) = ("LOSAT_X_TBNUNITXLATE", "0");

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812 (above). The
// default and `=0` give the same report at one and four threads. By default the long subject's
// plan has no residues and its six units (one chunk per frame) translate their frames, and the
// short subject's plan keeps its frames; with `=0` no plan is made so.
#[test]
fn the_report_does_not_depend_on_where_the_units_read_their_residues() {
    let long = inputs().2;
    let (expected, stderr) = tblastn(1, &[OFF]);
    assert!(
        String::from_utf8_lossy(&expected)
            .lines()
            .any(|line| line.split('\t').nth(1) == Some("long")),
        "the search must report hits in the long subject"
    );
    assert_eq!(stats(&stderr), None);
    for threads in [1, 4] {
        let (report, stderr) = tblastn(threads, &[]);
        assert!(
            report == expected,
            "{threads} threads: the default changes the report"
        );
        assert_eq!(
            stats(&stderr),
            Some([1, 6, unit_residues(long)]),
            "{threads} threads"
        );
    }
    let (report, _) = tblastn(4, &[OFF]);
    assert!(report == expected, "four threads with =0");
}
