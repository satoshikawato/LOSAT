//! BLAST XML (outfmt 5) of BLASTP and BLASTN: the one title that NCBI BLAST+ 2.17.0 cannot make (a
//! protein or nucleotide title of punctuation alone, such as `, ,`, on which `x_CleanAndCompress`
//! reads past the end of the string: NCBI crashes) is refused (PD-LOSAT-NCBI-DEFECTS) and the refusal
//! comes before the first byte (for BLASTN, of the query batch that has the title: after the
//! output of the batches before, where a later batch has it).
//!
//! The titles NCBI does write (a byte >= 0x80, a control byte, a title whose encoding NCBI cannot
//! guess; the `-db` name) are not refused: their documents, cut where NCBI's XML stream fails, and
//! their posts are the frozen rows `db.xb.*` of `tests/db_search_manifest/` and the cases of
//! `tests/xml_dbname.rs`.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:219-314
//! (`CDeflineGenerator::x_CleanAndCompress`, `left`, a `size_t`, wraps where a separator reaches
//! the end).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

mod temp_support;
use temp_support::{Finished, TempFile};

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt5_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// The binary with the extra environment variables `env` (NCBI's own `BATCH_SIZE`, for example).
fn run_with(args: &[&std::ffi::OsStr], env: &[(&str, &str)]) -> Finished<Output, TempFile> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
            command.env_remove(name);
        }
    }
    // An --all-features binary at 2+ threads keeps this variable (see the other tests).
    let log = TempFile::keep_on_failure("losat_outfmt5_worker", "worker.log");
    command.env("LOSAT_BLASTX_WORKER_LOG", &log);
    command.envs(env.iter().copied());
    let output = command.args(args).output().expect("run LOSAT");
    Finished::new(output, log)
}

fn fasta(path: &Path, records: &[(&[u8], &[u8])]) {
    let mut text = Vec::new();
    for (title, sequence) in records {
        text.extend_from_slice(b">");
        text.extend_from_slice(title);
        text.push(b'\n');
        text.extend_from_slice(sequence);
        text.push(b'\n');
    }
    fs::write(path, text).expect("write FASTA");
}

fn refused(output: &Output, program: &str, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with punctuation that NCBI BLAST+ 2.17.0 crashes on when it makes the definition line is not supported by LOSAT's {program}"
        )),
        "{stderr}"
    );
}

/// The title of punctuation that NCBI's `x_CleanAndCompress` reads past the end of.
const PUNCTUATION: &[u8] = b", ,";

/// 300 BLASTN queries of 400 letters from a fixed pseudo-random sequence generator; the titles
/// are `q<n> plain` except for the one at `bad_title`. Only the last query is `hit`.
fn many_queries(path: &Path, bad_title: Option<(usize, &[u8])>, hit: &[u8]) {
    let mut state = 12345u32;
    let mut text = Vec::new();
    for index in 0..300 {
        text.push(b'>');
        match bad_title {
            Some((at, title)) if at == index => text.extend_from_slice(title),
            _ => text.extend_from_slice(format!("q{index} plain").as_bytes()),
        }
        text.push(b'\n');
        if index == 299 {
            text.extend_from_slice(hit);
        } else {
            for _ in 0..400 {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                text.push(b"ACGT"[(state >> 24) as usize & 3]);
            }
        }
        text.push(b'\n');
    }
    fs::write(path, text).expect("write queries");
}

/// With more than one query batch (NCBI's `BATCH_SIZE` of 2,000 letters here, the default batch
/// sizes of BLASTN's application for 300 queries of 400 letters in the other) a refusal that only a
/// later batch finds comes after the output of the batches before, as NCBI's formatter writes each
/// batch's reports before it reads the next (the Owner's decision of 2026-11-16 on line RG's Q-X2c,
/// "the results of the query batches before are written"; remaining-gaps DECISIONS D-SO-2):
/// standard output (or the `-out` file) is the start of the document that the same search without
/// the refused title writes, to the end of an `Iteration` before the refused query's batch, without
/// the document's end. With the reports written after the loop (`LOSAT_X_RPTSTREAM=0`) the refusal
/// comes before the first byte. (A title NCBI writes does not make a later batch refuse: its
/// document is the one of `db.xb.*`.)
#[test]
fn a_refusal_in_a_later_batch_comes_after_the_batches_before() {
    let dir = temp_dir("batches");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let hit: Vec<u8> = {
        let mut state = 99u32;
        (0..400)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[(state >> 24) as usize & 3]
            })
            .collect()
    };
    let out = dir.join("o.xml");
    let args_of = |with_out: bool| {
        let mut args = vec![
            "blastn",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-task",
            "blastn",
            "-outfmt",
            "5",
        ];
        if with_out {
            args.extend(["-out", out.to_str().unwrap()]);
        }
        args.into_iter()
            .map(std::ffi::OsStr::new)
            .collect::<Vec<_>>()
    };
    for batch_size in [Some("2000"), None] {
        for streamed in [true, false] {
            let mut env: Vec<(&str, &str)> = batch_size
                .map(|size| ("BATCH_SIZE", size))
                .into_iter()
                .collect();
            if !streamed {
                env.push(("LOSAT_X_RPTSTREAM", "0"));
            }
            for what in ["a query title", "a subject title"] {
                // The search without the refused title: its whole document.
                many_queries(&query, None, &hit);
                fasta(&subject, &[(b"s1 subject", &hit)]);
                let whole = run_with(&args_of(false), &env);
                assert_eq!(whole.status.code(), Some(0), "{what}");
                // a query title: the 291st query; a subject title: the subject of the last
                // query's hit (the only query with one)
                if what == "a query title" {
                    many_queries(&query, Some((290, PUNCTUATION)), &hit);
                } else {
                    fasta(&subject, &[(PUNCTUATION, &hit)]);
                }
                for with_out in [false, true] {
                    let what = format!(
                        "{what}, BATCH_SIZE {batch_size:?}, streamed {streamed}, -out {with_out}"
                    );
                    fs::remove_file(&out).ok();
                    let output = run_with(&args_of(with_out), &env);
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
                    let program = "BLASTN";
                    let title = what.split(',').next().unwrap();
                    assert!(
                        stderr.contains(&format!(
                            "{title} with punctuation that NCBI BLAST+ 2.17.0 crashes on when it makes the definition line is not supported by LOSAT's {program}"
                        )),
                        "{what}: {stderr}"
                    );
                    let written = if with_out {
                        assert!(output.stdout.is_empty(), "{what}: output beside -out");
                        fs::read(&out).unwrap_or_default()
                    } else {
                        output.stdout.clone()
                    };
                    if !streamed {
                        assert!(written.is_empty(), "{what}: part of a document is written");
                        continue;
                    }
                    let iterations = written
                        .windows(b"</Iteration>".len())
                        .filter(|window| *window == b"</Iteration>")
                        .count();
                    assert!(
                        whole.stdout.starts_with(&written)
                            && written.ends_with(b"</Iteration>")
                            && (1..=290).contains(&iterations),
                        "{what}: {iterations} iterations of\n{}",
                        String::from_utf8_lossy(&written)
                    );
                }
            }
        }
    }
    fs::remove_dir_all(dir).ok();
}
