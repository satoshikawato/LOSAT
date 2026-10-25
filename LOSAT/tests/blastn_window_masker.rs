//! `blastn -window_masker_db` with counts files that NCBI cannot use, against the bytes of NCBI BLAST+
//! 2.17.0 (the official Linux x64 release, whose source tree the exception reports name; line NP
//! M/DECISIONS.md 14): a file that does not exist or cannot be read (exit 255, before `-query_loc` is
//! read and before the database is opened; an error archive with `-outfmt 11`), a format other than
//! the two binary ones (exit 1, after the `-dust` check), and damaged counts files, which stop the
//! batch that cannot be masked with one message per query of the batch (split or not: NCBI sets up the
//! whole batch first), after the reports of the batches before and the
//! outfmt 0 prolog (exit 3).
//!
//! The cases are `tests/fixtures/winmask_blastn/cases.tsv` with `<case>.out`, `.err` and `.rc`,
//! written by `tests/winmask/make_blastn_error_fixtures.py freeze`; the counts files under
//! `target/winmask_blastn/` are made here with the bytes that script writes. The searches that
//! succeed are the `wm.*` rows of `tests/blastn_regression_fixtures.py`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-458 (the check of the
//! counts file), c++/src/algo/blast/api/winmask_filter.cpp:146-151 (`s_BuildSeqMasker`),
//! c++/src/algo/blast/api/blast_setup_cxx.cpp:519-652 (the set-up of each query).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Once;

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixtures() -> PathBuf {
    engine().join("tests/fixtures/winmask_blastn")
}

fn words(values: &[u32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}

fn word_at(data: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(data[at..at + 4].try_into().unwrap())
}

/// Writes `content` to `path` unless the file already has it: through a file of its own and a
/// rename, so that another test process (nextest runs each test in one) never reads a file being
/// written.
fn put(path: &Path, content: &[u8]) {
    if std::fs::read(path).is_ok_and(|old| old == content) {
        return;
    }
    let part = path.with_extension(format!("part{}", std::process::id()));
    std::fs::write(&part, content).unwrap();
    std::fs::rename(&part, path).unwrap();
}

/// The counts files of the cases (`make_counts` of make_blastn_error_fixtures.py).
fn make_counts() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let work = engine().join("target/winmask_blastn");
        std::fs::create_dir_all(&work).unwrap();
        let data = std::fs::read(engine().join("tests/winmask/counts/mg1655_u9.obinary")).unwrap();
        let md = word_at(&data, 4) as usize;
        let (head, body) = data.split_at(8 + md);
        let (unit, m, k, bc) = (
            word_at(body, 4),
            word_at(body, 8),
            word_at(body, 12),
            word_at(body, 20),
        );
        let ht_at = 40 + 4 + 4 * ((1usize << (2 * unit)) / 32);
        let mut bad = body.to_vec();
        for index in 0..(1usize << k) {
            let at = ht_at + 4 * index;
            bad[at..at + 4].copy_from_slice(&((m << bc) | 3).to_le_bytes());
        }
        // The same with the entry of key 0 kept: a query without A, C, G or T only reads that
        // one, so a batch of such queries is masked and the next batch fails.
        let mut bad_after0 = bad.clone();
        bad_after0[ht_at..ht_at + 4].copy_from_slice(&body[ht_at..ht_at + 4]);
        let pct_md = b"##windowmasker-statistics-format-version:obinary 2.0.0\0##pct:\0";
        let unit17 = [head, &body[..4], &words(&[17]), &body[8..]].concat();
        let pct6 = [&words(&[3, pct_md.len() as u32])[..], pct_md, body].concat();
        let files: [(&str, Vec<u8>); 8] = [
            ("unit17.obinary", unit17),
            ("short_cba.obinary", [head, &body[..48]].concat()),
            ("bad_index.obinary", [head, &bad].concat()),
            ("bad_index_after0.obinary", [head, &bad_after0].concat()),
            ("pct6.obinary", pct6),
            ("binary_unit0.binary", words(&[0, 0, 5, 1, 1, 2, 3, 4])),
            (
                "ascii.counts",
                b"##windowmasker-statistics-format-version:ascii 1.1.0\n9\nAAAAAAAAA 5\n".to_vec(),
            ),
            ("empty.obinary", Vec::new()),
        ];
        for (name, content) in files {
            put(&work.join(name), &content);
        }
        // A file without read permission, made once (its bytes do not matter).
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let unreadable = work.join("unreadable.obinary");
            if !unreadable.exists() {
                let part = work.join(format!("unreadable.part{}", std::process::id()));
                std::fs::write(&part, &data).unwrap();
                std::fs::set_permissions(&part, std::fs::Permissions::from_mode(0o000)).unwrap();
                std::fs::rename(&part, &unreadable).unwrap();
            }
        }
    });
}

struct Case {
    name: String,
    env: Vec<(String, String)>,
    args: Vec<String>,
}

fn cases() -> Vec<Case> {
    std::fs::read_to_string(fixtures().join("cases.tsv"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            Case {
                name: cols[0].to_string(),
                env: cols[1]
                    .split_whitespace()
                    .map(|item| {
                        let (key, value) = item.split_once('=').unwrap();
                        (key.to_string(), value.to_string())
                    })
                    .collect(),
                args: cols[2]
                    .split(' ')
                    .map(|token| {
                        if token == "''" {
                            String::new()
                        } else {
                            token.replace("%20", " ")
                        }
                    })
                    .collect(),
            }
        })
        .collect()
}

/// Whether the case can run here: a file without read permission is readable to root, and on
/// Windows a directory does not open as a file (NCBI's Linux build reads it as an empty stream).
fn runs_here(case: &Case) -> bool {
    match case.name.as_str() {
        "unreadable" => {
            cfg!(unix)
                && std::fs::File::open(engine().join("target/winmask_blastn/unreadable.obinary"))
                    .is_err()
        }
        "directory" => cfg!(unix),
        _ => true,
    }
}

fn check(threads: &str) {
    make_counts();
    let mut failures = Vec::new();
    for case in cases().iter().filter(|case| runs_here(case)) {
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        command
            .current_dir(engine())
            .arg("blastn")
            .args(&case.args)
            .args(["-num_threads", threads]);
        for (name, _) in std::env::vars_os() {
            let text = name.to_string_lossy();
            if text.starts_with("LOSAT_")
                || text.starts_with("NCBI")
                || matches!(
                    text.as_ref(),
                    "BATCH_SIZE"
                        | "CHUNK_SIZE"
                        | "OVERLAP_CHUNK_SIZE"
                        | "BLASTDB"
                        | "ARCHIVE_FORMAT"
                )
            {
                command.env_remove(name);
            }
        }
        // The --all-features build's diagnostic probe needs its log at 2 threads and more.
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat_winmask_blastn_worker.log"),
        );
        command.envs(case.env.iter().map(|(key, value)| (key, value)));
        let output = command.output().expect("run losat");
        let read =
            |ext: &str| std::fs::read(fixtures().join(format!("{}.{ext}", case.name))).unwrap();
        let rc: i32 = String::from_utf8(read("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        if output.status.code() != Some(rc) {
            failures.push(format!(
                "{}: exit {:?}, NCBI {rc}",
                case.name,
                output.status.code()
            ));
        }
        if output.stdout != read("out") {
            failures.push(format!("{}: stdout differs", case.name));
        }
        if output.stderr != read("err") {
            failures.push(format!(
                "{}: stderr {:?}",
                case.name,
                String::from_utf8_lossy(&output.stderr)
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

// One test for both thread counts: the counts files are made by the test process.
#[test]
fn unusable_counts_files_give_ncbis_errors_at_one_and_four_threads() {
    check("1");
    check("4");
}

#[test]
fn every_case_has_its_frozen_files() {
    for case in cases() {
        for ext in ["out", "err", "rc"] {
            assert!(
                Path::new(&fixtures().join(format!("{}.{ext}", case.name))).exists(),
                "{} .{ext}",
                case.name
            );
        }
    }
}
