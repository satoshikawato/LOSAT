//! The discontiguous megablast tables filled by the search pool (`LOSAT_X_DISCLUTPAR`,
//! `src/algorithm/blastn/x_disc_lookup_par.rs`) on whole searches.
//!
//! The switch is in the strict set: on unless `LOSAT_X_DISCLUTPAR=0` (or `LOSAT_X_DEFAULT=0`). It
//! chooses between two LOSAT code paths with the same result: the tables of `s_FillDiscMBTable`
//! filled by one loop, or by the threads of the search pool, each owning whole blocks of cells of
//! both templates' tables. So the property under test is that the report does not depend on it:
//! for `-task dc-megablast` with one template and with two (`coding_and_optimal`), template
//! lengths 16, 18 and 21 and word sizes 11 and 12, with `-subject` and with `-db`, at 1 and 4
//! threads, the searches with the switch off (`=0`), by default, with it set to 1, and with the
//! shadow (`LOSAT_X_DISCLUTPARSHADOW`, which builds the tables the other way as well and compares
//! every element) give the same bytes. The shadow's counters show that the pool built the tables at 4 threads and not at 1,
//! and how many tables had two templates. The same holds with `LOSAT_X_PAIRPARDISC`, in
//! `LOSAT_X_QBATCHPAR` slots (private pools of two threads), for split query chunks
//! (`CHUNK_SIZE`), and for megablast with a template (`-word_size 11 -template_type ...`);
//! megablast without a template never builds this table. One search fills its table in several
//! chunks of the pool fill (a query of 0.66 Mb, both strands), with one template and with two.
//! That the parts run on several workers of the pool is checked by the unit test
//! `the_fill_runs_on_several_workers_of_the_pool`; that the bytes equal NCBI's is the business of
//! the blastn fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The switches that joined the default set on 2026-10-17 (perf/perfmb-default-bundle), and
/// `LOSAT_X_DCCSR` (2026-10-21; a one-template table above its size rule is built as CSR, not by
/// the fill these tests compare), each `0`, the default before, unless the search sets it: every search here
/// runs as it was written, with the switches off that it does not name (the default runs in
/// `default_bundle.rs`).
const PREVIOUS_DEFAULT: &[(&str, &str)] = &[
    ("LOSAT_X_BPDBCUT", "0"),
    ("LOSAT_X_BPLKPAHEAD", "0"),
    ("LOSAT_X_BPREDOLPT", "0"),
    ("LOSAT_X_BPREDOPAR", "0"),
    ("LOSAT_X_BPSUBJMETA", "0"),
    ("LOSAT_X_DCCSR", "0"),
    ("LOSAT_X_DIAGLV", "0"),
    ("LOSAT_X_PAIRPARDISC", "0"),
    ("LOSAT_X_PAIRPARFOLD", "0"),
    ("LOSAT_X_QBATCHPAR", "0"),
    ("LOSAT_X_QDUSTPAR", "0"),
    ("LOSAT_X_QPREPPAR", "0"),
    ("LOSAT_X_RPTSTREAM", "0"),
];

const SWITCH: (&str, &str) = ("LOSAT_X_DISCLUTPAR", "1");
const OFF: (&str, &str) = ("LOSAT_X_DISCLUTPAR", "0");
const SHADOW: (&str, &str) = ("LOSAT_X_DISCLUTPARSHADOW", "1");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of the first record of a FASTA file.
fn letters(name: &str) -> Vec<u8> {
    std::fs::read(fasta(name))
        .expect("read FASTA")
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A path in the temporary directory that no other test (of this or another process) uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_disclutpar_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A file or directory removed when dropped.
struct Temp(PathBuf);

impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

impl std::ops::Deref for Temp {
    type Target = Path;

    fn deref(&self) -> &Path {
        &self.0
    }
}

/// A query file: 14 pieces of LC738874 of 900 to 4,500 letters, every fourth with a run of `N`s
/// and every third with an IUPAC ambiguity code, then a query shorter than every template
/// (`short`, 12 letters), one of `N`s only (`allN`, for which NCBI warns on standard error) and
/// a poly-A query (`polyA`): 17 queries.
fn queries() -> Temp {
    let genome = letters("LC738874.fasta");
    let mut text = Vec::new();
    let mut put = |title: &str, seq: &[u8]| {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in seq.chunks(61) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    };
    for i in 0..14usize {
        let start = i * 18_000 + 211 * i;
        let len = 900 + (i * 1_531) % 3_600;
        let mut seq = genome[start..start + len].to_vec();
        if i % 4 == 0 {
            let at = len / 2;
            seq.splice(at..at, std::iter::repeat_n(b'N', 25));
        }
        if i % 3 == 0 {
            seq[len / 3] = b"RYKMSW"[i % 6];
        }
        put(&format!("piece{i}"), &seq);
    }
    put("short", b"ACGTACGTTTGA");
    put("allN", &[b'N'; 60]);
    put("polyA", &[b'A'; 700]);
    let path = temp_path("queries.fa");
    std::fs::write(&path, text).expect("write the queries");
    Temp(path)
}

/// The database made by `losat makedb` from LC738868, LC738869 and LC738870 (in a directory of
/// this test, removed when dropped), and its name.
fn database() -> (Temp, String) {
    let dir = Temp(temp_path("db"));
    std::fs::create_dir_all(&dir.0).unwrap();
    let set = dir.0.join("set.fna");
    let mut bytes = Vec::new();
    for name in ["LC738868", "LC738869", "LC738870"] {
        bytes.extend(std::fs::read(fasta(&format!("{name}.fasta"))).unwrap());
    }
    std::fs::write(&set, bytes).unwrap();
    let db = dir.0.join("set").to_str().unwrap().to_string();
    let output = losat()
        .arg("makedb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (dir, db)
}

/// The binary without any `LOSAT_*` variable (and `BATCH_SIZE` / `CHUNK_SIZE` /
/// `OVERLAP_CHUNK_SIZE`) that chooses a code path. `LOSAT_BLASTX_WORKER_LOG` stays (and is set to
/// a file of this process when the environment has none): a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join(format!(
                "losat_disclutpar_{}_worker.log",
                std::process::id()
            )),
        );
    }
    command
}

type Var = (&'static str, &'static str);

/// One search of `queries` against `subject` (`-subject FILE` or `-db NAME`), `-outfmt 6` on
/// standard output, with `vars` as its only variables that choose a code path. Returns the report
/// and the standard error.
fn blastn(
    queries: &Path,
    subject: &[&str],
    extra: &[&str],
    threads: u32,
    vars: &[Var],
) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(queries)
        .args(subject)
        .args(extra)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .envs(PREVIOUS_DEFAULT.iter().copied())
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "{extra:?} {vars:?} threads {threads}: {stderr}"
    );
    (output.stdout, stderr)
}

/// The shadow's line `[X_SHADOW] DISCLUTPAR tables=T pool=P two_templates=W cells=C diffs=0
/// identical` in `stderr`, whose other lines must be `expected` (the standard error of the run
/// without the switch, such as the warning for the query of `N`s only): `(T, P, W)`.
fn shadow_counts(stderr: &str, expected: &str) -> (u64, u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] DISCLUTPAR "))
        .unwrap_or_else(|| panic!("no shadow line:\n{stderr}"));
    assert!(line.ends_with(" diffs=0 identical"), "{line}");
    let rest: Vec<&str> = stderr.lines().filter(|l| *l != line).collect();
    assert_eq!(rest, expected.lines().collect::<Vec<_>>(), "{stderr}");
    let field = |name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
    };
    if field("tables=") > 0 {
        assert!(field("cells=") > 0, "{line}");
    }
    (field("tables="), field("pool="), field("two_templates="))
}

/// For one search (`subject` and `extra` arguments): the report with the switch off at one
/// thread, then off, by default (on), on (`=1`) and the shadow at 1 and 4 threads, all the same
/// bytes; the switch writes nothing, the shadow one line whose tables were built on the pool only
/// at 4 threads, with two templates exactly when `two_templates`.
fn check(subject: &[&str], extra: &[&str], two_templates: bool, env: &[Var]) {
    let queries = queries();
    let with = |more: &[Var]| -> Vec<Var> { env.iter().chain(more).copied().collect() };
    let (expected, expected_stderr) = blastn(&queries, subject, extra, 1, &with(&[OFF]));
    assert!(!expected.is_empty(), "{extra:?}: no hits");
    for threads in [1, 4] {
        let (report, stderr) = blastn(&queries, subject, extra, threads, &with(&[OFF]));
        assert!(
            report == expected,
            "{extra:?} threads {threads}: off differs"
        );
        assert_eq!(stderr, expected_stderr, "{extra:?} threads {threads}");
        for (on, what) in [(&[][..], "default"), (&[SWITCH][..], "on")] {
            let (report, stderr) = blastn(&queries, subject, extra, threads, &with(on));
            assert!(
                report == expected,
                "{extra:?} threads {threads}: {what} differs"
            );
            assert_eq!(
                stderr, expected_stderr,
                "{extra:?} threads {threads} {what}"
            );
        }
        let (report, stderr) = blastn(&queries, subject, extra, threads, &with(&[SHADOW]));
        assert!(
            report == expected,
            "{extra:?} threads {threads}: shadow differs"
        );
        let (tables, pool, two) = shadow_counts(&stderr, &expected_stderr);
        assert!(tables > 0, "{extra:?} threads {threads}: {stderr}");
        if threads == 1 {
            assert_eq!(pool, 0, "{extra:?}: {stderr}");
        } else {
            assert_eq!(pool, tables, "{extra:?}: {stderr}");
        }
        assert_eq!(two, if two_templates { tables } else { 0 }, "{stderr}");
    }
}

fn subject() -> Vec<String> {
    vec![
        "-subject".into(),
        fasta("LC738868.fasta").to_str().unwrap().into(),
    ]
}

fn args(v: &[String]) -> Vec<&str> {
    v.iter().map(String::as_str).collect()
}

#[test]
fn dc_megablast_default_template() {
    let s = subject();
    check(&args(&s), &["-task", "dc-megablast"], false, &[]);
}

#[test]
fn dc_megablast_two_templates_of_21() {
    let s = subject();
    check(
        &args(&s),
        &[
            "-task",
            "dc-megablast",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "21",
        ],
        true,
        &[],
    );
}

#[test]
fn dc_megablast_optimal_16_word_12() {
    let s = subject();
    check(
        &args(&s),
        &[
            "-task",
            "dc-megablast",
            "-word_size",
            "12",
            "-template_type",
            "optimal",
            "-template_length",
            "16",
        ],
        false,
        &[],
    );
}

#[test]
fn dc_megablast_two_templates_of_18_word_12() {
    let s = subject();
    check(
        &args(&s),
        &[
            "-task",
            "dc-megablast",
            "-word_size",
            "12",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "18",
        ],
        true,
        &[],
    );
}

#[test]
fn dc_megablast_database_one_and_two_templates() {
    let (_dir, db) = database();
    let subject = ["-db", db.as_str()];
    check(&subject, &["-task", "dc-megablast"], false, &[]);
    check(
        &subject,
        &[
            "-task",
            "dc-megablast",
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "16",
        ],
        true,
        &[],
    );
}

// With the two-hit seed stage on the pool as well (LOSAT_X_PAIRPARDISC).
#[test]
fn dc_megablast_with_pairpardisc() {
    let s = subject();
    check(
        &args(&s),
        &["-task", "dc-megablast"],
        false,
        &[("LOSAT_X_PAIRPARDISC", "1")],
    );
}

// The batches searched in LOSAT_X_QBATCHPAR slots (2 slots of 2 threads, each with a private
// pool, batches of 10,000 letters): the slots' tables are built on their pools, and the report is
// the same. Then a batch split into query chunks (`CHUNK_SIZE`), each chunk with its own table.
#[test]
fn qbatchpar_slots_and_split_chunks() {
    let queries = queries();
    let s = subject();
    let subject = args(&s);
    let extra = ["-task", "dc-megablast"];
    for (env, slots) in [
        (
            &[("BATCH_SIZE", "10000")][..],
            &[("LOSAT_X_QBATCHPAR", "1"), ("LOSAT_X_QBATCHPARSLOTS", "2")][..],
        ),
        (
            &[("BATCH_SIZE", "1000000"), ("CHUNK_SIZE", "7000")][..],
            &[][..],
        ),
    ] {
        let with =
            |more: &[Var]| -> Vec<Var> { env.iter().chain(slots).chain(more).copied().collect() };
        let (expected, expected_stderr) = blastn(&queries, &subject, &extra, 1, env);
        assert!(!expected.is_empty(), "{env:?}: no hits");
        let (report, stderr) = blastn(&queries, &subject, &extra, 4, &with(&[OFF]));
        assert!(report == expected, "{env:?} {slots:?}: off differs");
        assert_eq!(stderr, expected_stderr, "{env:?}");
        for on in [&[][..], &[SWITCH][..]] {
            let (report, stderr) = blastn(&queries, &subject, &extra, 4, &with(on));
            assert!(report == expected, "{env:?} {slots:?} {on:?}: on differs");
            assert_eq!(stderr, expected_stderr, "{env:?}");
        }
        let (report, stderr) = blastn(&queries, &subject, &extra, 4, &with(&[SHADOW]));
        assert!(report == expected, "{env:?} {slots:?}: shadow differs");
        let (tables, pool, _) = shadow_counts(&stderr, &expected_stderr);
        assert!(tables > 1 && pool == tables, "{env:?}: {stderr}");
    }
}

// Megablast builds this table only with a template (and word size 11 or 12); without one, the
// shadow compares nothing.
#[test]
fn megablast_with_and_without_a_template() {
    let s = subject();
    check(
        &args(&s),
        &[
            "-task",
            "megablast",
            "-word_size",
            "11",
            "-template_type",
            "coding",
            "-template_length",
            "18",
        ],
        false,
        &[],
    );
    let queries = queries();
    let (expected, expected_stderr) = blastn(&queries, &args(&s), &["-task", "megablast"], 1, &[]);
    let (report, stderr) = blastn(&queries, &args(&s), &["-task", "megablast"], 4, &[SHADOW]);
    assert!(report == expected, "megablast: shadow differs");
    assert_eq!(shadow_counts(&stderr, &expected_stderr), (0, 0, 0));
}

// With LOSAT_TIMING the builds print their summed wall time at exit, with or without the switch,
// and how many were built on a pool: all of them with the switch at 4 threads (by default and set
// to 1), none with it off (`=0`, or `LOSAT_X_DEFAULT=0`) or at 1 thread; and the split line: the
// chunks and phases of the pool fills (and its prefault rounds, at most one per fill), or the time
// of the original loop.
#[test]
fn timing_line_is_printed() {
    let queries = queries();
    let s = subject();
    let subject = args(&s);
    let extra = ["-task", "dc-megablast"];
    let (expected, _) = blastn(&queries, &subject, &extra, 1, &[]);
    for (threads, vars, on_pool) in [
        (4, &[("LOSAT_TIMING", "1"), SWITCH][..], true),
        (4, &[("LOSAT_TIMING", "1")][..], true),
        (4, &[("LOSAT_TIMING", "1"), OFF][..], false),
        (
            4,
            &[("LOSAT_TIMING", "1"), ("LOSAT_X_DEFAULT", "0")][..],
            false,
        ),
        (1, &[("LOSAT_TIMING", "1"), SWITCH][..], false),
    ] {
        let (report, stderr) = blastn(&queries, &subject, &extra, threads, vars);
        assert!(report == expected, "{vars:?} {threads}: the report differs");
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] disc_lookup_build: "))
            .unwrap_or_else(|| panic!("{vars:?}: no timing line:\n{stderr}"));
        let field = |name: &str| -> u64 {
            line.split([' ', '(', ')'])
                .find_map(|word| word.strip_prefix(name))
                .and_then(|value| value.parse().ok())
                .unwrap_or_else(|| panic!("no {name} in: {line}"))
        };
        let builds = field("builds=");
        assert!(builds > 0, "{line}");
        assert_eq!(
            field("parallel="),
            if on_pool { builds } else { 0 },
            "{line}"
        );
        // The split: the pool fills' chunks and phases, or the original loop's time.
        let chunks: u64 = split_field(&stderr, "chunks=").parse().unwrap();
        assert_eq!(chunks > 0, on_pool, "{stderr}");
        let nonzero = |name: &str| split_field(&stderr, name) != "0.0000s";
        assert_eq!(nonzero("pool_fill="), on_pool, "{stderr}");
        assert_eq!(nonzero("read_wall="), on_pool, "{stderr}");
        assert_eq!(nonzero("loop="), !on_pool, "{stderr}");
        // The prefault rounds: at most one per pool fill, none off the pool; they touch whole
        // 2 MiB regions of the cells only (512 pages each, one zero per page).
        let rounds: u64 = split_field(&stderr, "x_prefault_rounds=").parse().unwrap();
        let touches: u64 = split_field(&stderr, "x_prefault_touches=").parse().unwrap();
        assert!(rounds <= builds, "{stderr}");
        assert_eq!(touches % 512, 0, "{stderr}");
        if !on_pool {
            assert_eq!((rounds, touches), (0, 0), "{stderr}");
        }
    }
    // Megablast builds no discontiguous table: no line.
    let (_, stderr) = blastn(
        &queries,
        &subject,
        &["-task", "megablast"],
        4,
        &[("LOSAT_TIMING", "1"), SWITCH],
    );
    assert!(!stderr.contains("disc_lookup_build"), "{stderr}");
}

// The value of `name` (`name=VALUE`) in the `[TIMING] disc_lookup_split` line of `stderr`.
fn split_field(stderr: &str, name: &str) -> String {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[TIMING] disc_lookup_split: "))
        .unwrap_or_else(|| panic!("no split line:\n{stderr}"));
    line.split(' ')
        .find_map(|word| word.strip_prefix(name))
        .unwrap_or_else(|| panic!("no {name} in: {line}"))
        .to_string()
}

// A table filled in several chunks of the pool fill: AP027131 (0.66 Mb, so about 1.3 million
// positions on both strands, against chunks of 524,288) against NZ_CP006932, dc-megablast with
// one template and with two (`coding_and_optimal`, 18): off at 1 and 4 threads, on and the shadow at 4
// threads give the same bytes, the shadow compares pool-built tables with the original loop's,
// and LOSAT_TIMING's split counts more chunks than tables.
#[test]
fn dc_megablast_across_chunk_boundaries() {
    let query = fasta("AP027131.fasta");
    let subject_path = fasta("NZ_CP006932.fasta");
    let subject = ["-subject", subject_path.to_str().unwrap()];
    for (extra, two_templates) in [
        (&["-task", "dc-megablast"][..], false),
        (
            &[
                "-task",
                "dc-megablast",
                "-template_type",
                "coding_and_optimal",
                "-template_length",
                "18",
            ][..],
            true,
        ),
    ] {
        let (expected, expected_stderr) = blastn(&query, &subject, extra, 1, &[OFF]);
        assert!(!expected.is_empty(), "{extra:?}: no hits");
        let (report, stderr) = blastn(&query, &subject, extra, 4, &[OFF]);
        assert!(report == expected, "{extra:?}: off differs");
        assert_eq!(stderr, expected_stderr, "{extra:?}");
        let (report, stderr) = blastn(&query, &subject, extra, 4, &[SWITCH]);
        assert!(report == expected, "{extra:?}: on differs");
        assert_eq!(stderr, expected_stderr, "{extra:?}");
        let (report, stderr) = blastn(&query, &subject, extra, 4, &[SHADOW]);
        assert!(report == expected, "{extra:?}: shadow differs");
        let (tables, pool, two) = shadow_counts(&stderr, &expected_stderr);
        assert!(tables > 0 && pool == tables, "{extra:?}: {stderr}");
        assert_eq!(two, if two_templates { tables } else { 0 }, "{stderr}");
        let (report, stderr) = blastn(&query, &subject, extra, 4, &[SWITCH, ("LOSAT_TIMING", "1")]);
        assert!(
            report == expected,
            "{extra:?}: on with LOSAT_TIMING differs"
        );
        let chunks: u64 = split_field(&stderr, "chunks=").parse().unwrap();
        assert!(
            chunks > tables,
            "{extra:?}: {chunks} chunks for {tables} tables:\n{stderr}"
        );
    }
}

// LOSAT_X_DCCSR on a table above its size rule: MG1655 (4.6 Mb, one table of about 9.2 million
// words, more than a quarter of the 4,194,304 cells at word size 11; AP027131's table above has
// 818,022 words and keeps its chains) against AP027131, dc-megablast with one template, laid out
// as CSR, built on the pool at 4 threads and by one loop at 1 thread. The chains (`=0`), the CSR
// (`=1`) and the shadow (both, compared) give the same bytes at 1 and 4 threads; the shadow's
// line counts CSR tables built on the pool only at 4 threads, and LOSAT_TIMING counts the CSR
// builds the same way.
#[test]
fn dc_megablast_csr_of_a_large_table() {
    const CSR: (&str, &str) = ("LOSAT_X_DCCSR", "1");
    let query = fasta("MG1655.fna");
    let subject_path = fasta("AP027131.fasta");
    let subject = ["-subject", subject_path.to_str().unwrap()];
    let extra = &["-task", "dc-megablast"][..];
    let (expected, expected_stderr) = blastn(&query, &subject, extra, 1, &[]);
    assert!(!expected.is_empty(), "no hits");
    let field = |line: &str, name: &str| -> u64 {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
    };
    for threads in [1, 4] {
        let (report, stderr) = blastn(&query, &subject, extra, threads, &[CSR]);
        assert!(report == expected, "threads {threads}: on differs");
        assert_eq!(stderr, expected_stderr, "threads {threads}");
        let (report, stderr) = blastn(
            &query,
            &subject,
            extra,
            threads,
            &[("LOSAT_X_DCCSRSHADOW", "1")],
        );
        assert!(report == expected, "threads {threads}: shadow differs");
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[X_SHADOW] DCCSR "))
            .unwrap_or_else(|| panic!("no shadow line:\n{stderr}"));
        assert!(line.ends_with(" diffs=0 identical"), "{line}");
        let rest: Vec<&str> = stderr.lines().filter(|l| *l != line).collect();
        assert_eq!(
            rest,
            expected_stderr.lines().collect::<Vec<_>>(),
            "{stderr}"
        );
        let tables = field(line, "tables=");
        assert!(tables > 0, "{line}");
        let pool = if threads > 1 { tables } else { 0 };
        assert_eq!(field(line, "pool="), pool, "threads {threads}: {line}");
        let (report, stderr) = blastn(
            &query,
            &subject,
            extra,
            threads,
            &[CSR, ("LOSAT_TIMING", "1")],
        );
        assert!(
            report == expected,
            "threads {threads}: on with LOSAT_TIMING differs"
        );
        let line = stderr
            .lines()
            .find(|line| line.starts_with("[TIMING] disc_lookup_csr: "))
            .unwrap_or_else(|| panic!("no CSR timing line:\n{stderr}"));
        let builds = field(line, "builds=");
        assert!(builds > 0, "{line}");
        let parallel = if threads > 1 { builds } else { 0 };
        assert_eq!(
            field(line, "parallel="),
            parallel,
            "threads {threads}: {line}"
        );
    }
}
