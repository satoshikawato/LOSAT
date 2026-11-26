//! `losat blastdbcheck` against NCBI BLAST+ 2.17.0 blastdbcheck: each row of
//! `tests/fixtures/blastdbcheck/cases.tsv` runs in its own fresh run directory (`in/` copied from
//! `tests/fixtures/blastdb_aliastool/in`, an empty `home/`, and the row's case directory built by
//! the operations of `setup.tsv`) and gives NCBI's standard output, standard error (in LOSAT's
//! words, `losat::wording`; the run's directory replaced by `<RUN>`), exit status
//! (`ncbi/<row>.{out,err,rc}`) and, for each file the run made or changed, its size and FNV-1a
//! digest (`ncbi/<row>.files`), frozen from the oracle of row RG-D2b (2026-11-25; README.md of the
//! directory). NCBI is not run here.
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1-1560
//! The rows whose output depends on the time-seeded random sample (`seed` not `-`) are skipped
//! here; a crate-internal test covers them.

#![cfg(unix)]

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// Rows where LOSAT must refuse (exit status 1, standard error starting `Error: ` and containing
/// `not supported by LOSAT`) because NCBI's result is a crash or an endless loop that LOSAT does
/// not reproduce: a sequence file of no bytes, an empty `.rps` profile file, a cut string ISAM, a
/// garbage index found by `-dir` (status 139), `-stride` of 0 or less (killed after 5 s); NCBI's
/// frozen result stays beside it.
const LOSAT_REFUSES: [&str; 8] = [
    "d_p4_psq_zero",
    "d_p5_psq_zero",
    "c_rps0",
    "d_p4_psd_cut",
    "r_garbage",
    "r_garbage2",
    "stride0",
    "strideneg",
];

/// Rows covered by approved exception 32 (`Q-D2b-OV=A`, `PD-LOSAT-NCBI-DEFECTS` version 1.114,
/// Owner decision 2026-11-25 17:54 JST): NCBI's 32/64-bit signed integer arithmetic overflows and
/// wraps to produce false alias errors (`a_ov1`, `a_ov2`, and `-dir .` runs `a_dir`, `a_dir_v1`
/// that visit them) or a wrapped OID (`ends_intmin`), whereas LOSAT widens the calculation and
/// reports the exact result (NCBI's frozen output stays beside every row for provenance).
const APPROVED_WIDENED_OVERFLOW: [&str; 5] = ["a_ov1", "a_ov2", "a_dir", "a_dir_v1", "ends_intmin"];

/// `tests/fixtures/blastdbcheck`.
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdbcheck")
}

/// The base databases, shared with the blastdb_aliastool fixtures.
fn base_in() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/blastdb_aliastool/in")
}

/// One row of `cases.tsv`.
struct Row {
    id: String,
    dir: String,
    env: Vec<(String, String)>,
    seed: String,
    args: String,
}

fn rows() -> Vec<Row> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let mut part = line.splitn(5, '\t');
            let mut next = || part.next().unwrap_or("").to_string();
            let (id, dir, env, seed, args) = (next(), next(), next(), next(), next());
            let env = if env == "-" {
                Vec::new()
            } else {
                env.split(' ')
                    .map(|pair| {
                        let (name, value) = pair.split_once('=').unwrap();
                        (name.to_string(), value.to_string())
                    })
                    .collect()
            };
            Row {
                id,
                dir,
                env,
                seed,
                args,
            }
        })
        .collect()
}

/// The words of `text` as bash splits them: spaces separate words, single quotes keep their text
/// (the table's only quoting).
fn words(text: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut word: Option<String> = None;
    let mut quoted = false;
    for c in text.chars() {
        match c {
            '\'' => {
                quoted = !quoted;
                word.get_or_insert_with(String::new);
            }
            ' ' if !quoted => {
                if let Some(word) = word.take() {
                    words.push(word);
                }
            }
            _ => word.get_or_insert_with(String::new).push(c),
        }
    }
    words.extend(word);
    words
}

/// FNV-1a 64 of `bytes`.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// `text` of a `text` operation with `\n`, `\r`, `\t`, `\\` and `<RUN>` replaced.
fn unescape(text: &str, run: &str) -> String {
    let mut out = String::new();
    let mut rest = text;
    while let Some(c) = rest.chars().next() {
        if let Some(after) = rest.strip_prefix("<RUN>") {
            out.push_str(run);
            rest = after;
        } else if c == '\\' && rest.len() > 1 && "nrt\\".contains(&rest[1..2]) {
            out.push(match &rest[1..2] {
                "n" => '\n',
                "r" => '\r',
                "t" => '\t',
                _ => '\\',
            });
            rest = &rest[2..];
        } else {
            out.push(c);
            rest = &rest[c.len_utf8()..];
        }
    }
    out
}

fn hex_bytes(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[2 * i..2 * i + 2], 16).unwrap())
        .collect()
}

fn set_mode_644(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(path, fs::Permissions::from_mode(0o644)).unwrap();
}

/// Builds `run/in`, `run/home` and `run/cases/<dir>` of `row` by the operations of `setup.tsv`.
fn build(run: &Path, row: &Row) {
    let run_text = run.to_str().unwrap();
    copy_tree(&base_in(), &run.join("in"));
    fs::create_dir_all(run.join("home")).unwrap();
    let dir = run.join("cases").join(&row.dir);
    fs::create_dir_all(&dir).unwrap();
    let files = fixtures().join("files");
    let table = fs::read_to_string(fixtures().join("setup.tsv")).unwrap();
    for line in table
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let part: Vec<&str> = line.split('\t').collect();
        if part[0] != row.dir {
            continue;
        }
        let (op, a) = (part[1], &part[2..]);
        match op {
            "db" => {
                let sub = a.get(1).copied().unwrap_or(".");
                let prefix = a.get(2).copied().unwrap_or("");
                let target = dir.join(sub);
                fs::create_dir_all(&target).unwrap();
                let stem = format!("{}.", a[0]);
                for entry in fs::read_dir(base_in()).unwrap() {
                    let name = entry.unwrap().file_name().to_string_lossy().into_owned();
                    if name.starts_with(&stem) {
                        let to = target.join(format!("{prefix}{name}"));
                        fs::copy(base_in().join(&name), &to).unwrap();
                        set_mode_644(&to);
                    }
                }
            }
            "tax" => {
                for name in ["taxdb.bti", "taxdb.btd"] {
                    fs::copy(files.join("tax").join(a[0]).join(name), dir.join(name)).unwrap();
                }
            }
            "put" => {
                let to = dir.join(a[0]);
                fs::create_dir_all(to.parent().unwrap()).unwrap();
                fs::copy(files.join(a[1]), &to).unwrap();
                set_mode_644(&to);
            }
            "cp" => {
                fs::copy(dir.join(a[0]), dir.join(a[1])).unwrap();
            }
            "rm" => fs::remove_file(dir.join(a[0])).unwrap(),
            "trunc" => {
                let path = dir.join(a[0]);
                let mut bytes = fs::read(&path).unwrap();
                let size = if a[1] == "half" {
                    bytes.len() / 2
                } else if let Some(drop) = a[1].strip_prefix('-') {
                    bytes.len() - drop.parse::<usize>().unwrap()
                } else {
                    a[1].parse().unwrap()
                };
                bytes.truncate(size);
                fs::write(path, bytes).unwrap();
            }
            "poke" => {
                let path = dir.join(a[0]);
                let mut bytes = fs::read(&path).unwrap();
                let at: usize = a[1].parse().unwrap();
                let value = hex_bytes(a[2]);
                bytes[at..at + value.len()].copy_from_slice(&value);
                fs::write(path, bytes).unwrap();
            }
            "text" => {
                let text = unescape(a.get(1).copied().unwrap_or(""), run_text);
                fs::write(dir.join(a[0]), text).unwrap();
            }
            "mkdir" => fs::create_dir_all(dir.join(a[0])).unwrap(),
            "link" => std::os::unix::fs::symlink(a[1], dir.join(a[0])).unwrap(),
            other => panic!("unknown setup operation {other}"),
        }
    }
}

/// Every file below `dir` with its bytes, by its `/`-separated path relative to `dir`
/// (links to directories are not followed).
fn listing(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut found = BTreeMap::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            let kind = entry.file_type().unwrap();
            if kind.is_dir() {
                stack.push(path);
            } else if kind.is_symlink() && path.is_dir() {
                continue;
            } else if let Ok(bytes) = fs::read(&path) {
                let relative = path.strip_prefix(dir).unwrap();
                let parts: Vec<String> = relative
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                found.insert(parts.join("/"), bytes);
            }
        }
    }
    found
}

/// What a run gave: stdout, stderr, exit status, and the files it made, changed or removed
/// (`path\tsize\tdigest`).
type Run = (Vec<u8>, Vec<u8>, i32, String);

/// `bytes` with the run's directory `dir` replaced by `<RUN>`.
fn with_run(bytes: &[u8], dir: &Path) -> Vec<u8> {
    let run = dir.as_os_str().as_encoded_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i..].starts_with(run) {
            out.extend_from_slice(b"<RUN>");
            i += run.len();
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    out
}

/// A command of `binary` with no LOSAT, NCBI or BLAST variable, run in `cwd` in UTC with `home` as
/// its home directory.
fn command(binary: &Path, cwd: &Path, home: &Path, row: &Row) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", home)
        .env("HOME", home)
        .env("TZ", "UTC")
        .env("LC_ALL", "C")
        .env("BLAST_USAGE_REPORT", "false");
    for (name, value) in &row.env {
        command.env(name, value);
    }
    command.current_dir(cwd);
    command
}

/// The run of `losat blastdbcheck` for `row` in the fresh run directory `dir`.
fn run(dir: &Path, row: &Row) -> Run {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    let dir = dir.canonicalize().unwrap();
    build(&dir, row);
    let before = listing(&dir);
    let cwd = dir.join("cases").join(&row.dir);
    let arguments: Vec<String> = words(&row.args)
        .into_iter()
        .map(|word| word.replace("<RUN>", dir.to_str().unwrap()))
        .collect();
    let io = dir.with_extension("io");
    fs::create_dir_all(&io).unwrap();
    let mut child = command(
        Path::new(env!("CARGO_BIN_EXE_losat")),
        &cwd,
        &dir.join("home"),
        row,
    )
    .arg("blastdbcheck")
    .args(&arguments)
    .stdin(Stdio::null())
    .stdout(fs::File::create(io.join("out")).unwrap())
    .stderr(fs::File::create(io.join("err")).unwrap())
    .spawn()
    .unwrap();
    let started = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status.code().unwrap_or(-1);
        }
        if started.elapsed() > Duration::from_secs(120) {
            let _ = child.kill();
            let _ = child.wait();
            break -2;
        }
        std::thread::sleep(Duration::from_millis(10));
    };
    let read = |name: &str| {
        let mut bytes = Vec::new();
        fs::File::open(io.join(name))
            .unwrap()
            .read_to_end(&mut bytes)
            .unwrap();
        bytes
    };
    let (stdout, stderr) = (read("out"), read("err"));
    let _ = fs::remove_dir_all(&io);
    let after = listing(&dir);
    let mut files = String::new();
    for (path, bytes) in &after {
        if before.get(path) == Some(bytes) {
            continue;
        }
        let bytes = with_run(bytes, &dir);
        files.push_str(&format!(
            "{path}\t{}\t{:016x}\n",
            bytes.len(),
            fnv1a64(&bytes)
        ));
    }
    for path in before.keys().filter(|path| !after.contains_key(*path)) {
        files.push_str(&format!("{path}\tREMOVED\t-\n"));
    }
    (
        with_run(&stdout, &dir),
        with_run(&stderr, &dir),
        status,
        files,
    )
}

/// NCBI's run of the row `name`, its standard output and error in LOSAT's words
/// (`losat::wording`). A file the run made whose bytes are frozen (`ncbi/<name>.made/<path>`, the
/// `-logfile` rows: NCBI's exception texts) is given by the size and digest of those bytes in
/// LOSAT's words.
fn expected(name: &str) -> Run {
    let file = |extension: &str| fixtures().join("ncbi").join(format!("{name}.{extension}"));
    let err = losat::wording::losat_bytes(&fs::read(file("err")).unwrap()).into_owned();
    let out = losat::wording::losat_output_bytes(&fs::read(file("out")).unwrap()).into_owned();
    let made = file("made");
    let mut files = String::new();
    for line in fs::read_to_string(file("files")).unwrap().lines() {
        let path = line.split('\t').next().unwrap();
        match fs::read(made.join(path)) {
            Ok(bytes) => {
                let bytes = losat::wording::losat_output_bytes(&bytes);
                files.push_str(&format!(
                    "{path}\t{}\t{:016x}\n",
                    bytes.len(),
                    fnv1a64(&bytes)
                ));
            }
            Err(_) => {
                files.push_str(line);
                files.push('\n');
            }
        }
    }
    (
        out,
        err,
        fs::read_to_string(file("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap(),
        files,
    )
}

fn scratch(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("blastdbcheck-{name}-{}", std::process::id()))
}

/// The text of `bytes` for a report, cut to 800 characters.
fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).chars().take(800).collect()
}

/// Every list of rows that LOSAT must refuse.
const REFUSAL_LISTS: [&[&str]; 1] = [&LOSAT_REFUSES];

fn must_refuse(name: &str) -> bool {
    REFUSAL_LISTS.iter().any(|list| list.contains(&name))
}

/// Expected LOSAT output for an approved widened-overflow row (`Q-D2b-OV=A`,
/// `PD-LOSAT-NCBI-DEFECTS` exception 32): exact arithmetic avoids NCBI's signed overflow false
/// errors (`a_ov1`, `a_ov2`, and `-dir .` runs `a_dir`/`a_dir_v1` that visit them) and negative
/// wraparound (`ends_intmin`), keeping every other line and warning identical to NCBI.
fn expected_widened(name: &str) -> Run {
    if matches!(name, "a_dir" | "a_dir_v1") {
        let (out, err, rc, files) = expected(name);
        let text = String::from_utf8(out).unwrap();
        let adjusted = text
            .replace(
                "  <RUN>/cases/a/ov1.pal / AliasFileTest:   [ERROR] FIRST_OID not paired with LAST_OID\n",
                "",
            )
            .replace(
                "  <RUN>/cases/a/ov2.pal / AliasFileTest:   [ERROR] OIDLIST cannot be provided without MAXOID , LENGTH, or NSEQ\n",
                "",
            )
            .replace(
                " Result=FAILURE. 21 errors reported in 26 alias(es).\n",
                " Result=FAILURE. 19 errors reported in 26 alias(es).\n",
            )
            .replace("Total errors: 39\n", "Total errors: 37\n");
        return (adjusted.into_bytes(), err, rc, files);
    }
    let alias_ok = |alias: &str, warn: &str| -> Vec<u8> {
        let mut s = String::from(concat!(
            "Writing messages to <stdout> at verbosity (Detailed)\n",
            "ISAM testing is ENABLED.\n",
            "Legacy testing is DISABLED.\n",
            "TaxID testing is DISABLED.\n",
            "Testing every 100-th OID.\n",
            "\n",
            "Testing 1 volume(s).\n",
            " <RUN>/cases/a/p5\n",
            " Result=SUCCESS. No errors reported for 1 volume(s).\n",
            "Testing 1 alias(es).\n",
        ));
        s.push_str(&format!(" <RUN>/cases/a/{alias}\n{warn}"));
        s.push_str(" Result=SUCCESS. No errors reported for 1 alias(es).\n");
        s.into_bytes()
    };
    let (out, rc): (Vec<u8>, i32) = match name {
        "a_ov1" => (alias_ok("ov1.pal", ""), 0),
        "a_ov2" => (
            alias_ok(
                "ov2.pal",
                "  <RUN>/cases/a/ov2.pal / AliasFileTest:   [WARNING] oidlist file gi.gil does not have .msk extension\n",
            ),
            0,
        ),
        "ends_intmin" => (
            concat!(
                "Writing messages to <stdout> at verbosity (Minutiae)\n",
                "ISAM testing is ENABLED.\n",
                "Legacy testing is DISABLED.\n",
                "TaxID testing is DISABLED.\n",
                "Testing first -2147483648 and last -2147483648 OIDs.\n",
                "\n",
                "Testing 1 volume(s).\n",
                " <RUN>/in/p5\n",
                "  <RUN>/in/p5 / EndPoints: <testing -2147483648 OIDs at each end>\n",
                "  <RUN>/in/p5 / EndPoints:     Failed during headers, oid 2147483657: OID not found\n",
                "\n",
                "  <RUN>/in/p5 / EndPoints:       Status for OID 2147483657: FAIL\n",
                " Result=FAILURE. 1 errors reported in 1 volume(s).\n",
                "Testing 0 alias(es).\n",
                " Result=SUCCESS. No errors reported for 0 alias(es).\n",
                "\n",
                "Total errors: 1\n",
            )
            .as_bytes()
            .to_vec(),
            1,
        ),
        _ => panic!("unexpected approved widened row {name}"),
    };
    (out, Vec::new(), rc, String::new())
}

#[test]
fn blastdbcheck_gives_ncbis_bytes() {
    let mut differ = Vec::new();
    for row in rows() {
        if row.seed != "-" {
            continue;
        }
        let name = row.id.clone();
        let dir = scratch(&name);
        let got = run(&dir, &row);
        if must_refuse(&name) {
            let refused = got.2 == 1
                && got.1.starts_with(b"Error: ")
                && text(&got.1).contains("not supported by LOSAT");
            if !refused {
                eprintln!(
                    "--- {name}: expected LOSAT's refusal, got status {}\nstderr: {:?}\nfiles:\n{}",
                    got.2,
                    text(&got.1),
                    got.3
                );
                differ.push(name);
            } else {
                let _ = fs::remove_dir_all(&dir);
            }
            continue;
        }
        let want = if APPROVED_WIDENED_OVERFLOW.contains(&name.as_str()) {
            expected_widened(&name)
        } else {
            expected(&name)
        };
        if got != want {
            eprintln!(
                "--- {name}: LOSAT status {}, expected {}\nLOSAT stdout: {:?}\nexpected stdout:  {:?}\n\
                 LOSAT stderr: {:?}\nexpected stderr:  {:?}\nLOSAT files:\n{}expected files:\n{}",
                got.2,
                want.2,
                text(&got.0),
                text(&want.0),
                text(&got.1),
                text(&want.1),
                got.3,
                want.3,
            );
            differ.push(name);
        } else {
            let _ = fs::remove_dir_all(&dir);
        }
    }
    assert!(
        differ.is_empty(),
        "rows that differ from expected: {differ:?}"
    );
}

/// Every id in the refusal and approved-exception lists names a row of the table, and no row is in
/// two lists.
#[test]
fn refusal_lists_name_rows() {
    let ids: Vec<String> = rows().into_iter().map(|row| row.id).collect();
    let mut seen = Vec::new();
    for name in REFUSAL_LISTS
        .iter()
        .flat_map(|list| list.iter())
        .chain(APPROVED_WIDENED_OVERFLOW.iter())
    {
        assert!(ids.iter().any(|id| id == name), "no row {name}");
        assert!(!seen.contains(name), "{name} is in two lists");
        seen.push(name);
    }
}

/// D2b-3 cross-program parity: `blastdbcmd`, `blastp`, and `blastn` share the same 4 KiB
/// last-page zero-read behavior in `seqdb` and match NCBI BLAST+ 2.17.0 on the same truncated
/// inputs (`p4_pin_10`, `p4_pin_cut`, `p4_psq_half`, `p4_phr_half`, `n4_nsq_half`).
#[test]
fn zero_page_cross_program_parity() {
    let run_prog = |dir: &Path, sub: &str, args: &[&str]| -> (Vec<u8>, Vec<u8>, i32) {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_losat"));
        cmd.arg(sub)
            .args(args)
            .current_dir(dir)
            .stdin(Stdio::null())
            .env_clear()
            .env("PATH", "/usr/bin:/bin")
            .env("TZ", "UTC")
            .env("LC_ALL", "C");
        let out = cmd.output().unwrap();
        (out.stdout, out.stderr, out.status.code().unwrap_or(-1))
    };
    let cp_db = |dst: &Path, stem: &str| {
        fs::create_dir_all(dst).unwrap();
        let prefix = format!("{stem}.");
        for entry in fs::read_dir(base_in()).unwrap() {
            let entry = entry.unwrap();
            let name = entry.file_name().into_string().unwrap();
            if name.starts_with(&prefix) {
                fs::copy(entry.path(), dst.join(name)).unwrap();
            }
        }
    };

    let root = scratch("zp-xprog");
    fs::create_dir_all(&root).unwrap();
    let q_prot = root.join("q.faa");
    fs::write(&q_prot, b">q\nCSNFRHATLYNGCRFTMHYLCSNFRHATLYNGCRFTMHYL\n").unwrap();
    let q_nucl = root.join("q.fna");
    fs::write(
        &q_nucl,
        b">qn\nATGAAACCCTTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTT\n",
    )
    .unwrap();
    let q_prot_s = q_prot.to_str().unwrap();
    let q_nucl_s = q_nucl.to_str().unwrap();
    let blastp_args = ["-query", q_prot_s, "-db", "p4", "-outfmt", "6"];
    let prot_all_args = [
        "-db", "p4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %l %t",
    ];
    let run_ok = |dir: &Path, sub: &str, args: &[&str]| -> Vec<u8> {
        let (out, err, rc) = run_prog(dir, sub, args);
        assert_eq!((rc, err.as_slice()), (0, b"".as_slice()));
        out
    };
    let run_bp = |dir: &Path| run_prog(dir, "blastp", &blastp_args);
    let run_bp_ok = |dir: &Path| run_ok(dir, "blastp", &blastp_args);

    // 1. p4_pin_10: blastdbcmd -info sees 0 sequences; blastp exits 3 with memory map error.
    let d = root.join("p4_pin10");
    cp_db(&d, "p4");
    fs::OpenOptions::new()
        .write(true)
        .open(d.join("p4.pin"))
        .unwrap()
        .set_len(10)
        .unwrap();
    let out = run_ok(&d, "blastdbcmd", &["-db", "p4", "-dbtype", "prot", "-info"]);
    assert!(
        text(&out).contains("\t0 sequences; 0 total residues\n"),
        "{}",
        text(&out)
    );
    let (out, err, rc) = run_bp(&d);
    assert_eq!((rc, out.as_slice()), (3, b"".as_slice()));
    assert_eq!(err, b"LOSAT engine error: Database memory map file error\n");

    // 2. p4_pin_cut (-8 bytes): blastdbcmd prints signed negative lengths -364 and -1; blastp hits XP_000104.1.
    let d = root.join("p4_pincut");
    cp_db(&d, "p4");
    let pin_len = fs::metadata(d.join("p4.pin")).unwrap().len();
    fs::OpenOptions::new()
        .write(true)
        .open(d.join("p4.pin"))
        .unwrap()
        .set_len(pin_len - 8)
        .unwrap();
    let out = run_ok(&d, "blastdbcmd", &prot_all_args);
    assert_eq!(
        text(&out),
        concat!(
            "0 60 alpha protein\n",
            "1 45 beta predicted protein\n",
            "2 70 gamma protein\n",
            "3 50 delta predicted protein\n",
            "4 40 epsilon protein\n",
            "5 55 zeta protein\n",
            "6 35 pir protein\n",
            "7 -364 prf protein\n",
            "8 -1 local protein\n",
        )
    );
    let out = run_bp_ok(&d);
    assert_eq!(
        text(&out),
        "q\tXP_000104.1\t100.000\t40\t0\t0\t1\t40\t1\t40\t2.85e-28\t82.8\n\
         q\tXP_000104.1\t100.000\t30\t0\t0\t1\t30\t21\t50\t2.81e-20\t62.4\n"
    );

    // 3. p4_psq_half: zero-filled trailing protein bytes cut the alignment at residue 31 in blastp.
    let d = root.join("p4_psqhalf");
    cp_db(&d, "p4");
    let psq_len = fs::metadata(d.join("p4.psq")).unwrap().len();
    fs::OpenOptions::new()
        .write(true)
        .open(d.join("p4.psq"))
        .unwrap()
        .set_len(psq_len / 2)
        .unwrap();
    let out = run_bp_ok(&d);
    assert_eq!(
        text(&out),
        "q\tXP_000104.1\t100.000\t31\t0\t0\t1\t31\t1\t31\t2.62e-21\t65.1\n\
         q\tXP_000104.1\t100.000\t20\t0\t0\t21\t40\t1\t20\t5.39e-13\t44.3\n"
    );

    // 4. p4_phr_half: blastdbcmd outputs OIDs 0..3, posts 8 Error (803.7) lines for NUL bytes,
    // then stops with `byte 25: member seqid expected`; blastp exits 255 on zero-filled header.
    let d = root.join("p4_phrhalf");
    cp_db(&d, "p4");
    let phr_len = fs::metadata(d.join("p4.phr")).unwrap().len();
    fs::OpenOptions::new()
        .write(true)
        .open(d.join("p4.phr"))
        .unwrap()
        .set_len(phr_len / 2)
        .unwrap();
    let (out, err, rc) = run_prog(&d, "blastdbcmd", &prot_all_args);
    assert_eq!(rc, 1);
    assert_eq!(
        text(&out),
        "0 60 alpha protein\n1 45 beta predicted protein\n2 70 gamma protein\n3 50 delta predicted protein\n"
    );
    let err_s = String::from_utf8_lossy(&err);
    assert_eq!(
        err_s
            .matches("Error: (803.7) [blastdbcmd] Blast-def-line-set.E.title\n")
            .count(),
        8,
        "{err_s}"
    );
    assert!(
        err_s.ends_with("Error: [blastdbcmd] byte 25: member seqid expected\n"),
        "{err_s}"
    );
    let (out, err, rc) = run_bp(&d);
    assert_eq!((rc, out.as_slice()), (255, b"".as_slice()));
    assert_eq!(
        text(&err),
        "Error: byte 0: unexpected tag: None (0), should be constructed/None (32) ( at Blast-def-line-set.[])\n\n"
    );

    // 5. n4_nsq_half: zero-filled trailing 2-bit bytes yield lengths 120, 90, 148, 80 and blastn rc 0.
    let d = root.join("n4_nsqhalf");
    cp_db(&d, "n4");
    let nsq_len = fs::metadata(d.join("n4.nsq")).unwrap().len();
    fs::OpenOptions::new()
        .write(true)
        .open(d.join("n4.nsq"))
        .unwrap()
        .set_len(nsq_len / 2)
        .unwrap();
    let out = run_ok(
        &d,
        "blastdbcmd",
        &[
            "-db", "n4", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %l %t",
        ],
    );
    assert_eq!(
        text(&out),
        "0 120 n one\n1 90 n two predicted\n2 148 n three\n3 80 n four\n"
    );
    let out = run_ok(
        &d,
        "blastn",
        &["-query", q_nucl_s, "-db", "n4", "-outfmt", "6"],
    );
    assert!(out.is_empty(), "{}", text(&out));

    let _ = fs::remove_dir_all(&root);
}

/// Loads a frozen NCBI BLAST+ 2.17.0 cross-program fixture (`nonasc_xprog/<name>.{out,err,rc}`)
/// translated into LOSAT wording.
fn load_nonasc_expected(name: &str) -> (Vec<u8>, Vec<u8>, i32) {
    let base = fixtures().join("nonasc_xprog");
    let out =
        losat::wording::losat_output_bytes(&fs::read(base.join(format!("{name}.out"))).unwrap())
            .into_owned();
    let err = losat::wording::losat_bytes(&fs::read(base.join(format!("{name}.err"))).unwrap())
        .into_owned();
    let rc: i32 = fs::read_to_string(base.join(format!("{name}.rc")))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    (out, err, rc)
}

/// Runs a `losat <subcommand>` with clean environment and `<RUN>` substitution against `root`.
fn run_losat_sub(root: &Path, dir: &Path, sub: &str, args: &[&str]) -> (Vec<u8>, Vec<u8>, i32) {
    let home = root.join("home");
    fs::create_dir_all(&home).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg(sub)
        .args(args)
        .current_dir(dir)
        .stdin(Stdio::null())
        .env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("TZ", "UTC")
        .env("LC_ALL", "C")
        .env("HOME", &home)
        .env("NCBI_CONFIG_PATH", &home)
        .env("BLAST_USAGE_REPORT", "false")
        .output()
        .unwrap();
    (
        with_run(&out.stdout, root),
        with_run(&out.stderr, root),
        out.status.code().unwrap_or(-1),
    )
}

/// PR #1660 R2 cross-program parity for non-ascending index offsets (`d_p4_nonasc_seq`,
/// `d_p4_nonasc_hdr`, `n4_amb_eq_start`, `n4_amb_gt_next`) across `blastdbcmd`, `blastp`
/// (`-num_threads 1` and `4` in `[profile.test]` with overflow checks), `blastn`, and
/// `blastdbcheck` against frozen NCBI BLAST+ 2.17.0 outputs in `nonasc_xprog/`.
#[test]
fn nonasc_cross_program_parity() {
    let root = scratch("nonasc-xprog");
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let root = root.canonicalize().unwrap();

    let cp_db = |dst: &Path, stem: &str, with_tax: bool| {
        fs::create_dir_all(dst).unwrap();
        let prefix = format!("{stem}.");
        for entry in fs::read_dir(base_in()).unwrap() {
            let entry = entry.unwrap();
            let name = entry.file_name().into_string().unwrap();
            if name.starts_with(&prefix) {
                fs::copy(entry.path(), dst.join(name)).unwrap();
            }
        }
        if with_tax {
            let tax = fixtures().join("files/tax/good");
            for f in ["taxdb.bti", "taxdb.btd"] {
                fs::copy(tax.join(f), dst.join(f)).unwrap();
            }
        }
    };
    let poke = |path: PathBuf, at: usize, hex: &str| {
        let mut b = fs::read(&path).unwrap();
        let v = hex_bytes(hex);
        b[at..at + v.len()].copy_from_slice(&v);
        fs::write(path, b).unwrap();
    };
    let assert_refused = |label: &str, got: (Vec<u8>, Vec<u8>, i32)| {
        assert_eq!(got.2, 1, "{label}: expected rc=1, err={:?}", text(&got.1));
        assert!(
            got.1.starts_with(b"Error: ") && text(&got.1).contains("not supported by LOSAT"),
            "{label}: expected LOSAT refusal, got stderr={:?}",
            text(&got.1)
        );
    };

    let q3 = root.join("q3.faa");
    fs::write(&q3, b">q3\nCSNFRHATLYNGCRFTMHYLCSNFRHATLYNGCRFTMHYL\n").unwrap();
    let q2 = root.join("q2.faa");
    fs::write(&q2, b">q2\nRKFVIYQLDPFWRIVLEYPDRKFVIYQLDPFWRIVLEYPD\n").unwrap();
    let q4 = root.join("q4.faa");
    fs::write(&q4, b">q4\nKDWNCQIETHWPKCNEVQHTKDWNCQIETHWPKCNEVQHT\n").unwrap();
    let qn = root.join("qn.fna");
    fs::write(
        &qn,
        b">qn\nAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTG\n",
    )
    .unwrap();
    let (q3s, q2s, q4s, qns) = (
        q3.to_str().unwrap(),
        q2.to_str().unwrap(),
        q4.to_str().unwrap(),
        qn.to_str().unwrap(),
    );

    // 1. d_p4_nonasc_seq (p4.pin@106 = 01f4 -> seq=[1, 62, 500, 179, ...]).
    let d_seq = root.join("d_p4_nonasc_seq");
    cp_db(&d_seq, "p4", true);
    poke(d_seq.join("p4.pin"), 106, "01f4");
    for (name, sub, args) in [
        (
            "nonasc_seq_dbcmd_fasta",
            "blastdbcmd",
            vec!["-db", "p4", "-dbtype", "prot", "-entry", "all"],
        ),
        (
            "nonasc_seq_dbcmd_fmt",
            "blastdbcmd",
            vec![
                "-db", "p4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %l %t",
            ],
        ),
        (
            "nonasc_seq_blastp_t1",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ],
        ),
    ] {
        let got = run_losat_sub(&root, &d_seq, sub, &args);
        let want = load_nonasc_expected(name);
        assert_eq!(got, want, "{name}");
    }
    // Cross-thread invariance + debug-build overflow-check coverage: -num_threads 4 matches -num_threads 1
    // (NCBI -num_threads 4 drops OIDs 3..8 due to CSeqDBImpl::x_FillSeqBuffer returning empty chunk on len<0).
    let got_t4 = run_losat_sub(
        &root,
        &d_seq,
        "blastp",
        &[
            "-query",
            q3s,
            "-db",
            "p4",
            "-outfmt",
            "6",
            "-num_threads",
            "4",
        ],
    );
    assert_eq!(got_t4, load_nonasc_expected("nonasc_seq_blastp_t1"));

    // 2. d_p4_nonasc_hdr (p4.pin@68 = 00000000 -> hdr=[0, 79, 167, 0, 334, ...]).
    let d_hdr = root.join("d_p4_nonasc_hdr");
    cp_db(&d_hdr, "p4", true);
    poke(d_hdr.join("p4.pin"), 68, "00000000");
    for (name, sub, args) in [
        (
            "nonasc_hdr_dbcmd_fasta",
            "blastdbcmd",
            vec!["-db", "p4", "-dbtype", "prot", "-entry", "all"],
        ),
        (
            "nonasc_hdr_dbcmd_fmt",
            "blastdbcmd",
            vec![
                "-db", "p4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %l %t",
            ],
        ),
        (
            "nonasc_hdr_blastp_q2_t1",
            "blastp",
            vec![
                "-query",
                q2s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ],
        ),
        (
            "nonasc_hdr_blastp_q2_t1",
            "blastp",
            vec![
                "-query",
                q2s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_hdr_blastp_q4_t1",
            "blastp",
            vec![
                "-query",
                q4s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ],
        ),
        (
            "nonasc_hdr_blastp_q4_t1",
            "blastp",
            vec![
                "-query",
                q4s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
    ] {
        let got = run_losat_sub(&root, &d_hdr, sub, &args);
        let want = load_nonasc_expected(name);
        assert_eq!(got, want, "{name}");
    }
    // OID 3's header range [0, 334) decodes OID 0's gi|101 -> existing lookup-mismatch refusal.
    for t in ["1", "4"] {
        let got = run_losat_sub(
            &root,
            &d_hdr,
            "blastp",
            &[
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                t,
            ],
        );
        assert_refused(&format!("nonasc_hdr_blastp_q3_t{t}"), got);
    }

    // 3. Nucleotide ambiguity offset out of order: amb == start (OOB buffer[-1]) and amb > next
    // (reads past OID into OID 2; SIGABRT in NCBI blastn -num_threads 4). Uniform Rule-3 refusal.
    for (label, hex) in [
        ("n4_amb_eq_start", "00000020"),
        ("n4_amb_gt_next", "0000003b"),
    ] {
        let d = root.join(label);
        cp_db(&d, "n4", true);
        poke(d.join("n4.nin"), 100, hex);
        for (sub, args) in [
            (
                "blastdbcmd",
                vec!["-db", "n4", "-dbtype", "nucl", "-entry", "all"],
            ),
            (
                "blastn",
                vec![
                    "-query",
                    qns,
                    "-db",
                    "n4",
                    "-outfmt",
                    "6",
                    "-num_threads",
                    "1",
                ],
            ),
            (
                "blastn",
                vec![
                    "-query",
                    qns,
                    "-db",
                    "n4",
                    "-outfmt",
                    "6",
                    "-num_threads",
                    "4",
                ],
            ),
            (
                "blastdbcheck",
                vec!["-db", "n4", "-dbtype", "nucl", "-full"],
            ),
        ] {
            let got = run_losat_sub(&root, &d, sub, &args);
            assert_refused(&format!("{label}_{sub}"), got);
        }
    }

    let _ = fs::remove_dir_all(&root);
}
