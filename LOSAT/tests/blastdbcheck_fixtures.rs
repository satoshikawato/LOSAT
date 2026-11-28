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

/// Rows covered by approved exception 35 (`Q-D2b-OV=A`, `PD-LOSAT-NCBI-DEFECTS` version 1.124,
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
/// `PD-LOSAT-NCBI-DEFECTS` exception 35): exact arithmetic avoids NCBI's signed overflow false
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

const PROT_ALL_FASTA: [&str; 6] = ["-db", "p4", "-dbtype", "prot", "-entry", "all"];
const PROT_ALL_FMT: [&str; 8] = [
    "-db", "p4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %l %t",
];

fn copy_base_db(dst: &Path, stem: &str) {
    fs::create_dir_all(dst).unwrap();
    let prefix = format!("{stem}.");
    for entry in fs::read_dir(base_in()).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name().into_string().unwrap();
        if name.starts_with(&prefix) {
            fs::copy(entry.path(), dst.join(name)).unwrap();
        }
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
    let cp_db = copy_base_db;

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
    let prot_all_args = PROT_ALL_FMT;
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
    let out = run_ok(&d, "blastdbcmd", N4_DBCMD_FMT);
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

const N4_DBCMD: &[&str] = &["-db", "n4", "-dbtype", "nucl", "-entry", "all"];
const N4_DBCMD_FMT: &[&str] = &[
    "-db", "n4", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %l %t",
];
const N4_CHECK: &[&str] = &["-db", "n4", "-dbtype", "nucl", "-full"];
const N4_CHECK_V3: &[&str] = &["-db", "n4", "-dbtype", "nucl", "-full", "-verbosity", "3"];

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

/// Runs a `losat <subcommand>` with clean environment and `<RUN>` substitution against `root`,
/// optionally wrapped by an executable prefix (such as `["taskset", "-c", "0"]`).
fn run_losat_sub_prefixed(
    root: &Path,
    dir: &Path,
    prefix: &[&str],
    sub: &str,
    args: &[&str],
) -> (Vec<u8>, Vec<u8>, i32) {
    let home = root.join("home");
    fs::create_dir_all(&home).unwrap();
    let mut sub_cmd = if let Some((prog, rest)) = prefix.split_first() {
        let mut cmd = Command::new(prog);
        cmd.args(rest).arg(env!("CARGO_BIN_EXE_losat"));
        cmd
    } else {
        Command::new(env!("CARGO_BIN_EXE_losat"))
    };
    sub_cmd
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
        .env("BLAST_USAGE_REPORT", "false");
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        sub_cmd.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    let out = sub_cmd.output().unwrap();
    (
        with_run(&out.stdout, root),
        with_run(&out.stderr, root),
        out.status.code().unwrap_or(-1),
    )
}

/// Runs a `losat <subcommand>` with clean environment and `<RUN>` substitution against `root`.
fn run_losat_sub(root: &Path, dir: &Path, sub: &str, args: &[&str]) -> (Vec<u8>, Vec<u8>, i32) {
    run_losat_sub_prefixed(root, dir, &[], sub, args)
}

/// Cross-program parity for non-ascending index offsets (`d_p4_nonasc_seq`, `d_p4_nonasc_hdr`,
/// `n4_amb_eq_start`, `n4_amb_gt_next`, `n4_amb0_eq_0`) across `blastdbcmd`, `blastp`, `blastx`,
/// `psiblast`, `blastn`, and `blastdbcheck` (under `[profile.test]` with overflow checks enabled)
/// against frozen NCBI BLAST+ 2.17.0 outputs in `nonasc_xprog/`.
#[test]
fn nonasc_cross_program_parity() {
    let root = scratch("nonasc-xprog");
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let root = root.canonicalize().unwrap();

    let cp_db = |dst: &Path, stem: &str, with_tax: bool| {
        copy_base_db(dst, stem);
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
    let q3_dna = root.join("q3_dna.fna");
    fs::write(
        &q3_dna,
        b">q3_dna\n\
          TGCTCCAACTTCCGCCACGCAACGCTGTACAATGGCTGCCGCTTCACGATGCACTACCTG\
          TGCTCCAACTTCCGCCACGCAACGCTGTACAATGGCTGCCGCTTCACGATGCACTACCTG\n",
    )
    .unwrap();
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
    let (q3s, q3ds, q2s, q4s, qns) = (
        q3.to_str().unwrap(),
        q3_dna.to_str().unwrap(),
        q2.to_str().unwrap(),
        q4.to_str().unwrap(),
        qn.to_str().unwrap(),
    );

    let neg_gi = root.join("neg_gi.txt");
    fs::write(&neg_gi, b"103\n").unwrap();
    let pos_gi = root.join("pos_gi.txt");
    fs::write(&pos_gi, b"101\n102\n104\n105\n106\n").unwrap();
    let neg_si = root.join("neg_si.txt");
    fs::write(&neg_si, b"AAB00103.2\n").unwrap();
    let n4_neg_gi = root.join("n4_neg_gi.txt");
    fs::write(&n4_neg_gi, b"302\n").unwrap();
    let neg_gi2 = root.join("neg_gi2.txt");
    fs::write(&neg_gi2, b"103\n105\n").unwrap();
    let n4_neg_gi304 = root.join("n4_neg_gi304.txt");
    fs::write(&n4_neg_gi304, b"304\n").unwrap();
    let n4_neg_gi_all = root.join("n4_neg_gi_all.txt");
    fs::write(&n4_neg_gi_all, b"301\n302\n303\n304\n").unwrap();
    let (neg_gis, pos_gis, neg_sis, n4_neg_gis) = (
        neg_gi.to_str().unwrap(),
        pos_gi.to_str().unwrap(),
        neg_si.to_str().unwrap(),
        n4_neg_gi.to_str().unwrap(),
    );
    let (neg_gi2s, n4_neg_gi304s, n4_neg_gi_alls) = (
        neg_gi2.to_str().unwrap(),
        n4_neg_gi304.to_str().unwrap(),
        n4_neg_gi_all.to_str().unwrap(),
    );

    // 1. d_p4_nonasc_seq (p4.pin@106 = 01f4 -> seq=[1, 62, 500, 179, ...]).
    // At 2+ threads in DB-split mode, NCBI's CSeqDBImpl::x_FillSeqBuffer stops at OID 2 (len=-322)
    // before applying any GI or seqid list; at 1 thread or `-mt_mode 1` (split by queries), NCBI
    // scans all included OIDs.
    let d_seq = root.join("d_p4_nonasc_seq");
    cp_db(&d_seq, "p4", true);
    poke(d_seq.join("p4.pin"), 106, "01f4");
    for (name, sub, args) in [
        (
            "nonasc_seq_dbcmd_fasta",
            "blastdbcmd",
            PROT_ALL_FASTA.to_vec(),
        ),
        ("nonasc_seq_dbcmd_fmt", "blastdbcmd", PROT_ALL_FMT.to_vec()),
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
        (
            "nonasc_seq_blastp_t2",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "2",
            ],
        ),
        (
            "nonasc_seq_blastp_t4",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastp_mt1_t4",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-mt_mode",
                "1",
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastp_neggi_t1",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gis,
                "-num_threads",
                "1",
            ],
        ),
        (
            "nonasc_seq_blastp_neggi_t4",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gis,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastp_gi_t4",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-gilist",
                pos_gis,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastp_negsi_t4",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_seqidlist",
                neg_sis,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastp_neggi2_t4",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gi2s,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastp_neggi2_t1",
            "blastp",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gi2s,
                "-num_threads",
                "1",
            ],
        ),
        (
            "nonasc_seq_blastx_t1",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ],
        ),
        (
            "nonasc_seq_blastx_t4",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastx_mt1_t4",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-mt_mode",
                "1",
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastx_neggi_t1",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gis,
                "-num_threads",
                "1",
            ],
        ),
        (
            "nonasc_seq_blastx_neggi_t4",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gis,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastx_gi_t4",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-gilist",
                pos_gis,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastx_negsi_t4",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_seqidlist",
                neg_sis,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_blastx_neggi2_t4",
            "blastx",
            vec![
                "-query",
                q3ds,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gi2s,
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_psiblast_t1",
            "psiblast",
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
        (
            "nonasc_seq_psiblast_t4",
            "psiblast",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
        (
            "nonasc_seq_psiblast_neggi2_t4",
            "psiblast",
            vec![
                "-query",
                q3s,
                "-db",
                "p4",
                "-outfmt",
                "6",
                "-negative_gilist",
                neg_gi2s,
                "-num_threads",
                "4",
            ],
        ),
    ] {
        let got = run_losat_sub(&root, &d_seq, sub, &args);
        let want = load_nonasc_expected(name);
        assert_eq!(got, want, "{name}");
    }
    // Under Linux CPU affinity (`taskset -c 0`), NCBI's `CSystemInfo::GetCpuCount` still counts
    // online CPUs (`_SC_NPROCESSORS_ONLN`) and keeps the 4-thread DB-split stop on both unfiltered
    // and list-filtered searches.
    #[cfg(target_os = "linux")]
    if Command::new("taskset")
        .args(["-c", "0", "true"])
        .status()
        .is_ok_and(|s| s.success())
    {
        let base = [
            "-query",
            q3s,
            "-db",
            "p4",
            "-outfmt",
            "6",
            "-num_threads",
            "4",
        ];
        for (name, extra) in [
            ("nonasc_seq_blastp_t4", &[][..]),
            ("nonasc_seq_blastp_neggi_t4", &["-negative_gilist", neg_gis]),
        ] {
            let mut args = base.to_vec();
            args.extend_from_slice(extra);
            let got =
                run_losat_sub_prefixed(&root, &d_seq, &["taskset", "-c", "0"], "blastp", &args);
            let want = load_nonasc_expected(name);
            assert_eq!(got, want, "{name} under taskset -c 0");
        }
    }

    // 2. d_p4_nonasc_hdr (p4.pin@68 = 00000000 -> hdr=[0, 79, 167, 0, 334, ...]).
    let d_hdr = root.join("d_p4_nonasc_hdr");
    cp_db(&d_hdr, "p4", true);
    poke(d_hdr.join("p4.pin"), 68, "00000000");
    for (name, sub, args) in [
        (
            "nonasc_hdr_dbcmd_fasta",
            "blastdbcmd",
            PROT_ALL_FASTA.to_vec(),
        ),
        ("nonasc_hdr_dbcmd_fmt", "blastdbcmd", PROT_ALL_FMT.to_vec()),
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

    // 3. n4_amb_eq_start (n4.nin@100 = 00000020 -> amb[1] = seq[1] = 32; length = -4 inside mapping).
    // Deterministic in blastn, blastdbcmd, and blastdbcheck; blastn -t 4 stops at OID 1 via
    // x_FillSeqBuffer (even when a GI list excludes OID 1). tblastn and tblastx refuse out-of-order
    // ambiguity offsets (documented in refusals.md).
    let d_ambeq = root.join("n4_amb_eq_start");
    cp_db(&d_ambeq, "n4", true);
    poke(d_ambeq.join("n4.nin"), 100, "00000020");
    for (name, sub, args) in [
        ("n4_ambeq_dbcmd", "blastdbcmd", N4_DBCMD.to_vec()),
        ("n4_ambeq_dbcmd_fmt", "blastdbcmd", N4_DBCMD_FMT.to_vec()),
        (
            "n4_ambeq_blastn_t1",
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
            "n4_ambeq_blastn_t4",
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
            "n4_ambeq_blastn_mt1_t4",
            "blastn",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-mt_mode",
                "1",
                "-num_threads",
                "4",
            ],
        ),
        (
            "n4_ambeq_blastn_neggi_t1",
            "blastn",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gis,
                "-num_threads",
                "1",
            ],
        ),
        (
            "n4_ambeq_blastn_neggi_t4",
            "blastn",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gis,
                "-num_threads",
                "4",
            ],
        ),
        ("n4_ambeq_check", "blastdbcheck", N4_CHECK.to_vec()),
        ("n4_ambeq_check_v3", "blastdbcheck", N4_CHECK_V3.to_vec()),
    ] {
        let got = run_losat_sub(&root, &d_ambeq, sub, &args);
        let want = load_nonasc_expected(name);
        assert_eq!(got, want, "{name}");
    }

    // 4. n4_amb_gt_next (n4.nin@100 = 0000003b -> amb[1] = 59 > next = 55; length = 107).
    // Deterministic in blastdbcmd, blastn -t 1, and blastdbcheck; blastn at 2+ threads aborts in
    // NCBI (SIGABRT) and is refused.
    let d_ambgt = root.join("n4_amb_gt_next");
    cp_db(&d_ambgt, "n4", true);
    poke(d_ambgt.join("n4.nin"), 100, "0000003b");
    for (name, sub, args) in [
        ("n4_ambgt_dbcmd", "blastdbcmd", N4_DBCMD.to_vec()),
        ("n4_ambgt_dbcmd_fmt", "blastdbcmd", N4_DBCMD_FMT.to_vec()),
        (
            "n4_ambgt_blastn_t1",
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
        ("n4_ambgt_check", "blastdbcheck", N4_CHECK.to_vec()),
        ("n4_ambgt_check_v3", "blastdbcheck", N4_CHECK_V3.to_vec()),
    ] {
        let got = run_losat_sub(&root, &d_ambgt, sub, &args);
        let want = load_nonasc_expected(name);
        assert_eq!(got, want, "{name}");
    }
    for (label, args) in [
        (
            "n4_ambgt_blastn_t4",
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
            "n4_ambgt_blastn_mt1_t4",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-mt_mode",
                "1",
                "-num_threads",
                "4",
            ],
        ),
    ] {
        let got = run_losat_sub(&root, &d_ambgt, "blastn", &args);
        assert_refused(label, got);
    }

    // 5. n4_amb0_eq_0 (n4.nin@96 = 00000000 -> amb[0] = 0; NCBI reads file[-1] before its mapping
    // and crashes with SIGSEGV). Refused across all programs.
    let d_amb0 = root.join("n4_amb0_eq_0");
    cp_db(&d_amb0, "n4", true);
    poke(d_amb0.join("n4.nin"), 96, "00000000");
    for (sub, args) in [
        ("blastdbcmd", N4_DBCMD.to_vec()),
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
        ("blastdbcheck", N4_CHECK.to_vec()),
    ] {
        let got = run_losat_sub(&root, &d_amb0, sub, &args);
        assert_refused(&format!("n4_amb0_eq_0_{sub}"), got);
    }

    // 6. n4_amb_past_end (n4.nin@108 = 10000000 -> amb[3] = 2^28: OID 3's last byte lies 256 MiB
    // past the 114-byte .nsq; NCBI reads it and crashes with SIGSEGV). With `-negative_gilist
    // {304}` (OID 3 excluded) NCBI at 1 thread and `-mt_mode 1` never reads OID 3; at 2+ threads in
    // DB-split mode `x_FillSeqBuffer` reads the excluded OIDs after the last included one up to
    // the iteration range's end and crashes. A list that excludes every OID is scanned the same
    // way, and an alias file's OID range is read through the OID at its end. LOSAT refuses where
    // NCBI crashes.
    let d_past = root.join("n4_amb_past_end");
    cp_db(&d_past, "n4", true);
    poke(d_past.join("n4.nin"), 108, "10000000");
    // Alias files whose only filter is an OID range (OIDs 0..2 and 0..1): NCBI's iteration range
    // ends at OID 3 (2), and its loop reads that OID too before it compares it with the end.
    fs::write(
        d_past.join("n4r.nal"),
        b"TITLE r\nDBLIST n4\nFIRST_OID 1\nLAST_OID 3\n",
    )
    .unwrap();
    fs::write(
        d_past.join("n4r2.nal"),
        b"TITLE r2\nDBLIST n4\nFIRST_OID 1\nLAST_OID 2\n",
    )
    .unwrap();
    let mut past_refusal: Option<Vec<u8>> = None;
    for (name, args) in [
        (
            "n4_ambpast_blastn_neggi_t1",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gi304s,
                "-num_threads",
                "1",
            ],
        ),
        (
            "n4_ambpast_blastn_neggi_mt1_t4",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gi304s,
                "-mt_mode",
                "1",
                "-num_threads",
                "4",
            ],
        ),
        (
            "n4_ambpast_blastn_negall_t1",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gi_alls,
                "-num_threads",
                "1",
            ],
        ),
        (
            "n4_ambpast_blastn_range3_t1",
            vec![
                "-query",
                qns,
                "-db",
                "n4r",
                "-outfmt",
                "6",
                "-num_threads",
                "1",
            ],
        ),
        (
            "n4_ambpast_blastn_range2_t4",
            vec![
                "-query",
                qns,
                "-db",
                "n4r2",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
    ] {
        let got = run_losat_sub(&root, &d_past, "blastn", &args);
        let want = load_nonasc_expected(name);
        assert_eq!(got, want, "{name}");
    }
    for (name, args) in [
        (
            "n4_ambpast_blastn_t1",
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
            "n4_ambpast_blastn_t4",
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
            "n4_ambpast_blastn_neggi_t4",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gi304s,
                "-num_threads",
                "4",
            ],
        ),
        (
            "n4_ambpast_blastn_negall_t4",
            vec![
                "-query",
                qns,
                "-db",
                "n4",
                "-outfmt",
                "6",
                "-negative_gilist",
                n4_neg_gi_alls,
                "-num_threads",
                "4",
            ],
        ),
        (
            "n4_ambpast_blastn_range3_t4",
            vec![
                "-query",
                qns,
                "-db",
                "n4r",
                "-outfmt",
                "6",
                "-num_threads",
                "4",
            ],
        ),
    ] {
        // NCBI's frozen result (SIGSEGV, rc=-11, no output) stays beside the row for provenance.
        assert_eq!(load_nonasc_expected(name).2, -11, "{name}");
        let got = run_losat_sub(&root, &d_past, "blastn", &args);
        assert!(got.0.is_empty(), "{name}: stdout {:?}", text(&got.0));
        let err = got.1.clone();
        assert_refused(name, got);
        // Every row refuses as the first (1 thread, no list) does: the same read of OID 3's last
        // byte.
        match &past_refusal {
            None => {
                let past = "byte range 268435455..268435456 passes the end of the file";
                assert!(text(&err).contains(past), "{name}: {:?}", text(&err));
                past_refusal = Some(err);
            }
            Some(first) => assert_eq!(text(&err), text(first), "{name}"),
        }
    }

    let _ = fs::remove_dir_all(&root);
}
/// LOSAT's refusal where NCBI crashes (exit status 1, `Error: ...not supported by LOSAT`, no
/// standard output).
fn assert_tail_refusal(label: &str, got: &(Vec<u8>, Vec<u8>, i32)) {
    assert_eq!(text(&got.0), "", "{label}: stdout");
    assert_eq!(got.2, 1, "{label}: expected rc=1, err={:?}", text(&got.1));
    assert!(
        got.1.starts_with(b"Error: ") && text(&got.1).contains("not supported by LOSAT"),
        "{label}: expected LOSAT refusal, got stderr={:?}",
        text(&got.1)
    );
}

/// What a row of NCBI's frozen results asks of LOSAT: a crash (`SIGSEGV`, `rc=-11`) is refused,
/// anything else is byte-identical (standard output, standard error in LOSAT's words, status).
fn assert_tail_row(label: &str, got: &(Vec<u8>, Vec<u8>, i32), ncbi: (i32, &[u8], &[u8])) {
    if ncbi.0 == -11 {
        assert_tail_refusal(label, got);
        return;
    }
    let out = losat::wording::losat_output_bytes(ncbi.1).into_owned();
    let err = losat::wording::losat_bytes(ncbi.2).into_owned();
    assert_eq!(got, &(out, err, ncbi.0), "{label}");
}

/// The offsets of the arrays of a version 4 nucleotide `.nin`: (number of OIDs, header offsets,
/// sequence offsets, ambiguity offsets), each as a byte position in the file.
fn nin_v4_layout(nin: &[u8]) -> (usize, usize, usize, usize) {
    let be = |at: usize| u32::from_be_bytes(nin[at..at + 4].try_into().unwrap()) as usize;
    let mut at = 8;
    at += 4 + be(at);
    at += 4 + be(at);
    let count = be(at);
    at += 4 + 8 + 4;
    let seq = at + 4 * (count + 1);
    (count, at, seq, seq + 4 * (count + 1))
}

/// Sets the ambiguity offset of OID `oid` of `nin`: `2^28` makes NCBI read the last byte of the
/// sequence far past its file (`SIGSEGV`), the sequence's own start gives the signed length `-4`.
fn poke_amb_offset(nin: &Path, oid: usize, crash: bool) {
    let mut bytes = fs::read(nin).unwrap();
    let (count, _, seq, amb) = nin_v4_layout(&bytes);
    assert!(oid < count);
    let value: [u8; 4] = if crash {
        0x1000_0000u32.to_be_bytes()
    } else {
        bytes[seq + 4 * oid..seq + 4 * oid + 4].try_into().unwrap()
    };
    bytes[amb + 4 * oid..amb + 4 * oid + 4].copy_from_slice(&value);
    fs::write(nin, bytes).unwrap();
}

/// NCBI's frozen result `nonasc_xprog/<label>.{rc,out,err}` as NCBI wrote it: (exit status,
/// standard output, standard error). A label without a frozen row panics, so no row is skipped
/// for a typo.
fn frozen_tail_row(label: &str) -> (i32, Vec<u8>, Vec<u8>) {
    let base = fixtures().join("nonasc_xprog");
    let read = |ext: &str| {
        fs::read(base.join(format!("{label}.{ext}")))
            .unwrap_or_else(|error| panic!("{label}.{ext}: {error}"))
    };
    let rc = String::from_utf8(read("rc"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    (rc, read("out"), read("err"))
}

/// LOSAT's result of a row for the closing checks of a test: (exit status, lines of standard
/// output).
fn losat_status(got: &(Vec<u8>, Vec<u8>, i32)) -> (i32, usize) {
    (got.2, got.0.iter().filter(|&&byte| byte == b'\n').count())
}

/// The CPUs that this test may run on (`available_parallelism`: its affinity and cgroup quota).
/// LOSAT, like NCBI, caps `-num_threads` at the host's online CPUs (LOSAT's own count, which
/// falls back to this one), which are never fewer than these, so a row of NCBI's results at no
/// more threads than these runs as NCBI ran it (uncapped); a test skips a row at more threads.
fn usable_cpus() -> i64 {
    std::thread::available_parallelism().map_or(1, |n| n.get() as i64)
}

/// The guard of a test of rows at 1, 2 and 4 threads (`threads`: each row's), of which `ran`
/// ran: a row is skipped only for threads past the host's CPUs (`usable_cpus`), so on two CPUs
/// or more (every CI runner) every 1- and 2-thread row runs.
fn assert_rows_ran(test: &str, ran: usize, threads: &[i64]) {
    let cpus = usable_cpus();
    let runnable = threads.iter().filter(|&&t| t <= cpus).count();
    assert_eq!(
        ran,
        runnable,
        "{test}: {ran} of {} rows ran on {cpus} CPUs",
        threads.len()
    );
    assert!(ran > 0, "{test}: no row ran");
}

/// PR #1660 review R6 (`n6db/`: six sequences, OID 1 with the signed length `-4`, OID 5 with
/// its last byte `2^28` past the file): a block that the scan's own gap stops, before a last
/// block whose trailing excluded OID would crash. NCBI's results (`n6_<row>.{out,err,rc}`) are
/// `SIGSEGV` at one thread without a list, and rows without a crash otherwise (OID 0's rows only
/// at 2 and 4 threads); LOSAT is byte-identical, and refuses the crash. A row at more threads
/// than the host's CPUs is skipped (`usable_cpus`).
#[test]
fn nonasc_gap_stop_before_crash_parity() {
    let root = scratch("n6-gapstop");
    let _ = fs::remove_dir_all(&root);
    let db = root.join("db");
    fs::create_dir_all(&db).unwrap();
    let root = root.canonicalize().unwrap();
    let db = db.canonicalize().unwrap();
    let source = fixtures().join("nonasc_xprog/n6db");
    for entry in fs::read_dir(&source).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), db.join(entry.file_name())).unwrap();
    }
    let query = root.join("q.fna");
    fs::write(
        &query,
        b">q\nAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTG\n",
    )
    .unwrap();
    let lists = [
        ("neg135", "602\n604\n606\n"),
        ("neg35", "604\n606\n"),
        ("neg5", "606\n"),
    ];
    for (name, gis) in lists {
        fs::write(root.join(format!("{name}.txt")), gis).unwrap();
    }
    let query = query.to_str().unwrap().to_string();
    let mut rows = Vec::new();
    for (named, list) in [
        ("t", None),
        ("neg135_t", Some("neg135")),
        ("neg35_t", Some("neg35")),
        ("neg5_t", Some("neg5")),
    ] {
        for threads in [1i64, 2, 4] {
            rows.push((named, list, threads));
        }
    }
    let cpus = usable_cpus();
    let mut results = BTreeMap::new();
    for &(named, list, threads) in &rows {
        let label = format!("n6_{named}{threads}");
        let ncbi = frozen_tail_row(&label);
        if threads > cpus {
            continue;
        }
        let threads = threads.to_string();
        let mut args = vec!["-query", &query, "-db", "n6", "-outfmt", "6"];
        let list_path;
        if let Some(list) = list {
            list_path = root
                .join(format!("{list}.txt"))
                .to_str()
                .unwrap()
                .to_string();
            args.extend(["-negative_gilist", &list_path]);
        }
        args.extend(["-num_threads", &threads]);
        let got = run_losat_sub(&root, &db, "blastn", &args);
        assert_tail_row(&label, &got, (ncbi.0, &ncbi.1, &ncbi.2));
        results.insert(label, losat_status(&got));
    }
    let threads: Vec<i64> = rows.iter().map(|row| row.2).collect();
    assert_rows_ran(
        "nonasc_gap_stop_before_crash_parity",
        results.len(),
        &threads,
    );
    // LOSAT's results: one thread without a list is refused (NCBI's crash), the lists complete
    // at one thread, and at 2 threads the scan stops at OID 1 (OID 0's rows only).
    assert_eq!(results["n6_t1"], (1, 0));
    assert_eq!(results["n6_neg35_t1"], (0, 15));
    assert_eq!(results["n6_neg5_t1"], (0, 20));
    if cpus >= 2 {
        assert_eq!(results["n6_t2"], (0, 5));
        assert_eq!(results["n6_neg5_t2"], (0, 5));
    }
    let _ = fs::remove_dir_all(&root);
}

/// The sequence file of `pv1` in a row of `nonasc_protein_scan_stop_parity`.
#[derive(Clone, Copy, PartialEq)]
enum Pv1Psq {
    Kept,
    Empty,
    Missing,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:324-373,684-710
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1783-1838
/// The multithreaded scan of a protein database whose later volume's sequence file has no bytes
/// or is missing (RG-D2 tail R2 review finding 8, R4 review findings 1 and 2). `p2db/` has the
/// volumes `pv0` (GIs 7001-7003) and `pv1` (7004-7006) under an alias `p`; `p2_<row>` (blastp)
/// and `p2x_<row>` (blastx) hold NCBI's results at 1 and 2 threads. The scan reads every raw OID
/// of the iteration range, the excluded ones after the last included one too, up to the first
/// negative length; its first read in `pv1` opens `pv1.psq` (missing: NCBI's memory map error),
/// and a read of a positive length from the empty file crashes (`SIGSEGV`, refused by LOSAT).
/// `gapstop`: OID 1 of `pv0` has the length -1 and the list leaves it out before the block of
/// `pv1`, so the scan stops before `pv1`; `v1gap`: the negative OID is OID 1 of `pv1`;
/// `emptyv1_skipfirst`: no stop before the empty file. `trail`: the list holds only OID 0, so
/// only the excluded OIDs after it reach `pv1` (`trail_neg`: OID 1 of `pv0` of length -1 stops
/// them first). `none`: the list includes no OID; one thread warns that the filtering left no
/// database, two threads read the whole range first. `v1firstneg`: OID 0 of `pv1` has the
/// length -1 and `pv1.psq` is empty: a non-positive length reads nothing, so two threads stop
/// there and one thread crashes at the next OID; with the list of OIDs 0 and 3 nothing is read
/// from the file (blastx: NCBI's error of the length -1). A row at more threads than the host's
/// CPUs is skipped (`usable_cpus`).
#[test]
fn nonasc_protein_scan_stop_parity() {
    use Pv1Psq::{Empty, Kept, Missing};
    let root = scratch("p2-stop");
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let root = root.canonicalize().unwrap();
    let base = fixtures().join("nonasc_xprog");
    let protein_query = root.join("q.faa");
    fs::write(
        &protein_query,
        b">q\nYKNTARICGENSIPVEWIAHQKGPGDFYYRFFAAHHGGLMHVHGHPLANQGFKDMLYWAY\n",
    )
    .unwrap();
    // The codons of the protein query.
    let nucleotide_query = root.join("q.fna");
    fs::write(
        &nucleotide_query,
        concat!(
            ">qn\nTATAAAAATACTGCTCGTATTTGTGGTGAAAATTCTATTCCTGTTGAATGGATTGCTCATCAAAAAGGTCCTGG",
            "TGATTTTTATTATCGTTTTTTTGCTGCTCATCATGGTGGTCTTATGCATGTTCATGGTCATCCTCTTGCTAATCAAGGT",
            "TTTAAAGATATGCTTTATTGGGCTTAT\n"
        ),
    )
    .unwrap();
    // Protein `oid` of the volume index `pin` gets the length -1: its end (the next start) is
    // set to its start.
    let negative = |pin: &Path, oid: usize| {
        let mut bytes = fs::read(pin).unwrap();
        let (count, _, seq, _) = nin_v4_layout(&bytes);
        assert!(oid < count);
        let start: [u8; 4] = bytes[seq + 4 * oid..seq + 4 * oid + 4].try_into().unwrap();
        bytes[seq + 4 * (oid + 1)..seq + 4 * (oid + 1) + 4].copy_from_slice(&start);
        fs::write(pin, bytes).unwrap();
    };
    // (row, GI list, the OID of length -1, `pv1.psq`, blastx too)
    let cases: [(&str, Option<&str>, Option<(&str, usize)>, Pv1Psq, bool); 13] = [
        (
            "gapstop",
            Some("7001\n7004\n7005\n7006\n"),
            Some(("pv0", 1)),
            Kept,
            false,
        ),
        (
            "gapstop_emptyv1",
            Some("7001\n7004\n7005\n7006\n"),
            Some(("pv0", 1)),
            Empty,
            false,
        ),
        (
            "gapstop_missingv1",
            Some("7001\n7004\n7005\n7006\n"),
            Some(("pv0", 1)),
            Missing,
            false,
        ),
        (
            "v1gap_emptyv1",
            Some("7001\n7006\n"),
            Some(("pv1", 1)),
            Empty,
            false,
        ),
        (
            "v1gap_missingv1",
            Some("7001\n7006\n"),
            Some(("pv1", 1)),
            Missing,
            false,
        ),
        (
            "emptyv1_skipfirst",
            Some("7001\n7006\n"),
            None,
            Empty,
            false,
        ),
        ("trail_missingv1", Some("7001\n"), None, Missing, true),
        (
            "trail_neg_missingv1",
            Some("7001\n"),
            Some(("pv0", 1)),
            Missing,
            true,
        ),
        ("trail_emptyv1", Some("7001\n"), None, Empty, true),
        ("none_missingv1", Some("9999\n"), None, Missing, true),
        (
            "none_neg_missingv1",
            Some("9999\n"),
            Some(("pv0", 1)),
            Missing,
            true,
        ),
        ("v1firstneg_emptyv1", None, Some(("pv1", 0)), Empty, true),
        (
            "v1firstneg_emptyv1_list",
            Some("7001\n7004\n"),
            Some(("pv1", 0)),
            Empty,
            true,
        ),
    ];
    let mut rows = Vec::new();
    for case in cases {
        for program in ["blastp", "blastx"] {
            if program == "blastp" || case.4 {
                for threads in [1i64, 2] {
                    rows.push((case, program, threads));
                }
            }
        }
    }
    let cpus = usable_cpus();
    let mut results = BTreeMap::new();
    for &((named, list, negative_oid, psq, _), program, threads) in &rows {
        let prefix = if program == "blastp" { "p2" } else { "p2x" };
        let label = format!("{prefix}_{named}_t{threads}");
        let ncbi = frozen_tail_row(&label);
        if threads > cpus {
            continue;
        }
        let dir = root.join(format!("run_{label}"));
        fs::create_dir_all(&dir).unwrap();
        for entry in fs::read_dir(base.join("p2db")).unwrap() {
            let entry = entry.unwrap();
            fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
        }
        fs::write(dir.join("p.pal"), "TITLE p\nDBLIST pv0 pv1\n").unwrap();
        if let Some((volume, oid)) = negative_oid {
            negative(&dir.join(format!("{volume}.pin")), oid);
        }
        match psq {
            Kept => {}
            Empty => fs::write(dir.join("pv1.psq"), b"").unwrap(),
            Missing => fs::remove_file(dir.join("pv1.psq")).unwrap(),
        }
        let query = if program == "blastp" {
            &protein_query
        } else {
            &nucleotide_query
        };
        let threads = threads.to_string();
        let mut args = vec![
            "-query",
            query.to_str().unwrap(),
            "-db",
            "p",
            "-outfmt",
            "6",
        ];
        if let Some(list) = list {
            fs::write(dir.join("gl.txt"), list).unwrap();
            args.extend(["-gilist", "gl.txt"]);
        }
        args.extend(["-num_threads", &threads]);
        let got = run_losat_sub(&root, &dir, program, &args);
        assert_tail_row(&label, &got, (ncbi.0, &ncbi.1, &ncbi.2));
        results.insert(label, losat_status(&got));
    }
    let threads: Vec<i64> = rows.iter().map(|row| row.2).collect();
    assert_rows_ran("nonasc_protein_scan_stop_parity", results.len(), &threads);
    // LOSAT's results. One thread reads only the included OIDs: the missing file after them is
    // never opened, the list of no OID warns, and the empty file is refused where a positive
    // length is read from it (NCBI's crash).
    assert_eq!(results["p2_trail_missingv1_t1"], (0, 1));
    assert_eq!(results["p2x_trail_missingv1_t1"], (0, 1));
    assert_eq!(results["p2_none_missingv1_t1"], (0, 0));
    assert_eq!(results["p2_v1firstneg_emptyv1_t1"], (1, 0));
    assert_eq!(results["p2_v1firstneg_emptyv1_list_t1"], (0, 1));
    if cpus >= 2 {
        // Two threads: the stop before `pv1` completes; the excluded OIDs after the last
        // included one open `pv1.psq` (missing: the error, also where no OID is included)
        // unless a negative length stops them first; the length -1 of OID 0 of `pv1` stops the
        // scan before the empty file.
        assert_eq!(results["p2_gapstop_missingv1_t2"], (0, 1));
        assert_eq!(results["p2_emptyv1_skipfirst_t2"], (1, 0));
        assert_eq!(results["p2_trail_missingv1_t2"], (3, 0));
        assert_eq!(results["p2x_trail_missingv1_t2"], (3, 0));
        assert_eq!(results["p2_trail_neg_missingv1_t2"], (0, 1));
        assert_eq!(results["p2_none_missingv1_t2"], (3, 0));
        assert_eq!(results["p2x_none_missingv1_t2"], (3, 0));
        assert_eq!(results["p2_none_neg_missingv1_t2"], (0, 0));
        assert_eq!(results["p2_v1firstneg_emptyv1_t2"], (0, 3));
        assert_eq!(results["p2x_v1firstneg_emptyv1_t2"], (0, 3));
    }
    let _ = fs::remove_dir_all(&root);
}

/// A later volume whose sequence file has no bytes (RG-D2 tail, item 3). `mv3db/` has three
/// volumes `v0` (5 OIDs), `v1` (6) and `v2` (4); `mv3e_rows.tsv` holds NCBI's results for
/// four kinds of the database (`empty`: `v1.nsq` of 0 bytes; `neg`: and OID 2 of `v0` with the
/// signed length `-4`; `zero`: and the first sequence offset of `v1` 0; `negzero`: both), the
/// aliases `n` (the three volumes), `r_head`, `r_tail` and `r_pre` (OID ranges), with or without
/// a GI list, at 1, 2 and 4 threads. NCBI's read at an empty file's sequence offset 0 (a null
/// pointer) gives the length -1, which stops the multithreaded scan and which one thread skips
/// (`posv1first`: the list of OID 5 only, the first of `v1`), and it crashes where it reads
/// another offset; LOSAT matches both. A row at more threads than the host's CPUs is skipped
/// (`usable_cpus`).
#[test]
fn nonasc_empty_later_volume_sequence_file_parity() {
    let root = scratch("mv3-empty");
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let root = root.canonicalize().unwrap();
    let base = fixtures().join("nonasc_xprog");
    let query = root.join("qn.fna");
    fs::write(
        &query,
        b">qn\nAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTGAATCCAGGCTTG\n",
    )
    .unwrap();
    let query = query.to_str().unwrap().to_string();
    let gis = |oids: &[u32]| {
        oids.iter()
            .map(|o| format!("{}\n", 601 + o))
            .collect::<String>()
    };
    let all_v1: Vec<u32> = (5..11).collect();
    let all_v0: Vec<u32> = (0..5).collect();
    let all_v2: Vec<u32> = (11..15).collect();
    let v0_v2: Vec<u32> = all_v0.iter().chain(&all_v2).copied().collect();
    let list_of = |name: &str| -> Option<(&'static str, String)> {
        match name {
            "none" => None,
            "negv1" => Some(("-negative_gilist", gis(&all_v1))),
            "negv1first" => Some(("-negative_gilist", gis(&[5]))),
            "negv0" => Some(("-negative_gilist", gis(&[0]))),
            "posv0v2" => Some(("-gilist", gis(&v0_v2))),
            "posv0" => Some(("-gilist", gis(&all_v0))),
            "negoid2" => Some(("-negative_gilist", gis(&[2]))),
            "posv1first" => Some(("-gilist", gis(&[5]))),
            "posv1first2" => Some(("-gilist", gis(&[5, 6]))),
            other => panic!("unknown list {other}"),
        }
    };
    let setup = |kind: &str| -> PathBuf {
        let dir = root.join(format!("run_{kind}"));
        fs::create_dir_all(&dir).unwrap();
        for entry in fs::read_dir(base.join("mv3db")).unwrap() {
            let entry = entry.unwrap();
            fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
        }
        fs::write(dir.join("n.nal"), "TITLE n\nDBLIST v0 v1 v2\n").unwrap();
        for (name, first, last) in [("r_head", 1, 7), ("r_tail", 12, 15), ("r_pre", 1, 5)] {
            let head = format!("TITLE r\nDBLIST n\nFIRST_OID {first}\nLAST_OID {last}\n");
            fs::write(dir.join(format!("{name}.nal")), head).unwrap();
        }
        fs::write(dir.join("v1.nsq"), b"").unwrap();
        if kind == "neg" || kind == "negzero" {
            poke_amb_offset(&dir.join("v0.nin"), 2, false);
        }
        if kind == "zero" || kind == "negzero" {
            let mut nin = fs::read(dir.join("v1.nin")).unwrap();
            let (_, _, seq, _) = nin_v4_layout(&nin);
            nin[seq..seq + 4].copy_from_slice(&0u32.to_be_bytes());
            fs::write(dir.join("v1.nin"), nin).unwrap();
        }
        dir
    };
    let dirs: Vec<(String, PathBuf)> = ["empty", "neg", "zero", "negzero"]
        .iter()
        .map(|kind| (kind.to_string(), setup(kind)))
        .collect();
    let table = fs::read_to_string(base.join("mv3e_rows.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = table
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    let cpus = usable_cpus();
    let mut results = BTreeMap::new();
    let mut threads_of_rows = Vec::new();
    for cell in &rows {
        let label = cell[0];
        let (rest, threads) = label
            .strip_prefix("mv3e_")
            .unwrap()
            .rsplit_once("_t")
            .unwrap();
        let threads: i64 = threads.parse().unwrap();
        threads_of_rows.push(threads);
        if threads > cpus {
            continue;
        }
        let ncbi_rc: i32 = cell[1].parse().unwrap();
        let (kind, rest) = rest.split_once('_').unwrap();
        let threads = threads.to_string();
        let (db, list) = ["r_head_", "r_tail_", "r_pre_", "n_"]
            .iter()
            .find_map(|prefix| {
                rest.strip_prefix(prefix)
                    .map(|l| (prefix.trim_end_matches('_'), l))
            })
            .unwrap();
        let dir = &dirs.iter().find(|(k, _)| k == kind).unwrap().1;
        let mut args = vec![
            "-query",
            &query,
            "-db",
            db,
            "-outfmt",
            "6",
            "-num_threads",
            &threads,
        ];
        let list_path = dir.join(format!("{list}.txt"));
        if let Some((flag, content)) = list_of(list) {
            fs::write(&list_path, content).unwrap();
            args.extend([flag, list_path.to_str().unwrap()]);
        }
        let got = run_losat_sub(&root, dir, "blastn", &args);
        let file = |name: &str| {
            if name == "-" {
                Vec::new()
            } else {
                fs::read(base.join(name)).unwrap()
            }
        };
        assert_tail_row(label, &got, (ncbi_rc, &file(cell[2]), &file(cell[3])));
        results.insert(label.to_string(), losat_status(&got));
    }
    assert_rows_ran(
        "nonasc_empty_later_volume_sequence_file_parity",
        results.len(),
        &threads_of_rows,
    );
    // LOSAT's results: every row of the table at 4 threads or fewer, 59 refused (NCBI's crashes)
    // and 145 completed, among them the skip at one thread and the stop at two of a read at
    // offset 0 of the empty file.
    let refused = results.values().filter(|status| status.0 == 1).count();
    let completed = results.values().filter(|status| status.0 == 0).count();
    assert_eq!(refused + completed, results.len(), "{results:?}");
    if cpus >= 4 {
        assert_eq!((refused, completed), (59, 145));
    }
    assert_eq!(results["mv3e_zero_n_posv1first_t1"], (0, 0));
    assert_eq!(results["mv3e_zero_n_posv1first2_t1"], (1, 0));
    if cpus >= 2 {
        assert_eq!(results["mv3e_zero_n_posv1first2_t2"], (0, 0));
    }
    let _ = fs::remove_dir_all(&root);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:374-399
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
/// The search of the repeats of `-filtering_db` is NCBI's `CLocalBlast` of one thread, whatever
/// `-num_threads` gives the search (RG-D2 tail R4 item 5). `fdb/` holds `main` (`s0` holds R, the
/// first 150 bases of the query `fdb_q.fna`, `s1` holds U, its last 150) and the databases of
/// repeats `rep` (`r2` holds R) and `repa` (`r1` holds R), made with NCBI `makeblastdb` (v4);
/// `fdb_<case>_t<T>.{out,err,rc}` are NCBI's results at 1 and 2 threads. `clean`: R is a repeat,
/// so only `s1` is reported. `neg`: OID 1 of `rep` has the signed length -3; the one-thread scan
/// reads on past it and finds R (a multithreaded scan would stop there and leave R unmasked).
/// `repa`'s OID 1 has its ambiguity offset 4 bytes past the next sequence's start (`ambnext`):
/// the repeat search's traceback cannot read its ambiguity data, and each query's set-up posts
/// the error (exit status 3), as at one thread (a multithreaded traceback aborts). A row at more
/// threads than the host's CPUs is skipped (`usable_cpus`).
#[test]
fn nonasc_filtering_db_repeat_search_runs_one_thread() {
    let root = scratch("fdb-one-thread");
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let root = root.canonicalize().unwrap();
    let base = fixtures().join("nonasc_xprog");
    let query = base.join("fdb_q.fna");
    let query = query.to_str().unwrap().to_string();
    let cpus = usable_cpus();
    let mut results = BTreeMap::new();
    let mut threads_of_rows = Vec::new();
    for (named, repeats) in [("clean", "rep"), ("neg", "rep"), ("ambnext", "repa")] {
        for threads in [1i64, 2] {
            let label = format!("fdb_{named}_t{threads}");
            let ncbi = frozen_tail_row(&label);
            threads_of_rows.push(threads);
            if threads > cpus {
                continue;
            }
            let dir = root.join(format!("run_{label}"));
            fs::create_dir_all(&dir).unwrap();
            for entry in fs::read_dir(base.join("fdb")).unwrap() {
                let entry = entry.unwrap();
                fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
            }
            if named == "neg" {
                poke_amb_offset(&dir.join("rep.nin"), 1, false);
            } else if named == "ambnext" {
                let nin = dir.join("repa.nin");
                let mut bytes = fs::read(&nin).unwrap();
                let (_, _, seq, amb) = nin_v4_layout(&bytes);
                let next = u32::from_be_bytes(bytes[seq + 8..seq + 12].try_into().unwrap());
                bytes[amb + 4..amb + 8].copy_from_slice(&(next + 4).to_be_bytes());
                fs::write(&nin, bytes).unwrap();
            }
            let threads = threads.to_string();
            let args = [
                "-query",
                &query,
                "-db",
                "main",
                "-filtering_db",
                repeats,
                "-outfmt",
                "6",
                "-num_threads",
                &threads,
            ];
            let got = run_losat_sub(&root, &dir, "blastn", &args);
            assert_tail_row(&label, &got, (ncbi.0, &ncbi.1, &ncbi.2));
            results.insert(label, got);
        }
    }
    assert_rows_ran(
        "nonasc_filtering_db_repeat_search_runs_one_thread",
        results.len(),
        &threads_of_rows,
    );
    // LOSAT's results: R is masked past the negative length (the rows of `clean`), and the
    // traceback's error is posted, at 2 threads as at 1.
    let mut threads = vec![1];
    if cpus >= 2 {
        threads.push(2);
    }
    for t in threads {
        let got = |named: &str| &results[&format!("fdb_{named}_t{t}")];
        assert_eq!(losat_status(got("clean")), (0, 1));
        assert_eq!(got("neg"), got("clean"));
        assert_eq!(losat_status(got("ambnext")), (3, 0));
    }
    let _ = fs::remove_dir_all(&root);
}
