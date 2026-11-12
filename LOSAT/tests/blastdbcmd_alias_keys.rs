//! `losat blastdbcmd` on alias files with the keys of line RG, series D, item D1a
//! (`remaining-gaps-20261112/D/DESIGN_D1a.md`): `NSEQ` and `LENGTH` read as NCBI's walks of the alias
//! tree (values that are not the sums of the volumes included), `STATS_NSEQ`, `STATS_TOTLEN` and
//! `MIN_SEQ_LENGTH` (read, never printed), `FIRST_OID` / `LAST_OID` ranges (nested ranges ANDed, two
//! alias files over one volume ORed, the global OIDs of several volumes), `VDBLIST` (not read), a
//! volume that an alias file and a name both give, and NCBI's `CStringException` for a value that
//! does not convert. The expected bytes are NCBI BLAST+ 2.17.0 blastdbcmd's (the release build on
//! the developer machine, 2026-11-12, `remaining-gaps-20261112/D/fixtures/make_d1a_blastdbcmd.py`);
//! the alias files are in `tests/blastdbcmd/db/ak/`. `{R}` stands for the directory of the runs
//! (`tests/blastdbcmd`). NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// `fs::canonicalize` without the `\\\\?\\` prefix that Windows puts on its result (as
/// `blastdbcmd_retrieve.rs`).
fn canonical(path: &Path) -> PathBuf {
    let path = fs::canonicalize(path).unwrap();
    #[cfg(windows)]
    if let Some(plain) = path.to_str().and_then(|text| text.strip_prefix(r"\\?\")) {
        if plain.as_bytes().get(1) == Some(&b':') {
            return PathBuf::from(plain);
        }
    }
    path
}

fn run_dir() -> PathBuf {
    canonical(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("blastdbcmd"),
    )
}

/// A directory of its own for one test (the tests of a binary run in parallel), removed at
/// the end.
struct TestDir(PathBuf);

impl TestDir {
    fn new(stem: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "losat_blastdbcmd_alias_keys_{stem}_{}_{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        TestDir(canonical(&dir))
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One run of `losat blastdbcmd words...` in `tests/blastdbcmd`: stdout, stderr, exit status.
/// No `BLASTDB`, no registry file (`NCBI_CONFIG_PATH` names `config` alone), no LOSAT variable,
/// nothing on the standard input.
fn blastdbcmd(words: &[&str], config: &Path) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    let output = command
        .env("NCBI_CONFIG_PATH", config)
        .stdin(std::process::Stdio::null())
        .current_dir(run_dir())
        .arg("blastdbcmd")
        .args(words)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// The frozen file `name` of `tests/blastdbcmd/expected`, as it is.
fn frozen(name: &str) -> Vec<u8> {
    fs::read(run_dir().join("expected").join(name)).unwrap()
}

/// `text` with `{R}` placed.
fn placed(text: &[u8]) -> Vec<u8> {
    let root = run_dir();
    let root = root.to_str().unwrap().as_bytes();
    let mut placed = Vec::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.windows(3).position(|window| window == b"{R}") {
        placed.extend_from_slice(&rest[..at]);
        placed.extend_from_slice(root);
        rest = &rest[at + 3..];
    }
    placed.extend_from_slice(rest);
    placed
}

fn holds_root(text: &[u8]) -> bool {
    text.windows(3).any(|window| window == b"{R}")
}

/// The cases: id, words, exit status.
#[rustfmt::skip]
const CASES: &[(&str, &[&str], i32)] = &[
    ("A01", &["-db", "db/ak/p_nseq", "-info"], 0),
    ("A02", &["-db", "db/ak/p_nseq", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A03", &["-db", "db/ak/j_nseq", "-metadata"], 0),
    ("A04", &["-db", "db/ak/p_range", "-info"], 0),
    ("A05", &["-db", "db/ak/p_range", "-info", "-exact_length"], 0),
    ("A06", &["-db", "db/ak/p_range", "-entry", "all", "-outfmt", "%o %a %l %T"], 0),
    ("A07", &["-db", "db/ak/p_range", "-entry", "XP_000001.1", "-outfmt", "%o %a"], 1),
    ("A08", &["-db", "db/ak/p_range", "-entry", "Q00003,A0A000", "-outfmt", "%o %a"], 0),
    ("A09", &["-db", "db/ak/p_range", "-entry", "all"], 0),
    ("A10", &["-db", "db/ak/p_union", "-info"], 0),
    ("A11", &["-db", "db/ak/p_union", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A12", &["-db", "db/ak/p_twice", "-info"], 0),
    ("A13", &["-db", "db/ak/p_twice", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A14", &["-db", "db/ak/p_stats", "-info"], 0),
    ("A15", &["-db", "db/ak/p_nseq0", "-info"], 0),
    ("A16", &["-db", "db/ak/d_glob", "-info"], 0),
    ("A17", &["-db", "db/ak/d_glob", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A18", &["-db", "db/ak/d_rev", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A19", &["-db", "db/ak/n_vdb", "-info"], 0),
    ("A20", &["-db", "db/ak/n_vdbonly", "-info", "-dbtype", "nucl"], 2),
    ("A21", &["-db", "db/ak/p_bad_nseq", "-info"], 255),
    ("A22", &["-db", "db/ak/p_bad_nseq", "-info", "-dbtype", "prot"], 255),
    ("A23", &["-db", "db/ak/p_bad_first", "-entry", "all", "-dbtype", "prot"], 255),
    ("A24", &["-db", "db/ak/p_nest", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A25", &["-db", "db/ak/p_empty", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A26", &["-db", "db/ak/p_empty", "-info"], 0),
    ("A27", &["-db", "db/ak/p_nseq4g", "-info"], 0),
    ("A28", &["-db", "db/ak/p_min", "-info"], 0),
    ("A29", &["-db", "db/ak/j_range", "-metadata"], 0),
    ("A30", &["-db", "db/palias db/p5", "-info"], 0),
    ("A31", &["-db", "db/palias db/p5", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A32", &["-db", "db/palias db/p5", "-entry", "P99999.2", "-outfmt", "%o %a"], 0),
];

#[test]
fn alias_keys_are_ncbi_s_bytes() {
    let config = TestDir::new("cases");
    for &(id, words, status) in CASES {
        let (want_out, want_err) = (
            frozen(&format!("{id}.stdout")),
            frozen(&format!("{id}.stderr")),
        );
        // A path in the lines is written with the platform's separators, and `-metadata` cuts
        // a name at the platform's separator alone.
        if cfg!(windows)
            && (words.contains(&"-metadata") || holds_root(&want_out) || holds_root(&want_err))
        {
            continue;
        }
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        assert_eq!(
            String::from_utf8_lossy(&stdout),
            String::from_utf8_lossy(&placed(&want_out)),
            "{id} {words:?}: stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            String::from_utf8_lossy(&placed(&want_err)),
            "{id} {words:?}: stderr"
        );
        assert_eq!(code, status, "{id} {words:?}: exit status");
    }
}
