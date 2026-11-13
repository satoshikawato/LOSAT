//! `losat blastdbcmd` on alias files with the keys of line RG, series D, item D1a
//! (`remaining-gaps-20261112/D/DESIGN_D1a.md`): `NSEQ` and `LENGTH` read as NCBI's walks of the alias
//! tree (values that are not the sums of the volumes included), `STATS_NSEQ`, `STATS_TOTLEN` and
//! `MIN_SEQ_LENGTH` (read, never printed), `FIRST_OID` / `LAST_OID` ranges (nested ranges ANDed, two
//! alias files over one volume ORed, the global OIDs of several volumes), `VDBLIST` (not read), a
//! volume that an alias file and a name both give, and NCBI's `CStringException` for a value that
//! does not convert, and (item D1a-2) `-ipg` / `-ipg_batch` for a PIG whose OID an alias file's
//! range leaves out (NCBI's `PigToOid` asks the volumes, and its formatters read that OID). The
//! expected bytes are NCBI BLAST+ 2.17.0 blastdbcmd's (the release build on the developer machine,
//! 2026-11-12 and 2026-11-13, `remaining-gaps-20261112/D/fixtures/make_d1a_blastdbcmd.py` and
//! `make_d1a2_blastdbcmd.py`); the alias files are in `tests/blastdbcmd/db/ak/`, the batch files
//! in `tests/blastdbcmd/in/`. `{R}` stands for the directory of the runs (`tests/blastdbcmd`). NCBI
//! is not run here.

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
    // D1a-2: `-ipg` / `-ipg_batch` read the OID of a PIG whether the alias files' ranges include
    // it or not (`PigToOid` asks the volumes; probe p4c_c 0-25 of D1a-2). A52 is approved difference
    // 5 of PD-LOSAT-BLASTDBCMD-DEVIATIONS (a mask on a protein in FASTA without a range): its
    // expected bytes are NCBI's header with the letters of NCBI's `-outfmt %s` run of the request.
    ("A33", &["-db", "db/ak/m_r", "-ipg", "1001"], 0),
    ("A34", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%o|%a|%g|%t|%l|%P|%T|%L|%e|%n|%h|%i"], 0),
    ("A35", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%s"], 0),
    ("A36", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%m21"], 0),
    ("A39", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%f"], 0),
    ("A40", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%X|%K|%S|%B|%C|%N"], 0),
    ("A41", &["-db", "db/ak/m_r", "-ipg", "1002", "-outfmt", "%o|%a|%l"], 0),
    ("A42", &["-db", "db/m/pmg", "-ipg", "1001", "-outfmt", "%o|%a|%g|%t|%l|%P|%T|%L|%e|%n|%h|%i"], 0),
    ("A43", &["-db", "db/m/pmg", "-ipg", "1001"], 0),
    ("A44", &["-db", "db/ak/m_s", "-ipg", "1003", "-outfmt", "%o|%a|%l|%s"], 0),
    ("A45", &["-db", "db/ak/m_s", "-ipg", "1003"], 0),
    ("A46", &["-db", "db/ak/m_t", "-ipg", "1001", "-outfmt", "%o|%a|%l"], 0),
    ("A47", &["-db", "db/ak/m_r", "-ipg", "9999"], 1),
    ("A48", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%o|%l"], 1),
    ("A49", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%o"], 0),
    ("A50", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%s"], 1),
    ("A51", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%a"], 1),
    ("A52", &["-db", "db/ak/m_r", "-ipg", "1001", "-mask_sequence_with", "21"], 0),
    ("A53", &["-db", "db/ak/m_r", "-ipg", "1001", "-range", "2-10"], 0),
    ("A54", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%s", "-range", "2-10", "-mask_sequence_with", "21"], 0),
    ("A55", &["-db", "db/ak/m_r", "-ipg_batch", "in/ak_pig.txt"], 1),
    ("A56", &["-db", "db/ak/m_r", "-ipg_batch", "in/ak_pig.txt", "-outfmt", "%o %a %l"], 1),
    ("A57", &["-db", "db/ak/m_r", "-ipg_batch", "in/ak_pigm.txt", "-outfmt", "%o %s"], 0),
    ("A58", &["-db", "db/ak/m_s", "-ipg_batch", "in/ak_pig.txt", "-outfmt", "%o %a"], 1),
];

/// The cases of probe p4c_c that LOSAT's BLASTDBCMD refuses: `-outfmt %d` and `%b` (the ASN.1
/// texts, line BDC DECISIONS 9), which NCBI writes for the OID of a PIG outside the range too.
#[rustfmt::skip]
const REFUSALS: &[(&str, &[&str])] = &[
    ("A37", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%d"]),
    ("A38", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%b"]),
];

#[test]
fn asn1_texts_of_a_pig_outside_the_range_are_refused() {
    let config = TestDir::new("refusals");
    for &(id, words) in REFUSALS {
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert_eq!(stdout, b"", "{id} {words:?}");
        assert_eq!(code, 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n"),
            "{id} {words:?}: {stderr}"
        );
    }
}

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
