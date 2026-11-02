//! `losat blastdbcmd -db "a b"`: a `-db` of several names, read as one database, against NCBI
//! BLAST+ 2.17.0 blastdbcmd. The expected bytes are NCBI's, frozen by line BDC's probe N on ws
//! (`tests/blastdbcmd/README.md`), but for the differences 2, 3, 5 and 6 of
//! `PD-LOSAT-BLASTDBCMD-DEVIATIONS`. `{R}` stands for the directory of the probe's runs (here
//! `tests/blastdbcmd`, the current directory of each run). NCBI is not run here.

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
            "losat_blastdbcmd_names_{stem}_{}_{nanos}",
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

/// The probe's cases: id, words, exit status. The shapes of the list (quotes, spaces, one name
/// three times), a name not found, another molecule, versions 4 and 5 together, each output on
/// a list. `db/palias` is an alias of `p5` and `pm5`; `j/lst` holds the taxonomy files that `db/p5`
/// has not; `mask_nucl_ids` and `mask_nucl_y` give one algorithm the ids 101 and 100. Approved differences (expected bytes from NCBI's own neighbouring
/// run, line BDC's `probe/make_expected.py`): 2, `gnl|BL_ORD_ID|n` (`N846`); 3, a mask and
/// `-strand minus` in FASTA (`N1214`); 5, a mask on a protein in FASTA (`N1359`); 6, the mask of
/// the record before, on one accession of two sequences (`N1742`; this widening of item 6 was
/// approved with version 1.4 of the record).
#[rustfmt::skip]
const CASES: &[(&str, &[&str], i32)] = &[
    ("I37", &["-db", "db/p5 db/n5", "-info"], 2),
    ("I38", &["-db", "db/p5 db/pm5", "-info"], 0),
    ("F34", &["-db", "db/p5 db/pm5", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("N119", &["-db", "\"db/p5\" \"db/pm5\"", "-info"], 0),
    ("N2013", &["-db", "\"j/lst/p5\" \"j/lst/p5n\"", "-metadata"], 0),
    ("N128", &["-db", "\"db/p5 db/pm5\"", "-info"], 2),
    ("N146", &["-db", "\"\"", "-info"], 2),
    ("N194", &["-db", "db/p5 db/p5 db/p5", "-info"], 0),
    ("N2001", &["-db", "j/lst/p5 j/lst/p5n", "-metadata"], 0),
    ("N227", &["-db", "db/p5 db/none", "-info"], 2),
    ("N235", &["-db", "db/none db/p5", "-info", "-dbtype", "prot"], 2),
    ("N257", &["-db", "db/p4 db/p5", "-info"], 2),
    ("N259", &["-db", "db/p4 db/p5", "-info", "-dbtype", "prot"], 2),
    ("N249", &["-db", "db/p5 db/n5", "-info", "-dbtype", "nucl"], 2),
    ("N729", &["-db", "db/p5 db/pm5", "-entry", "P10001.1", "-outfmt", "%o %a %l"], 0),
    ("N735", &["-db", "db/p5 db/pm5", "-entry", "XP_000001.1,P10001.1", "-outfmt", "%o %a %l"], 0),
    ("N769", &["-db", "db/p5 db/p5n", "-entry", "XP_000001.1", "-outfmt", "%o %a %l"], 1),
    ("N846", &["-db", "db/p4 db/p4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%o %a %l"], 0),
    ("N743", &["-db", "db/p5 db/pm5", "-entry_batch", "in/n_p.txt", "-outfmt", "%o %a"], 1),
    ("N946", &["-db", "db/d5 db/d5b", "-entry", "XP_800001", "-get_dups", "-outfmt", "%o %a"], 0),
    ("N2003", &["-db", "j/lst/p5 j/lst/p5n", "-tax_info"], 0),
    ("N443", &["-db", "db/p4 db/p4n", "-tax_info"], 2),
    ("N2005", &["-db", "j/lst/p5 j/lst/p5n", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a %T"], 0),
    ("N2006", &["-db", "j/lst/p5 j/lst/p5n", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a %T"], 0),
    ("N1178", &["-db", "db/m/pmg db/p4", "-ipg", "1001", "-outfmt", "%o %a %P"], 0),
    ("N557", &["-db", "db/nm5 db/n5", "-info"], 0),
    ("N1211", &["-db", "db/nm5 db/n5", "-entry", "all", "-outfmt", "%o|%a|%l|%m11"], 0),
    ("N1212", &["-db", "db/nm5 db/n5", "-entry", "all", "-mask_sequence_with", "11"], 0),
    ("N1214", &["-db", "db/nm5 db/n5", "-entry", "all", "-mask_sequence_with", "11", "-strand", "minus"], 0),
    ("N1359", &["-db", "db/pm5 db/p5", "-entry", "all", "-mask_sequence_with", "21"], 0),
    ("N1742", &["-db", "db/nm4 db/m/nms4", "-entry", "NM_100001.1", "-mask_sequence_with", "12"], 0),
    ("N1451", &["-db", "../seqdb/db/mask_nucl_ids ../seqdb/db/mask_nucl_y", "-entry", "all", "-outfmt", "%o|%a|%l|%m101"], 0),
    ("N1446", &["-db", "../seqdb/db/mask_nucl_ids ../seqdb/db/mask_nucl_y", "-entry", "all", "-mask_sequence_with", "wmY"], 0),
    ("N1581", &["-db", "db/p5 db/pm5", "-list", "db/j/lst"], 1),
    ("N373", &["-db", "db/palias db/pv", "-info"], 0),
];

#[test]
fn several_names_are_ncbi_s_bytes() {
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

/// What LOSAT refuses on a list, with the lines written before the refusal: a volume that an alias
/// and a name both give, for every option, as the search programs do (NCBI's `-info` counts its
/// sequences twice, 39 for 21; its `-entry` answers); `-ipg`
/// without PIG files, where NCBI writes an OID it never set; an algorithm that the volume of a
/// sequence does not know, where NCBI reads the masks of another one.
#[rustfmt::skip]
const REFUSED: &[(&str, &[&str], usize)] = &[
    ("N389", &["-db", "db/palias db/p5", "-info"], 0),
    ("N1198", &["-db", "db/p4 db/p4n", "-ipg", "1001", "-outfmt", "%o %a %P"], 0),
    ("N1387", &["-db", "../seqdb/db/mask_nucl_v4 ../seqdb/db/mask_nucl_ids", "-entry", "all", "-outfmt", "%o|%a|%l|%m100"], 8),
];

#[test]
fn what_ncbi_does_not_answer_on_a_list_is_refused() {
    let config = TestDir::new("refused");
    for &(id, words, lines) in REFUSED {
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert_eq!(
            stdout.iter().filter(|&&byte| byte == b'\n').count(),
            lines,
            "{id} {words:?}"
        );
        assert_eq!(code, 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.starts_with("Error: ")
                && stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n")
                && stderr.lines().count() == 1,
            "{id} {words:?}: {stderr}"
        );
    }
}
