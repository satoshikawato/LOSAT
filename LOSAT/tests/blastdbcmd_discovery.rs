//! `losat blastdbcmd`: the taxonomy (`-taxids`, `-taxidlist`, `-no_taxid_expansion`,
//! `-target_only`, `-tax_info`), the metadata (`-metadata`, `-metadata_output_prefix`), the
//! lists of databases (`-list`, `-list_outfmt`, `-recursive`, `-remove_redundant_dbs`),
//! `-show_blastdb_search_path` and the IPGs (`-ipg`) against NCBI BLAST+ 2.17.0 blastdbcmd. The
//! expected bytes are NCBI's, frozen by line BDC's probe J and K on ws (`tests/blastdbcmd/README.md`);
//! `{R}` stands for the directory of the probe's runs (here `tests/blastdbcmd`, the current
//! directory of each run). NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// `fs::canonicalize` without the `\\?\` prefix that Windows puts on its result (as
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
            "losat_blastdbcmd_discovery_{stem}_{}_{nanos}",
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

/// One run of `losat blastdbcmd words...` in `cwd`: stdout, stderr, exit status. `BLASTDB` is
/// `blastdb` (a directory of `tests/blastdbcmd`) or unset; no registry file (`NCBI_CONFIG_PATH`
/// names `config` alone), no LOSAT variable.
fn blastdbcmd(
    words: &[&str],
    cwd: &Path,
    blastdb: Option<&str>,
    config: &Path,
) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    if let Some(dir) = blastdb {
        command.env("BLASTDB", run_dir().join(dir));
    }
    let output = command
        .env("NCBI_CONFIG_PATH", config)
        .current_dir(cwd)
        .arg("blastdbcmd")
        .args(words)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// The frozen file `name` of `tests/blastdbcmd/expected`, with `{R}` placed.
fn expected(name: &str) -> Vec<u8> {
    let text = fs::read(run_dir().join("expected").join(name)).unwrap();
    let text = String::from_utf8(text).unwrap();
    text.replace("{R}", run_dir().to_str().unwrap())
        .into_bytes()
}

/// The probe's cases that LOSAT answers as NCBI does: id, the taxonomy directory of `BLASTDB`,
/// whether the output names a path (the lines of `-list`, the search path, the `dbname` of
/// `-metadata`, read with the platform's separator: Unix only), words, NCBI's exit status.
/// The databases are copies of the probe's (`j/lst`: version 4 and 5, without parsed ids, in
/// three volumes, an alias of two; `j/jd`: `p4` with other dates and a title to escape;
/// `j/jtax5`: a taxonomy over the taxids of `p5`). `K13` (`-taxids 1`, which never ends in NCBI:
/// taxid 1 is its own parent) is the approved difference of DI-8 DECISIONS 10, its expected
/// output NCBI's run of the explicit list of every taxid that taxid 1 reaches (K25).
#[rustfmt::skip]
const CASES: &[(&str, Option<&str>, bool, &[&str], i32)] = &[
    ("J160", None, true, &["-list", "j/lst"], 0),
    ("J161", None, true, &["-list", "j/lst", "-list_outfmt", "%f|%p|%t|%d|%l|%n|%U|%v"], 0),
    ("J163", None, true, &["-list", "j/lst", "-dbtype", "prot", "-list_outfmt", "%f %t"], 0),
    ("J164", None, true, &["-list", "j/lst/", "-dbtype", "nucl"], 0),
    ("J169", None, false, &["-list", "j/lst", "-list_outfmt", "abc"], 1),
    ("J170", None, false, &["-list", "j/lst", "-list_outfmt", "%x"], 1),
    ("J172", None, false, &["-list", "j/lst", "-list_outfmt", "%f %"], 1),
    ("J173", None, false, &["-list", "nosuchdir", "-list_outfmt", "abc"], 1),
    ("J178", None, false, &["-list", "j/lst", "-list_outfmt", ""], 1),
    ("L06", None, false, &["-list", "nosuchdir"], 0),
    ("S38", None, false, &["-list", ""], 0),
    ("L07", None, true, &["-show_blastdb_search_path"], 0),
    ("J180", None, true, &["-db", "j/lst/p4", "-metadata"], 0),
    ("J181", None, true, &["-db", "j/lst/p5", "-metadata", "-metadata_output_prefix", "/x/y/"], 0),
    ("J184", None, true, &["-db", "j/lst/p5n", "-metadata"], 0),
    ("J185", None, true, &["-db", "j/jd/am12", "-metadata"], 0),
    ("J186", None, true, &["-db", "j/jd/pm12", "-metadata"], 0),
    ("J187", None, true, &["-db", "j/jd/d01", "-metadata"], 0),
    ("J188", None, true, &["-db", "j/jd/d1", "-metadata"], 0),
    ("J191", None, true, &["-db", "j/jd/esc", "-metadata"], 0),
    ("J192", None, true, &["-db", "j/lst/dv4", "-metadata"], 0),
    ("J197", None, true, &["-db", "j/lst/p5", "-metadata", "-metadata_output_prefix", ""], 0),
    ("J198", None, true, &["-db", "j/lst/pal", "-metadata"], 0),
    ("J146", None, false, &["-db", "j/lst/p5n", "-tax_info"], 0),
    ("J153", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "%f"], 0),
    ("J155", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "%n"], 0),
    ("J157", None, false, &["-db", "j/lst/n4", "-tax_info"], 2),
    ("J159", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "%T %%"], 0),
    ("J143", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "abc"], 1),
    ("T07", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a %T"], 0),
    ("T08", None, false, &["-db", "j/lst/p5n", "-taxids", "9606", "-no_taxid_expansion"], 1),
    ("S34", None, false, &["-db", "j/lst/p4", "-taxids", "9606", "-no_taxid_expansion"], 1),
    ("J106", None, false, &["-db", "j/lst/p5", "-taxids", "0", "-no_taxid_expansion", "-outfmt", "%o %a %T"], 0),
    ("J110", None, false, &["-db", "j/lst/p5", "-taxidlist", "in/nosuch.txt", "-no_taxid_expansion"], 1),
    ("J111", None, true, &["-db", "j/lst/p4", "-taxids", "9606", "-target_only"], 2),
    ("K09", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-dbtype", "prot", "-strand", "minus"], 0),
    ("K14", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-target_only"], 0),
    ("K17", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-target_only", "-dbtype", "prot"], 0),
    ("K10", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "9604", "-outfmt", "%o %a %T"], 0),
    ("K11", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "10066", "-target_only", "-outfmt", "%o %a %T"], 0),
    ("K12", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "9604,10066"], 0),
    ("K13", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "1", "-outfmt", "%o %T"], 0),
];

/// The probe's cases that LOSAT refuses (one line naming LOSAT's BLASTDBCMD, exit status 1,
/// nothing written): NCBI's exception that names its source files (`-taxids` with an empty id,
/// a date `CTime` does not read), NCBI's protein letters written reversed (a taxid list on a
/// protein database read with `-dbtype guess` and `-strand minus`), and the OID that NCBI's
/// `-ipg` never sets on a database without PIG files.
#[rustfmt::skip]
const REFUSED: &[(&str, &[&str])] = &[
    ("J101", &["-db", "j/lst/p5", "-taxids", "4840,,4841", "-no_taxid_expansion", "-outfmt", "%a %T"]),
    ("J114", &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-strand", "minus"]),
    ("J189", &["-db", "j/jd/feb30", "-metadata"]),
    ("J190", &["-db", "j/jd/bad", "-metadata"]),
    ("J204", &["-db", "j/lst/p5", "-ipg", "5"]),
];

#[test]
fn taxonomy_metadata_and_lists_are_ncbi_s_bytes() {
    let config = TestDir::new("cases");
    for &(id, blastdb, unix_only, words, status) in CASES {
        if unix_only && cfg!(windows) {
            continue;
        }
        let (stdout, stderr, code) = blastdbcmd(words, &run_dir(), blastdb, &config.0);
        assert_eq!(
            String::from_utf8_lossy(&stdout),
            String::from_utf8_lossy(&expected(&format!("{id}.stdout"))),
            "{id} {words:?}: stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            String::from_utf8_lossy(&expected(&format!("{id}.stderr"))),
            "{id} {words:?}: stderr"
        );
        assert_eq!(code, status, "{id} {words:?}: exit status");
    }
}

#[test]
fn what_ncbi_does_not_answer_is_refused() {
    let config = TestDir::new("refused");
    for &(id, words) in REFUSED {
        let (stdout, stderr, code) = blastdbcmd(words, &run_dir(), None, &config.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert!(stdout.is_empty(), "{id} {words:?}: stdout");
        assert_eq!(stderr.lines().count(), 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.contains("LOSAT's BLASTDBCMD"),
            "{id} {words:?}: {stderr}"
        );
        assert_eq!(code, 1, "{id} {words:?}: exit status");
    }
}

/// The entries of the probe's `db/j/lsx`: files the walk counts without opening them, a
/// directory named like a database, links to a file, to nowhere and to a directory, names that
/// sort apart from their bytes once quoted.
#[cfg(unix)]
fn make_lsx(dir: &Path) {
    use std::os::unix::fs::symlink;
    fs::create_dir_all(dir.join("sub").join("deep")).unwrap();
    fs::create_dir_all(dir.join("d.pin")).unwrap();
    for name in [
        "a b.pin", "a!.nin", ".h.pin", "X.PIN", "zz.pin", "n5.nin", "n5x.nin", "n55.nin",
    ] {
        fs::write(dir.join(name), b"").unwrap();
    }
    fs::write(dir.join("zz.pal"), b"DBLIST zz\n").unwrap();
    fs::write(dir.join("q.nal"), b"TITLE q\nDBLIST n5 \"n5x\"\n").unwrap();
    fs::write(dir.join("sub").join("y.nin"), b"").unwrap();
    fs::write(dir.join("sub").join("deep").join("z.pin"), b"").unwrap();
    fs::write(dir.join("d.pin").join("inner.pin"), b"").unwrap();
    symlink("zz.pin", dir.join("ln.pin")).unwrap();
    symlink("nowhere.nin", dir.join("br.nin")).unwrap();
    symlink("sub", dir.join("lsub")).unwrap();
}

// Probes J165-J167 and J176 of line BDC on a copy of `db/j/lsx` built here.
#[cfg(unix)]
#[test]
fn the_walk_lists_files_links_and_directories_as_ncbi_s() {
    let work = TestDir::new("lsx");
    let lsx = work.0.join("lsx");
    make_lsx(&lsx);
    let cases: [(&str, &[&str]); 4] = [
        ("J165", &["-list", "lsx"]),
        ("J166", &["-list", "lsx", "-recursive"]),
        ("J167", &["-list", "lsx", "-remove_redundant_dbs"]),
        (
            "J176",
            &[
                "-list",
                "lsx",
                "-recursive",
                "-dbtype",
                "nucl",
                "-remove_redundant_dbs",
                "-list_outfmt",
                "%p:%f",
            ],
        ),
    ];
    for (id, words) in cases {
        let (stdout, stderr, code) = blastdbcmd(words, &work.0, None, &work.0);
        let frozen = fs::read(run_dir().join("expected").join(format!("{id}.stdout"))).unwrap();
        let want = String::from_utf8(frozen)
            .unwrap()
            .replace("db/j/lsx", "lsx");
        assert_eq!(String::from_utf8_lossy(&stdout), want, "{id} {words:?}");
        assert!(
            stderr.is_empty(),
            "{id}: {}",
            String::from_utf8_lossy(&stderr)
        );
        assert_eq!(code, 0, "{id}");
    }
}
