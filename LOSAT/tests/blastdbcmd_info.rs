//! `losat blastdbcmd`: the `-info` of a database, the argument errors, `-out` and `-logfile`,
//! against NCBI BLAST+ 2.17.0 blastdbcmd. The expected bytes are NCBI's, frozen by line BDC's
//! probe on ws (`tests/blastdbcmd/README.md`); `{R}` stands for the directory of the probe's
//! runs (here `tests/blastdbcmd`, the current directory of each run) and `{T}` for `tests/`.
//! NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// `fs::canonicalize` without the `\\?\` prefix that Windows puts on its result: a verbatim path
/// takes `/` as a character of a name (`{O}/B23.out` below is not a file of the directory then),
/// and the product prints the paths it lists without the prefix.
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

fn tests_dir() -> PathBuf {
    canonical(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests"))
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
            "losat_blastdbcmd_{stem}_{}_{nanos}",
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

/// One run of `losat blastdbcmd words...` in `tests/blastdbcmd`, `{O}` in a word standing for
/// `out_dir`: stdout, stderr, exit status. No `BLASTDB`, no registry file (`NCBI_CONFIG_PATH`
/// names `out_dir` alone), no LOSAT variable.
fn blastdbcmd(words: &[&str], out_dir: &Path) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    let out_dir = out_dir.to_str().unwrap();
    let output = command
        .env("NCBI_CONFIG_PATH", out_dir)
        .current_dir(tests_dir().join("blastdbcmd"))
        .arg("blastdbcmd")
        .args(words.iter().map(|word| word.replace("{O}", out_dir)))
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// The frozen NCBI file `name` of `tests/blastdbcmd/expected`, with `{R}` and `{T}` placed.
fn expected(name: &str) -> Vec<u8> {
    let tests = tests_dir();
    let text = fs::read(tests.join("blastdbcmd/expected").join(name)).unwrap();
    let text = String::from_utf8(text).unwrap();
    native_paths(text, &tests)
}

/// The frozen text with `{R}` and `{T}` placed. NCBI's files come from Linux; on Windows the
/// paths in them are the ones NCBI's code makes there (`SeqDB_ConvertOSPath`: the separator of
/// the platform; the search path is split at `;`), so the rest of each placed path is written
/// with `\`, so is the name in `database [name]`, and the `:` of the search path `[{R}::]` is `;`.
#[cfg(not(windows))]
fn native_paths(text: String, tests: &Path) -> Vec<u8> {
    text.replace("{T}", tests.to_str().unwrap())
        .replace("{R}", tests.join("blastdbcmd").to_str().unwrap())
        .into_bytes()
}

#[cfg(windows)]
fn native_paths(text: String, tests: &Path) -> Vec<u8> {
    let text = text.replace("{R}::]", "{R};;]");
    // The database name of "No alias or index file found for ... database [name]" goes through
    // the same conversion (`SeqDB_ConvertOSPath`).
    let mut converted = String::new();
    let mut rest = text.as_str();
    while let Some(at) = rest.find("database [") {
        let name_at = at + "database [".len();
        converted.push_str(&rest[..name_at]);
        let end = rest[name_at..]
            .find(']')
            .map_or(rest.len(), |n| name_at + n);
        converted.push_str(&rest[name_at..end].replace('/', "\\"));
        rest = &rest[end..];
    }
    converted.push_str(rest);
    let text = converted;
    let mut out = String::new();
    let mut rest = text.as_str();
    while let Some(at) = rest.find('{') {
        out.push_str(&rest[..at]);
        let root = match &rest[at..at + 3] {
            "{R}" => tests.join("blastdbcmd"),
            "{T}" => tests.to_path_buf(),
            _ => unreachable!("a placeholder of the frozen files"),
        };
        out.push_str(root.to_str().unwrap());
        rest = &rest[at + 3..];
        let end = rest
            .find(|c: char| !(c.is_ascii_alphanumeric() || "_./-".contains(c)))
            .unwrap_or(rest.len());
        out.push_str(&rest[..end].replace('/', "\\"));
        rest = &rest[end..];
    }
    out.push_str(rest);
    out.into_bytes()
}

/// The probe's cases that LOSAT answers as NCBI does: id, words, NCBI's exit status. The
/// database `fx/LOSAT/tests/seqdb/db/X` of the probe is `../seqdb/db/X` here.
const MATCHES: &[(&str, &[&str], i32)] = &[
    ("I06", &["-db", "db/p5", "-info"], 0),
    ("I07", &["-db", "db/p4", "-info"], 0),
    ("I13", &["-db", "db/n5n", "-info"], 0),
    ("I15", &["-db", "db/nL", "-info"], 0),
    ("I22", &["-db", "../seqdb/db/alias_nested", "-info"], 0),
    ("I24", &["-db", "../seqdb/db/prot_mv", "-info"], 0),
    ("I25", &["-db", "../seqdb/db/empty", "-info"], 0),
    ("I26", &["-db", "db/p5", "-dbtype", "prot", "-info"], 0),
    ("I27", &["-db", "db/p5", "-dbtype", "nucl", "-info"], 2),
    ("I29", &["-db", "db/nosuch", "-info"], 2),
    ("I30", &["-info"], 2),
    (
        "I33",
        &["-db", "db/p5", "-info", "-entry", "XP_000001.1"],
        1,
    ),
    ("I34", &["-db", "./db/p5", "-info"], 0),
    ("I35", &["-db", "db/p5.pin", "-info"], 2),
    ("B01", &["-db", "db/p5", "-info", "-foo"], 1),
    ("B02", &["-db"], 1),
    ("B03", &["-db", "db/p5", "-info", "-dbtype", "protein"], 1),
    ("B04", &["-db", "db/p5", "-info", "-db", "db/n5"], 1),
    ("B05", &["-db", "db/p5", "-exact_length"], 1),
    (
        "B06",
        &["-db", "db/p5", "-info", "-out", "/nonexistent/dir/x.txt"],
        1,
    ),
    ("B09", &["-db", "db/p5", "-info", "extra"], 1),
    ("B10", &["-db", "db/p5", "-info=1"], 1),
    ("B11", &["-db=db/p5", "-info"], 0),
    ("B12", &["-db", "", "-info"], 2),
    ("B14", &["-db", "db/p5", "-info", "-outfmt", "%f"], 1),
    ("B16", &["-db", "db/p5", "-info", "-line_length", "abc"], 1),
    ("B17", &["-db", "db/p5", "-info", "-h"], 0),
    ("B26", &["-db", "  db/p5  ", "-info"], 0),
    (
        "B29",
        &["-db", "db/p5", "-info", "-exact_length", "-exact_length"],
        1,
    ),
    ("B30", &["-db", "db/p5", "-info", "-help"], 0),
    ("B32", &["-db", "db/p5", "-info", "--", "-x"], 1),
    (
        "B49",
        &["-db", "db/p5", "-info", "-metadata_output_prefix", "x"],
        1,
    ),
    ("B52", &["-db", "db/p5", "-info", "-line_length", "60"], 0),
    ("B57", &["-db", "\"db/p5\"", "-info"], 0),
    ("B58", &["-db", "db/p5", "-info", "-out", ""], 1),
    (
        "B61",
        &["-db", "db/p5", "-info", "-line_length", "2147483648"],
        1,
    ),
    // `-logfile -`: the diagnostics on stderr, no file; `-logfile ""` and `/dev/null`: none.
    ("C01", &["-db", "db/nosuch", "-info", "-logfile", "-"], 2),
    (
        "C02",
        &["-db", "db/p5", "-info", "-foo", "-logfile", "-"],
        1,
    ),
    ("C03", &["-db", "db/nosuch", "-info", "-logfile", ""], 2),
    (
        "C04",
        &["-db", "db/p5", "-info", "-foo", "-logfile", "/dev/null"],
        1,
    ),
    (
        "C19",
        &["-db", "db/nosuch", "-info", "-logfile=", "-foo"],
        1,
    ),
    ("C05", &["-db", "db/p5", "-info", "-help-full"], 0),
    ("C07", &["-db", "db/p5", "-info", "-version=1"], 1),
    ("C08", &["-db", "db/p5", "-info", "-foo.bar"], 1),
    ("C09", &["-db", "db/p5", "-info-file"], 0),
    // blastdbcmd without `-info` and without an entry: NCBI opens the database, then asks
    // for the query type (`blastdbcmd_retrieve.rs` has the retrieval).
    ("B40", &["-db", "db/p5"], 1),
    ("C10", &["-db-file", "db/p5", "-info"], 0),
    (
        "C11",
        &["-db", "db/p5", "-info", "-line_length", "\u{3042}"],
        1,
    ),
    ("C12", &["-db", "db/p5", "-info", "-line_length", "1\n2"], 1),
    ("C13", &["-db", "db/p5", "-info", "---x"], 1),
    ("C14", &["-db", "db/p5", "-info", "-=x"], 1),
    ("C15", &["-db", "db/p5", "-info", "--x"], 1),
    (
        "C18",
        &[
            "-db",
            "db/p5",
            "-info",
            "-dbtype",
            "prot",
            "-dbtype-file",
            "nucl",
        ],
        1,
    ),
];

#[test]
fn info_and_argument_errors_are_ncbi_s_bytes() {
    let dir = TestDir::new("match");
    for &(id, words, status) in MATCHES {
        // C04 sends the diagnostics to `/dev/null`; Windows has no such file, so there `-logfile`
        // is "not accessible" and NCBI's answer on Windows is not measured.
        if cfg!(windows) && id == "C04" {
            continue;
        }
        let (stdout, stderr, code) = blastdbcmd(words, &dir.0);
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
    // `-logfile -` is stderr, not a file named `-` in the current directory.
    assert!(!tests_dir().join("blastdbcmd/-").exists());
}

// NCBI opens `-out` before the database (B43: the file is there and empty after the database
// error), and appends its diagnostics to the file of `-logfile` (B08 then B38 into one file).
#[test]
fn out_and_logfile_take_what_ncbi_writes_there() {
    let dir = TestDir::new("files");
    let run = |words: &[&str]| blastdbcmd(words, &dir.0);
    let file = |name: &str| fs::read(dir.0.join(name)).unwrap();

    let (stdout, stderr, code) = run(&["-db", "db/p5", "-info", "-out", "{O}/B23.out"]);
    assert_eq!((stdout, stderr, code), (vec![], vec![], 0), "B23");
    assert_eq!(file("B23.out"), expected("I06.stdout"), "B23 -out");

    let (stdout, stderr, code) = run(&["-db", "db/nosuch", "-info", "-out", "{O}/B43.out"]);
    assert_eq!(
        (stdout, stderr, code),
        (vec![], expected("B43.stderr"), 2),
        "B43"
    );
    assert_eq!(file("B43.out"), b"", "B43 -out");

    let (stdout, stderr, code) = run(&[
        "-db",
        "db/p5",
        "-info",
        "-logfile",
        "{O}/B07.logfile",
        "-foo",
    ]);
    assert_eq!(
        (stdout, String::from_utf8_lossy(&stderr), code),
        (vec![], String::from_utf8_lossy(&expected("B07.stderr")), 1),
        "B07"
    );
    assert_eq!(file("B07.logfile"), expected("B07.logfile"), "B07 -logfile");

    let first = run(&["-db", "db/nosuch", "-info", "-logfile", "{O}/B08.logfile"]);
    assert_eq!(first, (vec![], vec![], 2), "B08");
    let second = run(&["-db", "db/nosuch2", "-info", "-logfile={O}/B08.logfile"]);
    assert_eq!(second, (vec![], vec![], 2), "B38");
    assert_eq!(
        file("B08.logfile"),
        expected("B08.logfile"),
        "B08 + B38 -logfile"
    );
}

/// The probe's cases that LOSAT's BLASTDBCMD does not serve (line BDC DECISIONS 6): several
/// names in `-db` (I38), the toolkit words LOSAT refuses.
const REFUSALS: &[(&str, &[&str])] = &[
    ("I38", &["-db", "db/p5 db/pm5", "-info"]),
    ("B20", &["-db", "db/p5", "-info", "-version-full"]),
    ("B21", &["-db", "db/p5", "-info", "-conffile", "x.ini"]),
    ("B22", &["-db", "db/p5", "-info", "-dryrun"]),
    ("C06", &["-db", "db/p5", "-info", "-xmlhelp"]),
];

#[test]
fn what_blastdbcmd_does_not_serve_is_refused_without_output() {
    let dir = TestDir::new("refuse");
    for &(id, words) in REFUSALS {
        let (stdout, stderr, code) = blastdbcmd(words, &dir.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert_eq!(stdout, b"", "{id} {words:?}");
        assert_eq!(code, 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.starts_with("Error: ")
                && stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n")
                && stderr.lines().count() == 1,
            "{id} {words:?}: {stderr}"
        );
    }
}

// line CLI DECISIONS 4: NCBI's first line, then LOSAT's package line; `-version` is read
// anywhere before `--`, also in the place of a value and before an argument error.
#[test]
fn version_is_read_before_the_arguments() {
    let dir = TestDir::new("version");
    let version = format!(
        "blastdbcmd: 2.17.0+\n Package: LOSAT {}, output compatible with blast 2.17.0\n",
        env!("CARGO_PKG_VERSION")
    );
    for words in [
        &["-version"][..],
        &["-db", "-version"],
        &["-foo", "-version"],
    ] {
        let (stdout, stderr, code) = blastdbcmd(words, &dir.0);
        assert_eq!(
            (String::from_utf8(stdout).unwrap(), stderr, code),
            (version.clone(), vec![], 0),
            "{words:?}"
        );
    }
    let (stdout, _, code) = blastdbcmd(&["-db", "db/p5", "-info", "--", "-version"], &dir.0);
    assert_eq!((stdout, code), (vec![], 1), "-version after --");
}
