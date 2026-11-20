//! `losat blastdbcmd`: the line ends of the batch files (`-entry_batch`, `-taxidlist`,
//! `-ipg_batch`) against NCBI BLAST+ 2.17.0 blastdbcmd on LOSAT's platform. NCBI's
//! `NcbiGetlineEOL` ends a line at `\n` on Linux, at `\n` without the `\r` before it on
//! Windows, whose C runtime also reads the file in text mode (`\r\n` as `\n`, the end at a
//! Ctrl-Z), and at `\r` or `\n` on macOS. The expected bytes are NCBI's on each platform,
//! frozen by line BDC's probe O (`tests/blastdbcmd/README.md`); the inputs are written here,
//! byte for byte. NCBI is not run here.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
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
            "losat_blastdbcmd_line_ends_{stem}_{}_{nanos}",
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

/// The directory of `tests/blastdbcmd/eol` that holds NCBI's outputs on LOSAT's platform.
fn platform() -> &'static str {
    if cfg!(windows) {
        "windows"
    } else if cfg!(target_os = "macos") {
        "macos"
    } else {
        "linux"
    }
}

/// The batch files of probe O. `P00019.3` and `AAA00018.1` (of `db/p5`) hold no `_` or `|`,
/// so a version 5 database reads them with `NStr::ToUpper` alone, which keeps a `\r`; the
/// Seq-id parser drops it from `XP_000001.1` (`e_us`). `e_ord` is `in/crlf_d.txt` (D65). `9606` and `562` are taxids of
/// `j/lst/p5`, `1001` and `1003` PIGs of `db/m/pmg`, whose mask algorithm 21 the second field
/// of a line names.
#[rustfmt::skip]
const INPUTS: &[(&str, &[u8])] = &[
    ("e_lf", b"P00019.3\nAAA00018.1\n"),
    ("e_crlf", b"P00019.3\r\nAAA00018.1\r\n"),
    ("e_cr", b"P00019.3\rAAA00018.1\r"),
    ("e_crcrlf", b"P00019.3\r\r\nAAA00018.1\r\n"),
    ("e_lfcr", b"P00019.3\n\rAAA00018.1\n"),
    ("e_mid", b"P00019.3\rAAA00018.1\n"),
    ("e_crlast", b"P00019.3\r\nAAA00018.1\r"),
    ("e_blank", b"\r\n\r\nP00019.3\r\n \r\n\t\r\nAAA00018.1\n"),
    ("e_sub", b"P00019.3\n\x1aAAA00018.1\n"),
    ("e_subin", b"P00019.3\x1a\nAAA00018.1\n"),
    ("e_crsub", b"P00019.3\r\x1a\nAAA00018.1\n"),
    ("e_subend", b"P00019.3\nAAA00018.1\n\x1a"),
    ("e_rng", b"AAA00018.1 2-10\r\nP00019.3 3-12\n"),
    ("e_rngcr", b"AAA00018.1 2-10\rP00019.3 3-12\r\n"),
    ("e_tab", b"AAA00018.1\t2-10\r\nP00019.3\t \r\n"),
    ("e_cr2", b"P00019.3\r\rAAA00018.1\r\n"),
    ("e_crlfcr", b"P00019.3\r\n\rAAA00018.1\r\n"),
    ("e_us", b"XP_000001.1\r\nXP_000020.1\r\n"),
    ("e_ord", b"gnl|BL_ORD_ID|3\r\ngnl|BL_ORD_ID|1 2-10\r\n"),
    ("t_lf", b"9606\n562\n"),
    ("t_crlf", b"9606\r\n562\r\n"),
    ("t_cr", b"9606\r562\r"),
    ("t_crcrlf", b"9606\r\r\n562\n"),
    ("t_blank", b"\r\n9606\r\n"),
    ("t_sub", b"9606\n\x1a562\n"),
    ("t_subend", b"9606\n562\n\x1a"),
    ("p_lf", b"1001 21\n1003\n"),
    ("p_crlf", b"1001 21\r\n1003\r\n"),
    ("p_cr", b"1001 21\r1003\r"),
    ("p_crcrlf", b"1001 21\r\r\n1003\n"),
    ("p_sub", b"1001 21\n\x1a1003\n"),
];

const ENTRY: &[&str] = &["-db", "db/p5", "-entry_batch"];
const TAX: &[&str] = &[
    "-db",
    "j/lst/p5",
    "-no_taxid_expansion",
    "-outfmt",
    "%a %T",
    "-taxidlist",
];
// `%s`, not the FASTA: without a range NCBI's FASTA masks only the first third of a protein
// (approved difference 5 of `PD-LOSAT-BLASTDBCMD-DEVIATIONS`).
const PIG: &[&str] = &["-db", "db/m/pmg", "-outfmt", "%o %s", "-ipg_batch"];
const PIG_FASTA: &[&str] = &["-db", "db/m/pmg", "-ipg_batch"];

/// A case of probe O: id, the words before the file, the input, whether it is read from the
/// standard input (`-`), the words after the file.
type Case = (
    &'static str,
    &'static [&'static str],
    &'static str,
    bool,
    &'static [&'static str],
);

#[rustfmt::skip]
const CASES: &[Case] = &[
    ("O01", ENTRY, "e_lf", false, &["-outfmt", "%a"]),
    ("O02", ENTRY, "e_crlf", false, &["-outfmt", "%a"]),
    ("O03", ENTRY, "e_cr", false, &["-outfmt", "%a"]),
    ("O04", ENTRY, "e_crcrlf", false, &["-outfmt", "%a"]),
    ("O05", ENTRY, "e_lfcr", false, &["-outfmt", "%a"]),
    ("O06", ENTRY, "e_mid", false, &["-outfmt", "%a"]),
    ("O07", ENTRY, "e_crlast", false, &["-outfmt", "%a"]),
    ("O08", ENTRY, "e_blank", false, &["-outfmt", "%a"]),
    ("O09", ENTRY, "e_sub", false, &["-outfmt", "%a"]),
    ("O10", ENTRY, "e_subin", false, &["-outfmt", "%a"]),
    ("O11", ENTRY, "e_crsub", false, &["-outfmt", "%a"]),
    ("O12", ENTRY, "e_subend", false, &["-outfmt", "%a"]),
    ("O13", ENTRY, "e_rng", false, &[]),
    ("O14", ENTRY, "e_rngcr", false, &[]),
    ("O15", ENTRY, "e_tab", false, &[]),
    ("O16", ENTRY, "e_cr2", false, &["-outfmt", "%a"]),
    ("O17", ENTRY, "e_crlfcr", false, &["-outfmt", "%a"]),
    ("O18", ENTRY, "e_crlf", false, &["-outfmt", "%a", "-get_dups"]),
    ("O19", ENTRY, "e_cr", false, &["-outfmt", "%a", "-get_dups"]),
    ("O20", ENTRY, "e_sub", false, &["-outfmt", "%a", "-get_dups"]),
    ("O21", ENTRY, "e_crlf", true, &["-outfmt", "%a"]),
    ("O22", ENTRY, "e_cr", true, &["-outfmt", "%a"]),
    ("O23", ENTRY, "e_sub", true, &["-outfmt", "%a"]),
    ("O24", ENTRY, "e_crcrlf", true, &["-outfmt", "%a"]),
    ("O25", ENTRY, "e_crlf", false, &["-outfmt", "%a", "-target_only"]),
    ("O26", ENTRY, "e_us", false, &["-outfmt", "%a"]),
    ("O27", &["-db", "db/p4n", "-entry_batch"], "e_ord", false, &[]),
    ("O31", TAX, "t_lf", false, &[]),
    ("O32", TAX, "t_crlf", false, &[]),
    ("O33", TAX, "t_cr", false, &[]),
    ("O34", TAX, "t_crcrlf", false, &[]),
    ("O35", TAX, "t_blank", false, &[]),
    ("O36", TAX, "t_sub", false, &[]),
    ("O37", TAX, "t_subend", false, &[]),
    ("O38", TAX, "t_crlf", true, &[]),
    ("O39", TAX, "t_cr", true, &[]),
    ("O41", PIG, "p_lf", false, &[]),
    ("O42", PIG, "p_crlf", false, &[]),
    ("O43", PIG, "p_cr", false, &[]),
    ("O44", PIG, "p_crcrlf", false, &[]),
    ("O45", PIG, "p_sub", false, &[]),
    ("O46", PIG, "p_crlf", true, &[]),
    ("O47", PIG_FASTA, "p_cr", true, &["-outfmt", "%o %a"]),
];

/// The frozen file `name` of NCBI's outputs on this platform. On Windows NCBI writes its
/// standard output and error in text mode, `\r\n` for each `\n`; LOSAT writes `\n` on every
/// platform (`PD-NCBI-PLATFORM-VARIANCE`), so each `\r` that the text mode put before a `\n`
/// is taken away.
fn ncbi(name: &str) -> Vec<u8> {
    let bytes = fs::read(run_dir().join("eol").join(platform()).join(name)).unwrap();
    if !cfg!(windows) {
        return bytes;
    }
    let mut kept = Vec::with_capacity(bytes.len());
    for (at, &byte) in bytes.iter().enumerate() {
        if !(byte == b'\r' && bytes.get(at + 1) == Some(&b'\n')) {
            kept.push(byte);
        }
    }
    kept
}

/// NCBI's exit status of each case on this platform (the probe's `manifest.tsv`).
fn statuses() -> HashMap<String, String> {
    let text = String::from_utf8(
        fs::read(run_dir().join("eol").join(platform()).join("manifest.tsv")).unwrap(),
    )
    .unwrap();
    text.lines()
        .skip(1)
        .filter(|line| *line != "RUNDONE")
        .map(|line| {
            let mut fields = line.split('\t');
            let id = fields.next().unwrap().to_string();
            (id, fields.next().unwrap().to_string())
        })
        .collect()
}

#[test]
fn batch_files_end_their_lines_as_ncbi_does_on_this_platform() {
    let dir = TestDir::new("cases");
    let inputs: HashMap<&str, &[u8]> = INPUTS.iter().copied().collect();
    let statuses = statuses();
    assert_eq!(statuses.len(), CASES.len());
    for &(id, head, input, on_stdin, tail) in CASES {
        let file = dir.0.join(format!("{input}.txt"));
        fs::write(&file, inputs[input]).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (name, _) in std::env::vars_os() {
            let text = name.to_string_lossy();
            if text.starts_with("LOSAT_")
                || text.starts_with("NCBI_")
                || text.starts_with("BLASTDB")
            {
                command.env_remove(&name);
            }
        }
        let file_word = if on_stdin {
            "-".to_string()
        } else {
            file.to_str().unwrap().to_string()
        };
        let stdin = if on_stdin {
            Stdio::from(fs::File::open(&file).unwrap())
        } else {
            Stdio::null()
        };
        let output = command
            .current_dir(run_dir())
            .stdin(stdin)
            .arg("blastdbcmd")
            .args(head)
            .arg(&file_word)
            .args(tail)
            .output()
            .unwrap();
        let code = output.status.code().unwrap();
        let status = &statuses[id];
        if status == "255" {
            // NCBI's exception that names its source files: a taxid that is not a number
            // (line BDC DECISIONS 31-3).
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(output.stdout.is_empty(), "{id}: stdout");
            assert!(stderr.contains("LOSAT's BLASTDBCMD"), "{id}: {stderr}");
            assert_eq!(code, 1, "{id}: exit status");
            continue;
        }
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&ncbi(&format!("{id}.stdout"))),
            "{id} {input}: stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            String::from_utf8_lossy(&ncbi(&format!("{id}.stderr"))),
            "{id} {input}: stderr"
        );
        assert_eq!(code.to_string(), *status, "{id} {input}: exit status");
    }
}
