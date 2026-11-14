//! Where dustmasker, segmasker and blastdbcmd write the texts of a damaged taxonomy database
//! (line RG item T5b). NCBI's `CSeqDB::GetBioseq` looks up every positive taxid of each sequence
//! it reads (the maskers' reader, blastdbcmd's FASTA formatter twice per sequence), and
//! `GetTaxInfo` (`%S`) and `-tax_info` look up directly; `GetTaxNames` posts its texts through
//! `ERR_POST`, and `cerr` is tied to `cout`, so each text comes after the output written so far
//! and before the sequence or line whose lookup wrote it, where both streams go to one file
//! (`2>&1`); with `-logfile` the texts go to that file. The expectations are NCBI 2.17.0's merged
//! stream of the same run, and its log file (`tests/fixtures/taxdb_tool_texts/<case>.{txt,log}`),
//! on `tests/db_search/t5b` (make_t5b.py there) with `BLASTDB` naming
//! `tests/seqdb/taxdb/<variant>`.
//!
//! Where a lookup reads a record past the mapped end of `taxdb.btd` (the empty `taxdb.btd` of
//! `btd_empty`), NCBI crashes (SIGSEGV, every format; PD-LOSAT-NCBI-DEFECTS rule 3): LOSAT refuses
//! the run there, with one line that names its program and exit status 1.

use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// The cases: name (`<case>.<variant>`), taxonomy database (under `tests/seqdb/taxdb`), program,
/// words, and whether the run writes its diagnostics to a `-logfile` file.
const CASES: &[(&str, &str, &str, &[&str], bool)] = &[
    ("dust.bad_magic", "bad_magic", "dustmasker", DUST, false),
    ("seg.bad_magic", "bad_magic", "segmasker", SEG, false),
    ("cmd_all.bad_magic", "bad_magic", "blastdbcmd", ALL, false),
    (
        "cmd_fmtS.bad_magic",
        "bad_magic",
        "blastdbcmd",
        FMT_S,
        false,
    ),
    (
        "cmd_taxinfo.bad_magic",
        "bad_magic",
        "blastdbcmd",
        TAX_INFO,
        false,
    ),
    (
        "dust.offset_last_past_end",
        OFFSET,
        "dustmasker",
        DUST,
        false,
    ),
    ("seg.offset_last_past_end", OFFSET, "segmasker", SEG, false),
    (
        "cmd_all.offset_last_past_end",
        OFFSET,
        "blastdbcmd",
        ALL,
        false,
    ),
    (
        "cmd_fmtS.offset_last_past_end",
        OFFSET,
        "blastdbcmd",
        FMT_S,
        false,
    ),
    (
        "cmd_taxinfo.offset_last_past_end",
        OFFSET,
        "blastdbcmd",
        TAX_INFO,
        false,
    ),
    (
        "dust_log.offset_last_past_end",
        OFFSET,
        "dustmasker",
        DUST,
        true,
    ),
    (
        "seg_log.offset_last_past_end",
        OFFSET,
        "segmasker",
        SEG,
        true,
    ),
    (
        "cmd_log.offset_last_past_end",
        OFFSET,
        "blastdbcmd",
        ALL,
        true,
    ),
];

const OFFSET: &str = "offset_last_past_end";
const DUST: &[&str] = &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "fasta"];
const SEG: &[&str] = &["-in", "pmv", "-infmt", "blastdb", "-outfmt", "fasta"];
const ALL: &[&str] = &["-db", "nmv", "-entry", "all"];
const FMT_S: &[&str] = &["-db", "nmv", "-entry", "all", "-outfmt", "%a %T %S"];
const TAX_INFO: &[&str] = &["-db", "nmv", "-tax_info", "-outfmt", "%T %S"];

/// The runs on `btd_empty` where NCBI crashes (probe of line RG T5b, exit status -11 in each):
/// program and words.
const REFUSED: &[(&str, &[&str])] = &[
    (
        "dustmasker",
        &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "fasta"],
    ),
    (
        "dustmasker",
        &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "interval"],
    ),
    (
        "dustmasker",
        &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "acclist"],
    ),
    (
        "dustmasker",
        &[
            "-in",
            "nmv",
            "-infmt",
            "blastdb",
            "-outfmt",
            "seqloc_asn1_text",
        ],
    ),
    (
        "dustmasker",
        &[
            "-in",
            "nmv",
            "-infmt",
            "blastdb",
            "-outfmt",
            "maskinfo_asn1_text",
        ],
    ),
    (
        "segmasker",
        &["-in", "pmv", "-infmt", "blastdb", "-outfmt", "fasta"],
    ),
    (
        "segmasker",
        &["-in", "pmv", "-infmt", "blastdb", "-outfmt", "interval"],
    ),
    (
        "segmasker",
        &[
            "-in",
            "pmv",
            "-infmt",
            "blastdb",
            "-outfmt",
            "seqloc_asn1_text",
        ],
    ),
    (
        "segmasker",
        &[
            "-in",
            "pmv",
            "-infmt",
            "blastdb",
            "-outfmt",
            "maskinfo_asn1_text",
        ],
    ),
    ("blastdbcmd", &["-db", "nmv", "-entry", "all"]),
    (
        "blastdbcmd",
        &["-db", "nmv", "-entry", "all", "-long_seqids"],
    ),
    (
        "blastdbcmd",
        &["-db", "pmv", "-entry", "all", "-long_seqids"],
    ),
    (
        "blastdbcmd",
        &["-db", "nmv", "-entry", "all", "-outfmt", "%a %T %S"],
    ),
    (
        "blastdbcmd",
        &["-db", "nmv", "-entry", "all", "-outfmt", "%a %C %N"],
    ),
    (
        "blastdbcmd",
        &["-db", "nmv", "-entry", "AB410003.1,AB410002.1"],
    ),
    (
        "blastdbcmd",
        &["-db", "nmv", "-taxids", "9606", "-no_taxid_expansion"],
    ),
    (
        "blastdbcmd",
        &["-db", "nmv", "-tax_info", "-outfmt", "%T %S"],
    ),
];

/// A fresh directory of the target's temporary space for one run, removed when dropped.
struct RunDir(std::path::PathBuf);

impl RunDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-taxdb-tool-texts-{}-{name}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        RunDir(dir)
    }
}

impl Drop for RunDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// `losat <program> words...` in `tests/db_search/t5b` with `BLASTDB` naming the taxonomy
/// database `taxdb`, standard output and standard error into one file of `dir`: that file's
/// bytes and the exit status. No LOSAT, NCBI or BLAST variable, no registry file.
fn merged(dir: &RunDir, taxdb: &str, program: &str, words: &[&str]) -> (Vec<u8>, i32) {
    let path = dir.0.join("merged");
    let file = std::fs::File::create(&path).unwrap();
    let mut command = Command::new(BIN);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    let status = command
        .current_dir(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("db_search")
                .join("t5b"),
        )
        .arg(program)
        .args(words)
        .stdin(Stdio::null())
        .stdout(Stdio::from(file.try_clone().unwrap()))
        .stderr(Stdio::from(file))
        .env("BLASTDB", format!("../../seqdb/taxdb/{taxdb}"))
        .env("NCBI_CONFIG_PATH", &dir.0)
        .env("HOME", &dir.0)
        .env("USERPROFILE", &dir.0)
        .status()
        .expect("run LOSAT CLI");
    (std::fs::read(&path).unwrap(), status.code().unwrap())
}

fn fixture(name: &str) -> Vec<u8> {
    std::fs::read(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("taxdb_tool_texts")
            .join(name),
    )
    .unwrap()
}

#[test]
fn the_texts_of_a_damaged_taxonomy_database_come_where_ncbi_reads_the_taxids() {
    for &(name, taxdb, program, words, logfile) in CASES {
        let dir = RunDir::new(name);
        let log = dir.0.join("log");
        let mut words = words.to_vec();
        let log_text = log.to_str().unwrap().to_string();
        if logfile {
            words.extend(["-logfile", &log_text]);
        }
        let (bytes, status) = merged(&dir, taxdb, program, &words);
        assert_eq!(status, 0, "{name}: {}", String::from_utf8_lossy(&bytes));
        assert!(
            bytes == fixture(&format!("{name}.txt")),
            "{name}:\n{}",
            String::from_utf8_lossy(&bytes)
        );
        if logfile {
            assert_eq!(
                String::from_utf8_lossy(&std::fs::read(&log).unwrap()),
                String::from_utf8_lossy(&fixture(&format!("{name}.log"))),
                "{name}: the log file"
            );
        }
    }
}

#[test]
fn a_lookup_ncbi_crashes_on_is_refused() {
    for &(program, words) in REFUSED {
        let dir = RunDir::new(&format!("refused-{program}"));
        let mut command = Command::new(BIN);
        for (name, _) in std::env::vars_os() {
            let text = name.to_string_lossy();
            if text.starts_with("LOSAT_") || text.starts_with("NCBI") || text.starts_with("BLAST") {
                command.env_remove(&name);
            }
        }
        let output = command
            .current_dir(
                std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("tests")
                    .join("db_search")
                    .join("t5b"),
            )
            .arg(program)
            .args(words)
            .stdin(Stdio::null())
            .env("BLASTDB", "../../seqdb/taxdb/btd_empty")
            .env("NCBI_CONFIG_PATH", &dir.0)
            .env("HOME", &dir.0)
            .env("USERPROFILE", &dir.0)
            .output()
            .expect("run LOSAT CLI");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(stderr.lines().count(), 1, "{program} {words:?}: {stderr}");
        assert!(
            stderr.contains(&format!(
                "is not supported by LOSAT's {}",
                program.to_uppercase()
            )) && stderr.contains("reaches past the mapped end"),
            "{program} {words:?}: {stderr}"
        );
        assert_eq!(output.status.code(), Some(1), "{program} {words:?}");
    }
}
