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
//! Where a lookup reads a record of an empty `taxdb.btd` (`btd_empty`, which NCBI maps as a null
//! pointer), NCBI crashes (SIGSEGV, every format that looks a taxid up; PD-LOSAT-NCBI-DEFECTS rule
//! 3): LOSAT refuses the run there, with one line that names its program and exit status 1. Where a
//! record of another `taxdb.btd` reaches past NCBI's mapping of the file
//! (`tests/db_search/t5b/taxdb/<variant>`, make_t5b.py --taxdb), NCBI reads memory past its mapping
//! and LOSAT refuses only a written name that comes from there (the same line and status).
//!
//! With `-out /dev/stdout` and standard error in the same pipe (line RG item T8a, Linux), the
//! report is NCBI's file stream (an 8,191-byte buffer, its writes at a flush or where an insertion
//! does not fit) and no tie makes it write before a text: each text lands after NCBI's last write
//! of the report, which can be inside a record or a line (`PIPE_CASES`).

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

/// The runs that write a name NCBI reads past its mapping of `taxdb.btd` (line RG item T5b fix 2;
/// NCBI 2.17.0 completes them, with bytes that depend on its memory map): taxonomy database under
/// `tests/db_search/t5b/taxdb`, program and words.
const PAST_MAPPING: &[(&str, &str, &[&str])] = &[
    (
        "kingdom_past_page",
        "blastdbcmd",
        &["-db", "nmv", "-entry", "all", "-outfmt", "%a %K"],
    ),
    (
        "kingdom_past_page",
        "blastdbcmd",
        &["-db", "nmv", "-tax_info", "-outfmt", "%T|%K"],
    ),
    (
        "kingdom_past_page",
        "blastdbcmd",
        &["-db", "nmv", "-tax_info"],
    ),
    (
        "names_past_page",
        "blastdbcmd",
        &["-db", "nmv", "-entry", "all", "-outfmt", "%a %T %S"],
    ),
    (
        "names_past_page",
        "blastdbcmd",
        &["-db", "nmv", "-entry", "all", "-outfmt", "%a %C %N"],
    ),
    (
        "names_past_page",
        "blastdbcmd",
        &["-db", "pmv", "-entry", "all", "-long_seqids"],
    ),
    (
        "names_past_page",
        "blastdbcmd",
        &["-db", "nmv", "-tax_info", "-outfmt", "%T|%S"],
    ),
    (
        "names_past_page",
        "dustmasker",
        &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "fasta"],
    ),
    (
        "names_past_page",
        "segmasker",
        &["-in", "pmv", "-infmt", "blastdb", "-outfmt", "fasta"],
    ),
];

/// One run of `losat program words...` in `tests/db_search/t5b` with `BLASTDB=blastdb`: one line of
/// standard error that names LOSAT's program and holds `reason`, exit status 1.
fn assert_refused(program: &str, words: &[&str], blastdb: &str, reason: &str) {
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
        .env("BLASTDB", blastdb)
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
        )) && stderr.contains(reason),
        "{program} {words:?}: {stderr}"
    );
    assert_eq!(output.status.code(), Some(1), "{program} {words:?}");
}

#[test]
fn a_name_read_past_ncbi_s_mapping_is_refused_where_it_is_written() {
    for &(variant, program, words) in PAST_MAPPING {
        assert_refused(
            program,
            words,
            &format!("taxdb/{variant}"),
            "NCBI BLAST+'s mapping",
        );
    }
}

#[test]
fn a_lookup_ncbi_crashes_on_is_refused() {
    for &(program, words) in REFUSED {
        assert_refused(
            program,
            words,
            "../../seqdb/taxdb/btd_empty",
            "reaches past the mapped end",
        );
    }
}

/// The cases of a report written to `-out /dev/stdout` (line RG item T8a, Linux): name
/// (`<mode>.<case>.<variant>`), directory under `tests/db_search`, taxonomy database (under
/// `tests/seqdb/taxdb`), program, words and NCBI's exit status. Mode `pipe`: standard output a
/// pipe, standard error into it, and the bytes read from it are NCBI's (`<name>.txt`); `split`: no
/// `-out`, standard output and standard error two pipes (`<name>.out`, `<name>.err`: windowmasker's
/// texts of a taxonomy database that does not open come once, W1); `closed`: standard output a pipe
/// whose reader closed before the start, where SIGPIPE ends NCBI's process at its first write of the
/// report: LOSAT's standard error is NCBI's (`<name>.err`, the texts before that write) and its exit
/// status 0 (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10). The inputs of `t8`:
/// `tests/db_search/t8/make_t8.py`; `long_*` read records longer than NCBI's buffer, so that the
/// offset error lands inside a record; `unit_*` lands where a piece of 1,024 letters of NCBI's
/// sequence reader does not fit the buffer; `long_cmd_ll1` where a line feed (a `put`) does not.
#[cfg(target_os = "linux")]
const PIPE_CASES: &[(&str, &str, &str, &str, &[&str], i32)] = &[
    (
        "pipe.dust_fasta.bad_magic",
        "t5b",
        "bad_magic",
        "dustmasker",
        DUST,
        0,
    ),
    (
        "pipe.dust_interval.bad_magic",
        "t5b",
        "bad_magic",
        "dustmasker",
        DUST_INTERVAL,
        0,
    ),
    (
        "pipe.dust_acclist.bad_magic",
        "t5b",
        "bad_magic",
        "dustmasker",
        DUST_ACCLIST,
        0,
    ),
    (
        "pipe.seg_fasta.bad_magic",
        "t5b",
        "bad_magic",
        "segmasker",
        SEG,
        0,
    ),
    (
        "pipe.cmd_all.bad_magic",
        "t5b",
        "bad_magic",
        "blastdbcmd",
        ALL,
        0,
    ),
    (
        "pipe.cmd_fmtS.bad_magic",
        "t5b",
        "bad_magic",
        "blastdbcmd",
        FMT_S,
        0,
    ),
    (
        "pipe.dust_fasta.offset_last_past_end",
        "t5b",
        OFFSET,
        "dustmasker",
        DUST,
        0,
    ),
    (
        "pipe.dust_interval.offset_last_past_end",
        "t5b",
        OFFSET,
        "dustmasker",
        DUST_INTERVAL,
        0,
    ),
    (
        "pipe.dust_acclist.offset_last_past_end",
        "t5b",
        OFFSET,
        "dustmasker",
        DUST_ACCLIST,
        0,
    ),
    (
        "pipe.seg_fasta.offset_last_past_end",
        "t5b",
        OFFSET,
        "segmasker",
        SEG,
        0,
    ),
    (
        "pipe.cmd_all.offset_last_past_end",
        "t5b",
        OFFSET,
        "blastdbcmd",
        ALL,
        0,
    ),
    (
        "pipe.cmd_fmtS.offset_last_past_end",
        "t5b",
        OFFSET,
        "blastdbcmd",
        FMT_S,
        0,
    ),
    (
        "pipe.cmd_nf_entries.good",
        "t5b",
        "good",
        "blastdbcmd",
        NOT_FOUND,
        1,
    ),
    (
        "pipe.cmd_nf_batch.good",
        "t5b",
        "good",
        "blastdbcmd",
        &["-db", "nmv", "-entry_batch", "../t8/bt.txt"],
        1,
    ),
    (
        "pipe.wm_counts.bad_magic",
        "t5b",
        "bad_magic",
        "windowmasker",
        MK_COUNTS,
        0,
    ),
    (
        "split.wm_counts.bad_magic",
        "t5b",
        "bad_magic",
        "windowmasker",
        MK_COUNTS,
        0,
    ),
    (
        "pipe.wm_counts.count_zero",
        "t5b",
        "count_zero",
        "windowmasker",
        MK_COUNTS,
        0,
    ),
    (
        "split.wm_counts.count_zero",
        "t5b",
        "count_zero",
        "windowmasker",
        MK_COUNTS,
        0,
    ),
    (
        "pipe.c2bm.none",
        "t8",
        "good",
        "convert2blastmask",
        &[
            "-in",
            "wq.fna",
            "-masking_algorithm",
            "x",
            "-masking_options",
            "y",
        ],
        0,
    ),
    (
        "closed.cmd_nf_entries.good",
        "t5b",
        "good",
        "blastdbcmd",
        NOT_FOUND,
        0,
    ),
    (
        "closed.dust_acclist.offset_last_past_end",
        "t5b",
        OFFSET,
        "dustmasker",
        DUST_ACCLIST,
        0,
    ),
    (
        "pipe.long_cmd_all.offset_last_past_end",
        "t8",
        OFFSET,
        "blastdbcmd",
        &["-db", "nl", "-entry", "all"],
        0,
    ),
    (
        "pipe.long_cmd_ll1000.offset_last_past_end",
        "t8",
        OFFSET,
        "blastdbcmd",
        &["-db", "nl", "-entry", "all", "-line_length", "1000"],
        0,
    ),
    (
        "pipe.long_dust_fasta.offset_last_past_end",
        "t8",
        OFFSET,
        "dustmasker",
        &["-in", "nl", "-infmt", "blastdb", "-outfmt", "fasta"],
        0,
    ),
    (
        "pipe.unit_cmd_ll1000.offset_last_past_end",
        "t8",
        OFFSET,
        "blastdbcmd",
        &["-db", "nu", "-entry", "all", "-line_length", "1000"],
        0,
    ),
    (
        "pipe.long_cmd_range.offset_last_past_end",
        "t8",
        OFFSET,
        "blastdbcmd",
        &[
            "-db", "nl", "-entry", "all", "-range", "100-9000", "-strand", "minus",
        ],
        0,
    ),
    (
        "pipe.long_cmd_ll1.offset_last_past_end",
        "t8",
        OFFSET,
        "blastdbcmd",
        &["-db", "nl", "-entry", "all", "-line_length", "1"],
        0,
    ),
];

#[cfg(target_os = "linux")]
const DUST_INTERVAL: &[&str] = &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "interval"];
#[cfg(target_os = "linux")]
const DUST_ACCLIST: &[&str] = &["-in", "nmv", "-infmt", "blastdb", "-outfmt", "acclist"];
#[cfg(target_os = "linux")]
const NOT_FOUND: &[&str] = &[
    "-db",
    "nmv",
    "-entry",
    "AB410002.1,XX000001.1,AB410003.1,YY000002.1,AB410004.1",
];
#[cfg(target_os = "linux")]
const MK_COUNTS: &[&str] = &["-mk_counts", "-in", "nmv", "-infmt", "blastdb"];

/// The run of a `PIPE_CASES` case `name` (its mode the first word): the bytes of standard output
/// (`pipe`: standard output and standard error), of standard error, and the exit status.
#[cfg(target_os = "linux")]
fn piped(
    dir: &RunDir,
    name: &str,
    sub: &str,
    taxdb: &str,
    program: &str,
    words: &[&str],
) -> (Vec<u8>, Vec<u8>, i32) {
    use std::io::Read;
    let mode = name.split('.').next().unwrap();
    let mut command = Command::new(BIN);
    for (variable, _) in std::env::vars_os() {
        let text = variable.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI") || text.starts_with("BLAST") {
            command.env_remove(&variable);
        }
    }
    command
        .current_dir(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("db_search")
                .join(sub),
        )
        .arg(program)
        .args(words)
        .stdin(Stdio::null())
        .env("BLASTDB", format!("../../seqdb/taxdb/{taxdb}"))
        .env("NCBI_CONFIG_PATH", &dir.0)
        .env("HOME", &dir.0)
        .env("USERPROFILE", &dir.0);
    match mode {
        "pipe" => {
            let (mut reader, writer) = std::io::pipe().unwrap();
            command
                .args(["-out", "/dev/stdout"])
                .stdout(writer.try_clone().unwrap())
                .stderr(writer);
            let mut child = command.spawn().expect("run LOSAT CLI");
            // The command holds the pipe's writing ends: without them the read ends with the run.
            drop(command);
            let mut bytes = Vec::new();
            reader.read_to_end(&mut bytes).unwrap();
            (bytes, Vec::new(), child.wait().unwrap().code().unwrap())
        }
        "closed" => {
            let (reader, writer) = std::io::pipe().unwrap();
            drop(reader);
            let output = command
                .args(["-out", "/dev/stdout"])
                .stdout(writer)
                .stderr(Stdio::piped())
                .output()
                .expect("run LOSAT CLI");
            (Vec::new(), output.stderr, output.status.code().unwrap())
        }
        _ => {
            let output = command.output().expect("run LOSAT CLI");
            (output.stdout, output.stderr, output.status.code().unwrap())
        }
    }
}

#[cfg(target_os = "linux")]
#[test]
fn texts_land_after_ncbi_s_last_write_of_the_out_file() {
    for &(name, sub, taxdb, program, words, status) in PIPE_CASES {
        let dir = RunDir::new(name);
        let (out, err, code) = piped(&dir, name, sub, taxdb, program, words);
        assert_eq!(code, status, "{name}: {}", String::from_utf8_lossy(&err));
        let (want_out, want_err) = match name.split('.').next().unwrap() {
            "pipe" => (fixture(&format!("{name}.txt")), Vec::new()),
            "closed" => (Vec::new(), fixture(&format!("{name}.err"))),
            _ => (
                fixture(&format!("{name}.out")),
                fixture(&format!("{name}.err")),
            ),
        };
        assert!(
            out == losat::wording::losat_output_bytes(&want_out).as_ref(),
            "{name}: standard output\n{}",
            String::from_utf8_lossy(&out)
        );
        assert_eq!(
            String::from_utf8_lossy(&err),
            String::from_utf8_lossy(&losat::wording::losat_bytes(&want_err)),
            "{name}: standard error"
        );
    }
}
