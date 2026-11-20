//! A `-db` name with a byte >= 0x80 in the XML family (outfmt 5, 13, 14, 15, 16): the name reaches the
//! document raw (`BlastOutput_db` of outfmt 5, `search-target` of XML2 and JSON), where NCBI BLAST+
//! 2.17.0's US-ASCII XML stream fails at the byte and its JSON stream copies it. LOSAT writes what NCBI
//! writes, for the five programs (`tests/fixtures/xml_dbname/<program>_<outfmt>.{stdout,stderr,rc}`, and
//! for outfmt 13 and 14 the master file `<program>_<outfmt>.out` and its side files `.out_<n>.<ext>`):
//!
//! - outfmt 5: the header of the first query fails, nothing on standard output, `Error: NCBI C++ Exception:`
//!   and its trace on standard error, exit status 255;
//! - outfmt 16: the document up to the byte, `Error: [prog]  Frame type= eFrameNamed, Object type=
//!   BlastOutput2NCBI C++ Exception:` and its trace, exit status 0;
//! - outfmt 14: the side file cut at the byte, the master file empty, the trace, exit status 255;
//! - outfmt 15 and 13: the name as it is (invalid UTF-8 or not), exit status 0.
//!
//! The runs are made from a directory that holds `q.fa` and the database `xb_caf\u{e9}` (a copy of
//! `tests/db_search/xb_p` or `xb_n`, the inputs of the rows `db.xb.*`), with the names as the command line
//! gives them, so that the document does not hold the path of a temporary directory. NCBI's bytes were
//! made by the same command with `-num_threads 1`. NCBI is the comparison oracle only.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:556-583 (`WriteEncodedChar`: a byte >= 0x80 of
//! the string goes through `CUtf8::SymbolToChar` into US-ASCII), c++/src/serial/objostrjson.cpp:756-769.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

mod temp_support;
use temp_support::{home_and_worker_log, unique_dir};

/// The program, the query file of `tests/db_search/xb`, and the database (`p` or `n`) it searches.
const PROGRAMS: [(&str, &str, &str); 5] = [
    ("blastp", "xb_qp_clean.faa", "p"),
    ("blastn", "xb_qn_clean.fna", "n"),
    ("blastx", "xb_qn_clean.fna", "p"),
    ("tblastn", "xb_qp_clean.faa", "n"),
    ("tblastx", "xb_qn_clean.fna", "n"),
];

fn engine() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

/// NCBI's frozen file `name`, in LOSAT's words (NCBI's bytes through `losat::wording`, AGENTS.md
/// rule 8).
fn frozen(name: &str) -> Vec<u8> {
    let path = engine().join("tests/fixtures/xml_dbname").join(name);
    let bytes = fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    losat::wording::losat_output_bytes(&bytes).into_owned()
}

/// The frozen files of `case` that a run of outfmt 13 or 14 leaves beside its master file: the master
/// `<case>.out` and `<case>.out_<n>.<ext>`, sorted.
fn frozen_files(case: &str) -> Vec<String> {
    let mut found: Vec<String> = fs::read_dir(engine().join("tests/fixtures/xml_dbname"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .filter(|name| name.starts_with(&format!("{case}.out")))
        .collect();
    found.sort();
    found
}

fn names(dir: &Path) -> Vec<String> {
    let mut found: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    found.sort();
    found
}

#[test]
fn a_database_name_outside_us_ascii_is_written_as_ncbi_writes_it() {
    for (program, query, database) in PROGRAMS {
        let extensions = if database == "p" {
            ["phr", "pin", "psq"]
        } else {
            ["nhr", "nin", "nsq"]
        };
        for format in ["5", "13", "14", "15", "16"] {
            let case = format!("{program}_{format}");
            for threads in ["1", "4"] {
                let dir = unique_dir("losat_xml_dbname");
                fs::copy(
                    engine().join("tests/db_search/xb").join(query),
                    dir.join("q.fa"),
                )
                .expect("copy the query");
                for extension in extensions {
                    fs::copy(
                        engine().join(format!("tests/db_search/xb_{database}.{extension}")),
                        dir.join(format!("xb_caf\u{e9}.{extension}")),
                    )
                    .expect("copy the database");
                }
                let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
                for (name, _) in std::env::vars_os() {
                    let name = name.to_string_lossy().into_owned();
                    if name.starts_with("LOSAT_") || name.starts_with("NCBI") || name == "BLASTDB" {
                        command.env_remove(name);
                    }
                }
                let files = home_and_worker_log(&mut command, "losat_xml_dbname");
                command
                    .current_dir(&*dir)
                    .args([program, "-query", "q.fa", "-db", "xb_caf\u{e9}"])
                    .args(["-num_threads", threads, "-outfmt", format]);
                if format == "13" || format == "14" {
                    command.arg("-out").arg(format!("{case}.out"));
                }
                let output = command.output().expect("run LOSAT");
                drop(files);
                let label = format!("{case} threads {threads}");
                let rc = String::from_utf8(frozen(&format!("{case}.rc"))).unwrap();
                assert_eq!(
                    output.status.code(),
                    Some(rc.trim().parse::<i32>().unwrap()),
                    "{label}: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&frozen(&format!("{case}.stdout"))),
                    "{label}: stdout"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    String::from_utf8_lossy(&frozen(&format!("{case}.stderr"))),
                    "{label}: stderr"
                );
                if format == "13" || format == "14" {
                    let mut expected = frozen_files(&case);
                    expected.extend(
                        names(&dir)
                            .into_iter()
                            .filter(|name| !name.starts_with(&case))
                            .collect::<Vec<_>>(),
                    );
                    expected.sort();
                    assert_eq!(names(&dir), expected, "{label}: the files");
                    for name in frozen_files(&case) {
                        assert_eq!(
                            fs::read(dir.join(&name)).unwrap(),
                            frozen(&name),
                            "{label}: {name}"
                        );
                    }
                }
            }
        }
    }
}
