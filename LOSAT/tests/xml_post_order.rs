//! Where NCBI BLAST+ 2.17.0's XML streams post what they replace or throw, in a merged stream
//! (`2>&1`, a terminal): item K4h, review 1 M3 and M2. NCBI's diagnostics go to `cerr`, which
//! flushes `cout` first, so a post lands after the bytes that NCBI's output buffer had written:
//! the outfmt 5 document is written a `Hit` at a time, the XML2 report of outfmt 16 and a side
//! file of outfmt 14 in buffers of 4 KiB (`report::serial::out_buffer::xml_written_at`), so a
//! control byte's 803.7 post after many hits comes in the middle of the report
//! (`c01at30.*`, `ent1at40.f5`, the message's post of `msg.c01.f16`), an exception's post after
//! the bytes before the failing one (`u2at60.f16`), a post and then a throw in one report
//! (`two.c01u2.*`), and an `Iteration` left out after earlier ones (`q3u2.f5`). With the
//! prefetch off (`PRE_FETCH_SEQS_LIMIT`), the `Iteration` left out loads none of its subjects: the
//! texts of the taxonomy database's lookups (`m2*`: `x16_taxid_99990`, whose two sequences share
//! taxid 99990, with `tests/seqdb/taxdb/offset_last_past_end`'s offset error or `bad_magic`'s text
//! of opening) come at NCBI's next load that makes the same lookup, of whichever sequence
//! (`report::query_warnings::ReadLoad`). Each case's expectation is NCBI 2.17.0's merged stream of
//! the same run, and for `-out` the files it wrote (`tests/fixtures/xml_post_order/<case>.*`;
//! inputs `tests/db_search/make_xpost_inputs.py`).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:548-554
//! ```c++
//! inline bool BAD_CHAR(char x) {
//!     return (x < 0x20 && x > 0x0 && x != 0x9 && x != 0xA && x != 0xD);
//! }
//! inline char CObjectOStreamXml::x_VerifyChar(char x) {
//!     return BAD_CHAR(x) ?
//!         ReplaceVisibleChar(x, x_FixCharsMethod(), this, kEmptyStr, x_FixCharsSubst()) : x;
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375
//! ```c++
//!     const CSeq_id& kSeqId = kAlignSet.Get().front()->GetSeq_id(1);
//!
//!     try {
//!         const CBioseq_Handle& kSubjBioseqHandle = scope->GetBioseqHandle(kSeqId);
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// A case: name, program, options (a `.faa` name is an input under `tests/db_search/xpost`; with
/// a taxonomy database the paths are relative to `tests/db_search`), the taxonomy database under
/// `tests/seqdb/taxdb` and `PRE_FETCH_SEQS_LIMIT` (`None`: a `-subject` run, in an empty
/// directory), whether the run writes `-out out`, and NCBI's exit status.
type Case = (
    &'static str,
    &'static str,
    &'static [&'static str],
    Option<(&'static str, &'static str)>,
    bool,
    i32,
);

const CASES: &[Case] = &[
    (
        "c01at30.f16",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s60_c01_at30.faa",
            "-outfmt",
            "16",
        ],
        None,
        false,
        0,
    ),
    (
        "c01at30.f5",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s60_c01_at30.faa",
            "-outfmt",
            "5",
        ],
        None,
        false,
        0,
    ),
    (
        "c01at30.f14",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s60_c01_at30.faa",
            "-outfmt",
            "14",
        ],
        None,
        true,
        0,
    ),
    (
        "msg.c01.f16",
        "blastp",
        &[
            "-query",
            "q1_c01first_O.faa",
            "-subject",
            "s60.faa",
            "-outfmt",
            "16",
        ],
        None,
        false,
        0,
    ),
    (
        "ent1at40.f5",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s60_ent1_at40.faa",
            "-outfmt",
            "5",
        ],
        None,
        false,
        0,
    ),
    (
        "u2at60.f16",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s60_u2_at60.faa",
            "-outfmt",
            "16",
        ],
        None,
        false,
        0,
    ),
    (
        "two.c01u2.f16",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s2_c01_u2.faa",
            "-outfmt",
            "16",
        ],
        None,
        false,
        0,
    ),
    (
        "two.c01u2.f5",
        "blastp",
        &[
            "-query",
            "q1.faa",
            "-subject",
            "s2_c01_u2.faa",
            "-outfmt",
            "5",
        ],
        None,
        false,
        0,
    ),
    (
        "q3u2.f5",
        "blastp",
        &[
            "-query",
            "q3_u2mid.faa",
            "-subject",
            "s60.faa",
            "-outfmt",
            "5",
        ],
        None,
        false,
        0,
    ),
    (
        "m2qa.off.f5",
        "blastp",
        &[
            "-query",
            "xpost/m2_qa.faa",
            "-db",
            "x16_taxid_99990",
            "-outfmt",
            "5",
        ],
        Some(("offset_last_past_end", "1")),
        false,
        0,
    ),
    (
        "m2qa.magic.f5",
        "blastp",
        &[
            "-query",
            "xpost/m2_qa.faa",
            "-db",
            "x16_taxid_99990",
            "-outfmt",
            "5",
        ],
        Some(("bad_magic", "0")),
        false,
        0,
    ),
    (
        "m2qd.off.f5",
        "blastp",
        &[
            "-query",
            "xpost/m2_qd.faa",
            "-db",
            "x16_taxid_99990",
            "-outfmt",
            "5",
        ],
        Some(("offset_last_past_end", "0")),
        false,
        0,
    ),
];

fn tests_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
}

/// Runs a case with standard output and standard error into one file, in a fresh directory of the
/// target's temporary space: its exit status, the merged stream's bytes and, for `-out`, the
/// files the run wrote (each `=== <name> <length>` and its bytes, by name).
fn run(case: &Case) -> (i32, Vec<u8>, Vec<u8>) {
    let (name, program, options, taxdb, to_file, _) = *case;
    let dir = std::env::temp_dir().join(format!(
        "losat-xml-post-order-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("merged");
    let file = std::fs::File::create(&path).unwrap();
    let db_search = tests_dir().join("db_search");
    let mut command = Command::new(BIN);
    command.arg(program);
    match taxdb {
        None => {
            command
                .current_dir(&dir)
                .args(options.iter().map(|&option| {
                    if option.ends_with(".faa") {
                        db_search.join("xpost").join(option).into_os_string()
                    } else {
                        option.into()
                    }
                }));
        }
        Some(_) => {
            command.current_dir(&db_search).args(options);
        }
    }
    if to_file {
        command.args(["-out", "out"]);
    }
    command
        .stdout(Stdio::from(file.try_clone().unwrap()))
        .stderr(Stdio::from(file));
    for (variable, _) in std::env::vars() {
        // `LOSAT_BLASTX_WORKER_LOG` stays: an `--all-features` binary writes its probe there.
        if variable.starts_with("LOSAT_X_")
            || variable.starts_with("NCBI")
            || matches!(
                variable.as_str(),
                "BATCH_SIZE" | "BLASTDB" | "CHUNK_SIZE" | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(variable);
        }
    }
    if let Some((taxdb, limit)) = taxdb {
        command
            .env("BLASTDB", format!("../seqdb/taxdb/{taxdb}"))
            .env("PRE_FETCH_SEQS_LIMIT", limit);
    }
    command.env("HOME", &dir).env("USERPROFILE", &dir);
    let status = command.status().expect("run LOSAT CLI");
    let merged = std::fs::read(&path).unwrap();
    let mut files = Vec::new();
    if to_file {
        let mut names: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .filter(|name| name != "merged")
            .collect();
        names.sort();
        for name in names {
            let data = std::fs::read(dir.join(&name)).unwrap();
            files.extend_from_slice(
                format!("=== {} {}\n", name.to_string_lossy(), data.len()).as_bytes(),
            );
            files.extend(data);
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
    (status.code().unwrap_or(-1), merged, files)
}

#[test]
fn the_xml_streams_post_where_ncbi_posts_in_a_merged_stream() {
    let fixtures = tests_dir().join("fixtures").join("xml_post_order");
    for case in CASES {
        let (name, .., to_file, rc) = *case;
        let (status, merged, files) = run(case);
        let expected = std::fs::read(fixtures.join(format!("{name}.merged"))).unwrap();
        assert!(
            status == rc && merged == expected,
            "{name}: exit {status} (NCBI {rc}):\n{}",
            String::from_utf8_lossy(&merged)
        );
        if to_file {
            let expected = std::fs::read(fixtures.join(format!("{name}.files"))).unwrap();
            assert!(
                files == expected,
                "{name}: the files:\n{}",
                String::from_utf8_lossy(&files)
            );
        }
    }
}
