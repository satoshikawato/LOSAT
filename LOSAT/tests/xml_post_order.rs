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
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:372-375

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// A case: name, program, options (a `.faa` name is an input under `tests/db_search/xpost`; with
/// a taxonomy database the paths are relative to `tests/db_search`), the taxonomy database under
/// `tests/seqdb/taxdb` and the run's other variables (`VAR=VALUE`, separated by spaces:
/// `PRE_FETCH_SEQS_LIMIT`, `BATCH_SIZE`) (`None`: a `-subject` run, in an empty directory),
/// whether the run writes `-out out` (in the run's own directory), and NCBI's exit status.
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
        Some(("offset_last_past_end", "PRE_FETCH_SEQS_LIMIT=1")),
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
        Some(("bad_magic", "PRE_FETCH_SEQS_LIMIT=0")),
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
        Some(("offset_last_past_end", "PRE_FETCH_SEQS_LIMIT=0")),
        false,
        0,
    ),
    (
        "r2m1.bp.mid.off.pre.f14",
        "blastp",
        &[
            "-query",
            "xpost/r2_mid.faa",
            "-db",
            "x16_taxid_99990",
            "-outfmt",
            "14",
        ],
        Some(("offset_last_past_end", "")),
        true,
        255,
    ),
    (
        "r2m1.tn.first.off.lim0.f14",
        "tblastn",
        &[
            "-query",
            "xpost/r2_first.faa",
            "-db",
            "xpost/db/nST/nST",
            "-outfmt",
            "14",
        ],
        Some(("offset_last_past_end", "PRE_FETCH_SEQS_LIMIT=0")),
        true,
        255,
    ),
    (
        "r2m2e.bp.mv.off.lim0.f5",
        "blastp",
        &[
            "-query",
            "xpost/r2_mv.faa",
            "-db",
            "xpost/db/pA/pA xpost/db/pB/pB",
            "-outfmt",
            "5",
            "-evalue",
            "1e-20",
        ],
        Some(("offset_last_past_end", "PRE_FETCH_SEQS_LIMIT=0")),
        false,
        0,
    ),
    (
        "r2m2.tn.qa.off.lim0.f5",
        "tblastn",
        &[
            "-query",
            "xpost/r2_qa.faa",
            "-db",
            "xpost/db/nST/nST",
            "-outfmt",
            "5",
        ],
        Some(("offset_last_past_end", "PRE_FETCH_SEQS_LIMIT=0")),
        false,
        0,
    ),
    (
        "r2ce.mq3.lim0.f16",
        "blastx",
        &[
            "-query",
            "xpost/r2_mq3.fna",
            "-db",
            "xpost/db/x2gtax/x2gtax",
            "-ungapped",
            "-window_size",
            "0",
            "-comp_based_stats",
            "0",
            "-seg",
            "no",
            "-word_size",
            "5",
            "-threshold",
            "2",
            "-evalue",
            "1e5",
            "-outfmt",
            "16",
        ],
        Some(("offset_last_past_end", "PRE_FETCH_SEQS_LIMIT=0")),
        false,
        0,
    ),
    (
        "m2b.bp.qm.off.mix.f5",
        "blastp",
        &[
            "-query",
            "xpost/m2b_bp_qm.faa",
            "-db",
            "x16_taxid_99990",
            "-outfmt",
            "5",
        ],
        Some((
            "offset_last_past_end",
            "PRE_FETCH_SEQS_LIMIT=3 BATCH_SIZE=1",
        )),
        false,
        0,
    ),
    (
        "m2b.bn.qm.off.mix.f5",
        "blastn",
        &[
            "-query",
            "xpost/m2b_bn_qm.fna",
            "-db",
            "x16_ntaxid_99990",
            "-outfmt",
            "5",
        ],
        Some((
            "offset_last_past_end",
            "PRE_FETCH_SEQS_LIMIT=3 BATCH_SIZE=1",
        )),
        false,
        0,
    ),
    (
        "m2b.tn.qm.off.mix.f5",
        "tblastn",
        &[
            "-query",
            "xpost/m2b_tn_qm.faa",
            "-db",
            "x16_ntaxid_99990",
            "-outfmt",
            "5",
        ],
        Some((
            "offset_last_past_end",
            "PRE_FETCH_SEQS_LIMIT=3 BATCH_SIZE=1",
        )),
        false,
        0,
    ),
    (
        "m2b.tx.qm.off.mix.f5",
        "tblastx",
        &[
            "-query",
            "xpost/m2b_tx_qm.fna",
            "-db",
            "x16_ntaxid_99990",
            "-outfmt",
            "5",
            "-evalue",
            "1e-20",
        ],
        Some((
            "offset_last_past_end",
            "PRE_FETCH_SEQS_LIMIT=1 BATCH_SIZE=1",
        )),
        false,
        0,
    ),
];

/// Removes the variables that change a run from `command`'s environment.
fn clear_env(command: &mut Command) {
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
}

fn tests_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
}

/// The path of the case's query file (`-query`), as `run` gives it to LOSAT.
fn query_path(case: &Case) -> Option<PathBuf> {
    let (_, _, options, taxdb, ..) = *case;
    let at = options.iter().position(|&option| option == "-query")?;
    let db_search = tests_dir().join("db_search");
    Some(match taxdb {
        None => db_search.join("xpost").join(options[at + 1]),
        Some(_) => db_search.join(options[at + 1]),
    })
}

/// Runs a case with standard output and standard error into one file, in a fresh directory of the
/// target's temporary space: its exit status, the merged stream's bytes and, for `-out`, the
/// files the run wrote (each `=== <name> <length>` and its bytes, by name). With `defline`, the
/// query is that file's text instead of the case's (`balanced_title`).
fn run(case: &Case, query: Option<&str>) -> (i32, Vec<u8>, Vec<u8>) {
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
    let balanced_query = query.map(|query| {
        // Beside the run's directory: the files of `-out` are read from it.
        let path = dir.with_extension("q.faa");
        std::fs::write(&path, query).unwrap();
        path
    });
    let options: Vec<std::ffi::OsString> = options
        .iter()
        .enumerate()
        .map(|(at, &option)| match &balanced_query {
            Some(query) if at > 0 && options[at - 1] == "-query" => query.clone().into_os_string(),
            _ => option.into(),
        })
        .collect();
    let mut command = Command::new(BIN);
    command.arg(program);
    match taxdb {
        None => {
            command.current_dir(&dir).args(options.iter().map(|option| {
                if option.to_string_lossy().ends_with(".faa") && Path::new(option).is_relative() {
                    db_search.join("xpost").join(option).into_os_string()
                } else {
                    option.clone()
                }
            }));
        }
        Some(_) => {
            command.current_dir(&db_search).args(&options);
        }
    }
    if to_file {
        // In the run's directory (a run with a taxonomy database runs in `tests/db_search`).
        command.arg("-out").arg(match taxdb {
            None => PathBuf::from("out"),
            Some(_) => dir.join("out"),
        });
    }
    command
        .stdout(Stdio::from(file.try_clone().unwrap()))
        .stderr(Stdio::from(file));
    clear_env(&mut command);
    if let Some((taxdb, variables)) = taxdb {
        command.env("BLASTDB", format!("../seqdb/taxdb/{taxdb}"));
        for variable in variables.split_whitespace() {
            let (name, value) = variable.split_once('=').expect("VAR=VALUE");
            command.env(name, value);
        }
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
    if let Some(query) = balanced_query {
        std::fs::remove_file(query).unwrap();
    }
    (status.code().unwrap_or(-1), merged, files)
}

/// NCBI's merged stream of the row `name` in LOSAT's words (NCBI's bytes through
/// `losat::wording`, AGENTS.md rule 8).
fn expected_merged(fixtures: &Path, name: &str) -> Vec<u8> {
    let frozen = std::fs::read(fixtures.join(format!("{name}.merged"))).unwrap();
    losat::wording::losat_output_bytes(&frozen).into_owned()
}

/// The merged run of `case` and what LOSAT writes then: where LOSAT's version adds bytes before
/// NCBI's first post (`Error: (` or `Warning: (`, NCBI's standard output before it), the query
/// with a shorter title (`balanced_title`); otherwise the case and [`expected_merged`].
fn merged_run(fixtures: &Path, case: &Case) -> (Option<String>, Vec<u8>) {
    let (name, _, options, ..) = *case;
    let frozen = std::fs::read(fixtures.join(format!("{name}.merged"))).unwrap();
    let first_post = (0..frozen.len()).find(|&at| {
        frozen[at..].starts_with(b"Error: (") || frozen[at..].starts_with(b"Warning: (")
    });
    let balanced = first_post.zip(query_path(case)).and_then(|(at, query)| {
        let query = std::fs::read_to_string(query).unwrap();
        let parse_deflines = options.contains(&"-parse_deflines");
        crate::balanced_title::balanced_title(&frozen, &frozen[..at], &query, parse_deflines)
    });
    match balanced {
        Some(balanced) => (Some(balanced.query), balanced.expected),
        None => (None, expected_merged(fixtures, name)),
    }
}

/// NCBI's files of the row `name` in LOSAT's words: the `=== <name> <length>` records of `run`,
/// each file's bytes through `losat::wording` (AGENTS.md rule 8) and its length that of the words.
fn expected_files(fixtures: &Path, name: &str) -> Vec<u8> {
    let frozen = std::fs::read(fixtures.join(format!("{name}.files"))).unwrap();
    let mut words = Vec::new();
    let mut rest = frozen.as_slice();
    while !rest.is_empty() {
        let line_end = rest.iter().position(|&byte| byte == b'\n').unwrap();
        let header = std::str::from_utf8(&rest[..line_end]).unwrap();
        let (file, length) = header
            .strip_prefix("=== ")
            .and_then(|header| header.rsplit_once(' '))
            .expect("=== <name> <length>");
        let length: usize = length.parse().unwrap();
        let data = &rest[line_end + 1..line_end + 1 + length];
        let data = losat::wording::losat_output_bytes(data);
        words.extend_from_slice(format!("=== {file} {}\n", data.len()).as_bytes());
        words.extend_from_slice(&data);
        rest = &rest[line_end + 1 + length..];
    }
    words
}

#[test]
fn the_xml_streams_post_where_ncbi_posts_in_a_merged_stream() {
    let fixtures = tests_dir().join("fixtures").join("xml_post_order");
    for case in CASES {
        let (name, .., to_file, rc) = *case;
        let (query, expected) = merged_run(&fixtures, case);
        let (status, merged, files) = run(case, query.as_deref());
        let first = merged
            .iter()
            .zip(&expected)
            .position(|(got, want)| got != want)
            .unwrap_or(merged.len().min(expected.len()));
        let around = |bytes: &[u8]| {
            String::from_utf8_lossy(&bytes[first.saturating_sub(80)..(first + 80).min(bytes.len())])
                .into_owned()
        };
        assert!(
            status == rc && merged == expected,
            "{name}: exit {status} (NCBI {rc}); the streams differ from byte {first} (balanced \
             title: {}):\n--- LOSAT\n{}\n--- NCBI in LOSAT's words\n{}",
            query
                .as_deref()
                .map_or("no", |query| query.lines().next().unwrap_or("")),
            around(&merged),
            around(&expected)
        );
        if to_file {
            let expected = expected_files(&fixtures, name);
            assert!(
                files == expected,
                "{name}: the files:\n{}",
                String::from_utf8_lossy(&files)
            );
        }
    }
}

/// Review 2 of K4h, M1: `-export_search_strategy` to a FIFO, which cannot be cut back, in a run
/// that query 2's 14 side file ends (`BATCH_SIZE=1`): NCBI BLAST+ 2.17.0 writes the requests of
/// the batches through query 2's (2 of 3), and so does LOSAT, which holds each request until the
/// run reaches its batch.
#[cfg(unix)]
#[test]
fn a_strategy_fifo_gets_the_requests_that_ncbi_writes() {
    let fixtures = tests_dir().join("fixtures").join("xml_post_order");
    let dir =
        std::env::temp_dir().join(format!("losat-xml-post-order-{}-fifo", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let fifo = dir.join("ss.fifo");
    assert!(Command::new("mkfifo")
        .arg(&fifo)
        .status()
        .expect("mkfifo")
        .success());
    let reader = {
        let fifo = fifo.clone();
        std::thread::spawn(move || std::fs::read(fifo).unwrap())
    };
    let path = dir.join("merged");
    let file = std::fs::File::create(&path).unwrap();
    let mut command = Command::new(BIN);
    command
        .arg("blastp")
        .current_dir(tests_dir().join("db_search"))
        .args([
            "-query",
            "xpost/r2_s3mid.faa",
            "-db",
            "xpost/db/pR3/pR3",
            "-outfmt",
            "14",
        ])
        .arg("-out")
        .arg(dir.join("out"))
        .arg("-export_search_strategy")
        .arg(&fifo)
        .stdout(Stdio::from(file.try_clone().unwrap()))
        .stderr(Stdio::from(file));
    clear_env(&mut command);
    command
        .env("BATCH_SIZE", "1")
        .env("HOME", &dir)
        .env("USERPROFILE", &dir);
    let status = command.status().expect("run LOSAT CLI");
    // A run that ends before it opens the FIFO leaves the reader waiting in its open: a writer
    // that opens the FIFO without blocking and closes it ends that wait (the reader reads nothing),
    // so the run's failure shows below instead of a hang. Where the reader is done, nothing has the
    // FIFO open for reading and this open fails (ENXIO), which changes nothing.
    {
        use std::os::unix::fs::OpenOptionsExt;
        let _ = std::fs::OpenOptions::new()
            .write(true)
            .custom_flags(rustix::fs::OFlags::NONBLOCK.bits() as i32)
            .open(&fifo);
    }
    let strategy = reader.join().unwrap();
    let merged = std::fs::read(&path).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    let name = "fifo.bp.s3mid.f14";
    assert_eq!(
        status.code(),
        Some(255),
        "{}",
        String::from_utf8_lossy(&merged)
    );
    assert!(
        merged == expected_merged(&fixtures, name),
        "{name}:\n{}",
        String::from_utf8_lossy(&merged)
    );
    assert!(
        strategy
            == losat::wording::losat_output_bytes(
                &std::fs::read(fixtures.join(format!("{name}.ss"))).unwrap()
            )
            .into_owned(),
        "{name}: the strategy:\n{}",
        String::from_utf8_lossy(&strategy)
    );
}
