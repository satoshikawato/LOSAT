//! psiblast's XML2 and JSON documents (outfmt 13 to 16, line PS item PS15b-4) in a run of LOSAT
//! Web: the format's stream holds NCBI BLAST+ 2.17.0's frozen bytes of the same argv in LOSAT's
//! words (the rows `db.ps15b4.*` of `LOSAT/tests/db_search_manifest/psiblast.ps15b4.tsv`, their
//! files under `LOSAT/tests/fixtures/db_search`), while the run's outfmt 0, 6 and 7, which a run
//! writes too in searches of their own (`search_groups`: psiblast's list cut of 13 to 16,
//! `psiblast_list_cut`), are those of the run without `-outfmt`. Review S2 of PS15b-4.

use super::*;
use std::collections::BTreeMap;

/// The repository's `LOSAT/` directory, where the paths of the rows start.
fn losat_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../LOSAT")
}

/// The file `path` of `LOSAT/`.
fn read(path: &str) -> Vec<u8> {
    std::fs::read(losat_dir().join(path)).unwrap_or_else(|error| panic!("{path}: {error}"))
}

/// NCBI's frozen file `name` of a row (`LOSAT/tests/fixtures/db_search`), in LOSAT's words.
fn frozen(name: &str) -> Vec<u8> {
    losat::web_adapter::losat_output_bytes(&read(&format!("tests/fixtures/db_search/{name}")))
        .into_owned()
}

/// The streams emitted since `EMITTED` was cleared, each stream's chunks joined.
fn emitted() -> BTreeMap<u32, Vec<u8>> {
    let mut streams = BTreeMap::new();
    crate::emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            streams
                .entry(*stream)
                .or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        }
    });
    streams
}

/// The result and streams of `run` of `words` with the FASTA query `query` and subject `subject`
/// registered as a host registers them.
fn subject_run(
    words: &[&str],
    query: &[u8],
    subject: &[u8],
) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let program = Program::parse(words[0]).unwrap();
    let (query, _) = crate::store::register(words[0], 0, true, false, query).unwrap();
    let (subject, _) = crate::store::register(words[0], 1, true, false, subject).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = crate::store::with_inputs(
        program,
        query,
        subject,
        words,
        Some(true),
        |queries, end, subjects| run(words, queries, end, subjects),
    );
    crate::store::release(query).unwrap();
    crate::store::release(subject).unwrap();
    (result.unwrap(), emitted())
}

/// The result and streams of `run` of `words` (`-db`), whose paths starting with `tests/` name
/// files of `LOSAT/` (the run's working directory, as the host's preopened `.` gives it), with the
/// query file `query` of `LOSAT/` and no subject (line WQA D12a).
fn database_run(words: &[&str], query: &str) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let dir = losat_dir();
    let words: Vec<String> = words
        .iter()
        .map(|word| {
            if word.starts_with("tests/") {
                dir.join(word).to_string_lossy().into_owned()
            } else {
                word.to_string()
            }
        })
        .collect();
    let words: Vec<&str> = words.iter().map(String::as_str).collect();
    let loaders = data_loaders(&words).unwrap();
    let (handle, _) = crate::store::register(words[0], 0, loaders, false, &read(query)).unwrap();
    crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let program = Program::parse(words[0]).unwrap();
    let result = crate::store::with_inputs(
        program,
        handle,
        0,
        &words,
        Some(loaders),
        |queries, end, subjects| {
            assert!(subjects.is_empty());
            run(&words, queries, end, subjects)
        },
    );
    crate::store::release(handle).unwrap();
    (result.unwrap(), emitted())
}

/// Bytes as text, so that a failed comparison prints readable text.
fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The streams of `with` (a run with `-outfmt`) that a run writes besides the format's own (0, 6,
/// 7 and the hits) equal those of `plain` (the run without it), and so does stream 3: NCBI's
/// standard error of the row is empty.
fn others_are_those_of_the_plain_run(
    row: &str,
    with: &BTreeMap<u32, Vec<u8>>,
    plain: &BTreeMap<u32, Vec<u8>>,
) {
    for stream in [0, 6, 7, STREAM_HITS, STREAM_DIAGNOSTICS] {
        assert!(
            with.get(&stream) == plain.get(&stream),
            "{row}: stream {stream}: {:?}",
            with.get(&stream).map(|bytes| text(bytes))
        );
    }
}

// Two FASTA queries, two rounds, `-subject` (rows `db.ps15b4.s.two.it2.15` and `.16`): one
// report per query holding a `Search` per round.
#[test]
fn psiblast_writes_ncbis_xml2_and_json_of_a_subject_search() {
    let query = read("tests/fasta/psiblast/two.faa");
    let subject = read("tests/fasta/PajaWSV.faa");
    let words = [
        "psiblast",
        "-query",
        "q",
        "-subject",
        "s",
        "-num_iterations",
        "2",
    ];
    let (plain_result, plain) = subject_run(&words, &query, &subject);
    plain_result.unwrap();
    for outfmt in ["15", "16"] {
        let row = format!("db.ps15b4.s.two.it2.{outfmt}");
        let words = [&words[..], &["-outfmt", outfmt]].concat();
        validate(&words).unwrap_or_else(|error| panic!("{row}: {error}"));
        let (result, with) = subject_run(&words, &query, &subject);
        result.unwrap_or_else(|error| panic!("{row}: {error}"));
        assert_eq!(
            text(with.get(&STREAM_CUSTOM).map_or(&[][..], Vec::as_slice)),
            text(&frozen(&format!("{row}.out"))),
            "{row}"
        );
        others_are_those_of_the_plain_run(&row, &with, &plain);
    }
}

// One FASTA query, three rounds, `-db tests/db_search/paja` (rows `db.ps15b4.f.tk.it3.13` to
// `.16`): the query's report is made after NCBI drops the query from its scope (its Seq-loc's id,
// no title, each Hsp without rows). Outfmt 13 and 14 are a ZIP archive on stream 4: NCBI's side
// files, named as the CLI names them for `-out results`, and NCBI's master file
// (`results.json` / `results.xml`), which lists them by those names.
#[test]
fn psiblast_writes_ncbis_xml2_and_json_of_a_database_search() {
    let query = "tests/fasta/psiblast/tk.faa";
    let words = [
        "psiblast",
        "-db",
        "tests/db_search/paja",
        "-query",
        "q",
        "-num_iterations",
        "3",
    ];
    let (plain_result, plain) = database_run(&words, query);
    plain_result.unwrap();
    for outfmt in ["13", "14", "15", "16"] {
        let row = format!("db.ps15b4.f.tk.it3.{outfmt}");
        let words = [&words[..], &["-outfmt", outfmt]].concat();
        validate(&words).unwrap_or_else(|error| panic!("{row}: {error}"));
        assert!(searches_database(&words), "{row}");
        let (result, with) = database_run(&words, query);
        result.unwrap_or_else(|error| panic!("{row}: {error}"));
        let stream = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
        let number: u32 = outfmt.parse().unwrap();
        if matches!(number, 15 | 16) {
            assert_eq!(text(&stream), text(&frozen(&format!("{row}.out"))), "{row}");
        } else {
            let extension = if number == 13 { "json" } else { "xml" };
            let mut expected = Vec::new();
            for n in 1.. {
                let side = format!("{row}.out_{n}.{extension}");
                if !losat_dir()
                    .join(format!("tests/fixtures/db_search/{side}"))
                    .exists()
                {
                    break;
                }
                expected.push((format!("results_{n}.{extension}"), frozen(&side)));
            }
            assert!(!expected.is_empty(), "{row}: side files");
            let master =
                text(&frozen(&format!("{row}.out"))).replace(&format!("{row}.out_"), "results_");
            expected.push((master_name(number), master.into_bytes()));
            let files = crate::zip::read(&stream);
            assert_eq!(files.len(), expected.len(), "{row}: entries");
            for ((name, bytes), (want_name, want)) in files.iter().zip(&expected) {
                assert_eq!(text(name), *want_name, "{row}");
                assert!(bytes == want, "{row}: {want_name}: {}", text(bytes));
            }
        }
        others_are_those_of_the_plain_run(&row, &with, &plain);
    }
}
