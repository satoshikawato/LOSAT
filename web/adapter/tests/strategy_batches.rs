//! LOSAT Web's search strategy file of a run of more than one query batch (line RG item F1w,
//! review 1 S3): the rows of `LOSAT/tests/fixtures/strategy/cases.tsv` whose `BATCH_SIZE` gives
//! their queries two or more batches, with `-subject`, with `-subject` and `-searchsp` (each
//! request holds the search space that the batch before it set), and with `-db` and `-outfmt 11`
//! (the archive and the strategy share one options handle). The file of each run (stream 8) holds
//! NCBI BLAST+ 2.17.0's bytes in LOSAT's words, and the argv's output NCBI's standard output.
//!
//! The browser gives the engine no environment, so its runs batch at NCBI's default sizes: a test
//! process of its own sets `BATCH_SIZE` before each run, one run at a time, as the command line's
//! rows set it, which reaches the same batch loops with the rows' small inputs and NCBI's files.
//! Its working directory is `LOSAT/`, where the rows' names are.

use losat_web_adapter::{emit, run, store, zip};
use std::collections::BTreeMap;

/// A row of `cases.tsv`: its program, its `BATCH_SIZE` and its arguments, with the strategy file
/// named `s.asn`.
struct Row {
    program: String,
    batch_size: String,
    words: Vec<String>,
}

fn row(name: &str) -> Row {
    let cases = std::fs::read_to_string("tests/fixtures/strategy/cases.tsv").unwrap();
    let cells: Vec<&str> = cases
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|cells| cells[0] == name)
        .unwrap_or_else(|| panic!("no row {name}"));
    assert_eq!((cells[4], cells[5]), ("ncbi", "0"), "{name}");
    let batch_size = cells[2]
        .strip_prefix("BATCH_SIZE=")
        .unwrap_or_else(|| panic!("{name}: {}", cells[2]));
    Row {
        program: cells[1].to_string(),
        batch_size: batch_size.to_string(),
        words: cells[3]
            .split(' ')
            .map(|word| word.replace("{ss}", "s.asn"))
            .collect(),
    }
}

/// NCBI's file `name.ext` of the fixtures, in LOSAT's words (empty when NCBI wrote none).
fn frozen(name: &str, ext: &str) -> Vec<u8> {
    let bytes = std::fs::read(format!("tests/fixtures/strategy/{name}.{ext}")).unwrap_or_default();
    losat::web_adapter::losat_output_bytes(&bytes).into_owned()
}

/// The value of `option` in `words`.
fn value<'a>(words: &'a [String], option: &str) -> Option<&'a str> {
    let at = words.iter().position(|word| word == option)?;
    words.get(at + 1).map(String::as_str)
}

/// The run of the row in LOSAT Web (`run::run`), with its query and subject registered as the
/// host registers them: its result and the bytes of each stream that it emitted, joined.
fn web_run(row: &Row) -> (Result<(), String>, BTreeMap<u32, Vec<u8>>) {
    let mut words: Vec<&str> = vec![row.program.as_str()];
    words.extend(row.words.iter().map(String::as_str));
    let program = run::Program::parse(words[0]).unwrap();
    let loaders = run::data_loaders(&words).unwrap();
    let read = |option| std::fs::read(value(&row.words, option).unwrap()).unwrap();
    let (query, _) = store::register(words[0], 0, loaders, &read("-query")).unwrap();
    let subject = match run::searches_database(&words) {
        true => 0,
        false => {
            store::register(words[0], 1, loaders, &read("-subject"))
                .unwrap()
                .0
        }
    };
    emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
    let result = store::with_inputs(
        program,
        query,
        subject,
        &words,
        Some(loaders),
        |queries, end, subjects| run::run(&words, queries, end, subjects),
    );
    store::release(query).unwrap();
    if subject != 0 {
        store::release(subject).unwrap();
    }
    let mut streams = BTreeMap::new();
    emit::EMITTED.with(|emitted| {
        for (stream, bytes) in emitted.borrow().iter() {
            streams
                .entry(*stream)
                .or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        }
    });
    (result.unwrap(), streams)
}

#[test]
fn a_run_of_several_batches_writes_ncbis_strategy_file() {
    std::env::set_current_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT")).unwrap();
    for program in ["blastp", "blastn", "blastx", "tblastn", "tblastx"] {
        for kind in ["sub.batches", "sub.searchsp.batches", "db.batches.f11"] {
            let name = format!("{program}.{kind}");
            let row = row(&name);
            // One run at a time in this process: the variable is read where each run batches.
            std::env::set_var("BATCH_SIZE", &row.batch_size);
            let (result, streams) = web_run(&row);
            std::env::remove_var("BATCH_SIZE");
            result.unwrap_or_else(|error| panic!("{name}: {error}"));
            let files: Vec<(Vec<u8>, Vec<u8>)> = zip::read(&streams[&emit::STREAM_RUN_FILES]);
            let want = frozen(&name, "ss");
            assert!(
                files.len() == 1 && files[0].0 == b"s.asn" && files[0].1 == want,
                "{name}: the file of the run ({} bytes) is not NCBI's ({} bytes)",
                files.first().map_or(0, |(_, bytes)| bytes.len()),
                want.len()
            );
            // NCBI's requests are one per batch: more than one here.
            let requests = String::from_utf8_lossy(&want)
                .matches("Blast4-request ::=")
                .count();
            assert!(requests > 1, "{name}: {requests} request(s)");
            // The argv's own output: outfmt 6 on its stream, outfmt 11 (the archive) on the
            // stream of the other formats.
            let stream = match value(&row.words, "-outfmt") {
                Some("6") => 6,
                _ => emit::STREAM_CUSTOM,
            };
            assert!(
                streams.get(&stream).map(Vec::as_slice).unwrap_or_default()
                    == frozen(&name, "out").as_slice(),
                "{name}: the output on stream {stream} is not NCBI's standard output"
            );
            assert_eq!(
                String::from_utf8_lossy(
                    streams
                        .get(&emit::STREAM_DIAGNOSTICS)
                        .map(Vec::as_slice)
                        .unwrap_or_default()
                ),
                String::from_utf8_lossy(&frozen(&name, "err")),
                "{name}: stream 3 is NCBI's standard error"
            );
        }
    }
}
