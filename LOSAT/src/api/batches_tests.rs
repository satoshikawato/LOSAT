use super::*;
use crate::api::{run_collect, validate, FastaRecord, Program};

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

#[test]
fn a_batch_ends_where_it_reaches_its_size_and_a_query_past_the_interval_adds_nothing() {
    let query: Vec<FastaRecord> = (0..4)
        .map(|index| FastaRecord::new(format!("Query_{}", index + 1), b"q", PROTEIN))
        .collect();
    let subject = vec![FastaRecord::new("Subject_1", b"s", PROTEIN)];
    let search = validate(Program::Blastp, &[]).unwrap();
    let batches = run_batches(&search, Inputs::new(&query, &subject)).unwrap();
    assert_eq!(batches.query_count(), 4);
    let mut batches = batches;
    // The default size is 10,000 letters: one batch.
    let only = batches.next().unwrap().unwrap();
    assert_eq!(only.query_indices, 0..4);
    assert_eq!(only.hsps.len(), 4);
    assert!(batches.next().is_none());
    assert!(batches.next().is_none());
    let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
    assert_eq!(only.hsps, outcome.hsps);
    assert_eq!(only.query_stats, outcome.query_stats);
}

/// The version 4 protein database `tests/db_search/t5/phd_utf8` (its OID 1 a UTF8String
/// title) copied to `dir` as `name`, with the header record of OID 2 cut to its first byte:
/// NCBI's searches end there with exit status 0 and no message (`run::RunEnd::WithoutPost`).
#[cfg(unix)]
fn utf8_then_cut_database(dir: &std::path::Path) -> String {
    let source = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/t5");
    let read = |ext: &str| std::fs::read(source.join(format!("phd_utf8.{ext}"))).unwrap();
    let mut index = read("pin");
    let header = read("phr");
    let word = |at: usize| u32::from_be_bytes(index[at..at + 4].try_into().unwrap()) as usize;
    let mut at = 8;
    for _ in 0..2 {
        at += 4 + word(at);
    }
    let count = word(at);
    at += 4 + 8 + 4;
    let offsets: Vec<usize> = (0..=count).map(|k| word(at + 4 * k)).collect();
    let mut records: Vec<Vec<u8>> = (0..count)
        .map(|k| header[offsets[k]..offsets[k + 1]].to_vec())
        .collect();
    records[2] = vec![0x30];
    let mut offset = 0u32;
    for (k, record) in std::iter::once(&Vec::new()).chain(&records).enumerate() {
        offset += record.len() as u32;
        index[at + 4 * k..at + 4 * k + 4].copy_from_slice(&offset.to_be_bytes());
    }
    std::fs::write(dir.join("u.pin"), index).unwrap();
    std::fs::write(dir.join("u.phr"), records.concat()).unwrap();
    std::fs::write(dir.join("u.psq"), read("psq")).unwrap();
    dir.join("u").to_str().unwrap().to_string()
}

/// The batch in which the run ends without a post and with exit status 0 is the last, and
/// carries its diagnostics: NCBI's warning of the UTF8String read before the end, which the
/// command line writes (review 2 of RG-T5d, D3), in a run of one batch (with the run's own
/// messages) and in the second of two.
#[cfg(unix)]
#[test]
fn the_batch_where_the_run_ends_carries_its_diagnostics() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-batches-end-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let db = utf8_then_cut_database(&dir);
    let fasta = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/t5/hd_q.faa"),
    )
    .unwrap();
    let letters: Vec<&str> = fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    let record = |name: &str, seq: &str| FastaRecord::new(name.to_string(), b"", seq.as_bytes());
    let warning = "Warning: (803.10) [blastp] CObjectIStreamAsnBinary: UTF8String data for \
                   VisibleString member Blast-def-line-set.E.title, ASN.1 specification may \
                   need an update\n";
    let search = validate(Program::Blastp, &["-evalue", "1e-10"]).unwrap();
    // One batch: q0 reports OID 0, q1 reads OID 1 (the warning), q2 OID 2 (the end).
    let one: Vec<FastaRecord> = (0..3)
        .map(|k| record(&format!("q{k}"), letters[k]))
        .collect();
    let mut batches = run_batches(&search, Inputs::database(&one, &db)).unwrap();
    let only = batches.next().unwrap().unwrap();
    assert_eq!(only.query_indices, 0..3);
    assert!(only.hsps.is_empty());
    assert_eq!(only.query_stats, vec![None; 3]);
    assert_eq!(String::from_utf8_lossy(&only.diagnostics), warning);
    assert!(batches.next().is_none());
    // Two batches: a first query of 10,000 letters fills the first; the second ends the run.
    let long: String = letters[0].chars().cycle().take(10_000).collect();
    let two = vec![
        record("q0", &long),
        record("q1", letters[1]),
        record("q2", letters[2]),
    ];
    let mut batches = run_batches(&search, Inputs::database(&two, &db)).unwrap();
    let first = batches.next().unwrap().unwrap();
    assert_eq!(first.query_indices, 0..1);
    assert!(first.diagnostics.is_empty());
    let last = batches.next().unwrap().unwrap();
    assert_eq!(last.query_indices, 1..3);
    assert!(last.hsps.is_empty());
    assert_eq!(String::from_utf8_lossy(&last.diagnostics), warning);
    assert!(batches.next().is_none());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_failure_ends_the_batches() {
    let query = vec![FastaRecord::new("Query_1", b"q", PROTEIN)];
    let search = validate(Program::Blastp, &[]).unwrap();
    let mut batches = run_batches(&search, Inputs::new(&query, &[])).unwrap();
    assert!(matches!(batches.next(), Some(Err(Error::Search(_)))));
    assert!(batches.next().is_none());
}

#[test]
fn queries_without_a_record_are_one_call() {
    let subject = vec![FastaRecord::new("Subject_1", b"s", PROTEIN)];
    let search = validate(Program::Blastp, &[]).unwrap();
    let mut batches = run_batches(&search, Inputs::new(&[], &subject)).unwrap();
    let only = batches.next().unwrap().unwrap();
    assert_eq!(only.query_indices, 0..0);
    assert_eq!(
        String::from_utf8_lossy(&only.diagnostics),
        "Warning: [blastp] Query is Empty!\n"
    );
    assert!(batches.next().is_none());
}
