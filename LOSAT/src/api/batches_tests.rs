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
