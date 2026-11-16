//! `LOSAT_X_RPTSTREAM` (`src/algorithm/blastn/blast_engine/x_rpt_stream.rs`) through the library
//! API (`run_collect`, which gives the entry `run_local` the Web adapter's outputs: hit records and
//! a format observer). Their HSP numbers count over the whole search and the hit records are
//! handed over once, so the switch does not act there (`off=outputs`) and the reports are written
//! after the loop.
//!
//! The property under test: with the switch set for this process (it is read once, at the first
//! search) and query batches of 4,000 letters (`BATCH_SIZE`), the formats, the warnings, the hit
//! records and the row and section ranges of the reports are those of a search of one batch (every
//! query in the first batch). Written batch by batch they would not be: the HSP numbers would
//! start again at every batch. This file holds one test, so that no other test of its process
//! reads the environment while it changes `BATCH_SIZE`.

mod run_local_support;

use losat::api::{Outcome, Program};
use run_local_support::{fixture_sequence, Files, TempFasta};

fn search(
    query: &TempFasta,
    subject: &TempFasta,
    formats: &[&str],
    batch_size: Option<&str>,
) -> Outcome {
    match batch_size {
        Some(size) => std::env::set_var("BATCH_SIZE", size),
        None => std::env::remove_var("BATCH_SIZE"),
    }
    Files::new(Program::Blastn, &query.0, &subject.0).collect(formats, &[])
}

#[test]
fn hit_records_and_an_observer_keep_the_reports_after_the_loop() {
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            std::env::remove_var(name);
        }
    }
    std::env::set_var("LOSAT_X_RPTSTREAM", "1");
    let genome = fixture_sequence("LC738874.fasta");
    let pieces: Vec<(String, Vec<u8>)> = (0..12usize)
        .map(|i| {
            let start = i * 13_000;
            (
                format!("piece{i}"),
                genome[start..start + 1_500 + 200 * i].to_vec(),
            )
        })
        .collect();
    let records: Vec<(&str, &[u8])> = pieces
        .iter()
        .map(|(name, seq)| (name.as_str(), seq.as_slice()))
        .collect();
    let query = TempFasta::new("rptstream_api_query.fna", &records);
    let subject_letters = fixture_sequence("LC738868.fasta");
    let subject = TempFasta::new("rptstream_api_subject.fna", &[("s1", &subject_letters)]);
    // One batch (the whole input in the first batch), then batches of 4,000 letters; one format
    // (the hit records and the observer alone keep the switch off) and three.
    for formats in [&["6"][..], &["0"], &["0", "6", "7"]] {
        let whole = search(&query, &subject, formats, None);
        let batched = search(&query, &subject, formats, Some("4000"));
        std::env::remove_var("BATCH_SIZE");
        assert!(
            !whole.hsps.is_empty()
                && whole
                    .reports
                    .iter()
                    .all(|report| report.spans.iter().flatten().count() > 0)
        );
        for (batched, whole) in batched.reports.iter().zip(&whole.reports) {
            assert_eq!(batched.bytes, whole.bytes, "{formats:?}: formats");
            assert_eq!(
                batched.spans, whole.spans,
                "{formats:?}: row and section ranges"
            );
            assert_eq!(
                batched.subject_headings, whole.subject_headings,
                "{formats:?}: subject headings"
            );
        }
        assert_eq!(
            batched.diagnostics, whole.diagnostics,
            "{formats:?}: warnings"
        );
        assert_eq!(
            format!("{:?}", batched.hsps),
            format!("{:?}", whole.hsps),
            "{formats:?}: hit records"
        );
        // The HSP numbers run over the whole search once in every format: every HSP has its
        // place, and the places follow the numbers in the report.
        for report in &batched.reports {
            assert_eq!(report.spans.len(), batched.hsps.len(), "{formats:?}");
            let mut end = 0;
            for (index, span) in report.spans.iter().enumerate() {
                let span = span
                    .clone()
                    .unwrap_or_else(|| panic!("{formats:?}: HSP {index}"));
                assert!(span.start >= end, "{formats:?}: HSP {index} in write order");
                end = span.end;
            }
        }
    }
}
