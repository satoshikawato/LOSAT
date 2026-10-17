//! `LOSAT_X_RPTSTREAM` (`src/algorithm/blastn/blast_engine/x_rpt_stream.rs`) through the library
//! entry `run_local` with the Web adapter's outputs: hit records and a format observer. Their HSP
//! numbers count over the whole search and the hit records are handed over once, so the switch
//! does not act there (`off=outputs`) and the reports are written after the loop.
//!
//! The property under test: with the switch set for this process (it is read once, at the first
//! search) and query batches of 4,000 letters (`BATCH_SIZE`), the formats, the warnings, the hit
//! records and the observer's row and section ranges are those of a search of one batch (every
//! query in the first batch). Written batch by batch they would not be: the observer's HSP
//! numbers would start again at every batch. This file holds one test, so that no other test of
//! its process reads the environment while it changes `BATCH_SIZE`.

#[allow(dead_code)]
mod run_local_support;

use run_local_support::{fixture_sequence, reader_records, run_formats, Run, TempFasta};
use LOSAT::algorithm::blastn::BlastnArgs;
use LOSAT::api::local_blast::run_local_blastn;
use LOSAT::blastinput::fasta_reader::ReaderConfig;
use LOSAT::cli::{try_parse_from, Cli, Commands};

fn args(query: &TempFasta, subject: &TempFasta) -> BlastnArgs {
    let (query, subject) = (
        query.0.display().to_string(),
        subject.0.display().to_string(),
    );
    let cli: Cli = try_parse_from([
        "LOSAT", "blastn", "-query", &query, "-subject", &subject, "-outfmt", "6",
    ])
    .expect("blastn argv");
    match cli.command {
        Commands::Blastn(args) => args,
        _ => unreachable!("blastn subcommand"),
    }
}

fn search(
    query: &TempFasta,
    subject: &TempFasta,
    formats: &[&str],
    batch_size: Option<&str>,
) -> Run {
    match batch_size {
        Some(size) => std::env::set_var("BATCH_SIZE", size),
        None => std::env::remove_var("BATCH_SIZE"),
    }
    let queries = reader_records(&query.0, ReaderConfig::query("BLASTN", false, true));
    let subjects = reader_records(&subject.0, ReaderConfig::subject("BLASTN", false, true));
    run_formats(formats, true, |outputs| {
        run_local_blastn(args(query, subject), &queries, &subjects, outputs)
    })
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
        assert!(!whole.hits.is_empty() && whole.ranges.iter().all(|r| !r.is_empty()));
        assert_eq!(batched.outputs, whole.outputs, "{formats:?}: formats");
        assert_eq!(
            batched.diagnostics, whole.diagnostics,
            "{formats:?}: warnings"
        );
        assert_eq!(
            format!("{:?}", batched.hits),
            format!("{:?}", whole.hits),
            "{formats:?}: hit records"
        );
        assert_eq!(batched.ranges, whole.ranges, "{formats:?}: observer ranges");
        assert_eq!(
            batched.subjects, whole.subjects,
            "{formats:?}: observer subject headings"
        );
        // The HSP numbers run over the whole search once in every format.
        for ranges in &batched.ranges {
            let indices: Vec<usize> = ranges.iter().map(|range| range.0).collect();
            assert_eq!(
                indices,
                (0..batched.hits.len()).collect::<Vec<_>>(),
                "{formats:?}"
            );
        }
    }
}
