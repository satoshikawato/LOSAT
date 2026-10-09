//! ABI v1's BLASTN search (plan TD-1): `bio` reads the inputs, ABI v1's checks of them come
//! in their frozen order, and the records that pass them enter the engine's `run_local` as
//! the records of NCBI's reader (`v1_bio::from_bio`). Moved unchanged from
//! `algorithm/blastn/blast_engine/run.rs` in session SFc, S10.

use anyhow::{Context, Result};

use super::v1_bio::from_bio;
use crate::algorithm::blastn::blast_engine::{check_subjects_not_empty, run_local};
use crate::algorithm::blastn::hsp::BlastnOutputFormat;
use crate::algorithm::blastn::scoring::{check_losat_limits, check_scoring_options};
use crate::algorithm::blastn::BlastnArgs;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;

fn fasta_records_from_bytes(bytes: &[u8]) -> Result<Vec<bio::io::fasta::Record>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:484-490
    // ```c
    // void
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  objects::ENa_strand strand_opt,
    //                  TSearchMessages& messages)
    // ```
    bio::io::fasta::Reader::new(bytes)
        .records()
        // NCBI reference: c++/src/objtools/readers/fasta.cpp:428-431
        // FASTA_ERROR(LineNumber(), "CFastaReader: Expected defline around line " << LineNumber(), ...);
        .collect::<std::result::Result<Vec<_>, _>>()
        .context("failed to parse in-memory FASTA")
}

/// ABI v1's reading of a BLASTN `-outfmt` value, which is frozen (plan TD-1): it keeps the
/// values that it accepted before `parse_blastn_output_format` followed NCBI.
fn v1_output_format(spec: &str) -> Result<BlastnOutputFormat> {
    let mut parts = spec.split_whitespace();
    let format = parts.next().unwrap_or("6");
    if parts.next().is_some() {
        anyhow::bail!("unsupported BLASTN custom outfmt specification: {spec:?}");
    }
    match format {
        "0" => Ok(BlastnOutputFormat::Pairwise),
        "6" => Ok(BlastnOutputFormat::Tabular),
        "7" => Ok(BlastnOutputFormat::TabularWithComments),
        _ => anyhow::bail!("unsupported BLASTN output format: {format}"),
    }
}

pub(super) fn run_web_pair(
    args: BlastnArgs,
    query_fasta: &str,
    subject_fasta: &str,
) -> Result<Vec<u8>> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // ```c
    // BlastSeqBlkSetSequence(subj, sequence.data.release(),
    //    ((sentinels == eSentinels) ? sequence.length - 2 :
    //     sequence.length));
    // ...
    // SBlastSequence compressed_seq =
    //     subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
    //                               eNa_strand_plus, eNoSentinels);
    // ```
    let queries = fasta_records_from_bytes(query_fasta.as_bytes()).context("query FASTA")?;
    let subjects = fasta_records_from_bytes(subject_fasta.as_bytes()).context("subject FASTA")?;
    // NCBI reference: ncbi-blast/c++/include/algo/blast/api/local_blast.hpp:76-78
    // ```c
    // CLocalBlast(CRef<IQueryFactory> query_factory,
    //             CRef<CBlastOptionsHandle> opts_handle,
    //             CRef<CLocalDbAdapter> db);
    // ```
    // Web ABI v1 runs the same local search as the CLI and keeps the report in memory.
    let mut output = Vec::new();
    let outfmt = args.outfmt.clone();
    // Plan TD-1: ABI v1 is frozen except for fail-fast fixes, so its BLASTN keeps the
    // formats that it had (6 and 7) and rejects outfmt 0 with the error that it gave
    // before the engine implemented outfmt 0.
    let outfmt = match v1_output_format(&outfmt)? {
        BlastnOutputFormat::Pairwise => anyhow::bail!("unsupported BLASTN output format: 0"),
        BlastnOutputFormat::Tabular => "6".to_string(),
        BlastnOutputFormat::TabularWithComments => "7".to_string(),
    };
    // ABI v1 is frozen (plan TD-1) except fail-fast fixes. Before S07+ it checked the
    // thread count first (in its search pool) and gave the empty report of an empty
    // query, as NCBI does after reading the subjects and checking the options; S07+ keeps
    // NCBI's errors there (fail-fast) and adds its checks of the deflines, the records and
    // LOSAT's limits only for a search.
    crate::utils::threading::validate_threads(args.num_threads)?;
    if queries.is_empty() {
        check_subjects_not_empty(&subjects)?;
        check_scoring_options(&args)?;
        return Ok(output);
    }
    use super::v1_bio::{check_deflines, check_records, check_sequence_lines};
    // The deflines that NCBI reads differently are rejected (a fail-fast fix, plan TD-1).
    check_deflines(subject_fasta.as_bytes(), "subject")?;
    check_deflines(query_fasta.as_bytes(), "query")?;
    check_sequence_lines(subject_fasta.as_bytes(), "subject")?;
    check_sequence_lines(query_fasta.as_bytes(), "query")?;
    if !subjects.is_empty() {
        check_scoring_options(&args)?;
        check_losat_limits(&args)?;
        check_records(&subjects, "subject")?;
        check_records(&queries, "query")?;
    }
    // ABI v1 keeps `bio` and its checks (plan TD-1, `AUTHORITY.md` §J-5): the records that
    // pass them are those that NCBI's reader reads alike, and enter the search as its
    // records (`FastaRecord::from_bio`, which reads `U` as `T`).
    let queries: Vec<FastaRecord> = queries
        .iter()
        .enumerate()
        .map(|(index, record)| from_bio(record, index + 1, "Query_", false))
        .collect();
    let subjects: Vec<FastaRecord> = subjects
        .iter()
        .enumerate()
        .map(|(index, record)| from_bio(record, index + 1, "Subject_", false))
        .collect();
    let mut stderr = std::io::stderr();
    run_local(
        args,
        &queries,
        &subjects,
        &mut ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr),
    )?;
    Ok(output)
}
