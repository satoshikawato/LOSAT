//! ABI v1's BLASTN search (plan TD-1): the checks of the request stay in their frozen order,
//! NCBI's reader reads the inputs (`v1_input`), and the records enter the engine's `run_local`.
//! Moved from `algorithm/blastn/blast_engine/run.rs` in session SFc, S10.

use anyhow::Result;

use super::v1_input::{read_queries, read_subjects};
use crate::algorithm::blastn::blast_engine::{check_subjects_not_empty, run_local};
use crate::algorithm::blastn::hsp::BlastnOutputFormat;
use crate::algorithm::blastn::scoring::{check_losat_limits, check_scoring_options};
use crate::algorithm::blastn::BlastnArgs;
use crate::api::local_blast::{OutputSink, ReportOutputs};

/// ABI v1's reading of a BLASTN `-outfmt` value, which is frozen (plan TD-1): it keeps the
/// values that it accepted before BLASTN's reading of `-outfmt` followed NCBI (now
/// `blastinput/app.rs` `parse_formatting_string`).
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
        // `v1_output_format` reads 0, 6 and 7 only (no field list).
        BlastnOutputFormat::Seqalign(_)
        | BlastnOutputFormat::CustomTabular(_)
        | BlastnOutputFormat::Xml
        | BlastnOutputFormat::Xml2(_)
        | BlastnOutputFormat::Archive
        | BlastnOutputFormat::Sam(_) => {
            anyhow::bail!("unsupported BLASTN output format")
        }
    };
    // ABI v1 is frozen (plan TD-1) except fail-fast fixes. Before S07+ it checked the
    // thread count first (in its search pool) and gave the empty report of an empty
    // query, as NCBI does after reading the subjects and checking the options; S07+ keeps
    // NCBI's errors there (fail-fast) and adds its checks of LOSAT's limits only for a
    // search.
    crate::utils::threading::validate_threads(args.num_threads)?;
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
    let subjects = read_subjects(subject_fasta.as_bytes(), "BLASTN", false)?;
    let queries = read_queries(query_fasta.as_bytes(), "BLASTN", false)?;
    if queries.is_empty() {
        check_subjects_not_empty(&subjects)?;
        check_scoring_options(&args)?;
        return Ok(output);
    }
    if !subjects.is_empty() {
        check_scoring_options(&args)?;
        check_losat_limits(&args)?;
    }
    let mut stderr = std::io::stderr();
    run_local(
        args,
        &queries,
        &subjects,
        &mut ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr),
    )?;
    Ok(output)
}
