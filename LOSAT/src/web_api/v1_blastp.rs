//! ABI v1's BLASTP search (plan TD-1): the checks of the request stay frozen, NCBI's reader
//! reads the inputs (`v1_input`), and the records enter the engine's `run_local`. Moved from
//! `algorithm/blastp/blast_engine.rs` in session SFc, S10.

use anyhow::{bail, Result};

use super::v1_input::{read_queries, read_subjects};
use crate::algorithm::blastp::blast_engine::{check_options, run_local};
use crate::algorithm::blastp::BlastpArgs;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::report::OutputFormat;

pub(super) fn run_web_pair(
    args: BlastpArgs,
    query_fasta: &str,
    subject_fasta: &str,
) -> Result<Vec<u8>> {
    run_web_pair_bytes(
        args,
        query_fasta.as_bytes(),
        subject_fasta.as_bytes(),
        "",
        "",
    )
}

/// The tabular fields that ABI v1's BLASTP wrote (before line CLI C2 moved BLASTP to the shared
/// tabular writer, which writes every one of NCBI's fields): its arguments are frozen (plan TD-1),
/// so a field outside them stays its error.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:47-198
/// ```c++
/// const SFormatSpec sc_FormatSpecifiers[kNumTabularOutputFormatSpecifiers] = {
///     SFormatSpec("qseqid",
///                 "Query Seq-id",
///                 eQuerySeqId),
/// ```
const V1_TABULAR_FIELDS: [&str; 30] = [
    "qseqid", "qacc", "qaccver", "qlen", "sseqid", "sacc", "saccver", "slen", "qstart", "qend",
    "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident",
    "mismatch", "positive", "gapopen", "gaps", "ppos", "qframe", "sframe", "frames", "btop",
    "stitle",
];

/// ABI v1's checks of a BLASTP request, whose arguments, formats and errors are frozen
/// (plan TD-1): the options are checked before the output format, and a tabular field that
/// ABI v1's BLASTP did not write (`V1_TABULAR_FIELDS`) is an error, as in v1 before
/// S08+ moved BLASTP to NCBI's application layer (which reads `-outfmt` first and ignores a
/// token that is not a field).
fn check_web_v1_request(args: &BlastpArgs) -> Result<()> {
    check_options(args)?;
    let (format, fields) = OutputFormat::parse(&args.outfmt).map_err(anyhow::Error::msg)?;
    if let Some(fields) = fields {
        if format == OutputFormat::Pairwise {
            bail!("blastp outfmt 0 does not accept custom field lists");
        }
        for token in fields.split_whitespace() {
            if !token.eq_ignore_ascii_case("std") && !V1_TABULAR_FIELDS.contains(&token) {
                bail!("unsupported blastp tabular field '{token}'");
            }
        }
    }
    Ok(())
}

pub(super) fn run_web_pair_bytes(
    args: BlastpArgs,
    query_fasta: &[u8],
    subject_fasta: &[u8],
    query_label: &str,
    subject_label: &str,
) -> Result<Vec<u8>> {
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
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    //     db_length = BlastSeqSrcGetTotLen(seq_src);
    //
    //     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    //
    //     /* iterate over all subject sequences */
    //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //            != BLAST_SEQSRC_EOF) {
    // ...
    //        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //            continue;
    //        }
    // ```
    // NCBI reference: ncbi-blast/c++/include/algo/blast/api/local_blast.hpp:76-78
    // ```c
    // CLocalBlast(CRef<IQueryFactory> query_factory,
    //             CRef<CBlastOptionsHandle> opts_handle,
    //             CRef<CLocalDbAdapter> db);
    // ```
    // Web ABI v1 runs the same local search as the CLI and keeps the report in memory.
    check_web_v1_request(&args)?;
    let subjects = read_subjects(subject_fasta, "BLASTP", true)?;
    let queries = read_queries(query_fasta, "BLASTP", true)?;
    let mut output = Vec::new();
    let outfmt = args.outfmt.clone();
    let mut stderr = std::io::stderr();
    let mut outputs = ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr);
    run_local(
        args,
        &queries,
        &subjects,
        query_label,
        subject_label,
        &mut outputs,
    )?;
    Ok(output)
}
