//! ABI v1's BLASTP search (plan TD-1): `bio` reads the inputs (or the FASTA store keeps
//! them), ABI v1's checks of them come where it made them (`V1Records`, the engine's
//! `V1Checks`), and the records that pass them enter the search as the records of NCBI's
//! reader (`v1_bio::from_bio`). Moved unchanged from `algorithm/blastp/blast_engine.rs` in
//! session SFc, S10.

use anyhow::{bail, Context, Result};
use bio::io::fasta;

use super::v1_bio::from_bio;
use crate::algorithm::blastp::blast_engine::{
    blastp_tabular_field, check_options, run_local_with, V1Checks,
};
use crate::algorithm::blastp::BlastpArgs;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::OutputFormat;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:609-616
// ```c
//                 sequence = queries.GetBlastSequence(index,
//                                                     encoding,
//                                                     eNa_strand_unknown,
//                                                     eSentinels,
//                                                     &warnings);
//
//                 int offset = qinfo->contexts[ctx_index].query_offset;
//                 memcpy(&buf.get()[offset], sequence.data.get(),
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
fn fasta_records_from_bytes(bytes: &[u8]) -> Result<Vec<fasta::Record>> {
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
    fasta::Reader::new(bytes)
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(anyhow::Error::from)
}

pub(super) fn run_web_pair(
    args: BlastpArgs,
    query_fasta: &str,
    subject_fasta: &str,
) -> Result<Vec<u8>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:609-616
    // ```c
    //                 sequence = queries.GetBlastSequence(index,
    //                                                     encoding,
    //                                                     eNa_strand_unknown,
    //                                                     eSentinels,
    //                                                     &warnings);
    //
    //                 int offset = qinfo->contexts[ctx_index].query_offset;
    //                 memcpy(&buf.get()[offset], sequence.data.get(),
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
    let queries = fasta_records_from_bytes(query_fasta.as_bytes())
        .context("failed to parse web query FASTA")?;
    let subjects = fasta_records_from_bytes(subject_fasta.as_bytes())
        .context("failed to parse web subject FASTA")?;
    run_web_pair_records(
        args,
        (&queries, query_fasta.as_bytes()),
        (&subjects, subject_fasta.as_bytes()),
        "",
        "",
    )
}

/// ABI v1's checks of a BLASTP request, whose arguments, formats and errors are frozen
/// (plan TD-1): the options are checked before the output format, and a tabular field that
/// LOSAT's BLASTP does not write is an error, as in v1 before S08+ moved BLASTP to NCBI's
/// application layer (which reads `-outfmt` first and ignores a token that is not a field).
fn check_web_v1_request(args: &BlastpArgs) -> Result<()> {
    check_options(args)?;
    let (format, fields) = OutputFormat::parse(&args.outfmt).map_err(anyhow::Error::msg)?;
    if let Some(fields) = fields {
        if format == OutputFormat::Pairwise {
            bail!("blastp outfmt 0 does not accept custom field lists");
        }
        for token in fields.split_whitespace() {
            if !token.eq_ignore_ascii_case("std")
                && !matches!(blastp_tabular_field(token), Ok(Some(_)))
            {
                bail!("unsupported blastp tabular field '{token}'");
            }
        }
    }
    Ok(())
}

/// ABI v1's BLASTP search of `bio`'s records of the query and subject FASTA bytes (each
/// role's records with the bytes that they were read from).
pub(super) fn run_web_pair_records(
    args: BlastpArgs,
    (query_records, query_fasta): (&[fasta::Record], &[u8]),
    (subject_records, subject_fasta): (&[fasta::Record], &[u8]),
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
    let shown = v1_shown(&args.outfmt);
    let mut output = Vec::new();
    let outfmt = args.outfmt.clone();
    let mut stderr = std::io::stderr();
    let mut outputs = ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr);
    // ABI v1 keeps `bio` and its checks (plan TD-1): the checks that it made in the search
    // come where they came (`V1Records`), and the records that pass them enter the search
    // as NCBI's reader's records (`FastaRecord::from_bio`, which skips the white space at the
    // start of a title as NCBI's defline parser does; ABI v1's BLASTP checks no defline
    // but those to which `bio` gives an empty ID, after the search).
    let v1_record = |index: usize, record: &fasta::Record, prefix: &str| {
        from_bio(record, index + 1, prefix, true)
    };
    let queries: Vec<FastaRecord> = query_records
        .iter()
        .enumerate()
        .map(|(index, record)| v1_record(index, record, "Query_"))
        .collect();
    let subjects: Vec<FastaRecord> = subject_records
        .iter()
        .enumerate()
        .map(|(index, record)| v1_record(index, record, "Subject_"))
        .collect();
    run_local_with(
        args,
        &queries,
        &subjects,
        query_label,
        subject_label,
        &mut outputs,
        Some(&V1Records {
            queries: query_records,
            subjects: subject_records,
        }),
    )?;
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-213
    // ```c++
    //     // trim leading whitespace from title (is this appropriate?)
    //     while (title_start < len
    //         &&  isspace((unsigned char)defline[title_start])) {
    //         ++title_start;
    //     }
    // ```
    // The deflines to which `bio` gives an empty ID and whose records NCBI's reader reads
    // otherwise are rejected (`check_bio_deflines_of`), subjects first, as ABI v1's
    // BLASTN rejects its deflines, and only for a search (a query input without records
    // gives NCBI's `Query is Empty!`). The rejection comes after ABI v1's other checks,
    // which plan TD-1 freezes with their order, so that it changes no other v1 error.
    if !query_records.is_empty() {
        use super::v1_bio::check_bio_deflines_of;
        check_bio_deflines_of(subject_fasta, "subject", "BLASTP", shown.subjects)?;
        check_bio_deflines_of(query_fasta, "query", "BLASTP", shown.queries)?;
    }
    Ok(output)
}

/// What each role of an ABI v1 BLASTP search shows (`check_bio_deflines_of`), by its
/// output format, which `check_web_v1_request` has accepted: the subjects' names (a subject
/// without a title is `unnamed`) and their titles in outfmt 0, and the queries' local IDs in
/// a tabular format with a query ID field (`std` has `qaccver`); the outfmt 0 and 7 query
/// lines show the title only.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-41
/// ```c++
/// const char* kDfltArgTabularOutputFmt =
///     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
///     "evalue bitscore";
/// const char* kDfltArgTabularOutputFmtTag("std");
/// ```
fn v1_shown(outfmt: &str) -> V1Shown {
    use super::v1_bio::Shown;
    let Ok((format, fields)) = OutputFormat::parse(outfmt) else {
        return V1Shown::default();
    };
    let query_ids = format != OutputFormat::Pairwise
        && fields.is_none_or(|fields| {
            fields.split_whitespace().any(|token| {
                token.eq_ignore_ascii_case("std") || matches!(token, "qseqid" | "qacc" | "qaccver")
            })
        });
    V1Shown {
        subjects: Shown {
            names: true,
            local_ids: false,
            outfmt0_titles: format == OutputFormat::Pairwise,
        },
        queries: Shown {
            names: query_ids,
            local_ids: query_ids,
            outfmt0_titles: false,
        },
    }
}

/// `v1_shown`'s answer for the subjects and the queries.
#[derive(Default)]
struct V1Shown {
    subjects: super::v1_bio::Shown,
    queries: super::v1_bio::Shown,
}

/// ABI v1's records as `bio` reads them (`run_web_pair_records`). Plan TD-1 freezes ABI
/// v1's accepted inputs, its messages and their order: the search makes ABI v1's checks of
/// these records where it made them when it searched them (`check_records`, the engine's
/// query-interval check, and its `check_shown_subject_titles` before the outfmt 0 report),
/// and searches the records of NCBI's reader made from them (`from_bio`), which those
/// checks guarantee to be read alike.
struct V1Records<'a> {
    queries: &'a [fasta::Record],
    subjects: &'a [fasta::Record],
}

impl V1Checks for V1Records<'_> {
    /// ABI v1's checks of the residues, after the warnings of the subjects without
    /// residues: a residue that NCBI's protein reader removes (subjects, then queries), and
    /// a query without residues. NCBI reads those inputs with a warning (or reports a query
    /// without residues), and ABI v1 rejects them as it did.
    fn check_records(&self) -> Result<()> {
        use super::v1_bio::{check_protein_residues_of, check_records_have_residues_of};
        check_protein_residues_of(self.subjects, "subject", "BLASTP")?;
        check_protein_residues_of(self.queries, "query", "BLASTP")?;
        check_records_have_residues_of(self.queries, "query", "BLASTP")
    }
}
