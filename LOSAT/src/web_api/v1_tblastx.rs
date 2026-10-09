//! ABI v1's TBLASTX search (plan TD-1): `bio` reads the inputs, ABI v1's checks of them come
//! where it made them (`V1Records`, the engine's `V1Checks`), and the records that pass
//! them enter the search as the records of NCBI's reader (`v1_bio::from_bio`). Moved
//! unchanged from `algorithm/tblastx/blast_engine/run_impl.rs` in session SFc, S10.

use anyhow::{Context, Result};
use bio::io::fasta;

use super::v1_bio::from_bio;
use crate::algorithm::tblastx::blast_engine::{
    check_losat_limits, report, run_local_with, V1Checks,
};
use crate::algorithm::tblastx::TblastxArgs;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, parse_optional_range, Placements, RangeRole, SequenceRange,
};

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
        // NCBI reference: c++/src/objtools/readers/fasta.cpp:428-431
        // FASTA_ERROR(LineNumber(), "CFastaReader: Expected defline around line " << LineNumber(), ...);
        .collect::<std::result::Result<Vec<_>, _>>()
        .context("failed to parse in-memory FASTA")
}

pub(super) fn run_web_pair(
    args: TblastxArgs,
    query_fasta: &str,
    subject_fasta: &str,
) -> Result<Vec<u8>> {
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
    let mut stderr = std::io::stderr();
    // ABI v1 keeps `bio` and its checks (plan TD-1): the checks that it made in the search
    // come where they came (`V1Records`), and the records that pass them, which NCBI's
    // reader reads alike, enter the search as its records (`FastaRecord::from_bio`, which
    // reads `U` as `T`).
    let query_records: Vec<FastaRecord> = queries
        .iter()
        .enumerate()
        .map(|(index, record)| from_bio(record, index + 1, "Query_", false))
        .collect();
    let subject_records: Vec<FastaRecord> = subjects
        .iter()
        .enumerate()
        .map(|(index, record)| from_bio(record, index + 1, "Subject_", false))
        .collect();
    run_local_with(
        args,
        &query_records,
        &subject_records,
        &mut ReportOutputs::single(&outfmt, OutputSink::Writer(&mut output), &mut stderr),
        Some(&V1Records {
            queries: &queries,
            subjects: &subjects,
        }),
    )?;
    Ok(output)
}

/// ABI v1's records as `bio` reads them (`run_web_pair`). Plan TD-1 freezes ABI v1's
/// accepted inputs, its messages and their order: the search makes ABI v1's checks of
/// these records where it made them when it searched them (`check_search`, and
/// `check_shown_subject_titles` before the reports), and searches the records of NCBI's
/// reader made from them (`from_bio`), which those checks guarantee to be
/// read alike.
struct V1Records<'a> {
    queries: &'a [fasta::Record],
    subjects: &'a [fasta::Record],
}

impl V1Checks for V1Records<'_> {
    /// ABI v1's checks where the search starts, in their order: the outfmt 0 and 7 titles,
    /// LOSAT's limits and environment, the residues, the records without residues and the
    /// intervals without letters (the subjects cut to `-subject_loc`, the queries to
    /// `-query_loc`). NCBI reads those inputs without a message (or with a warning about
    /// removed residues), and the records of NCBI's reader search them as NCBI does; ABI v1
    /// rejects them as it did.
    fn check_search(
        &self,
        args: &TblastxArgs,
        query_range: Option<&SequenceRange>,
        outputs: &ReportOutputs<'_>,
    ) -> Result<()> {
        use super::v1_bio::{check_records_have_residues_of, check_residues_of};
        use crate::blastinput::seq_range::check_no_empty_interval;
        // The ranges cut records of NCBI's reader (`seq_range`); these have the residues of
        // ABI v1's records and their `bio` ID as the title, which ABI v1's messages name them
        // by. The cut `bio` subjects keep their ID and description.
        let named_by_id = |records: &[fasta::Record]| -> Vec<FastaRecord> {
            records
                .iter()
                .map(|record| FastaRecord::new(String::new(), record.id().as_bytes(), record.seq()))
                .collect()
        };
        let subject_range =
            parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "TBLASTX")?;
        let whole_named_subjects = named_by_id(self.subjects);
        let ranged_subjects = cut_subjects(&whole_named_subjects, subject_range.as_ref())?;
        let whole_subjects = Placements::default();
        let cut_subjects_bio: Vec<fasta::Record>;
        let (subjects, named_subjects, subject_placements) = match &ranged_subjects {
            Some((cut, placements)) => {
                cut_subjects_bio = self
                    .subjects
                    .iter()
                    .zip(cut)
                    .map(|(record, cut)| {
                        fasta::Record::with_attrs(record.id(), record.desc(), cut.seq())
                    })
                    .collect();
                (cut_subjects_bio.as_slice(), cut.as_slice(), placements)
            }
            None => (
                self.subjects,
                whole_named_subjects.as_slice(),
                &whole_subjects,
            ),
        };
        check_report_titles(self.queries, subjects, outputs)?;
        check_losat_limits(args)?;
        crate::blastinput::app::check_unsupported_environment("TBLASTX")?;
        check_residues_of(subjects, "subject", "TBLASTX")?;
        check_no_empty_interval(
            named_subjects,
            subject_placements,
            &[],
            RangeRole::Subject,
            "TBLASTX",
        )?;
        for (records, role) in [(subjects, "subject"), (self.queries, "query")] {
            check_residues_of(records, role, "TBLASTX")?;
            check_records_have_residues_of(records, role, "TBLASTX")?;
        }
        if let Some(range) = query_range {
            let ranged = cut_queries(&named_by_id(self.queries), range);
            check_no_empty_interval(
                &ranged.records,
                &ranged.input.placements,
                &ranged.input.ordinals,
                RangeRole::Query,
                "TBLASTX",
            )?;
        }
        Ok(())
    }
}

/// ABI v1's limits on the titles that outfmt 0 and 7 print (`V1Records`).
///
/// ABI v1 reads the deflines with `bio`, which does not reproduce NCBI's `CFastaReader`
/// for every defline (the title ends at the first byte below a space,
/// fasta_reader_utils.cpp:215-225). Those reports therefore reject a defline that is
/// empty, starts with white space, or has a control character or a non-ASCII byte. The
/// outfmt 0 titles of the subjects are checked where the report shows them
/// (`check_shown_subject_titles`).
fn check_report_titles(
    queries: &[fasta::Record],
    subjects: &[fasta::Record],
    outputs: &ReportOutputs<'_>,
) -> Result<()> {
    if outputs
        .formats
        .iter()
        .all(|format| report::output_format(format.outfmt) == report::TblastxOutputFormat::Tabular)
    {
        return Ok(());
    }
    for (role, records) in [("query", queries), ("subject", subjects)] {
        for (index, record) in records.iter().enumerate() {
            let defline = match record.desc() {
                Some(desc) => format!("{} {desc}", record.id()),
                None => record.id().to_string(),
            };
            if defline.is_empty()
                || !defline.is_ascii()
                || defline.bytes().any(|byte| byte < b' ')
                || record.id().is_empty()
            {
                anyhow::bail!(
                    "{role} record {} has a defline that is empty, starts with white space or has a control character or a non-ASCII byte, which NCBI BLAST+ reads differently in the outfmt 0 and 7 titles; this is not supported by LOSAT's TBLASTX",
                    index + 1
                );
            }
        }
    }
    Ok(())
}
