//! The reports of each query batch of a TBLASTX search written when the batch loop takes the
//! batch, and its final hit lists dropped, as NCBI prints each batch's result sets in its batch
//! loop (line DBID, ID4t fix 2b); otherwise the loop keeps the final lists of every batch and
//! writes all reports after the loop.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_batch_pipe.rs` is), so that it reads that file's types.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-215
//!
//! When the stream acts (`check`), the loop writes a batch's reports right after it has taken
//! the batch (`-db`: its titles read and its prefetch done) and added its queries' warnings: it
//! calls the same report code (`BatchReports::write`) with the batch's final lists, numbered
//! from the batch's first query as the batch made them, its queries' statistics, warnings and
//! ids, and the batch's part of the query ranges (`batch_input`). After the loop it calls it once
//! more with no query and the epilog. The batch's lists and warnings are then dropped.
//!
//! Why the bytes are those of the reports written after the loop (the argument of blastn's
//! `x_rpt_stream.rs`, for TBLASTX):
//! - every report of a query depends only on that query's final list and on fixed inputs. The
//!   lists are final when the batch returns them (`HitCollector::finish`, or `culled_hit_order`
//!   for `-culling_limit`: `-max_hsps`, `-qcov_hsp_perc` and `-max_target_seqs` per query, in
//!   report order; nothing sorts across queries afterwards, ID4t fix 1). The coordinates move by
//!   the query's and the subject's own range start (`shift_to_records`), and the rows, the
//!   outfmt 7 headers and the outfmt 0 sections are written query by query from the query's
//!   list, record, statistics and warnings, the subjects and the options. The writers write the
//!   queries in order, so the reports of consecutive batches written one after the other are the
//!   reports of all of them (outfmt 0 writes a window of queries at a time already, ID4t fix
//!   2a). Formats 1 to 4 are outfmt 0 with another alignment section, made from the query's
//!   aligned hits and records alone (`report::query_anchored`).
//! - a database's reported sequences are numbered by rank among the sequences that the lists of
//!   the call name (`database::reported_subjects`); the rank only indexes the call's records,
//!   titles and ids, and every comparison is already made (the lists are final), so a rank among
//!   one batch's sequences writes the same bytes as one among all.
//! - what is not per query: the outfmt 0 prolog (written before the search, as without the
//!   stream: a file sink, which writes it with the reports, is not streamed), the outfmt 0
//!   epilog (the database lines and the options, fixed inputs) and outfmt 7's last line (the
//!   number of queries formatted, NCBI's `m_QueriesFormatted`: the last call counts the queries
//!   of every batch, `TblastxReportRun::queries_before`). Outfmt 6 has none.
//! - the errors: an error of NCBI's that a batch meets (a reader error, an `eEOF` batch, a batch
//!   without a query or without letters, an unknown matrix, the average subject length, a
//!   database that cannot be read, a header file that cannot be read) comes after the reports of
//!   the batches before, without the epilog, in both: here they were written when the loop took
//!   them, so the error paths write only what follows them (`BatchReports::queries_written`). So
//!   do LOSAT's refusals of what a batch's search meets (`ErrorAfterBatchReports`: a statistics
//!   count of 0 in the setup of the first batch with a valid context, `search_query_batch`, a
//!   search space that overflows NCBI's `Int8`, `context_cutoffs::search_space_overflow`, an
//!   ungapped extension that scores higher past the end of its query frame,
//!   `extension::query_batch`, a one-hit merge of lengths that add up to 0, the ungapped
//!   re-evaluation of HSPs on a subject (or its `-subject_loc` range) of fewer than 12 bases,
//!   `SubjectEmit::refuse_untranslated_frame`, and the re-evaluation's fetch of a subject whose
//!   header cannot be read under a negative seqid list, `subject_loop::x_fetch_keeps`; and since
//!   line RG series D's D-SO-5 a batch whose concatenated query is too long for the diag table,
//!   `check_diag_table_window`, and a `-db` whose sequences could not be read before the search,
//!   where the read error of the first batch with a valid context is LOSAT's refusal: both after
//!   the batch's reader messages too, as NCBI's errors of the batch's search). LOSAT's other
//!   refusals of a batch (the Owner's decision of 2026-11-16 on line RG's Q-X2c: "the results of
//!   the query batches before are written") come after the reports of the batches before and
//!   nothing of the batch, as in the first batch, where they come before its first byte: those of
//!   its reports, made by the batch's report call before it writes (a `-subject` record without a
//!   tabular id, `tabular_subject_ids`; with `-db` a sequence's stored data as the reports read
//!   it, `database::reported_subjects`: the header, the stored title, the tabular id, the parsed
//!   Seq-ids, two sequences that the reports would join, a load that NCBI's places give another
//!   sequence; the refusals of the pairwise formats, `anchored::refuse_pairwise_reports`), and
//!   those that the loop meets between the batch's search and its reports
//!   (`BatchReports::write_before_refusal`: the titles' reading, a `-db` sequence with the id of
//!   a query under `-parse_deflines`, `common::database::refuse_reported_ids_of_queries`, the
//!   reported sequences' reads, the prefetch's taxids and taxonomy database, NCBI's taxonomy
//!   cache, the reports' loads one by one). None is looked for before the search (before
//!   D-SO-5, the stream was off where one could come after the first batch: `Off::Window`,
//!   `SubjectIds`, `Unreadable`, `LargeParsed`, `TaxonomyCache`, `Database`, `QueryIds`, and a
//!   `-db`'s headers were all read before the search). The reports of a batch walk the loads of
//!   NCBI's data loader from those of the batches before (`database::batch_loads_before`: the
//!   run's record, which each call adds its batch to), as the one call's walk holds them. The
//!   refusals that depend on the arguments only (the matrix, the `-seg` value) come in the first
//!   batch, where nothing of the stream is written yet. Formats 1 to 4 (line RG T3c fix 1,
//!   `check`): the exception of NCBI's merge in a query's section (a query record that is also a
//!   `-subject` record, `query_anchored::write_section`) is NCBI's error, after the reports of the
//!   batches before and of the batch's queries before, and the query's report up to its
//!   description table, without the epilog, as NCBI's; the refusal of a `-db` sequence whose
//!   generated id (`lcl|Query_N`) is that of a query of its batch (NCBI's scope is the batch's,
//!   `ResetScopeHistory`) comes before that batch's first byte, after the batches before. With
//!   one batch (the first batch takes every query record) the one call is that batch's. LOSAT's
//!   internal invariant errors and panics are not argued; a write error of the output stops the
//!   run where it happens (at the end without the stream).
//! - what the report path writes to standard error (the query warnings) keeps its bytes and its
//!   place before each query's report. `LOSAT_TIMING`'s blocks are printed when the loop takes a
//!   batch, as before, so with the stream a block comes after the warnings of the batches before
//!   it (the same lines in another order). The debugging output of the search, which a batch
//!   writes while it is searched, keeps the reports at the end (`Off::Debug`), where the loop
//!   writes them one batch at a time (`Writes::Held`, `BatchReports::write_batches`): the calls
//!   of the stream, so a batch's refusal comes after the batches before in both.
//!
//! Off (the reports after the loop, in one call): outputs other than one format of the command
//! line or ABI v1 (the hit records and the format observer of the Web adapter, whose HSP numbers
//! count over the whole search; a file sink, created when its report is written; two or more
//! formats, which the library API allows on one stream); formats other than 0 to 4, 6 and 7 (the
//! BLAST XML documents, the Seq-annots and the archive, whose framing and refusals are not argued
//! here); the tabular reports with a field list or a delimiter and outfmt 10 and 20
//! (`Off::CustomTabular`: their refusals of a subject, a name of the taxonomy database or a
//! definition line's ids, must come before the first byte, DECISIONS 5 of line CLI C2); one batch.
//! Their refusals of a later batch come before the first byte, as before; NCBI's errors and
//! LOSAT's refusals of a batch's search (`ErrorAfterBatchReports`) after the batches before.
//!
//! Diagnostics: with `LOSAT_X_STATS`, `[X_STATS] TBXRPTSTREAM off= check_s=` before the search
//! (`off`: `none` when the loop writes each batch's reports, else the reason, `Off::name`;
//! `check_s`: the seconds of `check`).
//!
//! Memory: a batch's lists, statistics and warnings are dropped after its write. Kept for the
//! whole search, as without the stream: one (empty) warnings slot per query, and with `-db` the
//! header of every distinct sequence that a batch reported (`database::DatabaseTitles`, at most
//! the database's sequences and at most the queries times the hit list size).

use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::seq_range::{next_ranged_batch, Placements, QueryInput};

use super::{report, writes_debug_output};

/// Why the reports of a search are written after its batch loop (`check`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The first batch takes every query record: the one call after the loop writes them.
    OneBatch,
    /// Hit records, a format observer, a file sink, or two or more formats.
    Outputs,
    /// A format other than 0 to 4, 6 and 7.
    Format,
    /// A tabular report with a field list or a delimiter, or outfmt 10 or 20
    /// (`TblastxOutputFormat::CustomTabular`), which is written after the loop (the module
    /// doc).
    CustomTabular,
    /// The debugging output of the search: the loop holds the reports and writes them after it,
    /// one batch at a time (`held`).
    Debug,
}

impl Off {
    /// The name of `[X_STATS] TBXRPTSTREAM off=`.
    pub(super) fn name(self) -> &'static str {
        match self {
            Off::OneBatch => "one_batch",
            Off::Outputs => "outputs",
            Off::Format => "format",
            Off::CustomTabular => "custom_tabular",
            Off::Debug => "debug",
        }
    }
}

/// Whether the loop writes each batch's reports when it takes the batch (`Ok`), or why not
/// (`input`: the search's query input, `batch_size`: the first batch's size).
pub(super) fn check(
    outputs: &ReportOutputs<'_>,
    input: &QueryInput,
    batch_size: u32,
) -> Result<(), Off> {
    let (first_end, _) = next_ranged_batch(&input.input_lengths, &input.skipped, 0, batch_size);
    if first_end == input.input_lengths.len() {
        return Err(Off::OneBatch);
    }
    let [format] = outputs.formats.as_slice() else {
        return Err(Off::Outputs);
    };
    if outputs.hits.is_some()
        || outputs.observer.is_some()
        || matches!(format.sink, OutputSink::File(_))
    {
        return Err(Off::Outputs);
    }
    let format = report::output_format(format.outfmt);
    // The custom tabular reports (line CLI C2, PR 2): their writer refuses a name of the
    // taxonomy database for a subject with a taxid (DECISIONS 3) and reads the taxids and
    // definition lines of the reported sequences, which a later batch's report call would do
    // after the first batch's rows were written; and outfmt 20's header line comes with the
    // first report. They keep the one report call after the search, where every refusal comes
    // before the first byte, and never take the streamed path of the default 6 and 7.
    if let report::TblastxOutputFormat::CustomTabular(_) = format {
        return Err(Off::CustomTabular);
    }
    // Outfmt 1 to 4 are written per batch as outfmt 0 (line RG T3c fix 1): NCBI's exception of a
    // section's merge comes where NCBI's does (`query_anchored::write_section`), and the refusal of
    // a `-db` sequence whose generated id is that of a query of its batch
    // (`anchored::refuse_pairwise_reports`) before that batch's first byte, after the batches
    // before, as NCBI's error in a later batch's format (DECISIONS T3a-7 (2) of line RG, BLASTN's
    // and BLASTP's per batch and block).
    if !(format.is_pairwise() || format.is_tabular()) {
        return Err(Off::Format);
    }
    if writes_debug_output() {
        return Err(Off::Debug);
    }
    Ok(())
}

/// How the batch loop writes the reports of a search's batches (`check`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum Writes {
    /// Each batch's reports when the loop takes the batch (`check` is `Ok`).
    Each,
    /// After the loop, one batch at a time, as `Each` writes them
    /// (`BatchReports::write_batches`): a streamable output that only the debugging output keeps
    /// after the loop (`Off::Debug`).
    Held,
    /// After the loop, in one call (the other reasons of `Off`).
    #[default]
    Whole,
}

impl Writes {
    /// The writes of a search whose stream is `stream` (`check`).
    pub(super) fn of(stream: Result<(), Off>) -> Self {
        match stream {
            Ok(()) => Writes::Each,
            Err(Off::Debug) => Writes::Held,
            Err(_) => Writes::Whole,
        }
    }
}

/// The part of the query input of the searched queries `searched`, numbered from 0 as a batch's
/// lists number them (blastn's `x_rpt_stream::batch_input`).
pub(super) fn batch_input(input: &QueryInput, searched: std::ops::Range<usize>) -> QueryInput {
    QueryInput {
        placements: if input.placements.ranged() {
            Placements(input.placements.0[searched.clone()].to_vec())
        } else {
            Placements::default()
        },
        ordinals: input.ordinals[searched].to_vec(),
        input_lengths: Vec::new(),
        skipped: Vec::new(),
        range_length: input.range_length,
    }
}

#[cfg(test)]
#[path = "report_stream_tests.rs"]
mod tests;
