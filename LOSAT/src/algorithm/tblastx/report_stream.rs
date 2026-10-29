//! The reports of each query batch of a TBLASTX search written when the batch loop takes the
//! batch, and its final hit lists dropped, as NCBI prints each batch's result sets in its batch
//! loop (line DBID, ID4t fix 2b); otherwise the loop keeps the final lists of every batch and
//! writes all reports after the loop.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_batch_pipe.rs` is), so that it reads that file's types.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-215
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!     ...
//!             CRef<CSearchResultSet> results = lcl_blast->Run();
//!     ...
//!                 ITERATE(CSearchResultSet, result, *results) {
//!                     formatter.PrintOneResultSet(**result, query_batch);
//!                 }
//!     ...
//!         }
//!         formatter.PrintEpilog(opt);
//! ```
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
//!   2a).
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
//!   them, so the error paths write only what follows them (`BatchReports::queries_written`). An
//!   error that LOSAT raises where NCBI would not stop (LOSAT's refusals) is written without the
//!   reports of the batches before; so the stream acts only where no such refusal can come after
//!   the first batch (`check`): a query batch whose concatenated query is too long for the diag
//!   table (`check_diag_table_window`), a search space that overflows (`SearchSpaceOverflow`),
//!   a `-db` whose sequences could not be read before the search (some of those read errors are
//!   LOSAT's refusals, raised by the first batch with a valid context, `Off::Unreadable`),
//!   and the refusals of the reports: a `-subject` record without a tabular id
//!   (`tabular_subject_ids`), and with `-db` every sequence's stored data as the reports read it
//!   (`TblastxDatabase::reports_may_refuse`: the header, the stored title, the tabular id, the
//!   parsed Seq-ids, two sequences that the reports would join, the prefetch's taxids and
//!   taxonomy database). The refusals that depend on the arguments only (the matrix, the
//!   `-seg` value) come in the first batch, where nothing of the stream is written yet. With one
//!   batch (the first batch takes every query record) nothing can come after it, so the stream
//!   does not act and nothing is checked. LOSAT's internal invariant errors and panics are not
//!   argued; a write error of the output stops the run where it happens (at the end without the
//!   stream).
//! - what the report path writes to standard error (the query warnings) keeps its bytes and its
//!   place before each query's report. `LOSAT_TIMING`'s blocks are printed when the loop takes a
//!   batch, as before, so with the stream a block comes after the warnings of the batches before
//!   it (the same lines in another order). The debugging output of the search, which a batch
//!   writes while it is searched, keeps the reports at the end (`Off::Debug`).
//!
//! Off (the reports after the loop): outputs other than one format of the command line or ABI v1
//! (the hit records and the format observer of the Web adapter, whose HSP numbers count over the
//! whole search; a file sink, created when its report is written; two or more formats, which the
//! library API allows on one stream); formats other than 0, 6 and 7 (the BLAST XML documents,
//! the Seq-annots and the archive, whose framing and refusals are not argued here); ABI v1's
//! outfmt 0, which refuses a reported subject's title (`check_shown_subject_titles`); the
//! debugging output; and the refusals above.
//!
//! Memory: a batch's lists, statistics and warnings are dropped after its write. Kept for the
//! whole search, as without the stream: one (empty) warnings slot per query, and with `-db` the
//! header of every distinct sequence that a batch reported (`database::DatabaseTitles`, at most
//! the database's sequences and at most the queries times the hit list size).

use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::seq_range::{next_ranged_batch, Placements, QueryInput};
use crate::utils::threading::SearchPool;

use super::{
    report, subject_tabular_id, writes_debug_output, QueryRuns, RunSubjects, ShownSubjects,
    TblastxArgs,
};

/// Why the reports of a search are written after its batch loop (`check`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The first batch takes every query record: the one call after the loop writes them.
    OneBatch,
    /// Hit records, a format observer, a file sink, or two or more formats.
    Outputs,
    /// A format other than 0, 6 and 7.
    Format,
    /// ABI v1's outfmt 0, which refuses a reported subject's title.
    V1Titles,
    /// The debugging output of the search.
    Debug,
    /// A concatenated query of a batch that the diag table could refuse.
    Window,
    /// A search space that could overflow.
    SearchSpace,
    /// A `-subject` record without a tabular id.
    SubjectIds,
    /// A `-db` whose sequences could not be read before the search: LOSAT refuses some of those
    /// errors where NCBI's read error stops the run, at the first batch with a valid context,
    /// which may follow a batch without one (`run_batch_loop`).
    Unreadable,
    /// A database sequence whose report could be refused (`ReportRefusal`).
    Database(crate::algorithm::common::database::ReportRefusal),
}

/// Whether the loop writes each batch's reports when it takes the batch (`Ok`), or why not.
/// With `-db` it reads every header of the database in `pool`
/// (`TblastxDatabase::reports_may_refuse`), unless an earlier condition is off.
#[allow(clippy::too_many_arguments)]
pub(super) fn check(
    args: &TblastxArgs,
    outputs: &ReportOutputs<'_>,
    v1_titles: bool,
    query_runs: &QueryRuns<'_>,
    subjects: &RunSubjects<'_>,
    unreadable: bool,
    batch_size: u32,
    pool: &SearchPool<'_>,
) -> Result<(), Off> {
    let input = query_runs.input;
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
    let pairwise = format == report::TblastxOutputFormat::Pairwise;
    if !(pairwise || format.is_tabular()) {
        return Err(Off::Format);
    }
    if v1_titles && pairwise {
        return Err(Off::V1Titles);
    }
    if writes_debug_output() {
        return Err(Off::Debug);
    }
    if unreadable {
        return Err(Off::Unreadable);
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-315
    // A batch's concatenated query holds each query's six frames, each at most a third of the
    // query's letters, after a sentinel byte: at most `2 * letters + 7` per query. All queries
    // bound every batch (`check_diag_table_window`).
    let lengths = || {
        query_runs
            .records
            .iter()
            .map(|record| record.seq().len() as i128)
    };
    let concatenated: i128 = lengths().map(|letters| 2 * letters + 7).sum::<i128>() + 1;
    if concatenated + i128::from(args.window()) > 1 << 30 {
        return Err(Off::Window);
    }
    // The search space of a context (`search_query_batch`, `try_compute_eff_lengths_tblastx`):
    // the effective database length `db - n a`, at least 1, times the context's length minus
    // the length adjustment `a`. The adjustment is between 0 and the smaller root of
    // `n a^2 - (q n + db) a + c`, at most `2 c / (q n + db) <= 2 q` with `c <= q db`
    // (`compute_length_adjustment_ncbi`; the bound of blastn's `x_rpt_stream::check`), and a
    // context has at most a third of its query's letters, plus one.
    let db_length = i128::from(subjects.numbers.effective_length() / 3);
    let factor = 2 * (lengths().max().unwrap_or(0) / 3 + 1) + 1;
    if (db_length.abs() + i128::from(subjects.numbers.count) * factor) * factor
        > i128::from(i64::MAX)
    {
        return Err(Off::SearchSpace);
    }
    match &subjects.shown {
        ShownSubjects::Fasta { records, .. } => {
            if format.is_tabular()
                && records
                    .iter()
                    .any(|record| subject_tabular_id(record).is_none())
            {
                return Err(Off::SubjectIds);
            }
            Ok(())
        }
        ShownSubjects::Database { database, .. } => {
            let shown = crate::algorithm::common::database::ParsedShown {
                tabular: format.is_tabular(),
                pairwise,
                joined: pairwise,
                sam: false,
            };
            // The batches' prefetch reads the taxonomy database (`run_batch_loop`).
            let prefetch = format.prefetches_subjects() && database.taxdb_found();
            match database.reports_may_refuse(shown, prefetch, pool) {
                Some(refusal) => Err(Off::Database(refusal)),
                None => Ok(()),
            }
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
