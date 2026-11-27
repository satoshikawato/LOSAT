//! EXPERIMENT (LOSAT_X_RPTSTREAM): blastn writes the reports of
//! each query batch when the batch loop of `run_in_pool` takes the batch, and drops its final hit
//! lists, as NCBI prints each batch's result sets in its batch loop; without the switch the loop
//! keeps the final lists of every batch and writes all reports after the loop.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-318
//!
//! With the switch on, the loop writes a batch's reports right after it has searched the batch (in
//! batch order), finished its lists (`-db`: `finish_query_hit_list`, the titles, the HSPs of score
//! 0) and given its gie to the mixer. It calls the same report code
//! (`write_reports`) with the batch's queries, their final lists numbered from 0 (as the batch
//! made them), their warnings, lengths, masks, Karlin-Altschul blocks and search spaces, and the
//! batch's part of the query ranges; after the loop it calls it once more with no query and the
//! epilog. The lists are then dropped (the batch's other per-query values are not kept either).
//!
//! Why the bytes are those of the reports written after the loop:
//! - every report of a query depends only on that query's final list and on fixed inputs: the
//!   lists are finished per query (`finish_query_hit_list`: `-max_hsps`, `-qcov_hsp_perc`,
//!   `-subject_besthit`, the sort, `-max_target_seqs`; the HSPs of score 0), the coordinates move
//!   by the query's and the
//!   subject's own range start (`shift_to_records`), and the rows, outfmt 7 headers and outfmt 0
//!   sections are written query by query from the query's list, record, masks, blocks and search
//!   space, the subject records and the options. BLASTN has no cross-query sort, and its culling
//!   (`-culling_limit`, `culling.rs`) is NCBI's pipe, which culls each query's tree on its own and is the first
//!   step of `finish_query_hit_list`, as is its best hit pipe (`-best_hit_*`, `besthit.rs`, a list
//!   per query).
//!   The writers write the queries in order, so the reports of consecutive batches written one
//!   after the other are the reports of all of them. A query's warnings come before its report
//!   in both (`QueryWarnings`), with the first format (or all at once without a format).
//! - a database's reported sequences are numbered by rank among the sequences that the lists of
//!   the call name (`database::reported_subjects`); the rank only indexes the call's records and
//!   ids, and every comparison that the reports make is already done (the lists are finished), so
//!   a rank among one batch's sequences writes the same bytes as one among all.
//! - what is not per query: the outfmt 0 prolog (written before the search, as without the
//!   switch; a sink that writes it with the reports is a file sink, which is not streamed), the
//!   outfmt 0 epilog (the database lines and the options, fixed inputs), outfmt 7's last line
//!   (the number of queries formatted, NCBI's `m_QueriesFormatted`: the last call counts the
//!   queries of every batch, `BlastnReportInputs::queries_before`), and outfmt 6 has none. The
//!   formats 8, 9 and 12 (OF3) write per query like 6. BLAST XML (5) and the XML2 and JSON
//!   documents (16 and 15 on the output, 14 and 13 in side files) write per query as NCBI's
//!   formatter writes them (NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-748):
//!   outfmt 5's header comes with the run's first report (`queries_before` 0) and its iterations
//!   are numbered over the run (`BlastXmlWriter::new`); the document's state of 13 to 16 (its
//!   start, the number of its reports, which JSON's separator and the side files' names read, and
//!   the held report of a `-subject` run, which NCBI makes when the next query's result arrives, in
//!   the next batch or in the epilog: blast_format.cpp:2362-2403) is handed from call to call
//!   (`Xml2Carry`, `Stream::xml2`), as the BLASTP block search hands it; the footer, the end and
//!   the master file come with the epilog call. Under `-parse_deflines` each call shows its
//!   queries by the run's scopes (`QueryScopes::numbered_from`: which results join the held
//!   report is the run's, the labels the call's). The tabular reports with a field list or a
//!   delimiter and outfmt 10 and 20 (`algorithm::blastn::tabular`) write per query like 6 and 7
//!   (the delimiter and the kind change only how a row is written; `qcovs` and `qcovus` are a
//!   query's and a subject's; the texts of the rows' lookups of the taxonomy database are a
//!   taxid's), and outfmt 20's line of the field names comes with the first report
//!   (`queries_before` 0: the first batch's call), or, where an error of NCBI's stops the first
//!   batch, from `write_before_failed_batch` as without the switch. The texts of opening the
//!   taxonomy database that a batch's rows keep (`TaxdbRun::open_for_rows`) go to the first
//!   report call whose rows look a name up (`common::database::tabular_names`), as they go to the
//!   first such row of the one call.
//! - the errors: an error of NCBI's that a batch meets (a reader error, an `eEOF` batch, a batch
//!   of queries without letters, `after_batch_reports`, a header file that cannot be read)
//!   comes after the reports of the batches before, without the epilog, in both: here they were
//!   written when the loop took them, so the error paths write only what follows them (the
//!   batch's title and reader warnings). Search space overflows
//!   (`context_statistics::setup_failure` -> `ErrorAfterBatchReports`) follow the same path, so
//!   both stream on and stream off write the reports of the batches before and then the refusal,
//!   matching NCBI. Tabular row stops (`FailedLookupRefusal`, including a subject set's failed
//!   tabular id from `tabular_subject_ids`) are met while writing the rows in query order
//!   (`report::tabular::row_stop_end`), so both stream on and stream off write the reports of the
//!   batches before, the stopping batch's posts and report bytes up to the stopping field, and
//!   then the refusal, matching NCBI. Any other error that LOSAT
//!   raises where NCBI would not stop (LOSAT's other refusals: stored title, ambiguity data,
//!   prefetch, parsed Seq-id, taxonomy record, and the like) is written without the reports of the
//!   batches before by the off path, and after them with the switch (drop-in DECISIONS 44: as NCBI
//!   reads them, the stored data and the report ids of a batch's reported subjects are checked in
//!   that batch, nothing is read before the search). These are the refusals of a reported subject:
//!   with `-db` a definition line that does not read as the reports read it, holds a byte outside
//!   ' ' ..= '~' or gives no tabular subject id (`database::reported_subjects`), with parsed
//!   Seq-ids a header that the reports refuse and two sequences that outfmt 0 would join
//!   (`check_reported_parsed`, `refuse_joined_hits`), one that the batches' prefetch reads where
//!   no report reads it (the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`, a table of fields
//!   that read nothing of the subject) and cannot read as NCBI does (`check_prefetched_parsed`),
//!   ambiguity data that the traceback cannot decode (`traceback_listed`), a taxid or a lookup of
//!   the taxonomy database that the prefetch refuses (`batch_prefetch`), and what a custom tabular
//!   report reads (the taxonomy fields, the names of the taxonomy database, `sall*`); with a
//!   subject set a report id (`tabular_subject_ids`; the ids are made once where LOSAT refuses
//!   none, `SubjectCache`) and with `-parse_deflines` the parsed ids that a custom tabular report
//!   shows. Each is made by the report call (or the search) of the batch before that batch writes
//!   a byte, so the output is the reports of the batches before (without the epilog) and their
//!   query warnings, then the refusal and its exit status (the search's refusals, of the titles'
//!   reading, the traceback and the prefetch, come after the batch's title warnings, as without
//!   the switch; where several sequences are refused, the one named is the batch's first, not the
//!   search's). The parsed ids of `-parse_deflines` without `sall*` are refused while the rows are
//!   made (`TabularSubjects::values`), after the rows of the queries before in both paths, so
//!   their bytes do not change. An error of NCBI's of the taxonomy database or of reading the
//!   database that a later batch's prefetch or report meets likewise comes after the reports of
//!   the batches before, as NCBI's loop writes them (without the switch after none); any other
//!   input that LOSAT does not refuse gives the same bytes either way. Outfmt 17 writes its
//!   document only at the end (`SamDocument::finish`), so a refusal in a later batch finds nothing
//!   written there either.
//!   The XML family has refusals of its own (a `-db` name, query title or id, or a title of a
//!   subject with a hit that cannot be written, `blast_xml.rs`). Outfmt 5 and the `-db` XML2 and
//!   JSON documents make a batch's reports in its call, so their refusals come before the call's
//!   first byte: after the reports of the batches before, without the footer, the end or the
//!   master file (the Owner's decision of 2026-11-16 on line RG's Q-X2c, "the results of the
//!   query batches before are written"; remaining-gaps DECISIONS D-SO-2). A `-subject` XML2 or JSON
//!   document makes a report where NCBI makes it (when the result of another query arrives, in
//!   the same batch or a later one, or in the epilog: blast_format.cpp:2362-2403), and its refusal
//!   comes there, after the reports before it and only the bytes NCBI has flushed by then: a
//!   refused first report leaves nothing, not even the document's start, which waits for the
//!   output's first write or flush as NCBI's waits in its stream's buffer (`blast_xml::StartFirst`;
//!   remaining-gaps DECISIONS D-SO-6). One rule for any batch layout: the refused report is held
//!   with its refusal (`Xml2Document::hold_refused`), carried to the next call where a batch ends
//!   with it, and where an error stops the run first the report and its refusal are dropped, as NCBI
//!   never makes the report (the off path's `write_before_failed_batch` drops them alike; the start
//!   that no call wrote is written there, `Stream::start_xml2_at_stop`). The refusals of the options and the
//!   query lengths that a later batch can meet come after the reports of the batches before in
//!   both paths (`after_batch_reports`, through `context_statistics::setup_failure`): scoring
//!   options without Karlin-Altschul values and matrix-only scoring with the greedy extension (in
//!   `context_statistics`, which the first batch with a valid query meets), and the ungapped word
//!   cutoff of a context of 2^30 letters or more or of an average subject length whose search
//!   space overflows (`batch_setup::batch_cutoff_set`), so they need no check here (line RG
//!   series D, after D1c-1c-b PR 2b). The refusals of a split batch's query chunks
//!   (`search_query_chunks`) and of the fetch under a negative seqid list come after the reports of
//!   the batches before in both paths (`after_batch_reports`, line RG item X-short-bat-2), so a
//!   negative `CHUNK_SIZE` streams. With one batch (the first batch takes every query record)
//!   nothing can come after it, so none of these checks is made and no subject value is made ahead
//!   (the one call makes them as without the switch). LOSAT's internal invariant errors and panics
//!   are not argued; a write error of the output stops the run where it happens (at the end without
//!   the switch); a `-db` file that changes or fails after the preliminary stage has read it is not
//!   argued either.
//! - what the report path writes to standard error (the query warnings) keeps its bytes and its
//!   place before each query's report; `-verbose` and the debug and trace output, which a batch
//!   writes while it is searched and the reports write again, keep the reports at the end
//!   (`Off::Debug`). `LOSAT_TIMING`'s blocks are printed when the loop takes a batch, as before,
//!   so with the switch a block comes after the warnings of the batches before it (the same
//!   lines in another order).
//!
//! Off (the reports after the loop, as without the switch): the switch off; outputs other
//! than one format of the command line's kind (the hit records and the format observer of the Web
//! adapter, whose HSP numbers count over the whole search; file sinks, which are created when
//! their report is written; two or more formats, which the library API allows on one stream,
//! where the off path writes each format for all queries and the batch writes would interleave
//! them); and `-verbose` and the debug output (`x_debug_output`). The command line (`-out` file or standard output) and ABI
//! v1 are streamed (an ABI caller gets the reports of the batches before a refusal, then its
//! error).
//!
//! The taxonomy database (DBID item T): NCBI prefetches the subjects that each batch reports
//! before the batch's reports in outfmt 0, 5, 6, 7 and 17, and with `taxdb.bti` on the search path
//! the taxonomy database writes its messages there; the batch loop adds them to the warnings of
//! the batch's first query (`run_in_pool`), so both paths write them before that query's report.
//! Their lookups can be refused (a taxonomy database or a taxid that LOSAT cannot read as NCBI
//! does), in the batch that first reports the sequence, and so can the run for NCBI's taxonomy
//! cache (`TaxdbRun::cache_collision`), in a batch that fetches a subject (a report's reads at
//! the first, a prefetch where NCBI's order of its loads decides the texts); either writes the
//! refusal after the reports of the batches before (above).
//!
//! `finish_sam` at NCBI's five error sites writes the document of the batches before as in the off
//! path. Memory: the records of the batches before are rendered text, and each batch's lists are
//! dropped, as for any streamed search.
//!
//! Memory: a batch's lists, flags, masks, Karlin-Altschul blocks and search spaces are dropped
//! after its write, and so are its queries' warnings. Kept for the whole search, as without the
//! switch: one (empty) warnings slot per query, and with `-db` the header of every distinct
//! sequence that a batch reported (`database::DatabaseTitles`: its stored title or its parsed
//! Seq-ids and titles, `ReportedHeader`; at most the database's sequences and at most the queries
//! times the hit list size), which the next batches read instead of reading the definition line
//! again.
//!
//! With `LOSAT_X_STATS` and the switch,
//! `[X_STATS] RPTSTREAM batches= queries= off= one_batch= check_s=` (batches: the batches whose
//! reports were written by the loop; off: why the switch did not act, `none` when it did;
//! `one_batch`: 1 when the checks were skipped for a search of one batch; `check_s`: the seconds of
//! `check`).

use std::ops::Range;
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::Result;

use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::seq_range::{next_ranged_batch, Placements, QueryInput};
use crate::report::blast_xml2::Xml2Carry;
use crate::report::sam::SamDocument;
use crate::utils::dust::MaskedInterval;

use super::{
    record_title, x_debug_output, BlastnArgs, BlastnOutputFormat, QueryBatching, SubjectLetters,
    SubjectRanges, SubjectSource,
};

// No NCBI counterpart: reads the switch once (LOSAT_X_RPTSTREAM is a strict-set switch since
// 2026-10-17: on unless it is `0`). It chooses when LOSAT writes the same reports.
/// Whether the switch is on: the reports may be written batch by batch, so a batch's query
/// records need not outlive the batch (LOSAT_X_QSTREAM reads them per batch only then,
/// `x_query_stream`).
pub(super) fn switch_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_RPTSTREAM"))
}

/// Why the reports of a search are written after its batch loop (`Stream::start`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The switch is off (`LOSAT_X_RPTSTREAM=0`).
    NotSet,
    /// Hit records, a format observer, a file sink, or two or more formats (`ReportOutputs` other
    /// than one format of the command line's kind).
    Outputs,
    /// `-verbose`, the debug variables, the blastn trace or the thread diagnostics.
    Debug,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Outputs => "outputs",
            Off::Debug => "debug",
        }
    }
}

/// The report values of a subject set's subjects that every batch's reports would make again
/// (`post_process_hits_and_write`), made once: the tabular ids of every subject (`None` where
/// LOSAT refuses one: each batch's report then makes them and refuses the subjects it reports,
/// as without the switch), and for outfmt 0 the titles and, with `-lcase_masking`, the lowercase
/// masks. None of them depends on the hit lists.
pub(super) struct SubjectCache {
    pub(super) ids: Option<Vec<crate::algorithm::common::tabular::SubjectId>>,
    pub(super) titles: Option<Vec<Arc<[u8]>>>,
    pub(super) masks: Option<Vec<Vec<MaskedInterval>>>,
}

/// What the report calls of a search wrote before a call (`write_reports`): the number of
/// queries whose reports were written, the subject values made once, the SAM document that
/// the calls add to (outfmt 17 writes its document at the end: `SamDocument`), the framing state
/// of the XML2 and JSON documents that the calls hand on (outfmt 13 to 16: `Xml2Carry`), and
/// whether a later call follows this one. The default is the off path's single call.
#[derive(Default)]
pub(super) struct Written<'a> {
    pub(super) queries: usize,
    pub(super) subjects: Option<&'a SubjectCache>,
    pub(super) sam: Option<&'a mut SamDocument>,
    pub(super) xml2: Option<&'a Mutex<Xml2Carry>>,
    /// A later call of the stream writes more of the run's reports (a batch's call): a call
    /// without the epilog does not end the run (`BlastnReportInputs::stops`).
    pub(super) continues: bool,
}

/// What `check` found where the switch acts: the subject values made once (a subject set of a
/// search of two or more batches), and whether the search is of one batch (no check made).
struct Ready {
    cache: Option<SubjectCache>,
    one_batch: bool,
}

/// The streaming of one search's reports: whether it acts (`off`) and what it has written. Its
/// drop prints the summary line.
pub(super) struct Stream {
    off: Option<Off>,
    cache: Option<SubjectCache>,
    /// The SAM document of the run (outfmt 17 joins the streaming: its records are held as
    /// rendered bytes until the last call, or an NCBI error, writes it).
    sam: Option<SamDocument>,
    /// The framing state of the run's XML2 or JSON document (outfmt 13 to 16), which each call
    /// takes up from the call before and hands on (`Xml2Carry`: the reports written, the start,
    /// the held report of a `-subject` run and its refusal).
    xml2: Option<Mutex<Xml2Carry>>,
    batches: u64,
    queries: usize,
    one_batch: bool,
    check_seconds: f64,
}

impl Stream {
    /// The streaming of a search with these outputs and subjects, `query_input` its query input,
    /// `batching` with the first batch's size.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn start(
        args: &BlastnArgs,
        outputs: &ReportOutputs<'_>,
        output_formats: &[BlastnOutputFormat],
        source: SubjectSource<'_>,
        subject_ranges: SubjectRanges<'_>,
        query_input: &QueryInput,
        batching: QueryBatching,
    ) -> Self {
        let started = std::time::Instant::now();
        let checked = check(
            switch_on(),
            args,
            outputs,
            output_formats,
            source,
            subject_ranges,
            query_input,
            batching,
        );
        let check_seconds = started.elapsed().as_secs_f64();
        let (off, cache, one_batch) = match checked {
            Ok(Ready { cache, one_batch }) => (None, cache, one_batch),
            Err(off) => (Some(off), None, false),
        };
        if switch_on() {
            crate::utils::xreach::x_reach!("LOSAT_X_RPTSTREAM", off.is_none());
        }
        // NCBI writes the document when its formatter is destroyed, whatever the number of
        // batches (`CBlast_SAM_Formatter`): the stream holds it for the calls.
        let sam = output_formats
            .iter()
            .find_map(|format| format.sam_spec())
            .filter(|_| off.is_none())
            .map(|spec| sam_document(args, spec));
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2332-2417
        // NCBI's formatter keeps the document's state for the run (`m_XMLFileCount`, the
        // accumulated results of a `-subject` run): the stream holds it between the calls.
        let xml2 = output_formats
            .iter()
            .any(|format| format.xml2_kind().is_some())
            .then(Mutex::default)
            .filter(|_| off.is_none());
        Stream {
            off,
            cache,
            sam,
            xml2,
            batches: 0,
            queries: 0,
            one_batch,
            check_seconds,
        }
    }

    /// Whether the batch loop writes each batch's reports.
    pub(super) fn acting(&self) -> bool {
        self.off.is_none()
    }

    /// Writes the reports of one batch of `queries` queries (or, with `epilog`, the epilogs) to
    /// `outputs` with `write`, which gets what earlier calls wrote.
    pub(super) fn write<F>(
        &mut self,
        outputs: &mut ReportOutputs<'_>,
        queries: usize,
        epilog: bool,
        write: F,
    ) -> Result<()>
    where
        F: FnOnce(&mut ReportOutputs<'_>, Written<'_>) -> Result<()>,
    {
        debug_assert!(self.acting());
        let written = Written {
            queries: self.queries,
            subjects: self.cache.as_ref(),
            sam: self.sam.as_mut(),
            xml2: self.xml2.as_ref(),
            continues: !epilog,
        };
        let result = write(outputs, written);
        if !epilog {
            self.batches += 1;
        }
        self.queries += queries;
        result
    }

    /// Writes the SAM document of the batches taken so far, if it holds a record: NCBI's
    /// formatter is destroyed when an error leaves the `try` block that holds it, so the
    /// records of the batches before the error are written before the error's lines
    /// (`blastn_app.cpp:186-327`). The off path does the same in its `write_before_failed_batch`
    /// call.
    pub(super) fn finish_sam(&mut self, outputs: &mut ReportOutputs<'_>) -> Result<()> {
        if !self.sam.as_ref().is_some_and(SamDocument::has_records) {
            return Ok(());
        }
        self.write(outputs, 0, true, |outputs, written| {
            let document = written.sam.expect("the stream holds a SAM document");
            // One format: more are never streamed (`Off::Outputs`).
            let mut writer = outputs.formats[0].sink.open()?;
            document.finish(&mut writer)?;
            writer.flush()?;
            Ok(())
        })
    }
}

impl Stream {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:258
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:360-367
    /// Writes the start of the run's `-subject` XML2 or JSON document where an error of NCBI's
    /// stops the run before any call wrote it (the calls held their reports: `blast_xml::StartFirst`):
    /// NCBI wrote it before the first batch (`PrintProlog`), and its normal exit flushes it. The
    /// off path's call writes it the same way (`write_before_failed_batch`).
    pub(super) fn start_xml2_at_stop(
        &mut self,
        outputs: &mut ReportOutputs<'_>,
        output_formats: &[BlastnOutputFormat],
    ) -> Result<()> {
        let (Some(carry), Some(kind)) = (
            self.xml2.as_ref(),
            output_formats.first().and_then(|format| format.xml2_kind()),
        ) else {
            return Ok(());
        };
        let mut start = Vec::new();
        carry
            .lock()
            .expect("the carry")
            .start_at_stop(kind, &mut start);
        if start.is_empty() {
            return Ok(());
        }
        // One format: more are never streamed (`Off::Outputs`).
        let mut writer = outputs.formats[0].sink.open()?;
        writer.write_all(&start)?;
        writer.flush()?;
        Ok(())
    }
}

// No NCBI counterpart: the summary line of the switch (with LOSAT_X_STATS); it prints counters and
// changes no value NCBI computes.
impl Drop for Stream {
    fn drop(&mut self) {
        if !switch_on() {
            return;
        }
        let off = self.off.map_or("none", Off::name);
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!(
                "[X_STATS] RPTSTREAM batches={} queries={} off={off} one_batch={} check_s={:.4}",
                self.batches,
                self.queries,
                u8::from(self.one_batch),
                self.check_seconds
            );
        }
    }
}

/// The SAM document of an outfmt 17 report with `spec`: `args` give the command line that its
/// `@PG` line shows (`BlastnArgs::command_line`).
pub(super) fn sam_document(args: &BlastnArgs, spec: crate::report::sam::SamSpec) -> SamDocument {
    SamDocument::new(
        spec,
        crate::wording::LOSAT_VERSION,
        args.command_line.clone(),
    )
}

/// The query input of the searched queries `searched` (a batch's), numbered from 0 as the batch's
/// lists number them: their placements (with `-query_loc`; the reports read only these:
/// `shift_to_records`, `pairwise_hits`, `blastn_pairwise_report`) and input positions, and the
/// range length. The input records' lengths are not the reports' (left empty).
pub(super) fn batch_input(query_input: &QueryInput, searched: Range<usize>) -> QueryInput {
    QueryInput {
        placements: if query_input.placements.ranged() {
            Placements(query_input.placements.0[searched.clone()].to_vec())
        } else {
            Placements::default()
        },
        ordinals: query_input.ordinals[searched].to_vec(),
        input_lengths: Vec::new(),
        skipped: Vec::new(),
        range_length: query_input.range_length,
    }
}

/// The reasons of `check` that need neither the query records nor the subjects: the switch off,
/// the outputs, and the verbose and debug output. Where one
/// applies the reports are written after the loop whatever the query file holds, so
/// LOSAT_X_QSTREAM reads the records before the search (its `off=reports`) instead of reading them
/// twice to keep them all for those reports. `on` is `switch_on()`.
pub(super) fn refused_before_reading(
    on: bool,
    args: &BlastnArgs,
    outputs: &ReportOutputs<'_>,
) -> Option<Off> {
    if !on {
        return Some(Off::NotSet);
    }
    if outputs.hits.is_some()
        || outputs.observer.is_some()
        || outputs.formats.len() > 1
        || outputs
            .formats
            .iter()
            .any(|format| matches!(format.sink, OutputSink::File(_)))
    {
        return Some(Off::Outputs);
    }
    if x_debug_output(args) {
        return Some(Off::Debug);
    }
    None
}

/// Why the switch does not act for this search, or what it needs where it does (`Stream::start`;
/// the module doc has the reasons). `on` is `switch_on()` (a parameter for the unit tests).
#[allow(clippy::too_many_arguments)]
fn check(
    on: bool,
    args: &BlastnArgs,
    outputs: &ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    source: SubjectSource<'_>,
    subject_ranges: SubjectRanges<'_>,
    query_input: &QueryInput,
    batching: QueryBatching,
) -> std::result::Result<Ready, Off> {
    if let Some(off) = refused_before_reading(on, args, outputs) {
        return Err(off);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-170
    // The first batch takes every record (`next_ranged_batch` with the first batch's size, as
    // the loop reads it): no refusal can come after it (a refusal of that batch finds what its one
    // call writes before it in both paths: nothing, or for a `-subject` XML2 or JSON document the
    // reports before the refused one, `blast_xml::write_blastn_xml`; what ends the reading after the
    // records comes in that batch or in a batch of its own, NCBI's errors, which both paths write
    // alike), so nothing is checked and the one call makes the subject values as without the
    // switch.
    let (first_end, _) = next_ranged_batch(
        &query_input.input_lengths,
        &query_input.skipped,
        0,
        batching.batch_size,
    );
    if first_end == query_input.input_lengths.len() {
        return Ok(Ready {
            cache: None,
            one_batch: true,
        });
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:310-312
    // The stored data and the report ids of the reported subjects (a database's definition
    // lines, ambiguity data and taxonomy, a subject set's ids) are checked by the report call of
    // the batch that reports them, as NCBI's formatter reads them for each batch's result sets
    // (drop-in DECISIONS 44): a refusal there comes after the reports of the batches before (the
    // module doc).
    let cache = match source {
        SubjectSource::Prepared(prepared) => {
            // `tabular_subject_ids` of every subject, where every lookup finds its id; only the
            // tabular reports read them (`post_process_hits_and_write`), so none is made for the
            // others.
            let tabular = output_formats.iter().any(|format| format.is_tabular());
            let ids: Option<Vec<crate::algorithm::common::tabular::SubjectId>> = tabular
                .then(|| {
                    prepared
                        .records
                        .iter()
                        .map(|subject| {
                            subject
                                .subject_tabular_id(false)
                                .ok()
                                .flatten()
                                .map(|id| Some(Arc::from(id)))
                        })
                        .collect()
                })
                .flatten();
            let pairwise = output_formats.iter().any(|format| format.is_pairwise());
            let titles = pairwise.then(|| prepared.records.iter().map(record_title).collect());
            // The records as read with `-subject_loc` (never packed), else the searched ones.
            let masks = (pairwise && args.lcase_masking).then(|| {
                let subject_records = match subject_ranges.records {
                    Some(records) => SubjectLetters::new(records, None),
                    None => prepared.letters(),
                };
                (0..subject_records.records().len())
                    .map(|index| subject_records.lowercase_masks(index))
                    .collect()
            });
            Some(SubjectCache { ids, titles, masks })
        }
        SubjectSource::Database(_) => None,
    };
    Ok(Ready {
        cache,
        one_batch: false,
    })
}

#[cfg(test)]
#[path = "x_rpt_stream_tests.rs"]
mod tests;
