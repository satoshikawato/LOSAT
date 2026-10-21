//! EXPERIMENT (LOSAT_X_RPTSTREAM, shadow LOSAT_X_RPTSTREAMSHADOW): blastn writes the reports of
//! each query batch when the batch loop of `run_in_pool` takes the batch, and drops its final hit
//! lists, as NCBI prints each batch's result sets in its batch loop; without the switch the loop
//! keeps the final lists of every batch and writes all reports after the loop.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-318
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!     ...
//!                 results = lcl_blast.Run();
//!                 if (!batch_size)
//!                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
//!     ...
//!                 ITERATE(CSearchResultSet, result, *results) {
//!                     formatter.PrintOneResultSet(**result, query_batch);
//!                 }
//!     ...
//!         }
//!         BLAST_PROF_START( APP.POST );
//!         formatter.PrintEpilog(opt);
//! ```
//!
//! With the switch on, the loop writes a batch's reports right after it has taken the batch (in
//! batch order: a LOSAT_X_QBATCHPAR slot's result is taken under the loop's own key, in order,
//! as before), finished its lists (`-db`: `finish_query_hit_list`, the titles, the HSPs of score
//! 0) and given its gie to the mixer and the slots. It calls the same report code
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
//!   step of `finish_query_hit_list`.
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
//!   formats 8, 9 and 12 (OF3) write per query like 6; BLAST XML (5, OF4a) is not argued here
//!   and is never streamed (`Off::Xml`, below).
//! - the errors: an error of NCBI's that a batch meets (a reader error, an `eEOF` batch, a batch
//!   of queries without letters, `ErrorAfterInvalidBatches`, a header file that cannot be read)
//!   comes after the reports of the batches before, without the epilog, in both: here they were
//!   written when the loop took them, so the error paths write only what follows them (the
//!   batch's title and reader warnings). An error that LOSAT raises where NCBI would not stop
//!   (LOSAT's refusals) is written without the reports of the batches before; so with two or more
//!   batches the switch acts only where no such refusal can come after the first batch (`check`).
//!   BLAST XML has refusals of its own (a `-db` name, query title or id, or a title of a subject
//!   with a hit that cannot be written, `blast_xml.rs`), made per report call, so `-outfmt 5` is
//!   never streamed (`Off::Xml`, any number of batches): the report runs once after the search and
//!   every refusal comes before the first byte. Without XML the conditions are:
//!   scoring options with Karlin-Altschul values and greedy gap costs that LOSAT supports (the
//!   refusals of `search_query_batch` for scoring options without values or large greedy gap
//!   costs), no negative `CHUNK_SIZE` (the refusal in `search_query_chunks`), search spaces that
//!   cannot overflow (`SearchSpace`), subjects whose report ids are all supported (a subject set's
//!   `tabular_subject_ids`, made once here), and with `-db` stored data that the traceback and the
//!   reports read without a refusal: every definition line (as blastp's block search does) and
//!   the ambiguity data of every sequence, which the preliminary stage does not read
//!   (`database::x_streaming_refusal`; one pass over the database's sequences before the search,
//!   in the pool). With one batch (the first batch takes every query record) nothing can come
//!   after it, so none of these checks is made and no subject value is made ahead (the one call
//!   makes them as without the switch). LOSAT's internal invariant errors and panics are not
//!   argued; a write error of the output stops the run where it happens (at the end without the
//!   switch); a `-db` file that changes or fails after the preliminary stage has read it is not
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
//! them); `-verbose` and the debug output (`x_debug_output`, the list LOSAT_X_QBATCHPAR uses); and
//! the refusals above. The command line (`-out` file or standard output) and ABI v1 are streamed.
//!
//! Memory: a batch's lists, flags, masks, Karlin-Altschul blocks and search spaces are dropped
//! after its write, and so are its queries' warnings. Kept for the whole search, as without the
//! switch: one (empty) warnings slot per query, and with `-db` the stored title of every distinct
//! sequence that a batch reported (`database::DatabaseTitles`, at most the database's sequences
//! and at most the queries times the hit list size), which the next batches read instead of
//! reading the definition line again.
//!
//! LOSAT_X_QBATCHPAR: the writes happen while the slots search the batches ahead (they write
//! nothing; the batch loop alone writes), which is the point for the many-batch searches. The
//! error paths still stop every slot before they write. A slot's panic, which ends a release
//! build at once, then ends it after the reports of the batches taken so far, possibly in the
//! middle of one (what its buffered writer had flushed).
//!
//! The shadow (`LOSAT_X_RPTSTREAMSHADOW`, implies the switch) keeps the final lists and values of
//! every batch as the off path does, copies every byte that the batch writes send to each format
//! and to standard error, and at the end (and on the error paths that write the reports of the
//! batches before, and on LOSAT's errors of a batch, where the off path writes none) writes the
//! reports the off path's way into memory and compares the bytes: on a difference it prints
//! `[X_SHADOW] RPTSTREAM ... DIFFERENT` and panics; at the end of the search it prints
//! `[X_SHADOW] RPTSTREAM batches= bytes= checks= off= identical` (batches: the batches whose
//! reports were written by the loop; bytes: the bytes compared; checks: the comparisons; off: why
//! the switch did not act, `none` when it did). With `LOSAT_X_STATS` and the switch,
//! `[X_STATS] RPTSTREAM batches= queries= off= one_batch= check_s=` (`one_batch`: 1 when the
//! checks were skipped for a search of one batch; `check_s`: the seconds of `check`, the `-db`
//! pass included).

use std::io::{self, Write};
use std::ops::Range;
use std::sync::{Arc, OnceLock};

use anyhow::Result;

use crate::api::local_blast::{FormatOutput, OutputSink, ReportOutputs};
use crate::blastinput::seq_range::{next_ranged_batch, Placements, QueryInput};
use crate::config::NuclScoringSpec;
use crate::stats::NuclValues;
use crate::utils::dust::MaskedInterval;
use crate::utils::threading::SearchPool;

use super::{
    check_greedy_gap_costs, configure_task, database, record_title, x_debug_output, BlastnArgs,
    BlastnOutputFormat, QueryBatching, SubjectLetters, SubjectRanges, SubjectSource,
};

// No NCBI counterpart: reads the switches once (LOSAT_X_RPTSTREAM is a strict-set switch since
// 2026-10-17: on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set; the shadow, off
// unless set to any value, implies the switch). They choose when LOSAT writes the same reports.
/// 0 = off, 1 = LOSAT_X_RPTSTREAM, 2 = LOSAT_X_RPTSTREAMSHADOW.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_RPTSTREAMSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_RPTSTREAM") {
            1
        } else {
            0
        }
    })
}

/// Whether the switch (or its shadow) is set: the reports may be written batch by batch, so a
/// batch's query records need not outlive the batch (LOSAT_X_QSTREAM reads them per batch only
/// then, `x_query_stream`).
pub(super) fn switch_on() -> bool {
    mode() != 0
}

/// Why the reports of a search are written after its batch loop (`Stream::start`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// Neither switch is set.
    NotSet,
    /// Hit records, a format observer, a file sink, or two or more formats (`ReportOutputs` other
    /// than one format of the command line's kind).
    Outputs,
    /// `-verbose`, the debug variables, the blastn trace or the thread diagnostics.
    Debug,
    /// Scoring options without Karlin-Altschul values, or greedy gap costs that LOSAT refuses.
    Scoring,
    /// `-outfmt 5` (BLAST XML), which is written after the loop (the module doc has the reason).
    Xml,
    /// A negative `CHUNK_SIZE`.
    ChunkSize,
    /// A database length times a query length that could overflow NCBI's 64-bit search space.
    SearchSpace,
    /// A subject of the set whose tabular id LOSAT refuses.
    SubjectIds,
    /// A database whose definition lines do not all read as the reports read them.
    Titles,
    /// A database sequence whose ambiguity data the traceback cannot decode.
    Ambiguity,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Outputs => "outputs",
            Off::Debug => "debug",
            Off::Xml => "xml",
            Off::Scoring => "scoring",
            Off::ChunkSize => "chunk_size",
            Off::SearchSpace => "search_space",
            Off::SubjectIds => "subject_ids",
            Off::Titles => "titles",
            Off::Ambiguity => "ambiguity",
        }
    }
}

/// The report values of a subject set's subjects that every batch's reports would make again
/// (`post_process_hits_and_write`), made once: the tabular ids of every subject (all supported:
/// `Off::SubjectIds`), and for outfmt 0 the titles and, with `-lcase_masking`, the lowercase
/// masks. None of them depends on the hit lists.
pub(super) struct SubjectCache {
    pub(super) ids: Vec<Arc<[u8]>>,
    pub(super) titles: Option<Vec<Arc<[u8]>>>,
    pub(super) masks: Option<Vec<Vec<MaskedInterval>>>,
}

/// What the report calls of a search wrote before a call (`write_reports`): the number of
/// queries whose reports were written, and the subject values made once. The default is the
/// off path's single call.
#[derive(Clone, Copy, Default)]
pub(super) struct Written<'a> {
    pub(super) queries: usize,
    pub(super) subjects: Option<&'a SubjectCache>,
}

/// The bytes that the batch writes sent to each format and to standard error, and the
/// comparisons made (the shadow).
struct Shadow {
    outfmts: Vec<String>,
    formats: Vec<Vec<u8>>,
    diagnostics: Vec<u8>,
    checks: u64,
    compared: u64,
}

/// What `check` found where the switch acts: the subject values made once (a subject set of a
/// search of two or more batches), and whether the search is of one batch (no check made).
struct Ready {
    cache: Option<SubjectCache>,
    one_batch: bool,
}

/// The streaming of one search's reports: whether it acts (`off`), what it has written, and the
/// shadow's copies. Its drop prints the summary lines.
pub(super) struct Stream {
    off: Option<Off>,
    cache: Option<SubjectCache>,
    shadow: Option<Shadow>,
    batches: u64,
    queries: usize,
    one_batch: bool,
    check_seconds: f64,
}

impl Stream {
    /// The streaming of a search with these outputs and subjects, `query_input` and `lengths` its
    /// query input and the searched queries' lengths, `batching` with the first batch's size
    /// (`run_in_pool`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn start(
        args: &BlastnArgs,
        outputs: &ReportOutputs<'_>,
        output_formats: &[BlastnOutputFormat],
        source: SubjectSource<'_>,
        subject_ranges: SubjectRanges<'_>,
        query_input: &QueryInput,
        lengths: &[usize],
        batching: QueryBatching,
        pool: &SearchPool<'_>,
    ) -> Self {
        let started = std::time::Instant::now();
        let checked = check(
            mode(),
            args,
            outputs,
            output_formats,
            source,
            subject_ranges,
            query_input,
            lengths,
            batching,
            pool,
        );
        let check_seconds = started.elapsed().as_secs_f64();
        let (off, cache, one_batch) = match checked {
            Ok(Ready { cache, one_batch }) => (None, cache, one_batch),
            Err(off) => (Some(off), None, false),
        };
        let shadow = (mode() == 2 && off.is_none()).then(|| Shadow {
            outfmts: outputs
                .formats
                .iter()
                .map(|format| format.outfmt.to_string())
                .collect(),
            formats: vec![Vec::new(); outputs.formats.len()],
            diagnostics: Vec::new(),
            checks: 0,
            compared: 0,
        });
        Stream {
            off,
            cache,
            shadow,
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

    /// Whether the off path's values are kept for the shadow's comparison.
    pub(super) fn shadow(&self) -> bool {
        self.shadow.is_some()
    }

    /// Writes the reports of one batch of `queries` queries (or, with `epilog`, the epilogs) to
    /// `outputs` with `write`, which gets what earlier calls wrote. The shadow copies the bytes.
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
        };
        let result = match self.shadow.as_mut() {
            None => write(outputs, written),
            Some(shadow) => {
                let outfmts: Vec<&str> =
                    outputs.formats.iter().map(|format| format.outfmt).collect();
                let mut tees: Vec<SinkTee<'_, '_>> = outputs
                    .formats
                    .iter_mut()
                    .zip(shadow.formats.iter_mut())
                    .map(|(format, copy)| SinkTee {
                        sink: &mut format.sink,
                        copy,
                    })
                    .collect();
                let mut diagnostics = WriteTee {
                    inner: &mut *outputs.diagnostics,
                    copy: &mut shadow.diagnostics,
                };
                let mut teed = ReportOutputs {
                    formats: tees
                        .iter_mut()
                        .zip(outfmts)
                        .map(|(tee, outfmt)| FormatOutput {
                            outfmt,
                            sink: OutputSink::Writer(tee),
                        })
                        .collect(),
                    diagnostics: &mut diagnostics,
                    hits: None,
                    observer: None,
                };
                write(&mut teed, written)
            }
        };
        if !epilog {
            self.batches += 1;
        }
        self.queries += queries;
        result
    }

    /// The shadow: writes with `write_old` (the off path's call) into memory and compares the
    /// bytes with those that the batch writes sent (`what` names the place in a difference).
    /// Without the shadow nothing runs.
    pub(super) fn check<F>(&mut self, what: &str, write_old: F)
    where
        F: FnOnce(&mut ReportOutputs<'_>) -> Result<()>,
    {
        let Some(shadow) = self.shadow.as_mut() else {
            return;
        };
        let mut formats: Vec<Vec<u8>> = vec![Vec::new(); shadow.outfmts.len()];
        let mut diagnostics = Vec::new();
        let result = {
            let mut memory = ReportOutputs {
                formats: shadow
                    .outfmts
                    .iter()
                    .zip(formats.iter_mut())
                    .map(|(outfmt, bytes)| FormatOutput {
                        outfmt,
                        sink: OutputSink::Writer(bytes),
                    })
                    .collect(),
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            write_old(&mut memory)
        };
        shadow.checks += 1;
        let difference = match result {
            Err(error) => Some(format!(
                "the reports written after the loop fail: {error:#}"
            )),
            Ok(()) => formats
                .iter()
                .zip(&shadow.formats)
                .enumerate()
                .find_map(|(index, (old, streamed))| {
                    first_difference(old, streamed).map(|at| {
                        format!(
                            "format {index} ({}) differs at byte {at} of {} / {}",
                            shadow.outfmts[index],
                            old.len(),
                            streamed.len()
                        )
                    })
                })
                .or_else(|| {
                    first_difference(&diagnostics, &shadow.diagnostics).map(|at| {
                        format!(
                            "the standard error differs at byte {at} of {} / {}",
                            diagnostics.len(),
                            shadow.diagnostics.len()
                        )
                    })
                }),
        };
        if let Some(difference) = difference {
            eprintln!(
                "[X_SHADOW] RPTSTREAM batches={} checks={} DIFFERENT ({what}: {difference})",
                self.batches, shadow.checks
            );
            panic!("LOSAT_X_RPTSTREAMSHADOW: the reports written batch by batch differ from those written after the loop");
        }
        shadow.compared +=
            formats.iter().map(|bytes| bytes.len() as u64).sum::<u64>() + diagnostics.len() as u64;
    }
}

// No NCBI counterpart: the summary lines of the switch (with LOSAT_X_STATS) and of the shadow;
// they print counters and change no value NCBI computes.
impl Drop for Stream {
    fn drop(&mut self) {
        if mode() == 0 {
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
        if mode() == 2 {
            let (checks, compared) = self
                .shadow
                .as_ref()
                .map_or((0, 0), |shadow| (shadow.checks, shadow.compared));
            eprintln!(
                "[X_SHADOW] RPTSTREAM batches={} bytes={compared} checks={checks} off={off} identical",
                self.batches
            );
        }
    }
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

/// The first byte at which `a` and `b` differ (or the shorter one's length), `None` if equal.
fn first_difference(a: &[u8], b: &[u8]) -> Option<usize> {
    if a == b {
        return None;
    }
    Some(
        a.iter()
            .zip(b)
            .position(|(left, right)| left != right)
            .unwrap_or(a.len().min(b.len())),
    )
}

/// A format's sink during a batch write of the shadow: the bytes go to the sink and are copied.
/// A sink of the command line or of ABI v1 is a writer (the file sink is not streamed, `check`);
/// `Stdout` is written unbuffered, the same bytes in the same order.
struct SinkTee<'s, 'r> {
    sink: &'s mut OutputSink<'r>,
    copy: &'s mut Vec<u8>,
}

impl Write for SinkTee<'_, '_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = match self.sink {
            OutputSink::Writer(writer) => writer.write(buf)?,
            OutputSink::Stdout => io::stdout().write(buf)?,
            OutputSink::File(_) => {
                return Err(io::Error::other(
                    "LOSAT_X_RPTSTREAM does not stream a file sink",
                ))
            }
        };
        self.copy.extend_from_slice(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        match self.sink {
            OutputSink::Writer(writer) => writer.flush(),
            OutputSink::Stdout => io::stdout().flush(),
            OutputSink::File(_) => Ok(()),
        }
    }
}

/// Standard error during a batch write of the shadow: the bytes go to it and are copied.
struct WriteTee<'s> {
    inner: &'s mut (dyn Write + Send),
    copy: &'s mut Vec<u8>,
}

impl Write for WriteTee<'_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.copy.extend_from_slice(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Why the switch does not act for this search, or what it needs where it does (`Stream::start`;
/// the module doc has the reasons). `mode` is `mode()` (a parameter for the unit tests).
#[allow(clippy::too_many_arguments)]
fn check(
    mode: u8,
    args: &BlastnArgs,
    outputs: &ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    source: SubjectSource<'_>,
    subject_ranges: SubjectRanges<'_>,
    query_input: &QueryInput,
    lengths: &[usize],
    batching: QueryBatching,
    pool: &SearchPool<'_>,
) -> std::result::Result<Ready, Off> {
    if mode == 0 {
        return Err(Off::NotSet);
    }
    if outputs.hits.is_some()
        || outputs.observer.is_some()
        || outputs.formats.len() > 1
        || outputs
            .formats
            .iter()
            .any(|format| matches!(format.sink, OutputSink::File(_)))
    {
        return Err(Off::Outputs);
    }
    if x_debug_output(args) {
        return Err(Off::Debug);
    }
    // BLAST XML (DESIGN_OF4 Q8, Owner-delegated 2026-10-20): the XML writer refuses a `-db` name,
    // a query title or id, or the title of a subject with a hit that it cannot write as NCBI's
    // XML (DESIGN_OF4 Q2), and those checks are made by the report call of the batch that holds
    // the query or the subject. With several batches a later batch's refusal would come after
    // the first batch's Iterations were written, so outfmt 5 is never streamed: the one report
    // call after the search finds every refusal before the first byte.
    if output_formats.contains(&BlastnOutputFormat::Xml) {
        return Err(Off::Xml);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-170
    // ```c++
    //     TSeqPos size_read = 0;
    //
    //     while (size_read < GetBatchSize()) {
    //
    //         if (End())
    //             break;
    // ```
    // The first batch takes every record (`next_ranged_batch` with the first batch's size, as
    // the loop reads it): no refusal can come after it (a refusal of that batch finds nothing
    // written; what ends the reading after the records comes in that batch or in a batch of its
    // own, NCBI's errors, which both paths write alike), so nothing is checked and the one call
    // makes the subject values as without the switch.
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
    // The scoring of `search_query_batch` (`context_blocks`: `NuclValues::gapped` fails exactly
    // where `check_gaps` does, for every context) and its greedy gap costs: neither depends on
    // the queries.
    let config = configure_task(args);
    let spec = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    if NuclValues::new(spec.reward, spec.penalty)
        .and_then(|values| values.check_gaps(&spec))
        .is_err()
        || check_greedy_gap_costs(args).is_err()
    {
        return Err(Off::Scoring);
    }
    if batching.split.negative_chunk_size() {
        return Err(Off::ChunkSize);
    }
    // The search space of a context (`search_query_batch`): the effective database length
    // `db - n x a` (at least 1) times the effective query length `q - a` (at least 1), with the
    // length adjustment `a` between 0 and the smaller root of `n a^2 - (q n + db) a + c`, which is
    // at most `2 c / (q n + db) <= 2 q` (`compute_length_adjustment_ncbi`, `c <= q db`). So
    // `(db + n (2 q + 1)) (2 q + 1)` bounds every product that can overflow there.
    let metadata = source.metadata();
    let db_length = if args.db_length() > 0 {
        i128::from(args.db_length())
    } else {
        metadata.db_len_total as i128
    };
    let factor = 2 * lengths.iter().copied().max().unwrap_or(0) as i128 + 1;
    if (db_length.abs() + metadata.db_num_seqs as i128 * factor) * factor > i128::from(i64::MAX) {
        return Err(Off::SearchSpace);
    }
    match source {
        SubjectSource::Prepared(prepared) => {
            // `tabular_subject_ids` of every subject.
            let ids: Option<Vec<Arc<[u8]>>> = prepared
                .records
                .iter()
                .map(|subject| {
                    crate::report::defline::tabular_subject_id(
                        subject.local_id.as_bytes(),
                        &subject.title,
                        false,
                    )
                    .ok()
                    .map(Arc::from)
                })
                .collect();
            let ids = ids.ok_or(Off::SubjectIds)?;
            let pairwise = output_formats.contains(&BlastnOutputFormat::Pairwise);
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
            Ok(Ready {
                cache: Some(SubjectCache { ids, titles, masks }),
                one_batch: false,
            })
        }
        SubjectSource::Database(db) => match database::x_streaming_refusal(db, pool) {
            Some(off) => Err(off),
            None => Ok(Ready {
                cache: None,
                one_batch: false,
            }),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::seq_range::RecordPlacement;

    // EXPERIMENT (LOSAT_X_RPTSTREAM): a batch's part of the query input numbers its queries from
    // 0, as the batch's lists do; without `-query_loc` it has no placements either.
    #[test]
    fn a_batch_input_is_the_batch_part_numbered_from_zero() {
        let placements = (0..5)
            .map(|index| RecordPlacement {
                offset: 10 * index,
                length: 100 + index,
            })
            .collect();
        let input = QueryInput {
            placements: Placements(placements),
            ordinals: vec![0, 1, 3, 4, 6],
            input_lengths: vec![100; 7],
            skipped: vec![false, false, true, false, false, true, false],
            range_length: Some(50),
        };
        let batch = batch_input(&input, 2..4);
        assert!(batch.placements.ranged());
        assert_eq!(
            (batch.placements.offset(0), batch.placements.length(1, 7)),
            (20, 103)
        );
        assert_eq!(batch.placements.0.len(), 2);
        assert_eq!(batch.ordinals, [3, 4]);
        assert_eq!(batch.range_length, Some(50));
        let whole = QueryInput {
            placements: Placements::default(),
            ordinals: vec![0, 1, 2],
            input_lengths: vec![10; 3],
            skipped: vec![false; 3],
            range_length: None,
        };
        let batch = batch_input(&whole, 1..3);
        assert!(!batch.placements.ranged());
        assert_eq!(batch.placements.offset(1), 0);
        assert_eq!(batch.ordinals, [1, 2]);
    }

    // EXPERIMENT (LOSAT_X_RPTSTREAMSHADOW): the shadow's tees send every byte on and copy it.
    #[test]
    fn the_tees_send_and_copy_every_byte() {
        let mut target: Vec<u8> = Vec::new();
        let mut copy = Vec::new();
        {
            let mut sink = OutputSink::Writer(&mut target);
            let mut tee = SinkTee {
                sink: &mut sink,
                copy: &mut copy,
            };
            tee.write_all(b"row 1\n").unwrap();
            write!(tee, "row {}\n", 2).unwrap();
            tee.flush().unwrap();
        }
        assert_eq!(target, b"row 1\nrow 2\n");
        assert_eq!(copy, target);
        let mut errors: Vec<u8> = Vec::new();
        let mut copy = Vec::new();
        {
            let mut tee = WriteTee {
                inner: &mut errors,
                copy: &mut copy,
            };
            tee.write_all(b"Warning: one\n").unwrap();
        }
        assert_eq!(
            (errors.as_slice(), copy.as_slice()),
            (&b"Warning: one\n"[..], &b"Warning: one\n"[..])
        );
        let path = std::path::Path::new("never-created");
        let mut sink = OutputSink::File(path);
        let mut copy = Vec::new();
        let mut tee = SinkTee {
            sink: &mut sink,
            copy: &mut copy,
        };
        assert!(tee.write_all(b"x").is_err(), "a file sink is not streamed");
        assert!(!path.exists());
    }

    /// The arguments of a `blastn` command line (`extra` after `-query q -subject s`).
    fn args(extra: &[&str]) -> BlastnArgs {
        use clap::Parser;
        #[derive(Parser)]
        struct Options {
            #[command(flatten)]
            blastn: BlastnArgs,
        }
        let mut argv = vec!["test", "-query", "q", "-subject", "s"];
        argv.extend_from_slice(extra);
        crate::cli::try_parse_from::<Options, _, _>(argv)
            .unwrap()
            .blastn
    }

    fn record(
        index: usize,
        title: &[u8],
        letters: usize,
    ) -> crate::blastinput::fasta_reader::FastaRecord {
        crate::blastinput::fasta_reader::FastaRecord {
            local_id: format!("Subject_{}", index + 1),
            title: title.to_vec(),
            sequence: b"ACGTTGCA".iter().copied().cycle().take(letters).collect(),
            warnings: Vec::new(),
        }
    }

    /// `check` for queries of 1,000 letters each (`queries`), first batches of `batch_size`
    /// letters, the subjects `subjects` and the outputs `outputs`.
    fn checked(
        mode: u8,
        extra: &[&str],
        subjects: &[crate::blastinput::fasta_reader::FastaRecord],
        queries: usize,
        batch_size: u32,
        split: crate::algorithm::blastn::query_split::SplitSizes,
        outputs: &mut ReportOutputs<'_>,
    ) -> std::result::Result<Ready, Off> {
        checked_as(
            &[BlastnOutputFormat::Tabular],
            mode,
            extra,
            subjects,
            queries,
            batch_size,
            split,
            outputs,
        )
    }

    /// `checked` for the output formats `formats`.
    #[allow(clippy::too_many_arguments)]
    fn checked_as(
        formats: &[BlastnOutputFormat],
        mode: u8,
        extra: &[&str],
        subjects: &[crate::blastinput::fasta_reader::FastaRecord],
        queries: usize,
        batch_size: u32,
        split: crate::algorithm::blastn::query_split::SplitSizes,
        outputs: &mut ReportOutputs<'_>,
    ) -> std::result::Result<Ready, Off> {
        let args = args(extra);
        let query_records: Vec<_> = (0..queries)
            .map(|index| record(index, b"q", 1000))
            .collect();
        let query_input = QueryInput::whole(&query_records);
        let lengths: Vec<usize> = query_records.iter().map(|q| q.sequence.len()).collect();
        let prepared = super::super::PreparedSubjects {
            records: subjects,
            stores: None,
            metadata: super::super::subject_metadata_from_records(subjects, |index| {
                subjects[index].sequence.len()
            }),
            packed: Vec::new(),
        };
        let batching = QueryBatching {
            fixed_batch_size: batch_size as i32,
            batch_size,
            split,
        };
        // A one-thread pool (the `-db` pass is not reached for a subject set).
        let pool = SearchPool::x_serial();
        check(
            mode,
            &args,
            &*outputs,
            formats,
            SubjectSource::Prepared(&prepared),
            SubjectRanges {
                placements: &Placements::default(),
                records: None,
            },
            &query_input,
            &lengths,
            batching,
            &pool,
        )
    }

    // EXPERIMENT (LOSAT_X_RPTSTREAM): where the switch acts, and why not elsewhere: one format of
    // the command line's kind only; with one batch no refusal is checked (and nothing is made
    // ahead); with more, every refusal that could come after the first batch is.
    #[test]
    fn the_switch_acts_only_where_no_refusal_can_follow_the_first_batch() {
        use crate::algorithm::blastn::query_split::SplitSizes;
        let plain = SplitSizes::new(false, None, None);
        let good = vec![
            record(0, b"s1 a subject", 4000),
            record(1, b"s2 another", 3000),
        ];
        let refused = vec![
            record(0, b"s1 a subject", 4000),
            record(1, b"Subject_\x81 x", 3000),
        ];
        let (mut one, mut two) = (Vec::<u8>::new(), Vec::<u8>::new());
        let mut errors = Vec::<u8>::new();
        let single =
            |sink: &mut Vec<u8>, errors: &mut Vec<u8>| -> std::result::Result<Ready, Off> {
                let mut outputs = ReportOutputs::single("6", OutputSink::Writer(sink), errors);
                checked(1, &[], &good, 10, 4000, plain, &mut outputs)
            };
        // Not set.
        let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors);
        assert_eq!(
            checked(0, &[], &good, 10, 4000, plain, &mut outputs).err(),
            Some(Off::NotSet)
        );
        drop(outputs);
        // Hit records, an observer-like second output: two formats; a file sink.
        let mut on_hits = |_: &[crate::report::PairwiseHit]| {};
        let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors);
        outputs.hits = Some(&mut on_hits);
        assert_eq!(
            checked(1, &[], &good, 10, 4000, plain, &mut outputs).err(),
            Some(Off::Outputs)
        );
        drop(outputs);
        let mut errors2 = Vec::<u8>::new();
        let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors2);
        outputs.formats.push(FormatOutput {
            outfmt: "7",
            sink: OutputSink::Writer(&mut two),
        });
        assert_eq!(
            checked(1, &[], &good, 10, 4000, plain, &mut outputs).err(),
            Some(Off::Outputs)
        );
        drop(outputs);
        let path = std::path::Path::new("never-created");
        let mut outputs = ReportOutputs::single("6", OutputSink::File(path), &mut errors);
        assert_eq!(
            checked(1, &[], &good, 10, 4000, plain, &mut outputs).err(),
            Some(Off::Outputs)
        );
        drop(outputs);
        // Several batches (4,000 letters of 10,000): the subject values are made once.
        let ready = single(&mut one, &mut errors).expect("acts");
        assert!(!ready.one_batch);
        assert_eq!(ready.cache.as_ref().map(|cache| cache.ids.len()), Some(2));
        // One batch (20,000 letters): no check, nothing made ahead, even with a refused id.
        let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors);
        let ready = checked(1, &[], &refused, 10, 20_000, plain, &mut outputs).expect("acts");
        assert!(ready.one_batch && ready.cache.is_none());
        // Several batches: the refusals that could come after the first batch.
        assert_eq!(
            checked(1, &[], &refused, 10, 4000, plain, &mut outputs).err(),
            Some(Off::SubjectIds)
        );
        let scoring = [
            "-task",
            "blastn",
            "-reward",
            "2",
            "-penalty",
            "-3",
            "-gapopen",
            "1",
            "-gapextend",
            "1",
        ];
        assert_eq!(
            checked(1, &scoring, &good, 10, 4000, plain, &mut outputs).err(),
            Some(Off::Scoring)
        );
        assert_eq!(
            checked(
                1,
                &["-dbsize", "9000000000000000000"],
                &good,
                10,
                4000,
                plain,
                &mut outputs
            )
            .err(),
            Some(Off::SearchSpace)
        );
        assert_eq!(
            checked(
                1,
                &[],
                &good,
                10,
                4000,
                SplitSizes::new(false, Some(-5), Some(-10)),
                &mut outputs
            )
            .err(),
            Some(Off::ChunkSize)
        );
        // The shadow mode checks the same.
        assert!(checked(2, &[], &good, 10, 4000, plain, &mut outputs).is_ok());
        // BLAST XML is never streamed, with one batch (20,000 letters) or several.
        for batch_size in [20_000, 4000] {
            let xml = [BlastnOutputFormat::Xml];
            assert_eq!(
                checked_as(&xml, 1, &[], &good, 10, batch_size, plain, &mut outputs).err(),
                Some(Off::Xml)
            );
        }
    }

    #[test]
    fn the_first_difference_is_the_first_unequal_byte_or_the_shorter_length() {
        assert_eq!(first_difference(b"abc", b"abc"), None);
        assert_eq!(first_difference(b"abc", b"abd"), Some(2));
        assert_eq!(first_difference(b"ab", b"abc"), Some(2));
        assert_eq!(first_difference(b"", b"a"), Some(0));
    }
}
