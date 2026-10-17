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
//!   lists are finished per query (`finish_query_hit_list`: `-max_hsps`, `-subject_besthit`, the
//!   sort, `-max_target_seqs`; the HSPs of score 0), the coordinates move by the query's and the
//!   subject's own range start (`shift_to_records`), and the rows, outfmt 7 headers and outfmt 0
//!   sections are written query by query from the query's list, record, masks, blocks and search
//!   space, the subject records and the options. BLASTN has no culling and no cross-query sort.
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
//!   queries of every batch, `BlastnReportInputs::queries_before`), and outfmt 6 has none. LOSAT's
//!   BLASTN writes outfmt 0, 6 and 7 only (no XML, JSON, SAM or HTML).
//! - the errors: an error of NCBI's that a batch meets (a reader error, an `eEOF` batch, a batch
//!   of queries without letters, `ErrorAfterInvalidBatches`, a header file that cannot be read)
//!   comes after the reports of the batches before, without the epilog, in both: here they were
//!   written when the loop took them, so the error paths write only what follows them (the
//!   batch's title and reader warnings). An error that LOSAT raises where NCBI would not stop
//!   (LOSAT's refusals) is written without the reports of the batches before; so the switch acts
//!   only where no such refusal can come after the first batch (`off`): scoring options with
//!   Karlin-Altschul values and greedy gap costs that LOSAT supports (the refusals of
//!   `search_query_batch` for scoring options without values or large greedy gap costs), no
//!   negative `CHUNK_SIZE` (the refusal in `search_query_chunks`), search spaces that cannot
//!   overflow (`SearchSpace`), and subjects whose report ids are all supported (a subject set's
//!   `tabular_subject_ids`, made once here; with `-db` every definition line is read before the
//!   search, as blastp's block search does, `database::x_titles_allow_streaming`). LOSAT's
//!   internal invariant errors and panics are not argued; a write error of the output stops the
//!   run where it happens (at the end without the switch).
//! - what the report path writes to standard error (the query warnings) keeps its bytes and its
//!   place before each query's report; `-verbose` and the debug and trace output, which a batch
//!   writes while it is searched and the reports write again, keep the reports at the end
//!   (`Off::Debug`). `LOSAT_TIMING`'s blocks are printed when the loop takes a batch, as before,
//!   so with the switch a block comes after the warnings of the batches before it.
//!
//! Off (the reports after the loop, as without the switch): the switch not set; outputs other
//! than the command line's (the hit records and the format observer of the Web adapter, whose
//! HSP numbers count over the whole search, and file sinks, which are created when their report
//! is written); `-verbose` and the debug output (the same list as LOSAT_X_QBATCHPAR's); and the
//! refusals above. The command line (`-out` file or standard output) and ABI v1 are streamed.
//!
//! LOSAT_X_QBATCHPAR: the writes happen while the slots search the batches ahead (they write
//! nothing; the batch loop alone writes), which is the point for the many-batch searches. The
//! error paths still stop every slot before they write. A slot's panic, which ends a release
//! build at once, then ends it after the reports of the batches taken so far.
//!
//! The shadow (`LOSAT_X_RPTSTREAMSHADOW`, implies the switch) keeps the final lists and values of
//! every batch as the off path does, copies every byte that the batch writes send to each format
//! and to standard error, and at the end (and on the error paths that write the reports of the
//! batches before) writes the reports the off path's way into memory and compares the bytes: on
//! a difference it prints `[X_SHADOW] RPTSTREAM ... DIFFERENT` and panics; at the end of the
//! search it prints `[X_SHADOW] RPTSTREAM batches= bytes= checks= off= identical` (batches: the
//! batches whose reports were written by the loop; bytes: the bytes compared; checks: the
//! comparisons; off: why the switch did not act, `none` when it did). With `LOSAT_X_STATS` and the
//! switch, `[X_STATS] RPTSTREAM batches= queries= off=`.

use std::io::{self, Write};
use std::ops::Range;
use std::sync::{Arc, OnceLock};

use anyhow::Result;

use crate::api::local_blast::{FormatOutput, OutputSink, ReportOutputs};
use crate::blastinput::seq_range::{Placements, QueryInput};
use crate::config::NuclScoringSpec;
use crate::stats::NuclValues;
use crate::utils::dust::MaskedInterval;
use crate::utils::threading::SearchPool;

use super::{
    blastn_trace, check_greedy_gap_costs, configure_task, database, record_title, BlastnArgs,
    BlastnOutputFormat, QueryBatching, SubjectLetters, SubjectRanges, SubjectSource,
};

// No NCBI counterpart: reads the switches once (any value turns them on; the shadow implies the
// switch; neither is in the strict set). They choose when LOSAT writes the same reports.
/// 0 = off, 1 = LOSAT_X_RPTSTREAM, 2 = LOSAT_X_RPTSTREAMSHADOW.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_RPTSTREAMSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_RPTSTREAM").is_some() {
            1
        } else {
            0
        }
    })
}

/// Why the reports of a search are written after its batch loop (`Stream::start`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// Neither switch is set.
    NotSet,
    /// Hit records, a format observer or a file sink (`ReportOutputs` other than the command
    /// line's).
    Outputs,
    /// `-verbose`, the debug variables, the blastn trace or the thread diagnostics.
    Debug,
    /// Scoring options without Karlin-Altschul values, or greedy gap costs that LOSAT refuses.
    Scoring,
    /// A negative `CHUNK_SIZE`.
    ChunkSize,
    /// A database length times a query length that could overflow NCBI's 64-bit search space.
    SearchSpace,
    /// A subject of the set whose tabular id LOSAT refuses.
    SubjectIds,
    /// A database whose definition lines do not all read as the reports read them.
    Titles,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Outputs => "outputs",
            Off::Debug => "debug",
            Off::Scoring => "scoring",
            Off::ChunkSize => "chunk_size",
            Off::SearchSpace => "search_space",
            Off::SubjectIds => "subject_ids",
            Off::Titles => "titles",
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

/// The streaming of one search's reports: whether it acts (`off`), what it has written, and the
/// shadow's copies. Its drop prints the summary lines.
pub(super) struct Stream {
    off: Option<Off>,
    cache: Option<SubjectCache>,
    shadow: Option<Shadow>,
    batches: u64,
    queries: usize,
}

impl Stream {
    /// The streaming of a search with these outputs and subjects, `lengths` the searched queries'
    /// lengths (`run_in_pool`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn start(
        args: &BlastnArgs,
        outputs: &ReportOutputs<'_>,
        output_formats: &[BlastnOutputFormat],
        source: SubjectSource<'_>,
        subject_ranges: SubjectRanges<'_>,
        lengths: &[usize],
        batching: QueryBatching,
        pool: &SearchPool<'_>,
    ) -> Self {
        let (off, cache) = match check(
            args,
            outputs,
            output_formats,
            source,
            subject_ranges,
            lengths,
            batching,
            pool,
        ) {
            Ok(cache) => (None, cache),
            Err(off) => (Some(off), None),
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
                "[X_STATS] RPTSTREAM batches={} queries={} off={off}",
                self.batches, self.queries
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

/// Whether a search writes debug or trace output while its batches are searched (the list of
/// `x_qbatch_par_off`): the reports then stay after the loop, where that output comes before
/// them.
fn debug_output(args: &BlastnArgs) -> bool {
    args.verbose
        || [
            "BLEMIR_DEBUG",
            "LOSAT_DEBUG_BLASTN",
            "LOSAT_DEBUG_COORDS",
            "LOSAT_DEBUG_COORDS_START",
        ]
        .iter()
        .any(|name| std::env::var_os(name).is_some())
        || blastn_trace::enabled()
        || crate::utils::threading::diagnostics_enabled()
}

/// Why the switch does not act for this search, or the subject values made once for a subject
/// set (`Stream::start`; the module doc has the reasons).
#[allow(clippy::too_many_arguments)]
fn check(
    args: &BlastnArgs,
    outputs: &ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    source: SubjectSource<'_>,
    subject_ranges: SubjectRanges<'_>,
    lengths: &[usize],
    batching: QueryBatching,
    pool: &SearchPool<'_>,
) -> std::result::Result<Option<SubjectCache>, Off> {
    if mode() == 0 {
        return Err(Off::NotSet);
    }
    if outputs.hits.is_some()
        || outputs.observer.is_some()
        || outputs
            .formats
            .iter()
            .any(|format| matches!(format.sink, OutputSink::File(_)))
    {
        return Err(Off::Outputs);
    }
    if debug_output(args) {
        return Err(Off::Debug);
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
            Ok(Some(SubjectCache { ids, titles, masks }))
        }
        SubjectSource::Database(db) => {
            if database::x_titles_allow_streaming(db, pool) {
                Ok(None)
            } else {
                Err(Off::Titles)
            }
        }
    }
}
