//! BLASTX's BLAST XML (outfmt 5): the `CBlastOutput` document of one query batch's final hit
//! lists (`report::blast_xml`), one query at a time.
//!
//! NCBI makes the Iteration of each query from its Seq-aligns, in the order the Seq-annot
//! formats show them (`seqalign`): subject by subject in hit-list order, each subject's HSPs
//! sorted by e-value (`render_queries` sorts them). What `s_SeqAlignSetToXMLHsps` reads of an
//! HSP is its Seq-align (coordinates in the records, score, e-value and bit score) and the two
//! rows of the alignment, from which it counts identities and positives again; nothing of the
//! engine's counts and nothing of the outfmt 0 midline. The rows are the ones outfmt 0 shows
//! (`report::alignment`), the masks the context's nucleotide intervals
//! (`QueryContext::dna_masks`, after the core's protein-to-DNA conversion).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-744
//! (`PrintOneResultSet` is called once per query's result set; the accumulated lists hold that
//! one query.)
//!
//! The XML2 and JSON documents of outfmt 15 and 16, and 14 and 13 as side files with a master file
//! (`report::blast_xml2`), take the same path
//! with another writer: `write_blastx_xml` is generic over `StructuredXml`, and the run keeps
//! one `Xml2Writer` across the query batches (`xml2_writer`), as it passes `queries_before` to
//! outfmt 5's. Its document start is written with the prolog, before the first batch.
//!
//! Every refusal (`XmlRefusal`: a title NCBI's definition line generator reads past the end of,
//! `names::Refusal`, where NCBI crashes) comes before the first byte of the document, though the
//! query file is read batch by batch:
//!
//! - `prescan_queries` reads the query file once more before the search (only with outfmt 5,
//!   13 to 16) and refuses the first query whose title the writer would refuse;
//! - when some `-subject` title would be refused if its subject had a hit
//!   (`subjects_may_be_refused`), or some query's message would be (a title that is refused only
//!   in the message of a query NCBI did not search), or (with `-db`) the taxonomy database is
//!   on the search path and the lookup of the taxid of some sequence of the database would be
//!   refused (a damaged record, which outfmt 5 reads
//!   in each batch's prefetch, `ReportedSubjects::batch_prefetch`: checked once before the
//!   first byte, `taxdb_may_refuse`, with the names for XML2 and JSON only), or the run's
//!   fetches may be refused for NCBI's taxonomy cache (`BlastxDatabase::cache_may_refuse`: outfmt
//!   5's prefetch, and the XML2 and JSON reports' loads of a batch that reports a subject with a
//!   positive taxid, `read_reported_taxonomy`; after the prolog), so that a later batch could
//!   refuse it,
//!   the run holds the document in memory (`HeldOutput`; with outfmt 13 and 14 the side files'
//!   bytes too, `HeldSideFiles`) and writes it at the epilog; a refusal found while the batches
//!   are reported drops it, with nothing written (an empty master file, no side file).

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{self, Write};
use std::rc::Rc;
use std::sync::{Arc, Mutex, MutexGuard};

use anyhow::{ensure, Context, Result};

use super::args::ResolvedOptions;
use super::empty_rows::{
    row_range_exception, rows_empty, rows_range, std_seg_lengths_not_accurate, SubjectLoc,
};
use super::input::FastaRecord;
use super::report::{
    display_query_residue, has_seqalign, query_endpoints, reported, ShownSubjects,
};
use super::results::Hsp;
use super::runtime::BatchResults;
use crate::algorithm::common::database::TaxdbNotes;
use crate::algorithm::common::database::TaxonomyRefusal;
use crate::common::GapEditOp;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names::{self, QueryText};
use crate::report::blast_xml::subjects::{HitNames, SubjectNames};
use crate::report::blast_xml::CaughtException;
use crate::report::blast_xml::{
    filter, hit_ranges, written_queries, BlastXmlWriter, QueryMask, StructuredXml, XmlHit, XmlHsp,
    XmlQuery, XmlRun, XmlStats, REFERENCE_GAPPED_BLAST,
};
use crate::report::blast_xml2::framing::side_file_open_error;
use crate::report::blast_xml2::{
    RunSideFiles, SideFiles, Xml2Document, Xml2Kind, Xml2Run, Xml2Target, Xml2Writer,
};
use crate::report::query_warnings::{EMPTY_DATABASE_MESSAGE, INVALID_QUERY_MESSAGE};
use crate::report::serial::seqalign::{read_num_ident, seqalign_evalue};

const PROGRAM: &str = "BLASTX";

/// LOSAT's refusal of a title that NCBI crashes on (`names::Refusal`; not an error of NCBI's
/// engine: the run ends with exit status 1 and the message, `native::report_error` passes it
/// on).
#[derive(Debug)]
pub(super) struct XmlRefusal(String);

impl std::fmt::Display for XmlRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for XmlRefusal {}

impl XmlRefusal {
    /// The error of the taxonomy lookups of outfmt 15 and 16 (`ReportedSubjects::for_batch`): a
    /// [`TaxonomyRefusal`] becomes a refusal that carries its message (it drops what the run
    /// holds); any other error, such as a header or taxonomy file that cannot be read, is
    /// returned as it is, so that the run writes what it holds, as NCBI writes what it has
    /// before a reader error.
    pub(super) fn of_taxonomy(error: anyhow::Error) -> anyhow::Error {
        if error.is::<TaxonomyRefusal>() {
            XmlRefusal(format!("{error:#}")).into()
        } else {
            error
        }
    }

    /// The error of a query's writes: LOSAT's refusal of a load that the stream leaves out
    /// (`report::query_warnings::CutRefusal`) becomes a refusal that carries its message (it drops
    /// what the run holds, as the run's other refusals of the taxonomy cache's slot do); any other
    /// error is returned as it is.
    pub(super) fn of_cut(error: anyhow::Error) -> anyhow::Error {
        let refusal = error
            .downcast_ref::<io::Error>()
            .and_then(crate::report::query_warnings::CutRefusal::of)
            .map(str::to_string);
        match refusal {
            Some(text) => XmlRefusal(text).into(),
            None => error,
        }
    }
}

/// Whether `error` is one of LOSAT's refusals for which a run that holds its BLAST XML, XML2 or
/// JSON document drops what it holds, so that the refusal follows no byte: a title NCBI crashes
/// on (`XmlRefusal`), the database counts of the reports that NCBI's sums overflow
/// (`ReportSumRefusal`), a reported database sequence with parsed Seq-ids
/// (`ParsedRefusal`), or a taxonomy database name or record (`TaxonomyRefusal`). After any other
/// error the run writes what it holds, the bytes of the batches before it, as the run without a
/// hold does (the command line's `HeldOutput`; a library call's `HeldFormats`).
pub(super) fn drops_held_document(error: &anyhow::Error) -> bool {
    error.is::<XmlRefusal>()
        || error.is::<crate::algorithm::common::database::ReportSumRefusal>()
        || error.is::<crate::algorithm::common::database::ParsedRefusal>()
        || error.is::<TaxonomyRefusal>()
}

fn refuse(refusal: names::Refusal, what: &str) -> anyhow::Error {
    XmlRefusal(refusal.error(PROGRAM, what).to_string()).into()
}

/// The masked regions of a query as NCBI hands them to the XML writer
/// (`Blast_GetSeqLocInfoVector`, blast_aux.cpp:930-959): per context, in context order, the
/// context's nucleotide intervals moved into the record (`seq_range::record_mask`: other than an
/// empty one and the whole interval), with the context's frame, which counts from the ends of the
/// searched letters (`-query_loc`); none for a batch that NCBI did not search.
pub(super) fn query_masks(
    batch: &BatchResults,
    qidx: usize,
    query: &FastaRecord,
) -> Vec<QueryMask> {
    let (length, offset) = (query.sequence.len(), query.placement().offset);
    let mut masks = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-224
    // The result of a batch that NCBI did not search has no masks (RPSTBLASTN's merge of several
    // volumes takes volume 0's, rpsblast_local.cpp:174-175).
    if batch.search_skipped {
        return masks;
    }
    for context in &batch.prepared.contexts[qidx * 6..qidx * 6 + 6] {
        for &(from, to) in &context.dna_masks {
            let Some((from, to)) =
                crate::blastinput::seq_range::record_mask(from, to, length, offset)
            else {
                continue;
            };
            masks.push(QueryMask {
                from,
                to,
                frame: i32::from(context.frame),
            });
        }
    }
    masks
}

/// The `Statistics` of a query: the first valid context's length adjustment and search
/// space, with the gapped Karlin block (the table value; the ideal ungapped block of the
/// context for `-ungapped`); the unsearched values for a query without a valid context.
fn query_stats(
    batch: &BatchResults,
    qidx: usize,
    gapped: bool,
    gapped_karlin: &crate::stats::tables::KarlinParams,
) -> XmlStats {
    let Some(parameters) = batch.parameters[qidx * 6..qidx * 6 + 6]
        .iter()
        .find(|parameters| parameters.valid)
    else {
        return XmlStats::NONE;
    };
    let karlin = if gapped {
        gapped_karlin
    } else {
        &parameters.ungapped
    };
    XmlStats {
        hsp_len: parameters.length_adjustment,
        eff_space: parameters.search_space,
        // RPSTBLASTN's footer K of a batch's first query (`BatchResults::gapped_k`, 1.2 times the
        // profile's) is the XML's `Statistics_kappa` too: one `kbp_gap` block per query.
        kappa: match &batch.gapped_k {
            Some(gapped_k) if gapped => gapped_k[qidx],
            _ => karlin.k,
        },
        lambda: karlin.lambda,
        entropy: karlin.h,
    }
}

/// The gapped Karlin block of the search: the table value of the matrix and the gap costs
/// (`parameters::gapped_karlin`; an ungapped search reads its contexts' blocks instead).
fn gapped_karlin(options: &ResolvedOptions) -> crate::stats::tables::KarlinParams {
    super::parameters::gapped_karlin(options)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
/// Hands the observer the statistics of every query of `batches`, the values of the XML reports,
/// each query numbered among the run's queries (`queries_before` queries before `batches`).
pub(super) fn observe_query_stats(
    options: &ResolvedOptions,
    batches: &[BatchResults],
    queries_before: usize,
    observer: &mut dyn crate::api::local_blast::FormatObserver,
) {
    let karlin = gapped_karlin(options);
    for batch in batches {
        for qidx in 0..batch.queries.len() {
            observer.query_stats(
                queries_before + batch.query_ordinal + qidx,
                query_stats(batch, qidx, options.gapped, &karlin),
            );
        }
    }
}

/// The messages of a query's result (the warnings the other formats write to standard error,
/// `report::query_warning`): the Karlin-Altschul message of a batch NCBI did not search, and
/// the message of a record without letters.
fn query_messages(query: &FastaRecord, search_skipped: bool, filtered_empty: bool) -> Vec<String> {
    let mut messages = Vec::new();
    if search_skipped {
        messages.push(INVALID_QUERY_MESSAGE.to_string());
    }
    if query.sequence.is_empty() {
        messages.push("Sequence contains no data".to_string());
    }
    if filtered_empty && !search_skipped {
        messages.push(EMPTY_DATABASE_MESSAGE.to_string());
    }
    messages
}

/// The names of a `Hit`, for `-subject` (`Subject_N` and the generated definition line) or a
/// database's sequence (`gnl|BL_ORD_ID|n`, the stored title as stored).
fn subject_names(
    shown: &ShownSubjects<'_>,
    subject: &FastaRecord,
    s_idx: usize,
) -> Result<SubjectNames> {
    match shown {
        ShownSubjects::Fasta { .. } => match subject.believed() {
            Some(ids) => SubjectNames::of_believed(ids, &subject.title, true),
            None => SubjectNames::of_record(subject.internal_id.as_bytes(), &subject.title, true),
        }
        .map_err(|refusal| refuse(refusal, "a subject title")),
        ShownSubjects::Database { reported, .. } => Ok(match reported.parsed_subject(s_idx) {
            Some(parsed) => SubjectNames::of_parsed(parsed, reported.parsed_taxonomy(s_idx)),
            None => SubjectNames::of_database_sequence(
                reported.oid(s_idx),
                &reported.stored_titles()[s_idx],
                // XML2 and JSON alone show the taxid (`ReportedSubjects` reads it for them only).
                reported.taxid(s_idx),
                reported.taxonomy(s_idx),
            ),
        }),
    }
}

/// The rows of an HSP's alignment into `qrow` and `srow` (cleared first): the letters
/// `report::alignment` shows in outfmt 0 (the query's codons translated, the subject's
/// residues, `-` for a gap), without the rest of what it builds (the masked row, the midline
/// counts), which the XML counts again from the rows itself. Returns the HSP's endpoints
/// `(query_from, query_to)` as the Seq-align has them, lower first.
fn xml_rows(
    hsp: &Hsp,
    subject: &FastaRecord,
    query: &FastaRecord,
    code: &crate::utils::genetic_code::GeneticCode,
    qrow: &mut Vec<u8>,
    srow: &mut Vec<u8>,
) -> Result<(i32, i32)> {
    let raw = &hsp.hsp;
    let (q_end, s_end) = hsp.seqalign_ends();
    qrow.clear();
    srow.clear();
    let (mut qi, mut si) = (usize::try_from(raw.q_start)?, usize::try_from(raw.s_start)?);
    let whole = [GapEditOp::Sub(u32::try_from(raw.q_end - raw.q_start)?)];
    let script: &[GapEditOp] = if hsp.edit_script.is_empty() {
        &whole
    } else {
        &hsp.edit_script
    };
    for &op in script {
        let (count, consume_q, consume_s) = match op {
            GapEditOp::Sub(n) => (n, true, true),
            GapEditOp::Ins(n) => (n, true, false),
            GapEditOp::Del(n) => (n, false, true),
        };
        for _ in 0..count {
            qrow.push(if consume_q {
                display_query_residue(&query.sequence, raw.frame, qi, code)
            } else {
                b'-'
            });
            srow.push(if consume_s {
                subject.sequence[si]
            } else {
                b'-'
            });
            qi += usize::from(consume_q);
            si += usize::from(consume_s);
        }
    }
    ensure!(
        qi == q_end as usize && si == s_end as usize,
        "BLASTX edit script does not reach final endpoints"
    );
    let (qstart, qend) = query_endpoints(query.sequence.len(), raw.frame, raw.q_start, q_end);
    Ok((qstart.min(qend), qstart.max(qend)))
}

/// The writer of the XML2 and JSON documents of a BLASTX run, with the side files of outfmt 13 and
/// 14 (`RunSideFiles`) of the run.
pub(super) type RunXml2 = Xml2Writer<RunSideFiles>;

/// The output of a run with outfmt 5, 13, 14, 15 or 16: it passes the bytes on, or, once `hold`
/// is called, keeps them in memory until `release` (the epilog), so that a refusal found while
/// the batches are reported leaves nothing written. A run that ends with any other error writes
/// what it held (the bytes of the batches before it, as the run without a hold does); `discard`
/// is for the refusal alone.
///
/// Outfmt 13 and 14 write their reports to side files: `side` is the hold of those files, which
/// follows this output's `hold`, `discard` and `release` (the master file is this output).
#[derive(Clone)]
pub(super) struct HeldOutput(Rc<RefCell<HeldState>>);

struct HeldState {
    inner: Box<dyn Write>,
    holding: bool,
    held: Vec<u8>,
    discarded: bool,
    side: Option<HeldSideFiles>,
    /// The run's standard error (`HeldOutput::stderr`), and what the run writes there while the
    /// output is held: each text with the length of the held output before it.
    stderr: Option<(Box<dyn Write>, Vec<(usize, Vec<u8>)>)>,
    /// The command line's report stream behind `inner` (`HeldOutput::with_report`), which takes
    /// NCBI's regions and flushes (`ostream::ReportWrite`).
    report: Option<crate::report::ostream::ReportHandle>,
    /// While the output is held, in the order they were made: the end of each insertion (each
    /// write, one of NCBI's), the regions and flushes that the formatter declared to the report's
    /// stream (`ostream::Cut`), at the held output's length then, and the index in the texts
    /// before which each comes.
    events: Vec<(usize, usize, crate::report::ostream::Cut)>,
}

impl HeldState {
    /// Writes the held output with the texts written to standard error while it was held, each
    /// after the output before it (NCBI's `cerr` flushes `cout` first), and empties both: each of
    /// the held insertions by itself, with the regions and flushes declared between them, in the
    /// order they were made with the texts (`ostream::CutCursor::write_until`), so that the
    /// report's stream fails where NCBI's does. Where NCBI's formatter catches a failed insertion
    /// (`ostream::is_caught_report_write_failure`: an `Iteration` of BLAST XML, whose post the
    /// stream made), the rest of that region is not written, nor the texts of its reads (NCBI
    /// makes no later `Hit`), as the cursor skips it; after an exception that ends the run,
    /// nothing more.
    fn write_held(&mut self) -> std::io::Result<()> {
        use crate::report::ostream::{is_caught_report_write_failure, Cuts};
        let held = std::mem::take(&mut self.held);
        let texts = self
            .stderr
            .as_mut()
            .map(|(_, texts)| std::mem::take(texts))
            .unwrap_or_default();
        let events = std::mem::take(&mut self.events);
        let mut cuts = Cuts::default();
        cuts.reserve(events.len());
        // By event: the number of texts written before it.
        let mut befores = Vec::with_capacity(events.len());
        for (at, before, cut) in events {
            cuts.push(at, cut);
            befores.push(before);
        }
        let HeldState {
            inner,
            stderr,
            report,
            ..
        } = self;
        let mut target = HeldTarget {
            inner: &mut **inner,
            report: report.as_ref(),
        };
        let mut cursor = cuts.cursor();
        // A failure that NCBI's formatter catches ended its region (the cursor skips the rest).
        let caught = |written: std::io::Result<()>| match written {
            Err(error) if is_caught_report_write_failure(&error) => Ok(()),
            written => written,
        };
        let mut written = 0;
        for (k, (at, text)) in texts.into_iter().enumerate() {
            let cuts_end = befores.partition_point(|&before| before <= k);
            caught(cursor.write_until(&mut target, &held, written, at, cuts_end))?;
            written = written.max(at);
            if cursor.skipping() {
                continue;
            }
            // NCBI's writes before the text are made, then the text, whose stream (tied) makes
            // them first too; where they fail in a caught region, the region ends there and its
            // texts are not written.
            let flushed = target.inner.flush().and_then(|()| match stderr.as_mut() {
                Some((stderr, _)) => stderr.write_all(&text).and_then(|()| stderr.flush()),
                None => Ok(()),
            });
            match flushed {
                Err(error) if is_caught_report_write_failure(&error) => {
                    cursor.caught_outside(&mut target)?;
                }
                flushed => flushed?,
            }
        }
        caught(cursor.write(&mut target, &held, written, held.len()))?;
        caught(target.inner.flush())
    }
}

/// The held output's writer (`inner`) and its report stream (`report`), which takes NCBI's
/// regions and flushes: what `HeldState::write_held` writes the held output to.
struct HeldTarget<'h> {
    inner: &'h mut dyn Write,
    report: Option<&'h crate::report::ostream::ReportHandle>,
}

impl Write for HeldTarget<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.inner.write(bytes)
    }

    fn write_all(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        self.inner.write_all(bytes)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

impl crate::report::ostream::ReportWrite for HeldTarget<'_> {
    fn site(&mut self, site: crate::report::ostream::Site) -> std::io::Result<()> {
        self.report.map_or(Ok(()), |report| report.declare(site))
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        self.report.map_or(Ok(()), |report| report.ncbi_flush())
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        self.report.map_or(Ok(()), |report| report.resolve(pieces))
    }
}

impl HeldOutput {
    pub(super) fn new(inner: Box<dyn Write>, side: Option<HeldSideFiles>) -> Self {
        HeldOutput(Rc::new(RefCell::new(HeldState {
            inner,
            holding: false,
            held: Vec::new(),
            discarded: false,
            side,
            stderr: None,
            report: None,
            events: Vec::new(),
        })))
    }

    /// This output before the command line's report stream `report` (`native::run_inner`): the
    /// regions and flushes that the formatter declares (`ostream::ReportWrite`) go to it, kept in
    /// order with the bytes while the output is held.
    pub(super) fn with_report(self, report: crate::report::ostream::ReportHandle) -> Self {
        self.0.borrow_mut().report = Some(report);
        self
    }

    /// A region or flush declared to the report's stream: at once, or held with the output.
    fn declare(&self, event: crate::report::ostream::Cut) -> std::io::Result<()> {
        let mut state = self.0.borrow_mut();
        if state.discarded {
            return Ok(());
        }
        if state.holding {
            let at = state.held.len();
            let texts = state.stderr.as_ref().map_or(0, |(_, texts)| texts.len());
            state.events.push((at, texts, event));
            return Ok(());
        }
        match (event, &state.report) {
            (crate::report::ostream::Cut::Site(site), Some(report)) => report.declare(site),
            (crate::report::ostream::Cut::Flush, Some(report)) => report.ncbi_flush(),
            _ => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:493-500
    /// The run's standard error `inner` through this output (`HeldStderr`): while the output is
    /// held, what the run writes there is kept with the place in the output where NCBI, which
    /// writes its document as it goes, writes it (the taxonomy texts of a `Hit`'s read between
    /// the `Hit`s), and written there at `release`; at `discard` the texts alone.
    pub(super) fn stderr(&self, inner: Box<dyn Write>) -> HeldStderr {
        self.0.borrow_mut().stderr = Some((inner, Vec::new()));
        HeldStderr(Rc::clone(&self.0))
    }

    /// Keeps everything written from now on in memory.
    pub(super) fn hold(&self) {
        let mut state = self.0.borrow_mut();
        state.holding = true;
        if let Some(side) = &state.side {
            side.hold();
        }
    }

    /// Drops what is held, and everything written after (the run ends with a refusal).
    pub(super) fn discard(&self) {
        let mut state = self.0.borrow_mut();
        state.discarded = true;
        state.held.clear();
        if let Some(side) = &state.side {
            side.discard();
        }
        // What the run wrote to standard error in the meantime stays (before the refusal).
        if let Some((stderr, texts)) = state.stderr.as_mut() {
            for (_, text) in std::mem::take(texts) {
                let _ = stderr.write_all(&text);
            }
            let _ = stderr.flush();
        }
    }

    /// Writes what is held and passes the bytes on again: the side files first (in order), then
    /// the master file. A side file that cannot be opened stops the run there, as NCBI's does
    /// (`Cannot open output file`, exit status 1): the side files before it stay and the master
    /// file stays empty.
    pub(super) fn release(&self) -> Result<()> {
        let mut state = self.0.borrow_mut();
        state.holding = false;
        if let Some(side) = state.side.clone() {
            if let Err(error) = side.release() {
                state.discarded = true;
                state.held.clear();
                // What the run wrote to standard error in the meantime stays.
                if let Some((stderr, texts)) = state.stderr.as_mut() {
                    for (_, text) in std::mem::take(texts) {
                        let _ = stderr.write_all(&text);
                    }
                    let _ = stderr.flush();
                }
                return Err(error);
            }
        }
        state.write_held()?;
        Ok(())
    }
}

/// The run's standard error while its output may be held (`HeldOutput::stderr`).
pub(super) struct HeldStderr(Rc<RefCell<HeldState>>);

impl Write for HeldStderr {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut state = self.0.borrow_mut();
        let at = state.held.len();
        let holding = state.holding && !state.discarded;
        let Some((stderr, texts)) = state.stderr.as_mut() else {
            return Ok(bytes.len());
        };
        if holding {
            texts.push((at, bytes.to_vec()));
            Ok(bytes.len())
        } else {
            stderr.write(bytes)
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let mut state = self.0.borrow_mut();
        let holding = state.holding && !state.discarded;
        match state.stderr.as_mut() {
            Some((stderr, _)) if !holding => stderr.flush(),
            _ => Ok(()),
        }
    }
}

impl Write for HeldOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut state = self.0.borrow_mut();
        if state.discarded {
            Ok(bytes.len())
        } else if state.holding {
            state.held.extend_from_slice(bytes);
            let at = state.held.len();
            let texts = state.stderr.as_ref().map_or(0, |(_, texts)| texts.len());
            state
                .events
                .push((at, texts, crate::report::ostream::Cut::Insertion));
            Ok(bytes.len())
        } else {
            state.inner.write_all(bytes)?;
            Ok(bytes.len())
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let mut state = self.0.borrow_mut();
        if state.holding || state.discarded {
            Ok(())
        } else {
            state.inner.flush()
        }
    }
}

/// The formatter's regions and flushes, in order with the held bytes (`HeldOutput::declare`).
impl crate::report::ostream::ReportWrite for HeldOutput {
    fn site(&mut self, site: crate::report::ostream::Site) -> std::io::Result<()> {
        self.declare(crate::report::ostream::Cut::Site(site))
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        self.declare(crate::report::ostream::Cut::Flush)
    }
}

impl Drop for HeldState {
    fn drop(&mut self) {
        let texts = self
            .stderr
            .as_ref()
            .is_some_and(|(_, texts)| !texts.is_empty());
        if !self.discarded && (!self.held.is_empty() || texts) {
            let _ = self.write_held();
        }
    }
}

/// The side files of an outfmt 13 or 14 run (the files next to `-out`, `FsSideFiles`, or a
/// caller's receiver of them), which a hold keeps in memory: the n-th report is written at once,
/// or, while the run holds its output, at `release` (the epilog) in the order of the reports. A
/// refusal found while the batches are reported drops them, so that the run leaves no side file
/// (design OF5 4.4); any other error ends the run with the side files of the earlier reports
/// written (`SideState::drop`), as NCBI leaves them. The hold of the command line's document
/// (`HeldOutput`) and of a library call's formats (`HeldFormats`) lead it.
#[derive(Clone)]
pub(super) struct HeldSideFiles(Arc<Mutex<SideState>>);

struct SideState {
    files: RunSideFiles,
    holding: bool,
    discarded: bool,
    /// The reports held, as `(n, bytes)` in the order they were made.
    held: Vec<(u32, Vec<u8>)>,
}

impl HeldSideFiles {
    pub(super) fn new(files: RunSideFiles) -> Self {
        HeldSideFiles(Arc::new(Mutex::new(SideState {
            files,
            holding: false,
            discarded: false,
            held: Vec::new(),
        })))
    }

    fn state(&self) -> MutexGuard<'_, SideState> {
        self.0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    pub(super) fn hold(&self) {
        self.state().holding = true;
    }

    pub(super) fn discard(&self) {
        let mut state = self.state();
        state.discarded = true;
        state.held.clear();
    }

    /// Writes the held side files, each whole (a file is created, truncating one of that name).
    pub(super) fn release(&self) -> Result<()> {
        let mut state = self.state();
        state.holding = false;
        state.write_held()
    }

    /// The n-th report, whole: dropped after `discard`; kept while the run holds its output, once
    /// its file can be opened (where NCBI opens it, `SideFiles::check_openable`: a failure is
    /// NCBI's "Cannot open output file"); else written at once.
    fn put(&self, n: u32, report: &[u8]) -> Result<()> {
        let mut state = self.state();
        if state.discarded {
            return Ok(());
        }
        if state.holding {
            state
                .files
                .check_openable(n)
                .map_err(|_| side_file_open_error())?;
            state.held.push((n, report.to_vec()));
            return Ok(());
        }
        state.files.write_file(n, report)
    }
}

impl SideState {
    /// A file that cannot be created is NCBI's "Cannot open output file"; a write to one that was
    /// opened is the side-file write error of every program (`SideFiles::write_file`).
    fn write_held(&mut self) -> Result<()> {
        let held = std::mem::take(&mut self.held);
        for (n, bytes) in held {
            self.files.write_file(n, &bytes)?;
        }
        Ok(())
    }
}

impl Drop for SideState {
    fn drop(&mut self) {
        if !self.discarded {
            let _ = self.write_held();
        }
    }
}

/// A report from `HeldSideFiles::create`: its bytes go to the side files (`HeldSideFiles::put`) at
/// its first flush, or when it is dropped without one.
struct HeldReport {
    side: HeldSideFiles,
    n: u32,
    bytes: Vec<u8>,
    put: bool,
}

impl HeldReport {
    fn put(&mut self) -> io::Result<()> {
        if std::mem::replace(&mut self.put, true) {
            return Ok(());
        }
        let bytes = std::mem::take(&mut self.bytes);
        self.side.put(self.n, &bytes).map_err(io::Error::other)
    }
}

impl Write for HeldReport {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.put()
    }
}

impl Drop for HeldReport {
    fn drop(&mut self) {
        let _ = self.put();
    }
}

impl SideFiles for HeldSideFiles {
    /// A writer that hands the report over when it is flushed or dropped (the trait's other
    /// callers; a run writes through `write_file`). A run that holds its side files checks here
    /// that the file can be opened.
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>> {
        let state = self.state();
        if state.discarded {
            return Ok(Box::new(io::sink()));
        }
        if state.holding {
            state.files.check_openable(n)?;
        }
        Ok(Box::new(HeldReport {
            side: self.clone(),
            n,
            bytes: Vec::new(),
            put: false,
        }))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        self.state().files.listed(n)
    }

    fn write_file(&mut self, n: u32, report: &[u8]) -> Result<()> {
        self.put(n, report)
    }
}

/// No NCBI counterpart: the output streams of a library call in which a later batch may still be
/// refused by a BLAST XML, XML2 or JSON document (`local::report_batches`), held in the order in
/// which they are written: each write and each flush keeps its stream's index (`writer`,
/// `HeldWrite`), so that `release` makes them in that order, as the run would have made them
/// without the hold (the reports of every format, as a stream that they all go to reads them from
/// NCBI; a buffered stream's bytes reach their device before the text written after them on
/// another stream), and `discard` drops them. Outside `hold` the writes pass through. The side
/// files of an outfmt 13 or 14 document of the run follow it (`side_files`): held, written before
/// the streams and dropped as they are.
pub(super) struct HeldFormats<S: Write>(std::sync::Arc<std::sync::Mutex<FormatsState<S>>>);

struct FormatsState<S> {
    holding: bool,
    discarded: bool,
    sinks: Vec<S>,
    held: Vec<HeldWrite>,
    /// The side files of the run's outfmt 13 and 14 documents.
    sides: Vec<HeldSideFiles>,
}

/// What one stream of a `HeldFormats` was given while held, between the writes of the others:
/// `bytes`, then a flush where `flushed` (the writes of one stream with no flush between them
/// are one `HeldWrite`).
struct HeldWrite {
    sink: usize,
    bytes: Vec<u8>,
    flushed: bool,
}

/// A stream of a `HeldFormats`.
pub(crate) struct HeldFormatWriter<S: Write> {
    hold: HeldFormats<S>,
    sink: usize,
}

impl<S: Write> Clone for HeldFormats<S> {
    fn clone(&self) -> Self {
        HeldFormats(self.0.clone())
    }
}

impl<S: Write> HeldFormats<S> {
    pub(crate) fn new() -> Self {
        HeldFormats(std::sync::Arc::new(std::sync::Mutex::new(FormatsState {
            holding: false,
            discarded: false,
            sinks: Vec::new(),
            held: Vec::new(),
            sides: Vec::new(),
        })))
    }

    fn state(&self) -> std::sync::MutexGuard<'_, FormatsState<S>> {
        self.0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// A stream that writes to `sink`.
    pub(crate) fn writer(&self, sink: S) -> HeldFormatWriter<S> {
        let mut state = self.state();
        state.sinks.push(sink);
        HeldFormatWriter {
            hold: self.clone(),
            sink: state.sinks.len() - 1,
        }
    }

    /// The side files `files` of an outfmt 13 or 14 document of the run, which follow this hold
    /// (`blast_xml::HeldSideFiles`): the command line holds its document's side files the same way
    /// (`blast_xml::HeldOutput`).
    pub(crate) fn side_files(
        &self,
        files: crate::report::blast_xml2::RunSideFiles,
    ) -> crate::report::blast_xml2::RunSideFiles {
        let side = HeldSideFiles::new(files);
        let mut state = self.state();
        if state.discarded {
            side.discard();
        } else if state.holding {
            side.hold();
        }
        state.sides.push(side.clone());
        Box::new(side)
    }

    /// Keeps everything written from now on.
    pub(crate) fn hold(&self) {
        let mut state = self.state();
        state.holding = true;
        for side in &state.sides {
            side.hold();
        }
    }

    /// Writes what is held and passes the writes on again: the side files first (in order), then
    /// the streams' writes and flushes in their order, as the run without the hold made them. A
    /// side file that cannot be opened stops the run there, as NCBI's does (`Cannot open output
    /// file`, exit status 1): the side files before it stay, and the streams and the side files of
    /// the later documents get nothing of what was held. A write or flush of a stream that fails
    /// stops the release there with its error (as the run would have stopped at it), and nothing
    /// held after it is written then or later. The writes of a stream for which `refused` is true
    /// (a format that failed alone, `FormatRefusal`) are dropped: its sink gets none of them.
    pub(crate) fn release(&self, refused: &dyn Fn(usize) -> bool) -> Result<()> {
        let mut state = self.state();
        state.holding = false;
        let held = std::mem::take(&mut state.held);
        if state.discarded {
            return Ok(());
        }
        for (index, side) in state.sides.iter().enumerate() {
            if let Err(error) = side.release() {
                // Nothing more of what was held is written: neither the streams, nor the side
                // files of the documents after this one.
                for later in &state.sides[index + 1..] {
                    later.discard();
                }
                state.discarded = true;
                return Err(error);
            }
        }
        for HeldWrite {
            sink,
            bytes,
            flushed,
        } in held
        {
            if refused(sink) {
                continue;
            }
            let sink = &mut state.sinks[sink];
            let written =
                sink.write_all(&bytes)
                    .and_then(|()| if flushed { sink.flush() } else { Ok(()) });
            if let Err(error) = written {
                state.discarded = true;
                return Err(error.into());
            }
        }
        Ok(())
    }

    /// Drops what is held and everything written after.
    pub(crate) fn discard(&self) {
        let mut state = self.state();
        state.discarded = true;
        state.held.clear();
        for side in &state.sides {
            side.discard();
        }
    }

    /// The end of a run that stopped with `error` and held its output for its BLAST XML, XML2 or
    /// JSON documents: dropped after one of the refusals that such a hold is for
    /// (`blast_xml::drops_held_document`), written after any other error (the bytes of the batches
    /// before it, as the run without a hold writes them, a write that fails then being the run's
    /// error in place of `error`: the run without the hold stopped there, before it), as the
    /// command line's document (`blast_xml::HeldOutput`). `refused` as for `release`.
    pub(crate) fn end_document_with(
        &self,
        error: &anyhow::Error,
        refused: &dyn Fn(usize) -> bool,
    ) -> Result<()> {
        if drops_held_document(error) {
            self.discard();
            Ok(())
        } else {
            self.release(refused)
        }
    }
}

impl<S: Write> Write for HeldFormatWriter<S> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut state = self.hold.state();
        if state.discarded {
            return Ok(bytes.len());
        }
        if state.holding {
            match state.held.last_mut() {
                Some(held) if held.sink == self.sink && !held.flushed => {
                    held.bytes.extend_from_slice(bytes)
                }
                _ => state.held.push(HeldWrite {
                    sink: self.sink,
                    bytes: bytes.to_vec(),
                    flushed: false,
                }),
            }
            return Ok(bytes.len());
        }
        state.sinks[self.sink].write(bytes)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let mut state = self.hold.state();
        if state.discarded {
            return Ok(());
        }
        if state.holding {
            match state.held.last_mut() {
                Some(held) if held.sink == self.sink => held.flushed = true,
                _ => state.held.push(HeldWrite {
                    sink: self.sink,
                    bytes: Vec::new(),
                    flushed: true,
                }),
            }
            return Ok(());
        }
        state.sinks[self.sink].flush()
    }
}

/// The `QueryText` of a query, refused as the writer refuses it. `search_skipped`: NCBI did not
/// search the query's batch (its message then carries the query's id, made of its title).
fn query_text(
    program: &str,
    query: &FastaRecord,
    search_skipped: bool,
    filtered_empty: bool,
    has_hits: bool,
    merged_volumes: bool,
) -> Result<QueryText> {
    // RPSTBLASTN's results of several volumes keep no warning (`report::query_warning`).
    let messages = if merged_volumes {
        Vec::new()
    } else {
        query_messages(query, search_skipped, filtered_empty)
    };
    names::query_text(
        program,
        query.internal_id.as_bytes(),
        query.warning_id(),
        query.believed(),
        &query.title,
        false,
        &messages,
        has_hits,
    )
    .map_err(|error| anyhow::Error::new(XmlRefusal(error.to_string())))
}

/// What the read of the queries before the search found (`scan_queries`).
#[derive(Default)]
pub(super) struct QueryScan {
    /// The run stops at a batch after this many queries of the batches before it (a reader
    /// error, or a batch that NCBI fails at its set-up): the report of the last of them is never
    /// written where the document holds it back (`written_queries`). `None` when the run reads
    /// to the end of the file, or stops at its first batch.
    pub stops_after: Option<usize>,
    /// Under `-parse_deflines`, the scopes of the queries read, in the search's batches
    /// (`prescan_queries`; `QueryScopes::read_ahead`).
    pub query_scopes: Option<crate::blastinput::query_batch::QueryScopes>,
}

/// Reads the query file once before the search (`scan_queries`), with the search's reader
/// reading as the search reads (`settings`: same records and titles, cut to `range` as the search
/// cuts them: a record the range skips is never reported and is not checked) in the search's
/// batches of `batch_size` letters. A reader error ends the read (under `-parse_deflines`, a first
/// line that is not a defline; the run reads the file again and fails there as it always did,
/// after the reports of the batches before it), and so does a record it refuses for every format
/// (a title that is not UTF-8). A batch size of 0 reads no query (NCBI's first batch is empty and
/// the search stops at it): nothing is checked. Under `-parse_deflines` (`query_scopes`, without a
/// batch yet) the read adds the batches to the scopes of the queries, `checks_titles` as
/// `scan_queries` does.
///
/// `file` is the query file the search reads (`-query`, a regular file, named `query` in NCBI's
/// error when its position cannot be read or set): the read takes it and gives it back with the
/// scan, at the position it had, so the search's reading starts where it would have started. It
/// reads that open file itself, never a second handle on it, which the WASI targets cannot make.
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
// NCBI reads ahead of the search on the argument's one stream (`m_Ios`, c++/src/corelib/ncbiargs.cpp:727-735) and
// seeks it back to `orig_p`. This read has no NCBI counterpart (NCBI reads the query file once; this
// read is LOSAT's, for its refusals before the first byte), and it reads the same way: on the one
// open file, put back where it started.
#[allow(clippy::too_many_arguments)]
pub(super) fn prescan_queries(
    file: crate::blastinput::input_files::InputFile,
    query: &std::path::Path,
    settings: &super::input_source::QueryStreamSettings,
    options: &ResolvedOptions,
    lowercase: bool,
    batch_size: usize,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    holds_back: bool,
    checks_titles: bool,
    query_scopes: Option<crate::blastinput::query_batch::QueryScopes>,
) -> Result<(QueryScan, crate::blastinput::input_files::InputFile)> {
    use std::io::{Seek, SeekFrom};
    let inaccessible = || crate::cli::inaccessible("query", query);
    let mut file = file;
    let start = file.stream_position().map_err(|_| inaccessible())?;
    let scan = if batch_size == 0 {
        QueryScan::default()
    } else {
        let mut stream = settings.stream_before_search(file);
        let scan = scan_query_stream(
            &mut stream,
            options,
            lowercase,
            batch_size,
            range,
            holds_back,
            checks_titles,
            query_scopes,
        )?;
        file = stream.into_input();
        scan
    };
    file.seek(SeekFrom::Start(start))
        .map_err(|_| inaccessible())?;
    Ok((scan, file))
}

/// The read of `prescan_queries` on its stream of the query file.
#[allow(clippy::too_many_arguments)]
fn scan_query_stream<R: std::io::Read>(
    stream: &mut super::input_source::QueryStream<R>,
    options: &ResolvedOptions,
    lowercase: bool,
    batch_size: usize,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    holds_back: bool,
    checks_titles: bool,
    query_scopes: Option<crate::blastinput::query_batch::QueryScopes>,
) -> Result<QueryScan> {
    scan_queries(
        options,
        holds_back,
        checks_titles,
        query_scopes,
        &mut |deliver| stream.read_batches(lowercase, batch_size, range, deliver, &mut |_| Ok(())),
    )
}

/// Refuses the first query the writer would refuse (a title NCBI crashes on, `names::Refusal`)
/// in a report that the document writes, among the reader batches that `read` gives its
/// argument in order, so that no refusal can follow a written byte (`QueryScan`; the command
/// line's `prescan_queries`, the records of `local::run_local`).
///
/// The batches end where the run stops before a batch's search, whose records are never written
/// and need no check: the end of `read` with an error (a reader error, or another error of LOSAT's
/// reading, which the search meets in its place) and a batch that NCBI fails at its set-up
/// (`archive::written_batch`, which `options` decide). A document that holds each report back
/// until the next query's arrives (`holds_back`: a `-subject` XML2 or JSON document) never writes
/// the last query of the batches before such a stop: that query is checked only when the next
/// batch is written or `read` ends without an error. Without `checks_titles` (a run whose XML
/// format fails for its matrix name before any report, `xml_matrix_failure`) no title is checked.
///
/// Under `-parse_deflines` (`query_scopes`, without a batch yet) the read adds each batch that NCBI
/// does not fail to the scopes of the queries (`QueryScopes::add_batch_ids`; a batch that NCBI fails
/// reports nothing of its scope, review M3, 2026-11-20), which the search reads from
/// `QueryScan::query_scopes`. LOSAT's refusal of two queries of one scope with one id stops the run
/// at its batch, after the reports of the batches before, as a batch that NCBI fails does: the read
/// ends there, and the search gives the refusal where it adds that batch (`QueryScopes::read_ahead`;
/// line WQA D10 round 4, 2026-11-20). A refusal that the run's XML failure for its matrix name may
/// take (`ScopeRefusal::xml_takes`) does not end the read: the search's add of its batch decides
/// (`native::add_batch_scope`). The batches are the reader's (`read`), as NCBI counts them
/// (`query_batch::BatchCount`).
pub(super) fn scan_queries(
    options: &ResolvedOptions,
    holds_back: bool,
    checks_titles: bool,
    mut query_scopes: Option<crate::blastinput::query_batch::QueryScopes>,
    read: &mut dyn FnMut(&mut dyn FnMut(&[FastaRecord]) -> Result<()>) -> Result<()>,
) -> Result<QueryScan> {
    let mut scan = QueryScan::default();
    // The queries of the batches that the run writes, and the check of the last of them while
    // the document holds its report back.
    let mut written = 0usize;
    let mut held: Option<Result<()>> = None;
    // The batches added to the scopes, and the refusals of the batches: those that the run's XML
    // failure for its matrix name may take (`ScopeRefusal::xml_takes`), after which the read goes
    // on, and the one that ends it.
    let mut batches = 0usize;
    let mut refusals = Vec::new();
    let scanned = read(&mut |batch| {
        if !super::archive::written_batch(batch, options) {
            anyhow::bail!("NCBI stops the run at this batch");
        }
        if let Some(scopes) = query_scopes.as_mut() {
            let added = scopes.add_batch_records(
                super::input::scope_records(batch)
                    .into_iter()
                    .map(|(believed, length, place)| (believed, length, place.is_some())),
                written,
            );
            if let Err(error) = added {
                let xml_takes = crate::blastinput::query_batch::ScopeRefusal::of(&error)
                    .is_some_and(crate::blastinput::query_batch::ScopeRefusal::xml_takes);
                refusals.push((batches, error));
                if !xml_takes {
                    anyhow::bail!("LOSAT refuses this batch");
                }
            }
            batches += 1;
        }
        // A run whose XML format fails for its matrix name writes no report, and no title
        // (`xml_matrix_failure`): its batches are read for their scopes alone.
        if !checks_titles {
            written += batch.len();
            return Ok(());
        }
        // The next batch is written: so is the report held before it.
        if let Some(check) = held.take() {
            check?;
        }
        for (index, query) in batch.iter().enumerate() {
            let check = query_text(
                options.program().upper,
                query,
                false,
                false,
                true,
                options.merged_volumes(),
            )
            .map(drop);
            if holds_back && index + 1 == batch.len() {
                held = Some(check);
            } else {
                check?;
            }
        }
        written += batch.len();
        Ok(())
    });
    scan.query_scopes = query_scopes.map(|scopes| scopes.read_ahead(refusals));
    match scanned {
        Err(error) if error.is::<XmlRefusal>() => Err(error),
        // The end of the queries: the epilog writes the held report.
        Ok(()) => {
            if let Some(check) = held {
                check?;
            }
            Ok(scan)
        }
        Err(_) => {
            scan.stops_after = (written > 0).then_some(written);
            Ok(scan)
        }
    }
}

/// Whether the XML documents of a run (BLAST XML, outfmt 5; the XML2 and JSON documents, 13 to
/// 16; `outfmts` are the run's format numbers) may refuse after their first byte, the decision of
/// both the command line (`native::run_inner`, `HeldOutput::hold`) and the library
/// (`local::write_run`, `archive::ArchiveHold`), which then hold their output so that the refusal
/// follows no byte: with `-db` a later batch's taxonomy database or subject with parsed Seq-ids
/// (`BlastxDatabase::later_batch_may_refuse`) or the database counts the reports would print where
/// NCBI's sums overflow (`common::database::report_sum_refusal`, refused at the first report:
/// `StructuredXml::start` in `write_blastx_xml`), and without it a `-subject` title
/// (`subjects_may_be_refused`).
pub(super) fn documents_may_refuse(
    outfmts: &[i32],
    database: Option<&super::database::BlastxDatabase>,
    subjects: &[FastaRecord],
) -> bool {
    let documents = || {
        outfmts
            .iter()
            .filter(|&&outfmt| matches!(outfmt, 5 | 13..=16))
    };
    if documents().next().is_none() {
        return false;
    }
    match database {
        Some(database) => {
            let pairwise = database.pairwise_database();
            database.later_batch_may_refuse()
                || documents().any(|&outfmt| {
                    crate::algorithm::common::database::report_sum_refusal(
                        pairwise.num_sequences,
                        pairwise.total_letters_overflow,
                        outfmt,
                    )
                    .is_some()
                })
        }
        None => subjects_may_be_refused(subjects),
    }
}

/// Whether the title of some `-subject` would be refused if the subject had a hit (a title whose
/// `GenerateDefline` throws is the writers' catch, `SubjectNames::of_unknown`).
pub(super) fn subjects_may_be_refused(subjects: &[FastaRecord]) -> bool {
    subjects.iter().any(|subject| {
        names::xml_subject_names(subject.internal_id.as_bytes(), &subject.title, true).is_err()
    })
}

/// The document of outfmt 13, 14, 15 or 16 for `outfmt` (`None` for any other format).
pub(super) fn xml2_kind(outfmt: u8) -> Option<Xml2Kind> {
    match outfmt {
        13 => Some(Xml2Kind::Files13),
        14 => Some(Xml2Kind::Files14),
        15 => Some(Xml2Kind::Single15),
        16 => Some(Xml2Kind::Single16),
        _ => None,
    }
}

/// The writer of the XML2 or JSON document of a run (`kind`), which the run keeps across its
/// query batches: its count of reports and, with `-subject`, the report held back one query
/// (`Xml2Document`). `-subject` has one `Subject_N` per record of the file (`subjects`). `files`
/// are the side files of outfmt 13 and 14.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:544-554
pub(super) fn xml2_writer(
    kind: Xml2Kind,
    options: &ResolvedOptions,
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    files: RunSideFiles,
) -> RunXml2 {
    let params = Xml2Run {
        target: match shown {
            ShownSubjects::Fasta { .. } => Xml2Target::Subjects(subjects.len()),
            ShownSubjects::Database { .. } => Xml2Target::Database,
        },
        subject_names: match shown {
            ShownSubjects::Fasta { .. } => {
                crate::report::blast_xml2::subject_names(subjects.iter().map(FastaRecord::believed))
            }
            ShownSubjects::Database { .. } => None,
        },
        // `-ungapped` writes no gap costs; it needs `-comp_based_stats 0`, so no `cbs` either.
        gapped: options.gapped,
        cbs: i32::from(options.composition),
        query_gencode: Some(i32::from(options.query_gencode)),
        db_gencode: None,
    };
    let delayed = params.target.holds_back();
    Xml2Writer::with_document(Xml2Document::with_files(kind, delayed, files), params)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:121 (cited at
// `blast_xml2::reports_held`)
/// Where BLASTX's run of the formats `formats` (the command line's one, or the library's) stops
/// for a matrix name that NCBI's BLAST XML report data does not know
/// (`blast_xml::matrix::xml_formatting_failure`; the first of the formats to fail stops the run),
/// with the fact of the search: subject records (`subject_set`: `-subject`, `m_IsBl2Seq`); BLASTX
/// does not iterate. Decided once for the run and given to the batch whose report fails
/// (`write_blastx_xml`), to the end of a run of one query (`report::write_epilog`) and to the
/// reads of the titles before the search, which then check no title (`native::run_inner`,
/// `local::scan_xml_titles`).
pub(super) fn xml_matrix_failure(
    formats: &[ResolvedOptions],
    subject_set: bool,
) -> Option<crate::report::blast_xml::matrix::XmlMatrixFailure> {
    let outfmts: Vec<String> = formats
        .iter()
        .map(|options| options.outfmt.to_string())
        .collect();
    crate::report::blast_xml::matrix::xml_formatting_failure(
        outfmts.iter().map(String::as_str),
        &formats.first()?.matrix,
        subject_set,
        false,
    )
}

/// Writes the documents of every query of `batches` (one query batch, after the
/// `queries_before` queries of the batches before it) through `document` (BLAST XML, or the XML2
/// / JSON documents) to `writer`: the document's header before the first query of the run, each
/// query's Iteration or report, and nothing else; the end of the document is `write_epilog`'s.
/// What LOSAT refuses is refused before this call's first byte, for the reports that the document
/// writes: a query title or a title of a subject with hits that NCBI crashes on
/// (`names::Refusal`; with a hold, the whole run's bytes stay in memory, `HeldOutput`); every
/// other byte NCBI's stream writes or fails at by itself (`report::serial::xml`). A run that stops after
/// this batch (`stops`) never writes the `-subject` report of its last query, which NCBI holds
/// back (`written_queries`). NCBI's messages about the taxonomy database go to `diagnostics` at
/// the first `Hit` that needs them, those of the batch's prefetch where it begins.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_blastx_xml<W: StructuredXml>(
    writer: &mut impl crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
    document: &mut W,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
    queries_before: usize,
    stops: bool,
    xml_matrix: Option<crate::blastinput::query_batch::XmlFailingQuery>,
) -> Result<()> {
    let gapped_karlin = gapped_karlin(options);
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2389 (cited at
    // `blast_xml::matrix::XmlMatrixFailure`)
    // The formatting table is found by the name as typed (`strcmp`): a lower-case name stops the
    // report data of the run's first query, or where the report is held of NCBI's next query of
    // another Seq-id (`XmlMatrixFailure::SecondQueryOrEnd`; `xml_matrix`, the run's
    // `xml_matrix_failure` at its query of the read-ahead scopes, `QueryScopes::xml_failing_query`),
    // in the batch of that query (`BLAST engine error:`, exit 3), before any byte of the reports; a
    // batch before it writes nothing (its reader messages and request posts are the batch loop's).
    let Some(positives) = FormattingMatrix::new(&options.matrix) else {
        let failure =
            xml_matrix.context("LOSAT internal error: no run decision for an XML matrix")?;
        if failure.fails_in(queries_before..queries_before + queries.len()) {
            return Err(crate::report::blast_xml::matrix::unsupported_score_matrix());
        }
        return Ok(());
    };
    let (db_name, db_num, db_len, db_len_overflow) = match shown {
        ShownSubjects::Fasta { .. } => (String::new(), 0, 0, false),
        ShownSubjects::Database { name, pairwise, .. } => (
            name.to_string(),
            pairwise.num_sequences,
            pairwise.total_letters,
            pairwise.total_letters_overflow,
        ),
    };
    let filter = filter::filter_string(
        None,
        options.seg.filter_options(),
        None,
        None,
        options.soft_masking,
    );
    let run = XmlRun {
        program: options.program().lower,
        version: &format!(
            "{} {}",
            options.program().upper,
            crate::wording::report_version!()
        ),
        reference: REFERENCE_GAPPED_BLAST,
        db_name: &db_name,
        matrix: Some(&options.matrix),
        expect: options.evalue,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
        reward_penalty: None,
        filter: &filter,
        db_num,
        db_len,
        db_len_overflow,
        positives: Some(&positives),
        mask_char: b'X',
    };
    // The definition line and message of every query, and the names of the subjects with hits,
    // made (and refused) before the first byte, for the written reports.
    let reports = written_queries(
        batches.iter().map(|batch| batch.queries.len()).sum(),
        stops,
        document.holds_back(),
    );
    // The refusal of the database counts that the written reports would print, before any byte
    // of the batch (`StructuredXml::start`; the XML2 and JSON documents' start came with the
    // prolog, which a run that may refuse them holds, `native::run_inner`).
    let mut buf: Vec<u8> = Vec::new();
    document.start((reports > 0).then_some(&run), &mut buf)?;
    writer.write_all(&buf)?;
    buf.clear();
    let mut texts: Vec<QueryText> = Vec::new();
    let mut subject_cache: HashMap<usize, HitNames> = HashMap::new();
    // What the reports do at the `Hit` of each subject (`HitLoads`, the one answer of the batch:
    // a subject whose `Hit` throws, one that no report loads); a subject set's are all loaded.
    let every_loaded = crate::report::db_subject::HitLoads::default();
    let loads = match shown {
        ShownSubjects::Database { reported, .. } => reported.loads(),
        ShownSubjects::Fasta { .. } => &every_loaded,
    };
    // RPSTBLASTN: an entry of the profile database whose Seq-align id is the generated id
    // (`lcl|Query_N`) of a query of the batch that reports it is that query for NCBI's scope: its
    // Seq-align conversion raises an exception after half a document. Refused before the reports
    // (PS15a-1, followup `rps_lcl_query_title.md`), as rpsblast's XML is
    // (`report::rps_generated_scopes`, `common::database::refuse_generated_query_id`).
    let generated_scopes = super::report::rps_generated_scopes(options, queries, shown, batches);
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-711
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2358-2362
    // Every document writes the lists of each query that the cut to the hit list size keeps
    // (`TrimSeqAlign` counts following Seq-aligns of matching subject ids as one, as
    // `PruneSeqalign`): RPSTBLASTN's merged lists of several volumes, cut once per query
    // (`BatchResults::shown_lists`); every list otherwise.
    for batch in batches {
        for qidx in 0..batch.queries.len() {
            if texts.len() == reports {
                break;
            }
            let lists = batch.shown_lists(qidx, options);
            if let Some(scopes) = &generated_scopes {
                super::report::refuse_rps_generated_ids(
                    scopes,
                    shown,
                    lists,
                    batch.query_ordinal + qidx,
                    options,
                    |_, _| true,
                    crate::algorithm::common::database::XML_EFFECT,
                )?;
            }
            let query = &queries[batch.query_ordinal + qidx];
            let has_hits = lists.iter().any(|list| reported(list, options.gapped));
            texts.push(query_text(
                options.program().upper,
                query,
                batch.search_skipped,
                options.filtered_empty,
                has_hits,
                options.merged_volumes(),
            )?);
            for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                if let std::collections::hash_map::Entry::Vacant(slot) =
                    subject_cache.entry(list.oid)
                {
                    // A subject whose lookup throws where its `Hit` loads it has the names of
                    // NCBI's `catch` with the exception; one that no report loads gets no names
                    // (`HitNames::of`).
                    let parsed = match shown {
                        ShownSubjects::Database { reported, .. } => {
                            reported.parsed_subject(list.oid)
                        }
                        ShownSubjects::Fasta { .. } => None,
                    };
                    slot.insert(HitNames::of(loads, list.oid, parsed, || {
                        subject_names(shown, &subjects[list.oid], list.oid)
                    })?);
                }
            }
        }
    }
    let loc = SubjectLoc::of(shown);
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    let no_notes = TaxdbNotes::default();
    let mut taxdb_printer = match shown {
        ShownSubjects::Database { reported, .. } => reported.taxdb_notes(),
        ShownSubjects::Fasta { .. } => &no_notes,
    }
    .printer();
    // The rows of one HSP at a time, reused.
    let (mut qrow, mut srow): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut texts = texts.iter();
    for batch in batches {
        for qidx in 0..batch.queries.len() {
            // The lists that the cut keeps, as above.
            let lists = batch.shown_lists(qidx, options);
            let query = &queries[batch.query_ordinal + qidx];
            let Some(text) = texts.next() else {
                // The held report that is never written (`written_queries`); the one held before
                // it is written, or fails (`CaughtException`: its message follows).
                let post = match document.unwritten(&mut buf) {
                    Ok(()) => None,
                    Err(error) => match error.downcast::<CaughtException>() {
                        Ok(CaughtException(exception)) => {
                            Some(crate::ncbi_exception::error_post_line(
                                None,
                                Some(options.program().lower),
                                &document.caught_message(&exception),
                            ))
                        }
                        Err(error) => return Err(error),
                    },
                };
                write_with_posts(&buf, document.take_posts(), writer, diagnostics)?;
                buf.clear();
                if let Some(post) = post {
                    writer.flush()?;
                    diagnostics.write_all(&post)?;
                    diagnostics.flush()?;
                }
                break;
            };
            // NCBI's prefetch of the batch, where it begins (`TaxdbPrinter::batch_texts`).
            let prefetch = taxdb_printer.batch_texts(batch.query_ordinal + qidx);
            if !prefetch.is_empty() {
                writer.flush()?;
                diagnostics.write_all(&prefetch)?;
                diagnostics.flush()?;
            }
            let length = query.sequence.len();
            // A record cut to `-query_loc`: `query-len` is the interval's, the coordinates and
            // frames the record's (`RemapToQueryLoc`, `report::alignment`).
            let place = query.placement();
            let masks = query_masks(batch, qidx, query);
            let xml_query = XmlQuery {
                id: text.query_id(query.internal_id.as_bytes()),
                title: &text.title,
                len: i32::try_from(length)?,
                stats: query_stats(batch, qidx, options.gapped, &gapped_karlin),
                message: text.message.as_deref(),
                // Read by XML2 and JSON alone: a query without NCBI's alignment set
                // (`BatchResults::has_alignment_set`).
                searched: batch.has_alignment_set(options),
                masks: &masks,
            };
            // The texts of the subjects that the query's report reads one by one, written where
            // NCBI loads the subject of each `Hit` (`HitSink::position`, `ReadCursor`).
            let mut reads = crate::report::query_warnings::ReadCursor::new(std::slice::from_ref(
                shown.query_reads(batch.query_ordinal + qidx),
            ));
            // What NCBI posts where it loads the subject of each `Hit` (`HitPosts`): those texts
            // and the taxonomy database's messages.
            let mut posts = crate::report::query_warnings::HitPosts::new(reads.is_some());
            let written = document.query(&run, &xml_query, &mut buf, |sink| {
                // A copy of one record of the database that follows another is shown in its
                // `Hit`.
                let hits = hit_ranges(
                    lists,
                    |list| reported(list, options.gapped).then_some(list.oid),
                    |previous, oid| match shown {
                        ShownSubjects::Database { reported, .. } => {
                            match (
                                reported.parsed_subject(previous),
                                reported.parsed_subject(oid),
                            ) {
                                (Some(earlier), Some(later)) => {
                                    crate::report::db_subject::joins_as_copy(earlier, later)
                                }
                                _ => false,
                            }
                        }
                        ShownSubjects::Fasta { .. } => false,
                    },
                );
                for range in hits {
                    let first = lists[range.start].oid;
                    let names = subject_cache
                        .get(&first)
                        .context("names of a subject with hits")?
                        .names()?;
                    let hit = XmlHit::of(names, i32::try_from(subjects[first].placement().length)?);
                    // NCBI's messages about the taxonomy database: opening and the volume's
                    // lookup at the first `Hit` that needs them, the lookup of the `Hit` every
                    // time.
                    posts.hit(sink.position(), &hit, first as u32, |texts| {
                        taxdb_printer.at_names((batch.query_ordinal + qidx, first), names, texts)
                    })?;
                    sink.hit(&hit, |hsps| {
                        for list in lists[range]
                            .iter()
                            .filter(|list| reported(list, options.gapped))
                        {
                            let subject = &subjects[list.oid];
                            let subject_offset = subject.placement().offset as i32;
                            for hsp in list
                                .hsps
                                .iter()
                                .filter(|hsp| has_seqalign(hsp, options.gapped))
                            {
                                let identities = read_num_ident(hsp.identity as i64, None);
                                // The `catch` of an HSP whose subject's lookup throws comes before
                                // NCBI reads its rows (`HspSink::reads_alignment`).
                                if hsps.unaligned([XmlHsp::unaligned(
                                    hsp.hsp.score,
                                    hsp.bit_score,
                                    seqalign_evalue(hsp.evalue),
                                    identities,
                                )])? {
                                    continue;
                                }
                                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:173-175
                                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:243
                                // An HSP of length 0 (`empty_rows`): the Dense-seg of `-subject`'s
                                // mapped Std-seg throws; else the start of a row whose range the
                                // XML cannot take, the query's first (`empty_rows::rows_range`:
                                // empty rows with `-db` and the rows that `-subject`'s mapper
                                // leaves, `SubjectLoc::of_list`; RPSTBLASTN's row that starts with
                                // a gap). NCBI catches it in `s_BlastXMLAddIteration`
                                // (`CaughtException`).
                                let thrown = match loc.of_list(list, options.gapped) {
                                    SubjectLoc::Interval => {
                                        rows_empty(hsp).then(std_seg_lengths_not_accurate)
                                    }
                                    SubjectLoc::Whole => rows_range(hsp).map(row_range_exception),
                                };
                                if let Some(exception) = thrown {
                                    return Err(CaughtException(exception).into());
                                }
                                // The Seq-align's endpoints, lower first: the minus frame's
                                // outfmt 6 order is not the XML's.
                                let (query_from, query_to) =
                                    xml_rows(hsp, subject, query, &code, &mut qrow, &mut srow)?;
                                let (query_from, query_to) = (
                                    query_from + place.offset as i32,
                                    query_to + place.offset as i32,
                                );
                                let plus = hsp.hsp.frame > 0;
                                hsps.hsp(&XmlHsp {
                                    score: hsp.hsp.score,
                                    bit_score: hsp.bit_score,
                                    evalue: seqalign_evalue(hsp.evalue),
                                    query_from,
                                    query_to,
                                    hit_from: hsp.hsp.s_start + 1 + subject_offset,
                                    hit_to: hsp.seqalign_ends().1 + subject_offset,
                                    query_frame: translation_frame(
                                        plus,
                                        query_from,
                                        query_to,
                                        place.length as i32,
                                    ),
                                    hit_frame: 0,
                                    qrow: &qrow,
                                    srow: &srow,
                                    ops: Some(&hsp.edit_script),
                                    // The first residue of the row: the lowest nucleotide of the
                                    // interval on the plus strand, the highest on the minus
                                    // strand, where the row runs backwards.
                                    query_first: if plus { query_from - 1 } else { query_to - 1 },
                                    query_step: if plus { 3 } else { -3 },
                                    alignment_identities: identities,
                                })?;
                            }
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            });
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:745-760
            // A query whose HSP throws (`CaughtException`): BLAST XML leaves its `Iteration`
            // open after the `Hit`s before (each written as it is made), the XML2 and JSON
            // documents leave its report out (`Xml2Document::report`; for `-subject`, whose
            // report is held, where the next query's arrives); the document's message follows
            // the report and the texts of its loads, and the next query is written. Another
            // failed query leaves its partial bytes in `buf`, which the error return drops.
            let caught = match written {
                Ok(()) => None,
                Err(error) => match error.downcast::<CaughtException>() {
                    Ok(CaughtException(exception)) => Some(crate::ncbi_exception::error_post_line(
                        None,
                        Some(options.program().lower),
                        &document.caught_message(&exception),
                    )),
                    Err(error) => return Err(error),
                },
            };
            posts.document(document.take_posts());
            crate::report::query_warnings::write_posts_to(
                &buf,
                &posts,
                0,
                &mut *writer,
                &mut *diagnostics,
                reads.as_mut(),
            )
            .map_err(XmlRefusal::of_cut)?;
            if let Some(post) = caught {
                writer.flush()?;
                diagnostics.write_all(&post)?;
                diagnostics.flush()?;
            }
            // A text whose read no `Hit` placed comes after the query's report; one that a cut
            // report left out goes to NCBI's next load that makes its ask (`ReadCursor::cut`).
            if let Some(reads) = reads.as_mut() {
                reads.before_query(1, &mut *writer, &mut *diagnostics)?;
            }
            buf.clear();
        }
    }
    Ok(())
}

/// The end of the document (`PrintEpilog`): its trailer when `queries` queries were formatted,
/// else the epilog's end of line alone. For the XML2 and JSON documents `xml2` is the run's
/// writer: it writes the report held back by a `-subject` run and the document's end; `program`
/// names the run's program in its message (`ncbi_exception::error_post_line`).
pub(super) fn write_xml_epilog(
    writer: &mut impl crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
    program: &str,
    queries: usize,
    xml2: Option<&mut RunXml2>,
) -> Result<()> {
    let mut buf = Vec::new();
    let (finished, posts) = match xml2 {
        Some(xml2) => {
            let finished = xml2.finish(&mut buf).map_err(|error| {
                // The held report of a `-subject` run that fails (`CaughtException`): its
                // message, then the end of the document.
                match error.downcast::<CaughtException>() {
                    Ok(CaughtException(exception)) => Ok(crate::ncbi_exception::error_post_line(
                        None,
                        Some(program),
                        &xml2.caught_message(&exception),
                    )),
                    Err(error) => Err(error),
                }
            });
            (finished, xml2.take_posts())
        }
        None => {
            // The epilog's insertion, line feed and flush (`BlastXmlWriter::finish`'s cuts).
            let mut document = BlastXmlWriter::new(queries);
            let finished = document.finish(&mut buf).map_err(Err);
            (finished, document.take_posts())
        }
    };
    let post = match finished {
        Ok(()) => None,
        Err(Ok(post)) => Some(post),
        Err(Err(error)) => return Err(error),
    };
    if let Some(post) = post {
        writer.flush()?;
        diagnostics.write_all(&post)?;
        diagnostics.flush()?;
    }
    write_with_posts(&buf, posts, writer, diagnostics)
}

/// Writes `buf`, a part of the document that holds no `Hit` (a held report and the document's
/// end), with what the document's writer posted while it wrote it (`StructuredXml::take_posts`),
/// each to `diagnostics` after the bytes before it (`query_warnings::write_posts_to`).
fn write_with_posts(
    buf: &[u8],
    posts: crate::report::blast_xml::DocumentPosts,
    writer: &mut dyn crate::report::ostream::ReportWrite,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    let mut hits = crate::report::query_warnings::HitPosts::default();
    hits.document(posts);
    crate::report::query_warnings::write_posts_to(buf, &hits, 0, writer, diagnostics, None)
}

#[cfg(test)]
#[path = "blast_xml_taxonomy_error_tests.rs"]
mod taxonomy_error_tests;

#[cfg(test)]
#[path = "blast_xml_scan_tests.rs"]
mod scan_tests;

#[cfg(test)]
#[path = "blast_xml_held_side_files_tests.rs"]
mod held_side_files_tests;

#[cfg(test)]
#[path = "blast_xml_held_formats_tests.rs"]
mod held_formats_tests;
