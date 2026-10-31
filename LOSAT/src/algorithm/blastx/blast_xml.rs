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
//! ```c++
//!     } else if (m_FormatType == CFormattingArgs::eXml) {
//!         CRef<CSearchResults> res(const_cast<CSearchResults*>(&results));
//!         res->TrimSeqAlign(m_HitlistSize);
//!         m_AccumulatedResults.push_back(res);
//!         ...
//!         	CCmdLineBlastXMLReportData report_data(m_AccumulatedQueries,
//!                                                    m_AccumulatedResults,
//!                                                    *m_Options, m_DbName, m_DbIsAA,
//!                                                    m_QueryGenCode, m_DbGenCode,
//!                                                    m_IsRemoteSearch);
//!         	BlastXML_FormatReport(xml_output, &report_data, &m_Outfile,
//!                                   m_BlastXMLIncremental.GetPointer());
//! ```
//! (`PrintOneResultSet` is called once per query's result set; the accumulated lists hold that
//! one query.)
//!
//! The XML2 and JSON documents of outfmt 15 and 16, and 14 and 13 as side files with a master file
//! (`report::blast_xml2`), take the same path
//! with another writer: `write_blastx_xml` is generic over `StructuredXml`, and the run keeps
//! one `Xml2Writer` across the query batches (`xml2_writer`), as it passes `queries_before` to
//! outfmt 5's. Its document start is written with the prolog, before the first batch.
//!
//! Every refusal (`XmlRefusal`: a byte NCBI's XML stream cannot write in a title) comes before
//! the first byte of the document, though the query file is read batch by batch:
//!
//! - `prescan_queries` reads the query file once more before the search (only with outfmt 5,
//!   13 to 16) and refuses the first query whose title the writer would refuse;
//! - when some `-subject` title would be refused if its subject had a hit
//!   (`subjects_may_be_refused`), or some query's message would be (a title that is refused only
//!   in the message of a query NCBI did not search), or (with `-db`) the taxonomy database is
//!   on the search path and the lookup of the taxid of some sequence of the database would be
//!   refused (a scientific name NCBI's XML2 cannot write, a damaged record, which outfmt 5 reads
//!   in each batch's prefetch, `ReportedSubjects::prefetch_messages`: checked once before the
//!   first byte, `taxdb_may_refuse`, with the names for XML2 and JSON only), or outfmt 5's
//!   prefetch is refused for NCBI's taxonomy cache (`BlastxDatabase::prefetch_refused`, at the
//!   first batch, after the prolog), so that a later batch could refuse it,
//!   the run holds the document in memory (`HeldOutput`; with outfmt 13 and 14 the side files'
//!   bytes too, `HeldSideFiles`) and writes it at the epilog; a refusal found while the batches
//!   are reported drops it, with nothing written (an empty master file, no side file).

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;
use std::rc::Rc;
use std::sync::{Arc, Mutex, MutexGuard};

use anyhow::{ensure, Context, Result};

use super::args::ResolvedOptions;
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
use crate::report::blast_xml::subjects::SubjectNames;
use crate::report::blast_xml::{
    filter, BlastXmlWriter, QueryMask, StructuredXml, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats,
    REFERENCE_GAPPED_BLAST,
};
use crate::report::blast_xml2::framing::{side_file_open_error, side_file_write_error};
use crate::report::blast_xml2::{
    FsSideFiles, RunSideFiles, SideFiles, Xml2Document, Xml2Kind, Xml2Run, Xml2Target, Xml2Writer,
};
use crate::report::query_warnings::{EMPTY_DATABASE_MESSAGE, INVALID_QUERY_MESSAGE};
use crate::report::serial::seqalign::seqalign_evalue;

const PROGRAM: &str = "BLASTX";

/// LOSAT's refusal of a title or name that NCBI's XML stream cannot write (not an error of
/// NCBI's engine: the run ends with exit status 1 and the message, `native::report_error`
/// passes it on).
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
}

fn refuse(refusal: names::Refusal, what: &str) -> anyhow::Error {
    XmlRefusal(refusal.error(PROGRAM, what).to_string()).into()
}

/// The masked regions of a query as NCBI hands them to the XML writer
/// (`Blast_GetSeqLocInfoVector`, blast_aux.cpp:930-959): per context, in context order, the
/// context's nucleotide intervals other than an empty one and the whole query, with the
/// context's frame. LOSAT's BLASTX has no `-query_loc`, so the interval is the whole record.
fn query_masks(batch: &BatchResults, qidx: usize, length: usize) -> Vec<QueryMask> {
    let mut masks = Vec::new();
    for context in &batch.prepared.contexts[qidx * 6..qidx * 6 + 6] {
        for &(from, to) in &context.dna_masks {
            if from > to || (from == 0 && to == length as i32 - 1) {
                continue;
            }
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
        kappa: karlin.k,
        lambda: karlin.lambda,
        entropy: karlin.h,
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
    oid: usize,
) -> Result<SubjectNames> {
    match shown {
        ShownSubjects::Fasta { .. } => SubjectNames::of_record(
            subject.internal_id.as_bytes(),
            subject.title.as_bytes(),
            true,
        )
        .map_err(|refusal| refuse(refusal, "a subject title")),
        ShownSubjects::Database { reported, .. } => match reported.parsed_subject(oid) {
            Some(parsed) => SubjectNames::of_parsed(parsed, reported.parsed_taxonomy(oid)),
            None => SubjectNames::of_database_sequence(
                reported.ordinal(oid),
                &reported.stored_titles()[oid],
                // XML2 and JSON alone show the taxid (`ReportedSubjects` reads it for them only).
                reported.taxid(oid),
                reported.taxonomy(oid),
            ),
        }
        .map_err(|refusal| refuse(refusal, "a database title")),
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
        qi == raw.q_end as usize && si == raw.s_end as usize,
        "BLASTX edit script does not reach final endpoints"
    );
    let (qstart, qend) = query_endpoints(query.sequence.len(), raw.frame, raw.q_start, raw.q_end);
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
}

impl HeldOutput {
    pub(super) fn new(inner: Box<dyn Write>, side: Option<HeldSideFiles>) -> Self {
        HeldOutput(Rc::new(RefCell::new(HeldState {
            inner,
            holding: false,
            held: Vec::new(),
            discarded: false,
            side,
        })))
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
                return Err(error);
            }
        }
        let held = std::mem::take(&mut state.held);
        state.inner.write_all(&held)?;
        state.inner.flush()?;
        Ok(())
    }
}

impl Write for HeldOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut state = self.0.borrow_mut();
        if state.discarded {
            Ok(bytes.len())
        } else if state.holding {
            state.held.extend_from_slice(bytes);
            Ok(bytes.len())
        } else {
            state.inner.write(bytes)
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

impl Drop for HeldState {
    fn drop(&mut self) {
        if !self.discarded && !self.held.is_empty() {
            let _ = self.inner.write_all(&self.held);
            let _ = self.inner.flush();
        }
    }
}

/// The side files of an outfmt 13 or 14 run (`FsSideFiles`), which a hold keeps in memory: the
/// n-th report is written to its file at once, or, while the run holds its output, at `release`
/// (the epilog) in the order of the reports. A refusal found while the batches are reported
/// drops them, so that the run leaves no side file (design OF5 4.4); any other error ends the
/// run with the side files of the earlier reports written (`SideState::drop`), as NCBI leaves
/// them.
#[derive(Clone)]
pub(super) struct HeldSideFiles(Arc<Mutex<SideState>>);

struct SideState {
    files: FsSideFiles,
    holding: bool,
    discarded: bool,
    /// The reports held, as `(n, bytes)` in the order they were made.
    held: Vec<(u32, Vec<u8>)>,
}

impl HeldSideFiles {
    pub(super) fn new(files: FsSideFiles) -> Self {
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

    fn hold(&self) {
        self.state().holding = true;
    }

    fn discard(&self) {
        let mut state = self.state();
        state.discarded = true;
        state.held.clear();
    }

    /// Writes the held side files, each to a new file (truncating one of that name).
    fn release(&self) -> Result<()> {
        let mut state = self.state();
        state.holding = false;
        state.write_held()
    }
}

impl SideState {
    /// A file that cannot be created is NCBI's "Cannot open output file"; a write to one that was
    /// opened is the side-file write error of every program (`side_file_write_error`).
    fn write_held(&mut self) -> Result<()> {
        let held = std::mem::take(&mut self.held);
        for (n, bytes) in held {
            let mut file =
                std::fs::File::create(self.files.path(n)).map_err(|_| side_file_open_error())?;
            file.write_all(&bytes)
                .map_err(|_| side_file_write_error())?;
        }
        Ok(())
    }
}

/// Opens the side file of the n-th report as NCBI does when it writes the report, without
/// leaving a mark: a path that exists is opened for writing and not truncated (a directory or a
/// file without write permission fails), a path that does not exist is created and removed. The
/// run that holds its side files makes this check when the report is made, so that an
/// unopenable file stops it where NCBI's run stops, not at the epilog after every later batch
/// was searched (DECISIONS 35 S1). The bytes are written at `release`, which creates the file
/// again.
fn check_side_file_openable(path: &Path) -> io::Result<()> {
    match std::fs::OpenOptions::new().write(true).open(path) {
        Ok(_) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            // A dangling link is created through (`File::create` follows it); only a path that
            // did not exist at all is removed again.
            let existed = path.symlink_metadata().is_ok();
            std::fs::File::create(path)?;
            if !existed {
                std::fs::remove_file(path)?;
            }
            Ok(())
        }
        Err(error) => Err(error),
    }
}

impl Drop for SideState {
    fn drop(&mut self) {
        if !self.discarded {
            let _ = self.write_held();
        }
    }
}

/// A report held for its side file: it joins the held reports when the document is done with it.
struct HeldReport {
    side: HeldSideFiles,
    n: u32,
    bytes: Vec<u8>,
}

impl Write for HeldReport {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Drop for HeldReport {
    fn drop(&mut self) {
        let mut state = self.side.state();
        if !state.discarded {
            state.held.push((self.n, std::mem::take(&mut self.bytes)));
        }
    }
}

impl SideFiles for HeldSideFiles {
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>> {
        let state = self.state();
        if state.discarded {
            return Ok(Box::new(io::sink()));
        }
        if state.holding {
            check_side_file_openable(&state.files.path(n))?;
            return Ok(Box::new(HeldReport {
                side: self.clone(),
                n,
                bytes: Vec::new(),
            }));
        }
        Ok(Box::new(io::BufWriter::new(std::fs::File::create(
            state.files.path(n),
        )?)))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        self.state().files.listed(n)
    }
}

/// The `QueryText` of a query, refused as the writer refuses it. `search_skipped`: NCBI did not
/// search the query's batch (its message then carries the query's id, made of its title).
fn query_text(
    query: &FastaRecord,
    search_skipped: bool,
    filtered_empty: bool,
    has_hits: bool,
) -> Result<QueryText> {
    names::query_text(
        PROGRAM,
        query.warning_id(),
        query.believed(),
        query.title.as_bytes(),
        false,
        &query_messages(query, search_skipped, filtered_empty),
        has_hits,
    )
    .map_err(|error| anyhow::Error::new(XmlRefusal(error.to_string())))
}

/// Reads the query file once before the search and refuses the first query the writer would
/// refuse (a title with a byte NCBI's XML stream cannot write), so that no refusal can follow
/// a written byte. Returns whether some query would be refused if NCBI had not searched its
/// batch (a title that only its message cannot carry); the run then holds the document.
///
/// The reader is the search's own (same records and titles). A record it refuses for every
/// format (a title that is not UTF-8) ends the scan: the run reads the file again and fails
/// there as it always did.
pub(super) fn prescan_queries(file: std::fs::File, lowercase: bool) -> Result<bool> {
    let mut stream = super::input::FastaStream::from_file(file);
    let mut may_be_refused = false;
    let scanned = super::input::parse_fasta_stream_with_warnings(
        &mut stream,
        false,
        lowercase,
        super::query_setup::BATCH_SIZE,
        &mut |batch| {
            for query in batch {
                query_text(query, false, false, true)?;
                may_be_refused |= query_text(query, true, false, true).is_err();
            }
            Ok(())
        },
        &mut |_| Ok(()),
    );
    match scanned {
        Err(error) if error.is::<XmlRefusal>() => Err(error),
        _ => Ok(may_be_refused),
    }
}

/// Whether the title of some `-subject` would be refused if the subject had a hit.
pub(super) fn subjects_may_be_refused(subjects: &[FastaRecord]) -> bool {
    subjects.iter().any(|subject| {
        names::xml_subject_names(
            subject.internal_id.as_bytes(),
            subject.title.as_bytes(),
            true,
        )
        .is_err()
    })
}

/// Refuses a `-db` name that NCBI's XML stream cannot write (a byte outside US-ASCII); the
/// run checks it before the search, as the Seq-annot formats check theirs.
pub(super) fn check_database_name(name: &str) -> Result<()> {
    names::check_database_name(PROGRAM, name).map_err(|error| XmlRefusal(error.to_string()).into())
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
/// ```c++
///     // Composition-based statistics
///     if (data->GetCompositionBasedStats() > 0) {
///         params.SetCbs(data->GetCompositionBasedStats());
///     }
///     // Query genetic code
///     if (data->GetQueryGeneticCode() > 0) {
///         params.SetQuery_gencode(data->GetQueryGeneticCode());
///     }
/// ```
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
        // `-ungapped` writes no gap costs; it needs `-comp_based_stats 0`, so no `cbs` either.
        gapped: options.gapped,
        cbs: i32::from(options.composition),
        query_gencode: Some(i32::from(options.query_gencode)),
        db_gencode: None,
    };
    let delayed = matches!(params.target, Xml2Target::Subjects(_));
    Xml2Writer::with_document(Xml2Document::with_files(kind, delayed, files), params)
}

/// Writes the documents of every query of `batches` (one query batch) through `document` (BLAST
/// XML, or the XML2 / JSON documents) to `writer`: the document's header before the first query
/// of the run, each query's Iteration or report, and nothing else; the end of the document is
/// `write_epilog`'s. Everything NCBI cannot write is refused before this call's first byte: the
/// query titles and messages, the titles of the subjects with hits (and, with a hold, the whole
/// run's bytes stay in memory, `HeldOutput`). NCBI's messages about the taxonomy database go to
/// `diagnostics` at the first `Hit` that needs them.
pub(super) fn write_blastx_xml<W: StructuredXml>(
    writer: &mut impl Write,
    diagnostics: &mut dyn Write,
    document: &mut W,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
) -> Result<()> {
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    let gapped_karlin = crate::stats::tables::lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let positives = FormattingMatrix::new(&options.matrix)
        .with_context(|| format!("no formatting matrix for {}", options.matrix))?;
    let (db_name, db_num, db_len) = match shown {
        ShownSubjects::Fasta { .. } => (String::new(), 0, 0),
        ShownSubjects::Database { name, pairwise, .. } => (
            name.to_string(),
            i32::try_from(pairwise.num_sequences)?,
            i64::try_from(pairwise.total_letters)?,
        ),
    };
    check_database_name(&db_name)?;
    let filter = filter::filter_string(
        None,
        options
            .seg
            .enabled
            .then_some((options.seg.window, options.seg.locut, options.seg.hicut)),
        None,
        options.soft_masking,
    );
    let run = XmlRun {
        program: "blastx",
        version: "BLASTX 2.17.0+",
        reference: REFERENCE_GAPPED_BLAST,
        db_name: db_name.as_bytes(),
        matrix: Some(&options.matrix),
        expect: options.evalue,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
        reward_penalty: None,
        filter: &filter,
        db_num,
        db_len,
        positives: Some(&positives),
        mask_char: b'X',
    };
    // The definition line and message of every query, and the names of the subjects with hits,
    // made (and refused) before the first byte.
    let mut texts: Vec<QueryText> = Vec::new();
    let mut subject_cache: HashMap<usize, SubjectNames> = HashMap::new();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            let has_hits = lists.iter().any(|list| reported(list, options.gapped));
            texts.push(query_text(
                query,
                batch.search_skipped,
                options.filtered_empty,
                has_hits,
            )?);
            for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                if let std::collections::hash_map::Entry::Vacant(slot) =
                    subject_cache.entry(list.oid)
                {
                    slot.insert(subject_names(shown, &subjects[list.oid], list.oid)?);
                }
            }
        }
    }
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    let no_notes = TaxdbNotes::default();
    let mut taxdb_printer = match shown {
        ShownSubjects::Database { reported, .. } => reported.taxdb_notes(),
        ShownSubjects::Fasta { .. } => &no_notes,
    }
    .printer();
    let mut buf: Vec<u8> = Vec::new();
    // The rows of one HSP at a time, reused.
    let (mut qrow, mut srow): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut texts = texts.iter();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            let text = texts.next().context("one text per query")?;
            let length = query.sequence.len();
            let masks = query_masks(batch, qidx, length);
            let xml_query = XmlQuery {
                id: query.internal_id.as_bytes(),
                title: &text.title,
                len: i32::try_from(length)?,
                stats: query_stats(batch, qidx, options.gapped, &gapped_karlin),
                message: text.message.as_deref(),
                // Read by XML2 and JSON alone: the query of a batch NCBI did not search.
                searched: !batch.search_skipped,
                masks: &masks,
            };
            let written = document.query(&run, &xml_query, &mut buf, |sink| {
                for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                    let subject = &subjects[list.oid];
                    let names = subject_cache
                        .get(&list.oid)
                        .context("names of a subject with hits")?;
                    let hit = XmlHit::of(names, i32::try_from(subject.sequence.len())?);
                    // NCBI's messages about the taxonomy database: opening and the volume's
                    // lookup at the first `Hit` that needs them, the lookup of the `Hit` every
                    // time.
                    taxdb_printer.at_names(list.oid, names, &mut *diagnostics)?;
                    sink.hit(&hit, |hsps| {
                        for hsp in list
                            .hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                        {
                            // The Seq-align's endpoints, lower first: the minus frame's
                            // outfmt 6 order is not the XML's.
                            let (query_from, query_to) =
                                xml_rows(hsp, subject, query, &code, &mut qrow, &mut srow)?;
                            let plus = hsp.hsp.frame > 0;
                            hsps.hsp(&XmlHsp {
                                score: hsp.hsp.score,
                                bit_score: hsp.bit_score,
                                evalue: seqalign_evalue(hsp.evalue),
                                query_from,
                                query_to,
                                hit_from: hsp.hsp.s_start + 1,
                                hit_to: hsp.hsp.s_end,
                                query_frame: translation_frame(
                                    plus,
                                    query_from,
                                    query_to,
                                    length as i32,
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
                            })?;
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            });
            // A failed query leaves its partial bytes in `buf`, which the error return drops.
            written?;
            writer.write_all(&buf)?;
            buf.clear();
        }
    }
    Ok(())
}

/// The end of the document (`PrintEpilog`): its trailer when `queries` queries were formatted,
/// else the epilog's end of line alone. For the XML2 and JSON documents `xml2` is the run's
/// writer: it writes the report held back by a `-subject` run and the document's end.
pub(super) fn write_xml_epilog(
    writer: &mut impl Write,
    queries: usize,
    xml2: Option<&mut RunXml2>,
) -> Result<()> {
    let mut buf = Vec::new();
    match xml2 {
        Some(xml2) => xml2.finish(&mut buf)?,
        None => BlastXmlWriter::new(queries).finish(&mut buf)?,
    }
    writer.write_all(&buf)?;
    Ok(())
}

#[cfg(test)]
mod taxonomy_error_tests {
    //! `XmlRefusal::of_taxonomy`: only a refusal of the taxonomy lookups drops a held document;
    //! a read error passes through, and the run writes what it holds (`HeldState::drop`).
    use super::*;
    use crate::cli::NativeError;
    use crate::seqdb::SeqDbError;
    use std::path::PathBuf;

    fn through(error: SeqDbError) -> anyhow::Error {
        XmlRefusal::of_taxonomy(crate::algorithm::common::database::taxdb_error_for_tests(
            error,
        ))
    }

    #[test]
    fn a_refused_record_or_name_is_a_refusal_with_its_message() {
        let error = through(SeqDbError::Unsupported {
            path: PathBuf::from("taxdb.btd"),
            message: "the taxonomy record of taxid 9606 reaches past the mapped end".into(),
        });
        assert!(error.is::<XmlRefusal>());
        assert!(error.to_string().contains("taxid 9606"), "{error}");
    }

    #[test]
    fn a_read_error_is_not_a_refusal() {
        let io = SeqDbError::Io {
            path: PathBuf::from("taxdb.btd"),
            source: std::io::Error::other("Input/output error"),
        };
        let error = through(io);
        assert!(!error.is::<XmlRefusal>(), "{error}");
        let format = SeqDbError::Format {
            path: PathBuf::from("db.phr"),
            message: "byte range 0..9 passes the end of the file (4 bytes)".into(),
        };
        let error = through(format);
        assert!(!error.is::<XmlRefusal>(), "{error}");
        // An error of NCBI's own (exit status 2) is passed on whole.
        let native = anyhow::Error::new(NativeError {
            exit: 2,
            message: "BLAST Database error: x\n".into(),
        });
        assert!(XmlRefusal::of_taxonomy(native).is::<NativeError>());
    }
}

#[cfg(test)]
mod held_side_files_tests {
    //! The hold of the side files of outfmt 13 and 14 (`HeldSideFiles`, `HeldOutput`): nothing on
    //! the file system while the run holds its output, the files in report order at `release`, none
    //! after `discard`, and the files made so far when the run ends with another error.
    use super::*;
    use crate::cli::NativeError;
    use crate::report::blast_xml2::Xml2Kind;
    use std::path::{Path, PathBuf};

    fn work(stem: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("losat_held_side_{stem}_{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// A held output over `<dir>/o.out` and its side files `<dir>/o.out_<n>.xml`.
    fn run(dir: &Path) -> (HeldOutput, HeldSideFiles) {
        let master = std::fs::File::create(dir.join("o.out")).unwrap();
        let side = HeldSideFiles::new(FsSideFiles::new(&dir.join("o.out"), Xml2Kind::Files14));
        (HeldOutput::new(Box::new(master), Some(side.clone())), side)
    }

    fn report(side: &mut HeldSideFiles, n: u32, bytes: &[u8]) {
        let mut file = side.create(n).unwrap();
        file.write_all(bytes).unwrap();
        file.flush().unwrap();
    }

    fn names(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn without_a_hold_the_files_are_written_at_once() {
        let dir = work("nohold");
        let (held, mut side) = run(&dir);
        report(&mut side, 1, b"one");
        assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one");
        held.release().unwrap();
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn a_hold_keeps_the_files_in_memory_until_the_release_and_writes_them_in_order() {
        let dir = work("hold");
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"one");
        report(&mut side, 2, b"two");
        let mut master = held.clone();
        master.write_all(b"master").unwrap();
        assert_eq!(names(&dir), ["o.out"]);
        assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"");
        held.release().unwrap();
        assert_eq!(names(&dir), ["o.out", "o.out_1.xml", "o.out_2.xml"]);
        assert_eq!(std::fs::read(dir.join("o.out_2.xml")).unwrap(), b"two");
        assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"master");
        // after the release the files are written at once again
        report(&mut side, 3, b"three");
        assert_eq!(std::fs::read(dir.join("o.out_3.xml")).unwrap(), b"three");
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn a_discard_leaves_no_side_file_and_an_empty_master() {
        let dir = work("discard");
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"one");
        let mut master = held.clone();
        master.write_all(b"master").unwrap();
        held.discard();
        report(&mut side, 2, b"two");
        drop((held, side, master));
        assert_eq!(names(&dir), ["o.out"]);
        assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"");
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn any_other_end_of_the_run_writes_the_files_made_so_far() {
        let dir = work("dropped");
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"one");
        drop((held, side));
        assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one");
        std::fs::remove_dir_all(dir).ok();
    }

    /// DECISIONS 35 S1: the run that holds its side files opens each one when its report is made, as
    /// NCBI does: a file that cannot be opened stops the run there (the earlier reports are written
    /// by the drop, the master stays empty), and the check leaves no file of its own.
    #[test]
    fn a_side_file_that_cannot_be_opened_stops_at_its_report_while_holding() {
        let dir = work("blocked");
        std::fs::create_dir(dir.join("o.out_2.xml")).unwrap();
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"one");
        // nothing but the master file yet, and the check of the first report left no mark
        assert_eq!(names(&dir), ["o.out", "o.out_2.xml"]);
        let error = side.create(2).err().expect("a directory cannot be opened");
        assert!(error.raw_os_error().is_some(), "{error:?}");
        let mut master = held.clone();
        master.write_all(b"master").unwrap();
        drop((held, side, master));
        // the file before the blocked one stays; the master has nothing of the held bytes of the
        // failed report
        assert_eq!(names(&dir), ["o.out", "o.out_1.xml", "o.out_2.xml"]);
        assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one");
        std::fs::remove_dir_all(dir).ok();
    }

    /// The check opens an existing file for writing without truncating it, and a refusal after it
    /// (`discard`) leaves the file as it was.
    #[test]
    fn the_check_of_an_existing_side_file_does_not_truncate_it() {
        let dir = work("existing");
        std::fs::write(dir.join("o.out_1.xml"), b"old").unwrap();
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"new");
        assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"old");
        held.discard();
        drop((held, side));
        assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"old");
        std::fs::remove_dir_all(dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn a_read_only_side_file_stops_at_its_report_while_holding() {
        use std::os::unix::fs::PermissionsExt;
        let dir = work("readonly");
        std::fs::write(dir.join("o.out_1.xml"), b"old").unwrap();
        std::fs::set_permissions(
            dir.join("o.out_1.xml"),
            std::fs::Permissions::from_mode(0o444),
        )
        .unwrap();
        let (held, mut side) = run(&dir);
        held.hold();
        // a process with the right to write anyway (root) opens it, as it opens it for NCBI
        let writable = std::fs::OpenOptions::new()
            .write(true)
            .open(dir.join("o.out_1.xml"))
            .is_ok();
        assert_eq!(side.create(1).is_ok(), writable);
        drop((held, side));
        std::fs::remove_dir_all(dir).ok();
    }

    /// A write that fails after the open (a link to `/dev/full`) is the write error of every
    /// program, not NCBI's "Cannot open output file" (DECISIONS 35 S2).
    #[cfg(target_os = "linux")]
    #[test]
    fn a_side_file_write_error_at_the_release_is_the_write_error() {
        let dir = work("full");
        std::os::unix::fs::symlink("/dev/full", dir.join("o.out_1.xml")).unwrap();
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"one");
        let mut master = held.clone();
        master.write_all(b"master").unwrap();
        let error = held.release().unwrap_err();
        let native = error.downcast_ref::<NativeError>().expect("a NativeError");
        assert_eq!(native.exit, 6);
        assert_eq!(native.message, "BLAST failed to write output\n");
        drop((held, side, master));
        assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"");
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn a_side_file_that_cannot_be_created_at_the_release_is_the_open_error() {
        // The check at the report passed; a directory put in its place before the release is
        // what `File::create` then refuses.
        let dir = work("late");
        let (held, mut side) = run(&dir);
        held.hold();
        report(&mut side, 1, b"one");
        std::fs::create_dir(dir.join("o.out_1.xml")).unwrap();
        let error = held.release().unwrap_err();
        let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
        assert_eq!(native.exit, 1);
        assert_eq!(
            native.message,
            "Command line argument error: Cannot open output file\n"
        );
        drop((held, side));
        std::fs::remove_dir_all(dir).ok();
    }
}
