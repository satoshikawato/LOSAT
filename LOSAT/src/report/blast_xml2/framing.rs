//! The framings of an XML2 / JSON document: where a query's report goes and what comes before and
//! after the reports (DESIGN_OF5 1.1, 1.2).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2332-2356
//! There is one document per query and every document is complete (the report header, the
//! parameters, one `Search`). Outfmt 16 writes them one after the other after the XML2 root's
//! start tag, outfmt 15 as the elements of one JSON array, outfmts 14 and 13 each to a file of
//! its own (`<base>_<n>.xml`, `<base>_<n>.json`) with a master file at `-out` (OF5c).
//!
//! `-subject` searches keep the last query's result back until the next query's arrives
//! (`m_IsBl2Seq`, the accumulated results of blast_format.cpp:2363-2404) and write it in the
//! epilog. NCBI's epilog is not reached after a failed query batch, so that held report is
//! never written: the document of a failed run lacks the last formatted query's report. The
//! programs do not make that report (`Xml2Document::unwritten`) and refuse none of its texts
//! (`blast_xml::written_queries`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2385, :2383-2389
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2206-2219
//!
//! A result whose query's Seq-id matches that of the held report's first result joins the held
//! report (`XmlQuery::joins_held`): NCBI accumulates it, and writes the results of one Seq-id as
//! one report whose `bl2seq` list has a `Search` for each (`Xml2Writer`'s join, `HeldSearches`).
//!
//! A report that LOSAT refuses (a title that NCBI's definition line generator reads past the end of,
//! `blast_xml::names::Refusal`) is not made, and is refused where NCBI makes it
//! (`Xml2Document::hold_refused`): a `-db` report at its arrival, after the reports before it; a
//! held `-subject` report where the document would write it (the arrival of another query's result,
//! in the same report call or a later one of the run, or the epilog), since NCBI makes that report,
//! and crashes on it, only when it writes it, after the reports before it (each flushed when it is
//! written); a run that an error stops first never makes it (it is dropped with the report). The
//! side file of a refused report of outfmt 13 or 14 is left empty, as NCBI opens it before it makes
//! the report.

use std::io::{self, Write};
use std::path::Path;

use anyhow::Result;

use super::files::{run_side_files, SideFileSink};
use crate::report::blast_xml::writer::{DocumentPost, DocumentPosts};
use crate::report::serial::xml::{XmlThrow, XML2_ROOT_ATTRIBUTES};

/// What a report's writer hands the framing with the report's bytes (`Xml2Document::report`):
/// what NCBI posts while it writes the report (`DocumentPost`, at offsets of the report), and the
/// exception of a 14 side file, which ends the run when the report is written (blastxml2_format.cpp:
/// 768-778 has no `catch`; outfmt 16's report catches its exception and posts it,
/// `report::blast_xml2::writer`). Dropping it loses them, so it is `must_use`.
#[derive(Clone, Debug, Default)]
#[must_use = "a report's posts and the side file's exception go to the framing (Xml2Document::report)"]
pub struct MadeReport {
    pub posts: Vec<DocumentPost>,
    pub ends_run: Option<Box<XmlThrow>>,
}

/// Which of NCBI's four document forms a run writes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Xml2Kind {
    /// Outfmt 16 (`eXml2_S`): all reports on the output stream.
    Single16,
    /// Outfmt 15 (`eJson_S`): all reports in one JSON array on the output stream.
    Single15,
    /// Outfmt 14 (`eXml2`): one `<base>_<n>.xml` file per report and a master file.
    Files14,
    /// Outfmt 13 (`eJson`): one `<base>_<n>.json` file per report and a master file.
    Files13,
}

impl Xml2Kind {
    /// The document form of a `-outfmt` format (`None` for the formats that are not XML2 or JSON
    /// documents of `CBlastOutput2`).
    pub fn of(format: crate::blastinput::app::ReportFormat) -> Option<Self> {
        use crate::blastinput::app::ReportFormat;
        match format {
            ReportFormat::Xml2S => Some(Self::Single16),
            ReportFormat::JsonS => Some(Self::Single15),
            ReportFormat::Xml2 => Some(Self::Files14),
            ReportFormat::Json => Some(Self::Files13),
            _ => None,
        }
    }

    /// The `-outfmt` format of the document form (the inverse of `of`).
    pub fn report_format(self) -> crate::blastinput::app::ReportFormat {
        use crate::blastinput::app::ReportFormat;
        match self {
            Self::Single16 => ReportFormat::Xml2S,
            Self::Single15 => ReportFormat::JsonS,
            Self::Files14 => ReportFormat::Xml2,
            Self::Files13 => ReportFormat::Json,
        }
    }

    /// Whether the document is written as side files with a master file (outfmt 13 and 14).
    pub fn is_files(self) -> bool {
        matches!(self, Self::Files14 | Self::Files13)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:412-417, 750-760
    // (cited at `report::db_subject::HitLoads`: the `Hit` whose lookup throws, and the report
    // function that catches it)
    /// What the document does after a `Hit` whose subject's lookup throws
    /// (`report::db_subject::AfterThrow`): the stream forms (15 and 16) leave the query's report
    /// out with the exception's post, the forms of side files (13 and 14) end the run with it.
    pub fn after_throw(self) -> crate::report::db_subject::AfterThrow {
        if self.is_files() {
            crate::report::db_subject::AfterThrow::EndsRun
        } else {
            crate::report::db_subject::AfterThrow::EndsReport
        }
    }
}

/// A side file NCBI cannot open: `CArgException` "Cannot open output file", which the application
/// reports as "Command line argument error: ..." with exit status 1 (see `Xml2Document::write`).
pub fn side_file_open_error() -> anyhow::Error {
    crate::cli::NativeError {
        exit: 1,
        message: "Command line argument error: Cannot open output file\n".into(),
    }
    .into()
}

/// A write to a side file that fails after it was opened (a full disk) where no stream of NCBI's
/// is modelled: a caller's receiver that fails, and a library call's side files
/// (`ostream::write_side_file` writes the command line's as NCBI does: an 802.5 post and exit
/// status 0). The run stops with the message and the exit status of a failed write to `-out` in
/// outfmt 0 (`BLAST failed to write output`, exit status 6; the side files written before stay,
/// the master file stays empty).
pub fn side_file_write_error() -> anyhow::Error {
    crate::cli::NativeError::write_failure().into()
}

/// The files of a 13 / 14 run: `create(n)` opens the file of the n-th report (`<base>_<n>.xml` or
/// `<base>_<n>.json`), and `listed(n)` is the name the master file gives for it (the two differ:
/// the master drops the directory and always uses the XML rule for the base, `files::base_name`).
pub trait SideFiles {
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>>;
    fn listed(&self, n: u32) -> Vec<u8>;

    /// Writes the n-th report, whole, as its side file: opens it (`create`; a failure is NCBI's
    /// "Cannot open output file", `side_file_open_error`) and writes it (a failure is
    /// `side_file_write_error`). A caller's receiver takes the bytes as they are
    /// (`files::ReceivedSideFiles`).
    fn write_file(&mut self, n: u32, report: &[u8]) -> Result<()> {
        let mut file = self.create(n).map_err(|_| side_file_open_error())?;
        file.write_all(report)
            .and_then(|()| file.flush())
            .map_err(|_| side_file_write_error())
    }

    /// Whether the n-th side file can be opened, asked without leaving a mark by a run that holds
    /// its side files in memory, where NCBI opens the file (`FsSideFiles::check_openable`). Side
    /// files that are no files of their own (a caller's receiver) can always be.
    fn check_openable(&self, _n: u32) -> io::Result<()> {
        Ok(())
    }
}

/// The side files of a run that writes none (outfmt 15 and 16), and of the callers that have
/// neither a file name nor a receiver for them (`files::run_side_files`).
pub struct NoSideFiles;

impl SideFiles for NoSideFiles {
    fn create(&mut self, _n: u32) -> io::Result<Box<dyn Write + '_>> {
        Err(io::Error::other(
            "outfmt 13 and 14 need a file name for the side files",
        ))
    }

    fn listed(&self, _n: u32) -> Vec<u8> {
        Vec::new()
    }
}

/// The side files a run chooses at run time (`Xml2Document::for_run`).
pub type RunSideFiles = Box<dyn SideFiles + Send>;

impl SideFiles for RunSideFiles {
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>> {
        (**self).create(n)
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        (**self).listed(n)
    }

    fn write_file(&mut self, n: u32, report: &[u8]) -> Result<()> {
        (**self).write_file(n, report)
    }

    fn check_openable(&self, n: u32) -> io::Result<()> {
        (**self).check_openable(n)
    }
}

/// The start of an outfmt 16 document (`BlastXML2_PrintHeader`, blastxml2_format.cpp:781-802):
/// the declaration and the `BlastXML2` root's start tag with its namespace attributes, which ends
/// the line before `>`.
pub fn header16() -> String {
    format!("<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastXML2{XML2_ROOT_ATTRIBUTES}>")
}

/// The start of an outfmt 15 document (blastxml2_format.cpp:834-837).
pub const HEAD15: &[u8] = b"{\n\"BlastOutput2\": [\n";

/// The end of an outfmt 16 document (`PrintEpilog`, blast_format.cpp:2419-2423).
pub const TAIL16: &[u8] = b"</BlastXML2>\n";

/// The end of an outfmt 15 document (blast_format.cpp:2440-2445).
pub const TAIL15: &[u8] = b"]\n}\n";

/// The start of the master file of outfmt 14 (`x_GenerateXML2MasterFile`, blast_format.cpp:2425-2437):
/// one `<xi:include href="<file>"/>` line follows for each side file, then `</BlastXML2>`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2425-2437
pub const MASTER14_HEAD: &[u8] = b"<?xml version=\"1.0\"?>\n<BlastXML2\nxmlns=\"http://www.ncbi.nlm.nih.gov\"\nxmlns:xi=\"http://www.w3.org/2003/XInclude\"\nxmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\nxs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\">\n";

/// What `PrintProlog` writes for a document of `kind` (nothing for the side-file forms).
pub fn prolog_bytes(kind: Xml2Kind) -> Vec<u8> {
    match kind {
        Xml2Kind::Single16 => header16().into_bytes(),
        Xml2Kind::Single15 => HEAD15.to_vec(),
        Xml2Kind::Files14 | Xml2Kind::Files13 => Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xml2format.cpp:95-110
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:579-612
/// What the held report of a `-subject` run keeps for the results that join it
/// (`XmlQuery::joins_held`, `Xml2Writer`): where its last `Search` ends, the id, title and length
/// that its first `Search` shows and the masks that it lists and applies (NCBI's report shows
/// every `Search` with its query `m_Query` and the masks of its first result), and whether its
/// stream threw (no later `Search` is written). The buffers are reused from report to report.
#[derive(Clone, Debug, Default)]
pub struct HeldSearches {
    pub(super) end: usize,
    pub(super) id: Vec<u8>,
    pub(super) title: Vec<u8>,
    pub(super) len: i32,
    pub(super) masks: Vec<crate::report::blast_xml::QueryMask>,
    pub(super) cut: bool,
}

/// The held report of a `-subject` run as a result that joins it changes it
/// (`Xml2Document::held_report`): its bytes, what its writer handed with them, the exception of
/// its making, and its `Search`es.
pub(super) struct HeldReport<'d> {
    pub(super) report: &'d mut Vec<u8>,
    pub(super) made: &'d mut MadeReport,
    pub(super) failure: &'d mut Option<crate::ncbi_exception::NcbiException>,
    pub(super) searches: &'d mut HeldSearches,
}

/// The document of one run: the reports of the formatted queries and the state that the
/// framing needs (the count of reports, the held report of a `-subject` run).
pub struct Xml2Document<F: SideFiles = NoSideFiles> {
    kind: Xml2Kind,
    /// `m_XMLFileCount`: the reports written so far.
    written: u32,
    /// Whether the document's start has been written. NCBI writes it before the search
    /// (`PrintProlog`); LOSAT where the program writes it (`start`: BLASTX before the search, the
    /// others' reports before the warnings of their first query, `StructuredXml::start`), else
    /// with the first report or the epilog.
    started: bool,
    /// A `-subject` run (`m_IsBl2Seq`): a report is written when the next one arrives.
    delayed: bool,
    pending: Vec<u8>,
    pending_made: MadeReport,
    has_pending: bool,
    /// The held report of a `-subject` run failed with this exception (`CaughtException`): NCBI
    /// makes it, and fails, when it writes it.
    pending_failure: Option<crate::ncbi_exception::NcbiException>,
    /// The `Search`es of the held report (`HeldSearches`).
    held: HeldSearches,
    /// LOSAT's refusal of the held report (`hold_refused`), returned where it would be written.
    held_refusal: Option<anyhow::Error>,
    files: F,
    /// The posts of the reports written to the output buffers since `take_posts`, and the
    /// exception of a 14 side file that ends the run after them.
    posts: DocumentPosts,
}

impl Xml2Document<NoSideFiles> {
    /// A document of kind `kind` on the output stream; `delayed` for a `-subject` search. The
    /// unit tests make documents this way; a run uses `for_run`.
    #[cfg(test)]
    pub fn new(kind: Xml2Kind, delayed: bool) -> Self {
        Self::with_files(kind, delayed, NoSideFiles)
    }
}

impl Xml2Document<RunSideFiles> {
    /// The document of a run with `-out` `out`: the side files of outfmt 13 and 14 are written
    /// next to it (`FsSideFiles`), or handed to the caller's receiver when the format's output
    /// gives one (`sink`, `files::run_side_files`); outfmt 15 and 16 have none. A caller with
    /// outfmt 13 or 14 and neither gets an error when the first report is written.
    pub fn for_run(
        kind: Xml2Kind,
        delayed: bool,
        out: Option<&Path>,
        sink: Option<&SideFileSink>,
    ) -> Self {
        Self::with_files(kind, delayed, run_side_files(kind, out, sink))
    }
}

/// The framing state that one report call hands to the next when a run's reports are written by
/// several calls (the BLASTP block search, `x_query_blocks.rs`; BLASTN's report stream,
/// `blastn::blast_engine::x_rpt_stream`): the count of reports written, whether the start of the
/// document has been written, and the held report of a `-subject` run (with its refusal,
/// `Xml2Document::hold_refused`). (BLASTN's report stream hands the sizes of the writer's buffers
/// on beside it: `report::blast_xml::BufferSizes`.)
#[derive(Default)]
pub struct Xml2Carry {
    written: u32,
    started: bool,
    pending: Vec<u8>,
    pending_made: MadeReport,
    has_pending: bool,
    pending_failure: Option<crate::ncbi_exception::NcbiException>,
    held: HeldSearches,
    held_refusal: Option<anyhow::Error>,
}

impl Xml2Carry {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:360-367
    /// The start of a document of kind `kind` where the run stops with an error of NCBI's before any
    /// call wrote it (NCBI's exit flushes the start that `PrintProlog` wrote): appended to `out`
    /// once, nothing where it was written.
    pub fn start_at_stop(&mut self, kind: Xml2Kind, out: &mut Vec<u8>) {
        if !std::mem::replace(&mut self.started, true) {
            out.extend_from_slice(&prolog_bytes(kind));
        }
    }
}

impl<F: SideFiles> Xml2Document<F> {
    /// Takes the state of the run so far, leaving the document as new. A later call's document
    /// takes it up with [`Xml2Document::resume`].
    pub fn take_carry(&mut self) -> Xml2Carry {
        Xml2Carry {
            written: std::mem::take(&mut self.written),
            started: std::mem::take(&mut self.started),
            pending: std::mem::take(&mut self.pending),
            pending_made: std::mem::take(&mut self.pending_made),
            has_pending: std::mem::take(&mut self.has_pending),
            pending_failure: self.pending_failure.take(),
            held: std::mem::take(&mut self.held),
            held_refusal: self.held_refusal.take(),
        }
    }

    /// Continues the run whose state `carry` holds (from [`Xml2Document::take_carry`]): the document
    /// is then the one the earlier calls wrote, so its start is not written again, a report
    /// is numbered after the written ones, and a held report is written when the next arrives.
    pub fn resume(&mut self, carry: Xml2Carry) {
        self.written = carry.written;
        self.started = carry.started;
        self.pending = carry.pending;
        self.pending_made = carry.pending_made;
        self.has_pending = carry.has_pending;
        self.pending_failure = carry.pending_failure;
        self.held = carry.held;
        self.held_refusal = carry.held_refusal;
    }

    /// The posts of the reports written since the last call, at offsets of the output buffers of
    /// the calls that wrote them, and the exception that ends the run after them
    /// (`StructuredXml::take_posts`).
    pub fn take_posts(&mut self) -> DocumentPosts {
        std::mem::take(&mut self.posts)
    }

    pub fn with_files(kind: Xml2Kind, delayed: bool, files: F) -> Self {
        Self {
            kind,
            written: 0,
            started: false,
            delayed,
            pending: Vec::new(),
            pending_made: MadeReport::default(),
            has_pending: false,
            pending_failure: None,
            held: HeldSearches::default(),
            held_refusal: None,
            files,
            posts: DocumentPosts::default(),
        }
    }

    pub fn kind(&self) -> Xml2Kind {
        self.kind
    }

    /// Writes the document's start now, when a program writes it before its first report
    /// (`PrintProlog`, before the search: BLASTX, so that a run that stops at its first batch
    /// still has it; the reports of the other programs before the warnings of their first query,
    /// `StructuredXml::start`). The first report or the epilog then writes no start of its own.
    pub fn start(&mut self, out: &mut Vec<u8>) {
        self.prolog(out);
    }

    /// The start that `start` gave a caller was not written: the caller holds it until the output's
    /// first write or flush, as NCBI's stream holds the start that `PrintProlog` wrote until its
    /// first flush, and the call ended before one (a `-subject` report call that only held its
    /// report: BLASTN's report stream). The next call or the epilog writes it.
    pub fn start_not_written(&mut self) {
        self.started = false;
    }

    /// The document's start, once (`PrintProlog`, blast_format.cpp:360-367).
    fn prolog(&mut self, out: &mut Vec<u8>) {
        if std::mem::replace(&mut self.started, true) {
            return;
        }
        out.extend_from_slice(&prolog_bytes(self.kind));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:877-882
    /// Hands one query's report (`x_PrintXML2Report`): appends to `out` what the framing writes
    /// before the report that this call writes (the document's start if it is not written yet,
    /// outfmt 15's separator `,\n` before every report but the first), then `make` formats the
    /// query's report into `report`, given where the report written now starts in `out`, and the
    /// framing appends that report to `out`; `report` is left empty. NCBI writes the separator
    /// (`x_WriteXML2`, blast_format.cpp:2337-2342) before `s_FillBlastOutput` loads the report's
    /// subjects, and writes the report after those loads: what they post comes between the two.
    /// For a `-subject` run the report is held and the one held before it is written, with what
    /// `make` handed with it (`MadeReport`: its posts come when it is written, after what the
    /// call posted before). A failure of `make` returns it with the framing's bytes in `out` (the
    /// caller truncates `out`) and the document's count of reports unchanged.
    pub fn report(
        &mut self,
        report: &mut Vec<u8>,
        out: &mut Vec<u8>,
        make: impl FnOnce(usize, &mut Vec<u8>) -> Result<MadeReport>,
    ) -> Result<()> {
        // NCBI makes the held report when this result arrives, before it accumulates this one.
        self.release_held()?;
        self.prolog(out);
        // The report written now: this one, or for a `-subject` run the one held before it.
        if !self.delayed || self.has_pending {
            self.separator(out);
        }
        report.clear();
        // A report that raised an exception that NCBI's report function catches
        // (`CaughtException`) is a failed report (`failed`), held as the others are; it made no
        // bytes and posted nothing of its stream.
        let (mut made, failure) = match make(out.len(), report) {
            Ok(made) => (made, None),
            Err(error) => match error.downcast::<crate::report::blast_xml::CaughtException>() {
                Ok(caught) => (MadeReport::default(), Some(caught.0)),
                Err(error) => {
                    report.clear();
                    return Err(error);
                }
            },
        };
        // NCBI makes a report where it writes it: the making that failed now loaded no subject
        // after the one whose HSP raised the exception (`DocumentPosts::cut`); a `-subject` run,
        // whose report is held, reads no database.
        if failure.is_some() && !self.delayed {
            self.posts.cut = true;
        }
        let failure = if self.delayed {
            // The report just made is held; the one held before it is written now.
            let had_pending = std::mem::replace(&mut self.has_pending, true);
            std::mem::swap(report, &mut self.pending);
            std::mem::swap(&mut made, &mut self.pending_made);
            let held = std::mem::replace(&mut self.pending_failure, failure);
            if !had_pending {
                report.clear();
                return Ok(());
            }
            held
        } else {
            failure
        };
        if let Some(exception) = failure {
            report.clear();
            return self.failed(exception);
        }
        let result = self.write(report, made, out);
        report.clear();
        result
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2334-2354
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-778
    /// A report whose making raised `exception` where NCBI writes it: counted as written (the
    /// separator before it stays, so a JSON document can hold `,\n,\n`); the stream forms catch
    /// the exception (`CaughtException`, the program posts its message and goes on), the file
    /// forms (outfmt 13 and 14) have opened the report's file, which stays empty, and end the
    /// program with it (`CATCH_ALL`, exit status 255, 2 for a `CSeqDBException` of a subject's
    /// lookup, `cli::NativeError::exception`; the master file stays empty: no epilog)
    /// after what the making posted (the loads of its subjects): the exception is handed back
    /// after the posts (`DocumentPosts::ends_run`, by `query_warnings::write_hit_posts`).
    fn failed(&mut self, exception: crate::ncbi_exception::NcbiException) -> Result<()> {
        match self.kind {
            Xml2Kind::Single16 | Xml2Kind::Single15 => {
                self.written += 1;
                Err(crate::report::blast_xml::CaughtException(exception).into())
            }
            Xml2Kind::Files14 | Xml2Kind::Files13 => {
                self.file_left_empty()?;
                self.posts.ends_run = Some(Box::new(exception));
                Ok(())
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-778
    /// A report of outfmt 13 or 14 that is not made where NCBI opens its file: it fails (`failed`)
    /// or LOSAT refuses it (`refused`). The report is counted as written and its file is written
    /// empty.
    fn file_left_empty(&mut self) -> Result<()> {
        self.written += 1;
        self.files.write_file(self.written, &[])
    }

    /// Whether a report is held back until the next one arrives (a `-subject` run).
    pub fn holds_back(&self) -> bool {
        self.delayed
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2416
    /// The arrival of a result whose report LOSAT refuses (`refusal`). A document that writes each
    /// report at once (a `-db` run) refuses it now (`refused`): NCBI makes the report here, after
    /// the reports before it. A document that holds each report back (a `-subject` run): NCBI makes
    /// the accumulated report, and crashes on it, only where it writes it (when the result of another
    /// query arrives, in a later batch too, or in the epilog), so the refusal is held in the
    /// report's place and returned there (`release_held`); a run that stops before never makes it.
    /// As at a report's arrival, `out` gets the report held before it (`write_held`, with the
    /// document's start where it is the first), unless the result joins it (`joins_held`): then
    /// the held report is refused (its first refusal kept). The refusal is held whatever the held
    /// report's writing returns (its failure that NCBI's report function catches,
    /// `CaughtException`, is posted and the run goes on). A refused first report thus writes
    /// nothing, not even the start.
    pub fn hold_refused(
        &mut self,
        out: &mut Vec<u8>,
        joins_held: bool,
        refusal: anyhow::Error,
    ) -> Result<()> {
        if !self.delayed {
            return Err(self.refused(refusal));
        }
        if joins_held && self.has_pending {
            self.held_refusal.get_or_insert(refusal);
            return Ok(());
        }
        let written = self.write_held(out);
        self.has_pending = true;
        self.pending.clear();
        self.pending_made = MadeReport::default();
        self.pending_failure = None;
        self.held = HeldSearches::default();
        self.held_refusal = Some(refusal);
        written
    }

    /// Whether the held report is refused (`hold_refused`): a result that joins it adds nothing.
    pub fn held_refused(&self) -> bool {
        self.has_pending && self.held_refusal.is_some()
    }

    /// Where the document writes the held report: its refusal (`hold_refused`, `refused`), if it
    /// has one.
    fn release_held(&mut self) -> Result<()> {
        match self.held_refusal.take() {
            Some(refusal) if self.has_pending => Err(self.refused(refusal)),
            _ => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:768-778
    /// LOSAT's refusal of a report where NCBI makes it (`hold_refused`, `release_held`): the forms
    /// of side files (outfmt 13 and 14) have opened the report's file, which stays empty, as for a
    /// report that fails (`failed`; NCBI BLAST+ 2.17.0 leaves `o_2.xml` empty after `o_1.xml`, and
    /// `o_1.xml` empty for a refused first report, with `-subject` and `-db` alike: remaining-gaps
    /// DECISIONS D-SO-8). A file that cannot be opened fails first (`side_file_open_error`).
    fn refused(&mut self, refusal: anyhow::Error) -> anyhow::Error {
        if self.kind.is_files() {
            if let Err(error) = self.file_left_empty() {
                return error;
            }
        }
        refusal
    }

    /// The held report, written where the result of another query arrives (`hold_refused`) or the
    /// run stops (`unwritten`): its refusal where LOSAT refuses it (`release_held`), else the
    /// document's start where it is not written yet, the separator and the report (or its
    /// failure, `failed`). Nothing where no report is held.
    fn write_held(&mut self, out: &mut Vec<u8>) -> Result<()> {
        if !self.has_pending {
            return Ok(());
        }
        self.release_held()?;
        self.has_pending = false;
        let pending = std::mem::take(&mut self.pending);
        let made = std::mem::take(&mut self.pending_made);
        self.prolog(out);
        self.separator(out);
        if let Some(exception) = self.pending_failure.take() {
            return self.failed(exception);
        }
        self.write(&pending, made, out)
    }

    /// The `Search`es of the report that the next report call holds (`HeldSearches`), which its
    /// writer sets for the report it makes.
    pub(super) fn held_searches_mut(&mut self) -> &mut HeldSearches {
        &mut self.held
    }

    /// The held report of a `-subject` run, for a result that joins it (`HeldReport`); `None`
    /// where no report is held.
    pub(super) fn held_report(&mut self) -> Option<HeldReport<'_>> {
        (self.delayed && self.has_pending).then_some(HeldReport {
            report: &mut self.pending,
            made: &mut self.pending_made,
            failure: &mut self.pending_failure,
            searches: &mut self.held,
        })
    }

    /// The arrival of the report of the last query of a `-subject` run that stops, which NCBI
    /// holds back and never writes: the report is not made, and the one held before it is
    /// written now, as `report` would write it, unless the result joins it (`joins_held`,
    /// `XmlQuery::joins_held`): NCBI accumulates the result, and the held report is never written
    /// either.
    pub fn unwritten(&mut self, out: &mut Vec<u8>, joins_held: bool) -> Result<()> {
        debug_assert!(self.delayed, "only a held report is never written");
        if !self.has_pending || joins_held {
            // Never written, so never made: nor is its refusal. The run stops with NCBI's error,
            // whose exit flushes the start that `PrintProlog` wrote.
            self.has_pending = false;
            self.held_refusal = None;
            self.prolog(out);
            return Ok(());
        }
        // A refused held report crashes NCBI before its stream flushed the start.
        self.write_held(out)
    }

    /// `x_WriteXML2`'s separator of outfmt 15's reports: `,\n` before every report but the first.
    fn separator(&self, out: &mut Vec<u8>) {
        if self.kind == Xml2Kind::Single15 && self.written > 0 {
            out.extend_from_slice(b",\n");
        }
    }

    /// `x_WriteXML2`: one report into the stream (after its separator, `separator`) or its own
    /// file, its posts after the bytes of the output before them (a side file's after the output
    /// so far). A side file whose report threw is written up to the failing byte, and the
    /// exception ends the run after the report's posts (`DocumentPosts::ends_run`, handed back by
    /// the writer of the posts; the master file stays empty: `epilog` writes no end).
    fn write(&mut self, report: &[u8], made: MadeReport, out: &mut Vec<u8>) -> Result<()> {
        self.written += 1;
        let at = if self.kind.is_files() { 0 } else { 1 };
        let base = out.len();
        let MadeReport { posts, ends_run } = made;
        self.posts
            .posts
            .extend(posts.into_iter().map(|post| DocumentPost {
                at: base + at * post.at,
                text: post.text,
            }));
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:775-777
        // The subjects' loads (`s_FillBlastOutput`) and the posts come before the exception.
        if let Some(thrown) = ends_run {
            self.files.write_file(self.written, report)?;
            self.posts.ends_run = Some(Box::new(thrown.exception()));
            return Ok(());
        }
        match self.kind {
            Xml2Kind::Single16 | Xml2Kind::Single15 => {
                // The report is one object of NCBI's stream inside the report function's `try`
                // (`ostream::Site::Xml2Report`), then `Close`'s flush; what follows is outside.
                use crate::report::ostream::{Cut, Site};
                let json = self.kind == Xml2Kind::Single15;
                self.posts
                    .cuts
                    .push(base, Cut::Site(Site::Xml2Report { json }));
                out.extend_from_slice(report);
                self.posts.cuts.push(out.len(), Cut::Flush);
                self.posts.cuts.push(out.len(), Cut::Site(Site::Report));
            }
            Xml2Kind::Files14 | Xml2Kind::Files13 => {
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:770-774
                // The `CArgException` is caught by the application (blast_app_util.hpp:177-180):
                // "Command line argument error: " and the message, exit status 1. The side files
                // written before stay, and the master file stays empty (no epilog).
                self.files.write_file(self.written, report)?;
            }
        }
        Ok(())
    }

    /// `PrintEpilog` (blast_format.cpp:2203-2230): the held report of a `-subject` run, then
    /// the document's end. A run that stops at a failed batch does not call it.
    pub fn epilog(&mut self, out: &mut Vec<u8>) -> Result<()> {
        // A refused held report crashes NCBI before its stream flushed the start.
        if self.has_pending {
            self.release_held()?;
        }
        self.prolog(out);
        // The held report that fails: the end of a stream form follows its message
        // (`CaughtException`, returned after the end); a file form stops here.
        let mut caught = None;
        if self.has_pending {
            self.has_pending = false;
            let pending = std::mem::take(&mut self.pending);
            let made = std::mem::take(&mut self.pending_made);
            self.separator(out);
            match self.pending_failure.take() {
                Some(exception) => match self.failed(exception) {
                    Err(error) if error.is::<crate::report::blast_xml::CaughtException>() => {
                        caught = Some(error)
                    }
                    result => result?,
                },
                None => self.write(&pending, made, out)?,
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:775-777
            // No `catch`: the held report's exception leaves `PrintEpilog` before the document's
            // end (the master file stays empty).
            if self.posts.ends_run.is_some() {
                return Ok(());
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2203-2230
        // The held report's post comes before the document's end (`DocumentPosts::caught_at`).
        if caught.is_some() {
            self.posts.caught_at = Some(out.len());
        }
        match self.kind {
            Xml2Kind::Single16 => out.extend_from_slice(TAIL16),
            Xml2Kind::Single15 => out.extend_from_slice(TAIL15),
            Xml2Kind::Files14 => {
                out.extend_from_slice(MASTER14_HEAD);
                for n in 1..=self.written {
                    out.extend_from_slice(b"\t<xi:include href=\"");
                    out.extend_from_slice(&self.files.listed(n));
                    out.extend_from_slice(b"\"/>\n");
                }
                out.extend_from_slice(b"</BlastXML2>\n");
            }
            Xml2Kind::Files13 => {
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2447-2457
                // The master ends without an end of line.
                out.extend_from_slice(b"{\n\t\"BlastJSON\": [\n");
                for n in 1..=self.written {
                    out.extend_from_slice(b"\t\t{\"File\": \"");
                    out.extend_from_slice(&self.files.listed(n));
                    out.extend_from_slice(b"\" }");
                    if n != self.written {
                        out.push(b',');
                    }
                    out.push(b'\n');
                }
                out.extend_from_slice(b"\t]\n}");
            }
        }
        match caught {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }
}
