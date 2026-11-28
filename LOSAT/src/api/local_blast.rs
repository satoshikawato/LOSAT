//! Local BLAST Search
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp
//!
//! This module provides the main entry point for running local BLAST searches, and the
//! output routing shared by every entry that runs through it (the CLI, web ABI v1 and
//! web ABI v2). A program routed through the shared entry exposes `run_local`, which
//! searches once and hands the result to every requested output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::record::FastaRecord;
use crate::report::blast_xml::XmlStats;
use crate::report::PairwiseHit;

/// What every program's `run_local` returns when its arguments name `-db` (no NCBI counterpart:
/// the caller contract of LOSAT's entries). `run_local` searches the caller's subject records; a
/// `-db` search runs through the program's `run_local_database` (or the command line's `run`).
pub(crate) const RUN_LOCAL_DB_REFUSAL: &str =
    "-db (a database search) is run by run_local_database, not run_local";

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:106-111
pub use crate::algorithm::tblastx::blast_engine::run_local as run_local_tblastx;
// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:113-119
pub use crate::algorithm::tblastx::blast_engine::run_local_database as run_local_database_tblastx;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:128-133
pub use crate::algorithm::blastn::blast_engine::run_local as run_local_blastn;
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
pub use crate::algorithm::blastn::blast_engine::run_local_database as run_local_database_blastn;

// Re-export the pool of a run of several searches (`run_local` calls of one argv), which they
// share instead of building one each (no NCBI counterpart: NCBI's CLI searches once per run).
pub use crate::utils::threading::with_run_pool;

// Re-export BLASTP run function
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:167-299
pub use crate::algorithm::blastp::blast_engine::run_local as run_local_blastp;
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:194-199
pub use crate::algorithm::blastp::blast_engine::run_local_database as run_local_database_blastp;

// Re-export the TBLASTN entry
// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:289-301
pub use crate::algorithm::tblastn::run_local as run_local_tblastn;
// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
pub use crate::algorithm::tblastn::run_local_database as run_local_database_tblastn;

// Re-export BLASTX's shared entry (web ABI v2).
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
pub use crate::algorithm::blastx::run_local as run_local_blastx;
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-200
pub use crate::algorithm::blastx::run_local_database as run_local_database_blastx;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
/// Where one requested output format is written.
#[non_exhaustive]
pub enum OutputSink<'a> {
    /// Standard output, buffered (the CLI without `-out`, NCBI's default `-`).
    Stdout,
    /// A file named by `-out`, created when the report is written.
    File(&'a Path),
    /// A writer provided by the caller. It is not buffered here, so that every byte a
    /// formatter has written reaches it before each `FormatObserver` call.
    Writer(&'a mut (dyn Write + Send)),
    /// The `-out` file of the library's `run_to_file`, made where the command line makes its
    /// `-out` (`ReportOutputs::open_out`); otherwise as `Writer`.
    Out(&'a mut OutFile),
    /// A writer provided by the caller, as `Writer`, and the caller's receiver of the side files
    /// of outfmt 13 and 14 (a caller without files: LOSAT Web). The writer gets the master file,
    /// the receiver each side file, named as next to the receiver's `-out`
    /// (`SideFileSink::out`), which the format takes as its `-out`: its XInclude check passes.
    /// For any other format it is `Writer`.
    Files(
        &'a mut (dyn Write + Send),
        crate::report::blast_xml2::SideFileSink,
    ),
    /// The command line's report stream (`-out` or standard output), which writes as NCBI's
    /// stream does (`report::ostream::ReportStream`); otherwise as `Writer`.
    Report(&'a mut crate::report::ostream::ReportStream),
}

impl OutputSink<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
    // NCBI opens the `-out` stream while it processes the arguments, before the search.
    // LOSAT has always created the file when the report is written, after the search;
    // this routing keeps that existing timing unchanged.
    /// Opens the sink for writing the report.
    pub fn open(&mut self) -> io::Result<Box<dyn Write + '_>> {
        Ok(match self {
            OutputSink::Stdout => Box::new(BufWriter::new(io::stdout())),
            OutputSink::File(path) => Box::new(BufWriter::new(File::create(path)?)),
            OutputSink::Writer(writer) | OutputSink::Files(writer, _) => Box::new(&mut **writer),
            OutputSink::Out(out) => Box::new(&mut **out),
            OutputSink::Report(stream) => Box::new(&mut **stream),
        })
    }

    /// The command line's report stream behind the sink (`OutputSink::Report`), or the `-out`
    /// stream of `run_to_file` (`OutputSink::Out`, `OutFile::report_handle`), which takes NCBI's
    /// regions and flushes (`report::ostream::ReportWrite`).
    pub(crate) fn report_handle(&self) -> Option<crate::report::ostream::ReportHandle> {
        match self {
            OutputSink::Report(stream) => Some(stream.handle()),
            OutputSink::Out(out) => out.report_handle(),
            _ => None,
        }
    }

    /// The caller's receiver of the side files of outfmt 13 and 14 (`OutputSink::Files`).
    pub fn side_files(&self) -> Option<&crate::report::blast_xml2::SideFileSink> {
        match self {
            OutputSink::Files(_, sink) => Some(sink),
            _ => None,
        }
    }

    /// Whether the output is a writer of the caller (`Writer` or `Files`), or the `-out` file of
    /// `run_to_file`, which is written as one (`Out`).
    pub(crate) fn is_writer(&self) -> bool {
        matches!(
            self,
            OutputSink::Writer(_)
                | OutputSink::Files(..)
                | OutputSink::Out(_)
                | OutputSink::Report(_)
        )
    }

    /// The same sink, borrowed for a view of the outputs (`ViewParts::view`).
    fn reborrow(&mut self) -> OutputSink<'_> {
        match self {
            OutputSink::Stdout => OutputSink::Stdout,
            OutputSink::File(path) => OutputSink::File(path),
            OutputSink::Writer(writer) => OutputSink::Writer(&mut **writer),
            OutputSink::Out(out) => OutputSink::Out(out),
            OutputSink::Report(stream) => OutputSink::Report(stream),
            OutputSink::Files(writer, side) => OutputSink::Files(&mut **writer, side.clone()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:913-916 (cited at
// `ReportOutputs::finish_archive`): `PrintErrorArchive` opens `-out` for the error archive of an
// error that came before the argument handling opened it.
/// The `-out` file of the library's `run_to_file` (`OutputSink::Out`). It is made where the
/// command line makes `-out`: each program's entry calls `ReportOutputs::open_out` after it has
/// read the subjects (or the taxid or seqid list of `-db`) and before it processes the options,
/// so an error after that point leaves the file empty and an error before it leaves no file, as on
/// the command line. The first write makes it too: the error archive of an earlier error (outfmt
/// 11, `ReportOutputs::finish_archive`), as NCBI's `PrintErrorArchive` opens `-out`. A making that
/// failed is not tried again. The file is written through the command line's `-out` stream of the
/// format (`report::ostream::ReportStream::of_library_file`: NCBI's buffer, regions and flushes),
/// whose posts are the call's diagnostics in the command line's order (`Posts`), and the call's
/// run ends it as the command line's run ends its `-out` (`OutFile::finish`).
pub struct OutFile {
    path: std::path::PathBuf,
    /// The application's name in NCBI's posts.
    app: &'static str,
    /// The file's `-outfmt` value.
    outfmt: String,
    /// The posts of the file's stream, until the call's diagnostics take them.
    posts: Posts,
    state: OutState,
}

/// Whether `OutFile` has been made.
enum OutState {
    Unmade,
    Made {
        stream: crate::report::ostream::ReportStream,
        /// Whether `OutFile::finish` ended it.
        ended: bool,
    },
    Failed,
}

impl OutFile {
    /// The file `path` of a search (`cli::create_out_file`), the application `app`, that writes
    /// the report of `outfmt`, not made yet.
    pub(crate) fn new(path: &Path, app: &'static str, outfmt: &str) -> Self {
        OutFile {
            path: path.to_path_buf(),
            app,
            outfmt: outfmt.to_owned(),
            posts: Posts::default(),
            state: OutState::Unmade,
        }
    }

    /// The posts that the file's stream makes (`PostsFirst` writes them to the call's
    /// diagnostics).
    pub(crate) fn posts(&self) -> Posts {
        self.posts.clone()
    }

    /// Makes the file as the command line makes `-out` (`cli::create_out_file`), once, behind the
    /// command line's `-out` stream of the format (`ReportStream::of_library_file`).
    fn make(&mut self) -> anyhow::Result<()> {
        match self.state {
            OutState::Unmade => match crate::cli::create_out_file(&self.path) {
                Ok(file) => {
                    // No strategy file comes before the report: the library refuses
                    // `-export_search_strategy` (`search::file_options`).
                    self.state = OutState::Made {
                        stream: crate::report::ostream::ReportStream::of_library_file(
                            file,
                            self.app,
                            &self.outfmt,
                            Box::new(self.posts.clone()),
                        ),
                        ended: false,
                    };
                    Ok(())
                }
                Err(error) => {
                    self.state = OutState::Failed;
                    Err(error)
                }
            },
            OutState::Made { .. } => Ok(()),
            OutState::Failed => Err(anyhow::anyhow!("the -out file could not be made")),
        }
    }

    /// The file's stream, which takes NCBI's regions and flushes, once the file is made
    /// (`OutputSink::report_handle`).
    pub(crate) fn report_handle(&self) -> Option<crate::report::ostream::ReportHandle> {
        match &self.state {
            OutState::Made { stream, .. } => Some(stream.handle()),
            OutState::Unmade | OutState::Failed => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
    /// The end of the call's run, `result` being its result, as the command line's run ends its
    /// `-out` stream (`ReportStream::finish`). A file not made has nothing to end.
    /// The file ends once: the error archive after the run's end is written on the ended stream
    /// (`ReportOutputs::finish_archive`), as the command line's.
    pub(crate) fn finish(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        match &mut self.state {
            OutState::Made { stream, ended } if !*ended => {
                *ended = true;
                stream.finish(result)
            }
            OutState::Made { .. } | OutState::Unmade | OutState::Failed => result,
        }
    }
}

impl Write for OutFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if matches!(self.state, OutState::Unmade) {
            self.make().map_err(io::Error::other)?;
        }
        match &mut self.state {
            OutState::Made { stream, .. } => stream.write(buf),
            _ => Err(io::Error::other("the -out file could not be made")),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.state {
            OutState::Made { stream, .. } => stream.flush(),
            _ => Ok(()),
        }
    }
}

/// The posts of NCBI's that the `-out` stream of `run_to_file` makes while its report is written
/// (`OutFile`), which the command line writes to its standard error, until the call's diagnostics
/// take them (`PostsFirst`).
#[derive(Clone, Default)]
pub(crate) struct Posts(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

impl Posts {
    fn take(&self) -> Vec<u8> {
        std::mem::take(
            &mut *self
                .0
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner),
        )
    }
}

impl Write for Posts {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The call's diagnostics `inner`, which take the `-out` stream's posts so far (`Posts`) before
/// each of their own writes and at `drain`: the posts keep their place among the call's messages,
/// as on the command line's standard error.
pub(crate) struct PostsFirst<'a> {
    posts: Posts,
    inner: &'a mut (dyn Write + Send),
}

impl<'a> PostsFirst<'a> {
    pub(crate) fn new(posts: Posts, inner: &'a mut (dyn Write + Send)) -> Self {
        PostsFirst { posts, inner }
    }

    /// Writes the posts made so far.
    pub(crate) fn drain(&mut self) -> io::Result<()> {
        let posts = self.posts.take();
        if posts.is_empty() {
            return Ok(());
        }
        self.inner.write_all(&posts)
    }
}

impl Write for PostsFirst<'_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.drain()?;
        self.inner.write(buf)
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        self.drain()?;
        self.inner.write_all(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.drain()?;
        self.inner.flush()
    }
}

/// The report of the `-outfmt` value `outfmt` (`None` for a value that NCBI's reading refuses).
fn report_format_of(outfmt: impl AsRef<[u8]>) -> Option<crate::blastinput::app::ReportFormat> {
    use crate::blastinput::app::{parse_formatting_string, ReportFormat};
    parse_formatting_string(outfmt)
        .ok()
        .and_then(|choice| ReportFormat::from_number(choice.number))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
/// One requested output: an `-outfmt` value, parsed and validated by the program exactly
/// as on the command line, where to write it, and what LOSAT's refusal of its report does.
pub struct FormatOutput<'a> {
    pub outfmt: &'a [u8],
    pub sink: OutputSink<'a>,
    pub refusal: FormatRefusal,
}

/// What an error of one format's report does (`FormatRefusal::check`), the error that the command
/// line's run of that format alone stops with: LOSAT's refusals that a program's check of the
/// format makes before the report's rows (a tabular row at which NCBI stops,
/// `report::tabular::check_failed_lookups`, every program; the `-db` list's counts that the
/// prolog of outfmt 0 to 4 prints where NCBI's sums overflow,
/// `common::database::check_formats_report_sums`, every program; BLASTP's subject title of
/// punctuation in outfmt 0 to 4; BLASTX's subject of empty rows that the mapper rewrites,
/// `blastx::empty_rows::refuse_mapped_without_scores`), and every error that ends the format's
/// report in its middle (LOSAT's refusals there, such as the `-sorthits` total that NCBI pads
/// without end, `report::pairwise::TotalScoreRefusal`, and NCBI's errors of the format alone,
/// such as BLASTX's `Row is empty` in outfmt 0): each program passes the result of each format's
/// report, of each query batch where it writes them batch by batch, to its format's `check`.
#[derive(Debug, Default)]
pub enum FormatRefusal {
    /// The error fails the call: the command line's one format, the library's reports.
    #[default]
    FailsCall,
    /// The error fails this format alone: the call writes no more of its report and writes the
    /// other formats (LOSAT Web's formats that the argv does not ask for, docs/web/abi_v2.md §7).
    /// What the call wrote of it before (outfmt 0's prolog to a writer before the search, the
    /// rows before the error, those of BLASTX's earlier query batches) is the caller's to drop.
    FailsAlone,
    /// The error of a format that failed alone.
    Refused(anyhow::Error),
}

impl FormatRefusal {
    /// `checked`, a check of the format's report before its rows or the result of (a query batch
    /// of) its report (`Err`: the error that stops the format): `Ok(true)` where the report is
    /// written, `Ok(false)` where it is not (the format fails alone, now or earlier: `Refused`
    /// keeps the first error), and the error itself where it fails the call.
    pub(crate) fn check(&mut self, checked: anyhow::Result<()>) -> anyhow::Result<bool> {
        match (checked, &*self) {
            (_, FormatRefusal::Refused(_)) => Ok(false),
            (Ok(()), _) => Ok(true),
            (Err(error), FormatRefusal::FailsCall) => Err(error),
            (Err(error), FormatRefusal::FailsAlone) => {
                *self = FormatRefusal::Refused(error);
                Ok(false)
            }
        }
    }

    /// Whether the format failed alone (`Refused`): the call writes no more of it.
    pub(crate) fn refused(&self) -> bool {
        matches!(self, FormatRefusal::Refused(_))
    }

    /// The refusal of a format that failed alone, taken by the caller after the call.
    pub fn take_refused(&mut self) -> Option<anyhow::Error> {
        match std::mem::take(self) {
            FormatRefusal::Refused(error) => Some(error),
            other => {
                *self = other;
                None
            }
        }
    }
}

impl FormatOutput<'_> {
    /// The `-out` of this format: the one that the caller's receiver of side files names them by
    /// (`OutputSink::Files`), else the run's `out`. Outfmt 13 and 14 write their side files next
    /// to it and fail without one (`blastinput::app::xinclude_check_for`).
    pub(crate) fn side_out<'s>(&'s self, out: Option<&'s Path>) -> Option<&'s Path> {
        self.sink
            .side_files()
            .map(crate::report::blast_xml2::SideFileSink::out)
            .or(out)
    }

    /// The report of this output's `-outfmt` (`None` for a value that NCBI's reading refuses).
    pub(crate) fn report_format(&self) -> Option<crate::blastinput::app::ReportFormat> {
        report_format_of(self.outfmt)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1452
    /// Whether this output's report writes the queries' result warnings
    /// (`ReportFormat::writes_query_warnings`; a value that NCBI's reading refuses writes no report,
    /// and counts as one that does).
    pub(crate) fn writes_query_warnings(&self) -> bool {
        self.report_format()
            .is_none_or(crate::blastinput::app::ReportFormat::writes_query_warnings)
    }

    /// This format's receiver, among the caller's (`later`, `ReportOutputs::later_warnings`; the
    /// format is the call's `index`th), of the queries' result warnings that its reports write:
    /// `None` where they are not kept for the later formats (`keeps`,
    /// `ReportOutputs::keeps_later_query_warnings`) or this format's report writes none. The
    /// format's writer writes them there before each query's report
    /// (`QueryWarnings::with_later`).
    ///
    /// NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
    pub(crate) fn later_warnings<'s>(
        &self,
        keeps: bool,
        later: &'s mut Option<&mut [Vec<u8>]>,
        index: usize,
    ) -> Option<&'s mut Vec<u8>> {
        match keeps && self.writes_query_warnings() {
            true => later.as_deref_mut()?.get_mut(index),
            false => None,
        }
    }

    /// Whether this output is the BLAST archive (outfmt 11).
    fn is_archive(&self) -> bool {
        self.report_format() == Some(crate::blastinput::app::ReportFormat::ArchiveFormat)
    }
}

/// `ReportOutputs::open_out` of the outputs `formats` (for a caller that holds the diagnostics
/// stream of the same outputs).
pub(crate) fn open_out(formats: &mut [FormatOutput<'_>]) -> anyhow::Result<()> {
    for format in formats {
        if let OutputSink::Out(out) = &mut format.sink {
            out.make()?;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
// Every output format prints from the same final result set, so a position in that
// result identifies an HSP in every format.
/// Index of an HSP in the final, sorted hit list of a search: the list that
/// `ReportOutputs::hits` receives. It identifies the HSP in every output format.
pub type HspIndex = usize;

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1970-1973
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3956-3975
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3613-3632
// The observed units follow NCBI's: one printed row per HSP in outfmt 6 and 7, and in
// outfmt 0 the output of one x_DisplayAlnvecInfo call without the subject defline that
// it prints before the first HSP of each subject (the section starts at " Score =").

// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
/// A sequence of a database that a search reported (`FormatObserver::database_subjects`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportedSubject {
    /// The sequence's OID in the database (the `n` of `gnl|BL_ORD_ID|n`), the same in every
    /// search of the database, also one that a taxid or Seq-id list restricts.
    pub oid: u32,
    /// The `saccver` of the default tabular report, its bytes: where NCBI's places give another
    /// database sequence's Bioseq for the sequence's Seq-align id (its holder), the holder's, as
    /// the tables write it.
    pub id: crate::api::Id,
    /// The letters of the sequence whose `saccver` `id` is (the holder's length, the tables'
    /// `slen`, where `id` is the holder's).
    pub length: usize,
}

/// Reports which HSP a formatter is writing, so that a caller can relate each HSP to its
/// bytes in every output, and the statistics of each query of the search. The written bytes do
/// not change. The query-anchored reports (outfmt 1 to 4) write one section of all of a query's
/// HSPs merged, which has no part of one HSP or subject: they call none of the HSP and subject
/// events.
pub trait FormatObserver {
    /// The formatter for `ReportOutputs::formats[format]` starts the row (outfmt 6/7) or
    /// the alignment section (outfmt 0: the score lines and the alignment, without the
    /// subject defline) of `hsp`. Every byte written for that format before this call
    /// has reached its sink.
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex);
    /// The formatter for `ReportOutputs::formats[format]` finished the row or section of
    /// `hsp`. Every byte of it has reached the sink.
    fn hsp_end(&mut self, format: usize, hsp: HspIndex);

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3613-3632
    /// The outfmt 0 formatter for `ReportOutputs::formats[format]` starts the heading of
    /// a subject (its defline and length line), which precedes the first HSP of that
    /// subject, `first_hsp`, and the HSPs of the copies of one database record that follow it
    /// and that the heading also covers (`report::db_subject::joins_as_copy`). Every byte
    /// written for that format before this call has reached its sink. Formatters without
    /// subject headings never call it.
    fn subject_begin(&mut self, _format: usize, _first_hsp: HspIndex) {}
    /// The formatter finished the heading that `subject_begin` announced. Every byte of
    /// it has reached the sink.
    fn subject_end(&mut self, _format: usize, _first_hsp: HspIndex) {}

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-94
    // NCBI keeps one `CBlastAncillaryData` per query with the query's result; here each program
    // hands over what its XML reports print of it, whichever formats were asked for.
    /// The statistics of query `q_idx` of the search: the values that the `Statistics` of its
    /// `Iteration` in outfmt 5 shows (`XmlStats::NONE` for a query without a valid context).
    /// `q_idx` counts as `Hit::q_idx` of the final hit list. The call comes once per query,
    /// before the reports of the query are written, and only with an observer (the program does
    /// not make these values otherwise).
    fn query_stats(&mut self, _q_idx: usize, _stats: XmlStats) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // What NCBI writes before this loop is written once for the run; a call of `run_local`
    // per query batch writes it again, and this event tells where it ends.
    /// The checks of the run are done and the subjects are read: what the diagnostics receive
    /// from here on belongs to the queries.
    fn queries_begin(&mut self) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:299-300
    /// BLASTN's query batch just searched made `extensions` successful initial extensions, which
    /// size the next batch (`CBatchSizeMixer`). The other programs' batches have one size.
    fn batch_extensions(&mut self, _extensions: i32) {}

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
    // NCBI's reports read a subject's ids and length by the Seq-align's id; a database search's
    // hit list names its sequences by rank, and this event tells what the ranks are.
    /// The sequences that `Hit::s_idx` of a database search's final hit list (`ReportOutputs::hits`)
    /// names, by that index. The call comes once, with `ReportOutputs::hits` only (the program
    /// does not make these values otherwise).
    fn database_subjects(&mut self, _subjects: Vec<ReportedSubject>) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537-543
    // A PSSM search's query is the PSSM's Bioseq, which the inputs do not hold
    // (`api::Inputs::pssm`): the search reads it with the options.
    /// The query of a PSSM search (psiblast's and tblastn's `-in_pssm`): the record of the PSSM's
    /// sequence as the reports show it, the query 0 of the hit lists. The call comes once, before
    /// the search's queries (`queries_begin`), and only for a PSSM query.
    fn pssm_query(&mut self, _query: &FastaRecord) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:250-285
    // NCBI's psiblast searches one query per round and formats each round's result set; LOSAT's
    // round driver makes one search call per round of a query (`blastp::blast_engine::psi_round`).
    /// psiblast's round driver starts round `round` (1-based, NCBI's `GetIterationNumber`) of the
    /// searched query `query` (the index that `Hit::q_idx` counts in a search without rounds):
    /// the next hit list (`ReportOutputs::hits`), its HSP events, `query_stats` (whose `q_idx` is
    /// 0, the round's only query) and `database_subjects` belong to that round, and its HSP
    /// numbers start again at 0. A search without rounds never calls it.
    fn round_begin(&mut self, _query: usize, _round: u32) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:525-528
    /// The rounds of the searched query `query` converged (NCBI's `HasConverged`): every subject
    /// below `-inclusion_ethresh` in its last round was below it in the round before too. It
    /// comes after the query's last round, whatever the formats.
    fn rounds_converged(&mut self, _query: usize) {}
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1970-1973
/// The observer of one output format, as handed to that format's formatter.
pub struct FormatProbe<'o> {
    observer: &'o mut dyn FormatObserver,
    format: usize,
}

impl<'o> FormatProbe<'o> {
    pub fn new(observer: &'o mut dyn FormatObserver, format: usize) -> Self {
        Self { observer, format }
    }

    pub fn begin(&mut self, hsp: HspIndex) {
        self.observer.hsp_begin(self.format, hsp);
    }

    pub fn end(&mut self, hsp: HspIndex) {
        self.observer.hsp_end(self.format, hsp);
    }

    pub fn subject_begin(&mut self, first_hsp: HspIndex) {
        self.observer.subject_begin(self.format, first_hsp);
    }

    pub fn subject_end(&mut self, first_hsp: HspIndex) {
        self.observer.subject_end(self.format, first_hsp);
    }
}

/// The receiver of the final hit list (`ReportOutputs::hits`).
pub type HitsSink<'a> = dyn FnMut(std::sync::Arc<Vec<PairwiseHit>>) + Send + 'a;

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:195-295
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
// NCBI formats one CSearchResultSet without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// The outputs of one search. Every part is `Send`, because the search runs the
/// formatters inside its thread pool.
pub struct ReportOutputs<'a> {
    /// Every requested output format. The CLI requests exactly one.
    pub formats: Vec<FormatOutput<'a>>,
    /// Receives the warnings that the CLI writes to standard error, each once.
    pub diagnostics: &'a mut (dyn Write + Send),
    /// By format, as `formats`: receives what the command line's run of a format after the first
    /// writes to its standard error with its queries' reports and the first format's run does
    /// not. A structured first format (outfmt 5, 8, 9, 11 to 17) writes none of the queries'
    /// result warnings (`ReportFormat::writes_query_warnings`), and the run of a later format
    /// that writes them posts each when the query's report starts: such a format's writer writes
    /// them to `later_warnings[i]` there, `i` being its place in `formats`, so a format whose
    /// report stops (it failed alone, `FormatRefusal`) has given those of the queries that it
    /// reached, and one refused before its report none. BLASTX gives each such format's whole
    /// standard error of the reports, which may repeat lines of the diagnostics. The texts of
    /// several formats repeat each other: the caller keeps each line once (LOSAT Web writes
    /// outfmt 0, 6 and 7 besides the argv's `-outfmt`, and adds what their standard error adds;
    /// line WQA's sweep, D3 review r3, 2026-11-19). Nothing is written where the first format
    /// writes the warnings itself, nor to a place beyond the slice. `None` drops them, as the
    /// command line, whose one format is the first. In BLASTN, BLASTP, TBLASTN and TBLASTX, a
    /// later format that gets a receiver also writes its writer's own lines (the taxonomy
    /// lookups of a tabular field list, outfmt 20's prolog warning), and they go to the
    /// diagnostics, not to its receiver, so receivers can add lines to the diagnostics there;
    /// BLASTX gives them to the receiver. LOSAT Web's later formats are the plain outfmt 0, 6
    /// and 7, which write none (reviews r2 L1 and r3 L2 of 2026-11-20).
    pub later_warnings: Option<&'a mut [Vec<u8>]>,
    /// Receives the final, sorted hit list with rendered alignments, once per call. It may come
    /// before or after the reports are written (blastp's query blocks hand it over after the
    /// last block, `LOSAT_X_BPQBLOCK`). The list is shared, not copied: it holds an alignment's
    /// rows for every HSP, which the reports read while the caller keeps it for its records.
    pub hits: Option<&'a mut HitsSink<'a>>,
    /// Receives the row and section boundaries of each HSP in each format. The events of one
    /// format come in that format's write order; those of different formats may interleave in
    /// any order (blastp's query blocks write every format block by block).
    pub observer: Option<&'a mut (dyn FormatObserver + Send)>,
    /// The `-outfmt` values of the call's other searches (`api::run`'s searches, one for each
    /// hit list size, `search_groups`), which this search does not write: the refusal of a
    /// subject whose lookup throws where its batch's prefetch is off is the call's, decided over
    /// all of its formats (`report::db_subject::FormatThrow::of_search`; review 3 of line DBID
    /// pfoff-3, N2). Empty where the call is one search (the command line, every other caller).
    pub other_outfmts: &'a [&'a str],
    /// Receives the posts of the `-html` reports' reads of NCBI's formatter registries
    /// (`report::html_ids::HtmlPost`: `./.ncbirc` and `InitConfig`'s file) where `diagnostics`
    /// takes the repeats of an earlier search of the same call (`api::run`'s searches after the
    /// first, `search_groups`): such a report posts what a run of its format alone posts, and
    /// only its search writes it (line RG T4a2 fix round 2, review c C3 = review d D5). `None`
    /// writes them to `diagnostics` (the command line, every other caller).
    pub html_posts: Option<&'a mut (dyn Write + Send)>,
}

impl<'a> ReportOutputs<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1452 (cited at
    // `FormatOutput::writes_query_warnings`): a structured report returns before the warnings.
    /// Whether the first format writes the queries' result warnings with its reports (a call
    /// without formats too): a structured format does not (`ReportFormat::writes_query_warnings`).
    pub(crate) fn first_writes_query_warnings(&self) -> bool {
        self.formats
            .first()
            .is_none_or(FormatOutput::writes_query_warnings)
    }

    /// Whether the queries' result warnings are kept for the formats after the first
    /// (`Posts::later`), each of which writes them to its own receiver
    /// (`FormatOutput::later_warnings`): the first format writes none and the caller takes the
    /// later formats' (`later_warnings`).
    ///
    /// NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
    pub(crate) fn keeps_later_query_warnings(&self) -> bool {
        self.later_warnings.is_some() && !self.first_writes_query_warnings()
    }

    /// Makes the `-out` file of `run_to_file` (`OutputSink::Out`) where the command line opens its
    /// `-out`: each program's entry calls it after reading the subjects (or the taxid or seqid list
    /// of `-db`) and before processing the options. The other sinks are not opened here.
    pub(crate) fn open_out(&mut self) -> anyhow::Result<()> {
        open_out(&mut self.formats)
    }

    /// The end of a run whose result is `result`: the `-out` file of `run_to_file`
    /// (`OutputSink::Out`) ends as the command line's `-out` stream ends (`OutFile::finish`); the
    /// other sinks are the caller's.
    pub(crate) fn finish_out(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        self.formats
            .iter_mut()
            .fold(result, |result, format| match &mut format.sink {
                OutputSink::Out(out) => out.finish(result),
                _ => result,
            })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627 (cited at
    // `app::early_format_choice`)
    /// Each requested format as the command line's early reading of its `-outfmt` reads it before
    /// the option handlers (`app::early_format_choice`), with LOSAT's refusal of a format that
    /// `program` does not write (`EarlyFormat::report_format`, with the format's `-out`,
    /// `FormatOutput::side_out` of the run's `out`). `strategy`: the argv names a search strategy
    /// file (`-export_search_strategy`), with which a value that does not parse is left to the
    /// program's formatting handler, after the arguments' own errors, as on the command line (line
    /// RG F1w review 2, S-a: LOSAT Web gives such a value as the run's one format and as the
    /// arguments' `-outfmt`, which the handlers read again).
    pub(crate) fn early_choices(
        &self,
        program: crate::blastinput::app::Program,
        out: Option<&Path>,
        strategy: bool,
    ) -> anyhow::Result<Vec<crate::blastinput::app::EarlyFormat>> {
        self.formats
            .iter()
            .map(|format| {
                let choice = crate::blastinput::app::early_format_choice(format.outfmt, strategy)?;
                choice.report_format(program, format.side_out(out))?;
                Ok(choice)
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:887-899
    /// The first requested format that is an XInclude format (13, 14) whose side files have no
    /// file name to go by (`FormatOutput::side_out` of the run's `-out` `out`: no
    /// `api::run_to_file` file and no receiver of them), for which NCBI's `UseXInclude` fails
    /// (`blastinput::app::xinclude_check_for`, the one check of the command line and the library).
    pub(crate) fn xinclude_choice(
        &self,
        out: Option<&std::path::Path>,
    ) -> anyhow::Result<Option<crate::blastinput::app::FormatChoice>> {
        for format in &self.formats {
            let choice = crate::blastinput::app::parse_formatting_string(format.outfmt)?;
            if crate::blastinput::app::xinclude_check_for(&choice, format.side_out(out)).is_err() {
                return Ok(Some(choice));
            }
        }
        Ok(None)
    }

    /// NCBI's error of the first XInclude format without a file name (`xinclude_choice`). Each
    /// program's entry calls it where its command line makes the formatter: after `Query is
    /// Empty!` and the subjects' warnings.
    pub(crate) fn xinclude_check(&self, out: Option<&std::path::Path>) -> anyhow::Result<()> {
        match self.xinclude_choice(out)? {
            Some(choice) => crate::blastinput::app::xinclude_check_for(&choice, None),
            None => Ok(()),
        }
    }

    /// Whether every requested format is the BLAST archive (outfmt 11), which NCBI's apps but
    /// psiblast write without fetching the report's sequences ahead
    /// (`blastinput::app::check_search_environment`); a call without formats (hit records
    /// only) is not.
    pub(crate) fn only_archive(&self) -> bool {
        !self.formats.is_empty() && self.formats.iter().all(FormatOutput::is_archive)
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:913-916
    /// `handle` with the `-outfmt` text of the requested archive format: an error archive
    /// follows that text (`ArchiveHandle::with_outfmt_text`).
    pub(crate) fn with_archive_outfmt(
        &self,
        handle: crate::report::archive::ArchiveHandle,
    ) -> crate::report::archive::ArchiveHandle {
        match self.formats.iter().find(|format| format.is_archive()) {
            Some(format) => handle.with_outfmt_text(format.outfmt),
            None => handle,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
    /// The end of a run with a BLAST archive (outfmt 11) that `result` closes: NCBI writes an
    /// error archive after the archives of the batches before when the application stops with one
    /// of its own errors, and when the run ends with messages that no archive carried. `no_query`:
    /// the run ended at `Query is Empty!`, where NCBI returns before its `PrintErrorArchive`.
    ///
    /// The `-out` file of `run_to_file` ends first, as the command line's report stream ends
    /// (`~CBlastFormat`, `CATCH_ALL`) before `PrintErrorArchive`, and the archive is written on
    /// its stream (`OutFile`): what that stream posted before is the call's diagnostics, and its
    /// posts while it writes the archive come after the error's text on the command line's
    /// standard error (`archive::write_error_archive_after_post`), so they end the error's
    /// message.
    pub(crate) fn finish_archive(
        &mut self,
        handle: &crate::report::archive::ArchiveHandle,
        no_query: bool,
        result: anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let mut result = self.finish_out(result);
        let _ = self.diagnostics.flush();
        let lines = match &result {
            // The errors of `Run`, the files' among them (`run_local_database` opens no file
            // that the argument handling opens).
            Err(error) => crate::report::archive::error_lines(error, true)
                .or_else(|| crate::report::archive::silent_end_lines(handle, error)),
            Ok(()) if no_query => None,
            Ok(()) => handle.has_messages().then(Vec::new),
        };
        if let (Some(lines), Some(format)) = (
            lines,
            self.formats.iter_mut().find(|format| format.is_archive()),
        ) {
            // `PrintErrorArchive` opens `-out` before the archive's encoding is read
            // (`ARCHIVE_FORMAT`, which can write nothing), and ignores a file it cannot open: the
            // file of `run_to_file` is made here when the search ended before making it.
            if let (OutputSink::Out(out), true) = (&mut format.sink, handle.error_archive_enabled())
            {
                let _ = out.make();
            }
            let written = {
                let mut writer = crate::report::ostream::NcbiOstream::sink(&mut format.sink)?;
                crate::report::archive::write_error_archive(handle, &lines, &mut writer)
                    .and_then(|()| writer.flush().map_err(anyhow::Error::from))
            };
            if let (Err(error), OutputSink::Out(out)) = (&mut result, &format.sink) {
                let after = out.posts().take();
                if let Some(native) = error.downcast_mut::<crate::cli::NativeError>() {
                    native.message.extend_from_slice(&after);
                }
            }
            if result.is_ok() {
                written?;
            }
        }
        result
    }

    /// Runs `f` with the same outputs and a diagnostics stream that also saves what it
    /// writes as the messages of the BLAST archive (outfmt 11) of `handle`.
    pub(crate) fn with_archive_diagnostics<R>(
        &mut self,
        handle: crate::report::archive::ArchiveHandle,
        f: impl FnOnce(&mut ReportOutputs<'_>) -> R,
    ) -> R {
        let (formats, diagnostics, mut parts) = self.split();
        let mut tee = crate::report::archive::DiagTee::new(diagnostics, handle);
        let mut view = parts.view(
            formats
                .iter_mut()
                .map(|format| FormatOutput {
                    outfmt: format.outfmt,
                    sink: format.sink.reborrow(),
                    // Moved into the view and back after `f` (a format that fails alone).
                    refusal: std::mem::take(&mut format.refusal),
                })
                .collect(),
            &mut tee,
        );
        let result = f(&mut view);
        let refusals = view.into_refusals();
        FormatOutput::take_back(formats, refusals);
        result
    }

    /// The outputs as their formats, their diagnostics stream and the rest (`ViewParts`), for a
    /// view of them that sees the formats through other sinks or the diagnostics through another
    /// stream (`ViewParts::view`).
    fn split(
        &mut self,
    ) -> (
        &mut Vec<FormatOutput<'a>>,
        &mut (dyn Write + Send),
        ViewParts<'_, 'a>,
    ) {
        let ReportOutputs {
            other_outfmts,
            formats,
            diagnostics,
            later_warnings,
            hits,
            observer,
            html_posts,
        } = self;
        (
            formats,
            &mut **diagnostics,
            ViewParts {
                other_outfmts,
                later_warnings,
                hits,
                observer,
                html_posts,
            },
        )
    }

    /// The refusal of each format of a view (`ViewParts::view`), which `FormatOutput::take_back`
    /// gives back to the formats it was lent from.
    fn into_refusals(self) -> Vec<FormatRefusal> {
        self.formats
            .into_iter()
            .map(|format| format.refusal)
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:250-262 (cited at
    // `strategy::check_report_stream`): the report stream before NCBI's first request.
    /// A run of records (each program's `run_local` and `run_local_database`) as the command
    /// line's run of the argv ends: `run` searches and reports, with the BLAST archive (outfmt 11)
    /// of `archive`, whose diagnostics stream saves the messages (`with_archive_diagnostics`); the
    /// run ends with the search strategy file of the arguments (`strategy`,
    /// `strategy::end_run`, LOSAT Web's: `strategy::open_strategy_file_in`), and then the archive's
    /// error archive (`finish_archive`; `no_query`: the run has no query, NCBI's `Query is
    /// Empty!`). Each format's writer takes the run's bytes as the command line's report stream
    /// takes them: where LOSAT refuses the strategy file at NCBI's first request
    /// (`strategy::check_report_stream`), what a format was written before that point (NCBI's
    /// prolog) stays before the refusal, and a run that NCBI ends before it keeps what each format
    /// was written.
    pub(crate) fn run_records(
        &mut self,
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
        archive: Option<&crate::report::archive::ArchiveHandle>,
        no_query: bool,
        run: impl FnOnce(&mut ReportOutputs<'_>) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let result = match archive {
            Some(handle) => self.with_archive_diagnostics(handle.clone(), run),
            None => run(self),
        };
        let result = crate::report::archive::strategy::end_run(strategy, result);
        match archive {
            Some(handle) => self.finish_archive(handle, no_query, result),
            None => result,
        }
    }
}

/// The fields of a `ReportOutputs` other than its formats and its diagnostics stream
/// (`ReportOutputs::split`), which a view of the outputs borrows as they are.
struct ViewParts<'s, 'a> {
    other_outfmts: &'a [&'a str],
    later_warnings: &'s mut Option<&'a mut [Vec<u8>]>,
    hits: &'s mut Option<&'a mut HitsSink<'a>>,
    observer: &'s mut Option<&'a mut (dyn FormatObserver + Send)>,
    html_posts: &'s mut Option<&'a mut (dyn Write + Send)>,
}

impl ViewParts<'_, '_> {
    /// The outputs seen through `formats` (each format's `-outfmt` text, the sink that the step
    /// writes it through, and its refusal, moved into the view and back after it) and the
    /// diagnostics stream `diagnostics`, with the other fields as they are: the one rebuild of the
    /// outputs for a step that sees them through other sinks or another diagnostics stream (the
    /// BLAST archive's saved messages, `ReportOutputs::with_archive_diagnostics`). The view's
    /// refusals go back with `ReportOutputs::into_refusals` and `FormatOutput::take_back`.
    fn view<'v>(
        &'v mut self,
        formats: Vec<FormatOutput<'v>>,
        diagnostics: &'v mut (dyn Write + Send),
    ) -> ReportOutputs<'v> {
        ReportOutputs {
            other_outfmts: self.other_outfmts,
            formats,
            diagnostics,
            later_warnings: self.later_warnings.as_deref_mut(),
            // A later search's `-html` posts go where they went (`ReportOutputs::html_posts`).
            html_posts: self
                .html_posts
                .as_deref_mut()
                .map(|posts| posts as &mut (dyn Write + Send)),
            hits: self
                .hits
                .as_mut()
                .map(|hits| &mut **hits as &mut HitsSink<'_>),
            observer: self
                .observer
                .as_mut()
                .map(|observer| &mut **observer as &mut (dyn FormatObserver + Send)),
        }
    }
}

impl FormatOutput<'_> {
    /// Gives each of `formats` back the refusal of its view (`ReportOutputs::into_refusals`).
    fn take_back(formats: &mut [FormatOutput<'_>], refusals: Vec<FormatRefusal>) {
        for (format, refusal) in formats.iter_mut().zip(refusals) {
            format.refusal = refusal;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:225-226
// One command-line invocation constructs one formatter for its single `-outfmt`.
impl<'a> ReportOutputs<'a> {
    /// The command-line shape: one output format and warnings, no records or observer.
    pub fn single(
        outfmt: &'a (impl AsRef<[u8]> + ?Sized),
        sink: OutputSink<'a>,
        diagnostics: &'a mut (dyn Write + Send),
    ) -> Self {
        Self {
            other_outfmts: &[],
            html_posts: None,
            formats: vec![FormatOutput {
                outfmt: outfmt.as_ref(),
                sink,
                refusal: FormatRefusal::FailsCall,
            }],
            diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        }
    }
}

#[cfg(test)]
#[path = "local_blast_tests/mod.rs"]
mod local_blast_tests;
