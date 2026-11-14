//! Local BLAST Search
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp
//!
//! This module provides the main entry point for running local BLAST searches, and the
//! output routing shared by every entry that runs through it (the CLI, web ABI v1 and
//! web ABI v2). A program routed through the shared entry exposes `run_local`, which
//! searches once and hands the result to every requested output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::report::blast_xml::XmlStats;
use crate::report::PairwiseHit;

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:106-111
// ```c
// if(RecoverSearchStrategy(args, m_CmdLineArgs)) {
// 	opts_hndl.Reset(&*m_CmdLineArgs->SetOptionsForSavedStrategy(args));
// }
// else {
// 	opts_hndl.Reset(&*m_CmdLineArgs->SetOptions(args));
// }
// ```
pub use crate::algorithm::tblastx::blast_engine::run_local as run_local_tblastx;
// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:113-119
// ```c++
//         /*** Initialize the database/subject ***/
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope;
//         InitializeSubject(db_args, opts_hndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
pub use crate::algorithm::tblastx::blast_engine::run_local_database as run_local_database_tblastx;

// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:128-133
// ```c
// if(RecoverSearchStrategy(args, m_CmdLineArgs)) {
// 	m_OptsHndl.Reset(&*m_CmdLineArgs->SetOptionsForSavedStrategy(args));
// }
// else {
// 	m_OptsHndl.Reset(&*m_CmdLineArgs->SetOptions(args));
// }
// ```
pub use crate::algorithm::blastn::blast_engine::run_local as run_local_blastn;
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
// ```c++
//         /*** Initialize the database/subject ***/
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope;
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
pub use crate::algorithm::blastn::blast_engine::run_local_database as run_local_database_blastn;

// Re-export the pool of a run of several searches (`run_local` calls of one argv), which they
// share instead of building one each (no NCBI counterpart: NCBI's CLI searches once per run).
pub use crate::utils::threading::with_run_pool;

// Re-export BLASTP run function
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:1196-1218
// ```c
// CLocalBlast::Run()
// {
//     x_SetupSearch();
//     x_RunSearch();
// }
// ```
pub use crate::algorithm::blastp::blast_engine::run_local as run_local_blastp;
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:194-199
// ```c++
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope;
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
//         _ASSERT(db_adapter && scope);
// ```
pub use crate::algorithm::blastp::blast_engine::run_local_database as run_local_database_blastp;

// Re-export the TBLASTN entry
// NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:289-301
// ```c
// CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
// lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
// results = lcl_blast.Run();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     formatter.PrintOneResultSet(**result, query);
// }
// ```
pub use crate::algorithm::tblastn::run_local as run_local_tblastn;
// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
// ```c++
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
pub use crate::algorithm::tblastn::run_local_database as run_local_database_tblastn;

// Re-export BLASTX's shared entry (web ABI v2).
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
pub use crate::algorithm::blastx::run_local as run_local_blastx;
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-200
// ```c++
//         /*** Initialize the database/subject ***/
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope;
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
pub use crate::algorithm::blastx::run_local_database as run_local_database_blastx;

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
// ```c
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
//     : m_FormatType(format_type), ..., m_Outfile(outfile),
//       m_NumSummary(num_summary), ...
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
// ```c
// arg_desc.AddDefaultKey(kArgOutput, "output_file",
//                "Output file name",
//                CArgDescriptions::eOutputFile, "-");
// ```
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
    /// A writer provided by the caller, as `Writer`, and the caller's receiver of the side files
    /// of outfmt 13 and 14 (a caller without files: LOSAT Web). The writer gets the master file,
    /// the receiver each side file, named as next to the receiver's `-out`
    /// (`SideFileSink::out`), which the format takes as its `-out`: its XInclude check passes.
    /// For any other format it is `Writer`.
    Files(
        &'a mut (dyn Write + Send),
        crate::report::blast_xml2::SideFileSink,
    ),
}

impl OutputSink<'_> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
    // ```c
    // else {
    //     m_OutputStream = &args[kArgOutput].AsOutputFile();
    // }
    // ```
    // NCBI opens the `-out` stream while it processes the arguments, before the search.
    // LOSAT has always created the file when the report is written, after the search;
    // this routing keeps that existing timing unchanged.
    /// Opens the sink for writing the report.
    pub fn open(&mut self) -> io::Result<Box<dyn Write + '_>> {
        Ok(match self {
            OutputSink::Stdout => Box::new(BufWriter::new(io::stdout())),
            OutputSink::File(path) => Box::new(BufWriter::new(File::create(path)?)),
            OutputSink::Writer(writer) | OutputSink::Files(writer, _) => Box::new(&mut **writer),
        })
    }

    /// Opens the sink as the stream of the format's report (`cli::NcbiOstream`, which the organism
    /// report of outfmt 18 takes): `open`'s writer, passed on as it is.
    pub(crate) fn ostream(&mut self) -> io::Result<crate::cli::NcbiOstream<'_>> {
        self.open().map(crate::cli::NcbiOstream::new)
    }

    /// The caller's receiver of the side files of outfmt 13 and 14 (`OutputSink::Files`).
    pub fn side_files(&self) -> Option<&crate::report::blast_xml2::SideFileSink> {
        match self {
            OutputSink::Files(_, sink) => Some(sink),
            _ => None,
        }
    }

    /// Whether the output is a writer of the caller (`Writer` or `Files`).
    pub(crate) fn is_writer(&self) -> bool {
        matches!(self, OutputSink::Writer(_) | OutputSink::Files(..))
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
// ```c
// if (args[kArgOutputFormat]) {
//     string fmt_choice =
//         NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
// ```
/// One requested output: an `-outfmt` value, parsed and validated by the program exactly
/// as on the command line, and where to write it.
pub struct FormatOutput<'a> {
    pub outfmt: &'a str,
    pub sink: OutputSink<'a>,
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
        use crate::blastinput::app::{parse_formatting_string, ReportFormat};
        parse_formatting_string(self.outfmt)
            .ok()
            .and_then(|choice| ReportFormat::from_number(choice.number))
    }

    /// Whether this output is the BLAST archive (outfmt 11).
    fn is_archive(&self) -> bool {
        self.report_format() == Some(crate::blastinput::app::ReportFormat::ArchiveFormat)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
// ```c
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
// ```
// Every output format prints from the same final result set, so a position in that
// result identifies an HSP in every format.
/// Index of an HSP in the final, sorted hit list of a search: the list that
/// `ReportOutputs::hits` receives. It identifies the HSP in every output format.
pub type HspIndex = usize;

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// subid=&(avRef->GetSeqId(1));
// bool showDefLine = previousId.Empty() || !subid->Match(*previousId);
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     // Add tab in front of field, except for the first field.
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3956-3975
// ```c++
// void CDisplaySeqalign::x_DisplayAlnvecInfo(CNcbiOstream& out,
//                                            SAlnInfo* aln_vec_info,
//                                            bool show_defline)
// {
//     ...
// 	if(!m_AlignTemplates) {
// 		x_ShowAlnvecInfo(out,aln_vec_info,show_defline);
// 	}
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3613-3632
// ```c++
// void CDisplaySeqalign::x_ShowAlnvecInfo(CNcbiOstream& out,
//                                            SAlnInfo* aln_vec_info,
//                                            bool show_defline)
// {
// 	bool showSortControls = false;
//     if(show_defline) {
// 		...
// 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
// 				out<< deflines;
// 		...
// 			out << "\n";
// ```
// The observed units follow NCBI's: one printed row per HSP in outfmt 6 and 7, and in
// outfmt 0 the output of one x_DisplayAlnvecInfo call without the subject defline that
// it prints before the first HSP of each subject (the section starts at " Score =").

// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
// ```c++
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// ```
/// A sequence of a database that a search reported (`FormatObserver::database_subjects`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportedSubject {
    /// The sequence's ordinal id in the database (the `n` of `gnl|BL_ORD_ID|n`), the same in
    /// every search of the database.
    pub oid: u32,
    /// The `saccver` of the default tabular report.
    pub id: std::sync::Arc<str>,
    /// Its letters.
    pub length: usize,
}

/// Reports which HSP a formatter is writing, so that a caller can relate each HSP to its
/// bytes in every output, and the statistics of each query of the search. The written bytes do
/// not change.
pub trait FormatObserver {
    /// The formatter for `ReportOutputs::formats[format]` starts the row (outfmt 6/7) or
    /// the alignment section (outfmt 0: the score lines and the alignment, without the
    /// subject defline) of `hsp`. Every byte written for that format before this call
    /// has reached its sink.
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex);
    /// The formatter for `ReportOutputs::formats[format]` finished the row or section of
    /// `hsp`. Every byte of it has reached the sink.
    fn hsp_end(&mut self, format: usize, hsp: HspIndex);

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3613-3632
    // ```c++
    // void CDisplaySeqalign::x_ShowAlnvecInfo(CNcbiOstream& out,
    //                                            SAlnInfo* aln_vec_info,
    //                                            bool show_defline)
    // {
    // 	bool showSortControls = false;
    //     if(show_defline) {
    // 		...
    // 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
    // 				out<< deflines;
    // 		...
    // 			out << "\n";
    // ```
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
    // ```c++
    //     // find the first valid context corresponding to this query
    //     for (i = 0; i < context_per_query; i++) {
    //         BlastContextInfo *ctx = query_info->contexts +
    //                                 query_number * context_per_query + i;
    //         if (ctx->is_valid) {
    //             m_SearchSpace = ctx->eff_searchsp;
    //             m_LengthAdjustment = ctx->length_adjustment;
    //             break;
    //         }
    //     }
    //     if (i >= context_per_query) {
    //         return; // we didn't find a valid context :(
    //     }
    // ```
    // NCBI keeps one `CBlastAncillaryData` per query with the query's result; here each program
    // hands over what its XML reports print of it, whichever formats were asked for.
    /// The statistics of query `q_idx` of the search: the values that the `Statistics` of its
    /// `Iteration` in outfmt 5 shows (`XmlStats::NONE` for a query without a valid context).
    /// `q_idx` counts as `Hit::q_idx` of the final hit list. The call comes once per query,
    /// before the reports of the query are written, and only with an observer (the program does
    /// not make these values otherwise).
    fn query_stats(&mut self, _q_idx: usize, _stats: XmlStats) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // ```c++
    //         formatter.PrintProlog();
    //         ...
    //         /*** Process the input ***/
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //         ...
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // What NCBI writes before this loop is written once for the run; a call of `run_local`
    // per query batch writes it again, and this event tells where it ends.
    /// The checks of the run are done and the subjects are read: what the diagnostics receive
    /// from here on belongs to the queries.
    fn queries_begin(&mut self) {}

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:299-300
    // ```c++
    //                 if (!batch_size)
    //                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
    // ```
    /// BLASTN's query batch just searched made `extensions` successful initial extensions, which
    /// size the next batch (`CBatchSizeMixer`). The other programs' batches have one size.
    fn batch_extensions(&mut self, _extensions: i32) {}

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
    // ```c++
    //        		const CBioseq_Handle& subject_bh =
    //                 scope.GetBioseqHandle(align.GetSeq_id(1));
    //             if(setSubjectId) {
    //             	SetSubjectId(subject_bh);
    //             }
    //             m_SubjectLength = subject_bh.GetBioseqLength();
    // ```
    // NCBI's reports read a subject's ids and length by the Seq-align's id; a database search's
    // hit list names its sequences by rank, and this event tells what the ranks are.
    /// The sequences that `Hit::s_idx` of a database search's final hit list (`ReportOutputs::hits`)
    /// names, by that index. The call comes once, with `ReportOutputs::hits` only (the program
    /// does not make these values otherwise).
    fn database_subjects(&mut self, _subjects: Vec<ReportedSubject>) {}
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
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

// NCBI reference: ncbi-blast/c++/src/app/blast/blastp_app.cpp:195-295
// ```c
// CRef<CLocalDbAdapter> db_adapter;
// CRef<CScope> scope;
// InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                  db_adapter, scope);
// CBlastFormat formatter(opt, *db_adapter,
//                        fmt_args->GetFormattedOutputChoice(), ...);
// formatter.PrintProlog();
// for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//         CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//         results = lcl_blast.Run();
//         ITERATE(CSearchResultSet, result, *results) {
//             formatter.PrintOneResultSet(**result, query_batch);
//         }
// }
// formatter.PrintEpilog(opt);
// ```
// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one CSearchResultSet without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// The outputs of one search. Every part is `Send`, because the search runs the
/// formatters inside its thread pool.
pub struct ReportOutputs<'a> {
    /// Every requested output format. The CLI requests exactly one.
    pub formats: Vec<FormatOutput<'a>>,
    /// Receives the warnings that the CLI writes to standard error, each once.
    pub diagnostics: &'a mut (dyn Write + Send),
    /// Receives the final, sorted hit list with rendered alignments, once per call. It may come
    /// before or after the reports are written (blastp's query blocks hand it over after the
    /// last block, `LOSAT_X_BPQBLOCK`). The list is shared, not copied: it holds an alignment's
    /// rows for every HSP, which the reports read while the caller keeps it for its records.
    pub hits: Option<&'a mut HitsSink<'a>>,
    /// Receives the row and section boundaries of each HSP in each format. The events of one
    /// format come in that format's write order; those of different formats may interleave in
    /// any order (blastp's query blocks write every format block by block).
    pub observer: Option<&'a mut (dyn FormatObserver + Send)>,
}

impl ReportOutputs<'_> {
    /// Whether every requested format is the BLAST archive (outfmt 11), which NCBI's apps but
    /// psiblast write without fetching the report's sequences ahead
    /// (`blastinput::app::check_search_environment`); a call without formats (hit records
    /// only) is not.
    pub(crate) fn only_archive(&self) -> bool {
        !self.formats.is_empty() && self.formats.iter().all(FormatOutput::is_archive)
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:913-916
    // ```c++
    // void PrintErrorArchive(const CArgs & a, const list<CRef<CBlast4_error> > & msg)
    // {
    // 	try {
    // 		if(NStr::StringToInt(NStr::TruncateSpaces(a[align_format::kArgOutputFormat].AsString())) == CFormattingArgs::eArchiveFormat) {
    // ```
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
    // ```c++
    //     } CATCH_ALL(status)
    //
    //     if(!m_Bah.GetMessages().empty()) {
    //     	const CArgs & a = GetArgs();
    //     	PrintErrorArchive(a, m_Bah.GetMessages());
    //     }
    // ```
    /// The end of a run with a BLAST archive (outfmt 11) that `result` closes: NCBI writes an
    /// error archive after the archives of the batches before when the application stops with one
    /// of its own errors, and when the run ends with messages that no archive carried. `no_query`:
    /// the run ended at `Query is Empty!`, where NCBI returns before its `PrintErrorArchive`.
    pub(crate) fn finish_archive(
        &mut self,
        handle: &crate::report::archive::ArchiveHandle,
        no_query: bool,
        result: anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let lines = match &result {
            // The errors of `Run`, the files' among them (`run_local_database` opens no file
            // that the argument handling opens).
            Err(error) => crate::report::archive::error_lines(error, true),
            Ok(()) if no_query => None,
            Ok(()) => handle.has_messages().then(Vec::new),
        };
        if let (Some(lines), Some(format)) = (
            lines,
            self.formats.iter_mut().find(|format| format.is_archive()),
        ) {
            let mut writer = format.sink.open()?;
            let written = crate::report::archive::write_error_archive(handle, &lines, &mut writer)
                .and_then(|()| writer.flush().map_err(anyhow::Error::from));
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
        let ReportOutputs {
            formats,
            diagnostics,
            hits,
            observer,
        } = self;
        let mut tee = crate::report::archive::DiagTee::new(&mut **diagnostics, handle);
        let mut view = ReportOutputs {
            formats: formats
                .iter_mut()
                .map(|format| FormatOutput {
                    outfmt: format.outfmt,
                    sink: match &mut format.sink {
                        OutputSink::Stdout => OutputSink::Stdout,
                        OutputSink::File(path) => OutputSink::File(path),
                        OutputSink::Writer(writer) => OutputSink::Writer(&mut **writer),
                        OutputSink::Files(writer, sink) => {
                            OutputSink::Files(&mut **writer, sink.clone())
                        }
                    },
                })
                .collect(),
            diagnostics: &mut tee,
            hits: hits.as_mut().map(|hits| &mut **hits as &mut HitsSink<'_>),
            observer: observer
                .as_mut()
                .map(|observer| &mut **observer as &mut (dyn FormatObserver + Send)),
        };
        f(&mut view)
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blastp_app.cpp:225-226
// ```c
// CBlastFormat formatter(opt, *db_adapter,
//                        fmt_args->GetFormattedOutputChoice(),
// ```
// One command-line invocation constructs one formatter for its single `-outfmt`.
impl<'a> ReportOutputs<'a> {
    /// The command-line shape: one output format and warnings, no records or observer.
    pub fn single(
        outfmt: &'a str,
        sink: OutputSink<'a>,
        diagnostics: &'a mut (dyn Write + Send),
    ) -> Self {
        Self {
            formats: vec![FormatOutput { outfmt, sink }],
            diagnostics,
            hits: None,
            observer: None,
        }
    }
}

#[cfg(test)]
#[path = "local_blast_tests/mod.rs"]
mod local_blast_tests;
