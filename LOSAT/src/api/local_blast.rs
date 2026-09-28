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

use crate::report::PairwiseHit;

// Re-export TBLASTX run function
pub use crate::algorithm::tblastx::blast_engine::run as run_tblastx;

// Re-export BLASTN run function
pub use crate::algorithm::blastn::blast_engine::run as run_blastn;

// Re-export BLASTP run function
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:1196-1218
// ```c
// CLocalBlast::Run()
// {
//     x_SetupSearch();
//     x_RunSearch();
// }
// ```
pub use crate::algorithm::blastp::blast_engine::run as run_blastp;
pub use crate::algorithm::blastp::blast_engine::run_local as run_local_blastp;

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
pub enum OutputSink<'a> {
    /// Standard output, buffered (the CLI without `-out`, NCBI's default `-`).
    Stdout,
    /// A file named by `-out`, created when the report is written.
    File(&'a Path),
    /// A writer provided by the caller. It is not buffered here, so that every byte a
    /// formatter has written reaches it before each `FormatObserver` call.
    Writer(&'a mut (dyn Write + Send)),
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
            OutputSink::Writer(writer) => Box::new(&mut **writer),
        })
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
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3613-3630
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

/// Reports which HSP a formatter is writing, so that a caller can relate each HSP to its
/// bytes in every output. The written bytes do not change.
pub trait FormatObserver {
    /// The formatter for `ReportOutputs::formats[format]` starts the row (outfmt 6/7) or
    /// the alignment section (outfmt 0: the score lines and the alignment, without the
    /// subject defline) of `hsp`. Every byte written for that format before this call
    /// has reached its sink.
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex);
    /// The formatter for `ReportOutputs::formats[format]` finished the row or section of
    /// `hsp`. Every byte of it has reached the sink.
    fn hsp_end(&mut self, format: usize, hsp: HspIndex);
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
}

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
// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-465
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ITERATE(CSearchResultSet, result, *results) {
//         formatter.PrintOneResultSet(**result, queries);
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
    /// Receives the final, sorted hit list with rendered alignments.
    pub hits: Option<&'a mut (dyn FnMut(&[PairwiseHit]) + Send)>,
    /// Receives the row and section boundaries of each HSP in each format.
    pub observer: Option<&'a mut (dyn FormatObserver + Send)>,
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
