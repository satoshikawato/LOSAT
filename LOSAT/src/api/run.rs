//! Running a search from memory: `run` writes each report to the caller's writer, `run_collect`
//! gathers the reports, the HSP records, the statistics of each query and the diagnostics.
//!
//! NCBI's application searches once and prints the result set with one formatter per
//! requested output; LOSAT's `run_local` entries do the same for records in memory, and this
//! module hands them the caller's inputs and outputs.
//!
//! NCBI reference: ncbi-blast/c++/src/app/blast/blastp_app.cpp:195-295
//! ```c
//! CRef<CLocalDbAdapter> db_adapter;
//! CRef<CScope> scope;
//! InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//!                  db_adapter, scope);
//! CBlastFormat formatter(opt, *db_adapter,
//!                        fmt_args->GetFormattedOutputChoice(), ...);
//! formatter.PrintProlog();
//! for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!         CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//!         results = lcl_blast.Run();
//!         ITERATE(CSearchResultSet, result, *results) {
//!             formatter.PrintOneResultSet(**result, query_batch);
//!         }
//! }
//! formatter.PrintEpilog(opt);
//! ```
#![warn(missing_docs, unused)]

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::error::Error;
use super::input::{check_archive_inputs, is_protein, Inputs, Role};
use super::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_database_blastn,
    run_local_database_blastp, run_local_database_blastx, run_local_database_tblastn,
    run_local_database_tblastx, run_local_tblastn, run_local_tblastx, with_run_pool,
    FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs, ReportedSubject,
};
use super::program::Program;
use super::search::Search;
use super::start::{CallStart, RunLog, RunSink};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::{
    parse_optional_range, record_interval, RangeRole, RecordInterval, SequenceRange,
};
use crate::cli::Commands;
use crate::report::blast_xml::XmlStats;
use crate::report::PairwiseHit;

/// One report of a search: the bytes that the command line writes to `-out` for this
/// `-outfmt`, and where each HSP of `Outcome::hsps` lies in them.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Report {
    /// The `-outfmt` value, as requested.
    pub outfmt: String,
    /// The report's bytes.
    pub bytes: Vec<u8>,
    /// By `HspRecord::index`: the byte range in `bytes` of the HSP's row in a tabular format or
    /// of its section in the pairwise format (its score lines and alignment, without the
    /// subject's heading). `None` where the report does not show the HSP or came from another
    /// search than the records (`run_collect`).
    pub spans: Vec<Option<Range<usize>>>,
    /// By `HspRecord::index`: the byte range in `bytes` of the heading of the HSP's subject (the
    /// defline and `Length=` line that the pairwise format writes before the first HSP of each
    /// subject); every HSP of the subject in that query has the same range. `None` for the
    /// formats without headings.
    pub subject_headings: Vec<Option<Range<usize>>>,
}

/// A strand of a nucleotide sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// The sequence as read.
    Plus,
    /// Its reverse complement.
    Minus,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1061-1067
// ```c++
//         if(x_IsFieldRequested(eSubjectStrand))
//         {
//         	if(!subject_is_na)
//         		m_SubjectStrand = NA;
//         	else
//         		m_SubjectStrand = ((s_start - s_end) > 0 )? "minus":"plus";
//         }
// ```
/// The strand that NCBI's `sstrand` column gives a row of a nucleotide sequence, and no strand
/// for a protein: an end before the start is the minus strand.
fn strand(protein: bool, start: usize, end: usize) -> Option<Strand> {
    (!protein).then_some(if start > end {
        Strand::Minus
    } else {
        Strand::Plus
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
// ```c
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
// ```
// Every output format prints from the same final result set, so a position in that result
// identifies an HSP in every format (`HspIndex`).
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-132
// ```c
// typedef struct BlastHSP {
//    Int4 score;           /**< This HSP's raw score */
//    Int4 num_ident;       /**< Number of identical base pairs in this HSP */
//    double bit_score;     /**< Bit score, calculated from score */
//    double evalue;        /**< This HSP's e-value */
//    BlastSeg query;       /**< Query sequence info. */
//    BlastSeg subject;     /**< Subject sequence info. */
// ```
/// One HSP of the final, sorted hit list of a search, with the engine's values. Each field that
/// names a column of `-outfmt 6` holds the number that the column prints for the HSP's row, not
/// rounded (`tests/library_api.rs` holds them to the columns); the text of a report is its
/// tabular row (`Report::spans`), not the values here formatted again.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct HspRecord {
    /// 0-based position of the HSP in the final, sorted hit list; the same HSP has the same
    /// index in every report (`Report::spans`).
    pub index: usize,
    /// 0-based index of the query record in the inputs (records without residues count; with
    /// `-query_loc`, a query record that NCBI skips because its interval starts past its end
    /// keeps its place).
    pub q_idx: usize,
    /// 0-based index of the subject record in the inputs; for a database
    /// (`Inputs::database`), the sequence's OID (the `n` of `gnl|BL_ORD_ID|n`).
    pub s_idx: usize,
    /// 0-based position of the HSP among the HSPs of its query.
    pub rank: usize,
    /// The query's id (`qaccver`): the first word of its title, its generated `Query_N` without a
    /// title, and with `-parse_deflines` the accession and version of its id. A byte that is not
    /// UTF-8 shows as U+FFFD here and as itself in the reports. The HSPs of a query share it.
    pub query_id: Arc<str>,
    /// The subject's id (`saccver`), by the same rules. The HSPs of a subject share it.
    pub subject_id: Arc<str>,
    /// The raw alignment score (`score`).
    pub raw_score: i32,
    /// The bit score (`bitscore`), not rounded.
    pub bit_score: f64,
    /// The E-value (`evalue`), not rounded. It is a value of the whole search: the same pair
    /// has another E-value among other subjects (`QueryStats::eff_space`).
    pub e_value: f64,
    /// 1-based query start (`qstart`).
    pub q_start: usize,
    /// 1-based query end (`qend`): before `q_start` on the minus strand.
    pub q_end: usize,
    /// 1-based subject start (`sstart`).
    pub s_start: usize,
    /// 1-based subject end (`send`): before `s_start` on the minus strand.
    pub s_end: usize,
    /// The strand of the query row: `None` for a protein query. BLASTN shows every query row on
    /// the plus strand and turns the subject instead, as NCBI's reports do.
    pub query_strand: Option<Strand>,
    /// The strand of the subject row (`sstrand`): `None` for a protein subject. A row of one
    /// letter has its start at its end, and NCBI's rule calls it plus.
    pub subject_strand: Option<Strand>,
    /// The query's translation frame (`qframe`): `None` where the program does not translate
    /// the query.
    pub query_frame: Option<i32>,
    /// The subject's translation frame (`sframe`): `None` where the program does not translate
    /// the subject.
    pub subject_frame: Option<i32>,
    /// The query's length (`qlen`) in the letters of the record, nucleotides for a translated
    /// query: with `-query_loc` the length of the interval as given, even where it ends past
    /// the record.
    pub query_length: usize,
    /// The subject's length (`slen`) in the letters of the record or the database sequence,
    /// nucleotides for a translated subject: the whole record's, also with `-subject_loc`.
    pub subject_length: usize,
    /// The alignment's columns, gaps included (`length`).
    pub length: usize,
    /// The columns of identical letters (`nident`).
    pub num_ident: usize,
    /// The columns of a positive score (`positive`): the identities in a nucleotide search,
    /// which has no scoring matrix.
    pub num_positives: usize,
    /// The columns of different letters (`mismatch`): `length - num_ident - gaps`, as NCBI
    /// computes it. Below zero where a subject's stored gap character faces a gap of the query,
    /// which counts as an identity and as a gap.
    pub mismatch: i64,
    /// The gap openings (`gapopen`).
    pub gap_open: usize,
    /// The gap columns (`gaps`).
    pub gaps: usize,
    /// How many HSPs NCBI's sum statistics linked with this one (`Expect(n)` of the pairwise
    /// format), where the program links HSPs; `None` otherwise. 0 and 1 both mean an HSP on its
    /// own (`BlastHSP::num`).
    pub sum_n: Option<i32>,
    /// The composition-based statistics method that scored the HSP, where the program has one:
    /// 0 none, 1 `Composition-based stats`, 2 `Compositional matrix adjust`.
    pub comp_adjust_method: Option<u8>,
    /// The aligned query row, `-` for gaps, where the search renders it.
    pub query_aligned: Option<String>,
    /// The aligned subject row, `-` for gaps, where the search renders it.
    pub subject_aligned: Option<String>,
}

impl HspRecord {
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:437-442
    // ```c++
    // inline void CBlastTabularInfo::x_PrintPercentIdentical(void)
    // {
    //     double perc_ident =
    //         (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
    //     m_Ostream << NStr::DoubleToString(perc_ident, 3);
    // }
    // ```
    /// The percentage of identical columns (`pident`), not rounded.
    pub fn percent_identity(&self) -> f64 {
        crate::report::outfmt6::percent_identity_ncbi(self.num_ident, self.length)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1109-1120
    // ```c
    // Blast_HSPGetAdjustedOffsets(EBlastProgramType program, BlastHSP* hsp,
    //                             Int4 query_length, Int4 subject_length,
    //                             Int4* q_start, Int4* q_end,
    //                             Int4* s_start, Int4* s_end)
    // {
    //    if (!hsp->gap_info) {
    //       *q_start = hsp->query.offset + 1;
    //       *q_end = hsp->query.end;
    //       *s_start = hsp->subject.offset + 1;
    //       *s_end = hsp->subject.end;
    //       return;
    //    }
    // ```
    // The reports' coordinates are 1-based offsets of a `BlastSeg`, whose `offset` and `end` are
    // a 0-based half-open interval; the methods below give that form, on the plus strand.
    /// The query interval as a 0-based half-open range on the record's plus strand, whichever
    /// strand the row reads: `q_start` and `q_end` in the order of the record.
    pub fn query_range(&self) -> Range<usize> {
        self.q_start.min(self.q_end).saturating_sub(1)..self.q_start.max(self.q_end)
    }

    /// The subject interval as a 0-based half-open range on the record's plus strand, whichever
    /// strand the row reads: `s_start` and `s_end` in the order of the record.
    pub fn subject_range(&self) -> Range<usize> {
        self.s_start.min(self.s_end).saturating_sub(1)..self.s_start.max(self.s_end)
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_results.hpp:115-122
// ```c++
//     /// Retrieve the search space for this query sequence. If the
//     /// results correspond to a blastx search, the search space will
//     /// refer to protein letters
//     Int8 GetSearchSpace() const {
//         return m_SearchSpace;
//     }
//     /// Retrieve the length adjustment for boundary conditions
//     Int8 GetLengthAdjustment() const {
// ```
/// The statistics of one query of a search (NCBI's `CBlastAncillaryData`): the values that the
/// `Statistics` of the query's `Iteration` in outfmt 5 shows. The E-values of the query's HSPs
/// follow from them and the raw scores, so they say what an E-value of this search means.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct QueryStats {
    /// The effective search space (`Statistics_eff-space`), in protein letters where the query
    /// or the subject is translated.
    pub eff_space: i64,
    /// The length adjustment (`Statistics_hsp-len`).
    pub hsp_len: i32,
    /// Karlin-Altschul lambda (`Statistics_lambda`).
    pub lambda: f64,
    /// Karlin-Altschul K (`Statistics_kappa`).
    pub kappa: f64,
    /// Karlin-Altschul H (`Statistics_entropy`).
    pub entropy: f64,
}

impl QueryStats {
    /// The statistics of a query, or `None` for the values of a query without a valid context.
    pub(super) fn of_context(stats: XmlStats) -> Option<Self> {
        (stats != XmlStats::NONE).then_some(Self {
            eff_space: stats.eff_space,
            hsp_len: stats.hsp_len,
            lambda: stats.lambda,
            kappa: stats.kappa,
            entropy: stats.entropy,
        })
    }
}

/// What `run_collect` gathers from a search.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Outcome {
    /// One report per requested format, in the requested order.
    pub reports: Vec<Report>,
    /// The HSPs of the final hit list. Their places in the reports are in `Report::spans`.
    pub hsps: Vec<HspRecord>,
    /// One entry per query record of the inputs, in their order (`HspRecord::q_idx`). `None` for
    /// a query that has no statistics: no valid context (no residues, all `N` or all `X`), or
    /// not searched (its `-query_loc` interval starts past its end, or NCBI searches no query of
    /// its batch).
    pub query_stats: Vec<Option<QueryStats>>,
    /// What the command line writes to its standard error for the search: the warnings of the
    /// reader and of the search, each once. Empty for a search without warnings.
    pub diagnostics: Vec<u8>,
}

/// A writer into memory that tells its position, for the observer of its report.
pub(super) struct CountingWriter {
    pub(super) bytes: Vec<u8>,
    pub(super) position: Arc<AtomicUsize>,
}

impl CountingWriter {
    pub(super) fn new() -> Self {
        Self {
            bytes: Vec::new(),
            position: Arc::new(AtomicUsize::new(0)),
        }
    }
}

impl Write for CountingWriter {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.bytes.extend_from_slice(data);
        self.position.fetch_add(data.len(), Ordering::SeqCst);
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

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
/// Records, for every report of one search, the byte range of each HSP's row or section and
/// of the heading of each HSP's subject, from the formatter's observer events.
struct RangeRecorder {
    positions: Vec<Arc<AtomicUsize>>,
    /// Per report: the HSP whose row or section is being written, and where it started.
    open: Vec<Option<(HspIndex, usize)>>,
    /// Per report, per HSP: the row or section.
    spans: Vec<Vec<Option<Range<usize>>>>,
    /// Per report, per HSP: the heading of the HSP's subject.
    headings: Vec<Vec<Option<Range<usize>>>>,
    /// Per report: the start of the heading being written, and the last complete heading.
    heading_start: Vec<Option<usize>>,
    heading: Vec<Option<Range<usize>>>,
    /// By the engine's query number: the query's statistics.
    query_stats: Vec<Option<QueryStats>>,
    /// The sequences that a database search's hit list names.
    database_subjects: Option<Vec<ReportedSubject>>,
    /// The first observer event that does not fit the formatter contract; the search fails.
    fault: Option<String>,
}

impl RangeRecorder {
    fn new(positions: Vec<Arc<AtomicUsize>>) -> Self {
        let reports = positions.len();
        Self {
            positions,
            open: vec![None; reports],
            spans: vec![Vec::new(); reports],
            headings: vec![Vec::new(); reports],
            heading_start: vec![None; reports],
            heading: vec![None; reports],
            query_stats: Vec::new(),
            database_subjects: None,
            fault: None,
        }
    }

    fn position(&self, report: usize) -> usize {
        self.positions[report].load(Ordering::SeqCst)
    }

    fn fault(&mut self, message: String) {
        self.fault.get_or_insert(message);
    }
}

/// Stores `range` for `hsp`; a second range for the same HSP is a fault.
fn store(
    ranges: &mut Vec<Option<Range<usize>>>,
    hsp: HspIndex,
    range: Range<usize>,
) -> Result<(), String> {
    if ranges.len() <= hsp {
        ranges.resize(hsp + 1, None);
    }
    match ranges[hsp].replace(range) {
        None => Ok(()),
        Some(_) => Err(format!("HSP {hsp} was reported twice in one format")),
    }
}

impl FormatObserver for RangeRecorder {
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex) {
        let start = self.position(format);
        if let Some((open, _)) = self.open[format].replace((hsp, start)) {
            self.fault(format!("HSP {hsp} started before HSP {open} ended"));
        }
    }

    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        let Some((begun, start)) = self.open[format].take() else {
            self.fault(format!("HSP {hsp} ended without a start"));
            return;
        };
        if begun != hsp {
            self.fault(format!("HSP {begun} started but HSP {hsp} ended"));
            return;
        }
        let end = self.position(format);
        let stored = store(&mut self.spans[format], hsp, start..end).and_then(|()| {
            match self.heading[format].clone() {
                Some(heading) => store(&mut self.headings[format], hsp, heading),
                None => Ok(()),
            }
        });
        if let Err(fault) = stored {
            self.fault(fault);
        }
    }

    fn subject_begin(&mut self, format: usize, _first_hsp: HspIndex) {
        self.heading_start[format] = Some(self.position(format));
    }

    fn subject_end(&mut self, format: usize, _first_hsp: HspIndex) {
        if let Some(start) = self.heading_start[format].take() {
            self.heading[format] = Some(start..self.position(format));
        }
    }

    fn query_stats(&mut self, q_idx: usize, stats: XmlStats) {
        if self.query_stats.len() <= q_idx {
            self.query_stats.resize(q_idx + 1, None);
        }
        self.query_stats[q_idx] = QueryStats::of_context(stats);
    }

    fn database_subjects(&mut self, subjects: Vec<ReportedSubject>) {
        self.database_subjects = Some(subjects);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155
// ```c++
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
//         catch (const exception&) {
//             continue; //SB-2307. ignore well formed, not found accession
//         }
// ```
/// The index in the inputs of each query that the engine searches: NCBI's batch reader skips a
/// query whose `-query_loc` interval starts more than one letter past its end
/// (`seq_range::cut_queries`), so the engine's `q_idx` counts the searched queries, and the
/// records give the index in the inputs. `None` without `-query_loc` (every query is searched).
pub(super) fn searched_queries(
    range: Option<SequenceRange>,
    records: &[FastaRecord],
) -> Option<Vec<usize>> {
    let range = range?;
    Some(
        records
            .iter()
            .enumerate()
            .filter(|(_, record)| {
                record_interval(&range, record.seq().len()) != RecordInterval::PastEnd
            })
            .map(|(index, _)| index)
            .collect(),
    )
}

/// The `-query_loc` of the program's arguments.
fn query_loc(command: &Commands) -> Option<&str> {
    match command {
        Commands::Blastp(a) => a.query_loc.as_deref(),
        Commands::Tblastn(a) => a.query_loc.as_deref(),
        Commands::Blastn(a) => a.query_loc.as_deref(),
        Commands::Tblastx(a) => a.query_loc.as_deref(),
        Commands::Blastx(a) => a.query_loc.as_deref(),
        Commands::Psiblast(a) => a.query_loc.as_deref(),
        Commands::Rpsblast(a) => a.query_loc.as_deref(),
        Commands::Rpstblastn(a) => a.query_loc.as_deref(),
        Commands::Makedb(_) | Commands::Rpsbproc(_) => None,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// ```c
// int num_threads = args[kArgNumThreads].AsInteger();
// if (num_threads > kMaxValue) {
//     m_NumThreads = kMaxValue;
// } else {
//     m_NumThreads = num_threads;
// }
// ```
/// The `-query_loc` interval of the program's arguments. A value that the search will refuse
/// has no interval here.
pub(super) fn query_range_of(command: &Commands) -> Option<SequenceRange> {
    parse_optional_range(query_loc(command), RangeRole::Query)
        .ok()
        .flatten()
}

/// `-num_threads` of the program's arguments.
fn num_threads(command: &Commands) -> usize {
    match command {
        Commands::Blastp(a) => a.num_threads,
        Commands::Tblastn(a) => a.num_threads,
        Commands::Blastn(a) => a.num_threads,
        Commands::Tblastx(a) => a.num_threads,
        // A positive `i32` (the parser's `positive`).
        Commands::Blastx(a) => a.num_threads as usize,
        Commands::Psiblast(a) => a.num_threads,
        Commands::Rpsblast(a) => usize::try_from(a.num_threads.max(1)).unwrap_or(1),
        Commands::Rpstblastn(a) => usize::try_from(a.num_threads.max(1)).unwrap_or(1),
        Commands::Makedb(_) | Commands::Rpsbproc(_) => 1,
    }
}

/// The options that set the hit list size: `-max_target_seqs`, `-num_descriptions` and
/// `-num_alignments` of the program's arguments.
fn count_options(command: &Commands) -> [Option<usize>; 3] {
    match command {
        Commands::Blastp(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Tblastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Blastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Tblastx(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        // A positive `i32` (the parser's `positive`).
        Commands::Blastx(a) => [
            a.max_target_seqs.map(|n| n as usize),
            a.num_descriptions,
            a.num_alignments,
        ],
        Commands::Psiblast(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Rpsblast(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Rpstblastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Makedb(_) | Commands::Rpsbproc(_) => [None; 3],
    }
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973
// ```c
//     int hitlist_size = opt.GetHitlistSize();
//     if(m_OutputFormat <= eFlatQueryAnchoredNoIdentities) {
//     	 m_NumDescriptions = m_DfltNumDescriptions;
//     	 m_NumAlignments = m_DfltNumAlignments;
//     	 if (args.Exist(kArgNumDescriptions) && args[kArgNumDescriptions]) {
//     	    m_NumDescriptions = args[kArgNumDescriptions].AsInteger();
//     	 }
//          if (args.Exist(kArgNumAlignments) && args[kArgNumAlignments]) {
//     		m_NumAlignments = args[kArgNumAlignments].AsInteger();
//     	}
// ```
/// The hit list size that the command line derives for the output format `number` from the
/// options of `command` (`report_counts`).
fn hit_list_size(command: &Commands, number: i32) -> usize {
    let [max_target_seqs, descriptions, alignments] = count_options(command);
    crate::blastinput::report_counts::report_counts(
        number,
        max_target_seqs,
        descriptions,
        alignments,
    )
    .hitlist_size
}

/// The searches of a run: the requested output formats grouped by hit list size, each group a
/// list of positions in `numbers`, the tabular formats first. One group is one search
/// (`run_local` call), and every format of it is written from that search, as the command line
/// writes it alone: the command line's search depends on the hit list size, which
/// `-num_descriptions` and `-num_alignments` (without `-max_target_seqs`) can make differ
/// between the pairwise formats (0 to 4) and the tabular ones. Without those options the sizes
/// are equal and the one group is every format in the requested order. The tabular group goes
/// first because the warnings that the formats share belong to the first search's diagnostics,
/// and the command line writes the tabular formats' (`ReportCounts::warnings`).
///
/// LOSAT Web ABI v2 groups its searches the same way (`web/adapter/src/run.rs`).
pub(crate) fn search_groups(command: &Commands, numbers: &[i32]) -> Vec<Vec<usize>> {
    let [_, descriptions, alignments] = count_options(command);
    if descriptions.is_none() && alignments.is_none() {
        return vec![(0..numbers.len()).collect()];
    }
    let mut ordered: Vec<usize> = (0..numbers.len()).collect();
    ordered.sort_by_key(|&position| numbers[position] <= 4);
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for position in ordered {
        let size = hit_list_size(command, numbers[position]);
        match groups.iter_mut().find(|(group, _)| *group == size) {
            Some((_, group)) => group.push(position),
            None => groups.push((size, vec![position])),
        }
    }
    groups.into_iter().map(|(_, group)| group).collect()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// ```c++
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044
// ```c
//     // Setup some debugging features from environment variables.
//     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
//         SetDiagTrace(eDT_Enable, eDT_Enable);
//     }
//     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-220
// ```c++
// bool CSeqDBTaxInfo::GetTaxNames(TTaxId           tax_id,
//                                 SSeqDBTaxInfo  & info  )
// {
// 	static CTaxDBFileInfo t;
// ```
// The organism report (outfmt 18) looks the taxonomy database up itself, through the search
// path of NCBI's databases: a search of records reads the application settings for it, as the
// command line does (`common::database::SubjectTaxdb`).
/// The taxonomy database of `name`'s organism report when `outputs` has one (outfmt 18), found
/// with the call's settings (`start`).
fn subject_taxdb(
    name: &str,
    outputs: &mut ReportOutputs<'_>,
    start: &mut CallStart,
) -> anyhow::Result<Option<crate::algorithm::common::database::SubjectTaxdb>> {
    let tax_report = outputs.formats.iter().any(|format| {
        crate::blastinput::app::parse_formatting_string(format.outfmt)
            .is_ok_and(|choice| choice.number == 18)
    });
    if !tax_report {
        return Ok(None);
    }
    let settings = start.settings(&mut *outputs.diagnostics)?;
    crate::algorithm::common::database::SubjectTaxdb::of(&settings, &name.to_ascii_uppercase())
}

// NCBI's application layer reads its environment and registry files before its own code runs,
// and a database is found by them (the `BLASTDB` search path); LOSAT rejects the settings that
// change the output. A search of records reads neither, but for the directory of BLASTN's
// `-window_masker_taxid`, the taxonomy database of the organism report (outfmt 18) and the file
// of a `-matrix` that NCBI's tables do not have (TBLASTX, ungapped BLASTX), which NCBI's library
// looks for itself. The settings are the call's application start's
// (`start::CallStart::settings`).
/// One `run_local` call: the program's entry for the parsed arguments and `inputs`.
pub(super) fn run_local(
    command: Commands,
    inputs: Inputs<'_>,
    outputs: &mut ReportOutputs<'_>,
    start: &mut CallStart,
) -> anyhow::Result<()> {
    let (query, subject) = (inputs.query, inputs.subject);
    if inputs.database {
        let settings = start.settings(&mut *outputs.diagnostics)?;
        return match command {
            Commands::Blastp(args) => run_local_database_blastp(args, query, settings, outputs),
            Commands::Tblastn(args) => run_local_database_tblastn(args, query, settings, outputs),
            Commands::Blastn(args) => run_local_database_blastn(args, query, settings, outputs),
            Commands::Tblastx(args) => run_local_database_tblastx(args, query, settings, outputs),
            Commands::Blastx(args) => run_local_database_blastx(args, query, settings, outputs),
            Commands::Psiblast(args) => {
                let args = psiblast_args(args, outputs)?;
                run_local_database_blastp(args, query, settings, outputs)
            }
            Commands::Rpsblast(_)
            | Commands::Rpstblastn(_)
            | Commands::Makedb(_)
            | Commands::Rpsbproc(_) => {
                unreachable!("Program::parse admits the search programs only")
            }
        };
    }
    match command {
        Commands::Blastp(mut args) => {
            args.subject_taxdb = subject_taxdb("blastp", outputs, start)?;
            run_local_blastp(args, query, subject, "", "", outputs)
        }
        Commands::Tblastn(mut args) => {
            args.subject_taxdb = subject_taxdb("tblastn", outputs, start)?;
            args.tabular_taxdb = args.subject_taxdb.as_ref().map(|_| true);
            run_local_tblastn(args, query, subject, outputs)
        }
        Commands::Blastn(mut args) => {
            args.subject_taxdb = subject_taxdb("blastn", outputs, start)?;
            args.tabular_taxdb = args.subject_taxdb.as_ref().map(|_| true);
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-448
            // ```c++
            //     if (args.Exist(kArgWindowMaskerDatabase) &&
            //         args[kArgWindowMaskerDatabase]) {
            //         const string& stat_file = args[kArgWindowMaskerDatabase].AsString();
            // ```
            // The counts file of `-window_masker_db` is read from the file system, as on the
            // command line (LOSAT Web's host gives its own source).
            if args.window_masker_db.is_some() {
                args.window_masker_counts = Some(Arc::new(crate::utils::winmask::FileCounts));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:381-388
            // ```c++
            //     CRef<CNcbiRegistry> reg(new CNcbiRegistry(empty_stream,
            //                                               IRegistry::fWithNcbirc));
            //     ...
            //     reg->Add(*env_reg, CNcbiRegistry::ePriority_MaxUser);
            //     retval = reg->Get(kSection, kEnvVar);
            // ```
            // The directory of `-window_masker_taxid` is looked for by NCBI's library itself, in
            // the environment and `.ncbirc`; LOSAT reads it with the application settings, as
            // a database search of the library does.
            if args.window_masker_taxid.is_some_and(|taxid| taxid != 0) {
                let settings = start.settings(&mut *outputs.diagnostics)?;
                args.window_masker_registry = Some(settings.window_masker_registry);
                args.window_masker_counts = Some(Arc::new(crate::utils::winmask::FileCounts));
            }
            run_local_blastn(args, query, subject, outputs)
        }
        Commands::Tblastx(mut args) => {
            // The file of a `-matrix` that NCBI's tables do not have, first, as the command
            // line's `run` looks it up (with the call's settings; `TblastxArgs::find_matrix_file`).
            if args.scoring_matrix() == crate::algorithm::tblastx::args::TblastxMatrix::Unknown {
                let settings = start.settings(&mut *outputs.diagnostics)?;
                args.find_matrix_file(&settings.matrix_search)?;
            }
            args.subject_taxdb = subject_taxdb("tblastx", outputs, start)?;
            args.tabular_taxdb = args.subject_taxdb.as_ref().map(|_| true);
            run_local_tblastx(args, query, subject, outputs)
        }
        Commands::Blastx(mut args) => {
            args.subject_taxdb = subject_taxdb("blastx", outputs, start)?;
            // The file of a `-matrix` that NCBI's tables do not have, looked up after the options
            // are resolved, as on the command line (`ResolvedOptions::find_matrix_file`).
            if args.unknown_matrix_name() {
                args.matrix_search = Some(start.settings(&mut *outputs.diagnostics)?.matrix_search);
            }
            run_local_blastx(args, query, subject, outputs)
        }
        Commands::Psiblast(args) => {
            let args = psiblast_args(args, outputs)?;
            run_local_blastp(args, query, subject, "", "", outputs)
        }
        Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {
            unreachable!("Program::parse admits the search programs only")
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:233-241
// ```c++
//     CRef<IQueryFactory> query_factory(s_MakeQueryFactory(pssm, query));
//             CRef<CPsiBlast> psiblast;
//             if (run_token & 2)
//             {
//                 _ASSERT(psi_opts.NotEmpty());
//                 if (pssm.NotEmpty()) {
//                     psiblast.Reset(new CPsiBlast(pssm, db_adapter, psi_opts));
//                 } else {
//                     psiblast.Reset(new CPsiBlast(query_factory, db_adapter, psi_opts));
// ```
/// blastp's arguments for psiblast's first round (`algorithm::psiblast`), after LOSAT's refusals
/// of each requested format.
fn psiblast_args(
    args: crate::algorithm::psiblast::PsiblastArgs,
    outputs: &ReportOutputs<'_>,
) -> anyhow::Result<crate::algorithm::blastp::args::BlastpArgs> {
    for format in &outputs.formats {
        crate::algorithm::psiblast::PsiblastArgs::check_outfmt(format.outfmt)?;
    }
    args.into_blastp_args()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
// ```c
// if (args[kArgOutputFormat]) {
//     string fmt_choice =
//         NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
// ```
// No NCBI counterpart: NCBI's application reads its records itself, with the option.
/// The refusal of records read without their ids for a search with `-parse_deflines`: the
/// reports would show the local ids of a search without the option.
pub(super) fn check_inputs(search: &Search, inputs: Inputs<'_>) -> Result<(), Error> {
    let records = || inputs.query.iter().chain(inputs.subject);
    if search.parse_deflines()
        && records().next().is_some()
        && records().all(|record| record.believed().is_none())
    {
        return Err(Error::Argv(
            "-parse_deflines needs the records of read_inputs, which reads the ids".to_string(),
        ));
    }
    Ok(())
}

/// The searches of a call: `formats` grouped by hit list size (`search_groups`), and the threads
/// of the pool they share. A `-outfmt` value that the command line rejects, and a call without
/// formats (the engines refuse one), fail here, before any search; the command line's error of
/// a value goes to tblastx's `-logfile` file, as there (`RunLog::post`).
fn searches(
    search: &Search,
    inputs: Inputs<'_>,
    formats: &[&str],
    log: &mut RunLog,
) -> Result<(Vec<Vec<usize>>, usize), Error> {
    use crate::blastinput::app::output_format::ARCHIVE_FORMAT;
    use crate::blastinput::app::parse_formatting_string;
    if formats.is_empty() {
        return Err(Error::Argv("no output format was requested".to_string()));
    }
    check_inputs(search, inputs)?;
    let numbers: Vec<i32> = formats
        .iter()
        .map(|outfmt| {
            parse_formatting_string(outfmt)
                .map(|choice| choice.number)
                .map_err(|error| {
                    log.post(&error);
                    Error::Argv(format!("{error:#}"))
                })
        })
        .collect::<Result<_, _>>()?;
    if numbers.contains(&ARCHIVE_FORMAT) {
        check_archive_inputs(search.program(), inputs)?;
    }
    let command = search.command(inputs)?;
    Ok((search_groups(&command, &numbers), num_threads(&command)))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:286-288
// ```c++
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
// The formatters print to their streams as each result comes; here the streams are the caller's.
/// Runs the search and writes each report to its writer as the search makes it, and the command
/// line's warnings to `diagnostics`. `reports` holds an `-outfmt` value and a writer per report.
/// Each writer receives the bytes that the command line writes to `-out` for that `-outfmt` alone.
/// The formats are grouped as in `run_collect`, which says when that takes two searches and whose
/// warnings the diagnostics are. The writers are not buffered here, so a file's writer should be a
/// `BufWriter`.
///
/// The error is the command line's: the message of its parser or of the engine
/// (`Error::Search`). The search writes nothing to the process's standard output or error, and
/// never ends the process. What was written before an error stays in the writers.
pub fn run(
    search: &Search,
    inputs: Inputs<'_>,
    reports: &mut [(&str, &mut (dyn Write + Send))],
    diagnostics: &mut (dyn Write + Send),
) -> Result<(), Error> {
    let mut start = CallStart::new(search, inputs, &mut *diagnostics)?;
    if start.dry_run() {
        return Ok(());
    }
    let mut log = start.take_log();
    let formats: Vec<&str> = reports.iter().map(|(outfmt, _)| *outfmt).collect();
    let (groups, threads) = searches(search, inputs, &formats, &mut log)?;
    let mut writers: Vec<Option<&mut (dyn Write + Send)>> = reports
        .iter_mut()
        .map(|(_, writer)| Some(&mut **writer))
        .collect();
    with_run_pool(threads, search.program().name(), || -> Result<(), Error> {
        for (ordinal, group) in groups.iter().enumerate() {
            let command = search.command(inputs)?;
            // The warnings of the first search include the reader's and the queries' warnings
            // that the next searches repeat; the command line writes them once.
            let mut repeated = Vec::new();
            let mut sink = if ordinal == 0 {
                log.sink(&mut *diagnostics, true)
            } else {
                log.sink(&mut repeated, false)
            };
            let result = {
                let mut outputs = ReportOutputs {
                    formats: group
                        .iter()
                        .map(|&position| FormatOutput {
                            outfmt: formats[position],
                            sink: OutputSink::Writer(
                                writers[position].take().expect("one group per format"),
                            ),
                        })
                        .collect(),
                    diagnostics: &mut sink,
                    hits: None,
                    observer: None,
                };
                run_local(command, inputs, &mut outputs, &mut start)
            };
            result.map_err(|error| sink.failed(&error))?;
        }
        Ok(())
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
// ```c++
//     else {
//         m_OutputStream = &args[kArgOutput].AsOutputFile();
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:255-257
// ```c++
//         if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
//         	formatter.SetBaseFile(args[kArgOutput].AsString());
//         }
// ```
/// Runs the search and writes the report of one `-outfmt` to the file `path`, as the command
/// line's `-out path` does, and the command line's warnings to `diagnostics`. Outfmt 13 and 14
/// come through here: the file is the master file, which names the side files that the search
/// writes next to it (`path_1.json` or `path_1.xml`, ...). The error is the command line's, as in
/// [`run`].
pub fn run_to_file(
    search: &Search,
    inputs: Inputs<'_>,
    outfmt: &str,
    path: &Path,
    diagnostics: &mut (dyn Write + Send),
) -> Result<(), Error> {
    // `-` is the command line's standard output, which the library does not write.
    if path.as_os_str() == "-" {
        return Err(Error::Argv(
            "run_to_file writes a file: `run` gives a report to a writer".to_string(),
        ));
    }
    let mut start = CallStart::new(search, inputs, &mut *diagnostics)?;
    // The command line's dry run does not make the `-out` file (line RG item E1b).
    if start.dry_run() {
        return Ok(());
    }
    // tblastx writes the errors of an `-outfmt` value and of an `-out` it cannot make to its
    // `-logfile` file, as its other errors (`RunLog::failed`).
    let mut log = start.take_log();
    let threads = {
        let formats = [outfmt];
        searches(search, inputs, &formats, &mut log)?.1
    };
    let mut command = search.command(inputs)?;
    set_output_file(&mut command, path);
    let file =
        File::create(path).map_err(|_| log.failed(&crate::cli::inaccessible("out", path)))?;
    let mut file = BufWriter::new(file);
    with_run_pool(threads, search.program().name(), || -> Result<(), Error> {
        let mut sink = log.sink(&mut *diagnostics, true);
        let result = {
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(&mut file),
                }],
                diagnostics: &mut sink,
                hits: None,
                observer: None,
            };
            run_local(command, inputs, &mut outputs, &mut start)
        };
        result.map_err(|error| sink.failed(&error))
    })?;
    file.flush()
        .map_err(|error| Error::search(&anyhow::Error::from(error)))
}

/// The search's `-out` file: a report names its side files by it (`run_to_file`).
fn set_output_file(command: &mut Commands, path: &Path) {
    let out = Some(path.to_path_buf());
    match command {
        Commands::Blastp(args) => args.out = out,
        Commands::Tblastn(args) => args.out = out,
        Commands::Blastn(args) => args.out = out,
        Commands::Tblastx(args) => args.out = out,
        Commands::Blastx(args) => args.out = out,
        Commands::Psiblast(args) => args.out = out,
        Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {
            unreachable!("Program::parse admits the search programs only")
        }
    }
}

/// What one search of `run_collect` gathered.
struct Searched {
    /// Per format of the group: the report's bytes.
    writers: Vec<CountingWriter>,
    recorder: RangeRecorder,
    hits: Option<Arc<Vec<PairwiseHit>>>,
    /// The `-query_loc` interval.
    query_range: Option<SequenceRange>,
    ordinals: Option<Vec<usize>>,
}

/// One search of `run_collect`: the formats at `positions` of `formats`.
fn collect_one(
    search: &Search,
    inputs: Inputs<'_>,
    formats: &[&str],
    positions: &[usize],
    start: &mut CallStart,
    mut sink: RunSink<'_>,
) -> Result<Searched, Error> {
    let command = search.command(inputs)?;
    let query_range = query_range_of(&command);
    let ordinals = searched_queries(query_range, inputs.query);
    let mut writers: Vec<CountingWriter> =
        positions.iter().map(|_| CountingWriter::new()).collect();
    let mut recorder = RangeRecorder::new(writers.iter().map(|w| w.position.clone()).collect());
    let mut collected: Option<Arc<Vec<PairwiseHit>>> = None;
    let result = {
        let mut on_hits = |hits| collected = Some(hits);
        let mut outputs = ReportOutputs {
            formats: positions
                .iter()
                .zip(writers.iter_mut())
                .map(|(&position, writer)| FormatOutput {
                    outfmt: formats[position],
                    sink: OutputSink::Writer(writer),
                })
                .collect(),
            diagnostics: &mut sink,
            hits: Some(&mut on_hits),
            observer: Some(&mut recorder),
        };
        run_local(command, inputs, &mut outputs, start)
    };
    result.map_err(|error| sink.failed(&error))?;
    if let Some((hsp, _)) = recorder.open.iter().flatten().next() {
        recorder.fault(format!("HSP {hsp} started but did not end"));
    }
    if let Some(fault) = recorder.fault.take() {
        return Err(Error::Search(format!(
            "internal error in the output ranges: {fault}"
        )));
    }
    Ok(Searched {
        writers,
        recorder,
        hits: collected,
        query_range,
        ordinals,
    })
}

/// The index in the inputs of the engine's query `q_idx` (`searched_queries`).
pub(super) fn query_ordinal(ordinals: Option<&[usize]>, q_idx: usize) -> usize {
    ordinals
        .and_then(|ordinals| ordinals.get(q_idx).copied())
        .unwrap_or(q_idx)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1071-1077
// ```c++
//     int query_frame = 1, subject_frame = 1;
//     if (kTranslated) {
//         if (x_IsFieldRequested(eQueryFrame) || x_IsFieldRequested(eFrames)) {
//             query_frame = CAlignFormatUtil::
//                 GetFrame (q_start - 1, ds.GetSeqStrand(kQueryRow),
//                           scope.GetBioseqHandle(align.GetSeq_id(0)));
//         }
// ```
// The `-outfmt` columns print a frame for every row (1 for an untranslated row of BLASTP and
// BLASTN, the row's strand for BLASTN's subject, 0 for the protein row of BLASTX and TBLASTN); a
// translation frame exists where the program translates, the rest is the strand.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:880-887
// ```c++
//     if (x_IsFieldRequested(eQueryLength) && m_QueryLength == 0) {
//         ...
//     	if(m_QueryRange.NotEmpty())
//     		m_QueryLength = m_QueryRange.GetLength();
//     	else
//     		m_QueryLength = sequence::GetLength(align.GetSeq_id(0), &scope);
// ```
/// The HSP records of the final hit list of a search of `inputs`. `qlen` is the record's length
/// or the length of the `-query_loc` interval as given; `slen` is the record's, or the database
/// sequence's (`database_subjects`, which the hit list's subject numbers index).
// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
// ```c++
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// ```
// `shown_id` and `subject_tabular_id` are the first two columns of the default report.
/// The id of a record, read once: many HSPs show it.
fn shared_id(
    ids: &mut [Option<Arc<str>>],
    index: usize,
    read: impl FnOnce() -> Vec<u8>,
) -> Arc<str> {
    ids[index]
        .get_or_insert_with(|| Arc::from(String::from_utf8_lossy(&read())))
        .clone()
}

pub(super) fn records(
    program: Program,
    inputs: Inputs<'_>,
    hits: &[PairwiseHit],
    query_range: Option<SequenceRange>,
    ordinals: Option<&[usize]>,
    database_subjects: Option<&[ReportedSubject]>,
) -> Vec<HspRecord> {
    let query_protein = is_protein(program, Role::Query);
    let subject_protein = is_protein(program, Role::Subject);
    let translated_query = matches!(program, Program::Blastx | Program::Tblastx);
    let translated_subject = matches!(program, Program::Tblastn | Program::Tblastx);
    let mut ranks: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut query_ids = vec![None; inputs.query.len()];
    let mut subject_ids = vec![None; inputs.subject.len()];
    hits.iter()
        .enumerate()
        .map(|(index, record)| {
            let hit = &record.hit;
            let q_idx = query_ordinal(ordinals, hit.q_idx as usize);
            let s_idx = hit.s_idx as usize;
            let length = |records: &[FastaRecord], index: usize| {
                records.get(index).map_or(0, |record| record.seq().len())
            };
            let rank = ranks.entry(q_idx).or_default();
            let this_rank = *rank;
            *rank += 1;
            // A database's sequences are the hit list's (`FormatObserver::database_subjects`).
            let (subject_index, subject_id, subject_length) = match database_subjects {
                Some(subjects) => {
                    let subject = &subjects[s_idx];
                    (subject.oid as usize, subject.id.clone(), subject.length)
                }
                None => (
                    s_idx,
                    shared_id(&mut subject_ids, s_idx, || {
                        let subject = &inputs.subject[s_idx];
                        subject
                            .subject_tabular_id(subject_protein)
                            .unwrap_or_else(|_| subject.shown_id().to_vec())
                    }),
                    length(inputs.subject, s_idx),
                ),
            };
            HspRecord {
                index,
                q_idx,
                s_idx: subject_index,
                rank: this_rank,
                query_id: shared_id(&mut query_ids, q_idx, || {
                    inputs.query[q_idx].shown_id().to_vec()
                }),
                subject_id,
                raw_score: hit.raw_score,
                bit_score: hit.bit_score,
                e_value: hit.e_value,
                q_start: hit.q_start,
                q_end: hit.q_end,
                s_start: hit.s_start,
                s_end: hit.s_end,
                query_strand: strand(query_protein, hit.q_start, hit.q_end),
                subject_strand: strand(subject_protein, hit.s_start, hit.s_end),
                query_frame: record
                    .query_frame
                    .filter(|_| translated_query)
                    .map(i32::from),
                subject_frame: record
                    .subject_frame
                    .filter(|_| translated_subject)
                    .map(i32::from),
                query_length: query_range
                    .map_or_else(|| length(inputs.query, q_idx), |range| range.length()),
                subject_length,
                length: hit.length,
                num_ident: hit.num_ident,
                // BLASTN's hits have no positives: its reports show the identities.
                num_positives: record.positives.unwrap_or(hit.num_ident),
                mismatch: hit.mismatch,
                gap_open: hit.gapopen,
                gaps: record.gaps.unwrap_or_else(|| hit.gap_letters()),
                sum_n: record.sum_n,
                comp_adjust_method: record.comp_adjust_method,
                query_aligned: record.query_seq.clone(),
                subject_aligned: record.subject_seq.clone(),
            }
        })
        .collect()
}

/// Runs the search and gathers the report of each requested `-outfmt` value, the HSP records,
/// the statistics of each query and the diagnostics.
///
/// Each report is what the command line writes to `-out` for that `-outfmt` alone: the formats
/// are grouped by hit list size (`search_groups`), one search per group, so a run with
/// `-num_descriptions` or `-num_alignments` may search twice. The HSP records and the statistics
/// are those of the first search (the tabular formats'), with their rows and sections in the
/// reports of that search; the reports of another search give them no span. The diagnostics
/// are the first search's: the command line's standard error for its first format (the first
/// tabular format when the formats split, the first requested format otherwise), which holds
/// the warnings the command line writes once for the run.
///
/// Any failure fails the call with the command line's message, for the format whose search
/// failed. A `-outfmt` value that the command line rejects, and a call without formats, fail
/// before any search (`Error::Argv`).
pub fn run_collect(
    search: &Search,
    inputs: Inputs<'_>,
    formats: &[&str],
) -> Result<Outcome, Error> {
    let mut diagnostics: Vec<u8> = Vec::new();
    let mut start = CallStart::new(search, inputs, &mut diagnostics)?;
    if start.dry_run() {
        return Ok(Outcome {
            reports: formats
                .iter()
                .map(|outfmt| Report {
                    outfmt: outfmt.to_string(),
                    bytes: Vec::new(),
                    spans: Vec::new(),
                    subject_headings: Vec::new(),
                })
                .collect(),
            hsps: Vec::new(),
            query_stats: vec![None; inputs.query.len()],
            diagnostics,
        });
    }
    let mut log = start.take_log();
    let (groups, threads) = searches(search, inputs, formats, &mut log)?;
    let mut searched: Vec<Searched> = Vec::new();
    // The searches of one call share one thread pool (`with_run_pool`), as the command line's
    // one search has one.
    with_run_pool(threads, search.program().name(), || -> Result<(), Error> {
        for (ordinal, group) in groups.iter().enumerate() {
            // The warnings of the first search include the reader's and the queries' warnings
            // that the next searches repeat; the command line writes them once (the next
            // searches' do not go to the `-logfile` file either; a failure does, as in `run`).
            let mut repeated = Vec::new();
            let sink = if ordinal == 0 {
                log.sink(&mut diagnostics, true)
            } else {
                log.sink(&mut repeated, false)
            };
            searched.push(collect_one(
                search, inputs, formats, group, &mut start, sink,
            )?);
        }
        Ok(())
    })?;
    let first = searched.first();
    let ordinals = first.and_then(|first| first.ordinals.as_deref());
    let hsps = match first {
        Some(Searched {
            hits: Some(hits),
            query_range,
            recorder,
            ..
        }) => records(
            search.program(),
            inputs,
            hits,
            *query_range,
            ordinals,
            recorder.database_subjects.as_deref(),
        ),
        _ => Vec::new(),
    };
    let mut query_stats = vec![None; inputs.query.len()];
    if let Some(first) = first {
        for (q_idx, stats) in first.recorder.query_stats.iter().enumerate() {
            if let Some(slot) = query_stats.get_mut(query_ordinal(ordinals, q_idx)) {
                *slot = *stats;
            }
        }
    }
    let mut reports: Vec<Option<Report>> = (0..formats.len()).map(|_| None).collect();
    for (ordinal, (group, done)) in groups.iter().zip(&mut searched).enumerate() {
        for (slot, (&position, writer)) in group.iter().zip(done.writers.drain(..)).enumerate() {
            // The places of the records' HSPs: those of the records' search, by HSP.
            let places = |recorded: &mut Vec<Vec<Option<Range<usize>>>>| {
                let mut places = match ordinal {
                    0 => std::mem::take(&mut recorded[slot]),
                    _ => Vec::new(),
                };
                places.resize(hsps.len(), None);
                places
            };
            reports[position] = Some(Report {
                outfmt: formats[position].to_string(),
                bytes: writer.bytes,
                spans: places(&mut done.recorder.spans),
                subject_headings: places(&mut done.recorder.headings),
            });
        }
    }
    Ok(Outcome {
        reports: reports
            .into_iter()
            .map(|report| report.expect("every requested format belongs to one group"))
            .collect(),
        hsps,
        query_stats,
        diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::{validate, Program};

    const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

    fn one_protein() -> (Vec<FastaRecord>, Vec<FastaRecord>) {
        (
            vec![FastaRecord::new("Query_1", b"q", PROTEIN)],
            vec![FastaRecord::new("Subject_1", b"s", PROTEIN)],
        )
    }

    #[test]
    fn the_strand_is_ncbis_rule_and_the_ranges_are_half_open_on_the_plus_strand() {
        assert_eq!(strand(true, 9, 1), None);
        assert_eq!(strand(false, 1, 9), Some(Strand::Plus));
        assert_eq!(strand(false, 9, 1), Some(Strand::Minus));
        // A row of one letter: `(s_start - s_end) > 0` is false.
        assert_eq!(strand(false, 5, 5), Some(Strand::Plus));
        let (query, subject) = one_protein();
        let search = validate(Program::Blastp, &[]).unwrap();
        let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
        let mut hsp = outcome.hsps[0].clone();
        (hsp.q_start, hsp.q_end, hsp.s_start, hsp.s_end) = (3, 9, 9, 3);
        assert_eq!(hsp.query_range(), 2..9);
        assert_eq!(hsp.subject_range(), 2..9);
        (hsp.s_start, hsp.s_end) = (5, 5);
        assert_eq!(hsp.subject_range(), 4..5);
    }

    #[test]
    fn groups_follow_the_hit_list_sizes() {
        let inputs = Inputs::new(&[], &[]);
        let plain = validate(Program::Blastp, &[])
            .unwrap()
            .command(inputs)
            .unwrap();
        assert_eq!(search_groups(&plain, &[6, 0, 7]), vec![vec![0, 1, 2]]);
        // -num_alignments changes the pairwise formats' hit list size only: two searches, the
        // tabular formats first.
        let split = validate(Program::Blastp, &["-num_alignments", "1"])
            .unwrap()
            .command(inputs)
            .unwrap();
        assert_eq!(search_groups(&split, &[0, 6, 7]), vec![vec![1, 2], vec![0]]);
        // -max_target_seqs sets every format's size.
        let same = validate(Program::Blastp, &["-max_target_seqs", "5"])
            .unwrap()
            .command(inputs)
            .unwrap();
        assert_eq!(search_groups(&same, &[0, 6]), vec![vec![0, 1]]);
    }

    #[test]
    fn a_self_search_gives_one_hsp_with_its_rows_and_sections() {
        let (query, subject) = one_protein();
        let search = validate(Program::Blastp, &[]).unwrap();
        let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6", "0"]).unwrap();
        assert_eq!(outcome.reports.len(), 2);
        assert_eq!(outcome.reports[0].outfmt, "6");
        assert_eq!(outcome.reports[1].outfmt, "0");
        assert!(
            outcome.diagnostics.is_empty(),
            "{:?}",
            String::from_utf8_lossy(&outcome.diagnostics)
        );
        assert_eq!(outcome.hsps.len(), 1);
        let hsp = &outcome.hsps[0];
        assert_eq!((hsp.index, hsp.q_idx, hsp.s_idx, hsp.rank), (0, 0, 0, 0));
        assert_eq!(
            (hsp.q_start, hsp.q_end, hsp.s_start, hsp.s_end),
            (1, 66, 1, 66)
        );
        assert!(hsp.bit_score > 100.0);
        assert!(hsp.e_value < 1e-30);
        assert_eq!(
            hsp.query_aligned.as_deref(),
            Some(std::str::from_utf8(PROTEIN).unwrap())
        );
        // The row in outfmt 6 is the report's one line; the section in outfmt 0 starts at
        // the score line and the heading at the defline.
        let row = outcome.reports[0].spans[0].clone().unwrap();
        let tabular = &outcome.reports[0].bytes;
        assert_eq!(&tabular[row.clone()], tabular.as_slice());
        assert!(tabular.ends_with(b"\n"));
        assert_eq!(tabular.iter().filter(|&&b| b == b'\n').count(), 1);
        let section = outcome.reports[1].spans[0].clone().unwrap();
        let pairwise = &outcome.reports[1].bytes;
        assert!(pairwise[section].starts_with(b" Score ="));
        let heading = outcome.reports[1].subject_headings[0].clone().unwrap();
        assert!(pairwise[heading].starts_with(b">"));
        assert!(outcome.reports[0].subject_headings[0].is_none());
    }

    #[test]
    fn two_searches_give_each_format_its_own_report_and_records_from_the_first() {
        let (query, subject) = one_protein();
        let search = validate(Program::Blastp, &["-num_alignments", "1"]).unwrap();
        let outcome = run_collect(&search, Inputs::new(&query, &subject), &["0", "6"]).unwrap();
        assert_eq!(outcome.reports[0].outfmt, "0");
        assert_eq!(outcome.reports[1].outfmt, "6");
        assert_eq!(outcome.hsps.len(), 1);
        // The records come from the tabular search: a row in report 1, no section in report
        // 0 (another search).
        assert!(outcome.reports[1].spans[0].is_some());
        assert_eq!(outcome.reports[0].spans, [None]);
        // The same bytes as the formats searched alone.
        let alone = run_collect(&search, Inputs::new(&query, &subject), &["0"]).unwrap();
        assert_eq!(alone.reports[0].bytes, outcome.reports[0].bytes);
    }

    #[test]
    fn the_run_writes_to_the_callers_outputs_and_the_errors_are_the_command_lines() {
        let (query, subject) = one_protein();
        let search = validate(Program::Blastp, &[]).unwrap();
        let mut report = Vec::new();
        let mut diagnostics = Vec::new();
        run(
            &search,
            Inputs::new(&query, &subject),
            &mut [("6", &mut report)],
            &mut diagnostics,
        )
        .unwrap();
        let collected = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
        assert_eq!(report, collected.reports[0].bytes);
        assert!(diagnostics.is_empty());
        // A format the command line rejects fails before any search, as does a call without one.
        let error = run_collect(&search, Inputs::new(&query, &subject), &["99"]).unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{error}");
        let error = run_collect(&search, Inputs::new(&query, &subject), &[]).unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{error}");
        // An engine error is the command line's message.
        let empty: Vec<FastaRecord> = Vec::new();
        let error = run_collect(&search, Inputs::new(&query, &empty), &["6"]).unwrap_err();
        assert!(matches!(error, Error::Search(_)), "{error}");
        assert!(
            error.message().contains("Empty CBlastQueryVector"),
            "{error}"
        );
    }
}
