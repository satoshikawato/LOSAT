//! Running a search from memory: `run` writes to the caller's outputs, `run_collect` gathers
//! the reports, the HSP records and the diagnostics.
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

use std::io::{self, Write};
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::error::Error;
use super::input::Inputs;
use super::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    with_run_pool, FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs,
};
use super::search::Search;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::cli::Commands;
use crate::report::PairwiseHit;

/// One report of a search: the bytes that the command line writes to `-out` for this
/// `-outfmt`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Report {
    /// The `-outfmt` value, as requested.
    pub outfmt: String,
    /// The report's bytes.
    pub bytes: Vec<u8>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
// ```c
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
// ```
// Every output format prints from the same final result set, so a position in that result
// identifies an HSP in every format (`HspIndex`).
/// One HSP of the final, sorted hit list of a search, with the engine's values and its place in
/// each report. The fields are those of LOSAT Web ABI v2's HSP record (docs/web/abi_v2.md §8).
/// Displayed numbers (identities, the E-value as printed) are read from the HSP's row of a
/// tabular report (`spans`), not formatted again from the values here.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct HspRecord {
    /// 0-based position of the HSP in the final, sorted hit list; the same HSP has the same
    /// index in every report.
    pub index: usize,
    /// 0-based index of the query record in `Inputs::query` (records without residues count;
    /// with `-query_loc`, a query record that NCBI skips because its interval starts past its
    /// end keeps its place).
    pub q_idx: usize,
    /// 0-based index of the subject record in `Inputs::subject`.
    pub s_idx: usize,
    /// 0-based position of the HSP among the HSPs of its query.
    pub rank: usize,
    /// The raw alignment score.
    pub raw_score: i32,
    /// The bit score, not rounded.
    pub bit_score: f64,
    /// The E-value, not rounded.
    pub e_value: f64,
    /// 1-based query start, as printed in outfmt 6.
    pub q_start: usize,
    /// 1-based query end, as printed in outfmt 6.
    pub q_end: usize,
    /// 1-based subject start, as printed in outfmt 6 (start > end: minus strand).
    pub s_start: usize,
    /// 1-based subject end, as printed in outfmt 6.
    pub s_end: usize,
    /// The query's translation frame, for the translated programs.
    pub query_frame: Option<i32>,
    /// The subject's translation frame, for the translated programs.
    pub subject_frame: Option<i32>,
    /// The subject's length.
    pub subject_length: Option<usize>,
    /// The aligned query row, `-` for gaps, where the search renders it.
    pub query_aligned: Option<String>,
    /// The aligned subject row, `-` for gaps, where the search renders it.
    pub subject_aligned: Option<String>,
    /// For each report of the `Outcome`, in its order: the byte range of this HSP in the
    /// report's bytes, the row of a tabular format or the section of the pairwise format (its
    /// score lines and alignment, without the subject's heading). `None` where the report does
    /// not show the HSP, or where the report came from another search than the records
    /// (`run_collect`).
    pub spans: Vec<Option<Range<usize>>>,
    /// For each report, in its order: the byte range of the heading of this HSP's subject (the
    /// defline and `Length=` line that the pairwise format writes before the first HSP of each
    /// subject); every HSP of the subject in that query has the same range. `None` for the
    /// formats without headings.
    pub subject_headings: Vec<Option<Range<usize>>>,
}

/// What `run_collect` gathers from a search.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Outcome {
    /// One report per requested format, in the requested order.
    pub reports: Vec<Report>,
    /// The HSPs of the final hit list, with their places in the reports.
    pub hsps: Vec<HspRecord>,
    /// What the command line writes to its standard error for the search: the warnings of the
    /// reader and of the search, each once. Empty for a search without warnings.
    pub diagnostics: Vec<u8>,
}

/// A writer into memory that tells its position, for the observer of its report.
struct CountingWriter {
    bytes: Vec<u8>,
    position: Arc<AtomicUsize>,
}

impl CountingWriter {
    fn new() -> Self {
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
/// The index in `Inputs::query` of each query that the engine searches: NCBI's batch reader
/// skips a query whose `-query_loc` interval starts more than one letter past its end
/// (`seq_range::cut_queries`), so the engine's `q_idx` counts the searched queries, and the
/// records give the index in the inputs. `None` without `-query_loc` (every query is searched).
fn searched_queries(
    query_loc: Option<&str>,
    program: &str,
    records: &[FastaRecord],
) -> Option<Vec<usize>> {
    use crate::blastinput::seq_range::{
        parse_optional_range, record_interval, RangeRole, RecordInterval,
    };
    let range = parse_optional_range(query_loc, RangeRole::Query, program).ok()??;
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

/// The `-query_loc` of the program's arguments and the program's name in NCBI's messages.
fn query_loc(command: &Commands) -> (Option<&str>, &'static str) {
    match command {
        Commands::Blastp(a) => (a.query_loc.as_deref(), "BLASTP"),
        Commands::Tblastn(a) => (a.query_loc.as_deref(), "TBLASTN"),
        Commands::Blastn(a) => (a.query_loc.as_deref(), "BLASTN"),
        Commands::Tblastx(a) => (a.query_loc.as_deref(), "TBLASTX"),
        Commands::Blastx(a) => (a.query_loc.as_deref(), "BLASTX"),
        Commands::Psiblast(_) | Commands::Makedb(_) => (None, ""),
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
/// `-num_threads` of the program's arguments.
fn num_threads(command: &Commands) -> usize {
    match command {
        Commands::Blastp(a) => a.num_threads,
        Commands::Tblastn(a) => a.num_threads,
        Commands::Blastn(a) => a.num_threads,
        Commands::Tblastx(a) => a.num_threads,
        // A positive `i32` (the parser's `positive`).
        Commands::Blastx(a) => a.num_threads as usize,
        Commands::Psiblast(_) | Commands::Makedb(_) => 1,
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
        Commands::Psiblast(_) | Commands::Makedb(_) => [None; 3],
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
pub fn search_groups(command: &Commands, numbers: &[i32]) -> Vec<Vec<usize>> {
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
/// One `run_local` call: the program's entry for the parsed arguments.
fn run_local(
    command: Commands,
    query: &[FastaRecord],
    subject: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> anyhow::Result<()> {
    match command {
        Commands::Blastp(args) => run_local_blastp(args, query, subject, "", "", outputs),
        Commands::Tblastn(args) => run_local_tblastn(args, query, subject, outputs),
        Commands::Blastn(args) => run_local_blastn(args, query, subject, outputs),
        Commands::Tblastx(args) => run_local_tblastx(args, query, subject, outputs),
        Commands::Blastx(args) => run_local_blastx(args, query, subject, outputs),
        Commands::Psiblast(_) | Commands::Makedb(_) => {
            unreachable!("Program::parse admits the five search programs only")
        }
    }
}

/// Runs the search once and writes to the caller's outputs: every format of `outputs` from the
/// one search, the diagnostics, the hit list and the observer events as `ReportOutputs`
/// describes them. The formats share one search, so they must have one hit list size: formats
/// whose sizes differ (`search_groups`) are separate `run` calls, or one `run_collect`.
///
/// The error is the command line's: the message of its parser or of the engine
/// (`Error::Search`). The search writes nothing to the process's standard output or error, and
/// never ends the process.
pub fn run(
    search: &Search,
    inputs: Inputs<'_>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<(), Error> {
    let command = search.command(inputs.query_name, inputs.subject_name);
    run_local(command, inputs.query, inputs.subject, outputs).map_err(|error| Error::search(&error))
}

/// What one search of `run_collect` gathered.
struct Searched {
    /// Per format of the group: the report's bytes.
    writers: Vec<CountingWriter>,
    recorder: RangeRecorder,
    hits: Option<Arc<Vec<PairwiseHit>>>,
    ordinals: Option<Vec<usize>>,
}

/// One search of `run_collect`: the formats at `positions` of `formats`.
fn collect_one(
    search: &Search,
    inputs: Inputs<'_>,
    formats: &[&str],
    positions: &[usize],
    diagnostics: &mut (dyn Write + Send),
) -> Result<Searched, Error> {
    let command = search.command(inputs.query_name, inputs.subject_name);
    let (query_loc, program) = query_loc(&command);
    let ordinals = searched_queries(query_loc, program, inputs.query);
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
            diagnostics,
            hits: Some(&mut on_hits),
            observer: Some(&mut recorder),
        };
        run_local(command, inputs.query, inputs.subject, &mut outputs)
    };
    result.map_err(|error| Error::search(&error))?;
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
        ordinals,
    })
}

/// The HSP records of the hit list of the first search, with their places in the reports
/// (`None` for the reports of the other searches).
fn records(
    hits: &[PairwiseHit],
    recorder: &RangeRecorder,
    positions: &[usize],
    ordinals: Option<&[usize]>,
    report_count: usize,
) -> Vec<HspRecord> {
    let mut ranks: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    hits.iter()
        .enumerate()
        .map(|(index, record)| {
            let hit = &record.hit;
            let q_idx = ordinals
                .and_then(|ordinals| ordinals.get(hit.q_idx as usize).copied())
                .unwrap_or(hit.q_idx as usize);
            let rank = ranks.entry(q_idx).or_default();
            let this_rank = *rank;
            *rank += 1;
            let mut spans = vec![None; report_count];
            let mut subject_headings = vec![None; report_count];
            for (slot, &position) in positions.iter().enumerate() {
                spans[position] = recorder.spans[slot].get(index).cloned().flatten();
                subject_headings[position] = recorder.headings[slot].get(index).cloned().flatten();
            }
            HspRecord {
                index,
                q_idx,
                s_idx: hit.s_idx as usize,
                rank: this_rank,
                raw_score: hit.raw_score,
                bit_score: hit.bit_score,
                e_value: hit.e_value,
                q_start: hit.q_start,
                q_end: hit.q_end,
                s_start: hit.s_start,
                s_end: hit.s_end,
                query_frame: record.query_frame.map(i32::from),
                subject_frame: record.subject_frame.map(i32::from),
                subject_length: record.subject_length,
                query_aligned: record.query_seq.clone(),
                subject_aligned: record.subject_seq.clone(),
                spans,
                subject_headings,
            }
        })
        .collect()
}

/// Runs the search and gathers the report of each requested `-outfmt` value, the HSP records
/// and the diagnostics.
///
/// Each report is what the command line writes to `-out` for that `-outfmt` alone: the formats
/// are grouped by hit list size (`search_groups`), one search per group, so a run with
/// `-num_descriptions` or `-num_alignments` may search twice. The HSP records are the final hit
/// list of the first search (the tabular formats'), with their rows and sections in the
/// reports of that search; the reports of another search give them no span. The diagnostics
/// are the first search's: the command line's standard error for its first format (the first
/// tabular format when the formats split, the first requested format otherwise), which holds
/// the warnings the command line writes once for the run.
///
/// Any failure fails the call with the command line's message, for the format whose search
/// failed. A `-outfmt` value that the command line rejects fails before any search
/// (`Error::Argv`).
pub fn run_collect(
    search: &Search,
    inputs: Inputs<'_>,
    formats: &[&str],
) -> Result<Outcome, Error> {
    use crate::blastinput::app::parse_formatting_string;
    let numbers: Vec<i32> = formats
        .iter()
        .map(|outfmt| {
            parse_formatting_string(outfmt)
                .map(|choice| choice.number)
                .map_err(|error| Error::Argv(format!("{error:#}")))
        })
        .collect::<Result<_, _>>()?;
    let command = search.command(inputs.query_name, inputs.subject_name);
    let groups = search_groups(&command, &numbers);
    let threads = num_threads(&command);
    drop(command);
    let mut diagnostics: Vec<u8> = Vec::new();
    let mut searched: Vec<Searched> = Vec::new();
    // The searches of one call share one thread pool (`with_run_pool`), as the command line's
    // one search has one.
    with_run_pool(threads, search.program().name(), || -> Result<(), Error> {
        for (ordinal, group) in groups.iter().enumerate() {
            // The warnings of the first search include the reader's and the queries' warnings
            // that the next searches repeat; the command line writes them once.
            let mut scratch = Vec::new();
            let sink: &mut (dyn Write + Send) = if ordinal == 0 {
                &mut diagnostics
            } else {
                &mut scratch
            };
            searched.push(collect_one(search, inputs, formats, group, sink)?);
        }
        Ok(())
    })?;
    let mut reports: Vec<Option<Report>> = (0..formats.len()).map(|_| None).collect();
    for (group, done) in groups.iter().zip(&mut searched) {
        for (&position, writer) in group.iter().zip(done.writers.drain(..)) {
            reports[position] = Some(Report {
                outfmt: formats[position].to_string(),
                bytes: writer.bytes,
            });
        }
    }
    let hsps = match (groups.first(), searched.first()) {
        (Some(group), Some(first)) => match &first.hits {
            Some(hits) => records(
                hits,
                &first.recorder,
                group,
                first.ordinals.as_deref(),
                formats.len(),
            ),
            None => Vec::new(),
        },
        _ => Vec::new(),
    };
    Ok(Outcome {
        reports: reports
            .into_iter()
            .map(|report| report.expect("every requested format belongs to one group"))
            .collect(),
        hsps,
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
    fn groups_follow_the_hit_list_sizes() {
        let plain = validate(Program::Blastp, &[]).unwrap().command("q", "s");
        assert_eq!(search_groups(&plain, &[6, 0, 7]), vec![vec![0, 1, 2]]);
        // -num_alignments changes the pairwise formats' hit list size only: two searches, the
        // tabular formats first.
        let split = validate(Program::Blastp, &["-num_alignments", "1"])
            .unwrap()
            .command("q", "s");
        assert_eq!(search_groups(&split, &[0, 6, 7]), vec![vec![1, 2], vec![0]]);
        // -max_target_seqs sets every format's size.
        let same = validate(Program::Blastp, &["-max_target_seqs", "5"])
            .unwrap()
            .command("q", "s");
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
        let row = hsp.spans[0].clone().unwrap();
        let tabular = &outcome.reports[0].bytes;
        assert_eq!(&tabular[row.clone()], tabular.as_slice());
        assert!(tabular.ends_with(b"\n"));
        assert_eq!(tabular.iter().filter(|&&b| b == b'\n').count(), 1);
        let section = hsp.spans[1].clone().unwrap();
        let pairwise = &outcome.reports[1].bytes;
        assert!(pairwise[section].starts_with(b" Score ="));
        let heading = hsp.subject_headings[1].clone().unwrap();
        assert!(pairwise[heading].starts_with(b">"));
        assert!(hsp.subject_headings[0].is_none());
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
        assert!(outcome.hsps[0].spans[1].is_some());
        assert!(outcome.hsps[0].spans[0].is_none());
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
        {
            let mut outputs =
                ReportOutputs::single("6", OutputSink::Writer(&mut report), &mut diagnostics);
            run(&search, Inputs::new(&query, &subject), &mut outputs).unwrap();
        }
        let collected = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
        assert_eq!(report, collected.reports[0].bytes);
        assert!(diagnostics.is_empty());
        // A format the command line rejects fails before any search.
        let error = run_collect(&search, Inputs::new(&query, &subject), &["99"]).unwrap_err();
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
