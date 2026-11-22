//! Running a search from memory: `run` writes each report to the caller's writer, `run_collect`
//! gathers the reports, the HSP records, the statistics of each query and the diagnostics.
//!
//! NCBI's application searches once and prints the result set with one formatter per
//! requested output; LOSAT's `run_local` entries do the same for records in memory, and this
//! module hands them the caller's inputs and outputs.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:195-295
#![warn(missing_docs, unused)]

use std::io::{self, Write};
use std::ops::Range;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::error::Error;
use super::id::Id;
use super::input::{check_archive_inputs, is_protein, Inputs, Role};
use super::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_database_blastn,
    run_local_database_blastp, run_local_database_blastx, run_local_database_tblastn,
    run_local_database_tblastx, run_local_tblastn, run_local_tblastx, with_run_pool,
    FormatObserver, FormatOutput, HspIndex, OutFile, OutputSink, PostsFirst, ReportOutputs,
    ReportedSubject,
};
use super::parts::Parts;
use super::program::Program;
use super::search::Search;
use super::start::{CallStart, RunLog, RunSink};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::{
    parse_optional_range, record_interval, RangeRole, RecordInterval, SequenceRange,
};
use crate::cli::Commands;
use crate::report::archive::options::RequestKind;
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
/// The strand that NCBI's `sstrand` column gives a row of a nucleotide sequence, and no strand
/// for a protein: an end before the start is the minus strand.
fn strand(protein: bool, start: usize, end: usize) -> Option<Strand> {
    (!protein).then_some(if start > end {
        Strand::Minus
    } else {
        Strand::Plus
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
// Every output format prints from the same final result set, so a position in that result
// identifies an HSP in every format (`HspIndex`).
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-132
/// One HSP of the final, sorted hit list of a search, with the engine's values. Each field that
/// names a column of `-outfmt 6` holds the number that the column prints for the HSP's row, not
/// rounded (`tests/library_api.rs` holds them to the columns); the text of a report is its
/// tabular row (`Report::spans`), not the values here formatted again.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct HspRecord {
    /// 0-based position of the HSP in the final, sorted hit list; the same HSP has the same
    /// index in every report (`Report::spans`). With psiblast's rounds the hit lists of the
    /// rounds of each query follow one another, as the reports write them.
    pub index: usize,
    /// 0-based index of the query record in the inputs (records without residues count; with
    /// `-query_loc`, a query record that NCBI skips because its interval starts past its end
    /// keeps its place).
    pub q_idx: usize,
    /// 0-based index of the subject record in the inputs; for a database
    /// (`Inputs::database`), the sequence's OID (the `n` of `gnl|BL_ORD_ID|n`).
    pub s_idx: usize,
    /// 0-based position of the HSP among the HSPs of its query in its round (`round`).
    pub rank: usize,
    /// psiblast's round that found the HSP (1-based, NCBI's `Results from round N` and outfmt 7's
    /// `# Iteration: N`); 1 in a search without rounds. The HSPs of each round of a query follow
    /// those of the round before, each round's in the order of its own hit list.
    pub round: u32,
    /// The query's id (`qaccver`): the first word of its title, its generated `Query_N` without a
    /// title, and with `-parse_deflines` the accession and version of its id; the bytes of the
    /// reports, also those that are not UTF-8 ([`Id`]). Cloning is cheap; it does not copy the
    /// bytes.
    pub query_id: Id,
    /// The subject's id (`saccver`), by the same rules; for a database, the `saccver` of the
    /// sequence. Where NCBI writes no `saccver` for the subject (its lookup fails: a title whose
    /// encoding it cannot tell, or `-parse_deflines` ids that LOSAT refuses), the id by the
    /// query's rules; for a database sequence, the first word of its title. Cloning is cheap; it
    /// does not copy the bytes.
    pub subject_id: Id,
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
    /// How many HSPs NCBI's sum statistics linked with this one (the `n` of `Expect` in the
    /// pairwise format), where the program links HSPs; `None` otherwise. 0 and 1 both mean an HSP
    /// on its own (`BlastHSP::num`).
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
    /// The percentage of identical columns (`pident`), not rounded.
    pub fn percent_identity(&self) -> f64 {
        crate::report::outfmt6::percent_identity_ncbi(self.num_ident, self.length)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1109-1120
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_iteration.cpp:76-91
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_results.hpp:211-215
// Each round of psiblast is a result set of its own, with its own ancillary data: a PSSM round's
// Karlin-Altschul values are the PSSM's.
/// The rounds of one query of a search: psiblast's rounds (`-num_iterations` other than 1),
/// one round for the other searches.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct QueryRounds {
    /// The statistics of each round of the query, round 1 first: what fixed the E-values of the
    /// HSPs of that round (`HspRecord::round`). The first is `Outcome::query_stats`'s; a search
    /// without rounds, and a query that was not searched, have that one only.
    pub stats: Vec<Option<QueryStats>>,
    /// Whether the query's rounds converged (NCBI's `Search has CONVERGED!`, `HasConverged`):
    /// after two rounds or more, every subject below `-inclusion_ethresh` in the last round was
    /// below it in the round before too (`algorithm::psiblast::rounds::IterationState`). False
    /// without rounds.
    pub converged: bool,
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
    /// its batch). With psiblast's rounds, round 1's (`rounds`).
    pub query_stats: Vec<Option<QueryStats>>,
    /// One entry per query record of the inputs, in their order: the statistics of each of its
    /// rounds and whether they converged.
    pub rounds: Vec<QueryRounds>,
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3613-3632
/// Records, for every report of one search, the byte range of each HSP's row or section and
/// of the heading of each HSP's subject, from the formatter's observer events; and the parts of
/// the search's results (`Parts`: with psiblast's rounds, each round's call numbers its HSPs
/// from 0).
struct RangeRecorder {
    positions: Vec<Arc<AtomicUsize>>,
    /// Per report: the HSP whose row or section is being written, and where it started.
    open: Vec<Option<(HspIndex, usize)>>,
    /// Per part, per report, per HSP of the part: the row or section.
    spans: Vec<Vec<Vec<Option<Range<usize>>>>>,
    /// Per part, per report, per HSP of the part: the heading of the HSP's subject.
    headings: Vec<Vec<Vec<Option<Range<usize>>>>>,
    /// Per report: the start of the heading being written, and the last complete heading.
    heading_start: Vec<Option<usize>>,
    heading: Vec<Option<Range<usize>>>,
    /// The statistics, database ranks and rounds of each part.
    parts: Parts,
    /// The first observer event that does not fit the formatter contract; the search fails.
    fault: Option<String>,
}

impl RangeRecorder {
    fn new(positions: Vec<Arc<AtomicUsize>>) -> Self {
        let reports = positions.len();
        Self {
            positions,
            open: vec![None; reports],
            spans: vec![vec![Vec::new(); reports]],
            headings: vec![vec![Vec::new(); reports]],
            heading_start: vec![None; reports],
            heading: vec![None; reports],
            parts: Parts::new(),
            fault: None,
        }
    }

    /// The places of part `part`'s HSPs in report `slot`, `length` of them.
    fn places(
        recorded: &mut [Vec<Vec<Option<Range<usize>>>>],
        part: usize,
        slot: usize,
        length: usize,
    ) -> Vec<Option<Range<usize>>> {
        let mut places = recorded
            .get_mut(part)
            .map(|reports| std::mem::take(&mut reports[slot]))
            .unwrap_or_default();
        places.resize(length, None);
        places
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
        let part = self.parts.current();
        let stored =
            store(&mut self.spans[part][format], hsp, start..end).and_then(|()| {
                match self.heading[format].clone() {
                    Some(heading) => store(&mut self.headings[part][format], hsp, heading),
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
        self.parts.query_stats(q_idx, stats);
    }

    fn database_subjects(&mut self, subjects: Vec<ReportedSubject>) {
        self.parts.database_subjects(subjects);
    }

    fn round_begin(&mut self, query: usize, round: u32) {
        self.parts.round_begin(query, round);
        // A round's report starts its subjects' headings again.
        self.heading_start.fill(None);
        self.heading.fill(None);
        let reports = self.positions.len();
        self.spans
            .resize(self.parts.current() + 1, vec![Vec::new(); reports]);
        self.headings
            .resize(self.parts.current() + 1, vec![Vec::new(); reports]);
    }

    fn rounds_converged(&mut self, query: usize) {
        self.parts.rounds_converged(query);
    }

    fn pssm_query(&mut self, query: &FastaRecord) {
        self.parts.pssm_query(query);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973
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

/// The `-outfmt` values of `formats` that the search of the positions `group` (a group of
/// `search_groups`) does not write, in the call's order: those of the call's other searches, over
/// which, with its own, a search decides the call's refusals (`ReportOutputs::other_outfmts`).
/// Empty where the call is one search.
fn other_formats<'f>(formats: &[&'f str], group: &[usize]) -> Vec<&'f str> {
    formats
        .iter()
        .enumerate()
        .filter(|(position, _)| !group.contains(position))
        .map(|(_, outfmt)| *outfmt)
        .collect()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-220
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-133
// The formatters of the plural names (`x_CheckTaxDB`) and the organism report (outfmt 18) look
// the taxonomy database up themselves, through the search path of NCBI's databases: a search of
// records reads the application settings for them, as the command line does, and the search asks
// where the command line asks (`common::database::subject_taxonomy`).
/// The `[BLAST] BLASTDB` setting of the call (`start`) when a report of `outputs` (as `app`'s
/// formatter writes it, `TabularSpec::of_app`) looks every taxid up (the plural names of the
/// taxonomy database, the organism report); `None` otherwise.
fn taxdb_setting(
    outputs: &ReportOutputs<'_>,
    start: &CallStart,
    app: RequestKind,
) -> Option<crate::algorithm::common::database::TaxdbSetting> {
    use crate::blastinput::app::{parse_formatting_string, ReportFormat};
    let looks_up = outputs.formats.iter().any(|format| {
        parse_formatting_string(format.outfmt).is_ok_and(|choice| {
            ReportFormat::from_number(choice.number)
                .and_then(|report| {
                    crate::report::tabular::TabularSpec::of_app(report, &choice, app)
                })
                .is_some_and(|spec| spec.looks_up_every_taxid())
        })
    });
    looks_up.then(|| crate::algorithm::common::database::TaxdbSetting::of(&start.settings()))
}

// NCBI's application layer reads its environment and registry files before its own code runs,
// and a database is found by them (the `BLASTDB` search path); LOSAT rejects the settings that
// change the output. A search of records reads them where NCBI's library looks something up
// itself: the directory of BLASTN's `-window_masker_taxid`, the taxonomy database of the plural
// names and the organism report (outfmt 18) and the file of a `-matrix` that NCBI's tables do not
// have (TBLASTX, ungapped BLASTX). The settings are the call's application start's
// (`start::CallStart::settings`).
/// One `run_local` call: the program's entry for the parsed arguments and `inputs`
/// (`run_local_entry`), whose end is the command line's run's: a run that NCBI ends with exit
/// status 0 and no word ends the call without an error, and an NCBI error after NCBI's main-thread
/// diagnostic stream failed is the text that the command line writes for it.
pub(super) fn run_local(
    command: Commands,
    inputs: Inputs<'_>,
    outputs: &mut ReportOutputs<'_>,
    start: &mut CallStart,
) -> anyhow::Result<()> {
    // The call's code runs on this thread with the call's own diagnostic stream (one batch of
    // `run_batches` after another's of another call, too), and its searches' workers with it
    // (`utils::threading::with_search_pool`).
    let _run = start.enter();
    let result = run_local_entry(command, inputs, outputs, start);
    // The `-out` file of `run_to_file` ends as the command line's run ends its `-out` stream
    // (`ReportOutputs::finish_out`), before the run's error is written as `CATCH_ALL` writes it.
    let Err(error) = outputs.finish_out(result) else {
        return Ok(());
    };
    let Some(native) = error.downcast_ref::<crate::cli::NativeError>() else {
        return Err(error);
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    // `CATCH_ALL`'s silent end (`cli::NativeError::silent_end`: the `-in_pssm` reader's
    // `CEofException`, the object manager's failed write after its diagnostic stream failed): NCBI's
    // run returns 0 without a word, so the call succeeds with what it wrote (tblastn's error
    // archive of the saved messages, `report::archive::silent_end_lines`, or nothing).
    if crate::cli::is_silent_end(&error) {
        return Ok(());
    }
    // After NCBI's object manager failed the stream (`ncbi_exception::main_diag_failbit`, this
    // call's run, `CallStart::new`), `CATCH_ALL` writes its line feeds alone
    // (`cli::catch_all_text`), as the command line writes the error.
    if crate::ncbi_exception::main_diag_failbit() {
        let message = crate::cli::catch_all_text(&native.message).to_vec();
        return Err(crate::cli::NativeError {
            exit: native.exit,
            message,
        }
        .into());
    }
    Err(error)
}

/// The program's entry of a `run_local` call.
fn run_local_entry(
    command: Commands,
    inputs: Inputs<'_>,
    outputs: &mut ReportOutputs<'_>,
    start: &mut CallStart,
) -> anyhow::Result<()> {
    let (query, subject) = (inputs.query, inputs.subject);
    // How the reading of the query records ended, which the search meets in its batch
    // (`Inputs::query_end`).
    let end = inputs.end();
    if inputs.database {
        let settings = start.settings();
        return match command {
            Commands::Blastp(args) => {
                run_local_database_blastp(args, query, end, settings, outputs)
            }
            Commands::Tblastn(args) => {
                run_local_database_tblastn(args, query, end, settings, outputs)
            }
            Commands::Blastn(args) => {
                run_local_database_blastn(args, query, end, settings, outputs)
            }
            Commands::Tblastx(args) => {
                run_local_database_tblastx(args, query, end, settings, outputs)
            }
            Commands::Blastx(args) => {
                run_local_database_blastx(args, query, end, settings, outputs)
            }
            Commands::Psiblast(args) => {
                let args = psiblast_args(args, outputs)?;
                run_local_database_blastp(args, query, end, settings, outputs)
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
            args.taxdb_setting = taxdb_setting(outputs, start, RequestKind::Blastp);
            // The file of an ungapped search's `-matrix` that NCBI's tables do not have, looked up
            // where the options are resolved, as on the command line (`BlastpArgs::matrix_search`).
            if args.unknown_matrix_name() {
                args.matrix_search = Some(start.settings().matrix_search);
            }
            run_local_blastp(args, query, end, subject, "", "", outputs)
        }
        Commands::Tblastn(mut args) => {
            args.taxdb_setting = taxdb_setting(outputs, start, RequestKind::Tblastn);
            run_local_tblastn(args, query, end, subject, outputs)
        }
        Commands::Blastn(mut args) => {
            args.taxdb_setting = taxdb_setting(outputs, start, RequestKind::Blastn);
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-448
            // The counts file of `-window_masker_db` is read from the file system, as on the
            // command line (LOSAT Web's host gives its own source).
            if args.window_masker_db.is_some() {
                args.window_masker_counts = Some(Arc::new(crate::utils::winmask::FileCounts));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:381-388
            // The directory of `-window_masker_taxid` is looked for by NCBI's library itself, in
            // the environment and `.ncbirc`; LOSAT reads it with the application settings, as
            // a database search of the library does.
            if args.window_masker_taxid.is_some_and(|taxid| taxid != 0) {
                let settings = start.settings();
                args.window_masker_registry = Some(settings.window_masker_registry);
                args.window_masker_counts = Some(Arc::new(crate::utils::winmask::FileCounts));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:368-372
            // The database of `-filtering_db` is looked up as the database of `-db` (the call's
            // settings), by the repeat search of each query batch.
            if args.filtering_db.is_some() {
                let settings = start.settings();
                args.resolve_repeat_filter(&settings);
            }
            run_local_blastn(args, query, end, subject, outputs)
        }
        Commands::Tblastx(mut args) => {
            // The file of a `-matrix` that NCBI's tables do not have, first, as the command
            // line's `run` looks it up (with the call's settings; `TblastxArgs::find_matrix_file`).
            if args.scoring_matrix() == crate::algorithm::tblastx::args::TblastxMatrix::Unknown {
                let settings = start.settings();
                args.find_matrix_file(&settings.matrix_search)?;
            }
            args.taxdb_setting = taxdb_setting(outputs, start, RequestKind::Tblastx);
            run_local_tblastx(args, query, end, subject, outputs)
        }
        Commands::Blastx(mut args) => {
            args.taxdb_setting = taxdb_setting(outputs, start, RequestKind::Blastx);
            // The file of a `-matrix` that NCBI's tables do not have, looked up after the options
            // are resolved, as on the command line (`ResolvedOptions::find_matrix_file`).
            if args.unknown_matrix_name() {
                args.matrix_search = Some(start.settings().matrix_search);
            }
            run_local_blastx(args, query, end, subject, outputs)
        }
        Commands::Psiblast(args) => {
            let mut args = psiblast_args(args, outputs)?;
            args.taxdb_setting = taxdb_setting(outputs, start, RequestKind::Psiblast);
            run_local_blastp(args, query, end, subject, "", "", outputs)
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
/// blastp's arguments for psiblast's search and its rounds (`algorithm::psiblast`,
/// `PsiblastArgs::into_blastp_args`), after LOSAT's refusals of each requested format.
fn psiblast_args(
    args: crate::algorithm::psiblast::PsiblastArgs,
    outputs: &ReportOutputs<'_>,
) -> anyhow::Result<crate::algorithm::blastp::args::BlastpArgs> {
    for format in &outputs.formats {
        // The library writes no search strategy file (`search::validate`).
        crate::algorithm::psiblast::PsiblastArgs::check_outfmt(
            format.outfmt,
            format.side_out(args.out.as_deref()),
            false,
        )?;
    }
    args.into_blastp_args()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
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
            "-parse_deflines needs the records of read_inputs, which reads the ids".into(),
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
    start: &CallStart,
    log: &mut RunLog,
) -> Result<(Vec<Vec<usize>>, usize), Error> {
    use crate::blastinput::app::output_format::ARCHIVE_FORMAT;
    use crate::blastinput::app::parse_formatting_string;
    if formats.is_empty() {
        return Err(Error::Argv("no output format was requested".into()));
    }
    check_inputs(search, inputs)?;
    let numbers: Vec<i32> = formats
        .iter()
        .map(|outfmt| {
            parse_formatting_string(outfmt)
                .map(|choice| choice.number)
                .map_err(|error| {
                    log.post(&error);
                    Error::Argv(super::error::chain_bytes(&error))
                })
        })
        .collect::<Result<_, _>>()?;
    if numbers.contains(&ARCHIVE_FORMAT) {
        check_archive_inputs(search.program(), inputs)?;
    }
    let command = start.command_ref();
    let groups = search_groups(command, &numbers);
    if groups.len() > 1 {
        if let Some(option) = pssm_file_option(command) {
            return Err(Error::Argv(format!(
                "{option} is not accepted by the library API with formats of two hit list sizes, which take two searches: the PSSM files are one search's; search those formats in calls of their own"
            ).into()));
        }
    }
    Ok((groups, num_threads(command)))
}

// No NCBI counterpart: the command line writes one format per run, so one search writes its PSSM
// files; the library's call of formats of two hit list sizes searches twice (`search_groups`).
/// The PSSM file option (`-out_pssm`, else `-out_ascii_pssm`) of a psiblast search that can write
/// one (`PsiblastPssmOptions::writes_files`).
fn pssm_file_option(command: &Commands) -> Option<&'static str> {
    let Commands::Psiblast(args) = command else {
        return None;
    };
    let options = args.clone().into_blastp_args().ok()?.psiblast_pssm?;
    options.writes_files().then(|| {
        if options.out_pssm.is_some() {
            "-out_pssm"
        } else {
            "-out_ascii_pssm"
        }
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:286-288
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
    let (groups, threads) = searches(search, inputs, &formats, &start, &mut log)?;
    let mut writers: Vec<Option<&mut (dyn Write + Send)>> = reports
        .iter_mut()
        .map(|(_, writer)| Some(&mut **writer))
        .collect();
    with_run_pool(threads, search.program().name(), || -> Result<(), Error> {
        for (ordinal, group) in groups.iter().enumerate() {
            let command = start.command();
            // The warnings of the first search include the reader's and the queries' warnings
            // that the next searches repeat; the command line writes them once.
            let mut repeated = Vec::new();
            let mut sink = if ordinal == 0 {
                log.sink(&mut *diagnostics, true)
            } else {
                log.sink(&mut repeated, false)
            };
            // The formats of the other searches, of which the call's refusals are decided
            // (`ReportOutputs::other_outfmts`).
            let others = other_formats(&formats, group);
            let result = {
                let mut outputs = ReportOutputs {
                    other_outfmts: &others,
                    formats: group
                        .iter()
                        .map(|&position| FormatOutput {
                            outfmt: formats[position],
                            refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                            sink: OutputSink::Writer(
                                writers[position].take().expect("one group per format"),
                            ),
                        })
                        .collect(),
                    diagnostics: &mut sink,
                    later_warnings: None,
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
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:255-257
/// Runs the search and writes the report of one `-outfmt` to the file `path`, as the command
/// line's `-out path` does, and the command line's warnings to `diagnostics`. Outfmt 13 and 14
/// come through here: the file is the master file, which names the side files that the search
/// writes next to it (`path_1.json` or `path_1.xml`, ...). The error is the command line's, as in
/// [`run`].
///
/// `path` is the command line's `-out`: its parser checks the name (NCBI's `Illegal value` for a
/// file name of 256 bytes or more, dry run or not), and the search makes the file where the
/// command line makes `-out`, after it has read the subjects (or the taxid or seqid list of `-db`)
/// and before it processes the options: an error before that point leaves no file, an error after
/// it an empty one, and the subjects' warnings come before an error of a file that cannot be made.
/// A dry run makes no file.
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
            "run_to_file writes a file: `run` gives a report to a writer".into(),
        ));
    }
    let search = &search.writing_to(path);
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
        searches(search, inputs, &formats, &start, &mut log)?.1
    };
    let command = start.command();
    let app = search.program().name();
    let mut file = OutFile::new(path, out_file_program(search.program()), app, outfmt);
    let posts = file.posts();
    with_run_pool(threads, app, || -> Result<(), Error> {
        let mut sink = log.sink(&mut *diagnostics, true);
        let result = {
            // The file's posts, in their place among the run's messages (`PostsFirst`).
            let mut messages = PostsFirst::new(posts, &mut sink);
            let mut outputs = ReportOutputs {
                other_outfmts: &[],
                formats: vec![FormatOutput {
                    outfmt,
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Out(&mut file),
                }],
                diagnostics: &mut messages,
                later_warnings: None,
                hits: None,
                observer: None,
            };
            let result = run_local(command, inputs, &mut outputs, &mut start);
            drop(outputs);
            let _ = messages.drain();
            result
        };
        result.map_err(|error| sink.failed(&error))
    })
}

/// The program whose `-out` file names must be UTF-8, as its command line checks them
/// (`cli::create_out_file`): every search program but BLASTX.
fn out_file_program(program: Program) -> Option<&'static str> {
    match program {
        Program::Blastn => Some("BLASTN"),
        Program::Blastp => Some("BLASTP"),
        Program::Psiblast => Some("PSIBLAST"),
        Program::Tblastn => Some("TBLASTN"),
        Program::Tblastx => Some("TBLASTX"),
        Program::Blastx => None,
    }
}

/// What one search of `run_collect` gathered.
struct Searched {
    /// Per format of the group: the report's bytes.
    writers: Vec<CountingWriter>,
    recorder: RangeRecorder,
    /// The hit list of each search call (`Parts`).
    hits: Vec<Arc<Vec<PairwiseHit>>>,
    /// The `-query_loc` interval.
    query_range: Option<SequenceRange>,
    ordinals: Option<Vec<usize>>,
}

/// One search of `run_collect`: the formats at `positions` of `formats`.
fn collect_one(
    inputs: Inputs<'_>,
    formats: &[&str],
    positions: &[usize],
    start: &mut CallStart,
    mut sink: RunSink<'_>,
) -> Result<Searched, Error> {
    let command = start.command();
    let query_range = query_range_of(&command);
    let ordinals = searched_queries(query_range, inputs.query);
    let mut writers: Vec<CountingWriter> =
        positions.iter().map(|_| CountingWriter::new()).collect();
    let mut recorder = RangeRecorder::new(writers.iter().map(|w| w.position.clone()).collect());
    let mut collected: Vec<Arc<Vec<PairwiseHit>>> = Vec::new();
    let others = other_formats(formats, positions);
    let result = {
        let mut on_hits = |hits| collected.push(hits);
        let mut outputs = ReportOutputs {
            other_outfmts: &others,
            formats: positions
                .iter()
                .zip(writers.iter_mut())
                .map(|(&position, writer)| FormatOutput {
                    outfmt: formats[position],
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(writer),
                })
                .collect(),
            diagnostics: &mut sink,
            later_warnings: None,
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
        return Err(Error::Search(
            format!("internal error in the output ranges: {fault}").into(),
        ));
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
// The `-outfmt` columns print a frame for every row (1 for an untranslated row of BLASTP and
// BLASTN, the row's strand for BLASTN's subject, 0 for the protein row of BLASTX and TBLASTN); a
// translation frame exists where the program translates, the rest is the strand.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:880-887
// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
// `shown_id` and `subject_tabular_id` are the first two columns of the default report, whose
// bytes `Id` keeps (`api::id`).
/// The id of a record, read once: many HSPs show it.
fn shared_id(ids: &mut [Option<Id>], index: usize, read: impl FnOnce() -> Vec<u8>) -> Id {
    ids[index].get_or_insert_with(|| Id::new(&read())).clone()
}

/// The HSP records of the final hit list of a search of `inputs`. `qlen` is the record's length
/// or the length of the `-query_loc` interval as given; `slen` is the record's, or the database
/// sequence's (`database_subjects`, which the hit list's subject numbers index).
///
/// `query_base` is the searched query that the hit list numbers 0: a round of psiblast searches
/// one query (`parts::Parts`), a search without rounds numbers the searched queries (0).
pub(super) fn records(
    program: Program,
    inputs: Inputs<'_>,
    hits: &[PairwiseHit],
    query_range: Option<SequenceRange>,
    ordinals: Option<&[usize]>,
    database_subjects: Option<&[ReportedSubject]>,
    query_base: usize,
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
            let q_idx = query_ordinal(ordinals, query_base + hit.q_idx as usize);
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
                            .ok()
                            .flatten()
                            .unwrap_or_else(|| subject.shown_id().to_vec())
                    }),
                    length(inputs.subject, s_idx),
                ),
            };
            HspRecord {
                index,
                q_idx,
                s_idx: subject_index,
                rank: this_rank,
                round: 1,
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
            query_stats: vec![None; inputs.query_count()],
            rounds: vec![
                QueryRounds {
                    stats: vec![None],
                    converged: false,
                };
                inputs.query_count()
            ],
            diagnostics,
        });
    }
    let mut log = start.take_log();
    let (groups, threads) = searches(search, inputs, formats, &start, &mut log)?;
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
            searched.push(collect_one(inputs, formats, group, &mut start, sink)?);
        }
        Ok(())
    })?;
    // The records, statistics and rounds of the first search.
    let finished = match searched.first_mut() {
        Some(first) => std::mem::replace(&mut first.recorder.parts, Parts::new())
            .finish(
                search.program(),
                inputs,
                &first.hits,
                first.query_range,
                first.ordinals.as_deref(),
                0,
            )
            .map_err(|fault| {
                Error::Search(format!("internal error in the search's parts: {fault}").into())
            })?,
        None => Parts::new()
            .finish(search.program(), inputs, &[], None, None, 0)
            .expect("no hit list"),
    };
    let hsps = finished.hsps;
    let mut reports: Vec<Option<Report>> = (0..formats.len()).map(|_| None).collect();
    for (ordinal, (group, done)) in groups.iter().zip(&mut searched).enumerate() {
        for (slot, (&position, writer)) in group.iter().zip(done.writers.drain(..)).enumerate() {
            // The places of the records' HSPs: those of the records' search, by part and HSP.
            let places = |recorded: &mut Vec<Vec<Vec<Option<Range<usize>>>>>| {
                let mut places = Vec::with_capacity(hsps.len());
                if ordinal == 0 {
                    for (part, &length) in finished.lengths.iter().enumerate() {
                        places.extend(RangeRecorder::places(recorded, part, slot, length));
                    }
                }
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
        query_stats: finished.query_stats,
        rounds: finished.rounds,
        diagnostics,
    })
}

#[cfg(test)]
#[path = "run_tests.rs"]
mod tests;
