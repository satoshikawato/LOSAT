//! The search: argv parsing (`validate`) and `run` (docs/web/abi_v2.md §4, §7-8).

use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use LOSAT::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_tblastn, run_local_tblastx, FormatObserver,
    FormatOutput, HspIndex, OutputSink, ReportOutputs,
};
use LOSAT::blastinput::fasta_reader::FastaRecord;
use LOSAT::cli::{Cli, Commands};
use LOSAT::report::PairwiseHit;

use crate::emit::{emit, StreamWriter, STREAM_DIAGNOSTICS, STREAM_HITS};
use crate::json;

/// The programs that ABI v2 runs. BLASTX joins in session SX.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    Blastp,
    Tblastn,
    Blastn,
    Tblastx,
}

impl Program {
    pub fn parse(name: &str) -> Result<Self, String> {
        match name {
            "blastp" => Ok(Self::Blastp),
            "tblastn" => Ok(Self::Tblastn),
            "blastn" => Ok(Self::Blastn),
            "tblastx" => Ok(Self::Tblastx),
            "blastx" => Err("blastx is not available in LOSAT Web ABI v2 yet".into()),
            other => Err(format!("unknown program '{other}'")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Blastp => "blastp",
            Self::Tblastn => "tblastn",
            Self::Blastn => "blastn",
            Self::Tblastx => "tblastx",
        }
    }

    /// The output formats that the engine implements for this program; a run writes
    /// exactly these, each on the stream of the same number.
    pub fn formats(self) -> &'static [u32] {
        match self {
            Self::Blastp | Self::Tblastn | Self::Blastn | Self::Tblastx => &[0, 6, 7],
        }
    }
}

/// Parses an argv (program name first) with the engine's CLI parser.
pub fn parse(words: &[&str]) -> Result<(Program, Commands), String> {
    for word in words.iter().skip(1) {
        let key = word.split('=').next().unwrap_or("");
        if matches!(key, "-out" | "--out" | "-outfmt" | "--outfmt") {
            return Err(format!(
                "{key} is not accepted: a run writes every output format of the program"
            ));
        }
    }
    // The formats come from `Program::formats`, and every program accepts the CLI's
    // default `-outfmt 0`, so the host's words are parsed exactly as on the command line,
    // with the CLI's binary name: an error (an unknown program too) is the CLI's message.
    let argv = std::iter::once("losat").chain(words.iter().copied());
    let cli: Cli =
        LOSAT::cli::try_parse_from(argv).map_err(|error| LOSAT::cli::render_message(&error))?;
    let program = Program::parse(words.first().copied().unwrap_or(""))?;
    Ok((program, cli.command))
}

/// `validate`: `parse`, then the syntax of `-subject_loc` and `-query_loc`, BLASTN's
/// `-dust` value and the scoring options that NCBI rejects or that its Karlin-Altschul
/// tables do not support, with NCBI's message (for a batch of one query), and TBLASTX's
/// checks of its options alone (NCBI's `-evalue` check, the options that LOSAT's TBLASTX
/// rejects).
/// BLASTN's checks include the hit list size, which differs between the searches of a run
/// (`search_groups`): the run is valid when at least one search is (`run` lets the others
/// fail, see there), and the error is that of the first search, the tabular formats'.
/// `run` leaves them to the engine, which reports them as the CLI does.
pub fn validate(words: &[&str]) -> Result<(), String> {
    use LOSAT::blastinput::seq_range::{parse_optional_range, RangeRole};
    // The ranges' syntax, in NCBI's order: the subject range where the subjects are read
    // (before the other option handlers), the query range with the query options (the
    // engine's option checks read it). Whether a range fits the records is checked by
    // `run`, which has them.
    let subject_range = |subject_loc: Option<&str>, program: &str| {
        parse_optional_range(subject_loc, RangeRole::Subject, program).map(|_| ())
    };
    let (program, command) = parse(words)?;
    let groups = search_groups(program.formats(), &command);
    match command {
        Commands::Blastn(mut args) => subject_range(args.subject_loc.as_deref(), "BLASTN")
            .and_then(|()| args.resolve_dust())
            .and_then(|()| {
                parse_optional_range(args.query_loc.as_deref(), RangeRole::Query, "BLASTN")
                    .map(|_| ())
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| {
                let mut first_error = None;
                for group in &groups {
                    match blastn_group_check(&args, group) {
                        Ok(()) => return Ok(()),
                        Err(error) => {
                            first_error.get_or_insert(error);
                        }
                    }
                }
                first_error.map_or(Ok(()), Err)
            })?,
        Commands::Tblastx(args) => subject_range(args.subject_loc.as_deref(), "TBLASTX")
            .and_then(|()| LOSAT::algorithm::tblastx::blast_engine::check_options(&args))
            .map_err(|error| format!("{error:#}"))?,
        Commands::Blastp(args) => subject_range(args.subject_loc.as_deref(), "BLASTP")
            .and_then(|()| LOSAT::algorithm::blastp::blast_engine::check_options(&args))
            .map_err(|error| format!("{error:#}"))?,
        Commands::Tblastn(args) => subject_range(args.subject_loc.as_deref(), "TBLASTN")
            .and_then(|()| LOSAT::algorithm::tblastn::check_options(&args))
            .map_err(|error| format!("{error:#}"))?,
        _ => {}
    }
    Ok(())
}

/// BLASTN's checks (`check_scoring_for_formats`) for the search of the output formats `group`,
/// whose hit list size is the one the CLI derives for them. `validate` runs it for each search
/// of a run, and `run` asks it which search is rejected (`group_rejected`).
fn blastn_group_check(
    args: &LOSAT::algorithm::blastn::BlastnArgs,
    group: &[u32],
) -> Result<(), String> {
    let outfmts: Vec<i32> = group.iter().map(|&format| format as i32).collect();
    LOSAT::algorithm::blastn::scoring::check_scoring_for_formats(args, &outfmts)
        .map_err(|error| format!("{error:#}"))
}

/// Whether the options' checks reject the search of the output formats `group`, as the CLI
/// rejects those formats alone because of their hit list size: BLASTN's checks fail
/// (`blastn_group_check`: a hit list of 0, or a preliminary hit list size below 1), or, for
/// the other programs, the hit list size is 0. `run` tolerates a failed search of a
/// two-search run only for this reason (`fails_run`).
fn group_rejected(command: &Commands, group: &[u32]) -> bool {
    match command {
        Commands::Blastn(args) => blastn_group_check(args, group).is_err(),
        _ => hit_list_size(command, group[0]) == 0,
    }
}

/// Whether the failure of one search fails the run. A failed search of a run of two
/// searches is tolerated only when it wrote no byte of its formats and `rejected` says that
/// the options' checks reject its hit list size (`group_rejected`); any other error (for
/// example a thread pool that cannot start) fails the run. When the other search fails too,
/// `run` fails with the tabular search's error.
fn fails_run(failed: &Failed, rejected: bool) -> bool {
    failed.wrote || !rejected
}

/// Records the byte range of every HSP row (stream 6) and section (stream 0), and of
/// the subject heading (stream 0) that each section belongs to.
struct RangeRecorder {
    streams: Vec<u32>,
    positions: Vec<Arc<AtomicU64>>,
    open: Vec<Option<(HspIndex, u64)>>,
    out6: Vec<Option<(u64, u64)>>,
    out0: Vec<Option<(u64, u64)>>,
    out0_subject: Vec<Option<(u64, u64)>>,
    /// The start of the heading being written, and the last complete heading.
    subject_start: Option<u64>,
    subject: Option<(u64, u64)>,
    /// The first observer event that does not fit the formatter contract; the run fails.
    fault: Option<String>,
}

impl RangeRecorder {
    fn position(&self, format: usize) -> u64 {
        self.positions[format].load(Ordering::SeqCst)
    }
}

/// Stores the range of `hsp`; a second range for the same HSP is a fault.
fn store(
    ranges: &mut Vec<Option<(u64, u64)>>,
    hsp: HspIndex,
    range: (u64, u64),
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
            self.fault
                .get_or_insert(format!("HSP {hsp} started before HSP {open} ended"));
        }
    }

    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        let Some((begun, start)) = self.open[format].take() else {
            self.fault
                .get_or_insert(format!("HSP {hsp} ended without a start"));
            return;
        };
        if begun != hsp {
            self.fault
                .get_or_insert(format!("HSP {begun} started but HSP {hsp} ended"));
            return;
        }
        let range = (start, self.position(format));
        let stored = match self.streams[format] {
            6 => store(&mut self.out6, hsp, range),
            0 => store(&mut self.out0, hsp, range).and_then(|()| match self.subject {
                Some(subject) => store(&mut self.out0_subject, hsp, subject),
                None => Ok(()),
            }),
            _ => Ok(()),
        };
        if let Err(fault) = stored {
            self.fault.get_or_insert(fault);
        }
    }

    fn subject_begin(&mut self, format: usize, _first_hsp: HspIndex) {
        if self.streams[format] == 0 {
            self.subject_start = Some(self.position(format));
        }
    }

    fn subject_end(&mut self, format: usize, _first_hsp: HspIndex) {
        if self.streams[format] == 0 {
            if let Some(start) = self.subject_start.take() {
                self.subject = Some((start, self.position(format)));
            }
        }
    }
}

/// The registered index of each query that the engine searches: NCBI's batch reader skips
/// a query whose `-query_loc` interval starts more than one letter past its end
/// (`seq_range::cut_queries`), so the engine's `q_idx` counts the searched queries; the
/// HSP records give the registered index (docs/web/abi_v2.md §8, port plan Q8). `None`
/// without `-query_loc` (every query is searched); a range that does not parse fails the
/// run before any hit.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155
/// ```c++
///         try { q.Reset(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
///                 break;
///             }
///             throw;
///         }
///         catch (const exception&) {
///             continue; //SB-2307. ignore well formed, not found accession
///         }
/// ```
fn searched_queries(
    query_loc: Option<&str>,
    program: &str,
    records: &[FastaRecord],
) -> Option<Vec<usize>> {
    use LOSAT::blastinput::seq_range::{
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

/// The options that set the hit list size: `-max_target_seqs`, `-num_descriptions` and
/// `-num_alignments` of the program's arguments.
fn count_options(command: &Commands) -> [Option<usize>; 3] {
    match command {
        Commands::Blastp(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Tblastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Blastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Tblastx(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        // Rejected by `Program::parse`.
        Commands::Blastx(_) | Commands::Makedb(_) => [None; 3],
    }
}

/// The hit list size that the CLI derives for the output format `format` from the options
/// of `command` (`report_counts`, NCBI `blast_args.cpp:2894-2973`).
pub fn hit_list_size(command: &Commands, format: u32) -> usize {
    let [max_target_seqs, descriptions, alignments] = count_options(command);
    LOSAT::blastinput::report_counts::report_counts(
        format as i32,
        max_target_seqs,
        descriptions,
        alignments,
    )
    .hitlist_size
}

/// The searches of a run: the output formats grouped by hit list size, the tabular
/// formats first. One group is one `run_local` call, and every format of it is written from
/// that search, as the CLI writes it alone: the CLI's search depends on the hit list size,
/// which `-num_descriptions` and `-num_alignments` (without `-max_target_seqs`) can make
/// differ between the pairwise format 0 and the tabular formats. Without those options the
/// sizes are equal and the group is the program's formats in their order (the one search of
/// before). The tabular group goes first because the formatting handler's warnings belong to
/// the first format, and the CLI writes the tabular formats' (`ReportCounts::warnings`).
pub fn search_groups(formats: &[u32], command: &Commands) -> Vec<Vec<u32>> {
    let [_, descriptions, alignments] = count_options(command);
    if descriptions.is_none() && alignments.is_none() {
        return vec![formats.to_vec()];
    }
    let mut ordered: Vec<u32> = formats.to_vec();
    ordered.sort_by_key(|&format| format <= 4);
    let mut groups: Vec<(usize, Vec<u32>)> = Vec::new();
    for format in ordered {
        let size = hit_list_size(command, format);
        match groups.iter_mut().find(|(group, _)| *group == size) {
            Some((_, group)) => group.push(format),
            None => groups.push((size, vec![format])),
        }
    }
    groups.into_iter().map(|(_, group)| group).collect()
}

/// A search that failed.
struct Failed {
    /// The error of the run, when this failure fails it.
    message: String,
    /// What the CLI writes to standard error for the error: NCBI's message (a
    /// `NativeError`), or `Error: ` and the error, as `main` returns it.
    cli_stderr: Vec<u8>,
    /// Whether the search had written a byte of any of its formats. The host discards
    /// the streams of a failed run, but a search that wrote nothing leaves the other
    /// searches' streams as the CLI writes them.
    wrote: bool,
}

impl Failed {
    /// `error` is the engine's `anyhow::Error`, which this crate does not name.
    fn new<E>(error: &E, wrote: bool) -> Self
    where
        E: std::ops::Deref<Target = dyn std::error::Error + Send + Sync + 'static>
            + std::fmt::Display
            + std::fmt::Debug,
    {
        let cli_stderr = match error.downcast_ref::<LOSAT::cli::NativeError>() {
            Some(native) => native.message.clone().into_bytes(),
            None => format!("Error: {error:?}\n").into_bytes(),
        };
        Self {
            message: format!("{error:#}"),
            cli_stderr,
            wrote,
        }
    }
}

/// What one search wrote.
struct Searched {
    writers: Vec<StreamWriter>,
    recorder: RangeRecorder,
    hits: Option<Vec<PairwiseHit>>,
    ordinals: Option<Vec<usize>>,
}

/// One `run_local` call that writes the output formats `formats` (a group of
/// `search_groups`), with the warnings of the CLI to `diagnostics`.
fn search(
    command: Commands,
    formats: &[u32],
    query: &[FastaRecord],
    subject: &[FastaRecord],
    diagnostics: &mut (dyn std::io::Write + Send),
) -> Result<Searched, Failed> {
    let mut writers: Vec<StreamWriter> = formats.iter().map(|&f| StreamWriter::new(f)).collect();
    let mut recorder = RangeRecorder {
        streams: formats.to_vec(),
        positions: writers.iter().map(StreamWriter::position).collect(),
        open: vec![None; formats.len()],
        out6: Vec::new(),
        out0: Vec::new(),
        out0_subject: Vec::new(),
        subject_start: None,
        subject: None,
        fault: None,
    };
    let outfmts: Vec<String> = formats.iter().map(u32::to_string).collect();
    let mut collected: Option<Vec<PairwiseHit>> = None;
    let ordinals;
    let outcome = {
        let mut on_hits = |hits: &[PairwiseHit]| collected = Some(hits.to_vec());
        let mut outputs = ReportOutputs {
            formats: outfmts
                .iter()
                .zip(writers.iter_mut())
                .map(|(outfmt, writer)| FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(writer),
                })
                .collect(),
            diagnostics,
            hits: Some(&mut on_hits),
            observer: Some(&mut recorder),
        };
        let result = match command {
            Commands::Blastp(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), "BLASTP", query);
                run_local_blastp(args, query, subject, "", "", &mut outputs)
            }
            Commands::Tblastn(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), "TBLASTN", query);
                run_local_tblastn(args, query, subject, &mut outputs)
            }
            Commands::Blastn(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), "BLASTN", query);
                run_local_blastn(args, query, subject, &mut outputs)
            }
            Commands::Tblastx(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), "TBLASTX", query);
                run_local_tblastx(args, query, subject, &mut outputs)
            }
            Commands::Blastx(_) | Commands::Makedb(_) => {
                unreachable!("rejected by Program::parse")
            }
        };
        result
    };
    if let Err(error) = outcome {
        let wrote = writers
            .iter()
            .any(|writer| writer.position().load(Ordering::SeqCst) > 0);
        return Err(Failed::new(&error, wrote));
    }
    if let Some((hsp, _)) = recorder.open.iter().flatten().next() {
        recorder
            .fault
            .get_or_insert(format!("HSP {hsp} started but did not end"));
    }
    if let Some(fault) = recorder.fault.take() {
        let message = format!("internal error in the output ranges: {fault}");
        return Err(Failed {
            cli_stderr: message.clone().into_bytes(),
            message,
            wrote: true,
        });
    }
    Ok(Searched {
        writers,
        recorder,
        hits: collected,
        ordinals,
    })
}

/// The HSPs of one (query, subject) pair, in the list's order.
type Pairs = std::collections::HashMap<(u32, u32), Vec<usize>>;

fn pairs(hits: &[PairwiseHit], only: impl Fn(usize) -> bool) -> Pairs {
    let mut pairs = Pairs::new();
    for (index, record) in hits.iter().enumerate() {
        if only(index) {
            pairs
                .entry((record.hit.q_idx, record.hit.s_idx))
                .or_default()
                .push(index);
        }
    }
    pairs
}

/// Whether two searches found the same HSP: every field of the `PairwiseHit`, the numbers
/// of `f64` by their bits.
fn same_hsp(a: &PairwiseHit, b: &PairwiseHit) -> bool {
    let (x, y) = (&a.hit, &b.hit);
    x.identity.to_bits() == y.identity.to_bits()
        && x.length == y.length
        && x.mismatch == y.mismatch
        && x.gapopen == y.gapopen
        && (x.q_start, x.q_end, x.s_start, x.s_end) == (y.q_start, y.q_end, y.s_start, y.s_end)
        && x.e_value.to_bits() == y.e_value.to_bits()
        && x.bit_score.to_bits() == y.bit_score.to_bits()
        && x.num_ident == y.num_ident
        && x.query_frame == y.query_frame
        && x.query_length == y.query_length
        && (x.q_idx, x.s_idx, x.raw_score) == (y.q_idx, y.s_idx, y.raw_score)
        && (x.sort_query_offset, x.sort_query_end) == (y.sort_query_offset, y.sort_query_end)
        && (x.sort_subject_offset, x.sort_subject_end)
            == (y.sort_subject_offset, y.sort_subject_end)
        && x.has_sort_offsets == y.has_sort_offsets
        && x.gap_info == y.gap_info
        && x.num_positives == y.num_positives
        && a.query_seq == b.query_seq
        && a.subject_seq == b.subject_seq
        && a.query_frame == b.query_frame
        && a.subject_frame == b.subject_frame
        && a.positives == b.positives
        && a.gaps == b.gaps
        && a.subject_length == b.subject_length
        && a.subject_title == b.subject_title
        && a.comp_adjust_method == b.comp_adjust_method
        && a.sum_n == b.sum_n
}

/// Gives the HSPs of the tabular search `first` the sections of format 0 from the other
/// search, `second`, for the pairs (query, subject) whose HSPs in `second`'s format 0 are
/// exactly those of `first`: the same count and, in order, every field of the HSP equal
/// (`same_hsp`). A format 0 hit list
/// of another size can keep other subjects or other HSPs of a subject (the preliminary
/// list evicts a subject before its alignment is extended), and such a pair stays without
/// a section. Records, rows and sections that this leaves unlinked are still correct.
fn link_format0(
    first: &[PairwiseHit],
    first_recorder: &mut RangeRecorder,
    second: &[PairwiseHit],
    second_recorder: &RangeRecorder,
) {
    let headed = pairs(second, |index| {
        second_recorder.out0.get(index).copied().flatten().is_some()
    });
    let ours = pairs(first, |_| true);
    let same = same_hsp;
    for (pair, theirs) in headed {
        let Some(mine) = ours.get(&pair) else {
            continue;
        };
        if mine.len() != theirs.len()
            || !mine
                .iter()
                .zip(&theirs)
                .all(|(&a, &b)| same(&first[a], &second[b]))
        {
            continue;
        }
        for (&a, &b) in mine.iter().zip(&theirs) {
            for (target, source) in [
                (&mut first_recorder.out0, &second_recorder.out0),
                (
                    &mut first_recorder.out0_subject,
                    &second_recorder.out0_subject,
                ),
            ] {
                if target.len() <= a {
                    target.resize(a + 1, None);
                }
                target[a] = source.get(b).copied().flatten();
            }
        }
    }
}

/// Runs the searches of the run and emits its streams. On failure nothing more is emitted
/// and the host discards what it received.
///
/// Without `-num_descriptions` and `-num_alignments` this is one search that writes every
/// format. With them it is one search per distinct hit list size of the formats
/// (`search_groups`), so that each format's stream is the CLI's output for that format
/// alone. The HSP records follow the tabular search; a record gets the section of format 0
/// only where the search of format 0 has exactly its HSPs (`link_format0`)
/// (docs/web/abi_v2.md §5, §7, §8).
///
/// A search of a two-search run that fails before it writes a byte of its formats, and
/// that the options' checks reject because of its hit list size (`group_rejected`:
/// `-num_alignments 0` for the tabular formats, a BLASTN `-num_descriptions` of 2^30 or more
/// for format 0), does
/// not fail the run when the other search succeeds, as the CLI fails that format alone: the
/// streams of its formats stay empty (with no tabular search there are no records). Any other
/// error of a search (for example a thread pool that cannot start) fails the run, as does a
/// failure after a byte was written; so do both searches failing, as the CLI fails every
/// format.
///
/// Stream 3 is the CLI's standard error for the tabular formats (the tabular search's
/// diagnostics, with its error when it failed). When a search failed it is followed by the
/// CLI's standard error for format 0 (the diagnostics of the format 0 search, with its
/// error when it failed). The diagnostics of the format 0 search of a run without a failure
/// repeat the first search's and are dropped.
pub fn run(words: &[&str], query: &[FastaRecord], subject: &[FastaRecord]) -> Result<(), String> {
    use std::io::Write as _;
    let (program, command) = parse(words)?;
    let groups = search_groups(program.formats(), &command);
    let mut diagnostics = StreamWriter::new(STREAM_DIAGNOSTICS);
    let mut command = Some(command);
    // Each search: what it wrote or why it failed, and (after the first) its diagnostics.
    let mut outcomes: Vec<(Result<Searched, Failed>, Vec<u8>)> = Vec::new();
    for group in &groups {
        let command = match command.take() {
            Some(command) => command,
            None => parse(words)?.1,
        };
        let mut scratch = Vec::new();
        // Whether the options' checks reject this search, found before `search` takes the
        // command.
        let rejected = group_rejected(&command, group);
        // The warnings of the first search include the reader's and the queries' warnings
        // that the next searches repeat; the CLI writes them once.
        let result = if outcomes.is_empty() {
            search(command, group, query, subject, &mut diagnostics)
        } else {
            search(command, group, query, subject, &mut scratch)
        };
        if let Err(failed) = &result {
            if fails_run(failed, rejected) {
                return Err(failed.message.clone());
            }
        }
        outcomes.push((result, scratch));
    }
    if outcomes.iter().all(|(result, _)| result.is_err()) {
        // The CLI fails every format; the error is that of the tabular formats.
        return Err(match outcomes.remove(0).0 {
            Err(failed) => failed.message,
            Ok(_) => unreachable!("every search failed"),
        });
    }
    if outcomes.iter().any(|(result, _)| result.is_err()) {
        for (position, (result, scratch)) in outcomes.iter().enumerate() {
            if position > 0 {
                let _ = diagnostics.write_all(scratch);
            }
            if let Err(failed) = result {
                let _ = diagnostics.write_all(&failed.cli_stderr);
            }
        }
    }
    let mut searched: Vec<Searched> = Vec::new();
    let mut tabular_failed = false;
    for (position, (result, _)) in outcomes.into_iter().enumerate() {
        match result {
            Ok(done) => searched.push(done),
            Err(_) => tabular_failed |= position == 0,
        }
    }
    let mut writers = Vec::new();
    let mut records = None;
    if !tabular_failed {
        let mut first = searched.remove(0);
        if let Some(hits) = &first.hits {
            for other in &searched {
                if let Some(second) = &other.hits {
                    link_format0(hits, &mut first.recorder, second, &other.recorder);
                }
            }
        }
        writers.append(&mut first.writers);
        records = first
            .hits
            .map(|hits| (hits, first.recorder, first.ordinals));
    }
    for other in &mut searched {
        writers.append(&mut other.writers);
    }
    for writer in &mut writers {
        writer.finish();
    }
    diagnostics.finish();
    if let Some((hits, recorder, ordinals)) = records {
        emit_hit_records(&hits, &recorder, ordinals.as_deref());
    }
    Ok(())
}

/// Emits one JSON object per HSP on stream 1 (docs/web/abi_v2.md §8); `ordinals` maps the
/// engine's `q_idx` (searched queries) to the registered index (`searched_queries`).
fn emit_hit_records(hits: &[PairwiseHit], recorder: &RangeRecorder, ordinals: Option<&[usize]>) {
    let mut lines = String::new();
    let mut ranks: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for (index, record) in hits.iter().enumerate() {
        let hit = &record.hit;
        let q_idx = ordinals
            .and_then(|ordinals| ordinals.get(hit.q_idx as usize).copied())
            .unwrap_or(hit.q_idx as usize);
        let rank = ranks.entry(q_idx).or_default();
        let _ = write!(
            lines,
            "{{\"index\":{index},\"q_idx\":{q_idx},\"s_idx\":{},\"rank\":{rank},\"raw_score\":{},\"bit_score\":",
            hit.s_idx, hit.raw_score
        );
        *rank += 1;
        json::number(&mut lines, hit.bit_score);
        lines.push_str(",\"e_value\":");
        json::number(&mut lines, hit.e_value);
        let _ = write!(
            lines,
            ",\"q_start\":{},\"q_end\":{},\"s_start\":{},\"s_end\":{},\"query_frame\":",
            hit.q_start, hit.q_end, hit.s_start, hit.s_end
        );
        json::optional(&mut lines, record.query_frame);
        lines.push_str(",\"subject_frame\":");
        json::optional(&mut lines, record.subject_frame);
        lines.push_str(",\"subject_length\":");
        json::optional(&mut lines, record.subject_length);
        lines.push_str(",\"query_aligned\":");
        json::optional_string(&mut lines, record.query_seq.as_deref());
        lines.push_str(",\"subject_aligned\":");
        json::optional_string(&mut lines, record.subject_seq.as_deref());
        lines.push_str(",\"out6\":");
        json::range(&mut lines, recorder.out6.get(index).copied().flatten());
        lines.push_str(",\"out0\":");
        json::range(&mut lines, recorder.out0.get(index).copied().flatten());
        lines.push_str(",\"out0_subject\":");
        json::range(
            &mut lines,
            recorder.out0_subject.get(index).copied().flatten(),
        );
        lines.push_str("}\n");
        if lines.len() >= crate::emit::CHUNK {
            emit(STREAM_HITS, lines.as_bytes());
            lines.clear();
        }
    }
    emit(STREAM_HITS, lines.as_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn argv_is_parsed_as_on_the_command_line() {
        let words = ["blastp", "-query", "q", "-subject", "s", "-evalue", "0.1"];
        assert!(matches!(
            parse(&words),
            Ok((Program::Blastp, Commands::Blastp(_)))
        ));
        for owned in ["-outfmt", "-out", "-outfmt=6"] {
            let error = parse(&["blastp", "-query", "q", "-subject", "s", owned, "6"]).unwrap_err();
            assert!(error.contains("is not accepted"), "{owned}: {error}");
        }
        // -num_descriptions and -num_alignments are accepted (also as `-x=v`),
        // parsed as on the command line; -max_target_seqs excludes both (clap's conflict).
        for program in ["blastp", "tblastn", "blastn", "tblastx"] {
            for option in [
                &["-num_descriptions", "3"][..],
                &["-num_alignments", "0"],
                &["-num_alignments=3"],
            ] {
                let words = [&[program, "-query", "q", "-subject", "s"][..], option].concat();
                assert!(parse(&words).is_ok(), "{program} {option:?}");
                assert!(validate(&words).is_ok(), "{program} {option:?}");
            }
            let words = [
                program,
                "-query",
                "q",
                "-subject",
                "s",
                "-max_target_seqs",
                "3",
            ];
            assert!(parse(&words).is_ok(), "{program}");
            for option in ["-num_descriptions", "-num_alignments"] {
                let both = [&words[..], &[option, "3"]].concat();
                let cli = LOSAT::cli::try_parse_from::<Cli, _, _>(
                    std::iter::once("losat").chain(both.iter().copied()),
                )
                .unwrap_err();
                assert_eq!(
                    parse(&both).unwrap_err(),
                    LOSAT::cli::render_message(&cli),
                    "{program} {option}"
                );
                assert!(parse(&both).unwrap_err().contains("cannot be used with"));
            }
        }
        assert!(parse(&["blastx", "-query", "q", "-subject", "s"]).is_err());
        // An unknown program and a parser error are the CLI's messages.
        for words in [&["nosuch"][..], &["tblastx", "-query", "q"][..]] {
            let cli = LOSAT::cli::try_parse_from::<Cli, _, _>(
                std::iter::once("losat").chain(words.iter().copied()),
            )
            .map(|_| ())
            .map_err(|error| LOSAT::cli::render_message(&error));
            assert_eq!(parse(words).map(|_| ()), cli, "{words:?}");
        }
        assert_eq!(Program::Tblastx.formats(), &[0, 6, 7]);
        // validate rejects BLASTN scoring that NCBI rejects, with NCBI's message.
        // validate reads BLASTN's -dust as the run does.
        let error =
            validate(&["blastn", "-query", "q", "-subject", "s", "-dust", "20 64"]).unwrap_err();
        assert!(
            error.starts_with("BLAST query/options error: Invalid number of arguments"),
            "{error}"
        );
        let blastn = ["blastn", "-query", "q", "-subject", "s", "-task", "blastn"];
        assert!(validate(&blastn).is_ok());
        let unsupported = [&blastn[..], &["-reward", "1", "-penalty", "-6"]].concat();
        assert!(
            parse(&unsupported).is_ok(),
            "run leaves the check to the engine"
        );
        let error =
            validate(&[&blastn[..], &["-gapopen", "0", "-gapextend", "0"]].concat()).unwrap_err();
        assert!(
            error.starts_with("BLAST query/options error: Greedy extension"),
            "{error}"
        );
        let error = validate(&unsupported).unwrap_err();
        assert!(
            error.starts_with("BLAST engine error: Error: Substitution scores 1 and -6"),
            "{error}"
        );
        // dc-megablast and blastn-short (Session SD): the task defaults and NCBI's template
        // checks, with NCBI's messages (blast_options.c:1247-1261,1399-1413).
        for task in ["dc-megablast", "blastn-short"] {
            assert!(validate(&["blastn", "-query", "q", "-subject", "s", "-task", task]).is_ok());
        }
        let dc = [
            "blastn",
            "-query",
            "q",
            "-subject",
            "s",
            "-task",
            "dc-megablast",
        ];
        let template = [
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "21",
        ];
        assert!(validate(&[&dc[..], &template[..], &["-word_size", "12"]].concat()).is_ok());
        let error = validate(&[&dc[..], &["-word_size", "13"]].concat()).unwrap_err();
        assert!(
            error.starts_with(
                "BLAST query/options error: Invalid discontiguous template parameters: word size must be either 11 or 12"
            ),
            "{error}"
        );
        let error =
            validate(&[&blastn[..], &template[..], &["-word_size", "11"]].concat()).unwrap_err();
        assert!(
            error.starts_with(
                "BLAST query/options error: Invalid lookup table type for discontiguous Mega BLAST"
            ),
            "{error}"
        );
        // The options require each other: a parser error, the CLI's message.
        assert!(validate(&[&dc[..], &["-template_type", "coding"]].concat()).is_err());
        // validate reports TBLASTX's checks of its options alone, as the run does.
        let tblastx = ["tblastx", "-query", "q", "-subject", "s"];
        assert!(validate(&tblastx).is_ok());
        for (extra, start) in [
            (["-evalue", "0"], "BLAST query/options error: expect value"),
            (
                ["-window_size", "0"],
                "-window_size 0 (the one-hit word finder)",
            ),
        ] {
            let error = validate(&[&tblastx[..], &extra[..]].concat()).unwrap_err();
            assert!(error.starts_with(start), "{error}");
        }
        // An error of the CLI parser keeps the CLI's message.
        for words in [
            &["blastp", "-query", "q", "-subject", "s", "-evalue"][..],
            &[
                "blastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-nosuchoption",
                "1",
            ][..],
            &[
                "tblastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-db_gencode",
                "7",
            ][..],
        ] {
            let cli = LOSAT::cli::try_parse_from::<Cli, _, _>(
                ["losat"].into_iter().chain(words.iter().copied()),
            )
            .unwrap_err();
            assert_eq!(
                parse(words).unwrap_err(),
                LOSAT::cli::render_message(&cli),
                "{words:?}"
            );
        }
    }

    fn recorder() -> RangeRecorder {
        RangeRecorder {
            streams: vec![0, 6],
            positions: vec![Arc::new(AtomicU64::new(0)), Arc::new(AtomicU64::new(0))],
            open: vec![None, None],
            out6: Vec::new(),
            out0: Vec::new(),
            out0_subject: Vec::new(),
            subject_start: None,
            subject: None,
            fault: None,
        }
    }

    // LOSAT Web gives the engine subject records: BLASTP's `-db` (a database read from
    // files, CLI only until LOSAT Web takes a database) is rejected by `validate` and by
    // `run`, before any search.
    #[test]
    fn a_database_search_is_rejected() {
        let words = ["blastp", "-query", "q", "-db", "swissprot"];
        let error = validate(&words).unwrap_err();
        assert_eq!(
            error,
            "-db (a database search) is not supported by LOSAT Web yet"
        );
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [FastaRecord::new("Query_1", b"q", protein)];
        let subject = [FastaRecord::new("Subject_1", b"s", protein)];
        assert_eq!(
            run(&words, &query, &subject).unwrap_err(),
            "-db (a database search) is not supported by LOSAT Web yet"
        );
        // `-db` with `-subject` is the parser's conflict, as on the command line.
        let both = ["blastp", "-query", "q", "-subject", "s", "-db", "swissprot"];
        assert!(validate(&both).unwrap_err().contains("cannot be used with"));
        // `-dbsize` is an option of the search: accepted with the subject records.
        validate(&[
            "blastp", "-query", "q", "-subject", "s", "-dbsize", "1000000",
        ])
        .unwrap();
    }

    // The HSP records give the registered index of a query where `-query_loc` skips
    // an earlier record (its interval starts past its end), and records without residues
    // are registered records too.
    #[test]
    fn q_idx_is_the_registered_index() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [
            FastaRecord::new("Query_1", b"short", b"MK"),
            FastaRecord::new("Query_2", b"empty", b""),
            FastaRecord::new("Query_3", b"q", protein),
        ];
        let subject = [FastaRecord::new("Subject_1", b"s", protein)];
        let q_idx = |words: &[&str]| {
            crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
            run(words, &query[1..], &subject).unwrap();
            let hits = crate::emit::EMITTED.with(|emitted| {
                emitted
                    .borrow()
                    .iter()
                    .filter(|(stream, _)| *stream == crate::emit::STREAM_HITS)
                    .flat_map(|(_, bytes)| bytes.clone())
                    .collect::<Vec<u8>>()
            });
            let hits = String::from_utf8(hits).unwrap();
            assert!(!hits.is_empty(), "{words:?}");
            hits.lines()
                .map(|line| {
                    let at = line.find("\"q_idx\":").unwrap() + 8;
                    line[at..]
                        .split(',')
                        .next()
                        .unwrap()
                        .parse::<usize>()
                        .unwrap()
                })
                .collect::<Vec<_>>()
        };
        // The query without residues is searched as a query without data (engine q_idx 1);
        // with the range it is skipped (engine q_idx 0), and the record stays the second.
        let words = ["blastp", "-query", "q", "-subject", "s"];
        assert!(q_idx(&words).iter().all(|&index| index == 1));
        let ranged = [&words[..], &["-query_loc", "5-40"]].concat();
        assert!(q_idx(&ranged).iter().all(|&index| index == 1));
        let ordinals = searched_queries(Some("5-40"), "BLASTP", &query).unwrap();
        assert_eq!(ordinals, vec![2]);
        assert_eq!(searched_queries(None, "BLASTP", &query), None);
    }

    fn fasta(text: &[u8]) -> Vec<FastaRecord> {
        let text = std::str::from_utf8(text).unwrap();
        text.split('>')
            .filter(|entry| !entry.is_empty())
            .map(|entry| {
                let (name, residues) = entry.split_once('\n').unwrap();
                FastaRecord::new(name, name.as_bytes(), residues.trim_end().as_bytes())
            })
            .collect()
    }

    /// Bytes as text, so that a failed comparison prints readable text.
    fn text(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).into_owned()
    }

    /// The bytes of each stream that `run` emitted, joined.
    fn streams(
        words: &[&str],
        query: &[FastaRecord],
        subject: &[FastaRecord],
    ) -> Result<std::collections::BTreeMap<u32, Vec<u8>>, String> {
        crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
        run(words, query, subject)?;
        let mut streams = std::collections::BTreeMap::new();
        crate::emit::EMITTED.with(|emitted| {
            for (stream, bytes) in emitted.borrow().iter() {
                streams
                    .entry(*stream)
                    .or_insert_with(Vec::new)
                    .extend_from_slice(bytes);
            }
        });
        Ok(streams)
    }

    // The searches of a run: one per distinct hit list size of the formats 0, 6 and 7 (NCBI
    // `blast_args.cpp:2894-2973`), the tabular formats first; without -num_descriptions and
    // -num_alignments the one search of the program's formats in their order.
    #[test]
    fn a_run_searches_once_for_each_hit_list_size() {
        let groups = |program: &str, options: &[&str]| {
            let words = [&[program, "-query", "q", "-subject", "s"][..], options].concat();
            let (program, command) = parse(&words).unwrap();
            search_groups(program.formats(), &command)
        };
        let one = vec![vec![0, 6, 7]];
        let one_tabular_first = vec![vec![6, 7, 0]];
        let two = vec![vec![6, 7], vec![0]];
        for program in ["blastp", "tblastn", "blastn", "tblastx"] {
            for (options, expected) in [
                (&[][..], &one),
                (&["-max_target_seqs", "3"], &one),
                // The hit list of format 0 is max(-num_descriptions, -num_alignments), with
                // 500 and 250 for the one that is not given; the tabular formats' is
                // -num_alignments, or 500.
                (&["-num_descriptions", "500"], &one_tabular_first),
                (&["-num_descriptions", "3"], &two),
                (&["-num_alignments", "600"], &one_tabular_first),
                (&["-num_alignments", "3"], &two),
                (&["-num_alignments", "0"], &two),
                (
                    &["-num_descriptions", "2", "-num_alignments", "4"],
                    &one_tabular_first,
                ),
                (&["-num_descriptions", "4", "-num_alignments", "2"], &two),
                (
                    &["-num_descriptions", "0", "-num_alignments", "0"],
                    &one_tabular_first,
                ),
            ] {
                assert_eq!(&groups(program, options), expected, "{program} {options:?}");
            }
        }
    }

    const SPLIT_QUERY: &[u8] =
        include_bytes!("../../../LOSAT/tests/fasta/outfmt0/numdesc_split_q.fna");
    const SPLIT_SUBJECT: &[u8] =
        include_bytes!("../../../LOSAT/tests/fasta/outfmt0/numdesc_split_s.fna");
    const SPLIT_OUT0: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_na3.out");
    const SPLIT_OUT6: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_na3.6.out");
    const SPLIT_ERR6: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_na3.6.err");

    // The one input where outfmt 0 and outfmt 6/7 keep other subjects with -num_alignments 3
    // (the hit list of outfmt 0 is 500, the tabular formats' 3, and the subject A_split, whose
    // gap the preliminary X-drop cannot cross, is evicted from the small list): two searches.
    // Each stream is NCBI's frozen bytes of its format alone; the HSP records are the
    // tabular list, and only the HSPs that outfmt 0 shows identically have its sections.
    #[test]
    fn two_searches_write_each_format_as_the_cli_does() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        let words = [
            "blastn",
            "-query",
            "tests/fasta/outfmt0/numdesc_split_q.fna",
            "-subject",
            "tests/fasta/outfmt0/numdesc_split_s.fna",
            "-task",
            "blastn",
            "-num_alignments",
            "3",
        ];
        let streams = streams(&words, &query, &subject).unwrap();
        assert_eq!(streams[&0], SPLIT_OUT0);
        assert_eq!(streams[&6], SPLIT_OUT6);
        // The warning of the tabular formats, once (outfmt 0 has none).
        assert_eq!(streams[&3], SPLIT_ERR6);
        let out0 = String::from_utf8(streams[&0].clone()).unwrap();
        let headings: Vec<&str> = out0.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headings, ["> A_split", "> D00", "> D01"]);
        let records = String::from_utf8(streams[&1].clone()).unwrap();
        let field = |line: &str, name: &str| {
            let at = line.find(&format!("\"{name}\":")).unwrap() + name.len() + 3;
            line[at..]
                .split(&[',', '}'][..])
                .next()
                .unwrap()
                .to_string()
        };
        let lines: Vec<&str> = records.lines().collect();
        // D00 D01 D02 are the subjects 1, 2 and 3; D02 is not in outfmt 0's list.
        assert_eq!(
            lines.iter().map(|l| field(l, "s_idx")).collect::<Vec<_>>(),
            ["1", "2", "3"]
        );
        let pointed: Vec<bool> = lines.iter().map(|l| !l.contains("\"out0\":null")).collect();
        assert_eq!(pointed, [true, true, false]);
        assert!(lines[2].contains("\"out0_subject\":null"));
        let out0_subject = |line: &str| {
            let at = line.find("\"out0_subject\":[").unwrap() + 16;
            let range: Vec<usize> = line[at..]
                .split(']')
                .next()
                .unwrap()
                .split(',')
                .map(|n| n.parse().unwrap())
                .collect();
            out0[range[0]..range[1]].lines().next().unwrap().to_string()
        };
        assert_eq!(out0_subject(lines[0]), "> D00");
        assert_eq!(out0_subject(lines[1]), "> D01");
    }

    // The same command line without a difference of the hit list sizes is one search.
    #[test]
    fn equal_hit_list_sizes_are_one_search() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        fn words<'a>(extra: &[&'a str]) -> Vec<&'a str> {
            [
                &["blastn", "-query", "q", "-subject", "s", "-task", "blastn"][..],
                extra,
            ]
            .concat()
        }
        let plain = streams(&words(&[]), &query, &subject).unwrap();
        // -num_alignments 500: format 0 keeps the default hit list of 500 as the tabular
        // formats do (one search, and no warning).
        let same = streams(&words(&["-num_alignments", "500"]), &query, &subject).unwrap();
        for stream in [0, 6, 7, 1] {
            assert_eq!(plain[&stream], same[&stream], "stream {stream}");
        }
        assert!(plain.get(&3).is_none_or(|bytes| bytes.is_empty()));
        assert!(same.get(&3).is_none_or(|bytes| bytes.is_empty()));
        // -num_descriptions 500 is the other default that makes one search; outfmt 6 and 7
        // ignore it with NCBI's warning, which stream 3 carries once.
        let ignored = streams(&words(&["-num_descriptions", "500"]), &query, &subject).unwrap();
        for stream in [0, 6, 7, 1] {
            assert_eq!(plain[&stream], ignored[&stream], "stream {stream}");
        }
        assert_eq!(
            String::from_utf8(ignored[&3].clone()).unwrap(),
            "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n"
        );
    }

    // The CLI's standard error for `-outfmt 6` with `-num_descriptions 3`: the formatting
    // handler's two warnings, then NCBI's error for a hit list of 0 (`-num_alignments 0`).
    const ERR6_ND3_NA0: &str = "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\nWarning: [blastn] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n";
    const SPLIT_ND3_NA0_OUT0: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_nd3_na0.out");
    const SPLIT_ND3_NA0_ERR0: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_nd3_na0.err");

    // -num_alignments 0 gives the tabular formats a hit list of 0, which the CLI rejects with
    // NCBI's error for -outfmt 6 and 7; format 0 is written as the CLI writes it and the run
    // does not fail. Stream 0 is NCBI's frozen bytes, streams 6, 7 and 1 are empty, and stream 3
    // is the CLI's standard error for outfmt 6 followed by that for outfmt 0.
    #[test]
    fn a_hit_list_size_of_zero_fails_the_tabular_formats_only() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        let words = [
            "blastn",
            "-query",
            "tests/fasta/outfmt0/numdesc_split_q.fna",
            "-subject",
            "tests/fasta/outfmt0/numdesc_split_s.fna",
            "-task",
            "blastn",
            "-num_descriptions",
            "3",
            "-num_alignments",
            "0",
        ];
        let streams = streams(&words, &query, &subject).unwrap();
        assert_eq!(text(&streams[&0]), text(SPLIT_ND3_NA0_OUT0));
        for stream in [1, 6, 7] {
            assert!(
                streams.get(&stream).is_none_or(|bytes| bytes.is_empty()),
                "stream {stream}"
            );
        }
        assert_eq!(
            text(&streams[&3]),
            text(&[ERR6_ND3_NA0.as_bytes(), SPLIT_ND3_NA0_ERR0].concat())
        );
        // Both formats' hit lists of 0 fail every format, as the CLI does.
        let both = [
            &words[..7],
            &["-num_descriptions", "0", "-num_alignments", "0"],
        ]
        .concat();
        let error = run(&both, &query, &subject).unwrap_err();
        assert!(error.contains("No hits are being saved"), "{error}");
        assert!(validate(&both)
            .unwrap_err()
            .contains("No hits are being saved"));
        // validate accepts -num_alignments 0: the format 0 search is valid.
        validate(&words).unwrap();
    }

    // The reverse: a -num_descriptions of 2^30 or more overflows BLASTN's preliminary hit list
    // size for format 0 only (the tabular formats ignore -num_descriptions). The CLI's outfmt 6
    // and 7 succeed and its outfmt 0 fails with LOSAT's error (NCBI crashes there); the run
    // writes the tabular formats and the records as the CLI does, stream 0 stays empty, and
    // stream 3 is the CLI's standard error for outfmt 6 followed by that for outfmt 0.
    #[test]
    fn a_hit_list_size_that_overflows_fails_format_0_only() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        let plain_words = ["blastn", "-query", "q", "-subject", "s", "-task", "blastn"];
        let plain = streams(&plain_words, &query, &subject).unwrap();
        let words = [&plain_words[..], &["-num_descriptions", "2147483647"]].concat();
        validate(&words).unwrap();
        let streams = streams(&words, &query, &subject).unwrap();
        // The tabular formats are the plain run's (-num_descriptions is ignored by them).
        for stream in [6, 7] {
            assert_eq!(streams[&stream], plain[&stream], "stream {stream}");
        }
        assert!(streams.get(&0).is_none_or(|bytes| bytes.is_empty()));
        // The records follow the tabular list, with no section of format 0.
        let records = String::from_utf8(streams[&1].clone()).unwrap();
        assert!(records.lines().count() > 0);
        assert!(records.lines().all(|line| line.contains("\"out0\":null")));
        assert_eq!(
            String::from_utf8(streams[&3].clone()).unwrap(),
            "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\nError: a hit list size of 2147483647 (-num_descriptions / -num_alignments), whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTN\n"
        );
        // The same on the other options of the format 0 search: both searches rejected fails.
        let both = [&words[..], &["-num_alignments", "0"]].concat();
        assert!(run(&both, &query, &subject).is_err());
        assert!(validate(&both).is_err());
    }

    // The tolerance of a failed search is the options' checks' rejection of its hit list
    // size and nothing else: any other error (a thread pool that cannot start, which a real
    // run cannot trigger cheaply: the pool is built by the engine from the host's threads)
    // fails the run, as does a failure after a byte was written.
    #[test]
    fn only_a_rejected_hit_list_size_is_tolerated() {
        let failed = |wrote| Failed {
            message: "failed to build blastn pool with 4 threads".into(),
            cli_stderr: Vec::new(),
            wrote,
        };
        assert!(fails_run(&failed(false), false));
        assert!(fails_run(&failed(true), true));
        assert!(fails_run(&failed(true), false));
        assert!(!fails_run(&failed(false), true));
        // The predicate is the options' checks': rejected for a hit list of 0 and for the
        // overflow of N + 50, not for 2^30 or any other size the CLI accepts.
        let command = |extra: &[&str]| {
            let words = [
                &["blastn", "-query", "q", "-subject", "s", "-task", "blastn"][..],
                extra,
            ]
            .concat();
            parse(&words).unwrap().1
        };
        let rejected = |extra: &[&str], group: &[u32]| group_rejected(&command(extra), group);
        assert!(rejected(&["-num_alignments", "0"], &[6, 7]));
        assert!(!rejected(&["-num_alignments", "0"], &[0]));
        assert!(rejected(&["-num_descriptions", "2147483598"], &[0]));
        assert!(!rejected(&["-num_descriptions", "2147483597"], &[0]));
        assert!(!rejected(&["-num_descriptions", "1073741824"], &[0]));
        assert!(!rejected(&[], &[6, 7, 0]));
        let blastp = parse(&["blastp", "-num_alignments", "0"]).unwrap().1;
        assert!(group_rejected(&blastp, &[6]));
        assert!(!group_rejected(&blastp, &[0]));
    }

    // Observer events that break the formatter contract fail the run instead of
    // producing wrong ranges.
    #[test]
    fn unmatched_or_repeated_hsp_events_are_faults() {
        let mut unmatched = recorder();
        unmatched.hsp_end(1, 0);
        assert!(unmatched.fault.is_some());

        let mut crossed = recorder();
        crossed.hsp_begin(1, 0);
        crossed.hsp_end(1, 1);
        assert!(crossed.fault.is_some());

        let mut repeated = recorder();
        for _ in 0..2 {
            repeated.hsp_begin(1, 3);
            repeated.hsp_end(1, 3);
        }
        assert!(repeated.fault.is_some());

        let mut doubled = recorder();
        doubled.hsp_begin(1, 0);
        doubled.hsp_begin(1, 1);
        assert!(doubled.fault.is_some());

        let mut valid = recorder();
        valid.hsp_begin(1, 0);
        valid.hsp_end(1, 0);
        valid.subject_begin(0, 0);
        valid.subject_end(0, 0);
        valid.hsp_begin(0, 0);
        valid.hsp_end(0, 0);
        assert!(valid.fault.is_none());
        assert_eq!(
            (valid.out6.len(), valid.out0.len(), valid.out0_subject.len()),
            (1, 1, 1)
        );
    }
}
