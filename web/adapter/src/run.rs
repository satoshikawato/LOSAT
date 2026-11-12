//! The search: argv parsing (`validate`) and `run` (docs/web/abi_v2.md §4, §7-8).

use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use LOSAT::api::FastaRecord;
use LOSAT::web_adapter::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    Cli, Commands, CountsSource, FormatObserver, FormatOutput, HspIndex, OutputSink, PairwiseHit,
    ReportOutputs,
};

use crate::emit::{emit, StreamWriter, STREAM_CUSTOM, STREAM_DIAGNOSTICS, STREAM_HITS};
use crate::json;

/// The programs that ABI v2 runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    Blastp,
    Tblastn,
    Blastn,
    Tblastx,
    Blastx,
}

impl Program {
    pub fn parse(name: &str) -> Result<Self, String> {
        match name {
            "blastp" => Ok(Self::Blastp),
            "tblastn" => Ok(Self::Tblastn),
            "blastn" => Ok(Self::Blastn),
            "tblastx" => Ok(Self::Tblastx),
            "blastx" => Ok(Self::Blastx),
            "psiblast" => Err("psiblast is not available in the browser (ABI v2)".into()),
            "rpsblast" => Err("rpsblast is not available in the browser (ABI v2)".into()),
            other => Err(format!("unknown program '{other}'")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Blastp => "blastp",
            Self::Tblastn => "tblastn",
            Self::Blastn => "blastn",
            Self::Tblastx => "tblastx",
            Self::Blastx => "blastx",
        }
    }

    /// The output formats that the engine implements for this program; a run writes
    /// exactly these, each on the stream of the same number.
    pub fn formats(self) -> &'static [u32] {
        match self {
            Self::Blastp | Self::Tblastn | Self::Blastn | Self::Tblastx | Self::Blastx => {
                &[0, 6, 7]
            }
        }
    }
}

/// One output of a run: a format that every run of the program writes (`Program::formats`), on
/// the stream of its number, or the output of the argv's own `-outfmt` (any other format that the
/// CLI writes, or 6 and 7 with fields), on `STREAM_CUSTOM` (docs/web/abi_v2.md §5, §7).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Output {
    pub stream: u32,
    /// The format's number, which sets its hit list size.
    pub number: u32,
    /// The `-outfmt` value that the engine writes it for.
    pub outfmt: String,
}

impl Output {
    fn format(number: u32) -> Self {
        Self {
            stream: number,
            number,
            outfmt: number.to_string(),
        }
    }
}

/// The outputs of a run: the output of the argv's `-outfmt` (`parse_run`) if it has one, then
/// the program's formats. The first format of a search writes the warnings that the formats
/// share (the formatting handler's, the queries'), and the output of the `-outfmt` adds those of
/// its fields' writer (NCBI's taxdb warning), so that stream 3 is the CLI's standard error for
/// that `-outfmt` (BLASTX writes the warnings of its first format only).
pub fn outputs(program: Program, custom: Option<&str>) -> Vec<Output> {
    let mut outputs = Vec::new();
    if let Some(outfmt) = custom {
        let choice =
            LOSAT::web_adapter::parse_formatting_string(outfmt).expect("parse_run read the value");
        outputs.push(Output {
            stream: STREAM_CUSTOM,
            number: choice.number as u32,
            outfmt: outfmt.to_string(),
        });
    }
    outputs.extend(program.formats().iter().map(|&f| Output::format(f)));
    outputs
}

/// The `-outfmt` value of the program's arguments.
fn outfmt_of(command: &mut Commands) -> &mut String {
    match command {
        Commands::Blastp(a) => &mut a.outfmt,
        Commands::Tblastn(a) => &mut a.outfmt,
        Commands::Blastn(a) => &mut a.outfmt,
        Commands::Tblastx(a) => &mut a.outfmt,
        Commands::Blastx(a) => &mut a.outfmt,
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {
            unreachable!("rejected by Program::parse")
        }
    }
}

/// The engine's name of the program, for its reading of `-outfmt`.
fn engine_program(program: Program) -> LOSAT::web_adapter::EngineProgram {
    use LOSAT::web_adapter::EngineProgram as Engine;
    match program {
        Program::Blastp => Engine::Blastp,
        Program::Tblastn => Engine::Tblastn,
        Program::Blastn => Engine::Blastn,
        Program::Tblastx => Engine::Tblastx,
        Program::Blastx => Engine::Blastx,
    }
}

/// Parses an argv (program name first) with the engine's CLI parser (`parse_run`), without the
/// output of its `-outfmt`.
pub fn parse(words: &[&str]) -> Result<(Program, Commands), String> {
    parse_run(words).map(|(program, command, _)| (program, command))
}

/// Parses an argv (program name first) with the engine's CLI parser, and reads its `-outfmt` as
/// the CLI does (`parse_formatting_string`, `report_format`, with the CLI's messages; the
/// formatting handler's check of a custom delimiter comes later, `custom_handler_check`). A run
/// writes the program's formats 0, 6 and 7 whatever the value; every other format that the CLI
/// writes (5, 8 to 16, 20, BLASTN's 17), and 6 and 7 with words after the number (fields,
/// `delim=`), add an output of their own, the third value returned (docs/web/abi_v2.md §7). The
/// rule reads the words, not the fields they make, so that the application reads it from the argv
/// alike: `6 std` adds an output with the bytes of outfmt 6. 0 with words adds none (the CLI writes
/// outfmt 0). The output of outfmt 13 and 14 is a ZIP archive of the master file and the side
/// files (`search`). The arguments come back with
/// `-outfmt` at the CLI's default, `0`: the engine is given the formats of each search
/// (`ReportOutputs`), and the checks of the arguments alone (`check_options`) read the formats
/// as they do without `-outfmt`, but for an `-outfmt` whose custom delimiter the CLI rejects
/// (`validate`).
pub fn parse_run(words: &[&str]) -> Result<(Program, Commands, Option<String>), String> {
    for word in words.iter().skip(1) {
        let key = word.split('=').next().unwrap_or("");
        if matches!(key, "-out" | "--out") {
            return Err(format!(
                "{key} is not accepted: a run writes every output format of the program"
            ));
        }
        // The records are registered and read before the options of a run are known
        // (`store::register`), so their ids could not be parsed: the run would write the report
        // of a search without the option.
        if matches!(key, "-parse_deflines" | "--parse_deflines") {
            return Err(format!("{key} is not supported in the browser yet"));
        }
    }
    // The formats come from `Program::formats`, and every program accepts the CLI's
    // default `-outfmt 0`, so the host's words are parsed exactly as on the command line,
    // with the CLI's binary name: an error (an unknown program too) is the CLI's message.
    let argv = std::iter::once("losat").chain(words.iter().copied());
    let cli: Cli = LOSAT::web_adapter::try_parse_cli(argv)
        .map_err(|error| LOSAT::web_adapter::render_message(&error))?;
    let program = Program::parse(words.first().copied().unwrap_or(""))?;
    let mut command = cli.command;
    let outfmt = std::mem::replace(outfmt_of(&mut command), "0".to_string());
    let choice = LOSAT::web_adapter::parse_formatting_string(&outfmt)
        .map_err(|error| format!("{error:#}"))?;
    // The `@PG` line of SAM (outfmt 17) shows the command line, which `main` gives BLASTN: the
    // run's, `losat` and its argv, each word followed by one space (`cli::command_line_bytes`). Only
    // a run that writes SAM holds it.
    if let (Commands::Blastn(args), LOSAT::web_adapter::SAM_FORMAT) = (&mut command, choice.number)
    {
        let mut line = b"losat ".to_vec();
        for word in words {
            line.extend_from_slice(word.as_bytes());
            line.push(b' ');
        }
        args.command_line = Some(line);
    }
    // What the CLI rejects, with its message: a format that LOSAT's program does not write, and the
    // formats that NCBI's formatting handler rejects (`report_format` gives none: SAM of another
    // program, AIRR, FASTA), which a run never gives the engine; the CLI rejects those after the
    // ranges and the query options, so for an argv with an error of both the CLI reports the
    // other. The CLI is taken with `-out`, which the formats 13 and 14 need there (a run gives
    // them a receiver of their side files, `search`). The handler's
    // check of a custom delimiter comes with the program's checks of the options (`validate`) and
    // in the engine's `run_local`.
    let engine = engine_program(program);
    LOSAT::web_adapter::report_format(&choice, engine, Some(std::path::Path::new("out")))
        .and_then(|format| match format {
            Some(_) => Ok(()),
            None => LOSAT::web_adapter::formatting_handler_check(&choice, engine),
        })
        .map_err(|error| format!("{error:#}"))?;
    // NCBI's white space around the value (`NStr::TruncateSpaces`), as `parse_formatting_string`
    // reads it: a space inside is followed by words.
    let words_follow = outfmt
        .trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r'))
        .contains(' ');
    let custom = match choice.number {
        0 => None,
        6 | 7 => words_follow.then_some(outfmt),
        _ => Some(outfmt),
    };
    Ok((program, command, custom))
}

/// Whether the run of `words` writes the BLAST archive (outfmt 11), whose Bioseqs show the whole
/// records (`store::with_inputs`).
pub fn writes_archive(words: &[&str]) -> bool {
    parse_run(words).is_ok_and(|(_, _, custom)| {
        custom.is_some_and(|outfmt| {
            LOSAT::web_adapter::parse_formatting_string(&outfmt)
                .is_ok_and(|choice| choice.number == LOSAT::web_adapter::ARCHIVE_FORMAT)
        })
    })
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
    use LOSAT::web_adapter::{parse_optional_range, RangeRole};
    // The ranges' syntax, in NCBI's order: the subject range where the subjects are read
    // (before the other option handlers), the query range with the query options (the
    // engine's option checks read it). Whether a range fits the records is checked by
    // `run`, which has them.
    let subject_range = |subject_loc: Option<&str>, program: &str| {
        parse_optional_range(subject_loc, RangeRole::Subject, program).map(|_| ())
    };
    let (program, command, custom) = parse_run(words)?;
    let groups = search_groups(&outputs(program, custom.as_deref()), &command);
    let custom_check = || custom_handler_check(program, custom.as_deref());
    // A custom delimiter that a field uses makes the CLI reject the argv, whatever the hit list
    // sizes of the searches: its error is that of the CLI's checks with the argv's `-outfmt`, which
    // make the handler's check where the CLI does (the error may come before it).
    let custom_rejected = custom_check().is_err();
    let with_custom = |outfmt: &mut String| {
        if custom_rejected {
            outfmt.clone_from(custom.as_ref().expect("a rejected -outfmt"));
        }
    };
    // The CLI checks the thread count before the options (`validate_threads`), the checks of the
    // options alone last; with the argv's `-outfmt` they stop at its delimiter first.
    if custom_rejected {
        LOSAT::web_adapter::validate_threads(num_threads(&command))
            .map_err(|error| format!("{error:#}"))?;
    }
    match command {
        Commands::Blastn(mut args) => subject_range(args.subject_loc.as_deref(), "BLASTN")
            .and_then(|()| args.resolve_dust())
            .and_then(|()| {
                if validates_window_masker(&args) {
                    give_blastn_counts(&mut args, &crate::files::source());
                    args.resolve_window_masker()
                } else {
                    Ok(())
                }
            })
            .and_then(|()| {
                parse_optional_range(args.query_loc.as_deref(), RangeRole::Query, "BLASTN")
                    .map(|_| ())
            })
            .map_err(|error| format!("{error:#}"))
            // The formatting handler comes after the query options (`process_options`).
            .and_then(|()| custom_check())
            .and_then(|()| {
                let mut first_error = None;
                for group in &groups {
                    match blastn_group_check(&mut args, group) {
                        Ok(()) => return Ok(()),
                        Err(error) => {
                            first_error.get_or_insert(error);
                        }
                    }
                }
                first_error.map_or(Ok(()), Err)
            })?,
        // The checks of the options read the formats as without `-outfmt` (`parse_run`); with an
        // `-outfmt` whose custom delimiter the CLI rejects they read the argv's (`with_custom`),
        // and fail as the CLI's.
        Commands::Tblastx(mut args) => subject_range(args.subject_loc.as_deref(), "TBLASTX")
            .and_then(|()| {
                with_custom(&mut args.outfmt);
                LOSAT::web_adapter::check_tblastx_options(&args)
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| custom_check())?,
        Commands::Blastp(mut args) => subject_range(args.subject_loc.as_deref(), "BLASTP")
            .and_then(|()| {
                with_custom(&mut args.outfmt);
                LOSAT::web_adapter::check_blastp_options(&args)
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| custom_check())?,
        Commands::Tblastn(mut args) => subject_range(args.subject_loc.as_deref(), "TBLASTN")
            .and_then(|()| {
                with_custom(&mut args.outfmt);
                LOSAT::web_adapter::check_tblastn_options(&args)
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| custom_check())?,
        // The checks of each search's formats, as for BLASTN: BLASTX resolves its options for
        // each format (the hit list size). Its ranges are read among them, where its run reads
        // them: after LOSAT's limits of the options (`blastx::check_options`).
        Commands::Blastx(args) => {
            // The tabular search's formats hold the argv's `-outfmt`: its error is the run's.
            if custom_rejected {
                return blastx_group_check(&args, &groups[0]).and_then(|()| custom_check());
            }
            let mut first_error = None;
            for group in &groups {
                match blastx_group_check(&args, group) {
                    Ok(()) => return Ok(()),
                    Err(error) => {
                        first_error.get_or_insert(error);
                    }
                }
            }
            if let Some(error) = first_error {
                return Err(error);
            }
        }
        _ => {}
    }
    Ok(())
}

/// NCBI's formatting handler's check of the argv's own `-outfmt` (a custom delimiter that a field
/// uses to separate its entries), with the CLI's message. `validate` and `run` ask it whether the
/// CLI rejects the argv: `validate` then reports the error of the program's checks of the options
/// with that `-outfmt`, which make it where the CLI does (BLASTN's right after the query range),
/// and `run` tolerates no failed search. A run's engine makes it (`run_local`).
fn custom_handler_check(program: Program, custom: Option<&str>) -> Result<(), String> {
    let Some(outfmt) = custom else {
        return Ok(());
    };
    LOSAT::web_adapter::parse_formatting_string(outfmt)
        .and_then(|choice| {
            LOSAT::web_adapter::formatting_handler_check(&choice, engine_program(program))
        })
        .map_err(|error| format!("{error:#}"))
}

/// BLASTN's checks of the options for the search of the output formats `group`, whose hit list
/// size is the one the CLI derives for them: the checks that `run_local` makes before it reads
/// a query (`process_options`: the counts of the report, `check_scoring_options`; then
/// `check_losat_limits`). `validate` runs it for each search of a run, and `run` asks it which
/// search is rejected (`group_rejected`).
///
/// It does not check the Karlin-Altschul tables (scores and gap costs without a row, NCBI's
/// `Error: Substitution scores ...`) nor LOSAT's limit of greedy gap costs: a CLI run meets them
/// only when some query context has an ungapped Karlin-Altschul block (`context_blocks`,
/// `any_valid_context` in `blast_engine/run.rs`), which depends on the query's letters (a query
/// of N only has none, and the same scores give rc 0 with a warning). `validate` has no query,
/// so `run` makes these checks and stops with the CLI's message and exit code.
/// Making them for a query of any letters would reject what the CLI accepts (WA-FM-01).
fn blastn_group_check(
    args: &mut LOSAT::web_adapter::BlastnArgs,
    group: &[Output],
) -> Result<(), String> {
    use LOSAT::web_adapter::{check_losat_limits, check_scoring_options};
    // The host has no database files.
    if args.db.is_some() {
        return Err("-db (a database search) is not supported in the browser yet".to_string());
    }
    let outfmts: Vec<i32> = group.iter().map(|output| output.number as i32).collect();
    // `process_options` sets the counts that the checks read.
    let (counts, _) = LOSAT::web_adapter::run_report_counts(
        "blastn",
        &outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
    )
    .map_err(|error| format!("{error:#}"))?;
    args.report_counts = counts;
    check_scoring_options(args)
        .and_then(|()| check_losat_limits(args))
        .map_err(|error| format!("{error:#}"))
}

/// BLASTX's checks (`blastx::check_options`) for the search of the output formats `group`, as
/// `run_local` makes them before it searches. `validate` runs it for each search of a run, and
/// `run` asks it which search is rejected (`group_rejected`).
fn blastx_group_check(
    args: &LOSAT::web_adapter::BlastxArgs,
    group: &[Output],
) -> Result<(), String> {
    let outfmts: Vec<&str> = group.iter().map(|output| output.outfmt.as_str()).collect();
    LOSAT::web_adapter::check_blastx_options(args, &outfmts).map_err(|error| format!("{error:#}"))
}

/// Whether the options' checks reject the search of the output formats `group`, as the CLI
/// rejects those formats alone because of their hit list size: BLASTN's and BLASTX's checks
/// fail (`blastn_group_check`: a hit list of 0, or a preliminary hit list size below 1;
/// `blastx_group_check`), or, for the other programs, the hit list size is 0. `run` tolerates
/// a failed search of a two-search run only for this reason (`fails_run`). With the checks
/// themselves as the predicate, a search that fails them for another reason than its hit list
/// (for example `-evalue 0`) is rejected as the other search is, so `run` fails with the
/// tabular search's error, which `validate` reports.
fn group_rejected(command: &mut Commands, group: &[Output]) -> bool {
    match command {
        Commands::Blastn(args) => blastn_group_check(args, group).is_err(),
        Commands::Blastx(args) => blastx_group_check(args, group).is_err(),
        _ => hit_list_size(command, group[0].number) == 0,
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
    use LOSAT::web_adapter::{parse_optional_range, record_interval, RangeRole, RecordInterval};
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
        // A positive `i32` (the parser's `positive`).
        Commands::Blastx(a) => [
            a.max_target_seqs.map(|n| n as usize),
            a.num_descriptions,
            a.num_alignments,
        ],
        // Rejected by `Program::parse`.
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => [None; 3],
    }
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
        // Rejected by `Program::parse`.
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => 1,
    }
}

/// The hit list size that the CLI derives for the output format `format` from the options
/// of `command` (`report_counts`, NCBI `blast_args.cpp:2894-2973`).
pub fn hit_list_size(command: &Commands, format: u32) -> usize {
    let [max_target_seqs, descriptions, alignments] = count_options(command);
    LOSAT::web_adapter::report_counts(format as i32, max_target_seqs, descriptions, alignments)
        .hitlist_size
}

/// The searches of a run: the output formats grouped by hit list size, the tabular
/// formats first. One group is one `run_local` call, and every format of it is written from
/// that search, as the CLI writes it alone: the CLI's search depends on the hit list size,
/// which `-num_descriptions` and `-num_alignments` (without `-max_target_seqs`) can make
/// differ between the pairwise format 0 and the tabular formats. Without those options the
/// sizes are equal and the group is the program's formats with the tabular ones first (the one
/// search of before). The tabular formats go first because the formatting handler's warnings
/// belong to the first format, and the CLI writes the tabular formats' (`ReportCounts::warnings`;
/// `-line_length`, `-sorthits` and `-sorthsps` warn for formats above 4 or other than 0 only,
/// blast_args.cpp:2947-2989). The
/// output of the argv's `-outfmt` (a format above 4) has the hit list size of the tabular
/// formats, so it joins the tabular group, first (`outputs`): it adds no search.
pub fn search_groups(formats: &[Output], command: &Commands) -> Vec<Vec<Output>> {
    let [_, descriptions, alignments] = count_options(command);
    let mut ordered: Vec<Output> = formats.to_vec();
    ordered.sort_by_key(|format| format.number <= 4);
    if descriptions.is_none() && alignments.is_none() {
        return vec![ordered];
    }
    let mut groups: Vec<(usize, Vec<Output>)> = Vec::new();
    for format in ordered {
        let size = hit_list_size(command, format.number);
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
        let cli_stderr = match error.downcast_ref::<LOSAT::web_adapter::NativeError>() {
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

/// The name that the side files of outfmt 13 and 14 follow, as the CLI's `-out` (no directory and
/// no extension, so that the side files have the names that the master file lists: `results_1.json`
/// for 13 too, where an `-out` ending in `.json` lists other names than its side files have,
/// `report::blast_xml2::files`).
pub const SIDE_FILE_BASE: &str = "results";

/// The name of the master file in the archive of outfmt 13 or 14.
pub fn master_name(number: u32) -> String {
    let extension = if number == 13 { "json" } else { "xml" };
    format!("{SIDE_FILE_BASE}.{extension}")
}

/// The output of outfmt 13 or 14 (stream 4): a ZIP archive (`zip::ZipStream`) of the side files,
/// added as the engine hands each over, and of the master file, added after the run.
struct ArchivedFiles(std::sync::Mutex<Option<crate::zip::ZipStream<StreamWriter>>>);

impl ArchivedFiles {
    fn take(&self) -> crate::zip::ZipStream<StreamWriter> {
        self.0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .take()
            .expect("the archive is taken once, after the run")
    }
}

impl LOSAT::web_adapter::SideFileReceiver for ArchivedFiles {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
        match self
            .0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .as_mut()
        {
            Some(zip) => zip.add(name, bytes),
            None => Err(std::io::Error::other("the archive was closed")),
        }
    }
}

/// What one search wrote.
struct Searched {
    writers: Vec<StreamWriter>,
    recorder: RangeRecorder,
    hits: Option<std::sync::Arc<Vec<PairwiseHit>>>,
    ordinals: Option<Vec<usize>>,
}

/// One `run_local` call that writes the output formats `formats` (a group of
/// `search_groups`), with the warnings of the CLI to `diagnostics`.
fn search(
    command: Commands,
    formats: &[Output],
    query: &[FastaRecord],
    subject: &[FastaRecord],
    diagnostics: &mut (dyn std::io::Write + Send),
) -> Result<Searched, Failed> {
    let mut writers: Vec<StreamWriter> = formats
        .iter()
        .map(|format| StreamWriter::new(format.stream))
        .collect();
    let mut recorder = RangeRecorder {
        streams: formats.iter().map(|format| format.stream).collect(),
        positions: writers.iter().map(StreamWriter::position).collect(),
        open: vec![None; formats.len()],
        out6: Vec::new(),
        out0: Vec::new(),
        out0_subject: Vec::new(),
        subject_start: None,
        subject: None,
        fault: None,
    };
    let mut collected: Option<std::sync::Arc<Vec<PairwiseHit>>> = None;
    // Outfmt 13 and 14: the master file goes to memory and the side files to a ZIP archive on the
    // format's stream, to which the master file is added after the run (`ArchivedFiles`).
    let mut archives: Vec<Option<(Vec<u8>, Arc<ArchivedFiles>)>> = formats
        .iter()
        .zip(writers.iter_mut())
        .map(|(format, writer)| {
            matches!(format.number, 13 | 14).then(|| {
                let stream = std::mem::replace(writer, StreamWriter::new(format.stream));
                let zip = crate::zip::ZipStream::new(stream);
                (
                    Vec::new(),
                    Arc::new(ArchivedFiles(std::sync::Mutex::new(Some(zip)))),
                )
            })
        })
        .collect();
    let ordinals;
    let outcome = {
        let mut on_hits = |hits| collected = Some(hits);
        let mut outputs = ReportOutputs {
            formats: formats
                .iter()
                .zip(writers.iter_mut())
                .zip(archives.iter_mut())
                .map(|((format, writer), archive)| FormatOutput {
                    outfmt: &format.outfmt,
                    sink: match archive {
                        Some((master, files)) => OutputSink::Files(
                            master,
                            LOSAT::web_adapter::SideFileSink::new(SIDE_FILE_BASE, files.clone()),
                        ),
                        None => OutputSink::Writer(writer),
                    },
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
            Commands::Blastx(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), "BLASTX", query);
                run_local_blastx(args, query, subject, &mut outputs)
            }
            Commands::Psiblast(_)
            | Commands::Rpsblast(_)
            | Commands::Makedb(_)
            | Commands::Rpsbproc(_) => {
                unreachable!("rejected by Program::parse")
            }
        };
        result
    };
    // The archives: the master file after the side files and the central directory, on the
    // format's stream again. A run that failed leaves the side files written (the host discards
    // them).
    for ((format, writer), archive) in formats.iter().zip(writers.iter_mut()).zip(&mut archives) {
        let Some((master, files)) = archive.take() else {
            continue;
        };
        let zip = files.take();
        *writer = if outcome.is_ok() {
            let mut zip = zip;
            match zip
                .add(master_name(format.number).as_bytes(), &master)
                .and_then(|()| zip.finish())
            {
                Ok(stream) => stream,
                Err(error) => {
                    return Err(Failed {
                        cli_stderr: format!("{error}\n").into_bytes(),
                        message: format!("the archive of outfmt {}: {error}", format.number),
                        wrote: true,
                    })
                }
            }
        } else {
            zip.into_inner()
        };
    }
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

/// Whether `validate` reads the WindowMasker options where the command line reads them: after
/// `-dust`, before the query range and the other options (`resolve_window_masker`: the check of
/// the file of `-window_masker_db`, with the CLI's message for a format other than the two binary
/// ones, then "at most one of" the options, then the file of `-window_masker_taxid`). Not when
/// `-window_masker_db` names a file that is not registered: a host that validates an argv without
/// its files (the application's control reactor) leaves the file to `run`, where a name that is
/// not registered is the CLI's file that cannot be opened. A taxid whose file is not registered
/// is the CLI's taxid without a file, whose error comes with the search (`batch_error`). An argv
/// without the options is validated as before named input files.
fn validates_window_masker(args: &LOSAT::web_adapter::BlastnArgs) -> bool {
    match args.window_masker_db.as_deref() {
        Some(name) => crate::files::registered(name),
        None => args.window_masker_taxid.is_some(),
    }
}

/// The counts source of BLASTN's `-window_masker_db` and `-window_masker_taxid`: the named input
/// files registered when the run starts (`files::source`), which the engine opens by the option's
/// value, or by the taxid's name in `files::WINDOW_MASKER_PATH`, as the command line opens the
/// file of that name. `None` without the options, so a run without them gives the engine what it
/// gave before named input files.
fn window_masker_counts(command: &Commands) -> Option<Arc<dyn CountsSource>> {
    match command {
        Commands::Blastn(args)
            if args.window_masker_db.is_some() || args.window_masker_taxid.is_some() =>
        {
            Some(crate::files::source())
        }
        _ => None,
    }
}

/// Gives a search's command the run's counts source (`window_masker_counts`).
fn give_counts(command: &mut Commands, counts: Option<&Arc<dyn CountsSource>>) {
    if let (Commands::Blastn(args), Some(counts)) = (command, counts) {
        give_blastn_counts(args, counts);
    }
}

/// Gives BLASTN's arguments the counts source, and with `-window_masker_taxid` LOSAT Web's
/// registry of counts files, where the command line gives NCBI's registry: the directory
/// `files::WINDOW_MASKER_PATH`, and no posts of a `.ncbirc` (LOSAT Web reads none).
fn give_blastn_counts(args: &mut LOSAT::web_adapter::BlastnArgs, counts: &Arc<dyn CountsSource>) {
    args.window_masker_counts = Some(counts.clone());
    if args.window_masker_taxid.is_some() {
        args.window_masker_registry = Some(Ok(Default::default()));
        if let Some(Ok(registry)) = &mut args.window_masker_registry {
            registry.path = Some(crate::files::WINDOW_MASKER_PATH.to_string());
        }
    }
}

/// Runs the searches of the run and emits its streams. On failure nothing more is emitted
/// and the host discards what it received.
///
/// Without `-num_descriptions` and `-num_alignments` this is one search that writes every
/// format. With them it is one search per distinct hit list size of the formats
/// (`search_groups`), so that each format's stream is the CLI's output for that format
/// alone. The output of the argv's `-outfmt` (stream 4) is written by the tabular search, as the
/// CLI writes that `-outfmt`. The HSP records follow the tabular search;
/// a record gets the section of format 0 only where the search of format 0 has exactly its HSPs
/// (`link_format0`) (docs/web/abi_v2.md §5, §7, §8).
///
/// A search of a two-search run that fails before it writes a byte of its formats, and
/// that the options' checks reject because of its hit list size (`group_rejected`:
/// `-num_alignments 0` for the tabular formats, a BLASTN `-num_descriptions` of 2,147,483,598
/// or more for format 0, where N + 50 overflows the 32-bit preliminary hit list size), does
/// not fail the run when the other search succeeds, as the CLI fails that format alone: the
/// streams of its formats stay empty (with no tabular search there are no records). Any other
/// error of a search (for example a thread pool that cannot start) fails the run, as does a
/// failure after a byte was written; so do both searches failing, as the CLI fails every
/// format.
///
/// Stream 3 is the CLI's standard error for the tabular formats (the tabular search's
/// diagnostics, with its error when it failed): for `-outfmt 6`, or for the argv's `-outfmt`
/// when it has an output of its own (`outputs`). When a search failed it is followed by the
/// CLI's standard error for format 0 (the diagnostics of the format 0 search, with its
/// error when it failed). The diagnostics of the format 0 search of a run without a failure
/// repeat the first search's and are dropped.
pub fn run(words: &[&str], query: &[FastaRecord], subject: &[FastaRecord]) -> Result<(), String> {
    use std::io::Write as _;
    let (program, command, custom) = parse_run(words)?;
    let groups = search_groups(&outputs(program, custom.as_deref()), &command);
    let custom_accepted = custom_handler_check(program, custom.as_deref()).is_ok();
    let mut diagnostics = StreamWriter::new(STREAM_DIAGNOSTICS);
    let threads = num_threads(&command);
    let counts = window_masker_counts(&command);
    let mut command = Some(command);
    // Each search: what it wrote or why it failed, and (after the first) its diagnostics.
    let mut outcomes: Vec<(Result<Searched, Failed>, Vec<u8>)> = Vec::new();
    // The searches of the run share one thread pool (`with_run_pool`): a browser holds a thread
    // worker for every thread that may start, and the threads of a pool may not have returned when
    // the next pool starts, so each pool of its own needed thread workers of its own.
    LOSAT::web_adapter::with_run_pool(threads, program.name(), || -> Result<(), String> {
        for group in &groups {
            let mut command = match command.take() {
                Some(command) => command,
                None => parse(words)?.1,
            };
            give_counts(&mut command, counts.as_ref());
            let mut scratch = Vec::new();
            // Whether the options' checks reject this search, found before `search` takes the
            // command.
            let rejected = group_rejected(&mut command, group);
            // The warnings of the first search include the reader's and the queries' warnings
            // that the next searches repeat; the CLI writes them once.
            let result = if outcomes.is_empty() {
                search(command, group, query, subject, &mut diagnostics)
            } else {
                search(command, group, query, subject, &mut scratch)
            };
            // The CLI rejects an `-outfmt` whose custom delimiter a field uses: the failure of its
            // search (the tabular one) is never tolerated.
            if let Err(failed) = &result {
                if fails_run(failed, rejected && custom_accepted) {
                    return Err(failed.message.clone());
                }
            }
            outcomes.push((result, scratch));
        }
        Ok(())
    })?;
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
        for owned in ["-out", "-out=x"] {
            let error = parse(&["blastp", "-query", "q", "-subject", "s", owned, "6"]).unwrap_err();
            assert!(error.contains("is not accepted"), "{owned}: {error}");
        }
        // -num_descriptions and -num_alignments are accepted (also as `-x=v`),
        // parsed as on the command line; -max_target_seqs excludes both (clap's conflict).
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
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
                let cli = LOSAT::web_adapter::try_parse_cli(
                    std::iter::once("losat").chain(both.iter().copied()),
                )
                .unwrap_err();
                assert_eq!(
                    parse(&both).unwrap_err(),
                    LOSAT::web_adapter::render_message(&cli),
                    "{program} {option}"
                );
                assert!(parse(&both).unwrap_err().contains("cannot be used with"));
            }
        }
        assert!(matches!(
            parse(&["blastx", "-query", "q", "-subject", "s"]),
            Ok((Program::Blastx, Commands::Blastx(_)))
        ));
        assert!(parse(&["psiblast", "-query", "q", "-subject", "s"]).is_err());
        // An unknown program and a parser error are the CLI's messages.
        for words in [&["nosuch"][..], &["tblastx", "-query", "q"][..]] {
            let cli = LOSAT::web_adapter::try_parse_cli(
                std::iter::once("losat").chain(words.iter().copied()),
            )
            .map(|_| ())
            .map_err(|error| LOSAT::web_adapter::render_message(&error));
            assert_eq!(parse(words).map(|_| ()), cli, "{words:?}");
        }
        assert_eq!(Program::Tblastx.formats(), &[0, 6, 7]);
        assert_eq!(Program::Blastx.formats(), &[0, 6, 7]);
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
        // The tables depend on the query's letters, so `run` checks them (WA-FM-01,
        // `blastn_scores_without_table_row_are_left_to_the_run`).
        assert!(validate(&unsupported).is_ok());
        let letters = b"ACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA";
        let error = run(
            &unsupported,
            &[FastaRecord::new("Query_1", b"q", letters)],
            &[FastaRecord::new("Subject_1", b"s", letters)],
        )
        .unwrap_err();
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
        // validate reports BLASTX's checks of each search's options, as the run does (NCBI's
        // errors and LOSAT's refusals); the tabular search that -num_alignments 0 rejects alone
        // is tolerated, as `run` tolerates it.
        let blastx = ["blastx", "-query", "q", "-subject", "s"];
        assert!(validate(&blastx).is_ok());
        assert!(validate(&[&blastx[..], &["-num_alignments", "0"]].concat()).is_ok());
        for (extra, start) in [
            (
                &["-evalue", "0"][..],
                "BLAST query/options error: expect value",
            ),
            (
                &["-matrix", "IDENTITY"],
                "BLAST query/options error: IDENTITY is not a supported matrix",
            ),
            (
                &["-num_descriptions", "0", "-num_alignments", "0"],
                "BLAST query/options error: No hits are being saved",
            ),
        ] {
            let error = validate(&[&blastx[..], extra].concat()).unwrap_err();
            assert!(error.starts_with(start), "{extra:?}: {error}");
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
            let cli = LOSAT::web_adapter::try_parse_cli(
                ["losat"].into_iter().chain(words.iter().copied()),
            )
            .unwrap_err();
            assert_eq!(
                parse(words).unwrap_err(),
                LOSAT::web_adapter::render_message(&cli),
                "{words:?}"
            );
        }
    }

    /// The text a run wrote on `stream`.
    fn emitted(stream: u32) -> String {
        let bytes = crate::emit::EMITTED.with(|emitted| {
            emitted
                .borrow()
                .iter()
                .filter(|(number, _)| *number == stream)
                .flat_map(|(_, bytes)| bytes.clone())
                .collect::<Vec<u8>>()
        });
        String::from_utf8(bytes).unwrap()
    }

    // WA-FM-01. The CLI meets the Karlin-Altschul tables (scores and gap costs without a row,
    // LOSAT's limit of greedy gap costs) only when some query context has an ungapped
    // Karlin-Altschul block (`context_blocks`), which depends on the query's letters: a query
    // of N only has none, and the same scores give a warning and rc 0. `validate` has no query,
    // so it makes the checks of the options only, and `run` rejects with the CLI's message and
    // exit code (`validate` accepts a subset of what the CLI accepts).
    #[test]
    fn blastn_scores_without_table_row_are_left_to_the_run() {
        let megablast = ["blastn", "-query", "q", "-subject", "s"];
        let blastn = [&megablast[..], &["-task", "blastn"]].concat();
        let with = |base: &[&str], extra: &[&str]| validate(&[base, extra].concat());
        // The scores of the symptom: no table row, and a pair that the CLI rejects for some
        // queries.
        for base in [&megablast[..], &blastn[..]] {
            for pair in [
                &["-reward", "7"][..],
                &["-reward", "7", "-penalty", "-2"],
                &["-reward", "3", "-penalty", "-1"],
                &["-reward", "4", "-penalty", "-1"],
                &["-reward", "5", "-penalty", "-1"],
                &["-reward", "2", "-penalty", "-1"],
            ] {
                assert_eq!(with(base, pair), Ok(()), "{base:?} {pair:?}");
            }
        }
        // The checks of the options stand: NCBI's and LOSAT's refusals.
        let error = with(&megablast, &["-reward", "7", "-evalue", "0"]).unwrap_err();
        assert!(
            error.starts_with("BLAST query/options error: expect value"),
            "{error}"
        );
        let gaps = [
            "-reward",
            "3",
            "-penalty",
            "-1",
            "-gapopen",
            "5",
            "-gapextend",
            "0",
        ];
        let error = with(&blastn, &gaps).unwrap_err();
        assert!(
            error.contains("gap extension penalty cannot be 0"),
            "{error}"
        );
        let error = with(&megablast, &["-reward", "0"]).unwrap_err();
        assert!(error.contains("a reward of 0"), "{error}");

        // `run` meets the tables with the query's letters, as the CLI does.
        let search = |base: &[&str], extra: &[&str], letters: &[u8]| {
            crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
            let query = [FastaRecord::new("Query_1", b"mq1 query", letters)];
            let subject = [FastaRecord::new("Subject_1", b"ms1", ACGT)];
            let result = run(&[base, extra].concat(), &query, &subject);
            (
                result,
                emitted(crate::emit::STREAM_DIAGNOSTICS),
                emitted(crate::emit::STREAM_HITS),
            )
        };
        const ACGT: &[u8] = b"ACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA";
        const NNNN: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        const WARNING: &str =
            "Warning: [blastn] Query_1 mq1 query: Could not calculate ungapped Karlin-Altschul";
        // An expected score that is not negative (no block for ACGT): a warning, no hit.
        for base in [&megablast[..], &blastn[..]] {
            for pair in [
                &["-reward", "3", "-penalty", "-1"][..],
                &["-reward", "5", "-penalty", "-1"],
            ] {
                let (result, diagnostics, hits) = search(base, pair, ACGT);
                assert_eq!(result, Ok(()), "{base:?} {pair:?}");
                assert!(diagnostics.starts_with(WARNING), "{diagnostics}");
                assert_eq!(hits, "");
            }
        }
        let (result, diagnostics, _) = search(&megablast, &["-reward", "7"], ACGT);
        assert_eq!(result, Ok(()));
        assert!(diagnostics.starts_with(WARNING), "{diagnostics}");
        // A query of N only has no block for any scores: `-task blastn -reward 7` is a warning.
        let (result, diagnostics, hits) = search(&blastn, &["-reward", "7"], NNNN);
        assert_eq!(result, Ok(()));
        assert!(diagnostics.starts_with(WARNING), "{diagnostics}");
        assert_eq!(hits, "");
        // A query of A, C, G and T has a block for 7 and -3, which has no row: the CLI's error.
        let (result, _, _) = search(&blastn, &["-reward", "7"], ACGT);
        let error = result.unwrap_err();
        assert_eq!(
            error.trim_end(),
            "BLAST engine error: Error: Substitution scores 7 and -3 are not supported"
        );
        // Letters of ambiguity put a block under scores whose expected score is 0: (3, -1) with
        // a B (LOSAT's own limit of the gap costs, or NCBI's table error, as the CLI).
        let mut ambiguous = ACGT.to_vec();
        ambiguous.extend_from_slice(b"B");
        let (result, _, _) = search(&megablast, &["-reward", "3", "-penalty", "-1"], &ambiguous);
        assert!(result.is_err());
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
    // The CLI writes these records' ids from the option; a run here has read them without it.
    #[test]
    fn parse_deflines_is_rejected() {
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            let words = [program, "-query", "q", "-subject", "s", "-parse_deflines"];
            let expected = "-parse_deflines is not supported in the browser yet";
            assert_eq!(validate(&words).unwrap_err(), expected);
            assert_eq!(run(&words, &[], &[]).unwrap_err(), expected);
        }
    }

    // `run`, before any search.
    #[test]
    fn a_database_search_is_rejected() {
        let words = ["blastp", "-query", "q", "-db", "swissprot"];
        let error = validate(&words).unwrap_err();
        assert_eq!(
            error,
            "-db (a database search) is not supported in the browser yet"
        );
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [FastaRecord::new("Query_1", b"q", protein)];
        let subject = [FastaRecord::new("Subject_1", b"s", protein)];
        assert_eq!(
            run(&words, &query, &subject).unwrap_err(),
            "-db (a database search) is not supported in the browser yet"
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

    /// The program's formats `numbers` (`Program::formats`), as a search group holds them.
    fn formats(numbers: &[u32]) -> Vec<Output> {
        numbers
            .iter()
            .map(|&number| Output::format(number))
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
    // -num_alignments the one search of the program's formats, the tabular ones first.
    #[test]
    fn a_run_searches_once_for_each_hit_list_size() {
        let groups = |program: &str, options: &[&str]| {
            let words = [&[program, "-query", "q", "-subject", "s"][..], options].concat();
            let (program, command, custom) = parse_run(&words).unwrap();
            search_groups(&outputs(program, custom.as_deref()), &command)
                .iter()
                .map(|group| group.iter().map(|output| output.stream).collect())
                .collect::<Vec<Vec<u32>>>()
        };
        let one = vec![vec![6, 7, 0]];
        let one_tabular_first = vec![vec![6, 7, 0]];
        let two = vec![vec![6, 7], vec![0]];
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
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
            // The output of the argv's -outfmt (stream 4) has the hit list size of 6 and 7: it
            // joins the tabular search, first, and adds none.
            for (options, expected) in [
                (&["-outfmt", "6 qseqid"][..], vec![vec![4, 6, 7, 0]]),
                (
                    &["-outfmt", "7 std qcovs", "-num_alignments", "3"],
                    vec![vec![4, 6, 7], vec![0]],
                ),
                (
                    &["-outfmt=6 sseqid", "-num_descriptions", "500"],
                    vec![vec![4, 6, 7, 0]],
                ),
                (&["-outfmt", "7"], vec![vec![6, 7, 0]]),
            ] {
                assert_eq!(groups(program, options), expected, "{program} {options:?}");
            }
        }
    }

    // -outfmt is read as the CLI reads it, with the CLI's messages; 0, 6 and 7 add no output, 6 and
    // 7 with words after the number add one (WA2-OP-03, Owner 2026-11-05), and so does every other
    // format that the CLI writes (W11), 13 and 14 without the CLI's need of -out.
    #[test]
    fn outfmt_is_read_as_on_the_command_line() {
        let custom = |program: &str, outfmt: &str| {
            parse_run(&[program, "-query", "q", "-subject", "s", "-outfmt", outfmt])
                .map(|(_, _, custom)| custom)
        };
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            for outfmt in ["0", "6", "7", " 7 ", "0 qseqid", "+6"] {
                assert_eq!(custom(program, outfmt), Ok(None), "{program} {outfmt:?}");
            }
            for outfmt in [
                "6 qseqid sseqid pident",
                "7 std qcovs",
                "6 std",
                "6 delim=,",
                " 6 -qseqid ",
                "+7 sseqid",
            ] {
                assert_eq!(
                    custom(program, outfmt),
                    Ok(Some(outfmt.to_string())),
                    "{program} {outfmt:?}"
                );
            }
            // The -flag=value form, and the arguments come back with the CLI's default.
            let words = [program, "-query", "q", "-subject", "s", "-outfmt=6 qseqid"];
            let (_, mut command, given) = parse_run(&words).unwrap();
            assert_eq!(given.as_deref(), Some("6 qseqid"));
            assert_eq!(outfmt_of(&mut command), "0");
            // Every other format that the CLI writes adds an output, with or without words (W11).
            for outfmt in [
                "5",
                "8",
                "9",
                "10",
                "12",
                "15",
                "16",
                "20",
                " 5 ",
                "+010",
                "10 qseqid",
                "20 delim=| std",
                "5 qseqid",
                "11",
                "13",
                "14",
                "14 qseqid",
            ] {
                assert_eq!(
                    custom(program, outfmt),
                    Ok(Some(outfmt.to_string())),
                    "{program} {outfmt:?}"
                );
                let words = [program, "-query", "q", "-subject", "s", "-outfmt", outfmt];
                assert_eq!(validate(&words), Ok(()), "{program} {outfmt:?}");
            }
            // What the CLI rejects, with its message.
            let engine = engine_program(Program::parse(program).unwrap());
            for outfmt in ["abc", "99", "1", "6\tqseqid"] {
                let cli = LOSAT::web_adapter::parse_formatting_string(outfmt)
                    .and_then(|choice| {
                        LOSAT::web_adapter::report_format(&choice, engine, None)?;
                        LOSAT::web_adapter::formatting_handler_check(&choice, engine)
                    })
                    .unwrap_err();
                assert_eq!(
                    custom(program, outfmt),
                    Err(format!("{cli:#}")),
                    "{program} {outfmt:?}"
                );
            }
            // The formatting handler's check of a custom delimiter comes where the CLI makes it:
            // after the ranges (and BLASTN's query options), in validate and in the engine's run.
            let conflict = "6 delim=; qseqid staxids";
            assert_eq!(custom(program, conflict), Ok(Some(conflict.to_string())));
            let words = [program, "-query", "q", "-subject", "s", "-outfmt", conflict];
            let handler = LOSAT::web_adapter::parse_formatting_string(conflict)
                .and_then(|choice| LOSAT::web_adapter::formatting_handler_check(&choice, engine))
                .unwrap_err();
            assert_eq!(validate(&words), Err(format!("{handler:#}")), "{program}");
            assert!(validate(&words)
                .unwrap_err()
                .contains("Your custom record separator (;)"));
            let ranged = [&words[..], &["-subject_loc", "5-2"]].concat();
            assert!(
                validate(&ranged)
                    .unwrap_err()
                    .starts_with("BLAST engine error: Invalid specification of subject location"),
                "{program}: {:?}",
                validate(&ranged)
            );
            let query = [FastaRecord::new("q", b"q", b"ACGTACGTACGT")];
            assert_eq!(
                run(&words, &query, &query),
                Err(format!("{handler:#}")),
                "{program}"
            );
            // The CLI rejects the argv whatever the hit list sizes: a run of two searches does not
            // write the other one, and validate gives the error of the CLI's checks with the argv's
            // -outfmt (review 2).
            for extra in [
                &["-num_descriptions", "3"][..],
                &["-num_alignments", "3"],
                &["-num_alignments", "0"],
                &["-evalue", "0"],
            ] {
                let both = [&words[..], extra].concat();
                assert_eq!(
                    validate(&both),
                    Err(format!("{handler:#}")),
                    "{program} {extra:?}"
                );
                assert_eq!(
                    run(&both, &query, &query),
                    Err(format!("{handler:#}")),
                    "{program} {extra:?}"
                );
            }
            let twice = [
                program, "-query", "q", "-subject", "s", "-outfmt", "6", "-outfmt", "7",
            ];
            assert!(parse(&twice)
                .unwrap_err()
                .contains("cannot be used multiple times"));
            assert!(validate(&[
                program,
                "-query",
                "q",
                "-subject",
                "s",
                "-outfmt",
                "7 std qcovs"
            ])
            .is_ok());
        }
        assert_eq!(
            custom("blastn", "abc").unwrap_err(),
            "BLAST query/options error: 'abc' is not a valid output format\nPlease refer to the BLAST+ user manual.\n"
        );
        assert_eq!(
            custom("blastn", "99").unwrap_err(),
            "Error: Formatting choice is out of range\n"
        );
        assert_eq!(custom("blastp", "17").unwrap_err(), "BLAST query/options error: SAM format is only applicable to blastn\nPlease refer to the BLAST+ user manual.\n");
        assert_eq!(custom("blastn", "17 SQ"), Ok(Some("17 SQ".to_string())));
    }

    /// The streams of `run` over inputs that the program's reader reads (`store::register`).
    fn registered_streams(
        words: &[&str],
        query: &[u8],
        subject: &[u8],
    ) -> std::collections::BTreeMap<u32, Vec<u8>> {
        let program = Program::parse(words[0]).unwrap();
        let (query, _) = crate::store::register(words[0], 0, query).unwrap();
        let (subject, _) = crate::store::register(words[0], 1, subject).unwrap();
        crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
        let result = crate::store::with_inputs(
            program,
            query,
            subject,
            writes_archive(words),
            |queries, subjects| run(words, queries, subjects),
        );
        crate::store::release(query).unwrap();
        crate::store::release(subject).unwrap();
        result.unwrap().unwrap();
        let mut streams = std::collections::BTreeMap::new();
        crate::emit::EMITTED.with(|emitted| {
            for (stream, bytes) in emitted.borrow().iter() {
                streams
                    .entry(*stream)
                    .or_insert_with(Vec::new)
                    .extend_from_slice(bytes);
            }
        });
        streams
    }

    /// NCBI BLAST+ 2.17.0's frozen outputs of `-outfmt` with fields, for every program
    /// (`LOSAT/tests/outfmt0_manifest.tsv`, rows `tab.*`: id, program, query, subject, -outfmt).
    const CUSTOM_ROWS: [(&str, &str, &str, &str, &str); 12] = [
        ("tab.bn.delim_pipe.6", "blastn", "tests/db_search/nslice_q.fna", "tests/db_search/nslice.fna", "6 delim=| std sstrand"),
        ("tab.bn.align.7", "blastn", "tests/db_search/nslice_q.fna", "tests/db_search/nslice.fna", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        ("tab.bp.cov.6", "blastp", "tests/db_search/ids/ids_q.faa", "tests/fasta/outfmt11/arch_bp_s.faa", "6 std qcovs qcovhsp qcovus"),
        ("tab.bp.ids.7", "blastp", "tests/db_search/ids/ids_q.faa", "tests/db_search/ids/prot.faa", "7 std qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen"),
        ("tab.bx.qframe.6", "blastx", "tests/db_search/ids/ids_q.fna", "tests/fasta/outfmt11/arch_bp_s.faa", "6 qseqid qstart frames sframe qframe"),
        ("tab.bx.align.7", "blastx", "tests/db_search/ids/ids_q.fna", "tests/fasta/outfmt11/arch_bp_s.faa", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        ("tab.tbn.delim_pipe.6", "tblastn", "tests/db_search/tab_tbn_q.faa", "tests/db_search/tslice.fna", "6 delim=| std sframe sstrand"),
        ("tab.tbn.align.7", "tblastn", "tests/db_search/tab_tbn_q.faa", "tests/db_search/tslice.fna", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        ("tab.tbx.removed.6", "tblastx", "tests/db_search/tab_tbx_q.fna", "tests/db_search/xslice.fna", "6 -qseqid"),
        ("tab.tbx.align.7", "tblastx", "tests/db_search/tab_tbx_q.fna", "tests/db_search/xslice.fna", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        // Fields whose writer warns (taxdb, once per query block).
        ("tab.bn.tax.7", "blastn", "tests/db_search/nslice_q.fna", "tests/db_search/nslice.fna", "7 staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus"),
        ("tab.bx.delim_comma.6", "blastx", "tests/db_search/ids/ids_q.fna", "tests/fasta/outfmt11/arch_bp_s.faa", "6 delim=, qseqid sseqid scomnames"),
    ];

    // The output of the argv's -outfmt (6 or 7 with fields here) is stream 4, NCBI's bytes of the same
    // command line, in the five programs, and stream 3 NCBI's standard error of it (with the
    // warnings of its fields' writer); the formats 0, 6 and 7 and the HSP records are those of
    // the run without -outfmt.
    #[test]
    fn a_custom_outfmt_is_written_as_the_cli_writes_it() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        for (id, program, query, subject, outfmt) in CUSTOM_ROWS {
            let words = [program, "-query", query, "-subject", subject];
            let plain = registered_streams(&words, &read(query), &read(subject));
            let with = [&words[..], &["-outfmt", outfmt]].concat();
            let custom = registered_streams(&with, &read(query), &read(subject));
            let ncbi = read(&format!("tests/fixtures/outfmt0/{id}.out"));
            assert_eq!(text(&custom[&STREAM_CUSTOM]), text(&ncbi), "{id}");
            let stderr = std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                .unwrap_or_default();
            let diagnostics = custom.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(text(&diagnostics), text(&stderr), "{id}: stream 3");
            for stream in [0, 6, 7, STREAM_HITS] {
                assert_eq!(
                    plain.get(&stream),
                    custom.get(&stream),
                    "{id}: stream {stream}"
                );
            }
            assert!(!plain.contains_key(&STREAM_CUSTOM), "{id}");
        }
    }

    /// NCBI BLAST+ 2.17.0's frozen outputs of the formats besides 0, 6 and 7 (rows of
    /// `LOSAT/tests/outfmt0_manifest.tsv`: id and argv, run from LOSAT/): every program in 5, 8, 9,
    /// 12, 15 and 16, 10 and 20 with and without fields, the BLAST archive (11) of four programs
    /// with `-query_loc`, `-subject_loc` and deflines that the title does not give, and SAM (17).
    const FORMAT_ROWS: [(&str, &[&str]); 55] = [
        (
            "sm1.blastp.x16.soft.5",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/x16_seg_q.faa",
                "-subject",
                "tests/fasta/outfmt0/x16_seg_s.faa",
                "-seg",
                "yes",
                "-soft_masking",
                "true",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sm1.blastp.x16.soft.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/x16_seg_q.faa",
                "-subject",
                "tests/fasta/outfmt0/x16_seg_s.faa",
                "-seg",
                "yes",
                "-soft_masking",
                "true",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sm1.blastp.x16.soft.16",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/x16_seg_q.faa",
                "-subject",
                "tests/fasta/outfmt0/x16_seg_s.faa",
                "-seg",
                "yes",
                "-soft_masking",
                "true",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sm2.tblastx.mask.soft.5",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-soft_masking",
                "true",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sm2.tblastx.mask.soft.16",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-soft_masking",
                "true",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sm3.megablast.mask.hard.5",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-task",
                "megablast",
                "-soft_masking",
                "false",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sm3.megablast.mask.hard.16",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-task",
                "megablast",
                "-soft_masking",
                "false",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sm3.megablast.mask.hard.11",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-task",
                "megablast",
                "-soft_masking",
                "false",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.multi.megablast.8",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/multi_subject.fasta",
                "-task",
                "megablast",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.multi.megablast.12",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/multi_subject.fasta",
                "-task",
                "megablast",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.8",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.12",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.e2e.tblastn.xdrop_final.8",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-xdrop_gap_final",
                "1e9",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.e2d.tblastn.sloc.12",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-subject",
                "tests/fasta/outfmt0/e2d_n_sub8k.fa",
                "-subject_loc",
                "1001-6000",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.tblastx.nohit.8",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/tblastx_multi_subject.fasta",
                "-evalue",
                "1e-100",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.tblastx.nohit.12",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/tblastx_multi_subject.fasta",
                "-evalue",
                "1e-100",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.blastx.multi.8",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.blastx.multi.12",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.multi.megablast.9",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/multi_subject.fasta",
                "-task",
                "megablast",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.9",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.e2d.tblastn.sloc.9",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-subject",
                "tests/fasta/outfmt0/e2d_n_sub8k.fa",
                "-subject_loc",
                "1001-6000",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.tblastx.nohit.9",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/tblastx_multi_subject.fasta",
                "-evalue",
                "1e-100",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.blastx.multi.9",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-outfmt",
                "9",
            ],
        ),
        (
            "x5.tblastn.sloc2_470",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject_loc",
                "2-470",
                "-outfmt",
                "5",
            ],
        ),
        (
            "x5.blastx.ungapped_cbs0",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/x5_bx_q_mix.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "1",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.15",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "15",
            ],
        ),
        (
            "sa.x16.blastn.minus.15",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/x16n_q.fna",
                "-subject",
                "tests/fasta/outfmt0/x16n_s.fna",
                "-outfmt",
                "15",
            ],
        ),
        (
            "x16.blastx.ungapped_cbs0.16",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/x5_bx_q_mix.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "1",
                "-outfmt",
                "16",
            ],
        ),
        (
            "x16.blastx.ungapped_cbs0.15",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/x5_bx_q_mix.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "1",
                "-outfmt",
                "15",
            ],
        ),
        (
            "sa.x16.tblastn.sloc.16",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject_loc",
                "2-470",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sa.x16.tblastn.sloc.15",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject_loc",
                "2-470",
                "-outfmt",
                "15",
            ],
        ),
        (
            "sa.x16.tblastx.sub.15",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/x16x_q.fna",
                "-subject",
                "tests/fasta/outfmt0/x16x_s.fna",
                "-evalue",
                "0.01",
                "-max_target_seqs",
                "2",
                "-outfmt",
                "15",
            ],
        ),
        (
            "culling.tblastn.dom_cl2.11",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/culling_dom_q.faa",
                "-subject",
                "tests/fasta/outfmt0/culling_dom_s.fna",
                "-culling_limit",
                "2",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.sam.multi.17p",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/sam_multi_q.fna",
                "-subject",
                "tests/fasta/outfmt0/sam_multi_s.fna",
                "-task",
                "megablast",
                "-outfmt",
                "17",
            ],
        ),
        (
            "sa.sam.multi.17q",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/sam_multi_q.fna",
                "-subject",
                "tests/fasta/outfmt0/sam_multi_s.fna",
                "-task",
                "megablast",
                "-outfmt",
                "17 SQ",
            ],
        ),
        (
            "sa.arch.bn.g4.11",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt11/arch_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_s.fna",
                "-max_target_seqs",
                "1",
                "-dust",
                "20 64 1",
                "-query_loc",
                "10-300",
                "-subject_loc",
                "10-300",
                "-dbsize",
                "1000000",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.bn.qloc.11",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt11/arch_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_s.fna",
                "-query_loc",
                "10-300",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.bp.sloc.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt11/arch_bp_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-subject_loc",
                "10-100",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.bp.titles.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt11/arch_bp_titles.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbn.qloc.11",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt11/arch_tbn_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_tbn_s.fna",
                "-query_loc",
                "10-100",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbx.def.11",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt11/arch_tbx_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_tbx_s.fna",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbx.qloc.11",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt11/arch_tbx_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_tbx_s.fna",
                "-query_loc",
                "10-100",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbx.stitles.11",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt11/arch_tbx_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_tbx_s_titles.fna",
                "-outfmt",
                "11",
            ],
        ),
        (
            "tab.bn.csv.10",
            &[
                "blastn",
                "-query",
                "tests/db_search/nslice_q.fna",
                "-subject",
                "tests/db_search/nslice.fna",
                "-outfmt",
                "10 std sstrand",
            ],
        ),
        (
            "tab.bn.csvhdr.20",
            &[
                "blastn",
                "-query",
                "tests/db_search/nslice_q.fna",
                "-subject",
                "tests/db_search/nslice.fna",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.bp.csv.10",
            &[
                "blastp",
                "-query",
                "tests/db_search/ids/ids_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "10 std",
            ],
        ),
        (
            "tab.bp.csvhdr.20",
            &[
                "blastp",
                "-query",
                "tests/db_search/ids/ids_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.bx.csv.10",
            &[
                "blastx",
                "-query",
                "tests/db_search/ids/ids_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "10 std",
            ],
        ),
        (
            "tab.bx.csvhdr.20",
            &[
                "blastx",
                "-query",
                "tests/db_search/ids/ids_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.tbn.csv.10",
            &[
                "tblastn",
                "-query",
                "tests/db_search/tab_tbn_q.faa",
                "-subject",
                "tests/db_search/tslice.fna",
                "-outfmt",
                "10 std sstrand",
            ],
        ),
        (
            "tab.tbn.csvhdr.20",
            &[
                "tblastn",
                "-query",
                "tests/db_search/tab_tbn_q.faa",
                "-subject",
                "tests/db_search/tslice.fna",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.tbx.csv.10",
            &[
                "tblastx",
                "-query",
                "tests/db_search/tab_tbx_q.fna",
                "-subject",
                "tests/db_search/xslice.fna",
                "-outfmt",
                "10 std sstrand",
            ],
        ),
        (
            "tab.tbx.csvhdr.20",
            &[
                "tblastx",
                "-query",
                "tests/db_search/tab_tbx_q.fna",
                "-subject",
                "tests/db_search/xslice.fna",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "sxd.blastx.qloc.10",
            &[
                "blastx",
                "-query",
                "tests/db_search/bslice_q.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-query_loc",
                "1001-3000",
                "-outfmt",
                "10",
            ],
        ),
        (
            "sxd.blastx.qframe.20",
            &[
                "blastx",
                "-query",
                "tests/db_search/bslice_q.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-query_loc",
                "1002-3000",
                "-outfmt",
                "20",
            ],
        ),
    ];

    /// A SAM document without the CL tag of its `@PG` line (`db_search_fixtures.strip_sam_cl`, the
    /// rule of the SAM rows: the CLI writes its own command line there).
    fn without_cl(document: &[u8]) -> Vec<u8> {
        let text = String::from_utf8(document.to_vec()).unwrap();
        let head = "\n@PG\tID:0\tVN:2.17.0+";
        let Some(at) = text.find(&format!("{head}\tCL:")) else {
            return document.to_vec();
        };
        let start = at + head.len();
        let end = start + text[start..].find("\tPN:blastn\n").unwrap();
        format!("{}{}", &text[..start], &text[end..]).into_bytes()
    }

    // The output of an -outfmt besides 0, 6 and 7 is stream 4, NCBI's bytes of the same command
    // line (W11), and stream 3 NCBI's standard error of it; the formats 0, 6 and 7 and the HSP
    // records are those of the run without -outfmt. The archive's Bioseqs are the whole records
    // and the deflines as read; SAM's @PG line shows the run's command line.
    #[test]
    fn every_format_is_written_as_the_cli_writes_it() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        for (id, words) in FORMAT_ROWS {
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            let with = registered_streams(words, &query, &subject);
            let plain = registered_streams(&words[..words.len() - 2], &query, &subject);
            let mut ncbi = read(&format!("tests/fixtures/outfmt0/{id}.out"));
            let mut ours = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
            if at("-outfmt").starts_with("17") {
                let cl = format!("\tCL:losat {} \tPN:blastn\n", words.join(" "));
                assert!(text(&ours).contains(&cl), "{id}: {cl:?}");
                (ncbi, ours) = (without_cl(&ncbi), without_cl(&ours));
            }
            assert!(ours == ncbi, "{id}: stream 4");
            let stderr = std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                .unwrap_or_default();
            let diagnostics = with.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(text(&diagnostics), text(&stderr), "{id}: stream 3");
            for stream in [0, 6, 7, STREAM_HITS] {
                assert!(
                    plain.get(&stream) == with.get(&stream),
                    "{id}: stream {stream}"
                );
            }
        }
    }

    /// The rows of NCBI's frozen outputs (`LOSAT/tests/outfmt0_manifest.tsv`) for which `keep(program,
    /// format number, files column)` holds, as (fixture id, argv): the arguments of
    /// `run_oracle.search_argv` (none of the rows taken here quotes a value).
    fn manifest_rows(keep: impl Fn(&str, &str, &str) -> bool) -> Vec<(String, Vec<String>)> {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let manifest =
            std::fs::read_to_string(format!("{root}/tests/outfmt0_manifest.tsv")).unwrap();
        let mut lines = manifest.lines().filter(|line| !line.starts_with('#'));
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
        let mut rows = Vec::new();
        for line in lines {
            let row: Vec<&str> = line.split('\t').collect();
            let field = |name: &str| row.get(column(name)).copied().unwrap_or("");
            let outfmt = match field("outfmt") {
                "" => "0",
                outfmt => outfmt,
            };
            let number = outfmt
                .trim()
                .split(' ')
                .next()
                .unwrap()
                .trim_start_matches('+');
            let number = number.trim_start_matches('0');
            if !keep(field("program"), number, field("files")) {
                continue;
            }
            let extra = field("extra_args");
            assert!(
                !extra.contains(['"', '\'', '\\']),
                "{}",
                field("fixture_id")
            );
            let mut argv = vec![
                field("program").to_string(),
                "-query".to_string(),
                field("query").to_string(),
                "-subject".to_string(),
                field("subject").to_string(),
            ];
            if !field("task").is_empty() {
                argv.extend(["-task".to_string(), field("task").to_string()]);
            }
            argv.extend(extra.split_whitespace().map(str::to_string));
            argv.extend(["-outfmt".to_string(), outfmt.to_string()]);
            rows.push((field("fixture_id").to_string(), argv));
        }
        rows
    }

    // Outfmt 13 and 14 are a ZIP archive on stream 4: NCBI's side files, named as the CLI names them
    // for `-out results`, and NCBI's master file (`results.json` / `results.xml`), which lists them by
    // those names; stream 3 is NCBI's standard error, and the other streams are those of the run
    // without -outfmt (the 18 rows of the manifest with side files).
    #[test]
    fn outfmt_13_and_14_are_an_archive_of_ncbis_files() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        let rows = manifest_rows(|_, _, files| files == "xinclude");
        assert_eq!(rows.len(), 18);
        let programs: std::collections::BTreeSet<&str> =
            rows.iter().map(|(_, argv)| argv[0].as_str()).collect();
        assert_eq!(programs.len(), 5);
        for (id, argv) in &rows {
            let words: Vec<&str> = argv.iter().map(String::as_str).collect();
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let number: u32 = at("-outfmt").trim().parse().unwrap();
            let extension = if number == 13 { "json" } else { "xml" };
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            let with = registered_streams(&words, &query, &subject);
            let plain = registered_streams(&words[..words.len() - 2], &query, &subject);
            let mut expected = Vec::new();
            for n in 1.. {
                let path = format!("{root}/tests/fixtures/outfmt0/{id}.out_{n}.{extension}");
                let Ok(side) = std::fs::read(path) else {
                    break;
                };
                expected.push((format!("results_{n}.{extension}").into_bytes(), side));
            }
            let master = text(&read(&format!("tests/fixtures/outfmt0/{id}.out")))
                .replace(&format!("{id}.out_"), "results_");
            expected.push((master_name(number).into_bytes(), master.into_bytes()));
            let archive = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
            let files = crate::zip::read(&archive);
            assert_eq!(files.len(), expected.len(), "{id}: entries");
            for ((name, bytes), (want_name, want)) in files.iter().zip(&expected) {
                assert_eq!(text(name), text(want_name), "{id}");
                assert!(bytes == want, "{id}: {}", text(name));
            }
            let stderr = std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                .unwrap_or_default();
            let diagnostics = with.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(text(&diagnostics), text(&stderr), "{id}: stream 3");
            for stream in [0, 6, 7, STREAM_HITS] {
                assert!(
                    plain.get(&stream) == with.get(&stream),
                    "{id}: stream {stream}"
                );
            }
        }
    }

    // BLASTX writes the BLAST archive (outfmt 11) as NCBI does: the 57 rows of the manifest that a
    // run takes (not the one with -parse_deflines, nor the one with a `>?` gap line, which the
    // store refuses: maintainer decision 3), stream 4 and stream 3.
    #[test]
    fn blastx_writes_ncbis_archive() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        let gap_line = |path: &str| {
            read(path)
                .split(|b| *b == b'\n')
                .any(|l| l.starts_with(b">?"))
        };
        let rows = manifest_rows(|program, number, _| program == "blastx" && number == "11");
        let rows: Vec<_> = rows
            .into_iter()
            .filter(|(_, argv)| !argv.iter().any(|word| word == "-parse_deflines"))
            .filter(|(_, argv)| {
                let file = |flag: &str| &argv[argv.iter().position(|w| w == flag).unwrap() + 1];
                !gap_line(file("-query")) && !gap_line(file("-subject"))
            })
            .collect();
        assert_eq!(rows.len(), 57);
        for (id, argv) in &rows {
            let words: Vec<&str> = argv.iter().map(String::as_str).collect();
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            let with = registered_streams(&words, &query, &subject);
            let ncbi = read(&format!("tests/fixtures/outfmt0/{id}.out"));
            let ours = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
            assert!(ours == ncbi, "{id}: stream 4");
            let stderr = std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                .unwrap_or_default();
            let diagnostics = with.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(text(&diagnostics), text(&stderr), "{id}: stream 3");
        }
    }

    // The store keeps the deflines that the archive cannot make from the title (a control byte, white
    // space after the `>` or at the end; not a CRLF line end), aside from the records: only a run that
    // writes the archive gives the records a source, which they hold no longer after it.
    #[test]
    fn the_store_keeps_the_deflines_that_the_title_does_not_give() {
        let fasta: &[u8] = b">a one\r\nACGT\n>b\twith tab\nGG\n>  c\nTT\n>d \nAA\n>e\r\nCC\n";
        let (handle, _) = crate::store::register("blastn", 0, fasta).unwrap();
        let (subject, _) = crate::store::register("blastn", 1, b">s\nACGTACGT\n").unwrap();
        let sources = |whole: bool| {
            crate::store::with_inputs(
                Program::Blastn,
                handle,
                subject,
                whole,
                |queries, subjects| {
                    let source = |record: &FastaRecord| {
                        record.source.as_ref().map(|source| source.defline.clone())
                    };
                    (
                        queries.iter().map(source).collect::<Vec<_>>(),
                        subjects.iter().map(source).collect::<Vec<_>>(),
                    )
                },
            )
            .unwrap()
        };
        assert_eq!(sources(false), (vec![None; 5], vec![None]));
        let kept = |line: &[u8]| Some(line.to_vec());
        assert_eq!(
            sources(true),
            (
                vec![
                    kept(b""),
                    kept(b">b\twith tab"),
                    kept(b">  c"),
                    kept(b">d "),
                    kept(b"")
                ],
                vec![kept(b"")]
            )
        );
        assert_eq!(sources(false), (vec![None; 5], vec![None]));
        crate::store::release(handle).unwrap();
        crate::store::release(subject).unwrap();
        // A first record without a defline, a kept defline at the end of the input without a line
        // end, and CR line ends.
        for (fasta, kept_lines) in [
            (
                &b"ACGT\n>b\tx\nGG\n"[..],
                vec![b"".to_vec(), b">b\tx".to_vec()],
            ),
            (&b">a\nAC\n>b \t"[..], vec![b"".to_vec(), b">b \t".to_vec()]),
            (
                &b">a \rAC\r>b\rGG\r"[..],
                vec![b">a ".to_vec(), b"".to_vec()],
            ),
        ] {
            let (handle, _) = crate::store::register("blastn", 0, fasta).unwrap();
            let (subject, _) = crate::store::register("blastn", 1, b">s\nACGTACGT\n").unwrap();
            let defline = |whole: bool| {
                crate::store::with_inputs(Program::Blastn, handle, subject, whole, |queries, _| {
                    queries
                        .iter()
                        .map(|record| record.source.as_ref().map(|source| source.defline.clone()))
                        .collect::<Vec<_>>()
                })
                .unwrap()
            };
            assert_eq!(defline(false), vec![None; kept_lines.len()], "{fasta:?}");
            assert_eq!(
                defline(true),
                kept_lines.into_iter().map(Some).collect::<Vec<_>>(),
                "{fasta:?}"
            );
            crate::store::release(handle).unwrap();
            crate::store::release(subject).unwrap();
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
        // The output of -outfmt "6 std" comes from the tabular search: NCBI's outfmt 6 with
        // -num_alignments 3, without A_split, which the hit list of outfmt 0 keeps. The other
        // streams do not change.
        let custom = [&words[..], &["-outfmt", "6 std"]].concat();
        let with = self::streams(&custom, &query, &subject).unwrap();
        assert_eq!(with[&STREAM_CUSTOM], SPLIT_OUT6);
        for stream in [0, 6, 7, STREAM_HITS, STREAM_DIAGNOSTICS] {
            assert_eq!(with.get(&stream), streams.get(&stream), "stream {stream}");
        }
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

    // The reverse: a -num_descriptions of 2,147,483,598 or more (N + 50 overflows) makes BLASTN's preliminary hit list
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
        let rejected =
            |extra: &[&str], group: &[u32]| group_rejected(&mut command(extra), &formats(group));
        assert!(rejected(&["-num_alignments", "0"], &[6, 7]));
        assert!(!rejected(&["-num_alignments", "0"], &[0]));
        assert!(rejected(&["-num_descriptions", "2147483598"], &[0]));
        assert!(!rejected(&["-num_descriptions", "2147483597"], &[0]));
        assert!(!rejected(&["-num_descriptions", "1073741824"], &[0]));
        assert!(!rejected(&[], &[6, 7, 0]));
        let mut blastp = parse(&["blastp", "-num_alignments", "0"]).unwrap().1;
        assert!(group_rejected(&mut blastp, &formats(&[6])));
        assert!(!group_rejected(&mut blastp, &formats(&[0])));
    }

    // A coding sequence of the protein of `q_idx_is_the_registered_index`.
    const CODING: &[u8] = b"ATGAAAACCGCGTATATTGCGAAACAGCGTCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGTCAGCTGGAAGAACGTCTGGGCCTGATTGAAGTGCAGGCGCCGATTCTGAGCCGTGTGGGCGATGGCACCCAGGATAACCTGAGCGGCGCGGAAAAAGCGGTGCAGGTGAAAGTGAAAGCGCTGCCGGATGCGCAG";

    // BLASTX searches a nucleotide query against protein subjects and writes every format, the
    // HSP records and the warnings as the other programs do, one search for each hit list size:
    // with -num_alignments 0 the tabular search fails alone with NCBI's error on stream 3, and
    // format 0 is written.
    #[test]
    fn blastx_runs_as_the_other_programs_do() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [FastaRecord::new("q", b"q", CODING)];
        let subject = [FastaRecord::new("s", b"s", protein)];
        let words = ["blastx", "-query", "q.fna", "-subject", "s.faa"];
        let plain = streams(&words, &query, &subject).unwrap();
        let row = text(&plain[&6]);
        assert!(
            row.starts_with("q\ts\t100.000\t66\t0\t0\t1\t198\t1\t66\t"),
            "{row}"
        );
        assert_eq!(row.lines().count(), 1);
        assert!(text(&plain[&7]).contains("# BLASTX "));
        assert!(text(&plain[&0]).starts_with("BLASTX "));
        let records = text(&plain[&1]);
        assert_eq!(records.lines().count(), 1);
        assert!(records.contains("\"q_idx\":0"), "{records}");
        assert!(!records.contains("\"out0\":null"), "{records}");
        assert!(plain.get(&3).is_none_or(|bytes| bytes.is_empty()));
        let zero = streams(
            &[&words[..], &["-num_alignments", "0"]].concat(),
            &query,
            &subject,
        )
        .unwrap();
        assert!(text(&zero[&0]).contains("Sequences producing significant alignments"));
        for stream in [1, 6, 7] {
            assert!(
                zero.get(&stream).is_none_or(|bytes| bytes.is_empty()),
                "stream {stream}"
            );
        }
        // The CLI's standard error for `-outfmt 6`: the hit list size's warning, then NCBI's error.
        assert_eq!(
            text(&zero[&3]),
            "Warning: [blastx] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"
        );
        // Both searches fail, each for its own reason (the tabular formats' hit list is checked
        // before the e-value): the run fails with the tabular search's error, as validate.
        let both = [&words[..], &["-num_alignments", "0", "-evalue", "0"]].concat();
        let error = run(&both, &query, &subject).unwrap_err();
        assert!(error.contains("No hits are being saved"), "{error}");
        assert_eq!(validate(&both).unwrap_err(), error);
        let mut format0 = parse(&both).unwrap().1;
        assert!(
            group_rejected(&mut format0, &formats(&[0]))
                && group_rejected(&mut format0, &formats(&[6, 7]))
        );
    }

    // BLASTX reads -query_loc and -subject_loc as the CLI does (DW-39): the rows give record
    // coordinates, a query that the range skips keeps its registered index in the HSP records,
    // and validate reads the ranges where the run reads them, after LOSAT's limits of the
    // options and the subject range first.
    #[test]
    fn blastx_ranges_as_the_other_programs() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [
            FastaRecord::new("short", b"short", b"ACGTACGT"),
            FastaRecord::new("q", b"q", CODING),
        ];
        let subject = [FastaRecord::new("s", b"s", protein)];
        let words = [
            "blastx",
            "-query",
            "q.fna",
            "-subject",
            "s.faa",
            "-query_loc",
            "31-198",
            "-subject_loc",
            "5-66",
        ];
        assert!(validate(&words).is_ok());
        let ranged = streams(&words, &query, &subject).unwrap();
        let row = text(&ranged[&6]);
        assert!(
            row.starts_with("q\ts\t100.000\t56\t0\t0\t31\t198\t11\t66\t"),
            "{row}"
        );
        let records = text(&ranged[&1]);
        assert_eq!(records.lines().count(), 1);
        assert!(records.contains("\"q_idx\":1"), "{records}");
        assert_eq!(
            searched_queries(Some("31-198"), "BLASTX", &query),
            Some(vec![1])
        );
        // The errors of the argv, as the run reports them.
        for extra in [
            &["-subject_loc", "5-2", "-query_loc", "9-1"][..],
            &["-query_loc", "9-1", "-evalue", "0"],
            &["-subject_loc", "5-2", "-matrix", "IDENTITY"],
        ] {
            let argv = [&words[..5], extra].concat();
            let error = run(&argv, &query, &subject).unwrap_err();
            assert_eq!(validate(&argv).unwrap_err(), error, "{extra:?}");
        }
        let error =
            validate(&[&words[..5], &["-subject_loc", "5-2", "-query_loc", "9-1"]].concat())
                .unwrap_err();
        assert!(error.contains("subject location"), "{error}");
        // The subject range comes before NCBI's errors of the options (IDENTITY, which BLASTX
        // does not take).
        let error =
            validate(&[&words[..5], &["-subject_loc", "5-2", "-matrix", "IDENTITY"]].concat())
                .unwrap_err();
        assert!(error.contains("subject location"), "{error}");
        // A subject range that starts past the subject's end is an error of the records, which
        // validate does not see.
        let past = [&words[..5], &["-subject_loc", "70-80"]].concat();
        assert!(validate(&past).is_ok());
        let error = run(&past, &query, &subject).unwrap_err();
        assert!(error.contains("Invalid from coordinate"), "{error}");
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
