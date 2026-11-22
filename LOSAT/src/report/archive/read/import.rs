//! `-import_search_strategy <file>`: what NCBI's `RecoverSearchStrategy` and
//! `s_ImportSearchStrategy` (blast_app_util.cpp:532-628) make of a strategy file before the
//! running program takes it as its own options and inputs (`SetOptionsForSavedStrategy`,
//! blast_args.cpp:3800-3836; each program's `apply_strategy`).
//!
//! The file is read as `ExtractBlast4Request` reads it (remote_blast.cpp:2257-2331): first as a
//! `Blast4-get-search-strategy-reply`, then, from the start of the input, as a `Blast4-request`;
//! a serial error of either read is NCBI's `Failed to read search strategy`. The options are the
//! strategy's (`read::options::import_request`); the queries are the FASTA text NCBI writes and
//! reads back (`read::queries`), handed to the program as an unnamed temporary file (NCBI's
//! `CTmpFile`); the subject is the strategy's database or its Bioseqs.
//!
//! LOSAT refuses, as not yet ported, what it has no reader for: binary ASN.1 and XML strategies,
//! `seq-loc-list` queries and subjects, PSSM queries in the programs that take them, and GI lists.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:609-628

use std::io::{Read as IoRead, Seek};
use std::path::Path;

use anyhow::Result;

use super::options::{import_request, EProgram, ImportedOptions, SeqRange, Setter};
use super::queries::{render_queries, RenderError};
use super::request::{Queries, Request, Subject};
use crate::cli::{Commands, NativeError};
use crate::ncbi_exception::NcbiException;
use crate::report::archive::options::{Param, ParamValue};
use crate::report::serial::asn_types::read::ReadError;
use crate::report::serial::asn_types::stream::{guess_encoding, Encoding, GuessStop};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::Value;

/// What the running program's command line says before the strategy is read.
pub(crate) struct Command<'a> {
    /// `-import_search_strategy`'s value (`-`: the standard input).
    pub(crate) path: &'a Path,
    /// The strategy's bytes when the caller holds them (the library API and LOSAT Web, which have
    /// no file to name): read in place of the file `path`.
    pub(crate) bytes: Option<&'a [u8]>,
    /// The running program's name in LOSAT's refusals (`BLASTN`, ...).
    pub(crate) program: &'static str,
    /// `override_query`: `-query` has a value other than `-`.
    pub(crate) override_query: bool,
    /// `override_subject` (`CBlastDatabaseArgs::HasBeenSet`): `-db` or `-subject` given.
    pub(crate) override_subject: bool,
    /// `-use_index` given with a value other than its default (false).
    pub(crate) use_index_set: bool,
}

/// The strategy's subject when the command line does not override it.
pub(crate) enum ImportedSubject {
    /// `s_ImportDatabase`: the database name as stored.
    Database(String),
    /// `s_ImportSubjects`: the Bioseqs, each searched whole.
    Sequences(Vec<Value>),
}

/// What the import hands to the running program.
pub(crate) struct Imported {
    /// `CImportStrategy`'s options (the builder's setters, the database filters, the query
    /// range and PSI-BLAST's iteration count).
    pub(crate) options: ImportedOptions,
    /// The queries as NCBI's temporary FASTA file, `None` when `-query` overrides them.
    pub(crate) queries: Option<Vec<u8>>,
    /// The subject, `None` when `-db` or `-subject` overrides it.
    pub(crate) subject: Option<ImportedSubject>,
    /// The warnings posted so far, in order (the text after `Warning: [<program>] `).
    pub(crate) warnings: Vec<String>,
}

impl Imported {
    /// The strategy's program (`m_OptsHandle->GetOptions().GetProgramType()` and the task).
    pub(crate) fn program(&self) -> EProgram {
        self.options.options.program
    }

    /// `GetQueryRange()`, `None` when it is empty or `-query` overrides the queries (the range
    /// is set only with the strategy's queries, blast_app_util.cpp:594-598).
    pub(crate) fn query_range(&self) -> Option<SeqRange> {
        let range = self.options.options.query_range;
        (self.queries.is_some() && range != SeqRange::EMPTY).then_some(range)
    }
}

/// `s_ImportSearchStrategy`'s `CInputException` for a file that does not read as a request.
fn failed_to_read() -> anyhow::Error {
    crate::blastinput::app::options_error("Failed to read search strategy")
}

/// LOSAT's refusal of a part of a strategy that it does not read yet.
fn unported(what: &str, program: &str) -> anyhow::Error {
    anyhow::anyhow!("{what} is not supported by LOSAT's {program}")
}

/// What `CATCH_ALL` (blast_app_util.hpp:168-262) reports of an exception the import throws.
fn caught(exception: &NcbiException) -> anyhow::Error {
    NativeError::exception(exception).into()
}

/// The error of a read of `ExtractBlast4Request`'s second attempt: a `CSerialException` is
/// `Failed to read search strategy`; any other exception leaves the import for `CATCH_ALL`.
fn read_error(error: ReadError, program: &str) -> anyhow::Error {
    match error {
        ReadError::Serial(_) => failed_to_read(),
        ReadError::Exception(exception) => caught(&exception),
        // `CEofException` is a `CIOException` other than `eFlush`: `CATCH_ALL` ends the run
        // without a post and with its status, 0.
        ReadError::Eof => NativeError::silent_end().into(),
        ReadError::Std(what) => NativeError::std_exception(what).into(),
        ReadError::Unsupported(what) => unported(&what, program),
    }
}

/// The guess's answer for one read of `ExtractBlast4Request`: `Ok(true)` for text, `Ok(false)`
/// for an input that NCBI's switch does not read (`default`), LOSAT's refusal for binary ASN.1
/// and XML.
fn text_encoding(encoding: &Result<Encoding, GuessStop>, program: &str) -> Result<bool> {
    match encoding {
        Ok(Encoding::Text) => Ok(true),
        Ok(Encoding::Other) => Ok(false),
        Ok(Encoding::Binary) => Err(unported("a search strategy in binary ASN.1", program)),
        Ok(Encoding::Xml) => Err(unported("a search strategy in XML", program)),
        Err(GuessStop::Exception(exception)) => Err(caught(exception)),
        Err(GuessStop::Undecided(detail)) => Err(unported(
            &format!("a search strategy file whose format NCBI's format guess decides by reading memory past it ({detail})"),
            program,
        )),
    }
}

/// The posts of NCBI's reader: LOSAT refuses a strategy on which the reader posts, as it does
/// not write those posts in the import's order yet.
fn no_posts(posts: &[crate::report::serial::asn_types::Post], program: &str) -> Result<()> {
    if posts.is_empty() {
        Ok(())
    } else {
        Err(unported(
            "a search strategy whose reading posts messages",
            program,
        ))
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2258-2273
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2295-2309
/// `ExtractBlast4Request(in)`: the request of a strategy file, the reply type first.
fn extract_request<R: IoRead + Seek>(input: R, program: &str) -> Result<Value> {
    let guessed = guess_encoding(input);
    // The first attempt: its exceptions are caught, and the guess of an input that is not read
    // here is a failed attempt too.
    let first = match guessed.encoding {
        Ok(Encoding::Text) => {
            let mut stream = guessed.into_text();
            let (read, posts) = stream.read(named::BLAST4_GET_SEARCH_STRATEGY_REPLY);
            match read {
                Ok(value) => {
                    no_posts(&posts, program)?;
                    return Ok(value);
                }
                Err(_) => stream.rewind(),
            }
        }
        // Both types are a SEQUENCE (`Blast4-get-search-strategy-reply ::= Blast4-request`): a
        // binary stream whose first byte is not its tag (0x30) fails both reads, and NCBI's
        // `s_ImportSearchStrategy` reports `Failed to read search strategy` (a text file with a
        // byte order mark is guessed binary).
        Ok(Encoding::Binary) => {
            let (taken, _) = guessed.into_parts();
            if taken.first() != Some(&0x30) {
                return Err(failed_to_read());
            }
            return Err(unported("a search strategy in binary ASN.1", program));
        }
        Ok(Encoding::Xml) => {
            text_encoding(&guessed.encoding, program)?;
            unreachable!("binary and XML strategies are refused")
        }
        Ok(Encoding::Other) | Err(GuessStop::Exception(_)) => {
            let (_, mut input) = guessed.into_parts();
            input.seek(std::io::SeekFrom::Start(0)).map(|_| input)
        }
        Err(GuessStop::Undecided(_)) => {
            text_encoding(&guessed.encoding, program)?;
            unreachable!("an undecided guess is refused")
        }
    };
    // `in.seekg(0)` of an input that cannot seek (a pipe) fails its stream: NCBI's second guess
    // reads nothing and its `default` throws.
    let input = first.map_err(|_| failed_to_read())?;
    let guessed = guess_encoding(input);
    if !text_encoding(&guessed.encoding, program)? {
        // `NCBI_THROW(CSerialException, eInvalidData, "Unrecognized input format ")`.
        return Err(failed_to_read());
    }
    let mut stream = guessed.into_text();
    let (read, posts) = stream.read(named::BLAST4_REQUEST);
    no_posts(&posts, program)?;
    read.map_err(|error| read_error(error, program))
}

/// `GetImportSearchStrategyStream`: the file (`AsInputFile()`; `-` is the standard input), or
/// the argument layer's error for a file that cannot be opened.
fn open_strategy(path: &Path, program: &str) -> Result<crate::blastinput::input_files::InputFile> {
    if path.as_os_str() == "-" {
        return crate::blastinput::input_files::open_input(path, "import_search_strategy", program);
    }
    std::fs::File::open(path)
        .map(crate::blastinput::input_files::InputFile::from)
        .map_err(|_| crate::cli::inaccessible("import_search_strategy", path))
}

/// `Blast_ProgramNameFromType` of the options' program type (blast_program.c), for the PSSM
/// error of `s_ImportPssm`.
fn program_type_name(program: EProgram) -> &'static str {
    match program {
        EProgram::Blastn
        | EProgram::Megablast
        | EProgram::DiscMegablast
        | EProgram::VecScreen
        | EProgram::Mapper => "blastn",
        EProgram::PhiBlastp | EProgram::Blastp | EProgram::DeltaBlast => "blastp",
        EProgram::PhiBlastn => "phiblastn",
        other => other.task_name(),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:543-549
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:565-599
/// `RecoverSearchStrategy` up to `SetOptionsForSavedStrategy`: the strategy's options, queries
/// and subject, and the warnings posted on the way. An error comes after the warnings posted
/// before it, which the caller writes first ([`Imported::warnings`] of the error's
/// [`ImportFailure`]).
pub(crate) fn recover(command: &Command<'_>) -> std::result::Result<Imported, ImportFailure> {
    let mut warnings = Vec::new();
    let fail = |warnings: &Vec<String>, error: anyhow::Error| ImportFailure {
        warnings: warnings.clone(),
        error,
    };
    if command.use_index_set {
        warnings.push("Overriding megablast BLAST DB indexed options in saved strategy".into());
    }
    let program = command.program;
    let value = match command.bytes {
        Some(bytes) => extract_request(std::io::Cursor::new(bytes), program),
        None => {
            open_strategy(command.path, program).and_then(|input| extract_request(input, program))
        }
    }
    .map_err(|e| fail(&warnings, e))?;
    let request =
        Request::from_strategy(value).map_err(|exception| fail(&warnings, caught(&exception)))?;
    let options = import_request(&request).map_err(|e| fail(&warnings, caught(&e)))?;
    let subject = if command.override_subject {
        warnings.push("Overriding database/subject in saved strategy".into());
        None
    } else {
        Some(import_subject(request.subject, &options, program).map_err(|e| fail(&warnings, e))?)
    };
    let queries = if command.override_query {
        warnings.push("Overriding query in saved strategy".into());
        None
    } else {
        Some(import_queries(&request.queries, &options, program).map_err(|e| fail(&warnings, e))?)
    };
    Ok(Imported {
        options,
        queries,
        subject,
        warnings,
    })
}

/// An import that failed, after the warnings it posted.
pub(crate) struct ImportFailure {
    pub(crate) warnings: Vec<String>,
    pub(crate) error: anyhow::Error,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:452-465
/// `s_ImportDatabase` and `s_ImportSubjects`: the strategy's subject.
fn import_subject(
    subject: Subject,
    options: &ImportedOptions,
    program: &str,
) -> Result<ImportedSubject> {
    match subject {
        Subject::Database(name) => {
            let built = &options.options;
            if let Some(limit) = &built.entrez_query {
                return Err(NativeError::std_exception(&format!(
                    "Entrez query '{limit}' will not be processed locally.\nPlease use the -remote option."
                ))
                .into());
            }
            if built.gi_list.is_some() || built.negative_gi_list.is_some() {
                return Err(unported("a search strategy with a GI list", program));
            }
            Ok(ImportedSubject::Database(
                String::from_utf8_lossy(&name).into_owned(),
            ))
        }
        Subject::Sequences(bioseqs) => Ok(ImportedSubject::Sequences(bioseqs)),
        Subject::SeqLocList => Err(unported(
            "a search strategy whose subject is a list of sequence locations",
            program,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:372-392
/// `s_ImportPssm` and `s_ImportQueries`: the queries as the text of NCBI's temporary FASTA file
/// (`ImportRun::queries_input` makes the file where a run reads one).
fn import_queries(queries: &Queries, options: &ImportedOptions, program: &str) -> Result<Vec<u8>> {
    match queries {
        Queries::Pssm(_) => {
            if matches!(program, "BLASTN" | "BLASTP" | "BLASTX" | "TBLASTX") {
                Err(crate::blastinput::app::engine_error(&format!(
                    "PSSM found in saved strategy, but not supported for {}",
                    program_type_name(options.options.program)
                )))
            } else {
                Err(unported("a search strategy whose query is a PSSM", program))
            }
        }
        Queries::SeqLocList(_) => Err(unported(
            "a search strategy whose queries are a list of sequence locations",
            program,
        )),
        Queries::BioseqSet(set) => {
            let text =
                render_queries(V::new(named::BIOSEQ_SET, set)).map_err(|error| match error {
                    RenderError::Exception(exception) => caught(&exception),
                    RenderError::Unsupported(part) => unported(&part, program),
                })?;
            Ok(text)
        }
    }
}

/// The arguments of the running program's command line as NCBI's `CArgs` holds them: every
/// argument given, with its text, and every argument with a default that is not given, with the
/// default's text.
pub(crate) struct NcbiArguments {
    /// `(name, text)`, in `CArgs::GetAll`'s order (by name, `operator<` of the names,
    /// ncbiargs.cpp:983-986).
    all: Vec<(String, String)>,
}

impl NcbiArguments {
    /// The arguments of `program` (NCBI's name, `blastn` ...) given as `given` (name and text, a
    /// flag's text empty), over the defaults of NCBI's argument descriptions.
    pub(crate) fn new(program: &str, given: &[(String, String)]) -> Self {
        let mut all: std::collections::BTreeMap<String, String> =
            crate::ncbi_xmlhelp::ncbi_argument_defaults(program)
                .into_iter()
                .map(|(name, text)| (name.to_string(), text))
                .collect();
        for (name, text) in given {
            all.insert(name.clone(), text.clone());
        }
        Self {
            all: all.into_iter().collect(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3698-3726
/// `can_override` of `x_IssueWarningsForIgnoredOptions`.
const CAN_OVERRIDE: &[&str] = &[
    "query",
    "query_loc",
    "subject",
    "subject_loc",
    "lcase_masking",
    "db",
    "dbsize",
    "entrez_query",
    "db_soft_mask",
    "db_hard_mask",
    "use_index",
    "index_name",
    "strand",
    "parse_deflines",
    "out",
    "outfmt",
    "num_descriptions",
    "num_alignments",
    "max_target_seqs",
    "remote",
    "num_threads",
    "import_search_strategy",
    "remote_verbose",
    "verbose",
];

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3728-3773
/// `has_defaults` of `x_IssueWarningsForIgnoredOptions` for the strategy's program and task
/// (cmdline_flags.cpp: `kDfltArgCompBasedStats` `2`, `...Delta` `1`, `kDfltArgNoFiltering` `no`,
/// `kDfltArgSegFiltering` `12 2.2 2.5`, `kDfltArgDustFiltering` `20 64 1`, the lookup-table
/// masking `false`/`true`, `kDfltArgOldStyleIndex` `true`, `kDfltArgMaskLevel` `-1`,
/// `kDfltOffDiagonalRange` and `kDfltArgMaxIntronLength` 0, `kDfltArgPSINumIterations` 1;
/// `BLAST_EXPECT_VALUE` 10, the gap triggers 22 and 27, `BLAST_GENETIC_CODE` 1,
/// `PSI_PSEUDO_COUNT_CONST` 0 and `PSI_INCLUSION_ETHRESH` 0.002 as `NStr`'s texts).
fn has_default(name: &str, program: EProgram, task: &str) -> Option<&'static str> {
    let protein = !program.query_is_nucleotide();
    Some(match name {
        "comp_based_stats" => {
            if matches!(program, EProgram::RpsBlast | EProgram::RpsTblastn) {
                "1"
            } else {
                "2"
            }
        }
        "evalue" => "10",
        "old_style_index" => "true",
        "seg" if protein => {
            if task.contains("blastp") || task.contains("psiblast") {
                "no"
            } else {
                "12 2.2 2.5"
            }
        }
        "dust" if !protein => "20 64 1",
        "soft_masking" => {
            if protein {
                "false"
            } else {
                "true"
            }
        }
        "gap_trigger" => {
            if protein {
                "22"
            } else {
                "27"
            }
        }
        "off_diagonal_range" => "0",
        "mask_level" => "-1",
        "max_intron_length" => "0",
        "query_gencode" | "db_gencode" => "1",
        "pseudocount" => "0",
        "inclusion_ethresh" => "0.002",
        "num_iterations" => "1",
        _ => return None,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3775-3797
/// `x_IssueWarningsForIgnoredOptions`: the warnings (text after `Warning: [<program>] `) for the
/// arguments of the command line that the strategy's options replace, given the strategy's
/// program and task.
pub(crate) fn ignored_option_warnings(
    arguments: &NcbiArguments,
    program: EProgram,
    task: &str,
) -> Vec<String> {
    let mut warnings = Vec::new();
    for (name, text) in &arguments.all {
        let warning = format!("{name} cannot be overridden when using a search strategy");
        let default = if name == "task" {
            Some(task)
        } else {
            has_default(name, program, task)
        };
        if let Some(default) = default {
            if default == text {
                continue;
            }
            if name == "task" && text == "megablast" {
                continue;
            }
            warnings.push(warning.clone());
        }
        if !CAN_OVERRIDE.contains(&name.as_str()) {
            warnings.push(warning);
        }
    }
    warnings
}

/// What a run that imports a strategy carries into the program's own flow (each program's
/// `import` field): the files that stand for its query and subject inputs, and the warnings that
/// NCBI posts during the import.
#[derive(Clone, Debug, Default)]
pub(crate) struct ImportRun {
    /// The strategy's queries (the text of NCBI's temporary FASTA file), read in place of
    /// `-query` (`queries_input`).
    pub(crate) queries: Option<std::sync::Arc<Vec<u8>>>,
    /// The strategy's subject Bioseqs as FASTA text (`queries::render_subjects`), read in place
    /// of `-subject` and named in no report (`User specified sequence set.`).
    pub(crate) subjects: Option<std::sync::Arc<Vec<u8>>>,
    /// The strategy's bytes, which a caller without a file gives before the import reads them
    /// (`apply_strategy` reads them in place of the file of `-import_search_strategy`); `None` on
    /// the command line.
    pub(crate) source: Option<std::sync::Arc<Vec<u8>>>,
    /// `RecoverSearchStrategy`'s warnings, posted before the program's own.
    pub(crate) warnings: Vec<String>,
    /// `x_IssueWarningsForIgnoredOptions`' warnings, posted after the formatting handler's.
    pub(crate) option_warnings: Vec<String>,
    /// The running program's NCBI name (`Start::app`), the prefix of the warnings.
    pub(crate) app: &'static str,
    /// The effective search space is the strategy's `EffectiveSearchSpace`, which does not set
    /// `OLD_FSC` as `-searchsp` does (blast_args.cpp:303-307).
    pub(crate) searchsp_from_strategy: bool,
    /// The strategy's query range (`Start::query_range`), which the program's query options
    /// take in place of `-query_loc`'s.
    pub(crate) query_range: Option<crate::blastinput::seq_range::SequenceRange>,
    /// The strategy's parameters as its setters record them, in the file's order: the
    /// archive's options (`options::OptionInputs::imported`).
    pub(crate) recorded: Vec<crate::report::archive::options::Param>,
}

impl ImportRun {
    /// `RecoverSearchStrategy`'s warnings as the running program posts them.
    pub(crate) fn warning_lines(&self) -> Vec<u8> {
        Self::lines(&self.warnings, self.app)
    }

    /// `x_IssueWarningsForIgnoredOptions`' warnings as the running program posts them.
    pub(crate) fn option_warning_lines(&self) -> Vec<u8> {
        Self::lines(&self.option_warnings, self.app)
    }

    /// The lines of `warnings` as the program `program` (NCBI's name) posts them.
    pub(crate) fn lines(warnings: &[String], program: &str) -> Vec<u8> {
        let mut text = Vec::new();
        for warning in warnings {
            text.extend_from_slice(format!("Warning: [{program}] {warning}\n").as_bytes());
        }
        text
    }

    /// The query input of the run: the strategy's queries, else `-query`.
    pub(crate) fn open_query(
        import: Option<&ImportRun>,
        path: &Path,
        program: &str,
    ) -> Result<crate::blastinput::input_files::InputFile> {
        match import.and_then(|run| run.queries.as_ref()) {
            Some(text) => Self::input_of("query", text),
            None => crate::blastinput::input_files::open_input(path, "query", program),
        }
    }

    /// An input file of the strategy's `text` (an unnamed temporary file, which a run that reads
    /// its sequences from a file opens in place of `-query` or `-subject`).
    pub(crate) fn input_of(
        role: &str,
        text: &[u8],
    ) -> Result<crate::blastinput::input_files::InputFile> {
        Ok(crate::blastinput::input_files::InputFile::from(
            crate::blastinput::input_files::temporary_input(role, text)?,
        ))
    }

    /// `old_fsc` of a run: `-searchsp` sets it, the strategy's effective search space does not.
    pub(crate) fn old_fsc(import: Option<&ImportRun>, searchsp: Option<i64>) -> bool {
        (searchsp.is_some() && !import.is_some_and(|run| run.searchsp_from_strategy))
            || crate::blastinput::app::old_fsc_environment()
    }
}

/// The start of an import: the strategy read, and the program's defaults to apply it to.
pub(crate) struct Start {
    pub(crate) imported: Imported,
    /// The arguments of the command line as NCBI's `CArgs` holds them.
    pub(crate) arguments: NcbiArguments,
    /// The running program's NCBI name, the prefix of its warnings (`blastn`, ...): the program
    /// of the command line, also when it runs another program's strategy ([`RunAs`]).
    pub(crate) app: &'static str,
}

/// `RecoverSearchStrategy` for the running program `ncbi_name` (`blastp` ...): the strategy, or
/// the error after the warnings posted before it (written to `diagnostics` here, as NCBI posts
/// them before its `CATCH_ALL` reports the error; the command line's are its standard error).
pub(crate) fn start(
    command: &Command<'_>,
    ncbi_name: &'static str,
    diagnostics: &mut dyn std::io::Write,
) -> Result<Start> {
    let arguments = NcbiArguments::new(ncbi_name, &crate::cli::given_arguments());
    match recover(command) {
        Ok(imported) => Ok(Start {
            imported,
            arguments,
            app: ncbi_name,
        }),
        Err(failure) => {
            let _ = diagnostics.write_all(&ImportRun::lines(&failure.warnings, ncbi_name));
            Err(failure.error)
        }
    }
}

impl Start {
    /// The run state of the import, its queries and subjects in place (the subjects rendered for
    /// a reader with or without `-parse_deflines`), and the warnings of the options for the
    /// strategy's program and task.
    pub(crate) fn run(&self, parse_deflines: bool, program: &str) -> Result<ImportRun> {
        let subjects = match &self.imported.subject {
            Some(ImportedSubject::Sequences(bioseqs)) => {
                let values: Vec<V<'_>> = bioseqs
                    .iter()
                    .map(|bioseq| V::new(named::BIOSEQ, bioseq))
                    .collect();
                let text =
                    super::queries::render_subjects(&values, parse_deflines).map_err(|error| {
                        match error {
                            RenderError::Exception(exception) => caught(&exception),
                            RenderError::Unsupported(part) => unported(&part, program),
                        }
                    })?;
                Some(std::sync::Arc::new(text))
            }
            _ => None,
        };
        let task = self.imported.program().task_name();
        Ok(ImportRun {
            queries: self.imported.queries.clone().map(std::sync::Arc::new),
            subjects,
            source: None,
            warnings: self.imported.warnings.clone(),
            option_warnings: ignored_option_warnings(
                &self.arguments,
                self.imported.program(),
                task,
            ),
            app: self.app,
            searchsp_from_strategy: false,
            query_range: None,
            recorded: self
                .imported
                .options
                .options
                .setters
                .iter()
                .flat_map(recorded_params)
                .collect(),
        })
    }

    /// The strategy's database name, when its subject is a database the command line does not
    /// override.
    pub(crate) fn database(&self) -> Option<&str> {
        match &self.imported.subject {
            Some(ImportedSubject::Database(name)) => Some(name),
            _ => None,
        }
    }

    /// Posts the import's warnings for an error that `SetOptionsForSavedStrategy` raises before
    /// the run's diagnostics exist.
    pub(crate) fn post_warnings(&self) {
        use std::io::Write;
        let _ = std::io::stderr().write_all(&ImportRun::lines(&self.imported.warnings, self.app));
    }

    /// The strategy's taxid lists over its database (`s_ImportDatabase`, blast_app_util.cpp:468-495):
    /// the stored lists are the taxids themselves (the command line's expansion happened at the
    /// export), so a list is `-taxids` or `-negative_taxids` with `-no_taxid_expansion`. A
    /// strategy with both lists is not modelled.
    pub(crate) fn db_filter(
        &self,
        program: &str,
    ) -> Result<crate::blastinput::db_filter::DbFilterArgs> {
        let options = &self.imported.options.options;
        let join = |list: &[i32]| {
            list.iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(",")
        };
        let mut filter = crate::blastinput::db_filter::DbFilterArgs::default();
        // An empty list filters nothing (probe `E.blastp.taxid_empty`: the bytes of the strategy
        // without it).
        fn some(list: Option<&[i32]>) -> Option<&[i32]> {
            list.filter(|list| !list.is_empty())
        }
        match (
            some(options.taxid_list.as_deref()),
            some(options.negative_taxid_list.as_deref()),
        ) {
            (None, None) => return Ok(filter),
            (Some(list), None) => filter.taxids = Some(join(list)),
            (None, Some(list)) => filter.negative_taxids = Some(join(list)),
            (Some(_), Some(_)) => {
                return Err(unported("a search strategy with both taxid lists", program))
            }
        }
        filter.no_taxid_expansion = true;
        Ok(filter)
    }

    /// The strategy's database mask (`s_ImportDatabase`, blast_app_util.cpp:497-511): the
    /// algorithm's key before its id, and the masking type (soft when the strategy has none).
    pub(crate) fn db_mask(&self) -> Option<(DbMaskAlgorithm, i32)> {
        let options = &self.imported.options;
        let kind = if options.subject_masking_type == 0 {
            1
        } else {
            options.subject_masking_type
        };
        if !options.filtering_key.is_empty() {
            Some((DbMaskAlgorithm::Key(options.filtering_key.clone()), kind))
        } else if options.filtering_id != -1 {
            Some((DbMaskAlgorithm::Id(options.filtering_id), kind))
        } else {
            None
        }
    }

    /// LOSAT's refusal of a part of a strategy that the program does not take.
    pub(crate) fn unported(what: &str, program: &str) -> anyhow::Error {
        unported(what, program)
    }

    /// Whether the strategy's subject is sequences the command line does not override.
    pub(crate) fn has_subjects(&self) -> bool {
        matches!(self.imported.subject, Some(ImportedSubject::Sequences(_)))
    }

    /// The query range of the strategy (`strategy.GetQueryRange()`, set by
    /// `s_ImportSearchStrategy`, blast_app_util.cpp:594-598, unless empty), as the closed range
    /// of `CQueryOptionsArgs`: `None` when the command line gives `-query_loc` (its
    /// `ExtractAlgorithmOptions` sets the range again) or the strategy has none. A range without
    /// a start or an end keeps `kEmptyRange`'s bound, which `record_interval` reads as NCBI's
    /// reader does; a range whose end comes before its start (NCBI's `Invalid sequence range` at
    /// the first query) is not modelled.
    pub(crate) fn query_range(
        &self,
        given: Option<&String>,
        program: &str,
    ) -> Result<Option<crate::blastinput::seq_range::SequenceRange>> {
        use crate::blastinput::seq_range::SequenceRange;
        let Some(range) = self.imported.query_range() else {
            return Ok(None);
        };
        if given.is_some() {
            return Ok(None);
        }
        let range = SequenceRange {
            from: range.from as usize,
            to: range.to_open.wrapping_sub(1) as usize,
        };
        let from = if range.from == SequenceRange::EMPTY_FROM {
            0
        } else {
            range.from
        };
        let to = if range.to == SequenceRange::EMPTY_TO {
            0
        } else {
            range.to
        };
        if to > 0 && to < from {
            return Err(unported(
                "a search strategy whose query range ends before its start",
                program,
            ));
        }
        Ok(Some(range))
    }

    /// LOSAT's refusal of a setter that the program does not take from a strategy.
    pub(crate) fn unported_setter(setter: &super::options::Setter, program: &str) -> anyhow::Error {
        unported(
            &format!("a search strategy with the option {setter:?}"),
            program,
        )
    }

    /// LOSAT's refusal of a strategy of another program.
    /// A strategy of another program than the running one (`program`, LOSAT's name): the run of
    /// the strategy's program with the arguments the command line keeps (`carried`), where
    /// NCBI's output is that program's import ([`RunAs`]), else the refusal.
    pub(crate) fn run_as(&self, program: &str, carried: impl FnOnce() -> Carried) -> anyhow::Error {
        use EProgram::{Blastn, Blastx, DiscMegablast, Megablast, Tblastx};
        let target = self.imported.program();
        let command = match (self.app, target) {
            ("blastn" | "blastx", Tblastx) => {
                crate::algorithm::tblastx::archive::from_carried(carried()).map(Commands::Tblastx)
            }
            ("tblastx", Blastx) => {
                crate::algorithm::blastx::archive::from_carried(carried()).map(Commands::Blastx)
            }
            ("tblastx", Blastn | Megablast | DiscMegablast) => {
                crate::algorithm::blastn::archive::from_carried(carried()).map(Commands::Blastn)
            }
            _ => return self.unported_program(program),
        };
        match command {
            Ok(command) => anyhow::Error::new(RunAs(Box::new(command))),
            Err(error) => error,
        }
    }

    pub(crate) fn unported_program(&self, program: &str) -> anyhow::Error {
        unported(
            &format!(
                "a search strategy of {} (one program running another's strategy)",
                self.imported.program().task_name()
            ),
            program,
        )
    }
}

/// A database mask algorithm of a strategy (`SetFilteringAlgorithm`'s key or id).
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum DbMaskAlgorithm {
    Id(i32),
    Key(String),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2908
/// The formatting handler's early return for an imported hit list size above 500 (no
/// `-num_descriptions`, no `-num_alignments`, a pairwise format): it reads none of the options
/// after it, so `-line_length`, `-sorthits` and `-sorthsps` keep their defaults (probe
/// `I.blastn.mts600.ll`); the report counts follow it (`report_counts_from`).
pub(crate) fn formatting_early_return(
    formatting: &mut crate::blastinput::report_counts::FormattingArgs,
    num_descriptions_given: bool,
    num_alignments_given: bool,
    outfmt: &str,
) {
    let pairwise = crate::blastinput::app::parse_formatting_string(outfmt)
        .is_ok_and(|format| format.number <= 4);
    if formatting.hitlist_start.is_some_and(|size| size > 500)
        && !num_descriptions_given
        && !num_alignments_given
        && pairwise
    {
        formatting.line_length = None;
        formatting.sort_hits = None;
        formatting.sort_hsps = None;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:759-791,841-874
/// The DUST or SEG options as a strategy's setters leave them, in the file's order:
/// `SetDustFiltering`/`SetSegFiltering` drop the options and, when true, make them again with
/// their defaults; a parameter's setter makes them with their defaults when absent and sets its
/// value. Untouched (`None`), the task's default stays.
#[derive(Clone, Copy, Debug)]
pub(crate) struct FilterSetters<T: Copy> {
    defaults: T,
    state: Option<Option<T>>,
}

impl<T: Copy> FilterSetters<T> {
    fn on(&mut self, on: bool) {
        self.state = Some(on.then_some(self.defaults));
    }

    fn set(&mut self, change: impl FnOnce(&mut T)) {
        let mut options = self.state.flatten().unwrap_or(self.defaults);
        change(&mut options);
        self.state = Some(Some(options));
    }
}

/// DUST's level, window and linker (`SDustOptionsNew`: 20, 64, 1).
pub(crate) type DustSetters = FilterSetters<(i32, i32, i32)>;
/// SEG's window, low and high cut (`SSegOptionsNew`: 12, 2.2, 2.5).
pub(crate) type SegSetters = FilterSetters<(i32, f64, f64)>;

impl Default for DustSetters {
    fn default() -> Self {
        FilterSetters {
            defaults: (20, 64, 1),
            state: None,
        }
    }
}

impl Default for SegSetters {
    fn default() -> Self {
        FilterSetters {
            defaults: (12, 2.2, 2.5),
            state: None,
        }
    }
}

impl DustSetters {
    /// Takes a DUST setter (`true`), or leaves any other.
    pub(crate) fn take(&mut self, setter: &super::options::Setter) -> bool {
        use super::options::Setter;
        match *setter {
            Setter::DustFiltering(on) => self.on(on),
            Setter::DustFilteringLevel(n) => self.set(|o| o.0 = n),
            Setter::DustFilteringWindow(n) => self.set(|o| o.1 = n),
            Setter::DustFilteringLinker(n) => self.set(|o| o.2 = n),
            _ => return false,
        }
        true
    }

    /// The `-dust` value of the options the setters left, `None` when untouched.
    pub(crate) fn value(&self) -> Option<String> {
        self.state.map(|state| match state {
            None => "no".to_string(),
            Some((level, window, linker)) => format!("{level} {window} {linker}"),
        })
    }
}

impl SegSetters {
    /// Takes a SEG setter (`true`), or leaves any other.
    pub(crate) fn take(&mut self, setter: &super::options::Setter) -> bool {
        use super::options::Setter;
        match *setter {
            Setter::SegFiltering(on) => self.on(on),
            Setter::SegFilteringWindow(n) => self.set(|o| o.0 = n),
            Setter::SegFilteringLocut(v) => self.set(|o| o.1 = v),
            Setter::SegFilteringHicut(v) => self.set(|o| o.2 = v),
            _ => return false,
        }
        true
    }

    /// The `-seg` value of the options the setters left, `None` when untouched.
    pub(crate) fn value(&self) -> Option<String> {
        self.state.map(|state| match state {
            None => "no".to_string(),
            Some((window, locut, hicut)) => format!("{window} {locut} {hicut}"),
        })
    }
}

/// The parameters a setter of the builder records in the options handle (`CBlastOptions`'s
/// `Set*` calls `m_Remote->SetValue` with the field of the value it takes): the archive of an
/// imported run lists them where NCBI's handle does. `StrandOption`, whose `strand-type` value
/// LOSAT never writes, records nothing here (no program takes it from a strategy).
fn recorded_params(setter: &Setter) -> Vec<Param> {
    use crate::report::archive::fields as f;
    use ParamValue::{BigInteger, Boolean, Evalue, Integer, Real, Str};
    let one = |field, value| vec![Param::new(field, value)];
    match setter {
        Setter::BestHitScoreEdge(v) => one(f::BEST_HIT_SCORE_EDGE, Real(*v)),
        Setter::BestHitOverhang(v) => one(f::BEST_HIT_OVERHANG, Real(*v)),
        Setter::CompositionBasedStats(v) => one(f::COMPOSITION_BASED_STATS, Integer(*v)),
        Setter::CullingLimit(v) => one(f::CULLING_LIMIT, Integer(*v)),
        Setter::DbGeneticCode(v) => one(f::DB_GENETIC_CODE, Integer(*v)),
        Setter::DbLength(v) => one(f::DB_LENGTH, BigInteger(*v)),
        Setter::DustFiltering(v) => one(f::DUST_FILTERING, Boolean(*v)),
        Setter::DustFilteringLevel(v) => one(f::DUST_FILTERING_LEVEL, Integer(*v)),
        Setter::DustFilteringWindow(v) => one(f::DUST_FILTERING_WINDOW, Integer(*v)),
        Setter::DustFilteringLinker(v) => one(f::DUST_FILTERING_LINKER, Integer(*v)),
        Setter::DomainInclusionThreshold(v) => one(f::DOMAIN_INCLUSION_THRESHOLD, Real(*v)),
        Setter::EffectiveSearchSpace(v) => one(f::EFFECTIVE_SEARCH_SPACE, BigInteger(*v)),
        Setter::EvalueThreshold(v) => one(f::EVALUE_THRESHOLD, Evalue(*v)),
        Setter::FilterString(v) => one(f::FILTER_STRING, Str(v.clone())),
        Setter::GapExtensionCost(v) => one(f::GAP_EXTENSION_COST, Integer(*v)),
        Setter::GapOpeningCost(v) => one(f::GAP_OPENING_COST, Integer(*v)),
        Setter::GapTracebackAlgorithm(v) => one(f::GAP_TRACEBACK_ALGORITHM, Integer(*v)),
        Setter::GapTrigger(v) => one(f::GAP_TRIGGER, Real(*v)),
        Setter::GapXDropoff(v) => one(f::GAP_X_DROPOFF, Real(*v)),
        Setter::GapXDropoffFinal(v) => one(f::GAP_X_DROPOFF_FINAL, Real(*v)),
        Setter::GapExtnAlgorithm(v) => one(f::GAP_EXTN_ALGORITHM, Integer(*v)),
        Setter::HitlistSize(v) => one(f::HITLIST_SIZE, Integer(*v)),
        Setter::InclusionThreshold(v) => one(f::INCLUSION_THRESHOLD, Real(*v)),
        Setter::IgnoreMsaMaster(v) => one(f::IGNORE_MSA_MASTER, Boolean(*v)),
        Setter::LongestIntronLength(v) => one(f::LONGEST_INTRON_LENGTH, Integer(*v)),
        Setter::MbTemplateLength(v) => one(f::MB_TEMPLATE_LENGTH, Integer(i32::from(*v))),
        Setter::MbTemplateType(v) => one(f::MB_TEMPLATE_TYPE, Integer(i32::from(*v))),
        Setter::MatchReward(v) => one(f::MATCH_REWARD, Integer(*v)),
        Setter::MatrixName(v) => one(f::MATRIX_NAME, Str(v.clone())),
        Setter::MismatchPenalty(v) => one(f::MISMATCH_PENALTY, Integer(*v)),
        Setter::MaskAtHash(v) => one(f::MASK_AT_HASH, Boolean(*v)),
        Setter::MaxHspsPerSubject(v) => one(f::MAX_HSPS_PER_SUBJECT, Integer(*v)),
        Setter::OutOfFrameMode(v) => one(f::OUT_OF_FRAME_MODE, Boolean(*v)),
        Setter::PhiPattern { pattern, .. } => one(f::PHI_PATTERN, Str(pattern.clone())),
        Setter::PercentIdentity(v) => one(f::PERCENT_IDENTITY, Real(*v)),
        Setter::PseudoCount(v) => one(f::PSEUDO_COUNT, Integer(*v)),
        Setter::QueryGeneticCode(v) => one(f::QUERY_GENETIC_CODE, Integer(*v)),
        Setter::QueryCovHspPerc(v) => one(f::QUERY_COV_HSP_PERC, Real(*v)),
        Setter::RepeatFiltering(v) => one(f::REPEAT_FILTERING, Boolean(*v)),
        Setter::RepeatFilteringDb(v) => one(f::REPEAT_FILTERING_DB, Str(v.clone())),
        Setter::StrandOption(_) => Vec::new(),
        Setter::SegFiltering(v) => one(f::SEG_FILTERING, Boolean(*v)),
        Setter::SegFilteringWindow(v) => one(f::SEG_FILTERING_WINDOW, Integer(*v)),
        Setter::SegFilteringLocut(v) => one(f::SEG_FILTERING_LOCUT, Real(*v)),
        Setter::SegFilteringHicut(v) => one(f::SEG_FILTERING_HICUT, Real(*v)),
        Setter::SumStatisticsMode(v) => one(f::SUM_STATISTICS_MODE, Boolean(*v)),
        Setter::SmithWatermanMode(v) => one(f::SMITH_WATERMAN_MODE, Boolean(*v)),
        // `SetGappedMode(!v)` records `UngappedMode`, the value read.
        Setter::GappedMode(v) => one(f::GAPPED_MODE, Boolean(!*v)),
        Setter::UnifiedP(v) => one(f::UNIFIED_P, Integer(*v)),
        Setter::SubjectBestHit => one(f::SUBJECT_BEST_HIT, Boolean(true)),
        Setter::WindowMaskerTaxId(v) => one(f::WINDOW_MASKER_TAX_ID, Integer(*v)),
        Setter::WindowSize(v) => one(f::WINDOW_SIZE, Integer(*v)),
        Setter::WordSize(v) => one(f::WORD_SIZE, Integer(*v)),
        Setter::WordThreshold(v) => one(f::WORD_THRESHOLD, Integer(*v)),
        // `SetUseIndex(true, name, true)`, as `CMbIndexArgs` records it (`options::Class::MbIndex`).
        Setter::UseIndex { index_name } => vec![
            Param::new(f::FORCE_MB_INDEX, Boolean(false)),
            Param::new(f::MB_INDEX_NAME, Str(index_name.clone())),
        ],
    }
}

/// One program running another program's strategy (survey F2 2.6): the strategy decides the
/// engine, its options and the subjects; the command line's program its warnings (their prefix
/// and NCBI's argument table, [`Start::app`]). NCBI's output of blastn and blastx running a
/// TBLASTX strategy, and of tblastx running a BLASTX or BLASTN one, is that program's import
/// with the running program's warnings, and tblastx's report of a BLASTN strategy cites the
/// Gapped BLAST paper (probes `C.blastn.from.tblastx.db`, `C.blastx.from.tblastx.db`,
/// `C.tblastx.from.blastx.db`, `C.tblastx.from.blastn.db`). The error that carries the
/// strategy's program's arguments up to the dispatch (`cli::main_entry`), which runs them; no
/// byte is written before it.
#[derive(Debug)]
pub(crate) struct RunAs(pub(crate) Box<Commands>);

impl std::fmt::Display for RunAs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a search strategy of another program")
    }
}

impl std::error::Error for RunAs {}

/// The arguments that a run of an imported strategy keeps from the command line (the classes
/// `SetOptionsForSavedStrategy` extracts again: query, formatting, threads, database or
/// subjects), as one program hands them to the program whose strategy it runs ([`RunAs`]).
pub(crate) struct Carried {
    /// The running program's NCBI name ([`Start::app`]).
    pub(crate) app: &'static str,
    pub(crate) path: std::path::PathBuf,
    /// The strategy's bytes of a caller without a file ([`ImportRun::source`]).
    pub(crate) source: Option<std::sync::Arc<Vec<u8>>>,
    pub(crate) query: std::path::PathBuf,
    pub(crate) query_loc: Option<String>,
    /// `-strand` (`both` when not given).
    pub(crate) strand: String,
    pub(crate) lcase_masking: bool,
    pub(crate) parse_deflines: bool,
    pub(crate) id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    pub(crate) out: Option<std::path::PathBuf>,
    pub(crate) outfmt: String,
    pub(crate) show_gis: bool,
    pub(crate) num_descriptions: Option<usize>,
    pub(crate) num_alignments: Option<usize>,
    /// `-max_target_seqs` as given (blastx keeps values below 1 for its own error).
    pub(crate) max_target_seqs: Option<i64>,
    pub(crate) formatting: crate::blastinput::report_counts::FormattingArgs,
    pub(crate) num_threads: usize,
    pub(crate) mt_mode: Option<u8>,
    pub(crate) db: Option<String>,
    pub(crate) subject: Option<std::path::PathBuf>,
    pub(crate) subject_loc: Option<String>,
    pub(crate) db_filter: crate::blastinput::db_filter::DbFilterArgs,
    pub(crate) dbsize: Option<i64>,
}

impl Carried {
    /// The import run state that names the running program for the strategy's program
    /// (`ImportRun::app`, read by its `apply_strategy`).
    pub(crate) fn import_app(&self) -> Option<ImportRun> {
        Some(ImportRun {
            app: self.app,
            source: self.source.clone(),
            ..ImportRun::default()
        })
    }
}
