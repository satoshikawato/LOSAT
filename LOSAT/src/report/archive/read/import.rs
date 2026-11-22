//! `-import_search_strategy <file>`: what NCBI's `RecoverSearchStrategy` and
//! `s_ImportSearchStrategy` make of a strategy file before the running program takes it as its
//! own options and inputs (`SetOptionsForSavedStrategy`): [`apply`], the one import of the
//! command line, the library API and LOSAT Web, over each program's table of the strategy's
//! setters ([`StrategyArgs`]).
//!
//! The file is read as `ExtractBlast4Request` reads it: first as a
//! `Blast4-get-search-strategy-reply`, then, from the start of the input, as a `Blast4-request`;
//! a serial error of either read is NCBI's error that the strategy cannot be read
//! (`failed_to_read`). The options are the strategy's (`read::options::import_request`); the
//! queries are the FASTA text NCBI writes and reads back (`read::queries`), handed to the program
//! as an unnamed temporary file (NCBI's `CTmpFile`, [`StrategyText`]); the subject is the
//! strategy's database or its Bioseqs.
//!
//! LOSAT refuses, as not yet ported, what it has no reader for: binary ASN.1 and XML strategies,
//! `seq-loc-list` queries and subjects, PSSM queries in the programs that take them, and GI lists.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:532-628
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3837
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2257-2331

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

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:594-598
    /// `GetQueryRange()`, `None` when it is empty or `-query` overrides the queries (the import
    /// sets the range only with the strategy's queries).
    pub(crate) fn query_range(&self) -> Option<SeqRange> {
        let range = self.options.options.query_range;
        (self.queries.is_some() && range != SeqRange::EMPTY).then_some(range)
    }
}

/// `s_ImportSearchStrategy`'s `CInputException` for a file that does not read as a request
/// (its message in LOSAT's words is NCBI's).
fn failed_to_read() -> anyhow::Error {
    crate::blastinput::app::options_error("Failed to read search strategy")
}

/// LOSAT's refusal of a part of a strategy that it does not read yet.
fn unported(what: &str, program: &str) -> anyhow::Error {
    anyhow::anyhow!("{what} is not supported by LOSAT's {program}")
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-267
/// What `CATCH_ALL` reports of an exception the import throws.
fn caught(exception: &NcbiException) -> anyhow::Error {
    NativeError::exception(exception).into()
}

/// The error of a read of `ExtractBlast4Request`'s second attempt: a `CSerialException` is the
/// import's read failure (`failed_to_read`); any other exception leaves the import for
/// `CATCH_ALL`.
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
        // `s_ImportSearchStrategy` reports its read failure (`failed_to_read`; a text file with a
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
        // NCBI's switch over the guessed format throws a serial exception for a format it does
        // not read, which the import reports as its read failure.
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
        Queries::Pssm => {
            if matches!(program, "BLASTN" | "BLASTP" | "BLASTX" | "TBLASTX") {
                Err(crate::blastinput::app::engine_error(&format!(
                    "PSSM found in saved strategy, but not supported for {}",
                    program_type_name(options.options.program)
                )))
            } else {
                Err(unported("a search strategy whose query is a PSSM", program))
            }
        }
        Queries::SeqLocList => Err(unported(
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:983-986
/// The arguments of the running program's command line as NCBI's `CArgs` holds them: every
/// argument given, with its text, and every argument with a default that is not given, with the
/// default's text.
pub(crate) struct NcbiArguments {
    /// `(name, text)`, in `CArgs::GetAll`'s order (by name, as the descriptions compare their
    /// names).
    all: Vec<(String, String)>,
}

impl NcbiArguments {
    /// The arguments of `program` (NCBI's name, `blastn` ...) given as `given` (name and text, a
    /// flag's text empty), over the defaults of NCBI's argument descriptions (the tables of the
    /// program's `-xmlhelp` XML, `ncbi_app::args::arguments`).
    pub(crate) fn new(program: &str, given: &[(String, String)]) -> Self {
        let mut all: std::collections::BTreeMap<String, String> =
            crate::ncbi_app::args::arguments(program)
                .into_iter()
                .flat_map(crate::ncbi_app::args::Arguments::defaults)
                .map(|(name, text)| (name.to_string(), text.to_string()))
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
/// What a run that imports a strategy carries into the program's own flow (each program's
/// `import` field): the texts that stand for its query and subject inputs, and the warnings that
/// NCBI posts during the import.
#[derive(Clone, Debug, Default)]
pub(crate) struct ImportRun {
    /// The strategy's queries (the text of NCBI's temporary FASTA file), read in place of
    /// `-query` (`open_query`).
    pub(crate) queries: Option<std::sync::Arc<StrategyText>>,
    /// The strategy's subject Bioseqs as FASTA text (`queries::render_subjects`), read in place
    /// of `-subject` and named in no report (`User specified sequence set.`).
    pub(crate) subjects: Option<std::sync::Arc<StrategyText>>,
    /// `RecoverSearchStrategy`'s warnings, posted before the program's own.
    pub(crate) warnings: Vec<String>,
    /// `x_IssueWarningsForIgnoredOptions`' warnings, posted after the formatting handler's.
    pub(crate) option_warnings: Vec<String>,
    /// The running program's NCBI name (`Start::app`), the prefix of the warnings.
    pub(crate) app: &'static str,
    /// The effective search space is the strategy's `EffectiveSearchSpace`, which does not set
    /// `OLD_FSC` in the environment as the command line's `-searchsp` does.
    pub(crate) searchsp_from_strategy: bool,
    /// blastn's MegaBLAST database index of the import (`import_into`): the strategy's
    /// (`ForceMbIndex` with its name, a new-style index), or the one that the command line's
    /// `-use_index true` names over it; the name and whether the index is old-style.
    pub(crate) mb_index: Option<(String, bool)>,
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
        text: &StrategyText,
    ) -> Result<crate::blastinput::input_files::InputFile> {
        text.open(role)
    }

    /// `old_fsc` of a run: `-searchsp` sets it, the strategy's effective search space does not.
    pub(crate) fn old_fsc(import: Option<&ImportRun>, searchsp: Option<i64>) -> bool {
        (searchsp.is_some() && !import.is_some_and(|run| run.searchsp_from_strategy))
            || crate::blastinput::app::old_fsc_environment()
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:395-439
/// The strategy's queries or subjects as FASTA text, held once: in memory until a run first
/// opens it as an input, then only as the unnamed temporary file that the open wrote (NCBI writes
/// the queries to a temporary file at the import and reads the queries from it), whose handle
/// each later open reads from its start. A run holds no copy of the text beside the file, as the
/// command line's `-query` holds none beside its file. A caller that reads records from bytes
/// (the library API and LOSAT Web) takes the text itself (`take`).
#[derive(Debug)]
pub(crate) struct StrategyText(std::sync::Mutex<Held>);

/// Where a [`StrategyText`] is.
#[derive(Debug)]
enum Held {
    Text(Vec<u8>),
    File(std::fs::File),
    /// Taken by a caller that reads records (`take`).
    Taken,
}

impl StrategyText {
    pub(crate) fn new(text: Vec<u8>) -> Self {
        Self(std::sync::Mutex::new(Held::Text(text)))
    }

    /// The text as an input file of the run (`role`: `query` or `subject`), positioned at its
    /// start: the first open writes the text to the temporary file and drops it.
    pub(crate) fn open(&self, role: &str) -> Result<crate::blastinput::input_files::InputFile> {
        let mut held = self
            .0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Held::Text(text) = &*held {
            let file = crate::blastinput::input_files::temporary_input(role, text)?;
            *held = Held::File(file);
        }
        match &*held {
            Held::File(file) => {
                let mut copy = file.try_clone()?;
                copy.seek(std::io::SeekFrom::Start(0))?;
                Ok(crate::blastinput::input_files::InputFile::from(copy))
            }
            Held::Text(_) | Held::Taken => {
                anyhow::bail!("the search strategy's {role} sequences were taken as records")
            }
        }
    }

    /// The text, for a caller that reads records from it: what the strategy carries, read back
    /// from the temporary file once a run has opened it.
    pub(crate) fn take(&self) -> Result<Vec<u8>> {
        let mut held = self
            .0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        match std::mem::replace(&mut *held, Held::Taken) {
            Held::Text(text) => Ok(text),
            Held::File(file) => {
                let mut copy = file.try_clone()?;
                copy.seek(std::io::SeekFrom::Start(0))?;
                let mut text = Vec::new();
                copy.read_to_end(&mut text)?;
                *held = Held::File(file);
                Ok(text)
            }
            Held::Taken => Ok(Vec::new()),
        }
    }
}

/// The start of an import: the strategy read, and the program's defaults to apply it to.
pub(crate) struct Start {
    pub(crate) imported: Imported,
    /// The arguments of the command line as NCBI's `CArgs` holds them.
    pub(crate) arguments: NcbiArguments,
    /// The running program's NCBI name, the prefix of its warnings (`blastn`, ...): the program
    /// of the command line, also when it runs another program's strategy (`import`).
    pub(crate) app: &'static str,
}

/// `RecoverSearchStrategy` for the running program `ncbi_name` (`blastp` ...) whose command line
/// gave the arguments `given` ([`crate::cli::GivenArguments`]): the strategy, or the error after
/// the warnings posted before it (written to `diagnostics` here, as NCBI posts them before its
/// `CATCH_ALL` reports the error; the command line's are its standard error).
pub(crate) fn start(
    command: &Command<'_>,
    ncbi_name: &'static str,
    given: &[(String, String)],
    diagnostics: &mut dyn std::io::Write,
) -> Result<Start> {
    let arguments = NcbiArguments::new(ncbi_name, given);
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
    /// The run state of the import: the subjects rendered for a reader with or without
    /// `-parse_deflines`, the strategy's query range unless the command line gives `-query_loc`
    /// (`query_loc`), the queries moved in (no copy), and the warnings of the options for the
    /// strategy's program and task.
    pub(crate) fn run(
        &mut self,
        parse_deflines: bool,
        query_loc: Option<&String>,
        program: &str,
    ) -> Result<ImportRun> {
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
                Some(std::sync::Arc::new(StrategyText::new(text)))
            }
            _ => None,
        };
        let query_range = self.query_range(query_loc, program)?;
        let task = self.imported.program().task_name();
        Ok(ImportRun {
            queries: self
                .imported
                .queries
                .take()
                .map(|text| std::sync::Arc::new(StrategyText::new(text))),
            subjects,
            warnings: self.imported.warnings.clone(),
            option_warnings: ignored_option_warnings(
                &self.arguments,
                self.imported.program(),
                task,
            ),
            app: self.app,
            searchsp_from_strategy: false,
            mb_index: None,
            query_range,
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

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:468-496
    /// The strategy's taxid lists over its database (`s_ImportDatabase`):
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

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:498-511
    /// The strategy's database mask (`s_ImportDatabase`): the algorithm's key before its id, and
    /// the masking type (soft when the strategy has none).
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

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:594-598
    /// The query range of the strategy (`GetQueryRange`, which `s_ImportSearchStrategy` sets
    /// unless it is empty), as the closed range
    /// of `CQueryOptionsArgs`: `None` when the command line gives `-query_loc` (its
    /// `ExtractAlgorithmOptions` sets the range again) or the strategy has none. A range without
    /// a start or an end keeps `kEmptyRange`'s bound, which `record_interval` reads as NCBI's
    /// reader does; a range whose end comes before its start (NCBI's error at the first query)
    /// is not modelled.
    fn query_range(
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

    /// LOSAT's refusal of a setter that the program does not take from a strategy, the option
    /// named by its field in the strategy file and its value.
    pub(crate) fn unported_setter(setter: &Setter, program: &str) -> anyhow::Error {
        unported(
            &format!("a search strategy with the option {}", option_text(setter)),
            program,
        )
    }

    /// LOSAT's refusal of a strategy of another program than the running one, where NCBI's
    /// output is not that program's import.
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

/// A setter as the strategy file names it: its field and its value (`recorded_params`).
fn option_text(setter: &Setter) -> String {
    match setter {
        Setter::StrandOption(strand) => format!("StrandOption {strand}"),
        Setter::UseIndex { index_name } => format!("MbIndexName \"{index_name}\""),
        _ => match recorded_params(setter).first() {
            Some(param) => {
                let value = match &param.value {
                    ParamValue::Boolean(true) => "TRUE".to_string(),
                    ParamValue::Boolean(false) => "FALSE".to_string(),
                    ParamValue::Integer(n) => n.to_string(),
                    ParamValue::BigInteger(n) => n.to_string(),
                    ParamValue::Real(v) | ParamValue::Evalue(v) => v.to_string(),
                    ParamValue::Str(text) => format!("\"{text}\""),
                    ParamValue::IntegerList(list) => list
                        .iter()
                        .map(i32::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                };
                format!("{} {value}", param.name())
            }
            None => String::new(),
        },
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
        // The builder's `SetUseIndex(true, name, true)` records the force flag true; the run's
        // loader that uses no index sets it false again (`options::replayed_options`).
        Setter::UseIndex { index_name } => vec![
            Param::new(f::FORCE_MB_INDEX, Boolean(true)),
            Param::new(f::MB_INDEX_NAME, Str(index_name.clone())),
        ],
    }
}

/// One of the five search programs' arguments as the import of a search strategy sets them
/// (`SetOptionsForSavedStrategy`): the arguments that the program's command line keeps, and the
/// program's table of the strategy's setters over its own typed options, the options of its
/// command line. Everything else of the import is [`apply`]'s, the same for every program and
/// every caller (the command line, the library API, LOSAT Web).
pub(crate) trait StrategyArgs: Sized {
    /// The program, by LOSAT's name in its refusals.
    const PROGRAM: crate::blastinput::app::Program;
    /// NCBI's name of the program, the prefix of its warnings.
    const APP: &'static str;
    /// Whether a strategy of `program` is the program's own.
    fn runs(program: EProgram) -> bool;
    /// Whether the arguments give a strategy that is not imported yet.
    fn imports(&self) -> bool;
    /// The arguments of a command line that gives nothing (`losat <program>`).
    fn defaults() -> Result<Self>;
    /// What `SetOptionsForSavedStrategy` takes from the command line of these arguments.
    fn kept(self) -> Kept;
    /// The arguments that `kept` gives over these ones (the program's defaults).
    fn keep(&mut self, kept: Kept) -> Result<()>;
    /// The arguments that the import sets in the same way in every program.
    fn fields(&mut self) -> Fields<'_>;
    /// One search thread.
    fn one_thread(&mut self);
    /// `-task`.
    fn set_task(&mut self, task: &str);
    /// The program's table: `setter` set on the option of the command line that it stands for,
    /// `false` for a setter that the program does not take.
    fn set(&mut self, setter: &Setter) -> Result<bool>;
    /// The options of the setters that act together (`collected`), after the table, with the
    /// import's run state (`run`).
    fn finish(&mut self, collected: Collected, run: &ImportRun) -> Result<()>;
    fn into_command(self) -> Commands;
}

/// The arguments of the program that [`StrategyArgs::fields`] gives the import.
pub(crate) struct Fields<'a> {
    pub(crate) import_search_strategy: &'a mut Option<std::path::PathBuf>,
    pub(crate) db: &'a mut Option<String>,
    pub(crate) subject: &'a mut Option<std::path::PathBuf>,
    pub(crate) subject_loc: &'a mut Option<String>,
    pub(crate) db_filter: &'a mut crate::blastinput::db_filter::DbFilterArgs,
    pub(crate) db_mask: &'a mut crate::blastinput::db_mask::DbMaskArgs,
    pub(crate) dbsize: &'a mut Option<i64>,
    pub(crate) searchsp: &'a mut Option<i64>,
    pub(crate) formatting: &'a mut crate::blastinput::report_counts::FormattingArgs,
    pub(crate) num_descriptions: Option<usize>,
    pub(crate) num_alignments: Option<usize>,
    pub(crate) outfmt: &'a str,
    pub(crate) import: &'a mut Option<ImportRun>,
}

/// The arguments that the import keeps from the running program's command line: the classes
/// that `SetOptionsForSavedStrategy` extracts again (the query options, `-out`, the formatting
/// options, the threads, the database or subjects when given, the MegaBLAST index), and the
/// arguments the command line gave, for the warnings. One program hands them to the program
/// whose strategy it runs.
pub(crate) struct Kept {
    pub(crate) path: std::path::PathBuf,
    pub(crate) given: Vec<(String, String)>,
    pub(crate) query: std::path::PathBuf,
    pub(crate) query_loc: Option<String>,
    /// `-strand`, `None` for a program without it (a protein query).
    pub(crate) strand: Option<String>,
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
    /// `-num_threads` as given (blastx keeps its own type).
    pub(crate) num_threads: i64,
    pub(crate) mt_mode: Option<u8>,
    /// `CBlastDatabaseArgs`, when the command line gives `-db` or `-subject` (`HasBeenSet`).
    pub(crate) database: Option<KeptDatabase>,
    /// `-use_index true` given: NCBI's warning that the strategy's index options are overridden.
    pub(crate) use_index: bool,
    /// What `CMbIndexArgs` of the command line asks for.
    pub(crate) mb_index: MbIndexRequest,
}

/// `CBlastDatabaseArgs` of a command line that gives `-db` or `-subject`.
pub(crate) struct KeptDatabase {
    pub(crate) db: Option<String>,
    pub(crate) subject: Option<std::path::PathBuf>,
    pub(crate) subject_loc: Option<String>,
    pub(crate) db_filter: crate::blastinput::db_filter::DbFilterArgs,
    pub(crate) db_mask: crate::blastinput::db_mask::DbMaskArgs,
    pub(crate) dbsize: Option<i64>,
}

impl KeptDatabase {
    /// The database arguments of a command line, `None` when it gives neither `-db` nor
    /// `-subject`.
    pub(crate) fn given(
        db: Option<String>,
        subject: Option<std::path::PathBuf>,
        subject_loc: Option<String>,
        db_filter: crate::blastinput::db_filter::DbFilterArgs,
        db_mask: crate::blastinput::db_mask::DbMaskArgs,
        dbsize: Option<i64>,
    ) -> Option<Self> {
        (db.is_some() || subject.is_some()).then_some(Self {
            db,
            subject,
            subject_loc,
            db_filter,
            db_mask,
            dbsize,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3377-3417
/// What `CMbIndexArgs::ExtractAlgorithmOptions` of the command line does with an imported
/// strategy (`SetOptionsForSavedStrategy` runs it for blastn, whose `-use_index` always has a
/// value): with `-use_index true`, the command line's `-task` megablast and no `-subject`, it
/// names the index by `-index_name` (old style) or `-db` (new style), or throws; otherwise the
/// strategy's index options stay.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum MbIndexRequest {
    None,
    Index { name: String, old_style: bool },
    NoName,
}

impl MbIndexRequest {
    /// The request of blastn's command line.
    pub(crate) fn of(
        use_index: Option<bool>,
        task: &str,
        subject_given: bool,
        index_name: Option<&String>,
        db: Option<&String>,
    ) -> Self {
        if use_index != Some(true) || task != "megablast" || subject_given {
            return Self::None;
        }
        match (index_name, db) {
            (Some(name), _) => Self::Index {
                name: name.clone(),
                old_style: true,
            },
            (None, Some(db)) => Self::Index {
                name: db.clone(),
                old_style: false,
            },
            (None, None) => Self::NoName,
        }
    }
}

/// The repeat filter database of `SRepeatFilterOptionsNew` (`kDefaultRepeatFilterDb`).
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:246
const DEFAULT_REPEAT_FILTER_DB: &str = "repeat/repeat_9606";

/// The setters of a strategy that act together, collected in the file's order and applied by
/// the program after its table ([`StrategyArgs::finish`]).
#[derive(Default)]
pub(crate) struct Collected {
    pub(crate) dust: DustSetters,
    pub(crate) seg: SegSetters,
    pub(crate) composition: Option<i32>,
    pub(crate) unified_p: Option<i32>,
    pub(crate) extension: Option<i32>,
    pub(crate) traceback: Option<i32>,
    pub(crate) template_type: Option<u8>,
    pub(crate) template_length: Option<u8>,
    /// The repeat filter: untouched (`None`), off (`Some(None)`) or its database.
    pub(crate) repeat: Option<Option<String>>,
    /// The MegaBLAST index that the strategy forces (`ForceMbIndex`), by its name.
    pub(crate) index: Option<String>,
}

impl Collected {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:925-959
    /// Takes a setter that acts with others (`true`), or leaves it to the program's table. The
    /// repeat filter's switch drops its options and, when on, makes them again with the default
    /// database; its database setter makes them when absent and names the database.
    fn take(&mut self, setter: &Setter) -> bool {
        if self.dust.take(setter) || self.seg.take(setter) {
            return true;
        }
        match setter {
            Setter::CompositionBasedStats(n) => self.composition = Some(*n),
            Setter::UnifiedP(n) => self.unified_p = Some(*n),
            Setter::GapExtnAlgorithm(n) => self.extension = Some(*n),
            Setter::GapTracebackAlgorithm(n) => self.traceback = Some(*n),
            Setter::MbTemplateType(n) => self.template_type = Some(*n),
            Setter::MbTemplateLength(n) => self.template_length = Some(*n),
            Setter::RepeatFiltering(on) => {
                self.repeat = Some(on.then(|| DEFAULT_REPEAT_FILTER_DB.to_string()))
            }
            Setter::RepeatFilteringDb(name) => self.repeat = Some(Some(name.clone())),
            Setter::UseIndex { index_name } => self.index = Some(index_name.clone()),
            _ => return false,
        }
        true
    }

    /// The refusal of a program that has none of `option` (the strategy's field) unless the
    /// setters left it off or untouched (`off`).
    fn refuse_unless(off: bool, option: &str, program: &str) -> Result<()> {
        if off {
            Ok(())
        } else {
            Err(unported(
                &format!("a search strategy with the option {option}"),
                program,
            ))
        }
    }

    /// DUST left off or untouched (a program of protein queries).
    pub(crate) fn no_dust(&self, program: &str) -> Result<()> {
        let off = self.dust.value().is_none_or(|value| value == "no");
        Self::refuse_unless(off, "DustFiltering", program)
    }

    /// SEG left off or untouched (blastn).
    pub(crate) fn no_seg(&self, program: &str) -> Result<()> {
        let off = self.seg.value().is_none_or(|value| value == "no");
        Self::refuse_unless(off, "SegFiltering", program)
    }

    /// No composition-based statistics and no unified P-value (blastn).
    pub(crate) fn no_composition(&self, program: &str) -> Result<()> {
        Self::refuse_unless(self.composition.is_none(), "CompositionBasedStats", program)?;
        Self::refuse_unless(self.unified_p.is_none(), "UnifiedP", program)
    }

    /// The repeat filter left off or untouched (the programs other than blastn).
    pub(crate) fn no_repeat(&self, program: &str) -> Result<()> {
        let off = self.repeat.as_ref().is_none_or(Option::is_none);
        Self::refuse_unless(off, "RepeatFilteringDB", program)
    }

    /// No gapped extension or traceback algorithm, no discontiguous template (the programs other
    /// than blastn).
    pub(crate) fn no_nucleotide_algorithms(&self, program: &str) -> Result<()> {
        Self::refuse_unless(self.extension.is_none(), "GapExtnAlgorithm", program)?;
        Self::refuse_unless(self.traceback.is_none(), "GapTracebackAlgorithm", program)?;
        Self::refuse_unless(self.template_type.is_none(), "MBTemplateType", program)?;
        Self::refuse_unless(self.template_length.is_none(), "MBTemplateLength", program)
    }

    /// The `-comp_based_stats` value of the setters (`None` when untouched): a mode 0 to 3,
    /// with the unified P-value (`u`) where the program takes it (`unified`).
    pub(crate) fn composition_value(&self, unified: bool, program: &str) -> Result<Option<String>> {
        let unified_p = self.unified_p.unwrap_or(0);
        let Some(mode) = self.composition else {
            return Self::refuse_unless(unified_p == 0, "UnifiedP", program).map(|()| None);
        };
        if !(0..=3).contains(&mode) || (unified_p != 0 && !unified) {
            return Err(unported(
                &format!("a search strategy with composition-based statistics {mode} (unified P {unified_p})"),
                program,
            ));
        }
        Ok(Some(format!(
            "{mode}{}",
            if unified_p != 0 { "u" } else { "" }
        )))
    }
}

/// `-import_search_strategy` (`RecoverSearchStrategy`, `SetOptionsForSavedStrategy`): the
/// arguments of the strategy's program for the strategy file named in `command` (or its bytes,
/// `source`, for a caller without a file), which the run then searches as any other: the
/// program's defaults with the strategy's options set on them in the file's order, and from the
/// command line only what NCBI extracts again. A strategy of another program is that program's
/// arguments where NCBI's output is that program's import with the running program's warnings
/// (blastn and blastx running a TBLASTX strategy, tblastx a BLASTX or BLASTN one: probes
/// `C.blastn.from.tblastx.db`, `C.blastx.from.tblastx.db`, `C.tblastx.from.blastx.db`,
/// `C.tblastx.from.blastn.db`); other pairs are refused. Arguments without a strategy are
/// returned as they are. The one import of the command line, the library API and LOSAT Web.
pub(crate) fn apply(
    command: Commands,
    source: Option<&[u8]>,
    diagnostics: &mut dyn std::io::Write,
) -> Result<Commands> {
    use crate::algorithm::{blastn, blastp, blastx, tblastn, tblastx};
    match command {
        Commands::Blastn(args) if args.imports() => {
            import::<blastn::BlastnArgs>(args, source, diagnostics)
        }
        Commands::Blastp(args) if args.imports() => {
            import::<blastp::BlastpArgs>(args, source, diagnostics)
        }
        Commands::Blastx(args) if args.imports() => {
            import::<blastx::BlastxArgs>(args, source, diagnostics)
        }
        Commands::Tblastn(args) if args.imports() => {
            import::<tblastn::TblastnArgs>(args, source, diagnostics)
        }
        Commands::Tblastx(args) if args.imports() => {
            import::<tblastx::TblastxArgs>(args, source, diagnostics)
        }
        other => Ok(other),
    }
}

/// Whether the arguments `command` give a strategy that is not imported yet ([`apply`] imports
/// it).
pub(crate) fn imports(command: &Commands) -> bool {
    match command {
        Commands::Blastn(args) => args.imports(),
        Commands::Blastp(args) => args.imports(),
        Commands::Blastx(args) => args.imports(),
        Commands::Tblastn(args) => args.imports(),
        Commands::Tblastx(args) => args.imports(),
        _ => false,
    }
}

/// The import of the running program `A`'s arguments.
fn import<A: StrategyArgs>(
    args: A,
    source: Option<&[u8]>,
    diagnostics: &mut dyn std::io::Write,
) -> Result<Commands> {
    use crate::algorithm::{blastn::BlastnArgs, blastx::BlastxArgs, tblastx::TblastxArgs};
    use EProgram::{Blastn, Blastx, DiscMegablast, Megablast, Tblastx};
    let kept = args.kept();
    let program = A::PROGRAM.name();
    let start = start(
        &Command {
            path: &kept.path,
            bytes: source,
            program,
            override_query: kept.query.as_os_str() != "-",
            override_subject: kept.database.is_some(),
            use_index_set: kept.use_index,
        },
        A::APP,
        &kept.given,
        diagnostics,
    )?;
    let target = start.imported.program();
    if A::runs(target) {
        return import_as::<A>(start, kept, diagnostics);
    }
    match (A::APP, target) {
        ("blastn" | "blastx", Tblastx) => import_as::<TblastxArgs>(start, kept, diagnostics),
        ("tblastx", Blastx) => import_as::<BlastxArgs>(start, kept, diagnostics),
        ("tblastx", Blastn | Megablast | DiscMegablast) => {
            import_as::<BlastnArgs>(start, kept, diagnostics)
        }
        _ => {
            let _ = diagnostics.write_all(&ImportRun::lines(&start.imported.warnings, start.app));
            Err(start.unported_program(program))
        }
    }
}

/// The import as the arguments of the program `T` (the running program's or, for another
/// program's strategy, that program's). An error comes after the import's warnings, which NCBI
/// has posted by then.
fn import_as<T: StrategyArgs>(
    mut start: Start,
    kept: Kept,
    diagnostics: &mut dyn std::io::Write,
) -> Result<Commands> {
    let imported = import_into::<T>(&mut start, kept);
    if imported.is_err() {
        let _ = diagnostics.write_all(&ImportRun::lines(&start.imported.warnings, start.app));
    }
    imported
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3837
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2565-2572
fn import_into<T: StrategyArgs>(start: &mut Start, mut kept: Kept) -> Result<Commands> {
    let program = T::PROGRAM.name();
    let mut out = T::defaults()?;
    let path = kept.path.clone();
    let parse_deflines = kept.parse_deflines;
    let query_loc = kept.query_loc.clone();
    let database = kept.database.take();
    let mb_index = std::mem::replace(&mut kept.mb_index, MbIndexRequest::None);
    out.keep(kept)?;
    let mut run = start.run(parse_deflines, query_loc.as_ref(), program)?;
    let fields = out.fields();
    *fields.import_search_strategy = Some(path);
    // `CBlastDatabaseArgs` of the command line when it gives the database or the subjects, else
    // the strategy's database, its filters and its mask, or its subjects. The command line's
    // `-dbsize` is set after the strategy's options, unless the strategy's search space trumps
    // it.
    let mut cli_dbsize = None;
    let mut one_thread = false;
    match database {
        Some(given) => {
            *fields.db = given.db;
            *fields.subject = given.subject;
            *fields.subject_loc = given.subject_loc;
            *fields.db_filter = given.db_filter;
            *fields.db_mask = given.db_mask;
            cli_dbsize = given.dbsize;
        }
        None => {
            if let Some(name) = start.database() {
                *fields.db = Some(name.to_string());
                *fields.db_filter = start.db_filter(program)?;
                if let Some((algorithm, kind)) = start.db_mask() {
                    let value = match algorithm {
                        DbMaskAlgorithm::Id(id) => id.to_string(),
                        DbMaskAlgorithm::Key(key) => key,
                    };
                    match kind {
                        1 => fields.db_mask.db_soft_mask = Some(value),
                        2 => fields.db_mask.db_hard_mask = Some(value),
                        _ => {
                            return Err(unported(
                                "a search strategy with this database masking type",
                                program,
                            ))
                        }
                    }
                }
            } else if start.has_subjects() {
                // The Owner's decision F2-Q1 (2026-11-20 「Q1=C」, PD-LOSAT-NCBI-DEFECTS exception
                // 22): a strategy's subjects are searched once, as with `-num_threads 1`, whatever
                // `-num_threads` says (NCBI searches them once per thread or aborts: its thread
                // check reads only the command line's `-subject`). The one place of the rule, for
                // every program and caller.
                one_thread = true;
            }
        }
    }
    if one_thread {
        out.one_thread();
    }
    out.set_task(start.imported.program().task_name());
    let mut collected = Collected::default();
    for setter in &start.imported.options.options.setters {
        if collected.take(setter) {
            continue;
        }
        let fields = out.fields();
        match setter {
            Setter::EffectiveSearchSpace(n) => {
                *fields.searchsp = Some(*n);
                run.searchsp_from_strategy = true;
            }
            Setter::DbLength(n) => *fields.dbsize = Some(*n),
            Setter::HitlistSize(n) => fields.formatting.hitlist_start = usize::try_from(*n).ok(),
            other => {
                if !out.set(other)? {
                    return Err(Start::unported_setter(other, program));
                }
            }
        }
    }
    // `CMbIndexArgs` of the command line, after the strategy's options: the index it names, else
    // the strategy's (a new-style index by its name). The index is searched by blastn's loader
    // of a database (`blast_engine::mb_index`), which LOSAT runs for blastn's searches of a
    // database only.
    let index = match mb_index {
        MbIndexRequest::NoName => {
            return Err(crate::blastinput::app::options_error(
                "Can not deduce database index name",
            ))
        }
        MbIndexRequest::Index { name, old_style } => Some((name, old_style)),
        MbIndexRequest::None => collected.index.take().map(|name| (name, false)),
    };
    if index.is_some() && (T::APP != "blastn" || out.fields().db.is_none()) {
        return Err(unported(
            "a search strategy with a MegaBLAST database index (ForceMbIndex) outside a BLASTN search of a database",
            program,
        ));
    }
    run.mb_index = index;
    out.finish(collected, &run)?;
    let fields = out.fields();
    if fields.searchsp.is_none() {
        if let Some(size) = cli_dbsize {
            *fields.dbsize = Some(size);
        }
    }
    formatting_early_return(
        fields.formatting,
        fields.num_descriptions.is_some(),
        fields.num_alignments.is_some(),
        fields.outfmt,
    );
    *fields.import = Some(run);
    Ok(out.into_command())
}

/// The arguments of `losat <name>` with nothing else given (a program's defaults).
pub(crate) fn defaults_of(name: &str) -> Result<Commands> {
    let cli: crate::cli::Cli = crate::cli::try_parse_from(["losat", name])?;
    Ok(cli.command)
}
