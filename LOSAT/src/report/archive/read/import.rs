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
//! ```c++
//! bool
//! RecoverSearchStrategy(const CArgs& args, blast::CBlastAppArgs* cmdline_args)
//! {
//!     CNcbiIstream* in = cmdline_args->GetImportSearchStrategyStream(args);
//!     if ( !in ) {
//!         return false;
//!     }
//!     const bool is_remote_search =
//!         (args.Exist(kArgRemote) && args[kArgRemote].HasValue() && args[kArgRemote].AsBoolean());
//!     const bool override_query = (args[kArgQuery].HasValue() &&
//!                                  args[kArgQuery].AsString() != kDfltArgQuery);
//!     const bool override_subject = CBlastDatabaseArgs::HasBeenSet(args);
//!
//!     if (CMbIndexArgs::HasBeenSet(args)) {
//!     	if (args[kArgUseIndex].AsBoolean() != kDfltArgUseIndex)
//!     		ERR_POST(Warning << "Overriding megablast BLAST DB indexed options in saved strategy");
//!     }
//!
//!     s_ImportSearchStrategy(in, cmdline_args, is_remote_search, override_query,
//!                            override_subject);
//!
//!     return true;
//! }
//! ```

use std::io::{Read as IoRead, Seek};
use std::path::Path;

use anyhow::Result;

use super::options::{import_request, EProgram, ImportedOptions, SeqRange};
use super::queries::{render_queries, RenderError};
use super::request::{Queries, Request, Subject};
use crate::cli::NativeError;
use crate::ncbi_exception::NcbiException;
use crate::report::serial::asn_types::read::ReadError;
use crate::report::serial::asn_types::stream::{guess_encoding, Encoding, GuessStop};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::Value;

/// What the running program's command line says before the strategy is read.
pub(crate) struct Command<'a> {
    /// `-import_search_strategy`'s value (`-`: the standard input).
    pub(crate) path: &'a Path,
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
    pub(crate) queries: Option<std::fs::File>,
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
// ```c++
//     // First try to read a Blast4-get-search-strategy-reply...
//     CRef<CBlast4_get_search_strategy_reply> b4_ss_reply;
//     bool succeeded = false;
//     try {
//         switch (CFormatGuess().Format(in)) {
//         case CFormatGuess::eBinaryASN:
//             b4_ss_reply.Reset(new CBlast4_get_search_strategy_reply);
//             in >> MSerial_AsnBinary >> *b4_ss_reply;
//             succeeded = true;
//             break;
//         case CFormatGuess::eTextASN:
//             b4_ss_reply.Reset(new CBlast4_get_search_strategy_reply);
//             in >> MSerial_AsnText >> *b4_ss_reply;
//             succeeded = true;
//             break;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2295-2309
// ```c++
//     b4_ss_reply.Reset();
//     in.seekg(0);
//     // Go for broke and try the Blast4-request...
//     retval.Reset(new CBlast4_request);
//     switch (CFormatGuess().Format(in)) {
//     case CFormatGuess::eBinaryASN:
//         in >> MSerial_AsnBinary >> *retval;
//         break;
//     case CFormatGuess::eTextASN:
//         in >> MSerial_AsnText >> *retval;
//         break;
// ```
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
// ```c++
//     CRef<CBlast4_request> b4req;
//     try {
//         b4req = ExtractBlast4Request(*in);
//     } catch (const CSerialException&) {
//         NCBI_THROW(CInputException, eInvalidInput,
//                    "Failed to read search strategy");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:565-599
// ```c++
//     // Get the subject
//     if (override_subject) {
//         ERR_POST(Warning << "Overriding database/subject in saved strategy");
//     } else {
//         CRef<blast::CBlastDatabaseArgs> db_args;
//         CRef<CBlast4_subject> subj = strategy.GetSubject();
//         const bool subject_is_protein = Blast_SubjectIsProtein(prog) ? true : false;
//
//         if (subj->IsDatabase()) {
//             db_args = s_ImportDatabase(*subj, strategy.GetOptionsBuilder(),
//                                        subject_is_protein, is_remote_search);
//         } else {
//             db_args = s_ImportSubjects(*subj, subject_is_protein);
//         }
//         _ASSERT(db_args.NotEmpty());
//         cmdline_args->SetBlastDatabaseArgs(db_args);
//     }
//
//     // Get the query, queries, or pssm
//     if (override_query) {
//         ERR_POST(Warning << "Overriding query in saved strategy");
//     } else {
//         CRef<CBlast4_queries> queries = strategy.GetQueries();
//         if (queries->IsPssm()) {
//             s_ImportPssm(*queries, opts_hndl, cmdline_args);
//         } else {
//             s_ImportQueries(*queries, opts_hndl, cmdline_args);
//         }
// ```
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
    let input = open_strategy(command.path, program).map_err(|e| fail(&warnings, e))?;
    let value = extract_request(input, program).map_err(|e| fail(&warnings, e))?;
    let request =
        Request::from_strategy(value).map_err(|stop| fail(&warnings, stop_error(stop, program)))?;
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

/// The error of a request that `CImportStrategy` does not take.
fn stop_error(stop: super::Stop, program: &str) -> anyhow::Error {
    match stop {
        super::Stop::Exception(exception) => caught(&exception),
        other => unported(
            &format!("a search strategy that NCBI stops on ({other:?})"),
            program,
        ),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:452-465
// ```c++
//     const string dbname(subj.GetDatabase());
//     CRef<CSearchDatabase> search_db(new CSearchDatabase(dbname, mol));
//     if (opts_builder.HaveEntrezQuery()) {
//         string limit(opts_builder.GetEntrezQuery());
//         search_db->SetEntrezQueryLimitation(limit);
//         if ( !is_remote_search ) {
//             string msg("Entrez query '");
//             msg += limit + string("' will not be processed locally.\n");
//             msg += string("Please use the -remote option.");
//             throw runtime_error(msg);
//         }
//     }
// ```
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
        Subject::SeqLocList(_) => Err(unported(
            "a search strategy whose subject is a list of sequence locations",
            program,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:372-392
// ```c++
//     } else {
//         EBlastProgramType p = opts_hndl->GetOptions().GetProgramType();
//         string msg("PSSM found in saved strategy, but not supported ");
//         msg += "for " + Blast_ProgramNameFromType(p);
//         NCBI_THROW(CBlastException, eNotSupported, msg);
//     }
// ```
/// `s_ImportPssm` and `s_ImportQueries`: the queries as NCBI's temporary FASTA file.
fn import_queries(
    queries: &Queries,
    options: &ImportedOptions,
    program: &str,
) -> Result<std::fs::File> {
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
            crate::blastinput::input_files::temporary_input("query", &text).map_err(|error| {
                anyhow::anyhow!("cannot write the strategy's queries to a temporary file: {error}")
            })
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
// ```c++
// void
// CBlastAppArgs::x_IssueWarningsForIgnoredOptions(const CArgs& args)
// {
//     set<string> can_override;
//     can_override.insert(kArgQuery);
//     can_override.insert(kArgQueryLocation);
//     can_override.insert(kArgSubject);
//     can_override.insert(kArgSubjectLocation);
//     can_override.insert(kArgUseLCaseMasking);
//     can_override.insert(kArgDb);
//     can_override.insert(kArgDbSize);
//     can_override.insert(kArgEntrezQuery);
//     can_override.insert(kArgDbSoftMask);
//     can_override.insert(kArgDbHardMask);
//     can_override.insert(kArgUseIndex);
//     can_override.insert(kArgIndexName);
//     can_override.insert(kArgStrand);
//     can_override.insert(kArgParseDeflines);
//     can_override.insert(kArgOutput);
//     can_override.insert(kArgOutputFormat);
//     can_override.insert(kArgNumDescriptions);
//     can_override.insert(kArgNumAlignments);
//     can_override.insert(kArgMaxTargetSequences);
//     can_override.insert(kArgRemote);
//     can_override.insert(kArgNumThreads);
//     can_override.insert(kArgInputSearchStrategy);
//     can_override.insert(kArgRemote);
//     can_override.insert("remote_verbose");
//     can_override.insert("verbose");
// ```
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
// ```c++
//     // this stores the arguments (and their defaults) that cannot be overriden
//     map<string, string> has_defaults;
//     EBlastProgramType prog = m_OptsHandle->GetOptions().GetProgramType();
//     has_defaults[kArgCompBasedStats] =
//     		Blast_ProgramIsRpsBlast(prog) ? kDfltArgCompBasedStatsDelta:kDfltArgCompBasedStats;
//     // FIX the line below for igblast, and add igblast options
//     has_defaults[kArgEvalue] = NStr::DoubleToString(BLAST_EXPECT_VALUE);
//     has_defaults[kTask] = m_Task;
//     has_defaults[kArgOldStyleIndex] = kDfltArgOldStyleIndex;
//
//     if (Blast_QueryIsProtein(prog)) {
//         if (NStr::Find(m_Task, "blastp") != NPOS ||
//             NStr::Find(m_Task, "psiblast") != NPOS) {
//             has_defaults[kArgSegFiltering] = kDfltArgNoFiltering;
//         } else {
//             has_defaults[kArgSegFiltering] = kDfltArgSegFiltering;
//         }
//         has_defaults[kArgLookupTableMaskingOnly] =
//             kDfltArgLookupTableMaskingOnlyProt;
//         has_defaults[kArgGapTrigger] =
//             NStr::IntToString((int)BLAST_GAP_TRIGGER_PROT);
//     } else {
//         has_defaults[kArgDustFiltering] = kDfltArgDustFiltering;
//         has_defaults[kArgLookupTableMaskingOnly] =
//             kDfltArgLookupTableMaskingOnlyNucl;
//         has_defaults[kArgGapTrigger] =
//             NStr::IntToString((int)BLAST_GAP_TRIGGER_NUCL);
//     }
//     has_defaults[kArgOffDiagonalRange] =
//         NStr::IntToString(kDfltOffDiagonalRange);
//     has_defaults[kArgMaskLevel] = kDfltArgMaskLevel;
//     has_defaults[kArgMaxIntronLength] =
//         NStr::IntToString(kDfltArgMaxIntronLength);
//     has_defaults[kArgQueryGeneticCode] = NStr::IntToString((int)BLAST_GENETIC_CODE);
//     has_defaults[kArgDbGeneticCode] = NStr::IntToString((int)BLAST_GENETIC_CODE);
//     // pssm engine/psiblast default options
//     has_defaults[kArgPSIPseudocount] =
//         NStr::IntToString(PSI_PSEUDO_COUNT_CONST);
//     has_defaults[kArgPSIInclusionEThreshold] =
//         NStr::DoubleToString(PSI_INCLUSION_ETHRESH);
//     has_defaults[kArgPSINumIterations] =
//         NStr::IntToString(kDfltArgPSINumIterations);
// ```
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
// ```c++
//     typedef vector< CRef<CArgValue> > TArgs;
//     TArgs arguments = args.GetAll();
//     ITERATE(TArgs, a, arguments) {
//         const string& arg_name = (*a)->GetName();
//         const string& arg_value = (*a)->AsString();
//         // if it has a default value, ignore it if it's not different from the
//         // default, otherwise, issue a warning
//         if (has_defaults.find(arg_name) != has_defaults.end()) {
//             if (has_defaults[arg_name] == arg_value) {
//                 continue;
//             } else {
//                 if (arg_name == kTask && arg_value == "megablast") {
//                     // No need to issue warning here, as it's OK to change this
//                     continue;
//                 }
//                 ERR_POST(Warning << arg_name << " cannot be overridden when "
//                          "using a search strategy");
//             }
//         }
//         // if the argument cannot be overridden, issue a warning
//         if (can_override.find(arg_name) == can_override.end()) {
//             ERR_POST(Warning << arg_name << " cannot be overridden when "
//                      "using a search strategy");
//         }
//     }
// ```
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
    /// The strategy's queries (NCBI's temporary FASTA file), read in place of `-query`.
    pub(crate) queries: Option<std::sync::Arc<std::fs::File>>,
    /// The strategy's subject Bioseqs as FASTA text (`queries::render_subjects`), read in place
    /// of `-subject` and named in no report (`User specified sequence set.`).
    pub(crate) subjects: Option<std::sync::Arc<std::fs::File>>,
    /// `RecoverSearchStrategy`'s warnings, posted before the program's own.
    pub(crate) warnings: Vec<String>,
    /// `x_IssueWarningsForIgnoredOptions`' warnings, posted after the formatting handler's.
    pub(crate) option_warnings: Vec<String>,
    /// The effective search space is the strategy's `EffectiveSearchSpace`, which does not set
    /// `OLD_FSC` as `-searchsp` does (blast_args.cpp:303-307).
    pub(crate) searchsp_from_strategy: bool,
    /// The strategy's query range (`Start::query_range`), which the program's query options
    /// take in place of `-query_loc`'s.
    pub(crate) query_range: Option<crate::blastinput::seq_range::SequenceRange>,
}

impl ImportRun {
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
            Some(file) => Ok(file.try_clone()?.into()),
            None => crate::blastinput::input_files::open_input(path, "query", program),
        }
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
}

/// `RecoverSearchStrategy` for the running program `ncbi_name` (`blastp` ...): the strategy, or
/// the error after the warnings posted before it (written to the standard error here, as NCBI
/// posts them before its `CATCH_ALL` reports the error).
pub(crate) fn start(command: &Command<'_>, ncbi_name: &str) -> Result<Start> {
    let arguments = NcbiArguments::new(ncbi_name, &crate::cli::given_arguments());
    match recover(command) {
        Ok(imported) => Ok(Start {
            imported,
            arguments,
        }),
        Err(failure) => {
            use std::io::Write;
            let _ = std::io::stderr().write_all(&ImportRun::lines(&failure.warnings, ncbi_name));
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
                Some(std::sync::Arc::new(
                    crate::blastinput::input_files::temporary_input("subject", &text)?,
                ))
            }
            _ => None,
        };
        let task = self.imported.program().task_name();
        Ok(ImportRun {
            queries: match &self.imported.queries {
                Some(file) => Some(std::sync::Arc::new(file.try_clone()?)),
                None => None,
            },
            subjects,
            warnings: self.imported.warnings.clone(),
            option_warnings: ignored_option_warnings(
                &self.arguments,
                self.imported.program(),
                task,
            ),
            searchsp_from_strategy: false,
            query_range: None,
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
