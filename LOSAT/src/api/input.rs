//! The inputs of a search: FASTA records in memory, read by the program's own reader.
//!
//! The command line reads its `-query` and `-subject` files with NCBI's FASTA reader
//! (`CBlastFastaInputSource`, ported in `blastinput::fasta_reader`); the library API reads
//! bytes with the same reader and the same flags, so a search of the records it returns is the
//! search of the files that hold those bytes (LOSAT Web ABI v2 `register`, docs/web/abi_v2.md
//! §4, reads them the same way).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:316-361
//! ```c++
//! void
//! CBlastFastaInputSource::x_InitInputReader()
//! {
//!     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
//!                                     CFastaReader::fParseRawID:
//!                                     (CFastaReader::fNoParseID |
//!                                      CFastaReader::fDLOptional);
//!
//!     // Allow CFastaReader fSkipCheck flag to be set based
//!     // on new CBlastInputSourceConfig property - GetSkipSeqCheck() -RMH-
//!     flags += ( m_Config.GetSkipSeqCheck() ? CFastaReader::fSkipCheck : 0 );
//!
//!     flags += (m_ReadProteins
//!               ? CFastaReader::fAssumeProt
//!               : CFastaReader::fAssumeNuc);
//! ```
#![warn(missing_docs, unused)]

use std::io::Read;

use super::error::Error;
use super::program::Program;
use super::search::Search;
use crate::blastinput::fasta_reader::{
    read_all, read_queries, FastaInputSource, FastaRecord, ParseErrorCode, QueryRecords, ReadError,
    ReaderConfig,
};
use crate::blastinput::input_files::open_file_input;
use crate::blastinput::query_batch::QueryReading;
use crate::cli::Commands;

/// Which input of a search a FASTA holds. The reader's flags and the generated local IDs
/// (`Query_N`, `Subject_N`) depend on it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Role {
    /// The `-query` input.
    Query,
    /// The `-subject` input.
    Subject,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:202-216
// ```c++
//         SDataLoaderConfig dlconfig =
//             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
//                                                    db_adapter);
//         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
//                                      query_opts->UseLowercaseMasks(),
//                                      query_opts->GetParseDeflines(),
//                                      query_opts->GetRange());
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:78-80
// ```c++
//     static const string kDataLoadersConfig("DATA_LOADERS");
//
//     if (registry.HasEntry("BLAST", kDataLoadersConfig)) {
// ```
/// Whether the program reads the input of `role` as protein (`fAssumeProt`) or nucleotide
/// (`fAssumeNuc`) sequences.
pub(super) fn is_protein(program: Program, role: Role) -> bool {
    match (program, role) {
        (Program::Blastn, _) | (Program::Tblastx, _) => false,
        (Program::Blastp | Program::Psiblast, _) => true,
        (Program::Tblastn, Role::Query) => true,
        (Program::Tblastn, Role::Subject) => false,
        (Program::Blastx, Role::Query) => false,
        (Program::Blastx, Role::Subject) => true,
    }
}

/// The reader configuration of the command line for the input of `role`: the program's flags,
/// the data loaders on (the command line without `.ncbirc` or `DATA_LOADERS`). `writes_archive`
/// is whether the search writes the BLAST archive (outfmt 11), or `None` where its formats are not
/// known yet (`read_inputs`).
fn reader_config(program: Program, role: Role, writes_archive: Option<bool>) -> ReaderConfig {
    let name = program.reader_name();
    let protein = is_protein(program, role);
    let config = match role {
        Role::Query => ReaderConfig::query(name, protein, true),
        Role::Subject => ReaderConfig::subject(name, protein, true),
    };
    match writes_archive {
        // The archive (outfmt 11) shows each record's defline as it is written, and the formats of
        // the searches of the records are not known here (`check_archive_inputs`).
        None => config.with_defline(true).with_archive_unknown(),
        // The command line's reader, which keeps the deflines for the archive and refuses as it
        // reads what the archive cannot write (`refuse_range_in_archive`).
        Some(writes) => config.with_defline(writes),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:633-640
// ```c++
//     if (has_range) {
//         auto old_id = FindBestChoice(GetIDs(), CSeq_id::BestRank);
//         // generate a new ID, and record its relation to the given one(s).
//         SetIDs().clear();
//         GenerateID();
//
//         CRef<CSeq_align> align = xCreateAlignment(old_id,
//             GetIDs().front(), range_start, range_end);
// ```
// LOSAT's archive does not write this Seq-hist alignment of a record with a range; the command
// line's reader refuses such a record as it reads it.
/// The refusal of a search that writes the BLAST archive of records with a range in their first
/// word under `-parse_deflines` (`refuse_range_in_archive`).
pub(super) fn check_archive_inputs(program: Program, inputs: Inputs<'_>) -> Result<(), Error> {
    use crate::blastinput::fasta_reader::believed::{defline_has_range, refuse_range_in_archive};
    let name = program.reader_name();
    for (role, records, protein) in [
        ("query", inputs.query, is_protein(program, Role::Query)),
        (
            "subject",
            inputs.subject,
            is_protein(program, Role::Subject),
        ),
    ] {
        let ranged = records.iter().any(|record| {
            record.believed().is_some()
                && record
                    .source
                    .as_ref()
                    .is_some_and(|source| defline_has_range(&source.defline, protein, name))
        });
        if ranged {
            return Err(Error::search(&refuse_range_in_archive(role, name)));
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-152
// ```c++
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
// ```
/// Reads every record of a FASTA input, as the command line reads the file named by `-query`
/// or `-subject` of `program`: the same reader, the same flags, the same generated local IDs.
///
/// Every record the reader returns is kept, records without residues and a first record
/// without a defline included. The reader's messages about a record stay with it
/// (`FastaRecord::warnings`) and a search writes them where the command line writes them, so
/// a search of the records is the command line's search of the file. An input of white space
/// only has no record. A reader error is the command line's message (`Error::Input`).
pub fn read_fasta(program: Program, role: Role, bytes: &[u8]) -> Result<Vec<FastaRecord>, Error> {
    read_source(
        FastaInputSource::from_bytes(bytes, reader_config(program, role, None)),
        role,
    )
}

/// Reads every record of an input source; a reader error is the command line's message: the
/// subjects end as `GetAllSeqs` ends (`ReadError::into_all_seqs_error`), the queries as a batch of
/// queries ends (`ReadError::into_app_error`, `GetNextSeqBatch` rethrows the reader's exception).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:205-214
/// ```c++
///         catch (const CObjReaderParseException& e) {
///             auto err = e.GetErrCode();
///             if (err == CObjReaderParseException::eEOF) {
///                 break;
///             } else if (err == CObjReaderParseException::eNoDefline) {
///                 CNcbiStrstream ss;
///                 ss << "Query input doesn't start with "
///                     "a defline or comment, line " << e.GetPos() << ends;
///                 NCBI_THROW(CInputException, eInvalidInput, ss.str());
///             }
/// ```
fn read_source<R: std::io::Read>(
    mut source: FastaInputSource<R>,
    role: Role,
) -> Result<Vec<FastaRecord>, Error> {
    read_all(&mut source, &mut |_| Ok(())).map_err(|error| {
        let error = match role {
            Role::Query => error.into_app_error(),
            Role::Subject => error.into_all_seqs_error(),
        };
        Error::Input(super::error::input_bytes(&error))
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
// ```c++
// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
// {
//     if ( !m_Ios ) {
//         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
//             "File is not accessible",AsString()));
// ```
/// The bytes of the file `path`, as the command line opens its `-query` or `-subject` argument:
/// the same message for a file that does not open (`Error::Input`), and a file named `.gz`,
/// `.bz2` or `.zst` is inflated (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, exception 8). `read_fasta`
/// and `read_inputs` read the bytes.
pub fn read_input_file(role: Role, path: &std::path::Path) -> Result<Vec<u8>, Error> {
    let input_error = |error: anyhow::Error| Error::Input(super::error::input_bytes(&error));
    let name = match role {
        Role::Query => "query",
        Role::Subject => "subject",
    };
    let mut file = open_file_input(path, name).map_err(input_error)?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .map_err(|error| input_error(error.into()))?;
    Ok(bytes)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:188-207
// ```c++
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                           db_adapter, scope);
//         ...
//         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
// ```
/// Reads the query and subject FASTA bytes of `search`, as the command line reads its `-query`
/// and `-subject` files: the records that `read_fasta` gives, and with `-parse_deflines` in the
/// search's argv the ids of their deflines. The subjects come first, in one scope of ids with the
/// queries, as NCBI reads them: a query and a subject with one id are one sequence to the reports.
///
/// A search with `-parse_deflines` takes its records from here; a search of a database gives no
/// subject bytes. Returns `(query, subject)`: every query record, or the reader's error
/// (`Error::Input`) where a query record does not read, before any search. [`read_input_records`]
/// keeps the query records read before that error, and the error, which a search of them meets
/// where the command line meets it: after the reports of the query batches before.
pub fn read_inputs(
    search: &Search,
    query: &[u8],
    subject: &[u8],
) -> Result<(Vec<FastaRecord>, Vec<FastaRecord>), Error> {
    let records = read_input_records(search, query, subject)?;
    match records.query_end {
        QueryEnd::Error(error) => Err(Error::Input(error.message)),
        QueryEnd::Input | QueryEnd::BlankLines => Ok((records.query, records.subject)),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
/// Reads the query and subject FASTA bytes of `search` as [`read_inputs`] does, and keeps how the
/// reading of the query ended ([`QueryEnd`]): the command line reads its queries one batch at a
/// time, so a query record that its reader cannot read ends the run in the batch that reads it,
/// after the reports of the batches before (line WQA's D10, 2026-11-19). The records are those
/// read before the end; a search of [`InputRecords::inputs`] (or of [`Inputs::query_end`]) meets
/// the end where the command line meets it. A subject's error fails here, as the command line
/// reads every subject before the first query.
pub fn read_input_records(
    search: &Search,
    query: &[u8],
    subject: &[u8],
) -> Result<InputRecords, Error> {
    read_search_inputs(
        search.program(),
        &search.any_subject(),
        None,
        query,
        subject,
    )
}

/// The records of a search's inputs as the command line reads its `-query` and `-subject` files
/// ([`read_input_records`]).
#[derive(Debug, Clone, Default)]
pub struct InputRecords {
    /// The query records read before the end of the query input (`query_end`).
    pub query: Vec<FastaRecord>,
    /// How the reading of the query ended after `query`.
    pub query_end: QueryEnd,
    /// The subject records.
    pub subject: Vec<FastaRecord>,
}

impl InputRecords {
    /// The query records against the subject records ([`Inputs::new`]), with the query's end.
    pub fn inputs(&self) -> Inputs<'_> {
        Inputs::new(&self.query, &self.subject).query_end(&self.query_end)
    }

    /// The query records against the BLAST database `name` ([`Inputs::database`]), with the
    /// query's end.
    pub fn database<'a>(&'a self, name: &'a str) -> Inputs<'a> {
        Inputs::database(&self.query, name).query_end(&self.query_end)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-152
/// How the reading of a query input ended after its records ([`InputRecords::query_end`]): what
/// the command line's query batches meet after the last record read. A search meets it in the
/// batch being read when the records run out before the batch reaches its size, and else in a
/// batch of its own.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum QueryEnd {
    /// The end of the input (`CBlastInput::End()`); an input of white space only has no record
    /// and is NCBI's `Query is Empty!` (`IsIStreamEmpty`).
    #[default]
    Input,
    /// Blank and comment lines where the reader expected a record (`eEOF`): the batch being read
    /// ends there. An input that has nothing else gives a first batch without a query, NCBI's
    /// `BLAST engine error: Empty CBlastQueryVector`.
    BlankLines,
    /// The reader's error on the record after the records read: the run ends with it in that
    /// record's batch.
    Error(QueryError),
}

/// The reader's error on a query record ([`QueryEnd::Error`]): the command line's message, and
/// the reader's messages about the record that it wrote before the error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryError {
    /// The command line's message without its last newline (`Error::Input`'s text).
    message: String,
    /// The reader's messages about the record, before the error.
    warnings: Vec<u8>,
    failure: ReadFailure,
}

/// A reader error as data, so that each search of the records meets the same error
/// (`QueryError::read_error`).
#[derive(Debug, Clone, PartialEq, Eq)]
enum ReadFailure {
    /// `ReadError::Parse`.
    Parse {
        code: ParseErrorCode,
        message: String,
        line: u64,
    },
    /// `ReadError::Unsupported`, LOSAT's refusal: its text, and its exit status where it is a
    /// `NativeError`.
    Refused { exit: Option<i32>, text: String },
}

impl QueryError {
    fn new(error: ReadError, warnings: Vec<u8>) -> Self {
        let failure = match &error {
            ReadError::Parse {
                code,
                message,
                line,
            } => ReadFailure::Parse {
                code: *code,
                message: message.clone(),
                line: *line,
            },
            ReadError::Unsupported(refusal) => {
                match refusal.downcast_ref::<crate::cli::NativeError>() {
                    Some(native) => ReadFailure::Refused {
                        exit: Some(native.exit),
                        text: native.message.clone(),
                    },
                    None => ReadFailure::Refused {
                        exit: None,
                        text: format!("{refusal:#}"),
                    },
                }
            }
            // The messages of the records are kept, not written: no write fails here.
            ReadError::Write(error) => ReadFailure::Refused {
                exit: None,
                text: error.to_string(),
            },
        };
        let message = format!("{:#}", error.into_app_error())
            .trim_end_matches('\n')
            .to_string();
        Self {
            message,
            warnings,
            failure,
        }
    }

    /// The command line's message, without its last newline: what [`read_inputs`] fails with.
    pub fn message(&self) -> &str {
        &self.message
    }

    /// The reader's error, as the reader gave it.
    fn read_error(&self) -> ReadError {
        match &self.failure {
            ReadFailure::Parse {
                code,
                message,
                line,
            } => ReadError::Parse {
                code: *code,
                message: message.clone(),
                line: *line,
            },
            ReadFailure::Refused {
                exit: Some(exit),
                text,
            } => ReadError::Unsupported(
                crate::cli::NativeError {
                    exit: *exit,
                    message: text.clone(),
                }
                .into(),
            ),
            ReadFailure::Refused { exit: None, text } => {
                ReadError::Unsupported(anyhow::anyhow!("{text}"))
            }
        }
    }

    /// The reader's messages about the record before the error, which the search writes before
    /// the error.
    pub(crate) fn warnings(&self) -> &[u8] {
        &self.warnings
    }

    /// The error as the programs end with it (`ReadError::into_app_error`).
    pub(crate) fn app_error(&self) -> anyhow::Error {
        self.read_error().into_app_error()
    }
}

impl std::fmt::Display for QueryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl QueryEnd {
    /// The end of the command line's reading of the queries (`fasta_reader::read_queries`).
    fn of(end: crate::blastinput::fasta_reader::QueryEnd) -> Self {
        use crate::blastinput::fasta_reader::QueryEnd as Read;
        match end {
            Read::Input => QueryEnd::Input,
            Read::EofError => QueryEnd::BlankLines,
            Read::Error { error, warnings } => QueryEnd::Error(QueryError::new(error, warnings)),
        }
    }

    /// The reading that a program's search meets after the records (`QueryReading`): the records
    /// of an input read to its end are `QueryReading::Records` (no record: `Query is Empty!`).
    pub(crate) fn reading(&self) -> QueryReading {
        match self {
            QueryEnd::Input => QueryReading::Records,
            QueryEnd::BlankLines => QueryReading::BlankLines,
            QueryEnd::Error(error) => QueryReading::Error {
                error: error.read_error(),
                warnings: error.warnings.clone(),
            },
        }
    }

    /// Whether the end is the end of the input: records without one are NCBI's `Query is Empty!`
    /// where a program checks it before its batches (`run_local` of TBLASTN and BLASTX).
    pub(crate) fn is_input(&self) -> bool {
        matches!(self, QueryEnd::Input)
    }
}

/// The query records of `bytes` as the command line reads its `-query` file with `config`, and
/// how the reading ended: NCBI's `IsIStreamEmpty` first (white space only: no record, the end of
/// the input), then the batches' reading (`fasta_reader::read_queries`), which stops at the
/// reader's first error. LOSAT Web's `register` reads a query so (`web_adapter`).
pub fn read_query_records(bytes: &[u8], config: ReaderConfig) -> (Vec<FastaRecord>, QueryEnd) {
    read_query_source(FastaInputSource::from_bytes(bytes, config), bytes)
}

/// `read_query_records` of a source made of `bytes`.
fn read_query_source<R: std::io::Read>(
    mut source: FastaInputSource<R>,
    bytes: &[u8],
) -> (Vec<FastaRecord>, QueryEnd) {
    if crate::blastinput::input_files::is_blank(bytes) {
        return (Vec::new(), QueryEnd::Input);
    }
    let QueryRecords { records, end } = read_queries(&mut source);
    (records, QueryEnd::of(end))
}

/// Whether the arguments `command` have `-parse_deflines`: a search that reads the ids of the
/// deflines, whose query and subjects are read together for its arguments (`read_search_inputs`).
pub fn parse_deflines_of(command: &Commands) -> bool {
    match command {
        Commands::Blastp(args) => args.parse_deflines,
        Commands::Tblastn(args) => args.parse_deflines,
        Commands::Blastn(args) => args.parse_deflines,
        Commands::Tblastx(args) => args.parse_deflines,
        Commands::Blastx(args) => args.parse_deflines,
        Commands::Psiblast(args) => args.parse_deflines,
        Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => false,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:67-72
// ```c++
//     // Set an appropriate default for the strand
//     if (m_Strand == eNa_strand_other) {
//         m_Strand = (m_DLConfig.m_IsLoadingProteins)
//             ? eNa_strand_unknown
//             : eNa_strand_both;
//     }
// ```
/// The `-strand` value of the arguments `command`, for the programs that have the option.
fn strand_of(command: &Commands) -> Option<&str> {
    match command {
        Commands::Blastn(args) => Some(args.strand.as_str()),
        Commands::Blastx(args) => Some(args.strand.as_str()),
        _ => None,
    }
}

/// The query and subject FASTA bytes of a search of `program` with the parsed arguments `command`
/// (the command line's, of that program), read as the command line reads its `-query` and
/// `-subject` files: the reader of `read_fasta`, with `-parse_deflines` the ids of the deflines,
/// read in one scope of ids, the subjects first, and the query's `-strand`. `writes_archive` is
/// whether the search writes the BLAST archive (outfmt 11), whose Bioseqs show the deflines as
/// read: the reader keeps them then, and refuses an id with a range as it reads it, as the command
/// line's reader does (`refuse_range_in_archive`). `None` where the formats are not known yet
/// (`read_inputs`): every defline is kept, and the search refuses such a record when it writes the
/// archive (`check_archive_inputs`). A subject's reader error is the command line's message
/// (`read_source`); the query keeps the records read before its reader's error and the error
/// (`QueryEnd`), which the search meets in its batch, as the command line's batches read the
/// queries (`read_query_records`). The one reading of `read_input_records` and LOSAT Web's
/// `register_inputs` (`web/adapter`, docs/web/abi_v2.md §4).
pub fn read_search_inputs(
    program: Program,
    command: &Commands,
    writes_archive: Option<bool>,
    query: &[u8],
    subject: &[u8],
) -> Result<InputRecords, Error> {
    let parse_deflines = parse_deflines_of(command);
    let scope = parse_deflines.then(Default::default);
    let source = |role: Role, bytes| {
        let config =
            reader_config(program, role, writes_archive).with_parse_deflines(parse_deflines);
        let config = match strand_of(command) {
            Some(strand) if role == Role::Query => config.with_strand(strand),
            _ => config,
        };
        let mut source = FastaInputSource::from_bytes(bytes, config);
        source.set_id_scope(scope.clone());
        source
    };
    let subject = read_source(source(Role::Subject, subject), Role::Subject)?;
    let (query, query_end) = read_query_source(source(Role::Query, query), query);
    Ok(InputRecords {
        query,
        query_end,
        subject,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:797-799
// ```c++
// 	    string dbname;
// 	    if (m_IsDbScan)
// 		dbname = string("User specified sequence set (Input: ") + m_SubjectTag + string(")");
// ```
/// The inputs of one search: the query and subject records, and the names that the reports
/// show for them where the command line shows its `-query` and `-subject` file names (the
/// `Database:` line of a `-subject` search, `User specified sequence set (Input: NAME)`).
///
/// The records come from [`read_fasta`] with the search's program and the input's role, or
/// from [`FastaRecord::new`]. A database subject comes from [`Inputs::database`]: the engine
/// opens it by its name, as the command line opens `-db`. A PSSM query of psiblast and tblastn
/// comes from [`Inputs::pssm`] and [`Inputs::pssm_database`]: the engine reads the file, as the
/// command line reads `-in_pssm`.
#[derive(Debug, Clone, Copy)]
pub struct Inputs<'a> {
    /// The query records (`-query`); none for a PSSM query (`pssm`).
    pub(crate) query: &'a [FastaRecord],
    /// The subject records (`-subject`); none for a database.
    pub(crate) subject: &'a [FastaRecord],
    /// What the reports show as the query's file name. The command line's default is the
    /// `-query` argument; the API's is `query`.
    pub(crate) query_name: &'a str,
    /// What the reports show as the subject set's file name, or the name of the database. The
    /// command line's default is the `-subject` argument; the API's is `subject`.
    pub(crate) subject_name: &'a str,
    /// Whether the subject is the database `subject_name` (`-db`).
    pub(crate) database: bool,
    /// The PSSM file whose sequence is the one query (`-in_pssm`), in place of query records.
    pub(crate) pssm: Option<&'a std::path::Path>,
    /// How the reading of the query records ended (`Inputs::query_end`); `None`: at the end of
    /// the input.
    pub(crate) query_end: Option<&'a QueryEnd>,
}

impl<'a> Inputs<'a> {
    /// The query and subject records, shown as `query` and `subject`.
    pub fn new(query: &'a [FastaRecord], subject: &'a [FastaRecord]) -> Self {
        Self {
            query,
            subject,
            query_name: "query",
            subject_name: "subject",
            database: false,
            pssm: None,
            query_end: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    // ```c++
    //     if (args.Exist(kArgDb) && args[kArgDb]) {
    // 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
    //
    //         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
    //                                              mol_type));
    // ```
    /// The query records against the BLAST database `name`: the command line's `-db name`, one
    /// or more space-separated names, resolved as there (the `BLASTDB` search path). The search
    /// opens the database.
    pub fn database(query: &'a [FastaRecord], name: &'a str) -> Self {
        Self {
            query,
            subject: &[],
            query_name: "query",
            subject_name: name,
            database: true,
            pssm: None,
            query_end: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
    // ```c++
    //     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
    //         CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
    //         _ASSERT(m_Pssm.Empty());
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:450
    // ```c++
    //         CRef<CPssmWithParameters> pssm = m_CmdLineArgs->GetInputPssm();
    // ```
    /// The query of the PSSM file `path` (psiblast's checkpoint file, as `-out_pssm` writes it)
    /// against the subject records: the command line's `-in_pssm path` of psiblast and of tblastn
    /// (PSI-TBLASTN), in place of `-query`. The search reads the file as the command line reads
    /// it, with its messages and errors (not `-`, the command line's standard input); the one
    /// query is the PSSM's sequence (query 0 of the results, `HspRecord::query_id` as the reports
    /// show it).
    pub fn pssm(path: &'a std::path::Path, subject: &'a [FastaRecord]) -> Self {
        Self {
            pssm: Some(path),
            ..Self::new(&[], subject)
        }
    }

    /// The query of the PSSM file `path` against the BLAST database `name`: the command line's
    /// `-in_pssm path -db name` (`Inputs::pssm`, `Inputs::database`).
    pub fn pssm_database(path: &'a std::path::Path, name: &'a str) -> Self {
        Self {
            pssm: Some(path),
            ..Self::database(&[], name)
        }
    }

    /// The same inputs, whose query records end with `end` ([`read_input_records`]): a search
    /// meets it where the command line's query batches meet it, after the reports of the batches
    /// before ([`QueryEnd`]). Without it the query records are those of a whole input. A PSSM
    /// query has no records and no end.
    pub fn query_end(self, end: &'a QueryEnd) -> Self {
        Self {
            query_end: Some(end),
            ..self
        }
    }

    /// How the reading of the query records ended (`Inputs::query_end`).
    pub(crate) fn end(&self) -> &'a QueryEnd {
        static INPUT: QueryEnd = QueryEnd::Input;
        self.query_end.unwrap_or(&INPUT)
    }

    /// The number of queries: the records, or the one sequence of a PSSM file.
    pub(crate) fn query_count(&self) -> usize {
        match self.pssm {
            Some(_) => 1,
            None => self.query.len(),
        }
    }

    /// The same inputs, shown with these names (the file names of a command line whose output
    /// the search is to reproduce). A database keeps its name, a PSSM query its file.
    pub fn named(self, query_name: &'a str, subject_name: &'a str) -> Self {
        Self {
            query_name,
            subject_name: if self.database {
                self.subject_name
            } else {
                subject_name
            },
            ..self
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
    const PROTEIN: &[u8] =
        b">p1 a protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n";

    #[test]
    fn records_carry_the_readers_local_ids_and_titles() {
        let query = read_fasta(Program::Blastp, Role::Query, PROTEIN).unwrap();
        let subject = read_fasta(Program::Blastp, Role::Subject, PROTEIN).unwrap();
        assert_eq!(query.len(), 1);
        assert_eq!(query[0].local_id, "Query_1");
        assert_eq!(subject[0].local_id, "Subject_1");
        assert_eq!(query[0].title(), b"p1 a protein");
        assert_eq!(query[0].seq().len(), 66);
        // An input of white space only has no record; an empty input neither.
        assert!(read_fasta(Program::Blastn, Role::Query, b" \n\t\n")
            .unwrap()
            .is_empty());
        assert!(read_fasta(Program::Blastn, Role::Query, b"")
            .unwrap()
            .is_empty());
    }

    #[test]
    fn the_readers_flags_follow_the_programs_molecules() {
        assert!(!is_protein(Program::Blastn, Role::Query));
        assert!(!is_protein(Program::Blastn, Role::Subject));
        assert!(is_protein(Program::Blastp, Role::Query));
        assert!(is_protein(Program::Blastp, Role::Subject));
        assert!(!is_protein(Program::Blastx, Role::Query));
        assert!(is_protein(Program::Blastx, Role::Subject));
        assert!(is_protein(Program::Tblastn, Role::Query));
        assert!(!is_protein(Program::Tblastn, Role::Subject));
        assert!(!is_protein(Program::Tblastx, Role::Query));
        assert!(!is_protein(Program::Tblastx, Role::Subject));
    }

    #[test]
    fn a_reader_error_is_the_command_lines_message() {
        // A first line that is neither a defline nor plausible data (`!`, `#` and `;` open
        // comment lines): NCBI's CFastaReader error, with the CLI's prefix.
        let error = read_fasta(Program::Blastn, Role::Query, b"xyz{}[]|~`^&*()=+\n");
        match error {
            Err(error @ Error::Input(_)) => {
                let message = error.to_string();
                assert!(
                    message.starts_with("BLAST query error: CFastaReader: Near line 1"),
                    "{message}"
                );
                assert!(!message.ends_with('\n'));
            }
            other => panic!("expected a reader error, got {other:?}"),
        }
    }

    #[test]
    fn a_file_reads_as_the_command_line_opens_it() {
        use std::io::Write;
        let directory = crate::utils::test_dir::temp_dir().join(format!(
            "losat-read-input-file-{}",
            crate::utils::test_dir::process_id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let plain = directory.join("q.fa");
        std::fs::write(&plain, PROTEIN).unwrap();
        assert_eq!(read_input_file(Role::Query, &plain).unwrap(), PROTEIN);
        let packed = directory.join("q.fa.gz");
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&packed).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(PROTEIN).unwrap();
        encoder.finish().unwrap();
        assert_eq!(read_input_file(Role::Subject, &packed).unwrap(), PROTEIN);
        let missing = directory.join("none.fa");
        let error = read_input_file(Role::Subject, &missing).unwrap_err();
        assert_eq!(
            error,
            Error::Input(
                format!(
                "Command line argument error: Argument \"subject\". File is not accessible:  `{}'",
                missing.display()
            )
                .into()
            )
        );
        // A name that is not UTF-8 keeps its bytes, as on the command line (line LIB-B).
        #[cfg(unix)]
        {
            use std::os::unix::ffi::OsStrExt;
            let missing = directory.join(std::ffi::OsStr::from_bytes(b"none\xE9.fa"));
            let error = read_input_file(Role::Query, &missing).unwrap_err();
            let expected = [
                b"Command line argument error: Argument \"query\". File is not accessible:  `"
                    .as_slice(),
                missing.as_os_str().as_bytes(),
                b"'",
            ]
            .concat();
            assert_eq!(error, Error::Input(expected));
        }
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn inputs_default_names_are_query_and_subject() {
        let records = read_fasta(Program::Blastp, Role::Query, PROTEIN).unwrap();
        let inputs = Inputs::new(&records, &records);
        assert_eq!(
            (inputs.query_name, inputs.subject_name),
            ("query", "subject")
        );
        let named = inputs.named("q.faa", "s.faa");
        assert_eq!((named.query_name, named.subject_name), ("q.faa", "s.faa"));
    }
}

#[cfg(test)]
#[path = "reader_end_tests.rs"]
mod reader_end_tests;
