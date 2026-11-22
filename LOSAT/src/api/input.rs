//! The inputs of a search: FASTA records in memory, read by the program's own reader.
//!
//! The command line reads its `-query` and `-subject` files with NCBI's FASTA reader
//! (`CBlastFastaInputSource`, ported in `blastinput::fasta_reader`); the library API reads
//! bytes with the same reader and the same flags, so a search of the records it returns is the
//! search of the files that hold those bytes (LOSAT Web ABI v2 `register`, docs/web/abi_v2.md
//! §4, reads them the same way).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:316-361
#![warn(missing_docs, unused)]

use super::error::Error;
use super::program::Program;
use super::search::Search;
use crate::blastinput::fasta_reader::{
    read_all, read_queries, FastaInputSource, FastaRecord, ParseErrorCode, QueryRecords, ReadError,
    ReaderConfig,
};
use crate::blastinput::input_files::{open_file_input, read_to_failure};
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
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:78-80
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

/// A query that a file holds in place of FASTA records (`Inputs::pssm`, `Inputs::msa`): the
/// search reads the file as the command line reads its `-in_pssm` or `-in_msa`, and the one query
/// is the sequence it gives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum QueryFile<'a> {
    /// A PSSM file (`-in_pssm` of psiblast and tblastn).
    Pssm(&'a std::path::Path),
    /// A multiple sequence alignment (`-in_msa` of psiblast), made into the query's PSSM.
    Alignment(&'a std::path::Path),
}

impl<'a> QueryFile<'a> {
    /// The file.
    pub(crate) fn path(self) -> &'a std::path::Path {
        match self {
            QueryFile::Pssm(path) | QueryFile::Alignment(path) => path,
        }
    }
}

/// The reader configuration of the command line for the input of `role`: the program's flags,
/// and NCBI's data loaders on (the command line without `.ncbirc` or `DATA_LOADERS`) or as the
/// application's settings give them (`data_loaders`, `ApplicationSettings::data_loaders`).
/// `writes_archive` is whether the search writes the BLAST archive (outfmt 11), or `None` where its
/// formats are not known yet (`read_inputs`).
fn reader_config(
    program: Program,
    role: Role,
    writes_archive: Option<bool>,
    data_loaders: bool,
) -> ReaderConfig {
    let name = program.reader_name();
    let protein = is_protein(program, role);
    let config = match role {
        Role::Query => ReaderConfig::query(name, protein, data_loaders),
        Role::Subject => ReaderConfig::subject(name, protein, data_loaders),
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
        FastaInputSource::from_bytes(bytes, reader_config(program, role, None, true)),
        role,
    )
}

/// Reads every record of an input source; a reader error is the command line's message: the
/// subjects end as `GetAllSeqs` ends (`ReadError::into_all_seqs_error`), the queries as a batch of
/// queries ends (`ReadError::into_app_error`, `GetNextSeqBatch` rethrows the reader's exception).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:205-214
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
/// The bytes of the file `path`, as the command line opens and reads its `-query` or `-subject`
/// argument: the same message for a file that does not open (`Error::Input`), a file named `.gz`,
/// `.bz2` or `.zst` is inflated (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, exception 8), and a read
/// that fails ends the bytes, as it ends the command line's input: a directory gives no byte (the
/// command line's query of a directory is empty). `read_fasta` and `read_inputs` read the bytes.
pub fn read_input_file(role: Role, path: &std::path::Path) -> Result<Vec<u8>, Error> {
    let name = match role {
        Role::Query => "query",
        Role::Subject => "subject",
    };
    let mut file = open_file_input(path, name)
        .map_err(|error| Error::Input(super::error::input_bytes(&error)))?;
    Ok(read_to_failure(&mut file))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:188-213
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
        true,
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
    /// `Empty CBlastQueryVector` (`LOSAT engine error: Empty CBlastQueryVector`).
    BlankLines,
    /// The reader's error on the record after the records read: the run ends with it in that
    /// record's batch.
    Error(QueryError),
}

/// The reader's error on a query record ([`QueryEnd::Error`]): the command line's message, and
/// the reader's messages about the record that it wrote before the error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryError {
    /// The command line's message without its last newline (`Error::Input`'s bytes).
    message: Vec<u8>,
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
    /// `ReadError::Unsupported`, LOSAT's refusal: its text (a `NativeError`'s bytes), and its
    /// exit status where it is a `NativeError`.
    Refused { exit: Option<i32>, text: Vec<u8> },
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
                        text: format!("{refusal:#}").into_bytes(),
                    },
                }
            }
            // The messages of the records are kept, not written: no write fails here.
            ReadError::Write(error) => ReadFailure::Refused {
                exit: None,
                text: error.to_string().into_bytes(),
            },
        };
        let message = super::error::input_bytes(&error.into_app_error());
        Self {
            message,
            warnings,
            failure,
        }
    }

    /// The command line's message, without its last newline: what [`read_inputs`] fails with
    /// (`Error::Input`).
    pub fn message(&self) -> &[u8] {
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
                ReadError::Unsupported(anyhow::anyhow!("{}", String::from_utf8_lossy(text)))
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
        f.write_str(&String::from_utf8_lossy(&self.message))
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
/// The `-strand` value of the arguments `command`, for the programs that have the option.
fn strand_of(command: &Commands) -> Option<&str> {
    match command {
        Commands::Blastn(args) => Some(args.strand.as_str()),
        Commands::Blastx(args) => Some(args.strand.as_str()),
        Commands::Tblastx(args) => Some(args.strand.as_str()),
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
/// archive (`check_archive_inputs`). `data_loaders` is whether NCBI's data loaders are on
/// (`ApplicationSettings::data_loaders` of the search's start; the library reads with them on,
/// `CallStart::settings`), under which the reader refuses a first line that NCBI would fetch as a
/// Seq-id and without which it reads that line as data. A subject's reader error is the command
/// line's message (`read_source`); the query keeps the records read before its reader's error and
/// the error (`QueryEnd`), which the search meets in its batch, as the command line's batches read
/// the queries (`read_query_records`). The one reading of `read_input_records` and LOSAT Web's
/// `register_inputs` (`web/adapter`, docs/web/abi_v2.md §4).
pub fn read_search_inputs(
    program: Program,
    command: &Commands,
    writes_archive: Option<bool>,
    data_loaders: bool,
    query: &[u8],
    subject: &[u8],
) -> Result<InputRecords, Error> {
    let parse_deflines = parse_deflines_of(command);
    let scope = parse_deflines.then(Default::default);
    let source = |role: Role, bytes| {
        let config = reader_config(program, role, writes_archive, data_loaders)
            .with_parse_deflines(parse_deflines);
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
/// The inputs of one search: the query and subject records, and the names that the reports
/// show for them where the command line shows its `-query` and `-subject` file names (the
/// `Database:` line of a `-subject` search, `User specified sequence set (Input: NAME)`).
///
/// The records come from [`read_fasta`] with the search's program and the input's role, or
/// from [`FastaRecord::new`]. A database subject comes from [`Inputs::database`]: the engine
/// opens it by its name, as the command line opens `-db`. A PSSM query of psiblast and tblastn
/// comes from [`Inputs::pssm`] and [`Inputs::pssm_database`]: the engine reads the file, as the
/// command line reads `-in_pssm`. A query alignment of psiblast comes from [`Inputs::msa`] and
/// [`Inputs::msa_database`]: the engine reads the file as the command line reads `-in_msa`.
#[derive(Debug, Clone, Copy)]
pub struct Inputs<'a> {
    /// The query records (`-query`); none for a query file (`query_file`).
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
    /// The file whose sequence is the one query (`-in_pssm`, `-in_msa`), in place of query
    /// records.
    pub(crate) query_file: Option<QueryFile<'a>>,
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
            query_file: None,
            query_end: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
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
            query_file: None,
            query_end: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:450
    /// The query of the PSSM file `path` (psiblast's checkpoint file, as `-out_pssm` writes it)
    /// against the subject records: the command line's `-in_pssm path` of psiblast and of tblastn
    /// (PSI-TBLASTN), in place of `-query`. The search reads the file as the command line reads
    /// it, with its messages and errors (not `-`, the command line's standard input); the one
    /// query is the PSSM's sequence (query 0 of the results, `HspRecord::query_id` as the reports
    /// show it).
    pub fn pssm(path: &'a std::path::Path, subject: &'a [FastaRecord]) -> Self {
        Self {
            query_file: Some(QueryFile::Pssm(path)),
            ..Self::new(&[], subject)
        }
    }

    /// The query of the PSSM file `path` against the BLAST database `name`: the command line's
    /// `-in_pssm path -db name` (`Inputs::pssm`, `Inputs::database`).
    pub fn pssm_database(path: &'a std::path::Path, name: &'a str) -> Self {
        Self {
            query_file: Some(QueryFile::Pssm(path)),
            ..Self::database(&[], name)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1184-1229
    /// The query of the alignment file `path` against the subject records: the command line's
    /// `-in_msa path` of psiblast, in place of `-query`. The search reads the file as the command
    /// line reads it, with its messages and errors (not `-`, the command line's standard input);
    /// the one query is the alignment's master (query 0 of the results, `HspRecord::query_id` as
    /// the reports show it), whose row `-msa_master_idx` and `-ignore_msa_master` of the argv
    /// choose, as on the command line. Only psiblast searches an alignment.
    pub fn msa(path: &'a std::path::Path, subject: &'a [FastaRecord]) -> Self {
        Self {
            query_file: Some(QueryFile::Alignment(path)),
            ..Self::new(&[], subject)
        }
    }

    /// The query of the alignment file `path` against the BLAST database `name`: the command
    /// line's `-in_msa path -db name` (`Inputs::msa`, `Inputs::database`).
    pub fn msa_database(path: &'a std::path::Path, name: &'a str) -> Self {
        Self {
            query_file: Some(QueryFile::Alignment(path)),
            ..Self::database(&[], name)
        }
    }

    /// The same inputs, whose query records end with `end` ([`read_input_records`]): a search
    /// meets it where the command line's query batches meet it, after the reports of the batches
    /// before ([`QueryEnd`]). Without it the query records are those of a whole input. A PSSM or
    /// alignment query has no records and no end.
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

    /// The number of queries: the records, or the one sequence of a PSSM or alignment file.
    pub(crate) fn query_count(&self) -> usize {
        match self.query_file {
            Some(_) => 1,
            None => self.query.len(),
        }
    }

    /// The same inputs, shown with these names (the file names of a command line whose output
    /// the search is to reproduce). A database keeps its name, a PSSM or alignment query its file.
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
#[path = "input_tests.rs"]
mod tests;
