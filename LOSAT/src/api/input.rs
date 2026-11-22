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

use std::io::Read;

use super::error::Error;
use super::program::Program;
use super::search::Search;
use crate::blastinput::fasta_reader::{read_all, FastaInputSource, FastaRecord, ReaderConfig};
use crate::blastinput::input_files::open_file_input;
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

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:188-213
/// Reads the query and subject FASTA bytes of `search`, as the command line reads its `-query`
/// and `-subject` files: the records that `read_fasta` gives, and with `-parse_deflines` in the
/// search's argv the ids of their deflines. The subjects come first, in one scope of ids with the
/// queries, as NCBI reads them: a query and a subject with one id are one sequence to the reports.
///
/// A search with `-parse_deflines` takes its records from here; a search of a database gives no
/// subject bytes. Returns `(query, subject)`.
pub fn read_inputs(
    search: &Search,
    query: &[u8],
    subject: &[u8],
) -> Result<(Vec<FastaRecord>, Vec<FastaRecord>), Error> {
    read_search_inputs(
        search.program(),
        &search.any_subject(),
        None,
        true,
        query,
        subject,
    )
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
/// Seq-id and without which it reads that line as data. A reader error is the command line's
/// message (`read_source`). The one reading of `read_inputs` and LOSAT Web's `register_inputs`
/// (`web/adapter`, docs/web/abi_v2.md §4). Returns `(query, subject)`.
pub fn read_search_inputs(
    program: Program,
    command: &Commands,
    writes_archive: Option<bool>,
    data_loaders: bool,
    query: &[u8],
    subject: &[u8],
) -> Result<(Vec<FastaRecord>, Vec<FastaRecord>), Error> {
    let parse_deflines = parse_deflines_of(command);
    let scope = parse_deflines.then(Default::default);
    let read = |role: Role, bytes: &[u8]| {
        let config = reader_config(program, role, writes_archive, data_loaders)
            .with_parse_deflines(parse_deflines);
        let config = match strand_of(command) {
            Some(strand) if role == Role::Query => config.with_strand(strand),
            _ => config,
        };
        let mut source = FastaInputSource::from_bytes(bytes, config);
        source.set_id_scope(scope.clone());
        read_source(source, role)
    };
    let subject = read(Role::Subject, subject)?;
    let query = read(Role::Query, query)?;
    Ok((query, subject))
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
    /// The search strategy whose queries and subject (or database) these are
    /// (`Inputs::strategy`).
    pub(crate) strategy: Option<&'a super::strategy::Strategy>,
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
            strategy: None,
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
            strategy: None,
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

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:565-599
    /// The queries and the subjects (or the database) that the search strategy `strategy` carries
    /// (`Search::read_strategy`): the command line's `-import_search_strategy FILE` without
    /// `-query`, `-subject` or `-db`. The search is the one that read the strategy; its argv
    /// supplies the options that the import takes from the command line, and the strategy the
    /// others. The reports name the subject set as the command line names an imported one
    /// (`User specified sequence set.`).
    pub fn strategy(strategy: &'a super::strategy::Strategy) -> Self {
        Self {
            query: strategy.queries(),
            subject: strategy.subjects(),
            query_name: "query",
            subject_name: strategy.database().unwrap_or("subject"),
            database: strategy.database().is_some(),
            query_file: None,
            strategy: Some(strategy),
        }
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
