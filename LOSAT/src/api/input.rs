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
use crate::blastinput::fasta_reader::{read_all, FastaInputSource, FastaRecord, ReaderConfig};
use crate::blastinput::input_files::open_file_input;

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
        (Program::Blastp, _) => true,
        (Program::Tblastn, Role::Query) => true,
        (Program::Tblastn, Role::Subject) => false,
        (Program::Blastx, Role::Query) => false,
        (Program::Blastx, Role::Subject) => true,
    }
}

/// The reader configuration of the command line for the input of `role`: the program's flags,
/// the data loaders on (the command line without `.ncbirc` or `DATA_LOADERS`).
fn reader_config(program: Program, role: Role) -> ReaderConfig {
    let name = program.reader_name();
    let protein = is_protein(program, role);
    match role {
        Role::Query => ReaderConfig::query(name, protein, true),
        Role::Subject => ReaderConfig::subject(name, protein, true),
    }
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
    read_source(FastaInputSource::from_bytes(
        bytes,
        reader_config(program, role),
    ))
}

/// Reads every record of an input source; a reader error is the command line's message.
fn read_source<R: std::io::Read>(
    mut source: FastaInputSource<R>,
) -> Result<Vec<FastaRecord>, Error> {
    read_all(&mut source, &mut |_| Ok(())).map_err(|error| {
        Error::Input(
            format!("{:#}", error.into_app_error())
                .trim_end_matches('\n')
                .to_string(),
        )
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
    let input_error = |error: anyhow::Error| {
        Error::Input(format!("{error:#}").trim_end_matches('\n').to_string())
    };
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
/// A search with `-parse_deflines` takes its records from here. Returns `(query, subject)`.
pub fn read_inputs(
    search: &Search,
    query: &[u8],
    subject: &[u8],
) -> Result<(Vec<FastaRecord>, Vec<FastaRecord>), Error> {
    let parse_deflines = search.parse_deflines();
    let scope = parse_deflines.then(Default::default);
    let read = |role: Role, bytes: &[u8]| {
        let config = reader_config(search.program(), role).with_parse_deflines(parse_deflines);
        let mut source = FastaInputSource::from_bytes(bytes, config);
        source.set_id_scope(scope.clone());
        read_source(source)
    };
    let subject = read(Role::Subject, subject)?;
    let query = read(Role::Query, query)?;
    Ok((query, subject))
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
/// from [`FastaRecord::new`]. The library API has no database (`-db`) input yet
/// (docs/losat_library_surface_plan.md §3.3): it will be another constructor, so the fields are
/// not part of the API.
#[derive(Debug, Clone, Copy)]
pub struct Inputs<'a> {
    /// The query records (`-query`).
    pub(crate) query: &'a [FastaRecord],
    /// The subject records (`-subject`).
    pub(crate) subject: &'a [FastaRecord],
    /// What the reports show as the query's file name. The command line's default is the
    /// `-query` argument; the API's is `query`.
    pub(crate) query_name: &'a str,
    /// What the reports show as the subject set's file name. The command line's default is the
    /// `-subject` argument; the API's is `subject`.
    pub(crate) subject_name: &'a str,
}

impl<'a> Inputs<'a> {
    /// The query and subject records, shown as `query` and `subject`.
    pub fn new(query: &'a [FastaRecord], subject: &'a [FastaRecord]) -> Self {
        Self {
            query,
            subject,
            query_name: "query",
            subject_name: "subject",
        }
    }

    /// The same inputs, shown with these names (the file names of a command line whose output
    /// the search is to reproduce).
    pub fn named(self, query_name: &'a str, subject_name: &'a str) -> Self {
        Self {
            query_name,
            subject_name,
            ..self
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            Err(Error::Input(message)) => {
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
        let directory =
            std::env::temp_dir().join(format!("losat-read-input-file-{}", std::process::id()));
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
            Error::Input(format!(
                "Command line argument error: Argument \"subject\". File is not accessible:  `{}'",
                missing.display()
            ))
        );
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
