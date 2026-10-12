//! `LOSAT makedb`: the application layer of makeblastdb
//! (`c++/src/app/blastdb/makeblastdb.cpp`, `CMakeBlastDBApp`) for a FASTA
//! input and a version 4 database without parsed Seq-ids.

use std::io::{Read, Write};

use clap::Args;

use super::build_db::BuildDatabase;
use super::data_size::{string_to_uint8_data_size, uint8_to_string_data_size};
use super::format_guess::appears_fasta;
use super::time::LocalTime;
use super::WriteDbError;
use crate::seqdb::{SeqDb, SeqDbError, SeqType};

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:224-244
// ```c++
//     arg_desc->SetCurrentGroup("Input options");
//     arg_desc->AddDefaultKey(kInput, "input_file",
//                             "Input file/database name",
//                             CArgDescriptions::eInputFile, "-");
//     arg_desc->AddDefaultKey("input_type", "type",
//                              "Type of the data specified in input_file",
//                              CArgDescriptions::eString, "fasta");
//     arg_desc->SetConstraint("input_type", &(*new CArgAllow_Strings,
//                                             "fasta", "blastdb",
//                                             "asn1_bin",
//                                             "asn1_txt"));
//
//     arg_desc->AddKey(kArgDbType, "molecule_type",
//                      "Molecule type of target db", CArgDescriptions::eString);
//     arg_desc->SetConstraint(kArgDbType, &(*new CArgAllow_Strings,
//                                         "nucl", "prot"));
//
//     arg_desc->SetCurrentGroup("Configuration options");
//     arg_desc->AddOptionalKey(kArgDbTitle, "database_title",
//                              "Title for BLAST database\n" + dflt,
//                              CArgDescriptions::eString);
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:292-306
// ```c++
//     arg_desc->SetCurrentGroup("Output options");
//     arg_desc->AddOptionalKey(kOutput, "database_name",
//                              "Name of BLAST database to be created\n" + dflt +
//                              "Required if multiple file(s)/database(s) are "
//                              "provided as input",
//                              CArgDescriptions::eString);
//     arg_desc->AddDefaultKey("blastdb_version", "version",
//                              "Version of BLAST database to be created",
//                              CArgDescriptions::eInteger,
//                              NStr::NumericToString(static_cast<int>(eBDB_Version5)));
//     arg_desc->SetConstraint("blastdb_version",
//                             new CArgAllow_Integers(eBDB_Version4, eBDB_Version5));
//     arg_desc->AddDefaultKey("max_file_sz", "number_of_bytes",
//                             "Maximum file size for BLAST database files",
//                             CArgDescriptions::eString, "3GB");
// ```
// LOSAT writes format version 4 only (PD-LOSAT-DATABASE-SEARCH, P2), so 4 is
// the default of `-blastdb_version` here and 5 is refused. The options of
// makeblastdb that LOSAT does not support are not registered; `crate::cli`
// names them in its error.
/// Arguments of `LOSAT makedb`, under makeblastdb's names.
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct MakedbArgs {
    /// Input FASTA file (`-` is standard input)
    #[arg(long = "in", default_value = "-", value_name = "input_file")]
    pub input: String,
    /// Type of the data specified in input_file (only fasta is supported)
    #[arg(long, default_value = "fasta", value_name = "type",
          value_parser = ["fasta", "blastdb", "asn1_bin", "asn1_txt"])]
    pub input_type: String,
    /// Molecule type of target db
    #[arg(long, value_name = "molecule_type", value_parser = ["nucl", "prot"])]
    pub dbtype: String,
    /// Title for BLAST database [default: the input file name given to -in]
    #[arg(long, value_name = "database_title")]
    pub title: Option<String>,
    /// Name of BLAST database to be created [default: the input file name given to -in]
    #[arg(long, value_name = "database_name")]
    pub out: Option<String>,
    /// Version of BLAST database to be created (only 4 is supported)
    #[arg(long, default_value_t = 4, value_name = "version",
          value_parser = clap::value_parser!(i32).range(4..=5))]
    pub blastdb_version: i32,
    /// Maximum file size for BLAST database files
    #[arg(long, default_value = "3GB", value_name = "number_of_bytes")]
    pub max_file_sz: String,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:76-77
// ```c++
// static const char* s_AutoHelp     = "h";
// static const char* s_AutoHelpFull = "help";
// ```
/// The options of makeblastdb 2.17.0 that `LOSAT makedb` refuses
/// (makeblastdb.cpp:246-333), and `-h`, the short usage that the NCBI C++
/// Toolkit gives every program (`-help` prints LOSAT's own usage, as for the
/// search programs).
pub const UNSUPPORTED_OPTIONS: [&str; 15] = [
    "h",
    "parse_seqids",
    "hash_index",
    "randomize",
    "mask_data",
    "mask_id",
    "mask_desc",
    "gi_mask",
    "gi_mask_name",
    "metadata_output_prefix",
    "logfile",
    "taxid",
    "taxid_map",
    "oid_masks",
    "version",
];

/// How `LOSAT makedb` ends when it does not build a database, or cannot
/// open the one it built.
#[derive(Debug)]
pub enum MakedbError {
    /// An error of makeblastdb, with makeblastdb's exit code.
    Build(WriteDbError),
    /// The database was written and cannot be opened: NCBI's
    /// `CSeqDBException`, the text makeblastdb prints after "BLAST Database
    /// error: ".
    Database(SeqDbError),
    /// A request makeblastdb serves and LOSAT does not.
    Unsupported(String),
}

impl MakedbError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:149-163
    // ```c++
    // #define BLAST_EXIT_SUCCESS          0
    // /// Command line binary exit code: error in input query/options
    // #define BLAST_INPUT_ERROR           1
    // /// Command line binary exit code: error in database/subject
    // #define BLAST_DATABASE_ERROR        2
    // /// Command line binary exit code: error in BLAST engine
    // #define BLAST_ENGINE_ERROR          3
    // /// Command line binary exit code: BLAST run out of memory
    // #define BLAST_OUT_OF_MEMORY         4
    // /// Command line binary exit code: Network error encountered
    // #define BLAST_NETWORK_ERROR         5
    // /// Command line binary exit code: error when writing output
    // #define BLAST_OUTPUT_ERROR          6
    // /// Command line binary exit code: unknown error
    // #define BLAST_UNKNOWN_ERROR         255
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
    // ```c++
    // #define CATCH_ALL(exit_code)                                                \
    //     catch (const CInvalidDataException& e) {                                \
    //         LOG_POST(Error << "BLAST options error: " << e.GetMsg());           \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:199-203
    // ```c++
    //         else {                                                              \
    //         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
    //         	exit_code = BLAST_DATABASE_ERROR;                               \
    //         }                                                                   \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:212-215
    // ```c++
    //     catch (const CWriteDBException& e) {                                    \
    //         LOG_POST(Error << "BLAST Database creation error: " << e.GetMsg()); \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:256-259
    // ```c++
    //     catch (const std::bad_alloc&) {                                         \
    //         LOG_POST(Error << "BLAST ran out of memory");                       \
    //         exit_code = BLAST_OUT_OF_MEMORY;                                    \
    //     }                                                                       \
    // ```
    /// The exit code: 1 (`BLAST_INPUT_ERROR`) for an error of the options or
    /// of the database creation, 2 (`BLAST_DATABASE_ERROR`) when the new
    /// database cannot be opened, 4 (`BLAST_OUT_OF_MEMORY`) for a failed
    /// allocation, and 255 (`BLAST_UNKNOWN_ERROR`) otherwise, which includes
    /// the failures of the file system that LOSAT reports itself.
    pub fn exit_code(&self) -> i32 {
        match self {
            MakedbError::Build(
                WriteDbError::Other(_) | WriteDbError::Io { .. } | WriteDbError::Rename { .. },
            ) => 255,
            MakedbError::Build(WriteDbError::OutOfMemory) => 4,
            MakedbError::Database(_) => 2,
            MakedbError::Build(WriteDbError::Creation(_) | WriteDbError::Options(_))
            | MakedbError::Unsupported(_) => 1,
        }
    }
}

impl std::fmt::Display for MakedbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MakedbError::Build(error) => error.fmt(f),
            MakedbError::Database(error) => write!(f, "BLAST Database error: {error}"),
            MakedbError::Unsupported(message) => write!(f, "LOSAT makedb: {message}"),
        }
    }
}

impl From<WriteDbError> for MakedbError {
    fn from(error: WriteDbError) -> Self {
        MakedbError::Build(error)
    }
}

fn options_error(message: impl Into<String>) -> MakedbError {
    MakedbError::Build(WriteDbError::Options(message.into()))
}

fn unsupported(message: impl Into<String>) -> MakedbError {
    MakedbError::Unsupported(message.into())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1762-1807
// ```c++
// void SeqDB_SplitQuoted(const string             & dbname,
//                        vector<CSeqDB_Substring> & dbs,
//                        bool keep_quote)
// {
//     // split names
//
//     const char * sp = dbname.data();
//
//     bool quoted = false;
//     unsigned begin = 0;
//
//     for(unsigned i = 0; i < dbname.size(); i++) {
//         char ch = dbname[i];
//
//         if (quoted) {
//             // Quoted mode sees '"' as the only actionable token.
//             if (ch == '"') {
//                 if (begin < i) {
//                 	if(keep_quote) i++;
//                     dbs.push_back(CSeqDB_Substring(sp + begin, sp + i));
//                 }
//                 begin = i + 1;
//                 quoted = false;
//             }
//         } else {
//             // Non-quote mode: Space or quote starts the next string.
//
//             if (ch == ' ') {
//                 if (begin < i) {
//                     dbs.push_back(CSeqDB_Substring(sp + begin, sp + i));
//                 }
//                 begin = i + 1;
//             } else if (ch == '"') {
//                 if (begin < i) {
//                     dbs.push_back(CSeqDB_Substring(sp + begin, sp + i));
//                 }
//                 begin = keep_quote ? i : i + 1;
//                 quoted = true;
//             }
//         }
//     }
//
//     if (begin < dbname.size()) {
//         dbs.push_back(CSeqDB_Substring(sp + begin, sp + dbname.size()));
//     }
// }
// ```
// makeblastdb splits `-in` without keeping the quotes (`keep_quote` false).
/// `SeqDB_SplitQuoted`: the input names of `-in` (spaces separate names,
/// double quotes keep a name with spaces together).
fn split_quoted(names: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut quoted = false;
    let mut begin = 0;
    for (i, ch) in names.bytes().enumerate() {
        if quoted {
            if ch == b'"' {
                if begin < i {
                    out.push(&names[begin..i]);
                }
                begin = i + 1;
                quoted = false;
            }
        } else if ch == b' ' || ch == b'"' {
            if begin < i {
                out.push(&names[begin..i]);
            }
            begin = i + 1;
            quoted = ch == b'"';
        }
    }
    if begin < names.len() {
        out.push(&names[begin..]);
    }
    out
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1128-1136
// ```c++
// void CMakeBlastDBApp::x_BuildDatabase()
// {
//     const string env_skip =
//             GetEnvironment().Get("NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS");
//     m_SkipUnver = (env_skip.empty() == false);
//
//     const string dont_scan_bioseq =
//         GetEnvironment().Get("NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS");
//     const bool scan_bioseq_4_cfastareader_usrobj = static_cast<bool>(dont_scan_bioseq.empty());
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:562-566
// ```c++
// 	try {
//         while(!data.eof())
//         {
//             CSeqEntrySource src(data, fmt, m_SkipUnver);
//             found = found || m_DB->AddSequences(src);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1733-1738,1764-1766
// ```c++
//     string fasta;
//
//     // Scan the CBioseq for the CFastaReader user object.
//
//     if (scan_bioseq_4_cfastareader_usrobj) {
//         ITERATE(list< CRef< CSeqdesc > >, iter, bioseq.GetDescr().Get()) {
// ...
//
//     }
//     if (fasta.empty())
// ```
// makeblastdb reads two variables of its own. LOSAT reproduces makeblastdb
// with both unset (or empty, which `Get` cannot tell from unset).
//
// * `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS`: with a
//   value, the definition line that the FASTA reader keeps in its user
//   object is not used, and the title stored for a sequence comes from the
//   title descriptor instead. The titles differ for a definition line with a
//   tab, a Ctrl-A or white space at its ends, so the variable changes the
//   database of a FASTA input; LOSAT refuses to run with it.
// * `NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS`: `m_SkipUnver` is read by
//   `CSeqEntrySource` only, the source of ASN.1 input. A FASTA input goes
//   through `x_AddFasta`, which does not look at it, so the variable changes
//   nothing here and is not refused.
/// Refuses the environment variable of makeblastdb that changes the database
/// of a FASTA input. `env` gives the value of a variable.
pub fn check_environment(env: &dyn Fn(&str) -> Option<std::ffi::OsString>) -> Result<(), String> {
    const DONT_SCAN: &str = "NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS";
    if env(DONT_SCAN).is_some_and(|value| !value.is_empty()) {
        return Err(format!(
            "the environment variable {DONT_SCAN}, with which NCBI makeblastdb stores other \
             titles for some FASTA definition lines, is not supported by LOSAT's makedb"
        ));
    }
    Ok(())
}

/// The checks of the arguments that LOSAT makes before it touches a file:
/// what LOSAT does not support, and the value of `-max_file_sz`.
fn check_arguments(args: &MakedbArgs) -> Result<u64, MakedbError> {
    if args.input_type != "fasta" {
        return Err(unsupported(format!(
            "-input_type {} is not supported; only fasta input is",
            args.input_type
        )));
    }
    if args.blastdb_version != 4 {
        return Err(unsupported(
            "-blastdb_version 5 is not supported; LOSAT writes format version 4 \
             (-blastdb_version 4, the default here)",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1267-1272
    // ```c++
    //     Uint8 bytes = NStr::StringToUInt8_DataSize(args["max_file_sz"].AsString());
    //     static const Uint8 MAX_VOL_FILE_SIZE = 0x100000000;
    //     if (bytes >= MAX_VOL_FILE_SIZE) {
    //         NCBI_THROW(CInvalidDataException, eInvalidInput,
    //                 "max_file_sz must be < 4 GiB");
    //     }
    // ```
    // makeblastdb reads the value after it has written the head of the log
    // and removed an existing database of the same name; LOSAT reads it
    // first, so that a wrong value costs nothing. A value NCBI cannot convert
    // throws a `CStringException`, which `CATCH_ALL` reports through its
    // `CException` clause (`MakedbError::exit_code`): "Error: NCBI C++
    // Exception:" and a line with the source path of NCBI's build, the
    // exception's class and its sentence, exit status 255. LOSAT prints the
    // sentence after "Error: -max_file_sz: ", with the same status.
    let bytes = string_to_uint8_data_size(&args.max_file_sz).map_err(|error| {
        MakedbError::Build(WriteDbError::Other(format!(
            "-max_file_sz: {}",
            error.message
        )))
    })?;
    if bytes >= 0x1_0000_0000 {
        return Err(options_error("max_file_sz must be < 4 GiB"));
    }
    Ok(bytes)
}

/// The messages of the FASTA reader go to the standard error stream as they
/// arrive.
fn warn_to(err: &mut dyn Write) -> impl FnMut(&[u8]) + '_ {
    move |line: &[u8]| {
        let _ = err.write_all(line);
        let _ = err.write_all(b"\n");
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1146-1197
// ```c++
//     string dbname = (args[kOutput].HasValue()
//                      ? args[kOutput]
//                      : args[kInput]).AsString();
//     ...
//     vector<string> input_files;
//     NStr::Split(dbname, kInputSeparators, input_files);
//     if (dbname == "-" || input_files.size() > 1) {
//         NCBI_THROW(CInvalidDataException, eInvalidInput,
//             "Please provide a database name using -" + kOutput);
//     }
//     ...
//     string title = (args[kArgDbTitle].HasValue()
//                     ? args[kArgDbTitle]
//                     : args[kInput]).AsString();
//     ...
//     if (title == "-") {
//         NCBI_THROW(CInvalidDataException, eInvalidInput,
//                          "Please provide a title using -title");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1273-1276
// ```c++
//     *m_LogFile << "Maximum file size: "
//                << Uint8ToString_DataSize(bytes) << endl;
//
//     m_DB->SetMaxFileSize(bytes);
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1294-1296
// ```c++
//     x_ProcessInputData(args[kInput].AsString(), is_protein);
//
//     bool success = m_DB->EndBuild();
// ```
// After a successful build makeblastdb opens the new database with its
// reader and writes the JSON metadata file (`.pjs` / `.njs`); an error of
// that open ends makeblastdb with "BLAST Database error" and exit status 2
// although the database was written. LOSAT opens the database for the same
// outcome (`open_new_database`); it does not write the JSON file, which
// NCBI's programs do not read when they search.
/// `CMakeBlastDBApp::x_BuildDatabase`. The log goes to `out` (makeblastdb's
/// standard output), warnings and diagnostics to `err`; `stdin` is read when
/// `-in` is `-`. `now` is the local time for the database's date (`None`
/// where the platform has no clock).
pub fn build_database(
    args: &MakedbArgs,
    now: Option<&LocalTime>,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let max_file_size = check_arguments(args)?;
    let is_protein = args.dbtype == "prot";

    let dbname = args.out.as_deref().unwrap_or(&args.input);
    // NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-275,282-301
    // ```c++
    //     void Do(TContainer&     target,
    //             TPosContainer&  token_pos,
    //             const TString&  empty_str = TString())
    //     {
    //         auto target_initial_size = target.size();
    //
    //         // Special cases
    //         if (m_Str.empty()) {
    //             return;
    //         } else if (m_Delim.empty()) {
    //             target.push_back(m_Str);
    //             token_pos.push_back(0);
    // ...
    //
    //         CTempStringList part_collector(m_Storage);
    //         SIZE_TYPE prev_pos;
    //         SIZE_TYPE delim_pos = NPOS;
    //         m_Pos = 0;
    //         do {
    //             Advance(&part_collector, &prev_pos, &delim_pos);
    //             target.push_back(empty_str); // reserve space for value added next line
    //             part_collector.Join(&target.back());
    //             part_collector.Clear();
    //             token_pos.push_back(prev_pos);
    //         } while ( !AtEnd() );
    //
    //         if ( (m_Flags & NStr::fSplit_Truncate_End) == 0 ) {
    //             // account training delimiter
    //             if (delim_pos != NPOS) {
    //                 // add empty token after last delimiter
    //                 target.push_back(empty_str);
    //                 token_pos.push_back(delim_pos + 1);
    //             }
    // ```
    // `NStr::Split(dbname, " ", input_files)` runs with no flag
    // (ncbistr.hpp:2551-2555): adjacent delimiters are not merged and none is
    // truncated, so every space ends a token, an empty one included ("db "
    // is "db" and "", " db" is "" and "db"). A name with a space anywhere
    // gives more than one token; the empty name gives none.
    if dbname == "-" || dbname.contains(' ') {
        return Err(options_error("Please provide a database name using -out"));
    }
    let title = args.title.as_deref().unwrap_or(&args.input);
    if title == "-" {
        return Err(options_error("Please provide a title using -title"));
    }

    // One FASTA input: makeblastdb concatenates several (`-in "a b"`).
    let names = split_quoted(&args.input);
    let [name] = names.as_slice() else {
        return Err(unsupported(
            "-in names several input files (or none); LOSAT makedb reads one FASTA file",
        ));
    };
    let name = *name;
    let lower = name.to_ascii_lowercase();
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1041-1050
    // ```c++
    // 		if( (is_gz  = NStr::EndsWith( seq_file, ".gz", NStr::eNocase)) ||
    // 		    (is_bz2 = NStr::EndsWith( seq_file, ".bz2",NStr::eNocase)) ||
    // 		    (is_zst = NStr::EndsWith( seq_file, ".zst",NStr::eNocase))
    // 		  )
    // 		{
    // 		    if( is_bz2  ) method = CDecompressIStream::eBZip2;
    // 		    if( is_zst  ) method = CDecompressIStream::eZstd;
    // 		    comp_seq_istream.reset(new CDecompressIStream(f , method ));
    // 		    data_in =  comp_seq_istream.get();
    // 		}
    // ```
    if name != "-"
        && [".gz", ".bz2", ".zst"]
            .iter()
            .any(|ext| lower.ends_with(ext))
    {
        return Err(unsupported(format!(
            "{name}: makeblastdb reads a file of this name as compressed; LOSAT makedb reads \
             uncompressed FASTA only"
        )));
    }

    let Some(now) = now else {
        return Err(unsupported(
            "this platform has no clock for the database's date",
        ));
    };
    let mut db = BuildDatabase::new(dbname, title, is_protein, out, now)?;
    let _ = writeln!(
        db.log(),
        "Maximum file size: {}",
        uint8_to_string_data_size(max_file_size)
    );
    db.set_max_file_size(max_file_size);

    let result = if name == "-" {
        add_fasta(&mut db, stdin, err)
    } else {
        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:804-812
        // ```c++
        //         CFile input_file(seq_file);
        //         if ( !input_file.Exists() ) {
        //             string error_msg = "File " + seq_file + " does not exist";
        //             NCBI_THROW(CInvalidDataException, eInvalidInput, error_msg);
        //         }
        //         if (input_file.GetLength() == 0) {
        //             string error_msg = "File " + seq_file + " is empty";
        //             NCBI_THROW(CInvalidDataException, eInvalidInput, error_msg);
        //         }
        // ```
        match std::fs::metadata(name) {
            Ok(meta) if meta.is_file() && meta.len() == 0 => {
                Err(options_error(format!("File {name} is empty")))
            }
            // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:832-836
            // ```c++
            //         if(input_type == eFasta && x_ConvertToSupportedType(x_GuessFileType( *data_in ) ) != eFasta) {
            //             string msg = "\nInput file " + seq_file + " does NOT appear to be FASTA (processing anyway).\n" \
            //                  + "Advise validating database with 'blastdbcheck -dbtype [prot|nucl] -db ${DBNAME}'\n";
            //             ERR_POST(Warning << msg);
            //         }
            // ```
            // makeblastdb opens the file a second time for this test, so the
            // build reads it from its start. `ERR_POST` ends the message
            // with a line feed.
            Ok(meta) if meta.is_file() => match std::fs::File::open(name) {
                Ok(mut file) => {
                    let looks_like_fasta = std::fs::File::open(name).is_ok_and(appears_fasta);
                    if !looks_like_fasta {
                        let _ = write!(
                            err,
                            "Warning: [makeblastdb] \nInput file {name} does NOT appear to be \
                             FASTA (processing anyway).\nAdvise validating database with \
                             'blastdbcheck -dbtype [prot|nucl] -db ${{DBNAME}}'\n\n"
                        );
                    }
                    add_fasta(&mut db, &mut file, err)
                }
                // A file that exists and cannot be opened. makeblastdb does
                // not test the stream it opened: `AddFasta` skips a stream
                // that is not good, no sequence is added, the build ends
                // with "No volumes were created." and exit status 0 (after
                // the FASTA warning above, for an empty sample). LOSAT does
                // not reproduce that: it reports the file, with the exit
                // status of its file-system errors.
                // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1398-1408
                // ```c++
                // bool CBuildDatabase::AddFasta(CNcbiIstream & fasta_file)
                // {
                //     // Add any fasta sequences as well.
                //     bool success = true;
                //
                //     if (fasta_file) {
                //         CFastaBioseqSource fbs(fasta_file,
                //                                m_IsProtein,
                //                                m_ParseIDs,
                //                                m_LongIDs);
                //
                // ```
                Err(source) => Err(MakedbError::Build(WriteDbError::Io {
                    path: name.into(),
                    source,
                })),
            },
            _ => Err(options_error(format!("File {name} does not exist"))),
        }
    };
    result?;

    let closed = {
        let mut warn = warn_to(err);
        db.end_build(false, &mut warn)
    };
    if closed? {
        open_new_database(db.output_dbname(), is_protein)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1296-1300,1306-1310
// ```c++
//     bool success = m_DB->EndBuild();
//     if(success) {
//     	string new_db = m_DB->GetOutputDbName();
//     	CSeqDB::ESeqType t = is_protein? CSeqDB::eProtein: CSeqDB::eNucleotide;
//     	CSeqDB sdb(new_db, t);
// ...
//     	CRef<CBlast_db_metadata> m = sdb.GetDBMetaData(output_prefix);
//     	string extn (kEmptyStr);
//     	SeqDB_GetMetadataFileExtension(is_protein, extn);
//     	string metadata_filename = new_db + "." + extn;
//         ofstream out(metadata_filename.c_str());
// ```
// `CSeqDB sdb(new_db, t)` throws a `CSeqDBException` when the reader cannot
// open what was written; the usual cause is an alias file of an earlier
// database of the same name that `DeleteBlastDb` could not resolve and
// therefore left in place, which now hides the new single volume. LOSAT
// resolves the database here for that outcome alone, and it only resolves:
// the alias files are read and every volume is located, but no volume is
// opened. The volumes were written a moment ago, and LOSAT's reader holds the
// index of an open volume in memory (12 bytes for each sequence). So the
// error `CSeqDB` gives for more volumes than the process may hold open files
// ("Too many open files, please raise the open file limit", exit status 2
// after a complete build) is not reproduced.
//
// The path is one database name (`SeqDb::locate_one_name`), which departs
// from NCBI on purpose: `CSeqDB` splits the name at spaces and double quotes,
// so for an output path with such a character makeblastdb ends here with
// exit status 2 ("No alias or index file found ...") after every successful
// build. Apart from that the resolution follows NCBI's rules, the volume of
// the same name in the current directory included.
//
// LOSAT does not write the JSON file that makeblastdb writes next. A file of
// that name left by an earlier database would describe a database that is
// gone (the removal of an existing database takes the JSON file of a single
// volume, but not `NAME.pjs` of a database of several volumes), so it is
// removed where makeblastdb would overwrite it: after the database resolved.
/// The end of `x_BuildDatabase` after a successful build: resolves the new
/// database and removes the JSON metadata file of an earlier one.
fn open_new_database(output_dbname: &str, is_protein: bool) -> Result<(), MakedbError> {
    let seqtype = if is_protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    SeqDb::locate_one_name(output_dbname, seqtype).map_err(MakedbError::Database)?;
    let metadata_filename = format!("{output_dbname}.{}js", seqtype.letter());
    let _ = std::fs::remove_file(metadata_filename);
    Ok(())
}

fn add_fasta(
    db: &mut BuildDatabase<'_>,
    input: &mut dyn Read,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    // The reader's warnings and the serializer's diagnostics arrive in the
    // order makeblastdb posts them.
    let mut warn = warn_to(err);
    db.add_fasta(input, &mut warn).map_err(MakedbError::Build)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1167-1171
// ```c++
//     m_LogFile << "\n\nBuilding a new DB, current time: "
//               << CTime(CTime::eCurrent).AsString() << endl;
//
//     m_LogFile << "New DB name:   " << output_dbname << endl;
//     m_LogFile << "New DB title:  " << title << endl;
// ```
// Every line of the log ends with `endl`, which flushes the stream, so the
// log and the messages of the standard error stream reach a common
// destination (`2>&1`) in the order of the events. The standard output of
// the Rust library is line buffered whatever it is connected to; the log
// goes to it without another buffer.
/// Runs `LOSAT makedb` and returns the exit code.
pub fn run(args: &MakedbArgs) -> i32 {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    let now = LocalTime::now();
    let status = run_with(args, now.as_ref(), &mut stdin.lock(), &mut out, &mut err);
    let _ = out.flush();
    status
}

/// One run with the given time and streams: the build, then the error
/// message and the exit code of a failed build.
pub(crate) fn run_with(
    args: &MakedbArgs,
    now: Option<&LocalTime>,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    match build_database(args, now, stdin, out, err) {
        Ok(()) => 0,
        Err(error) => {
            let _ = writeln!(err, "{error}");
            error.exit_code()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn input_names_split_at_spaces_and_quotes() {
        assert_eq!(split_quoted("a.fa"), ["a.fa"]);
        assert_eq!(split_quoted("a.fa b.fa"), ["a.fa", "b.fa"]);
        assert_eq!(split_quoted("\"a b.fa\""), ["a b.fa"]);
        assert_eq!(split_quoted("\"a b.fa\" c.fa"), ["a b.fa", "c.fa"]);
        assert_eq!(split_quoted("  a.fa  "), ["a.fa"]);
        assert!(split_quoted("").is_empty());
    }

    #[test]
    fn the_variable_that_changes_stored_titles_is_refused() {
        const DONT_SCAN: &str = "NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS";
        let with = |name: &'static str, value: &'static str| {
            move |asked: &str| (asked == name).then(|| std::ffi::OsString::from(value))
        };
        assert!(check_environment(&|_: &str| None).is_ok());
        let error = check_environment(&with(DONT_SCAN, "1")).unwrap_err();
        assert!(error.contains(DONT_SCAN), "{error}");
        // An empty value is no value for makeblastdb.
        assert!(check_environment(&with(DONT_SCAN, "")).is_ok());
        // Read by the source of ASN.1 input only.
        assert!(check_environment(&with("NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS", "1")).is_ok());
    }

    #[test]
    fn errors_end_with_makeblastdb_s_exit_codes_and_prefixes() {
        let io = || std::io::Error::other("no");
        for (error, code, text) in [
            (
                MakedbError::Build(WriteDbError::Options("o".into())),
                1,
                "BLAST options error: o",
            ),
            (
                MakedbError::Build(WriteDbError::Creation("c".into())),
                1,
                "BLAST Database creation error: c",
            ),
            (
                MakedbError::Database(SeqDbError::NotFound("d".into())),
                2,
                "BLAST Database error: d",
            ),
            (
                MakedbError::Build(WriteDbError::OutOfMemory),
                4,
                "BLAST ran out of memory",
            ),
            (
                MakedbError::Build(WriteDbError::Other("x".into())),
                255,
                "Error: x",
            ),
            (
                MakedbError::Build(WriteDbError::Io {
                    path: "p".into(),
                    source: io(),
                }),
                255,
                "Error: p: no",
            ),
            (
                MakedbError::Build(WriteDbError::Rename {
                    from: "a".into(),
                    to: "b".into(),
                    source: io(),
                }),
                255,
                "Error: a -> b: no",
            ),
        ] {
            assert_eq!(error.exit_code(), code, "{error}");
            assert_eq!(error.to_string(), text);
        }
    }
}
