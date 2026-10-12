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

/// The options of makeblastdb 2.17.0 that `LOSAT makedb` refuses
/// (makeblastdb.cpp:246-333).
pub const UNSUPPORTED_OPTIONS: [&str; 14] = [
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

/// How `LOSAT makedb` ends when it does not build a database.
#[derive(Debug)]
pub enum MakedbError {
    /// An error of makeblastdb, with makeblastdb's exit code.
    Build(WriteDbError),
    /// A request makeblastdb serves and LOSAT does not.
    Unsupported(String),
}

impl MakedbError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
    // ```c++
    // #define CATCH_ALL(exit_code)                                                \
    //     catch (const CInvalidDataException& e) {                                \
    //         LOG_POST(Error << "BLAST options error: " << e.GetMsg());           \
    //         exit_code = BLAST_INPUT_ERROR;                                      \
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
    /// The exit code: 1 (`BLAST_INPUT_ERROR`) for an error of the options or
    /// of the database creation, 255 (`BLAST_UNKNOWN_ERROR`) otherwise.
    pub fn exit_code(&self) -> i32 {
        match self {
            MakedbError::Build(WriteDbError::Other(_)) => 255,
            MakedbError::Build(_) | MakedbError::Unsupported(_) => 1,
        }
    }
}

impl std::fmt::Display for MakedbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MakedbError::Build(error) => error.fmt(f),
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
    // makeblastdb reads the value after it has removed an existing database
    // of the same name; LOSAT reads it first, so that a wrong value costs
    // nothing. A value NCBI cannot convert ends makeblastdb with an uncaught
    // exception (exit 255); here it is an error of the options.
    let bytes = string_to_uint8_data_size(&args.max_file_sz)
        .map_err(|error| options_error(format!("-max_file_sz: {}", error.message)))?;
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
// makeblastdb then writes the JSON metadata file (`.pjs` / `.njs`,
// makeblastdb.cpp:1297-1318). NCBI's programs do not read it when they
// search, and LOSAT does not write it.
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
    if dbname == "-" || dbname.split(' ').filter(|part| !part.is_empty()).count() > 1 {
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
                Err(source) => Err(MakedbError::Build(WriteDbError::Io {
                    path: name.into(),
                    source,
                })),
            },
            _ => Err(options_error(format!("File {name} does not exist"))),
        }
    };
    result?;

    let closed = db.end_build(false);
    let _ = err.write_all(&db.take_diagnostics());
    closed?;
    Ok(())
}

fn add_fasta(
    db: &mut BuildDatabase<'_>,
    input: &mut dyn Read,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    // The reader's warnings and the serializer's diagnostics arrive in the
    // order makeblastdb posts them; what the last sequence leaves follows.
    let result = {
        let mut warn = warn_to(err);
        db.add_fasta(input, &mut warn)
    };
    let _ = err.write_all(&db.take_diagnostics());
    result.map_err(MakedbError::Build)
}

/// Runs `LOSAT makedb` and returns the exit code.
pub fn run(args: &MakedbArgs) -> i32 {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = std::io::BufWriter::new(stdout.lock());
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
}
