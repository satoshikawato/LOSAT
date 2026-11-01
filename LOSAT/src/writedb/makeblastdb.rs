//! `losat makeblastdb`: the application layer of makeblastdb
//! (`c++/src/app/blastdb/makeblastdb.cpp`, `CMakeBlastDBApp`) for FASTA
//! input and a database of format version 5 or 4 without parsed Seq-ids.

use std::io::{BufRead, Read, Write};

use clap::Args;

use super::build_db::BuildDatabase;
use super::data_size::{string_to_uint8_data_size, uint8_to_string_data_size};
use super::format_guess::appears_fasta;
use super::metadata::{db_files_metadata, last_updated, write_metadata_file, DbMetadata};
use super::platform::Platform;
use super::taxid_set::{TaxIdMapError, TaxIdSet};
use super::time::LocalTime;
use super::writer::{build_lmdb_file_name, DbVersion, Indices};
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
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:309-312
// ```c++
//     arg_desc->AddOptionalKey("logfile", "File_Name",
//                              "File to which the program log should be redirected",
//                              CArgDescriptions::eOutputFile,
//                              CArgDescriptions::fAppend);
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:246-251
// ```c++
//     arg_desc->AddFlag("parse_seqids",
//                       "Option to parse seqid for FASTA input if set, for all other input types seqids are parsed automatically", true);
//
//     arg_desc->AddFlag("hash_index",
//                       "Create index of sequence hash values.",
//                       true);
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:319-329
// ```c++
//     arg_desc->AddOptionalKey("taxid", "TaxID",
//                              "Taxonomy ID to assign to all sequences",
//                              CArgDescriptions::eInteger);
//     arg_desc->SetConstraint("taxid", new CArgAllowValuesGreaterThanOrEqual(0));
//     arg_desc->SetDependency("taxid", CArgDescriptions::eExcludes, "taxid_map");
//
//     arg_desc->AddOptionalKey("taxid_map", "TaxIDMapFile",
//              "Text file mapping sequence IDs to taxonomy IDs.\n"
//              "Format:<SequenceId> <TaxonomyId><newline>",
//              CArgDescriptions::eInputFile);
//     arg_desc->SetDependency("taxid_map", CArgDescriptions::eRequires, "parse_seqids");
// ```
// LOSAT's parser checks the range of `-taxid`, the exclusion and the requirement; its errors
// are the parser's (exit status 2, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES) where NCBI's argument
// reader prints its usage (exit status 1). `check_arguments` repeats the two dependencies for
// the library. The options of makeblastdb that LOSAT does not support are not registered;
// `crate::cli` names them in its error.
/// Arguments of `losat makeblastdb`, under makeblastdb's names.
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
    /// Version of BLAST database to be created
    #[arg(long, default_value_t = 5, value_name = "version",
          value_parser = clap::value_parser!(i32).range(4..=5))]
    pub blastdb_version: i32,
    /// Maximum file size for BLAST database files
    #[arg(long, default_value = "3GB", value_name = "number_of_bytes")]
    pub max_file_sz: String,
    /// File to which the program log should be redirected
    #[arg(long, value_name = "File_Name")]
    pub logfile: Option<String>,
    /// Create index of sequence hash values.
    #[arg(long)]
    pub hash_index: bool,
    /// Taxonomy ID to assign to all sequences
    #[arg(long, value_name = "TaxID", value_parser = clap::value_parser!(i32).range(0..),
          conflicts_with = "taxid_map")]
    pub taxid: Option<i32>,
    /// Text file mapping sequence IDs to taxonomy IDs. Format:<SequenceId> <TaxonomyId><newline>
    #[arg(long, value_name = "TaxIDMapFile", requires = "parse_seqids")]
    pub taxid_map: Option<String>,
    /// Option to parse seqid for FASTA input if set, for all other input types seqids are
    /// parsed automatically
    #[arg(long)]
    pub parse_seqids: bool,
}

/// The options of makeblastdb 2.17.0 that `losat makeblastdb` refuses
/// (makeblastdb.cpp:246-333). `-h` and `-version` are served (`crate::cli`);
/// `-logfile`, `-parse_seqids`, `-hash_index`, `-taxid` and `-taxid_map` are read.
pub const UNSUPPORTED_OPTIONS: [&str; 8] = [
    "randomize",
    "mask_data",
    "mask_id",
    "mask_desc",
    "gi_mask",
    "gi_mask_name",
    "metadata_output_prefix",
    "oid_masks",
];

/// How `losat makeblastdb` ends when it does not build a database, or cannot
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
    /// NCBI's `CArgException` raised after the arguments were parsed: the text
    /// makeblastdb prints after "Command line argument error: ".
    Argument(String),
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
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
    // ```c++
    //     catch (const CArgException& e) {                                        \
    //         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
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
            MakedbError::Build(
                WriteDbError::Creation(_)
                | WriteDbError::CreationBytes(_)
                | WriteDbError::Options(_),
            )
            | MakedbError::Unsupported(_)
            | MakedbError::Argument(_) => 1,
        }
    }
}

impl std::fmt::Display for MakedbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MakedbError::Build(error) => error.fmt(f),
            MakedbError::Database(error) => write!(f, "BLAST Database error: {error}"),
            MakedbError::Unsupported(message) => write!(f, "losat makeblastdb: {message}"),
            MakedbError::Argument(message) => write!(f, "Command line argument error: {message}"),
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
             titles for some FASTA definition lines, is not supported by LOSAT's makeblastdb"
        ));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:739-747
// ```c++
//     CNcbiApplication* app = CNcbiApplication::Instance();
//     if (app) {
//         const CNcbiRegistry& registry = app->GetConfig();
//         const string& value = registry.Get("BLAST", "MAX_SEQID_LENGTH");
//         if (!value.empty()) {
//             try { m_FastaReader->SetMaxIDLength(NStr::StringToUInt(value)); }
//             catch (const exception&) {}  // Ignore errors
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1211-1216
// ```c++
//     bool long_seqids = false;
//     CNcbiApplication* app = CNcbiApplication::Instance();
//     if (app) {
//         const CNcbiRegistry& registry = app->GetConfig();
//         long_seqids = (registry.Get("BLAST", "LONG_SEQID") == "1");
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:485-491,647-653
// ```c++
// void CFastaReader::SetMaxIDLength(Uint4 max_len)
// {
//     CFastaDeflineReader::s_MaxLocalIDLength =
//     CFastaDeflineReader::s_MaxGeneralTagLength =
//     CFastaDeflineReader::s_MaxAccessionLength = m_MaxIDLength = max_len;
//     m_bModifiedMaxIdLength=true;
// }
// ...
// void CFastaReader::x_SetDeflineParseInfo(SDefLineParseInfo& info)
// {
//     info.fBaseFlags = m_iFlags;
//     info.fFastaFlags = GetFlags();
//     info.maxIdLength = m_bModifiedMaxIdLength ?
//                        m_MaxIDLength :
//                        0;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:430-434
// ```c++
//     if (info.maxIdLength) {
//         s_IdValidate.SetMaxLocalIDLength(info.maxIdLength);
//         s_IdValidate.SetMaxGeneralTagLength(info.maxIdLength);
//         s_IdValidate.SetMaxAccessionLength(info.maxIdLength);
//     }
// ```
// The registry's limit of a Seq-id's length changes which FASTA ids makeblastdb parses with
// `-parse_seqids` when `NStr::StringToUInt` converts the value and the value is not 0 (0 is
// the reader's "not set"; the static `s_Max*` copies are read nowhere). A value it cannot
// convert is ignored (`catch`), as is an empty one. LOSAT's FASTA reader keeps NCBI's default
// limits (makedb-v5 DECISIONS 12), so it refuses the entry with `-parse_seqids` exactly when
// makeblastdb would use it. `[BLAST] LONG_SEQID` is refused for every program by
// `check_ncbi_application_settings` (it is not one of the entries that change nothing).
/// Refuses `[BLAST] MAX_SEQID_LENGTH` (`max_seqid_length`, from
/// `ApplicationSettings::max_seqid_length`: the value of the first layer that has the entry,
/// `Err` for a file value LOSAT cannot read) for a build with `-parse_seqids` when makeblastdb
/// would use the value.
pub fn check_registry_settings(
    args: &MakedbArgs,
    max_seqid_length: Result<Option<&str>, &str>,
) -> Result<(), String> {
    if !args.parse_seqids {
        return Ok(());
    }
    match max_seqid_length {
        Err(reason) => Err(format!(
            "{reason} ([BLAST] MAX_SEQID_LENGTH changes the Seq-ids NCBI makeblastdb parses \
             with -parse_seqids)"
        )),
        Ok(Some(value)) if string_to_uint(value).is_some_and(|limit| limit != 0) => Err(format!(
            "the registry entry [BLAST] MAX_SEQID_LENGTH ({value}), which changes the Seq-ids \
             NCBI makeblastdb parses with -parse_seqids, is not supported by LOSAT's makeblastdb"
        )),
        Ok(_) => Ok(()),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:646-656
// ```c++
// unsigned int
// NStr::StringToUInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Uint8 value = StringToUInt8(str, flags, base);
//     if ( value > kMax_UInt ) {
//         S2N_CONVERT_ERROR(unsigned int, "overflow", ERANGE, 0);
//     }
//     return (unsigned int) value;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:884-918
// ```c++
//     if ( base == 10  &&  (flags & slow_flags) == 0 ) {
//         // fast conversion
//
//         // Current position in the string
//         CTempString::const_iterator ptr = str.begin(), end = str.end();
//
//         // Determine sign
//         if ( ptr != end && *ptr == '+' ) {
//             ++ptr;
//         }
//         if ( ptr == end ) {
//             S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
//         }
//
//         // Begin conversion
//         Uint8 n = 0;
//
//         const Uint8 limdiv = kMax_UI8/10;
//         const int   limoff = int(kMax_UI8 % 10);
//
//         do {
//             char ch = *ptr;
//             int  delta = ch - '0';
//             if ( unsigned(delta) >= 10 ) {
//                 S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
//             }
//             // Overflow check
//             if ( n >= limdiv && (n > limdiv || delta > limoff) ) {
//                 S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, ptr-str.begin());
//             }
//             n = n*10+delta;
//         } while ( ++ptr != end );
//
//         return n;
//     }
// ```
/// `NStr::StringToUInt(value)` with no flags: `None` where it throws.
fn string_to_uint(value: &str) -> Option<u32> {
    let digits = value.strip_prefix('+').unwrap_or(value);
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let mut n: u64 = 0;
    for c in digits.bytes() {
        n = n.checked_mul(10)?.checked_add(u64::from(c - b'0'))?;
    }
    u32::try_from(n).ok()
}

/// The checks of the arguments that LOSAT makes first: what LOSAT does not
/// support, and the format version.
fn check_arguments(args: &MakedbArgs) -> Result<DbVersion, MakedbError> {
    if args.input_type != "fasta" {
        return Err(unsupported(format!(
            "-input_type {} is not supported; only fasta input is",
            args.input_type
        )));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1218-1219
    // ```c++
    //     const EBlastDbVersion dbver =
    //         static_cast<EBlastDbVersion>(args["blastdb_version"].AsInteger());
    // ```
    // The parser allows 4 and 5 only (`CArgAllow_Integers(eBDB_Version4, eBDB_Version5)`).
    let Some(version) = DbVersion::from_number(args.blastdb_version) else {
        return Err(unsupported(format!(
            "-blastdb_version {} is not a version makeblastdb writes (4 or 5)",
            args.blastdb_version
        )));
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:322-329
    // ```c++
    //     arg_desc->SetConstraint("taxid", new CArgAllowValuesGreaterThanOrEqual(0));
    //     arg_desc->SetDependency("taxid", CArgDescriptions::eExcludes, "taxid_map");
    //
    //     arg_desc->AddOptionalKey("taxid_map", "TaxIDMapFile",
    //              "Text file mapping sequence IDs to taxonomy IDs.\n"
    //              "Format:<SequenceId> <TaxonomyId><newline>",
    //              CArgDescriptions::eInputFile);
    //     arg_desc->SetDependency("taxid_map", CArgDescriptions::eRequires, "parse_seqids");
    // ```
    // The parser checks these on the command line (`-taxid` has its range there); a caller of
    // the library that sets the fields is checked here with the parser's sentences.
    if args.taxid_map.is_some() {
        if !args.parse_seqids {
            return Err(MakedbError::Argument(
                "Argument \"parse_seqids\". Must be specified, as it is required by argument:  \
                 `taxid_map'"
                    .to_string(),
            ));
        }
        if args.taxid.is_some() {
            return Err(MakedbError::Argument(
                "Argument \"taxid\". Incompatible with argument:  `taxid_map'".to_string(),
            ));
        }
    }
    Ok(version)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1278-1287
// ```c++
//     if (args["taxid"].HasValue()) {
//         _ASSERT( !args["taxid_map"].HasValue() );
//         CRef<CTaxIdSet> taxids(new CTaxIdSet(TAX_ID_FROM(int, args["taxid"].AsInteger())));
//         m_DB->SetTaxids(*taxids);
//     } else if (args["taxid_map"].HasValue()) {
//         _ASSERT( !args["taxid"].HasValue() );
//         CRef<CTaxIdSet> taxids(new CTaxIdSet());
//         taxids->SetMappingFromFile(args["taxid_map"].AsInputFile());
//         m_DB->SetTaxids(*taxids);
//     }
// ```
// `AsInputFile` opens the file here: one that cannot be opened is the parser's error
// "File is not accessible" (makeblastdb 2.17.0: "Command line argument error: Argument
// "taxid_map". File is not accessible:  `<name>'", exit status 1, after the head of the
// log). A number the map cannot convert ends makeblastdb with its exception report
// (`CATCH_ALL`, "Error: ", exit status 255; the build tree in its path is the official
// release's, makedb-v5 DECISIONS 29).
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-737
// ```c++
//     if (AsString() == "-") {
// #if defined(NCBI_OS_MSWIN)
//         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
// #endif
//         m_Ios  = &cin;
//     } else if ( !AsString().empty() ) {
//         if (!fstrm) {
//             fstrm = new CNcbiIfstream;
//         }
//         if (fstrm) {
//             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
//                 fstrm = NULL;
//             } else {
//                 m_DeleteFlag = true;
//             }
//         }
//         m_Ios = fstrm;
//     }
//     CArg_Ios::x_Open(flags);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-621,95-99
// ```c++
// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
// {
//     if ( !m_Ios ) {
//         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
//             "File is not accessible",AsString()));
//     }
// }
// string s_ArgExptMsg(const string& name, const string& what, const string& attr)
// {
//     return string("Argument \"") + (name.empty() ? s_ExtraName : name) +
//         "\". " + what + (attr.empty() ? attr : ":  `" + attr + "'");
// }
// ```
// `-` is the standard input: the map is read before the FASTA input, so with `-in -` as well
// the map takes all of it and the FASTA reader finds none ("No volumes were created.", as
// makeblastdb). The empty name opens nothing and its message has no name part.
/// The taxids of `-taxid` or `-taxid_map` (`SetTaxids`), read where makeblastdb reads them;
/// `stdin` is read for `-taxid_map -`, and the second value is then `Some(failed)`: whether
/// the standard input has failed (its state when `-in -` reads it next).
fn read_taxids(
    args: &MakedbArgs,
    stdin: &mut dyn Read,
    err: &mut dyn Write,
) -> Result<(Option<TaxIdSet>, Option<bool>), MakedbError> {
    if let Some(taxid) = args.taxid {
        return Ok((Some(TaxIdSet::with_global(taxid)), None));
    }
    let Some(name) = args.taxid_map.as_deref() else {
        return Ok((None, None));
    };
    let not_accessible = || {
        let attr = if name.is_empty() {
            String::new()
        } else {
            format!(":  `{name}'")
        };
        MakedbError::Argument(format!(
            "Argument \"taxid_map\". File is not accessible{attr}"
        ))
    };
    let mut file;
    let mut standard_input;
    let reader: &mut dyn BufRead = if name == "-" {
        standard_input = std::io::BufReader::new(stdin);
        &mut standard_input
    } else if name.is_empty() {
        return Err(not_accessible());
    } else {
        file = std::io::BufReader::new(std::fs::File::open(name).map_err(|_| not_accessible())?);
        &mut file
    };
    let mut taxids = TaxIdSet::default();
    let mut warn = warn_to(err);
    let failed = taxids
        .set_mapping_from_file(reader, &mut |warning| {
            super::fasta_source::defline_ids::post_seq_id_warning(warning, &mut warn)
        })
        .map_err(|error| match error {
            TaxIdMapError::Exception(report) => MakedbError::Build(WriteDbError::Other(report)),
            TaxIdMapError::BlankLine(line) => unsupported(format!(
                "line {line} of the -taxid_map file {name} holds only white space, which \
                 makeblastdb 2.17.0 cannot read (it crashes)"
            )),
        })?;
    Ok((Some(taxids), (name == "-").then_some(failed)))
}

/// The value of `-max_file_sz`, read before the build touches a file.
fn max_file_size(args: &MakedbArgs) -> Result<u64, MakedbError> {
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
    // before (after the checks of the names, the title and `-logfile`), so
    // that a wrong value costs nothing. A value NCBI cannot convert
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

/// A name of `-in` that makeblastdb reads through a decompressor.
fn is_compressed_name(name: &str) -> bool {
    let lower = name.to_ascii_lowercase();
    [".gz", ".bz2", ".zst"]
        .iter()
        .any(|ext| lower.ends_with(ext))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:792-836
// ```c++
// void
// CMakeBlastDBApp::x_VerifyInputFilesType(const vector<CTempString>& filenames,
//                                    CMakeBlastDBApp::ESupportedInputFormats input_type)
// {
//     //Let other part of the program deal with blastdb input
//     if(eBlastDb == input_type)
//         return;
//
//     // Guess the input data type
//     for (size_t i = 0; i < filenames.size(); i++) {
//         string seq_file = filenames[i];
//
//         CFile input_file(seq_file);
//         if ( !input_file.Exists() ) {
//             string error_msg = "File " + seq_file + " does not exist";
//             NCBI_THROW(CInvalidDataException, eInvalidInput, error_msg);
//         }
//         if (input_file.GetLength() == 0) {
//             string error_msg = "File " + seq_file + " is empty";
//             NCBI_THROW(CInvalidDataException, eInvalidInput, error_msg);
//         }
//
//         CNcbiIfstream f(seq_file.c_str(), ios::binary);
//         ...
//         if(input_type == eFasta && x_ConvertToSupportedType(x_GuessFileType( *data_in ) ) != eFasta) {
//             string msg = "\nInput file " + seq_file + " does NOT appear to be FASTA (processing anyway).\n" \
//                  + "Advise validating database with 'blastdbcheck -dbtype [prot|nucl] -db ${DBNAME}'\n";
//             ERR_POST(Warning << msg);
//         }
//     }
//     return;
// }
// ```
// Every file is checked, in the order of `-in`, before the first one is read.
// `CFile::Exists` is true for a regular file, also through a link; anything
// else "does not exist". The test of the format opens the file again, so the
// build reads it from its start. `ERR_POST` ends the message with a line feed.
/// `x_VerifyInputFilesType` for FASTA input.
fn verify_input_files(names: &[&str], err: &mut dyn Write) -> Result<(), MakedbError> {
    for &seq_file in names {
        match std::fs::metadata(seq_file) {
            Ok(meta) if meta.is_file() && meta.len() == 0 => {
                return Err(options_error(format!("File {seq_file} is empty")));
            }
            Ok(meta) if meta.is_file() => {}
            _ => return Err(options_error(format!("File {seq_file} does not exist"))),
        }
        let looks_like_fasta = std::fs::File::open(seq_file).is_ok_and(appears_fasta);
        if !looks_like_fasta {
            let _ = write!(
                err,
                "Warning: [makeblastdb] \nInput file {seq_file} does NOT appear to be \
                 FASTA (processing anyway).\nAdvise validating database with \
                 'blastdbcheck -dbtype [prot|nucl] -db ${{DBNAME}}'\n\n"
            );
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1146-1201
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
//
//     m_LogFile = & (args["logfile"].HasValue()
//                    ? args["logfile"].AsOutputFile()
//                    : cout);
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
// The title defaults to the whole value of `-in`, several names included. After a
// successful build makeblastdb opens the new database with its reader and writes
// the JSON metadata file (`open_new_database`); an error of that open ends
// makeblastdb with "BLAST Database error" and exit status 2 although the database
// was written.
/// `CMakeBlastDBApp::x_BuildDatabase`. The log goes to `out` (makeblastdb's
/// standard output, or the file of `-logfile`), warnings and diagnostics to
/// `err`; `stdin` is read when `-in` is `-`. `now` is the local time for the
/// database's date (`None` where the platform has no clock), `platform` the
/// platform whose makeblastdb the files follow ([`Platform::host`] for a build).
pub fn build_database(
    args: &MakedbArgs,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let mut end_message = None;
    let result =
        build_database_with_log(args, now, platform, stdin, out, err, None, &mut end_message);
    if let Some(message) = end_message {
        let _ = writeln!(err, "{message}");
    }
    result
}

/// [`build_database`]; `logfile_error` is the error of a `-logfile` that could
/// not be opened, raised where makeblastdb opens the file. `end_message` gets the
/// message of `~CBuildDatabase` (`-taxid_map` matched no sequence), which makeblastdb
/// posts after everything else, the error of the build included.
#[allow(clippy::too_many_arguments)]
fn build_database_with_log(
    args: &MakedbArgs,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
    logfile_error: Option<MakedbError>,
    end_message: &mut Option<&'static str>,
) -> Result<(), MakedbError> {
    let version = check_arguments(args)?;
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
    if let Some(error) = logfile_error {
        return Err(error);
    }
    let max_file_size = max_file_size(args)?;

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1016-1028
    // ```c++
    // void CMakeBlastDBApp::x_ProcessInputData(const string & paths,
    //                                          bool           is_protein)
    // {
    //     vector<CTempString> names;
    //     SeqDB_SplitQuoted(paths, names);
    //     ...
    //     if(eBlastDb != input_fmt) {
    //         if (names[0] == "-") {
    //             x_AddSequenceData(cin, build_fmt);
    //         }
    // ```
    // A first name `-` reads the standard input and nothing else. With no name at
    // all (`-in ""`) NCBI reads `names[0]` of an empty vector; LOSAT refuses that.
    let names = split_quoted(&args.input);
    if names.is_empty() {
        return Err(unsupported(
            "-in names no input file; losat makeblastdb reads FASTA files or the standard input",
        ));
    }
    let from_stdin = names[0] == "-";
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050
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
    if !from_stdin {
        if let Some(name) = names.iter().find(|name| is_compressed_name(name)) {
            return Err(unsupported(format!(
                "{name}: makeblastdb reads a file of this name as compressed; losat makeblastdb \
                 reads uncompressed FASTA only"
            )));
        }
    }

    let Some(now) = now else {
        return Err(unsupported(
            "this platform has no clock for the database's date",
        ));
    };
    let mut db = BuildDatabase::new(
        dbname,
        title,
        is_protein,
        out,
        now,
        version,
        platform,
        Indices {
            full: args.parse_seqids,
            hash: args.hash_index,
        },
    )?;
    let _ = writeln!(
        db.log(),
        "Maximum file size: {}",
        uint8_to_string_data_size(max_file_size)
    );
    db.set_max_file_size(max_file_size);
    let (taxids, stdin_map) = read_taxids(args, stdin, err)?;
    if let Some(taxids) = taxids {
        db.set_taxids(taxids);
    }
    let result = add_input_and_close(
        &mut db,
        &names,
        from_stdin.then_some(stdin_map),
        title,
        is_protein,
        now,
        stdin,
        err,
    );
    *end_message = db.taxid_message();
    result
}

/// The rest of `x_BuildDatabase`: the input files, `EndBuild` and the metadata file.
/// `from_stdin` is `Some(map)` for `-in -`, `map` being `Some(failed)` when `-taxid_map -`
/// read the standard input first (`failed`: whether it left the stream failed).
#[allow(clippy::too_many_arguments)]
fn add_input_and_close(
    db: &mut BuildDatabase<'_>,
    names: &[&str],
    from_stdin: Option<Option<bool>>,
    title: &str,
    is_protein: bool,
    now: &LocalTime,
    stdin: &mut dyn Read,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    if let Some(stdin_map) = from_stdin {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1398-1420
        // ```c++
        // bool CBuildDatabase::AddFasta(CNcbiIstream & fasta_file)
        // {
        //     // Add any fasta sequences as well.
        //     bool success = true;
        //
        //     if (fasta_file) {
        //         CFastaBioseqSource fbs(fasta_file,
        // ...
        //     }
        //     return success;
        // }
        // ```
        // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:619-636
        // ```c++
        // bool CBufferedLineReader::x_ReadBuffer()
        // {
        //     _ASSERT(m_Reader);
        //
        //     if ( m_Eof ) {
        //         return false;
        //     }
        //
        //     m_InputPos += CT_OFF_TYPE(m_End - m_Buffer.get());
        //     m_Pos = m_End = m_Buffer.get();
        //     for (bool flag = true; flag; ) {
        //         size_t size;
        //         ERW_Result result =
        //             m_Reader->Read(m_Buffer.get(), m_BufferSize, &size);
        //         switch (result) {
        //         case eRW_NotImplemented:
        //         case eRW_Error:
        //             NCBI_THROW(CIOException, eRead, "Read error");
        // ```
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
        // ```c++
        //     catch (const CIOException& e) {                                         \
        //         if (e.GetErrCode() == CIOException::eFlush) {                       \
        //             LOG_POST(Error << "BLAST failed to write output: " << e.GetMsg());\
        //             exit_code = BLAST_OUTPUT_ERROR;                                 \
        //         }                                                                   \
        //     }                                                                       \
        // ```
        // After `-taxid_map -` has read the standard input to its end (a read that found
        // nothing) or to an error, `cin` has failed: makeblastdb adds nothing from it and
        // reports no error ("No volumes were created.", exit status 0). A map whose last
        // line has no line feed leaves `cin` at its end with only `eofbit` set: the line
        // reader of `CFastaBioseqSource`, made before `AddFasta`'s `try`, cannot read it and
        // throws `CIOException::eRead`, which `CATCH_ALL` swallows: the build stops there,
        // without `EndBuild` and without a message, exit status 0 (makeblastdb 2.17.0, case
        // T_map_only_stdin_noeol; a deterministic NCBI result, reproduced).
        match stdin_map {
            None => add_fasta(db, stdin, err)?,
            Some(false) => return Ok(()),
            // A failed stream: `AddFasta` skips it.
            Some(true) => {}
        }
    } else {
        verify_input_files(names, err)?;
        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1030-1054
        // ```c++
        //              x_VerifyInputFilesType(names, input_fmt);
        //             for (size_t i = 0; i < names.size(); i++) {
        //                 string seq_file = names[i];
        //                 CNcbiIfstream f(seq_file.c_str(), ios::binary);
        // 		CNcbiIstream *data_in = &f;  // uncompressed by default
        //         ...
        // 		// common part
        // 		x_AddSequenceData( *data_in , build_fmt);
        //             }
        // ```
        // Each file goes to its own `AddFasta` (`x_AddSequenceData` ->
        // `x_AddFasta`): its own FASTA reader (line numbers and local ids start
        // again), its own "Adding sequences" line, and "No sequences added" for a
        // file without a sequence.
        for &name in names {
            match std::fs::File::open(name) {
                Ok(mut file) => add_fasta(db, &mut file, err)?,
                // A file that exists and cannot be opened. makeblastdb does
                // not test the stream it opened: `AddFasta` skips a stream
                // that is not good and adds nothing from it (a build of that
                // file alone ends with "No volumes were created." and exit
                // status 0). LOSAT does not reproduce that: it reports the
                // file, with the exit status of its file-system errors, and
                // removes what earlier files added.
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
                Err(source) => {
                    let mut warn = warn_to(err);
                    let _ = db.end_build(true, &mut warn);
                    return Err(MakedbError::Build(WriteDbError::Io {
                        path: name.into(),
                        source,
                    }));
                }
            }
        }
    }

    let closed = {
        let mut warn = warn_to(err);
        db.end_build(false, &mut warn)
    };
    if closed? {
        open_new_database(db, title, is_protein, now)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1296-1328
// ```c++
//     bool success = m_DB->EndBuild();
//     if(success) {
//     	string new_db = m_DB->GetOutputDbName();
//     	CSeqDB::ESeqType t = is_protein? CSeqDB::eProtein: CSeqDB::eNucleotide;
//     	CSeqDB sdb(new_db, t);
//         string output_prefix = args["metadata_output_prefix"]
//                 ? args["metadata_output_prefix"].AsString()
//                 : kEmptyStr;
//         ...
//     	CRef<CBlast_db_metadata> m = sdb.GetDBMetaData(output_prefix);
//     	string extn (kEmptyStr);
//     	SeqDB_GetMetadataFileExtension(is_protein, extn);
//     	string metadata_filename = new_db + "." + extn;
//         ofstream out(metadata_filename.c_str());
//         ...
//         out << NcbiEndl;
//     }
// ```
// `CSeqDB sdb(new_db, t)` throws a `CSeqDBException` when the reader cannot
// open what was written; the usual cause is an alias file of an earlier
// database of the same name that `DeleteBlastDb` could not resolve and
// therefore left in place, which now hides the new single volume. LOSAT
// resolves the database here for that outcome and for the values of the
// metadata, and it only resolves: the alias files are read and every volume is
// located, but no volume is opened. The volumes were written a moment ago, and
// LOSAT's reader holds the index of an open volume in memory (12 bytes for each
// sequence). So the error `CSeqDB` gives for more volumes than the process may
// hold open files ("Too many open files, please raise the open file limit",
// exit status 2 after a complete build) is not reproduced. The totals of the
// metadata are those of the volumes written (what the reader would read from
// their index files), and their date is the build's.
//
// The path is one database name (`SeqDb::locate_one_name_metadata`), which
// departs from NCBI on purpose: `CSeqDB` splits the name at spaces and double
// quotes, so for an output path with such a character makeblastdb ends here with
// exit status 2 ("No alias or index file found ...") after every successful
// build. Apart from that the resolution follows NCBI's rules, the volume of the
// same name in the current directory included. For the same reason the `dbname`
// of the metadata is the file name of the one path.
/// The end of `x_BuildDatabase` after a successful build: resolves the new
/// database and writes its metadata file.
fn open_new_database(
    db: &BuildDatabase<'_>,
    title: &str,
    is_protein: bool,
    now: &LocalTime,
) -> Result<(), MakedbError> {
    let new_db = db.output_dbname();
    let seqtype = if is_protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let located =
        SeqDb::locate_one_name_metadata(new_db, seqtype, &|_| Some(title.as_bytes().to_vec()))
            .map_err(MakedbError::Database)?;
    let output_db = db.output_db();
    let version = output_db.version();
    let totals = output_db.volume_totals();
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1705-1747
    // ```c++
    // 	GetTotals(CSeqDB::eFilteredAll, &num_seqs, &total_length, true);
    // 	vector<string> dblist;
    //     NStr::Split(GetDBNameList(), " ", dblist, NStr::fSplit_Tokenize);
    //     NON_CONST_ITERATE(vector<string>, itr, dblist) {
    //         size_t off = (*itr).find_last_of(CFile::GetPathSeparator());
    //         if (off != string::npos ) {
    //         	(*itr).erase(0, off+1);
    //         }
    //     }
    //
    // 	string dbnames = NStr::Join(dblist, " ");
    // 	m->SetDbname(dbnames);
    //
    //     m->SetDbtype(GetSequenceType() == CSeqDB::eProtein ? "Protein" : "Nucleotide" );
    // 	m->SetDb_version(GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5?5:4);
    // 	m->SetDescription(GetTitle());
    // 	m->SetNumber_of_letters(total_length);
    // 	m->SetNumber_of_sequences(num_seqs);
    //     ...
    // 	Int8 disk_bytes(0), cached_bytes(0);
    // 	x_GetDBFilesMetaData(disk_bytes, cached_bytes, m->SetFiles(), user_path);
    // 	m->SetBytes_total(disk_bytes);
    // 	m->SetBytes_to_cache(cached_bytes);
    //
    // 	m->SetNumber_of_volumes(m_Impl->GetNumOfVols());
    //
    // 	if (GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) {
    // 		set<TTaxId> tax_ids;
    // 		GetDBTaxIds(tax_ids);
    // 		if((tax_ids.size() > 1) || ((tax_ids.size() == 1) && (0 != *tax_ids.begin()))){
    // 			m->SetNumber_of_taxids(static_cast<int>(tax_ids.size()));
    // 		}
    // 	}
    // ```
    // Without `-taxid` and `-taxid_map` every OID has tax id 0, so
    // `number-of-taxids` is not set. `GetDBTaxIds` reads the taxids of the tax id files,
    // which are the ones the writer gave them (`WriteDb::db_tax_ids`).
    let dbname = match new_db.rfind(std::path::MAIN_SEPARATOR) {
        Some(off) => &new_db[off + 1..],
        None => new_db,
    };
    let lmdb: Vec<String> = (version == DbVersion::V5)
        .then(|| build_lmdb_file_name(new_db, is_protein))
        .into_iter()
        .collect();
    let (bytes_total, bytes_to_cache, files) = db_files_metadata(
        &located.volumes,
        &located.alias_files,
        &lmdb,
        is_protein,
        "",
    );
    let metadata = DbMetadata {
        dbname: dbname.as_bytes().to_vec(),
        dbtype: if is_protein { "Protein" } else { "Nucleotide" },
        db_version: version.number(),
        description: located.title,
        number_of_letters: totals.iter().map(|(_, _, letters)| letters).sum(),
        // NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.def:11-12
        // ```
        // number-of-letters._type = Uint8
        // number-of-sequences._type = Int8
        // ```
        // NCBI's `GetDBMetaData` counts in an `int num_seqs` (seqdb.cpp:1707), which a
        // database of more than 2^31 - 1 sequences would overflow; LOSAT sums in `Int8`.
        number_of_sequences: totals.iter().map(|(_, oids, _)| i64::from(*oids)).sum(),
        last_updated: last_updated(now),
        number_of_volumes: located.volumes.len() as i32,
        number_of_taxids: (version == DbVersion::V5)
            .then(|| output_db.db_tax_ids())
            .filter(|tax_ids| {
                tax_ids.len() > 1 || (tax_ids.len() == 1 && tax_ids.first() != Some(&0))
            })
            .map(|tax_ids| tax_ids.len() as i32),
        bytes_total,
        bytes_to_cache,
        files,
    };
    write_metadata_file(new_db, is_protein, &metadata).map_err(MakedbError::Build)
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
/// Runs `losat makeblastdb` and returns the exit code.
pub fn run(args: &MakedbArgs) -> i32 {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    let now = LocalTime::now();
    let status = run_with(
        args,
        now.as_ref(),
        Platform::host(),
        &mut stdin.lock(),
        &mut out,
        &mut err,
    );
    let _ = out.flush();
    status
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:843-871
// ```c++
//     // Check if logfile is set in the args.
//     m_LogFile = 0;
//     if (!m_DisableArgDesc && argc > 1  &&  argv  &&  diag != eDS_User) {
//         for (int i = 1;  i < argc;  i++) {
//             ...
//             if ( NStr::strcmp(argv[i], s_ArgLogFile) == 0 ) {
//                 if (!argv[++i]) {
//                     continue;
//                 }
//                 m_LogFile = argv[i];
//             } else if (NStr::StartsWith(argv[i], s_ArgLogFile)) {
//                 const char *a = argv[i] + strlen(s_ArgLogFile);
//                 if (*a == '=') {
//                     m_LogFile = ++a;
//                 }
//             }
//         }
//     }
//
//     // Setup logging as soon as possible.
//     // Setup for diagnostics
//     try {
//         CDiagContext::SetupDiag(diag, 0, eDCM_NoChange, m_LogFile);
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-620
// ```c++
// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
// {
//     ...
//         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
//                    "File is not accessible", m_Path));
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7593-7601 (SetLogFile)
// ```c++
//         // Check special filenames
//         if ( file_name.empty()  ||  file_name == "/dev/null" ) {
//             // no output
//             SetDiagStream(0, quick_flush, 0, 0, kLogName_None);
//         }
//         else if (file_name == "-") {
//             // output to stderr
//             SetDiagStream(&NcbiCerr, quick_flush, 0, 0, kLogName_Stderr);
//         }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:776-792 (CArg_OutputFile::x_Open)
// ```c++
//     if (AsString() == "-") {
// #if defined(NCBI_OS_MSWIN)
//         NcbiSys_setmode(NcbiSys_fileno(stdout), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
// #endif
//         m_Ios = &cout;
//     } else if ( !AsString().empty() ) {
//         if (!fstrm) {
//             fstrm = new CNcbiOfstream;
//         }
//         if (fstrm) {
//             if (x_CreatePath(m_CurrentFlags)) {
//                 fstrm->open(AsString().c_str(),IOS_BASE::out | mode);
//             }
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
//                 fstrm = NULL;
// ```
// With `-logfile F` makeblastdb sends its diagnostics (warnings, errors and the
// final error message) to F from the start, and its log to F as well
// (`AsOutputFile` with `fAppend`, makeblastdb.cpp:1199-1201): both append, in the
// order they are written, and nothing goes to the standard streams (two runs
// append to one file). A file that cannot be opened leaves the diagnostics on the
// standard error stream, and the build stops where makeblastdb opens the log,
// after the checks of the database name and the title, with "Command line
// argument error: Argument "logfile". File is not accessible:  `F'" and exit
// status 1. LOSAT opens the file once and writes both streams through the one
// handle. Three names are special. `-` keeps the diagnostics on the standard
// error stream and opens the standard output for the log: makeblastdb runs as
// without `-logfile`. `/dev/null` (that name, on every platform) and the empty
// name discard the diagnostics; the log then goes to the file `/dev/null`, and
// the empty name is no file, so its build stops with the error above, which is a
// diagnostic and so is not shown (makeblastdb 2.17.0: `-logfile ""` prints
// nothing, exit status 1). LOSAT reads `-logfile` as an option of makeblastdb
// only, not anywhere in the command line as NCBI's application does (for example
// as the value of `-title`).
/// One run with the given time, platform and streams: the build, then the error
/// message and the exit code of a failed build. With `-logfile` the log and the
/// messages go where makeblastdb sends them instead of `out` and `err`.
pub(crate) fn run_with(
    args: &MakedbArgs,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    let mut end_message = None;
    let name = match args.logfile.as_deref() {
        None | Some("-") => {
            let result = build_database_with_log(
                args,
                now,
                platform,
                stdin,
                out,
                err,
                None,
                &mut end_message,
            );
            return finish(result, end_message, err);
        }
        Some(name) => name,
    };
    let silent = name.is_empty() || name == "/dev/null";
    let opened = if name.is_empty() {
        None
    } else {
        std::fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(name)
            .ok()
    };
    let mut discard = std::io::sink();
    let diag: &mut dyn Write = if silent { &mut discard } else { err };
    match opened {
        Some(file) if silent => {
            let mut log = &file;
            let result = build_database_with_log(
                args,
                now,
                platform,
                stdin,
                &mut log,
                diag,
                None,
                &mut end_message,
            );
            finish(result, end_message, diag)
        }
        Some(file) => {
            let (mut log, mut diag) = (&file, &file);
            let result = build_database_with_log(
                args,
                now,
                platform,
                stdin,
                &mut log,
                &mut diag,
                None,
                &mut end_message,
            );
            finish(result, end_message, &mut diag)
        }
        None => {
            let error = MakedbError::Argument(format!(
                "Argument \"logfile\". File is not accessible:  `{name}'"
            ));
            let result = build_database_with_log(
                args,
                now,
                platform,
                stdin,
                out,
                diag,
                Some(error),
                &mut end_message,
            );
            finish(result, end_message, diag)
        }
    }
}

/// The error message and the exit code of a build, then the message of
/// `~CBuildDatabase` (`end_message`).
fn finish(result: Result<(), MakedbError>, end_message: Option<&str>, err: &mut dyn Write) -> i32 {
    let status = match result {
        Ok(()) => 0,
        Err(MakedbError::Build(WriteDbError::CreationBytes(message))) => {
            let _ = err.write_all(b"BLAST Database creation error: ");
            let _ = err.write_all(&message);
            let _ = err.write_all(b"\n");
            1
        }
        Err(error) => {
            let _ = writeln!(err, "{error}");
            error.exit_code()
        }
    };
    if let Some(message) = end_message {
        let _ = writeln!(err, "{message}");
    }
    status
}

#[cfg(test)]
mod tests {
    use super::*;

    // makedb-v5 DECISIONS 12: `[BLAST] MAX_SEQID_LENGTH` is refused with `-parse_seqids`
    // (build_db.cpp:739-747, quoted at `check_registry_settings`), the library's argument
    // rules are makeblastdb's (makeblastdb.cpp:322-329, quoted at `check_arguments`).
    #[test]
    fn registry_and_argument_rules_of_parse_seqids() {
        use clap::Parser;
        #[derive(Parser)]
        struct Wrap {
            #[command(flatten)]
            args: MakedbArgs,
        }
        let parse = |words: &[&str]| {
            let mut argv = vec!["makeblastdb"];
            argv.extend(words);
            Wrap::try_parse_from(argv).unwrap().args
        };
        let plain = parse(&["--dbtype", "prot"]);
        let parsed = parse(&["--dbtype", "prot", "--parse_seqids"]);
        assert_eq!(check_registry_settings(&plain, Ok(Some("100"))), Ok(()));
        assert_eq!(check_registry_settings(&plain, Err("unread")), Ok(()));
        assert_eq!(check_registry_settings(&parsed, Ok(None)), Ok(()));
        let error = check_registry_settings(&parsed, Ok(Some("100"))).unwrap_err();
        assert!(error.contains("[BLAST] MAX_SEQID_LENGTH (100)"), "{error}");
        // The values makeblastdb uses (StringToUInt converts them, not 0) are refused; the
        // ones it ignores are accepted (makeblastdb 2.17.0: pr2/fix/n4probe).
        for used in ["5", "+5", "05", "4294967295"] {
            assert!(
                check_registry_settings(&parsed, Ok(Some(used))).is_err(),
                "{used}"
            );
        }
        for ignored in [
            "",
            "0",
            "+0",
            "abc",
            " 5",
            "5 ",
            "-5",
            "+",
            "4294967296",
            "5x",
        ] {
            assert_eq!(
                check_registry_settings(&parsed, Ok(Some(ignored))),
                Ok(()),
                "{ignored:?}"
            );
        }
        let error =
            check_registry_settings(&parsed, Err("the registry file x sets ...")).unwrap_err();
        assert!(error.starts_with("the registry file x sets ..."), "{error}");
        let mut map = parse(&["--dbtype", "prot"]);
        map.taxid_map = Some("m.txt".to_string());
        let error = check_arguments(&map).unwrap_err();
        assert_eq!(error.exit_code(), 1);
        assert_eq!(
            error.to_string(),
            "Command line argument error: Argument \"parse_seqids\". Must be specified, as it \
             is required by argument:  `taxid_map'"
        );
        map.parse_seqids = true;
        assert!(check_arguments(&map).is_ok());
        map.taxid = Some(1);
        assert_eq!(
            check_arguments(&map).unwrap_err().to_string(),
            "Command line argument error: Argument \"taxid\". Incompatible with argument:  \
             `taxid_map'"
        );
    }

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
