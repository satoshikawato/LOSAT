//! `losat makeblastdb`: the application layer of makeblastdb
//! (`c++/src/app/blastdb/makeblastdb.cpp`, `CMakeBlastDBApp`) for FASTA input and for the
//! sequences of existing databases (`-input_type blastdb`), and a database of format version 5
//! or 4.

use crate::wording::DATABASE_ERROR;
use std::io::{BufRead, Read, Write};

use clap::Args;

use super::build_db::BuildDatabase;
use super::data_size::uint8_to_string_data_size;
use super::format_guess::appears_fasta;
use super::mask_data::process::{process_mask_data, MaskDataOptions};
use super::metadata::{db_files_metadata, metadata_user_path, write_metadata_file, DbMetadata};
use super::platform::Platform;
use super::raw_source::{scan_leaf_tax_ids, RawCopyError, RawSeqDbSource};
use super::taxid_set::{TaxIdMapError, TaxIdSet};
use super::time::LocalTime;
use super::writer::{build_lmdb_file_name, DbVersion, Indices};
use super::WriteDbError;
use crate::blastinput::compressed::{compression_of_name, Compression};
use crate::blastinput::ncbi_environment::{diag_target, ApplicationSettings, DiagStream};
use crate::ncbi_str::string_to_uint8_data_size;
use crate::seqdb::date::iso8601_date_time_exception;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:224-244
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:292-306
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:309-312
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:246-251
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:319-329
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:257-277
// LOSAT's parser checks the range of `-taxid`, the exclusion and the requirements (`-taxid_map`
// needs `-parse_seqids`, `-mask_id` needs `-mask_data`, `-mask_desc` needs `-mask_id`); its errors
// are the parser's (exit status 2, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES) where NCBI's argument
// reader prints its usage (exit status 1). `check_arguments` repeats the two dependencies for
// the library. The options of makeblastdb that LOSAT does not support are not registered;
// `crate::cli` names them in its error.
/// Arguments of `losat makeblastdb`, under makeblastdb's names.
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct MakedbArgs {
    /// Input FASTA file (`-` is standard input), or the input databases of -input_type blastdb
    #[arg(long = "in", default_value = "-", value_name = "input_file")]
    pub input: String,
    /// Type of the data specified in input_file (fasta and blastdb are supported)
    #[arg(long, default_value = "fasta", value_name = "type",
          value_parser = ["fasta", "blastdb", "asn1_bin", "asn1_txt"])]
    pub input_type: String,
    /// Molecule type of target db
    #[arg(long, value_name = "molecule_type", value_parser = ["nucl", "prot"])]
    pub dbtype: String,
    /// Title for BLAST database [default: the input file name given to -in; for -input_type
    /// blastdb, the title of the input database]
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
    /// Path prefix for location of database files in metadata
    #[arg(long, value_name = "String")]
    pub metadata_output_prefix: Option<String>,
    /// File to which the program log should be redirected
    #[arg(long, value_name = "File_Name")]
    pub logfile: Option<String>,
    /// Create index of sequence hash values.
    #[arg(long)]
    pub hash_index: bool,
    /// Taxonomy ID to assign to all sequences
    // NCBI's integer and its constraint (`nonnegative_ncbi_integer`): `0x10` is 16 (NCBI BLAST+
    // 2.17.0, line RG item E1b probe `makeblastdb.k.taxid.12`: rc 0 with `-dryrun`).
    #[arg(long, value_name = "TaxID",
          value_parser = crate::blastinput::value_parsers::nonnegative_ncbi_integer,
          conflicts_with = "taxid_map")]
    pub taxid: Option<i32>,
    /// Text file mapping sequence IDs to taxonomy IDs. Format:<SequenceId> <TaxonomyId><newline>
    #[allow(rustdoc::invalid_html_tags, reason = "NCBI's help, printed by clap")]
    #[arg(long, value_name = "TaxIDMapFile", requires = "parse_seqids")]
    pub taxid_map: Option<String>,
    /// Comma-separated list of input files containing masking data as produced by NCBI masking
    /// applications (e.g. dustmasker, segmasker, windowmasker)
    #[arg(long, value_name = "mask_data_files")]
    pub mask_data: Option<String>,
    /// Comma-separated list of strings to uniquely identify the masking algorithm
    #[arg(long, value_name = "mask_algo_ids", requires = "mask_data")]
    pub mask_id: Option<String>,
    /// Comma-separated list of free form strings to describe the masking algorithm details
    #[arg(long, value_name = "mask_algo_descriptions", requires = "mask_id")]
    pub mask_desc: Option<String>,
    /// Option to parse seqid for FASTA input if set, for all other input types seqids are
    /// parsed automatically
    #[arg(long)]
    pub parse_seqids: bool,
    /// 0x01 Exclude Model
    // NCBI's integer (`AsInteger`, an `int`: `0x01` is 1, `-1` sets every bit), made an `Uint8`
    // (makeblastdb.cpp:1227-1230).
    #[arg(long, value_name = "oid_masks",
          value_parser = crate::blastinput::value_parsers::ncbi_integer)]
    pub oid_masks: Option<i32>,
}

/// The options of makeblastdb 2.17.0 that `losat makeblastdb` refuses
/// (makeblastdb.cpp:246-333). `-h` and `-version` are served (`crate::cli`);
/// `-logfile`, `-parse_seqids`, `-hash_index`, `-taxid`, `-taxid_map`, `-mask_data`, `-mask_id`,
/// `-mask_desc`, `-metadata_output_prefix` and `-oid_masks` are read. `-gi_mask` and
/// `-gi_mask_name` write GI-indexed mask files (item M4c), which LOSAT does not write yet.
pub const UNSUPPORTED_OPTIONS: [&str; 3] = ["randomize", "gi_mask", "gi_mask_name"];

/// How `losat makeblastdb` ends when it does not build a database, or cannot
/// open the one it built.
#[derive(Debug)]
pub enum MakedbError {
    /// An error of makeblastdb, with makeblastdb's exit code.
    Build(WriteDbError),
    /// A database cannot be opened or read (the input databases of `-input_type blastdb`, or
    /// the database written): NCBI's `CSeqDBException`, the text makeblastdb prints after
    /// "BLAST Database error: ".
    Database(SeqDbError),
    /// A request makeblastdb serves and LOSAT does not.
    Unsupported(String),
    /// NCBI's `CIOException::eFlush` of the metadata file's stream: the text makeblastdb
    /// prints after "BLAST failed to write output: ", exit status 6.
    Output(&'static str),
    /// NCBI's `CArgException` raised after the arguments were parsed: the text
    /// makeblastdb prints after "Command line argument error: ".
    Argument(String),
}

impl MakedbError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:149-163
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:167-171
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:199-203
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:212-215
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:256-259
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    /// The exit code: 1 (`BLAST_INPUT_ERROR`) for an error of the options or
    /// of the database creation, 2 (`BLAST_DATABASE_ERROR`) when the new
    /// database cannot be opened, 4 (`BLAST_OUT_OF_MEMORY`) for a failed
    /// allocation, 6 (`BLAST_OUTPUT_ERROR`) when the metadata file's stream cannot be
    /// flushed, and 255 (`BLAST_UNKNOWN_ERROR`) otherwise, which includes LOSAT's own
    /// refusals and an input file that cannot be read.
    pub fn exit_code(&self) -> i32 {
        match self {
            MakedbError::Build(
                WriteDbError::Exception(_) | WriteDbError::Other(_) | WriteDbError::Io { .. },
            ) => 255,
            MakedbError::Build(WriteDbError::OutOfMemory) => 4,
            MakedbError::Output(_) => 6,
            MakedbError::Database(error) if !is_refusal(error) => 2,
            MakedbError::Build(WriteDbError::Database(_)) => 2,
            MakedbError::Database(_)
            | MakedbError::Build(
                WriteDbError::Creation(_)
                | WriteDbError::CreationBytes(_)
                | WriteDbError::Options(_)
                | WriteDbError::Unsupported(_),
            )
            | MakedbError::Unsupported(_)
            | MakedbError::Argument(_) => 1,
        }
    }
}

impl MakedbError {
    /// The message makeblastdb prints for the error (without the line feed that ends it), as
    /// bytes ([`WriteDbError::message`]).
    pub fn message(&self) -> Vec<u8> {
        match self {
            MakedbError::Build(error) => error.message(),
            // LOSAT's refusal of a database where NCBI crashes or reads memory outside its
            // mapping, in LOSAT's words.
            MakedbError::Database(error) if is_refusal(error) => {
                format!("losat makeblastdb: {error}").into_bytes()
            }
            MakedbError::Database(SeqDbError::Ncbi(exception)) => [
                DATABASE_ERROR.as_bytes(),
                b" ",
                exception.what().trim_ascii_end(),
            ]
            .concat(),
            // NCBI prints `e.GetMsg()` of the `CSeqDBException`; the text of an error that is
            // NCBI's own is what `SeqDbError::ncbi_message` gives (the file's path is not part
            // of it; `SeqDbError::NotFound` is bytes), and any other is LOSAT's error in its
            // own words.
            MakedbError::Database(error) => match error.ncbi_message() {
                Some(text) => [DATABASE_ERROR.as_bytes(), b" ", &text].concat(),
                None => format!("{DATABASE_ERROR} {error}").into_bytes(),
            },
            MakedbError::Unsupported(message) => format!("losat makeblastdb: {message}").into(),
            MakedbError::Output(message) => {
                format!("BLAST failed to write output: {message}").into()
            }
            MakedbError::Argument(message) => {
                format!("Command line argument error: {message}").into()
            }
        }
    }
}

impl std::fmt::Display for MakedbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.message()))
    }
}

/// Whether `error` is LOSAT's refusal of a database where NCBI crashes or gives a result that
/// depends on its memory ([`SeqDbError::NcbiCrashes`], [`SeqDbError::Unsupported`]), which is
/// printed as the other refusals of makeblastdb are, with exit status 1.
fn is_refusal(error: &SeqDbError) -> bool {
    matches!(
        error,
        SeqDbError::NcbiCrashes { .. } | SeqDbError::Unsupported { .. }
    )
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

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
/// An error of LOSAT's reader of the input databases (`-input_type blastdb`), as the other
/// programs that read a database report it: NCBI's own (`SeqDbError::ncbi_message`) is a
/// database error ("BLAST Database error: ", exit status 2), another exception of NCBI's
/// (`SeqDbError::ncbi_exception`) is reported as "Error: " and its report (exit status 255), and
/// anything else is LOSAT's refusal.
fn input_error(error: SeqDbError) -> MakedbError {
    match error {
        SeqDbError::Ncbi(exception) => MakedbError::Build(WriteDbError::Exception(exception)),
        error if error.ncbi_message().is_some() => MakedbError::Database(error),
        error => unsupported(format!(
            "{error}; this database is not supported as input by LOSAT's makeblastdb"
        )),
    }
}

/// The error of the copy of the input databases.
fn raw_copy_error(error: RawCopyError) -> MakedbError {
    match error {
        RawCopyError::Read(error) => input_error(error),
        RawCopyError::Build(error) => MakedbError::Build(error),
    }
}

/// Whether `error` is LOSAT's own (a refusal, a failed write of LOSAT's), which removes the
/// files of the build, rather than an error of NCBI's, after which makeblastdb's writer is
/// closed when it is destroyed (`BuildDatabase::close_at_destruction`).
fn is_losat_error(error: &MakedbError) -> bool {
    matches!(
        error,
        MakedbError::Unsupported(_)
            | MakedbError::Build(
                WriteDbError::Unsupported(_) | WriteDbError::Other(_) | WriteDbError::Io { .. }
            )
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1762-1807
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1717-1741
/// `SeqDB_CombineAndQuote`: the names as one database list, separated by spaces, a name with a
/// space between double quotes.
fn combine_and_quote(names: &[&str]) -> String {
    let mut out = String::new();
    for name in names {
        if !out.is_empty() {
            out.push(' ');
        }
        if name.contains(' ') {
            out.push('"');
            out.push_str(name);
            out.push('"');
        } else {
            out.push_str(name);
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1128-1136
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:562-566
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1733-1738,1764-1766
// makeblastdb reads two variables of its own. LOSAT reproduces makeblastdb
// with both unset (or empty, which `Get` cannot tell from unset). Neither is read for the
// sequences of a database (`-input_type blastdb`, which `AddSequences` copies as they are
// stored), so neither is refused for that input.
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
/// of a FASTA input (`args`). `env` gives the value of a variable.
pub fn check_environment(
    args: &MakedbArgs,
    env: &dyn Fn(&str) -> Option<std::ffi::OsString>,
) -> Result<(), String> {
    const DONT_SCAN: &str = "NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS";
    if args.input_type != "fasta" {
        return Ok(());
    }
    if env(DONT_SCAN).is_some_and(|value| !value.is_empty()) {
        return Err(format!(
            "the environment variable {DONT_SCAN}, with which NCBI makeblastdb stores other \
             titles for some FASTA definition lines, is not supported by LOSAT's makeblastdb"
        ));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:739-747
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1211-1216
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:485-491,647-653
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:430-434
// The registry's limit of a Seq-id's length changes which FASTA ids makeblastdb parses with
// `-parse_seqids` when `NStr::StringToUInt` converts the value and the value is not 0 (0 is
// the reader's "not set"; the static `s_Max*` copies are read nowhere). A value it cannot
// convert is ignored (`catch`), as is an empty one. LOSAT's FASTA reader keeps NCBI's default
// limits (makedb-v5 DECISIONS 12), so it refuses the entry with `-parse_seqids` exactly when
// makeblastdb would use it. `[BLAST] LONG_SEQID` is refused for every program by the
// application settings (`ApplicationStart::settings`; it is not one of the entries that change
// nothing). The FASTA reader reads only a FASTA input: the entry changes nothing for
// `-input_type blastdb`, whose Seq-ids are copied as they are stored.
/// Refuses `[BLAST] MAX_SEQID_LENGTH` (`max_seqid_length`, from
/// `ApplicationSettings::max_seqid_length`: the value of the first layer that has the entry,
/// `Err` for a file value LOSAT cannot read) for a build with `-parse_seqids` when makeblastdb
/// would use the value.
pub fn check_registry_settings(
    args: &MakedbArgs,
    max_seqid_length: Result<Option<&str>, &str>,
) -> Result<(), String> {
    if !args.parse_seqids || args.input_type != "fasta" {
        return Ok(());
    }
    match max_seqid_length {
        Err(reason) => Err(format!(
            "{reason} ([BLAST] MAX_SEQID_LENGTH changes the Seq-ids NCBI makeblastdb parses \
             with -parse_seqids)"
        )),
        // `NStr::StringToUInt(value)` with no flags: no limit where it throws.
        Ok(Some(value))
            if crate::ncbi_str::string_to_uint(value.as_bytes()).is_ok_and(|limit| limit != 0) =>
        {
            Err(format!(
                "the registry entry [BLAST] MAX_SEQID_LENGTH ({value}), which changes the Seq-ids \
                 NCBI makeblastdb parses with -parse_seqids, is not supported by LOSAT's makeblastdb"
            ))
        }
        Ok(_) => Ok(()),
    }
}

/// The checks of the arguments that LOSAT makes first: what LOSAT does not
/// support, and the format version.
fn check_arguments(args: &MakedbArgs) -> Result<DbVersion, MakedbError> {
    if args.input_type != "fasta" && args.input_type != "blastdb" {
        return Err(unsupported(format!(
            "-input_type {} is not supported; only fasta and blastdb input are",
            args.input_type
        )));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1218-1219
    // The parser allows 4 and 5 only (`CArgAllow_Integers(eBDB_Version4, eBDB_Version5)`).
    let Some(version) = DbVersion::from_number(args.blastdb_version) else {
        return Err(unsupported(format!(
            "-blastdb_version {} is not a version makeblastdb writes (4 or 5)",
            args.blastdb_version
        )));
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:322-329
    // The parser checks these on the command line (`-taxid` has its range there); a caller of
    // the library that sets the fields is checked here with the parser's sentences.
    // The same for the masking options: `-mask_id` requires `-mask_data`, `-mask_desc` requires
    // `-mask_id` (makeblastdb.cpp:276-277, cited at `MakedbArgs`).
    for (name, given, needed, needed_given) in [
        (
            "mask_id",
            args.mask_id.is_some(),
            "mask_data",
            args.mask_data.is_some(),
        ),
        (
            "mask_desc",
            args.mask_desc.is_some(),
            "mask_id",
            args.mask_id.is_some(),
        ),
    ] {
        if given && !needed_given {
            return Err(MakedbError::Argument(format!(
                "Argument \"{needed}\". Must be specified, as it is required by argument:  \
                 `{name}'"
            )));
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1005-1014
    // `-taxid_map` requires `-parse_seqids` on the command line for every input type, although
    // an input other than FASTA always has its Seq-ids parsed.
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
// `AsInputFile` opens the file here: one that cannot be opened is the parser's error
// "File is not accessible" (makeblastdb 2.17.0: "Command line argument error: Argument
// "taxid_map". File is not accessible:  `<name>'", exit status 1, after the head of the
// log). A number the map cannot convert ends makeblastdb with its exception report
// (`CATCH_ALL`, "Error: ", exit status 255; the build tree in its path is the official
// release's, makedb-v5 DECISIONS 29).
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-737
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-621,95-99
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
            TaxIdMapError::Exception(e) => MakedbError::Build(WriteDbError::Exception(e)),
            TaxIdMapError::BlankLine(line) => unsupported(format!(
                "line {line} of the -taxid_map file {name} holds only white space, which \
                 makeblastdb 2.17.0 cannot read (it crashes)"
            )),
        })?;
    Ok((Some(taxids), (name == "-").then_some(failed)))
}

/// The value of `-max_file_sz`, read where makeblastdb reads it: after the constructor of
/// `CBuildDatabase` (the checks of the output directory, then the head of the log and the
/// removal of an existing database of the same name, build_db.cpp:1165-1177), before the log
/// line "Maximum file size".
fn max_file_size(args: &MakedbArgs) -> Result<u64, MakedbError> {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1267-1272
    // A value NCBI cannot convert throws a `CStringException`, which `CATCH_ALL` reports
    // through its `CException` clause (`MakedbError::exit_code`): "Error: " and the
    // exception's report, exit status 255.
    let bytes = string_to_uint8_data_size(&args.max_file_sz)
        .map_err(|error| MakedbError::Build(WriteDbError::Exception(error.exception())))?;
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

/// The input of `-in x.gz`, `x.bz2` or `x.zst`, inflated while it is read
/// (`compressed::Compression::decoder`: every member, stream or frame, as pbzip2 and pzstd write
/// several), as makeblastdb's `CDecompressIStream` reads it (which reads every gzip member but
/// only the first bzip2 stream or zstd frame: `PD-LOSAT-NCBI-DEFECTS`). A fault of the compressed
/// data (not that format, cut short) is kept apart from a failed read of the file itself, so that
/// the build can name it (`fault`).
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050
struct CompressedInput<'a> {
    inflater: Box<dyn Read + 'a>,
    /// Set by `ReadWatch` when a read of the file itself failed.
    file_failed: std::rc::Rc<std::cell::Cell<bool>>,
    fault: Option<String>,
}

/// A reader that remembers whether a read of its own failed.
struct ReadWatch<R: Read> {
    inner: R,
    failed: std::rc::Rc<std::cell::Cell<bool>>,
}

impl<R: Read> Read for ReadWatch<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let result = self.inner.read(buf);
        if let Err(error) = &result {
            if error.kind() != std::io::ErrorKind::Interrupted {
                self.failed.set(true);
            }
        }
        result
    }
}

impl<'a> CompressedInput<'a> {
    fn new(input: impl Read + 'a, compression: Compression) -> Self {
        let file_failed = std::rc::Rc::new(std::cell::Cell::new(false));
        let watch = ReadWatch {
            inner: input,
            failed: std::rc::Rc::clone(&file_failed),
        };
        Self {
            inflater: compression.decoder(watch),
            file_failed,
            fault: None,
        }
    }

    /// The fault of the compressed data that ended the reading, if the file itself read without
    /// error.
    fn fault(&self) -> Option<&str> {
        if self.file_failed.get() {
            return None;
        }
        self.fault.as_deref()
    }
}

impl Read for CompressedInput<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let result = self.inflater.read(buf);
        if let Err(error) = &result {
            if error.kind() != std::io::ErrorKind::Interrupted && self.fault.is_none() {
                self.fault = Some(error.to_string());
            }
        }
        result
    }
}

/// The sentence of a `-in` name ending in `.gz`, `.bz2` or `.zst` whose bytes are not that
/// format, or are cut short (line GZ, DECISIONS 5; for `.gz`: where makeblastdb's decompressor
/// fails, `compressed::GzipMembers`): exit status 1, the files of the build removed.
fn compressed_fault_error(name: &str, compression: Compression, fault: &str) -> MakedbError {
    unsupported(format!("-in: {}", compression.not_data(name, fault)))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:792-836
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
        // A compressed name is guessed on its inflated bytes, as `x_VerifyInputFilesType` guesses
        // on its `CDecompressIStream` (makeblastdb.cpp:815-830); a fault of the compressed data
        // reads as "not FASTA".
        let looks_like_fasta = std::fs::File::open(seq_file).is_ok_and(|file| {
            match compression_of_name(std::path::Path::new(seq_file)) {
                Some(compression) => appears_fasta(CompressedInput::new(file, compression)),
                None => appears_fasta(file),
            }
        });
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
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1273-1276
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1294-1296
// The title defaults to the whole value of `-in`, several names included. After a
// successful build makeblastdb opens the new database with its reader and writes
// the JSON metadata file (`open_new_database`); an error of that open ends
// makeblastdb with "BLAST Database error" and exit status 2 although the database
// was written.
//
// For `-input_type blastdb` the output must be named apart from `-in` (the same string is an
// error even when `-out` is given), the title defaults to the title of the one input database,
// which is opened for it before the log is opened (several names need `-title`), and the
// indices are always full (`x_ShouldParseSeqIds`: every input but FASTA has its Seq-ids
// parsed).
/// `CMakeBlastDBApp::x_BuildDatabase`. The log goes to `out` (makeblastdb's
/// standard output, or the file of `-logfile`), warnings and diagnostics to
/// `err`; `stdin` is read when `-in` is `-`. `now` is the local time for the
/// database's date (`None` where the platform has no clock), `platform` the
/// platform whose makeblastdb the files follow ([`Platform::host`] for a build).
///
/// `logfile_error` is the error of a `-logfile` that could not be opened, raised where
/// makeblastdb opens the file. `end_messages` gets the messages of `~CBuildDatabase` (no
/// sequence matched a mask, `-taxid_map` matched no sequence), which makeblastdb posts after
/// everything else, the error of the build included.
#[allow(clippy::too_many_arguments)]
fn build_database_with_log(
    args: &MakedbArgs,
    registry_blastdb: Option<&str>,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
    logfile_error: Option<MakedbError>,
    end_messages: &mut Vec<Vec<u8>>,
) -> Result<(), MakedbError> {
    let version = check_arguments(args)?;
    let is_protein = args.dbtype == "prot";
    let blastdb_input = args.input_type == "blastdb";

    let dbname = args.out.as_deref().unwrap_or(&args.input);
    if blastdb_input && dbname == args.input {
        return Err(options_error(
            "Cannot create a BLAST database from an existing one without changing the output \
             name, please provide a (different) database name using -out",
        ));
    }
    // NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-275,282-301
    // `NStr::Split(dbname, " ", input_files)` runs with no flag
    // (ncbistr.hpp:2551-2555): adjacent delimiters are not merged and none is
    // truncated, so every space ends a token, an empty one included ("db "
    // is "db" and "", " db" is "" and "db"). A name with a space anywhere
    // gives more than one token; the empty name gives none.
    if dbname == "-" || dbname.contains(' ') {
        return Err(options_error("Please provide a database name using -out"));
    }
    let title = match args.title.as_deref() {
        Some(title) => title.as_bytes().to_vec(),
        None if blastdb_input => input_database_title(&args.input, is_protein, registry_blastdb)?,
        None => args.input.as_bytes().to_vec(),
    };
    if title == b"-" {
        return Err(options_error("Please provide a title using -title"));
    }
    if let Some(error) = logfile_error {
        return Err(error);
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1016-1028
    // A first name `-` reads the standard input and nothing else. With no name at
    // all (`-in ""`) NCBI reads `names[0]` of an empty vector (makeblastdb 2.17.0: a
    // segmentation fault after the log line "Maximum file size"); LOSAT refuses that where
    // NCBI reads it, after `-max_file_sz` and `-taxid_map`, whose errors come first.
    let names = split_quoted(&args.input);
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050
    // A name ending in `.gz`, `.bz2` or `.zst` is read through its decoder (`CompressedInput`).

    let Some(now) = now else {
        return Err(unsupported(
            "this platform has no clock for the database's date",
        ));
    };
    let mut db = BuildDatabase::new(
        dbname,
        &title,
        is_protein,
        out,
        now,
        version,
        platform,
        Indices {
            full: args.parse_seqids || blastdb_input,
            hash: args.hash_index,
            trace: false,
        },
    )?;
    if let Some(oid_masks) = args.oid_masks {
        db.set_oid_masks(oid_masks as u64);
    }
    let max_file_size = max_file_size(args)?;
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
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1288-1294
    // The mask files are read before the input is looked at. A failure ends the run at once
    // (`CATCH_ALL`), and the destructor's messages follow it (`end_messages`).
    let mask_options = MaskDataOptions {
        mask_data: args.mask_data.as_deref(),
        mask_id: args.mask_id.as_deref(),
        mask_desc: args.mask_desc.as_deref(),
        root: None,
    };
    let masked = {
        let mut warn = warn_to(err);
        process_mask_data(&mut db, &mask_options, &mut warn)
    };
    if let Err(error) = masked {
        let error = error.into_makedb_error();
        *end_messages = destructor_messages(&db, error.as_ref());
        return error.map_or(Ok(()), Err);
    }
    let user_path = metadata_user_path(args.metadata_output_prefix.as_deref());
    if blastdb_input {
        let result = copy_input_databases(
            &mut db,
            &names,
            &user_path,
            is_protein,
            registry_blastdb,
            err,
        );
        *end_messages = destructor_messages(&db, result.as_ref().err());
        if result.as_ref().is_err_and(|error| !is_losat_error(error)) {
            end_messages.extend(db.close_at_destruction());
        }
        return result;
    }
    if names.is_empty() {
        return Err(unsupported(
            "-in names no input file; losat makeblastdb reads FASTA files or the standard input",
        ));
    }
    let from_stdin = names[0] == "-";
    let result = add_input_and_close(
        &mut db,
        &names,
        from_stdin.then_some(stdin_map),
        &user_path,
        is_protein,
        registry_blastdb,
        stdin,
        err,
    );
    *end_messages = destructor_messages(&db, result.as_ref().err());
    result
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1175-1188
// The names are split as the reader splits a list; with no name at all NCBI opens the first of
// none, which LOSAT refuses.
/// The title of `-input_type blastdb` without `-title`: the title of the one database `input`
/// names (`GetTitle`).
fn input_database_title(
    input: &str,
    is_protein: bool,
    registry_blastdb: Option<&str>,
) -> Result<Vec<u8>, MakedbError> {
    let names = split_quoted(input);
    if names.len() > 1 {
        return Err(options_error("Please provide a title using -title"));
    }
    let Some(&name) = names.first() else {
        return Err(unsupported(
            "-in names no database and -title is not given; NCBI makeblastdb then reads the \
             title of the first database of none (no defined result)",
        ));
    };
    let search = SearchPath::from_environment(registry_blastdb).map_err(input_error)?;
    let db = SeqDb::open_in(name, seq_type(is_protein), &search).map_err(input_error)?;
    Ok(db.title().to_vec())
}

fn seq_type(is_protein: bool) -> SeqType {
    if is_protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1077-1124
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:792-798
// The names are not checked as files. Each is opened alone first: one NCBI's reader cannot
// open is posted and kept, and the open of all of them together then fails with the reader's
// error. A name LOSAT's reader refuses, or one that raises another exception of NCBI's (which
// the check does not catch), ends the build there. The leaf taxids of all the input's
// definition lines are read before the copy; the count of sequences is not checked (a copy of
// none ends with "No volumes were created.").
//
// The copy is not inside a `try` that ends the build: an error of NCBI's leaves the output
// open, and `build_database_with_log` closes it after the messages, as makeblastdb's
// destructors do. LOSAT's refusals and its own errors end the build and remove its files.
/// `x_ProcessInputData` for `-input_type blastdb`, then `EndBuild` and the metadata file
/// (`open_new_database`).
#[allow(clippy::too_many_arguments)]
fn copy_input_databases(
    db: &mut BuildDatabase<'_>,
    names: &[&str],
    user_path: &str,
    is_protein: bool,
    registry_blastdb: Option<&str>,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let seqtype = seq_type(is_protein);
    let search = SearchPath::from_environment(registry_blastdb).map_err(input_error)?;
    for &name in names {
        match SeqDb::open_in(name, seqtype, &search) {
            Ok(_) => {}
            Err(error) if error.ncbi_message().is_some() => {
                let _ = writeln!(
                    err,
                    "Error: [makeblastdb] Unable to open input {name} as BLAST db"
                );
            }
            Err(error) => return Err(input_error(error)),
        }
    }
    if names.is_empty() {
        return Err(options_error(
            "No valid input FASTA file or BLAST db is found.",
        ));
    }
    let input = SeqDb::open_in(&combine_and_quote(names), seqtype, &search).map_err(input_error)?;
    db.set_leaf_tax_ids(scan_leaf_tax_ids(&input).map_err(raw_copy_error)?);
    let mut source = RawSeqDbSource::new(&input, db).map_err(raw_copy_error)?;
    let copied = {
        let mut warn = warn_to(err);
        db.add_raw_sequences(&mut source, &mut warn)
    };
    if let Err(error) = copied {
        let error = raw_copy_error(error);
        if is_losat_error(&error) {
            let mut warn = warn_to(err);
            let _ = db.end_build(true, &mut warn);
        }
        return Err(error);
    }
    let closed = {
        let mut warn = warn_to(err);
        db.end_build(false, &mut warn)
    };
    if closed? {
        open_new_database(db, user_path, is_protein, registry_blastdb, err)?;
    }
    Ok(())
}

/// What `~CBuildDatabase` posts (build_db.cpp:1204-1213), in its order: the message about the
/// masks, then the one about the taxids. Each entry is one message, without its final newline.
/// Nothing after a refusal of LOSAT (`error`): those messages are NCBI's words about the state
/// of NCBI's run, which LOSAT did not follow (the refusal of `-in ""` posts none either).
fn destructor_messages(db: &BuildDatabase<'_>, error: Option<&MakedbError>) -> Vec<Vec<u8>> {
    if matches!(
        error,
        Some(MakedbError::Unsupported(_) | MakedbError::Build(WriteDbError::Unsupported(_)))
    ) || matches!(error, Some(MakedbError::Database(error)) if is_refusal(error))
    {
        return Vec::new();
    }
    let mut messages = Vec::new();
    db.post_destruction_diagnostics(&mut |line: &[u8]| {
        messages.push(line.to_vec());
    });
    messages.extend(
        db.taxid_message()
            .map(|message| message.as_bytes().to_vec()),
    );
    messages
}

/// The rest of `x_BuildDatabase`: the input files, `EndBuild` and the metadata file.
/// `from_stdin` is `Some(map)` for `-in -`, `map` being `Some(failed)` when `-taxid_map -`
/// read the standard input first (`failed`: whether it left the stream failed). `user_path`
/// comes before the names of the metadata's files (`metadata_user_path`).
#[allow(clippy::too_many_arguments)]
fn add_input_and_close(
    db: &mut BuildDatabase<'_>,
    names: &[&str],
    from_stdin: Option<Option<bool>>,
    user_path: &str,
    is_protein: bool,
    registry_blastdb: Option<&str>,
    stdin: &mut dyn Read,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    if let Some(stdin_map) = from_stdin {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1398-1420
        // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:619-636
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
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
        // Each file goes to its own `AddFasta` (`x_AddSequenceData` ->
        // `x_AddFasta`): its own FASTA reader (line numbers and local ids start
        // again), its own "Adding sequences" line, and "No sequences added" for a
        // file without a sequence.
        for &name in names {
            // A file that exists and cannot be opened (`x_VerifyInputFilesType` has warned that
            // it does "NOT appear to be FASTA"). makeblastdb does not test the stream it
            // opened: `AddFasta` skips a stream that is not good and adds nothing from it, not
            // even its "Adding sequences" line; a build of such files alone ends with "No
            // volumes were created." and exit status 0 (makeblastdb 2.17.0, a compressed name
            // too). A deterministic NCBI result, reproduced.
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1398-1408
            let Ok(mut file) = std::fs::File::open(name) else {
                continue;
            };
            // A compressed name is inflated while it is read (`CompressedInput`); a fault of
            // its data is reported by name, after the build removed its files as for any
            // failed input.
            match compression_of_name(std::path::Path::new(name)) {
                Some(compression) => {
                    let mut input = CompressedInput::new(file, compression);
                    let result = add_fasta(db, &mut input, err);
                    if let Some(fault) = input.fault() {
                        return Err(compressed_fault_error(name, compression, fault));
                    }
                    result?
                }
                None => add_fasta(db, &mut file, err)?,
            }
        }
    }

    let closed = {
        let mut warn = warn_to(err);
        db.end_build(false, &mut warn)
    };
    if closed? {
        open_new_database(db, user_path, is_protein, registry_blastdb, err)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1296-1328
// `CSeqDB sdb(new_db, t)` throws a `CSeqDBException` when the reader cannot
// open what was written; the usual causes are an alias file of an earlier
// database of the same name that `DeleteBlastDb` could not resolve and
// therefore left in place, which now hides the new single volume, and an
// index file or alias file that is not a file (one that could not be created
// or renamed, a link to a device: `CFile::Exists` asks for a regular file), or
// an index file that a failed write left short (below).
// The other files of a volume are not opened then: one that is missing is
// left out of the metadata (`db_files_metadata`). LOSAT
// resolves the database here for that outcome and for the values of the
// metadata: the alias files are read, every volume is located and the header of
// its index file is read as NCBI's reader reads it (`SeqDb::locate_one_name_metadata`),
// but no volume is opened. The volumes were written a moment ago, and
// LOSAT's reader holds the index of an open volume in memory (12 bytes for each
// sequence). So the error `CSeqDB` gives for more volumes than the process may
// hold open files ("Too many open files, please raise the open file limit",
// exit status 2 after a complete build) is not reproduced. The totals, the
// title, the date and the version of the metadata are those the headers hold,
// which are the writer's unless an index file's write failed: NCBI's reader
// takes a file that ends early to be zeros up to the end of its last page
// (a missing string is empty, a missing count is 0, a missing version is not 4
// or 5 and ends the build with "Not a valid version 4 or 5 database."), and a
// file of no bytes is mapped as a null pointer, and NCBI's process is killed
// by SIGSEGV (no metadata file). LOSAT refuses it ("NCBI crashes on it"; the
// files of the build stay as NCBI leaves them). A read past the last page
// reads memory that is not the file's: what NCBI does there depends on its
// memory (it may go on, or fault), so LOSAT refuses that too, in other words.
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
    user_path: &str,
    is_protein: bool,
    registry_blastdb: Option<&str>,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let new_db = db.output_dbname();
    let seqtype = seq_type(is_protein);
    // `CSeqDB` names its search path in the message of a database that is not there: the
    // current directory, `BLASTDB` and the registry entry `[BLAST] BLASTDB` (`registry_blastdb`,
    // of the application settings that the callers of `run_with` read before the build; a value
    // that LOSAT does not read is left out).
    let search = SearchPath::from_environment(registry_blastdb).map_err(MakedbError::Database)?;
    let located = match SeqDb::locate_one_name_metadata(new_db, seqtype, &search) {
        Ok(located) => located,
        Err(error) => return Err(MakedbError::Database(error)),
    };
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1649-1656
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:1638-1644,1662-1668
    // The totals are those of the volumes `CSeqDB` opened: the ones the alias file names (a
    // cut alias file names fewer than the writer made), or the one volume. Each volume's
    // count and length are the ones its index file holds in its header
    // (`LocatedDatabase::volume_totals`); no volume is opened.
    let totals = &located.volume_totals;
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1705-1747
    // Without `-taxid` and `-taxid_map` every OID has tax id 0, so
    // `number-of-taxids` is not set. `GetDBTaxIds` reads the taxids of the tax id files,
    // which are the ones the writer gave them (`WriteDb::db_tax_ids`).
    let dbname = match new_db.rfind(std::path::MAIN_SEPARATOR) {
        Some(off) => &new_db[off + 1..],
        None => new_db,
    };
    let version5 = located.version5;
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1730-1732
    // The date is the one the volumes' indexes hold; `CTime` throws where a date that a failed
    // write cut short does not read, and makeblastdb ends with the exception (exit status
    // 255), after the files of the build were written.
    let last_updated = match &located.date {
        Ok(date) => iso8601_date_time_exception(date),
        Err(error) => Err((*error.exception).clone()),
    }
    .map_err(|exception| MakedbError::Build(WriteDbError::Exception(exception)))?;
    let lmdb: Vec<std::path::PathBuf> = version5
        .then(|| build_lmdb_file_name(new_db, is_protein).into())
        .into_iter()
        .collect();
    let (bytes_total, bytes_to_cache, files) = db_files_metadata(
        &located.volumes,
        &located.alias_files,
        &lmdb,
        is_protein,
        user_path,
    );
    let metadata = DbMetadata {
        dbname: dbname.as_bytes().to_vec(),
        dbtype: if is_protein { "Protein" } else { "Nucleotide" },
        db_version: if version5 { 5 } else { 4 },
        description: located.title,
        number_of_letters: totals.iter().map(|(_, letters)| letters).sum(),
        // NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.def:11-12
        // NCBI's `GetDBMetaData` counts in an `int num_seqs` (c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1707), which a
        // database of more than 2^31 - 1 sequences would overflow; LOSAT sums in `Int8`.
        number_of_sequences: totals.iter().map(|(oids, _)| i64::from(*oids)).sum(),
        last_updated,
        number_of_volumes: located.volumes.len() as i32,
        number_of_taxids: version5
            .then(|| db.output_db().db_tax_ids())
            .filter(|tax_ids| {
                tax_ids.len() > 1 || (tax_ids.len() == 1 && tax_ids.first() != Some(&0))
            })
            .map(|tax_ids| tax_ids.len() as i32),
        bytes_total,
        bytes_to_cache,
        files,
    };
    match write_metadata_file(new_db, is_protein, &metadata, err) {
        Some(message) => Err(MakedbError::Output(message)),
        None => Ok(()),
    }
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
// Every line of the log ends with `endl`, which flushes the stream, so the
// log and the messages of the standard error stream reach a common
// destination (`2>&1`) in the order of the events. The standard output of
// the Rust library is line buffered whatever it is connected to; the log
// goes to it without another buffer.
/// Runs `losat makeblastdb` with the settings of its application start and returns the exit code.
pub fn run(args: &MakedbArgs, settings: &ApplicationSettings) -> i32 {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    let now = LocalTime::now();
    let status = run_with(
        args,
        registry_blastdb(settings),
        now.as_ref(),
        Platform::host(),
        &mut stdin.lock(),
        &mut out,
        &mut err,
    );
    let _ = out.flush();
    status
}

/// The registry entry `[BLAST] BLASTDB` of `settings` that `run_with` takes: `None` when the
/// registry has no such entry or LOSAT does not read its value (`ApplicationSettings::blastdb`).
pub(crate) fn registry_blastdb(settings: &ApplicationSettings) -> Option<&str> {
    settings.blastdb.as_ref().ok().and_then(Option::as_deref)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:843-871
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-620
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7593-7601 (SetLogFile)
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:776-792 (CArg_OutputFile::x_Open)
// With `-logfile F` makeblastdb sends its diagnostics (warnings, errors and the
// final error message) to F from the start, and its log to F as well
// (`AsOutputFile` with `fAppend`, makeblastdb.cpp:1199-1201): both append, in the
// order they are written, and nothing goes to the standard streams (two runs
// append to one file). A file that cannot be opened leaves the diagnostics on the
// standard error stream, and the build stops where makeblastdb opens the log,
// after the checks of the database name and the title, with "Command line
// argument error: Argument "logfile". File is not accessible:  `F'" and exit
// status 1. The diagnostics' file is the one that NCBI's diagnostics open for the
// name (`ncbi_environment::diag_target`) and the log's the one that the system opens
// under the name as written (`log_of`): one file for most names, two where NCBI's
// diagnostics read the name otherwise (` mk.log` is `mk.log` there, `a/../mk.log`
// the `mk.log` beside the link `a`), and no log for `logs/`, which the
// diagnostics read as the file `logs` (NCBI BLAST+ 2.17.0; line WQA's D11 round 2,
// 2026-11-20). Three names are special. `-` keeps the diagnostics on the standard
// error stream and opens the standard output for the log: makeblastdb runs as
// without `-logfile`. `/dev/null` (that name, on every platform) and the empty
// name discard the diagnostics; the log then goes to the file `/dev/null`, and
// the empty name is no file, so its build stops with the error above, which is a
// diagnostic and so is not shown (makeblastdb 2.17.0: `-logfile ""` prints
// nothing, exit status 1). LOSAT reads `-logfile` as an option of makeblastdb
// only, not anywhere in the command line as NCBI's application does (for example
// as the value of `-title`).
/// One run with the registry entry `[BLAST] BLASTDB` of the application settings
/// (`registry_blastdb`) and the given time, platform and streams: the build, then the error
/// message and the exit code of a failed build. With `-logfile` the log and the
/// messages go where makeblastdb sends them instead of `out` and `err`.
pub(crate) fn run_with(
    args: &MakedbArgs,
    registry_blastdb: Option<&str>,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    let mut end_messages = Vec::new();
    let name = match args.logfile.as_deref() {
        None | Some("-") => {
            let result = build_database_with_log(
                args,
                registry_blastdb,
                now,
                platform,
                stdin,
                out,
                err,
                None,
                &mut end_messages,
            );
            return finish(result, &end_messages, err);
        }
        Some(name) => name,
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1199-1201
    // The diagnostics go to the target of `-logfile` (`ncbi_environment::diag_target`, which
    // owns the special names), and the log to `AsOutputFile` of the same name (`log_of`).
    let target = diag_target(Some(std::ffi::OsStr::new(name)));
    let log = log_of(name);
    let mut diag = DiagStream {
        target: &target,
        stderr: err,
    };
    match log {
        Some(mut log) => {
            let result = build_database_with_log(
                args,
                registry_blastdb,
                now,
                platform,
                stdin,
                &mut log,
                &mut diag,
                None,
                &mut end_messages,
            );
            finish(result, &end_messages, &mut diag)
        }
        None => {
            let error = MakedbError::Argument(format!(
                "Argument \"logfile\". File is not accessible:  `{name}'"
            ));
            let result = build_database_with_log(
                args,
                registry_blastdb,
                now,
                platform,
                stdin,
                out,
                &mut diag,
                Some(error),
                &mut end_messages,
            );
            finish(result, &end_messages, &mut diag)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1199-1201
/// The log of `-logfile name` (`AsOutputFile`, append): the file that the system opens under
/// the name as written; none for the empty name and a name that it cannot open (`logs/`). The
/// diagnostics open the file of NCBI's reading of the name (`diag_target`), which is this file
/// for most names and another one where that reading differs (` mk.log`, `a/../mk.log` beside
/// a link `a`; NCBI BLAST+ 2.17.0 makeblastdb writes its log to ` mk.log` and its
/// diagnostics to `mk.log`; line WQA's D11 round 2, 2026-11-20). Both append, so where they are
/// one file the two streams land in it in the order they are written, as NCBI's two appending
/// streams on the file do. (A second open, not a copy of the diagnostics' handle, `try_clone`,
/// which the WASI targets cannot make: the run stopped with "File is not accessible".)
fn log_of(name: &str) -> Option<std::fs::File> {
    if name.is_empty() {
        return None;
    }
    std::fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(name)
        .ok()
}

/// The error message and the exit code of a build, then the messages of
/// `~CBuildDatabase` (`end_messages`).
fn finish(result: Result<(), MakedbError>, end_messages: &[Vec<u8>], err: &mut dyn Write) -> i32 {
    let status = match result {
        Ok(()) => 0,
        Err(error) => {
            let _ = err.write_all(&[error.message().as_slice(), b"\n"].concat());
            error.exit_code()
        }
    };
    for message in end_messages {
        let _ = err.write_all(&[message.as_slice(), b"\n"].concat());
    }
    status
}

#[cfg(test)]
#[path = "makeblastdb_tests.rs"]
mod tests;
