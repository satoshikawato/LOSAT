//! `losat makeblastdb`: the application layer of makeblastdb
//! (`c++/src/app/blastdb/makeblastdb.cpp`, `CMakeBlastDBApp`) for FASTA input and for the
//! sequences of existing databases (`-input_type blastdb`), and a database of format version 5
//! or 4.

use crate::wording::DATABASE_ERROR;
use std::ffi::OsString;
use std::io::{BufRead, Read, Write};

use super::asn_source::AsnFormat;
use super::build_db::{endl, BuildDatabase};
use super::data_size::uint8_to_string_data_size;
use super::format_guess::appears_fasta;
use super::mask_data::process::{process_mask_data, MaskDataOptions};
use super::metadata::{db_files_metadata, metadata_user_path, write_metadata_file, DbMetadata};
use super::ofstream::StreamBuf;
use super::platform::Platform;
use super::raw_source::{scan_leaf_tax_ids, RawCopyError, RawSeqDbSource};
use super::taxid_set::{TaxIdMapError, TaxIdSet};
use super::time::LocalTime;
use super::writer::{build_lmdb_file_name, DbVersion, Indices};
use super::WriteDbError;
use crate::blastinput::compressed::{compression_of_name, Compression};
use crate::blastinput::ncbi_environment::{
    diag_target, ApplicationSettings, DiagStream, ReadVerification,
};
use crate::ncbi_str::string_to_uint8_data_size_bytes;
use crate::seqdb::date::iso8601_date_time_exception;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:224-244
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:292-306
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:309-312
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:246-251
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:319-329
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:257-288
// NCBI's argument reader (`cli::reader`) checks the constraints, the exclusions and the
// requirements (`-taxid_map` needs `-parse_seqids`, `-mask_id` needs `-mask_data`, `-mask_desc`
// needs `-mask_id`, `-gi_mask` needs `-parse_seqids` and excludes `-mask_id`, `-gi_mask_name`
// needs `-mask_data` and `-gi_mask`) of the command line and of a library call, in NCBI's order
// and words (its usage, exit status 1); `from_args` then reads the arguments it read. The options
// of makeblastdb that LOSAT does not support are refused by name after the reader
// (`ncbi_app::args::refusal`).
/// Arguments of `losat makeblastdb`, under makeblastdb's names.
#[derive(Debug, Clone)]
pub struct MakedbArgs {
    /// Input FASTA or ASN.1 file (`-` is standard input), or the input databases of -input_type
    /// blastdb
    pub input: OsString,
    /// Type of the data specified in input_file (fasta, blastdb, asn1_txt and asn1_bin are
    /// supported)
    pub input_type: String,
    /// Molecule type of target db
    pub dbtype: String,
    /// Title for BLAST database [default: the input file name given to -in; for -input_type
    /// blastdb, the title of the input database]
    pub title: Option<OsString>,
    /// Name of BLAST database to be created [default: the input file name given to -in]
    pub out: Option<OsString>,
    /// Version of BLAST database to be created
    pub blastdb_version: i32,
    /// Maximum file size for BLAST database files
    pub max_file_sz: OsString,
    /// Path prefix for location of database files in metadata
    pub metadata_output_prefix: Option<OsString>,
    /// File to which the program log should be redirected
    pub logfile: Option<OsString>,
    /// Create index of sequence hash values.
    pub hash_index: bool,
    /// Taxonomy ID to assign to all sequences
    // NCBI's integer and its constraint, at least 0 (NCBI's reader, `ncbi_app::constraint`): `0x10`
    // is 16 (NCBI BLAST+ 2.17.0, line RG item E1b probe `makeblastdb.k.taxid.12`: rc 0 with
    // `-dryrun`).
    pub taxid: Option<i32>,
    /// Text file mapping sequence IDs to taxonomy IDs. Format: `<SequenceId> <TaxonomyId><newline>`
    pub taxid_map: Option<OsString>,
    /// Comma-separated list of input files containing masking data as produced by NCBI masking
    /// applications (e.g. dustmasker, segmasker, windowmasker)
    pub mask_data: Option<OsString>,
    /// Comma-separated list of strings to uniquely identify the masking algorithm
    pub mask_id: Option<OsString>,
    /// Comma-separated list of free form strings to describe the masking algorithm details
    pub mask_desc: Option<OsString>,
    /// Create GI indexed masking data.
    pub gi_mask: bool,
    /// Comma-separated list of masking data output files.
    pub gi_mask_name: Option<OsString>,
    /// Option to parse seqid for FASTA input if set, for all other input types seqids are
    /// parsed automatically
    pub parse_seqids: bool,
    /// 0x01 Exclude Model
    // NCBI's integer (`AsInteger`, an `int`: `0x01` is 1, `-1` sets every bit), made an `Uint8`
    // (makeblastdb.cpp:1227-1230).
    pub oid_masks: Option<i32>,
    /// Copy the sequences of -input_type blastdb in a random order, seeded by the clock's second
    // NCBI's argument is hidden (`fHidden`, makeblastdb.cpp:253-254): its help lists it with
    // `-help-full` only. LOSAT's help is NCBI's text (`cli::reader`): `-help-full` lists it, `-help` does not.
    pub randomize: bool,
}

impl MakedbArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, makeblastdb's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(MakedbArgs {
            input: args.option_or_default("in")?,
            input_type: args.option_or_default("input_type")?,
            dbtype: args.option_or_default("dbtype")?,
            title: args.given_option("title")?,
            out: args.given_option("out")?,
            blastdb_version: args.option_or_default("blastdb_version")?,
            max_file_sz: args.option_or_default("max_file_sz")?,
            metadata_output_prefix: args.given_option("metadata_output_prefix")?,
            logfile: args.given_option("logfile")?,
            hash_index: args.flag("hash_index"),
            taxid: args.given_option("taxid")?,
            taxid_map: args.given_option("taxid_map")?,
            mask_data: args.given_option("mask_data")?,
            mask_id: args.given_option("mask_id")?,
            mask_desc: args.given_option("mask_desc")?,
            gi_mask: args.flag("gi_mask"),
            gi_mask_name: args.given_option("gi_mask_name")?,
            parse_seqids: args.flag("parse_seqids"),
            oid_masks: args.given_option("oid_masks")?,
            randomize: args.flag("randomize"),
        })
    }
}

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
    /// makeblastdb prints after "Command line argument error: ", as bytes.
    Argument(Vec<u8>),
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
            MakedbError::Build(WriteDbError::Native(native)) => native.exit,
            MakedbError::Build(
                WriteDbError::Exception(_)
                | WriteDbError::Other(_)
                | WriteDbError::Lmdb(_)
                | WriteDbError::Std(_)
                | WriteDbError::Io { .. },
            ) => 255,
            MakedbError::Build(WriteDbError::OutOfMemory) => 4,
            MakedbError::Build(WriteDbError::ClosedPipe) => 0,
            MakedbError::Output(_) => 6,
            MakedbError::Database(error) if !is_refusal(error) => 2,
            MakedbError::Build(WriteDbError::Database(_)) => 2,
            MakedbError::Database(_)
            | MakedbError::Build(
                WriteDbError::Creation(_)
                | WriteDbError::CreationBytes(_)
                | WriteDbError::GiOnlyBetweenKeys { .. }
                | WriteDbError::CloseRefusal(_)
                | WriteDbError::Options(_)
                | WriteDbError::OptionsBytes(_)
                | WriteDbError::Unsupported(_)
                | WriteDbError::OrdinalIdWithoutIds,
            )
            | MakedbError::Unsupported(_)
            | MakedbError::Argument(_) => 1,
        }
    }
}

impl MakedbError {
    /// Whether the run ended at a closed pipe ([`WriteDbError::ClosedPipe`]).
    fn is_closed_pipe(&self) -> bool {
        matches!(self, MakedbError::Build(error) if error.is_closed_pipe())
    }

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
                [b"Command line argument error: ".as_slice(), message].concat()
            }
        }
    }
}

impl std::fmt::Display for MakedbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.message()))
    }
}

/// Whether `error` is LOSAT's refusal of a database where NCBI crashes, throws what LOSAT does
/// not report or gives a result that depends on its memory ([`SeqDbError::NcbiCrashes`],
/// [`SeqDbError::NcbiThrows`], [`SeqDbError::Unsupported`]), which is printed as the other
/// refusals of makeblastdb are, with exit status 1.
fn is_refusal(error: &SeqDbError) -> bool {
    matches!(
        error,
        SeqDbError::NcbiCrashes { .. }
            | SeqDbError::NcbiThrows { .. }
            | SeqDbError::Unsupported { .. }
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

fn options_error_bytes(message: Vec<u8>) -> MakedbError {
    MakedbError::Build(WriteDbError::OptionsBytes(message))
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
        SeqDbError::Setting(refused) => {
            unsupported(format!("{refused} is not supported by LOSAT's makeblastdb"))
        }
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

/// Whether `error` is LOSAT's own (a refusal, a failed write of LOSAT's), which removes every
/// file of the build (`BuildDatabase::erase_after_refusal`), rather than an error of NCBI's,
/// after which makeblastdb's writer is closed when it is destroyed
/// (`BuildDatabase::close_at_destruction`) and the files that NCBI leaves stay.
fn is_losat_error(error: &MakedbError) -> bool {
    match error {
        MakedbError::Unsupported(_) => true,
        MakedbError::Build(error) => error.is_losat_own(),
        _ => false,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1762-1807
// makeblastdb splits `-in` without keeping the quotes (`keep_quote` false).
/// `SeqDB_SplitQuoted`: the input names of `-in` (spaces separate names,
/// double quotes keep a name with spaces together).
fn split_quoted(names: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut quoted = false;
    let mut begin = 0;
    for (i, &ch) in names.iter().enumerate() {
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
fn combine_and_quote(names: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    for name in names {
        if !out.is_empty() {
            out.push(b' ');
        }
        if name.contains(&b' ') {
            out.push(b'"');
            out.extend_from_slice(name);
            out.push(b'"');
        } else {
            out.extend_from_slice(name);
        }
    }
    out
}

/// Makeblastdb's own environment variables (`makeblastdb.cpp:1130-1136`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MakedbEnv {
    /// `NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS` non-empty.
    pub(crate) skip_unverified: bool,
    /// `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` unset or empty.
    pub(crate) scan_cfastareader: bool,
}

impl Default for MakedbEnv {
    fn default() -> Self {
        Self {
            skip_unverified: false,
            scan_cfastareader: true,
        }
    }
}

impl MakedbEnv {
    /// Reads makeblastdb's two environment variables via `env`.
    pub(crate) fn from_env(env: &dyn Fn(&str) -> Option<std::ffi::OsString>) -> Self {
        Self {
            skip_unverified: env("NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS")
                .is_some_and(|value| !value.is_empty()),
            scan_cfastareader: env(
                "NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS",
            )
            .is_none_or(|value| value.is_empty()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1128-1136
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:562-566
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1733-1738,1764-1766
// makeblastdb reads two variables of its own. ASN.1 inputs (`asn1_txt`, `asn1_bin`) honor both
// variables directly (`MakedbEnv`). Neither is read for `-input_type blastdb`.
//
// * `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS`: for FASTA input, LOSAT
//   reproduces makeblastdb with this variable unset and refuses to run with it set.
// * `NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS`: read by `CSeqEntrySource` only (ASN.1 input),
//   so it changes nothing for FASTA and is not refused there.
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
        // `NStr::StringToUInt` of `value` with no flags: no limit where it throws.
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

/// The format version of `-blastdb_version`, which NCBI's argument reader allows only as 4 or 5
/// (`CArgAllow_Integers` from `eBDB_Version4` to `eBDB_Version5`); a library caller that set the
/// field itself gets LOSAT's refusal of another number.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1218-1219
fn db_version(args: &MakedbArgs) -> Result<DbVersion, MakedbError> {
    DbVersion::from_number(args.blastdb_version).ok_or_else(|| {
        unsupported(format!(
            "-blastdb_version {} is not a version makeblastdb writes (4 or 5)",
            args.blastdb_version
        ))
    })
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
    let name_bytes = crate::utils::os_bytes::os_str_bytes(name);
    let not_accessible = || {
        let mut msg = b"Argument \"taxid_map\". File is not accessible".to_vec();
        if !name_bytes.is_empty() {
            msg.extend_from_slice(b":  `");
            msg.extend_from_slice(&name_bytes);
            msg.push(b'\'');
        }
        MakedbError::Argument(msg)
    };
    let mut file;
    let mut standard_input;
    let is_stdin = name_bytes.as_ref() == b"-";
    let reader: &mut dyn BufRead = if is_stdin {
        standard_input = std::io::BufReader::new(stdin);
        &mut standard_input
    } else if name_bytes.is_empty() {
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
                "line {line} of the -taxid_map file {} holds only white space, which \
                 makeblastdb 2.17.0 cannot read (it crashes)",
                String::from_utf8_lossy(&name_bytes)
            )),
        })?;
    Ok((Some(taxids), is_stdin.then_some(failed)))
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
    let raw = crate::utils::os_bytes::os_str_bytes(&args.max_file_sz);
    let bytes = string_to_uint8_data_size_bytes(&raw)
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
fn compressed_fault_error(name: &[u8], compression: Compression, fault: &str) -> MakedbError {
    unsupported(format!(
        "-in: {}",
        compression.not_data(&String::from_utf8_lossy(name), fault)
    ))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:792-836
// Every file is checked, in the order of `-in`, before the first one is read.
// `CFile::Exists` is true for a regular file, also through a link; anything
// else "does not exist". The test of the format opens the file again, so the
// build reads it from its start. `ERR_POST` ends the message with a line feed.
/// `x_VerifyInputFilesType` for FASTA and ASN.1 inputs (`check_fasta` is true only for FASTA).
fn verify_input_files(
    names: &[&[u8]],
    check_fasta: bool,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    for &seq_file in names {
        let path = crate::utils::os_bytes::path_of_bytes(seq_file);
        match std::fs::metadata(&path) {
            Ok(meta) if meta.is_file() && meta.len() == 0 => {
                return Err(options_error_bytes(
                    [b"File ".as_slice(), seq_file, b" is empty"].concat(),
                ));
            }
            Ok(meta) if meta.is_file() => {}
            _ => {
                return Err(options_error_bytes(
                    [b"File ".as_slice(), seq_file, b" does not exist"].concat(),
                ));
            }
        }
        if !check_fasta {
            continue;
        }
        // A compressed name is guessed on its inflated bytes, as `x_VerifyInputFilesType` guesses
        // on its `CDecompressIStream` (makeblastdb.cpp:815-830); a fault of the compressed data
        // reads as "not FASTA".
        let looks_like_fasta =
            std::fs::File::open(&path).is_ok_and(|file| match compression_of_name(&path) {
                Some(compression) => appears_fasta(CompressedInput::new(file, compression)),
                None => appears_fasta(file),
            });
        if !looks_like_fasta {
            let _ = err.write_all(
                &[
                    b"Warning: [makeblastdb] \nInput file ".as_slice(),
                    seq_file,
                    b" does NOT appear to be FASTA (processing anyway).\n\
                      Advise validating database with 'blastdbcheck -dbtype [prot|nucl] -db ${DBNAME}'\n\n",
                ]
                .concat(),
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
    env: MakedbEnv,
    settings: &ApplicationSettings,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    cout: Option<&mut dyn Write>,
    err: &mut dyn Write,
    logfile_error: Option<MakedbError>,
    end_messages: &mut Vec<Vec<u8>>,
) -> Result<(), MakedbError> {
    let registry_blastdb = registry_blastdb(settings);
    let version = db_version(args)?;
    let is_protein = args.dbtype == "prot";
    let blastdb_input = args.input_type == "blastdb";
    let asn_format = match args.input_type.as_str() {
        "asn1_txt" => Some(AsnFormat::Text),
        "asn1_bin" => Some(AsnFormat::Binary),
        _ => None,
    };

    let input_bytes = crate::utils::os_bytes::os_str_bytes(&args.input);
    let out_bytes = args
        .out
        .as_deref()
        .map(crate::utils::os_bytes::os_str_bytes);
    let dbname: &[u8] = out_bytes.as_deref().unwrap_or(&input_bytes);
    if blastdb_input && dbname == input_bytes.as_ref() {
        return Err(options_error(
            "Cannot create a BLAST database from an existing one without changing the output \
             name, please provide a (different) database name using -out",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1156-1159
    if args.randomize && !blastdb_input {
        return Err(options_error(
            "Option randomize is only applicable to blastdb input type",
        ));
    }
    // The input types LOSAT does not read are refused after NCBI's checks of the input type.
    if args.input_type != "fasta" && !blastdb_input && asn_format.is_none() {
        return Err(unsupported(format!(
            "-input_type {} is not supported; only fasta, blastdb, asn1_txt and asn1_bin input are",
            args.input_type
        )));
    }
    // NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:264-275,282-301
    // `NStr::Split` of `dbname` at `" "` into `input_files` runs with no flag
    // (ncbistr.hpp:2551-2555): adjacent delimiters are not merged and none is
    // truncated, so every space ends a token, an empty one included ("db "
    // is "db" and "", " db" is "" and "db"). A name with a space anywhere
    // gives more than one token; the empty name gives none.
    if dbname == b"-" || dbname.contains(&b' ') {
        return Err(options_error("Please provide a database name using -out"));
    }
    let title = match args.title.as_deref() {
        Some(title) => crate::utils::os_bytes::os_str_bytes(title).into_owned(),
        None if blastdb_input => input_database_title(&input_bytes, is_protein, registry_blastdb)?,
        None => input_bytes.as_ref().to_vec(),
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
    let names = split_quoted(&input_bytes);
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
            // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1005-1014
            full: args.parse_seqids || args.input_type != "fasta",
            hash: args.hash_index,
            trace: false,
        },
    )?;
    if let Some(cout) = cout {
        db.set_cout(cout);
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1205,1231-1241
    // makeblastdb gives `-gi_mask` to the writer when it makes it, before any masking algorithm is
    // registered.
    if args.gi_mask {
        db.use_gi_masks();
    }
    if let Some(oid_masks) = args.oid_masks {
        db.set_oid_masks(oid_masks as u64);
    }
    let max_file_size = max_file_size(args)?;
    let _ = write!(
        db.log(),
        "Maximum file size: {}",
        uint8_to_string_data_size(max_file_size)
    );
    endl(db.log())?;
    db.set_max_file_size(max_file_size);
    let (taxids, stdin_map) = read_taxids(args, stdin, err)?;
    if let Some(taxids) = taxids {
        db.set_taxids(taxids);
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1288-1294
    // The mask files are read before the input is looked at. A failure ends the run at once
    // (`CATCH_ALL`), and the destructor's messages follow it (`end_messages`).
    let mask_data = args
        .mask_data
        .as_deref()
        .map(crate::utils::os_bytes::os_str_bytes);
    let mask_id = args
        .mask_id
        .as_deref()
        .map(crate::utils::os_bytes::os_str_bytes);
    let mask_desc = args
        .mask_desc
        .as_deref()
        .map(crate::utils::os_bytes::os_str_bytes);
    let gi_mask_name = args
        .gi_mask_name
        .as_deref()
        .map(crate::utils::os_bytes::os_str_bytes);
    let mask_options = MaskDataOptions {
        mask_data: mask_data.as_deref(),
        mask_id: mask_id.as_deref(),
        mask_desc: mask_desc.as_deref(),
        gi_mask_name: gi_mask_name.as_deref(),
        root: None,
    };
    let masked = {
        let mut warn = warn_to(err);
        process_mask_data(&mut db, &mask_options, &mut warn)
    };
    if let Err(error) = masked {
        let error = error.into_makedb_error(&settings.read_verification);
        *end_messages = destructor_messages(&db, error.as_ref());
        return error.map_or(Ok(()), Err);
    }
    let metadata_prefix = args
        .metadata_output_prefix
        .as_deref()
        .map(crate::utils::os_bytes::os_str_bytes);
    let user_path = metadata_user_path(metadata_prefix.as_deref());
    let result = if blastdb_input {
        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:665-668
        // The seed of `-randomize` is the seconds of the clock (`time` with a null pointer), which
        // `mt19937` takes modulo 2^32. NCBI reads the clock again when the copy starts, after the
        // input's headers are read; LOSAT reads it once for the build (`now`), so that a build
        // whose clock is fixed (the tests) has one instant. Under a clock that runs, the two may
        // take different seconds.
        let randomize = args.randomize.then_some(now.unix_seconds as u32);
        copy_input_databases(
            &mut db,
            &names,
            &user_path,
            is_protein,
            registry_blastdb,
            &settings.read_verification,
            randomize,
            err,
        )
    } else if names.is_empty() {
        let what = if asn_format.is_some() {
            "ASN.1"
        } else {
            "FASTA"
        };
        Err(unsupported(format!(
            "-in names no input file; losat makeblastdb reads {what} files or the standard input"
        )))
    } else if let Some(format) = asn_format {
        add_asn_and_close(
            &mut db,
            &names,
            (names[0] == b"-").then_some(stdin_map),
            format,
            env,
            &user_path,
            is_protein,
            registry_blastdb,
            &settings.read_verification,
            stdin,
            err,
        )
    } else {
        add_input_and_close(
            &mut db,
            &names,
            (names[0] == b"-").then_some(stdin_map),
            &user_path,
            is_protein,
            registry_blastdb,
            stdin,
            err,
        )
    };
    // One rule for every input and every point of a refusal: a build that LOSAT refuses leaves
    // no file of its own.
    if result.as_ref().is_err_and(is_losat_error) {
        db.erase_after_refusal();
    }
    *end_messages = destructor_messages(&db, result.as_ref().err());
    // NCBI's error ends a copy or an ASN.1 build without `EndBuild`: the writer is closed when
    // it is destroyed. A FASTA input's error has ended the build (`AddFasta`).
    if (blastdb_input || asn_format.is_some())
        && result
            .as_ref()
            .is_err_and(|error| !is_losat_error(error) && !error.is_closed_pipe())
    {
        end_messages.extend(db.close_at_destruction());
    }
    result
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1175-1188
// The names are split as the reader splits a list; with no name at all NCBI opens the first of
// none, which LOSAT refuses.
/// The title of `-input_type blastdb` without `-title`: the title of the one database `input`
/// names (`GetTitle`).
fn input_database_title(
    input: &[u8],
    is_protein: bool,
    registry_blastdb: Option<&[u8]>,
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
    let search = SearchPath::from_environment(registry_blastdb);
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
/// (`open_new_database`); `read_verification` is the run's `SERIAL_VERIFY_DATA_READ`
/// (`ApplicationSettings::read_verification`), `randomize` the seed of `-randomize`.
#[allow(clippy::too_many_arguments)]
fn copy_input_databases(
    db: &mut BuildDatabase<'_>,
    names: &[&[u8]],
    user_path: &[u8],
    is_protein: bool,
    registry_blastdb: Option<&[u8]>,
    read_verification: &ReadVerification,
    randomize: Option<u32>,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let seqtype = seq_type(is_protein);
    let search = SearchPath::from_environment(registry_blastdb);
    for &name in names {
        match SeqDb::open_in(name, seqtype, &search) {
            Ok(_) => {}
            Err(error) if error.ncbi_message().is_some() => {
                let _ = err.write_all(
                    &[
                        b"Error: [makeblastdb] Unable to open input ".as_slice(),
                        name,
                        b" as BLAST db\n",
                    ]
                    .concat(),
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
    let mut input =
        SeqDb::open_in(&combine_and_quote(names), seqtype, &search).map_err(input_error)?;
    // The definition lines are read as the copy reaches them: a record that fails NCBI's check
    // of a mandatory member under a `SERIAL_VERIFY_DATA_READ` that turns it off is refused there
    // (`SeqDb::refuse_unchecked_reads`, `input_error`).
    input.refuse_unchecked_reads(read_verification);
    let leafs = {
        let mut warn = warn_to(err);
        scan_leaf_tax_ids(&input, &mut warn, db.utf8_once())
    };
    db.set_leaf_tax_ids(leafs.map_err(raw_copy_error)?)?;
    let mut source = RawSeqDbSource::new(&input, db, randomize).map_err(raw_copy_error)?;
    let copied = {
        let mut warn = warn_to(err);
        db.add_raw_sequences(&mut source, &mut warn)
    };
    if let Err(error) = copied {
        let error = raw_copy_error(error);
        if is_losat_error(&error) {
            let mut warn = warn_to(err);
            // A closed pipe at the log of that `EndBuild` ends the run there (`endl`); the
            // refused build leaves no file all the same (`build_database_with_log`).
            if let Err(stop) = db.end_build(true, &mut warn) {
                if stop.is_closed_pipe() {
                    db.erase_after_refusal();
                    return Err(stop.into());
                }
            }
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
    messages.extend(db.taxid_message());
    messages
}

/// The rest of `x_BuildDatabase`: the input files, `EndBuild` and the metadata file.
/// `from_stdin` is `Some(map)` for `-in -`, `map` being `Some(failed)` when `-taxid_map -`
/// read the standard input first (`failed`: whether it left the stream failed). `user_path`
/// comes before the names of the metadata's files (`metadata_user_path`).
#[allow(clippy::too_many_arguments)]
fn add_input_and_close(
    db: &mut BuildDatabase<'_>,
    names: &[&[u8]],
    from_stdin: Option<Option<bool>>,
    user_path: &[u8],
    is_protein: bool,
    registry_blastdb: Option<&[u8]>,
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
        verify_input_files(names, true, err)?;
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
            let path = crate::utils::os_bytes::path_of_bytes(name);
            let Ok(mut file) = std::fs::File::open(&path) else {
                continue;
            };
            // A compressed name is inflated while it is read (`CompressedInput`); a fault of
            // its data is reported by name, after the build removed its files as for any
            // failed input.
            match compression_of_name(&path) {
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

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1026-1055,559-578
/// `x_ProcessInputData` + `EndBuild` + `open_new_database` for `-input_type asn1_txt` and
/// `-input_type asn1_bin`. Unlike FASTA (`AddFasta`), an error from NCBI leaves the build open so
/// `build_database_with_log` runs `close_at_destruction` (writing any sequences already added);
/// a refusal of LOSAT's (including `read_verification` refusals for `SERIAL_VERIFY_DATA_READ` and
/// `SERIAL_VERIFY_DATA_GET`) writes no `EndBuild` log line, and `build_database_with_log` removes
/// every file of the build (`erase_after_refusal`).
#[allow(clippy::too_many_arguments)]
fn add_asn_and_close(
    db: &mut BuildDatabase<'_>,
    names: &[&[u8]],
    from_stdin: Option<Option<bool>>,
    format: AsnFormat,
    env: MakedbEnv,
    user_path: &[u8],
    is_protein: bool,
    registry_blastdb: Option<&[u8]>,
    read_verification: &ReadVerification,
    stdin: &mut dyn Read,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let feed_one = |db: &mut BuildDatabase<'_>,
                    reader: &mut dyn Read,
                    err: &mut dyn Write|
     -> Result<(), MakedbError> {
        let res = {
            let mut warn = warn_to(err);
            db.add_seq_entries(
                reader,
                format,
                env.skip_unverified,
                env.scan_cfastareader,
                &mut warn,
            )
        };
        res.map_err(|error| {
            let error = match error {
                WriteDbError::Exception(exception) => match read_verification
                    .refused(&exception, "ASN.1 input (-input_type asn1_txt|asn1_bin)")
                {
                    Some(refused) => {
                        unsupported(format!("{refused} is not supported by LOSAT's makeblastdb"))
                    }
                    None => MakedbError::Build(WriteDbError::Exception(exception)),
                },
                other => MakedbError::Build(other),
            };
            error
        })
    };

    if let Some(stdin_map) = from_stdin {
        // After `-taxid_map -` drained stdin, NCBI's `x_AddSeqEntries` sees EOF, posts
        // `No sequences written`, and ends with `No volumes were created.` (exit 0).
        let mut empty = std::io::empty();
        let reader: &mut dyn Read = match stdin_map {
            None => stdin,
            Some(_) => &mut empty,
        };
        feed_one(db, reader, err)?;
    } else {
        verify_input_files(names, false, err)?;
        for &name in names {
            let path = crate::utils::os_bytes::path_of_bytes(name);
            let Ok(mut file) = std::fs::File::open(&path) else {
                // File exists and is non-empty, but unreadable (e.g. mode 000): NCBI's
                // `CNcbiIfstream` fails and `x_AddSeqEntries` hits EOF (`No sequences written`).
                let mut empty = std::io::empty();
                let _ = feed_one(db, &mut empty, err);
                continue;
            };
            match compression_of_name(&path) {
                Some(compression) => {
                    let mut input = CompressedInput::new(file, compression);
                    let res = feed_one(db, &mut input, err);
                    if let Some(fault) = input.fault() {
                        return Err(compressed_fault_error(name, compression, fault));
                    }
                    res?;
                }
                None => feed_one(db, &mut file, err)?,
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
// The `CSeqDB` `sdb` made from `new_db` and `t` throws a `CSeqDBException` when the reader cannot
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
    user_path: &[u8],
    is_protein: bool,
    registry_blastdb: Option<&[u8]>,
    err: &mut dyn Write,
) -> Result<(), MakedbError> {
    let new_db = db.output_dbname();
    let seqtype = seq_type(is_protein);
    // `CSeqDB` names its search path in the message of a database that is not there: the
    // current directory, `BLASTDB` and the registry entry `[BLAST] BLASTDB` (`registry_blastdb`,
    // of the application settings that the callers of `run_with` read before the build).
    let search = SearchPath::from_environment(registry_blastdb);
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
    let dbname = match new_db
        .iter()
        .rposition(|&byte| byte == std::path::MAIN_SEPARATOR as u8)
    {
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
    // The name of a database makeblastdb built is not empty (`build_lmdb_file_name`'s error).
    let lmdb: Vec<std::path::PathBuf> = version5
        .then(|| build_lmdb_file_name(new_db, is_protein).ok())
        .flatten()
        .map(crate::utils::os_bytes::path_from_bytes)
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
        dbname: dbname.to_vec(),
        dbtype: if is_protein { "Protein" } else { "Nucleotide" },
        db_version: if version5 { 5 } else { 4 },
        description: located.title,
        number_of_letters: totals.iter().map(|(_, letters)| letters).sum(),
        // NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.def:11-12
        // NCBI's `GetDBMetaData` counts in `num_seqs`, an `int`
        // (c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1707), which a database of more than
        // 2^31 - 1 sequences would overflow; LOSAT sums in `Int8`.
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
// The log is `cout`, and every line of it ends with a new line and a flush (`build_db::endl`),
// so the log and the messages of the standard error stream reach a common destination (`2>&1`)
// in the order of the events. `CoutLog` writes it as NCBI's `cout` (`writedb::cout`) over the
// process's standard output with no buffer of Rust's (`cli::standard_output`). The buffer
// writes a line where NCBI's does, a write that fails loses what the buffer held, and the
// stream writes nothing after it (`badbit`; line RG item KC-11b-2).
/// Runs `losat makeblastdb` with the settings of its application start and returns the exit code.
pub fn run(args: &MakedbArgs, settings: &ApplicationSettings) -> i32 {
    let stdin = std::io::stdin();
    let stderr = std::io::stderr();
    let mut device = crate::cli::standard_output();
    let mut out = CoutLog(super::cout::cout(
        &mut device,
        crate::cli::NcbiStdio::cout(),
    ));
    let mut err = stderr.lock();
    let now = LocalTime::now();
    let status = run_with(
        args,
        settings,
        now.as_ref(),
        Platform::host(),
        &mut stdin.lock(),
        &mut out,
        &mut err,
    );
    // The exit's flush of `cout`.
    let _ = out.flush();
    status
}

/// makeblastdb's log on NCBI's `cout`: a write is an insertion, a flush is the stream's flush
/// (the end of each line, `build_db::endl`). The stream's state is not reported (NCBI's
/// makeblastdb does not read it), but for a write that found the pipe closed where SIGPIPE ends
/// NCBI's process (`StreamBuf::closed_pipe`): from that write on, every write and flush returns
/// the error of a closed pipe, and the build ends at the end of the line (`build_db::endl`).
struct CoutLog<'a>(super::cout::Cout<'a>);

impl CoutLog<'_> {
    /// The error of a closed pipe once a write found one.
    fn closed_pipe(&self) -> std::io::Result<()> {
        if self.0.rdbuf().closed_pipe() {
            return Err(std::io::ErrorKind::BrokenPipe.into());
        }
        Ok(())
    }
}

impl Write for CoutLog<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.insert(buf);
        self.closed_pipe().map(|()| buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush();
        self.closed_pipe()
    }
}

/// The registry entry `[BLAST] BLASTDB` of `settings` that `run_with` takes (its bytes): `None`
/// when the registry has no such entry (`ApplicationSettings::blastdb`).
pub(crate) fn registry_blastdb(settings: &ApplicationSettings) -> Option<&[u8]> {
    settings.blastdb.as_deref()
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
/// One run with the application settings `settings` (the registry entry `[BLAST] BLASTDB`,
/// `registry_blastdb`, and the checks of NCBI's reads, `ApplicationSettings::read_verification`)
/// and the given time, platform and streams: the build, then the error message and the exit code
/// of a failed build. With `-logfile` the log and the messages go where makeblastdb sends them
/// instead of `out` and `err`.
pub(crate) fn run_with(
    args: &MakedbArgs,
    settings: &ApplicationSettings,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    run_with_env(
        args,
        MakedbEnv::from_env(&|name| std::env::var_os(name)),
        settings,
        now,
        platform,
        stdin,
        out,
        err,
    )
}

/// [`run_with`] with explicit [`MakedbEnv`] flags (so tests can cover both environment variables
/// without mutating process environment).
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_with_env(
    args: &MakedbArgs,
    env: MakedbEnv,
    settings: &ApplicationSettings,
    now: Option<&LocalTime>,
    platform: Platform,
    stdin: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    // The run's own diagnostic stream of the main thread, which the object manager fails on the
    // ASN.1 input (`objmgr::annot_index`): its later posts are their headers alone
    // (`ncbi_exception::post_text`), and its end the line feed of its report (`catch_all_report`).
    let _diag = crate::ncbi_exception::MainDiag::enter(Some(std::sync::Arc::new(
        crate::ncbi_exception::MainDiag::new(),
    )));
    let mut end_messages = Vec::new();
    let name = match args.logfile.as_deref() {
        None => {
            let result = build_database_with_log(
                args,
                env,
                settings,
                now,
                platform,
                stdin,
                out,
                None,
                err,
                None,
                &mut end_messages,
            );
            return finish(result, &end_messages, err);
        }
        Some(name) if crate::utils::os_bytes::os_str_bytes(name).as_ref() == b"-" => {
            let result = build_database_with_log(
                args,
                env,
                settings,
                now,
                platform,
                stdin,
                out,
                None,
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
    let target = diag_target(Some(name));
    let log = log_of(name);
    let mut diag = DiagStream {
        target: &target,
        stderr: err,
    };
    match log {
        Some(log) => {
            let mut log = LogLines {
                file: log,
                line: Vec::new(),
            };
            let result = build_database_with_log(
                args,
                env,
                settings,
                now,
                platform,
                stdin,
                &mut log,
                Some(out),
                &mut diag,
                None,
                &mut end_messages,
            );
            let _ = log.flush();
            finish(result, &end_messages, &mut diag)
        }
        None => {
            let mut msg = b"Argument \"logfile\". File is not accessible:  `".to_vec();
            msg.extend_from_slice(&crate::utils::os_bytes::os_str_bytes(name));
            msg.push(b'\'');
            let error = MakedbError::Argument(msg);
            let result = build_database_with_log(
                args,
                env,
                settings,
                now,
                platform,
                stdin,
                out,
                None,
                &mut diag,
                Some(error),
                &mut end_messages,
            );
            finish(result, &end_messages, &mut diag)
        }
    }
}

/// The file of `-logfile` as NCBI's log stream writes it: each line, ended by its `NcbiEndl`, in
/// one write (line RG KD-11b2: the pieces of a line are not written one by one); the bytes of a
/// line not ended yet at the end.
struct LogLines {
    file: std::fs::File,
    /// The bytes of the line being written.
    line: Vec<u8>,
}

impl Write for LogLines {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let mut rest = buf;
        while let Some(at) = rest.iter().position(|&byte| byte == b'\n') {
            self.line.extend_from_slice(&rest[..=at]);
            rest = &rest[at + 1..];
            let line = std::mem::take(&mut self.line);
            self.file.write_all(&line)?;
        }
        self.line.extend_from_slice(rest);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if !self.line.is_empty() {
            let line = std::mem::take(&mut self.line);
            self.file.write_all(&line)?;
        }
        self.file.flush()
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
fn log_of(name: &std::ffi::OsStr) -> Option<std::fs::File> {
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
/// `~CBuildDatabase` (`end_messages`); nothing after a closed pipe.
fn finish(result: Result<(), MakedbError>, end_messages: &[Vec<u8>], err: &mut dyn Write) -> i32 {
    let status = match result {
        Ok(()) => 0,
        // A closed pipe ends the run where SIGPIPE ends NCBI's process: no message, and none of
        // `~CBuildDatabase`'s (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10).
        Err(error) if error.is_closed_pipe() => return error.exit_code(),
        Err(error) => {
            let _ = err.write_all(&catch_all_report(&error));
            error.exit_code()
        }
    };
    for message in end_messages {
        let _ = err.write_all(&[message.as_slice(), b"\n"].concat());
    }
    status
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-263
/// What the run writes for `error` at its end: its message and a line feed. NCBI's end is its
/// general handler's `LOG_POST`, of which a failed diagnostic stream of the main thread writes the
/// line feed alone (`cli::catch_all_text`; the object manager fails it on the ASN.1 input,
/// `objmgr::annot_index`); the end of the failed write of an object's text writes nothing
/// (`WriteDbError::Native`). LOSAT's refusals keep their words.
fn catch_all_report(error: &MakedbError) -> Vec<u8> {
    if let MakedbError::Build(WriteDbError::Native(native)) = error {
        return crate::cli::catch_all_text(&native.message).to_vec();
    }
    let report = [error.message().as_slice(), b"\n"].concat();
    if is_losat_error(error) {
        return report;
    }
    crate::cli::catch_all_text(&report).to_vec()
}

#[cfg(test)]
#[path = "makeblastdb_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "makeblastdb_asn1_tests.rs"]
mod asn1_tests;
