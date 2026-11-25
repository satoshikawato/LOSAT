//! makeprofiledb's `Init` (`x_InitProgramParameters`), `Run` (`x_Run`, `x_MakeVol`, the alias
//! file of several volumes) and the clean-up of its destructor, for `-dbtype rps`.

use std::io::{Read, Write};
use std::path::Path;

use crate::blastinput::ncbi_environment::{ApplicationSettings, ReadVerification};
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{self, Cin, OutputFile};
use crate::ncbi_app::Diag;
use crate::report::run_files::OpenedFile;
use crate::report::serial::asn_types::{Post, ReadError};
use crate::report::serial::pssm::{PssmWithParameters, SeqData, SeqEntry};
use crate::seqid::{BlastDefline, DeflineSet};
use crate::utils::os_bytes::{path_from_bytes, path_of_bytes};
use crate::writedb::alias_file::{self, AliasError, FilterType, CLOCK_ERROR};
use crate::writedb::bioseq::{self, Bioseq, ParsedIds};
use crate::writedb::build_db::{create_directories, delete_blast_db};
use crate::writedb::convert::{iupacna_high_byte, short_nucleotide_data, ShortNucleotideData};
use crate::writedb::deflines::{bioseq_deflines, first_org_taxid, BioseqDeflinesError};
use crate::writedb::files::{make_short_name, DbVersion};
use crate::writedb::taxid_set::{TaxIdMapError, TaxIdSet};
use crate::writedb::time::LocalTime;
use crate::writedb::volume::Indices;
use crate::writedb::writer::WriteDb;
use crate::writedb::{Platform, WriteDbError};

use super::profile::{
    check_input_scoremat, fill_in_rps_db_parameters, query_length, update_rps_db_info,
    update_statistics, ProfileError, Scoremat,
};
use super::rps_files::{self, ProfileData, RpsError, RpsVolume, VolumeParams};

/// The limit of each file of a volume's database (`SetMaxFileSize` of 4,000,000,000 bytes).
const MAX_FILE_SIZE: u64 = 4_000_000_000;

/// What makeprofiledb's `Run` reports, by the clause of its handler that catches it.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1879-1904
#[derive(Debug)]
enum RunError {
    /// `CInputException`: "INPUT ERROR: " and the message, `BLAST_INPUT_ERROR` (1).
    Input(Vec<u8>),
    /// `CSeqDBException`: "ERROR: ", `BLAST_DATABASE_ERROR` (2).
    SeqDb(Vec<u8>),
    /// `CBlastException`: "ERROR: ", `BLAST_INPUT_ERROR` (1).
    Blast(Vec<u8>),
    /// Another `CException`: "ERROR: ", `BLAST_UNKNOWN_ERROR` (255).
    Exception(Vec<u8>),
    /// Any other exception: "Error: Unknown exception", 255.
    Unknown,
    /// LOSAT's refusal of `what` (a noun phrase), in the form of LOSAT's other programs
    /// ([`unsupported`]): a line on the standard error stream, exit status 1.
    Refused(String),
    /// A failure of LOSAT's own (a file that could not be written, a clock that could not be
    /// read), where NCBI has no message: the text on the standard error stream, exit status 1.
    Failed(String),
}

/// LOSAT's refusal of `what` (a noun phrase) in makeprofiledb: `<what> is not supported by
/// LOSAT's MAKEPROFILEDB`, the form of LOSAT's other programs.
fn unsupported(what: &str) -> String {
    crate::seqmasks_io::run::unsupported(what, super::PROGRAM)
}

impl RunError {
    /// The text that the report of this error writes after "Error: " (NCBI's prefix and message,
    /// LOSAT's refusal), for a refusal that names it.
    fn text(&self) -> String {
        let bytes =
            |prefix: &str, message: &[u8]| format!("{prefix}{}", String::from_utf8_lossy(message));
        match self {
            RunError::Input(message) => bytes("INPUT ERROR: ", message),
            RunError::SeqDb(message) | RunError::Blast(message) | RunError::Exception(message) => {
                bytes("ERROR: ", message)
            }
            RunError::Unknown => "Unknown exception".to_string(),
            RunError::Refused(what) => unsupported(what),
            RunError::Failed(text) => text.clone(),
        }
    }
}

impl From<RpsError> for RunError {
    fn from(error: RpsError) -> Self {
        match error {
            RpsError::SeqDb(message) => RunError::SeqDb(message.into()),
            RpsError::Blast(message) => RunError::Blast(message.into()),
            RpsError::Unknown => RunError::Unknown,
            RpsError::Refused(reason) => RunError::Refused(reason),
            RpsError::Io(error) => RunError::Failed(format!(
                "writing a file of the RPS database failed: {error}"
            )),
        }
    }
}

impl From<ProfileError> for RunError {
    fn from(error: ProfileError) -> Self {
        match error {
            ProfileError::Input(message) => RunError::Input(message),
            ProfileError::Rps(error) => error.into(),
            ProfileError::Blast(message) => RunError::Blast(message.into_bytes()),
            ProfileError::Exception(exception) => RunError::Exception(exception.get_msg().to_vec()),
            ProfileError::Unknown => RunError::Unknown,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:69-80 (the
// exceptions of CWriteDB: `CWriteDBException`, a `CException`; `CSeqDBException` of LMDB)
/// The writer's errors as makeprofiledb's handler reports them: its exceptions are
/// `CException`s ("ERROR: " and the message, 255), those of LMDB's environment
/// `CSeqDBException`s (2).
// A `std::exception` (`std::invalid_argument` of an empty name's LMDB file, `std::bad_alloc`) is
// the handler's unknown exception.
impl From<WriteDbError> for RunError {
    fn from(error: WriteDbError) -> Self {
        match error {
            // A GI-only query between queries with other ids in a version 5 database: NCBI's
            // writer aborts at its close (survey I01g, I01x), as makeblastdb's does, which
            // LOSAT refuses (makedb-v5 DECISIONS 10).
            WriteDbError::GiOnlyBetweenKeys { oid, name } => RunError::Refused(format!(
                "a version 5 database with -index true in which the query {} (OID {oid}) has \
                 only GI ids and lies between queries with other ids (NCBI's makeprofiledb \
                 aborts there; -blastdb_version 4 builds it)",
                String::from_utf8_lossy(&name)
            )),
            WriteDbError::Creation(message) | WriteDbError::Options(message) => {
                RunError::Exception(message.into_bytes())
            }
            WriteDbError::CreationBytes(message) => RunError::Exception(message),
            WriteDbError::Exception(exception) => RunError::Exception(exception.get_msg().to_vec()),
            WriteDbError::Database(message) => RunError::SeqDb(message.into_bytes()),
            WriteDbError::OutOfMemory
            | WriteDbError::InvalidArgument(_)
            | WriteDbError::LengthError(_) => RunError::Unknown,
            WriteDbError::OrdinalIdWithoutIds => RunError::Refused(
                "a query whose first definition line has no Seq-id, with -index false (NCBI's \
                 makeprofiledb writes its gnl|BL_ORD_ID id past the end of an empty list)"
                    .to_string(),
            ),
            WriteDbError::Other(message) | WriteDbError::Unsupported(message) => {
                RunError::Refused(format!("the database error '{message}'"))
            }
            WriteDbError::Io { path, source } => {
                RunError::Failed(format!("{}: {source}", path.display()))
            }
            // Only makeblastdb's log raises it (`writedb::build_db::endl`); makeprofiledb's
            // writer never does.
            WriteDbError::ClosedPipe => unreachable!("makeblastdb's closed log in makeprofiledb"),
        }
    }
}

/// The parameters of `x_InitProgramParameters` that the build reads.
struct Params<'a> {
    /// The standard output (`cout`): the database writer's lines (`WriteDb::take_stdout`), and
    /// `m_LogFile` without a file of `-logfile`.
    stdout: &'a mut dyn Write,
    /// `m_LogFile` when `-logfile` names a file (`None`: the standard output).
    log_file: Option<OpenedFile>,
    /// The list of profile files (`m_InPssmList`), read whole.
    in_list: Vec<u8>,
    /// `-binary` (`m_binary_scoremat`): read profile files as ASN.1 binary (`MSerial_AsnBinary`).
    binary: bool,
    /// `-title` (its bytes; `-in` without it).
    title: Vec<u8>,
    word_default_score_threshold: f64,
    out_db_name: Vec<u8>,
    max_smp_files_per_vol: i32,
    create_index_file: bool,
    gap_open: i32,
    gap_extend: i32,
    pssm_scale_factor: f64,
    matrix: String,
    use_model_threshold: bool,
    db_ver: DbVersion,
    taxids: TaxIdSet,
    user_taxids: bool,
    /// The run's `SERIAL_VERIFY_DATA_READ` and `SERIAL_VERIFY_DATA_GET`
    /// (`ApplicationSettings::read_verification`), which the read of each profile checks
    /// (`read_profile`).
    read_verification: ReadVerification,
}

impl Params<'_> {
    /// `m_LogFile`: the file of `-logfile`, or the standard output.
    fn log(&mut self) -> &mut dyn Write {
        match &mut self.log_file {
            Some(file) => file,
            None => &mut *self.stdout,
        }
    }
}

/// What `Init` ends with instead of the parameters: an exit status, its report written.
type InitEnd = i32;

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:485-583
// An exception of `Init` is the application's: a `CArgException` of a file argument
// (`AsOutputFile`, `AsInputFile`) is reported as the reader's errors are, exit status 1; another
// (`SetMappingFromFile`'s conversions) with "Application's initialization failed", exit status 2.
/// `x_InitProgramParameters`.
fn init<'a>(
    args: &Args,
    input: &mut Cin<'_>,
    out: &'a mut dyn Write,
    diag: &mut Diag,
) -> Result<Params<'a>, InitEnd> {
    let log_file = if args.given("logfile") {
        match args.output_file("logfile", &mut *out) {
            Ok(OutputFile::File(file)) => Some(file),
            Ok(OutputFile::Stdout(_)) => None,
            // LOSAT Web's refusal of a name that it does not write (`OpenFiles::create`), which
            // the files of the command line and the library never give.
            Ok(OutputFile::Refused(reason)) => return Err(diag.refuse(&unsupported(&reason))),
            Err(error) => return Err(diag.init_arg_exception(error)),
        }
    } else {
        None
    };
    let in_name = args
        .os_value("in")
        .ok_or_else(|| crate::ncbi_app::args::no_value("in"))
        .map_err(|error| diag.init_arg_exception(error))?;
    let mut in_file = files::open_input_file("in", in_name, input)
        .map_err(|error| diag.init_arg_exception(error))?;
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:585-600
    // `x_GetSMPFilenames` reads the stream in `Run`; LOSAT reads it here. A stream whose reads
    // fail (a directory) never reaches its end, and NCBI's loop does not end.
    let mut in_list = Vec::new();
    if let Err(error) = in_file.read_to_end(&mut in_list) {
        return Err(diag.refuse(&unsupported(&format!(
            "an -in list that cannot be read (`{}': {error}; NCBI's makeprofiledb reads it \
             without end)",
            in_name.to_string_lossy()
        ))));
    }
    drop(in_file);
    let in_bytes = crate::utils::os_bytes::os_str_bytes(in_name);
    let binary = args.flag("binary");
    let title = args
        .bytes("title")
        .map_or_else(|| in_bytes.to_vec(), std::borrow::Cow::into_owned);
    let word_default_score_threshold = args.double("threshold").expect("a default");
    let out_db_name = args
        .bytes("out")
        .map_or_else(|| in_bytes.to_vec(), std::borrow::Cow::into_owned);
    let max_smp_files_per_vol = args.integer("max_smp_vol").expect("a default");
    let create_index_file = args.boolean("index").expect("a default");
    let matrix = args.value("matrix").expect("a default").to_string();
    let (default_gap_open, default_gap_extend) =
        crate::stats::protein_options::protein_gap_existence_extend_params(&matrix)
            .unwrap_or((0, 0));
    let gap_open = args.integer("gapopen").unwrap_or(default_gap_open);
    let gap_extend = args.integer("gapextend").unwrap_or(default_gap_extend);
    let pssm_scale_factor = args.double("scale").expect("a default");
    let use_model_threshold = !args.flag("force");
    let db_ver = match args.integer("blastdb_version") {
        Some(4) => DbVersion::V4,
        _ => DbVersion::V5,
    };
    let (taxids, user_taxids) = if let Some(taxid) = args.integer("taxid") {
        (TaxIdSet::with_global(taxid), true)
    } else if let Some(name) = args.os_value("taxid_map") {
        let file = files::open_input_file("taxid_map", name, input)
            .map_err(|error| diag.init_arg_exception(error))?;
        let mut taxids = TaxIdSet::default();
        let mut warnings = Vec::new();
        let read = taxids
            .set_mapping_from_file(&mut std::io::BufReader::new(file), &mut |warning| {
                warnings.push(warning)
            });
        for warning in warnings {
            post_seq_id_warning(diag, &warning);
        }
        match read {
            Ok(_) => {}
            Err(TaxIdMapError::Exception(exception)) => {
                diag.write_text(&exception.report_default_ex(
                    106,
                    15,
                    "Application's initialization failed",
                ));
                return Err(2);
            }
            Err(TaxIdMapError::BlankLine(line)) => {
                return Err(diag.refuse(&unsupported(&format!(
                    "a -taxid_map line of white space only (line {line} of `{}'; NCBI's \
                     makeprofiledb crashes on it)",
                    name.to_string_lossy()
                ))))
            }
        }
        (taxids, true)
    } else {
        (TaxIdSet::default(), false)
    };
    Ok(Params {
        stdout: out,
        log_file,
        in_list,
        binary,
        title,
        word_default_score_threshold,
        out_db_name,
        max_smp_files_per_vol,
        create_index_file,
        gap_open,
        gap_extend,
        pssm_scale_factor,
        matrix,
        use_model_threshold,
        db_ver,
        taxids,
        user_taxids,
        read_verification: ReadVerification::default(),
    })
}

/// A warning of the Seq-id parser (a warning that `CSeq_id` posts with `ERR_POST_X`), as
/// makeprofiledb posts it (no prefix of its own).
fn post_seq_id_warning(diag: &mut Diag, warning: &crate::seqid::accguide::SeqIdWarning) {
    let mut line = format!("Warning: (1306.{}) ", warning.subcode).into_bytes();
    line.extend_from_slice(&warning.message);
    diag.post_bytes(&line);
}

/// The state that makeprofiledb's destructor reads: the volumes made (`m_VolNames`), whether the
/// build ended (`m_Done`) and the output's name.
struct Built {
    vol_names: Vec<Vec<u8>>,
    done: bool,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1866-1904
/// makeprofiledb's `Run` after `Init` (`args` read): the build, its report, and the clean-up of
/// its destructor; the exit status.
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    input: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let mut params = match init(args, input, out, diag) {
        Ok(params) => params,
        Err(status) => return status,
    };
    params.read_verification = settings.read_verification.clone();
    if let Some(refusal) = refusal(args) {
        return diag.refuse(&refusal);
    }
    let mut built = Built {
        vol_names: Vec::new(),
        done: false,
    };
    let status = match x_run(&mut params, &mut built, settings, diag) {
        Ok(()) => 0,
        Err(error) => report(error, diag),
    };
    let _ = params.log().flush();
    destroy(&built, &params.out_db_name);
    status
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1304-1308
/// makeprofiledb with `-dryrun`: its `Init` (the parameters, the files of `-logfile`, `-in` and
/// `-taxid_map`) and no `Run`; the exit status.
pub(super) fn dry_run(
    args: &Args,
    _settings: &ApplicationSettings,
    input: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    match init(args, input, out, diag) {
        Ok(mut params) => {
            let _ = params.log().flush();
            0
        }
        Err(status) => status,
    }
}

/// LOSAT's refusals of a command line that NCBI's makeprofiledb reads: database types other than
/// `rps` (row PS14).
fn refusal(args: &Args) -> Option<String> {
    match args.value("dbtype") {
        Some("rps") | None => None,
        Some(dbtype) => Some(unsupported(&format!(
            "-dbtype {dbtype} (only -dbtype rps is built)"
        ))),
    }
}

/// The report of `Run`'s handler; the exit status.
fn report(error: RunError, diag: &mut Diag) -> i32 {
    let post = |diag: &mut Diag, prefix: &[u8], message: &[u8]| {
        diag.err_post_bytes(&[prefix, message].concat());
    };
    match error {
        RunError::Input(message) => {
            post(diag, b"INPUT ERROR: ", &message);
            1
        }
        RunError::SeqDb(message) => {
            post(diag, b"ERROR: ", &message);
            2
        }
        RunError::Blast(message) => {
            post(diag, b"ERROR: ", &message);
            1
        }
        RunError::Exception(message) => {
            post(diag, b"ERROR: ", &message);
            255
        }
        RunError::Unknown => {
            diag.err_post("Error: Unknown exception");
            255
        }
        RunError::Refused(what) => diag.refuse(&unsupported(&what)),
        RunError::Failed(text) => diag.refuse(&text),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:334-375
/// makeprofiledb's destructor: after a build that did not end, the RPS files of every volume
/// begun (and the alias file of several); after one that ended, the `.pog` file of each volume.
fn destroy(built: &Built, out_db_name: &[u8]) {
    let remove = |name: &[u8], extension: &str| {
        let path = [name, b".", extension.as_bytes()].concat();
        crate::writedb::dir_entry::remove(&path_of_bytes(&path));
    };
    if !built.done {
        for name in &built.vol_names {
            for extension in ["rps", "loo", "aux", "freq"] {
                remove(name, extension);
            }
        }
        if built.vol_names.len() > 1 {
            remove(out_db_name, "pal");
        }
    } else {
        for name in &built.vol_names {
            remove(name, "pog");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:680-703
/// `s_DeleteMakeprofileDb`: the RPS files of a database of the name (each extension of the
/// single volume, else of the volumes `.00`, `.01`, ... while they are removed), then the
/// database itself (`DeleteBlastDb`); whether the database was removed.
fn delete_makeprofile_db(name: &[u8]) -> Result<bool, RunError> {
    let named = |suffix: String| path_of_bytes(&[name, suffix.as_bytes()].concat()).into_owned();
    for extension in [
        ".rps", ".loo", ".aux", ".freq", ".blocks", ".wcounts", ".obsr",
    ] {
        if !crate::writedb::dir_entry::remove(&named(extension.to_string())) {
            let mut index = 0u32;
            while crate::writedb::dir_entry::remove(&named(format!(
                ".{}{}{extension}",
                index / 10,
                index % 10
            ))) {
                index += 1;
            }
        }
    }
    Ok(delete_blast_db(name, true)?)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:585-600
/// `x_GetSMPFilenames`: the names of the list, split at spaces, tabs and carriage returns.
fn smp_filenames(list: &[u8]) -> Result<Vec<&[u8]>, RunError> {
    let names: Vec<&[u8]> = list
        .split(|&c| matches!(c, b' ' | b'\t' | b'\r' | b'\n'))
        .filter(|name| !name.is_empty())
        .collect();
    if names.is_empty() {
        return Err(RunError::Input(
            b"Input file contains no smp filnames".to_vec(),
        ));
    }
    Ok(names)
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1228-1263
// `m_NumOfVols` is one more than the quotient of the count and the maximum, with C's division: a
// maximum of 0 divides by zero (NCBI 2.17.0: SIGFPE), and so does a count of volumes of 0
// (`num_seqs` is the count divided by `m_NumOfVols`); a negative count
// makes no volume, no alias file and exit status 0.
/// `x_Run`.
fn x_run(
    params: &mut Params<'_>,
    built: &mut Built,
    settings: &ApplicationSettings,
    diag: &mut Diag,
) -> Result<(), RunError> {
    create_directories(&params.out_db_name)?;
    if delete_makeprofile_db(&params.out_db_name)? {
        let _ = params
            .log()
            .write_all(b"Deleted existing BLAST database with identical name.\n");
        let _ = params.log().flush();
    }
    let list = std::mem::take(&mut params.in_list);
    let smps = smp_filenames(&list)?;
    let num_smps = i32::try_from(smps.len())
        .map_err(|_| RunError::Refused("a list of 2^31 profiles or more".to_string()))?;
    let max = params.max_smp_files_per_vol;
    let undefined = || {
        RunError::Refused(format!(
            "-max_smp_vol {max} for {num_smps} profiles (NCBI's makeprofiledb divides by zero)"
        ))
    };
    if max == 0 {
        return Err(undefined());
    }
    let num_of_vols = num_smps / max + 1;
    if num_of_vols == 0 {
        return Err(undefined());
    }
    let num_seqs = num_smps / num_of_vols;
    let mut residue_seqs = num_smps % num_of_vols;
    if num_of_vols == 1 {
        make_vol(params, built, diag, None, &smps)?;
        built.done = true;
        return Ok(());
    }
    let mut b = 0usize;
    let mut r = num_seqs.max(0) as usize;
    for i in 0..num_of_vols.max(0) {
        make_vol(
            params,
            built,
            diag,
            Some(i as u32),
            &smps[b..r.min(smps.len())],
        )?;
        b = r;
        r = b + num_seqs as usize;
        if residue_seqs > 0 {
            r += 1;
            residue_seqs -= 1;
        }
    }
    if usize::try_from(num_of_vols).is_ok_and(|vols| vols == built.vol_names.len()) {
        create_alias_file(params, built, settings, diag)?;
        built.done = true;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1865-1876
// The volumes are named by their file names (`SeqDB_RemoveDirName`), without a list
// (`eNoAliasFilterType`); the reader opens the new alias file on the database search path (the
// current directory, `BLASTDB`, the registry's `[BLAST] BLASTDB`). The names and the title are
// their bytes.
/// `x_CreateAliasFile`: the alias file of the volumes, its totals, and the log line
/// (`writedb::alias_file::create_alias_file`).
fn create_alias_file(
    params: &Params<'_>,
    built: &Built,
    settings: &ApplicationSettings,
    diag: &mut Diag,
) -> Result<(), RunError> {
    let names: Vec<&[u8]> = built
        .vol_names
        .iter()
        .map(|name| crate::seqdb::alias::base_name(name))
        .collect();
    let path = crate::seqdb::SearchPath::from_environment(settings.blastdb.as_deref());
    let created = alias_file::create_alias_file(
        &params.out_db_name,
        &names,
        crate::seqdb::SeqType::Protein,
        b"",
        &params.title,
        FilterType::None,
        &path,
    );
    match created {
        Ok(line) => {
            diag.post_bytes(&line);
            Ok(())
        }
        Err(AliasError::SeqDb(text)) => Err(RunError::SeqDb(text)),
        Err(AliasError::Reader(error)) => Err(RunError::Failed(
            crate::algorithm::common::database::losat_reader_error(&error, super::PROGRAM),
        )),
        Err(AliasError::Clock) => Err(RunError::Failed(CLOCK_ERROR.to_string())),
    }
}

/// Posts the messages emitted by the shared ASN.1 reader (`Post`) as makeprofiledb posts them:
/// makeprofiledb leaves NCBI's diagnostic post level at its default (`eDiag_Error`), so warnings
/// (`Post::Utf8Tag`, 803.10) are suppressed while errors (`Post::BadChar`, `Post::Error`,
/// `Post::UnendedByteBlock`) are written with no application prefix.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1381-1395
fn post_asn_posts(diag: &mut Diag, posts: &[Post]) {
    for post in posts {
        match post {
            Post::BadChar(replaced) => diag.post_replaced(std::slice::from_ref(replaced)),
            Post::Error(text) => diag.err_post(text),
            Post::Utf8Tag(_) => {}
            other => diag.write_text(&other.line(None)),
        }
    }
}

/// Maps an ASN.1 [`ReadError`] of root type `root_type` to makeprofiledb's `Run` error: a
/// `CSerialException` or `CException`'s newest `GetMsg()`, the root object's `Frame type=` for
/// `CEofException`, `Unknown` for `std::exception`, and `Refused` for LOSAT's reader limits or
/// a check turned off by `read_verification` (`ReadVerification::refused`).
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1878-1904
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
fn read_error_to_run_error(
    error: ReadError,
    root_type: &str,
    input: &str,
    read_verification: &ReadVerification,
    refuse: impl FnOnce(String) -> String,
) -> RunError {
    match error {
        ReadError::Serial(exception) => match read_verification.refused(&exception, input) {
            Some(refused) => RunError::Refused(refused),
            None => RunError::Exception(exception.get_msg().to_vec()),
        },
        ReadError::Exception(exception) => RunError::Exception(exception.get_msg().to_vec()),
        ReadError::Eof => RunError::Exception(
            format!(" Frame type= eFrameNamed, Object type= {root_type}").into_bytes(),
        ),
        ReadError::Std(_) | ReadError::OutOfMemory => RunError::Unknown,
        ReadError::Unsupported(what) => RunError::Refused(refuse(what)),
    }
}

/// A profile file read: the reader's posts written, its errors as makeprofiledb's handler
/// reports them.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1381-1395
// The read of `pssm` from `in_stream` in ASN.1 text (`MSerial_AsnText`) or binary
// (`MSerial_AsnBinary` when `-binary` is set): a `CSerialException`'s message (`line N: ...` or
// `byte N: ...`), the frame of the root object for an input that ends in it (`CEofException`),
// the message of the newest exception of the reader's stack; a `std::exception` is the handler's
// unknown one.
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
// A mandatory member that the profile skips is NCBI's read error unless `SERIAL_VERIFY_DATA_READ`
// turns the check off, which LOSAT refuses there (`read_verification`, `ReadVerification::refused`).
fn read_profile(
    path: &Path,
    filename: &[u8],
    binary: bool,
    read_verification: &ReadVerification,
    diag: &mut Diag,
) -> Result<PssmWithParameters, RunError> {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1386-1395
    // The file is read through a `CNcbiIfstream` whose open is not checked: a file that does not
    // open (no permission) reads as an empty stream, whose reader stops at once.
    let bytes = std::fs::read(path).unwrap_or_default();
    let read = crate::blastinput::pssm_text::read_profile(&bytes, binary);
    post_asn_posts(diag, &read.posts);
    use crate::blastinput::pssm_text::PssmReadError;
    let refuse = |what: String| {
        format!(
            "a profile that LOSAT's PSSM reader does not read ({}: {what})",
            String::from_utf8_lossy(filename)
        )
    };
    match read.result {
        Ok(pssm) => Ok(pssm),
        Err(PssmReadError::Ncbi(error)) => Err(read_error_to_run_error(
            error,
            "PssmWithParameters",
            "a profile (-in)",
            read_verification,
            refuse,
        )),
        Err(PssmReadError::Refused(what)) => Err(RunError::Refused(refuse(what))),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:464-466
// The writer's line on `cout` is written with `endl`: at once, before what the run writes after
// it.
/// The writer's diagnostics since the last call, posted, and its text on the standard output
/// written there.
fn post_writer(diag: &mut Diag, stdout: &mut dyn Write, output_db: &mut WriteDb) {
    let text = output_db.take_diagnostics();
    if !text.is_empty() {
        diag.write_text(&text);
    }
    let text = output_db.take_stdout();
    if !text.is_empty() {
        let _ = stdout.write_all(&text);
        let _ = stdout.flush();
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1366-1492
// On an error the volume's writer is destroyed, which closes it (`~CWriteDB`): the profiles
// added stay in the database's files.
/// `x_MakeVol`: the volume `vol` (`None`: the only one) of the profiles `smps`.
fn make_vol(
    params: &mut Params<'_>,
    built: &mut Built,
    diag: &mut Diag,
    vol: Option<u32>,
    smps: &[&[u8]],
) -> Result<(), RunError> {
    let db_name = match vol {
        None => params.out_db_name.clone(),
        Some(index) => make_short_name(&params.out_db_name, index),
    };
    let volume_params = VolumeParams {
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        matrix: params.matrix.as_bytes().to_vec(),
        scale_factor: rps_files::int4_of(params.pssm_scale_factor.ceil()).ok_or_else(|| {
            RunError::Refused(format!(
                "-scale {} (NCBI's conversion of its ceiling to a 32-bit integer has no defined \
                 result)",
                params.pssm_scale_factor
            ))
        })?,
    };
    let mut rps = rps_files::open_volume(&db_name, smps.len(), volume_params)?;
    built.vol_names.push(db_name.clone());
    let now = LocalTime::now().ok_or_else(|| {
        RunError::Failed("the clock cannot be read for the database's date".to_string())
    })?;
    let indices = if params.create_index_file {
        Indices {
            full: true,
            hash: false,
            trace: true,
        }
    } else {
        Indices::default()
    };
    let mut output_db = WriteDb::new(
        &db_name,
        true,
        &params.title,
        now.index_date().as_bytes(),
        MAX_FILE_SIZE,
        params.db_ver,
        Platform::host(),
        indices,
    );
    output_db.set_diag_prefix(None);
    output_db.set_max_file_size(MAX_FILE_SIZE);
    let result = add_profiles(params, diag, &mut rps, &mut output_db, smps);
    // `output_db->Close()` after the loop, or the writer's destructor in the unwinding of an
    // error (`~CWriteDB_Impl` closes it), then the destructors of the writer's members.
    let closed = output_db.close();
    if let Err(RunError::Refused(_)) = result {
        // LOSAT's refusal: NCBI's run has no defined result there, and nothing after it is
        // reported.
        drop(output_db.take_diagnostics());
        drop(output_db.take_stdout());
        return result;
    }
    post_writer(diag, params.stdout, &mut output_db);
    match (result, closed) {
        (Ok(()), Ok(())) => {}
        // The `Close` after the loop threw: its error leaves `x_MakeVol`, the destructor's
        // `Close` returns at once and the members close (the last volume under its numbered
        // name, the LMDB files).
        (Ok(()), Err(close_error)) => {
            let destroyed = output_db.destroy_after_failed_close();
            post_writer(diag, params.stdout, &mut output_db);
            let error = RunError::from(close_error);
            return Err(match destroyed {
                Ok(()) => error,
                Err(destructor) => aborts(error, &destructor),
            });
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:100-108
        // In the unwinding the destructor posts a `CWriteDBException` of its `Close` and the first
        // error goes on; any other exception leaves the destructor and ends NCBI's process
        // (`std::terminate`), which LOSAT refuses.
        (Err(error), Err(close_error)) => {
            if !matches!(
                close_error,
                WriteDbError::Creation(_) | WriteDbError::CreationBytes(_)
            ) {
                return Err(aborts(error, &close_error));
            }
            diag.err_post_bytes(&close_error.message());
            let destroyed = output_db.destroy_after_failed_close();
            post_writer(diag, params.stdout, &mut output_db);
            return Err(match destroyed {
                Ok(()) => error,
                Err(destructor) => aborts(error, &destructor),
            });
        }
        (Err(error), Ok(())) => {
            let destroyed = output_db.destroy_after_close();
            post_writer(diag, params.stdout, &mut output_db);
            return Err(match destroyed {
                Ok(()) => error,
                Err(destructor) => aborts(error, &destructor),
            });
        }
    }
    rps.close()?;
    Ok(())
}

/// LOSAT's refusal where the writer, closed or destroyed after the error `first`, throws
/// `then` out of a destructor, where NCBI's process ends (`std::terminate`, SIGABRT); `first`
/// itself when it is LOSAT's refusal (NCBI's run left its defined behaviour there).
fn aborts(first: RunError, then: &WriteDbError) -> RunError {
    if let RunError::Refused(_) = first {
        return first;
    }
    let message = then.message();
    let message = message.strip_suffix(b"\n").unwrap_or(&message);
    RunError::Refused(format!(
        "an error ({}) after which the database's writer stops with another ({}; NCBI's \
         makeprofiledb aborts there)",
        first.text(),
        String::from_utf8_lossy(message)
    ))
}

/// The profiles of a volume added to its writer and its RPS files (the loop of `x_MakeVol`).
fn add_profiles(
    params: &mut Params<'_>,
    diag: &mut Diag,
    rps: &mut RpsVolume,
    output_db: &mut WriteDb,
    smps: &[&[u8]],
) -> Result<(), RunError> {
    for &filename in smps {
        let path = path_from_bytes(filename);
        if !std::fs::metadata(&path).is_ok_and(|meta| meta.is_file()) {
            return Err(RunError::Input([filename, b" does not exists"].concat()));
        }
        let mut pssm_w_parameters = read_profile(
            &path,
            filename,
            params.binary,
            &params.read_verification,
            diag,
        )?;
        let sm = check_input_scoremat(&pssm_w_parameters, filename)?;
        // `seq_size`, the `int` of `GetQueryLength()`: negative past `INT_MAX` (`query_length`).
        let seq_size = query_length(&pssm_w_parameters);
        add_to_database(params, diag, output_db, &pssm_w_parameters)?;
        if !rps.has_lookup() {
            update_rps_db_info(&mut rps.params, &pssm_w_parameters);
            fill_in_rps_db_parameters(&rps.params, &pssm_w_parameters)?;
            let mut word_score_threshold = params.word_default_score_threshold;
            if sm == Scoremat::HasPssm {
                let final_data = pssm_w_parameters
                    .pssm
                    .final_data
                    .as_ref()
                    .expect("has a PSSM");
                rps.params.scale_factor = final_data.scaling_factor.unwrap_or(1);
                if params.use_model_threshold {
                    if let Some(threshold) = final_data.word_score_threshold {
                        word_score_threshold = threshold;
                    }
                }
            } else {
                update_statistics(&rps.params, &mut pssm_w_parameters, seq_size)?;
            }
            rps.add_first(f64::from(rps.params.scale_factor) * word_score_threshold)?;
        } else {
            fill_in_rps_db_parameters(&rps.params, &pssm_w_parameters)?;
            if sm == Scoremat::FreqOnly {
                update_statistics(&rps.params, &mut pssm_w_parameters, seq_size)?;
            }
            let final_data = pssm_w_parameters
                .pssm
                .final_data
                .as_ref()
                .expect("final data");
            if final_data.scaling_factor.unwrap_or(1) != rps.params.scale_factor {
                return Err(RunError::Blast(b"Scaling factors do not match".to_vec()));
            }
            let threshold = match final_data.word_score_threshold {
                Some(threshold) if params.use_model_threshold => threshold,
                _ => params.word_default_score_threshold,
            };
            rps.set_threshold(f64::from(rps.params.scale_factor) * threshold)?;
        }
        let pssm = &pssm_w_parameters.pssm;
        let final_data = pssm.final_data.as_ref().expect("final data");
        let by_row = pssm.by_row.unwrap_or(false);
        let Ok(seq_size) = usize::try_from(seq_size) else {
            return Err(
                rps_files::negative_length_error(seq_size, by_row, &final_data.scores).into(),
            );
        };
        let freq_ratios = pssm
            .intermediate_data
            .as_ref()
            .map_or(&[][..], |data| &data.freq_ratios[..]);
        rps.add_profile(&ProfileData {
            seq_size,
            num_rows: pssm.num_rows as usize,
            by_row,
            scores: &final_data.scores,
            freq_ratios,
            kappa: final_data.kappa,
        })?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1310-1317,1406-1437
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:865-874
// A query with descriptors gets the definition lines of `ExtractBioseqDeflines` (its binary
// ones, else one line of its ids and title: none without ids, "No deflines provided."), one
// without them a line of its ids. Without `-taxid` and `-taxid_map` the taxid of the query's
// first `org` descriptor maps its first id (`AddTaxId`; a query without ids maps none), and
// `FixTaxId` then gives each definition line the map's taxid of its ids (0 for none). With
// `-index false` the writer puts `gnl|BL_ORD_ID|<oid>` in place of the first id of the first
// definition line: NCBI writes it past the end of an empty list (undefined behaviour), which
// LOSAT refuses. With `-index true` a sequence without ids stops at the writer's duplicate check.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiobj.cpp:1001-1012
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1878-1901
/// A query descriptor stated as `NULL` that a loop over the descriptors dereferences
/// (`BioseqDeflinesError::NullDescriptor`): the `CCoreException` of a NULL `CRef`, which `Run`'s
/// handler reports as a `CException` (`Error: ERROR: Attempt to access NULL pointer.`, exit status
/// 255); the volume begun is closed as at any error (makeprofiledb 2.17.0, probe D3c-3c `probeB`:
/// nothing for the first profile, the closed volume of the profiles before for a later one). The
/// message is the one owner's (`ncbi_exception::NULL_POINTER`).
fn null_descriptor() -> RunError {
    RunError::Exception(crate::ncbi_exception::NULL_POINTER.as_bytes().to_vec())
}

/// The profile's query added to the database: its definition lines (the query's descriptors,
/// else one line of its ids), their taxids, the sequence.
fn add_to_database(
    params: &mut Params<'_>,
    diag: &mut Diag,
    output_db: &mut WriteDb,
    pssm_w_parameters: &PssmWithParameters,
) -> Result<(), RunError> {
    let Some(SeqEntry::Seq(query)) = &pssm_w_parameters.pssm.query else {
        unreachable!("check_input_scoremat found the query's Bioseq");
    };
    let (mut deflines, null_seqids) = match &query.descr {
        Some(descr) => {
            let read = bioseq_deflines(&query.ids, descr);
            post_asn_posts(diag, &read.posts);
            let extracted = read.result.map_err(|error| match error {
                BioseqDeflinesError::NoDeflines => {
                    RunError::Exception(b"Error: No deflines provided.".to_vec())
                }
                BioseqDeflinesError::BinaryHeader(error) => read_error_to_run_error(
                    error,
                    "Blast-def-line-set",
                    "binary definition lines (ASN1_BlastDefLine)",
                    &params.read_verification,
                    |what| {
                        format!(
                            "binary definition lines (a user object ASN1_BlastDefLine) that \
                             LOSAT's reader does not read ({what})"
                        )
                    },
                ),
                BioseqDeflinesError::NullDescriptor => null_descriptor(),
            })?;
            // The binary definition lines are not encoded again here (`bin_hdr` is kept).
            if !extracted.binary {
                output_db.post_extracted_deflines(&extracted.set);
                post_writer(diag, params.stdout, output_db);
            }
            (extracted.set, extracted.null_seqids)
        }
        None => (
            DeflineSet {
                deflines: vec![BlastDefline {
                    ids: query.ids.clone(),
                    ..BlastDefline::default()
                }],
                utf8: false,
            },
            Vec::new(),
        ),
    };
    if !params.user_taxids {
        let first_id = query.ids.first();
        let taxid = match query.descr.as_deref().map(first_org_taxid) {
            Some(Err(_)) => return Err(null_descriptor()),
            Some(Ok(taxid)) => taxid,
            None => None,
        };
        if let (Some(taxid), Some(first_id)) = (taxid, first_id) {
            let mut warnings = Vec::new();
            params
                .taxids
                .add_tax_id(first_id, taxid, &mut |warning| warnings.push(warning))
                .map_err(|exception| RunError::Exception(exception.get_msg().to_vec()))?;
            for warning in warnings {
                post_seq_id_warning(diag, &warning);
            }
        }
    }
    if params
        .taxids
        .fix_tax_id_checked(&mut deflines, &null_seqids)
        .is_err()
    {
        return Err(null_descriptor());
    }
    let inst = &query.inst;
    let seq_data = inst.seq_data.as_ref().map(|data| match data {
        SeqData::Iupacna(letters) => bioseq::SeqData::Iupacna(letters.clone()),
        SeqData::Iupacaa(letters) => bioseq::SeqData::Iupacaa(letters.clone()),
        SeqData::Ncbi2na(bases) => bioseq::SeqData::Ncbi2na(bases.clone()),
        SeqData::Ncbi4na(bases) => bioseq::SeqData::Ncbi4na(bases.clone()),
        SeqData::Ncbieaa(letters) => bioseq::SeqData::Ncbieaa(letters.clone()),
        SeqData::Ncbistdaa(codes) => bioseq::SeqData::Ncbistdaa(codes.clone()),
        SeqData::Other(_) => bioseq::SeqData::Other,
    });
    let length = inst.length.unwrap_or(0);
    let iupacna_hi = seq_data.as_ref().and_then(iupacna_high_byte);
    let short = seq_data
        .as_ref()
        .and_then(|data| short_nucleotide_data(data, length));
    let bioseq = Bioseq {
        local_id: String::new(),
        title: Vec::new(),
        is_aa: inst.mol.1 == 3,
        length: u64::from(length),
        seq_data,
        parsed: Some(ParsedIds {
            ids: query.ids.clone(),
            hist_old_id: None,
            defline_field: Vec::new(),
            title_descr: None,
        }),
        delta_fix: None,
    };
    let added = output_db.add_sequence(bioseq);
    post_writer(diag, params.stdout, output_db);
    added?;
    output_db.set_deflines_with_null_seqids(deflines, null_seqids);
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:590-604
    // NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
    // Before `WriteDB_Ncbi4naToBinary`, `WriteDB_IupacnaToBinary` runs `x_ConvertIupacnaTo4na` on
    // all `v.size()` bytes; its pair loop indexes the conversion table with a signed `char`, so
    // any byte >= 0x80 in `0..(v.len() & !1)` reads before the table (undefined behaviour).
    if let Some((index, byte)) = iupacna_hi {
        return Err(RunError::Refused(format!(
            "a profile whose query has an iupacna byte 0x{byte:02X} at position {index} (NCBI's \
             makeprofiledb indexes its iupacna-to-ncbi4na table with a signed char and reads \
             before the table)"
        )));
    }
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:158-180
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:100-108,324-331
    // A sequence that the writer keeps (its molecule checked) is published at the next
    // `AddSequence`, at `Close`, or at the destructor's `Close` after an error: NCBI's
    // conversion of short nucleotide data then reads past it (`short_nucleotide_data`), whatever
    // the profile and whatever follows.
    if let Some(short) = short {
        let ShortNucleotideData {
            coding,
            text,
            held,
            read,
        } = short;
        let count = |n: usize, unit: &str| match n {
            1 => format!("1 {unit}"),
            _ => format!("{n} {unit}s"),
        };
        let (held, read) = if text {
            (
                count(held, "letter"),
                format!("{} of their ncbi4na form", count(read, "byte")),
            )
        } else {
            (count(held, "byte"), count(read, "byte"))
        };
        return Err(RunError::Refused(format!(
            "a profile whose query has {held} of {coding} data for a length of {length} (NCBI's \
             makeprofiledb reads {read} there, past the data)"
        )));
    }
    Ok(())
}

#[cfg(test)]
#[path = "work_tests.rs"]
mod tests;
