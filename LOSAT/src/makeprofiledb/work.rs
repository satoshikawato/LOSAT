//! makeprofiledb's `Init` (`x_InitProgramParameters`), `Run` (`x_Run`, `x_MakeVol`, the alias
//! file of several volumes) and the clean-up of its destructor, for `-dbtype rps`.

use std::io::{Read, Write};
use std::path::Path;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{self, Cin, OutputFile};
use crate::ncbi_app::Diag;
use crate::report::serial::asn_types::{Post, ReadError};
use crate::report::serial::pssm::{PssmWithParameters, SeqData, SeqEntry};
use crate::seqid::{BlastDefline, DeflineSet};
use crate::utils::os_bytes::path_from_bytes;
use crate::writedb::bioseq::{self, Bioseq, ParsedIds};
use crate::writedb::build_db::{create_directories, delete_blast_db};
use crate::writedb::deflines::{bioseq_deflines, first_org_taxid, BioseqDeflinesError};
use crate::writedb::files::{make_short_name, DbVersion};
use crate::writedb::taxid_set::{TaxIdMapError, TaxIdSet};
use crate::writedb::time::LocalTime;
use crate::writedb::volume::Indices;
use crate::writedb::writer::WriteDb;
use crate::writedb::{Platform, WriteDbError};

use super::profile::{
    check_input_scoremat, fill_in_rps_db_parameters, update_rps_db_info, update_statistics,
    ProfileError, Scoremat,
};
use super::rps_files::{self, ProfileData, RpsError, RpsVolume, VolumeParams};

/// The limit of each file of a volume's database (`SetMaxFileSize(4000000000)`).
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
    /// LOSAT's refusal (a line on the standard error stream, exit status 1).
    Refused(String),
}

impl From<RpsError> for RunError {
    fn from(error: RpsError) -> Self {
        match error {
            RpsError::SeqDb(message) => RunError::SeqDb(message.into()),
            RpsError::Blast(message) => RunError::Blast(message.into()),
            RpsError::Unknown => RunError::Unknown,
            RpsError::Refused(reason) => RunError::Refused(reason),
            RpsError::Io(error) => RunError::Refused(format!(
                "a file of the RPS database could not be written ({error}); NCBI makeprofiledb \
                 does not test its writes and leaves the file short"
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
            ProfileError::Exception(message) => RunError::Exception(message.into_bytes()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:69-80 (the
// exceptions of CWriteDB: `CWriteDBException`, a `CException`; `CSeqDBException` of LMDB)
/// The writer's errors as makeprofiledb's handler reports them: its exceptions are
/// `CException`s ("ERROR: " and the message, 255), those of LMDB's environment
/// `CSeqDBException`s (2).
impl From<WriteDbError> for RunError {
    fn from(error: WriteDbError) -> Self {
        // A GI-only query between queries with other ids in a version 5 database: NCBI's writer
        // aborts at its close (survey I01g, I01x), as makeblastdb's does; LOSAT's refusal is
        // makeblastdb's (makedb-v5 DECISIONS 10).
        let gi_only = crate::writedb::writer::is_gi_only_between_keys(&error);
        match error {
            WriteDbError::CreationBytes(message) if gi_only => {
                RunError::Refused(String::from_utf8_lossy(&message).into_owned())
            }
            WriteDbError::Creation(message) | WriteDbError::Options(message) => {
                RunError::Exception(message.into_bytes())
            }
            WriteDbError::CreationBytes(message) => RunError::Exception(message),
            WriteDbError::Exception(exception) => RunError::Exception(exception.get_msg().to_vec()),
            WriteDbError::Database(message) => RunError::SeqDb(message.into_bytes()),
            WriteDbError::OutOfMemory => RunError::Unknown,
            WriteDbError::Other(message) | WriteDbError::Unsupported(message) => {
                RunError::Refused(message)
            }
            WriteDbError::Io { path, source } => {
                RunError::Refused(format!("{}: {source}", path.display()))
            }
        }
    }
}

/// The parameters of `x_InitProgramParameters` that the build reads.
struct Params<'a> {
    /// `m_LogFile`: the file of `-logfile`, or the standard output.
    log: OutputFile<'a>,
    /// The list of profile files (`m_InPssmList`), read whole.
    in_list: Vec<u8>,
    title: String,
    word_default_score_threshold: f64,
    out_db_name: String,
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
    let log = if args.given("logfile") {
        args.output_file("logfile", out)
            .map_err(|error| diag.init_arg_exception(error))?
    } else {
        OutputFile::Stdout(out)
    };
    let in_name = args
        .as_string("in")
        .map_err(|error| diag.init_arg_exception(error))?;
    let mut in_file = files::open_input_file("in", in_name, input)
        .map_err(|error| diag.init_arg_exception(error))?;
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:585-600
    // `x_GetSMPFilenames` reads the stream in `Run`; LOSAT reads it here. A stream whose reads
    // fail (a directory) never reaches its end, and NCBI's loop does not end.
    let mut in_list = Vec::new();
    if let Err(error) = in_file.read_to_end(&mut in_list) {
        return Err(diag.refuse(&format!(
            "makeprofiledb: the list of profiles `{in_name}' cannot be read ({error}); NCBI \
             makeprofiledb 2.17.0 reads such a stream without end"
        )));
    }
    drop(in_file);
    let in_string = in_name.to_string();
    let title = args
        .value("title")
        .map_or_else(|| in_string.clone(), str::to_string);
    let word_default_score_threshold = args.double("threshold").expect("a default");
    let out_db_name = args
        .value("out")
        .map_or_else(|| in_string.clone(), str::to_string);
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
    } else if let Some(name) = args.value("taxid_map") {
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
                return Err(diag.refuse(&format!(
                    "makeprofiledb: line {line} of -taxid_map `{name}' holds white space only, \
                     on which NCBI makeprofiledb 2.17.0 crashes"
                )))
            }
        }
        (taxids, true)
    } else {
        (TaxIdSet::default(), false)
    };
    Ok(Params {
        log,
        in_list,
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
    })
}

/// A warning of the Seq-id parser (`CSeq_id`'s `ERR_POST_X(Warning ...)`), as makeprofiledb
/// posts it (no prefix of its own).
fn post_seq_id_warning(diag: &mut Diag, warning: &crate::seqid::accguide::SeqIdWarning) {
    let mut line = format!("Warning: (1306.{}) ", warning.subcode).into_bytes();
    line.extend_from_slice(&warning.message);
    diag.post_bytes(&line);
}

/// The state that makeprofiledb's destructor reads: the volumes made (`m_VolNames`), whether the
/// build ended (`m_Done`) and the output's name.
struct Built {
    vol_names: Vec<String>,
    done: bool,
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1866-1904
/// makeprofiledb's `Run` after `Init` (`args` read): the build, its report, and the clean-up of
/// its destructor; the exit status.
pub(super) fn run(
    args: &Args,
    _settings: &ApplicationSettings,
    input: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let mut params = match init(args, input, out, diag) {
        Ok(params) => params,
        Err(status) => return status,
    };
    if let Some(reason) = refusal(args, &params) {
        return diag.refuse(&reason);
    }
    let mut built = Built {
        vol_names: Vec::new(),
        done: false,
    };
    let status = match x_run(&mut params, &mut built, diag) {
        Ok(()) => 0,
        Err(error) => report(error, diag),
    };
    let _ = params.log.flush();
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
            let _ = params.log.flush();
            0
        }
        Err(status) => status,
    }
}

/// LOSAT's refusals of a command line that NCBI's makeprofiledb reads: the binary profiles of
/// `-binary` (their reader, row PS10b), the database types other than `rps` (row PS14), and an
/// empty database name.
fn refusal(args: &Args, params: &Params<'_>) -> Option<String> {
    if args.flag("binary") {
        return Some(
            "makeprofiledb -binary (profiles in ASN.1 binary) is not supported by LOSAT's \
             makeprofiledb yet; give the profiles in ASN.1 text"
                .to_string(),
        );
    }
    match args.value("dbtype") {
        Some("rps") | None => {}
        Some(dbtype) => {
            return Some(format!(
                "makeprofiledb -dbtype {dbtype} is not supported by LOSAT's makeprofiledb yet; \
                 it builds -dbtype rps"
            ))
        }
    }
    if params.out_db_name.is_empty() {
        return Some(
            "makeprofiledb with an empty database name, with which NCBI makeprofiledb 2.17.0 \
             stops with an unknown exception or aborts"
                .to_string(),
        );
    }
    None
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
        RunError::Refused(reason) => diag.refuse(&format!("makeprofiledb: {reason}")),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:334-375
/// makeprofiledb's destructor: after a build that did not end, the RPS files of every volume
/// begun (and the alias file of several); after one that ended, the `.pog` file of each volume.
fn destroy(built: &Built, out_db_name: &str) {
    let remove = |name: String| {
        crate::writedb::dir_entry::remove(Path::new(&name));
    };
    if !built.done {
        for name in &built.vol_names {
            for extension in ["rps", "loo", "aux", "freq"] {
                remove(format!("{name}.{extension}"));
            }
        }
        if built.vol_names.len() > 1 {
            remove(format!("{out_db_name}.pal"));
        }
    } else {
        for name in &built.vol_names {
            remove(format!("{name}.pog"));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:680-703
/// `s_DeleteMakeprofileDb`: the RPS files of a database of the name (each extension of the
/// single volume, else of the volumes `.00`, `.01`, ... while they are removed), then the
/// database itself (`DeleteBlastDb`); whether the database was removed.
fn delete_makeprofile_db(name: &str) -> Result<bool, RunError> {
    for extension in [
        ".rps", ".loo", ".aux", ".freq", ".blocks", ".wcounts", ".obsr",
    ] {
        if !crate::writedb::dir_entry::remove(Path::new(&format!("{name}{extension}"))) {
            let mut index = 0u32;
            while crate::writedb::dir_entry::remove(Path::new(&format!(
                "{name}.{}{}{extension}",
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
// `m_NumOfVols = n / max + 1` with C's division: a maximum of 0 divides by zero (NCBI 2.17.0:
// SIGFPE), and so does a count of volumes of 0 (`num_seqs = n / m_NumOfVols`); a negative count
// makes no volume, no alias file and exit status 0.
/// `x_Run`.
fn x_run(params: &mut Params<'_>, built: &mut Built, diag: &mut Diag) -> Result<(), RunError> {
    create_directories(&params.out_db_name)?;
    if delete_makeprofile_db(&params.out_db_name)? {
        let _ = params
            .log
            .write_all(b"Deleted existing BLAST database with identical name.\n");
        let _ = params.log.flush();
    }
    let list = std::mem::take(&mut params.in_list);
    let smps = smp_filenames(&list)?;
    let num_smps = i32::try_from(smps.len())
        .map_err(|_| RunError::Refused("a list of 2^31 profiles or more".to_string()))?;
    let max = params.max_smp_files_per_vol;
    let undefined = || {
        RunError::Refused(format!(
            "-max_smp_vol {max} for {num_smps} profiles, with which NCBI makeprofiledb 2.17.0 \
             divides by zero"
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
        create_alias_file(params, built, diag)?;
        built.done = true;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1865-1876
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:284-303,305-328,358-414
// The totals are those of NCBI's reader opening the new alias file (`CSeqDB`): a database name
// that the reader splits (a space) does not open, and the alias file is removed.
/// `x_CreateAliasFile`: the alias file of the volumes, its totals, and the log line.
fn create_alias_file(params: &Params<'_>, built: &Built, diag: &mut Diag) -> Result<(), RunError> {
    let now = LocalTime::now().ok_or_else(|| {
        RunError::Refused("the clock cannot be read for the alias file's date".to_string())
    })?;
    let fname = format!("{}.pal", params.out_db_name);
    let mut text = format!("#\n# Alias file created {}\n#\n", now.log_time()).into_bytes();
    if !params.title.is_empty() {
        text.extend_from_slice(format!("TITLE {}\n", params.title).as_bytes());
    }
    text.extend_from_slice(b"DBLIST ");
    for name in &built.vol_names {
        // `SeqDB_RemoveDirName`: the name after the platform's last path separator.
        let base = name
            .rsplit(std::path::MAIN_SEPARATOR)
            .next()
            .unwrap_or(name);
        text.extend_from_slice(format!("\"{base}\" ").as_bytes());
    }
    text.extend_from_slice(b"\n");
    let _ = std::fs::write(&fname, &text);
    let totals = std::env::current_dir()
        .ok()
        .and_then(|cwd| crate::seqdb::alias::SearchPath::new(&cwd, None, None).ok())
        .and_then(|search| {
            crate::seqdb::SeqDb::open_in(
                &params.out_db_name,
                crate::seqdb::SeqType::Protein,
                &search,
            )
            .ok()
        })
        .map(|db| (db.num_seqs(), db.total_length()));
    let Some((num_seqs, length)) = totals else {
        crate::writedb::dir_entry::remove(Path::new(&fname));
        return Err(RunError::SeqDb(
            b"BLASTDB alias file creation failed.  Some referenced files may be missing".to_vec(),
        ));
    };
    if num_seqs == 0 {
        crate::writedb::dir_entry::remove(Path::new(&fname));
        return Err(RunError::SeqDb(
            b"No seqs in  were found in BLAST database".to_vec(),
        ));
    }
    if let Ok(mut file) = std::fs::OpenOptions::new().append(true).open(&fname) {
        let _ = write!(file, "NSEQ {num_seqs}\nLENGTH {length}\n");
    }
    diag.post(&format!(
        "Created protein BLAST (alias) database {} with {num_seqs} sequences",
        params.out_db_name
    ));
    Ok(())
}

/// A profile file read: the reader's posts written, its errors as makeprofiledb's handler
/// reports them.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1381-1393
// `in_stream >> MSerial_AsnText >> pssm`: a `CSerialException`'s message (`line N: ...`), the
// frame of the root object for a text that ends in it (`CEofException`), the message of the
// newest exception of the reader's stack; a `std::exception` is the handler's unknown one.
fn read_profile(
    path: &Path,
    filename: &[u8],
    diag: &mut Diag,
) -> Result<PssmWithParameters, RunError> {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1390-1394
    // The file is read through a `CNcbiIfstream` whose open is not checked: a file that does not
    // open (no permission) reads as an empty stream, whose reader stops at once.
    let text = std::fs::read(path).unwrap_or_default();
    let read = crate::blastinput::pssm_text::read_pssm(&text);
    for post in &read.posts {
        match post {
            Post::BadChar(replaced) => diag.post_replaced(std::slice::from_ref(replaced)),
            Post::Error(text) => diag.err_post(text),
        }
    }
    use crate::blastinput::pssm_text::PssmReadError;
    match read.result {
        Ok(pssm) => Ok(pssm),
        Err(PssmReadError::Ncbi(ReadError::Serial(message))) => {
            Err(RunError::Exception(message.into_bytes()))
        }
        Err(PssmReadError::Ncbi(ReadError::Eof)) => Err(RunError::Exception(
            b" Frame type= eFrameNamed, Object type= PssmWithParameters".to_vec(),
        )),
        Err(PssmReadError::Ncbi(ReadError::Exception(exception))) => {
            Err(RunError::Exception(exception.get_msg().to_vec()))
        }
        Err(PssmReadError::Ncbi(ReadError::Std(_))) => Err(RunError::Unknown),
        Err(PssmReadError::Refused(what))
        | Err(PssmReadError::Ncbi(ReadError::Unsupported(what))) => {
            Err(RunError::Refused(format!(
                "a profile that LOSAT's reader does not read ({}: {what})",
                String::from_utf8_lossy(filename)
            )))
        }
    }
}

/// The writer's diagnostics since the last call, posted.
fn post_writer(diag: &mut Diag, output_db: &mut WriteDb) {
    let text = output_db.take_diagnostics();
    if !text.is_empty() {
        diag.write_text(&text);
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
                "-scale {}, whose ceiling NCBI makeprofiledb converts to a 32-bit integer \
                 without a defined result",
                params.pssm_scale_factor
            ))
        })?,
    };
    let mut rps = rps_files::open_volume(&db_name, smps.len(), volume_params)?;
    built.vol_names.push(db_name.clone());
    let now = LocalTime::now().ok_or_else(|| {
        RunError::Refused("the clock cannot be read for the database's date".to_string())
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
        params.title.as_bytes(),
        now.index_date().as_bytes(),
        MAX_FILE_SIZE,
        params.db_ver,
        Platform::host(),
        indices,
    );
    output_db.set_max_file_size(MAX_FILE_SIZE);
    let result = add_profiles(params, diag, &mut rps, &mut output_db, smps);
    // `output_db->Close()` after the loop, or the writer's destructor in the unwinding of an
    // error (`~CWriteDB_Impl` closes it).
    let closed = output_db.close();
    post_writer(diag, &mut output_db);
    if let Err(close_error) = closed {
        // A `Close` that threw: the destructor's `Close` returns at once and the writer's members
        // close (the last volume under its numbered name, the LMDB files).
        output_db.destroy_after_failed_close();
        post_writer(diag, &mut output_db);
        let Err(error) = result else {
            return Err(close_error.into());
        };
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:100-108
        // In the unwinding the destructor posts a `CWriteDBException` and the first error goes
        // on; any other exception leaves the destructor and ends NCBI's process
        // (`std::terminate`), which LOSAT refuses.
        let creation = matches!(
            close_error,
            WriteDbError::Creation(_) | WriteDbError::CreationBytes(_)
        ) && !crate::writedb::writer::is_gi_only_between_keys(&close_error);
        if !creation {
            let close = RunError::from(close_error);
            return Err(RunError::Refused(format!(
                "the error {error:?} left a database whose close fails ({close:?}), where NCBI \
                 makeprofiledb 2.17.0 aborts"
            )));
        }
        diag.err_post_bytes(&close_error.message());
        return Err(error);
    }
    result?;
    rps.close()?;
    Ok(())
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
        let mut pssm_w_parameters = read_profile(&path, filename, diag)?;
        let sm = check_input_scoremat(&pssm_w_parameters, filename)?;
        let seq_size = pssm_w_parameters.pssm.num_columns as usize;
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
        let freq_ratios = pssm
            .intermediate_data
            .as_ref()
            .map_or(&[][..], |data| &data.freq_ratios[..]);
        rps.add_profile(&ProfileData {
            seq_size,
            num_rows: pssm.num_rows as usize,
            by_row: pssm.by_row.unwrap_or(false),
            scores: &final_data.scores,
            freq_ratios,
            kappa: final_data.kappa,
        })?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1395-1430
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:1210-1226
// Without `-taxid` and `-taxid_map` the taxid of the query's first `org` descriptor maps its
// first id (`AddTaxId`), and `FixTaxId` then gives each definition line the map's taxid of its
// ids (0 for none). With `-index false` the writer puts `gnl|BL_ORD_ID|<oid>` in place of the
// first id of the first definition line, which NCBI does past an empty list of ids; LOSAT
// refuses a query without ids in both modes.
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
    if query.ids.is_empty() {
        return Err(RunError::Refused(
            "a profile whose query has no Seq-id, whose definition line NCBI makeprofiledb \
             writes past an empty list"
                .to_string(),
        ));
    }
    let mut deflines = match &query.descr {
        Some(descr) => {
            let deflines = bioseq_deflines(&query.ids, descr).map_err(|error| match error {
                BioseqDeflinesError::NoDeflines => {
                    RunError::Exception(b"Error: No deflines provided.".to_vec())
                }
                BioseqDeflinesError::BinaryHeader => RunError::Refused(
                    "a profile whose query holds binary definition lines (a user object \
                     ASN1_BlastDefLine), which LOSAT does not read yet"
                        .to_string(),
                ),
            })?;
            output_db.post_extracted_deflines(&deflines);
            post_writer(diag, output_db);
            deflines
        }
        None => DeflineSet {
            deflines: vec![BlastDefline {
                ids: query.ids.clone(),
                ..BlastDefline::default()
            }],
        },
    };
    if !params.user_taxids {
        if let Some(taxid) = query.descr.as_deref().and_then(first_org_taxid) {
            let mut warnings = Vec::new();
            params
                .taxids
                .add_tax_id(&query.ids[0], taxid, &mut |warning| warnings.push(warning))
                .map_err(|exception| RunError::Exception(exception.get_msg().to_vec()))?;
            for warning in warnings {
                post_seq_id_warning(diag, &warning);
            }
        }
    }
    params.taxids.fix_tax_id(&mut deflines);
    let inst = &query.inst;
    let seq_data = inst.seq_data.as_ref().map(|data| match data {
        SeqData::Ncbieaa(letters) => bioseq::SeqData::Ncbieaa(letters.clone()),
        SeqData::Iupacaa(letters) => bioseq::SeqData::Iupacaa(letters.clone()),
        SeqData::Ncbistdaa(codes) => bioseq::SeqData::Ncbistdaa(codes.clone()),
        SeqData::Other(_) => bioseq::SeqData::Other,
    });
    let bioseq = Bioseq {
        local_id: String::new(),
        title: Vec::new(),
        is_aa: inst.mol.1 == 3,
        length: u64::from(inst.length.unwrap_or(0)),
        seq_data,
        parsed: Some(ParsedIds {
            ids: query.ids.clone(),
            hist_old_id: None,
            defline_field: Vec::new(),
            title_descr: None,
        }),
        past_seq_pos: None,
    };
    let added = output_db.add_sequence(bioseq);
    post_writer(diag, output_db);
    added?;
    output_db.set_deflines(deflines);
    Ok(())
}
