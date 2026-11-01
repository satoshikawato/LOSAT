//! `losat blastdbcmd`: NCBI BLAST+ 2.17.0 blastdbcmd, of which LOSAT serves `-info` (the
//! information of a database: `-db`, `-dbtype`, `-exact_length`, `-out`, `-logfile`), the
//! retrieval of sequences (`retrieve`: `-entry all`, `-entry`, `-entry_batch`, every
//! `-outfmt` but the ASN.1 texts, `-range`, `-strand`, `-line_length`, `-ctrl_a`,
//! `-target_only`, `-get_dups`, `-long_seqids`), `-h`, `-help` and `-version`. Its command
//! line is NCBI's (`args`, NCBI's argument reader), not clap's: `main` hands it the words
//! after `blastdbcmd`. The other uses of the program (`-list`, taxonomy, metadata, IPGs)
//! are refused, and so are the lookups of an accession or GI in the tables of a version 4
//! database built with `-parse_seqids` (line BDC's BDC2c) and a database that LOSAT does
//! not read (the masks of `makeblastdb -mask_data`, several names in `-db`).

mod args;
mod fasta;
mod format;
mod letters;
mod range;
mod retrieve;

use std::ffi::OsString;
use std::fs::{File, OpenOptions};
use std::io::Write;

use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType};

use args::{Args, ArgsError};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "BLASTDBCMD";

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:907-909
// ```c++
//     arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
//                   "BLAST database client, version " + CBlastVersion().Print());
// ```
/// The text of `blastdbcmd -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions in `args`). Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the
/// argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3829-3834
// ```c++
//     if (detailed) {
//         x.AddDetails(arr);
//     } else {
//         arr.push_back(kEmptyStr);
//         arr.push_back("Use '-help' to print detailed descriptions of command line arguments");
//     }
// ```
/// The text of `blastdbcmd -h`: the usage and the description, then NCBI's pointer to `-help`.
fn short_usage() -> String {
    let end = HELP.find("OPTIONAL ARGUMENTS").unwrap_or(HELP.len());
    format!(
        "{}Use '-help' to print detailed descriptions of command line arguments\n",
        &HELP[..end]
    )
}

/// Runs `blastdbcmd` with `words` (the command line after `losat blastdbcmd`) on the process's
/// standard streams; the exit status.
pub fn run(words: &[OsString]) -> i32 {
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    let status = run_with(words, &mut out, &mut err);
    let _ = out.flush();
    let _ = err.flush();
    status
}

/// Where NCBI's diagnostics go (`SetLogFile` of the `-logfile` value).
enum DiagTarget {
    /// The standard stream `cerr`: no `-logfile`, `-logfile -`, or a file that cannot be
    /// opened.
    Stderr,
    /// The file of `-logfile`, appended to.
    File(File),
    /// Nowhere: `-logfile ""` and `-logfile /dev/null` (a stream handler without a stream).
    Discard,
}

struct Diag<'a> {
    target: DiagTarget,
    stderr: &'a mut dyn Write,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608
// ```c++
//     // Check if a non-existing dir is specified
//     if ( !s_IsSpecialLogName(file_name) ) {
//         string dir = CFile(file_name).GetDir();
//         if ( !dir.empty()  &&  !CDir(dir).Exists() ) {
//             return false;
//         }
//     }
//     ...
//         // Check special filenames
//         if ( file_name.empty()  ||  file_name == "/dev/null" ) {
//             // no output
//             SetDiagStream(0, quick_flush, 0, 0, kLogName_None);
//         }
//         else if (file_name == "-") {
//             // output to stderr
//             SetDiagStream(&NcbiCerr, quick_flush, 0, 0, kLogName_Stderr);
//         }
//         else {
//             // output to file
//             unique_ptr<CFileDiagHandler> fhandler(new CFileDiagHandler());
//             if ( !fhandler->SetLogFile(file_name, eDiagFile_All, quick_flush) ) {
//                 ERR_POST_X(9, Info << "Failed to initialize log: " << file_name);
//                 return false;
//             }
// ```
/// The target of the `-logfile` value `name` (a file NCBI fails to open leaves `cerr`).
fn diag_target(name: Option<&str>) -> DiagTarget {
    match name {
        None | Some("-") => DiagTarget::Stderr,
        Some("" | "/dev/null") => DiagTarget::Discard,
        Some(name) => match OpenOptions::new().create(true).append(true).open(name) {
            Ok(file) => DiagTarget::File(file),
            Err(_) => DiagTarget::Stderr,
        },
    }
}

impl Diag<'_> {
    /// Posts one diagnostic line (`LOG_POST`, `ERR_POST`).
    fn post(&mut self, line: &str) {
        match &mut self.target {
            DiagTarget::File(file) => {
                let _ = writeln!(file, "{line}");
            }
            DiagTarget::Stderr => {
                let _ = writeln!(self.stderr, "{line}");
                let _ = self.stderr.flush();
            }
            DiagTarget::Discard => {}
        }
    }

    /// LOSAT's refusal: a diagnostic line, exit status 1.
    fn refuse(&mut self, reason: &str) -> i32 {
        self.post(&format!("Error: {reason}"));
        1
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1175-1177
    // ```c++
    //     SetDiagPostLevel(eDiag_Warning);
    //     SetDiagPostPrefix("blastdbcmd");
    // ```
    /// `ERR_POST(Error << message)` inside `Run`: the severity, the prefix of the
    /// application, the message.
    fn err_post(&mut self, message: &str) {
        self.post(&format!("Error: [blastdbcmd] {message}"));
    }

    /// `err_post` of a message that is not UTF-8.
    fn err_post_bytes(&mut self, message: &[u8]) {
        let mut line = b"Error: [blastdbcmd] ".to_vec();
        line.extend_from_slice(message);
        line.push(b'\n');
        match &mut self.target {
            DiagTarget::File(file) => {
                let _ = file.write_all(&line);
            }
            DiagTarget::Stderr => {
                let _ = self.stderr.write_all(&line);
                let _ = self.stderr.flush();
            }
            DiagTarget::Discard => {}
        }
    }
}

/// What NCBI's application reads before the arguments (`CNcbiApplicationAPI::AppMain`).
enum Preparsed {
    /// `-version`: print the version, exit status 0.
    Version,
    /// The words left for the argument reader, and whether there were words that LOSAT
    /// refuses although NCBI goes on reading: `-conffile FILE` (NCBI reads that registry file
    /// before the arguments) and `-dryrun` (NCBI reads the arguments, then runs nothing).
    Words {
        words: Vec<String>,
        conffile: bool,
        dryrun: bool,
    },
    /// A word LOSAT refuses where NCBI stops reading (`-version-full`, ...).
    Refused(&'static str),
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-866
// ```c++
//     // Check if logfile is set in the args.
//     m_LogFile = 0;
//     if (!m_DisableArgDesc && argc > 1  &&  argv  &&  diag != eDS_User) {
//         for (int i = 1;  i < argc;  i++) {
//             if ( !argv[i] ) {
//                 continue;
//             }
//             if ( NStr::strcmp(argv[i], s_ArgDelimiter) == 0 ) {
//                 break;
//             }
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
// ```
/// The file of `-logfile` (`m_LogFile`): the value of its last `-logfile FILE` or
/// `-logfile=FILE` before `--`.
fn log_file(words: &[String]) -> Option<&str> {
    let mut log = None;
    let mut at = 0;
    while at < words.len() {
        let word = words[at].as_str();
        if word == "--" {
            break;
        }
        if word == "-logfile" {
            at += 1;
            if let Some(name) = words.get(at) {
                log = Some(name.as_str());
            }
        } else if let Some(name) = word.strip_prefix("-logfile=") {
            log = Some(name);
        }
        at += 1;
    }
    log
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1015
// ```c++
//         for (int i = 1;  i < argc;  i++) {
//             if ( !argv[i] ) {
//                 continue;
//             }
//             if ( NStr::strcmp(argv[i], s_ArgDelimiter) == 0 ) {
//                 skip = true;
//             }
//             if (skip) {
//                 v[real_arg_index++] = argv[i];
//                 continue;
//             }
//             // Log file - ignore if diag is eDS_User - the user wants to
//             // take care about logging.
//             if ( diag != eDS_User  &&
//                 NStr::strcmp(argv[i], s_ArgLogFile) == 0 ) {
//                 if (!argv[++i]) {
//                     continue;
//                 }
//                 v[real_arg_index++] = argv[i - 1];
//                 v[real_arg_index++] = argv[i];
//                 // Configuration file
//             } else if ( NStr::strcmp(argv[i], s_ArgCfgFile) == 0 ) {
//                 if (!argv[++i]) {
//                     continue;
//                 }
//                 v[real_arg_index++] = argv[i - 1];
//                 v[real_arg_index++] = argv[i];
//                 conf_arg = argv[i];
//
//             }
//             else if (NStr::StartsWith(argv[i], s_ArgCfgFile)) {
//                 v[real_arg_index++] = argv[i];
//                 const char* a = argv[i] + strlen(s_ArgCfgFile);
//                 if (*a == '=') {
//                     conf_arg = ++a;
//                 }
//
//                 // Version
//             } else if ( NStr::strcmp(argv[i], s_ArgVersion) == 0 ) {
//                 delete[] v;
//                 // Print VERSION
//                 cout << GetFullVersion().Print( appname, CVersionAPI::fVersionInfo | CVersionAPI::fPackageShort );
//                 diag_context.DiscardMessages();
//                 return 0;
//
//                 // Full version
//             } else if ( NStr::strcmp(argv[i], s_ArgFullVersion) == 0 ) {
//                 ...
//             } else if ( NStr::strcmp(argv[i], s_ArgDryRun) == 0 ) {
//                 m_DryRun = true;
//
//                 // Save real argument
//             } else {
//                 v[real_arg_index++] = argv[i];
//             }
//         }
// ```
/// The words NCBI's application takes before its argument reader: `-version` anywhere before
/// `--` (but not as the value of `-logfile` or `-conffile`), the full versions, `-dryrun` and
/// `-conffile`. A `-logfile` or `-conffile` without a value is dropped. LOSAT refuses
/// `-conffile` (another registry file) and `-dryrun` (NCBI then runs nothing), as the search
/// programs do (line CLI DECISIONS 4), and the full versions (NCBI's build details).
fn preparse(words: &[String]) -> Preparsed {
    let mut kept = Vec::with_capacity(words.len());
    let mut conffile = false;
    let mut dryrun = false;
    let mut skip = false;
    let mut at = 0;
    while at < words.len() {
        let word = words[at].as_str();
        if word == "--" {
            skip = true;
        }
        if skip {
            kept.push(words[at].clone());
            at += 1;
            continue;
        }
        match word {
            "-logfile" | "-conffile" => {
                at += 1;
                if let Some(value) = words.get(at) {
                    kept.push(word.to_string());
                    kept.push(value.clone());
                    conffile |= word == "-conffile";
                }
            }
            "-version" => return Preparsed::Version,
            "-version-full" => return Preparsed::Refused("-version-full"),
            "-version-full-xml" => return Preparsed::Refused("-version-full-xml"),
            "-version-full-json" => return Preparsed::Refused("-version-full-json"),
            "-dryrun" => dryrun = true,
            _ => {
                conffile |= word.starts_with("-conffile=");
                kept.push(words[at].clone());
            }
        }
        at += 1;
    }
    Preparsed::Words {
        words: kept,
        conffile,
        dryrun,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1086-1112
// ```c++
//     catch (const CArgException& e) {
//         // Print USAGE and the exception error message
//         if ( e.GetErrCode() != CArgException::eNoValue &&  m_ArgDesc.get() ) {
//             x_AddDefaultArgs();
//             string str;
//             if ( !m_ArgDesc->IsSetMiscFlag(CArgDescriptions::fNoUsage) ) {
//                 m_ArgDesc->PrintUsage(str);
//                 cerr << str;
//             }
//             if ( !m_ArgDesc->IsSetMiscFlag(CArgDescriptions::fDupErrToCerr) ) {
//                 CStreamDiagHandler* errh =
//                     dynamic_cast<CStreamDiagHandler*>(GetDiagHandler());
//                 if (!errh  ||  errh->GetStream() != &cerr) {
//                     cerr << "Error in command-line arguments. "
//                         "See error logs for more details." << endl;
//                 }
//             }
//             ...
//             cerr << string(72, '=') << endl << endl;
//         }
//         SetDiagPostAllFlags(eDPF_Severity);
//         NCBI_REPORT_EXCEPTION_X(18, "", e);
//         got_exception = true;
//         exit_code = 1;
//     }
// ```
/// The run of one command line, its output on `out` and its errors on `err`; the exit status.
///
/// The order is NCBI's: the diagnostic file of `-logfile` is opened first (append; stderr when
/// it cannot be opened), then `-version` is read, then the arguments; an argument error prints
/// the usage and a line of `=` on stderr and NCBI's report of the exception as diagnostics.
/// LOSAT's refusals come where NCBI would start the work they stand for, and before anything is
/// written: the program's other uses after the arguments are read, the database where it is
/// opened (after `-out` is opened, as NCBI's `Run` opens it first).
pub(crate) fn run_with(words: &[OsString], out: &mut dyn Write, err: &mut dyn Write) -> i32 {
    let Some(words) = words
        .iter()
        .map(|word| word.to_str().map(str::to_string))
        .collect::<Option<Vec<String>>>()
    else {
        let _ = writeln!(
            err,
            "Error: a word of the command line is not UTF-8; NCBI BLAST+ reads its bytes, which \
             is not supported by LOSAT's {PROGRAM}"
        );
        return 1;
    };
    let mut diag = Diag {
        target: diag_target(log_file(&words)),
        stderr: err,
    };
    let (words, dryrun) = match preparse(&words) {
        Preparsed::Version => {
            let _ = write!(
                out,
                "blastdbcmd: 2.17.0+\n Package: LOSAT {}, output compatible with blast 2.17.0\n",
                env!("CARGO_PKG_VERSION")
            );
            return 0;
        }
        Preparsed::Refused(word) => return diag.refuse(&unsupported_option(word)),
        Preparsed::Words { conffile: true, .. } => {
            return diag.refuse(&unsupported_option("-conffile"))
        }
        Preparsed::Words { words, dryrun, .. } => (words, dryrun),
    };
    let args = match args::parse(&words) {
        Ok(args) => args,
        Err(ArgsError::Usage) => {
            let _ = out.write_all(short_usage().as_bytes());
            return 0;
        }
        Err(ArgsError::Help) => {
            let _ = out.write_all(HELP.as_bytes());
            return 0;
        }
        Err(ArgsError::Invalid { code, messages }) => {
            // `eNoValue` (which prints no usage) is not raised by this reader. The diagnostics
            // of a `-logfile` go to a stream that is not `cerr`: NCBI then says where they are.
            let _ = diag.stderr.write_all(short_usage().as_bytes());
            if !matches!(diag.target, DiagTarget::Stderr) {
                let _ = writeln!(
                    diag.stderr,
                    "Error in command-line arguments. See error logs for more details."
                );
            }
            let _ = write!(diag.stderr, "{}\n\n", "=".repeat(72));
            let _ = diag.stderr.flush();
            for message in &messages {
                diag.post(&format!("Error: {message}"));
            }
            if let Some(last) = messages.last() {
                diag.post(&format!("Error:  (CArgException::{code}) {last}"));
            }
            return 1;
        }
        Err(ArgsError::Refused(reason)) => return diag.refuse(&reason),
    };
    if dryrun {
        return diag.refuse(&unsupported_option("-dryrun"));
    }
    if let Some(reason) = unserved_use(&args) {
        return diag.refuse(&reason);
    }
    let settings =
        match crate::blastinput::ncbi_environment::check_ncbi_application_settings("blastdbcmd") {
            Ok(settings) => settings,
            Err(reason) => return diag.refuse(&reason),
        };
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => return diag.refuse(&refusal(reason)),
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1179-1180
    // ```c++
    //     try {
    //         CNcbiOstream& out = args["out"].AsOutputFile();
    // ```
    let mut file;
    let out: &mut dyn Write = match args.value("out").unwrap_or("-") {
        "-" => out,
        name => match File::create(name) {
            Ok(created) => {
                file = created;
                &mut file
            }
            Err(_) => {
                // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-620
                // ```c++
                // void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
                // {
                //     if ( !m_Ios ) {
                //         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
                //             "File is not accessible",AsString()));
                //     }
                // }
                // ```
                // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
                // ```c++
                //     catch (const CArgException& e) {                                        \
                //         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
                //         exit_code = BLAST_INPUT_ERROR;                                      \
                //     }                                                                       \
                // ```
                diag.post(&format!(
                    "Command line argument error: {}",
                    args::arg_message("out", "File is not accessible", name)
                ));
                return 1;
            }
        },
    };
    let path = match SearchPath::from_environment(registry) {
        Ok(path) => path,
        Err(error) => return database_error(&mut diag, error),
    };
    let name = args.value("db").unwrap_or("nr");
    let (db, volumes) = match open_database(name, args.value("dbtype").unwrap_or("guess"), &path) {
        Ok(opened) => opened,
        Err(error) => return database_error(&mut diag, error),
    };
    if !args.given("info") {
        return retrieve::process(&db, &args, &path, &mut diag, out);
    }
    // NCBI prints the date bytes of the index file; the reader keeps them as text
    // (`SeqDb::date`), with U+FFFD for bytes that are not UTF-8.
    if db.date().contains('\u{fffd}') {
        return diag.refuse(&refusal(&format!(
            "{name}: the date of the database is not UTF-8 text"
        )));
    }
    let mut text = Vec::new();
    write_info(&mut text, &db, &volumes);
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does this.
    let _ = out.write_all(&text);
    let _ = out.flush();
    0
}

/// LOSAT's refusal of the option `word`.
fn unsupported_option(word: &str) -> String {
    format!("the NCBI BLAST+ option {word} is not supported by LOSAT's {PROGRAM}")
}

/// LOSAT's refusal of a database (`algorithm::common::database::refusal`).
fn refusal(reason: &str) -> String {
    crate::algorithm::common::database::refusal(reason, PROGRAM).to_string()
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1180-1232
// ```c++
//         if (args["show_blastdb_search_path"]) {
//             ...
//         } else if (args["list"]) {
//             ...
//         }
//
//         if (args["info"]) {
//         	x_InitBlastDB();
//             x_PrintBlastDatabaseInformation();
//         }
//         else if (args["metadata"]) {
//             ...
//         else if (args["tax_info"]) {
//             ...
//         else if(args[kArgTaxIdList].HasValue() ||
//                 args[kArgTaxIdListFile].HasValue()) {
//             ...
//         else {
//         	x_InitBlastDB();
//        		status = x_ProcessSearchRequest();
//     	}
// ```
/// LOSAT's refusal of a use of blastdbcmd other than `-info` and the retrieval of sequences
/// (`-entry`, `-entry_batch`), named by the option that selects it in NCBI's `Run` (line
/// BDC's BDC3: the lists, the taxonomy, the metadata, the IPGs). With `-info`, NCBI's reader
/// refuses every option the information does not use but `-line_length`,
/// `-mask_sequence_with`, `-ipg` and `-ipg_batch`, which NCBI then ignores, and so does
/// LOSAT. Without a selector NCBI opens the database and then asks for one
/// (`retrieve`).
fn unserved_use(args: &Args) -> Option<String> {
    if args.given("info") {
        return None;
    }
    [
        "show_blastdb_search_path",
        "list",
        "metadata",
        "tax_info",
        "taxids",
        "taxidlist",
        "ipg",
        "ipg_batch",
    ]
    .into_iter()
    .find(|name| args.given(name))
    .map(|name| unsupported_option(&format!("-{name}")))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
// ```c++
// static CSeqDBImpl *
// s_SeqDBInit(const string       & dbname,
//             char                 prot_nucl,
//             ...
// {
//     CSeqDBImpl * impl = 0;
//
//     if (prot_nucl == '-') {
//         try {
//             prot_nucl = 'p';
//             impl = new CSeqDBImpl(dbname,
//                                   prot_nucl,
//                                   ...
//         }
//         catch(CSeqDBException &) {
//             prot_nucl = 'n';
//         }
//     }
//
//     if (! impl) {
//         impl = new CSeqDBImpl(dbname,
//                               prot_nucl,
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:524-531
// ```c++
//     CSeqDB::ESeqType seqtype = ParseMoleculeTypeString(args[kArgDbType].AsString());
//     m_BlastDb.Reset(new CSeqDBExpert(args[kArgDb].AsString(), seqtype));
// ```
/// Opens the database `name` (the `-db` value as given: blastdbcmd does not trim it) of the
/// type `dbtype` (`prot`, `nucl`, or `guess`: protein first, then nucleotide), with the list of
/// its volumes for `-info`. With `guess`, an error that NCBI raises too (`ncbi_message`) moves
/// on to the nucleotide database, and the nucleotide error is the one reported; an error of
/// LOSAT's is reported at once, since NCBI may open that protein database.
fn open_database(
    name: &str,
    dbtype: &str,
    path: &SearchPath,
) -> Result<(SeqDb, Vec<String>), SeqDbError> {
    // The search programs search the databases of several names as one; LOSAT's BLASTDBCMD
    // keeps one database (line BDC DECISIONS 6).
    if crate::seqdb::alias::lists_several_names(name) {
        return Err(SeqDbError::Unsupported {
            path: std::path::PathBuf::from(name),
            message: "the name lists several databases; LOSAT's BLASTDBCMD reads one database \
                      (a name with spaces goes between double quotes)"
                .to_string(),
        });
    }
    let types: &[SeqType] = match dbtype {
        "prot" => &[SeqType::Protein],
        "nucl" => &[SeqType::Nucleotide],
        _ => &[SeqType::Protein, SeqType::Nucleotide],
    };
    let mut types = types.iter().peekable();
    loop {
        let Some(&seqtype) = types.next() else {
            unreachable!("at least one sequence type is tried");
        };
        let opened = SeqDb::open_in(name, seqtype, path).and_then(|db| {
            let volumes = SeqDb::info_volumes(name, seqtype, path)?;
            Ok((db, volumes))
        });
        match opened {
            Err(error) if types.peek().is_some() && error.ncbi_message().is_some() => {}
            opened => return opened,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// ```c++
//     catch (const CSeqDBException& e) {                                      \
//         ...
//         else {                                                              \
//         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
//         	exit_code = BLAST_DATABASE_ERROR;                               \
//         }                                                                   \
//     }                                                                       \
// ```
/// A reader error: NCBI's own (`SeqDbError::ncbi_message`) as `BLAST Database error: <text>`
/// with exit status 2, else LOSAT's refusal.
fn database_error(diag: &mut Diag, error: SeqDbError) -> i32 {
    match error.ncbi_message() {
        Some(text) => {
            diag.post(&format!("BLAST Database error: {text}"));
            2
        }
        None => diag.refuse(&refusal(&error.to_string())),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1841-1850
// ```c++
//     if ((flags & NStr::fWithCommas)) {
//         int cnt = -1;
//         do {
//             if (++cnt == 3) {
//                 *--pos = ',';
//                 cnt = 0;
//             }
//             *--pos = '0' + value % 10;
//             value /= 10;
//         } while (value);
// ```
/// `n` in decimal with a comma between groups of three digits (`NStr::fWithCommas`).
fn with_commas(n: u64) -> String {
    let digits = n.to_string();
    let mut text = String::with_capacity(digits.len() + digits.len() / 3);
    for (at, digit) in digits.chars().enumerate() {
        if at > 0 && (digits.len() - at) % 3 == 0 {
            text.push(',');
        }
        text.push(digit);
    }
    text
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:533-573
// ```c++
//     static const NStr::TNumToStringFlags kFlags = NStr::fWithCommas;
//     const string kLetters = m_DbIsProtein ? "residues" : "bases";
//     const string kVersion = (m_BlastDb->GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) ? "5":"4";
//     ...
//     // Print basic database information
//     out << "Database: " << m_BlastDb->GetTitle() << endl
//         << "\t" << NStr::IntToString(m_BlastDb->GetNumSeqs(), kFlags)
//         << " sequences; ";
//         if(args["exact_length"])
//         	out << NStr::UInt8ToString(m_BlastDb->GetExactTotalLength(), kFlags);
//         else
//         	out << NStr::UInt8ToString(m_BlastDb->GetTotalLength(), kFlags);
//     out << " total " << kLetters << endl << endl
//         << "Date: " << m_BlastDb->GetDate()
//         << "\tLongest sequence: "
//         << NStr::IntToString(m_BlastDb->GetMaxLength(), kFlags) << " "
//         << kLetters << endl << endl;
//
//     out << "BLASTDB Version: " << kVersion << endl;
//     ...
//     out << m_BlastDb->GetAvailableMaskAlgorithmDescriptions();
//     ...
//     out << endl << "Volumes:" << endl;
//     ITERATE(vector<string>, file_name, volumes) {
//         out << "\t" << *file_name << endl;
//     }
// ```
/// The text of `-info` for `db`, whose `Volumes:` are `volumes`.
///
/// `-exact_length` changes nothing here: NCBI's exact total differs from the total only when
/// its totals are counted over a filtered list of sequences (an alias file with `GILIST`,
/// `SEQIDLIST`, `OIDLIST`, ..., seqdbimpl.cpp:104-131), which LOSAT's reader refuses. The
/// filtering algorithms are never listed: LOSAT's reader refuses the volumes that have them
/// (column files).
fn write_info(out: &mut Vec<u8>, db: &SeqDb, volumes: &[String]) {
    let letters = match db.seqtype() {
        SeqType::Protein => "residues",
        SeqType::Nucleotide => "bases",
    };
    out.extend_from_slice(b"Database: ");
    out.extend_from_slice(db.title());
    out.push(b'\n');
    let _ = write!(
        out,
        "\t{} sequences; {} total {letters}\n\nDate: {}\tLongest sequence: {} {letters}\n\n\
         BLASTDB Version: {}\n\nVolumes:\n",
        with_commas(u64::from(db.num_seqs())),
        with_commas(db.total_length()),
        db.date(),
        with_commas(u64::from(db.max_length())),
        if db.is_version5() { "5" } else { "4" },
    );
    for volume in volumes {
        let _ = writeln!(out, "\t{volume}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(line: &[&str]) -> Vec<String> {
        line.iter().map(|word| word.to_string()).collect()
    }

    #[test]
    fn commas_group_three_digits() {
        assert_eq!(with_commas(0), "0");
        assert_eq!(with_commas(999), "999");
        assert_eq!(with_commas(1_000), "1,000");
        assert_eq!(with_commas(186_063_409), "186,063,409");
    }

    #[test]
    fn the_short_usage_is_the_head_of_the_help() {
        let usage = short_usage();
        assert_eq!(usage.len(), 780);
        assert!(usage.starts_with("USAGE\n  blastdbcmd [-h] [-help] [-db dbname]"));
        assert!(usage.ends_with(
            "   BLAST database client, version 2.17.0+\n\nUse '-help' to print detailed \
             descriptions of command line arguments\n"
        ));
        assert_eq!(HELP.len(), 9450);
    }

    #[test]
    fn the_log_file_is_the_last_one_before_the_delimiter() {
        assert_eq!(
            log_file(&words(&["-logfile", "a", "-logfile=b"])),
            Some("b")
        );
        assert_eq!(log_file(&words(&["-logfile"])), None);
        assert_eq!(log_file(&words(&["--", "-logfile", "a"])), None);
        assert_eq!(log_file(&words(&["-logfilex", "a"])), None);
    }

    #[test]
    fn preparse_takes_ncbi_s_words() {
        let kept = |line: &[&str]| match preparse(&words(line)) {
            Preparsed::Words {
                words,
                conffile,
                dryrun,
            } => (words, conffile, dryrun),
            Preparsed::Version => panic!("version"),
            Preparsed::Refused(word) => panic!("refused {word}"),
        };
        assert!(matches!(
            preparse(&words(&["-db", "-version"])),
            Preparsed::Version
        ));
        // The value of `-logfile` is not read as `-version`; a last `-logfile` is dropped.
        assert_eq!(
            kept(&["-info", "-logfile", "-version"]).0,
            words(&["-info", "-logfile", "-version"])
        );
        assert_eq!(kept(&["-info", "-logfile"]).0, words(&["-info"]));
        assert_eq!(kept(&["--", "-version"]).0, words(&["--", "-version"]));
        assert_eq!(
            kept(&["-info", "-dryrun"]),
            (words(&["-info"]), false, true)
        );
        assert!(kept(&["-dryrun", "-conffile", "x"]).1);
        assert!(kept(&["-conffile=x"]).1);
        assert!(!kept(&["-info", "-conffile"]).1);
        assert!(matches!(
            preparse(&words(&["-version-full", "-version"])),
            Preparsed::Refused("-version-full")
        ));
    }
}
