//! `losat blastdbcmd`: NCBI BLAST+ 2.17.0 blastdbcmd, of which LOSAT serves `-info` (the
//! information of a database: `-db`, `-dbtype`, `-exact_length`, `-out`, `-logfile`), the
//! retrieval of sequences (`retrieve`: `-entry all`, `-entry`, `-entry_batch`, `-ipg`,
//! `-ipg_batch`, `-taxids`, `-taxidlist`, every `-outfmt` but the ASN.1 texts, `-range`,
//! `-strand`, `-line_length`, `-ctrl_a`, `-target_only`, `-get_dups`, `-long_seqids`), the
//! taxids of a database (`taxonomy`: `-tax_info`), its metadata (`metadata`: `-metadata`,
//! `-metadata_output_prefix`), the databases of a directory (`list`: `-list`,
//! `-list_outfmt`, `-recursive`, `-remove_redundant_dbs`), `-show_blastdb_search_path`, `-h`,
//! `-help` and `-version`. Its command line is NCBI's (`args`, NCBI's argument reader), not
//! clap's: `main` hands it the words after `blastdbcmd`. The masks of `makeblastdb -mask_data`
//! are read (`-info`, `-mask_sequence_with`, `%m`). A `-db` of several names is read as NCBI's
//! reader reads it, as one database (`SeqDb::open_in`).

mod args;
mod fasta;
mod format;
mod letters;
mod list;
mod metadata;
mod range;
mod retrieve;
mod taxonomy;

use std::ffi::OsString;
use std::fs::File;
use std::io::Write;

use crate::blastinput::ncbi_environment::{
    application_start_of_process, diag_target, log_file, preparse, DiagTarget, Preparsed,
};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType, TaxidList};

use args::ArgsError;

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

struct Diag<'a> {
    target: DiagTarget,
    stderr: &'a mut dyn Write,
}

impl Diag<'_> {
    /// Posts one diagnostic line (`LOG_POST`, `ERR_POST`).
    fn post(&mut self, line: &str) {
        self.target
            .write_text(format!("{line}\n").as_bytes(), self.stderr);
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
        self.post_bytes(&line);
    }

    /// `post` of a line that is not UTF-8.
    fn post_bytes(&mut self, message: &[u8]) {
        let mut line = message.to_vec();
        line.push(b'\n');
        self.target.write_text(&line, self.stderr);
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
    // NCBI's application start (`ncbi_environment::ApplicationStart`), which the command line's
    // entry has made: the file of `-logfile` gets the usage report's messages when NCBI opens it,
    // before `-version`.
    let start = application_start_of_process();
    if let Some(start) = start {
        start.write_log_setup(&mut diag.target);
    }
    let (words, dryrun) = match preparse(&words) {
        Preparsed::Version => {
            let _ = write!(
                out,
                "blastdbcmd: {}\n Package: LOSAT {}, output compatible with blast {}\n",
                crate::NCBI_BLAST_VERSION,
                env!("CARGO_PKG_VERSION"),
                crate::NCBI_BLAST_RELEASE
            );
            return 0;
        }
        Preparsed::Refused(word) => return diag.refuse(&unsupported_option(word)),
        // The file of `-conffile` is the configuration of the command line's application start
        // (`start_application`); the library API runs none (`cli::toolkit_file_word_refusal`).
        Preparsed::Words {
            conffile: Some(_), ..
        } if start.is_none() => {
            return diag.refuse(&crate::cli::toolkit_file_word_refusal(PROGRAM, "-conffile"))
        }
        Preparsed::Words { words, dryrun, .. } => (words, dryrun),
    };
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit` loads the
    // configuration before the arguments are read; quoted at `ApplicationStart`)
    // ```c++
    // void CNcbiApplicationAPI::x_TryInit(EAppDiagStream diag, const char* conf)
    // {
    //     // Load registry from the config file
    //     if ( conf ) {
    //         string x_conf(conf);
    //         LoadConfig(*m_Config, &x_conf);
    // ```
    // The configuration's messages go where the diagnostics go; a `blastdbcmd.ini` that NCBI
    // reports as malformed ends the run with exit status 2 (NCBI BLAST+ 2.17.0, g_reg probes
    // `o2`, `inibad_null`).
    if let Some(start) = start {
        match start.write_config_load(&mut diag.target, diag.stderr) {
            Ok(None) => {}
            Ok(Some(status)) => return status,
            Err(reason) => return diag.refuse(&reason),
        }
    }
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
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1181-1184
    // ```c++
    //         if (args["show_blastdb_search_path"]) {
    //             out << CSeqDB::GenerateSearchPath() << NcbiEndl;
    //             return status;
    // ```
    if args.given("show_blastdb_search_path") {
        let _ = writeln!(out, "{}", path.as_str());
        let _ = out.flush();
        return 0;
    }
    if args.given("list") {
        return list::list(&args, &path, &mut diag, out);
    }
    let name = args.value("db").unwrap_or("nr");
    let dbtype = args.value("dbtype").unwrap_or("guess");
    let selected = ["info", "metadata", "tax_info"]
        .into_iter()
        .any(|selector| args.given(selector));
    if !selected && (args.given("taxids") || args.given("taxidlist")) {
        // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:358-369
        // ```c++
        //     CSeqDB::ESeqType seqtype = ParseMoleculeTypeString(args[kArgDbType].AsString());
        //     m_DbIsProtein = static_cast<bool>(seqtype == CSeqDB::eProtein);
        //     m_TargetOnly = args["target_only"];
        //    	if(m_TargetOnly) {
        //     	CRef<CSeqDBGiList> taxid_list(new CSeqDBGiList());
        //     	taxid_list->AddTaxIds(m_TaxIdList);
        //    		m_BlastDb.Reset(new CSeqDBExpert(args[kArgDb].AsString(), seqtype, taxid_list.GetPointer()));
        //     }
        //    	else {
        //    		m_BlastDb.Reset(new CSeqDBExpert(args[kArgDb].AsString(), seqtype));
        //    	}
        // ```
        // With `-target_only` the list's errors are those of the open, inside NCBI's guess
        // of the type: a protein database the list fails on is then looked for as a
        // nucleotide one (probe J111 of line BDC: `No alias or index file found for
        // nucleotide database [db/p4]`).
        let tax_ids = match taxonomy::tax_id_list(&args, &path, &mut diag) {
            Ok(tax_ids) => tax_ids,
            Err(taxonomy::Ended(status)) => return status,
        };
        let target_only = args.given("target_only");
        let list = TaxidList::Positive(tax_ids.clone());
        let restrict = |db: &mut SeqDb| -> Result<(), SeqDbError> {
            if target_only {
                db.restrict_taxids(&list)?;
            }
            Ok(())
        };
        let (db, _) = match open_database(name, dbtype, &path, &restrict) {
            Ok(opened) => opened,
            Err(error) => return database_error(&mut diag, error),
        };
        let dbtype_protein = dbtype == "prot";
        return retrieve::process(
            &db,
            &args,
            &path,
            &mut diag,
            out,
            Some(&tax_ids),
            dbtype_protein,
        );
    }
    let (db, volumes) = match open_database(name, dbtype, &path, &|_| Ok(())) {
        Ok(opened) => opened,
        Err(error) => return database_error(&mut diag, error),
    };
    if !args.given("info") {
        if args.given("metadata") {
            return metadata::metadata(&db, name, &args, &path, &mut diag, out);
        }
        if args.given("tax_info") {
            return taxonomy::tax_info(&db, &args, &path, &mut diag, out);
        }
        let protein = db.seqtype() == SeqType::Protein;
        return retrieve::process(&db, &args, &path, &mut diag, out, None, protein);
    }
    // NCBI prints the date bytes of the index file; the reader keeps them as text
    // (`SeqDb::date`), with U+FFFD for bytes that are not UTF-8.
    if db.date().contains('\u{fffd}') {
        return diag.refuse(&refusal(&format!(
            "{name}: the date of the database is not UTF-8 text"
        )));
    }
    let mut text = Vec::new();
    write_info_head(&mut text, &db);
    let masks = db
        .mask_algorithms()
        .and_then(|algorithms| algorithms.descriptions());
    let tail = match &masks {
        Ok(masks) => {
            text.extend_from_slice(masks.as_bytes());
            write_info_volumes(&mut text, &volumes);
            None
        }
        Err(error) => Some(error),
    };
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does this.
    // The lines before a bad mask column are written when the error is NCBI's own.
    if tail.is_none_or(|error| error.ncbi_message().is_some()) {
        let _ = out.write_all(&text);
        let _ = out.flush();
    }
    match masks {
        Ok(_) => 0,
        Err(error) => database_error(&mut diag, error),
    }
}

/// LOSAT's refusal of the option `word`.
fn unsupported_option(word: &str) -> String {
    format!("the NCBI BLAST+ option {word} is not supported by LOSAT's {PROGRAM}")
}

/// LOSAT's refusal of a database (`algorithm::common::database::refusal`).
fn refusal(reason: &str) -> String {
    crate::algorithm::common::database::refusal(reason, PROGRAM).to_string()
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
/// Opens the database `name` (the `-db` value as given, one name or several: blastdbcmd does not
/// trim it) of the type `dbtype` (`prot`, `nucl`, or `guess`: protein first, then nucleotide),
/// with the list of its volumes for `-info`, and `prepare`s it (the taxid list of `-target_only`, which NCBI's
/// constructor applies). With `guess`, an error that NCBI raises too (`ncbi_message`) moves
/// on to the nucleotide database, and the nucleotide error is the one reported; an error of
/// LOSAT's is reported at once, since NCBI may open that protein database.
fn open_database(
    name: &str,
    dbtype: &str,
    path: &SearchPath,
    prepare: &dyn Fn(&mut SeqDb) -> Result<(), SeqDbError>,
) -> Result<(SeqDb, Vec<String>), SeqDbError> {
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
        // The alias files' counts fail where NCBI's constructor reaches them, which the list of
        // `-target_only` moves (`SeqDb::check_alias_values`).
        let opened = SeqDb::open_search_in(name, seqtype, path).and_then(|mut db| {
            let volumes = SeqDb::info_volumes(name, seqtype, path)?;
            prepare(&mut db)?;
            db.check_alias_values()?;
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
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
// ```c++
//     catch (const CException& e) {                                           \
//         LOG_POST(Error << "Error: " << e.what());                           \
//         exit_code = BLAST_UNKNOWN_ERROR;                                    \
//     }                                                                       \
// ```
/// A reader error: NCBI's own (`SeqDbError::ncbi_message`) as `BLAST Database error: <text>`
/// with exit status 2, another exception of NCBI's (`SeqDbError::ncbi_exception`: an alias
/// file's value that does not convert) as `Error: <what()>` with exit status 255, else LOSAT's
/// refusal.
fn database_error(diag: &mut Diag, error: SeqDbError) -> i32 {
    if let Some(exception) = error.ncbi_exception() {
        diag.post(&format!("Error: {}", exception.what()));
        return 255;
    }
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
        if at > 0 && (digits.len() - at).is_multiple_of(3) {
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
/// The text of `-info` for `db` up to its filtering algorithms.
///
/// The counts are `GetNumSeqs` and `GetTotalLength` of the database opened without a list:
/// the walks of the alias tree (`SeqDb::printed_num_seqs`). `-exact_length` changes nothing
/// here: NCBI's exact total differs from the total only when its totals are counted over a
/// filtered list of sequences (an alias file with `GILIST`, `SEQIDLIST`, `OIDLIST`, ...,
/// seqdbimpl.cpp:104-131), which LOSAT's reader refuses; the OID ranges (`FIRST_OID`,
/// `LAST_OID`) do not make NCBI scan.
fn write_info_head(out: &mut Vec<u8>, db: &SeqDb) {
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
         BLASTDB Version: {}\n",
        with_commas(u64::from(db.printed_num_seqs())),
        with_commas(db.printed_total_length()),
        db.date(),
        with_commas(u64::from(db.max_length())),
        if db.is_version5() { "5" } else { "4" },
    );
}

/// The text of `-info` after the filtering algorithms: the `Volumes:` are `volumes`.
fn write_info_volumes(out: &mut Vec<u8>, volumes: &[String]) {
    out.extend_from_slice(b"\nVolumes:\n");
    for volume in volumes {
        let _ = writeln!(out, "\t{volume}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn commas_group_three_digits() {
        assert_eq!(with_commas(0), "0");
        assert_eq!(with_commas(999), "999");
        assert_eq!(with_commas(1_000), "1,000");
        assert_eq!(with_commas(186_063_409), "186,063,409");
    }

    // help.txt is NCBI's frozen `-help` text and is not edited or split: this keeps the version
    // that it names and the version that LOSAT writes elsewhere from drifting apart.
    #[test]
    fn the_frozen_help_names_the_version_that_losat_writes() {
        assert!(HELP.contains(&format!(
            "   BLAST database client, version {}\n",
            crate::NCBI_BLAST_VERSION
        )));
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
}
