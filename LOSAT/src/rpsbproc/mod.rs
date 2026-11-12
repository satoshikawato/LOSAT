//! `losat rpsbproc`: NCBI's post-RPSBLAST processing utility rpsbproc 0.5.1 (the package of
//! 2025-05-13 on `ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/rpsbproc/`), which reads the ASN.1 text
//! archives of `rpsblast -outfmt 11` and the CDD annotation data and writes the conserved domain
//! hits of each query. LOSAT gives the domain hits (`-t doms`, with `-m rep`, `std` or `full`
//! and `-f`) and the sites and motifs (`-t feats`, `-t both`, the default) of protein queries;
//! nucleotide queries (rpstblastn's archives; item PS16c) and the deprecated XML input (item
//! PS16d) are refused explicitly.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1853-1871
//! ```c++
//! int main(int argc, char * argv[])
//! {
//! 	string launchAlias = CDirEntry::NormalizePath(argv[0]);
//! 	string base, ext;
//! 	CDirEntry::SplitPath(launchAlias, &CRpsbProcApp::m_st_launch_path, &base, &ext);
//! 	...
//! 	return CRpsbProcApp().AppMain(argc, argv);
//! }
//! ```
//!
//! NCBI's rpsbproc is not part of BLAST+ and its source is not in the pinned BLAST tree; the
//! references `RpsbProc 0.5.1, 2025-05-13` name the files of its source package
//! (`RpsbProc-src.tar.gz` of 2025-05-13, sha256 a20f4088...6086), the same source that bioconda
//! builds on the BLAST+ 2.17.0 toolkit.

mod archive;
mod args;
mod data;
mod domains;
mod prosite;
mod query;
mod segset;
mod sites;

pub(crate) use args::alias_word;
pub use args::RpsbprocArgs;

use std::ffi::{OsStr, OsString};
use std::io::{BufReader, Read, Write};

use crate::ncbi_exception::{err_post_x, NcbiException};

use archive::{ArchiveStream, ParamValue, Queries, RequestBody, SeqEntry, Subject};
use args::TargetData;
use data::{CddData, DataError, DataFiles};
use domains::{
    domain_records, feature_records, map_seg_set, translated_pos_map, DataMode, FeatureRecord,
    Record,
};
use query::{check_masks, choice_selection, process_blast_results, DomQuery, ProcessError};
use segset::SegSet;
use sites::{MotifError, MotifState};

/// rpsbproc's version (`VERDEF`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:74-87
// ```c++
// #define VERDEF 0, 5, 1
//
// static constexpr const struct TVER
// {
// 	int major;
// 	int minor;
// 	int patch_level;
// } VER = {VERDEF};
//
// std::ostream & operator << (std::ostream &os, const struct TVER &v)
// {
// 	os << 'v' << v.major << '.' << v.minor << '.' << v.patch_level;
// 	return os;
// }
// ```
pub const VERSION: &str = "0.5.1";

/// `EErrorSubCode`: the exit codes.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:91-101
// ```c++
// enum EErrorSubCode
// {
// 	e_OK = 0,
// 	e_Info,
// 	e_DataLoadError,
// 	e_InputFileOpenError,
// 	e_OutputFileOpenError,
// 	e_InputDataError,
// 	//-----------------------
// 	eErr_TotalErrorCodes
// };
// ```
const E_OK: i32 = 0;
const E_DATA_LOAD_ERROR: i32 = 2;
const E_INPUT_FILE_OPEN_ERROR: i32 = 3;
const E_OUTPUT_FILE_OPEN_ERROR: i32 = 4;
const E_INPUT_DATA_ERROR: i32 = 5;

/// `OVERLAPLEADING` as rpsbproc stores it (`unsigned int` of 10000000000).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:338
// ```c++
// 	constexpr static const unsigned int OVERLAPLEADING = 10000000000;
// ```
const OVERLAPLEADING: u32 = (10_000_000_000u64 % (1u64 << 32)) as u32;

/// How `losat rpsbproc` ends: an exit status, or a refusal (its message; exit 1).
pub enum Outcome {
    Exit(i32),
    Refused(String),
}

/// glibc's `%g` (an `ostream`'s default output of a double), `-0` for negative zero.
fn format_g(value: f64) -> String {
    if value == 0.0 && value.is_sign_negative() {
        return "-0".to_string();
    }
    crate::report::pairwise::cpp_default_double(value)
}

/// The writer of the report (`m_ostream`): standard output or the `-o` file. NCBI's `cout` is
/// flushed by each `endl` and its `cerr` is unbuffered; LOSAT writes the report in blocks and
/// flushes it before each line it writes to standard error, so the two streams interleave as
/// NCBI's do.
///
/// A failed write stops LOSAT (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 9): NCBI's
/// `ostream` sets its bad bit and the run goes on to exit 0 (a closed pipe ends it by SIGPIPE).
struct Report {
    out: Box<dyn Write>,
    quiet: bool,
    /// The first write error; nothing is written after it.
    failed: Option<std::io::Error>,
}

impl Report {
    fn write(&mut self, bytes: &[u8]) {
        if self.failed.is_none() {
            if let Err(error) = self.out.write_all(bytes) {
                self.failed = Some(error);
            }
        }
    }

    fn line(&mut self, parts: &[&[u8]]) {
        for part in parts {
            self.write(part);
        }
        self.write(b"\n");
    }

    fn flush(&mut self) {
        if self.failed.is_none() {
            if let Err(error) = self.out.flush() {
                self.failed = Some(error);
            }
        }
    }

    /// A progress line on standard error (not with `-q`, nor after a failed write).
    fn progress(&mut self, text: &str) {
        if !self.quiet {
            self.flush();
            if self.failed.is_none() {
                log_post(text);
            }
        }
    }

    /// Whether a write failed: the run stops (`run` reports it).
    fn stopped(&mut self) -> bool {
        self.flush();
        self.failed.is_some()
    }
}

/// `LOG_POST_X`: rpsbproc's messages are its text alone on standard error.
fn log_post(text: &str) {
    write_stderr(&format!("{text}\n"));
}

/// Writes `text` to standard error. A failed write (a reader that has closed its pipe) ends
/// LOSAT at once with exit status 4 and no message, where NCBI's rpsbproc dies of SIGPIPE
/// (exit status 141 from the shell; PS16a review oracle M1): PD-LOSAT-CLI-NONSEARCH-DIFFERENCES
/// exception 9. The report is flushed before every line that rpsbproc writes to standard error
/// (`Report::progress`, `Report::stopped`), so it holds no bytes then.
pub(crate) fn write_stderr(text: &str) {
    let mut stderr = std::io::stderr().lock();
    if stderr.write_all(text.as_bytes()).is_err() {
        std::process::exit(E_OUTPUT_FILE_OPEN_ERROR);
    }
}

/// Writes `text` (help or version) to standard output; a failed write ends LOSAT as a failed
/// write of the report does (exception 9).
pub(crate) fn write_stdout(text: &str) {
    let mut stdout = std::io::stdout().lock();
    if let Err(error) = stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.flush())
    {
        write_stderr(&format!("rpsbproc failed to write output: {error}\n"));
        std::process::exit(E_OUTPUT_FILE_OPEN_ERROR);
    }
}

/// `NCBI_USE_ERRCODE_X` of ncbiapp.cpp (`Corelib_App`), the error code of the application's
/// reports.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:65
// ```c++
// #define NCBI_USE_ERRCODE_X   Corelib_App
// ```
// NCBI reference (598d8ae6): c++/include/corelib/error_codes.hpp:49
// ```c++
// NCBI_DEFINE_ERRCODE_X(Corelib_App,        106, 24);
// ```
const CORELIB_APP: u32 = 106;

/// The report of an exception that leaves rpsbproc's `Run` (`CNcbiApplication::x_TryMain`):
/// the exception's own line, then the application's, exit 3. For the exceptions of `Init`
/// (`initialization`) the application's line has subcode 15 and the exit is 2.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:669-674,711-728
// ```c++
//             catch (const CException& e) {
//                 NCBI_REPORT_EXCEPTION_X(15,
//                                         "Application's initialization failed", e);
//                 *got_exception = true;
//                 *exit_code = 2;
//             }
//             ...
//             catch (const CException& e) {
//                 ...
//                 NCBI_REPORT_EXCEPTION_X(16,
//                                         "Application's execution failed", e);
//                 *got_exception = true;
//                 *exit_code = 3;
//             }
// ```
// (rpsbproc sets no diagnostic flags, so the report is `NcbiException::report_default_ex`'s; the
// lines as NCBI's rpsbproc printed them: `Error: (CInvalidChoiceSelection::eFail) <message>` and
// `Error: (106.16) Application's execution failed (CInvalidChoiceSelection::eFail) <message>`
// (probe 2, a database with `lcl|` ids); `Error: (CStringException::eConvert) <message>` and
// `Error: (106.15) Application's initialization failed (CStringException::eConvert) <message>`
// (probe 4, data files with numbers that do not convert).)
fn report_exception(exception: &NcbiException, initialization: bool) -> i32 {
    let (subcode, title, status) = if initialization {
        (15, "Application's initialization failed", 2)
    } else {
        (16, "Application's execution failed", 3)
    };
    write_stderr(&exception.report_default_ex(CORELIB_APP, subcode, title));
    status
}

/// `CDir(path).Exists()`.
fn dir_exists(path: &OsStr) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_dir())
}

/// `CFile(path).Exists()`.
fn file_exists(path: &OsStr) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_file())
}

/// `a + "/" + b` of `std::string`s.
fn join(a: &OsStr, b: &str) -> OsString {
    let mut out = a.to_os_string();
    out.push("/");
    out.push(b);
    out
}

/// `CDirEntry::NormalizePath(path)` (`eIgnoreLinks`): the path's `.` and empty components
/// dropped and `..` taken back where possible, without its trailing `/`
/// (`ncbi_environment::normalize_path`, the Unix form, on every target).
// NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:501-502
// ```c++
//     static string NormalizePath(const string& path,
//                                 EFollowLinks  follow_links = eIgnoreLinks);
// ```
fn normalize_path(path: &[u8]) -> Vec<u8> {
    crate::blastinput::ncbi_environment::normalize_path(path, None)
        .expect("a path normalized without links has no link limit")
}

/// The paths where rpsbproc looks for its data without `-d`: the launch directory of the
/// program (`argv[0]` normalized) and its real directory (links followed, made absolute).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1855-1868
// ```c++
// 	string launchAlias = CDirEntry::NormalizePath(argv[0]);
// 	string base, ext;
// 	CDirEntry::SplitPath(launchAlias, &CRpsbProcApp::m_st_launch_path, &base, &ext);
// 	if (!CRpsbProcApp::m_st_launch_path.empty())
// 		CRpsbProcApp::m_st_launch_path.pop_back();	//get rid of the tailing '/'
// 	// -- this is deliberately leave alone
// 	CRpsbProcApp::m_st_exec_name = base + ext;
//
// 	CDir launchPathDir(launchAlias);
// 	launchPathDir.DereferencePath();
//
// 	CDirEntry::SplitPath(CDirEntry::CreateAbsolutePath(launchPathDir.GetPath(), CDirEntry::eRelativeToCwd), &CRpsbProcApp::m_st_real_path);
// 	if (!CRpsbProcApp::m_st_real_path.empty())
// 		CRpsbProcApp::m_st_real_path.pop_back();
// ```
// The binary is `losat` (or a link `rpsbproc` to it, `--install-shims`). The launch path is
// ported as NCBI makes it. For the real path NCBI follows the links of each component of the
// path (`s_DereferencePath`, ncbifile.cpp) relative to the current directory, normalizes it
// lexically and puts the current directory before a relative one; LOSAT takes the canonical
// path of the launch path (`fs::canonicalize`), else the current directory with the path. They
// differ only where a `..` follows a linked directory (the canonical path goes up from the
// link's target, NCBI's lexical normalization from the link).
fn program_paths() -> (OsString, OsString) {
    let argv0 = std::env::args_os().next().unwrap_or_default();
    let alias = normalize_path(argv0.as_encoded_bytes());
    let dir_of = |path: &[u8]| -> Vec<u8> {
        match path.iter().rposition(|&b| b == b'/') {
            Some(at) => path[..at].to_vec(),
            None => Vec::new(),
        }
    };
    // SAFETY: the bytes are an `OsStr`'s (`argv[0]`, or a path made of its parts and `/`).
    let to_os = |bytes: Vec<u8>| unsafe { OsString::from_encoded_bytes_unchecked(bytes) };
    let launch = dir_of(&alias);
    let alias_path = std::path::PathBuf::from(to_os(alias.clone()));
    let real_path = match std::fs::canonicalize(&alias_path) {
        Ok(path) => path.into_os_string().into_encoded_bytes(),
        Err(_) => {
            let mut absolute = std::env::current_dir()
                .map(|dir| dir.into_os_string().into_encoded_bytes())
                .unwrap_or_default();
            if alias.first() == Some(&b'/') {
                absolute = alias.clone();
            } else {
                absolute.push(b'/');
                absolute.extend_from_slice(&alias);
            }
            normalize_path(&absolute)
        }
    };
    (to_os(launch), to_os(dir_of(&real_path)))
}

/// The data directory (`#Annotation data at:`), searched as rpsbproc searches it.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:572-643
// ```c++
// 	m_init_msg << "#Annotation data at:";
//
// 	if (AV_Data.HasValue())	//take command line designated data path as first priority
// 	{
// 		const string &datapath = AV_Data.AsString();
// 		if (!datapath.empty())
// 		{
// 			actual_path = datapath + "/" + _datadir;
// 			if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())
// 				goto labelDataPathFound;
//
// 			actual_path = datapath;
// 			if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())
// 				goto labelDataPathFound;
//
// 			LOG_POST_X(e_DataLoadError, "Specified data location (" << datapath << ") does not exist or does not contain data files");
// 			m_app_status = e_DataLoadError;
// 			m_init_msg << "(error)";
// 			actual_path.clear();
// 		}
// 		...
// 	}
// 	else	//not specified in command line, so try home directory (where the binary located) first
// 	{
// 		// -- first look for current directory for a "data" subdirectory
// 		actual_path = "./" + _datadir;
//
// 		if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())	//found
// 			goto labelDataPathFound;
//
// 		actual_path = "./";
// 		if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())	//found
// 			goto labelDataPathFound;
//
// 		actual_path = CRpsbProcApp::m_st_launch_path + "/" + _datadir;
// 		if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())	//found
// 			goto labelDataPathFound;
//
// 		actual_path = CRpsbProcApp::m_st_launch_path;
//
// 		actual_path = CRpsbProcApp::m_st_launch_path + "/" + _datadir;
// 		if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())	//found
// 			goto labelDataPathFound;
//
// 		if (CRpsbProcApp::m_st_real_path != CRpsbProcApp::m_st_launch_path)
// 		{
// 			actual_path = CRpsbProcApp::m_st_real_path + "/" + _datadir;
// 			if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())	//not there
// 				goto labelDataPathFound;
//
// 			actual_path = CRpsbProcApp::m_st_real_path;
// 			if (CDir(actual_path).Exists() && CFile(actual_path + "/" + _cdd).Exists())	//not there
// 				goto labelDataPathFound;
// 		}
// 		actual_path.clear();
// 	}
// ```
// (The registry's `[datafiles]` names, main:535-566, are rejected with every other registry
// entry that changes the output: `check_ncbi_application_settings`.)
fn find_data_path(data_path: Option<&OsStr>) -> Option<OsString> {
    let found = |path: &OsStr| dir_exists(path) && file_exists(&join(path, "cddid.tbl"));
    match data_path {
        Some(datapath) => {
            if datapath.is_empty() {
                return None;
            }
            let candidate = join(datapath, "data");
            if found(&candidate) {
                return Some(candidate);
            }
            if found(datapath) {
                return Some(datapath.to_os_string());
            }
            log_post(&format!(
                "Specified data location ({}) does not exist or does not contain data files",
                datapath.to_string_lossy()
            ));
            None
        }
        None => {
            for candidate in [OsString::from("./data"), OsString::from("./")] {
                if found(&candidate) {
                    return Some(candidate);
                }
            }
            let (launch, real) = program_paths();
            let candidate = join(&launch, "data");
            if found(&candidate) {
                return Some(candidate);
            }
            if real != launch {
                let candidate = join(&real, "data");
                if found(&candidate) {
                    return Some(candidate);
                }
                if found(&real) {
                    return Some(real);
                }
            }
            None
        }
    }
}

/// The 9 lines that rpsbproc writes to standard error first (also with `-q`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:397-406
// ```c++
// 	if (!m_silent)
// 	{
// 		cerr << PROGRAM_TITLE << ' ' << VER << endl;
// 		cerr << "-------------------------------------------------------" << endl;
// 		cerr << "Download:" << endl;
// 		cerr << "\tftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/" << endl << endl;
// 		cerr << "Blast applications download:" << endl;
// 		cerr << "\tftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/" << endl << endl;
// 		cerr << "For more information please refer to the README file" << endl;
// 	}
// ```
// (`m_silent` is still false here: `-q` is read at main:677.)
fn banner() {
    write_stderr(&format!(
        "Post-RPSBLAST Processing Utility v{VERSION}\n\
         -------------------------------------------------------\n\
         Download:\n\
         \tftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/\n\
         \n\
         Blast applications download:\n\
         \tftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n\
         \n\
         For more information please refer to the README file\n"
    ));
}

/// The header of the report.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:679-716
// ```c++
// 	if (e_OK == m_app_status)
// 	{
//
// 		m_init_msg << "#E-Value cutoff:\t" << m_evalue << endl;
// 		m_init_msg << "#Redundancy:\t" << GetDisplay<TDataModes> (m_dmode) << endl;
// 		m_init_msg << "#Data requested:\t" << GetDisplay<TTargetData>(m_tdata) << endl;
// 		m_init_msg << "#Output format -- tab-delimited table" << endl;
// 		m_init_msg << "#Show superfamilies: " << (m_show_clusters ? "YES" : "NO") << endl;
// 		m_init_msg << '#' << DATASTART << endl;
// 		m_init_msg << '#' << SESSIONSTART << DELIMIT << "<session-ordinal>" << DELIMIT << "<program>" << DELIMIT << "<database>" << DELIMIT << "<score-matrix>" << DELIMIT << "<evalue-threshold>" << endl;
// 		m_init_msg << '#' << QUERYSTART << DELIMIT << "<query-id>" << DELIMIT << "<seq-type>" << DELIMIT << "<seq-length>" << DELIMIT << "<definition-line>" << endl;
//
//
// 		if (TTargetData::e_feats != m_tdata)
// 		{
// 			m_init_msg << '#' << DOMSTART << endl;
// 			m_init_msg << '#' << "<session-ordinal>" << DELIMIT << "<query-id[readingframe]>" << DELIMIT << "<hit-type>" << DELIMIT << "<PSSM-ID>" << DELIMIT << "<from>" << DELIMIT << "<to>" << DELIMIT << "<E-Value>" << DELIMIT << "<bitscore>" << DELIMIT << "<accession>" << DELIMIT << "<short-name>" << DELIMIT << "<incomplete>" << DELIMIT << "<superfamily PSSM-ID>" << endl;
// 			m_init_msg << "#more such lines......" << endl;
// 			m_init_msg << '#' << DOMEND << endl;
// 		}
// 		...
// 		m_init_msg << '#' << QUERYEND << DELIMIT << "<query-id>" << endl;
// 		m_init_msg << "#more query sections.." << endl;
// 		m_init_msg << '#' << SESSIONEND << DELIMIT << "<session-ordinal>" << endl;
// 		m_init_msg << "#more session sections.." << endl;
// 		m_init_msg << '#' << DATAEND << endl;
// 		m_init_msg << "#=====================================================================" << endl;
// 	}
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:699-711
// ```c++
// 		if (TTargetData::e_doms != m_tdata)
// 		{
// 			m_init_msg << '#' << FEATSTART << endl;
// 			m_init_msg << '#' << "<session-ordinal>" << DELIMIT << "<query-id[readingframe]>" << DELIMIT << "<annot-type>" << DELIMIT << "<title>" << DELIMIT << "<residue(coordinates)>" << DELIMIT << "<complete-size>" << DELIMIT << "<mapped-size>" << DELIMIT << "<source-domain>" << endl;
// 			m_init_msg << "#more such lines......" << endl;
// 			m_init_msg << '#' << FEATEND << endl;
// 			m_init_msg << '#' << MOTIFSTART << endl;
// 			m_init_msg << '#' << "<session-ordinal>" << DELIMIT << "<query-id[readingframe]>" << DELIMIT << "<title>" << DELIMIT << "<from>" << DELIMIT << "<to>" << DELIMIT << "<source-domain>" << endl;
// 			m_init_msg << "#more such lines......" << endl;
// 			m_init_msg << '#' << MOTIFEND << endl;
// 		}
// ```
fn header_tail(evalue: f64, mode: DataMode, show_families: bool, target: TargetData) -> String {
    let mode = match mode {
        DataMode::Rep => "Concise",
        DataMode::Std => "Standard",
        DataMode::Full => "Full",
    };
    let (requested, domains, sites) = match target {
        TargetData::Doms => ("Domain hits", true, false),
        TargetData::Feats => ("Site annotations", false, true),
        TargetData::Both => ("Domain hits and site annotations", true, true),
    };
    let mut text = format!(
        "#E-Value cutoff:\t{}\n\
         #Redundancy:\t{mode}\n\
         #Data requested:\t{requested}\n\
         #Output format -- tab-delimited table\n\
         #Show superfamilies: {}\n\
         #DATA\n\
         #SESSION\t<session-ordinal>\t<program>\t<database>\t<score-matrix>\t<evalue-threshold>\n\
         #QUERY\t<query-id>\t<seq-type>\t<seq-length>\t<definition-line>\n",
        format_g(evalue),
        if show_families { "YES" } else { "NO" }
    );
    if domains {
        text.push_str(
            "#DOMAINS\n\
             #<session-ordinal>\t<query-id[readingframe]>\t<hit-type>\t<PSSM-ID>\t<from>\t<to>\t<E-Value>\t<bitscore>\t<accession>\t<short-name>\t<incomplete>\t<superfamily PSSM-ID>\n\
             #more such lines......\n\
             #ENDDOMAINS\n",
        );
    }
    if sites {
        text.push_str(
            "#SITES\n\
             #<session-ordinal>\t<query-id[readingframe]>\t<annot-type>\t<title>\t<residue(coordinates)>\t<complete-size>\t<mapped-size>\t<source-domain>\n\
             #more such lines......\n\
             #ENDSITES\n\
             #MOTIFS\n\
             #<session-ordinal>\t<query-id[readingframe]>\t<title>\t<from>\t<to>\t<source-domain>\n\
             #more such lines......\n\
             #ENDMOTIFS\n",
        );
    }
    text.push_str(
        "#ENDQUERY\t<query-id>\n\
         #more query sections..\n\
         #ENDSESSION\t<session-ordinal>\n\
         #more session sections..\n\
         #ENDDATA\n\
         #=====================================================================\n",
    );
    text
}

/// The first variable of `variables` that changes how NCBI's object stream reads an archive or
/// what its objects give, with a value that changes it; `None` when there is none.
// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:57-68
// ```c++
// NCBI_PARAM_DECL(bool, OBJECTS, SEQ_DESCR_ALLOW_EMPTY);
// NCBI_PARAM_DEF_EX(bool, OBJECTS, SEQ_DESCR_ALLOW_EMPTY, false,
//                   eParam_NoThread, OBJECTS_SEQ_DESCR_ALLOW_EMPTY);
// static CSafeStatic<NCBI_PARAM_TYPE(OBJECTS, SEQ_DESCR_ALLOW_EMPTY)> s_SeqDescrAllowEmpty;
//
// void CSeq_descr::PostRead(void) const
// {
//     if ( !s_SeqDescrAllowEmpty->Get() && Get().empty() ) {
//         NCBI_THROW(CSerialException, eInvalidData,
//                    "empty Seq-descr is not allowed");
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:257-287
// ```c++
// ESerialVerifyData CObjectIStream::x_GetVerifyDataDefault(void)
// {
//     ESerialVerifyData now = TSerialVerifyData::GetThreadDefault();
//     if (now == eSerialVerifyData_Default) {
//         now = TSerialVerifyData::GetDefault();
//         if (now == eSerialVerifyData_Default) {
// // this is to provide compatibility with old implementation
//             const char* str = getenv(SERIAL_VERIFY_DATA_READ);
//             if (str) {
//                 if (NStr::CompareNocase(str,"YES") == 0) {
//                     now = eSerialVerifyData_Yes;
//                 } else if (NStr::CompareNocase(str,"NO") == 0) {
//                     now = eSerialVerifyData_No;
//                 } else if (NStr::CompareNocase(str,"NEVER") == 0) {
//                     now = eSerialVerifyData_Never;
//                 } else  if (NStr::CompareNocase(str,"ALWAYS") == 0) {
//                     now = eSerialVerifyData_Always;
//                 } else  if (NStr::CompareNocase(str,"DEFVALUE") == 0) {
//                     now = eSerialVerifyData_DefValue;
//                 } else  if (NStr::CompareNocase(str,"DEFVALUE_ALWAYS") == 0) {
//                     now = eSerialVerifyData_DefValueAlways;
//                 }
//             }
//         }
//     }
//     if (now != eSerialVerifyData_Default) {
//         return now;
//     }
//     // change the default here, if you like
//     return eSerialVerifyData_Yes;
// }
// ```
// (`CSerialObject::x_GetVerifyData`, serialobject.cpp:180-221, reads `SERIAL_VERIFY_DATA_GET`
// the same way; `NO`, `NEVER`, `DEFVALUE` and `DEFVALUE_ALWAYS` turn the checks off: a missing
// member is no longer a read error, and the getter of an unset member gives its default instead
// of throwing `CUnassignedMember`. `YES`, `ALWAYS` and other values keep NCBI's default.)
// `OBJECTS_SEQ_DESCR_ALLOW_EMPTY` lets an empty `Seq-descr` read; LOSAT rejects it with any value.
// The other parameters of the readers are registry entries (`NCBI_CONFIG__SERIAL__...`, which
// `check_ncbi_application_settings` rejects) or change no result (`SERIAL_READ_MMAPBYTESOURCE`,
// `SERIAL_DISABLE_DELAY_BUFFERS` and the binary ASN.1 ones); the gi parameters of the Seq-ids
// are refused in `run_sessions`, where a query has a gi.
fn reader_variable(variables: impl Iterator<Item = (OsString, OsString)>) -> Option<String> {
    for (name, value) in variables {
        let (Some(name), value) = (name.to_str(), value.to_string_lossy()) else {
            continue;
        };
        let changes = match name {
            "OBJECTS_SEQ_DESCR_ALLOW_EMPTY" => true,
            "SERIAL_VERIFY_DATA_READ" | "SERIAL_VERIFY_DATA_GET" => {
                ["NO", "NEVER", "DEFVALUE", "DEFVALUE_ALWAYS"]
                    .iter()
                    .any(|off| value.eq_ignore_ascii_case(off))
            }
            _ => false,
        };
        if changes {
            return Some(format!(
                "the environment variable {name}={value}, which changes how NCBI's object stream reads an archive, is not supported by LOSAT's rpsbproc"
            ));
        }
    }
    None
}

/// Runs `losat rpsbproc` with `args` (NCBI's `Init` and `Run`).
pub fn run(args: RpsbprocArgs) -> Outcome {
    // File names that are not UTF-8 are refused as the other programs refuse them (NCBI writes
    // their bytes in its messages and the report's header).
    for (role, name) in [
        ("i", &args.infile),
        ("o", &args.outfile),
        ("d", &args.data_path),
    ] {
        if let Some(name) = name {
            if let Err(error) = crate::blastinput::input_files::check_utf8_file_name(
                std::path::Path::new(name),
                role,
                "rpsbproc",
            ) {
                return Outcome::Refused(error.to_string());
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044 (quoted at the BLAST
    // programs in `cli::main_entry`): the application layer reads its environment and its
    // registry files (`rpsbproc.ini`, `.ncbirc`) before `Init`.
    let settings =
        match crate::blastinput::ncbi_environment::check_ncbi_application_settings("rpsbproc") {
            Ok(settings) => settings,
            Err(reason) => return Outcome::Refused(reason),
        };
    if let Some(reason) = reader_variable(std::env::vars_os()) {
        return Outcome::Refused(reason);
    }
    banner();
    let mut status = E_OK;
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:447-498
    // ```c++
    // 	m_init_msg << "#Input data file:";
    // 	if (AV_InFile.HasValue())
    // 	{
    // 		const string &infname = AV_InFile.AsString();
    // 		if (CFile(infname).Exists())
    // 		{
    // 			ifstream ifs(infname.c_str(), ios::binary | ios::in);
    // 			if (ifs.good())
    // 			{
    // 				m_infiles.AppendFile(infname);
    // 				m_init_msg << '\t' << infname;
    // 			}
    // 			else
    // 			{
    // 				LOG_POST_X(e_InputFileOpenError, "Unable to open input file: " << infname);
    // 				m_app_status = e_InputFileOpenError;
    // 			}
    // 			ifs.close();
    // 		}
    // 		else
    // 		{
    // 			LOG_POST_X(e_InputFileOpenError, "Input file does not exist: " << infname);
    // 			m_app_status = e_InputFileOpenError;
    // 		}
    // 	}
    //
    // 	if (m_infiles.NLeft() > 0)
    // 		m_istream = m_infiles.OpenStream();
    // 	else
    // 		m_init_msg << "\tstdin";
    //
    // 	m_init_msg << endl;
    //
    // 	m_init_msg << "#Output data file:";
    //
    // 	if (AV_OutFile.HasValue())
    // 	{
    // 		const string & ofilename = AV_OutFile.AsString();
    // 		m_outfile.open(ofilename.c_str(), ios::binary | ios::out);	//append! not write.
    // 		if (m_outfile.good())
    // 		{
    // 			m_ostream = &m_outfile;
    // 			m_init_msg << '\t' << ofilename;
    // 		}
    // 		else
    // 		{
    // 			LOG_POST_X(e_OutputFileOpenError, "Unable to open output file: " << ofilename);
    // 			m_app_status = e_OutputFileOpenError;
    // 		}
    // 	}
    // 	else
    // 		m_init_msg << "\t(stdout)";
    // ```
    let mut header: Vec<u8> = b"#Post-RPSBLAST Processing Utility\n\n#Input data file:".to_vec();
    let mut input: Option<Box<dyn Read>> = None;
    if let Some(infile) = &args.infile {
        let name = infile.to_string_lossy();
        if file_exists(infile) {
            match std::fs::File::open(infile) {
                Ok(file) => {
                    input = Some(Box::new(file));
                    header.push(b'\t');
                    header.extend_from_slice(infile.as_encoded_bytes());
                }
                Err(_) => {
                    log_post(&format!("Unable to open input file: {name}"));
                    status = E_INPUT_FILE_OPEN_ERROR;
                }
            }
        } else {
            log_post(&format!("Input file does not exist: {name}"));
            status = E_INPUT_FILE_OPEN_ERROR;
        }
    }
    let input: Box<dyn Read> = match input {
        Some(file) => file,
        None => {
            header.extend_from_slice(b"\tstdin");
            Box::new(std::io::stdin())
        }
    };
    header.extend_from_slice(b"\n#Output data file:");
    let mut out: Box<dyn Write> = Box::new(std::io::stdout());
    if let Some(outfile) = &args.outfile {
        match std::fs::File::create(outfile) {
            Ok(file) => {
                out = Box::new(file);
                header.push(b'\t');
                header.extend_from_slice(outfile.as_encoded_bytes());
            }
            Err(_) => {
                log_post(&format!(
                    "Unable to open output file: {}",
                    outfile.to_string_lossy()
                ));
                status = E_OUTPUT_FILE_OPEN_ERROR;
            }
        }
    } else {
        header.extend_from_slice(b"\t(stdout)");
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:501-522
    // ```c++
    // 	// -- check stream type
    // 	// -- test if the stream contains xml or asn
    // 	char c = 0;
    // 	while (m_istream->good() && c <= ' ')
    // 		m_istream->get(c);
    //
    // 	if (!m_istream->good())
    // 	{
    // 		m_app_status = e_InputDataError;
    // 		return;
    // 	}
    //
    // 	if ('<' == c)	//xml
    // 	{
    // 		LOG_POST_X(e_Info, "Warning: XML format support is deprecated!");
    // 		m_sfmt = eSerial_Xml;
    // 	}
    //
    // 	m_istream->putback(c);
    //
    //
    // 	m_init_msg << endl;
    // ```
    // (`c` is a signed `char`: the bytes from 0x80 are skipped as white space.)
    let mut input = BufReader::with_capacity(1 << 16, input);
    let first = loop {
        let mut byte = [0u8];
        match input.read(&mut byte) {
            Ok(1) => {
                if (byte[0] as i8) > b' ' as i8 {
                    break Some(byte[0]);
                }
            }
            Ok(_) => break None,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(_) => break None,
        }
    };
    let Some(first) = first else {
        return Outcome::Exit(E_INPUT_DATA_ERROR);
    };
    // XML is refused only when NCBI's run would read it: after a failed `-i` or `-o` NCBI posts
    // its warning, loads the data and ends with the status (PS16a review oracle M4).
    if first == b'<' {
        if status == E_OK {
            return Outcome::Refused(
                "the XML input of rpsblast -outfmt 5 (which NCBI's rpsbproc reads as deprecated) is not supported by LOSAT's rpsbproc, which reads the ASN.1 text archives of -outfmt 11".to_string(),
            );
        }
        log_post("Warning: XML format support is deprecated!");
    }
    header.push(b'\n');
    header.extend_from_slice(b"#Annotation data at:");
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:645-669
    // ```c++
    // labelDataPathFound:
    //
    // 	if (!actual_path.empty())
    // 	{
    // 		m_init_msg << '\t' << actual_path;
    // 		try
    // 		{
    // 			m_dom_cluster_src.LoadData_real(actual_path + "/" + _cdd, actual_path + "/" + _feats, actual_path + "/" + _genfeats, actual_path + "/" + _spthr, actual_path + "/" + _cdt, actual_path + "/" + _clst);
    // 		}
    // 		catch (CSimpleException e)
    // 		{
    // 			LOG_POST_X(e_DataLoadError, "Annotation data loading unsuccessful: " << e.what());
    // 			m_app_status = e_DataLoadError;
    // 			m_init_msg << "(error)";
    // 		}
    // 	}
    // 	else
    // 	{
    // 		m_init_msg << "\tNot found!";
    // 		LOG_POST_X(e_DataLoadError, "Cannot find data files");
    // 		m_app_status = e_DataLoadError;
    // 		m_init_msg << "(error)";
    // 	}
    // ```
    let mut data = CddData::default();
    match find_data_path(args.data_path.as_deref()) {
        Some(path) => {
            header.push(b'\t');
            header.extend_from_slice(path.as_encoded_bytes());
            let base = path.to_string_lossy().into_owned();
            let file = |name: &str| format!("{base}/{name}");
            let files = DataFiles {
                cddids: file("cddid.tbl"),
                spfeats: file("cddannot.dat"),
                genfeats: file("cddannot_generic.dat"),
                bscores: file("bitscore_specific.txt"),
                cdtrack: file("cdtrack.txt"),
                famlinks: file("family_superfamily_links"),
            };
            if path.to_str().is_none() {
                return Outcome::Refused(
                    "a data path that is not UTF-8 is not supported by LOSAT's rpsbproc"
                        .to_string(),
                );
            }
            match CddData::load(&files) {
                Ok(loaded) => data = loaded,
                Err(DataError::Missing(path)) => {
                    log_post(&format!(
                        "Annotation data loading unsuccessful: Data file {path} does not exist"
                    ));
                    status = E_DATA_LOAD_ERROR;
                }
                Err(DataError::Conversion(error)) => {
                    return Outcome::Exit(report_exception(&error, true));
                }
                // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:675-679
                // ```c++
                //             catch (const exception& e) {
                //                 ERR_POST_X(6, "Application's initialization failed: " << e.what());
                //                 *got_exception = true;
                //                 *exit_code = 2;
                //             }
                // ```
                // (A `std::exception` of the loader is no `CSimpleException`, so it leaves
                // `Init`; NCBI's rpsbproc printed `Error: (106.6) ...` with the `what()` of
                // libstdc++, PS16a review oracle A4-1 cases `annot10`, `gen10` and `ct3`.)
                Err(DataError::Std(what)) => {
                    write_stderr(&err_post_x(
                        CORELIB_APP,
                        6,
                        &format!("Application's initialization failed: {what}"),
                    ));
                    return Outcome::Exit(2);
                }
                Err(DataError::Unsupported(message)) => return Outcome::Refused(message),
            }
        }
        None => {
            log_post("Cannot find data files");
            status = E_DATA_LOAD_ERROR;
        }
    }
    header.push(b'\n');
    if status != E_OK {
        return Outcome::Exit(status);
    }
    header.extend_from_slice(
        header_tail(args.evalue, args.mode(), args.show_families, args.target()).as_bytes(),
    );
    let mut report = Report {
        out: Box::new(std::io::BufWriter::with_capacity(1 << 16, out)),
        quiet: args.quiet,
        failed: None,
    };
    let outcome = run_sessions(
        &mut report,
        &header,
        ArchiveStream::new(input, &[first]),
        &data,
        args.mode(),
        args.show_families,
        args.target(),
        settings.seq_id_rank_variable,
    );
    // PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 9 (version 1.10): the first failed write of
    // the report stops LOSAT with this message and exit 4, rpsbproc's code of the output file
    // (`e_OutputFileOpenError`); NCBI prints nothing for it.
    report.flush();
    if let Some(error) = report.failed {
        write_stderr(&format!("rpsbproc failed to write output: {error}\n"));
        return Outcome::Exit(E_OUTPUT_FILE_OPEN_ERROR);
    }
    outcome
}

/// The name of a `Blast4-value` variant.
fn value_variant(value: &ParamValue) -> &'static str {
    match value {
        ParamValue::String(_) => "string",
        ParamValue::Cutoff(_) => "cutoff",
        ParamValue::Other(name) => name,
    }
}

/// The sessions: one per archive (`Run`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:727-739
// ```c++
// int CRpsbProcApp::Run(void)
// {
// 	if (e_OK == m_app_status)
// 	{
// 		// -- output initial message
// 		(*m_ostream) << m_init_msg.str();
// 		(*m_ostream) << DATASTART << endl;
// 		m_init_msg.clear();
//
// 		//CNcbiCdAlignProcessor proc(&m_dom_cluster_src);
// 		COfflCdAlignProcessor proc(&m_dom_cluster_src);
// 		//list<_TOfflDomQuery> qobjs;
// 		size_t objCounter = 0, qCounter = 0;
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:793-917
// ```c++
// 				if (!m_silent) cerr << "RPSBlast session " << objCounter << " start..." << endl;
// 				...
// 					switch (reqbody.Which())
// 					{
// 					case CBlast4_request_body::e_Queue_search:
// 						{
// 							const CBlast4_queue_search_request &qsreq = reqbody.GetQueue_search();
//
// 							// -- start session output
// 							(*m_ostream) << SESSIONSTART << DELIMIT << ++objCounter << DELIMIT << qsreq.GetProgram();
//
// 							if (req.IsSetIdent() && req.CanGetIdent())
// 								(*m_ostream) << DELIMIT << req.GetIdent();
//
// 							(*m_ostream) << DELIMIT << qsreq.GetSubject().GetDatabase();
//
// 							if (qsreq.IsSetAlgorithm_options() && qsreq.CanGetAlgorithm_options())
// 							{
// 								const CBlast4_parameters &algoParams = qsreq.GetAlgorithm_options();
// 								CRef<CBlast4_parameter> refGCode = algoParams.GetParamByName("MatrixName");
// 								if (!refGCode.IsNull())
// 									(*m_ostream) << DELIMIT << refGCode->GetValue().GetString();
//
// 								refGCode = algoParams.GetParamByName("EvalueThreshold");
// 								if (!refGCode.IsNull())
// 									(*m_ostream) << DELIMIT << refGCode->GetValue().GetCutoff().GetE_value();
// 							}
//
// 							(*m_ostream) << endl;
//
//
// 							const CBlast4_queries &blast_queries = qsreq.GetQueries();
//
// 							switch (blast_queries.Which())
// 							{
// 							case CBlast4_queries::e_Bioseq_set:
// 								pQuerySeqs = &(blast_queries.GetBioseq_set());
// 								break;
// 							default:
// 								LOG_POST_X(e_InputDataError, "Wrong type of request data: CBlast4_queries::e_Bioseq_set expected");
// 								m_app_status = e_InputDataError;
// 								goto labelReturn;
// 							}	//switch (blast_queries.Which())
// 						}
// 						break;
// 					default:
// 						LOG_POST_X(e_InputDataError, "Wrong type of request data: CBlast4_request_body::e_Queue_search expected");
// 						m_app_status = e_InputDataError;
// 						goto labelReturn;
// 					}	//switch (reqbody.Which())
//
//
// 					if (objBlastOutput.IsSetResults() && objBlastOutput.CanGetResults())
// 					{
// 						const CBlast4_archive::TResults & resobj = objBlastOutput.GetResults();	//
//
// 						if (resobj.IsSetMasks() && resobj.CanGetMasks())
// 						{
// 							const list< CRef< CBlast4_mask > > & blmasks = resobj.GetMasks();
//
// 							if (!blmasks.empty())
// 								ConvertBlastMaskListToSeqLocInfoVector(blmasks, vecMasks);
// 						}
//
// 						if (resobj.IsSetAlignments() && resobj.CanGetAlignments())
// 						{
// 							pAligns = &(resobj.GetAlignments());
// 						}
// 						else
// 						{
// 							LOG_POST_X(e_InputDataError, "Wrong type of results data: CBlast4_get_search_results_reply::e_alignments expected");
//
// 							goto labelNextObj;
// 						}
//
// 					}
// 					...
// 				if (!m_silent) cerr << "Calculating...." << endl;
// 				qobjs.clear();
// 				proc.ProcessBlastResults(qobjs, pQuerySeqs->GetSeq_set(), vecMasks, pAligns->Get(), gcode);
//
// 				for (const auto & v: qobjs)
// 				{
// 					PrintDomQuery(v, objCounter);
// 					++qCounter;
// 				}
//
// 			labelNextObj:
// 				(*m_ostream) << SESSIONEND << DELIMIT << objCounter << endl;
// 				if (!m_silent) cerr << "RPSBlast session " << objCounter << " done..." << endl;
// 			}
// 		}
//
//
// 		(*m_ostream) << DATAEND << endl;
// 		if (!m_silent)
// 			cerr << "All done, " << qCounter << " queries processed" << endl << endl;
// 		// -- do real job, if anything goes wrong, set m_app_status and goto labelReturn;
// 	}
// labelReturn:
// 	return m_app_status;
// }
// ```
// The `request` and `results` members are mandatory in the module, so an object that reads has
// both (main:886-896 are not reached). The refusals that a session's request decides (a
// translated search, a nucleotide query: item PS16c) come before anything of that session.
fn run_sessions<R: Read>(
    report: &mut Report,
    header: &[u8],
    mut stream: ArchiveStream<R>,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    seq_id_rank_variable: Option<&'static str>,
) -> Outcome {
    // NCBI's function-static motif positions (`MotifCheck`): one list for the whole run.
    let mut motif_state = MotifState::default();
    report.write(header);
    report.line(&[b"DATA"]);
    let mut session = 0usize;
    let mut queries_done = 0usize;
    loop {
        if report.stopped() {
            return Outcome::Exit(E_OK);
        }
        let archive = match stream.next_archive() {
            Ok(Some(archive)) => archive,
            Ok(None) => break,
            Err(refusal) => return Outcome::Refused(refusal.message),
        };
        if let RequestBody::QueueSearch(search) = &archive.body {
            if search.program == b"tblastn" {
                return Outcome::Refused(
                    "an archive of rpstblastn (program tblastn, nucleotide queries) is not supported by LOSAT's rpsbproc (item PS16c)".to_string(),
                );
            }
            if let Queries::BioseqSet(entries) = &search.queries {
                if entries
                    .iter()
                    .any(|entry| matches!(entry, SeqEntry::Seq(query) if query.mol != 3))
                {
                    return Outcome::Refused(
                        "a nucleotide query is not supported by LOSAT's rpsbproc (item PS16c)"
                            .to_string(),
                    );
                }
                // `SEQ_ID_PREFER_ACCESSION_OVER_GI` and `SEQ_ID_AVOID_GI` change the rank of a gi
                // (`ApplicationSettings::seq_id_rank_variable`, Seq_id.cpp:3064-3174), so the
                // label of a query with a gi among its ids (`query_accession`).
                if let Some(name) = seq_id_rank_variable {
                    if entries.iter().any(|entry| {
                        matches!(entry, SeqEntry::Seq(query) if query.ids.iter().any(|id| matches!(id, crate::seqid::model::ParsedSeqId::Gi(_))))
                    }) {
                        return Outcome::Refused(format!(
                            "the environment variable {name}, which changes the Seq-id that NCBI's rpsbproc prints for a query with a gi, is not supported by LOSAT's rpsbproc"
                        ));
                    }
                }
            }
        }
        report.progress(&format!("RPSBlast session {session} start..."));
        let RequestBody::QueueSearch(search) = &archive.body else {
            log_post("Wrong type of request data: CBlast4_request_body::e_Queue_search expected");
            return Outcome::Exit(E_INPUT_DATA_ERROR);
        };
        session += 1;
        let session_text = session.to_string();
        report.write(b"SESSION\t");
        report.write(session_text.as_bytes());
        report.write(b"\t");
        report.write(&search.program);
        if let Some(ident) = &archive.ident {
            report.write(b"\t");
            report.write(ident);
        }
        let process_failure = |report: &mut Report, error: ProcessError| -> Outcome {
            if report.stopped() {
                return Outcome::Exit(E_OK);
            }
            match error {
                ProcessError::Exception(exception) => {
                    Outcome::Exit(report_exception(&exception, false))
                }
                ProcessError::Unsupported(message) => Outcome::Refused(message),
            }
        };
        match &search.subject {
            Subject::Database(database) => {
                report.write(b"\t");
                report.write(database);
            }
            Subject::Other(name) => {
                return process_failure(
                    report,
                    choice_selection("NCBI-Blast4", "Blast4-subject", "database", name),
                );
            }
        }
        if search.has_algorithm_options {
            match &search.matrix_name {
                None => {}
                Some(ParamValue::String(name)) => {
                    report.write(b"\t");
                    report.write(name);
                }
                Some(other) => {
                    return process_failure(
                        report,
                        choice_selection(
                            "NCBI-Blast4",
                            "Blast4-value",
                            "string",
                            value_variant(other),
                        ),
                    );
                }
            }
            match &search.evalue_threshold {
                None => {}
                Some(ParamValue::Cutoff(archive::Cutoff::EValue(evalue))) => {
                    report.write(b"\t");
                    report.write(format_g(*evalue).as_bytes());
                }
                Some(ParamValue::Cutoff(archive::Cutoff::RawScore)) => {
                    return process_failure(
                        report,
                        choice_selection("NCBI-Blast4", "Blast4-cutoff", "e-value", "raw-score"),
                    );
                }
                Some(other) => {
                    return process_failure(
                        report,
                        choice_selection(
                            "NCBI-Blast4",
                            "Blast4-value",
                            "cutoff",
                            value_variant(other),
                        ),
                    );
                }
            }
        }
        report.write(b"\n");
        let Queries::BioseqSet(entries) = &search.queries else {
            if report.stopped() {
                return Outcome::Exit(E_OK);
            }
            log_post("Wrong type of request data: CBlast4_queries::e_Bioseq_set expected");
            return Outcome::Exit(E_INPUT_DATA_ERROR);
        };
        if let Some(masks) = &archive.results.masks {
            if let Err(error) = check_masks(masks) {
                return process_failure(report, error);
            }
        }
        match &archive.results.alignments {
            None => {
                if report.stopped() {
                    return Outcome::Exit(E_OK);
                }
                log_post(
                    "Wrong type of results data: CBlast4_get_search_results_reply::e_alignments expected",
                );
            }
            Some(aligns) => {
                report.progress("Calculating....");
                let queries = match process_blast_results(entries, aligns, data) {
                    Ok(queries) => queries,
                    Err(error) => return process_failure(report, error),
                };
                for query in &queries {
                    if let Err(error) = print_dom_query(
                        report,
                        query,
                        &session_text,
                        data,
                        mode,
                        show_families,
                        target,
                        &mut motif_state,
                    ) {
                        if report.stopped() {
                            return Outcome::Exit(E_OK);
                        }
                        return match error {
                            // `THROW_SIMPLE` leaves `Run`: `CNcbiApplication` reports a
                            // `std::exception` of the execution (probe 3 of PS16b: the line
                            // below, exit 3, the report ends before the site's line).
                            PrintError::MotifParse { motif, position } => {
                                write_stderr(&err_post_x(
                                    CORELIB_APP,
                                    7,
                                    &format!(
                                        "Application's execution failed: Motif string parse error -- Motif = {}, error position: {position}",
                                        String::from_utf8_lossy(&motif)
                                    ),
                                ));
                                Outcome::Exit(3)
                            }
                            PrintError::Unsupported(message) => Outcome::Refused(message),
                        };
                    }
                    queries_done += 1;
                    if report.stopped() {
                        return Outcome::Exit(E_OK);
                    }
                }
            }
        }
        report.line(&[b"ENDSESSION\t", session_text.as_bytes()]);
        report.progress(&format!("RPSBlast session {session} done..."));
    }
    report.line(&[b"ENDDATA"]);
    report.progress(&format!("All done, {queries_done} queries processed\n"));
    Outcome::Exit(E_OK)
}

/// `PrintDomQuery` of a protein query with `-t doms`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:947-1031
// ```c++
// void CRpsbProcApp::PrintDomQuery(const TDomQuery &q, size_t idxBlObj) const
// {
// 	static constexpr const SeqLen_t codonEnd = READINGFRAME::RF_SIZE - 1;
//
// 	(*m_ostream) << QUERYSTART << DELIMIT << q.m_strAccession << DELIMIT;
// 	...
// 	if (!q.m_bIsNa)
// 	{
// 		(*m_ostream) << QUERY_TYPE_PEPTIDE << DELIMIT << q.m_uiSeqLen << DELIMIT << q.m_strTitle << endl;
// 		vector<TDomSeqAlignIndex::__TCdAlignRecord> vecDomAligns, vecFeatAligns;	//vecDomFams for output of superfamilies of domain hits
// 		q.m_dimSplitAligns[0].CreateRecordSets(q.m_vecAlignments, m_dom_cluster_src, vecDomAligns, vecFeatAligns, m_dmode);
//
// 		// - so that we don't output redundent domain families
// 		map<unsigned int, TDomSeqAlignIndex::__TCdAlignRecord > mapDomFams;
//
// 		if (TTargetData::e_feats != m_tdata)
// 		{
// 			if (!vecDomAligns.empty())
// 			{
// 				(*m_ostream) << DOMSTART << endl;
// 				if (TDataModes::e_rep == m_dmode)
// 				{
// 					for (vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterRec = vecDomAligns.begin(), iterRecEnd = vecDomAligns.end(); iterRecEnd != iterRec; ++iterRec)
// 					{
// 						if (NULL == iterRec->pCdInfo) continue;
// 						(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << DELIMIT;
// 						if (1 == iterRec->pAlign->m_iRepClass)	//multi, should be no more
// 							PrintDomLine(HITTYPE_MULTIDOM, *iterRec);
// 						else if (iterRec->pAlign->m_bSpecQualified)
// 						{
// 							PrintDomLine(HITTYPE_SPECIFIC, *iterRec);
// 							if (NULL != iterRec->pClst)
// 							{
// 								unsigned int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCdInfo->m_uiClusterPssmId;
// 								//if (mapDomFams.end() == mapDomFams.find(key))
// 								mapDomFams.emplace(key, *iterRec);
// 							}
// 						}
// 						else if (NULL != iterRec->pClst)
// 							PrintClstLine(*iterRec);
// 						else	//not qualified for specific, but no cluster -- non-specific
// 							PrintDomLine(HITTYPE_NONSPECIFIC, *iterRec);
// 						(*m_ostream) << endl;
// 					}	//dom loops end
// 				}	//rep mode
// 				else	//nonrep mode
// 				{
// 					for (vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterRec = vecDomAligns.begin(), iterRecEnd = vecDomAligns.end(); iterRecEnd != iterRec; ++iterRec)
// 					{
// 						if (NULL == iterRec->pCdInfo) continue;
// 						(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << DELIMIT;
// 						if (1 == iterRec->pAlign->m_iRepClass)	//multi
// 							PrintDomLine(HITTYPE_MULTIDOM, *iterRec);
// 						else if (iterRec->pAlign->m_bRep && iterRec->pAlign->m_bSpecQualified)
// 							PrintDomLine(HITTYPE_SPECIFIC, *iterRec);
// 						else
// 							PrintDomLine(HITTYPE_NONSPECIFIC, *iterRec);
// 						if (NULL != iterRec->pClst)
// 						{
// 							unsigned int key = iterRec->pAlign->m_iRegionIdx * OVERLAPLEADING + iterRec->pCdInfo->m_uiClusterPssmId;
// 							//if (mapDomFams.end() == mapDomFams.find(key))
// 							mapDomFams.emplace(key, *iterRec);
// 						}
// 						(*m_ostream) << endl;
// 					}
// 				}
//
// 				(*m_ostream) << DOMEND << endl;
//
// 				if (m_show_clusters)	//user want to see all superfams
// 				{
// 					(*m_ostream) << FAMSTART << endl;
// 					for (const auto & v : mapDomFams)
// 					{
// 						(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << DELIMIT;
// 						PrintClstLine(v.second);
// 						(*m_ostream) << endl;
// 					}
// 					(*m_ostream) << FAMEND << endl;
// 				}
// 			}
// 		}
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1801-1802
// ```c++
// 	(*m_ostream) << QUERYEND << endl;
// 	if (!m_silent) cerr << "End processing Query " << q.m_strTitle << endl;
// ```
// `m_uiClusterPssmId` is never set by the offline data (`TDomain` keeps its 0): the key of a
// superfamily is its region times `OVERLAPLEADING`, an `unsigned int` product that wraps, and
// `emplace` keeps the first record of a key. `Multidom` (`m_iRepClass` 1) is never made
// (`domains::calculate`), so it is never printed.
#[allow(clippy::too_many_arguments)]
fn print_dom_query(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    show_families: bool,
    target: TargetData,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let length = query.length.to_string();
    report.line(&[
        b"QUERY\t",
        &query.accession,
        b"\tPeptide\t",
        length.as_bytes(),
        b"\t",
        &query.title,
    ]);
    let records = if target == TargetData::Feats {
        Vec::new()
    } else {
        domain_records(&query.index, &query.aligns, data, mode)
    };
    if !records.is_empty() {
        let mut families: std::collections::BTreeMap<u32, Record> = Default::default();
        report.line(&[b"DOMAINS"]);
        for record in &records {
            let align = &query.aligns[record.align];
            report.write(session.as_bytes());
            report.write(b"\t");
            report.write(&query.accession);
            report.write(b"\t");
            let key = (align.region as u32).wrapping_mul(OVERLAPLEADING);
            match mode {
                DataMode::Rep => {
                    if align.spec_qualified {
                        print_dom_line(report, b"Specific", record, query, data);
                        if record.cluster.is_some() {
                            families.entry(key).or_insert(*record);
                        }
                    } else if record.cluster.is_some() {
                        print_cluster_line(report, record, query, data);
                    } else {
                        print_dom_line(report, b"Non-specific", record, query, data);
                    }
                }
                DataMode::Std | DataMode::Full => {
                    let kind: &[u8] = if align.rep && align.spec_qualified {
                        b"Specific"
                    } else {
                        b"Non-specific"
                    };
                    print_dom_line(report, kind, record, query, data);
                    if record.cluster.is_some() {
                        families.entry(key).or_insert(*record);
                    }
                }
            }
            report.write(b"\n");
        }
        report.line(&[b"ENDDOMAINS"]);
        if show_families {
            report.line(&[b"SUPERFAMILIES"]);
            for record in families.values() {
                report.write(session.as_bytes());
                report.write(b"\t");
                report.write(&query.accession);
                report.write(b"\t");
                print_cluster_line(report, record, query, data);
                report.write(b"\n");
            }
            report.line(&[b"ENDSUPERFAMILIES"]);
        }
    }
    if target != TargetData::Doms {
        print_features(report, query, session, data, mode, motif_state)?;
    }
    report.line(&[b"ENDQUERY"]);
    if report.quiet || report.stopped() {
        return Ok(());
    }
    let mut stderr = std::io::stderr().lock();
    let _ = stderr.write_all(b"End processing Query ");
    let _ = stderr.write_all(&query.title);
    let _ = stderr.write_all(b"\n");
    Ok(())
}

/// Why the report of a query stops.
enum PrintError {
    /// A motif that NCBI's parser refuses (`THROW_SIMPLE`).
    MotifParse { motif: Vec<u8>, position: usize },
    /// A case LOSAT does not reproduce.
    Unsupported(String),
}

impl From<MotifError> for PrintError {
    fn from(error: MotifError) -> Self {
        match error {
            MotifError::Parse { motif, position } => PrintError::MotifParse { motif, position },
            MotifError::ReadPast => PrintError::Unsupported(
                "a motif match that reads the residues of a site beyond their end (undefined in NCBI's rpsbproc) is not supported by LOSAT's rpsbproc".to_string(),
            ),
        }
    }
}

/// `(char)toupper(q.m_strSeqData[res])`: the residue at `res`. `std::string::operator[]` gives
/// the terminating 0 at its size; beyond that, and before 0, NCBI reads bytes outside the
/// string (undefined), which LOSAT refuses.
fn upper_residue(query: &DomQuery, res: i64) -> Result<u8, PrintError> {
    let len = query.residues.len() as i64;
    if (0..len).contains(&res) {
        Ok(query.residues[res as usize].to_ascii_uppercase())
    } else if res == len {
        Ok(0)
    } else {
        Err(PrintError::Unsupported(format!(
            "a site or motif at position {} of a query with {len} residues (NCBI reads outside the sequence) is not supported by LOSAT's rpsbproc",
            res + 1
        )))
    }
}

/// A site of a structural motif type, collected for the MOTIFS part
/// (`__MotifType_base`: the site, the record, the source PSSM-ID).
struct MotifRecord<'a> {
    site: &'a data::Site,
    record: usize,
    source: u32,
}

/// The source column: the PSSM-ID, or the accession of a model with a PSSM-ID of
/// `INVALIDPSSMID` or more that the data has.
fn write_source(report: &mut Report, data: &CddData, source: u32) {
    if source < data::INVALID_PSSM_ID {
        report.write(source.to_string().as_bytes());
    } else if let Some(domain) = data.find_domain(source) {
        report.write(data.text(data.domains[domain].accession));
    } else {
        report.write(source.to_string().as_bytes());
    }
}

/// A site's positions as a set (`CSegSet featsegs(*iterFeat)`: `ori_from` is `from`).
fn site_segset(data: &CddData, site: &data::Site) -> SegSet {
    let mut set = SegSet::new();
    for seg in data.site_segs(site) {
        set.add(i64::from(seg[0]), i64::from(seg[1]));
    }
    set
}

/// `featsegs.GetTotalResidues() / (double)GetCompleteSize() > 0.8` of a mapped site.
fn mostly_mapped(featsegs: &SegSet, site: &data::Site) -> bool {
    !featsegs.is_empty() && (featsegs.total_residues() as f64 / f64::from(site.complete_size)) > 0.8
}

/// The SITES and MOTIFS parts of a protein query (`PrintDomQuery`, `-t feats` and `both`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1033-1266
// ```c++
// 		if (TTargetData::e_doms != m_tdata)
// 		{
// 			vector<TDomSeqAlignIndex::__MotifType_base> vecMotifRecs;
//
// 			vector<TDomSeqAlignIndex::__TCdAlignRecord> :: const_iterator iterFeatAlign = vecFeatAligns.begin(), iterFeatAlignEnd = vecFeatAligns.end();
//
// 			if (iterFeatAlignEnd != iterFeatAlign)	//not empty())
// 			{
// 				if (!iterFeatAlign->pCdInfo->m_bIsStructDom)	//has regular domains
// 				{
// 					(*m_ostream) << FEATSTART << endl;
// 					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
// 					{
// 						if (iterFeatAlign->pCdInfo->m_bIsStructDom) break;
//
// 						bool bIsSpecific = iterFeatAlign->pAlign->m_bRep && iterFeatAlign->pAlign->m_bSpecQualified;
// 						const char * pType = ANNOTTYPE_SPECIFIC;
// 						const list <TDomSite> * pFeats = &(iterFeatAlign->pCdInfo->m_lstSpecificFeatures);
// 						unsigned int uiSrcPssmId = iterFeatAlign->pAlign->m_uiPssmId;
//
// 						if (!bIsSpecific && NULL != iterFeatAlign->pRootCdInfo)
// 						{
// 							pType = ANNOTTYPE_GENERIC;
// 							pFeats = &(iterFeatAlign->pCdInfo->m_lstGenericFeatures);
// 							uiSrcPssmId = iterFeatAlign->pCdInfo->m_uiHierarchyRoot;
// 						}
//
// 						for (list<TDomSite> :: const_iterator iterFeat = pFeats->begin(); iterFeat != pFeats->end(); ++iterFeat)
// 						{
// 							if (TDomSite::eType_StructMotif == iterFeat->m_iType)
// 							{
// 								vecMotifRecs.emplace_back(iterFeat, iterFeatAlign, bIsSpecific, uiSrcPssmId);
// 								continue;
// 							}
//
// 							CSegSet featsegs(*iterFeat);
// 							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
// 							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
// 							{
// 								vector<TSeg_base::TResiduePos> vecMappedPos;
// 								iterFeatAlign->pAlign->GetTranslatedPosMap(featsegs, q.m_uiSeqLen, vecMappedPos);
//
// 								if (iterFeat->MotifCheck(vecMappedPos, q.m_strSeqData) > 0) continue;	//failed motif check
//
// 								(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << DELIMIT << pType << DELIMIT << iterFeat->m_strTitle << DELIMIT;
//
// 								const CSegSet::TSegs& segs = featsegs.GetSegs();
//
// 								size_t mapped = 0;
//
// 								char dimDelimit[2] = {0, 0};
// 								for (CSegSet::TSegs::const_iterator iterSeg = segs.begin(); iterSeg != segs.end(); ++iterSeg)
// 								{
// 									for (int res = iterSeg->from; res <= iterSeg->to; ++res)
// 									{
// 										(*m_ostream) << dimDelimit << (char)toupper(q.m_strSeqData[res]) << (res + COORDSBASE);
// 										dimDelimit[0] = COORDELIMIT;
// 										++mapped;
// 									}
// 								}
// 								(*m_ostream) << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT;
// 								if (uiSrcPssmId < INVALIDPSSMID)
// 									(*m_ostream) << uiSrcPssmId;
// 								else
// 								{
// 									const TDomain * pSrcCd = m_dom_cluster_src.FindCdInfo(uiSrcPssmId);
// 									if (nullptr != pSrcCd)
// 										(*m_ostream) << pSrcCd->m_strAccession;
// 									else
// 										(*m_ostream) << uiSrcPssmId;
// 								}
// 								(*m_ostream) << endl;
// 							}
// 						}
// 					}
//
// 					// -- non-motif features from SD
//
// 					for ( ; iterFeatAlignEnd != iterFeatAlign; ++iterFeatAlign)
// 					{
// 						for (list<TDomSite> :: const_iterator iterFeat = iterFeatAlign->pCdInfo->m_lstSpecificFeatures.begin(); iterFeat != iterFeatAlign->pCdInfo->m_lstSpecificFeatures.end(); ++iterFeat)
// 						{
// 							if (TDomSite::eType_StructMotif == iterFeat->m_iType)
// 							{
// 								vecMotifRecs.emplace_back(iterFeat, iterFeatAlign, true, 0);
// 								continue;
// 							}
//
// 							CSegSet featsegs(*iterFeat);
// 							iterFeatAlign->pAlign->MapSegSet(featsegs, false);
//
// 							if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterFeat->GetCompleteSize() > 0.8))
// 							{
// 								... (the MotifCheck and the line as above, with ANNOTTYPE_SPECIFIC)
// 								(*m_ostream) << DELIMIT << iterFeat->GetTotalResidues() << DELIMIT << mapped << DELIMIT << 0/*iterFeatAlign->pAlign->m_uiPssmId*/ << endl;
// 							}
// 						}
// 					}
//
// 					(*m_ostream) << FEATEND << endl;
// 				}
//
// 				// -- Motif part
// 				if (!vecMotifRecs.empty())
// 				{
// 					// -- start SD
// 					// -- features
// 					(*m_ostream) << MOTIFSTART << endl;
// 					vector<TDomSeqAlignIndex::__MotifType_base> :: const_iterator iterM = vecMotifRecs.begin(), iterMEnd = vecMotifRecs.end();
// 					for ( ; iterMEnd != iterM; ++iterM)
// 					{
// 						if (iterM->iterAlignRec->pCdInfo->m_bIsStructDom) break;
//
// 						CSegSet featsegs(*iterM->iterMotifFeat);
// 						iterM->iterAlignRec->pAlign->MapSegSet(featsegs);
//
// 						if (!featsegs.IsEmpty() && ((double)featsegs.GetTotalResidues() / (double)iterM->iterMotifFeat->GetCompleteSize() > 0.8))
// 						{
// 							(*m_ostream) << idxBlObj << DELIMIT << q.m_strAccession << DELIMIT << iterM->iterMotifFeat->m_strTitle << DELIMIT;
// 							int res0 = featsegs.GetLeft(), res1 = featsegs.GetRight();
// 							(*m_ostream) << (char)toupper(q.m_strSeqData[res0]) << res0 + COORDSBASE << DELIMIT << res1 + COORDSBASE << (char)toupper(q.m_strSeqData[res1]) << DELIMIT;
//
// 							if (iterM->uiSrcPssmId < INVALIDPSSMID)
// 								(*m_ostream) << iterM->uiSrcPssmId;
// 							else
// 							{
// 								const TDomain *pSrcDom = m_dom_cluster_src.FindCdInfo(iterM->uiSrcPssmId);
// 								if (nullptr != pSrcDom)
// 									(*m_ostream) << pSrcDom->m_strAccession;
// 								else
// 									(*m_ostream) << iterM->uiSrcPssmId;
// 							}
//
// 							(*m_ostream) << endl;
// 						}
// 					}
//
// 					// Motif from SD
// 					for ( ; iterMEnd != iterM; ++iterM)
// 					{
// 						CSegSet featsegs(*iterM->iterMotifFeat);
// 						if (!iterM->iterAlignRec->pAlign->m_ClipSet.IsEmpty())
// 							featsegs.Cross(iterM->iterAlignRec->pAlign->m_ClipSet);
//
// 						iterM->iterAlignRec->pAlign->MapSegSet(featsegs);
// 						... (the same line)
// 					}
// 					(*m_ostream) << MOTIFEND << endl;
// 				}
// 			}
// 		}
// ```
// (An `sd` record first in the list prints no SITES part and no motif is collected then.)
fn print_features(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    mode: DataMode,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let records: Vec<FeatureRecord> = feature_records(&query.index, &query.aligns, data, mode);
    let Some(first) = records.first() else {
        return Ok(());
    };
    let mut motifs: Vec<MotifRecord> = Vec::new();
    if !data.domains[first.domain].is_struct_dom {
        report.line(&[b"SITES"]);
        let mut at = 0usize;
        while at < records.len() {
            let record = &records[at];
            if data.domains[record.domain].is_struct_dom {
                break;
            }
            let align = &query.aligns[record.align];
            let is_specific = align.rep && align.spec_qualified;
            let (kind, specific, source): (&[u8], bool, u32) =
                if !is_specific && record.root.is_some() {
                    (
                        b"Generic",
                        false,
                        data.domains[record.domain].hierarchy_root,
                    )
                } else {
                    (b"Specific", true, align.pssm_id)
                };
            for site in data.sites(record.domain, specific) {
                if site.site_type == data::TYPE_STRUCT_MOTIF {
                    motifs.push(MotifRecord {
                        site,
                        record: at,
                        source,
                    });
                    continue;
                }
                print_site(
                    report,
                    query,
                    session,
                    data,
                    align,
                    site,
                    kind,
                    Some(source),
                    motif_state,
                )?;
            }
            at += 1;
        }
        // Sites of the structural domains: source 0.
        while at < records.len() {
            let record = &records[at];
            let align = &query.aligns[record.align];
            for site in data.sites(record.domain, true) {
                if site.site_type == data::TYPE_STRUCT_MOTIF {
                    motifs.push(MotifRecord {
                        site,
                        record: at,
                        source: 0,
                    });
                    continue;
                }
                print_site(
                    report,
                    query,
                    session,
                    data,
                    align,
                    site,
                    b"Specific",
                    None,
                    motif_state,
                )?;
            }
            at += 1;
        }
        report.line(&[b"ENDSITES"]);
    }
    if !motifs.is_empty() {
        report.line(&[b"MOTIFS"]);
        // The first loop stops at the first `sd` record; the second takes the rest.
        let mut structural = false;
        for motif in &motifs {
            let record = &records[motif.record];
            structural = structural || data.domains[record.domain].is_struct_dom;
            let align = &query.aligns[record.align];
            let mut featsegs = site_segset(data, motif.site);
            if structural && !align.clip_set.is_empty() {
                featsegs.cross(&align.clip_set);
            }
            map_seg_set(&mut featsegs, align);
            if mostly_mapped(&featsegs, motif.site) {
                let (res0, res1) = (featsegs.left(), featsegs.right());
                report.write(session.as_bytes());
                report.write(b"\t");
                report.write(&query.accession);
                report.write(b"\t");
                report.write(data.text(motif.site.title));
                report.write(b"\t");
                report.write(&[upper_residue(query, res0)?]);
                report.write((res0 + 1).to_string().as_bytes());
                report.write(b"\t");
                report.write((res1 + 1).to_string().as_bytes());
                report.write(&[upper_residue(query, res1)?]);
                report.write(b"\t");
                write_source(report, data, motif.source);
                report.write(b"\n");
            }
        }
        report.line(&[b"ENDMOTIFS"]);
    }
    Ok(())
}

/// One site line of the SITES part. `source`: the PSSM-ID column (`None`: the structural
/// domains' `0`).
#[allow(clippy::too_many_arguments)]
fn print_site(
    report: &mut Report,
    query: &DomQuery,
    session: &str,
    data: &CddData,
    align: &domains::DomAlign,
    site: &data::Site,
    kind: &[u8],
    source: Option<u32>,
    motif_state: &mut MotifState,
) -> Result<(), PrintError> {
    let mut featsegs = site_segset(data, site);
    map_seg_set(&mut featsegs, align);
    if !mostly_mapped(&featsegs, site) {
        return Ok(());
    }
    let mapped_positions = translated_pos_map(&featsegs);
    if motif_state.motif_check(data, site, &mapped_positions, &query.residues)? > 0 {
        return Ok(());
    }
    report.write(session.as_bytes());
    report.write(b"\t");
    report.write(&query.accession);
    report.write(b"\t");
    report.write(kind);
    report.write(b"\t");
    report.write(data.text(site.title));
    report.write(b"\t");
    let mut mapped = 0usize;
    let mut delimit: &[u8] = b"";
    for seg in featsegs.segs() {
        for res in seg.from..=seg.to {
            report.write(delimit);
            report.write(&[upper_residue(query, res)?]);
            report.write((res + 1).to_string().as_bytes());
            delimit = b",";
            mapped += 1;
        }
    }
    let complete: i64 = data
        .site_segs(site)
        .iter()
        .map(|seg| i64::from(seg[1]) - i64::from(seg[0]) + 1)
        .sum();
    report.write(format!("\t{complete}\t{mapped}\t").as_bytes());
    match source {
        Some(source) => write_source(report, data, source),
        None => report.write(b"0"),
    }
    report.write(b"\n");
    Ok(())
}

/// `N`, `C`, both or `-`: more than 20 % of the model missing at that end.
fn print_incomplete(report: &mut Report, n_missing: f64, c_missing: f64) {
    let mut missing = false;
    if n_missing > 0.2 {
        report.write(b"N");
        missing = true;
    }
    if c_missing > 0.2 {
        report.write(b"C");
        missing = true;
    }
    if !missing {
        report.write(b"-");
    }
}

/// `PrintDomLine`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1807-1832
// ```c++
// void CRpsbProcApp::PrintDomLine(const char *pType, const TDomSeqAlignIndex::__TCdAlignRecord &rec) const
// {
// 	(*m_ostream) << pType << DELIMIT << (rec.pAlign->m_uiPssmId >= INVALIDPSSMID ? 0 : rec.pAlign->m_uiPssmId) << DELIMIT << rec.pAlign->m_iFrom + COORDSBASE << DELIMIT << rec.pAlign->m_iTo + COORDSBASE << DELIMIT << rec.pAlign->m_dEValue << DELIMIT << rec.pAlign->m_dBitScore << DELIMIT << rec.pCdInfo->m_strAccession << DELIMIT << rec.pCdInfo->m_strShortName << DELIMIT;
// 	bool bTmMiss = false;
// 	if (rec.pAlign->m_dNMissing > 0.2)
// 	{
// 		(*m_ostream) << 'N';
// 		bTmMiss = true;
// 	}
// 	if (rec.pAlign->m_dCMissing > 0.2)
// 	{
// 		(*m_ostream) << 'C';
// 		bTmMiss = true;
// 	}
// 	if (!bTmMiss)
// 		(*m_ostream) << '-';
// 	(*m_ostream) << DELIMIT;
// 	if (NULL == rec.pClst)
// 		(*m_ostream) << '-';
// 	else if (rec.pClst->m_uiPssmId >= INVALIDPSSMID)
// 		(*m_ostream) << rec.pClst->m_strAccession;
// 	else
// 		(*m_ostream) << rec.pClst->m_uiPssmId;
// 	//<< (NULL == rec.pClst ? '-' : rec.pClst->m_uiPssmId);
//
// }
// ```
// (`m_iFrom + COORDSBASE` is `int` arithmetic, which overflows at `INT_MAX` (undefined in C++);
// LOSAT adds in `i64`, as `segset::Seg` keeps its sums.)
fn print_dom_line(
    report: &mut Report,
    kind: &[u8],
    record: &Record,
    query: &DomQuery,
    data: &CddData,
) {
    let align = &query.aligns[record.align];
    let domain = &data.domains[record.domain];
    let pssm_id = if align.pssm_id >= data::INVALID_PSSM_ID {
        0
    } else {
        align.pssm_id
    };
    let numbers = format!(
        "\t{pssm_id}\t{}\t{}\t{}\t{}\t",
        i64::from(align.from) + 1,
        i64::from(align.to) + 1,
        format_g(align.e_value),
        format_g(align.bit_score)
    );
    report.write(kind);
    report.write(numbers.as_bytes());
    report.write(data.text(domain.accession));
    report.write(b"\t");
    report.write(data.text(domain.short_name));
    report.write(b"\t");
    print_incomplete(report, align.n_missing, align.c_missing);
    report.write(b"\t");
    match record.cluster {
        None => report.write(b"-"),
        Some(cluster) => {
            let cluster = &data.clusters[cluster];
            if cluster.pssm_id >= data::INVALID_PSSM_ID {
                report.write(data.text(cluster.accession));
            } else {
                report.write(cluster.pssm_id.to_string().as_bytes());
            }
        }
    }
}

/// `PrintClstLine`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1834-1851
// ```c++
// void CRpsbProcApp::PrintClstLine(const TDomSeqAlignIndex::__TCdAlignRecord &rec) const
// {
// 	(*m_ostream) << HITTYPE_CLUSTER << DELIMIT << (rec.pCdInfo->m_uiClusterPssmId >= INVALIDPSSMID ? 0 : rec.pCdInfo->m_uiClusterPssmId) << DELIMIT << rec.pAlign->m_iFrom + COORDSBASE << DELIMIT << rec.pAlign->m_iTo + COORDSBASE << DELIMIT << rec.pAlign->m_dEValue << DELIMIT << rec.pAlign->m_dBitScore << DELIMIT << rec.pClst->m_strAccession << DELIMIT << rec.pClst->m_strShortName << DELIMIT;
// 	bool bTmMiss = false;
// 	...
// 	(*m_ostream) << DELIMIT << '-';
// }
// ```
// (`m_uiClusterPssmId` is 0: the offline loader never sets it.)
fn print_cluster_line(report: &mut Report, record: &Record, query: &DomQuery, data: &CddData) {
    let align = &query.aligns[record.align];
    let cluster = &data.clusters[record.cluster.expect("a superfamily line has its cluster")];
    let numbers = format!(
        "Superfamily\t0\t{}\t{}\t{}\t{}\t",
        i64::from(align.from) + 1,
        i64::from(align.to) + 1,
        format_g(align.e_value),
        format_g(align.bit_score)
    );
    report.write(numbers.as_bytes());
    report.write(data.text(cluster.accession));
    report.write(b"\t");
    report.write(data.text(cluster.short_name));
    report.write(b"\t");
    print_incomplete(report, align.n_missing, align.c_missing);
    report.write(b"\t-");
}

#[cfg(test)]
mod tests {
    use super::*;

    // ncbifile.cpp:820-971 (quoted at `ncbi_environment::normalize_path`).
    #[test]
    fn normalize_path_as_ncbi() {
        for (path, normalized) in [
            ("", ""),
            ("/", "/"),
            ("///", "/"),
            ("losat", "losat"),
            ("./losat", "losat"),
            ("bin/./losat", "bin/losat"),
            ("bin//losat", "bin/losat"),
            ("bin/../losat", "losat"),
            ("../bin/losat", "../bin/losat"),
            ("/../bin/losat", "/bin/losat"),
            ("bin/losat/", "bin/losat"),
            (".", "."),
            ("./", "."),
        ] {
            assert_eq!(
                normalize_path(path.as_bytes()),
                normalized.as_bytes(),
                "{path}"
            );
        }
    }

    // `ostream << double` of -0.0 is glibc's `%g`, "-0".
    #[test]
    fn format_g_keeps_the_sign_of_zero() {
        assert_eq!(format_g(-0.0), "-0");
        assert_eq!(format_g(0.0), "0");
        assert_eq!(format_g(3.06789e-55), "3.06789e-55");
        assert_eq!(format_g(179.393), "179.393");
        assert_eq!(OVERLAPLEADING, 1_410_065_408);
    }
}
