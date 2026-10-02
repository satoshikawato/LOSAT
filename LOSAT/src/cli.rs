//! Public CLI v2 boundary. Clap's double-dash representation is internal only.

use std::{ffi::OsString, fmt};

use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
// ```c++
//     static const string kProgram("blastx");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Translated Query-Protein Subject BLAST"));
// ```
use crate::algorithm::{blastn, blastp, blastx, tblastn, tblastx};

// NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
// ```c++
// const string kArgQuery("query");
// const string kArgSubject("subject");
// const string kArgNumThreads("num_threads");
// ```
// LOSAT retains its program subcommand; search parameter names are NCBI names.
#[derive(Parser, Debug)]
#[command(
    name = "losat",
    version,
    about = "Pure-Rust pairwise alignment with NCBI-verified fixtures and genetic-code exceptions"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
    // ```c++
    //     static const string kProgram("blastx");
    //     arg.Reset(new CProgramDescriptionArgs(kProgram,
    //                                   "Translated Query-Protein Subject BLAST"));
    // ```
    /// Translated nucleotide query vs protein subject (local blastx)
    Blastx(blastx::BlastxArgs),
    /// Pairwise nucleotide alignment (megablast [default], blastn)
    Blastn(blastn::BlastnArgs),
    /// Pairwise protein alignment (blastp)
    Blastp(blastp::BlastpArgs),
    /// Pairwise 6-frame translated nucleotide alignment (tblastx)
    Tblastx(tblastx::TblastxArgs),
    // NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp:47-52:
    // static const string kProgram("tblastn"); SetTask("tblastn");
    /// Protein query vs translated nucleotide subject (local tblastn)
    Tblastn(tblastn::TblastnArgs),
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-170,203-207,332-349
// ```c++
// arg_desc.AddOptionalKey(kArgWordSize, "int_value", description, CArgDescriptions::eInteger);
// arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value", ...);
// arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options", ..., CArgDescriptions::eString, ...);
// ```
// Translate registered canonical keys only. Consume each value once so a path
// or a quoted filtering specification can never be interpreted as another key.
pub fn try_parse_from<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: Parser + CommandFactory,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let command = T::command();
    let mut input = argv.into_iter().map(Into::into);
    let mut translated = vec![input.next().unwrap_or_else(|| "losat".into())];
    let mut scope = &command;
    while let Some(token) = input.next() {
        let text = token.to_str().unwrap_or("");
        if let Some(subcommand) = scope.find_subcommand(text) {
            scope = subcommand;
            translated.push(token);
            continue;
        }
        if matches!(text, "-help" | "--help") {
            translated.push("--help".into());
            continue;
        }
        if text == "--version" && std::ptr::eq(scope, &command) {
            translated.push(token);
            continue;
        }
        let key = text.strip_prefix('-').filter(|key| !key.starts_with('-'));
        let (name, inline) = key
            .unwrap_or("")
            .split_once('=')
            .map_or((key.unwrap_or(""), None), |(k, v)| (k, Some(v)));
        let Some(arg) = scope
            .get_arguments()
            .find(|arg| arg.get_long() == Some(name))
        else {
            // NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp:64-129:
            // m_Args.push_back(arg) registers shared, formatting, DB, and PSI groups.
            // The named NCBI options below have no Stage B Rust behavior.
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
            // ```c++
            //     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
            //     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
            //     m_BlastDbArgs->SetIPGFilteringSupport(true);
            //     arg.Reset(m_BlastDbArgs);
            //     m_Args.push_back(arg);
            // ```
            // Product scope explicitly refuses unported NCBI capabilities.
            if scope.get_name() == "blastx" && is_unported_blastx_arg(name) {
                return Err(clap::Error::raw(ErrorKind::InvalidValue,
                    format!("unsupported BLASTX option '-{name}': outside the declared local FASTA scope")));
            }
            if scope.get_name() == "blastn" && is_unported_blastn_arg(name) {
                return Err(clap::Error::raw(
                    ErrorKind::InvalidValue,
                    format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's BLASTN"),
                ));
            }
            if scope.get_name() == "tblastn" && is_unported_tblastn_arg(name) {
                return Err(clap::Error::raw(
                    ErrorKind::InvalidValue,
                    format!("unsupported TBLASTN option '-{name}': Rust behavior is unimplemented"),
                ));
            }
            return Err(clap::Error::raw(
                ErrorKind::UnknownArgument,
                format!("unknown option or argument '{text}'; use -help for CLI v2 syntax"),
            ));
        };
        if arg.get_action().takes_values() {
            let value = match inline {
                Some(value) => OsString::from(value),
                None => input.next().ok_or_else(|| {
                    clap::Error::raw(
                        ErrorKind::InvalidValue,
                        format!("-{name} requires one value"),
                    )
                })?,
            };
            let mut internal = OsString::from(format!("--{name}="));
            internal.push(value);
            translated.push(internal);
        } else {
            if inline.is_some() {
                return Err(clap::Error::raw(
                    ErrorKind::TooManyValues,
                    format!("-{name} is a flag and does not take a value"),
                ));
            }
            translated.push(format!("--{name}").into());
        }
    }
    T::try_parse_from(translated)
}

// NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
// ```c++
// const string kArgWordSize("word_size");
// const string kArgCompBasedStats("comp_based_stats");
// ```
// Render the same canonical names in help, usage and parser diagnostics.
pub fn render_message(error: &clap::Error) -> String {
    let mut message = error
        .to_string()
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    let command = Cli::command();
    for scope in std::iter::once(&command).chain(command.get_subcommands()) {
        for arg in scope.get_arguments() {
            if let Some(name) = arg
                .get_long()
                .filter(|name| !matches!(*name, "help" | "version"))
            {
                message = message.replace(&format!("--{name}"), &format!("-{name}"));
            }
        }
    }
    message
        .replace("'-h'", "'-help'")
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version")
}

// NCBI reference: c++/src/algo/blast/blastinput/tblastn_args.cpp:64-129
// ```c++
// m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
// arg.Reset(new CGenericSearchArgs(kQueryIsProtein));
// m_HspFilteringArgs.Reset(new CHspFilteringArgs);
// m_FormattingArgs.Reset(new CFormattingArgs);
// m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eNucleotideDb));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
// ```c++
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
//     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
//     arg.Reset(m_BlastDbArgs);
//     m_Args.push_back(arg);
//
//     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
//     arg.Reset(m_StdCmdLineArgs);
//     m_Args.push_back(arg);
// ```
// The options of NCBI blastn 2.17.0+ (-help) that LOSAT's BLASTN does not implement
// (AGENTS.md rule 2: explicit unsupported errors).
fn is_unported_blastn_arg(name: &str) -> bool {
    matches!(
        name,
        "best_hit_overhang"
            | "best_hit_score_edge"
            | "culling_limit"
            | "db"
            | "db_hard_mask"
            | "db_soft_mask"
            | "dbsize"
            | "entrez_query"
            | "export_search_strategy"
            | "filtering_db"
            | "gilist"
            | "h"
            | "html"
            | "import_search_strategy"
            | "index_name"
            | "line_length"
            | "min_raw_gapped_score"
            | "mt_mode"
            | "negative_gilist"
            | "negative_seqidlist"
            | "negative_taxidlist"
            | "negative_taxids"
            | "no_greedy"
            | "no_taxid_expansion"
            | "num_alignments"
            | "num_descriptions"
            | "off_diagonal_range"
            | "parse_deflines"
            | "qcov_hsp_perc"
            | "query_loc"
            | "remote"
            | "searchsp"
            | "seqidlist"
            | "show_gis"
            | "soft_masking"
            | "sorthits"
            | "sorthsps"
            | "strand"
            | "subject_loc"
            | "taxidlist"
            | "taxids"
            | "template_length"
            | "template_type"
            | "ungapped"
            | "use_index"
            | "version"
            | "window_masker_db"
            | "window_masker_taxid"
            | "window_size"
            | "xdrop_gap"
            | "xdrop_gap_final"
            | "xdrop_ungap"
    )
}

// Names are from the pinned 2.17.0+ -help and have no implemented Rust path yet.
fn is_unported_tblastn_arg(name: &str) -> bool {
    matches!(
        name,
        "query_loc"
            | "show_gis"
            | "num_descriptions"
            | "num_alignments"
            | "line_length"
            | "html"
            | "sorthits"
            | "sorthsps"
            | "lcase_masking"
            | "gilist"
            | "seqidlist"
            | "negative_gilist"
            | "negative_seqidlist"
            | "taxids"
            | "negative_taxids"
            | "taxidlist"
            | "negative_taxidlist"
            | "no_taxid_expansion"
            | "entrez_query"
            | "db_soft_mask"
            | "db_hard_mask"
            | "qcov_hsp_perc"
            | "max_hsps"
            | "culling_limit"
            | "best_hit_overhang"
            | "best_hit_score_edge"
            | "subject_besthit"
            | "dbsize"
            | "searchsp"
            | "import_search_strategy"
            | "export_search_strategy"
            | "xdrop_ungap"
            | "xdrop_gap"
            | "xdrop_gap_final"
            | "parse_deflines"
            | "mt_mode"
            | "use_sw_tback"
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
// ```c++
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
//     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
//     m_BlastDbArgs->SetIPGFilteringSupport(true);
//     arg.Reset(m_BlastDbArgs);
//     m_Args.push_back(arg);
// ```
// Named capabilities exist in NCBI but are explicitly outside the BLASTX product scope.
fn is_unported_blastx_arg(name: &str) -> bool {
    matches!(
        name,
        "db" | "query_loc"
            | "subject_loc"
            | "show_gis"
            | "num_descriptions"
            | "num_alignments"
            | "line_length"
            | "html"
            | "sorthits"
            | "sorthsps"
            | "gilist"
            | "seqidlist"
            | "negative_gilist"
            | "negative_seqidlist"
            | "taxids"
            | "negative_taxids"
            | "taxidlist"
            | "negative_taxidlist"
            | "no_taxid_expansion"
            | "entrez_query"
            | "db_soft_mask"
            | "db_hard_mask"
            | "ipglist"
            | "negative_ipglist"
            | "qcov_hsp_perc"
            | "best_hit_overhang"
            | "best_hit_score_edge"
            | "dbsize"
            | "searchsp"
            | "import_search_strategy"
            | "export_search_strategy"
            | "xdrop_ungap"
            | "xdrop_gap"
            | "xdrop_gap_final"
            | "parse_deflines"
            | "mt_mode"
            | "remote"
            | "use_sw_tback"
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
// ```c++
//     catch (const CArgException& e) {                                        \
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-230
// ```c++
//     catch (const blast::CBlastException& e) {                               \
//         const string& msg = e.GetMsg();                                     \
//         if (e.GetErrCode() == CBlastException::eInvalidOptions) {           \
//             LOG_POST(Error << "BLAST options error: " << e.GetMsg());       \
//             exit_code = BLAST_INPUT_ERROR;                                  \
//         } else if ((NStr::Find(msg, "Out of memory") != NPOS) ||            \
//             (NStr::Find(msg, "Failed to allocate") != NPOS)) {              \
//             LOG_POST(Error << "BLAST ran out of memory: " << e.GetMsg());   \
//             exit_code = BLAST_OUT_OF_MEMORY;                                \
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
//         }                                                                   \
//     }                                                                       \
//     catch (const blast::CBlastSystemException& e) {                         \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// ```c++
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
#[derive(Debug)]
pub struct NativeError {
    pub exit: i32,
    pub message: String,
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:178-180
// ```c++
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
impl fmt::Display for NativeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
// ```c++
//     catch (const CArgException& e) {                                        \
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
impl std::error::Error for NativeError {}

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:118-119
// ```c
// {
//     m_Outfile.exceptions(NcbiBadbit);
// ```
/// The report's output stream (the `-out` file or standard output), which records whether
/// a write to it failed: NCBI's outfmt 0 stream throws there, which the application
/// reports as "BLAST failed to write output" (blast_app_util.hpp:252-255).
pub struct ReportStream {
    pub inner: Box<dyn std::io::Write + Send>,
    pub failed: bool,
}

impl std::io::Write for ReportStream {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf);
        self.failed |= written.is_err();
        written
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let flushed = self.inner.flush();
        self.failed |= flushed.is_err();
        flushed
    }
}

/// Standard output as the report's stream (`ReportStream::inner`).
///
/// NCBI writes the report to `cout`; when the process starts with its standard output
/// closed (`>&-`), the first write fails (outfmt 0: "BLAST failed to write output", exit
/// 6; outfmt 6 and 7 abort). Rust's runtime opens `/dev/null` on a closed standard
/// descriptor before `main`, read and write, so the writes would succeed: on Linux, that
/// descriptor is recognised (`/proc/self/fdinfo/1`; a shell's `> /dev/null` opens it write
/// only) and every write fails, as the approved CLI differences expect
/// (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES).
pub fn report_standard_output() -> Box<dyn std::io::Write + Send> {
    if standard_output_was_closed() {
        Box::new(ClosedStandardOutput)
    } else {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    }
}

/// Whether standard output is the `/dev/null` that the runtime opened read and write in
/// place of a closed descriptor (O_ACCMODE is 3 and O_RDWR is 2 on Linux).
#[cfg(target_os = "linux")]
fn standard_output_was_closed() -> bool {
    let is_dev_null = std::fs::read_link("/proc/self/fd/1")
        .is_ok_and(|target| target == std::path::Path::new("/dev/null"));
    is_dev_null
        && std::fs::read_to_string("/proc/self/fdinfo/1")
            .ok()
            .and_then(|info| {
                info.lines()
                    .find_map(|line| line.strip_prefix("flags:"))
                    .and_then(|flags| u32::from_str_radix(flags.trim(), 8).ok())
            })
            .is_some_and(|flags| flags & 3 == 2)
}

#[cfg(not(target_os = "linux"))]
fn standard_output_was_closed() -> bool {
    false
}

/// A standard output that was closed when the process started: every write fails.
struct ClosedStandardOutput;

impl std::io::Write for ClosedStandardOutput {
    fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("standard output is closed"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
// ```c++
// string s_ArgExptMsg(const string& name, const string& what, const string& attr)
// {
//     return string("Argument \"") + (name.empty() ? s_ExtraName : name) +
//         "\". " + what + (attr.empty() ? attr : ":  `" + attr + "'");
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
// ```c++
// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
// {
//     if ( !m_Ios ) {
//         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
//             "File is not accessible",AsString()));
// ```
pub fn inaccessible(name: &str, path: &std::path::Path) -> anyhow::Error {
    let value = if path.as_os_str().is_empty() {
        String::new()
    } else {
        format!(":  `{}'", path.display())
    };
    NativeError {
        exit: 1,
        message: format!(
            "Command line argument error: Argument \"{name}\". File is not accessible{value}\n"
        ),
    }
    .into()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
// ```c++
//     catch (const CArgException& e) {                                        \
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// ```c++
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
// ```c++
// CBlastTabularInfo::~CBlastTabularInfo()
// {
//     m_Ostream.flush();
// }
// ```
pub fn exit_on_native_error(error: &anyhow::Error) {
    if let Some(e) = error.downcast_ref::<NativeError>() {
        eprint!("{}", e.message);
        if e.exit < 0 {
            std::process::abort();
        }
        std::process::exit(e.exit);
    }
}
