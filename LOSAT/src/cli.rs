//! Public CLI v2 boundary. Clap's double-dash representation is internal only.

use std::ffi::OsString;

use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};

use crate::algorithm::{blastn, blastp, tblastn, tblastx};

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
    about = "Pure-Rust pairwise alignment with NCBI BLAST+ parity on certified fixtures"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
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
