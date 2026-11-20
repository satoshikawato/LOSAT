//! `losat windowmasker`: NCBI BLAST+ 2.17.0 windowmasker 1.0.0 (`CWinMaskApplication`), the
//! masker of over-represented and low-complexity regions of nucleotide sequences by the counts
//! of their units (WindowMasker, `utils::winmask`) and, with `-dust`, the symmetric DUST
//! algorithm. Its command line is NCBI's (`ncbi_app`), not clap's: `main` hands it the words
//! after `windowmasker`. LOSAT masks with a counts file (`-ustat`); the making of counts
//! (`-mk_counts`) and their conversion (`-convert`) are refused. The work of its `Run` is in
//! `app` (the configuration, `CWinMaskConfig`) and `work` (the masking).

mod app;
mod ids;
mod sdust;
mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, Constraint, Dep, Desc, Kind};
use crate::ncbi_app::Tool;
use crate::ncbi_version::Component;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "WINDOWMASKER";

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:76,99-110
// ```c++
// const char * const CWinMaskApplication::USAGE_LINE
//     = "Window based sequence masker";
// ...
// void CWinMaskApplication::Init(void)
// {
//     HideStdArgs(fHideLogfile | fHideConffile | fHideVersion | fHideDryRun);
//     unique_ptr< CArgDescriptions > arg_desc( new CArgDescriptions );
//
//     // Set the program description
//     arg_desc->SetUsageContext( GetArguments().GetProgramBasename(),
//                                USAGE_LINE );
//
//     CWinMaskConfig::AddWinMaskArgs(*arg_desc);
// ```
/// The text of `windowmasker -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions in `DESCS` after `HideStdArgs`), frozen. Its part before `OPTIONAL ARGUMENTS` is
/// the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

const INT_MAX: i64 = i32::MAX as i64;

// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:50-213
// ```c++
//     arg_desc.SetCurrentGroup("WindowMasker options");
//
//     // Adding command line arguments descriptions
//     if(type == eAny)
//         arg_desc.AddOptionalKey( "ustat", "unit_counts",
//                                  "file with unit counts",
//                                  CArgDescriptions::eString);
// ...
//     if(determine_input)
//         arg_desc.AddDefaultKey( kInput, "input_file_name", ...,
//                                  CArgDescriptions::eInputFile, "-" );
//     arg_desc.AddDefaultKey( kOutput, "output_file_name", ...,
//                              CArgDescriptions::eOutputFile, "-" );
//     if(type == eAny || type == eComputeCounts){
//         arg_desc.AddDefaultKey( "checkdup", ..., CArgDescriptions::eBoolean, "false" );
//         if(determine_input)
//             arg_desc.AddDefaultKey( "fa_list", ..., CArgDescriptions::eBoolean, "false" );
//         arg_desc.AddDefaultKey( "mem", ..., CArgDescriptions::eInteger, "1536" );
//         arg_desc.AddDefaultKey( "meta", ..., CArgDescriptions::eString, "" );
//         arg_desc.AddOptionalKey( "unit", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "genome_size", ..., CArgDescriptions::eInteger );
//         arg_desc.SetConstraint( "mem", new CArgAllow_Integers( 1, kMax_Int ) );
//         arg_desc.SetConstraint( "unit", new CArgAllow_Integers( 1, 16 ) );
//     }
//     if(type == eAny || type >= eGenerateMasks){
//         arg_desc.AddOptionalKey( "window", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "t_extend", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "t_thres", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "set_t_high", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "set_t_low", ..., CArgDescriptions::eInteger );
//         arg_desc.SetConstraint( "window", new CArgAllow_Integers( 1, kMax_Int ) );
//         arg_desc.SetConstraint( "t_extend", new CArgAllow_Integers( 0, kMax_Int ) );
//         arg_desc.SetConstraint( "t_thres", new CArgAllow_Integers( 1, kMax_Int ) );
//         arg_desc.SetConstraint( "set_t_high", new CArgAllow_Integers( 1, kMax_Int ) );
//         arg_desc.SetConstraint( "set_t_low", new CArgAllow_Integers( 1, kMax_Int ) );
//         arg_desc.AddFlag      ( "parse_seqids", "Parse Seq-ids in FASTA input", true );
//         arg_desc.AddDefaultKey( kOutputFormat, ..., CArgDescriptions::eString, *kOutputFormats );
//         ...
//         arg_desc.SetConstraint( kOutputFormat, strings_allowed );
//     }
//     if(type != eConvertCounts){
//         arg_desc.AddOptionalKey( "t_high", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "t_high_pct", ..., CArgDescriptions::eDouble );
//         arg_desc.AddOptionalKey( "t_thres_pct", ..., CArgDescriptions::eDouble );
//         arg_desc.AddOptionalKey( "t_extend_pct", ..., CArgDescriptions::eDouble );
//         arg_desc.AddOptionalKey( "t_low", ..., CArgDescriptions::eInteger );
//         arg_desc.AddOptionalKey( "t_low_pct", ..., CArgDescriptions::eDouble );
//         arg_desc.SetConstraint( "t_high", new CArgAllow_Integers( 1, kMax_Int ) );
//         arg_desc.SetConstraint( "t_high_pct", new CArgAllow_Doubles( 0.0, 100.0 ) );
//         ... (the same for t_thres_pct, t_extend_pct; t_low as t_high; t_low_pct as t_high_pct)
//         arg_desc.AddDefaultKey( kInputFormat, ..., CArgDescriptions::eString, *kInputFormats );
//         arg_desc.AddDefaultKey( "exclude_ids", ..., CArgDescriptions::eString, "" );
//         arg_desc.AddDefaultKey( "ids", ..., CArgDescriptions::eString, "" );
//         arg_desc.AddDefaultKey( "text_match", ..., CArgDescriptions::eBoolean, "T" );
//         ...
//         strings_allowed->Allow("seqids");
//         arg_desc.SetConstraint( kInputFormat, strings_allowed );
//     }
//     if(type < eGenerateMasks){
//         arg_desc.AddDefaultKey( "sformat", ..., CArgDescriptions::eString, "ascii" );
//         arg_desc.SetConstraint( "sformat", (new CArgAllow_Strings())
//                                  ->Allow( "ascii" )->Allow( "binary" )
//                                  ->Allow( "oascii" )->Allow( "obinary" ) );
//         arg_desc.AddDefaultKey( "smem", ..., CArgDescriptions::eInteger, "512" );
//     }
//     if(type == eAny || type >= eGenerateMasksWithDuster) {
//         arg_desc.AddDefaultKey( "dust", ..., CArgDescriptions::eBoolean, "F" );
//         arg_desc.AddDefaultKey( "dust_level", ..., CArgDescriptions::eInteger, "20" );
//     }
//
//     if(type == eAny){
//         arg_desc.AddFlag( "mk_counts", "generate frequency counts for a database" );
//         arg_desc.AddFlag( "convert", "convert counts between different formats" );
// ```
// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_cmdline_args.cpp:40-51
// ```c++
// const std::string kInput("in");
// const std::string kInputFormat("infmt");
// const std::string kOutput("out");
// const std::string kOutputFormat("outfmt");
// const size_t kNumInputFormats = 2;
// const size_t kNumOutputFormats = 8;
// const char* kInputFormats[] =  { "fasta", "blastdb" };
// const char* kOutputFormats[] =
//     { "interval", "fasta",
//       "seqloc_asn1_bin", "seqloc_asn1_text", "seqloc_xml",
//       "maskinfo_asn1_bin", "maskinfo_asn1_text", "maskinfo_xml",
//     };
// ```
/// The arguments of `CWinMaskConfig::AddWinMaskArgs(arg_desc)` (`eAny`, the input determined)
/// in the order it describes them; NCBI's reader adds the help flags before them and `-logfile`
/// and `-conffile` after them (`ncbi_app::args`). The values that a `CArgAllow_Strings` allows
/// are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("ustat", Kind::Text, None, None),
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key("checkdup", Kind::Boolean, Some("false"), None),
    key("fa_list", Kind::Boolean, Some("false"), None),
    key(
        "mem",
        Kind::Integer,
        Some("1536"),
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("meta", Kind::Text, Some(""), None),
    key("unit", Kind::Integer, None, Some(Constraint::Integers(1, 16))),
    key("genome_size", Kind::Integer, None, None),
    key(
        "window",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "t_extend",
        Kind::Integer,
        None,
        Some(Constraint::Integers(0, INT_MAX)),
    ),
    key(
        "t_thres",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "set_t_high",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "set_t_low",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    flag("parse_seqids"),
    key(
        "outfmt",
        Kind::Text,
        Some("interval"),
        Some(Constraint::Strings(&[
            "fasta",
            "interval",
            "maskinfo_asn1_bin",
            "maskinfo_asn1_text",
            "maskinfo_xml",
            "seqloc_asn1_bin",
            "seqloc_asn1_text",
            "seqloc_xml",
        ])),
    ),
    key(
        "t_high",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("t_high_pct", Kind::Double, None, Some(PCT)),
    key("t_thres_pct", Kind::Double, None, Some(PCT)),
    key("t_extend_pct", Kind::Double, None, Some(PCT)),
    key(
        "t_low",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("t_low_pct", Kind::Double, None, Some(PCT)),
    key(
        "infmt",
        Kind::Text,
        Some("fasta"),
        Some(Constraint::Strings(&["blastdb", "fasta", "seqids"])),
    ),
    key("exclude_ids", Kind::Text, Some(""), None),
    key("ids", Kind::Text, Some(""), None),
    key("text_match", Kind::Boolean, Some("T"), None),
    key(
        "sformat",
        Kind::Text,
        Some("ascii"),
        Some(Constraint::Strings(&["ascii", "binary", "oascii", "obinary"])),
    ),
    key("smem", Kind::Integer, Some("512"), None),
    key("dust", Kind::Boolean, Some("F"), None),
    key("dust_level", Kind::Integer, Some("20"), None),
    flag("mk_counts"),
    flag("convert"),
];

/// `CArgAllow_Doubles( 0.0, 100.0 )` of the `-t_*_pct` keys.
const PCT: Constraint = Constraint::Doubles(0.0, 100.0, "0..100");

// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:212-254
// ```c++
//         arg_desc.CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "outfmt" );
//         ... (in the order below)
// ```
/// The `SetDependency` calls of `AddWinMaskArgs`, in their order.
fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    let excludes = |of: &'static str, names: &[&'static str]| {
        names
            .iter()
            .map(move |&name| (of, Dep::Excludes, name))
            .collect::<Vec<_>>()
    };
    let mut list = excludes(
        "mk_counts",
        &[
            "outfmt",
            "ustat",
            "window",
            "t_thres",
            "t_extend",
            "set_t_low",
            "set_t_high",
            "dust",
            "dust_level",
            "convert",
        ],
    );
    list.extend(excludes(
        "ustat",
        &[
            "checkdup",
            "fa_list",
            "mem",
            "unit",
            "genome_size",
            "sformat",
            "smem",
            "convert",
        ],
    ));
    list.extend(excludes(
        "convert",
        &[
            "checkdup",
            "window",
            "t_extend",
            "t_thres",
            "t_high",
            "t_low",
            "set_t_low",
            "set_t_high",
            "infmt",
            "outfmt",
            "parse_seqids",
            "fa_list",
            "mem",
            "unit",
            "genome_size",
            "dust",
            "dust_level",
            "exclude_ids",
            "ids",
            "text_match",
        ],
    ));
    list
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:79-84
// ```c++
// CWinMaskApplication::CWinMaskApplication() {
//     CRef<CVersion> version(new CVersion());
//     version->SetVersionInfo(1, 0, 0);
//     ...
//     SetFullVersion(version);
// }
// ```
/// windowmasker as NCBI's application runs it (`ncbi_app`): its own version and no prefix on
/// its diagnostics.
static TOOL: Tool = Tool {
    program: "windowmasker",
    refusal_name: PROGRAM,
    version: "1.0.0",
    help: HELP,
    descs: DESCS,
    dependencies,
    post_prefix: None,
    components: COMPONENTS,
};

/// The components that `CWinMaskApplication`'s constructor adds: the algorithm, the counts
/// algorithm and the four counts formats, which LOSAT's masker and counts readers implement.
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:84-95
/// ```c++
///     version->AddComponentVersion( new CSeqMaskerVersion(
///                 CSeqMasker::AlgoVersion ) );
///     version->AddComponentVersion( new CSeqMaskerVersion(
///                 CSeqMaskerOstat::StatAlgoVersion ) );
///     version->AddComponentVersion( new CSeqMaskerVersion(
///                 CSeqMaskerOstatAscii::FormatVersion ) );
///     version->AddComponentVersion( new CSeqMaskerVersion(
///                 CSeqMaskerOstatBin::FormatVersion ) );
///     version->AddComponentVersion( new CSeqMaskerVersion(
///                 CSeqMaskerOstatOptAscii::FormatVersion ) );
///     version->AddComponentVersion( new CSeqMaskerVersion(
///                 CSeqMaskerOstatOptBin::FormatVersion ) );
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:57-60,
/// seq_masker_ostat.cpp:41-47 (`STAT_ALGO_*`), seq_masker_ostat_opt_bin.cpp:39-43 (and the other
/// formats' `STAT_FMT_*`)
/// ```c++
/// #define WIN_MASK_ALGO_NAME "window-masker-algorithm"
/// #define WIN_MASK_ALGO_VER_MAJOR 1
/// #define WIN_MASK_ALGO_VER_MINOR 0
/// #define WIN_MASK_ALGO_VER_PATCH 0
/// ...
/// #define STAT_FMT_COMPONENT_NAME "windowmasker-statistics-format-version"
/// #define STAT_FMT_VER_MAJOR 1
/// #define STAT_FMT_VER_MINOR 1
/// #define STAT_FMT_VER_PATCH 0
/// #define STAT_FMT_VER_PFX "obinary "
/// ```
const COMPONENTS: &[Component] = &[
    component("window-masker-algorithm", "", 0),
    component("windowmasker-statistics-algorithm-version", "", 0),
    component("windowmasker-statistics-format-version", "ascii ", 1),
    component("windowmasker-statistics-format-version", "binary ", 1),
    component("windowmasker-statistics-format-version", "oascii ", 1),
    component("windowmasker-statistics-format-version", "obinary ", 1),
];

const fn component(name: &'static str, prefix: &'static str, minor: u32) -> Component {
    Component {
        name,
        prefix,
        major: 1,
        minor,
        patch: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/main.cpp
// ```c++
//     return CWinMaskApplication().AppMain(argc, argv);
// ```
/// Runs `windowmasker` with `words` (the command line after `losat windowmasker`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, app::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around windowmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_windowmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, app::run)
}

#[cfg(test)]
mod tests;
