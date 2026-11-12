//! `losat dustmasker`: NCBI BLAST+ 2.17.0 dustmasker 1.0.0 (`CDustMaskApplication`), the masker
//! of low-complexity regions of nucleotide sequences by the symmetric DUST algorithm. Its
//! command line is NCBI's (`ncbi_app`, NCBI's application and argument reader), not clap's:
//! `main` hands it the words after `dustmasker`. The work of its `Run` is in `app`.

mod app;
mod fasta;
mod masker;
mod work;

use std::ffi::OsString;
use std::io::Write;

use crate::ncbi_app::args::{flag, key, Constraint, Desc, Kind};
use crate::ncbi_app::{no_dependencies, Tool};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "DUSTMASKER";

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:76-85
// ```c++
// const char * const CDustMaskApplication::USAGE_LINE
//     = "Low complexity region masker based on Symmetric DUST algorithm";
//
// //-------------------------------------------------------------------------
// void CDustMaskApplication::Init(void)
// {
//     HideStdArgs(fHideLogfile | fHideConffile | fHideVersion | fHideDryRun);
//     unique_ptr< CArgDescriptions > arg_desc( new CArgDescriptions );
//     arg_desc->SetUsageContext( GetArguments().GetProgramBasename(),
//                                USAGE_LINE );
// ```
/// The text of `dustmasker -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions in `DESCS` after `HideStdArgs`: `-logfile`, `-conffile`, `-version` and `-dryrun`
/// are not shown). Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the argument
/// errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:86-117
// ```c++
//     arg_desc->AddDefaultKey( kInput, "input_file_name",
//                              "input file name",
//                              CArgDescriptions::eInputFile, "-" );
//     arg_desc->AddDefaultKey( kOutput, "output_file_name",
//                              "output file name",
//                              CArgDescriptions::eOutputFile, "-");
//     arg_desc->AddDefaultKey( "window", "window_size",
//                              "DUST window length",
//                              CArgDescriptions::eInteger, "64" );
//     arg_desc->AddDefaultKey( "level", "level",
//                              "DUST level (score threshold for subwindows)",
//                              CArgDescriptions::eInteger, "20" );
//     arg_desc->AddDefaultKey( "linker", "linker",
//                              "DUST linker (how close masked intervals "
//                              "should be to get merged together).",
//                              CArgDescriptions::eInteger, "1" );
//     arg_desc->AddDefaultKey( kInputFormat, "input_format",
//                              "input format (possible values: fasta, blastdb)",
//                              CArgDescriptions::eString, *kInputFormats );
//     arg_desc->AddDefaultKey( kOutputFormat, "output_format",
//                              "output format",
//                              CArgDescriptions::eString, *kOutputFormats );
//     arg_desc->AddFlag      ( "parse_seqids",
//                              "Parse Seq-ids in FASTA input", true );
//     arg_desc->AddFlag      ( "hard_masking",
//                              "Use hard masking for fasta outfmt", true );
//     CArgAllow_Strings* strings_allowed = new CArgAllow_Strings();
//     for (size_t i = 0; i < kNumOutputFormats; i++) {
//         strings_allowed->Allow(kOutputFormats[i]);
//     }
//     strings_allowed->Allow("acclist");
//     arg_desc->SetConstraint( kOutputFormat, strings_allowed );
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
/// The arguments of `CDustMaskApplication::Init` in the order it describes them; NCBI's reader
/// adds the help flags before them and `-logfile` and `-conffile` after them (`ncbi_app::args`).
/// The values that `-outfmt` allows are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key("window", Kind::Integer, Some("64"), None),
    key("level", Kind::Integer, Some("20"), None),
    key("linker", Kind::Integer, Some("1"), None),
    key("infmt", Kind::Text, Some("fasta"), None),
    key(
        "outfmt",
        Kind::Text,
        Some("interval"),
        Some(Constraint::Strings(&[
            "acclist",
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
    flag("parse_seqids"),
    flag("hard_masking"),
];

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.hpp:51-55
// ```c++
//     CDustMaskApplication() {
//         CRef<CVersion> version(new CVersion());
//         version->SetVersionInfo(1, 0, 0);
//         SetFullVersion(version);
//     }
// ```
/// dustmasker as NCBI's application runs it (`ncbi_app`): its own version, no dependencies
/// between its arguments and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "dustmasker",
    refusal_name: PROGRAM,
    version: "1.0.0",
    help: HELP,
    descs: DESCS,
    dependencies: no_dependencies,
    post_prefix: None,
};

// NCBI reference (598d8ae6): c++/src/app/dustmasker/main.cpp:39-42
// ```c++
// int main( int argc, char * argv[] )
// {
//     return CDustMaskApplication().AppMain(argc, argv);
// }
// ```
/// Runs `dustmasker` with `words` (the command line after `losat dustmasker`) on the process's
/// standard streams; the exit status.
pub(crate) fn run(words: &[OsString]) -> i32 {
    crate::ncbi_app::run(&TOOL, words, app::run)
}

/// The run of one command line, its output on `out` and its errors on `err`; the exit status
/// (NCBI's application around dustmasker's `Run`: `ncbi_app::run_with`).
#[cfg_attr(
    not(test),
    expect(dead_code, reason = "the library API does not run dustmasker yet")
)]
pub(crate) fn run_with(words: &[OsString], out: &mut dyn Write, err: &mut dyn Write) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, out, err, app::run)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_words(line: &[&str]) -> (i32, String, String) {
        let words: Vec<OsString> = line.iter().map(OsString::from).collect();
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = run_with(&words, &mut out, &mut err);
        (
            status,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    #[test]
    fn the_short_usage_is_the_head_of_the_help() {
        let usage = crate::ncbi_app::short_usage(HELP);
        assert_eq!(usage.len(), 443);
        assert!(
            usage.starts_with("USAGE\n  dustmasker [-h] [-help] [-xmlhelp] [-in input_file_name]")
        );
        assert_eq!(HELP.len(), 1673);
    }

    #[test]
    fn the_version_is_dustmasker_s() {
        assert_eq!(
            run_words(&["-version"]),
            (
                0,
                format!(
                    "dustmasker: 1.0.0\n Package: LOSAT {}, output compatible with blast {}\n",
                    env!("CARGO_PKG_VERSION"),
                    crate::NCBI_BLAST_RELEASE
                ),
                String::new()
            )
        );
    }

    // The NCBI C++ Toolkit options that LOSAT's dustmasker refuses, as blastdbcmd does.
    #[test]
    fn toolkit_options_are_refused_as_blastdbcmd_refuses_them() {
        for (line, word) in [
            (&["-xmlhelp"][..], "-xmlhelp"),
            (&["-version-full"], "-version-full"),
            (&["-version-full-xml"], "-version-full-xml"),
            (&["-version-full-json"], "-version-full-json"),
            (&["-dryrun"], "-dryrun"),
        ] {
            assert_eq!(
                run_words(line),
                (
                    1,
                    String::new(),
                    format!(
                        "Error: the NCBI BLAST+ option {word} is not supported by LOSAT's \
                         DUSTMASKER\n"
                    )
                ),
                "{line:?}"
            );
        }
    }

    // The file of `-conffile` is read by the command line's application start; the library API,
    // which these runs are (no start), refuses it as blastdbcmd's does.
    #[test]
    fn conffile_is_refused_without_the_command_line_s_application_start() {
        assert_eq!(
            run_words(&["-conffile", "x.ini"]),
            (
                1,
                String::new(),
                format!(
                    "Error: {}\n",
                    crate::cli::toolkit_file_word_refusal("dustmasker", "-conffile")
                )
            )
        );
    }

    // `s_StringToInt8` reads a hexadecimal value of a key without a constraint (NCBI 2.17.0:
    // `-window 0x10` runs, probe `big_win_0x10`); `-window 0x` is 0.
    #[test]
    fn integer_keys_take_hexadecimal_values() {
        let words =
            |line: &[&str]| -> Vec<String> { line.iter().map(|word| word.to_string()).collect() };
        let args = crate::ncbi_app::args::parse(&TOOL, &words(&["-window", "0x10"])).unwrap();
        assert_eq!(args.integer("window"), Some(16));
        assert_eq!(args.integer("level"), Some(20));
        assert_eq!(args.integer("linker"), Some(1));
        assert_eq!(args.value("outfmt"), Some("interval"));
        assert_eq!(args.value("infmt"), Some("fasta"));
        assert!(!args.flag("parse_seqids") && !args.given("window-x"));
        let args = crate::ncbi_app::args::parse(&TOOL, &words(&["-window=0X"])).unwrap();
        assert_eq!(args.integer("window"), Some(0));
    }
}
