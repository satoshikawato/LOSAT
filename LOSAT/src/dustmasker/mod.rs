//! `losat dustmasker`: NCBI BLAST+ 2.17.0 dustmasker 1.0.0 (`CDustMaskApplication`), the masker
//! of low-complexity regions of nucleotide sequences by the symmetric DUST algorithm. Its
//! command line is NCBI's (`ncbi_app`, NCBI's application and argument reader), not clap's:
//! `main` hands it the words after `dustmasker`. The work of its `Run` is in `app`.

mod app;
mod fasta;
mod masker;
mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

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
/// Runs `dustmasker` with `words` (the command line after `losat dustmasker`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, app::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around dustmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_dustmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, app::run)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_words(line: &[&str]) -> (i32, String, String) {
        let words: Vec<OsString> = line.iter().map(OsString::from).collect();
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = run_with(&words, &mut std::io::empty(), &mut out, &mut err);
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

    // The library run refuses `-xmlhelp`, which the command line answers with NCBI's XML, as
    // blastdbcmd's (line RG DECISIONS E2x-1).
    #[test]
    fn the_library_refuses_xmlhelp_as_blastdbcmd_does() {
        assert_eq!(
            run_words(&["-xmlhelp"]),
            (
                1,
                String::new(),
                format!("Error: {}\n", crate::cli::xml_help_refusal("DUSTMASKER"))
            )
        );
    }

    // The full versions are NCBI's forms with LOSAT's values, printed before the arguments are
    // read (line RG item E3; NCBI BLAST+ 2.17.0 `dustmasker -version-full`: `dustmasker: 1.0.0`).
    #[test]
    fn the_full_versions_are_dustmasker_s() {
        use crate::ncbi_version::{print, VersionRequest};
        for (word, request) in [
            ("-version-full", VersionRequest::Full),
            ("-version-full-xml", VersionRequest::FullXml),
            ("-version-full-json", VersionRequest::FullJson),
        ] {
            assert_eq!(
                run_words(&["-level", "x", word, "-version"]),
                (0, print(request, "dustmasker", "1.0.0"), String::new()),
                "{word}"
            );
        }
    }

    /// An output whose every write fails (a full disk).
    struct FullDisk;

    impl Write for FullDisk {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("no space left"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("no space left"))
        }
    }

    // A closed pipe ends NCBI's process by SIGPIPE outside Windows (exit status 141, no
    // message), whatever the format: the run ends at once, silently (PD-LOSAT-CLI-NONSEARCH-
    // DIFFERENCES exception 10; review of line RG item K1b, F2).
    #[test]
    #[cfg(not(windows))]
    fn a_closed_pipe_ends_the_run_silently() {
        /// A pipe whose reader has closed.
        struct ClosedPipe;

        impl Write for ClosedPipe {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
            }

            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }

        let input = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/dustmasker/in/basic.fa")
            .into_os_string();
        for outfmt in ["interval", "fasta", "seqloc_asn1_text", "maskinfo_asn1_bin"] {
            let words = [
                OsString::from("-in"),
                input.clone(),
                OsString::from("-outfmt"),
                OsString::from(outfmt),
            ];
            let mut err = Vec::new();
            let status = run_with(&words, &mut std::io::empty(), &mut ClosedPipe, &mut err);
            assert_eq!((status, err), (0, Vec::new()), "{outfmt}");
        }
    }

    // NCBI 2.17.0 with `-out /dev/full` (review 2 of K1, F1): `interval`, `acclist` and `fasta`
    // write nothing more, silently, exit status 0; the serial writers' `CObjectOStream` reports
    // the failure (exit status 3, 0 or 134 by the case), which LOSAT refuses.
    #[test]
    fn a_failed_write_is_silent_for_the_text_writers_and_refused_for_the_serial_ones() {
        let input = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("dustmasker")
            .join("in")
            .join("basic.fa");
        for (outfmt, status) in [
            ("interval", 0),
            ("acclist", 0),
            ("fasta", 0),
            ("seqloc_asn1_text", 1),
            ("seqloc_xml", 1),
            ("maskinfo_asn1_bin", 1),
            ("maskinfo_xml", 1),
        ] {
            let words = [
                OsString::from("-in"),
                input.clone().into_os_string(),
                OsString::from("-outfmt"),
                OsString::from(outfmt),
            ];
            let mut err = Vec::new();
            assert_eq!(
                run_with(&words, &mut std::io::empty(), &mut FullDisk, &mut err),
                status,
                "{outfmt}"
            );
            let err = String::from_utf8(err).unwrap();
            let expected = match status {
                0 => String::new(),
                _ => format!(
                    "Error: a failed write of the -outfmt {outfmt} output (NCBI's serial stream \
                     reports it) is not supported by LOSAT's DUSTMASKER\n"
                ),
            };
            assert_eq!(err, expected, "{outfmt}");
        }
    }

    // The file of `-conffile` is read, and `-dryrun` run, by the command line's application
    // start; the library API, which these runs are (no start), refuses them as blastdbcmd's does.
    #[test]
    fn conffile_and_dryrun_are_refused_without_the_command_line_s_application_start() {
        for (line, word) in [
            (&["-conffile", "x.ini"][..], "-conffile"),
            (&["-dryrun"], "-dryrun"),
        ] {
            assert_eq!(
                run_words(line),
                (
                    1,
                    String::new(),
                    format!(
                        "Error: {}\n",
                        crate::cli::toolkit_file_word_refusal("dustmasker", word)
                    )
                ),
                "{line:?}"
            );
        }
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
