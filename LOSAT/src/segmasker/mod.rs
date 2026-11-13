//! `losat segmasker`: NCBI BLAST+ 2.17.0 segmasker 1.0.0 (`SegMaskerApplication`,
//! c++/src/app/segmasker/segmasker.cpp), the masker of low-complexity regions of protein
//! sequences by the SEG algorithm. Its command line is NCBI's (`ncbi_app`), as dustmasker's;
//! its input and output are the masking programs' (`seqmasks_io`), and its SEG is BLAST's
//! (`utils::seg::SegMasker`, without merging overlapping segments). The work of its `Run` is in
//! `work`.

mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, Constraint, Desc, Kind};
use crate::ncbi_app::{no_dependencies, Tool};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "SEGMASKER";

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:96-104
// ```c++
// const char * const SegMaskerApplication::USAGE_LINE
//     = "Low complexity region masker based on the SEG algorithm";
//
// void SegMaskerApplication::Init(void)
// {
//     HideStdArgs(fHideLogfile | fHideConffile | fHideVersion | fHideDryRun);
// ```
/// The text of `segmasker -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`, `-version`
/// and `-dryrun` are hidden. Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the
/// argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:112-150
// ```c++
//     arg_desc->SetCurrentGroup("Input/output options");
//     arg_desc->AddDefaultKey(kInput, "input_file_name",
//                             "input file name",
//                             CArgDescriptions::eInputFile, "-");
//     arg_desc->AddDefaultKey(kOutput, "output_file_name",
//                             "output file name",
//                             CArgDescriptions::eOutputFile, "-");
//     arg_desc->AddDefaultKey(kInputFormat, "input_format",
//                             "controls the format of the masker input",
//                             CArgDescriptions::eString, *kInputFormats);
//     CArgAllow_Strings* strings_allowed = new CArgAllow_Strings();
//     for (size_t i = 0; i < kNumInputFormats; i++) {
//         strings_allowed->Allow(kInputFormats[i]);
//     }
//     arg_desc->SetConstraint(kInputFormat, strings_allowed);
//     arg_desc->AddFlag      ( "parse_seqids",
//                              "Parse Seq-ids in FASTA input", true );
//
//     arg_desc->AddDefaultKey(kOutputFormat, "output_format",
//                             "controls the format of the masker output",
//                             CArgDescriptions::eString, *kOutputFormats);
//     strings_allowed = new CArgAllow_Strings();
//     for (size_t i = 0; i < kNumOutputFormats; i++) {
//         strings_allowed->Allow(kOutputFormats[i]);
//     }
//     arg_desc->SetConstraint(kOutputFormat, strings_allowed);
//
//     arg_desc->SetCurrentGroup("SEG algorithm options");
//     arg_desc->AddDefaultKey("window", "integer_value", "SEG window",
//                             CArgDescriptions::eInteger,
//                             NStr::IntToString(kSegWindow));
//     arg_desc->AddDefaultKey("locut", "float_value", "SEG locut",
//                             CArgDescriptions::eDouble,
//                             NStr::DoubleToString(kSegLocut));
//     arg_desc->AddDefaultKey("hicut", "float_value", "SEG hicut",
//                             CArgDescriptions::eDouble,
//                             NStr::DoubleToString(kSegHicut));
// ```
// (`kInputFormats` and the eight `kOutputFormats`: mask_cmdline_args.cpp:40-51, quoted at
// dustmasker's descriptions; segmasker allows no `acclist`. The defaults `12`, `2.2`, `2.5` are
// those of segmasker 2.17.0's `-help`.)
/// The arguments of `SegMaskerApplication::Init` in the order it describes them; the values that
/// `-infmt` and `-outfmt` allow are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key(
        "infmt",
        Kind::Text,
        Some("fasta"),
        Some(Constraint::Strings(&["blastdb", "fasta"])),
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
    key("window", Kind::Integer, Some("12"), None),
    key("locut", Kind::Double, Some("2.2"), None),
    key("hicut", Kind::Double, Some("2.5"), None),
];

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:69-74
// ```c++
//     SegMaskerApplication() {
//         CRef<CVersion> version(new CVersion());
//         version->SetVersionInfo(1, 0, 0);
//         SetFullVersion(version);
//     }
// ```
/// segmasker as NCBI's application runs it (`ncbi_app`): its own version, no dependencies
/// between its arguments and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "segmasker",
    refusal_name: PROGRAM,
    version: "1.0.0",
    help: HELP,
    descs: DESCS,
    dependencies: no_dependencies,
    post_prefix: None,
};

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:265-272
// ```c++
// int main(int argc, const char* argv[])
// {
//     // Execute main application function
//     return SegMaskerApplication().AppMain(argc, argv);
// }
// ```
/// Runs `segmasker` with `words` (the command line after `losat segmasker`) on the process's
/// standard streams; the exit status.
pub(crate) fn run(words: &[OsString]) -> i32 {
    crate::ncbi_app::run(&TOOL, words, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around segmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_segmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, work::run)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_words(line: &[&str], input: &[u8]) -> (i32, String, String) {
        let words: Vec<OsString> = line.iter().map(OsString::from).collect();
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = run_with(&words, &mut &input[..], &mut out, &mut err);
        (
            status,
            String::from_utf8_lossy(&out).into_owned(),
            String::from_utf8(err).unwrap(),
        )
    }

    #[test]
    fn the_short_usage_is_the_head_of_the_help() {
        let usage = crate::ncbi_app::short_usage(HELP);
        assert_eq!(usage.len(), 431);
        assert!(
            usage.starts_with("USAGE\n  segmasker [-h] [-help] [-xmlhelp] [-in input_file_name]")
        );
        assert_eq!(HELP.len(), 1582);
    }

    #[test]
    fn the_version_is_segmasker_s() {
        assert_eq!(
            run_words(&["-version"], b""),
            (
                0,
                format!(
                    "segmasker: 1.0.0\n Package: LOSAT {}, output compatible with blast {}\n",
                    env!("CARGO_PKG_VERSION"),
                    crate::NCBI_BLAST_RELEASE
                ),
                String::new()
            )
        );
    }

    // NCBI 2.17.0 segmasker on in/basic.fa of the fixtures, from the standard input.
    #[test]
    fn the_standard_input_is_masked() {
        let input = b">seq1 low complexity test\n\
            MKVLAAGIVALLLAAQQQQQQQQQQQQQQQQPPPPPPPPGSGSGSGSGSGSAEDKLMNVWY\n\
            RTSHCFIPKLEDGHKLMNRSTVWY\n";
        assert_eq!(
            run_words(&[], input),
            (
                0,
                ">seq1 low complexity test\n2 - 14\n15 - 50\n".to_string(),
                String::new()
            )
        );
    }

    // The keys as `CArgs` gives them: `-window 0x10` is 16 (`s_StringToInt8`), `-locut 2e` is 2
    // (`StringToDoublePosix` reads the mark), a comma is no decimal point.
    #[test]
    fn the_keys_read_as_ncbis() {
        let words =
            |line: &[&str]| -> Vec<String> { line.iter().map(|word| word.to_string()).collect() };
        let args = crate::ncbi_app::args::parse(&TOOL, &words(&[])).unwrap();
        assert_eq!(args.integer("window"), Some(12));
        assert_eq!(args.double("locut"), Some(Ok(2.2)));
        assert_eq!(args.double("hicut"), Some(Ok(2.5)));
        assert_eq!(args.value("infmt"), Some("fasta"));
        assert_eq!(args.value("outfmt"), Some("interval"));
        let args = crate::ncbi_app::args::parse(
            &TOOL,
            &words(&["-window", "0x10", "-locut", "2e", "-hicut=-nan"]),
        )
        .unwrap();
        assert_eq!(args.integer("window"), Some(16));
        assert_eq!(args.double("locut"), Some(Ok(2.0)));
        let hicut = args.double("hicut").unwrap().unwrap();
        assert!(hicut.is_nan() && hicut.is_sign_negative());
        assert!(crate::ncbi_app::args::parse(&TOOL, &words(&["-locut", "2,2"])).is_err());
    }

    // A `-locut` that NCBI computes in extended precision is read as NCBI's reader reads it (the
    // usage of `-h`, rc 0), and refused where `Run` reads it (review of line RG item K2, L1).
    #[test]
    fn an_extended_precision_value_is_refused_by_the_run() {
        let value = "1.23456789012345678901e";
        let (status, out, err) = run_words(&["-locut", value, "-h"], b"");
        assert_eq!(
            (status, out, err),
            (0, crate::ncbi_app::short_usage(HELP), String::new())
        );
        assert_eq!(
            run_words(&["-locut", value], b">a\nQ\n"),
            (
                1,
                String::new(),
                format!(
                    "Error: the value '{value}' of -locut, which NCBI computes in extended \
                     precision, is not supported by LOSAT's SEGMASKER\n"
                )
            )
        );
    }

    // A database byte above 127 (tests/segmasker/db/r201: the byte 201 for two X) is undefined in
    // NCBI's SEG, and 201 stops its window: refused (review of line RG item K2, M1).
    #[test]
    fn a_database_byte_above_127_is_refused() {
        let db = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/segmasker/db/r201");
        let db = db.to_str().unwrap();
        assert_eq!(
            run_words(&["-infmt", "blastdb", "-in", db], b""),
            (
                1,
                String::new(),
                "Error: a database residue byte 201, above 127 (NCBI's SEG reads it as a negative \
                 char, and stops its window at 201), is not supported by LOSAT's SEGMASKER\n"
                    .to_string()
            )
        );
    }

    // `-in` and `-out` that do not open: `Run` catches the `CArgException` and writes
    // `e.what()`, exit status 1; the reader is made before the writer, so `-in` is reported first
    // (NCBI 2.17.0, fixtures `out_bad_noin`, `out_bad`).
    #[test]
    fn files_that_do_not_open_are_reported_by_what() {
        let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
        let line = |name: &str, value: &str| {
            format!(
                "NCBI C++ Exception:\n    T0 \"{root}/src/corelib/ncbiargs.cpp\", line 618: Error: \
                 (CArgException::eNoFile) Argument \"{name}\". File is not accessible:  \
                 `{value}'\n\n"
            )
        };
        assert_eq!(
            run_words(&["-in", "/nonexistent/x.fa", "-out", "/nonexistent/y"], b""),
            (1, String::new(), line("in", "/nonexistent/x.fa"))
        );
        assert_eq!(
            run_words(&["-out", "/nonexistent/y"], b">a\nQQQ\n"),
            (1, String::new(), line("out", "/nonexistent/y"))
        );
    }

    // NCBI's toolkit options that LOSAT's segmasker refuses, as dustmasker and blastdbcmd do.
    #[test]
    fn toolkit_options_are_refused_as_dustmasker_refuses_them() {
        for word in [
            "-xmlhelp",
            "-version-full",
            "-version-full-xml",
            "-version-full-json",
        ] {
            assert_eq!(
                run_words(&[word], b""),
                (
                    1,
                    String::new(),
                    format!(
                        "Error: the NCBI BLAST+ option {word} is not supported by LOSAT's \
                         SEGMASKER\n"
                    )
                ),
                "{word}"
            );
        }
    }

    /// A pipe whose reader has closed.
    #[cfg(not(windows))]
    struct ClosedPipe;

    #[cfg(not(windows))]
    impl Write for ClosedPipe {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // A closed pipe ends NCBI's process by SIGPIPE outside Windows: the run ends at once,
    // silently, exit status 0, as dustmasker's (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10).
    #[test]
    #[cfg(not(windows))]
    fn a_closed_pipe_ends_the_run_silently() {
        for outfmt in ["interval", "fasta", "seqloc_asn1_text", "maskinfo_asn1_bin"] {
            let words = [OsString::from("-outfmt"), OsString::from(outfmt)];
            let mut err = Vec::new();
            let status = run_with(
                &words,
                &mut &b">a\nQQQQQQQQQQQQQQQQQQQQ\n"[..],
                &mut ClosedPipe,
                &mut err,
            );
            assert_eq!((status, err), (0, Vec::new()), "{outfmt}");
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

    // As dustmasker's (review 2 of K1, F1): the text writers go on silently, exit status 0; the
    // serial writers' failure, which NCBI's `CObjectOStream` reports, is LOSAT's refusal.
    #[test]
    fn a_failed_write_is_silent_for_the_text_writers_and_refused_for_the_serial_ones() {
        for (outfmt, status) in [
            ("interval", 0),
            ("fasta", 0),
            ("seqloc_asn1_text", 1),
            ("maskinfo_xml", 1),
        ] {
            let words = [OsString::from("-outfmt"), OsString::from(outfmt)];
            let mut err = Vec::new();
            assert_eq!(
                run_with(
                    &words,
                    &mut &b">a\nQQQQQQQQQQQQQQQQQQQQ\n"[..],
                    &mut FullDisk,
                    &mut err
                ),
                status,
                "{outfmt}"
            );
            let expected = match status {
                0 => String::new(),
                _ => format!(
                    "Error: a failed write of the -outfmt {outfmt} output (NCBI's serial stream \
                     reports it) is not supported by LOSAT's SEGMASKER\n"
                ),
            };
            assert_eq!(String::from_utf8(err).unwrap(), expected, "{outfmt}");
        }
    }

    // The file of `-conffile` is read, and `-dryrun` run, by the command line's application
    // start; the library API, which these runs are (no start), refuses them as dustmasker's does.
    #[test]
    fn conffile_and_dryrun_are_refused_without_the_command_line_s_application_start() {
        for (line, word) in [
            (&["-conffile", "x.ini"][..], "-conffile"),
            (&["-dryrun"], "-dryrun"),
        ] {
            assert_eq!(
                run_words(line, b""),
                (
                    1,
                    String::new(),
                    format!(
                        "Error: {}\n",
                        crate::cli::toolkit_file_word_refusal("segmasker", word)
                    )
                ),
                "{line:?}"
            );
        }
    }
}
