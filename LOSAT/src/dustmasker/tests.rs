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
    assert!(usage.starts_with("USAGE\n  dustmasker [-h] [-help] [-xmlhelp] [-in input_file_name]"));
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

// The library API's run (`run_with`) runs NCBI's application start as the command line's
// does (line RG item E1w-2): the file of `-conffile` is loaded in place of `dustmasker.ini` (a
// name that is not a file is NCBI's critical post and abort, exit status 255: the table's row
// `conf_missing`), and `-dryrun` ends the run after the arguments, reading no input.
#[test]
fn conffile_and_dryrun_run_as_on_the_command_line() {
    let missing = "losat-dustmasker-unit-missing.ini";
    assert_eq!(
        run_words(&["-conffile", missing]),
        (
            255,
            String::new(),
            format!("Critical: (106.23) Registry file \"{missing}\" not found, aborting.\n")
        )
    );
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-dustmasker-unit-conffile-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let conffile = dir.join("c.ini");
    std::fs::write(&conffile, "[BLAST]\n").unwrap();
    let conffile = conffile.to_str().unwrap();
    let missing_input = dir.join("missing.fa");
    let missing_input = missing_input.to_str().unwrap();
    for line in [
        &["-dryrun"][..],
        &["-in", missing_input, "-dryrun"],
        &["-conffile", conffile, "-dryrun"],
        &["-dryrun", "-outfmt", "fasta", "-in", missing_input],
    ] {
        assert_eq!(
            run_words(line),
            (0, String::new(), String::new()),
            "{line:?}"
        );
    }
    // A file that loads changes nothing: the empty input's run, as without it.
    assert_eq!(run_words(&["-conffile", conffile]), run_words(&[]));
    let _ = std::fs::remove_dir_all(&dir);
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
