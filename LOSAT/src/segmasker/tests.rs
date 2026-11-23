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
    assert!(usage.starts_with("USAGE\n  segmasker [-h] [-help] [-xmlhelp] [-in input_file_name]"));
    assert_eq!(HELP.len(), 1582);
}

#[test]
fn the_version_is_losat_s() {
    assert_eq!(
        run_words(&["-version"], b""),
        (
            0,
            format!(
                "segmasker: LOSAT {}\n Package: LOSAT {}, output compatible with blast {}\n",
                env!("CARGO_PKG_VERSION"),
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
    let words = |line: &[&str]| -> Vec<OsString> { line.iter().map(OsString::from).collect() };
    let args = crate::ncbi_app::args::parse(&TOOL, &words(&[])).unwrap();
    assert_eq!(args.integer("window"), Some(12));
    assert_eq!(args.double("locut"), Some(2.2));
    assert_eq!(args.double("hicut"), Some(2.5));
    assert_eq!(args.value("infmt"), Some("fasta"));
    assert_eq!(args.value("outfmt"), Some("interval"));
    let args = crate::ncbi_app::args::parse(
        &TOOL,
        &words(&["-window", "0x10", "-locut", "2e", "-hicut=-nan"]),
    )
    .unwrap();
    assert_eq!(args.integer("window"), Some(16));
    assert_eq!(args.double("locut"), Some(2.0));
    let hicut = args.double("hicut").unwrap();
    assert!(hicut.is_nan() && hicut.is_sign_negative());
    assert!(crate::ncbi_app::args::parse(&TOOL, &words(&["-locut", "2,2"])).is_err());
}

// A `-locut` of more than 18 digits before a final `e`, which NCBI's `StringToDoublePosix`
// adds up in x87 extended precision (refused after the review of line RG item K2, L1, until
// line SW, SW8): NCBI 2.17.0's bits (its corelib, `swebok-improve-20261105/sw8/nfix/probe/x87`)
// and its masks of in/basic.fa with `-hicut 3.5`.
#[test]
fn an_extended_precision_value_is_read_as_ncbi_reads_it() {
    let input = b">seq1 low complexity test\n\
            MKVLAAGIVALLLAAQQQQQQQQQQQQQQQQPPPPPPPPGSGSGSGSGSGSAEDKLMNVWY\n\
            RTSHCFIPKLEDGHKLMNRSTVWY\n";
    let words = |line: &[&str]| -> Vec<OsString> { line.iter().map(OsString::from).collect() };
    for (locut, bits, masks) in [
        (
            "1.23456789012345678901e",
            0x3ff3_c0ca_428c_59fb_u64,
            "15 - 50\n",
        ),
        (
            "3.30000000000000000000000000001e",
            0x400a_6666_6666_6666,
            "2 - 14\n15 - 50\n65 - 81\n",
        ),
    ] {
        let args = crate::ncbi_app::args::parse(&TOOL, &words(&["-locut", locut])).unwrap();
        assert_eq!(args.double("locut").unwrap().to_bits(), bits);
        assert_eq!(
            run_words(&["-locut", locut, "-hicut", "3.5"], input),
            (
                0,
                format!(">seq1 low complexity test\n{masks}"),
                String::new()
            )
        );
    }
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
    // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
    let line = |name: &str, value: &str| {
        crate::wording::losat_text(&format!(
            "NCBI C++ Exception:\n    T0 \"{root}/src/corelib/ncbiargs.cpp\", line 618: Error: \
                 (CArgException::eNoFile) Argument \"{name}\". File is not accessible:  \
                 `{value}'\n\n"
        ))
        .into_owned()
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

// The library run refuses `-xmlhelp`, which the command line answers with NCBI's XML, as
// dustmasker's and blastdbcmd's (line RG DECISIONS E2x-1).
#[test]
fn the_library_refuses_xmlhelp_as_dustmasker_does() {
    assert_eq!(
        run_words(&["-xmlhelp"], b""),
        (
            1,
            String::new(),
            format!("Error: {}\n", crate::cli::xml_help_refusal("SEGMASKER"))
        )
    );
}

// The full versions, as dustmasker's (line RG item E3; NCBI BLAST+ 2.17.0
// `segmasker -version-full`: `segmasker: 1.0.0`).
#[test]
fn the_full_versions_are_segmasker_s() {
    use crate::ncbi_version::{print, VersionRequest};
    for (word, request) in [
        ("-version-full", VersionRequest::Full),
        ("-version-full-xml", VersionRequest::FullXml),
        ("-version-full-json", VersionRequest::FullJson),
    ] {
        assert_eq!(
            run_words(&["-out", word], b""),
            (0, print(request, "segmasker"), String::new()),
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

// The library API's run (`run_with`) runs NCBI's application start as the command line's
// does (line RG item E1w-2): the file of `-conffile` is loaded in place of `segmasker.ini` (a
// name that is not a file is NCBI's critical post and abort, exit status 255: the table's row
// `conf_missing`), and `-dryrun` ends the run after the arguments, reading no input.
#[test]
fn conffile_and_dryrun_run_as_on_the_command_line() {
    let missing = "losat-segmasker-unit-missing.ini";
    assert_eq!(
        run_words(&["-conffile", missing], b""),
        (
            255,
            String::new(),
            format!("Critical: (106.23) Registry file \"{missing}\" not found, aborting.\n")
        )
    );
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-segmasker-unit-conffile-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let conffile = dir.join("c.ini");
    std::fs::write(&conffile, "[BLAST]\n").unwrap();
    let conffile = conffile.to_str().unwrap();
    let missing_input = dir.join("missing.fa");
    let missing_input = missing_input.to_str().unwrap();
    let record = b">p\nMKKLLPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPGGGGGHHHH\n";
    for line in [
        &["-dryrun"][..],
        &["-in", missing_input, "-dryrun"],
        &["-conffile", conffile, "-dryrun"],
    ] {
        assert_eq!(
            run_words(line, record),
            (0, String::new(), String::new()),
            "{line:?}"
        );
    }
    // A file that loads changes nothing.
    let plain = run_words(&[], record);
    assert!(!plain.1.is_empty());
    assert_eq!(run_words(&["-conffile", conffile], record), plain);
    let _ = std::fs::remove_dir_all(&dir);
}
