use super::*;

// A run of dustmasker's arguments (`include_str!` of its help) whose body writes one line into its
// `-out`: the status, the standard output and the standard error.
#[cfg(unix)]
fn run_dustmasker_words(words: Vec<OsString>) -> (i32, Vec<u8>, Vec<u8>) {
    static TOOL: Tool = Tool {
        program: "dustmasker",
        refusal_name: "DUSTMASKER",
        help: include_str!("../dustmasker/help.txt"),
        post_prefix: None,
        components: &[],
        help_full: None,
        short_description: None,
        dry_run: None,
    };
    fn body(
        args: &Args,
        _: &ApplicationSettings,
        _: &mut files::Cin<'_>,
        out: &mut dyn Write,
        _: &mut Diag,
    ) -> i32 {
        let mut output = args.output_file("out", out).unwrap();
        output.write_all(b"ran\n").unwrap();
        0
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let status = run_with(
        &TOOL,
        &words,
        &mut std::io::empty(),
        &mut out,
        &mut err,
        body,
    );
    (status, out, err)
}

// NCBI prints the versions before it reads the arguments, whatever their bytes (review 1 F4 of
// line RG item E3: NCBI BLAST+ 2.17.0 `blastdbcmd -out $'\xff' -version-full`, rc 0), and its
// reader takes a word's bytes as they are: its errors and its help words come first, and quote
// the bytes (line RG item K4i-2a, probe `K/k4i2/nonutf8`). The body then reads each value by
// its bytes (`Args::bytes`), as NCBI's strings hold them (line RG item K10a, probe `K/kd`).
#[cfg(unix)]
#[test]
fn a_word_that_is_not_utf8_is_read_by_its_bytes() {
    use std::os::unix::ffi::OsStringExt;
    let run = run_dustmasker_words;
    let bytes = || OsString::from_vec(vec![0xff]);
    for (word, request) in [
        ("-version", crate::ncbi_version::VersionRequest::Version),
        (
            "-version-full-json",
            crate::ncbi_version::VersionRequest::FullJson,
        ),
    ] {
        let printed = crate::ncbi_version::print(request, "dustmasker");
        assert_eq!(
            run(vec!["-infmt".into(), bytes(), word.into()]),
            (0, printed.into_bytes(), Vec::new())
        );
    }
    // The body runs with such values.
    for words in [
        vec!["-infmt".into(), bytes()],
        vec!["-in".into(), bytes(), "-infmt".into(), bytes()],
    ] {
        assert_eq!(run(words), (0, b"ran\n".to_vec(), Vec::new()));
    }
    // `-dryrun` runs no body: status 0 (NCBI BLAST+ 2.17.0 dustmasker `-infmt $'\xff' -dryrun`).
    assert_eq!(
        run(vec!["-infmt".into(), bytes(), "-dryrun".into()]),
        (0, Vec::new(), Vec::new())
    );
    // A help word stops the reader: the help, status 0.
    let (status, out, err) = run(vec!["-infmt".into(), bytes(), "-h".into()]);
    assert_eq!((status, err), (0, Vec::new()));
    assert!(out.starts_with(b"USAGE\n  dustmasker"), "{out:?}");
    // The reader's errors come first, also after the word: an unknown argument.
    let (status, out, err) = run(vec!["-infmt".into(), bytes(), "-bogus".into()]);
    let err = String::from_utf8(err).unwrap();
    assert_eq!((status, out), (1, Vec::new()));
    assert!(err.contains("Unknown argument: \"bogus\""), "{err}");
    // A refused value is quoted with its bytes.
    let (status, _, err) = run(vec!["-outfmt".into(), bytes()]);
    assert_eq!(status, 1);
    let quoted = b"Illegal value, expected `acclist', ";
    assert!(err.windows(quoted.len()).any(|at| at == quoted), "{err:?}");
    assert!(err.windows(3).any(|at| at == b"`\xff'"), "{err:?}");
    // A name that is not UTF-8 is NCBI's invalid name, with its bytes.
    let (status, _, err) = run(vec![OsString::from_vec(b"-o\xff".to_vec())]);
    assert_eq!(status, 1);
    let quoted = b"Invalid argument name: o\xff";
    assert!(err.windows(quoted.len()).any(|at| at == quoted), "{err:?}");
}

// `-out` is opened by its bytes, as NCBI's `AsOutputFile` opens it. Not on macOS, whose APFS
// refuses a file name that is not UTF-8 (as `tests/non_utf8_names.rs`).
#[cfg(all(unix, not(target_os = "macos")))]
#[test]
fn an_out_name_that_is_not_utf8_is_opened_by_its_bytes() {
    use std::os::unix::ffi::OsStringExt;
    let run = run_dustmasker_words;
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-ncbi-app-out-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let mut path = dir.clone().into_os_string().into_vec();
    path.extend_from_slice(b"/o\xff");
    let path = OsString::from_vec(path);
    assert_eq!(
        run(vec!["-out".into(), path.clone()]),
        (0, Vec::new(), Vec::new())
    );
    assert_eq!(std::fs::read(&path).unwrap(), b"ran\n");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_run_s_arg_exception_is_reported_after_the_usage_with_its_backlog() {
    static TOOL: Tool = Tool {
        program: "p",
        refusal_name: "P",
        help: "USAGE\n  p [-h]\n\nOPTIONAL ARGUMENTS\n -h\n",
        post_prefix: None,
        components: &[],
        help_full: None,
        short_description: None,
        dry_run: None,
    };
    let mut stderr = Vec::new();
    let mut diag = Diag::new(&DiagTarget::Stderr, &mut stderr, &TOOL);
    let error = ArgException::new("eNoFile", "Argument \"in\". X".to_string());
    assert_eq!(diag.run_arg_exception(error), 1);
    assert_eq!(
        String::from_utf8(stderr).unwrap(),
        // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
        crate::wording::losat_text(&format!(
            "USAGE\n  p [-h]\n\nUse '-help' to print detailed descriptions of command line \
                 arguments\n{}\n\nError: Argument \"in\". X\nError: Application's execution \
                 failed\nError:  (CArgException::eNoFile) Argument \"in\". X\n",
            "=".repeat(72)
        ))
    );
    // `eNoValue` prints no usage.
    let mut stderr = Vec::new();
    let mut diag = Diag::new(&DiagTarget::Stderr, &mut stderr, &TOOL);
    diag.run_arg_exception(ArgException::new("eNoValue", "m".to_string()));
    assert_eq!(
        String::from_utf8(stderr).unwrap(),
        crate::wording::losat_text(
            "Error: m\nError: Application's execution failed\nError:  (CArgException::eNoValue) \
                 m\n"
        )
    );
}
