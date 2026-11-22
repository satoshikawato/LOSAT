use super::*;

// NCBI prints the versions before it reads the arguments, whatever their bytes (review 1 F4 of
// line RG item E3: NCBI BLAST+ 2.17.0 `blastdbcmd -out $'\xff' -version-full`, rc 0), and its
// reader takes a word's bytes as they are: its errors and its help words come first, and quote
// the bytes (line RG item K4i-2a, probe `K/k4i2/nonutf8`). LOSAT refuses a value that is not
// UTF-8 after the read, but that of `-logfile`, which the application opens by its bytes.
#[cfg(unix)]
#[test]
fn a_word_that_is_not_utf8_is_read_first_and_refused_after_the_read() {
    use std::os::unix::ffi::OsStringExt;
    static TOOL: Tool = Tool {
        program: "dustmasker",
        refusal_name: "DUSTMASKER",
        help: include_str!("../dustmasker/help.txt"),
        post_prefix: None,
        components: &[],
    };
    fn body(
        _: &Args,
        _: &ApplicationSettings,
        _: &mut files::Cin<'_>,
        _: &mut dyn Write,
        _: &mut Diag,
    ) -> i32 {
        unreachable!("the run")
    }
    let run = |words: Vec<OsString>| {
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
    };
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
            run(vec!["-out".into(), bytes(), word.into()]),
            (0, printed.into_bytes(), Vec::new())
        );
    }
    let refusal = b"Error: a word of the command line is not UTF-8; NCBI BLAST+ reads its bytes, \
                    which is not supported by LOSAT's DUSTMASKER\n";
    assert_eq!(
        run(vec!["-out".into(), bytes()]),
        (1, Vec::new(), refusal.to_vec())
    );
    // A help word stops the reader: the help, status 0.
    let (status, out, err) = run(vec!["-out".into(), bytes(), "-h".into()]);
    assert_eq!((status, err), (0, Vec::new()));
    assert!(out.starts_with(b"USAGE\n  dustmasker"), "{out:?}");
    // The reader's errors come first, also after the word: an unknown argument.
    let (status, out, err) = run(vec!["-out".into(), bytes(), "-bogus".into()]);
    let err = String::from_utf8(err).unwrap();
    assert_eq!((status, out), (1, Vec::new()));
    assert!(err.contains("Unknown argument: \"bogus\""), "{err}");
    assert!(!err.contains("not UTF-8"), "{err}");
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

#[test]
fn a_run_s_arg_exception_is_reported_after_the_usage_with_its_backlog() {
    static TOOL: Tool = Tool {
        program: "p",
        refusal_name: "P",
        help: "USAGE\n  p [-h]\n\nOPTIONAL ARGUMENTS\n -h\n",
        post_prefix: None,
        components: &[],
    };
    let mut stderr = Vec::new();
    let mut diag = Diag {
        target: DiagTarget::Stderr,
        stderr: &mut stderr,
        tool: &TOOL,
    };
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
    let mut diag = Diag {
        target: DiagTarget::Stderr,
        stderr: &mut stderr,
        tool: &TOOL,
    };
    diag.run_arg_exception(ArgException::new("eNoValue", "m".to_string()));
    assert_eq!(
        String::from_utf8(stderr).unwrap(),
        crate::wording::losat_text(
            "Error: m\nError: Application's execution failed\nError:  (CArgException::eNoValue) \
                 m\n"
        )
    );
}
