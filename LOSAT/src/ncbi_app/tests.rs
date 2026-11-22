use super::*;

// NCBI prints the versions before it reads the arguments, whatever their bytes (review 1 F4 of
// line RG item E3: NCBI BLAST+ 2.17.0 `blastdbcmd -out $'\xff' -version-full`, rc 0); a word
// that is not UTF-8 is refused otherwise.
#[cfg(unix)]
#[test]
fn the_versions_come_before_the_refusal_of_a_word_that_is_not_utf8() {
    use std::os::unix::ffi::OsStringExt;
    static TOOL: Tool = Tool {
        program: "p",
        refusal_name: "P",
        help: "USAGE\n  p [-h]\n\nOPTIONAL ARGUMENTS\n -h\n",
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
        (
            status,
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    };
    let bytes = OsString::from_vec(vec![0xff]);
    for (word, request) in [
        ("-version", crate::ncbi_version::VersionRequest::Version),
        (
            "-version-full-json",
            crate::ncbi_version::VersionRequest::FullJson,
        ),
    ] {
        assert_eq!(
            run(vec!["-out".into(), bytes.clone(), word.into()]),
            (0, crate::ncbi_version::print(request, "p"), String::new())
        );
    }
    assert_eq!(
        run(vec!["-out".into(), bytes]),
        (
            1,
            String::new(),
            "Error: a word of the command line is not UTF-8; NCBI BLAST+ reads its bytes, \
                 which is not supported by LOSAT's P\n"
                .to_string()
        )
    );
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
        record: None,
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
        record: None,
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
