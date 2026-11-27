use super::*;

/// A chain of 1,000,000 rethrows (the frames that an error at the bottom of a PSSM nested
/// 500,000 deep leaves, `NCBI_RETHROW_SAME` in each) is dropped, and one of 100,000 is
/// cloned, compared, shown and reported, on a thread of 256 KiB (line RG D3-1 review 7 B1).
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_chain_of_any_length_is_handled_link_by_link() {
    let chain = |links: usize| {
        let mut chain = NcbiException::new(ThrowSite::generated("CSerialException"), "");
        for _ in 0..links {
            chain = chain.rethrow_same("src/serial/objistrasn.cpp", 1, "f", String::new());
        }
        chain
    };
    /// The bytes that a format writes.
    struct Count(usize);
    impl std::fmt::Write for Count {
        fn write_str(&mut self, s: &str) -> std::fmt::Result {
            self.0 += s.len();
            Ok(())
        }
    }
    std::thread::Builder::new()
        .stack_size(256 << 10)
        .spawn(move || {
            drop(chain(1_000_000));
            let long = chain(100_000);
            let copy = long.clone();
            assert!(copy == long);
            let other = copy.clone().rethrow_same("x", 2, "g", "one more");
            assert!(other != long);
            assert_eq!(
                long.what().iter().filter(|&&byte| byte == b'\n').count(),
                1 + 100_001
            );
            let mut shown = Count(0);
            std::fmt::write(&mut shown, format_args!("{long:?}")).expect("shown");
            assert!(shown.0 > 100_000);
        })
        .expect("a thread")
        .join()
        .expect("no overflow");
}

// Each run's stream is its own (line RG D3-1 review 4 B1): a library call entered on a thread
// reads its stream there, another call entered in turn on the same thread reads its own, and
// leaving a call gives the thread back the stream it had; a thread that runs no call reads
// the command line's.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn each_run_has_a_diagnostic_stream_of_its_own() {
    let a = std::sync::Arc::new(MainDiag::new());
    let b = std::sync::Arc::new(MainDiag::new());
    {
        let _a = MainDiag::enter(Some(a.clone()));
        assert!(!main_diag_failbit());
        set_main_diag_failbit();
        assert!(main_diag_failbit());
        {
            let _b = MainDiag::enter(Some(b.clone()));
            assert!(!main_diag_failbit());
        }
        assert!(main_diag_failbit());
    }
    assert!(a.failed());
    assert!(!b.failed());
    assert!(!COMMAND_LINE.failed());
    let entered = || CURRENT.with(|current| current.borrow().is_some());
    assert!(!entered());
    let _b = MainDiag::enter(Some(b.clone()));
    assert!(!main_diag_failbit());
    // Another thread runs no call until it enters one.
    std::thread::scope(|scope| {
        scope.spawn(|| assert!(!entered()));
    });
}

// makeblastdb 2.17.0 after the object manager failed the stream on its ASN.1 input (line RG D3c-6a,
// rows d3c6a_y08, y09, P_x49 of the makeblastdb ASN.1 manifest): each later post is its header
// and its line feed.
#[test]
fn a_post_after_the_stream_failed_is_its_header() {
    let diag = std::sync::Arc::new(MainDiag::new());
    let _run = MainDiag::enter(Some(diag));
    let header = b"Warning: [makeblastdb] ";
    assert_eq!(
        post_text(header, b"No sequences written"),
        b"Warning: [makeblastdb] No sequences written"
    );
    set_main_diag_failbit();
    assert_eq!(post_text(header, b"No sequences written"), header);
}

// makeblastdb 2.17.0 (the official Linux x64 release) with `MAX_LMDB_TXN_ENTRY=-1x` and two
// sequences (probe of 2026-11-08, case txn_endptr_2seq): a report with the position.
#[test]
fn a_string_exception_reports_its_position() {
    let site = ThrowSite {
        file: "src/corelib/ncbistr.cpp",
        line: 870,
        kind: "CStringException::eConvert",
        function: Some("ncbi::NStr::StringToInt8"),
    };
    let e = NcbiException::with_pos(site, "Cannot convert string '-1x' to Int8", 2);
    assert_eq!(e.get_msg(), b"Cannot convert string '-1x' to Int8");
    assert_eq!(e.site(), site);
    assert_eq!(
        String::from_utf8(e.what()).unwrap(),
        crate::wording::losat_text(
            "NCBI C++ Exception:\n    T0 \"/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../../src/corelib/ncbistr.cpp\", line 870: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string '-1x' to Int8 (m_Pos = 2)\n",
        )
    );
}

// blastn 2.17.0 (the official Linux x64 release) `-window_masker_db` with a counts file whose
// header is cut short (line NP, M/STATE.md "M2 の実機確認"): the reader's exception, without a
// function, and the factory's rethrow, oldest first.
#[test]
fn a_rethrow_reports_the_chain_oldest_first() {
    let reader = NcbiException::new(
        ThrowSite {
            file: "src/algo/winmask/seq_masker_istat_obinary.cpp",
            line: 63,
            kind: "Exception::format error",
            function: None,
        },
        "file too short",
    );
    let factory = reader.rethrow(
        ThrowSite {
            file: "src/algo/winmask/seq_masker_istat_factory.cpp",
            line: 358,
            kind: "Exception::creation failure",
            function: Some("ncbi::CSeqMaskerIstatFactory::create"),
        },
        "could not create a unit counts container",
    );
    assert_eq!(
        factory.get_msg(),
        b"could not create a unit counts container"
    );
    assert_eq!(
        String::from_utf8(factory.what()).unwrap(),
        crate::wording::losat_text(
            "NCBI C++ Exception:\n    T0 \"/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../../src/algo/winmask/seq_masker_istat_obinary.cpp\", line 63: Error: (Exception::format error) file too short\n    T0 \"/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../../src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n",
        )
    );
}

// blastn 2.17.0 (the official Linux x64 release) with a `blastn.ini` of the line `garbage`
// (env parity g_reg probe `x_ini`, 2026-11-11): `NCBI_REPORT_EXCEPTION_X` (subcode 15, title
// "Application's initialization failed") of the application layer (error code 106), and of a chain.
#[test]
fn a_report_with_a_title_has_a_line_per_exception_then_the_title() {
    let site = ThrowSite {
        file: "src/corelib/ncbireg.cpp",
        line: 722,
        kind: "CRegistryException::eEntry",
        function: Some("ncbi::IRWRegistry::x_Read"),
    };
    let e = NcbiException::with_pos(site, "Invalid registry entry format: 'garbage'", 1);
    assert_eq!(
        String::from_utf8(e.report_default_ex(106, 15, "Application's initialization failed")).unwrap(),
        crate::wording::losat_text(
            "Error: (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\n",
        )
    );
    let outer = NcbiException::new(site, "first").rethrow(
        ThrowSite {
            kind: "CAppException::eSetupDiag",
            ..site
        },
        "second",
    );
    assert_eq!(
        String::from_utf8(outer.report_default_ex(106, 15, "T")).unwrap(),
        crate::wording::losat_text(
            "Error: (CRegistryException::eEntry) first\nError: (CAppException::eSetupDiag) second\nError: (106.15) T (CAppException::eSetupDiag) second\n",
        )
    );
}

// NCBI's rpsbproc 0.5.1 (tests/fixtures/rpsbproc/oracle/d_wsline.err and lcl.err): the
// application's report of an exception with the default diagnostic flags, without the place
// or the function, the position of a `CStringException` included; and that of an exception
// of generated code.
#[test]
fn the_default_report_names_the_class_and_the_message() {
    let site = ThrowSite {
        file: "src/corelib/ncbistr.cpp",
        line: 895,
        kind: "CStringException::eConvert",
        function: Some("ncbi::NStr::StringToUInt8"),
    };
    let e = NcbiException::with_pos(site, "Cannot convert string '' to Uint8", 0);
    assert_eq!(
        String::from_utf8(e.report_default_ex(106, 15, "Application's initialization failed")).unwrap(),
        crate::wording::losat_text(
            "Error: (CStringException::eConvert) Cannot convert string '' to Uint8 (m_Pos = 0)\nError: (106.15) Application's initialization failed (CStringException::eConvert) Cannot convert string '' to Uint8 (m_Pos = 0)\n",
        )
    );
    let e = NcbiException::new(
        ThrowSite::generated("CInvalidChoiceSelection::eFail"),
        "CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local",
    );
    assert_eq!(
        String::from_utf8(e.report_default_ex(106, 16, "Application's execution failed")).unwrap(),
        crate::wording::losat_text(
            "Error: (CInvalidChoiceSelection::eFail) CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local\nError: (106.16) Application's execution failed (CInvalidChoiceSelection::eFail) CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local\n",
        )
    );
}

/// The source files of the crate below `dir`, with their paths.
// The text of an exception report has one owner: a second format string for it (as four
// modules had until 2026-11-08: winmask, the Seq-id parser, the taxid map and the seqid list
// filter) fails here. Test files (`tests.rs`, `*_tests.rs`, files under a `tests/` directory)
// and test modules hold NCBI's bytes as expected values and are not read. The needles are
// the escaped forms in ordinary string literals: a raw string (`r#"..."#`) or a text put
// together with `concat!` or from pieces is not caught.
#[test]
fn exception_reports_are_formatted_only_here() {
    use crate::wording::tests::{code_lines, rust_files};
    let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut files = Vec::new();
    rust_files(&src, &mut files);
    assert!(
        files.len() > 100,
        "{} files below {}",
        files.len(),
        src.display()
    );
    let mut offenders = Vec::new();
    for path in files {
        let relative = path.strip_prefix(&src).unwrap();
        let stem = path.file_stem().unwrap().to_string_lossy();
        let test_file = stem == "tests"
            || stem.ends_with("_tests")
            || relative
                .components()
                .any(|part| part.as_os_str() == "tests");
        if relative == std::path::Path::new("ncbi_exception.rs") || test_file {
            continue;
        }
        let text = std::fs::read_to_string(&path).unwrap();
        for (number, line) in code_lines(&text) {
            if ["NCBI C++ Exception:\\n", "T0 \\\"", ": Error: ("]
                .iter()
                .any(|needle| line.contains(needle))
            {
                offenders.push(format!("{}:{number}: {}", relative.display(), line.trim()));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "format an NCBI exception report with ncbi_exception::NcbiException:\n{}",
        offenders.join("\n")
    );
}

/// `add_frame_infos` of a list of names is `CSerialException::AddFrameInfo` of each in turn
/// (serial/exception.cpp: a name that is not empty goes in front of the stack, after a `.` when the
/// stack has one; an empty one changes nothing), also on a stack that a rethrow kept.
#[test]
fn frame_names_are_added_as_ncbi_adds_them_one_by_one() {
    fn one_by_one(stack: &mut String, name: &str) {
        if !name.is_empty() && !stack.is_empty() {
            stack.insert(0, '.');
        }
        stack.insert_str(0, name);
    }
    let site = ThrowSite::generated("CSerialException");
    let names = [
        "E",
        "",
        "strand",
        "",
        "[]",
        "[]",
        "mix",
        "",
        "masks",
        "Blast-db-mask-info",
    ];
    for kept in [None, Some("x.y")] {
        for split in 0..=names.len() {
            let mut e = NcbiException::new(site, "m");
            let mut stack = String::new();
            if let Some(kept) = kept {
                e.add_frame_infos(&[kept]);
                one_by_one(&mut stack, kept);
            }
            e.add_frame_infos(&names[..split]);
            e.add_frame_infos(&names[split..]);
            for name in names {
                one_by_one(&mut stack, name);
            }
            let expected = (!stack.is_empty()).then(|| format!(" at {stack}").into_boxed_str());
            assert_eq!(e.extra, expected, "{kept:?} {split}");
        }
    }
    let mut e = NcbiException::new(site, "m");
    e.add_frame_infos(&["", ""]);
    assert_eq!(e.extra, None);
}
