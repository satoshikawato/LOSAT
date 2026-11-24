use super::*;

#[test]
fn diag_table_windows_over_2_to_the_30_are_rejected() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-57: the doubling loop ends only when qlen + window <= 2^30.
    assert!(check_diag_table_window(1000, (1 << 30) - 1000, "BLASTP").is_ok());
    for window in [(1 << 30) - 999, i32::MAX] {
        let error = format!(
            "{:#}",
            check_diag_table_window(1000, window, "TBLASTX").unwrap_err()
        );
        assert!(
            error.contains("not supported by LOSAT's TBLASTX"),
            "{window}: {error}"
        );
    }
}

#[test]
fn seg_values_that_strtod_does_not_read_are_ncbi_errors() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:405-406,423-428 (StringToDouble with default flags).
    for value in [
        "12 2.5x 2.5",
        "12 2,2 2.5",
        "12 1e 2.5",
        "12 0x 2.5",
        "12 . 2.5",
        "12 nan(1) 2.5",
    ] {
        let error = format!("{:#}", parse_seg_option(value).unwrap_err());
        assert!(
            error.contains("Invalid input for filtering parameters"),
            "{value}: {error}"
        );
    }
    // NCBI 2.17.0 searches with an infinite or NaN locut or hicut (sw8/nfix/probe/seg2), as
    // LOSAT does: the values as read (the filter string prints them), defaults for those not
    // above 0 in the query filter (`SegSpec::params`).
    for (value, locut, hicut) in [
        ("12 +inf 2.5", f64::INFINITY, 2.5),
        ("12 1e400 2.5", f64::INFINITY, 2.5),
        ("12 2.2 -nan", 2.2, -f64::NAN),
        ("12 +nan(1) 2.5", f64::from_bits(0x7ff8_0000_0000_0001), 2.5),
        ("12 -inf -inf", f64::NEG_INFINITY, f64::NEG_INFINITY),
    ] {
        match parse_seg_option(value).unwrap() {
            SegSpec::WindowLocutHicut {
                window: 12,
                locut: read_locut,
                hicut: read_hicut,
            } => assert_eq!(
                (read_locut.to_bits(), read_hicut.to_bits()),
                (locut.to_bits(), hicut.to_bits()),
                "{value}"
            ),
            other => panic!("{value}: {other:?}"),
        }
    }
    let params = |value: &str| parse_seg_option(value).unwrap().params().unwrap();
    let (inf_low, nan_high) = (params("12 +inf 2.5"), params("12 2.2 -nan"));
    assert_eq!(
        (inf_low.locut, inf_low.hicut),
        (f64::INFINITY, f64::INFINITY)
    );
    assert_eq!((nan_high.locut, nan_high.hicut), (2.2, 2.5));
    // NCBI's error of the hicut, with an infinite locut.
    let error = format!("{:#}", parse_seg_option("12 +inf 2.5x").unwrap_err());
    assert!(
        error.contains("Invalid input for filtering parameters"),
        "{error}"
    );
    assert!(parse_seg_option("12 2.2 2.5").is_ok());
    // `strtod` reads hexadecimal numbers: NCBI 2.17.0's `blastp -seg "12 0x1.2p1 0x1.4p1"`
    // searches as `-seg "12 2.25 2.5"`, `-seg "12 0x10 2.5"` with the locut 16.
    for (value, locut, hicut) in [
        ("12 0x1.2p1 0x1.4p1", 2.25, 2.5),
        ("12 0x10 2.5", 16.0, 2.5),
    ] {
        match parse_seg_option(value).unwrap() {
            SegSpec::WindowLocutHicut {
                window: 12,
                locut: read_locut,
                hicut: read_hicut,
            } => assert_eq!((read_locut, read_hicut), (locut, hicut), "{value}"),
            other => panic!("{value}: {other:?}"),
        }
    }
}

fn message(error: anyhow::Error) -> String {
    match error.downcast::<NativeError>() {
        Ok(native) => format!(
            "{} {}",
            native.exit,
            String::from_utf8(native.message).unwrap()
        ),
        Err(other) => format!("losat {other}"),
    }
}

// NCBI BLAST+ 2.17.0 -outfmt spellings (c++/src/algo/blast/blastinput/blast_args.cpp:2801-2851).
#[test]
fn formatting_strings_are_read_as_ncbi_reads_them() {
    assert_eq!(
        parse_formatting_string(" 6 ").unwrap(),
        FormatChoice {
            number: 6,
            spec: String::new(),
            delimiter: String::new(),
        }
    );
    assert_eq!(parse_formatting_string("+6").unwrap().number, 6);
    assert_eq!(parse_formatting_string("06").unwrap().normalized(), "6");
    assert!(parse_formatting_string("6 std").unwrap().custom());
    assert!(!parse_formatting_string("0 std").unwrap().custom());
    assert_eq!(parse_formatting_string("0 std").unwrap().normalized(), "0");
    let delimited = parse_formatting_string("7 delim=, qseqid  sseqid ").unwrap();
    assert_eq!(
        (delimited.spec.as_str(), delimited.delimiter.as_str()),
        ("qseqid  sseqid", ",")
    );
    assert_eq!(delimited.normalized(), "7 qseqid  sseqid");
    assert!(!parse_formatting_string("6 delim=").unwrap().custom());
    assert!(message(parse_formatting_string("6 delim").unwrap_err())
        .contains("Delimiter format is invalid"));
    assert_eq!(
        message(parse_formatting_string("abc").unwrap_err()),
        "1 LOSAT query/options error: 'abc' is not a valid output format\nPlease refer to the LOSAT user manual.\n"
    );
    assert_eq!(
        message(parse_formatting_string("22").unwrap_err()),
        "255 Error: Formatting choice is out of range\n"
    );
    assert_eq!(
        message(parse_formatting_string("-1").unwrap_err()),
        "255 Error: Formatting choice is out of range\n"
    );
}

const PROGRAMS: [Program; 5] = [
    Program::Blastn,
    Program::Blastp,
    Program::Blastx,
    Program::Tblastn,
    Program::Tblastx,
];

// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1023-1074: the numbers of `EOutputFormat`; 19 (AIRR) and 21 (FASTA)
// are igblastn's and magicblast's, 22 is `eEndValue`.
#[test]
fn report_formats_are_ncbis_format_numbers() {
    for number in (0..=18).chain([20]) {
        let format = ReportFormat::from_number(number).unwrap();
        assert_eq!(format.number(), number);
    }
    for number in [-1, AIRR_REARRANGEMENT, FASTA, END_VALUE, 99] {
        assert_eq!(ReportFormat::from_number(number), None, "{number}");
    }
    assert_eq!(ReportFormat::from_number(17), Some(ReportFormat::Sam));
    assert_eq!(
        ReportFormat::from_number(20),
        Some(ReportFormat::CommaSeparatedValuesWithHeader)
    );
}

// The one table of the formats LOSAT writes (`Program::writes`): 0, 6, 7, 8, 9 and 12 in
// every program (8 and 12: line OF track S, OF3b; 9: OF3c), and 5 in every program
// (BLASTN, BLASTP: OF4a; BLASTX, TBLASTN, TBLASTX: OF4b); 15 and 16 in BLASTP and BLASTN
// (OF5a), BLASTX (OF5b part 1), TBLASTN (OF5b part 2) and TBLASTX (OF5b part 3): every program;
// 13 and 14 in BLASTP and BLASTN (OF5c part 1), BLASTX, TBLASTN and TBLASTX (part 2): every
// program writes 0, 5 to 9 and 12 to 16.
// 17 (SAM) in BLASTN alone (OF6); 10 and 20 in every program with the shared tabular
// writer (line CLI C2, PRs 1 to 3).
#[test]
fn every_program_writes_its_formats_only() {
    for program in PROGRAMS {
        let written: Vec<i32> = (0..END_VALUE)
            .filter_map(ReportFormat::from_number)
            .filter(|&format| program.writes(format))
            .map(ReportFormat::number)
            .collect();
        // 11 (the BLAST archive): every program (OF8a parts 1 and 2: BLASTN, BLASTP; OF8b
        // parts 1 to 3: BLASTX, TBLASTN, TBLASTX).
        // 10 and 20 with the shared tabular writer (line CLI C2: every program); 18, the
        // organism report (line RG T2: every program); 1 to 4, the query-anchored views (line RG
        // T3a: BLASTP; T3b: BLASTN, TBLASTN; T3c: BLASTX, TBLASTX: every program).
        let expected: &[i32] = match program {
            Program::Blastn => &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
            ],
            _ => &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20,
            ],
        };
        assert_eq!(written, expected, "{program:?}");
    }
    assert_eq!(
        PROGRAMS.map(Program::name),
        ["BLASTN", "BLASTP", "BLASTX", "TBLASTN", "TBLASTX"]
    );
    assert_eq!(
        PROGRAMS.map(Program::is_sam),
        [true, false, false, false, false]
    );
    assert_eq!(
        PROGRAMS.map(Program::custom_fields),
        [true, true, true, true, true]
    );
    assert_eq!(
        PROGRAMS.map(Program::shared_tabular),
        [true, true, true, true, true]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452: the structured reports return before the result
// warnings.
#[test]
fn structured_formats_write_no_query_warnings() {
    let silent: Vec<i32> = (0..END_VALUE)
        .filter_map(ReportFormat::from_number)
        .filter(|format| !format.writes_query_warnings())
        .map(ReportFormat::number)
        .collect();
    assert_eq!(silent, [5, 8, 9, 11, 12, 13, 14, 15, 16, 17]);
}

#[test]
fn report_formats_defer_ncbis_own_errors() {
    let choice = |number| FormatChoice {
        number,
        spec: String::new(),
        delimiter: String::new(),
    };
    let file = Some(Path::new("o.json"));
    for program in PROGRAMS {
        let name = program.name();
        assert_eq!(
            report_format(&choice(0), program, None).unwrap(),
            Some(ReportFormat::Pairwise)
        );
        // NCBI rejects these later (`formatting_handler_check`, `xinclude_check`).
        for number in [AIRR_REARRANGEMENT, FASTA, JSON, XML2] {
            assert_eq!(
                report_format(&choice(number), program, None).unwrap(),
                None,
                "{name} {number}"
            );
        }
        assert_eq!(
            report_format(&choice(JSON), program, Some(Path::new("-"))).unwrap(),
            None
        );
        // Formats that NCBI writes and LOSAT does not, before the inputs are read.
        for number in (1..=16)
            .filter(|&number| !matches!(number, PAIRWISE | TABULAR | TABULAR_WITH_COMMENTS))
            .filter(|&number| {
                let format = ReportFormat::from_number(number).unwrap();
                !program.writes(format)
            })
            .chain([18, 20])
            .filter(|&number| {
                !ReportFormat::from_number(number).is_some_and(|format| program.writes(format))
            })
        {
            assert_eq!(
                report_format(&choice(number), program, file)
                    .unwrap_err()
                    .to_string(),
                format!("output format {number} is not supported by LOSAT's {name}")
            );
        }
        let sam = report_format(&choice(SAM), program, None);
        if program.is_sam() {
            assert_eq!(sam.unwrap(), Some(ReportFormat::Sam));
        } else {
            assert_eq!(sam.unwrap(), None);
        }
        // A delimiter: the shared tabular writer's (BLASTN, TBLASTN, TBLASTX), else LOSAT's
        // rejection.
        let delimited = parse_formatting_string("6 delim=, std").unwrap();
        if program.shared_tabular() {
            assert_eq!(
                report_format(&delimited, program, None).unwrap(),
                Some(ReportFormat::Tabular)
            );
        } else {
            assert!(report_format(&delimited, program, None)
                .unwrap_err()
                .to_string()
                .ends_with(&format!("is not supported by LOSAT's {name}")));
        }
    }
    // A number outside NCBI's enum (no `-outfmt` reads so) is LOSAT's rejection.
    assert!(report_format(&choice(END_VALUE), Program::Blastp, None).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886 and c++/src/app/blast/blast_app_util.cpp:888-899.
#[test]
fn ncbis_later_format_errors_are_its_messages() {
    let choice = |number| FormatChoice {
        number,
        spec: String::new(),
        delimiter: String::new(),
    };
    let options =
        |text: &str| format!("1 LOSAT query/options error: {text}\n{USER_MANUAL_FOOTER}\n");
    for program in PROGRAMS {
        let sam = formatting_handler_check(&choice(SAM), program);
        if program.is_sam() {
            assert!(sam.is_ok());
        } else {
            assert_eq!(
                message(sam.unwrap_err()),
                options("SAM format is only applicable to blastn")
            );
        }
        assert_eq!(
            message(formatting_handler_check(&choice(AIRR_REARRANGEMENT), program).unwrap_err()),
            options("AIRR rearrangement format is only applicable to igblastn")
        );
        assert_eq!(
            message(formatting_handler_check(&choice(FASTA), program).unwrap_err()),
            options("FASTA output format is only applicable to magicblast")
        );
        for number in (0..=16).chain([18, 20]) {
            assert!(formatting_handler_check(&choice(number), program).is_ok());
        }
    }
    assert_eq!(
        message(xinclude_check(&choice(JSON)).unwrap_err()),
        options("Please provide a file name for outfmt 13.")
    );
    assert_eq!(
        message(xinclude_check(&choice(XML2)).unwrap_err()),
        options("Please provide a file name for outfmt 14.")
    );
    for number in (0..=12).chain(15..=21) {
        assert!(xinclude_check(&choice(number)).is_ok());
    }
}

#[test]
fn seg_values_follow_ncbis_filtering_handler() {
    assert_eq!(parse_seg_option("no").unwrap(), SegSpec::No);
    assert!(message(parse_seg_option("12 2.2").unwrap_err())
        .contains("Invalid number of arguments to filtering option"));
    assert!(message(parse_seg_option("12  2.2 2.5").unwrap_err())
        .contains("Invalid number of arguments to filtering option"));
    assert!(message(parse_seg_option("x 2.2 2.5").unwrap_err())
        .contains("Invalid input for filtering parameters"));
    assert!(message(parse_seg_option("12 inf 2.5").unwrap_err())
        .contains("Invalid input for filtering parameters"));
    assert_eq!(
        parse_seg_option("+12 2.2 .5").unwrap(),
        SegSpec::WindowLocutHicut {
            window: 12,
            locut: 2.2,
            hicut: 0.5
        }
    );
}

#[test]
fn composition_modes_read_the_first_character() {
    for (value, mode) in [
        ("0", CompositionMode::NoCompositionBasedStats),
        ("F", CompositionMode::NoCompositionBasedStats),
        ("", CompositionMode::NoCompositionBasedStats),
        ("x", CompositionMode::NoCompositionBasedStats),
        ("4", CompositionMode::NoCompositionBasedStats),
        ("1", CompositionMode::CompositionBasedStats),
        ("D", CompositionMode::CompositionMatrixAdjust),
        ("t", CompositionMode::CompositionMatrixAdjust),
        ("2xyz", CompositionMode::CompositionMatrixAdjust),
        ("3", CompositionMode::CompoForceFullMatrixAdjust),
    ] {
        assert_eq!(
            parse_comp_based_stats(value, true, false).unwrap().0,
            mode,
            "{value}"
        );
    }
    assert!(parse_comp_based_stats("2U", true, false).unwrap().1);
    assert!(!parse_comp_based_stats("2U", false, false).unwrap().1);
    assert!(!parse_comp_based_stats("0u", true, false).unwrap().1);
    assert!(parse_comp_based_stats("2", true, true).is_err());
    assert!(parse_comp_based_stats("x", true, true).is_ok());
}

// NCBI BLAST+ 2.17.0 (the official Linux x64 release) with `CHUNK_SIZE`,
// `OVERLAP_CHUNK_SIZE`, `BATCH_SIZE` and `PRE_FETCH_SEQS_LIMIT` set to each value (env
// parity probe h_probe of 2026-11-12, REPORT.md "Values"): `NStr::StringToInt` with flags 0
// takes a sign and no spaces; a value that it cannot convert ends the program with
// `CATCH_ALL`'s report of its `CStringException`, exit status 255.
#[test]
fn environment_values_are_converted_as_ncbis_string_to_int() {
    use std::ffi::OsStr;
    for (text, value) in [("5", 5), ("+5", 5), ("-1", -1), ("0", 0)] {
        assert_eq!(
            environment_int(OsStr::new(text)).unwrap(),
            value,
            "{text:?}"
        );
    }
    let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
    for (text, line, function, kind, pos) in [
        ("abc", 861, "StringToInt8", "Int8", 0),
        (" 5", 861, "StringToInt8", "Int8", 0),
        ("", 861, "StringToInt8", "Int8", 0),
        ("1.5", 870, "StringToInt8", "Int8", 1),
        ("5x", 870, "StringToInt8", "Int8", 1),
        ("5 ", 870, "StringToInt8", "Int8", 1),
        (
            "99999999999999999999",
            852,
            "StringToInt8",
            "Int8, overflow",
            18,
        ),
        (
            "-99999999999999999999",
            852,
            "StringToInt8",
            "Int8, overflow",
            19,
        ),
        ("2147483648", 640, "StringToInt", "int, overflow", 0),
        ("-2147483649", 640, "StringToInt", "int, overflow", 0),
    ] {
        let error = environment_int(OsStr::new(text)).unwrap_err();
        let native = error.downcast_ref::<NativeError>().expect("NCBI's report");
        assert_eq!(native.exit, 255, "{text:?}");
        assert_eq!(
            native.message,
            // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
            crate::wording::losat_text(&format!(
                "Error: NCBI C++ Exception:\n    T0 \"{root}/src/corelib/ncbistr.cpp\", line {line}: Error: (CStringException::eConvert) ncbi::NStr::{function}() - Cannot convert string '{text}' to {kind} (m_Pos = {pos})\n\n"
            ))
            .as_bytes(),
            "{text:?}"
        );
    }
}

/// A format NCBI cannot read is left to the formatting handler when a strategy file is asked
/// for (`early_format_choice`), which reads it again and fails there; a readable one is as it
/// was read, and without a strategy file the error comes at once.
#[test]
fn an_unreadable_format_is_deferred_with_a_strategy_file() {
    for outfmt in ["99", "abc", "6 delimx", " -1"] {
        let message = |error: anyhow::Error| error.to_string();
        assert!(early_format_choice(outfmt, false).is_err(), "{outfmt}");
        let choice = early_format_choice(outfmt, true).unwrap();
        assert_eq!(choice, EarlyFormat::Deferred, "{outfmt}");
        assert_eq!(choice.number(), None, "{outfmt}");
        assert!(choice
            .report_format(Program::Blastp, None)
            .unwrap()
            .is_none());
        assert_eq!(
            message(choice.at_formatting_handler(outfmt).unwrap_err()),
            message(parse_formatting_string(outfmt).unwrap_err()),
            "{outfmt}"
        );
    }
    let choice = early_format_choice("6 qseqid", true).unwrap();
    let parsed = choice.parsed().unwrap().clone();
    assert_eq!(choice.number(), Some(6));
    assert_eq!(choice.at_formatting_handler("99").unwrap(), parsed);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1535-1553
#[test]
fn the_reports_of_the_formats_0_to_4_read_the_html_links() {
    // NCBI calls `CAlignFormatUtil::InitConfig` only from the id links of the description table,
    // which formats 0 to 4 write alike (the query-anchored ones too, item T4b of line RG); the
    // organism report (18) and the formats that ignore -html read no registry file.
    for number in 0..END_VALUE {
        let links = (0..=4).contains(&number);
        assert_eq!(formats_link_ids(&[number]), links, "outfmt {number}");
        assert_eq!(
            formats_link_ids(&[6, number, 18]),
            links,
            "outfmt 6, {number}, 18"
        );
    }
    for formats in [&[18][..], &[6, 18], &[5], &[]] {
        assert!(html_links(true, Program::Blastp, None, formats)
            .unwrap()
            .is_none());
    }
    for formats in [&[0][..], &[1], &[2], &[3], &[4]] {
        assert!(html_links(false, Program::Blastn, Some("db"), formats)
            .unwrap()
            .is_none());
    }
}
