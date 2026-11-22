//! `extract_formatting` and `sam_check` with NCBI BLAST+ 2.17.0 blast_formatter's texts (line RG
//! item F3: probes P2.149-P2.200 of `-outfmt`, P2.001-P2.121 of the counts and the warnings).

use super::*;

fn options(outfmt: &str) -> Options {
    Options {
        outfmt: outfmt.to_string(),
        show_gis: false,
        num_descriptions: None,
        num_alignments: None,
        line_length: None,
        html: false,
        sort_hits: None,
        sort_hsps: None,
        max_target_seqs: None,
        out: "-".to_string(),
        parse_deflines: false,
    }
}

/// The exit status and text of `CATCH_ALL` for an error of `extract_formatting`.
fn caught(error: anyhow::Error) -> (i32, String) {
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .expect("NCBI's error");
    (
        native.exit,
        String::from_utf8(native.message.clone()).unwrap(),
    )
}

fn error_of(outfmt: &str) -> (i32, String) {
    caught(extract_formatting(&options(outfmt), 500).unwrap_err())
}

fn options_error_text(message: &str) -> (i32, String) {
    caught(crate::blastinput::app::options_error(message))
}

#[test]
fn every_format_of_ncbi_s_list_is_read() {
    for number in (0..=18).chain([20]) {
        let formatting = extract_formatting(&options(&number.to_string()), 500).unwrap();
        assert_eq!(formatting.choice.number, number);
        assert!(!formatting.choice.custom(), "{number}");
        assert!(formatting.warnings.is_empty(), "{number}");
    }
    // `NStr::StringToInt` after `TruncateSpaces` (P2.171, P2.173, P2.185, P2.187).
    for text in [" 6", "6 ", "06", "+6"] {
        assert_eq!(
            extract_formatting(&options(text), 500)
                .unwrap()
                .choice
                .number,
            6,
            "{text:?}"
        );
    }
}

#[test]
fn the_errors_of_outfmt_are_ncbi_s() {
    // P2.149, P2.183, P2.189, P2.191.
    for text in ["abc", "1e1", "6.0", ""] {
        assert_eq!(
            error_of(text),
            options_error_text(&format!("'{text}' is not a valid output format")),
            "{text:?}"
        );
    }
    // P2.151, P2.157, P2.159: `std::out_of_range`, exit status 255.
    for text in ["-1", "22", "100"] {
        assert_eq!(
            error_of(text),
            (
                255,
                "Error: Formatting choice is out of range\n".to_string()
            ),
            "{text:?}"
        );
    }
    // P2.153, P2.155: the formats of igblastn and magicblast.
    assert_eq!(
        error_of("19"),
        options_error_text("AIRR rearrangement format is only applicable to igblastn")
    );
    assert_eq!(
        error_of("21"),
        options_error_text("FASTA output format is only applicable to magicblast")
    );
    // P2.181, P2.163, P2.167: the delimiter.
    assert_eq!(
        error_of("6 delim"),
        options_error_text("Delimiter format is invalid. Valid format is delim=<delimiter value>")
    );
    for (outfmt, delim, field) in [
        ("6 delim=; staxids", ";", "staxids"),
        ("6 delim=<> salltitles", "<>", "salltitles"),
    ] {
        assert_eq!(
            error_of(outfmt),
            options_error_text(&format!(
                "Your custom record separator ({delim}) is also used by the format specifier \
                 ({field}) to separate multiple entries. Please use a different record separator \
                 (delim keyword)."
            )),
            "{outfmt}"
        );
    }
}

#[test]
fn the_custom_specifications_are_the_search_programs() {
    let choice = |text: &str| extract_formatting(&options(text), 500).unwrap().choice;
    // P2.161, P2.175, P2.177, P2.179.
    let custom = choice("6 delim=");
    assert_eq!((custom.spec.as_str(), custom.delimiter.as_str()), ("", ""));
    assert_eq!(choice("6 qseqid foo").spec, "qseqid foo");
    let custom = choice("7 delim=, qseqid");
    assert_eq!(
        (
            custom.number,
            custom.spec.as_str(),
            custom.delimiter.as_str()
        ),
        (7, "qseqid", ",")
    );
    assert_eq!(choice("6 delim=: qseqid sseqid").delimiter, ":");
    // P2.169: a format that takes no fields drops them.
    assert!(!choice("0 qseqid").custom());
    // SAM (eIsSAM) keeps its words (P2.193, P2.195).
    assert_eq!(choice("17").number, 17);
    assert_eq!(choice("17 SQ").spec, "SQ");
}

// blast_formatter.cpp:350-355: SAM for blastn's archives only (P2.197, P2.199: blastp).
#[test]
fn sam_is_for_blastn_archives() {
    let sam = extract_formatting(&options("17 SQ"), 500).unwrap().choice;
    assert!(sam_check(&sam, true).is_ok());
    assert_eq!(
        caught(sam_check(&sam, false).unwrap_err()),
        options_error_text("SAM format is only applicable to blastn results")
    );
    let tabular = extract_formatting(&options("6"), 500).unwrap().choice;
    assert!(sam_check(&tabular, false).is_ok());
}

/// The counts of `options` from the hit list size `hitlist_size`.
fn counts(options: &Options, hitlist_size: usize) -> (usize, usize, usize) {
    let counts = extract_formatting(options, hitlist_size).unwrap().counts;
    (
        counts.hitlist_size,
        counts.num_descriptions,
        counts.num_alignments,
    )
}

fn warnings(options: &Options, hitlist_size: usize) -> String {
    String::from_utf8(extract_formatting(options, hitlist_size).unwrap().warnings).unwrap()
}

#[test]
fn the_counts_start_from_the_object_s_hit_list_size() {
    let mut pairwise = options("0");
    assert_eq!(counts(&pairwise, 500), (500, 500, 250));
    assert_eq!(counts(&pairwise, 100), (100, 500, 250));
    pairwise.max_target_seqs = Some(1);
    assert_eq!(counts(&pairwise, 500), (1, 1, 1));
    // P2.016: the few matches, without a prefix.
    assert_eq!(
        warnings(&pairwise, 500),
        "Warning: Examining 5 or more matches is recommended\n"
    );
    let mut tabular = options("6");
    assert_eq!(counts(&tabular, 100), (100, 100, 100));
    tabular.num_descriptions = Some(1);
    // P2.002: -num_descriptions is ignored above 4.
    assert_eq!(counts(&tabular, 500), (500, 500, 500));
    assert_eq!(
        warnings(&tabular, 500),
        "Warning: The parameter -num_descriptions is ignored for output formats > 4 . Use \
         -max_target_seqs to control output\n"
    );
    tabular.num_descriptions = None;
    tabular.num_alignments = Some(0);
    assert_eq!(counts(&tabular, 500), (0, 0, 0));
}

// blast_args.cpp:2894-2908 (quoted at `report_counts::kept_hitlist_counts`): above 500, the
// formats up to 4 keep the object's size and read nothing else.
#[test]
fn a_hit_list_above_500_is_kept_for_the_pairwise_formats() {
    let mut pairwise = options("0");
    pairwise.max_target_seqs = Some(2);
    pairwise.line_length = Some(30);
    pairwise.sort_hits = Some(1);
    pairwise.sort_hsps = Some(1);
    let formatting = extract_formatting(&pairwise, 1000).unwrap();
    assert_eq!(
        (
            formatting.counts.hitlist_size,
            formatting.counts.num_descriptions,
            formatting.counts.num_alignments
        ),
        (1000, 1000, 500)
    );
    assert_eq!(formatting.formatting, FormattingArgs::default());
    assert!(formatting.warnings.is_empty());
    // -num_alignments reads the command line again.
    pairwise.max_target_seqs = None;
    pairwise.num_alignments = Some(3);
    assert_eq!(counts(&pairwise, 1000), (500, 500, 3));
    // Above 4 the size is read as below 500.
    assert_eq!(counts(&options("6"), 1000), (1000, 1000, 1000));
}

// P2.088-P2.117 and blast_args.cpp:2947-2989 (quoted at `ReportCounts::warnings_with`).
#[test]
fn the_formatting_options_are_kept_where_ncbi_reads_them() {
    let mut given = options("0");
    given.line_length = Some(30);
    given.sort_hits = Some(2);
    given.sort_hsps = Some(1);
    let formatting = extract_formatting(&given, 500).unwrap();
    assert_eq!(
        formatting.formatting,
        FormattingArgs {
            line_length: Some(30),
            sort_hits: Some(2),
            sort_hsps: Some(1),
            hitlist_start: None,
        }
    );
    assert!(formatting.warnings.is_empty());
    given.outfmt = "1".to_string();
    let formatting = extract_formatting(&given, 500).unwrap();
    assert_eq!(formatting.formatting.sort_hsps, None);
    assert_eq!(
        String::from_utf8(formatting.warnings).unwrap(),
        "Warning: The parameter -sorthsps is ignored for output formats != 0.\n"
    );
    given.outfmt = "7".to_string();
    given.num_descriptions = Some(3);
    let formatting = extract_formatting(&given, 500).unwrap();
    assert_eq!(formatting.formatting, FormattingArgs::default());
    assert_eq!(
        String::from_utf8(formatting.warnings).unwrap(),
        "Warning: The parameter -num_descriptions is ignored for output formats > 4 . Use \
         -max_target_seqs to control output\n\
         Warning: The parameter -line_length is not applicable for output formats > 4 .\n\
         Warning: The parameter -sorthits is ignored for output formats > 4.\n\
         Warning: The parameter -sorthsps is ignored for output formats != 0.\n"
    );
    let mut flags = options("8");
    flags.show_gis = true;
    flags.html = true;
    let formatting = extract_formatting(&flags, 500).unwrap();
    assert!(formatting.show_gis && formatting.html);
}
