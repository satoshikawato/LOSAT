use super::*;

fn counts(
    outfmt: i32,
    mts: Option<usize>,
    nd: Option<usize>,
    na: Option<usize>,
) -> (usize, usize, usize) {
    let c = report_counts(outfmt, mts, nd, na);
    (c.hitlist_size, c.num_descriptions, c.num_alignments)
}

#[test]
fn format_up_to_4_takes_the_defaults_and_the_larger_given_number() {
    assert_eq!(counts(0, None, None, None), (500, 500, 250));
    assert_eq!(counts(0, Some(7), None, None), (7, 7, 7));
    assert_eq!(counts(0, Some(500), None, None), (500, 500, 500));
    // -num_descriptions 3: the hit list is max(3, 250).
    assert_eq!(counts(0, None, Some(3), None), (250, 3, 250));
    // -num_alignments 2: max(500, 2).
    assert_eq!(counts(0, None, None, Some(2)), (500, 500, 2));
    assert_eq!(counts(0, None, Some(2), Some(4)), (4, 2, 4));
    assert_eq!(counts(0, None, Some(4), Some(2)), (4, 4, 2));
    assert_eq!(counts(0, None, Some(600), None), (600, 600, 250));
    assert_eq!(counts(0, None, None, Some(600)), (600, 500, 600));
    assert_eq!(counts(0, None, Some(0), None), (250, 0, 250));
    assert_eq!(counts(0, None, Some(0), Some(0)), (0, 0, 0));
    assert_eq!(counts(4, None, Some(1), Some(1)), (1, 1, 1));
}

#[test]
fn format_above_4_takes_the_target_count_or_the_alignment_count() {
    for outfmt in [5, 6, 7, 8, 12, 20] {
        assert_eq!(counts(outfmt, None, None, None), (500, 500, 500));
        assert_eq!(counts(outfmt, Some(9), None, None), (9, 9, 9));
        // -num_descriptions is ignored (and warned about).
        assert_eq!(counts(outfmt, None, Some(3), None), (500, 500, 500));
        assert_eq!(counts(outfmt, None, None, Some(2)), (2, 2, 2));
        assert_eq!(counts(outfmt, None, Some(3), Some(2)), (2, 2, 2));
        assert_eq!(counts(outfmt, None, None, Some(0)), (0, 0, 0));
    }
}

#[test]
fn warnings_come_in_ncbis_order() {
    let c = report_counts(6, None, Some(3), Some(3));
    assert_eq!(
        String::from_utf8(c.warnings("blastp", 6, true)).unwrap(),
        "Warning: [blastp] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\nWarning: [blastp] Examining 5 or more matches is recommended\n"
    );
    // The hit list size, not -num_descriptions, decides the second warning.
    let c = report_counts(0, None, Some(3), None);
    assert!(c.warnings("blastn", 0, true).is_empty());
    let c = report_counts(7, None, Some(3), None);
    assert_eq!(
        String::from_utf8(c.warnings("tblastx", 7, true)).unwrap(),
        "Warning: [tblastx] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n"
    );
    let c = report_counts(0, None, None, Some(3));
    assert!(c.warnings("blastp", 0, false).is_empty());
    let c = report_counts(6, None, None, Some(3));
    assert_eq!(
        String::from_utf8(c.warnings("blastp", 6, false)).unwrap(),
        "Warning: [blastp] Examining 5 or more matches is recommended\n"
    );
    let c = report_counts(0, Some(4), None, None);
    assert_eq!(
        String::from_utf8(c.warnings("blastx", 0, false)).unwrap(),
        "Warning: [blastx] Examining 5 or more matches is recommended\n"
    );
}

// blast_args.cpp:2947-2989 (cited at `warnings_with`) and NCBI 2.17.0's rpsblast
// (DESIGN_ps11c `f.allfmt.c0.7.err`): -num_descriptions, -line_length, -sorthits, the few
// matches, then -sorthsps; outfmt 0 reads -line_length and -sorthits and warns of nothing.
#[test]
fn the_formatting_options_warn_in_ncbis_order() {
    let all = FormattingArgs {
        line_length: Some(80),
        sort_hits: Some(1),
        sort_hsps: Some(1),
        hitlist_start: None,
    };
    let (counts, warnings) =
        run_report_counts_with("rpsblast", &[7], Some(2), Some(3), None, &all).unwrap();
    assert_eq!(counts.hitlist_size, 2);
    assert_eq!(
        String::from_utf8(warnings).unwrap(),
        "Warning: [rpsblast] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n\
             Warning: [rpsblast] The parameter -line_length is not applicable for output formats > 4 .\n\
             Warning: [rpsblast] The parameter -sorthits is ignored for output formats > 4.\n\
             Warning: [rpsblast] Examining 5 or more matches is recommended\n\
             Warning: [rpsblast] The parameter -sorthsps is ignored for output formats != 0.\n"
    );
    // Outfmt 0: only the few matches.
    let (_, warnings) =
        run_report_counts_with("rpsblast", &[0], Some(2), None, None, &all).unwrap();
    assert_eq!(
        String::from_utf8(warnings).unwrap(),
        "Warning: [rpsblast] Examining 5 or more matches is recommended\n"
    );
    // -sorthsps alone with outfmt 6, and -line_length 1 with outfmt 10.
    let hsps = FormattingArgs {
        sort_hsps: Some(4),
        ..FormattingArgs::default()
    };
    assert_eq!(
        String::from_utf8(
            report_counts(6, None, None, None).warnings_with("rpsblast", 6, false, &hsps)
        )
        .unwrap(),
        "Warning: [rpsblast] The parameter -sorthsps is ignored for output formats != 0.\n"
    );
    let length = FormattingArgs {
        line_length: Some(1),
        ..FormattingArgs::default()
    };
    assert_eq!(
        String::from_utf8(report_counts(10, None, None, None).warnings_with("rpsblast", 10, false, &length))
            .unwrap(),
        "Warning: [rpsblast] The parameter -line_length is not applicable for output formats > 4 .\n"
    );
    assert_eq!(length.line_length(), 1);
    assert_eq!(FormattingArgs::default().line_length(), 60);
}

#[test]
fn a_hit_list_of_zero_is_the_no_hits_error() {
    assert!(report_counts(0, None, Some(1), None).validate().is_ok());
    let err = report_counts(0, None, Some(0), Some(0))
        .validate()
        .unwrap_err();
    assert_eq!(
        format!("{err:#}"),
        "LOSAT query/options error: No hits are being saved\nPlease refer to the LOSAT user manual.\n"
    );
    assert!(report_counts(6, None, None, Some(0)).validate().is_err());
    assert!(report_counts(0, None, Some(0), None).validate().is_ok());
}

#[test]
fn formats_share_a_run_only_with_one_hit_list_size() {
    // -num_descriptions 3: the hit list is 250 for format 0 and 500 for format 6.
    assert!(run_report_counts("blastp", &[0, 6], None, Some(3), None).is_err());
    assert!(run_report_counts("blastp", &[6, 0], None, Some(3), None).is_err());
    // Neither option: format 0 shows 500 and 250, the tabular formats read the hit list.
    let (counts, warnings) = run_report_counts("blastp", &[6, 0], None, None, None).unwrap();
    assert_eq!(counts, ReportCounts::default());
    assert!(warnings.is_empty());
    // -num_descriptions 3 -num_alignments 500: one hit list size of 500.
    let (counts, warnings) =
        run_report_counts("blastp", &[6, 0], None, Some(3), Some(500)).unwrap();
    assert_eq!(
        (
            counts.hitlist_size,
            counts.num_descriptions,
            counts.num_alignments
        ),
        (500, 3, 500)
    );
    assert_eq!(
        String::from_utf8(warnings).unwrap(),
        "Warning: [blastp] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n"
    );
    let (counts, warnings) = run_report_counts("blastn", &[6, 7], Some(2), None, None).unwrap();
    assert_eq!(counts.hitlist_size, 2);
    assert_eq!(
        String::from_utf8(warnings).unwrap(),
        "Warning: [blastn] Examining 5 or more matches is recommended\n"
    );
    assert!(run_report_counts("blastn", &[], None, None, None).is_err());
}

/// A hit list size above 500 in the options is kept by the pairwise formats when the command line
/// gives neither `-num_descriptions` nor `-num_alignments` (the descriptions the size, the
/// alignments half of it), whatever `-max_target_seqs` says; otherwise the counts go on as for
/// the default start (probes `Y.*.imp.*`).
#[test]
fn a_hit_list_above_500_in_the_options_is_kept_by_the_pairwise_formats() {
    let kept = ReportCounts {
        hitlist_size: 600,
        num_descriptions: 600,
        num_alignments: 300,
    };
    for outfmt in 0..=4 {
        assert_eq!(kept_hitlist_counts(600, outfmt, false, false), Some(kept));
        assert_eq!(report_counts_from(600, outfmt, Some(501), None, None), kept);
    }
    assert_eq!(kept_hitlist_counts(500, 0, false, false), None);
    assert_eq!(kept_hitlist_counts(600, 5, false, false), None);
    assert_eq!(kept_hitlist_counts(600, 0, true, false), None);
    assert_eq!(kept_hitlist_counts(600, 0, false, true), None);
    assert_eq!(
        report_counts_from(600, 0, None, None, Some(300)),
        ReportCounts {
            hitlist_size: 500,
            num_descriptions: 500,
            num_alignments: 300,
        }
    );
}
