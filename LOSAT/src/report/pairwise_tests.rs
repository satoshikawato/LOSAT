use super::*;

// No NCBI counterpart: the size of LOSAT's record of an outfmt 0 HSP, one per reported HSP
// until the report is written. It is not larger than on dev ed9178f49 (376 bytes on 64-bit
// targets): the masks of the rows are not kept per HSP (line RG T4a1; the rows hold them in
// lower case, `query_anchored::output_seq`). A boxed list per masked HSP had made it 384 bytes
// and raised max RSS by up to 2 % (3,000 BLASTX HSPs). On 32-bit targets (wasm32) nothing is
// pinned.
#[cfg(target_pointer_width = "64")]
#[test]
fn pairwise_hit_is_no_larger_than_on_dev() {
    let size = std::mem::size_of::<PairwiseHit>();
    assert!(size <= 376, "{size}");
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:533-538
// (cited at `write_blastp_database_header_spacing`)
// The `Database:` lines of a FASTA subject set print NCBI's `Int8` counts whole: a set of
// 2^32 letters or more (two records of 3,000,000,000 letters, each below 2^32 as on wasm32) is
// not cut to a 32-bit `usize`.
#[test]
fn subject_set_database_lines_print_totals_of_2_to_the_32_letters_or_more() {
    let mut prolog = Vec::new();
    write_blastp_pairwise_prolog(
        &mut prolog,
        "2.17.0+",
        "BLASTP",
        b"s",
        2,
        6_000_000_000,
        None,
        None,
    )
    .unwrap();
    let prolog = String::from_utf8(prolog).unwrap();
    assert!(
        prolog.contains("Database: s.\n           2 sequences; 6,000,000,000 total letters\n"),
        "{prolog}"
    );
    let mut epilog = Vec::new();
    write_final_database_report(&mut epilog, b"s", 2, 6_000_000_000).unwrap();
    let epilog = String::from_utf8(epilog).unwrap();
    assert!(
        epilog.contains("  Number of letters in database: 6,000,000,000\n  Number of sequences in database:  2\n"),
        "{epilog}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:395 (cited at
// `write_prolog`): the prolog's two `NcbiEndl`, each a newline and a flush.
// Every program's prolog makes NCBI's flush after each of these two newlines
// (`ReportWrite::ncbi_flush`, line RG item T6a; line RG T4a1's review note N3 had kept each
// program's flushes of dev): blastp's driver (blastp, psiblast, rpsblast), megablast and the
// other blastn tasks, blastx, rpstblastn, tblastn (also with `-in_pssm`) and tblastx, the plain
// prolog and the `-html` one (line RG T4a2 review B4: blastp's `-html -out /dev/full` with
// `GETURL_DEBUG` wrote `REGISTRY: .ncbirc` before NCBI's write error).
#[test]
fn every_prolog_flushes_after_its_two_newlines() {
    /// Keeps what is written and the number of bytes written at each flush of NCBI's, and at
    /// each of LOSAT's (`std::io::Write::flush`).
    #[derive(Default)]
    struct Flushes(Vec<u8>, Vec<usize>, Vec<usize>);
    impl Write for Flushes {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.0.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            self.2.push(self.0.len());
            Ok(())
        }
    }
    impl crate::report::ostream::ReportWrite for Flushes {
        fn ncbi_flush(&mut self) -> io::Result<()> {
            self.1.push(self.0.len());
            Ok(())
        }
    }
    type Html<'h> = Option<(&'h crate::report::html_ids::HtmlLinks, &'h mut dyn Write)>;
    type Prolog<'a> = Box<dyn for<'h> Fn(&mut Flushes, Html<'h>) -> io::Result<()> + 'a>;
    let version = "2.17.0+";
    let report = BlastpPairwiseReport {
        version: version.into(),
        database_name: "d".into(),
        database_num_sequences: 1,
        database_total_letters: 10,
        matrix_name: Some("BLOSUM62".into()),
        gap_open: 11,
        gap_extend: 1,
        word_threshold: 12.0,
        window_size: 40,
        gapped_karlin: Default::default(),
        gumbel: None,
        gapped: true,
        show_linked_set_size: false,
        num_descriptions: 500,
        num_alignments: 250,
    };
    let blastp = |program: &'static str| -> Prolog<'_> {
        Box::new(move |out: &mut Flushes, html: Html<'_>| {
            write_blastp_pairwise_prolog(out, version, program, b"d", 1, 10, None, html)
        })
    };
    let blastn = |megablast: bool| -> Prolog<'_> {
        Box::new(move |out: &mut Flushes, html: Html<'_>| {
            let intro = BlastnIntro {
                megablast,
                indexed_megablast: false,
                matrix_only: false,
            };
            write_blastn_pairwise_prolog(out, version, intro, b"d", 1, 10, None, html)
        })
    };
    let tblastn = |program: &'static str| -> Prolog<'_> {
        Box::new(move |out: &mut Flushes, html: Html<'_>| {
            write_tblastn_pairwise_prolog(out, program, version, b"d", 1, 10, None, html)
        })
    };
    let blastx = |program: &'static str| -> Prolog<'_> {
        let report = &report;
        Box::new(move |out: &mut Flushes, html: Html<'_>| {
            write_blastx_prolog(out, program, report, None, html)
        })
    };
    let prologs: [(&str, Prolog<'_>); 10] = [
        ("blastp", blastp("blastp")),
        ("psiblast", blastp("psiblast")),
        ("rpsblast", blastp("rpsblast")),
        ("megablast", blastn(true)),
        ("blastn", blastn(false)),
        (
            "tblastx",
            Box::new(|out: &mut Flushes, html: Html<'_>| {
                write_tblastx_pairwise_prolog(out, version, b"d", 1, 10, None, html)
            }),
        ),
        ("tblastn", tblastn("TBLASTN")),
        ("psitblastn", tblastn("PSITBLASTN")),
        ("blastx", blastx("BLASTX")),
        ("rpstblastn", blastx("RPSTBLASTN")),
    ];
    let links = html_links_with(None);
    for (program, prolog) in &prologs {
        let mut out = Flushes::default();
        prolog(&mut out, None).unwrap();
        // The version line, then the two newlines, each flushed.
        let line = out.0.iter().position(|&byte| byte == b'\n').unwrap() + 1;
        assert_eq!(&out.0[line..line + 2], b"\n\n", "{program}");
        assert_eq!(out.1, vec![line + 1, line + 2], "{program}");
        // `-html`: the bold version line after the HTML head, then the two newlines, each flushed
        // (no `./.ncbirc`: the references post nothing and flush nowhere).
        let mut out = Flushes::default();
        let mut diagnostics = Vec::new();
        prolog(&mut out, Some((&links, &mut diagnostics))).unwrap();
        let line = HTML_PREFIX.len()
            + out.0[HTML_PREFIX.len()..]
                .iter()
                .position(|&byte| byte == b'\n')
                .unwrap()
            + 1;
        assert!(out.0[..line].ends_with(b"</b>\n"), "{program}");
        assert_eq!(&out.0[line..line + 2], b"\n\n", "{program}");
        assert_eq!(out.1, [line + 1, line + 2], "{program}");
        assert!(out.2.is_empty(), "{program}");
        assert!(diagnostics.is_empty());
    }
    // A reference that posts hands the report written so far to its stream first (LOSAT's flush,
    // `query_warnings::write_html_post`), where NCBI's `cerr`, tied to `cout`, flushes it
    // (`ostream::cerr_post`), so the posts follow the bytes before them; NCBI's own flushes are
    // the `NcbiEndl` pair alone.
    let posting = html_links_with(Some(b"[BLASTFMTUTIL]\nPROTOCOL=https:\nPROTOCOL=https:\n"));
    let mut out = Flushes::default();
    let mut diagnostics = Vec::new();
    (prologs[0].1)(&mut out, Some((&posting, &mut diagnostics))).unwrap();
    assert_eq!(
        String::from_utf8(diagnostics).unwrap(),
        "Warning: (110.8) [blastp] Found multiple [BLASTFMTUTIL] PROTOCOL settings; using the one from line 3\n"
            .repeat(2)
    );
    assert_eq!(out.1.len(), 2);
    assert_eq!(out.2.len(), 2);
}

/// `write_flatfile_wrapped_with` finds the break of a line without newlines or backspaces from the window's
/// end; it equals the scan of every column (the form it takes for a text that may hold a
/// backspace) on texts of every scored class, at every width.
#[test]
fn the_break_found_from_the_end_is_the_scans() {
    let alphabet = b"ab ,-\n\t\x0bxyz  \x08";
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for _ in 0..20_000 {
        let len = (next() % 160) as usize;
        let backspace = next() % 4 == 0;
        let text: Vec<u8> = (0..len)
            .map(|_| alphabet[(next() % alphabet.len() as u64) as usize])
            .map(|c| if c == 8 && !backspace { b'q' } else { c })
            .collect();
        let width = 1 + (next() % 70) as usize;
        let mut found = Vec::new();
        write_flatfile_wrapped(&mut found, &text, width).unwrap();
        let mut scanned = Vec::new();
        let every_column = TextFacts {
            newlines: true,
            backspaces: true,
        };
        write_flatfile_wrapped_with(&mut scanned, &text, width, every_column).unwrap();
        assert_eq!(
            found,
            scanned,
            "{:?} at width {width}",
            String::from_utf8_lossy(&text)
        );
    }
}

/// The plain description label is the id text and the title as one string, cut to 68 columns
/// or padded to them, wherever the parts split it (line RG, item T4a0).
#[test]
fn a_description_label_in_parts_writes_the_joined_line() {
    let cut = |id: &[u8], title: &[u8]| {
        let mut out = Vec::new();
        DescriptionLabel {
            id: id.to_vec(),
            title: title.to_vec(),
            caught: false,
        }
        .write_cut(&mut out)
        .unwrap();
        String::from_utf8(out).unwrap()
    };
    let joined = |label: &[u8]| cut(b"", label);
    assert_eq!(cut(b"P1 ", b"title"), format!("P1 title{}", " ".repeat(60)));
    assert_eq!(joined(b"P1 title"), cut(b"P1 ", b"title"));
    let long = [b'x'; 70];
    for split in [0, 3, 64, 65, 66, 70] {
        let line = cut(&long[..split], &long[split..]);
        assert_eq!(line, format!("{}...", "x".repeat(65)), "split {split}");
    }
    let exact = [b'y'; 68];
    assert_eq!(cut(&exact[..10], &exact[10..]), "y".repeat(68));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2266-2268 (cited at
// `write_matrix_line`)
// The matrix line of a named matrix; without one (RPS-BLAST before a query batch's set-up) the
// line stops after `Matrix: ` with the stream's failure, which a program reports as NCBI's
// "BLAST failed to write output" (NCBI BLAST+ 2.17.0: rpsblast and rpstblastn, an empty pipe,
// outfmt 0: stdout ends in `Matrix: `, exit status 6; rows `db.rps.stdin.empty.0`,
// `db.rpstbn.stdin.empty.*`).
#[test]
fn a_null_matrix_name_stops_the_epilog_after_its_label() {
    let mut named = Vec::new();
    write_matrix_line(&mut named, Some("BLOSUM62")).unwrap();
    assert_eq!(named, b"Matrix: BLOSUM62\n");
    let mut unnamed = Vec::new();
    let error = write_matrix_line(&mut unnamed, None).unwrap_err();
    assert_eq!(unnamed, b"Matrix: ");
    assert!(crate::report::ostream::is_report_write_failure(&error));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:642-649 (`s_SetFlags`,
// cited at `report::query_anchored`) and c++/src/algo/blast/blastinput/blast_args.cpp:2980-2987
// (cited at `PairwiseView::report_sort`).
#[test]
fn the_views_of_formats_0_to_4_and_their_sorts() {
    use crate::blastinput::app::ReportFormat;
    let sort = crate::report::hit_sort::ReportSort {
        hits: Some(3),
        hsps: Some(1),
    };
    for number in 0..=4 {
        let view = PairwiseView::of(ReportFormat::from_number(number).unwrap()).unwrap();
        assert_eq!(view.outfmt(), number);
        let expected = match number {
            0 => sort,
            _ => crate::report::hit_sort::ReportSort {
                hits: Some(3),
                hsps: None,
            },
        };
        assert_eq!(view.report_sort(sort), expected);
    }
    assert_eq!(
        PairwiseView::of(ReportFormat::QueryAnchoredIdentities),
        Some(PairwiseView::QueryAnchored {
            identity: true,
            anchored: true
        })
    );
    assert_eq!(
        PairwiseView::of(ReportFormat::FlatQueryAnchoredNoIdentities),
        Some(PairwiseView::QueryAnchored {
            identity: false,
            anchored: false
        })
    );
    assert_eq!(PairwiseView::of(ReportFormat::Tabular), None);
}

fn long_titled_database() -> PairwiseDatabase {
    PairwiseDatabase::new(vec![DatabaseInfo {
        definition: [
            &b"a database title that is much longer than sixty-eight letters, with commas, "[..],
            b"hyphen-ated words and a-very-long-unbroken-token-",
            &[b'z'; 80],
            b" at its end",
        ]
        .concat(),
        date: "Oct 13, 2026  1:18 PM".to_string(),
        num_sequences: 17,
        total_letters: 2287,
        filt_algorithm_name: String::new(),
        filt_algorithm_options: String::new(),
    }])
}

// NCBI BLAST+ 2.17.0 `blastp -db longtitle -query hba.faa -num_threads 1` (campaign
// dbsearch-20261012, scripts/db3_compare.py case shapes.longtitle.0): the prolog's
// `Database:` lines. The title is wrapped at 68 without the period that the text of
// a FASTA subject set gets.
#[test]
fn database_header_is_the_stored_title_wrapped_and_the_real_counts() {
    let mut written = Vec::new();
    write_blast_database_header(&mut written, &long_titled_database()).unwrap();
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "Database: a database title that is much longer than sixty-eight letters, with\n\
             commas, hyphen-ated words and\n\
             a-very-long-unbroken-token-\n\
             zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz\n\
             zzzzzzzzzzzz at its end\n           17 sequences; 2,287 total letters\n\n"
    );
    let mut written = Vec::new();
    write_blast_database_header(
        &mut written,
        &PairwiseDatabase::new(vec![DatabaseInfo {
            definition: b"Swiss-Prot".to_vec(),
            date: "Oct 12, 2026  1:01 AM".to_string(),
            num_sequences: 575_748,
            total_letters: 209_017_843,
            filt_algorithm_name: String::new(),
            filt_algorithm_options: String::new(),
        }]),
    )
    .unwrap();
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "Database: Swiss-Prot\n           575,748 sequences; 209,017,843 total letters\n\n"
    );
}

// NCBI BLAST+ 2.17.0 `blastn -db db/b5\ db/a5 -db_soft_mask 11` (probe m04 of the M4a-2 survey):
// the prolog shows the first name's algorithm after all titles, the epilog each name's own;
// `, options:` only for an algorithm with options (probe m06c, `'windowmasker'`).
#[test]
fn masked_databases_show_their_algorithms() {
    let database = |definition: &[u8], name: &str, options: &str| DatabaseInfo {
        definition: definition.to_vec(),
        date: "Nov 14, 2026  4:29 AM".to_string(),
        num_sequences: 5,
        total_letters: 7360,
        filt_algorithm_name: name.to_string(),
        filt_algorithm_options: options.to_string(),
    };
    let shown = PairwiseDatabase::new(vec![
        database(b"dust40 only", "dust", "window=32; level=40; linker=1"),
        database(b"mask probe v5", "windowmasker", ""),
    ]);
    let mut written = Vec::new();
    write_blast_database_header(&mut written, &shown).unwrap();
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "Database: dust40 only; mask probe v5\n\
             Masked using: 'dust', options: 'window=32; level=40; linker=1'\n\
             \x20          10 sequences; 14,720 total letters\n\n"
    );
    let mut written = Vec::new();
    write_final_database_block(&mut written, &shown).unwrap();
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "  Database: dust40 only\n\
             \x20 Masked using: 'dust', options: 'window=32; level=40; linker=1'\n\
             \x20   Posted date:  Nov 14, 2026  4:29 AM\n\
             \x20 Number of letters in database: 7,360\n\
             \x20 Number of sequences in database:  5\n\n\
             \x20 Database: mask probe v5\n\
             \x20 Masked using: 'windowmasker'\n\
             \x20   Posted date:  Nov 14, 2026  4:29 AM\n\
             \x20 Number of letters in database: 7,360\n\
             \x20 Number of sequences in database:  5\n\n\n\n"
    );
}

// The same report's epilog: the posted date of the database after its title.
#[test]
fn database_footer_has_the_posted_date() {
    let mut written = Vec::new();
    write_final_database_block(&mut written, &long_titled_database()).unwrap();
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "  Database: a database title that is much longer than sixty-eight letters, with\n\
             commas, hyphen-ated words and\n\
             a-very-long-unbroken-token-\n\
             zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz\n\
             zzzzzzzzzzzz at its end\n    Posted date:  Oct 13, 2026  1:18 PM\n  \
             Number of letters in database: 2,287\n  \
             Number of sequences in database:  17\n\n\n\n"
    );
}

fn heading_of(title: &[u8]) -> String {
    let mut hit = PairwiseHit::from(Hit {
        identity: 100.0,
        length: 10,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 10,
        s_start: 1,
        s_end: 10,
        e_value: 1e-5,
        bit_score: 30.0,
        num_ident: 10,
        query_frame: 0,
        query_length: 10,
        q_idx: 0,
        s_idx: 0,
        raw_score: 60,
        sort_query_offset: 0,
        sort_query_end: 10,
        sort_subject_offset: 0,
        sort_subject_end: 10,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 10,
    });
    hit.subject_length = Some(142);
    let mut written = Vec::new();
    write_database_subject_alignments_part(
        &mut written,
        None,
        DatabaseHeading::Title(title, &crate::report::defline::UNNAMED_SOURCE),
        true,
        &[],
        None,
        0,
        &[&hit],
        &[0],
        0..1,
        None,
        &mut crate::report::query_warnings::no_places,
        |writer: &mut Vec<u8>, _, _| writer.write_all(b"<hsp>\n"),
    )
    .unwrap();
    String::from_utf8(written).unwrap()
}

// NCBI BLAST+ 2.17.0 `blastp -db titles -query hba.faa -num_threads 1` (case
// shapes.titles.0 of the same script): the heading of a database subject is `>` and
// the stored title with no space between them, wrapped after the first white space
// from every sixtieth letter; an empty title is the generated one.
#[test]
fn database_subject_heading_is_the_stored_title_after_the_mark() {
    assert_eq!(
        heading_of(b"s_plain plain copy"),
        ">s_plain plain copy\nLength=142\n\n<hsp>\n"
    );
    assert_eq!(
        heading_of(
            b"a title that is longer than sixty letters and has spaces so that the heading wraps at a space after the sixtieth"
        ),
        ">a title that is longer than sixty letters and has spaces so that \n\
             the heading wraps at a space after the sixtieth\nLength=142\n\n<hsp>\n"
    );
    assert_eq!(
        heading_of(b"a title of exactly sixty letters long 123456789 123456789 1234"),
        ">a title of exactly sixty letters long 123456789 123456789 1234\nLength=142\n\n<hsp>\n"
    );
    let unbroken = b"a_title_without_any_space_that_is_longer_than_sixty_letters_so_that_the_heading_cannot_wrap_at_a_space_anywhere";
    assert_eq!(
        heading_of(unbroken),
        format!(
            ">{}\nLength=142\n\n<hsp>\n",
            String::from_utf8_lossy(unbroken)
        )
    );
    let seventy = [
        &[b'x'; 70][..],
        b" then words after seventy letters without a space before them at all in the title",
    ]
    .concat();
    assert_eq!(
        heading_of(&seventy),
        format!(
            ">{} \nthen words after seventy letters without a space before \n\
                 them at all in the title\nLength=142\n\n<hsp>\n",
            "x".repeat(70)
        )
    );
    // The prefixes that the generated title loses stay in the stored one.
    assert_eq!(
        heading_of(b"TPA: a third party annotation prefix"),
        ">TPA: a third party annotation prefix\nLength=142\n\n<hsp>\n"
    );
    // case shapes.punct00.0: an empty stored title.
    assert_eq!(
        heading_of(b""),
        ">unnamed protein product\nLength=142\n\n<hsp>\n"
    );
}

#[test]
fn k_bits_follows_ctoolkit_compatible() {
    assert_eq!(ncbi_k_bits_for(true), "(bits)");
    assert_eq!(ncbi_k_bits_for(false), "(Bits)");
    assert_eq!(ncbi_k_bits_for(true).len(), ncbi_k_bits_for(false).len());
}

#[test]
fn ka_value_matches_c_percent_hash_8_3g() {
    // Expected strings from C's `sprintf` `"%#8.3g "` formatting semantics
    // (Python's '%#.3g' % value gives the same strings).
    let cases: [(f64, &str); 32] = [
        (0.0, "0.00"),
        (-1.0, "-1.00"),
        (0.625, "0.625"),
        (0.41, "0.410"),
        (0.78, "0.780"),
        (1.37, "1.37"),
        (1.28, "1.28"),
        (0.46, "0.460"),
        (0.85, "0.850"),
        (0.99996, "1.00"),
        (9.9996, "10.0"),
        (999.5, "1.00e+03"),
        (1.23e-05, "1.23e-05"),
        (0.0001, "0.000100"),
        (9.999e-05, "0.000100"),
        (9.9994e-05, "0.000100"),
        (1000.0, "1.00e+03"),
        (99.95, "100."),
        (99.94, "99.9"),
        (0.125, "0.125"),
        (0.375, "0.375"),
        (12.5, "12.5"),
        (125.0, "125."),
        (1.5, "1.50"),
        (2.25e-07, "2.25e-07"),
        (123456.0, "1.23e+05"),
        (-0.00042, "-0.000420"),
        (5e-324, "4.94e-324"),
        (f64::MAX, "1.80e+308"),
        (0.0009995, "0.000999"),
        (0.00099949, "0.000999"),
        (-999.49, "-999."),
    ];
    for (value, expected) in cases {
        assert_eq!(format_ncbi_ka_value(value), expected, "value {value:e}");
        let mut field = Vec::new();
        write_ncbi_ka_field(&mut field, value).unwrap();
        assert_eq!(String::from_utf8(field).unwrap(), format!("{expected:>8} "));
    }
    assert_eq!(format_ncbi_ka_value(f64::NAN), "nan");
    assert_eq!(format_ncbi_ka_value(f64::INFINITY), "inf");
    assert_eq!(format_ncbi_ka_value(f64::NEG_INFINITY), "-inf");
}

fn make_hit() -> Hit {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    Hit {
        identity: 95.0,
        length: 100,
        mismatch: 5,
        gapopen: 0,
        q_start: 1,
        q_end: 100,
        s_start: 1,
        s_end: 100,
        e_value: 1e-50,
        bit_score: 185.5,
        num_ident: 95,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
        query_frame: 1,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
        raw_score: 200,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 95,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5088-5340
// Break-priority tiers (`EScore`): `eForced`, `ePunct`, `eComma`, `eSpace`, `eNewline`.
#[test]
fn flatfile_wrapping_preserves_ncbi_break_boundaries() {
    for (input, width, expected) in [
        ("abcd", 4, "abcd\n"),
        ("abc-def", 4, "abc-\ndef\n"),
        ("a b c", 4, "a b\nc\n"),
        ("abcd ef", 4, "abcd\nef\n"),
        ("a\n\nb", 4, "a\n\nb\n"),
        ("abc\u{8}d", 4, "abd\n"),
        ("ab\u{8}cd\nz", 68, "ab\u{8}cd\nz\n"),
        ("abc,def", 4, "abc,\ndef\n"),
    ] {
        let mut output = Vec::new();
        write_flatfile_wrapped(&mut output, input.as_bytes(), width).unwrap();
        assert_eq!(output, expected.as_bytes(), "{input:?}");
    }
}

#[test]
fn test_write_subject_header() {
    let mut output = Vec::new();
    write_subject_header(&mut output, b"seq1", Some(&b"Test sequence"[..]), Some(500)).unwrap();
    let output_str = String::from_utf8(output).unwrap();

    assert!(output_str.contains("> seq1 Test sequence"));
    assert!(output_str.contains("Length=500"));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1626 (cited at
// `write_sequence_row`): each row of `m_LineLen` columns; rpsblast's `-line_length` (1 or
// more) sets the width. A width of 1 prints one column per row with its own coordinates,
// the coordinates' column as wide as the HSP's widest coordinate.
#[test]
fn alignment_rows_of_one_column_step_one_letter() {
    let hit = Hit {
        length: 5,
        q_start: 8,
        q_end: 12,
        s_start: 98,
        s_end: 102,
        ..make_hit()
    };
    let row_of = |line_length: usize| {
        let config = PairwiseConfig {
            line_length,
            program: "rpsblast".to_string(),
            show_frame: false,
            ..PairwiseConfig::default()
        };
        let mut output = Vec::new();
        write_alignment_with_sequences(
            &mut output,
            &hit,
            "ACDeF",
            "ACDEW",
            Default::default(),
            None,
            &config,
            &PlainRows,
        )
        .unwrap();
        String::from_utf8(output).unwrap()
    };
    let one = row_of(1);
    let queries: Vec<&str> = one.lines().filter(|l| l.starts_with("Query")).collect();
    assert_eq!(
        queries,
        [
            "Query  8    A  8",
            "Query  9    C  9",
            "Query  10   D  10",
            "Query  11   e  11",
            "Query  12   F  12"
        ]
    );
    let subjects: Vec<&str> = one.lines().filter(|l| l.starts_with("Sbjct")).collect();
    assert_eq!(subjects[0], "Sbjct  98   A  98");
    assert_eq!(subjects[4], "Sbjct  102  W  102");
    // The middle line reads the letters without case (the lower-case mask is a display).
    assert_eq!(one.lines().nth(13).unwrap(), "            E");
    // Widths at and above the HSP's length give the one row of the default.
    assert_eq!(row_of(5), row_of(60));
    assert_eq!(row_of(1000), row_of(60));
    // A first row of one residue at 1 shows no coordinates (`j == 0 && start == 1 && end ==
    // 1`), as NCBI's `-line_length 1` rows of a query that aligns from its first letter.
    let hit = Hit {
        length: 2,
        q_start: 1,
        q_end: 2,
        s_start: 1,
        s_end: 2,
        ..make_hit()
    };
    let config = PairwiseConfig {
        line_length: 1,
        program: "rpsblast".to_string(),
        show_frame: false,
        ..PairwiseConfig::default()
    };
    let mut output = Vec::new();
    write_alignment_with_sequences(
        &mut output,
        &hit,
        "MQ",
        "MQ",
        Default::default(),
        None,
        &config,
        &PlainRows,
    )
    .unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "Query     M  \n          M\nSbjct     M  \n\nQuery  2  Q  2\n          Q\nSbjct  2  Q  2\n\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1626 (cited at
// `write_sequence_row`): the rule `j == 0 && start == 1 && end == 1` hides the coordinates of
// a first row of one residue at 1 in BLASTX's subject and TBLASTN's query, the rows of letters
// (NCBI 2.17.0 `-line_length 1`: "Sbjct       M  ", "Query        M  "); a translated row's
// first residue spans the nucleotides 1 to 3, so it keeps them, in TBLASTX's rows too
// ("Query  1    F  3").
#[test]
fn translated_alignment_rows_of_one_column_hide_only_a_letter_row_at_1() {
    let config = PairwiseConfig {
        line_length: 1,
        ..PairwiseConfig::default()
    };
    let hit = |q_end: usize, s_end: usize, query_frame, subject_frame| {
        let mut hit = PairwiseHit::from(Hit {
            length: 2,
            q_start: 1,
            q_end,
            s_start: 1,
            s_end,
            gap_info: None,
            ..make_hit()
        });
        hit.query_seq = Some("MR".to_string());
        hit.subject_seq = Some("MR".to_string());
        hit.query_frame = query_frame;
        hit.subject_frame = subject_frame;
        hit
    };
    let text = |write: &dyn Fn(&mut Vec<u8>) -> io::Result<()>| {
        let mut output = Vec::new();
        write(&mut output).unwrap();
        String::from_utf8(output).unwrap()
    };
    let blastx = hit(6, 2, Some(1), None);
    assert_eq!(
        text(&|out| write_blastx_alignment(out, &blastx, None, &config)),
        "Query  1  M  3\n          M\nSbjct     M  \n\n\
             Query  4  R  6\n          R\nSbjct  2  R  2\n\n"
    );
    let tblastn = hit(2, 6, None, Some(1));
    assert_eq!(
        text(&|out| write_tblastn_alignment(out, &tblastn, None, &config)),
        "Query     M  \n          M\nSbjct  1  M  3\n\n\
             Query  2  R  2\n          R\nSbjct  4  R  6\n\n"
    );
    let tblastx = hit(6, 6, Some(1), Some(1));
    assert_eq!(
        text(&|out| write_tblastx_alignment(out, &tblastx, None, 6, &config)),
        "Query  1  M  3\n          M\nSbjct  1  M  3\n\n\
             Query  4  R  6\n          R\nSbjct  4  R  6\n\n"
    );
    // On a minus frame the first residue of a translated row at 1 to 3 runs from 3 to 1.
    // The query record's length is the HSP's (frame -1 ends at its last nucleotide).
    let mut minus = hit(1, 2, Some(-1), None);
    (minus.hit.q_start, minus.hit.query_length) = (6, 6);
    assert_eq!(
        text(&|out| write_blastx_alignment(out, &minus, None, &config)),
        "Query  6  M  4\n          M\nSbjct     M  \n\n\
             Query  3  R  1\n          R\nSbjct  2  R  2\n\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1599
// NCBI BLAST+ 2.17.0 (line RG T1, cases d3.*): with a wrapped screen only the end of a
// translated row moves, to the end of its first screen; the letters are the whole row's
// and the rows of letters do not change. blastx (minus frame) `Query  2707 ... 2706` at
// 1431655766, `... 2703` at 1431655767, the row's `... 1325` at 1431655765; tblastn's
// subject (plus frame) 574..575, 574..578 and 574..1431; tblastx's query 2..3, 2..6 and
// 2..598 and minus-frame rows 598..597 and 598..594.
#[test]
fn a_translated_row_ends_where_ncbi_s_wrapped_screen_ends() {
    let config = |line_length| PairwiseConfig {
        line_length,
        ..PairwiseConfig::default()
    };
    let hit = |q: (usize, usize), s: (usize, usize), query_frame, subject_frame| {
        let mut hit = PairwiseHit::from(Hit {
            length: 5,
            q_start: q.0,
            q_end: q.1,
            s_start: s.0,
            s_end: s.1,
            gap_info: None,
            ..make_hit()
        });
        hit.query_seq = Some("MRKVL".to_string());
        hit.subject_seq = Some("MRKVL".to_string());
        hit.query_frame = query_frame;
        hit.subject_frame = subject_frame;
        hit
    };
    let text = |write: &dyn Fn(&mut Vec<u8>) -> io::Result<()>| {
        let mut output = Vec::new();
        write(&mut output).unwrap();
        String::from_utf8(output).unwrap()
    };
    // The records' lengths match the frames: the BLASTX query's of frame -1 at 2707, the
    // TBLASTX subject's of frame -1 at 598 (and the TBLASTX query's frame at 2 is +2).
    let mut blastx = hit((2707, 2693), (4, 8), Some(-1), None);
    blastx.hit.query_length = 2707;
    for (line_length, end) in [(1431655765, 2693), (1431655766, 2706), (1431655767, 2703)] {
        assert_eq!(
            text(&|out| write_blastx_alignment(out, &blastx, None, &config(line_length))),
            format!("Query  2707  MRKVL  {end}\n             MRKVL\nSbjct  4     MRKVL  8\n\n"),
            "-line_length {line_length}"
        );
    }
    let tblastn = hit((1, 5), (574, 588), None, Some(1));
    for (line_length, end) in [(1431655765, 588), (1431655766, 575), (1431655767, 578)] {
        assert_eq!(
            text(&|out| write_tblastn_alignment(out, &tblastn, None, &config(line_length))),
            format!("Query  1    MRKVL  5\n            MRKVL\nSbjct  574  MRKVL  {end}\n\n"),
            "-line_length {line_length}"
        );
    }
    let mut tblastx = hit((2, 16), (598, 584), Some(2), Some(-1));
    tblastx.subject_length = Some(598);
    for (line_length, ends) in [
        (1431655765, (16, 584)),
        (1431655766, (3, 597)),
        (1431655767, (6, 594)),
        (2147483647, (16, 584)),
    ] {
        assert_eq!(
            text(&|out| write_tblastx_alignment(out, &tblastx, None, 16, &config(line_length))),
            format!(
                "Query  2    MRKVL  {}\n            MRKVL\nSbjct  598  MRKVL  {}\n\n",
                ends.0, ends.1
            ),
            "-line_length {line_length}"
        );
    }
    // Rows of letters keep their ends (BLASTP's rows).
    let blastp = hit((1, 5), (11, 15), None, None);
    let blastp_text = |line_length| {
        text(&|out| {
            write_alignment_with_sequences(
                out,
                &blastp.hit,
                "MRKVL",
                "MRKVL",
                Default::default(),
                None,
                &PairwiseConfig {
                    program: "blastp".into(),
                    ..config(line_length)
                },
                &PlainRows,
            )
        })
    };
    assert_eq!(blastp_text(1431655766), blastp_text(60));
    assert!(blastp_text(1431655766).starts_with("Query  1   MRKVL  5\n"));
}

// `write_buffered` writes every letter once and in order across its buffer's ends, `char`s of
// two to four bytes included.
#[test]
fn buffered_letters_are_written_in_order_across_buffers() {
    let chars: String = (0..1000).map(|i| ['A', 'é', '€', '𝄞'][i % 7 % 4]).collect();
    let mut out = Vec::new();
    write_buffered(&mut out, chars.chars()).unwrap();
    assert_eq!(out, chars.as_bytes());
    let bytes: Vec<u8> = (0..1000u32).map(|i| (i % 251) as u8).collect();
    let mut out = Vec::new();
    write_buffered(&mut out, bytes.iter().copied()).unwrap();
    assert_eq!(out, bytes);
}

// The parallel writer's estimate of an HSP's text (`x_pairwise_par::hit_estimate`) holds at
// every line length: the rows of a 300-column HSP whose coordinates have 4 digits (55 bytes
// besides the letters per line) and about 300 bytes of its header lines fit in it, and it is
// at most a tenth more than the rows and its 400 bytes of header lines.
#[test]
fn the_writer_s_estimate_of_an_hsp_holds_at_every_line_length() {
    let hit = Hit {
        length: 300,
        q_start: 1001,
        q_end: 1300,
        s_start: 2001,
        s_end: 2300,
        gap_info: None,
        ..make_hit()
    };
    let letters = "ACDEFGHIKLMNPQRSTVWY".repeat(15);
    for line_length in [1, 2, 3, 7, 59, 60, 61, 299, 300, 1000, 2147483647] {
        let mut rows = Vec::new();
        write_alignment_with_sequences(
            &mut rows,
            &hit,
            &letters,
            &letters,
            Default::default(),
            None,
            &PairwiseConfig {
                line_length,
                program: "blastp".into(),
                ..PairwiseConfig::default()
            },
            &PlainRows,
        )
        .unwrap();
        let estimate = crate::report::x_pairwise_par::hit_estimate(300, line_length);
        assert!(
            rows.len() + 300 <= estimate && estimate <= rows.len() * 11 / 10 + 400,
            "-line_length {line_length}: rows {} B, estimate {estimate} B",
            rows.len()
        );
    }
}

/// NCBI's `GetWholeAlnSeqString` (c++/src/objtools/alnmgr/alnvec.cpp:199-397, cited at
/// `query_anchored::aln_vec::AlnMap::aln_seq_string`), a second transcription, for one row of
/// `segments` (a 0-based Dense-seg start, -1 for a gap, and a length in columns), as
/// written, unsigned arithmetic included: the first entry of `scrn_rights`.
fn ncbi_first_screen_right(
    segments: &[(i64, u32)],
    plus: bool,
    width: u32,
    scrn_width: u32,
) -> i64 {
    let (mut aln_pos, mut scrn_pos) = (0u32, 0u32);
    let (mut lft, mut rgt) = (-1i64, -1i64);
    let mut rights = Vec::new();
    let left_seg = segments.iter().position(|seg| seg.0 >= 0).unwrap();
    let toward = if plus { -1 } else { 1 };
    for (seg, &(start, seg_len)) in segments.iter().enumerate() {
        let len = seg_len.wrapping_mul(width);
        if start >= 0 {
            let stop = start + i64::from(len) - 1;
            let first = if plus { start } else { stop };
            if lft < 0 {
                lft = first;
                if rgt < 0 {
                    rgt = lft;
                }
            }
            let nscrns = (aln_pos.wrapping_sub(scrn_pos) / scrn_width) as i32;
            for _ in 0..nscrns {
                rights.push(rgt);
                scrn_pos = scrn_pos.wrapping_add(scrn_width);
            }
            if nscrns > 0 {
                lft = first;
            }
            let nscrns = (aln_pos.wrapping_add(len).wrapping_sub(scrn_pos) / scrn_width) as i32;
            let mut curr_pos = aln_pos;
            for _ in 0..nscrns {
                let delta = if plus {
                    scrn_width.wrapping_sub(curr_pos.wrapping_sub(scrn_pos))
                } else {
                    curr_pos.wrapping_sub(scrn_pos).wrapping_sub(scrn_width)
                } as i32 as i64;
                if if plus { lft < start } else { lft > stop } {
                    lft = first + delta;
                    rgt = lft + toward;
                } else {
                    rgt = lft + toward + delta;
                    lft += delta;
                }
                if seg == left_seg && lft == rgt {
                    rgt += toward;
                }
                rights.push(rgt);
                scrn_pos = scrn_pos.wrapping_add(scrn_width);
                curr_pos = scrn_pos;
            }
            if aln_pos.wrapping_add(len) <= scrn_pos {
                lft = -1;
            }
            rgt = if plus { stop } else { start };
        }
        aln_pos = aln_pos.wrapping_add(len);
    }
    let pos_diff = aln_pos.wrapping_sub(scrn_pos);
    if pos_diff > 0 {
        let nscrns = pos_diff / scrn_width + u32::from(pos_diff % scrn_width != 0);
        for _ in 0..nscrns {
            rights.push(rgt);
        }
    }
    rights[0]
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:310-358 (cited at
// `query_anchored::aln_vec::AlnMap::aln_seq_string`): subject rows of residues and gaps, on
// both strands, of width 1 or 3, with screens narrower than the row (wrapped ones: 3 * L -
// 2^32 nucleotides, 2, 5, 8, ...) and wider, and residue runs split into segments at random
// (the other row's gaps split them): the end of the first line that `write_hsp_rows` prints
// (`query_anchored::HspAlnVec`, the owner of the rows' coordinates) is the first screen's
// right of a second transcription of NCBI's computation (`ncbi_first_screen_right`).
#[test]
fn the_first_line_ends_at_ncbi_s_first_screen() {
    use crate::report::serial::seqalign::Strand;
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut random = |bound: u64| {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        state.wrapping_mul(0x2545_f491_4f6c_dd1d) % bound
    };
    let mut rows = crate::report::query_anchored::HspAlnVec::default();
    for case in 0..3000 {
        let letters: u32 = if case % 5 == 0 { 1 } else { 3 };
        let columns = 1 + random(60) as usize;
        // Residues, with runs of gaps after the first one (an HSP's row begins with one).
        let mut row: Vec<u8> = Vec::with_capacity(columns);
        while row.len() < columns {
            if row.is_empty() || random(4) != 0 {
                row.push(b'A');
            } else {
                let gaps = 1 + random(4) as usize;
                row.resize(row.len() + gaps, b'-');
            }
        }
        row.truncate(columns);
        let plus = random(2) == 0;
        let start = 1 + random(5000) as i64;
        let start = if plus {
            start
        } else {
            start + i64::from(letters) * columns as i64
        };
        let span = u64::from(letters) * columns as u64;
        let line_length = match (letters, random(3)) {
            // A wrapped screen narrower than the row: (W + 2^32) / 3 for a W of 2 mod 3.
            (3, 0 | 1) if span > 2 => {
                let screen = 2 + 3 * random((span - 2).div_ceil(3));
                ((screen + (1u64 << 32)) / 3) as usize
            }
            (3, 2) => [1431655765, 2147483647][random(2) as usize],
            _ => columns + random(100) as usize,
        };
        // The row's segments, runs of residues split at random, under a query row of
        // residues (width 1, plus strand).
        let mut segments = Vec::new();
        let mut next = start - 1;
        let mut column = 0;
        while column < row.len() {
            let gap = row[column] == b'-';
            let mut run = 1;
            while column + run < row.len()
                && (row[column + run] == b'-') == gap
                && (gap || random(3) != 0)
            {
                run += 1;
            }
            let nucleotides = i64::from(letters) * run as i64;
            if gap {
                segments.push((-1, run as u32));
            } else if plus {
                segments.push((next, run as u32));
                next += nucleotides;
            } else {
                segments.push((next - nucleotides + 1, run as u32));
                next -= nucleotides;
            }
            column += run;
        }
        let (mut starts, mut lens, mut strands) = (Vec::new(), Vec::new(), Vec::new());
        let mut query = 0;
        for &(subject, len) in &segments {
            starts.extend([query, subject as i32]);
            lens.push(len as i32);
            strands.extend([
                Strand::Plus,
                if plus { Strand::Plus } else { Strand::Minus },
            ]);
            query += len as i32;
        }
        rows.set(&starts, &lens, &strands, [1, letters], line_length)
            .unwrap();
        let last = line_length.min(columns) as i64 - 1;
        let (_, right) = rows.line(1, 0, 0, last).unwrap().unwrap();
        let scrn_width = (line_length as u32).wrapping_mul(letters);
        assert_eq!(
            i64::from(right) + 1,
            ncbi_first_screen_right(&segments, plus, letters, scrn_width) + 1,
            "plus={plus} letters={letters} L={line_length} segments={segments:?}"
        );
    }
}

#[test]
fn test_write_hsp_info() {
    let hit = PairwiseHit::from(make_hit());
    let config = PairwiseConfig::default();

    let mut output = Vec::new();
    write_hsp_info(&mut output, &hit, &config).unwrap();
    let output_str = String::from_utf8(output).unwrap();

    assert!(output_str.contains("Score ="));
    assert!(output_str.contains("Expect ="));
    assert!(output_str.contains("Identities ="));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3586-3591 (quoted
// at `write_hsp_info_raw_score`): under matrix-only scoring the header of an HSP is its raw
// score and a blank line, then the identities as in every report.
#[test]
fn matrix_only_hsp_header_is_the_raw_score() {
    let hit = PairwiseHit::from(make_hit());
    let config = PairwiseConfig {
        program: "blastn".to_string(),
        show_frame: false,
        ..PairwiseConfig::default()
    };
    let mut raw = Vec::new();
    write_hsp_info_raw_score(&mut raw, &hit, &config).unwrap();
    let raw = String::from_utf8(raw).unwrap();
    let mut full = Vec::new();
    write_hsp_info(&mut full, &hit, &config).unwrap();
    let full = String::from_utf8(full).unwrap();
    let header = format!(" Score = {}\n\n", hit.hit.raw_score);
    assert!(raw.starts_with(&header), "{raw}");
    assert!(!raw.contains("Expect") && !raw.contains("bits"), "{raw}");
    // The rest is the full header's from its identities on.
    let identities = full.find(" Identities = ").expect("identities");
    assert_eq!(&raw[header.len()..], &full[identities..]);
}

// NCBI 2.17.0 (line RG X6 probe `f_rm_sub_0`, ws:~/losat-bench/work/rgx6/s2): the outfmt 0
// prolog of matrix-only scoring (`write_prolog` with `disable_ka_stats`, blast_format.cpp:376-395),
// and megablast's reference after the RMBlast lines for `-task megablast`.
#[test]
fn matrix_only_prolog_is_rmblastns() {
    let intro = |megablast| {
        let mut prolog = Vec::new();
        write_blastn_pairwise_prolog(
            &mut prolog,
            "2.17.0+",
            BlastnIntro {
                megablast,
                indexed_megablast: false,
                matrix_only: true,
            },
            b"User specified sequence set (Input: s1.fa)",
            1,
            100_000,
            None,
            None,
        )
        .unwrap();
        String::from_utf8(prolog).unwrap()
    };
    let blastn = intro(false);
    assert!(
        blastn.starts_with(
            "RMBLASTN 2.17.0+\n\n\nReference: Robert M. Hubley, Arian Smit\nRMBlast - RepeatMasker Search Engine\n2010 <http://www.repeatmasker.org>\n\nReference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.\n"
        ),
        "{blastn}"
    );
    assert!(
        blastn.contains(
            "25:3389-3402.\n\n\n\nDatabase: User specified sequence set (Input: s1.fa).\n           1 sequences; 100,000 total letters\n"
        ),
        "{blastn}"
    );
    let megablast = intro(true);
    assert!(
        megablast.starts_with("RMBLASTN 2.17.0+\n\n\nReference: Robert M. Hubley"),
        "{megablast}"
    );
    assert!(
        megablast.contains("2010 <http://www.repeatmasker.org>\n\nReference: Zheng Zhang"),
        "{megablast}"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4017-4018
// A BLASTN HSP of one letter has the same start and end; its strand is its frame's.
#[test]
fn blastn_strand_comes_from_the_hsp_frame() {
    let config = PairwiseConfig {
        program: "blastn".to_string(),
        show_frame: false,
        ..PairwiseConfig::default()
    };
    for (query_frame, strand) in [(1, " Strand=Plus/Plus\n"), (-1, " Strand=Plus/Minus\n")] {
        let mut hit = PairwiseHit::from(Hit {
            length: 1,
            s_start: 10,
            s_end: 10,
            query_frame,
            ..make_hit()
        });
        hit.query_seq = Some("G".to_string());
        hit.subject_seq = Some("G".to_string());
        let mut output = Vec::new();
        write_hsp_info(&mut output, &hit, &config).unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.contains(strand), "{output}");
    }
}

// Line DBID item ID4f: BLASTN's outfmt 0 makes its rows a window of queries at a time and
// writes each window's report (`algorithm::blastn::blast_engine::run`); the reports of
// consecutive query ranges, each given only its own hits, are the report of all queries,
// on one thread and on a pool.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn blastn_reports_of_consecutive_query_ranges_are_the_whole_report() {
    let config = PairwiseConfig {
        program: "blastn".to_string(),
        show_frame: false,
        ..PairwiseConfig::default()
    };
    let karlin = KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 0.0,
        beta: 0.0,
    };
    let queries: Vec<BlastnPairwiseQuery> = (0..4)
        .map(|q_idx| BlastnPairwiseQuery {
            query_name: format!("q{q_idx} query {q_idx}").into_bytes(),
            query_length: 100,
            karlin: Some((karlin, karlin)),
            effective_search_space: 10_000,
        })
        .collect();
    let report = BlastnPairwiseReport {
        intro: BlastnIntro {
            megablast: true,
            indexed_megablast: false,
            matrix_only: false,
        },
        database_name: b"subjects".to_vec(),
        database_num_sequences: 3,
        database_total_letters: 300,
        database: None,
        reward: 1,
        penalty: -2,
        gap_open: 0,
        gap_extend: 0,
        zero_gap_extension_formula: true,
        window_size: 0,
        gapped: true,
        num_descriptions: 500,
        num_alignments: 250,
        unsearched: vec![false; queries.len()],
        epilog: true,
        parsed_subjects: Vec::new(),
        subject_sources: Vec::new(),
        subject_oids: None,
    };
    // The prolog, which the search writes before the first query batch (`run`'s
    // `write_pairwise_prologs`), then the reports.
    let mut prolog = Vec::new();
    write_blastn_pairwise_prolog(
        &mut prolog,
        "2.17.0+",
        report.intro,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
        None,
        None,
    )
    .unwrap();
    // Query 0: subjects 1 (2 HSPs) and 2; query 1: no hits; queries 2 and 3: subject 1.
    let hits: Vec<PairwiseHit> = [(0, 1), (0, 1), (0, 2), (2, 1), (3, 1)]
        .into_iter()
        .map(|(q_idx, s_idx)| {
            let mut hit = PairwiseHit::from(Hit {
                q_idx,
                s_idx,
                ..make_hit()
            });
            hit.query_seq = Some("A".repeat(100));
            hit.subject_seq = Some("A".repeat(100));
            hit.gaps = Some(0);
            hit.subject_length = Some(100);
            hit
        })
        .collect();
    let titles: Vec<Arc<[u8]>> = (0..3)
        .map(|s_idx| Arc::from(format!("s{s_idx} subject {s_idx}").as_bytes()))
        .collect();
    crate::utils::threading::with_search_pool(4, "test", |pool| {
        for x_pool in [None, Some(pool)] {
            let write = |range: std::ops::Range<usize>| {
                let range_hits: Vec<PairwiseHit> = hits
                    .iter()
                    .filter(|hit| range.contains(&(hit.hit.q_idx as usize)))
                    .cloned()
                    .collect();
                let mut written = Vec::new();
                write_blastn_pairwise_report(
                    &range_hits,
                    &mut written,
                    &config,
                    &queries,
                    range,
                    &titles,
                    &report,
                    None,
                    None,
                    x_pool,
                    None,
                    None,
                    &crate::report::db_subject::HitLoads::default(),
                )
                .unwrap();
                written
            };
            let whole = write(0..queries.len());
            let text = String::from_utf8([prolog.as_slice(), &whole].concat()).unwrap();
            assert!(text.starts_with("BLASTN 2.17.0+"), "{text}");
            assert_eq!(text.matches("Query= ").count(), queries.len(), "{text}");
            assert!(text.contains("No hits found"), "{text}");
            for cuts in [vec![1], vec![2], vec![3], vec![1, 2, 3], vec![1, 3]] {
                let mut joined = prolog.clone();
                let mut start = 0;
                for end in cuts.into_iter().chain([queries.len()]) {
                    joined.extend(write(start..end));
                    start = end;
                }
                assert_eq!(String::from_utf8(joined).unwrap(), text);
            }
        }
        Ok(())
    })
    .unwrap();
}

/// The run's `-html` data with `./.ncbirc` `ncbirc` (`None`: no file) for the references and
/// no file for the id links.
fn html_links_with(ncbirc: Option<&[u8]>) -> crate::report::html_ids::HtmlLinks {
    crate::report::html_ids::HtmlLinks {
        registry: crate::report::formatter_registry::InitConfig::of("blastp", None, false).unwrap(),
        db_name: Vec::new(),
        nucleotide: false,
        ncbirc: crate::report::formatter_registry::FormatterNcbirc::of(
            "blastp",
            ncbirc.map(<[u8]>::to_vec),
        )
        .unwrap(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:371-373 (cited at
// `HTML_PREFIX`)
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2290-2292 (cited at
// `HTML_SUFFIX`)
/// The `-html` prolog and epilog end of formats 0-4 are NCBI BLAST+ 2.17.0's in LOSAT's words
/// (line RG T4 probe `bp.s.f0.h`, `T/t4/out`: its bytes before `Database:` and its last five
/// lines): the HTML head and `<PRE>`, the version line in bold, the linked references, and after
/// the epilog's last line the end of `<PRE>`, body and document. The plain forms are unchanged.
#[test]
fn the_html_prolog_and_epilog_end_are_ncbis() {
    const NCBI_PROLOG: &str = "<HTML>\n\
<HEAD><TITLE>BLAST Search Results</TITLE></HEAD>\n\
<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" VLINK=\"#660099\" ALINK=\"#660099\">\n\
<PRE>\n\
\n\
<b>BLASTP 2.17.0+</b>\n\
\n\
\n\
<b><a\n\
href=\"https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=R\n\
etrieve&list_uids=9254694&dopt=Citation\">Reference</a>:</b>\n\
Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer,\n\
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),\n\
\"Gapped BLAST and PSI-BLAST: a new generation of protein database\n\
search programs\", Nucleic Acids Res. 25:3389-3402.\n\
\n\
\n\
<b><a\n\
href=\"https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=R\n\
etrieve&list_uids=11452024&dopt=Citation\">Reference for\n\
composition-based statistics</a>:</b>\n\
Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, Sergei\n\
Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and\n\
Stephen F. Altschul (2001), \"Improving the accuracy of PSI-BLAST\n\
protein database searches with composition-based statistics and\n\
other refinements\", Nucleic Acids Res. 29:2994-3005.\n\
\n\
\n\
\n\
";
    const NCBI_EPILOG_END: &str = "Neighboring words threshold: 11\n\
Window for multiple hits: 40\n\
</PRE>\n\
</BODY>\n\
</HTML>\n\
";
    let links = html_links_with(None);
    let mut out = Vec::new();
    write_prolog(
        &mut out,
        "blastp",
        crate::wording::REPORT_VERSION,
        false,
        false,
        false,
        Some((&links, &mut Vec::new())),
    )
    .unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        crate::wording::losat_text(NCBI_PROLOG)
    );
    let mut plain = Vec::new();
    write_prolog(
        &mut plain,
        "blastp",
        crate::wording::REPORT_VERSION,
        false,
        false,
        false,
        None,
    )
    .unwrap();
    assert!(plain.starts_with(
        format!("BLASTP {}\n\n\nReference: ", crate::wording::REPORT_VERSION).as_bytes()
    ));
    let mut end = Vec::new();
    write_epilog_end(&mut end, 11.0, 40, true).unwrap();
    assert_eq!(String::from_utf8(end).unwrap(), NCBI_EPILOG_END);
    let mut end = Vec::new();
    write_epilog_end(&mut end, 11.0, 40, false).unwrap();
    assert_eq!(
        end,
        b"Neighboring words threshold: 11\nWindow for multiple hits: 40\n"
    );
    // Matrix-only scoring names `rmblastn` in the same bold version line
    // (`BlastPrintVersionInfo` with `m_IsHTML`, blast_format.cpp:376-381).
    let mut rmblastn = Vec::new();
    write_prolog(
        &mut rmblastn,
        "blastn",
        crate::wording::REPORT_VERSION,
        false,
        false,
        true,
        Some((&links, &mut Vec::new())),
    )
    .unwrap();
    let head = format!(
        "{}<b>RMBLASTN {}</b>\n\n\nReference: Robert M. Hubley",
        std::str::from_utf8(HTML_PREFIX).unwrap(),
        crate::wording::REPORT_VERSION
    );
    assert!(rmblastn.starts_with(head.as_bytes()));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:96-99 (cited at
// `formatter_registry::FormatterNcbirc`)
/// A malformed `./.ncbirc` ends the `-html` prolog after the version line and its two newlines
/// with NCBI's exception (NCBI BLAST+ 2.17.0, line RG T4 probe `bp.di.f0.rc_bad.h`: the first 159
/// bytes of `bp.s.f0.h`'s prolog, exit status 255); the plain prolog does not read the file.
#[test]
fn a_malformed_ncbirc_ends_the_html_prolog_after_the_version_line() {
    let links = html_links_with(Some(b"this line has no equals sign\n"));
    let mut out = Vec::new();
    let error = write_prolog(
        &mut out,
        "blastp",
        crate::wording::REPORT_VERSION,
        false,
        false,
        false,
        Some((&links, &mut Vec::new())),
    )
    .unwrap_err();
    assert!(crate::ncbi_exception::NcbiException::of_io_error(&error).is_some());
    let mut expected = HTML_PREFIX.to_vec();
    expected.extend_from_slice(
        format!("<b>BLASTP {}</b>\n\n\n", crate::wording::REPORT_VERSION).as_bytes(),
    );
    assert_eq!(out, expected);
    assert_eq!(
        expected.len() - crate::wording::REPORT_VERSION.len() + crate::wording::NCBI_VERSION.len(),
        159
    );
    let mut plain = Vec::new();
    write_prolog(
        &mut plain,
        "blastp",
        crate::wording::REPORT_VERSION,
        false,
        false,
        false,
        None,
    )
    .unwrap();
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2485-2573 (cited at
// `query_anchored::output_seq`): under `-html` the rows' masks (`HtmlRowMasks`, from the rows'
// builder) put each mask's columns on a line in a `<font color="#000000">` span. NCBI BLAST+
// 2.17.0's rows of `blastp -seg yes -lcase_masking -soft_masking true -outfmt 0 -html` of a query
// of 70 Q with the letters 70-79 in lower case (line RG survey T4 §1.9, `mk.bp.segadj.f0.h`): the
// SEG mask 0-69 starts at the HSP's first column and runs past the first line, so that line has
// no span; the lower-case mask that touches it is a mask of its own (NCBI's list is not merged),
// with its span back to back. Without `-html` the same rows hold the masks in lower case only.
#[test]
fn html_rows_put_each_mask_of_the_rows_list_in_a_span() {
    const TAIL: &str = "AEVCLNSHCYPLCPDTSPYSITYRIVVIAWPQTAFFYFIE\
                        ELLDLEQNTLFHRGEDEWKVQLNYADVVIRPGINGINCHAVEQREHIFTPVCFLNKNDCG\
                        YKLMAHGSIYIWIGHHDINCSVYDASQLAPYWCWVYARRVLEFSWCHLCDGGRVYYETAC\
                        DPPCGRPATYTNNETYPWRV";
    let query = format!("{}vplqwyrgpk{TAIL}", "q".repeat(70));
    let subject = format!("{}VPLQWYRGPK{TAIL}", "Q".repeat(70))
        .replacen("SITY", "SVTY", 1)
        .replacen("CHAVE", "CHAFE", 1)
        .replacen("RRVLE", "RRVVE", 1);
    assert_eq!(query.len(), 260);
    let hit = Hit {
        length: 260,
        q_start: 1,
        q_end: 260,
        s_start: 1,
        s_end: 260,
        ..make_hit()
    };
    let config = PairwiseConfig {
        program: "blastp".to_string(),
        show_frame: false,
        ..PairwiseConfig::default()
    };
    let mut masks = crate::report::query_anchored::HspMasks::default();
    masks.set_row(
        0,
        0,
        query.as_bytes(),
        |_| true,
        (0, 259),
        true,
        1,
        [(0, 69, 0), (70, 79, 0)],
    );
    let rows = |masks: Option<&crate::report::query_anchored::HspMasks>| {
        let mut output = Vec::new();
        write_alignment_with_sequences(
            &mut output,
            &hit,
            &query,
            &subject,
            Default::default(),
            masks,
            &config,
            &PlainRows,
        )
        .unwrap();
        String::from_utf8(output).unwrap()
    };
    let html = rows(Some(&masks));
    const OPEN: &str = "<font color=\"#000000\">";
    const CLOSE: &str = "</font>";
    let lines: Vec<&str> = html.lines().collect();
    assert_eq!(lines[0], format!("Query  1    {}  60", "q".repeat(60)));
    assert_eq!(
        lines[4],
        format!(
            "Query  61   {OPEN}qqqqqqqqqq{CLOSE}{OPEN}vplqwyrgpk{CLOSE}\
             AEVCLNSHCYPLCPDTSPYSITYRIVVIAWPQTAFFYFIE  120"
        )
    );
    assert_eq!(
        lines[5],
        "            QQQQQQQQQQVPLQWYRGPKAEVCLNSHCYPLCPDTSPYS+TYRIVVIAWPQTAFFYFIE"
    );
    assert_eq!(html.matches(OPEN).count(), 2);
    // The spans are all that `-html` adds.
    assert_eq!(html.replace(OPEN, "").replace(CLOSE, ""), rows(None));
    // The table of the spans gives the masks of an HSP by its index in the run's hits.
    let table = [masks.copy_rows()];
    let of = HtmlRowMasks::new(Some(&table[..]), 4);
    assert!(of.of(3).is_none() && of.of(5).is_none());
    assert_eq!(rows(of.of(4)), html);
    assert!(HtmlRowMasks::new(None, 0).of(0).is_none());
}

// Line PS item PS7c (review B3): the HSP whose values the description table and `-sorthits`'
// identity show is the first that raises the highest bit score from 0 (`GetSeqAlignSetCalcParams`):
// a first HSP of a NaN or negative bit score does not hide a later one above 0, and a subject
// without one shows 0.0 and 0.0 and an identity of 0 over 1.
#[test]
fn the_shown_hsp_raises_the_highest_bit_score_from_zero() {
    let hit = |bit_score: f64, num_ident: usize| {
        PairwiseHit::from(Hit {
            bit_score,
            num_ident,
            ..make_hit()
        })
    };
    let (nan, negative, five, five_again) =
        (hit(f64::NAN, 1), hit(-1.0, 2), hit(5.0, 3), hit(5.0, 4));
    let shown = |hits: &[&PairwiseHit]| shown_hsp(hits).map(|hit| hit.hit.num_ident);
    assert_eq!(shown(&[&nan, &five]), Some(3));
    assert_eq!(shown(&[&negative, &five, &five_again]), Some(3));
    assert_eq!(shown(&[&nan, &negative]), None);
    assert_eq!(table_scores(shown_hsp(&[&nan, &negative])), (0.0, 0.0));
    assert_eq!(table_scores(shown_hsp(&[&nan, &five])), (5.0, 1e-50));
}

// Review R2 finding 4: with `-html` and a skipped database heading, a later part of a subject
// (`LOSAT_X_PAIRWISEPAR`) that opens with class N HSPs must write their `alignment skipped` lines
// before retrying the heading and placing its HTML posts at the part's first non-class-N HSP; a
// later part holding only class N HSPs must place no heading posts at all.
#[test]
fn later_part_html_posts_wait_for_first_shown_hsp_with_class_n() {
    use crate::report::defline::SourceTaxname;
    use crate::report::html_ids::{HtmlHeading, HtmlPost, HtmlReport};
    use crate::report::query_warnings::ReportPlace;

    let bad_source = SourceTaxname::Named(Arc::from(&b"Bad \x81 name"[..]))
        .of_stored(|| b"BL_ORD_ID:0".to_vec());
    let links = Arc::new(html_links_with(None));
    let report = HtmlReport::new(&links, [crate::seqdb::Oid::of_unrestricted(0)]);

    let make = |left_of_record: bool| {
        let mut hit = PairwiseHit::from(make_hit());
        hit.subject_length = Some(100);
        hit.edge = if left_of_record {
            PairwiseEdge::QueryLeftOfRecord
        } else {
            PairwiseEdge::Normal
        };
        hit
    };
    // 0: shown (starts heading, skipped + 1 InitConfig post)
    // 1: class N (skipped, no post)
    // 2: class N (skipped, no post)
    // 3: shown (retries heading, skipped + 1 InitConfig post)
    // 4: class N (skipped, no post)
    // 5: class N (skipped, no post)
    let owned = [
        make(false),
        make(true),
        make(true),
        make(false),
        make(true),
        make(true),
    ];
    let hits: Vec<&PairwiseHit> = owned.iter().collect();
    let indices: Vec<HspIndex> = (0..hits.len()).collect();

    let run_parts = |parts: &[std::ops::Range<usize>]| -> (Vec<u8>, Vec<(usize, HtmlPost)>) {
        let mut out = Vec::new();
        let mut placed = Vec::new();
        for part in parts {
            let posts = std::cell::RefCell::new(Vec::new());
            let html = HtmlHeading {
                report: &report,
                s_idx: 0,
                rank: 1,
                posts: &posts,
            };
            let mut hook = |place: ReportPlace, writer: &mut Vec<u8>| -> io::Result<()> {
                if let ReportPlace::Html(post) = place {
                    placed.push((writer.len(), post));
                }
                Ok(())
            };
            write_database_subject_alignments_part(
                &mut out,
                None,
                DatabaseHeading::Title(&[], &bad_source),
                false,
                &[],
                Some(&[0]),
                0,
                &hits,
                &indices,
                part.clone(),
                Some(html),
                &mut hook,
                |writer, _, _| writer.write_all(b"<hsp>\n"),
            )
            .unwrap();
        }
        (out, placed)
    };

    let (serial_out, serial_posts) = run_parts(&[0..6]);
    // Part 0..1 has first_shown; part 1..4 opens with 2 class N HSPs before shown idx 3;
    // part 4..6 has ONLY class N HSPs (must emit no posts).
    let (split_out, split_posts) = run_parts(&[0..1, 1..4, 4..6]);
    assert_eq!(split_out, serial_out);
    assert_eq!(split_posts, serial_posts);
    // Exactly two InitConfig posts: before idx 0 (byte 0) and before idx 3 (after 3 lines).
    assert_eq!(split_posts.len(), 2);
    assert_eq!(split_posts[0].0, 0);
    assert!(split_posts[1].0 > 0);
    let prefix = std::str::from_utf8(&split_out[..split_posts[1].0]).unwrap();
    assert_eq!(prefix.lines().count(), 3);
}

/// The `Gaps =` of an HSP's header, the tables' counts and the library's `HspRecord` count the
/// columns inside the query row's aligned range (`PairwiseHit::query_range_lengths`, NCBI's
/// `x_GetNumGaps` and `GetAlignLengths`): a Smith-Waterman alignment's end gaps in the query's
/// row are left out, the subject row's and an X-drop script's inner gaps stay.
#[test]
fn header_gaps_leave_out_the_query_rows_end_gaps() {
    use crate::common::GapEditOp::{Del, Ins, Sub};
    use crate::report::serial::seqalign::AlignLengths;
    let counts = |length: usize, gapopen: usize, ops: Vec<crate::common::GapEditOp>| {
        let AlignLengths {
            length,
            gaps,
            gap_opens,
        } = PairwiseHit::from(Hit {
            length,
            mismatch: 5,
            gapopen,
            gap_info: Some(ops),
            ..make_hit()
        })
        .query_range_lengths();
        (length, gaps, gap_opens)
    };
    // 95 identities and 5 mismatches in each: the rest of `length` are gap columns.
    assert_eq!(counts(195, 1, vec![Del(95), Sub(100)]), (100, 0, 0));
    assert_eq!(
        counts(105, 2, vec![Sub(60), Ins(2), Sub(40), Del(3)]),
        (102, 2, 1)
    );
    assert_eq!(counts(104, 1, vec![Ins(4), Sub(100)]), (104, 4, 1));
    assert_eq!(counts(103, 1, vec![Sub(50), Del(3), Sub(50)]), (103, 3, 1));
    let plain = PairwiseHit::from(make_hit());
    assert_eq!(plain.query_range_lengths().gaps, 0);
    assert_eq!(plain.query_range_lengths().length, plain.hit.length);
}

#[test]
fn tblastn_nullb_skip_split_parts_match_serial_for_record_and_database_html() {
    use crate::report::html_ids::{HtmlPost, HtmlReport};
    use crate::report::query_warnings::ReportPlace;

    let links = Arc::new(html_links_with(None));
    let db_report = HtmlReport::new(&links, [crate::seqdb::Oid::of_unrestricted(0)]);
    let subj_report = HtmlReport::new(&links, []);

    let make = |edge: PairwiseEdge| {
        let mut hit = PairwiseHit::from(make_hit());
        hit.subject_length = Some(59);
        hit.edge = edge;
        hit
    };
    // 0: skipped (NullbSkipTranslate)
    // 1: shown (emits heading once)
    // 2: skipped (NullbSkipTranslate)
    // 3: shown (no duplicate heading)
    let owned = [
        make(PairwiseEdge::NullbSkipTranslate),
        make(PairwiseEdge::Normal),
        make(PairwiseEdge::NullbSkipTranslate),
        make(PairwiseEdge::Normal),
    ];
    let hits: Vec<&PairwiseHit> = owned.iter().collect();
    let indices: Vec<HspIndex> = (0..hits.len()).collect();

    let run_parts = |report: &HtmlReport,
                     parts: &[std::ops::Range<usize>],
                     is_db: bool|
     -> (Vec<u8>, Vec<(usize, HtmlPost)>) {
        let source = crate::report::defline::SourceTaxname::default();
        let mut out = Vec::new();
        let mut placed = Vec::new();
        for part in parts {
            let posts = std::cell::RefCell::new(Vec::new());
            let html = HtmlHeading {
                report,
                s_idx: 0,
                rank: 1,
                posts: &posts,
            };
            let mut hook = |place: ReportPlace, writer: &mut Vec<u8>| -> io::Result<()> {
                if let ReportPlace::Html(post) = place {
                    placed.push((writer.len(), post));
                }
                Ok(())
            };
            if is_db {
                write_database_subject_alignments_part(
                    &mut out,
                    None,
                    DatabaseHeading::Title(b"a0_r2", &source),
                    false,
                    &[],
                    Some(&[0]),
                    0,
                    &hits,
                    &indices,
                    part.clone(),
                    Some(html),
                    &mut hook,
                    |writer, _, _| writer.write_all(b"<hsp>\n"),
                )
                .unwrap();
            } else {
                write_subject_alignments_part(
                    &mut out,
                    None,
                    0,
                    &[],
                    b"a0_r2",
                    false,
                    &hits,
                    &indices,
                    part.clone(),
                    Some(html),
                    &mut hook,
                    |writer, _, _| writer.write_all(b"<hsp>\n"),
                )
                .unwrap();
            }
        }
        (out, placed)
    };

    // Record (-subject) heading across split parts (including -html):
    let (rec_serial, rec_posts_s) = run_parts(&subj_report, &[0..4], false);
    let (rec_split, rec_posts_p) = run_parts(&subj_report, &[0..1, 1..2, 2..4], false);
    assert_eq!(rec_split, rec_serial);
    assert_eq!(rec_posts_p, rec_posts_s);
    let rec_str = std::str::from_utf8(&rec_serial).unwrap();
    assert_eq!(rec_str.matches("alignment skipped").count(), 2);
    assert_eq!(rec_str.matches("<a name=Subject_1></a>").count(), 1);
    assert_eq!(rec_str.matches("<hsp>").count(), 2);

    // Database (-db) heading across split parts (including -html):
    let (db_serial, db_posts_s) = run_parts(&db_report, &[0..4], true);
    let (db_split, db_posts_p) = run_parts(&db_report, &[0..1, 1..2, 2..4], true);
    assert_eq!(db_split, db_serial);
    assert_eq!(db_posts_p, db_posts_s);
    let db_str = std::str::from_utf8(&db_serial).unwrap();
    assert_eq!(db_str.matches("BL_ORD_ID:0").count(), 3);
    assert_eq!(db_str.matches("alignment skipped").count(), 2);
    assert_eq!(db_str.matches("<a name=BL_ORD_ID:0></a>").count(), 1);
    assert_eq!(db_str.matches("<hsp>").count(), 2);
}

#[test]
fn tblastn_nullb_empty_subject_row_html_preprocess_and_id_texts() {
    use crate::report::x_pairwise_par::{QueryGroup, QuerySections};

    assert_eq!(std::mem::size_of::<PairwiseEdge>(), 1);
    assert_eq!(
        subject_alignment_id_texts(&[], Some(&[7]), 0),
        [b"BL_ORD_ID:7".to_vec(), b"BL_ORD_ID:7".to_vec()]
    );
    assert_eq!(
        subject_alignment_id_texts(&[], None, 0),
        [b"Subject_1".to_vec(), b"Subject_1".to_vec()]
    );

    let check = |view: PairwiseView, left: bool| -> io::Result<()> {
        let mut hit = PairwiseHit::from(make_hit());
        hit.edge = PairwiseEdge::NullbEmptySubjectRow { left };
        let queries = [BlastpPairwiseQuery {
            query_name: Vec::new(),
            query_length: 30,
            ungapped_karlin: Default::default(),
            effective_search_space: 0,
            valid: true,
            batch_skipped: false,
        }];
        let query_titles = [Arc::<[u8]>::from(&b"q"[..])];
        let subject_titles = [Arc::<[u8]>::from(&b"s"[..])];
        let report = BlastpPairwiseReport {
            version: String::new(),
            database_name: Vec::new(),
            database_num_sequences: 1,
            database_total_letters: 57,
            matrix_name: None,
            gap_open: 11,
            gap_extend: 1,
            word_threshold: 1.0,
            window_size: 0,
            gapped_karlin: Default::default(),
            gumbel: None,
            gapped: true,
            show_linked_set_size: false,
            num_descriptions: 0,
            num_alignments: 1,
        };
        let config = PairwiseConfig {
            view,
            ..Default::default()
        };
        let loads = crate::report::db_subject::HitLoads::default();
        let oids = [0u32];
        let sections = TblastnSections {
            config: &config,
            queries: &queries,
            query_titles: &query_titles,
            query_validity: &[true],
            query_batch_skipped: &[false],
            parsed_subjects: &[],
            subject_sources: &[],
            report: &report,
            database: None,
            anchored: None,
            row_masks: HtmlRowMasks::new(None, 0),
            loads: &loads,
            subject_oids: Some(&oids),
            has_nullb_denseg_edge: true,
        };
        let group = QueryGroup::for_sections(&[(0, &hit)], &subject_titles, &sections).unwrap();
        let mut out = Vec::new();
        sections.subject_densegs(&mut out, &group, 0)
    };

    // Outfmt 0 `-html` (`!merges_alignments`): throws alnmap.cpp:753 (left=false) or :732 (left=true).
    let err_right = check(PairwiseView::Pairwise, false).unwrap_err();
    let exc_right = crate::ncbi_exception::NcbiException::of_io_error(&err_right).unwrap();
    assert_eq!(exc_right.site().line, 753);
    assert!(exc_right
        .to_string()
        .contains("CAlnVec::x_GetSeqRightSeg(): Invalid Dense-seg: Row 1 contains gaps only."));

    let err_left = check(PairwiseView::Pairwise, true).unwrap_err();
    let exc_left = crate::ncbi_exception::NcbiException::of_io_error(&err_left).unwrap();
    assert_eq!(exc_left.site().line, 732);
    assert!(exc_left
        .to_string()
        .contains("CAlnVec::x_GetSeqLeftSeg(): Invalid Dense-seg: Row 1 contains gaps only."));

    // Outfmt 1-4 (`merges_alignments`): `CreateDensegFromStdseg` succeeds without reading row 1's start.
    check(
        PairwiseView::QueryAnchored {
            identity: true,
            anchored: true,
        },
        false,
    )
    .unwrap();
    check(
        PairwiseView::QueryAnchored {
            identity: true,
            anchored: true,
        },
        true,
    )
    .unwrap();
}
