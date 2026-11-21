//! Tests of the alignment reader. The expected answers are NCBI BLAST+ 2.17.0's (psiblast
//! `-in_msa`, the probes of line PS item PS6 `ps6/probe`, named after each case): the read of an
//! alignment succeeds or fails as `CAlnReader` does, whatever the message of its report.
use super::*;

fn guess(text: &[u8]) -> AlignFormat {
    guess_format(&mut PeekAheadStream::new(text))
}

fn rows(text: &[u8]) -> Result<Vec<Vec<u8>>, AlnError> {
    read_alignment(text).map(|alignment| alignment.rows)
}

/// A well-formed Clustal file of one block: the header, two blank lines, the rows, a conservation
/// line and a blank line (what ClustalW writes).
fn clustal(header: &str, rows: &[(&str, &str)]) -> Vec<u8> {
    let mut text = format!("{header}\n\n\n");
    for (id, letters) in rows {
        text.push_str(&format!("{id:<8}{letters}\n"));
    }
    text.push_str(&format!("        {}\n\n", " ".repeat(rows[0].1.len())));
    text.into_bytes()
}

const BASIC: &[u8] = b"CLUSTAL W (1.83) multiple sequence alignment\n\n\nq1      MKVLA-GT\ns1      MKVIAAGT\ns2      -KVLAAGT\n         *:** ***\n";

#[test]
fn the_format_guess_picks_the_formats_in_order() {
    assert_eq!(guess(BASIC), AlignFormat::Clustal);
    assert_eq!(guess(b"#NEXUS\nbegin data;\n"), AlignFormat::Nexus);
    assert_eq!(guess(b"#nexus\n"), AlignFormat::Nexus);
    assert_eq!(guess(b">a\nMKV-\n>b\nMKVL\n"), AlignFormat::FastaGap);
    assert_eq!(
        guess(b" 2 5\nseq1 MKVLA\nseq2 MKVLA\n"),
        AlignFormat::Phylip
    );
    assert_eq!(guess(b""), AlignFormat::Unknown);
    assert_eq!(guess(b"\x00\x01\x02binary\n"), AlignFormat::Unknown);
}

// `xSampleIsClustal`: only line 1 can be the header ("clustalw" or "clustal w", any case, after the
// white space is cut); without one, three blocks that each end with a conservation line are.
#[test]
fn clustal_headers_and_the_three_block_rule() {
    assert_eq!(
        guess(b"CLUSTAL W\n\nq1 MKV\ns1 MKV\n"),
        AlignFormat::Clustal
    );
    assert_eq!(
        guess(b"clustal w\n\nq1 MKV\ns1 MKV\n"),
        AlignFormat::Clustal
    );
    // The sample's lines are trimmed before the test: a leading space does not matter to the
    // guess (the read of such a file still fails, `r.hdr.06`).
    assert_eq!(
        guess(b" CLUSTAL W\n\nq1 MKV\ns1 MKV\n"),
        AlignFormat::Clustal
    );
    assert_ne!(
        guess(b"CLUSTAL O(1.2.4) multiple sequence alignment\n\nq1 MKV\ns1 MKV\n"),
        AlignFormat::Clustal
    );
    // "CLUSTAL" alone is guessed as Clustal but its read fails (`r.hdr.00`, `r.hdr3.00`).
    assert_eq!(
        read_alignment(&clustal("CLUSTAL", &[("q1", "MKVLA"), ("s1", "MKVLA")])).unwrap_err(),
        AlnError::NotAvailable
    );
    for header in ["CLUSTAL W", "clustalw", "clustal w", "Clustal W"] {
        assert!(read_alignment(&clustal(header, &[("q1", "MKVLA"), ("s1", "MKVLA")])).is_ok());
    }
    // Without a header: three blocks, each ending with a conservation line, then a blank line.
    let blocks = b"q1 MKV\ns1 MKV\n *:*\n\nq1 LAG\ns1 LAG\n *:*\n\nq1 TTT\ns1 TTT\n *:*\n\n";
    assert_eq!(guess(blocks), AlignFormat::Clustal);
    let two = b"q1 MKV\ns1 MKV\n *:*\n\nq1 LAG\ns1 LAG\n *:*\n\n";
    assert_ne!(guess(two), AlignFormat::Clustal);
    // A block whose last line has a letter is no conservation line.
    assert_ne!(
        guess(b"q1 MKV\ns1 MKV\n\nq1 LAG\ns1 LAG\n"),
        AlignFormat::Clustal
    );
}

#[test]
fn a_clustal_alignment_is_read_with_its_case_and_gaps() {
    let aln = read_alignment(BASIC).unwrap();
    assert_eq!(aln.rows.len(), 3);
    assert_eq!(aln.rows[0], b"MKVLA-GT");
    assert_eq!(aln.rows[2], b"-KVLAAGT");
    // The Bioseq of each row has no gap and upper-case letters.
    assert_eq!(aln.bioseqs[0].letters, b"MKVLAGT");
    assert_eq!(aln.bioseqs[2].letters, b"KVLAAGT");
    assert!(aln.bioseqs.iter().all(|bioseq| bioseq.mol == Mol::Aa));
    let lower = rows(&clustal("CLUSTAL W", &[("q1", "mkVLA"), ("s1", "MKVLA")])).unwrap();
    assert_eq!(lower[0], b"mkVLA");
}

#[test]
fn the_lines_of_a_block_may_be_numbered_and_split_over_blocks() {
    let text = b"CLUSTAL W\n\n\nq1      MKVLA 5\ns1      MKVIA 5\n         ***:*\n\nq1      AGT 8\ns1      AGT 8\n         ***\n\n";
    let aln = rows(text).unwrap();
    assert_eq!(aln[0], b"MKVLAAGT");
    assert_eq!(aln[1], b"MKVIAAGT");
}

// What NCBI's reader rejects, in the order of its scanner: every report throws, so each of
// these is one read failure ("Seq_entry is not available until after Read()").
#[test]
fn a_read_that_fails_is_not_available() {
    let failures: &[(&str, &[u8])] = &[
        ("one sequence", b"CLUSTAL W\n\nq1 MKVLA\n"),
        ("no sequence", b"CLUSTAL W\n\n"),
        ("unequal rows", b"CLUSTAL W\n\nq1 MKVLA\ns1 MKV\n"),
        ("a duplicate id", b"CLUSTAL W\n\nq1 MKVLA\nq1 MKVLA\n"),
        (
            "ids that differ by case",
            b"CLUSTAL W\n\nq1 MKVLA\nQ1 MKVLA\n",
        ),
        ("a digit in a row", b"CLUSTAL W\n\nq1 MKV1A\ns1 MKVLA\n"),
        ("a dot in a row", b"CLUSTAL W\n\nq1 MKV.A\ns1 MKVLA\n"),
        (
            "a missing row in a later block",
            b"CLUSTAL W\n\nq1 MKV\ns1 MKV\n\nq1 LAG\n",
        ),
        ("not an alignment", b"\x00\x01\x02binary\n"),
    ];
    for (name, text) in failures {
        assert_eq!(
            read_alignment(text).unwrap_err(),
            AlnError::NotAvailable,
            "{name}"
        );
    }
}

// `x_GetSequenceMolType`: a row with both T and U is the one error that is thrown from
// `GetSeqEntry`, and printed.
#[test]
fn a_nucleotide_row_with_t_and_u_is_an_invalid_mol_type() {
    let both = "ACGUACGTACGUACGTACGUACGTACGU";
    let text = clustal("CLUSTAL W", &[("q1", both), ("s1", both)]);
    assert_eq!(read_alignment(&text).unwrap_err(), AlnError::InvalidMolType);
    let dna = read_alignment(&clustal(
        "CLUSTAL W",
        &[
            ("q1", "ACGTACGTACGTACGTACGT"),
            ("s1", "ACGTACGTACGTACGTACGT"),
        ],
    ))
    .unwrap();
    assert!(dna
        .bioseqs
        .iter()
        .all(|bioseq| bioseq.mol == Mol::Nucleotide));
}

/// Id strings of every shape that the Seq-id parser tells apart: one to four pieces of FASTA
/// tags (known ones of each field count, unknown ones of two and three letters), accessions,
/// numbers and empty pieces, joined by `|` (41,370 strings).
fn id_strings() -> Vec<Vec<u8>> {
    const PIECES: [&str; 14] = [
        "gi", "sp", "lcl", "gnl", "pdb", "pat", "ref", "xx", "abc", "", "123", "P12345", "1ABC",
        "A",
    ];
    let mut strings: Vec<Vec<u8>> = Vec::new();
    let mut last: Vec<String> = vec![String::new()];
    for length in 1..=4 {
        let mut next = Vec::new();
        for head in &last {
            for piece in PIECES {
                let joined = if length == 1 {
                    piece.to_string()
                } else {
                    format!("{head}|{piece}")
                };
                strings.push(joined.clone().into_bytes());
                next.push(joined);
            }
        }
        last = next;
    }
    strings
}

// S3 of the PS6 review, precondition 1 (`NO_SEQ_ID_WARNING`): under the reader's flags the Seq-id
// parser posts no warning, whatever the id (it fails instead); with `fParse_PartialOK`, which the
// reader does not set, the same strings do post some.
#[test]
fn the_seq_id_parser_posts_no_warning_under_the_readers_flags() {
    let strings = id_strings();
    assert_eq!(
        strings.len(),
        14 + 14 * 14 + 14 * 14 * 14 + 14 * 14 * 14 * 14
    );
    let mut with_partial = 0;
    for id in &strings {
        let mut warned = false;
        let _ = parse_ids(id, ID_PARSE_FLAGS, &mut |_| warned = true);
        assert!(!warned, "{}", String::from_utf8_lossy(id));
        let mut warned = false;
        let _ = parse_ids(id, ID_PARSE_FLAGS | parse_flags::PARTIAL_OK, &mut |_| {
            warned = true
        });
        with_partial += usize::from(warned);
        assert!(
            !matches!(
                row_ids(id),
                Err(AlnError::Internal(what)) if what == NO_SEQ_ID_WARNING
            ),
            "{}",
            String::from_utf8_lossy(id)
        );
    }
    assert!(with_partial > 0);
}

// S3 of the PS6 review, precondition 2 (`AN_ID_OF_EVERY_ROW`): the parser gives a list without
// an id only for a string of white space, and the scanner gives no such id (the id is the first
// token of a line whose ends were cut).
#[test]
fn every_scanned_id_parses_to_an_id_or_fails() {
    for id in id_strings() {
        let blank = crate::ncbi_str::truncate_spaces(&id).is_empty();
        match parse_ids(&id, ID_PARSE_FLAGS, &mut |_| {}) {
            Ok(ids) => assert_eq!(
                ids.iter().flatten().next().is_none(),
                blank,
                "{}",
                String::from_utf8_lossy(&id)
            ),
            Err(_) => assert!(!blank),
        }
        assert_eq!(
            matches!(row_ids(&id), Err(AlnError::Internal(_))),
            blank,
            "{}",
            String::from_utf8_lossy(&id)
        );
    }
    // The state itself: a blank id string, which only a caller other than the scanner can give.
    for blank in [&b""[..], b" ", b"\t \x0b\x0c\r\n"] {
        assert_eq!(
            row_ids(blank),
            Err(AlnError::Internal(AN_ID_OF_EVERY_ROW)),
            "{blank:?}"
        );
    }
    // The scanner's ids: never blank, also where a row's line starts or ends with white space of
    // any kind or holds other white space than the two delimiters.
    let texts: &[&[u8]] = &[
        BASIC,
        b"CLUSTAL W\n\n  q1 MKV\n\ts1\tMKV  \n\n",
        b"CLUSTAL W\n\n\x0bq1 MKV\n\x0c\x0bs1 MKV\x0b\n\n",
        b"CLUSTAL W\n\nq1\x0b MKV\n\x0bs1\x0c MKV\n\n",
        b"CLUSTAL W\n\n\x0b MKV 3\n\x0c LKV 3\n\n",
        b"CLUSTAL W\r\n\r\nq1 MKV\r\ns1 MKV\r\n\r\n",
    ];
    let mut scanned_ids = 0;
    for text in texts {
        let Ok(scanned) = clustal::scan(&mut PeekAheadStream::new(text)) else {
            continue;
        };
        for id in &scanned.ids {
            assert!(
                !crate::ncbi_str::truncate_spaces(id).is_empty(),
                "{:?}",
                String::from_utf8_lossy(text)
            );
            assert!(!matches!(row_ids(id), Err(AlnError::Internal(_))));
            scanned_ids += 1;
        }
    }
    assert!(scanned_ids >= 9, "{scanned_ids}");
}

// The two formats that the guess tells last (review probes `g.sequin`, `g.multalin` of line PS
// item PS6: NCBI BLAST+ 2.17.0 fails to read both: its scanners find no sequence in them).
#[test]
fn sequin_and_multalin_samples_of_the_guess_are_not_read() {
    let sequin = b"\n        10        20        30        40        50        60\nBDT63543.1  --MRFIAFLAVVMLAAFIVAALFTANSNYLDGNINKNLDLTKSLQLHGSFTEDDIAANNT\nBDT63093.1  --MSFLAFLIVALLAAYVFSLLFTLNSNYLDHKINTSLDLTSSLKLQGQFTDKDRAENNR\n\n";
    assert_eq!(guess(sequin), AlignFormat::Sequin);
    assert_eq!(read_alignment(sequin).unwrap_err(), AlnError::NotAvailable);
    let multalin = b"//\n\n       1          50\nBDT63543.1  --MRFIAFLAVVMLAAFIVAALFTANSNYLDGNINKNLDLTKSLQLHGSF\nBDT63093.1  --MSFLAFLIVALLAAYVFSLLFTLNSNYLDHKINTSLDLTSSLKLQGQF\n\n";
    assert_eq!(guess(multalin), AlignFormat::MultAlign);
    assert_eq!(
        read_alignment(multalin).unwrap_err(),
        AlnError::NotAvailable
    );
}
