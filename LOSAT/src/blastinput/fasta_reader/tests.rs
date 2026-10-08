//! Unit tests of the FASTA reader against NCBI BLAST+ 2.17.0 observations
//! (docs/evidence/losat_web_e2h/inventory: RD, LR and BI notes).

use super::*;

fn config(protein: bool) -> ReaderConfig {
    ReaderConfig::query("BLASTN", protein, false)
}

/// The records and the messages of a nucleotide (or protein) query input.
fn read(bytes: &[u8], protein: bool) -> (Vec<FastaRecord>, String, Option<String>) {
    let mut source = FastaInputSource::from_bytes(bytes, config(protein));
    let mut messages = Vec::new();
    let mut records = Vec::new();
    let mut error = None;
    while !source.end() {
        match source.next_sequence(&mut |message: &[u8]| {
            messages.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => records.push(record),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(other) => {
                error = Some(other.to_string());
                break;
            }
        }
    }
    (
        records,
        String::from_utf8_lossy(&messages).into_owned(),
        error,
    )
}

fn titles(records: &[FastaRecord]) -> Vec<&[u8]> {
    records
        .iter()
        .map(|record| record.title.as_slice())
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-225
// (ParseDefline: white space after `>` skipped, the title ends at the first byte below
// 0x20 after its first byte, trailing white space trimmed by x_ApplyMods).
#[test]
fn titles_follow_parse_defline() {
    let (records, messages, error) = read(
        b">   lead title  \nACGT\n>q1\ttab title\nACGT\n>q1 ctl\x01tail\nACGT\n>\x01q1 title\nACGT\n>q1\x00zz title\nACGT\n>q1 caf\xc3\xa9 \xff\nACGT\n>\nACGT\n>   \nACGT\n>q [organism=x] d\nACGT\n",
        false,
    );
    assert!(error.is_none() && messages.is_empty());
    assert_eq!(
        titles(&records),
        vec![
            &b"lead title"[..],
            b"q1",
            b"q1 ctl",
            b"\x01q1 title",
            b"q1",
            b"q1 caf\xc3\xa9 \xff",
            b"",
            b"",
            b"q [organism=x] d",
        ]
    );
    // CSeqIdGenerator: one ID per record, titled or not.
    assert_eq!(records[6].local_id, "Query_7");
    assert_eq!(records[6].shown_id(), b"Query_7");
    assert_eq!(records[0].shown_id(), b"lead");
    assert_eq!(records[1].local_id, "Query_2");
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1616-1679 (the title test
// runs at the end of the record, after the data-line warnings, on the untrimmed title).
#[test]
fn title_warnings_come_at_the_end_of_the_record() {
    let nucleotide_warning = "FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?\n";
    let (_, messages, _) = read(b">q1 ACGTACGTACGTACGTACGT\nAC-GT1\n", false);
    assert_eq!(
        messages,
        format!(
            "CFastaReader: Hyphens are invalid and will be ignored around line 2\nFASTA-Reader: Ignoring invalid residues at position(s): On line 2: 6\n{nucleotide_warning}"
        )
    );
    for (title, warns) in [
        (&b">q1 ACGTACGTACGTACGTACGT  \nA\n"[..], false),
        (b">ACGTACGTACGTACGTACGT\nA\n", false),
        (b">ACGTACGTACGTACGTACGTa\nA\n", true),
        (b">q acgtACGTacgtACGTacgt\nA\n", true),
        (b">q ACGUACGTACGTACGTACGT\nA\n", false),
        (b">q ACGTACGTACGTACGTACGT\r\nA\r\n", true),
    ] {
        let (_, messages, _) = read(title, false);
        assert_eq!(messages == nucleotide_warning, warns, "{title:?}");
    }
    let letters = "A".repeat(50);
    let (_, messages, _) = read(format!(">{letters}\nA\n").as_bytes(), true);
    assert!(messages.is_empty());
    let (_, messages, _) = read(format!(">q{letters}\nA\n").as_bytes(), true);
    assert_eq!(messages, "FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?\n");
    let (_, messages, _) = read(format!(">q{letters} \nA\n").as_bytes(), true);
    assert!(messages.is_empty());
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-1016 (ParseDataLine)
#[test]
fn data_lines_follow_parse_data_line() {
    let (records, messages, error) = read(
        b">q\n  AC GT\tac ; ignored XYZ\n!comment\n#c\n;c\n\nuU-*E1\n",
        false,
    );
    assert!(error.is_none());
    assert_eq!(records[0].sequence, b"ACGTactT");
    assert_eq!(
        messages,
        "CFastaReader: Hyphens are invalid and will be ignored around line 7\nFASTA-Reader: Ignoring invalid residues at position(s): On line 7: 4-6\n"
    );
    let (records, messages, _) = read(b">p\nACDEFGHIKLMNPQRSTVWYUOJBZX*\nac-1\n", true);
    assert_eq!(records[0].sequence, b"ACDEFGHIKLMNPQRSTVWYUOJBZX*ac");
    assert_eq!(
        messages,
        "CFastaReader: Hyphens are invalid and will be ignored around line 3\nFASTA-Reader: Ignoring invalid residues at position(s): On line 3: 4\n"
    );
    // U+00A0 is two bad bytes (not white space for `isspace` in the C locale); a byte
    // order mark before the sequence is three. The lines are long enough to pass
    // CheckDataLine, as in the oracle inputs scratch_RD/s_nbsp (50 letters, U+00A0,
    // 50 letters: "On line 2: 51-52") and p_bom_nodef ("On line 1: 1-3"); shorter
    // lines fail it (`check_data_line_errors`).
    let fifty = b"GTGTGAATCG".repeat(5);
    let mut text = b">q\n".to_vec();
    text.extend_from_slice(&fifty);
    text.extend_from_slice(b"\xc2\xa0");
    text.extend_from_slice(&fifty);
    text.push(b'\n');
    let (records, messages, error) = read(&text, false);
    assert!(error.is_none());
    assert_eq!(records[0].sequence, [&fifty[..], &fifty[..]].concat());
    assert_eq!(
        messages,
        "FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 51-52\n"
    );
    let mut text = b"\xef\xbb\xbf".to_vec();
    text.extend_from_slice(&fifty);
    text.push(b'\n');
    let (records, messages, error) = read(&text, false);
    assert!(error.is_none());
    assert_eq!(records[0].sequence, fifty);
    assert_eq!(records[0].title, b"");
    assert_eq!(
        messages,
        "FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 1-3\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:713-772 (CheckDataLine runs
// while the record has no residue; hyphens are neutral).
#[test]
fn check_data_line_errors() {
    for (text, line) in [
        (&b">q\n12345\n"[..], 2),
        (b">q\n----\n", 2),
        (b">q\nAC;rest\n", 2),
        (b">a\nACGT\n>b\nAC\n>c\n\n123456\n", 7),
        (b"\xef\xbb\xbf>q1\nACGT\n", 1),
        // `bad >= good / 3` with `len_to_check > 3` (fasta.cpp:751-758): 4 letters and
        // the two bytes of U+00A0; 8 letters after a byte order mark.
        (b">q\nACGT\xc2\xa0\n", 2),
        (b"\xef\xbb\xbfACGTACGT\n", 1),
    ] {
        let (_, _, error) = read(text, false);
        assert_eq!(
            error.unwrap(),
            format!("CFastaReader: Near line {line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."),
            "{text:?}"
        );
    }
    // Only before the first residue of a record.
    let (records, _, error) = read(b">q\nACGT\n12345\n", false);
    assert!(error.is_none());
    assert_eq!(records[0].sequence, b"ACGT");
    let (_, _, error) = read(b">q\nAC\n", false);
    assert!(error.is_none());
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-372,1094-1180,1385-1416
#[test]
fn gap_lines_become_runs_of_n_or_x() {
    let (records, messages, error) = read(b">q\nAC\n>?3\nGT\n>?unk 2\n>?\n>?abc\n >?2\n", false);
    assert!(error.is_none());
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].sequence, b"ACNNNGTNNNNNN");
    // `>?abc` has no digits ("Bad gap size", size 1), and the `abc` left after them is
    // not a `[key=value]` modifier, so the modifier loop warns once and gives up
    // (fasta.cpp:1118-1166; oracle scratch_RD/g_abc: both warnings for line 3).
    assert_eq!(
        messages,
        "CFastaReader: Bad gap size at line 6\nCFastaReader: Bad gap size at line 7\nCFastaReader: Problem parsing gap mods at line 7\n"
    );
    let (records, _, _) = read(b">p\nAC\n>?2\nW\n", true);
    assert_eq!(records[0].sequence, b"ACXXW");
    // `>?_x` is a defline `>x`.
    let (records, _, _) = read(b">a\nAC\n>?_x y\nGT\n", false);
    assert_eq!(titles(&records), vec![&b"a"[..], b"x y"]);
    let (_, messages, _) = read(
        b">q\nAC\n>?5 [gap-type=centromere] [linkage-evidence=map] [foo=1] [bar=2]\n>?5 [gap-type=x]\n>?5 [linkage-evidence=pcr;bogus]\n>?5 [gap-type=within-clone\n>?5 [gap-type=between-scaffolds]\n>?5 [gap-type=unknown][gap-type=telomere]\nAC\n",
        false,
    );
    assert_eq!(
        messages,
        "Unknown gap modifier name(s): bar\nUnknown gap modifier name(s): foo\nUnknown gap-type: x\nUnknown linkage-evidence: pcr;bogus\nFASTA-Reader: Unknown gap-type can have linkage-evidence of type 'unspecified' only.\nCFastaReader: Problem parsing gap mods at line 6\nCFastaReader: This gap-type should have at least one specified linkage-evidence.\nThere were conflicting gap-types around line 8\n"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1062-1075 (masks stay
// open over white space, bad bytes, comments, line ends and gap lines).
#[test]
fn lowercase_masks_follow_open_and_close_mask() {
    let (records, _, _) = read(b">q\nACgt a;x\nc1-c\n>?3\nacGT\n", false);
    assert_eq!(records[0].sequence, b"ACgtaccnnnacGT");
    let (records, _, _) = read(b">q\nAC\n>?2\nacGT\n", false);
    assert_eq!(records[0].sequence, b"ACNNacGT");
    let (records, _, _) = read(b">p\nac*dE\n", true);
    assert_eq!(records[0].sequence, b"ac*dE");
    let (records, _, _) = read(b">q\nacgu\n>?2\n", false);
    assert_eq!(records[0].sequence, b"acgtnn");
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:375-397 (blank and comment
// lines, records without a defline, records without residues).
#[test]
fn records_without_defline_or_residues() {
    let (records, messages, error) = read(b"\n;c\n ACGT\n>b\n>c\n;x\n\n>d\nA\n", false);
    assert!(error.is_none() && messages.is_empty());
    assert_eq!(
        records
            .iter()
            .map(|record| (record.local_id.as_str(), record.sequence.len()))
            .collect::<Vec<_>>(),
        vec![
            ("Query_1", 4),
            ("Query_2", 0),
            ("Query_3", 0),
            ("Query_4", 1)
        ]
    );
    // ` >q1` is a data line (trimmed to `>q1`, which passes CheckDataLine: 1 good, 1
    // bad, 3 bytes). `>`, `q` (not a nucleotide letter) and `1` are all bad residues
    // (fasta.cpp:892-898,933-935; oracle scratch_RD/p_lead_ws_gt: "On line 1: 1-3").
    // Only comment lines: no record (eEOF, `Expected defline`).
    let (records, messages, error) = read(b" >q1\nACGT\n", false);
    assert!(error.is_none());
    assert_eq!(records[0].sequence, b"ACGT");
    assert_eq!(
        messages,
        "FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 1-3\n"
    );
    let mut source = FastaInputSource::from_bytes(b"; only\n\n", config(false));
    assert!(matches!(
        source.next_sequence(&mut |_| Ok(())),
        Err(ReadError::Parse {
            code: ParseErrorCode::Eof,
            ..
        })
    ));
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:155-296 (line ends LF, CRLF, CR
// and mixed give the same lines and line numbers).
#[test]
fn line_ends_give_the_same_records() {
    let text = b">q1 a\nAC-GT\nAC\n>q2\nGG\n";
    let expected = read(text, false);
    for other in [
        text.iter()
            .flat_map(|&byte| {
                if byte == b'\n' {
                    vec![b'\r', b'\n']
                } else {
                    vec![byte]
                }
            })
            .collect::<Vec<u8>>(),
        text.iter()
            .map(|&byte| if byte == b'\n' { b'\r' } else { byte })
            .collect(),
    ] {
        let got = read(&other, false);
        assert_eq!(got.0, expected.0, "{other:?}");
        assert_eq!(got.1, expected.1, "{other:?}");
    }
}

/// The CheckDataLine error of `line` (fasta.cpp:751-758).
fn is_check_data_line_error(result: &Result<FastaRecord, ReadError>, line: u64) -> bool {
    matches!(
        result,
        Err(ReadError::Parse {
            code: ParseErrorCode::Format,
            line: error_line,
            ..
        }) if *error_line == line
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-159
// (with data loaders the first line that each ReadOneSeq call reads is tried as a
// Seq-id; a blank, `>` or comment line is not).
#[test]
fn seq_id_lines_are_rejected_with_data_loaders() {
    let with_loaders = ReaderConfig::query("BLASTN", false, true);
    for text in [&b"AB123456\n"[..], b"gb|AB123456|\n", b"12345\n"] {
        let mut source = FastaInputSource::from_bytes(text, with_loaders);
        assert!(
            matches!(
                source.next_sequence(&mut |_| Ok(())),
                Err(ReadError::Unsupported(_))
            ),
            "{text:?}"
        );
        assert!(source.end(), "{text:?}");
    }
    // After a rejection the source is past the one line that NCBI reads as the Seq-id
    // (`*++GetLineReader()` with no `UngetLine` when the CSeq_id parse succeeds,
    // blast_fasta_input.cpp:132-144), and the local-ID counter has not moved: the next
    // call reads from the following line, as NCBI's next ReadOneSeq does after a Seq-id
    // record (oracle BI net2: `AB123456` then `>q1`; net5: five Seq-id lines, then the
    // first FASTA record is Query_1). LOSAT's callers stop at the rejection (`read_all`,
    // `read_queries`).
    let mut source =
        FastaInputSource::from_bytes(b"AB123456\n12345\n>q1\nACGT\n>\nGG\n", with_loaders);
    for _ in 0..2 {
        assert!(matches!(
            source.next_sequence(&mut |_| Ok(())),
            Err(ReadError::Unsupported(_))
        ));
    }
    let q1 = source.next_sequence(&mut |_| Ok(())).unwrap();
    assert_eq!(
        (q1.local_id.as_str(), q1.title(), q1.seq()),
        ("Query_1", &b"q1"[..], &b"ACGT"[..])
    );
    let untitled = source.next_sequence(&mut |_| Ok(())).unwrap();
    assert_eq!(
        (untitled.local_id.as_str(), untitled.title(), untitled.seq()),
        ("Query_2", &b""[..], &b"GG"[..])
    );
    assert!(source.end());
    // A line of letters only throws "Malformatted ID" (Seq_id.cpp:2496-2503 on the retry
    // without fParse_ValidLocal): it is FASTA data, and the record takes the next lines.
    let mut source = FastaInputSource::from_bytes(b"ACGTACGT\nACGT\n", with_loaders);
    assert_eq!(
        source.next_sequence(&mut |_| Ok(())).unwrap().sequence,
        b"ACGTACGTACGT"
    );
    // A blank first line is not tried (`line.empty()`): CFastaReader reads `AB123456` as
    // data, and two letters in eight bytes fail CheckDataLine.
    let mut source = FastaInputSource::from_bytes(b"\nAB123456\n", with_loaders);
    assert!(is_check_data_line_error(
        &source.next_sequence(&mut |_| Ok(())),
        2
    ));
    // Without data loaders (`DATA_LOADERS=none`: CCustomizedFastaReader,
    // blast_fasta_input.cpp:353-356) no line is tried: oracle BI fl_acc, fl_gb and fl_dig
    // fail CheckDataLine at line 1, fl_text (`ACGT ACGT`) is one record.
    let without_loaders = ReaderConfig::query("BLASTN", false, false);
    for first in [&b"AB123456"[..], b"gb|AB123456.1|", b"1234"] {
        let text = [first, b"\nGCTAAAGACAATTACATAACATACACGTCAGC\n"].concat();
        let mut source = FastaInputSource::from_bytes(&text, without_loaders);
        assert!(
            is_check_data_line_error(&source.next_sequence(&mut |_| Ok(())), 1),
            "{first:?}"
        );
    }
    let mut source = FastaInputSource::from_bytes(b"ACGT ACGT\nGCTA\n", without_loaders);
    let record = source.next_sequence(&mut |_| Ok(())).unwrap();
    assert_eq!(
        (record.local_id.as_str(), record.seq()),
        ("Query_1", &b"ACGTACGTGCTA"[..])
    );
}
