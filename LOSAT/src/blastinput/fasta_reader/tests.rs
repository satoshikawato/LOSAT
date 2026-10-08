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

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2457-2551 (CSeq_id::Set with
// fParse_AnyRaw | fParse_ValidLocal), 632-642 (s_CheckForFastaTag), 1634-1830
// (IdentifyAccession), 96-124 (kSupportedRawDbtags); blast_fasta_input.cpp:132-151.
#[test]
fn seq_id_lines_follow_cseq_id() {
    use super::seq_id::{classify, SeqIdLine::*};
    for (line, kind) in [
        // Oracle BI net1-net6 (data loaders on): fetched, or skipped / fatal as Seq-ids.
        (&b"AB123456"[..], GuideFormat),
        (b"AB000000", GuideFormat),
        (b"AB123456.1", GuideFormat),
        (b"ab123456", GuideFormat),
        (b"gb|AB123456.1|", FastaTag),
        (b"lcl|foo", FastaTag),
        (b"P01308", GuideFormat),
        // ... and read as FASTA (net4, net5).
        (b"MKVLAAGIVGLLLAQ", Fasta),
        (b"GCTAAAGACAATTACATAACATACACGTCAGCACGAAACT", Fasta),
        // FASTA tags: the position of the bar first, then the tag in any case.
        (b"gi|abc|", FastaTag),
        (b"GB|x", FastaTag),
        (b"gnl|db|x", FastaTag),
        (b"ref|NC_000001|", FastaTag),
        (b"tr|x", FastaTag),
        (b"gb|", Fasta),
        (b"xx|abc", Fasta),
        (b"abc|x", Fasta),
        (b"AB|123456", Fasta),
        // GI, PDB, PRF and Swiss-Prot shapes (no table).
        (b"1", Gi),
        (b"1234", Gi),
        (b"2024", Gi),
        (b"0123", Fasta),
        (b"12.5", Fasta),
        (b"1234.1", Fasta),
        (b"1ABC", Pdb),
        (b"1abc_B", Pdb),
        (b"1ABC-A", Pdb),
        (b"1ABC|A", Pdb),
        (b"1ABC\x00", Pdb),
        (b"1ABCD", Fasta),
        (b"1AB", Fasta),
        (b"1 acgt acgt", Fasta),
        (b"123456A", Prf),
        (b"123456AB", Prf),
        (b"123456AB:x", Prf),
        (b"123456ABC", Fasta),
        (b"Q9XYZ1", Swissprot),
        (b"Q9XYZ1.2", Swissprot),
        (b"A0A023GPI8", Swissprot),
        // Accession-guide formats (letters + digits), rejected without the guide's rules:
        // `ZZ123456` and `N12345` are FASTA for NCBI (unreserved prefixes).
        (b"NC_000001.11", GuideFormat),
        (b"NP_000001", GuideFormat),
        (b"P12345", GuideFormat),
        (b"ZZ123456", GuideFormat),
        (b"N12345", GuideFormat),
        (b"A?B12345", GuideFormat),
        (b"ABCD01P000001", GuideFormat),
        (b"ABCD01S00000", Fasta),
        (b"ABC_01P000001", Fasta),
        (b"XP_123", Fasta),
        (b"ACGT1234", Fasta),
        (b"contig1", Fasta),
        (b"AB 123456", Fasta),
        (b"A-B12345", Fasta),
        (b"AB123456.1.2", Fasta),
        (b"AB123456.", Fasta),
        (b"AB123456.x", Fasta),
        // General IDs of the whitelist (`dbGSS` and `dbSTS` never match).
        (b"SRA:SRR000001", General),
        (b"sra:x", General),
        (b"TIGR:abc", General),
        (b"DBGSS:x", Fasta),
        (b"dbGSS:x", Fasta),
        (b"FOO:bar", Fasta),
        (b"SRA", Fasta),
        // Not tried (empty, or not alphanumeric first), and trimmed.
        (b"", Fasta),
        (b">q", Fasta),
        (b";c", Fasta),
        (b"\xef\xbb\xbfAB123456", Fasta),
        (b"-AB123456", Fasta),
        (b" \tAB123456 \x0b", GuideFormat),
        (b"ACGT ACGT", Fasta),
        (b"ACGTACGT", Fasta),
        (b"A*0", Fasta),
    ] {
        assert_eq!(classify(line), kind, "{:?}", String::from_utf8_lossy(line));
    }
}

/// The records and rejections of an input read to its end, as the oracle's batches see
/// them (each rejected Seq-id line is one step).
fn steps(bytes: &[u8], config: ReaderConfig) -> Vec<Result<(String, Vec<u8>, Vec<u8>), String>> {
    let mut source = FastaInputSource::from_bytes(bytes, config);
    let mut steps = Vec::new();
    while !source.end() {
        match source.next_sequence(&mut |_| Ok(())) {
            Ok(record) => steps.push(Ok((record.local_id, record.title, record.sequence))),
            Err(ReadError::Unsupported(error)) => steps.push(Err(format!("{error:#}"))),
            Err(other) => panic!("{other}"),
        }
    }
    steps
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:132-159
// (oracle BI net4 and net5: each Seq-id line is one record, the next line starts a new
// one; a malformatted line is FASTA and takes the lines up to the next `>`).
#[test]
fn seq_id_records_of_the_network_oracle() {
    let blastp = ReaderConfig::query("BLASTP", true, true);
    let got = steps(
        b"P01308\nAB123456\nMKVLAAGIVGLLLAQ\nMKVLAAGIVGLL\n>pq2 two\nAYMLPDDDHWIAYNNYRFWSWFPGFLIYIHGVHPSYDQCEECTKIPKQYSLTGDFISVYY\nDKHHADKQRICCGLNTWFDNFRMTWFWCASNLCAYDSDLMFDDRVNCFMNDSDQASWAYP\n",
        blastp,
    );
    assert_eq!(got.len(), 4);
    assert!(got[0]
        .as_ref()
        .unwrap_err()
        .contains("(\"P01308\") is not a defline"));
    assert!(got[0]
        .as_ref()
        .unwrap_err()
        .contains("not supported by LOSAT's BLASTP"));
    assert!(got[1].is_err());
    assert_eq!(
        got[2],
        Ok((
            "Query_1".to_string(),
            Vec::new(),
            b"MKVLAAGIVGLLLAQMKVLAAGIVGLL".to_vec()
        ))
    );
    let (id, title, sequence) = got[3].clone().unwrap();
    assert_eq!(
        (id.as_str(), &title[..], sequence.len()),
        ("Query_2", &b"pq2 two"[..], 120)
    );

    let blastn = ReaderConfig::query("BLASTN", false, true);
    let got = steps(
        b"AB123456.1\ngb|AB123456.1|\nab123456\nlcl|foo\nP01308\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACT\nTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTG\nCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATT\nTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCT\nGATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGTAACCGAATA\nATGCGTTCGCTCTATTGACT\n>q3 three\nTTCGTACCTTGGGGGTCGTTACCACTCTGTTCCCACGAGCGGCATTTCTGGATGGCCAGC\n",
        blastn,
    );
    assert_eq!(got.len(), 7);
    assert!(got[..5].iter().all(Result::is_err));
    let (id, title, sequence) = got[5].clone().unwrap();
    assert_eq!(
        (id.as_str(), title.len(), sequence.len()),
        ("Query_1", 0, 300)
    );
    assert!(sequence.starts_with(b"GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGG"));
    let (id, title, _) = got[6].clone().unwrap();
    assert_eq!((id.as_str(), &title[..]), ("Query_2", &b"q3 three"[..]));

    // A subject (net6) is rejected the same way, naming its role.
    let subject = ReaderConfig::subject("BLASTN", false, true);
    let got = steps(b"AB000000\n", subject);
    assert!(got[0]
        .as_ref()
        .unwrap_err()
        .starts_with("the first line of the subject (\"AB000000\")"));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:352-356
// (DATA_LOADERS=none: CCustomizedFastaReader tries no line as a Seq-id; oracle BI fl_*).
#[test]
fn first_lines_without_data_loaders() {
    let sequence_line = b"\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG\n";
    for (first, messages) in [
        // fl_colon, fl_lcl: data with bad residues.
        (
            &b"SRA:SRR000001"[..],
            "FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 4, 8-13\n",
        ),
        (
            b"lcl|foo",
            "FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 1, 3-7\n",
        ),
        // fl_long: 60 letters.
        (&[b'A'; 60][..], ""),
    ] {
        let (records, got, error) = read(&[first, &sequence_line[..]].concat(), false);
        assert!(error.is_none(), "{first:?}");
        assert_eq!(got, messages, "{first:?}");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].local_id, "Query_1");
    }
    // fl_dig0, fl_gi, fl_prot: CheckDataLine fails at line 1.
    for first in [&b"0123"[..], b"gi|abc|", b"P01308"] {
        let (_, _, error) = read(&[first, &sequence_line[..]].concat(), false);
        assert_eq!(
            error.unwrap(),
            "CFastaReader: Near line 1, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.",
            "{first:?}"
        );
    }
}

/// The lines of the hyphen warnings in `messages`.
fn hyphen_lines(messages: &str) -> Vec<u64> {
    messages
        .lines()
        .filter_map(|line| {
            line.strip_prefix("CFastaReader: Hyphens are invalid and will be ignored around line ")
        })
        .map(|number| number.parse().unwrap())
        .collect()
}

/// The probe of the LR inventory (`scratch_LR/in/probe_*.fa`): eight lines, five of them
/// data lines with a hyphen.
fn probe(eol: &[u8], final_eol: bool) -> Vec<u8> {
    let lines: [&[u8]; 8] = [
        b">p1 one",
        b"TTTCCTCATGCAATTCAAAACCATGTCCGT-AATGTAGGCGAAATAGTAAA",
        b"CCATTTTACGGAGGATACCAAATTC-CTCCTTATTCAGGACCTAAC",
        b"",
        b"CTGAGGTAAACCAGGTCTCTCC-GCCCCCTTATAAAAGCTGTT",
        b">p2 two",
        b"GCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGG-CAATGACGGATATATATTAAA",
        b"AAGTGTTTTAAGATACATTG-AGGCCCGTTCGTGCTCCTCGC",
    ];
    let mut text = lines.join(eol);
    if final_eol {
        text.extend_from_slice(eol);
    }
    text
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-296 (the first line decides
// the end-of-line style; `\n\r` switches the reader to the mixed style, where it is two
// breaks). Oracle LR n_blastn_probe_{lf,cr,crlf,lfcr} and n_blastn_probe_{lf,cr}_nofinal
// (the other two files without a final line end were not run).
#[test]
fn line_end_style_is_decided_by_the_first_line() {
    // The records without their messages (whose line numbers differ).
    let records = |records: Vec<FastaRecord>| {
        records
            .into_iter()
            .map(|record| (record.local_id, record.title, record.sequence))
            .collect::<Vec<_>>()
    };
    let expected = read(&probe(b"\n", true), false);
    assert_eq!(hyphen_lines(&expected.1), vec![2, 3, 5, 7, 8]);
    let expected = records(expected.0);
    for (eol, lines) in [
        (&b"\n"[..], vec![2, 3, 5, 7, 8]),
        (b"\r", vec![2, 3, 5, 7, 8]),
        (b"\r\n", vec![2, 3, 5, 7, 8]),
        (b"\n\r", vec![3, 4, 7, 11, 13]),
    ] {
        for final_eol in [true, false] {
            let got = read(&probe(eol, final_eol), false);
            assert_eq!(hyphen_lines(&got.1), lines, "{eol:?} {final_eol}");
            assert_eq!(records(got.0), expected, "{eol:?} {final_eol}");
            assert!(got.2.is_none());
        }
    }
    // A CR-style file reads a CRLF as one break, an LF-style file a CR before the LF
    // (n_blastn_w_cr_crlf, w_lf_then_crlf, w_crlf_then_lf): the same records.
    let a = b"GGCCCAGTCCAGATCCTCGGAAGTCC";
    let b = b"ATTGGGTCATAAACAAACATCATG";
    for text in [
        [&b">q1 x\r"[..], a, b"\r\n>q2 y\r", b, b"\r\n"].concat(),
        [&b">q1 x\n"[..], a, b"\r\n>q2 y\r\n", b, b"\n"].concat(),
        [&b">q1 x\r\n"[..], a, b"\n>q2 y\n", b, b"\n"].concat(),
        [&b">q1 x\r\r\n"[..], a, b"\r\r\n>q2 y\r\r\n", b, b"\r\r\n"].concat(),
        [&b">q1 x\n"[..], a, b"\n>q2 y\n", b, b"\r"].concat(),
    ] {
        let (records, messages, error) = read(&text, false);
        assert!(error.is_none() && messages.is_empty(), "{text:?}");
        assert_eq!(
            records
                .iter()
                .map(|record| (record.title.as_slice(), record.sequence.as_slice()))
                .collect::<Vec<_>>(),
            vec![(&b"q1 x"[..], &a[..]), (b"q2 y", &b[..])],
            "{text:?}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:245-267 (x_AdvanceEOLSimple: the
// tail after an embedded alternate end of line is pushed back without the line's own
// terminator, so it is glued to the next line). Oracle LR n_blastn_two_mixed, n_lost_c.
#[test]
fn a_lone_cr_glues_the_tail_to_the_next_line() {
    let second = b"TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTGGACACTCGCT";
    let text = [
        &b">q1 first\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG\r\nCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC\r>q2 second\n"[..],
        second,
        b"\n",
    ]
    .concat();
    let (records, messages, error) = read(&text, false);
    assert!(error.is_none());
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].title, b"q1 first");
    assert_eq!(records[0].sequence.len(), 120);
    assert_eq!(records[1].title, [&b"q2 second"[..], second].concat());
    assert!(records[1].sequence.is_empty());
    assert_eq!(
        messages,
        "FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?\n"
    );
    // lost_c: CR style; the LF inside line 3 pushes `ATTAC-...` back, which is then
    // glued to `TT-GG...` (one line 4, no line 5).
    let (records, messages, _) = read(
        b">f\rAC-GTACGTAGGCTAGCTAGGATCGATCG\rGG-ACACGTAGGCTAGCTAGGATCGATCG\nATTAC-CAGACGTAGGCTAGCTAGGATCGATCG\rTT-GGACGTAGGCTAGCTAGGATCGATCG\r",
        false,
    );
    assert_eq!(hyphen_lines(&messages), vec![2, 3, 4]);
    assert_eq!(
        records[0].sequence,
        b"ACGTACGTAGGCTAGCTAGGATCGATCGGGACACGTAGGCTAGCTAGGATCGATCGATTACCAGACGTAGGCTAGCTAGGATCGATCGTTGGACGTAGGCTAGCTAGGATCGATCG"
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:413-443 (a pushed-back tail
// of at most 256 bytes steps back into the current buffer, which keeps the stream's
// end-of-file state) and c++/src/util/line_reader.cpp:100-104 (AtEOF). Oracle LR
// n_lost_g55 and n_lost_g55_nl: the tail after the second embedded end of line of the
// last line is never read (no warning for its hyphen); n_lost_a, _b, _d (one embedded
// end of line): nothing is lost.
#[test]
fn the_tail_after_a_second_embedded_line_end_can_be_lost() {
    let head = b">first\r\nTC-TTGGCTCAATCCTAGGTGGGCATGTTTCCTAATGCCC\rTT-TTTAACGTGAGGGTTCGCGTTTTTATCCCACCTAGC\r\r\nACGTACGT\n-ACGTACGTACGTACGTACGTAC";
    for final_lf in [false, true] {
        let mut text = head.to_vec();
        if final_lf {
            text.push(b'\n');
        }
        let (records, messages, error) = read(&text, false);
        assert!(error.is_none());
        assert_eq!(hyphen_lines(&messages), vec![2, 3], "{final_lf}");
        assert_eq!(
            records[0].sequence,
            b"TCTTGGCTCAATCCTAGGTGGGCATGTTTCCTAATGCCCTTTTTAACGTGAGGGTTCGCGTTTTTATCCCACCTAGCACGTACGT"
        );
    }
    for text in [
        &b">f\rAC-GTACGTAGGCTAGCTAGGATCGATCG\rGG-ACACGTAGGCTAGCTAGGATCGATCG\nATTAC-CAGACGTAGGCTAGCTAGGATCGATCG"[..],
        b">f\rAC-GTACGTAGGCTAGCTAGGATCGATCG\rGG-ACACGTAGGCTAGCTAGGATCGATCG\nATTAC-CAGACGTAGGCTAGCTAGGATCGATCG\n",
        b">f\nAC-GTACGTAGGCTAGCTAGGATCGATCG\nGG-ACACGTAGGCTAGCTAGGATCGATCG\rATTAC-CAGACGTAGGCTAGCTAGGATCGATCG",
    ] {
        let (records, messages, _) = read(text, false);
        assert_eq!(hyphen_lines(&messages), vec![2, 3, 4], "{text:?}");
        assert_eq!(
            records[0].sequence,
            b"ACGTACGTAGGCTAGCTAGGATCGATCGGGACACGTAGGCTAGCTAGGATCGATCGATTACCAGACGTAGGCTAGCTAGGATCGATCG"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:155-163 (every line read counts,
// blank, white-space and comment lines too). Oracle LR n_cmt_probe, n_cmt_probe_cr and
// n_near_{lf,cr,crlf}.
#[test]
fn line_numbers_count_skipped_lines() {
    for text in [
        &b"#c1\n;c2\n  \n>p\nAC-GTACGTAGGCTAGCTAGGATCGATCG\n"[..],
        b"#c1\r;c2\r  \r>p\rAC-GTACGTAGGCTAGCTAGGATCGATCG\r",
    ] {
        let (records, messages, error) = read(text, false);
        assert!(error.is_none());
        assert_eq!(records[0].title, b"p");
        assert_eq!(hyphen_lines(&messages), vec![5], "{text:?}");
    }
    for eol in [&b"\n"[..], b"\r", b"\r\n"] {
        let text = [
            &b">q1 a"[..],
            b"ACGTACGTAGCTAGCTAGCTAGCATCGATCGACTAGC",
            b">q2 b",
            b"@@@@@@@@@@@@@@@@@@@@@@@@@@",
            b"",
        ]
        .join(eol);
        let (records, _, error) = read(&text, false);
        assert_eq!(records.len(), 1);
        assert_eq!(
            error.unwrap(),
            "CFastaReader: Near line 4, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.",
            "{eol:?}"
        );
    }
}

/// What reading an input gives, step by step: records, rejected lines and the error that
/// ends it, with the messages written on the way.
#[derive(Debug, PartialEq)]
enum Step {
    Record(FastaRecord),
    Rejected(String),
    Failed {
        code: ParseErrorCode,
        message: String,
        line: u64,
    },
}

fn outcome<R: Read>(mut source: FastaInputSource<R>) -> (Vec<Step>, Vec<u8>) {
    let mut steps = Vec::new();
    let mut messages = Vec::new();
    while !source.end() {
        match source.next_sequence(&mut |message: &[u8]| {
            messages.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => steps.push(Step::Record(record)),
            Err(ReadError::Unsupported(error)) => steps.push(Step::Rejected(format!("{error:#}"))),
            Err(ReadError::Parse {
                code,
                message,
                line,
            }) => {
                steps.push(Step::Failed {
                    code,
                    message,
                    line,
                });
                break;
            }
            Err(ReadError::Write(error)) => panic!("{error}"),
        }
    }
    (steps, messages)
}

/// splitmix64.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.below(items.len())]
    }
}

const LINE_ENDS: [&[u8]; 4] = [b"\n", b"\r\n", b"\r", b"\n\r"];

/// One generated line (without its end): deflines, gap lines, comments, blank lines,
/// Seq-id-like lines, odd bytes and data lines with the reader's edge cases.
fn generated_line(rng: &mut Rng, out: &mut Vec<u8>) {
    const DEFLINES: [&[u8]; 14] = [
        b">q1 title",
        b">",
        b">   ",
        b"> lead",
        b">q\tx",
        b">q ACGTACGTACGTACGTACGTAC",
        b">?3",
        b">?",
        b">?abc",
        b">?unk 2",
        b">?_x y",
        b">?5 [gap-type=centromere]",
        b">q\x01c",
        b">\xff\xfe t",
    ];
    const OTHER: [&[u8]; 19] = [
        b";c",
        b"#c",
        b"!c",
        b" ;x",
        b"",
        b"  ",
        b"\t",
        b"\x0b",
        b"AB123456",
        b"1234",
        b"gb|X|",
        b"lcl|foo",
        b"SRA:x",
        b"ACGT ACGT",
        b"0123",
        b"\xef\xbb\xbfACGT",
        b"\x00",
        b"AC\xc2\xa0GT",
        b"\xff",
    ];
    const RESIDUES: &[u8] = b"ACGTACGTACGTNacgtnURYE*-;  1x\t";
    match rng.below(10) {
        0 | 1 => out.extend_from_slice(*rng.pick(&DEFLINES)),
        2 | 3 => out.extend_from_slice(*rng.pick(&OTHER)),
        _ => {
            for _ in 0..rng.below(90) {
                out.push(*rng.pick(RESIDUES));
            }
        }
    }
}

/// A generated input: a few to a few thousand lines in one line-end style, with lines
/// in other styles, an embedded line end or none at the end now and then.
fn generated_input(rng: &mut Rng, large: bool) -> Vec<u8> {
    let style = *rng.pick(&LINE_ENDS);
    let lines = if large {
        100 + rng.below(400)
    } else {
        rng.below(12)
    };
    let mut text = Vec::new();
    for index in 0..lines {
        if large && index % 40 != 0 {
            // Long runs of plain data lines, to cross the 8191-byte windows.
            for _ in 0..40 + rng.below(80) {
                text.push(*rng.pick(b"ACGTACGTacgtN-"));
            }
        } else {
            generated_line(rng, &mut text);
        }
        if rng.below(if large { 60 } else { 8 }) == 0 {
            // An embedded line end of another style inside the line.
            text.extend_from_slice(*rng.pick(&LINE_ENDS));
            for _ in 0..rng.below(40) {
                text.push(*rng.pick(b"ACGT-"));
            }
        }
        let last = index + 1 == lines;
        if !(last && rng.below(2) == 0) {
            let end = if rng.below(if large { 80 } else { 6 }) == 0 {
                *rng.pick(&LINE_ENDS)
            } else {
                style
            };
            text.extend_from_slice(end);
        }
    }
    text
}

/// A file in which an LF sits inside a CR-style line `offset` bytes before the end of
/// window `window` (8191 bytes each), followed by `tail` bytes and `end`. The CRLF style
/// of line 1 switches to CR at the embedded CR of line 2, whose tail is pushed back, so
/// every later byte comes through the pushback buffers; `CPushback_Streambuf` refills
/// them with `in_avail`, the rest of the file for an `ifstream`, when a window is used
/// up. With `offset` 1 and no `end`, a tail of 2 to 256 bytes is lost from a file
/// (`file_rows_through_from_file_and_from_bytes` has the oracle's rows of this kind).
fn window_boundary_input(window: usize, offset: usize, tail: usize, end: &[u8]) -> Vec<u8> {
    let lf = 8191 * window - offset;
    let mut text = b">first\r\nAC-GT\rACGT\n".to_vec();
    while text.len() + 53 < lf - 10 {
        text.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\r");
    }
    while text.len() < lf {
        text.push(b'C');
    }
    text.push(b'\n');
    text.extend(std::iter::repeat_n(b'G', tail));
    text.extend_from_slice(end);
    text
}

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236 (the pushback
// buffer refills with `in_avail()`). The CLI reads a file; the ABI reads bytes as a file
// with the same bytes: both must give the same records, messages and errors, and so must
// the reader with its bulk paths off (one byte at a time).
#[test]
fn bytes_read_as_a_file_with_the_same_bytes() {
    let path = std::env::temp_dir().join(format!(
        "losat-fasta-reader-bytes-vs-file-{}",
        std::process::id()
    ));
    let configs = [
        ReaderConfig::query("BLASTN", false, false),
        ReaderConfig::subject("BLASTP", true, false),
        ReaderConfig::query("TBLASTX", false, true),
    ];
    let mut check = |bytes: &[u8], what: &str| {
        std::fs::write(&path, bytes).unwrap();
        for config in configs {
            let from_bytes = outcome(FastaInputSource::from_bytes(bytes, config));
            let from_file = outcome(FastaInputSource::from_file(
                std::fs::File::open(&path).unwrap(),
                config,
            ));
            assert_eq!(from_bytes, from_file, "{what} {config:?} {bytes:?}");
            // LOSAT's bulk paths (port plan R1) against one byte at a time.
            let reference = outcome(FastaInputSource::from_bytes_byte_at_a_time(bytes, config));
            assert_eq!(from_bytes, reference, "{what} {config:?} bulk {bytes:?}");
        }
    };
    for window in 1..=2 {
        for offset in 1..=4 {
            for tail in [0, 1, 2, 23, 255, 256, 257, 4100] {
                for end in [&b""[..], b"\r", b"\n", b"\r\n"] {
                    check(
                        &window_boundary_input(window, offset, tail, end),
                        &format!("window {window} offset {offset} tail {tail} end {end:?}"),
                    );
                }
            }
        }
    }
    let mut rng = Rng(0x5f0e_2026_1008);
    for case in 0..3000 {
        check(&generated_input(&mut rng, false), &format!("small {case}"));
    }
    for case in 0..150 {
        check(&generated_input(&mut rng, true), &format!("large {case}"));
    }
    std::fs::remove_file(&path).unwrap();
}

/// Reader-only speed against `bio::io::fasta` (port plan R1), median of three runs per
/// file. Run with `LOSAT_FASTA_READER_BENCH=<file>[,<file>...] cargo test --release --lib
/// fasta_reader::tests::reader_speed -- --ignored --nocapture`.
#[test]
#[ignore]
fn reader_speed_against_bio() {
    let Ok(files) = std::env::var("LOSAT_FASTA_READER_BENCH") else {
        return;
    };
    let median = |mut times: Vec<f64>| {
        times.sort_by(f64::total_cmp);
        times[times.len() / 2]
    };
    for path in files.split(',') {
        let (mut reader_times, mut bio_times) = (Vec::new(), Vec::new());
        for _ in 0..3 {
            let start = std::time::Instant::now();
            let mut source = FastaInputSource::from_file(
                std::fs::File::open(path).unwrap(),
                ReaderConfig::subject("BLASTN", false, false),
            );
            let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
            reader_times.push(start.elapsed().as_secs_f64());
            let reader_letters: usize = records.iter().map(|record| record.sequence.len()).sum();
            let reader_records = records.len();
            drop(records);
            let start = std::time::Instant::now();
            let records: Vec<bio::io::fasta::Record> = bio::io::fasta::Reader::from_file(path)
                .unwrap()
                .records()
                .map(Result::unwrap)
                .collect();
            bio_times.push(start.elapsed().as_secs_f64());
            let bio_letters: usize = records.iter().map(|record| record.seq().len()).sum();
            assert_eq!(
                (reader_records, reader_letters),
                (records.len(), bio_letters)
            );
        }
        let (reader, bio) = (median(reader_times), median(bio_times));
        println!(
            "reader_speed\t{path}\treader {reader:.3} s\tbio {bio:.3} s\tratio {:.2}",
            reader / bio
        );
    }
}

/// The exit code and message of a program error (`crate::cli::NativeError`).
fn native(error: &anyhow::Error) -> (i32, String) {
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .unwrap_or_else(|| panic!("not a NativeError: {error:#}"));
    (native.exit, native.message.clone())
}

const NUC_TITLE_WARNING: &str = "FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?\n";

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:198-219 and
// blast_fasta_input.cpp:450-453 (GetAllSeqs reads a record whole, its messages written as
// it is read, then x_FastaToSeqLoc checks the range; the range's CInputException is not
// caught, so no later record is read).
#[test]
fn read_subjects_writes_messages_per_record_and_stops_at_a_range_past_the_end() {
    let input = b">s1 ACGTACGTACGTACGTACGTA\nACGTACGTAC-GT\n>s2 GGGGGGGGGGGGGGGGGGGGG\nACG\n>s3 TTTTTTTTTTTTTTTTTTTTTT\nACGTACGT\n";
    let config = ReaderConfig::subject("BLASTN", false, false);
    let mut written = Vec::new();
    let mut source = FastaInputSource::from_bytes(input, config);
    let range = crate::blastinput::seq_range::SequenceRange { from: 5, to: 9 };
    let error = read_subjects(&mut source, Some(&range), &mut |message: &[u8]| {
        written.extend_from_slice(message);
        Ok(())
    })
    .unwrap_err();
    assert_eq!(
        native(&error),
        (
            1,
            "BLAST query/options error: Invalid from coordinate (greater than sequence length)\nPlease refer to the BLAST+ user manual.\n".to_string()
        )
    );
    // s1's hyphen and title warnings and s2's title warning; s3 is never read.
    assert_eq!(
        String::from_utf8(written).unwrap(),
        format!(
            "CFastaReader: Hyphens are invalid and will be ignored around line 2\n{NUC_TITLE_WARNING}{NUC_TITLE_WARNING}"
        )
    );

    // A range that starts just past a record's end (from == length) is not an error, and
    // the records come back whole, each with its messages.
    let mut written = Vec::new();
    let mut source = FastaInputSource::from_bytes(input, config);
    let range = crate::blastinput::seq_range::SequenceRange { from: 3, to: 9 };
    let records = read_subjects(&mut source, Some(&range), &mut |message: &[u8]| {
        written.extend_from_slice(message);
        Ok(())
    })
    .unwrap();
    assert_eq!(
        records
            .iter()
            .map(|record| (record.local_id.as_str(), record.sequence.as_slice()))
            .collect::<Vec<_>>(),
        vec![
            ("Subject_1", &b"ACGTACGTACGT"[..]),
            ("Subject_2", b"ACG"),
            ("Subject_3", b"ACGTACGT")
        ]
    );
    assert_eq!(written.iter().filter(|&&byte| byte == b'\n').count(), 4);
    assert_eq!(
        records[0].warnings,
        format!("CFastaReader: Hyphens are invalid and will be ignored around line 2\n{NUC_TITLE_WARNING}").into_bytes()
    );
    let without_range = read_subjects(
        &mut FastaInputSource::from_bytes(input, config),
        None,
        &mut |_| Ok(()),
    )
    .unwrap();
    assert_eq!(without_range, records);
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184 (a reader
// exception while the subjects are read is `BLAST query error: <msg>`, exit 1; the
// messages of the records before it were written).
#[test]
fn read_subjects_reader_errors_end_the_program() {
    let mut written = Vec::new();
    let error = read_subjects(
        &mut FastaInputSource::from_bytes(
            b">s1 ACGTACGTACGTACGTACGTA\nACGT\n>s2\n12345\n>s3\nACGT\n",
            ReaderConfig::subject("TBLASTX", false, false),
        ),
        None,
        &mut |message: &[u8]| {
            written.extend_from_slice(message);
            Ok(())
        },
    )
    .unwrap_err();
    assert_eq!(
        native(&error),
        (
            1,
            "BLAST query error: CFastaReader: Near line 4, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n".to_string()
        )
    );
    assert_eq!(written, NUC_TITLE_WARNING.as_bytes());
    // A first line that NCBI may fetch as a Seq-id is LOSAT's explicit rejection (NCBI's
    // GetAllSeqs does not skip a subject it cannot fetch).
    let input = b"SRA:SRR000001\nACGT\n>s1\nACGT\n\n\n";
    let error = read_subjects(
        &mut FastaInputSource::from_bytes(input, ReaderConfig::subject("TBLASTN", false, true)),
        None,
        &mut |_| Ok(()),
    )
    .unwrap_err();
    assert!(error.downcast_ref::<crate::cli::NativeError>().is_none());
    assert!(
        format!("{error:#}").contains("not supported by LOSAT's TBLASTN"),
        "{error:#}"
    );
    // Without data loaders the same line is FASTA data (oracle BI fl_colon), and only
    // blank lines follow the last record: the reader ends there (`End()`).
    let mut written = Vec::new();
    let records = read_subjects(
        &mut FastaInputSource::from_bytes(input, ReaderConfig::subject("TBLASTN", false, false)),
        None,
        &mut |message: &[u8]| {
            written.extend_from_slice(message);
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(
        records
            .iter()
            .map(|record| (record.local_id.as_str(), record.sequence.as_slice()))
            .collect::<Vec<_>>(),
        vec![("Subject_1", &b"SRASRRACGT"[..]), ("Subject_2", b"ACGT")]
    );
    assert_eq!(
        written,
        b"FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 4, 8-13\n"
    );
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:855-873 (IsIStreamEmpty,
// non-Windows branch: a stream without a position is never empty; otherwise skipws `>>`
// either fails (empty) or the position is restored).
#[test]
fn stream_is_empty_reads_files_and_never_pipes() {
    let path = std::env::temp_dir().join(format!(
        "losat-fasta-reader-stream-is-empty-{}",
        std::process::id()
    ));
    let config = ReaderConfig::query("BLASTN", false, false);
    for (bytes, empty) in [
        (&b""[..], true),
        (b" \t\n\r\n\x0b\x0c", true),
        (b"\n\n>q1\nACGT\n", false),
        (b"  ACGT", false),
        (b"\x00", false),
    ] {
        std::fs::write(&path, bytes).unwrap();
        let mut source = FastaInputSource::from_file(std::fs::File::open(&path).unwrap(), config);
        assert_eq!(source.stream_is_empty(), empty, "{bytes:?}");
        if !empty && bytes.contains(&b'A') {
            // The position is restored: the reader sees the whole input.
            let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
            assert_eq!(records.len(), 1, "{bytes:?}");
            assert_eq!(records[0].sequence, b"ACGT", "{bytes:?}");
        }
    }
    std::fs::remove_file(&path).unwrap();
    #[cfg(unix)]
    for bytes in [&b""[..], b"  \n\n", b">q1\nACGT\n"] {
        let (reader, mut writer) = std::io::pipe().unwrap();
        std::io::Write::write_all(&mut writer, bytes).unwrap();
        drop(writer);
        let file = std::fs::File::from(std::os::fd::OwnedFd::from(reader));
        let mut source = FastaInputSource::from_file(file, config);
        assert!(!source.stream_is_empty(), "pipe {bytes:?}");
        if !bytes.is_empty() && bytes[0] == b'>' {
            let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
            assert_eq!(records[0].sequence, b"ACGT");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447, 1616-1679,
// 2038-2043 and fasta_reader_utils.cpp:466-483 (the bridge gives the reader's record for
// the inputs that today's checks accept: the title, `U` stored as `T` for nucleotides,
// the lowercase letters, the local ID and the title warning).
#[test]
fn from_bio_gives_the_readers_record_for_accepted_inputs() {
    let fifty = "A".repeat(50);
    let inputs: [(&[u8], bool); 6] = [
        (b">q1 a description\nACGTacgtNNRY\n", false),
        (b">q1\nACGUacguTT\nAC\n", false),
        (b">q1 x ACGTACGTACGTACGTACGTA\nACGT\n", false),
        (b">q1  two  spaces\nAC\n>q2 second\nGG\n", false),
        (b">p1 desc\nMKVLUuX*\n", true),
        (format!(">p1 q{fifty}\nMKV\n").leak().as_bytes(), true),
    ];
    for (bytes, protein) in inputs {
        let config = ReaderConfig::query("BLASTN", protein, false);
        let reader = read_all(
            &mut FastaInputSource::from_bytes(bytes, config),
            &mut |_| Ok(()),
        )
        .unwrap();
        let bio: Vec<FastaRecord> = bio::io::fasta::Reader::new(bytes)
            .records()
            .enumerate()
            .map(|(index, record)| {
                FastaRecord::from_bio(&record.unwrap(), index + 1, "Query_", protein)
            })
            .collect();
        assert_eq!(bio, reader, "{:?}", String::from_utf8_lossy(bytes));
    }
    let record = FastaRecord::from_bio(
        &bio::io::fasta::Record::with_attrs("s", None, b"ACGU"),
        7,
        "Subject_",
        false,
    );
    assert_eq!(
        (record.local_id.as_str(), record.sequence.as_slice()),
        ("Subject_7", &b"ACGT"[..])
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-466
// (a range keeps the record and its ID; the search sees the interval). The generic range,
// warning and metadata code gives the same results for a `bio` record and its bridge.
#[test]
fn cut_and_the_input_record_bridge_agree() {
    let mut record = FastaRecord::new("Query_2", b"q2 title", b"ACGTacgtAC");
    record.warnings = b"w\n".to_vec();
    let cut = record.cut(2, 6);
    assert_eq!(cut, FastaRecord::new("Query_2", b"q2 title", b"GTac"));
    assert_eq!(InputRecord::cut(&record, 2, 6), cut);

    let bio_records = vec![
        bio::io::fasta::Record::with_attrs("q1", Some("first query"), b"ACGTACGTAC"),
        bio::io::fasta::Record::with_attrs("q2", None, b"AC"),
        bio::io::fasta::Record::with_attrs("q3", Some("x"), b"acgtACGTacgtAC"),
    ];
    let records: Vec<FastaRecord> = bio_records
        .iter()
        .enumerate()
        .map(|(index, record)| FastaRecord::from_bio(record, index + 1, "Query_", false))
        .collect();
    for (bio_record, record) in bio_records.iter().zip(&records) {
        assert_eq!(bio_record.title_bytes(), record.title_bytes());
        assert_eq!(InputRecord::seq(bio_record), InputRecord::seq(record));
    }
    let range = crate::blastinput::seq_range::SequenceRange { from: 3, to: 7 };
    let from_bio = crate::blastinput::seq_range::cut_queries(&bio_records, &range);
    let from_reader = crate::blastinput::seq_range::cut_queries(&records, &range);
    assert_eq!(from_bio.input.ordinals, from_reader.input.ordinals);
    assert_eq!(from_bio.input.skipped, vec![false, true, false]);
    for (a, b) in from_bio.records.iter().zip(&from_reader.records) {
        assert_eq!(
            (InputRecord::seq(a), a.title_bytes()),
            (InputRecord::seq(b), b.title_bytes())
        );
    }
    // As subjects, q2's interval starts past its end: NCBI's range error for both.
    assert!(crate::blastinput::seq_range::cut_subjects(&bio_records, Some(&range)).is_err());
    assert!(crate::blastinput::seq_range::cut_subjects(&records, Some(&range)).is_err());
    let first_two = crate::blastinput::seq_range::SequenceRange { from: 0, to: 1 };
    let (cut_bio, _) = crate::blastinput::seq_range::cut_subjects(&bio_records, Some(&first_two))
        .unwrap()
        .unwrap();
    let (cut_reader, placements) =
        crate::blastinput::seq_range::cut_subjects(&records, Some(&first_two))
            .unwrap()
            .unwrap();
    assert_eq!(placements.length(2, 2), 14);
    for (a, b) in cut_bio.iter().zip(&cut_reader) {
        assert_eq!(InputRecord::seq(a), InputRecord::seq(b));
    }
    for index in 0..3 {
        assert_eq!(
            crate::report::query_warnings::invalid_query_warning(
                "blastn",
                index,
                &bio_records[index]
            ),
            crate::report::query_warnings::invalid_query_warning("blastn", index, &records[index])
        );
    }
    assert_eq!(
        crate::algorithm::blastn::coordination::subject_metadata_from_records(&bio_records)
            .subject_ids,
        crate::algorithm::blastn::coordination::subject_metadata_from_records(&records).subject_ids
    );
}
