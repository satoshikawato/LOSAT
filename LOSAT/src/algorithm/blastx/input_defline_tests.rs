use super::*;
use std::io::Cursor;

/// The records of `data` read as a nucleotide query file, with or without `-parse_deflines`.
fn read(data: &[u8], parse_deflines: bool) -> Result<Vec<FastaRecord>> {
    let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
    stream.parse_deflines = parse_deflines;
    parse_fasta_stream_keeping(
        &mut stream,
        false,
        false,
        100,
        false,
        None,
        &mut |_| Ok(()),
        &mut |_| Ok(()),
    )
}

/// `read` with NCBI's data loaders on (`FastaStream::seq_id_role`): the titles and letters.
fn read_with_loaders(data: &[u8], parse_deflines: bool) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    let mut stream = FastaStream::new(Cursor::new(data.to_vec()));
    stream.parse_deflines = parse_deflines;
    stream.seq_id_role = Some("query");
    let records = parse_fasta_stream_keeping(
        &mut stream,
        false,
        false,
        100,
        false,
        None,
        &mut |_| Ok(()),
        &mut |_| Ok(()),
    )?;
    Ok(records
        .into_iter()
        .map(|record| (record.title, record.sequence))
        .collect())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-133
// With the data loaders on, BLASTX's reader refuses the first lines that the other programs'
// reader refuses, with that reader's text, also under `-parse_deflines` (the Seq-id is tried
// before `CFastaReader` reads); a first line that NCBI reads as FASTA, and a Seq-id after the
// first line, are read as with the data loaders off.
#[test]
fn a_first_line_ncbi_tries_as_a_seq_id_is_refused_as_the_shared_reader_refuses_it() {
    use crate::blastinput::fasta_reader::{FastaInputSource, ReadError, ReaderConfig};
    // The first lines, each before a record with LF line ends, and whole inputs: CR-only and
    // CRLF line ends, and a file of the line alone without a line end.
    let first_lines: [&[u8]; 9] = [
        b"lcl|abc",
        b"12345",
        b"AB123456",
        b"P12345",
        b"gnl|db|x",
        b"  lcl|abc  ",
        b"lcl|abc\r",
        b"\x0blcl|abc",
        b"\x0c12345",
    ];
    let inputs = first_lines
        .iter()
        .map(|first| [*first, b"\n>q1 x\nACGTACGT\n"].concat())
        .chain([
            b"lcl|abc\r>q1 x\rACGTACGT\r".to_vec(),
            b"12345\r\n>q1 x\r\nACGTACGT\r\n".to_vec(),
            b"lcl|abc".to_vec(),
        ]);
    for data in inputs {
        let config = ReaderConfig::query("BLASTX", false, true);
        let expected =
            match FastaInputSource::from_bytes(&data, config).next_sequence(&mut |_| Ok(())) {
                Err(ReadError::Unsupported(refusal)) => refusal.to_string(),
                other => panic!("{data:?}: {other:?}"),
            };
        assert!(
            expected.contains("not supported by LOSAT's BLASTX"),
            "{expected}"
        );
        for parse_deflines in [false, true] {
            let error = read_with_loaders(&data, parse_deflines).unwrap_err();
            assert_eq!(error.to_string(), expected, "{data:?}");
        }
        // With the data loaders off the line is data: the letters of a first record without a
        // defline, or `CheckDataLine`'s error for a line that is not plausible data (`12345`).
        match read(&data, false) {
            Ok(records) => assert!(records[0].title.is_empty(), "{data:?}"),
            Err(error) => assert!(error.is::<FastaParseError>(), "{data:?}: {error}"),
        }
    }
    // The same with CR-only and CRLF line ends and without a line end: the line handed back
    // keeps the reader's line ends and line numbers.
    let fasta: [&[u8]; 10] = [
        b"ACGTACGTAC\n>q1 x\nACGT\n",
        b"ACGT ACGT\n>q1 x\nACGT\n",
        b"\nlcl|abc\n>q1 x\nACGT\n",
        b";c\nlcl|abc\n>q1 x\nACGT\n",
        b">lcl|abc x\nACGT\n>q1 x\nACGT\n",
        b"ACGTACGTAC\r>q1 x\rACGT\r",
        b"ACGTACGTAC\r\n>q1 x\r\nACGT\r\n",
        b"\r\nlcl|abc\r\n>q1 x\r\nACGT\r\n",
        b"ACGTACGTAC",
        b">q1 x",
    ];
    for data in fasta {
        let without: Vec<_> = read(data, false)
            .unwrap()
            .into_iter()
            .map(|record| (record.title, record.sequence))
            .collect();
        assert_eq!(read_with_loaders(data, false).unwrap(), without, "{data:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:383-398
// Under -parse_deflines a gap line before any defline stops the read at its line, whatever
// its size or modifiers; without it, the gap line starts a record without a defline. A gap
// line inside a record is that record's letters either way.
#[test]
fn a_gap_line_before_any_defline_stops_a_parse_deflines_read() {
    for (data, line) in [
        (&b">?5\nACGT\n>r1\nACGT\n"[..], 1),
        (b">?5", 1),
        (b">?0\nACGT\n", 1),
        (b">?abc\nACGT\n", 1),
        (b">?unk100\n>?7\nACGT\n", 1),
        (b"  >?5\r\nACGT\r\n", 1),
        (b";c\n\n#x\n>?5\nACGT\n", 4),
    ] {
        let error = read(data, true).unwrap_err();
        assert_eq!(
            error
                .downcast_ref::<NoDeflineError>()
                .map(|e| (e.0, e.to_string())),
            Some((
                line,
                format!(
                    "CFastaReader: Input doesn't start with a defline or comment around line {line}"
                )
            )),
            "{}",
            String::from_utf8_lossy(data)
        );
    }
    let records = read(b">?5\nACGT\n>r1\nACGT\n", false).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence, b"NNNNNACGT");
    let records = read(b">r0\nAC\n>?5\nGT\n", true).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].sequence, b"ACNNNNNGT");
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-780
// BLASTX's reader reads a gap line as the reader of the other programs reads it
// (`fasta_reader::reader::gap_line_size`): the same letters, lower-case runs and warnings,
// whatever the size and the modifiers, before any defline, inside a record, as its last line
// and inside a lower-case run. NCBI 2.17.0 reads them so in all five programs (probe
// `ncbi-parity-20261019/SM/gapline`).
#[test]
fn a_gap_line_reads_as_the_reader_of_the_other_programs_reads_it() {
    use crate::blastinput::fasta_reader::{
        FastaInputSource, ParseErrorCode, ReadError, ReaderConfig,
    };
    let gaps: [&[u8]; 23] = [
        b">?5",
        b">? 5 ",
        b">?0",
        b">?abc",
        b">?-5",
        b">?+5",
        b">?5abc",
        b">?",
        b">?unk",
        b">?unk7",
        b">?unk 7",
        b">?99999999999",
        b">?4294967296",
        b">?5 [gap-type=short-arm]",
        b">?5 [gap-type=bogus]",
        b">?5 [linkage-evidence=paired-ends]",
        b">?5 [gap-type=within scaffold][linkage-evidence=paired-ends]",
        b">?5 [gap-type=short-arm][gap-type=centromere]",
        b">?5 [foo=bar]",
        b">?5 [gap-type",
        b">?0 [gap-type=short-arm] x",
        b">?5 [gap-type=caf\xe9]",
        b">?5 [gap-type=a\0b] [linkage-evidence=pcr]",
    ];
    // BLASTX's records of `input` (a nucleotide query file, `-lcase_masking`): the letters with
    // the lower-case runs in lower case, as the reader of the other programs gives them, and
    // the warnings as BLASTX writes them.
    let blastx = |input: &[u8]| {
        let mut warnings = Vec::new();
        let records = parse_fasta_stream_keeping(
            &mut FastaStream::new(Cursor::new(input.to_vec())),
            false,
            true,
            100,
            false,
            None,
            &mut |_| Ok(()),
            &mut |warning| crate::algorithm::blastx::native::reader_warning(warning, &mut warnings),
        )
        .unwrap();
        let letters: Vec<Vec<u8>> = records
            .iter()
            .map(|record| {
                let mut letters = record.sequence.to_ascii_uppercase();
                for &(left, right) in &record.lowercase_masks {
                    letters[left as usize..=right as usize].make_ascii_lowercase();
                }
                letters
            })
            .collect();
        (letters, warnings)
    };
    for gap in gaps {
        for input in [
            [gap, b"\nACGT\n>r1\nAC\n"].concat(),
            [b">r0\nAC\n", gap, b"\nGT\n>r1\nAC\n"].concat(),
            [b">r0\nAC\n", gap, b"\n>r1\nAC\n"].concat(),
            [b">r0\nac\n", gap, b"\ngtAC\n>r1\nac\n>?2\n"].concat(),
        ] {
            let context = String::from_utf8_lossy(&input).into_owned();
            let (letters, warnings) = blastx(&input);
            let mut source =
                FastaInputSource::from_bytes(&input, ReaderConfig::query("BLASTX", false, false));
            let (mut messages, mut expected) = (Vec::new(), Vec::new());
            while !source.end() {
                match source.next_sequence(&mut |message: &[u8]| {
                    messages.extend_from_slice(message);
                    Ok(())
                }) {
                    Ok(record) => expected.push(record.sequence),
                    Err(ReadError::Parse {
                        code: ParseErrorCode::Eof,
                        ..
                    }) => break,
                    Err(error) => panic!("{context}: {error}"),
                }
            }
            assert_eq!(letters, expected, "{context}");
            assert_eq!(warnings, messages, "{context}");
        }
    }
    // NCBI 2.17.0's warnings and letters (`blastx -query` with `>?abc` inside a record, probe
    // case `blastx.q.mid_abc.f6.sub`; `>?5 [gap-type=caf\xe9]` keeps its byte).
    let (letters, warnings) = blastx(b">r0\nAC\n>?abc\nGT\n>r1\nAC\n>?5 [gap-type=caf\xe9]\n");
    assert_eq!(letters, vec![b"ACNGT".to_vec(), b"ACNNNNN".to_vec()]);
    assert_eq!(
        warnings,
        b"CFastaReader: Bad gap size at line 3\nCFastaReader: Problem parsing gap mods at line 3\nUnknown gap-type: caf\xe9\n"
    );
}

// `xSetSeqMol` + `AssignMolType` + `CSeqportUtil::Pack` through BLASTX's reader (folded into
// `fasta_reader::reader::scan_data_line`, line RG-XT item X7-5): a protein subject with a
// nucleotide-named id (`gnl|TI|`) drops non-IUPACna letters with NCBI's nucleotide warning; a
// nucleotide query with a protein-named id (`sp|`) keeps protein letters without warning, maps
// U->T and non-IUPACna letters to N, and preserves lowercase masks.
#[test]
fn blastx_other_molecule_accession_follows_ncbis_reading() {
    let read_pd = |input: &[u8], protein: bool| {
        let mut warnings = Vec::new();
        let mut stream = FastaStream::new(Cursor::new(input.to_vec()));
        stream.parse_deflines = true;
        let records = parse_fasta_stream_keeping(
            &mut stream,
            protein,
            true,
            100,
            false,
            None,
            &mut |_| Ok(()),
            &mut |warning| crate::algorithm::blastx::native::reader_warning(warning, &mut warnings),
        )
        .unwrap();
        (records, String::from_utf8(warnings).unwrap())
    };
    let (records, warnings) = read_pd(b">gnl|TI|12347 trace\nMKLVEFu\n", true);
    assert_eq!(records[0].sequence, b"MKVU");
    assert_eq!(records[0].lowercase_masks, vec![(3, 3)]);
    assert_eq!(
        warnings,
        "FASTA-Reader: Ignoring invalid nucleotide residues at position(s): On line 2: 3, 5-6\n"
    );

    // Nucleotide queries under protein-named accessions (including ranged ids and `>?` gaps).
    let nuc_cases: [(&[u8], &[u8], &[(i32, i32)]); 2] = [
        (
            b">sp|P99995.1|NAME_W\nACGUefILPQXZ*\n",
            b"ACGTNNNNNNNNN",
            &[(4, 5)],
        ),
        (b">sp|P99995.1:1-20\nae\n>?2\nfgA\n", b"ANNNNGA", &[(0, 5)]),
    ];
    for (input, seq, masks) in nuc_cases {
        let (records, warnings) = read_pd(input, false);
        assert_eq!(records[0].sequence, seq);
        assert_eq!(records[0].lowercase_masks, masks);
        assert!(warnings.is_empty(), "{warnings}");
    }

    // Lowercase mask across a `>?` gap line on a ranged nucleotide-named protein subject.
    let (records, warnings) = read_pd(b">gnl|TI|12347:1-20 trace\naeu\n>?2\nfgcA\n", true);
    assert_eq!(records[0].sequence, b"AUXXGCA");
    assert_eq!(records[0].lowercase_masks, vec![(0, 5)]);
    assert!(
        warnings.contains("On line 2: 2") && warnings.contains("On line 4: 1"),
        "{warnings}"
    );
}
