//! Tests of the mask FASTA reader against dustmasker and segmasker 2.17.0.
//!
//! The fixtures are NCBI's outputs for the inputs of `tests/dustmasker/in` (`../cases.tsv`)
//! and of the reader's own probes (`tests/dustmasker/reader`). For every case that runs a
//! program with `-in <file>` and at most `-outfmt fasta` or `-outfmt interval` and
//! `-parse_seqids` (an environment `DIAG_POST_LEVEL=Warning` lets the warnings through),
//! the reader must give as many records as NCBI wrote, NCBI's letters (`-outfmt fasta`,
//! upper-cased: the masked letters are in lower case), NCBI's ids (`-parse_seqids`, the id
//! part of the interval format's definition lines) and NCBI's standard error (warnings, and
//! the report of the error that ends the reading).

use std::path::{Path, PathBuf};

use super::*;
use crate::seqid::label;
use crate::writedb::fasta_source::reader::ReaderWarning;

/// What the reader gives for an input: the records, the warnings (each a line as NCBI
/// writes it at the post level Warning, ending with a line feed) and the error.
struct Run {
    records: Vec<MaskRecord>,
    warnings: String,
    error: Option<MaskReadError>,
}

fn run(input: &[u8], protein: bool, parse_seqids: bool) -> Run {
    let mut records = Vec::new();
    let mut warnings = String::new();
    let mut reader = match MaskFastaReader::new(input, protein, parse_seqids) {
        Ok(reader) => reader,
        Err(error) => {
            return Run {
                records,
                warnings,
                error: Some(MaskReadError::Read(error)),
            }
        }
    };
    let mut warn = |warning: ReaderWarning<'_>| match warning {
        ReaderWarning::Line(line) => {
            warnings.push_str(&String::from_utf8_lossy(line));
            warnings.push('\n');
        }
        // `ERR_POST_X` of an application without a post prefix.
        ReaderWarning::SeqId(warning) => {
            warnings.push_str(&format!(
                "Warning: (1306.{}) {}\n",
                warning.subcode,
                String::from_utf8_lossy(&warning.message)
            ));
        }
    };
    let error = loop {
        match reader.next(&mut warn) {
            Ok(Some(record)) => records.push(record),
            Ok(None) => break None,
            Err(error) => break Some(error),
        }
    };
    assert!(matches!(reader.next(&mut |_| ()), Ok(None)));
    Run {
        records,
        warnings,
        error,
    }
}

/// `CSeq_id::GetStringDescr` of the Bioseq with `eFormat_FastA`, for a record with the ids `ids`,
/// as the masking programs write it in their id lines (`label::bioseq_fasta`: its gi and best id).
fn fasta_ids(ids: &[ParsedSeqId], protein: bool) -> String {
    let mut out = Vec::new();
    label::bioseq_fasta(ids, protein, &mut out);
    String::from_utf8(out).unwrap()
}

/// One case of a `cases.tsv` that the reader decides.
struct Case {
    name: String,
    dir: PathBuf,
    segmasker: bool,
    input: String,
    fasta: bool,
    parse_seqids: bool,
    warnings: bool,
}

/// The cases of `dir/cases.tsv` that run a masking program with only the options the reader
/// decides.
fn cases(dir: &Path) -> Vec<Case> {
    let text = std::fs::read_to_string(dir.join("cases.tsv")).unwrap();
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (name, program, args) = match fields[..] {
            [name, args] => (name, "dustmasker", args),
            [name, program, args] => (name, program, args),
            _ => panic!("{line}"),
        };
        let mut words = args.split_whitespace().peekable();
        let mut warnings = false;
        while let Some(word) = words.peek() {
            if *word == "DIAG_POST_LEVEL=Warning" {
                warnings = true;
                words.next();
            } else {
                break;
            }
        }
        let (mut input, mut fasta, mut parse_seqids, mut decided) = (None, false, false, true);
        while let Some(word) = words.next() {
            match (word, words.peek().copied()) {
                ("-in", Some(file)) if file != "-" && input.is_none() => {
                    input = Some(file.to_string());
                    words.next();
                }
                ("-outfmt", Some("fasta")) => {
                    fasta = true;
                    words.next();
                }
                ("-outfmt", Some("interval")) => {
                    words.next();
                }
                ("-parse_seqids", next) if next.is_none_or(|next| next.starts_with('-')) => {
                    parse_seqids = true;
                }
                _ => decided = false,
            }
        }
        if !decided || !matches!(program, "dustmasker" | "segmasker") {
            continue;
        }
        let Some(input) = input else { continue };
        if !dir.join(&input).is_file() {
            continue;
        }
        out.push(Case {
            name: name.to_string(),
            dir: dir.to_path_buf(),
            segmasker: program == "segmasker",
            input,
            fasta,
            parse_seqids,
            warnings,
        });
    }
    out
}

/// The records of a masking program's output: each definition line (without `>`) and the
/// lines after it.
fn ncbi_records(out: &str) -> Vec<(String, String)> {
    let mut records: Vec<(String, String)> = Vec::new();
    for line in out.lines() {
        if let Some(defline) = line.strip_prefix('>') {
            records.push((defline.to_string(), String::new()));
        } else if let Some(record) = records.last_mut() {
            record.1.push_str(line);
        }
    }
    records
}

/// Compares a case with NCBI's output; returns the differences.
fn check(case: &Case) -> Vec<String> {
    let ncbi = |ext: &str| {
        std::fs::read(case.dir.join("ncbi").join(format!("{}.{ext}", case.name)))
            .map(|bytes| String::from_utf8_lossy(&bytes).into_owned())
            .unwrap()
    };
    let (out, err, rc) = (ncbi("out"), ncbi("err"), ncbi("rc"));
    let input = std::fs::read(case.dir.join(&case.input)).unwrap();
    let run = run(&input, case.segmasker, case.parse_seqids);
    let mut problems = Vec::new();
    let expected = ncbi_records(&out);
    if expected.len() != run.records.len() {
        problems.push(format!(
            "{} records, NCBI wrote {}",
            run.records.len(),
            expected.len()
        ));
    }
    for (index, (record, (defline, body))) in run.records.iter().zip(&expected).enumerate() {
        if case.fasta {
            let letters = String::from_utf8_lossy(&record.letters);
            if *letters != body.to_ascii_uppercase() {
                problems.push(format!(
                    "record {index}: letters {letters}, NCBI {}",
                    body.to_ascii_uppercase()
                ));
            }
        }
        if case.parse_seqids {
            let ncbi_ids = defline.split(' ').next().unwrap_or_default();
            let ids = fasta_ids(&record.ids, case.segmasker);
            if ids != ncbi_ids {
                problems.push(format!("record {index}: ids {ids}, NCBI {ncbi_ids}"));
            }
        }
    }
    let mut stderr = if case.warnings {
        run.warnings.clone()
    } else {
        String::new()
    };
    let mut status = 0;
    match &run.error {
        None => {}
        Some(MaskReadError::Ncbi(error)) if case.segmasker => {
            stderr.push_str(&String::from_utf8(error.what()).unwrap());
            stderr.push('\n');
            status = 1;
        }
        Some(MaskReadError::Ncbi(error)) => {
            stderr.push_str(
                &String::from_utf8(error.report_default_ex(
                    106,
                    16,
                    "Application's execution failed",
                ))
                .unwrap(),
            );
            status = 3;
        }
        Some(other) => problems.push(format!("error {other:?}")),
    }
    // NCBI's standard error (it names NCBI's own build trees) in LOSAT's words
    // (`crate::wording::losat_text`).
    if stderr != crate::wording::losat_text(&err) {
        problems.push(format!("standard error\n{stderr}NCBI\n{err}"));
    }
    if rc.trim() != status.to_string() {
        problems.push(format!("status {status}, NCBI {}", rc.trim()));
    }
    problems
}

fn fixture_dir(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/dustmasker")
        .join(name)
}

/// `Path::canonicalize`, which WASI lacks (it returns "unsupported" there): the paths here are
/// built from one absolute directory without links, so their components name them.
fn canonical(path: &Path) -> PathBuf {
    #[cfg(not(target_os = "wasi"))]
    {
        path.canonicalize().unwrap()
    }
    #[cfg(target_os = "wasi")]
    {
        path.components().collect()
    }
}

/// Every case of both fixture sets that the reader decides, and every input of
/// `tests/dustmasker/in` and of the reader's probes among them.
#[test]
fn every_input_reads_as_the_masking_programs_read_it() {
    let dirs = [fixture_dir(""), fixture_dir("reader")];
    let mut failures = Vec::new();
    let mut covered = std::collections::BTreeSet::new();
    let mut checked = 0;
    for dir in &dirs {
        for case in cases(dir) {
            let problems = check(&case);
            checked += 1;
            if problems.is_empty() {
                covered.insert(canonical(&dir.join(&case.input)));
            } else {
                failures.push(format!("{}: {}", case.name, problems.join("; ")));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert!(checked > 150, "{checked} cases");
    for dir in &dirs {
        for entry in std::fs::read_dir(dir.join("in")).unwrap() {
            let path = canonical(&entry.unwrap().path());
            assert!(
                covered.contains(&path),
                "no case decides {}",
                path.display()
            );
        }
    }
}

// The cases that the task names (line RG K1): letters (`<input>_fasta`), the error outcomes
// of short.fa, noheader.fa and empty.fa, and the ids under `-parse_seqids`.
#[test]
fn the_named_cases_are_decided() {
    let dir = fixture_dir("");
    let names: Vec<String> = cases(&dir).into_iter().map(|case| case.name).collect();
    for name in [
        "basic_fasta",
        "deflines_fasta",
        "ambig_fasta",
        "gaps_fasta",
        "protein_fasta",
        "badchar_fasta",
        "width80_fasta",
        "nruns_fasta",
        "short_int",
        "noheader_int",
        "empty_int",
        "deflines_int_ps",
        "mdefl_interval_ps",
        "dup_int_ps",
        "basic_int_ps",
        "short_int_ps",
    ] {
        assert!(names.iter().any(|n| n == name), "{name} is not decided");
    }
}

#[test]
fn letters_of_nucleotide_and_protein_records() {
    // A `U` is `T` (NCBI2na and NCBI4na hold no `U`), a hyphen and a `>?` gap are `N` in a
    // nucleotide sequence; lower case is upper case; bad letters are dropped.
    let read = run(
        b">a\nACGTU-acgtuRYKMSWBDHVN\nAC*EFGT\n>?3\nTT\n",
        false,
        false,
    );
    assert!(read.error.is_none());
    assert_eq!(read.records[0].letters, b"ACGTTNACGTTRYKMSWBDHVNACGTNNNTT");
    // In a protein `*` and `-` are `X`, as is a `>?` gap; `J`, `O`, `U` stay.
    let read = run(b">p\nMKV*-JOUBZXacd\n>?3\nAC\n", true, false);
    assert_eq!(read.records[0].letters, b"MKVXXJOUBZXACDXXXAC");
}

#[test]
fn ids_and_titles() {
    let input = b">\nACGT\n>  x  two words  \nACGT\n>lcl|r:1-5 ranged\nACGT\n>123456 t\nACGT\n";
    let read = run(input, false, false);
    let local = |n| vec![ParsedSeqId::Local(ObjectId::Id(n))];
    assert_eq!(
        read.records
            .iter()
            .map(|r| r.ids.clone())
            .collect::<Vec<_>>(),
        [local(1), local(2), local(3), local(4)]
    );
    assert_eq!(read.records[0].title, None);
    assert_eq!(read.records[1].title.as_deref(), Some(&b"x  two words"[..]));
    let read = run(input, false, true);
    assert_eq!(read.records[0].ids, local(1));
    assert_eq!(read.records[1].title.as_deref(), Some(&b"two words"[..]));
    assert_eq!(read.records[2].ids, local(2));
    assert_eq!(read.records[3].ids, local(123456));
}

#[test]
fn blank_and_comment_lines_only_end_without_error() {
    for input in [&b""[..], b"\n\n", b"\n;only\n#x\n\n"] {
        let read = run(input, false, false);
        assert!(read.records.is_empty() && read.error.is_none());
    }
}

#[test]
fn the_error_of_a_record_without_residues() {
    let read = run(b">a\nACGT\n>b\n>c\nACGT\n", false, false);
    assert_eq!(read.records.len(), 1);
    let Some(MaskReadError::Ncbi(error)) = &read.error else {
        panic!("{:?}", read.error)
    };
    assert_eq!(error.get_msg(), b"FASTA-Reader: No residues given");
    assert_eq!(
        String::from_utf8(error.report_default_ex(106, 16, "Application's execution failed"))
            .unwrap(),
        // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
        crate::wording::losat_text(
            "Error: (CObjReaderParseException::eNoResidues) FASTA-Reader: No residues given (m_Pos = 3)\n\
             Error: (106.16) Application's execution failed (CObjReaderParseException::eNoResidues) FASTA-Reader: No residues given (m_Pos = 3)\n"
        )
    );
}

/// A record of 4,294,967,296 letters or more, counting its gaps (NCBI's length is a 32-bit
/// `TSeqPos`; line RG DECISIONS K3-6, K3-8, KC-11d), reads as NCBI's object manager reads it:
/// the first (length modulo 2^32) letters, the gap lines as runs of `N` or `X`; gaps that sum to
/// a multiple of 2^32 stay pieces, which read as `N` or `X` where the map resolves them and
/// throw where a gap is 4294967295 or a piece ends past 32 bits. Bytes of dustmasker and
/// segmasker 2.17.0 (the official release, remaining-gaps-20261112/K/kc/p11d/p2/local).
#[test]
fn a_record_past_32_bits() {
    let letters = |run: &Run| -> Vec<String> {
        run.records
            .iter()
            .map(|record| String::from_utf8(record.letters.clone()).unwrap())
            .collect()
    };
    // The run of `N` is made, then cut to the length, 0.
    let read = run(b">a\nACGT\n>b\nA\n>?4294967295\n", false, false);
    assert!(read.error.is_none(), "{:?}", read.error);
    assert_eq!(letters(&read), ["ACGT", ""]);
    // Cut to 1 (the protein's `X` after it are not read).
    let read = run(b">a\nMKV\n>b\nM\n>?4294967295\n>?4294967295\nKK\n", true, false);
    assert!(read.error.is_none(), "{:?}", read.error);
    assert_eq!(letters(&read), ["MKV", "M"]);
    // Gaps that sum to 2^32: pieces, read as N or X (dustmasker -outfmt fasta:
    // `ACGTACGTACGTnnnnnnnn` and `nnnnnn`; segmasker: `XXXXXX`).
    let read = run(
        b">a\nACGTACGTACGT\n>?2147483648\n>?2147483648\nACGTaaaa\n>b\nacgtAC\n",
        false,
        false,
    );
    assert!(read.error.is_none(), "{:?}", read.error);
    assert_eq!(letters(&read), ["ACGTACGTACGTNNNNNNNN", "ACGTAC"]);
    let read = run(b">a\n>?2147483648\n>?2147483648\nACGTAC\n", false, false);
    assert_eq!(letters(&read), ["NNNNNN"]);
    let read = run(b">a\n>?2147483648\n>?2147483648\nMKVLAG\n", true, false);
    assert_eq!(letters(&read), ["XXXXXX"]);
    let read = run(b">a\nMKVL\n>?4294967283\n>?13\nQQQQ\n", true, false);
    assert_eq!(letters(&read), ["MKVLXXXX"]);
    // A gap of 4294967295 within the length, and a piece that ends past 32 bits, throw
    // (`Invalid sequence length`, `Sequence position overflow`); the records before are read.
    for (input, protein, error, length) in [
        (
            &b">z\nACGT\n>a\nACGTacgtACGT\n>?4294967295\n>?1\nACGTaaaa\n"[..],
            false,
            SeqMapError::InvalidLength,
            20,
        ),
        (
            b">z\nACGT\n>a\nACGTACGT\n>?1\n>?4294967295\nACGT\n",
            false,
            SeqMapError::InvalidLength,
            12,
        ),
        (
            b">z\nACGT\n>a\nACGTACGTAC\n>?4294967290\n>?6\nACGTAAAA\n",
            false,
            SeqMapError::PositionOverflow,
            18,
        ),
        (
            b">z\nMKVL\n>a\nMKVLMKVLMK\n>?4294967290\n>?6\nQQQQ\n",
            true,
            SeqMapError::PositionOverflow,
            14,
        ),
    ] {
        let read = run(input, protein, false);
        assert_eq!(read.records.len(), 1);
        let Some(MaskReadError::SeqMap {
            error: found,
            length: found_length,
            id,
        }) = read.error
        else {
            panic!("{:?}", read.error)
        };
        assert_eq!((found, found_length, id.as_str()), (error, length, "lcl|2"));
    }
    // A raw sequence of 4294967295 letters (`kInvalidSeqPos`) throws at once; one letter
    // fewer is read (its letters are not built here: 4 GiB).
    let raw = |length| SeqInst::Raw {
        length,
        seq_data: SeqData::Ncbi4na(Vec::new()),
    };
    assert_eq!(
        seq_map_error(&raw(4_294_967_295)),
        Some(SeqMapError::InvalidRawLength)
    );
    assert_eq!(seq_map_error(&raw(4_294_967_294)), None);
    assert_eq!(
        String::from_utf8(seq_map_exception(SeqMapError::InvalidRawLength).ncbi_what()).unwrap(),
        format!(
            "NCBI C++ Exception:\n    T0 \"{}/src/objmgr/seq_map.cpp\", line 367: Error: \
             (CSeqMapException::eSeqDataError) Invalid sequence length\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        )
    );
    assert_eq!(
        String::from_utf8(seq_vector_exception(14).ncbi_what()).unwrap(),
        format!(
            "NCBI C++ Exception:\n    T0 \"{}/src/objmgr/seq_vector_ci.cpp\", line 731: Error: \
             (CSeqVectorException::eSeqDataError) ncbi::objects::CSeqVector_CI::GetSeqData() - \
             CSeqVector_CI::GetSeqData: cannot get seq-data in range: 0-14\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        )
    );
}

/// Past 2^32 residues (here past 64) the record keeps those since its position last wrapped,
/// and an append of `AssembleSeq` from past them is `std::out_of_range` (line RG DECISIONS
/// KC-11d).
#[test]
fn residues_past_32_bits() {
    let read = |input: &[u8]| {
        let mut reader = MaskFastaReader::new(input, false, false)
            .unwrap()
            .with_seq_pos_modulus(64);
        let mut records = Vec::new();
        let error = loop {
            match reader.next(&mut |_| {}) {
                Ok(Some(record)) => records.push(String::from_utf8(record.letters).unwrap()),
                Ok(None) => break None,
                Err(error) => break Some(error),
            }
        };
        (records, error)
    };
    let mut input = b">a\n".to_vec();
    input.extend(b"A".repeat(60).chunks(10).flat_map(|line| [line, b"\n"].concat()));
    input.extend_from_slice(b"ACGTAC\nGGGG\n>b\nAC\n");
    let (records, error) = read(&input);
    assert!(error.is_none());
    assert_eq!(records, ["ACGGGG", "AC"]);
    // A gap at 20, past the 12 kept, and one at 71 (7): the runs after the residues.
    let mut input = b">a\n".to_vec();
    input.extend(b"A".repeat(20).chunks(10).flat_map(|line| [line, b"\n"].concat()));
    input.extend_from_slice(b">?3\n");
    input.extend(b"C".repeat(51).chunks(10).flat_map(|line| [line, b"\n"].concat()));
    input.extend_from_slice(b">?2\nGGGGG\n");
    assert_eq!(read(&input).0, ["CCCCCCCGGGGGNNNNN"]);
    // Gaps at 20 and 30, past the 6 kept.
    let mut input = b">a\n".to_vec();
    input.extend(b"A".repeat(20).chunks(10).flat_map(|line| [line, b"\n"].concat()));
    input.extend_from_slice(b">?3\nCCCCCCCCCC\n>?4\n");
    input.extend(b"G".repeat(40).chunks(10).flat_map(|line| [line, b"\n"].concat()));
    let (records, error) = read(&input);
    assert!(records.is_empty());
    let Some(MaskReadError::Std(what)) = error else {
        panic!("{error:?}")
    };
    assert_eq!(
        what,
        "basic_string::append: __pos (which is 20) > this->size() (which is 6)"
    );
}

struct FailingRead<'a>(&'a [u8]);

impl Read for FailingRead<'_> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if self.0.is_empty() {
            return Err(std::io::Error::other("the disk is gone"));
        }
        let count = buffer.len().min(self.0.len());
        buffer[..count].copy_from_slice(&self.0[..count]);
        self.0 = &self.0[count..];
        Ok(count)
    }
}

#[test]
fn a_failed_read_is_the_read_error() {
    assert_eq!(
        MaskFastaReader::new(FailingRead(b""), false, false)
            .err()
            .map(|error| error.to_string()),
        Some("the disk is gone".to_string())
    );
    let outcome = MaskFastaReader::new(FailingRead(b">a\nACGT\n>b\nGG\n"), false, false)
        .map_err(MaskReadError::Read)
        .and_then(|mut reader| {
            let first = reader.next(&mut |_| ());
            assert!(matches!(reader.next(&mut |_| ()), Ok(None)));
            first
        });
    match outcome {
        Err(MaskReadError::Read(error)) => assert_eq!(error.to_string(), "the disk is gone"),
        other => panic!("{other:?}"),
    }
}

/// NCBI's packing chooses NCBI2na for `A`, `C`, `G`, `T`, `U` and NCBI4na otherwise; both
/// give the same IUPACna letter for those five (so the letters do not depend on the
/// arrangement of literals). The reader stores the residues in upper case.
#[test]
fn both_codings_give_the_same_letters() {
    for &base in b"ACGTU" {
        assert_eq!(
            IUPACNA_FROM[0][usize::from(base)],
            IUPACNA_FROM[1][usize::from(base)]
        );
    }
    let mut odd = *b"-N@";
    read_back(&mut odd, Coding::Ncbi4na);
    assert_eq!(&odd, b"NNN");
}

/// Each residue reads back in place as it reads back from the data that `pack` makes of a
/// literal: every byte, at every place of a packed byte and in a last, partial byte.
#[test]
fn letters_read_back_in_place_as_from_the_packed_data() {
    let residues: Vec<u8> = (0..=255u8).chain(0..=255u8).skip(1).take(509).collect();
    for shift in 0..4 {
        let literal = &residues[shift..];
        let packed = pack::iupacna_to_ncbi2na(literal);
        let unpacked: Vec<u8> = (0..literal.len())
            .map(|i| b"ACGT"[usize::from(packed[i / 4] >> (6 - 2 * (i % 4)) & 3)])
            .collect();
        let mut letters = literal.to_vec();
        read_back(&mut letters, Coding::Ncbi2na);
        assert_eq!(letters, unpacked, "NCBI2na from {shift}");
        let packed = pack::iupacna_to_ncbi4na(literal);
        let unpacked: Vec<u8> = (0..literal.len())
            .map(|i| b"NACMGRSVTWYHKDBN"[usize::from(packed[i / 2] >> (4 - 4 * (i % 2)) & 0xf)])
            .collect();
        let mut letters = literal.to_vec();
        read_back(&mut letters, Coding::Ncbi4na);
        assert_eq!(letters, unpacked, "NCBI4na from {shift}");
    }
}

/// A raw sequence (one literal) and a delta (literals of both codings) give the letters in
/// the buffer of the residues; the reader then reads the next record as before.
#[test]
fn raw_and_delta_records_and_the_next_one() {
    let read = run(
        b">r\nACGTACGTAC\n>d\nACGTNNNNNNNNNNNNNNNNACGTU\n>r2\nRYKM\n",
        false,
        false,
    );
    assert!(read.error.is_none());
    let letters: Vec<&[u8]> = read.records.iter().map(|r| &r.letters[..]).collect();
    assert_eq!(
        letters,
        [&b"ACGTACGTAC"[..], b"ACGTNNNNNNNNNNNNNNNNACGTT", b"RYKM"]
    );
}
