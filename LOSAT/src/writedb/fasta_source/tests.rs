//! Tests of the FASTA source against makeblastdb 2.17.0 (`-blastdb_version 4`, no
//! `-parse_seqids`).
//!
//! The expected values are what makeblastdb printed and what its databases hold
//! (`makeblastdb_cases`, and the values quoted in the tests). The tests add the part of
//! `CBuildDatabase::AddSequences` that decides whether a sequence is stored (`run`), so
//! that a whole input compares with makeblastdb's messages and database.

mod makeblastdb_cases;

use std::io::Read;

use super::line_reader::BUFFER_SIZE;
use super::*;
use crate::seqdb::{SeqDb, SeqType};

/// A sequence of a makeblastdb database.
pub(super) struct Rec {
    /// The title of the header file.
    pub(super) title: &'static [u8],
    pub(super) length: u64,
    /// The residues as letters (NCBIstdaa or the bases with their ambiguity codes).
    pub(super) residues: &'static str,
    /// Whether the sequence has an ambiguity block (possibly an empty one).
    pub(super) ambiguity: bool,
}

/// One input and what makeblastdb made of it.
pub(super) struct Case {
    pub(super) name: &'static str,
    pub(super) protein: bool,
    pub(super) input: &'static [u8],
    /// The reader's warnings, in order.
    pub(super) warnings: &'static [&'static str],
    /// `Ignoring sequence '<id>' as it has no sequence data`.
    pub(super) ignored: &'static [&'static str],
    /// The text after `BLAST Database creation error: `.
    pub(super) error: Option<&'static str>,
    /// The sequences of the database; `None` when makeblastdb wrote none.
    pub(super) records: Option<&'static [Rec]>,
}

/// What makeblastdb reports for an input, and the sequences it stores.
#[derive(Debug, PartialEq, Eq)]
struct Outcome {
    warnings: Vec<Vec<u8>>,
    ignored: Vec<String>,
    error: Option<String>,
    added: Vec<Bioseq>,
}

/// What `CBuildDatabase::AddSequences` and the writer do with a sequence of the source,
/// as far as it decides what is stored.
enum Fate {
    Stored,
    /// `Ignoring sequence '<id>' as it has no sequence data`.
    Ignored(String),
    /// The text after `BLAST Database creation error: `.
    Error(String),
}

/// The tests of `AddSequences` (build_db.cpp:864-876: the length limit, then
/// `Ignoring sequence` for a length of 0 or when `x_EditAndAddBioseq` finds no sequence
/// data) and the writer's error for sequence data it cannot read
/// (writedb_impl.cpp:1032-1038).
fn fate(bioseq: &Bioseq) -> Fate {
    if bioseq.length > 0x7fff_ffff {
        return Fate::Error(format!(
            "Sequences longer than 2,147,483,647 bases are not supported.  Offending sequence is lcl|{}",
            bioseq.local_id
        ));
    }
    if bioseq.length == 0 || bioseq.seq_data.is_none() {
        return Fate::Ignored(format!("lcl|{}", bioseq.local_id));
    }
    if bioseq.seq_data == Some(SeqData::Other) {
        return Fate::Error(format!(
            "Unable to process sequence for entry [{}].",
            bioseq.local_id
        ));
    }
    Fate::Stored
}

/// Reads `input` to its end as `CBuildDatabase::AddFasta` does: the source, `fate` for
/// each sequence, and `No sequences added` (build_db.cpp:1410-1412).
fn run<R: Read>(input: R, protein: bool) -> Outcome {
    let mut source = FastaBioseqSource::new(input, protein);
    let mut outcome = Outcome {
        warnings: Vec::new(),
        ignored: Vec::new(),
        error: None,
        added: Vec::new(),
    };
    loop {
        let warnings = &mut outcome.warnings;
        let bioseq = match source.next_bioseq(&mut |line: &[u8]| warnings.push(line.to_vec())) {
            Ok(Some(bioseq)) => bioseq,
            Ok(None) => break,
            Err(error) => {
                outcome.error = Some(error.message);
                return outcome;
            }
        };
        assert_eq!(bioseq.is_aa, protein);
        match fate(&bioseq) {
            Fate::Stored => outcome.added.push(bioseq),
            Fate::Ignored(id) => outcome.ignored.push(id),
            Fate::Error(message) => {
                outcome.error = Some(message);
                return outcome;
            }
        }
    }
    assert_eq!(source.next_bioseq(&mut |_| ()), Ok(None));
    if outcome.added.is_empty() {
        outcome.error = Some("No sequences added".to_string());
    }
    outcome
}

/// A reader that returns at most a few bytes for each call.
struct ShortReads<'a> {
    rest: &'a [u8],
    step: usize,
}

impl Read for ShortReads<'_> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        self.step = self.step % 7 + 1;
        let count = self.step.min(buffer.len()).min(self.rest.len());
        buffer[..count].copy_from_slice(&self.rest[..count]);
        self.rest = &self.rest[count..];
        Ok(count)
    }
}

/// The residues of a stored sequence as letters, with the checks of the packed forms: the
/// number of bytes, and zero bits after the last base.
fn residues(bioseq: &Bioseq) -> String {
    let length = bioseq.length as usize;
    let letters: Vec<u8> = match bioseq.seq_data.as_ref().expect("sequence data") {
        SeqData::Ncbieaa(bytes) => {
            assert_eq!(bytes.len(), length);
            bytes.clone()
        }
        SeqData::Ncbi2na(bytes) => {
            assert_eq!(bytes.len(), length.div_ceil(4));
            if length % 4 != 0 {
                assert_eq!(bytes[length / 4] << (2 * (length % 4)), 0);
            }
            (0..length)
                .map(|i| b"ACGT"[usize::from(bytes[i / 4] >> (6 - 2 * (i % 4)) & 3)])
                .collect()
        }
        SeqData::Ncbi4na(bytes) => {
            assert_eq!(bytes.len(), length.div_ceil(2));
            if length % 2 != 0 {
                assert_eq!(bytes[length / 2] & 15, 0);
            }
            (0..length)
                .map(|i| b"-ACMGRSVTWYHKDBN"[usize::from(bytes[i / 2] >> (4 - 4 * (i % 2)) & 15)])
                .collect()
        }
        other => panic!("unexpected sequence data {other:?}"),
    };
    String::from_utf8(letters).expect("letters")
}

/// Whether the writer gives the sequence an ambiguity block (`WriteDB_Ncbi4naToBinary`).
fn has_ambiguity_block(bioseq: &Bioseq) -> bool {
    matches!(bioseq.seq_data, Some(SeqData::Ncbi4na(_)))
}

/// The title as the header file has it: the ASN.1 writer replaces the bytes that a
/// VisibleString cannot hold by `#`.
fn header_title(title: &[u8]) -> Vec<u8> {
    title
        .iter()
        .map(|&byte| {
            if (0x20..=0x7e).contains(&byte) {
                byte
            } else {
                b'#'
            }
        })
        .collect()
}

/// The differences between an outcome and makeblastdb's.
fn differences(case: &Case, outcome: &Outcome) -> Vec<String> {
    let mut found = Vec::new();
    let warnings: Vec<String> = outcome
        .warnings
        .iter()
        .map(|line| String::from_utf8_lossy(line).into_owned())
        .collect();
    if warnings != case.warnings {
        found.push(format!(
            "warnings {warnings:?}, makeblastdb {:?}",
            case.warnings
        ));
    }
    if outcome.ignored != case.ignored {
        found.push(format!(
            "ignored {:?}, makeblastdb {:?}",
            outcome.ignored, case.ignored
        ));
    }
    if outcome.error.as_deref() != case.error {
        found.push(format!(
            "error {:?}, makeblastdb {:?}",
            outcome.error, case.error
        ));
    }
    // A database is written only without an error.
    let Some(records) = case.records else {
        return found;
    };
    if outcome.added.len() != records.len() {
        found.push(format!(
            "{} sequences, makeblastdb {}",
            outcome.added.len(),
            records.len()
        ));
        return found;
    }
    for (index, (bioseq, record)) in outcome.added.iter().zip(records).enumerate() {
        if header_title(&bioseq.title) != record.title {
            found.push(format!(
                "sequence {index}: title {:?}, makeblastdb {:?}",
                String::from_utf8_lossy(&bioseq.title),
                String::from_utf8_lossy(record.title)
            ));
        }
        if bioseq.length != record.length {
            found.push(format!(
                "sequence {index}: length {}, makeblastdb {}",
                bioseq.length, record.length
            ));
            continue;
        }
        let letters = residues(bioseq);
        if letters != record.residues {
            found.push(format!(
                "sequence {index}: residues {letters}, makeblastdb {}",
                record.residues
            ));
        }
        if has_ambiguity_block(bioseq) != record.ambiguity {
            found.push(format!(
                "sequence {index}: ambiguity block {}, makeblastdb {}",
                has_ambiguity_block(bioseq),
                record.ambiguity
            ));
        }
    }
    found
}

/// Every small case of the probe (`db2_edge_probe.py` and `db2_fs_cases.py`): the
/// warnings, the ignored sequences, the error, and each sequence of the database.
#[test]
fn makeblastdb_cases_match() {
    let mut failures = Vec::new();
    for case in makeblastdb_cases::CASES {
        let outcome = run(case.input, case.protein);
        for difference in differences(case, &outcome) {
            failures.push(format!("{}: {difference}", case.name));
        }
        let short = run(
            ShortReads {
                rest: case.input,
                step: 0,
            },
            case.protein,
        );
        if short != outcome {
            failures.push(format!("{}: short reads change the result", case.name));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert!(makeblastdb_cases::CASES.len() > 150);
}

fn one(input: &[u8], protein: bool) -> Bioseq {
    let mut source = FastaBioseqSource::new(input, protein);
    let bioseq = source
        .next_bioseq(&mut |_| ())
        .expect("no error")
        .expect("a sequence");
    assert_eq!(source.next_bioseq(&mut |_| ()), Ok(None));
    bioseq
}

/// The `Seq-data` that reaches the writer: its choice and bytes.
#[test]
fn seq_data_choice_and_bytes() {
    let data = |input: &[u8], protein| {
        let bioseq = one(input, protein);
        (bioseq.length, bioseq.seq_data)
    };
    // One or two unambiguous bases are NCBI4na (the writer adds an empty ambiguity block:
    // makeblastdb's `G` is `81 00000000`), three or more NCBI2na.
    assert_eq!(
        data(b">a\nG\n", false),
        (1, Some(SeqData::Ncbi4na(vec![0x40])))
    );
    assert_eq!(
        data(b">a\nGA\n", false),
        (2, Some(SeqData::Ncbi4na(vec![0x41])))
    );
    assert_eq!(
        data(b">a\nGAT\n", false),
        (3, Some(SeqData::Ncbi2na(vec![0x8c])))
    );
    assert_eq!(
        data(b">a\nacgt\n", false),
        (4, Some(SeqData::Ncbi2na(vec![0x1b])))
    );
    // `U` is `T`.
    assert_eq!(
        data(b">a\nACGUUu\n", false),
        (6, Some(SeqData::Ncbi2na(vec![0x1b, 0xf0])))
    );
    // Any ambiguity code or gap makes the whole sequence NCBI4na.
    assert_eq!(
        data(b">a\nACGTN\n", false),
        (5, Some(SeqData::Ncbi4na(vec![0x12, 0x48, 0xf0])))
    );
    assert_eq!(
        data(b">a\nAC-GT\n", false),
        (5, Some(SeqData::Ncbi4na(vec![0x12, 0xf4, 0x80])))
    );
    assert_eq!(
        data(b">a\nAC\n>?3\nGT\n", false),
        (7, Some(SeqData::Ncbi4na(vec![0x12, 0xff, 0xf4, 0x80])))
    );
    // A lone hyphen at the end of a line is a gap of unknown length: one `N`, and the
    // sequence is one longer than the reader's.
    assert_eq!(
        data(b">a\nACGT-\n", false),
        (5, Some(SeqData::Ncbi4na(vec![0x12, 0x48, 0xf0])))
    );
    // ... but first in the record it leaves a delta that `s_FixBioseqDeltas` does not
    // read: no sequence data (makeblastdb: "Ignoring sequence").
    assert_eq!(data(b">a\n-\nACGT\n", false), (4, None));
    assert_eq!(
        data(b">a\n-\n ACGT\n", false),
        (5, Some(SeqData::Ncbi4na(vec![0xf1, 0x24, 0x80])))
    );
    // A record that is one gap has sequence data without a choice (makeblastdb: "Unable
    // to process sequence for entry"), of length 0 for a lone hyphen (ignored).
    assert_eq!(data(b">a\n>?10\n", false), (10, Some(SeqData::Other)));
    assert_eq!(data(b">a\n--\n", false), (2, Some(SeqData::Other)));
    assert_eq!(data(b">a\n-\n", false), (0, Some(SeqData::Other)));
    assert_eq!(
        data(b">a\n>?1\n>?1\n", false),
        (2, Some(SeqData::Ncbi4na(vec![0xff])))
    );
    // No residues: no sequence data.
    assert_eq!(data(b">a\n", false), (0, None));
    assert_eq!(data(b">a\n", true), (0, None));
    // Protein: NCBIeaa as read, in upper case; `-` and `*` are residues, a `>?` gap a
    // run of `X`.
    assert_eq!(
        data(b">a\nmk*V-l\n", true),
        (6, Some(SeqData::Ncbieaa(b"MK*V-L".to_vec())))
    );
    assert_eq!(
        data(b">a\nMK\n>?3\nVL\n", true),
        (7, Some(SeqData::Ncbieaa(b"MKXXXVL".to_vec())))
    );
    assert_eq!(
        data(b">a\n>?2\n", true),
        (2, Some(SeqData::Ncbieaa(b"XX".to_vec())))
    );
}

/// The local IDs count the records of the input from 1, the empty ones too.
#[test]
fn local_ids_and_titles() {
    let mut source = FastaBioseqSource::new(&b">a one\nMKV\n>b two\n>\nGG\n"[..], true);
    let mut read = Vec::new();
    while let Some(bioseq) = source.next_bioseq(&mut |_| ()).unwrap() {
        read.push((bioseq.local_id, bioseq.title, bioseq.is_aa, bioseq.length));
    }
    assert_eq!(
        read,
        [
            ("1".to_string(), b"a one".to_vec(), true, 3),
            ("2".to_string(), b"b two".to_vec(), true, 0),
            ("3".to_string(), b"".to_vec(), true, 2),
        ]
    );
    // The title is the definition line after `>`: its bytes, with a tab as three spaces
    // and a ctrl-A as one (`title`).
    assert_eq!(
        one(b">  a\ttab \x01ctl \xc3\xa9 \x7f  \nACGT\n", false).title,
        b"  a   tab  ctl \xc3\xa9 \x7f  "
    );
    assert_eq!(one(b">?_gap t\nACGT\n", false).title, b"gap t");
    assert_eq!(one(b">a", false).title, b"a");
}

/// `prefix`, then lines of 60 `A` and a shorter last line, so that the input has `total`
/// bytes and ends with `tail_eol`.
fn pad_lines(prefix: &[u8], total: usize, eol: &[u8], tail_eol: &[u8]) -> Vec<u8> {
    let mut out = prefix.to_vec();
    while out.len() + 60 + eol.len() + 1 + tail_eol.len() <= total {
        out.extend_from_slice(&[b'A'; 60]);
        out.extend_from_slice(eol);
    }
    let rest = total - out.len() - tail_eol.len();
    assert!(rest >= 1);
    out.extend(std::iter::repeat_n(b'A', rest));
    out.extend_from_slice(tail_eol);
    assert_eq!(out.len(), total);
    out
}

/// The outcome of `input`, read at once and in short reads.
fn outcome(input: &[u8], protein: bool) -> Outcome {
    let outcome = run(input, protein);
    let short = run(
        ShortReads {
            rest: input,
            step: 0,
        },
        protein,
    );
    assert_eq!(short, outcome);
    outcome
}

/// The lengths of the stored sequences and the error.
fn lengths_and_error(input: &[u8], protein: bool) -> (Vec<u64>, Option<String>) {
    let outcome = outcome(input, protein);
    assert!(outcome.warnings.is_empty(), "{:?}", outcome.warnings);
    (
        outcome.added.iter().map(|bioseq| bioseq.length).collect(),
        outcome.error,
    )
}

/// Lines and the 32 KiB buffer of `CBufferedLineReader` (the probe cases `fs_p_buf_*`:
/// the lengths are those of makeblastdb's databases).
#[test]
fn lines_across_buffers() {
    let buffer = BUFFER_SIZE;
    assert_eq!(buffer, 32768);
    // A line feed that is the last byte of a buffer, and the first of the next.
    assert_eq!(
        lengths_and_error(&pad_lines(b">a\n", buffer, b"\n", b"\n"), true),
        (vec![32227], None)
    );
    assert_eq!(
        lengths_and_error(&pad_lines(b">a\n", buffer + 1, b"\n", b"\n"), true),
        (vec![32228], None)
    );
    // CR LF inside the end of a buffer; CR as its last byte.
    assert_eq!(
        lengths_and_error(&pad_lines(b">>a\r\n", buffer, b"\r\n", b"\r\n"), true),
        (vec![31705], None)
    );
    assert_eq!(
        lengths_and_error(&pad_lines(b">>a\r", buffer, b"\r", b"\r"), true),
        (vec![32226], None)
    );
    // CR LF across two buffers with more input after it.
    let mut more = pad_lines(b">>a\r\n", buffer + 1, b"\r\n", b"\r\n");
    more.extend_from_slice(b">b\r\nGG\r\n");
    assert_eq!(lengths_and_error(&more, true), (vec![31706, 2], None));
    // A definition line across two buffers.
    let mut straddle = pad_lines(b">a\n", buffer - 2, b"\n", b"\n");
    straddle.extend_from_slice(b">bcdef g\nGG\n");
    let read = outcome(&straddle, true);
    assert_eq!(read.added.len(), 2);
    assert_eq!(read.added[0].length, 32225);
    assert_eq!(read.added[1].title, b"bcdef g");
    // Lines longer than a buffer (`x_LoadLong`).
    let mut long = b">a\n".to_vec();
    long.extend(b"MKVLAAGG".repeat(20000));
    long.extend_from_slice(b"\n>b\nGG\n");
    assert_eq!(lengths_and_error(&long, true), (vec![160000, 2], None));
    let mut long = b">a\n".to_vec();
    long.extend(std::iter::repeat_n(b'M', 2 * buffer - 4));
    long.push(b'\n');
    assert_eq!(lengths_and_error(&long, true), (vec![65532], None));
    let mut no_eol = b">a\n".to_vec();
    no_eol.extend(std::iter::repeat_n(b'M', buffer - 3));
    assert_eq!(lengths_and_error(&no_eol, true), (vec![32765], None));
    // A title longer than a buffer.
    let mut title = b">".to_vec();
    title.extend(std::iter::repeat_n(b't', buffer + 100));
    title.extend_from_slice(b"\nMKV\n>b\nGG\n");
    let read = run(&title[..], true);
    assert_eq!(read.added[0].title.len(), buffer + 100);
    assert_eq!(read.warnings.len(), 1);
    assert_eq!(read.added[1].title, b"b");
}

/// After a CR that is the last byte of a buffer, when the next read brings only the LF
/// that follows it and the input ends there, NCBI's `PeekChar` returns the second byte
/// of the buffer before (`line_reader`): a `>` there makes the reader take the end of
/// the input for an empty line that starts a record, and the next `GetNext` fails.
/// makeblastdb (probe cases `fs_p_buf_crlf_straddle_eof*`, `fs_p_buf_cr_then_lf_only`,
/// `fs_p_buf_long_line_cr_straddle*`, `fs_n_buf_crlf_straddle_eof_gt`).
#[test]
fn stale_peek_after_cr_lf_across_buffers() {
    let buffer = BUFFER_SIZE;
    // The second byte of the input is `a`: nothing happens.
    assert_eq!(
        lengths_and_error(&pad_lines(b">a\r\n", buffer + 1, b"\r\n", b"\r\n"), true),
        (vec![31707], None)
    );
    // The second byte of the input is `>`.
    let input = pad_lines(b">>a\r\n", buffer + 1, b"\r\n", b"\r\n");
    for protein in [true, false] {
        let read = outcome(&input, protein);
        assert_eq!(
            read.error.as_deref(),
            Some("Expected defline around line 531")
        );
        assert_eq!(read.added.len(), 1);
        assert_eq!(read.added[0].length, 31706);
        assert_eq!(read.added[0].title, b">a");
    }
    // CR line ends, and CR LF at the end.
    assert_eq!(
        lengths_and_error(&pad_lines(b">>a\r", buffer + 1, b"\r", b"\r\n"), true).1,
        Some("Expected defline around line 540".to_string())
    );
    // Two buffers: the byte is the second one of the second buffer, an `A`.
    assert_eq!(
        lengths_and_error(
            &pad_lines(b">>a\r\n", 2 * buffer + 1, b"\r\n", b"\r\n"),
            true
        ),
        (vec![63418], None)
    );
    // The same at the end of a line longer than a buffer (`x_LoadLong`).
    let mut long = b">>a\r\n".to_vec();
    long.extend(std::iter::repeat_n(b'M', 2 * buffer - 6));
    long.extend_from_slice(b"\r\n");
    assert_eq!(long.len(), 2 * buffer + 1);
    assert_eq!(lengths_and_error(&long, true), (vec![65530], None));
    let mut long = b">>a\r\n".to_vec();
    long.extend(std::iter::repeat_n(b'M', buffer + 1 - 5));
    long.push(b'>');
    long.extend(std::iter::repeat_n(b'M', buffer - 3));
    long.extend_from_slice(b"\r\n");
    assert_eq!(long.len(), 2 * buffer + 1);
    let read = outcome(&long, true);
    assert_eq!(
        read.warnings,
        [b"FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 32765".to_vec()]
    );
    assert_eq!(
        read.error.as_deref(),
        Some("Expected defline around line 3")
    );
    // With a record after it the input does not end there.
    long.extend_from_slice(b">b\r\nGG\r\n");
    let read = outcome(&long, true);
    assert_eq!(read.error, None);
    assert_eq!(
        read.added
            .iter()
            .map(|bioseq| bioseq.length)
            .collect::<Vec<_>>(),
        [65529, 2]
    );
}

fn random_bases(count: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..count)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect()
}

fn fasta_lines(sequence: &[u8], width: usize) -> Vec<u8> {
    let mut out = b">r\n".to_vec();
    for line in sequence.chunks(width) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
    out
}

/// Long sequences: unambiguous bases are NCBI2na whatever their number (makeblastdb:
/// no ambiguity block for 1 Mb and 10 Mb of random bases, probe cases
/// `fs_n_random_*`), and one `N` makes the whole sequence NCBI4na.
#[test]
fn long_nucleotide_sequences() {
    let bases = random_bases(1_000_003, 0x9e37_79b9_7f4a_7c15);
    for input in [fasta_lines(&bases, 80), fasta_lines(&bases, bases.len())] {
        let bioseq = one(&input, false);
        assert_eq!(bioseq.length, 1_000_003);
        assert!(matches!(bioseq.seq_data, Some(SeqData::Ncbi2na(_))));
        assert_eq!(residues(&bioseq).as_bytes(), bases);
    }
    let mut one_n = bases.clone();
    one_n[500_000] = b'N';
    let bioseq = one(&fasta_lines(&one_n, 80), false);
    assert!(has_ambiguity_block(&bioseq));
    assert_eq!(residues(&bioseq).as_bytes(), one_n);
    // Runs of `N` between long unambiguous runs: `Pack` makes several literals, and
    // `s_FixBioseqDeltas` one NCBI4na sequence of them (probe cases `n_n_runs`, `fs_n_mix`).
    let mut runs = Vec::new();
    for run in [15, 16, 4095, 4096, 4097] {
        runs.extend_from_slice(&bases[..3000]);
        runs.extend(std::iter::repeat_n(b'N', run));
    }
    runs.extend_from_slice(&bases[..7]);
    let bioseq = one(&fasta_lines(&runs, 70), false);
    assert_eq!(bioseq.length, runs.len() as u64);
    assert!(has_ambiguity_block(&bioseq));
    assert_eq!(residues(&bioseq).as_bytes(), runs);
}

/// A record above 2,147,483,647 residues comes with its length and no sequence data
/// (makeblastdb, probe case `fs_n_gap_too_long`: "Sequences longer than 2,147,483,647
/// bases are not supported.  Offending sequence is lcl|1").
#[test]
fn record_above_the_length_limit() {
    for (input, protein) in [
        (&b">a\nA\n>?2147483647\n"[..], false),
        (&b">a\nM\n>?2147483647\n"[..], true),
        (&b">a\n>?4294967295\n>?4294967295\n"[..], false),
    ] {
        let bioseq = one(input, protein);
        assert!(bioseq.length > 0x7fff_ffff);
        assert_eq!(bioseq.seq_data, None);
        assert!(matches!(fate(&bioseq), Fate::Error(_)));
    }
    assert_eq!(one(b">a\nA\n>?2147483647\n", false).length, 2_147_483_648);
    // At the limit a record that is one gap is not built either way.
    let bioseq = one(b">a\n>?2147483647\n", false);
    assert_eq!(
        (bioseq.length, bioseq.seq_data),
        (2_147_483_647, Some(SeqData::Other))
    );
}

/// A reader that fails after its bytes.
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
fn read_error_ends_the_reading() {
    let mut source = FastaBioseqSource::new(FailingRead(b">a\nMKV\n>b\nGG\n"), true);
    let first = source.next_bioseq(&mut |_| ());
    assert_eq!(
        first,
        Err(SourceError {
            message: "Read error: the disk is gone".to_string(),
            out_of_memory: false,
        })
    );
    assert_eq!(source.next_bioseq(&mut |_| ()), Ok(None));
}

/// A delta whose gap is longer than any allocation is the error that makeblastdb reports
/// as "BLAST ran out of memory", not an abort. (The reader hands over no length above
/// 2,147,483,647; the length here only makes the reservation fail on every machine.)
#[test]
fn a_delta_that_cannot_be_allocated_is_an_error() {
    let delta = |gap: u64| SeqInst::Delta {
        length: 1 + gap,
        residues: b"A".to_vec(),
        pieces: vec![
            DeltaSeq::Literal(pack::Segment {
                coding: pack::Coding::Ncbi4na,
                start: 0,
                length: 1,
            }),
            DeltaSeq::Gap { length: gap },
        ],
    };
    assert_eq!(
        fix_bioseq_deltas(delta(2)).unwrap(),
        (3, Some(SeqData::Ncbi4na(vec![0x1f, 0xf0])))
    );
    assert!(fix_bioseq_deltas(delta(u64::MAX - 1)).is_err());
    let error = SourceError::out_of_memory();
    assert!(error.out_of_memory);
    assert_eq!(error.message, "BLAST ran out of memory");
}

/// NCBIstdaa: the residue at position `i` has the code `i`
/// (c++/src/algo/blast/core/blast_encoding.c:115-118, `NCBISTDAA_TO_AMINOACID`).
const NCBISTDAA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
/// BLASTNA: the base at position `i` has the code `i`
/// (c++/src/algo/blast/core/blast_encoding.c:89-106, `BLASTNA_TO_IUPACNA`).
const BLASTNA: &[u8; 16] = b"ACGTRYMKWSBDHVN-";

/// Compares FASTA files outside the repository with the databases that makeblastdb
/// (`-blastdb_version 4`, no `-parse_seqids`) built from them.
///
/// `LOSAT_FASTA_SOURCE_CHECK` names a list: one line for each input, with the
/// tab-separated columns name, `prot` or `nucl`, the FASTA file, the database (`-` when
/// makeblastdb wrote none) and a file for the transcript. The source reads each FASTA
/// file; every sequence that makeblastdb stores (`fate`) is compared with the sequence
/// of the same position in the database: the title (with `#` for the bytes outside
/// 0x20..0x7E), the length, the residues, and whether there is an ambiguity block. The
/// transcript gets the warnings, the ignored sequences and the error, for the
/// comparison with makeblastdb's messages (`db2_fs_check.py` of the task folder).
#[test]
#[ignore = "needs LOSAT_FASTA_SOURCE_CHECK (FASTA files and databases outside the repository)"]
fn fasta_files_match_their_makeblastdb_databases() {
    use std::io::Write;

    let list = std::env::var("LOSAT_FASTA_SOURCE_CHECK").expect("LOSAT_FASTA_SOURCE_CHECK");
    let text = std::fs::read_to_string(&list).expect("the list");
    let mut failed = Vec::new();
    for line in text.lines().filter(|line| !line.is_empty()) {
        let columns: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            columns.len(),
            5,
            "name, prot|nucl, FASTA, database, transcript"
        );
        let (name, protein) = (columns[0], columns[1] == "prot");
        let started = std::time::Instant::now();
        let file = std::fs::File::open(columns[2]).unwrap_or_else(|e| panic!("{name}: {e}"));
        let mut source = FastaBioseqSource::new(file, protein);
        let seqtype = if protein {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        };
        let db = (columns[3] != "-")
            .then(|| SeqDb::open(columns[3], seqtype).unwrap_or_else(|e| panic!("{name}: {e}")));
        let mut transcript = std::io::BufWriter::new(
            std::fs::File::create(columns[4]).unwrap_or_else(|e| panic!("{name}: {e}")),
        );
        let mut warnings: Vec<Vec<u8>> = Vec::new();
        let mut mismatches = Vec::new();
        let mut block = None;
        let mut block_end = 0u32;
        let (mut added, mut letters, mut with_ambiguity) = (0u32, 0u64, 0u32);
        let mut error = None;
        loop {
            let next = source.next_bioseq(&mut |line: &[u8]| warnings.push(line.to_vec()));
            for warning in warnings.drain(..) {
                transcript.write_all(b"warning\t").unwrap();
                transcript.write_all(&warning).unwrap();
                transcript.write_all(b"\n").unwrap();
            }
            let bioseq = match next {
                Ok(Some(bioseq)) => bioseq,
                Ok(None) => break,
                Err(source_error) => {
                    error = Some(source_error.message);
                    break;
                }
            };
            match fate(&bioseq) {
                Fate::Stored => {}
                Fate::Ignored(id) => {
                    writeln!(transcript, "ignored\t{id}").unwrap();
                    continue;
                }
                Fate::Error(message) => {
                    error = Some(message);
                    break;
                }
            }
            let ordinal = added;
            added += 1;
            letters += bioseq.length;
            with_ambiguity += u32::from(has_ambiguity_block(&bioseq));
            let Some(db) = &db else {
                continue;
            };
            if ordinal >= db.num_seqs() {
                continue;
            }
            if ordinal >= block_end {
                let range = db.next_block(ordinal, 16 << 20).expect("a block");
                block_end = range.end_ordinal();
                block = Some(db.read_block(range).expect("read block"));
            }
            let block = block.as_ref().expect("a block");
            let i = (ordinal - block.first_ordinal()) as usize;
            let deflines = db.deflines(ordinal).expect("deflines");
            if deflines.len() != 1 || deflines[0].title != header_title(&bioseq.title) {
                mismatches.push(format!("sequence {ordinal}: title"));
            }
            let (stored, stored_ambiguity): (Vec<u8>, bool) = if protein {
                (
                    block
                        .protein(i)
                        .iter()
                        .map(|&code| NCBISTDAA[usize::from(code)])
                        .collect(),
                    false,
                )
            } else {
                let blastna = block.nucleotide_blastna(i).expect("bases");
                (
                    blastna[1..blastna.len() - 1]
                        .iter()
                        .map(|&code| BLASTNA[usize::from(code)])
                        .collect(),
                    !block.ambiguity_bytes(i).is_empty(),
                )
            };
            if stored.len() as u64 != bioseq.length {
                mismatches.push(format!(
                    "sequence {ordinal}: length {}, database {}",
                    bioseq.length,
                    stored.len()
                ));
            } else if residues(&bioseq).as_bytes() != stored {
                mismatches.push(format!("sequence {ordinal}: residues"));
            }
            if has_ambiguity_block(&bioseq) != stored_ambiguity {
                mismatches.push(format!(
                    "sequence {ordinal}: ambiguity block {}, database {stored_ambiguity}",
                    has_ambiguity_block(&bioseq)
                ));
            }
        }
        if error.is_none() && added == 0 {
            error = Some("No sequences added".to_string());
        }
        if let Some(error) = &error {
            writeln!(transcript, "error\t{error}").unwrap();
        }
        match &db {
            // makeblastdb leaves a database only when it ends without an error.
            Some(db) => {
                if db.num_seqs() != added {
                    mismatches.push(format!("{added} sequences, database {}", db.num_seqs()));
                }
                if db.total_length() != letters {
                    mismatches.push(format!("{letters} letters, database {}", db.total_length()));
                }
                if error.is_some() {
                    mismatches.push(format!("error {error:?} with a database"));
                }
            }
            None => {
                if error.is_none() {
                    mismatches.push("no error without a database".to_string());
                }
            }
        }
        for mismatch in mismatches.iter().take(20) {
            writeln!(transcript, "mismatch\t{mismatch}").unwrap();
        }
        writeln!(
            transcript,
            "summary\t{added} sequences\t{letters} letters\t{with_ambiguity} with ambiguity block\t{} mismatches\t{:.2} s",
            mismatches.len(),
            started.elapsed().as_secs_f64()
        )
        .unwrap();
        transcript.flush().unwrap();
        if !mismatches.is_empty() {
            failed.push(format!("{name}: {}", mismatches[0]));
        }
    }
    assert!(failed.is_empty(), "{}", failed.join("\n"));
}
