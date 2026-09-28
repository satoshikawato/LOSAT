//! Property tests of the index scan against `bio::io::fasta` 1.6.0 (plan TD-8): the
//! records, IDs and lengths agree with the parser; the offsets, the line layout and the
//! residue counts agree with a reference walk that records where every residue is; and
//! scanning in chunks of any size gives the same result as scanning at once.
//!
//! Inputs cover LF and CRLF, uneven lines, blank lines, interior and trailing
//! whitespace (including Unicode whitespace), non-ASCII headers, empty records, invalid
//! UTF-8 and records larger than the checkpoint spacing.

use bio::io::fasta;
use losat_web_adapter::scan::{scan, LineLayout, ScanRecord, Scanner, CHECKPOINT_EVERY};

/// A small deterministic generator (xorshift64*).
struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }
    fn chance(&mut self, percent: u64) -> bool {
        self.below(100) < percent
    }
    fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.below(items.len() as u64) as usize]
    }
}

/// One record as the reference walk sees it.
#[derive(Debug)]
struct Expected {
    id: String,
    header_offset: u64,
    sequence_offset: u64,
    end_offset: u64,
    /// The input offset of every sequence byte.
    residue_offsets: Vec<u64>,
    residues: Vec<u8>,
}

/// bio's rules (see `src/scan.rs`), walked line by line with the offsets kept.
fn reference(input: &[u8]) -> Result<Vec<Expected>, String> {
    let mut lines = Vec::new();
    let mut start = 0;
    while start < input.len() {
        let end = input[start..]
            .iter()
            .position(|&byte| byte == b'\n')
            .map_or(input.len(), |at| start + at + 1);
        lines.push((start, &input[start..end]));
        start = end;
    }
    let mut records: Vec<Expected> = Vec::new();
    let mut current: Option<(Expected, bool)> = None;
    for (index, &(offset, line)) in lines.iter().enumerate() {
        let text = std::str::from_utf8(line).map_err(|_| "stream did not contain valid UTF-8")?;
        if index == 0 && !text.starts_with('>') {
            return Err("Expected > at record start.".into());
        }
        if let Some(header) = text.strip_prefix('>') {
            if let Some((record, has_description)) = current.take() {
                let mut record = record;
                record.end_offset = offset as u64;
                if record.id.is_empty() && !has_description && record.residues.is_empty() {
                    return Ok(records);
                }
                records.push(record);
            }
            let mut fields = header.trim_end().splitn(2, char::is_whitespace);
            let id = fields.next().unwrap_or("").to_string();
            let has_description = fields.next().is_some();
            current = Some((
                Expected {
                    id,
                    header_offset: offset as u64,
                    sequence_offset: (offset + line.len()) as u64,
                    end_offset: 0,
                    residue_offsets: Vec::new(),
                    residues: Vec::new(),
                },
                has_description,
            ));
        } else if let Some((record, _)) = current.as_mut() {
            let kept = text.trim_end().len();
            for (at, &byte) in line[..kept].iter().enumerate() {
                record.residue_offsets.push((offset + at) as u64);
                record.residues.push(byte);
            }
        }
    }
    if let Some((mut record, has_description)) = current.take() {
        record.end_offset = input.len() as u64;
        if !(record.id.is_empty() && !has_description && record.residues.is_empty()) {
            records.push(record);
        }
    }
    Ok(records)
}

fn bio_records(input: &[u8]) -> Result<Vec<fasta::Record>, String> {
    fasta::Reader::new(input)
        .records()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| error.to_string())
}

fn residue_offset(record: &ScanRecord, input: &[u8], residue: u64) -> u64 {
    match &record.layout {
        LineLayout::Uniform { width, eol } => {
            record.sequence_offset + residue / width * (width + eol) + residue % width
        }
        LineLayout::Checkpoints { offsets } => {
            // The application's forward read: from the checkpoint, count the bytes that
            // bio keeps (everything but each line's trailing whitespace).
            let checkpoint = residue / CHECKPOINT_EVERY;
            let mut offset = offsets[checkpoint as usize] as usize;
            let mut remaining = residue - checkpoint * CHECKPOINT_EVERY;
            loop {
                // The kept bytes end where the trimmed line ends; a checkpoint may fall
                // inside a character, so the line is trimmed from its start.
                let line_start = input[..offset]
                    .iter()
                    .rposition(|&byte| byte == b'\n')
                    .map_or(0, |at| at + 1);
                let line_end = input[offset..]
                    .iter()
                    .position(|&byte| byte == b'\n')
                    .map_or(input.len(), |at| offset + at);
                let kept_end = line_start
                    + std::str::from_utf8(&input[line_start..line_end])
                        .expect("valid line")
                        .trim_end()
                        .len();
                let kept = kept_end.saturating_sub(offset) as u64;
                if remaining < kept {
                    return (offset as u64) + remaining;
                }
                remaining -= kept;
                offset = line_end + 1;
            }
        }
    }
}

/// What the generated inputs exercised.
#[derive(Default, Debug)]
struct Coverage {
    errors: usize,
    uniform: usize,
    checkpoints: usize,
    large: usize,
    crlf: usize,
    non_ascii: usize,
}

fn check(input: &[u8], random: &mut Random, context: &str, coverage: &mut Coverage) {
    if input.windows(2).any(|pair| pair == b"\r\n") {
        coverage.crlf += 1;
    }
    if !input.is_ascii() {
        coverage.non_ascii += 1;
    }
    let expected = reference(input);
    let parsed = bio_records(input);
    let scanned = scan(input);
    assert_eq!(
        expected.is_err(),
        parsed.is_err(),
        "{context}: reference vs bio"
    );
    match (&expected, &scanned) {
        (Err(error), Err(scan_error)) => {
            assert_eq!(error, scan_error, "{context}: error");
            coverage.errors += 1;
            return;
        }
        (Ok(_), Ok(_)) => {}
        _ => panic!("{context}: reference {expected:?} vs scan {scanned:?}"),
    }
    let (expected, parsed, scanned) = (expected.unwrap(), parsed.unwrap(), scanned.unwrap());
    assert_eq!(parsed.len(), expected.len(), "{context}: bio record count");
    assert_eq!(
        scanned.len(),
        expected.len(),
        "{context}: scan record count"
    );
    for ((record, bio), want) in scanned.iter().zip(&parsed).zip(&expected) {
        assert_eq!(record.id, bio.id(), "{context}: id");
        assert_eq!(
            bio.seq(),
            want.residues.as_slice(),
            "{context}: reference residues"
        );
        assert_eq!(record.length, bio.seq().len() as u64, "{context}: length");
        assert_eq!(
            record.header_offset, want.header_offset,
            "{context}: header offset"
        );
        assert_eq!(
            record.sequence_offset, want.sequence_offset,
            "{context}: sequence offset"
        );
        assert_eq!(record.end_offset, want.end_offset, "{context}: end offset");
        let mut counts = [0u64; 256];
        for &byte in bio.seq() {
            counts[byte as usize] += 1;
        }
        assert_eq!(*record.residue_counts, counts, "{context}: residue counts");
        match &record.layout {
            LineLayout::Uniform { width, .. } if *width > 0 => coverage.uniform += 1,
            LineLayout::Checkpoints { .. } => coverage.checkpoints += 1,
            _ => {}
        }
        if record.length > 2 * CHECKPOINT_EVERY {
            coverage.large += 1;
        }
        if let LineLayout::Checkpoints { offsets } = &record.layout {
            let expected_checkpoints: Vec<u64> = want
                .residue_offsets
                .iter()
                .step_by(CHECKPOINT_EVERY as usize)
                .copied()
                .collect();
            assert_eq!(offsets, &expected_checkpoints, "{context}: checkpoints");
        }
        // Every residue of small records, a sample of large ones.
        let length = record.length;
        let samples: Vec<u64> = if length <= 4096 {
            (0..length).collect()
        } else {
            (0..512)
                .map(|_| random.below(length))
                .chain([0, length - 1])
                .collect()
        };
        for residue in samples {
            assert_eq!(
                residue_offset(record, input, residue),
                want.residue_offsets[residue as usize],
                "{context}: offset of residue {residue} ({:?})",
                record.layout
            );
        }
    }
    // Chunked scans give the same records.
    for chunking in 0..3 {
        let mut scanner = Scanner::new();
        let mut at = 0;
        while at < input.len() {
            let size = match chunking {
                0 => 1,
                1 => 1 + random.below(7) as usize,
                _ => 1 + random.below(4096) as usize,
            };
            let end = (at + size).min(input.len());
            scanner.feed(&input[at..end]);
            at = end;
        }
        assert_eq!(
            scanner.finish().as_ref(),
            Ok(&scanned),
            "{context}: chunking {chunking}"
        );
    }
}

const IDS: [&str; 8] = [
    "seq1",
    "chr_2",
    "gi|123|ref|NC_1.1|",
    "é-query",
    "日本",
    "a",
    "",
    "x\u{a0}y",
];
const DESCRIPTIONS: [&str; 6] = ["", " desc", " two words", "\tTAB", " ümlaut ✓", "  "];
const RESIDUES: &[u8] = b"ACGTNacgtnRYKMSWBDHV*-XU";
const WHITESPACE: [&str; 6] = [" ", "\t", "  ", "\u{a0}", "\u{3000}", "\r"];

fn generate(random: &mut Random) -> Vec<u8> {
    let crlf = random.chance(30);
    let eol = if crlf { "\r\n" } else { "\n" };
    let mut text = String::new();
    let records = random.below(6);
    for _ in 0..records {
        text.push('>');
        text.push_str(random.pick::<&str>(&IDS));
        text.push_str(random.pick::<&str>(&DESCRIPTIONS));
        text.push_str(eol);
        let length = match random.below(10) {
            0 => 0,
            1 => CHECKPOINT_EVERY * 2 + random.below(CHECKPOINT_EVERY),
            _ => random.below(400),
        };
        let uniform = random.chance(50);
        let width = 1 + random.below(80);
        let mut written = 0;
        if random.chance(10) {
            text.push_str(eol); // a blank line before the sequence
        }
        while written < length {
            let line = if uniform { width } else { 1 + random.below(90) }.min(length - written);
            for position in 0..line {
                if !uniform && position > 0 && random.chance(2) {
                    text.push_str(random.pick::<&str>(&WHITESPACE)); // interior whitespace
                }
                text.push(*random.pick(RESIDUES) as char);
            }
            written += line;
            if !uniform && random.chance(10) {
                text.push_str(random.pick::<&str>(&WHITESPACE)); // trailing whitespace
            }
            text.push_str(if !uniform && random.chance(10) {
                "\r\n"
            } else {
                eol
            });
            if !uniform && random.chance(5) {
                text.push_str(eol); // a blank line inside the sequence
            }
        }
        if random.chance(10) {
            text.push_str(eol); // a blank line after the sequence
        }
    }
    let mut bytes = text.into_bytes();
    if random.chance(10) && !bytes.is_empty() {
        bytes.pop(); // no final newline, or a cut CRLF
    }
    if random.chance(3) {
        let at = random.below(bytes.len() as u64 + 1) as usize;
        bytes.insert(at, 0xff); // invalid UTF-8
    }
    if random.chance(3) {
        bytes.insert(0, b'\n'); // not a record start
    }
    bytes
}

#[test]
fn scan_agrees_with_the_parser_on_generated_inputs() {
    let mut random = Random(0x5eed_1234_abcd_0001);
    let mut coverage = Coverage::default();
    for case in 0..2000 {
        let input = generate(&mut random);
        check(&input, &mut random, &format!("case {case}"), &mut coverage);
    }
    // The generator reaches every kind of input the scan distinguishes.
    eprintln!("{coverage:?}");
    let Coverage {
        errors,
        uniform,
        checkpoints,
        large,
        crlf,
        non_ascii,
    } = coverage;
    assert!(
        errors >= 50
            && uniform >= 500
            && checkpoints >= 500
            && large >= 50
            && crlf >= 300
            && non_ascii >= 500,
        "{coverage:?}"
    );
}

#[test]
fn scan_agrees_with_the_parser_on_edge_cases() {
    let mut random = Random(7);
    let cases: [&[u8]; 17] = [
        b"",
        b">",
        b">\n",
        b">\n>b\nAC\n",
        b"> x\nAC\n",
        b">a\n",
        b">a",
        b">a\r\nACGT\r\nAC\r\n",
        b">a\nACGT\nACGT\nAC",
        b">a\nACGT\n\nACGT\n",
        b">a\nAC GT  \nAC\t\n",
        b">a\nACGT\nACGTA\n",
        b"\n>a\nAC\n",
        b"AC\n>a\n",
        b">a\n>b\n>c\nA\n",
        b">\xc3\xa9t\xc3\xa9 d\n\xe3\x80\x80AC\xe3\x80\x80\n",
        b">a\nAC\n>\n>b\xff\n",
    ];
    for (index, input) in cases.iter().enumerate() {
        check(
            input,
            &mut random,
            &format!("edge {index}"),
            &mut Coverage::default(),
        );
    }
}
