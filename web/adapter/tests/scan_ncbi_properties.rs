//! Property tests of the index scan's parser kinds 1 (nucleotide flags) and 2 (protein
//! flags) against the engine's reader (`FastaInputSource::from_bytes`, NCBI BLAST+'s
//! `CFastaReader` with data loaders on): the records, IDs, lengths, residue counts and
//! reader errors agree with the reader; the offsets and the line layout agree with a
//! reference walk that records where every residue is (a pull-based model of NCBI's line
//! reader over the same bytes, itself checked against the reader's records); a record's
//! residues are found again from its layout (the uniform formula, or the read-forward
//! rules from every checkpoint); and scanning in chunks of any size gives the same result
//! as scanning at once.
//!
//! Inputs cover LF, CRLF, CR only and mixed line ends (a CR inside an LF file, an LF
//! inside a CR file, the tail that a file loses after a second pushback), no final line
//! end, comment and blank lines, a first record without a defline, `>?` gap lines,
//! `-`, `;`, interior white space, invalid residues, non-ASCII and non-UTF-8 bytes, empty
//! records, a BOM, Seq-id first lines, implausible first data lines, and records longer
//! than the checkpoint spacing.

use losat_web_adapter::scan::{scan_ncbi, LineLayout, NcbiScanner, ScanRecord, CHECKPOINT_EVERY};
use LOSAT::blastinput::fasta_reader::{
    FastaInputSource, FastaRecord, ParseErrorCode, ReadError, ReaderConfig,
};

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
    fn bytes(&mut self, items: &[&'static [u8]]) -> &'static [u8] {
        items[self.below(items.len() as u64) as usize]
    }
}

// ---------------------------------------------------------------------------------------
// The engine's reader.

#[derive(Debug)]
enum EngineEnd {
    Done,
    Parse { message: String, line: u64 },
    Unsupported(String),
}

fn engine(input: &[u8], protein: bool) -> (Vec<FastaRecord>, EngineEnd) {
    let mut source =
        FastaInputSource::from_bytes(input, ReaderConfig::subject("SCAN", protein, true));
    let mut records = Vec::new();
    loop {
        if source.end() {
            return (records, EngineEnd::Done);
        }
        match source.next_sequence(&mut |_| Ok(())) {
            Ok(record) => records.push(record),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => return (records, EngineEnd::Done),
            Err(ReadError::Parse { message, line, .. }) => {
                return (records, EngineEnd::Parse { message, line })
            }
            Err(ReadError::Unsupported(error)) => {
                return (records, EngineEnd::Unsupported(format!("{error:#}")))
            }
            Err(ReadError::Write(error)) => panic!("{error}"),
        }
    }
}

// ---------------------------------------------------------------------------------------
// The reference walk: NCBI's line reader (CStreamLineReader over the pushback streambuf of
// an ifstream, line_reader.cpp:154-281, stream_utils.cpp:224-443) and record reader
// (fasta.cpp:312-440), one byte at a time over offsets into the input.

const WINDOW: usize = 8191;
const MIN_BUF: usize = 4096;

struct Pushback {
    offsets: Vec<usize>,
    position: usize,
    capacity: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Eol {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}

struct Lines<'a> {
    input: &'a [u8],
    /// The file's window: the bytes `cursor..filled` are buffered; `filled` is the
    /// position of the file.
    cursor: usize,
    filled: usize,
    pushback: Vec<Pushback>,
    failed: bool,
    eol: Eol,
    line: Vec<usize>,
    line_number: u64,
    unget: bool,
    /// Bytes were pushed back (NCBI joins lines).
    pushed: bool,
    /// Pushed-back bytes stepped back into the buffer at the end of the input.
    lost: bool,
}

impl<'a> Lines<'a> {
    fn new(input: &'a [u8]) -> Self {
        Self {
            input,
            cursor: 0,
            filled: 0,
            pushback: Vec::new(),
            failed: false,
            eol: Eol::Unknown,
            line: Vec::new(),
            line_number: 0,
            unget: false,
            pushed: false,
            lost: false,
        }
    }

    fn fill(&mut self) {
        if self.cursor == self.filled {
            self.filled = (self.filled + WINDOW).min(self.input.len());
        }
    }

    fn raw_peek(&mut self) -> Option<usize> {
        loop {
            let Some(top) = self.pushback.last() else {
                self.fill();
                return (self.cursor < self.filled).then_some(self.cursor);
            };
            if top.position < top.offsets.len() {
                return Some(top.offsets[top.position]);
            }
            if self.pushback.len() > 1 {
                let below = self.pushback.remove(self.pushback.len() - 2);
                if below.position < below.offsets.len() {
                    *self.pushback.last_mut().unwrap() = below;
                }
                continue;
            }
            // x_FillBuffer(in_avail()): at most the buffer's size (at least 4096).
            let buffered = self.filled - self.cursor;
            let available = if buffered != 0 {
                buffered
            } else {
                self.input.len() - self.filled
            };
            let capacity = top.capacity.max(MIN_BUF);
            let requested = capacity.min(available.max(1));
            let mut offsets = Vec::new();
            while offsets.len() < requested {
                let want = requested - offsets.len();
                if self.cursor == self.filled && want >= WINDOW {
                    // A read of a whole window or more bypasses the empty window.
                    let n = want.min(self.input.len() - self.filled);
                    if n == 0 {
                        break;
                    }
                    offsets.extend(self.filled..self.filled + n);
                    self.filled += n;
                    self.cursor = self.filled;
                } else {
                    self.fill();
                    let n = want.min(self.filled - self.cursor);
                    if n == 0 {
                        break;
                    }
                    offsets.extend(self.cursor..self.cursor + n);
                    self.cursor += n;
                }
            }
            if offsets.is_empty() {
                return None;
            }
            let top = self.pushback.last_mut().unwrap();
            top.capacity = capacity;
            top.offsets = offsets;
            top.position = 0;
        }
    }

    fn raw_take(&mut self) -> Option<usize> {
        let offset = self.raw_peek()?;
        match self.pushback.last_mut() {
            Some(top) => top.position += 1,
            None => self.cursor += 1,
        }
        Some(offset)
    }

    fn push_back(&mut self, offsets: &[usize]) {
        let mut remaining = offsets.len();
        if let Some(top) = self.pushback.last_mut() {
            if offsets.len() <= MIN_BUF >> 4 {
                let take = top.position.min(remaining);
                top.position -= take;
                remaining -= take;
                top.offsets[top.position..top.position + take]
                    .copy_from_slice(&offsets[remaining..]);
            }
        }
        if remaining != 0 {
            self.pushback.push(Pushback {
                offsets: offsets[..remaining].to_vec(),
                position: 0,
                capacity: remaining,
            });
            self.failed = false;
        } else if !offsets.is_empty() && self.failed {
            self.lost = true;
        }
        self.pushed |= !offsets.is_empty();
    }

    fn stream_at_eof(&mut self) -> bool {
        if !self.failed && self.raw_peek().is_none() {
            self.failed = true;
        }
        self.failed
    }

    fn next_byte(&mut self) -> Option<usize> {
        if self.failed {
            return None;
        }
        let offset = self.raw_take();
        if offset.is_none() {
            self.failed = true;
        }
        offset
    }

    fn peek(&mut self) -> Option<usize> {
        if self.failed {
            return None;
        }
        let offset = self.raw_peek();
        if offset.is_none() {
            self.failed = true;
        }
        offset
    }

    fn byte(&self, offset: usize) -> u8 {
        self.input[offset]
    }

    /// NcbiGetline with one or two delimiters; the last byte read.
    fn getline(&mut self, delimiters: &[u8]) -> Option<u8> {
        self.line.clear();
        loop {
            let offset = if delimiters.len() == 1 {
                self.next_byte()
            } else {
                self.raw_take()
            };
            let Some(offset) = offset else {
                self.failed = true;
                return self.line.last().map(|&offset| self.byte(offset));
            };
            let byte = self.byte(offset);
            if let Some(position) = delimiters.iter().position(|&d| d == byte) {
                let mut last = Some(byte);
                if delimiters.len() > 1 {
                    if let Some(next) = self.raw_peek() {
                        if delimiters[position + 1..].contains(&self.byte(next)) {
                            last = self.raw_take().map(|offset| self.byte(offset));
                        }
                    }
                }
                return last;
            }
            self.line.push(offset);
        }
    }

    fn advance_simple(&mut self, eol: u8, alternate: u8) -> Eol {
        self.line.clear();
        let mut watched = None;
        loop {
            let Some(offset) = self.next_byte() else {
                self.failed = true;
                break;
            };
            let byte = self.byte(offset);
            if byte == eol {
                break;
            }
            if byte == alternate && watched.is_none() {
                watched = Some(self.line.len());
            }
            self.line.push(offset);
        }
        if let Some(position) = watched {
            let position = position + 1;
            if eol != b'\n' || position != self.line.len() {
                let rest = self.line[position..].to_vec();
                self.push_back(&rest);
                self.eol = Eol::Mixed;
            }
            self.line.truncate(position - 1);
            return if self.eol == Eol::Mixed {
                Eol::Mixed
            } else {
                Eol::CrLf
            };
        }
        if eol == b'\r' && !self.failed {
            if let Some(next) = self.raw_peek() {
                if self.byte(next) == alternate {
                    self.raw_take();
                    return Eol::CrLf;
                }
            }
        }
        if eol == b'\r' {
            Eol::Cr
        } else {
            Eol::Lf
        }
    }

    fn advance(&mut self) {
        match self.eol {
            Eol::Unknown => {
                let last = self.getline(b"\r\n");
                if self.failed && !self.line.is_empty() {
                    self.failed = false;
                }
                match last {
                    Some(b'\r') => self.eol = Eol::Cr,
                    Some(b'\n') => self.eol = Eol::CrLf,
                    _ => {}
                }
            }
            Eol::Cr => {
                self.advance_simple(b'\r', b'\n');
            }
            Eol::Lf => {
                self.advance_simple(b'\n', b'\r');
            }
            Eol::CrLf => {
                let style = self.advance_simple(b'\n', b'\r');
                if style == Eol::Mixed {
                    self.eol = Eol::Cr;
                } else if style != Eol::CrLf {
                    self.eol = Eol::Lf;
                }
            }
            Eol::Mixed => {
                self.getline(b"\r\n");
            }
        }
    }

    fn at_eof(&mut self) -> bool {
        !self.unget && self.stream_at_eof()
    }

    fn peek_char(&mut self) -> Option<u8> {
        if self.at_eof() {
            return self.peek().map(|offset| self.byte(offset));
        }
        if self.unget {
            return Some(self.line.first().map_or(0, |&offset| self.byte(offset)));
        }
        match self.peek().map(|offset| self.byte(offset)) {
            Some(b'\n' | b'\r') => Some(0),
            other => other,
        }
    }

    fn unget_line(&mut self) {
        if !self.unget && self.line_number != 0 {
            self.line_number -= 1;
            self.unget = true;
        }
    }

    fn next_line(&mut self) -> Vec<usize> {
        if self.at_eof() {
            self.line.clear();
            return Vec::new();
        }
        self.line_number += 1;
        if self.unget {
            self.unget = false;
        } else {
            self.advance();
        }
        self.line.clone()
    }
}

fn space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

fn residue_type(byte: u8, protein: bool) -> bool {
    let upper = byte.to_ascii_uppercase();
    if protein {
        byte.is_ascii_alphabetic() || byte == b'*'
    } else {
        b"ABCDGHKMRSTUVWYN".contains(&upper)
    }
}

fn stored(byte: u8, protein: bool) -> u8 {
    let upper = byte.to_ascii_uppercase();
    if !protein && upper == b'U' {
        b'T'
    } else {
        upper
    }
}

#[derive(Debug, Clone)]
struct RefRecord {
    title: Vec<u8>,
    header_offset: usize,
    sequence_offset: usize,
    end_offset: usize,
    /// Stored residues (upper-cased) and where they are (none for a gap's residues).
    residues: Vec<(u8, Option<usize>)>,
    /// The number of the line that ended the record (`u64::MAX` at the end of the input).
    close_line: u64,
}

#[derive(Debug, Default)]
struct Reference {
    records: Vec<RefRecord>,
    /// The first gap line.
    gap_line: Option<u64>,
    joins: bool,
    lost: bool,
}

fn trim(line: &[usize], input: &[u8]) -> Vec<usize> {
    let start = line
        .iter()
        .position(|&o| !space(input[o]))
        .unwrap_or(line.len());
    let end = line
        .iter()
        .rposition(|&o| !space(input[o]))
        .map_or(start, |p| p + 1);
    line[start..end.max(start)].to_vec()
}

/// `>?_` at the start becomes `>` (its offset is the `>`'s).
fn modify(line: Vec<usize>, input: &[u8]) -> Vec<usize> {
    if line.len() >= 3
        && line[..3]
            .iter()
            .map(|&o| input[o])
            .eq(b">?_".iter().copied())
    {
        let mut modified = vec![line[0]];
        modified.extend_from_slice(&line[3..]);
        modified
    } else {
        line
    }
}

fn starts_with(line: &[usize], input: &[u8], prefix: &[u8]) -> bool {
    line.len() >= prefix.len()
        && line[..prefix.len()]
            .iter()
            .map(|&o| input[o])
            .eq(prefix.iter().copied())
}

/// ParseGapLine's length (fasta.cpp:1094-1349): leading digits after `>?` (and `unk`).
fn gap_length(line: &[u8]) -> usize {
    let mut rest: &[u8] = &line[2..];
    while rest.first().is_some_and(|&b| space(b)) {
        rest = &rest[1..];
    }
    if let Some(after) = rest.strip_prefix(b"unk") {
        rest = after;
        while rest.first().is_some_and(|&b| space(b)) {
            rest = &rest[1..];
        }
    }
    let digits = rest.iter().take_while(|b| b.is_ascii_digit()).count();
    std::str::from_utf8(&rest[..digits])
        .ok()
        .and_then(|text| text.parse::<u32>().ok())
        .filter(|&size| size > 0)
        .map_or(1, |size| size as usize)
}

/// The title of a defline (ParseDefline with fNoParseID, then x_ApplyMods' trim).
fn title(defline: &[u8]) -> Vec<u8> {
    let len = defline.len();
    if len <= 1 || defline[1..].iter().all(|&b| space(b)) {
        return Vec::new();
    }
    let start = 1 + defline[1..].iter().position(|&b| !space(b)).unwrap();
    let end = start
        + 1
        + defline[start + 1..]
            .iter()
            .position(|&b| b < b' ')
            .unwrap_or(len - start - 1);
    let text = &defline[start..end];
    let last = text.iter().rposition(|&b| !space(b)).map_or(0, |p| p + 1);
    text[..last].to_vec()
}

fn reference(input: &[u8], protein: bool) -> Reference {
    let mut lines = Lines::new(input);
    let mut out = Reference::default();
    let total = input.len();
    let bytes = |line: &[usize]| line.iter().map(|&o| input[o]).collect::<Vec<u8>>();
    while !lines.at_eof() {
        // ReadOneSeq of CBlastInputReader: the first line is tried as a Seq-id.
        lines.next_line();
        lines.unget_line();
        let mut need_defline = true;
        let mut record = RefRecord {
            title: Vec::new(),
            header_offset: 0,
            sequence_offset: 0,
            end_offset: total,
            residues: Vec::new(),
            close_line: u64::MAX,
        };
        while !lines.at_eof() {
            let c = lines.peek_char();
            assert!(!lines.at_eof(), "unexpected end of file");
            if c == Some(b'>') {
                let raw = lines.next_line();
                let modified = modify(raw.clone(), input);
                if starts_with(&modified, input, b">?") {
                    lines.unget_line();
                } else if need_defline {
                    record.title = title(&bytes(&modified));
                    record.header_offset = raw[0];
                    record.sequence_offset = lines.peek().unwrap_or(total);
                    need_defline = false;
                    continue;
                } else {
                    record.end_offset = raw[0];
                    record.close_line = lines.line_number;
                    lines.unget_line();
                    break;
                }
            }
            let raw = lines.next_line();
            let line = trim(&raw, input);
            if line.is_empty() || matches!(input[line[0]], b'!' | b'#' | b';') {
                continue;
            }
            if need_defline {
                need_defline = false;
            }
            let line = modify(line, input);
            if starts_with(&line, input, b">?") {
                out.gap_line.get_or_insert(lines.line_number);
                let gap = if protein { b'X' } else { b'N' };
                for _ in 0..gap_length(&bytes(&line)) {
                    record.residues.push((gap, None));
                }
                continue;
            }
            for &offset in &line {
                let byte = input[offset];
                if byte == b';' {
                    break;
                }
                if residue_type(byte, protein) {
                    record.residues.push((stored(byte, protein), Some(offset)));
                }
            }
        }
        if need_defline {
            break;
        }
        out.records.push(record);
    }
    out.joins = lines.pushed;
    out.lost = lines.lost;
    out
}

// ---------------------------------------------------------------------------------------
// The layout as the application uses it.

/// The uniform layout of residue offsets, if they have one.
fn uniform(offsets: &[usize], sequence_offset: usize) -> Option<(u64, u64)> {
    if offsets.is_empty() {
        return Some((0, 1));
    }
    if offsets[0] != sequence_offset {
        return None;
    }
    let Some(width) = (1..offsets.len()).find(|&i| offsets[i] != offsets[i - 1] + 1) else {
        return Some((offsets.len() as u64, 1));
    };
    let eol = offsets[width] - offsets[width - 1] - 1;
    offsets
        .iter()
        .enumerate()
        .all(|(i, &o)| o == sequence_offset + i / width * (width + eol) + i % width)
        .then_some((width as u64, eol as u64))
}

/// The read-forward rules (docs/web/abi_v2.md §9, kinds 1 and 2) from a checkpoint:
/// the offsets of at most `count` residues, the checkpoint's first.
fn read_forward(input: &[u8], from: usize, count: usize, protein: bool) -> Vec<usize> {
    #[derive(PartialEq)]
    enum State {
        LineStart,
        Data,
        Skip,
    }
    let mut state = State::Data;
    let mut found = Vec::new();
    let mut offset = from;
    while found.len() < count && offset < input.len() {
        let byte = input[offset];
        if byte == b'\r' || byte == b'\n' {
            state = State::LineStart;
        } else {
            if state == State::LineStart && !matches!(byte, b' ' | b'\t' | 0x0b | 0x0c) {
                state = if matches!(byte, b'!' | b'#' | b';') {
                    State::Skip
                } else {
                    State::Data
                };
            }
            if state == State::Data {
                if byte == b';' {
                    state = State::Skip;
                } else if residue_type(byte, protein) {
                    found.push(offset);
                }
            }
        }
        offset += 1;
    }
    found
}

/// Whether the read-forward rules find the residues from every checkpoint.
fn readable(input: &[u8], offsets: &[usize], protein: bool) -> bool {
    let every = CHECKPOINT_EVERY as usize;
    (0..offsets.len()).step_by(every).all(|start| {
        let count = every.min(offsets.len() - start);
        read_forward(input, offsets[start], count, protein) == offsets[start..start + count]
    })
}

// ---------------------------------------------------------------------------------------
// The expected scan.

#[derive(Debug)]
enum Expected {
    Records(Vec<(RefRecord, LineLayout)>),
    SeqId,
    Gap(u64),
    Parse(String),
    NotReadable(usize),
}

/// What the generated inputs exercised.
#[derive(Default, Debug)]
struct Coverage {
    records: usize,
    uniform: usize,
    checkpoints: usize,
    large: usize,
    empty: usize,
    headerless: usize,
    seq_id: usize,
    gap: usize,
    parse: usize,
    not_readable: usize,
    joins: usize,
    lost: usize,
    crlf: usize,
    cr_only: usize,
    non_utf8: usize,
}

fn expected(input: &[u8], protein: bool, context: &str, coverage: &mut Coverage) -> Expected {
    let (records, end) = engine(input, protein);
    let reference = reference(input, protein);
    // The reference reads the reader's records (titles and residues).
    for (index, record) in records.iter().enumerate() {
        let want = &reference.records[index];
        assert_eq!(record.title, want.title, "{context}: title {index}");
        let residues: Vec<u8> = want.residues.iter().map(|&(r, _)| r).collect();
        assert_eq!(
            record.sequence.to_ascii_uppercase(),
            residues,
            "{context}: residues {index}"
        );
    }
    if matches!(end, EngineEnd::Done) {
        assert_eq!(records.len(), reference.records.len(), "{context}: count");
    }
    let stop_line = match &end {
        EngineEnd::Unsupported(message) => {
            assert!(
                message.contains("may be a sequence identifier"),
                "{context}: {message}"
            );
            coverage.seq_id += 1;
            return Expected::SeqId;
        }
        EngineEnd::Parse { line, .. } => *line,
        EngineEnd::Done => u64::MAX,
    };
    // The first rejection, by line: a gap line, or a record that the layout cannot locate
    // (rejected when the line that ends it starts).
    let mut rejection: Option<(u64, Expected)> = reference
        .gap_line
        .filter(|&line| line < stop_line)
        .map(|line| (line, Expected::Gap(line)));
    coverage.joins += reference.joins as usize;
    coverage.lost += reference.lost as usize;
    let mut layouts = Vec::new();
    for (index, record) in reference.records.iter().enumerate() {
        // Records that close before the reader's error (or at the end of the input).
        let closes = record.close_line < stop_line
            || (stop_line == u64::MAX && record.close_line == u64::MAX);
        if !closes {
            break;
        }
        if rejection
            .as_ref()
            .is_some_and(|(line, _)| *line <= record.close_line)
        {
            break;
        }
        if record.residues.iter().any(|(_, offset)| offset.is_none()) {
            break; // a gap: rejected at the gap line
        }
        let offsets: Vec<usize> = record.residues.iter().map(|&(_, o)| o.unwrap()).collect();
        let layout = match uniform(&offsets, record.sequence_offset) {
            Some((width, eol)) => LineLayout::Uniform { width, eol },
            None if readable(input, &offsets, protein) => LineLayout::Checkpoints {
                offsets: offsets
                    .iter()
                    .step_by(CHECKPOINT_EVERY as usize)
                    .map(|&o| o as u64)
                    .collect(),
            },
            None => {
                rejection = Some((record.close_line, Expected::NotReadable(index)));
                break;
            }
        };
        layouts.push((record.clone(), layout));
        if record.close_line == u64::MAX {
            break;
        }
    }
    if let Some((_, rejection)) = rejection {
        match rejection {
            Expected::Gap(_) => coverage.gap += 1,
            _ => coverage.not_readable += 1,
        }
        return rejection;
    }
    if let EngineEnd::Parse { message, .. } = end {
        coverage.parse += 1;
        return Expected::Parse(message);
    }
    Expected::Records(layouts)
}

fn check(input: &[u8], random: &mut Random, context: &str, coverage: &mut Coverage) {
    for protein in [false, true] {
        let context = format!("{context} kind {}", 1 + protein as u8);
        let want = expected(input, protein, &context, coverage);
        let scanned = scan_ncbi(input, protein);
        match (&want, &scanned) {
            (Expected::Records(records), Ok(scanned)) => {
                assert_eq!(scanned.len(), records.len(), "{context}: record count");
                for (index, ((record, layout), got)) in records.iter().zip(scanned).enumerate() {
                    check_record(input, protein, record, layout, got, random, &context, index);
                    coverage.records += 1;
                    match layout {
                        LineLayout::Uniform { width, .. } if *width > 0 => coverage.uniform += 1,
                        LineLayout::Checkpoints { .. } => coverage.checkpoints += 1,
                        _ => {}
                    }
                    if got.length == 0 {
                        coverage.empty += 1;
                    }
                    if got.length > 2 * CHECKPOINT_EVERY {
                        coverage.large += 1;
                    }
                    if index == 0 && got.header_offset == got.sequence_offset {
                        coverage.headerless += 1;
                    }
                }
            }
            (Expected::SeqId, Err(error)) => assert!(
                error.contains("may be a sequence identifier")
                    && error.contains("not supported by LOSAT Web"),
                "{context}: {error}"
            ),
            (Expected::Gap(line), Err(error)) => assert!(
                error.starts_with(&format!("line {line} is a gap line"))
                    && error.contains("not supported by LOSAT Web"),
                "{context}: {error}"
            ),
            (Expected::Parse(message), Err(error)) => {
                assert_eq!(error, message, "{context}: reader error")
            }
            (Expected::NotReadable(index), Err(error)) => assert!(
                error.starts_with(&format!("record {} ", index + 1))
                    && error.contains("reads two of its lines as one")
                    && error.contains("not supported by LOSAT Web"),
                "{context}: {error}"
            ),
            _ => panic!("{context}: expected {want:?}, scanned {scanned:?}"),
        }
        // Chunked scans give the same result.
        for chunking in 0..3 {
            let mut scanner = NcbiScanner::new(protein);
            let mut at = 0;
            while at < input.len() {
                let size = match chunking {
                    0 => 1,
                    1 => 1 + random.below(7) as usize,
                    _ => 1 + random.below(9000) as usize,
                };
                let end = (at + size).min(input.len());
                scanner.feed(&input[at..end]);
                at = end;
            }
            assert_eq!(scanner.finish(), scanned, "{context}: chunking {chunking}");
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn check_record(
    input: &[u8],
    protein: bool,
    record: &RefRecord,
    layout: &LineLayout,
    got: &ScanRecord,
    random: &mut Random,
    context: &str,
    index: usize,
) {
    let first_word = record
        .title
        .split(|&b| b == b' ')
        .next()
        .unwrap_or_default();
    assert_eq!(
        got.id,
        String::from_utf8_lossy(first_word),
        "{context}: id {index}"
    );
    assert_eq!(
        got.header_offset, record.header_offset as u64,
        "{context}: header offset {index}"
    );
    assert_eq!(
        got.sequence_offset, record.sequence_offset as u64,
        "{context}: sequence offset {index}"
    );
    assert_eq!(
        got.end_offset, record.end_offset as u64,
        "{context}: end offset {index}"
    );
    assert_eq!(
        got.length,
        record.residues.len() as u64,
        "{context}: length {index}"
    );
    let mut counts = [0u64; 256];
    for &(residue, _) in &record.residues {
        counts[residue as usize] += 1;
    }
    assert_eq!(
        *got.residue_counts, counts,
        "{context}: residue counts {index}"
    );
    assert_eq!(&got.layout, layout, "{context}: layout {index}");
    // The layout locates the residues: the formula, or reading forward from checkpoints.
    let offsets: Vec<usize> = record.residues.iter().map(|&(_, o)| o.unwrap()).collect();
    let length = offsets.len();
    let samples: Vec<usize> = if length <= 4096 {
        (0..length).collect()
    } else {
        (0..512)
            .map(|_| random.below(length as u64) as usize)
            .chain([0, length - 1])
            .collect()
    };
    for residue in samples {
        let at = match &got.layout {
            LineLayout::Uniform { width, eol } => {
                let r = residue as u64;
                (got.sequence_offset + r / width * (width + eol) + r % width) as usize
            }
            LineLayout::Checkpoints { offsets: points } => {
                let every = CHECKPOINT_EVERY as usize;
                let from = points[residue / every] as usize;
                let read = read_forward(input, from, residue % every + 1, protein);
                read[residue % every]
            }
        };
        assert_eq!(
            at, offsets[residue],
            "{context}: residue {residue} of {index}"
        );
        assert_eq!(
            stored(input[at], protein),
            record.residues[residue].0,
            "{context}: residue byte"
        );
    }
}

// ---------------------------------------------------------------------------------------
// Inputs.

const TITLES: [&[u8]; 14] = [
    b"seq1",
    b"chr_2 two words",
    b"gi|123|ref|NC_1.1| desc",
    b" leading",
    b"\ttab\tinside",
    b"\xc3\xa9t\xc3\xa9 caf\xc3\xa9",
    b"bad\xff\xfeutf8 x",
    b"",
    b"   ",
    b"x\x01ctl",
    b"ACGTACGTACGTACGTACGTACGTACGT",
    b"?_odd",
    b"q;semi",
    b"\xef\xbb\xbfbom",
];
const NUCLEOTIDES: &[u8] = b"ACGTACGTACGTNacgtnRYKMSWBDHVuU";
const PROTEINS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwy*UOBZJ";
const NOISE: &[&[u8]] = &[
    b" ",
    b"\t",
    b"-",
    b"--",
    b"1",
    b"9",
    b".",
    b"@",
    b"\x0b",
    b"\x0c",
    b"\xc3\xa9",
    b"\xff",
    b"*",
    b"e",
    b"E",
    b"x",
    b"X",
    b"\xa0",
    b">",
    b"!",
    b"#",
];
const COMMENTS: &[&[u8]] = &[
    b"# comment",
    b"!bang",
    b";semi colon",
    b"  # indented",
    b"\t;tab",
    b"",
    b"   ",
];
const FIRST_LINES: &[&[u8]] = &[
    b"AB123456",
    b"gb|X|",
    b"NC_000001.1",
    b"ACGTACGTAC",
    b"contig1",
    b"0123",
    b"1abc",
    b"\xef\xbb\xbfACGT",
];

#[derive(Clone, Copy)]
enum Style {
    Lf,
    CrLf,
    Cr,
    Mixed,
}

fn eol(style: Style, random: &mut Random) -> &'static [u8] {
    match style {
        Style::Lf => b"\n",
        Style::CrLf => b"\r\n",
        Style::Cr => b"\r",
        Style::Mixed => random.bytes(&[b"\n", b"\r\n", b"\r", b"\n\r", b"\r\r\n"]),
    }
}

fn generate(random: &mut Random) -> Vec<u8> {
    let style = *random.pick(&[Style::Lf, Style::Lf, Style::CrLf, Style::Cr, Style::Mixed]);
    let protein_letters = random.chance(40);
    let letters = if protein_letters {
        PROTEINS
    } else {
        NUCLEOTIDES
    };
    let noisy = random.chance(40);
    let mut text: Vec<u8> = Vec::new();
    if random.chance(3) {
        text.extend_from_slice(b"\xef\xbb\xbf");
    }
    if random.chance(4) {
        text.extend_from_slice(random.bytes(FIRST_LINES));
        text.extend_from_slice(eol(style, random));
    }
    if random.chance(10) {
        text.extend_from_slice(random.bytes(COMMENTS));
        text.extend_from_slice(eol(style, random));
    }
    let records = random.below(6);
    let headerless = random.chance(10);
    // At most one gap line, in a random record.
    let gap_record = random.chance(8).then(|| random.below(records.max(1)));
    for record in 0..records {
        if !(record == 0 && headerless) {
            text.push(b'>');
            if random.chance(3) {
                text.extend_from_slice(b"?_");
            }
            text.extend_from_slice(random.bytes(&TITLES));
            text.extend_from_slice(eol(style, random));
        }
        let length = match random.below(12) {
            0 => 0,
            1 => CHECKPOINT_EVERY * 2 + random.below(CHECKPOINT_EVERY),
            _ => random.below(400),
        };
        let uniform = random.chance(50);
        let width = 1 + random.below(80);
        let mut written = 0;
        if random.chance(10) {
            text.extend_from_slice(eol(style, random)); // a blank line before the sequence
        }
        while written < length {
            let line = if uniform { width } else { 1 + random.below(90) }.min(length - written);
            for position in 0..line {
                if noisy && position > 0 && random.chance(3) {
                    text.extend_from_slice(random.bytes(NOISE));
                }
                text.push(*random.pick(letters));
            }
            written += line;
            if noisy && random.chance(5) {
                text.push(b';');
                text.extend_from_slice(b"ACGT note");
            }
            if !uniform && random.chance(3) {
                // A line end of another style inside the line (NCBI joins lines there),
                // sometimes before a comment mark or after a `;`.
                text.extend_from_slice(random.bytes(&[
                    b"\r", b"\n", b"\r#c", b"\n#c", b"\r!", b"\r;c", b";c\r", b"\r  #",
                ]));
            }
            text.extend_from_slice(eol(style, random));
            if !uniform && random.chance(4) {
                text.extend_from_slice(random.bytes(COMMENTS));
                text.extend_from_slice(eol(style, random));
            }
            if gap_record == Some(record) && random.chance(20) {
                // A gap line.
                text.extend_from_slice(random.bytes(&[
                    b">?5",
                    b">?",
                    b">?unk100",
                    b" >?3",
                    b">?_?2",
                ]));
                text.extend_from_slice(eol(style, random));
            }
        }
        if random.chance(10) {
            text.extend_from_slice(eol(style, random)); // a blank line after the sequence
        }
        if noisy && random.chance(5) {
            text.extend_from_slice(b"@@@@ junk");
            text.extend_from_slice(eol(style, random));
        }
    }
    if random.chance(15) && !text.is_empty() {
        text.pop(); // no final line end, or a cut CRLF
    }
    if random.chance(10) {
        // A tail after a late line end of the other style.
        text.extend_from_slice(random.bytes(&[
            b"\nACGT",
            b"\rACGT",
            b"\nAC-GT\n",
            b"\r\nGG",
            b"\rAC\nGGGG",
            b"\rAC\nGG\nTT",
        ]));
    }
    text
}

/// A file whose line 2 has a CR inside a CRLF-style line, whose last CR-style line has an
/// LF `offset` bytes before the end of 8191-byte window `window`, followed by `tail` bytes
/// and `end` (the engine's `window_boundary_input`: a tail of 2 to 256 bytes can be lost).
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

fn note(input: &[u8], coverage: &mut Coverage) {
    let crlf = input.windows(2).any(|pair| pair == b"\r\n");
    if crlf {
        coverage.crlf += 1;
    }
    if input.contains(&b'\r') && !input.contains(&b'\n') {
        coverage.cr_only += 1;
    }
    if std::str::from_utf8(input).is_err() {
        coverage.non_utf8 += 1;
    }
}

#[test]
fn ncbi_scan_agrees_with_the_reader_on_generated_inputs() {
    let mut random = Random(0x5eed_2026_1009_0009);
    let mut coverage = Coverage::default();
    let cases = std::env::var("SCAN_NCBI_CASES")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(1500);
    for case in 0..cases {
        let input = generate(&mut random);
        note(&input, &mut coverage);
        check(&input, &mut random, &format!("case {case}"), &mut coverage);
    }
    eprintln!("{coverage:?}");
    let Coverage {
        records,
        uniform,
        checkpoints,
        large,
        empty,
        headerless,
        seq_id,
        gap,
        parse,
        not_readable,
        joins,
        lost,
        crlf,
        cr_only,
        non_utf8,
    } = coverage;
    assert!(
        records >= 3000
            && uniform >= 500
            && checkpoints >= 500
            && large >= 50
            && empty >= 50
            && headerless >= 50
            && seq_id >= 20
            && gap >= 20
            && parse >= 50
            && not_readable >= 50
            && joins >= 300
            && lost >= 5
            && crlf >= 300
            && cr_only >= 100
            && non_utf8 >= 100,
        "{coverage:?}"
    );
}

#[test]
fn ncbi_scan_agrees_with_the_reader_on_edge_cases() {
    let mut random = Random(9);
    let mut coverage = Coverage::default();
    let mut cases: Vec<Vec<u8>> = [
        &b""[..],
        b"\n",
        b"   \n\t\n",
        b"# only a comment\n;another\n",
        b">",
        b">\n",
        b"> \n",
        b">\n>b\nAC\n",
        b">a\nACGT\nAC",
        b">a\r\nACGT\r\nAC\r\n",
        b">a\rACGT\rAC\r",
        b">a\nACGT\r\nAC\rGT\n",
        b">a\nAC GT  \n-AC;GT\nAC\t\n",
        b"ACGT\nACGT\n>b\nAC\n",
        b"# c\nACGT\n",
        b"\n\nACGT\n",
        b">a\n>?5\nACGT\n",
        b">a\nACGT\n>?\n",
        b">a\nACGT\n  >?3\n",
        b">a\nACGT\n>?_?2\n",
        b">?_x desc\nACGT\n",
        b">?_\nACGT\n",
        b">a\n >x\nACGT\n",
        b">a\n@@@@\n",
        b">a\nAC;@@@@@@@@\n",
        b">a\n---\n",
        b">a\n12345\n",
        b">a\nACGT\n@@@@\n",
        b"AB123456\nACGT\n",
        b"gb|X|\n",
        b"  NC_000001.1\n",
        b"contig1\nACGT\n",
        b"\xef\xbb\xbfACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        b"\xef\xbb\xbf>q1\nACGT\n",
        b">\xc3\xa9t\xc3\xa9 d\nACGT\n",
        b">bad\xff id\nACGT\n",
        b">a\nAC\xc3\xa9GT\xffAC\n",
        // Joins: a CR inside an LF file, an LF inside a CR file.
        b">q1 first\nACGTACGT\r\nCTTAAG\r>q2 second\nTGGCATTTTTATTACACTCAGAAACAG\n",
        b">f\rAC-GTACGTAGGC\rGG-ACACGTAGG\nATTAC-CAGACG\rTT-GGACGT\r",
        b">first\r\nTC-TTGGCTCAAT\rTT-TTTAACGTGAGG\r\r\nACGTACGT\n-ACGTACGTACGTAC",
        b">first\r\nTC-TTGGCTCAAT\rTT-TTTAACGTGAGG\r\r\nACGTACGT\n-ACGTACGTACGTAC\n",
        b">f\nAC\r#c\nACGT\nGG\n",
        b">f\nAC\rGT;x\nACGT\nGG\n",
        b">f\r\nA\r#cmt\nACGT\n",
        b">f\nAC\r>g\nGT\n",
        b">f\r\nAC\rGT\nTT\nAA",
        b">f\r\nAC\rGT\nTT\nAA\n",
        b">f\r\nAC\rGT\nTT\n>?5",
        b">f\r\nAC\rGT\nTT\n@@@@",
        b">f\nAC\r#c\nACGT\nGG\nG\n",
        b">f\nAC\rGT;x\nACGT\nGG\nG\n",
        b">f\r\nA\r#cmt\nACGT\nGG\nG\n",
        b">f\nAC\nGT\r!x\nACGT\nG\n",
        b">f\rAC\rGT\n#x\rACGT\rG\r",
    ]
    .iter()
    .map(|case| case.to_vec())
    .collect();
    for window in 1..=2 {
        for offset in 1..=4 {
            for tail in [0, 1, 2, 23, 255, 256, 257, 4100] {
                for end in [&b""[..], b"\r", b"\n", b"\r\n"] {
                    cases.push(window_boundary_input(window, offset, tail, end));
                }
            }
        }
    }
    // A record past the checkpoint spacing with irregular lines, and a uniform one.
    let mut long = b">long\n".to_vec();
    for line in 0..3000 {
        long.extend(std::iter::repeat_n(b'A', 40 + line % 50));
        long.extend_from_slice(if line % 7 == 0 { b"\r\n" } else { b"\n" });
    }
    cases.push(long);
    let mut long = b">long\r\n".to_vec();
    for _ in 0..3000 {
        long.extend_from_slice(&[b'C'; 60]);
        long.extend_from_slice(b"\r\n");
    }
    cases.push(long);
    for (index, input) in cases.iter().enumerate() {
        note(input, &mut coverage);
        check(input, &mut random, &format!("edge {index}"), &mut coverage);
    }
    eprintln!("{coverage:?}");
}
