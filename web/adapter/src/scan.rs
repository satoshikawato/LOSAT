//! The index scan of FASTA input (docs/web/abi_v2.md §9, plan TD-8). It exists only so
//! that the application can extract original records; every search reads its input
//! with the program's own parser, and `register` checks that the two agree.
//!
//! Parser kind 0 reproduces `bio::io::fasta::Reader` 1.6.0 (`src/io/fasta.rs:316-344`),
//! the reader of BLASTP, TBLASTN, BLASTN and TBLASTX:
//! - input is read line by line as UTF-8 (`read_line`); invalid UTF-8 is an error;
//! - the first line must start with `>` ("Expected > at record start.");
//! - the ID is the header text before its first whitespace character, after
//!   `trim_end`;
//! - every other line until the next line starting with `>` is a sequence line whose
//!   `trim_end` is appended, so interior whitespace is kept and trailing whitespace
//!   (including `\r`) is not;
//! - reading stops at the first empty record (no ID, no description, no sequence).
//!
//! The scan works on chunks of any size and never holds a sequence line, only a header
//! line and a run of whitespace whose fate (interior or trailing) is still open.

/// Residues between two checkpoints of a record whose lines are not regular.
pub const CHECKPOINT_EVERY: u64 = 65_536;

/// How to find the byte offset of a residue of a record in the scanned input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LineLayout {
    /// Residue `i` is at `sequence_offset + (i / width) * (width + eol) + i % width`.
    Uniform { width: u64, eol: u64 },
    /// `offsets[k]` is the byte offset of residue `k * CHECKPOINT_EVERY`; the residues
    /// in between are found by reading forward under the parser's rules.
    Checkpoints { offsets: Vec<u64> },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScanRecord {
    pub id: String,
    /// Offset of the `>` of the header line.
    pub header_offset: u64,
    /// Offset of the first byte after the header line.
    pub sequence_offset: u64,
    /// Offset of the next header line, or the end of the input.
    pub end_offset: u64,
    /// Sequence length in bytes, as the parser reports it.
    pub length: u64,
    pub layout: LineLayout,
    /// Count of every byte value in the sequence.
    pub residue_counts: Box<[u64; 256]>,
}

/// Per-record state while its sequence lines are read.
struct RecordBuilder {
    record: ScanRecord,
    /// The header had a description (the text after the ID).
    has_description: bool,
    /// Residues committed on the current line.
    line_residues: u64,
    /// The uniform layout candidate: (width, eol) of the first residue line.
    uniform: Option<(u64, u64)>,
    uniform_possible: bool,
    residue_lines: u64,
    /// A residue line shorter than the width was seen (it must be the last one).
    short_line_seen: bool,
    /// A line without residues followed a residue line.
    blank_after_residues: bool,
    checkpoints: Vec<u64>,
}

impl RecordBuilder {
    fn new(id: String, has_description: bool, header_offset: u64, sequence_offset: u64) -> Self {
        Self {
            record: ScanRecord {
                id,
                header_offset,
                sequence_offset,
                end_offset: sequence_offset,
                length: 0,
                layout: LineLayout::Uniform { width: 0, eol: 1 },
                residue_counts: Box::new([0; 256]),
            },
            has_description,
            line_residues: 0,
            uniform: None,
            uniform_possible: true,
            residue_lines: 0,
            short_line_seen: false,
            blank_after_residues: false,
            checkpoints: Vec::new(),
        }
    }

    /// Commits `bytes`, which start at input offset `offset`, as residues.
    fn commit(&mut self, offset: u64, bytes: &[u8]) {
        let first = self.record.length;
        let count = bytes.len() as u64;
        // Checkpoints at every multiple of CHECKPOINT_EVERY in [first, first + count).
        let mut next = first.div_ceil(CHECKPOINT_EVERY) * CHECKPOINT_EVERY;
        while next < first + count {
            self.checkpoints.push(offset + (next - first));
            next += CHECKPOINT_EVERY;
        }
        for &byte in bytes {
            self.record.residue_counts[byte as usize] += 1;
        }
        self.record.length += count;
        self.line_residues += count;
    }

    /// Ends a sequence line whose trailing whitespace (the uncommitted run) is `tail`.
    /// `tail` includes the `\n` when the line has one.
    fn end_line(&mut self, tail: &[u8], line_start: u64) {
        let residues = std::mem::take(&mut self.line_residues);
        if residues == 0 {
            if self.residue_lines > 0 {
                self.blank_after_residues = true;
            } else if !tail.is_empty() {
                // A line before the first residue line: residue 0 is not at the
                // start of the sequence.
                self.uniform_possible = false;
            }
            return;
        }
        let eol = match tail {
            b"\n" => Some(1),
            b"\r\n" => Some(2),
            b"" => None, // the last line of the input
            _ => {
                self.uniform_possible = false;
                None
            }
        };
        if self.blank_after_residues || self.short_line_seen {
            self.uniform_possible = false;
        }
        match self.uniform {
            None => {
                if line_start != self.record.sequence_offset {
                    self.uniform_possible = false;
                }
                self.uniform = Some((residues, eol.unwrap_or(1)));
            }
            Some((width, first_eol)) => {
                if residues > width || (eol.is_some() && eol != Some(first_eol)) {
                    self.uniform_possible = false;
                }
                if residues < width {
                    self.short_line_seen = true;
                }
            }
        }
        self.residue_lines += 1;
    }

    fn finish(mut self, end_offset: u64) -> (ScanRecord, bool) {
        self.record.end_offset = end_offset;
        self.record.layout = match (self.uniform_possible, self.uniform) {
            (true, Some((width, eol))) => LineLayout::Uniform { width, eol },
            (true, None) => LineLayout::Uniform { width: 0, eol: 1 },
            (false, _) => LineLayout::Checkpoints {
                offsets: self.checkpoints,
            },
        };
        let empty = self.record.id.is_empty() && !self.has_description && self.record.length == 0;
        (self.record, empty)
    }
}

enum Line {
    /// Nothing has been read yet.
    First,
    /// The first line does not start with `>`. bio reads (and validates) the whole
    /// line before it reports that.
    NotARecord,
    /// Inside a header line; its bytes so far.
    Header(Vec<u8>),
    /// Inside a sequence line.
    Sequence,
}

/// A streaming scan of one input.
pub struct Scanner {
    offset: u64,
    line_start: u64,
    at_line_start: bool,
    line: Line,
    current: Option<RecordBuilder>,
    records: Vec<ScanRecord>,
    /// Bytes of an incomplete UTF-8 character and its first offset.
    partial: Vec<u8>,
    partial_offset: u64,
    /// Whitespace on the current sequence line not yet known to be interior.
    pending: Vec<u8>,
    pending_offset: u64,
    /// The last record was empty: bio has read (and validated) the header line that
    /// ended it, and stops after that line.
    draining: bool,
    stopped: bool,
    error: Option<String>,
}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}

impl Scanner {
    pub fn new() -> Self {
        Self {
            offset: 0,
            line_start: 0,
            at_line_start: true,
            line: Line::First,
            current: None,
            records: Vec::new(),
            partial: Vec::new(),
            partial_offset: 0,
            pending: Vec::new(),
            pending_offset: 0,
            draining: false,
            stopped: false,
            error: None,
        }
    }

    /// Scans the next chunk of the input.
    pub fn feed(&mut self, chunk: &[u8]) {
        for &byte in chunk {
            if self.stopped || self.error.is_some() {
                return;
            }
            let offset = self.offset;
            self.offset += 1;
            if self.partial.is_empty() && byte.is_ascii() {
                self.character(byte as char, &[byte], offset);
                continue;
            }
            if self.partial.is_empty() {
                self.partial_offset = offset;
            }
            self.partial.push(byte);
            match std::str::from_utf8(&self.partial) {
                Ok(text) => {
                    let character = text.chars().next().expect("one character");
                    let bytes = std::mem::take(&mut self.partial);
                    self.character(character, &bytes, self.partial_offset);
                }
                Err(error) if error.error_len().is_none() && self.partial.len() < 4 => {}
                Err(_) => self.error = Some("stream did not contain valid UTF-8".into()),
            }
        }
    }

    /// Ends the input and returns the records, or the parser's error.
    pub fn finish(mut self) -> Result<Vec<ScanRecord>, String> {
        if let Some(error) = self.error {
            return Err(error);
        }
        if !self.partial.is_empty() {
            return Err("stream did not contain valid UTF-8".into());
        }
        if self.stopped || self.draining {
            return Ok(self.records);
        }
        let end = self.offset;
        match std::mem::replace(&mut self.line, Line::Sequence) {
            Line::First => {}
            Line::NotARecord => return Err("Expected > at record start.".into()),
            Line::Header(bytes) => {
                self.end_header(&bytes, end)?;
            }
            Line::Sequence => {
                if !self.at_line_start {
                    let tail = std::mem::take(&mut self.pending);
                    if let Some(current) = self.current.as_mut() {
                        current.end_line(&tail, self.line_start);
                    }
                }
            }
        }
        self.close_record(end);
        Ok(self.records)
    }

    fn character(&mut self, character: char, bytes: &[u8], offset: u64) {
        if self.at_line_start {
            self.at_line_start = false;
            self.line_start = offset;
            if character == '>' {
                if matches!(self.line, Line::Sequence) {
                    self.close_record(offset);
                }
                self.line = Line::Header(Vec::new());
            } else if matches!(self.line, Line::First) {
                self.line = Line::NotARecord;
            }
        }
        match &mut self.line {
            Line::NotARecord => {
                if character == '\n' {
                    self.error = Some("Expected > at record start.".into());
                }
            }
            Line::Header(header) => {
                if character == '\n' && self.draining {
                    self.stopped = true;
                } else if character == '\n' {
                    let header = std::mem::take(header);
                    if let Err(error) = self.end_header(&header, offset + 1) {
                        self.error = Some(error);
                        return;
                    }
                    self.line = Line::Sequence;
                    self.at_line_start = true;
                } else {
                    header.extend_from_slice(bytes);
                }
            }
            Line::Sequence => {
                if character == '\n' {
                    let mut tail = std::mem::take(&mut self.pending);
                    tail.push(b'\n');
                    if let Some(current) = self.current.as_mut() {
                        current.end_line(&tail, self.line_start);
                    }
                    self.at_line_start = true;
                } else if character.is_whitespace() {
                    if self.pending.is_empty() {
                        self.pending_offset = offset;
                    }
                    self.pending.extend_from_slice(bytes);
                } else {
                    let pending = std::mem::take(&mut self.pending);
                    if let Some(current) = self.current.as_mut() {
                        if !pending.is_empty() {
                            current.commit(self.pending_offset, &pending);
                        }
                        current.commit(offset, bytes);
                    }
                }
            }
            Line::First => unreachable!("the first line starts with '>' or is an error"),
        }
    }

    /// A complete header line (without its `\n`); `sequence_offset` follows it.
    fn end_header(&mut self, header: &[u8], sequence_offset: u64) -> Result<(), String> {
        let text = std::str::from_utf8(header).map_err(|_| "stream did not contain valid UTF-8")?;
        // bio: self.line[1..].trim_end().splitn(2, char::is_whitespace)
        let mut fields = text[1..].trim_end().splitn(2, char::is_whitespace);
        let id = fields.next().unwrap_or("").to_string();
        let has_description = fields.next().is_some();
        self.current = Some(RecordBuilder::new(
            id,
            has_description,
            self.line_start,
            sequence_offset,
        ));
        Ok(())
    }

    fn close_record(&mut self, end_offset: u64) {
        if let Some(current) = self.current.take() {
            let (record, empty) = current.finish(end_offset);
            if empty {
                // bio's `records()` ends at the first empty record, after the line
                // that ended it.
                self.draining = true;
            } else {
                self.records.push(record);
            }
        }
    }
}

/// Scans a whole input at once.
pub fn scan(input: &[u8]) -> Result<Vec<ScanRecord>, String> {
    let mut scanner = Scanner::new();
    scanner.feed(input);
    scanner.finish()
}

/// The *scan* response (docs/web/abi_v2.md §9).
pub fn to_json(records: &[ScanRecord]) -> String {
    use std::fmt::Write as _;
    let mut out = String::from("{\"records\":[");
    for (index, record) in records.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "{{\"index\":{index},\"id\":");
        crate::json::string(&mut out, &record.id);
        let _ = write!(
            out,
            ",\"header_offset\":{},\"sequence_offset\":{},\"end_offset\":{},\"length\":{},\"line_layout\":",
            record.header_offset, record.sequence_offset, record.end_offset, record.length
        );
        match &record.layout {
            LineLayout::Uniform { width, eol } => {
                let _ = write!(
                    out,
                    "{{\"kind\":\"uniform\",\"width\":{width},\"eol\":{eol}}}"
                );
            }
            LineLayout::Checkpoints { offsets } => {
                let _ = write!(
                    out,
                    "{{\"kind\":\"checkpoints\",\"every\":{CHECKPOINT_EVERY},\"offsets\":["
                );
                for (index, offset) in offsets.iter().enumerate() {
                    if index > 0 {
                        out.push(',');
                    }
                    let _ = write!(out, "{offset}");
                }
                out.push_str("]}");
            }
        }
        out.push_str(",\"residue_counts\":{");
        let mut first = true;
        for (byte, &count) in record.residue_counts.iter().enumerate() {
            if count == 0 {
                continue;
            }
            if !first {
                out.push(',');
            }
            first = false;
            // Printable ASCII residues are their own keys; other bytes are "0xNN".
            let key = if (0x21..=0x7e).contains(&byte) {
                (byte as u8 as char).to_string()
            } else {
                format!("0x{byte:02x}")
            };
            crate::json::string(&mut out, &key);
            let _ = write!(out, ":{count}");
        }
        out.push_str("}}");
    }
    out.push_str("]}");
    out
}
