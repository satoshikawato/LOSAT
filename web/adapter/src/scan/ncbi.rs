//! The index scan (docs/web/abi_v2.md §9): the records that NCBI BLAST+'s FASTA input source reads,
//! as the engine reads them (the engine's `fasta_reader`, the port of `CBlastFastaInputSource`,
//! `CFastaReader` and `CStreamLineReader`), from the bytes of a file
//! (`losat::web_adapter::FastaInputSource::from_bytes`), with the reader of a program's input
//! (`ReaderConfig`: the program, the role, its flags and whether NCBI's data loaders are on).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-345
//! Kind 1 has the flags of nucleotide input (`fAssumeNuc`: BLASTN, TBLASTX, the TBLASTN
//! subject, the BLASTX query), kind 2 those of protein input (`fAssumeProt`: BLASTP, the TBLASTN
//! query, the BLASTX subject).
//!
//! The scan takes the input in chunks of any size and holds only a header's first word, the
//! first line of an input whose first line may be a Seq-id (with the data loaders on), the line
//! of a `>?` gap line, at most 70 bytes of a record's first data line and at most 257 bytes whose
//! loss is not yet decided. It writes no messages. Each record is reported with every record that
//! the reader returns (records without residues too): the ID is the title up to its first space
//! (empty without a title), the length is the number of letters that the reader stores (a gap
//! line's letters too, counted by the reader's own parser of gap lines,
//! `losat::web_adapter::gap_line_letters`), and `residue_counts` counts them upper-cased (`U` is
//! `T` for kind 1, as the reader stores it; a gap's letters are `N` for kind 1, `X` for kind 2).
//!
//! The scan fails with the reader's errors, in the reader's order and words: `CheckDataLine`, and
//! the refusals of the reader of the program's input, which the command line makes too: a first
//! line that NCBI's data loaders would read as a Seq-id (with the data loaders on), and a record
//! longer than 2147483647 letters (`losat::web_adapter::record_too_long`, which the scan gives
//! because such a record does not fit in the reader's memory in the browser). Until line WQA's
//! D13 (2026-11-20) the scan also reported where each residue was (`line_layout`), which no
//! application read, and refused a `>?` gap line and lines that NCBI's line reader joins for it:
//! inputs that the command line reads.

use losat::web_adapter::{
    defline_named_molecule, gap_line_letters, record_too_long, DataLineAlphabet, FastaInputSource,
    ReadError, ReaderConfig, MAX_RECORD_LETTERS,
};

use super::{ScanError, ScanRecord};

/// Bytes of a regular file that an `ifstream` buffers at once (`FastaStream::new` in
/// `LOSAT/src/blastinput/fasta_reader/stream.rs`, checked against an independent C++ oracle).
const FILE_WINDOW: u64 = 8191;

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:117
const MIN_BUF_SIZE: u64 = 4096;

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:413-417
const STEP_BACK_MAX: usize = (MIN_BUF_SIZE >> 4) as usize;

/// `CheckDataLine`'s window.
const CHECK_WINDOW: usize = 70;

/// Sentinel values in the per-byte residue mapping tables (`residue_table`): non-residue bytes
/// map to `BYTE_OTHER` (`0`), `;` maps to `BYTE_COMMENT` (`1`), and every valid residue maps
/// directly to its stored upper-case ASCII byte (`b'A'..=b'Z'` or `b'*'`, all `> 1`).
const BYTE_OTHER: u8 = 0;
const BYTE_COMMENT: u8 = 1;

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:854-935
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1448
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:392-457
// `bIsNuc` is `fAssumeNuc`; `bAllowLetterGaps` is false (no `fParseGaps`); `-`, white space
// and bad bytes are not stored (`fHyphensIgnoreAndWarn`), so they map to `BYTE_OTHER`.
// When `program_protein` is false: `U`/`u` maps to `b'T'`, and if `parse_protein` is true
// (`CFastaReader::xSetSeqMol` named protein), non-IUPACna protein letters map to `b'N'`
// (`pack_protein_read_nucleotide`). When `program_protein` is true, kept letters are stored
// upper-cased (`U` stays `b'U'`).
const fn residue_table(parse_protein: bool, program_protein: bool) -> [u8; 256] {
    let mut table = [BYTE_OTHER; 256];
    let nucleotide = b"ABCDGHKMRSTUVWYN";
    let mut index = 0;
    while index < nucleotide.len() {
        let upper = nucleotide[index];
        let stored = if !program_protein && upper == b'U' {
            b'T'
        } else {
            upper
        };
        table[upper as usize] = stored;
        table[upper.to_ascii_lowercase() as usize] = stored;
        index += 1;
    }
    if parse_protein {
        let protein_only = b"EFIJLOPQZX";
        let mut index = 0;
        while index < protein_only.len() {
            let upper = protein_only[index];
            let stored = if program_protein { upper } else { b'N' };
            table[upper as usize] = stored;
            table[upper.to_ascii_lowercase() as usize] = stored;
            index += 1;
        }
        table[b'*' as usize] = if program_protein { b'*' } else { b'N' };
    }
    table[b';' as usize] = BYTE_COMMENT;
    table
}

const NUC_FOR_NUC: [u8; 256] = residue_table(false, false);
const PROT_FOR_NUC: [u8; 256] = residue_table(true, false);
const NUC_FOR_PROT: [u8; 256] = residue_table(false, true);
const PROT_FOR_PROT: [u8; 256] = residue_table(true, true);

#[inline]
const fn select_residue_table(parse_protein: bool, program_protein: bool) -> &'static [u8; 256] {
    match (parse_protein, program_protein) {
        (false, false) => &NUC_FOR_NUC,
        (true, false) => &PROT_FOR_NUC,
        (false, true) => &NUC_FOR_PROT,
        (true, true) => &PROT_FOR_PROT,
    }
}

/// `isspace` of the C locale (`TruncateSpaces_Unsafe`, `ParseDefline`); CR and LF never
/// reach a line's content.
const fn space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// What a byte of the input is to NCBI's line reader.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Role {
    /// A byte of a line.
    Content,
    /// A byte of a line's end of line.
    Delimiter,
    /// An end-of-line byte that the line reader consumes without ending a line (the LF
    /// after a CR that it pushed back, the CR after an LF that it pushed back).
    Dropped,
    /// A byte that the line reader never returns (the tail lost to the pushback's
    /// end-of-file state).
    Lost,
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:167-173
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Eol {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}

/// The pushback of a line's tail after a CR in the CRLF style (the first pushback of a
/// file that can be followed by a second one).
struct FirstPush {
    /// Offset of the first pushed-back byte.
    start: u64,
    /// Offset of the LF that ended the physical line (consumed, not pushed back).
    lf: Option<u64>,
}

/// NCBI's `CStreamLineReader` over the bytes of a file, pushed in chunks: each byte goes to
/// the record reader with its role, in the order of the input. The line reader's states
/// are those of `FastaStream` and `LineReader` in the engine; a pushed-back tail is read
/// again in the input's order, so the pushback is the end-of-line byte that it skips
/// (`Role::Dropped`) and the style that reads the tail.
struct Lines {
    eol: Eol,
    /// Offset of the next byte.
    offset: u64,
    /// A CR ended the line (`Unknown`, `Cr`, `Mixed`): an LF next is part of its end of line.
    cr_end: bool,
    /// The line has a CR (`Lf`, `CrLf`): an LF next makes it a CRLF.
    lf_cr: bool,
    drop_lf: bool,
    drop_cr: bool,
    first_push: Option<FirstPush>,
    /// The bytes after an LF of a CR-style line read through the pushback buffer (the
    /// second pushback, whose bytes may be lost) and the offset of the first.
    held: Option<(u64, Vec<u8>)>,
}

impl Lines {
    fn new() -> Self {
        Self {
            eol: Eol::Unknown,
            offset: 0,
            cr_end: false,
            lf_cr: false,
            drop_lf: false,
            drop_cr: false,
            first_push: None,
            held: None,
        }
    }

    fn byte(&mut self, byte: u8, records: &mut Records) {
        let offset = self.offset;
        self.offset += 1;
        if let Some((_, held)) = self.held.as_mut() {
            if byte != b'\r' {
                held.push(byte);
                if held.len() > STEP_BACK_MAX {
                    // More than 256 bytes are pushed back into a new buffer, which clears
                    // the end-of-file state: nothing is lost.
                    self.release(records);
                }
                return;
            }
            // The CR that ends the line: the second pushback is not at the end of the input.
            self.release(records);
        }
        self.visible(byte, offset, records);
    }

    fn release(&mut self, records: &mut Records) {
        if let Some((start, held)) = self.held.take() {
            for (index, &byte) in held.iter().enumerate() {
                self.visible(byte, start + index as u64, records);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-243
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:245-268
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:271-281
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:90-102
    /// One byte that the line reader reads (not held). A CR in the LF or CRLF style that
    /// is not followed by its LF ends the line and pushes the rest of the physical line
    /// back: the LF that ended the physical line was consumed (`drop_lf`), and the rest is
    /// read in the mixed style (from LF) or in the CR style (from CRLF). An LF inside a
    /// CR-style line ends the line and pushes the rest back: the CR that ended the read was
    /// consumed (`drop_cr`), and the rest is read in the mixed style.
    fn visible(&mut self, byte: u8, offset: u64, records: &mut Records) {
        if byte == b'\n' && self.drop_lf {
            self.drop_lf = false;
            if let Some(push) = self.first_push.as_mut() {
                push.lf = Some(offset);
            }
            records.byte(byte, offset, Role::Dropped);
            return;
        }
        if byte == b'\r' && self.drop_cr {
            self.drop_cr = false;
            records.byte(byte, offset, Role::Dropped);
            return;
        }
        if std::mem::take(&mut self.cr_end) {
            if self.eol == Eol::Unknown {
                self.eol = if byte == b'\n' { Eol::CrLf } else { Eol::Cr };
            }
            if byte == b'\n' {
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                return;
            }
            records.line_end();
        }
        if std::mem::take(&mut self.lf_cr) {
            if byte == b'\n' {
                // CR LF: the style does not change.
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                return;
            }
            // A CR inside the line: it ends the line; this byte starts the pushed-back rest.
            records.line_end();
            self.drop_lf = true;
            if self.eol == Eol::CrLf {
                self.eol = Eol::Cr;
                self.first_push = Some(FirstPush {
                    start: offset,
                    lf: None,
                });
            } else {
                self.eol = Eol::Mixed;
            }
        }
        match (self.eol, byte) {
            (Eol::Unknown | Eol::Mixed | Eol::Cr, b'\r') => {
                records.byte(byte, offset, Role::Delimiter);
                self.cr_end = true;
            }
            (Eol::Unknown | Eol::Mixed, b'\n') => {
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                if self.eol == Eol::Unknown {
                    self.eol = Eol::CrLf;
                }
            }
            (Eol::Cr, b'\n') => {
                // The first LF of a CR-style line.
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                self.eol = Eol::Mixed;
                self.drop_cr = true;
                if self.first_push.is_some() {
                    // A second pushback: its bytes go back into the first one's buffer when
                    // they fit, and the end-of-file state stays if the read ended there.
                    self.held = Some((offset + 1, Vec::new()));
                }
            }
            (Eol::Lf | Eol::CrLf, b'\r') => {
                records.byte(byte, offset, Role::Delimiter);
                self.lf_cr = true;
            }
            (Eol::Lf | Eol::CrLf, b'\n') => {
                records.byte(byte, offset, Role::Delimiter);
                records.line_end();
                if self.eol == Eol::CrLf {
                    self.eol = Eol::Lf;
                }
            }
            _ => records.byte(byte, offset, Role::Content),
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:418-434
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    /// The end of the input. A second pushback at the end of the input (the CR-style line
    /// ended there) whose bytes all step back into the current pushback buffer keeps the
    /// stream's end-of-file state: `AtEOF` is true and those bytes are lost. A new buffer
    /// (more than 256 bytes, or more than the buffer has read) clears the state.
    fn finish(&mut self, records: &mut Records) {
        if let Some((start, held)) = self.held.take() {
            let lost = match &self.first_push {
                Some(FirstPush {
                    start: first,
                    lf: Some(lf),
                }) => {
                    !held.is_empty()
                        && last_refill(lf + 1, lf - first, self.offset) >= held.len() as u64
                }
                _ => false,
            };
            if lost {
                for (index, &byte) in held.iter().enumerate() {
                    records.byte(byte, start + index as u64, Role::Lost);
                }
            } else {
                self.held = Some((start, held));
                self.release(records);
            }
        }
        // A CR at the end ends the line (in the LF and CRLF styles too: `pos` is the size of
        // `m_Line`, nothing is pushed back).
        self.cr_end = false;
        self.lf_cr = false;
        records.line_end();
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-327
/// The length of the last refill of the first pushback buffer before the end of a file of
/// `total` bytes: the buffer holds the `pushed` bytes of the first pushback, which ended
/// with the LF before `resume`, and then refills with `in_avail()` bytes, at most its size
/// (at least 4096). `in_avail()` of the file's stream is the rest of its 8191-byte window,
/// or the rest of the file when the window is used up; a read of a whole window or more
/// with the window used up bypasses it (`FastaStream::raw_peek` in the engine, whose
/// windows start at multiples of 8191 until the first pushback).
fn last_refill(resume: u64, pushed: u64, total: u64) -> u64 {
    let capacity = pushed.max(MIN_BUF_SIZE);
    let mut position = resume;
    let mut window = if resume.is_multiple_of(FILE_WINDOW) {
        0
    } else {
        ((resume / FILE_WINDOW + 1) * FILE_WINDOW).min(total) - resume
    };
    let mut last = pushed;
    loop {
        let read = if window > 0 {
            let read = capacity.min(window);
            window -= read;
            read
        } else {
            let rest = total - position;
            if rest == 0 {
                return last;
            }
            let read = capacity.min(rest);
            if read < FILE_WINDOW {
                window = FILE_WINDOW.min(rest) - read;
            }
            read
        };
        position += read;
        last = read;
    }
}

/// How a line that starts with `>` begins (at most 4 bytes).
#[derive(Clone, Copy)]
struct Prefix {
    bytes: [u8; 4],
    len: usize,
}

impl Prefix {
    fn new(byte: u8) -> Self {
        Self {
            bytes: [byte; 4],
            len: 1,
        }
    }

    fn push(&mut self, byte: u8) {
        self.bytes[self.len] = byte;
        self.len += 1;
    }

    /// The bytes after the `>` and its `?_` (`skip` 3), or after the `>` alone (`skip` 1).
    fn after(&self, skip: usize) -> &[u8] {
        &self.bytes[skip..self.len]
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-374
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:402-409
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-779
    /// What the line is once enough of it is known (`at_end`: the line has no more bytes):
    /// `Some(None)` a gap line, `Some(Some(skip))` not a gap line whose `>?_` (when
    /// `skip` is 3) becomes `>`, `None` not known yet.
    fn kind(&self, at_end: bool) -> Option<Option<usize>> {
        let byte = |index: usize| (index < self.len).then(|| self.bytes[index]);
        match (byte(1), byte(2), byte(3)) {
            (None, ..) => at_end.then_some(Some(1)),
            (Some(b'?'), None, _) => at_end.then_some(None),
            (Some(b'?'), Some(b'_'), None) => at_end.then_some(Some(3)),
            (Some(b'?'), Some(b'_'), Some(b'?')) => Some(None),
            (Some(b'?'), Some(b'_'), Some(_)) => Some(Some(3)),
            (Some(b'?'), Some(_), _) => Some(None),
            (Some(_), ..) => Some(Some(1)),
        }
    }
}

/// The ID of a defline being read: the title's first word.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-165
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:168-180
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:209-225
/// `x_ApplyMods` then drops the title's white space at both ends (its first byte is not
/// white space, and a space ends the first word anyway). The ID is the title up to its
/// first space (`FastaRecord::shown_id`, plan §1).
#[derive(Default)]
struct Title {
    started: bool,
    word_done: bool,
    title_done: bool,
}

impl Title {
    fn byte(&mut self, byte: u8, id: &mut Vec<u8>) {
        if self.title_done {
            return;
        }
        if !self.started {
            if !space(byte) {
                self.started = true;
                id.push(byte);
            }
            return;
        }
        if byte < b' ' {
            self.title_done = true;
        } else if !self.word_done {
            if byte == b' ' {
                self.word_done = true;
            } else {
                id.push(byte);
            }
        }
    }
}

/// `CheckDataLine` of a record's first data lines (no residue stored yet): the first 70
/// bytes of the line (after `TruncateSpaces_Unsafe` and the `>?_` change) and how long
/// the line is up to its last byte other than white space.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:710-760
/// The warning about ambiguous nucleotides is never written (`bIsNuc` is false before
/// `AssembleSeq`).
struct Check {
    bytes: [u8; CHECK_WINDOW],
    len: usize,
    /// One past the last byte other than white space, at most 71.
    end: usize,
}

impl Check {
    fn new() -> Self {
        Self {
            bytes: [0; CHECK_WINDOW],
            len: 0,
            end: 0,
        }
    }

    fn byte(&mut self, byte: u8) {
        let index = self.len;
        if index < CHECK_WINDOW {
            self.bytes[index] = byte;
            self.len += 1;
            if !space(byte) {
                self.end = index + 1;
            }
        } else if !space(byte) {
            self.end = CHECK_WINDOW + 1;
        }
    }

    fn plausible(&self) -> bool {
        let len_to_check = self.end.min(CHECK_WINDOW);
        let (mut good, mut bad) = (0usize, 0usize);
        for &c in &self.bytes[..len_to_check] {
            if c.is_ascii_alphabetic() || c == b'*' {
                good += 1;
            } else if c == b'-' || space(c) || c.is_ascii_digit() {
            } else if c == b';' {
                break;
            } else {
                bad += 1;
            }
        }
        !(bad >= good / 3 && (len_to_check > 3 || good == 0 || bad > good))
    }
}

/// The state of the line being read.
enum LineState {
    /// No byte of the line yet.
    Fresh,
    /// The line starts with `>`: a defline or a gap line.
    Angle(Prefix),
    /// The defline of the current record.
    Title(Title),
    /// White space at the start of a line that does not start with `>`.
    Leading,
    /// A data line whose first byte other than white space is `>` (a gap line or not),
    /// with its `CheckDataLine`.
    DataAngle(Prefix, Option<Check>),
    /// A data line: whether a `;` ended its data, and its `CheckDataLine`.
    Data {
        semicolon: bool,
        check: Option<Check>,
    },
    /// A `>?` gap line: the line as `ParseGapLine` reads it (after `>?_` became `>`), whose
    /// letters are counted at its end.
    Gap(Vec<u8>),
    /// A comment line.
    Skip,
}

/// A record being read.
struct Builder {
    id: Vec<u8>,
    header_offset: u64,
    /// Whether the record starts with a defline (`>`), not the input's first record without one.
    defline: bool,
    sequence_offset: u64,
    /// The letters that the reader stores, a gap line's letters too.
    length: u64,
    /// The residues of the data lines (the reader's `m_SeqData`, without the gaps): while there
    /// is none, a data line gets `CheckDataLine`.
    residues: u64,
    residue_counts: Box<[u64; 256]>,
}

impl Builder {
    fn new(header_offset: u64, sequence_offset: u64, defline: bool) -> Self {
        Self {
            id: Vec::new(),
            header_offset,
            defline,
            sequence_offset,
            length: 0,
            residues: 0,
            residue_counts: Box::new([0; 256]),
        }
    }

    /// Stores a residue of a data line (upper-cased; a nucleotide `U` as `T`).
    fn store(&mut self, residue: u8) {
        self.residue_counts[residue as usize] += 1;
        self.residues += 1;
        self.length += 1;
    }

    /// Stores the `letters` of a gap line: the run of `gap` (`N` or `X`) that `AssembleSeq` puts
    /// in the record's place of the line.
    fn gap(&mut self, gap: u8, letters: u32) {
        self.residue_counts[gap as usize] += u64::from(letters);
        self.length += u64::from(letters);
    }

    fn finish(self, end_offset: u64) -> ScanRecord {
        ScanRecord {
            id: String::from_utf8_lossy(&self.id).into_owned(),
            header_offset: self.header_offset,
            sequence_offset: self.sequence_offset,
            end_offset,
            length: self.length,
            residue_counts: self.residue_counts,
        }
    }
}

/// The record reader (`CFastaReader::ReadOneSeq` as `CBlastInputReader` calls it), fed
/// with the line reader's bytes.
struct Records {
    /// The reader of the program's input whose records these are.
    config: ReaderConfig,
    /// Residue mapping table for the current record, decided once at `Records::new` for the
    /// default path and updated once per record at defline end when `-parse_deflines` is active.
    residues: &'static [u8; 256],
    line_number: u64,
    in_line: bool,
    line_start: u64,
    state: LineState,
    /// The first line from its first byte other than white space, held while it may be a
    /// Seq-id.
    first_line: Option<Vec<u8>>,
    /// No record has started (`need_defline` of the first `ReadOneSeq`).
    need_defline: bool,
    current: Option<Builder>,
    /// The current record's `sequence_offset` is the next byte that the line reader
    /// reads after its defline.
    awaiting_sequence: bool,
    records: Vec<ScanRecord>,
    error: Option<ScanError>,
    /// The reader's own error that ended the records.
    stop: Option<ReaderStop>,
}

impl Records {
    fn new(config: ReaderConfig) -> Self {
        Self {
            config,
            residues: select_residue_table(config.protein, config.protein),
            line_number: 0,
            in_line: false,
            line_start: 0,
            state: LineState::Fresh,
            first_line: None,
            need_defline: true,
            current: None,
            awaiting_sequence: false,
            records: Vec::new(),
            error: None,
            stop: None,
        }
    }

    /// Whether the scan has ended: at a refusal, or at the reader's error.
    fn halted(&self) -> bool {
        self.error.is_some() || self.stop.is_some()
    }

    fn fail(&mut self, error: ScanError) {
        if !self.halted() {
            self.error = Some(error);
        }
    }

    /// The reader's error `message` on the current record: the records read before it are the
    /// input's records, and the error ends them (`ReaderStop`).
    fn stop(&mut self, message: String) {
        if !self.halted() {
            let (offset, defline) = self.current.as_ref().map_or((0, false), |current| {
                (current.header_offset, current.defline)
            });
            self.stop = Some(ReaderStop {
                offset,
                defline,
                message,
            });
        }
    }

    /// One byte of the input, in the input's order, with its role.
    fn byte(&mut self, byte: u8, offset: u64, role: Role) {
        if self.halted() {
            return;
        }
        if matches!(role, Role::Content | Role::Delimiter) {
            if std::mem::take(&mut self.awaiting_sequence) {
                if let Some(current) = self.current.as_mut() {
                    current.sequence_offset = offset;
                }
            }
            if !self.in_line {
                // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-161
                self.in_line = true;
                self.line_number += 1;
                self.line_start = offset;
                self.state = LineState::Fresh;
            }
        }
        if role == Role::Content {
            self.content(byte);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:339-340
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:376-389
    /// A byte of the current line.
    fn content(&mut self, byte: u8) {
        if let Some(line) = self.first_line.as_mut() {
            line.push(byte);
        }
        match std::mem::replace(&mut self.state, LineState::Skip) {
            LineState::Fresh if byte == b'>' => self.state = LineState::Angle(Prefix::new(byte)),
            LineState::Fresh | LineState::Leading => self.leading(byte),
            LineState::Angle(mut prefix) => {
                prefix.push(byte);
                self.angle(prefix, false);
            }
            LineState::Title(mut title) => {
                if let Some(current) = self.current.as_mut() {
                    title.byte(byte, &mut current.id);
                }
                self.state = LineState::Title(title);
            }
            LineState::DataAngle(mut prefix, check) => {
                prefix.push(byte);
                self.data_angle(prefix, check, false);
            }
            LineState::Data { semicolon, check } => {
                self.state = LineState::Data { semicolon, check };
                if let LineState::Data {
                    check: Some(check), ..
                } = &mut self.state
                {
                    check.byte(byte);
                }
                self.residue(byte);
            }
            LineState::Gap(mut line) => {
                line.push(byte);
                self.state = LineState::Gap(line);
            }
            LineState::Skip => {}
        }
    }

    /// The first byte other than white space of a line that does not start with `>`.
    fn leading(&mut self, byte: u8) {
        if space(byte) {
            self.state = LineState::Leading;
            return;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:132-133
        if self.config.data_loaders && self.line_number == 1 && byte.is_ascii_alphanumeric() {
            self.first_line = Some(vec![byte]);
        }
        if matches!(byte, b'!' | b'#' | b';') {
            return;
        }
        self.start_headerless();
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:715-717
        let check = self
            .current
            .as_ref()
            .is_some_and(|current| current.residues == 0)
            .then(Check::new);
        if byte == b'>' {
            self.state = LineState::DataAngle(Prefix::new(byte), check);
            return;
        }
        let mut check = check;
        if let Some(check) = check.as_mut() {
            check.byte(byte);
        }
        self.state = LineState::Data {
            semicolon: false,
            check,
        };
        self.residue(byte);
    }

    /// A data line that comes before any defline: the reader parses an empty defline there
    /// (`ParseDefLine`), so the input's first record has no title (`fDLOptional`). Every later
    /// record starts at the defline that ended the one before.
    fn start_headerless(&mut self) {
        if self.need_defline {
            self.need_defline = false;
            self.current = Some(Builder::new(0, 0, false));
        }
    }

    /// A line that starts with `>`, once its kind is known: a defline starts a record (and ends
    /// the current one), a gap line is a data line.
    fn angle(&mut self, prefix: Prefix, at_end: bool) {
        match prefix.kind(at_end) {
            None => self.state = LineState::Angle(prefix),
            Some(None) => {
                self.start_headerless();
                self.state = LineState::Gap(gap_line(&prefix));
            }
            Some(Some(skip)) => {
                if let Some(current) = self.current.take() {
                    self.close(current, self.line_start);
                }
                self.need_defline = false;
                let mut current = Builder::new(self.line_start, self.line_start, true);
                let mut title = Title::default();
                for &byte in prefix.after(skip) {
                    title.byte(byte, &mut current.id);
                }
                self.current = Some(current);
                self.state = LineState::Title(title);
            }
        }
    }

    /// A data line whose first byte other than white space is `>`, once its kind is known.
    fn data_angle(&mut self, prefix: Prefix, check: Option<Check>, at_end: bool) {
        match prefix.kind(at_end) {
            None => self.state = LineState::DataAngle(prefix, check),
            Some(None) => self.state = LineState::Gap(gap_line(&prefix)),
            Some(Some(skip)) => {
                let mut check = check;
                if let Some(check) = check.as_mut() {
                    check.byte(prefix.bytes[0]);
                    for &byte in prefix.after(skip) {
                        check.byte(byte);
                    }
                }
                self.state = LineState::Data {
                    semicolon: false,
                    check,
                };
                // `>`, `?` and `_` are not residues: only the line's last byte, when this
                // byte decided the kind, can be one.
                match (at_end, prefix.len) {
                    (false, len) if len > 1 => self.residue(prefix.bytes[len - 1]),
                    _ => {}
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:809-814
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:938-950
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:970-973
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1448
    /// A byte of a data line's data: looked up once in the current record's precomputed
    /// `self.residues` table (which encodes the record's alphabet, upper-casing, `U`→`T`, and
    /// nucleotide packing for protein-named records); `;` ends the data.
    #[inline]
    fn residue(&mut self, byte: u8) {
        let LineState::Data { semicolon, .. } = &mut self.state else {
            return;
        };
        if *semicolon {
            return;
        }
        let mapped = self.residues[byte as usize];
        if mapped <= BYTE_COMMENT {
            if mapped == BYTE_COMMENT {
                *semicolon = true;
            }
            return;
        }
        if let Some(current) = self.current.as_mut() {
            current.store(mapped);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-779
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1386-1389
    // `ParseDataLine` hands a line that starts with `>?` to the parser of gap lines before any
    // check of data; without `fParseGaps` (BLAST's reader) the record's gaps are then stored as
    // runs of `X` in protein and `N` in nucleotide data (line WQA's D13 round 3, L2; license
    // rule 1 of 2026-11-20).
    /// The end of a gap line: its letters, as the reader's parser of gap lines counts them
    /// (`gap_line_letters`), are a run of `N` (kind 1) or `X` (kind 2) in the record. A gap line
    /// gets no `CheckDataLine` and stores no residue in `m_SeqData`, so the next data line is
    /// checked as the record's first.
    fn gap(&mut self, line: &[u8]) {
        let letters = gap_line_letters(line).expect("a gap line starts with >?");
        let gap = if self.config.protein { b'X' } else { b'N' };
        if let Some(current) = self.current.as_mut() {
            current.gap(gap, letters);
        }
    }

    /// The end of the current line, if one is open.
    fn line_end(&mut self) {
        if !std::mem::take(&mut self.in_line) || self.halted() {
            return;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-162
        // `CBlastInputReader` reads only with the data loaders on (`first_line` is held only
        // then). Every later record starts at the `>` line that ended the one before, so only the
        // input's first line is tried.
        if let Some(line) = self.first_line.take() {
            if let Some(message) = seq_id_refusal(&line, self.config) {
                self.fail(ScanError::Refused { message, record: 0 });
                return;
            }
        }
        match std::mem::replace(&mut self.state, LineState::Fresh) {
            LineState::Angle(prefix) => self.angle(prefix, true),
            LineState::DataAngle(prefix, check) => self.data_angle(prefix, check, true),
            state => self.state = state,
        }
        match std::mem::replace(&mut self.state, LineState::Fresh) {
            LineState::Title(_) => {
                if self.config.parse_deflines {
                    if let Some(current) = self.current.as_ref() {
                        let mut defline = Vec::with_capacity(1 + current.id.len());
                        defline.push(b'>');
                        defline.extend_from_slice(&current.id);
                        let named = defline_named_molecule(&defline, self.config.protein);
                        let parse_protein =
                            DataLineAlphabet::new(named, self.config.protein).is_protein();
                        self.residues = select_residue_table(parse_protein, self.config.protein);
                    }
                }
                self.awaiting_sequence = true;
            }
            LineState::Data {
                check: Some(check), ..
            } if !check.plausible() => {
                let line = self.line_number;
                self.stop(format!(
                    "CFastaReader: Near line {line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
                ));
                return;
            }
            LineState::Gap(line) => self.gap(&line),
            _ => {}
        }
        // The reader checks the record's length after each data line and gap line, after
        // `CheckDataLine` (`FastaReader::check_length`).
        if let Some(current) = self
            .current
            .as_ref()
            .filter(|current| current.length > MAX_RECORD_LETTERS)
        {
            let record = current.header_offset;
            self.fail(ScanError::Refused {
                message: reader_error(record_too_long(&self.config)),
                record,
            });
        }
    }

    fn close(&mut self, current: Builder, end_offset: u64) {
        self.records.push(current.finish(end_offset));
    }

    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:428-432
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:203-208
    /// The end of the input (`total` bytes): the current record ends there; an input with
    /// only white space and comment lines has no record. The reader's own error ends the records
    /// before it (`ReaderStop`); the reader's refusals fail the scan.
    fn finish(mut self, total: u64) -> Result<Scanned, ScanError> {
        if !self.halted() {
            if let Some(mut current) = self.current.take() {
                if self.awaiting_sequence {
                    current.sequence_offset = total;
                }
                self.close(current, total);
            }
        }
        match self.error {
            Some(error) => Err(error),
            None => Ok(Scanned {
                records: self.records,
                stop: self.stop,
            }),
        }
    }
}

/// The line of a gap line that starts with `prefix`, as `ParseGapLine` reads it: `>?_` (at the
/// start of a line, after its white space) becomes `>`.
fn gap_line(prefix: &Prefix) -> Vec<u8> {
    let skip = if prefix.after(1).starts_with(b"?_") {
        3
    } else {
        1
    };
    let mut line = vec![b'>'];
    line.extend_from_slice(prefix.after(skip));
    line
}

/// An error of the engine's reader as the command line reports it, without its line end
/// (`store::reader_error`).
fn reader_error(error: ReadError) -> String {
    crate::store::reader_error(error)
}

/// The refusal of the reader of `config` (with NCBI's data loaders on) of an input whose first
/// line is `line` (from its first byte other than white space), when NCBI would read that line as
/// a Seq-id: the engine's reader decides and words it (`ReadError::Unsupported` for such a line,
/// `seq_id.rs`), as `register` and the command line refuse the input.
fn seq_id_refusal(line: &[u8], config: ReaderConfig) -> Option<String> {
    let mut bytes = Vec::with_capacity(line.len() + 1);
    bytes.extend_from_slice(line);
    bytes.push(b'\n');
    let mut source = FastaInputSource::from_bytes(&bytes, config);
    match source.next_sequence(&mut |_| Ok(())) {
        Err(error @ ReadError::Unsupported(_)) => Some(reader_error(error)),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-152
/// The records of a scan, and the reader's own error that ended them (`CheckDataLine`): the
/// command line's query batches search the records before it and end with the error in the batch
/// that reads the record (line WQA's D10, 2026-11-19).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Scanned {
    pub records: Vec<ScanRecord>,
    pub stop: Option<ReaderStop>,
}

/// The reader's error on the record after the records of a scan (`Scanned::stop`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReaderStop {
    /// The offset of the record that the reader fails on: its defline's `>` (0 for a first record
    /// without a defline). The bytes from there on are read by no record.
    pub offset: u64,
    /// Whether that record starts with a defline. One without (the input's first line, read as
    /// letters: a byte-order mark, a line of other text) cannot follow another input, whose last
    /// record would read it (review M2, 2026-11-20).
    pub defline: bool,
    /// The reader's message (`CFastaReader: Near line …`).
    pub message: String,
}

/// A streaming scan with the reader of a program's input (`config`: kind 1 for nucleotide
/// input, 2 for protein input, NCBI's data loaders on or off).
pub struct NcbiScanner {
    lines: Lines,
    records: Records,
}

impl NcbiScanner {
    pub fn new(config: ReaderConfig) -> Self {
        Self {
            lines: Lines::new(),
            records: Records::new(config),
        }
    }

    /// Scans the next chunk of the input.
    pub fn feed(&mut self, chunk: &[u8]) {
        for &byte in chunk {
            if self.records.halted() {
                return;
            }
            self.lines.byte(byte, &mut self.records);
        }
    }

    /// Ends the input and returns the records and the reader's error that ended them, or the
    /// reader's refusal.
    pub fn finish(mut self) -> Result<Scanned, ScanError> {
        if !self.records.halted() {
            self.lines.finish(&mut self.records);
        }
        let total = self.lines.offset;
        self.records.finish(total)
    }
}
