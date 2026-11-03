//! The index scan of FASTA input (docs/web/abi_v2.md §9, plan TD-8). It exists only so
//! that the application can extract original records; every search reads its input
//! with the program's own reader, and `register` checks that the scan of the program's
//! kind (1 or 2) agrees with it.
//!
//! Parser kinds 1 (nucleotide input) and 2 (protein input) give the records of NCBI
//! BLAST+'s reader as the engine reads them (`ncbi.rs`). Kind 0, the `bio::io::fasta`
//! 1.6.0 reader that BLASTP, TBLASTN, BLASTN and TBLASTX used before session SF, was
//! removed after the application switched to kinds 1 and 2 (plan DW-38).

mod ncbi;

pub use ncbi::NcbiScanner;

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
    /// The title up to its first space (empty without a title), with bytes that are not
    /// UTF-8 replaced.
    pub id: String,
    /// Offset of the defline's `>` (0 for a first record without a defline).
    pub header_offset: u64,
    /// Offset of the first byte after the defline (0 for a first record without a defline).
    pub sequence_offset: u64,
    /// Offset of the next defline, or the end of the input.
    pub end_offset: u64,
    /// The number of residues that the reader stores.
    pub length: u64,
    pub layout: LineLayout,
    /// Count of every stored residue, upper-cased.
    pub residue_counts: Box<[u64; 256]>,
}

/// A streaming scan of one input with a parser kind.
pub struct Scanner(NcbiScanner);

impl Scanner {
    /// A scan with parser kind `kind` (1 or 2).
    pub fn of_kind(kind: u32) -> Option<Self> {
        match kind {
            1 => Some(Self(NcbiScanner::new(false))),
            2 => Some(Self(NcbiScanner::new(true))),
            _ => None,
        }
    }

    /// Scans the next chunk of the input.
    pub fn feed(&mut self, chunk: &[u8]) {
        self.0.feed(chunk);
    }

    /// Ends the input and returns the records, or the scan's error.
    pub fn finish(self) -> Result<Vec<ScanRecord>, String> {
        self.0.finish()
    }
}

/// Scans a whole input at once with parser kind 1 (`protein` false) or 2.
pub fn scan_ncbi(input: &[u8], protein: bool) -> Result<Vec<ScanRecord>, String> {
    let mut scanner = NcbiScanner::new(protein);
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
