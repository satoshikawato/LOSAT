//! The index scan of FASTA input (docs/web/abi_v2.md §9, plan TD-8). It exists only so
//! that the application can extract original records; every search reads its input
//! with the program's own reader, and `register` checks that the scan agrees with it.
//!
//! The scan reads the input with the reader of a program's input (`ReaderConfig`: kind 1 for
//! nucleotide input, kind 2 for protein input, NCBI's data loaders on or off) and gives the
//! records of NCBI BLAST+'s reader as the engine reads them (`ncbi.rs`), with the reader's errors
//! in the command line's words. Kind 0, the `bio::io::fasta` 1.6.0 reader that BLASTP, TBLASTN,
//! BLASTN and TBLASTX used before session SF, was removed after the application switched to kinds
//! 1 and 2 (plan DW-38).

mod ncbi;

use losat::web_adapter::ReaderConfig;

pub use ncbi::NcbiScanner;

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
    /// The number of letters that the reader stores (a gap line's too).
    pub length: u64,
    /// Count of every stored letter, upper-cased.
    pub residue_counts: Box<[u64; 256]>,
}

/// Why a scan failed: the reader's error or its refusal, which the reading of the input with the
/// same reader gives too, as the command line's reading of the file does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScanError {
    /// The reader's error (`CheckDataLine`): NCBI's text, without the command line's
    /// `BLAST query error: `.
    Reader(String),
    /// The reader's refusal of the input in the command line's words: a first line that NCBI's
    /// data loaders would fetch as a Seq-id, or a record longer than LOSAT reads. `record` is the
    /// offset where the refused record starts; the records before it are whole.
    Refused { message: String, record: u64 },
}

impl ScanError {
    /// The message, as `scan_end` reports it.
    pub fn message(&self) -> &str {
        match self {
            ScanError::Reader(message) | ScanError::Refused { message, .. } => message,
        }
    }
}

impl std::fmt::Display for ScanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.message())
    }
}

/// A streaming scan of one input with the reader of a program's input.
pub struct Scanner(NcbiScanner);

impl Scanner {
    /// A scan with the reader `config` (`store::reader_config`).
    pub fn new(config: ReaderConfig) -> Self {
        Self(NcbiScanner::new(config))
    }

    /// Scans the next chunk of the input.
    pub fn feed(&mut self, chunk: &[u8]) {
        self.0.feed(chunk);
    }

    /// Ends the input and returns the records, or the scan's error.
    pub fn finish(self) -> Result<Vec<ScanRecord>, ScanError> {
        self.0.finish()
    }
}

/// Scans a whole input at once with the reader `config`.
pub fn scan_ncbi(input: &[u8], config: ReaderConfig) -> Result<Vec<ScanRecord>, ScanError> {
    let mut scanner = NcbiScanner::new(config);
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
            ",\"header_offset\":{},\"sequence_offset\":{},\"end_offset\":{},\"length\":{}",
            record.header_offset, record.sequence_offset, record.end_offset, record.length
        );
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
