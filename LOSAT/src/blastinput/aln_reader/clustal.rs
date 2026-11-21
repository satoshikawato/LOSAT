//! NCBI's Clustal scanner (`CAlnScannerClustal`) and the verification of the base scanner
//! (`CAlnScanner::xVerifyAlignmentData`): the ids and the rows of a Clustal file, or the failure
//! that NCBI's scanner throws (every throw of the scanner is a failed read, `super::read`).

use super::base::{verify_rows, ScanFailure, Scanned, SeqInfo};
use super::lines::PeekAheadStream;
use crate::ncbi_str::{string_to_int, truncate_spaces};
use crate::utils::format_guess::tokens;

struct BlockInfo {
    seq_count: usize,
    line_length: usize,
    first: bool,
}

/// `sTerminateBlock`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:80-98
fn terminate_block(num_seqs: &mut usize, block: &mut BlockInfo) -> Result<(), ScanFailure> {
    if block.first {
        *num_seqs = block.seq_count;
        block.first = false;
    } else if *num_seqs != block.seq_count {
        return Err(ScanFailure::Failed);
    }
    block.seq_count = 0;
    block.line_length = 0;
    Ok(())
}

/// `sIsConservationLine`: the result of `line.find_first_not_of(":.* \t\r\n")` taken as a `bool`,
/// so true unless that position is 0: a line that starts with `:`, `.` or `*` is a conservation
/// line also when other text follows, and a line of nothing else is one too (`npos` is not 0).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:200-206
fn is_conservation_line(line: &[u8]) -> bool {
    line.iter().position(|c| !b":.* \t\r\n".contains(c)) != Some(0)
}

/// `CAlnScannerClustal::xImportAlignmentData` and `sProcessClustalDataLine`, then the base
/// scanner's `xVerifyAlignmentData`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:117-196
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:219-303
/// The base scanner (`xGetExistingSeqIdInfo`: the same id, or one that differs by case only, is
/// a duplicate):
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner.cpp:205-223
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Result<Scanned, ScanFailure> {
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut sequences: Vec<Vec<&[u8]>> = Vec::new();
    let mut in_block = false;
    let mut num_seqs = 0usize;
    let mut block = BlockInfo {
        seq_count: 0,
        line_length: 0,
        first: true,
    };
    while let Some((raw, line_count)) = input.read_line() {
        if line_count == 1 && raw.len() >= 7 && raw[..7].eq_ignore_ascii_case(b"clustal") {
            continue;
        }
        let line = truncate_spaces(raw);
        if line.is_empty() {
            if in_block {
                terminate_block(&mut num_seqs, &mut block)?;
                in_block = false;
            }
            continue;
        }
        if is_conservation_line(line) {
            if !in_block {
                return Err(ScanFailure::Failed);
            }
            terminate_block(&mut num_seqs, &mut block)?;
            in_block = false;
            continue;
        }
        // `NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize)`: the repository's one port
        // (`utils::format_guess::tokens`), two or three tokens: the id, the data and a count.
        let tokens = tokens(line, b" \t");
        let (Some(id), Some(data)) = (tokens.get(0), tokens.get(1)) else {
            return Err(ScanFailure::Failed);
        };
        if tokens.len() > 3 {
            return Err(ScanFailure::Failed);
        }
        if let Some(count) = tokens.get(2) {
            // `StringToInt(..., fConvErr_NoThrow)`: 0 for a token that is no number.
            if string_to_int(count, 0).unwrap_or(0) == 0 {
                return Err(ScanFailure::Failed);
            }
        }
        in_block = true;
        if block.first {
            if ids.iter().any(|existing| existing.eq_ignore_ascii_case(id)) {
                return Err(ScanFailure::Failed);
            }
            ids.push(id.to_vec());
            sequences.push(Vec::new());
        } else {
            if block.seq_count >= num_seqs {
                return Err(ScanFailure::Failed);
            }
            if id != ids[block.seq_count].as_slice() {
                return Err(ScanFailure::Failed);
            }
        }
        if block.seq_count == 0 {
            block.line_length = data.len();
        } else if data.len() != block.line_length {
            return Err(ScanFailure::Failed);
        }
        sequences[block.seq_count].push(data);
        block.seq_count += 1;
    }
    if in_block {
        return Err(ScanFailure::Failed);
    }
    let rows: Vec<Vec<u8>> = sequences.iter().map(|parts| parts.concat()).collect();
    let info = SeqInfo::protein();
    verify_rows(&info, &rows)?;
    Ok(Scanned {
        ids,
        rows,
        deflines: Vec::new(),
        info,
    })
}
