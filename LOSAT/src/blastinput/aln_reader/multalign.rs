//! NCBI's MultAlin scanner (`CAlnScannerMultAlign`): an optional `//` line, then blocks of an
//! offsets line, one data line per sequence and an optional `Consensus` line, the gap `.`.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_multalign.cpp:50-128

use super::base::{
    process_data_line, split_merge, Scan, ScanFailure, Scanned, SeqInfo, PROTEIN_ALPHABET,
};
use super::lines::PeekAheadStream;
use crate::ncbi_str::{string_to_int, truncate_spaces};

/// `xGetExpectedDataSize`: the offsets line holds one or two numbers, and with two the size is
/// `end - start + 1`, else 0 (the size of the first data line then sets it).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_multalign.cpp:213-243
fn expected_data_size(line: &[u8]) -> Scan<usize> {
    let tokens = split_merge(line, b" ");
    if tokens.len() > 2 {
        return Err(ScanFailure::Failed);
    }
    // `tokens[0]` of a line that has no token (an empty string cannot get here: the caller
    // skips empty lines).
    let start = string_to_int(tokens[0], 0).map_err(|_| ScanFailure::Failed)?;
    if tokens.len() == 2 {
        let end = string_to_int(tokens[1], 0).map_err(|_| ScanFailure::Failed)?;
        // `endOffset - startOffset + 1` in `int`, stored as `size_t`.
        let size = end.wrapping_sub(start).wrapping_add(1);
        return Ok(size as i64 as usize);
    }
    Ok(0)
}

/// `xImportAlignmentData`, `xAdjustSequenceInfo` (the gap `.`, no match or missing character) and
/// `xVerifySingleSequenceData` (alphabet and `.` anywhere).
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Scan<Scanned> {
    let Some((first, _)) = input.read_line() else {
        return Err(ScanFailure::Failed);
    };
    if first.starts_with(b"//") {
        let Some((second, _)) = input.read_line() else {
            return Err(ScanFailure::Failed);
        };
        if !second.is_empty() {
            return Err(ScanFailure::Failed);
        }
    }
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut expecting_offsets = true;
    let mut in_first_block = true;
    let mut expected_size = 0usize;
    let mut expected_num_sequences = 0usize;
    let mut line_in_block = 0usize;
    while let Some((raw, _)) = input.read_line() {
        let line = truncate_spaces(raw);
        if expecting_offsets {
            if line.is_empty() {
                continue;
            }
            expected_size = expected_data_size(line)?;
            expecting_offsets = false;
            line_in_block = 0;
            continue;
        }
        if line.is_empty() || line.starts_with(b"Consensus") {
            if in_first_block {
                expected_num_sequences = ids.len();
                in_first_block = false;
            } else if line_in_block != expected_num_sequences {
                return Err(ScanFailure::Failed);
            }
            expecting_offsets = true;
            continue;
        }
        if !in_first_block && line_in_block == expected_num_sequences {
            return Err(ScanFailure::Failed);
        }
        let (id, data) = process_data_line(line)?;
        if expected_size == 0 {
            expected_size = data.len();
        } else if data.len() != expected_size {
            return Err(ScanFailure::Failed);
        }
        if in_first_block {
            ids.push(id);
            sequences.push(data);
        } else {
            if ids.get(line_in_block).map(Vec::as_slice) != Some(id.as_slice()) {
                return Err(ScanFailure::Failed);
            }
            sequences[line_in_block].extend_from_slice(&data);
        }
        line_in_block += 1;
    }
    let info = SeqInfo {
        gap: b".".to_vec(),
        match_chars: Vec::new(),
        missing: Vec::new(),
    };
    for row in &sequences {
        if row
            .iter()
            .any(|c| !PROTEIN_ALPHABET.contains(c) && *c != b'.')
        {
            return Err(ScanFailure::Failed);
        }
    }
    Ok(Scanned {
        ids,
        rows: sequences,
        deflines: Vec::new(),
        info,
    })
}
