//! NCBI's Sequin scanner (`CAlnScannerSequin`): blocks of data lines (an id and chunks of
//! letters, or `id > start end chunks..` ) separated by blank lines and offsets lines, ended by a
//! `//` line; the match character `.` stands for the letter of the block's first row.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_sequin.cpp:50-259

use super::base::{
    existing_id, split_merge, verify_rows, IdMatch, Scan, ScanFailure, Scanned, SeqInfo,
};
use super::lines::PeekAheadStream;
use crate::ncbi_str::truncate_spaces;

/// `xIsSequinOffsetsLine`: at most five tokens (split at spaces only), each ending in `0`.
fn is_offsets_line(line: &[u8]) -> bool {
    let tokens = split_merge(line, b" ");
    tokens.len() <= 5 && tokens.iter().all(|token| token.ends_with(b"0"))
}

/// `xExtractSequinSequenceData`: the id and the letters of a data line, or `None` for a line with
/// fewer than two tokens (or fewer than five when the second token is `>`).
fn extract_data(line: &[u8]) -> Option<(Vec<u8>, Vec<u8>)> {
    let tokens = split_merge(line, b" ");
    if tokens.len() < 2 {
        return None;
    }
    let id = tokens[0].to_vec();
    let mut data = Vec::new();
    if tokens[1] == b">" {
        if tokens.len() < 5 {
            return None;
        }
        for token in &tokens[3..tokens.len() - 1] {
            data.extend_from_slice(token);
        }
    } else {
        for token in &tokens[1..] {
            data.extend_from_slice(token);
        }
    }
    Some((id, data))
}

/// `xImportAlignmentData`, `xAdjustSequenceInfo` (gap `-`, match `.`) and the base verification.
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Scan<Scanned> {
    let mut processing_data = true;
    let mut in_first_block = false;
    let mut line_in_block = 0usize;
    let mut ref_data: Vec<u8> = Vec::new();
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    while let Some((raw, _)) = input.read_line() {
        let line = truncate_spaces(raw);
        if line.is_empty() {
            if processing_data {
                if line_in_block != ids.len() {
                    return Err(ScanFailure::Failed);
                }
                processing_data = false;
                if !ids.is_empty() {
                    in_first_block = false;
                }
            }
            continue;
        }
        if is_offsets_line(line) {
            if processing_data {
                return Err(ScanFailure::Failed);
            }
            continue;
        }
        if line == b"//" {
            if processing_data {
                return Err(ScanFailure::Failed);
            }
            break;
        }
        if !processing_data {
            processing_data = true;
            line_in_block = 0;
            ref_data.clear();
            if ids.is_empty() {
                in_first_block = true;
            }
        }
        let Some((id, mut data)) = extract_data(line) else {
            return Err(ScanFailure::Failed);
        };
        if in_first_block {
            if existing_id(&ids, &id) != IdMatch::Different {
                return Err(ScanFailure::Failed);
            }
            ids.push(id);
            sequences.push(Vec::new());
        } else {
            if line_in_block >= ids.len() {
                return Err(ScanFailure::Failed);
            }
            if ids[line_in_block] != id {
                return Err(ScanFailure::Failed);
            }
        }
        if ref_data.is_empty() {
            ref_data = data.clone();
        } else {
            if data.len() != ref_data.len() {
                return Err(ScanFailure::Failed);
            }
            for (letter, reference) in data.iter_mut().zip(&ref_data) {
                if *letter == b'.' {
                    *letter = *reference;
                }
            }
        }
        sequences[line_in_block].extend_from_slice(&data);
        line_in_block += 1;
    }
    let info = SeqInfo {
        gap: b"-".to_vec(),
        match_chars: b".".to_vec(),
        missing: Vec::new(),
    };
    verify_rows(&info, &sequences)?;
    Ok(Scanned {
        ids,
        rows: sequences,
        deflines: Vec::new(),
        info,
    })
}
