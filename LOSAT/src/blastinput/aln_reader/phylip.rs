//! NCBI's PHYLIP scanner (`CAlnScannerPhylip`): a first line with the sequence count and the
//! length, then blocks of one data line per sequence (the first block starts each line with the
//! id), blank lines skipped, `>`/`[` lines kept as deflines, and the lines `begin ncbi;`,
//! `sequin` and `end;` ignored.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_phylip.cpp:80-232

use super::base::{
    existing_id, process_defline, split_merge, strip_blanks, verify_rows, IdMatch, Scan,
    ScanFailure, Scanned, SeqInfo,
};
use super::lines::PeekAheadStream;
use crate::ncbi_str::{string_to_int, truncate_spaces};
use crate::utils::format_guess::split_in_two;

/// The cause of the refusal of a count of 0 (the remainder of the data line count by the sequence count).
pub(crate) const ZERO_COUNT: &str =
    "NCBI BLAST+ 2.17.0 stops with an integer division by zero in the alignment reader";

/// `xImportAlignmentData`, `xVerifyAlignmentData` and the base verification.
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Scan<Scanned> {
    // The first line is the guessed one: two numbers.
    let Some((first, _)) = input.read_line() else {
        return Err(ScanFailure::Failed);
    };
    let first = truncate_spaces(first);
    let tokens = split_merge(first, b" \t");
    let count = string_to_int(tokens[0], 0).map_err(ScanFailure::Exception)?;
    let length = string_to_int(tokens[1], 0).map_err(ScanFailure::Exception)?;
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut deflines: Vec<Vec<u8>> = Vec::new();
    // The sequence of the first row, to its length: `mSequences[0]` summed over its data lines.
    let mut data_line_count = 0usize;
    let mut block_line_length = 0usize;
    // The sequence count is an `int`, and the remainder of the data line count by it needs it non-zero.
    let count_u = count as usize;
    while let Some((raw, _)) = input.read_line() {
        let line = truncate_spaces(raw);
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' || line[0] == b'[' {
            let (dummy, def) = process_defline(line)?;
            if !dummy.is_empty() {
                return Err(ScanFailure::Failed);
            }
            deflines.push(def);
            continue;
        }
        let lower = line.to_ascii_lowercase();
        if matches!(
            lower.as_slice(),
            b"begin ncbi;" | b"begin ncbi" | b"sequin" | b"end;" | b"end"
        ) {
            continue;
        }
        if count == 0 {
            return Err(ScanFailure::Undefined(ZERO_COUNT));
        }
        let new_block = data_line_count.is_multiple_of(count_u);
        let data: Vec<u8> = if data_line_count < count_u {
            let (id, rest, _) = split_in_two(line, b" \t", true);
            if rest.is_empty() {
                return Err(ScanFailure::Failed);
            }
            if existing_id(&ids, id) != IdMatch::Different {
                return Err(ScanFailure::Failed);
            }
            ids.push(id.to_vec());
            sequences.push(Vec::new());
            strip_blanks(rest)
        } else {
            strip_blanks(line)
        };
        if new_block {
            block_line_length = data.len();
        } else if data.len() != block_line_length {
            return Err(ScanFailure::Failed);
        }
        sequences[data_line_count % count_u].extend_from_slice(&data);
        data_line_count += 1;
    }
    if count == 0 {
        return Err(ScanFailure::Undefined(ZERO_COUNT));
    }
    if !data_line_count.is_multiple_of(count_u) {
        return Err(ScanFailure::Failed);
    }
    // `xVerifyAlignmentData`.
    if count_u != ids.len() {
        return Err(ScanFailure::Failed);
    }
    if length as usize != sequences[0].len() {
        return Err(ScanFailure::Failed);
    }
    let info = SeqInfo::protein();
    verify_rows(&info, &sequences)?;
    Ok(Scanned {
        ids,
        rows: sequences,
        deflines,
        info,
    })
}
