//! NCBI's FASTA-gap scanner (`CAlnScannerFastaGap`): a `>id defline` line, then the data lines
//! of the row (blanks inside a line are dropped), a blank line or the next `>` ends the row.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_fastagap.cpp:66-176

use super::base::{
    existing_id, strip_blanks, verify_rows, IdMatch, Scan, ScanFailure, Scanned, SeqInfo,
};
use super::lines::PeekAheadStream;
use crate::ncbi_str::truncate_spaces;
use crate::utils::format_guess::split_in_two;

/// `xImportAlignmentData` and the base verification. A data line outside a row is read as a
/// defline too (its first character is cut off: `sSplitFastaDef` does `substr(1)`).
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Scan<Scanned> {
    let mut waiting_for_data = false;
    let mut expected_sizes: Vec<usize> = Vec::new();
    let mut data_line_index = 0usize;
    let mut first_sequence = true;
    let mut expected_num_lines = 0usize;
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut deflines: Vec<Vec<u8>> = Vec::new();
    let mut rows: Vec<Vec<u8>> = Vec::new();
    while let Some((raw, _line_number)) = input.read_line() {
        let line = truncate_spaces(raw);
        if waiting_for_data {
            if line.is_empty() || line.starts_with(b">") {
                waiting_for_data = false;
                first_sequence = false;
            } else {
                let data = strip_blanks(line);
                if first_sequence {
                    expected_sizes.push(data.len());
                } else {
                    if expected_num_lines == 0 {
                        expected_num_lines = expected_sizes.len();
                    }
                    let expected = if data_line_index < expected_num_lines {
                        expected_sizes[data_line_index]
                    } else {
                        0
                    };
                    if data.len() != expected {
                        return Err(ScanFailure::Failed);
                    }
                }
                rows.last_mut()
                    .expect("a row is open while waiting for its data")
                    .extend_from_slice(&data);
                data_line_index += 1;
                continue;
            }
        }
        if line.is_empty() {
            continue;
        }
        let def = truncate_spaces(&line[1..]);
        let (id, defline, _) = split_in_two(def, b" \t", true);
        if id.is_empty() || existing_id(&ids, id) != IdMatch::Different {
            return Err(ScanFailure::Failed);
        }
        ids.push(id.to_vec());
        deflines.push(defline.to_vec());
        rows.push(Vec::new());
        waiting_for_data = true;
        data_line_index = 0;
    }
    let info = SeqInfo::protein();
    verify_rows(&info, &rows)?;
    Ok(Scanned {
        ids,
        rows,
        deflines,
        info,
    })
}
