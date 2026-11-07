//! Mask helpers shared by the programs: the lower-case runs of a record's letters.
//!
//! `lowercase_runs` moved from `tblastn/archive.rs` (SW6 item 2) so that `algorithm/common`
//! depends on no program.

/// The runs of lower-case letters of a record (inclusive 0-based positions).
pub(crate) fn lowercase_runs(letters: &[u8]) -> Vec<(i32, i32)> {
    let mut runs = Vec::new();
    let mut start = None;
    for (index, letter) in letters.iter().enumerate() {
        if letter.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(from) = start.take() {
            runs.push((from as i32, index as i32 - 1));
        }
    }
    if let Some(from) = start {
        runs.push((from as i32, letters.len() as i32 - 1));
    }
    runs
}
