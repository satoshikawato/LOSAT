//! Mask helpers shared by the programs: the lower-case runs of a record's letters, and the mask
//! coordinate conversions between a DNA query and its translated frames (`dna_to_protein_masks`,
//! `protein_to_dna_masks`) and the merge of overlapping masks (`combine_masks`).
//!
//! `lowercase_runs` moved from `tblastn/archive.rs` and the translated-mask functions from
//! `blastx/query_setup.rs` (SW6 items 2 and 3) so that `algorithm/common`
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1604-1635
pub fn dna_to_protein_masks(masks: &[(i32, i32)], dna_length: usize, frame: i8) -> Vec<(i32, i32)> {
    let length = dna_length as i32;
    let frame = frame as i32;
    let aa_length = (length - (frame.abs() - 1)).max(0) / 3;
    masks
        .iter()
        .map(|&(left, right)| {
            let (from, to) = if frame < 0 {
                ((length + frame - right) / 3, (length + frame - left) / 3)
            } else {
                ((left - frame + 1) / 3, (right - frame + 1) / 3)
            };
            (from.max(0).min(aa_length - 1), to.max(0).min(aa_length - 1))
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:928-954
pub fn protein_to_dna_masks(masks: &[(i32, i32)], dna_length: usize, frame: i8) -> Vec<(i32, i32)> {
    let length = dna_length as i32;
    let frame = frame as i32;
    masks
        .iter()
        .map(|&(left, right)| {
            let (from, to) = if frame < 0 {
                (length - 3 * right + frame + 1, length - 3 * left + frame)
            } else {
                (3 * left + frame - 1, 3 * right + frame - 1)
            };
            (from.max(0).min(length - 1), to.max(0).min(length - 1))
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:983-1004
pub fn combine_masks(masks: Vec<(i32, i32)>) -> Vec<(i32, i32)> {
    combine_linked(masks, 0)
}

/// `BlastSeqLocCombine(mask_loc, link_value)` of inclusive `(left, right)` intervals: sorted by
/// their left ends, an interval joins the one before it when that one's right end plus
/// `link_value` is past its left end (`combine_masks` with 0; the repeat masks of blastn's
/// `-filtering_db` with `REPEAT_MASK_LINK_VALUE`, `blastn::repeat_filter`). NCBI's `qsort`
/// compares the left ends only; the result does not depend on the order of equal ones.
pub(crate) fn combine_linked<T>(mut masks: Vec<(T, T)>, link_value: T) -> Vec<(T, T)>
where
    T: Copy + Ord + std::ops::Add<Output = T>,
{
    masks.sort_by_key(|r| r.0);
    let mut result: Vec<(T, T)> = Vec::new();
    for mask in masks {
        if let Some(last) = result.last_mut().filter(|r| r.1 + link_value > mask.0) {
            last.1 = last.1.max(mask.1);
        } else {
            result.push(mask);
        }
    }
    result
}
