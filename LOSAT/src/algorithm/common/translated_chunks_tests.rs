use super::*;
use crate::core::subject_split::{complement_masks, HardSplitState};

/// The hard chunks of frame `frame` of a sequence of `length` bases with the stored masks `masks`:
/// (offset, length, scan ranges) per chunk.
fn hard_chunks(
    masks: &[(u32, u32)],
    length: usize,
    frame: i8,
) -> Vec<(usize, usize, Vec<(i32, i32)>)> {
    let mut ranges = Vec::new();
    complement_masks(masks, length, &mut ranges);
    let frame_length = (length - (frame.unsigned_abs() as usize - 1)) / 3;
    let mut split = HardSplitState::default();
    translated_db_mask_chunks(&ranges, true, frame, frame_length, length / 3, &mut split)
        .unwrap()
        .into_iter()
        .map(|(chunk, scan)| (chunk.offset, chunk.length, scan))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:150-152, 162, 221-311
#[test]
fn a_whole_frame_scans_to_the_frame_s_end_and_a_tail_masked_one_its_chunk() {
    // No mask, 1,200 bases: frame +2 has 399 residues; its one range maps to [0, 400), a chunk
    // that ends on the frame's sentinel and takes the "no chunking" branch, which scans the
    // frame's letters only.
    assert_eq!(hard_chunks(&[], 1200, 2), vec![(0, 400, vec![(0, 399)])]);
    assert_eq!(hard_chunks(&[], 1200, 1), vec![(0, 400, vec![(0, 400)])]);
    // A tail mask [1000, 1200): frame +1's chunk from 0 takes the same branch while it ends at the
    // mask (333 of 400 residues): it scans its own letters (PD-LOSAT-NCBI-DEFECTS exception 21).
    assert_eq!(
        hard_chunks(&[(1000, 1200)], 1200, 1),
        vec![(0, 333, vec![(0, 333)])]
    );
    // Frame +2 (399 residues): the range after the mask starts at 400, past the frame; the chunk
    // from 0 is a normal one and scans its letters.
    assert_eq!(
        hard_chunks(&[(1000, 1200)], 1200, 2),
        vec![(0, 333, vec![(0, 333)])]
    );
    // The near input, the mask 4 bases shorter: a normal chunk, then one residue after the mask.
    assert_eq!(
        hard_chunks(&[(1000, 1196)], 1200, 1),
        vec![(0, 333, vec![(0, 333)]), (398, 2, vec![(0, 2)])]
    );
    // A head mask [0, 100) in frame -1 (400 residues): its first part, from 0, ends at the mask.
    assert_eq!(
        hard_chunks(&[(0, 100)], 1200, -1),
        vec![(0, 367, vec![(0, 367)])]
    );
}
