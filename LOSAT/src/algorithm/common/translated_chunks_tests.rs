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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:286-307 (cited at
// `soft_chunk_scan_ranges` and `scan_lead`)
// The soft masks of probe M4b-T2 p1 (2026-11-23), stored out of order on 15,300,000 bases: frame
// +1 (5,100,000 residues) has two chunks. The second's ranges are NCBI's (its scan found the
// probe's copies there): an inverted one, and one from 3,899,900 residues before the chunk, which
// its scan reads from there.
#[test]
fn masks_stored_out_of_order_give_a_range_before_its_chunk() {
    let mut ranges = Vec::new();
    complement_masks(
        &[
            (14_999_400, 14_999_850),
            (14_999_880, 15_180_000),
            (3_000_000, 3_300_000),
        ],
        15_300_000,
        &mut ranges,
    );
    let mut split = HardSplitState::default();
    let chunks = translated_db_mask_chunks(&ranges, false, 1, 5_100_000, 5_100_000, &mut split)
        .unwrap()
        .into_iter()
        .map(|(chunk, scan)| {
            let lead = scan_lead(&scan);
            (chunk.offset, chunk.length, scan, lead)
        })
        .collect::<Vec<_>>();
    assert_eq!(
        chunks,
        vec![
            (
                0,
                5_000_000,
                vec![(0, 4_999_800), (4_999_950, 4_999_960)],
                0
            ),
            (
                4_999_900,
                100_100,
                vec![(50, 60), (60_100, -3_999_900), (-3_899_900, 100_100)],
                3_899_900
            ),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
// The letters around a chunk of 4 letters with 3 before it: offsets below 0 read the letters before
// the chunk, the left walks' subject starts at the chunk's first letter and the right walks' ends at
// its last.
#[test]
fn chunk_letters_read_before_the_chunk_and_walk_inside_it() {
    let letters = [10, 11, 12, 20, 21, 22, 23, 30];
    let chunk = ChunkLetters {
        letters: &letters,
        lead: 3,
        length: 4,
    };
    assert_eq!(
        (-3..5).map(|s| chunk.at(s)).collect::<Vec<_>>(),
        letters.to_vec()
    );
    assert_eq!(chunk.chunk_onward(), &[20, 21, 22, 23, 30]);
    assert_eq!(chunk.to_chunk_end(), &[10, 11, 12, 20, 21, 22, 23]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:145-215, 221-276, 805-830
#[test]
fn a_frame_after_an_empty_first_part_scans_the_whole_frame() {
    // [0, 1) of 1,200 bases: the plus frames' ranges are (0, 0), (0, 400). The empty part's chunk
    // scans nothing; the chunk from 0 after it takes the "no chunking" branch (NCBI then frees
    // memory that it did not allocate, PD-LOSAT-NCBI-DEFECTS exception 27) and scans the whole
    // frame: 400 residues in frame +1, the 399 of frame +2 (its chunk ends on the sentinel).
    assert_eq!(
        hard_chunks(&[(0, 1)], 1200, 1),
        vec![(0, 0, vec![(0, 0)]), (0, 400, vec![(0, 400)])]
    );
    assert_eq!(
        hard_chunks(&[(1, 2)], 1200, 2),
        vec![(0, 0, vec![(0, 0)]), (0, 400, vec![(0, 399)])]
    );
    // The same chunks without the mask, as in the near input (one chunk of the whole frame).
    assert_eq!(hard_chunks(&[], 1200, 1), vec![(0, 400, vec![(0, 400)])]);
    // The other strand keeps an empty last part: frame -1's one chunk is the whole frame.
    assert_eq!(
        hard_chunks(&[(0, 1)], 1200, -1),
        vec![(0, 400, vec![(0, 400)])]
    );
    // [1200, 1202) of 1,202 bases: the minus frames start with the empty part.
    assert_eq!(
        hard_chunks(&[(1200, 1202)], 1202, -1),
        vec![(0, 0, vec![(0, 0)]), (0, 400, vec![(0, 400)])]
    );
    // [0, 3) and [1200, 1201) of 1,201 bases: frame -1's chunk after the empty part ends at the
    // mask [0, 3) before the frame's end (399 of 400 residues) and scans its own letters
    // (exception 21).
    assert_eq!(
        hard_chunks(&[(0, 3), (1200, 1201)], 1201, -1),
        vec![(0, 0, vec![(0, 0)]), (0, 399, vec![(0, 399)])]
    );
    // An empty mask at 0 before a tail mask: frame +1's chunk after the empty part ends at the
    // tail mask (exception 21).
    assert_eq!(
        hard_chunks(&[(0, 0), (1100, 1200)], 1200, 1),
        vec![(0, 0, vec![(0, 0)]), (0, 366, vec![(0, 366)])]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-310,478-483,804-830,
// 1318-1331,1368-1369,1429
#[test]
fn translated_soft_mask_exits_match_chunks_and_skip_short_mask_reads() {
    let mut split = HardSplitState::default();
    let mut exits = Vec::new();

    // Subjects shorter than 3 * word + 2 bases produce no exits and do not read masks.
    let lead = translated_soft_mask_exits(10, 3, &mut exits, &mut split, || -> Result<_> {
        panic!("short subject must not read masks")
    })
    .unwrap();
    assert!(!lead);
    assert!(exits.is_empty());

    // Subjects whose frame +1 has at most MAX_DBSEQ_LEN residues emit the six frame lengths and
    // do not read masks.
    let lead =
        translated_soft_mask_exits(15_000_000, 3, &mut exits, &mut split, || -> Result<_> {
            panic!("unchunked subject must not read masks")
        })
        .unwrap();
    assert!(!lead);
    assert_eq!(
        exits,
        vec![5_000_000, 4_999_999, 4_999_999, 5_000_000, 4_999_999, 4_999_999]
    );

    // 15,300,000 bases with T2 out-of-order masks: exits equal translated_db_mask_chunks across
    // all six frames in NCBI's order and lead is true.
    let t2_masks = [
        (14_999_400, 14_999_850),
        (14_999_880, 15_180_000),
        (3_000_000, 3_300_000),
    ];
    let length = 15_300_000usize;
    let mut ranges = Vec::new();
    complement_masks(&t2_masks, length, &mut ranges);
    let mut expected_exits = Vec::new();
    let mut expected_lead = false;
    let mut check_split = HardSplitState::default();
    for frame in NCBI_FRAMES {
        let fl = frame_length(length, frame).unwrap();
        let calls =
            translated_db_mask_chunks(&ranges, false, frame, fl, length / 3, &mut check_split)
                .unwrap();
        for (chunk, scan) in calls {
            expected_exits.push(chunk.length);
            expected_lead |= scan_lead(&scan) > 0;
        }
    }
    let lead = translated_soft_mask_exits(length, 3, &mut exits, &mut split, || {
        Ok::<_, anyhow::Error>(Cow::Borrowed(ranges.as_slice()))
    })
    .unwrap();
    assert!(expected_lead);
    assert_eq!(lead, expected_lead);
    assert_eq!(exits, expected_exits);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276,805-830
#[test]
fn overlapping_hard_masks_record_negative_ncbi_length_and_full_range() {
    let mut ranges = Vec::new();
    complement_masks(&[(120, 600), (180, 714)], 720, &mut ranges);
    let mut split = HardSplitState::default();
    let plus1 = translated_db_mask_chunks(&ranges, true, 1, 240, 240, &mut split).unwrap();
    let summary: Vec<_> = plus1
        .iter()
        .map(|(c, _)| (c.offset, c.length, c.ncbi_length, c.full_range))
        .collect();
    assert_eq!(
        summary,
        vec![
            (0, 40, 40, false),
            (200, 0, -140, false),
            (238, 2, 2, false)
        ]
    );

    // Overlapping tail masks trigger `full_range` on frame +1.
    complement_masks(&[(120, 720), (180, 720)], 720, &mut ranges);
    let tail = translated_db_mask_chunks(&ranges, true, 1, 240, 240, &mut split).unwrap();
    assert_eq!(tail.len(), 1);
    assert!(tail[0].0.full_range);
    assert_eq!(tail[0].0.length, 40);
}
