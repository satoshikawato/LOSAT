use super::*;

fn chunks(masks: &[(u32, u32)], length: usize) -> Vec<HardChunk> {
    let mut ranges = Vec::new();
    complement_masks(masks, length, &mut ranges);
    let mut split = HardSplitState::default();
    split.start(length, ranges);
    std::iter::from_fn(|| split.next_chunk(DBSEQ_CHUNK_OVERLAP)).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 579-584 (quoted
// at `HardSplitState::next_chunk`)
// The fields that the refusals read: the "no chunking" branch of a chunk from 0 before a mask
// that runs to the end (after empty masks at 0 as well, review F1), the range of each chunk,
// and the overlap read past the list after the last range's chunk (a short last range, and the
// last full chunk of a long one, review F3).
#[test]
fn hard_chunks_name_their_range_and_the_overlap_read_past_the_list() {
    let found = chunks(&[(0, 0), (0, 0), (1000, 1200)], 1200);
    assert_eq!(found.len(), 3);
    assert_eq!((found[2].offset, found[2].length), (0, 1000));
    assert!(found[2].full_range && found[2].range == 2 && !found[2].overlap_past_list);
    let found = chunks(&[(1001, 1002)], 1200);
    assert_eq!(
        found
            .iter()
            .map(|c| (c.offset, c.length, c.residual, c.range, c.overlap_past_list))
            .collect::<Vec<_>>(),
        [(0, 1001, 0, 0, false), (1000, 200, 2, 1, true)]
    );
    let found = chunks(&[(1001, 1002)], 5_300_000);
    assert_eq!(
        found
            .iter()
            .map(|c| (c.offset, c.length, c.range, c.overlap, c.overlap_past_list))
            .collect::<Vec<_>>(),
        [
            (0, 1001, 0, 100, false),
            (1000, 5_000_000, 1, 100, false),
            (5_000_900, 299_100, 1, 100, true),
        ]
    );
}
