use super::*;

fn chunks(masks: &[(u32, u32)], length: usize) -> Vec<HardChunk> {
    chunks_of(masks, length, true)
}

fn chunks_of(masks: &[(u32, u32)], length: usize, nucleotide: bool) -> Vec<HardChunk> {
    let mut ranges = Vec::new();
    complement_masks(masks, length, &mut ranges);
    let mut split = HardSplitState::default();
    split.start(length, ranges, nucleotide);
    std::iter::from_fn(|| split.next_chunk(DBSEQ_CHUNK_OVERLAP)).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 579-584 (cited
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 578-584 (cited
// at `HardSplitState::next_chunk`)
// A protein program's chunks start at their hard range's left (`residual` 0), so the overlap
// that the merge after a chunk reads is the next range's left against that offset; the first
// mask running to the end takes the "no chunking" branch for the chunk before it.
#[test]
fn protein_chunks_start_at_their_range() {
    let found = chunks_of(&[(150, 190)], 400, false);
    assert_eq!(
        found
            .iter()
            .map(|c| (c.offset, c.length, c.residual, c.range, c.overlap))
            .collect::<Vec<_>>(),
        [(0, 150, 0, 0, 100), (190, 210, 0, 1, 100)]
    );
    assert!(!found[0].full_range && !found[1].full_range);
    let found = chunks_of(&[(300, 400)], 400, false);
    assert_eq!(found.len(), 1);
    assert!(found[0].full_range && found[0].length == 300);
    let found = chunks_of(&[(0, 400)], 400, false);
    assert_eq!(found.len(), 1);
    assert!(found[0].full_range && found[0].length == 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276 (cited at
// `HardSplitState::next_chunk`)
// Masks that overlap or nest, or are out of the order of their start (series NP, M4b-3b): the range
// between them runs back from the end of one mask to the start of the next, and its chunk's length
// is NCBI's `Int4` difference, negative (no letters); touching masks give a chunk of length 0. The
// chunks stop once the next range starts at the end, also where a later range would start before
// it.
#[test]
fn protein_chunks_of_overlapping_masks_have_signed_lengths() {
    let of = |masks: &[(u32, u32)], length: usize| {
        chunks_of(masks, length, false)
            .iter()
            .map(|c| (c.offset, c.length, c.letters()))
            .collect::<Vec<_>>()
    };
    assert_eq!(
        of(&[(40, 60), (55, 70)], 240),
        [(0, 40, 40), (60, -5, 0), (70, 170, 170)]
    );
    assert_eq!(
        of(&[(40, 100), (50, 60)], 240),
        [(0, 40, 40), (100, -50, 0), (60, 180, 180)]
    );
    assert_eq!(
        of(&[(40, 60), (60, 70)], 240),
        [(0, 40, 40), (60, 0, 0), (70, 170, 170)]
    );
    assert_eq!(of(&[(40, 240), (60, 100)], 240), [(0, 40, 40)]);
}
