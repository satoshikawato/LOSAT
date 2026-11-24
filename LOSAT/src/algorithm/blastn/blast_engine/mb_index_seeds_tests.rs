use super::*;

#[test]
fn a_chunks_number_counts_ncbis_chunk_offsets() {
    let step = MAX_DBSEQ_LEN - 100;
    assert_eq!(chunk_number(0, 100), 0);
    assert_eq!(chunk_number(step, 100), 1);
    assert_eq!(chunk_number(2 * step, 100), 2);
    // An overlap whose step is not a whole byte: each offset rounds down to one.
    let step = MAX_DBSEQ_LEN - 101;
    assert_eq!(chunk_number(step - step % 4, 101), 1);
    let next = step - step % 4 + step;
    assert_eq!(chunk_number(next - next % 4, 101), 2);
}

#[test]
fn the_diagonal_hash_keeps_ncbis_entries() {
    let mut hash = IrHash::new();
    let d1 = 0x1000_0000u32 + 5;
    let d2 = d1 + IR_HASH_SIZE;
    // An empty key gives its first entry.
    let e = hash.locate(d1);
    assert_eq!((e.diag, e.qend), (0, 0));
    e.diag = d1;
    e.qend = 10;
    // Another diagonal of the key: a new entry after the first.
    let e = hash.locate(d2);
    assert_eq!((e.diag, e.qend), (d2, 0));
    e.qend = 20;
    assert_eq!(hash.pool.len(), 1);
    // Found after the first: its values move to the first entry.
    let e = hash.locate(d2);
    assert_eq!((e.diag, e.qend), (d2, 20));
    let key = (d1 % IR_HASH_SIZE) as usize;
    assert_eq!((hash.heads[key].diag, hash.heads[key].qend), (d2, 20));
    assert_eq!((hash.pool[0].diag, hash.pool[0].qend), (d1, 10));
    let e = hash.locate(d1);
    assert_eq!((e.diag, e.qend), (d1, 10));
    assert_eq!(hash.pool.len(), 1);
    // A first entry whose end is 0 is taken by any diagonal of its key.
    let d3 = 0x1000_0000u32 + 6;
    let e = hash.locate(d3);
    assert_eq!((e.diag, e.qend), (0, 0));
    hash.clear();
    assert_eq!(hash.heads[key].next, 0);
    assert_eq!((hash.heads[key].diag, hash.heads[key].qend), (0, 0));
    assert!(hash.pool.is_empty() && hash.used.is_empty());
}

#[test]
fn the_lookup_segments_are_the_unmasked_intervals_of_valid_contexts() {
    let context = |query_offset, len, masks: &[(usize, usize)]| QueryContext {
        query_idx: 0,
        frame: 1,
        query_offset,
        len,
        masks: masks
            .iter()
            .map(|&(start, end)| crate::utils::dust::MaskedInterval { start, end })
            .collect(),
    };
    let contexts = [
        context(0, 50, &[(10, 20)]),
        context(51, 50, &[]),
        context(102, 30, &[(0, 30)]),
        context(133, 20, &[(15, 20)]),
    ];
    assert_eq!(
        lookup_segments(&contexts, &[true, false, true, true]),
        [(0, 10), (20, 50), (133, 148)]
    );
}

/// The seeds that `seeds` (left end in the sequence, query offset) give the hard chunks of a
/// sequence of `length` letters with the hard ranges `ranges`: each chunk's offset, residual,
/// length and seeds (query offset, offset in the chunk).
fn placed(
    length: usize,
    ranges: &[(i32, i32)],
    seeds: &[(usize, u32)],
) -> Vec<(usize, usize, usize, Vec<(u32, u32)>)> {
    let mut hard = HardSeeds::chunks_of(length, ranges, 100, 0);
    for &(start, q_off) in seeds {
        hard.place(Seed { q_off, s_off: 0 }, start, ranges);
    }
    hard.chunks
        .iter()
        .map(|chunk| {
            (
                chunk.offset,
                chunk.residual,
                chunk.length,
                chunk
                    .seeds
                    .iter()
                    .map(|seed| (seed.q_off, seed.s_off))
                    .collect(),
            )
        })
        .collect()
}

/// Each seed goes to the range that holds its left end, in the order given; a range after a mask
/// starts a chunk at a whole byte before it; a seed at a masked letter is not placed.
#[test]
fn a_seed_goes_to_the_range_that_holds_it() {
    // Masks [100, 150) and [210, 230).
    let ranges = [(0, 100), (150, 210), (230, 1000)];
    let got = placed(
        1000,
        &ranges,
        &[(60, 5), (180, 9), (110, 0), (150, 3), (999, 1), (40, 2)],
    );
    assert_eq!(
        got,
        vec![
            (0, 0, 100, vec![(5, 60), (2, 40)]),
            (148, 2, 62, vec![(9, 32), (3, 2)]),
            (228, 2, 772, vec![(1, 771)]),
        ]
    );
}

/// Each seed of a range over `MAX_DBSEQ_LEN` letters goes to the last chunk that starts at or
/// before it, its overlap with the chunk before included, as the index places it.
#[test]
fn a_long_range_places_a_seed_in_the_later_chunk() {
    let step = MAX_DBSEQ_LEN - 100;
    let length = MAX_DBSEQ_LEN + 1000;
    let ranges = [(0, 40), (64, length as i32)];
    let got = placed(length, &ranges, &[(64 + step - 10, 0), (64 + step + 10, 1)]);
    assert_eq!(
        got,
        vec![
            (0, 0, 40, vec![]),
            (64, 0, MAX_DBSEQ_LEN, vec![(0, (step - 10) as u32)]),
            (64 + step, 0, length - 64 - step, vec![(1, 10)]),
        ]
    );
}

/// The chunk loop's chunk finds its seeds by offset, residual and length; chunks of hard ranges of
/// no letters at one place hold none.
#[test]
fn a_hard_chunk_finds_its_seeds() {
    // Masks [0, 9), [9, 9) and [9, 20): two empty ranges at 9.
    let ranges = [(0, 0), (9, 9), (9, 9), (20, 300)];
    let mut hard = HardSeeds::chunks_of(300, &ranges, 100, 3);
    hard.word_size = 28;
    hard.place(Seed { q_off: 7, s_off: 0 }, 30, &ranges);
    let offsets: Vec<_> = hard
        .chunks
        .iter()
        .map(|chunk| (chunk.offset, chunk.residual, chunk.length))
        .collect();
    assert_eq!(offsets, vec![(0, 0, 0), (8, 1, 1), (8, 1, 1), (20, 0, 280)]);
    assert!(hard.chunk(8, 1, 1).is_none());
    let (seeds, word_size, at) = hard.chunk(20, 0, 280).unwrap();
    assert_eq!((word_size, at), (28, 3));
    assert_eq!(
        seeds,
        &[Seed {
            q_off: 7,
            s_off: 10
        }][..]
    );
}
