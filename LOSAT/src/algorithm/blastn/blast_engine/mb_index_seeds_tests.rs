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
