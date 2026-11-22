use super::super::MaskedInterval;
use super::*;

fn words<'a>(
    refs: &'a [&'a [u8]],
    offsets: &'a [i32],
    masks: &'a [Vec<MaskedInterval>],
    lut: usize,
) -> Words<'a> {
    Words {
        queries_blastna: refs,
        query_offsets: offsets,
        query_masks: masks,
        word_length: 16,
        lut_word_length: lut,
        kmer_mask: (1u64 << (2 * lut)) - 1,
        ascending_cells: false,
    }
}

// The scope: 4^12 cells and `lut_word_length` 12, newest first, no database word counts or
// BLEMIR_DEBUG, `next_pos` of at most 2^24 elements, query offsets that ascend without overlap
// (queries shorter than a word do not count) and no inverted mask.
#[test]
fn the_scope_takes_ascending_queries_of_4_12_tables_only() {
    let queries: Vec<Vec<u8>> = vec![vec![0; 100], vec![1; 5], vec![2; 50]];
    let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let plain = vec![Vec::new(); 3];
    let ok = words(&refs, &[0, 101, 107], &plain, 12);
    assert!(in_scope(&ok, CELLS, false, false, 158));
    assert!(in_scope(&ok, CELLS, false, false, CELLS));
    assert!(!in_scope(&ok, CELLS, false, false, CELLS + 1));
    assert!(!in_scope(&ok, CELLS, true, false, 158));
    assert!(!in_scope(&ok, CELLS, false, true, 158));
    let lut11 = words(&refs, &[0, 101, 107], &plain, 11);
    assert!(!in_scope(&lut11, 1 << 22, false, false, 158));
    let ascending = Words {
        ascending_cells: true,
        ..words(&refs, &[0, 101, 107], &plain, 12)
    };
    assert!(!in_scope(&ascending, CELLS, false, false, 158));
    // a query shorter than a word does not count; overlapping or descending ones do
    let short_first = words(&refs, &[0, 0, 101], &plain, 12);
    assert!(in_scope(&short_first, CELLS, false, false, 152));
    let overlapping = words(&refs, &[0, 101, 99], &plain, 12);
    assert!(!in_scope(&overlapping, CELLS, false, false, 150));
    let descending = words(&refs, &[60, 0, 0], &plain, 12);
    assert!(!in_scope(&descending, CELLS, false, false, 161));
    let inverted = vec![
        vec![MaskedInterval { start: 70, end: 20 }],
        Vec::new(),
        Vec::new(),
    ];
    let inverted = words(&refs, &[0, 101, 107], &inverted, 12);
    assert!(!in_scope(&inverted, CELLS, false, false, 158));
}

#[test]
fn the_words_hold_the_elements_and_the_last_read() {
    for n in 0..64 {
        let w = words_for(n);
        assert!(4 * w >= 3 * n + 1, "{n}");
        assert!(capacity(w) >= n, "{n}");
        assert!(words_for(capacity(w)) <= w, "{n}");
    }
    assert_eq!(words_for(CELLS), 12_582_913);
}

#[test]
fn every_store_writes_its_element_only() {
    let n = 9;
    let values = [
        0,
        1,
        MASK,
        0x00AB_CDEF,
        0x0001_0000,
        0x0000_FFFF,
        7,
        MASK - 1,
        0x0080_0000,
    ];
    // set (and the exact store) keep the neighbours; set_spill writes in increasing order
    let mut a = vec![0u32; words_for(n)];
    let mut b = vec![0u32; words_for(n)];
    let mut c = vec![0u32; words_for(n)];
    for (i, &v) in values.iter().enumerate().rev() {
        set(&mut a, i, v);
        // SAFETY: element `i` lies in `b`, which no other thread uses.
        unsafe { store_exact(b.as_mut_ptr(), i, v) };
    }
    for (i, &v) in values.iter().enumerate() {
        set_spill(&mut c, i, v);
    }
    for (i, &v) in values.iter().enumerate() {
        assert_eq!(get(&a, i), v);
        assert_eq!(get(&b, i), v);
        assert_eq!(get(&c, i), v);
        // SAFETY: as above.
        unsafe {
            assert_eq!(load_exact(b.as_ptr(), i), v);
            assert_eq!(load_word(b.as_ptr(), i) & MASK, v);
        }
    }
    assert_eq!(a, b);
    assert_eq!(a, c);
    // a merged store keeps the next element's byte
    // SAFETY: as above.
    unsafe {
        let w = load_word(a.as_ptr(), 3);
        store_word(a.as_mut_ptr(), 3, (w & !MASK) | 5);
    }
    assert_eq!(get(&a, 3), 5);
    assert_eq!(get(&a, 4), values[4]);
    assert_eq!(widen(&a, n)[4..], values[4..]);
}
