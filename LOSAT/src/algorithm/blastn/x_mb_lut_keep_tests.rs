use super::super::x_mb_lookup_par::tests::{random_case, Lcg};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::super::x_mb_narrow;
use super::super::{build_mb_lookup, MbLookupTable};
use super::*;

/// The kept table's length and the number of its first elements that may be dirty.
fn kept(keep: &Keep) -> Option<(usize, usize)> {
    let kept = keep.kept.lock().unwrap();
    kept.as_ref().map(|kept| (kept.table.len(), kept.dirty))
}

/// The build's table (its first cells) and the other arrays and fields equal the fresh
/// build's (both of 3-byte elements, LOSAT_X_MB24, or both not).
fn assert_same(
    built: &MbLookupTable,
    fresh: &MbLookupTable,
    cells: usize,
    next_pos_len: usize,
    what: &str,
) {
    assert_eq!(built.x_narrow.is_some(), fresh.x_narrow.is_some(), "{what}");
    let (built_heads, built_next) = built.x_wide(cells, next_pos_len);
    let (fresh_heads, fresh_next) = fresh.x_wide(cells, next_pos_len);
    assert!(built_heads == fresh_heads, "{what}: hashtable");
    assert!(built_next == fresh_next, "{what}: next_pos");
    assert!(built.pv_array == fresh.pv_array, "{what}: pv_array");
    assert_eq!(built.pv_array_bts, fresh.pv_array_bts, "{what}");
    assert_eq!(built.longest_chain, fresh.longest_chain, "{what}");
}

// A keep takes back the longest table with the start that builds wrote, gives it out to
// builds of no more cells with their cells set to 0 (the written elements past them stay and
// stay counted), and frees it for a longer build.
#[test]
fn the_keep_holds_the_longest_table_and_gives_it_out_zeroed() {
    let keep = Keep::default();
    let (table, from_keep, dirty) = keep.take(8, false);
    assert_eq!((from_keep, dirty), (false, 8));
    assert_eq!(table, vec![0; 8]);
    keep.give_back(vec![7; 8], 8);
    keep.give_back(vec![5; 4], 4);
    assert_eq!(kept(&keep), Some((8, 8)), "a shorter table is freed");
    let (mut table, from_keep, dirty) = keep.take(4, false);
    assert!(from_keep, "a shorter build takes the kept table");
    assert_eq!(
        dirty, 8,
        "the written elements past the build's cells stay counted"
    );
    assert_eq!(
        table,
        [0, 0, 0, 0, 7, 7, 7, 7],
        "only the build's cells are set to 0"
    );
    assert_eq!(kept(&keep), None);
    table[..2].fill(3);
    keep.give_back(table, dirty);
    assert_eq!(kept(&keep), Some((8, 8)));
    let (table, from_keep, dirty) = keep.take(8, false);
    assert_eq!((from_keep, dirty), (true, 8));
    assert_eq!(table, vec![0; 8]);
    keep.give_back(table, 4);
    assert_eq!(kept(&keep), Some((8, 4)));
    let (table, from_keep, dirty) = keep.take(2, false);
    assert_eq!((from_keep, dirty), (true, 4));
    keep.give_back(table, dirty);
    let (table, from_keep, _) = keep.take(16, false);
    assert!(!from_keep, "a longer build allocates");
    assert_eq!(table, vec![0; 16]);
    assert_eq!(kept(&keep), None, "the shorter kept table is freed first");
    keep.give_back(vec![9; 1 << 20], 1 << 20);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let (table, from_keep, dirty) = pool.install(|| keep.take(1 << 12, true));
        assert_eq!((from_keep, dirty), (true, 1 << 20));
        assert!(
            table[..1 << 12].iter().all(|&cell| cell == 0),
            "zeroed on the pool"
        );
        assert!(
            table[1 << 12..].iter().all(|&cell| cell == 9),
            "past the build's cells"
        );
    }
}

// After the loop's last preliminary stage the keep frees its table (`release`); a release
// of an empty keep frees nothing, and a later build allocates its table as the first build.
#[test]
fn a_release_frees_the_kept_table() {
    let keep = Keep::default();
    assert!(!keep.release(), "nothing kept");
    let (table, from_keep, dirty) = keep.take(8, false);
    assert!(!from_keep);
    keep.give_back(table, dirty);
    assert_eq!(kept(&keep), Some((8, 8)));
    assert!(keep.release(), "the kept table is freed");
    assert_eq!(kept(&keep), None);
    let (table, from_keep, dirty) = keep.take(8, false);
    assert_eq!(
        (from_keep, dirty),
        (false, 8),
        "a fresh table after a release"
    );
    assert_eq!(table, vec![0; 8]);
}

// LOSAT_X_MBCSR: the batches of a search can take either side of its size rule at 4^12 cells, so
// a chains build can follow a CSR build on one keep; it is the fresh build.
#[test]
fn chains_after_a_csr_of_4_12_cells_are_the_fresh_build() {
    let mut rng = Lcg(0x6b65_6571);
    let keep = Arc::new(Keep::default());
    let (queries, offsets, masks, word_length, _) = random_case(&mut rng);
    let word_length = word_length.max(12);
    let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let build = |keep: Option<&Arc<Keep>>, letters: Option<u64>| {
        build_mb_lookup(
            &refs,
            &offsets,
            word_length,
            12,
            &masks,
            None,
            0,
            1_000_000,
            false,
            keep,
            letters,
        )
    };
    drop(build(Some(&keep), None));
    let csr = build(Some(&keep), Some(u64::MAX));
    assert_eq!(
        csr.x_csr.is_some(),
        super::super::x_bn_csr::x_will_pack(false) && super::super::x_bn_csr::x_mb12_on()
    );
    drop(csr);
    let fresh = build(None, None);
    let built = build(Some(&keep), None);
    let next_pos_len = queries
        .iter()
        .zip(&offsets)
        .map(|(q, &o)| o.max(0) as usize + q.len())
        .max()
        .unwrap_or(0)
        + 1;
    assert_same(&built, &fresh, 1 << 24, next_pos_len, "chains after a CSR");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1262
// Random queries (ambiguity codes, masks, both cell orders) of tables of 4^4 to 4^12 cells in
// a row through one keep, one after another and on pools of 2 and 5 threads (the pool fill
// with the `resident` flag; these queries are too short for the prefault round, which the
// whole-search tests cover): every build equals the same build on a fresh table, every build
// of no more cells than the kept table takes it (filled by earlier builds), and the elements
// past the start that builds wrote stay 0. With LOSAT_X_MB24 the newest-first tables of 4^12
// cells take the words of their 3-byte heads, fewer than 4^12. With LOSAT_X_BNCSRPACK a
// newest-first table of 4^8 to 4^11 cells is CSR at any size, so those tables are built in
// reading order here, which keeps their chains (the switch off builds them newest first too).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn builds_on_a_kept_table_are_the_fresh_builds() {
    const LUTS: [usize; 8] = [11, 8, 11, 11, 4, 12, 9, 12];
    let mut rng = Lcg(0x6b65_6570);
    for threads in [1usize, 2, 5] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let keep = Arc::new(Keep::default());
            let (mut reused, mut reused_dirty) = (0, 0);
            for case in 0..LUTS.len() * 2 {
                let (queries, offsets, masks, word_length, _) = random_case(&mut rng);
                let lut = LUTS[case % LUTS.len()];
                let word_length = word_length.max(lut);
                let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                let ascending = case % 3 == 1
                    || (super::super::x_bn_csr::x_will_pack(false) && (8..=11).contains(&lut));
                let fresh = build_mb_lookup(
                    &refs,
                    &offsets,
                    word_length,
                    lut,
                    &masks,
                    None,
                    0,
                    1_000_000,
                    ascending,
                    None,
                    None,
                );
                let next_pos_len = queries
                    .iter()
                    .zip(&offsets)
                    .map(|(q, &o)| o.max(0) as usize + q.len())
                    .max()
                    .unwrap_or(0)
                    + 1;
                let cells = if fresh.x_narrow.is_some() {
                    x_mb_narrow::words_for(1 << (2 * lut))
                } else {
                    1 << (2 * lut)
                };
                let before = keep.kept.lock().unwrap().as_ref().map(|kept| {
                    assert!(
                        kept.table[kept.dirty..].iter().all(|&cell| cell == 0),
                        "case {case}: an element past the dirty cells is not 0"
                    );
                    (
                        kept.table.len(),
                        kept.dirty,
                        kept.table[..kept.dirty.min(cells)].iter().any(|&c| c != 0),
                    )
                });
                let built = build_mb_lookup(
                    &refs,
                    &offsets,
                    word_length,
                    lut,
                    &masks,
                    None,
                    0,
                    1_000_000,
                    ascending,
                    Some(&keep),
                    None,
                );
                let dirty = match before {
                    Some((len, dirty, written)) if len >= cells => {
                        assert_eq!(kept(&keep), None, "case {case}: the kept table taken");
                        assert_eq!(built.x_table().len(), len, "case {case}");
                        reused += 1;
                        reused_dirty += usize::from(written);
                        dirty.max(cells)
                    }
                    _ => {
                        assert_eq!(built.x_table().len(), cells, "case {case}: fresh");
                        cells
                    }
                };
                assert_same(
                    &built,
                    &fresh,
                    1 << (2 * lut),
                    next_pos_len,
                    &format!("case {case} lut {lut} threads {threads}"),
                );
                let len = built.x_table().len();
                drop(built);
                assert_eq!(kept(&keep), Some((len, dirty)), "case {case}: given back");
            }
            // Every build but the first (4^11) and the first 4^12, which is longer than the
            // kept table: 14 of 16. With LOSAT_X_MB24 the first 4^12 (case 5) is a 3-byte
            // table, shorter than a 4-byte one, so the next 4^12 (case 7, ascending: 4-byte)
            // is longer than the kept table too: 13.
            let expected = if x_mb_narrow::on() { 13 } else { 14 };
            assert_eq!(
                reused, expected,
                "threads {threads}: builds on a kept table"
            );
            assert!(
                reused_dirty >= 7,
                "threads {threads}: {reused_dirty} filled tables kept"
            );
        });
    }
}
