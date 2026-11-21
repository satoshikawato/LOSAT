use super::super::build_two_stage_lookup;
use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
// Random queries (with ambiguity codes and masks) through the original loop and through
// the partitioned fill with 1..8 workers: every table element must be the same.
pub(crate) struct Lcg(pub(crate) u64);
impl Lcg {
    pub(crate) fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    pub(crate) fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

pub(crate) fn random_case(
    rng: &mut Lcg,
) -> (
    Vec<Vec<u8>>,
    Vec<i32>,
    Vec<Vec<MaskedInterval>>,
    usize,
    usize,
) {
    let lut_word_length = [4usize, 8, 9, 11, 12][rng.below(5) as usize];
    let word_length = lut_word_length + rng.below(17) as usize;
    let n_queries = 1 + rng.below(3) as usize;
    let mut queries = Vec::new();
    let mut offsets = Vec::new();
    let mut masks = Vec::new();
    let mut offset = 0i32;
    for _ in 0..n_queries {
        let len = 1 + rng.below(3000) as usize;
        // a small alphabet of repeats makes long chains
        let motif: Vec<u8> = (0..1 + rng.below(40)).map(|_| rng.below(4) as u8).collect();
        let seq: Vec<u8> = (0..len)
            .map(|i| {
                let r = rng.below(100);
                if r < 2 {
                    14 // N
                } else if r < 50 {
                    motif[i % motif.len()]
                } else {
                    rng.below(4) as u8
                }
            })
            .collect();
        let mut m = Vec::new();
        let mut cursor = 0usize;
        while cursor < len && rng.below(3) != 0 {
            let start = cursor + rng.below(400) as usize;
            let end = start + 1 + rng.below(200) as usize;
            if start >= len {
                break;
            }
            m.push(MaskedInterval {
                start,
                end: end.min(len),
            });
            cursor = end;
        }
        offsets.push(offset);
        offset += len as i32 + 1;
        queries.push(seq);
        masks.push(m);
    }
    (queries, offsets, masks, word_length, lut_word_length)
}

/// The fill of this module without the private copies and without the read-once fill (the
/// loop's statements, part by part; checked against the loop below), for the tests of
/// `x_mb_lookup_once.rs`.
pub(crate) fn fill_old(parts: usize, words: &Words<'_>, tables: Tables<'_>) {
    fill_mode(false, false, parts, words, tables, false, false);
}

#[test]
fn partitioned_fill_builds_the_same_tables() {
    let mut rng = Lcg(0x5EED_1234);
    for case in 0..60 {
        let (queries, offsets, masks, word_length, lut_word_length) = random_case(&mut rng);
        let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let ascending_cells = case % 2 == 1;
        let approx = queries.iter().map(Vec::len).sum::<usize>();
        let reference = build_two_stage_lookup(
            &query_refs,
            &offsets,
            word_length,
            lut_word_length,
            &masks,
            None,
            0,
            approx,
            ascending_cells,
        );
        let mb = &reference.mb_lookup;
        // (a table of 3-byte elements, LOSAT_X_MB24, gives its chains back as 4-byte arrays)
        let (cells, next_len) = (
            1usize << (2 * lut_word_length),
            next_len_of(&queries, &offsets),
        );
        let (mb_hashtable, mb_next_pos) = reference.x_chain_arrays(next_len);
        for parts in 1..=8usize {
            let mut hashtable = vec![0u32; cells];
            let mut next_pos = vec![0u32; next_len];
            let mut pv_array = vec![0 as PvArrayType; mb.pv_array.len()];
            let mut helper_array = vec![0u32; (cells >> 11).max(1)];
            let words = Words {
                queries_blastna: &query_refs,
                query_offsets: &offsets,
                query_masks: &masks,
                word_length,
                lut_word_length,
                kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                ascending_cells,
            };
            // The fill of this module (the read-once fill has its own tests in
            // `x_mb_lookup_once.rs`).
            fill_mode(
                false,
                false,
                parts,
                &words,
                Tables {
                    hashtable: &mut hashtable,
                    next_pos: &mut next_pos,
                    pv_array: &mut pv_array,
                    pv_array_bts: mb.pv_array_bts,
                    helper_array: &mut helper_array,
                    x_h24: false,
                },
                false,
                false,
            );
            let mut longest_chain = 2usize;
            for &value in &helper_array {
                longest_chain = longest_chain.max(value as usize);
            }
            assert_eq!(hashtable, mb_hashtable, "case {case} parts {parts}");
            assert_eq!(next_pos, mb_next_pos, "case {case} parts {parts}");
            assert_eq!(pv_array, mb.pv_array, "case {case} parts {parts}");
            assert_eq!(longest_chain, mb.longest_chain, "case {case} parts {parts}");
            // LOSAT_X_MB24: the same fill into 3-byte elements (newest first only)
            if !ascending_cells {
                let (h, n, pv, helper) = fill_h24(
                    false,
                    parts,
                    false,
                    &words,
                    next_len,
                    mb.pv_array.len(),
                    mb.pv_array_bts,
                );
                assert_eq!(h, mb_hashtable, "case {case} parts {parts} h24");
                assert_eq!(n, mb_next_pos, "case {case} parts {parts} h24");
                assert_eq!(pv, mb.pv_array, "case {case} parts {parts} h24");
                assert_eq!(helper, helper_array, "case {case} parts {parts} h24");
            }
        }
    }
}

/// The elements of `next_pos` for these queries (the largest offset plus length, plus 1).
pub(crate) fn next_len_of(queries: &[Vec<u8>], offsets: &[i32]) -> usize {
    queries
        .iter()
        .zip(offsets)
        .map(|(q, &o)| o.max(0) as usize + q.len())
        .max()
        .unwrap_or(0)
        + 1
}

/// The fill of this module (`once`: the read-once fill) into tables of 3-byte elements
/// (LOSAT_X_MB24) of `words.kmer_mask + 1` cells, widened: `(hashtable, next_pos, pv_array,
/// helper_array)`.
pub(crate) fn fill_h24(
    once: bool,
    parts: usize,
    in_pool: bool,
    words: &Words<'_>,
    next_len: usize,
    pv_len: usize,
    pv_array_bts: usize,
) -> (Vec<u32>, Vec<u32>, Vec<PvArrayType>, Vec<u32>) {
    use super::super::x_mb_narrow::{widen, words_for};
    let cells = words.kmer_mask as usize + 1;
    let mut hashtable = vec![0u32; words_for(cells)];
    let mut next_pos = vec![0u32; words_for(next_len)];
    let mut pv_array = vec![0 as PvArrayType; pv_len];
    let mut helper_array = vec![0u32; (cells >> 11).max(1)];
    fill_mode(
        once,
        false,
        parts,
        words,
        Tables {
            hashtable: &mut hashtable,
            next_pos: &mut next_pos,
            pv_array: &mut pv_array,
            pv_array_bts,
            helper_array: &mut helper_array,
            x_h24: true,
        },
        in_pool,
        false,
    );
    (
        widen(&hashtable, cells),
        widen(&next_pos, next_len),
        pv_array,
        helper_array,
    )
}

// The same with real threads (a pool of 1..8 threads; the owner is inside the pool).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn partitioned_fill_on_a_pool_builds_the_same_tables() {
    let mut rng = Lcg(0xFACE_0FF5);
    for case in 0..12 {
        let (queries, offsets, masks, word_length, lut_word_length) = random_case(&mut rng);
        let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let ascending_cells = case % 2 == 0;
        let reference = build_two_stage_lookup(
            &query_refs,
            &offsets,
            word_length,
            lut_word_length,
            &masks,
            None,
            0,
            1_000_000,
            ascending_cells,
        );
        let mb = &reference.mb_lookup;
        let (cells, next_len) = (
            1usize << (2 * lut_word_length),
            next_len_of(&queries, &offsets),
        );
        let (mb_hashtable, mb_next_pos) = reference.x_chain_arrays(next_len);
        for threads in 1..=8usize {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let words = Words {
                queries_blastna: &query_refs,
                query_offsets: &offsets,
                query_masks: &masks,
                word_length,
                lut_word_length,
                kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                ascending_cells,
            };
            // LOSAT_X_MB24: the fill of this module and the read-once fill into 3-byte
            // elements on the pool (newest first only)
            if !ascending_cells {
                for once in [false, true] {
                    let (h, n, pv, _) = pool.install(|| {
                        fill_h24(
                            once,
                            threads,
                            true,
                            &words,
                            next_len,
                            mb.pv_array.len(),
                            mb.pv_array_bts,
                        )
                    });
                    let what = format!("case {case} threads {threads} h24 once {once}");
                    assert_eq!(h, mb_hashtable, "{what}");
                    assert_eq!(n, mb_next_pos, "{what}");
                    assert_eq!(pv, mb.pv_array, "{what}");
                }
            }
            let (hashtable, next_pos, pv_array) = pool.install(|| {
                let mut hashtable = vec![0u32; cells];
                let mut next_pos = vec![0u32; next_len];
                let mut pv_array = vec![0 as PvArrayType; mb.pv_array.len()];
                let mut helper_array = vec![0u32; (cells >> 11).max(1)];
                let words = Words {
                    queries_blastna: &query_refs,
                    query_offsets: &offsets,
                    query_masks: &masks,
                    word_length,
                    lut_word_length,
                    kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                    ascending_cells,
                };
                fill_mode(
                    false,
                    false,
                    threads,
                    &words,
                    Tables {
                        hashtable: &mut hashtable,
                        next_pos: &mut next_pos,
                        pv_array: &mut pv_array,
                        pv_array_bts: mb.pv_array_bts,
                        helper_array: &mut helper_array,
                        x_h24: false,
                    },
                    true,
                    false,
                );
                (hashtable, next_pos, pv_array)
            });
            assert_eq!(hashtable, mb_hashtable, "case {case} threads {threads}");
            assert_eq!(next_pos, mb_next_pos, "case {case} threads {threads}");
            assert_eq!(pv_array, mb.pv_array, "case {case} threads {threads}");
        }
    }
}

// The tables of one case, built by the fill of this module in one part (the loop's
// statements in reading order; the tests above check it against `build_two_stage_lookup`),
// from start values `helper0` / `pv0` of the counters and presence-vector words.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Built {
    hashtable: Vec<u32>,
    next_pos: Vec<u32>,
    pv_array: Vec<PvArrayType>,
    helper_array: Vec<u32>,
    stats: FillStats,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn build_with(
    parts: usize,
    in_pool: bool,
    words: &Words<'_>,
    shape: (usize, usize, usize, usize),
    helper0: &[u32],
    pv0: &[PvArrayType],
) -> Built {
    let (cells, next_len, pv_len, pv_array_bts) = shape;
    let mut b = Built {
        hashtable: vec![0u32; cells],
        next_pos: vec![0u32; next_len],
        pv_array: pv0.to_vec(),
        helper_array: helper0.to_vec(),
        stats: FillStats::default(),
    };
    assert_eq!(b.pv_array.len(), pv_len);
    b.stats = fill_mode(
        false,
        false,
        parts,
        words,
        Tables {
            hashtable: &mut b.hashtable,
            next_pos: &mut b.next_pos,
            pv_array: &mut b.pv_array,
            pv_array_bts,
            helper_array: &mut b.helper_array,
            x_h24: false,
        },
        in_pool,
        false,
    );
    b
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn assert_same(a: &Built, b: &Built, what: &str) {
    assert_eq!(a.hashtable, b.hashtable, "{what}: hashtable");
    assert_eq!(a.next_pos, b.next_pos, "{what}: next_pos");
    assert_eq!(a.pv_array, b.pv_array, "{what}: pv_array");
    assert_eq!(a.helper_array, b.helper_array, "{what}: helper_array");
}

// A query with long unmasked repeats (many words per cell, so counters grow) and random bases.
pub(crate) fn repeat_case(rng: &mut Lcg, len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| {
            if (i / 700) % 3 == 0 {
                [0u8, 0, 1, 0, 2][i % 5]
            } else if rng.below(200) == 0 {
                14
            } else {
                rng.below(4) as u8
            }
        })
        .collect()
}

// The fill on real pools of 2, 3, 5 and 16 threads (and off the pool with the same
// parts): every element of the four arrays and the longest chain equal the one-part fill's
// and `build_two_stage_lookup`'s, for random cases (ambiguity codes, masks, several queries,
// `lut_word_length` 4..12, both `ascending_cells`) and for `lut_word_length` 11 and 12 tables
// (`pv_array_bts` 8 to 10).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn fill_gives_the_loops_tables_on_pools() {
    let mut rng = Lcg(0x0B1D_5EED);
    let mut cases = Vec::new();
    for case in 0..8 {
        let (queries, offsets, masks, word_length, lut_word_length) = random_case(&mut rng);
        cases.push((
            queries,
            offsets,
            masks,
            word_length,
            lut_word_length,
            case % 2 == 1,
        ));
    }
    for (lut, asc) in [(11usize, false), (12, false), (12, true)] {
        let q = repeat_case(&mut rng, 40_000);
        cases.push((vec![q], vec![0], vec![Vec::new()], lut + 5, lut, asc));
    }
    let mut seen_large_tables = 0;
    for (case, (queries, offsets, masks, word_length, lut_word_length, asc)) in
        cases.iter().enumerate()
    {
        let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let reference = build_two_stage_lookup(
            &query_refs,
            offsets,
            *word_length,
            *lut_word_length,
            masks,
            None,
            0,
            1_000_000,
            *asc,
        );
        let mb = &reference.mb_lookup;
        let (cells, next_len) = (
            1usize << (2 * lut_word_length),
            next_len_of(queries, offsets),
        );
        // (a table of 3-byte elements, LOSAT_X_MB24, gives its chains back as 4-byte arrays)
        let (mb_hashtable, mb_next_pos) = reference.x_chain_arrays(next_len);
        let shape = (cells, next_len, mb.pv_array.len(), mb.pv_array_bts);
        if cells >= 1 << 22 && mb.pv_array_bts > 7 {
            seen_large_tables += 1;
        }
        let helper0 = vec![0u32; (cells >> 11).max(1)];
        let pv0 = vec![0 as PvArrayType; mb.pv_array.len()];
        let words = Words {
            queries_blastna: &query_refs,
            query_offsets: offsets,
            query_masks: masks,
            word_length: *word_length,
            lut_word_length: *lut_word_length,
            kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
            ascending_cells: *asc,
        };
        let loop_built = build_with(1, false, &words, shape, &helper0, &pv0);
        assert_eq!(loop_built.hashtable, mb_hashtable, "case {case}");
        assert_eq!(loop_built.next_pos, mb_next_pos, "case {case}");
        assert_eq!(loop_built.pv_array, mb.pv_array, "case {case}");
        let chain = loop_built
            .helper_array
            .iter()
            .fold(2usize, |m, &v| m.max(v as usize));
        assert_eq!(chain, mb.longest_chain, "case {case}");
        for threads in [2usize, 3, 5, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let pooled = pool.install(|| build_with(threads, true, &words, shape, &helper0, &pv0));
            assert_same(&loop_built, &pooled, &format!("case {case} pool {threads}"));
            let serial = build_with(threads, false, &words, shape, &helper0, &pv0);
            assert_same(
                &loop_built,
                &serial,
                &format!("case {case} parts {threads}"),
            );
        }
    }
    assert!(
        seen_large_tables >= 3,
        "the lut 11 / 12 cases ({seen_large_tables})"
    );
}

// Each worker's copies are a whole number of 64-byte lines starting on a line boundary, the
// workers' copies are disjoint and in part order, and zeroed.
#[test]
fn private_copies_have_lines_of_their_own() {
    assert_eq!(std::mem::size_of::<Line>(), 64);
    assert_eq!(std::mem::align_of::<Line>(), 64);
    for (parts, helper_len, pv_len) in [
        (1usize, 1usize, 1usize),
        (2, 2048, 16384),
        (3, 8192, 32768),
        (5, 15, 17),
        (16, 2048, 32768),
        (16, 8192, 16384),
    ] {
        let mut copies = Copies::new(parts, helper_len, pv_len);
        let stride = helper_len.div_ceil(16) + pv_len.div_ceil(16);
        assert!(copies.lines.capacity() >= stride * parts);
        let base = copies.lines.as_ptr() as usize;
        assert_eq!(base % 64, 0);
        let mut seen = 0;
        let mut ends = Vec::new();
        for (p, x) in copies.parts().enumerate() {
            let mut x = x.zeroed();
            let c0 = x.counters.lines.as_ptr() as usize;
            let v0 = x.pv.lines.as_ptr() as usize;
            assert_eq!(c0, base + p * stride * 64, "parts {parts} part {p}");
            assert_eq!(
                v0,
                c0 + helper_len.div_ceil(16) * 64,
                "parts {parts} part {p}"
            );
            assert_eq!(x.counters.lines.len(), helper_len.div_ceil(16));
            assert_eq!(x.pv.lines.len(), pv_len.div_ceil(16));
            assert!((0..helper_len).all(|i| x.counters.get(i) == 0));
            assert!((0..pv_len).all(|i| x.pv.get(i) == 0));
            assert!(!x.has_values());
            *x.counters.slot(helper_len - 1) = 1;
            *x.pv.slot(pv_len - 1) = 1;
            assert!(x.has_values());
            ends.push((c0, v0 + x.pv.lines.len() * 64));
            seen += 1;
        }
        assert_eq!(seen, parts);
        for pair in ends.windows(2) {
            assert!(pair[0].1 <= pair[1].0, "parts {parts}: copies overlap");
        }
    }
}

// The kept buffers, on a list of this test's own (the module's other tests run at
// the same time on `SPARE_LINES`). On Linux with glibc a fill takes the largest kept buffer when
// it is large enough; one too small is freed and the fill gets a new buffer of
// `lines_capacity(n)` lines (32 MiB + 4 KiB or more), not the old one grown (std's aligned
// realloc frees the old block, which moves glibc's mmap threshold); `Drop` gives the buffer back,
// and frees it without waiting when the lock is taken. Elsewhere nothing is kept.
#[test]
fn kept_buffers_are_taken_largest_first_never_grown_and_given_back() {
    static SPARE: SpareLines = std::sync::Mutex::new(Vec::new());
    let kept = || -> Vec<(usize, usize)> {
        SPARE
            .lock()
            .unwrap()
            .iter()
            .map(|l| (l.as_ptr() as usize, l.capacity()))
            .collect()
    };
    let at = |c: &Copies| (c.lines.as_ptr() as usize, c.lines.capacity());
    // 4 parts of 2 + 16 lines.
    let (parts, helper_len, pv_len) = (4usize, 32usize, 256usize);
    let n = parts * (helper_len / 16 + pv_len / 16);
    let mapped = ((32 << 20) + 4096) / 64;
    assert_eq!(lines_capacity(n), if GLIBC { mapped } else { n });
    assert_eq!(lines_capacity(mapped + 7), mapped + 7);

    // The first fill makes a new buffer of `lines_capacity(n)` lines.
    let copies = Copies::with_spare(&SPARE, parts, helper_len, pv_len);
    let first = at(&copies);
    assert!(first.1 >= lines_capacity(n));
    drop(copies);
    if !GLIBC {
        // Nothing is kept, and a buffer in the list is not taken.
        assert!(kept().is_empty(), "kept off glibc");
        SPARE.lock().unwrap().push(Vec::with_capacity(4 * n));
        let before = kept();
        drop(Copies::with_spare(&SPARE, parts, helper_len, pv_len));
        assert_eq!(kept(), before, "taken or kept off glibc");
        SPARE.lock().unwrap().clear();
        return;
    }
    assert_eq!(kept(), [first], "Drop gives the buffer back");
    // The next fill takes it, and gives it back.
    let copies = Copies::with_spare(&SPARE, parts, helper_len, pv_len);
    assert_eq!(at(&copies), first);
    assert!(kept().is_empty());
    drop(copies);
    assert_eq!(kept(), [first]);
    // With a smaller buffer kept too, the fill takes the largest.
    SPARE.lock().unwrap().push(Vec::with_capacity(n));
    let copies = Copies::with_spare(&SPARE, parts, helper_len, pv_len);
    assert_eq!(at(&copies), first, "not the largest kept buffer");
    let small = kept();
    assert_eq!(small.len(), 1);
    assert!(small[0].1 >= n && small[0].1 < first.1);
    // The lock taken at the end of a fill: the buffer is freed, not kept, and nothing waits.
    {
        let _taken = SPARE.lock().unwrap();
        drop(copies);
    }
    assert_eq!(kept(), small, "kept although the lock was taken");
    // A fill larger than the only kept buffer: that one is freed, and the fill gets a new
    // buffer of `lines_capacity` lines, not the old one grown to its size.
    let copies = Copies::with_spare(&SPARE, 2 * parts, helper_len, pv_len);
    assert!(
        copies.lines.capacity() >= lines_capacity(2 * n),
        "grown to {} lines instead of a new buffer of {}",
        copies.lines.capacity(),
        lines_capacity(2 * n)
    );
    assert!(kept().is_empty(), "the buffer too small is kept");
    drop(copies);
    assert_eq!(kept().len(), 1);
    // More than 32 MiB of copies (205 parts of 160 KiB, not touched here): the kept buffer is
    // too small, and the fill gets a new one of its own size.
    let copies = Copies::with_spare(&SPARE, 205, 8192, 32768);
    assert_eq!(copies.stride, 2560);
    assert!(copies.lines.capacity() >= 205 * 2560);
    assert!(kept().is_empty());
    drop(copies);
    assert_eq!(kept().len(), 1);
    SPARE.lock().unwrap().clear();
}

// Many parts: off the pool 64 / 254 / 255 / 256 / 300 parts (the last two clamp to 255), on a
// 64-thread pool 64 and 255 parts (more parts than threads, and at lut 8 more parts than
// blocks); lut 8 / 11 / 12; the fill equals the loop.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn fill_with_many_parts_gives_the_loops_tables() {
    let mut rng = Lcg(0x00C0_FFEE);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(64)
        .build()
        .unwrap();
    for (lut, asc, approx) in [
        (8usize, true, 30_000usize),
        (11, false, 1_000_000),
        (12, true, 400_000),
    ] {
        let q = repeat_case(&mut rng, 25_000);
        let queries = [q];
        let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let offsets = vec![0i32];
        let masks = vec![Vec::new()];
        let reference =
            build_two_stage_lookup(&refs, &offsets, lut + 4, lut, &masks, None, 0, approx, asc);
        let mb = &reference.mb_lookup;
        // (a table laid out as CSR, LOSAT_X_BNCSR, gives its chains back as the fill's arrays)
        let (hashtable, next_pos) = reference.x_chain_arrays(queries[0].len() + 1);
        let shape = (
            hashtable.len(),
            next_pos.len(),
            mb.pv_array.len(),
            mb.pv_array_bts,
        );
        let words = Words {
            queries_blastna: &refs,
            query_offsets: &offsets,
            query_masks: &masks,
            word_length: lut + 4,
            lut_word_length: lut,
            kmer_mask: (1u64 << (2 * lut)) - 1,
            ascending_cells: asc,
        };
        let helper0 = vec![0u32; (hashtable.len() >> 11).max(1)];
        let pv0 = vec![0 as PvArrayType; mb.pv_array.len()];
        let loop_built = build_with(1, false, &words, shape, &helper0, &pv0);
        assert_eq!(loop_built.hashtable, hashtable);
        assert_eq!(loop_built.next_pos, next_pos);
        assert_eq!(loop_built.pv_array, mb.pv_array);
        for parts in [64usize, 254, 255, 256, 300] {
            let b = build_with(parts, false, &words, shape, &helper0, &pv0);
            assert_same(&loop_built, &b, &format!("lut {lut} parts {parts}"));
        }
        for parts in [64usize, 255] {
            let b = pool.install(|| build_with(parts, true, &words, shape, &helper0, &pv0));
            assert_same(&loop_built, &b, &format!("lut {lut} pool 64 parts {parts}"));
        }
    }
}

// Every counter and presence-vector word is flushed by exactly one part, its block's owner,
// for every reachable table shape (lut 4..12; `approx_table_entries` below, between and above
// the cutoffs that halve the presence vector; with and without the database filter) and
// P = 1..=17, 31, 64, 127, 128, 254, 255, 256 and 300 (the last two clamp): each part's
// copies hold `part + 1` everywhere, so after all flushes a counter
// holds its owner + 1 (the flush adds to 0) and a word its owner + 1 (one OR). A double flush,
// a missed one or a wrong owner fails. Also: no presence-vector word spans two blocks, and
// the copies are `P x (ceil(hashsize / 32768) + ceil(pv_size / 16))` lines, at most 160 KiB per
// part (the module documentation's memory sentence).
#[test]
fn every_counter_and_word_is_flushed_once_by_its_owner() {
    let mut shapes = 0;
    for lut in 4usize..=12 {
        let mut worst = 0usize;
        for approx in [1_000usize, 20_000, 1_000_000] {
            for db in [false, true] {
                let table = 1usize << (2 * lut);
                let (pv_len, bts) = super::super::compute_mb_pv_params(table, approx, db, lut);
                let helper_len = table.div_ceil(2048).max(1);
                let per_part = table.div_ceil(32768) + pv_len.div_ceil(16);
                worst = worst.max(per_part * 64);
                for parts_in in (1usize..=17).chain([31, 64, 127, 128, 254, 255, 256, 300]) {
                    let (block_shift, owner) = plan(table, bts, parts_in);
                    let parts = parts_in.clamp(1, 255);
                    assert!(block_shift >= bts && block_shift >= 11);
                    let mut helper = vec![0u32; helper_len];
                    let mut pv = vec![0 as PvArrayType; pv_len];
                    let mut hashtable = vec![0u32; 0];
                    let mut next_pos = vec![0u32; 0];
                    let t = SharedTables {
                        hashtable: Shared::new(&mut hashtable),
                        next_pos: Shared::new(&mut next_pos),
                        pv_array: Shared::new(&mut pv),
                        pv_array_bts: bts,
                        helper_array: Shared::new(&mut helper),
                    };
                    let mut copies = Copies::new(parts, helper_len, pv_len);
                    assert!(copies.lines.capacity() >= parts * per_part, "lut {lut}");
                    assert_eq!(copies.stride, per_part, "lut {lut}");
                    for (p, x) in copies.parts().enumerate() {
                        let x = x.zeroed();
                        for line in x.counters.lines.iter_mut().chain(x.pv.lines.iter_mut()) {
                            line.0 = [p as u32 + 1; 16];
                        }
                        flush_private(p as u8, &owner, block_shift, t, &x);
                    }
                    let what = format!("lut {lut} approx {approx} db {db} P {parts_in}");
                    for (h, &value) in helper.iter().enumerate() {
                        let block = (h << 11) >> block_shift;
                        let o = owner[block.min(owner.len() - 1)] as u32;
                        assert_eq!(value, o + 1, "{what} counter {h}");
                    }
                    for (w, &value) in pv.iter().enumerate() {
                        let block = (w << bts) >> block_shift;
                        let o = owner[block.min(owner.len() - 1)] as u32;
                        assert_eq!(value, o + 1, "{what} word {w}");
                        assert_eq!(block, (((w + 1) << bts) - 1) >> block_shift, "{what}");
                    }
                    shapes += 1;
                }
            }
        }
        if lut >= 11 {
            assert!(worst <= 160 * 1024, "lut {lut}: {worst} bytes per part");
        }
    }
    assert!(shapes > 1000, "{shapes}");
}

// LOSAT_X_MB24 at its edges. Words of the cells on both sides of block boundaries (`plan`'s
// blocks of 2^11 to 2^20 cells, whose owners differ) and of the first and last cells, 40 of
// each (long chains), through the
// serial build and both pool fills on pools of 2, 5 and 16 threads; and one query placed so
// that `next_pos` has 2^24 elements (its last offset 2^24 - 1: 3-byte elements) or 2^24 + 1
// (4-byte, as before). The chains, presence vector and counters equal the loop's.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn three_byte_tables_at_block_edges_and_at_the_last_offset() {
    let lut = 12usize;
    let cells = 1usize << (2 * lut);
    let kmer = |code: usize| -> Vec<u8> {
        (0..lut)
            .rev()
            .map(|k| ((code >> (2 * k)) & 3) as u8)
            .collect()
    };
    // and the first and last cells (the last one's 4-byte read takes the slack byte)
    let mut edges = vec![0, cells - 1];
    for s in 11..=20 {
        for k in [1usize, 3, 7] {
            if k << s < cells {
                edges.extend([(k << s) - 1, k << s]);
            }
        }
    }
    let mut seq = Vec::new();
    for _ in 0..40 {
        for &cell in &edges {
            seq.extend(kmer(cell));
            // an N, so that the only words are those of the edge cells
            seq.push(14);
        }
    }
    let mut rng = Lcg(0x0024_ED6E);
    let tail: Vec<u8> = (0..200).map(|_| rng.below(4) as u8).collect();
    for (queries, offsets) in [
        (vec![seq.clone()], vec![0i32]),
        (vec![tail.clone()], vec![(1i32 << 24) - 201]),
        (vec![tail.clone()], vec![(1i32 << 24) - 200]),
    ] {
        let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let masks = vec![Vec::new()];
        let next_len = next_len_of(&queries, &offsets);
        let what = format!("next_pos of {next_len}");
        let reference =
            build_two_stage_lookup(&refs, &offsets, 16, lut, &masks, None, 0, 1_000_000, false);
        let mb = &reference.mb_lookup;
        assert_eq!(
            mb.x_narrow.is_some(),
            super::super::x_mb_narrow::on() && next_len <= 1 << 24,
            "{what}"
        );
        let words = Words {
            queries_blastna: &refs,
            query_offsets: &offsets,
            query_masks: &masks,
            word_length: 16,
            lut_word_length: lut,
            kmer_mask: (cells - 1) as u64,
            ascending_cells: false,
        };
        let shape = (cells, next_len, mb.pv_array.len(), mb.pv_array_bts);
        let helper0 = vec![0u32; cells >> 11];
        let pv0 = vec![0 as PvArrayType; mb.pv_array.len()];
        let loop_built = build_with(1, false, &words, shape, &helper0, &pv0);
        let (hashtable, next_pos) = reference.x_chain_arrays(next_len);
        assert!(hashtable == loop_built.hashtable, "{what}: serial heads");
        assert!(next_pos == loop_built.next_pos, "{what}: serial chains");
        assert_eq!(mb.pv_array, loop_built.pv_array, "{what}: serial");
        if next_len > 1 << 24 {
            continue;
        }
        for threads in [2usize, 5, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for once in [false, true] {
                let (h, n, pv, helper) = pool
                    .install(|| fill_h24(once, threads, true, &words, next_len, shape.2, shape.3));
                let what = format!("{what} pool {threads} once {once}");
                assert!(h == loop_built.hashtable, "{what}: heads");
                assert!(n == loop_built.next_pos, "{what}: chains");
                assert_eq!(pv, loop_built.pv_array, "{what}");
                assert_eq!(helper, loop_built.helper_array, "{what}");
            }
        }
    }
}
