use super::super::x_disc_lookup_par::PAGE_BYTES;
use super::super::x_mb_lookup_par::tests::{random_case, repeat_case, Lcg};
use super::super::{build_two_stage_lookup, MaskedInterval};
use super::*;

/// The tables of one fill.
#[derive(Debug, PartialEq, Eq)]
struct Built {
    hashtable: Vec<u32>,
    next_pos: Vec<u32>,
    pv_array: Vec<PvArrayType>,
    helper_array: Vec<u32>,
}

/// One case: the queries, their offsets and masks, `word_length`, `lut_word_length` and
/// `ascending_cells`.
struct Case {
    queries: Vec<Vec<u8>>,
    offsets: Vec<i32>,
    masks: Vec<Vec<MaskedInterval>>,
    word_length: usize,
    lut: usize,
    ascending: bool,
}

impl Case {
    fn random(rng: &mut Lcg, ascending: bool) -> Self {
        let (queries, offsets, masks, word_length, lut) = random_case(rng);
        Self {
            queries,
            offsets,
            masks,
            word_length,
            lut,
            ascending,
        }
    }

    fn one(query: Vec<u8>, word_length: usize, lut: usize, ascending: bool) -> Self {
        Self {
            queries: vec![query],
            offsets: vec![0],
            masks: vec![Vec::new()],
            word_length,
            lut,
            ascending,
        }
    }

    fn refs(&self) -> Vec<&[u8]> {
        self.queries.iter().map(Vec::as_slice).collect()
    }

    fn words<'a>(&'a self, refs: &'a [&'a [u8]]) -> Words<'a> {
        Words {
            queries_blastna: refs,
            query_offsets: &self.offsets,
            query_masks: &self.masks,
            word_length: self.word_length,
            lut_word_length: self.lut,
            kmer_mask: (1u64 << (2 * self.lut)) - 1,
            ascending_cells: self.ascending,
        }
    }

    /// `build_two_stage_lookup`'s table (the loop), with its counters (from 0) recomputed by
    /// the loop's statements in one part of `x_mb_lookup_par`'s fill (the tables there are
    /// checked against the loop) and its `pv_array_bts`.
    fn reference(&self, approx: usize) -> (Built, usize, usize) {
        let refs = self.refs();
        let r = build_two_stage_lookup(
            &refs,
            &self.offsets,
            self.word_length,
            self.lut,
            &self.masks,
            None,
            0,
            approx,
            self.ascending,
        );
        let mb = &r.mb_lookup;
        // (a table laid out as CSR, LOSAT_X_BNCSR, gives its chains back as the fill's arrays)
        let next_pos_len = self
            .queries
            .iter()
            .zip(&self.offsets)
            .map(|(seq, &off)| off.max(0) as usize + seq.len())
            .max()
            .unwrap_or(0)
            + 1;
        let (hashtable, next_pos) = r.x_chain_arrays(next_pos_len);
        let helper = self.helper_of(hashtable.len(), mb.pv_array.len(), mb.pv_array_bts);
        (
            Built {
                hashtable,
                next_pos,
                pv_array: mb.pv_array.clone(),
                helper_array: helper,
            },
            mb.pv_array_bts,
            mb.longest_chain,
        )
    }

    /// The loop's counters: the fill of `x_mb_lookup_par` in one part, which applies
    /// the loop's statements in reading order.
    fn helper_of(&self, cells: usize, pv_len: usize, bts: usize) -> Vec<u32> {
        let refs = self.refs();
        let words = self.words(&refs);
        let mut hashtable = vec![0u32; cells];
        let mut next_pos = vec![0u32; self.next_len()];
        let mut pv_array = vec![0 as PvArrayType; pv_len];
        let mut helper_array = vec![0u32; (cells >> 11).max(1)];
        super::super::x_mb_lookup_par::tests::fill_old(
            1,
            &words,
            Tables {
                hashtable: &mut hashtable,
                next_pos: &mut next_pos,
                pv_array: &mut pv_array,
                pv_array_bts: bts,
                helper_array: &mut helper_array,
                x_h24: false,
            },
        );
        helper_array
    }

    fn next_len(&self) -> usize {
        self.queries
            .iter()
            .zip(&self.offsets)
            .map(|(q, &o)| o.max(0) as usize + q.len())
            .max()
            .unwrap_or(0)
            + 1
    }

    /// The read-once fill, from start values `helper0` / `pv0`, with the prefault round
    /// (LOSAT_X_MBLUTPREFAULT), which `check_prefault` checks.
    fn once(
        &self,
        shape: (usize, usize),
        parts: usize,
        in_pool: bool,
        chunk: usize,
        helper0: &[u32],
        pv0: &[PvArrayType],
    ) -> (Built, OnceStats) {
        let (b, stats) = self.once_mode(shape, parts, in_pool, chunk, helper0, pv0, true);
        check_prefault(&b.hashtable, &stats, self.first_chunk(chunk), true);
        (b, stats)
    }

    /// `once` with LOSAT_X_MBLUTPREFAULT given (`prefault`), unchecked.
    #[allow(clippy::too_many_arguments)]
    fn once_mode(
        &self,
        shape: (usize, usize),
        parts: usize,
        in_pool: bool,
        chunk: usize,
        helper0: &[u32],
        pv0: &[PvArrayType],
        prefault: bool,
    ) -> (Built, OnceStats) {
        let (cells, bts) = shape;
        let refs = self.refs();
        let words = self.words(&refs);
        let mut b = Built {
            hashtable: vec![0u32; cells],
            next_pos: vec![0u32; self.next_len()],
            pv_array: pv0.to_vec(),
            helper_array: helper0.to_vec(),
        };
        let stats = fill_chunked(
            parts,
            &words,
            Tables {
                hashtable: &mut b.hashtable,
                next_pos: &mut b.next_pos,
                pv_array: &mut b.pv_array,
                pv_array_bts: bts,
                helper_array: &mut b.helper_array,
                x_h24: false,
            },
            in_pool,
            false,
            chunk,
            prefault,
        );
        (b, stats)
    }

    /// LOSAT_X_MB24: the read-once fill from zero start values into 3-byte elements of `cells`
    /// cells (newest first), with LOSAT_X_MBLUTPREFAULT given (`prefault`) and checked on the
    /// words of the heads (`check_prefault`); the tables widened, and the pages touched.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn once_h24(
        &self,
        reference: &Built,
        bts: usize,
        parts: usize,
        in_pool: bool,
        chunk: usize,
        prefault: bool,
    ) -> (Built, OnceStats, usize) {
        use super::super::x_mb_narrow::{widen, words_for};
        assert!(!self.ascending);
        let cells = reference.hashtable.len();
        let refs = self.refs();
        let words = self.words(&refs);
        let mut hashtable = vec![0u32; words_for(cells)];
        let mut next_pos = vec![0u32; words_for(self.next_len())];
        let mut b = Built {
            hashtable: Vec::new(),
            next_pos: Vec::new(),
            pv_array: vec![0 as PvArrayType; reference.pv_array.len()],
            helper_array: vec![0u32; reference.helper_array.len()],
        };
        let stats = fill_chunked(
            parts,
            &words,
            Tables {
                hashtable: &mut hashtable,
                next_pos: &mut next_pos,
                pv_array: &mut b.pv_array,
                pv_array_bts: bts,
                helper_array: &mut b.helper_array,
                x_h24: true,
            },
            in_pool,
            false,
            chunk,
            prefault,
        );
        let pages = check_prefault(&hashtable, &stats, self.first_chunk(chunk), prefault);
        b.hashtable = widen(&hashtable, cells);
        b.next_pos = widen(&next_pos, self.next_len());
        (b, stats, pages)
    }

    /// The read-once fill from zero start values.
    fn once0(
        &self,
        reference: &Built,
        bts: usize,
        parts: usize,
        in_pool: bool,
        chunk: usize,
    ) -> (Built, OnceStats) {
        self.once(
            (reference.hashtable.len(), bts),
            parts,
            in_pool,
            chunk,
            &vec![0u32; reference.helper_array.len()],
            &vec![0 as PvArrayType; reference.pv_array.len()],
        )
    }

    /// The number of stream positions (the loop's reads).
    fn positions(&self) -> usize {
        let refs = self.refs();
        segments(&self.words(&refs)).1
    }

    /// The positions of the first chunk of the processing order for chunks of `chunk` (the
    /// last stream chunk for `ascending_cells`).
    fn first_chunk(&self, chunk: usize) -> usize {
        let (total, chunk) = (self.positions(), chunk.max(1));
        if self.ascending && total > 0 {
            total - (total.div_ceil(chunk) - 1) * chunk
        } else {
            chunk.min(total)
        }
    }
}

/// The prefault round of a fill of `v` (`hashtable` as filled in place) with
/// LOSAT_X_MBLUTPREFAULT `prefault` and a first chunk of `first` positions, from its log of
/// writes (`OnceStats::touch_log`): each wrote 0 over 0; the round ran once exactly when the
/// switch is on and the first chunk has at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per
/// 1,024 cells of `hashtable`, never otherwise; and the pages it wrote to are exactly the pages
/// the first chunk wrote (`OnceStats::first_pages`), each once, when it ran. So it maps no page
/// the fill does not write, and misses none of the first chunk's. Returns the number of pages
/// touched.
fn check_prefault(v: &[u32], stats: &OnceStats, first: usize, prefault: bool) -> usize {
    let per_page = PREFAULT_MIN_WORDS_PER_PAGE * (v.len() / 1024).max(1);
    let round = prefault && first > 0 && first >= per_page;
    assert_eq!(stats.prefault_rounds, usize::from(round), "{stats:?}");
    assert_eq!(stats.touches, stats.touch_log.len());
    let mut touched = std::collections::BTreeSet::new();
    for t in &stats.touch_log {
        assert!(
            t.array == HASHTABLE && t.before == 0 && t.after == 0 && t.index < v.len(),
            "the prefault wrote {} over {} at element {} of array {}",
            t.after,
            t.before,
            t.index,
            t.array
        );
        let page = (v.as_ptr() as usize + 4 * t.index) / PAGE_BYTES;
        assert!(touched.insert(page), "page {page} touched twice");
    }
    if round {
        assert_eq!(
            touched, stats.first_pages,
            "touched / the first chunk's pages"
        );
    } else {
        assert!(touched.is_empty());
    }
    touched.len()
}

fn longest(helper: &[u32]) -> usize {
    helper.iter().fold(2usize, |m, &v| m.max(v as usize))
}

fn assert_built(reference: &Built, chain: usize, b: &Built, what: &str) {
    assert_eq!(b.hashtable, reference.hashtable, "{what}: hashtable");
    assert_eq!(b.next_pos, reference.next_pos, "{what}: next_pos");
    assert_eq!(b.pv_array, reference.pv_array, "{what}: pv_array");
    assert_eq!(
        b.helper_array, reference.helper_array,
        "{what}: helper_array"
    );
    assert_eq!(longest(&b.helper_array), chain, "{what}: longest_chain");
}

// Random cases (ambiguity codes, masks, repeats, `lut_word_length` 4 / 8 / 9 / 11 / 12,
// `word_length` above it, 1-3 queries, both `ascending_cells`), off the pool with 1..8 parts
// and chunks of 1 (small cases), 7, 1000 positions and the default: every element of the four
// arrays and the longest chain equal the loop's. The sub-ranges are in reading order, the
// chunks are `ceil(positions / chunk)`.
#[test]
fn once_fill_builds_the_loops_tables() {
    let mut rng = Lcg(0x0ECE_5EED);
    let mut multi_chunk = 0;
    for case in 0..40 {
        let c = Case::random(&mut rng, case % 2 == 1);
        let (reference, bts, chain) = c.reference(1_000_000);
        assert_eq!(longest(&reference.helper_array), chain, "case {case}");
        let positions = c.positions();
        for parts in 1..=8usize {
            for chunk in [1usize, 7, 1000, CHUNK_POSITIONS] {
                if chunk == 1 && (positions > 3000 || parts > 3) {
                    continue;
                }
                let (b, stats) = c.once0(&reference, bts, parts, false, chunk);
                let what = format!("case {case} parts {parts} chunk {chunk}");
                assert_built(&reference, chain, &b, &what);
                assert_eq!(stats.chunks, positions.div_ceil(chunk), "{what}");
                if positions > 0 {
                    assert_eq!(
                        stats.arena_lines,
                        Shape::new(chunk.min(positions), parts.clamp(1, 255)).lines(),
                        "{what}"
                    );
                }
                assert_eq!(stats.subranges, stats.subranges_in_order, "{what}");
                multi_chunk += usize::from(stats.chunks >= 2 && stats.subranges > 0);
            }
        }
    }
    assert!(multi_chunk > 100, "{multi_chunk}");
}

// Edge cases: ranges shorter than `word_length` (skipped as the loop does), queries shorter
// than `lut_word_length`, all-N queries, no words at all, a word ending exactly at a span or
// chunk boundary after an ambiguity, and queries placed at offsets with gaps.
#[test]
fn short_ranges_short_queries_and_no_words() {
    let n = |len: usize| vec![14u8; len];
    let acgt = |len: usize, seed: u64| {
        let mut rng = Lcg(seed);
        (0..len).map(|_| rng.below(4) as u8).collect::<Vec<u8>>()
    };
    let mut cases = vec![
        // all N
        Case::one(n(500), 12, 8, false),
        // shorter than lut
        Case::one(acgt(6, 1), 8, 8, true),
        // empty
        Case::one(Vec::new(), 12, 8, false),
    ];
    // ranges shorter than word_length between masks; one N every 9 bases
    let q = acgt(2000, 2);
    let mut masks = Vec::new();
    for start in (0..2000).step_by(40) {
        masks.push(MaskedInterval {
            start: start + 15 + (start % 7),
            end: start + 40,
        });
    }
    cases.push(Case {
        queries: vec![q.clone(), acgt(5, 3), n(30), q.clone()],
        offsets: vec![0, 2001, 2010, 2100],
        masks: vec![masks.clone(), Vec::new(), Vec::new(), Vec::new()],
        word_length: 16,
        lut: 11,
        ascending: false,
    });
    let mut q9 = acgt(3000, 4);
    for i in (0..3000).step_by(9) {
        q9[i] = 14;
    }
    cases.push(Case::one(q9, 8, 8, true));
    cases.push(Case {
        queries: vec![q.clone(), q],
        offsets: vec![0, 2001],
        masks: vec![masks.clone(), masks],
        word_length: 12,
        lut: 12,
        ascending: true,
    });
    // Queries of exactly `lut_word_length` letters (one word each, with `word_length` equal to
    // it) and one letter shorter (no word), at `lut_word_length` 8 and 12, both orders; with
    // other queries around them so the word's offset is not the first.
    let first_exact = cases.len();
    for (lut, asc, seed) in [
        (8usize, false, 5u64),
        (12, true, 6),
        (12, false, 7),
        (8, true, 8),
    ] {
        cases.push(Case {
            queries: vec![
                acgt(lut, seed),
                acgt(lut - 1, seed + 10),
                acgt(lut, seed + 20),
            ],
            offsets: vec![0, lut as i32 + 1, 2 * lut as i32 + 1],
            masks: vec![Vec::new(), Vec::new(), Vec::new()],
            word_length: lut,
            lut,
            ascending: asc,
        });
    }
    for (i, c) in cases.iter().enumerate() {
        let (reference, bts, chain) = c.reference(30_000);
        if i >= first_exact {
            // The two queries of `lut_word_length` letters have one word each (offsets 1 and
            // `2 lut + 2`), the one between them none.
            let mut heads: Vec<u32> = reference
                .hashtable
                .iter()
                .chain(&reference.next_pos)
                .copied()
                .filter(|&v| v != 0)
                .collect();
            heads.sort_unstable();
            let lut = c.lut as u32;
            assert_eq!(heads.len(), 2, "edge {i}: {heads:?}");
            assert!(
                heads.contains(&1) && heads.contains(&(2 * lut + 2)),
                "edge {i}: {heads:?}"
            );
        }
        for parts in [1usize, 2, 3, 5, 8] {
            for chunk in [1usize, 2, 7, 64, CHUNK_POSITIONS] {
                let (b, _) = c.once0(&reference, bts, parts, false, chunk);
                assert_built(
                    &reference,
                    chain,
                    &b,
                    &format!("edge {i} parts {parts} chunk {chunk}"),
                );
            }
        }
    }
}

// On real pools of 2, 3, 5, 16 and 64 threads (the fill inside the pool, as in a search):
// random cases and long unmasked repeats at `lut_word_length` 11 and 12 (many words per cell,
// counters that grow, presence-vector words shared between blocks), default chunks and chunks
// of 1000 positions; at least two owners' copies held values.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn once_fill_on_pools_builds_the_loops_tables() {
    let mut rng = Lcg(0x00C0_FFEE_0001);
    let mut cases = Vec::new();
    for case in 0..6 {
        cases.push(Case::random(&mut rng, case % 2 == 0));
    }
    for (lut, asc) in [(11usize, false), (12, true), (12, false)] {
        cases.push(Case::one(repeat_case(&mut rng, 40_000), lut + 5, lut, asc));
    }
    for (i, c) in cases.iter().enumerate() {
        let (reference, bts, chain) = c.reference(1_000_000);
        for threads in [2usize, 3, 5, 16, 64] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for chunk in [1000usize, CHUNK_POSITIONS] {
                let (b, stats) = pool.install(|| c.once0(&reference, bts, threads, true, chunk));
                let what = format!("case {i} pool {threads} chunk {chunk}");
                assert_built(&reference, chain, &b, &what);
                assert_eq!(stats.subranges, stats.subranges_in_order, "{what}");
                if i >= 6 && threads >= 3 {
                    assert!(stats.owners_with_values >= 2, "{what}: {stats:?}");
                }
            }
        }
    }
}

// Many parts: 64 / 254 / 255 / 256 / 300 off the pool (the last two clamp to 255) and 255 on
// a 64-thread pool, at `lut_word_length` 8 (32 blocks: most owners have none), 11 and 12.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn once_fill_with_many_parts_builds_the_loops_tables() {
    let mut rng = Lcg(0x0255_0255);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(64)
        .build()
        .unwrap();
    for (lut, asc, approx) in [
        (8usize, true, 30_000usize),
        (11, false, 1_000_000),
        (12, true, 400_000),
    ] {
        let c = Case::one(repeat_case(&mut rng, 25_000), lut + 4, lut, asc);
        let (reference, bts, chain) = c.reference(approx);
        for parts in [64usize, 254, 255, 256, 300] {
            let (b, _) = c.once0(&reference, bts, parts, false, CHUNK_POSITIONS);
            assert_built(&reference, chain, &b, &format!("lut {lut} parts {parts}"));
        }
        for chunk in [999usize, CHUNK_POSITIONS] {
            let (b, _) = pool.install(|| c.once0(&reference, bts, 255, true, chunk));
            assert_built(
                &reference,
                chain,
                &b,
                &format!("lut {lut} pool 64 parts 255 chunk {chunk}"),
            );
        }
    }
}

// From start values other than 0 (production starts at 0): counters at u32::MAX - k and
// others, presence-vector words with bits set; the flush's `saturating_add` and OR give the
// loop's saturating increments and kept bits (the loop's values: `x_mb_lookup_par`'s
// fill in one part from the same start values).
#[test]
fn once_flush_matches_the_loop_from_any_start() {
    let mut rng = Lcg(0x5A7_0F1D_0002);
    for (lut, asc) in [(8usize, false), (11, true), (12, false)] {
        let c = Case::one(repeat_case(&mut rng, 30_000), lut + 4, lut, asc);
        let (reference, bts, _) = c.reference(30_000);
        let helper_len = reference.helper_array.len();
        let helper0: Vec<u32> = (0..helper_len)
            .map(|h| match h % 4 {
                0 => u32::MAX,
                1 => u32::MAX - (h as u32 % 7),
                2 => 0,
                _ => (h as u32).wrapping_mul(2_654_435_761),
            })
            .collect();
        let pv0: Vec<PvArrayType> = (0..reference.pv_array.len())
            .map(|w| {
                if w % 3 == 0 {
                    (w as u32).wrapping_mul(40_503)
                } else {
                    0
                }
            })
            .collect();
        let saturated = (0..helper_len)
            .filter(|&h| {
                helper0[h] != u32::MAX
                    && helper0[h].checked_add(reference.helper_array[h]).is_none()
            })
            .count();
        if lut >= 11 {
            assert!(saturated > 0, "lut {lut}: no counter saturates");
        }
        let expect_helper: Vec<u32> = (0..helper_len)
            .map(|h| helper0[h].saturating_add(reference.helper_array[h]))
            .collect();
        let expect_pv: Vec<PvArrayType> = (0..reference.pv_array.len())
            .map(|w| pv0[w] | reference.pv_array[w])
            .collect();
        for parts in [2usize, 3, 5, 16] {
            let (b, _) = c.once(
                (reference.hashtable.len(), bts),
                parts,
                false,
                CHUNK_POSITIONS,
                &helper0,
                &pv0,
            );
            let what = format!("lut {lut} parts {parts}");
            assert_eq!(b.helper_array, expect_helper, "{what}");
            assert_eq!(b.pv_array, expect_pv, "{what}");
            assert_eq!(b.hashtable, reference.hashtable, "{what}");
            assert_eq!(b.next_pos, reference.next_pos, "{what}");
        }
    }
}

// `next_pos` chain order: for megablast a cell's chain lists its offsets newest first, for
// `ascending_cells` in reading order; and the longest chain of a long repeat (counters in the
// thousands) is the loop's.
#[test]
fn chain_order_and_longest_chain() {
    // One 8-mer repeated: "ACGTACGT" period 4, so every 8-mer start at 4k is the same cell.
    let q: Vec<u8> = (0..4000).map(|i| [0u8, 1, 2, 3][i % 4]).collect();
    for asc in [false, true] {
        let c = Case::one(q.clone(), 8, 8, asc);
        let (reference, bts, chain) = c.reference(30_000);
        assert!(chain > 900, "{chain}");
        let (b, _) = c.once0(&reference, bts, 4, false, 100);
        assert_built(&reference, chain, &b, &format!("asc {asc}"));
        // The chain of the cell of the first word (offset 1): walk it.
        let cell = b.hashtable.iter().position(|&v| v != 0).unwrap();
        let mut walk = Vec::new();
        let mut at = b.hashtable[cell];
        while at != 0 {
            walk.push(at);
            at = b.next_pos[at as usize];
        }
        assert!(walk.len() > 900);
        if asc {
            assert!(walk.windows(2).all(|w| w[0] < w[1]), "ascending");
        } else {
            assert!(walk.windows(2).all(|w| w[0] > w[1]), "newest first");
        }
    }
}

// LOSAT_X_MBLUTPREFAULT on pools of 2, 5 and 16 threads and off the pool: long repeats and
// random letters at `lut_word_length` 8, 11 and 12, both `ascending_cells`, first chunks on
// both sides of `PREFAULT_MIN_WORDS_PER_PAGE` positions per page: the tables are the loop's
// with the round and without it, the round touches exactly the first chunk's pages
// (`check_prefault`) and nothing without the switch; some rounds touched thousands of pages,
// and some fills ran without a round because their first chunk was too short.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn prefault_touches_the_first_chunks_pages_and_changes_no_table() {
    let mut rng = Lcg(0x9F_A017_0003);
    let random = |rng: &mut Lcg, len: usize| (0..len).map(|_| rng.below(4) as u8).collect();
    let mut cases = Vec::new();
    for (lut, asc) in [(8usize, false), (11, true), (12, false), (12, true)] {
        cases.push((Case::one(random(&mut rng, 150_000), lut + 4, lut, asc), lut));
        cases.push((
            Case::one(repeat_case(&mut rng, 40_000), lut + 4, lut, asc),
            lut,
        ));
    }
    let (mut most, mut skipped) = (0usize, 0usize);
    for (i, (c, lut)) in cases.iter().enumerate() {
        let (reference, bts, chain) = c.reference(if *lut == 8 { 30_000 } else { 1_000_000 });
        let shape = (reference.hashtable.len(), bts);
        let zeros = (
            vec![0u32; reference.helper_array.len()],
            vec![0 as PvArrayType; reference.pv_array.len()],
        );
        for threads in [1usize, 2, 5, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for chunk in [4_096usize, 70_000, CHUNK_POSITIONS] {
                let what = format!("case {i} lut {lut} pool {threads} chunk {chunk}");
                for prefault in [true, false] {
                    let (b, stats) = pool.install(|| {
                        c.once_mode(
                            shape,
                            threads,
                            threads > 1,
                            chunk,
                            &zeros.0,
                            &zeros.1,
                            prefault,
                        )
                    });
                    assert_built(&reference, chain, &b, &format!("{what} {prefault}"));
                    let pages =
                        check_prefault(&b.hashtable, &stats, c.first_chunk(chunk), prefault);
                    if prefault {
                        most = most.max(pages);
                        skipped += usize::from(stats.prefault_rounds == 0);
                    }
                }
            }
        }
    }
    assert!(most > 4_000, "{most}");
    assert!(skipped > 0);
}

// LOSAT_X_MB24: the read-once fill into 3-byte elements, newest first: random cases (`lut_word_length`
// 4 to 12) off the pool with 1..8 parts and chunks of 7, 1000 positions and the default, then
// long repeats and random letters at 11 and 12 on pools of 2, 5 and 16 threads (with
// LOSAT_X_MBLUTPREFAULT and without; first chunks on both sides of the gate): every element of
// the four arrays and the longest chain equal the loop's, and the round touches exactly the
// pages of the words that the first chunk wrote (`check_prefault` on the words of the heads).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn once_fill_of_3_byte_cells_builds_the_loops_tables() {
    let mut rng = Lcg(0x0003_B17E);
    for case in 0..24 {
        let c = Case::random(&mut rng, false);
        let (reference, bts, chain) = c.reference(1_000_000);
        for parts in 1..=8usize {
            for chunk in [7usize, 1000, CHUNK_POSITIONS] {
                let (b, stats, _) = c.once_h24(&reference, bts, parts, false, chunk, true);
                let what = format!("case {case} parts {parts} chunk {chunk}");
                assert_built(&reference, chain, &b, &what);
                assert_eq!(stats.subranges, stats.subranges_in_order, "{what}");
            }
        }
    }
    let random = |rng: &mut Lcg, len: usize| (0..len).map(|_| rng.below(4) as u8).collect();
    let mut cases = Vec::new();
    for lut in [11usize, 12] {
        cases.push(Case::one(random(&mut rng, 150_000), lut + 4, lut, false));
        cases.push(Case::one(
            repeat_case(&mut rng, 40_000),
            lut + 4,
            lut,
            false,
        ));
    }
    let (mut most, mut skipped) = (0usize, 0usize);
    for (i, c) in cases.iter().enumerate() {
        let (reference, bts, chain) = c.reference(1_000_000);
        for threads in [2usize, 5, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for chunk in [4_096usize, 70_000, CHUNK_POSITIONS] {
                for prefault in [true, false] {
                    let what = format!("case {i} pool {threads} chunk {chunk} {prefault}");
                    let (b, stats, pages) = pool
                        .install(|| c.once_h24(&reference, bts, threads, true, chunk, prefault));
                    assert_built(&reference, chain, &b, &what);
                    if prefault {
                        most = most.max(pages);
                        skipped += usize::from(stats.prefault_rounds == 0);
                    }
                    if threads >= 5 && c.positions() >= 40_000 {
                        assert!(stats.owners_with_values >= 2, "{what}: {stats:?}");
                    }
                }
            }
        }
    }
    assert!(most > 4_000, "{most}");
    assert!(skipped > 0);
}

// LOSAT_X_MBLUTPREFAULT with `hashtable` at element offsets 0, 4 (glibc's mapped blocks) and
// 1,023 (a first page of one element) in its first page, at `lut_word_length` 11 and 12, on a
// pool of 5 threads, with a first chunk just at the gate (8 positions per 1,024 cells) and one
// below it: the tables are the loop's, the round touches exactly the first chunk's pages, and it
// runs on both sides of the gate as `prefault_gate` says, whatever the offset.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn prefault_pages_do_not_depend_on_the_tables_offset() {
    fn aligned(len: usize, offset: usize) -> (Vec<u32>, usize) {
        let v = vec![0u32; len + offset + 1024];
        let skip = (PAGE_BYTES - v.as_ptr() as usize % PAGE_BYTES) % PAGE_BYTES / 4;
        (v, skip + offset)
    }
    let mut rng = Lcg(0x9F_A017_0004);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(5)
        .build()
        .unwrap();
    let mut rounds = 0;
    for lut in [11usize, 12] {
        let q: Vec<u8> = (0..200_000).map(|_| rng.below(4) as u8).collect();
        let c = Case::one(q, lut + 4, lut, false);
        let (reference, bts, chain) = c.reference(1_000_000);
        let cells = reference.hashtable.len();
        let gate = PREFAULT_MIN_WORDS_PER_PAGE * (cells >> 10);
        for chunk in [gate, gate - 1] {
            for offset in [0usize, 4, 1023] {
                let refs = c.refs();
                let words = c.words(&refs);
                let (mut cells_v, c0) = aligned(cells, offset);
                let mut b = Built {
                    hashtable: Vec::new(),
                    next_pos: vec![0u32; c.next_len()],
                    pv_array: vec![0 as PvArrayType; reference.pv_array.len()],
                    helper_array: vec![0u32; reference.helper_array.len()],
                };
                let stats = pool.install(|| {
                    fill_chunked(
                        5,
                        &words,
                        Tables {
                            hashtable: &mut cells_v[c0..c0 + cells],
                            next_pos: &mut b.next_pos,
                            pv_array: &mut b.pv_array,
                            pv_array_bts: bts,
                            helper_array: &mut b.helper_array,
                            x_h24: false,
                        },
                        true,
                        false,
                        chunk,
                        true,
                    )
                });
                let what = format!("lut {lut} chunk {chunk} offset {offset}");
                let table = &cells_v[c0..c0 + cells];
                check_prefault(table, &stats, c.first_chunk(chunk), true);
                assert_eq!(stats.prefault_rounds, usize::from(chunk == gate), "{what}");
                rounds += stats.prefault_rounds;
                b.hashtable = table.to_vec();
                assert_built(&reference, chain, &b, &what);
            }
        }
    }
    assert_eq!(rounds, 6);
}

// The arena: each reader's region and row on lines of their own, the regions within
// `2 x (chunk + 8 P)` entries and the rows `P x ceil((P + 1) / 16)` lines; 8.0 MiB at 16 and 64
// parts with the default chunks.
#[test]
fn arena_regions_and_rows_have_lines_of_their_own() {
    for parts in (1usize..=17).chain([31, 64, 128, 255]) {
        for chunk in [1usize, 7, 1000, CHUNK_POSITIONS] {
            let shape = Shape::new(chunk, parts);
            assert!(
                shape.half >= chunk.div_ceil(parts),
                "P {parts} chunk {chunk}"
            );
            assert_eq!(shape.half % LINE_ENTRIES, 0);
            assert_eq!(shape.region_lines * LINE_ENTRIES, 2 * shape.half);
            assert!(shape.row_lines * LINE_WORDS > parts);
            assert!(
                parts * 2 * shape.half <= 2 * (chunk + 8 * parts),
                "P {parts} chunk {chunk}"
            );
            let mut lines: Vec<Line> = Vec::with_capacity(shape.lines());
            let arena = SharedArena {
                base: lines.spare_capacity_mut().as_mut_ptr() as *mut Line,
                shape,
            };
            let base = arena.base as usize;
            assert_eq!(base % 64, 0);
            let end = base + shape.lines() * 64;
            let mut spans = Vec::new();
            for a in 0..parts {
                let s = arena.scratch(a) as usize;
                let t = arena.sorted(a) as usize;
                let r = arena.row(a) as usize;
                assert_eq!(s % 64, 0);
                assert_eq!(t, s + shape.half * ENTRY_BYTES);
                assert_eq!(r % 64, 0);
                spans.push((s, s + shape.region_lines * 64));
                spans.push((r, r + shape.row_lines * 64));
            }
            spans.sort();
            for w in spans.windows(2) {
                assert!(w[0].1 <= w[1].0, "P {parts} chunk {chunk}: overlap");
            }
            assert!(spans.last().unwrap().1 <= end);
        }
    }
    for parts in [16usize, 64] {
        let bytes = Shape::new(CHUNK_POSITIONS, parts).lines() * 64;
        assert!(bytes <= (8 << 20) + 64 * 1024, "P {parts}: {bytes}");
    }
    assert!(Shape::new(CHUNK_POSITIONS, 255).lines() * 64 <= (8 << 20) + 400 * 1024);
    // The chunk follows the workers up to 2^19 positions: about P x 0.5 MiB of arena.
    assert_eq!(chunk_positions(0), 1 << 15);
    assert_eq!(chunk_positions(1), 1 << 15);
    assert_eq!(chunk_positions(4), 1 << 17);
    assert_eq!(chunk_positions(5), 5 << 15);
    assert_eq!(chunk_positions(16), CHUNK_POSITIONS);
    assert_eq!(chunk_positions(64), CHUNK_POSITIONS);
    assert_eq!(chunk_positions(usize::MAX), CHUNK_POSITIONS);
    for parts in [2usize, 4, 5, 8, 16] {
        let bytes = Shape::new(chunk_positions(parts), parts).lines() * 64;
        assert!(bytes <= parts * (1 << 19) + 64 * 1024, "P {parts}: {bytes}");
    }
}

// EXPERIMENT (LOSAT_X_BNCSR): the CSR built on a pool, in chunks of any size and with owners
// without cells (more threads than blocks), is the serial build's (`x_bn_csr::serial`: every
// run, the presence vector and the counters) and lists the loop's chains, newest first and in
// reading order; `csr` builds on the pool there and on this thread off it.
//
// Cost: each chunk takes about four passes of the pool's `SpinBarrier`, where a thread spins 2^16
// times before it yields its CPU, so a pool of more threads than CPUs (8 and 40 threads on a CI
// runner's 4) waits at every pass for threads that are not running. The one-position chunks of
// those pools therefore read small cases (`random_case(.., 2, 40)`: 85 to 165 letters, no masked
// interval), the pool of 8 threads with one cell of 256 words (one past a byte's count) instead
// of `x_with_repeats`; the repeat cells are read at the other chunk sizes (on the pool of 40
// threads, in round 2 at chunks of 5,000 positions), and the pools of 2 and 3 threads read every
// case at every chunk size.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn csr_on_a_pool_is_the_serial_build() {
    use super::super::x_bn_csr::tests::{
        chains, pv_params, random_case, x_push_query, x_with_repeats,
    };
    use super::super::x_bn_csr::{serial, x_assert_runs_are_the_chains};
    let mut rng = Lcg(0xB4C5_0003);
    for threads in [2usize, 3, 8, 40] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            let rounds = if threads == 40 { 3 } else { 8 };
            for round in 0..rounds {
                // lut 8 has 32 blocks of 2,048 cells: 40 threads leave owners without cells
                // one table of LOSAT_X_MBCSR's 4^12 cells, newest first, with repeats
                let lut = if threads == 40 {
                    8
                } else if threads == 3 && round == 0 {
                    12
                } else {
                    [8usize, 9, 10, 11][round % 4]
                };
                // both orders for every lut and every chunk size
                let ascending = (round + round / 4) % 2 == 1;
                let chunk = if lut == 12 {
                    chunk_positions(threads)
                } else if threads == 40 && round == 2 {
                    // owners without cells with the repeat cells (lut 8, descending: the
                    // packed head and its map of the words past 255), in a few chunks
                    assert!(lut == 8 && !ascending, "round 2 at 40 threads is packed");
                    5000
                } else {
                    [1usize, 97, 5000, chunk_positions(threads)][(round / 2) % 4]
                };
                // one-position chunks on a pool of more threads than a CI runner's 4 CPUs (above)
                let small = chunk == 1 && threads > 3;
                let mut case = if small {
                    random_case(&mut rng, 2, 40)
                } else {
                    random_case(&mut rng, 2 + round, 2_000 + 900 * round)
                };
                if round % 2 == 0 && !small {
                    // cells of more words than a byte of LOSAT_X_BNCSRPACK's counts holds
                    x_with_repeats(&mut case, lut);
                } else if round % 2 == 0 && threads == 8 {
                    // one cell of 256 words (letter G), each in a chunk of its own
                    x_push_query(&mut case, vec![2u8; 256 + lut - 1]);
                }
                let refs = case.refs();
                let words = case.words(&refs, lut + rng.below(4) as usize, lut, ascending);
                let cells = 1usize << (2 * lut);
                let (pv_len, bts) = pv_params(&case, lut);
                let what = format!("threads {threads} round {round} lut {lut} chunk {chunk}");
                let mut pv = vec![0 as PvArrayType; pv_len];
                let mut helper = vec![0u32; cells >> 11];
                let csr = csr_on_pool(
                    threads,
                    &words,
                    cells,
                    &mut pv,
                    bts,
                    &mut helper,
                    chunk,
                    Vec::new(),
                    None,
                );
                // EXPERIMENT (LOSAT_X_CSRQSKEEP): on a kept `qs` of stale values (a packed head's
                // counts in its first bytes), the same CSR.
                super::super::x_dc_csr::x_assert_same_on_a_kept_qs(
                    &csr,
                    |x_qs| {
                        let mut pv = vec![0 as PvArrayType; pv_len];
                        let mut helper = vec![0u32; cells >> 11];
                        csr_on_pool(
                            threads,
                            &words,
                            cells,
                            &mut pv,
                            bts,
                            &mut helper,
                            chunk,
                            Vec::new(),
                            x_qs,
                        )
                    },
                    &what,
                );
                let mut pv1 = vec![0 as PvArrayType; pv_len];
                let mut helper1 = vec![0u32; cells >> 11];
                let one = serial(&words, cells, &mut pv1, bts, &mut helper1, Vec::new(), None);
                assert_eq!(csr.words(), one.words(), "{what}");
                for cell in 0..cells {
                    assert_eq!(
                        csr.x_run_any(cell),
                        one.x_run_any(cell),
                        "{what} cell {cell}"
                    );
                }
                assert!(pv == pv1 && helper == helper1, "{what}");
                let (hashtable, next_pos, fill_pv, fill_helper) =
                    chains(&case, &words, cells, pv_len, bts);
                x_assert_runs_are_the_chains(
                    &csr,
                    &pv,
                    &helper,
                    &hashtable,
                    &next_pos,
                    &fill_pv,
                    &fill_helper,
                    &what,
                );
            }
            let case = random_case(&mut rng, 2, 3_000);
            let refs = case.refs();
            let words = case.words(&refs, 12, 11, false);
            let (pv_len, bts) = pv_params(&case, 11);
            let mut pv = vec![0 as PvArrayType; pv_len];
            let mut helper = vec![0u32; (1 << 22) >> 11];
            let (_, on_pool) = csr(
                threads,
                &words,
                1 << 22,
                &mut pv,
                bts,
                &mut helper,
                Vec::new(),
                None,
            );
            assert!(on_pool, "threads {threads}");
            Ok(())
        })
        .unwrap();
    }
    let case = random_case(&mut rng, 2, 3_000);
    let refs = case.refs();
    let words = case.words(&refs, 12, 11, false);
    let (pv_len, bts) = pv_params(&case, 11);
    let mut pv = vec![0 as PvArrayType; pv_len];
    let mut helper = vec![0u32; (1 << 22) >> 11];
    let (_, on_pool) = csr(
        4,
        &words,
        1 << 22,
        &mut pv,
        bts,
        &mut helper,
        Vec::new(),
        None,
    );
    assert!(!on_pool, "off a pool");
}
