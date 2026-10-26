//! EXPERIMENT (LOSAT_X_MBFILL): the one-thread word loop of
//! `build_mb_lookup` (the contiguous megablast table and the two-stage table of blastn) with the
//! words read run by run and the update of each applied `X_DELAY` words after a prefetch of its
//! cell, as `x_fill_mb_queued`'s `Delay` queue does (x86_64), or in batches of `X_DELAY` words
//! whose cells are read first, as its `Batch` queue does (wasm32).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1091
//! ```c
//!       for (index = from; index <= last_offset; index++) {
//!          val = *++seq;
//!          /* if an ambiguity is encountered, do not add
//!             any words that would contain it */
//!          if ((val & BLAST2NA_MASK) != 0) {
//!             ecode = 0;
//!             pos = seq + kLutWordLength;
//!             continue;
//!          }
//!
//!          /* get next base */
//!          ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
//!          if (seq < pos)
//!             continue;
//! ...
//!          if (mb_lt->hashtable[ecode] == 0) {
//!             PV_SET(pv_array, ecode, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
//!          mb_lt->hashtable[ecode] = index;
//! ```
//!
//! Where it runs. `build_mb_lookup` calls `fill` where it would call `x_fill_mb_queued`: the
//! table is built on one thread (no pool of two threads or more: PAIRPAR fills the others), without
//! database word counts (`limit_lookup`, which the command line does not offer), on x86_64 and
//! wasm32. On x86_64 the queue prefetches the cell of each word and applies the word `X_DELAY`
//! words later (`run_words`). wasm32 has no prefetch instruction, so there the queue is a batch
//! (`run_words_batch`): `X_DELAY` words are collected, their cells read in one loop (the cache
//! misses overlap), then the words applied in order. Other targets (aarch64) keep the queued loop
//! and its `Batch` queue: the run loop has not been timed there.
//!
//! Why the table is the same. A word is the `lut_word_length` bases that end at a position, and
//! NCBI counts a word only after that many unambiguous bases in a row (`pos` is moved past each
//! ambiguity). This loop reads each run of unambiguous bases of a range once: its first
//! `lut_word_length - 1` bases only build the code, and every later base of the run ends one word,
//! whose code is the last `lut_word_length` bases (`kLutMask` keeps `2 * lut_word_length` bits)
//! and whose `index` is `query_offset + pos - lut_word_length + 2` for the base at `pos`. Those
//! are NCBI's words in NCBI's order. Each update is the statements above on the same cell and
//! index, applied in the same order (the queue is first in, first out: it only prefetches, or
//! reads, the cell of a word before it applies it). `fill` checks once that every index into the
//! tables is in range, so the loop reads them without a check:
//!
//! - `lut_word_length >= 1` and the cell is the code masked with `table_size - 1 = 4^lut_word_length - 1`;
//! - `query_offset + seq.len() < next_pos.len()` for every query (`next_pos` has one more element
//!   than the longest `query_offset + seq.len()`), and an index is at most
//!   `query_offset + seq.len()`;
//! - `(table_size - 1) >> pv_array_bts` is a presence-vector word (NCBI sizes the array so; the
//!   bound check of `pv_set_shift` never fails) and `(table_size - 1) / kCompressionFactor` a
//!   helper counter;
//! - `next_pos` has at most `u32::MAX` elements, so an index fits its `u32`, and the queries have
//!   at most `u32::MAX` letters, so there are fewer than 2^32 words: a helper counter stays below
//!   2^32 and `+ 1` equals `x_mb_apply`'s `saturating_add(1)`.
//!
//! Why it is faster. The queued loop counted the unambiguous bases of every position, carried the
//! database-count test and the bound checks of every index, and kept the queue's state in memory:
//! about 87 instructions a word (line R, R10). Here the loop over a run is a function of its own
//! (`run_words`), so the tables, the queue position and the code stay in registers; the presence
//! bit comes from a table, so the only shift by a variable count (`pv_array_bts`) keeps `cl`; the
//! queue is 64 words deep (32 was 1 % slower, 16 the old depth). On wasm32 the batch keeps the
//! queued loop's order and depth; the gain is the loop over a run and the checks made once. There,
//! applying each word at once, or batching only the first `X_DELAY` words of a run, was faster on
//! long runs on a server core but slower than the batch on a client CPU (short runs, long runs).
//!
//! The variable is read once. `LOSAT_X_MBFILL` (a switch of the strict set,
//! `utils/xswitch.rs`: on unless it is `0`) runs this module.
//!
//! EXPERIMENT (LOSAT_X_MB24): a table of 3-byte elements (`x_mb_narrow.rs`, `XMbTables::x_h24`)
//! runs its own instance of the loops (`H24`): a head is read as the 4 bytes at its element and
//! written back with the fourth byte unchanged, a link is one 4-byte store whose fourth byte is the
//! next element's first byte, still 0 (the words come in increasing order of offset), and the
//! prefetch and the batch read take the element's address. The 4-byte instances are unchanged.

use super::x_mb_lookup_par::Words;
use super::{build_unmasked_ranges, x_mb_narrow, XMbTables, BLAST2NA_MASK, PV_ARRAY_MASK};
#[cfg(test)]
use super::{x_fill_mb_queued, x_mb_queue, XMbQueue};
use std::sync::OnceLock;

/// `kCompressionFactor` (2048), the number of cells that share one helper counter, as a shift.
const K_COMPRESSION_FACTOR_SHIFT: usize = 11;

/// The words between the prefetch of a word's cell and its update (x86_64), or of a batch
/// (`X_BATCH`).
const X_DELAY: usize = 64;

/// The queue is a batch (`run_words_batch`): every target but x86_64, whose queue prefetches
/// (`run_words`); `fill` takes x86_64 and wasm32 only.
const X_BATCH: bool = !cfg!(target_arch = "x86_64");

/// No NCBI counterpart: reads LOSAT_X_MBFILL once (a switch of the strict set: on unless it is
/// `0`; off runs `x_fill_mb_queued`). It chooses a LOSAT code path and changes no value NCBI
/// computes.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBFILL"))
}

/// Dispatch of `build_mb_lookup`: fills `t` as `x_fill_mb_queued` would and returns true, or
/// returns false (switch off, or a build this module does not take) and leaves `t` as it was.
pub(super) fn fill(
    words: &Words<'_>,
    db_word_counts: Option<&[u8]>,
    t: &mut XMbTables<'_>,
) -> bool {
    if !on() {
        return false;
    }
    let taken = cfg!(any(target_arch = "x86_64", target_arch = "wasm32"))
        && db_word_counts.is_none()
        && in_range(words, t);
    crate::utils::xreach::x_reach!("LOSAT_X_MBFILL", taken);
    if !taken {
        return false;
    }
    // SAFETY: `in_range(words, t)` holds.
    let (n_words, ambiguous) = unsafe {
        if t.x_h24 {
            // (`in_range`: a 3-byte table is newest first)
            fill_runs::<false, X_BATCH, true>(words, t)
        } else if words.ascending_cells {
            fill_runs::<true, X_BATCH, false>(words, t)
        } else {
            fill_runs::<false, X_BATCH, false>(words, t)
        }
    };
    *t.total_positions += n_words;
    *t.ambiguous_skipped += ambiguous;
    true
}

/// The bounds the loop relies on (module documentation). A 3-byte table (`x_h24`) holds
/// `x_mb_narrow::capacity` elements in its words, values below 2^24, and is newest first.
fn in_range(words: &Words<'_>, t: &XMbTables<'_>) -> bool {
    let lut_word_length = words.lut_word_length;
    let last_cell = words.kmer_mask as usize;
    let letters = words
        .queries_blastna
        .iter()
        .fold(0usize, |sum, seq| sum.saturating_add(seq.len()));
    // (3-byte elements: indexes and values below 2^24 as well)
    let (cells, next_len) = if t.x_h24 {
        (
            x_mb_narrow::capacity(t.hashtable.len()),
            x_mb_narrow::capacity(t.next_pos.len()).min(x_mb_narrow::CELLS),
        )
    } else {
        (t.hashtable.len(), t.next_pos.len())
    };
    (1..=16).contains(&lut_word_length)
        && words.kmer_mask == (1u64 << (2 * lut_word_length)) - 1
        && if t.x_h24 {
            cells > last_cell && !words.ascending_cells
        } else {
            cells == last_cell.wrapping_add(1)
        }
        && (last_cell >> t.pv_array_bts) < t.pv_array.len()
        && (last_cell >> K_COMPRESSION_FACTOR_SHIFT) < t.helper_array.len()
        && next_len <= u32::MAX as usize
        && letters <= u32::MAX as usize
        && words
            .queries_blastna
            .iter()
            .zip(words.query_offsets)
            .all(|(seq, &offset)| (offset.max(0) as usize).saturating_add(seq.len()) < next_len)
}

/// `1 << k`, the bit of cell `k` (mod 32) in its presence-vector word, read from a table: a shift
/// by a variable count needs `cl` on x86-64 without BMI2, and the loop keeps `pv_array_bts` there.
static X_PV_BIT: [u32; 32] = {
    let mut bits = [0u32; 32];
    let mut k = 0;
    while k < 32 {
        bits[k] = 1 << k;
        k += 1;
    }
    bits
};

/// The tables of `XMbTables` as pointers, so the word loop keeps them in registers.
#[derive(Clone, Copy)]
struct Raw {
    hashtable: *mut u32,
    next_pos: *mut u32,
    pv_array: *mut u32,
    pv_array_bts: usize,
    helper_array: *mut u32,
}

impl Raw {
    fn new(t: &mut XMbTables<'_>) -> Self {
        Self {
            hashtable: t.hashtable.as_mut_ptr(),
            next_pos: t.next_pos.as_mut_ptr(),
            pv_array: t.pv_array.as_mut_ptr(),
            pv_array_bts: t.pv_array_bts,
            helper_array: t.helper_array.as_mut_ptr(),
        }
    }

    /// The update of one word (`x_mb_apply` with the indexes `in_range` checked); `cells`: the
    /// words in the order applied (`ASC`). `H24`: on 3-byte elements (module documentation).
    /// SAFETY: `in_range` held for the tables of `self`, and `bucket` / `q_off_1` are a word, the
    /// words in increasing order of `q_off_1` with `H24`.
    #[inline(always)]
    unsafe fn apply<const ASC: bool, const H24: bool>(
        self,
        cells: &mut Vec<(u32, u32)>,
        bucket: usize,
        q_off_1: u32,
    ) {
        unsafe {
            let word = if H24 {
                x_mb_narrow::load_word(self.hashtable, bucket)
            } else {
                0
            };
            let old = if H24 {
                word & x_mb_narrow::MASK
            } else {
                *self.hashtable.add(bucket)
            };
            if old == 0 {
                *self.pv_array.add(bucket >> self.pv_array_bts) |= X_PV_BIT[bucket & PV_ARRAY_MASK];
            } else {
                let h = self.helper_array.add(bucket >> K_COMPRESSION_FACTOR_SHIFT);
                *h = (*h).wrapping_add(1);
            }
            if H24 {
                // the link's fourth byte is the next element's first, not written yet
                x_mb_narrow::store_word(self.next_pos, q_off_1 as usize, old);
                x_mb_narrow::store_word(
                    self.hashtable,
                    bucket,
                    (word & !x_mb_narrow::MASK) | q_off_1,
                );
            } else {
                *self.next_pos.add(q_off_1 as usize) = old;
                *self.hashtable.add(bucket) = q_off_1;
            }
        }
        if ASC {
            cells.push((bucket as u32, q_off_1));
        }
    }
}

/// The last `X_DELAY` words read, each applied when it leaves (`n` words entered so far).
struct Queue {
    ring: [(u32, u32); X_DELAY],
    n: usize,
}

/// The fill: the words of `build_mb_lookup` in its order, read run by run (module documentation),
/// each update applied `X_DELAY` words later, or in batches of `X_DELAY` (`BATCH`). Returns the
/// number of words and of ambiguous bases.
/// SAFETY: `in_range(words, t)` holds.
#[inline(never)]
unsafe fn fill_runs<const ASC: bool, const BATCH: bool, const H24: bool>(
    words: &Words<'_>,
    t: &mut XMbTables<'_>,
) -> (usize, usize) {
    let raw = Raw::new(t);
    let cells = &mut *t.words;
    let mut queue = Queue {
        ring: [(0u32, 0u32); X_DELAY],
        n: 0,
    };
    // `BATCH`: the cells read ahead, xor-ed (kept so that the reads stay)
    let mut touched = 0u32;
    let lut_word_length = words.lut_word_length;
    let (mut n_words, mut ambiguous) = (0usize, 0usize);
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        if seq.len() < lut_word_length {
            continue;
        }
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let query_offset = words.query_offsets[q_idx].max(0) as usize;
        for (range_start, range_end) in build_unmasked_ranges(seq.len(), masks) {
            if words.word_length > range_end.saturating_sub(range_start) {
                continue;
            }
            let s = &seq[range_start..range_end];
            // the index of the word that ends at `s[i]` is `index_base + i`
            let index_base = (query_offset + range_start + 2).wrapping_sub(lut_word_length);
            let mut i = 0usize;
            'run: while i < s.len() {
                // the first `lut_word_length - 1` bases of a run
                let mut code = 0u64;
                let mut got = 1usize;
                while got < lut_word_length {
                    let Some(&base) = s.get(i) else {
                        break 'run;
                    };
                    i += 1;
                    if (base & BLAST2NA_MASK) != 0 {
                        code = 0;
                        got = 1;
                        ambiguous += 1;
                        continue;
                    }
                    code = (code << 2) | base as u64;
                    got += 1;
                }
                // every further base of the run ends a word
                let first = i;
                // SAFETY: `in_range` holds; the words are those of `words`.
                i = unsafe {
                    if BATCH {
                        run_words_batch::<ASC, H24>(
                            s,
                            i,
                            code,
                            words.kmer_mask,
                            index_base,
                            &mut queue,
                            &mut touched,
                            raw,
                            cells,
                        )
                    } else {
                        run_words::<ASC, H24>(
                            s,
                            i,
                            code,
                            words.kmer_mask,
                            index_base,
                            &mut queue,
                            raw,
                            cells,
                        )
                    }
                };
                n_words += i - first;
                if i < s.len() {
                    // `s[i]` is ambiguous
                    ambiguous += 1;
                    i += 1;
                }
            }
        }
    }
    // (with `BATCH`, `queue.n < X_DELAY`: the words of the last batch, in order)
    for k in queue.n - queue.n.min(X_DELAY)..queue.n {
        let (b, q) = queue.ring[k % X_DELAY];
        // SAFETY: as above.
        unsafe { raw.apply::<ASC, H24>(cells, b as usize, q) };
    }
    if BATCH {
        std::hint::black_box(touched);
    }
    (n_words, ambiguous)
}

/// The words that end at `s[i..]` up to the first ambiguous base (whose position it returns, or
/// `s.len()`), `code` holding the `lut_word_length - 1` bases before `s[i]`: each is queued, and the
/// word it pushes out of the queue applied. Its own function, so its loop has the registers.
/// SAFETY: `in_range` held for the tables of `raw`; the words are words of the fill.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn run_words<const ASC: bool, const H24: bool>(
    s: &[u8],
    mut i: usize,
    mut code: u64,
    kmer_mask: u64,
    index_base: usize,
    queue: &mut Queue,
    raw: Raw,
    cells: &mut Vec<(u32, u32)>,
) -> usize {
    let mut ring = queue.ring;
    let mut n = queue.n;
    // until the queue is full no word leaves it
    while n < X_DELAY {
        let Some(&base) = s.get(i) else {
            break;
        };
        if (base & BLAST2NA_MASK) != 0 {
            break;
        }
        code = ((code << 2) | base as u64) & kmer_mask;
        x_prefetch::<H24>(raw.hashtable, code as usize);
        ring[n] = (code as u32, index_base.wrapping_add(i) as u32);
        n += 1;
        i += 1;
    }
    if n >= X_DELAY {
        while let Some(&base) = s.get(i) {
            if (base & BLAST2NA_MASK) != 0 {
                break;
            }
            code = ((code << 2) | base as u64) & kmer_mask;
            x_prefetch::<H24>(raw.hashtable, code as usize);
            let slot = n % X_DELAY;
            let (b, q) = ring[slot];
            // SAFETY: `(b, q)` is a word of the fill.
            unsafe { raw.apply::<ASC, H24>(cells, b as usize, q) };
            ring[slot] = (code as u32, index_base.wrapping_add(i) as u32);
            n += 1;
            i += 1;
        }
    }
    queue.ring = ring;
    queue.n = n;
    i
}

/// `run_words` with the queue as a batch (`X_BATCH`): each word goes to `queue.ring[queue.n]`;
/// when `X_DELAY` are there, their cells are read in one loop (so the cache misses overlap without
/// a prefetch instruction; the values go to `touched`), then the words are applied in order and the
/// batch starts again. On return `queue.n < X_DELAY`. The ring stays in `queue`: on wasm32 a local
/// copy would be in memory as well.
/// SAFETY: as `run_words`.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn run_words_batch<const ASC: bool, const H24: bool>(
    s: &[u8],
    mut i: usize,
    mut code: u64,
    kmer_mask: u64,
    index_base: usize,
    queue: &mut Queue,
    touched: &mut u32,
    raw: Raw,
    cells: &mut Vec<(u32, u32)>,
) -> usize {
    let mut n = queue.n;
    while let Some(&base) = s.get(i) {
        if (base & BLAST2NA_MASK) != 0 {
            break;
        }
        code = ((code << 2) | base as u64) & kmer_mask;
        queue.ring[n] = (code as u32, index_base.wrapping_add(i) as u32);
        n += 1;
        i += 1;
        if n == X_DELAY {
            let mut read = *touched;
            for &(b, _) in &queue.ring {
                // SAFETY: `b` is the cell of a word of the fill.
                read ^= unsafe {
                    if H24 {
                        x_mb_narrow::load_word(raw.hashtable, b as usize)
                    } else {
                        *raw.hashtable.add(b as usize)
                    }
                };
            }
            *touched = read;
            for &(b, q) in &queue.ring {
                // SAFETY: `(b, q)` is a word of the fill.
                unsafe { raw.apply::<ASC, H24>(cells, b as usize, q) };
            }
            n = 0;
        }
    }
    queue.n = n;
    i
}

/// A prefetch of cell `bucket` of the table at `hashtable` (x86_64; nothing elsewhere); `H24`: of
/// its 3-byte element.
#[inline(always)]
fn x_prefetch<const H24: bool>(hashtable: *mut u32, bucket: usize) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not fault.
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            if H24 {
                x_mb_narrow::at(hashtable, bucket)
            } else {
                hashtable.wrapping_add(bucket) as *const i8
            },
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (hashtable, bucket);
}

#[cfg(test)]
mod tests {
    use super::super::compute_mb_pv_params;
    use super::super::x_mb_lookup_par::tests::{random_case, repeat_case, Lcg};
    use super::*;

    #[derive(PartialEq, Debug)]
    struct Built {
        hashtable: Vec<u32>,
        next_pos: Vec<u32>,
        pv_array: Vec<u32>,
        helper_array: Vec<u32>,
        words: Vec<(u32, u32)>,
        total_positions: usize,
        ambiguous_skipped: usize,
    }

    /// The arrays a fill makes; `h24`: in 3-byte elements (LOSAT_X_MB24), returned as 4-byte values.
    fn build(words: &Words<'_>, h24: bool, fill: impl FnOnce(&mut XMbTables<'_>)) -> Built {
        let table_size = words.kmer_mask as usize + 1;
        let entries: usize = words.queries_blastna.iter().map(|q| q.len()).sum();
        let (pv_size, pv_array_bts) =
            compute_mb_pv_params(table_size, entries, false, words.lut_word_length);
        let next_len = words
            .queries_blastna
            .iter()
            .zip(words.query_offsets)
            .map(|(q, &o)| o.max(0) as usize + q.len())
            .max()
            .unwrap_or(0)
            + 1;
        let width = |n: usize| if h24 { x_mb_narrow::words_for(n) } else { n };
        let mut b = Built {
            hashtable: vec![0; width(table_size)],
            next_pos: vec![0; width(next_len)],
            pv_array: vec![0; pv_size],
            helper_array: vec![0; (table_size >> K_COMPRESSION_FACTOR_SHIFT).max(1)],
            words: Vec::new(),
            total_positions: 0,
            ambiguous_skipped: 0,
        };
        fill(&mut XMbTables {
            hashtable: &mut b.hashtable,
            next_pos: &mut b.next_pos,
            pv_array: &mut b.pv_array,
            pv_array_bts,
            helper_array: &mut b.helper_array,
            words: &mut b.words,
            total_positions: &mut b.total_positions,
            ambiguous_skipped: &mut b.ambiguous_skipped,
            x_h24: h24,
        });
        if h24 {
            b.hashtable = x_mb_narrow::widen(&b.hashtable, table_size);
            b.next_pos = x_mb_narrow::widen(&b.next_pos, next_len);
        }
        b
    }

    /// Both queues (prefetch and batch) and both cell orders build the arrays and counters of
    /// `x_fill_mb_queued` on random queries (ambiguity codes, masks, several queries,
    /// `lut_word_length` 4..12, runs shorter and longer than the queue), on single queries with
    /// fewer words than the queue holds (`lut_word_length` 6, 7, 10, 12) and on long repeats (long
    /// chains, helper counts); newest first, so do the instances for 3-byte elements (LOSAT_X_MB24)
    /// of both and of the queued loop.
    #[test]
    fn the_run_loop_builds_the_queued_loops_arrays() {
        let mut rng = Lcg(0x4D42_F111);
        for case in 0..60 {
            let (queries, offsets, masks, word_length, lut_word_length) = if case % 6 == 4 {
                // fewer words than the queue holds, at the table widths `random_case` leaves out
                let lut = [6usize, 7, 10, 12][(case / 6) % 4];
                let len = lut + rng.below(60) as usize;
                let seq: Vec<u8> = (0..len)
                    .map(|_| {
                        if rng.below(30) == 0 {
                            14
                        } else {
                            rng.below(4) as u8
                        }
                    })
                    .collect();
                (vec![seq], vec![0], vec![Vec::new()], lut, lut)
            } else if case % 6 == 5 {
                let lut = [11usize, 12][case % 2];
                (
                    vec![repeat_case(&mut rng, 9000)],
                    vec![0],
                    vec![Vec::new()],
                    lut + 16,
                    lut,
                )
            } else {
                random_case(&mut rng)
            };
            let refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
            for ascending_cells in [false, true] {
                let words = Words {
                    queries_blastna: &refs,
                    query_offsets: &offsets,
                    query_masks: &masks,
                    word_length,
                    lut_word_length,
                    kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                    ascending_cells,
                };
                let reference = build(&words, false, |t| {
                    x_fill_mb_queued::<false>(
                        x_mb_queue(),
                        words.queries_blastna,
                        words.query_offsets,
                        words.query_masks,
                        word_length,
                        lut_word_length,
                        words.kmer_mask,
                        None,
                        0,
                        ascending_cells,
                        t,
                    )
                });
                assert!(reference.total_positions > 0 || queries.iter().all(|q| q.len() < 80));
                for (batch, h24) in [(false, false), (true, false), (false, true), (true, true)] {
                    if h24 && ascending_cells {
                        continue;
                    }
                    let got = build(&words, h24, |t| {
                        assert!(in_range(&words, t));
                        // SAFETY: checked above.
                        let (n, a) = unsafe {
                            match (ascending_cells, batch, h24) {
                                (true, false, _) => fill_runs::<true, false, false>(&words, t),
                                (false, false, false) => {
                                    fill_runs::<false, false, false>(&words, t)
                                }
                                (true, true, _) => fill_runs::<true, true, false>(&words, t),
                                (false, true, false) => fill_runs::<false, true, false>(&words, t),
                                (false, false, true) => fill_runs::<false, false, true>(&words, t),
                                (false, true, true) => fill_runs::<false, true, true>(&words, t),
                            }
                        };
                        *t.total_positions += n;
                        *t.ambiguous_skipped += a;
                    });
                    assert!(
                        got == reference,
                        "case {case} ascending {ascending_cells} batch {batch} h24 {h24}"
                    );
                }
                if !ascending_cells {
                    for queue in [XMbQueue::Delay, XMbQueue::Batch(64)] {
                        let got = build(&words, true, |t| {
                            x_fill_mb_queued::<true>(
                                queue,
                                words.queries_blastna,
                                words.query_offsets,
                                words.query_masks,
                                word_length,
                                lut_word_length,
                                words.kmer_mask,
                                None,
                                0,
                                false,
                                t,
                            )
                        });
                        assert!(got == reference, "case {case} queued h24");
                    }
                }
            }
        }
    }
}
