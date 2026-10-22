//! EXPERIMENT (LOSAT_X_MBFILL, shadow LOSAT_X_MBFILLSHADOW): the one-thread word loop of
//! `build_mb_lookup` (the contiguous megablast table and the two-stage table of blastn) with the
//! words read by a run-based loop and the update queued as `x_fill_mb_queued`'s `Delay` does.
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
//! table is built on one thread (no pool of two threads or more, or LOSAT_X_PAIRPAR off), without
//! database word counts (`-limit_lookup`), on x86_64 (the queue needs a prefetch instruction).
//!
//! Why the table is the same. A word is the `lut_word_length` bases that end at a position, and
//! NCBI counts a word only after that many unambiguous bases in a row (`pos` is moved past each
//! ambiguity). This loop reads each run of unambiguous bases of a range once: its first
//! `lut_word_length - 1` bases only build the code, and every later base of the run ends one word,
//! whose code is the last `lut_word_length` bases (`kLutMask` keeps `2 * lut_word_length` bits)
//! and whose `index` is `query_offset + pos - lut_word_length + 2` for the base at `pos`. Those
//! are NCBI's words in NCBI's order. Each update is the statements above on the same cell and
//! index, applied in the same order (the queue only prefetches the cell of a word `X_DELAY` words
//! before it applies it). `fill` checks once that every index into the tables is in range, so the
//! loop reads them without a check:
//!
//! - the cell is the code masked with `table_size - 1`;
//! - `query_offset + seq.len() < next_pos.len()` for every query (`next_pos` has one more element
//!   than the longest `query_offset + seq.len()`), and an index is at most
//!   `query_offset + seq.len()`;
//! - `(table_size - 1) >> pv_array_bts` is a presence-vector word (NCBI sizes the array so; the
//!   bound check of `pv_set_shift` never fails) and `(table_size - 1) / kCompressionFactor` a
//!   helper counter;
//! - fewer than 2^32 words, since `next_pos` has at most `u32::MAX` elements, so a helper counter
//!   stays below 2^32 and `+ 1` equals `x_mb_apply`'s `saturating_add(1)`.
//!
//! The variables are read once. `LOSAT_X_MBFILL` (a switch of the strict set,
//! `utils/xswitch.rs`: on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set) runs this
//! module. `LOSAT_X_MBFILLSHADOW` (set to any value; not in the strict set; implies it) runs this
//! module, builds the arrays again with `x_fill_mb_queued` and panics unless every array and
//! counter is equal; `main` prints one `[X_SHADOW] MBFILL compared=N` line.

use super::x_mb_lookup_par::Words;
use super::{
    build_unmasked_ranges, x_fill_mb_queued, x_mb_apply, x_mb_queue, XMbTables, BLAST2NA_MASK,
    PV_ARRAY_MASK,
};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// `kCompressionFactor` (2048), the number of cells that share one helper counter, as a shift.
const K_COMPRESSION_FACTOR_SHIFT: usize = 11;

// No NCBI counterpart: the number of builds the shadow compared (relaxed; it orders nothing).
static X_SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: reads LOSAT_X_MBFILL / LOSAT_X_MBFILLSHADOW once. 0 = `x_fill_mb_queued`,
/// 1 = this module, 2 = both and compare; SHADOW wins. It chooses a LOSAT code path and changes no
/// value NCBI computes.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_MBFILLSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_MBFILL") {
            1
        } else {
            0
        }
    })
}

/// PROBE (not for merge): the loop shape, LOSAT_X_MBFILLARM. 1 = the word loop of
/// `x_fill_mb_queued` with a closure and a delay of 32; 2 = the run loop, delay 32; 3 = the run
/// loop, delay 32, update without the branch on the old head; 4 = the run loop, delay 64;
/// 5 = the run loop, delay 64, without the branch.
fn arm() -> u8 {
    static ARM: OnceLock<u8> = OnceLock::new();
    *ARM.get_or_init(|| {
        std::env::var("LOSAT_X_MBFILLARM")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(2)
    })
}

/// The summary line `main` prints at exit with LOSAT_X_MBFILLSHADOW.
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    if mode() != 2 {
        return;
    }
    let compared = X_SHADOW_COMPARED.load(Ordering::Relaxed);
    eprintln!(
        "[X_SHADOW] MBFILL compared={compared} ({})",
        if compared == 0 {
            "nothing compared: no one-thread megablast table was built"
        } else {
            "every table identical"
        }
    );
}

/// Dispatch of `build_mb_lookup`: fills `t` as `x_fill_mb_queued` would and returns true, or
/// returns false (switch off, or a build this module does not take) and leaves `t` as it was.
pub(super) fn fill(
    words: &Words<'_>,
    db_word_counts: Option<&[u8]>,
    t: &mut XMbTables<'_>,
) -> bool {
    let mode = mode();
    if mode == 0 {
        return false;
    }
    let taken = cfg!(target_arch = "x86_64") && db_word_counts.is_none() && in_range(words, t);
    crate::utils::xreach::x_reach!("LOSAT_X_MBFILL", taken);
    if !taken {
        return false;
    }
    let (n_words, ambiguous) = fill_arm(arm(), words, t);
    *t.total_positions += n_words;
    *t.ambiguous_skipped += ambiguous;
    if mode == 2 {
        shadow_check(words, t);
    }
    true
}

/// PROBE: the fill of arm `arm` (see `arm`); returns the number of words and of ambiguous bases.
/// `in_range(words, t)` must hold.
fn fill_arm(arm: u8, words: &Words<'_>, t: &mut XMbTables<'_>) -> (usize, usize) {
    // SAFETY (every arm): the caller checked `in_range(words, t)`.
    unsafe {
        match (arm, words.ascending_cells) {
            (1, false) => fill_delay::<32, false>(words, t),
            (1, true) => fill_delay::<32, true>(words, t),
            (3, false) => fill_runs::<32, true, false>(words, t),
            (3, true) => fill_runs::<32, true, true>(words, t),
            (4, false) => fill_runs::<64, false, false>(words, t),
            (4, true) => fill_runs::<64, false, true>(words, t),
            (5, false) => fill_runs::<64, true, false>(words, t),
            (5, true) => fill_runs::<64, true, true>(words, t),
            (_, false) => fill_runs::<32, false, false>(words, t),
            (_, true) => fill_runs::<32, false, true>(words, t),
        }
    }
}

/// The bounds the loop relies on (module documentation).
fn in_range(words: &Words<'_>, t: &XMbTables<'_>) -> bool {
    let last_cell = words.kmer_mask as usize;
    t.hashtable.len() == last_cell.wrapping_add(1)
        && (last_cell >> t.pv_array_bts) < t.pv_array.len()
        && (last_cell >> K_COMPRESSION_FACTOR_SHIFT) < t.helper_array.len()
        && t.next_pos.len() <= u32::MAX as usize
        && words
            .queries_blastna
            .iter()
            .zip(words.query_offsets)
            .all(|(seq, &offset)| {
                (offset.max(0) as usize).saturating_add(seq.len()) < t.next_pos.len()
            })
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

    /// The update of one word (`x_mb_apply` with the indexes `in_range` checked). `BF`: the
    /// presence bit and the helper counter are always written, with 0 where the C branch skips
    /// them. `cells`: the words in the order applied (`ASC`).
    /// SAFETY: `in_range` held for the tables of `self`, and `bucket` / `q_off_1` are a word.
    #[inline(always)]
    unsafe fn apply<const BF: bool, const ASC: bool>(
        self,
        cells: &mut Vec<(u32, u32)>,
        bucket: usize,
        q_off_1: u32,
    ) {
        let bit = X_PV_BIT[bucket & PV_ARRAY_MASK];
        unsafe {
            let old = *self.hashtable.add(bucket);
            if BF {
                let empty = (old == 0) as u32;
                *self.pv_array.add(bucket >> self.pv_array_bts) |= bit & empty.wrapping_neg();
                let h = self.helper_array.add(bucket >> K_COMPRESSION_FACTOR_SHIFT);
                *h = (*h).wrapping_add(empty ^ 1);
            } else if old == 0 {
                *self.pv_array.add(bucket >> self.pv_array_bts) |= bit;
            } else {
                let h = self.helper_array.add(bucket >> K_COMPRESSION_FACTOR_SHIFT);
                *h = (*h).wrapping_add(1);
            }
            *self.next_pos.add(q_off_1 as usize) = old;
            *self.hashtable.add(bucket) = q_off_1;
        }
        if ASC {
            cells.push((bucket as u32, q_off_1));
        }
    }
}

/// The last `D` words read, each applied when it leaves (`n` words entered so far).
struct Queue<const D: usize> {
    ring: [(u32, u32); D],
    n: usize,
}

/// PROBE (arm 1): the word loop of `x_fill_mb_queued` with the update of each word applied `D`
/// words later, after a prefetch of its cell. Returns the number of words and of ambiguous bases.
#[inline(never)]
fn fill_delay<const D: usize, const ASC: bool>(
    words: &Words<'_>,
    t: &mut XMbTables<'_>,
) -> (usize, usize) {
    let mut ring = [(0u32, 0u32); D];
    let mut n = 0usize;
    let mut queue = |t: &mut XMbTables<'_>, bucket: usize, q_off_1: u32| {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: `bucket` < table size (masked code); a prefetch does not fault.
        unsafe {
            core::arch::x86_64::_mm_prefetch(
                t.hashtable.as_ptr().wrapping_add(bucket) as *const i8,
                core::arch::x86_64::_MM_HINT_T0,
            );
        }
        let slot = n % D;
        if n >= D {
            let (b, q) = ring[slot];
            x_mb_apply(t, ASC, b as usize, q);
        }
        ring[slot] = (bucket as u32, q_off_1);
        n += 1;
    };
    let counts = for_each_word_loop(words, |bucket, q_off_1| queue(t, bucket, q_off_1));
    for k in n - n.min(D)..n {
        let (b, q) = ring[k % D];
        x_mb_apply(t, ASC, b as usize, q);
    }
    counts
}

/// The fill: the words of `build_mb_lookup` in its order, read run by run (module documentation),
/// each update applied `D` words later, after a prefetch of its cell. Returns the number of words
/// and of ambiguous bases.
/// SAFETY: `in_range(words, t)` holds.
#[inline(never)]
unsafe fn fill_runs<const D: usize, const BF: bool, const ASC: bool>(
    words: &Words<'_>,
    t: &mut XMbTables<'_>,
) -> (usize, usize) {
    let raw = Raw::new(t);
    let cells = &mut *t.words;
    let mut queue = Queue {
        ring: [(0u32, 0u32); D],
        n: 0,
    };
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
                    run_words::<D, BF, ASC>(
                        s,
                        i,
                        code,
                        words.kmer_mask,
                        index_base,
                        &mut queue,
                        raw,
                        cells,
                    )
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
    for k in queue.n - queue.n.min(D)..queue.n {
        let (b, q) = queue.ring[k % D];
        // SAFETY: as above.
        unsafe { raw.apply::<BF, ASC>(cells, b as usize, q) };
    }
    (n_words, ambiguous)
}

/// The words that end at `s[i..]` up to the first ambiguous base (whose position it returns, or
/// `s.len()`), `code` holding the `lut_word_length - 1` bases before `s[i]`: each is queued, and the
/// word it pushes out of the queue applied. Its own function, so its loop has the registers.
/// SAFETY: `in_range` held for the tables of `raw`; the words are words of the fill.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn run_words<const D: usize, const BF: bool, const ASC: bool>(
    s: &[u8],
    mut i: usize,
    mut code: u64,
    kmer_mask: u64,
    index_base: usize,
    queue: &mut Queue<D>,
    raw: Raw,
    cells: &mut Vec<(u32, u32)>,
) -> usize {
    let mut ring = queue.ring;
    let mut n = queue.n;
    // until the queue is full no word leaves it
    while n < D {
        let Some(&base) = s.get(i) else {
            break;
        };
        if (base & BLAST2NA_MASK) != 0 {
            break;
        }
        code = ((code << 2) | base as u64) & kmer_mask;
        x_prefetch(raw.hashtable, code as usize);
        ring[n] = (code as u32, index_base.wrapping_add(i) as u32);
        n += 1;
        i += 1;
    }
    if n >= D {
        while let Some(&base) = s.get(i) {
            if (base & BLAST2NA_MASK) != 0 {
                break;
            }
            code = ((code << 2) | base as u64) & kmer_mask;
            x_prefetch(raw.hashtable, code as usize);
            let slot = n % D;
            let (b, q) = ring[slot];
            // SAFETY: `(b, q)` is a word of the fill.
            unsafe { raw.apply::<BF, ASC>(cells, b as usize, q) };
            ring[slot] = (code as u32, index_base.wrapping_add(i) as u32);
            n += 1;
            i += 1;
        }
    }
    queue.ring = ring;
    queue.n = n;
    i
}

/// A prefetch of cell `bucket` of the table at `hashtable` (x86_64; nothing elsewhere).
#[inline(always)]
fn x_prefetch(hashtable: *mut u32, bucket: usize) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not fault.
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            hashtable.wrapping_add(bucket) as *const i8,
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (hashtable, bucket);
}

/// PROBE (arm 1): the word loop of `x_fill_mb_queued` without the database word counts.
#[inline(always)]
fn for_each_word_loop(words: &Words<'_>, mut f: impl FnMut(usize, u32)) -> (usize, usize) {
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
            let mut current_kmer: u64 = 0;
            let mut valid_bases = 0usize;
            for pos in range_start..range_end {
                let base = seq[pos];
                if (base & BLAST2NA_MASK) != 0 {
                    current_kmer = 0;
                    valid_bases = 0;
                    ambiguous += 1;
                    continue;
                }
                current_kmer = ((current_kmer << 2) | base as u64) & words.kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length {
                    continue;
                }
                n_words += 1;
                f(
                    current_kmer as usize,
                    (query_offset + (pos + 1 - lut_word_length) + 1) as u32,
                );
            }
        }
    }
    (n_words, ambiguous)
}

/// LOSAT_X_MBFILLSHADOW: builds the arrays again with `x_fill_mb_queued` on fresh tables and
/// panics at the first one that differs from what `fill` built.
#[cold]
#[inline(never)]
fn shadow_check(words: &Words<'_>, t: &XMbTables<'_>) {
    let mut hashtable = vec![0u32; t.hashtable.len()];
    let mut next_pos = vec![0u32; t.next_pos.len()];
    let mut pv_array = vec![0u32; t.pv_array.len()];
    let mut helper_array = vec![0u32; t.helper_array.len()];
    let mut ref_words = Vec::new();
    let (mut total_positions, mut ambiguous_skipped) = (0usize, 0usize);
    x_fill_mb_queued(
        x_mb_queue(),
        words.queries_blastna,
        words.query_offsets,
        words.query_masks,
        words.word_length,
        words.lut_word_length,
        words.kmer_mask,
        None,
        0,
        words.ascending_cells,
        &mut XMbTables {
            hashtable: &mut hashtable,
            next_pos: &mut next_pos,
            pv_array: &mut pv_array,
            pv_array_bts: t.pv_array_bts,
            helper_array: &mut helper_array,
            words: &mut ref_words,
            total_positions: &mut total_positions,
            ambiguous_skipped: &mut ambiguous_skipped,
        },
    );
    let first_difference = |a: &[u32], b: &[u32]| a.iter().zip(b).position(|(x, y)| x != y);
    for (name, new, old) in [
        ("hashtable", &*t.hashtable, &hashtable[..]),
        ("next_pos", &*t.next_pos, &next_pos[..]),
        ("pv_array", &*t.pv_array, &pv_array[..]),
        ("helper_array", &*t.helper_array, &helper_array[..]),
    ] {
        if let Some(i) = first_difference(new, old) {
            panic!(
                "[X_SHADOW] MBFILL {name}[{i}] = {} differs from {}",
                new[i], old[i]
            );
        }
    }
    assert!(
        *t.words == ref_words
            && *t.total_positions == total_positions
            && *t.ambiguous_skipped == ambiguous_skipped,
        "[X_SHADOW] MBFILL words or counters differ"
    );
    X_SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
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

    fn build(words: &Words<'_>, fill: impl FnOnce(&mut XMbTables<'_>)) -> Built {
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
        let mut b = Built {
            hashtable: vec![0; table_size],
            next_pos: vec![0; next_len],
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
        });
        b
    }

    /// Every arm and both cell orders build the arrays and counters of `x_fill_mb_queued` on random
    /// queries (ambiguity codes, masks, several queries, `lut_word_length` 4..12) and on long
    /// repeats (long chains, helper counts).
    #[test]
    fn every_arm_builds_the_queued_loops_arrays() {
        let mut rng = Lcg(0x4D42_F111);
        for case in 0..60 {
            let (queries, offsets, masks, word_length, lut_word_length) = if case % 6 == 5 {
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
                let reference = build(&words, |t| {
                    x_fill_mb_queued(
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
                assert!(reference.total_positions > 0 || queries.iter().all(|q| q.len() < 40));
                for arm in 1..=5 {
                    let got = build(&words, |t| {
                        assert!(in_range(&words, t));
                        let (n, a) = fill_arm(arm, &words, t);
                        *t.total_positions += n;
                        *t.ambiguous_skipped += a;
                    });
                    assert!(
                        got == reference,
                        "case {case} arm {arm} ascending {ascending_cells}"
                    );
                }
            }
        }
    }
}
