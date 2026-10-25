//! PROBE (LOSAT_X_MBHRANK=1, branch exp/mb-head-rank-probe; v4-A10 s20, MF E2; not a PR): the
//! head table of a contiguous megablast table of 4^12 cells (`build_mb_lookup`) with one element
//! per non-empty cell instead of one per cell. NCBI's `hashtable` (4^12 x 4 B = 64 MiB) becomes
//! `ranks` (one u64 per 32 cells: the number of non-empty cells before them in the high half,
//! their presence bits in the low half; 4 MiB) and `heads` (one u32 per non-empty cell, in cell
//! order). `head(cell)` is `hashtable[cell]`: 0 for an empty cell, otherwise
//! `heads[rank(cell)]`. `next_pos` is NCBI's.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
//! ```c
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
//! The fill (`serial`, one thread). Pass 1 reads the words in the fill's order (`each_word`) and
//! stores `cell + 1` at `next_pos[index]`, sets the cell's presence bit and adds 1 to its helper
//! counter. Pass 2 turns the bits into running counts and makes the presence vector (a cell's PV
//! bit is set at its first word: every non-empty cell) and the helper counters (the fill adds 1 at
//! every further word: words minus non-empty cells of each counter's block). Pass 3 walks
//! `next_pos` in index order, which is the fill's order when the queries' offsets ascend (`plan`),
//! turns each `cell + 1` into `slot + 1` `2 * D` places ahead of the update (prefetching the
//! cell's `ranks` element `D` places before that and the slot's `heads` element after), and
//! applies NCBI's update to `heads[slot]`. So the chains are the fill's, and `head` returns the
//! fill's `hashtable`.
//!
//! With LOSAT_X_MBHRANKTRACE=1 each build prints `[X_HRANK] words= present= cells= heads_kib=
//! ranks_kib= pass1_ms= pass2_ms= pass3_ms=`; a chain table built without the probe prints its
//! present cells as `[X_HRANK] chains words= present= cells=`.

use super::x_mb_lookup_par::Words;
use super::{PvArrayType, PV_ARRAY_BTS};

const K_COMPRESSION_FACTOR: usize = 2048;
/// The places between the prefetch of a `ranks` element and its read, and between the prefetch
/// of a `heads` element and its update.
const D: usize = 32;

pub(super) fn on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("LOSAT_X_MBHRANK").is_ok_and(|v| v == "1"))
}

pub(super) fn trace() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("LOSAT_X_MBHRANKTRACE").is_ok_and(|v| v == "1"))
}

/// The probe's scope: 4^12 cells, newest first, no database word counts, no BLEMIR_DEBUG, a PV
/// word of at least 32 cells, and queries whose offsets ascend without overlap (index order is
/// then the fill's order).
pub(super) fn plan(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    pv_array_bts: usize,
) -> bool {
    if !on() || cells != 1 << 24 || words.ascending_cells || db_word_counts || debug {
        return false;
    }
    if pv_array_bts < PV_ARRAY_BTS {
        return false;
    }
    let mut end = 0usize;
    for (seq, &offset) in words.queries_blastna.iter().zip(words.query_offsets) {
        let start = offset.max(0) as usize;
        if start < end {
            return false;
        }
        end = start + seq.len();
    }
    true
}

pub(super) struct HeadRank {
    ranks: Vec<u64>,
    heads: Vec<u32>,
}

impl HeadRank {
    /// `hashtable[cell]` of the fill.
    #[inline(always)]
    pub(super) fn head(&self, cell: usize) -> u32 {
        let Some(&r) = self.ranks.get(cell >> 5) else {
            return 0;
        };
        let bits = r as u32;
        let bit = 1u32 << (cell & 31);
        if bits & bit == 0 {
            return 0;
        }
        let slot = (r >> 32) as usize + (bits & (bit - 1)).count_ones() as usize;
        self.heads[slot]
    }

    pub(super) fn present(&self) -> usize {
        self.heads.len()
    }
}

/// The table's `ranks` and `heads`, `next_pos`, presence vector and helper counters (zeroed, of
/// the fill's lengths).
pub(super) fn serial(
    words: &Words<'_>,
    cells: usize,
    next_pos: &mut [u32],
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
) -> HeadRank {
    let t0 = std::time::Instant::now();
    let mut ranks = vec![0u64; cells >> 5];
    let mut n_words = 0usize;
    super::x_bn_csr::each_word(words, |cell, q_off_1| {
        next_pos[q_off_1 as usize] = cell as u32 + 1;
        ranks[cell >> 5] |= 1u64 << (cell & 31);
        let h = &mut helper_array[cell / K_COMPRESSION_FACTOR];
        *h = h.wrapping_add(1);
        n_words += 1;
    });
    let t1 = std::time::Instant::now();
    let present = sums(&mut ranks, pv_array, pv_array_bts, helper_array);
    let t2 = std::time::Instant::now();
    let mut heads = vec![0u32; present];
    apply(&ranks, &mut heads, next_pos);
    let t3 = std::time::Instant::now();
    if trace() {
        eprintln!(
            "[X_HRANK] words={} present={} cells={} heads_kib={} ranks_kib={} pass1_ms={:.2} pass2_ms={:.2} pass3_ms={:.2}",
            n_words,
            present,
            cells,
            present * 4 / 1024,
            ranks.len() * 8 / 1024,
            (t1 - t0).as_secs_f64() * 1e3,
            (t2 - t1).as_secs_f64() * 1e3,
            (t3 - t2).as_secs_f64() * 1e3,
        );
    }
    HeadRank { ranks, heads }
}

/// Pass 2: running counts into the high halves, the presence vector and the helper counters.
fn sums(
    ranks: &mut [u64],
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
) -> usize {
    let pv_shift = pv_array_bts - PV_ARRAY_BTS;
    let mut acc = 0u64;
    for (w, r) in ranks.iter_mut().enumerate() {
        let bits = *r as u32;
        *r = (acc << 32) | u64::from(bits);
        if bits != 0 {
            let n = bits.count_ones();
            acc += u64::from(n);
            if let Some(pv) = pv_array.get_mut(w >> pv_shift) {
                *pv |= bits;
            }
            let h = &mut helper_array[(w << 5) / K_COMPRESSION_FACTOR];
            *h = h.wrapping_sub(n);
        }
    }
    acc as usize
}

#[inline(always)]
fn slot(ranks: &[u64], cell: usize) -> usize {
    let r = ranks[cell >> 5];
    let bits = r as u32;
    (r >> 32) as usize + (bits & ((1u32 << (cell & 31)) - 1)).count_ones() as usize
}

/// Pass 3 (module comment).
fn apply(ranks: &[u64], heads: &mut [u32], next_pos: &mut [u32]) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("popcnt") {
        // SAFETY: popcnt is available.
        return unsafe { apply_popcnt(ranks, heads, next_pos) };
    }
    apply_any(ranks, heads, next_pos)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "popcnt")]
unsafe fn apply_popcnt(ranks: &[u64], heads: &mut [u32], next_pos: &mut [u32]) {
    apply_any(ranks, heads, next_pos)
}

#[inline(always)]
fn apply_any(ranks: &[u64], heads: &mut [u32], next_pos: &mut [u32]) {
    let n = next_pos.len();
    for q in 1..(1 + 2 * D).min(n) {
        let c1 = next_pos[q];
        if c1 != 0 {
            next_pos[q] = slot(ranks, (c1 - 1) as usize) as u32 + 1;
        }
    }
    for q in 1..n {
        if let Some(&c1) = next_pos.get(q + 3 * D) {
            if c1 != 0 {
                prefetch(ranks.as_ptr().wrapping_add((c1 - 1) as usize >> 5).cast());
            }
        }
        let qb = q + 2 * D;
        if qb < n {
            let c1 = next_pos[qb];
            if c1 != 0 {
                let s = slot(ranks, (c1 - 1) as usize);
                next_pos[qb] = s as u32 + 1;
                prefetch(heads.as_ptr().wrapping_add(s).cast());
            }
        }
        let s1 = next_pos[q];
        if s1 != 0 {
            let s = (s1 - 1) as usize;
            next_pos[q] = heads[s];
            heads[s] = q as u32;
        }
    }
}

#[inline(always)]
fn prefetch(_p: *const i8) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not fault.
    unsafe {
        core::arch::x86_64::_mm_prefetch(_p, core::arch::x86_64::_MM_HINT_T0)
    };
}

/// The non-empty cells of a chain table (LOSAT_X_MBHRANKTRACE only).
pub(super) fn trace_chains(hashtable: &[u32], next_pos: &[u32]) {
    let present = hashtable.iter().filter(|&&h| h != 0).count();
    let mut words = 0usize;
    for &h in hashtable {
        let mut q = h;
        while q != 0 {
            words += 1;
            q = next_pos[q as usize];
        }
    }
    eprintln!(
        "[X_HRANK] chains words={} present={} cells={}",
        words,
        present,
        hashtable.len()
    );
}
