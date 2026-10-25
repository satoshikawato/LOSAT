//! PROBE (LOSAT_X_MB24=1|2, branch exp/mb-head-rank-probe; v4-A10 s20; not a PR): the chains
//! (`next_pos`, =1) or the chains and the head table (`hashtable`, =2) of a contiguous megablast
//! table of 4^12 cells (`build_mb_lookup`) with 3-byte elements instead of 4-byte ones. An element
//! holds a 1-based query offset, so a table whose `next_pos` has fewer than 2^24 elements loses
//! nothing. Element `i` is the low 24 bits of the little-endian u32 at byte `3 * i` (each array has
//! one byte more than its elements, so that read stays inside it).
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
//! The fill (`serial`, one thread) is `x_mb_fill::run_words`'s: the words read run by run, the
//! cell of each prefetched, the update applied `X_DELAY` words later. The chain element of a word
//! is written as a whole u32 (its top byte, the next element's first byte, is still 0: the words
//! come in increasing index order, `plan`); a head is written by merging its 3 bytes into the u32.

use super::x_mb_lookup_par::Words;
use super::{PvArrayType, BLAST2NA_MASK, PV_ARRAY_MASK};

const X_DELAY: usize = 64;

pub(super) fn mode() -> u8 {
    static MODE: std::sync::OnceLock<u8> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| match std::env::var("LOSAT_X_MB24").as_deref() {
        Ok("1") => 1,
        Ok("2") => 2,
        _ => 0,
    })
}

/// The scope of LOSAT_X_MBNPKEEP for the 4-byte chains: `plan`'s without a mode.
pub(super) fn np_scope(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    next_pos_len: usize,
    helper_len: usize,
) -> bool {
    np_keep()
        && scope(
            words,
            cells,
            db_word_counts,
            debug,
            next_pos_len,
            helper_len,
        )
}

/// The probe's scope: 4^12 cells, newest first, no database word counts, no BLEMIR_DEBUG,
/// `next_pos` of fewer than 2^24 elements, and queries whose offsets ascend without overlap.
pub(super) fn plan(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    next_pos_len: usize,
    helper_len: usize,
) -> u8 {
    let m = mode();
    if m == 0
        || !scope(
            words,
            cells,
            db_word_counts,
            debug,
            next_pos_len,
            helper_len,
        )
    {
        return 0;
    }
    m
}

fn scope(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    next_pos_len: usize,
    helper_len: usize,
) -> bool {
    if cells != 1 << 24
        || words.kmer_mask != (1u64 << 24) - 1
        || words.lut_word_length != 12
        || words.ascending_cells
        || db_word_counts
        || debug
        || next_pos_len >= 1 << 24
        || (words.kmer_mask as usize >> 11) >= helper_len
    {
        return false;
    }
    let mut end = 0usize;
    for (seq, &offset) in words.queries_blastna.iter().zip(words.query_offsets) {
        let start = offset.max(0) as usize;
        if start < end || start + seq.len() >= next_pos_len {
            return false;
        }
        end = start + seq.len();
    }
    true
}

/// LOSAT_X_MBNPKEEP=1: the chains' buffer of a build in this module's scope goes back to a spare
/// when the table is dropped and the next build takes it without zeroing it: a chain is read only
/// at the indexes of words (a head, then `next_pos` of each word), which the fill writes, so the
/// other elements may hold anything.
pub(super) fn np_keep() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("LOSAT_X_MBNPKEEP").is_ok_and(|v| v == "1"))
}

static SPARE_NP: std::sync::Mutex<Vec<u32>> = std::sync::Mutex::new(Vec::new());
static SPARE_CHAIN: std::sync::Mutex<Vec<u8>> = std::sync::Mutex::new(Vec::new());

/// A `next_pos` of `len` elements: the spare's first `len` (not zeroed) when it is long enough,
/// otherwise a zeroed one.
pub(super) fn take_np(len: usize) -> Vec<u32> {
    let mut v = SPARE_NP
        .lock()
        .map(|mut s| std::mem::take(&mut *s))
        .unwrap_or_default();
    if v.len() >= len {
        v.truncate(len);
        v
    } else {
        drop(v);
        vec![0u32; len]
    }
}

pub(super) fn give_np(v: Vec<u32>) {
    if let Ok(mut s) = SPARE_NP.lock() {
        *s = v;
    }
}

fn take_chain(len: usize) -> Vec<u8> {
    if np_keep() {
        let mut v = SPARE_CHAIN
            .lock()
            .map(|mut s| std::mem::take(&mut *s))
            .unwrap_or_default();
        if v.len() >= len {
            v.truncate(len);
            return v;
        }
    }
    vec![0u8; len]
}

impl Drop for Narrow {
    fn drop(&mut self) {
        if np_keep() {
            if let Ok(mut s) = SPARE_CHAIN.lock() {
                *s = std::mem::take(&mut self.chain);
            }
        }
    }
}

pub(super) struct Narrow {
    /// 3-byte heads (=2), empty with =1 (the table's `hashtable` holds them).
    heads: Vec<u8>,
    chain: Vec<u8>,
}

#[inline(always)]
fn read3(v: &[u8], i: usize) -> u32 {
    let b = &v[3 * i..3 * i + 4];
    u32::from_le_bytes([b[0], b[1], b[2], b[3]]) & 0x00FF_FFFF
}

impl Narrow {
    #[inline(always)]
    pub(super) fn head(&self, hashtable: &[u32], cell: usize) -> u32 {
        if self.heads.is_empty() {
            hashtable[cell]
        } else {
            read3(&self.heads, cell)
        }
    }

    #[inline(always)]
    pub(super) fn next(&self, q_off: u32) -> u32 {
        read3(&self.chain, q_off as usize)
    }
}

/// The table's 3-byte arrays (and with =1 `hashtable`), presence vector and helper counters
/// (zeroed, of the fill's lengths).
pub(super) fn serial(
    m: u8,
    words: &Words<'_>,
    cells: usize,
    hashtable: &mut [u32],
    next_pos_len: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
) -> Narrow {
    let mut heads = if m == 2 {
        vec![0u8; 3 * cells + 1]
    } else {
        Vec::new()
    };
    let mut chain = take_chain(3 * next_pos_len + 1);
    assert!(m == 2 || hashtable.len() == cells);
    assert!((cells - 1) >> pv_array_bts < pv_array.len());
    let raw = Raw {
        heads: heads.as_mut_ptr(),
        hashtable: hashtable.as_mut_ptr(),
        chain: chain.as_mut_ptr(),
        pv_array: pv_array.as_mut_ptr(),
        pv_array_bts,
        helper_array: helper_array.as_mut_ptr(),
    };
    // SAFETY: `plan` held: every cell is below 2^24 = `cells`, every index below `next_pos_len`.
    unsafe {
        if m == 2 {
            fill::<true>(words, raw)
        } else {
            fill::<false>(words, raw)
        }
    };
    Narrow { heads, chain }
}

#[derive(Clone, Copy)]
struct Raw {
    heads: *mut u8,
    hashtable: *mut u32,
    chain: *mut u8,
    pv_array: *mut PvArrayType,
    pv_array_bts: usize,
    helper_array: *mut u32,
}

impl Raw {
    #[inline(always)]
    unsafe fn cell_ptr<const H24: bool>(self, bucket: usize) -> *const i8 {
        if H24 {
            unsafe { self.heads.add(3 * bucket).cast() }
        } else {
            unsafe { self.hashtable.add(bucket).cast() }
        }
    }

    /// SAFETY: as `serial`'s.
    #[inline(always)]
    unsafe fn apply<const H24: bool>(self, bucket: usize, q_off_1: u32) {
        unsafe {
            let old = if H24 {
                let p = self.heads.add(3 * bucket);
                let w = (p as *const u32).read_unaligned();
                (p as *mut u32).write_unaligned((w & 0xFF00_0000) | q_off_1);
                w & 0x00FF_FFFF
            } else {
                let p = self.hashtable.add(bucket);
                let old = *p;
                *p = q_off_1;
                old
            };
            if old == 0 {
                *self.pv_array.add(bucket >> self.pv_array_bts) |= 1u32 << (bucket & PV_ARRAY_MASK);
            } else {
                let h = self.helper_array.add(bucket >> 11);
                *h = (*h).wrapping_add(1);
            }
            (self.chain.add(3 * q_off_1 as usize) as *mut u32).write_unaligned(old);
        }
    }
}

/// SAFETY: as `serial`'s.
unsafe fn fill<const H24: bool>(words: &Words<'_>, raw: Raw) {
    let mut ring = [(0u32, 0u32); X_DELAY];
    let mut n = 0usize;
    let lut_word_length = words.lut_word_length;
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
        for (range_start, range_end) in super::build_unmasked_ranges(seq.len(), masks) {
            if words.word_length > range_end.saturating_sub(range_start) {
                continue;
            }
            let s = &seq[range_start..range_end];
            let index_base = (query_offset + range_start + 2).wrapping_sub(lut_word_length);
            let mut i = 0usize;
            'run: while i < s.len() {
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
                        continue;
                    }
                    code = (code << 2) | base as u64;
                    got += 1;
                }
                // SAFETY: as this function's.
                i = unsafe {
                    run_words::<H24>(
                        s,
                        i,
                        code,
                        words.kmer_mask,
                        index_base,
                        &mut ring,
                        &mut n,
                        raw,
                    )
                };
                if i < s.len() {
                    i += 1;
                }
            }
        }
    }
    for k in n - n.min(X_DELAY)..n {
        let (b, q) = ring[k % X_DELAY];
        // SAFETY: as above.
        unsafe { raw.apply::<H24>(b as usize, q) };
    }
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn run_words<const H24: bool>(
    s: &[u8],
    mut i: usize,
    mut code: u64,
    kmer_mask: u64,
    index_base: usize,
    ring_ref: &mut [(u32, u32); X_DELAY],
    n_ref: &mut usize,
    raw: Raw,
) -> usize {
    let mut ring = *ring_ref;
    let mut n = *n_ref;
    while n < X_DELAY {
        let Some(&base) = s.get(i) else {
            break;
        };
        if (base & BLAST2NA_MASK) != 0 {
            break;
        }
        code = ((code << 2) | base as u64) & kmer_mask;
        prefetch(unsafe { raw.cell_ptr::<H24>(code as usize) });
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
            prefetch(unsafe { raw.cell_ptr::<H24>(code as usize) });
            let slot = n % X_DELAY;
            let (b, q) = ring[slot];
            // SAFETY: `(b, q)` is a word of the fill.
            unsafe { raw.apply::<H24>(b as usize, q) };
            ring[slot] = (code as u32, index_base.wrapping_add(i) as u32);
            n += 1;
            i += 1;
        }
    }
    *ring_ref = ring;
    *n_ref = n;
    i
}

#[inline(always)]
fn prefetch(_p: *const i8) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not fault.
    unsafe {
        core::arch::x86_64::_mm_prefetch(_p, core::arch::x86_64::_MM_HINT_T0)
    };
}
