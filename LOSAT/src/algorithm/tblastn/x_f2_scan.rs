//! PROBE (line F s3, p3; not for merging): arms of tblastn's lookup scan that give the same
//! offset pairs and `s_range` as `x_scan_one_range_bf`, to cost the scan's floor on a real input.
//!
//! `LOSAT_X_F2SCAN` picks the arm (read once; default 0, which is `x_scan_one_range_bf` as it is):
//!   1  V1: the presence test reads a one-byte-per-word table (`pvb`, 32 KB, built from `pv` with
//!      the plan) instead of the 4 KB bit vector; the rolling index is kept as in the base loop.
//!   2  V2: V1's table, and the indexes of a block are made 16 at a time with AVX2
//!      (`ix_k = l[k] << 10 | l[k+1] << 5 | l[k+2]`, three shifted loads, no dependency from word
//!      to word); the scalar loop reads the index and the byte, and writes the word's position
//!      without a branch (the list's length advances by the byte).
//!   3  the base loop (`x_candidates_w3c5` as it is) in this file's block structure, with the
//!      candidate loop in a function of its own: the control for the structure, and the way to
//!      split the base's cycles into candidate loop and copy loop.
//! `LOSAT_X_F2CHECK` (set) counts words, candidates and hits of every call and sums a hash of every
//! call's offset pairs and `s_range[1]`; `main` prints them as `[X_STATS] F2SCAN ...`.

use super::search_seed::XScanPlan;
use super::x_scan_bf::{x_scan_one_range_bf, X_BLOCK, X_SCAN_BF_SLACK};
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
use crate::algorithm::tblastx::lookup::{BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS};
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicU64, Ordering};

const MASK: usize = (1 << 15) - 1;

pub(super) fn x_f2_arm() -> u8 {
    static ARM: std::sync::OnceLock<u8> = std::sync::OnceLock::new();
    *ARM.get_or_init(|| {
        std::env::var("LOSAT_X_F2SCAN")
            .ok()
            .and_then(|v| v.parse::<u8>().ok())
            .filter(|v| *v <= 3)
            .unwrap_or(0)
    })
}

pub(super) fn x_f2_check() -> bool {
    static CHECK: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *CHECK.get_or_init(|| std::env::var_os("LOSAT_X_F2CHECK").is_some())
}

/// Whether `x_scan_unit` leaves the base path.
pub(super) fn x_f2_active() -> bool {
    static ACTIVE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ACTIVE.get_or_init(|| x_f2_arm() != 0 || x_f2_check())
}

fn w3c5(lookup: &BlastAaLookupTable) -> bool {
    lookup.word_length == 3 && lookup.charsize == 5 && lookup.mask as usize == MASK
}

/// The one-byte-per-word presence table (arms 1 and 2); empty otherwise.
pub(super) fn x_f2_build_pvb(lookup: &BlastAaLookupTable) -> Vec<u8> {
    if x_f2_check() {
        static ONCE: std::sync::Once = std::sync::Once::new();
        ONCE.call_once(|| {
            eprintln!(
                "[X_STATS] F2SCAN lookup word_length={} charsize={} mask={:#x} pv_len={} backbone_len={} arm={}",
                lookup.word_length,
                lookup.charsize,
                lookup.mask,
                lookup.pv.len(),
                lookup.backbone.len(),
                x_f2_arm()
            )
        });
    }
    if x_f2_arm() == 0 || !w3c5(lookup) || lookup.pv.len() < ((MASK + 1) >> PV_ARRAY_BTS) {
        return Vec::new();
    }
    (0..=MASK)
        .map(|ix| ((lookup.pv[ix >> PV_ARRAY_BTS] >> (ix & 31)) & 1) as u8)
        .collect()
}

static CALLS: AtomicU64 = AtomicU64::new(0);
static UNITS: AtomicU64 = AtomicU64::new(0);
static RANGES: AtomicU64 = AtomicU64::new(0);
static WORDS_USED: AtomicU64 = AtomicU64::new(0);
static WORDS_LOOKED: AtomicU64 = AtomicU64::new(0);
static CANDS: AtomicU64 = AtomicU64::new(0);
static HITS: AtomicU64 = AtomicU64::new(0);
static STOPS: AtomicU64 = AtomicU64::new(0);
static HASH: AtomicU64 = AtomicU64::new(0);

pub(super) fn print_summary() {
    if x_f2_check() {
        eprintln!(
            "[X_STATS] F2SCAN arm={} units={} ranges={} calls={} stops={} words_used={} words_looked={} cands={} hits={} hash={:016x}",
            x_f2_arm(),
            UNITS.load(Ordering::Relaxed),
            RANGES.load(Ordering::Relaxed),
            CALLS.load(Ordering::Relaxed),
            STOPS.load(Ordering::Relaxed),
            WORDS_USED.load(Ordering::Relaxed),
            WORDS_LOOKED.load(Ordering::Relaxed),
            CANDS.load(Ordering::Relaxed),
            HITS.load(Ordering::Relaxed),
            HASH.load(Ordering::Relaxed),
        );
    }
}

/// `x_scan_unit` for the probe arms; same ranges, same calls, same `on_pairs` sequence.
pub(super) fn x_scan_unit_f2(
    plan: &XScanPlan,
    unit: usize,
    pairs: &mut Vec<BlastOffsetPair>,
    mut on_pairs: impl FnMut(&[BlastOffsetPair]),
) {
    let arm = x_f2_arm();
    let check = x_f2_check();
    let (frame_index, chunk, ranges) = &plan.subject.units[unit];
    let frame = &plan.subject.frames[*frame_index];
    let pairs_len = plan.pair_capacity + X_SCAN_BF_SLACK;
    if pairs.len() < pairs_len {
        pairs.resize(pairs_len, BlastOffsetPair::default());
    }
    let pair_capacity = plan.pair_capacity as i32;
    let (mut calls, mut stops, mut used, mut looked, mut cands, mut hits_sum, mut hash) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    for (range_index, &(left, right)) in ranges.iter().enumerate() {
        let mut end = right - plan.lookup.word_length as i32;
        if range_index == 0 && end < left {
            end = left;
        }
        let mut scan_range = [0, left, end];
        while scan_range[1] <= scan_range[2] {
            let first = scan_range[1];
            let hits = x_scan_range_f2(
                arm,
                &plan.lookup,
                &plan.pvb,
                &frame.aa_seq[1 + chunk.offset..],
                pairs,
                pair_capacity,
                &mut scan_range,
            );
            if check {
                calls += 1;
                let next = scan_range[1];
                let last = scan_range[2];
                if next > last {
                    used += (last - first + 1) as u64;
                    looked += (last - first + 1) as u64;
                } else {
                    stops += 1;
                    used += (next - first) as u64;
                    let block = (next - first) as usize / X_BLOCK;
                    let block_end = (first as usize + X_BLOCK * (block + 1) - 1).min(last as usize);
                    looked += (block_end - first as usize + 1) as u64;
                }
                let mut h: u64 = 0xcbf2_9ce4_8422_2325 ^ hits as u64;
                let mut previous = u32::MAX;
                for pair in &pairs[..hits as usize] {
                    if pair.s_off != previous {
                        cands += 1;
                        previous = pair.s_off;
                    }
                    h = (h ^ (((pair.q_off as u64) << 32) | pair.s_off as u64))
                        .wrapping_mul(0x0000_0100_0000_01b3)
                        .rotate_left(17);
                }
                h = (h ^ next as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
                hash = hash.wrapping_add(h);
                hits_sum += hits as u64;
            }
            on_pairs(&pairs[..hits as usize]);
        }
    }
    if check {
        UNITS.fetch_add(1, Ordering::Relaxed);
        RANGES.fetch_add(ranges.len() as u64, Ordering::Relaxed);
        CALLS.fetch_add(calls, Ordering::Relaxed);
        STOPS.fetch_add(stops, Ordering::Relaxed);
        WORDS_USED.fetch_add(used, Ordering::Relaxed);
        WORDS_LOOKED.fetch_add(looked, Ordering::Relaxed);
        CANDS.fetch_add(cands, Ordering::Relaxed);
        HITS.fetch_add(hits_sum, Ordering::Relaxed);
        HASH.fetch_add(hash, Ordering::Relaxed);
    }
}

fn x_scan_range_f2(
    arm: u8,
    lookup: &BlastAaLookupTable,
    pvb: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    if arm == 0 || !w3c5(lookup) || (arm != 3 && pvb.len() != MASK + 1) {
        return x_scan_one_range_bf(lookup, subject, offset_pairs, array_size, s_range);
    }
    #[cfg(target_arch = "x86_64")]
    let avx2 = std::is_x86_feature_detected!("avx2");
    #[cfg(not(target_arch = "x86_64"))]
    let avx2 = false;
    match arm {
        1 => scan::<1>(lookup, pvb, subject, offset_pairs, array_size, s_range),
        2 if avx2 => scan::<2>(lookup, pvb, subject, offset_pairs, array_size, s_range),
        2 => scan::<1>(lookup, pvb, subject, offset_pairs, array_size, s_range),
        _ => scan::<3>(lookup, pvb, subject, offset_pairs, array_size, s_range),
    }
}

/// Copies the hits of one candidate word as `x_scan_one_range_bf` does; `None` when the offset
/// array has no room (the scan stops at this word).
#[inline(always)]
unsafe fn copy_cell(
    backbone: *const crate::algorithm::tblastx::lookup::BackboneCell,
    overflow: *const i32,
    dest: *mut BlastOffsetPair,
    ix: usize,
    s_off: u32,
    totalhits: i32,
    array_size: i32,
) -> Option<i32> {
    let cell = &*backbone.add(ix);
    let numhits = cell.num_used;
    if numhits > array_size - totalhits {
        return None;
    }
    let d = dest.add(totalhits as usize);
    if numhits as usize <= AA_HITS_PER_CELL {
        for i in 0..AA_HITS_PER_CELL {
            *d.add(i) = BlastOffsetPair {
                q_off: cell.entries[i] as u32,
                s_off,
            };
        }
    } else {
        let src = overflow.add(cell.entries[0] as usize);
        for i in 0..numhits as usize {
            *d.add(i) = BlastOffsetPair {
                q_off: *src.add(i) as u32,
                s_off,
            };
        }
    }
    Some(totalhits + numhits)
}

fn scan<const ARM: u8>(
    lookup: &BlastAaLookupTable,
    pvb: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    if s_range[1] > s_range[2] {
        return 0;
    }
    assert!(
        array_size >= 0 && offset_pairs.len() >= array_size as usize + X_SCAN_BF_SLACK,
        "offset array of {} entries for array_size {array_size}",
        offset_pairs.len()
    );
    let pv = lookup.pv.as_ptr();
    let pvb = pvb.as_ptr();
    let backbone = lookup.backbone.as_ptr();
    let overflow = lookup.overflow.as_ptr();
    let dest = offset_pairs.as_mut_ptr();
    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;

    let mut index: usize = 0;
    for i in 0..2 {
        let ch = unsafe { *subject.get_unchecked(s_first + i) } as usize;
        index = (index << 5) | ch;
    }
    let mut cand = [MaybeUninit::<u64>::uninit(); X_BLOCK];
    let mut cand16 = [MaybeUninit::<u16>::uninit(); X_BLOCK];
    let mut ixbuf = [0u16; X_BLOCK + 16];
    let mut totalhits: i32 = 0;
    let mut block_first = s_first;
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        let block = (block_first, block_last);
        if ARM == 2 {
            #[cfg(target_arch = "x86_64")]
            {
                // SAFETY: AVX2 was detected by the caller.
                let n = unsafe { v2_candidates(subject, pvb, block, &mut ixbuf, &mut cand16) };
                for c in &cand16[..n] {
                    // SAFETY: the first `n` entries were written.
                    let k = unsafe { c.assume_init() } as usize;
                    let s_off = (block_first + k) as u32;
                    let ix = ixbuf[k] as usize;
                    // SAFETY: as in `x_scan_one_range_bf`.
                    match unsafe {
                        copy_cell(backbone, overflow, dest, ix, s_off, totalhits, array_size)
                    } {
                        Some(total) => totalhits = total,
                        None => {
                            s_range[1] = s_off as i32;
                            return totalhits;
                        }
                    }
                }
            }
        } else {
            let n = if ARM == 1 {
                v1_candidates(subject, pvb, block, &mut index, &mut cand)
            } else {
                base_candidates(subject, pv, block, &mut index, &mut cand)
            };
            for c in &cand[..n] {
                // SAFETY: the first `n` candidates were written.
                let c = unsafe { c.assume_init() };
                let s_off = (c >> 32) as u32;
                let ix = (c & 0xffff_ffff) as usize;
                // SAFETY: as in `x_scan_one_range_bf`.
                match unsafe { copy_cell(backbone, overflow, dest, ix, s_off, totalhits, array_size) }
                {
                    Some(total) => totalhits = total,
                    None => {
                        s_range[1] = s_off as i32;
                        return totalhits;
                    }
                }
            }
        }
        block_first = block_last + 1;
    }
    s_range[1] = (s_last + 1) as i32;
    totalhits
}

/// Arm 3: `x_candidates_w3c5` of `x_scan_bf.rs`, unchanged, in a function of its own.
#[inline(never)]
fn base_candidates(
    subject: &[u8],
    pv: *const u32,
    block: (usize, usize),
    index: &mut usize,
    cand: &mut [MaybeUninit<u64>; X_BLOCK],
) -> usize {
    let words = block.1 - block.0 + 1;
    let letters = unsafe { subject.as_ptr().add(block.0 + 2) };
    let out = cand.as_mut_ptr().cast::<u64>();
    let mut ix = *index;
    let mut n = 0usize;
    let mut position = (block.0 as u64) << 32;
    for k in 0..words {
        let new_char = unsafe { *letters.add(k) } as usize;
        ix = ((ix << 5) | new_char) & MASK;
        let bit = (unsafe { *pv.add(ix >> PV_ARRAY_BTS) } >> (ix & 31)) & 1;
        unsafe { out.add(n).write(position | ix as u64) };
        n += bit as usize;
        position += 1 << 32;
    }
    *index = ix;
    n
}

/// Arm 1: the same loop with the presence test on the byte table.
#[inline(never)]
fn v1_candidates(
    subject: &[u8],
    pvb: *const u8,
    block: (usize, usize),
    index: &mut usize,
    cand: &mut [MaybeUninit<u64>; X_BLOCK],
) -> usize {
    let words = block.1 - block.0 + 1;
    let letters = unsafe { subject.as_ptr().add(block.0 + 2) };
    let out = cand.as_mut_ptr().cast::<u64>();
    let mut ix = *index;
    let mut n = 0usize;
    let mut position = (block.0 as u64) << 32;
    for k in 0..words {
        let new_char = unsafe { *letters.add(k) } as usize;
        ix = ((ix << 5) | new_char) & MASK;
        let bit = unsafe { *pvb.add(ix) };
        unsafe { out.add(n).write(position | ix as u64) };
        n += bit as usize;
        position += 1 << 32;
    }
    *index = ix;
    n
}

/// Arm 2: indexes 16 at a time (AVX2), then a scalar loop on the byte table that lists the
/// positions (relative to the block) of the words whose byte is 1.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn v2_candidates(
    subject: &[u8],
    pvb: *const u8,
    block: (usize, usize),
    ixbuf: &mut [u16; X_BLOCK + 16],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    use std::arch::x86_64::*;
    let words = block.1 - block.0 + 1;
    let base = subject.as_ptr().add(block.0);
    let available = subject.len() - block.0;
    let ixp = ixbuf.as_mut_ptr();
    let mut k = 0usize;
    let m = _mm256_set1_epi16(MASK as i16);
    // Word k needs the letters k, k + 1, k + 2; a group of 16 reads bytes k .. k + 17.
    while k + 16 <= words && k + 18 <= available {
        let a = _mm256_cvtepu8_epi16(_mm_loadu_si128(base.add(k) as *const __m128i));
        let b = _mm256_cvtepu8_epi16(_mm_loadu_si128(base.add(k + 1) as *const __m128i));
        let c = _mm256_cvtepu8_epi16(_mm_loadu_si128(base.add(k + 2) as *const __m128i));
        let v = _mm256_or_si256(
            _mm256_or_si256(_mm256_slli_epi16(a, 10), _mm256_slli_epi16(b, 5)),
            c,
        );
        _mm256_storeu_si256(ixp.add(k) as *mut __m256i, _mm256_and_si256(v, m));
        k += 16;
    }
    while k < words {
        let v = ((*base.add(k) as u16) << 10)
            | ((*base.add(k + 1) as u16) << 5)
            | (*base.add(k + 2) as u16);
        *ixp.add(k) = v & MASK as u16;
        k += 1;
    }
    let out = cand.as_mut_ptr().cast::<u16>();
    let mut n = 0usize;
    for k in 0..words {
        let ix = *ixp.add(k) as usize;
        let bit = *pvb.add(ix) as usize;
        out.add(n).write(k as u16);
        n += bit;
    }
    n
}
