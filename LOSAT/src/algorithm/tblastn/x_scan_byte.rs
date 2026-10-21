//! EXPERIMENT (LOSAT_X_TBNSCANBYTE): tblastn's lookup scan with the table indexes of a block of
//! words made at once and their presence read from one byte per index, and the hits of a block
//! copied with one capacity test.
//!
//! NCBI's scan (c++/src/algo/blast/core/blast_aascan.c:75-127, `s_BlastAaScanSubject`, kept as
//! `tblastx/blast_aascan.rs` `s_blast_aa_scan_subject_one_range`) rolls the index of the next word
//! from the one before (`ComputeTableIndexIncremental`), tests its presence-vector bit and copies
//! the word's hits when the offset array has room. LOSAT_X_TBNSCANBF (`x_scan_bf.rs`) takes the
//! branch on the bit away, in blocks of `X_BLOCK` words; on a genome searched with one protein
//! (HsG_p1 x hsap_dna `-db`, n1) the scan was still 55-58 % of the cycles, 5.38 cycles per word,
//! with the rolling index a chain from word to word (line F probe P3, v4-F/REPORT.md).
//!
//! For a table of word size 3 and 5-bit letters (tblastn's, `fits`) the index of the word at `s`
//! is `(subject[s] << 10 | subject[s + 1] << 5 | subject[s + 2]) & 0x7fff`: the rolled index is
//! `((index << 5) | letter) & mask`, and the letters before `s` are shifted past bit 14 and masked
//! away, so both give the same number for any bytes. Here each block's indexes are made from the
//! letters directly, 16 at a time with AVX2 (or by a plain loop), with no chain between words. Its
//! candidates are listed from `presence`, the presence vector with one byte per index (1 where
//! `pv` has the bit, 0 elsewhere, 32 KB for the 15-bit index space), so the test is one load
//! instead of a load, a shift and a mask. Indexes past the backbone read 0 and are never
//! candidates.
//!
//! The copy then takes the block's candidates in position order, as NCBI's loop does. While three
//! hits for every candidate left in the block fit in the offset array (one test for the rest of the
//! block), the hits of a cell of up to `AA_HITS_PER_CELL` are written without the capacity test, as
//! four pairs (with AVX2 one 32-byte store: the cell broadcast, one `vpermd`, one `vpblendd` with
//! the subject offset) of which the first `num_used` count. A cell with more hits, and every
//! candidate of a block that could fill the array, takes NCBI's test: when its hits do not fit,
//! `s_range[1]` is that word's position and the scan returns, so the next call resumes there. The
//! offset pairs and the `s_range` values of every call are NCBI's (the words after the stop were
//! looked at but not used; the next call looks at them again, as with LOSAT_X_TBNSCANBF).
//!
//! The offset array must have `X_SCAN_BYTE_SLACK` entries after `array_size`; the entries past the
//! returned count are not part of the result. `presence` depends on the lookup table only; each
//! thread builds it once per table (`presence_bytes`).
//!
//! Probe and source of the form: line F s3 P3 V2 (candidates) and s5 P3b arm 7 (copy),
//! `exp/f2-tbn-scan` ff918e61: HsG_p1 x hsap_dna `-db` n1 whole-run cycles x0.703 of dev, with
//! the offset pairs and `s_range` of every call and the output equal.
//!
//! The switch is read once. It is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_TBNSCANBYTE=0`, or `LOSAT_X_DEFAULT=0` and it is not set. On x86-64 with AVX2 the
//! scan uses the AVX2 forms of both stages; elsewhere, and with `LOSAT_X_TBNSCANBYTE=2`, the plain
//! loops (same values). A table of another shape is scanned as before (LOSAT_X_TBNSCANBF or NCBI's
//! loop). With `LOSAT_X_TBNSCANBYTESHADOW` set (any value; it turns the path on) every scan call
//! also runs NCBI's loop on a second offset array and the count, the pairs and `s_range` must be
//! equal; `main` prints `[X_SHADOW] TBNSCANBYTE calls_compared= pairs_compared=` at exit when the
//! unit path scanned any range this way. With `LOSAT_X_STATS` set, `main` prints
//! `[X_STATS] TBNSCANBYTE units= calls= avx2=` (the units scanned this way, their scan calls, and
//! whether the AVX2 forms ran).

use crate::algorithm::tblastx::blast_aascan::{s_blast_aa_scan_subject_one_range, BlastOffsetPair};
use crate::algorithm::tblastx::lookup::{
    BackboneCell, BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): subject words whose indexes are made before their hits are
/// copied.
const X_BLOCK: usize = 256;

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the index mask of a table of word size 3 and 5-bit letters.
const X_MASK: usize = (1 << 15) - 1;

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): offset-array entries needed after `array_size` (a cell of up
/// to `AA_HITS_PER_CELL` hits is written as `AA_HITS_PER_CELL + 1` pairs).
pub(super) const X_SCAN_BYTE_SLACK: usize = AA_HITS_PER_CELL + 1;

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): how the unit path scans a table that `fits`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum XScanByte {
    /// As before the switch (LOSAT_X_TBNSCANBF or NCBI's loop).
    Off,
    /// `x_scan_one_range_byte`; `avx2` = the AVX2 forms of both stages.
    On { avx2: bool },
    /// `x_scan_one_range_byte` and NCBI's loop on every call, compared
    /// (LOSAT_X_TBNSCANBYTESHADOW).
    Shadow { avx2: bool },
}

impl XScanByte {
    /// Whether `x_scan_one_range_byte` runs.
    pub(super) fn on(self) -> bool {
        self != XScanByte::Off
    }
}

// No NCBI counterpart: reads LOSAT_X_TBNSCANBYTE and LOSAT_X_TBNSCANBYTESHADOW once and detects
// AVX2; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the scan of the unit path. LOSAT_X_TBNSCANBYTE is a switch of
/// the strict set, on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set; `2` takes the
/// plain loops also where AVX2 is there. The shadow wins.
pub(super) fn x_scan_byte() -> XScanByte {
    static MODE: std::sync::OnceLock<XScanByte> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| {
        let shadow = std::env::var_os("LOSAT_X_TBNSCANBYTESHADOW").is_some();
        let value = crate::utils::xswitch::value("LOSAT_X_TBNSCANBYTE", 1..=2, 1);
        if !shadow && value.is_none() {
            return XScanByte::Off;
        }
        let avx2 = value != Some(2) && avx2_detected();
        if shadow {
            XScanByte::Shadow { avx2 }
        } else {
            XScanByte::On { avx2 }
        }
    })
}

/// Whether this machine has AVX2 (always false off x86-64).
fn avx2_detected() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): whether `x_scan_one_range_byte` can scan this table (word
/// size 3, 5-bit letters, a 15-bit mask).
pub(super) fn fits(lookup: &BlastAaLookupTable) -> bool {
    lookup.word_length == 3 && lookup.charsize == 5 && lookup.mask as usize == X_MASK
}

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence vector of `lookup` with one byte per index of
/// the 15-bit space: 1 where `pv` has the bit and the backbone has the cell, 0 elsewhere (a bit is
/// set only for a cell with hits, so this is `pv`; the second test keeps every candidate's cell
/// inside the backbone whatever `pv` holds).
fn build_presence(lookup: &BlastAaLookupTable) -> Arc<[u8]> {
    (0..=X_MASK)
        .map(|ix| {
            let bit = lookup
                .pv
                .get(ix >> PV_ARRAY_BTS)
                .map_or(0, |word| ((word >> (ix & PV_ARRAY_MASK)) & 1) as u8);
            if ix < lookup.backbone.len() {
                bit
            } else {
                0
            }
        })
        .collect()
}

thread_local! {
    /// The table `presence_bytes` built last on this thread and its bytes. The `Weak` keeps the
    /// table's allocation, so no other table can have its address while the entry is kept.
    static PRESENCE: std::cell::RefCell<(std::sync::Weak<BlastAaLookupTable>, Arc<[u8]>)> =
        std::cell::RefCell::new((std::sync::Weak::new(), Arc::from(Vec::new())));
}

// No NCBI counterpart: the presence vector in another layout, made once per table and thread; it
// does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence bytes of `lookup` (see `build_presence`), built
/// on this thread when it last saw another table. The table must `fit`.
pub(super) fn presence_bytes(lookup: &Arc<BlastAaLookupTable>) -> Arc<[u8]> {
    PRESENCE.with(|kept| {
        let mut kept = kept.borrow_mut();
        if kept.0.as_ptr() != Arc::as_ptr(lookup) {
            *kept = (Arc::downgrade(lookup), build_presence(lookup));
        }
        Arc::clone(&kept.1)
    })
}

thread_local! {
    /// The offset array of NCBI's loop under LOSAT_X_TBNSCANBYTESHADOW.
    static SHADOW_PAIRS: std::cell::RefCell<Vec<BlastOffsetPair>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

// No NCBI counterpart: LOSAT_X_STATS and shadow counters; they do not change any value NCBI
// computes.
/// Units scanned by `x_scan_one_range_byte` and their scan calls (`[0]`, `[1]`); calls and pairs
/// the shadow compared (`[2]`, `[3]`).
static COUNTS: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): counts one unit and its `calls` for the
/// `[X_STATS] TBNSCANBYTE` line (`LOSAT_X_STATS` set).
pub(super) fn x_count_unit(calls: u64) {
    if super::x_scan_bf::x_stats() {
        COUNTS[0].fetch_add(1, Ordering::Relaxed);
        COUNTS[1].fetch_add(calls, Ordering::Relaxed);
    }
}

/// The `[X_STATS] TBNSCANBYTE` line with `LOSAT_X_STATS` and the `[X_SHADOW] TBNSCANBYTE` line
/// with LOSAT_X_TBNSCANBYTESHADOW (printed by `main` at exit), when the unit path scanned any unit
/// this way.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    let [units, calls, compared, pairs] = COUNTS.each_ref().map(|c| c.load(Ordering::Relaxed));
    let mode = x_scan_byte_if_read();
    if units > 0 {
        let avx2 = matches!(
            mode,
            Some(XScanByte::On { avx2: true } | XScanByte::Shadow { avx2: true })
        );
        eprintln!("[X_STATS] TBNSCANBYTE units={units} calls={calls} avx2={avx2}");
    }
    if compared > 0 {
        eprintln!(
            "[X_SHADOW] TBNSCANBYTE calls_compared={compared} pairs_compared={pairs} (offset pairs and scan ranges identical)"
        );
    }
}

/// The mode if a search read it (the summary must not read the switch for programs that never
/// scanned).
fn x_scan_byte_if_read() -> Option<XScanByte> {
    (COUNTS[0].load(Ordering::Relaxed) + COUNTS[2].load(Ordering::Relaxed) > 0).then(x_scan_byte)
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (see `x_scan_bf.rs`).
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): `s_blast_aa_scan_subject_one_range` for a table that `fits`
/// (see the module text). Same arguments, same result, same `s_range`; `presence` is
/// `presence_bytes(lookup)` and `offset_pairs` must hold `array_size + X_SCAN_BYTE_SLACK` entries.
/// With `mode` the shadow, NCBI's loop runs on a second array and must give the same.
pub(super) fn x_scan_one_range_byte(
    mode: XScanByte,
    lookup: &BlastAaLookupTable,
    presence: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    match mode {
        XScanByte::On { avx2 } => scan(
            avx2,
            lookup,
            presence,
            subject,
            offset_pairs,
            array_size,
            s_range,
        ),
        XScanByte::Shadow { avx2 } => {
            let mut reference_range = *s_range;
            let (want, same_pairs) = SHADOW_PAIRS.with(|reference| {
                let mut reference = reference.borrow_mut();
                reference.clear();
                reference.resize(offset_pairs.len(), BlastOffsetPair::default());
                let want = s_blast_aa_scan_subject_one_range(
                    lookup,
                    subject,
                    &mut reference,
                    array_size,
                    &mut reference_range,
                );
                let got = scan(
                    avx2,
                    lookup,
                    presence,
                    subject,
                    offset_pairs,
                    array_size,
                    s_range,
                );
                let same = got == want
                    && offset_pairs[..got as usize]
                        .iter()
                        .zip(&reference[..want as usize])
                        .all(|(a, b)| a.q_off == b.q_off && a.s_off == b.s_off);
                (want, (got, same))
            });
            let (got, same) = same_pairs;
            let same = same && *s_range == reference_range;
            assert!(
                same,
                "LOSAT_X_TBNSCANBYTESHADOW: {got} pairs and s_range {s_range:?} vs {want} and {reference_range:?}"
            );
            COUNTS[2].fetch_add(1, Ordering::Relaxed);
            COUNTS[3].fetch_add(got as u64, Ordering::Relaxed);
            got
        }
        XScanByte::Off => unreachable!("x_scan_one_range_byte with the switch off"),
    }
}

/// The scan of one call (module text), with the AVX2 forms when `avx2`.
fn scan(
    avx2: bool,
    lookup: &BlastAaLookupTable,
    presence: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    if s_range[1] > s_range[2] {
        return 0;
    }
    assert!(
        fits(lookup) && presence.len() == X_MASK + 1,
        "LOSAT_X_TBNSCANBYTE: table of word size {} with {} presence bytes",
        lookup.word_length,
        presence.len()
    );
    assert!(
        array_size >= 0 && offset_pairs.len() >= array_size as usize + X_SCAN_BYTE_SLACK,
        "offset array of {} entries for array_size {array_size}",
        offset_pairs.len()
    );
    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;
    // NCBI's loop reads `subject[s + 2]` for every `s <= s_last`.
    assert!(
        s_last + 2 < subject.len(),
        "scan to {s_last} in {}",
        subject.len()
    );
    let table = Table {
        backbone: lookup.backbone.as_ptr(),
        overflow: lookup.overflow.as_ptr(),
    };
    let dest = offset_pairs.as_mut_ptr();
    let mut index = [0u16; X_BLOCK];
    let mut cand = [MaybeUninit::<u16>::uninit(); X_BLOCK];
    let mut totalhits: i32 = 0;
    let mut block_first = s_first;
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        let letters = &subject[block_first..block_last + 3];
        // SAFETY: `presence` has an entry for every 15-bit index and is 1 only at indexes of the
        // backbone (`build_presence`), so every candidate's cell is inside it; `dest` has
        // `array_size + X_SCAN_BYTE_SLACK` entries; AVX2 only when the caller detected it.
        #[cfg(target_arch = "x86_64")]
        let (total, stop) = unsafe {
            if avx2 {
                let n = candidates_avx2(letters, presence, &mut index, &mut cand);
                copy_avx2(
                    &table,
                    dest,
                    &index,
                    &cand,
                    n,
                    block_first,
                    totalhits,
                    array_size,
                )
            } else {
                let n = candidates(letters, presence, &mut index, &mut cand);
                copy(
                    &table,
                    dest,
                    &index,
                    &cand,
                    n,
                    block_first,
                    totalhits,
                    array_size,
                )
            }
        };
        // SAFETY: as above.
        #[cfg(not(target_arch = "x86_64"))]
        let (total, stop) = unsafe {
            let _ = avx2;
            let n = candidates(letters, presence, &mut index, &mut cand);
            copy(
                &table,
                dest,
                &index,
                &cand,
                n,
                block_first,
                totalhits,
                array_size,
            )
        };
        if let Some(s_off) = stop {
            // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:118-123
            s_range[1] = s_off as i32;
            return total;
        }
        totalhits = total;
        block_first = block_last + 1;
    }
    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:125-127
    s_range[1] = (s_last + 1) as i32;
    totalhits
}

/// The lookup table's backbone and overflow, for the copy stage.
struct Table {
    backbone: *const BackboneCell,
    overflow: *const i32,
}

// NCBI reference: c++/include/algo/blast/core/blast_lookup.h:96-107,121-127
// (ComputeTableIndex, ComputeTableIndexIncremental) and 55-57 (PV_TEST).
/// The plain-loop candidate stage: writes the index of every word of the block (word `k` has the
/// letters `letters[k..k + 3]`; `letters` holds the block's words and two letters more) to `index`
/// and lists the positions (from the block's first word) of the words whose presence byte is 1;
/// returns how many.
#[inline(never)]
fn candidates(
    letters: &[u8],
    presence: &[u8],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    let words = letters.len() - 2;
    let index = &mut index[..words];
    let (a, b, c) = (
        &letters[..words],
        &letters[1..words + 1],
        &letters[2..words + 2],
    );
    for k in 0..words {
        index[k] = ((a[k] as u16) << 10 | (b[k] as u16) << 5 | c[k] as u16) & X_MASK as u16;
    }
    list(index, presence, cand)
}

/// The candidates of `index` by their presence bytes, without a branch.
#[inline(always)]
fn list(index: &[u16], presence: &[u8], cand: &mut [MaybeUninit<u16>; X_BLOCK]) -> usize {
    let presence: &[u8; X_MASK + 1] = presence.try_into().expect("presence bytes");
    let out = cand.as_mut_ptr().cast::<u16>();
    let mut n = 0usize;
    for (k, &ix) in index.iter().enumerate() {
        let bit = presence[ix as usize & X_MASK] as usize;
        // SAFETY: `n <= k < index.len() <= X_BLOCK`.
        unsafe { out.add(n).write(k as u16) };
        n += bit;
    }
    n
}

/// The AVX2 candidate stage: `candidates` with the indexes made 16 at a time (three shifted loads
/// of 16 letters each); the same values.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn candidates_avx2(
    letters: &[u8],
    presence: &[u8],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    use std::arch::x86_64::*;
    let words = letters.len() - 2;
    let l = letters.as_ptr();
    let ixp = index.as_mut_ptr();
    let mask = _mm256_set1_epi16(X_MASK as i16);
    let mut k = 0usize;
    // Words k .. k + 16 read the letters k .. k + 18 (the last load ends at k + 18).
    while k + 16 <= words {
        let a = _mm256_cvtepu8_epi16(_mm_loadu_si128(l.add(k) as *const __m128i));
        let b = _mm256_cvtepu8_epi16(_mm_loadu_si128(l.add(k + 1) as *const __m128i));
        let c = _mm256_cvtepu8_epi16(_mm_loadu_si128(l.add(k + 2) as *const __m128i));
        let v = _mm256_or_si256(
            _mm256_or_si256(_mm256_slli_epi16::<10>(a), _mm256_slli_epi16::<5>(b)),
            c,
        );
        _mm256_storeu_si256(ixp.add(k) as *mut __m256i, _mm256_and_si256(v, mask));
        k += 16;
    }
    while k < words {
        let v = (*l.add(k) as u16) << 10 | (*l.add(k + 1) as u16) << 5 | *l.add(k + 2) as u16;
        *ixp.add(k) = v & X_MASK as u16;
        k += 1;
    }
    list(&index[..words], presence, cand)
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:99-123
// ```c
//         numhits = bbc[index].num_used;
//         if (numhits <= (array_size - totalhits)) {
//             if (numhits <= AA_HITS_PER_CELL) src = bbc[index].payload.entries;
//             else src = ovfl + bbc[index].payload.overflow_cursor;
//             for (i = 0; i < numhits; i++) {
//                 offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//                 offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
//             }
//             totalhits += numhits;
//         } else {
//             s_range[1] = (Int4)(s - subject->sequence);
//             return totalhits;
//         }
// ```
/// The tested copy of one candidate (NCBI's statements above): `None` when the word's hits do not
/// fit (the scan stops at this word), else the new total.
///
/// # Safety
/// `ix` is a candidate's index; `dest` has `array_size + X_SCAN_BYTE_SLACK` entries.
#[inline(always)]
unsafe fn copy_tested(
    table: &Table,
    dest: *mut BlastOffsetPair,
    ix: usize,
    s_off: u32,
    totalhits: i32,
    array_size: i32,
) -> Option<i32> {
    let cell = &*table.backbone.add(ix);
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
        let src = table.overflow.add(cell.entries[0] as usize);
        for i in 0..numhits as usize {
            *d.add(i) = BlastOffsetPair {
                q_off: *src.add(i) as u32,
                s_off,
            };
        }
    }
    Some(totalhits + numhits)
}

/// Whether three hits for each of the `left` candidates fit after `totalhits` (then no candidate
/// of up to `AA_HITS_PER_CELL` hits can stop the scan).
#[inline(always)]
fn room_for_rest(totalhits: i32, left: usize, array_size: i32) -> bool {
    totalhits as i64 + (AA_HITS_PER_CELL * left) as i64 <= array_size as i64
}

/// The plain-loop copy stage for the candidates `cand[..n]` of the block that starts at
/// `block_first`: the new total, and the subject offset of the word at which the offset array has
/// no room, if one has none.
///
/// # Safety
/// As `copy_tested`, for every candidate.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn copy(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    block_first: usize,
    mut totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    let cand = cand.as_ptr().cast::<u16>();
    let mut i = 0usize;
    while i < n {
        if room_for_rest(totalhits, n - i, array_size) {
            while i < n {
                let k = *cand.add(i) as usize;
                let cell = &*table.backbone.add(index[k] as usize);
                let numhits = cell.num_used;
                if numhits as usize > AA_HITS_PER_CELL {
                    break;
                }
                let s_off = (block_first + k) as u32;
                let d = dest.add(totalhits as usize);
                for j in 0..AA_HITS_PER_CELL {
                    *d.add(j) = BlastOffsetPair {
                        q_off: cell.entries[j] as u32,
                        s_off,
                    };
                }
                totalhits += numhits;
                i += 1;
            }
            if i == n {
                break;
            }
        }
        let k = *cand.add(i) as usize;
        let s_off = (block_first + k) as u32;
        match copy_tested(table, dest, index[k] as usize, s_off, totalhits, array_size) {
            Some(total) => totalhits = total,
            None => return (totalhits, Some(s_off)),
        }
        i += 1;
    }
    (totalhits, None)
}

/// The AVX2 copy stage: `copy` with the pairs of a cell of up to `AA_HITS_PER_CELL` hits written by
/// one 32-byte store (`[e0, s, e1, s, e2, s, e2, s]` from the cell `[num_used, e0, e1, e2]`); the
/// same values in the first `num_used` pairs.
///
/// # Safety
/// As `copy`; AVX2 must be there.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn copy_avx2(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    block_first: usize,
    mut totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    use std::arch::x86_64::*;
    const _: () = assert!(std::mem::size_of::<BackboneCell>() == 16 && AA_HITS_PER_CELL == 3);
    const _: () = assert!(std::mem::size_of::<BlastOffsetPair>() == 8);
    let order = _mm256_setr_epi32(1, 0, 2, 0, 3, 0, 3, 0);
    let cand = cand.as_ptr().cast::<u16>();
    let index = index.as_ptr();
    let mut i = 0usize;
    while i < n {
        if room_for_rest(totalhits, n - i, array_size) {
            let mut d = dest.add(totalhits as usize);
            while i < n {
                let k = *cand.add(i) as usize;
                let cell = table.backbone.add(*index.add(k) as usize);
                let numhits = (*cell).num_used;
                if numhits as usize > AA_HITS_PER_CELL {
                    break;
                }
                let c = _mm256_broadcastsi128_si256(_mm_loadu_si128(cell as *const __m128i));
                let s = _mm256_set1_epi32((block_first + k) as i32);
                let v = _mm256_blend_epi32::<0xAA>(_mm256_permutevar8x32_epi32(c, order), s);
                _mm256_storeu_si256(d as *mut __m256i, v);
                d = d.add(numhits as usize);
                totalhits += numhits;
                i += 1;
            }
            if i == n {
                break;
            }
        }
        let k = *cand.add(i) as usize;
        let s_off = (block_first + k) as u32;
        match copy_tested(
            table,
            dest,
            *index.add(k) as usize,
            s_off,
            totalhits,
            array_size,
        ) {
            Some(total) => totalhits = total,
            None => return (totalhits, Some(s_off)),
        }
        i += 1;
    }
    (totalhits, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::search_seed::encode_tblastn_lookup_query;
    use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;

    struct XorShift(u64);

    impl XorShift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

    fn lookup_for(queries: &[Vec<u8>], threshold: i32) -> BlastAaLookupTable {
        let frames: Vec<_> = queries
            .iter()
            .map(|query| vec![encode_tblastn_lookup_query(query, None, false)])
            .collect();
        let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        build_ncbi_lookup_for_profile(
            &frames,
            threshold,
            &karlin,
            false,
            ScoringMatrix::Blosum62,
            3,
        )
        .0
    }

    /// The scan calls of one range as `x_scan_unit` makes them: (offset pairs, s_range[1]) after
    /// each call; `avx2` = `None` for NCBI's loop.
    fn calls(
        lookup: &BlastAaLookupTable,
        presence: &[u8],
        subject: &[u8],
        range: (i32, i32),
        array_size: usize,
        avx2: Option<bool>,
    ) -> Vec<(Vec<(u32, u32)>, i32)> {
        let mut pairs = vec![BlastOffsetPair::default(); array_size + X_SCAN_BYTE_SLACK];
        let mut scan_range = [0, range.0, range.1];
        let mut out = Vec::new();
        while scan_range[1] <= scan_range[2] {
            let hits = match avx2 {
                Some(avx2) => scan(
                    avx2,
                    lookup,
                    presence,
                    subject,
                    &mut pairs,
                    array_size as i32,
                    &mut scan_range,
                ),
                None => s_blast_aa_scan_subject_one_range(
                    lookup,
                    subject,
                    &mut pairs,
                    array_size as i32,
                    &mut scan_range,
                ),
            };
            let got = pairs[..hits as usize]
                .iter()
                .map(|pair| (pair.q_off, pair.s_off))
                .collect();
            out.push((got, scan_range[1]));
            assert!(out.len() < 1_000_000, "scan does not advance");
        }
        out
    }

    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (the function compared
    // against, as `s_blast_aa_scan_subject_one_range`). Not a port of an NCBI unit case: the pairs
    // and the `s_range[1]` of every call must equal the reference's, with offset arrays that stop
    // the scan after a few words (as aascan_unit_test.cpp's ScanOffsetTest scans), cells past
    // `AA_HITS_PER_CELL` hits, every letter value of ncbistdaa, and block ends at every alignment.
    #[test]
    fn x_scan_byte_equals_ncbi_scan_on_random_subjects() {
        let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
        let forms: Vec<bool> = if avx2_detected() {
            vec![false, true]
        } else {
            vec![false]
        };
        let mut compared = 0usize;
        let mut stops = 0usize;
        let mut big_cells = 0usize;
        for case in 0..24 {
            let num_queries = [1, 1, 3, 12][case % 4];
            let queries: Vec<Vec<u8>> = (0..num_queries)
                .map(|_| {
                    let len = 20 + rng.below(380);
                    (0..len).map(|_| AMINO[rng.below(AMINO.len())]).collect()
                })
                .collect();
            let threshold = [11, 13, 15][case % 3];
            let lookup = lookup_for(&queries, threshold);
            assert!(fits(&lookup), "tblastn's table must fit");
            let presence = build_presence(&lookup);
            big_cells += lookup
                .backbone
                .iter()
                .filter(|cell| cell.num_used as usize > AA_HITS_PER_CELL)
                .count();
            let longest = lookup.longest_chain.max(1) as usize;
            // A frame in ncbistdaa (0..=27; a translated subject has 1..=25), sentinel 0 after it.
            let len = [3, 4, 255, 256, 257, 271, 1_000, 20_000][case % 8];
            let mut subject: Vec<u8> = (0..len).map(|_| rng.below(28) as u8).collect();
            subject.push(0);
            let last = len as i32 - 3;
            let ranges = [
                (0, last.max(0)),
                (1.min(last.max(0)), last.max(0)),
                (len as i32 / 3, last.max(len as i32 / 3)),
                (5.min(last.max(0)), (last - 7).max(5.min(last.max(0)))),
            ];
            for &range in &ranges {
                for array_size in [
                    longest,
                    longest + 1,
                    longest + 2,
                    3 * longest,
                    3 * longest + 1,
                    longest * 1024,
                ] {
                    let want = calls(&lookup, &presence, &subject, range, array_size, None);
                    for &avx2 in &forms {
                        let got =
                            calls(&lookup, &presence, &subject, range, array_size, Some(avx2));
                        assert_eq!(
                            got, want,
                            "case {case} range {range:?} array_size {array_size} avx2 {avx2}"
                        );
                    }
                    compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                    stops += want.len().saturating_sub(1);
                }
            }
        }
        // The cases above must reach the capacity stop, many hits and cells of the overflow.
        assert!(compared > 100_000, "{compared} pairs compared");
        assert!(stops > 100, "{stops} capacity stops");
        assert!(big_cells > 0, "no cell past {AA_HITS_PER_CELL} hits");
    }

    // No NCBI counterpart: the index made from the letters is the rolled one
    // (`ComputeTableIndexIncremental` from `ComputeTableIndex`) for every byte value, and the
    // presence bytes are the bits of `pv`.
    #[test]
    fn indexes_and_presence_equal_the_rolled_index_and_pv() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let letters: Vec<u8> = (0..X_BLOCK + 2).map(|_| rng.next() as u8).collect();
        let lookup = lookup_for(&[b"MKTAYIAKQRQISFVKSHFSRQ".to_vec()], 11);
        let presence = build_presence(&lookup);
        assert_eq!(presence.len(), X_MASK + 1);
        for ix in 0..=X_MASK {
            let bit = lookup
                .pv
                .get(ix >> PV_ARRAY_BTS)
                .is_some_and(|w| (w >> (ix & PV_ARRAY_MASK)) & 1 == 1);
            assert_eq!(presence[ix] == 1, bit, "index {ix}");
            assert!(
                presence[ix] == 0 || ix < lookup.backbone.len(),
                "index {ix}"
            );
        }
        for words in [1, 15, 16, 17, 31, 32, 33, 100, X_BLOCK] {
            let block = &letters[..words + 2];
            let mut forms = vec![false];
            if avx2_detected() {
                forms.push(true);
            }
            for avx2 in forms {
                let mut index = [0u16; X_BLOCK];
                let mut cand = [MaybeUninit::<u16>::uninit(); X_BLOCK];
                let n = run_candidates(avx2, block, &presence, &mut index, &mut cand);
                let mut rolled = ((block[0] as usize) << 5 | block[1] as usize) & X_MASK;
                let mut want = Vec::new();
                for k in 0..words {
                    rolled = ((rolled << 5) | block[k + 2] as usize) & X_MASK;
                    assert_eq!(
                        index[k] as usize, rolled,
                        "word {k} of {words}, avx2 {avx2}"
                    );
                    if presence[rolled] == 1 {
                        want.push(k as u16);
                    }
                }
                // SAFETY: the first `n` entries were written.
                let got: Vec<u16> = cand[..n]
                    .iter()
                    .map(|c| unsafe { c.assume_init() })
                    .collect();
                assert_eq!(got, want, "{words} words, avx2 {avx2}");
            }
        }
    }

    /// `candidates_avx2` when `avx2` (detected by the caller), else `candidates`.
    fn run_candidates(
        avx2: bool,
        letters: &[u8],
        presence: &[u8],
        index: &mut [u16; X_BLOCK],
        cand: &mut [MaybeUninit<u16>; X_BLOCK],
    ) -> usize {
        #[cfg(target_arch = "x86_64")]
        if avx2 {
            // SAFETY: AVX2 was detected.
            return unsafe { candidates_avx2(letters, presence, index, cand) };
        }
        assert!(!avx2, "AVX2 off x86-64");
        candidates(letters, presence, index, cand)
    }

    // No NCBI counterpart: the kept presence bytes follow the table (a new table, also one at a
    // freed table's address, gets its own).
    #[test]
    fn presence_bytes_are_rebuilt_for_another_table() {
        let one = Arc::new(lookup_for(&[b"MKTAYIAKQRQISFVKSHFSRQ".to_vec()], 11));
        let two = Arc::new(lookup_for(&[b"WWWWCCCCHHHHMMMM".to_vec()], 11));
        assert_eq!(*presence_bytes(&one), *build_presence(&one));
        assert_eq!(*presence_bytes(&two), *build_presence(&two));
        assert_eq!(*presence_bytes(&one), *build_presence(&one));
        drop(one);
        let three = Arc::new(lookup_for(&[b"PPPPGGGGAAAA".to_vec()], 11));
        assert_eq!(*presence_bytes(&three), *build_presence(&three));
    }
}
