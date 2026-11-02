//! EXPERIMENT (LOSAT_X_TBNSCANBYTE): tblastn's lookup scan with the table indexes of a block of
//! words made at once and their presence read from one byte per index, and the hits of a block
//! copied with one capacity test.
//!
//! NCBI's scan (c++/src/algo/blast/core/blast_aascan.c:75-127, `s_BlastAaScanSubject`, kept as
//! `tblastx/blast_aascan.rs` `s_blast_aa_scan_subject_one_range`) rolls the index of the next word
//! from the one before (`ComputeTableIndexIncremental`), tests its presence-vector bit and copies
//! the word's hits when the offset array has room. The branch-free scan (`x_scan_bf.rs`) takes the
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
//! The copy then takes the block's candidates in position order, as NCBI's loop does, each once.
//! While three hits for every candidate left in the block fit in the offset array (one test, made
//! again only after a cell of more hits), the hits of a cell of up to `AA_HITS_PER_CELL` are
//! written without the capacity test, as three pairs of one 8-byte store each (with AVX2 four pairs
//! in one 32-byte store: the cell broadcast, one `vpermd`, one `vpblendd` with the subject offset)
//! of which the first `num_used` count. A cell of more hits, and every candidate met while the
//! block could fill the array, takes NCBI's test; a larger cell's hits are then written four pairs
//! at a time from the overflow array. When a word's hits do not fit, `s_range[1]` is that word's
//! position and the scan returns, so the next call resumes there. The offset pairs and the
//! `s_range` values of every call are NCBI's (the words after the stop were looked at but not
//! used; the next call looks at them again, as with the branch-free scan).
//!
//! EXPERIMENT (LOSAT_X_TBNSCANRESUME): the unit path keeps the block where a call stops (its
//! indexes, its candidates and the word that did not fit, `XScanHeld`), and the next call of the
//! same range copies on from that candidate instead of making the block's candidates again
//! (`x_scan_one_range_held`; the copy loop is the one above). A stop then costs a call, not a
//! block of words, and the unit path gives the offset array NCBI's size, 4096 + the longest chain
//! (`search_seed.rs` `x_scan_unit`; `lookup_wrap.c:264-266`; at most the plan's size), instead of
//! `LOSAT_X_TBNPAIRCAP`'s 32768 + it. The pairs and the `s_range` of every call are still NCBI's.
//!
//! The offset array must have `X_SCAN_BYTE_SLACK` entries after `array_size`; the entries past the
//! returned count are not part of the result. `presence` depends on the lookup table only; each
//! thread builds it once per table (`presence_bytes`).
//!
//! Probe and source of the form: line F s3 P3 V2 (candidates) and s5 P3b arm 7 (copy),
//! `exp/f2-tbn-scan` ff918e61: HsG_p1 x hsap_dna `-db` n1 whole-run cycles x0.703 of dev, with
//! the offset pairs and `s_range` of every call and the output equal.
//!
//! Both stages have an AVX2 form and a plain-loop form (the same values; the plain index loop is
//! left to the compiler's vectorizer). On the workstation (one core, the scan before this switch =
//! the branch-free scan) the AVX2 forms take x0.75 of the time of one protein x the human genome and
//! x0.98-0.99 where most words are candidates (many proteins x a genome, 55 proteins x 16S); the
//! plain loops take x0.79 and x0.83 (the LOSAT Web serial reactor, one protein x E. coli) of the
//! time where few words are candidates, but x1.01-1.06 where most are: their copy of a cell's hits
//! costs more than the copy before. So by default the switch takes the AVX2 forms where the CPU
//! has AVX2, and elsewhere the plain loops for sparse tables only: those where at most three tenths
//! of the backbone cells have hits (`X_SPARSE_TENTHS`; known when the presence bytes are built).
//! The first k proteins of AP027131 x E. coli (one core): the plain loops against the scan before
//! took x0.84-0.95 in the LOSAT Web serial reactor and cycles x0.84-0.90 on x86-64 at a fill of
//! 0.09-0.24 (k = 1-8), x1.00 and x0.95 at 0.34 (k = 16), and x1.02-1.03 and x0.99-1.00 at 0.42-0.46.
//!
//! The switch is read once. It is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_TBNSCANBYTE=0`; on, it takes the AVX2 forms
//! of both stages on x86-64 with AVX2, the plain loops for sparse tables on x86-64 without AVX2 and
//! in WebAssembly, and changes nothing elsewhere (aarch64: not timed). `LOSAT_X_TBNSCANBYTE=2`
//! takes the plain loops on any CPU for every table, `=3` the plain loops for sparse tables on any
//! CPU. A table of another shape, and without AVX2 a table that is not sparse, is scanned as before
//! (the branch-free scan of `x_scan_bf.rs`). With `LOSAT_X_STATS` set, `main` prints
//! `[X_STATS] TBNSCANBYTE units= calls= avx2=` (the units scanned this way, their scan calls, and
//! whether the AVX2 forms ran).
//!
//! EXPERIMENT (LOSAT_X_BXSCANBYTE): blastx's direct seed loop (`blastx/x_seed_scan.rs`) scans its
//! standard table (word size 3, 5-bit letters: every one `fits`) with `x_scan_one_range_byte` too,
//! under its own switch of the same values and CPU policy (`x_mode`). Its presence bytes are built
//! once per table and kept with it (`x_table_presence`; blastx's `Lookup`).

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
pub(crate) const X_SCAN_BYTE_SLACK: usize = AA_HITS_PER_CELL + 1;

/// A table is sparse when at most `X_SPARSE_TENTHS` tenths of its backbone cells have hits
/// (presence bytes of 1): the plain loops then win (`XScanByte::Sparse`).
const X_SPARSE_TENTHS: usize = 3;

/// Whether CPUs without AVX2 take the plain loops for sparse tables by default: on the targets they
/// were timed on (x86-64 without AVX2, WebAssembly); elsewhere (aarch64) only with `=3`.
const X_SPARSE_DEFAULT: bool = cfg!(any(target_arch = "x86_64", target_arch = "wasm32"));

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): how the unit path scans a table that `fits`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XScanByte {
    /// As before the switch (the branch-free scan of `x_scan_bf.rs`).
    Off,
    /// `x_scan_one_range_byte`; `avx2` = the AVX2 forms of both stages.
    On { avx2: bool },
    /// `x_scan_one_range_byte` with the plain loops for a sparse table (`presence_bytes`), as
    /// before for the others.
    Sparse,
}

impl XScanByte {
    /// Whether `x_scan_one_range_byte` runs.
    pub(crate) fn on(self) -> bool {
        self != XScanByte::Off
    }
}

// No NCBI counterpart: reads LOSAT_X_TBNSCANBYTE once and detects AVX2; it does not change any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the scan of the unit path. LOSAT_X_TBNSCANBYTE is a switch of
/// the strict set, on unless it is `0`. On, it takes the
/// AVX2 forms where the CPU has AVX2, and elsewhere the plain loops for sparse tables only (the
/// plain loops lose where most words are candidates, see the module text; on aarch64 nothing);
/// `2` takes the plain loops on any CPU and table, `3` the plain loops for sparse tables on any
/// CPU.
pub(super) fn x_scan_byte() -> XScanByte {
    static MODE: std::sync::OnceLock<XScanByte> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| x_mode("LOSAT_X_TBNSCANBYTE"))
}

// No NCBI counterpart: reads a switch and detects AVX2; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE, LOSAT_X_BXSCANBYTE): the scan that the strict-set switch `name`
/// (one of these two) asks for on this CPU, read now (the callers keep it): `x_scan_byte`'s values.
pub(crate) fn x_mode(name: &str) -> XScanByte {
    let value = crate::utils::xswitch::value(name, 1..=3, 1);
    let avx2 = !matches!(value, Some(2 | 3)) && avx2_detected();
    if value == Some(2) || (value.is_some() && avx2) {
        XScanByte::On { avx2 }
    } else if value == Some(3) || (value == Some(1) && X_SPARSE_DEFAULT) {
        XScanByte::Sparse
    } else {
        XScanByte::Off
    }
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
    /// The table `presence_bytes` built last on this thread, its bytes and whether it is sparse.
    /// The `Weak` keeps the table's allocation, so no other table can have its address while the
    /// entry is kept.
    #[allow(clippy::type_complexity)]
    static PRESENCE: std::cell::RefCell<(std::sync::Weak<BlastAaLookupTable>, Arc<[u8]>, bool)> =
        std::cell::RefCell::new((std::sync::Weak::new(), Arc::from(Vec::new()), false));
}

/// Frees the presence bytes this thread keeps, at the end of its loop over the chunks of the
/// subject flow (`x_oid_flow::x_loop_free`, LOSAT_X_TBNLOOPFREE). The next table this thread scans
/// with has them built again, as every query batch's table has.
#[cfg(feature = "parallel")]
pub(super) fn x_free_presence() {
    let _ = PRESENCE.try_with(|kept| {
        if let Ok(mut kept) = kept.try_borrow_mut() {
            *kept = (std::sync::Weak::new(), Arc::from(Vec::new()), false);
        }
    });
}

// No NCBI counterpart: the presence vector in another layout, made once per table and thread; it
// does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence bytes of `lookup` (see `build_presence`), built
/// on this thread when it last saw another table. The table must `fit`.
pub(super) fn presence_bytes(lookup: &Arc<BlastAaLookupTable>) -> (Arc<[u8]>, bool) {
    PRESENCE.with(|kept| {
        let mut kept = kept.borrow_mut();
        if kept.0.as_ptr() != Arc::as_ptr(lookup) {
            let bytes = build_presence(lookup);
            let sparse = sparse(&bytes, lookup);
            *kept = (Arc::downgrade(lookup), bytes, sparse);
        }
        (Arc::clone(&kept.1), kept.2)
    })
}

/// Whether at most `X_SPARSE_TENTHS` tenths of the backbone cells of `lookup` have hits.
fn sparse(presence: &[u8], lookup: &BlastAaLookupTable) -> bool {
    let ones = presence.iter().filter(|&&b| b == 1).count();
    ones * 10 <= lookup.backbone.len() * X_SPARSE_TENTHS
}

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence bytes for `x_scan_one_range_byte` when `mode`
/// scans `lookup` this way (it is on, the table `fits`, and it is sparse if `mode` asks), else
/// `None` (the scan before).
pub(super) fn x_presence(mode: XScanByte, lookup: &Arc<BlastAaLookupTable>) -> Option<Arc<[u8]>> {
    if !mode.on() || !fits(lookup) {
        return None;
    }
    let (bytes, sparse) = presence_bytes(lookup);
    (mode != XScanByte::Sparse || sparse).then_some(bytes)
}

// No NCBI counterpart: the presence vector in another layout, made once per table; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXSCANBYTE): the presence bytes of `lookup` (see `build_presence`) and
/// whether it is sparse, for a caller that keeps them with the table (blastx's `Lookup`); `None`
/// for a table that does not `fit`.
pub(crate) fn x_table_presence(lookup: &BlastAaLookupTable) -> Option<(Arc<[u8]>, bool)> {
    fits(lookup).then(|| {
        let bytes = build_presence(lookup);
        let sparse = sparse(&bytes, lookup);
        (bytes, sparse)
    })
}

/// EXPERIMENT (LOSAT_X_BXSCANBYTE): the bytes of `kept` (from `x_table_presence`) when `mode`
/// scans that table on presence bytes (it is on, and the table is sparse if `mode` asks), else
/// `None` (the scan before).
pub(crate) fn x_kept_presence(mode: XScanByte, kept: &Option<(Arc<[u8]>, bool)>) -> Option<&[u8]> {
    match kept {
        Some((bytes, sparse)) if mode.on() && (mode != XScanByte::Sparse || *sparse) => Some(bytes),
        _ => None,
    }
}

// No NCBI counterpart: LOSAT_X_STATS counters; they do not change any value NCBI computes.
/// Units scanned by `x_scan_one_range_byte` and their scan calls (`[0]`, `[1]`).
static COUNTS: [AtomicU64; 2] = [const { AtomicU64::new(0) }; 2];

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): counts one unit and its `calls` for the
/// `[X_STATS] TBNSCANBYTE` line (`LOSAT_X_STATS` set).
pub(super) fn x_count_unit(calls: u64) {
    if super::x_scan_bf::x_stats() {
        COUNTS[0].fetch_add(1, Ordering::Relaxed);
        COUNTS[1].fetch_add(calls, Ordering::Relaxed);
    }
}

/// The `[X_STATS] TBNSCANBYTE` line (printed by `main` at exit) with `LOSAT_X_STATS`, when the unit
/// path scanned any unit this way.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    let [units, calls] = COUNTS.each_ref().map(|c| c.load(Ordering::Relaxed));
    let mode = x_scan_byte_if_read();
    if units > 0 {
        let avx2 = matches!(mode, Some(XScanByte::On { avx2: true }));
        eprintln!("[X_STATS] TBNSCANBYTE units={units} calls={calls} avx2={avx2}");
    }
}

/// The mode if a search read it (the summary must not read the switch for programs that never
/// scanned).
fn x_scan_byte_if_read() -> Option<XScanByte> {
    (COUNTS[0].load(Ordering::Relaxed) > 0).then(x_scan_byte)
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (see `x_scan_bf.rs`).
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): `s_blast_aa_scan_subject_one_range` for a table that `fits`
/// (see the module text). Same arguments, same result, same `s_range`; `presence` is
/// `presence_bytes(lookup)` and `offset_pairs` must hold `array_size + X_SCAN_BYTE_SLACK` entries.
/// Not inlined: the switch-off scan of `x_scan_unit` keeps its code. Blastx's direct loop calls it
/// too (LOSAT_X_BXSCANBYTE).
#[inline(never)]
pub(crate) fn x_scan_one_range_byte(
    mode: XScanByte,
    lookup: &BlastAaLookupTable,
    presence: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    x_scan_one_range_held(
        mode,
        lookup,
        presence,
        subject,
        offset_pairs,
        array_size,
        s_range,
        &mut XScanHeld::new(),
    )
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (see `x_scan_bf.rs`).
/// EXPERIMENT (LOSAT_X_TBNSCANRESUME): `x_scan_one_range_byte` that keeps in `held` the block of a
/// call that stops on a full offset array, and copies the rest of it when the next call is that
/// call's resumption (the same subject and `s_range[2]`, `s_range[1]` the word that did not fit).
/// Same result and `s_range` as NCBI's scan. A new `held` holds nothing.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
pub(crate) fn x_scan_one_range_held(
    mode: XScanByte,
    lookup: &BlastAaLookupTable,
    presence: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
    held: &mut XScanHeld,
) -> i32 {
    match mode {
        XScanByte::On { .. } | XScanByte::Sparse => scan(
            mode == XScanByte::On { avx2: true },
            lookup,
            presence,
            subject,
            offset_pairs,
            array_size,
            s_range,
            held,
        ),
        XScanByte::Off => unreachable!("x_scan_one_range_byte with the switch off"),
    }
}

/// EXPERIMENT (LOSAT_X_TBNSCANRESUME): the block of the scan call that stopped on a full offset
/// array, for its resumption (`x_scan_one_range_held`): the indexes and candidates of the block,
/// the first candidate not copied, and the call that resumes it.
pub(crate) struct XScanHeld {
    index: [u16; X_BLOCK],
    cand: [MaybeUninit<u16>; X_BLOCK],
    /// The block's candidates (`cand[..n]`); a block is held while `next < n`.
    n: usize,
    /// The first candidate not copied.
    next: usize,
    /// The subject offset of the block's first word.
    block_first: usize,
    /// The resuming call: the subject's address and length, `s_range[1]` and `s_range[2]`.
    at: (usize, usize, i32, i32),
}

impl XScanHeld {
    /// One that holds nothing.
    pub(crate) fn new() -> Self {
        XScanHeld {
            index: [0; X_BLOCK],
            cand: [MaybeUninit::uninit(); X_BLOCK],
            n: 0,
            next: 0,
            block_first: 0,
            at: (0, 0, 0, 0),
        }
    }
}

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANRESUME): whether the unit path's scan on presence bytes resumes a
/// stopped call in its block (`x_scan_one_range_held`) and takes NCBI's offset-array size
/// (`x_scan_unit`). In the strict set: on unless `LOSAT_X_TBNSCANRESUME=0`; read once.
pub(super) fn x_scan_resume() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNSCANRESUME"))
}

/// The scan of one call (module text), with the AVX2 forms when `avx2`; the block where it stops is
/// kept in `held`, and copied on from there when this call resumes it.
#[allow(clippy::too_many_arguments)]
fn scan(
    avx2: bool,
    lookup: &BlastAaLookupTable,
    presence: &[u8],
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
    held: &mut XScanHeld,
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
    #[cfg(target_arch = "x86_64")]
    debug_assert!(!avx2 || avx2_detected(), "AVX2 forms without AVX2");
    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;
    // NCBI's loop reads `subject[s + 2]` for every `s <= s_last`. A first range of one letter at
    // the end of the slice (`x_scan_unit` sets its end to its left edge, aa_ungapped.c:500-501)
    // reads one byte past it: that call takes NCBI's loop, as without this switch.
    if s_last + 2 >= subject.len() {
        held.n = 0;
        return s_blast_aa_scan_subject_one_range(
            lookup,
            subject,
            offset_pairs,
            array_size,
            s_range,
        );
    }
    let table = Table {
        backbone: lookup.backbone.as_ptr(),
        overflow: lookup.overflow.as_ptr(),
        overflow_len: lookup.overflow.len(),
    };
    let dest = offset_pairs.as_mut_ptr();
    let mut totalhits: i32 = 0;
    let mut block_first = s_first;
    // The held block, from its first candidate not copied, when this call resumes the one that
    // stopped in it (its first candidate is the word at `s_range[1]`).
    let subject_at = (subject.as_ptr() as usize, subject.len());
    let mut first = 0;
    if held.next < held.n && held.at == (subject_at.0, subject_at.1, s_range[1], s_range[2]) {
        first = held.next;
        block_first = held.block_first;
    }
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        if first == 0 {
            let letters = &subject[block_first..block_last + 3];
            held.n = candidates_any(avx2, letters, presence, &mut held.index, &mut held.cand);
        }
        // SAFETY: `presence` has an entry for every 15-bit index and is 1 only at indexes of the
        // backbone (`build_presence`), so every candidate's cell is inside it; `dest` has
        // `array_size + X_SCAN_BYTE_SLACK` entries; AVX2 only when the caller detected it.
        let (total, stop) = unsafe {
            copy_any(
                avx2,
                &table,
                dest,
                held,
                first,
                block_first,
                totalhits,
                array_size,
            )
        };
        if let Some(s_off) = stop {
            // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:118-123
            s_range[1] = s_off as i32;
            // The word's place among the block's candidates (which ascend), for the next call.
            // SAFETY: `cand[..n]` was written by the candidate stage.
            let cand =
                unsafe { &*(&held.cand[..held.n] as *const [MaybeUninit<u16>] as *const [u16]) };
            held.next = cand.partition_point(|&k| (k as usize) < s_off as usize - block_first);
            held.block_first = block_first;
            held.at = (subject_at.0, subject_at.1, s_range[1], s_range[2]);
            return total;
        }
        totalhits = total;
        first = 0;
        block_first = block_last + 1;
    }
    held.n = 0;
    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:125-127
    s_range[1] = (s_last + 1) as i32;
    totalhits
}

/// The candidate stage of a block (`candidates`, or `candidates_avx2` when `avx2`).
#[inline(always)]
fn candidates_any(
    avx2: bool,
    letters: &[u8],
    presence: &[u8],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: AVX2 only when the caller detected it.
        return unsafe { candidates_avx2(letters, presence, index, cand) };
    }
    let _ = avx2;
    candidates(letters, presence, index, cand)
}

/// The copy stage of `held`'s block from its candidate `first` (`copy`, or `copy_avx2` when
/// `avx2`).
///
/// # Safety
/// As `copy_block`; AVX2 must be there when `avx2`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn copy_any(
    avx2: bool,
    table: &Table,
    dest: *mut BlastOffsetPair,
    held: &XScanHeld,
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        return copy_avx2(
            table,
            dest,
            &held.index,
            &held.cand,
            held.n,
            first,
            block_first,
            totalhits,
            array_size,
        );
    }
    let _ = avx2;
    copy(
        table,
        dest,
        &held.index,
        &held.cand,
        held.n,
        first,
        block_first,
        totalhits,
        array_size,
    )
}

/// The lookup table's backbone and overflow, for the copy stage.
struct Table {
    backbone: *const BackboneCell,
    overflow: *const i32,
    overflow_len: usize,
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

/// Whether three hits for each of the `left` candidates fit after `totalhits` (then no candidate
/// of up to `AA_HITS_PER_CELL` hits can stop the scan).
#[inline(always)]
fn room_for_rest(totalhits: i32, left: usize, array_size: i32) -> bool {
    totalhits as i64 + (AA_HITS_PER_CELL * left) as i64 <= array_size as i64
}

/// The plain-loop copy stage for the candidates `cand[..n]` of the block that starts at
/// `block_first` (`copy_block`).
///
/// # Safety
/// As `copy_block`.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn copy(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    copy_block::<false>(
        table,
        dest,
        index,
        cand,
        n,
        first,
        block_first,
        totalhits,
        array_size,
    )
}

/// The AVX2 copy stage (`copy_block` with the AVX2 writes).
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
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    copy_block::<true>(
        table,
        dest,
        index,
        cand,
        n,
        first,
        block_first,
        totalhits,
        array_size,
    )
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
/// The copy stage for the candidates `cand[first..n]` of the block that starts at `block_first`
/// (NCBI's statements above, in NCBI's order): the new total, and the subject offset of the word
/// at which the offset array has no room, if one has none. While three hits for every candidate left fit, a cell of up to `AA_HITS_PER_CELL` hits is
/// written without NCBI's test (after it they still fit); a larger cell takes the test, and the
/// room is computed again after it. A candidate met without that room takes the test whatever its
/// size. `AVX2` = the AVX2 writes (`write`).
///
/// # Safety
/// The candidates' indexes are those of backbone cells of `table`; `dest` has `array_size +
/// X_SCAN_BYTE_SLACK` entries; AVX2 must be there when `AVX2`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn copy_block<const AVX2: bool>(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    first: usize,
    block_first: usize,
    mut totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    let cand = cand.as_ptr().cast::<u16>();
    let index = index.as_ptr();
    let mut i = first;
    while i < n {
        if room_for_rest(totalhits, n - i, array_size) {
            while i < n {
                let k = *cand.add(i) as usize;
                let cell = table.backbone.add(*index.add(k) as usize);
                let numhits = (*cell).num_used;
                let s_off = (block_first + k) as u32;
                i += 1;
                if numhits as usize <= AA_HITS_PER_CELL {
                    write::<AVX2>(table, dest.add(totalhits as usize), cell, s_off);
                    totalhits += numhits;
                } else {
                    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:100,118-123
                    if numhits > array_size - totalhits {
                        return (totalhits, Some(s_off));
                    }
                    write::<AVX2>(table, dest.add(totalhits as usize), cell, s_off);
                    totalhits += numhits;
                    if !room_for_rest(totalhits, n - i, array_size) {
                        break;
                    }
                }
            }
        } else {
            let k = *cand.add(i) as usize;
            let cell = table.backbone.add(*index.add(k) as usize);
            let numhits = (*cell).num_used;
            let s_off = (block_first + k) as u32;
            // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:100,118-123
            if numhits > array_size - totalhits {
                return (totalhits, Some(s_off));
            }
            write::<AVX2>(table, dest.add(totalhits as usize), cell, s_off);
            totalhits += numhits;
            i += 1;
        }
    }
    (totalhits, None)
}

/// One offset pair as the 8 bytes it is in memory (`q_off`, then `s_off`: `BlastOffsetPair` is
/// `repr(C)`), to write it with one store.
#[inline(always)]
fn pair_bits(q_off: u32, s_off: u32) -> u64 {
    if cfg!(target_endian = "little") {
        q_off as u64 | (s_off as u64) << 32
    } else {
        (q_off as u64) << 32 | s_off as u64
    }
}

/// Writes the hits of `cell` (subject offset `s_off`) at `d`: a cell of up to `AA_HITS_PER_CELL`
/// hits as `AA_HITS_PER_CELL` pairs (with AVX2 `AA_HITS_PER_CELL + 1`), a larger one from the
/// overflow array four pairs at a time. The pairs past the cell's hits are not part of the result;
/// a larger cell writes at most three of them, and reads past its query offsets only while they
/// stay inside the overflow array (else its last pairs are written one by one).
///
/// # Safety
/// `cell` is a backbone cell of `table`; `d` has room for the cell's hits and
/// `X_SCAN_BYTE_SLACK` entries more; AVX2 must be there when `AVX2`.
#[inline(always)]
unsafe fn write<const AVX2: bool>(
    table: &Table,
    d: *mut BlastOffsetPair,
    cell: *const BackboneCell,
    s_off: u32,
) {
    #[cfg(target_arch = "x86_64")]
    if AVX2 {
        return write_avx2(table, d, cell, s_off);
    }
    let numhits = (*cell).num_used as usize;
    let d = d.cast::<u64>();
    if numhits <= AA_HITS_PER_CELL {
        for j in 0..AA_HITS_PER_CELL {
            d.add(j)
                .write_unaligned(pair_bits((*cell).entries[j] as u32, s_off));
        }
    } else {
        let cursor = (*cell).entries[0] as usize;
        let src = table.overflow.add(cursor);
        let whole = chain_whole(table, cursor, numhits);
        // SAFETY: the first `whole` query offsets are inside the overflow array (`chain_whole`) and
        // `d` has room for `whole <= numhits + 3` pairs; the two do not overlap.
        chain_fours(
            std::slice::from_raw_parts_mut(d.cast::<u8>(), 8 * whole),
            std::slice::from_raw_parts(src, whole),
            s_off,
        );
        for j in whole..numhits {
            d.add(j)
                .write_unaligned(pair_bits(*src.add(j) as u32, s_off));
        }
    }
}

/// Writes the pairs of the query offsets `src` (a multiple of four) with the subject offset
/// `s_off` to `dst` (8 bytes per pair), four at a time. The two slices are separate arguments so
/// that the compiler knows they do not overlap and can widen the four without a runtime test.
#[inline(always)]
fn chain_fours(dst: &mut [u8], src: &[i32], s_off: u32) {
    for (out, q) in dst.chunks_exact_mut(32).zip(src.chunks_exact(4)) {
        for l in 0..4 {
            out[8 * l..8 * l + 8].copy_from_slice(&pair_bits(q[l] as u32, s_off).to_ne_bytes());
        }
    }
}

/// The pairs of an overflow chain of `numhits` hits at `cursor` written four at a time: all of
/// them, rounded up to four, when those query offsets are inside the overflow array, else the
/// whole fours.
#[inline(always)]
fn chain_whole(table: &Table, cursor: usize, numhits: usize) -> usize {
    if cursor + numhits.next_multiple_of(4) <= table.overflow_len {
        numhits.next_multiple_of(4)
    } else {
        numhits & !3
    }
}

/// `write` with AVX2: a cell of up to `AA_HITS_PER_CELL` hits as one 32-byte store
/// (`[e0, s, e1, s, e2, s, e2, s]` from the cell `[num_used, e0, e1, e2]`), a larger one as 32-byte
/// stores of four query offsets widened to 64 bits with `s_off` in the high halves.
///
/// # Safety
/// As `write`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn write_avx2(
    table: &Table,
    d: *mut BlastOffsetPair,
    cell: *const BackboneCell,
    s_off: u32,
) {
    use std::arch::x86_64::*;
    const _: () = assert!(std::mem::size_of::<BackboneCell>() == 16 && AA_HITS_PER_CELL == 3);
    const _: () =
        assert!(std::mem::size_of::<BlastOffsetPair>() == 8 && cfg!(target_endian = "little"));
    let numhits = (*cell).num_used as usize;
    if numhits <= AA_HITS_PER_CELL {
        let order = _mm256_setr_epi32(1, 0, 2, 0, 3, 0, 3, 0);
        let c = _mm256_broadcastsi128_si256(_mm_loadu_si128(cell as *const __m128i));
        let s = _mm256_set1_epi32(s_off as i32);
        let v = _mm256_blend_epi32::<0xAA>(_mm256_permutevar8x32_epi32(c, order), s);
        _mm256_storeu_si256(d as *mut __m256i, v);
    } else {
        let cursor = (*cell).entries[0] as usize;
        let src = table.overflow.add(cursor);
        let high = _mm256_set1_epi64x((s_off as i64) << 32);
        let whole = chain_whole(table, cursor, numhits);
        let mut j = 0;
        while j < whole {
            let q = _mm256_cvtepu32_epi64(_mm_loadu_si128(src.add(j) as *const __m128i));
            _mm256_storeu_si256(d.add(j) as *mut __m256i, _mm256_or_si256(q, high));
            j += 4;
        }
        while j < numhits {
            *d.add(j) = BlastOffsetPair {
                q_off: *src.add(j) as u32,
                s_off,
            };
            j += 1;
        }
    }
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
    /// each call; `avx2` = `None` for NCBI's loop. `RESUMED` counts the calls that resume the
    /// block held by the call before: a block is held and was kept for this call's subject and
    /// `s_range` (`resume`: one `XScanHeld` for the range, as `LOSAT_X_TBNSCANRESUME` keeps it;
    /// else a new one per call).
    fn calls(
        lookup: &BlastAaLookupTable,
        presence: &[u8],
        subject: &[u8],
        range: (i32, i32),
        array_size: usize,
        avx2: Option<bool>,
    ) -> Vec<(Vec<(u32, u32)>, i32)> {
        calls_held(lookup, presence, subject, range, array_size, avx2, false)
    }

    thread_local! {
        static RESUMED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    }

    /// `calls`, with one `XScanHeld` for the range when `resume`.
    fn calls_held(
        lookup: &BlastAaLookupTable,
        presence: &[u8],
        subject: &[u8],
        range: (i32, i32),
        array_size: usize,
        avx2: Option<bool>,
        resume: bool,
    ) -> Vec<(Vec<(u32, u32)>, i32)> {
        let mut held = XScanHeld::new();
        // Guard entries after the slack: no form may write past `array_size + X_SCAN_BYTE_SLACK`.
        let guard = BlastOffsetPair {
            q_off: 0xdead_beef,
            s_off: 0xfeed_f00d,
        };
        let end = array_size + X_SCAN_BYTE_SLACK;
        let mut pairs = vec![BlastOffsetPair::default(); end];
        pairs.resize(end + 2 * X_SCAN_BYTE_SLACK, guard);
        let mut scan_range = [0, range.0, range.1];
        let mut out = Vec::new();
        while scan_range[1] <= scan_range[2] {
            let hits = match avx2 {
                Some(avx2) => {
                    let at = (subject.as_ptr() as usize, subject.len());
                    if !resume {
                        held = XScanHeld::new();
                    } else if held.next < held.n
                        && held.at == (at.0, at.1, scan_range[1], scan_range[2])
                    {
                        RESUMED.with(|r| r.set(r.get() + 1));
                    }
                    scan(
                        avx2,
                        lookup,
                        presence,
                        subject,
                        &mut pairs,
                        array_size as i32,
                        &mut scan_range,
                        &mut held,
                    )
                }
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
            assert!(
                pairs[end..]
                    .iter()
                    .all(|pair| (pair.q_off, pair.s_off) == (guard.q_off, guard.s_off)),
                "write past the slack (array_size {array_size}, avx2 {avx2:?})"
            );
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
        RESUMED.with(|r| r.set(0));
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
                    4096 + longest,
                    32768 + longest,
                    longest * 1024,
                ] {
                    let want = calls(&lookup, &presence, &subject, range, array_size, None);
                    for &avx2 in &forms {
                        for resume in [false, true] {
                            let got = calls_held(
                                &lookup,
                                &presence,
                                &subject,
                                range,
                                array_size,
                                Some(avx2),
                                resume,
                            );
                            assert_eq!(
                                got, want,
                                "case {case} range {range:?} array_size {array_size} avx2 {avx2} \
                                 resume {resume}"
                            );
                        }
                    }
                    compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                    stops += want.len().saturating_sub(1);
                }
            }
        }
        // The cases above must reach the capacity stop (each one resumed in its block with
        // `resume`), many hits and cells of the overflow.
        assert!(compared > 100_000, "{compared} pairs compared");
        assert!(stops > 100, "{stops} capacity stops");
        assert_eq!(
            RESUMED.with(|r| r.get()),
            stops * forms.len(),
            "every stop resumed in its block"
        );
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

    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:99-115 (the copy of an overflow chain).
    // Not a port of an NCBI unit case: the chain that ends the overflow array, of a length that is
    // not a multiple of four, is copied exactly (the AVX2 form then writes its last pairs one by
    // one instead of reading past the array), with offset arrays from the longest chain up (as
    // `x_scan_unit`'s: a smaller one would never hold the longest cell, and NCBI's loop would stop
    // at that word on every call).
    #[test]
    fn the_last_overflow_chain_is_copied_exactly() {
        let mut rng = XorShift(0x5851_f42d_4c95_7f2d);
        let mut tails = 0usize;
        for case in 0..64 {
            let queries: Vec<Vec<u8>> = (0..1 + case % 5)
                .map(|_| {
                    let len = 30 + rng.below(600);
                    (0..len).map(|_| AMINO[rng.below(AMINO.len())]).collect()
                })
                .collect();
            let lookup = lookup_for(&queries, [11, 12, 13][case % 3]);
            let Some((ix, cell)) = lookup
                .backbone
                .iter()
                .enumerate()
                .filter(|(_, cell)| cell.num_used as usize > AA_HITS_PER_CELL)
                .max_by_key(|(_, cell)| cell.entries[0])
            else {
                continue;
            };
            let hits = cell.num_used as usize;
            let cursor = cell.entries[0] as usize;
            if cursor + hits != lookup.overflow.len() || hits % 4 == 0 {
                continue;
            }
            tails += 1;
            let presence = build_presence(&lookup);
            // The chain's word three times, then the sentinel.
            let word = [(ix >> 10) as u8, ((ix >> 5) & 31) as u8, (ix & 31) as u8];
            let mut subject: Vec<u8> = word.repeat(3);
            subject.push(0);
            let longest = lookup.longest_chain.max(1) as usize;
            for array_size in [
                longest,
                longest + 1,
                longest + hits,
                3 * longest,
                4096 + longest,
                32768 + longest,
                longest * 1024,
            ] {
                let range = (0, subject.len() as i32 - 4);
                let want = calls(&lookup, &presence, &subject, range, array_size, None);
                assert!(
                    want[0].0.len() >= hits,
                    "case {case}: the chain was not copied"
                );
                let mut forms = vec![false];
                if avx2_detected() {
                    forms.push(true);
                }
                for avx2 in forms {
                    for resume in [false, true] {
                        let got = calls_held(
                            &lookup,
                            &presence,
                            &subject,
                            range,
                            array_size,
                            Some(avx2),
                            resume,
                        );
                        assert_eq!(
                            got, want,
                            "case {case} array_size {array_size} avx2 {avx2} resume {resume}"
                        );
                    }
                }
            }
        }
        assert!(tails >= 4, "{tails} tables whose last chain ends the array");
    }

    // No NCBI counterpart: the kept presence bytes follow the table (a new table, also one at a
    // freed table's address, gets its own).
    #[test]
    fn presence_bytes_are_rebuilt_for_another_table() {
        let one = Arc::new(lookup_for(&[b"MKTAYIAKQRQISFVKSHFSRQ".to_vec()], 11));
        let two = Arc::new(lookup_for(&[b"WWWWCCCCHHHHMMMM".to_vec()], 11));
        assert_eq!(*presence_bytes(&one).0, *build_presence(&one));
        assert_eq!(*presence_bytes(&two).0, *build_presence(&two));
        assert_eq!(*presence_bytes(&one).0, *build_presence(&one));
        drop(one);
        let three = Arc::new(lookup_for(&[b"PPPPGGGGAAAA".to_vec()], 11));
        assert_eq!(*presence_bytes(&three).0, *build_presence(&three));
    }

    // No NCBI counterpart: a table is sparse when at most three tenths of its backbone cells have
    // hits; the table of one random protein of 400 letters is (0.21 at threshold 13), that of 64 is
    // not (0.54), and `XScanByte::Sparse` scans only the sparse one this way.
    #[test]
    fn one_protein_makes_a_sparse_table_and_many_do_not() {
        let mut rng = XorShift(0x0bad_5eed_1234_5678);
        let mut protein =
            || -> Vec<u8> { (0..400).map(|_| AMINO[rng.below(AMINO.len())]).collect() };
        let one = Arc::new(lookup_for(&[protein()], 13));
        let many = Arc::new(lookup_for(
            &(0..64).map(|_| protein()).collect::<Vec<_>>(),
            11,
        ));
        for (table, want) in [(&one, true), (&many, false)] {
            let bytes = build_presence(table);
            let ones = bytes.iter().filter(|&&b| b == 1).count();
            let cells = table.backbone.len();
            assert_eq!(sparse(&bytes, table), ones * 10 <= cells * X_SPARSE_TENTHS);
            assert_eq!(presence_bytes(table).1, want, "{ones} of {cells} cells");
            assert_eq!(x_presence(XScanByte::Sparse, table).is_some(), want);
            let mode = XScanByte::On { avx2: false };
            assert!(x_presence(mode, table).is_some(), "{mode:?}");
            assert!(x_presence(XScanByte::Off, table).is_none());
        }
    }
}
