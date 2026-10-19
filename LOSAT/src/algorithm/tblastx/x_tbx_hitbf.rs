//! EXPERIMENT (LOSAT_X_TBXHITBF / LOSAT_X_TBXHITBFSHADOW): the hit loop of a unit of the unit path
//! (x_pair_par.rs `scan_unit`, x_tbx_extbatch.rs `scan_unit_batched`) with the too-far and overlap
//! outcomes of the two-hit test taken without a branch. A child module of x_pair_par.rs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-552
//! ```c
//! if (diag_array[diag_coord].flag) {
//! ...
//! }
//! else {
//!     last_hit = diag_array[diag_coord].last_hit - diag_offset;
//!     diff = subject_offset - last_hit;
//!
//!     if (diff >= window) {
//!         diag_array[diag_coord].last_hit =
//!             subject_offset + diag_offset;
//!         continue;
//!     }
//!
//!     if (diff < wordsize) {
//!         continue;
//!     }
//! ```
//! In the reference loop most hits end in one of these two `continue`s (on AP027131 x AP027133,
//! 80 % too far and 16 % overlap), and which one is taken depends on the data: the too-far test
//! was mispredicted 0.14 times per hit there (v4-R6 REPORT, workstation). `x_hits_bf` computes both tests from one
//! subtraction and writes the cell's next value, chosen without a branch: `s + window` (too far) or
//! the value it read (overlap: the same value back). A flag-1 cell or a hit that passes both tests
//! (a candidate) goes to the per-hit body of the scan (`x_hit`, or `x_hit_batched`), one branch that
//! is rarely taken. The loop walks the batch by pointer, and the prefetch of
//! LOSAT_X_TBXPREFETCH (`PF`) reads the entry `X_PREFETCH_AHEAD` ahead without a bound test (the
//! offset-pair vector has that many entries more; `scan_unit_hits`).
//!
//! Why the result is that of the reference loop. Let L be a cell's `last_hit` (the 31-bit field,
//! sign-extended) and s the hit's subject offset. A unit's table starts at raw 0 with
//! `diag_offset == window` (x_pair_par.rs module doc, item 1), so in the flag-0 branch
//! `diff = s - L + window`, and in exact integers
//! - too far: `diff >= window` <=> `s - L >= 0`;
//! - overlap: `diff < wordsize` <=> `s - L < wordsize - window`;
//! - otherwise `wordsize - window <= s - L < 0`, the candidate.
//!
//! The rule `fits` (checked per unit): `1 <= wordsize < window` and chunk length + window < 2^29.
//! Every s is below the chunk length. Every value a unit stores with flag 0 is 0 (the fresh cell)
//! or `s + window` (too far, the flag reset, the context check, an extension without a right
//! extension), so a flag-0 cell holds L in [0, 2^29) in its raw word (bit 31 clear), with nothing
//! cut by the 31-bit field. The values stored with flag 1 (a right extension; LOSAT_X_TBXEXTBATCH's
//! pending marks) are never decoded here. So for a flag-0 cell `d2 = 2s - 2 raw` (wrapping i32) is
//! exactly 2 (s - L), with |2 (s - L)| < 2^30, and with `ow2 = 2 (wordsize - window)`, which lies in
//! (-2^30, 0):
//! - too far <=> `d2 >= 0`;
//! - overlap <=> `d2 < ow2` (both sides even);
//! - candidate <=> `ow2 <= d2 < 0` <=> `(d2 - ow2) as u32 < -ow2`: d2 - ow2 is exact in i32 (it
//!   lies in (-2^30, 2^31)), and outside the range it is negative (2^31 or more as u32) or at least
//!   -ow2.
//!
//! The reference's too-far write stores `s + diag_offset` truncated to 31 bits with the flag kept
//! (0); `s + window` is below 2^29, so the raw word `s + window` is that value. The overlap path
//! writes back the word it read. A flag-1 cell or a candidate goes to the per-hit body unchanged,
//! which then reads and writes the cell as in the reference. Each hit reads and writes only its own
//! cell, in scan order, so every cell goes through the reference's states and every call of the
//! per-hit body sees the cell the reference loop gives it; the init HSPs are saved by the same calls
//! in the same order. A unit outside `fits` runs the reference loop, which keeps the 31-bit field's
//! truncation (memory `losat-no-ncbi-bug-reproduction`, 2026-10-12); the batch guard `last_hit_fits`
//! already keeps the units of the unit path inside it.
//!
//! The rule: the switch is in the strict set (on unless `LOSAT_X_TBXHITBF` is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set); the probes of LOSAT_TIMING / LOSAT_X_STATS run the
//! reference loop (they count every outcome). `LOSAT_X_TBXHITBFSHADOW=1` (turns the path on) runs
//! the reference loop and this one on each unit inside `fits`, each on a fresh table, and stops at
//! the first difference of the saved records or of the final table.

use super::x_tbx_extbatch::{XShadowList, XShadowRef};
use super::*;

// No NCBI counterpart: counters of LOSAT_X_TBXHITBFSHADOW (units compared, records compared) and of
// the units outside `fits` while the switch is on; they do not change any value NCBI computes.
static HB_SHADOW_UNITS: AtomicU64 = AtomicU64::new(0);
static HB_SHADOW_RECORDS: AtomicU64 = AtomicU64::new(0);
static HB_OUTSIDE: AtomicU64 = AtomicU64::new(0);

/// `fits`: chunk length + window stays below this (2^29; module doc).
const X_HITBF_LIMIT: i64 = 1 << 29;

// No NCBI counterpart: which hit loop the units of the unit path use (scheduling only). Read once
// per process, as the other modes of x_pair_par.rs.
/// 0 = off, 1 = LOSAT_X_TBXHITBF (a strict-set switch, `utils/xswitch.rs`), 2 =
/// LOSAT_X_TBXHITBFSHADOW (set to any value; turns the path on).
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_TBXHITBFSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_TBXHITBF") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: the range in which `x_hits_bf` takes the reference's outcomes (module doc).
/// Whether a unit with a subject chunk of `chunk_length` letters may run `x_hits_bf`:
/// `1 <= wordsize < window` and `chunk_length + window < 2^29`.
pub(super) fn fits(chunk_length: usize, window: i32, wordsize: i32) -> bool {
    wordsize >= 1
        && window > wordsize
        && chunk_length < X_HITBF_LIMIT as usize
        && (chunk_length as i64) + (window as i64) < X_HITBF_LIMIT
}

// No NCBI counterpart: the switch's mode `m` for one unit (module doc).
/// The mode a unit runs: `m`, or 0 (the reference loop) when the unit is outside `fits`.
fn unit_mode(m: u8, chunk_length: usize, window: i32, wordsize: i32) -> u8 {
    if m == 0 {
        return 0;
    }
    if !fits(chunk_length, window, wordsize) {
        HB_OUTSIDE.fetch_add(1, AtomicOrdering::Relaxed);
        return 0;
    }
    m
}

/// The summary line of the shadow.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] TBXHITBF units_compared={} records_compared={} outside={} (records and unit tables identical)",
            HB_SHADOW_UNITS.load(AtomicOrdering::Relaxed),
            HB_SHADOW_RECORDS.load(AtomicOrdering::Relaxed),
            HB_OUTSIDE.load(AtomicOrdering::Relaxed)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:510-552 (module doc)
// ```c
// for (i = 0; i < hits; ++i) {
//     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ...
//     diag_coord = (query_offset - subject_offset) & diag_mask;
// ```
// The hit loop of `scan_unit` over one batch: the cell of each hit as there (the raw diagonal minus
// the unit's lowest one), the too-far and overlap outcomes of the flag-0 branch from one
// subtraction with the cell's next value selected (module doc), and `rare` for a flag-1 cell or a
// candidate. `ow2` and `cand_len` are locals so that they stay in registers.
/// The two-hit test of the `n` hits at `pairs` on the unit table at `diag_ptr` (cell of a hit
/// (q, s): q - s - lo), with `rare(cell, q, s)` for a hit whose cell has flag 1 or that passes the
/// test; `PF` prefetches the cell of the hit `X_PREFETCH_AHEAD` entries ahead (x86-64).
///
/// SAFETY: `pairs` points to `n` entries, and with `PF` to `X_PREFETCH_AHEAD` initialized entries
/// after them; every hit's cell is below `width`, the length of the table `diag_ptr` points into;
/// `fits(chunk length, window, wordsize)` holds, every subject offset is below the chunk length and
/// every flag-0 cell holds a value the unit stored (module doc); `rare` writes at most the cell it
/// is given.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) unsafe fn x_hits_bf<const PF: bool>(
    diag_ptr: *mut DiagStruct,
    pairs: *const OffsetPair,
    n: usize,
    lo: u32,
    width: usize,
    window: i32,
    wordsize: i32,
    mut rare: impl FnMut(usize, u32, u32),
) {
    let window_u = window as u32;
    // 2 (wordsize - window), below 0, and the length of the candidate range [ow2, 0).
    let ow2 = (wordsize - window).wrapping_mul(2);
    let cand_len = ow2.wrapping_neg() as u32;
    let mut p = pairs;
    // SAFETY: the caller's bound (n entries).
    let end = unsafe { pairs.add(n) };
    while p < end {
        // EXPERIMENT (LOSAT_X_TBXPREFETCH): the cell of the hit X_PREFETCH_AHEAD entries ahead,
        // computed as below; past the batch the entry is a stale one and the address is only
        // prefetched (no architectural effect, `wrapping_add`).
        #[cfg(target_arch = "x86_64")]
        if PF {
            // SAFETY: the caller's bound (X_PREFETCH_AHEAD initialized entries after the n).
            let ahead = unsafe { p.add(X_PREFETCH_AHEAD).read() };
            let cell = ahead.q_off.wrapping_sub(ahead.s_off).wrapping_sub(lo) as usize;
            // SAFETY: a prefetch never faults.
            unsafe {
                core::arch::x86_64::_mm_prefetch(
                    diag_ptr.wrapping_add(cell) as *const i8,
                    core::arch::x86_64::_MM_HINT_T0,
                )
            };
        }
        // SAFETY: p < end.
        let OffsetPair {
            q_off: query_offset,
            s_off: subject_offset,
        } = unsafe { p.read() };
        // SAFETY: p < end, so p + 1 <= end.
        p = unsafe { p.add(1) };
        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
        // Unit cell: the raw diagonal minus the unit's lowest one (no mask).
        let cell = query_offset.wrapping_sub(subject_offset).wrapping_sub(lo) as usize;
        debug_assert!(cell < width, "unit diagonal outside its table");
        // SAFETY: cell < width (the caller's bound).
        let entry = unsafe { diag_ptr.add(cell) };
        // SAFETY: as above.
        let raw = unsafe { (*entry).raw_bits() };
        // 2 (s - L) for a flag-0 cell (module doc); any value for a flag-1 cell, which goes to
        // `rare` below.
        let d2 = (subject_offset.wrapping_add(subject_offset) as i32)
            .wrapping_sub(raw.wrapping_add(raw) as i32);
        let candidate = (d2.wrapping_sub(ow2) as u32) < cand_len;
        // [C] if (diag_array[diag_coord].flag), or a hit past both tests.
        if ((raw as i32) < 0) | candidate {
            rare(cell, query_offset, subject_offset);
            continue;
        }
        // [C] if (diff >= window) last_hit = subject_offset + diag_offset;
        // [C] if (diff < wordsize) continue;
        let next = if d2 >= 0 {
            subject_offset.wrapping_add(window_u)
        } else {
            raw
        };
        // SAFETY: as above.
        unsafe { (*entry).set_raw_bits(next) };
    }
}

// No NCBI counterpart: compares two LOSAT scans of one unit (LOSAT_X_TBXHITBFSHADOW).
/// Panics with the unit at the first difference of the records or of the tables.
fn shadow_compare<S: XShadowList>(
    scan: &str,
    inp: &XSubject<'_>,
    unit: &XUnit,
    records: &S,
    reference: &S,
    diag: &[DiagStruct],
    reference_diag: &[DiagStruct],
) {
    if records != reference || diag != reference_diag {
        let cell = diag.iter().zip(reference_diag).position(|(a, b)| a != b);
        panic!(
            "LOSAT_X_TBXHITBFSHADOW: {scan}: subject {} frame {} chunk {} group {}: records {} vs {} (first difference at {}), first table difference at {:?}",
            inp.s_idx,
            unit.f,
            unit.k,
            unit.g,
            records.len_of(),
            reference.len_of(),
            records.first_difference(reference),
            cell
        );
    }
    HB_SHADOW_UNITS.fetch_add(1, AtomicOrdering::Relaxed);
    HB_SHADOW_RECORDS.fetch_add(records.len_of() as u64, AtomicOrdering::Relaxed);
}

impl XPairPar {
    // EXPERIMENT (LOSAT_X_TBXHITBF): dispatch of the units without probes (`run_unit`,
    // `run_unit_compact`) by the switch and `fits`.
    /// `scan_unit` without probes, or its branch-free hit loop (LOSAT_X_TBXHITBFSHADOW: both,
    /// compared).
    #[inline(always)]
    pub(super) fn scan_unit_or_bf<S: XShadowList, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        let m = unit_mode(mode(), unit.chunk.length, inp.window, inp.wordsize);
        self.scan_unit_mode::<S, PF>(m, unit, inp, scratch, init_hsps);
    }

    /// `scan_unit_or_bf` in mode `m` (0, 1 or 2 as `mode`; the caller has applied `fits`).
    #[inline(always)]
    fn scan_unit_mode<S: XShadowList, const PF: bool>(
        &self,
        m: u8,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        match m {
            0 => self.scan_unit::<S, XNoProbe, PF>(unit, inp, scratch, init_hsps, &mut XNoProbe),
            1 => self.scan_unit_bf::<S, PF>(unit, inp, scratch, init_hsps),
            _ => self.scan_unit_bf_shadow::<S, PF>(unit, inp, scratch, init_hsps),
        }
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): its own function, so that the loop gets its own registers and
    // the units outside the switch keep the code of `scan_unit`.
    /// `scan_unit` with the branch-free hit loop.
    #[inline(never)]
    fn scan_unit_bf<S: XInitSink, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        self.scan_unit_hits::<S, XNoProbe, PF, true>(unit, inp, scratch, init_hsps, &mut XNoProbe);
    }

    // EXPERIMENT (LOSAT_X_TBXHITBFSHADOW): the reference loop (with its own probe type, so that the
    // shadow adds no caller to the instance the units run) on a fresh table into its own list, then
    // the branch-free loop, compared.
    /// The branch-free scan of one unit, compared with `scan_unit`.
    #[cold]
    #[inline(never)]
    fn scan_unit_bf_shadow<S: XShadowList, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        let mut reference = S::default();
        self.scan_unit::<S, XShadowRef, PF>(unit, inp, scratch, &mut reference, &mut XShadowRef);
        let reference_diag = scratch.diag.clone();
        self.scan_unit_bf::<S, PF>(unit, inp, scratch, init_hsps);
        shadow_compare(
            "scan_unit",
            inp,
            unit,
            init_hsps,
            &reference,
            &scratch.diag,
            &reference_diag,
        );
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): dispatch of the batched scan (LOSAT_X_TBXEXTBATCH) by the
    // switch and `fits`.
    /// `scan_unit_batched` with full batches, or its branch-free hit loop
    /// (LOSAT_X_TBXHITBFSHADOW: both, compared).
    #[inline(always)]
    pub(super) fn scan_unit_batched_or_bf<S: XShadowList, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        let m = unit_mode(mode(), unit.chunk.length, inp.window, inp.wordsize);
        self.scan_unit_batched_mode::<S, PF>(m, unit, inp, scratch, init_hsps);
    }

    /// `scan_unit_batched_or_bf` in mode `m` (as `scan_unit_mode`).
    #[inline(always)]
    fn scan_unit_batched_mode<S: XShadowList, const PF: bool>(
        &self,
        m: u8,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        match m {
            0 => {
                self.scan_unit_batched::<S, PF>(
                    unit,
                    inp,
                    scratch,
                    init_hsps,
                    x_tbx_extbatch::BATCH,
                );
            }
            1 => self.scan_unit_batched_bf::<S, PF>(unit, inp, scratch, init_hsps),
            _ => self.scan_unit_batched_bf_shadow::<S, PF>(unit, inp, scratch, init_hsps),
        }
    }

    // EXPERIMENT (LOSAT_X_TBXHITBF): its own function, as `scan_unit_bf`.
    /// `scan_unit_batched` with full batches and the branch-free hit loop.
    #[inline(never)]
    fn scan_unit_batched_bf<S: XInitSink, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        self.scan_unit_batched_hits::<S, PF, true>(
            unit,
            inp,
            scratch,
            init_hsps,
            x_tbx_extbatch::BATCH,
        );
    }

    // EXPERIMENT (LOSAT_X_TBXHITBFSHADOW): the batched scan with the loop before the switch, then
    // with the branch-free loop, each on a fresh table, compared.
    /// The branch-free batched scan of one unit, compared with `scan_unit_batched`.
    #[cold]
    #[inline(never)]
    fn scan_unit_batched_bf_shadow<S: XShadowList, const PF: bool>(
        &self,
        unit: &XUnit,
        inp: &XSubject<'_>,
        scratch: &mut XScratch,
        init_hsps: &mut S,
    ) {
        let mut reference = S::default();
        self.scan_unit_batched::<S, PF>(unit, inp, scratch, &mut reference, x_tbx_extbatch::BATCH);
        let reference_diag = scratch.diag.clone();
        self.scan_unit_batched_bf::<S, PF>(unit, inp, scratch, init_hsps);
        shadow_compare(
            "scan_unit_batched",
            inp,
            unit,
            init_hsps,
            &reference,
            &scratch.diag,
            &reference_diag,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A frame of amino acids (NCBISTDAA codes) with its sentinels.
    fn frame(frame: i8, residues: &[u8]) -> QueryFrame {
        let mut aa_seq = Vec::with_capacity(residues.len() + 2);
        aa_seq.push(0);
        aa_seq.extend_from_slice(residues);
        aa_seq.push(0);
        QueryFrame {
            frame,
            aa_seq,
            aa_seq_nomask: None,
            aa_len: residues.len(),
            orig_len: residues.len() * 3 + 2,
            seg_masks: Vec::new(),
        }
    }

    /// xorshift64*, for reproducible random values.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Residues from `alphabet`, with stretches copied from `source` (shared diagonals: long
    /// extensions, right extensions and later hits on their cells).
    fn residues(rng: &mut Rng, len: usize, alphabet: &[u8], source: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(len);
        while out.len() < len {
            if !source.is_empty() && rng.below(3) == 0 {
                let at = rng.below(source.len() as u64) as usize;
                let take = (rng.below(60) as usize + 5).min(source.len() - at);
                out.extend_from_slice(&source[at..at + take]);
            } else {
                out.push(alphabet[rng.below(alphabet.len() as u64) as usize]);
            }
        }
        out.truncate(len);
        out
    }

    /// A probe that counts the outcomes of the reference loop's hits.
    #[derive(Default)]
    struct XCount([u64; 7]);

    impl XProbe for XCount {
        type Mark = ();
        fn unit_begin(&mut self, _width: usize) {}
        fn batch(&mut self, _hits: usize) {}
        fn at(&mut self, _cell: usize, _s: u32) {}
        fn hit(&mut self, outcome: XHit, _window: i32) {
            self.0[outcome as usize] += 1;
        }
        fn batch_end(&mut self) {}
        fn ext_begin(&mut self) {}
        fn ext_end(
            &mut self,
            _mark: (),
            _q: &[u8],
            _s: &[u8],
            _s_left: usize,
            _s_right: usize,
            _q_right: usize,
            _x_drop: i32,
        ) {
        }
        fn unit_end(&mut self) {}
    }

    /// What the tests' `rare` does to a cell (the same in both loops): a flag-1 cell is masked
    /// below its `last_hit` and reset above it, a candidate sets the flag with `last_hit = s +
    /// window` on odd s and stores `s + window` on even s, as `x_hit` stores one of those.
    fn rare_apply(cell: &mut DiagStruct, s: u32, window: i32) {
        let next = s.wrapping_add(window as u32) as i32;
        if cell.flag() != 0 {
            if (next as u32) >= cell.last_hit() as u32 {
                cell.set_last_hit(next);
                cell.set_flag(0);
            }
        } else {
            cell.set_last_hit(next);
            cell.set_flag(s & 1);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:518-552 (module doc)
    // `x_hits_bf` against the reference's statements on random tables and hits over the range of
    // `fits`, with the last_hit values on both sides of each boundary of the two tests.
    #[test]
    fn x_hits_bf_takes_the_reference_outcomes() {
        let mut rng = Rng(0x0b1f_5eed_7777_0001);
        // too far, overlap, flag 1, candidate; s - L at 0, -1, wordsize - window, one below it.
        let mut taken = [0u64; 4];
        let mut edges = [0u64; 4];
        for round in 0..3000 {
            let wordsize = [2, 3, 3, 5][round % 4];
            let window = match round % 5 {
                0 => wordsize + 1,
                1 | 4 => 40,
                2 => 64,
                _ => 1 << 20,
            };
            // Every seventh round the longest chunk `fits` takes.
            let chunk_len: usize = if round % 7 == 0 {
                ((X_HITBF_LIMIT - 1) - window as i64) as usize
            } else {
                1 + rng.below(5000) as usize
            };
            assert!(fits(chunk_len, window, wordsize));
            assert!(!fits(chunk_len + 1, window, wordsize) || round % 7 != 0);
            let width = 1 + rng.below(48) as usize;
            let lo = rng.next() as u32;
            let n = rng.below(400) as usize;
            let mut table = vec![DiagStruct::default(); width];
            let mut touched = vec![false; width];
            let mut pairs: Vec<OffsetPair> = Vec::with_capacity(n + X_PREFETCH_AHEAD);
            for _ in 0..n {
                let cell = rng.below(width as u64) as usize;
                let s = match rng.below(8) {
                    0 => 0,
                    1 => chunk_len as u32 - 1,
                    _ => rng.below(chunk_len as u64) as u32,
                };
                if !touched[cell] {
                    touched[cell] = true;
                    // The first hit on a cell: a start value on an edge of a test, or any value a
                    // unit stores.
                    let delta: i64 = match rng.below(9) {
                        0 => 0,
                        1 => -1,
                        2 => 1,
                        3 => (wordsize - window) as i64,
                        4 => (wordsize - window) as i64 - 1,
                        5 => (wordsize - window) as i64 + 1,
                        _ => i64::MAX,
                    };
                    let last_hit = s as i64 - delta;
                    let c = &mut table[cell];
                    if delta != i64::MAX && last_hit >= 0 {
                        c.set_raw_bits(last_hit as u32);
                    } else if rng.below(4) == 0 {
                        c.set_raw_bits(0);
                    } else {
                        let v = rng.below(chunk_len as u64 + window as u64) as u32;
                        c.set_raw_bits(v | if rng.below(5) == 0 { 0x8000_0000 } else { 0 });
                    }
                }
                let q = (cell as u32).wrapping_add(s).wrapping_add(lo);
                pairs.push(OffsetPair { q_off: q, s_off: s });
            }
            // Stale entries after the batch, as the vector keeps them.
            for _ in 0..X_PREFETCH_AHEAD {
                pairs.push(OffsetPair {
                    q_off: rng.next() as u32,
                    s_off: rng.next() as u32,
                });
            }

            // The reference's statements (diag_offset = window).
            let mut reference = table.clone();
            let mut reference_log = Vec::new();
            for pair in &pairs[..n] {
                let (q, s) = (pair.q_off, pair.s_off);
                let cell = q.wrapping_sub(s).wrapping_sub(lo) as usize;
                let c = &mut reference[cell];
                if c.flag() != 0 {
                    taken[2] += 1;
                    reference_log.push((cell, q, s));
                    rare_apply(c, s, window);
                    continue;
                }
                let d = s as i64 - c.last_hit() as i64;
                for (k, edge) in [
                    0,
                    -1,
                    (wordsize - window) as i64,
                    (wordsize - window) as i64 - 1,
                ]
                .into_iter()
                .enumerate()
                {
                    edges[k] += u64::from(d == edge);
                }
                let last_hit = c.last_hit() - window;
                let diff = s.wrapping_sub(last_hit as u32) as i32;
                if diff >= window {
                    taken[0] += 1;
                    c.set_last_hit(s.wrapping_add(window as u32) as i32);
                    continue;
                }
                if diff < wordsize {
                    taken[1] += 1;
                    continue;
                }
                taken[3] += 1;
                reference_log.push((cell, q, s));
                rare_apply(c, s, window);
            }

            for pf in [false, true] {
                let mut got = table.clone();
                let mut log = Vec::new();
                let ptr = got.as_mut_ptr();
                let rare = |cell: usize, q: u32, s: u32| {
                    log.push((cell, q, s));
                    // SAFETY: cell < width, the table's length (the loop's bound).
                    rare_apply(unsafe { &mut *ptr.add(cell) }, s, window);
                };
                // SAFETY: n entries and X_PREFETCH_AHEAD more; every cell is below width; `fits`
                // holds and every flag-0 value is below chunk_len + window.
                unsafe {
                    if pf {
                        x_hits_bf::<true>(ptr, pairs.as_ptr(), n, lo, width, window, wordsize, rare)
                    } else {
                        x_hits_bf::<false>(
                            ptr,
                            pairs.as_ptr(),
                            n,
                            lo,
                            width,
                            window,
                            wordsize,
                            rare,
                        )
                    }
                };
                assert_eq!(log, reference_log, "round {round} pf {pf}: rare calls");
                assert_eq!(got, reference, "round {round} pf {pf}: table");
            }
        }
        assert!(taken.iter().all(|&t| t > 1000), "outcomes {taken:?}");
        assert!(edges.iter().all(|&e| e > 100), "edges {edges:?}");
    }

    // No NCBI counterpart: compares the branch-free scans with the scans before the switch on the
    // same units (module doc), and the batched ones (LOSAT_X_TBXEXTBATCH) with batches of 1, 3 and
    // 16 candidates.
    #[test]
    fn bf_scans_equal_the_reference_scans() {
        let alphabet: Vec<u8> = "ACDEKLMNPQRSTVWY"
            .bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect();
        let narrow: Vec<u8> = "AKLE"
            .bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let mut rng = Rng(0x41b0_f00d_2026_1019);
        let mut outcomes = [0u64; 7];
        let mut batched_units = 0usize;
        for case in 0..16 {
            let abc = if case % 3 == 2 { &narrow } else { &alphabet };
            let q_len = 150 + rng.below(350) as usize;
            let frames: Vec<QueryFrame> = [1i8, 2, 3, -1, -2, -3]
                .iter()
                .map(|&f| frame(f, &residues(&mut rng, q_len, abc, &[])))
                .collect();
            let source: Vec<u8> = frames
                .iter()
                .flat_map(|f| f.aa_seq[1..f.aa_len + 1].to_vec())
                .collect();
            let s_len = 300 + rng.below(500) as usize;
            // Every fourth case scans the query against itself (long extensions, flag-1 cells hit
            // again inside and past their extensions).
            let s_frames: Vec<QueryFrame> = if case % 4 == 3 {
                frames.clone()
            } else {
                [1i8, 2, 3, -1, -2, -3]
                    .iter()
                    .map(|&f| frame(f, &residues(&mut rng, s_len, abc, &source)))
                    .collect()
            };
            let (lookup, contexts) = build_ncbi_lookup(&[frames.clone()], 11, &params, true);
            let cutoff_scores: Vec<i32> = (0..contexts.len())
                .map(|i| 10 + (i as i32 % 4) * 3)
                .collect();
            let x_dropoff_per_context: Vec<i32> = (0..contexts.len())
                .map(|i| [7, 16, 30, 200][i % 4])
                .collect();
            let pp = XPairPar {
                groups: build_groups(&lookup, &contexts),
                shadow: false,
                stub: false,
            };
            // The default window, the smallest `fits` takes, a wide one, and the smallest of the
            // batched scan.
            let window = [40, 4, 64, 3 + x_tbx_extbatch::BATCH as i32][case % 4];
            let inp = XSubject {
                s_frames_preliminary: &s_frames,
                s_frames: &s_frames,
                contexts: &contexts,
                cutoff_scores: &cutoff_scores,
                x_dropoff_per_context: &x_dropoff_per_context,
                window,
                wordsize: 3,
                s_idx: 0,
                s_len: 999,
            };
            let mut scratch = XScratch::new(false);
            for (f, s_frame) in s_frames.iter().enumerate() {
                // Every other case scans a chunk that starts inside the frame (offset > 0).
                let offset = if case % 2 == 1 { s_frame.aa_len / 3 } else { 0 };
                let chunk = SubjectChunk {
                    offset,
                    length: s_frame.aa_len - offset,
                    overlap: 0,
                };
                assert!(fits(chunk.length, window, 3));
                for (g, group) in pp.groups.iter().enumerate() {
                    if group.lookup.longest_chain <= 0 {
                        continue;
                    }
                    let unit = XUnit { f, k: 0, chunk, g };
                    let mut count = XCount::default();
                    let mut reference: Vec<InitHSP> = Vec::new();
                    pp.scan_unit::<_, _, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut reference,
                        &mut count,
                    );
                    let reference_diag = scratch.diag.clone();
                    for (o, c) in outcomes.iter_mut().zip(count.0) {
                        *o += c;
                    }
                    for pf in [false, true] {
                        let mut got: Vec<InitHSP> = Vec::new();
                        if pf {
                            pp.scan_unit_bf::<_, true>(&unit, &inp, &mut scratch, &mut got);
                        } else {
                            pp.scan_unit_bf::<_, false>(&unit, &inp, &mut scratch, &mut got);
                        }
                        assert_eq!(got, reference, "case {case} frame {f} group {g} pf {pf}");
                        assert!(
                            scratch.diag == reference_diag,
                            "case {case} frame {f} group {g} pf {pf}: table"
                        );
                    }
                    let mut compact: Vec<XUnitHsp> = Vec::new();
                    let mut compact_ref: Vec<XUnitHsp> = Vec::new();
                    pp.scan_unit::<_, _, false>(
                        &unit,
                        &inp,
                        &mut scratch,
                        &mut compact_ref,
                        &mut XNoProbe,
                    );
                    pp.scan_unit_bf::<_, false>(&unit, &inp, &mut scratch, &mut compact);
                    assert_eq!(
                        compact, compact_ref,
                        "case {case} frame {f} group {g}: compact"
                    );
                    // The batched scan's windows leave room for its pending marks.
                    if window < 3 + x_tbx_extbatch::BATCH as i32 {
                        continue;
                    }
                    batched_units += 1;
                    for len in [1, 3, 16] {
                        let mut got: Vec<InitHSP> = Vec::new();
                        if len == 3 {
                            pp.scan_unit_batched_hits::<_, true, true>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                                len,
                            );
                        } else {
                            pp.scan_unit_batched_hits::<_, false, true>(
                                &unit,
                                &inp,
                                &mut scratch,
                                &mut got,
                                len,
                            );
                        }
                        assert_eq!(
                            got, reference,
                            "case {case} frame {f} group {g} batch {len}"
                        );
                        assert!(
                            scratch.diag == reference_diag,
                            "case {case} frame {f} group {g} batch {len}: table"
                        );
                    }
                }
            }
        }
        // Every outcome of the reference loop but the context check (one query: the unit's frames
        // are contexts of one group only where they are apart) was reached.
        let [masked, reset, too_far, overlap, _context, extended, saved] = outcomes;
        assert!(
            masked > 100 && reset > 100 && too_far > 10_000 && overlap > 1000,
            "outcomes {outcomes:?}"
        );
        assert!(extended > 1000 && saved > 100, "outcomes {outcomes:?}");
        assert!(batched_units > 20, "batched units {batched_units}");
    }

    // No NCBI counterpart: the rule (module doc) at its boundaries.
    #[test]
    fn fits_at_its_boundaries_and_units_outside_run_the_reference() {
        let top = (X_HITBF_LIMIT - 1) as usize;
        assert!(fits(top - 40, 40, 3));
        assert!(!fits(top - 39, 40, 3));
        assert!(fits(0, 40, 3));
        assert!(fits(100, 4, 3));
        assert!(!fits(100, 3, 3));
        assert!(!fits(100, 0, 3));
        assert!(!fits(100, 40, 0));
        assert!(!fits(usize::MAX, 40, 3));
        assert!(!fits(100, i32::MAX, 3));
        for m in [0u8, 1, 2] {
            assert_eq!(unit_mode(m, 100, 40, 3), m);
            assert_eq!(unit_mode(m, top - 39, 40, 3), 0);
            assert_eq!(unit_mode(m, 100, 3, 3), 0);
        }
    }

    // No NCBI counterpart: the shadow's two scans of each unit (LOSAT_X_TBXHITBFSHADOW), for the
    // unit scan and the batched one.
    #[test]
    fn the_shadow_compares_both_scans() {
        let alphabet: Vec<u8> = "ACDEKLMNPQRSTVWY"
            .bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let mut rng = Rng(0x5a5a_0bf0_0000_0003);
        let frames: Vec<QueryFrame> = [1i8, 2, -1]
            .iter()
            .map(|&f| frame(f, &residues(&mut rng, 300, &alphabet, &[])))
            .collect();
        let s_frames = frames.clone();
        let (lookup, contexts) = build_ncbi_lookup(&[frames], 11, &params, true);
        let cutoff_scores = vec![12; contexts.len()];
        let x_dropoff_per_context = vec![16; contexts.len()];
        let pp = XPairPar {
            groups: build_groups(&lookup, &contexts),
            shadow: false,
            stub: false,
        };
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &contexts,
            cutoff_scores: &cutoff_scores,
            x_dropoff_per_context: &x_dropoff_per_context,
            window: 40,
            wordsize: 3,
            s_idx: 0,
            s_len: 999,
        };
        let mut scratch = XScratch::new(false);
        let before = HB_SHADOW_UNITS.load(AtomicOrdering::Relaxed);
        let mut units = 0u64;
        for (f, s_frame) in s_frames.iter().enumerate() {
            let chunk = SubjectChunk {
                offset: 0,
                length: s_frame.aa_len,
                overlap: 0,
            };
            for g in 0..pp.groups.len() {
                if pp.groups[g].lookup.longest_chain <= 0 {
                    continue;
                }
                let unit = XUnit { f, k: 0, chunk, g };
                let mut reference: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit::<_, _, false>(
                    &unit,
                    &inp,
                    &mut scratch,
                    &mut reference,
                    &mut XNoProbe,
                );
                let mut shadowed: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit_mode::<_, false>(2, &unit, &inp, &mut scratch, &mut shadowed);
                assert_eq!(shadowed, reference);
                let mut batched: Vec<XUnitHsp> = Vec::new();
                pp.scan_unit_batched_mode::<_, false>(2, &unit, &inp, &mut scratch, &mut batched);
                assert_eq!(batched, reference);
                units += 2;
            }
        }
        assert!(units > 0);
        assert!(HB_SHADOW_UNITS.load(AtomicOrdering::Relaxed) >= before + units);
    }

    // No NCBI counterpart: the shadow stops at a difference.
    #[test]
    #[should_panic(expected = "LOSAT_X_TBXHITBFSHADOW: scan_unit")]
    fn the_shadow_stops_at_a_difference() {
        let s_frames = vec![frame(1, &[1, 2, 3])];
        let inp = XSubject {
            s_frames_preliminary: &s_frames,
            s_frames: &s_frames,
            contexts: &[],
            cutoff_scores: &[],
            x_dropoff_per_context: &[],
            window: 40,
            wordsize: 3,
            s_idx: 0,
            s_len: 999,
        };
        let unit = XUnit {
            f: 0,
            k: 0,
            chunk: SubjectChunk {
                offset: 0,
                length: 3,
                overlap: 0,
            },
            g: 0,
        };
        let diag = vec![DiagStruct::default(); 4];
        let mut other = diag.clone();
        other[2].set_raw_bits(41);
        let records: Vec<XUnitHsp> = Vec::new();
        shadow_compare("scan_unit", &inp, &unit, &records, &records, &diag, &other);
    }
}
