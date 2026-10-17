//! EXPERIMENT (LOSAT_X_BPSEED): the blastp two-hit seed stage with the subject scan and the
//! diagonal test in one loop.
//!
//! The reference path copies every word hit of the subject into the offset-pair buffer
//! (`s_blast_aa_scan_subject_one_range`, NCBI's `s_BlastAaScanSubject`) and then reads the buffer
//! back for the two-hit test (NCBI's `s_BlastAaWordFinder_TwoHit`). Here the test reads the hit
//! list of the lookup cell in place, in the same order (subject offset ascending, the cell's list
//! in its stored order), and applies the outcomes that only touch the hit's diagonal entry. A hit
//! that needs the rest of the test (the context check and the ungapped extension) goes to the
//! caller's closure, in scan order, with the values the reference computes at that point.
//!
//! Why the outcome is the reference's: the test of a hit reads and writes only the entry of its
//! diagonal `(q_off - s_off) & diag_mask`, and the closure of a candidate writes only that entry
//! and the HSP list. The hits are taken in the reference's order, so every entry sees the same
//! sequence of reads and writes and the HSP list gets the same appends in the same order.
//! The buffer's capacity (`array_size`) only cuts the reference's scan into rounds; it changes
//! neither the order nor the set of hits.
//!
//! `LOSAT_X_BPSEED` is a switch of the strict set (`utils/xswitch.rs`): on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set.
//!
//! EXPERIMENT (LOSAT_X_BPDIAG): with the fused loop on, the hits of one subject offset (one lookup
//! cell) are tested in batches of eight instead of one after the other: for every hit of a batch
//! the entry's new value and whether the hit is a candidate are computed without branches in
//! AVX2 registers (eight lanes), every entry of the batch is stored, and then the candidates go to
//! the closure in the cell's order.
//! Why the outcome is the per-hit loop's:
//! - The query offsets of one cell are distinct: NCBI's `s_AddWordHits` adds a query offset to a
//!   cell once, as a neighbour of the offset's word or (when the word's self score is below the
//!   threshold, or the threshold is 0) as the word itself, which the neighbour enumeration then
//!   does not return. The diagonal table is at least the query's length plus the window
//!   (`BlastDiagTableNew`), so the hits of one subject offset fall on distinct entries.
//! - The test of a hit reads and writes only its entry, and the closure writes only the
//!   candidate's entry and the HSP list. So the hits of one subject offset are independent, and
//!   only the order of the closure's calls matters, which is the cell's order.
//! - A candidate's entry keeps its value in the test (the batch stores its old value back), so
//!   storing the batch before the calls changes nothing the closure reads, and a call's write is
//!   never overwritten by the batch.
//! - The values are the bit-field's (31-bit `last_hit`, sign-extended on read, truncated on
//!   write) and the i32 arithmetic of `x_test_hit`, which wraps (see there).
//!
//! `LOSAT_X_BPDIAG` is a switch of the strict set: on unless it is `0`, or `LOSAT_X_DEFAULT=0`
//! and it is not set. It needs `LOSAT_X_BPSEED` (or `LOSAT_X_BPSEEDSHADOW`, which then compares
//! the batched loop with the reference). Without AVX2
//! (another architecture, or an x86-64 CPU without it) the switch keeps the per-hit loop: the
//! scalar batch was slower than the per-hit loop on both machines it was timed on (2026-10-14,
//! human block 1 x Species4, one thread, scan_ungapped: workstation 12.56 s against 11.59 s,
//! local 7.53 / 7.08 s against 6.78 / 6.88 s). The scalar batch remains as the AVX2 path's handler
//! of a cell's last batch of fewer than eight hits and as the cell-by-cell version the tests
//! compare with; `LOSAT_X_BPDIAG_SIMD=0` forces it (diagnostic, used by the tests; not in the
//! strict set). With `LOSAT_TIMING=1` each search of a query set prints
//! `[TIMING] x_bpdiag: mode=<avx2|avx2-defer|scalar|per-hit> cells=<n> batched_hits=<m>
//! candidates=<k>` (the path that ran, the cells tested in batches, their hits, and the hits passed
//! to the closure; zeros with `per-hit`).
//!
//! EXPERIMENT (LOSAT_X_BPDEFER): with the AVX2 batches on, a batch's candidates are not called at
//! the end of the batch: their query offsets are packed into a buffer (the candidate lanes moved to
//! the front of the register by a permutation from a table of the 256 lane masks, the register
//! stored, the buffer's length advanced by the number of candidates), and the closure is called for
//! the buffer in order after the cell's last batch of eight. The batch loop then has no branch on
//! the candidates (on the workstation 60 % of the cell test's mispredicts were the per-batch "no
//! candidate" branch, 2026-10-16). Why the outcome is the per-batch calls': the hits of one cell lie
//! on distinct entries and a call writes only its candidate's entry (above), so the batches after a
//! candidate neither read nor write what its call reads or writes, and the calls keep their order
//! (the batches' order, then the lanes'); a candidate's entry still holds the value its test read.
//! A cell's last batch of fewer than eight hits is the scalar batch, called after the buffer.
//! `LOSAT_X_BPDEFER` is a switch of the strict set: on unless it is `0`, or `LOSAT_X_DEFAULT=0`
//! and it is not set. It acts only where the AVX2 batch runs (`LOSAT_X_BPDIAG` with AVX2, not
//! with `LOSAT_X_BPDIAG_SIMD=0`). `LOSAT_X_BPSEEDSHADOW` compares it with the reference like the
//! other modes.
//!
//! EXPERIMENT (LOSAT_X_BXSCANBF, blastx's switch): `x_scan_two_hit_with(.., branch_free: true, ..)`
//! finds the cells with hits by the candidate loop of LOSAT_X_TBNSCANBF
//! (`tblastn/x_scan_bf.rs`, `x_block_candidates`): the positions of a block of words are written
//! to a list whose length advances by the presence bit, without a branch on it, and the listed
//! cells then go to the cell test in position order. The cells, their order and their hit lists
//! are the presence-test loop's, so every outcome above is unchanged. blastp's own loop
//! (`x_scan_two_hit`) takes the same cells under its own switch, LOSAT_X_BPSCANBF (off unless set,
//! not in the strict set; [`x_bp_scan_bf`]).

use crate::algorithm::tblastn::x_scan_bf::{x_block_candidates, X_BLOCK};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::{
    BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use crate::utils::matrix::BLASTAA_SIZE;
use std::mem::MaybeUninit;

// No NCBI counterpart: reads the LOSAT_X_BPSEED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSEED): a strict-set switch, on unless it is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set.
pub(crate) fn x_bp_seed() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPSEED"))
}

// No NCBI counterpart: how many hits ahead the fused loop reads the diagonal entries
// (LOSAT_X_BPSEED_PREFETCH, default 16; 0 = no prefetch). Tuning only.
pub(super) fn x_prefetch_distance() -> usize {
    use std::sync::OnceLock;
    static D: OnceLock<usize> = OnceLock::new();
    *D.get_or_init(|| {
        std::env::var("LOSAT_X_BPSEED_PREFETCH")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(16)
    })
}

// No NCBI counterpart: a cache hint for one diagonal entry; it reads nothing and changes no value.
#[inline(always)]
pub(super) fn x_prefetch(diag_array: &[DiagStruct], index: usize) {
    #[cfg(target_arch = "x86_64")]
    if let Some(entry) = diag_array.get(index) {
        // SAFETY: a prefetch of a valid address has no effect on program state.
        unsafe {
            std::arch::x86_64::_mm_prefetch(
                (entry as *const DiagStruct).cast::<i8>(),
                std::arch::x86_64::_MM_HINT_T0,
            );
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (diag_array, index);
}

// No NCBI counterpart: reads the LOSAT_X_BPDIAG switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDIAG): a strict-set switch, on unless it is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set. It acts in the fused loop only.
pub(crate) fn x_bp_diag() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDIAG"))
}

// No NCBI counterpart: reads the LOSAT_X_BPDEFER switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDEFER): a strict-set switch, on unless it is `0`, or
/// `LOSAT_X_DEFAULT=0` and it is not set. It acts on the AVX2 batches only.
pub(crate) fn x_bp_defer() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDEFER"))
}

// No NCBI counterpart: reads the LOSAT_X_BPSCANBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSCANBF): blastp's fused loop finds the cells with [`x_scan_cells_bf`]
/// (module header). Off unless set, to any value (not in the strict set).
pub(crate) fn x_bp_scan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BPSCANBF").is_some())
}

/// How the fused loop tests the hits of one lookup cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XDiagMode {
    /// One hit after the other ([`x_test_hit`]): the fused loop without LOSAT_X_BPDIAG.
    PerHit,
    /// LOSAT_X_BPDIAG with LOSAT_X_BPDIAG_SIMD=0 (diagnostic): batches of [`X_DIAG_LANES`] hits in
    /// scalar code ([`x_test_batch_scalar`]).
    Scalar,
    /// LOSAT_X_BPDIAG with AVX2 and LOSAT_X_BPDEFER off: batches of eight hits in AVX2 registers
    /// ([`x_test_batch_avx2`]); a cell's last batch of fewer than eight hits is a scalar batch.
    #[cfg(target_arch = "x86_64")]
    Avx2,
    /// LOSAT_X_BPDIAG with AVX2 and LOSAT_X_BPDEFER: the AVX2 batches with the candidates of a
    /// cell's batches of eight called after its last one ([`x_cell_avx2_defer`]).
    #[cfg(target_arch = "x86_64")]
    Avx2Defer,
}

impl XDiagMode {
    /// The name in the `[TIMING] x_bpdiag` line (and blastx's `[TIMING] blastx_diag`).
    pub(crate) fn name(self) -> &'static str {
        match self {
            XDiagMode::PerHit => "per-hit",
            XDiagMode::Scalar => "scalar",
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2 => "avx2",
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2Defer => "avx2-defer",
        }
    }
}

// No NCBI counterpart: whether the CPU has AVX2; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
fn x_avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

// No NCBI counterpart: chooses LOSAT's code path once (LOSAT_X_BPDIAG; the AVX2 batch where the
// CPU has AVX2, with LOSAT_X_BPDEFER the deferred calls, otherwise the per-hit loop;
// LOSAT_X_BPDIAG_SIMD=0 forces the scalar batch); it does not change any value NCBI computes.
fn x_diag_mode() -> XDiagMode {
    use std::sync::OnceLock;
    static MODE: OnceLock<XDiagMode> = OnceLock::new();
    *MODE.get_or_init(|| {
        if !x_bp_diag() {
            return XDiagMode::PerHit;
        }
        if std::env::var_os("LOSAT_X_BPDIAG_SIMD").is_some_and(|v| v.as_os_str() == "0") {
            return XDiagMode::Scalar;
        }
        #[cfg(target_arch = "x86_64")]
        if x_avx2_available() {
            return if x_bp_defer() {
                XDiagMode::Avx2Defer
            } else {
                XDiagMode::Avx2
            };
        }
        // Without AVX2 the scalar batch would be slower than the per-hit loop (module header).
        XDiagMode::PerHit
    })
}

// No NCBI counterpart: chooses LOSAT's code path once per call site (LOSAT_X_BXDIAG); it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): the cell test of blastx's fused seed loop
/// (`blastx/x_seed_diag.rs`): `"per-hit"` or `"scalar"` when asked by name (diagnostics), otherwise
/// the choice LOSAT_X_BPDIAG makes for blastp: the AVX2 batch where the CPU has AVX2, else the
/// per-hit loop (the scalar batch was slower than the per-hit loop; module header).
pub(crate) fn x_cell_test(requested: &str) -> XDiagMode {
    match requested {
        "per-hit" => XDiagMode::PerHit,
        "scalar" => XDiagMode::Scalar,
        _ => {
            #[cfg(target_arch = "x86_64")]
            if x_avx2_available() {
                return XDiagMode::Avx2;
            }
            XDiagMode::PerHit
        }
    }
}

// No NCBI counterpart: counters of LOSAT_X_BPDIAG for its timing line.
static X_DIAG_CELLS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_DIAG_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_DIAG_CANDIDATES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BPDIAG): the counts of one subject's scan (per cell, not per hit), added to
/// the totals once.
#[derive(Default)]
struct XDiagCounts {
    cells: u64,
    hits: u64,
    candidates: u64,
}

impl XDiagCounts {
    #[inline(always)]
    fn add_cell(&mut self, hits: usize, candidates: u32) {
        self.cells += 1;
        self.hits += hits as u64;
        self.candidates += u64::from(candidates);
    }

    fn add_to_totals(&self) {
        if self.cells != 0 {
            use std::sync::atomic::Ordering::Relaxed;
            X_DIAG_CELLS.fetch_add(self.cells, Relaxed);
            X_DIAG_HITS.fetch_add(self.hits, Relaxed);
            X_DIAG_CANDIDATES.fetch_add(self.candidates, Relaxed);
        }
    }
}

// No NCBI counterpart: a LOSAT timing line.
/// EXPERIMENT (LOSAT_X_BPDIAG): with the switch on, the timing line of the cells tested in batches
/// since the last one.
pub(crate) fn x_bp_diag_print_timing() {
    if x_bp_diag() {
        use std::sync::atomic::Ordering::Relaxed;
        eprintln!(
            "[TIMING] x_bpdiag: mode={} cells={} batched_hits={} candidates={}",
            x_diag_mode().name(),
            X_DIAG_CELLS.swap(0, Relaxed),
            X_DIAG_HITS.swap(0, Relaxed),
            X_DIAG_CANDIDATES.swap(0, Relaxed)
        );
    }
}

// No NCBI counterpart: NCBI indexes `diag_array` directly; this keeps the table with its mask so
// that the loops below read entries without a bounds check on a range checked once.
/// The diagonal table of the fused loop and its mask; `new` checks that every
/// `(q_off - s_off) & mask` names an entry.
struct XDiagTable<'a> {
    entries: &'a mut [DiagStruct],
    mask: i32,
}

impl<'a> XDiagTable<'a> {
    fn new(entries: &'a mut [DiagStruct], mask: i32) -> Self {
        assert!(mask >= 0 && entries.len() > mask as usize);
        Self { entries, mask }
    }

    /// The entry index of a hit.
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:516
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    // ```
    #[inline(always)]
    fn coord(&self, query_offset: i32, subject_offset: i32) -> usize {
        (query_offset.wrapping_sub(subject_offset) & self.mask) as usize
    }

    #[inline(always)]
    fn entry(&mut self, query_offset: i32, subject_offset: i32) -> &mut DiagStruct {
        let coord = self.coord(query_offset, subject_offset);
        // SAFETY: coord <= mask < entries.len() (checked in `new`).
        unsafe { self.entries.get_unchecked_mut(coord) }
    }

    #[inline(always)]
    fn prefetch(&self, query_offset: i32, subject_offset: i32) {
        x_prefetch(self.entries, self.coord(query_offset, subject_offset));
    }
}

/// The values of the diagonal test that stay the same for a subject.
#[derive(Clone, Copy, Debug)]
struct XDiagTest {
    diag_offset: i32,
    window: i32,
    word_size: i32,
    /// LOSAT_X_BPSEED_PREFETCH (the per-hit loop only).
    prefetch: usize,
}

// No NCBI counterpart: reads the LOSAT_X_BPSEEDSHADOW switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSEEDSHADOW): run the fused loop and the reference for every subject from the
/// same diagonal table and stop (panic) if the initial HSPs or the table differ. Off unless set.
pub(crate) fn x_bp_seed_shadow() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BPSEEDSHADOW").is_some())
}

// No NCBI counterpart: counters of the shadow mode and of the fused loop, for the line at exit.
static X_SHADOW_SUBJECTS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_SHADOW_HSPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_SHADOW_SHORT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BPSEEDSHADOW): compares one subject's two runs; panics on a difference.
///
/// `compare_table` is false for a subject shorter than a word less one letter (one residue with
/// word size 3): the reference's only word there ends one byte past the subject's buffer, whose
/// unknown value can make the word index name a cell with hits (NCBI reads the same byte). Such
/// a hit is the subject's only scan position, so it is never extended, and the entry it leaves
/// is `window` plus the subject's length behind the next subject's first offset, which every
/// later hit on the diagonal reads as no previous hit, like any other entry of an earlier
/// subject. The initial HSPs are compared; the tables then may differ in that entry only.
pub(crate) fn x_shadow_compare(
    fused_hsps: &[String],
    reference_hsps: &[String],
    fused_diag: &[DiagStruct],
    reference_diag: &[DiagStruct],
    compare_table: bool,
) {
    if fused_hsps != reference_hsps {
        let at = fused_hsps
            .iter()
            .zip(reference_hsps)
            .position(|(a, b)| a != b)
            .unwrap_or(fused_hsps.len().min(reference_hsps.len()));
        panic!(
            "LOSAT_X_BPSEEDSHADOW: initial HSPs differ at {at} (fused {} / reference {}): fused {:?} reference {:?}",
            fused_hsps.len(),
            reference_hsps.len(),
            fused_hsps.get(at),
            reference_hsps.get(at)
        );
    }
    if !compare_table {
        X_SHADOW_SHORT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    } else if let Some(at) = fused_diag
        .iter()
        .zip(reference_diag)
        .position(|(a, b)| a.raw_bits() != b.raw_bits())
    {
        panic!(
            "LOSAT_X_BPSEEDSHADOW: diagonal entry {at} differs: fused {:#x} reference {:#x}",
            fused_diag[at].raw_bits(),
            reference_diag[at].raw_bits()
        );
    }
    assert_eq!(fused_diag.len(), reference_diag.len());
    X_SHADOW_SUBJECTS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    X_SHADOW_HSPS.fetch_add(
        fused_hsps.len() as u64,
        std::sync::atomic::Ordering::Relaxed,
    );
}

/// EXPERIMENT (LOSAT_X_BPSEEDSHADOW): one line at exit with what the shadow mode compared.
pub fn x_bp_seed_print_summary() {
    if x_bp_seed_shadow() {
        let subjects = X_SHADOW_SUBJECTS.load(std::sync::atomic::Ordering::Relaxed);
        let hsps = X_SHADOW_HSPS.load(std::sync::atomic::Ordering::Relaxed);
        let short = X_SHADOW_SHORT.load(std::sync::atomic::Ordering::Relaxed);
        if subjects == 0 {
            eprintln!("[X_SHADOW] BPSEED nothing compared (no two-hit scan of an amino-acid lookup table ran)");
        } else if short == 0 {
            eprintln!("[X_SHADOW] BPSEED subjects={subjects} initial_hsps={hsps} identical");
        } else {
            eprintln!(
                "[X_SHADOW] BPSEED subjects={subjects} initial_hsps={hsps} identical (tables not compared for {short} subjects shorter than a word)"
            );
        }
    }
}

/// The diagonal test of one word hit as NCBI applies it before an extension, on the hit's own
/// entry. Returns `Some((last_hit, diff))` when the hit needs the rest of the test.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
//             /* If the reset bit is set, an extension just happened. */
//             if (diag_array[diag_coord].flag) {
//                 /* If we've already extended past this hit, skip it. */
//                 if ((Int4) (subject_offset + diag_offset) <
//                     diag_array[diag_coord].last_hit) {
//                     continue;
//                 }
//                 /* Otherwise, start a new hit. */
//                 else {
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     diag_array[diag_coord].flag = 0;
//                 }
//             }
//             /* If the reset bit is cleared, try to start an extension. */
//             else {
//                 /* find the distance to the last hit on this diagonal */
//                 last_hit = diag_array[diag_coord].last_hit - diag_offset;
//                 diff = subject_offset - last_hit;
//
//                 if (diff >= window) {
//                     /* We are beyond the window for this diagonal; start a
//                        new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
//
//                 /* If the difference is less than the wordsize (i.e. last
//                    hit and this hit overlap), give up */
//
//                 if (diff < wordsize) {
//                     continue;
//                 }
// ```
//
// The sums and differences wrap (two's complement i32): that is what `+` and `-` do in release
// builds (no overflow checks), and what NCBI's arithmetic on its `Uint4 subject_offset` gives.
// LOSAT_X_BPDIAG's batches compute the same wrapping values per lane, and their tests compare
// with this function on offsets whose sum passes 2^31 - 1.
#[inline(always)]
pub(super) fn x_test_hit(
    entry: &mut DiagStruct,
    subject_offset: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
) -> Option<(i32, i32)> {
    if entry.flag() != 0 {
        let subject_plus_offset = subject_offset.wrapping_add(diag_offset);
        if subject_plus_offset >= entry.last_hit() {
            entry.set_last_hit(subject_plus_offset);
            entry.set_flag(0);
        }
        return None;
    }
    let last_hit = entry.last_hit().wrapping_sub(diag_offset);
    let diff = subject_offset.wrapping_sub(last_hit);
    if diff >= window {
        entry.set_last_hit(subject_offset.wrapping_add(diag_offset));
        return None;
    }
    if diff < word_size {
        return None;
    }
    Some((last_hit, diff))
}

/// Scans `subject` with `lookup` (two-hit mode) and applies the diagonal test to every word hit
/// in scan order; `candidate(query_offset, subject_offset, last_hit, diff, entry)` receives each
/// hit that needs the rest of the test, with its diagonal entry.
///
/// `subject` starts at the first residue and keeps the sentinel after the last one, as the
/// reference scan reads it; `subject_length` is the number of residues.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:495-512
// ```c
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
//
//         totalhits += hits;
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ```
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit<F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let test = XDiagTest {
        diag_offset,
        window,
        word_size,
        prefetch: x_prefetch_distance(),
    };
    let table = &mut XDiagTable::new(diag_array, diag_mask);
    // EXPERIMENT (LOSAT_X_BPSCANBF): the same cells found without a branch on the presence bit.
    let counts = if x_bp_scan_bf() {
        x_scan_two_hit_bf(
            x_diag_mode(),
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    } else {
        x_scan_two_hit_in(
            x_diag_mode(),
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    };
    counts.add_to_totals();
}

/// EXPERIMENT (LOSAT_X_BPSCANBF): [`x_scan_two_hit_in`] with the cells found by [`x_scan_cells_bf`].
/// Out of line, and with a closure of its own type, so that the default instance in
/// [`x_scan_two_hit`] and the cell tests it calls are compiled as without the switch: with the
/// same closure type the cell test had two callers, and the default instance ran 2-6 % more
/// instructions.
#[inline(never)]
fn x_scan_two_hit_bf<F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut own = |query_offset, subject_offset, last_hit, diff, entry: &mut DiagStruct| {
        candidate(query_offset, subject_offset, last_hit, diff, entry)
    };
    x_scan_two_hit_in_with::<true, _>(mode, lookup, subject, subject_length, table, test, &mut own)
}

/// EXPERIMENT (LOSAT_X_BXDIAG): [`x_scan_two_hit`] with the cell test `mode` (blastx's choice,
/// [`x_cell_test`], not LOSAT_X_BPDIAG's); returns the cells tested in batches, their hits and the
/// hits passed to `candidate` (zeros for [`XDiagMode::PerHit`]). Nothing is added to the totals of
/// the `[TIMING] x_bpdiag` line. EXPERIMENT (LOSAT_X_BXSCANBF): `branch_free` finds the cells with
/// [`x_scan_cells_bf`] (module header).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:495-512 (quoted at `x_scan_two_hit`)
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit_with<F>(
    mode: XDiagMode,
    branch_free: bool,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) -> [u64; 3]
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let test = XDiagTest {
        diag_offset,
        window,
        word_size,
        prefetch: x_prefetch_distance(),
    };
    let table = &mut XDiagTable::new(diag_array, diag_mask);
    // A constant per instantiation, so that each scan loop has one call of the cell test and the
    // compiler inlines it, as in blastp's loop.
    let counts = if branch_free {
        x_scan_two_hit_in_with::<true, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    } else {
        x_scan_two_hit_in_with::<false, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    };
    [counts.cells, counts.hits, counts.candidates]
}

/// [`x_scan_two_hit`] with the cell test of `mode`; returns the counts of LOSAT_X_BPDIAG (zero
/// for [`XDiagMode::PerHit`]).
fn x_scan_two_hit_in<F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_scan_two_hit_in_with::<false, F>(
        mode,
        lookup,
        subject,
        subject_length,
        table,
        test,
        candidate,
    )
}

/// [`x_scan_two_hit_in`], with the cells found by [`x_scan_cells_bf`] when `BRANCH_FREE`
/// (EXPERIMENT, LOSAT_X_BXSCANBF).
fn x_scan_two_hit_in_with<const BRANCH_FREE: bool, F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut counts = XDiagCounts::default();
    match mode {
        XDiagMode::PerHit => {
            x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
                x_cell_per_hit(hits, s, table, test, candidate)
            })
        }
        XDiagMode::Scalar => {
            x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
                let candidates = x_cell_scalar(hits, s, table, test, candidate);
                counts.add_cell(hits.len(), candidates);
            })
        }
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2 => {
            assert!(x_avx2_available(), "the AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                x_scan_cells_avx2::<BRANCH_FREE, F>(
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                    &mut counts,
                )
            }
        }
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2Defer => {
            assert!(x_avx2_available(), "the deferred AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                x_scan_cells_avx2_defer::<BRANCH_FREE, F>(
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                    &mut counts,
                )
            }
        }
    }
    counts
}

/// [`x_scan_cells`] with [`x_cell_avx2`], compiled for AVX2 as a whole so that the cell test is
/// inlined into the scan.
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_scan_cells_avx2<const BRANCH_FREE: bool, F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates = unsafe { x_cell_avx2(hits, s, table, test, candidate) };
        counts.add_cell(hits.len(), candidates);
    });
}

thread_local! {
    /// EXPERIMENT (LOSAT_X_BPDEFER): the candidate buffer of a thread's scans, kept between them.
    static X_DEFER_BUFFER: std::cell::Cell<Vec<i32>> = const { std::cell::Cell::new(Vec::new()) };
}

/// [`x_scan_cells_avx2`] with [`x_cell_avx2_defer`] (LOSAT_X_BPDEFER).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_scan_cells_avx2_defer<const BRANCH_FREE: bool, F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut buffer = X_DEFER_BUFFER.take();
    x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates = unsafe { x_cell_avx2_defer(hits, s, table, test, candidate, &mut buffer) };
        counts.add_cell(hits.len(), candidates);
    });
    X_DEFER_BUFFER.set(buffer);
}

/// EXPERIMENT (LOSAT_X_BXSCANBF): [`x_scan_cells`], or [`x_scan_cells_bf`] when `BRANCH_FREE` (the
/// same calls). One loop per instantiation, so `on_cell` has one call site and is inlined.
#[inline(always)]
fn x_scan_cells_with<const BRANCH_FREE: bool, C>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    on_cell: C,
) where
    C: FnMut(i32, &[i32]),
{
    if BRANCH_FREE {
        x_scan_cells_bf(lookup, subject, subject_length, on_cell);
    } else {
        x_scan_cells(lookup, subject, subject_length, on_cell);
    }
}

/// The subject scan of [`x_scan_two_hit`]: `on_cell(subject_offset, hits)` for every scan position
/// whose lookup cell has hits, in scan order, with the cell's list in its stored order.
// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99
// ```c
//    for (s = s_first; s <= s_last; s++) {
//       /* compute the index value */
//       index = ComputeTableIndexIncremental(word_length, lut->charsize,
//                                            lut->mask, s, index);
//       /* if there are hits... */
//       if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
//          Int4 numhits = bbcell[index].num_used;
//          ...
//          if (numhits <= AA_HITS_PER_CELL)
//             /* hits live in thick_backbone */
//             src = bbcell[index].payload.entries;
//          else
//             /* hits live in overflow array */
//             src = & (lut->overflow[bbcell[index].payload.overflow_cursor]);
//
//          /* copy the hits. */
//          for (i = 0; i < numhits; i++) {
//             offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//             offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
//          }
// ```
#[inline(always)]
fn x_scan_cells<C>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    mut on_cell: C,
) where
    C: FnMut(i32, &[i32]),
{
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let s_first = 0i32;
    let mut s_last =
        i32::try_from(subject_length).expect("subject length must fit Int4") - lookup.word_length;
    if s_last < s_first {
        s_last = s_first;
    }
    // A subject shorter than a word has one scan position whose word runs past the sentinel; the
    // reference reads that byte unchecked (one byte past the buffer for a one-residue subject).
    // It is read as a sentinel here. With one scan position no two hits share a diagonal, and the
    // next subject starts window + length further on the diagonals, so the byte reaches no output.
    let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0) as usize;
    let pv = &lookup.pv;
    let backbone = &lookup.backbone;
    let overflow = &lookup.overflow;

    let mut index = 0usize;
    for i in 0..word_length - 1 {
        index = (index << charsize) | byte_at(i);
    }
    for s in s_first..=s_last {
        let new_char = byte_at(s as usize + word_length - 1);
        index = ((index << charsize) | new_char) & mask;
        if pv[index >> PV_ARRAY_BTS] & (1u32 << (index & PV_ARRAY_MASK)) == 0 {
            continue;
        }
        let cell = &backbone[index];
        let numhits = cell.num_used as usize;
        let hits: &[i32] = if numhits <= AA_HITS_PER_CELL {
            &cell.entries[..numhits]
        } else {
            let cursor = cell.entries[0] as usize;
            &overflow[cursor..cursor + numhits]
        };
        on_cell(s, hits);
    }
}

/// EXPERIMENT (LOSAT_X_BXSCANBF): [`x_scan_cells`] with the cells that have hits found without a
/// branch on the presence bit: for each block of [`X_BLOCK`] positions [`x_block_candidates`] lists
/// `position << 32 | index`, the list's length advanced by the bit, and the listed cells then go to
/// `on_cell` in position order. The positions, indices and presence bits are those of
/// [`x_scan_cells`]'s loop (the same rolling index over the same bytes), so the calls are the same.
// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99 (quoted at `x_scan_cells`)
#[inline(always)]
fn x_scan_cells_bf<C>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    mut on_cell: C,
) where
    C: FnMut(i32, &[i32]),
{
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    i32::try_from(subject_length).expect("subject length must fit Int4");
    // The bytes x_scan_cells reads: positions 0..=s_last, s_last = max(subject_length -
    // word_length, 0), each with its word. A subject shorter than a word has one position whose
    // word runs past the subject, and x_scan_cells reads the missing bytes as 0: here they come
    // from a zero-padded copy.
    let mut padded = [0u8; 8];
    let (bytes, s_last) = if subject_length >= word_length {
        (&subject[..subject_length], subject_length - word_length)
    } else {
        assert!(word_length <= padded.len(), "word size {word_length}");
        let n = subject.len().min(word_length);
        padded[..n].copy_from_slice(&subject[..n]);
        (&padded[..word_length], 0)
    };
    // The backbone has NCBI's `backbone_size` cells (the index of a word of letters below
    // BLASTAA_SIZE, plus one) and the presence vector a bit for each, so the unchecked presence
    // reads of `x_block_candidates` stay inside it when every byte read is below BLASTAA_SIZE
    // (x_scan_cells would panic on an index past the vector instead).
    let largest = bytes.iter().copied().max().unwrap_or(0);
    assert!(
        usize::from(largest) < BLASTAA_SIZE
            && lookup.pv.len() << PV_ARRAY_BTS >= lookup.backbone.len(),
        "letter {largest} outside the lookup table's alphabet"
    );
    let backbone = &lookup.backbone;
    let overflow = &lookup.overflow;

    let mut index = 0usize;
    for &byte in &bytes[..word_length - 1] {
        index = (index << charsize) | byte as usize;
    }
    let mut cand = [MaybeUninit::<u64>::uninit(); X_BLOCK];
    let mut block_first = 0usize;
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        // Reads bytes[block_last + word_length - 1] at most, the last byte of `bytes`.
        let n = x_block_candidates(
            bytes,
            lookup.pv.as_ptr(),
            (word_length, charsize, mask),
            (block_first, block_last),
            &mut index,
            &mut cand,
        );
        for c in &cand[..n] {
            // SAFETY: the first `n` candidates were written by `x_block_candidates`.
            let c = unsafe { c.assume_init() };
            let cell = &backbone[(c & 0xffff_ffff) as usize];
            let numhits = cell.num_used as usize;
            let hits: &[i32] = if numhits <= AA_HITS_PER_CELL {
                &cell.entries[..numhits]
            } else {
                let cursor = cell.entries[0] as usize;
                &overflow[cursor..cursor + numhits]
            };
            on_cell((c >> 32) as i32, hits);
        }
        block_first = block_last + 1;
    }
}

/// The diagonal test of a cell's hits one after the other (the fused loop without
/// LOSAT_X_BPDIAG).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[inline(always)]
fn x_cell_per_hit<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    for (i, &query_offset) in hits.iter().enumerate() {
        // Read ahead the entry of a later hit of this cell: the hits of one subject offset
        // fall on distinct diagonals spread over the table, so each is a cache miss that the
        // hardware cannot predict. A prefetch reads nothing the test depends on.
        if let Some(&ahead) = hits.get(i + test.prefetch) {
            table.prefetch(ahead, s);
        }
        let entry = table.entry(query_offset, s);
        if let Some((last_hit, diff)) =
            x_test_hit(entry, s, test.diag_offset, test.window, test.word_size)
        {
            candidate(query_offset, s, last_hit, diff, entry);
        }
    }
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the hits of one batch, the lanes of an AVX2 register.
const X_DIAG_LANES: usize = 8;

/// EXPERIMENT (LOSAT_X_BPDIAG): the 31-bit `last_hit` field of a packed entry.
const X_LAST_HIT_MASK: u32 = 0x7fff_ffff;

// No NCBI counterpart: a check of the property the batches rest on (module header).
/// EXPERIMENT (LOSAT_X_BPDIAG): whether the entries of a batch's hits are distinct.
fn x_batch_coords_distinct(coords: &[usize]) -> bool {
    coords
        .iter()
        .enumerate()
        .all(|(i, coord)| !coords[i + 1..].contains(coord))
}

/// EXPERIMENT (LOSAT_X_BPDIAG): calls `candidate` for the lanes set in `lanes` (bit `i` = hit `i`
/// of `hits`) in lane order, which is the cell's order, with the values [`x_test_hit`] returns.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[inline(always)]
fn x_call_candidates<F>(
    mut lanes: u32,
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    while lanes != 0 {
        let lane = lanes.trailing_zeros() as usize;
        lanes &= lanes - 1;
        let query_offset = hits[lane];
        let entry = table.entry(query_offset, s);
        // The batch stored a candidate's old value back, so this is the entry the test read.
        let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
        let diff = s.wrapping_sub(last_hit);
        candidate(query_offset, s, last_hit, diff, entry);
    }
}

/// EXPERIMENT (LOSAT_X_BPDIAG): [`x_test_hit`] for up to [`X_DIAG_LANES`] hits of one cell: for
/// each hit the new entry and whether it is a candidate, without branches on the entry; every
/// entry is stored, then the candidates are called in order. Returns the number of candidates.
///
/// With `flag` the entry's bit 31, `last` its sign-extended 31-bit `last_hit`,
/// `sum = s + diag_offset` and `diff = s - (last - diag_offset)` (wrapping, as in `x_test_hit`):
/// - write = flag ? sum >= last : diff >= window; the new entry is then `sum` truncated to 31
///   bits with the flag cleared (`set_last_hit` + `set_flag(0)`, or `set_last_hit` on an entry
///   without the flag), otherwise the old entry;
/// - candidate = !flag && word_size <= diff < window (its entry is not written).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[inline(always)]
fn x_test_batch_scalar<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    debug_assert!(hits.len() <= X_DIAG_LANES);
    let subject_plus_offset = s.wrapping_add(test.diag_offset);
    let written = (subject_plus_offset as u32) & X_LAST_HIT_MASK;
    let mut lanes = 0u32;
    for (lane, &query_offset) in hits.iter().enumerate() {
        let entry = table.entry(query_offset, s);
        let old = entry.raw_bits();
        let flag = old >> 31 != 0;
        let last = entry.last_hit();
        let diff = s.wrapping_sub(last.wrapping_sub(test.diag_offset));
        let write = if flag {
            subject_plus_offset >= last
        } else {
            diff >= test.window
        };
        entry.set_raw_bits(if write { written } else { old });
        let is_candidate = !flag & (diff >= test.word_size) & (diff < test.window);
        lanes |= u32::from(is_candidate) << lane;
    }
    debug_assert!(x_batch_coords_distinct(
        &hits
            .iter()
            .map(|&query_offset| table.coord(query_offset, s))
            .collect::<Vec<_>>()
    ));
    x_call_candidates(lanes, hits, s, table, test, candidate);
    lanes.count_ones()
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the cell's hits in batches of [`X_DIAG_LANES`], each by
/// [`x_test_batch_scalar`]; returns the number of candidates. The batches read no entry ahead
/// (LOSAT_X_BPSEED_PREFETCH is the per-hit loop's): prefetching the batch 8 or 16 hits ahead was
/// not faster than none on 2026-10-14 (local, human block 1 x Species4 and the OrthoFinder
/// sample job, one thread).
#[inline(always)]
fn x_cell_scalar<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut candidates = 0;
    for batch in hits.chunks(X_DIAG_LANES) {
        candidates += x_test_batch_scalar(batch, s, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPDIAG): [`x_test_batch_scalar`] for exactly eight hits in AVX2 registers.
/// The lane arithmetic is the scalar code's: `_mm256_sub_epi32` wraps like `wrapping_sub`, and
/// `_mm256_cmpgt_epi32` is the signed i32 comparison. All eight entries are stored (the old value
/// where the test writes nothing): storing only the written lanes, from the movemask, was slower
/// on 2026-10-14 (local, human block 1 x Species4, one thread).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_test_batch_avx2<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let (_, lanes) = x_batch_avx2(hits, s, table, test);
    x_call_candidates(lanes, hits, s, table, test, candidate);
    lanes.count_ones()
}

/// The test of [`x_test_batch_avx2`] without the calls: stores the batch's entries and returns the
/// hits' register and the candidate lanes (bit `i` = hit `i`).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_batch_avx2(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> (std::arch::x86_64::__m256i, u32) {
    use std::arch::x86_64::*;
    assert!(hits.len() == X_DIAG_LANES);
    let subject_plus_offset = s.wrapping_add(test.diag_offset);
    let written = (subject_plus_offset as u32) & X_LAST_HIT_MASK;
    let query = _mm256_loadu_si256(hits.as_ptr().cast());
    // (q - s) & diag_mask per lane
    let coords = _mm256_and_si256(
        _mm256_sub_epi32(query, _mm256_set1_epi32(s)),
        _mm256_set1_epi32(table.mask),
    );
    let mut lane_coords = [0u32; X_DIAG_LANES];
    _mm256_storeu_si256(lane_coords.as_mut_ptr().cast(), coords);
    // The eight entries with scalar loads, not `_mm256_i32gather_epi32`: on the workstation (Xeon
    // 6338, Gather Data Sampling microcode mitigation on) the gather was slower, 2026-10-14, one
    // pinned core, scan_ungapped: human block 1 per-hit 11.59 s, gather 10.17 s, scalar loads
    // 9.13 s; C. elegans block 5 32.68, 28.55, 26.45 s. Locally (no mitigation) both were equal.
    let mut lane_old = [0u32; X_DIAG_LANES];
    for (old, &coord) in lane_old.iter_mut().zip(&lane_coords) {
        // coord <= mask < entries.len() (XDiagTable::new)
        *old = table.entries.get_unchecked(coord as usize).raw_bits();
    }
    let old = _mm256_loadu_si256(lane_old.as_ptr().cast());
    // flag: all ones where bit 31 is set; last: the 31-bit field sign-extended (last_hit()).
    let flag = _mm256_srai_epi32::<31>(old);
    let last = _mm256_srai_epi32::<1>(_mm256_slli_epi32::<1>(old));
    let diff = _mm256_sub_epi32(
        _mm256_set1_epi32(s),
        _mm256_sub_epi32(last, _mm256_set1_epi32(test.diag_offset)),
    );
    // keep (no write) = flag ? sum < last : diff < window
    let sum_below_last = _mm256_cmpgt_epi32(last, _mm256_set1_epi32(subject_plus_offset));
    let in_window = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.window), diff);
    let keep = _mm256_blendv_epi8(in_window, sum_below_last, flag);
    // candidate = !flag && diff >= word_size && diff < window
    let below_word = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.word_size), diff);
    let candidates = _mm256_andnot_si256(_mm256_or_si256(flag, below_word), in_window);
    let new = _mm256_blendv_epi8(_mm256_set1_epi32(written as i32), old, keep);
    let mut lane_new = [0u32; X_DIAG_LANES];
    _mm256_storeu_si256(lane_new.as_mut_ptr().cast(), new);
    for (&coord, &raw) in lane_coords.iter().zip(&lane_new) {
        // coord <= mask < entries.len()
        table
            .entries
            .get_unchecked_mut(coord as usize)
            .set_raw_bits(raw);
    }
    debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
    let lanes = _mm256_movemask_ps(_mm256_castsi256_ps(candidates)) as u32;
    (query, lanes)
}

/// EXPERIMENT (LOSAT_X_BPDEFER): for each mask of eight lanes, the set lanes in order, then zeros:
/// the permutation that moves a batch's candidates to the front of its register.
const X_PACK_LANES: [[u8; X_DIAG_LANES]; 256] = {
    let mut table = [[0u8; X_DIAG_LANES]; 256];
    let mut mask = 0;
    while mask < 256 {
        let mut n = 0;
        let mut lane = 0;
        while lane < X_DIAG_LANES {
            if mask >> lane & 1 != 0 {
                table[mask][n] = lane as u8;
                n += 1;
            }
            lane += 1;
        }
        mask += 1;
    }
    table
};

/// EXPERIMENT (LOSAT_X_BPDEFER): [`x_cell_avx2`] with the calls of the batches of eight after the
/// last of them (module header): each batch packs its candidates' query offsets into `buffer`, then
/// the closure is called for the buffer in order, as [`x_call_candidates`] calls it, and the last
/// batch of fewer than eight hits is a scalar batch. Returns the number of candidates.
///
/// # Safety
/// The CPU must have AVX2.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2_defer<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    buffer: &mut Vec<i32>,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let mut batches = hits.chunks_exact(X_DIAG_LANES);
    // Room for a full register past the last candidate (each batch stores eight lanes).
    buffer.clear();
    buffer.reserve(hits.len() + X_DIAG_LANES);
    let out = buffer.as_mut_ptr();
    let mut packed = 0usize;
    for batch in &mut batches {
        let (query, lanes) = x_batch_avx2(batch, s, table, test);
        let order = _mm256_cvtepu8_epi32(_mm_loadl_epi64(
            X_PACK_LANES[lanes as usize].as_ptr().cast(),
        ));
        // packed <= the hits of the batches before this one, so the eight stores stay inside
        // the reserved capacity.
        _mm256_storeu_si256(
            out.add(packed).cast(),
            _mm256_permutevar8x32_epi32(query, order),
        );
        packed += lanes.count_ones() as usize;
    }
    // The first `packed` values were stored above.
    buffer.set_len(packed);
    for &query_offset in buffer.iter() {
        let entry = table.entry(query_offset, s);
        // The batch stored a candidate's old value back, and no later batch of the cell touches
        // its entry, so this is the entry the test read.
        let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
        let diff = s.wrapping_sub(last_hit);
        candidate(query_offset, s, last_hit, diff, entry);
    }
    let mut candidates = packed as u32;
    let rest = batches.remainder();
    if !rest.is_empty() {
        candidates += x_test_batch_scalar(rest, s, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the cell's hits in batches of eight by [`x_test_batch_avx2`], the
/// last batch of fewer than eight by [`x_test_batch_scalar`]; returns the number of candidates.
/// No entry is read ahead (see [`x_cell_scalar`]).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut candidates = 0;
    let mut batches = hits.chunks_exact(X_DIAG_LANES);
    for batch in &mut batches {
        candidates += x_test_batch_avx2(batch, s, table, test, candidate);
    }
    let rest = batches.remainder();
    if !rest.is_empty() {
        candidates += x_test_batch_scalar(rest, s, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPSEED): the context of an offset of the concatenated query without the
/// binary search over the contexts: a coarse table gives the context at the start of the offset's
/// block, and a walk over the (small, contiguous) frame bases moves on to the last context starting
/// at or before the offset. The result is the largest context index whose `frame_base` is at most
/// the offset (0 below the first), which is what the binary search of
/// `blastp_context_idx_for_query_offset_with_bounds` returns. The coarse table has at most
/// `COARSE_ENTRIES` entries so that both arrays stay in the first cache levels: a table with one
/// entry per offset (as large as the diagonal table) cost a cache miss for every candidate hit.
// NCBI reference: c++/src/algo/blast/core/blast_query_info.c:211-232
// ```c
// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
// {
//     Int4 m=0, b=0, e=0, size=0;
//
//     size = A->last_context+1;
//     ...
//     while (b < e - 1) {
//         m = (b + e) / 2;
//         if (A->contexts[m].query_offset > n)
//             e = m;
//         else
//             b = m;
//     }
//     return b;
// }
// ```
pub(crate) struct XContextTable {
    /// Frame bases in context order, then `i32::MAX`.
    starts: Vec<i32>,
    /// Context index at offset `b << shift`.
    coarse: Vec<u32>,
    shift: u32,
    span: usize,
}

impl XContextTable {
    const COARSE_ENTRIES: usize = 1 << 14;

    /// `frame_bases` in context order (non-decreasing); `span` = the offsets to cover.
    pub(crate) fn new(frame_bases: &[i32], span: usize) -> Self {
        let mut shift = 0u32;
        while (span >> shift) >= Self::COARSE_ENTRIES {
            shift += 1;
        }
        Self::with_shift(frame_bases, span, shift)
    }

    pub(crate) fn with_shift(frame_bases: &[i32], span: usize, shift: u32) -> Self {
        let mut starts = Vec::with_capacity(frame_bases.len() + 2);
        starts.extend_from_slice(frame_bases);
        starts.push(i32::MAX);
        if frame_bases.is_empty() {
            // no context: every offset gives index 0, as the binary search does
            starts.push(i32::MAX);
        }
        let blocks = (span >> shift) + 1;
        let mut coarse = Vec::with_capacity(blocks);
        let mut idx = 0usize;
        for b in 0..blocks {
            let offset = i64::try_from(b << shift).expect("query offset must fit i64");
            while idx + 1 < frame_bases.len() && i64::from(starts[idx + 1]) <= offset {
                idx += 1;
            }
            coarse.push(u32::try_from(idx).expect("context index must fit u32"));
        }
        Self {
            starts,
            coarse,
            shift,
            span,
        }
    }

    /// The context of `query_offset`, or `None` outside the table.
    #[inline(always)]
    pub(crate) fn get(&self, query_offset: i32) -> Option<usize> {
        let q = usize::try_from(query_offset)
            .ok()
            .filter(|&q| q < self.span)?;
        let mut idx = self.coarse[q >> self.shift] as usize;
        // starts[last] = i32::MAX ends the walk; the coarse entry is a context index, so idx + 1
        // is always inside `starts`.
        while self.starts[idx + 1] <= query_offset {
            idx += 1;
        }
        Some(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(last_hits: &[i32]) -> Vec<DiagStruct> {
        last_hits
            .iter()
            .map(|&last_hit| {
                let mut entry = DiagStruct::default();
                entry.set_last_hit(last_hit);
                entry
            })
            .collect()
    }

    #[test]
    fn the_shadow_compares_the_table_unless_the_subject_is_shorter_than_a_word() {
        let hsps = vec!["hsp".to_string()];
        // A subject shorter than a word less one letter: only the initial HSPs are compared.
        x_shadow_compare(&hsps, &hsps, &table(&[0, 40]), &table(&[0, 0]), false);
        x_shadow_compare(&hsps, &hsps, &table(&[0, 40]), &table(&[0, 40]), true);
        let differs = std::panic::catch_unwind(|| {
            x_shadow_compare(&hsps, &hsps, &table(&[0, 40]), &table(&[0, 0]), true)
        });
        assert!(differs.is_err(), "a table difference stops the run");
        let hsps_differ = std::panic::catch_unwind(|| {
            x_shadow_compare(&hsps, &[], &table(&[0]), &table(&[0]), false)
        });
        assert!(
            hsps_differ.is_err(),
            "an initial HSP difference stops the run"
        );
    }

    // ---- EXPERIMENT (LOSAT_X_BPDIAG): the batched cell test against the per-hit test ----

    /// One call of the candidate closure: query offset, subject offset, last_hit, diff, and the
    /// entry before the call.
    type Call = (i32, i32, i32, i32, u32);

    fn next(seed: &mut u64) -> u32 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*seed >> 32) as u32
    }

    fn entry_of(raw: u32) -> DiagStruct {
        let mut entry = DiagStruct::default();
        entry.set_raw_bits(raw);
        entry
    }

    /// The modes this machine runs: the per-hit test (the reference), the scalar batch and, with
    /// AVX2, the AVX2 batch.
    fn cell_modes() -> Vec<XDiagMode> {
        #[allow(unused_mut)]
        let mut modes = vec![XDiagMode::PerHit, XDiagMode::Scalar];
        #[cfg(target_arch = "x86_64")]
        if x_avx2_available() {
            modes.push(XDiagMode::Avx2);
            modes.push(XDiagMode::Avx2Defer);
        } else {
            eprintln!("no AVX2 on this CPU: the AVX2 batch is not compared");
        }
        modes
    }

    /// The candidate closure of the tests: writes the candidate's entry as the extension does
    /// (with the flag for some), so a batch that stored its entries after the call would leave
    /// another table.
    fn extend(query_offset: i32, subject_offset: i32, entry: &mut DiagStruct, diag_offset: i32) {
        if (query_offset ^ subject_offset) & 1 == 0 {
            entry.set_flag(1);
            entry.set_last_hit(subject_offset.wrapping_add(diag_offset).wrapping_add(5));
        } else {
            entry.set_last_hit(subject_offset.wrapping_add(diag_offset));
        }
    }

    /// Tests the cells `(subject_offset, hits)` in order with `mode`; returns the closure's calls.
    /// A batch mode's candidate count must be its number of calls.
    fn run_cells(
        mode: XDiagMode,
        entries: &mut [DiagStruct],
        mask: i32,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
    ) -> Vec<Call> {
        let mut calls = Vec::new();
        let mut counted = 0u32;
        #[cfg(target_arch = "x86_64")]
        let mut buffer = Vec::new();
        let mut table = XDiagTable::new(entries, mask);
        let mut candidate = |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
            calls.push((q, s, last_hit, diff, entry.raw_bits()));
            extend(q, s, entry, test.diag_offset);
        };
        for (s, hits) in cells {
            match mode {
                XDiagMode::PerHit => x_cell_per_hit(hits, *s, &mut table, test, &mut candidate),
                XDiagMode::Scalar => {
                    counted += x_cell_scalar(hits, *s, &mut table, test, &mut candidate);
                }
                #[cfg(target_arch = "x86_64")]
                XDiagMode::Avx2 => {
                    assert!(x_avx2_available());
                    // SAFETY: the CPU has AVX2 (asserted above).
                    counted += unsafe { x_cell_avx2(hits, *s, &mut table, test, &mut candidate) };
                }
                #[cfg(target_arch = "x86_64")]
                XDiagMode::Avx2Defer => {
                    assert!(x_avx2_available());
                    // SAFETY: the CPU has AVX2 (asserted above).
                    counted += unsafe {
                        x_cell_avx2_defer(hits, *s, &mut table, test, &mut candidate, &mut buffer)
                    };
                }
            }
        }
        if mode != XDiagMode::PerHit {
            assert_eq!(counted as usize, calls.len(), "{mode:?}: candidate count");
        }
        calls
    }

    /// Every mode against the per-hit test from the same table: the same calls in the same order
    /// with the same arguments, and the same table bits. Returns the number of calls.
    fn same_in_every_mode(
        entries: &[DiagStruct],
        mask: i32,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
    ) -> usize {
        let mut reference = entries.to_vec();
        let reference_calls = run_cells(XDiagMode::PerHit, &mut reference, mask, test, cells);
        for mode in cell_modes() {
            let mut got = entries.to_vec();
            let calls = run_cells(mode, &mut got, mask, test, cells);
            assert_eq!(calls, reference_calls, "{mode:?} {test:?}: candidate calls");
            let at = got
                .iter()
                .zip(&reference)
                .position(|(a, b)| a.raw_bits() != b.raw_bits());
            assert!(
                at.is_none(),
                "{mode:?} {test:?}: entry {at:?}: {:#x} reference {:#x}",
                got[at.unwrap_or(0)].raw_bits(),
                reference[at.unwrap_or(0)].raw_bits()
            );
        }
        reference_calls.len()
    }

    /// `n` distinct query offsets from `base .. base + size` (so on distinct entries of a table of
    /// `size` entries), in random order.
    fn distinct_offsets(seed: &mut u64, size: usize, n: usize, base: i32) -> Vec<i32> {
        let mut all: Vec<i32> = (0..size as i32).collect();
        for i in 0..n {
            let j = i + next(seed) as usize % (size - i);
            all.swap(i, j);
        }
        all.truncate(n);
        all.into_iter().map(|q| q + base).collect()
    }

    /// Diagonal offsets with sums past 2^30 (the 31-bit field truncates) and past 2^31 - 1 (the
    /// i32 sum wraps) for subject offsets up to a few thousand.
    const DIAG_OFFSETS: [i32; 6] = [
        40,
        1_000,
        (1 << 30) - 600,
        (1 << 30) + 77,
        i32::MAX / 4,
        i32::MAX - 700,
    ];

    const HIT_COUNTS: [usize; 13] = [0, 1, 2, 7, 8, 9, 15, 16, 17, 24, 64, 333, 1024];

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[test]
    fn batched_cells_equal_the_per_hit_test_on_random_tables() {
        let mut seed = 0x5eed_u64;
        let size = 1024usize;
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for (window, word_size) in [(40, 3), (16, 2), (3, 5), (1, 0)] {
            for diag_offset in DIAG_OFFSETS {
                for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: [0, 3, 16][k % 3],
                    };
                    let s0 = (next(&mut seed) % 4000) as i32;
                    // Half the entries any bit pattern, half near the cells' subject offsets
                    // (flag set or not), so every outcome of the test occurs.
                    let entries: Vec<DiagStruct> = (0..size)
                        .map(|_| {
                            let r = next(&mut seed);
                            if r & 1 == 0 {
                                entry_of(next(&mut seed))
                            } else {
                                let back = (next(&mut seed) % (window as u32 + 20)) as i32 - 8;
                                let last = s0.wrapping_add(diag_offset).wrapping_sub(back);
                                entry_of((last as u32 & X_LAST_HIT_MASK) | (r & 2) << 30)
                            }
                        })
                        .collect();
                    // Later cells read what earlier cells (and their candidates) wrote.
                    let cells: Vec<(i32, Vec<i32>)> = [0, 1, 1, 3, 10, 41]
                        .iter()
                        .map(|&step| {
                            let base = (next(&mut seed) % 3) as i32 * 7;
                            (s0 + step, distinct_offsets(&mut seed, size, n, base))
                        })
                        .collect();
                    calls += same_in_every_mode(&entries, mask, test, &cells);
                }
            }
        }
        assert!(calls > 1000, "the tables must give candidates: {calls}");
    }

    /// The entry that gives the hit `(query_offset, s)` the outcome `case` of the test.
    fn adversarial_entry(case: usize, s: i32, test: XDiagTest) -> u32 {
        let sum = s.wrapping_add(test.diag_offset);
        // last_hit for a hit without the flag at distance `diff`: s - diff + diag_offset
        let at_distance = |diff: i32| s.wrapping_sub(diff).wrapping_add(test.diag_offset);
        let flagged = |last: i32| 0x8000_0000 | (last as u32 & X_LAST_HIT_MASK);
        let clear = |last: i32| last as u32 & X_LAST_HIT_MASK;
        match case % 11 {
            // flag set: s + diag_offset above, equal to, below last_hit
            0 => flagged(sum.wrapping_add(1)),
            1 => flagged(sum),
            2 => flagged(sum.wrapping_sub(1)),
            // flag clear: diff at the window and word size edges, 0, negative, far
            3 => clear(at_distance(test.window)),
            4 => clear(at_distance(test.window - 1)),
            5 => clear(at_distance(test.word_size)),
            6 => clear(at_distance(test.word_size - 1)),
            7 => clear(at_distance(0)),
            8 => clear(at_distance(-3)),
            9 => clear(at_distance(1 << 29)),
            _ => flagged(at_distance(test.window)),
        }
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[test]
    fn batched_cells_equal_the_per_hit_test_at_every_edge_of_the_test() {
        let mut seed = 0xed9e_u64;
        let size = 2048usize;
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for (window, word_size) in [(40, 3), (4, 3), (2, 2)] {
            for diag_offset in DIAG_OFFSETS {
                for s in [0, 5, 599, 600, 601, 699, 700, 701, 3000] {
                    for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                        let test = XDiagTest {
                            diag_offset,
                            window,
                            word_size,
                            prefetch: [16, 0, 1][k % 3],
                        };
                        let hits = distinct_offsets(&mut seed, size, n, 0);
                        let mut entries: Vec<DiagStruct> =
                            (0..size).map(|_| entry_of(next(&mut seed))).collect();
                        let rotate = next(&mut seed) as usize;
                        for (j, &q) in hits.iter().enumerate() {
                            let coord = ((q - s) & mask) as usize;
                            entries[coord] = entry_of(adversarial_entry(j + rotate, s, test));
                        }
                        // The same hits again one and window offsets later: the entries the first
                        // cell and its candidates wrote are tested again.
                        let cells =
                            vec![(s, hits.clone()), (s + 1, hits.clone()), (s + window, hits)];
                        calls += same_in_every_mode(&entries, mask, test, &cells);
                    }
                }
            }
        }
        assert!(calls > 1000, "the edges must give candidates: {calls}");
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:495-512 and blast_aascan.c:60-99
    // (quoted at `x_scan_two_hit` and `x_scan_cells`)
    #[test]
    fn batched_scans_of_real_lookup_tables_equal_the_per_hit_scan() {
        use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
        use crate::algorithm::tblastx::translation::QueryFrame;
        use crate::config::ScoringMatrix;
        use crate::stats::lookup_protein_params_ungapped;

        // NCBISTDAA codes of the 20 standard amino acids
        const STANDARD: [u8; 20] = [
            1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
        ];
        let mut seed = 0x1007_u64;
        let mut protein = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| STANDARD[next(&mut seed) as usize % STANDARD.len()])
                .collect()
        };
        let mut queries: Vec<Vec<u8>> = (0..6).map(|i| protein(80 + 90 * i)).collect();
        // Low complexity: cells far above AA_HITS_PER_CELL (the overflow array).
        let mut low = vec![11u8; 150];
        low.extend(protein(30));
        low.extend(std::iter::repeat_n([17u8, 18], 60).flatten());
        queries.push(low);
        let frames: Vec<Vec<QueryFrame>> = queries
            .iter()
            .map(|seq| {
                let mut aa_seq = vec![0u8];
                aa_seq.extend_from_slice(seq);
                aa_seq.push(0);
                vec![QueryFrame {
                    frame: 1,
                    aa_seq,
                    aa_seq_nomask: None,
                    aa_len: seq.len(),
                    orig_len: seq.len(),
                    seg_masks: Vec::new(),
                }]
            })
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, _) = build_ncbi_lookup(&frames, 11, &params, false);
        assert!(
            lookup
                .backbone
                .iter()
                .any(|cell| cell.num_used as usize > 3 * X_DIAG_LANES),
            "the query set must give long cells"
        );
        // Subjects: unrelated proteins, copies of the queries with every 5th letter changed (two
        // hits on a diagonal, extensions), the low-complexity query, and short ones.
        let mut subjects: Vec<Vec<u8>> = queries
            .iter()
            .map(|q| {
                q.iter()
                    .enumerate()
                    .map(|(i, &a)| if i % 5 == 4 { STANDARD[i % 20] } else { a })
                    .collect()
            })
            .collect();
        subjects.extend((0..4).map(|i| protein(50 + 200 * i)));
        subjects.push(queries[6].clone());
        subjects.push(protein(1));
        subjects.push(protein(3));
        let window = 40;
        let mut size = 1usize;
        while size < lookup.query_length as usize + window as usize {
            size <<= 1;
        }
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for start in [window, (1 << 30) - 2_000, i32::MAX - 3_000] {
            for prefetch in [16, 0, 5] {
                let mut results = Vec::new();
                // EXPERIMENT (LOSAT_X_BXSCANBF): every cell test also with the branch-free scan.
                let runs = cell_modes()
                    .into_iter()
                    .flat_map(|mode| [(mode, false), (mode, true)]);
                for (mode, branch_free) in runs {
                    let mut entries = vec![DiagStruct::clear(window); size];
                    let mut log: Vec<(usize, Call)> = Vec::new();
                    let mut diag_offset: i32 = start;
                    let mut counts = (0u64, 0u64, 0u64);
                    for (k, subject) in subjects.iter().enumerate() {
                        let mut buffer = subject.clone();
                        buffer.push(0);
                        let test = XDiagTest {
                            diag_offset,
                            window,
                            word_size: 3,
                            prefetch,
                        };
                        let mut candidate =
                            |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
                                log.push((k, (q, s, last_hit, diff, entry.raw_bits())));
                                extend(q, s, entry, diag_offset);
                            };
                        let table = &mut XDiagTable::new(&mut entries, mask);
                        let c = if branch_free {
                            x_scan_two_hit_in_with::<true, _>(
                                mode,
                                &lookup,
                                &buffer,
                                subject.len(),
                                table,
                                test,
                                &mut candidate,
                            )
                        } else {
                            x_scan_two_hit_in(
                                mode,
                                &lookup,
                                &buffer,
                                subject.len(),
                                table,
                                test,
                                &mut candidate,
                            )
                        };
                        counts = (
                            counts.0 + c.cells,
                            counts.1 + c.hits,
                            counts.2 + c.candidates,
                        );
                        diag_offset = diag_offset.wrapping_add(subject.len() as i32 + window);
                    }
                    if mode == XDiagMode::PerHit {
                        assert_eq!(counts, (0, 0, 0), "the per-hit loop counts nothing");
                    } else {
                        assert!(counts.0 > 0 && counts.1 >= counts.0, "{mode:?}: {counts:?}");
                        assert_eq!(counts.2 as usize, log.len(), "{mode:?}: candidates");
                    }
                    let bits: Vec<u32> = entries.iter().map(|e| e.raw_bits()).collect();
                    results.push(((mode, branch_free), log, bits));
                }
                let (_, reference_log, reference_bits) = &results[0];
                for (mode, log, bits) in &results[1..] {
                    assert!(
                        log == reference_log,
                        "{mode:?} start {start} prefetch {prefetch}: candidate calls differ"
                    );
                    assert!(
                        bits == reference_bits,
                        "{mode:?} start {start} prefetch {prefetch}: tables differ"
                    );
                }
                calls += reference_log.len();
            }
        }
        assert!(calls > 100, "the subjects must give candidates: {calls}");
    }

    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99 (quoted at `x_scan_cells`)
    // EXPERIMENT (LOSAT_X_BXSCANBF): the branch-free scan calls `on_cell` with the positions and
    // hit lists of the presence-test loop, on subjects across the block boundaries, with the
    // letters of translated frames (X, stop; all of 1-27) and with the word sizes 2-4 (the general
    // loop).
    #[test]
    fn the_branch_free_scan_finds_the_cells_of_the_presence_test() {
        use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
        use crate::algorithm::tblastx::translation::QueryFrame;
        use crate::config::ScoringMatrix;
        use crate::stats::lookup_protein_params_ungapped;

        let mut seed = 0x5ca7_u64;
        let mut letters = |n: usize, top: u64| -> Vec<u8> {
            (0..n)
                .map(|_| 1 + (u64::from(next(&mut seed)) % top) as u8)
                .collect()
        };
        let mut compared = 0usize;
        for (case, word_length) in [3, 3, 2, 4, 3].into_iter().enumerate() {
            let frames: Vec<Vec<QueryFrame>> = (0..1 + 3 * case)
                .map(|i| {
                    let seq = letters(60 + 70 * i, 22);
                    let mut aa_seq = vec![0u8];
                    aa_seq.extend_from_slice(&seq);
                    aa_seq.push(0);
                    vec![QueryFrame {
                        frame: 1,
                        aa_seq,
                        aa_seq_nomask: None,
                        aa_len: seq.len(),
                        orig_len: seq.len(),
                        seg_masks: Vec::new(),
                    }]
                })
                .collect();
            let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
            let (lookup, _) = build_ncbi_lookup_for_profile(
                &frames,
                [11, 13][case % 2],
                &params,
                false,
                ScoringMatrix::Blosum62,
                word_length,
            );
            for length in [0, 1, 2, 3, 4, 255, 256, 257, 258, 511, 512, 513, 2_000] {
                let mut subject = letters(length, 27);
                subject.push(0);
                let mut want: Vec<(i32, Vec<i32>)> = Vec::new();
                let mut got: Vec<(i32, Vec<i32>)> = Vec::new();
                x_scan_cells(&lookup, &subject, length, |s, hits| {
                    want.push((s, hits.to_vec()))
                });
                x_scan_cells_with::<true, _>(&lookup, &subject, length, |s, hits| {
                    got.push((s, hits.to_vec()))
                });
                assert_eq!(got, want, "case {case} length {length}");
                compared += want.len();
            }
        }
        assert!(compared > 1_000, "{compared} cells compared");
    }

    #[test]
    fn the_pack_table_lists_the_set_lanes_in_order() {
        for (mask, order) in X_PACK_LANES.iter().enumerate() {
            let lanes: Vec<u8> = (0..X_DIAG_LANES as u8)
                .filter(|&lane| mask >> lane & 1 != 0)
                .collect();
            assert_eq!(&order[..lanes.len()], &lanes[..], "mask {mask:#x}");
            assert!(order[lanes.len()..].iter().all(|&lane| lane == 0));
        }
    }

    #[test]
    fn the_batch_check_sees_a_repeated_entry() {
        assert!(x_batch_coords_distinct(&[]));
        assert!(x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 4]));
        assert!(!x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 3]));
    }
}
