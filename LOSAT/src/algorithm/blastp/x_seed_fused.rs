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
//! Off unless `LOSAT_X_BPSEED` is set (any value). Not in the strict set.
//!
//! EXPERIMENT (LOSAT_X_BPDIAG): with the fused loop on, the hits of one subject offset (one lookup
//! cell) are tested in batches of eight instead of one after the other: for every hit of a batch
//! the entry's new value and whether the hit is a candidate are computed without branches (in
//! AVX2 registers where the CPU has AVX2, eight lanes, the entries read with one gather), every
//! entry of the batch is stored, and then the candidates go to the closure in the cell's order.
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
//! Off unless `LOSAT_X_BPDIAG` is set (any value); it needs `LOSAT_X_BPSEED` (or
//! `LOSAT_X_BPSEEDSHADOW`, which then compares the batched loop with the reference).
//! `LOSAT_X_BPDIAG_SIMD=0` keeps the scalar batch where the CPU has AVX2, and
//! `LOSAT_X_BPDIAG_GATHER=0` reads a batch's entries with eight scalar loads instead of the AVX2
//! gather (tuning and tests; the gather can be slow under the Gather Data Sampling microcode
//! mitigation). Not in the strict set. With `LOSAT_TIMING=1` each search of a query set prints
//! `[TIMING] x_bpdiag: cells=<n> batched_hits=<m> candidates=<k>` (the cells tested in batches,
//! their hits, and the hits passed to the closure).

use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::{
    BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};

// No NCBI counterpart: reads the LOSAT_X_BPSEED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSEED): off unless set to any value.
pub(crate) fn x_bp_seed() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BPSEED").is_some())
}

// No NCBI counterpart: how many hits ahead the fused loop reads the diagonal entries
// (LOSAT_X_BPSEED_PREFETCH, default 16; 0 = no prefetch). Tuning only.
fn x_prefetch_distance() -> usize {
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
fn x_prefetch(diag_array: &[DiagStruct], index: usize) {
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
/// EXPERIMENT (LOSAT_X_BPDIAG): off unless set to any value. It acts in the fused loop only.
pub(crate) fn x_bp_diag() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BPDIAG").is_some())
}

/// How the fused loop tests the hits of one lookup cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum XDiagMode {
    /// One hit after the other ([`x_test_hit`]): the fused loop without LOSAT_X_BPDIAG.
    PerHit,
    /// LOSAT_X_BPDIAG: batches of [`X_DIAG_LANES`] hits in scalar code ([`x_test_batch_scalar`]).
    Scalar,
    /// LOSAT_X_BPDIAG with AVX2: batches of eight hits in AVX2 registers
    /// ([`x_test_batch_avx2`]); a cell's last batch of fewer than eight hits is a scalar batch.
    /// `gather`: the eight entries are read with one AVX2 gather; without it
    /// (LOSAT_X_BPDIAG_GATHER=0) with eight scalar loads into the register.
    #[cfg(target_arch = "x86_64")]
    Avx2 { gather: bool },
}

// No NCBI counterpart: whether the CPU has AVX2; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
fn x_avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

// No NCBI counterpart: chooses LOSAT's code path once (LOSAT_X_BPDIAG, LOSAT_X_BPDIAG_SIMD=0 for
// the scalar batch on a CPU with AVX2, LOSAT_X_BPDIAG_GATHER=0 for scalar loads instead of the AVX2
// gather); it does not change any value NCBI computes.
fn x_diag_mode() -> XDiagMode {
    use std::sync::OnceLock;
    static MODE: OnceLock<XDiagMode> = OnceLock::new();
    *MODE.get_or_init(|| {
        if !x_bp_diag() {
            return XDiagMode::PerHit;
        }
        let simd = std::env::var_os("LOSAT_X_BPDIAG_SIMD").is_none_or(|v| v.as_os_str() != "0");
        #[cfg(target_arch = "x86_64")]
        if simd && x_avx2_available() {
            // Tuning only: a CPU with the Gather Data Sampling microcode mitigation (Xeon 6338 on
            // the workstation) may run vpgatherdd slowly.
            let gather =
                std::env::var_os("LOSAT_X_BPDIAG_GATHER").is_none_or(|v| v.as_os_str() != "0");
            return XDiagMode::Avx2 { gather };
        }
        let _ = simd;
        XDiagMode::Scalar
    })
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
            "[TIMING] x_bpdiag: cells={} batched_hits={} candidates={}",
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
fn x_test_hit(
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
    x_scan_two_hit_in(
        x_diag_mode(),
        lookup,
        subject,
        subject_length,
        &mut XDiagTable::new(diag_array, diag_mask),
        test,
        &mut candidate,
    )
    .add_to_totals();
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
    let mut counts = XDiagCounts::default();
    match mode {
        XDiagMode::PerHit => x_scan_cells(lookup, subject, subject_length, |s, hits| {
            x_cell_per_hit(hits, s, table, test, candidate)
        }),
        XDiagMode::Scalar => x_scan_cells(lookup, subject, subject_length, |s, hits| {
            let candidates = x_cell_scalar(hits, s, table, test, candidate);
            counts.add_cell(hits.len(), candidates);
        }),
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2 { gather } => {
            assert!(x_avx2_available(), "the AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                if gather {
                    x_scan_cells_avx2::<true, F>(
                        lookup,
                        subject,
                        subject_length,
                        table,
                        test,
                        candidate,
                        &mut counts,
                    )
                } else {
                    x_scan_cells_avx2::<false, F>(
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
unsafe fn x_scan_cells_avx2<const GATHER: bool, F>(
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
    x_scan_cells(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates = unsafe { x_cell_avx2::<GATHER, F>(hits, s, table, test, candidate) };
        counts.add_cell(hits.len(), candidates);
    });
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
/// on 2026-10-14 (local, human block 1 x Species4, one thread). `GATHER`: the entries are read
/// with one gather; otherwise with eight scalar loads from the lanes' indices
/// (LOSAT_X_BPDIAG_GATHER=0, for CPUs where the gather is slow), the same values.
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_test_batch_avx2<const GATHER: bool, F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
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
    // Every lane is <= mask < entries.len() (XDiagTable::new), and an entry is a
    // repr(transparent) u32.
    let old = if GATHER {
        _mm256_i32gather_epi32::<4>(table.entries.as_ptr().cast::<i32>(), coords)
    } else {
        let mut lane_old = [0u32; X_DIAG_LANES];
        for (old, &coord) in lane_old.iter_mut().zip(&lane_coords) {
            *old = table.entries.get_unchecked(coord as usize).raw_bits();
        }
        _mm256_loadu_si256(lane_old.as_ptr().cast())
    };
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
    x_call_candidates(lanes, hits, s, table, test, candidate);
    lanes.count_ones()
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
unsafe fn x_cell_avx2<const GATHER: bool, F>(
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
        candidates += x_test_batch_avx2::<GATHER, F>(batch, s, table, test, candidate);
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
            modes.push(XDiagMode::Avx2 { gather: true });
            modes.push(XDiagMode::Avx2 { gather: false });
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
                XDiagMode::Avx2 { gather } => {
                    assert!(x_avx2_available());
                    // SAFETY: the CPU has AVX2 (asserted above).
                    counted += unsafe {
                        if gather {
                            x_cell_avx2::<true, _>(hits, *s, &mut table, test, &mut candidate)
                        } else {
                            x_cell_avx2::<false, _>(hits, *s, &mut table, test, &mut candidate)
                        }
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
                for mode in cell_modes() {
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
                        let c = x_scan_two_hit_in(
                            mode,
                            &lookup,
                            &buffer,
                            subject.len(),
                            &mut XDiagTable::new(&mut entries, mask),
                            test,
                            &mut candidate,
                        );
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
                    results.push((mode, log, bits));
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

    #[test]
    fn the_batch_check_sees_a_repeated_entry() {
        assert!(x_batch_coords_distinct(&[]));
        assert!(x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 4]));
        assert!(!x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 3]));
    }
}
