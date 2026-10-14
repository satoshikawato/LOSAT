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
#[inline(always)]
pub(super) fn x_test_hit(
    entry: &mut DiagStruct,
    subject_offset: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
) -> Option<(i32, i32)> {
    if entry.flag() != 0 {
        let subject_plus_offset = subject_offset + diag_offset;
        if subject_plus_offset >= entry.last_hit() {
            entry.set_last_hit(subject_plus_offset);
            entry.set_flag(0);
        }
        return None;
    }
    let last_hit = entry.last_hit() - diag_offset;
    let diff = subject_offset - last_hit;
    if diff >= window {
        entry.set_last_hit(subject_offset + diag_offset);
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
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let s_first = 0i32;
    let mut s_last =
        i32::try_from(subject_length).expect("subject length must fit Int4") - lookup.word_length;
    if s_last < s_first {
        s_last = s_first;
    }
    assert!(diag_mask >= 0 && diag_array.len() > diag_mask as usize);
    // A subject shorter than a word has one scan position whose word runs past the sentinel; the
    // reference reads that byte unchecked (one byte past the buffer for a one-residue subject).
    // It is read as a sentinel here. With one scan position no two hits share a diagonal, and the
    // next subject starts window + length further on the diagonals, so the byte reaches no output.
    let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0) as usize;
    let prefetch_distance = x_prefetch_distance();
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
        for (i, &query_offset) in hits.iter().enumerate() {
            // Read ahead the entry of a later hit of this cell: the hits of one subject offset
            // fall on distinct diagonals spread over the table, so each is a cache miss that the
            // hardware cannot predict. A prefetch reads nothing the test depends on.
            if let Some(&ahead) = hits.get(i + prefetch_distance) {
                x_prefetch(diag_array, ((ahead - s) & diag_mask) as usize);
            }
            let diag_coord = ((query_offset - s) & diag_mask) as usize;
            // SAFETY: diag_coord <= diag_mask < diag_array.len() (asserted above).
            let entry = unsafe { diag_array.get_unchecked_mut(diag_coord) };
            if let Some((last_hit, diff)) = x_test_hit(entry, s, diag_offset, window, word_size) {
                candidate(query_offset, s, last_hit, diff, entry);
            }
        }
    }
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
}
