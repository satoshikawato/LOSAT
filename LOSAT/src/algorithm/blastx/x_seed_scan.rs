//! The BLASTX seed loop reads the scan buffer directly, and the buffer is kept per thread.
//!
//! NCBI allocates the offset-pair array once per search thread and its word finder reads the
//! pairs straight from that array:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
//! ```c
//!     aux_struct->offset_pairs =
//!       (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-511
//! ```c
//!     while (scan_range[1] <= scan_range[2]) {
//!         /* scan the subject sequence for hits */
//!         hits = scansub(lookup_wrap, subject,
//!                                   offset_pairs, array_size, scan_range);
//!
//!         totalhits += hits;
//!         /* for each hit, */
//!         for (i = 0; i < hits; ++i) {
//!             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ```
//!
//! The reference (`Lookup::scan` + `word_finder`, which the diagnostic stages run) allocates and
//! zero-fills a buffer of `capacity()` pairs per (chunk, subject) and copies every batch into a
//! new `(i32, i32)` vector before the diagonal loop reads it. `x_word_finder_scan` is
//! `word_finder` with the scan loop of `Lookup::scan` inlined: the same `scan_range` state
//! machine, the same `array_size` (`capacity()`), the scanner writing into a per-thread buffer
//! (grown when a table needs more, never shrunk), and the diagonal statements reading the pairs in
//! place.
//!
//! Identity argument: the scanners (`s_blast_aa_scan_subject_one_range`, the compressed
//! `scan_subject`) only write the buffer and return the number of pairs written, and only
//! `pairs[0..hits)` is read afterwards, so stale contents beyond `hits` are never read. The
//! pairs are visited in the same order with the same values (`u32 as i32`, as the copy did),
//! and the diagonal statements below are those of `word_finder`. The caller uses this
//! path only when nothing traces the seeds (not in the diagnostic stages).
//!
//! EXPERIMENT (LOSAT_X_BXSCANBF): with this switch (a switch of the strict set: on unless it is
//! `0`) the scan of the standard lookup table is
//! tblastn's branch-free scan (`tblastn/x_scan_bf.rs`, `x_scan_one_range_bf`, LOSAT_X_TBNSCANBF's
//! loop): the presence bits of a block of words advance a candidate list instead of a branch, and
//! the listed words' hits are copied in position order with NCBI's stop and resume when the offset
//! array is full. Its pairs
//! and `scan_range` values are those of `s_blast_aa_scan_subject_one_range`, so the loop above
//! reads the same pairs; it writes up to `X_SCAN_BF_SLACK` pairs past `array_size`, which the
//! per-thread buffer then holds and nothing reads. With LOSAT_X_BXDIAG the fused loop takes the
//! same candidate loop (`x_seed_diag.rs`). The compressed lookup table keeps its scan.
//!
//! EXPERIMENT (LOSAT_X_BXHITBF): the direct loop takes the two outcomes most two-hit tests end in,
//! too far and overlap, without a branch (tblastn's LOSAT_X_TBNHITBF and TBLASTX's
//! LOSAT_X_TBXHITBF do the same in their loops):
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-551
//! ```c
//!                 if (diff >= window) {
//!                     /* We are beyond the window for this diagonal; start a
//!                        new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//!                 }
//!
//!                 /* If the difference is less than the wordsize (i.e. last
//!                    hit and this hit overlap), give up */
//!
//!                 if (diff < wordsize) {
//!                     continue;
//!                 }
//! ```
//! Where the direct loop scans every subject (`LOSAT_X_BXDIAG=0`, and the default without AVX2:
//! aarch64, WebAssembly) 79-81 % of the hits of the three blastx cases of v4-A10 probe bx1 are too
//! far and 14.5-17.4 % overlap, which of the two is hard to predict, and 1.2-1.4 % meet a flag-1
//! cell and 2.1-2.8 % are candidates for the extension.
//!
//! The test is blastp's (`blastp::x_seed_fused::XHitBf`, whose module header gives the arithmetic):
//! it reads the packed cell `raw` (`last_hit` in bits 0-30, the flag in bit 31) and works with
//! twice the values, `d2 = 2 (s + diag_offset) - (raw << 1)` in wrapping `u32` arithmetic read as
//! `i32`, which is twice the reference's `diff = s - (last_hit - diag_offset)` for a flag-0 cell
//! whose `last_hit` is in `[-window, diag_offset + length)` where `XHitBf::fits` holds (checked
//! once per subject): `0 <= diag_offset` and `diag_offset + length + 2 window < 2^30` (`length` the
//! subject's letters). This table holds no other flag-0 value: the cleared value, 0 in a new table,
//! or an `s + diag_offset` of this subject or of an earlier one. Every store with flag 0 is such an
//! `s + diag_offset` (the extension's `s_last_off - (wordsize - 1) + diag_offset` is stored with
//! flag 1, and the flag goes back to 0 only with a store of `s + diag_offset`); the offset only
//! grows between clears (`Diagonals::finish_subject`), so an earlier subject's stores are below
//! this subject's offset, and none of them wrapped in the 31-bit field, as they were below 2^30.
//! Then too far (`diff >= window`) is `d2 >= 2 window`, and a candidate (`wordsize <= diff <
//! window`) is `(d2 - 2 wordsize) as u32 < 2 (window - wordsize)`, one unsigned compare. The rule
//! also asks `1 <= wordsize < window`, a two-hit search (the one-hit search, window 0, keeps the
//! statements). A flag-1 cell (`raw as i32 < 0`) or a candidate goes on to the statements of
//! `word_finder`, unchanged; any other hit stores `s + diag_offset` (too far; the flag stays
//! 0) or what it read (overlap), selected without a branch. So the table holds after each hit what
//! the statements leave there, the statements see the same table for every hit they run, and the
//! init HSPs are the same, in the same order.
//!
//! The switch is read once; it is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_BXHITBF=0`. The loop with the test is its own
//! function (`x_word_finder_direct_bf`, never inlined), so the loop without it keeps its code. The
//! fused loop of LOSAT_X_BXDIAG keeps its own tests.

use super::*;
use crate::algorithm::blastp::x_seed_fused::XHitBf;
use crate::algorithm::tblastn::x_scan_bf::{x_scan_one_range_bf, X_SCAN_BF_SLACK};
use std::cell::RefCell;

// No NCBI counterpart: reads the LOSAT_X_BXSCANBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXSCANBF): whether the seed scans (the direct loop and BXDIAG's fused loop)
/// find the words with hits by the branch-free candidate loop (module header); a switch of the
/// strict set, on unless it is `0`.
pub(crate) fn x_bxscan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXSCANBF"))
}

// No NCBI counterpart: reads the LOSAT_X_BXHITBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXHITBF): whether the direct loop takes the too-far and overlap outcomes
/// without a branch where [`XHitBf::fits`] holds (module header); a switch of the strict set, on
/// unless it is `0`.
pub(crate) fn x_bxhitbf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXHITBF"))
}

/// The `[X_SHADOW]` summaries of the BXDIAG, BXSUBJFLOW, BXDBSTREAM and BXEVALEMPTY switches
/// (printed by `main` at exit, one line per switch that ran in shadow mode).
// No NCBI counterpart: prints the shadow-mode counters; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    super::x_seed_diag::x_print_shadow_summary();
    #[cfg(feature = "parallel")]
    super::super::search::x_subject_flow::x_print_shadow_summary();
    crate::algorithm::blastx::x_db_stream::x_print_shadow_summary();
    crate::algorithm::blastx::statistics::x_print_eval_empty_summary();
}

thread_local! {
    /// The per-thread offset-pair array (blast_engine.c:1040-1041).
    static X_PAIRS: RefCell<Vec<OffsetPair>> = const { RefCell::new(Vec::new()) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
// ...
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
// ...
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
// ```
/// `word_finder` without the per-call buffer and the batch copy (the seed loop of a search that
/// traces nothing).
pub(crate) fn x_word_finder_scan(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
) -> Result<Vec<InitHsp>> {
    // EXPERIMENT (LOSAT_X_BXDIAG): the fused loop where it acts (x_seed_diag.rs); where it does
    // not, the direct loop.
    if super::x_seed_diag::x_bxdiag_mode() != 0 {
        if let Some(hits) = super::x_seed_diag::x_word_finder_diag(
            batch, parameters, options, subject, lookup, diagonals,
        )? {
            return Ok(hits);
        }
    }
    x_word_finder_direct(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        x_bxscan_bf(),
        x_bxhitbf(),
    )
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's hits and table.
/// EXPERIMENT (LOSAT_X_BXHITBF): whether the direct loop for `subject` is the branch-free one:
/// `hit_bf` (the switch) and [`XHitBf::fits`] at the table's offset before the subject.
fn x_hitbf_acts(
    hit_bf: bool,
    options: &ResolvedOptions,
    subject: &[u8],
    diagonals: &Diagonals,
) -> bool {
    hit_bf
        && XHitBf::fits(
            diagonals.offset,
            subject.len().saturating_sub(1),
            options.window_size,
            options.word_size,
        )
}

// No NCBI counterpart: chooses the instance of LOSAT's direct loop (`x_word_finder_direct_with`, which
// ports the NCBI code); both give the reference's hits and table.
// The direct loop: `Lookup::scan` and `word_finder` in one function (`x_word_finder_direct_with`).
// EXPERIMENT (LOSAT_X_BXSCANBF): `branch_free` scans the standard table with `x_scan_one_range_bf`.
// EXPERIMENT (LOSAT_X_BXHITBF): with `hit_bf`, where `XHitBf::fits` holds for the subject, the
// instance with the branch-free two-hit test, `x_word_finder_direct_bf` (module header).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_direct(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    branch_free: bool,
    hit_bf: bool,
) -> Result<Vec<InitHsp>> {
    let acts = x_hitbf_acts(hit_bf, options, subject, diagonals);
    if hit_bf {
        crate::utils::xreach::x_reach!("LOSAT_X_BXHITBF", acts);
    }
    if acts {
        return x_word_finder_direct_bf(
            batch,
            parameters,
            options,
            subject,
            lookup,
            diagonals,
            branch_free,
        );
    }
    x_word_finder_direct_with::<false>(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        branch_free,
    )
}

// No NCBI counterpart: an instance of `x_word_finder_direct_with`, kept out of line; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXHITBF): the direct loop with the branch-free two-hit test, in a function
/// of its own so that the loop without the test keeps its code (module header).
#[inline(never)]
fn x_word_finder_direct_bf(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    branch_free: bool,
) -> Result<Vec<InitHsp>> {
    x_word_finder_direct_with::<true>(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        branch_free,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-511
// ```c
//         wordsize = lookup->word_length;
//     }
//
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
// ...
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ```
// `Lookup::scan` (its range set-up and loop) and `word_finder` (the statements per pair)
// in one function; the pairs are read from the scanner's buffer. EXPERIMENT (LOSAT_X_BXSCANBF):
// `branch_free` scans the standard table with `x_scan_one_range_bf` (module header).
// EXPERIMENT (LOSAT_X_BXHITBF): with `HIT_BF` (only where the caller found `XHitBf::fits` to hold)
// each hit first goes through the branch-free two-hit test (`XHitBf::test`), and only a flag-1 cell
// or a candidate reaches the statements of `word_finder`.
#[inline(always)]
fn x_word_finder_direct_with<const HIT_BF: bool>(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    branch_free: bool,
) -> Result<Vec<InitHsp>> {
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let mask = diagonals.entries.len() as i32 - 1;
    let mut hits = Vec::new();
    // Lookup::scan's set-up, with length = subject.len() - 1 as word_finder passes it.
    let length = subject.len().saturating_sub(1);
    let window = options.window_size;
    let lookup_word = match lookup {
        Lookup::Standard(l) => l.word_length,
        Lookup::Compressed(l) => l.word_length,
    };
    let ranges = [(0, length as i32)];
    let mut range = [0, 0, length as i32 - lookup_word];
    if window > 0 && range[2] < range[1] {
        range[2] = range[1];
    }
    // The scanners read `scanned`: the subject, or a padded copy when its one scan position's
    // word runs past the sentinel (seed.rs `x_scan_bytes`); the statements below read `subject`.
    let mut pad = Vec::new();
    let scanned = x_scan_bytes(subject, range[2], lookup_word, &mut pad);
    let capacity = lookup.capacity();
    // EXPERIMENT (LOSAT_X_BXSCANBF): room for the pairs the branch-free scan writes past
    // `array_size` (never read).
    let room = capacity + if branch_free { X_SCAN_BF_SLACK } else { 0 };
    // EXPERIMENT (LOSAT_X_BXHITBF): the constants of the branch-free test (the caller checked
    // XHitBf::fits).
    let bf = XHitBf::in_range(diagonals.offset, window, word);
    X_PAIRS.with(|buffer| {
        let mut buffer = buffer.borrow_mut();
        if buffer.len() < room {
            buffer.resize(room, OffsetPair::default());
        }
        let pairs = &mut buffer[..room];
        let array_size = capacity as i32;
        while range[1] <= range[2] {
            let n = match lookup {
                Lookup::Standard(l) if branch_free => {
                    x_scan_one_range_bf(l, scanned, pairs, array_size, &mut range)
                }
                Lookup::Standard(l) => {
                    s_blast_aa_scan_subject_one_range(l, scanned, pairs, array_size, &mut range)
                }
                Lookup::Compressed(l) => {
                    l.scan_subject(scanned, &ranges, pairs, array_size, &mut range)
                }
            };
            for pair in &pairs[..n as usize] {
                let (q, s) = (pair.q_off as i32, pair.s_off as i32);
                // EXPERIMENT (LOSAT_X_BXHITBF): a flag-0 cell's too-far and overlap outcomes
                // without a branch; a flag-1 cell or a candidate goes on to the statements below.
                if HIT_BF {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516
                    // ```c
                    //             diag_coord = (query_offset - subject_offset) & diag_mask;
                    // ```
                    let index = ((q - s) & mask) as usize;
                    debug_assert!(index < diagonals.entries.len());
                    // SAFETY: 0 <= index <= mask = entries.len() - 1 (a power of two less one,
                    // `Diagonals::new`).
                    let cell = unsafe { diagonals.entries.get_unchecked_mut(index) };
                    if bf.test(cell, bf.here(s)) {
                        continue;
                    }
                }
                // From here to the end of the loop body: the statements of word_finder.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,783-784
                // ```c
                //                 curr_context = BSearchContextInfo(query_offset, query_info);
                // ...
                //                 Int4 curr_context = BSearchContextInfo(query_offset,
                //                                                        query_info);
                // ```
                let context_of = |q: i32| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                };
                let context;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516,776
                // ```c
                //             diag_coord = (query_offset - subject_offset) & diag_mask;
                // ...
                //             diag_coord = (subject_offset - query_offset) & diag_mask;
                // ```
                let index = if options.window_size == 0 {
                    (s - q) & mask
                } else {
                    (q - s) & mask
                } as usize;
                let entry = &mut diagonals.entries[index];
                let (u, end, extended) = if options.window_size == 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:777-781
                    // ```c
                    //             diff = subject_offset -
                    //                 (diag_array[diag_coord].last_hit - diag_offset);
                    //
                    //             /* do an extension, but only if we have not already extended this
                    //                far */
                    // ```
                    if s.wrapping_sub(entry.last_hit().wrapping_sub(diagonals.offset)) < 0 {
                        continue;
                    }
                    context = context_of(q);
                    let p = &parameters[context];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:787-790
                    // ```c
                    //                 score = s_BlastAaExtendOneHit(matrix, subject, query,
                    //                                               subject_offset, query_offset,
                    //                                               cutoffs->x_dropoff,
                    //                                               &hsp_q, &hsp_s, &hsp_len,
                    // ```
                    let Some(result) = extend_one_hit_blosum62(
                        query,
                        subject,
                        q as usize,
                        s as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, true)
                } else {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551
                    // ```c
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
                    if entry.flag() != 0 {
                        if s.wrapping_add(diagonals.offset) < entry.last_hit() {
                            continue;
                        }
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        entry.set_flag(0);
                        continue;
                    }
                    let previous = entry.last_hit().wrapping_sub(diagonals.offset);
                    let diff = s.wrapping_sub(previous);
                    if diff >= options.window_size {
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        continue;
                    }
                    if diff < word {
                        continue;
                    }
                    context = context_of(q);
                    let p = &parameters[context];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
                    // ```c
                    //                 if (query_offset - diff <
                    //                     query_info->contexts[curr_context].query_offset) {
                    //
                    //                     /* there was no last hit for this diagnol; start a new hit */
                    //                     diag_array[diag_coord].last_hit =
                    //                         subject_offset + diag_offset;
                    //                     continue;
                    //                 }
                    // ```
                    if q - diff < batch.contexts[context].offset as i32 {
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        continue;
                    }
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                    // ```c
                    //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
                    //                                               last_hit + wordsize,
                    //                                               subject_offset, query_offset,
                    //                                               cutoffs->x_dropoff,
                    //                                               &hsp_q, &hsp_s,
                    //                                               &hsp_len, use_pssm,
                    //                                               wordsize, &right_extend,
                    // ```
                    let Some(result) = extend_two_hit_blosum62(
                        query,
                        subject,
                        (previous + word) as usize,
                        s as usize,
                        q as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, result.right_extend)
                };
                let p = &parameters[context];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-606
                // ```c
                //                 if (score >= cutoffs->cutoff_score)
                //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                //                                      query_offset, subject_offset, hsp_len,
                //                                      score);
                //
                //                 /* If an extension to the right happened, reset the last hit
                //                    so that future hits to this diagonal must start over. */
                //
                //                 if (right_extend) {
                //                     diag_array[diag_coord].flag = 1;
                //                     diag_array[diag_coord].last_hit =
                //                         s_last_off - (wordsize - 1) + diag_offset;
                //                 }
                //                 /* Otherwise, make the present hit into the previous hit for
                //                    this diagonal */
                //                 else {
                //                     diag_array[diag_coord].last_hit =
                //                         subject_offset + diag_offset;
                //                 }
                // ```
                if u.score >= p.word_cutoff {
                    hits.push(InitHsp {
                        q_seed: q,
                        s_seed: s,
                        q_start: u.q_start,
                        s_start: u.s_start,
                        length: u.length,
                        score: u.score,
                    });
                }
                if options.window_size == 0 {
                    entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
                } else if extended {
                    entry.set_flag(1);
                    entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
                } else {
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                }
            }
        }
    });
    // The rest of word_finder, unchanged.
    diagonals.finish_subject(subject.len().saturating_sub(1))?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // The frozen pinned-NCBI INIT rows of the session F fixture (one-hit and two-hit windows,
    // fresh, retained and reset diagonal tables): the direct loop gives the NCBI hits and the
    // same diagonal table as the reference loop.
    /// An entry with `last_hit` and the flag set.
    fn entry(last_hit: i32) -> DiagStruct {
        let mut entry = DiagStruct::default();
        entry.set_last_hit(last_hit);
        entry.set_flag(1);
        entry
    }

    /// The direct loop's arms: (the branch-free scan, the branch-free two-hit test).
    const ARMS: [(bool, bool); 4] = [(false, false), (true, false), (false, true), (true, true)];

    #[test]
    fn direct_scan_loop_matches_ncbi_hits_and_reference_diagonal_table() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let q = parse_fasta(&std::fs::read(dir.join("query.fna")).unwrap(), false, false).unwrap();
        let s = parse_fasta(
            &std::fs::read(dir.join("subject.faa")).unwrap(),
            true,
            false,
        )
        .unwrap();
        let mut subject = s[0]
            .sequence
            .iter()
            .copied()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect::<Vec<_>>();
        subject.push(0);
        for window in [0, 1, 39, 40, 41] {
            let args = [
                "losat",
                "blastx",
                "-query",
                "query",
                "-subject",
                "subject",
                "-seg",
                "no",
                "-comp_based_stats",
                "0",
                "-window_size",
                &window.to_string(),
            ];
            let crate::cli::Commands::Blastx(args) =
                crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                    .unwrap()
                    .command
            else {
                unreachable!()
            };
            let options = args.resolve().unwrap();
            let batch = prepare_queries(&q, &options).unwrap();
            let mut params = parameters::score_block(&batch);
            parameters::effective_lengths(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                None,
            )
            .unwrap();
            let lookup = Lookup::new(&batch, &params, &options).unwrap();
            parameters::subject_parameters(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                s[0].sequence.len(),
                None,
            )
            .unwrap();
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            let expected =
                std::fs::read_to_string(dir.join(format!("window{window}.stage"))).unwrap();
            let initial = expected
                .lines()
                .filter(|x| x.starts_with("INIT\t"))
                .map(|x| {
                    let f = x.split('\t').collect::<Vec<_>>();
                    InitHsp {
                        q_seed: f[3].parse().unwrap(),
                        s_seed: f[4].parse().unwrap(),
                        q_start: f[5].parse().unwrap(),
                        s_start: f[6].parse().unwrap(),
                        length: f[7].parse().unwrap(),
                        score: f[8].parse().unwrap(),
                    }
                })
                .collect::<Vec<_>>();
            let mut reference = Diagonals::new(length, window).unwrap();
            // EXPERIMENT (LOSAT_X_BXSCANBF, LOSAT_X_BXHITBF): the direct loop with and without the
            // branch-free scan and the branch-free two-hit test, each on its own table, compared
            // alike.
            let mut tables = ARMS.map(|_| Diagonals::new(length, window).unwrap());
            for state in 0..4 {
                for d in tables.iter_mut().chain([&mut reference]) {
                    if state == 2 {
                        d.offset = 536_870_910;
                        d.entries.fill(entry(123));
                    }
                    if state == 3 {
                        d.offset = 536_870_911;
                        d.entries.fill(entry(536_870_900));
                        d.finish_subject(83).unwrap();
                    }
                }
                // The branch-free test acts on every two-hit window above the word size here.
                assert_eq!(
                    x_hitbf_acts(true, &options, &subject, &reference),
                    window > options.word_size,
                    "window={window} state={state}"
                );
                let plain = word_finder(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap();
                assert_eq!(plain, initial, "window={window} state={state}");
                for (&(branch_free, hit_bf), d) in ARMS.iter().zip(&mut tables) {
                    let hits = x_word_finder_direct(
                        &batch,
                        &params,
                        &options,
                        &subject,
                        &lookup,
                        d,
                        branch_free,
                        hit_bf,
                    )
                    .unwrap();
                    assert_eq!(
                        hits, plain,
                        "window={window} state={state} scan_bf={branch_free} hit_bf={hit_bf}"
                    );
                    assert!(
                        d.offset == reference.offset && d.entries == reference.entries,
                        "window={window} state={state} scan_bf={branch_free} hit_bf={hit_bf}"
                    );
                }
            }
            // EXPERIMENT (LOSAT_X_BXHITBF): the edge of the branch-free test's range. At the most
            // offset `XHitBf::fits` allows for this subject (cleared cells) the test acts, and one
            // more it does not; the loops agree at both.
            if window > options.word_size {
                let most = (1 << 30) - 1 - (subject.len() as i32 - 1) - 2 * window;
                for (offset, acts) in [(most, true), (most + 1, false)] {
                    for d in tables.iter_mut().chain([&mut reference]) {
                        d.offset = offset;
                        d.entries.fill(DiagStruct::clear(window));
                    }
                    assert_eq!(
                        x_hitbf_acts(true, &options, &subject, &reference),
                        acts,
                        "window={window} offset={offset}"
                    );
                    let plain = word_finder(
                        &batch,
                        &params,
                        &options,
                        &subject,
                        &lookup,
                        &mut reference,
                        |_| {},
                    )
                    .unwrap();
                    // Only differences of offsets count below 2^30: the NCBI rows of a new table.
                    assert_eq!(plain, initial, "window={window} offset={offset}");
                    for (&(branch_free, hit_bf), d) in ARMS.iter().zip(&mut tables) {
                        let hits = x_word_finder_direct(
                            &batch,
                            &params,
                            &options,
                            &subject,
                            &lookup,
                            d,
                            branch_free,
                            hit_bf,
                        )
                        .unwrap();
                        let what = format!(
                            "window={window} offset={offset} scan_bf={branch_free} hit_bf={hit_bf}"
                        );
                        assert_eq!(hits, plain, "{what}");
                        assert!(
                            d.offset == reference.offset && d.entries == reference.entries,
                            "{what}"
                        );
                    }
                }
            }
            // A subject of one residue (shorter than a word less one letter) and one of two: with a
            // two-hit window the one scan position reads a word from position 0, which runs past
            // the sentinel of the one-residue subject. The loops scan it from a padded copy
            // (seed.rs `x_scan_bytes`; the slices here hold exactly the residues and the
            // sentinel, so an unchecked read past them trips the debug checks), and they agree.
            for n in [1, 2] {
                let short = subject[..n].iter().copied().chain([0]).collect::<Vec<u8>>();
                let mut reference = Diagonals::new(length, window).unwrap();
                let plain = word_finder(
                    &batch,
                    &params,
                    &options,
                    &short,
                    &lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap();
                assert!(plain.is_empty(), "window={window} {n} residues: {plain:?}");
                for (branch_free, hit_bf) in ARMS {
                    let mut d = Diagonals::new(length, window).unwrap();
                    let hits = x_word_finder_direct(
                        &batch,
                        &params,
                        &options,
                        &short,
                        &lookup,
                        &mut d,
                        branch_free,
                        hit_bf,
                    )
                    .unwrap();
                    assert_eq!(
                        hits, plain,
                        "window={window} {n} residues, scan_bf={branch_free} hit_bf={hit_bf}"
                    );
                    assert!(
                        d.offset == reference.offset && d.entries == reference.entries,
                        "window={window} {n} residues, scan_bf={branch_free} hit_bf={hit_bf}"
                    );
                }
            }
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
            // ```c
            // typedef struct DiagStruct {
            //    signed int last_hit   : 31; /**< Offset of the last hit */
            //    unsigned int flag      : 1 ; /**< Reset the next extension? */
            // } DiagStruct;
            // ```
            // An offset just below 2^30 (NCBI reaches such offsets at a -window_size near 2^29,
            // from the second subject on): the stores past 2^30 wrap in the 31-bit field, and the
            // loops read the wrapped values back alike. tests/blastx_diag_bitfield.rs compares the
            // rows with NCBI 2.17.0 at such a window. (The offset is past INT4_MAX / 4, so the
            // tables are cleared at the end of the subject.) The offset is outside the range of
            // the branch-free test, so its arms run the loop without it.
            for d in tables.iter_mut().chain([&mut reference]) {
                d.offset = (1 << 30) - 64;
                d.entries.fill(DiagStruct::clear(window));
            }
            assert!(!x_hitbf_acts(true, &options, &subject, &reference));
            let plain = word_finder(
                &batch,
                &params,
                &options,
                &subject,
                &lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            for (&(branch_free, hit_bf), d) in ARMS.iter().zip(&mut tables) {
                let hits = x_word_finder_direct(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    d,
                    branch_free,
                    hit_bf,
                )
                .unwrap();
                assert_eq!(
                    hits, plain,
                    "window={window} wrapped, scan_bf={branch_free} hit_bf={hit_bf}"
                );
                assert!(
                    d.offset == reference.offset && d.entries == reference.entries,
                    "window={window} wrapped, scan_bf={branch_free} hit_bf={hit_bf}"
                );
            }
            if !initial.is_empty() {
                assert_ne!(
                    plain, initial,
                    "window={window}: the wrapped stores change the seeds"
                );
            }
        }
    }

    /// The statements of `word_finder` up to the two `continue`s of a flag-0 cell
    /// (aa_ungapped.c:519-551): the cell after a hit that is too far or overlaps, or `None` for a
    /// flag-1 cell or a candidate for the extension.
    fn two_hit_reference(
        cell: DiagStruct,
        s: i32,
        offset: i32,
        window: i32,
        word: i32,
    ) -> Option<DiagStruct> {
        if cell.flag() != 0 {
            return None;
        }
        let diff = s.wrapping_sub(cell.last_hit().wrapping_sub(offset));
        if diff >= window {
            let mut next = cell;
            next.set_last_hit(s.wrapping_add(offset));
            return Some(next);
        }
        (diff < word).then_some(cell)
    }

    /// `XHitBf::test` on a copy of `cell`: the cell after it and whether it took the hit.
    fn bf_outcome(
        cell: DiagStruct,
        s: i32,
        offset: i32,
        window: i32,
        word: i32,
    ) -> (DiagStruct, bool) {
        let mut next = cell;
        let bf = XHitBf::in_range(offset, window, word);
        let taken = bf.test(&mut next, bf.here(s));
        (next, taken)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551 (module header).
    // EXPERIMENT (LOSAT_X_BXHITBF): the branch-free test against the reference statements at the
    // edges of its range: offsets from 0 to the most `XHitBf::fits` allows for the subject's length
    // (and that it allows no more), word sizes from 1 to window - 1, windows from 2 to 2^28 - 1,
    // subject positions at both ends, and cells around the outcomes' edges (diff = wordsize - 1,
    // wordsize, window - 1, window), the cleared cell (`-window`), a new table's 0, the lowest and
    // highest `last_hit` a flag-0 cell can hold, and flag-1 cells. Past the range a cleared cell
    // can be read wrong, which is why the range is checked.
    #[test]
    fn the_branch_free_test_gives_the_reference_outcomes_in_its_range() {
        // too far, overlap, candidate, flag 1
        let mut outcomes = [0u64; 4];
        for window in [2, 3, 4, 40, 41, 1000, (1 << 28) - 1] {
            let mut words = vec![1, 2, 3, window / 2, window - 1];
            words.retain(|&w| 1 <= w && w < window);
            words.sort_unstable();
            words.dedup();
            for word in words {
                for length in [1usize, 2, 100, 1 << 20] {
                    let most = (1i64 << 30) - 1 - length as i64 - 2 * i64::from(window);
                    let most = i32::try_from(most).unwrap();
                    assert!(XHitBf::fits(most, length, window, word));
                    assert!(!XHitBf::fits(most + 1, length, window, word));
                    assert!(!XHitBf::fits(-1, length, window, word));
                    for offset in [0, most / 2, most] {
                        let high = offset + length as i32;
                        for s in [0, 1, length as i32 / 2, length as i32 - 1] {
                            if !(0..length as i32).contains(&s) {
                                continue;
                            }
                            let here = s + offset;
                            let mut last_hits = vec![-window, 0, high - 1, offset, here + 1];
                            for d in [word - 1, word, window - 1, window, window + 1, 0] {
                                last_hits.push(here - d);
                            }
                            for last_hit in last_hits {
                                for flag in [0, 1] {
                                    // A flag-0 cell holds a last_hit in [-window, high).
                                    if flag == 0 && !(-window..high).contains(&last_hit) {
                                        continue;
                                    }
                                    let mut cell = DiagStruct::default();
                                    cell.set_last_hit(last_hit);
                                    cell.set_flag(flag);
                                    let (next, taken) = bf_outcome(cell, s, offset, window, word);
                                    let what = format!(
                                        "window={window} word={word} length={length} offset={offset} s={s} last_hit={last_hit} flag={flag}"
                                    );
                                    match two_hit_reference(cell, s, offset, window, word) {
                                        Some(expected) => {
                                            assert!(taken, "{what}");
                                            assert_eq!(next, expected, "{what}");
                                            outcomes[usize::from(next == cell)] += 1;
                                        }
                                        None => {
                                            assert!(!taken, "{what}");
                                            assert_eq!(next, cell, "{what}");
                                            outcomes[2 + flag as usize] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        assert!(outcomes.iter().all(|&n| n > 0), "{outcomes:?}");
        assert!(!XHitBf::fits(0, 10, 3, 3));
        assert!(!XHitBf::fits(0, 10, 3, 0));
        assert!(!XHitBf::fits(0, 10, 0, 3));
        assert!(XHitBf::fits(0, 10, 3, 2));
        // Just past the range: a cleared cell at the end of a subject is too far for the
        // reference (diff = 2^30 + 499), but twice that wraps in an i32 and the branch-free test
        // reads it as an overlap.
        let (window, word, length) = (1000, 3, 1000usize);
        let offset = (1 << 30) - 1500;
        assert!(!XHitBf::fits(offset, length, window, word));
        let cleared = DiagStruct::clear(window);
        let expected = two_hit_reference(cleared, 999, offset, window, word).unwrap();
        assert_ne!(expected, cleared);
        assert_ne!(
            bf_outcome(cleared, 999, offset, window, word),
            (expected, true)
        );
    }
}
