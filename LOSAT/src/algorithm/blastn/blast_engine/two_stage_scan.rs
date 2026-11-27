//! The scan of a blastn subject chunk with a two-stage lookup table (`TwoStageLookup`: every table
//! of megablast and dc-megablast, and the tables of `-task blastn` whose lookup width is below the
//! word size), as `prelim_subject::prelim_chunk_hits` runs it when PAIRPAR
//! (`x_pair_seed.rs`; for a discontiguous table `x_pair_seed_disc.rs`) does not take the chunk, on
//! one thread (`scan_two_stage_lookup`): the subject's words in scan order, through the
//! discontiguous scan (`scan_subject_disc_words_with_ranges`; with LOSAT_X_DCCSR, a table laid out
//! as CSR through `x_dc_csr::Ring`), the scan of a contiguous table laid out as CSR (LOSAT_X_BNCSR,
//! one scan per form of its head, `x_scan_csr!`) or the contiguous scan
//! (`scan_subject_kmers_with_ranges`). The lookup hits of each word, in the table's order, become
//! offset pairs of a batch, which goes to the steps of the pairs (`pair_extend::extend_offset_pairs`:
//! the word extension, the diagonal, the word type, the ungapped extension and the save of the
//! initial hit) when it is full and once after the scan.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684 (the loop of
//! `BlastNaWordFinder`, :1597-1695, with the word lengths of a megablast table, :1625-1637), where
//! `scansub` is the megablast scanner that `s_MBChooseScanSubject`
//! (c++/src/algo/blast/core/blast_nascan.c:2602-2680) picks and `extend` is the variant
//! `BlastChooseNaExtend` (c++/src/algo/blast/core/na_ungapped.c:1791-1836) picks:
//! `s_BlastNaExtendDirect` for a discontiguous table or one whose width is the word size
//! (:1798-1799), else `s_BlastNaExtendAligned` or `s_BlastNaExtend` (:1052-1174). For the small-query table the
//! scanner is `s_SmallNaChooseScanSubject`'s (blast_nascan.c:1298) and the extension the one
//! `BlastChooseNaExtend` picks at na_ungapped.c:1806-1818 (`s_BlastSmallNaExtendAlignedOneByte` or
//! `s_BlastSmallNaExtend`); its word lengths are at na_ungapped.c:1617-1624
//! NCBI's scanner returns when its batch is full and the loop calls it again where it stopped;
//! this scan hands the batch to the steps of the pairs at that word and goes on (the same pairs in
//! the same batches, in the same order).
//!
//! This module is a child of `run` for the items it keeps private (the types `DiagHashTable`,
//! `OffsetPair`, `SubjectChunk` and `UngappedHit`; `scan_subject_kmers_with_ranges`,
//! `scan_subject_disc_words_with_ranges`, `x_mb_hit_prefetch`, `x_bn_pv_once`, `x_bn_pv_once_of`,
//! `x_bn_rank_once`, `x_bn_rank_once_on` and `OFFSET_ARRAY_SIZE`). Its function is the statements
//! of the branch of `prelim_chunk_hits` that were there, in their order; the strips scans of
//! PAIRPAR (`x_pair_seed.rs` `scan_strip`, `x_pair_seed_disc.rs` `scan_strip_disc`) walk the lookup
//! hits in the order `on_word` does. What the scan reads is `TwoStageScanInputs` (the batch's word
//! finder, `prelim_subject::ScanBatch` or its copy per subject, and the chunk's own values), and
//! what it writes is `TwoStageScanState`, both built once per chunk. Nothing here runs on the pool.

use super::neg_window::NegSteps;
use super::pair_extend;
use super::prelim_subject::ScanBatch;
use super::{
    scan_subject_disc_words_with_ranges, scan_subject_kmers_with_ranges, x_bn_pv_once,
    x_bn_pv_once_of, x_bn_rank_once, x_bn_rank_once_on, x_dc_csr, x_mb_hit_prefetch, DiagHashTable,
    DiagStruct, OffsetPair, SmallNaWord, SubjectChunk, UngappedHit, OFFSET_ARRAY_SIZE,
};
use crate::algorithm::blastn::lookup::{TwoStageLookup, XHeadFifo};
use std::mem::MaybeUninit;

/// What the scan of one subject chunk with a two-stage table reads: the batch's word finder
/// (`ScanBatch`, or its copy per subject: the query block, the cutoffs, the scoring, word lengths,
/// the scan step and window, the diagonal table's form and sizes, the diagnostics), its two-stage
/// table, and the chunk's own values: the packed subject from the chunk's start, the chunk's
/// length, ranges and whether the subject is masked, the chunk (its offset in the subject), the
/// diagonal table's offset for this chunk, whether the hash's retrieve takes the level bound of
/// LOSAT_X_DIAGLV, and the subject's index and trace id (the trace lines). Built once per chunk by
/// `prelim_chunk_hits`, after PAIRPAR declines the chunk.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1597-1602,1625-1634
/// (`BlastNaWordFinder`, whose parameters are pointers to the search's shared structs and the
/// chunk) and c++/src/algo/blast/core/blast_engine.c:451 (the batch's size)
/// The level bound and the diagnostics have no NCBI counterpart.
pub(super) struct TwoStageScanInputs<'a> {
    pub(super) scan: &'a ScanBatch<'a>,
    pub(super) two_stage: &'a TwoStageLookup,
    pub(super) search_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    pub(super) subject_seq_ranges: &'a [(i32, i32)],
    pub(super) subject_masked: bool,
    pub(super) chunk: &'a SubjectChunk,
    pub(super) diag_offset: isize,
    pub(super) x_lv_ref: bool,
    pub(super) s_idx: usize,
    pub(super) s_id: &'a str,
}

/// What the scan writes: the chunk's initial hits, the diagonal table (the array's levels and
/// lengths, or the hash), the batch of offset pairs, and the debug counters (the subject positions
/// scanned, the lookup hits, the seeds in the debug window, the hits under the cutoff and the hits
/// saved). Built once per chunk by `prelim_chunk_hits` with `TwoStageScanInputs`, from
/// its scratch and its counters (references only, the types the scan had as a branch of that
/// fn); its borrows end with the call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1603-1606
/// (`BlastNaWordFinder`)
/// The counters have no NCBI counterpart (diagnostics only).
pub(super) struct TwoStageScanState<'s> {
    /// The chunk's initial hits (NCBI's `init_hitlist`).
    pub(super) ungapped_hits: &'s mut Vec<UngappedHit>,
    /// The diagonal array (a query block of at most 8000 letters) or the diagonal hash.
    pub(super) hit_level_array: &'s mut Vec<DiagStruct>,
    pub(super) hit_len_array: &'s mut Vec<u8>,
    pub(super) diag_hash: &'s mut DiagHashTable,
    /// The batch of offset pairs (NCBI's `offset_pairs`, allocated for `offset_array_size` pairs).
    pub(super) offset_pairs: &'s mut Vec<OffsetPair>,
    /// The debug counters (written only when `debug_enabled`).
    pub(super) dbg_total_s_positions: &'s mut usize,
    pub(super) dbg_seeds_found: &'s mut usize,
    pub(super) dbg_window_seeds: &'s mut usize,
    pub(super) dbg_ungapped_low: &'s mut usize,
    pub(super) dbg_gapped_attempted: &'s mut usize,
}

/// The scan of one subject chunk with a two-stage table and the steps of the offset pairs it
/// collects, in the scan's order: the chunk's constants (the word and lookup word lengths, the
/// small-query table's word, the diagonal mask, the hash's insert window, two hits or one, the
/// level bound's constants), the inputs and the state of the steps of the pairs
/// (`pair_extend::PairExtendInputs`, `PairExtendState`), then the subject's words (`on_word`: the
/// lookup hits of the word, and of the second template's word, appended to the batch; the batch
/// goes to `pair_extend::extend_offset_pairs` before a word when it holds `OFFSET_ARRAY_SIZE` pairs
/// and after it when it holds `offset_array_size`) through the discontiguous, CSR or contiguous
/// scan, the scan's time (LOSAT_TIMING, without the time of the ungapped extensions,
/// `scan_pause_ns`), the last batch, and the debug `[PERF]` line. `neg` is where the batches'
/// steps go (`neg_window::NegSteps`): `()` for a window of 0 or more, the table of a window below 0
/// (`neg_window::neg_chunk_hits`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1504-1506,1523-1537
/// (`s_MBScanSubject_Any`, :1482-1626, one of the scanners that `s_MBChooseScanSubject`,
/// :2602-2680, picks for the loop of `BlastNaWordFinder`, na_ungapped.c:1673-1684), whose batch
/// ends before a word with hits once it holds `max_hits - longest_chain` pairs; the loop then hands
/// the batch to `extend`, `s_BlastNaExtend` (na_ungapped.c:1052-1174)
/// `#[inline(always)]`: the branch of `prelim_chunk_hits` that this was stays in that
/// fn's code (one copy in each of its two instances, over `UNGAPPED`). Out of line, its
/// dc-megablast scan loop (LOSAT_X_BNDCONE=0, then inlined in this fn) took another register
/// allocation and block layout: one more instruction and five executed instructions of
/// loop-alignment padding per subject byte (+0.84 % for a 3,000-letter query, callgrind, n1).
/// `on_word` and the other callbacks stay functions of their own, as they were. Their closures no
/// longer sit in the generic over `UNGAPPED`, so they, and the generic scans called with them, have
/// one instance for both instances of `prelim_chunk_hits` and run out of line: the discontiguous
/// scan's loop (13.46 M -> 13.39 M instructions for that query) and the LOSAT_X_DCCSR ring's scan
/// (26.67 M -> 27.04 M for a 600 kb x 300 kb pair).
#[inline(always)]
pub(super) fn scan_two_stage_lookup<N: NegSteps>(
    x: &TwoStageScanInputs<'_>,
    state: &mut TwoStageScanState<'_>,
    neg: N,
) {
    let &TwoStageScanInputs {
        scan,
        two_stage,
        search_seq_packed,
        s_len,
        subject_seq_ranges,
        subject_masked,
        chunk,
        diag_offset,
        x_lv_ref,
        s_idx,
        s_id,
    } = x;
    // What the scan reads of the batch's word finder (the steps of the pairs read the rest there).
    let &ScanBatch {
        scan_step,
        offset_array_size,
        small_na_compressed_query,
        encoded_query_concat_blastna,
        word_length,
        lut_word_length,
        debug_enabled,
        debug_mode,
        blastn_debug,
        timing_enabled,
        timing_ref,
        ..
    } = scan;
    let ungapped_hits = &mut *state.ungapped_hits;
    let hit_level_array = &mut *state.hit_level_array;
    let hit_len_array = &mut *state.hit_len_array;
    let diag_hash = &mut *state.diag_hash;
    let dbg_total_s_positions = &mut *state.dbg_total_s_positions;
    let dbg_seeds_found = &mut *state.dbg_seeds_found;
    let dbg_window_seeds = &mut *state.dbg_window_seeds;
    let dbg_ungapped_low = &mut *state.dbg_ungapped_low;
    let dbg_gapped_attempted = &mut *state.dbg_gapped_attempted;
    // For two-stage lookup, use lut_word_length (8) for scanning (`scan.lut_word_length` and
    // `scan.word_length` are the table's, `ScanBatch::new`).
    // NCBI's small-query table always keeps `masked_locations`, so
    // `s_TypeOfWord` checks the lookup words of every word hit (the
    // compressed query reads ambiguity codes as bases, and the lookup
    // table has no words with ambiguity codes).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:407-411
    #[expect(
        clippy::unnecessary_lazy_evaluations,
        reason = "the word is built for the small-query table only, as before the move (its query was a `Vec` there, whose borrow the lint does not count as cheap)"
    )]
    let small_na_word = (!small_na_compressed_query.is_empty()).then(|| SmallNaWord {
        compressed_query: small_na_compressed_query,
        subject_packed: search_seq_packed,
        query_length: encoded_query_concat_blastna.len(),
        word_length,
        lut_word_length,
    });

    // The diagonal table's length and mask for the steps of the pairs (0 with the hash), the
    // hash's insert window and the level bound's constants (LOSAT_X_DIAGLV) are the batch's
    // (`ScanBatch::new`).

    // DEBUG: Count loop iterations
    let mut dbg_left_ext_iters = 0usize;
    let mut dbg_right_ext_iters = 0usize;
    let mut dbg_ungapped_ext_calls = 0usize;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
    let scan_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    let mut scan_pause_ns: u64 = 0;
    let dbg_start_time = if debug_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:451-497
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1482-1534
    let offset_pairs = state.offset_pairs.as_mut_slice();
    let mut offset_pairs_len = 0usize;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:451-497
    let x_hit_prefetch = scan.x_hit_prefetch;
    // What the steps of the offset pairs of this chunk
    // (`pair_extend::extend_offset_pairs`) read besides the batch's word finder (`scan`, a
    // parameter of their own): the chunk's values, built once per chunk.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1052-1060
    let pair_inputs = pair_extend::PairExtendInputs {
        two_stage,
        small_na_word: &small_na_word,
        subject_masked,
        search_seq_packed,
        x_hit_prefetch,
        s_len,
        chunk,
        diag_offset,
        x_lv_ref,
        s_idx,
        s_id,
    };
    // What the steps of the offset pairs of this chunk write, held by the scan's
    // callbacks for the whole scan.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1058-1059
    let mut pair_state = pair_extend::PairExtendState {
        ungapped_hits: &mut *ungapped_hits,
        hit_level_array: &mut *hit_level_array,
        hit_len_array: &mut *hit_len_array,
        diag_hash: &mut *diag_hash,
        dbg_window_seeds: &mut *dbg_window_seeds,
        dbg_left_ext_iters: &mut dbg_left_ext_iters,
        dbg_right_ext_iters: &mut dbg_right_ext_iters,
        dbg_ungapped_ext_calls: &mut dbg_ungapped_ext_calls,
        dbg_ungapped_low: &mut *dbg_ungapped_low,
        dbg_gapped_attempted: &mut *dbg_gapped_attempted,
        neg,
    };

    // The presence vector and its shift, held by the scans' per-word closures by value
    // (`lookup::XPvView`), which test each word there before they call out.
    let pv_view = two_stage.x_pv_view();
    // The sides of the forms this chunk's scan takes, recorded while the run records them
    // (`ScanBatch::x_recording`, the batch's; `record_forms`).
    if scan.x_recording {
        record_forms(scan, two_stage);
    }
    // The batch's form of a contiguous table laid out as CSR (`ScanBatch::x_mb_head`, chosen
    // once per batch by `MbCsrHead::choose`), tested before `on_word` is made: the scans of
    // that form do not call it, and made before the test it was stored in each of their chunks.
    if let Some(x_head) = scan.x_mb_head {
        // EXPERIMENT (LOSAT_X_BNCSR): a contiguous table laid out as CSR has no
        // chains; its words with hits (the presence test of `for_each_hit` first)
        // go through `x_on_word_csr`, `on_word`'s statements with the cell's run
        // (the chain's offsets in the chain's order) instead of the chain walk, in
        // the same order, each `RING` words late through the ring that prefetches
        // their cells' offsets and runs (`lookup::x_dc_csr`). Every word still
        // counts for the debug summary, as in `on_word`. A batch of pairs that
        // reached `OFFSET_ARRAY_SIZE` goes before the next word with hits (`on_word`:
        // before the next word), with the same pairs, and before that word's warning,
        // so the debug output keeps its order. `on_word` keeps its two callers, so its
        // code is the one without the switch.
        // EXPERIMENT (LOSAT_X_BNCSRPACK): the form of the table's head is chosen
        // once (`x_dc_csr::XHead`), and `x_scan_csr!` expands to each form's scan,
        // so the scan of a head of one element per cell is the one without the
        // switch.
        // EXPERIMENT (LOSAT_X_BNPVONCE): a packed head's presence bits hold one bit
        // per cell, so its scan tests the word there only (`XPacked::x_present`), not
        // first in the coarser presence vector, whose bit stands for several cells;
        // for a head whose bits fit the cache (`x_bn_pv_once`). A larger packed head keeps the
        // presence vector in front and tests its own bit after it, in the closure too, so a
        // word that passes the vector without a run of its cell is not kept either (the
        // vector's test first: a word it drops reads no bit of the head).
        // EXPERIMENT (LOSAT_X_BNRANKONCE): with `rank_once` a word's run is the one of the
        // head's element the ring hands over with it (`XHead::run_at`; a packed head's word is
        // ranked once, where it enters the ring), else the cell's run read again (`XHead::run`).
        // The ring's methods are called as `XRing`'s, not as those of `Ring` itself (a head of
        // one element per cell), whose `emit` gets the word only.
        // `x_scan_csr!(head, |k| test)`: the scan of `head` whose per-word closure
        // keeps a word `k` for the ring when `test` holds (`XWordBatch`).
        let rank_once = scan.x_rank_once;
        macro_rules! x_scan_csr {
            ($head:expr, |$k:ident| $has:expr) => {{
                use x_dc_csr::XHead;
                let head = $head;
                let mut x_on_word_csr = |(kmer_start, s_range, index): (usize, usize, u64),
                                         e: usize| {
                    if offset_pairs_len >= OFFSET_ARRAY_SIZE {
                        pair_extend::extend_offset_pairs(
                            scan,
                            &pair_inputs,
                            &mut pair_state,
                            offset_pairs,
                            &mut offset_pairs_len,
                            &mut scan_pause_ns,
                            true,
                        );
                    }
                    let run: &[u32] = if rank_once {
                        head.run_at(e)
                    } else {
                        head.run(index as usize)
                    };
                    if (debug_mode || blastn_debug) && run.len() > 1000 {
                        eprintln!(
                            "[WARN] Large matches_slice: len={} for kmer at position {}",
                            run.len(),
                            kmer_start
                        );
                    }
                    for &q_off_1 in run {
                        if debug_enabled {
                            *dbg_seeds_found += 1;
                        }
                        offset_pairs[offset_pairs_len] = OffsetPair {
                            q_off: q_off_1 as usize - 1,
                            s_off: kmer_start,
                            s_range,
                        };
                        offset_pairs_len += 1;
                    }
                    if offset_pairs_len >= offset_array_size {
                        pair_extend::extend_offset_pairs(
                            scan,
                            &pair_inputs,
                            &mut pair_state,
                            offset_pairs,
                            &mut offset_pairs_len,
                            &mut scan_pause_ns,
                            true,
                        );
                    }
                };
                let mut ring = head.ring::<(usize, usize, u64)>();
                let mb_scan_kind = scan.mb_scan_kind(subject_masked);
                // The per-word closure holds what it reads by value (`move`: the debug switch,
                // the presence vector and its shift, or the packed head) and the places it
                // writes as its own references, so the scan keeps them in registers across the
                // calls instead of reading them again through references for every word.
                let mut words = XWordBatch::new(&mut ring, &mut x_on_word_csr);
                let words_ref = &mut words;
                let total_s_positions = &mut *dbg_total_s_positions;
                scan_subject_kmers_with_ranges(
                    search_seq_packed,
                    s_len,
                    word_length,
                    lut_word_length,
                    scan_step,
                    &subject_seq_ranges,
                    subject_masked,
                    mb_scan_kind,
                    move |kmer_start, s_range, current_lut_kmer| {
                        if debug_enabled {
                            *total_s_positions += 1;
                        }
                        // A word that passes `test` goes to its batch out of line (`x_keep_word`),
                        // so this closure is the test and that call and the scan inlines it: a
                        // word the test drops makes no call.
                        let $k = current_lut_kmer;
                        if $has {
                            x_keep_word(&mut *words_ref, kmer_start, s_range, current_lut_kmer);
                        }
                    },
                );
                words.finish();
                x_dc_csr::XRing::flush(&mut ring, |w| w.2 as usize, &mut x_on_word_csr);
            }};
        }
        match x_head {
            MbCsrHead::PackedOnce(packed) => {
                x_scan_csr!(packed, |k| packed.x_present(k as usize))
            }
            MbCsrHead::Packed(packed) => {
                x_scan_csr!(packed, |k| pv_view.has_hits(k)
                    && packed.x_present(k as usize))
            }
            MbCsrHead::Csr(csr) => x_scan_csr!(csr, |k| pv_view.has_hits(k)),
        }
    } else {
        let mut on_word = |kmer_start: usize,
                           s_range: usize,
                           current_lut_kmer: u64,
                           second_template_index: Option<u64>| {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
            if debug_enabled {
                *dbg_total_s_positions += 1;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
            if debug_mode || blastn_debug {
                // DEBUG: Track processing time
                let _dbg_start = std::time::Instant::now();
            }

            // Lookup using lut_word_length k-mer
            let matches_len = if debug_mode || blastn_debug {
                two_stage.count_hits(current_lut_kmer)
            } else {
                0
            };

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
            // DEBUG: Log max matches
            if (debug_mode || blastn_debug) && matches_len > 1000 {
                eprintln!(
                    "[WARN] Large matches_slice: len={} for kmer at position {}",
                    matches_len, kmer_start
                );
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1504-1533
            if offset_pairs_len >= OFFSET_ARRAY_SIZE {
                pair_extend::extend_offset_pairs(
                    scan,
                    &pair_inputs,
                    &mut pair_state,
                    offset_pairs,
                    &mut offset_pairs_len,
                    &mut scan_pause_ns,
                    true,
                );
            }

            // For each match, add to the offset_pairs buffer.
            two_stage.for_each_hit(current_lut_kmer, |q_off_1| {
                if debug_enabled {
                    *dbg_seeds_found += 1;
                }
                if q_off_1 != 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
                    offset_pairs[offset_pairs_len] = OffsetPair {
                        q_off: q_off_1 as usize - 1,
                        s_off: kmer_start,
                        s_range,
                    };
                    offset_pairs_len += 1;
                }
            });
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1456-1469
            // The second template's hits of the same subject offset follow
            // the first template's (`longest_chain` counts both chains).
            if let Some(second_template_index) = second_template_index {
                two_stage.for_each_hit2(second_template_index, |q_off_1| {
                    if debug_enabled {
                        *dbg_seeds_found += 1;
                    }
                    offset_pairs[offset_pairs_len] = OffsetPair {
                        q_off: q_off_1 as usize - 1,
                        s_off: kmer_start,
                        s_range,
                    };
                    offset_pairs_len += 1;
                });
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
            if offset_pairs_len >= offset_array_size {
                pair_extend::extend_offset_pairs(
                    scan,
                    &pair_inputs,
                    &mut pair_state,
                    offset_pairs,
                    &mut offset_pairs_len,
                    &mut scan_pause_ns,
                    true,
                );
            }
        };
        if let Some(disc) = two_stage.disc() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1684
            // EXPERIMENT (LOSAT_X_DCCSR): a table laid out as CSR has no chains; its
            // words go through `x_on_word_csr`, `on_word`'s statements for one template
            // with the cell's run (the chain's offsets in the chain's order) instead of
            // the chain walk, in the same order, each `RING` words late through the ring
            // that prefetches their cells' offsets and runs (`lookup::x_dc_csr`). `on_word`
            // keeps its two callers, so its code is the one without the switch.
            if let Some(csr) = disc.x_csr() {
                let mut x_on_word_csr = |(kmer_start, s_range, index): (usize, usize, u64)| {
                    if debug_enabled {
                        *dbg_total_s_positions += 1;
                    }
                    let run: &[u32] = if two_stage.has_hits(index) {
                        csr.run(index as usize)
                    } else {
                        &[]
                    };
                    if (debug_mode || blastn_debug) && run.len() > 1000 {
                        eprintln!(
                            "[WARN] Large matches_slice: len={} for kmer at position {}",
                            run.len(),
                            kmer_start
                        );
                    }
                    if offset_pairs_len >= OFFSET_ARRAY_SIZE {
                        pair_extend::extend_offset_pairs(
                            scan,
                            &pair_inputs,
                            &mut pair_state,
                            offset_pairs,
                            &mut offset_pairs_len,
                            &mut scan_pause_ns,
                            true,
                        );
                    }
                    for &q_off_1 in run {
                        if debug_enabled {
                            *dbg_seeds_found += 1;
                        }
                        offset_pairs[offset_pairs_len] = OffsetPair {
                            q_off: q_off_1 as usize - 1,
                            s_off: kmer_start,
                            s_range,
                        };
                        offset_pairs_len += 1;
                    }
                    if offset_pairs_len >= offset_array_size {
                        pair_extend::extend_offset_pairs(
                            scan,
                            &pair_inputs,
                            &mut pair_state,
                            offset_pairs,
                            &mut offset_pairs_len,
                            &mut scan_pause_ns,
                            true,
                        );
                    }
                };
                let mut ring: x_dc_csr::Ring<'_, (usize, usize, u64)> = x_dc_csr::Ring::new(csr);
                // A word without hits (its presence bit, which `x_on_word_csr` tests first) does
                // not enter the ring: tested here, inline in the scan, it costs neither the ring's
                // prefetches nor a place in it, and is only counted for the debug summary. The
                // words with hits leave in the same order, so the same pairs go out in the same
                // order (a batch that reached `OFFSET_ARRAY_SIZE` goes before the next word with
                // hits, with the same pairs and state). With the debug output on, every word still
                // enters, as the warning for a long run comes before that batch's window lines.
                let every_word = debug_mode || blastn_debug;
                let mut words_without_hits = 0usize;
                scan_subject_disc_words_with_ranges(
                    search_seq_packed,
                    s_len,
                    disc,
                    subject_seq_ranges,
                    subject_masked,
                    |s_off, s_range, index, index2| {
                        debug_assert!(index2.is_none());
                        if !every_word && !pv_view.has_hits(index) {
                            words_without_hits += 1;
                            return;
                        }
                        ring.push(
                            (s_off, s_range, index),
                            |w| w.2 as usize,
                            &mut x_on_word_csr,
                        );
                    },
                );
                ring.flush(&mut x_on_word_csr);
                if debug_enabled {
                    *dbg_total_s_positions += words_without_hits;
                }
            } else {
                // A word whose chains are both empty (the presence vector, which `for_each_hit` and
                // `for_each_hit2` test first) adds no pair: tested here, inline in the scan, so
                // such a word makes no call to `on_word` and is only counted for the debug summary.
                // A batch of pairs that reached `OFFSET_ARRAY_SIZE` then goes before the next word
                // with hits (`on_word`: before the next word), with the same pairs and the same
                // state, so every value is the same. With the debug output on, every word still
                // goes to `on_word`, whose warning for a long chain comes before that batch's
                // window lines. A word that calls `on_word` has its chains' heads prefetched and
                // goes through the FIFO (`lookup::XHeadFifo`) to `on_word`, in the same order, so
                // the same pairs and batches go out in the same order (with the debug output on,
                // its lines too).
                let every_word = debug_mode || blastn_debug;
                let mut words_without_hits = 0usize;
                let mut emit =
                    |(s_off, s_range, index, index2): (usize, usize, u64, Option<u64>)| {
                        on_word(s_off, s_range, index, index2)
                    };
                let mut fifo: XHeadFifo<(usize, usize, u64, Option<u64>)> = XHeadFifo::new();
                // (the closure holds the flag, the presence vector and the table by value, and the
                // rest by reference: its values stay in registers over the scan's loop instead of
                // being read again every word through references to the locals here)
                let (fifo_ref, emit_ref, without_ref) =
                    (&mut fifo, &mut emit, &mut words_without_hits);
                scan_subject_disc_words_with_ranges(
                    search_seq_packed,
                    s_len,
                    disc,
                    subject_seq_ranges,
                    subject_masked,
                    move |s_off, s_range, index, index2| {
                        if every_word
                            || pv_view.has_hits(index)
                            || index2.is_some_and(|i2| pv_view.has_hits(i2))
                        {
                            two_stage.x_prefetch_heads(index, index2);
                            fifo_ref.push((s_off, s_range, index, index2), emit_ref);
                        } else {
                            *without_ref += 1;
                        }
                    },
                );
                fifo.flush(&mut emit);
                if debug_enabled {
                    *dbg_total_s_positions += words_without_hits;
                }
            }
        } else {
            let mb_scan_kind = scan.mb_scan_kind(subject_masked);
            // A word whose chain is empty (the presence vector, which `for_each_hit` tests first)
            // adds no pair: tested here, inline in the scan, so such a word makes no call to
            // `on_word` and is only counted for the debug summary (as the disc scans above). A
            // batch of pairs that reached `OFFSET_ARRAY_SIZE` then goes before the next word with
            // hits, with the same pairs and the same state, so every value is the same. With the
            // debug output on, every word still goes to `on_word`, whose warning for a long chain
            // comes before that batch's window lines.
            let every_word = debug_mode || blastn_debug;
            let mut words_without_hits = 0usize;
            scan_subject_kmers_with_ranges(
                search_seq_packed,
                s_len,
                word_length,
                lut_word_length,
                scan_step,
                subject_seq_ranges,
                subject_masked,
                mb_scan_kind,
                |kmer_start, s_range, current_lut_kmer| {
                    if every_word || pv_view.has_hits(current_lut_kmer) {
                        on_word(kmer_start, s_range, current_lut_kmer, None);
                    } else {
                        words_without_hits += 1;
                    }
                },
            );
            if debug_enabled {
                *dbg_total_s_positions += words_without_hits;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
    if let Some(scan_start) = scan_start {
        let elapsed_ns = scan_start.elapsed().as_nanos() as u64;
        let scan_ns = elapsed_ns.saturating_sub(scan_pause_ns);
        if let Some(timing) = timing_ref {
            timing
                .scan_ns
                .fetch_add(scan_ns, std::sync::atomic::Ordering::Relaxed);
            timing
                .scan_calls
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1675-1683
    if offset_pairs_len > 0 {
        pair_extend::extend_offset_pairs(
            scan,
            &pair_inputs,
            &mut pair_state,
            offset_pairs,
            &mut offset_pairs_len,
            &mut scan_pause_ns,
            false,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
    // DEBUG: Log stats for this subject
    if debug_mode || blastn_debug {
        if let Some(dbg_start_time) = dbg_start_time {
            let elapsed = dbg_start_time.elapsed();
            eprintln!("[PERF] Subject scan took {:?}: left_ext={}, right_ext={}, ungapped={}, seeds={}, valid_pos={}",
                    elapsed, dbg_left_ext_iters, dbg_right_ext_iters, dbg_ungapped_ext_calls,
                    dbg_seeds_found, dbg_total_s_positions);
        }
    }
}

/// The head of a contiguous two-stage table laid out as CSR (LOSAT_X_BNCSR) in the form its scan
/// takes (`scan_two_stage_lookup`'s `x_scan_csr!`), chosen once per batch (`MbCsrHead::choose`,
/// for `ScanBatch::x_mb_head`), as NCBI chooses a table's scanner once: a packed head
/// (LOSAT_X_BNCSRPACK) whose presence bits are each word's only test (LOSAT_X_BNPVONCE,
/// `PackedOnce`), a packed head behind the presence vector (`Packed`), or a head of one element per
/// cell (`Csr`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1022-1026 (the scanner chosen
/// once for the search)
/// The forms of the head have no NCBI counterpart (layout only).
#[derive(Clone, Copy)]
pub(super) enum MbCsrHead<'a> {
    PackedOnce(x_dc_csr::XPacked<'a>),
    Packed(x_dc_csr::XPacked<'a>),
    Csr(&'a x_dc_csr::Csr),
}

impl<'a> MbCsrHead<'a> {
    /// The form of the batch's table `two_stage` (`None` without a two-stage table, for a
    /// discontiguous one and for a table with chains), and whether a word's run is the element
    /// that the ring hands over with it (LOSAT_X_BNRANKONCE, for a table with a CSR head: `false`
    /// otherwise, where no scan reads it). The sides are recorded by the chunks' scans
    /// (`record_forms`), as they were recorded when each chunk chose them.
    pub(super) fn choose(two_stage: Option<&'a TwoStageLookup>) -> (Option<Self>, bool) {
        let Some(csr) = two_stage
            .filter(|table| table.disc().is_none())
            .and_then(|table| table.x_mb_csr())
        else {
            return (None, false);
        };
        let head = match csr.x_packed() {
            Some(packed) if x_bn_pv_once_of(&packed) => MbCsrHead::PackedOnce(packed),
            Some(packed) => MbCsrHead::Packed(packed),
            None => MbCsrHead::Csr(csr),
        };
        (Some(head), x_bn_rank_once_on())
    }
}

/// No NCBI counterpart: the sides of the switches whose forms a chunk's scan with the two-stage
/// table `two_stage` takes (`LOSAT_X_REACH`), called by `scan_two_stage_lookup` only while the run
/// records them (`ScanBatch::x_recording`): LOSAT_X_MBHITPF; for a contiguous table
/// LOSAT_X_BNRANKONCE (a CSR head or chains) and, for a packed head, LOSAT_X_BNPVONCE; the sides
/// each chunk recorded where it chose its forms (`MbCsrHead::choose` chooses them once now).
#[cold]
#[inline(never)]
fn record_forms(scan: &ScanBatch<'_>, two_stage: &TwoStageLookup) {
    x_mb_hit_prefetch();
    if two_stage.disc().is_some() {
        return;
    }
    x_bn_rank_once(scan.x_mb_head.is_some());
    if let Some(packed) = two_stage.x_mb_csr().and_then(|csr| csr.x_packed()) {
        x_bn_pv_once(&packed);
    }
}

/// No NCBI counterpart (call layout only): how many words of a CSR scan (`x_scan_csr!`) its
/// closure gathers before they go to the head's ring together (`XWordBatch`). A power of two.
const X_WORD_BATCH: usize = 32;

/// No NCBI counterpart (scheduling only): the words of the scan of a CSR table (`x_scan_csr!`)
/// that pass its test, on their way to the head's ring. While the ring fills (no push has handed a
/// word over yet) each word goes to it at once (`push_one`), so its prefetches start when it is
/// found, as when every word went straight to the ring; from then on the words gather here, in
/// order, and go to the ring `X_WORD_BATCH` at a time (`push`), so the ring's pushes and the
/// hand-overs of the words that leave it (`x_on_word_csr`, the seeding) run as one loop instead of
/// one word between two of the scan's. BlobToolKit's and MOB-suite's megablast searches run 3-4 %
/// fewer cycles (ws node 1, one core, 10 rounds against dev 9943c4db7, each round faster). A
/// subject whose words never fill the ring (mlst's alleles: about 7 words each) keeps the old
/// schedule: the batch only waits for a full ring, as a batch from the first word made mlst's
/// searches 1.6-2.3 % slower (their prefetches came just before the flush that used them). The
/// words, their order and every value are the same: the ring gets the words that pass in the
/// scan's order, the last ones (`finish`) before it is flushed.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1523-1536
/// (`s_MBScanSubject_Any`: the presence test of each word in the scan's loop, the retrieve only for
/// a word with hits)
struct XWordBatch<'a, R, E> {
    /// `words[..n]`: the words that passed and have not gone to the ring, in the scan's order; a
    /// place no word has reached since the last push is uninitialised (never read). (Not zeroed
    /// for each subject: mlst's 225,777 alleles spent 29 M instructions, 1.2 %, zeroing them.)
    words: [MaybeUninit<(usize, usize, u64)>; X_WORD_BATCH],
    n: usize,
    /// Whether a push has handed a word over (the ring is full); until then the words go to the
    /// ring one at a time.
    ring_full: bool,
    ring: &'a mut R,
    emit: &'a mut E,
}

impl<'a, R, E> XWordBatch<'a, R, E>
where
    R: x_dc_csr::XRing<(usize, usize, u64)>,
    E: FnMut((usize, usize, u64), usize),
{
    fn new(ring: &'a mut R, emit: &'a mut E) -> Self {
        Self {
            words: [MaybeUninit::uninit(); X_WORD_BATCH],
            n: 0,
            ring_full: false,
            ring,
            emit,
        }
    }

    /// A word of the scan that passed its test: to the ring while it fills, else to the batch; a
    /// full batch goes to the ring.
    #[inline(always)]
    fn keep(&mut self, kmer_start: usize, s_range: usize, kmer: u64) {
        if !self.ring_full {
            self.push_one(kmer_start, s_range, kmer);
            return;
        }
        // (`n < X_WORD_BATCH` here: a full batch goes at once, so the mask never wraps a word)
        self.words[self.n & (X_WORD_BATCH - 1)] = MaybeUninit::new((kmer_start, s_range, kmer));
        self.n += 1;
        if self.n == X_WORD_BATCH {
            self.push();
        }
    }

    /// A word while the ring fills goes to it at once; whether the ring is full now (out of line,
    /// so `keep` stays a test and a store, and `x_keep_word` jumps here).
    #[inline(never)]
    fn push_one(&mut self, kmer_start: usize, s_range: usize, kmer: u64) {
        self.ring_full = self.ring.push(
            (kmer_start, s_range, kmer),
            |w| w.2 as usize,
            &mut *self.emit,
        );
    }

    /// The words gathered go to the ring, in order (out of line: once per `X_WORD_BATCH` words).
    #[inline(never)]
    fn push(&mut self) {
        for m in 0..self.n {
            // SAFETY: `m < n`, and `keep` wrote the places `0..n` since the last push.
            let w = unsafe { self.words[m].assume_init() };
            self.ring.push(w, |w| w.2 as usize, &mut *self.emit);
        }
        self.n = 0;
    }

    /// The words gathered after the last full batch go to the ring.
    fn finish(&mut self) {
        if self.n > 0 {
            self.push();
        }
    }
}

/// No NCBI counterpart (call layout only): a word of a CSR scan that passes its test, to its
/// batch (`XWordBatch::keep`), out of line with the word's three values as its arguments, so the
/// scan's closure is the test and this call and the scan inlines it: a word the test drops makes no
/// call (the ring's push inlined in the closure made it too large for the scan to inline, and the
/// batch's store inlined there cost mlst's scan two more reads per word, 1.3 % more instructions).
#[inline(never)]
fn x_keep_word<R, E>(words: &mut XWordBatch<'_, R, E>, kmer_start: usize, s_range: usize, kmer: u64)
where
    R: x_dc_csr::XRing<(usize, usize, u64)>,
    E: FnMut((usize, usize, u64), usize),
{
    words.keep(kmer_start, s_range, kmer);
}

#[cfg(test)]
#[path = "two_stage_scan_tests.rs"]
mod tests;
