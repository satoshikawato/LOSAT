//! EXPERIMENT (LOSAT_X_PAIRPAR, shadow LOSAT_X_PAIRPARSHADOW): the one-hit seed stage of one
//! subject chunk (scan, lookup chains, word extension, diagonal hash, ungapped extension) with
//! the scan of the subject cut into strips that the threads of the search pool scan, and the
//! diagonal-hash fold done by the owner thread in the original scan order.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
//! ```c
//!     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//!
//!         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
//!
//!         if (hitsfound == 0)
//!             continue;
//!
//!         total_hits += hitsfound;
//!         hits_extended += extend(offset_pairs, hitsfound, word_params,
//!                                 lookup_wrap, query, subject, matrix,
//!                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-838,939-941
//! ```c
//!     diag = s_off - q_off;
//!     s_end = s_off + word_length;
//!     s_off_pos = s_off + hash_table->offset;
//!     s_end_pos = s_end + hash_table->offset;
//!
//!     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
//! ...
//!     if (s_off_pos < last_hit) return 0;
//! ...
//!     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
//!                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
//!                           hit_ready, s_off_pos, window_size + Delta + 1);
//! ```
//!
//! NCBI scans the subject in pieces already: `scansub` stops when the offset array is full
//! and the loop above calls it again from `scan_range[1]`, so a scan kernel restarts at any
//! position of its grid `start + k * scan_step` and derives its state from that position.
//! This module cuts the grid of each scan range into strips (`plan_strips`). Phase 1, on
//! any thread: the reference kernel scans one strip with `s_range` of the whole range (the
//! value `scan_range[2] + lut_word_length` NCBI passes), walks the lookup chains as
//! `on_word` does, and runs the steps of `process_offset_pairs` that read nothing but the
//! inputs (the context of the query offset, the bounds test, the word extension of
//! `s_BlastNaExtend`); a pair that passes becomes a `Seed`, in the order of the pairs.
//! Phase 2, on the owner thread only: the seeds of strip 0, 1, 2, ... in that order go
//! through a statement-for-statement copy of the rest of `process_offset_pairs` for one-hit
//! mode (diagonal-hash retrieve, the `last_hit` test, `type_of_word`, the ungapped
//! extension, the cutoff, `ungapped_hits.push`, the diagonal-hash insert) on the real
//! `diag_hash` and `ungapped_hits`. The owner therefore sees the reference sequence of
//! pairs (strip order = scan order; within a strip the kernel order and the chain order),
//! and every state change is made by the same statements in the same order: the diagonal
//! hash (including the cells that `insert` takes over from stale diagonals of the same
//! bucket) and the hit list are the ones of the reference. Seeds are never reordered by
//! diagonal. The strip boundaries depend on the number of threads, the result does not:
//! a strip boundary only splits the reference sequence of pairs, which the owner rejoins
//! in order (so `-num_threads 1` and `-num_threads 8` give the same bytes).
//!
//! Scope (anything else returns `false`, and the reference code runs): the diagonal hash
//! (query longer than 8000), `window_size == 0` (one hit), `scan_range == 0`, a contiguous
//! megablast table, an unmasked subject, no small-query word extension, no debug, trace or
//! window diagnostics, one subject searched by the owner of a pool of two threads or more.

use super::{
    determine_scanning_offsets, extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi,
    mask_lookup_index, packed_kmer_at_seed_mask, scan_subject_kmers_range,
    scan_subject_kmers_range_mb_10_1, scan_subject_kmers_range_mb_10_2,
    scan_subject_kmers_range_mb_10_3, scan_subject_kmers_range_mb_11_1mod4,
    scan_subject_kmers_range_mb_11_2mod4, scan_subject_kmers_range_mb_11_3mod4,
    scan_subject_kmers_range_mb_9_1, scan_subject_kmers_range_mb_9_2,
    scan_subject_kmers_range_mb_any, select_mb_scan_kind, type_of_word, BlastnTiming,
    DiagHashTable, MbScanSubjectKind, QueryContext, QueryContextIndex, UngappedHit,
    COMPRESSION_RATIO,
};
use crate::algorithm::blastn::lookup::TwoStageLookup;
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};

/// No NCBI counterpart: reads the switches once (LOSAT_X_PAIRPAR, LOSAT_X_PAIRPARSHADOW); they
/// do not change any value NCBI computes.
fn switches() -> (bool, bool) {
    static SWITCHES: OnceLock<(bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        (
            std::env::var_os("LOSAT_X_PAIRPAR").is_some(),
            std::env::var_os("LOSAT_X_PAIRPARSHADOW").is_some(),
        )
    })
}

// No NCBI counterpart: counters of the shadow summary and of LOSAT_TIMING (scheduling only).
static SHADOW_CHUNKS: AtomicU64 = AtomicU64::new(0);
static SHADOW_HITS: AtomicU64 = AtomicU64::new(0);
static SHADOW_CELLS: AtomicU64 = AtomicU64::new(0);
static X_CHUNKS: AtomicU64 = AtomicU64::new(0);
static X_STRIPS: AtomicU64 = AtomicU64::new(0);
static X_SEEDS: AtomicU64 = AtomicU64::new(0);
static X_SKIPPED: AtomicU64 = AtomicU64::new(0);
static X_FOLD_NS: AtomicU64 = AtomicU64::new(0);
static X_OWNER_SCAN_NS: AtomicU64 = AtomicU64::new(0);
static X_WAIT_NS: AtomicU64 = AtomicU64::new(0);

/// The inputs of the seed stage of one subject chunk, borrowed from
/// `collect_prelim_hits_for_chunk` (the same names as there).
pub(super) struct Inputs<'a> {
    pub(super) two_stage: &'a TwoStageLookup,
    pub(super) search_seq_packed: &'a [u8],
    pub(super) s_len: usize,
    pub(super) subject_seq_ranges: &'a [(i32, i32)],
    pub(super) subject_masked: bool,
    pub(super) scan_step: usize,
    pub(super) query_context_index: &'a QueryContextIndex,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) encoded_query_concat_blastna: &'a [u8],
    pub(super) encoded_query_concat_blastna_with_sentinels: &'a [u8],
    pub(super) query_four_base: &'a [u8],
    pub(super) cutoff_scores: &'a [i32],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) score_matrix: &'a [i32; 256],
    pub(super) nucl_score_table: &'a [i32; 256],
    pub(super) diag_offset: isize,
    pub(super) diag_hash_window: i32,
    pub(super) use_array_indexing: bool,
    pub(super) window_size: usize,
    pub(super) scan_range: usize,
    pub(super) small_na_word: bool,
    pub(super) diagnostics: bool,
    pub(super) one_subject_owner: bool,
    pub(super) timing: Option<&'a BlastnTiming>,
}

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:135-149
/// ```c
/// typedef union BlastOffsetPair {
///     struct {
///         Uint4 q_off;  /**< Query offset */
///         Uint4 s_off;  /**< Subject offset */
///     } qs_offsets;     /**< Query/subject offset pair */
/// ```
/// A pair after the word extension of `s_BlastNaExtend` (`q_offset -= ext_left`,
/// `s_offset -= ext_left`), with its context and the `s_range` of its scan range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Seed {
    context: u32,
    q: u32,
    s: u32,
    s_range: u32,
}

/// One strip of the scan grid: the kernel runs from `start` to `end` with `s_range`.
#[derive(Clone, Copy, Debug)]
struct StripJob {
    start: usize,
    end: usize,
    s_range: usize,
}

/// No NCBI counterpart: the number of threads of the search pool this thread runs in, when
/// there is one with two threads or more (scheduling only).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn pool_threads() -> Option<usize> {
    rayon::current_thread_index()?;
    let threads = rayon::current_num_threads();
    (threads > 1).then_some(threads)
}

#[cfg(any(
    not(feature = "parallel"),
    all(target_arch = "wasm32", not(feature = "wasm-threads"))
))]
fn pool_threads() -> Option<usize> {
    None
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// ```c
///     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
///
///         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
/// ```
/// The scan ranges `scan_subject_kmers_with_ranges` visits for an unmasked subject
/// (`scan_range[1] = 0`, `scan_range[2] = subject->length - lut_word_length`), each cut into
/// strips of at most `per` grid points; the last strip of a range ends at the range's end.
/// `s_range` is the range's `scan_range[2] + lut_word_length`.
fn plan_strips(
    seq_ranges: &[(i32, i32)],
    s_len: usize,
    word_length: usize,
    lut_word_length: usize,
    scan_step: usize,
    target_strips: usize,
) -> Vec<StripJob> {
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    if seq_ranges.is_empty() || scan_step == 0 {
        return Vec::new();
    }
    let mut scan_range = [0i32, 0i32, s_len as i32 - lut_word_length as i32];
    while determine_scanning_offsets(
        seq_ranges,
        word_length as i32,
        lut_word_length as i32,
        &mut scan_range,
    ) {
        let start = scan_range[1];
        let end = scan_range[2];
        if start >= 0 && end >= 0 {
            ranges.push((start as usize, end as usize));
        }
        scan_range[1] = scan_range[2] + 1;
    }
    let grid_points = |&(start, end): &(usize, usize)| {
        if start > end {
            0
        } else {
            (end - start) / scan_step + 1
        }
    };
    let total: usize = ranges.iter().map(grid_points).sum();
    const MIN_STRIP: usize = 512;
    let per = total.div_ceil(target_strips.max(1)).max(MIN_STRIP);
    let mut jobs = Vec::new();
    for range in &ranges {
        let (start, end) = *range;
        let s_range = end.saturating_add(lut_word_length);
        let n = grid_points(range);
        if n == 0 {
            // the kernel is called and returns at once (start > end)
            jobs.push(StripJob {
                start,
                end,
                s_range,
            });
            continue;
        }
        let mut k = 0usize;
        while k < n {
            let k2 = (k + per).min(n);
            jobs.push(StripJob {
                start: start + k * scan_step,
                end: if k2 == n {
                    end
                } else {
                    start + (k2 - 1) * scan_step
                },
                s_range,
            });
            k = k2;
        }
    }
    jobs
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1504-1533
/// ```c
///    max_hits -= mb_lt->longest_chain;
/// ...
///             if (total_hits >= max_hits)
///                 break;
/// ...
///             total_hits += s_BlastMBLookupRetrieve(mb_lt,
///                 index, offset_pairs + total_hits, s_off);
/// ```
/// Phase 1 for one strip: the scan kernel `scan_subject_kmers_with_ranges` picks for this
/// lookup width and step (unmasked subject), from `job.start` to `job.end`; every lookup
/// chain hit as in `on_word` (the chain order); every pair through `seed_of_pair`.
fn scan_strip(x: &Inputs<'_>, job: &StripJob, out: &mut Vec<Seed>) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let packed = x.search_seq_packed;
    let subject_len = x.s_len;
    let scan_step = x.scan_step;
    let s_range = job.s_range;
    let (start, end) = (job.start, job.end);
    let mut on_kmer = |kmer_start: usize, current_lut_kmer: u64| {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
        // ```c
        // offset_pairs[i].qs_offsets.q_off   = q_off - 1;
        // offset_pairs[i++].qs_offsets.s_off = s_off;
        // ```
        two_stage.for_each_hit(current_lut_kmer, |q_off_1| {
            if q_off_1 != 0 {
                seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
            }
        });
    };
    let mb_scan_kind = select_mb_scan_kind(lut_word_length, scan_step, x.subject_masked);
    match mb_scan_kind {
        Some(MbScanSubjectKind::Any) => scan_subject_kmers_range_mb_any(
            packed,
            subject_len,
            lut_word_length,
            scan_step,
            x.subject_masked,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan9_1) => scan_subject_kmers_range_mb_9_1(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan9_2) => scan_subject_kmers_range_mb_9_2(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan10_1) => scan_subject_kmers_range_mb_10_1(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan10_2) => scan_subject_kmers_range_mb_10_2(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan10_3) => scan_subject_kmers_range_mb_10_3(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan11_1Mod4) => scan_subject_kmers_range_mb_11_1mod4(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan11_2Mod4) => scan_subject_kmers_range_mb_11_2mod4(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        Some(MbScanSubjectKind::Scan11_3Mod4) => scan_subject_kmers_range_mb_11_3mod4(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            &mut on_kmer,
        ),
        None => scan_subject_kmers_range(
            packed,
            subject_len,
            lut_word_length,
            scan_step,
            x.subject_masked,
            start,
            end,
            &mut on_kmer,
        ),
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1106-1144
/// ```c
///         for (; ext_left < MIN(ext_to, s_offset); ++ext_left) {
///             s_off--;
///             q--;
///             if (s_off % COMPRESSION_RATIO == 3)
///                 s--;
///             if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6)
///                 != *q)
///                 break;
///         }
///
///         /* do the right extension if the left extension did not find all
///            the bases required */
///
///         if (ext_left < ext_to) {
///             Int4 ext_right = 0;
///             s_off = s_offset + lut_word_length;
///             if (s_off + ext_to - ext_left > s_range)
///                 continue;
/// ...
///             if (ext_left + ext_right < ext_to)
///                 continue;
///         }
///
///         q_offset -= ext_left;
///         s_offset -= ext_left;
/// ```
/// Phase 1 for one pair: the statements of `process_offset_pairs` before the diagonal is
/// read (context, bounds test, word extension), which read only the inputs. A pair that
/// the reference drops with `continue` here is dropped; any other becomes a `Seed`.
#[inline(always)]
fn seed_of_pair(
    x: &Inputs<'_>,
    q_off0: usize,
    kmer_start: usize,
    s_range: usize,
    out: &mut Vec<Seed>,
) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let word_length = two_stage.word_length();
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let encoded_query_concat_blastna_with_sentinels = x.encoded_query_concat_blastna_with_sentinels;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:730-733
    // ```c
    // Int4 context = BSearchContextInfo(q_off, query_info);
    // ```
    let context_idx = x.query_context_index.context_for_offset(q_off0);
    let ctx = &x.query_contexts[context_idx];
    let q_seq = ctx.seq.as_slice();
    let q_context_start = ctx.query_offset as usize;
    let q_context_end = q_context_start + q_seq.len();

    if q_off0 + lut_word_length > q_context_end || kmer_start + lut_word_length > s_len {
        return;
    }

    let (q_ext_start, s_ext_start) = if word_length > lut_word_length {
        let ext_to = word_length - lut_word_length;
        let max_ext_left = ext_to.min(kmer_start).min(q_off0.saturating_add(1));
        let mut ext_left = 0usize;
        let mut q_left = q_off0 + 1;
        let mut s_off = kmer_start;
        let mut s_idx = s_off / COMPRESSION_RATIO;
        while ext_left < max_ext_left {
            // LOSAT's sentinel query buffer stores logical query offset `n` at `n + 1`.
            q_left -= 1;
            s_off -= 1;
            if s_off % COMPRESSION_RATIO == COMPRESSION_RATIO - 1 {
                s_idx -= 1;
            }
            // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by max_ext_left.
            let s_byte = unsafe { *search_seq_packed.get_unchecked(s_idx) };
            let s_base = ((s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6) as u8;
            // SAFETY: q_left is bounded by the outer-sentinel query buffer.
            let q_base =
                unsafe { *encoded_query_concat_blastna_with_sentinels.get_unchecked(q_left) };
            if s_base != q_base {
                break;
            }
            ext_left += 1;
        }

        let mut ext_right = 0usize;
        if ext_left < ext_to {
            let mut s_off = kmer_start + lut_word_length;
            if s_off + (ext_to - ext_left) > s_range {
                return;
            }
            let mut q_right = q_off0 + lut_word_length + 1;
            if q_right + (ext_to - ext_left) > encoded_query_concat_blastna_with_sentinels.len() {
                return;
            }
            let mut s_idx = s_off / COMPRESSION_RATIO;
            while ext_right < (ext_to - ext_left) {
                // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by s_len check above.
                let s_byte = unsafe { *search_seq_packed.get_unchecked(s_idx) };
                let s_base = ((s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6) as u8;
                // SAFETY: q_right is bounded by the outer-sentinel query buffer.
                let q_base =
                    unsafe { *encoded_query_concat_blastna_with_sentinels.get_unchecked(q_right) };
                if s_base != q_base {
                    break;
                }
                ext_right += 1;
                q_right += 1;
                s_off += 1;
                if s_off % COMPRESSION_RATIO == 0 {
                    s_idx += 1;
                }
            }
            if ext_left + ext_right < ext_to {
                return;
            }
        }
        (q_off0 - ext_left, kmer_start - ext_left)
    } else {
        (q_off0, kmer_start)
    };
    out.push(Seed {
        context: context_idx as u32,
        q: q_ext_start as u32,
        s: s_ext_start as u32,
        s_range: s_range as u32,
    });
}

/// No NCBI counterpart: counts of one fold (LOSAT_TIMING summary).
#[derive(Default)]
struct FoldCounts {
    seeds: u64,
    skipped: u64,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-838,887-895,923-933,939-941
/// ```c
///     diag = s_off - q_off;
///     s_end = s_off + word_length;
///     s_off_pos = s_off + hash_table->offset;
///     s_end_pos = s_end + hash_table->offset;
///
///     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
///
///     /* if there is no record in hashtable, we set last_hit to be a very negative number */
///     if(!rc)  last_hit = 0;
///
///     /* hit within the explored area should be rejected*/
///     if (s_off_pos < last_hit) return 0;
/// ...
///     } else if (check_masks) {
///         /* check the masks for the word */
///         if(!s_TypeOfWord(query, subject, &q_off, &s_off,
///                         query_mask, query_info, s_range,
///                         word_length, lut_word_length, lut, FALSE, &extended)) return 0;
///         /* update the right end*/
///         s_end += extended;
///         s_end_pos += extended;
///     }
/// ...
///             if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
/// ...
///                 BLAST_SaveInitialHit(init_hitlist, q_off, s_off, final_data);
///                 s_end_pos = ungapped_data->length + ungapped_data->s_start
///                           + hash_table->offset;
///             } else {
///                 hit_ready = 0;
///             }
/// ...
///     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
///                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
///                           hit_ready, s_off_pos, window_size + Delta + 1);
/// ```
/// Phase 2 (owner only): the rest of `process_offset_pairs` for the seeds of one strip, in
/// order, statement for statement for the diagonal hash in one-hit mode (`two_hits` false,
/// so `off_found` stays false and `hit_ready` true until the cutoff test); the debug, trace
/// and diagonal-array branches are outside the scope and left out.
fn fold(
    x: &Inputs<'_>,
    seeds: &[Seed],
    diag_hash: &mut DiagHashTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    counts: &mut FoldCounts,
) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let word_length = two_stage.word_length();
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let diag_offset = x.diag_offset;
    let diag_hash_window = x.diag_hash_window;
    let timing_enabled = timing.is_some();
    counts.seeds += seeds.len() as u64;
    for seed in seeds {
        let context_idx = seed.context as usize;
        let ctx = &x.query_contexts[context_idx];
        let q_idx = context_idx as u32;
        let q_context_start = ctx.query_offset as usize;
        let q_context_end = q_context_start + ctx.seq.len();
        let cutoff_score = x.cutoff_scores[context_idx];
        let q_ext_start = seed.q as usize;
        let s_ext_start = seed.s as usize;
        let s_range = seed.s_range as usize;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:828-839
        // ```c
        // diag = s_off - q_off;
        // s_end = s_off + word_length;
        // s_off_pos = s_off + hash_table->offset;
        // rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
        // if(!rc)  last_hit = 0;
        // if (s_off_pos < last_hit) return 0;
        // ```
        let diag = s_ext_start as isize - q_ext_start as isize;
        let (last_hit, _hit_saved) = {
            let (level, _hit_len, hit_saved) =
                diag_hash.retrieve(diag as i32).unwrap_or((0, 0, false));
            (level, hit_saved)
        };
        let s_off_pos = s_ext_start + diag_offset as usize;
        let s_off_pos_i32 = s_off_pos as i32;

        // Hit within explored area should be rejected
        if s_off_pos_i32 < last_hit {
            counts.skipped += 1;
            continue;
        }

        let mut q_off = q_ext_start;
        let mut s_off = s_ext_start;
        let mut s_end = s_ext_start + word_length;
        let mut s_end_pos = s_end + diag_offset as usize;
        let off_found = false;
        let mut hit_ready = true;
        let query_mask = ctx.masks.as_slice();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:41-69 (s_MBLookup)
        // ```c
        // if (! PV_TEST(pv, index, mb_lt->pv_array_bts)) {
        //     return FALSE;
        // }
        // q_off = mb_lt->hashtable[index];
        // while (q_off) {
        //     if (q_off == q_pos) return TRUE;
        //     q_off = mb_lt->next_pos[q_off];
        // }
        // ```
        let mut is_seed_masked = |s_pos: usize, q_pos: usize| -> bool {
            if s_pos + lut_word_length > s_len {
                return true;
            }
            let kmer = mask_lookup_index(
                packed_kmer_at_seed_mask(search_seq_packed, s_pos, lut_word_length),
                lut_word_length,
            );
            if !two_stage.has_hits(kmer) {
                return true;
            }
            let q_off_1 = (q_pos + 1) as u32;
            !two_stage.contains_hit(kmer, q_off_1)
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:887-895
        // ```c
        //     }  else if (check_masks) {
        // 	 /* check the masks for the word */
        //         if (!s_TypeOfWord(query, subject, &q_off, &s_off,
        //                           query_mask, query_info, s_range,
        //                           word_length, lut_word_length, lut, FALSE, &extended)) return 0;
        //         /* update the right end*/
        //         s_end += extended;
        //         s_end_pos += extended;
        //     }
        // ```
        // (`two_hits` is false; the small-query word extension is outside the scope, so the
        // locations are the query masks.)
        let (wt, ext, q_off_adj, s_off_adj) = type_of_word(
            q_off,
            s_off,
            !query_mask.is_empty(),
            q_context_end,
            s_range,
            word_length,
            lut_word_length,
            false, // check_double = FALSE (not in two-hit block)
            &mut is_seed_masked,
        );
        if wt == 0 {
            // Non-word, skip this hit
            continue;
        }
        q_off = q_off_adj;
        s_off = s_off_adj;
        s_end += ext;
        s_end_pos += ext;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
        // ```c
        // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
        //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
        //                       hit_ready, s_off_pos, window_size + Delta + 1);
        // ```
        // (`hit_ready` is still 1 in one-hit mode; the reference keeps this test.)
        if !hit_ready {
            diag_hash.insert(
                diag as i32,
                s_end_pos as i32,
                (s_end_pos - s_off_pos) as i32,
                false,
                s_off_pos as i32,
                diag_hash_window,
            );
            continue;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:908-919
        // ```c
        //             if ( word_params->options->program_number == eBlastTypeBlastn &&
        //                  (word_params->matrix_only_scoring || word_length < 11))
        //             {
        //                 s_NuclUngappedExtendExact(query, subject, matrix, q_off,
        //                                   s_off, -(cutoffs->x_dropoff), ungapped_data);
        //             }else {
        //                 s_NuclUngappedExtend(query, subject, matrix, q_off, s_end,
        //                                  s_off, -(cutoffs->x_dropoff),
        //                                  ungapped_data,
        //                                  word_params->nucl_score_table,
        //                                  cutoffs->reduced_nucl_cutoff_score);
        //             }
        // ```
        let x_dropoff = x.x_dropoff_scores[context_idx];
        let reduced_cutoff = x.reduced_cutoff_scores[context_idx];
        let ungapped_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let ungapped = if word_length < 11 {
            extend_hit_ungapped_exact_ncbi(
                x.encoded_query_concat_blastna,
                search_seq_packed,
                q_off,
                s_off,
                s_len,
                x_dropoff,
                x.score_matrix,
            )
        } else {
            extend_hit_ungapped_approx_ncbi(
                x.encoded_query_concat_blastna,
                x.query_four_base,
                search_seq_packed,
                q_off,
                s_off,
                s_end,
                s_len,
                x_dropoff,
                x.nucl_score_table,
                reduced_cutoff,
                x.score_matrix,
            )
        };
        if let Some(ungapped_start) = ungapped_start {
            let elapsed_ns = ungapped_start.elapsed().as_nanos() as u64;
            if let Some(timing) = timing {
                timing
                    .ungapped_ns
                    .fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
                timing
                    .ungapped_calls
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
        }
        debug_assert!(ungapped.q_start >= q_context_start);
        let qs = ungapped.q_start - q_context_start;
        let qe = qs + ungapped.length;
        let ss = ungapped.s_start;
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:928-929
        // ```c
        //                 s_end_pos = ungapped_data->length + ungapped_data->s_start
        //                           + hash_table->offset;
        // ```
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:923-932
        // ```c
        //             if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
        // ...
        //             } else {
        //                 hit_ready = 0;
        //             }
        // ```
        if !(off_found || ungapped_score >= cutoff_score) {
            hit_ready = false;
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
            // ```c
            // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
            //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
            //                       hit_ready, s_off_pos, window_size + Delta + 1);
            // ```
            diag_hash.insert(
                diag as i32,
                s_end_pos as i32,
                if hit_ready {
                    0
                } else {
                    (s_end_pos - s_off_pos) as i32
                },
                hit_ready,
                s_off_pos as i32,
                diag_hash_window,
            );
            continue;
        }

        // NCBI architecture: Collect ungapped hit for batch processing
        ungapped_hits.push(UngappedHit {
            context_idx: q_idx,
            query_idx: ctx.query_idx,
            query_frame: ctx.frame,
            query_context_offset: ctx.query_offset,
            seed_q_off: q_off - q_context_start,
            seed_s_off: s_off,
            qs,
            qe,
            ss,
            se: ungapped_se,
            score: ungapped_score,
        });

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
        // ```c
        // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
        //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
        //                       hit_ready, s_off_pos, window_size + Delta + 1);
        // ```
        diag_hash.insert(
            diag as i32,
            ungapped_s_end_pos as i32,
            0,
            true,
            s_off_pos as i32,
            diag_hash_window,
        );
    }
}

/// No NCBI counterpart: the strips of one scan and their results (scheduling only). Strip
/// `i` is scanned once, by whichever thread claims it; the owner takes the results in
/// strip order.
struct Strips<T> {
    cursor: AtomicUsize,
    slots: Vec<Mutex<Option<Vec<T>>>>,
    ready: Vec<AtomicBool>,
    failed: AtomicBool,
}

// No NCBI counterpart: marks the strips as failed if a scan unwinds, so that the owner stops
// waiting for its result.
struct Producing<'a, T> {
    strips: &'a Strips<T>,
    done: bool,
}

impl<T> Drop for Producing<'_, T> {
    fn drop(&mut self) {
        if !self.done {
            self.strips.failed.store(true, Ordering::Release);
        }
    }
}

impl<T> Strips<T> {
    fn new(n: usize) -> Self {
        Self {
            cursor: AtomicUsize::new(0),
            slots: (0..n).map(|_| Mutex::new(None)).collect(),
            ready: (0..n).map(|_| AtomicBool::new(false)).collect(),
            failed: AtomicBool::new(false),
        }
    }

    fn claim(&self) -> Option<usize> {
        let i = self.cursor.fetch_add(1, Ordering::Relaxed);
        (i < self.slots.len()).then_some(i)
    }

    fn produce(&self, i: usize, produce: &impl Fn(usize, &mut Vec<T>)) {
        let mut guard = Producing {
            strips: self,
            done: false,
        };
        let mut out = Vec::new();
        produce(i, &mut out);
        *self.slots[i].lock().unwrap_or_else(|e| e.into_inner()) = Some(out);
        self.ready[i].store(true, Ordering::Release);
        guard.done = true;
    }

    fn helper(&self, produce: &impl Fn(usize, &mut Vec<T>)) {
        while let Some(i) = self.claim() {
            self.produce(i, produce);
        }
    }

    /// Owner: the result of strip `i`; while it is not there, scan a later strip.
    fn take(
        &self,
        i: usize,
        produce: &impl Fn(usize, &mut Vec<T>),
        scan_ns: &mut u64,
        wait_ns: &mut u64,
        timed: bool,
    ) -> Vec<T> {
        let mut idle = 0u32;
        let mut wait_start: Option<std::time::Instant> = None;
        loop {
            if self.ready[i].load(Ordering::Acquire) {
                if let Some(t0) = wait_start {
                    *wait_ns += t0.elapsed().as_nanos() as u64;
                }
                return self.slots[i]
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .take()
                    .unwrap_or_default();
            }
            assert!(
                !self.failed.load(Ordering::Acquire),
                "LOSAT_X_PAIRPAR: a strip scan failed"
            );
            match self.claim() {
                Some(j) => {
                    idle = 0;
                    if let Some(t0) = wait_start.take() {
                        *wait_ns += t0.elapsed().as_nanos() as u64;
                    }
                    let t0 = timed.then(std::time::Instant::now);
                    self.produce(j, produce);
                    if let Some(t0) = t0 {
                        *scan_ns += t0.elapsed().as_nanos() as u64;
                    }
                }
                None => {
                    if timed && wait_start.is_none() {
                        wait_start = Some(std::time::Instant::now());
                    }
                    idle = idle.saturating_add(1);
                    if idle < 128 {
                        std::hint::spin_loop();
                    } else {
                        std::thread::yield_now();
                    }
                }
            }
        }
    }
}

/// No NCBI counterpart (scheduling only): runs `produce` for strips `0..n` on this thread
/// and `helpers` other threads of its pool, and `consume` on this thread for strips
/// `0, 1, ..., n - 1` in that order, each as soon as its result is there.
fn run_strips<T: Send>(
    n: usize,
    helpers: usize,
    produce: &(impl Fn(usize, &mut Vec<T>) + Sync),
    mut consume: impl FnMut(usize, Vec<T>),
    scan_ns: &mut u64,
    wait_ns: &mut u64,
    timed: bool,
) {
    let strips = Strips::new(n);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if helpers > 0 && rayon::current_thread_index().is_some() {
        let strips = &strips;
        rayon::in_place_scope(|scope| {
            for _ in 0..helpers.min(n) {
                scope.spawn(move |_| strips.helper(produce));
            }
            for i in 0..n {
                let out = strips.take(i, produce, scan_ns, wait_ns, timed);
                consume(i, out);
            }
        });
        return;
    }
    let _ = helpers;
    for i in 0..n {
        let out = strips.take(i, produce, scan_ns, wait_ns, timed);
        consume(i, out);
    }
}

/// No NCBI counterpart: what one seed stage did (LOSAT_TIMING summary).
#[derive(Default)]
struct StageStats {
    strips: usize,
    counts: FoldCounts,
    scan_pause_ns: u64,
    fold_ns: u64,
    owner_scan_ns: u64,
    wait_ns: u64,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// ```c
///     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
///
///         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
/// ...
///         hits_extended += extend(offset_pairs, hitsfound, word_params,
/// ```
/// The seed stage of one chunk on `threads` threads (this one and `threads - 1` helpers of
/// its pool): strips scanned by any thread, folded here in strip order into `hash` and
/// `hits`.
fn seed_stage(
    x: &Inputs<'_>,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
) -> StageStats {
    let jobs = plan_strips(
        x.subject_seq_ranges,
        x.s_len,
        x.two_stage.word_length(),
        x.two_stage.lut_word_length(),
        x.scan_step,
        8 * threads.max(4),
    );
    let produce = |i: usize, out: &mut Vec<Seed>| scan_strip(x, &jobs[i], out);
    let mut stats = StageStats {
        strips: jobs.len(),
        ..StageStats::default()
    };
    let StageStats {
        counts,
        scan_pause_ns,
        fold_ns,
        owner_scan_ns,
        wait_ns,
        ..
    } = &mut stats;
    run_strips(
        jobs.len(),
        threads.saturating_sub(1),
        &produce,
        |_, seeds| {
            let t0 = timed.then(std::time::Instant::now);
            fold(x, &seeds, hash, hits, timing, scan_pause_ns, counts);
            if let Some(t0) = t0 {
                *fold_ns += t0.elapsed().as_nanos() as u64;
            }
        },
        owner_scan_ns,
        wait_ns,
        timed,
    );
    stats
}

/// The result of the shadow computation of one chunk, compared after the reference ran.
struct Stash {
    hits: Vec<UngappedHit>,
    hash: DiagHashTable,
}

thread_local! {
    static STASH: RefCell<Option<Stash>> = const { RefCell::new(None) };
}

// No NCBI counterpart: a copy of the diagonal hash for the shadow computation.
fn copy_hash(h: &DiagHashTable) -> DiagHashTable {
    DiagHashTable {
        num_buckets: h.num_buckets,
        occupancy: h.occupancy,
        capacity: h.capacity,
        backbone: h.backbone.clone(),
        chain: h.chain.clone(),
        offset: h.offset,
        window: h.window,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// ```c
///     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
///
///         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
/// ...
///         hits_extended += extend(offset_pairs, hitsfound, word_params,
///                                 lookup_wrap, query, subject, matrix,
///                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
///     }
/// ```
/// Dispatch of `collect_prelim_hits_for_chunk`: the seed stage of this chunk (the loop above
/// with its one-hit extension) when the switch is on and the chunk is in scope; `true` means
/// it has run and the reference scan is skipped. With LOSAT_X_PAIRPARSHADOW the stage runs on
/// copies of the diagonal hash and the hit list, `false` is returned so the reference runs,
/// and `shadow_check` compares the two.
pub(super) fn run(
    x: &Inputs<'_>,
    diag_hash: &mut DiagHashTable,
    ungapped_hits: &mut Vec<UngappedHit>,
) -> bool {
    let (on, shadow) = switches();
    if !on && !shadow {
        return false;
    }
    if shadow {
        STASH.with(|s| s.borrow_mut().take());
    }
    if x.use_array_indexing
        || x.window_size != 0
        || x.scan_range != 0
        || x.two_stage.disc().is_some()
        || x.subject_masked
        || x.small_na_word
        || x.diagnostics
        || !x.one_subject_owner
    {
        return false;
    }
    let threads = match pool_threads() {
        Some(threads) => threads,
        // the shadow also runs without a pool (strips on this thread)
        None if shadow => 1,
        None => return false,
    };
    let timing = if shadow { None } else { x.timing };
    let timed = timing.is_some() || std::env::var_os("LOSAT_TIMING").is_some();
    let t_start = timed.then(std::time::Instant::now);
    let mut shadow_hash = shadow.then(|| copy_hash(diag_hash));
    let mut shadow_hits: Vec<UngappedHit> = Vec::new();
    let stats = {
        let (hash, hits): (&mut DiagHashTable, &mut Vec<UngappedHit>) = match shadow_hash.as_mut() {
            Some(hash) => (hash, &mut shadow_hits),
            None => (diag_hash, ungapped_hits),
        };
        seed_stage(x, threads, hash, hits, timing, timed)
    };
    if let (Some(t_start), Some(timing)) = (t_start, timing) {
        // `[TIMING] scan_lookup`: the time of the stage minus the ungapped extensions, as
        // the reference counts its scan (which also leaves out the last flush).
        let elapsed_ns = t_start.elapsed().as_nanos() as u64;
        let scan_ns = elapsed_ns.saturating_sub(stats.scan_pause_ns);
        timing.scan_ns.fetch_add(scan_ns, Ordering::Relaxed);
        timing.scan_calls.fetch_add(1, Ordering::Relaxed);
    }
    X_CHUNKS.fetch_add(1, Ordering::Relaxed);
    X_STRIPS.fetch_add(stats.strips as u64, Ordering::Relaxed);
    X_SEEDS.fetch_add(stats.counts.seeds, Ordering::Relaxed);
    X_SKIPPED.fetch_add(stats.counts.skipped, Ordering::Relaxed);
    X_FOLD_NS.fetch_add(stats.fold_ns, Ordering::Relaxed);
    X_OWNER_SCAN_NS.fetch_add(stats.owner_scan_ns, Ordering::Relaxed);
    X_WAIT_NS.fetch_add(stats.wait_ns, Ordering::Relaxed);
    if let Some(hash) = shadow_hash {
        STASH.with(|s| {
            *s.borrow_mut() = Some(Stash {
                hits: shadow_hits,
                hash,
            })
        });
        return false;
    }
    true
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:307-309
/// ```c
///     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
///           sizeof(BlastInitHSP), score_compare_match);
/// ```
/// LOSAT_X_PAIRPARSHADOW, before the sort above: the hit list and the diagonal hash the
/// reference built for this chunk against the ones of `run`, element by element; panics on
/// the first difference.
pub(super) fn shadow_check(ungapped_hits: &[UngappedHit], diag_hash: &DiagHashTable) {
    let Some(stash) = STASH.with(|s| s.borrow_mut().take()) else {
        return;
    };
    let key = |h: &UngappedHit| {
        (
            h.context_idx,
            h.query_idx,
            h.query_frame,
            h.query_context_offset,
            h.seed_q_off,
            h.seed_s_off,
            h.qs,
            h.qe,
            h.ss,
            h.se,
            h.score,
        )
    };
    for (i, (r, p)) in ungapped_hits.iter().zip(&stash.hits).enumerate() {
        assert!(
            key(r) == key(p),
            "LOSAT_X_PAIRPARSHADOW: ungapped hit {i} differs: reference {:?} strips {:?}",
            key(r),
            key(p)
        );
    }
    assert_eq!(
        ungapped_hits.len(),
        stash.hits.len(),
        "LOSAT_X_PAIRPARSHADOW: number of ungapped hits differs"
    );
    let h = &stash.hash;
    assert!(
        (h.num_buckets, h.occupancy, h.capacity, h.offset, h.window)
            == (
                diag_hash.num_buckets,
                diag_hash.occupancy,
                diag_hash.capacity,
                diag_hash.offset,
                diag_hash.window
            ),
        "LOSAT_X_PAIRPARSHADOW: diagonal hash header differs: reference occupancy {} capacity {} strips occupancy {} capacity {}",
        diag_hash.occupancy,
        diag_hash.capacity,
        h.occupancy,
        h.capacity
    );
    if let Some(b) = (0..diag_hash.backbone.len()).find(|&b| diag_hash.backbone[b] != h.backbone[b])
    {
        panic!(
            "LOSAT_X_PAIRPARSHADOW: diagonal hash bucket {b} differs: reference head {} strips head {}",
            diag_hash.backbone[b], h.backbone[b]
        );
    }
    let cell = |c: &super::DiagHashCell| (c.diag, c.level, c.hit_saved, c.hit_len, c.next);
    let used = diag_hash.occupancy as usize;
    for (i, (r, p)) in diag_hash.chain[..used]
        .iter()
        .zip(&h.chain[..used])
        .enumerate()
    {
        let (r, p) = (cell(r), cell(p));
        assert!(
            r == p,
            "LOSAT_X_PAIRPARSHADOW: diagonal hash cell {i} differs: reference {r:?} strips {p:?}"
        );
    }
    SHADOW_CHUNKS.fetch_add(1, Ordering::Relaxed);
    SHADOW_HITS.fetch_add(ungapped_hits.len() as u64, Ordering::Relaxed);
    SHADOW_CELLS.fetch_add(diag_hash.occupancy as u64, Ordering::Relaxed);
}

/// No NCBI counterpart: the one summary line of LOSAT_X_PAIRPARSHADOW, and with LOSAT_TIMING
/// the counters of LOSAT_X_PAIRPAR, on stderr at the end of the search.
pub(super) fn summary() {
    let (on, shadow) = switches();
    if shadow {
        let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
        let lookup = &crate::algorithm::blastn::lookup::x_mb_lookup_par::SHADOW_TABLES;
        let lookup_cells = &crate::algorithm::blastn::lookup::x_mb_lookup_par::SHADOW_CELLS;
        eprintln!(
            "[X_PAIRPARSHADOW] blastn: lookup tables equal={} (elements {}), seed chunks equal={} (ungapped hits {}, hash cells {})",
            load(lookup),
            load(lookup_cells),
            load(&SHADOW_CHUNKS),
            load(&SHADOW_HITS),
            load(&SHADOW_CELLS)
        );
    }
    if (on || shadow) && std::env::var_os("LOSAT_TIMING").is_some() {
        let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
        let s = |c: &AtomicU64| c.load(Ordering::Relaxed) as f64 / 1e9;
        eprintln!(
            "[TIMING] x_pairpar_seed: chunks={} strips={} seeds={} skipped={} fold={:.3}s owner_scan={:.3}s owner_wait={:.3}s",
            load(&X_CHUNKS),
            load(&X_STRIPS),
            load(&X_SEEDS),
            load(&X_SKIPPED),
            s(&X_FOLD_NS),
            s(&X_OWNER_SCAN_NS),
            s(&X_WAIT_NS)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::super::{
        build_blastna_matrix, build_nucl_score_table, build_query_blastna_concat_buffers,
        build_query_four_base_bytes, diag_hash_insert_window, scan_subject_kmers_with_ranges,
    };
    use super::*;
    use crate::algorithm::blastn::lookup::{build_two_stage_lookup, reverse_complement};
    use crate::core::blast_encoding::{encode_iupac_to_blastna, encode_subject_ncbi2na_packed};
    use crate::utils::dust::MaskedInterval;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    type HashState = (Vec<u32>, Vec<(i32, i32, bool, i32, u32)>, u32, u32, i32);

    fn hash_state(h: &DiagHashTable) -> HashState {
        (
            h.backbone.clone(),
            h.chain[..h.occupancy as usize]
                .iter()
                .map(|c| (c.diag, c.level, c.hit_saved, c.hit_len, c.next))
                .collect(),
            h.occupancy,
            h.capacity,
            h.offset,
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-425
    // ```c
    //     while (index) {
    //         if (table->chain[index].diag == diag) {
    // ...
    //         } else {
    //             if (s_off - table->chain[index].level > window_size) {
    //                 table->chain[index].diag = diag;
    // ```
    // A random stream of seed outcomes on the diagonal hash: retrieve, the `last_hit` test,
    // then an insert as `fold` does (below the cutoff: `hit_len` = s_end_pos - s_off_pos,
    // not saved; saved: the ungapped end, `hit_len` 0). Diagonals d + 512 k share a bucket,
    // so cells of other diagonals are taken over (the stale rule); some cells start at
    // `level = -window` (a table reset with a window, the round-1 sign-extension case) and
    // the subject offsets start below zero around them. The stream folded from strips
    // scanned on 1..8 threads must give the decisions and the table of the serial fold.
    #[derive(Clone, Copy, Debug)]
    struct Op {
        diag: i32,
        s_off_pos: i32,
        s_end_pos: i32,
        pass: bool,
        ungapped_end: i32,
    }

    fn apply(hash: &mut DiagHashTable, op: &Op, window: i32, decisions: &mut Vec<u8>) {
        let last_hit = hash.retrieve(op.diag).map_or(0, |(level, _, _)| level);
        if op.s_off_pos < last_hit {
            decisions.push(0);
            return;
        }
        if op.pass {
            decisions.push(2);
            hash.insert(op.diag, op.ungapped_end, 0, true, op.s_off_pos, window);
        } else {
            decisions.push(1);
            hash.insert(
                op.diag,
                op.s_end_pos,
                op.s_end_pos - op.s_off_pos,
                false,
                op.s_off_pos,
                window,
            );
        }
    }

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn diag_hash_fold_from_strips_equals_serial_fold() {
        let mut rng = Lcg(0xD1A6_0001);
        for case in 0..40 {
            let word_length = if case % 2 == 0 { 28 } else { 11 };
            let window = diag_hash_insert_window(0, 0, word_length);
            let table_window = 16;
            let mut start = DiagHashTable::new(table_window);
            let base = rng.below(1000) as i32 - 500;
            let diags: Vec<i32> = (0..12)
                .map(|k| {
                    if k < 8 {
                        base + 512 * k
                    } else {
                        rng.below(4000) as i32 - 2000
                    }
                })
                .collect();
            // cells left at `last_hit = -window` by a reset with a window
            for &d in diags.iter().take(3) {
                start.insert(d, -table_window, 0, false, -table_window, window);
            }
            let n = 200 + rng.below(3000) as usize;
            let mut s = -20i32;
            let ops: Vec<Op> = (0..n)
                .map(|_| {
                    s += rng.below(9) as i32 - 2;
                    let s_end = s + word_length as i32 + rng.below(4) as i32;
                    Op {
                        diag: diags[rng.below(diags.len() as u64) as usize],
                        s_off_pos: s,
                        s_end_pos: s_end,
                        pass: rng.below(3) == 0,
                        ungapped_end: s_end + rng.below(200) as i32,
                    }
                })
                .collect();
            let mut serial = copy_hash(&start);
            let mut serial_decisions = Vec::new();
            for op in &ops {
                apply(&mut serial, op, window, &mut serial_decisions);
            }
            for threads in 1..=8usize {
                let mut cuts = vec![0usize];
                while *cuts.last().unwrap() < n {
                    let next = cuts.last().unwrap() + 1 + rng.below(64) as usize;
                    cuts.push(next.min(n));
                }
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let (state, decisions) = pool.install(|| {
                    let mut hash = copy_hash(&start);
                    let mut decisions = Vec::new();
                    let produce = |i: usize, out: &mut Vec<Op>| {
                        out.extend_from_slice(&ops[cuts[i]..cuts[i + 1]]);
                    };
                    let (mut a, mut b) = (0u64, 0u64);
                    run_strips(
                        cuts.len() - 1,
                        threads - 1,
                        &produce,
                        |_, strip| {
                            for op in &strip {
                                apply(&mut hash, op, window, &mut decisions);
                            }
                        },
                        &mut a,
                        &mut b,
                        false,
                    );
                    (hash_state(&hash), decisions)
                });
                assert_eq!(decisions, serial_decisions, "case {case} threads {threads}");
                assert_eq!(state, hash_state(&serial), "case {case} threads {threads}");
            }
        }
    }

    fn random_dna(rng: &mut Lcg, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
    // ```c
    //     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
    //
    //         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
    // ```
    // The whole stage on random sequences (the subject made of mutated copies of query
    // pieces on both strands, query masks on one context) for the scan kernels of every
    // lookup width and step: the pairs of the reference wrapper
    // `scan_subject_kmers_with_ranges` folded in one go, against `seed_stage` on 1..8
    // threads (strips planned for each thread count). Hit lists and diagonal hashes must be
    // the same.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn seed_stage_on_strips_equals_one_scan() {
        let mut rng = Lcg(0x5EED_5EED);
        let configs: [(usize, usize); 14] = [
            (28, 12),
            (19, 12),
            (11, 11),
            (16, 11),
            (18, 11),
            (21, 11),
            (9, 9),
            (10, 9),
            (10, 10),
            (11, 10),
            (12, 10),
            (16, 8),
            (15, 8),
            (24, 12),
        ];
        for (case, &(word_length, lut_word_length)) in configs.iter().enumerate() {
            let scan_step = word_length - lut_word_length + 1;
            let q_len = 3000 + rng.below(3000) as usize;
            let mut query = random_dna(&mut rng, q_len);
            // a repeat inside the query (long lookup chains)
            let motif = random_dna(&mut rng, 40);
            for k in 0..20 {
                let at = (k * 97) % (q_len - 40);
                query[at..at + 40].copy_from_slice(&motif);
            }
            let minus = reverse_complement(&query);
            let mut subject = Vec::new();
            while subject.len() < 40_000 {
                if rng.below(3) == 0 {
                    let len = 50 + rng.below(400) as usize;
                    subject.extend(random_dna(&mut rng, len));
                } else {
                    let src = if rng.below(2) == 0 { &query } else { &minus };
                    let len = 30 + rng.below(1500) as usize;
                    let at = rng.below((q_len - len.min(q_len - 1)) as u64) as usize;
                    let end = (at + len).min(q_len);
                    for &b in &src[at..end] {
                        let r = rng.below(100);
                        if r < 3 {
                            subject.push(b"ACGT"[rng.below(4) as usize]);
                        } else if r < 4 {
                            // an indel
                        } else {
                            subject.push(b);
                        }
                    }
                }
            }
            let s_len = subject.len();
            let packed = encode_subject_ncbi2na_packed(&subject);
            let ctx_seqs = vec![
                encode_iupac_to_blastna(&query),
                encode_iupac_to_blastna(&minus),
            ];
            let offsets = vec![0i32, q_len as i32 + 1];
            let masks = vec![
                vec![
                    MaskedInterval {
                        start: 100,
                        end: 180,
                    },
                    MaskedInterval {
                        start: 1000,
                        end: 1013,
                    },
                ],
                Vec::new(),
            ];
            let concat_len = 2 * q_len + 1;
            let (concat, concat_sentinels) =
                build_query_blastna_concat_buffers(&ctx_seqs, &offsets, concat_len);
            let lookup = build_two_stage_lookup(
                &ctx_seqs,
                &offsets,
                word_length,
                lut_word_length,
                &masks,
                None,
                0,
                2 * q_len,
                false,
            );
            let contexts: Vec<QueryContext> = (0..2)
                .map(|c| QueryContext {
                    query_idx: 0,
                    frame: if c == 0 { 1 } else { -1 },
                    query_offset: offsets[c],
                    seq: ctx_seqs[c].clone(),
                    masks: masks[c].clone(),
                })
                .collect();
            let index = QueryContextIndex::new(&contexts);
            let four_base = build_query_four_base_bytes(&concat);
            let score_matrix = build_blastna_matrix(1, -2);
            let nucl_score_table = build_nucl_score_table(1, -2);
            let seq_ranges = [(0i32, s_len as i32)];
            let diag_offset = 37 + rng.below(1000) as isize;
            let x = Inputs {
                two_stage: &lookup,
                search_seq_packed: &packed,
                s_len,
                subject_seq_ranges: &seq_ranges,
                subject_masked: false,
                scan_step,
                query_context_index: &index,
                query_contexts: &contexts,
                encoded_query_concat_blastna: &concat,
                encoded_query_concat_blastna_with_sentinels: &concat_sentinels,
                query_four_base: &four_base,
                cutoff_scores: &[22, 22],
                x_dropoff_scores: &[16, 16],
                reduced_cutoff_scores: &[18, 18],
                score_matrix: &score_matrix,
                nucl_score_table: &nucl_score_table,
                diag_offset,
                diag_hash_window: diag_hash_insert_window(0, 0, word_length),
                use_array_indexing: false,
                window_size: 0,
                scan_range: 0,
                small_na_word: false,
                diagnostics: false,
                one_subject_owner: true,
                timing: None,
            };
            let new_hash = || {
                let mut h = DiagHashTable::new(0);
                h.offset = diag_offset as i32;
                h
            };
            // the reference wrapper visits the pairs in one scan
            let mut seeds = Vec::new();
            let kind = select_mb_scan_kind(lut_word_length, scan_step, false);
            scan_subject_kmers_with_ranges(
                &packed,
                s_len,
                word_length,
                lut_word_length,
                scan_step,
                &seq_ranges,
                false,
                kind,
                |kmer_start, s_range, kmer| {
                    lookup.for_each_hit(kmer, |q_off_1| {
                        seed_of_pair(&x, q_off_1 as usize - 1, kmer_start, s_range, &mut seeds);
                    });
                },
            );
            let mut serial_hash = new_hash();
            let mut serial_hits = Vec::new();
            let mut counts = FoldCounts::default();
            let mut pause = 0u64;
            fold(
                &x,
                &seeds,
                &mut serial_hash,
                &mut serial_hits,
                None,
                &mut pause,
                &mut counts,
            );
            assert!(
                serial_hits.len() > 5 && counts.skipped > 0,
                "case {case}: the test sequences give too few hits ({} hits, {} skipped)",
                serial_hits.len(),
                counts.skipped
            );
            let key = |h: &UngappedHit| {
                (
                    h.context_idx,
                    h.seed_q_off,
                    h.seed_s_off,
                    h.qs,
                    h.qe,
                    h.ss,
                    h.se,
                    h.score,
                )
            };
            let serial_keys: Vec<_> = serial_hits.iter().map(key).collect();
            for threads in 1..=8usize {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let (state, keys, stats_seeds) = pool.install(|| {
                    let mut hash = new_hash();
                    let mut hits = Vec::new();
                    let stats = seed_stage(&x, threads, &mut hash, &mut hits, None, false);
                    (
                        hash_state(&hash),
                        hits.iter().map(key).collect::<Vec<_>>(),
                        stats.counts.seeds,
                    )
                });
                assert_eq!(
                    stats_seeds,
                    seeds.len() as u64,
                    "case {case} threads {threads}"
                );
                assert_eq!(keys, serial_keys, "case {case} threads {threads}");
                assert_eq!(
                    state,
                    hash_state(&serial_hash),
                    "case {case} threads {threads}"
                );
            }
        }
    }
}
