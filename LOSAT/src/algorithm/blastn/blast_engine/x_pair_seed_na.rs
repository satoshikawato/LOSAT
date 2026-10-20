//! EXPERIMENT (LOSAT_X_PAIRPARNA, shadow LOSAT_X_PAIRPARNASHADOW): the one-hit seed stage of
//! one subject chunk searched with a table that is not two-stage (the direct and Na tables of
//! words 7 to 10, `-task blastn-short`), with the subject scan cut into strips that the
//! threads of the search pool scan, and the diagonal-hash fold done in the original pair
//! order (by the owner, or split by diagonal-hash bucket over the pool with
//! LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`). The scheduling (`plan_strips`, `run_strips`,
//! `part::seed_stage_part`) and the scan kernels (`scan_kernel`) are those of
//! `LOSAT_X_PAIRPAR` in the parent module.
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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
//! ```c
//!     for (; s <= s_end; s += scan_step) {
//!         num_hits = s_BlastLookupGetNumHits(lookup, index);
//!         if (num_hits == 0)
//!             continue;
//!         s_BlastLookupRetrieve(lookup,
//!                               index,
//!                               offset_pairs + total_hits,
//!                               s_off);
//!         total_hits += num_hits;
//!     }
//! ```
//!
//! The reference for these tables is the branch of `run.rs` after `if
//! two_stage_lookup_ref.is_none()` ("Original lookup method"): one call of
//! `scan_subject_kmers_with_ranges` with `word_length == lut_word_length`, and per pair of a
//! lookup hit the context of the query offset, the diagonal, the `last_hit` test,
//! `s_TypeOfWord`, the ungapped extension on the context's own query slice, the cutoff,
//! `ungapped_hits.push` and the diagonal-hash insert. Phase 1, on any thread
//! (`scan_strip_na`): the same scan kernel on one strip of a scan range with the `s_range` of
//! the whole range; every lookup hit of a word in the table's order; for each pair the context
//! of the query offset (which reads only the inputs), giving a `Seed` (no bounds test and no
//! word extension: the reference has none for these tables). Phase 2 (`fold_na`): the seeds of
//! strip 0, 1, 2, ... in that order go through a statement-for-statement copy of the rest of
//! the reference for one-hit mode with the diagonal hash, on the real `diag_hash` and
//! `ungapped_hits` (or, with LOSAT_X_PAIRPARFOLD, a partition's seeds on that partition's
//! table; the parent's argument in `x_pair_seed_part.rs` holds for this fold unchanged: it
//! reads and changes only its seed's bucket, appends to the hit list, and reads the inputs
//! otherwise). So every decision and every state change is the reference's, whatever the
//! number of threads and strips.
//!
//! Scope (anything else returns `false`, and the reference code runs): no two-stage table,
//! the diagonal hash (query longer than 8000), `window_size == 0` (one hit), `scan_range ==
//! 0`, an unmasked subject, no debug, trace or window diagnostics, one subject searched by the
//! owner of a pool of two threads or more (the shadow also runs without a pool).
//!
//! Strip size: a lookup hit is a seed here (no word extension filters it), about 100 per
//! subject offset for 30,000 24-mers against a 3 Mb genome (297.7 M pairs), so the strips are
//! the parent's smallest (`plan_strips`'s 512 offsets, about 51 k seeds or 0.8 MB on that
//! input) and the helpers may scan at most `AHEAD_PER_THREAD` strips per thread past the
//! slowest fold.

use super::super::x_diag_lv;
use super::super::{
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, mask_lookup_index,
    packed_kmer_at_seed_mask, type_of_word, BlastnTiming, DiagHashTable, QueryContext,
    QueryContextIndex, UngappedHit,
};
use super::{
    copy_hash, part, plan_strips, pool_threads, run_strips, scan_kernel, FoldCounts, PartScratch,
    Seed, StageStats, Stash, StashPath, StripJob, STASH,
};
use crate::algorithm::blastn::lookup::{NaLookupTable, PvDirectLookup};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads the switches once (LOSAT_X_PAIRPARNA, on only when set to a
/// value other than `0` for now; LOSAT_X_PAIRPARNASHADOW, off unless set to any value); they
/// do not change any value NCBI computes.
fn switches() -> (bool, bool) {
    static SWITCHES: OnceLock<(bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        (
            std::env::var_os("LOSAT_X_PAIRPARNA").is_some_and(|v| v != "0"),
            std::env::var_os("LOSAT_X_PAIRPARNASHADOW").is_some(),
        )
    })
}

/// No NCBI counterpart (scheduling only): the strips per thread asked of `plan_strips` (its
/// 512-offset minimum decides for any subject over 512 x this x threads offsets), and the
/// strips per thread that the helpers may scan past the slowest fold.
const STRIPS_PER_THREAD: usize = 512;
const AHEAD_PER_THREAD: usize = 2;

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

/// The inputs of the seed stage of one subject chunk with a table that is not two-stage,
/// borrowed from `collect_prelim_hits_for_chunk` (the same names as there; `word_length` is
/// its `safe_k`, the word length and the lookup word length of these tables).
pub(in super::super) struct NaInputs<'a> {
    pub(in super::super) pv_direct_lookup: Option<&'a PvDirectLookup>,
    pub(in super::super) na_lookup: Option<&'a NaLookupTable>,
    pub(in super::super) use_direct_lookup: bool,
    pub(in super::super) word_length: usize,
    pub(in super::super) search_seq_packed: &'a [u8],
    pub(in super::super) s_len: usize,
    pub(in super::super) subject_seq_ranges: &'a [(i32, i32)],
    pub(in super::super) subject_masked: bool,
    pub(in super::super) scan_step: usize,
    pub(in super::super) query_context_index: &'a QueryContextIndex,
    pub(in super::super) query_contexts: &'a [QueryContext],
    pub(in super::super) encoded_queries_blastna: &'a [&'a [u8]],
    pub(in super::super) cutoff_scores: &'a [i32],
    pub(in super::super) x_dropoff_scores: &'a [i32],
    pub(in super::super) reduced_cutoff_scores: &'a [i32],
    pub(in super::super) score_matrix: &'a [i32; 256],
    pub(in super::super) nucl_score_table: &'a [i32; 256],
    pub(in super::super) diag_offset: isize,
    pub(in super::super) diag_hash_window: i32,
    pub(in super::super) use_array_indexing: bool,
    pub(in super::super) window_size: usize,
    pub(in super::super) scan_range: usize,
    pub(in super::super) diagnostics: bool,
    pub(in super::super) one_subject_owner: bool,
    pub(in super::super) timing: Option<&'a BlastnTiming>,
}

impl NaInputs<'_> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:41-79
    /// ```c
    /// num_hits = s_BlastLookupGetNumHits(lookup, index);
    /// if (num_hits == 0)
    ///     continue;
    /// s_BlastLookupRetrieve(lookup,
    ///                       index,
    ///                       offset_pairs + total_hits,
    ///                       s_off);
    /// ```
    /// The reference's lookup of a word: the presence-vector direct table for word sizes up to
    /// 13, the Na table otherwise.
    #[inline(always)]
    fn hits(&self, kmer: u64) -> &[u32] {
        if self.use_direct_lookup {
            self.pv_direct_lookup
                .map(|pv_dl| pv_dl.get_hits_checked(kmer))
                .unwrap_or(&[])
        } else {
            self.na_lookup
                .map(|na_lt| na_lt.get_hits_checked(kmer))
                .unwrap_or(&[])
        }
    }
}

/// No NCBI counterpart: whether a chunk is in the scope of this module (see the module
/// comment).
fn in_scope(x: &NaInputs<'_>) -> bool {
    !x.use_array_indexing
        && x.window_size == 0
        && x.scan_range == 0
        && !x.subject_masked
        && !x.diagnostics
        && x.one_subject_owner
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:193-207
/// ```c
///     for (; s <= s_end; s += scan_step) {
///         num_hits = s_BlastLookupGetNumHits(lookup, index);
///         if (num_hits == 0)
///             continue;
///         s_BlastLookupRetrieve(lookup,
///                               index,
///                               offset_pairs + total_hits,
///                               s_off);
///         total_hits += num_hits;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:900-903
/// ```c
///         Int4 context = BSearchContextInfo(q_off, query_info);
/// ```
/// Phase 1 for one strip: the scan kernel `scan_subject_kmers_with_ranges` picks for this
/// word length and step (unmasked subject), from `job.start` to `job.end`; every lookup hit of
/// a word in the table's order, `q_off_1 == 0` skipped as the reference skips it; each pair
/// becomes a `Seed` with the context of its query offset.
fn scan_strip_na(x: &NaInputs<'_>, job: &StripJob, out: &mut Vec<Seed>) {
    let s_range = job.s_range;
    let mut on_kmer = |kmer_start: usize, current_kmer: u64| {
        for &q_off_1 in x.hits(current_kmer) {
            if q_off_1 == 0 {
                continue;
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1412
            // ```c
            // offset_pairs[i].qs_offsets.q_off = q_off - 1;
            // ```
            let q_off0 = q_off_1 as usize - 1;
            let context_idx = x.query_context_index.context_for_offset(q_off0);
            out.push(Seed {
                context: context_idx as u32,
                q: q_off0 as u32,
                s: kmer_start as u32,
                s_range: s_range as u32,
            });
        }
    };
    scan_kernel(
        x.search_seq_packed,
        x.s_len,
        x.word_length,
        x.scan_step,
        x.subject_masked,
        job.start,
        job.end,
        &mut on_kmer,
    );
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:656-672
/// ```c
///     diag = s_off + diag_table->diag_array_length - q_off;
///     real_diag = diag & diag_table->diag_mask;
///     last_hit = hit_level_array[real_diag].last_hit;
///     s_off_pos = s_off + diag_table->offset;
///     if (s_off_pos < last_hit) return 0;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-839,887-895,923-933,939-941
/// ```c
///     diag = s_off - q_off;
///     s_end = s_off + word_length;
///     s_off_pos = s_off + hash_table->offset;
///     s_end_pos = s_end + hash_table->offset;
///
///     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
/// ...
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
/// Phase 2 for the seeds of one strip, in order: the reference's statements after the context
/// for one-hit mode with the diagonal hash (`two_hits` false, so `off_found` stays false and
/// `hit_ready` true until the cutoff test), with the reference's coordinates (the query offset
/// within its context for `s_TypeOfWord` and the extension, which runs on the context's own
/// query slice); the debug, trace and diagonal-array branches are outside the scope and left
/// out.
///
/// With `PT` (LOSAT_X_PAIRPARFOLD) the seeds are one partition's seeds of one strip,
/// `diag_hash` and `ungapped_hits` are the partition's, and `out` gets the key buckets the
/// inserts change, and the time of the ungapped extensions (instead of the shared LOSAT_TIMING
/// counters); without it `out` and `pos` are not read.
#[allow(clippy::too_many_arguments)]
fn fold_na<const PT: bool>(
    x: &NaInputs<'_>,
    seeds: &[Seed],
    diag_hash: &mut DiagHashTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    counts: &mut FoldCounts,
    out: &mut part::PartOut,
    pos: &[u32],
) {
    let safe_k = x.word_length;
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let diag_offset = x.diag_offset;
    let diag_hash_window = x.diag_hash_window;
    let timing_enabled = timing.is_some();
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the table has level bounds; the bound's
    // constants of these seeds (`s_end` adds `safe_k`).
    let x_lv = diag_hash.x_lv_on();
    let x_lv_seeds = x_diag_lv::LvSeeds::new(x.window_size > 0, safe_k, x.window_size);
    counts.seeds += seeds.len() as u64;
    for (x_k, seed) in seeds.iter().enumerate() {
        // EXPERIMENT (LOSAT_X_PAIRPARFOLD): the seed's place in the stage.
        let x_key = if PT { out.key(pos, x_k) } else { 0 };
        let q_off0 = seed.q as usize;
        let kmer_start = seed.s as usize;
        let s_range = seed.s_range as usize;
        let context_idx = seed.context as usize;
        let ctx = &x.query_contexts[context_idx];
        let q_idx = context_idx as u32;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:268-269
        // ```c
        // Uint1 *q_start = query->sequence;
        // ```
        let q_seq_blastna = x.encoded_queries_blastna[context_idx];
        let q_pos_usize = q_off0 - ctx.query_offset as usize;

        // Use pre-computed cutoff score (computed once per query-subject pair)
        let cutoff_score = x.cutoff_scores[context_idx];

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:828-829
        // ```c
        // diag = s_off - q_off;
        // s_end = s_off + word_length;
        // ```
        let diag = kmer_start as isize - q_off0 as isize;

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:833-836
        // ```c
        // rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
        // if(!rc)  last_hit = 0;
        // ```
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same decisions without the walk
        // of a bucket below the bound (`s_off_pos` below).
        let (last_hit, _hit_saved) = {
            let (level, _hit_len, hit_saved) = if x_lv {
                diag_hash.retrieve_below(
                    diag as i32,
                    &x_lv_seeds,
                    kmer_start + diag_offset as usize,
                )
            } else {
                diag_hash.retrieve(diag as i32)
            }
            .unwrap_or((0, 0, false));
            (level, hit_saved)
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:656-672
        // ```c
        // s_off_pos = s_off + diag_table->offset;
        // if (s_off_pos < last_hit) return 0;
        // ```
        let s_off_pos = kmer_start + diag_offset as usize;
        let s_off_pos_i32 = s_off_pos as i32;

        // Hit within explored area should be rejected
        if s_off_pos_i32 < last_hit {
            counts.skipped += 1;
            continue;
        }

        let mut q_off = q_pos_usize;
        let mut s_off = kmer_start;
        let mut s_end = kmer_start + safe_k;
        let mut s_end_pos = s_end + diag_offset as usize;
        let off_found = false;
        let mut hit_ready = true;
        let query_mask = ctx.masks.as_slice();

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:459-489
        // ```c
        // static NCBI_INLINE Boolean s_IsSeedMasked(...)
        // {
        //     ...
        //     return !(((T_Lookup_Callback)(lookup_wrap->lookup_callback))
        //                                      (lookup_wrap, index, q_pos));
        // }
        // ```
        let mut is_seed_masked = |s_pos: usize, q_pos: usize| -> bool {
            if s_pos + safe_k > s_len {
                return true;
            }
            let kmer = mask_lookup_index(
                packed_kmer_at_seed_mask(search_seq_packed, s_pos, safe_k),
                safe_k,
            );
            let hits = x.hits(kmer);
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
            // ```c
            // /* Also add 1 to all indices, because lookup table indices count
            //    from 1. */
            // ```
            let q_off_1 = (ctx.query_offset as usize + q_pos + 1) as u32;
            !hits.iter().any(|&hit_q_off| hit_q_off == q_off_1)
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:887-895
        // ```c
        //     } else if (check_masks) {
        //         if(!s_TypeOfWord(query, subject, &q_off, &s_off,
        //                         query_mask, query_info, s_range,
        //                         word_length, lut_word_length, lut, FALSE, &extended)) return 0;
        // ```
        let (wt, ext, q_off_adj, s_off_adj) = type_of_word(
            q_off,
            s_off,
            !query_mask.is_empty(),
            ctx.len,
            s_range,
            safe_k, // word_length
            safe_k, // lut_word_length (same for non-two-stage)
            false,  // check_double = FALSE (not in two-hit block)
            &mut is_seed_masked,
        );
        if wt == 0 {
            continue;
        }
        q_off = q_off_adj;
        s_off = s_off_adj;
        s_end += ext;
        s_end_pos += ext;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:908-919
        // ```c
        //         if (word_params->ungapped_extension) {
        //             /* perform ungapped extension */
        //             s_NuclUngappedExtend(query, subject, matrix, q_off, s_off, ...
        // ```
        let x_dropoff = x.x_dropoff_scores[context_idx];
        let reduced_cutoff = x.reduced_cutoff_scores[context_idx];
        let ungapped_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let ungapped = if safe_k < 11 {
            extend_hit_ungapped_exact_ncbi(
                q_seq_blastna,
                search_seq_packed,
                q_off,
                s_off,
                s_len,
                x_dropoff,
                x.score_matrix,
            )
        } else {
            extend_hit_ungapped_approx_ncbi(
                q_seq_blastna,
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
            if PT {
                out.ungapped_ns += elapsed_ns;
                out.ungapped_calls += 1;
            } else if let Some(timing) = timing {
                timing.ungapped_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                timing.ungapped_calls.fetch_add(1, Ordering::Relaxed);
            }
            *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
        }
        let qs = ungapped.q_start;
        let qe = ungapped.q_start + ungapped.length;
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
            part::insert::<PT>(
                diag_hash,
                out,
                x_key,
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
            seed_q_off: q_off,
            seed_s_off: s_off,
            qs,
            qe,
            ss,
            se: ungapped_se,
            score: ungapped_score,
        });
        if PT {
            out.hit_keys.push(x_key);
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:939-941
        // ```c
        // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
        //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
        //                       hit_ready, s_off_pos, window_size + Delta + 1);
        // ```
        part::insert::<PT>(
            diag_hash,
            out,
            x_key,
            diag as i32,
            ungapped_s_end_pos as i32,
            0,
            true,
            s_off_pos as i32,
            diag_hash_window,
        );
    }
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
/// its pool): the scan ranges of `scan_subject_kmers_with_ranges` (an unmasked subject,
/// `word_length == lut_word_length`) cut into strips scanned by any thread, folded here in
/// strip order into `hash` and `hits`. With `fold_part` (LOSAT_X_PAIRPARFOLD,
/// `x_pair_seed_part.rs`) the fold is split by diagonal-hash bucket over the threads instead.
fn seed_stage(
    x: &NaInputs<'_>,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
    fold_part: Option<part::FoldPart<'_>>,
) -> StageStats {
    let jobs = plan_strips(
        x.subject_seq_ranges,
        x.s_len,
        x.word_length,
        x.word_length,
        x.scan_step,
        STRIPS_PER_THREAD * threads.max(4),
    );
    let produce = |i: usize, out: &mut Vec<Seed>| scan_strip_na(x, &jobs[i], out);
    if let Some(fold_part) = fold_part {
        let parts = fold_part.parts;
        let fold_one = |st: &mut part::PartState, seeds: &[Seed], pos: &[u32]| {
            let part::PartState {
                hash,
                hits,
                out,
                counts,
                scan_pause_ns,
                ..
            } = st;
            fold_na::<true>(
                x,
                seeds,
                hash,
                hits,
                timing,
                scan_pause_ns,
                counts,
                out,
                pos,
            );
        };
        let r = part::seed_stage_part(
            jobs.len(),
            threads,
            AHEAD_PER_THREAD * threads,
            hash,
            hits,
            fold_part,
            &produce,
            &fold_one,
            timed,
        );
        part::record(part::NA_ONE_HIT, parts, &r);
        return StageStats::of_part(jobs.len(), r, timing);
    }
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
        AHEAD_PER_THREAD * threads,
        &produce,
        |_, seeds| {
            let t0 = timed.then(std::time::Instant::now);
            fold_na::<false>(
                x,
                seeds,
                hash,
                hits,
                timing,
                scan_pause_ns,
                counts,
                &mut part::PartOut::default(),
                &[],
            );
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
/// Dispatch of `collect_prelim_hits_for_chunk` for a table that is not two-stage: the seed
/// stage of this chunk (the loop above with its one-hit extension) when LOSAT_X_PAIRPARNA is
/// on and the chunk is in scope; `true` means it has run and the reference scan is skipped.
/// With LOSAT_X_PAIRPARNASHADOW the stage runs on copies of the diagonal hash and the hit
/// list, `false` is returned so the reference runs, and the parent's `shadow_check` compares
/// the two. `fold_scratch` holds the allocations of LOSAT_X_PAIRPARFOLD's partition tables
/// and buffers (`x_pair_seed_part.rs`).
pub(in super::super) fn run(
    x: &NaInputs<'_>,
    diag_hash: &mut DiagHashTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    fold_scratch: &mut PartScratch,
) -> bool {
    let (on, shadow) = switches();
    if !on && !shadow {
        return false;
    }
    if shadow {
        STASH.with(|s| s.borrow_mut().take());
    }
    if !in_scope(x) {
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
    // EXPERIMENT (LOSAT_X_PAIRPARFOLD): the fold split over the pool's threads (with the shadow,
    // on the copies, which the parent's `shadow_check` compares with the reference's table and
    // hits).
    let fold_part = part::partitions(threads).map(|parts| part::FoldPart {
        scratch: fold_scratch,
        parts,
    });
    let stats = {
        let (hash, hits): (&mut DiagHashTable, &mut Vec<UngappedHit>) = match shadow_hash.as_mut() {
            Some(hash) => (hash, &mut shadow_hits),
            None => (diag_hash, ungapped_hits),
        };
        seed_stage(x, threads, hash, hits, timing, timed, fold_part)
    };
    if let (Some(t_start), Some(timing)) = (t_start, timing) {
        // `[TIMING] scan_lookup`: the time of the stage minus the ungapped extensions, as
        // the reference counts its scan. With LOSAT_X_PAIRPARFOLD only the owner's extensions
        // are taken off (the partitions' extension times summed over the threads exceed the
        // stage's wall time).
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
                path: StashPath::Na,
            })
        });
        return false;
    }
    true
}

/// No NCBI counterpart: one chunk compared equal by the parent's `shadow_check`.
pub(super) fn count_shadow(hits: usize, cells: u32) {
    SHADOW_CHUNKS.fetch_add(1, Ordering::Relaxed);
    SHADOW_HITS.fetch_add(hits as u64, Ordering::Relaxed);
    SHADOW_CELLS.fetch_add(u64::from(cells), Ordering::Relaxed);
}

/// No NCBI counterpart: the one summary line of LOSAT_X_PAIRPARNASHADOW, and with
/// LOSAT_TIMING the counters of LOSAT_X_PAIRPARNA, on stderr at the end of the search.
pub(super) fn summary() {
    let (on, shadow) = switches();
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if shadow {
        eprintln!(
            "[X_PAIRPARNASHADOW] blastn: seed chunks equal={} (ungapped hits {}, hash cells {}){}",
            load(&SHADOW_CHUNKS),
            load(&SHADOW_HITS),
            load(&SHADOW_CELLS),
            part::shadow_fields(part::NA_ONE_HIT)
        );
    }
    if (on || shadow) && std::env::var_os("LOSAT_TIMING").is_some() {
        let s = |c: &AtomicU64| c.load(Ordering::Relaxed) as f64 / 1e9;
        eprintln!(
            "[TIMING] x_pairparna_seed: chunks={} strips={} seeds={} skipped={} fold={:.3}s owner_scan={:.3}s owner_wait={:.3}s{}",
            load(&X_CHUNKS),
            load(&X_STRIPS),
            load(&X_SEEDS),
            load(&X_SKIPPED),
            s(&X_FOLD_NS),
            s(&X_OWNER_SCAN_NS),
            s(&X_WAIT_NS),
            part::timing_fields(part::NA_ONE_HIT)
        );
    }
}
