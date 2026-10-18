//! EXPERIMENT (LOSAT_X_PAIRPARDISC, shadow LOSAT_X_PAIRPARDISCSHADOW): the two-hit seed
//! stage of discontiguous megablast (`-task dc-megablast`) for one subject chunk with the
//! diagonal hash, with the subject scan cut into strips that the threads of the search pool
//! scan, and the diagonal-hash fold done by the owner thread in the original pair order.
//! The scheduling (`plan_strips`, `run_strips`) and the per-pair word step (`seed_of_pair`)
//! are those of `LOSAT_X_PAIRPAR` in the parent module.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1624-1634,1673-1684
//! ```c
//!     else if (lookup_wrap->lut_type == eMBLookupTable) {
//!         BlastMBLookupTable *lookup =
//!                                 (BlastMBLookupTable *) lookup_wrap->lut;
//!         if (lookup->discontiguous) {
//!             word_length = lookup->template_length;
//!             lut_word_length = lookup->template_length;
//! ...
//!     while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
//!
//!         hitsfound = scansub(lookup_wrap, subject, offset_pairs, max_hits, &scan_range[1]);
//! ...
//!         hits_extended += extend(offset_pairs, hitsfound, word_params,
//!                                 lookup_wrap, query, subject, matrix,
//!                                 query_info, ewp, init_hitlist, scan_range[2] + lut_word_length);
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1651-1655
//! ```c
//!     if (subject->mask_type != eNoSubjMasking) {
//!         if (lookup_wrap->lut_type == eMBLookupTable &&
//!             ((BlastMBLookupTable *) lookup_wrap->lut)->discontiguous) {
//!             /* discontiguous scan subs assumes any (non-aligned starting offset */
//! ```
//! The discontiguous scan routines (`s_MB_DiscWordScanSubject_*`, stride 1) start at any
//! subject offset and derive their accumulator from it; NCBI restarts them at
//! `scan_range[1]` whenever its offset array is full. Phase 1, on any thread: the reference
//! kernel `disc_word_scan_subject` scans one strip of a scan range with the `s_range` of
//! the whole range (`scan_range[2] + lut_word_length`), every word's chains go out in
//! NCBI's order (per subject offset the first template's chain, then the second's,
//! `MB_ACCESS_HITS2`), and each pair goes through `seed_of_pair` (context of the query
//! offset and the bounds test; a discontiguous word has no word extension, since
//! `word_length == lut_word_length`). Phase 2, on the owner thread only: the seeds of strip
//! 0, 1, 2, ... in that order go through `fold_two_hit`, a statement-for-statement copy of
//! the two-hit diagonal-hash branch of `process_offset_pairs` (`run.rs`): the retrieve with
//! `last_hit` and `hit_saved`, the `last_hit` test, the window test, `s_TypeOfWord` with and
//! without `check_double`, the off-diagonal search, the ungapped extension, the cutoff,
//! `ungapped_hits.push` and every diagonal-hash insert with its arguments. The owner sees
//! the reference sequence of pairs and makes every state change with the same statements in
//! the same order, so the hit list and the diagonal hash are the reference's whatever the
//! number of threads and strips.
//!
//! Scope (anything else returns `false`, and the reference code runs): a discontiguous
//! megablast table (one or two templates), `window_size > 0` (two hits), the diagonal hash
//! (query block longer than 8000), `scan_range == 0`, an unmasked subject, no small-query
//! word extension, no debug, trace or window diagnostics, one subject searched by the owner
//! of a pool of two threads or more. A contiguous table in two-hit mode stays with the
//! reference (the parent module takes one-hit contiguous tables only).

use super::super::{
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, mask_lookup_index,
    packed_kmer_at_seed_mask, type_of_word, BlastnTiming, DiagHashTable, UngappedHit,
};
use super::{
    copy_hash, part, plan_strips, pool_threads, run_strips, seed_of_pair, FoldCounts, Inputs,
    PartScratch, Seed, StageStats, Stash, StripJob, STASH,
};
use crate::algorithm::blastn::disc_lookup::{choose_disc_scan_subject, disc_word_scan_subject};
use super::super::x_r4_probe::{self, prefetch, CellPv, CtxTab, Probe, RING};
use crate::algorithm::blastn::lookup::DiscTemplates;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads the switches once (LOSAT_X_PAIRPARDISC, a strict-set switch since
/// 2026-10-17: on unless it is `0`, or `LOSAT_X_DEFAULT=0` and it is not set;
/// LOSAT_X_PAIRPARDISCSHADOW, off unless set to any value); they do not change any value NCBI
/// computes.
fn switches() -> (bool, bool) {
    static SWITCHES: OnceLock<(bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        (
            crate::utils::xswitch::on("LOSAT_X_PAIRPARDISC"),
            std::env::var_os("LOSAT_X_PAIRPARDISCSHADOW").is_some(),
        )
    })
}

/// No NCBI counterpart (scheduling only): strips per thread of a chunk's scan, and strips
/// per thread that the helpers may scan past the one the owner folds. The owner's fold is
/// the slower side here (dc-megablast reaches the diagonal hash with most pairs), so the
/// helpers would otherwise hold the pairs of the whole chunk; small strips and a bounded
/// look-ahead keep that to a few strips (q5000 x 3 Mb at 4 threads: +7 MB instead of
/// +115 MB of seeds).
const STRIPS_PER_THREAD: usize = 64;
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

/// No NCBI counterpart: whether a chunk with a discontiguous table is in the scope of this
/// module (see the module comment); the parent module's scope is the one-hit contiguous
/// table.
fn in_scope(x: &Inputs<'_>) -> bool {
    x.two_stage.disc().is_some()
        && x.window_size > 0
        && !x.use_array_indexing
        && x.scan_range == 0
        && !x.subject_masked
        && !x.small_na_word
        && !x.diagnostics
        && x.one_subject_owner
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1456-1469
/// ```c
/// /** access the megablast lookup table with two templates */
/// #define MB_ACCESS_HITS2()                               \
///    if (total_hits >= max_hits)                          \
///        break;                                           \
///    if (s_BlastMBLookupHasHits(mb_lt, index)) {          \
///        total_hits += s_BlastMBLookupRetrieve(mb_lt,     \
///                   index, offset_pairs + total_hits,     \
///                   scan_range[0]);                               \
///    }                                                    \
///    if (s_BlastMBLookupHasHits(mb_lt, index2)) {         \
///        total_hits += s_BlastMBLookupRetrieve2(mb_lt,    \
///                   index2, offset_pairs + total_hits,    \
///                   scan_range[0]);                               \
///    }
/// ```
/// Phase 1 for one strip: the discontiguous scan kernel `scan_subject_disc_words_with_ranges`
/// picks (`choose_disc_scan_subject`), from `job.start` to `job.end` (stride 1); every word's
/// first-template chain, then its second-template chain, as `on_word` walks them; every pair
/// through `seed_of_pair`.
fn scan_strip_disc(x: &Inputs<'_>, disc: &DiscTemplates, job: &StripJob, out: &mut Vec<Seed>) {
    let two_stage = x.two_stage;
    let s_range = job.s_range;
    let kind = choose_disc_scan_subject(disc.template_type, disc.two_templates);
    disc_word_scan_subject(
        kind,
        x.search_seq_packed,
        disc.template_type,
        disc.second_template_type,
        disc.template_length as i32,
        [job.start as i32, job.end as i32],
        &mut |s_off, index, index2| {
            let kmer_start = s_off as usize;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
            // ```c
            //     while (q_off) {
            //         offset_pairs[i].qs_offsets.q_off   = q_off - 1;
            //         offset_pairs[i++].qs_offsets.s_off = s_off;
            //         q_off = lookup->next_pos[q_off];
            //     }
            // ```
            two_stage.for_each_hit(index as u64, |q_off_1| {
                if q_off_1 != 0 {
                    seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
                }
            });
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1428-1444
            // ```c
            //     Int4 q_off = lookup->hashtable2[index];
            //
            //     while (q_off) {
            //         offset_pairs[i].qs_offsets.q_off   = q_off - 1;
            //         offset_pairs[i++].qs_offsets.s_off = s_off;
            //         q_off = lookup->next_pos2[q_off];
            //     }
            // ```
            if let Some(index2) = index2 {
                two_stage.for_each_hit2(u64::from(index2), |q_off_1| {
                    seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
                });
            }
        },
    );
}

/// EXPERIMENT (exp/r4-scan-probes, v4-R R4; evidence only): `scan_strip_disc` with the R4_PROBE
/// arms (`x_r4_probe.rs`) for one template: the same words in the same order, `pf` words late
/// through a FIFO ring with the backbone prefetch (and with `ST2` the first chain element's), the
/// empty cells dropped with `PV`, and with `CTX` the seeds' contexts from the table.
fn scan_strip_disc_r4<const ST2: bool, const PV: bool, const CTX: bool>(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    job: &StripJob,
    out: &mut Vec<Seed>,
    pf: usize,
    pv: Option<&CellPv>,
    ctx: Option<&CtxTab>,
) {
    let two_stage = x.two_stage;
    let s_range = job.s_range;
    let kind = choose_disc_scan_subject(disc.template_type, disc.two_templates);
    let (hashtable, next_pos) = two_stage.x_r4_raw();
    let mut emit = |s_off: i32, index: u32| {
        let kmer_start = s_off as usize;
        two_stage.for_each_hit(index as u64, |q_off_1| {
            if q_off_1 != 0 {
                if CTX {
                    seed_of_pair_ctx(x, ctx.unwrap(), q_off_1 as usize - 1, kmer_start, s_range, out);
                } else {
                    seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
                }
            }
        });
    };
    let mut ring = [(0i32, 0u32); RING];
    let mut n = 0usize;
    let mask = pf.wrapping_sub(1);
    disc_word_scan_subject(
        kind,
        x.search_seq_packed,
        disc.template_type,
        disc.second_template_type,
        disc.template_length as i32,
        [job.start as i32, job.end as i32],
        &mut |s_off, index, _index2| {
            if PV && !pv.unwrap().test(index as usize) {
                return;
            }
            if pf == 0 {
                emit(s_off, index);
                return;
            }
            prefetch(hashtable.as_ptr().wrapping_add(index as usize));
            if ST2 && n >= pf / 2 {
                let mid = ring[(n - pf / 2) & mask].1 as usize;
                let v = hashtable[mid] as usize;
                if v != 0 {
                    prefetch(next_pos.as_ptr().wrapping_add(v));
                }
            }
            let slot = n & mask;
            if n >= pf {
                let (so, ix) = ring[slot];
                emit(so, ix);
            }
            ring[slot] = (s_off, index);
            n += 1;
        },
    );
    if pf > 0 {
        for m in n.saturating_sub(pf)..n {
            let (so, ix) = ring[m & mask];
            emit(so, ix);
        }
    }
}

/// EXPERIMENT (exp/r4-scan-probes): `seed_of_pair` for a discontiguous table
/// (`word_length == lut_word_length`, no word extension) with the context from `CtxTab`.
#[inline(always)]
fn seed_of_pair_ctx(
    x: &Inputs<'_>,
    ctx: &CtxTab,
    q_off0: usize,
    kmer_start: usize,
    s_range: usize,
    out: &mut Vec<Seed>,
) {
    let lut_word_length = x.two_stage.lut_word_length();
    let context_idx = ctx.context(x.query_context_index, q_off0);
    let c = &x.query_contexts[context_idx];
    let q_context_end = c.query_offset as usize + c.len;
    if q_off0 + lut_word_length > q_context_end || kmer_start + lut_word_length > x.s_len {
        return;
    }
    out.push(Seed {
        context: context_idx as u32,
        q: q_off0 as u32,
        s: kmer_start as u32,
        s_range: s_range as u32,
    });
}

/// EXPERIMENT (exp/r4-scan-probes): the instance of `scan_strip_disc_r4` for the arms.
fn scan_strip_disc_r4_dispatch(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    job: &StripJob,
    out: &mut Vec<Seed>,
    r4: Probe,
    pv: Option<&CellPv>,
    ctx: Option<&CtxTab>,
) {
    let pf = r4.pf;
    match (r4.st2, r4.pv, r4.ctx) {
        (false, false, false) => scan_strip_disc_r4::<false, false, false>(x, disc, job, out, pf, pv, ctx),
        (true, false, false) => scan_strip_disc_r4::<true, false, false>(x, disc, job, out, pf, pv, ctx),
        (false, true, false) => scan_strip_disc_r4::<false, true, false>(x, disc, job, out, pf, pv, ctx),
        (true, true, false) => scan_strip_disc_r4::<true, true, false>(x, disc, job, out, pf, pv, ctx),
        (false, false, true) => scan_strip_disc_r4::<false, false, true>(x, disc, job, out, pf, pv, ctx),
        (true, false, true) => scan_strip_disc_r4::<true, false, true>(x, disc, job, out, pf, pv, ctx),
        (false, true, true) => scan_strip_disc_r4::<false, true, true>(x, disc, job, out, pf, pv, ctx),
        (true, true, true) => scan_strip_disc_r4::<true, true, true>(x, disc, job, out, pf, pv, ctx),
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-841
/// ```c
///     Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
///     Int4 rc;
///
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
///
///     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:882-895,923-941
/// ```c
///             if (!off_found) {
///                 /* This is a new hit */
///                 hit_ready = 0;
///             }
///         }
///     }  else if (check_masks) {
/// 	 /* check the masks for the word */
///         if (!s_TypeOfWord(query, subject, &q_off, &s_off,
///                           query_mask, query_info, s_range,
///                           word_length, lut_word_length, lut, FALSE, &extended)) return 0;
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
/// order, statement for statement for the diagonal hash in two-hit mode; the debug, trace
/// and diagonal-array branches are outside the scope and left out, and the small-query
/// word extension (`small_na_word`) is out of the scope, so `s_TypeOfWord` checks the query
/// masks only.
///
/// With `PT` (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`) the seeds are one partition's seeds
/// of one strip, `diag_hash` and `ungapped_hits` are the partition's, and `out` gets the key
/// (`out.strip`, `pos[k]`) of every pushed hit and of every insert that appends a cell, the
/// buckets the inserts change, and the time of the ungapped extensions (instead of the shared
/// LOSAT_TIMING counters); without it `out` and `pos` are not read.
#[allow(clippy::too_many_arguments)]
fn fold_two_hit<const PT: bool>(
    x: &Inputs<'_>,
    seeds: &[Seed],
    diag_hash: &mut DiagHashTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    counts: &mut FoldCounts,
    out: &mut part::PartOut,
    pos: &[u32],
) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let word_length = two_stage.word_length();
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let diag_offset = x.diag_offset;
    let diag_hash_window = x.diag_hash_window;
    let window_size = x.window_size;
    let scan_range = x.scan_range;
    let two_hits = window_size > 0;
    let timing_enabled = timing.is_some();
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the table has level bounds; the bound's
    // constants of these seeds.
    let x_lv = diag_hash.x_lv_on();
    let x_lv_seeds = super::super::x_diag_lv::LvSeeds::new(two_hits, word_length, window_size);
    counts.seeds += seeds.len() as u64;
    for (x_k, seed) in seeds.iter().enumerate() {
        // EXPERIMENT (LOSAT_X_PAIRPARFOLD): the seed's place in the stage.
        let x_key = if PT { out.key(pos, x_k) } else { 0 };
        let context_idx = seed.context as usize;
        let ctx = &x.query_contexts[context_idx];
        let q_idx = context_idx as u32;
        let q_context_start = ctx.query_offset as usize;
        let q_context_end = q_context_start + ctx.len;
        let cutoff_score = x.cutoff_scores[context_idx];
        let q_ext_start = seed.q as usize;
        let s_ext_start = seed.s as usize;
        let s_range = seed.s_range as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-839
        // ```c
        //     diag = s_off - q_off;
        //     s_end = s_off + word_length;
        //     s_off_pos = s_off + hash_table->offset;
        //     s_end_pos = s_end + hash_table->offset;
        //
        //     rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
        //
        //     /* if there is no record in hashtable, we set last_hit to be a very negative number */
        //     if(!rc)  last_hit = 0;
        //
        //     /* hit within the explored area should be rejected*/
        //     if (s_off_pos < last_hit) return 0;
        // ```
        let diag = s_ext_start as isize - q_ext_start as isize;
        let (last_hit, hit_saved) = {
            // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same decisions without the walk
            // of a bucket below the bound (the `last_hit` test and the window test).
            let (level, _hit_len, hit_saved) = if x_lv {
                diag_hash.retrieve_below(
                    diag as i32,
                    &x_lv_seeds,
                    s_ext_start + diag_offset as usize,
                )
            } else {
                diag_hash.retrieve(diag as i32)
            }
            .unwrap_or((0, 0, false));
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
        let mut off_found = false;
        let mut hit_ready = true;
        let query_mask = ctx.masks.as_slice();

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:41-69 (s_MBLookup)
        // ```c
        //     if (! PV_TEST(pv, index, mb_lt->pv_array_bts)) {
        //         return FALSE;
        //     }
        //
        //     q_off = mb_lt->hashtable[index];
        //     while (q_off) {
        //         if (q_off == q_pos) return TRUE;
        //         q_off = mb_lt->next_pos[q_off];
        //     }
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:841-851
        // ```c
        //     if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
        //         /* check the masks for the word; also check to see if this
        //            first hit qualifies for a double-hit */
        //         word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
        //                                  query_mask, query_info, s_range,
        //                                  word_length, lut_word_length, lut, TRUE, &extended);
        //         if (!word_type) return 0;
        //         /* update the right end*/
        //         s_end += extended;
        //         s_end_pos += extended;
        // ```
        let s_end_pos_i32 = s_end_pos as i32;
        let window_end = last_hit + window_size as i32;
        if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
            let (word_type, extended, q_off_adj, s_off_adj) = type_of_word(
                q_off,
                s_off,
                !query_mask.is_empty(),
                q_context_end,
                s_range,
                word_length,
                lut_word_length,
                true, // check_double = TRUE
                &mut is_seed_masked,
            );
            q_off = q_off_adj;
            s_off = s_off_adj;

            // NCBI: if (!word_type) return 0;
            if word_type == 0 {
                continue;
            }

            // NCBI: s_end += extended; s_end_pos += extended;
            s_end += extended;
            s_end_pos += extended;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:853-886
            // ```c
            //         if (word_type == 1) {
            //             /* try off-diagonal */
            //             Int4 s_a = s_off_pos + word_length - window_size;
            //             Int4 s_b = s_end_pos - 2 * word_length;
            //             Int4 delta;
            //             if (Delta < 0) Delta = 0;
            //             for (delta = 1; delta <= Delta; ++delta) {
            //                 Int4 off_s_end = 0;
            //                 Int4 off_s_l = 0;
            //                 Int4 off_hit_saved = 0;
            //                 Int4 off_rc = s_BlastDiagHashRetrieve(hash_table, diag + delta,
            //                               &off_s_end, &off_s_l, &off_hit_saved);
            //                 if ( off_rc
            //                   && off_s_l
            //                   && off_s_end - delta >= s_a
            //                   && off_s_end - off_s_l <= s_b) {
            //                      off_found = TRUE;
            //                      break;
            //                 }
            //                 off_rc = s_BlastDiagHashRetrieve(hash_table, diag - delta,
            //                               &off_s_end, &off_s_l, &off_hit_saved);
            //                 if ( off_rc
            //                   && off_s_l
            //                   && off_s_end >= s_a
            //                   && off_s_end - off_s_l + delta <= s_b) {
            //                      off_found = TRUE;
            //                      break;
            //                 }
            //             }
            //             if (!off_found) {
            //                 /* This is a new hit */
            //                 hit_ready = 0;
            //             }
            //         }
            // ```
            // (`scan_range` is 0 in the scope, so `Delta` is 0 and the loop runs no
            // iteration; the statements are kept as the reference has them.)
            if word_type == 1 {
                let delta_calc = window_size as isize - word_length as isize;
                let delta_max = if delta_calc < 0 {
                    0
                } else {
                    scan_range.min(delta_calc as usize) as isize
                };
                let s_a = s_off_pos as isize + word_length as isize - window_size as isize;
                let s_b = s_end_pos as isize - 2 * word_length as isize;
                let diag_i32 = diag as i32;
                for delta in 1..=delta_max {
                    if let Some((off_s_end, off_s_l, _)) =
                        diag_hash.retrieve(diag_i32 + delta as i32)
                    {
                        if off_s_l > 0
                            && (off_s_end as isize - delta) >= s_a
                            && (off_s_end as isize - off_s_l as isize) <= s_b
                        {
                            off_found = true;
                            break;
                        }
                    }
                    if let Some((off_s_end, off_s_l, _)) =
                        diag_hash.retrieve(diag_i32 - delta as i32)
                    {
                        if off_s_l > 0
                            && (off_s_end as isize) >= s_a
                            && (off_s_end as isize - off_s_l as isize + delta) <= s_b
                        {
                            off_found = true;
                            break;
                        }
                    }
                }
                if !off_found {
                    hit_ready = false;
                }
            }
        } else {
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
                continue;
            }
            q_off = q_off_adj;
            s_off = s_off_adj;
            s_end += ext;
            s_end_pos += ext;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:897,939-941
        // ```c
        //     if (hit_ready) {
        // ...
        //     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
        //                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
        //                           hit_ready, s_off_pos, window_size + Delta + 1);
        // ```
        if !hit_ready {
            part::insert::<PT>(
                diag_hash,
                out,
                x_key,
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
                // a partition's own counters, added to the shared ones once per stage
                out.ungapped_ns += elapsed_ns;
                out.ungapped_calls += 1;
            } else if let Some(timing) = timing {
                timing.ungapped_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                timing.ungapped_calls.fetch_add(1, Ordering::Relaxed);
            }
            *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
        }
        let ss = ungapped.s_start;
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:928-929
        // ```c
        //                 s_end_pos = ungapped_data->length + ungapped_data->s_start
        //                           + hash_table->offset;
        // ```
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:923-932,939-941
        // ```c
        //             if (off_found || ungapped_data->score >= cutoffs->cutoff_score) {
        // ...
        //             } else {
        //                 hit_ready = 0;
        //             }
        // ...
        //     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
        //                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
        //                           hit_ready, s_off_pos, window_size + Delta + 1);
        // ```
        if !(off_found || ungapped_score >= cutoff_score) {
            hit_ready = false;
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
        // ```c
        //     init_hsp->offsets.qs_offsets.q_off -= query_start;
        //     if (init_hsp->ungapped_data) {
        //         init_hsp->ungapped_data->q_start -= query_start;
        //     }
        //     ASSERT(init_hsp->ungapped_data == NULL ||
        //            init_hsp->ungapped_data->q_start >= 0);
        // ```
        // NCBI checks the start against the context only for saved hits. A rejected
        // approximate result may start left of its context: the left loop of
        // s_NuclUngappedExtend runs over the concatenated query and reads the sentinel as
        // bases (na_ungapped.c:284-305, `len = MIN(q_ext, s_ext) / COMPRESSION_RATIO;`),
        // and that result is kept only below reduced_nucl_cutoff_score = 0.8 * cutoff_score
        // (blast_parameters.c:412), so it is rejected above and never reaches this point.
        debug_assert!(ungapped.q_start >= q_context_start);
        let qs = ungapped.q_start - q_context_start;
        let qe = qs + ungapped.length;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:924-927
        // ```c
        //                 BlastUngappedData *final_data =
        //                     (BlastUngappedData *) malloc(sizeof(BlastUngappedData));
        //                 *final_data = *ungapped_data;
        //                 BLAST_SaveInitialHit(init_hitlist, q_off, s_off, final_data);
        // ```
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
        if PT {
            out.hit_keys.push(x_key);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
        // ```c
        //     s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
        //                           (hit_ready) ? 0 : s_end_pos - s_off_pos,
        //                           hit_ready, s_off_pos, window_size + Delta + 1);
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

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1327-1331
/// ```c
///    if (lookup_options->mb_template_length > 0) {
///         /* discontiguous megablast */
///         mb_lt->scan_step = 1;
///         status = s_FillDiscMBTable(query, location, mb_lt, lookup_options);
///    }
/// ```
/// The seed stage of one chunk on `threads` threads (this one and `threads - 1` helpers of
/// its pool): the scan ranges of `scan_subject_disc_words_with_ranges` (an unmasked
/// subject: `scan_range[2] = subject->length - template_length`, stride 1) cut into strips
/// scanned by any thread, folded here in strip order into `hash` and `hits`. With `fold_part`
/// (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`) the fold is split by diagonal-hash bucket over
/// the threads instead.
#[allow(clippy::too_many_arguments)]
fn seed_stage(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
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
        disc.template_length,
        disc.template_length,
        1,
        STRIPS_PER_THREAD * threads.max(4),
    );
    // EXPERIMENT (exp/r4-scan-probes, v4-R R4): the R4_PROBE arms (one template only).
    let r4 = x_r4_probe::probe();
    let r4_on = r4.on() && !disc.two_templates;
    let r4_pv = (r4_on && r4.pv).then(|| CellPv::build(x.two_stage.x_r4_raw().0));
    let r4_ctx = (r4_on && r4.ctx).then(|| {
        assert_eq!(x.two_stage.word_length(), x.two_stage.lut_word_length());
        CtxTab::build(x.query_context_index)
    });
    let produce = |i: usize, out: &mut Vec<Seed>| {
        if r4_on {
            scan_strip_disc_r4_dispatch(x, disc, &jobs[i], out, r4, r4_pv.as_ref(), r4_ctx.as_ref());
        } else {
            scan_strip_disc(x, disc, &jobs[i], out);
        }
    };
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
            fold_two_hit::<true>(
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
        part::record(part::TWO_HITS, parts, &r);
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
            fold_two_hit::<false>(
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
/// The parent's `run` for a discontiguous table: the seed stage of this chunk (the loop
/// above with its two-hit extension) when LOSAT_X_PAIRPARDISC is on and the chunk is in
/// scope; `true` means it has run and the reference scan is skipped. With
/// LOSAT_X_PAIRPARDISCSHADOW the stage runs on copies of the diagonal hash and the hit list,
/// `false` is returned so the reference runs, and the parent's `shadow_check` compares the
/// two. `fold_scratch` holds the allocations of LOSAT_X_PAIRPARFOLD's partition tables and
/// buffers (`x_pair_seed_part.rs`).
pub(super) fn run(
    x: &Inputs<'_>,
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
    let Some(disc) = x.two_stage.disc() else {
        return false;
    };
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
        seed_stage(x, disc, threads, hash, hits, timing, timed, fold_part)
    };
    if let (Some(t_start), Some(timing)) = (t_start, timing) {
        // `[TIMING] scan_lookup`: the time of the stage minus the ungapped extensions, as
        // the reference counts its scan (which also leaves out the last flush). With
        // LOSAT_X_PAIRPARFOLD only the owner's extensions are taken off (the partitions'
        // extension times summed over the threads exceed the stage's wall time).
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
                disc: true,
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

/// No NCBI counterpart: the one summary line of LOSAT_X_PAIRPARDISCSHADOW, and with
/// LOSAT_TIMING the counters of LOSAT_X_PAIRPARDISC, on stderr at the end of the search.
pub(super) fn summary() {
    let (on, shadow) = switches();
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if shadow {
        eprintln!(
            "[X_PAIRPARDISCSHADOW] blastn: seed chunks equal={} (ungapped hits {}, hash cells {}){}",
            load(&SHADOW_CHUNKS),
            load(&SHADOW_HITS),
            load(&SHADOW_CELLS),
            part::shadow_fields(part::TWO_HITS)
        );
    }
    if (on || shadow) && std::env::var_os("LOSAT_TIMING").is_some() {
        let s = |c: &AtomicU64| c.load(Ordering::Relaxed) as f64 / 1e9;
        eprintln!(
            "[TIMING] x_pairpardisc_seed: chunks={} strips={} seeds={} skipped={} fold={:.3}s owner_scan={:.3}s owner_wait={:.3}s{}",
            load(&X_CHUNKS),
            load(&X_STRIPS),
            load(&X_SEEDS),
            load(&X_SKIPPED),
            s(&X_FOLD_NS),
            s(&X_OWNER_SCAN_NS),
            s(&X_WAIT_NS),
            part::timing_fields(part::TWO_HITS)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::{
        build_blastna_matrix, build_nucl_score_table, build_query_blastna_buffer,
        diag_hash_insert_window, scan_subject_disc_words_with_ranges, QueryContext,
        QueryContextIndex,
    };
    use super::*;
    use crate::algorithm::blastn::disc_lookup::DiscWordType;
    use crate::algorithm::blastn::lookup::{
        build_disc_mb_lookup, build_two_stage_lookup, reverse_complement, TwoStageLookup,
    };
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

    fn random_dna(rng: &mut Lcg, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
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

    /// A query of two contexts (plus and minus strand, masks on the plus strand) and a
    /// subject made of mutated copies of query pieces on both strands and random pieces.
    struct Case {
        contexts: Vec<QueryContext>,
        index: QueryContextIndex,
        ctx_seqs: Vec<Vec<u8>>,
        offsets: Vec<i32>,
        masks: Vec<Vec<MaskedInterval>>,
        concat: Vec<u8>,
        concat_sentinels: Vec<u8>,
        packed: Vec<u8>,
        s_len: usize,
        seq_ranges: [(i32, i32); 1],
        score_matrix: [i32; 256],
        nucl_score_table: [i32; 256],
    }

    fn make_case(rng: &mut Lcg, subject_len: usize) -> Case {
        let q_len = 4100 + rng.below(3000) as usize;
        let mut query = random_dna(rng, q_len);
        // a repeat inside the query (long lookup chains)
        let motif = random_dna(rng, 40);
        for k in 0..20 {
            let at = (k * 97) % (q_len - 40);
            query[at..at + 40].copy_from_slice(&motif);
        }
        let minus = reverse_complement(&query);
        let mut subject = Vec::new();
        while subject.len() < subject_len {
            if rng.below(3) == 0 {
                let len = 50 + rng.below(400) as usize;
                subject.extend(random_dna(rng, len));
            } else {
                let src = if rng.below(2) == 0 { &query } else { &minus };
                let len = 30 + rng.below(1500) as usize;
                let at = rng.below((q_len - len.min(q_len - 1)) as u64) as usize;
                let end = (at + len).min(q_len);
                for &b in &src[at..end] {
                    // dc-megablast finds diverged copies: 10 % substitutions
                    let r = rng.below(100);
                    if r < 10 {
                        subject.push(b"ACGT"[rng.below(4) as usize]);
                    } else if r < 11 {
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
        let concat_sentinels = build_query_blastna_buffer(&[&query], &offsets, concat_len);
        let concat = concat_sentinels[1..concat_len + 1].to_vec();
        let contexts: Vec<QueryContext> = (0..2)
            .map(|c| QueryContext {
                query_idx: 0,
                frame: if c == 0 { 1 } else { -1 },
                query_offset: offsets[c],
                len: ctx_seqs[c].len(),
                masks: masks[c].clone(),
            })
            .collect();
        let index = QueryContextIndex::new(&contexts);
        Case {
            contexts,
            index,
            ctx_seqs,
            offsets,
            masks,
            concat,
            concat_sentinels,
            packed,
            s_len,
            seq_ranges: [(0, s_len as i32)],
            score_matrix: build_blastna_matrix(2, -3),
            nucl_score_table: build_nucl_score_table(2, -3),
        }
    }

    fn ctx_refs(c: &Case) -> Vec<&[u8]> {
        c.ctx_seqs.iter().map(Vec::as_slice).collect()
    }

    fn inputs<'a>(
        c: &'a Case,
        lookup: &'a TwoStageLookup,
        diag_offset: isize,
        window_size: usize,
    ) -> Inputs<'a> {
        Inputs {
            two_stage: lookup,
            search_seq_packed: &c.packed,
            s_len: c.s_len,
            subject_seq_ranges: &c.seq_ranges,
            subject_masked: false,
            scan_step: 1,
            query_context_index: &c.index,
            query_contexts: &c.contexts,
            encoded_query_concat_blastna: &c.concat,
            encoded_query_concat_blastna_with_sentinels: &c.concat_sentinels,
            cutoff_scores: &[30, 30],
            x_dropoff_scores: &[20, 20],
            reduced_cutoff_scores: &[24, 24],
            score_matrix: &c.score_matrix,
            nucl_score_table: &c.nucl_score_table,
            diag_offset,
            diag_hash_window: diag_hash_insert_window(window_size, 0, lookup.word_length()),
            use_array_indexing: false,
            window_size,
            scan_range: 0,
            small_na_word: false,
            diagnostics: false,
            one_subject_owner: true,
            timing: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1651-1655
    // ```c
    //     if (subject->mask_type != eNoSubjMasking) {
    //         if (lookup_wrap->lut_type == eMBLookupTable &&
    //             ((BlastMBLookupTable *) lookup_wrap->lut)->discontiguous) {
    //             /* discontiguous scan subs assumes any (non-aligned starting offset */
    // ```
    // The whole two-hit stage on random sequences for every discontiguous template (word
    // size 11 and 12, coding, optimal and both, lengths 16, 18 and 21; each one-template
    // type picks one of the scan kernels, two templates the two-template kernel): the pairs
    // of the reference wrapper `scan_subject_disc_words_with_ranges` (one scan of the whole
    // range) folded in one go, against `seed_stage` on 1..8 threads (strips planned for each
    // thread count, so the kernels restart at many offsets). Hit lists and diagonal hashes
    // must be the same, and the serial fold must have taken every branch: the `last_hit`
    // rejection, the first hit of a window (insert with `hit_ready` 0), the second hit
    // (ungapped extension) both below and above the cutoff.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn disc_two_hit_stage_on_strips_equals_one_scan() {
        let mut rng = Lcg(0xD15C_5EED);
        let mut case_no = 0;
        for word_size in [11usize, 12] {
            for word_type in [
                DiscWordType::Coding,
                DiscWordType::Optimal,
                DiscWordType::TwoTemplates,
            ] {
                for template_length in [16usize, 18, 21] {
                    case_no += 1;
                    let c = make_case(&mut rng, 30_000);
                    let lookup = build_disc_mb_lookup(
                        &ctx_refs(&c),
                        &c.offsets,
                        word_size,
                        template_length,
                        word_type,
                        &c.masks,
                        c.concat.len(),
                    );
                    let disc = lookup.disc().expect("a discontiguous table");
                    assert_eq!(lookup.word_length(), template_length);
                    let diag_offset = 37 + rng.below(1000) as isize;
                    let window_size = 40;
                    let x = inputs(&c, &lookup, diag_offset, window_size);
                    assert!(in_scope(&x));
                    let new_hash = || {
                        let mut h = DiagHashTable::new(window_size as i32);
                        h.offset = diag_offset as i32;
                        h
                    };
                    let mut seeds = Vec::new();
                    scan_subject_disc_words_with_ranges(
                        &c.packed,
                        c.s_len,
                        disc,
                        &c.seq_ranges,
                        false,
                        |s_off, s_range, index, index2| {
                            lookup.for_each_hit(index, |q_off_1| {
                                if q_off_1 != 0 {
                                    seed_of_pair(
                                        &x,
                                        q_off_1 as usize - 1,
                                        s_off,
                                        s_range,
                                        &mut seeds,
                                    );
                                }
                            });
                            if let Some(index2) = index2 {
                                lookup.for_each_hit2(index2, |q_off_1| {
                                    seed_of_pair(
                                        &x,
                                        q_off_1 as usize - 1,
                                        s_off,
                                        s_range,
                                        &mut seeds,
                                    );
                                });
                            }
                        },
                    );
                    let mut serial_hash = new_hash();
                    let mut serial_hits = Vec::new();
                    let mut counts = FoldCounts::default();
                    let mut pause = 0u64;
                    fold_two_hit::<false>(
                        &x,
                        &seeds,
                        &mut serial_hash,
                        &mut serial_hits,
                        None,
                        &mut pause,
                        &mut counts,
                        &mut part::PartOut::default(),
                        &[],
                    );
                    let saved = serial_hash.chain[1..serial_hash.occupancy as usize]
                        .iter()
                        .filter(|cell| cell.hit_saved)
                        .count();
                    let first_hits = serial_hash.chain[1..serial_hash.occupancy as usize]
                        .iter()
                        .filter(|cell| !cell.hit_saved && cell.hit_len > 0)
                        .count();
                    assert!(
                        serial_hits.len() > 5 && counts.skipped > 0 && saved > 0 && first_hits > 0,
                        "case {case_no}: the test sequences give too few hits ({} hits, {} skipped, {saved} saved cells, {first_hits} first-hit cells)",
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
                        let (state, keys, stats_seeds, strips) = pool.install(|| {
                            let mut hash = new_hash();
                            let mut hits = Vec::new();
                            let stats = seed_stage(
                                &x, disc, threads, &mut hash, &mut hits, None, false, None,
                            );
                            (
                                hash_state(&hash),
                                hits.iter().map(key).collect::<Vec<_>>(),
                                stats.counts.seeds,
                                stats.strips,
                            )
                        });
                        assert!(strips > 1, "case {case_no} threads {threads}");
                        assert_eq!(
                            stats_seeds,
                            seeds.len() as u64,
                            "case {case_no} threads {threads}"
                        );
                        assert_eq!(keys, serial_keys, "case {case_no} threads {threads}");
                        assert_eq!(
                            state,
                            hash_state(&serial_hash),
                            "case {case_no} threads {threads}"
                        );
                    }
                    // LOSAT_X_DIAGLV (`x_diag_lv.rs`): on a table with level bounds (and with
                    // the shadow, which walks and checks every skipped bucket) the serial fold
                    // and the strips on 1, 2, 3 and 7 threads skip walks and give the hits and
                    // the table of the fold without bounds, and the same bounds.
                    for shadow in [false, true] {
                        let new_lv_hash = || {
                            let mut h = new_hash();
                            h.x_lv_start(shadow);
                            h
                        };
                        let mut lv_hash = new_lv_hash();
                        let mut lv_hits = Vec::new();
                        let before = super::super::super::x_diag_lv::test_skips();
                        fold_two_hit::<false>(
                            &x,
                            &seeds,
                            &mut lv_hash,
                            &mut lv_hits,
                            None,
                            &mut pause,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                        let skipped = super::super::super::x_diag_lv::test_skips() - before;
                        assert!(skipped > 0, "case {case_no}: no walk skipped");
                        assert_eq!(lv_hits.iter().map(key).collect::<Vec<_>>(), serial_keys);
                        assert_eq!(
                            hash_state(&lv_hash),
                            hash_state(&serial_hash),
                            "case {case_no}"
                        );
                        for threads in [1usize, 2, 3, 7] {
                            let pool = rayon::ThreadPoolBuilder::new()
                                .num_threads(threads)
                                .build()
                                .unwrap();
                            let (state, keys, max_level, skipped) = pool.install(|| {
                                let mut hash = new_lv_hash();
                                let mut hits = Vec::new();
                                let before = super::super::super::x_diag_lv::test_skips();
                                seed_stage(
                                    &x, disc, threads, &mut hash, &mut hits, None, false, None,
                                );
                                (
                                    hash_state(&hash),
                                    hits.iter().map(key).collect::<Vec<_>>(),
                                    hash.max_level.clone(),
                                    super::super::super::x_diag_lv::test_skips() - before,
                                )
                            });
                            let at = format!("case {case_no} threads {threads} shadow {shadow}");
                            assert!(skipped > 0, "{at}: no walk skipped");
                            assert_eq!(keys, serial_keys, "{at}");
                            assert_eq!(state, hash_state(&serial_hash), "{at}");
                            assert_eq!(max_level, lv_hash.max_level, "{at}");
                        }
                    }
                    // LOSAT_X_PAIRPARFOLD (`x_pair_seed_part.rs`): the two-hit fold split into
                    // 2, 3, 7, 8 and 16 partitions on 1..8 threads, on tables with and without
                    // level bounds, for two chunks in a row on one scratch (the second with the
                    // hash offset moved on): the hits and the table of the serial fold.
                    let x2 = Inputs {
                        diag_offset: diag_offset + c.s_len as isize + window_size as isize,
                        ..x
                    };
                    for levels in [false, true] {
                        let new_table = || {
                            let mut h = new_hash();
                            if levels {
                                h.x_lv_start(false);
                            }
                            h
                        };
                        let mut expected_hash = new_table();
                        let mut expected = Vec::new();
                        for xc in [&x, &x2] {
                            expected_hash.offset = xc.diag_offset as i32;
                            let mut hits = Vec::new();
                            fold_two_hit::<false>(
                                xc,
                                &seeds,
                                &mut expected_hash,
                                &mut hits,
                                None,
                                &mut pause,
                                &mut FoldCounts::default(),
                                &mut part::PartOut::default(),
                                &[],
                            );
                            expected.push((
                                hits.iter().map(key).collect::<Vec<_>>(),
                                hash_state(&expected_hash),
                                expected_hash.max_level.clone(),
                            ));
                        }
                        for (threads, parts) in [
                            (1usize, 2usize),
                            (2, 2),
                            (3, 3),
                            (3, 7),
                            (4, 4),
                            (5, 16),
                            (8, 8),
                            (8, 16),
                        ] {
                            let pool = rayon::ThreadPoolBuilder::new()
                                .num_threads(threads)
                                .build()
                                .unwrap();
                            pool.install(|| {
                                let mut hash = new_table();
                                let mut scratch = PartScratch::default();
                                for (chunk, xc) in [&x, &x2].into_iter().enumerate() {
                                    let at = format!(
                                        "case {case_no} threads {threads} parts {parts} levels {levels} chunk {chunk}"
                                    );
                                    hash.offset = xc.diag_offset as i32;
                                    let mut hits = Vec::new();
                                    let stats = seed_stage(
                                        xc,
                                        disc,
                                        threads,
                                        &mut hash,
                                        &mut hits,
                                        None,
                                        false,
                                        Some(part::FoldPart {
                                            scratch: &mut scratch,
                                            parts,
                                        }),
                                    );
                                    assert_eq!(stats.counts.seeds, seeds.len() as u64, "{at}");
                                    let (keys, state, max_level) = &expected[chunk];
                                    assert_eq!(
                                        &hits.iter().map(key).collect::<Vec<_>>(),
                                        keys,
                                        "{at}"
                                    );
                                    assert_eq!(&hash_state(&hash), state, "{at}");
                                    assert_eq!(&hash.max_level, max_level, "{at}");
                                }
                            });
                        }
                    }
                }
            }
        }
    }

    // The scope: a discontiguous table in two-hit mode with the diagonal hash only. A
    // contiguous table in two-hit mode (megablast or blastn with `-window_size`, which the
    // CLI rejects) is in neither this scope nor the parent's: the parent's `run` refuses it
    // (window_size != 0) and the reference runs.
    #[test]
    fn scope_is_the_two_hit_discontiguous_stage() {
        let mut rng = Lcg(0x5C0_9E);
        let c = make_case(&mut rng, 5_000);
        let disc = build_disc_mb_lookup(
            &ctx_refs(&c),
            &c.offsets,
            11,
            18,
            DiscWordType::Coding,
            &c.masks,
            c.concat.len(),
        );
        let x = inputs(&c, &disc, 40, 40);
        assert!(in_scope(&x));
        assert!(!in_scope(&Inputs {
            window_size: 0,
            ..inputs(&c, &disc, 40, 40)
        }));
        assert!(!in_scope(&Inputs {
            use_array_indexing: true,
            ..inputs(&c, &disc, 40, 40)
        }));
        assert!(!in_scope(&Inputs {
            scan_range: 1,
            ..inputs(&c, &disc, 40, 40)
        }));
        assert!(!in_scope(&Inputs {
            subject_masked: true,
            ..inputs(&c, &disc, 40, 40)
        }));
        assert!(!in_scope(&Inputs {
            diagnostics: true,
            ..inputs(&c, &disc, 40, 40)
        }));
        assert!(!in_scope(&Inputs {
            one_subject_owner: false,
            ..inputs(&c, &disc, 40, 40)
        }));
        let contiguous = build_two_stage_lookup(
            &ctx_refs(&c),
            &c.offsets,
            28,
            12,
            &c.masks,
            None,
            0,
            c.concat.len(),
            false,
        );
        assert!(!in_scope(&inputs(&c, &contiguous, 40, 40)));
        // The parent's `run` in a pool of two threads (where its one-hit stage would run):
        // a contiguous table with a window is refused and nothing is touched.
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(2)
                .build()
                .unwrap();
            pool.install(|| {
                let mut hash = DiagHashTable::new(40);
                let mut hits = Vec::new();
                let x = inputs(&c, &contiguous, 40, 40);
                let mut fold_scratch = PartScratch::default();
                assert!(!super::super::run(
                    &x,
                    &mut hash,
                    &mut hits,
                    &mut fold_scratch
                ));
                assert!(hits.is_empty());
                assert_eq!(hash.occupancy, 1);
            });
        }
    }
}
