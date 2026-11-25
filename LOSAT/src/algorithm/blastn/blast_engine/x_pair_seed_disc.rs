//! EXPERIMENT (LOSAT_X_PAIRPARDISC): the two-hit seed
//! stage of discontiguous megablast (`-task dc-megablast`) for one subject chunk with the
//! diagonal hash, with the subject scan cut into strips that the threads of the search pool
//! scan, and the diagonal-hash fold done by the owner thread in the original pair order.
//! The scheduling (`plan_strips`, `run_strips`) and the per-pair word step (`seed_of_pair`)
//! are those of PAIRPAR in the parent module.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1624-1634,1673-1684
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1651-1655
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
//! the two-hit diagonal-hash branch of `extend_offset_pairs` (`pair_extend.rs`): the retrieve with
//! `last_hit` and `hit_saved`, the `last_hit` test, the window test, `s_TypeOfWord` with and
//! without `check_double`, the off-diagonal search (`pair_extend::off_diagonal_found`), the
//! ungapped extension, the cutoff, `ungapped_hits.push` and every diagonal-hash insert with its
//! arguments. The owner sees
//! the reference sequence of pairs and makes every state change with the same statements in
//! the same order, so the hit list and the diagonal hash are the reference's whatever the
//! number of threads and strips.
//!
//! Scope (anything else returns `false`, and the reference code runs): a discontiguous
//! megablast table (one or two templates), `window_size > 0` (two hits), the diagonal hash
//! (query block longer than 8000), `scan_range == 0`, no small-query word extension, no debug,
//! trace or window diagnostics, a subject searched by the owner of a pool of two threads or
//! more by itself (one subject, or several one at a time under LOSAT_X_BNFEWSUBJ), or with
//! LOSAT_X_BNDCONE by a search on one thread (no pool: the strips are scanned and folded in
//! order on that thread). A masked subject is scanned in its `seq_ranges` as the reference
//! scans it, the words of each range with that range's `s_range`. A contiguous table in
//! two-hit mode stays with the reference (the parent module takes one-hit contiguous tables
//! only).
//!
//! EXPERIMENT (LOSAT_X_PAIRPARDISCARR): the same scope
//! with the diagonal array instead of the hash (a query block of at most 8000 letters) runs
//! on the same strips, folded by the owner thread through `fold_two_hit_array`, a
//! statement-for-statement copy of the two-hit diagonal-array branch of
//! `extend_offset_pairs` (`run_array`, `seed_stage_array`).

use super::super::x_ctx_o1::{XCtxGet, XNone, XView};
use super::super::{
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, mask_lookup_index,
    packed_kmer_at_seed_mask, type_of_word, BlastnTiming, DiagHashTable, DiagStruct, QueryContext,
    QueryContextIndex, UngappedHit,
};
use super::latest::{self, DiagLatest, Found};
use super::{
    part, plan_strips, pool_threads, run_strips, seed_of_pair, Cut, DiagArray, FoldCounts, Inputs,
    PartScratch, Resume, Seed, SeedOut, StageStats, StripJob,
};
use crate::algorithm::blastn::disc_lookup::{
    choose_disc_scan_subject, disc_word_scan_subject_until,
};
use crate::algorithm::blastn::lookup::x_dc_csr::{Csr, Ring};
use crate::algorithm::blastn::lookup::{DiscTemplates, XHeadFifo};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads LOSAT_X_PAIRPARDISC once (a strict-set switch since 2026-10-17: on
/// unless it is `0`), for the parent's `takes_owner_chunks`; it does not change any value NCBI
/// computes.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARDISC"))
}

/// No NCBI counterpart: reads LOSAT_X_PAIRPARDISCARR once (a strict-set switch since 2026-10-20:
/// on unless it is `0`), for the parent's `takes_owner_chunks`; it does not change any value NCBI
/// computes.
pub(super) fn arr_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_PAIRPARDISCARR"))
}

/// No NCBI counterpart: reads LOSAT_X_FOLDLEAN once (a strict-set switch since 2026-10-22: on
/// unless it is `0`). With it the strips folds
/// (`fold_two_hit`, `fold_two_hit_array`) run their `LEAN` instances: no call of
/// `type_of_word` where it returns at once (`word_length == lut_word_length`, always
/// so in the scope: a discontiguous word has no word extension), and the diagonal-hash insert
/// inlined (`DiagHashTable::insert_inl`). The same statements in the same order, so it does not
/// change any value NCBI computes.
fn lean_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_FOLDLEAN"))
}

/// No NCBI counterpart: reads LOSAT_X_DCSEEDLEAN once (a switch of the strict set: on unless
/// `LOSAT_X_DCSEEDLEAN=0`). With it the partitioned stage of the strips (`seed_stage_of` with a
/// `part::FoldPart`) folds without recording the buckets it changes and, in LOSAT_X_PAIRPARHOME's
/// home loop, scans each strip into its partitions' lists (`seed_stage_lists`). The same seeds folded in the same order into the
/// same tables, so it does not change any value NCBI computes.
fn seed_lean_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCSEEDLEAN"))
}

/// No NCBI counterpart: reads LOSAT_X_BNDCONE once (a strict-set switch since 2026-10-24: on unless
/// it is `0`). With it a chunk in the scope of either stage (`in_scope`, `in_scope_array`) that
/// has no pool (a search on one thread) runs the stage on this thread: its strips are scanned and
/// folded in order. The owner folds the pairs in
/// the reference order with the reference statements whatever the number of threads (module
/// comment), so it does not change any value NCBI computes.
fn one_thread_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNDCONE"))
}

/// No NCBI counterpart: reads LOSAT_X_DCLATEST once (a strict-set switch since 2026-11-23: on
/// unless it is `0`). With it the strips stage of the diagonal hash folds into maps of each
/// diagonal's latest insert (`x_dc_latest.rs`: the owner's map without partitions, each
/// partition's with LOSAT_X_PAIRPARFOLD) instead of the diagonal hash, and the partitioned merge
/// writes back no table. Each seed makes the decisions of the reference fold (the proof in the
/// module comment of `x_dc_latest.rs`), so it does not change any value NCBI computes.
fn latest_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCLATEST"))
}

/// No NCBI counterpart: reads LOSAT_X_DCEMITHOIST once (a strict-set switch since 2026-11-26: on
/// unless it is `0`). With it the strips scan of a table laid out as CSR (`scan_strip_disc_csr`,
/// LOSAT_X_DCCSR) hands each word's run to `DcPairs::run` (its `HOIST` instance): the values
/// `seed_of_pair` reads for every pair read once per strip, the context form chosen once per run,
/// the subject's bound tested once per word, and a strip's list given room for the run at once.
/// The same pairs give the same seeds in the same order (the argument at `DcPairs`), so it does not
/// change any value NCBI computes.
fn emit_hoist_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCEMITHOIST"))
}

// No NCBI counterpart: counter of LOSAT_TIMING (the strips scanned through the hoisted CSR emit of
// LOSAT_X_DCEMITHOIST).
static EMIT_HOIST_STRIPS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: LOSAT_X_REACH for LOSAT_X_DCLATEST, which acts in the strips stage of the
/// diagonal hash only: `own_path` where it folds into the maps, false where the chunk goes to the
/// reference fold. Nothing while the switch is off.
fn reach_latest(own_path: bool) {
    if latest_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_DCLATEST", own_path);
    }
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_DCLATEST) the level bound (`x_diag_lv::LvSeeds`) of
/// the first position of the chunk's stage, where every position of the chunk has one (module
/// comment of `x_dc_latest.rs`, (ii) and (vii)). Every chunk has: its positions lie between the
/// hash offset and the offset plus its length, the offset at least the window and below
/// a quarter of the largest Int4 plus a chunk (NCBI's clear), and a chunk at most twice `MAX_DBSEQ_LEN` letters,
/// far below `i32::MAX` less a word. `None` would send the chunk to the chains.
fn latest_bound(x: &Inputs<'_>) -> Option<i32> {
    let seeds =
        super::super::x_diag_lv::LvSeeds::new(true, x.two_stage.word_length(), x.window_size);
    let first = x.diag_offset as usize;
    let bound0 = seeds.bound(first);
    (bound0 != i32::MIN && seeds.bound(first + x.s_len) != i32::MIN).then_some(bound0)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// EXPERIMENT (LOSAT_X_DCLATEST) whether the chunk's scan ranges, in the order the stage scans
/// them (`plan_strips`: the reference's walk, `first_scan_range` and the ranges after it), never
/// step back: each range with positions starts at or past the last position of the one before.
/// The maps need S never decreasing (module comment of `x_dc_latest.rs`, (vi)). A subject of
/// `-db` with `-db_soft_mask` has the complement of its stored masks in their stored order
/// (`subject_split.rs`), which can step back (masks stored 600-699 then 100-199); NCBI scans them
/// in that order, and which cells its chains still hold after the step back depends on the order
/// of their take-overs, so such a chunk goes to the chains. One pass over the chunk's ranges.
fn latest_ranges_ascend(x: &Inputs<'_>, template_length: usize) -> bool {
    if x.subject_seq_ranges.is_empty() {
        return true;
    }
    let first = super::super::first_scan_range(
        x.subject_seq_ranges,
        x.s_len,
        template_length,
        template_length,
        x.subject_masked,
    );
    let mut last: Option<usize> = None;
    let mut ascend = true;
    super::super::scan_subject_kmers_with_offsets(
        x.subject_seq_ranges,
        template_length,
        template_length,
        first,
        |start, end| {
            if last.is_some_and(|last| start < last) {
                ascend = false;
            }
            last = Some(end);
        },
    );
    ascend
}

/// No NCBI counterpart: whether LOSAT_TIMING is set (read once). The stages time themselves when
/// it is; a chunk without a pool, one of many on one thread with LOSAT_X_BNDCONE, does not read the
/// environment again.
fn timing_env() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

/// No NCBI counterpart: LOSAT_X_REACH for LOSAT_X_FOLDLEAN, which acts in the strips folds only:
/// `own_path` where a chunk's strips stage runs, false where the chunk goes to the reference
/// fold. Nothing while the switch is off.
fn reach_lean(own_path: bool) {
    if lean_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_FOLDLEAN", own_path);
    }
}

/// No NCBI counterpart (scheduling only): strips per thread of a chunk's scan, and strips
/// per thread that the helpers may scan past the one the owner folds. The owner's fold is
/// the slower side here (dc-megablast reaches the diagonal hash with most pairs), so the
/// helpers would otherwise hold the pairs of the whole chunk; small strips and a bounded
/// look-ahead keep that to a few strips (q5000 x 3 Mb at 4 threads: +7 MB instead of
/// +115 MB of seeds).
const STRIPS_PER_THREAD: usize = 64;
const AHEAD_PER_THREAD: usize = 2;

// No NCBI counterpart: counters of LOSAT_TIMING (scheduling only).
static X_CHUNKS: AtomicU64 = AtomicU64::new(0);
static X_STRIPS: AtomicU64 = AtomicU64::new(0);
static X_SEEDS: AtomicU64 = AtomicU64::new(0);
static X_SKIPPED: AtomicU64 = AtomicU64::new(0);
static X_FOLD_NS: AtomicU64 = AtomicU64::new(0);
static X_OWNER_SCAN_NS: AtomicU64 = AtomicU64::new(0);
static X_WAIT_NS: AtomicU64 = AtomicU64::new(0);
// (the same for the diagonal-array stage, LOSAT_X_PAIRPARDISCARR)
static ARR_CHUNKS: AtomicU64 = AtomicU64::new(0);
static ARR_STRIPS: AtomicU64 = AtomicU64::new(0);
static ARR_SEEDS: AtomicU64 = AtomicU64::new(0);
static ARR_SKIPPED: AtomicU64 = AtomicU64::new(0);
static ARR_FOLD_NS: AtomicU64 = AtomicU64::new(0);
static ARR_OWNER_SCAN_NS: AtomicU64 = AtomicU64::new(0);
static ARR_WAIT_NS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: whether a chunk with a discontiguous table is in the scope of this
/// module (see the module comment); the parent module's scope is the one-hit contiguous
/// table.
fn in_scope(x: &Inputs<'_>) -> bool {
    x.two_stage.disc().is_some()
        && super::chunk_free_scope(
            true,
            x.use_array_indexing,
            x.window_size,
            x.scan_range,
            x.small_na_word,
            x.diagnostics,
        )
        && x.one_subject_owner
}

/// No NCBI counterpart: whether a chunk with a discontiguous table and the diagonal array (a
/// query block of at most 8000 letters) is in the scope of the stage of LOSAT_X_PAIRPARDISCARR:
/// the scope of `in_scope` with the array instead of the hash.
fn in_scope_array(x: &Inputs<'_>) -> bool {
    x.two_stage.disc().is_some()
        && x.use_array_indexing
        && super::chunk_free_scope(
            true,
            false,
            x.window_size,
            x.scan_range,
            x.small_na_word,
            x.diagnostics,
        )
        && x.one_subject_owner
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1456-1469
/// Phase 1 for one strip: the discontiguous scan kernel `scan_subject_disc_words_with_ranges`
/// picks (`choose_disc_scan_subject`), from `job.start` to `job.end` (stride 1); every word's
/// first-template chain, then its second-template chain, as `on_word` of
/// `two_stage_scan::scan_two_stage_lookup` walks them; every pair
/// through `seed_of_pair`, into `out` (a list, or with LOSAT_X_DCSEEDLEAN the strip's
/// partitions' lists).
fn scan_strip_disc<O: SeedOut>(x: &Inputs<'_>, disc: &DiscTemplates, job: &StripJob, out: &mut O) {
    scan_strip_disc_from::<O, false>(x, disc, job, Resume::at(job.start), out);
}

/// `scan_strip_disc` from the word at `from.s` (`job.start`, or where a piece stopped), that
/// word from where `from` says (`Resume`: the element of the first or the second template's
/// chain). EXPERIMENT (LOSAT_X_SEEDPIECE): with `PIECE` the scan stops at the first pair `out`
/// has no room for (`Cut`) and returns where it goes on; `None` when the strip is done.
pub(super) fn scan_strip_disc_from<O: SeedOut, const PIECE: bool>(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    job: &StripJob,
    from: Resume,
    out: &mut O,
) -> Option<Resume> {
    // EXPERIMENT (LOSAT_X_DCCSR): a table laid out as CSR has no chains to walk.
    if let Some(csr) = disc.x_csr() {
        // EXPERIMENT (LOSAT_X_DCEMITHOIST): its runs through `DcPairs::run`.
        if emit_hoist_on() {
            return scan_strip_disc_csr::<O, PIECE, true>(x, disc, csr, job, from, out);
        }
        return scan_strip_disc_csr::<O, PIECE, false>(x, disc, csr, job, from, out);
    }
    let two_stage = x.two_stage;
    // (the presence vector by value: `lookup::XPvView`)
    let pv_view = two_stage.x_pv_view();
    let s_range = job.s_range;
    let kind = choose_disc_scan_subject(disc.template_type, disc.two_templates);
    let cut = Cut::new(&from);
    let start = from.s;
    let mut emit = |(s_off, index, index2): (i32, u32, Option<u32>)| {
        x_disc_word_pairs::<O, PIECE>(x, &cut, s_range, s_off, index, index2, out);
    };
    let mut fifo: XHeadFifo<(i32, u32, Option<u32>)> = XHeadFifo::new();
    disc_word_scan_subject_until(
        kind,
        x.search_seq_packed,
        disc.template_type,
        disc.second_template_type,
        disc.template_length as i32,
        [start as i32, job.end as i32],
        &mut |s_off, index, index2| {
            // A word whose chains are both empty (the presence vector, which `for_each_hit` and
            // `for_each_hit2` test first) has no pair: tested here, so the scan inlines the test
            // and such a word makes no call. A word with a chain has its chains' heads prefetched
            // and goes through the FIFO (`lookup::XHeadFifo`) to `x_disc_word_pairs`, out of
            // line, in the same order; a cut stops only there, which runs only in a push, and is
            // seen after that push as before. (A piece that stopped inside a word's chain resumes
            // at that word, which has a chain, so `Cut::chain_start` still meets it first.)
            if !pv_view.has_hits(u64::from(index))
                && index2.is_none_or(|index2| !pv_view.has_hits(u64::from(index2)))
            {
                return !PIECE || cut.go();
            }
            two_stage.x_prefetch_heads(u64::from(index), index2.map(u64::from));
            let left = fifo.push((s_off, index, index2), &mut emit);
            !PIECE || !left || cut.go()
        },
    );
    // (a stopped scan drops the words still in the FIFO: the next piece scans them again)
    if !PIECE || cut.go() {
        fifo.flush(&mut emit);
    }
    cut.stop.get()
}

/// No NCBI counterpart (call layout only): the pairs of one word of `scan_strip_disc_from` that
/// has a chain (its presence bit is set for either template), out of line, so the scan's closure
/// is the presence test and this call, and the scan inlines it: NCBI's `MB_ACCESS_HITS2` tests
/// each word's presence in the scan's loop and retrieves a chain only for a word with hits. The
/// words, their order, the pairs and every value are the same as with the walk in the closure.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1456-1469
/// (`MB_ACCESS_HITS2`, above `scan_strip_disc`).
#[inline(never)]
fn x_disc_word_pairs<O: SeedOut, const PIECE: bool>(
    x: &Inputs<'_>,
    cut: &Cut,
    s_range: usize,
    s_off: i32,
    index: u32,
    index2: Option<u32>,
    out: &mut O,
) {
    let two_stage = x.two_stage;
    let kmer_start = s_off as usize;
    // (where a resumed word's walk starts, and whether its first chain is passed over)
    let (node, second) = if PIECE { cut.chain_start() } else { (0, false) };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
    let mut pair = |q_off_1: u32| {
        if q_off_1 != 0 {
            if PIECE && !cut.take(kmer_start, q_off_1, false, out) {
                return;
            }
            seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
        }
    };
    if !PIECE {
        two_stage.for_each_hit(index as u64, &mut pair);
    } else if !second {
        two_stage.for_each_hit_from(index as u64, node, &mut pair);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1428-1444
    if let Some(index2) = index2 {
        let mut pair2 = |q_off_1: u32| {
            if PIECE && !cut.take(kmer_start, q_off_1, true, out) {
                return;
            }
            seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
        };
        if PIECE {
            let node2 = if second { node } else { 0 };
            two_stage.for_each_hit2_from(u64::from(index2), node2, &mut pair2);
        } else {
            two_stage.for_each_hit2(u64::from(index2), &mut pair2);
        }
    }
}

/// EXPERIMENT (LOSAT_X_DCCSR): `scan_strip_disc` over a table laid out as CSR
/// (`lookup::x_dc_csr`, one template): the same words in the same order, each `RING` words
/// late through `x_dc_csr::Ring` (whose prefetches bring in the cell's offsets and the ends of
/// its run), then the same presence test and the cell's run, which lists the chain's offsets in
/// the chain's order; so the same pairs reach `seed_of_pair` in the same order.
///
/// EXPERIMENT (LOSAT_X_DCEMITHOIST): with `HOIST` each word's run goes to `DcPairs::run` instead
/// of `seed_of_pair` pair by pair: for each word the same seeds in the same order (`DcPairs`); the
/// same cut of a piece (`Cut::fit`, and `Cut::run_cold` with `seed_of_pair` itself). A word in the
/// ring has passed the presence test (the scan's closure lets no other word in, and the ring hands
/// out only the words it took), so `emit` does not test it again.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
#[inline(never)]
fn scan_strip_disc_csr<O: SeedOut, const PIECE: bool, const HOIST: bool>(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    csr: &Csr,
    job: &StripJob,
    from: Resume,
    out: &mut O,
) -> Option<Resume> {
    debug_assert!(!disc.two_templates);
    if HOIST {
        EMIT_HOIST_STRIPS.fetch_add(1, Ordering::Relaxed);
    }
    let two_stage = x.two_stage;
    // (the presence vector by value: `lookup::XPvView`)
    let pv_view = two_stage.x_pv_view();
    let s_range = job.s_range;
    // EXPERIMENT (LOSAT_X_DCEMITHOIST): the strip's values of `DcPairs` and the context index's
    // form, made once
    let hoist = if HOIST {
        Some((DcPairs::new(x, s_range), x.query_context_index.x_view()))
    } else {
        None
    };
    let kind = choose_disc_scan_subject(disc.template_type, disc.two_templates);
    let cut = Cut::new(&from);
    let start = from.s;
    let mut emit = |(s_off, index): (i32, u32)| {
        if !HOIST && !two_stage.has_hits(u64::from(index)) {
            return;
        }
        debug_assert!(two_stage.has_hits(u64::from(index)));
        let kmer_start = s_off as usize;
        let mut run = csr.run(index as usize);
        if PIECE {
            match cut.fit(run, out) {
                Ok(fit) => run = fit,
                Err((skip, rest)) => {
                    cut.run_cold(kmer_start, skip, rest, out, |q_off_1, out| {
                        seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out)
                    });
                    return;
                }
            }
        }
        if let (true, Some((pairs, form))) = (HOIST, hoist) {
            pairs.run(form, run, kmer_start, out);
            return;
        }
        for &q_off_1 in run {
            seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
        }
    };
    let mut ring: Ring<'_, (i32, u32)> = Ring::new(csr);
    disc_word_scan_subject_until(
        kind,
        x.search_seq_packed,
        disc.template_type,
        disc.second_template_type,
        disc.template_length as i32,
        [start as i32, job.end as i32],
        // (a cut stops only in `emit`, so only a word that let one out can have stopped it)
        &mut |s_off, index, _index2| {
            // A word without hits (its presence bit, which `emit` tests first) does not enter the
            // ring: tested here, inline in the scan, so it costs neither the ring's prefetches nor
            // a place in it. The words with hits leave in the same order, so the same pairs reach
            // `seed_of_pair` in the same order; a cut stops only in `emit`, which runs only in a
            // push, and is seen after that push as before.
            if !pv_view.has_hits(u64::from(index)) {
                return true;
            }
            let left = ring.push((s_off, index), |w| w.1 as usize, &mut emit);
            !PIECE || !left || cut.go()
        },
    );
    // (a stopped scan drops the words still in the ring: the next piece scans them again)
    if !PIECE || cut.go() {
        ring.flush(&mut emit);
    }
    cut.stop.get()
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1106-1144
/// EXPERIMENT (LOSAT_X_DCEMITHOIST): `seed_of_pair` for the pairs of a strip of a discontiguous
/// table laid out as CSR, with what it reads for every pair read once per strip. For a pair
/// (`q_off_1`, `kmer_start`) of the strip it gives `seed_of_pair`'s seed, or none where that drops
/// the pair:
/// - `word` is `lut_word_length` and `word_length` at once: both are the template's length for a
///   table with templates (`TwoStageLookup::lut_word_length`, `word_length`), and the strips scan
///   only such a table (`x.two_stage.disc()`, whose CSR this is); so `seed_of_pair`'s word extension
///   (`word_length > lut_word_length`) never runs, and the seed is the pair itself (`q_off0`,
///   `kmer_start`), as there. `new` checks the two lengths once.
/// - the bounds test is `seed_of_pair`'s `||` of two tests without side effects; the subject's
///   (`kmer_start + word > s_len`) depends on the word alone and is made once per word, before any
///   pair of it (`run`).
/// - the context is `context_for_offset`'s (`QueryContextIndex::x_context_with` with the form made
///   once per strip: the same answer for every offset), and `contexts[context]` is read without the
///   bounds check: `new` checks that the index was made of as many contexts as there are, at least
///   one, and every answer of the index is below that number (`QueryContextIndex::x_contexts`).
/// - `s_range as u32` is made once per strip.
#[derive(Clone, Copy)]
struct DcPairs<'a> {
    word: usize,
    s_len: usize,
    s_range: u32,
    contexts: &'a [QueryContext],
    index: &'a QueryContextIndex,
}

impl<'a> DcPairs<'a> {
    /// The values for the pairs of a strip of `x` scanned with `s_range`.
    fn new(x: &Inputs<'a>, s_range: usize) -> Self {
        let word = x.two_stage.lut_word_length();
        assert_eq!(word, x.two_stage.word_length(), "a discontiguous word");
        let contexts = x.query_contexts;
        let index = x.query_context_index;
        assert!(!contexts.is_empty() && index.x_contexts() == contexts.len());
        DcPairs {
            word,
            s_len: x.s_len,
            s_range: s_range as u32,
            contexts,
            index,
        }
    }

    /// The seeds of the pairs `run` of the word at `kmer_start`, in order, into `out` (which has
    /// room for all of them), with the context form `form` (`QueryContextIndex::x_view`): one
    /// loop per form, its constants in registers.
    #[inline(always)]
    fn run<O: SeedOut>(self, form: Option<XView<'_>>, run: &[u32], kmer_start: usize, out: &mut O) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:730-733
        // (the subject's half of the bounds test, the same for every pair of the word)
        if kmer_start + self.word > self.s_len {
            return;
        }
        match form {
            Some(XView::Div(form)) => self.seeds(&form, run, kmer_start, out),
            Some(XView::Blocks(form)) => self.seeds(&form, run, kmer_start, out),
            Some(XView::BlocksRanks(form)) => self.seeds(&form, run, kmer_start, out),
            None => self.seeds(&XNone, run, kmer_start, out),
        }
    }

    /// `run` with the form `form`.
    #[inline(always)]
    fn seeds<G: XCtxGet, O: SeedOut>(self, form: &G, run: &[u32], kmer_start: usize, out: &mut O) {
        let DcPairs {
            word,
            s_range,
            contexts,
            index,
            ..
        } = self;
        let s = kmer_start as u32;
        out.x_push_run(run, |q_off_1| {
            let q_off0 = q_off_1 as usize - 1;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:730-733
            let context_idx = index.x_context_with(form, q_off0);
            debug_assert!(context_idx < contexts.len());
            // SAFETY: `context_idx` is an answer of the index, below the number of contexts it
            // was made of, which `new` checked to be `contexts.len()` (struct comment).
            let ctx = unsafe { contexts.get_unchecked(context_idx) };
            let q_context_start = ctx.query_offset as usize;
            let q_context_end = q_context_start + ctx.len;
            if q_off0 + word > q_context_end {
                return None;
            }
            Some(Seed {
                context: context_idx as u32,
                q: q_off0 as u32,
                s,
                s_range,
            })
        });
    }
}

/// No NCBI counterpart (the layout of the table only): where `fold_two_hit_on` keeps the
/// diagonals' inserts. `Latest` is the map of each diagonal's latest insert of LOSAT_X_DCLATEST
/// (`x_dc_latest.rs`); `Chains` (tests only) is NCBI's diagonal hash, on which the tests run
/// `fold_two_hit_on` against `fold_two_hit_t`. The fold makes the same statements in the same order
/// on either.
trait FoldTable {
    /// Whether the fold gives `ahead` the seed `latest::AHEAD` places after the one being folded
    /// (not for the chains, whose fold then has no code for it).
    const AHEAD: bool = false;

    /// A hint for the seed `latest::AHEAD` places after the one being folded.
    #[inline(always)]
    fn ahead(&self, _seed: &Seed) {}

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-836
    /// The retrieve of the seed's own diagonal at `s_off_pos`: `last_hit` and `hit_saved` (0 and
    /// false with no cell).
    fn retrieve(&mut self, diag: i32, s_off_pos: usize) -> (i32, bool);

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:861-885
    /// The retrieve of another diagonal in the off-diagonal loop (`last_hit`, `hit_len`,
    /// `hit_saved`), which runs no iteration in the scope (`scan_range == 0`).
    fn retrieve_off(&self, diag: i32) -> Option<(i32, i32, bool)>;

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
    /// The insert at the end of the seed (`part::insert_t`'s arguments).
    #[allow(clippy::too_many_arguments)]
    fn insert(
        &mut self,
        out: &mut part::PartOut,
        key: u64,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    );
}

/// The diagonal hash of a fold (NCBI's table) for the tests of `fold_two_hit_on`: with `PT` a
/// partition's table (LOSAT_X_PAIRPARFOLD), `LEAN` and `MEMO` as for `fold_two_hit`, `TRACK` as for
/// `fold_two_hit_t`; with LOSAT_X_DIAGLV's bounds (`x_lv`) the retrieve skips the walks below the
/// seed's bound.
#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Chains<'a, const PT: bool, const LEAN: bool, const MEMO: bool, const TRACK: bool> {
    hash: &'a mut DiagHashTable,
    x_lv: bool,
    x_lv_seeds: super::super::x_diag_lv::LvSeeds,
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<'a, const PT: bool, const LEAN: bool, const MEMO: bool, const TRACK: bool>
    Chains<'a, PT, LEAN, MEMO, TRACK>
{
    /// The table of a fold of the seeds of `x`.
    #[inline(always)]
    fn new(x: &Inputs<'_>, hash: &'a mut DiagHashTable) -> Self {
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the table has level bounds; the bound's
        // constants of these seeds.
        let x_lv = hash.x_lv_on();
        let x_lv_seeds = super::super::x_diag_lv::LvSeeds::new(
            x.window_size > 0,
            x.two_stage.word_length(),
            x.window_size,
        );
        Chains {
            hash,
            x_lv,
            x_lv_seeds,
        }
    }
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<const PT: bool, const LEAN: bool, const MEMO: bool, const TRACK: bool> FoldTable
    for Chains<'_, PT, LEAN, MEMO, TRACK>
{
    #[inline(always)]
    fn retrieve(&mut self, diag: i32, s_off_pos: usize) -> (i32, bool) {
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same decisions without the walk
        // of a bucket below the bound (the `last_hit` test and the window test).
        let (level, _hit_len, hit_saved) = if self.x_lv {
            self.hash
                .retrieve_below::<MEMO>(diag, &self.x_lv_seeds, s_off_pos)
        } else {
            self.hash.retrieve(diag)
        }
        .unwrap_or((0, 0, false));
        (level, hit_saved)
    }

    #[inline(always)]
    fn retrieve_off(&self, diag: i32) -> Option<(i32, i32, bool)> {
        self.hash.retrieve(diag)
    }

    #[inline(always)]
    fn insert(
        &mut self,
        out: &mut part::PartOut,
        key: u64,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) {
        part::insert_t::<PT, LEAN, MEMO, TRACK>(
            self.hash,
            out,
            key,
            diag,
            level,
            hit_len,
            hit_saved,
            s_off,
            window_size,
        );
    }
}

/// EXPERIMENT (LOSAT_X_DCLATEST, `x_dc_latest.rs`): the map of the latest insert of each diagonal
/// as a fold's table (the owner's, or a partition's), with the level bounds of the fold's seeds
/// (`x_diag_lv::LvSeeds`, the bound below which an entry is dead) and the place `find` gave the
/// seed's insert.
struct Latest<'a> {
    map: &'a mut DiagLatest,
    seeds: super::super::x_diag_lv::LvSeeds,
    found: Found,
}

impl<'a> Latest<'a> {
    /// The map of a fold of the seeds of `x` (in the scope: two hits and `scan_range == 0`, which
    /// the map needs, module comment of `x_dc_latest.rs` (i)).
    #[inline(always)]
    fn new(x: &Inputs<'_>, map: &'a mut DiagLatest) -> Self {
        debug_assert!(x.window_size > 0 && x.scan_range == 0);
        Latest {
            map,
            seeds: super::super::x_diag_lv::LvSeeds::new(
                true,
                x.two_stage.word_length(),
                x.window_size,
            ),
            found: Found::NONE,
        }
    }
}

impl FoldTable for Latest<'_> {
    const AHEAD: bool = true;

    #[inline(always)]
    fn ahead(&self, seed: &Seed) {
        self.map
            .prefetch((seed.s as isize - seed.q as isize) as i32);
    }

    #[inline(always)]
    fn retrieve(&mut self, diag: i32, s_off_pos: usize) -> (i32, bool) {
        self.found = self.map.find(diag, s_off_pos as i32);
        (self.found.level, self.found.saved)
    }

    fn retrieve_off(&self, _diag: i32) -> Option<(i32, i32, bool)> {
        unreachable!("LOSAT_X_DCLATEST: an off-diagonal retrieve with scan_range 0")
    }

    #[inline(always)]
    fn insert(
        &mut self,
        _out: &mut part::PartOut,
        _key: u64,
        diag: i32,
        level: i32,
        _hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        _window_size: i32,
    ) {
        let bound = self.seeds.bound(s_off as usize);
        self.map
            .put(self.found, diag, level, hit_saved, s_off, bound);
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-841
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:882-895,923-941
/// Phase 2 (owner only): the rest of `pair_extend::extend_offset_pairs` for the seeds of one
/// strip, in order, statement for statement for the diagonal hash in two-hit mode; the debug, trace
/// and diagonal-array branches are outside the scope and left out, and the small-query
/// word extension (`small_na_word`) is out of the scope, so `s_TypeOfWord` checks the query
/// masks only.
///
/// With `PT` (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`) the seeds are one partition's seeds
/// of one strip, `diag_hash` and `ungapped_hits` are the partition's, and `out` gets the key
/// (`out.strip`, `pos[k]`) of every pushed hit and of every insert that appends a cell, the
/// buckets the inserts change, and the time of the ungapped extensions (instead of the shared
/// LOSAT_TIMING counters); without it `out` and `pos` are not read.
///
/// With `LEAN` (LOSAT_X_FOLDLEAN, `lean_on`) the fold takes `type_of_word`'s early return
/// itself (NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:531, which
/// returns 1 when the word length equals the lookup table's word length, with `*extended`
/// set to 0 and the offsets unchanged) when the word lengths are equal, and inlines the insert
/// (`insert_inl`). `MEMO` says the table has the memo of LOSAT_X_DIAGLV2 (`x_last`); without it
/// the retrieve and the insert do not test for one, so with both off this is the fold before
/// those switches.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn fold_two_hit<const PT: bool, const LEAN: bool, const MEMO: bool>(
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
    fold_two_hit_t::<PT, LEAN, MEMO, true>(
        x,
        seeds,
        diag_hash,
        ungapped_hits,
        timing,
        scan_pause_ns,
        counts,
        out,
        pos,
    );
}

/// `fold_two_hit`; EXPERIMENT (LOSAT_X_DCSEEDLEAN): without `TRACK` the inserts of a partition
/// do not record the buckets they change (`part::insert_t`; the merge then writes back every
/// bucket of the partition). The fold on the chains, the code before LOSAT_X_DCLATEST, written out
/// here and not through `FoldTable`: every form of `fold_two_hit_on` on the chains that was tried
/// (a table struct holding the hash; the trait on `DiagHashTable` itself, with the fold's
/// constants in a struct or in the locals of this function) compiled to other code than this, and
/// the struct form cost about 4 % more cycles at 8 and 16 threads with LOSAT_X_DCLATEST=0 (FD-3,
/// 2026-11-24). The tests run `fold_two_hit_on` on the chains against this function (hits and
/// table).
#[allow(clippy::too_many_arguments)]
fn fold_two_hit_t<const PT: bool, const LEAN: bool, const MEMO: bool, const TRACK: bool>(
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
    let exact_ungapped = x.exact_ungapped;
    // EXPERIMENT (LOSAT_X_FOLDLEAN): `type_of_word` returns at its early return (na_ungapped.c:531).
    let x_tow = LEAN && word_length == lut_word_length;
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
        let diag = s_ext_start as isize - q_ext_start as isize;
        let (last_hit, hit_saved) = {
            // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same decisions without the walk
            // of a bucket below the bound (the `last_hit` test and the window test).
            let (level, _hit_len, hit_saved) = if x_lv {
                diag_hash.retrieve_below::<MEMO>(
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
        let s_end_pos_i32 = s_end_pos as i32;
        let window_end = last_hit + window_size as i32;
        if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
            let (word_type, extended, q_off_adj, s_off_adj) = if x_tow {
                (1, 0, q_off, s_off)
            } else {
                type_of_word(
                    q_off,
                    s_off,
                    !query_mask.is_empty(),
                    q_context_end,
                    s_range,
                    word_length,
                    lut_word_length,
                    true, // check for double-word extension
                    &mut is_seed_masked,
                )
            };
            q_off = q_off_adj;
            s_off = s_off_adj;

            // Skip this seed when `word_type` is zero (non-word).
            if word_type == 0 {
                continue;
            }

            // Advance `s_end` and `s_end_pos` by `extended`.
            s_end += extended;
            s_end_pos += extended;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:853-886
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
            let (wt, ext, q_off_adj, s_off_adj) = if x_tow {
                (1, 0, q_off, s_off)
            } else {
                type_of_word(
                    q_off,
                    s_off,
                    !query_mask.is_empty(),
                    q_context_end,
                    s_range,
                    word_length,
                    lut_word_length,
                    false, // skip double-word check (not in two-hit block)
                    &mut is_seed_masked,
                )
            };
            if wt == 0 {
                continue;
            }
            q_off = q_off_adj;
            s_off = s_off_adj;
            s_end += ext;
            s_end_pos += ext;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:897,939-941
        if !hit_ready {
            part::insert_t::<PT, LEAN, MEMO, TRACK>(
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:908-920
        let x_dropoff = x.x_dropoff_scores[context_idx];
        let reduced_cutoff = x.reduced_cutoff_scores[context_idx];
        let ungapped_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let ungapped = if exact_ungapped {
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
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:928-929
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:923-932,939-941
        if !(off_found || ungapped_score >= cutoff_score) {
            hit_ready = false;
            part::insert_t::<PT, LEAN, MEMO, TRACK>(
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:924-927
        super::super::x_push_ungapped(
            ungapped_hits,
            UngappedHit::saved(q_idx, ctx, q_off, s_off, &ungapped),
        );
        if PT {
            out.hit_keys.push(x_key);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
        part::insert_t::<PT, LEAN, MEMO, TRACK>(
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

/// EXPERIMENT (LOSAT_X_DCLATEST): `fold_two_hit` on the map of the latest inserts `map` (the
/// owner's, or with `PT` a partition's): the inserts change no `DiagHashTable`, and the fold brings
/// in the map's bucket of the seed `latest::AHEAD` places ahead. The function of the fold on a
/// map, `map` its argument (`fold_two_hit_on` inlined with `Latest`).
#[allow(clippy::too_many_arguments)]
fn fold_two_hit_latest<const PT: bool, const LEAN: bool>(
    x: &Inputs<'_>,
    seeds: &[Seed],
    map: &mut DiagLatest,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    counts: &mut FoldCounts,
    out: &mut part::PartOut,
    pos: &[u32],
) {
    fold_two_hit_on::<PT, LEAN, _>(
        x,
        seeds,
        &mut Latest::new(x, map),
        ungapped_hits,
        timing,
        scan_pause_ns,
        counts,
        out,
        pos,
    );
}

/// The body of `fold_two_hit_latest` on the table `table`, inlined into it (so that the table's
/// fields and the place `find` gave a seed stay in registers): the statements of `fold_two_hit_t`
/// in the same order, through `FoldTable` (the tests run it on the chains, `Chains`, against
/// `fold_two_hit_t`).
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn fold_two_hit_on<const PT: bool, const LEAN: bool, T: FoldTable>(
    x: &Inputs<'_>,
    seeds: &[Seed],
    table: &mut T,
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
    let exact_ungapped = x.exact_ungapped;
    // EXPERIMENT (LOSAT_X_FOLDLEAN): `type_of_word` returns at its early return (na_ungapped.c:531).
    let x_tow = LEAN && word_length == lut_word_length;
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let diag_offset = x.diag_offset;
    let diag_hash_window = x.diag_hash_window;
    let window_size = x.window_size;
    let scan_range = x.scan_range;
    let two_hits = window_size > 0;
    let timing_enabled = timing.is_some();
    counts.seeds += seeds.len() as u64;
    for (x_k, seed) in seeds.iter().enumerate() {
        // EXPERIMENT (LOSAT_X_DCLATEST): the map's bucket of a later seed (a hint only).
        if T::AHEAD {
            if let Some(later) = seeds.get(x_k + latest::AHEAD) {
                table.ahead(later);
            }
        }
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
        let diag = s_ext_start as isize - q_ext_start as isize;
        let s_off_pos = s_ext_start + diag_offset as usize;
        let s_off_pos_i32 = s_off_pos as i32;
        let (last_hit, hit_saved) = table.retrieve(diag as i32, s_off_pos);

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
        let s_end_pos_i32 = s_end_pos as i32;
        let window_end = last_hit + window_size as i32;
        if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
            let (word_type, extended, q_off_adj, s_off_adj) = if x_tow {
                (1, 0, q_off, s_off)
            } else {
                type_of_word(
                    q_off,
                    s_off,
                    !query_mask.is_empty(),
                    q_context_end,
                    s_range,
                    word_length,
                    lut_word_length,
                    true, // check for double-word extension
                    &mut is_seed_masked,
                )
            };
            q_off = q_off_adj;
            s_off = s_off_adj;

            // Skip this seed when `word_type` is zero (non-word).
            if word_type == 0 {
                continue;
            }

            // Advance `s_end` and `s_end_pos` by `extended`.
            s_end += extended;
            s_end_pos += extended;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:853-886
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
                        table.retrieve_off(diag_i32 + delta as i32)
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
                        table.retrieve_off(diag_i32 - delta as i32)
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
            let (wt, ext, q_off_adj, s_off_adj) = if x_tow {
                (1, 0, q_off, s_off)
            } else {
                type_of_word(
                    q_off,
                    s_off,
                    !query_mask.is_empty(),
                    q_context_end,
                    s_range,
                    word_length,
                    lut_word_length,
                    false, // skip double-word check (not in two-hit block)
                    &mut is_seed_masked,
                )
            };
            if wt == 0 {
                continue;
            }
            q_off = q_off_adj;
            s_off = s_off_adj;
            s_end += ext;
            s_end_pos += ext;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:897,939-941
        if !hit_ready {
            table.insert(
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:908-920
        let x_dropoff = x.x_dropoff_scores[context_idx];
        let reduced_cutoff = x.reduced_cutoff_scores[context_idx];
        let ungapped_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let ungapped = if exact_ungapped {
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
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:928-929
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:923-932,939-941
        if !(off_found || ungapped_score >= cutoff_score) {
            hit_ready = false;
            table.insert(
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:924-927
        super::super::x_push_ungapped(
            ungapped_hits,
            UngappedHit::saved(q_idx, ctx, q_off, s_off, &ungapped),
        );
        if PT {
            out.hit_keys.push(x_key);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
        table.insert(
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
/// The seed stage of one chunk on `threads` threads (this one and `threads - 1` helpers of
/// its pool): the scan ranges of `scan_subject_disc_words_with_ranges` (the whole subject, or
/// a masked subject's `seq_ranges` one by one; stride 1) cut into strips scanned by any
/// thread, folded here in strip order into `hash` and `hits`. With `fold_part`
/// (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`) the fold is split by diagonal-hash bucket over
/// the threads instead.
///
/// EXPERIMENT (LOSAT_X_FOLDLEAN, LOSAT_X_DIAGLV2): one instance of the stage per pair of
/// switches (`seed_stage_of`), the table's memo deciding `MEMO`, so that with both off the
/// stage is the code before them.
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
    // EXPERIMENT (LOSAT_X_DIAGLV2): the memo, kept by the strips folds only, made again from
    // the bounds (the reference fold may have changed the table since the last stage).
    let memo = !hash.x_last.is_empty();
    if memo {
        hash.x_lv2_reset();
    }
    let lists = seed_lean_on();
    match (lean_on(), memo) {
        (false, false) => seed_stage_of::<false, false, false>(
            x, disc, threads, hash, hits, timing, timed, fold_part, lists, None,
        ),
        (false, true) => seed_stage_of::<false, true, false>(
            x, disc, threads, hash, hits, timing, timed, fold_part, lists, None,
        ),
        (true, false) => seed_stage_of::<true, false, false>(
            x, disc, threads, hash, hits, timing, timed, fold_part, lists, None,
        ),
        (true, true) => seed_stage_of::<true, true, false>(
            x, disc, threads, hash, hits, timing, timed, fold_part, lists, None,
        ),
    }
}

/// EXPERIMENT (LOSAT_X_DCLATEST): `seed_stage` on maps of the latest inserts (`x_dc_latest.rs`)
/// instead of the diagonal hash `hash`, which the stage neither reads nor changes: with `parts`
/// (LOSAT_X_PAIRPARFOLD) each partition's map (`part::PartState::latest`), else the owner's
/// (`PartScratch::latest`); each starts the stage at `bound0`, the level bound of the stage's
/// first position (`latest_bound`). No memo (LOSAT_X_DIAGLV2) is made. `lean`: LOSAT_X_FOLDLEAN,
/// `lists`: LOSAT_X_DCSEEDLEAN's partitioned stage (`run` passes the switches).
#[allow(clippy::too_many_arguments)]
fn seed_stage_latest(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
    scratch: &mut PartScratch,
    parts: Option<usize>,
    bound0: i32,
    lean: bool,
    lists: bool,
) -> StageStats {
    let (fold_part, owner) = match parts {
        Some(parts) => (
            Some(part::FoldPart {
                scratch,
                parts,
                latest: Some(bound0),
            }),
            None,
        ),
        None => {
            scratch.latest.begin(bound0);
            (None, Some(&mut scratch.latest))
        }
    };
    if lean {
        seed_stage_of::<true, false, true>(
            x, disc, threads, hash, hits, timing, timed, fold_part, lists, owner,
        )
    } else {
        seed_stage_of::<false, false, true>(
            x, disc, threads, hash, hits, timing, timed, fold_part, lists, owner,
        )
    }
}

/// `seed_stage` for one pair of LOSAT_X_FOLDLEAN (`LEAN`) and the table's memo of
/// LOSAT_X_DIAGLV2 (`MEMO`); `lists`: LOSAT_X_DCSEEDLEAN's partitioned stage (`seed_stage_lists`).
/// EXPERIMENT (LOSAT_X_DCLATEST): with `MAP` the folds keep the maps of the latest inserts
/// (`seed_stage_latest`: `owner` without partitions), not `hash`.
#[allow(clippy::too_many_arguments)]
fn seed_stage_of<const LEAN: bool, const MEMO: bool, const MAP: bool>(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
    fold_part: Option<part::FoldPart<'_>>,
    lists: bool,
    mut owner: Option<&mut DiagLatest>,
) -> StageStats {
    let jobs = plan_strips(
        x.subject_seq_ranges,
        x.s_len,
        disc.template_length,
        disc.template_length,
        1,
        x.subject_masked,
        STRIPS_PER_THREAD * threads.max(4),
    );
    let produce = |i: usize, out: &mut Vec<Seed>| scan_strip_disc(x, disc, &jobs[i], out);
    if let Some(fold_part) = fold_part {
        let parts = fold_part.parts;
        // EXPERIMENT (LOSAT_X_DCSEEDLEAN): every partitioned stage takes the partitions' lists and
        // the folds that do not record their buckets.
        if lists {
            crate::utils::xreach::x_reach!("LOSAT_X_DCSEEDLEAN", true);
            let r = seed_stage_lists::<LEAN, MEMO, MAP>(
                x, disc, &jobs, threads, hash, hits, timing, timed, fold_part,
            );
            part::record(part::TWO_HITS, parts, &r);
            return StageStats::of_part(jobs.len(), r, timing);
        }
        let fold_one = |st: &mut part::PartState, seeds: &[Seed], pos: &[u32]| {
            let part::PartState {
                hash,
                hits,
                out,
                counts,
                scan_pause_ns,
                latest: map,
                ..
            } = st;
            if MAP {
                // EXPERIMENT (LOSAT_X_DCLATEST): the partition's map.
                fold_two_hit_latest::<true, LEAN>(
                    x,
                    seeds,
                    map,
                    hits,
                    timing,
                    scan_pause_ns,
                    counts,
                    out,
                    pos,
                );
            } else {
                fold_two_hit::<true, LEAN, MEMO>(
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
            }
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
            if MAP {
                // EXPERIMENT (LOSAT_X_DCLATEST): the owner's map.
                let map = owner
                    .as_deref_mut()
                    .expect("LOSAT_X_DCLATEST: a stage without the owner's map");
                fold_two_hit_latest::<false, LEAN>(
                    x,
                    seeds,
                    map,
                    hits,
                    timing,
                    scan_pause_ns,
                    counts,
                    &mut part::PartOut::default(),
                    &[],
                );
            } else {
                fold_two_hit::<false, LEAN, MEMO>(
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
            }
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

/// No NCBI counterpart (scheduling only): EXPERIMENT (LOSAT_X_DCSEEDLEAN) `seed_stage_of`'s
/// partitioned stage (LOSAT_X_PAIRPARFOLD) on the strips `jobs`: the home loop's strips are
/// scanned into their partitions' lists (`scan_strip_disc` into a `part::PartSink`) instead of
/// into one list cut by partition, and the folds do not record the buckets they change
/// (`fold_two_hit_t` without `TRACK`), the merge writing back every bucket of each partition
/// (`part::seed_stage_part_lists`). Each partition folds the same seeds in the same order.
/// EXPERIMENT (LOSAT_X_DCLATEST): with `MAP` into the partitions' maps (`seed_stage_latest`).
#[allow(clippy::too_many_arguments)]
fn seed_stage_lists<const LEAN: bool, const MEMO: bool, const MAP: bool>(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    jobs: &[StripJob],
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
    fold_part: part::FoldPart<'_>,
) -> part::PartStats {
    let produce = |i: usize, out: &mut Vec<Seed>| scan_strip_disc(x, disc, &jobs[i], out);
    let produce_lists =
        |i: usize, out: &mut part::PartSink<'_>| scan_strip_disc(x, disc, &jobs[i], out);
    // EXPERIMENT (LOSAT_X_SEEDPIECE): a piece of strip i, from its start or a stop.
    let produce_piece = |i: usize, from: Option<Resume>, out: &mut part::PieceSink<'_>| {
        let job = &jobs[i];
        scan_strip_disc_from::<_, true>(x, disc, job, from.unwrap_or(Resume::at(job.start)), out)
    };
    let fold_one = |st: &mut part::PartState, seeds: &[Seed], pos: &[u32]| {
        let part::PartState {
            hash,
            hits,
            out,
            counts,
            scan_pause_ns,
            latest: map,
            ..
        } = st;
        if MAP {
            // EXPERIMENT (LOSAT_X_DCLATEST): the partition's map.
            fold_two_hit_latest::<true, LEAN>(
                x,
                seeds,
                map,
                hits,
                timing,
                scan_pause_ns,
                counts,
                out,
                pos,
            );
        } else {
            fold_two_hit_t::<true, LEAN, MEMO, false>(
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
        }
    };
    part::seed_stage_part_lists(
        jobs.len(),
        threads,
        AHEAD_PER_THREAD * threads,
        hash,
        hits,
        fold_part,
        &produce,
        &produce_lists,
        &produce_piece,
        &fold_one,
        timed,
    )
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// The parent's `run` for a discontiguous table: the seed stage of this chunk (the loop
/// above with its two-hit extension) when LOSAT_X_PAIRPARDISC is on and the chunk is in
/// scope; `true` means it has run and the reference scan is skipped. `fold_scratch` holds the
/// allocations of LOSAT_X_PAIRPARFOLD's partition tables and buffers (`x_pair_seed_part.rs`).
pub(super) fn run(
    x: &Inputs<'_>,
    diag_hash: &mut DiagHashTable,
    diag_array: DiagArray<'_>,
    ungapped_hits: &mut Vec<UngappedHit>,
    fold_scratch: &mut PartScratch,
) -> bool {
    // EXPERIMENT (LOSAT_X_PAIRPARDISCARR): a query block of at most 8000 letters has the
    // diagonal array, which the stage below never takes.
    if x.use_array_indexing {
        return run_array(x, diag_array, ungapped_hits);
    }
    if !on() {
        return false;
    }
    let Some(disc) = x.two_stage.disc() else {
        return false;
    };
    if !in_scope(x) {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARDISC", false);
        reach_lean(false);
        reach_latest(false);
        return false;
    }
    let threads = match pool_threads() {
        Some(threads) => threads,
        // LOSAT_X_BNDCONE: the stage runs without a pool
        None if one_thread_on() => {
            crate::utils::xreach::x_reach!("LOSAT_X_BNDCONE", true);
            1
        }
        None => {
            crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARDISC", false);
            reach_lean(false);
            reach_latest(false);
            return false;
        }
    };
    crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARDISC", true);
    reach_lean(true);
    // EXPERIMENT (LOSAT_X_DCLATEST): the folds keep maps of the latest inserts (every chunk of
    // the stage has the bounds, `latest_bound`), unless the chunk's scan ranges step back
    // (`latest_ranges_ascend`: those go to the chains).
    let latest = if latest_on() {
        latest_bound(x).filter(|_| {
            let ascend = latest_ranges_ascend(x, disc.template_length);
            latest::count_stage(ascend);
            ascend
        })
    } else {
        None
    };
    reach_latest(latest.is_some());
    if super::super::x_diag_lv::lv2_on() {
        // (the maps keep no memo)
        crate::utils::xreach::x_reach!(
            "LOSAT_X_DIAGLV2",
            latest.is_none() && !diag_hash.x_last.is_empty()
        );
    }
    let timing = x.timing;
    let timed = timing.is_some() || timing_env();
    let t_start = timed.then(std::time::Instant::now);
    // EXPERIMENT (LOSAT_X_PAIRPARFOLD): the fold split over the pool's threads.
    let parts = part::partitions(threads);
    let stats = match latest {
        Some(bound0) => seed_stage_latest(
            x,
            disc,
            threads,
            diag_hash,
            ungapped_hits,
            timing,
            timed,
            fold_scratch,
            parts,
            bound0,
            lean_on(),
            seed_lean_on(),
        ),
        None => seed_stage(
            x,
            disc,
            threads,
            diag_hash,
            ungapped_hits,
            timing,
            timed,
            parts.map(|parts| part::FoldPart {
                scratch: fold_scratch,
                parts,
                latest: None,
            }),
        ),
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
    true
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:660-677
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:768-772
/// EXPERIMENT (LOSAT_X_PAIRPARDISCARR), phase 2 on the owner thread: the seeds of one strip
/// through the two-hit diagonal-array branch of `extend_offset_pairs` (`pair_extend.rs`,
/// `s_BlastnDiagTableExtendInitialHit`), statement for statement: `fold_two_hit` with the
/// table reads and writes of the array (the diagonal masked by `diag_mask`, the bounds test
/// against `diag_array_size`, `last_hit` and `flag` in `hit_level_array`, the length in
/// `hit_len_array` as a byte) instead of the diagonal hash. The owner folds the strips in scan
/// order, so the array and the hit list are the reference's whatever the number of threads.
///
/// With `LEAN` (LOSAT_X_FOLDLEAN) it takes `type_of_word`'s early return itself as
/// `fold_two_hit` does.
fn fold_two_hit_array<const LEAN: bool>(
    x: &Inputs<'_>,
    seeds: &[Seed],
    arr: &mut DiagArray<'_>,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    counts: &mut FoldCounts,
) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let word_length = two_stage.word_length();
    let exact_ungapped = x.exact_ungapped;
    // EXPERIMENT (LOSAT_X_FOLDLEAN): `type_of_word` returns at its early return (na_ungapped.c:531).
    let x_tow = LEAN && word_length == lut_word_length;
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let diag_offset = x.diag_offset;
    let window_size = x.window_size;
    let scan_range = x.scan_range;
    let two_hits = window_size > 0;
    let timing_enabled = timing.is_some();
    let diag_array_length = arr.length;
    let diag_mask = arr.mask;
    let diag_array_size = arr.size;
    let hit_level_array: &mut [DiagStruct] = arr.levels;
    let hit_len_array: &mut [u8] = arr.lens;
    counts.seeds += seeds.len() as u64;
    for seed in seeds {
        let context_idx = seed.context as usize;
        let ctx = &x.query_contexts[context_idx];
        let q_idx = context_idx as u32;
        let q_context_start = ctx.query_offset as usize;
        let q_context_end = q_context_start + ctx.len;
        let cutoff_score = x.cutoff_scores[context_idx];
        let q_ext_start = seed.q as usize;
        let s_ext_start = seed.s as usize;
        let s_range = seed.s_range as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:663-668
        let diag = (s_ext_start as isize + diag_array_length - q_ext_start as isize) & diag_mask;
        let diag_idx = diag as usize;
        let (last_hit, hit_saved) = if diag_idx < diag_array_size {
            let diag_entry = &hit_level_array[diag_idx];
            (diag_entry.last_hit, diag_entry.flag != 0)
        } else {
            (0, false)
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:674-686
        let s_end_pos_i32 = s_end_pos as i32;
        let window_end = last_hit + window_size as i32;
        if two_hits && (hit_saved || s_end_pos_i32 > window_end) {
            let (word_type, extended, q_off_adj, s_off_adj) = if x_tow {
                (1, 0, q_off, s_off)
            } else {
                type_of_word(
                    q_off,
                    s_off,
                    !query_mask.is_empty(),
                    q_context_end,
                    s_range,
                    word_length,
                    lut_word_length,
                    true, // check for double-word extension
                    &mut is_seed_masked,
                )
            };
            q_off = q_off_adj;
            s_off = s_off_adj;

            // Skip this seed when `word_type` is zero (non-word).
            if word_type == 0 {
                continue;
            }

            // Advance `s_end` and `s_end_pos` by `extended`.
            s_end += extended;
            s_end_pos += extended;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:686-720
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
                let orig_diag = diag + diag_array_length;
                for delta in 1..=delta_max {
                    let off_diag = (orig_diag + delta) & diag_mask;
                    let (off_s_end, off_s_l) = {
                        let off_diag_idx = off_diag as usize;
                        if off_diag_idx < diag_array_size {
                            (
                                hit_level_array[off_diag_idx].last_hit,
                                hit_len_array[off_diag_idx],
                            )
                        } else {
                            (0, 0)
                        }
                    };
                    if off_s_l > 0
                        && (off_s_end as isize - delta) >= s_a
                        && (off_s_end as isize - off_s_l as isize) <= s_b
                    {
                        off_found = true;
                        break;
                    }
                    let off_diag = (orig_diag - delta) & diag_mask;
                    let (off_s_end, off_s_l) = {
                        let off_diag_idx = off_diag as usize;
                        if off_diag_idx < diag_array_size {
                            (
                                hit_level_array[off_diag_idx].last_hit,
                                hit_len_array[off_diag_idx],
                            )
                        } else {
                            (0, 0)
                        }
                    };
                    if off_s_l > 0
                        && (off_s_end as isize) >= s_a
                        && (off_s_end as isize - off_s_l as isize + delta) <= s_b
                    {
                        off_found = true;
                        break;
                    }
                }
                if !off_found {
                    hit_ready = false;
                }
            }
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:718-730
            let (wt, ext, q_off_adj, s_off_adj) = if x_tow {
                (1, 0, q_off, s_off)
            } else {
                type_of_word(
                    q_off,
                    s_off,
                    !query_mask.is_empty(),
                    q_context_end,
                    s_range,
                    word_length,
                    lut_word_length,
                    false, // skip double-word check (not in two-hit block)
                    &mut is_seed_masked,
                )
            };
            if wt == 0 {
                continue;
            }
            q_off = q_off_adj;
            s_off = s_off_adj;
            s_end += ext;
            s_end_pos += ext;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:728-772
        if !hit_ready {
            if diag_idx < diag_array_size {
                hit_level_array[diag_idx].last_hit = s_end_pos as i32;
                hit_level_array[diag_idx].flag = 0;
                if two_hits {
                    hit_len_array[diag_idx] = (s_end_pos - s_off_pos) as u8;
                }
            }
            continue;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:740-756
        let x_dropoff = x.x_dropoff_scores[context_idx];
        let reduced_cutoff = x.reduced_cutoff_scores[context_idx];
        let ungapped_start = if timing_enabled {
            Some(std::time::Instant::now())
        } else {
            None
        };
        let ungapped = if exact_ungapped {
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
            if let Some(timing) = timing {
                timing.ungapped_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
                timing.ungapped_calls.fetch_add(1, Ordering::Relaxed);
            }
            *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
        }
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:757-758
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:752-772
        if !(off_found || ungapped_score >= cutoff_score) {
            hit_ready = false;
            if diag_idx < diag_array_size {
                hit_level_array[diag_idx].last_hit = s_end_pos as i32;
                hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                if two_hits {
                    hit_len_array[diag_idx] = if hit_ready {
                        0
                    } else {
                        (s_end_pos - s_off_pos) as u8
                    };
                }
            }
            continue;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:753-756
        super::super::x_push_ungapped(
            ungapped_hits,
            UngappedHit::saved(q_idx, ctx, q_off, s_off, &ungapped),
        );

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:768-772
        if diag_idx < diag_array_size {
            hit_level_array[diag_idx].last_hit = ungapped_s_end_pos as i32;
            hit_level_array[diag_idx].flag = 1;
            if two_hits {
                hit_len_array[diag_idx] = 0;
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// EXPERIMENT (LOSAT_X_PAIRPARDISCARR): the seed stage of one chunk with the diagonal array on
/// `threads` threads: the strips of `seed_stage` scanned by any thread, folded here in strip
/// order into the array and `hits` (`fold_two_hit_array`). EXPERIMENT (LOSAT_X_FOLDLEAN): one
/// instance of the stage per state of the switch (`seed_stage_array_of`), so that with it off
/// the stage is the code before it.
fn seed_stage_array(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    threads: usize,
    arr: &mut DiagArray<'_>,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
) -> StageStats {
    if lean_on() {
        seed_stage_array_of::<true>(x, disc, threads, arr, hits, timing, timed)
    } else {
        seed_stage_array_of::<false>(x, disc, threads, arr, hits, timing, timed)
    }
}

/// `seed_stage_array` for one state of LOSAT_X_FOLDLEAN (`LEAN`).
fn seed_stage_array_of<const LEAN: bool>(
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    threads: usize,
    arr: &mut DiagArray<'_>,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
) -> StageStats {
    let jobs = plan_strips(
        x.subject_seq_ranges,
        x.s_len,
        disc.template_length,
        disc.template_length,
        1,
        x.subject_masked,
        STRIPS_PER_THREAD * threads.max(4),
    );
    let produce = |i: usize, out: &mut Vec<Seed>| scan_strip_disc(x, disc, &jobs[i], out);
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
            fold_two_hit_array::<LEAN>(x, seeds, arr, hits, timing, scan_pause_ns, counts);
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

/// EXPERIMENT (LOSAT_X_PAIRPARDISCARR): `run` for a chunk
/// with the diagonal array (a query block of at most 8000 letters): the seed stage of the
/// chunk on the strips with the owner's array fold when the switch is on and the chunk is in
/// scope (`in_scope_array`); `true` means it has run and the reference scan is skipped.
fn run_array(
    x: &Inputs<'_>,
    mut diag_array: DiagArray<'_>,
    ungapped_hits: &mut Vec<UngappedHit>,
) -> bool {
    if !arr_on() {
        return false;
    }
    let Some(disc) = x.two_stage.disc() else {
        return false;
    };
    if !in_scope_array(x) {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARDISCARR", false);
        reach_lean(false);
        return false;
    }
    let threads = match pool_threads() {
        Some(threads) => threads,
        // LOSAT_X_BNDCONE: the stage runs without a pool
        None if one_thread_on() => {
            crate::utils::xreach::x_reach!("LOSAT_X_BNDCONE", true);
            1
        }
        None => {
            crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARDISCARR", false);
            reach_lean(false);
            return false;
        }
    };
    crate::utils::xreach::x_reach!("LOSAT_X_PAIRPARDISCARR", true);
    reach_lean(true);
    let timing = x.timing;
    let timed = timing.is_some() || timing_env();
    let t_start = timed.then(std::time::Instant::now);
    let stats = seed_stage_array(
        x,
        disc,
        threads,
        &mut diag_array,
        ungapped_hits,
        timing,
        timed,
    );
    if let (Some(t_start), Some(timing)) = (t_start, timing) {
        // `[TIMING] scan_lookup`: the time of the stage minus the ungapped extensions, as the
        // reference counts its scan.
        let elapsed_ns = t_start.elapsed().as_nanos() as u64;
        let scan_ns = elapsed_ns.saturating_sub(stats.scan_pause_ns);
        timing.scan_ns.fetch_add(scan_ns, Ordering::Relaxed);
        timing.scan_calls.fetch_add(1, Ordering::Relaxed);
    }
    ARR_CHUNKS.fetch_add(1, Ordering::Relaxed);
    ARR_STRIPS.fetch_add(stats.strips as u64, Ordering::Relaxed);
    ARR_SEEDS.fetch_add(stats.counts.seeds, Ordering::Relaxed);
    ARR_SKIPPED.fetch_add(stats.counts.skipped, Ordering::Relaxed);
    ARR_FOLD_NS.fetch_add(stats.fold_ns, Ordering::Relaxed);
    ARR_OWNER_SCAN_NS.fetch_add(stats.owner_scan_ns, Ordering::Relaxed);
    ARR_WAIT_NS.fetch_add(stats.wait_ns, Ordering::Relaxed);
    true
}

/// No NCBI counterpart: with LOSAT_TIMING the counters of LOSAT_X_PAIRPARDISC and
/// LOSAT_X_PAIRPARDISCARR, on stderr at the end of the search.
pub(super) fn summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if on() && std::env::var_os("LOSAT_TIMING").is_some() {
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
        if latest_on() {
            eprintln!("[TIMING] x_dc_latest:{}", latest::timing_fields());
        }
    }
    if (on() || arr_on()) && std::env::var_os("LOSAT_TIMING").is_some() {
        eprintln!(
            "[TIMING] x_dc_emit_hoist: strips={}",
            load(&EMIT_HOIST_STRIPS)
        );
    }
    if arr_on() && std::env::var_os("LOSAT_TIMING").is_some() {
        let s = |c: &AtomicU64| c.load(Ordering::Relaxed) as f64 / 1e9;
        eprintln!(
            "[TIMING] x_pairpardiscarr_seed: chunks={} strips={} seeds={} skipped={} fold={:.3}s owner_scan={:.3}s owner_wait={:.3}s",
            load(&ARR_CHUNKS),
            load(&ARR_STRIPS),
            load(&ARR_SEEDS),
            load(&ARR_SKIPPED),
            s(&ARR_FOLD_NS),
            s(&ARR_OWNER_SCAN_NS),
            s(&ARR_WAIT_NS)
        );
    }
}

#[cfg(test)]
#[path = "x_pair_seed_disc_tests.rs"]
mod tests;
