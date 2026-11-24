//! EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable): the one-hit seed stage of one
//! subject chunk (scan, lookup chains, word extension, diagonal hash, ungapped extension) with
//! the scan of the subject cut into strips that the threads of the search pool scan, and the
//! diagonal-hash fold done by the owner thread in the original scan order.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-838,939-941
//!
//! NCBI scans the subject in pieces already: `scansub` stops when the offset array is full
//! and the loop above calls it again from `scan_range[1]`, so a scan kernel restarts at any
//! position of its grid `start + k * scan_step` and derives its state from that position.
//! This module cuts the grid of each scan range into strips (`plan_strips`). Phase 1, on
//! any thread: the reference kernel scans one strip with `s_range` of the whole range (the
//! value `scan_range[2] + lut_word_length` NCBI passes), walks the lookup chains as `on_word`
//! of the reference scan (`two_stage_scan::scan_two_stage_lookup`) does, and runs the steps of
//! `pair_extend::extend_offset_pairs` that read nothing
//! but the inputs (the context of the query offset, the bounds test, the word extension of
//! `s_BlastNaExtend`, which `extend_offset_pairs` calls as `word_extend_pair`); a pair that
//! passes becomes a `Seed`, in the order of the pairs.
//! Phase 2, on the owner thread only: the seeds of strip 0, 1, 2, ... in that order go
//! through a statement-for-statement copy of the rest of `extend_offset_pairs` for one-hit
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
//! megablast table, no small-query word extension, no debug, trace or window diagnostics, a
//! subject (masked or not: the strips cover its scan ranges, `plan_strips`) searched by the owner of a pool of two threads or more by itself
//! (one subject, or several one at a time under LOSAT_X_BNFEWSUBJ).
//! A chunk with a discontiguous megablast table goes to the child module `disc`
//! (`x_pair_seed_disc.rs`, LOSAT_X_PAIRPARDISC: the two-hit stage of dc-megablast on the same
//! strips, with its own switches; LOSAT_X_PAIRPARDISCARR: the same with the diagonal array).
//! A chunk searched with a table that is not two-stage (the direct table of word sizes 4 to
//! 8, `-task blastn-short`) goes to the child module `na` (`x_pair_seed_na.rs`,
//! LOSAT_X_PAIRPARNA, called by `run.rs` as `run_na`).

use super::DiagStruct;
use super::{
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, first_scan_range,
    mask_lookup_index, packed_kmer_at_seed_mask, scan_subject_kmers_range,
    scan_subject_kmers_range_mb_10_1, scan_subject_kmers_range_mb_10_2,
    scan_subject_kmers_range_mb_10_3, scan_subject_kmers_range_mb_11_1mod4,
    scan_subject_kmers_range_mb_11_2mod4, scan_subject_kmers_range_mb_11_3mod4,
    scan_subject_kmers_range_mb_9_1, scan_subject_kmers_range_mb_9_2,
    scan_subject_kmers_range_mb_any, scan_subject_kmers_with_offsets, select_mb_scan_kind,
    type_of_word, BlastnTiming, DiagHashTable, MbScanSubjectKind, QueryContext, QueryContextIndex,
    UngappedHit, COMPRESSION_RATIO,
};
use crate::algorithm::blastn::lookup::TwoStageLookup;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::atomic::{AtomicU64, Ordering};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::sync::Mutex;
use std::sync::OnceLock;

#[path = "x_pair_seed_disc.rs"]
mod disc;

// EXPERIMENT (LOSAT_X_DCLATEST): the map of each diagonal's latest insert that dc-megablast's
// strips folds keep instead of the diagonal hash (`disc`, `part`).
#[path = "x_dc_latest.rs"]
mod latest;

#[path = "x_pair_seed_na.rs"]
mod na;
pub(super) use na::{run as run_na, wants as wants_na, NaInputs};

// EXPERIMENT (LOSAT_X_PAIRPARFOLD): the fold of both strips paths split by diagonal-hash bucket
// over the pool (`x_pair_seed_part.rs`); its scratch is a field of `SubjectScratch`.
#[path = "x_pair_seed_part.rs"]
mod part;
pub(super) use part::PartScratch;

// No NCBI counterpart: counters of LOSAT_TIMING (scheduling only).
static X_CHUNKS: AtomicU64 = AtomicU64::new(0);
static X_STRIPS: AtomicU64 = AtomicU64::new(0);
static X_SEEDS: AtomicU64 = AtomicU64::new(0);
static X_SKIPPED: AtomicU64 = AtomicU64::new(0);
static X_FOLD_NS: AtomicU64 = AtomicU64::new(0);
static X_OWNER_SCAN_NS: AtomicU64 = AtomicU64::new(0);
static X_WAIT_NS: AtomicU64 = AtomicU64::new(0);

/// The inputs of the seed stage of one subject chunk, borrowed from
/// `prelim_subject::prelim_chunk_hits` (the same names as there).
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
    pub(super) cutoff_scores: &'a [i32],
    pub(super) x_dropoff_scores: &'a [i32],
    pub(super) reduced_cutoff_scores: &'a [i32],
    pub(super) score_matrix: &'a [i32; 256],
    pub(super) nucl_score_table: &'a [i32; 256],
    /// `ScanBatch::exact_ungapped` (matrix-only scoring or a word shorter than 11 letters).
    pub(super) exact_ungapped: bool,
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

/// The diagonal array of `prelim_subject::prelim_chunk_hits` (a query block of at most 8000
/// letters): `hit_level_array`, `hit_len_array`, `diag_table_length` (NCBI's
/// `diag_array_length`), `diag_table_mask` (`diag_mask`) and `diag_array_size`, the names
/// there, for the two-hit stage of LOSAT_X_PAIRPARDISCARR (child module `disc`). With the
/// diagonal hash the arrays are empty and nothing reads them.
pub(super) struct DiagArray<'a> {
    pub(super) levels: &'a mut Vec<DiagStruct>,
    pub(super) lens: &'a mut Vec<u8>,
    pub(super) length: isize,
    pub(super) mask: isize,
    pub(super) size: usize,
}

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:135-149
/// A pair after the word extension of `s_BlastNaExtend` (`q_offset` and `s_offset` each moved left
/// by `ext_left`), with its context and the `s_range` of its scan range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Seed {
    context: u32,
    q: u32,
    s: u32,
    s_range: u32,
}

/// No NCBI counterpart (scheduling only): where `seed_of_pair` puts a strip's seeds, in scan
/// order: a list (`Vec<Seed>`), or EXPERIMENT (LOSAT_X_DCSEEDLEAN, LOSAT_X_BNSEEDLISTS) the
/// lists of the strip's partitions (`part::PartSink`).
trait SeedOut {
    fn push(&mut self, seed: Seed);

    /// EXPERIMENT (LOSAT_X_SEEDPIECE): how many more seeds `push` takes (a seed piece's room,
    /// `x_pair_seed_piece.rs`); unbounded for a strip's list.
    #[inline(always)]
    fn room(&self) -> usize {
        usize::MAX
    }

    /// EXPERIMENT (LOSAT_X_SEEDPIECE): no room any more (a piece whose scan stopped, so that the
    /// words the scan still hands out go to `Cut::run_cold`, which leaves at once).
    #[inline(always)]
    fn close(&mut self) {}
}

/// EXPERIMENT (LOSAT_X_SEEDPIECE): where the scan of a strip goes on after a seed piece: the
/// subject offset of the word (`kmer_start`) of the first pair the piece had no room for, and
/// where that pair is in the word: for a word whose pairs are a run (a CSR cell, the direct
/// table's list), how many pairs of the run the scan read before it (every pair the reference
/// reads, whether `seed_of_pair` keeps it or not); for a word whose pairs are chains, the pair's
/// element of its chain (`node`, its `q_off`) and whether that is the second template's chain. A
/// scan started again from that word's offset gives the same words from that one on (NCBI's scan
/// resumes at `scan_range[0]`), and its first word with pairs is that word.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Resume {
    pub(super) s: usize,
    pub(super) skip: u32,
    pub(super) node: u32,
    pub(super) second: bool,
}

impl Resume {
    /// A strip's scan from its start.
    pub(super) fn at(s: usize) -> Self {
        Resume {
            s,
            ..Resume::default()
        }
    }
}

/// EXPERIMENT (LOSAT_X_SEEDPIECE): the cut of one piece's scan: where the first word goes on
/// (`Resume`: the pairs of its run to pass over, or the element of its chains to start from) and
/// where the scan stopped (`None` while it goes on). Cells, so that the closures of a scan (the
/// word's, the ring's `emit`) share it.
#[derive(Default)]
struct Cut {
    skip: std::cell::Cell<u32>,
    node: std::cell::Cell<u32>,
    second: bool,
    stop: std::cell::Cell<Option<Resume>>,
}

impl Cut {
    fn new(from: &Resume) -> Self {
        Cut {
            skip: std::cell::Cell::new(from.skip),
            node: std::cell::Cell::new(from.node),
            second: from.second,
            stop: std::cell::Cell::new(None),
        }
    }

    /// Whether the scan goes on (no stop yet).
    #[inline(always)]
    fn go(&self) -> bool {
        self.stop.get().is_none()
    }

    /// The pairs `run` of a word, after those the pieces before took (the first word passes
    /// over them): `Ok` with all of them when `out` has room for all (the caller hands each to
    /// its pair in its own loop, the loop of the strips' scan), else `Err` with the pairs passed
    /// over and the rest, for `run_cold`.
    ///
    /// (The caller's loop is the strips' loop, whose body runs for every pair, so that the
    /// loads of the inputs `seed_of_pair` reads stay out of the loop: a loop that tests the room
    /// before each pair, or that calls its pair from a second loop too, kept them in or kept the
    /// pair a call, about 30 instructions more per seed; results/s14p and s14q of line EB.)
    #[inline(always)]
    fn fit<'r, O: SeedOut>(
        &self,
        run: &'r [u32],
        out: &O,
    ) -> Result<&'r [u32], (usize, &'r [u32])> {
        // (a skip only at the word a piece starts at, so the store is out of the loop's way)
        let skip = self.skip.get() as usize;
        let run = if skip == 0 {
            run
        } else {
            self.skip.set(0);
            &run[skip..]
        };
        if run.len() <= out.room() {
            Ok(run)
        } else {
            Err((skip, run))
        }
    }

    /// The pairs `run` (the rest of the run of the word at subject offset `s` after `skip`
    /// pairs) through `pair`, while `out` has room: those before the first pair it has no room
    /// for, where the cut stops (and `out` is closed, so the later words leave at the test of
    /// `go`). `pair` pushes at most one seed.
    #[cold]
    #[inline(never)]
    fn run_cold<O: SeedOut>(
        &self,
        s: usize,
        skip: usize,
        run: &[u32],
        out: &mut O,
        mut pair: impl FnMut(u32, &mut O),
    ) {
        if !self.go() {
            return;
        }
        for (j, &q) in run.iter().enumerate() {
            if out.room() == 0 {
                self.stop.set(Some(Resume {
                    skip: (skip + j) as u32,
                    ..Resume::at(s)
                }));
                out.close();
                return;
            }
            pair(q, out);
        }
    }

    /// For a word whose pairs are chains: where the first word's walk of its chain starts (the
    /// element the piece before stopped at, or 0: the head), and whether its first template's
    /// chain is passed over (the piece before stopped in the second's); 0 and `false` for every
    /// later word.
    #[inline(always)]
    fn chain_start(&self) -> (u32, bool) {
        let node = self.node.replace(0);
        (node, node != 0 && self.second)
    }

    /// For a word whose pairs are chains: whether the pair at element `q_off_1` of the word at
    /// `s` (of the second template's chain where `second`) goes to `seed_of_pair`. The first pair
    /// `out` has no room for stops the cut there.
    #[inline(always)]
    fn take<O: SeedOut>(&self, s: usize, q_off_1: u32, second: bool, out: &mut O) -> bool {
        if !self.go() {
            return false;
        }
        if out.room() == 0 {
            self.stop.set(Some(Resume {
                node: q_off_1,
                second,
                ..Resume::at(s)
            }));
            out.close();
            return false;
        }
        true
    }
}

impl SeedOut for Vec<Seed> {
    #[inline(always)]
    fn push(&mut self, seed: Seed) {
        Vec::push(self, seed);
    }
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

/// No NCBI counterpart (scheduling only): the part of the scope of `run` (`disc`: of
/// `disc::in_scope`) that does not depend on the chunk: the diagonal hash, one hit for a
/// contiguous table and two hits for a discontiguous one, `scan_range == 0`, no small-query word
/// extension, no diagnostics.
fn chunk_free_scope(
    disc: bool,
    use_array_indexing: bool,
    window_size: usize,
    scan_range: usize,
    small_na_word: bool,
    diagnostics: bool,
) -> bool {
    (if disc {
        window_size > 0
    } else {
        window_size == 0
    }) && !use_array_indexing
        && scan_range == 0
        && !small_na_word
        && !diagnostics
}

/// No NCBI counterpart (scheduling only): the test of `run` for a chunk with a contiguous table
/// (the names of `Inputs`): `chunk_free_scope` and a chunk the subject's owner searches by
/// itself. `run` reads it; `wants` reads its parts (the batch's `owner_scope` and the chunk's
/// `one_subject_owner`).
#[inline(always)]
fn owner_chunk(
    use_array_indexing: bool,
    window_size: usize,
    scan_range: usize,
    small_na_word: bool,
    diagnostics: bool,
    one_subject_owner: bool,
) -> bool {
    chunk_free_scope(
        false,
        use_array_indexing,
        window_size,
        scan_range,
        small_na_word,
        diagnostics,
    ) && one_subject_owner
}

/// No NCBI counterpart (scheduling only): the part of `run`'s test of a chunk with a contiguous
/// table that is the same for every chunk and thread of a batch (`chunk_free_scope`, the names of
/// `Inputs`), taken once per batch (`prelim_subject::ScanBatch::new`, its `x_pair_scope`) for
/// `wants`.
pub(super) fn owner_scope(
    use_array_indexing: bool,
    window_size: usize,
    scan_range: usize,
    small_na_word: bool,
    diagnostics: bool,
) -> bool {
    chunk_free_scope(
        false,
        use_array_indexing,
        window_size,
        scan_range,
        small_na_word,
        diagnostics,
    )
}

/// No NCBI counterpart (scheduling only): whether `prelim_subject::prelim_chunk_hits` calls `run`
/// for a chunk (`disc`: the table is discontiguous), asked before it builds `Inputs`: always for a
/// discontiguous table (`disc::run` decides with its own switches), else `run`'s tests for a
/// contiguous table (`owner_chunk`): the batch's part (`scope`, `owner_scope`), a chunk the
/// subject's owner searches by itself (`one_subject_owner`) and a pool of two threads or more, so a
/// chunk that `run` would send back at once (every chunk of a database search: one thread has no
/// pool, and a pool searches the subjects in parallel) costs neither the inputs nor the call. Such
/// a chunk records the reference side of LOSAT_X_PAIRPAR here, as `run` recorded it, while the run
/// records the sides (`recording`, the batch's `xreach::recording`).
#[inline(always)]
pub(super) fn wants(disc: bool, scope: bool, one_subject_owner: bool, recording: bool) -> bool {
    if disc {
        return true;
    }
    let wanted = scope && one_subject_owner && pool_threads().is_some();
    if !wanted && recording {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", false);
    }
    wanted
}

/// No NCBI counterpart (scheduling only): whether the strips take a chunk that the owner of the
/// pool this thread runs in (two threads or more) searches by itself, in a search with these
/// inputs (the names of `Inputs`; `two_stage`: the search's two-stage table, without which the
/// stage never runs): `chunk_free_scope` for a contiguous table, with LOSAT_X_PAIRPARDISC on for
/// a discontiguous table, or for a discontiguous
/// table with the diagonal array the scope of `disc::in_scope_array` with LOSAT_X_PAIRPARDISCARR
/// on. LOSAT_X_BNFEWSUBJ (`x_few_subj.rs`) reads it before it gives a search the shape of a
/// one-subject search. (A chunk without a pool is not asked here: for a discontiguous table
/// LOSAT_X_BNDCONE decides in `disc::run` and `disc::run_array`.)
pub(super) fn takes_owner_chunks(
    two_stage: Option<&TwoStageLookup>,
    use_array_indexing: bool,
    window_size: usize,
    scan_range: usize,
    small_na_word: bool,
    diagnostics: bool,
) -> bool {
    let Some(two_stage) = two_stage else {
        return false;
    };
    let disc = two_stage.disc().is_some();
    // LOSAT_X_PAIRPARDISCARR: the array's stage has the hash stage's scope with the array.
    let disc_array = disc && use_array_indexing;
    let on = if disc_array {
        disc::arr_on()
    } else if disc {
        disc::on()
    } else {
        true
    };
    on && chunk_free_scope(
        disc,
        use_array_indexing && !disc_array,
        window_size,
        scan_range,
        small_na_word,
        diagnostics,
    ) && pool_threads().is_some()
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// The scan ranges `scan_subject_kmers_with_ranges` visits (from `first_scan_range`: the whole
/// subject when it is unmasked, else its `seq_ranges` one by one), each cut into strips of at
/// most `per` grid points; the last strip of a range ends at the range's end. `s_range` is the
/// range's `scan_range[2] + lut_word_length`.
fn plan_strips(
    seq_ranges: &[(i32, i32)],
    s_len: usize,
    word_length: usize,
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
    target_strips: usize,
) -> Vec<StripJob> {
    if seq_ranges.is_empty() || scan_step == 0 {
        return Vec::new();
    }
    let first = first_scan_range(
        seq_ranges,
        s_len,
        word_length,
        lut_word_length,
        subject_masked,
    );
    // the reference's loop over the ranges, walked again for each pass (no list of them)
    let ranges = |f: &mut dyn FnMut(usize, usize)| {
        scan_subject_kmers_with_offsets(seq_ranges, word_length, lut_word_length, first, |s, e| {
            f(s, e)
        })
    };
    let grid_points = |start: usize, end: usize| {
        if start > end {
            0
        } else {
            (end - start) / scan_step + 1
        }
    };
    let mut total = 0usize;
    ranges(&mut |start, end| total += grid_points(start, end));
    const MIN_STRIP: usize = 512;
    let per = total.div_ceil(target_strips.max(1)).max(MIN_STRIP);
    // a range without grid points still has its strip (the kernel is called and returns at
    // once, start > end)
    let mut count = 0usize;
    ranges(&mut |start, end| count += grid_points(start, end).div_ceil(per).max(1));
    let mut jobs = Vec::with_capacity(count);
    ranges(&mut |start, end| {
        let s_range = end.saturating_add(lut_word_length);
        let n = grid_points(start, end);
        if n == 0 {
            jobs.push(StripJob {
                start,
                end,
                s_range,
            });
            return;
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
    });
    debug_assert_eq!(jobs.len(), count);
    jobs
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1504-1533
/// Phase 1 for one strip: the scan kernel `scan_subject_kmers_with_ranges` picks for this
/// lookup width and step (`select_mb_scan_kind`: NCBI's Any scanner for a masked subject), from
/// `job.start` to `job.end`; every lookup chain hit as in `on_word` (the chain order); every pair through `seed_of_pair`.
fn scan_strip<O: SeedOut>(x: &Inputs<'_>, job: &StripJob, out: &mut O) {
    scan_strip_from::<O, false>(x, job, Resume::at(job.start), out);
}

/// `scan_strip` from the word at `from.s` (a word of the strip's scan: `job.start`, or where a
/// piece stopped), that word from where `from` says (`Resume`). EXPERIMENT (LOSAT_X_SEEDPIECE):
/// with `PIECE` the scan stops at the first pair `out` has no room for (`Cut`) and returns where
/// it goes on; `None` when the strip is done.
fn scan_strip_from<O: SeedOut, const PIECE: bool>(
    x: &Inputs<'_>,
    job: &StripJob,
    from: Resume,
    out: &mut O,
) -> Option<Resume> {
    // EXPERIMENT (LOSAT_X_BNCSR): a contiguous table laid out as CSR has no chains.
    // EXPERIMENT (LOSAT_X_BNCSRPACK): the form of its head chosen once, each form's scan its own.
    if let Some(csr) = x.two_stage.x_mb_csr() {
        // EXPERIMENT (LOSAT_X_BNPVONCE): one presence test per word on a head that fits the cache.
        return match csr.x_packed() {
            Some(packed) if super::x_bn_pv_once(&packed) => {
                scan_strip_csr::<O, _, false, PIECE>(x, &packed, job, from, out)
            }
            Some(packed) => scan_strip_csr::<O, _, true, PIECE>(x, &packed, job, from, out),
            None => scan_strip_csr::<O, _, true, PIECE>(x, csr, job, from, out),
        };
    }
    super::x_bn_rank_once(false);
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let packed = x.search_seq_packed;
    let subject_len = x.s_len;
    let scan_step = x.scan_step;
    let s_range = job.s_range;
    let end = job.end;
    let cut = Cut::new(&from);
    let mut on_kmer = |kmer_start: usize, current_lut_kmer: u64| {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
        let mut pair = |q_off_1: u32| {
            if q_off_1 != 0 {
                if PIECE && !cut.take(kmer_start, q_off_1, false, out) {
                    return;
                }
                seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
            }
        };
        if PIECE {
            let (node, _) = cut.chain_start();
            two_stage.for_each_hit_from(current_lut_kmer, node, &mut pair);
            cut.go()
        } else {
            two_stage.for_each_hit(current_lut_kmer, &mut pair);
            true
        }
    };
    scan_kernel(
        packed,
        subject_len,
        lut_word_length,
        scan_step,
        x.subject_masked,
        from.s,
        end,
        &mut on_kmer,
    );
    cut.stop.get()
}

/// EXPERIMENT (LOSAT_X_BNCSR): `scan_strip` over a contiguous table laid out as CSR
/// (`lookup::x_bn_csr`): the words with hits (the presence test of `for_each_hit` first), in
/// order, each `RING` words late through `x_dc_csr::Ring` (whose prefetches bring in the cell's
/// offsets and the ends of its run), then the cell's run, which holds the chain's offsets in the
/// chain's order; so the same pairs reach `seed_of_pair` in the same order. EXPERIMENT
/// (LOSAT_X_BNCSRPACK): over a head of either form (`x_dc_csr::XHead`).
#[inline(never)]
fn scan_strip_csr<
    O: SeedOut,
    H: crate::algorithm::blastn::lookup::x_dc_csr::XHead,
    const PV: bool,
    const PIECE: bool,
>(
    x: &Inputs<'_>,
    csr: &H,
    job: &StripJob,
    from: Resume,
    out: &mut O,
) -> Option<Resume> {
    use crate::algorithm::blastn::lookup::x_dc_csr::XRing;
    let two_stage = x.two_stage;
    let s_range = job.s_range;
    let cut = Cut::new(&from);
    // EXPERIMENT (LOSAT_X_BNRANKONCE): with `rank_once` a word's run is the one of the head's
    // element the ring hands over with it (a packed head's word is ranked once), else the cell's
    // run read again.
    let rank_once = super::x_bn_rank_once(true);
    let mut emit = |(kmer_start, index): (usize, u64), e: usize| {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
        let mut run = if rank_once {
            csr.run_at(e)
        } else {
            csr.run(index as usize)
        };
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
        for &q_off_1 in run {
            seed_of_pair(x, q_off_1 as usize - 1, kmer_start, s_range, out);
        }
    };
    let mut ring: H::Ring<'_, (usize, u64)> = csr.ring();
    // (a cut stops only in `emit`, so only a word that let one out can have stopped it)
    let mut on_kmer = |kmer_start: usize, current_lut_kmer: u64| {
        if !PV || two_stage.has_hits(current_lut_kmer) {
            let left = ring.push((kmer_start, current_lut_kmer), |w| w.1 as usize, &mut emit);
            !PIECE || !left || cut.go()
        } else {
            true
        }
    };
    scan_kernel(
        x.search_seq_packed,
        x.s_len,
        two_stage.lut_word_length(),
        x.scan_step,
        x.subject_masked,
        from.s,
        job.end,
        &mut on_kmer,
    );
    // (a stopped scan drops the words still in the ring: the next piece scans them again)
    if !PIECE || cut.go() {
        ring.flush(|w| w.1 as usize, &mut emit);
    }
    cut.stop.get()
}

/// The scan kernel `scan_subject_kmers_with_ranges` picks for this lookup width and step,
/// from `start` to `end` (`scan_strip`, and the child module `na`).
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn scan_kernel<F: FnMut(usize, u64) -> bool>(
    packed: &[u8],
    subject_len: usize,
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) {
    match select_mb_scan_kind(lut_word_length, scan_step, subject_masked) {
        Some(MbScanSubjectKind::Any) => scan_subject_kmers_range_mb_any(
            packed,
            subject_len,
            lut_word_length,
            scan_step,
            subject_masked,
            start,
            end,
            on_kmer,
        ),
        Some(MbScanSubjectKind::Scan9_1) => {
            scan_subject_kmers_range_mb_9_1(packed, subject_len, scan_step, start, end, on_kmer)
        }
        Some(MbScanSubjectKind::Scan9_2) => {
            scan_subject_kmers_range_mb_9_2(packed, subject_len, scan_step, start, end, on_kmer)
        }
        Some(MbScanSubjectKind::Scan10_1) => {
            scan_subject_kmers_range_mb_10_1(packed, subject_len, scan_step, start, end, on_kmer)
        }
        Some(MbScanSubjectKind::Scan10_2) => {
            scan_subject_kmers_range_mb_10_2(packed, subject_len, scan_step, start, end, on_kmer)
        }
        Some(MbScanSubjectKind::Scan10_3) => {
            scan_subject_kmers_range_mb_10_3(packed, subject_len, scan_step, start, end, on_kmer)
        }
        Some(MbScanSubjectKind::Scan11_1Mod4) => scan_subject_kmers_range_mb_11_1mod4(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            on_kmer,
        ),
        Some(MbScanSubjectKind::Scan11_2Mod4) => scan_subject_kmers_range_mb_11_2mod4(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            on_kmer,
        ),
        Some(MbScanSubjectKind::Scan11_3Mod4) => scan_subject_kmers_range_mb_11_3mod4(
            packed,
            subject_len,
            scan_step,
            start,
            end,
            on_kmer,
        ),
        None => scan_subject_kmers_range(
            packed,
            subject_len,
            lut_word_length,
            scan_step,
            subject_masked,
            start,
            end,
            on_kmer,
        ),
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1106-1144
/// Phase 1 for one pair: the statements of `pair_extend::extend_offset_pairs` before the
/// diagonal is read (context, bounds test, word extension), which read only the inputs. A pair that
/// the reference drops with `continue` here is dropped; any other becomes a `Seed`, pushed to
/// `out` (a strip's list, or with LOSAT_X_BNSEEDLISTS / LOSAT_X_DCSEEDLEAN its partition's list,
/// `part::PartSink`).
#[inline(always)]
fn seed_of_pair<O: SeedOut>(
    x: &Inputs<'_>,
    q_off0: usize,
    kmer_start: usize,
    s_range: usize,
    out: &mut O,
) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let word_length = two_stage.word_length();
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let encoded_query_concat_blastna_with_sentinels = x.encoded_query_concat_blastna_with_sentinels;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:730-733
    let context_idx = x.query_context_index.context_for_offset(q_off0);
    let ctx = &x.query_contexts[context_idx];
    let q_context_start = ctx.query_offset as usize;
    let q_context_end = q_context_start + ctx.len;

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
            let s_base = (s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6;
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
                let s_base = (s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6;
                // SAFETY: q_right is bounded by the outer-sentinel query buffer.
                let q_base =
                    unsafe { *encoded_query_concat_blastna_with_sentinels.get_unchecked(q_right) };
                if s_base != q_base {
                    break;
                }
                ext_right += 1;
                q_right += 1;
                s_off += 1;
                if s_off.is_multiple_of(COMPRESSION_RATIO) {
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
/// Phase 2 (owner only): the rest of `pair_extend::extend_offset_pairs` for the seeds of one
/// strip, in order, statement for statement for the diagonal hash in one-hit mode (`two_hits`
/// false, so `off_found` stays false and `hit_ready` true until the cutoff test); the debug, trace
/// and diagonal-array branches are outside the scope and left out.
///
/// With `PT` (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`) the seeds are one partition's seeds
/// of one strip, `diag_hash` and `ungapped_hits` are the partition's, and `out` gets the key
/// (`out.strip`, `pos[k]`) of every pushed hit and of every insert that appends a cell, the
/// buckets the inserts change, and the time of the ungapped extensions (instead of the shared
/// LOSAT_TIMING counters); without it `out` and `pos` are not read.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn fold<const PT: bool>(
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
    fold_t::<PT, true>(
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

/// `fold`; EXPERIMENT (LOSAT_X_BNSEEDLISTS): without `TRACK` the inserts of a partition do not
/// record the buckets they change (`part::insert_t`; the merge then writes back every bucket of
/// the partition).
#[allow(clippy::too_many_arguments)]
fn fold_t<const PT: bool, const TRACK: bool>(
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
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let diag_offset = x.diag_offset;
    let diag_hash_window = x.diag_hash_window;
    let timing_enabled = timing.is_some();
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the table has level bounds; the bound's
    // constants of these seeds (one hit: `window_size` is 0 in the scope).
    let x_lv = diag_hash.x_lv_on();
    let x_lv_seeds = super::x_diag_lv::LvSeeds::new(x.window_size > 0, word_length, x.window_size);
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
        let (last_hit, _hit_saved) = {
            // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): the same decisions without the walk
            // of a bucket below the bound. (A one-hit table keeps no memo, LOSAT_X_DIAGLV2.)
            let (level, _hit_len, hit_saved) = if x_lv {
                diag_hash.retrieve_below::<false>(
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
        let off_found = false;
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:887-895
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
            false, // skip double-word check (not in two-hit block)
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
        // (`hit_ready` is still 1 in one-hit mode; the reference keeps this test.)
        if !hit_ready {
            part::insert_t::<PT, false, false, TRACK>(
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
                timing
                    .ungapped_ns
                    .fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
                timing
                    .ungapped_calls
                    .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
        }
        let ungapped_se = ungapped.s_start + ungapped.length;
        let ungapped_score = ungapped.score;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:928-929
        let ungapped_s_end_pos = ungapped_se + diag_offset as usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:923-932
        if !(off_found || ungapped_score >= cutoff_score) {
            hit_ready = false;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
            part::insert_t::<PT, false, false, TRACK>(
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
        super::x_push_ungapped(
            ungapped_hits,
            UngappedHit::saved(q_idx, ctx, q_off, s_off, &ungapped),
        );
        if PT {
            out.hit_keys.push(x_key);
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:939-941
        part::insert_t::<PT, false, false, TRACK>(
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

/// No NCBI counterpart: the strips of one scan and their results (scheduling only). Strip
/// `i` is scanned once, by whichever thread claims it; the owner takes the results in
/// strip order. `taken` counts the strips the owner has taken; with `ahead > 0` a helper
/// claims a strip only while fewer than `ahead` strips past `taken` are claimed (a bound on
/// the results waiting for the owner, up to one strip per helper over it), and the buffers
/// the owner has consumed go back to `spare` for the next strips. Only a pool's stage has
/// them: on one thread `run_strips` scans and consumes each strip in turn into one buffer.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Strips<T> {
    cursor: AtomicUsize,
    taken: AtomicUsize,
    ahead: usize,
    spare: Mutex<Vec<Vec<T>>>,
    slots: Vec<Mutex<Option<Vec<T>>>>,
    ready: Vec<AtomicBool>,
    failed: AtomicBool,
}

// No NCBI counterpart: marks the strips as failed if a scan unwinds, so that the owner stops
// waiting for its result.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Producing<'a, T> {
    strips: &'a Strips<T>,
    done: bool,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<T> Drop for Producing<'_, T> {
    fn drop(&mut self) {
        if !self.done {
            self.strips.failed.store(true, Ordering::Release);
        }
    }
}

// No NCBI counterpart: when the owner leaves its loop (at the end, or unwinding), the helpers
// that wait for it to take a strip (`ahead > 0`) stop waiting.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Closing<'a, T>(&'a Strips<T>);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<T> Drop for Closing<'_, T> {
    fn drop(&mut self) {
        self.0.failed.store(true, Ordering::Release);
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<T> Strips<T> {
    fn new(n: usize, ahead: usize) -> Self {
        Self {
            cursor: AtomicUsize::new(0),
            taken: AtomicUsize::new(0),
            ahead,
            spare: Mutex::new(Vec::new()),
            slots: (0..n).map(|_| Mutex::new(None)).collect(),
            ready: (0..n).map(|_| AtomicBool::new(false)).collect(),
            failed: AtomicBool::new(false),
        }
    }

    /// A consumed buffer, kept for a later strip when the look-ahead is bounded.
    fn recycle(&self, mut out: Vec<T>) {
        if self.ahead > 0 {
            out.clear();
            self.spare
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .push(out);
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
        // The spare buffers exist only while strips run ahead (`recycle` fills them then); the
        // one-hit path (`ahead == 0`) takes no lock here.
        let mut out = if self.ahead > 0 {
            self.spare
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .pop()
                .unwrap_or_default()
        } else {
            Vec::new()
        };
        produce(i, &mut out);
        *self.slots[i].lock().unwrap_or_else(|e| e.into_inner()) = Some(out);
        self.ready[i].store(true, Ordering::Release);
        guard.done = true;
    }

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn helper(&self, produce: &impl Fn(usize, &mut Vec<T>)) {
        loop {
            if self.ahead > 0 {
                // spin a little, then sleep 20 us, doubling up to 1 ms (a helper that the
                // owner does not need costs a few wake-ups per millisecond); on wasm32 only
                // yield: `std::thread::sleep` panics under the browser's WASI shim (as in
                // `x_pair_seed_part.rs`'s `pause`)
                let mut idle = 0u32;
                #[cfg(not(target_arch = "wasm32"))]
                let mut nap = 20u64;
                loop {
                    let cursor = self.cursor.load(Ordering::Relaxed);
                    if cursor >= self.slots.len() || self.failed.load(Ordering::Acquire) {
                        return;
                    }
                    if cursor < self.taken.load(Ordering::Acquire) + self.ahead {
                        break;
                    }
                    idle = idle.saturating_add(1);
                    if idle < 128 {
                        std::hint::spin_loop();
                    } else {
                        #[cfg(not(target_arch = "wasm32"))]
                        {
                            std::thread::sleep(std::time::Duration::from_micros(nap));
                            nap = (nap * 2).min(1000);
                        }
                        #[cfg(target_arch = "wasm32")]
                        std::thread::yield_now();
                    }
                }
            }
            let Some(i) = self.claim() else {
                return;
            };
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
                let out = self.slots[i]
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .take()
                    .unwrap_or_default();
                self.taken.store(i + 1, Ordering::Release);
                return out;
            }
            assert!(
                !self.failed.load(Ordering::Acquire),
                "PAIRPAR / LOSAT_X_PAIRPARDISC: a strip scan failed"
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
/// `0, 1, ..., n - 1` in that order, each as soon as its result is there. `ahead` (0: no
/// bound) bounds how far the helpers scan past the strip being consumed. Without helpers (or
/// outside a pool) this thread scans and consumes the strips in turn through one buffer, with
/// no per-strip slots: a masked subject's chunk has a strip for every scan range.
#[allow(clippy::too_many_arguments)]
fn run_strips<T: Send>(
    n: usize,
    helpers: usize,
    ahead: usize,
    produce: &(impl Fn(usize, &mut Vec<T>) + Sync),
    mut consume: impl FnMut(usize, &[T]),
    scan_ns: &mut u64,
    wait_ns: &mut u64,
    timed: bool,
) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if helpers > 0 && rayon::current_thread_index().is_some() {
        let strips = Strips::new(n, ahead);
        let strips = &strips;
        rayon::in_place_scope(|scope| {
            for _ in 0..helpers.min(n) {
                scope.spawn(move |_| strips.helper(produce));
            }
            let _closing = Closing(strips);
            for i in 0..n {
                let out = strips.take(i, produce, scan_ns, wait_ns, timed);
                consume(i, &out);
                strips.recycle(out);
            }
        });
        return;
    }
    let _ = (helpers, ahead, &wait_ns);
    let mut out = Vec::new();
    for i in 0..n {
        out.clear();
        let t0 = timed.then(std::time::Instant::now);
        produce(i, &mut out);
        if let Some(t0) = t0 {
            *scan_ns += t0.elapsed().as_nanos() as u64;
        }
        consume(i, &out);
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

impl StageStats {
    /// No NCBI counterpart: a partitioned stage's counters (LOSAT_X_PAIRPARFOLD); its
    /// partitions' time in ungapped extensions goes to the LOSAT_TIMING counters here, once.
    fn of_part(strips: usize, r: part::PartStats, timing: Option<&BlastnTiming>) -> Self {
        if let Some(timing) = timing {
            timing
                .ungapped_ns
                .fetch_add(r.ungapped_ns, Ordering::Relaxed);
            timing
                .ungapped_calls
                .fetch_add(r.ungapped_calls, Ordering::Relaxed);
        }
        StageStats {
            strips,
            counts: r.counts,
            scan_pause_ns: r.scan_pause_ns,
            fold_ns: r.fold_ns,
            owner_scan_ns: r.owner_scan_ns,
            wait_ns: r.wait_ns,
        }
    }
}

/// No NCBI counterpart (scheduling only): strips per thread of a chunk's scan with the
/// partitioned fold (LOSAT_X_PAIRPARFOLD), and strips per thread that the participants may scan
/// past the slowest partition (`seed_stage`).
const PART_STRIPS_PER_THREAD: usize = 64;
const PART_AHEAD_PER_THREAD: usize = 2;

/// No NCBI counterpart: reads LOSAT_X_BNSEEDLISTS once (a strict-set switch: on unless it is
/// `0`). With it the partitioned stages of this module and of `na` (LOSAT_X_PAIRPARFOLD) take
/// LOSAT_X_DCSEEDLEAN's form (`part::seed_stage_part_lists`): the home loop scans each strip
/// into its partitions' lists, and the folds do not record the buckets they change, the merge
/// writing back every bucket of each partition. Each partition folds the same seeds in the same
/// order, so the hits and the table are those of the cut and the tracked folds.
pub(super) fn seed_lists_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNSEEDLISTS"))
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// The seed stage of one chunk on `threads` threads (this one and `threads - 1` helpers of
/// its pool): strips scanned by any thread, folded here in strip order into `hash` and
/// `hits`. With `fold_part` (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`) the fold is split by
/// diagonal-hash bucket over the threads instead, in the form LOSAT_X_BNSEEDLISTS says.
fn seed_stage(
    x: &Inputs<'_>,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
    fold_part: Option<part::FoldPart<'_>>,
) -> StageStats {
    seed_stage_of(
        x,
        threads,
        hash,
        hits,
        timing,
        timed,
        fold_part,
        seed_lists_on(),
    )
}

/// `seed_stage`; `lists`: LOSAT_X_BNSEEDLISTS's partitioned stage (`seed_lists_on`).
#[allow(clippy::too_many_arguments)]
fn seed_stage_of(
    x: &Inputs<'_>,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    timed: bool,
    fold_part: Option<part::FoldPart<'_>>,
    lists: bool,
) -> StageStats {
    // The partitioned fold holds the seeds of the strips within its look-ahead of the slowest
    // partition: `PART_AHEAD_PER_THREAD` of `PART_STRIPS_PER_THREAD` strips per thread, 1/32 of a
    // chunk's seeds (as `disc`), where 2 of 8 per thread held a quarter (blastn of the Sakai
    // genome against MG1655, 30 M seeds, at 4 threads: about +70 MiB over the search without
    // strips; LOSAT_X_BNFEWSUBJ gives that shape to searches of two subjects or more).
    let strips_per_thread = if fold_part.is_some() {
        PART_STRIPS_PER_THREAD
    } else {
        8
    };
    let jobs = plan_strips(
        x.subject_seq_ranges,
        x.s_len,
        x.two_stage.word_length(),
        x.two_stage.lut_word_length(),
        x.scan_step,
        x.subject_masked,
        strips_per_thread * threads.max(4),
    );
    let produce = |i: usize, out: &mut Vec<Seed>| scan_strip(x, &jobs[i], out);
    if let Some(fold_part) = fold_part {
        let parts = fold_part.parts;
        // EXPERIMENT (LOSAT_X_BNSEEDLISTS): every partitioned stage takes the partitions' lists and
        // the folds that do not record their buckets.
        if lists {
            crate::utils::xreach::x_reach!("LOSAT_X_BNSEEDLISTS", true);
            let produce_lists =
                |i: usize, out: &mut part::PartSink<'_>| scan_strip(x, &jobs[i], out);
            // EXPERIMENT (LOSAT_X_SEEDPIECE): a piece of strip i, from its start or a stop.
            let produce_piece = |i: usize, from: Option<Resume>, out: &mut part::PieceSink<'_>| {
                let job = &jobs[i];
                scan_strip_from::<_, true>(x, job, from.unwrap_or(Resume::at(job.start)), out)
            };
            let fold_one = |st: &mut part::PartState, seeds: &[Seed], pos: &[u32]| {
                let part::PartState {
                    hash,
                    hits,
                    out,
                    counts,
                    scan_pause_ns,
                    ..
                } = st;
                fold_t::<true, false>(
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
            let r = part::seed_stage_part_lists(
                jobs.len(),
                threads,
                PART_AHEAD_PER_THREAD * threads,
                hash,
                hits,
                fold_part,
                &produce,
                &produce_lists,
                &produce_piece,
                &fold_one,
                timed,
            );
            part::record(part::ONE_HIT, parts, &r);
            return StageStats::of_part(jobs.len(), r, timing);
        }
        let fold_one = |st: &mut part::PartState, seeds: &[Seed], pos: &[u32]| {
            let part::PartState {
                hash,
                hits,
                out,
                counts,
                scan_pause_ns,
                ..
            } = st;
            fold::<true>(
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
            PART_AHEAD_PER_THREAD * threads,
            hash,
            hits,
            fold_part,
            &produce,
            &fold_one,
            timed,
        );
        part::record(part::ONE_HIT, parts, &r);
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
        0,
        &produce,
        |_, seeds| {
            let t0 = timed.then(std::time::Instant::now);
            fold::<false>(
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

// No NCBI counterpart: a copy of the diagonal hash, for the unit tests.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[cfg(test)]
fn copy_hash(h: &DiagHashTable) -> DiagHashTable {
    DiagHashTable {
        num_buckets: h.num_buckets,
        occupancy: h.occupancy,
        capacity: h.capacity,
        backbone: h.backbone.clone(),
        chain: h.chain.clone(),
        offset: h.offset,
        window: h.window,
        max_level: h.max_level.clone(),
        x_last: h.x_last.clone(),
        #[cfg(debug_assertions)]
        x_last_fresh: h.x_last_fresh,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
/// Dispatch of `prelim_subject::prelim_chunk_hits`: the seed stage of this chunk (the loop above
/// with its one-hit extension) when the chunk is in scope (`chunk_free_scope`, a pool of two
/// threads or more, one subject owner); `true` means it has run and the
/// reference scan is skipped. `fold_scratch` holds the allocations of
/// LOSAT_X_PAIRPARFOLD's partition tables and buffers (`x_pair_seed_part.rs`; every
/// partitioned stage copies `diag_hash` into them).
pub(super) fn run(
    x: &Inputs<'_>,
    diag_hash: &mut DiagHashTable,
    diag_array: DiagArray<'_>,
    ungapped_hits: &mut Vec<UngappedHit>,
    fold_scratch: &mut PartScratch,
) -> bool {
    // A discontiguous megablast table: LOSAT_X_PAIRPARDISC (the diagonal hash) or
    // LOSAT_X_PAIRPARDISCARR (the diagonal array) decide.
    if x.two_stage.disc().is_some() {
        return disc::run(x, diag_hash, diag_array, ungapped_hits, fold_scratch);
    }
    if !owner_chunk(
        x.use_array_indexing,
        x.window_size,
        x.scan_range,
        x.small_na_word,
        x.diagnostics,
        x.one_subject_owner,
    ) {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", false);
        return false;
    }
    let Some(threads) = pool_threads() else {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", false);
        return false;
    };
    crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", true);
    let timing = x.timing;
    let timed = timing.is_some() || std::env::var_os("LOSAT_TIMING").is_some();
    let t_start = timed.then(std::time::Instant::now);
    // EXPERIMENT (LOSAT_X_PAIRPARFOLD): the fold split over the pool's threads.
    let fold_part = part::partitions(threads).map(|parts| part::FoldPart {
        scratch: fold_scratch,
        parts,
        latest: None,
    });
    let stats = seed_stage(
        x,
        threads,
        diag_hash,
        ungapped_hits,
        timing,
        timed,
        fold_part,
    );
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

/// No NCBI counterpart: with LOSAT_TIMING the counters of PAIRPAR, on stderr at the end of the
/// search; then those of
/// LOSAT_X_PAIRPARDISC (child module `disc`) and LOSAT_X_PAIRPARNA (child module `na`).
pub(super) fn summary() {
    if std::env::var_os("LOSAT_TIMING").is_some() {
        let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
        let s = |c: &AtomicU64| c.load(Ordering::Relaxed) as f64 / 1e9;
        eprintln!(
            "[TIMING] x_pairpar_seed: chunks={} strips={} seeds={} skipped={} fold={:.3}s owner_scan={:.3}s owner_wait={:.3}s{}",
            load(&X_CHUNKS),
            load(&X_STRIPS),
            load(&X_SEEDS),
            load(&X_SKIPPED),
            s(&X_FOLD_NS),
            s(&X_OWNER_SCAN_NS),
            s(&X_WAIT_NS),
            part::timing_fields(part::ONE_HIT)
        );
    }
    disc::summary();
    na::summary();
}

#[cfg(test)]
#[path = "x_pair_seed_tests.rs"]
mod tests;
