//! EXPERIMENT (LOSAT_X_DISCLUTPAR, shadow LOSAT_X_DISCLUTPARSHADOW): the discontiguous megablast
//! tables of `build_disc_mb_lookup` (`s_FillDiscMBTable`) filled by the threads of the search
//! pool, each thread owning a fixed set of table cells of both templates' tables. The plan of the
//! cells and the pool test are those of LOSAT_X_PAIRPAR's megablast table (`x_mb_lookup_par.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:754-767,775-786,788-804
//! ```c
//!       for (index = from; index <= to; index++) {
//!          val = *++seq;
//!          /* if an ambiguity is encountered, do not add
//!             any words that would contain it */
//!          if ((val & BLAST2NA_MASK) != 0) {
//!             accum = 0;
//!             pos = seq + template_length;
//!             continue;
//!          }
//!
//!          /* get next base */
//!          accum = (accum << BITS_PER_NUC) | val;
//!          if (seq < pos)
//!             continue;
//! ...
//!          ecode1 = ComputeDiscontiguousIndex(accum, template_type);
//!          if (mb_lt->hashtable[ecode1] == 0) {
//! ...
//!             PV_SET(pv_array, ecode1, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode1/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
//!          mb_lt->hashtable[ecode1] = index;
//!
//!          if (!kTwoTemplates)
//!             continue;
//! ...
//!          ecode2 = ComputeDiscontiguousIndex(accum, second_template_type);
//!          if (mb_lt->hashtable2[ecode2] == 0) {
//! ...
//!             PV_SET(pv_array, ecode2, pv_array_bts);
//!          }
//!          else {
//!             helper_array2[ecode2/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos2[index] = mb_lt->hashtable2[ecode2];
//!          mb_lt->hashtable2[ecode2] = index;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:808-813
//! ```c
//!    longest_chain = 2;
//!    for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
//!        longest_chain = MAX(longest_chain, helper_array[index]);
//!     /* +1 because helper_array is not incremented for the first position of a
//!       word */
//!    mb_lt->longest_chain = longest_chain + 1;
//! ```
//!
//! Why the tables are the same. The update of a word for the first template reads and writes
//! only the cell `ecode1` of `hashtable`, the presence-vector word `ecode1 >> pv_array_bts`,
//! the counter `ecode1 / kCompressionFactor` of `helper_array` and `next_pos[index]`; the
//! update for the second template only the cell `ecode2` of `hashtable2`, the presence-vector
//! word `ecode2 >> pv_array_bts` (of the same `pv_array`), the counter of `helper_array2` and
//! `next_pos2[index]`. `index` is the word's own 1-based query offset, unique within a template,
//! and `next_pos[index]` is the cell's head before the word's update. The cells (`hashsize` of
//! them, the same for both templates) are cut into blocks of `2^block_shift` cells with
//! `2^block_shift` a multiple of both `kCompressionFactor` (2048) and `2^pv_array_bts`
//! (`x_mb_lookup_par::plan`, blocks dealt round robin to P owners), so a block owns whole
//! presence-vector words and whole counters of either template.
//!
//! The words are read in the loop's order: the stream of all positions of all unmasked ranges of
//! all queries (`segments`), cut into chunks of `CHUNK_POSITIONS` positions, each chunk done in
//! three phases, each a set of P jobs that returns when all have returned:
//! - A (read): the chunk's positions are cut into P spans in stream order; reader `a` reads the
//!   words that end in its span with the loop's own statements (the same ambiguity resets per
//!   range, `accum` and `index`), entering a range up to `template_length - 1` positions before
//!   its span (never before the range's start). The template indices use only the last
//!   `template_length` bases of `accum` (the low `2 * template_length` bits; checked for all
//!   twelve templates), and a word is read iff the last `template_length` positions are inside the
//!   range and unambiguous, so the words of the span and their `index`, `ecode1` and `ecode2` are
//!   the loop's. It appends `(index, ecode1)` to bucket `(a, owner of ecode1's block)` and
//!   `(index, ecode2)` to the second template's bucket `(a, owner of ecode2's block)`, in reading
//!   order.
//! - B (update): owner `o` takes the buckets `(0, o)`, `(1, o)`, ..., `(P-1, o)`: every word of the
//!   chunk whose cell lies in its blocks, in reading order (the spans follow each other). For each
//!   it makes the loop's update of the cell (the presence bit when the cell is empty, else the
//!   counter; then the cell takes `index`) and keeps the cell's previous head in the entry. The
//!   chunks are done in order, so each cell of `hashtable` (and of `hashtable2`) receives the same
//!   words in the same order as in the loop: the same final head, the same counters, and each
//!   word's kept head is the value the loop writes to `next_pos[index]`. Both templates set bits
//!   in the same `pv_array`, but a presence-vector word lies in the block of the cells whose bits
//!   it holds, whichever template sets them (same `hashsize`, same plan), so it has one writer; a
//!   bit is only ever set, never cleared or read during the fill, so the order of the two
//!   templates' sets does not matter.
//! - C (link): reader `a` writes `next_pos[index]` (and `next_pos2[index]`) of the words it read
//!   from their kept heads. Each index is written once; a span's words have a contiguous run of
//!   query offsets, so the readers write disjoint parts of `next_pos`.
//!
//! `longest_chain` is computed afterwards from the same counters by the original code
//! (`lookup::disc_longest_chain`). Jobs write disjoint elements, so the number of workers changes
//! the running time only. No floating-point operation is involved.
//!
//! Why three phases (scheduling only): in a first version every worker read every word and applied
//! those of its cells, writing `next_pos[index]` itself; each worker then paid the whole read
//! (about 0.09 s for 10 million positions at the workstation's clock), and neighbouring
//! `next_pos` entries were written by different workers (a cache line moving between cores for
//! most writes): at 4 threads that build was no faster than the loop (ws, `logs/disclutpar_report.md`).
//! Here each position is read once, each owner reads only its own words, and `next_pos` is
//! written by spans. The buckets hold 8 bytes per word of one chunk and template (4 MiB per
//! template for `CHUNK_POSITIONS`), reused from chunk to chunk; phase B prefetches the cells 16
//! entries ahead on x86_64.
//!
//! The pool. The switch acts only where a LOSAT search pool of two threads or more is current
//! (`x_mb_lookup_par::pool_threads`: the calling thread is a rayon worker, and its pool has two
//! threads or more). `build_disc_mb_lookup` is called by `coordination::build_lookup_tables`,
//! which `run.rs` `search_query_batch` calls for a query batch (and for each query chunk of a
//! split batch). That function runs inside the search's pool: `utils::threading::with_search_pool`
//! builds a pool of `-num_threads` threads when there are two or more and runs the whole search
//! with `ThreadPool::install` on the pool's slot zero (`-subject` and `-db` alike); a
//! LOSAT_X_QBATCHPAR slot of `m >= 2` threads searches its batch the same way on a private pool
//! of `m` threads (`with_search_pool` again). With `-num_threads 1`, in a QBATCHPAR slot of one
//! thread, in a build without the `parallel` feature or WebAssembly without `wasm-threads`, the
//! calling thread is no rayon worker (or there is no pool) and the original loop runs.
//! `rayon::current_thread_index` reads a thread-local and starts no pool, and
//! `rayon::current_num_threads` is only asked on a worker, where it is the size of that worker's
//! pool, so rayon's global pool is never started or used (the command line never builds it; a
//! program that embeds LOSAT and calls a one-thread search from a thread of its own rayon pool
//! would see that pool used, as `x_mb_lookup_par` and `x_link_sort` say). Each phase is one
//! `rayon::scope` whose jobs use only their own buckets and cells and never wait for anything (no
//! lock, spin or sleep), so no job can block a thread that steals it.
//!
//! Paths that reach it: `blastn -task dc-megablast`, `-subject` and `-db`, with one template
//! (`coding` or `optimal`) or two (`coding_and_optimal`), template lengths 16, 18 and 21 and
//! word sizes 11 and 12, at two threads or more. Megablast and blastn build the contiguous
//! tables (`build_two_stage_lookup`, `build_na_lookup`, `build_pv_direct_lookup`) unless
//! `-template_type` and `-template_length` are given (NCBI applies them to any task,
//! `coordination::determine_template`); then they build this table too. The callers in
//! `x_pair_seed_disc.rs` are unit tests (on the test thread, so the original loop).
//!
//! The shadow implies the switch: after the tables are built, it builds them again the other
//! way (serially with the original loop when the pool built them; with four parts one after
//! another when no pool was there) and compares every element of every array and the longest
//! chain. A difference prints `[X_SHADOW] DISCLUTPAR ... DIFFERENT ...` and panics; at exit
//! `main` prints `[X_SHADOW] DISCLUTPAR tables=<compared> pool=<built on a pool>
//! two_templates=<of them with two templates> cells=<elements compared> diffs=0 identical`.
//! With LOSAT_TIMING (any switch state) `main` prints `[TIMING] disc_lookup_build: <s>s
//! (builds=<n> parallel=<n on a pool>)` when at least one table was built: `<s>` is the sum of
//! the wall times of every build of the process (allocation, fill and longest chain, not the
//! shadow's second build); with LOSAT_X_QBATCHPAR it includes the slots' builds, which run
//! concurrently.

use super::x_mb_lookup_par::{plan, pool_threads, Shared, K_COMPRESSION_FACTOR_SHIFT};
use super::{
    build_unmasked_ranges, disc_longest_chain, fill_disc_mb_tables, DiscTables, DiscTablesRef,
    DiscWords, PvArrayType, BLAST2NA_MASK, PV_ARRAY_MASK,
};
use crate::algorithm::blastn::disc_lookup::compute_discontiguous_index;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads the switches once (LOSAT_X_DISCLUTPAR, LOSAT_X_DISCLUTPARSHADOW;
/// off unless set to any value, not in the strict set; the shadow implies the switch). They
/// choose between two LOSAT code paths with the same result and change no value NCBI computes.
fn switches() -> (bool, bool) {
    static SWITCHES: OnceLock<(bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        let shadow = std::env::var_os("LOSAT_X_DISCLUTPARSHADOW").is_some();
        (
            shadow || std::env::var_os("LOSAT_X_DISCLUTPAR").is_some(),
            shadow,
        )
    })
}

/// No NCBI counterpart: the number of workers of the fill when `on`, i.e. the size of the
/// current search pool when there is one of two threads or more; `None` runs the original loop.
fn parts_for(on: bool) -> Option<usize> {
    if !on {
        return None;
    }
    pool_threads()
}

/// No NCBI counterpart: dispatch test of `build_disc_mb_lookup` (the switch, or its shadow, and
/// a search pool of two threads or more).
pub(super) fn parts() -> Option<usize> {
    parts_for(switches().0)
}

/// No NCBI counterpart: LOSAT_X_DISCLUTPARSHADOW is set.
pub(super) fn shadow() -> bool {
    switches().1
}

// No NCBI counterpart: counters of the shadow's summary line.
static SHADOW_TABLES: AtomicU64 = AtomicU64::new(0);
static SHADOW_POOL_TABLES: AtomicU64 = AtomicU64::new(0);
static SHADOW_TWO_TEMPLATE_TABLES: AtomicU64 = AtomicU64::new(0);
static SHADOW_CELLS: AtomicU64 = AtomicU64::new(0);
static SHADOW_DIFFS: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: the summed wall time of the builds, their number and the number built on a
// pool, for LOSAT_TIMING (any switch state); they change no value NCBI computes.
static BUILD_NANOS: AtomicU64 = AtomicU64::new(0);
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: LOSAT_TIMING is set (read once).
fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

/// No NCBI counterpart: the start of a build's `[TIMING]` interval (only with LOSAT_TIMING).
pub(super) fn timer_start() -> Option<std::time::Instant> {
    timing_on().then(std::time::Instant::now)
}

/// No NCBI counterpart: the end of a build's `[TIMING]` interval.
pub(super) fn record_build(start: Option<std::time::Instant>, on_pool: bool) {
    if let Some(start) = start {
        BUILD_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        BUILDS.fetch_add(1, Ordering::Relaxed);
        BUILDS_ON_POOL.fetch_add(u64::from(on_pool), Ordering::Relaxed);
    }
}

/// The `[TIMING]` line (with LOSAT_TIMING, when a table was built) and the shadow's summary line
/// (with LOSAT_X_DISCLUTPARSHADOW), printed by `main` at exit.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if timing_on() && load(&BUILDS) > 0 {
        eprintln!(
            "[TIMING] disc_lookup_build: {:.3}s (builds={} parallel={})",
            load(&BUILD_NANOS) as f64 / 1e9,
            load(&BUILDS),
            load(&BUILDS_ON_POOL)
        );
    }
    if shadow() {
        let diffs = load(&SHADOW_DIFFS);
        eprintln!(
            "[X_SHADOW] DISCLUTPAR tables={} pool={} two_templates={} cells={} diffs={diffs} {}",
            load(&SHADOW_TABLES),
            load(&SHADOW_POOL_TABLES),
            load(&SHADOW_TWO_TEMPLATE_TABLES),
            load(&SHADOW_CELLS),
            if diffs == 0 { "identical" } else { "DIFFERENT" }
        );
    }
}

/// The tables of one template, shared by the workers.
#[derive(Clone, Copy)]
struct SharedTemplate {
    hashtable: Shared<u32>,
    next_pos: Shared<u32>,
    helper_array: Shared<u32>,
}

#[derive(Clone, Copy)]
struct SharedTables {
    first: SharedTemplate,
    second: SharedTemplate,
    pv_array: Shared<PvArrayType>,
    pv_array_bts: usize,
}

impl SharedTables {
    fn new(t: DiscTables<'_>) -> Self {
        Self {
            first: SharedTemplate {
                hashtable: Shared::new(t.hashtable),
                next_pos: Shared::new(t.next_pos),
                helper_array: Shared::new(t.helper_array),
            },
            second: SharedTemplate {
                hashtable: Shared::new(t.hashtable2),
                next_pos: Shared::new(t.next_pos2),
                helper_array: Shared::new(t.helper_array2),
            },
            pv_array: Shared::new(t.pv_array),
            pv_array_bts: t.pv_array_bts,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:776-786,794-804
    // ```c
    //          if (mb_lt->hashtable[ecode1] == 0) {
    //             PV_SET(pv_array, ecode1, pv_array_bts);
    //          }
    //          else {
    //             helper_array[ecode1/kCompressionFactor]++;
    //          }
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
    //          mb_lt->hashtable[ecode1] = index;
    // ...
    //          if (mb_lt->hashtable2[ecode2] == 0) {
    //             PV_SET(pv_array, ecode2, pv_array_bts);
    //          }
    //          else {
    //             helper_array2[ecode2/kCompressionFactor]++;
    //          }
    //          mb_lt->next_pos2[index] = mb_lt->hashtable2[ecode2];
    //          mb_lt->hashtable2[ecode2] = index;
    // ```
    // The update of one word for one template (`second`), as in `fill_disc_mb_tables` (same
    // statements; the counters wrap as there), except that the value `next_pos[index]` gets (the
    // cell's previous head) is returned, for phase C to write; the presence bit is set in the
    // presence-vector word this worker owns.
    /// SAFETY: the worker owns the block of `ecode` (of the template's tables and of `pv_array`).
    #[inline(always)]
    unsafe fn update(self, second: bool, ecode: usize, index: u32) -> u32 {
        let t = if second { self.second } else { self.first };
        unsafe {
            let head = t.hashtable.get(ecode);
            if head == 0 {
                // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:49-50
                // ```c
                // #define PV_SET(lookup, index, shift) \
                //     lookup[(index) >> (shift)] |= (PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)
                // ```
                // As `pv_set_shift`, on the presence-vector word this worker owns.
                let array_idx = ecode >> self.pv_array_bts;
                if array_idx < self.pv_array.len {
                    let bit = (1 as PvArrayType) << (ecode & PV_ARRAY_MASK);
                    self.pv_array
                        .set(array_idx, self.pv_array.get(array_idx) | bit);
                }
            } else {
                let h = ecode >> K_COMPRESSION_FACTOR_SHIFT;
                t.helper_array.set(h, t.helper_array.get(h).wrapping_add(1));
            }
            t.hashtable.set(ecode, index);
            head
        }
    }
}

/// No NCBI counterpart (scheduling only): a word of one template in a bucket: its 1-based query
/// offset `index` and, until phase B updates it, its cell `ecode`; after that, the cell's previous
/// head, which is the value of `next_pos[index]`.
#[derive(Clone, Copy)]
struct Entry {
    index: u32,
    value: u32,
}

/// No NCBI counterpart (scheduling only): the buckets of one template, `parts x parts` lists of
/// words: bucket `(a, o)` holds the words that reader `a` read whose cell lies in owner `o`'s
/// blocks, in reading order. Shared by the jobs of a phase, each of which uses only its own
/// buckets: phase A's job `a` the row `(a, *)`, phase B's job `o` the column `(*, o)`, phase C's
/// job `a` the row `(a, *)`; the phases are separated by joins.
#[derive(Clone, Copy)]
struct SharedBuckets {
    ptr: *mut Vec<Entry>,
    parts: usize,
}

// SAFETY: the jobs of a phase use disjoint buckets (see `SharedBuckets`), the phases are separated
// by joins, and the buckets outlive the phases.
unsafe impl Send for SharedBuckets {}
// SAFETY: as above.
unsafe impl Sync for SharedBuckets {}

impl SharedBuckets {
    fn new(buckets: &mut [Vec<Entry>], parts: usize) -> Self {
        debug_assert!(buckets.is_empty() || buckets.len() == parts * parts);
        Self {
            ptr: buckets.as_mut_ptr(),
            parts,
        }
    }

    /// SAFETY: no other job of the current phase uses the bucket `(a, o)`, and the buckets outlive
    /// the returned borrow (the phase).
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    unsafe fn at<'b>(self, a: usize, o: usize) -> &'b mut Vec<Entry> {
        debug_assert!(a < self.parts && o < self.parts);
        unsafe { &mut *self.ptr.add(a * self.parts + o) }
    }
}

/// No NCBI counterpart (scheduling only): one unmasked range of one query (a location of
/// `s_FillDiscMBTable`) and the offset of its first position in the stream of all positions the
/// word loop reads, in reading order.
struct Segment<'a> {
    seq: &'a [u8],
    start: usize,
    end: usize,
    query_offset: usize,
    stream: usize,
}

/// The ranges of `fill_disc_mb_tables`'s loop, in its order (the same `build_unmasked_ranges`
/// per query), and the number of positions it reads.
fn segments<'a>(words: &DiscWords<'a>) -> (Vec<Segment<'a>>, usize) {
    let mut segments = Vec::new();
    let mut stream = 0usize;
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let query_offset = words.query_offsets[q_idx].max(0) as usize;
        for (start, end) in build_unmasked_ranges(seq.len(), masks) {
            if end <= start {
                continue;
            }
            segments.push(Segment {
                seq,
                start,
                end,
                query_offset,
                stream,
            });
            stream += end - start;
        }
    }
    (segments, stream)
}

/// No NCBI counterpart (scheduling only): the plan and the shared state of one fill.
struct Fill<'a, 'w> {
    words: &'a DiscWords<'w>,
    segments: &'a [Segment<'w>],
    owner: &'a [u8],
    block_shift: usize,
    parts: usize,
    tables: SharedTables,
    first: SharedBuckets,
    second: SharedBuckets,
}

impl Fill<'_, '_> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:749-767,775,793
    /// ```c
    ///       from = loc->ssr->left - (template_length - 2);
    ///       to = loc->ssr->right - (template_length - 2);
    ///       seq = query->sequence_start + loc->ssr->left;
    ///       pos = seq + template_length;
    ///
    ///       for (index = from; index <= to; index++) {
    ///          val = *++seq;
    ///          if ((val & BLAST2NA_MASK) != 0) {
    ///             accum = 0;
    ///             pos = seq + template_length;
    ///             continue;
    ///          }
    ///          accum = (accum << BITS_PER_NUC) | val;
    ///          if (seq < pos)
    ///             continue;
    /// ...
    ///          ecode1 = ComputeDiscontiguousIndex(accum, template_type);
    /// ...
    ///          ecode2 = ComputeDiscontiguousIndex(accum, second_template_type);
    /// ```
    /// Phase A, reader `a`: the words that end at the stream positions `[s0, s1)`, read by the
    /// loop of `fill_disc_mb_tables` (the same resets, `accum`, `index` and template indices),
    /// into the buckets `(a, owner of the cell)` of each template, in reading order. A range is
    /// entered up to `template_length - 1` positions before `s0` (not before its start), so that
    /// the template window of `accum` and the test `valid_bases < template_length` at the
    /// positions from `s0` on are those of the loop.
    ///
    /// SAFETY: no other job of the phase uses the buckets `(a, *)`.
    unsafe fn read(&self, a: usize, s0: usize, s1: usize) {
        let words = self.words;
        let two_templates = words.two_templates;
        let template_length = words.template_length;
        for o in 0..self.parts {
            // SAFETY: row `a` is this job's.
            unsafe {
                self.first.at(a, o).clear();
                if two_templates {
                    self.second.at(a, o).clear();
                }
            }
        }
        if s0 >= s1 {
            return;
        }
        let mut g = self
            .segments
            .partition_point(|seg| seg.stream + (seg.end - seg.start) <= s0);
        while g < self.segments.len() && self.segments[g].stream < s1 {
            let seg = &self.segments[g];
            let len = seg.end - seg.start;
            let first_pos = seg.start + (s0.max(seg.stream) - seg.stream);
            let end_pos = seg.start + (s1.min(seg.stream + len) - seg.stream);
            let warm_start = seg.start.max(first_pos.saturating_sub(template_length - 1));
            let mut accum: u64 = 0;
            let mut valid_bases = 0usize;
            for pos in warm_start..end_pos {
                let val = seg.seq[pos];
                if (val & BLAST2NA_MASK) != 0 {
                    accum = 0;
                    valid_bases = 0;
                    continue;
                }
                accum = (accum << 2) | val as u64;
                valid_bases += 1;
                if valid_bases < template_length || pos < first_pos {
                    continue;
                }
                // The 1-based start of the word in the concatenated query.
                let index = (seg.query_offset + pos + 2 - template_length) as u32;
                let ecode1 = compute_discontiguous_index(accum, words.template_type);
                let o1 = self.owner[(ecode1 as usize) >> self.block_shift] as usize;
                // SAFETY: row `a` is this job's.
                unsafe {
                    self.first.at(a, o1).push(Entry {
                        index,
                        value: ecode1,
                    })
                };
                if !two_templates {
                    continue;
                }
                let ecode2 = compute_discontiguous_index(accum, words.second_template_type);
                let o2 = self.owner[(ecode2 as usize) >> self.block_shift] as usize;
                // SAFETY: as above.
                unsafe {
                    self.second.at(a, o2).push(Entry {
                        index,
                        value: ecode2,
                    })
                };
            }
            g += 1;
        }
    }

    /// Phase B, owner `o`: the updates of the words of its cells, reader by reader (so in reading
    /// order), first template then second; each entry keeps the cell's previous head.
    ///
    /// SAFETY: no other job of the phase uses the buckets `(*, o)` or the cells of `o`'s blocks.
    unsafe fn update(&self, o: usize) {
        for a in 0..self.parts {
            // SAFETY: column `o` and the cells of `o`'s blocks are this job's.
            unsafe { update_all(self.tables, false, self.first.at(a, o)) };
        }
        if self.words.two_templates {
            for a in 0..self.parts {
                // SAFETY: as above.
                unsafe { update_all(self.tables, true, self.second.at(a, o)) };
            }
        }
    }

    /// Phase C, reader `a`: `next_pos[index]` (and `next_pos2[index]`) of the words it read.
    ///
    /// SAFETY: no other job of the phase uses the buckets `(a, *)`; the indices are the words' own.
    unsafe fn link(&self, a: usize) {
        for o in 0..self.parts {
            // SAFETY: row `a` is this job's; each index is written once.
            unsafe {
                for e in self.first.at(a, o).iter() {
                    self.tables.first.next_pos.set(e.index as usize, e.value);
                }
            }
        }
        if self.words.two_templates {
            for o in 0..self.parts {
                // SAFETY: as above.
                unsafe {
                    for e in self.second.at(a, o).iter() {
                        self.tables.second.next_pos.set(e.index as usize, e.value);
                    }
                }
            }
        }
    }
}

/// No NCBI counterpart (scheduling only): how far ahead phase B prefetches the cells of a bucket
/// on x86_64 (as the delay of `x_mb_lookup_par`'s fill).
#[cfg(target_arch = "x86_64")]
const X_AHEAD: usize = 16;

/// The updates of the entries of one bucket, in order; each entry's `value` becomes the cell's
/// previous head.
/// SAFETY: the cells of the entries lie in blocks of the calling job.
#[inline]
unsafe fn update_all(t: SharedTables, second: bool, entries: &mut [Entry]) {
    #[cfg(target_arch = "x86_64")]
    let cells = if second {
        t.second.hashtable
    } else {
        t.first.hashtable
    };
    for j in 0..entries.len() {
        #[cfg(target_arch = "x86_64")]
        if let Some(ahead) = entries.get(j + X_AHEAD) {
            debug_assert!((ahead.value as usize) < cells.len);
            // SAFETY: a cell index < the table size; a prefetch does not fault.
            unsafe {
                core::arch::x86_64::_mm_prefetch(
                    cells.ptr.add(ahead.value as usize) as *const i8,
                    core::arch::x86_64::_MM_HINT_T0,
                );
            }
        }
        let e = entries[j];
        // SAFETY: the cell is in the caller's blocks.
        entries[j].value = unsafe { t.update(second, e.value as usize, e.index) };
    }
}

/// No NCBI counterpart (scheduling only): the positions read per chunk (each chunk is read, updated
/// and linked before the next; it bounds the buckets to 8 bytes per word of a chunk per template).
const CHUNK_POSITIONS: usize = 1 << 19;

/// No NCBI counterpart (scheduling only): `job(0)`, ..., `job(parts - 1)`, on the pool of the
/// calling worker (`in_pool`) or one after another; returns when all have returned. The jobs never
/// wait for anything.
fn run_jobs<F: Fn(usize) + Sync>(parts: usize, in_pool: bool, job: F) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if in_pool && parts > 1 {
        // The scope runs on the pool of the calling worker; outside a pool it would start
        // rayon's global pool, which `parts` never asks for.
        debug_assert!(rayon::current_thread_index().is_some());
        let job = &job;
        rayon::scope(|scope| {
            for i in 1..parts {
                scope.spawn(move |_| job(i));
            }
            job(0);
        });
        return;
    }
    let _ = in_pool;
    for i in 0..parts {
        job(i);
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:733-806
/// ```c
///    for (loc = location; loc; loc = loc->next) {
/// ...
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
///          mb_lt->hashtable[ecode1] = index;
/// ...
///          mb_lt->next_pos2[index] = mb_lt->hashtable2[ecode2];
///          mb_lt->hashtable2[ecode2] = index;
/// ```
/// Fills the tables as the word loop of `fill_disc_mb_tables` does, with `parts` workers on the
/// search pool (or one after another, `in_pool == false`, for the shadow without a pool).
pub(super) fn fill(parts: usize, words: &DiscWords<'_>, tables: DiscTables<'_>, in_pool: bool) {
    fill_with(parts, words, tables, in_pool, CHUNK_POSITIONS, &|_| {});
}

/// `fill` with the chunk size and a `probe` that each phase-B job calls first, on its thread.
fn fill_with<P: Fn(u8) + Sync>(
    parts: usize,
    words: &DiscWords<'_>,
    tables: DiscTables<'_>,
    in_pool: bool,
    chunk: usize,
    probe: &P,
) {
    let table_size = tables.hashtable.len();
    debug_assert!(!words.two_templates || tables.hashtable2.len() == table_size);
    let (block_shift, owner) = plan(table_size, tables.pv_array_bts, parts);
    debug_assert!(
        block_shift >= K_COMPRESSION_FACTOR_SHIFT && block_shift >= tables.pv_array_bts,
        "a block must own whole counters and presence-vector words"
    );
    let parts = parts.clamp(1, 255);
    let chunk = chunk.max(1);
    let (segments, total) = segments(words);
    let mut first: Vec<Vec<Entry>> = (0..parts * parts).map(|_| Vec::new()).collect();
    let mut second: Vec<Vec<Entry>> = if words.two_templates {
        (0..parts * parts).map(|_| Vec::new()).collect()
    } else {
        Vec::new()
    };
    let fill = Fill {
        words,
        segments: &segments,
        owner: &owner,
        block_shift,
        parts,
        tables: SharedTables::new(tables),
        first: SharedBuckets::new(&mut first, parts),
        second: SharedBuckets::new(&mut second, parts),
    };
    let mut start = 0usize;
    while start < total {
        let end = total.min(start.saturating_add(chunk));
        let span = |a: usize| {
            (
                start + (end - start) * a / parts,
                start + (end - start) * (a + 1) / parts,
            )
        };
        // SAFETY (all three phases): job `i` of a phase uses only its own buckets and cells (see
        // `SharedBuckets`, `Fill::read`, `Fill::update`, `Fill::link`); `run_jobs` returns after
        // every job of the phase has returned.
        run_jobs(parts, in_pool, |a| {
            let (s0, s1) = span(a);
            unsafe { fill.read(a, s0, s1) }
        });
        run_jobs(parts, in_pool, |o| {
            probe(o as u8);
            unsafe { fill.update(o) }
        });
        run_jobs(parts, in_pool, |a| unsafe { fill.link(a) });
        start = end;
    }
}

/// No NCBI counterpart (scheduling only): the parts of the shadow's partitioned build when no
/// pool built the tables (as `x_mb_lookup_par::shadow_check`).
const SHADOW_PARTS_OFF_POOL: usize = 4;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:785-786,803-804,808-824
/// ```c
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
///          mb_lt->hashtable[ecode1] = index;
/// ...
///          mb_lt->next_pos2[index] = mb_lt->hashtable2[ecode2];
///          mb_lt->hashtable2[ecode2] = index;
/// ...
///    longest_chain = 2;
///    for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
///        longest_chain = MAX(longest_chain, helper_array[index]);
/// ```
/// LOSAT_X_DISCLUTPARSHADOW: builds the tables again the other way (with the original loop when
/// the pool built `built`, `pool_parts` = `Some`; otherwise partitioned in four parts one after
/// another) and compares every element of every array and the longest chain; panics at the first
/// difference after printing it.
pub(super) fn shadow_check(
    words: &DiscWords<'_>,
    pool_parts: Option<usize>,
    built: &DiscTablesRef<'_>,
    longest_chain: usize,
) {
    let mut hashtable = vec![0u32; built.hashtable.len()];
    let mut next_pos = vec![0u32; built.next_pos.len()];
    let mut helper_array = vec![0u32; built.helper_array.len()];
    let mut hashtable2 = vec![0u32; built.hashtable2.len()];
    let mut next_pos2 = vec![0u32; built.next_pos2.len()];
    let mut helper_array2 = vec![0u32; built.helper_array2.len()];
    let mut pv_array = vec![0 as PvArrayType; built.pv_array.len()];
    let tables = DiscTables {
        hashtable: &mut hashtable,
        next_pos: &mut next_pos,
        helper_array: &mut helper_array,
        hashtable2: &mut hashtable2,
        next_pos2: &mut next_pos2,
        helper_array2: &mut helper_array2,
        pv_array: &mut pv_array,
        pv_array_bts: built.pv_array_bts,
    };
    let other = match pool_parts {
        Some(_) => {
            fill_disc_mb_tables(words, tables);
            "serial"
        }
        None => {
            fill(SHADOW_PARTS_OFF_POOL, words, tables, false);
            "partitioned"
        }
    };
    let other_chain = disc_longest_chain(&helper_array, &helper_array2, words.two_templates);
    let mut difference = None;
    for (name, have, again) in [
        ("hashtable", built.hashtable, hashtable.as_slice()),
        ("next_pos", built.next_pos, next_pos.as_slice()),
        ("helper_array", built.helper_array, helper_array.as_slice()),
        ("hashtable2", built.hashtable2, hashtable2.as_slice()),
        ("next_pos2", built.next_pos2, next_pos2.as_slice()),
        (
            "helper_array2",
            built.helper_array2,
            helper_array2.as_slice(),
        ),
        ("pv_array", built.pv_array, pv_array.as_slice()),
    ] {
        if let Some(i) = (0..have.len()).find(|&i| have[i] != again[i]) {
            difference = Some(format!(
                "{name}[{i}]: built {} {other} {}",
                have[i], again[i]
            ));
            break;
        }
    }
    if difference.is_none() && longest_chain != other_chain {
        difference = Some(format!(
            "longest_chain: built {longest_chain} {other} {other_chain}"
        ));
    }
    if let Some(what) = difference {
        SHADOW_DIFFS.fetch_add(1, Ordering::Relaxed);
        eprintln!(
            "[X_SHADOW] DISCLUTPAR DIFFERENT {what} (pool parts {pool_parts:?}, compared with the {other} build; template_length {}, two_templates {})",
            words.template_length, words.two_templates
        );
        panic!("LOSAT_X_DISCLUTPARSHADOW: the discontiguous megablast tables differ: {what}");
    }
    SHADOW_TABLES.fetch_add(1, Ordering::Relaxed);
    SHADOW_POOL_TABLES.fetch_add(u64::from(pool_parts.is_some()), Ordering::Relaxed);
    SHADOW_TWO_TEMPLATE_TABLES.fetch_add(u64::from(words.two_templates), Ordering::Relaxed);
    SHADOW_CELLS.fetch_add(
        [
            built.hashtable.len(),
            built.next_pos.len(),
            built.helper_array.len(),
            built.hashtable2.len(),
            built.next_pos2.len(),
            built.helper_array2.len(),
            built.pv_array.len(),
        ]
        .iter()
        .sum::<usize>() as u64,
        Ordering::Relaxed,
    );
}

#[cfg(test)]
mod tests {
    use super::super::{build_disc_mb_lookup, compute_mb_pv_params};
    use super::*;
    use crate::algorithm::blastn::disc_lookup::{get_disc_template_type, DiscWordType};
    use crate::utils::dust::MaskedInterval;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:785-786,803-804
    // ```c
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
    //          mb_lt->hashtable[ecode1] = index;
    // ...
    //          mb_lt->next_pos2[index] = mb_lt->hashtable2[ecode2];
    //          mb_lt->hashtable2[ecode2] = index;
    // ```
    // Random queries (with runs of ambiguity codes, masks, queries shorter than the template and
    // queries without a mask list) through the original loop and through the partitioned fill:
    // every element of every table must be the same.
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

    struct Case {
        queries: Vec<Vec<u8>>,
        offsets: Vec<i32>,
        masks: Vec<Vec<MaskedInterval>>,
    }

    /// One to four queries in BLASTNA: lengths 1..3,000 (one in four shorter than any template),
    /// half of the letters from a short repeated motif (long chains), 2 % ambiguity codes (some in
    /// runs), random masks or none; the last query may have no mask list at all.
    fn random_case(rng: &mut Lcg) -> Case {
        let n_queries = 1 + rng.below(4) as usize;
        let mut queries = Vec::new();
        let mut offsets = Vec::new();
        let mut masks = Vec::new();
        let mut offset = 0i32;
        for _ in 0..n_queries {
            let len = if rng.below(4) == 0 {
                1 + rng.below(20) as usize
            } else {
                1 + rng.below(3000) as usize
            };
            let motif: Vec<u8> = (0..1 + rng.below(40)).map(|_| rng.below(4) as u8).collect();
            let mut seq = Vec::with_capacity(len);
            while seq.len() < len {
                let r = rng.below(100);
                if r < 1 {
                    let run = 1 + rng.below(30) as usize;
                    seq.extend(std::iter::repeat_n(14u8, run)); // N
                } else if r < 2 {
                    seq.push(4 + rng.below(12) as u8); // another ambiguity code
                } else if r < 50 {
                    seq.push(motif[seq.len() % motif.len()]);
                } else {
                    seq.push(rng.below(4) as u8);
                }
            }
            seq.truncate(len);
            let mut m = Vec::new();
            let mut cursor = 0usize;
            while cursor < len && rng.below(3) != 0 {
                let start = cursor + rng.below(400) as usize;
                let end = start + 1 + rng.below(200) as usize;
                if start >= len {
                    break;
                }
                m.push(MaskedInterval {
                    start,
                    end: end.min(len),
                });
                cursor = end;
            }
            offsets.push(offset);
            offset += len as i32 + 1;
            queries.push(seq);
            masks.push(m);
        }
        if rng.below(4) == 0 {
            masks.pop();
        }
        Case {
            queries,
            offsets,
            masks,
        }
    }

    /// Every (word size, template type, template length) that a discontiguous table has.
    fn combos() -> Vec<(usize, DiscWordType, usize)> {
        let mut v = Vec::new();
        for word_size in [11usize, 12] {
            for word_type in [
                DiscWordType::Coding,
                DiscWordType::Optimal,
                DiscWordType::TwoTemplates,
            ] {
                for template_length in [16usize, 18, 21] {
                    v.push((word_size, word_type, template_length));
                }
            }
        }
        v
    }

    #[derive(Debug, PartialEq, Eq)]
    struct Built {
        hashtable: Vec<u32>,
        next_pos: Vec<u32>,
        helper_array: Vec<u32>,
        hashtable2: Vec<u32>,
        next_pos2: Vec<u32>,
        helper_array2: Vec<u32>,
        pv_array: Vec<PvArrayType>,
        longest_chain: usize,
    }

    enum How {
        Serial,
        /// `fill` (the default chunk size), on the pool when `true`.
        Parts(usize, bool),
        /// `fill_with` with that many parts and positions per chunk, on the pool when `true`.
        Chunked(usize, usize, bool),
    }

    /// The tables of `case` allocated as `build_disc_mb_lookup` allocates them and filled `how`.
    fn build(case: &Case, combo: (usize, DiscWordType, usize), how: How) -> Built {
        let (word_size, word_type, template_length) = combo;
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let hashsize = 1usize << (2 * word_size);
        let approx: usize = case.queries.iter().map(Vec::len).sum();
        let (pv_size, pv_array_bts) = compute_mb_pv_params(hashsize, approx, false, word_size);
        let template_type =
            get_disc_template_type(word_size as i32, template_length as u8, word_type);
        let two_templates = word_type == DiscWordType::TwoTemplates;
        let second_template_type = if two_templates {
            template_type.next()
        } else {
            crate::algorithm::blastn::disc_lookup::DiscTemplateType::Contiguous
        };
        let query_length = case
            .queries
            .iter()
            .zip(&case.offsets)
            .map(|(q, &o)| o as usize + q.len())
            .max()
            .unwrap_or(0);
        let second = |n: usize| {
            if two_templates {
                vec![0u32; n]
            } else {
                Vec::new()
            }
        };
        let mut b = Built {
            hashtable: vec![0; hashsize],
            next_pos: vec![0; query_length + 1],
            helper_array: vec![0; hashsize / 2048],
            hashtable2: second(hashsize),
            next_pos2: second(query_length + 1),
            helper_array2: second(hashsize / 2048),
            pv_array: vec![0; pv_size],
            longest_chain: 0,
        };
        let words = DiscWords {
            queries_blastna: &query_refs,
            query_offsets: &case.offsets,
            query_masks: &case.masks,
            template_length,
            template_type,
            second_template_type,
            two_templates,
        };
        let tables = DiscTables {
            hashtable: &mut b.hashtable,
            next_pos: &mut b.next_pos,
            helper_array: &mut b.helper_array,
            hashtable2: &mut b.hashtable2,
            next_pos2: &mut b.next_pos2,
            helper_array2: &mut b.helper_array2,
            pv_array: &mut b.pv_array,
            pv_array_bts,
        };
        match how {
            How::Serial => fill_disc_mb_tables(&words, tables),
            How::Parts(parts, in_pool) => fill(parts, &words, tables, in_pool),
            How::Chunked(parts, chunk, in_pool) => {
                fill_with(parts, &words, tables, in_pool, chunk, &|_| {})
            }
        }
        b.longest_chain = disc_longest_chain(&b.helper_array, &b.helper_array2, two_templates);
        b
    }

    /// The serial tables have work in them: chains (counters above zero), and for two templates
    /// both tables.
    fn assert_not_trivial(b: &Built, two_templates: bool) {
        assert!(b.hashtable.iter().any(|&c| c != 0));
        assert!(b.pv_array.iter().any(|&w| w != 0));
        if two_templates {
            assert!(b.hashtable2.iter().any(|&c| c != 0));
        }
    }

    // The helper above allocates as `build_disc_mb_lookup` does: its serial tables are the
    // function's (which runs the original loop on this test thread, which is no pool worker).
    #[test]
    fn helper_builds_what_build_disc_mb_lookup_builds() {
        let mut rng = Lcg(0xD15C_0001);
        for combo in combos() {
            let case = random_case(&mut rng);
            let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
            let approx: usize = case.queries.iter().map(Vec::len).sum();
            let lookup = build_disc_mb_lookup(
                &query_refs,
                &case.offsets,
                combo.0,
                combo.2,
                combo.1,
                &case.masks,
                approx,
            );
            let b = build(&case, combo, How::Serial);
            let disc = lookup.disc().unwrap();
            assert_eq!(lookup.mb_lookup.hashtable, b.hashtable, "{combo:?}");
            assert_eq!(lookup.mb_lookup.next_pos, b.next_pos, "{combo:?}");
            assert_eq!(lookup.mb_lookup.pv_array, b.pv_array, "{combo:?}");
            assert_eq!(disc.hashtable2, b.hashtable2, "{combo:?}");
            assert_eq!(disc.next_pos2, b.next_pos2, "{combo:?}");
            assert_eq!(lookup.longest_chain(), b.longest_chain, "{combo:?}");
        }
    }

    // Partitioned in P = 1, 2, 3, 7 and 16 parts one after another (every template type and
    // length, one and two templates), with the default chunk and with chunks of 1 to 4,096
    // positions (one per part count, rotating by round): the same tables as the original loop.
    #[test]
    fn partitioned_fill_builds_the_same_tables() {
        let mut rng = Lcg(0xD15C_5EED);
        let mut chains = 0;
        for round in 0..3 {
            for combo in combos() {
                let case = random_case(&mut rng);
                let reference = build(&case, combo, How::Serial);
                if round == 0 && case.queries.iter().map(Vec::len).sum::<usize>() > 200 {
                    assert_not_trivial(&reference, combo.1 == DiscWordType::TwoTemplates);
                }
                chains += usize::from(reference.longest_chain > 3);
                for (k, parts) in [1usize, 2, 3, 7, 16].into_iter().enumerate() {
                    let got = build(&case, combo, How::Parts(parts, false));
                    assert!(
                        got == reference,
                        "round {round} {combo:?} parts {parts}: tables differ"
                    );
                    let chunk = [1usize, 17, 21, 333, 4096][(k + round) % 5];
                    let got = build(&case, combo, How::Chunked(parts, chunk, false));
                    assert!(
                        got == reference,
                        "round {round} {combo:?} parts {parts} chunk {chunk}: tables differ"
                    );
                }
            }
        }
        assert!(chains > 10, "the cases must have long chains ({chains})");
    }

    // Chunk and span boundaries anywhere relative to the words (chunks of 1 to 1,000 positions, so
    // spans start inside words, at ambiguity codes and at range starts, and many spans are empty):
    // the same tables as the original loop, for 1, 3 and 16 parts.
    #[test]
    fn chunk_boundaries_anywhere() {
        let mut rng = Lcg(0xC4A2_B0DE);
        for case_no in 0..6 {
            let case = random_case(&mut rng);
            let combo = [
                (11usize, DiscWordType::Coding, 18usize),
                (11, DiscWordType::TwoTemplates, 21),
                (11, DiscWordType::Optimal, 16),
            ][case_no % 3];
            let reference = build(&case, combo, How::Serial);
            for chunk in [1usize, 2, 5, 17, 18, 19, 20, 21, 22, 64, 1000] {
                for parts in [1usize, 3, 16] {
                    let got = build(&case, combo, How::Chunked(parts, chunk, false));
                    assert!(
                        got == reference,
                        "case {case_no} {combo:?} chunk {chunk} parts {parts}: tables differ"
                    );
                }
            }
        }
    }

    // Queries shorter than every template and empty mask lists (no mask list at all, and lists
    // without intervals), and a query that is all ambiguity codes: no word, and the same (empty)
    // tables either way.
    #[test]
    fn short_and_unmasked_queries() {
        let case = Case {
            queries: vec![vec![0, 1, 2, 3, 0, 1, 2], vec![14; 40], vec![3; 15]],
            offsets: vec![0, 8, 49],
            masks: vec![Vec::new()],
        };
        for combo in combos() {
            let reference = build(&case, combo, How::Serial);
            assert!(reference.hashtable.iter().all(|&c| c == 0), "{combo:?}");
            for parts in [2usize, 3, 7, 16] {
                assert!(build(&case, combo, How::Parts(parts, false)) == reference);
            }
        }
        // One long query without masks (`masks` empty), then the same with an empty list.
        let mut rng = Lcg(0xD15C_0E11);
        let seq: Vec<u8> = (0..5000).map(|_| rng.below(4) as u8).collect();
        for masks in [Vec::new(), vec![Vec::new()]] {
            let case = Case {
                queries: vec![seq.clone()],
                offsets: vec![0],
                masks,
            };
            for combo in combos() {
                let reference = build(&case, combo, How::Serial);
                assert_not_trivial(&reference, combo.1 == DiscWordType::TwoTemplates);
                for parts in [2usize, 3, 7, 16] {
                    assert!(build(&case, combo, How::Parts(parts, false)) == reference);
                }
            }
        }
    }

    // On real pools (`with_search_pool`, the search's pool) of 1, 2, 3 and 8 threads: the
    // dispatch gives the pool's size from 2 threads on and nothing on 1, and the fill on the pool
    // (also with more parts than threads, and with chunks of 97 positions) builds the same tables
    // as the original loop.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn fill_on_a_pool_builds_the_same_tables() {
        let mut rng = Lcg(0xFACE_D15C);
        let cases: Vec<Case> = (0..6).map(|_| random_case(&mut rng)).collect();
        for threads in [1usize, 2, 3, 8] {
            crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                let parts = parts_for(true);
                assert_eq!(parts, (threads > 1).then_some(threads), "threads {threads}");
                for (i, combo) in combos().into_iter().enumerate() {
                    let case = &cases[i % cases.len()];
                    let reference = build(case, combo, How::Serial);
                    let got = match parts {
                        Some(parts) => build(case, combo, How::Parts(parts, true)),
                        None => build(case, combo, How::Serial),
                    };
                    assert!(got == reference, "threads {threads} {combo:?}");
                    if threads > 1 {
                        let more = build(case, combo, How::Parts(16, true));
                        assert!(more == reference, "threads {threads} {combo:?} 16 parts");
                        let small = build(case, combo, How::Chunked(threads, 97, true));
                        assert!(
                            small == reference,
                            "threads {threads} {combo:?} chunks of 97"
                        );
                    }
                }
                Ok(())
            })
            .unwrap();
        }
    }

    // The dispatch: nothing with the switch off, nothing off a pool (this test thread is no
    // rayon worker), nothing on a search "pool" of one thread, the pool's size from two threads.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn dispatch_needs_a_search_pool_of_two_threads() {
        assert_eq!(rayon::current_thread_index(), None);
        assert_eq!(parts_for(false), None);
        assert_eq!(parts_for(true), None);
        for (threads, want) in [(1usize, None), (2, Some(2)), (5, Some(5))] {
            crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                assert_eq!(parts_for(true), want, "threads {threads}");
                assert_eq!(parts_for(false), None, "threads {threads}");
                Ok(())
            })
            .unwrap();
        }
    }

    /// The fill on a pool of three threads runs its parts on the workers of that pool, and on
    /// more than one of them: the probe that each phase-B job (owner) calls first sees a rayon
    /// worker of a pool of three threads, and the owners are spread over at least two workers. (Each probe call
    /// waits until a second worker has called it (the first call that waits 10 s in vain ends
    /// all waiting), so that the spread does not depend on the scheduler; a fill that ran its
    /// parts one after another on one thread has no second worker and fails after those 10 s.)
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn the_fill_runs_on_several_workers_of_the_pool() {
        use std::collections::BTreeSet;
        use std::sync::atomic::{AtomicBool, AtomicUsize};
        use std::sync::Mutex;
        use std::time::{Duration, Instant};

        let mut rng = Lcg(0x0B5E_12FE);
        let case = Case {
            queries: vec![(0..20_000).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        };
        let combo = (11usize, DiscWordType::TwoTemplates, 18usize);
        let reference = build(&case, combo, How::Serial);
        crate::utils::threading::with_search_pool(3, "blastn", |_pool| {
            let parts = parts_for(true).expect("a pool of three threads");
            assert_eq!(parts, 3);
            let calls = AtomicUsize::new(0);
            let in_pool_of_three = AtomicUsize::new(0);
            let workers = Mutex::new(BTreeSet::new());
            let gave_up = AtomicBool::new(false);
            let probe = |_part: u8| {
                calls.fetch_add(1, Ordering::SeqCst);
                if let Some(index) = rayon::current_thread_index() {
                    workers.lock().unwrap().insert(index);
                    if rayon::current_num_threads() == 3 {
                        in_pool_of_three.fetch_add(1, Ordering::SeqCst);
                    }
                }
                let start = Instant::now();
                while !gave_up.load(Ordering::SeqCst) && workers.lock().unwrap().len() < 2 {
                    if start.elapsed() > Duration::from_secs(10) {
                        gave_up.store(true, Ordering::SeqCst);
                    }
                    std::thread::yield_now();
                }
            };
            let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
            let mut got = build(&case, combo, How::Serial);
            got.hashtable.fill(0);
            got.next_pos.fill(0);
            got.helper_array.fill(0);
            got.hashtable2.fill(0);
            got.next_pos2.fill(0);
            got.helper_array2.fill(0);
            got.pv_array.fill(0);
            let (_, pv_array_bts) = compute_mb_pv_params(1 << 22, 20_000, false, 11);
            let template_type = get_disc_template_type(11, 18, DiscWordType::TwoTemplates);
            fill_with(
                parts,
                &DiscWords {
                    queries_blastna: &query_refs,
                    query_offsets: &case.offsets,
                    query_masks: &case.masks,
                    template_length: 18,
                    template_type,
                    second_template_type: template_type.next(),
                    two_templates: true,
                },
                DiscTables {
                    hashtable: &mut got.hashtable,
                    next_pos: &mut got.next_pos,
                    helper_array: &mut got.helper_array,
                    hashtable2: &mut got.hashtable2,
                    next_pos2: &mut got.next_pos2,
                    helper_array2: &mut got.helper_array2,
                    pv_array: &mut got.pv_array,
                    pv_array_bts,
                },
                true,
                CHUNK_POSITIONS,
                &probe,
            );
            got.longest_chain = disc_longest_chain(&got.helper_array, &got.helper_array2, true);
            assert!(got == reference, "the tables differ");
            assert_eq!(calls.load(Ordering::SeqCst), 3);
            assert_eq!(in_pool_of_three.load(Ordering::SeqCst), 3);
            let workers = workers.lock().unwrap();
            assert!(
                workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                "workers {workers:?}"
            );
            Ok(())
        })
        .unwrap();
    }

    /// The tables of `b` as the shadow reads them.
    fn tables_ref(b: &Built, pv_array_bts: usize) -> DiscTablesRef<'_> {
        DiscTablesRef {
            hashtable: &b.hashtable,
            next_pos: &b.next_pos,
            helper_array: &b.helper_array,
            hashtable2: &b.hashtable2,
            next_pos2: &b.next_pos2,
            helper_array2: &b.helper_array2,
            pv_array: &b.pv_array,
            pv_array_bts,
        }
    }

    /// `DiscWords` of `case` with the templates of `combo`.
    fn words_of<'a>(
        case: &'a Case,
        query_refs: &'a [&'a [u8]],
        combo: (usize, DiscWordType, usize),
    ) -> DiscWords<'a> {
        let template_type = get_disc_template_type(combo.0 as i32, combo.2 as u8, combo.1);
        let two_templates = combo.1 == DiscWordType::TwoTemplates;
        DiscWords {
            queries_blastna: query_refs,
            query_offsets: &case.offsets,
            query_masks: &case.masks,
            template_length: combo.2,
            template_type,
            second_template_type: if two_templates {
                template_type.next()
            } else {
                crate::algorithm::blastn::disc_lookup::DiscTemplateType::Contiguous
            },
            two_templates,
        }
    }

    // The shadow passes equal tables both ways (built on a pool: compared with the original loop;
    // built serially: compared with the partitioned fill) and counts them.
    #[test]
    fn shadow_check_passes_equal_tables() {
        let mut rng = Lcg(0x5AAD_0001);
        let case = random_case(&mut rng);
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let approx: usize = case.queries.iter().map(Vec::len).sum();
        for combo in [
            (11usize, DiscWordType::Coding, 18usize),
            (12, DiscWordType::TwoTemplates, 21),
        ] {
            let (_, bts) = compute_mb_pv_params(1 << (2 * combo.0), approx, false, combo.0);
            let words = words_of(&case, &query_refs, combo);
            let before = SHADOW_TABLES.load(Ordering::Relaxed);
            let two_before = SHADOW_TWO_TEMPLATE_TABLES.load(Ordering::Relaxed);
            let serial = build(&case, combo, How::Serial);
            shadow_check(
                &words,
                None,
                &tables_ref(&serial, bts),
                serial.longest_chain,
            );
            let parted = build(&case, combo, How::Parts(3, false));
            shadow_check(
                &words,
                Some(3),
                &tables_ref(&parted, bts),
                parted.longest_chain,
            );
            assert!(SHADOW_TABLES.load(Ordering::Relaxed) >= before + 2);
            if combo.1 == DiscWordType::TwoTemplates {
                assert!(SHADOW_TWO_TEMPLATE_TABLES.load(Ordering::Relaxed) >= two_before + 2);
            }
        }
    }

    // A table with one element changed (a second-template chain link, the last element of
    // `next_pos2`) makes the shadow panic.
    #[test]
    #[should_panic(expected = "LOSAT_X_DISCLUTPARSHADOW")]
    fn shadow_check_panics_on_a_difference() {
        let mut rng = Lcg(0x5AAD_0002);
        let case = Case {
            queries: vec![(0..3000).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        };
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let combo = (11usize, DiscWordType::TwoTemplates, 16usize);
        let (_, bts) = compute_mb_pv_params(1 << 22, 3000, false, 11);
        let words = words_of(&case, &query_refs, combo);
        let mut b = build(&case, combo, How::Serial);
        let last = b.next_pos2.len() - 1;
        b.next_pos2[last] ^= 1;
        shadow_check(&words, Some(2), &tables_ref(&b, bts), b.longest_chain);
    }
}
