//! EXPERIMENT (LOSAT_X_DISCLUTPAR, in the strict set: default ON since 2026-10-17, off with
//! `LOSAT_X_DISCLUTPAR=0`): the
//! discontiguous megablast tables of `build_disc_mb_lookup` (`s_FillDiscMBTable`) filled by the
//! threads of the search pool, each thread owning a fixed set of table cells of both templates'
//! tables. The plan of the cells and the pool test are those of PAIRPAR's megablast table
//! (`x_mb_lookup_par.rs`).
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
//! three phases and a prefault round, each a set of P jobs that returns when all have returned:
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
//! - P (prefault, see "The prefault round" below; in the first chunk only, when it is large
//!   enough): the jobs write 0 to one element of each page of the cells and of the presence
//!   vector that B of the first chunk writes.
//! - B (update): owner `o` takes the buckets `(0, o)`, `(1, o)`, ..., `(P-1, o)`: every word of the
//!   chunk whose cell lies in its blocks, in reading order (the spans follow each other). For each
//!   it makes the loop's update of the cell (the presence bit when the cell is empty, else the
//!   counter, in owner `o`'s private copy of the counters; then the cell takes `index`) and keeps
//!   the cell's previous head in the entry. The chunks are done in order, so each cell of
//!   `hashtable` (and of `hashtable2`) receives the same words in the same order as in the loop:
//!   the same final head, the same increments of the same counters, and each word's kept head is
//!   the value the loop writes to `next_pos[index]`. Both templates set bits in the same
//!   `pv_array`, but a presence-vector word lies in the block of the cells whose bits it holds,
//!   whichever template sets them (same `hashsize`, same plan), so it has one writer; a bit is
//!   only ever set, never cleared or read during the fill, so the order of the two templates' sets
//!   does not matter.
//! - C (link): reader `a` writes `next_pos[index]` (and `next_pos2[index]`) of the words it read
//!   from their kept heads. Each index is written once; a span's words have a contiguous run of
//!   query offsets, so the readers write disjoint parts of `next_pos`.
//!
//! After the last chunk the calling thread adds each owner's private counters into `helper_array`
//! (and `helper_array2`) (`OwnerCounters::flush`), and `longest_chain` is computed afterwards from
//! them by the original code (`lookup::disc_longest_chain`). Jobs write disjoint elements, so the
//! number of workers changes the running time only. No floating-point operation is involved.
//!
//! Why the counters are the loop's. Counter `h` of `helper_array` counts the cells
//! `h * 2048 .. (h + 1) * 2048`, which lie in one block (`2^block_shift` is a multiple of 2048), so
//! one owner `o` makes all its increments: one wrapping `+ 1` for each word of the template whose
//! cell lies there and is not empty when the word comes, in phase B, into its private copy, which
//! starts at 0 and which no other job touches. Those are the loop's words in the loop's order of
//! each cell (above), so after the last chunk the copy holds their number `n` modulo 2^32 and the
//! flush makes the counter `v.wrapping_add(n)`: the loop's `v + 1 + ... + 1` modulo 2^32 from the
//! same start `v` (0 in every build: the arrays are allocated zeroed, as NCBI's `calloc`). No job
//! reads or writes `helper_array` or `helper_array2` (the jobs' `SharedTables` holds neither), the
//! flush runs after the last join, and their only readers are after the fill
//! (`build_disc_mb_lookup`: `disc_longest_chain`; the arrays are dropped there).
//!
//! Why three phases (scheduling only): in a first version every worker read every word and applied
//! those of its cells, writing `next_pos[index]` itself; each worker then paid the whole read
//! (about 0.09 s for 10 million positions at the workstation's clock), and neighbouring
//! `next_pos` entries were written by different workers (a cache line moving between cores for
//! most writes): at 4 threads that build was no faster than the loop (ws, `logs/disclutpar_report.md`).
//! Here each position is read once, each owner reads only its own words, and `next_pos` is
//! written by spans; phase B prefetches the cells 16 entries ahead on x86_64.
//!
//! Why padded bucket headers and private counters (scheduling only). In the version of #252 the
//! jobs' summed time was 3.8 times (4 threads) and 4.9 times (16) the loop's fill, and the cause
//! was false sharing (workstation, perf counters inside the fill; `logs/disclutcause_report.md` of
//! perf-mb-20261014): (1) in phase B the 4-byte counters of `helper_array`, one per block of 2048
//! cells with the blocks dealt round robin, put the counters of 16 owners (P = 16) or 4 (P = 4)
//! on each 64-byte line, and every word on a non-empty cell wrote one; (2) in phase A the 24-byte
//! `Vec` headers of the buckets, a row of P per reader, shared lines with the neighbouring
//! readers' rows, and every push stores the length. That gave, per fill of 1,775 cDNAs, 3.2-4.5
//! million HITM loads (loads served by a line another core had modified) and as many RFOs
//! (stores that had to take a line from another core), 6.4-9.0 million of each for 3,492 cDNAs,
//! against at most a few thousand in the loop. The jobs' summed time (3.8 and 4.9 times the
//! loop's) is stall time, not work: the instructions were 1.46-1.6 times the loop's. Now each
//! bucket is a `Bucket` on a 64-byte line of its own (layout only: the same pushes into the same
//! buckets in the same order), and each owner counts into its own copy of the counters on lines
//! of their own (`OwnerCounters`), added into `helper_array` once after the last chunk. So the
//! counter lines and the bucket-header lines are no longer shared. The presence-vector words of
//! two neighbouring blocks still share a line (eight words, 32 bytes, per block at the
//! workstation's `pv_array_bts` 8); that kept 1.1-2.4 million HITM loads and was measured to
//! cost nothing (the report's `pvlate` test), so they were left alone. Measured as the probe of
//! that report: the fill 0.205 -> 0.088 s at 4 threads and 0.0835 -> 0.047 s at 16 on 1,775
//! cDNAs, and faster than the loop at 2 threads too, so the switch needs no thread threshold.
//!
//! Bucket memory. The buckets hold one entry of 8 bytes per word of the current chunk and
//! template, and a bucket's capacity is the power of two above its largest fill. On words that
//! spread over the cells (random letters, P = 4, 16 and 64 parts) the buckets of one template
//! reserve 8 MiB (twice a chunk's 4 MiB of entries), and nothing is reallocated from chunk to
//! chunk. Before each chunk a reader gives back the room above `bucket_keep` (four times a
//! bucket's share of a chunk, plus 64 entries) that one of its buckets grew for an earlier chunk,
//! so that the buckets of one template reserve at most `6 * CHUNK_POSITIONS + 68 * P^2` entries
//! (24 MiB plus 544 bytes per bucket) whatever the input (LOSAT_X_DCSEEDLEAN's
//! `csr_on_pool_steps`: two sets of buckets for chunks of half the size, at most
//! `6 * CHUNK_POSITIONS + 136 * P^2` entries, and twice the bucket headers). Measured on blocks of one chunk, each a
//! different unmasked repeat of period 1 to 4 (each chunk's words on the cells of one or two
//! owners, which change from chunk to chunk): 16, 15 and 12 MiB per template at P = 4, 16 and 64,
//! where buckets that kept their room reserved 16, 32 and 48 MiB (up to 2 x P x 4 MiB). The
//! bucket headers take 64 bytes each, `64 * P^2` bytes per template (16 KiB at P = 16).
//!
//! Counter memory. The owners' private counters take, per template, one u32 per counter of the
//! table per owner, each owner's copy on whole cache lines: `P * hashsize / 512` bytes, that is
//! P x 8 KiB at word size 11 and P x 32 KiB at word size 12 (128 KiB and 512 KiB at P = 16; at
//! most 255 x 32 KiB = 8 MiB), freed at the end of the fill.
//!
//! The prefault round (phase P, scheduling only). Between phases A and B of the first chunk, when
//! that chunk has at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per page of the cells, a fourth
//! set of P jobs writes 0 to one element of each page that phase B of that chunk writes, in the
//! whole 2 MiB regions of `hashtable` and `hashtable2` that it writes in full (`Fill::prefault`,
//! which gives the argument; `select_pages`; the presence vector, at most 128 KiB, holds no whole
//! region): each page has one toucher, the owner of the block of its middle element
//! (`page_touchers`); the readers mark in phase A the page of every word's cell and
//! presence-vector word, in bits of their own per array. The round runs once, before the first
//! phase B, so each write stores 0 over the 0 the element was allocated with; the later chunks
//! run without marks or round. At 4 KiB pages it touches only pages the fill writes; with
//! transparent huge pages (`enabled=always`) a first write can map a whole 2 MiB page where the
//! loop's first read maps the huge zero page, so the round touches only whole 2 MiB regions of
//! the array that the first chunk writes in full (`select_pages`, `PREFAULT_RULE`): a huge page
//! it maps holds only pages the fill writes anyway. Why: phase B reads a
//! cell before writing it, so on a page the process has not touched yet the kernel first maps its
//! shared zero page and then, at the write, a page of the process's own (two faults, and a
//! copy-on-write that flushes the page from the TLBs of the other CPUs of the search); written
//! first, the page takes one fault and no flush (on the workstation, Linux, 1,775 cDNAs at 16
//! threads: 108,000 TLB shootdowns received by the search's CPUs without the round, 100 with the
//! round of 1a05a8a2, which touched every marked page; the timing evidence is Linux only, the
//! saving rests on Linux's zero page and copy-on-write).
//! `next_pos` is written without a read and is not touched (touching it cost more than it saved).
//! Memory: the readers' marks take, per reader, one bit per page of the cells of each template and
//! of the presence vector, on whole cache lines (576 bytes per template at word size 11, 2,112 at
//! 12; 9 KiB and 33 KiB at P = 16), and the touchers' lists 4 bytes per page; freed at the end of
//! the fill.
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
//! lock, spin or sleep), so no job can block a thread that steals it. The CSR build of
//! LOSAT_X_DCCSR (`csr`) is the exception: one `rayon::broadcast` whose jobs spin at a barrier,
//! which needs every thread of the pool free to join it (its doc gives the conditions).
//!
//! Paths that reach it: `blastn -task dc-megablast`, `-subject` and `-db`, with one template
//! (`coding` or `optimal`) or two (`coding_and_optimal`), template lengths 16, 18 and 21 and
//! word sizes 11 and 12, at two threads or more. Megablast and blastn build the contiguous
//! tables (`build_two_stage_lookup`, `build_na_lookup`, `build_pv_direct_lookup`) unless
//! `-template_type` and `-template_length` are given (NCBI applies them to any task,
//! `coordination::determine_template`); then they build this table too. The callers in
//! `x_pair_seed_disc.rs` are unit tests (on the test thread, so the original loop).
//!
//! With LOSAT_TIMING (any switch state) `main` prints `[TIMING] disc_lookup_build: <s>s
//! (builds=<n> parallel=<n on a pool>)` when at least one table was built: `<s>` is the sum of
//! the wall times of every build of the process (allocation, fill and longest chain); with
//! LOSAT_X_QBATCHPAR it includes the slots' builds, which run concurrently. A second line,
//! `[TIMING] disc_lookup_split: alloc= loop= pool_fill= chain= chunks= read_wall= read_max=
//! read_sum= read_late= update_... link_... buckets_peak= x_prefault_wall= x_prefault_max=
//! x_prefault_sum= x_prefault_late= x_prefault_rounds= x_prefault_touches=`, splits the same
//! builds: the allocation of the tables, the fill by the original loop, the fill on the pool (the
//! prefault rounds included) and the longest chain (seconds, summed); then, for the pool fills
//! only, their chunks, and per phase (read = A, update = B, link = C), summed over the chunks: the
//! calling worker's wall time from before the phase's first job to the return of its last (`wall`),
//! the longest job (`max`), the jobs' summed running time (`sum`) and the latest start of a job
//! after the phase began (`late`); the most the buckets of one fill reserved (both templates); and
//! the same four times of the prefault rounds (phase P), the number of rounds and the zeros they
//! wrote (one per page touched). The jobs' start and end are taken only with LOSAT_TIMING.

use super::x_mb_lookup_par::{plan, pool_threads, Shared, K_COMPRESSION_FACTOR_SHIFT};
use super::{
    build_unmasked_ranges, DiscTables, DiscWords, PvArrayType, BLAST2NA_MASK, PV_ARRAY_MASK,
};
#[cfg(test)]
use super::{disc_longest_chain, fill_disc_mb_tables};
use crate::algorithm::blastn::disc_lookup::compute_discontiguous_index;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastn::disc_lookup::DiscTemplateType;
#[cfg(test)]
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads LOSAT_X_DISCLUTPAR once (a switch of the strict set: on unless
/// `LOSAT_X_DISCLUTPAR=0`). It chooses between two LOSAT code paths with the same result and
/// changes no value NCBI computes.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DISCLUTPAR"))
}

/// No NCBI counterpart: reads LOSAT_X_DCSEEDLEAN once (a switch of the strict set: on unless
/// `LOSAT_X_DCSEEDLEAN=0`): the CSR of a one-template table built on a pool (`csr`) is read and
/// placed in steps over two sets of buckets (`csr_on_pool_steps`) instead of `csr_on_pool`'s
/// chunks. It chooses between two LOSAT code paths with the same result and changes no value NCBI
/// computes.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn seed_lean_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCSEEDLEAN"))
}

/// No NCBI counterpart: the number of workers of the fill when `on`, i.e. the size of the
/// current search pool when there is one of two threads or more; `None` runs the original loop.
fn parts_for(on: bool) -> Option<usize> {
    if !on {
        return None;
    }
    pool_threads()
}

/// No NCBI counterpart: dispatch test of `build_disc_mb_lookup` (the switch and a search pool of
/// two threads or more).
pub(super) fn parts() -> Option<usize> {
    let parts = parts_for(on());
    if on() {
        crate::utils::xreach::x_reach!("LOSAT_X_DISCLUTPAR", parts.is_some());
    }
    parts
}

// No NCBI counterpart: the summed wall time of the builds, their number and the number built on a
// pool, for LOSAT_TIMING (any switch state); they change no value NCBI computes.
static BUILD_NANOS: AtomicU64 = AtomicU64::new(0);
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: LOSAT_TIMING's split of the builds (any switch state), summed over the
// builds of the process, in ns: the allocation of the tables, the fill by the original loop, the
// fill on the pool, and the longest chain.
static ALLOC_NANOS: AtomicU64 = AtomicU64::new(0);
static LOOP_NANOS: AtomicU64 = AtomicU64::new(0);
static POOL_FILL_NANOS: AtomicU64 = AtomicU64::new(0);
static CHAIN_NANOS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: LOSAT_TIMING's split of one phase of the pool fills (read, update, link),
/// summed over the pool fills and their chunks, in ns: the calling worker's wall time from before
/// the phase's first job to the return of its last (`wall`), the longest job (`max`), the jobs'
/// summed running time (`sum`) and the latest start of a job after the phase began (`late`).
struct PhaseTimes {
    wall: AtomicU64,
    max: AtomicU64,
    sum: AtomicU64,
    late: AtomicU64,
}

impl PhaseTimes {
    const fn new() -> Self {
        Self {
            wall: AtomicU64::new(0),
            max: AtomicU64::new(0),
            sum: AtomicU64::new(0),
            late: AtomicU64::new(0),
        }
    }
}

// No NCBI counterpart (LOSAT_TIMING only): the phases' split (read, update, link, and the
// prefault round, `PREFAULT_PHASE`), the chunks of the pool fills, the largest bucket memory of one
// pool fill (bytes reserved after a phase A, both templates), and the prefault rounds and their
// page touches (the zeros written, `Fill::prefault`).
static PHASE_TIMES: [PhaseTimes; 4] = [const { PhaseTimes::new() }; 4];
const PHASE_NAMES: [&str; 3] = ["read", "update", "link"];
/// The index of the prefault round in `PHASE_TIMES` (printed with an `x_prefault_` prefix).
const PREFAULT_PHASE: usize = 3;
static POOL_CHUNKS: AtomicU64 = AtomicU64::new(0);
static BUCKET_PEAK_BYTES: AtomicU64 = AtomicU64::new(0);
static PREFAULT_ROUNDS: AtomicU64 = AtomicU64::new(0);
static PREFAULT_TOUCHES: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: LOSAT_TIMING is set (read once).
fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

/// No NCBI counterpart: the marks of one build's `[TIMING]` intervals (only with LOSAT_TIMING): the
/// start of the build, of the fill and of the longest chain.
pub(super) struct BuildClock {
    start: Option<std::time::Instant>,
    fill: Option<std::time::Instant>,
    chain: Option<std::time::Instant>,
}

impl BuildClock {
    /// The start of a build (before the allocation of the tables).
    pub(super) fn start() -> Self {
        Self {
            start: timing_on().then(std::time::Instant::now),
            fill: None,
            chain: None,
        }
    }

    /// The start of the fill (the original loop or the pool fill).
    pub(super) fn fill(&mut self) {
        if self.start.is_some() {
            self.fill = Some(std::time::Instant::now());
        }
    }

    /// The start of the longest chain (the end of the fill).
    pub(super) fn chain(&mut self) {
        if self.start.is_some() {
            self.chain = Some(std::time::Instant::now());
        }
    }

    /// The end of the build: adds its intervals to the `[TIMING]` counters.
    pub(super) fn record(self, on_pool: bool) {
        let (Some(start), Some(fill), Some(chain)) = (self.start, self.fill, self.chain) else {
            return;
        };
        let end = std::time::Instant::now();
        let ns = |from: std::time::Instant, to: std::time::Instant| {
            to.saturating_duration_since(from).as_nanos() as u64
        };
        BUILD_NANOS.fetch_add(ns(start, end), Ordering::Relaxed);
        ALLOC_NANOS.fetch_add(ns(start, fill), Ordering::Relaxed);
        if on_pool {
            &POOL_FILL_NANOS
        } else {
            &LOOP_NANOS
        }
        .fetch_add(ns(fill, chain), Ordering::Relaxed);
        CHAIN_NANOS.fetch_add(ns(chain, end), Ordering::Relaxed);
        BUILDS.fetch_add(1, Ordering::Relaxed);
        BUILDS_ON_POOL.fetch_add(u64::from(on_pool), Ordering::Relaxed);
    }
}

/// The `[TIMING]` lines (with LOSAT_TIMING, when a table was built), printed by `main` at exit.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if timing_on() && load(&BUILDS) > 0 {
        let s = |c: &AtomicU64| load(c) as f64 / 1e9;
        eprintln!(
            "[TIMING] disc_lookup_build: {:.3}s (builds={} parallel={})",
            s(&BUILD_NANOS),
            load(&BUILDS),
            load(&BUILDS_ON_POOL)
        );
        let mut split = format!(
            "[TIMING] disc_lookup_split: alloc={:.4}s loop={:.4}s pool_fill={:.4}s chain={:.4}s chunks={}",
            s(&ALLOC_NANOS),
            s(&LOOP_NANOS),
            s(&POOL_FILL_NANOS),
            s(&CHAIN_NANOS),
            load(&POOL_CHUNKS)
        );
        for (name, t) in PHASE_NAMES.iter().zip(&PHASE_TIMES) {
            split.push_str(&format!(
                " {name}_wall={:.4}s {name}_max={:.4}s {name}_sum={:.4}s {name}_late={:.4}s",
                s(&t.wall),
                s(&t.max),
                s(&t.sum),
                s(&t.late)
            ));
        }
        split.push_str(&format!(
            " buckets_peak={:.1}MiB",
            load(&BUCKET_PEAK_BYTES) as f64 / 1048576.0
        ));
        let t = &PHASE_TIMES[PREFAULT_PHASE];
        split.push_str(&format!(
            " x_prefault_wall={:.4}s x_prefault_max={:.4}s x_prefault_sum={:.4}s x_prefault_late={:.4}s x_prefault_rounds={} x_prefault_touches={}",
            s(&t.wall),
            s(&t.max),
            s(&t.sum),
            s(&t.late),
            load(&PREFAULT_ROUNDS),
            load(&PREFAULT_TOUCHES)
        ));
        eprintln!("{split}");
    }
}

/// The tables of one template, shared by the workers (not its counters: phase B counts into the
/// owners' private copies, `OwnerCounters`).
#[derive(Clone, Copy)]
struct SharedTemplate {
    hashtable: Shared<u32>,
    next_pos: Shared<u32>,
}

/// The tables the jobs write, shared by the workers; `helper_array` and `helper_array2` are not
/// among them (only `OwnerCounters::flush` writes those, after the last join).
#[derive(Clone, Copy)]
struct SharedTables {
    first: SharedTemplate,
    second: SharedTemplate,
    pv_array: Shared<PvArrayType>,
    pv_array_bts: usize,
}

impl SharedTables {
    fn new(
        hashtable: &mut [u32],
        next_pos: &mut [u32],
        hashtable2: &mut [u32],
        next_pos2: &mut [u32],
        pv_array: &mut [PvArrayType],
        pv_array_bts: usize,
    ) -> Self {
        Self {
            first: SharedTemplate {
                hashtable: Shared::new(hashtable),
                next_pos: Shared::new(next_pos),
            },
            second: SharedTemplate {
                hashtable: Shared::new(hashtable2),
                next_pos: Shared::new(next_pos2),
            },
            pv_array: Shared::new(pv_array),
            pv_array_bts,
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
    // presence-vector word this worker owns; and the counter is incremented in `counters`, this
    // worker's private copy of the template's counters (`OwnerCounters::of`, same index).
    /// SAFETY: the worker owns the block of `ecode` (of the template's tables and of `pv_array`),
    /// and `counters` is its own copy.
    #[inline(always)]
    unsafe fn update(self, second: bool, counters: Shared<u32>, ecode: usize, index: u32) -> u32 {
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
                counters.set(h, counters.get(h).wrapping_add(1));
            }
            t.hashtable.set(ecode, index);
            head
        }
    }
}

/// No NCBI counterpart (scheduling only): 16 counters on a 64-byte cache line of their own.
#[repr(C, align(64))]
#[derive(Clone, Copy)]
struct CounterLine([u32; 16]);

/// The counters of a `CounterLine`.
const LINE_COUNTERS: usize = 16;

/// No NCBI counterpart (scheduling only): the owners' private copies of one template's counters
/// (`helper_array` or `helper_array2`, `len` counters each): owner `o`'s copy is the `len`
/// counters from counter `o * stride`, `stride` a whole number of lines, so no two owners' copies
/// share a cache line. In phase B owner `o` increments only the counters of its own blocks in its
/// copy (`SharedTables::update`); `flush` adds them into the template's counters after the last
/// chunk. `P * stride * 4` bytes (see "Counter memory" above).
struct OwnerCounters {
    lines: Vec<CounterLine>,
    stride: usize,
    len: usize,
}

/// No NCBI counterpart (scheduling only): `OwnerCounters` shared by the jobs of phase B, each of
/// which uses only its own copy (`of`).
#[derive(Clone, Copy)]
struct SharedCounters {
    all: Shared<u32>,
    stride: usize,
    len: usize,
}

impl OwnerCounters {
    /// Zeroed copies of `len` counters for `parts` owners.
    fn new(len: usize, parts: usize) -> Self {
        let stride = len.div_ceil(LINE_COUNTERS) * LINE_COUNTERS;
        Self {
            lines: vec![CounterLine([0; LINE_COUNTERS]); parts * stride / LINE_COUNTERS],
            stride,
            len,
        }
    }

    fn shared(&mut self) -> SharedCounters {
        SharedCounters {
            all: Shared {
                ptr: self.lines.as_mut_ptr().cast::<u32>(),
                len: self.lines.len() * LINE_COUNTERS,
            },
            stride: self.stride,
            len: self.len,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:781-783
    // ```c
    //          else {
    //             helper_array[ecode1/kCompressionFactor]++;
    //          }
    // ```
    /// After the last chunk (no job running): `helper[h] = helper[h].wrapping_add(n)`, `n` the
    /// count of counter `h` in the copy of the owner of its block (`owner[h >> counter_shift]`,
    /// `counter_shift = block_shift - 11`). Every increment of counter `h` was made in that copy
    /// (the module's "Why the counters are the loop's"), and the other copies hold 0 there.
    ///
    /// Serially on the calling thread: 2,048 (word size 11) or 8,192 (12) additions per template;
    /// a job per owner would cost a pool round (spawns and wake-ups) for that, and its owners
    /// would again write neighbouring counters of `helper`.
    fn flush(&self, helper: &mut [u32], owner: &[u8], counter_shift: usize) {
        debug_assert_eq!(helper.len(), self.len);
        for (h, counter) in helper.iter_mut().enumerate() {
            let i = owner[h >> counter_shift] as usize * self.stride + h;
            *counter = counter.wrapping_add(self.lines[i / LINE_COUNTERS].0[i % LINE_COUNTERS]);
        }
    }

    /// Tests only: how many owners' copies hold a non-zero counter (before the flush). A fill that
    /// counted into one shared copy would give 1 here whatever the number of owners.
    #[cfg(test)]
    fn owners_with_counts(&self, parts: usize) -> usize {
        (0..parts)
            .filter(|&o| {
                (0..self.len).any(|h| {
                    let i = o * self.stride + h;
                    self.lines[i / LINE_COUNTERS].0[i % LINE_COUNTERS] != 0
                })
            })
            .count()
    }
}

impl SharedCounters {
    /// Owner `o`'s copy, indexed as the template's counters.
    #[inline(always)]
    fn of(self, o: usize) -> Shared<u32> {
        debug_assert!(self.all.len == 0 || (o + 1) * self.stride <= self.all.len);
        Shared {
            // In bounds of the copies (`o < parts`); nothing is dereferenced here.
            ptr: self.all.ptr.wrapping_add(o * self.stride),
            len: self.len,
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

/// No NCBI counterpart (scheduling only): one bucket's list on a 64-byte cache line of its own, so
/// that the pushes of one reader (each stores the `Vec`'s length) never write a line that holds
/// another reader's bucket header (layout only).
#[repr(align(64))]
#[derive(Default)]
struct Bucket(Vec<Entry>);

/// No NCBI counterpart (scheduling only): the buckets of one template, `parts x parts` lists of
/// words: bucket `(a, o)` holds the words that reader `a` read whose cell lies in owner `o`'s
/// blocks, in reading order. Shared by the jobs of a phase, each of which uses only its own
/// buckets: phase A's job `a` the row `(a, *)`, phase B's job `o` the column `(*, o)`, phase C's
/// job `a` the row `(a, *)`; the phases are separated by joins.
#[derive(Clone, Copy)]
struct SharedBuckets {
    ptr: *mut Bucket,
    parts: usize,
    len: usize,
}

// SAFETY: the jobs of a phase use disjoint buckets (see `SharedBuckets`), the phases are separated
// by joins, and the buckets outlive the phases.
unsafe impl Send for SharedBuckets {}
// SAFETY: as above.
unsafe impl Sync for SharedBuckets {}

impl SharedBuckets {
    fn new(buckets: &mut [Bucket], parts: usize) -> Self {
        debug_assert!(buckets.is_empty() || buckets.len() == parts * parts);
        Self {
            ptr: buckets.as_mut_ptr(),
            parts,
            len: buckets.len(),
        }
    }

    /// The entries the buckets have room for (their summed capacity).
    ///
    /// SAFETY: no job is running (between two phases).
    unsafe fn capacity(self) -> usize {
        (0..self.len)
            .map(|i| unsafe { (*self.ptr.add(i)).0.capacity() })
            .sum()
    }

    /// SAFETY: no other job of the current phase uses the bucket `(a, o)`, and the buckets outlive
    /// the returned borrow (the phase).
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    unsafe fn at<'b>(self, a: usize, o: usize) -> &'b mut Vec<Entry> {
        debug_assert!(a < self.parts && o < self.parts);
        unsafe { &mut (*self.ptr.add(a * self.parts + o)).0 }
    }
}

/// No NCBI counterpart (scheduling only): the page size the prefault round assumes, 4 KiB (the
/// smallest page of the targets; where pages are larger, some get more than one touch, which writes
/// the same 0). The argument that the round maps no page the fill does not write holds at this
/// granularity; for 2 MiB transparent huge pages see `Fill::prefault` and the size rule.
pub(super) const PAGE_BYTES: usize = 4096;

/// log2 of the 4-byte elements of a page (every array the round touches holds 4-byte elements).
pub(super) const PAGE_ELEMS_SHIFT: usize = 10;

const _: () = assert!(PAGE_BYTES == 4 << PAGE_ELEMS_SHIFT);
const _: () = assert!(std::mem::size_of::<PvArrayType>() == 4);

/// The element offset of `ptr` (4-byte elements) in its page: element `i` of an array that starts
/// there lies on the array's page `(offset + i) >> PAGE_ELEMS_SHIFT` (its first page is 0).
pub(super) fn page_offset<T>(ptr: *const T) -> usize {
    (ptr as usize % PAGE_BYTES) / 4
}

/// The pages of an array of `len` 4-byte elements from page offset `offset`.
pub(super) fn pages_of(offset: usize, len: usize) -> usize {
    if len == 0 {
        0
    } else {
        ((offset + len - 1) >> PAGE_ELEMS_SHIFT) + 1
    }
}

/// The arrays whose pages the prefault round touches (`FaultJob::pages`, `Faults::marks`, and
/// the tests' log of touches): `HASHTABLE` + template and `PV_ARRAY`.
pub(super) const HASHTABLE: usize = 0;
const PV_ARRAY: usize = 2;

/// No NCBI counterpart (scheduling only): 512 page bits on a 64-byte line of their own.
#[repr(C, align(64))]
#[derive(Clone, Copy)]
struct BitLine([u64; 8]);

/// No NCBI counterpart (scheduling only): one bit per page of one array for each reader: reader
/// `r`'s bits are the `stride` words from word `r * stride`, a whole number of lines, so no two
/// readers' bits share a line. A set bit: in phase A of the current chunk, the reader read a word
/// whose cell (or the cell's presence-vector word) lies on that page.
pub(super) struct PageBits {
    lines: Vec<BitLine>,
    stride: usize,
}

impl PageBits {
    pub(super) fn new(pages: usize, parts: usize) -> Self {
        let stride = pages.div_ceil(64).div_ceil(8) * 8;
        Self {
            lines: vec![BitLine([0; 8]); parts * stride / 8],
            stride,
        }
    }

    pub(super) fn shared(&mut self) -> SharedBits {
        SharedBits {
            ptr: self.lines.as_mut_ptr().cast::<u64>(),
            stride: self.stride,
        }
    }
}

/// No NCBI counterpart (scheduling only): `PageBits` shared by the jobs: in phase A each reader
/// writes only its own bits; in the prefault round the jobs only read them.
#[derive(Clone, Copy)]
pub(super) struct SharedBits {
    ptr: *mut u64,
    stride: usize,
}

// SAFETY: in phase A each job writes only its own reader's bits, in the prefault round they are
// only read, the phases are separated by joins, and the bits outlive the phases.
unsafe impl Send for SharedBits {}
// SAFETY: as above.
unsafe impl Sync for SharedBits {}

impl SharedBits {
    /// Clears reader `r`'s bits.
    ///
    /// SAFETY: phase A, job `r`.
    pub(super) unsafe fn clear(self, r: usize) {
        // SAFETY: reader `r`'s bits (above).
        unsafe { std::ptr::write_bytes(self.ptr.add(r * self.stride), 0, self.stride) };
    }

    /// Sets reader `r`'s bit of page `page`.
    ///
    /// SAFETY: phase A, job `r`; `page < 64 * stride`.
    #[inline(always)]
    pub(super) unsafe fn mark(self, r: usize, page: usize) {
        debug_assert!(page < 64 * self.stride);
        // SAFETY: in reader `r`'s bits (above).
        unsafe { *self.ptr.add(r * self.stride + (page >> 6)) |= 1u64 << (page & 63) };
    }

    /// Reader `r`'s bits of pages `64 * w .. 64 * w + 64`.
    ///
    /// SAFETY: no job writes the bits; `w < stride`.
    #[inline(always)]
    pub(super) unsafe fn word(self, r: usize, w: usize) -> u64 {
        debug_assert!(w < self.stride);
        // SAFETY: as above.
        unsafe { *self.ptr.add(r * self.stride + w) }
    }

    /// Whether reader `r`'s bit of page `page` is set.
    ///
    /// SAFETY: no job writes the bits (the prefault round); `page < 64 * stride`.
    #[inline(always)]
    unsafe fn marked(self, r: usize, page: usize) -> bool {
        debug_assert!(page < 64 * self.stride);
        // SAFETY: as above.
        unsafe { *self.ptr.add(r * self.stride + (page >> 6)) & (1u64 << (page & 63)) != 0 }
    }
}

/// Tests only: one write of the prefault round: the array (`HASHTABLE` + template or `PV_ARRAY`),
/// the element, and its value before and after the write.
#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(super) struct Touch {
    pub(super) array: usize,
    pub(super) index: usize,
    pub(super) before: u32,
    pub(super) after: u32,
}

/// The log of the prefault round's writes: in tests a list of `Touch`es, else nothing.
#[cfg(test)]
pub(super) type TouchLog = Vec<Touch>;
#[cfg(not(test))]
pub(super) type TouchLog = ();

/// No NCBI counterpart (scheduling only): what toucher `j` of a fill keeps for the prefault round,
/// on lines of its own (it writes them in the round; the calling thread reads them after it): its
/// pages of `hashtable`, `hashtable2` and `pv_array` (`page_touchers`) and the number of zeros it
/// wrote.
#[repr(C, align(64))]
#[derive(Default)]
struct FaultJob {
    pages: [Vec<u32>; 3],
    touches: usize,
    log: TouchLog,
}

/// No NCBI counterpart (scheduling only): the `FaultJob`s of a fill shared by the jobs, each of
/// which uses only its own; the calling thread reads them between two phases.
#[derive(Clone, Copy)]
struct SharedJobs {
    ptr: *mut FaultJob,
    len: usize,
}

// SAFETY: each job of a phase uses only its own `FaultJob`, the phases are separated by joins,
// and the `FaultJob`s outlive the phases.
unsafe impl Send for SharedJobs {}
// SAFETY: as above.
unsafe impl Sync for SharedJobs {}

impl SharedJobs {
    /// SAFETY: no other job uses `FaultJob` `j` during the borrow, which ends with the phase.
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    unsafe fn at<'b>(self, j: usize) -> &'b mut FaultJob {
        debug_assert!(j < self.len);
        unsafe { &mut *self.ptr.add(j) }
    }
}

/// No NCBI counterpart (scheduling only): the prefault's state of one fill: the page offsets
/// (`page_offset`) and lengths of `hashtable`, `hashtable2` (empty for one template) and
/// `pv_array`, the readers' marks of their pages (`PageBits`), the pages chosen to be touched
/// (`select_pages`), and the jobs' `FaultJob`s.
#[derive(Clone, Copy)]
struct Faults {
    offsets: [usize; 3],
    lens: [usize; 3],
    marks: [SharedBits; 3],
    /// The pages to touch (one row each; `select_pages`), read only in the round.
    touch: [SharedBits; 3],
    jobs: SharedJobs,
}

/// No NCBI counterpart (scheduling only): the middle element of page `page` of an array of `len`
/// 4-byte elements from page offset `offset` (`page < pages_of(offset, len)`).
pub(super) fn middle(offset: usize, len: usize, page: usize) -> usize {
    let lo = (page << PAGE_ELEMS_SHIFT).saturating_sub(offset);
    let hi = (((page + 1) << PAGE_ELEMS_SHIFT) - offset).min(len);
    debug_assert!(lo < hi);
    (lo + hi - 1) / 2
}

/// No NCBI counterpart (scheduling only): gives each page of `hashtable`, `hashtable2` and
/// `pv_array` (`Faults::offsets`, `lens`) one toucher, the owner of the block (`plan`) that holds
/// its middle element (for `pv_array`, the block of the cells of that presence-vector word), in
/// that job's `FaultJob::pages`, in rising order.
fn page_touchers(
    jobs: &mut [FaultJob],
    owner: &[u8],
    block_shift: usize,
    pv_array_bts: usize,
    offsets: [usize; 3],
    lens: [usize; 3],
) {
    for array in [HASHTABLE, HASHTABLE + 1, PV_ARRAY] {
        let shift = if array == PV_ARRAY {
            block_shift - pv_array_bts
        } else {
            block_shift
        };
        for page in 0..pages_of(offsets[array], lens[array]) {
            let m = middle(offsets[array], lens[array], page);
            // `m < len`, and the blocks cover the table (`len >> shift <= owner.len()`).
            debug_assert!((m >> shift) < owner.len());
            let j = owner[m >> shift] as usize;
            jobs[j].pages[array].push(page as u32);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:717-718
// ```c
//       mb_lt->hashtable2 = (Int4*)calloc(mb_lt->hashsize, sizeof(Int4));
// ```
/// No NCBI counterpart (scheduling only): writes 0 to element `i` of `t`, which holds 0 (the
/// caller's argument): a page touch that leaves the value as it is. A volatile write, so that it
/// is made even though it stores the value the element has. In tests it records the write with
/// the element's value before and after it (`Touch`).
///
/// SAFETY: `i < t.len`, and no other job reads or writes element `i` during the phase.
#[inline(always)]
pub(super) unsafe fn touch(t: Shared<u32>, i: usize, array: usize, log: &mut TouchLog) {
    debug_assert!(i < t.len);
    #[cfg(test)]
    // SAFETY: as the caller's.
    let before = unsafe { t.get(i) };
    // SAFETY: as the caller's.
    unsafe { std::ptr::write_volatile(t.ptr.add(i), 0) };
    #[cfg(test)]
    log.push(Touch {
        array,
        index: i,
        before,
        // SAFETY: as the caller's.
        after: unsafe { t.get(i) },
    });
    #[cfg(not(test))]
    let _ = (array, log);
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
    /// The entries a bucket keeps room for from one chunk to the next (`bucket_keep`).
    keep: usize,
    tables: SharedTables,
    first: SharedBuckets,
    second: SharedBuckets,
    /// The owners' private counters of the first and the second template.
    counters: [SharedCounters; 2],
    /// The prefault round's state (`Fill::prefault`).
    faults: Faults,
    /// Tests only: how many reader jobs ran `read::<false>` (without the marking code).
    #[cfg(test)]
    unmarked_reads: AtomicUsize,
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
    /// The job first empties its buckets, giving back the room above `keep` entries that a bucket
    /// grew in an earlier chunk (`empty`). For the prefault round (`Fill::prefault`), for each
    /// array `k` with `mark[k]` (`HASHTABLE` + template, `PV_ARRAY`), it also marks in its own bits
    /// (`Faults::marks`) the page of each word's cell and of the cell's presence-vector word (one
    /// that the update writes: `array_idx < pv_array.len()`), after clearing them. `MARK` is
    /// whether any `mark[k]` is set: the chunks without marks (every chunk but the first, and the
    /// first too when the round does not run) run a copy of the loop without the marking code (on
    /// the workstation the loop with the tests in it but no marks to make took 7-11 % longer in
    /// phase A).
    ///
    /// SAFETY: no other job of the phase uses the buckets `(a, *)` or reader `a`'s marks.
    unsafe fn read<const MARK: bool>(&self, a: usize, s0: usize, s1: usize, mark: [bool; 3]) {
        let words = self.words;
        let two_templates = words.two_templates;
        let template_length = words.template_length;
        for o in 0..self.parts {
            // SAFETY: row `a` is this job's.
            unsafe {
                empty(self.first.at(a, o), self.keep);
                if two_templates {
                    empty(self.second.at(a, o), self.keep);
                }
            }
        }
        let faults = self.faults;
        debug_assert_eq!(MARK, mark.contains(&true));
        #[cfg(test)]
        if !MARK {
            self.unmarked_reads.fetch_add(1, Ordering::Relaxed);
        }
        for (k, &on) in mark.iter().enumerate() {
            if MARK && on {
                // SAFETY: reader `a`'s marks are this job's.
                unsafe { faults.marks[k].clear(a) };
            }
        }
        let [mark1, mark2, mark_pv] = mark.map(|on| MARK && on);
        let [cells1, cells2, pv] = faults.offsets;
        let (pv_len, bts) = (faults.lens[PV_ARRAY], self.tables.pv_array_bts);
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
                // SAFETY (the marks): reader `a`'s bits; the pages are in the arrays.
                if mark1 {
                    unsafe {
                        faults.marks[0].mark(a, (cells1 + ecode1 as usize) >> PAGE_ELEMS_SHIFT)
                    };
                }
                if mark_pv && ((ecode1 as usize) >> bts) < pv_len {
                    let page = (pv + ((ecode1 as usize) >> bts)) >> PAGE_ELEMS_SHIFT;
                    unsafe { faults.marks[PV_ARRAY].mark(a, page) };
                }
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
                if mark2 {
                    unsafe {
                        faults.marks[1].mark(a, (cells2 + ecode2 as usize) >> PAGE_ELEMS_SHIFT)
                    };
                }
                if mark_pv && ((ecode2 as usize) >> bts) < pv_len {
                    let page = (pv + ((ecode2 as usize) >> bts)) >> PAGE_ELEMS_SHIFT;
                    unsafe { faults.marks[PV_ARRAY].mark(a, page) };
                }
            }
            g += 1;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:776-786
    // ```c
    //          if (mb_lt->hashtable[ecode1] == 0) {
    //             PV_SET(pv_array, ecode1, pv_array_bts);
    //          }
    // ...
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
    //          mb_lt->hashtable[ecode1] = index;
    // ```
    /// Phase P (the prefault round, between phases A and B of the first chunk only, and only when
    /// that chunk has at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per page of the cells,
    /// `fill_with`), toucher `j`: writes 0 to one element of each page chosen by `select_pages`
    /// (pages of `hashtable` or `hashtable2` that phase B of the first chunk writes, in whole
    /// 2 MiB regions of the array that it writes in full) of this job's share, so that phase B
    /// finds the page mapped. It changes no value, so phase B, and the tables, are what they are
    /// without it (no NCBI counterpart; scheduling only).
    ///
    /// Each page of the three arrays has one toucher, the owner of the block that holds its middle
    /// element (`page_touchers`). In phase A of the first chunk the readers mark (`Fill::read`) the
    /// page of every word's cell and of the cell's presence-vector word; `select_pages` keeps the
    /// marked pages of the regions it accepts; toucher `j` writes 0 to the middle element (in a
    /// block of `j`) of each of its pages so chosen. The round
    /// runs before the first phase B, when no element of the three arrays has been written yet:
    /// each still holds the 0 it was allocated with, and no other job uses it during the round.
    /// So the round writes 0 over 0. A marked page is one the first chunk's phase B writes: it
    /// writes the cell of every word, and the presence-vector word of a cell when the cell's first
    /// word comes (here, the cell is empty before the first chunk). The later chunks run without
    /// marks (`read::<false>`) and without the round; their first writes to other pages are the
    /// loop's (a read, then the write). This needs the three arrays zeroed at the start of the
    /// fill, as their callers allocate them (NCBI's `calloc`: `build_disc_mb_lookup`).
    ///
    /// Pages and max RSS: at the granularity of 4 KiB pages the round touches only pages the fill
    /// writes, so the 4 KiB pages it maps are pages the fill maps anyway. Under transparent huge
    /// pages (`enabled=always`, or `madvise` with glibc's `glibc.malloc.hugetlb=1`) a first write
    /// can map a whole 2 MiB page where the loop's first read maps the shared huge zero page and
    /// its writes then map only their 4 KiB pages: touching every marked page raised max RSS up
    /// to the whole table for small queries on the WSL machine (`PREFAULT_RULE`). So
    /// `select_pages` keeps only whole regions of 512 pages that lie within the array and that
    /// the first chunk writes in full: a huge page mapped by the round then holds only pages the
    /// fill writes anyway, at either page size. (The region's first page may hold another
    /// allocation's bytes before the array; the region is then not within the array's pages and
    /// is not touched.)
    ///
    /// Why. Phase B reads a cell (and a presence-vector word) before it writes it. On a page the
    /// process has not touched, the read maps the kernel's shared zero page and the write then
    /// replaces it with a page of its own: a second fault, and a copy-on-write that flushes the
    /// page's entry from the TLBs of the other CPUs running the search. Written first, the page
    /// takes one fault and no flush. On the workstation (Linux 6.8, THP `madvise`) one thread
    /// touching 4,096 fresh pages took 8,193 faults and 13-18 ms reading first, 4,063 faults and
    /// 7 ms writing first; in the pool fill of 1,775 cDNAs at 16 threads the CPUs of the search
    /// received 108,000 TLB shootdowns without the round and 100 with it. `next_pos` is not
    /// touched: phase C writes it without reading, so its faults take no zero page and no flush,
    /// and touching its pages in the round cost more than it saved (fill 0.89 of the fill without
    /// any round at 16 threads, against 0.79 with the cells alone; 1.01 against 0.91 at 4).
    ///
    /// SAFETY: no other job of the phase uses `FaultJob` `j`; the readers' marks are only read;
    /// the elements written are as above.
    unsafe fn prefault(&self, j: usize) {
        let faults = self.faults;
        // SAFETY: `FaultJob` `j` is this job's.
        let FaultJob {
            pages,
            touches,
            log,
        } = unsafe { faults.jobs.at(j) };
        let tables = self.tables;
        for (k, pages) in pages.iter().enumerate() {
            let array = match k {
                HASHTABLE => tables.first.hashtable,
                PV_ARRAY => tables.pv_array,
                _ => tables.second.hashtable,
            };
            let (offset, len, chosen) = (faults.offsets[k], faults.lens[k], faults.touch[k]);
            for &page in pages {
                let page = page as usize;
                // SAFETY: the bits are only read in this round; `page` is in the array.
                if unsafe { chosen.marked(0, page) } {
                    // SAFETY: an element of a page of toucher `j`, before the first phase B.
                    unsafe { touch(array, middle(offset, len, page), k, log) };
                    *touches += 1;
                }
            }
        }
    }

    /// Phase B, owner `o`: the updates of the words of its cells, reader by reader (so in reading
    /// order), first template then second, counting into its own copies of the counters; each
    /// entry keeps the cell's previous head.
    ///
    /// SAFETY: no other job of the phase uses the buckets `(*, o)`, the cells of `o`'s blocks or
    /// `o`'s copies of the counters.
    unsafe fn update(&self, o: usize) {
        let counters = self.counters[0].of(o);
        for a in 0..self.parts {
            // SAFETY: column `o`, the cells of `o`'s blocks and `o`'s counters are this job's.
            unsafe { update_all(self.tables, false, counters, self.first.at(a, o)) };
        }
        if self.words.two_templates {
            let counters = self.counters[1].of(o);
            for a in 0..self.parts {
                // SAFETY: as above.
                unsafe { update_all(self.tables, true, counters, self.second.at(a, o)) };
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

/// The updates of the entries of one bucket, in order, counting into `counters`; each entry's
/// `value` becomes the cell's previous head.
/// SAFETY: the cells of the entries lie in blocks of the calling job, and `counters` is its own
/// copy of the template's counters.
#[inline]
unsafe fn update_all(t: SharedTables, second: bool, counters: Shared<u32>, entries: &mut [Entry]) {
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
        // SAFETY: the cell is in the caller's blocks; `counters` is the caller's.
        entries[j].value = unsafe { t.update(second, counters, e.value as usize, e.index) };
    }
}

/// No NCBI counterpart (scheduling only): the positions read per chunk (each chunk is read, updated
/// and linked before the next; the buckets hold the words of one chunk).
const CHUNK_POSITIONS: usize = 1 << 19;

/// No NCBI counterpart (scheduling only): the entries a bucket keeps room for from one chunk to
/// the next: four times a bucket's share of a chunk (`chunk / parts^2` words), plus 64. On input
/// whose words spread over the cells, a bucket's capacity (the power of two above its largest
/// fill) stays below this, so nothing is given back and nothing is reallocated from chunk to
/// chunk; a bucket that once held many more words (a chunk whose words fall on the cells of one
/// owner, such as an unmasked short repeat) gives the rest back before the next chunk. So the
/// buckets of one template reserve at most `parts^2 * keep + 2 * chunk` entries at any time
/// (each bucket at most `keep`, or twice what it holds after growing in the current chunk):
/// `6 * chunk + 68 * parts^2` entries of 8 bytes, whatever the input.
fn bucket_keep(chunk: usize, parts: usize) -> usize {
    4 * chunk.div_ceil(parts * parts) + 64
}

/// No NCBI counterpart (scheduling only): a bucket emptied for the next chunk, its capacity cut
/// back to `keep` entries when it is larger (`bucket_keep`).
#[inline]
fn empty(bucket: &mut Vec<Entry>, keep: usize) {
    bucket.clear();
    if bucket.capacity() > keep {
        bucket.shrink_to(keep);
    }
}

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

/// No NCBI counterpart (LOSAT_TIMING only): when each job of the current phase started and ended,
/// in ns after the phase began.
struct JobClock {
    begun: std::time::Instant,
    start: Vec<AtomicU64>,
    end: Vec<AtomicU64>,
}

impl JobClock {
    fn new(parts: usize) -> Self {
        Self {
            begun: std::time::Instant::now(),
            start: (0..parts).map(|_| AtomicU64::new(0)).collect(),
            end: (0..parts).map(|_| AtomicU64::new(0)).collect(),
        }
    }
}

/// No NCBI counterpart (scheduling only): one phase, `run_jobs`; with a clock (LOSAT_TIMING, pool
/// fills only) each job's start and end are taken and the phase's wall time, longest job, summed
/// job time and latest start are added to `PHASE_TIMES[phase]`. The jobs are the same either way.
fn run_phase<F: Fn(usize) + Sync>(
    clock: &mut Option<JobClock>,
    phase: usize,
    parts: usize,
    in_pool: bool,
    job: F,
) {
    let Some(clock) = clock.as_mut() else {
        run_jobs(parts, in_pool, job);
        return;
    };
    clock.begun = std::time::Instant::now();
    let c = &*clock;
    let ns = |at: std::time::Instant| at.saturating_duration_since(c.begun).as_nanos() as u64;
    run_jobs(parts, in_pool, |i| {
        let start = ns(std::time::Instant::now());
        job(i);
        let end = ns(std::time::Instant::now());
        c.start[i].store(start, Ordering::Relaxed);
        c.end[i].store(end, Ordering::Relaxed);
    });
    let wall = ns(std::time::Instant::now());
    let (mut max, mut sum, mut late) = (0u64, 0u64, 0u64);
    for i in 0..parts {
        let start = c.start[i].load(Ordering::Relaxed);
        let took = c.end[i].load(Ordering::Relaxed).saturating_sub(start);
        max = max.max(took);
        sum += took;
        late = late.max(start);
    }
    let t = &PHASE_TIMES[phase];
    t.wall.fetch_add(wall, Ordering::Relaxed);
    t.max.fetch_add(max, Ordering::Relaxed);
    t.sum.fetch_add(sum, Ordering::Relaxed);
    t.late.fetch_add(late, Ordering::Relaxed);
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
/// search pool (or one after another, `in_pool == false`, for the unit tests without a pool).
pub(super) fn fill(parts: usize, words: &DiscWords<'_>, tables: DiscTables<'_>, in_pool: bool) {
    fill_with(
        parts,
        words,
        tables,
        in_pool,
        CHUNK_POSITIONS,
        PREFAULT_RULE,
        &|_| {},
    );
}

/// EXPERIMENT (LOSAT_X_DCCSR, `x_dc_csr.rs`): the CSR of a one-template table built from its
/// words on the search pool of `parts` threads instead of the chains (`fill`), with the
/// presence bits and helper counters into `pv_array` and `helper_array` (zeroed, of the fill's
/// lengths); returns it and whether the pool built it (otherwise `x_dc_csr::build` did, on this
/// thread: no pool of `parts` threads around it, or a build without threads).
///
/// The words are read by the fill's readers (`Fill::read::<false>`, which touches no table,
/// counter or prefault state; the same chunks of `CHUNK_POSITIONS` positions and the same
/// spans) into the same `parts x parts` buckets, once to count and once to place. Owner `o`
/// holds a contiguous run of the cell blocks (block `b` goes to owner `b * parts / blocks`, so
/// the owners' cells follow each other in cell order) and makes every write to its cells:
/// - it zeroes `offsets` over its cells; per chunk, after the readers, for `a` in `0..parts`,
///   for each entry of bucket `(a, o)` in order: `offsets[ecode] += 1` (the element of the entry
///   16 ahead prefetched), adding up its entries;
/// - after the last chunk, its base is the number of entries of the owners before it, and it
///   runs the sum pass (`x_dc_csr::sum_cells`) over its cells from that base, with the
///   presence-vector words and counters of its blocks (a block owns whole ones);
/// - per chunk, after the readers, for `a` in `0..parts`, for each entry of bucket `(a, o)` in
///   order: `p = offsets[ecode] - 1; offsets[ecode] = p; qs[p] = index`.
///
/// A cell gets its words in the fill's order (chunks in order, the readers' spans in stream
/// order, a bucket in reading order), so the counts, sums and runs are those of
/// `x_dc_csr::build`. The whole build is one `rayon::broadcast` (one job on each of the pool's
/// `parts` threads) whose threads pass an `x_dc_csr::SpinBarrier` after each chunk's reading and
/// after each chunk's owner work: one `rayon::scope` per phase and chunk, as `fill` runs them,
/// started its last job about 0.15 ms late each time on the workstation (0.35 s over the 37
/// tables of 31,168 cDNAs at 4 threads per slot; v4-R R16). The barrier waits for every thread
/// of the pool, so it needs them free while the table is built: `build_disc_mb_lookup` runs on
/// the search thread before its batch is scanned, with nothing else queued on its pool, and the
/// broadcast's jobs wait for nothing but each other.
///
/// EXPERIMENT (LOSAT_X_DCSEEDLEAN): with the switch on, the pool's build is `csr_on_pool_steps`
/// (the same CSR: the reading and the owners' work of the chunks in steps over two sets of
/// buckets, one barrier per step).
pub(super) fn csr(
    parts: usize,
    words: &DiscWords<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
) -> (super::x_dc_csr::Csr, bool) {
    debug_assert!(!words.two_templates);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if (2..=255).contains(&parts)
        && rayon::current_thread_index().is_some()
        && rayon::current_num_threads() == parts
    {
        // EXPERIMENT (LOSAT_X_DCSEEDLEAN): every build on a pool takes the steps.
        let csr = if seed_lean_on() {
            crate::utils::xreach::x_reach!("LOSAT_X_DCSEEDLEAN", true);
            csr_on_pool_steps(
                parts,
                words,
                cells,
                pv_array,
                pv_array_bts,
                helper_array,
                STEP_POSITIONS,
                x_offsets,
            )
        } else {
            csr_on_pool(
                parts,
                words,
                cells,
                pv_array,
                pv_array_bts,
                helper_array,
                CHUNK_POSITIONS,
                x_offsets,
            )
        };
        return (csr, true);
    }
    let _ = parts;
    let csr = super::x_dc_csr::build(
        words,
        cells,
        pv_array,
        pv_array_bts,
        helper_array,
        x_offsets,
    );
    (csr, false)
}

/// `csr` on the calling worker's pool of `parts` threads (see there), in chunks of `chunk`
/// positions (`CHUNK_POSITIONS`; the tests take others).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn csr_on_pool(
    parts: usize,
    words: &DiscWords<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    chunk: usize,
    x_offsets: Vec<u32>,
) -> super::x_dc_csr::Csr {
    use super::x_dc_csr::{assemble, prefetch, sum_cells, Poison, SpinBarrier};
    let chunk = chunk.max(1);
    let (block_shift, mut owner) = plan(cells, pv_array_bts, parts);
    debug_assert!(
        block_shift >= K_COMPRESSION_FACTOR_SHIFT && block_shift >= pv_array_bts,
        "a block must own whole counters and presence-vector words"
    );
    let blocks = owner.len();
    for (b, o) in owner.iter_mut().enumerate() {
        *o = (b * parts / blocks) as u8;
    }
    // Owner `o`'s cells `own[o].0 .. own[o].1` (empty when it has no block).
    let mut own = vec![(0usize, 0usize); parts];
    for (b, &o) in owner.iter().enumerate().rev() {
        let r = (b << block_shift)..((b + 1) << block_shift).min(cells);
        let e = &mut own[o as usize];
        *e = (r.start, if e.1 == 0 { r.end } else { e.1 });
    }
    // Once per table: the owners' cells cover the table (every element of `offsets` is written),
    // and the positions, which bound the words, fit the `u32` offsets.
    assert_eq!(own.iter().map(|&(c0, c1)| c1 - c0).sum::<usize>(), cells);
    assert!(pv_array_bts >= 5 && helper_array.len() == cells >> 11);
    let (segments, total) = segments(words);
    assert!(
        total <= u32::MAX as usize,
        "DCCSR: more positions than u32 offsets"
    );
    // The reading-only `Fill`: `read::<false>` uses only the words, the segments, the owners and
    // the buckets.
    let mut first: Vec<Bucket> = (0..parts * parts).map(|_| Bucket::default()).collect();
    let mut second: Vec<Bucket> = Vec::new();
    let mut counters = [OwnerCounters::new(0, parts), OwnerCounters::new(0, parts)];
    let mut marks = [0, 1, 2].map(|_| PageBits::new(0, parts));
    let mut touch = [0, 1, 2].map(|_| PageBits::new(0, 1));
    let mut jobs: Vec<FaultJob> = Vec::new();
    let fill = Fill {
        words,
        segments: &segments,
        owner: &owner,
        block_shift,
        parts,
        keep: bucket_keep(chunk, parts),
        tables: SharedTables::new(&mut [], &mut [], &mut [], &mut [], &mut [], pv_array_bts),
        first: SharedBuckets::new(&mut first, parts),
        second: SharedBuckets::new(&mut second, parts),
        counters: [counters[0].shared(), counters[1].shared()],
        faults: Faults {
            offsets: [0; 3],
            lens: [0; 3],
            marks: [marks[0].shared(), marks[1].shared(), marks[2].shared()],
            touch: [touch[0].shared(), touch[1].shared(), touch[2].shared()],
            jobs: SharedJobs {
                ptr: jobs.as_mut_ptr(),
                len: 0,
            },
        },
        #[cfg(test)]
        unmarked_reads: AtomicUsize::new(0),
    };
    // `offsets` is written in full by its owners (zeroed, then counted, summed and placed) and
    // `offsets[cells]` below; `qs` gets one element per word (`total` positions bound them).
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): `x_offsets` is the batch loop's kept table or an empty vector.
    let mut offsets = super::x_dc_csr::x_offsets_buffer(x_offsets, cells + 1);
    let mut qs: Vec<u32> = Vec::with_capacity(total);
    let ob = Shared {
        ptr: offsets.as_mut_ptr(),
        len: cells + 1,
    };
    let qb = Shared {
        ptr: qs.as_mut_ptr(),
        len: total,
    };
    let pv = Shared::new(pv_array);
    let helper = Shared::new(helper_array);
    let mut own_words = vec![0u64; parts];
    let own_w = Shared::new(&mut own_words);
    let chunks = total.div_ceil(chunk);
    let span = |c: usize, a: usize| {
        let start = c * chunk;
        let end = total.min(start + chunk);
        (
            start + (end - start) * a / parts,
            start + (end - start) * (a + 1) / parts,
        )
    };
    let bar = SpinBarrier::new(parts);
    let (fill, own, bar) = (&fill, &own, &bar);
    rayon::broadcast(|ctx| {
        let o = ctx.index();
        let _poison = Poison(bar);
        let (ob, qb, pv, helper, own_w) = (ob, qb, pv, helper, own_w);
        let (c0, c1) = own[o];
        // SAFETY (all of the job's accesses): reader `o` uses only the buckets `(o, *)` and owner
        // `o` only the buckets `(*, o)`, between barriers that separate the readers from the
        // owners; owner `o` alone writes `offsets[c0 .. c1]`, the presence-vector words and
        // counters of its blocks, its entry of the owners' words (read by the later owners after
        // a barrier) and `qs` at its places, which lie in `base .. base + its words` (disjoint
        // from the other owners'); every element read was written before by the same thread or
        // before a barrier.
        unsafe {
            std::ptr::write_bytes(ob.ptr.add(c0), 0, c1 - c0);
            let mut got = 0u64;
            for c in 0..chunks {
                let (s0, s1) = span(c, o);
                fill.read::<false>(o, s0, s1, [false; 3]);
                bar.wait();
                for a in 0..parts {
                    let es = fill.first.at(a, o);
                    got += es.len() as u64;
                    for j in 0..es.len() {
                        if let Some(ahead) = es.get(j + 16) {
                            prefetch(ob.ptr.wrapping_add(ahead.value as usize));
                        }
                        let cell = es[j].value as usize;
                        ob.set(cell, ob.get(cell) + 1);
                    }
                }
                if c + 1 == chunks {
                    own_w.set(o, got);
                }
                bar.wait();
            }
            let base: u64 = (0..o).map(|p| own_w.get(p)).sum();
            if c0 < c1 {
                sum_cells(
                    pv.ptr,
                    pv.len,
                    helper.ptr,
                    ob.ptr,
                    base as u32,
                    c0,
                    c1,
                    pv_array_bts,
                );
            }
            for c in 0..chunks {
                let (s0, s1) = span(c, o);
                fill.read::<false>(o, s0, s1, [false; 3]);
                bar.wait();
                for a in 0..parts {
                    for e in fill.first.at(a, o).iter() {
                        let cell = e.value as usize;
                        let p = ob.get(cell) - 1;
                        ob.set(cell, p);
                        qb.set(p as usize, e.index);
                    }
                }
                if c + 1 < chunks {
                    bar.wait();
                }
            }
        }
    });
    let n_words: u64 = own_words.iter().sum();
    assert!(
        n_words <= total as u64 && n_words <= u32::MAX as u64,
        "DCCSR: more words than positions"
    );
    // SAFETY: no job is running; `offsets[cells]` is no cell's end. Every element of `offsets`
    // was written (the owners' cells cover `0 .. cells`) and the runs partition `0 .. n_words`,
    // each place written once.
    unsafe {
        ob.set(cells, n_words as u32);
        assemble(offsets, qs, cells, n_words as usize)
    }
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCSEEDLEAN): the positions of a chunk of
/// `csr_on_pool_steps`, half of `CHUNK_POSITIONS`, so that its two sets of buckets keep room
/// (`bucket_keep`) for about as many entries as `csr_on_pool`'s one set (module comment). Chunks
/// of `CHUNK_POSITIONS` took more cycles on the workstation (v4-A15 s18 screening).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const STEP_POSITIONS: usize = CHUNK_POSITIONS / 2;

/// NCBI reference: `Fill::read`'s (blast_nalookup.c:749-767,775,793), first template only.
/// EXPERIMENT (LOSAT_X_DCSEEDLEAN): the words that end at the stream positions `[s0, s1)`, read
/// by `Fill::read`'s loop (the same ranges entered up to `template_length - 1` positions before
/// `s0`, resets, `accum`, `index` and first-template index), `f(ecode1, index)` for each, in
/// reading order. The template is the const `TT` (`DiscTemplateType::from_value`), so the index is computed
/// without the per-word dispatch of `compute_discontiguous_index`, and a range is read through
/// one slice (one bounds check per range).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[inline(always)]
fn lean_span<const TT: u8>(
    segments: &[Segment<'_>],
    template_length: usize,
    s0: usize,
    s1: usize,
    f: &mut impl FnMut(u32, u32),
) {
    let template = const { DiscTemplateType::from_value(TT as i32) };
    if s0 >= s1 {
        return;
    }
    let mut g = segments.partition_point(|seg| seg.stream + (seg.end - seg.start) <= s0);
    while g < segments.len() && segments[g].stream < s1 {
        let seg = &segments[g];
        let len = seg.end - seg.start;
        let first_pos = seg.start + (s0.max(seg.stream) - seg.stream);
        let end_pos = seg.start + (s1.min(seg.stream + len) - seg.stream);
        let warm_start = seg.start.max(first_pos.saturating_sub(template_length - 1));
        // The 1-based start of the word that ends at `warm_start + i` in the concatenated query
        // is `base + i` (`Fill::read`'s `index`; formed only for words, where it is >= 1).
        let base = (seg.query_offset + warm_start + 2).wrapping_sub(template_length);
        let first = first_pos - warm_start;
        let mut accum: u64 = 0;
        let mut valid_bases = 0usize;
        for (i, &val) in seg.seq[warm_start..end_pos].iter().enumerate() {
            if (val & BLAST2NA_MASK) != 0 {
                accum = 0;
                valid_bases = 0;
                continue;
            }
            accum = (accum << 2) | val as u64;
            valid_bases += 1;
            if valid_bases < template_length || i < first {
                continue;
            }
            f(
                compute_discontiguous_index(accum, template),
                base.wrapping_add(i) as u32,
            );
        }
        g += 1;
    }
}

/// No NCBI counterpart (LOSAT_X_DCSEEDLEAN): `lean_span` with the words' template (one instance
/// per template).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[inline(always)]
fn lean_words(
    words: &DiscWords<'_>,
    segments: &[Segment<'_>],
    s0: usize,
    s1: usize,
    mut f: impl FnMut(u32, u32),
) {
    let tl = words.template_length;
    let f = &mut f;
    match words.template_type as u8 {
        1 => lean_span::<1>(segments, tl, s0, s1, f),
        2 => lean_span::<2>(segments, tl, s0, s1, f),
        3 => lean_span::<3>(segments, tl, s0, s1, f),
        4 => lean_span::<4>(segments, tl, s0, s1, f),
        5 => lean_span::<5>(segments, tl, s0, s1, f),
        6 => lean_span::<6>(segments, tl, s0, s1, f),
        7 => lean_span::<7>(segments, tl, s0, s1, f),
        8 => lean_span::<8>(segments, tl, s0, s1, f),
        9 => lean_span::<9>(segments, tl, s0, s1, f),
        10 => lean_span::<10>(segments, tl, s0, s1, f),
        11 => lean_span::<11>(segments, tl, s0, s1, f),
        12 => lean_span::<12>(segments, tl, s0, s1, f),
        _ => lean_span::<0>(segments, tl, s0, s1, f),
    }
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCSEEDLEAN): what one thread of
/// `csr_on_pool_steps` reads with.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct StepCtx<'a, 'w> {
    t: usize,
    chunks: usize,
    /// The entries a bucket keeps room for from one chunk to the next (`bucket_keep`).
    keep: usize,
    block_shift: usize,
    words: &'a DiscWords<'w>,
    segments: &'a [Segment<'w>],
    owner: &'a [u8],
    /// Chunk `k`'s buckets are set `k % 2`.
    sets: [SharedBuckets; 2],
    bar: &'a super::x_dc_csr::SpinBarrier,
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCSEEDLEAN): `Fill::read::<false>` of a
/// one-template table through `lean_words`: reader `c.t`'s words of `[s0, s1)` into row `c.t` of
/// `set`, each in the bucket of its cell's owner, after emptying the row (`empty`).
///
/// SAFETY: row `c.t` of `set` is this reader's (no other thread uses it during the step).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[inline(never)]
unsafe fn read_lean_row(c: &StepCtx<'_, '_>, set: SharedBuckets, s0: usize, s1: usize) {
    let (t, owner, shift) = (c.t, c.owner, c.block_shift);
    for o in 0..set.parts {
        // SAFETY: row `t` is this reader's.
        unsafe { empty(set.at(t, o), c.keep) };
    }
    lean_words(c.words, c.segments, s0, s1, |cell, index| {
        debug_assert!(((cell as usize) >> shift) < owner.len());
        // SAFETY: the cell is in the table, whose blocks `owner` covers; row `t` is this reader's.
        unsafe {
            let o = *owner.get_unchecked((cell as usize) >> shift) as usize;
            set.at(t, o).push(Entry { index, value: cell });
        }
    });
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCSEEDLEAN): one pass of `csr_on_pool_steps`
/// on thread `c.t`. Step `k` (`0 ..= chunks`) reads chunk `k` (if any) into row `c.t` of set
/// `k % 2`, then applies `op` to each entry of column `c.t` of set `(k - 1) % 2` (chunk `k - 1`,
/// if any: bucket `(0, t)` to `(parts - 1, t)`, each in order), with the cell of the entry 16
/// places after it in its bucket (if any); a barrier ends each step.
///
/// SAFETY: row `t` of a set is written only by reader `t` and column `t` read only by owner `t`;
/// a set read at step `k` was written at step `k - 1` and is written again at step `k + 1`, with
/// barriers between; `op` touches only owner `t`'s cells and places.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[inline(always)]
unsafe fn steps_pass(
    c: &StepCtx<'_, '_>,
    span: &impl Fn(usize, usize) -> (usize, usize),
    mut op: impl FnMut(Entry, Option<u32>),
) {
    let (t, parts) = (c.t, c.sets[0].parts);
    for k in 0..=c.chunks {
        if k < c.chunks {
            let (s0, s1) = span(k, t);
            // SAFETY: row `t` of set `k % 2` is this reader's during step `k`.
            unsafe { read_lean_row(c, c.sets[k % 2], s0, s1) };
        }
        if k >= 1 {
            let set = c.sets[(k - 1) % 2];
            for a in 0..parts {
                // SAFETY: column `t` of set `(k - 1) % 2` is this owner's during step `k`.
                let es = unsafe { set.at(a, t) };
                for j in 0..es.len() {
                    op(es[j], es.get(j + 16).map(|e| e.value));
                }
            }
        }
        c.bar.wait();
    }
}

/// No NCBI counterpart (scheduling only): EXPERIMENT (LOSAT_X_DCSEEDLEAN) `csr` on the calling
/// worker's pool of `parts` threads, as `csr_on_pool` builds it (the same owners, buckets, count
/// pass, sum pass and place pass, so the same CSR) but in steps: two sets of buckets, chunk `k`'s
/// in set `k % 2`; step `k` of a pass reads chunk `k` with the lean loop (`read_lean_row`) and
/// then does the owner's work of chunk `k - 1`, and one barrier ends the step (`steps_pass`).
/// `csr_on_pool` waits twice per chunk (after the reading, after the owners' work), so each
/// chunk costs the slowest reading plus the slowest owner's work; here a step costs the slowest
/// thread's reading and owner's work together (about as many barriers per pass: half-size
/// chunks, one barrier each). The place pass prefetches the `offsets` element of the entry 16
/// ahead, as the count pass does. Chunks of `chunk` positions (`STEP_POSITIONS`; the tests take
/// others).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[allow(clippy::too_many_arguments)]
fn csr_on_pool_steps(
    parts: usize,
    words: &DiscWords<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    chunk: usize,
    x_offsets: Vec<u32>,
) -> super::x_dc_csr::Csr {
    use super::x_dc_csr::{assemble, prefetch, sum_cells, Poison, SpinBarrier};
    let chunk = chunk.max(1);
    let (block_shift, mut owner) = plan(cells, pv_array_bts, parts);
    debug_assert!(
        block_shift >= K_COMPRESSION_FACTOR_SHIFT && block_shift >= pv_array_bts,
        "a block must own whole counters and presence-vector words"
    );
    let blocks = owner.len();
    for (b, o) in owner.iter_mut().enumerate() {
        *o = (b * parts / blocks) as u8;
    }
    // Owner `o`'s cells `own[o].0 .. own[o].1` (empty when it has no block).
    let mut own = vec![(0usize, 0usize); parts];
    for (b, &o) in owner.iter().enumerate().rev() {
        let r = (b << block_shift)..((b + 1) << block_shift).min(cells);
        let e = &mut own[o as usize];
        *e = (r.start, if e.1 == 0 { r.end } else { e.1 });
    }
    assert_eq!(own.iter().map(|&(c0, c1)| c1 - c0).sum::<usize>(), cells);
    assert!(pv_array_bts >= 5 && helper_array.len() == cells >> 11);
    let (segments, total) = segments(words);
    assert!(
        total <= u32::MAX as usize,
        "DCCSR: more positions than u32 offsets"
    );
    let mut set0: Vec<Bucket> = (0..parts * parts).map(|_| Bucket::default()).collect();
    let mut set1: Vec<Bucket> = (0..parts * parts).map(|_| Bucket::default()).collect();
    let sets = [
        SharedBuckets::new(&mut set0, parts),
        SharedBuckets::new(&mut set1, parts),
    ];
    // As `csr_on_pool`'s (EXPERIMENT LOSAT_X_MBLUTKEEP: `x_offsets` is the batch loop's kept table
    // or an empty vector).
    let mut offsets = super::x_dc_csr::x_offsets_buffer(x_offsets, cells + 1);
    let mut qs: Vec<u32> = Vec::with_capacity(total);
    let ob = Shared {
        ptr: offsets.as_mut_ptr(),
        len: cells + 1,
    };
    let qb = Shared {
        ptr: qs.as_mut_ptr(),
        len: total,
    };
    let pv = Shared::new(pv_array);
    let helper = Shared::new(helper_array);
    let mut own_words = vec![0u64; parts];
    let own_w = Shared::new(&mut own_words);
    let chunks = total.div_ceil(chunk);
    let span = |c: usize, a: usize| {
        let start = c * chunk;
        let end = total.min(start + chunk);
        (
            start + (end - start) * a / parts,
            start + (end - start) * (a + 1) / parts,
        )
    };
    let bar = SpinBarrier::new(parts);
    let keep = bucket_keep(chunk, parts);
    let (segments, own, owner, bar) = (&segments, &own, &owner, &bar);
    rayon::broadcast(|ctx| {
        let t = ctx.index();
        let _poison = Poison(bar);
        let (ob, qb, pv, helper, own_w) = (ob, qb, pv, helper, own_w);
        let (c0, c1) = own[t];
        let c = StepCtx {
            t,
            chunks,
            keep,
            block_shift,
            words,
            segments,
            owner,
            sets,
            bar,
        };
        // SAFETY (all of the job's accesses): `steps_pass`'s for the sets; owner `t` alone writes
        // `offsets[c0 .. c1]`, the presence-vector words and counters of its blocks, its entry of
        // the owners' words (read by the later owners after a barrier) and `qs` at its places,
        // which lie in `base .. base + its words` (disjoint from the other owners'), as in
        // `csr_on_pool`.
        unsafe {
            std::ptr::write_bytes(ob.ptr.add(c0), 0, c1 - c0);
            let mut got = 0u64;
            steps_pass(&c, &span, |e, ahead| {
                if let Some(cell) = ahead {
                    prefetch(ob.ptr.wrapping_add(cell as usize));
                }
                let cell = e.value as usize;
                ob.set(cell, ob.get(cell) + 1);
                got += 1;
            });
            own_w.set(t, got);
            bar.wait();
            let base: u64 = (0..t).map(|p| own_w.get(p)).sum();
            if c0 < c1 {
                sum_cells(
                    pv.ptr,
                    pv.len,
                    helper.ptr,
                    ob.ptr,
                    base as u32,
                    c0,
                    c1,
                    pv_array_bts,
                );
            }
            steps_pass(&c, &span, |e, ahead| {
                if let Some(cell) = ahead {
                    prefetch(ob.ptr.wrapping_add(cell as usize));
                }
                let cell = e.value as usize;
                let p = ob.get(cell) - 1;
                ob.set(cell, p);
                qb.set(p as usize, e.index);
            });
        }
    });
    let n_words: u64 = own_words.iter().sum();
    assert!(
        n_words <= total as u64 && n_words <= u32::MAX as u64,
        "DCCSR: more words than positions"
    );
    // SAFETY: as `csr_on_pool`'s.
    unsafe {
        ob.set(cells, n_words as u32);
        assemble(offsets, qs, cells, n_words as usize)
    }
}

/// No NCBI counterpart (scheduling only): when the prefault round runs and which pages it touches
/// (`Fill::prefault`): only when the first chunk has at least `min_words_per_page` positions per
/// page of the cells of one template (a cost gate: below it the readers' marks rarely complete a
/// region), and only the marked pages of the regions of `2^region_shift` pages (by address) whose
/// pages in the array the first chunk all writes (`select_pages`).
#[derive(Clone, Copy, Debug)]
struct PrefaultRule {
    min_words_per_page: usize,
    region_shift: u32,
}

/// No NCBI counterpart (scheduling only): the rule of `fill`. Regions of 512 pages, the 2 MiB of a
/// transparent huge page on x86_64 and aarch64 with 4 KiB pages, so that where a first write maps
/// a whole huge page, every page of it is a page of the array that the fill writes anyway. On the
/// WSL machine (THP `always`, free 2 MiB blocks) touching every marked page instead raised max
/// RSS above #263's fill by up to the whole table for small queries (1,000 letters: 23 against 17
/// MiB at word size 11, 134 against 18 MiB at 12 with two templates) and still by 2.8 MiB at word
/// size 12 with 160 positions per page; with these regions max RSS was #263's at every size
/// measured. And at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per page of the cells in the
/// first chunk: below 8 no region was complete in any case measured (the round touched nothing),
/// and the marks and the round's own allocations are skipped there.
const PREFAULT_RULE: PrefaultRule = PrefaultRule {
    min_words_per_page: PREFAULT_MIN_WORDS_PER_PAGE,
    region_shift: 9,
};

/// See `PrefaultRule::min_words_per_page` and `PREFAULT_RULE`.
pub(super) const PREFAULT_MIN_WORDS_PER_PAGE: usize = 8;

/// No NCBI counterpart (scheduling only): after phase A of the first chunk, between two phases
/// (no job running), sets in `faults.touch[k]` the pages the round touches: for each array, the
/// pages of the regions of `2^region_shift` pages (by address: the page number of an address,
/// shifted) that lie within the array's pages and are all marked by some reader. `bases[k]` is
/// the address of the array's element 0. So where a first write maps a whole region (a 2 MiB
/// transparent huge page), every page of it is a page of the array that the first chunk writes;
/// the presence vector (128 KiB at most) holds no whole region and is not touched.
///
/// SAFETY: no job is running; the bits are the fill's.
unsafe fn select_pages(faults: &Faults, parts: usize, bases: [usize; 3], region_shift: u32) {
    for k in [HASHTABLE, HASHTABLE + 1, PV_ARRAY] {
        let pages = pages_of(faults.offsets[k], faults.lens[k]);
        if pages == 0 {
            continue;
        }
        let words = pages.div_ceil(64);
        let mut all = vec![0u64; words];
        for r in 0..parts {
            for (w, word) in all.iter_mut().enumerate() {
                // SAFETY: no job is running; `w < stride`.
                *word |= unsafe { faults.marks[k].word(r, w) };
            }
        }
        let marked = |p: usize| all[p >> 6] & (1u64 << (p & 63)) != 0;
        let base = bases[k] / PAGE_BYTES;
        let mut p0 = 0;
        while p0 < pages {
            let region = (base + p0) >> region_shift;
            let p1 = (((region + 1) << region_shift) - base).min(pages);
            // A region qualifies when all its pages are pages of the array (none beyond its first
            // or last page, which may be another allocation's) and all are marked.
            let inside =
                region << region_shift >= base && (region + 1) << region_shift <= base + pages;
            if inside && (p0..p1).all(marked) {
                for p in p0..p1 {
                    // SAFETY: no job is running; row 0 of the fill's own bits.
                    unsafe { faults.touch[k].mark(0, p) };
                }
            }
            p0 = p1;
        }
    }
}

/// No NCBI counterpart (scheduling only): what one fill did, for LOSAT_TIMING and the tests: its
/// chunks and the largest summed capacity of the buckets of both templates after a phase A (the
/// most they reserve: phases B and C add nothing), in entries; 0 when not measured. In tests also
/// `counter_owners`, the number of owners' private copies of the counters (of both templates) that
/// hold a non-zero counter before the flush; 0 otherwise.
/// Also the prefault rounds run (0 or 1) and the zeros they wrote (`touches`); in tests every one
/// of those writes (`touch_log`), the pages of `hashtable`, `hashtable2` and `pv_array` that the
/// first chunk wrote (`first_pages`) and the reader jobs that ran without the marking code
/// (`unmarked_reads`).
#[derive(Clone, Debug, Default)]
struct FillStats {
    chunks: usize,
    peak_entries: usize,
    prefault_rounds: usize,
    touches: usize,
    #[cfg(test)]
    counter_owners: usize,
    #[cfg(test)]
    touch_log: Vec<Touch>,
    #[cfg(test)]
    first_pages: [std::collections::BTreeSet<usize>; 3],
    #[cfg(test)]
    unmarked_reads: usize,
    #[cfg(test)]
    region_shift: u32,
}

/// `fill` with the chunk size and a `probe` that each phase-B job calls first, on its thread.
#[allow(clippy::too_many_arguments)]
fn fill_with<P: Fn(u8) + Sync>(
    parts: usize,
    words: &DiscWords<'_>,
    tables: DiscTables<'_>,
    in_pool: bool,
    chunk: usize,
    rule: PrefaultRule,
    probe: &P,
) -> FillStats {
    let DiscTables {
        hashtable,
        next_pos,
        helper_array,
        hashtable2,
        next_pos2,
        helper_array2,
        pv_array,
        pv_array_bts,
    } = tables;
    let table_size = hashtable.len();
    debug_assert!(!words.two_templates || hashtable2.len() == table_size);
    let (block_shift, owner) = plan(table_size, pv_array_bts, parts);
    debug_assert!(
        block_shift >= K_COMPRESSION_FACTOR_SHIFT && block_shift >= pv_array_bts,
        "a block must own whole counters and presence-vector words"
    );
    let parts = parts.clamp(1, 255);
    let chunk = chunk.max(1);
    let (segments, total) = segments(words);
    let mut first: Vec<Bucket> = (0..parts * parts).map(|_| Bucket::default()).collect();
    let mut second: Vec<Bucket> = if words.two_templates {
        (0..parts * parts).map(|_| Bucket::default()).collect()
    } else {
        Vec::new()
    };
    // The owners' private counters (empty for the second template of one).
    let mut counters = [
        OwnerCounters::new(helper_array.len(), parts),
        OwnerCounters::new(
            if words.two_templates {
                helper_array2.len()
            } else {
                0
            },
            parts,
        ),
    ];
    // The prefault's state (`Faults`): the arrays' page offsets and lengths (no `hashtable2` for
    // one template), the readers' marks, and each job's `FaultJob` with its pages to touch.
    let offsets = [
        page_offset(hashtable.as_ptr()),
        page_offset(hashtable2.as_ptr()),
        page_offset(pv_array.as_ptr()),
    ];
    let lens = [
        hashtable.len(),
        if words.two_templates {
            hashtable2.len()
        } else {
            0
        },
        pv_array.len(),
    ];
    let mut jobs: Vec<FaultJob> = (0..parts).map(|_| FaultJob::default()).collect();
    // The prefault round runs in the first chunk only, and only when that chunk has at least
    // `rule.min_words_per_page` positions per page of the cells (`PrefaultRule`); without it
    // the marks, the pages to touch and the touchers' lists are empty.
    let first_chunk = total.min(chunk.max(1));
    let round_on = first_chunk > 0
        && first_chunk
            >= rule
                .min_words_per_page
                .saturating_mul(pages_of(offsets[0], lens[0]));
    let round_pages = |k: usize| {
        if round_on {
            pages_of(offsets[k], lens[k])
        } else {
            0
        }
    };
    let mut marks = [0, 1, 2].map(|k| PageBits::new(round_pages(k), parts));
    let mut touch = [0, 1, 2].map(|k| PageBits::new(round_pages(k), 1));
    if round_on {
        page_touchers(&mut jobs, &owner, block_shift, pv_array_bts, offsets, lens);
    }
    let bases = [
        hashtable.as_ptr() as usize,
        hashtable2.as_ptr() as usize,
        pv_array.as_ptr() as usize,
    ];
    let faults = Faults {
        offsets,
        lens,
        marks: [marks[0].shared(), marks[1].shared(), marks[2].shared()],
        touch: [touch[0].shared(), touch[1].shared(), touch[2].shared()],
        jobs: SharedJobs {
            ptr: jobs.as_mut_ptr(),
            len: jobs.len(),
        },
    };
    let fill = Fill {
        words,
        segments: &segments,
        owner: &owner,
        block_shift,
        parts,
        keep: bucket_keep(chunk, parts),
        tables: SharedTables::new(
            hashtable,
            next_pos,
            hashtable2,
            next_pos2,
            pv_array,
            pv_array_bts,
        ),
        first: SharedBuckets::new(&mut first, parts),
        second: SharedBuckets::new(&mut second, parts),
        counters: [counters[0].shared(), counters[1].shared()],
        faults,
        #[cfg(test)]
        unmarked_reads: AtomicUsize::new(0),
    };
    // LOSAT_TIMING's split covers the pool fills only.
    let mut clock = (in_pool && timing_on()).then(|| JobClock::new(parts));
    let measure = clock.is_some() || cfg!(test);
    let mut stats = FillStats::default();
    #[cfg(test)]
    {
        stats.region_shift = rule.region_shift;
    }
    let mut start = 0usize;
    while start < total {
        let end = total.min(start.saturating_add(chunk));
        let span = |a: usize| {
            (
                start + (end - start) * a / parts,
                start + (end - start) * (a + 1) / parts,
            )
        };
        // In the first chunk (with `round_on`) the readers mark the pages of the arrays that have
        // any, and the prefault round runs (`Fill::prefault`); later chunks run neither.
        let mark = [0, 1, 2].map(|k| round_on && stats.chunks == 0 && lens[k] > 0);
        // SAFETY (all four phases): job `i` of a phase uses only its own buckets, cells, marks
        // and `FaultJob` (see `SharedBuckets`, `Fill::read`, `Fill::prefault`, `Fill::update`,
        // `Fill::link`); `run_jobs` returns after every job of the phase has returned.
        run_phase(&mut clock, 0, parts, in_pool, |a| {
            let (s0, s1) = span(a);
            if mark.contains(&true) {
                unsafe { fill.read::<true>(a, s0, s1, mark) }
            } else {
                unsafe { fill.read::<false>(a, s0, s1, mark) }
            }
        });
        if measure {
            // SAFETY: between two phases.
            let reserved = unsafe { fill.first.capacity() + fill.second.capacity() };
            stats.peak_entries = stats.peak_entries.max(reserved);
        }
        // The prefault round, before the first phase B (`Fill::prefault`), on the pages of the
        // complete regions (`select_pages`).
        let round = mark.contains(&true);
        if round {
            // SAFETY: between two phases.
            unsafe { select_pages(&fill.faults, parts, bases, rule.region_shift) };
            run_phase(&mut clock, PREFAULT_PHASE, parts, in_pool, |j| unsafe {
                fill.prefault(j)
            });
            stats.prefault_rounds += 1;
        }
        run_phase(&mut clock, 1, parts, in_pool, |o| {
            probe(o as u8);
            unsafe { fill.update(o) }
        });
        run_phase(&mut clock, 2, parts, in_pool, |a| unsafe { fill.link(a) });
        // Tests: the pages of the three arrays the first chunk wrote (every written value is
        // other than 0), which the round must have touched.
        #[cfg(test)]
        if stats.chunks == 0 {
            let t = fill.tables;
            // SAFETY: between two phases; the arrays outlive the fill.
            stats.first_pages = [t.first.hashtable, t.second.hashtable, t.pv_array]
                .map(|a| nonzero_pages(unsafe { std::slice::from_raw_parts(a.ptr, a.len) }));
        }
        stats.chunks += 1;
        start = end;
    }
    #[cfg(test)]
    {
        stats.unmarked_reads = fill.unmarked_reads.load(Ordering::Relaxed);
    }
    #[cfg(test)]
    {
        stats.counter_owners = counters[0].owners_with_counts(parts)
            + if words.two_templates {
                counters[1].owners_with_counts(parts)
            } else {
                0
            };
    }
    // After the last join: the owners' counts into the template's counters.
    let counter_shift = block_shift - K_COMPRESSION_FACTOR_SHIFT;
    counters[0].flush(helper_array, &owner, counter_shift);
    if words.two_templates {
        counters[1].flush(helper_array2, &owner, counter_shift);
    }
    stats.touches = jobs.iter().map(|job| job.touches).sum();
    #[cfg(test)]
    {
        stats.touch_log = jobs.iter_mut().flat_map(|job| job.log.drain(..)).collect();
    }
    if clock.is_some() {
        POOL_CHUNKS.fetch_add(stats.chunks as u64, Ordering::Relaxed);
        BUCKET_PEAK_BYTES.fetch_max(
            (stats.peak_entries * std::mem::size_of::<Entry>()) as u64,
            Ordering::Relaxed,
        );
        PREFAULT_ROUNDS.fetch_add(stats.prefault_rounds as u64, Ordering::Relaxed);
        PREFAULT_TOUCHES.fetch_add(stats.touches as u64, Ordering::Relaxed);
    }
    stats
}

/// Tests only (`FillStats::first_pages`, `check_prefault`): the pages of `v` that hold an element
/// other than 0 (a fast test per page: most pages are all zero or all written).
#[cfg(test)]
pub(super) fn nonzero_pages(v: &[u32]) -> std::collections::BTreeSet<usize> {
    static ZEROS: [u32; 1024] = [0; 1024];
    let base = v.as_ptr() as usize;
    let mut pages = std::collections::BTreeSet::new();
    let mut i = 0;
    while i < v.len() {
        let page = (base + 4 * i) / PAGE_BYTES;
        let end = (((page + 1) * PAGE_BYTES - base) / 4).min(v.len());
        if v[i..end] != ZEROS[..end - i] {
            pages.insert(page);
        }
        i = end;
    }
    pages
}

#[cfg(test)]
mod tests {
    use super::super::{build_disc_mb_lookup, compute_mb_pv_params};
    use super::*;
    use crate::algorithm::blastn::disc_lookup::{get_disc_template_type, DiscWordType};
    use crate::utils::dust::MaskedInterval;
    use std::collections::BTreeSet;

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

    /// The prefault rule of the tests' partitioned builds: a round whatever the first chunk's size,
    /// and regions of one page (every marked page is touched).
    const EVERY_PAGE: PrefaultRule = PrefaultRule {
        min_words_per_page: 0,
        region_shift: 0,
    };

    #[derive(Clone, Debug, PartialEq, Eq)]
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
        /// `fill_with` off the pool with that many parts, positions per chunk and prefault rule
        /// (the other partitioned variants use `EVERY_PAGE`).
        Rule(usize, usize, PrefaultRule),
        /// `fill` (the default chunk size), on the pool when `true`.
        Parts(usize, bool),
        /// `fill_with` with that many parts and positions per chunk, on the pool when `true`.
        Chunked(usize, usize, bool),
    }

    /// The tables of `case` allocated as `build_disc_mb_lookup` allocates them and filled `how`.
    fn build(case: &Case, combo: (usize, DiscWordType, usize), how: How) -> Built {
        build_stats(case, combo, how).0
    }

    /// `build`, and what the fill did (`How::Chunked` only; else the default).
    fn build_stats(
        case: &Case,
        combo: (usize, DiscWordType, usize),
        how: How,
    ) -> (Built, FillStats) {
        build_from(case, combo, how, None)
    }

    /// `build_stats`, with the counters set to `start(h)` (counter `h` of `helper_array`) and
    /// `start(len + h)` (of `helper_array2`) before the fill when `start` is given; else zero, as
    /// `build_disc_mb_lookup` allocates them.
    fn build_from(
        case: &Case,
        combo: (usize, DiscWordType, usize),
        how: How,
        start: Option<&dyn Fn(usize) -> u32>,
    ) -> (Built, FillStats) {
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
        if let Some(start) = start {
            let len = b.helper_array.len();
            for (h, c) in b.helper_array.iter_mut().enumerate() {
                *c = start(h);
            }
            for (h, c) in b.helper_array2.iter_mut().enumerate() {
                *c = start(len + h);
            }
        }
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
        let stats = match how {
            How::Serial => {
                fill_disc_mb_tables(&words, tables);
                FillStats::default()
            }
            // `fill` is `fill_with` in the default chunks without a probe; here the prefault round
            // runs whatever the first chunk's size and touches every marked page (`EVERY_PAGE`),
            // so that every partitioned build checks it (`fill`'s rule has tests of its own).
            How::Parts(parts, in_pool) => fill_with(
                parts,
                &words,
                tables,
                in_pool,
                CHUNK_POSITIONS,
                EVERY_PAGE,
                &|_| {},
            ),
            How::Chunked(parts, chunk, in_pool) => {
                fill_with(parts, &words, tables, in_pool, chunk, EVERY_PAGE, &|_| {})
            }
            How::Rule(parts, chunk, rule) => {
                fill_with(parts, &words, tables, false, chunk, rule, &|_| {})
            }
        };
        if !matches!(how, How::Serial) {
            check_prefault(&b, &stats);
        }
        b.longest_chain = disc_longest_chain(&b.helper_array, &b.helper_array2, two_templates);
        (b, stats)
    }

    /// The prefault round of a partitioned fill of `b` (the arrays it filled in place):
    /// `check_prefault_arrays` on `hashtable`, `hashtable2` and `pv_array`.
    fn check_prefault(b: &Built, stats: &FillStats) -> [BTreeSet<usize>; 3] {
        check_prefault_arrays([&b.hashtable, &b.hashtable2, &b.pv_array], stats)
    }

    /// The prefault round of a fill of `arrays` (`hashtable`, `hashtable2`, `pv_array`, as filled
    /// in place), from its log of writes (`FillStats::touch_log`): each wrote 0 over 0; there was
    /// at most one round; and the pages it wrote to are, array by array, exactly the pages the
    /// first chunk wrote (`FillStats::first_pages`) in the regions of `2^region_shift` pages that lie
    /// within the array's pages and that the first chunk wrote in full, when the round ran, none
    /// otherwise. So it maps no page (and no region) the fill does not write, and misses none it
    /// should touch; and those pages hold a value after the fill. Returns the pages touched per
    /// array.
    fn check_prefault_arrays(arrays: [&[u32]; 3], stats: &FillStats) -> [BTreeSet<usize>; 3] {
        let mut touched: [BTreeSet<usize>; 3] = Default::default();
        for t in &stats.touch_log {
            assert!(
                t.before == 0 && t.after == 0,
                "the prefault wrote {} over {} at element {} of array {}",
                t.after,
                t.before,
                t.index,
                t.array
            );
            let v = arrays[t.array];
            assert!(t.index < v.len(), "{t:?}");
            touched[t.array].insert((v.as_ptr() as usize + 4 * t.index) / PAGE_BYTES);
        }
        assert_eq!(stats.touches, stats.touch_log.len());
        assert!(
            stats.prefault_rounds <= 1,
            "{} rounds",
            stats.prefault_rounds
        );
        for k in 0..3 {
            // The first chunk's pages in the regions (`2^region_shift` pages by address) that lie
            // within the array's pages and that the first chunk wrote in full.
            let first = if stats.prefault_rounds == 1 && !arrays[k].is_empty() {
                let (a0, a1) = (
                    arrays[k].as_ptr() as usize / PAGE_BYTES,
                    (arrays[k].as_ptr() as usize + 4 * arrays[k].len() - 1) / PAGE_BYTES,
                );
                let shift = stats.region_shift;
                stats.first_pages[k]
                    .iter()
                    .copied()
                    .filter(|&p| {
                        let region = p >> shift;
                        let (lo, hi) = (region << shift, ((region + 1) << shift) - 1);
                        lo >= a0 && hi <= a1 && (lo..=hi).all(|q| stats.first_pages[k].contains(&q))
                    })
                    .collect()
            } else {
                BTreeSet::new()
            };
            assert!(
                touched[k] == first,
                "array {k}: touched {} pages, the first chunk wrote {}; touched but not written {:?}, written but not touched {:?}",
                touched[k].len(),
                first.len(),
                touched[k].difference(&first).take(5).collect::<Vec<_>>(),
                first.difference(&touched[k]).take(5).collect::<Vec<_>>()
            );
            assert!(touched[k].is_subset(&nonzero_pages(arrays[k])), "array {k}");
        }
        touched
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
                EVERY_PAGE,
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

    // The premise of the warm-up of phase A: every index function of the twelve templates (the
    // second template of two included) reads only the last `template_length` bases of `accum`
    // (its low `2 * template_length` bits), and the index is a cell of the table.
    #[test]
    fn template_indices_use_only_the_last_template_length_bases() {
        let mut rng = Lcg(0x1234_5678);
        let mut uses = 0;
        for word_size in [11i32, 12] {
            for template_length in [16u8, 18, 21] {
                for word_type in [
                    DiscWordType::Coding,
                    DiscWordType::Optimal,
                    DiscWordType::TwoTemplates,
                ] {
                    let first = get_disc_template_type(word_size, template_length, word_type);
                    let mut templates = vec![first];
                    if word_type == DiscWordType::TwoTemplates {
                        templates.push(first.next());
                    }
                    let mask = (1u64 << (2 * template_length as u64)) - 1;
                    for t in templates {
                        assert!(
                            t != crate::algorithm::blastn::disc_lookup::DiscTemplateType::Contiguous
                        );
                        for _ in 0..20_000 {
                            let accum = (rng.next() << 33) ^ rng.next() ^ (rng.next() << 62);
                            let index = compute_discontiguous_index(accum, t);
                            assert_eq!(
                                index,
                                compute_discontiguous_index(accum & mask, t),
                                "{t:?} template_length {template_length}"
                            );
                            assert!((index as u64) < (1u64 << (2 * word_size as u64)), "{t:?}");
                        }
                        uses += 1;
                    }
                }
            }
        }
        assert_eq!(uses, 24);
    }

    /// Random layouts: up to 12 queries with gaps of 0-2 letters between their offsets, ranges of
    /// a few letters, masks that touch or overlap, short motifs, 2 % ambiguity codes.
    fn random_layout(rng: &mut Lcg) -> Case {
        let n_queries = 1 + rng.below(12) as usize;
        let mut queries = Vec::new();
        let mut offsets = Vec::new();
        let mut masks = Vec::new();
        let mut offset = rng.below(3) as i32;
        for _ in 0..n_queries {
            let len = match rng.below(5) {
                0 => rng.below(25) as usize,
                1 => 1 + rng.below(200) as usize,
                _ => 1 + rng.below(4000) as usize,
            };
            let motif: Vec<u8> = (0..1 + rng.below(6)).map(|_| rng.below(4) as u8).collect();
            let mut seq = Vec::with_capacity(len);
            while seq.len() < len {
                let r = rng.below(100);
                if r < 2 {
                    seq.push(4 + rng.below(12) as u8);
                } else if r < 60 {
                    seq.push(motif[seq.len() % motif.len()]);
                } else {
                    seq.push(rng.below(4) as u8);
                }
            }
            let mut m = Vec::new();
            let mut cursor = 0usize;
            while cursor < len && rng.below(4) != 0 {
                let start = cursor + rng.below(30) as usize;
                let end = start + rng.below(25) as usize;
                if start >= len {
                    break;
                }
                m.push(MaskedInterval {
                    start,
                    end: end.min(len),
                });
                cursor = if rng.below(3) == 0 { start } else { end };
            }
            offsets.push(offset);
            offset += len as i32 + 1 + rng.below(3) as i32;
            queries.push(seq);
            masks.push(m);
        }
        Case {
            queries,
            offsets,
            masks,
        }
    }

    // Random layouts (gaps between query offsets, touching and overlapping masks, many short
    // ranges) partitioned in 1-20 parts with chunks of 1-5,000 positions: the original loop's
    // tables.
    #[test]
    fn random_layouts_partitioned_equal_the_loop() {
        let mut rng = Lcg(0xBEEF_0001);
        let combos = combos();
        let mut words = 0usize;
        for round in 0..200 {
            let case = random_layout(&mut rng);
            let combo = combos[rng.below(combos.len() as u64) as usize];
            let reference = build(&case, combo, How::Serial);
            words += reference.hashtable.iter().filter(|&&c| c != 0).count();
            for _ in 0..4 {
                let parts = 1 + rng.below(20) as usize;
                let chunk = 1 + rng.below(5000) as usize;
                let got = build(&case, combo, How::Chunked(parts, chunk, false));
                assert!(
                    got == reference,
                    "round {round} {combo:?} parts {parts} chunk {chunk}"
                );
            }
        }
        assert!(words > 10_000, "{words}");
    }

    // On search pools of 2, 5 and 16 threads: random layouts in 1-24 parts and chunks of 1-3,000
    // positions, then a query of 1.3 million letters and its reverse complement in the default
    // chunks (five chunk boundaries, one inside a mask), one and two templates: the original
    // loop's tables.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn pool_fill_in_default_chunks() {
        let mut rng = Lcg(0xBEEF_0002);
        let combos = combos();
        for threads in [2usize, 5, 16] {
            crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                for round in 0..30 {
                    let case = random_layout(&mut rng);
                    let combo = combos[rng.below(combos.len() as u64) as usize];
                    let reference = build(&case, combo, How::Serial);
                    let chunk = 1 + rng.below(3000) as usize;
                    let parts = 1 + rng.below(24) as usize;
                    let got = build(&case, combo, How::Chunked(parts, chunk, true));
                    assert!(
                        got == reference,
                        "threads {threads} round {round} {combo:?} parts {parts} chunk {chunk}"
                    );
                }
                let seq: Vec<u8> = (0..1_300_000)
                    .map(|i| if i % 9973 < 3 { 14 } else { rng.below(4) as u8 })
                    .collect();
                let reverse: Vec<u8> = seq
                    .iter()
                    .rev()
                    .map(|&b| if b < 4 { 3 - b } else { b })
                    .collect();
                let case = Case {
                    queries: vec![seq, reverse],
                    offsets: vec![0, 1_300_001],
                    masks: vec![
                        vec![MaskedInterval {
                            start: 524_280,
                            end: 524_300,
                        }],
                        Vec::new(),
                    ],
                };
                for combo in [
                    (11usize, DiscWordType::TwoTemplates, 21usize),
                    (12, DiscWordType::Coding, 16),
                ] {
                    let reference = build(&case, combo, How::Serial);
                    let (got, stats) =
                        build_stats(&case, combo, How::Chunked(threads, CHUNK_POSITIONS, true));
                    assert_eq!(stats.chunks, 5, "threads {threads}");
                    assert!(got == reference, "threads {threads} long {combo:?}");
                }
                Ok(())
            })
            .unwrap();
        }
    }

    /// Blocks of `chunk` letters, each a different repeat of period 1 to 4 (unmasked): each
    /// chunk's words fall on a few cells, of one owner or two, which change from block to block.
    fn repeat_blocks(chunk: usize, blocks: usize) -> Case {
        let mut seq = Vec::with_capacity(chunk * blocks);
        for b in 0..blocks {
            let period = 1 + b % 4;
            let motif: Vec<u8> = (0..period)
                .map(|k| ((b * 7 + k * 3 + b / 4) % 4) as u8)
                .collect();
            seq.extend((0..chunk).map(|i| motif[i % period]));
        }
        Case {
            queries: vec![seq],
            offsets: vec![0],
            masks: Vec::new(),
        }
    }

    // The buckets' memory has a bound that does not grow with the parts (`bucket_keep`): at most
    // `6 * chunk + 68 * parts^2` entries per template, also when the chunks' words fall on the
    // cells of one owner after another (repeat blocks, where buckets that kept their room grew to
    // `2 * parts` times a chunk); and about twice a chunk per template on random letters. In the
    // default chunks, 4 and 16 parts, one and two templates.
    #[test]
    fn bucket_memory_is_bounded() {
        let chunk = CHUNK_POSITIONS;
        let blocks = 20;
        let mut rng = Lcg(0xBEEF_0003);
        let random = Case {
            queries: vec![(0..chunk * blocks).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        };
        let repeats = repeat_blocks(chunk, blocks);
        for (name, case) in [("random", &random), ("repeats", &repeats)] {
            for word_type in [DiscWordType::Coding, DiscWordType::TwoTemplates] {
                let templates = if word_type == DiscWordType::TwoTemplates {
                    2
                } else {
                    1
                };
                for parts in [4usize, 16] {
                    let (_, stats) =
                        build_stats(case, (11, word_type, 18), How::Chunked(parts, chunk, false));
                    let bound = templates * (6 * chunk + 68 * parts * parts);
                    eprintln!(
                        "bucket peak {name} {word_type:?} parts {parts}: {} entries = {:.1} MiB (bound {:.1} MiB)",
                        stats.peak_entries,
                        stats.peak_entries as f64 * 8.0 / 1048576.0,
                        bound as f64 * 8.0 / 1048576.0
                    );
                    assert_eq!(stats.chunks, blocks);
                    assert!(
                        stats.peak_entries <= bound,
                        "{name} {word_type:?} parts {parts}"
                    );
                    if name == "random" {
                        assert!(
                            stats.peak_entries <= templates * (2 * chunk + 64 * parts * parts),
                            "{name} {word_type:?} parts {parts}: {}",
                            stats.peak_entries
                        );
                    }
                }
            }
        }
    }

    /// One query of 60,000 letters: random letters with ten runs of 600 letters, each a repeat of
    /// period 1 to 5 (unmasked), so a few cells get hundreds of words (counters far above zero)
    /// and most blocks get no repeated cell (counters that stay where they start).
    fn counter_case() -> Case {
        let mut rng = Lcg(0xC0DE_0C01);
        let mut seq: Vec<u8> = (0..60_000).map(|_| rng.below(4) as u8).collect();
        for run in 0..10usize {
            let period = 1 + run % 5;
            let motif: Vec<u8> = (0..period).map(|_| rng.below(4) as u8).collect();
            let at = 3_000 + run * 5_600;
            for (i, letter) in seq[at..at + 600].iter_mut().enumerate() {
                *letter = motif[i % period];
            }
        }
        Case {
            queries: vec![seq],
            offsets: vec![0],
            masks: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:781-783,799-801
    // ```c
    //          else {
    //             helper_array[ecode1/kCompressionFactor]++;
    //          }
    // ...
    //          else {
    //             helper_array2[ecode2/kCompressionFactor]++;
    //          }
    // ```
    // The owners' private counters, added into the counters after the last chunk, give the loop's
    // `helper_array` and `helper_array2` on search pools of 2, 3, 5 and 16 threads (the default
    // chunk and chunks of 1,000 positions, one and two templates, word sizes 11 and 12), also
    // when the counters start where the loop's increments wrap them past u32::MAX (start values
    // u32::MAX, u32::MAX - 6 .. u32::MAX - 0, 0 and others): the case has blocks whose counters
    // wrap, blocks whose counters stay zero, and blocks no increment reaches in every class.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn private_counters_give_the_loops_counters() {
        let case = counter_case();
        let start = |h: usize| -> u32 {
            match h % 4 {
                0 => u32::MAX - (h % 7) as u32,
                1 => 0,
                2 => u32::MAX,
                _ => (h as u32).wrapping_mul(2_654_435_761),
            }
        };
        let combos = [
            (11usize, DiscWordType::Coding, 18usize),
            (11, DiscWordType::TwoTemplates, 21),
            (12, DiscWordType::TwoTemplates, 16),
            (12, DiscWordType::Optimal, 18),
        ];
        // The loop's tables from these starts, and the classes of blocks they have.
        let references: Vec<Built> = combos
            .iter()
            .map(|&combo| {
                let reference = build_from(&case, combo, How::Serial, Some(&start)).0;
                let approx = case.queries[0].len();
                let (_, bts) = compute_mb_pv_params(1 << (2 * combo.0), approx, false, combo.0);
                let (block_shift, _) = plan(1 << (2 * combo.0), bts, 2);
                let per_block = 1usize << (block_shift - K_COMPRESSION_FACTOR_SHIFT);
                let arrays: Vec<&[u32]> = if combo.1 == DiscWordType::TwoTemplates {
                    vec![&reference.helper_array, &reference.helper_array2]
                } else {
                    vec![&reference.helper_array]
                };
                for (t, got) in arrays.into_iter().enumerate() {
                    let base = t * reference.helper_array.len();
                    let (mut wrapped, mut zero, mut untouched_high) = (0, 0, 0);
                    for (b, block) in got.chunks(per_block).enumerate() {
                        let starts: Vec<u32> =
                            (0..block.len()).map(|i| start(base + b * per_block + i)).collect();
                        if block.iter().zip(&starts).any(|(&v, &s)| v < s) {
                            wrapped += 1;
                        }
                        if block.iter().zip(&starts).all(|(&v, &s)| v == s) {
                            if starts.iter().all(|&s| s == 0) {
                                zero += 1;
                            } else if starts.iter().any(|&s| s >= u32::MAX - 6) {
                                untouched_high += 1;
                            }
                        }
                    }
                    assert!(
                        wrapped > 0 && zero > 0 && untouched_high > 0,
                        "{combo:?} template {t}: blocks wrapped {wrapped}, zero {zero}, untouched near the wrap {untouched_high}"
                    );
                }
                reference
            })
            .collect();
        for threads in [2usize, 3, 5, 16] {
            crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                assert_eq!(parts_for(true), Some(threads));
                for (&combo, reference) in combos.iter().zip(&references) {
                    let got = build_from(&case, combo, How::Parts(threads, true), Some(&start)).0;
                    assert!(got == *reference, "threads {threads} {combo:?}");
                    let (got, stats) = build_from(
                        &case,
                        combo,
                        How::Chunked(threads, 1000, true),
                        Some(&start),
                    );
                    assert!(
                        got == *reference,
                        "threads {threads} {combo:?} chunks of 1,000"
                    );
                    // The owners count into copies of their own (not all into one shared copy,
                    // which is just as exact but brings back the false sharing): the case has
                    // repeats in several blocks, so at least two owners (of each template) have
                    // counts before the flush.
                    let templates = 1 + usize::from(combo.1 == DiscWordType::TwoTemplates);
                    assert!(
                        stats.counter_owners >= 2 * templates,
                        "threads {threads} {combo:?}: owners with counts {} (templates {templates})",
                        stats.counter_owners
                    );
                }
                Ok(())
            })
            .unwrap();
        }
        // From zero, as `build_disc_mb_lookup` allocates them: the loop's counters, not all zero.
        for &combo in &combos {
            let reference = build(&case, combo, How::Serial);
            assert!(reference.helper_array.iter().any(|&c| c > 100), "{combo:?}");
            for parts in [2usize, 3, 5, 16] {
                assert!(build(&case, combo, How::Parts(parts, false)) == reference);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:785-786
    // ```c
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
    //          mb_lt->hashtable[ecode1] = index;
    // ```
    // The prefault round (`Fill::prefault`, here without the size rule) on search pools of 2, 5
    // and 16 threads and off a pool in 1, 3 and 64 parts, word sizes 11 and 12, one and two
    // templates: a query of 400 letters (chunks of 97 positions and the default chunk) and one of
    // 1.2 million letters (the default chunk, and chunks of 5,000 on the pool of 5). The tables are
    // the loop's, and `check_prefault` (which every partitioned build of these tests runs) finds 0
    // written over 0, one round, and the pages touched equal to the pages the first chunk wrote,
    // array by array (`next_pos` is never touched). The short query leaves most pages of the cells
    // untouched by both; with the long one in the default chunks every page of the cells is
    // written but perhaps the end pages (with the allocator's 16-byte offset the last page holds 4
    // cells), and each page touched is touched exactly once (one toucher per page). Every reader
    // job of the later chunks ran without the marking code (`read::<false>`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn prefault_touches_the_pages_the_fill_writes() {
        let mut rng = Lcg(0xFA17_0001);
        let small = Case {
            queries: vec![(0..400).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        };
        let long = Case {
            queries: vec![(0..1_200_000)
                .map(|i| if i % 7919 < 2 { 14 } else { rng.below(4) as u8 })
                .collect()],
            offsets: vec![0],
            masks: vec![vec![MaskedInterval {
                start: 600_000,
                end: 610_000,
            }]],
        };
        let combos = [
            (11usize, DiscWordType::Coding, 18usize),
            (12, DiscWordType::TwoTemplates, 21),
        ];
        let cell_pages = |b: &Built, k: usize| {
            let v = if k == 0 { &b.hashtable } else { &b.hashtable2 };
            pages_of(page_offset(v.as_ptr()), v.len())
        };
        let check = |case: &Case, combo, how: How, reference: &Built, long: bool| {
            let (parts, chunk) = match how {
                How::Parts(p, _) => (p, CHUNK_POSITIONS),
                How::Chunked(p, c, _) => (p, c),
                _ => unreachable!(),
            };
            let what = format!("{combo:?} parts {parts} chunk {chunk}");
            let (got, stats) = build_stats(case, combo, how);
            assert!(got == *reference, "{what}: tables differ");
            let touched = check_prefault(&got, &stats);
            let templates = 1 + usize::from(combo.1 == DiscWordType::TwoTemplates);
            for t in 0..templates {
                let pages = cell_pages(&got, t);
                if long {
                    let n = touched[HASHTABLE + t].len();
                    if chunk == CHUNK_POSITIONS {
                        assert!(n + 2 >= pages, "{what} template {t}: {n} of {pages} pages");
                    }
                    let writes = stats.touch_log.iter().filter(|w| w.array == HASHTABLE + t);
                    assert_eq!(writes.count(), n, "{what} template {t}: one touch per page");
                } else {
                    let n = touched[HASHTABLE + t].len();
                    assert!(
                        n > 0 && n < pages / 2,
                        "{what} template {t}: {n} of {pages} pages"
                    );
                }
            }
            assert!(!touched[PV_ARRAY].is_empty(), "{what}");
            assert_eq!(stats.prefault_rounds, 1, "{what}");
            assert!(stats.chunks > 1 || !long, "{what}");
            assert_eq!(stats.unmarked_reads, (stats.chunks - 1) * parts, "{what}");
        };
        for &combo in &combos {
            let small_reference = build(&small, combo, How::Serial);
            let long_reference = build(&long, combo, How::Serial);
            for parts in [1usize, 3, 64] {
                for chunk in [97usize, CHUNK_POSITIONS] {
                    check(
                        &small,
                        combo,
                        How::Chunked(parts, chunk, false),
                        &small_reference,
                        false,
                    );
                }
            }
            check(&long, combo, How::Parts(3, false), &long_reference, true);
            for threads in [2usize, 5, 16] {
                crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                    check(
                        &small,
                        combo,
                        How::Chunked(threads, 97, true),
                        &small_reference,
                        false,
                    );
                    check(
                        &small,
                        combo,
                        How::Parts(threads, true),
                        &small_reference,
                        false,
                    );
                    check(
                        &long,
                        combo,
                        How::Parts(threads, true),
                        &long_reference,
                        true,
                    );
                    if threads == 5 {
                        check(
                            &long,
                            combo,
                            How::Chunked(threads, 5000, true),
                            &long_reference,
                            true,
                        );
                    }
                    Ok(())
                })
                .unwrap();
            }
        }
    }

    // Each page of the cells (one and two templates, word sizes 11 and 12) and of the presence
    // vector has one toucher, the owner of the block of its middle element, for 1, 2, 3, 16, 64
    // and 255 parts and page offsets 0, 4 and 1,023 elements; the pages come in rising order.
    #[test]
    fn every_page_has_one_toucher() {
        for word_size in [11usize, 12] {
            let table_size = 1usize << (2 * word_size);
            for approx in [100usize, 100_000] {
                let (pv_size, bts) = compute_mb_pv_params(table_size, approx, false, word_size);
                for parts in [1usize, 2, 3, 16, 64, 255] {
                    let (block_shift, owner) = plan(table_size, bts, parts);
                    for offset in [0usize, 4, 1023] {
                        let offsets = [offset, (offset + 7) % 1024, (offset + 512) % 1024];
                        let lens = [table_size, table_size, pv_size];
                        let mut jobs: Vec<FaultJob> =
                            (0..parts).map(|_| FaultJob::default()).collect();
                        page_touchers(&mut jobs, &owner, block_shift, bts, offsets, lens);
                        for k in [HASHTABLE, HASHTABLE + 1, PV_ARRAY] {
                            let pages = pages_of(offsets[k], lens[k]);
                            let mut seen = vec![0u8; pages];
                            for (j, job) in jobs.iter().enumerate() {
                                let list = &job.pages[k];
                                assert!(list.windows(2).all(|w| w[0] < w[1]));
                                for &page in list {
                                    seen[page as usize] += 1;
                                    let m = middle(offsets[k], lens[k], page as usize);
                                    let cell = if k == PV_ARRAY { m << bts } else { m };
                                    assert_eq!(owner[cell >> block_shift] as usize, j);
                                }
                            }
                            assert!(
                                seen.iter().all(|&n| n == 1),
                                "word size {word_size} parts {parts} offset {offset} array {k}"
                            );
                        }
                    }
                }
            }
        }
    }

    // `fill`'s size rule (`PREFAULT_MIN_WORDS_PER_PAGE`): the round runs only when the first chunk
    // has at least `k` positions per page of the cells. Word size 11 (4,096 or 4,097 pages), one
    // and two templates, 3 parts: a query of `k * 4,097 + 100` letters gets one round, one of
    // `k * 4,096 - 100` none, and so does the longer query in chunks of `k * 4,096 - 100`; with no
    // round nothing is touched and every reader job runs without the marking code. The tables
    // are the loop's either way.
    #[test]
    fn the_round_runs_only_when_the_first_chunk_is_dense() {
        let k = PREFAULT_MIN_WORDS_PER_PAGE;
        assert!(k >= 1);
        let mut rng = Lcg(0xFA17_0002);
        let letters = |n: usize, rng: &mut Lcg| -> Case {
            Case {
                queries: vec![(0..n).map(|_| rng.below(4) as u8).collect()],
                offsets: vec![0],
                masks: Vec::new(),
            }
        };
        let dense = letters(k * 4097 + 100, &mut rng);
        let sparse = letters(k * 4096 - 100, &mut rng);
        for combo in [
            (11usize, DiscWordType::Coding, 18usize),
            (11, DiscWordType::TwoTemplates, 16),
        ] {
            for (case, chunk, rounds) in [
                (&dense, CHUNK_POSITIONS, 1usize),
                (&sparse, CHUNK_POSITIONS, 0),
                (&dense, k * 4096 - 100, 0),
            ] {
                let reference = build(case, combo, How::Serial);
                let rule = PrefaultRule {
                    min_words_per_page: k,
                    region_shift: 0,
                };
                let (got, stats) = build_stats(case, combo, How::Rule(3, chunk, rule));
                let what = format!("{combo:?} chunk {chunk}");
                assert!(got == reference, "{what}: tables differ");
                assert_eq!(stats.prefault_rounds, rounds, "{what}");
                if rounds == 0 {
                    assert_eq!(stats.touches, 0, "{what}");
                    assert_eq!(stats.unmarked_reads, stats.chunks * 3, "{what}");
                }
            }
        }
    }

    // `fill`'s regions (`PREFAULT_RULE.region_shift`, 512 pages = 2 MiB) and smaller ones: the round
    // touches the first chunk's pages only in the regions that lie within the array's pages and
    // that the first chunk wrote in full (`check_prefault`, run by `build_from`, compares the
    // exact set). 300,000 random letters write every page of a word size 11 table, so all of its
    // whole regions are touched (whole multiples of 512 pages; the presence vector, 128 KiB,
    // holds none); 400 letters complete no region, so nothing is touched; with regions of 8
    // pages, 5,000 letters (word size 11, one template) and 40,000 (12, two templates) complete
    // some regions and not others. The tables are the loop's.
    #[test]
    fn the_round_touches_whole_regions_only() {
        let mut rng = Lcg(0xFA17_0004);
        let letters = |n: usize, rng: &mut Lcg| -> Case {
            Case {
                queries: vec![(0..n).map(|_| rng.below(4) as u8).collect()],
                offsets: vec![0],
                masks: Vec::new(),
            }
        };
        let (big, small, medium) = (
            letters(300_000, &mut rng),
            letters(400, &mut rng),
            letters(5_000, &mut rng),
        );
        let rule = |shift: u32| PrefaultRule {
            min_words_per_page: 0,
            region_shift: shift,
        };
        let coding = (11usize, DiscWordType::Coding, 18usize);
        let shift = PREFAULT_RULE.region_shift;
        let (got, stats) = build_stats(&big, coding, How::Rule(3, CHUNK_POSITIONS, rule(shift)));
        assert!(got == build(&big, coding, How::Serial));
        let pages = pages_of(page_offset(got.hashtable.as_ptr()), got.hashtable.len());
        let touched = check_prefault(&got, &stats);
        let cells = touched[HASHTABLE].len();
        assert!(
            cells + 2 * 512 >= pages && cells % 512 == 0 && cells > 0,
            "{cells} of {pages}"
        );
        assert!(touched[PV_ARRAY].is_empty());
        let (got, stats) = build_stats(&small, coding, How::Rule(3, CHUNK_POSITIONS, rule(shift)));
        assert!(got == build(&small, coding, How::Serial));
        assert_eq!(stats.prefault_rounds, 1);
        assert_eq!(stats.touches, 0);
        let medium12 = letters(40_000, &mut rng);
        for (combo, medium) in [
            (coding, &medium),
            ((12, DiscWordType::TwoTemplates, 21), &medium12),
        ] {
            let reference = build(medium, combo, How::Serial);
            let (got, stats) = build_stats(medium, combo, How::Rule(3, CHUNK_POSITIONS, rule(3)));
            assert!(got == reference, "{combo:?}");
            let touched = check_prefault(&got, &stats)[HASHTABLE].len();
            let first = stats.first_pages[HASHTABLE].len();
            assert!(
                touched > 0 && touched < first,
                "{combo:?}: touched {touched} of the first chunk's {first} pages"
            );
        }
    }

    // The round at other page offsets of the arrays than the allocator's: `hashtable`,
    // `hashtable2` and `pv_array` as slices of page-aligned buffers starting 0, 4 and 1,023
    // elements into a page (one and two templates, word size 11, 3 parts, a query of 30,000
    // letters in chunks of 7,000): the tables are the loop's, and `check_prefault_arrays` finds 0
    // over 0, one round and exactly the first chunk's pages.
    #[test]
    fn the_round_at_any_page_offset() {
        let mut rng = Lcg(0xFA17_0003);
        let case = Case {
            queries: vec![(0..30_000).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        };
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        /// A zeroed buffer with room for `len` elements from `offset` elements after a page start;
        /// returns it and that start.
        fn aligned(len: usize, offset: usize) -> (Vec<u32>, usize) {
            let v = vec![0u32; len + offset + 1024];
            let skip = (PAGE_BYTES - v.as_ptr() as usize % PAGE_BYTES) % PAGE_BYTES / 4;
            (v, skip + offset)
        }
        for combo in [
            (11usize, DiscWordType::Coding, 18usize),
            (11, DiscWordType::TwoTemplates, 21),
        ] {
            let reference = build(&case, combo, How::Serial);
            let two = combo.1 == DiscWordType::TwoTemplates;
            let (_, bts) = compute_mb_pv_params(1 << 22, 30_000, false, 11);
            for offset in [0usize, 4, 1023] {
                let mut got = reference.clone();
                let (mut cells, c0) = aligned(got.hashtable.len(), offset);
                let (mut cells2, c2) = aligned(got.hashtable2.len(), (offset + 7) % 1024);
                let (mut pv, p0) = aligned(got.pv_array.len(), (offset + 512) % 1024);
                got.next_pos.fill(0);
                got.helper_array.fill(0);
                got.next_pos2.fill(0);
                got.helper_array2.fill(0);
                let (n1, n2, np) = (
                    got.hashtable.len(),
                    got.hashtable2.len(),
                    got.pv_array.len(),
                );
                let stats = fill_with(
                    3,
                    &words_of(&case, &query_refs, combo),
                    DiscTables {
                        hashtable: &mut cells[c0..c0 + n1],
                        next_pos: &mut got.next_pos,
                        helper_array: &mut got.helper_array,
                        hashtable2: &mut cells2[c2..c2 + n2],
                        next_pos2: &mut got.next_pos2,
                        helper_array2: &mut got.helper_array2,
                        pv_array: &mut pv[p0..p0 + np],
                        pv_array_bts: bts,
                    },
                    false,
                    7_000,
                    EVERY_PAGE,
                    &|_| {},
                );
                let what = format!("{combo:?} offset {offset}");
                assert_eq!(page_offset(cells[c0..].as_ptr()), offset, "{what}");
                assert_eq!(&cells[c0..c0 + n1], &reference.hashtable[..], "{what}");
                assert_eq!(&cells2[c2..c2 + n2], &reference.hashtable2[..], "{what}");
                assert_eq!(&pv[p0..p0 + np], &reference.pv_array[..], "{what}");
                assert_eq!(got.next_pos, reference.next_pos, "{what}");
                assert_eq!(got.next_pos2, reference.next_pos2, "{what}");
                got.longest_chain = disc_longest_chain(&got.helper_array, &got.helper_array2, two);
                assert_eq!(got.helper_array, reference.helper_array, "{what}");
                assert_eq!(got.longest_chain, reference.longest_chain, "{what}");
                assert_eq!(stats.prefault_rounds, 1, "{what}");
                check_prefault_arrays(
                    [&cells[c0..c0 + n1], &cells2[c2..c2 + n2], &pv[p0..p0 + np]],
                    &stats,
                );
            }
        }
    }

    // The bucket headers and the owners' copies of the counters lie on 64-byte lines of their own:
    // a `Bucket` is one line (size and alignment 64), so the P x P headers of the grid are 64
    // bytes apart on line boundaries; a `CounterLine` is one line, and each owner's copy starts on
    // a line boundary and spans whole lines (also for a number of counters that is no multiple of
    // 16).
    #[test]
    fn bucket_headers_and_counter_copies_have_lines_of_their_own() {
        assert_eq!(std::mem::align_of::<Bucket>(), 64);
        assert_eq!(std::mem::size_of::<Bucket>(), 64);
        assert_eq!(std::mem::align_of::<CounterLine>(), 64);
        assert_eq!(std::mem::size_of::<CounterLine>(), 64);
        // The prefault's per-job state and the readers' page marks too.
        assert_eq!(std::mem::align_of::<FaultJob>(), 64);
        assert_eq!(std::mem::size_of::<FaultJob>() % 64, 0);
        assert_eq!(std::mem::align_of::<BitLine>(), 64);
        assert_eq!(std::mem::size_of::<BitLine>(), 64);
        for (pages, parts) in [(1usize, 1usize), (4097, 16), (16385, 3), (33, 255)] {
            let mut marks = PageBits::new(pages, parts);
            let shared = marks.shared();
            assert!(shared.stride % 8 == 0 && shared.stride * 64 >= pages);
            for r in 0..parts {
                assert_eq!(shared.ptr.wrapping_add(r * shared.stride) as usize % 64, 0);
            }
        }
        for parts in [1usize, 2, 3, 4, 16] {
            let grid: Vec<Bucket> = (0..parts * parts).map(|_| Bucket::default()).collect();
            for (i, bucket) in grid.iter().enumerate() {
                let at = bucket as *const Bucket as usize;
                assert_eq!(at % 64, 0, "parts {parts} bucket {i}");
                assert_eq!(
                    at - grid.as_ptr() as usize,
                    64 * i,
                    "parts {parts} bucket {i}"
                );
            }
            for len in [2048usize, 8192, 100, 1] {
                let mut counters = OwnerCounters::new(len, parts);
                let shared = counters.shared();
                assert_eq!(shared.stride % 16, 0);
                assert!(shared.stride >= len && shared.stride < len + 16);
                for o in 0..parts {
                    let copy = shared.of(o);
                    assert_eq!(
                        copy.ptr as usize % 64,
                        0,
                        "parts {parts} len {len} owner {o}"
                    );
                    assert_eq!(copy.len, len);
                    assert!(
                        copy.ptr as usize + 4 * shared.stride
                            <= shared.all.ptr as usize + 4 * shared.all.len
                    );
                }
            }
        }
    }

    // EXPERIMENT (LOSAT_X_DCCSR): the CSR built on the pool (`csr_on_pool`, and LOSAT_X_DCSEEDLEAN's
    // `csr_on_pool_steps`), in chunks of any size, is the serial build's and lists the chains of
    // the original loop (`x_dc_csr::x_assert_runs_are_the_chains` compares every cell's run with
    // its chain, the presence vector, the counters, the longest chain and the two builds; it
    // panics on a difference).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn check_csr(case: &Case, combo: (usize, DiscWordType, usize), parts: usize, chunk: usize) {
        let reference = build(case, combo, How::Serial);
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let words = words_of(case, &query_refs, combo);
        let cells = reference.hashtable.len();
        let approx: usize = case.queries.iter().map(Vec::len).sum();
        let (_, pv_array_bts) = compute_mb_pv_params(cells, approx, false, combo.0);
        // EXPERIMENT (LOSAT_X_DCSEEDLEAN): `csr_on_pool` and `csr_on_pool_steps`.
        for steps in [false, true] {
            let mut pv = vec![0; reference.pv_array.len()];
            let mut helper = vec![0; reference.helper_array.len()];
            let build = if steps {
                csr_on_pool_steps
            } else {
                csr_on_pool
            };
            let csr = build(
                parts,
                &words,
                cells,
                &mut pv,
                pv_array_bts,
                &mut helper,
                chunk,
                Vec::new(),
            );
            super::super::x_dc_csr::x_assert_runs_are_the_chains(
                &words,
                &csr,
                true,
                &pv,
                &helper,
                pv_array_bts,
                &reference.hashtable,
                &reference.next_pos,
                &reference.pv_array,
                &reference.helper_array,
                reference.longest_chain,
            );
        }
    }

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn csr_on_a_pool_lists_the_chains() {
        let mut rng = Lcg(0xC5A0_0001);
        let one_template: Vec<_> = combos()
            .into_iter()
            .filter(|c| c.1 != DiscWordType::TwoTemplates)
            .collect();
        for threads in [2usize, 3, 5] {
            crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                for round in 0..16 {
                    let case = if round % 2 == 0 {
                        random_case(&mut rng)
                    } else {
                        random_layout(&mut rng)
                    };
                    let combo = one_template[rng.below(one_template.len() as u64) as usize];
                    let chunk = [1usize, 97, 5000, CHUNK_POSITIONS][round % 4];
                    check_csr(&case, combo, threads, chunk);
                }
                // `csr` takes the pool there, and builds on this thread off it.
                let case = random_case(&mut rng);
                let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
                let words = words_of(&case, &query_refs, (11, DiscWordType::Coding, 16));
                let mut pv = vec![0; 1 << 15];
                let mut helper = vec![0; (1 << 22) / 2048];
                let (_, on_pool) = csr(
                    threads,
                    &words,
                    1 << 22,
                    &mut pv,
                    7,
                    &mut helper,
                    Vec::new(),
                );
                assert!(on_pool, "threads {threads}");
                Ok(())
            })
            .unwrap();
        }
        let case = random_case(&mut rng);
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let words = words_of(&case, &query_refs, (11, DiscWordType::Coding, 16));
        let mut pv = vec![0; 1 << 15];
        let mut helper = vec![0; (1 << 22) / 2048];
        let (_, on_pool) = csr(4, &words, 1 << 22, &mut pv, 7, &mut helper, Vec::new());
        assert!(!on_pool, "off a pool");
    }

    // EXPERIMENT (LOSAT_X_DCCSR): a table of more positions than a quarter of its cells, built on
    // a pool, is laid out as CSR (with the switch on) and answers `count_hits`, `contains_hit`
    // and `longest_chain` as the chains of the original loop; its runs are the chains.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn a_large_table_on_a_pool_answers_from_its_csr() {
        let mut rng = Lcg(0xC5A0_0002);
        let seq: Vec<u8> = (0..1_200_000)
            .map(|i| {
                if i % 7919 < 2 {
                    14
                } else if i % 3 == 0 {
                    (i / 5 % 4) as u8
                } else {
                    rng.below(4) as u8
                }
            })
            .collect();
        let case = Case {
            queries: vec![seq.clone(), seq[..300_000].to_vec()],
            offsets: vec![0, 1_200_001],
            masks: vec![
                vec![MaskedInterval {
                    start: 524_280,
                    end: 524_300,
                }],
                Vec::new(),
            ],
        };
        let combo = (11usize, DiscWordType::Coding, 16usize);
        let reference = build(&case, combo, How::Serial);
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let approx: usize = case.queries.iter().map(Vec::len).sum();
        crate::utils::threading::with_search_pool(4, "blastn", |_pool| {
            let lookup = build_disc_mb_lookup(
                &query_refs,
                &case.offsets,
                combo.0,
                combo.2,
                combo.1,
                &case.masks,
                approx,
            );
            let disc = lookup.disc().expect("a discontiguous table");
            assert_eq!(disc.x_csr().is_some(), super::super::x_dc_csr::on());
            assert_eq!(lookup.longest_chain(), reference.longest_chain);
            let mut chain = Vec::new();
            for cell in (0..reference.hashtable.len()).step_by(3) {
                chain.clear();
                let mut q = reference.hashtable[cell];
                while q != 0 {
                    chain.push(q);
                    q = reference.next_pos[q as usize];
                }
                assert_eq!(lookup.count_hits(cell as u64), chain.len(), "cell {cell}");
                for &q in chain.iter().take(3) {
                    assert!(lookup.contains_hit(cell as u64, q), "cell {cell} {q}");
                }
                assert!(!lookup.contains_hit(cell as u64, 1_500_002), "cell {cell}");
                if let Some(csr) = disc.x_csr() {
                    assert_eq!(csr.run(cell), chain.as_slice(), "cell {cell}");
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
