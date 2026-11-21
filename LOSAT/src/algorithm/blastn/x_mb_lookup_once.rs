//! EXPERIMENT (LOSAT_X_MBLUTONCE; in the strict set since 2026-10-18: on unless it is `0`)
//! : PAIRPAR's pool fill of the megablast lookup table of
//! `build_mb_lookup` (`x_mb_lookup_par.rs`) with each query position read once. In the fill of
//! `x_mb_lookup_par` every worker reads every word of the query and applies those of its own
//! cells, so each worker pays the whole read (0.65 s of each worker's time at 16 threads for
//! megablast of 31,168 D. melanogaster cDNAs on the workstation). Here the words are read by spans
//! and handed to the owners of their cells, as LOSAT_X_DISCLUTPAR does for the discontiguous
//! tables (`x_disc_lookup_par.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1033,1043-1045,1081-1091
//!
//! The plan. The cells are cut into blocks of `2^block_shift` cells dealt round robin to P owners
//! (`x_mb_lookup_par::plan`, the same plan as the fill there); P is the size of the search pool.
//! The positions the loop reads, all positions of the unmasked ranges of all queries in the loop's
//! order (`segments`: the loop's `build_unmasked_ranges` per query, without the queries shorter
//! than `lut_word_length` and the ranges shorter than `word_length`, which the loop skips), form
//! one stream, cut into chunks of `chunk_positions(P)` positions; each chunk is cut into P spans of
//! equal length in stream order, one per reader. Per chunk:
//! - A (read, job = reader `a`): reader `a` reads the words that end in its span with the loop's
//!   statements (the words of `x_range_words`, with the same ambiguity resets and k-mer, and
//!   query offset `q_off_1`), entering a range up to `lut_word_length - 1` positions before its
//!   span (never before the range's start). It writes each word `(q_off_1, cell)` to its scratch
//!   region in reading order, counting the words per owner, and then copies them, stably, into its
//!   sorted region grouped by owner: sub-range `(a, o)` holds the words of owner `o`'s cells in
//!   reading order, and reader `a`'s row of the arena holds the sub-ranges' bounds.
//! - B (update, job = owner `o`): owner `o` takes the sub-ranges `(0, o)`, `(1, o)`, ...,
//!   `(P-1, o)` and makes the loop's update for each word: `head = hashtable[cell]`; the presence
//!   bit in its private copy of `pv_array` when `head == 0`, otherwise a saturating `+1` of the
//!   cell's counter in its private copy of `helper_array` (the copies of
//!   `x_mb_lookup_par::Copies`, zeroed by the owner in its first job); then `hashtable[cell] =
//!   q_off_1`. It keeps `head` in the entry. On x86_64 it prefetches the cell 16 entries ahead.
//! - C (link, job = reader `a`): reader `a` writes `next_pos[q_off_1] = head` for the entries of
//!   its sorted region.
//!
//! C of one chunk and A of the next run in the same job (reader `a` links its region, then reads
//! the next chunk into it), so a chunk takes two `rayon::scope`s; the last scope links the last
//! chunk and has job `o` add owner `o`'s copies into `helper_array` and `pv_array` for its own
//! blocks (`x_mb_lookup_par::flush_private`). `longest_chain` is then computed by the original
//! code.
//!
//! Why the tables are the loop's. The words of a span are the loop's words ending there, with the
//! loop's cell and offset: the k-mer of a position depends only on the last `lut_word_length` bases
//! (the mask), and a word is taken iff its last `lut_word_length` positions lie in the range and
//! are unambiguous; a reader that starts `lut_word_length - 1` positions before its span (or at the
//! range's start) sees those positions, so `x_range_words` gives it a word at the same positions
//! of the span as the loop's. Each word is read by the one reader whose span holds its
//! last position. A block holds whole counters and whole presence-vector words (`2^block_shift` is
//! a multiple of 2048 and of `2^pv_array_bts`), so owner `o` is the only job that touches the cells
//! of its blocks, their counters and their presence-vector words. B takes the chunks in order, the
//! readers in order (their spans follow each other in the stream) and each sub-range in reading
//! order, and the copy into the sorted region is stable, so every cell receives the loop's words in
//! the loop's order. That fixes the cell's final head, each word's head before its update (exactly
//! the value the loop writes to `next_pos[q_off_1]`, which C writes once per word: the offsets are
//! unique), and the increments of each counter and the bits of each presence-vector word, which the
//! copies hold starting from 0; the flush gives `min(s + n, u32::MAX)` and `s | bits` from the
//! start values `s`, the loop's saturating increments and ORed bits (0 in every build;
//! `x_mb_lookup_par`'s "Why the counters and presence bits are the loop's"). So `hashtable`,
//! `next_pos`, `pv_array` and `helper_array`, and with them `longest_chain`, equal the loop's. Each
//! job of a phase touches only its own elements (the table above in the design note
//! `DESIGN_mbreadonce.md`: cells and copies by owner, `next_pos` and the regions by reader, the
//! kept heads of sub-range `(a, o)` by owner `o`); the phases are separated by the joins of
//! `rayon::scope`, and no job takes a lock, spins, sleeps or waits. No floating-point operation is
//! involved, and the number of workers changes the running time only.
//!
//! `ascending_cells` (the blastn task): the loop builds that table by applying the same update to
//! the words in reverse reading order (its re-threading, blast_lookup.c:74-76), so the chunks run
//! last to first, and B takes the readers `P-1 .. 0` and each sub-range back to front (with a
//! backward prefetch); A and C are unchanged. Each cell then receives its words in reverse reading
//! order, the order of the loop's re-threading, so the chains and heads are the loop's; the
//! counters and presence bits do not depend on the direction (a cell with n words gives its
//! counter n - 1 saturating increments and sets its presence bit either way). No per-word list is
//! kept for the re-threading.
//!
//! Memory. One arena per fill, of the calling thread, not zeroed (every element is written before
//! it is read): P reader regions of `2 x ceil(ceil(c / P) / 8)` lines of 64 bytes (scratch and
//! sorted, 8-byte entries), with `c` the smaller of the chunk (`chunk_positions`: `min(2^19, P <<
//! 15)` positions) and the positions the loop reads, and P rows of `ceil((P + 1) / 16)` lines (the
//! sub-range bounds), each region and row on lines of its own: at most `16 x (c + 8 P)` bytes of
//! regions, about P x 0.5 MiB up to 16 threads (2.0 MiB at 4, 8.0 MiB at 16) and 8.0 MiB at 16 to
//! 64 threads (8.3 MiB with the rows at 255) for a query of 2^19 positions or more; it does not
//! grow. That is per fill, not per process: fills that run at the same time each hold their own
//! arena. Until 2026-10-18 the chunk was 2^19 positions at every P, and the k query-batch slots of
//! LOSAT_X_QBATCHPAR (removed on 2026-10-28) held about k x 8 MiB: megablast of three E. coli
//! genomes (15.9 Mb) against LC738868 at 16 threads, 3 slots of 5 threads, max RSS +25 MiB against
//! the loop's fill (v4-D session 4, ws window w1). On Linux with
//! glibc its capacity is at least 32 MiB + 4 KiB (`x_mb_lookup_par::lines_capacity`, the rule of
//! `x_frame_capacity` in `tblastn/x_subject_translate.rs`), so glibc maps it and unmaps it when it
//! is freed at the end of the fill without moving its mmap threshold, and only the touched part is
//! resident; elsewhere (and on wasm32, whose linear memory never shrinks) it has its exact size and
//! is freed at the end of the fill. The ranges of the stream (`segments`) are a `Vec` of the
//! calling thread with one `Segment` (48 bytes on 64-bit targets) per unmasked range read (0.1-0.3
//! MB for the megablast tables of 31,168 cDNAs; a query masked into tens of thousands of ranges
//! takes more), freed at the end of the fill. The private copies are those of
//! `x_mb_lookup_par::Copies` (kept from fill to fill on glibc), one set per worker of each fill.
//! The per-range list of the loop's re-threading of `ascending_cells` (8 bytes a word) is not
//! needed.
//!
//! The pool. As PAIRPAR's fill (`x_mb_lookup_par::parts`): a LOSAT search pool of two
//! threads or more, no database word counts, no `BLEMIR_DEBUG`. Without the `parallel`
//! feature, or on wasm32 without `wasm-threads`, the loop runs; with `wasm-threads` the scopes run
//! on the search pool, which this module never builds (memory losat-wasm-threads-main-only).
//!
//! Timing: with LOSAT_TIMING `[TIMING] x_mb_lookup_build` counts `once=<n>` (the builds this fill
//! made) and `[TIMING] x_mb_lookup_once: fills= chunks= read_wall= read_max= read_sum= read_late=
//! update_... last_... prefault_... arena= prefault_rounds= prefault_touches=` gives, per phase
//! summed over the chunks of the search's fills: the calling worker's wall time from before
//! the phase's first job to the return of its last, the longest job, the jobs' summed time and the
//! latest start of a job; `read` is the scope of A (with C of the chunk before), `update` that of B,
//! `last` the scope that links the last chunk and flushes the copies and `prefault` the prefault
//! round (below); `arena` the largest arena of one fill (MiB); `prefault_rounds` the fills that ran
//! the round and `prefault_touches` the zeros those rounds wrote.
//!
//! EXPERIMENT (LOSAT_X_MBLUTPREFAULT; in the strict set since 2026-10-19: on unless it is `0`
//! , which turns this fill off anyway): the prefault round
//! (scheduling only), the form of LOSAT_X_DISCLUTPAR's round (`x_disc_lookup_par.rs`, "The
//! prefault round") for the one array here. Between phases A and B of the first chunk (in
//! processing order), when that chunk has at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per
//! 1,024 cells (4 KiB) of `hashtable` (`prefault_gate`), a set of P jobs writes 0 to one element of
//! each page of `hashtable` that
//! holds the cell of a word of that chunk: the readers mark in phase A the page of every word's
//! cell, in bits of their own (`PageBits`); the calling thread ORs them; toucher `j` writes 0 to
//! the middle element of each marked page whose middle element lies in a block of owner `j` (one
//! job per page). Phase B writes the cell of every word, so a marked page is a page the fill
//! writes; and before the first phase B every element of `hashtable` still holds the 0 it was
//! allocated with (NCBI's `calloc`: `build_mb_lookup` allocates the table zeroed),
//! so the round writes 0 over 0 and every array is what it is without it. Why: phase B reads a cell
//! before it writes it (`head = hashtable[cell]`), so on a page the process has not touched yet
//! the read maps the kernel's shared zero page and the write then replaces it with a page of the
//! process's own: a second fault, and a copy-on-write that flushes the page's entry from the TLBs
//! of every other CPU that runs the search. Written first, the page takes one fault and no flush.
//! The saving rests on Linux's zero page and copy-on-write, and the timing evidence is Linux only;
//! where memory comes zeroed and mapped (wasm32's linear memory) the round only costs its writes.
//! On the workstation (Linux 6.8, THP `madvise`; line R's R5, `v4-R/REPORT_R5.md` §2-5) the 7
//! tables of megablast of 65 records of three E. coli genomes against LC738868 at 64 threads made
//! 102,000 such copies, which sent 6.4 million TLB shootdowns (to 63 CPUs each); a probe round that
//! wrote every page of the table before the fill removed them and took whole runs to 0.73-0.94 of
//! the read-once fill's time; this round took whole runs to 0.61-0.93 of the time without it at
//! 4-32 threads on an i9 (WSL2, 2026-10-19; `docs/agents/diagnostics.md`). `next_pos` is not touched: phase C writes it without reading it (no
//! zero page, no flush), and touching it did not pay there (0.97-1.04). Pages and max RSS: at 4 KiB
//! pages the round touches only pages the fill writes, so the pages it maps are pages the fill maps
//! anyway. It does not keep to whole 2 MiB regions (LOSAT_X_DISCLUTPAR's round does, for its
//! discontiguous tables of fixed width): the width of this table grows with the query
//! (`choose_na_lookup_table`: 4^11 cells from 8,500 table entries, 4^12 from 300,000 for word sizes
//! above 12), so where the round runs nearly every page is written, and three E. coli genomes
//! complete no 2 MiB region of a 4^12 table in a first chunk of 163,840 positions (6 of 32 in one
//! of 2^19) while they leave 9 of its 16,385 pages unwritten in 2,000,000. Under transparent huge
//! pages (`enabled=always`) a first write can map a whole 2 MiB page where the loop's first read
//! maps the huge zero page; such a page holds marked pages and at most the pages of the table that
//! no word of the whole fill writes (on this WSL machine, THP `always`, megablast of 5 kb to 4.6 Mb of
//! MG1655 against LC738868 at 8 threads kept its max RSS within 0.4 MiB of the fill without the
//! round). Memory: with the round, the readers' marks, one bit per page per reader on whole cache
//! lines (576 bytes per reader at 4^11 cells and 2,112 at 4^12 for a table that starts 16 bytes into
//! a page, as glibc's mapped blocks do, 512 and 2,048 for a page-aligned one: 33 KiB at P = 16, 526
//! KiB at 255), freed at the end of the fill, and their OR (one row), freed after the round; in
//! every fill, P counters of the zeros written.
//!
//! EXPERIMENT (LOSAT_X_MB24): a table of 3-byte elements (`Tables::x_h24`, `x_mb_narrow.rs`) runs
//! its own instances of phases B and C, of the marking and of the prefault round (`H24`): phase B
//! reads and writes each head as exactly its 3 bytes (the next cell may be another owner's) and
//! phase C each link (the next offset may be another reader's, or linked already); the marks are
//! the pages of the words that hold a cell's first and last byte, and a page's toucher is the
//! owner of the cell of its middle word's first byte (the round's jobs write different words, so a
//! word with bytes of two owners' cells is still written by one job). The plan, the phases and the
//! order of the updates are the same, so the tables are the same.

use super::x_disc_lookup_par::{
    middle, page_offset, pages_of, touch, PageBits, SharedBits, TouchLog, HASHTABLE,
    PAGE_ELEMS_SHIFT, PREFAULT_MIN_WORDS_PER_PAGE,
};
use super::x_mb_lookup_par::{
    flush_private, lines_capacity, plan, timing_on, Copies, Line, Private, Shared, SharedTables,
    Tables, Unready, Words, K_COMPRESSION_FACTOR_SHIFT,
};
use super::{build_unmasked_ranges, x_range_words, PvArrayType, PV_ARRAY_MASK};
use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart (scheduling only): the most stream positions per chunk (each chunk is
/// read, updated and linked before the next; the arena holds the words of one chunk), as
/// LOSAT_X_DISCLUTPAR.
const CHUNK_POSITIONS: usize = 1 << 19;

/// No NCBI counterpart (scheduling only): the stream positions per chunk for each worker of the
/// fill, up to `CHUNK_POSITIONS` (`chunk_positions`).
const CHUNK_POSITIONS_PER_PART: usize = 1 << 15;

/// No NCBI counterpart (scheduling only): the positions per chunk of a fill of `parts` workers,
/// `min(2^19, parts << 15)`: the arena (16 bytes a position) is about `parts x 0.5 MiB` up to 16
/// workers and stops at 8 MiB beyond (the size was chosen for the query-batch slots of
/// LOSAT_X_QBATCHPAR, removed on 2026-10-28, whose fills ran at the same time). The chunk size
/// changes no element of the tables (the module documentation).
fn chunk_positions(parts: usize) -> usize {
    CHUNK_POSITIONS.min(parts.clamp(1, 255) * CHUNK_POSITIONS_PER_PART)
}

/// No NCBI counterpart (scheduling only): how far ahead phase B prefetches the cells of a
/// sub-range on x86_64 (as `x_mb_lookup_par`'s fill and LOSAT_X_DISCLUTPAR's phase B).
#[cfg(target_arch = "x86_64")]
const X_AHEAD: usize = 16;

/// No NCBI counterpart (LOSAT_TIMING only): one phase's times, summed over the chunks of the
/// search's fills, in ns (see the module documentation).
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

// No NCBI counterpart (LOSAT_TIMING only): the phases' split (the prefault round last,
// `PREFAULT_PHASE`), the fills and chunks of the search's read-once fills, the largest arena of one
// fill (bytes), and the prefault rounds and their page touches (the zeros written).
static PHASE_TIMES: [PhaseTimes; 4] = [const { PhaseTimes::new() }; 4];
const PHASE_NAMES: [&str; 4] = ["read", "update", "last", "prefault"];
/// The index of the prefault round in `PHASE_TIMES`.
const PREFAULT_PHASE: usize = 3;
static FILLS: AtomicU64 = AtomicU64::new(0);
static CHUNKS: AtomicU64 = AtomicU64::new(0);
static ARENA_PEAK_BYTES: AtomicU64 = AtomicU64::new(0);
static PREFAULT_ROUNDS: AtomicU64 = AtomicU64::new(0);
static PREFAULT_TOUCHES: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart (scheduling only): LOSAT_X_MBLUTPREFAULT, the prefault round (module
/// documentation), read once.
fn prefault_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBLUTPREFAULT"))
}

/// No NCBI counterpart: with LOSAT_TIMING, when a read-once fill of the search ran, the `[TIMING]
/// x_mb_lookup_once` line (module documentation), printed at exit after `x_mb_lookup_par`'s lines.
/// It prints counters and changes no value NCBI computes.
pub(super) fn print_timing() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if !timing_on() || load(&FILLS) == 0 {
        return;
    }
    let s = |c: &AtomicU64| load(c) as f64 / 1e9;
    let mut line = format!(
        "[TIMING] x_mb_lookup_once: fills={} chunks={}",
        load(&FILLS),
        load(&CHUNKS)
    );
    for (name, t) in PHASE_NAMES.iter().zip(&PHASE_TIMES) {
        line.push_str(&format!(
            " {name}_wall={:.4}s {name}_max={:.4}s {name}_sum={:.4}s {name}_late={:.4}s",
            s(&t.wall),
            s(&t.max),
            s(&t.sum),
            s(&t.late)
        ));
    }
    line.push_str(&format!(
        " arena={:.1}MiB prefault_rounds={} prefault_touches={}",
        load(&ARENA_PEAK_BYTES) as f64 / 1048576.0,
        load(&PREFAULT_ROUNDS),
        load(&PREFAULT_TOUCHES)
    ));
    eprintln!("{line}");
}

/// No NCBI counterpart (scheduling only): one unmasked range of one query that the loop of
/// `build_mb_lookup` reads, and the offset of its first position in the stream of all positions the
/// loop reads, in reading order.
struct Segment<'a> {
    seq: &'a [u8],
    start: usize,
    end: usize,
    query_offset: usize,
    stream: usize,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1019-1022
/// The ranges the loop of `build_mb_lookup` reads, in its order (the same `build_unmasked_ranges`
/// per query; the queries shorter than `lut_word_length` and the ranges shorter than `word_length`
/// skipped as there), and the number of positions it reads. In debug builds it checks the premise
/// of phase C that each query offset belongs to one word (`distinct_query_offsets`).
fn segments<'a>(words: &Words<'a>) -> (Vec<Segment<'a>>, usize) {
    debug_assert!(
        distinct_query_offsets(words),
        "the query intervals [offset, offset + length) of a megablast table overlap"
    );
    let mut segments = Vec::new();
    let mut stream = 0usize;
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        let seq: &'a [u8] = seq;
        if seq.len() < words.lut_word_length {
            continue;
        }
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let query_offset = words.query_offsets[q_idx].max(0) as usize;
        for (start, end) in build_unmasked_ranges(seq.len(), masks) {
            let range_len = end.saturating_sub(start);
            if words.word_length > range_len || range_len == 0 {
                continue;
            }
            segments.push(Segment {
                seq,
                start,
                end,
                query_offset,
                stream,
            });
            stream += range_len;
        }
    }
    (segments, stream)
}

/// No NCBI counterpart (debug builds only): whether the intervals `[offset, offset + length)` of
/// the queries that the loop reads (those of `lut_word_length` letters or more) are disjoint. The
/// words of a query have the 1-based offsets
/// `offset + 1 ..= offset + length - lut_word_length + 1`, inside that interval, so disjoint
/// intervals give every word its own offset and every `next_pos` element one writer in phase C.
/// The callers place the queries one after another with a separator (`offset + length + 1`),
/// which is disjoint.
fn distinct_query_offsets(words: &Words<'_>) -> bool {
    let mut spans: Vec<(usize, usize)> = words
        .queries_blastna
        .iter()
        .zip(words.query_offsets)
        .filter(|(seq, _)| seq.len() >= words.lut_word_length)
        .map(|(seq, &offset)| (offset.max(0) as usize, seq.len()))
        .collect();
    spans.sort_unstable();
    spans.windows(2).all(|w| w[0].0 + w[0].1 <= w[1].0)
}

/// No NCBI counterpart (scheduling only): a word in the arena: its 1-based query offset `index`
/// and, until phase B updates it, its cell; after that, the cell's head before the word's update,
/// the value of `next_pos[index]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(C)]
struct Entry {
    index: u32,
    value: u32,
}

const ENTRY_BYTES: usize = std::mem::size_of::<Entry>();
const LINE_BYTES: usize = std::mem::size_of::<Line>();
const _: () = assert!(ENTRY_BYTES == 8 && LINE_BYTES.is_multiple_of(ENTRY_BYTES));
/// Entries per 64-byte line.
const LINE_ENTRIES: usize = LINE_BYTES / ENTRY_BYTES;
/// Row values (u32) per 64-byte line.
const LINE_WORDS: usize = LINE_BYTES / std::mem::size_of::<u32>();

/// No NCBI counterpart (scheduling only): the layout of one fill's arena, in 64-byte lines: P
/// reader regions of `region_lines` lines (scratch, then sorted, `half` entries each), then P rows
/// of `row_lines` lines (the bounds of a reader's sub-ranges: P + 1 values).
#[derive(Clone, Copy, Debug)]
struct Shape {
    parts: usize,
    half: usize,
    region_lines: usize,
    row_lines: usize,
}

impl Shape {
    /// The layout for chunks of `chunk` positions and `parts` readers: a reader's span has at most
    /// `ceil(chunk / parts)` positions, so at most that many words.
    fn new(chunk: usize, parts: usize) -> Self {
        let half = chunk.div_ceil(parts).div_ceil(LINE_ENTRIES) * LINE_ENTRIES;
        Self {
            parts,
            half,
            region_lines: 2 * half / LINE_ENTRIES,
            row_lines: (parts + 1).div_ceil(LINE_WORDS),
        }
    }

    /// The lines of the arena.
    fn lines(&self) -> usize {
        self.parts * (self.region_lines + self.row_lines)
    }
}

/// No NCBI counterpart (scheduling only): the arena of one fill, shared by the jobs of a phase,
/// each of which uses only its own parts: phase A's job `a` reader `a`'s region and row, phase
/// B's job `o` the sub-ranges `(*, o)` (and reads every row), phase C's job `a` reader `a`'s
/// region (and row); the phases are separated by joins.
#[derive(Clone, Copy)]
struct SharedArena {
    base: *mut Line,
    shape: Shape,
}

// SAFETY: the jobs of a phase use disjoint parts of the arena (see `SharedArena`), the phases are
// separated by joins, and the arena outlives the phases.
unsafe impl Send for SharedArena {}
// SAFETY: as above.
unsafe impl Sync for SharedArena {}

impl SharedArena {
    /// Reader `a`'s scratch region (`shape.half` entries). Nothing is dereferenced here.
    #[inline(always)]
    fn scratch(self, a: usize) -> *mut Entry {
        debug_assert!(a < self.shape.parts);
        self.base.wrapping_add(a * self.shape.region_lines) as *mut Entry
    }

    /// Reader `a`'s sorted region (`shape.half` entries).
    #[inline(always)]
    fn sorted(self, a: usize) -> *mut Entry {
        self.scratch(a).wrapping_add(self.shape.half)
    }

    /// Reader `a`'s row: the start of each owner's sub-range, then the number of words.
    #[inline(always)]
    fn row(self, a: usize) -> *mut u32 {
        debug_assert!(a < self.shape.parts);
        self.base
            .wrapping_add(self.shape.parts * self.shape.region_lines + a * self.shape.row_lines)
            as *mut u32
    }
}

/// No NCBI counterpart (scheduling only): one owner's private copies of `helper_array` and
/// `pv_array`, zeroed by the owner in its first phase-B job.
struct OwnerCopies<'c> {
    unready: Option<Unready<'c>>,
    private: Option<Private<'c>>,
}

impl<'c> OwnerCopies<'c> {
    /// The copies, zeroed (on the calling thread) the first time.
    fn private(&mut self) -> &mut Private<'c> {
        let unready = &mut self.unready;
        self.private
            .get_or_insert_with(|| unready.take().expect("every owner has its copies").zeroed())
    }
}

/// No NCBI counterpart (scheduling only): one value per job, used only by that job.
struct Slots<T>(Vec<UnsafeCell<T>>);

// SAFETY: job `i` of a phase uses only slot `i` (`Slots::at`), and the phases are separated by
// joins.
unsafe impl<T: Send> Sync for Slots<T> {}

impl<T> Slots<T> {
    /// SAFETY: no other job of the current phase uses slot `i`.
    #[allow(clippy::mut_from_ref)]
    #[inline(always)]
    unsafe fn at(&self, i: usize) -> &mut T {
        // SAFETY: the caller's guarantee.
        unsafe { &mut *self.0[i].get() }
    }
}

/// No NCBI counterpart (scheduling only, tests only): what the fill checked on its way.
#[cfg(test)]
#[derive(Default)]
struct Checks {
    owners_with_values: std::sync::atomic::AtomicUsize,
    subranges: std::sync::atomic::AtomicUsize,
    subranges_in_order: std::sync::atomic::AtomicUsize,
}

/// No NCBI counterpart (scheduling only): the prefault round's state of one fill (module
/// documentation): the page offset (`page_offset`) of `hashtable`, the readers' marks of the pages
/// of their words' cells in the first chunk (one row per reader; no rows without the round), and
/// each job's count of the zeros it wrote (in tests also its log of them).
#[derive(Clone, Copy)]
struct Faults<'a> {
    offset: usize,
    marks: SharedBits,
    touched: &'a Slots<(usize, TouchLog)>,
}

/// No NCBI counterpart (scheduling only): the plan and the shared state of one fill.
struct Fill<'a, 'w, 'c> {
    words: &'a Words<'w>,
    segments: &'a [Segment<'w>],
    owner: &'a [u8],
    block_shift: usize,
    parts: usize,
    tables: SharedTables,
    arena: SharedArena,
    copies: &'a Slots<OwnerCopies<'c>>,
    faults: Faults<'a>,
    #[cfg(test)]
    checks: &'a Checks,
}

impl Fill<'_, '_, '_> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1045
    /// Phase A, reader `a`: the words that end at the stream positions `[s0, s1)`, read with the
    /// statements of `build_mb_lookup`'s loop (the same resets, k-mer and offsets), into its
    /// scratch region in reading order, then copied stably into its sorted region grouped by the
    /// owner of the cell; its row gets the start of each owner's sub-range and the number of words.
    /// A range is entered up to `lut_word_length - 1` positions before `s0` (not before its start),
    /// so that the k-mer and the words (`x_range_words`) at the positions from `s0` on are the
    /// loop's. With `MARK` (the first chunk, when the prefault round runs) it also marks in its own
    /// bits (`Faults::marks`), after clearing them, the page of `hashtable` of each word's cell; the
    /// other chunks run the loop without the marking code. `H24` (LOSAT_X_MB24, with `MARK`): the
    /// cells are 3-byte elements, and the pages marked those of the words that hold the first and
    /// the last byte of the cell's element.
    ///
    /// SAFETY: `a < parts`, no other job of the phase uses reader `a`'s region, row or marks, and
    /// the span has at most `shape.half` positions (`s1 - s0 <= half`): each caller's `span` cuts a
    /// chunk of `L <= chunk.min(total)` positions into `parts` spans of at most `ceil(L / parts)`
    /// positions (`floor(L(a+1)/P) - floor(La/P) <= ceil(L/P)`), and the `half` of
    /// `Shape::new(chunk.min(total), parts)` is `ceil(chunk.min(total) / parts)` rounded up.
    unsafe fn read<const MARK: bool, const H24: bool>(&self, a: usize, s0: usize, s1: usize) {
        let words = self.words;
        let lut_word_length = words.lut_word_length;
        let (marks, cells) = (self.faults.marks, self.faults.offset);
        if MARK {
            // SAFETY: reader `a`'s marks are this job's.
            unsafe { marks.clear(a) };
        }
        let scratch = self.arena.scratch(a);
        let mut counts = [0u32; 256];
        let mut n = 0usize;
        let mut g = self
            .segments
            .partition_point(|seg| seg.stream + (seg.end - seg.start) <= s0);
        while s0 < s1 && g < self.segments.len() && self.segments[g].stream < s1 {
            let seg = &self.segments[g];
            let len = seg.end - seg.start;
            let first_pos = seg.start + (s0.max(seg.stream) - seg.stream);
            let end_pos = seg.start + (s1.min(seg.stream + len) - seg.stream);
            let warm_start = seg.start.max(first_pos.saturating_sub(lut_word_length - 1));
            x_range_words(
                &seg.seq[warm_start..end_pos],
                lut_word_length,
                words.kmer_mask,
                || {},
                |current_kmer, i| {
                    let pos = warm_start + i;
                    // A word ends `lut_word_length - 1` letters or more after `warm_start`, so at
                    // or after `first_pos` (`warm_start` is at most that many letters before it).
                    debug_assert!(pos >= first_pos);
                    let cell = current_kmer as usize;
                    let q_off_1 = (seg.query_offset + (pos + 1 - lut_word_length) + 1) as u32;
                    debug_assert!(n < self.arena.shape.half);
                    // SAFETY: `n < half`, the entries of reader `a`'s scratch region: a word is
                    // written only at a stream position of `[s0, s1)` (the segments' stream ranges
                    // follow each other without overlap, `segments`), one at most per position, and
                    // `s1 - s0 <= half` (this function's contract).
                    unsafe {
                        scratch.add(n).write(Entry {
                            index: q_off_1,
                            value: cell as u32,
                        })
                    };
                    n += 1;
                    counts[self.owner[cell >> self.block_shift] as usize] += 1;
                    if MARK && H24 {
                        // SAFETY: reader `a`'s bits; the pages are `hashtable`'s (the element's
                        // bytes lie in its words).
                        unsafe {
                            marks.mark(a, (cells + ((3 * cell) >> 2)) >> PAGE_ELEMS_SHIFT);
                            marks.mark(a, (cells + ((3 * cell + 2) >> 2)) >> PAGE_ELEMS_SHIFT);
                        }
                    } else if MARK {
                        // SAFETY: reader `a`'s bits; the page is one of `hashtable`'s (`cell` < its
                        // length).
                        unsafe { marks.mark(a, (cells + cell) >> PAGE_ELEMS_SHIFT) };
                    }
                },
            );
            g += 1;
        }
        let row = self.arena.row(a);
        let mut at = [0u32; 256];
        let mut sum = 0u32;
        for o in 0..self.parts {
            at[o] = sum;
            // SAFETY: reader `a`'s row has `parts + 1` values.
            unsafe { row.add(o).write(sum) };
            sum += counts[o];
        }
        // SAFETY: as above.
        unsafe { row.add(self.parts).write(sum) };
        // SAFETY: the first `n` entries of the scratch region were written above.
        let read = unsafe { std::slice::from_raw_parts(scratch, n) };
        let sorted = self.arena.sorted(a);
        for e in read {
            let o = self.owner[(e.value as usize) >> self.block_shift] as usize;
            // SAFETY: `at[o]` < the end of owner `o`'s sub-range <= `n` <= `half`.
            unsafe { sorted.add(at[o] as usize).write(*e) };
            at[o] += 1;
        }
    }

    /// The prefault round (module documentation; between phases A and B of the first chunk only, and
    /// only when that chunk has at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per page of
    /// `hashtable`, `fill_chunked`), toucher `j`: writes 0 to the middle element of each page in
    /// `marked` (the OR of the readers' marks, by page of `hashtable`) whose middle element lies in
    /// a block of owner `j`, so that phase B finds the page mapped. A marked page holds the cell of a
    /// word of the first chunk, which phase B writes; no phase B has run, so the element holds the 0
    /// it was allocated with, and each page has one toucher: the round writes 0 over 0 and changes
    /// no value (no NCBI counterpart; scheduling only). `H24` (LOSAT_X_MB24): `hashtable` holds
    /// 3-byte cells in its words, and the middle word goes to the owner of the cell of its first
    /// byte (the jobs of the round write different words, so a word that holds bytes of two
    /// owners' cells is still written by one job only).
    ///
    /// SAFETY: no other job of the phase uses slot `j` of `Faults::touched`; the round runs before
    /// the first phase B, on a table allocated zeroed; `marked` holds pages of `hashtable` only.
    unsafe fn prefault<const H24: bool>(&self, j: usize, marked: &[u64]) {
        let (offset, len) = (self.faults.offset, self.tables.hashtable.len);
        // SAFETY: slot `j` is this job's.
        let (touches, log) = unsafe { self.faults.touched.at(j) };
        for (w, &bits) in marked.iter().enumerate() {
            let mut bits = bits;
            while bits != 0 {
                let page = (w << 6) + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                let m = middle(offset, len, page);
                let cell = if H24 {
                    (4 * m / 3).min((self.owner.len() << self.block_shift) - 1)
                } else {
                    m
                };
                if self.owner[cell >> self.block_shift] as usize == j {
                    // SAFETY: as above: element `m` of a page of this job, before the first phase B.
                    unsafe { touch(self.tables.hashtable, m, HASHTABLE, log) };
                    *touches += 1;
                }
            }
        }
    }

    /// Phase B, owner `o`: the updates of the words of its cells, reader by reader, each sub-range
    /// in reading order (both reversed for `ascending_cells`), counting into its own copies; each
    /// entry keeps the cell's head before the update. `H24` (LOSAT_X_MB24): on 3-byte cells.
    ///
    /// SAFETY: no other job of the phase uses the sub-ranges `(*, o)`, the cells of `o`'s blocks or
    /// `o`'s copies.
    unsafe fn update<const H24: bool>(&self, o: usize) {
        // SAFETY: slot `o` is this job's.
        let x = unsafe { self.copies.at(o) }.private();
        let ascending = self.words.ascending_cells;
        for k in 0..self.parts {
            let a = if ascending { self.parts - 1 - k } else { k };
            let row = self.arena.row(a);
            // SAFETY: reader `a` wrote its row in phase A (before the join); nothing writes it now.
            let (lo, hi) = unsafe { (row.add(o).read() as usize, row.add(o + 1).read() as usize) };
            // SAFETY: sub-range `(a, o)` was written in phase A and is this job's alone.
            let entries =
                unsafe { std::slice::from_raw_parts_mut(self.arena.sorted(a).add(lo), hi - lo) };
            #[cfg(test)]
            if !entries.is_empty() {
                self.checks.subranges.fetch_add(1, Ordering::Relaxed);
                if entries.windows(2).all(|w| w[0].index < w[1].index) {
                    self.checks
                        .subranges_in_order
                        .fetch_add(1, Ordering::Relaxed);
                }
            }
            // SAFETY: the cells of the entries lie in `o`'s blocks; `x` is `o`'s.
            unsafe {
                if ascending {
                    update_backward::<H24>(self.tables, x, entries);
                } else {
                    update_forward::<H24>(self.tables, x, entries);
                }
            }
        }
    }

    /// Phase C, reader `a`: `next_pos[index]` of the words it read, from their kept heads. `H24`
    /// (LOSAT_X_MB24): as 3-byte elements, each written as exactly its 3 bytes (the next element may
    /// be another reader's, or written already).
    ///
    /// SAFETY: no other job of the phase uses reader `a`'s region; the indices are the words' own.
    unsafe fn link<const H24: bool>(&self, a: usize) {
        // SAFETY: reader `a` wrote its row in phase A.
        let n = unsafe { self.arena.row(a).add(self.parts).read() } as usize;
        // SAFETY: the first `n` entries of the sorted region were written in phase A.
        let entries = unsafe { std::slice::from_raw_parts(self.arena.sorted(a), n) };
        for e in entries {
            // SAFETY: each index is written once, by the reader of its word.
            unsafe {
                if H24 {
                    let next_pos = self.tables.next_pos;
                    debug_assert!(3 * e.index as usize + 3 <= 4 * next_pos.len);
                    super::x_mb_narrow::store_exact(next_pos.ptr, e.index as usize, e.value);
                } else {
                    self.tables.next_pos.set(e.index as usize, e.value)
                }
            };
        }
    }

    /// `link` for a table of 3-byte cells (`h24`, LOSAT_X_MB24) or not.
    ///
    /// SAFETY: as `link`'s.
    #[inline(always)]
    unsafe fn link_any(&self, h24: bool, a: usize) {
        // SAFETY: the caller's.
        unsafe {
            if h24 {
                self.link::<true>(a)
            } else {
                self.link::<false>(a)
            }
        }
    }

    /// After the last phase B, owner `o`: its copies added into `helper_array` and ORed into
    /// `pv_array` for its own blocks (`x_mb_lookup_par::flush_private`).
    ///
    /// SAFETY: no other job of the phase uses `o`'s copies or the counters and presence-vector
    /// words of `o`'s blocks.
    unsafe fn flush(&self, o: usize) {
        // SAFETY: slot `o` is this job's.
        let slot = unsafe { self.copies.at(o) };
        if let Some(x) = slot.private.as_ref() {
            #[cfg(test)]
            if x.has_values() {
                self.checks
                    .owners_with_values
                    .fetch_add(1, Ordering::Relaxed);
            }
            flush_private(o as u8, self.owner, self.block_shift, self.tables, x);
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
/// The update of one word, as in the loop of `build_mb_lookup`, with the counter and the presence
/// bit in the owner's private copies; `next_pos[index]` is written in phase C from the head kept in
/// the entry. `H24` (LOSAT_X_MB24): on a 3-byte cell, read and written as exactly its 3 bytes (the
/// next cell may be another owner's).
///
/// SAFETY: the cell lies in a block of the calling owner, and `x` is its copies.
#[inline(always)]
unsafe fn update_one<const H24: bool>(t: SharedTables, x: &mut Private<'_>, e: &mut Entry) {
    let cell = e.value as usize;
    // SAFETY: the caller's guarantee (the cell < the table size: a masked k-mer).
    unsafe {
        let head = if H24 {
            debug_assert!(3 * cell + 3 <= 4 * t.hashtable.len);
            super::x_mb_narrow::load_exact(t.hashtable.ptr, cell)
        } else {
            t.hashtable.get(cell)
        };
        if head == 0 {
            let w = cell >> t.pv_array_bts;
            if w < t.pv_array.len {
                *x.pv.slot(w) |= (1 as PvArrayType) << (cell & PV_ARRAY_MASK);
            }
        } else {
            let c = x.counters.slot(cell >> K_COMPRESSION_FACTOR_SHIFT);
            *c = c.saturating_add(1);
        }
        if H24 {
            super::x_mb_narrow::store_exact(t.hashtable.ptr, cell, e.index);
        } else {
            t.hashtable.set(cell, e.index);
        }
        e.value = head;
    }
}

/// The updates of a sub-range's entries in order (phase B of a descending table).
///
/// SAFETY: as `update_one`, for every entry.
#[inline]
unsafe fn update_forward<const H24: bool>(
    t: SharedTables,
    x: &mut Private<'_>,
    entries: &mut [Entry],
) {
    for j in 0..entries.len() {
        #[cfg(target_arch = "x86_64")]
        if let Some(ahead) = entries.get(j + X_AHEAD) {
            prefetch::<H24>(t.hashtable, ahead.value as usize);
        }
        // SAFETY: the caller's guarantee.
        unsafe { update_one::<H24>(t, x, &mut entries[j]) };
    }
}

/// The updates of a sub-range's entries back to front (phase B of `ascending_cells`).
///
/// SAFETY: as `update_one`, for every entry.
#[inline]
unsafe fn update_backward<const H24: bool>(
    t: SharedTables,
    x: &mut Private<'_>,
    entries: &mut [Entry],
) {
    for j in (0..entries.len()).rev() {
        #[cfg(target_arch = "x86_64")]
        if j >= X_AHEAD {
            prefetch::<H24>(t.hashtable, entries[j - X_AHEAD].value as usize);
        }
        // SAFETY: the caller's guarantee.
        unsafe { update_one::<H24>(t, x, &mut entries[j]) };
    }
}

/// A prefetch of cell `cell` of `hashtable` (x86_64); `H24`: of its 3-byte element.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn prefetch<const H24: bool>(hashtable: Shared<u32>, cell: usize) {
    debug_assert!(H24 || cell < hashtable.len);
    // SAFETY: a cell index < the table size; a prefetch does not fault.
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            if H24 {
                super::x_mb_narrow::at(hashtable.ptr, cell)
            } else {
                hashtable.ptr.add(cell) as *const i8
            },
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
}

/// No NCBI counterpart (scheduling only): `job(0)`, ..., `job(parts - 1)`, on the pool of the
/// calling worker (`in_pool`) or one after another; returns when all have returned. The jobs never
/// wait for anything. (A copy of LOSAT_X_DISCLUTPAR's `run_jobs` in `x_disc_lookup_par.rs`, kept
/// here until that module's functions are shared.)
fn run_jobs<F: Fn(usize) + Sync>(parts: usize, in_pool: bool, job: F) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if in_pool && parts > 1 {
        // The scope runs on the pool of the calling worker; outside a pool it would start
        // rayon's global pool, which `x_mb_lookup_par::parts` never asks for.
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
/// in ns after the phase began (a copy of LOSAT_X_DISCLUTPAR's `JobClock`).
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

/// No NCBI counterpart (scheduling only): one phase, `run_jobs`; with a clock (LOSAT_TIMING, the
/// search's fills only) each job's start and end are taken and the phase's wall time, longest job,
/// summed job time and latest start are added to `PHASE_TIMES[phase]`. The jobs are the same
/// either way. (A copy of LOSAT_X_DISCLUTPAR's `run_phase`.)
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

/// No NCBI counterpart (scheduling only): what one fill did: its chunks, the lines of its arena,
/// its prefault rounds (0 or 1) and the zeros they wrote; in tests also the owners whose copies
/// held a value before the flush, the non-empty sub-ranges phase B took and how many of them were
/// in increasing query offset (reading) order, the prefault round's writes (`Touch`) and the pages
/// of `hashtable` (by address) that the first chunk wrote.
#[derive(Clone, Debug, Default)]
pub(super) struct OnceStats {
    pub(super) chunks: usize,
    pub(super) arena_lines: usize,
    pub(super) prefault_rounds: usize,
    pub(super) touches: usize,
    #[cfg(test)]
    pub(super) owners_with_values: usize,
    #[cfg(test)]
    subranges: usize,
    #[cfg(test)]
    subranges_in_order: usize,
    #[cfg(test)]
    touch_log: TouchLog,
    #[cfg(test)]
    first_pages: std::collections::BTreeSet<usize>,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1033,1081-1091
/// Fills the tables as the word loop of `build_mb_lookup` does (with its re-threading for
/// `ascending_cells`), each query position read once, with `parts` workers on the search pool (or
/// one after another, `in_pool == false`, where the unit tests run it without a pool); the module
/// documentation gives the phases and why every element of the four arrays is the loop's. `timed`:
/// a build of the search (LOSAT_TIMING's `[TIMING] x_mb_lookup_once` counts it). `hashtable`
/// must be zeroed (as `build_mb_lookup` allocates it, NCBI's `calloc`): the loop's
/// empty cells, and the premise of the prefault round (LOSAT_X_MBLUTPREFAULT). `resident`: a table
/// kept by LOSAT_X_MBLUTKEEP and zeroed, whose pages are mapped: no prefault round
/// (`x_mb_lut_keep.rs`, "The prefault round").
pub(super) fn fill(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    timed: bool,
    resident: bool,
) -> OnceStats {
    if resident && prefault_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_MBLUTPREFAULT", false);
    }
    fill_chunked(
        parts,
        words,
        tables,
        in_pool,
        timed,
        chunk_positions(parts),
        prefault_on() && !resident,
    )
}

/// No NCBI counterpart (scheduling only): whether a first chunk of `first` positions runs the prefault
/// round (LOSAT_X_MBLUTPREFAULT) on a table of `cells` cells: at least `PREFAULT_MIN_WORDS_PER_PAGE`
/// positions per 1,024 cells (4 KiB), counted on the table's size, not on the pages it spans, so
/// that the allocation's offset in its first page (16 bytes for glibc's mapped blocks, 0 for a
/// page-aligned one) does not move the boundary: 131,072 positions for 4^12 cells, the first chunk
/// of 4 workers (`chunk_positions`). With LOSAT_X_MB24 `cells` is the words of the 3-byte cells
/// (1,024 to a page): 98,304 positions for 4^12 cells.
fn prefault_gate(first: usize, cells: usize) -> bool {
    first > 0
        && first >= PREFAULT_MIN_WORDS_PER_PAGE.saturating_mul((cells >> PAGE_ELEMS_SHIFT).max(1))
}

/// `fill` with the chunk size and LOSAT_X_MBLUTPREFAULT (`prefault`) given (the tests take others).
fn fill_chunked(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    timed: bool,
    chunk: usize,
    prefault: bool,
) -> OnceStats {
    // LOSAT_X_MB24: a table of 3-byte cells (`x_mb_narrow.rs`) has `kmer_mask + 1` cells in fewer
    // words; its phases run their own instances (`H24`).
    let h24 = tables.x_h24;
    let table_size = if h24 {
        words.kmer_mask as usize + 1
    } else {
        tables.hashtable.len()
    };
    let (block_shift, owner) = plan(table_size, tables.pv_array_bts, parts);
    debug_assert!(
        block_shift >= K_COMPRESSION_FACTOR_SHIFT && block_shift >= tables.pv_array_bts,
        "a block must own whole counters and presence-vector words"
    );
    debug_assert!(!(h24 && words.ascending_cells));
    let parts = parts.clamp(1, 255);
    let chunk = chunk.max(1);
    let (segments, total) = segments(words);
    let mut stats = OnceStats::default();
    if total == 0 {
        // No position: the loop reads nothing either.
        return stats;
    }
    let shared = SharedTables {
        hashtable: Shared::new(tables.hashtable),
        next_pos: Shared::new(tables.next_pos),
        pv_array: Shared::new(tables.pv_array),
        pv_array_bts: tables.pv_array_bts,
        helper_array: Shared::new(tables.helper_array),
    };
    let mut copies = Copies::new(parts, shared.helper_array.len, shared.pv_array.len);
    let slots = Slots(
        copies
            .parts()
            .map(|unready| {
                UnsafeCell::new(OwnerCopies {
                    unready: Some(unready),
                    private: None,
                })
            })
            .collect(),
    );
    // A chunk holds at most `min(chunk, total)` positions: a small query gets a small arena.
    let shape = Shape::new(chunk.min(total), parts);
    let mut arena: Vec<Line> = Vec::with_capacity(lines_capacity(shape.lines()));
    let chunks = total.div_ceil(chunk);
    // Chunk `k` of the processing order: the chunks in stream order, last to first for
    // `ascending_cells`.
    let bounds = |k: usize| {
        let c = if words.ascending_cells {
            chunks - 1 - k
        } else {
            k
        };
        (c * chunk, total.min((c + 1) * chunk))
    };
    // The prefault round runs in the first chunk (of the processing order) only, and only when
    // that chunk has at least `PREFAULT_MIN_WORDS_PER_PAGE` positions per 4 KiB of `hashtable`
    // (`prefault_gate`); without it the readers' marks have no rows. (The words of `hashtable`:
    // its cells, or with LOSAT_X_MB24 the words that hold them.)
    let offset = page_offset(shared.hashtable.ptr);
    let pages = pages_of(offset, shared.hashtable.len);
    let (first_start, first_end) = bounds(0);
    let round_on = prefault && prefault_gate(first_end - first_start, shared.hashtable.len);
    if prefault {
        crate::utils::xreach::x_reach!("LOSAT_X_MBLUTPREFAULT", round_on);
    }
    let mut marks = PageBits::new(if round_on { pages } else { 0 }, parts);
    let touched = Slots(
        (0..parts)
            .map(|_| UnsafeCell::new((0usize, TouchLog::default())))
            .collect(),
    );
    #[cfg(test)]
    let checks = Checks::default();
    let fill = Fill {
        words,
        segments: &segments,
        owner: &owner,
        block_shift,
        parts,
        tables: shared,
        arena: SharedArena {
            base: arena.spare_capacity_mut().as_mut_ptr() as *mut Line,
            shape,
        },
        copies: &slots,
        faults: Faults {
            offset,
            marks: marks.shared(),
            touched: &touched,
        },
        #[cfg(test)]
        checks: &checks,
    };
    let span = |(start, end): (usize, usize), a: usize| {
        (
            start + (end - start) * a / parts,
            start + (end - start) * (a + 1) / parts,
        )
    };
    let mut clock = (timed && in_pool && timing_on()).then(|| JobClock::new(parts));
    // SAFETY (all phases): job `i` of a phase uses only reader `i`'s region, row and marks, owner
    // `i`'s sub-ranges, cells and copies, the `next_pos` entries of reader `i`'s words, and toucher
    // `i`'s pages and slot of `touched` (see `SharedArena`, `Fill::read`, `Fill::prefault`,
    // `Fill::update`, `Fill::link`, `Fill::flush`); `run_jobs` returns after every job of the phase
    // has returned.
    for k in 0..chunks {
        let current = bounds(k);
        // In the first chunk (with `round_on`) the readers mark their words' pages and the prefault
        // round runs before phase B (`Fill::prefault`); later chunks run neither.
        let mark = round_on && k == 0;
        run_phase(&mut clock, 0, parts, in_pool, |a| {
            if k > 0 {
                // C of the chunk before.
                // SAFETY: phase C of job `a < parts` links reader `a`'s words only (all phases).
                unsafe { fill.link_any(h24, a) };
            }
            let (s0, s1) = span(current, a);
            // Reader `a < parts` reads its span of chunk `current` (at most `chunk.min(total)`).
            // SAFETY: the reader is this job's (all phases); the span fits `shape.half` (`read`).
            match (mark, h24) {
                (true, true) => unsafe { fill.read::<true, true>(a, s0, s1) },
                (true, false) => unsafe { fill.read::<true, false>(a, s0, s1) },
                // (the words are read alike; only the marks depend on the cells' width)
                (false, _) => unsafe { fill.read::<false, false>(a, s0, s1) },
            }
        });
        if mark {
            // Between two phases: the OR of the readers' marks.
            let words = pages.div_ceil(64);
            let marked: Vec<u64> = (0..words)
                .map(|w| {
                    (0..parts).fold(0, |bits, r| {
                        // SAFETY: no job is running; `w` < the rows' stride.
                        bits | unsafe { fill.faults.marks.word(r, w) }
                    })
                })
                .collect();
            // SAFETY: toucher `j < parts` alone uses slot `j` of `touched` (all phases); the round
            // runs between phase A and the first phase B, on the zeroed `hashtable`, and `marked`
            // (the OR of the readers' marks) holds pages of `hashtable` only (`Fill::prefault`).
            run_phase(&mut clock, PREFAULT_PHASE, parts, in_pool, |j| unsafe {
                if h24 {
                    fill.prefault::<true>(j, &marked)
                } else {
                    fill.prefault::<false>(j, &marked)
                }
            });
            stats.prefault_rounds += 1;
        }
        // SAFETY: phase B of owner `o < parts` uses only the sub-ranges `(*, o)`, the cells of
        // `o`'s blocks and `o`'s copies (all phases); every reader of the chunk has returned.
        run_phase(&mut clock, 1, parts, in_pool, |o| unsafe {
            if h24 {
                fill.update::<true>(o)
            } else {
                fill.update::<false>(o)
            }
        });
        // Tests: the pages of `hashtable` the first chunk wrote (every written cell holds a 1-based
        // query offset, other than 0), which the round must have touched; for 3-byte cells
        // (LOSAT_X_MB24) the pages of the words of their first and last bytes (a cell's last byte
        // may be 0 on a page of its own).
        #[cfg(test)]
        if k == 0 {
            let t = fill.tables.hashtable;
            // SAFETY: between two phases; the table outlives the fill.
            let v = unsafe { std::slice::from_raw_parts(t.ptr, t.len) };
            stats.first_pages = if h24 {
                super::x_mb_narrow::written_pages(v, table_size)
            } else {
                super::x_disc_lookup_par::nonzero_pages(v)
            };
        }
    }
    // SAFETY: after the last phase B, job `i < parts` links reader `i`'s words and flushes
    // owner `i`'s copies into its own blocks' counters and presence words (all phases).
    run_phase(&mut clock, 2, parts, in_pool, |i| unsafe {
        fill.link_any(h24, i);
        fill.flush(i);
    });
    stats.chunks = chunks;
    stats.arena_lines = shape.lines();
    let mut touched = touched.0;
    stats.touches = touched.iter_mut().map(|slot| slot.get_mut().0).sum();
    #[cfg(test)]
    {
        stats.owners_with_values = checks.owners_with_values.into_inner();
        stats.subranges = checks.subranges.into_inner();
        stats.subranges_in_order = checks.subranges_in_order.into_inner();
        stats.touch_log = touched
            .iter_mut()
            .flat_map(|slot| slot.get_mut().1.drain(..))
            .collect();
    }
    if clock.is_some() {
        FILLS.fetch_add(1, Ordering::Relaxed);
        CHUNKS.fetch_add(chunks as u64, Ordering::Relaxed);
        ARENA_PEAK_BYTES.fetch_max((shape.lines() * LINE_BYTES) as u64, Ordering::Relaxed);
        PREFAULT_ROUNDS.fetch_add(stats.prefault_rounds as u64, Ordering::Relaxed);
        PREFAULT_TOUCHES.fetch_add(stats.touches as u64, Ordering::Relaxed);
    }
    stats
}

impl Fill<'_, '_, '_> {
    /// Sub-range `(a, o)` of the arena: the words reader `a` read in its span (`Fill::read`) whose
    /// cells lie in owner `o`'s blocks, in reading order (LOSAT_X_BNCSR's pool build, `csr`).
    ///
    /// SAFETY: reader `a` wrote its row and region in a phase that has ended, and no job writes
    /// them until the slice is dropped.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[inline(always)]
    unsafe fn subrange(&self, a: usize, o: usize) -> &[Entry] {
        let row = self.arena.row(a);
        // SAFETY: the caller's; a row has `parts + 1` values and `o < parts`.
        unsafe {
            let (lo, hi) = (row.add(o).read() as usize, row.add(o + 1).read() as usize);
            std::slice::from_raw_parts(self.arena.sorted(a).add(lo), hi - lo)
        }
    }
}

/// EXPERIMENT (LOSAT_X_BNCSR, `x_bn_csr.rs`): the CSR of a contiguous table built from its words
/// on the search pool of `parts` threads instead of the chains (`fill`), with the presence bits and
/// helper counters into `pv_array` and `helper_array` (zeroed, of the fill's lengths); returns it
/// and whether the pool built it (otherwise `x_bn_csr::serial` did, on this thread: no pool of
/// `parts` threads around it, or a build without threads).
///
/// The words are read by this module's readers (`Fill::read::<false>`, which touches no table,
/// copy or prefault state; chunks of `chunk_positions(parts)` positions, each cut into the
/// readers' spans in stream order) into the arena, once to count and once to place. Owner `o`
/// holds a contiguous run of the cell blocks (block `b` goes to owner `b * parts / blocks`, so the
/// owners' cells follow each other in cell order, as in LOSAT_X_DCCSR's `x_disc_lookup_par::csr`)
/// and makes every write to its cells' counts and cursors (`x_bn_csr.rs`: the count of `cell` is
/// `offsets[cell + count_shift]`):
/// - it zeroes its cells' counts; per chunk, after the readers, for `a` in `0..parts`, for each
///   entry of sub-range `(a, o)` in order: the cell's count + 1 (the count of the entry 16 ahead
///   prefetched), adding up its entries;
/// - after the last chunk, its base is the number of entries of the owners before it, and it runs
///   the sum pass (`x_dc_csr::sum_cells`) over its cells' counts from that base, with the
///   presence-vector words and counters of its blocks (a block owns whole ones);
/// - per chunk, after the readers, for `a` in `0..parts`, for each entry of sub-range `(a, o)` in
///   order, the place of `x_bn_csr::serial`: `qs[--offsets[cell]] = index` (newest first) or
///   `qs[offsets[cell + 1]++] = index` (`ascending_cells`).
///
/// A cell gets its words in the loop's order (chunks in order, the readers' spans in stream order,
/// a sub-range in reading order), so the counts, sums and runs are those of `x_bn_csr::serial`.
/// With `ascending_cells` the cursor of an owner's first cell, `offsets[c0 + 1]`, is the running
/// sum of the cell before it, which the owner before writes in its sum pass; every thread runs its
/// sum pass before it reads the first chunk of the place pass and waits at the barrier after that
/// reading, so no cursor is used before its sum is written, and the owners' cursors
/// (`offsets[c0 + 1 .. c1 + 1]`) are disjoint. The whole build is one `rayon::broadcast` (one job on
/// each of the pool's `parts` threads) whose threads pass an `x_dc_csr::SpinBarrier` after each
/// chunk's reading and after each chunk's owner work (one `rayon::scope` per phase and chunk, as
/// `fill` runs them, started its last job about 0.15 ms late each time on the workstation for the
/// discontiguous tables; v4-R R16). The barrier waits for every thread of the pool, so it needs
/// them free while the table is built: `build_mb_lookup` runs on the search thread before its
/// batch is scanned (`build_lookup_tables`, called by `search_query_batch`), with nothing else
/// queued on its pool, and the broadcast's jobs wait for nothing but each other.
///
/// Memory: the arena of `fill` (16 bytes a position of one chunk, `Shape`), not one list of the
/// whole query's words (8 bytes a word: 15 MiB more than the chains' fill at 8 threads for the
/// first form, v4-A15 s10).
pub(super) fn csr(
    parts: usize,
    words: &Words<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
) -> (super::x_dc_csr::Csr, bool) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if (2..=255).contains(&parts)
        && rayon::current_thread_index().is_some()
        && rayon::current_num_threads() == parts
    {
        let csr = csr_on_pool(
            parts,
            words,
            cells,
            pv_array,
            pv_array_bts,
            helper_array,
            chunk_positions(parts),
            x_offsets,
        );
        return (csr, true);
    }
    let _ = parts;
    let csr = super::x_bn_csr::serial(
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
/// positions (`chunk_positions`; the tests take others).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn csr_on_pool(
    parts: usize,
    words: &Words<'_>,
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
    // Once per table: the owners' cells cover the table (every count is written), and the
    // positions, which bound the words, fit the `u32` offsets.
    assert_eq!(own.iter().map(|&(c0, c1)| c1 - c0).sum::<usize>(), cells);
    assert!(cells.is_multiple_of(2048) && pv_array_bts >= 5 && helper_array.len() == cells >> 11);
    let (segments, total) = segments(words);
    assert!(
        total <= u32::MAX as usize,
        "BNCSR: more positions than u32 offsets"
    );
    // The reading-only `Fill`: `read::<false>` and `subrange` use only the words, the segments,
    // the owners and the arena.
    let shape = Shape::new(chunk.min(total), parts);
    let mut arena: Vec<Line> = Vec::with_capacity(lines_capacity(shape.lines()));
    let copies = Slots(Vec::new());
    let touched = Slots(Vec::new());
    let mut marks = PageBits::new(0, parts);
    #[cfg(test)]
    let checks = Checks::default();
    let fill = Fill {
        words,
        segments: &segments,
        owner: &owner,
        block_shift,
        parts,
        tables: SharedTables {
            hashtable: Shared::new(&mut []),
            next_pos: Shared::new(&mut []),
            pv_array: Shared::new(&mut []),
            pv_array_bts,
            helper_array: Shared::new(&mut []),
        },
        arena: SharedArena {
            base: arena.spare_capacity_mut().as_mut_ptr() as *mut Line,
            shape,
        },
        copies: &copies,
        faults: Faults {
            offset: 0,
            marks: marks.shared(),
            touched: &touched,
        },
        #[cfg(test)]
        checks: &checks,
    };
    // EXPERIMENT (LOSAT_X_BNCSRPACK): a packed head (`x_bn_csr::x_serial_packed`) on the pool, with
    // the two phases of the chunks' readings below. Owner `o` counts the words of its cells into
    // their bytes (`x_count8`; the words past 255 in its own map), then sets the presence bits of
    // its blocks (they own whole words of `bits`) and publishes its number of non-empty cells;
    // after a barrier it writes the ranks of its blocks from the number of the owners before it and
    // its cells' counts by rank (`x_rank_counts8`), and after another the sum and place phases
    // below run on the counts of its ranks instead of its cells.
    if super::x_bn_csr::x_packed(words.ascending_cells) {
        use super::x_dc_csr::{
            x_assemble_packed, x_count8, x_presence8, x_rank, x_rank_counts8, x_sum_packed, XPop,
        };
        // EXPERIMENT (LOSAT_X_BNPOPCNT): the form of the bit counts, read once per build.
        let pop = XPop::get();
        let blocks = cells >> 6;
        let mut bits = vec![0u64; 2 * blocks];
        // The non-empty cells are at most the cells and the words (bounded by the positions).
        let cap = cells.min(total) + 1;
        // EXPERIMENT (LOSAT_X_MBLUTKEEP): `x_offsets` is the batch loop's kept table or an empty vector.
        let mut starts = super::x_dc_csr::x_offsets_buffer(x_offsets, cap);
        // The counts take the first `cells` bytes of `qs`'s capacity (each owner zeroes its cells'
        // bytes first, and the place phase writes `qs` only after the barrier that follows the last
        // read of a count), so the build frees nothing in the middle of a batch (a freed block
        // there is trimmed and faulted in again, each trim interrupting every thread of the
        // process).
        let mut qs: Vec<u32> = Vec::with_capacity(total.max(cells.div_ceil(4)));
        let bb = Shared {
            ptr: bits.as_mut_ptr(),
            len: 2 * blocks,
        };
        let cb = Shared {
            ptr: qs.as_mut_ptr() as *mut u8,
            len: cells,
        };
        let sb = Shared {
            ptr: starts.as_mut_ptr(),
            len: cap,
        };
        let qb = Shared {
            ptr: qs.as_mut_ptr(),
            len: total,
        };
        let pv = Shared::new(pv_array);
        let helper = Shared::new(helper_array);
        let mut own_words = vec![0u64; parts];
        let own_w = Shared::new(&mut own_words);
        let mut own_present = vec![0u64; parts];
        let own_p = Shared::new(&mut own_present);
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
            let (bb, cb, sb, qb, pv, helper, own_w, own_p) =
                (bb, cb, sb, qb, pv, helper, own_w, own_p);
            let (c0, c1) = own[o];
            let (b0, b1) = (c0 >> 6, c1 >> 6);
            let mut more = std::collections::HashMap::new();
            // SAFETY (all of the job's accesses): as in the unpacked build below; owner `o` alone
            // writes its cells' bytes of the counts (before the barrier after which the place phase
            // writes `qs`), the words of its blocks in `bits` (presence, then
            // ranks, before the barrier after which every thread only reads `bits`), the counts of
            // its ranks (`starts[p0 .. p0 + its non-empty cells]`, disjoint from the other owners'),
            // its entries of the owners' words and non-empty cells (read by the later owners after a
            // barrier) and `qs` at its places.
            // SAFETY (`o`, `read`): `o = ctx.index() < parts` (the caller runs this only on a pool
            // of `parts` threads), and `span(c, o)` has at most `shape.half` positions (`read`).
            unsafe {
                std::ptr::write_bytes(cb.ptr.add(c0), 0, c1 - c0);
                let mut got = 0u64;
                for c in 0..chunks {
                    let (s0, s1) = span(c, o);
                    fill.read::<false, false>(o, s0, s1);
                    bar.wait();
                    for a in 0..parts {
                        let es = fill.subrange(a, o);
                        got += es.len() as u64;
                        for j in 0..es.len() {
                            if let Some(ahead) = es.get(j + 16) {
                                prefetch(cb.ptr.wrapping_add(ahead.value as usize));
                            }
                            x_count8(cb.ptr, &mut more, es[j].value as usize);
                        }
                    }
                    if c + 1 == chunks {
                        own_w.set(o, got);
                    }
                    bar.wait();
                }
                let mine = x_presence8(bb.ptr, cb.ptr, b0, b1, pop);
                own_p.set(o, mine as u64);
                bar.wait();
                let p0: u64 = (0..o).map(|p| own_p.get(p)).sum();
                let ranked = x_rank_counts8(bb.ptr, cb.ptr, &more, sb.ptr, b0, b1, p0 as usize);
                debug_assert_eq!(ranked, p0 as usize + mine);
                drop(more);
                bar.wait();
                let bits_r: &[u64] = std::slice::from_raw_parts(bb.ptr as *const u64, bb.len);
                let rank = |cell: u32| {
                    x_rank(bits_r, cell as usize, pop).expect("BNCSRPACK: a word of an empty cell")
                };
                let base: u64 = (0..o).map(|p| own_w.get(p)).sum();
                if c0 < c1 {
                    x_sum_packed(
                        pv.ptr,
                        pv.len,
                        helper.ptr,
                        bits_r,
                        sb.ptr,
                        base as u32,
                        b0,
                        b1,
                        pv_array_bts,
                        pop,
                    );
                }
                for c in 0..chunks {
                    let (s0, s1) = span(c, o);
                    fill.read::<false, false>(o, s0, s1);
                    bar.wait();
                    for a in 0..parts {
                        for e in fill.subrange(a, o) {
                            let r = rank(e.value);
                            let p = sb.get(r) - 1;
                            sb.set(r, p);
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
        let present = own_present.iter().sum::<u64>() as usize;
        assert!(
            n_words <= total as u64 && n_words <= u32::MAX as u64 && present < cap,
            "BNCSR: more words than positions"
        );
        // SAFETY: no job is running. Every element of `starts[..present]` was written (each owner
        // zeroed, counted, summed and placed its ranks, which cover `0 .. present`), `starts[present]`
        // (no cell's end) now, and the runs partition `0 .. n_words`, each place written once.
        unsafe {
            sb.set(present, n_words as u32);
            return x_assemble_packed(starts, bits, qs, present, n_words as usize);
        }
    }
    // `offsets` is written by its owners (their counts zeroed, then counted, summed and placed)
    // and by this thread (the elements before the first count, and `offsets[cells]` newest first);
    // `qs` gets one element per word (`total` positions bound them).
    let ascending = words.ascending_cells;
    let shift = super::x_bn_csr::count_shift(ascending);
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): `x_offsets` is the batch loop's kept table or an empty vector.
    let mut offsets = super::x_dc_csr::x_offsets_buffer(x_offsets, cells + 2);
    let mut qs: Vec<u32> = Vec::with_capacity(total);
    let ob = Shared {
        ptr: offsets.as_mut_ptr(),
        len: cells + 2,
    };
    // SAFETY: the first `shift` elements of the capacity, which no owner writes.
    unsafe { std::ptr::write_bytes(ob.ptr, 0, shift) };
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
        let cnt = ob.ptr.wrapping_add(shift);
        // SAFETY (all of the job's accesses): reader `o` uses only its own region and row of the
        // arena and owner `o` only the sub-ranges `(*, o)`, between barriers that separate the
        // readers from the owners; owner `o` alone writes its cells' counts and cursors
        // (`cnt[c0 .. c1]`, i.e. `offsets[c0 + shift .. c1 + shift]`, then `offsets[c0 + 1 .. c1 +
        // 1]` or `offsets[c0 .. c1]`), the presence-vector words and counters of its blocks, its
        // entry of the owners' words (read by the later owners after a barrier) and `qs` at its
        // places, which lie in `base .. base + its words` (disjoint from the other owners'); every
        // element read was written before by the same thread or before a barrier (the function's
        // documentation for the cursor of `c0`).
        // SAFETY (`o`, `read`): `o = ctx.index() < parts` (the caller runs this only on a pool of
        // `parts` threads), and `span(c, o)` has at most `shape.half` positions (`read`).
        unsafe {
            std::ptr::write_bytes(cnt.add(c0), 0, c1 - c0);
            let mut got = 0u64;
            for c in 0..chunks {
                let (s0, s1) = span(c, o);
                fill.read::<false, false>(o, s0, s1);
                bar.wait();
                for a in 0..parts {
                    let es = fill.subrange(a, o);
                    got += es.len() as u64;
                    for j in 0..es.len() {
                        if let Some(ahead) = es.get(j + 16) {
                            prefetch(cnt.wrapping_add(ahead.value as usize));
                        }
                        *cnt.add(es[j].value as usize) += 1;
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
                    cnt,
                    base as u32,
                    c0,
                    c1,
                    pv_array_bts,
                );
            }
            for c in 0..chunks {
                let (s0, s1) = span(c, o);
                fill.read::<false, false>(o, s0, s1);
                bar.wait();
                for a in 0..parts {
                    for e in fill.subrange(a, o) {
                        let cell = e.value as usize;
                        if ascending {
                            let cursor = ob.ptr.add(cell + 1);
                            qb.set(*cursor as usize, e.index);
                            *cursor += 1;
                        } else {
                            let p = ob.get(cell) - 1;
                            ob.set(cell, p);
                            qb.set(p as usize, e.index);
                        }
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
        "BNCSR: more words than positions"
    );
    // SAFETY: no job is running. Every element of `offsets[..cells + 1]` was written (the elements
    // before the counts here, the owners' counts covering the rest; newest first `offsets[cells]`,
    // no cell's end, now), and the runs partition `0 .. n_words`, each place written once.
    unsafe {
        if !ascending {
            ob.set(cells, n_words as u32);
        }
        assemble(offsets, qs, cells, n_words as usize)
    }
}

#[cfg(test)]
#[path = "x_mb_lookup_once_tests.rs"]
mod tests;
