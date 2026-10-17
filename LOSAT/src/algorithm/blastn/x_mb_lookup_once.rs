//! EXPERIMENT (LOSAT_X_MBLUTONCE; off unless set to any value, not in the strict set: its default
//! waits for the workstation's paired timings and max RSS, as LOSAT_X_MBLUTPAD's private copies,
//! which it always uses, raised max RSS there; `LOSAT_X_DEFAULT=0` turns `LOSAT_X_PAIRPAR` off and
//! with it this fill): `LOSAT_X_PAIRPAR`'s pool fill of the megablast lookup table of
//! `build_mb_lookup` (`x_mb_lookup_par.rs`) with each query position read once. In the fill of
//! `x_mb_lookup_par` every worker reads every word of the query and applies those of its own
//! cells, so each worker pays the whole read (0.65 s of each worker's time at 16 threads for
//! megablast of 31,168 D. melanogaster cDNAs on the workstation). Here the words are read by spans
//! and handed to the owners of their cells, as LOSAT_X_DISCLUTPAR does for the discontiguous
//! tables (`x_disc_lookup_par.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1033,1043-1045,1081-1091
//! ```c
//!       for (index = from; index <= last_offset; index++) {
//!          val = *++seq;
//! ...
//!          ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
//!          if (seq < pos)
//!             continue;
//! ...
//!          if (mb_lt->hashtable[ecode] == 0) {
//! ...
//!             PV_SET(pv_array, ecode, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
//!          mb_lt->hashtable[ecode] = index;
//! ```
//!
//! The plan. The cells are cut into blocks of `2^block_shift` cells dealt round robin to P owners
//! (`x_mb_lookup_par::plan`, the same plan as the fill there); P is the size of the search pool.
//! The positions the loop reads, all positions of the unmasked ranges of all queries in the loop's
//! order (`segments`: the loop's `build_unmasked_ranges` per query, without the queries shorter
//! than `lut_word_length` and the ranges shorter than `word_length`, which the loop skips), form
//! one stream, cut into chunks of `CHUNK_POSITIONS` positions; each chunk is cut into P spans of
//! equal length in stream order, one per reader. Per chunk:
//! - A (read, job = reader `a`): reader `a` reads the words that end in its span with the loop's
//!   statements (the same ambiguity resets, `current_kmer`, `valid_bases < lut_word_length` and
//!   query offset `q_off_1`), entering a range up to `lut_word_length - 1` positions before its
//!   span (never before the range's start). It writes each word `(q_off_1, cell)` to its scratch
//!   region in reading order, counting the words per owner, and then copies them, stably, into its
//!   sorted region grouped by owner: sub-range `(a, o)` holds the words of owner `o`'s cells in
//!   reading order, and reader `a`'s row of the arena holds the sub-ranges' bounds.
//! - B (update, job = owner `o`): owner `o` takes the sub-ranges `(0, o)`, `(1, o)`, ...,
//!   `(P-1, o)` and makes the loop's update for each word: `head = hashtable[cell]`; the presence
//!   bit in its private copy of `pv_array` when `head == 0`, otherwise a saturating `+1` of the
//!   cell's counter in its private copy of `helper_array` (the copies of LOSAT_X_MBLUTPAD,
//!   `x_mb_lookup_par::Copies`, zeroed by the owner in its first job); then `hashtable[cell] =
//!   q_off_1`. It keeps `head` in the entry. On x86_64 it prefetches the cell 16 entries ahead.
//! - C (link, job = reader `a`): reader `a` writes `next_pos[q_off_1] = head` for the entries of
//!   its sorted region.
//!
//! C of one chunk and A of the next run in the same job (reader `a` links its region, then reads
//! the next chunk into it), so a chunk takes two `rayon::scope`s; the last scope links the last
//! chunk and has job `o` add owner `o`'s copies into `helper_array` and `pv_array` for its own
//! blocks (`x_mb_lookup_par::flush_private`). `longest_chain` is then computed by the original code.
//!
//! Why the tables are the loop's. The words of a span are the loop's words ending there, with the
//! loop's cell and offset: the k-mer of a position depends only on the last `lut_word_length`
//! bases (the mask), and a word is taken iff its last `lut_word_length` positions lie in the range
//! and are unambiguous; a reader that starts `lut_word_length - 1` positions before its span (or
//! at the range's start) sees those positions, so its `valid_bases` reaches `lut_word_length` at the
//! same positions of the span as the loop's. Each word is read by the one reader whose span holds
//! its last position. A block holds whole counters and whole presence-vector words (`2^block_shift`
//! is a multiple of 2048 and of `2^pv_array_bts`), so owner `o` is the only job that touches the
//! cells of its blocks, their counters and their presence-vector words. B takes the chunks in
//! order, the readers in order (their spans follow each other in the stream) and each sub-range in
//! reading order, and the copy into the sorted region is stable, so every cell receives the loop's
//! words in the loop's order. That fixes the cell's final head, each word's head before its update
//! (exactly the value the loop writes to `next_pos[q_off_1]`, which C writes once per word: the
//! offsets are unique), and the increments of each counter and the bits of each presence-vector
//! word, which the copies hold starting from 0; the flush gives `min(s + n, u32::MAX)` and
//! `s | bits` from the start values `s`, the loop's saturating increments and ORed bits (0 in every
//! build; `x_mb_lookup_par`'s "Why the counters and presence bits are the loop's"). So
//! `hashtable`, `next_pos`, `pv_array` and `helper_array`, and with them `longest_chain`, equal the
//! loop's. Each job of a phase touches only its own elements (the table above in the design note
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
//! sorted, 8-byte entries), with `c` the smaller of `CHUNK_POSITIONS` and the positions the loop
//! reads, and P rows of `ceil((P + 1) / 16)` lines (the sub-range bounds), each region and row on
//! lines of its own: at most `16 x (c + 8 P)` bytes of regions, 8.0 MiB at 16 and at 64 threads and
//! 8.3 MiB with the rows at 255 for a query of 2^19 positions or more; it does not grow. On Linux with glibc its capacity is at least 32 MiB + 4 KiB
//! (`x_mb_lookup_par::lines_capacity`, the rule of `x_frame_capacity` in
//! `tblastn/x_subject_translate.rs`), so glibc maps it and unmaps it when it is freed at the end of
//! the fill without moving its mmap threshold, and only the touched 8 MiB are resident; elsewhere
//! (and on wasm32, whose linear memory never shrinks) it has its exact size. The private copies
//! are those of LOSAT_X_MBLUTPAD (`x_mb_lookup_par::Copies`, kept from fill to fill on glibc). The
//! per-range list of the loop's re-threading of `ascending_cells` (8 bytes a word) is not needed.
//!
//! The pool. As LOSAT_X_PAIRPAR's fill (`x_mb_lookup_par::parts`): a LOSAT search pool of two
//! threads or more (a LOSAT_X_QBATCHPAR slot of `m >= 2` threads fills on its own pool with P = m),
//! no `LOSAT_X_PAIRPARSHADOW`, no database word counts, no `BLEMIR_DEBUG`. Without the `parallel`
//! feature, or on wasm32 without `wasm-threads`, the loop runs; with `wasm-threads` the scopes run
//! on the search pool, which this module never builds (memory losat-wasm-threads-main-only).
//!
//! Shadow and timing: with `LOSAT_X_PAIRPARSHADOW` the second build of `x_mb_lookup_par`'s shadow
//! is this fill (`[X_SHADOW] MBLUTONCE tables= cells= diffs=0 identical` at exit); with LOSAT_TIMING
//! `[TIMING] x_mb_lookup_once: fills= chunks= read_wall= read_max= read_sum= read_late= update_...
//! last_... arena=`, per phase summed over the chunks of the search's fills (not the shadow's): the
//! calling worker's wall time from before the phase's first job to the return of its last, the
//! longest job, the jobs' summed time and the latest start of a job; `read` is the scope of A (with
//! C of the chunk before), `update` that of B and `last` the scope that links the last chunk and
//! flushes the copies; `arena` the largest arena of one fill (MiB).

use super::x_mb_lookup_par::{
    flush_private, lines_capacity, plan, timing_on, Copies, Line, Private, Shared, SharedTables,
    Tables, Unready, Words, K_COMPRESSION_FACTOR_SHIFT,
};
use super::{build_unmasked_ranges, PvArrayType, BLAST2NA_MASK, PV_ARRAY_MASK};
use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicU64, Ordering};

/// No NCBI counterpart (scheduling only): the stream positions per chunk (each chunk is read,
/// updated and linked before the next; the arena holds the words of one chunk), as LOSAT_X_DISCLUTPAR.
const CHUNK_POSITIONS: usize = 1 << 19;

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

// No NCBI counterpart (LOSAT_TIMING only): the phases' split, the fills and chunks of the search's
// read-once fills, and the largest arena of one fill (bytes).
static PHASE_TIMES: [PhaseTimes; 3] = [const { PhaseTimes::new() }; 3];
const PHASE_NAMES: [&str; 3] = ["read", "update", "last"];
static FILLS: AtomicU64 = AtomicU64::new(0);
static CHUNKS: AtomicU64 = AtomicU64::new(0);
static ARENA_PEAK_BYTES: AtomicU64 = AtomicU64::new(0);

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
        " arena={:.1}MiB",
        load(&ARENA_PEAK_BYTES) as f64 / 1048576.0
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
/// ```c
///       /* case of unmasked region >= kLutWordLength but < full_word_size,
///          so no hits should be generated. */
///       if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
///          continue;
/// ```
/// The ranges the loop of `build_mb_lookup` reads, in its order (the same `build_unmasked_ranges`
/// per query; the queries shorter than `lut_word_length` and the ranges shorter than `word_length`
/// skipped as there), and the number of positions it reads.
fn segments<'a>(words: &Words<'a>) -> (Vec<Segment<'a>>, usize) {
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
const _: () = assert!(ENTRY_BYTES == 8 && LINE_BYTES % ENTRY_BYTES == 0);
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
/// `pv_array` (LOSAT_X_MBLUTPAD's), zeroed by the owner in its first phase-B job.
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
    #[cfg(test)]
    checks: &'a Checks,
}

impl Fill<'_, '_, '_> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1045
    /// ```c
    ///       for (index = from; index <= last_offset; index++) {
    ///          val = *++seq;
    ///          /* if an ambiguity is encountered, do not add
    ///             any words that would contain it */
    ///          if ((val & BLAST2NA_MASK) != 0) {
    ///             ecode = 0;
    ///             pos = seq + kLutWordLength;
    ///             continue;
    ///          }
    ///
    ///          /* get next base */
    ///          ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
    ///          if (seq < pos)
    ///             continue;
    /// ```
    /// Phase A, reader `a`: the words that end at the stream positions `[s0, s1)`, read with the
    /// statements of `build_mb_lookup`'s loop (the same resets, k-mer and offsets), into its scratch
    /// region in reading order, then copied stably into its sorted region grouped by the owner of
    /// the cell; its row gets the start of each owner's sub-range and the number of words. A range
    /// is entered up to `lut_word_length - 1` positions before `s0` (not before its start), so that
    /// the k-mer and `valid_bases >= lut_word_length` at the positions from `s0` on are the loop's.
    ///
    /// SAFETY: no other job of the phase uses reader `a`'s region or row.
    unsafe fn read(&self, a: usize, s0: usize, s1: usize) {
        let words = self.words;
        let lut_word_length = words.lut_word_length;
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
            let mut current_kmer: u64 = 0;
            let mut valid_bases = 0usize;
            for pos in warm_start..end_pos {
                let base = seg.seq[pos];
                if (base & BLAST2NA_MASK) != 0 {
                    current_kmer = 0;
                    valid_bases = 0;
                    continue;
                }
                current_kmer = ((current_kmer << 2) | base as u64) & words.kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length || pos < first_pos {
                    continue;
                }
                let cell = current_kmer as usize;
                let q_off_1 = (seg.query_offset + (pos + 1 - lut_word_length) + 1) as u32;
                debug_assert!(n < self.arena.shape.half);
                // SAFETY: reader `a`'s scratch region has `half` entries, at least the positions
                // of its span (`Shape::new`), and one word per position at most.
                unsafe {
                    scratch.add(n).write(Entry {
                        index: q_off_1,
                        value: cell as u32,
                    })
                };
                n += 1;
                counts[self.owner[cell >> self.block_shift] as usize] += 1;
            }
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

    /// Phase B, owner `o`: the updates of the words of its cells, reader by reader, each sub-range
    /// in reading order (both reversed for `ascending_cells`), counting into its own copies; each
    /// entry keeps the cell's head before the update.
    ///
    /// SAFETY: no other job of the phase uses the sub-ranges `(*, o)`, the cells of `o`'s blocks or
    /// `o`'s copies.
    unsafe fn update(&self, o: usize) {
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
                    update_backward(self.tables, x, entries);
                } else {
                    update_forward(self.tables, x, entries);
                }
            }
        }
    }

    /// Phase C, reader `a`: `next_pos[index]` of the words it read, from their kept heads.
    ///
    /// SAFETY: no other job of the phase uses reader `a`'s region; the indices are the words' own.
    unsafe fn link(&self, a: usize) {
        // SAFETY: reader `a` wrote its row in phase A.
        let n = unsafe { self.arena.row(a).add(self.parts).read() } as usize;
        // SAFETY: the first `n` entries of the sorted region were written in phase A.
        let entries = unsafe { std::slice::from_raw_parts(self.arena.sorted(a), n) };
        for e in entries {
            // SAFETY: each index is written once, by the reader of its word.
            unsafe { self.tables.next_pos.set(e.index as usize, e.value) };
        }
    }

    /// After the last phase B, owner `o`: its copies added into `helper_array` and ORed into
    /// `pv_array` for its own blocks (`x_mb_lookup_par::flush_private`).
    ///
    /// SAFETY: no other job of the phase uses `o`'s copies or the counters and presence-vector words
    /// of `o`'s blocks.
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
/// ```c
///          if (mb_lt->hashtable[ecode] == 0) {
///             PV_SET(pv_array, ecode, pv_array_bts);
///          }
///          else {
///             helper_array[ecode/kCompressionFactor]++;
///          }
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
///          mb_lt->hashtable[ecode] = index;
/// ```
/// The update of one word, as in the loop of `build_mb_lookup`, with the counter and the presence
/// bit in the owner's private copies; `next_pos[index]` is written in phase C from the head kept in
/// the entry.
///
/// SAFETY: the cell lies in a block of the calling owner, and `x` is its copies.
#[inline(always)]
unsafe fn update_one(t: SharedTables, x: &mut Private<'_>, e: &mut Entry) {
    let cell = e.value as usize;
    // SAFETY: the caller's guarantee (the cell < the table size: a masked k-mer).
    unsafe {
        let head = t.hashtable.get(cell);
        if head == 0 {
            let w = cell >> t.pv_array_bts;
            if w < t.pv_array.len {
                *x.pv.slot(w) |= (1 as PvArrayType) << (cell & PV_ARRAY_MASK);
            }
        } else {
            let c = x.counters.slot(cell >> K_COMPRESSION_FACTOR_SHIFT);
            *c = c.saturating_add(1);
        }
        t.hashtable.set(cell, e.index);
        e.value = head;
    }
}

/// The updates of a sub-range's entries in order (phase B of a descending table).
///
/// SAFETY: as `update_one`, for every entry.
#[inline]
unsafe fn update_forward(t: SharedTables, x: &mut Private<'_>, entries: &mut [Entry]) {
    for j in 0..entries.len() {
        #[cfg(target_arch = "x86_64")]
        if let Some(ahead) = entries.get(j + X_AHEAD) {
            prefetch(t.hashtable, ahead.value as usize);
        }
        // SAFETY: the caller's guarantee.
        unsafe { update_one(t, x, &mut entries[j]) };
    }
}

/// The updates of a sub-range's entries back to front (phase B of `ascending_cells`).
///
/// SAFETY: as `update_one`, for every entry.
#[inline]
unsafe fn update_backward(t: SharedTables, x: &mut Private<'_>, entries: &mut [Entry]) {
    for j in (0..entries.len()).rev() {
        #[cfg(target_arch = "x86_64")]
        if j >= X_AHEAD {
            prefetch(t.hashtable, entries[j - X_AHEAD].value as usize);
        }
        // SAFETY: the caller's guarantee.
        unsafe { update_one(t, x, &mut entries[j]) };
    }
}

/// A prefetch of cell `cell` of `hashtable` (x86_64).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn prefetch(hashtable: Shared<u32>, cell: usize) {
    debug_assert!(cell < hashtable.len);
    // SAFETY: a cell index < the table size; a prefetch does not fault.
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            hashtable.ptr.add(cell) as *const i8,
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

/// No NCBI counterpart (scheduling only): what one fill did: its chunks and the lines of its arena;
/// in tests also the owners whose copies held a value before the flush, the non-empty sub-ranges
/// phase B took and how many of them were in increasing query offset (reading) order.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct OnceStats {
    pub(super) chunks: usize,
    pub(super) arena_lines: usize,
    #[cfg(test)]
    pub(super) owners_with_values: usize,
    #[cfg(test)]
    subranges: usize,
    #[cfg(test)]
    subranges_in_order: usize,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1033,1081-1091
/// ```c
///       for (index = from; index <= last_offset; index++) {
///          val = *++seq;
/// ...
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
///          mb_lt->hashtable[ecode] = index;
/// ```
/// Fills the tables as the word loop of `build_mb_lookup` does (with its re-threading for
/// `ascending_cells`), each query position read once, with `parts` workers on the search pool (or
/// one after another, `in_pool == false`, for the shadow without a pool); the module documentation
/// gives the phases and why every element of the four arrays is the loop's. `timed`: a build of
/// the search (LOSAT_TIMING's `[TIMING] x_mb_lookup_once` counts it), not the shadow's.
pub(super) fn fill(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    timed: bool,
) -> OnceStats {
    fill_chunked(parts, words, tables, in_pool, timed, CHUNK_POSITIONS)
}

/// `fill` with the chunk size (the tests take others).
fn fill_chunked(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    timed: bool,
    chunk: usize,
) -> OnceStats {
    let table_size = tables.hashtable.len();
    let (block_shift, owner) = plan(table_size, tables.pv_array_bts, parts);
    debug_assert!(
        block_shift >= K_COMPRESSION_FACTOR_SHIFT && block_shift >= tables.pv_array_bts,
        "a block must own whole counters and presence-vector words"
    );
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
    let mut copies = Copies::new(true, parts, shared.helper_array.len, shared.pv_array.len);
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
        #[cfg(test)]
        checks: &checks,
    };
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
    let span = |(start, end): (usize, usize), a: usize| {
        (
            start + (end - start) * a / parts,
            start + (end - start) * (a + 1) / parts,
        )
    };
    let mut clock = (timed && in_pool && timing_on()).then(|| JobClock::new(parts));
    // SAFETY (all phases): job `i` of a phase uses only reader `i`'s region and row, owner `i`'s
    // sub-ranges, cells and copies, and the `next_pos` entries of reader `i`'s words (see
    // `SharedArena`, `Fill::read`, `Fill::update`, `Fill::link`, `Fill::flush`); `run_jobs` returns
    // after every job of the phase has returned.
    for k in 0..chunks {
        let current = bounds(k);
        run_phase(&mut clock, 0, parts, in_pool, |a| {
            if k > 0 {
                // C of the chunk before.
                unsafe { fill.link(a) };
            }
            let (s0, s1) = span(current, a);
            unsafe { fill.read(a, s0, s1) };
        });
        run_phase(&mut clock, 1, parts, in_pool, |o| unsafe { fill.update(o) });
    }
    run_phase(&mut clock, 2, parts, in_pool, |i| unsafe {
        fill.link(i);
        fill.flush(i);
    });
    stats.chunks = chunks;
    stats.arena_lines = shape.lines();
    #[cfg(test)]
    {
        stats.owners_with_values = checks.owners_with_values.into_inner();
        stats.subranges = checks.subranges.into_inner();
        stats.subranges_in_order = checks.subranges_in_order.into_inner();
    }
    if clock.is_some() {
        FILLS.fetch_add(1, Ordering::Relaxed);
        CHUNKS.fetch_add(chunks as u64, Ordering::Relaxed);
        ARENA_PEAK_BYTES.fetch_max((shape.lines() * LINE_BYTES) as u64, Ordering::Relaxed);
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::super::x_mb_lookup_par::tests::{random_case, repeat_case, Lcg};
    use super::super::{build_two_stage_lookup, MaskedInterval};
    use super::*;

    /// The tables of one fill.
    #[derive(Debug, PartialEq, Eq)]
    struct Built {
        hashtable: Vec<u32>,
        next_pos: Vec<u32>,
        pv_array: Vec<PvArrayType>,
        helper_array: Vec<u32>,
    }

    /// One case: the queries, their offsets and masks, `word_length`, `lut_word_length` and
    /// `ascending_cells`.
    struct Case {
        queries: Vec<Vec<u8>>,
        offsets: Vec<i32>,
        masks: Vec<Vec<MaskedInterval>>,
        word_length: usize,
        lut: usize,
        ascending: bool,
    }

    impl Case {
        fn random(rng: &mut Lcg, ascending: bool) -> Self {
            let (queries, offsets, masks, word_length, lut) = random_case(rng);
            Self {
                queries,
                offsets,
                masks,
                word_length,
                lut,
                ascending,
            }
        }

        fn one(query: Vec<u8>, word_length: usize, lut: usize, ascending: bool) -> Self {
            Self {
                queries: vec![query],
                offsets: vec![0],
                masks: vec![Vec::new()],
                word_length,
                lut,
                ascending,
            }
        }

        fn refs(&self) -> Vec<&[u8]> {
            self.queries.iter().map(Vec::as_slice).collect()
        }

        fn words<'a>(&'a self, refs: &'a [&'a [u8]]) -> Words<'a> {
            Words {
                queries_blastna: refs,
                query_offsets: &self.offsets,
                query_masks: &self.masks,
                word_length: self.word_length,
                lut_word_length: self.lut,
                kmer_mask: (1u64 << (2 * self.lut)) - 1,
                ascending_cells: self.ascending,
            }
        }

        /// `build_two_stage_lookup`'s table (the loop), with its counters (from 0) recomputed by
        /// the loop's statements in one part of `x_mb_lookup_par`'s fill (the tables there are
        /// checked against the loop) and its `pv_array_bts`.
        fn reference(&self, approx: usize) -> (Built, usize, usize) {
            let refs = self.refs();
            let r = build_two_stage_lookup(
                &refs,
                &self.offsets,
                self.word_length,
                self.lut,
                &self.masks,
                None,
                0,
                approx,
                self.ascending,
            );
            let mb = &r.mb_lookup;
            let helper = self.helper_of(mb.hashtable.len(), mb.pv_array.len(), mb.pv_array_bts);
            (
                Built {
                    hashtable: mb.hashtable.clone(),
                    next_pos: mb.next_pos.clone(),
                    pv_array: mb.pv_array.clone(),
                    helper_array: helper,
                },
                mb.pv_array_bts,
                mb.longest_chain,
            )
        }

        /// The loop's counters: the unpadded fill of `x_mb_lookup_par` in one part, which applies
        /// the loop's statements in reading order.
        fn helper_of(&self, cells: usize, pv_len: usize, bts: usize) -> Vec<u32> {
            let refs = self.refs();
            let words = self.words(&refs);
            let mut hashtable = vec![0u32; cells];
            let mut next_pos = vec![0u32; self.next_len()];
            let mut pv_array = vec![0 as PvArrayType; pv_len];
            let mut helper_array = vec![0u32; (cells >> 11).max(1)];
            super::super::x_mb_lookup_par::tests::fill_old(
                1,
                &words,
                Tables {
                    hashtable: &mut hashtable,
                    next_pos: &mut next_pos,
                    pv_array: &mut pv_array,
                    pv_array_bts: bts,
                    helper_array: &mut helper_array,
                },
            );
            helper_array
        }

        fn next_len(&self) -> usize {
            self.queries
                .iter()
                .zip(&self.offsets)
                .map(|(q, &o)| o.max(0) as usize + q.len())
                .max()
                .unwrap_or(0)
                + 1
        }

        /// The read-once fill, from start values `helper0` / `pv0`.
        fn once(
            &self,
            shape: (usize, usize),
            parts: usize,
            in_pool: bool,
            chunk: usize,
            helper0: &[u32],
            pv0: &[PvArrayType],
        ) -> (Built, OnceStats) {
            let (cells, bts) = shape;
            let refs = self.refs();
            let words = self.words(&refs);
            let mut b = Built {
                hashtable: vec![0u32; cells],
                next_pos: vec![0u32; self.next_len()],
                pv_array: pv0.to_vec(),
                helper_array: helper0.to_vec(),
            };
            let stats = fill_chunked(
                parts,
                &words,
                Tables {
                    hashtable: &mut b.hashtable,
                    next_pos: &mut b.next_pos,
                    pv_array: &mut b.pv_array,
                    pv_array_bts: bts,
                    helper_array: &mut b.helper_array,
                },
                in_pool,
                false,
                chunk,
            );
            (b, stats)
        }

        /// The read-once fill from zero start values.
        fn once0(
            &self,
            reference: &Built,
            bts: usize,
            parts: usize,
            in_pool: bool,
            chunk: usize,
        ) -> (Built, OnceStats) {
            self.once(
                (reference.hashtable.len(), bts),
                parts,
                in_pool,
                chunk,
                &vec![0u32; reference.helper_array.len()],
                &vec![0 as PvArrayType; reference.pv_array.len()],
            )
        }

        /// The number of stream positions (the loop's reads).
        fn positions(&self) -> usize {
            let refs = self.refs();
            segments(&self.words(&refs)).1
        }
    }

    fn longest(helper: &[u32]) -> usize {
        helper.iter().fold(2usize, |m, &v| m.max(v as usize))
    }

    fn assert_built(reference: &Built, chain: usize, b: &Built, what: &str) {
        assert_eq!(b.hashtable, reference.hashtable, "{what}: hashtable");
        assert_eq!(b.next_pos, reference.next_pos, "{what}: next_pos");
        assert_eq!(b.pv_array, reference.pv_array, "{what}: pv_array");
        assert_eq!(
            b.helper_array, reference.helper_array,
            "{what}: helper_array"
        );
        assert_eq!(longest(&b.helper_array), chain, "{what}: longest_chain");
    }

    // Random cases (ambiguity codes, masks, repeats, `lut_word_length` 4 / 8 / 9 / 11 / 12,
    // `word_length` above it, 1-3 queries, both `ascending_cells`), off the pool with 1..8 parts
    // and chunks of 1 (small cases), 7, 1000 positions and the default: every element of the four
    // arrays and the longest chain equal the loop's. The sub-ranges are in reading order, the
    // chunks are `ceil(positions / chunk)`.
    #[test]
    fn once_fill_builds_the_loops_tables() {
        let mut rng = Lcg(0x0E_CE_5EED);
        let mut multi_chunk = 0;
        for case in 0..40 {
            let c = Case::random(&mut rng, case % 2 == 1);
            let (reference, bts, chain) = c.reference(1_000_000);
            assert_eq!(longest(&reference.helper_array), chain, "case {case}");
            let positions = c.positions();
            for parts in 1..=8usize {
                for chunk in [1usize, 7, 1000, CHUNK_POSITIONS] {
                    if chunk == 1 && (positions > 3000 || parts > 3) {
                        continue;
                    }
                    let (b, stats) = c.once0(&reference, bts, parts, false, chunk);
                    let what = format!("case {case} parts {parts} chunk {chunk}");
                    assert_built(&reference, chain, &b, &what);
                    assert_eq!(stats.chunks, positions.div_ceil(chunk), "{what}");
                    if positions > 0 {
                        assert_eq!(
                            stats.arena_lines,
                            Shape::new(chunk.min(positions), parts.clamp(1, 255)).lines(),
                            "{what}"
                        );
                    }
                    assert_eq!(stats.subranges, stats.subranges_in_order, "{what}");
                    multi_chunk += usize::from(stats.chunks >= 2 && stats.subranges > 0);
                }
            }
        }
        assert!(multi_chunk > 100, "{multi_chunk}");
    }

    // Edge cases: ranges shorter than `word_length` (skipped as the loop does), queries shorter
    // than `lut_word_length`, all-N queries, no words at all, a word ending exactly at a span or
    // chunk boundary after an ambiguity, and queries placed at offsets with gaps.
    #[test]
    fn short_ranges_short_queries_and_no_words() {
        let n = |len: usize| vec![14u8; len];
        let acgt = |len: usize, seed: u64| {
            let mut rng = Lcg(seed);
            (0..len).map(|_| rng.below(4) as u8).collect::<Vec<u8>>()
        };
        let mut cases = vec![
            // all N
            Case::one(n(500), 12, 8, false),
            // shorter than lut
            Case::one(acgt(6, 1), 8, 8, true),
            // empty
            Case::one(Vec::new(), 12, 8, false),
        ];
        // ranges shorter than word_length between masks; one N every 9 bases
        let q = acgt(2000, 2);
        let mut masks = Vec::new();
        for start in (0..2000).step_by(40) {
            masks.push(MaskedInterval {
                start: start + 15 + (start % 7),
                end: start + 40,
            });
        }
        cases.push(Case {
            queries: vec![q.clone(), acgt(5, 3), n(30), q.clone()],
            offsets: vec![0, 2001, 2010, 2100],
            masks: vec![masks.clone(), Vec::new(), Vec::new(), Vec::new()],
            word_length: 16,
            lut: 11,
            ascending: false,
        });
        let mut q9 = acgt(3000, 4);
        for i in (0..3000).step_by(9) {
            q9[i] = 14;
        }
        cases.push(Case::one(q9, 8, 8, true));
        cases.push(Case {
            queries: vec![q.clone(), q],
            offsets: vec![0, 2001],
            masks: vec![masks.clone(), masks],
            word_length: 12,
            lut: 12,
            ascending: true,
        });
        for (i, c) in cases.iter().enumerate() {
            let (reference, bts, chain) = c.reference(30_000);
            for parts in [1usize, 2, 3, 5, 8] {
                for chunk in [1usize, 2, 7, 64, CHUNK_POSITIONS] {
                    let (b, _) = c.once0(&reference, bts, parts, false, chunk);
                    assert_built(
                        &reference,
                        chain,
                        &b,
                        &format!("edge {i} parts {parts} chunk {chunk}"),
                    );
                }
            }
        }
    }

    // On real pools of 2, 3, 5, 16 and 64 threads (the fill inside the pool, as in a search):
    // random cases and long unmasked repeats at `lut_word_length` 11 and 12 (many words per cell,
    // counters that grow, presence-vector words shared between blocks), default chunks and chunks
    // of 1000 positions; at least two owners' copies held values.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn once_fill_on_pools_builds_the_loops_tables() {
        let mut rng = Lcg(0x00C0_FFEE_0001);
        let mut cases = Vec::new();
        for case in 0..6 {
            cases.push(Case::random(&mut rng, case % 2 == 0));
        }
        for (lut, asc) in [(11usize, false), (12, true), (12, false)] {
            cases.push(Case::one(repeat_case(&mut rng, 40_000), lut + 5, lut, asc));
        }
        for (i, c) in cases.iter().enumerate() {
            let (reference, bts, chain) = c.reference(1_000_000);
            for threads in [2usize, 3, 5, 16, 64] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                for chunk in [1000usize, CHUNK_POSITIONS] {
                    let (b, stats) =
                        pool.install(|| c.once0(&reference, bts, threads, true, chunk));
                    let what = format!("case {i} pool {threads} chunk {chunk}");
                    assert_built(&reference, chain, &b, &what);
                    assert_eq!(stats.subranges, stats.subranges_in_order, "{what}");
                    if i >= 6 && threads >= 3 {
                        assert!(stats.owners_with_values >= 2, "{what}: {stats:?}");
                    }
                }
            }
        }
    }

    // Many parts: 64 / 254 / 255 / 256 / 300 off the pool (the last two clamp to 255) and 255 on
    // a 64-thread pool, at `lut_word_length` 8 (32 blocks: most owners have none), 11 and 12.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn once_fill_with_many_parts_builds_the_loops_tables() {
        let mut rng = Lcg(0x0255_0255);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(64)
            .build()
            .unwrap();
        for (lut, asc, approx) in [
            (8usize, true, 30_000usize),
            (11, false, 1_000_000),
            (12, true, 400_000),
        ] {
            let c = Case::one(repeat_case(&mut rng, 25_000), lut + 4, lut, asc);
            let (reference, bts, chain) = c.reference(approx);
            for parts in [64usize, 254, 255, 256, 300] {
                let (b, _) = c.once0(&reference, bts, parts, false, CHUNK_POSITIONS);
                assert_built(&reference, chain, &b, &format!("lut {lut} parts {parts}"));
            }
            for chunk in [999usize, CHUNK_POSITIONS] {
                let (b, _) = pool.install(|| c.once0(&reference, bts, 255, true, chunk));
                assert_built(
                    &reference,
                    chain,
                    &b,
                    &format!("lut {lut} pool 64 parts 255 chunk {chunk}"),
                );
            }
        }
    }

    // From start values other than 0 (production starts at 0): counters at u32::MAX - k and
    // others, presence-vector words with bits set; the flush's `saturating_add` and OR give the
    // loop's saturating increments and kept bits (the loop's values: `x_mb_lookup_par`'s unpadded
    // fill in one part from the same start values).
    #[test]
    fn once_flush_matches_the_loop_from_any_start() {
        let mut rng = Lcg(0x5A7_0F1D_0002);
        for (lut, asc) in [(8usize, false), (11, true), (12, false)] {
            let c = Case::one(repeat_case(&mut rng, 30_000), lut + 4, lut, asc);
            let (reference, bts, _) = c.reference(30_000);
            let helper_len = reference.helper_array.len();
            let helper0: Vec<u32> = (0..helper_len)
                .map(|h| match h % 4 {
                    0 => u32::MAX,
                    1 => u32::MAX - (h as u32 % 7),
                    2 => 0,
                    _ => (h as u32).wrapping_mul(2_654_435_761),
                })
                .collect();
            let pv0: Vec<PvArrayType> = (0..reference.pv_array.len())
                .map(|w| {
                    if w % 3 == 0 {
                        (w as u32).wrapping_mul(40_503)
                    } else {
                        0
                    }
                })
                .collect();
            let saturated = (0..helper_len)
                .filter(|&h| {
                    helper0[h] != u32::MAX
                        && helper0[h].checked_add(reference.helper_array[h]).is_none()
                })
                .count();
            if lut >= 11 {
                assert!(saturated > 0, "lut {lut}: no counter saturates");
            }
            let expect_helper: Vec<u32> = (0..helper_len)
                .map(|h| helper0[h].saturating_add(reference.helper_array[h]))
                .collect();
            let expect_pv: Vec<PvArrayType> = (0..reference.pv_array.len())
                .map(|w| pv0[w] | reference.pv_array[w])
                .collect();
            for parts in [2usize, 3, 5, 16] {
                let (b, _) = c.once(
                    (reference.hashtable.len(), bts),
                    parts,
                    false,
                    CHUNK_POSITIONS,
                    &helper0,
                    &pv0,
                );
                let what = format!("lut {lut} parts {parts}");
                assert_eq!(b.helper_array, expect_helper, "{what}");
                assert_eq!(b.pv_array, expect_pv, "{what}");
                assert_eq!(b.hashtable, reference.hashtable, "{what}");
                assert_eq!(b.next_pos, reference.next_pos, "{what}");
            }
        }
    }

    // `next_pos` chain order: for megablast a cell's chain lists its offsets newest first, for
    // `ascending_cells` in reading order; and the longest chain of a long repeat (counters in the
    // thousands) is the loop's.
    #[test]
    fn chain_order_and_longest_chain() {
        // One 8-mer repeated: "ACGTACGT" period 4, so every 8-mer start at 4k is the same cell.
        let q: Vec<u8> = (0..4000).map(|i| [0u8, 1, 2, 3][i % 4]).collect();
        for asc in [false, true] {
            let c = Case::one(q.clone(), 8, 8, asc);
            let (reference, bts, chain) = c.reference(30_000);
            assert!(chain > 900, "{chain}");
            let (b, _) = c.once0(&reference, bts, 4, false, 100);
            assert_built(&reference, chain, &b, &format!("asc {asc}"));
            // The chain of the cell of the first word (offset 1): walk it.
            let cell = b.hashtable.iter().position(|&v| v != 0).unwrap();
            let mut walk = Vec::new();
            let mut at = b.hashtable[cell];
            while at != 0 {
                walk.push(at);
                at = b.next_pos[at as usize];
            }
            assert!(walk.len() > 900);
            if asc {
                assert!(walk.windows(2).all(|w| w[0] < w[1]), "ascending");
            } else {
                assert!(walk.windows(2).all(|w| w[0] > w[1]), "newest first");
            }
        }
    }

    // The arena: each reader's region and row on lines of their own, the regions within
    // `2 x (chunk + 8 P)` entries and the rows `P x ceil((P + 1) / 16)` lines; 8.0 MiB at 16 and 64
    // parts with the default chunks.
    #[test]
    fn arena_regions_and_rows_have_lines_of_their_own() {
        for parts in (1usize..=17).chain([31, 64, 128, 255]) {
            for chunk in [1usize, 7, 1000, CHUNK_POSITIONS] {
                let shape = Shape::new(chunk, parts);
                assert!(
                    shape.half >= chunk.div_ceil(parts),
                    "P {parts} chunk {chunk}"
                );
                assert_eq!(shape.half % LINE_ENTRIES, 0);
                assert_eq!(shape.region_lines * LINE_ENTRIES, 2 * shape.half);
                assert!(shape.row_lines * LINE_WORDS > parts);
                assert!(
                    parts * 2 * shape.half <= 2 * (chunk + 8 * parts),
                    "P {parts} chunk {chunk}"
                );
                let mut lines: Vec<Line> = Vec::with_capacity(shape.lines());
                let arena = SharedArena {
                    base: lines.spare_capacity_mut().as_mut_ptr() as *mut Line,
                    shape,
                };
                let base = arena.base as usize;
                assert_eq!(base % 64, 0);
                let end = base + shape.lines() * 64;
                let mut spans = Vec::new();
                for a in 0..parts {
                    let s = arena.scratch(a) as usize;
                    let t = arena.sorted(a) as usize;
                    let r = arena.row(a) as usize;
                    assert_eq!(s % 64, 0);
                    assert_eq!(t, s + shape.half * ENTRY_BYTES);
                    assert_eq!(r % 64, 0);
                    spans.push((s, s + shape.region_lines * 64));
                    spans.push((r, r + shape.row_lines * 64));
                }
                spans.sort();
                for w in spans.windows(2) {
                    assert!(w[0].1 <= w[1].0, "P {parts} chunk {chunk}: overlap");
                }
                assert!(spans.last().unwrap().1 <= end);
            }
        }
        for parts in [16usize, 64] {
            let bytes = Shape::new(CHUNK_POSITIONS, parts).lines() * 64;
            assert!(bytes <= (8 << 20) + 64 * 1024, "P {parts}: {bytes}");
        }
        assert!(Shape::new(CHUNK_POSITIONS, 255).lines() * 64 <= (8 << 20) + 400 * 1024);
    }
}
