//! EXPERIMENT (LOSAT_X_BPCLUT at more than one thread): the flat compressed lookup table of
//! `x_clut_direct.rs` built by the search's thread pool. The table is the serial builder's
//! (`x_try_build_flat_with_capacity`; for another matrix than BLOSUM62, `x_try_build_flat_of` with
//! the matrix's words, item SC2b), slot for slot, at every thread count.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 (`add_word_hit`, in
//! the header of `x_clut_direct.rs`)
//!
//! The serial builder makes the reference's appends (a query offset's cells, in the order of
//! `add_word_hits`, or of `LOSAT_X_CLUTPAIRS`' pair rows within an offset: `x_clut_direct.rs` header)
//! for the query offsets of its count ranges (`x_count_ranges`) in order. Here:
//! - chunks (`x_chunks`): each count range's query offsets, in order, cut into contiguous runs, a run
//!   given as ranges of lookup segments (`XChunk::pieces`: a range `(a, b + word_length - 1)` holds
//!   the offsets a..=b of its lookup segment, which `x_for_each_offset` visits as it visits them in
//!   the lookup segment, with the same `end <= query.len()` guard; not the pieces of appends below). The chunks in order hold the count ranges' offsets in order, so the
//!   chunks' appends one after the other are the serial count pass's appends in its order.
//! - stripes: the backbone cut into ranges of 2^shift cells. Each chunk sorts its appends by stripe,
//!   stably (counting sort), so a stripe's part of a chunk holds the chunk's appends to the stripe's
//!   cells in the chunk's order. A cell lies in one stripe; reading the stripe's parts of the chunks in
//!   chunk order gives the cell's appends in the serial order: its offsets o_0, ..., o_{n-1}, the
//!   append of rank k being the k-th one read.
//! - pieces (item SC2d-2): a chunk sorts its appends in pieces (`x_chunk_appends`: the appends of a
//!   run of its offsets, at most `X_PIECE_APPENDS` past the run's last offset's), each piece sorted
//!   by stripe as above; the pieces of a chunk in order hold the chunk's offsets in order, so reading
//!   a stripe's parts of the pieces in order (chunk by chunk) reads the stripe's part of each chunk in
//!   the chunk's order, and the above holds with pieces in place of chunks.
//! - each stripe's hits: the sum of the chunks' parts for it; a stripe's first hit is the hits of the
//!   stripes before it, so the stripes' slices of the hits are disjoint and in the serial layout's
//!   order (cells in index order).
//! - count (per stripe): each cell's count n (the serial count), and per count range the overflow
//!   cells of the reference's table of all batches that its appends take: an append that finds the
//!   count c takes one when `x_takes_overflow_cell(c)`, the serial count pass's statement, with c its
//!   rank, the count before it in the serial order. So the appends to a cell that holds n offsets at
//!   some point of the serial order have taken `x_overflow_cells_of(n)` cells by then (item SC9i of
//!   line NP: read from the counts, not counted append by append): a stripe's parts are read in
//!   chunk order, which holds the count ranges in order, and where they pass from one range to the
//!   next (`x_next_count_range`) the cells of the counts so far less those at the range before's
//!   start are that range's; the last range's are those of the final counts (at the layout) less
//!   those at its start. Summed over the stripes, a range's cells are those the serial pass counts
//!   for it (`st.cells - cells_before` at the range's end).
//! - layout (per stripe): the cells' starts are the counts' prefix sums from the stripe's first
//!   hit, so the serial starts; the longest chain and the count of non-empty cells are functions of
//!   the counts. The PV shift (`finalize`: `PV_ARRAY_BTS`, or more when at most 1 % of the cells
//!   are not empty) is `PV_ARRAY_BTS` whatever the counts when the base-2 logarithm (`ilog2`) of
//!   `backbone_size` divided by the product of 8 and `kTargetPVBytes` is 0 (word size 5): then each
//!   stripe sets the PV bits of its cells, in its own words (a stripe is 2^shift cells,
//!   shift >= `PV_ARRAY_BTS`); otherwise the PV array is made after the stripes from the starts (a
//!   cell is not empty when its start is below the next one), by `finalize`'s statements.
//! - fill (per stripe): the serial fill on the stripe's cells (each append, in the serial order, at
//!   its cell's cursor; the cursors start at the starts), then each cell of n > 5 offsets puts o_k at
//!   `x_emission_position(k, n)`: the serial fill's statements on the same cells.
//! - refusal: the serial builder tests a table when its cells pass the capacity (after a chunk of
//!   appends or at the table's end; the cells only grow, so a table passes the bound at some test
//!   exactly when it passes it at its end), once, and refuses when no segment comes before the
//!   table or the table's own appends fill the banks (`x_table_fills_banks`); the outcome does not
//!   depend on where the test is made. Here, after the stripes, each table's cells at its end are
//!   tested by the same rule, in the tables' order, so the build refuses exactly when the serial one
//!   does and names the same table (the first that refuses, `Refused(BanksFilled(i))`; the table
//!   made is dropped).
//! - a packed append holds the cell's index in its stripe (below 2^16) and the offset less the chunk's
//!   first offset (a chunk's offsets lie within 2^16 of it): nothing is lost.
//! - `TooLarge` (the reference builder runs) when the appends pass `x_max_appends`, as the serial
//!   builder's bound; the reference then decides the table or the refusal.
//!
//! The cuts (chunks, stripes) depend on the thread count; by the above the table does not. The chunks'
//! sizes are guided (item SC2d-2, `XChunkSizes::Guided`): a chunk takes the offsets left over
//! `X_GUIDED_SHARE` times the threads, at least `X_MIN_CHUNK_OFFSETS`, so that there are few chunks
//! (each chunk's sort and each of its stripes' parts costs a fixed amount) and the last ones are small
//! (the tasks end together); by the above any sizes give the same table.
//!
//! With `LOSAT_X_BPLKPFLAT` the same passes are scheduled otherwise and the table's vectors may come
//! from an earlier table (`XParallelCuts::passes`, `XParallelCuts::spare`); the table is the same:
//! - one-latch passes (`x_one_latch_pass`): the appends pass and the stripes pass each run in one
//!   `rayon::scope` with one spawned task per thread, instead of `par_iter`'s tree of joins; a task
//!   takes the index of the next item (a chunk, a stripe) from one counter and stores the item's
//!   result at that index. Each chunk's appends and each stripe's part are made by the same
//!   functions from the same inputs (the chunks, the stripes and the shift are chosen before the
//!   pass), the stripes write disjoint slices, and the results are read by index. The only state the
//!   chunks share, `appended` and `too_large`, gives `TooLarge` exactly when the appends of all
//!   chunks pass the bound, in any order (`x_chunk_appends`).
//! - recycled vectors (`XFlatSpare`, the previous table of a block search): `starts` and `hits` are
//!   that table's, cut or grown to their lengths (new entries 0), not zeroed. No entry of either is
//!   read before the stripes write it: each stripe's layout writes every entry of its slice of
//!   `starts`, the slices partition the cells, and `starts[backbone_size]` is set before the pass;
//!   the stripes' slices of `hits` partition it (their lengths are the stripes' hits, which sum to
//!   its length), each stripe's fill writes each slot of its slice once (a cell's cursor runs from
//!   its first slot over its count, the next cell's first slot), and the move to the scan order reads
//!   only slots the fill wrote. `pv` is zeroed as before (the stripes OR their bits into it). The
//!   build takes and sizes the vectors itself, before the stripes pass, never in a pass's tasks; a
//!   build that ends `TooLarge` drops the spare's vectors (the reference builder makes that table,
//!   and they are not kept beside it). The unit tests fill the spare's vectors with `u32::MAX` /
//!   `i32::MIN` first (a spare of garbage first gets room for a word-size-5 table), so that an
//!   entry no stripe writes differs from the fresh build's where the spare has the room.
//! - lent vectors (line NP, SC9d: `XParallelCuts::lent`, BLASTX's builds in parts with
//!   LOSAT_X_BXCLUTPOOL): `starts` and `hits` are those of a set the pool of BLASTX's tables lends,
//!   taken as the spare's above (no entry is read before the stripes write it), and `pv` that set's,
//!   cleared and zero-filled. The hits are of the table's offset type `H` (`u16` for BLASTX's query
//!   batches of at most `X_U16_LETTERS` letters, whose offsets are below 2^16): each stripe writes
//!   `H::from_offset` of the offset it wrote as `i32`, the same number.
//! - kept segments (item SC2d-2, `XSegmentPool`, the block search's): the appends pass takes a kept
//!   segment before it makes a new one, and the stripes pass keeps at most as many segments as threads
//!   for the next build, in place of what the pool still holds. A segment's values are not read: each piece writes every entry of its range
//!   before the stripes pass reads it.
//! - extra tasks (`XExtraTasks`: the contexts' Karlin blocks of a blastp block,
//!   `blastp/x_query_blocks.rs`): the one-latch appends pass runs them after its chunks (the indices
//!   past the chunks), each once, also when the build then ends `TooLarge` or refused. They read
//!   nothing the build writes and the build reads nothing they write.
//! - the two passes apart (item SC2d): the appends pass (`x_flat_appends_on`) hands the stripes pass
//!   (`x_flat_stripes_on`) owned values (`XFlatAppends`: the header, the count ranges, the chunks,
//!   the shift, the chunks' appends, each stripe's hits), which the stripes pass reads as it reads
//!   them when it follows at once (`x_flat_parallel_on`); a build ahead of its block runs them apart,
//!   possibly on other threads and later (`x_flat_ahead_appends_in_pool`,
//!   `x_flat_ahead_stripes_in_pool`), with the same values. The words of the appends pass (`XWords`,
//!   a function of the header, the threshold and the matrix) are made again for the capacity test.
//! - builds ahead (LOSAT_X_BPLKPAHEAD, `XFlatArgs::x_ahead`: the next query block's table, built by
//!   pool tasks during the block before it, its two passes apart (item SC2d),
//!   `blastp/blast_engine/x_lookup_ahead.rs`): the same build,
//!   whose `starts`, `hits` and `pv` all come from the spare that the block search's thread made ready
//!   (`XFlatSpare::x_make_room`), so that no vector of the table is allocated by a pool thread (a
//!   pool thread's allocations stay in its own glibc arena). `starts` and `hits` are sized and written
//!   as above (the room is there, so `x_reuse` takes the buffers); `pv` is cleared and filled with
//!   zeros, the values of a fresh zeroed vector. Where one of the three lacks the room the build stops
//!   after the appends pass and makes no table (`None`, the spare keeps the vectors); the block then
//!   builds its table itself.
//!
//! With `LOSAT_X_BPSTRIPEFILL` (a strict-set switch: on unless `=0`) the one-latch stripes pass
//! makes each stripe's part another way (`x_stripe_fill`, `XParallelCuts::stripe_fill`); the table
//! is the same: the stripe reads its parts of the chunks once, in chunk order, into a copy while it
//! counts the cells (the copy holds the cells' appends in the serial order, chunk by chunk, each
//! packed as in its chunk), then lays out and fills from the copy as `x_stripe_build` does, the
//! copy's part of chunk c taking chunk c's first offset.
//!
//! Waves (item SC9i of line NP): the backbone's cells are cut by their last letter (the index's
//! most significant digit, weight `A^(word_length - 1)`) into waves of `x_wave_letters` letters
//! each (`x_wave`: about `X_WAVE_CELLS` cells; one wave for a build ahead of its block, whose
//! passes run apart), and the build makes the appends to one wave's cells, sorts them by stripe
//! (the wave's stripes, 2^shift cells each from its first cell) and runs the wave's stripes before
//! it makes the next wave's appends, so that the packed appends held at once are one wave's instead
//! of the table's (4 bytes per hit; for `-word_size 7` they were twice the table's hits, 182 MB of a
//! 50-cDNA query batch at its end). The table is the same:
//! - a wave's appends are the chunks' appends to its cells, in order (`XWords::for_each_word_in`
//!   gives the cells of `for_each_word` whose last letter is in a range, each once, less the wave's
//!   first cell; the offsets of a chunk in order). The first of more than one wave keeps what it
//!   reads of each offset's word before the last letter's loop (`XWords::word_start`: the bound of
//!   the positions before the last, the word's own cell, the last residue; 16 bytes per offset),
//!   and the later waves read it (`XStarts`) instead of the query's letters, the same values. They
//!   differ from the one-wave build's only in their order within one query
//!   offset, and an offset appends each cell at most once (`for_each_word_in`), so each cell's
//!   appends, read by stripe as above, are its offsets in the serial order: the count, layout, fill
//!   and scan order of each stripe are those above.
//! - the waves' cells are the backbone's in order and their stripes partition them, each stripe in
//!   one wave; a wave's hits follow the waves' before it (its first stripe's first hit is their
//!   hits), so the stripes' slices of the hits are in the serial layout's order as above, and each
//!   cell's start is the serial start. A packed append's cell is its index less its stripe's first
//!   cell (below 2^shift).
//! - `TooLarge` when the appends of the waves so far pass the bound (`appended` across the waves):
//!   the build's appends pass it exactly when the serial build's do. The overflow cells of each count
//!   range are summed over every wave's stripes (each cell in one wave, its appends counted in their
//!   order), counted in every wave where there is more than one (a range's appends are known only at
//!   the end), and tested once at the end as above. With one wave they are summed only when a
//!   table's appends pass the capacity (`need_cells`; an append takes at most one cell).
//! - with more than one wave the passes after the first (each wave's stripes and the next wave's
//!   appends) run on one gang of the pool's threads per build (`XGang`): helpers spawned once take
//!   each pass's items and wait between passes, so that a pass does not wake the pool's idle threads
//!   again (scheduling only: each item is made once, by the same call). A build of one wave has one
//!   pass after its first (its stripes), a one-latch pass of its own (`XPasses::latched`).
//! - the PV array is set by the stripes only with one wave (its stripes are whole PV words from cell
//!   0); otherwise it is made from the starts as above.
//!
//! The stripes write the hits in place (`MaybeUninit`, each slot once; the length is set after each
//! wave), so no page of the hits is touched before its stripe writes it; with more than one wave the
//! vector is a block glibc maps (`x_mapped_entries`) and grows by each wave's hits (`mremap`); the
//! segments of packed appends go from one wave to the next (`XSegmentPool`).
//!
//! On with `LOSAT_X_BPCLUT` (a strict-set switch: on unless `=0`)
//! where the search's pool has more than one thread; at one thread the serial builder runs.

use super::*;
use crate::utils::threading::SearchPool;
use rayon::prelude::*;
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::{Condvar, Mutex, OnceLock, PoisonError};

// No NCBI counterpart: how the build is divided among the pool's threads. Scheduling only.
/// The fewest stripes per thread: many small stripes, so that the tasks of the stripes pass end
/// together and the stripe fill's copy (`x_stripe_fill`: a task's copy of one stripe's appends)
/// stays small; the copies of all tasks hold about `1 / X_PARTS_PER_THREAD` of the table's hits
/// (item SC2d-2).
const X_PARTS_PER_THREAD: usize = 64;
/// Guided chunk sizes (`XChunkSizes::Guided`, item SC2d-2): a chunk takes the offsets left over
/// `X_GUIDED_SHARE` times the threads, so that the first chunks are large (few chunks and parts: each
/// chunk's sort and each (chunk, stripe) part costs a fixed amount) and the last ones small (the tasks
/// of the appends pass end together).
const X_GUIDED_SHARE: usize = 2;
// No NCBI counterpart: where the chunks' packed appends are stored. Memory placement only.
/// The fewest entries of a segment of packed appends (`XSegmentOut`): 40 MiB, above glibc's largest
/// mmap threshold (32 MiB), so that each segment is a mapping of its own, returned whole when it is
/// dropped, whichever thread made it (a pool thread's smaller blocks stay in its arena). The pages
/// past the entries written are not touched.
#[cfg(not(target_arch = "wasm32"))]
const X_APPENDS_SEGMENT: usize = 10 << 20;
/// The fewest entries of a segment on wasm32: a piece's. There is no mapping to return there (the
/// linear memory only grows, and the allocator clears what it hands out), so a segment holds about
/// what its task writes.
#[cfg(target_arch = "wasm32")]
const X_APPENDS_SEGMENT: usize = X_PIECE_APPENDS;
/// The fewest query offsets of a chunk (about 65 appends each with blastp-fast's threshold; about
/// 12,800 with BLOSUM45 at word size 7).
const X_MIN_CHUNK_OFFSETS: usize = 256;
/// The most query offsets of a guided chunk (item SC2d-2): a task's scratch holds at most the appends
/// of this many offsets or a piece (`X_PIECE_APPENDS`), whichever is fewer (about 66,000 appends,
/// 260 KiB, with blastp-fast's threshold).
const X_MAX_GUIDED_OFFSETS: usize = 4 * X_MIN_CHUNK_OFFSETS;
// No NCBI counterpart: how a task's scratch is bounded (item SC2d-2). Memory only.
/// The appends of a piece (`x_chunk_appends`): a chunk's appends are sorted into the task's segment
/// in pieces of this many (the last one fewer), so that a task's scratch holds at most a piece and
/// the appends of the offset that passes it (512 KiB of targets, in a vector of 1 MiB once an offset
/// passes it) whatever a chunk's offsets add (about 380 appends per offset with BLOSUM45 and
/// blastp-fast's threshold, 12,800 at word size 7).
/// A task's scratch stays in its pool thread's malloc arena after the pass, so that the threads keep
/// a piece each between builds: 2^20 appends (4 MiB) per thread raised the peak RSS of BLOSUM45
/// blastp-fast at 8 threads by about 10 % (SC2d-2).
const X_PIECE_APPENDS: usize = 1 << 17;
/// The most (chunk, stripe) parts the stripe count allows: a stripe reads its part of every piece of
/// every chunk, so the (piece, stripe) parts are this many times the pieces per chunk (one for
/// blastp-fast's chunks; tens with BLOSUM45 at word size 7).
const X_MAX_PARTS: usize = 1 << 20;
// No NCBI counterpart: the 16 bits of a packed append that hold the offset less the chunk's first.
/// A chunk's offsets lie in `base..base + X_OFFSET_SPAN`.
const X_OFFSET_SPAN: i64 = 1 << 16;
// No NCBI counterpart: the 16 bits of a packed append that hold the cell's index in its stripe.
/// The widest stripe: 2^16 cells.
const X_MAX_STRIPE_SHIFT: usize = 16;

// No NCBI counterpart: the order of the build's work (item SC9i of line NP, module header).
/// The cells of a wave, about: a backbone's waves are the most of `ceil(backbone_size /
/// X_WAVE_CELLS)` with as many last letters each (`x_wave_letters`): one wave for word size 5's
/// 15^5 + 1 cells, 8 waves of 2 letters (1,518,750 cells) for word size 6's 15^6 + 1, 10 of one
/// letter (10^6 cells) for word size 7's 10^7 + 1.
const X_WAVE_CELLS: usize = 1 << 20;

// No NCBI counterpart: the overflow cells of the counts so far (`x_overflow_cells_of`), read where a
// stripe's appends pass from one count range to the next (module header).
/// The overflow cells that the appends counted in `counts` take.
fn x_overflow_cells_in(counts: &[u32]) -> usize {
    counts.iter().map(|&n| x_overflow_cells_of(n)).sum()
}

/// A contiguous run of the query offsets of one count range.
struct XChunk {
    /// The count range (index into `x_count_ranges`).
    range: usize,
    /// Segments whose offsets are the chunk's, in order.
    pieces: Vec<(i32, i32)>,
    /// The chunk's first offset; its offsets lie in `base..base + X_OFFSET_SPAN`.
    base: i32,
}

/// A piece of a chunk's appends to the cells of a wave (`XWave`; the appends of a run of its offsets,
/// in order) sorted by stripe, stably: its packed appends `packed` are
/// `segments[segment][at..at + len]` of the pass's `XAppends`; the wave's stripe j's (`local` j,
/// first cell `first = wave.cells.start + (j << shift)`) are
/// `packed[stripe_starts[j]..stripe_starts[j + 1]]`, each `(index - first) << 16 | (query_offset -
/// base)` with the chunk's `base`.
struct XChunkAppends {
    /// The chunk (index into the pass's chunks).
    chunk: usize,
    segment: usize,
    at: usize,
    len: usize,
    stripe_starts: Vec<u32>,
}

// No NCBI counterpart: the order of the build's work (item SC9i of line NP, module header).
/// A wave of the build: the cells of a run of last letters (the index's most significant digit),
/// whose appends the build makes, sorts by stripe and lays out and fills before it makes those of
/// the next wave.
struct XWave {
    /// The last letters of the wave's cells.
    letters: std::ops::Range<u32>,
    /// The wave's cells: those of its last letters (in the last wave, the backbone's last cell too,
    /// which no word takes).
    cells: std::ops::Range<u32>,
    /// The wave's stripes: its cells from its first, 2^shift each (the last one fewer).
    stripes: usize,
}

/// The appends pass's result: the pieces of the chunks' appends (in chunk order, each chunk's in the
/// order of its offsets), their packed appends in the segments that the pass's tasks filled.
struct XAppends {
    segments: Vec<Vec<u32>>,
    pieces: Vec<XChunkAppends>,
}

impl XAppends {
    /// A piece's packed appends.
    fn packed(&self, piece: &XChunkAppends) -> &[u32] {
        if piece.len == 0 {
            return &[];
        }
        &self.segments[piece.segment][piece.at..piece.at + piece.len]
    }

    /// The pieces and their packed appends, in order.
    fn iter(&self) -> impl Iterator<Item = (&XChunkAppends, &[u32])> {
        self.pieces.iter().map(|piece| (piece, self.packed(piece)))
    }
}

// No NCBI counterpart: memory placement of the packed appends (`X_APPENDS_SEGMENT`).
/// A task's segment of packed appends: chunks take ranges of it in turn; a full segment (and the
/// last one, when the task ends) goes to its slot of `registry`, whose index the chunks keep.
struct XSegmentOut<'r> {
    registry: &'r Mutex<Vec<Vec<u32>>>,
    /// The fewest entries of a segment (`XParallelCuts::appends_segment`).
    size: usize,
    /// Item SC2d-2: the block search's kept segments, taken before a new one is made.
    pool: Option<&'r XSegmentPool>,
    segment: Vec<u32>,
    /// The segment's slot in `registry` (`usize::MAX`: none yet).
    index: usize,
    used: usize,
}

impl<'r> XSegmentOut<'r> {
    fn new(
        registry: &'r Mutex<Vec<Vec<u32>>>,
        size: usize,
        pool: Option<&'r XSegmentPool>,
    ) -> Self {
        Self {
            registry,
            size,
            pool,
            segment: Vec::new(),
            index: usize::MAX,
            used: 0,
        }
    }

    /// `n` entries for one piece (each range is taken once; the piece writes every entry of it
    /// before any is read): the segment's slot, the range's start, the range.
    fn take(&mut self, n: usize) -> (usize, usize, &mut [u32]) {
        if self.segment.len() - self.used < n {
            self.finish();
            let mut registry = self.registry.lock().unwrap_or_else(PoisonError::into_inner);
            self.index = registry.len();
            registry.push(Vec::new());
            drop(registry);
            // A kept segment's values are not read: the range is written before it is read.
            self.segment = self
                .pool
                .and_then(|pool| pool.take(n))
                .unwrap_or_else(|| vec![0u32; n.max(self.size)]);
            self.used = 0;
        }
        let at = self.used;
        self.used += n;
        (self.index, at, &mut self.segment[at..at + n])
    }

    /// Moves the segment to its slot, cut to the entries the chunks took (item SC2d-2: a kept
    /// segment holds no more than they wrote, `XSegmentPool::keep`).
    fn finish(&mut self) {
        if self.index != usize::MAX {
            self.segment.truncate(self.used);
            self.registry.lock().unwrap_or_else(PoisonError::into_inner)[self.index] =
                std::mem::take(&mut self.segment);
        }
    }
}

impl Drop for XSegmentOut<'_> {
    fn drop(&mut self) {
        self.finish();
    }
}

/// A task's scratch in the appends pass: a chunk's targets and appends per offset, and the task's
/// segment of packed appends.
struct XAppendsScratch<'r> {
    targets: Vec<u32>,
    per_offset: Vec<(i32, u32)>,
    out: XSegmentOut<'r>,
}

impl<'r> XAppendsScratch<'r> {
    fn new(
        registry: &'r Mutex<Vec<Vec<u32>>>,
        segment: usize,
        pool: Option<&'r XSegmentPool>,
    ) -> Self {
        Self {
            targets: Vec::new(),
            per_offset: Vec::new(),
            out: XSegmentOut::new(registry, segment, pool),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// NCBI builds the table on one thread before the search threads start
// (c++/src/algo/blast/api/prelim_stage.cpp:145-188 starts them for the search); here the pool that
// will search builds it, the same table (module header).
/// EXPERIMENT (LOSAT_X_BPCLUT): `x_try_build_flat_compressed_lookup` on `pool`'s threads; `pool` has
/// more than one thread. With LOSAT_X_BPLKPFLAT (`x_flat_args.mode`) the one-latch passes on the
/// spare's vectors, with the extra tasks (module header).
pub(super) fn x_try_build_flat_in_pool(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<XFlatBuild> {
    let threads = pool.threads();
    let build = |passes, spare, segments, x_ahead, stripe_fill| {
        pool.install(|| {
            x_try_build_flat_parallel_with(
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                XParallelCuts {
                    overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
                    max_appends: None,
                    threads,
                    chunk_sizes: &[],
                    stripe_shift: None,
                    wave_cells: None,
                    appends_segment: X_APPENDS_SEGMENT,
                    piece_appends: X_PIECE_APPENDS,
                    segments,
                    passes,
                    spare,
                    lent: None,
                    x_ahead,
                    stripe_fill,
                },
            )
        })
    };
    if x_flat_args.mode == XLkpFlat::Off {
        // LOSAT_X_BPLKPAHEAD: a build ahead of its block takes only the switch's passes.
        if x_flat_args.x_ahead {
            return None;
        }
        return build(XPasses::Rayon, None, None, false, false);
    }
    let built = build(
        XPasses::OneLatch(x_flat_args.extra),
        x_flat_args.spare,
        x_flat_args.segments,
        x_flat_args.x_ahead,
        x_bpstripefill_on(),
    );
    if matches!(built, Some(XFlatBuild::Table(_))) {
        X_LKPFLAT_TABLES.fetch_add(1, Ordering::Relaxed);
    }
    built
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_try_build_flat_in_pool`)
/// `x_try_build_flat_in_pool` with the words of the search matrix `matrix`, any but BLOSUM62 (items
/// SC2b, SC2d-2): the same passes, spare vectors and build ahead (`x_flat_args`; another matrix's
/// build has no extra tasks, its contexts' Karlin blocks are computed before it); `pool` has more
/// than one thread.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_try_build_flat_in_pool_of(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<XFlatBuild> {
    let threads = pool.threads();
    let build = |passes, spare, segments, x_ahead, stripe_fill| {
        pool.install(|| {
            x_try_build_flat_parallel_of(
                matrix,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                XParallelCuts {
                    overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
                    max_appends: None,
                    threads,
                    chunk_sizes: &[],
                    stripe_shift: None,
                    wave_cells: None,
                    appends_segment: X_APPENDS_SEGMENT,
                    piece_appends: X_PIECE_APPENDS,
                    segments,
                    passes,
                    spare,
                    lent: None,
                    x_ahead,
                    stripe_fill,
                },
            )
        })
    };
    if x_flat_args.mode == XLkpFlat::Off {
        // LOSAT_X_BPLKPAHEAD: a build ahead of its block takes only the switch's passes.
        if x_flat_args.x_ahead {
            return None;
        }
        return build(XPasses::Rayon, None, None, false, false);
    }
    let built = build(
        XPasses::OneLatch(x_flat_args.extra),
        x_flat_args.spare,
        x_flat_args.segments,
        x_flat_args.x_ahead,
        x_bpstripefill_on(),
    );
    if matches!(built, Some(XFlatBuild::Table(_))) {
        X_LKPFLAT_TABLES.fetch_add(1, Ordering::Relaxed);
    }
    built
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_try_build_flat_in_pool`)
// NCBI builds a query chunk's table on one thread before the search threads start
// (c++/src/algo/blast/api/prelim_stage.cpp:145-188 starts them for the search); here a BLASTX unit
// builds it on the threads the unit cap leaves it, the same table (module header).
/// Item SC9b (2) of line NP: BLASTX's flat table of a query batch with the words of `matrix`, where
/// the bucket fill declines (`x_clut_direct.rs` `x_flat_in_parts`), built by `parts` tasks of the pool
/// the build runs on (`LOSAT_X_BXBUILDPAR`, `blastx/x_bx_mem.rs`): the one-latch passes
/// (`x_one_latch_pass`, the stripes by `x_stripe_fill` as LOSAT_X_BPSTRIPEFILL chooses), without
/// spare vectors, kept segments or extra tasks, the table's vectors `lent` by BLASTX's pool
/// (`LOSAT_X_BXCLUTPOOL`; else fresh ones), its hits of offset type `H` written by the stripes.
/// `parts` is above 1.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_try_build_flat_in_parts<H: XQueryOffset>(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    parts: usize,
    lent: Option<XFlatTableVectors<H>>,
) -> Option<XFlatBuild<H>> {
    let cuts = XParallelCuts {
        overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
        max_appends: None,
        threads: parts,
        chunk_sizes: &[],
        stripe_shift: None,
        wave_cells: None,
        appends_segment: X_APPENDS_SEGMENT,
        piece_appends: X_PIECE_APPENDS,
        segments: None,
        passes: XPasses::OneLatch(None),
        spare: None,
        lent,
        x_ahead: false,
        stripe_fill: x_bpstripefill_on(),
    };
    if matrix == ScoringMatrix::Blosum62 {
        x_try_build_flat_parallel_with(word_size, scaled_threshold, query, locations, tables, cuts)
    } else {
        x_try_build_flat_parallel_of(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            cuts,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_try_build_flat_in_pool`)
/// Item SC2d (LOSAT_X_BPLKPAHEAD): the first part of `x_try_build_flat_in_pool` (of
/// `x_try_build_flat_in_pool_of` for another search matrix `matrix`, item SC2d-2) for a build ahead
/// of its block (`x_flat_args.x_ahead`, the switch's mode on): the one-latch appends pass with the
/// extra tasks, on `pool`'s threads; the second part (`x_flat_ahead_stripes_in_pool`) runs once the
/// block search's thread has made the table's vectors ready (module header). `None` where the build
/// ends at the appends pass (`TooLarge`, or a backbone past the stripes, whose serial builder a build
/// ahead does not run: it would allocate the table on a pool thread) or the mode is `Off`: the block
/// then builds its table itself. `pool` has more than one thread.
#[allow(clippy::too_many_arguments)]
pub(super) fn x_flat_ahead_appends_in_pool(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<XFlatAppends> {
    if x_flat_args.mode == XLkpFlat::Off || !x_flat_args.x_ahead {
        return None;
    }
    let cuts = XParallelCuts::<i32> {
        overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
        max_appends: None,
        threads: pool.threads(),
        chunk_sizes: &[],
        stripe_shift: None,
        wave_cells: None,
        appends_segment: X_APPENDS_SEGMENT,
        piece_appends: X_PIECE_APPENDS,
        segments: x_flat_args.segments,
        passes: XPasses::OneLatch(x_flat_args.extra),
        // The spare's vectors are the second part's.
        spare: None,
        lent: None,
        x_ahead: true,
        stripe_fill: x_bpstripefill_on(),
    };
    let stage = pool.install(|| {
        if matrix == ScoringMatrix::Blosum62 {
            x_flat_appends_on(
                Blosum62Compressed,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                cuts,
                |_| None,
            )
        } else {
            x_flat_appends_on(
                matrix,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                cuts,
                |_| None,
            )
        }
    });
    match stage {
        XFlatStage::Appends(appends) => Some(appends),
        XFlatStage::Built(_) => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_try_build_flat_in_pool`)
/// Item SC2d (LOSAT_X_BPLKPAHEAD): the second part of a build ahead of its block, after
/// `x_flat_ahead_appends_in_pool` gave `appends` (of `matrix`, `query` and `locations`): the stripes
/// pass on `pool`'s threads, on the vectors that `spare` holds, all of them (without their room no
/// table is made and the spare keeps them, `XFlatArgs::x_ahead`).
pub(super) fn x_flat_ahead_stripes_in_pool(
    matrix: ScoringMatrix,
    appends: XFlatAppends,
    query: &[u8],
    locations: &[(i32, i32)],
    pool: &SearchPool<'_>,
    spare: &XFlatSpare,
    segments: Option<&XSegmentPool>,
) -> Option<XFlatBuild> {
    let built = pool.install(|| {
        if matrix == ScoringMatrix::Blosum62 {
            x_flat_stripes_on(
                Blosum62Compressed,
                appends,
                query,
                locations,
                Some(spare),
                None,
                true,
                segments,
            )
        } else {
            x_flat_stripes_on(
                matrix,
                appends,
                query,
                locations,
                Some(spare),
                None,
                true,
                segments,
            )
        }
    });
    if matches!(built, Some(XFlatBuild::Table(_))) {
        X_LKPFLAT_TABLES.fetch_add(1, Ordering::Relaxed);
    }
    built
}

// No NCBI counterpart: reads the LOSAT_X_BPSTRIPEFILL switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSTRIPEFILL): a strict-set switch, on unless it is `0`: the one-latch
/// stripes pass reads each stripe's parts once and fills from the copy (`x_stripe_fill`, module
/// header).
fn x_bpstripefill_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPSTRIPEFILL"))
}

/// EXPERIMENT (LOSAT_X_BPLKPFLAT): how the build's two passes reach the pool's threads.
#[derive(Clone, Copy)]
enum XPasses<'a> {
    /// `par_iter` (rayon's tree of joins per pass).
    Rayon,
    /// One-latch passes (`x_one_latch_pass`), the extra tasks (if any) in the appends pass after the
    /// chunks.
    OneLatch(Option<&'a dyn XExtraTasks>),
    /// Item SC9i of line NP: one-latch passes on a gang's helpers (`XGang::pass`), without extra
    /// tasks (the later waves' appends).
    Gang(&'a XGang),
}

impl XPasses<'_> {
    // No NCBI counterpart: scheduling only (item SC9i of line NP).
    /// A pass of `items` without extra tasks, the results of `item(scratch, i)` in order: on the gang's
    /// helpers (`XGang::pass`, `last`: the gang's last pass), else a one-latch pass of its own
    /// (`x_one_latch_pass`; the Rayon passes' stripes are `par_iter`'s, their callers').
    fn latched<S, R: Send + Sync>(
        self,
        threads: usize,
        items: usize,
        last: bool,
        init: impl Fn() -> S + Sync,
        item: impl Fn(&mut S, usize) -> R + Sync,
    ) -> Vec<R> {
        match self {
            XPasses::Gang(gang) => gang.pass(items, last, init, item),
            XPasses::OneLatch(_) | XPasses::Rayon => {
                x_one_latch_pass(threads, items, 0, init, item, |_| {})
            }
        }
    }
}

impl<H> XFlatSpare<H> {
    // No NCBI counterpart: buffer reuse only (`XFlatSpare`).
    /// The vectors kept (empty ones when none are); the spare then keeps none.
    fn take(&self) -> (Vec<u32>, Vec<H>) {
        std::mem::take(&mut *self.0.lock().unwrap_or_else(PoisonError::into_inner))
    }

    // No NCBI counterpart: which thread allocates the table's vectors (LOSAT_X_BPLKPAHEAD); it
    // changes no value NCBI computes.
    /// LOSAT_X_BPLKPAHEAD: whether the vectors kept have the room for `starts`, `hits` and `pv`
    /// entries (`XFlatSpare::x_make_room`); if so, the PV array's vector is taken from the spare (the
    /// others are taken with `take`).
    fn x_room_for(&self, starts: usize, hits: usize, pv: usize) -> Option<Vec<u32>> {
        let kept = self.0.lock().unwrap_or_else(PoisonError::into_inner);
        let mut kept_pv = self.1.lock().unwrap_or_else(PoisonError::into_inner);
        (kept.0.capacity() >= starts && kept.1.capacity() >= hits && kept_pv.capacity() >= pv)
            .then(|| std::mem::take(&mut *kept_pv))
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1217-1219 (`finalize`)
// Buffer reuse (LOSAT_X_BPLKPAHEAD): the PV array is zeroed as `calloc` zeroes NCBI's.
/// `len` zeroed entries in `pv`'s buffer, which has the room (`XFlatSpare::x_room_for`).
fn x_zeroed_in(mut pv: Vec<u32>, len: usize) -> Vec<u32> {
    debug_assert!(pv.capacity() >= len);
    pv.clear();
    pv.resize(len, 0);
    pv
}

// No NCBI counterpart: buffer reuse (LOSAT_X_BPLKPFLAT; module header: no entry is read before the
// stripes write it). It changes no value NCBI computes.
/// `len` entries in `spare`'s buffer when it has room for them (its values kept, new entries 0),
/// else in a fresh zeroed vector; and whether `spare`'s buffer was used.
fn x_reuse<T: Copy + Default>(spare: Vec<T>, len: usize) -> (Vec<T>, bool) {
    if spare.capacity() == 0 || spare.capacity() < len {
        drop(spare);
        return (vec![T::default(); len], false);
    }
    let mut reused = spare;
    reused.resize(len, T::default());
    (reused, true)
}

// No NCBI counterpart: how a pass's items reach the pool's threads (LOSAT_X_BPLKPFLAT). Scheduling
// only: each item is made once, by the same call, and its result is stored at its index.
/// The results of `item(scratch, i)` for i in 0..items, in order, made by up to `threads` tasks (at
/// most one per item or extra task) in one `rayon::scope` of the current pool: the calling thread and
/// the tasks it spawns each keep one `scratch` (`init`) and take the next index from one counter
/// until none is left; the indices past the items run `extra(i - items)`, `extras` of them.
///
/// The calling thread then waits for the items the spawned tasks took (item SC9i of line NP), not in
/// the scope's wait, which runs other jobs of the pool meanwhile: a job taken there (a subject of
/// another search) held the pass's end, by a few milliseconds each time, until it was done. It spins
/// a little (`X_PASS_END_SPINS`; an item left is short) and then sleeps until the task that ends the
/// last item wakes it, so that a pool of more threads than cores does not spin. The wait cannot
/// deadlock: an item runs no job of the pool (no nested parallel work) and depends on nothing but its
/// inputs, so every item taken ends; a task spawned and not yet started takes no item (the counter
/// is spent) and the scope's wait runs it at once.
fn x_one_latch_pass<S, R: Send + Sync>(
    threads: usize,
    items: usize,
    extras: usize,
    init: impl Fn() -> S + Sync,
    item: impl Fn(&mut S, usize) -> R + Sync,
    extra: impl Fn(usize) + Sync,
) -> Vec<R> {
    let total = items + extras;
    let results: Vec<OnceLock<R>> = (0..items).map(|_| OnceLock::new()).collect();
    if total > 0 {
        let next = AtomicUsize::new(0);
        let ended = XPassEnd {
            items: AtomicUsize::new(0),
            total,
            sleeping: Mutex::new(()),
            woken: Condvar::new(),
        };
        let (results, next, ended, init, item, extra) =
            (&results, &next, &ended, &init, &item, &extra);
        let take_items = move || {
            if next.load(Ordering::Relaxed) >= total {
                return;
            }
            let mut scratch = init();
            loop {
                let index = next.fetch_add(1, Ordering::Relaxed);
                if index >= total {
                    break;
                }
                // Counts the item as ended also when it unwinds, so that the calling thread's wait
                // ends and the scope passes the panic on.
                let _ended = XItemEnded(ended);
                if index < items {
                    let stored = results[index].set(item(&mut scratch, index)).is_ok();
                    assert!(stored, "item {index} taken twice");
                } else {
                    extra(index - items);
                }
            }
        };
        rayon::scope(|scope| {
            for _ in 1..threads.clamp(1, total) {
                scope.spawn(move |_| take_items());
            }
            take_items();
            ended.wait();
        });
    }
    results
        .into_iter()
        .map(|result| result.into_inner().expect("every item taken"))
        .collect()
}

// No NCBI counterpart: scheduling only (item SC9i of line NP, `x_one_latch_pass`).
/// How long the caller of a one-latch pass spins for the items left before it sleeps.
const X_PASS_END_SPINS: u32 = 1 << 12;

// No NCBI counterpart: scheduling only (item SC9i of line NP, `x_one_latch_pass`).
/// The end of a one-latch pass: the items ended, of `total`; the caller sleeps on `woken` (under
/// `sleeping`) until the task that ends the last one wakes it.
struct XPassEnd {
    items: AtomicUsize,
    total: usize,
    sleeping: Mutex<()>,
    woken: Condvar,
}

impl XPassEnd {
    /// Counts one more item as ended; the last one wakes the caller (under the lock, so that the
    /// caller, which tests the count under it before it sleeps, cannot miss it).
    fn end_one(&self) {
        if self.items.fetch_add(1, Ordering::AcqRel) + 1 == self.total {
            let _sleeping = self.sleeping.lock().unwrap_or_else(PoisonError::into_inner);
            self.woken.notify_all();
        }
    }

    /// Returns once every item has ended.
    fn wait(&self) {
        for _ in 0..X_PASS_END_SPINS {
            if self.items.load(Ordering::Acquire) >= self.total {
                return;
            }
            std::hint::spin_loop();
        }
        let mut sleeping = self.sleeping.lock().unwrap_or_else(PoisonError::into_inner);
        while self.items.load(Ordering::Acquire) < self.total {
            sleeping = self
                .woken
                .wait(sleeping)
                .unwrap_or_else(PoisonError::into_inner);
        }
    }
}

/// Counts one item of a one-latch pass as ended when dropped.
struct XItemEnded<'a>(&'a XPassEnd);

impl Drop for XItemEnded<'_> {
    fn drop(&mut self) {
        self.0.end_one();
    }
}

// No NCBI counterpart: scheduling only (item SC9i of line NP, `XGang`).
/// How long a helper of a gang spins for the next pass, and the caller for the helpers to leave a
/// pass or the gang, before they sleep.
const X_GANG_SPINS: u32 = 1 << 10;

// No NCBI counterpart: scheduling only (item SC9i of line NP). Each item is made once, by the same
// call, and its result stored at its index, as in `x_one_latch_pass`.
/// The helpers of a build's passes after its first (`x_flat_stripes_on`: each wave's stripes and the
/// next wave's appends): `threads - 1` tasks of the current pool spawned once, in one `rayon::scope`
/// (`x_with_gang`), that take the items of each pass the calling thread posts (`XGang::pass`) and
/// between passes wait for the next one (a short spin, then `posted`) until the build closes the
/// gang. A pass in a scope of its own (`x_one_latch_pass`) woke the pool's idle threads each time
/// (rayon's spawn), each of which then searched every thread's queue for jobs: with two passes per
/// wave, 8 or 10 waves and 32 threads, 1.4 to 1.6 G more cycles in rayon's steal and 5 % more user
/// time than one wave (ws node 1, 800 cDNAs).
///
/// A helper joins a pass only while it is open (`XGangSlot::take`), and is counted in `inside`
/// until it leaves the pass's item loop; once its own item loop ends (every item is taken), the
/// caller closes the pass to new helpers and waits for those inside, so that the loop, on the
/// caller's stack, outlives every call of it. The waits cannot deadlock: an item runs no job of the
/// pool and depends only on its inputs, and a helper waits for nothing but the caller's next pass or
/// the end; a helper that starts after the end does nothing.
struct XGang {
    threads: usize,
    slot: Mutex<XGangSlot>,
    /// The helpers sleep on it for a pass or the end.
    posted: Condvar,
    /// The caller sleeps on it for the helpers to leave a pass, or to end.
    left: Condvar,
    /// `XGangSlot::generation`, for the helpers' spin.
    generation: AtomicUsize,
    /// The helpers in the item loop of the open or closing pass (raised under `slot`, while the pass
    /// is open).
    inside: AtomicUsize,
    /// The helpers that have started and not yet ended (raised under `slot`, before the end).
    alive: AtomicUsize,
}

/// The pass a gang's helpers may join.
struct XGangSlot {
    /// Raised by each pass posted and by the end.
    generation: usize,
    /// The open pass's item loop; `None` between passes.
    take: Option<XTake>,
    closed: bool,
}

/// A pass's item loop with its lifetime erased: `XGang::pass` keeps the loop alive while a helper
/// can call it.
#[derive(Clone, Copy)]
struct XTake(*const (dyn Fn() + Sync));

// SAFETY: the loop is `Sync` (several threads call it at once), and the pointer is dereferenced only
// by a helper that joined the pass while it was open (`XGang::serve`), which `XGang::pass` waits for.
unsafe impl Send for XTake {}

impl XGang {
    fn new(threads: usize) -> Self {
        Self {
            threads: threads.max(1),
            slot: Mutex::new(XGangSlot {
                generation: 0,
                take: None,
                closed: false,
            }),
            posted: Condvar::new(),
            left: Condvar::new(),
            generation: AtomicUsize::new(0),
            inside: AtomicUsize::new(0),
            alive: AtomicUsize::new(0),
        }
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, XGangSlot> {
        self.slot.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Waits (a short spin, then `left`) until `count` is 0.
    fn wait_for_none(&self, count: &AtomicUsize) {
        for _ in 0..X_GANG_SPINS {
            if count.load(Ordering::Acquire) == 0 {
                return;
            }
            std::hint::spin_loop();
        }
        let mut slot = self.lock();
        while count.load(Ordering::Acquire) != 0 {
            slot = self.left.wait(slot).unwrap_or_else(PoisonError::into_inner);
        }
    }

    /// Lowers `count` by one; the last wakes the caller (under the lock, so that the caller, which
    /// tests the count under it before it sleeps, cannot miss it).
    fn leave(&self, count: &AtomicUsize) {
        if count.fetch_sub(1, Ordering::AcqRel) == 1 {
            let _slot = self.lock();
            self.left.notify_all();
        }
    }

    /// A helper: joins each pass posted until the gang is closed.
    fn serve(&self) {
        {
            // After the end, or after the last pass has closed; a helper that comes while the last
            // pass is open joins it.
            let slot = self.lock();
            if slot.closed && slot.take.is_none() {
                return;
            }
            self.alive.fetch_add(1, Ordering::AcqRel);
        }
        let _ended = XGangLeave(self, &self.alive);
        let mut seen = 0usize;
        loop {
            for _ in 0..X_GANG_SPINS {
                if self.generation.load(Ordering::Acquire) != seen {
                    break;
                }
                std::hint::spin_loop();
            }
            let take = {
                let mut slot = self.lock();
                while slot.generation == seen && !slot.closed {
                    slot = self
                        .posted
                        .wait(slot)
                        .unwrap_or_else(PoisonError::into_inner);
                }
                seen = slot.generation;
                match slot.take {
                    Some(take) => {
                        self.inside.fetch_add(1, Ordering::AcqRel);
                        (take, slot.closed)
                    }
                    // The pass closed before this helper came.
                    None if slot.closed => return,
                    None => continue,
                }
            };
            let (take, last) = take;
            {
                // Leaves the pass also when an item unwinds, so that the caller's wait ends and the
                // scope passes the panic on.
                let _left = XGangLeave(self, &self.inside);
                // SAFETY: the helper joined the pass while it was open (above, under the lock), and
                // `XGang::pass` keeps the loop alive until `inside` is 0.
                unsafe { (*take.0)() };
            }
            if last {
                return;
            }
        }
    }

    /// `x_one_latch_pass` without its extra tasks, on the caller and the gang's helpers: the
    /// results of `item(scratch, i)` for i in 0..items, in order, each task with its own `scratch`
    /// (`init`). After the `last` pass the helpers end at once (they wait for no next pass, and the
    /// end does not wake them).
    fn pass<S, R: Send + Sync>(
        &self,
        items: usize,
        last: bool,
        init: impl Fn() -> S + Sync,
        item: impl Fn(&mut S, usize) -> R + Sync,
    ) -> Vec<R> {
        let results: Vec<OnceLock<R>> = (0..items).map(|_| OnceLock::new()).collect();
        if items > 0 {
            let next = AtomicUsize::new(0);
            let take_items = || {
                if next.load(Ordering::Relaxed) >= items {
                    return;
                }
                let mut scratch = init();
                loop {
                    let index = next.fetch_add(1, Ordering::Relaxed);
                    if index >= items {
                        break;
                    }
                    let stored = results[index].set(item(&mut scratch, index)).is_ok();
                    assert!(stored, "item {index} taken twice");
                }
            };
            let take: *const (dyn Fn() + Sync + '_) = &take_items;
            // SAFETY: only the lifetime is erased; `_open` (also when the caller's items unwind)
            // closes the pass and waits for every helper in it before `take_items` is dropped.
            let take = XTake(unsafe {
                std::mem::transmute::<*const (dyn Fn() + Sync + '_), *const (dyn Fn() + Sync)>(take)
            });
            {
                let mut slot = self.lock();
                slot.generation += 1;
                slot.take = Some(take);
                slot.closed |= last;
                self.generation.store(slot.generation, Ordering::Release);
            }
            let _open = XGangPassOpen(self);
            self.posted.notify_all();
            take_items();
        }
        results
            .into_iter()
            .map(|result| result.into_inner().expect("every item taken"))
            .collect()
    }
}

/// Lowers a count of a gang when dropped (`XGang::leave`).
struct XGangLeave<'g>(&'g XGang, &'g AtomicUsize);

impl Drop for XGangLeave<'_> {
    fn drop(&mut self) {
        self.0.leave(self.1);
    }
}

/// An open pass of a gang: closes it to new helpers and waits for those in it when dropped.
struct XGangPassOpen<'g>(&'g XGang);

impl Drop for XGangPassOpen<'_> {
    fn drop(&mut self) {
        self.0.lock().take = None;
        self.0.wait_for_none(&self.0.inside);
    }
}

/// Closes a gang when dropped: its helpers end, and the caller waits for those that started.
struct XGangClose<'g>(&'g XGang);

impl Drop for XGangClose<'_> {
    fn drop(&mut self) {
        {
            let mut slot = self.0.lock();
            slot.closed = true;
            slot.generation += 1;
            self.0.generation.store(slot.generation, Ordering::Release);
        }
        self.0.posted.notify_all();
        self.0.wait_for_none(&self.0.alive);
    }
}

// No NCBI counterpart: scheduling only (item SC9i of line NP, `XGang`).
/// `body` with a gang of `threads - 1` helpers of the current pool, closed when `body` returns or
/// unwinds.
fn x_with_gang<R: Send>(threads: usize, body: impl FnOnce(&XGang) -> R + Send) -> R {
    let gang = XGang::new(threads);
    let gang = &gang;
    rayon::scope(move |scope| {
        for _ in 1..gang.threads {
            scope.spawn(move |_| gang.serve());
        }
        let _close = XGangClose(gang);
        body(gang)
    })
}

// No NCBI counterpart: scheduling only (line NP, SC9e).
/// `items` on `threads` tasks of the current pool (the one-latch pass), each task with its own
/// scratch from `init`: the passes of BLASTX's bucket fill in parts (`x_clut_direct.rs`
/// `x_bx_buckets_with`). The results in item order.
pub(super) fn x_items_in_parts<S, R: Send + Sync>(
    threads: usize,
    items: usize,
    init: impl Fn() -> S + Sync,
    item: impl Fn(&mut S, usize) -> R + Sync,
) -> Vec<R> {
    x_one_latch_pass(threads, items, 0, init, item, |_| {})
}

/// The parameters of `x_try_build_flat_parallel_with` that the module header's argument leaves free
/// (and the capacity and the bound of the serial builder, for the tests).
struct XParallelCuts<'a, H = i32> {
    overflow_cell_capacity: usize,
    /// `None`: `x_max_appends`, the serial builder's bound.
    max_appends: Option<usize>,
    threads: usize,
    /// The chunks' sizes in query offsets, in turn; empty: guided sizes (`XChunkSizes::Guided`).
    chunk_sizes: &'a [usize],
    /// The stripes' width 2^shift; `None`: the widest that makes `X_PARTS_PER_THREAD` stripes per
    /// thread, or fewer so that the chunks' parts are at most `X_MAX_PARTS`.
    stripe_shift: Option<usize>,
    /// Item SC9i of line NP: the cells of a wave (`x_wave`); `None`: `X_WAVE_CELLS` (smaller in the
    /// tests).
    wave_cells: Option<usize>,
    /// The fewest entries of a segment of packed appends (`XSegmentOut`): `X_APPENDS_SEGMENT`
    /// (smaller in the tests).
    appends_segment: usize,
    /// The appends of a piece (`x_chunk_appends`): `X_PIECE_APPENDS` (smaller in the tests).
    piece_appends: usize,
    /// Item SC2d-2: the block search's kept segments of packed appends (`XSegmentPool`): the appends
    /// pass takes them before it makes new ones, the stripes pass keeps the build's.
    segments: Option<&'a XSegmentPool>,
    /// LOSAT_X_BPLKPFLAT: how the passes are scheduled, and the extra tasks.
    passes: XPasses<'a>,
    /// LOSAT_X_BPLKPFLAT: the vectors of an earlier table for `starts` and `hits` (else fresh ones).
    spare: Option<&'a XFlatSpare<H>>,
    /// LOSAT_X_BXCLUTPOOL: the vectors BLASTX's pool lent for `starts`, `hits` and `pv` (else fresh
    /// ones; not with `spare`).
    lent: Option<XFlatTableVectors<H>>,
    /// LOSAT_X_BPLKPAHEAD: `starts`, `hits` and `pv` all come from `spare`, and without their room
    /// the build makes no table (`None`); nothing of the table is allocated here
    /// (`XFlatArgs::x_ahead`).
    x_ahead: bool,
    /// LOSAT_X_BPSTRIPEFILL: with one-latch passes, the stripes made by `x_stripe_fill` (else by
    /// `x_stripe_build`).
    stripe_fill: bool,
}

// No NCBI counterpart: the largest PV shift that `finalize` (blast_aalookup.c:1219-1222, below)
// chooses for a backbone of this size; the stripes are whole PV words at any shift up to it.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1213-1214
fn x_largest_pv_shift(backbone_size: usize) -> usize {
    PV_ARRAY_BTS + ilog2(backbone_size / (8 * X_TARGET_PV_BYTES))
}

// No NCBI counterpart: chooses the stripes' width. Scheduling only.
/// The widest stripes (2^shift cells, shift from `smallest` to `X_MAX_STRIPE_SHIFT`) that make at
/// least `stripes` of them, or the narrowest.
fn x_stripe_shift(backbone_size: usize, stripes: usize, smallest: usize) -> usize {
    let mut shift = X_MAX_STRIPE_SHIFT;
    while shift > smallest && backbone_size.div_ceil(1 << shift) < stripes {
        shift -= 1;
    }
    shift
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
// No NCBI counterpart: how the appends pass is divided among the pool's threads. Scheduling only.
/// The chunks' sizes in query offsets (`x_chunks`).
#[derive(Clone, Copy)]
enum XChunkSizes<'a> {
    /// Chunk i takes `sizes[i % sizes.len()]` (the tests' cuts).
    List(&'a [usize]),
    /// Guided (item SC2d-2): a chunk takes the offsets not yet in a chunk (of `offsets`) over
    /// `X_GUIDED_SHARE * threads`, from `X_MIN_CHUNK_OFFSETS` to `X_MAX_GUIDED_OFFSETS`.
    Guided { threads: usize, offsets: usize },
}

impl XChunkSizes<'_> {
    /// The size of chunk `chunk`, opened when `placed` offsets are in the chunks before it.
    fn size(self, chunk: usize, placed: usize) -> usize {
        match self {
            Self::List(sizes) => sizes[chunk % sizes.len()].max(1),
            Self::Guided { threads, offsets } => (offsets.saturating_sub(placed)
                / (X_GUIDED_SHARE * threads.max(1)))
            .clamp(X_MIN_CHUNK_OFFSETS, X_MAX_GUIDED_OFFSETS),
        }
    }
}

/// The count ranges' offsets cut into chunks (module header) of the sizes of `sizes`, fewer at a
/// range's end and where the offsets would pass `X_OFFSET_SPAN`.
fn x_chunks(
    word_size: usize,
    ranges: &[(std::ops::Range<usize>, bool)],
    locations: &[(i32, i32)],
    sizes: XChunkSizes<'_>,
) -> Vec<XChunk> {
    let word = word_size as i64;
    let mut chunks: Vec<XChunk> = Vec::new();
    // The offsets in the chunks so far (the open one's included).
    let mut placed = 0usize;
    for (range_index, (range, _)) in ranges.iter().enumerate() {
        // The chunk being filled (it is chunk `chunks.len()`), its offsets and its size.
        let mut open: Option<(XChunk, usize, usize)> = None;
        for &(left, right) in &locations[range.clone()] {
            let to = i64::from(right) - word + 1;
            let mut from = i64::from(left);
            while from <= to {
                let joins = matches!(&open, Some((chunk, count, size))
                    if *count < *size
                        && from >= i64::from(chunk.base)
                        && from < i64::from(chunk.base) + X_OFFSET_SPAN);
                if !joins {
                    if let Some((chunk, _, _)) = open.take() {
                        chunks.push(chunk);
                    }
                    open = Some((
                        XChunk {
                            range: range_index,
                            pieces: Vec::new(),
                            base: from as i32,
                        },
                        0,
                        sizes.size(chunks.len(), placed),
                    ));
                }
                let (chunk, count, size) = open.as_mut().expect("an open chunk");
                let size = *size;
                // `size - count` is at least 1 (`joins`); past `X_OFFSET_SPAN` it does not matter.
                let take = ((size - *count).min(X_OFFSET_SPAN as usize) as i64)
                    .min(to - from + 1)
                    .min(i64::from(chunk.base) + X_OFFSET_SPAN - from);
                // Offsets from..=from + take - 1: the segment's right end less what follows them.
                chunk
                    .pieces
                    .push((from as i32, (from + take - 1 + word - 1) as i32));
                *count += take as usize;
                placed += take as usize;
                from += take;
            }
        }
        if let Some((chunk, _, _)) = open.take() {
            chunks.push(chunk);
        }
    }
    chunks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (`x_chunks`) and
// :1075-1121 (`s_CompressedAddWordHits`, `XWords::for_each_word`)
/// A chunk's appends (the serial count pass's statements) in pieces of `piece_appends` appends (the
/// last one fewer; none for a chunk of no appends), each sorted by stripe into the task's segment as
/// soon as it is full, so that the scratch holds at most a piece (item SC2d-2); the pieces in order
/// hold the chunk's appends in order. None when the appends of all chunks pass `max_appends` (then
/// `too_large`). The pieces' `chunk` is left 0 (the pass sets it).
///
/// Item SC9i of line NP: the words' starts (`XWords::word_start`) are read from `starts` where the
/// chunk's are kept (the later waves), else made from the query; with `XStarts::Keep` (the first of
/// more than one wave) they are also returned, those of the chunk's offsets in order.
#[allow(clippy::too_many_arguments)]
fn x_chunk_appends(
    words: &XWords<'_>,
    query: &[u8],
    chunk: &XChunk,
    wave: &XWave,
    shift: usize,
    max_appends: usize,
    piece_appends: usize,
    appended: &AtomicUsize,
    too_large: &AtomicBool,
    starts: XStarts<'_>,
    scratch: &mut XAppendsScratch<'_>,
) -> (Vec<XChunkAppends>, Vec<XWordStart>) {
    let XAppendsScratch {
        targets,
        per_offset,
        out,
    } = scratch;
    targets.clear();
    per_offset.clear();
    let mut pieces = Vec::new();
    let mut kept = Vec::new();
    if too_large.load(Ordering::Relaxed) {
        return (pieces, kept);
    }
    // The appends not yet added to `appended`.
    let mut unreported = 0usize;
    let report = |appends: usize| {
        let total = appended.fetch_add(appends, Ordering::Relaxed) + appends;
        if total > max_appends {
            too_large.store(true, Ordering::Relaxed);
        }
        !too_large.load(Ordering::Relaxed)
    };
    // The bookkeeping after an offset's appends (`targets[before..]`); `false` past `max_appends`.
    let mut after = |offset: i32, before: usize, targets: &mut Vec<u32>| {
        let appends = targets.len() - before;
        if appends != 0 {
            per_offset.push((offset, appends as u32));
        }
        unreported += appends;
        if unreported >= X_COUNT_CHUNK && !report(std::mem::take(&mut unreported)) {
            return false;
        }
        if targets.len() >= piece_appends {
            pieces.push(x_sorted_piece(
                targets, per_offset, chunk.base, shift, wave, out,
            ));
            targets.clear();
            per_offset.clear();
        }
        true
    };
    match starts {
        // One wave: every cell, as `XWords::for_each_word` gives them.
        XStarts::Make => {
            x_for_each_offset(words.word_length, query, &chunk.pieces, |offset, word| {
                let before = targets.len();
                words.for_each_word(word, targets);
                after(offset, before, targets)
            });
        }
        XStarts::Keep => {
            x_for_each_offset(words.word_length, query, &chunk.pieces, |offset, word| {
                let Some(start) = words.word_start(word, offset) else {
                    return true;
                };
                kept.push(start);
                let before = targets.len();
                words.for_each_word_from(&start, word, &wave.letters, targets);
                after(offset, before, targets)
            });
        }
        XStarts::Kept(kept_starts) => {
            for start in kept_starts {
                let at = start.offset as usize;
                let before = targets.len();
                words.for_each_word_from(
                    start,
                    &query[at..at + words.word_length],
                    &wave.letters,
                    targets,
                );
                if !after(start.offset, before, targets) {
                    break;
                }
            }
        }
    }
    if !report(unreported) {
        return (Vec::new(), Vec::new());
    }
    if !targets.is_empty() {
        pieces.push(x_sorted_piece(
            targets, per_offset, chunk.base, shift, wave, out,
        ));
    }
    (pieces, kept)
}

// No NCBI counterpart: the order of the build's work (item SC9i of line NP).
/// Where `x_chunk_appends` takes the words' starts (`XWords::word_start`) of a chunk's offsets.
#[derive(Clone, Copy)]
enum XStarts<'a> {
    /// None: a build of one wave, whose appends are `XWords::for_each_word`'s (every cell of an
    /// offset, as it gives them: a wave of fewer letters chooses its last letter first, and
    /// testing each of the alphabet's letters against the bound one by one missed the branch
    /// predictor's guess often: blastp-fast +4.7 % of the wall time at 8 threads with +0.4 % of the
    /// instructions).
    Make,
    /// From the query, and kept for the later waves (the first of more than one wave).
    Keep,
    /// The chunk's kept starts (the later waves).
    Kept(&'a [XWordStart]),
}

// No NCBI counterpart: the layout of a piece of a chunk's appends (module header); it changes no
// value NCBI computes.
/// The appends `targets` of the offsets `per_offset` (each offset with its number of appends, in
/// order; offsets of a chunk whose first offset is `base`), sorted by stripe, stably (counting sort),
/// into a range of the task's segment `out`. A target is its cell less the wave's first cell
/// (`XWords::for_each_word_in`).
fn x_sorted_piece(
    targets: &[u32],
    per_offset: &[(i32, u32)],
    base: i32,
    shift: usize,
    wave: &XWave,
    out: &mut XSegmentOut<'_>,
) -> XChunkAppends {
    // The wave's stripes, from its first cell (the targets lie in its cells, from it).
    let stripes = wave.stripes;
    let mut stripe_starts = vec![0u32; stripes + 1];
    for &index in targets {
        stripe_starts[(index >> shift) as usize + 1] += 1;
    }
    for j in 0..stripes {
        stripe_starts[j + 1] += stripe_starts[j];
    }
    let mut cursors = stripe_starts[..stripes].to_vec();
    let len = targets.len();
    let (segment, at_segment, packed) = out.take(len);
    let mask = (1u32 << shift) - 1;
    let mut at = 0usize;
    for &(query_offset, appends) in per_offset {
        // Below X_OFFSET_SPAN (`x_chunks`).
        let delta = (query_offset - base) as u32;
        let end = at + appends as usize;
        for &index in &targets[at..end] {
            let cursor = &mut cursors[(index >> shift) as usize];
            packed[*cursor as usize] = (index & mask) << 16 | delta;
            *cursor += 1;
        }
        at = end;
    }
    XChunkAppends {
        chunk: 0,
        segment,
        at: at_segment,
        len,
        stripe_starts,
    }
}

/// One stripe's part of the table.
struct XStripe<'t, H> {
    /// The stripe's first cell.
    first: usize,
    /// The stripe's index among its wave's (`XChunkAppends::stripe_starts`).
    local: usize,
    /// The stripe's first hit.
    base: u32,
    starts: &'t mut [u32],
    /// The stripe's slots of the hits, each written once by the fill.
    hits: &'t mut [MaybeUninit<H>],
    /// The stripe's PV words when the PV shift is known before the counts (`x_largest_pv_shift` is
    /// `PV_ARRAY_BTS`), else empty.
    pv: &'t mut [u32],
}

// No NCBI counterpart: memory handling only (the hits written in place, item SC9i of line NP).
/// `slots` as the values written in them.
///
/// # Safety
/// Every slot of `slots` has been written.
unsafe fn x_written<H>(slots: &mut [MaybeUninit<H>]) -> &mut [H] {
    // SAFETY: `MaybeUninit<H>` has the layout of `H`, and the caller wrote every slot.
    unsafe { &mut *(slots as *mut [MaybeUninit<H>] as *mut [H]) }
}

/// What a stripe's counts give the whole table.
struct XStripeCounts {
    nonempty: usize,
    longest: u32,
    /// The overflow cells of the reference's table of all batches that each count range's appends
    /// to the stripe's cells take.
    cells: Vec<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849,1224-1230 (`add_word_hit`,
// module header of `x_clut_direct.rs`), :1197-1235 (`finalize`) and
// c++/src/algo/blast/core/blast_aascan.c:348-396 (the scan order, module header of
// `x_clut_direct.rs`)
/// The count, layout and fill of one stripe (module header).
#[allow(clippy::too_many_arguments)]
fn x_stripe_build<H: XQueryOffset>(
    part: XStripe<'_, H>,
    ranges: usize,
    chunks: &[XChunk],
    appends: &XAppends,
    need_cells: bool,
) -> XStripeCounts {
    let XStripe {
        first,
        local,
        base,
        starts,
        hits,
        pv,
    } = part;
    let parts = || {
        appends.iter().map(move |(piece, packed)| {
            (
                &chunks[piece.chunk],
                &packed
                    [piece.stripe_starts[local] as usize..piece.stripe_starts[local + 1] as usize],
            )
        })
    };

    // Count: each cell's appends in the serial order; with `need_cells`, the overflow cells of the
    // count ranges before the last (module header).
    let mut counts = vec![0u32; starts.len()];
    let mut cells = vec![0usize; ranges];
    let mut range = None;
    let mut cells_before = 0usize;
    for (chunk, part) in parts() {
        if part.is_empty() {
            continue;
        }
        if need_cells {
            x_next_count_range(
                &mut range,
                chunk.range,
                &counts,
                &mut cells,
                &mut cells_before,
            );
        }
        for &packed in part {
            counts[(packed >> 16) as usize] += 1;
        }
    }

    // Layout: the starts (from the stripe's first hit) and the PV bits; each count becomes its
    // cell's cursor, the cell's first slot in the stripe's hits. The overflow cells of the counts.
    let mut nonempty = 0usize;
    let mut longest = 0u32;
    let mut taken = 0usize;
    let mut total = base;
    let mut slot = 0u32;
    for (cell, (start, cursor)) in starts.iter_mut().zip(counts.iter_mut()).enumerate() {
        let count = *cursor;
        *start = total;
        if count > 0 {
            nonempty += 1;
            longest = longest.max(count);
            taken += x_overflow_cells_of(count);
            if !pv.is_empty() {
                // `compressed_pv_set` on the stripe's words (PV shift `PV_ARRAY_BTS`; the stripe
                // starts a word).
                pv[cell >> PV_ARRAY_BTS] |= 1u32 << ((first + cell) & PV_ARRAY_MASK);
            }
        }
        *cursor = slot;
        total += count;
        slot += count;
    }
    if let (true, Some(range)) = (need_cells, range) {
        cells[range] += taken - cells_before;
    }

    // Fill: the serial fill on the stripe's cells, the appends in order, each at its cell's cursor.
    for (chunk, part) in parts() {
        for &packed in part {
            let cursor = &mut counts[(packed >> 16) as usize];
            hits[*cursor as usize].write(H::from_offset(
                (chunk.base + (packed & 0xffff) as i32) as u32,
            ));
            *cursor += 1;
        }
    }
    // SAFETY: the cursors start at the cells' first slots in the stripe (the counts' prefix sums,
    // from 0) and each moves on once per append to its cell, which the count counted: every slot
    // of the stripe's hits, whose number is the stripe's appends, is written once above.
    let hits = unsafe { x_written(hits) };
    // A cursor ends at the next cell's first slot; the cells of more than 5 offsets put o_k at
    // `x_emission_position(k, n)` (the serial builder's statements).
    let mut begin = 0usize;
    for &end in &counts {
        let end = end as usize;
        if end - begin > COMPRESSED_HITS_PER_BACKBONE_CELL + 1 {
            x_to_scan_order(&mut hits[begin..end]);
        }
        begin = end;
    }
    XStripeCounts {
        nonempty,
        longest,
        cells,
    }
}

// No NCBI counterpart: where the overflow cells of the serial count pass are summed per count range
// (module header).
/// A stripe's count passes from the appends of `range` (`None`: none yet) to those of `next`, a
/// count range not before it: the overflow cells that the counts so far take, less those of the
/// ranges before (`cells_before`), are `range`'s.
fn x_next_count_range(
    range: &mut Option<usize>,
    next: usize,
    counts: &[u32],
    cells: &mut [usize],
    cells_before: &mut usize,
) {
    match *range {
        Some(current) if current != next => {
            let so_far = x_overflow_cells_in(counts);
            cells[current] += so_far - *cells_before;
            *cells_before = so_far;
        }
        _ => {}
    }
    *range = Some(next);
}

/// LOSAT_X_BPSTRIPEFILL: a thread's scratch for `x_stripe_fill`: the copy of the stripe's parts,
/// the cells' counts (then cursors), the cells of more than 5 offsets.
type XFillScratch = (Vec<u32>, Vec<u32>, Vec<u32>);

// NCBI reference (598d8ae6): as `x_stripe_build` (c++/src/algo/blast/core/blast_aalookup.c:787-849, :1197-1235)
/// LOSAT_X_BPSTRIPEFILL: `x_stripe_build` with the stripe's parts of the chunks read once, in chunk
/// order, into a copy while the cells are counted; the layout and the fill from the copy (module
/// header).
#[allow(clippy::too_many_arguments)]
fn x_stripe_fill<H: XQueryOffset>(
    part: XStripe<'_, H>,
    ranges: usize,
    chunks: &[XChunk],
    appends: &XAppends,
    need_cells: bool,
    scratch: &mut XFillScratch,
) -> XStripeCounts {
    let XStripe {
        first,
        local,
        base,
        starts,
        hits,
        pv,
    } = part;
    let (copy, counts, overflowed) = scratch;
    counts.clear();
    counts.resize(starts.len(), 0);
    copy.clear();
    copy.reserve(hits.len());
    let part_of = |piece: &XChunkAppends| {
        piece.stripe_starts[local] as usize..piece.stripe_starts[local + 1] as usize
    };

    // Count: each cell's appends in the serial order; with `need_cells`, the overflow cells of the
    // count ranges before the last, as `x_stripe_build`'s.
    let mut cells = vec![0usize; ranges];
    let mut range = None;
    let mut cells_before = 0usize;
    for (piece, packed) in appends.iter() {
        let chunk = &chunks[piece.chunk];
        let part = &packed[part_of(piece)];
        if part.is_empty() {
            continue;
        }
        if need_cells {
            x_next_count_range(
                &mut range,
                chunk.range,
                counts,
                &mut cells,
                &mut cells_before,
            );
        }
        copy.extend(part.iter().map(|&packed| {
            counts[(packed >> 16) as usize] += 1;
            packed
        }));
    }

    // Layout, as `x_stripe_build`'s; the cells of more than 5 offsets (those that take overflow
    // cells) noted.
    let mut nonempty = 0usize;
    let mut longest = 0u32;
    let mut taken = 0usize;
    let mut slot = 0u32;
    overflowed.clear();
    for (cell, (start, cursor)) in starts.iter_mut().zip(counts.iter_mut()).enumerate() {
        let count = *cursor;
        *start = base + slot;
        if count > 0 {
            nonempty += 1;
            longest = longest.max(count);
            if !pv.is_empty() {
                pv[cell >> PV_ARRAY_BTS] |= 1u32 << ((first + cell) & PV_ARRAY_MASK);
            }
            if count > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1 {
                overflowed.push(cell as u32);
                taken += x_overflow_cells_of(count);
            }
        }
        *cursor = slot;
        slot += count;
    }
    if let (true, Some(range)) = (need_cells, range) {
        cells[range] += taken - cells_before;
    }

    // Fill: the serial fill on the stripe's cells, the appends in order, each at its cell's cursor;
    // the copy holds the chunks' parts one after the other.
    let mut at = 0usize;
    for piece in &appends.pieces {
        let chunk = &chunks[piece.chunk];
        let end = at + part_of(piece).len();
        for &packed in &copy[at..end] {
            let cursor = &mut counts[(packed >> 16) as usize];
            hits[*cursor as usize].write(H::from_offset(
                (chunk.base + (packed & 0xffff) as i32) as u32,
            ));
            *cursor += 1;
        }
        at = end;
    }
    // SAFETY: as `x_stripe_build`'s: the copy holds the stripe's appends, each written once at its
    // cell's cursor, so every slot of the stripe's hits is written.
    let hits = unsafe { x_written(hits) };
    // A cursor ends at its cell's end; the cells of more than 5 offsets put o_k at
    // `x_emission_position(k, n)`.
    for &cell in overflowed.iter() {
        let cell = cell as usize;
        let begin = (starts[cell] - base) as usize;
        let end = counts[cell] as usize;
        x_to_scan_order(&mut hits[begin..end]);
    }
    XStripeCounts {
        nonempty,
        longest,
        cells,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (the passes) and
// :1197-1235 (`finalize`, the statements of `x_try_build_flat_with_capacity`'s layout)
/// The parallel build (module header) in the current rayon pool, with the cuts of `cuts`.
fn x_try_build_flat_parallel_with<H: XQueryOffset>(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    cuts: XParallelCuts<'_, H>,
) -> Option<XFlatBuild<H>> {
    x_flat_parallel_on(
        Blosum62Compressed,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        cuts,
        |capacity| {
            x_try_build_flat_with_capacity(
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                capacity,
                X_COUNT_CHUNK,
                0,
            )
            .map(XFlatBuild::x_offsets)
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310 (see
// `x_clut_direct::x_tbn_build_compressed_lookup_of`)
/// `x_try_build_flat_parallel_with` with the words of the search matrix `matrix`, any but BLOSUM62
/// (`build_compressed_score_matrix`, item SC2b): the same passes, and the serial builder of the
/// matrix's words (`x_try_build_flat_of`) where the parallel build takes the serial one.
fn x_try_build_flat_parallel_of<H: XQueryOffset>(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    cuts: XParallelCuts<'_, H>,
) -> Option<XFlatBuild<H>> {
    x_flat_parallel_on(
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        cuts,
        |capacity| {
            x_try_build_flat_of(
                matrix,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                capacity,
                X_COUNT_CHUNK,
            )
            .map(XFlatBuild::x_offsets)
        },
    )
}

// The body of both builders above, so that BLOSUM62's keeps its constant matrix (as
// `x_flat_with_capacity_of`); `serial` is the serial builder of the same words, with the capacity.
// The appends pass (`x_flat_appends_on`) and then the stripes pass (`x_flat_stripes_on`), one after
// the other here; a build ahead of its block may run them apart (item SC2d, module header).
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn x_flat_parallel_on<S: CompressedScores, H: XQueryOffset>(
    scores: S,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    mut cuts: XParallelCuts<'_, H>,
    serial: impl FnOnce(usize) -> Option<XFlatBuild<H>>,
) -> Option<XFlatBuild<H>> {
    let (spare, x_ahead, segments) = (cuts.spare, cuts.x_ahead, cuts.segments);
    let lent = cuts.lent.take();
    match x_flat_appends_on(
        scores,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        cuts,
        serial,
    ) {
        XFlatStage::Built(built) => built,
        XFlatStage::Appends(appends) => x_flat_stripes_on(
            scores, appends, query, locations, spare, lent, x_ahead, segments,
        ),
    }
}

/// What the appends pass of the parallel build gave (`x_flat_appends_on`).
#[expect(
    clippy::large_enum_variant,
    reason = "a value moved once per build, as `XFlatBuild`; boxing a variant adds an allocation"
)]
enum XFlatStage<H: XQueryOffset = i32> {
    /// The build's outcome without a stripes pass: the serial builder's (a backbone past the
    /// stripes), `TooLarge`, or none (no score matrix for the alphabet).
    Built(Option<XFlatBuild<H>>),
    /// The appends, for the stripes pass.
    Appends(XFlatAppends),
}

// No NCBI counterpart: the values the build's first pass hands to its second (module header, item
// SC2d); it changes no value NCBI computes.
/// The parallel build after its appends pass (`x_flat_appends_on`), all owned, so that the stripes
/// pass (`x_flat_stripes_on`) can run later, on another thread: the table's header and cuts, the
/// first wave's appends and each of its stripes' hits. The words (`XWords`) are not kept: the stripes
/// pass needs them only for the later waves' appends and the capacity test of a table that passes
/// the capacity, and makes them again there.
pub(crate) struct XFlatAppends {
    header: XHeader,
    largest_pv_shift: usize,
    ranges: Vec<(std::ops::Range<usize>, bool)>,
    chunks: Vec<XChunk>,
    /// Item SC9i of line NP: with more than one wave, each chunk's words' starts
    /// (`XWords::word_start`), which the later waves read (else none).
    word_starts: Vec<Vec<XWordStart>>,
    appends: XAppends,
    shift: usize,
    /// The first wave's stripes' hits.
    stripe_hits: Vec<usize>,
    /// Item SC9i of line NP: the waves (`x_wave`) and the last letters of each, and the appends of
    /// the first.
    waves: usize,
    wave_letters: usize,
    appended: usize,
    max_appends: usize,
    piece_appends: usize,
    appends_segment: usize,
    scaled_threshold: i32,
    threads: usize,
    overflow_cell_capacity: usize,
    /// The one-latch passes (else rayon's).
    one_latch: bool,
    /// LOSAT_X_BPSTRIPEFILL with the one-latch passes.
    fill: bool,
}

// The first pass of `x_flat_parallel_on` (module header): the chunks' appends, sorted by stripe,
// with the extra tasks of the one-latch passes; `XFlatStage::Built` where the build ends here.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn x_flat_appends_on<S: CompressedScores, H: XQueryOffset>(
    scores: S,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    cuts: XParallelCuts<'_, H>,
    serial: impl FnOnce(usize) -> Option<XFlatBuild<H>>,
) -> XFlatStage<H> {
    let Some(header) = XHeader::new(word_size) else {
        return XFlatStage::Built(None);
    };
    let backbone_size = header.backbone_size;
    let largest_pv_shift = x_largest_pv_shift(backbone_size);
    if largest_pv_shift > X_MAX_STRIPE_SHIFT {
        // A backbone past 2^21 * 2^11 cells: the serial builder. It ignores `cuts.x_ahead` (the
        // spare's vectors): a build ahead would allocate its table here, on the pool thread. This
        // is unreachable for blastp-fast's compressed backbones (759,375 cells for word size 5,
        // far below 2^32), and the table would have the same values.
        return XFlatStage::Built(serial(cuts.overflow_cell_capacity));
    }
    let Some(matrix) = scores.score_matrix(header.alphabet) else {
        return XFlatStage::Built(None);
    };
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let threads = cuts.threads.max(1);
    let max_appends = cuts
        .max_appends
        .unwrap_or_else(|| x_max_appends(backbone_size));

    let ranges = x_count_ranges(tables, locations.len());
    let sizes = if cuts.chunk_sizes.is_empty() {
        let word = word_size as i64;
        let offsets: i64 = ranges
            .iter()
            .flat_map(|(range, _)| &locations[range.clone()])
            .map(|&(left, right)| (i64::from(right) - word + 2 - i64::from(left)).max(0))
            .sum();
        XChunkSizes::Guided {
            threads,
            offsets: offsets as usize,
        }
    } else {
        XChunkSizes::List(cuts.chunk_sizes)
    };
    let chunks = x_chunks(word_size, &ranges, locations, sizes);
    // LOSAT_X_BPSTRIPEFILL (module header).
    let fill = cuts.stripe_fill && matches!(cuts.passes, XPasses::OneLatch(_));
    // Item SC9i of line NP: the waves (module header); one for a build ahead of its block, whose two
    // passes run apart.
    let wave_letters = if cuts.x_ahead {
        header.alphabet
    } else {
        x_wave_letters(&header, cuts.wave_cells.unwrap_or(X_WAVE_CELLS))
    };
    let waves = header.alphabet.div_ceil(wave_letters);
    // The stripes of a wave of `wave_letters` letters (the backbone with one wave): as many per wave
    // as the one-wave build has (wider ones, a share per wave, left the threads waiting at each
    // wave's end: x1.07 of the wall time at 32 threads, item SC9i).
    let wave_size = if waves == 1 {
        backbone_size
    } else {
        wave_letters * x_top_cells(&header)
    };
    let shift = cuts.stripe_shift.unwrap_or_else(|| {
        let per_wave = (X_PARTS_PER_THREAD * threads).min(X_MAX_PARTS / chunks.len().max(1));
        x_stripe_shift(wave_size, per_wave, largest_pv_shift)
    });
    assert!(
        (largest_pv_shift..=X_MAX_STRIPE_SHIFT).contains(&shift),
        "stripe shift {shift} outside {largest_pv_shift}..={X_MAX_STRIPE_SHIFT}"
    );

    // The first wave's appends, sorted by stripe.
    let appended = AtomicUsize::new(0);
    let wave = x_wave(0, wave_letters, &header, shift);
    let Some((appends, stripe_hits, word_starts)) = x_wave_appends(
        &words,
        query,
        &chunks,
        &wave,
        shift,
        threads,
        max_appends,
        &appended,
        cuts.piece_appends,
        cuts.appends_segment,
        cuts.segments,
        cuts.passes,
        waves > 1,
        None,
    ) else {
        // LOSAT_X_BPLKPFLAT: the reference builder makes this table; the spare's vectors are not
        // kept beside it (module header).
        if let Some(spare) = cuts.spare {
            drop(spare.take());
        }
        return XFlatStage::Built(Some(XFlatBuild::TooLarge));
    };
    XFlatStage::Appends(XFlatAppends {
        header,
        largest_pv_shift,
        ranges,
        chunks,
        word_starts,
        appends,
        shift,
        stripe_hits,
        waves,
        wave_letters,
        appended: appended.into_inner(),
        max_appends,
        piece_appends: cuts.piece_appends,
        appends_segment: cuts.appends_segment,
        scaled_threshold,
        threads,
        overflow_cell_capacity: cuts.overflow_cell_capacity,
        one_latch: matches!(cuts.passes, XPasses::OneLatch(_)),
        fill,
    })
}

// No NCBI counterpart: the order of the build's work (item SC9i of line NP, module header).
/// The cells of one last letter: `A^(word_length - 1)`.
fn x_top_cells(header: &XHeader) -> usize {
    (header.backbone_size - 1) / header.alphabet
}

// No NCBI counterpart: the order of the build's work (item SC9i of line NP, module header).
/// The last letters of a wave: as many as make `ceil(backbone_size / wave_cells)` waves, or the
/// fewest waves with that many letters each.
fn x_wave_letters(header: &XHeader, wave_cells: usize) -> usize {
    header
        .alphabet
        .div_ceil(header.backbone_size.div_ceil(wave_cells.max(1)))
}

// No NCBI counterpart: the order of the build's work (item SC9i of line NP, module header).
/// Wave `wave` of the waves of `letters` last letters each (the last wave fewer), its stripes
/// 2^shift cells from its first.
fn x_wave(wave: usize, letters: usize, header: &XHeader, shift: usize) -> XWave {
    let first = wave * letters;
    let end = (first + letters).min(header.alphabet);
    let top = x_top_cells(header);
    let cells_end = if end == header.alphabet {
        header.backbone_size
    } else {
        end * top
    };
    XWave {
        letters: first as u32..end as u32,
        cells: (first * top) as u32..cells_end as u32,
        stripes: (cells_end - first * top).div_ceil(1 << shift),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (`x_chunks`)
/// The appends pass of one wave: the chunks' appends to the wave's cells (`x_chunk_appends`), sorted
/// by stripe, on `threads` tasks (`passes`: the extra tasks of the one-latch passes run in the first
/// wave's pass), and each of the wave's stripes' hits. `None` when the appends of all chunks of all
/// waves so far (`appended`) pass `max_appends`. The chunks' words' starts as `x_chunk_appends`
/// takes them: from `kept` (each chunk's) where given, else from the query, and with `keep` returned
/// (each chunk's; else none).
#[allow(clippy::too_many_arguments)]
fn x_wave_appends(
    words: &XWords<'_>,
    query: &[u8],
    chunks: &[XChunk],
    wave: &XWave,
    shift: usize,
    threads: usize,
    max_appends: usize,
    appended: &AtomicUsize,
    piece_appends: usize,
    appends_segment: usize,
    segments: Option<&XSegmentPool>,
    passes: XPasses<'_>,
    keep: bool,
    kept: Option<&[Vec<XWordStart>]>,
) -> Option<(XAppends, Vec<usize>, Vec<Vec<XWordStart>>)> {
    let too_large = AtomicBool::new(false);
    let registry: Mutex<Vec<Vec<u32>>> = Mutex::new(Vec::new());
    let chunk_appends = |scratch: &mut XAppendsScratch<'_>, index: usize| {
        let starts = match (kept, keep) {
            (Some(kept), _) => XStarts::Kept(&kept[index]),
            (None, true) => XStarts::Keep,
            (None, false) => XStarts::Make,
        };
        x_chunk_appends(
            words,
            query,
            &chunks[index],
            wave,
            shift,
            max_appends,
            piece_appends,
            appended,
            &too_large,
            starts,
            scratch,
        )
    };
    let chunk_list: Vec<(Vec<XChunkAppends>, Vec<XWordStart>)> = match passes {
        XPasses::Rayon => (0..chunks.len())
            .into_par_iter()
            .map_init(
                || XAppendsScratch::new(&registry, appends_segment, segments),
                chunk_appends,
            )
            .collect(),
        // LOSAT_X_BPLKPFLAT: the extra tasks after the chunks (module header).
        XPasses::OneLatch(extra) => x_one_latch_pass(
            threads,
            chunks.len(),
            extra.map_or(0, |extra| extra.count()),
            || XAppendsScratch::new(&registry, appends_segment, segments),
            chunk_appends,
            |index| {
                if let Some(extra) = extra {
                    extra.run(index);
                }
            },
        ),
        XPasses::Gang(gang) => gang.pass(
            chunks.len(),
            false,
            || XAppendsScratch::new(&registry, appends_segment, segments),
            chunk_appends,
        ),
    };
    if too_large.load(Ordering::Relaxed) {
        return None;
    }

    // Each stripe's hits: its appends in all chunks.
    // Every task has ended: its segments are in their slots.
    let mut starts_kept = Vec::with_capacity(if keep { chunks.len() } else { 0 });
    let mut pieces = Vec::new();
    for (index, (chunk_pieces, chunk_starts)) in chunk_list.into_iter().enumerate() {
        pieces.extend(chunk_pieces.into_iter().map(|mut piece| {
            piece.chunk = index;
            piece
        }));
        if keep {
            starts_kept.push(chunk_starts);
        }
    }
    let appends = XAppends {
        segments: registry
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner),
        pieces,
    };
    let mut stripe_hits = vec![0usize; wave.stripes];
    for piece in &appends.pieces {
        for (hits, bounds) in stripe_hits.iter_mut().zip(piece.stripe_starts.windows(2)) {
            *hits += (bounds[1] - bounds[0]) as usize;
        }
    }
    Some((appends, stripe_hits, starts_kept))
}

// The second pass of `x_flat_parallel_on` (module header): the stripes (count, layout and fill) of
// the appends that `x_flat_appends_on` made, on `spare`'s vectors (with `x_ahead`, a build ahead of
// its block: all of them, and no table without their room), then each later wave's appends and
// stripes (item SC9i of line NP), the capacity test and the PV array. The words of `scores` are
// those of the appends pass (`XWords`, made again for the later waves and the capacity test).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn x_flat_stripes_on<S: CompressedScores, H: XQueryOffset>(
    scores: S,
    appends: XFlatAppends,
    query: &[u8],
    locations: &[(i32, i32)],
    spare: Option<&XFlatSpare<H>>,
    lent: Option<XFlatTableVectors<H>>,
    x_ahead: bool,
    segments: Option<&XSegmentPool>,
) -> Option<XFlatBuild<H>> {
    let XFlatAppends {
        header,
        largest_pv_shift,
        ranges,
        chunks,
        word_starts,
        appends,
        shift,
        stripe_hits,
        waves,
        wave_letters,
        appended,
        max_appends,
        piece_appends,
        appends_segment,
        scaled_threshold,
        threads,
        overflow_cell_capacity,
        one_latch,
        fill,
    } = appends;
    let backbone_size = header.backbone_size;
    let first_total: usize = stripe_hits.iter().sum();
    // The PV shift is `PV_ARRAY_BTS` whatever the counts when `x_largest_pv_shift` is (`finalize`,
    // below): with one wave (whose stripes are whole PV words from cell 0) the stripes set their PV
    // words; else the PV array is made from the starts below.
    let pv_in_stripes = largest_pv_shift == PV_ARRAY_BTS && waves == 1;
    // LOSAT_X_BPLKPAHEAD: a build ahead of its block takes every vector of its table from the spare
    // (its PV array's here) and allocates none; without their room it makes no table and leaves them
    // in the spare (module header). The PV array has at most `(backbone_size >> PV_ARRAY_BTS) + 1`
    // words at any PV shift. Such a build has one wave (`x_flat_appends_on`).
    let mut x_ahead_pv = None;
    if x_ahead {
        assert_eq!(waves, 1, "a build ahead of its block has one wave");
        let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
        x_ahead_pv =
            spare.and_then(|spare| spare.x_room_for(backbone_size + 1, first_total, pv_words));
        x_ahead_pv.as_ref()?;
    }
    // LOSAT_X_BXCLUTPOOL: the vectors BLASTX's pool lent (buffer reuse only, as the spare's).
    let (lent_kept, mut lent_pv) = match lent {
        Some(lent) => (Some((lent.starts, lent.hits)), lent.pv),
        None => (None, Vec::new()),
    };
    let mut pv = if pv_in_stripes {
        match x_ahead_pv.take() {
            Some(buffer) => x_zeroed_in(buffer, (backbone_size >> PV_ARRAY_BTS) + 1),
            None => x_zeroed_on(
                std::mem::take(&mut lent_pv),
                (backbone_size >> PV_ARRAY_BTS) + 1,
            ),
        }
    } else {
        Vec::new()
    };
    // LOSAT_X_BPLKPFLAT: an earlier table's vectors when the spare holds them (module header).
    let (spare_starts, spare_hits) = lent_kept
        .or_else(|| spare.map(XFlatSpare::take))
        .unwrap_or_default();
    let (mut starts, starts_reused) = x_reuse(spare_starts, backbone_size + 1);
    // The hits, written by the stripes in the vector's room past its length (each wave's slots once,
    // the length set after them): no slot is zeroed or read before its stripe writes it, and the
    // stripes, on the pool's threads, touch the pages first. A spare's vector is used when it has
    // the first wave's room (its values are not read). With more than one wave a fresh vector has
    // room for a block that glibc maps whatever its dynamic mmap threshold (`x_mapped_entries`), so
    // that each wave's growth moves the block's pages (`mremap`) instead of copying the hits so far
    // (a block in an arena is copied when the block after it is taken; the threshold rises to
    // 32 MiB once a mapped block is freed).
    let hits_reused = spare_hits.capacity() != 0 && spare_hits.capacity() >= first_total;
    let mut hits = if hits_reused {
        let mut reused = spare_hits;
        reused.clear();
        reused
    } else {
        drop(spare_hits);
        Vec::with_capacity(if waves > 1 {
            first_total.max(x_mapped_entries::<H>())
        } else {
            first_total
        })
    };
    X_LKPFLAT_RECYCLED.fetch_add(u64::from(starts_reused), Ordering::Relaxed);
    X_LKPFLAT_RECYCLED_HITS.fetch_add(u64::from(hits_reused), Ordering::Relaxed);

    // Item SC9i of line NP: the later waves' words (the appends pass's) and the segments of packed
    // appends kept from one wave to the next where the build has no pool of them.
    let words_of_later_waves = if waves > 1 {
        let matrix = scores.score_matrix(header.alphabet)?;
        Some((
            matrix,
            CompressedNeighborInfo::new(&matrix, header.alphabet),
        ))
    } else {
        None
    };
    let words = words_of_later_waves
        .as_ref()
        .map(|(matrix, info)| XWords::new(&header, scaled_threshold, matrix, info));
    let between_waves = XSegmentPool::default();
    let wave_segments = segments.or((waves > 1).then_some(&between_waves));
    let appended = AtomicUsize::new(appended);

    // No NCBI counterpart: whether the stripes sum the overflow cells per count range, which only the
    // capacity test after the waves reads (module header); it changes no value NCBI computes. An
    // append takes at most one cell, so a table of at most `overflow_cell_capacity` appends cannot
    // pass the capacity. With one wave the tables' appends are the appends pass's; with more, a
    // table's appends are known only after the last wave, so the stripes of every wave sum them.
    let need_cells = waves > 1 || {
        let mut range_appends = vec![0usize; ranges.len()];
        for piece in &appends.pieces {
            range_appends[chunks[piece.chunk].range] += piece.len;
        }
        ranges
            .iter()
            .zip(&range_appends)
            .any(|((_, is_table), &total)| *is_table && total > overflow_cell_capacity)
    };
    let mut stripe_counts: Vec<XStripeCounts> = Vec::new();
    let mut first_wave = Some((appends, stripe_hits));
    // The hits of the waves before.
    let mut hits_base = 0usize;
    // The waves' stripes and later appends (`passes`: with more than one wave the one-latch passes are
    // on a gang of the pool's threads, `XGang`, item SC9i of line NP); `false` when the appends pass
    // `max_appends`.
    let mut wave_passes = |passes: XPasses<'_>| -> bool {
        for wave_index in 0..waves {
            let wave = x_wave(wave_index, wave_letters, &header, shift);
            let (appends, stripe_hits) = match first_wave.take() {
                Some(first) => first,
                None => {
                    let words = words.as_ref().expect("the words of the later waves");
                    // The extra tasks (if any) ran in the first wave's appends pass.
                    let passes = match passes {
                        XPasses::OneLatch(_) => XPasses::OneLatch(None),
                        passes => passes,
                    };
                    let Some((made, made_hits, _)) = x_wave_appends(
                        words,
                        query,
                        &chunks,
                        &wave,
                        shift,
                        threads,
                        max_appends,
                        &appended,
                        piece_appends,
                        appends_segment,
                        wave_segments,
                        passes,
                        false,
                        Some(&word_starts),
                    ) else {
                        return false;
                    };
                    // The room for the wave's hits after those of the waves before (the vector grown in
                    // place where the allocator can: glibc moves a mapped block's pages, `mremap`).
                    hits.reserve_exact(made_hits.iter().sum());
                    (made, made_hits)
                }
            };
            let wave_total: usize = stripe_hits.iter().sum();
            debug_assert_eq!(hits.len(), hits_base);

            // The stripes: count, layout and fill.
            let mut parts = Vec::with_capacity(stripe_hits.len());
            let mut starts_rest = &mut starts[wave.cells.start as usize..wave.cells.end as usize];
            let mut hits_rest = &mut hits.spare_capacity_mut()[..wave_total];
            let mut pv_rest = &mut pv[..];
            let mut base = hits_base as u32;
            for (local, &stripe_total) in stripe_hits.iter().enumerate() {
                let first = wave.cells.start as usize + (local << shift);
                let (starts_part, rest) =
                    starts_rest.split_at_mut((wave.cells.end as usize - first).min(1 << shift));
                starts_rest = rest;
                let (hits_part, rest) = hits_rest.split_at_mut(stripe_total);
                hits_rest = rest;
                // A stripe's PV words (one wave: its stripes from cell 0); the last stripe takes the
                // words that are left (the `+ 1` word among them).
                let words = if !pv_in_stripes {
                    0
                } else if local + 1 < wave.stripes {
                    1usize << (shift - PV_ARRAY_BTS)
                } else {
                    pv_rest.len()
                };
                let (pv_part, rest) = pv_rest.split_at_mut(words);
                pv_rest = rest;
                parts.push(XStripe {
                    first,
                    local,
                    base,
                    starts: starts_part,
                    hits: hits_part,
                    pv: pv_part,
                });
                base += stripe_total as u32;
            }
            // The last pass of the gang (its helpers end after it).
            let last = wave_index + 1 == waves;
            let wave_counts: Vec<XStripeCounts> = match (passes, fill) {
                (XPasses::Rayon, _) => parts
                    .into_par_iter()
                    .map(|part| x_stripe_build(part, ranges.len(), &chunks, &appends, need_cells))
                    .collect(),
                (latched, true) => {
                    let parts: Vec<Mutex<Option<XStripe<'_, H>>>> = parts
                        .into_iter()
                        .map(|part| Mutex::new(Some(part)))
                        .collect();
                    latched.latched(
                        threads,
                        parts.len(),
                        last,
                        XFillScratch::default,
                        |scratch, stripe| {
                            let part = parts[stripe]
                                .lock()
                                .unwrap_or_else(PoisonError::into_inner)
                                .take()
                                .expect("each stripe taken once");
                            x_stripe_fill(
                                part,
                                ranges.len(),
                                &chunks,
                                &appends,
                                need_cells,
                                scratch,
                            )
                        },
                    )
                }
                (latched, false) => {
                    // Each stripe's part, taken by the task that takes its index.
                    let parts: Vec<Mutex<Option<XStripe<'_, H>>>> = parts
                        .into_iter()
                        .map(|part| Mutex::new(Some(part)))
                        .collect();
                    latched.latched(
                        threads,
                        parts.len(),
                        last,
                        || (),
                        |_, stripe| {
                            let part = parts[stripe]
                                .lock()
                                .unwrap_or_else(PoisonError::into_inner)
                                .take()
                                .expect("each stripe taken once");
                            x_stripe_build(part, ranges.len(), &chunks, &appends, need_cells)
                        },
                    )
                }
            };
            stripe_counts.extend(wave_counts);
            // SAFETY: the wave's stripes split its slots (their lengths are the stripes' hits, which sum
            // to `wave_total`) and each wrote every slot of its own (`x_stripe_build`, `x_stripe_fill`).
            unsafe { hits.set_len(hits_base + wave_total) };
            // Item SC2d-2: at most as many segments as threads kept for the next build
            // (`XSegmentPool`); item SC9i: and whole for the next wave (`XSegmentPool::keep_whole`).
            match wave_segments {
                Some(pool) if wave_index + 1 < waves => pool.keep_whole(appends.segments, threads),
                Some(pool) => pool.keep(appends.segments, threads),
                None => drop(appends),
            }
            hits_base += wave_total;
        }
        true
    };
    // The gang serves a run of passes (its helpers wait between them). A build of one wave has one pass
    // after its first (its stripes), a one-latch pass of its own as each single pass is: a gang around
    // it held the pool's threads across the build (blastp-fast's maximum RssAnon at 8 threads +3.4 MB
    // over 30 rounds, item SC9i of line NP).
    let fits = match (one_latch, waves > 1) {
        (true, true) => x_with_gang(threads, |gang| wave_passes(XPasses::Gang(gang))),
        (true, false) => wave_passes(XPasses::OneLatch(None)),
        (false, _) => wave_passes(XPasses::Rayon),
    };
    if !fits {
        return Some(XFlatBuild::TooLarge);
    }
    drop(between_waves);
    starts[backbone_size] = hits_base as u32;

    // NCBI's tables (`x_try_build_flat_with_capacity`'s test, at each table's end).
    for (range_index, (range, is_table)) in ranges.iter().enumerate() {
        if !is_table {
            continue;
        }
        let cells: usize = stripe_counts
            .iter()
            .map(|count| count.cells[range_index])
            .sum();
        if cells > overflow_cell_capacity {
            let words_before = range.start > 0;
            // The appends pass's words (made from the same matrix, which had the alphabet's).
            let refused = !words_before || {
                let matrix = scores.score_matrix(header.alphabet)?;
                let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
                let words = XWords::new(&header, scaled_threshold, &matrix, &info);
                x_table_fills_banks(
                    &words,
                    backbone_size,
                    query,
                    &locations[range.clone()],
                    overflow_cell_capacity,
                )
            };
            if refused {
                return Some(XFlatBuild::Refused(x_table_of_range(range_index)));
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1197-1235
    let nonempty_count: usize = stripe_counts.iter().map(|count| count.nonempty).sum();
    let longest_chain = stripe_counts
        .iter()
        .map(|count| count.longest)
        .max()
        .unwrap_or(0) as i32;
    let mut pv_array_bts = PV_ARRAY_BTS;
    if (nonempty_count as f64) <= 0.01 * (backbone_size as f64) {
        pv_array_bts += ilog2(backbone_size / (8 * X_TARGET_PV_BYTES));
    }
    if !pv_in_stripes {
        // LOSAT_X_BPLKPAHEAD: in the spare's buffer for a build ahead of its block (above).
        pv = match x_ahead_pv.take() {
            Some(buffer) => x_zeroed_in(buffer, (backbone_size >> pv_array_bts) + 1),
            None => x_zeroed_on(lent_pv, (backbone_size >> pv_array_bts) + 1),
        };
        for (index, bounds) in starts.windows(2).enumerate() {
            if bounds[1] > bounds[0] {
                compressed_pv_set(&mut pv, index, pv_array_bts);
            }
        }
    }
    assert!(!pv_in_stripes || pv_array_bts == PV_ARRAY_BTS);

    Some(XFlatBuild::Table(XFlatCompressedLookup {
        word_length: header.word_length as i32,
        compressed_alphabet_size: header.alphabet as i32,
        reciprocal_alphabet_size: header.reciprocal_alphabet_size,
        longest_chain,
        pv,
        pv_array_bts: pv_array_bts as i32,
        scaled_compress_table_u8: build_scaled_compress_table_u8(&header.scaled_compress_table),
        starts,
        hits,
    }))
}

#[cfg(test)]
#[path = "x_clut_par_tests.rs"]
mod tests;
