//! EXPERIMENT (LOSAT_X_BPCLUT at more than one thread): the flat compressed lookup table of
//! `x_clut_direct.rs` built by the search's thread pool. The table is the serial builder's
//! (`x_try_build_flat_with_capacity`), slot for slot, at every thread count.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
//! ```c
//!     for (loc = location; loc; loc = loc->next){
//!         Int4 from = loc->ssr->left;
//!         Int4 to = loc->ssr->right - lookup->word_length + 1;
//!
//!         for (offset = from; offset <= to; offset++){
//!             s_CompressedAddWordHits(&info, query->sequence, offset);
//!         }
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 (`add_word_hit`, in
//! the header of `x_clut_direct.rs`)
//! ```c
//!     Int4 num_entries = backbone_cell->num_used;
//!     ...
//!     backbone_cell->num_used++;
//! ```
//!
//! The serial builder makes the reference's appends (a query offset's cells, in the order of
//! `add_word_hits`) for the query offsets of its count ranges (`x_count_ranges`) in order. Here:
//! - chunks (`x_chunks`): each count range's query offsets, in order, cut into contiguous runs, a run
//!   given as pieces of segments (a piece `(a, b + word_length - 1)` holds the offsets a..=b of its
//!   segment, which `x_for_each_offset` visits as it visits them in the segment, with the same
//!   `end <= query.len()` guard). The chunks in order hold the count ranges' offsets in order, so the
//!   chunks' appends one after the other are the serial count pass's appends in its order.
//! - stripes: the backbone cut into ranges of 2^shift cells. Each chunk sorts its appends by stripe,
//!   stably (counting sort), so a stripe's part of a chunk holds the chunk's appends to the stripe's
//!   cells in the chunk's order. A cell lies in one stripe; reading the stripe's parts of the chunks in
//!   chunk order gives the cell's appends in the serial order: its offsets o_0, ..., o_{n-1}, the
//!   append of rank k being the k-th one read.
//! - each stripe's hits: the sum of the chunks' parts for it; a stripe's first hit is the hits of the
//!   stripes before it, so the stripes' slices of the hits are disjoint and in the serial layout's
//!   order (cells in index order).
//! - count (per stripe): each cell's count n (the serial count), and per count range the overflow
//!   cells of the reference's table of all batches that its appends take: an append that finds the
//!   count c takes one when `x_takes_overflow_cell(c)`, the serial count pass's statement, with c its
//!   rank, the count before it in the serial order. Summed over the stripes, a range's cells are
//!   those the serial pass counts for it (`st.cells - cells_before` at the range's end).
//! - layout (per stripe): the cells' starts are the counts' prefix sums from the stripe's first hit,
//!   so the serial starts; the longest chain and the count of non-empty cells are functions of the
//!   counts. The PV shift (`finalize`: `PV_ARRAY_BTS`, or more when at most 1 % of the cells are not
//!   empty) is `PV_ARRAY_BTS` whatever the counts when `ilog2(backbone_size / (8 * kTargetPVBytes))`
//!   is 0 (word size 5): then each stripe sets the PV bits of its cells, in its own words (a stripe
//!   is 2^shift cells, shift >= `PV_ARRAY_BTS`); otherwise the PV array is made after the stripes
//!   from the starts (a cell is not empty when its start is below the next one), by `finalize`'s
//!   statements.
//! - fill (per stripe): the serial fill on the stripe's cells (each append, in the serial order, at
//!   its cell's cursor; the cursors start at the starts), then each cell of n > 5 offsets puts o_k at
//!   `x_emission_position(k, n)`: the serial fill's statements on the same cells.
//! - refusal: the serial builder tests a table when its cells pass the capacity (after a chunk of
//!   appends or at the table's end; the cells only grow, so a table passes the bound at some test
//!   exactly when it passes it at its end), once, and refuses when no segment comes before the
//!   table or the table's own appends fill the banks (`x_table_fills_banks`); the outcome does not
//!   depend on where the test is made. Here, after the stripes, each table's cells at its end are
//!   tested by the same rule, so the build refuses exactly when the serial one does (any table that
//!   refuses makes the outcome `Refused`; the table made is dropped).
//! - a packed append holds the cell's index in its stripe (below 2^16) and the offset less the chunk's
//!   first offset (a chunk's offsets lie within 2^16 of it): nothing is lost.
//! - `TooLarge` (the reference builder runs) when the appends pass `x_max_appends`, as the serial
//!   builder's bound; the reference then decides the table or the refusal.
//! The cuts (chunks, stripes) depend on the thread count; by the above the table does not.
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
//! - extra tasks (`XExtraTasks`: the contexts' Karlin blocks of a blastp block,
//!   `blastp/x_query_blocks.rs`): the one-latch appends pass runs them after its chunks (the indices
//!   past the chunks), each once, also when the build then ends `TooLarge` or refused. They read
//!   nothing the build writes and the build reads nothing they write.
//! - builds ahead (LOSAT_X_BPLKPAHEAD, `XFlatArgs::x_ahead`: the next query block's table, built by a
//!   pool task during the block before it, `blastp/blast_engine/x_lookup_ahead.rs`): the same build,
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
//! copy's part of chunk c taking chunk c's first offset. The overflow cells are counted only when
//! some table's appends pass the capacity: an append takes at most one, so otherwise no table's test
//! can refuse it, and each stripe's counts of them are left 0.
//!
//! On with `LOSAT_X_BPCLUT` (a strict-set switch: on unless `=0`)
//! where the search's pool has more than one thread; at one thread the serial builder runs.

use super::*;
use crate::utils::threading::SearchPool;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::{Mutex, OnceLock, PoisonError};

// No NCBI counterpart: how the build is divided among the pool's threads. Scheduling only.
/// Chunks per thread, and the fewest stripes per thread.
const X_PARTS_PER_THREAD: usize = 8;
/// The fewest query offsets of a chunk (about 65 appends each with blastp-fast's threshold).
const X_MIN_CHUNK_OFFSETS: usize = 256;
/// The most (chunk, stripe) parts: a stripe reads its part of every chunk.
const X_MAX_PARTS: usize = 1 << 16;
// No NCBI counterpart: the 16 bits of a packed append that hold the offset less the chunk's first.
/// A chunk's offsets lie in `base..base + X_OFFSET_SPAN`.
const X_OFFSET_SPAN: i64 = 1 << 16;
// No NCBI counterpart: the 16 bits of a packed append that hold the cell's index in its stripe.
/// The widest stripe: 2^16 cells.
const X_MAX_STRIPE_SHIFT: usize = 16;

/// A contiguous run of the query offsets of one count range.
struct XChunk {
    /// The count range (index into `x_count_ranges`).
    range: usize,
    /// Segments whose offsets are the chunk's, in order.
    pieces: Vec<(i32, i32)>,
    /// The chunk's first offset; its offsets lie in `base..base + X_OFFSET_SPAN`.
    base: i32,
}

/// A chunk's appends sorted by stripe, stably: stripe j's are `packed[stripe_starts[j]..
/// stripe_starts[j + 1]]`, each `(index - (j << shift)) << 16 | (query_offset - base)`.
struct XChunkAppends {
    packed: Vec<u32>,
    stripe_starts: Vec<u32>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
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
    let build = |passes, spare, x_ahead, stripe_fill| {
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
                    passes,
                    spare,
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
        return build(XPasses::Rayon, None, false, false);
    }
    let built = build(
        XPasses::OneLatch(x_flat_args.extra),
        x_flat_args.spare,
        x_flat_args.x_ahead,
        x_bpstripefill_on(),
    );
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
}

impl XFlatSpare {
    // No NCBI counterpart: buffer reuse only (`XFlatSpare`).
    /// The vectors kept (empty ones when none are); the spare then keeps none.
    fn take(&self) -> (Vec<u32>, Vec<i32>) {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1224-1226 (`finalize`)
// ```c
//     pv = lookup->pv = (PV_ARRAY_TYPE *)calloc(
//                               (lookup->backbone_size >> pv_array_bts) + 1,
//                               sizeof(PV_ARRAY_TYPE));
// ```
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
/// The results of `item(scratch, i)` for i in 0..items, in order, made by `threads` tasks (at most
/// one per item or extra task) spawned in one `rayon::scope` of the current pool: each task keeps one
/// `scratch` (`init`) and takes the next index from one counter until none is left; the indices past
/// the items run `extra(i - items)`, `extras` of them.
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
        let (results, next, init, item, extra) = (&results, &next, &init, &item, &extra);
        rayon::scope(|scope| {
            for _ in 0..threads.clamp(1, total) {
                scope.spawn(move |_| {
                    let mut scratch = init();
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= total {
                            break;
                        }
                        if index < items {
                            let stored = results[index].set(item(&mut scratch, index)).is_ok();
                            assert!(stored, "item {index} taken twice");
                        } else {
                            extra(index - items);
                        }
                    }
                });
            }
        });
    }
    results
        .into_iter()
        .map(|result| result.into_inner().expect("every item taken"))
        .collect()
}

/// The parameters of `x_try_build_flat_parallel_with` that the module header's argument leaves free
/// (and the capacity and the bound of the serial builder, for the tests).
struct XParallelCuts<'a> {
    overflow_cell_capacity: usize,
    /// `None`: `x_max_appends`, the serial builder's bound.
    max_appends: Option<usize>,
    threads: usize,
    /// The chunks' sizes in query offsets, in turn; empty: the offsets over `X_PARTS_PER_THREAD`
    /// chunks per thread, at least `X_MIN_CHUNK_OFFSETS` each.
    chunk_sizes: &'a [usize],
    /// The stripes' width 2^shift; `None`: the widest that makes `X_PARTS_PER_THREAD` stripes per
    /// thread, or fewer so that the chunks' parts are at most `X_MAX_PARTS`.
    stripe_shift: Option<usize>,
    /// LOSAT_X_BPLKPFLAT: how the passes are scheduled, and the extra tasks.
    passes: XPasses<'a>,
    /// LOSAT_X_BPLKPFLAT: the vectors of an earlier table for `starts` and `hits` (else fresh ones).
    spare: Option<&'a XFlatSpare>,
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
// ```c
//     if (count <= 0.01 * lookup->backbone_size) {
//         pv_array_bts += ilog2(lookup->backbone_size / (8 * kTargetPVBytes));
//     }
// ```
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
// ```c
//         Int4 from = loc->ssr->left;
//         Int4 to = loc->ssr->right - lookup->word_length + 1;
//
//         for (offset = from; offset <= to; offset++){
// ```
/// The count ranges' offsets cut into chunks (module header) of `sizes[i % sizes.len()]` offsets
/// (chunk i), fewer at a range's end and where the offsets would pass `X_OFFSET_SPAN`.
fn x_chunks(
    word_size: usize,
    ranges: &[(std::ops::Range<usize>, bool)],
    locations: &[(i32, i32)],
    sizes: &[usize],
) -> Vec<XChunk> {
    let word = word_size as i64;
    let mut chunks: Vec<XChunk> = Vec::new();
    for (range_index, (range, _)) in ranges.iter().enumerate() {
        // The chunk being filled (it is chunk `chunks.len()`) and its offsets.
        let mut open: Option<(XChunk, usize)> = None;
        for &(left, right) in &locations[range.clone()] {
            let to = i64::from(right) - word + 1;
            let mut from = i64::from(left);
            while from <= to {
                let size = sizes[chunks.len() % sizes.len()].max(1);
                let joins = matches!(&open, Some((chunk, count))
                    if *count < size
                        && from >= i64::from(chunk.base)
                        && from < i64::from(chunk.base) + X_OFFSET_SPAN);
                if !joins {
                    if let Some((chunk, _)) = open.take() {
                        chunks.push(chunk);
                    }
                    open = Some((
                        XChunk {
                            range: range_index,
                            pieces: Vec::new(),
                            base: from as i32,
                        },
                        0,
                    ));
                }
                let size = sizes[chunks.len() % sizes.len()].max(1);
                let (chunk, count) = open.as_mut().expect("an open chunk");
                // `size - count` is at least 1 (`joins`); past `X_OFFSET_SPAN` it does not matter.
                let take = ((size - *count).min(X_OFFSET_SPAN as usize) as i64)
                    .min(to - from + 1)
                    .min(i64::from(chunk.base) + X_OFFSET_SPAN - from);
                // Offsets from..=from + take - 1: the segment's right end less what follows them.
                chunk
                    .pieces
                    .push((from as i32, (from + take - 1 + word - 1) as i32));
                *count += take as usize;
                from += take;
            }
        }
        if let Some((chunk, _)) = open.take() {
            chunks.push(chunk);
        }
    }
    chunks
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (`x_chunks`) and
// :1075-1121 (`s_CompressedAddWordHits`, `XWords::for_each_word`)
// ```c
//         for (offset = from; offset <= to; offset++){
//             s_CompressedAddWordHits(&info, query->sequence, offset);
//         }
// ```
/// A chunk's appends (the serial count pass's statements), sorted by stripe; empty when the appends
/// of all chunks pass `max_appends` (then `too_large`).
#[allow(clippy::too_many_arguments)]
fn x_chunk_appends(
    words: &XWords<'_>,
    query: &[u8],
    chunk: &XChunk,
    shift: usize,
    stripes: usize,
    max_appends: usize,
    appended: &AtomicUsize,
    too_large: &AtomicBool,
    targets: &mut Vec<u32>,
    per_offset: &mut Vec<(i32, u32)>,
) -> XChunkAppends {
    let empty = || XChunkAppends {
        packed: Vec::new(),
        stripe_starts: vec![0; stripes + 1],
    };
    targets.clear();
    per_offset.clear();
    if too_large.load(Ordering::Relaxed) {
        return empty();
    }
    // The appends not yet added to `appended`.
    let mut reported = 0usize;
    let report = |appends: usize| {
        let total = appended.fetch_add(appends, Ordering::Relaxed) + appends;
        if total > max_appends {
            too_large.store(true, Ordering::Relaxed);
        }
        !too_large.load(Ordering::Relaxed)
    };
    x_for_each_offset(words.word_length, query, &chunk.pieces, |offset, word| {
        let before = targets.len();
        words.for_each_word(word, targets);
        let appends = targets.len() - before;
        if appends != 0 {
            per_offset.push((offset, appends as u32));
        }
        if targets.len() - reported < X_COUNT_CHUNK {
            return true;
        }
        let new = targets.len() - reported;
        reported = targets.len();
        report(new)
    });
    if !report(targets.len() - reported) {
        return empty();
    }

    // Counting sort by stripe, stable.
    let mut stripe_starts = vec![0u32; stripes + 1];
    for &index in targets.iter() {
        stripe_starts[(index as usize >> shift) + 1] += 1;
    }
    for j in 0..stripes {
        stripe_starts[j + 1] += stripe_starts[j];
    }
    let mut cursors = stripe_starts[..stripes].to_vec();
    let mut packed = vec![0u32; targets.len()];
    let mask = (1u32 << shift) - 1;
    let mut at = 0usize;
    for &(query_offset, appends) in per_offset.iter() {
        // Below X_OFFSET_SPAN (`x_chunks`).
        let delta = (query_offset - chunk.base) as u32;
        let end = at + appends as usize;
        for &index in &targets[at..end] {
            let cursor = &mut cursors[index as usize >> shift];
            packed[*cursor as usize] = (index & mask) << 16 | delta;
            *cursor += 1;
        }
        at = end;
    }
    XChunkAppends {
        packed,
        stripe_starts,
    }
}

/// One stripe's part of the table.
struct XStripe<'t> {
    stripe: usize,
    /// The stripe's first hit.
    base: u32,
    starts: &'t mut [u32],
    hits: &'t mut [i32],
    /// The stripe's PV words when the PV shift is known before the counts (`x_largest_pv_shift` is
    /// `PV_ARRAY_BTS`), else empty.
    pv: &'t mut [u32],
}

/// What a stripe's counts give the whole table.
struct XStripeCounts {
    nonempty: usize,
    longest: u32,
    /// The overflow cells of the reference's table of all batches that each count range's appends
    /// to the stripe's cells take (0 where `x_stripe_fill` does not count them: no table's appends
    /// pass the capacity).
    cells: Vec<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 (`add_word_hit`,
// module header of `x_clut_direct.rs`), :1197-1235 (`finalize`) and
// c++/src/algo/blast/core/blast_aascan.c:348-396 (the scan order, module header of
// `x_clut_direct.rs`)
// ```c
//     Int4 num_entries = backbone_cell->num_used;
//     ...
//     backbone_cell->num_used++;
//     ...
//     for (i = 0; i < lookup->backbone_size; i++) {
//         count = lookup->backbone[i].num_used;
//
//         if (count > 0) {
//             /* set the corresponding bit in the pv_array */
//             PV_SET(pv, i, pv_array_bts);
//             longest_chain = MAX(count, longest_chain);
// ```
/// The count, layout and fill of one stripe (module header).
#[allow(clippy::too_many_arguments)]
fn x_stripe_build(
    part: XStripe<'_>,
    shift: usize,
    backbone_size: usize,
    ranges: usize,
    chunks: &[XChunk],
    appends: &[XChunkAppends],
) -> XStripeCounts {
    let XStripe {
        stripe,
        base,
        starts,
        hits,
        pv,
    } = part;
    let first = stripe << shift;
    let parts = || {
        chunks
            .iter()
            .zip(appends)
            .map(move |(chunk, chunk_appends)| {
                (
                    chunk,
                    &chunk_appends.packed[chunk_appends.stripe_starts[stripe] as usize
                        ..chunk_appends.stripe_starts[stripe + 1] as usize],
                )
            })
    };

    // Count: each cell's appends in the serial order; the overflow cells they take.
    let mut counts = vec![0u32; (backbone_size - first).min(1 << shift)];
    let mut cells = vec![0usize; ranges];
    for (chunk, part) in parts() {
        let mut taken = 0usize;
        for &packed in part {
            let count = &mut counts[(packed >> 16) as usize];
            taken += usize::from(x_takes_overflow_cell(*count));
            *count += 1;
        }
        cells[chunk.range] += taken;
    }

    // Layout: the starts (from the stripe's first hit) and the PV bits; each count becomes its
    // cell's cursor, the cell's first slot in the stripe's hits.
    let mut nonempty = 0usize;
    let mut longest = 0u32;
    let mut total = base;
    let mut slot = 0u32;
    for (cell, (start, cursor)) in starts.iter_mut().zip(counts.iter_mut()).enumerate() {
        let count = *cursor;
        *start = total;
        if count > 0 {
            nonempty += 1;
            longest = longest.max(count);
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

    // Fill: the serial fill on the stripe's cells, the appends in order, each at its cell's cursor.
    for (chunk, part) in parts() {
        for &packed in part {
            let cursor = &mut counts[(packed >> 16) as usize];
            hits[*cursor as usize] = chunk.base + (packed & 0xffff) as i32;
            *cursor += 1;
        }
    }
    // A cursor ends at the next cell's first slot; the cells of more than 5 offsets put o_k at
    // `x_emission_position(k, n)` (the serial builder's statements).
    let mut appended: Vec<i32> = Vec::new();
    let mut begin = 0usize;
    for &end in &counts {
        let end = end as usize;
        let n = (end - begin) as u32;
        if n > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1 {
            let cell = &mut hits[begin..end];
            appended.clear();
            appended.extend_from_slice(cell);
            for (k, &query_offset) in appended.iter().enumerate() {
                cell[x_emission_position(k as u32, n) as usize] = query_offset;
            }
        }
        begin = end;
    }
    XStripeCounts {
        nonempty,
        longest,
        cells,
    }
}

/// LOSAT_X_BPSTRIPEFILL: a thread's scratch for `x_stripe_fill`: the copy of the stripe's parts,
/// the cells' counts (then cursors), the cells of more than 5 offsets, a cell's offsets.
type XFillScratch = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<i32>);

// NCBI reference (598d8ae6): as `x_stripe_build` (blast_aalookup.c:787-849, :1197-1235)
// ```c
//     Int4 num_entries = backbone_cell->num_used;
//     ...
//     backbone_cell->num_used++;
// ```
/// LOSAT_X_BPSTRIPEFILL: `x_stripe_build` with the stripe's parts of the chunks read once, in chunk
/// order, into a copy while the cells are counted; the layout and the fill from the copy (module
/// header). The overflow cells are counted only with `need_cells`.
#[allow(clippy::too_many_arguments)]
fn x_stripe_fill(
    part: XStripe<'_>,
    shift: usize,
    backbone_size: usize,
    ranges: usize,
    chunks: &[XChunk],
    appends: &[XChunkAppends],
    need_cells: bool,
    scratch: &mut XFillScratch,
) -> XStripeCounts {
    let XStripe {
        stripe,
        base,
        starts,
        hits,
        pv,
    } = part;
    let (copy, counts, overflowed, appended) = scratch;
    let first = stripe << shift;
    counts.clear();
    counts.resize((backbone_size - first).min(1 << shift), 0);
    copy.clear();
    copy.reserve(hits.len());
    let part_of = |chunk_appends: &XChunkAppends| {
        chunk_appends.stripe_starts[stripe] as usize
            ..chunk_appends.stripe_starts[stripe + 1] as usize
    };

    // Count: each cell's appends in the serial order; the overflow cells they take.
    let mut cells = vec![0usize; ranges];
    for (chunk, chunk_appends) in chunks.iter().zip(appends) {
        let part = &chunk_appends.packed[part_of(chunk_appends)];
        if need_cells {
            let mut taken = 0usize;
            copy.extend(part.iter().map(|&packed| {
                let count = &mut counts[(packed >> 16) as usize];
                taken += usize::from(x_takes_overflow_cell(*count));
                *count += 1;
                packed
            }));
            cells[chunk.range] += taken;
        } else {
            copy.extend(part.iter().map(|&packed| {
                counts[(packed >> 16) as usize] += 1;
                packed
            }));
        }
    }

    // Layout, as `x_stripe_build`'s; the cells of more than 5 offsets noted.
    let mut nonempty = 0usize;
    let mut longest = 0u32;
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
            }
        }
        *cursor = slot;
        slot += count;
    }

    // Fill: the serial fill on the stripe's cells, the appends in order, each at its cell's cursor;
    // the copy holds the chunks' parts one after the other.
    let mut at = 0usize;
    for (chunk, chunk_appends) in chunks.iter().zip(appends) {
        let end = at + part_of(chunk_appends).len();
        for &packed in &copy[at..end] {
            let cursor = &mut counts[(packed >> 16) as usize];
            hits[*cursor as usize] = chunk.base + (packed & 0xffff) as i32;
            *cursor += 1;
        }
        at = end;
    }
    // A cursor ends at its cell's end; the cells of more than 5 offsets put o_k at
    // `x_emission_position(k, n)`.
    for &cell in overflowed.iter() {
        let cell = cell as usize;
        let begin = (starts[cell] - base) as usize;
        let end = counts[cell] as usize;
        let n = (end - begin) as u32;
        let cell_hits = &mut hits[begin..end];
        appended.clear();
        appended.extend_from_slice(cell_hits);
        for (k, &query_offset) in appended.iter().enumerate() {
            cell_hits[x_emission_position(k as u32, n) as usize] = query_offset;
        }
    }
    XStripeCounts {
        nonempty,
        longest,
        cells,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (the passes) and
// :1197-1235 (`finalize`, the statements of `x_try_build_flat_with_capacity`'s layout)
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
/// The parallel build (module header) in the current rayon pool, with the cuts of `cuts`.
fn x_try_build_flat_parallel_with(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    cuts: XParallelCuts<'_>,
) -> Option<XFlatBuild> {
    let header = XHeader::new(word_size)?;
    let backbone_size = header.backbone_size;
    let largest_pv_shift = x_largest_pv_shift(backbone_size);
    if largest_pv_shift > X_MAX_STRIPE_SHIFT {
        // A backbone past 2^21 * 2^11 cells: the serial builder. It ignores `cuts.x_ahead` (the
        // spare's vectors): a build ahead would allocate its table here, on the pool thread. This
        // is unreachable for blastp-fast's compressed backbones (759,375 cells for word size 5,
        // far below 2^32), and the table would have the same values.
        return x_try_build_flat_with_capacity(
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            cuts.overflow_cell_capacity,
            X_COUNT_CHUNK,
        );
    }
    let matrix = build_blosum62_compressed_score_matrix(header.alphabet)?;
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let threads = cuts.threads.max(1);
    let max_appends = cuts
        .max_appends
        .unwrap_or_else(|| x_max_appends(backbone_size));

    let ranges = x_count_ranges(tables, locations.len());
    let default_size;
    let sizes = if cuts.chunk_sizes.is_empty() {
        let word = word_size as i64;
        let offsets: i64 = ranges
            .iter()
            .flat_map(|(range, _)| &locations[range.clone()])
            .map(|&(left, right)| (i64::from(right) - word + 2 - i64::from(left)).max(0))
            .sum();
        default_size = [(offsets as usize)
            .div_ceil(X_PARTS_PER_THREAD * threads)
            .max(X_MIN_CHUNK_OFFSETS)];
        &default_size[..]
    } else {
        cuts.chunk_sizes
    };
    let chunks = x_chunks(word_size, &ranges, locations, sizes);
    // LOSAT_X_BPSTRIPEFILL (module header).
    let fill = cuts.stripe_fill && matches!(cuts.passes, XPasses::OneLatch(_));
    let shift = cuts.stripe_shift.unwrap_or_else(|| {
        let stripes = (X_PARTS_PER_THREAD * threads).min(X_MAX_PARTS / chunks.len().max(1));
        x_stripe_shift(backbone_size, stripes, largest_pv_shift)
    });
    assert!(
        (largest_pv_shift..=X_MAX_STRIPE_SHIFT).contains(&shift),
        "stripe shift {shift} outside {largest_pv_shift}..={X_MAX_STRIPE_SHIFT}"
    );
    let stripes = backbone_size.div_ceil(1 << shift);

    // The chunks' appends, sorted by stripe.
    let appended = AtomicUsize::new(0);
    let too_large = AtomicBool::new(false);
    let chunk_appends = |(targets, per_offset): &mut (Vec<u32>, Vec<(i32, u32)>), chunk| {
        x_chunk_appends(
            &words,
            query,
            chunk,
            shift,
            stripes,
            max_appends,
            &appended,
            &too_large,
            targets,
            per_offset,
        )
    };
    let appends: Vec<XChunkAppends> = match cuts.passes {
        XPasses::Rayon => chunks
            .par_iter()
            .map_init(
                || (Vec::<u32>::new(), Vec::<(i32, u32)>::new()),
                chunk_appends,
            )
            .collect(),
        // LOSAT_X_BPLKPFLAT: the extra tasks after the chunks (module header).
        XPasses::OneLatch(extra) => x_one_latch_pass(
            threads,
            chunks.len(),
            extra.map_or(0, |extra| extra.count()),
            || (Vec::<u32>::new(), Vec::<(i32, u32)>::new()),
            |scratch, index| chunk_appends(scratch, &chunks[index]),
            |index| {
                if let Some(extra) = extra {
                    extra.run(index);
                }
            },
        ),
    };
    if too_large.load(Ordering::Relaxed) {
        // LOSAT_X_BPLKPFLAT: the reference builder makes this table; the spare's vectors are not
        // kept beside it (module header).
        if let Some(spare) = cuts.spare {
            drop(spare.take());
        }
        return Some(XFlatBuild::TooLarge);
    }

    // Each stripe's hits: its appends in all chunks.
    let mut stripe_hits = vec![0usize; stripes];
    for chunk_appends in &appends {
        for (hits, bounds) in stripe_hits
            .iter_mut()
            .zip(chunk_appends.stripe_starts.windows(2))
        {
            *hits += (bounds[1] - bounds[0]) as usize;
        }
    }
    let total: usize = stripe_hits.iter().sum();
    // The PV shift is `PV_ARRAY_BTS` whatever the counts when `x_largest_pv_shift` is (`finalize`,
    // below): the stripes set their PV words; else the PV array is made from the starts below.
    let pv_in_stripes = largest_pv_shift == PV_ARRAY_BTS;
    // LOSAT_X_BPLKPAHEAD: a build ahead of its block takes every vector of its table from the spare
    // (its PV array's here) and allocates none; without their room it makes no table and leaves them
    // in the spare (module header). The PV array has at most `(backbone_size >> PV_ARRAY_BTS) + 1`
    // words at any PV shift.
    let mut x_ahead_pv = None;
    if cuts.x_ahead {
        let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
        x_ahead_pv = cuts
            .spare
            .and_then(|spare| spare.x_room_for(backbone_size + 1, total, pv_words));
        if x_ahead_pv.is_none() {
            return None;
        }
    }
    let mut pv = if pv_in_stripes {
        match x_ahead_pv.take() {
            Some(buffer) => x_zeroed_in(buffer, (backbone_size >> PV_ARRAY_BTS) + 1),
            None => vec![0u32; (backbone_size >> PV_ARRAY_BTS) + 1],
        }
    } else {
        Vec::new()
    };
    // LOSAT_X_BPLKPFLAT: an earlier table's vectors when the spare holds them (module header).
    let (spare_starts, spare_hits) = cuts.spare.map(XFlatSpare::take).unwrap_or_default();
    let (mut starts, starts_reused) = x_reuse(spare_starts, backbone_size + 1);
    let (mut hits, hits_reused) = x_reuse(spare_hits, total);
    X_LKPFLAT_RECYCLED.fetch_add(u64::from(starts_reused), Ordering::Relaxed);
    X_LKPFLAT_RECYCLED_HITS.fetch_add(u64::from(hits_reused), Ordering::Relaxed);
    starts[backbone_size] = total as u32;

    // The stripes: count, layout and fill.
    let mut parts = Vec::with_capacity(stripes);
    let mut starts_rest = &mut starts[..backbone_size];
    let mut hits_rest = &mut hits[..];
    let mut pv_rest = &mut pv[..];
    let mut base = 0u32;
    for (stripe, &stripe_total) in stripe_hits.iter().enumerate() {
        let first = stripe << shift;
        let (starts_part, rest) = starts_rest.split_at_mut((backbone_size - first).min(1 << shift));
        starts_rest = rest;
        let (hits_part, rest) = hits_rest.split_at_mut(stripe_total);
        hits_rest = rest;
        // A stripe's PV words; the last stripe takes the words that are left (the `+ 1` word among
        // them).
        let words = if !pv_in_stripes {
            0
        } else if stripe + 1 < stripes {
            1usize << (shift - PV_ARRAY_BTS)
        } else {
            pv_rest.len()
        };
        let (pv_part, rest) = pv_rest.split_at_mut(words);
        pv_rest = rest;
        parts.push(XStripe {
            stripe,
            base,
            starts: starts_part,
            hits: hits_part,
            pv: pv_part,
        });
        base += stripe_total as u32;
    }
    let stripe_counts: Vec<XStripeCounts> = match cuts.passes {
        XPasses::Rayon => parts
            .into_par_iter()
            .map(|part| x_stripe_build(part, shift, backbone_size, ranges.len(), &chunks, &appends))
            .collect(),
        XPasses::OneLatch(_) if fill => {
            // No NCBI counterpart: whether the stripes count the overflow cells that the capacity test
            // below reads (module header); it changes no value NCBI computes. An append takes at most
            // one: when no table's appends pass the capacity, no table's test can refuse it.
            let mut range_appends = vec![0usize; ranges.len()];
            for (chunk, chunk_appends) in chunks.iter().zip(&appends) {
                range_appends[chunk.range] += chunk_appends.packed.len();
            }
            let need_cells =
                ranges
                    .iter()
                    .zip(&range_appends)
                    .any(|((_, is_table), &range_total)| {
                        *is_table && range_total > cuts.overflow_cell_capacity
                    });
            let parts: Vec<Mutex<Option<XStripe<'_>>>> = parts
                .into_iter()
                .map(|part| Mutex::new(Some(part)))
                .collect();
            x_one_latch_pass(
                threads,
                parts.len(),
                0,
                XFillScratch::default,
                |scratch, stripe| {
                    let part = parts[stripe]
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner)
                        .take()
                        .expect("each stripe taken once");
                    x_stripe_fill(
                        part,
                        shift,
                        backbone_size,
                        ranges.len(),
                        &chunks,
                        &appends,
                        need_cells,
                        scratch,
                    )
                },
                |_| {},
            )
        }
        XPasses::OneLatch(_) => {
            // Each stripe's part, taken by the task that takes its index.
            let parts: Vec<Mutex<Option<XStripe<'_>>>> = parts
                .into_iter()
                .map(|part| Mutex::new(Some(part)))
                .collect();
            x_one_latch_pass(
                threads,
                parts.len(),
                0,
                || (),
                |_, stripe| {
                    let part = parts[stripe]
                        .lock()
                        .unwrap_or_else(PoisonError::into_inner)
                        .take()
                        .expect("each stripe taken once");
                    x_stripe_build(part, shift, backbone_size, ranges.len(), &chunks, &appends)
                },
                |_| {},
            )
        }
    };
    drop(appends);

    // NCBI's tables (`x_try_build_flat_with_capacity`'s test, at each table's end).
    for (range_index, (range, is_table)) in ranges.iter().enumerate() {
        if !is_table {
            continue;
        }
        let cells: usize = stripe_counts
            .iter()
            .map(|count| count.cells[range_index])
            .sum();
        if cells > cuts.overflow_cell_capacity {
            let words_before = range.start > 0;
            let refused = !(words_before
                && !x_table_fills_banks(
                    &words,
                    backbone_size,
                    query,
                    &locations[range.clone()],
                    cuts.overflow_cell_capacity,
                ));
            if refused {
                return Some(XFlatBuild::Refused);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1197-1235
    // ```c
    //     for (i = count = 0; i < lookup->backbone_size; i++) {
    //         if (lookup->backbone[i].num_used)
    //             count++;
    //     }
    //     ...
    //     pv_array_bts = PV_ARRAY_BTS;
    //     if (count <= 0.01 * lookup->backbone_size) {
    //         pv_array_bts += ilog2(lookup->backbone_size / (8 * kTargetPVBytes));
    //     }
    //
    //     pv = lookup->pv = (PV_ARRAY_TYPE *)calloc(
    //                               (lookup->backbone_size >> pv_array_bts) + 1,
    //                               sizeof(PV_ARRAY_TYPE));
    //     ...
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         count = lookup->backbone[i].num_used;
    //
    //         if (count > 0) {
    //             /* set the corresponding bit in the pv_array */
    //             PV_SET(pv, i, pv_array_bts);
    // ```
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
            None => vec![0u32; (backbone_size >> pv_array_bts) + 1],
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
        used_cells: nonempty_count,
        pv,
        pv_array_bts: pv_array_bts as i32,
        scaled_compress_table_u8: build_scaled_compress_table_u8(&header.scaled_compress_table),
        starts,
        hits,
    }))
}

#[cfg(test)]
mod tests {
    use super::super::tests::{queries, XorShift, FEW_LETTERS, LETTERS};
    use super::*;

    fn in_pool<R: Send>(threads: usize, work: impl FnOnce() -> R + Send) -> R {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("test pool")
            .install(work)
    }

    #[allow(clippy::too_many_arguments)]
    fn parallel(
        threads: usize,
        threshold: i32,
        query: &[u8],
        locations: &[(i32, i32)],
        tables: &[std::ops::Range<usize>],
        capacity: usize,
        max_appends: Option<usize>,
        chunk_sizes: &[usize],
        stripe_shift: Option<usize>,
        passes: XPasses<'_>,
        stripe_fill: bool,
        spare: Option<(&XFlatSpare, bool)>,
    ) -> XFlatBuild {
        if let Some((spare, true)) = spare {
            poison(spare);
        }
        in_pool(threads, || {
            x_try_build_flat_parallel_with(
                5,
                threshold,
                query,
                locations,
                tables,
                XParallelCuts {
                    overflow_cell_capacity: capacity,
                    max_appends,
                    threads,
                    chunk_sizes,
                    stripe_shift,
                    passes,
                    spare: spare.map(|(spare, _)| spare),
                    x_ahead: false,
                    stripe_fill,
                },
            )
        })
        .expect("word size 5")
    }

    /// Fills `spare`'s vectors, to their whole capacity, with `u32::MAX` / `i32::MIN`: an entry
    /// that no stripe writes then keeps a value that no build makes (the build takes the vectors'
    /// values as they are, module header).
    fn poison(spare: &XFlatSpare) {
        let mut kept = spare.0.lock().unwrap();
        let (starts, hits) = &mut *kept;
        starts.resize(starts.capacity(), u32::MAX);
        starts.fill(u32::MAX);
        hits.resize(hits.capacity(), i32::MIN);
        hits.fill(i32::MIN);
    }

    // The pass implementations (LOSAT_X_BPLKPFLAT): `par_iter`, and one-latch passes with the stripes
    // of `x_stripe_build` or (LOSAT_X_BPSTRIPEFILL) of `x_stripe_fill`.
    const FORMS: [(XPasses<'static>, bool); 3] = [
        (XPasses::Rayon, false),
        (XPasses::OneLatch(None), false),
        (XPasses::OneLatch(None), true),
    ];

    /// Extra tasks that count how often each one ran.
    struct Counted(Vec<AtomicUsize>);

    impl Counted {
        fn new(count: usize) -> Self {
            Self((0..count).map(|_| AtomicUsize::new(0)).collect())
        }

        fn assert_each_ran(&self, times: usize, what: &str) {
            for (index, runs) in self.0.iter().enumerate() {
                assert_eq!(
                    runs.load(Ordering::Relaxed),
                    times,
                    "{what}: extra task {index}"
                );
            }
        }
    }

    impl XExtraTasks for Counted {
        fn count(&self) -> usize {
            self.0.len()
        }

        fn run(&self, index: usize) {
            self.0[index].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// A spare of garbage: `starts` with room for a word-size-5 table (759,376 cells and one), `hits`
    /// of random length.
    fn garbage_spare(rng: &mut XorShift) -> XFlatSpare {
        let spare = XFlatSpare::default();
        *spare.0.lock().unwrap() = (
            vec![u32::MAX - 1; 759_377 + rng.below(9)],
            vec![i32::MIN + 1; rng.below(4000)],
        );
        spare
    }

    /// Reserves room in a spare of garbage for a word-size-5 table's `starts` and for 65,536 hits,
    /// so that `poison` fills whole vectors and the build never allocates a fresh (zeroed) vector
    /// where a spare is poisoned.
    fn give_room(spare: &XFlatSpare) {
        let mut kept = spare.0.lock().unwrap();
        let (starts, hits) = &mut *kept;
        starts.reserve(759_377);
        hits.reserve(1 << 16);
    }

    /// The capacities of the spare's vectors.
    fn spare_capacities(spare: &XFlatSpare) -> (usize, usize) {
        let kept = spare.0.lock().unwrap();
        (kept.0.capacity(), kept.1.capacity())
    }

    // The cuts of the tests: thread counts, chunk sizes (even and uneven), stripe widths (the
    // narrowest, 2^5 cells, a middle one and the widest).
    const THREADS: [usize; 4] = [1, 2, 3, 7];
    const CHUNK_SIZES: [&[usize]; 4] = [&[], &[1], &[13, 1, 97, 2, 5, 400], &[usize::MAX]];
    const SHIFTS: [Option<usize>; 4] = [None, Some(5), Some(11), Some(16)];

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
    // ```c
    //     for (loc = location; loc; loc = loc->next){
    //         Int4 from = loc->ssr->left;
    //         Int4 to = loc->ssr->right - lookup->word_length + 1;
    // ```
    // Random queries (few letters and low thresholds, so that cells overflow; all letters, with
    // letters outside the compressed alphabet), SEG-masked runs (several segments per query), several
    // NCBI tables and split batches (segments outside every table), threshold 0: at 1, 2, 3 and 7
    // threads, with even and uneven chunks and narrow to wide stripes, the parallel table is the
    // serial one slot for slot (and the reference's, cell by cell in scan order).
    #[test]
    fn parallel_table_equals_serial_on_random_queries() {
        let mut rng = XorShift(0x510e_527f_ade6_82d1);
        let mut overflow_seen = false;
        for round in 0..10 {
            let (letters, threshold): (&[u8], i32) = match round % 5 {
                0 => (&FEW_LETTERS[..], 1100),
                1 => (&FEW_LETTERS[..], 1600),
                2 => (&LETTERS[..], 1100),
                3 => (&[1u8, 11, 17, 21, 24, 10, 3][..], 2000),
                _ => (&FEW_LETTERS[..], 0),
            };
            let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(300)).collect();
            let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
            let mut tables = Vec::new();
            let mut q = 0;
            while q < lengths.len() {
                let end = (q + 1 + rng.below(2)).min(lengths.len());
                if rng.below(3) != 0 {
                    tables.push(first[q]..first[end]);
                }
                q = end;
            }
            let serial = flat(x_try_build_flat_compressed_lookup(
                5, threshold, &query, &locations, &tables,
            ));
            let reference =
                try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                    .expect("word size 5")
                    .expect("not refused");
            x_assert_same_compressed_table(&serial, &reference);
            overflow_seen |= serial.longest_chain > 5;
            // Every thread count with every chunk size and every form (`FORMS`); the stripe
            // widths in turn.
            for (i, threads) in THREADS.into_iter().enumerate() {
                for (j, chunk_sizes) in CHUNK_SIZES.into_iter().enumerate() {
                    let stripe_shift = SHIFTS[(i + j + round) % SHIFTS.len()];
                    for (passes, stripe_fill) in FORMS {
                        let built = flat(Some(parallel(
                            threads,
                            threshold,
                            &query,
                            &locations,
                            &tables,
                            NCBI_OVERFLOW_CELL_CAPACITY,
                            None,
                            chunk_sizes,
                            stripe_shift,
                            passes,
                            stripe_fill,
                            None,
                        )));
                        x_assert_same_flat_tables(&built, &serial);
                    }
                }
            }
        }
        assert!(overflow_seen, "cells with overflow lists");
    }

    fn flat(build: Option<XFlatBuild>) -> XFlatCompressedLookup {
        match build.expect("word size 5") {
            XFlatBuild::Table(table) => table,
            XFlatBuild::Refused => panic!("refused"),
            XFlatBuild::TooLarge => panic!("too large"),
        }
    }

    // A query longer than a chunk's offset span (2^16) in one segment, and in several: the chunks
    // are cut where the span ends, and the table is the serial one.
    #[test]
    fn parallel_table_equals_serial_past_the_offset_span() {
        let mut rng = XorShift(0x9b05_688c_2b3e_6c1f);
        for masked in [0, 3] {
            let (query, locations, first) = queries(&mut rng, &[70_000, 900], &LETTERS, masked);
            let tables = [first[0]..first[2]];
            let serial = flat(x_try_build_flat_compressed_lookup(
                5, 2000, &query, &locations, &tables,
            ));
            for (threads, chunk_sizes) in [(2, &[usize::MAX][..]), (3, &[][..]), (7, &[40_000][..])]
            {
                for (passes, stripe_fill) in FORMS {
                    let built = flat(Some(parallel(
                        threads,
                        2000,
                        &query,
                        &locations,
                        &tables,
                        NCBI_OVERFLOW_CELL_CAPACITY,
                        None,
                        chunk_sizes,
                        None,
                        passes,
                        stripe_fill,
                        None,
                    )));
                    x_assert_same_flat_tables(&built, &serial);
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
    // ```c
    //     if (lookup->curr_overflow_cell ==
    //                         COMPRESSED_OVERFLOW_CELLS_IN_BANK) {
    //         /* need a new bank */
    //         Int4 bank_idx = lookup->curr_overflow_bank + 1;
    // ```
    // The refusal at small capacities (around each table's own cells and its cells of the table of
    // all batches, as the serial builder's test): at 2, 3 and 7 threads and uneven chunks the
    // parallel builder refuses exactly when the reference does, and otherwise builds its table.
    // The bound on the appends: `TooLarge` exactly when the appends pass it. Both pass
    // implementations; the one-latch passes with extra tasks, each of which runs exactly once, also
    // when the build is refused or `TooLarge`, and on a spare of garbage (filled with `poison`'s
    // pattern or not), which a refused build consumes and a `TooLarge` build drops; the next build on
    // a spare of garbage is still exact.
    #[test]
    fn parallel_refusal_and_bound_are_the_serial_ones() {
        let mut rng = XorShift(0x1f83_d9ab_fb41_bd6b);
        let mut refused_seen = 0usize;
        let mut built_seen = 0usize;
        for round in 0..8 {
            let threshold = [1100, 1600, 2000][round % 3];
            let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 5 + rng.below(350)).collect();
            let (query, locations, first) = queries(&mut rng, &lengths, &FEW_LETTERS, round % 2);
            let mut tables = Vec::new();
            let mut q = 0;
            while q < lengths.len() {
                let end = (q + 1 + rng.below(3)).min(lengths.len());
                if rng.below(4) != 0 {
                    tables.push(first[q]..first[end]);
                }
                q = end;
            }
            let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
            let mut capacities = vec![1usize];
            for table in &tables {
                let mut own = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
                own.threshold = threshold;
                own.overflow_cell_capacity = usize::MAX;
                own.add_neighboring_words(&matrix, &query, &locations[table.clone()]);
                let cells = own.overflow_cells.len();
                capacities.extend([cells.saturating_sub(1), cells, cells + 1, 2 * cells]);
            }
            capacities.retain(|&capacity| capacity > 0);
            for capacity in capacities {
                let reference = try_build_blosum62_compressed_lookup_with_capacity(
                    5, threshold, &query, &locations, &tables, capacity,
                )
                .expect("word size 5");
                for (threads, chunk_sizes) in [(2, &[][..]), (3, &[7, 1, 50][..]), (7, &[1][..])] {
                    for (one_latch, stripe_fill) in [(false, false), (true, false), (true, true)] {
                        let extra = Counted::new(1 + rng.below(40));
                        let passes = if one_latch {
                            XPasses::OneLatch(Some(&extra))
                        } else {
                            XPasses::Rayon
                        };
                        let spare = garbage_spare(&mut rng);
                        let poison = rng.below(2) == 0;
                        if poison {
                            give_room(&spare);
                        }
                        match (
                            &reference,
                            parallel(
                                threads,
                                threshold,
                                &query,
                                &locations,
                                &tables,
                                capacity,
                                None,
                                chunk_sizes,
                                None,
                                passes,
                                stripe_fill,
                                one_latch.then_some((&spare, poison)),
                            ),
                        ) {
                            (Err(()), XFlatBuild::Refused) => refused_seen += 1,
                            (Ok(reference), XFlatBuild::Table(built)) => {
                                x_assert_same_compressed_table(&built, reference);
                                built_seen += 1;
                            }
                            (reference, _) => panic!(
                                "round {round}, capacity {capacity}, threads {threads}: reference refused {}",
                                reference.is_err()
                            ),
                        }
                        extra.assert_each_ran(
                            usize::from(one_latch),
                            &format!("round {round}, capacity {capacity}, threads {threads}"),
                        );
                        assert_eq!(spare_capacities(&spare) == (0, 0), one_latch);
                    }
                }
            }
            // The bound: below, at and above the number of appends (the serial table's hits).
            let serial = flat(x_try_build_flat_compressed_lookup(
                5, threshold, &query, &locations, &tables,
            ));
            let appends = serial.hits.len();
            for (max_appends, too_large) in [
                (appends.saturating_sub(1), appends > 0),
                (appends, false),
                (appends + 1, false),
            ] {
                for (one_latch, stripe_fill) in [(false, false), (true, false), (true, true)] {
                    let extra = Counted::new(1 + rng.below(40));
                    let passes = if one_latch {
                        XPasses::OneLatch(Some(&extra))
                    } else {
                        XPasses::Rayon
                    };
                    let spare = garbage_spare(&mut rng);
                    let poison = rng.below(2) == 0;
                    if poison {
                        give_room(&spare);
                    }
                    let built = parallel(
                        3,
                        threshold,
                        &query,
                        &locations,
                        &tables,
                        NCBI_OVERFLOW_CELL_CAPACITY,
                        Some(max_appends),
                        &[5],
                        None,
                        passes,
                        stripe_fill,
                        one_latch.then_some((&spare, poison)),
                    );
                    assert_eq!(
                        matches!(built, XFlatBuild::TooLarge),
                        too_large,
                        "round {round}, bound {max_appends}, appends {appends}"
                    );
                    if let XFlatBuild::Table(built) = built {
                        x_assert_same_flat_tables(&built, &serial);
                    }
                    extra.assert_each_ran(
                        usize::from(one_latch),
                        &format!("round {round}, bound {max_appends}"),
                    );
                    // A `TooLarge` build drops the spare's vectors, a table takes them.
                    assert_eq!(spare_capacities(&spare) == (0, 0), one_latch);
                }
            }
        }
        assert!(refused_seen > 0 && built_seen > 0);
    }

    // LOSAT_X_BPLKPFLAT: a build on vectors of garbage (`u32::MAX`, `i32::MIN`, random values) longer
    // or shorter than the table's, with and without the room for it, and then on the vectors that
    // build gives back, is the fresh build slot for slot, with every form and at 1, 2,
    // 3 and 7 threads; where the spare's `hits` has the room, the table holds that buffer.
    #[test]
    fn builds_on_recycled_vectors_equal_fresh_builds() {
        let mut rng = XorShift(0x6a09_e667_f3bc_c908);
        let spare = XFlatSpare::default();
        let mut reused_seen = 0usize;
        for round in 0..6 {
            let (letters, threshold): (&[u8], i32) = match round % 3 {
                0 => (&FEW_LETTERS[..], 1100),
                1 => (&LETTERS[..], 1100),
                _ => (&FEW_LETTERS[..], 1600),
            };
            let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(400)).collect();
            let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
            let tables = [first[0]..first[lengths.len()]];
            let serial = flat(x_try_build_flat_compressed_lookup(
                5, threshold, &query, &locations, &tables,
            ));
            let cells = serial.starts.len();
            let total = serial.hits.len();
            for case in 0..3 {
                for threads in THREADS {
                    for (passes, stripe_fill) in FORMS {
                        let mut random = |len: usize| -> Vec<u32> {
                            (0..len)
                                .map(|_| rng.below(u32::MAX as usize) as u32)
                                .collect()
                        };
                        let garbage: (Vec<u32>, Vec<i32>) = match case {
                            // Longer than the table's.
                            0 => (vec![u32::MAX; cells + 7], vec![i32::MIN; total + 100]),
                            // Random values; `hits` shorter, with room to grow in place.
                            1 => {
                                let starts = random(cells);
                                let mut hits = Vec::with_capacity(2 * total + 1);
                                hits.extend(
                                    random(total / 2).into_iter().map(|value| value as i32),
                                );
                                (starts, hits)
                            }
                            // No room: fresh vectors.
                            _ => (Vec::new(), vec![i32::MIN; total.saturating_sub(1)]),
                        };
                        *spare.0.lock().unwrap() = garbage;
                        // The garbage, then the vectors the first build gives back.
                        for build in 0..2 {
                            let (given, room) = {
                                let kept = spare.0.lock().unwrap();
                                (kept.1.as_ptr() as usize, kept.1.capacity() >= total)
                            };
                            let built = flat(Some(parallel(
                                threads,
                                threshold,
                                &query,
                                &locations,
                                &tables,
                                NCBI_OVERFLOW_CELL_CAPACITY,
                                None,
                                &[],
                                None,
                                passes,
                                stripe_fill,
                                Some((&spare, build == 1)),
                            )));
                            x_assert_same_flat_tables(&built, &serial);
                            if room && total > 0 {
                                assert_eq!(
                                    built.hits.as_ptr() as usize,
                                    given,
                                    "round {round}, case {case}, build {build}: the spare's buffer"
                                );
                                reused_seen += 1;
                            }
                            let kept = spare.0.lock().unwrap();
                            assert!(
                                kept.0.is_empty() && kept.1.is_empty(),
                                "the build took them"
                            );
                            drop(kept);
                            spare.give(built);
                        }
                    }
                }
            }
        }
        assert!(reused_seen > 0);
    }

    // LOSAT_X_BPLKPFLAT: the extra tasks run in the one-latch appends pass, each exactly once, at
    // every thread count and with any number of chunks (none included), and not in the serial
    // builder that a one-thread pool runs (`x_build_flat_in_pool`); the table does not change.
    #[test]
    fn extra_tasks_run_once_in_the_parallel_build_only() {
        let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
        let lengths = [120, 7, 300];
        let (query, locations, first) = queries(&mut rng, &lengths, &FEW_LETTERS, 1);
        let tables = [first[0]..first[3]];
        let serial = flat(x_try_build_flat_compressed_lookup(
            5, 1100, &query, &locations, &tables,
        ));
        for threads in THREADS {
            for chunk_sizes in CHUNK_SIZES {
                for (count, stripe_fill) in [0, 1, 3, 64]
                    .into_iter()
                    .flat_map(|count| [(count, false), (count, true)])
                {
                    let extra = Counted::new(count);
                    let built = flat(Some(parallel(
                        threads,
                        1100,
                        &query,
                        &locations,
                        &tables,
                        NCBI_OVERFLOW_CELL_CAPACITY,
                        None,
                        chunk_sizes,
                        None,
                        XPasses::OneLatch(Some(&extra)),
                        stripe_fill,
                        None,
                    )));
                    x_assert_same_flat_tables(&built, &serial);
                    extra.assert_each_ran(1, &format!("threads {threads}, {count} tasks"));
                }
            }
        }
        // No offsets: no chunk, the tasks still run once.
        for stripe_fill in [false, true] {
            let extra = Counted::new(5);
            let none = flat(Some(parallel(
                3,
                1100,
                &query,
                &[],
                &[],
                NCBI_OVERFLOW_CELL_CAPACITY,
                None,
                &[],
                None,
                XPasses::OneLatch(Some(&extra)),
                stripe_fill,
                None,
            )));
            assert!(none.hits.is_empty());
            extra.assert_each_ran(1, "no chunks");
        }
        // The entry of the search: the pool's build at 3 threads, the serial builder at one thread
        // (no task runs; the caller runs them), and nothing with the switch off.
        for (threads, mode, runs) in [
            (1, XLkpFlat::On, 0),
            (3, XLkpFlat::On, 1),
            (3, XLkpFlat::Off, 0),
        ] {
            let extra = Counted::new(9);
            let built = crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                Ok(super::super::x_build_flat_in_pool(
                    5,
                    1100,
                    &query,
                    &locations,
                    &tables,
                    pool,
                    XFlatArgs {
                        mode,
                        spare: None,
                        extra: Some(&extra),
                        x_ahead: false,
                    },
                ))
            })
            .expect("pool");
            x_assert_same_flat_tables(&flat(built), &serial);
            extra.assert_each_ran(runs, &format!("threads {threads}, {mode:?}"));
        }
    }

    /// LOSAT_X_BPLKPAHEAD: the one-latch build ahead of its block (`x_ahead`) on `spare`.
    #[allow(clippy::too_many_arguments)]
    fn ahead(
        threads: usize,
        threshold: i32,
        query: &[u8],
        locations: &[(i32, i32)],
        tables: &[std::ops::Range<usize>],
        spare: Option<&XFlatSpare>,
        poisoned: bool,
        stripe_fill: bool,
    ) -> Option<XFlatBuild> {
        if let (Some(spare), true) = (spare, poisoned) {
            poison(spare);
        }
        in_pool(threads, || {
            x_try_build_flat_parallel_with(
                5,
                threshold,
                query,
                locations,
                tables,
                XParallelCuts {
                    overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
                    max_appends: None,
                    threads,
                    chunk_sizes: &[],
                    stripe_shift: None,
                    passes: XPasses::OneLatch(None),
                    stripe_fill,
                    spare,
                    x_ahead: true,
                },
            )
        })
    }

    /// The addresses of the spare's three vectors (`starts`, `hits`, `pv`).
    fn spare_addresses(spare: &XFlatSpare) -> [usize; 3] {
        let kept = spare.0.lock().unwrap();
        let pv = spare.1.lock().unwrap();
        [
            kept.0.as_ptr() as usize,
            kept.1.as_ptr() as usize,
            pv.as_ptr() as usize,
        ]
    }

    // LOSAT_X_BPLKPAHEAD: a build ahead of its block takes `starts`, `hits` and the PV array from the
    // spare that `x_make_room` made ready (garbage in them, the hits shorter than the table's) and
    // allocates none of them: the table holds the spare's three buffers and equals the fresh serial
    // build slot for slot, at 1, 2, 3 and 7 threads, poisoned or not. A spare one entry short of the
    // room for any of the three (or no spare) gives no table and keeps its vectors untouched.
    #[test]
    fn builds_ahead_take_every_vector_from_the_spare() {
        let mut rng = XorShift(0x3c6e_f372_fe94_f82b);
        let mut built_seen = 0usize;
        for round in 0..4 {
            let (letters, threshold): (&[u8], i32) = match round % 2 {
                0 => (&FEW_LETTERS[..], 1100),
                _ => (&LETTERS[..], 1100),
            };
            let lengths: Vec<usize> = (0..2 + rng.below(4)).map(|_| 1 + rng.below(400)).collect();
            let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
            let tables = [first[0]..first[lengths.len()]];
            let serial = flat(x_try_build_flat_compressed_lookup(
                5, threshold, &query, &locations, &tables,
            ));
            let backbone_size = serial.starts.len() - 1;
            let total = serial.hits.len();
            let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
            for threads in THREADS {
                for (poison, stripe_fill) in
                    [(false, false), (false, true), (true, false), (true, true)]
                {
                    // The room, with garbage in it.
                    let spare = XFlatSpare::default();
                    let mut hits = Vec::with_capacity(total + 5);
                    hits.extend((0..total / 2).map(|_| rng.below(1 << 20) as i32));
                    *spare.0.lock().unwrap() = (vec![u32::MAX - 3; backbone_size + 1], hits);
                    *spare.1.lock().unwrap() = vec![u32::MAX; pv_words];
                    let block = XFlatSpare::default();
                    assert_eq!(spare.x_make_room(&block, backbone_size, total, total), 0);
                    let given = spare_addresses(&spare);
                    let built = flat(ahead(
                        threads,
                        threshold,
                        &query,
                        &locations,
                        &tables,
                        Some(&spare),
                        poison,
                        stripe_fill,
                    ));
                    x_assert_same_flat_tables(&built, &serial);
                    let held = [
                        built.starts.as_ptr() as usize,
                        built.hits.as_ptr() as usize,
                        built.pv.as_ptr() as usize,
                    ];
                    assert_eq!(held, given, "round {round}, threads {threads}");
                    assert_eq!(spare_capacities(&spare), (0, 0));
                    assert_eq!(spare.1.lock().unwrap().capacity(), 0);
                    built_seen += 1;
                    // One entry short of the room for each vector in turn, and no spare.
                    for short in 0..3 {
                        if short == 1 && total == 0 {
                            continue;
                        }
                        let spare = XFlatSpare::default();
                        *spare.0.lock().unwrap() = (
                            Vec::with_capacity(backbone_size + usize::from(short != 0)),
                            Vec::with_capacity(total + usize::from(short != 1) - 1),
                        );
                        *spare.1.lock().unwrap() =
                            Vec::with_capacity(pv_words - usize::from(short == 2));
                        let kept = spare_addresses(&spare);
                        let none = ahead(
                            threads,
                            threshold,
                            &query,
                            &locations,
                            &tables,
                            Some(&spare),
                            poison,
                            stripe_fill,
                        );
                        assert!(none.is_none(), "round {round}, short {short}: no table");
                        assert_eq!(spare_addresses(&spare), kept, "the spare keeps them");
                    }
                    let none = ahead(
                        threads,
                        threshold,
                        &query,
                        &locations,
                        &tables,
                        None,
                        poison,
                        stripe_fill,
                    );
                    assert!(none.is_none());
                }
            }
        }
        assert!(built_seen > 0);
    }

    // LOSAT_X_BPLKPAHEAD: `x_make_room` keeps the spare's own vectors when they have the room, else
    // moves the block spare's (the same buffers), and replaces each vector short of its room by an
    // empty one of the asked capacity; it writes nothing into them.
    #[test]
    fn make_room_keeps_or_moves_the_vectors_and_allocates_only_what_is_short() {
        let backbone_size = 759_375;
        let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
        // Nothing kept anywhere: three new vectors.
        let ahead_spare = XFlatSpare::default();
        let block = XFlatSpare::default();
        assert_eq!(
            ahead_spare.x_make_room(&block, backbone_size, 1000, 1500),
            3
        );
        {
            let kept = ahead_spare.0.lock().unwrap();
            assert!(kept.0.is_empty() && kept.0.capacity() >= backbone_size + 1);
            assert!(kept.1.is_empty() && kept.1.capacity() >= 1500);
            assert!(ahead_spare.1.lock().unwrap().capacity() >= pv_words);
        }
        // Enough room: nothing allocated, the same buffers.
        let before = spare_addresses(&ahead_spare);
        assert_eq!(
            ahead_spare.x_make_room(&block, backbone_size, 1400, 3000),
            0
        );
        assert_eq!(spare_addresses(&ahead_spare), before);
        // `hits` short: only it is replaced (empty, with the grown capacity).
        assert_eq!(
            ahead_spare.x_make_room(&block, backbone_size, 1600, 3200),
            1
        );
        let after = spare_addresses(&ahead_spare);
        assert_eq!([after[0], after[2]], [before[0], before[2]]);
        assert!(ahead_spare.0.lock().unwrap().1.capacity() >= 3200);
        // An empty spare takes the block spare's vectors (the same buffers, contents kept).
        let empty = XFlatSpare::default();
        *block.0.lock().unwrap() = (vec![7; backbone_size + 1], vec![-7; 5000]);
        let moved = spare_addresses(&block);
        assert_eq!(empty.x_make_room(&block, backbone_size, 5000, 7500), 1);
        let taken = spare_addresses(&empty);
        assert_eq!([taken[0], taken[1]], [moved[0], moved[1]]);
        assert_eq!(spare_capacities(&block), (0, 0));
        assert_eq!(empty.0.lock().unwrap().1.len(), 5000);
    }
}
