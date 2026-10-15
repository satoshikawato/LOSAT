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
//! Off unless `LOSAT_X_BPCLUT` is set and the search's pool has more than one thread; at one thread
//! the serial builder runs. `LOSAT_X_BPCLUTSHADOW` builds this table too (with more than one thread)
//! and compares it with the serial one.

use super::*;
use crate::utils::threading::SearchPool;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicUsize};

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
/// more than one thread.
pub(super) fn x_try_build_flat_in_pool(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
) -> Option<XFlatBuild> {
    let threads = pool.threads();
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
            },
        )
    })
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
    /// to the stripe's cells take.
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
        // A backbone past 2^21 * 2^11 cells: the serial builder.
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
    let appends: Vec<XChunkAppends> = chunks
        .par_iter()
        .map_init(
            || (Vec::<u32>::new(), Vec::<(i32, u32)>::new()),
            |(targets, per_offset), chunk| {
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
            },
        )
        .collect();
    if too_large.load(Ordering::Relaxed) {
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
    let mut pv = if pv_in_stripes {
        vec![0u32; (backbone_size >> PV_ARRAY_BTS) + 1]
    } else {
        Vec::new()
    };
    let mut starts = vec![0u32; backbone_size + 1];
    starts[backbone_size] = total as u32;
    let mut hits = vec![0i32; total];

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
    let stripe_counts: Vec<XStripeCounts> = parts
        .into_par_iter()
        .map(|part| x_stripe_build(part, shift, backbone_size, ranges.len(), &chunks, &appends))
        .collect();
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
        pv = vec![0u32; (backbone_size >> pv_array_bts) + 1];
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
    ) -> XFlatBuild {
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
                },
            )
        })
        .expect("word size 5")
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
            // Every thread count with every chunk size; the stripe widths in turn.
            for (i, threads) in THREADS.into_iter().enumerate() {
                for (j, chunk_sizes) in CHUNK_SIZES.into_iter().enumerate() {
                    let stripe_shift = SHIFTS[(i + j + round) % SHIFTS.len()];
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
                    )));
                    x_assert_same_flat_tables(&built, &serial);
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
                )));
                x_assert_same_flat_tables(&built, &serial);
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
    // The bound on the appends: `TooLarge` exactly when the appends pass it.
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
                );
                assert_eq!(
                    matches!(built, XFlatBuild::TooLarge),
                    too_large,
                    "round {round}, bound {max_appends}, appends {appends}"
                );
            }
        }
        assert!(refused_seen > 0 && built_seen > 0);
    }
}
