//! EXPERIMENT (LOSAT_X_BPCLUT; BPCLUT is a strict-set switch): the compressed lookup table of
//! blastp-fast (`CompressedAaLookupTable`) built flat, in three passes (count, layout, fill): each
//! cell's hits lie in one contiguous slice, already in the order NCBI's compressed scan copies them,
//! and the scans read a cell as that slice.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:348-396
//!
//! The reference (`try_build_blosum62_compressed_lookup`) runs `add_word_hits` for every query offset
//! of every lookup segment, in order: each call appends the offset to the cells of its words (the word
//! itself, then its neighbors in the order of the `add_word_hits_core` recursion) with `add_word_hit`,
//! which keeps a cell's first offsets in the backbone and the rest in a list of overflow cells. Here:
//! - count pass: the same appends in the same order (with `LOSAT_X_CLUTPAIRS`, a query offset's
//!   cells in another order, below) (`XWords::for_each_word`: the statements of
//!   `add_word_hits` and `add_word_hits_core`, the index of a neighbor formed as its letters are
//!   chosen, `sum c[i] * A^i` as `compute_encoded_compressed_index` forms it), recorded as the cell of
//!   each append and the number of appends of each query offset; then each cell's count goes up by
//!   one per recorded append, in order.
//! - layout: the cells' slices one after the other in index order (`starts`, the prefix sums of the
//!   counts), and the PV array and `longest_chain` from the counts by `finalize`'s statements.
//! - fill: the recorded appends again in order, each written at its cell's cursor (the cell's start,
//!   moved on by one per append), so that a cell's slice holds its offsets in the order of their
//!   appends; then the offset of rank k (0-based) of a cell whose count n is above 5 is moved to slot
//!   `x_emission_position(k, n)` (for n <= 5 that slot is k), in place (`x_to_scan_order`).
//!
//! The permutation. Let a cell take the offsets o_0, ..., o_{n-1} in this order (o_k is appended when
//! `num_used` is k). By `add_word_hit` (above): o_0 is `query_offset` and o_1..o_4 are
//! `payload.query_offsets[0..4]`; when n >= 6, o_5 makes the first overflow cell, which takes o_3, o_4
//! and o_5 in slots 0, 1, 2 while o_1, o_2 stay as `overflow_list.query_offsets[0..2]`; each o_k with
//! k >= 6 goes to slot (k - 3) & 3 of the head, and a new head is put in front when that slot is 0.
//! So for n >= 6 the offset o_k with k >= 3 lies in the g-th overflow cell of the cell, g = (k - 3) >> 2,
//! at slot s = (k - 3) & 3, the list runs from the newest cell (g = G - 1) to the first (g = 0, whose
//! `next` is NULL), G = ceil((n - 3) / 4) = floor(n / 4) cells, and the newest holds
//! m = n - 3 - 4 (G - 1) offsets, 1 <= m <= 4. The scan (above) writes o_0 at `dest[0]`; for n <= 5 it
//! writes o_k at `dest[k]`. For n >= 6 it writes o_1, o_2 at `dest[1]`, `dest[2]`, then
//! f = (n - 3) & 3 slots of the head; n - 3 = m (mod 4), so f = m when f != 0, and the walk then moves
//! to the next cell, while for f = 0 (m = 4) the walk starts at the head: either way the head's m
//! slots come first, then cells G - 2, ..., 0 with 4 slots each. Hence the emission position of o_k:
//! - k when n <= 5 or k < 3;
//! - 3 + s when g = G - 1 (the head);
//! - 3 + m + 4 (G - 2 - g) + s = n - 4 (g + 1) + s when g < G - 1.
//!
//! The head's positions are 3..3 + m and cell g's are n - 4 (g + 1)..n - 4 g, so the positions are a
//! permutation of 0..n. For n = 6 and 7 there is one cell and every position is k.
//!
//! Identity. The count pass makes the reference's appends in the reference's order:
//! `for_each_word` returns at a letter without a compressed letter, appends the word itself when
//! `threshold == 0 || self_score < threshold` at the index of `compute_compressed_index` (when its
//! `skip` is 0, as `add_unencoded`), stops when `threshold == 0`, and otherwise enumerates the
//! neighbors with the recursion's statements and break tests (its test of `overflow_banks_exhausted`
//! never fires on the reference's table of all batches, which has no limit); the loop over segments and
//! offsets is `add_location_words`' (the same `end <= query.len()` guard). So each cell gets the
//! reference's o_0..o_{n-1}, `num_used` is n, and its slice holds o_k at the scan's position of o_k:
//! the slice is the sequence the reference's scan copies for the cell. The PV bits, `pv_array_bts`
//! and `longest_chain` are `finalize`'s functions of the counts, computed by its statements.
//!
//! `LOSAT_X_CLUTPAIRS` (a strict-set switch, on unless it is `0`) takes the cells of the two
//! innermost positions of the recursion from a row per pair of query letters (`XPairRow`,
//! `XWords::neighbors_paired`): the same cells of each query offset, those of the last two positions
//! in the row's order instead of the loops'. A cell gets at most one append per query offset (its
//! cells are distinct, `XWords::for_each_word_in`), so each cell still gets o_0..o_{n-1}, and every
//! test of the passes (the count chunks, the refusal, the bucket and piece cuts) runs at the end of
//! a query offset: the order within an offset changes no value.
//!
//! The refusal. The reference tests NCBI's tables (`tables`, ranges of the segments) while it adds the
//! words: after each query offset of a table, while the table's appends have taken no more cells of the
//! table of all batches than NCBI's banks hold, the table passes (a batch's n offsets of a cell take
//! past that, once, it refuses when no words come before the table (the table so far is the batch's
//! own, and its cells only grow) or when the batch's own table fills the banks
//! (`compressed_table_fills_banks`, built from all the batch's segments). So it refuses exactly when
//! some table's own need, the sum of cells(n) over its cells, passes the capacity: such a table's cells
//! of the table of all batches pass the capacity too, and the test then refuses; and a refusal needs
//! that need past the capacity. The count pass counts the same cells of the table of all batches (one
//! when an append finds the count c = 5, or c >= 6 with (c - 3) & 3 = 0: the two calls of
//! `compressed_list_get_new_cell` in `add_word_hit`) and applies the same rule after each chunk of
//! appends and at the end of each table (the cells only grow, so a table that passes the bound passes
//! it at its end); `x_table_fills_banks` counts the batch's own appends per cell with the same rule and
//! is true exactly when they take more cells than the capacity, as the reference's own table is marked
//! exhausted exactly when a request finds the capacity taken. So the build refuses exactly when the
//! reference does, and the caller gives the same message.
//!
//! The scans. `scan_subject` is the reference's `scan_subject` (NCBI's `s_BlastCompressedAaScanSubject`)
//! statement for statement, with a cell's hits copied from its slice: the same `numhits` (the slice's
//! length, `num_used`), the same order. `for_each_cell` is the scan of the fused loop of
//! `blastp/x_seed_fused_compressed.rs` (`x_scan_cells`, exactness argument in that file's header) with
//! the slice in place of `x_cell_hits`.
//!
//! The bucket fill (`LOSAT_X_BPFFILL`, `x_try_build_flat_by_buckets`, a strict-set switch: on
//! unless `=0`) makes the same table with the count pass's random updates of the counts and the fill
//! pass's random writes replaced by passes within buckets of 2^10 cells:
//! - appends: the count pass's appends, in its order, recorded as it records them (the cell of each
//!   append and the appends of each query offset), not counted. An append takes at most one overflow
//!   cell (one call of `compressed_list_get_new_cell` at most), so a table whose appends do not pass
//!   the capacity takes no more cells than the capacity, and the count pass's test of it passes; the
//!   bucket fill declines a table whose appends pass it (the count-and-fill builder then runs, and
//!   refuses or builds as above). The bound `x_max_appends` is tested after each query offset, as
//!   the count pass tests it, so `TooLarge` comes at the same offset; no refusal can come before it
//!   (every table so far stayed within the capacity).
//! - partition: the appends sorted by bucket (the cell's index over 2^10), stably (a counting sort),
//!   each packed as `query_offset << 10 | index & (2^10 - 1)` (the bucket fill declines a query of
//!   2^22 letters or more, whose offsets would not fit). A bucket's cells are consecutive, so its
//!   appends take the slots of `hits` that its cells' slices will take, and they are written there.
//! - per bucket, in order: each cell's count n (the number of its appends), the starts as the counts'
//!   prefix sums from the bucket's first slot (the count-and-fill builder's starts, which are the
//!   prefix sums in index order), the PV bits (shift `PV_ARRAY_BTS` whatever the counts for word
//!   size 5, as in `x_clut_par.rs`; the bucket fill declines the other word sizes, which
//!   `-task blastp-fast` does not take) and the count of non-empty cells and longest chain; then the
//!   fill from a copy of the bucket's
//!   packed appends (the record of the appends is not read again) in their order, each at its cell's
//!   cursor, so a cell's slice holds o_0..o_{n-1} as the fill pass leaves it; then o_k goes to slot
//!   `x_emission_position(k, n)` in each cell of more than 5 offsets, the fill pass's statements.
//!
//! The table is the count-and-fill builder's slot for slot; only the order of the memory accesses
//! differs.
//!
//! BLASTX's tables hold `u16` query offsets (`XBxCompressedLookup::Flat`; its query batches of more
//! than 2^16 letters keep `i32` ones), so an offset above a cell's index in its bucket does not fit
//! a slot of `hits`. Its bucket fill (`x_bx_buckets_on`) keeps no record and no partition: each
//! append goes, as the count pass makes it, to the chunks of its bucket (`XBucketChunks`) as
//! `query_offset << 16 | index & (2^10 - 1)`, and each bucket is filled from its chunks in the same
//! order and by the same statements.
//!
//! With more than one thread the search's pool builds the same table (`x_clut_par.rs`, exactness
//! argument in its header); at one thread (and in builds without threads) the serial builder runs:
//! the bucket fill, or the count-and-fill builder where it declines.
//!
//! A strict-set switch (`utils/xswitch.rs`): on unless `LOSAT_X_BPCLUT=0`.
//!
//! `LOSAT_X_BPLKPFLAT` (also a strict-set switch, and it acts only where the flat table is built) changes only how the pool's build is
//! scheduled and where its vectors come from (`x_clut_par.rs`, exactness argument in its header): one
//! latch per pass, the table's `starts` and `hits` taken from the previous table of a block search
//! (`XFlatSpare`), and tasks of another computation (the contexts' Karlin blocks of
//! `blastp/x_query_blocks.rs`) run by the appends pass after its chunks (`XExtraTasks`).
//! `LOSAT_X_BPLKPAHEAD` (a strict-set switch since 2026-10-17) runs
//! the same build for the next query block ahead of that block, on vectors that the block search's
//! thread allocated (`XFlatArgs::x_ahead`, `XFlatSpare::x_make_room`).

use super::*;
use crate::algorithm::tblastx::scan::OffsetPair;
use crate::utils::threading::SearchPool;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, PoisonError};

// EXPERIMENT (LOSAT_X_BPCLUT at more than one thread): the same table built by the search's pool
// (exactness argument in the file's header).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_clut_par.rs"]
mod x_clut_par;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) use x_clut_par::XFlatAppends;

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPCLUT): a strict-set switch, on unless it is `0`: the flat table
/// (`true`) or the reference builder (`false`).
pub(crate) fn x_bpclut_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPCLUT"))
}

/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the mode the variable selects.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum XLkpFlat {
    #[default]
    Off,
    On,
}

// No NCBI counterpart: reads the LOSAT_X_BPLKPFLAT switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPLKPFLAT): a strict-set switch (on unless it is `0`).
pub(crate) fn x_bplkpflat_mode() -> XLkpFlat {
    use std::sync::OnceLock;
    static MODE: OnceLock<XLkpFlat> = OnceLock::new();
    *MODE.get_or_init(|| {
        if crate::utils::xswitch::on("LOSAT_X_BPLKPFLAT") {
            XLkpFlat::On
        } else {
            XLkpFlat::Off
        }
    })
}

// No NCBI counterpart: counters of the `[X_STATS] BPLKPFLAT` line at exit (LOSAT_X_STATS); they do
// not change any value NCBI computes.
static X_LKPFLAT_TABLES: AtomicU64 = AtomicU64::new(0);
static X_LKPFLAT_RECYCLED: AtomicU64 = AtomicU64::new(0);
static X_LKPFLAT_RECYCLED_HITS: AtomicU64 = AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BPLKPFLAT): with `LOSAT_X_STATS`, the line at exit (printed by `main`) when
/// the pool built a table: `tables` the flat tables the switch's build made; `recycled` and
/// `recycled_hits` the builds whose `starts`, and whose `hits`, are the buffer of an earlier table
/// (`hits` only where it had the room).
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn x_bplkpflat_print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let load = |counter: &AtomicU64| counter.load(Ordering::Relaxed);
        let tables = load(&X_LKPFLAT_TABLES);
        let recycled = load(&X_LKPFLAT_RECYCLED);
        let recycled_hits = load(&X_LKPFLAT_RECYCLED_HITS);
        if tables != 0 || recycled != 0 || recycled_hits != 0 {
            eprintln!(
                "[X_STATS] BPLKPFLAT tables={tables} recycled={recycled} recycled_hits={recycled_hits}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1318-1320,1359-1373
// NCBI allocates a zeroed table for each search and frees it with the search; a block search makes
// a table per block. Buffer reuse only: the parallel build writes every entry of the vectors it keeps
// before any entry is read (`x_clut_par.rs`), so the values left in them do not matter.
/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the `starts` and `hits` of a flat table that its search no longer
/// reads, kept by the block search (`blastp/x_query_blocks.rs` `XBlockSubjects`) for the next
/// parallel build: the table's vectors are allocated by the caller of the builds, once, and then
/// recycled instead of being zeroed again for each block. The second field (LOSAT_X_BPLKPAHEAD) holds
/// a PV array's vector for a build ahead of its block (`x_make_room`); it is empty otherwise.
#[derive(Default)]
pub(crate) struct XFlatSpare<H = i32>(Mutex<(Vec<u32>, Vec<H>)>, Mutex<Vec<u32>>);

impl XFlatSpare {
    // No NCBI counterpart: buffer reuse only (above).
    /// Keeps `table`'s `starts` and `hits` (in place of any vectors kept before); the rest of the
    /// table is dropped.
    pub(crate) fn give(&self, table: XFlatCompressedLookup) {
        *self.0.lock().unwrap_or_else(PoisonError::into_inner) = (table.starts, table.hits);
    }

    // No NCBI counterpart: buffer reuse only (above).
    /// Item SC2d (LOSAT_X_BPLKPAHEAD): moves the `starts` and `hits` kept here, if any, to `other`
    /// (in place of any it keeps); this spare then keeps none. The second spare gives the vectors of
    /// a block's table to the block spare where the build ahead of the next block made no table, so
    /// that the next block's own build takes them (`blastp/blast_engine/x_lookup_ahead.rs`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(crate) fn x_give_to(&self, other: &XFlatSpare) {
        // Each lock is held only around a move.
        let kept = std::mem::take(&mut *self.0.lock().unwrap_or_else(PoisonError::into_inner));
        if kept.0.capacity() != 0 || kept.1.capacity() != 0 {
            *other.0.lock().unwrap_or_else(PoisonError::into_inner) = kept;
        }
    }

    // No NCBI counterpart: buffer reuse, and the thread that allocates the table's vectors
    // (LOSAT_X_BPLKPAHEAD); it changes no value NCBI computes (the build writes every entry it keeps
    // before any is read, `x_clut_par.rs`).
    /// EXPERIMENT (LOSAT_X_BPLKPAHEAD): on the block search's thread, before the build of the next
    /// block's table ahead of that block (`blastp/blast_engine/x_lookup_ahead.rs`) takes this spare,
    /// which then allocates none of the table's vectors (`XFlatArgs::x_ahead`): keeps here vectors
    /// with the room for a table of `backbone_size` cells and `hits` hits and for its PV array at any
    /// PV shift. The vectors are this spare's own, else those `block` keeps (the block spare, which
    /// the block before gave back); `starts` and the PV array short of their room are replaced here
    /// by empty ones of `backbone_size + 1` and `(backbone_size >> PV_ARRAY_BTS) + 1` entries of
    /// capacity, and `hits` short of its room is grown in place to `hits_grown` (at least `hits`)
    /// entries of capacity (item SC2d; its entries are kept). Nothing is written into them. Returns
    /// the number of vectors allocated or grown.
    pub(crate) fn x_make_room(
        &self,
        block: &XFlatSpare,
        backbone_size: usize,
        hits: usize,
        hits_grown: usize,
    ) -> usize {
        // Each lock is held only around a move.
        let take = |spare: &XFlatSpare| {
            std::mem::take(&mut *spare.0.lock().unwrap_or_else(PoisonError::into_inner))
        };
        let mut kept = take(self);
        if kept.0.capacity() == 0 && kept.1.capacity() == 0 {
            kept = take(block);
        }
        let mut pv = std::mem::take(&mut *self.1.lock().unwrap_or_else(PoisonError::into_inner));
        let pv_words = (backbone_size >> PV_ARRAY_BTS) + 1;
        let mut allocated = 0;
        if kept.0.capacity() < backbone_size + 1 {
            kept.0 = Vec::with_capacity(backbone_size + 1);
            allocated += 1;
        }
        if kept.1.capacity() < hits {
            // Item SC2d: grown in place of a new vector, so that the pages the vector's last table
            // touched stay with it (the room past them is not resident).
            let grown = hits_grown.max(hits);
            kept.1.reserve_exact(grown - kept.1.len().min(grown));
            allocated += 1;
        }
        if pv.capacity() < pv_words {
            pv = Vec::with_capacity(pv_words);
            allocated += 1;
        }
        *self.0.lock().unwrap_or_else(PoisonError::into_inner) = kept;
        *self.1.lock().unwrap_or_else(PoisonError::into_inner) = pv;
        allocated
    }
}

// No NCBI counterpart: buffer reuse between the builds of a block search (memory placement only); it
// changes no value NCBI computes.
/// Item SC2d-2: segments of packed appends (`x_clut_par.rs` `XSegmentOut`) kept by the block search
/// (`blastp/x_query_blocks.rs` `XBlockSubjects`) from one build's stripes pass to the next build's
/// appends pass, at most as many as the build's threads, each cut to the entries that build wrote (the
/// pages past them returned), so that the next build writes their pages, which the last build
/// touched, instead of faulting the pages of fresh mappings. A kept segment's values are not read: a
/// piece writes every entry of its range before any is read. What the pool holds is at most the last
/// build's appends (as many as the hits of the table that build made).
#[derive(Default)]
pub(crate) struct XSegmentPool(Mutex<Vec<Vec<u32>>>);

impl XSegmentPool {
    // No NCBI counterpart: buffer reuse only (above).
    /// A kept segment of at least `n` entries, if any (the pool then keeps it no more).
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    fn take(&self, n: usize) -> Option<Vec<u32>> {
        // The lock is held only around a move.
        let mut kept = self.0.lock().unwrap_or_else(PoisonError::into_inner);
        let at = kept.iter().position(|segment| segment.len() >= n)?;
        Some(kept.swap_remove(at))
    }

    // No NCBI counterpart: buffer reuse only (above).
    /// Keeps at most `at_most` of `segments` (a build's, each cut to its entries) in place of what the
    /// pool holds: segments of earlier builds that this build did not take (none long enough for its
    /// pieces, or more than it needed) are dropped, so the pool holds at most the last build's.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    fn keep(&self, segments: Vec<Vec<u32>>, at_most: usize) {
        let new: Vec<Vec<u32>> = segments
            .into_iter()
            .filter(|segment| !segment.is_empty())
            .take(at_most)
            .map(|mut segment| {
                // The capacity past the entries written is returned (their pages too).
                segment.shrink_to_fit();
                segment
            })
            .collect();
        self.replace(new);
    }

    // No NCBI counterpart: buffer reuse only (item SC9i of line NP).
    /// Keeps at most `at_most` of `segments` (a wave's) for the next wave of the same build, in place
    /// of what the pool holds, each whole: at the length it was made with, so that the next wave
    /// writes the pages this one touched and no mapping is cut, freed or made between the waves of a
    /// build. A cut or freed mapping stops every thread of the process for its TLB: cut after each
    /// wave (`keep`), BLASTX `-word_size 6` with 800 cDNAs at 32 threads made 16,318 `mremap` and
    /// 13,784 `munmap` calls instead of dev's 71 and 2,412. The last wave's segments are cut
    /// (`keep`).
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    fn keep_whole(&self, segments: Vec<Vec<u32>>, at_most: usize) {
        let new: Vec<Vec<u32>> = segments
            .into_iter()
            .filter(|segment| segment.capacity() != 0)
            .take(at_most)
            .map(|mut segment| {
                // SAFETY: a segment is made as `vec![0u32; n]` (capacity n, every entry written) and
                // cut since only by `truncate` (`XSegmentOut::finish`), which leaves the entries past
                // its length written, or by `keep` (`shrink_to_fit`: capacity = length).
                unsafe { segment.set_len(segment.capacity()) };
                segment
            })
            .collect();
        self.replace(new);
    }

    /// Puts `new` in place of what the pool holds.
    fn replace(&self, new: Vec<Vec<u32>>) {
        // The lock is held only around a move; the old segments are dropped after it.
        let old = std::mem::replace(
            &mut *self.0.lock().unwrap_or_else(PoisonError::into_inner),
            new,
        );
        drop(old);
    }
}

/// EXPERIMENT (LOSAT_X_BPLKPFLAT): tasks of another computation that the parallel build's appends
/// pass runs after its chunks (`x_clut_par.rs`), each once; the caller of the build runs those the
/// build did not run (the serial builder runs none). A task does not read the table.
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    allow(dead_code)
)]
pub(crate) trait XExtraTasks: Sync {
    /// The number of tasks.
    fn count(&self) -> usize;
    /// Runs task `index` (below `count`).
    fn run(&self, index: usize);
}

/// EXPERIMENT (LOSAT_X_BPLKPFLAT): what the pool's build of a flat table gets besides its inputs:
/// the switch's mode (`Off`: the `par_iter` passes on fresh vectors; `spare` and `extra` are not
/// used), the vectors of an earlier table, and the extra tasks.
#[derive(Clone, Copy, Default)]
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    allow(dead_code)
)]
pub(crate) struct XFlatArgs<'a> {
    pub(crate) mode: XLkpFlat,
    /// The vectors of an earlier table (not read by the first part of a split build ahead,
    /// `x_bp_flat_ahead_appends`, whose second part takes the second spare's).
    pub(crate) spare: Option<&'a XFlatSpare>,
    pub(crate) extra: Option<&'a dyn XExtraTasks>,
    /// Item SC2d-2: the block search's kept segments of packed appends (`XSegmentPool`).
    pub(crate) segments: Option<&'a XSegmentPool>,
    /// LOSAT_X_BPLKPAHEAD: the build of a block's table ahead of that block, on a pool task
    /// (`blastp/blast_engine/x_lookup_ahead.rs`): it allocates none of the table's vectors, takes
    /// them all from `spare` (made ready by `XFlatSpare::x_make_room`), and makes no table where they
    /// lack the room, the appends pass the bound (`TooLarge`) or the switch's mode is `Off`; the
    /// block then builds its table itself.
    pub(crate) x_ahead: bool,
    /// The most appends the build makes, as its caller counted them (0 where it did not): the
    /// serial builders reserve their record of the appends once (`x_try_build_flat_reserving`).
    pub(crate) appends: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1189
const X_TARGET_PV_BYTES: usize = 262_144;

// No NCBI counterpart: the callers of `XFlatCompressedLookup::scan_subject` and
// `XFlatCompressedLookup::x_list_pairs_lean`, one copy each.
/// BLASTP's caller of `XFlatCompressedLookup::scan_subject` and `x_list_pairs_lean` (BLOSUM62's
/// scans, `blastp/x_seed_fused_clut.rs`; and the tests').
pub(crate) const X_BLASTP: u8 = 0;
/// TBLASTN's caller of `XFlatCompressedLookup::scan_subject` (`tblastn/search_seed.rs`).
pub(crate) const X_TBLASTN: u8 = 1;
/// BLASTX's caller of `XFlatCompressedLookup::scan_subject` (`blastx/seed.rs`, `x_seed_scan.rs`).
pub(crate) const X_BLASTX: u8 = 2;
/// BLASTP's scan of a `-subject` with lower-case letters (`-lcase_masking`), a copy of its own so
/// that `X_BLASTP` keeps its one caller.
pub(crate) const X_BLASTP_SOFT: u8 = 3;
/// BLASTP's scans with a matrix other than BLOSUM62 (item SC2b; `scan_subject` and
/// `x_list_pairs_lean`), copies of their own so that `X_BLASTP` keeps its one caller.
pub(crate) const X_BLASTP_MATRIX: u8 = 4;

/// EXPERIMENT (LOSAT_X_BPCLUT): the compressed lookup table with each cell's hits in one slice of
/// `hits` (`starts[index]..starts[index + 1]`), in NCBI's scan order. The other fields are those of
/// `BlastCompressedAaLookupTable` that the scans read.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:236-267
pub(crate) struct XFlatCompressedLookup<H: XQueryOffset = i32> {
    pub(crate) word_length: i32,
    pub(crate) compressed_alphabet_size: i32,
    pub(crate) reciprocal_alphabet_size: i32,
    pub(crate) longest_chain: i32,
    pub(crate) pv: Vec<u32>,
    pub(crate) pv_array_bts: i32,
    pub(crate) scaled_compress_table_u8: [i32; 256],
    /// `backbone_size + 1` entries: cell `index`'s hits are `hits[starts[index]..starts[index + 1]]`.
    pub(crate) starts: Vec<u32>,
    pub(crate) hits: Vec<H>,
}

// No NCBI counterpart: NCBI's cells hold Int4 query offsets.
/// The type of a flat table's query offsets: `i32` as NCBI's (BLASTP, TBLASTN), or `u16` for
/// BLASTX, whose query batches of at most 2 x 10002 letters have fewer than 2^16 letters in their
/// contexts (`X_U16_LETTERS`); the offsets are the same numbers.
pub(crate) trait XQueryOffset: Copy + Default + Ord + Send + Sync {
    fn from_offset(offset: u32) -> Self;
    fn offset(self) -> u32;
    /// `table` (of `i32` offsets) with offsets of this type, the same numbers.
    fn x_from_wide(table: XFlatCompressedLookup) -> XFlatCompressedLookup<Self>;
}

impl XQueryOffset for i32 {
    #[inline(always)]
    fn from_offset(offset: u32) -> Self {
        offset as i32
    }

    #[inline(always)]
    fn offset(self) -> u32 {
        self as u32
    }

    fn x_from_wide(table: XFlatCompressedLookup) -> XFlatCompressedLookup<Self> {
        table
    }
}

impl XQueryOffset for u16 {
    #[inline(always)]
    fn from_offset(offset: u32) -> Self {
        offset as u16
    }

    #[inline(always)]
    fn offset(self) -> u32 {
        u32::from(self)
    }

    fn x_from_wide(table: XFlatCompressedLookup) -> XFlatCompressedLookup<Self> {
        x_narrow_flat(table)
    }
}

/// The most letters a query of a `u16` flat table has (its offsets are below 2^16).
pub(crate) const X_U16_LETTERS: usize = 1 << 16;

/// What the flat builder made.
#[expect(
    clippy::large_enum_variant,
    reason = "boxing a variant changes the layout and the memory; that needs an RSS and time measurement"
)]
pub(crate) enum XFlatBuild<H: XQueryOffset = i32> {
    Table(XFlatCompressedLookup<H>),
    /// One of NCBI's tables fills its overflow banks (the reference's `Err(())`).
    Refused,
    /// More appends than the slices' `u32` offsets take: the reference builder runs instead.
    TooLarge,
}

impl XFlatBuild {
    /// The same outcome, a table's offsets of type `H` (`XQueryOffset::x_from_wide`).
    pub(crate) fn x_offsets<H: XQueryOffset>(self) -> XFlatBuild<H> {
        match self {
            XFlatBuild::Table(table) => XFlatBuild::Table(H::x_from_wide(table)),
            XFlatBuild::Refused => XFlatBuild::Refused,
            XFlatBuild::TooLarge => XFlatBuild::TooLarge,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 and
// c++/src/algo/blast/core/blast_aascan.c:348-396 (module header)
/// The position in NCBI's scan order of the offset appended with rank `k` (0-based) to a cell that
/// ends with `n` offsets (module header). The builders put a cell in that order with
/// `x_to_scan_order`; this is the statement of the order that the unit tests check it against.
#[cfg(test)]
pub(crate) fn x_emission_position(k: u32, n: u32) -> u32 {
    if n <= (COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1) || k < 3 {
        return k;
    }
    let cells = n >> 2;
    let g = (k - 3) >> 2;
    let s = (k - 3) & COMPRESSED_HITS_CELL_MASK as u32;
    if g + 1 == cells {
        3 + s
    } else {
        n - 4 * (g + 1) + s
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 and
// c++/src/algo/blast/core/blast_aascan.c:348-396 (module header, `x_emission_position`)
/// A cell's offsets, given in the order of their appends (o_0..o_{n-1}), put in NCBI's scan order in
/// place: o_k at `x_emission_position(k, n)`. That order keeps o_0..o_2 and takes the n >> 2 groups
/// after them (groups of four, the last one of (n & 3) + 1 offsets) last group first, each group in
/// its own order; so the slots after the third are reversed, and then the first group of the
/// reversed slots (the last group's) and each group of four after it are reversed back. No copy of
/// the cell is made (the unit test `scan_order_in_place_is_the_emission_positions` compares it with
/// `x_emission_position` for 0 to 300 offsets).
#[inline(always)]
pub(crate) fn x_to_scan_order<T>(cell: &mut [T]) {
    // `x_emission_position`'s `cells`: with one group (n <= 7, the cells of 5 or fewer offsets
    // among them) the order is the order of the appends.
    let groups = cell.len() >> 2;
    if groups < 2 {
        return;
    }
    let tail = &mut cell[3..];
    tail.reverse();
    let last_len = tail.len() - 4 * (groups - 1);
    let (last, full) = tail.split_at_mut(last_len);
    last.reverse();
    for group in full.chunks_exact_mut(4) {
        group.reverse();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:800-806,825-836
/// Whether the append to a cell that holds `c` offsets takes an overflow cell.
#[inline(always)]
fn x_takes_overflow_cell(c: u32) -> bool {
    c == COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1
        || (c > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1
            && (c - 3) & COMPRESSED_HITS_CELL_MASK as u32 == 0)
}

/// Where `XWords` puts the cells of the appends, in order.
trait XSink {
    fn one(&mut self, index: u32);
    /// The cells `partial + last[i]`, i in order.
    fn many(&mut self, partial: u32, last: &[u32]) {
        for &letter in last {
            self.one(partial + letter);
        }
    }
    /// The cells `partial + codes[i] * scale`, i in order (the pairs of `XPairRow`).
    fn many_coded(&mut self, partial: u32, scale: u32, codes: &[u8]) {
        for &code in codes {
            self.one(partial + u32::from(code) * scale);
        }
    }
}

// No NCBI counterpart: the count pass's record of the appends.
impl XSink for Vec<u32> {
    #[inline(always)]
    fn one(&mut self, index: u32) {
        self.push(index);
    }

    #[inline(always)]
    fn many(&mut self, partial: u32, last: &[u32]) {
        self.extend(last.iter().map(|&letter| partial + letter));
    }

    #[inline(always)]
    fn many_coded(&mut self, partial: u32, scale: u32, codes: &[u8]) {
        self.extend(codes.iter().map(|&code| partial + u32::from(code) * scale));
    }
}

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_CLUTPAIRS): a strict-set switch, on unless it is `0`: `XWords` takes the
/// cells of the two innermost positions it chooses from the pair table of their query letters
/// (`XPairTable`; `true`), or from the two loops of `add_word_hits` (`false`).
fn x_clutpairs_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_CLUTPAIRS"))
}

/// The pairs of a group of a row of `XPairTable` (`XPairGroup`), which `XPairRow::emit` compares at
/// a time (32 bytes of sums).
const X_PAIR_GROUP: usize = 16;
/// The largest compressed alphabet (word sizes 5 and 6), and the groups of its 15 x 15 pairs.
const X_PAIR_ALPHABET: usize = 15;
const X_PAIR_GROUPS: usize = (X_PAIR_ALPHABET * X_PAIR_ALPHABET).div_ceil(X_PAIR_GROUP);

// No NCBI counterpart: the pairs of letters of `add_word_hits`' two innermost loops (`XWords`).
/// The pairs of ranks g * 16 .. g * 16 + 15 of a row (`XPairRow`), as their sums and codes; the
/// slots past the alphabet's pairs hold `i16::MIN` and 0.
struct XPairGroup {
    sums: [i16; X_PAIR_GROUP],
    codes: [u8; X_PAIR_GROUP],
}

// No NCBI counterpart: the pairs of letters of `add_word_hits`' two innermost loops (`XWords`).
/// The row of `XPairTable` of query letters (a, b): the pairs (i, j) of ranks in their sorted rows
/// (`row_a = matrix_sorted[a]`, `row_b = matrix_sorted[b]`) in the order of `row_a[i] + row_b[j]`
/// descending, then of k = i * alphabet + j ascending (the stable sort by sum of the pairs in k
/// order), with the codes `matrix_sorted_char[a][i] + alphabet * matrix_sorted_char[b][j]` of their
/// two compressed letters (below alphabet^2), in groups of 16 (`XPairGroup`), made in order as they
/// are first read: the pairs a node reaches are a prefix of the row (`XWords::neighbors_paired`),
/// and most nodes read the first group only.
/// The row holds what its merge reads (`alphabet`, the sorted rows of a and b and their letters):
/// `make` and `emit` take no table, so nothing of the table is loaded on `emit`'s path (with the
/// table as an argument, the link-time optimizer moved `make`'s loads of its fields into `emit`'s
/// loop).
struct XPairRow {
    a: u8,
    b: u8,
    alphabet: u8,
    row_a: [i32; X_PAIR_ALPHABET],
    row_b: [i32; X_PAIR_ALPHABET],
    char_a: [u8; X_PAIR_ALPHABET],
    char_b: [u8; X_PAIR_ALPHABET],
    groups: [std::sync::OnceLock<XPairGroup>; X_PAIR_GROUPS],
    /// The merge of `XPairRow::make` after the groups made (`next`, see there; `made` groups).
    merge: Mutex<XPairMerge>,
}

/// The state of a row's merge between its groups: `next[i]` the rank j of row i's next pair,
/// `keys[i]` the key of that pair while row i is a candidate of the merge (`XPairRow::make`), else
/// `i32::MIN`, `made` the groups made.
struct XPairMerge {
    next: [u8; X_PAIR_ALPHABET],
    keys: [i32; X_PAIR_ALPHABET + 1],
    made: u8,
}

/// The merge's key of pair (i, j) of sum `sum` (`XPairRow::make`): the sum, then the smaller i, in
/// one i32 (the sums lie in i16, below 15 rows).
#[inline(always)]
fn x_pair_key(sum: i32, i: usize) -> i32 {
    (sum << 4) | (X_PAIR_ALPHABET - i) as i32
}

impl XPairRow {
    /// Group `g` (and the groups before it that are not made yet): the merge takes the row's pairs in
    /// order from the two sorted rows (the rows i of the pairs (i, j), j ascending), continuing from
    /// the last group made. `next[i]` is the rank j of row i's next pair. Row i's next pair is the
    /// next of the merge only when row i - 1 has gone past j (`next[i - 1] > next[i]`): else (i - 1,
    /// next[i - 1]) has a sum at least as large and comes first on a tie (the smaller k), so `next`
    /// never increases with i. These candidates hold their pair's key (`x_pair_key`: the sum, then
    /// the smaller i) in `keys`, the other rows `i32::MIN`, and the largest key is the next pair.
    /// Taking (i, j) changes only rows i (its next pair, or no candidate once row i - 1's next is
    /// j + 1 or its pairs end) and i + 1 (a candidate once row i passes its next rank j). One thread
    /// at a time (`merge`); a group is read only once it is set.
    #[cold]
    #[inline(never)]
    fn make(&self, g: usize) -> &XPairGroup {
        let mut merge = self.merge.lock().unwrap_or_else(PoisonError::into_inner);
        let merge = &mut *merge;
        // Copied once: the row's other fields share it with the merge's (a `Mutex`) and groups,
        // which the merge writes.
        let alphabet = usize::from(self.alphabet);
        let (a, b) = (usize::from(self.a), usize::from(self.b));
        let (row_a, row_b) = (self.row_a, self.row_b);
        let (row_a, row_b) = (&row_a[..alphabet], &row_b[..alphabet]);
        let (char_a, char_b) = (self.char_a, self.char_b);
        let made = merge.made;
        if made == 0 {
            // Every sum of the row lies between its first and last pair's (the rows are sorted):
            // `x_pair_table` checked them for the letters inside the compressed alphabet, the only
            // rows `XWords` reads; here for this row, so the keys and sums below fit.
            let (high, low) = (
                row_a[0] + row_b[0],
                row_a[alphabet - 1] + row_b[alphabet - 1],
            );
            assert!(
                low > i32::from(i16::MIN) && high <= i32::from(i16::MAX),
                "the pair sums {low}..={high} of letters ({a}, {b}) leave the pair table's i16"
            );
            merge.keys[0] = x_pair_key(high, 0);
        }
        while usize::from(merge.made) <= g {
            let first = usize::from(merge.made) * X_PAIR_GROUP;
            let end = (first + X_PAIR_GROUP).min(alphabet * alphabet);
            let mut group = XPairGroup {
                sums: [i16::MIN; X_PAIR_GROUP],
                codes: [0; X_PAIR_GROUP],
            };
            for slot in 0..end.saturating_sub(first) {
                let key = merge.keys.iter().fold(i32::MIN, |max, &key| max.max(key));
                let i = X_PAIR_ALPHABET - (key & 15) as usize;
                let j = usize::from(merge.next[i]);
                merge.next[i] = j as u8 + 1;
                merge.keys[i] =
                    if j + 1 < alphabet && (i == 0 || usize::from(merge.next[i - 1]) > j + 1) {
                        x_pair_key(row_a[i] + row_b[j + 1], i)
                    } else {
                        i32::MIN
                    };
                if i + 1 < alphabet && usize::from(merge.next[i + 1]) == j {
                    merge.keys[i + 1] = x_pair_key(row_a[i + 1] + row_b[j], i + 1);
                }
                group.sums[slot] = (key >> 4) as i16;
                group.codes[slot] = char_a[i] + alphabet as u8 * char_b[j];
            }
            let at = usize::from(merge.made);
            assert!(
                self.groups[at].set(group).is_ok(),
                "a pair group made twice"
            );
            merge.made += 1;
        }
        X_PAIR_GROUPS_MADE.fetch_add(u64::from(merge.made - made), Ordering::Relaxed);
        self.groups[g].get().expect("made above")
    }

    /// Gives `sink` the cells `partial + code * scale` of the row's pairs with a sum of at least
    /// `key`, a prefix of the row, in the row's order. The sums lie in (i16::MIN, i16::MAX]
    /// (`x_pair_table`), so a key above i16::MAX takes none, and one at or below i16::MIN + 1 takes
    /// all and no slot past them. A group at a time, without a branch inside a group; the last
    /// group of an alphabet's pairs is not full, so the loop ends there.
    #[inline(always)]
    fn emit(&self, key: i32, partial: u32, scale: u32, sink: &mut impl XSink) {
        if key > i32::from(i16::MAX) {
            return;
        }
        let key = key.max(i32::from(i16::MIN) + 1) as i16;
        for (g, group) in self.groups.iter().enumerate() {
            let group = match group.get() {
                Some(group) => group,
                None => self.make(g),
            };
            let reached: usize = group.sums.iter().map(|&sum| usize::from(sum >= key)).sum();
            sink.many_coded(partial, scale, &group.codes[..reached]);
            if reached < X_PAIR_GROUP {
                break;
            }
        }
    }
}

// No NCBI counterpart: the pairs of letters of `add_word_hits`' two innermost loops (`XWords`),
// for one compressed matrix; it changes no value NCBI computes.
/// The rows (`XPairRow`) of every pair of query letters for a compressed matrix (`alphabet`,
/// `matrix_sorted` and `matrix_sorted_char` of its `CompressedNeighborInfo`: the table's key),
/// each row made when first read (`XPairTable::row`, about 1,020 bytes on x86-64, its groups as
/// they are read). One table per matrix for the builds that run at once (`x_pair_table`: BLASTX's
/// query batches build on the pool's tasks side by side), freed when the last of them ends: a
/// table of each build kept the rows of each (blastx-fast, 800 queries, 16 threads: max RSS +5.2
/// MiB), and a table kept for the process the rows after its builds (UniRef50 at one thread:
/// +292 KiB). At most the 23 x 23 rows of the letters inside the 15-letter alphabet (0.52 MiB),
/// and 12.25 KiB for the table's `OnceLock`s.
struct XPairTable {
    alphabet: usize,
    matrix_sorted: [[i32; BLASTAA_SIZE]; BLASTAA_SIZE],
    matrix_sorted_char: [[u8; BLASTAA_SIZE]; BLASTAA_SIZE],
    rows: Vec<std::sync::OnceLock<Box<XPairRow>>>,
}

impl XPairTable {
    /// The row of query letters (a, b).
    #[inline(always)]
    fn row(&self, a: usize, b: usize) -> &XPairRow {
        self.rows[a * BLASTAA_SIZE + b].get_or_init(|| {
            X_PAIR_ROWS.fetch_add(1, Ordering::Relaxed);
            let alphabet = self.alphabet;
            let mut row = Box::new(XPairRow {
                a: a as u8,
                b: b as u8,
                alphabet: alphabet as u8,
                row_a: [0; X_PAIR_ALPHABET],
                row_b: [0; X_PAIR_ALPHABET],
                char_a: [0; X_PAIR_ALPHABET],
                char_b: [0; X_PAIR_ALPHABET],
                groups: [const { std::sync::OnceLock::new() }; X_PAIR_GROUPS],
                merge: Mutex::new(XPairMerge {
                    next: [0; X_PAIR_ALPHABET],
                    keys: [i32::MIN; X_PAIR_ALPHABET + 1],
                    made: 0,
                }),
            });
            row.row_a[..alphabet].copy_from_slice(&self.matrix_sorted[a][..alphabet]);
            row.row_b[..alphabet].copy_from_slice(&self.matrix_sorted[b][..alphabet]);
            row.char_a[..alphabet].copy_from_slice(&self.matrix_sorted_char[a][..alphabet]);
            row.char_b[..alphabet].copy_from_slice(&self.matrix_sorted_char[b][..alphabet]);
            row
        })
    }
}

// No NCBI counterpart: the pair tables of the builds that run (`x_pair_table`); it changes no
// value NCBI computes.
static X_PAIR_SHARED: Mutex<Vec<std::sync::Weak<XPairTable>>> = Mutex::new(Vec::new());

// No NCBI counterpart: counters of the `[X_STATS] CLUTPAIRS` line at exit (LOSAT_X_STATS); they
// do not change any value NCBI computes.
static X_PAIR_TABLES: AtomicU64 = AtomicU64::new(0);
static X_PAIR_ROWS: AtomicU64 = AtomicU64::new(0);
static X_PAIR_GROUPS_MADE: AtomicU64 = AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_CLUTPAIRS): with `LOSAT_X_STATS`, the line at exit (printed by `main`) when
/// `XWords` read a pair table: `tables` the pair tables made (one per matrix for the builds that
/// run at once, `x_pair_table`), `rows` their rows made (one per pair of query letters read), `groups` the rows'
/// groups of 16 pairs made.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn x_clutpairs_print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let tables = X_PAIR_TABLES.load(Ordering::Relaxed);
        let rows = X_PAIR_ROWS.load(Ordering::Relaxed);
        let groups = X_PAIR_GROUPS_MADE.load(Ordering::Relaxed);
        if tables != 0 || rows != 0 || groups != 0 {
            eprintln!("[X_STATS] CLUTPAIRS tables={tables} rows={rows} groups={groups}");
        }
    }
}

/// The pair table of the compressed matrix of `info` for the alphabet of `header` (`XWords::new`):
/// a running build's table of the same key (`X_PAIR_SHARED`), else a new one. Every sum of two
/// letters that `XWords` reads (letters inside the compressed alphabet) lies in (i16::MIN,
/// i16::MAX] for every matrix of NCBI's (test `pair_tables_hold_every_matrix`): checked here.
fn x_pair_table(header: &XHeader, info: &CompressedNeighborInfo) -> Arc<XPairTable> {
    let alphabet = header.alphabet;
    let mut shared = X_PAIR_SHARED.lock().unwrap_or_else(PoisonError::into_inner);
    shared.retain(|table| table.strong_count() > 0);
    let same = |table: &XPairTable| {
        table.alphabet == alphabet
            && table.matrix_sorted == info.matrix_sorted
            && table.matrix_sorted_char == info.matrix_sorted_char
    };
    if let Some(table) = shared
        .iter()
        .filter_map(std::sync::Weak::upgrade)
        .find(|table| same(table))
    {
        return table;
    }
    // A row's groups, and a code below alphabet^2 (one byte).
    assert!(
        alphabet <= X_PAIR_ALPHABET,
        "alphabet {alphabet} too large for the pair table"
    );
    let (mut high, mut low) = (i64::MIN, i64::MAX);
    for (c, &compressed) in header.compress_table.iter().enumerate() {
        if usize::from(compressed) < alphabet {
            high = high.max(i64::from(info.matrix_sorted[c][0]));
            low = low.min(i64::from(info.matrix_sorted[c][alphabet - 1]));
        }
    }
    assert!(
        2 * high <= i64::from(i16::MAX) && 2 * low > i64::from(i16::MIN),
        "the compressed matrix's pair sums ({}..={}) leave the pair table's i16",
        2 * low,
        2 * high
    );
    X_PAIR_TABLES.fetch_add(1, Ordering::Relaxed);
    let table = Arc::new(XPairTable {
        alphabet,
        matrix_sorted: info.matrix_sorted,
        matrix_sorted_char: info.matrix_sorted_char,
        rows: (0..BLASTAA_SIZE * BLASTAA_SIZE)
            .map(|_| std::sync::OnceLock::new())
            .collect(),
    });
    shared.push(Arc::downgrade(&table));
    table
}

/// No NCBI counterpart: a query word as `XWords` reads it: a slice of the table's word length
/// (`[u8]`, read at run time) or an array of a word length fixed by the build's instance
/// (`[u8; W]`, BLASTX's `x_bx_buckets_on`, item SC9 of line NP), whose length and checks the
/// compiler then folds, as it did for `-word_size 5` when that was BLASTX's only compressed one.
trait XQueryWord {
    /// The word length of the type, `None` for a slice.
    const LEN: Option<usize>;
    /// The letters.
    fn letters(&self) -> &[u8];
}

impl XQueryWord for [u8] {
    const LEN: Option<usize> = None;
    #[inline(always)]
    fn letters(&self) -> &[u8] {
        self
    }
}

impl<const W: usize> XQueryWord for [u8; W] {
    const LEN: Option<usize> = Some(W);
    #[inline(always)]
    fn letters(&self) -> &[u8] {
        self
    }
}

/// The words of one query offset, as `add_word_hits` adds them (module header).
struct XWords<'a> {
    word_length: usize,
    alphabet: usize,
    threshold: i32,
    compress_table: &'a [u8; BLASTAA_SIZE],
    scaled_compress_table: &'a [i32; BLASTAA_SIZE],
    matrix: &'a CompressedScoreMatrix,
    info: &'a CompressedNeighborInfo,
    /// `A^i` for the letter at position i of a word.
    power: [u32; 8],
    /// `matrix_sorted_char[c][i] * A^(word_length - 1)`: the last letter's term of a neighbor's index.
    last: [[u32; BLASTAA_SIZE]; BLASTAA_SIZE],
    /// `matrix_sorted_char[c][i] * A^(word_length - 2)`: the term of the letter before the last, the
    /// one that `for_each_word_in` chooses last (the parallel builder's, `x_clut_par.rs`).
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    below_last: [[u32; BLASTAA_SIZE]; BLASTAA_SIZE],
    /// EXPERIMENT (LOSAT_X_CLUTPAIRS): the pair table of the matrix (`XPairTable`), from which
    /// `neighbors_paired` and `neighbors_below_last_paired` take the cells of the two innermost
    /// positions; `None` with the switch off (`neighbors`, `neighbors_below_last`).
    /// A pointer: `XWords`' other fields keep the places they have without it.
    pairs: Option<Arc<XPairTable>>,
}

// No NCBI counterpart: the order of the build's work (item SC9i of line NP).
/// What `XWords::word_start` reads of a query offset's word: the offset, `row_max` of the positions
/// before the last (`bound`: the sum that the last letter's loop adds a letter's score to), the
/// word's own cell (`NO_CELL` when `add_word_hits` adds none: the word reaches the threshold or is
/// skipped), and the query's residue at the last position.
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    allow(dead_code)
)]
#[derive(Clone, Copy, Debug)]
struct XWordStart {
    offset: i32,
    bound: i32,
    own: u32,
    last: u8,
}

impl XWordStart {
    /// `own` when the word has no cell of its own (a backbone has fewer than 2^32 cells).
    const NO_CELL: u32 = u32::MAX;
}

impl<'a> XWords<'a> {
    /// No NCBI counterpart: the word length of a word of type `Q` (`XQueryWord`): the type's, or
    /// the table's.
    #[inline(always)]
    fn length_of<Q: XQueryWord + ?Sized>(&self) -> usize {
        debug_assert!(Q::LEN.is_none_or(|w| w == self.word_length));
        Q::LEN.unwrap_or(self.word_length)
    }

    fn new(
        header: &'a XHeader,
        threshold: i32,
        matrix: &'a CompressedScoreMatrix,
        info: &'a CompressedNeighborInfo,
    ) -> Self {
        Self::with_pairs(header, threshold, matrix, info, x_clutpairs_on())
    }

    /// `new`, with the pair table (`pairs`, LOSAT_X_CLUTPAIRS) or without it.
    fn with_pairs(
        header: &'a XHeader,
        threshold: i32,
        matrix: &'a CompressedScoreMatrix,
        info: &'a CompressedNeighborInfo,
        pairs: bool,
    ) -> Self {
        let mut power = [0u32; 8];
        let mut scale = 1u32;
        for p in power.iter_mut().take(header.word_length) {
            *p = scale;
            scale *= header.alphabet as u32;
        }
        let mut last = [[0u32; BLASTAA_SIZE]; BLASTAA_SIZE];
        let mut below_last = [[0u32; BLASTAA_SIZE]; BLASTAA_SIZE];
        for (c, (row, below_row)) in last.iter_mut().zip(below_last.iter_mut()).enumerate() {
            for (i, (term, below_term)) in row
                .iter_mut()
                .zip(below_row.iter_mut())
                .enumerate()
                .take(header.alphabet)
            {
                *term = u32::from(info.matrix_sorted_char[c][i]) * power[header.word_length - 1];
                *below_term =
                    u32::from(info.matrix_sorted_char[c][i]) * power[header.word_length - 2];
            }
            // `neighbors` takes the last letters as the prefix of the sorted row whose scores reach
            // the threshold; the recursion's loop stops at the first that does not, the same prefix
            // when the row does not increase (the sort of `CompressedNeighborInfo::new`).
            assert!(
                info.matrix_sorted[c][..header.alphabet]
                    .windows(2)
                    .all(|pair| pair[0] >= pair[1]),
                "the sorted compressed matrix rows do not increase"
            );
        }
        Self {
            word_length: header.word_length,
            alphabet: header.alphabet,
            threshold,
            compress_table: &header.compress_table,
            scaled_compress_table: &header.scaled_compress_table,
            matrix,
            info,
            power,
            last,
            below_last,
            pairs: pairs.then(|| x_pair_table(header, info)),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1075-1121
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:919-929
    /// Gives `sink` every cell that `add_word_hits` appends the offset of `word` to, in its order
    /// (with `pairs`, the cells of the last two positions in their pair row's order: module
    /// header).
    #[inline]
    fn for_each_word<Q: XQueryWord + ?Sized>(&self, query_word: &Q, sink: &mut impl XSink) {
        let word = query_word.letters();
        let mut self_score = 0i32;
        for &residue in word {
            let compressed = self.compress_table[residue as usize] as usize;
            if compressed >= self.alphabet {
                return;
            }
            self_score += self.matrix[residue as usize][compressed];
        }
        if self.threshold == 0 || self_score < self.threshold {
            let (index, skip) = compute_compressed_index(
                self.length_of::<Q>(),
                word,
                self.alphabet,
                self.scaled_compress_table,
            );
            if skip == 0 {
                sink.one(index as u32);
            }
        }
        if self.threshold == 0 {
            return;
        }
        let mut score = self.info.row_max[word[0] as usize];
        for &residue in &word[1..] {
            score += self.info.row_max[residue as usize];
        }
        match &self.pairs {
            Some(pairs) => {
                // The row of the query letters at the last two positions, read once for the word.
                let last = self.length_of::<Q>() - 1;
                let row = pairs.row(word[last - 1] as usize, word[last] as usize);
                self.neighbors_paired(row, query_word, score, 0, 0, sink);
            }
            None => self.neighbors(query_word, score, 0, 0, sink),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:894-896
    // `partial` is the index of the letters chosen so far, sum c[i] * A^i over i < current_pos (the
    // terms of `compute_encoded_compressed_index`). At the last position the loop's letters are the
    // prefix of the sorted row up to the first score below the threshold (`partition_point`; the row
    // does not increase, checked in `new`).
    fn neighbors<Q: XQueryWord + ?Sized>(
        &self,
        word: &Q,
        mut score: i32,
        current_pos: usize,
        partial: u32,
        sink: &mut impl XSink,
    ) {
        let query_word = word.letters();
        let curr_query_char = query_word[current_pos] as usize;
        score -= self.info.row_max[curr_query_char];
        let row_sorted = &self.info.matrix_sorted[curr_query_char][..self.alphabet];
        let threshold = self.threshold;

        if current_pos == self.length_of::<Q>() - 1 {
            let count =
                row_sorted.partition_point(|&letter_score| score + letter_score >= threshold);
            sink.many(partial, &self.last[curr_query_char][..count]);
            return;
        }

        let char_sorted = &self.info.matrix_sorted_char[curr_query_char];
        let scale = self.power[current_pos];
        if current_pos + 2 == self.length_of::<Q>() {
            // The call for the last position, inline: the same statements as the branch above.
            let last_char = query_word[current_pos + 1] as usize;
            let last_row = &self.info.matrix_sorted[last_char][..self.alphabet];
            let last = &self.last[last_char];
            let last_max = self.info.row_max[last_char];
            for (i, &letter_score) in row_sorted.iter().enumerate() {
                if score + letter_score < threshold {
                    break;
                }
                let next_score = score + letter_score - last_max;
                let count =
                    last_row.partition_point(|&last_score| next_score + last_score >= threshold);
                sink.many(partial + u32::from(char_sorted[i]) * scale, &last[..count]);
            }
            return;
        }
        for (i, &letter_score) in row_sorted.iter().enumerate() {
            if score + letter_score < threshold {
                break;
            }
            self.neighbors(
                word,
                score + letter_score,
                current_pos + 1,
                partial + u32::from(char_sorted[i]) * scale,
                sink,
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1075-1121 and :1019-1053
    // (`for_each_word`, `neighbors`)
    /// Item SC9i of line NP: the cells of `for_each_word` whose last letter (the index's most
    /// significant digit, of weight `A^(word_length - 1)`) is in `letters`, each once, in another
    /// order: the last letter first, only among `letters`, then the others as `neighbors` chooses
    /// them (from position 0 to the one before the last, whose loop is the innermost, so that
    /// consecutive cells lie `A^(word_length - 2)` apart, in different stripes of the parallel
    /// builder as `neighbors`' consecutive cells are). Each is given less the first cell of
    /// `letters`, `letters.start * A^(word_length - 1)` (the index in the parallel builder's wave,
    /// `x_clut_par.rs`).
    ///
    /// The same cells. `neighbors` gives every compressed word `c` with `sum_i matrix[w[i]][c_i] >=
    /// threshold` once: each position's loop runs the letters of a row sorted by score down to the
    /// first whose score with the best of the positions after it (`row_max`) misses the threshold, so
    /// the bound never drops a word that reaches it, and the words are the leaves of a tree of
    /// distinct letter choices. Choosing the positions in another order gives the same set (the sum
    /// and the bound do not depend on the order), and so does testing the last position's letters
    /// one by one, in any order, against the same bound (the sorted loop's break only skips letters
    /// that fail it). The word itself is appended (`threshold == 0 || self_score < threshold`) at
    /// its index, which no neighbor has (a neighbor's sum reaches the threshold, the word's is
    /// `self_score`).
    ///
    /// The order does not matter to the table. For one query offset the cells are distinct (above),
    /// so a cell's appends are one per offset, in the order of the offsets, whatever the order of
    /// the cells within an offset; the overflow cells an append takes depend only on its cell's count
    /// before it (`x_takes_overflow_cell`), its rank among the cell's appends.
    #[cfg(test)]
    fn for_each_word_in<Q: XQueryWord + ?Sized>(
        &self,
        query_word: &Q,
        letters: &std::ops::Range<u32>,
        sink: &mut impl XSink,
    ) {
        if let Some(start) = self.word_start(query_word, 0) {
            self.for_each_word_from(&start, query_word, letters, sink);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1075-1121 (as
    // `for_each_word`)
    /// What `for_each_word_in` reads of the word at `offset` before its last letter's loop, the same
    /// for every range of last letters (item SC9i of line NP: the parallel builder's first wave keeps
    /// it for the later waves, `x_clut_par.rs`); `None` when the word has a letter outside the
    /// alphabet (`add_word_hits` adds nothing for it).
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    #[inline]
    fn word_start<Q: XQueryWord + ?Sized>(
        &self,
        query_word: &Q,
        offset: i32,
    ) -> Option<XWordStart> {
        let word = query_word.letters();
        let last = self.length_of::<Q>() - 1;
        let mut self_score = 0i32;
        for &residue in word {
            let compressed = self.compress_table[residue as usize] as usize;
            if compressed >= self.alphabet {
                return None;
            }
            self_score += self.matrix[residue as usize][compressed];
        }
        let mut own = XWordStart::NO_CELL;
        if self.threshold == 0 || self_score < self.threshold {
            let (index, skip) = compute_compressed_index(
                self.length_of::<Q>(),
                word,
                self.alphabet,
                self.scaled_compress_table,
            );
            if skip == 0 {
                own = index as u32;
            }
        }
        let mut bound = 0i32;
        if self.threshold != 0 {
            for &residue in &word[..last] {
                bound += self.info.row_max[residue as usize];
            }
        }
        Some(XWordStart {
            offset,
            bound,
            own,
            last: word[last],
        })
    }

    /// `for_each_word_in` from the word's `start` (`word_start`): the word's own cell when its last
    /// letter (the index's most significant digit) is in `letters`, then the last position's loop
    /// over `letters` (`neighbors`' statements) and the positions below it.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    #[inline]
    fn for_each_word_from<Q: XQueryWord + ?Sized>(
        &self,
        start: &XWordStart,
        query_word: &Q,
        letters: &std::ops::Range<u32>,
        sink: &mut impl XSink,
    ) {
        let last = self.length_of::<Q>() - 1;
        let last_char = start.last as usize;
        if start.own != XWordStart::NO_CELL
            && letters.contains(&u32::from(self.compress_table[last_char]))
        {
            debug_assert_eq!(
                start.own / self.power[last],
                u32::from(self.compress_table[last_char])
            );
            sink.one(start.own - letters.start * self.power[last]);
        }
        if self.threshold == 0 {
            return;
        }
        let score = start.bound;
        // EXPERIMENT (LOSAT_X_CLUTPAIRS): the row of the query letters at the two innermost
        // positions of `neighbors_below_last` (w - 3, w - 2), read once for the word.
        let row = self.pairs.as_ref().map(|pairs| {
            let word = query_word.letters();
            pairs.row(word[last - 2] as usize, word[last - 1] as usize)
        });
        for letter in letters.clone() {
            let letter_score = self.matrix[last_char][letter as usize];
            if score + letter_score < self.threshold {
                continue;
            }
            let partial = (letter - letters.start) * self.power[last];
            match row {
                Some(row) => self.neighbors_below_last_paired(
                    row,
                    query_word,
                    score + letter_score,
                    0,
                    partial,
                    sink,
                ),
                None => {
                    self.neighbors_below_last(query_word, score + letter_score, 0, partial, sink)
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053 (`neighbors`)
    /// `for_each_word_in`'s recursion: `neighbors` with the last letter chosen (in `partial`, the
    /// index of the letters chosen, and in `score`, the bound of the word: the letters chosen, and
    /// `row_max` of the positions left), its positions from `pos` to the one before the last.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    fn neighbors_below_last<Q: XQueryWord + ?Sized>(
        &self,
        word: &Q,
        mut score: i32,
        pos: usize,
        partial: u32,
        sink: &mut impl XSink,
    ) {
        let query_word = word.letters();
        let below_last = self.length_of::<Q>() - 2;
        let curr_query_char = query_word[pos] as usize;
        score -= self.info.row_max[curr_query_char];
        let row_sorted = &self.info.matrix_sorted[curr_query_char][..self.alphabet];
        let threshold = self.threshold;

        if pos == below_last {
            let count =
                row_sorted.partition_point(|&letter_score| score + letter_score >= threshold);
            sink.many(partial, &self.below_last[curr_query_char][..count]);
            return;
        }

        let char_sorted = &self.info.matrix_sorted_char[curr_query_char];
        let scale = self.power[pos];
        if pos + 1 == below_last {
            // The call for the position before the last, inline: the same statements as the branch
            // above.
            let next_char = query_word[below_last] as usize;
            let next_row = &self.info.matrix_sorted[next_char][..self.alphabet];
            let terms = &self.below_last[next_char];
            let next_max = self.info.row_max[next_char];
            for (i, &letter_score) in row_sorted.iter().enumerate() {
                if score + letter_score < threshold {
                    break;
                }
                let next_score = score + letter_score - next_max;
                let count =
                    next_row.partition_point(|&next_letter| next_score + next_letter >= threshold);
                sink.many(partial + u32::from(char_sorted[i]) * scale, &terms[..count]);
            }
            return;
        }
        for (i, &letter_score) in row_sorted.iter().enumerate() {
            if score + letter_score < threshold {
                break;
            }
            self.neighbors_below_last(
                word,
                score + letter_score,
                pos + 1,
                partial + u32::from(char_sorted[i]) * scale,
                sink,
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053 (`neighbors`)
    /// EXPERIMENT (LOSAT_X_CLUTPAIRS): `neighbors` with the cells of the last two positions (w - 2
    /// and w - 1) from `row`, the row of their query letters (a, b) (`XPairTable::row`, read by
    /// `for_each_word` once for the word). `neighbors`' two loops
    /// there give, for the bound `score` (less `row_max[a]`), the letters i of a's sorted row with
    /// `score + row_a[i] >= threshold`, and for each the letters j of b's with `score + row_a[i] -
    /// row_max[b] + row_b[j] >= threshold`: the pairs with `row_a[i] + row_b[j] >= threshold -
    /// (score - row_max[b])` (the second test implies the first, `row_b[j] <= row_max[b]`), a
    /// prefix of the row. Their cells `partial + char_a[i] * A^(w-2) + char_b[j] * A^(w-1)` are
    /// `partial + code * A^(w-2)`. The same cells in another order within the offset, which does
    /// not change the table (`for_each_word_in`).
    fn neighbors_paired<Q: XQueryWord + ?Sized>(
        &self,
        row: &XPairRow,
        word: &Q,
        mut score: i32,
        current_pos: usize,
        partial: u32,
        sink: &mut impl XSink,
    ) {
        let query_word = word.letters();
        // The compressed word sizes are 5 to 7 (`XHeader::new`): the recursion stops at w - 2.
        debug_assert!(current_pos + 2 <= self.length_of::<Q>());
        let curr_query_char = query_word[current_pos] as usize;
        score -= self.info.row_max[curr_query_char];
        let threshold = self.threshold;
        let scale = self.power[current_pos];
        if current_pos + 2 == self.length_of::<Q>() {
            let last_char = query_word[current_pos + 1] as usize;
            let key = threshold.saturating_sub(score - self.info.row_max[last_char]);
            row.emit(key, partial, scale, sink);
            return;
        }
        let row_sorted = &self.info.matrix_sorted[curr_query_char][..self.alphabet];
        let char_sorted = &self.info.matrix_sorted_char[curr_query_char];
        for (i, &letter_score) in row_sorted.iter().enumerate() {
            if score + letter_score < threshold {
                break;
            }
            self.neighbors_paired(
                row,
                word,
                score + letter_score,
                current_pos + 1,
                partial + u32::from(char_sorted[i]) * scale,
                sink,
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053 (`neighbors`)
    /// EXPERIMENT (LOSAT_X_CLUTPAIRS): `neighbors_below_last` with the cells of its two innermost
    /// positions (w - 3 and w - 2) from `row`, the row of their query letters (read by
    /// `for_each_word_from` once for the word), as `neighbors_paired` takes the last two: the
    /// cells `partial + code * A^(w-3)`.
    #[cfg_attr(
        not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )),
        allow(dead_code)
    )]
    fn neighbors_below_last_paired<Q: XQueryWord + ?Sized>(
        &self,
        row: &XPairRow,
        word: &Q,
        mut score: i32,
        pos: usize,
        partial: u32,
        sink: &mut impl XSink,
    ) {
        let query_word = word.letters();
        let below_last = self.length_of::<Q>() - 2;
        // The compressed word sizes are 5 to 7 (`XHeader::new`): the recursion stops at w - 3.
        debug_assert!(pos < below_last);
        let curr_query_char = query_word[pos] as usize;
        score -= self.info.row_max[curr_query_char];
        let threshold = self.threshold;
        let scale = self.power[pos];
        if pos + 1 == below_last {
            let next_char = query_word[below_last] as usize;
            let key = threshold.saturating_sub(score - self.info.row_max[next_char]);
            row.emit(key, partial, scale, sink);
            return;
        }
        let row_sorted = &self.info.matrix_sorted[curr_query_char][..self.alphabet];
        let char_sorted = &self.info.matrix_sorted_char[curr_query_char];
        for (i, &letter_score) in row_sorted.iter().enumerate() {
            if score + letter_score < threshold {
                break;
            }
            self.neighbors_below_last_paired(
                row,
                word,
                score + letter_score,
                pos + 1,
                partial + u32::from(char_sorted[i]) * scale,
                sink,
            );
        }
    }
}

/// The parameters of `BlastCompressedAaLookupTable::new` without its backbone.
struct XHeader {
    word_length: usize,
    alphabet: usize,
    reciprocal_alphabet_size: i32,
    backbone_size: usize,
    compress_table: [u8; BLASTAA_SIZE],
    scaled_compress_table: [i32; BLASTAA_SIZE],
}

impl XHeader {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1294-1317
    // The expressions of `BlastCompressedAaLookupTable::new`.
    fn new(word_size: usize) -> Option<Self> {
        let (alphabet, reciprocal_alphabet_size) = compressed_lookup_parameters(word_size)?;
        let compress_table = build_compressed_translation(alphabet)?;
        let scaled_compress_table =
            build_scaled_compress_table(&compress_table, alphabet, word_size);
        let backbone_size = alphabet.checked_pow(word_size as u32)?.checked_add(1)?;
        Some(Self {
            word_length: word_size,
            alphabet,
            reciprocal_alphabet_size,
            backbone_size,
            compress_table,
            scaled_compress_table,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
/// The loop of `add_location_words`: calls `visit(query_offset, word)` for every query offset that it
/// passes to `add_word_hits`, in order, and stops when `visit` returns `false` (then `false`).
#[inline]
fn x_for_each_offset(
    word_length: usize,
    query: &[u8],
    locations: &[(i32, i32)],
    mut visit: impl FnMut(i32, &[u8]) -> bool,
) -> bool {
    let word_size = word_length as i32;
    for &(left, right) in locations {
        let to = right - word_size + 1;
        for offset in left..=to {
            let start = usize::try_from(offset)
                .expect("NCBI BLAST compressed query offset must be non-negative");
            let end = start + word_length;
            if end <= query.len() && !visit(offset, &query[start..end]) {
                return false;
            }
        }
    }
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
/// `compressed_table_fills_banks` without the table: whether the appends of the words of `locations`,
/// in one table of their own, take more than `overflow_cell_capacity` overflow cells (module header).
fn x_table_fills_banks(
    words: &XWords<'_>,
    backbone_size: usize,
    query: &[u8],
    locations: &[(i32, i32)],
    overflow_cell_capacity: usize,
) -> bool {
    struct XOwnTable {
        counts: Vec<u32>,
        cells: usize,
    }
    // No NCBI counterpart: the batch's own table's counts (`num_used`) and overflow cells.
    impl XSink for XOwnTable {
        fn one(&mut self, index: u32) {
            let c = self.counts[index as usize];
            self.cells += usize::from(x_takes_overflow_cell(c));
            self.counts[index as usize] = c + 1;
        }
    }
    let mut own = XOwnTable {
        counts: vec![0u32; backbone_size],
        cells: 0,
    };
    let fits = x_for_each_offset(words.word_length, query, locations, |_, word| {
        words.for_each_word(word, &mut own);
        own.cells <= overflow_cell_capacity
    });
    !fits
}

// No NCBI counterpart: the most appends the flat table takes; past it the reference builder runs
// (`TooLarge`, the same table). Two facts bound it. A query offset appends to each cell at most
// twice (its word and a neighbor), and the bound is tested after each query offset, so the
// counts and the slices' `u32` offsets stay below `i32::MAX` (`num_used` is an Int4) and the
// record holds at most the bound plus `2 * backbone_size` appends. And a `Vec` of `u32` grows
// by doubling: a full capacity `c` asks for `2c`, which panics ("capacity overflow") when
// `2c * 4` bytes pass `isize::MAX`, so a 32-bit target (isize::MAX = 2^31 - 1) must take the
// fallback before the record gets there. `X_APPENDS_LIMIT` is the smaller of the two bounds
// on the record's length: `i32::MAX` on a 64-bit target, 268,435,455 on a 32-bit one.
const X_APPENDS_LIMIT: usize = {
    let counts = i32::MAX as usize;
    let growth = isize::MAX as usize / (2 * core::mem::size_of::<u32>());
    if counts < growth {
        counts
    } else {
        growth
    }
};

fn x_max_appends(backbone_size: usize) -> usize {
    X_APPENDS_LIMIT.saturating_sub(2 * backbone_size)
}

/// EXPERIMENT (LOSAT_X_BPCLUT): `try_build_blosum62_compressed_lookup` as a flat table. With
/// LOSAT_X_BPFFILL the bucket fill builds it (the count-and-fill builder where it declines).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// The builds use `x_try_build_flat_reserving`; the tests build without a count.
#[cfg(test)]
pub(crate) fn x_try_build_flat_compressed_lookup(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
) -> Option<XFlatBuild> {
    x_try_build_flat_reserving(word_size, scaled_threshold, query, locations, tables, 0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_try_build_flat_compressed_lookup`)
/// `x_try_build_flat_compressed_lookup` with `appends`, the most appends the build makes as its
/// caller counted them (0 where it did not): the builders reserve their record of the appends
/// once instead of growing it (the same table).
fn x_try_build_flat_reserving(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    appends: usize,
) -> Option<XFlatBuild> {
    if x_bpffill_on() {
        match x_try_build_flat_by_buckets(
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            appends,
        ) {
            XBucketBuild::Built(built) => {
                crate::utils::xreach::x_reach!("LOSAT_X_BPFFILL", true);
                return built;
            }
            XBucketBuild::Declined => crate::utils::xreach::x_reach!("LOSAT_X_BPFFILL", false),
        }
    }
    x_try_build_flat_with_capacity(
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        NCBI_OVERFLOW_CELL_CAPACITY,
        X_COUNT_CHUNK,
        appends,
    )
}

// No NCBI counterpart: how many appends ahead the count pass reads the count of a cell, and the fill
// pass the cursor of a cell and then the slot it points to (cache hints only, `x_prefetch`).
const X_COUNT_AHEAD: usize = 32;
const X_FILL_STARTS_AHEAD: usize = 32;
const X_FILL_HITS_AHEAD: usize = 16;
// The fill reads a cursor for the hits hint after the starts hint has fetched its line.
const _: () = assert!(X_FILL_HITS_AHEAD <= X_FILL_STARTS_AHEAD);

// No NCBI counterpart: how many appends the count pass records before it counts them (and tests NCBI's
// table, as the reference does after each query offset); it bounds the record kept before a refusal.
const X_COUNT_CHUNK: usize = 1 << 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop over the
// segments, see `x_for_each_offset`)
/// The segments (ranges of `locations`) before, inside and between NCBI's tables, in order, each with
/// whether it is one of NCBI's tables (the words between tables, a split query batch, are added
/// without a test of the overflow banks).
fn x_count_ranges(
    tables: &[std::ops::Range<usize>],
    locations: usize,
) -> Vec<(std::ops::Range<usize>, bool)> {
    let mut ranges: Vec<(std::ops::Range<usize>, bool)> = Vec::with_capacity(2 * tables.len() + 1);
    let mut next = 0usize;
    for table in tables {
        ranges.push((next..table.start, false));
        ranges.push((table.clone(), true));
        next = table.end;
    }
    ranges.push((next..locations, false));
    ranges
}

/// The count pass's state: the cell of each append in order, (query offset, appends) of each query
/// offset with any, each cell's count of the appends counted so far, and the overflow cells that
/// those appends take in the reference's table of all batches.
struct XCount {
    targets: Vec<u32>,
    per_offset: Vec<(i32, u32)>,
    counts: Vec<u32>,
    counted: usize,
    cells: usize,
}

impl XCount {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 (module header)
    /// Counts the appends recorded since the last call, in order.
    fn count_new(&mut self) {
        let new = &self.targets[self.counted..];
        for (k, &index) in new.iter().enumerate() {
            if let Some(&ahead) = new.get(k + X_COUNT_AHEAD) {
                x_prefetch(self.counts.as_ptr().wrapping_add(ahead as usize));
            }
            let c = self.counts[index as usize];
            self.cells += usize::from(x_takes_overflow_cell(c));
            self.counts[index as usize] = c + 1;
        }
        self.counted = self.targets.len();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop, see
// `x_for_each_offset`), :787-849 (the appends, module header), :1182-1242 (`finalize`, below)
// The passes of the module header; the tests of NCBI's tables are those of
// `try_build_compressed_lookup_with_capacity`, made after each chunk of appends and at the end
// of a table instead of after each query offset (the reference tests a table once, at the first query
// offset where its appends have taken more cells of the table of all batches than the capacity, and its
// outcome, the batch's own table, does not depend on that offset; the cells taken only grow, so a table
// that passes that bound passes it at the end of the table).
#[allow(clippy::too_many_arguments)]
fn x_try_build_flat_with_capacity(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    count_chunk: usize,
    appends: usize,
) -> Option<XFlatBuild> {
    x_flat_with_capacity_of(
        Blosum62Compressed,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        overflow_cell_capacity,
        count_chunk,
        appends,
    )
}

/// `x_try_build_flat_with_capacity` with the words of the search matrix `matrix`
/// (`build_compressed_score_matrix`). `None` also for IDENTITY, whose table NCBI cannot build.
#[allow(clippy::too_many_arguments)]
fn x_try_build_flat_of(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    count_chunk: usize,
) -> Option<XFlatBuild> {
    x_flat_with_capacity_of(
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        overflow_cell_capacity,
        count_chunk,
        0,
    )
}

// The body of both builders above, so that BLOSUM62's keeps its constant matrix.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn x_flat_with_capacity_of<S: CompressedScores>(
    scores: S,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    count_chunk: usize,
    appends: usize,
) -> Option<XFlatBuild> {
    let header = XHeader::new(word_size)?;
    let matrix = scores.score_matrix(header.alphabet)?;
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let backbone_size = header.backbone_size;
    let max_appends = x_max_appends(backbone_size);

    // Count pass, over the segments before, inside and between NCBI's tables (the words between
    // tables, a split query batch, are added without a test).
    let ranges = x_count_ranges(tables, locations.len());
    // With `appends` (the caller's count), the record is reserved once; an offset has at most one
    // entry of `per_offset`.
    let mut st = XCount {
        targets: Vec::with_capacity(appends.min(max_appends + 1)),
        per_offset: Vec::with_capacity(if appends > 0 { query.len() } else { 0 }),
        // One more entry: the counts become the slices' starts.
        counts: vec![0u32; backbone_size + 1],
        counted: 0,
        cells: 0,
    };
    for (range, is_table) in ranges {
        st.count_new();
        let cells_before = st.cells;
        let words_before = range.start > 0;
        let range_locations = &locations[range];
        let mut tested = false;
        let mut refused = false;
        let mut too_large = false;
        let mut test = |st: &mut XCount| -> bool {
            st.count_new();
            if !is_table || tested || st.cells - cells_before <= overflow_cell_capacity {
                return true;
            }
            tested = true;
            refused = !words_before
                || x_table_fills_banks(
                    &words,
                    backbone_size,
                    query,
                    range_locations,
                    overflow_cell_capacity,
                );
            !refused
        };
        let completed = x_for_each_offset(word_size, query, range_locations, |offset, word| {
            let before = st.targets.len();
            words.for_each_word(word, &mut st.targets);
            let appends = st.targets.len() - before;
            if appends != 0 {
                st.per_offset.push((offset, appends as u32));
            }
            if st.targets.len() > max_appends {
                too_large = true;
                return false;
            }
            st.targets.len() - st.counted < count_chunk || test(&mut st)
        });
        if completed {
            test(&mut st);
        }
        if refused {
            return Some(XFlatBuild::Refused);
        }
        if too_large {
            return Some(XFlatBuild::TooLarge);
        }
    }
    st.count_new();
    let XCount {
        targets,
        per_offset,
        counts: mut starts,
        ..
    } = st;

    // Layout.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1197-1235
    // The statements of `finalize` on the counts (`num_used`); then each cell's start (the sum of the
    // counts before it) goes to the entry after the cell's, where the fill pass uses it as the cell's
    // cursor, and the first entry becomes 0.
    let nonempty_count = starts[..backbone_size]
        .iter()
        .filter(|&&count| count > 0)
        .count();
    let mut pv_array_bts = PV_ARRAY_BTS;
    if (nonempty_count as f64) <= 0.01 * (backbone_size as f64) {
        pv_array_bts += ilog2(backbone_size / (8 * X_TARGET_PV_BYTES));
    }
    let mut pv = vec![0u32; (backbone_size >> pv_array_bts) + 1];
    let mut longest_chain = 0i32;
    // The cells whose scan order is not the order of their appends (more than 5 offsets).
    let mut overflowed: Vec<u32> = Vec::new();
    let mut total = 0u32;
    // Entry `index + 1` holds cell `index + 1`'s count until cell `index`'s start replaces it.
    let mut count = std::mem::take(&mut starts[0]);
    for index in 0..backbone_size {
        let next = std::mem::replace(&mut starts[index + 1], total);
        if count > 0 {
            compressed_pv_set(&mut pv, index, pv_array_bts);
            longest_chain = longest_chain.max(count as i32);
            if count > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1 {
                overflowed.push(index as u32);
            }
        }
        total += count;
        count = next;
    }

    // Fill pass: the appends again in order, each at the next slot of its cell (the cursor of cell
    // `index` is entry `index + 1`), so each cell holds o_0..o_{n-1}; a cursor ends at the next cell's
    // start, so then entry `index` is cell `index`'s start and the last entry the number of hits. Then
    // each cell of more than 5 offsets puts o_k at `x_emission_position(k, n)` (cells of 5 or fewer:
    // position k).
    let mut hits = vec![0i32; total as usize];
    let mut at = 0usize;
    for &(query_offset, appends) in &per_offset {
        let end = at + appends as usize;
        for (k, &index) in targets[at..end].iter().enumerate() {
            if let Some(&ahead) = targets.get(at + k + X_FILL_STARTS_AHEAD) {
                x_prefetch(starts.as_ptr().wrapping_add(ahead as usize + 1));
            }
            if let Some(&near) = targets.get(at + k + X_FILL_HITS_AHEAD) {
                if let Some(&slot) = starts.get(near as usize + 1) {
                    x_prefetch(hits.as_ptr().wrapping_add(slot as usize));
                }
            }
            let cursor = &mut starts[index as usize + 1];
            hits[*cursor as usize] = query_offset;
            *cursor += 1;
        }
        at = end;
    }
    drop(targets);
    drop(per_offset);
    for &index in &overflowed {
        let index = index as usize;
        x_to_scan_order(&mut hits[starts[index] as usize..starts[index + 1] as usize]);
    }

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

// No NCBI counterpart: reads the LOSAT_X_BPFFILL switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPFFILL): a strict-set switch, on unless it is `0`: the serial builder of the
/// flat table fills it by buckets of cells (`x_try_build_flat_by_buckets`).
fn x_bpffill_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPFFILL"))
}

// No NCBI counterpart: the width of the bucket fill's buckets (LOSAT_X_BPFFILL), 2^10 cells: a
// bucket's counts and cursors (4 KiB) stay in L1 and the partition writes 742 streams for word
// size 5 (2^8 to 2^14 timed the same or slower, v4-BPF K2-2). Scheduling only.
const X_BUCKET_SHIFT: u32 = 10;

// No NCBI counterpart: whether a query of `letters` letters fits the bucket fill's packed append,
// the query offset above the cell's index in its bucket (offsets below 2^(32 - X_BUCKET_SHIFT)).
fn x_bucket_packs(letters: usize) -> bool {
    usize::BITS - letters.leading_zeros() + X_BUCKET_SHIFT <= u32::BITS
}

/// What the bucket fill made.
#[expect(
    clippy::large_enum_variant,
    reason = "boxing a variant changes the layout and the memory; that needs an RSS and time measurement"
)]
enum XBucketBuild<H: XQueryOffset = i32> {
    /// `x_try_build_flat_with_capacity`'s outcome: its table, its `TooLarge`, or `None` where the
    /// table's header cannot be made.
    Built(Option<XFlatBuild<H>>),
    /// The count-and-fill builder decides (`x_try_build_flat_by_buckets`).
    Declined,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop, see
// `x_for_each_offset`), :787-849 (the appends, module header), :1182-1242 (`finalize`, below)
/// EXPERIMENT (LOSAT_X_BPFFILL): `x_try_build_flat_with_capacity`'s outcome with the appends sorted
/// by buckets of 2^`X_BUCKET_SHIFT` cells before they are counted and placed (module header).
/// Declined where a table's appends pass `overflow_cell_capacity`, where the query has 2^22 letters
/// or more (the packed append) and where the PV shift depends on the counts (word sizes other than
/// 5, which `-task blastp-fast` does not take). Out of line, so that the count-and-fill builder
/// keeps its code.
#[inline(never)]
fn x_try_build_flat_by_buckets(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    appends: usize,
) -> XBucketBuild {
    x_buckets_on(
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        overflow_cell_capacity,
        false,
        appends,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1318-1320,1359-1373
// NCBI allocates a zeroed table for each query chunk and frees it with the chunk. Buffer reuse only:
// the kept pv is cleared and zero-filled as calloc gives it; the kept starts and hits, every entry of
// which the bucket fill writes before any is read, are only sized (`x_sized_on`).
/// EXPERIMENT (LOSAT_X_BXCLUTPOOL): whether BLASTX builds its flat compressed tables on vectors kept
/// between builds. Every thread that built a table kept a table's worth and its build's record of
/// the appends in its glibc arena after they were freed, about 20 MiB per thread on a split query.
pub(crate) fn x_bx_clut_pool() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCLUTPOOL"))
}

/// No NCBI counterpart: the vectors of a flat table kept for the next build (LOSAT_X_BXCLUTPOOL).
#[derive(Default)]
struct XFlatTableVectors<H = i32> {
    starts: Vec<u32>,
    hits: Vec<H>,
    pv: Vec<u32>,
}

// No NCBI counterpart: the kept sets and the number of sets lent to pooled builds whose tables have
// not come back. A set comes back only against one lent, so a table built without the pool
// (`LOSAT_X_BPFFILL=0`, a fill declined before it took a set) adds nothing, and the pool never holds
// more sets than pooled builds and their tables were alive at once.
struct XFlatPool<H = i32> {
    kept: Vec<XFlatTableVectors<H>>,
    lent: usize,
}

impl<H: Default> XFlatPool<H> {
    fn lend(&mut self) -> XFlatTableVectors<H> {
        self.lent += 1;
        self.kept.pop().unwrap_or_default()
    }

    /// Keeps `vectors` against a set lent; with none lent they are refused (and dropped by the caller).
    fn take_back(&mut self, vectors: XFlatTableVectors<H>) -> Result<(), XFlatTableVectors<H>> {
        if self.lent == 0 {
            return Err(vectors);
        }
        self.lent -= 1;
        self.kept.push(vectors);
        Ok(())
    }
}

// No NCBI counterpart: one pool per offset type: BLASTX's `u16` tables, and its `i32` ones of query
// batches of more than `X_U16_LETTERS` letters.
static X_FLAT_TABLES: Mutex<XFlatPool> = Mutex::new(XFlatPool {
    kept: Vec::new(),
    lent: 0,
});
static X_NARROW_FLAT_TABLES: Mutex<XFlatPool<u16>> = Mutex::new(XFlatPool {
    kept: Vec::new(),
    lent: 0,
});

/// The pool of the tables of query offsets `Self`.
trait XPooledOffset: Sized + 'static {
    fn x_pool() -> &'static Mutex<XFlatPool<Self>>;
}

impl XPooledOffset for i32 {
    fn x_pool() -> &'static Mutex<XFlatPool<Self>> {
        &X_FLAT_TABLES
    }
}

impl XPooledOffset for u16 {
    fn x_pool() -> &'static Mutex<XFlatPool<Self>> {
        &X_NARROW_FLAT_TABLES
    }
}

fn x_flat_pool<H: XPooledOffset>() -> std::sync::MutexGuard<'static, XFlatPool<H>> {
    H::x_pool().lock().unwrap_or_else(PoisonError::into_inner)
}

// No NCBI counterpart for the size: the record of a bucket fill's appends, one entry per hit of the
// table, is needed only during the build. Kept between builds, three builds side by side held three
// records through the search (ABRicate's call at 8 threads: 48 MiB instead of 40); freed from a
// thread's arena, it stayed resident there, as the tables did. glibc maps an allocation of at least
// its mmap threshold, but it raises the threshold to the size of each mapped block that is freed (up
// to 32 MiB), so its capacity is above that ceiling: it is always mapped and returned when dropped,
// and the pages beyond its length are never touched and are not resident (as
// `tblastn/x_subject_translate.rs` `x_frame_capacity`). Other allocators return large blocks on their
// own, and a wasm memory would keep the capacity. BLASTX's `u16` fill keeps no record; the chunks of
// its appends (`XBucketChunks`, at most 63 entries more per bucket) are cut from an arena of this
// capacity for the same reason.
/// EXPERIMENT (LOSAT_X_BXCLUTPOOL): the capacity, in entries, of the record of a pooled bucket fill.
const fn x_mapped_record_entries() -> usize {
    x_mapped_entries::<u32>()
}

/// The fewest entries of type `T` that glibc maps whatever its dynamic mmap threshold (above), 0
/// with another allocator.
const fn x_mapped_entries<T>() -> usize {
    if cfg!(all(target_os = "linux", target_env = "gnu")) {
        // glibc's largest dynamic mmap threshold (DEFAULT_MMAP_THRESHOLD_MAX on 64-bit), and one
        // page more
        ((32 << 20) + 4096) / std::mem::size_of::<T>()
    } else {
        0
    }
}

// No NCBI counterpart: `vector` cleared and zero-filled to `n` entries (as calloc gives them); one
// short of room is replaced by a zeroed one of `n` entries.
fn x_zeroed_on<T: Copy + Default>(mut vector: Vec<T>, n: usize) -> Vec<T> {
    if vector.capacity() < n {
        return vec![T::default(); n];
    }
    vector.clear();
    vector.resize(n, T::default());
    vector
}

// No NCBI counterpart: `vector` with `n` entries for a fill that writes every entry before any is
// read (the bucket fill's starts and hits): the values the table before left are not cleared, only
// entries added past them are zeroed. One short of room is replaced by a zeroed one (calloc'd, so
// that a large one is mapped and its pages are touched only when the fill writes them) of `n`
// entries, or with `slack` (the hits, whose number changes from table to table) of `n + n / 4`
// entries cut to `n`. Zero-filling them again cost the small searches 1-5 % of their user cycles.
// `S` is the bucket fill's matrix type (`x_buckets_with`), so that BLOSUM62's builders keep the
// callers, and the inlining, of their own instances.
#[expect(
    clippy::extra_unused_type_parameters,
    reason = "S only keeps one instance per bucket-fill matrix type (see above)"
)]
fn x_sized_on<S, T: Copy + Default>(mut vector: Vec<T>, n: usize, slack: bool) -> Vec<T> {
    if vector.capacity() < n {
        let mut fresh = vec![T::default(); if slack { n + n / 4 } else { n }];
        fresh.truncate(n);
        return fresh;
    }
    vector.resize(n, T::default());
    vector
}

// No NCBI counterpart: buffer reuse (LOSAT_X_BXCLUTPOOL); NCBI frees the vectors.
/// Keeps the vectors of a flat `table` for the next pooled build (the table is left empty), if a set
/// is lent ([`XFlatPool`]); otherwise the table keeps them.
pub(crate) fn x_give_back_flat(table: &mut XFlatCompressedLookup<u16>) {
    x_give_back_to_pool(table);
}

/// `x_give_back_flat` for a table of `i32` offsets.
pub(crate) fn x_give_back_wide_flat(table: &mut XFlatCompressedLookup) {
    x_give_back_to_pool(table);
}

fn x_give_back_to_pool<H: XQueryOffset + XPooledOffset>(table: &mut XFlatCompressedLookup<H>) {
    let mut pool = x_flat_pool::<H>();
    if pool.lent == 0 {
        return;
    }
    // Only the vectors an arena would keep; the table keeps (and drops) the others.
    let vectors = XFlatTableVectors {
        starts: x_take_arena_sized(&mut table.starts),
        hits: x_take_arena_sized(&mut table.hits),
        pv: x_take_arena_sized(&mut table.pv),
    };
    let kept = pool.take_back(vectors);
    debug_assert!(kept.is_ok());
}

// No NCBI counterpart: memory placement only (line NP, SC9d); NCBI frees the vectors.
/// `vector`'s buffer taken for the pool (`vector` left empty), unless glibc maps a block of its
/// capacity whatever its dynamic mmap threshold (more than its largest one,
/// DEFAULT_MMAP_THRESHOLD_MAX on 64-bit): then nothing is taken (an empty vector), and the buffer,
/// dropped with its table, goes back to the system. The pool is for the vectors that a pool thread's
/// arena would keep after they are freed; kept, the larger ones stayed resident between builds
/// (word size 7's 40 MB of starts per table: 800 cDNAs x scer_pep at 8 threads, 1.41-1.46 GB of
/// maximum RSS with them in the pool, 1.06-1.15 GB without the pool).
fn x_take_arena_sized<T>(vector: &mut Vec<T>) -> Vec<T> {
    let bytes = vector.capacity().saturating_mul(std::mem::size_of::<T>());
    if cfg!(all(target_os = "linux", target_env = "gnu")) && bytes > 32 << 20 {
        return Vec::new();
    }
    std::mem::take(vector)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (see
// `x_try_build_flat_by_buckets`)
/// `x_try_build_flat_by_buckets`, with `pooled` on vectors kept between builds
/// (`LOSAT_X_BXCLUTPOOL`); the table is the same.
#[allow(clippy::too_many_arguments)]
fn x_buckets_on(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    pooled: bool,
    appends: usize,
) -> XBucketBuild {
    x_buckets_with(
        Blosum62Compressed,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        overflow_cell_capacity,
        pooled,
        appends,
    )
}

/// `x_try_build_flat_by_buckets` with the words of the search matrix `matrix`
/// (`build_compressed_score_matrix`), as `x_try_build_flat_of` builds them.
fn x_buckets_of(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
) -> XBucketBuild {
    x_buckets_with(
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        overflow_cell_capacity,
        false,
        0,
    )
}

// The body of both builders above, so that BLOSUM62's keeps its constant matrix.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn x_buckets_with<S: CompressedScores>(
    scores: S,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    pooled: bool,
    appends: usize,
) -> XBucketBuild {
    if !x_bucket_packs(query.len()) {
        return XBucketBuild::Declined;
    }
    let Some(header) = XHeader::new(word_size) else {
        return XBucketBuild::Built(None);
    };
    let backbone_size = header.backbone_size;
    // `finalize`'s larger PV shift is `PV_ARRAY_BTS` (word size 5): the layout sets the PV bits.
    if ilog2(backbone_size / (8 * X_TARGET_PV_BYTES)) != 0 {
        return XBucketBuild::Declined;
    }
    let Some(matrix) = scores.score_matrix(header.alphabet) else {
        return XBucketBuild::Built(None);
    };
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let max_appends = x_max_appends(backbone_size);

    // The vectors in the count-and-fill builder's order: the starts (its counts) first, the hits
    // after the appends; the record of the appends is dropped at the end, the offsets' appends
    // after the partition (they bound the list of a cell's offsets moved to the scan order: a query
    // offset appends to a cell at most twice).
    let kept = if pooled {
        x_flat_pool::<i32>().lend()
    } else {
        XFlatTableVectors::default()
    };
    // One more entry: the number of hits.
    let mut starts = x_sized_on::<S, _>(kept.starts, backbone_size + 1, false);

    // Appends, over the segments before, inside and between NCBI's tables, recorded as the
    // count-and-fill builder records them; not counted.
    let ranges = x_count_ranges(tables, locations.len());
    // With `appends` (the caller's count), the record is reserved once; an offset has at most one
    // entry of `per_offset`.
    let mut targets: Vec<u32> = if pooled {
        Vec::with_capacity(x_mapped_record_entries())
    } else {
        Vec::with_capacity(appends.min(max_appends + 1))
    };
    let mut per_offset: Vec<(i32, u32)> =
        Vec::with_capacity(if appends > 0 { query.len() } else { 0 });
    // The kept vectors taken go back to the pool at a refusal.
    let give_back = |starts: Vec<u32>, hits: Vec<i32>, pv: Vec<u32>| {
        if pooled {
            let refused = x_flat_pool::<i32>()
                .take_back(XFlatTableVectors { starts, hits, pv })
                .err();
            drop(refused);
        }
    };
    for (range, is_table) in ranges {
        let range_start = targets.len();
        let mut too_large = false;
        let mut declined = false;
        x_for_each_offset(word_size, query, &locations[range], |offset, word| {
            let before = targets.len();
            words.for_each_word(word, &mut targets);
            let appends = targets.len() - before;
            if appends != 0 {
                per_offset.push((offset, appends as u32));
            }
            if targets.len() > max_appends {
                too_large = true;
                return false;
            }
            // An append takes at most one overflow cell: within the capacity the table's test
            // passes.
            if is_table && targets.len() - range_start > overflow_cell_capacity {
                declined = true;
                return false;
            }
            true
        });
        if too_large || declined {
            give_back(starts, kept.hits, kept.pv);
            return if too_large {
                XBucketBuild::Built(Some(XFlatBuild::TooLarge))
            } else {
                XBucketBuild::Declined
            };
        }
    }
    let total = targets.len();

    // Partition: the buckets' appends in order (a counting sort, stable), bucket b's at
    // `bucket_starts[b]..bucket_starts[b + 1]` of `hits`, the slots its cells' hits take.
    let buckets = backbone_size.div_ceil(1 << X_BUCKET_SHIFT);
    let mask = (1u32 << X_BUCKET_SHIFT) - 1;
    let mut bucket_starts = vec![0u32; buckets + 1];
    for &index in &targets {
        bucket_starts[(index >> X_BUCKET_SHIFT) as usize + 1] += 1;
    }
    for b in 0..buckets {
        bucket_starts[b + 1] += bucket_starts[b];
    }
    let mut pv = x_zeroed_on(kept.pv, (backbone_size >> PV_ARRAY_BTS) + 1);
    let mut hits = x_sized_on::<S, _>(kept.hits, total, pooled);
    {
        let mut cursors = bucket_starts[..buckets].to_vec();
        let out = hits.as_mut_ptr();
        let mut at = 0usize;
        for &(query_offset, appends) in &per_offset {
            debug_assert!((query_offset as u32) < 1 << (u32::BITS - X_BUCKET_SHIFT));
            let high = (query_offset as u32) << X_BUCKET_SHIFT;
            let end = at + appends as usize;
            for &index in &targets[at..end] {
                let cursor = &mut cursors[(index >> X_BUCKET_SHIFT) as usize];
                // SAFETY: the offsets' appends are `targets` in order, and bucket b's cursor starts
                // at `bucket_starts[b]` and moves on once per append to the bucket, which the
                // histogram counted, so it stays below `bucket_starts[b + 1] <= total = hits.len()`.
                unsafe {
                    out.add(*cursor as usize)
                        .write((high | (index & mask)) as i32)
                };
                *cursor += 1;
            }
            at = end;
        }
    }
    drop(per_offset);

    // Per bucket: its appends copied into `targets` (read no more), the counts of its cells, their
    // starts, the PV bits and `finalize`'s statistics, the fill, then the scan order of the cells
    // of more than 5 offsets.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1197-1235
    // `pv_array_bts` is `PV_ARRAY_BTS` whatever `count` (the larger shift is the same, tested
    // above), so the bits are set as the counts are read.
    let mut cell_slots = vec![0u32; 1 << X_BUCKET_SHIFT];
    let mut overflowed: Vec<u32> = Vec::new();
    let mut nonempty_count = 0usize;
    let mut longest_chain = 0i32;
    for b in 0..buckets {
        let first = b << X_BUCKET_SHIFT;
        let cells = (backbone_size - first).min(1 << X_BUCKET_SHIFT);
        let lo = bucket_starts[b] as usize;
        let hi = bucket_starts[b + 1] as usize;
        let part = &mut targets[..hi - lo];
        for (to, &from) in part.iter_mut().zip(&hits[lo..hi]) {
            *to = from as u32;
        }
        // Each cell's count, then its cursor.
        let slots = &mut cell_slots[..cells];
        slots.fill(0);
        for &packed in part.iter() {
            slots[(packed & mask) as usize] += 1;
        }
        overflowed.clear();
        let mut slot = lo as u32;
        for (cell, (start, count)) in starts[first..first + cells]
            .iter_mut()
            .zip(slots.iter_mut())
            .enumerate()
        {
            let n = *count;
            *start = slot;
            if n > 0 {
                nonempty_count += 1;
                longest_chain = longest_chain.max(n as i32);
                compressed_pv_set(&mut pv, first + cell, PV_ARRAY_BTS);
                if n > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1 {
                    overflowed.push(cell as u32);
                }
            }
            *count = slot;
            slot += n;
        }
        let out = hits.as_mut_ptr();
        for &packed in part.iter() {
            let at = &mut slots[(packed & mask) as usize];
            // SAFETY: a cell's cursor starts at its start and moves on once per append to the
            // cell, which the counts counted, so it stays within the bucket's slots `lo..hi`.
            unsafe {
                out.add(*at as usize)
                    .write((packed >> X_BUCKET_SHIFT) as i32)
            };
            *at += 1;
        }
        // A cursor ends at the cell's end.
        for &cell in &overflowed {
            let cell = cell as usize;
            let begin = starts[first + cell] as usize;
            let end = slots[cell] as usize;
            x_to_scan_order(&mut hits[begin..end]);
        }
    }
    starts[backbone_size] = total as u32;
    drop(targets);
    let mut pv_array_bts = PV_ARRAY_BTS;
    if (nonempty_count as f64) <= 0.01 * (backbone_size as f64) {
        pv_array_bts += ilog2(backbone_size / (8 * X_TARGET_PV_BYTES));
    }

    XBucketBuild::Built(Some(XFlatBuild::Table(XFlatCompressedLookup {
        word_length: header.word_length as i32,
        compressed_alphabet_size: header.alphabet as i32,
        reciprocal_alphabet_size: header.reciprocal_alphabet_size,
        longest_chain,
        pv,
        pv_array_bts: pv_array_bts as i32,
        scaled_compress_table_u8: build_scaled_compress_table_u8(&header.scaled_compress_table),
        starts,
        hits,
    })))
}

// No NCBI counterpart: the entries of a chunk of `XBucketChunks` (256 bytes).
const X_BX_CHUNK: u32 = 64;

// No NCBI counterpart: BLASTX's bucket fill (`x_bx_buckets_on`) keeps no record of the appends.
/// The appends of BLASTX's bucket fill, each written as it is made into the chunks of its bucket
/// of 2^`X_BUCKET_SHIFT` cells, as `query offset << 16 | the cell's index in its bucket`. The
/// chunks of all buckets are cut from one arena; a bucket's chunks are linked in the order of
/// their entries, which is the order of the appends.
struct XBucketChunks {
    arena: Vec<u32>,
    /// Per bucket: the arena index its next entry takes. An index at a chunk's start (0 before
    /// the bucket's first entry) opens a chunk: a bucket's next entry is never the first of the
    /// chunk it fills.
    next_entry: Vec<u32>,
    /// Per bucket: its first and last chunks (`u32::MAX` before the first) and its number of
    /// chunks.
    first: Vec<u32>,
    last: Vec<u32>,
    chunks: Vec<u32>,
    /// Per chunk: the next chunk of its bucket.
    link: Vec<u32>,
    /// The query offset of the appends being made, above the cells' indices.
    high: u32,
    total: usize,
}

impl XBucketChunks {
    fn new(buckets: usize, arena: Vec<u32>) -> Self {
        Self {
            arena,
            next_entry: vec![0; buckets],
            first: vec![u32::MAX; buckets],
            last: vec![u32::MAX; buckets],
            chunks: vec![0; buckets],
            link: Vec::new(),
            high: 0,
            total: 0,
        }
    }

    /// A new chunk for bucket `b`: the arena index of its first entry.
    #[cold]
    #[inline(never)]
    fn open(&mut self, b: usize) -> u32 {
        let chunk = self.link.len() as u32;
        self.arena.resize(self.arena.len() + X_BX_CHUNK as usize, 0);
        self.link.push(u32::MAX);
        match self.last[b] {
            u32::MAX => self.first[b] = chunk,
            last => self.link[last as usize] = chunk,
        }
        self.last[b] = chunk;
        self.chunks[b] += 1;
        chunk * X_BX_CHUNK
    }

    /// The entries of bucket `b`.
    fn count(&self, b: usize) -> u32 {
        match self.chunks[b] {
            0 => 0,
            n => (n - 1) * X_BX_CHUNK + ((self.next_entry[b] - 1) & (X_BX_CHUNK - 1)) + 1,
        }
    }

    /// `visit(entry)` for the entries of bucket `b`, in order.
    #[inline(always)]
    fn for_each(&self, b: usize, mut visit: impl FnMut(u32)) {
        let mut left = self.count(b);
        let mut chunk = self.first[b];
        while left != 0 {
            let n = left.min(X_BX_CHUNK);
            let at = (chunk * X_BX_CHUNK) as usize;
            for &entry in &self.arena[at..at + n as usize] {
                visit(entry);
            }
            left -= n;
            chunk = self.link[chunk as usize];
        }
    }
}

impl XSink for XBucketChunks {
    #[inline(always)]
    fn one(&mut self, index: u32) {
        let b = (index >> X_BUCKET_SHIFT) as usize;
        let mut at = self.next_entry[b];
        if at & (X_BX_CHUNK - 1) == 0 {
            at = self.open(b);
        }
        self.arena[at as usize] = self.high | (index & ((1 << X_BUCKET_SHIFT) - 1));
        self.next_entry[b] = at + 1;
        self.total += 1;
    }

    // The trait's loop, inlined into `XWords`' neighbour lists: with the bucket fill's several
    // instances (word sizes, offset types, runs) the compiler kept it a call per list, 1.3 % more
    // instructions of BLASTX word size 5 at 1 thread than dev's one instance (line NP, SC9e).
    #[inline(always)]
    fn many(&mut self, partial: u32, last: &[u32]) {
        for &letter in last {
            self.one(partial + letter);
        }
    }

    // As `many`, inlined (LOSAT_X_CLUTPAIRS).
    #[inline(always)]
    fn many_coded(&mut self, partial: u32, scale: u32, codes: &[u8]) {
        for &code in codes {
            self.one(partial + u32::from(code) * scale);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (see
// `x_try_build_flat_by_buckets`)
/// BLASTX's bucket fill: `x_buckets_on`'s table with `u16` offsets, with `pooled` on vectors kept
/// between builds (`LOSAT_X_BXCLUTPOOL`), in `parts` (`x_bx_buckets_with`). `WORD` is the word
/// length (`x_bx_build_compressed_lookup`'s), a constant of the instance down to the neighbouring
/// words (`XQueryWord`).
#[allow(clippy::too_many_arguments)]
fn x_bx_buckets_on<const WORD: usize>(
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    pooled: bool,
    parts: usize,
) -> XBucketBuild<u16> {
    x_bx_buckets_with::<_, u16, [u8; WORD]>(
        Blosum62Compressed,
        WORD,
        scaled_threshold,
        query,
        locations,
        tables,
        overflow_cell_capacity,
        pooled,
        parts,
        |word| word.try_into().expect("a word of the word length"),
    )
}

// No NCBI counterpart: a run of the query offsets of a build in parts (`x_bx_buckets_with`); it
// changes no value NCBI computes.
/// The query offsets of the count ranges `ranges` (over `locations`) cut into `parts` runs of about
/// as many offsets each (fewer runs where there are fewer offsets; one empty run where there are
/// none), in order: each run's segments with the index of their count range. One run is the
/// offsets of `locations` in the serial builders' order.
fn x_offset_parts(
    word_size: usize,
    ranges: &[(std::ops::Range<usize>, bool)],
    locations: &[(i32, i32)],
    parts: usize,
) -> Vec<Vec<(usize, (i32, i32))>> {
    let word = word_size as i64;
    let offsets = |&(left, right): &(i32, i32)| {
        usize::try_from(i64::from(right) - word + 2 - i64::from(left)).unwrap_or(0)
    };
    let total: usize = ranges
        .iter()
        .flat_map(|(range, _)| &locations[range.clone()])
        .map(offsets)
        .sum();
    let parts = parts.clamp(1, total.max(1));
    // The first offset (in order) of run `part`.
    let bound = |part: usize| (part * total).div_ceil(parts);
    let mut runs: Vec<Vec<(usize, (i32, i32))>> = vec![Vec::new(); parts];
    let (mut part, mut placed) = (0usize, 0usize);
    for (index, (range, _)) in ranges.iter().enumerate() {
        for &(left, right) in &locations[range.clone()] {
            let to = i64::from(right) - word + 1;
            let mut from = i64::from(left);
            while from <= to {
                while placed >= bound(part + 1) {
                    part += 1;
                }
                let take = ((bound(part + 1) - placed) as i64).min(to - from + 1);
                // Offsets from..=from + take - 1: the segment's right end less what follows them.
                runs[part].push((index, (from as i32, (from + take - 1 + word - 1) as i32)));
                placed += take as usize;
                from += take;
            }
        }
    }
    runs
}

// No NCBI counterpart: why a run of the bucket fill stopped (`x_bx_buckets_with`).
#[derive(Clone, Copy)]
enum XBucketStop {
    /// The appends passed the table's most (`x_max_appends`): the serial builders' `TooLarge`.
    TooLarge,
    /// A table's appends passed the overflow capacity: the bucket fill declines.
    Declined,
}

// No NCBI counterpart: one run's appends (`x_bx_buckets_with`).
struct XBucketRun {
    appends: XBucketChunks,
    /// The appends of each count range the run has offsets of, in order.
    range_appends: Vec<(usize, usize)>,
    stop: Option<XBucketStop>,
}

// No NCBI counterpart: runs `items` on this thread (`parts` 1, or a build without threads), or on
// `parts` tasks of the pool the build runs on (`x_clut_par::x_items_in_parts`), each with its own
// scratch from `init`; the results in item order.
fn x_bucket_items<T, R: Send + Sync>(
    parts: usize,
    items: usize,
    init: impl Fn() -> T + Sync,
    item: impl Fn(&mut T, usize) -> R + Sync,
) -> Vec<R> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parts > 1 {
        return x_clut_par::x_items_in_parts(parts, items, init, item);
    }
    let _ = parts;
    let mut scratch = init();
    (0..items).map(|index| item(&mut scratch, index)).collect()
}

/// No NCBI counterpart: the runs of BLASTX's bucket fill in parts per part (line NP, SC9e), and
/// its bucket ranges per part in the layout and fill; a few per part so that the tasks end close
/// together.
const X_BUCKET_RANGES_PER_PART: usize = 4;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop),
// :787-849 (the appends), :1182-1242 (`finalize`; see `x_buckets_with`)
// NCBI builds a query chunk's table on one thread before the search threads start
// (c++/src/algo/blast/api/prelim_stage.cpp:145-188 starts them for the search).
/// BLASTX's bucket fill: `x_buckets_with`'s table with offsets of type `H`, the words of `scores`,
/// with `pooled` on vectors kept between builds (`LOSAT_X_BXCLUTPOOL`), its chunks cut from a mapped
/// arena (`x_mapped_record_entries`). An offset above a cell's
/// index in its bucket no longer fits 32 bits of a `u16` table's hits, so there is no packed copy of
/// the appends in `hits`: each append goes as it is made to its bucket's chunks (`XBucketChunks`),
/// in place of the record of the appends and its partition, and each bucket's cells are filled from
/// them. Declined as `x_buckets_on`, and where the query has more than `X_U16_LETTERS` letters.
///
/// In `parts` (line NP, SC9e: a BLASTX unit under the cap, `LOSAT_X_BXBUILDPAR`, `blastx/x_bx_mem.rs`;
/// 1: on this thread), the same table slot for slot: the query offsets in order cut into runs
/// (`x_offset_parts`), each run's appends in chunks of its own (the run's appends in order); a
/// bucket's appends are then its runs' in run order, the appends of the serial fill in order, and
/// its slots follow those of the buckets before. The buckets, in ranges of whole buckets (whole PV
/// words: 2^`X_BUCKET_SHIFT` cells), are counted, laid out and filled as on one thread, each range
/// into its own slices of `starts`, `hits` and `pv`. The work is the serial fill's but for a pass
/// over each bucket's runs. Where the runs stop (too many appends, or a table's past the overflow
/// capacity), the serial fill decides at one part (it stops at the first; the runs only show that
/// it stops). `word_of` gives `scores`' query word type of a word of the word length.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn x_bx_buckets_with<S, H, Q>(
    scores: S,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    pooled: bool,
    parts: usize,
    word_of: impl Fn(&[u8]) -> &Q + Sync,
) -> XBucketBuild<H>
where
    S: CompressedScores,
    H: XQueryOffset + XPooledOffset,
    Q: XQueryWord + ?Sized,
{
    if query.len() > X_U16_LETTERS {
        return XBucketBuild::Declined;
    }
    let Some(header) = XHeader::new(word_size) else {
        return XBucketBuild::Built(None);
    };
    let backbone_size = header.backbone_size;
    // `finalize`'s larger PV shift is `PV_ARRAY_BTS` (word size 5): the layout sets the PV bits.
    if ilog2(backbone_size / (8 * X_TARGET_PV_BYTES)) != 0 {
        return XBucketBuild::Declined;
    }
    let Some(matrix) = scores.score_matrix(header.alphabet) else {
        return XBucketBuild::Built(None);
    };
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let max_appends = x_max_appends(backbone_size);
    let buckets = backbone_size.div_ceil(1 << X_BUCKET_SHIFT);
    let mask = (1u32 << X_BUCKET_SHIFT) - 1;
    let ranges = x_count_ranges(tables, locations.len());

    // Appends, over the segments before, inside and between NCBI's tables, in the count-and-fill
    // builder's order, into their buckets' chunks; not counted. A run stops where the serial fill
    // would: its appends past `max_appends` (the runs' together, `appended`, at more than one), or
    // a table's past the overflow capacity.
    let runs = |parts: usize| -> Result<Vec<XBucketRun>, XBucketStop> {
        let cuts = x_offset_parts(word_size, &ranges, locations, parts);
        let parts = cuts.len();
        let appended = std::sync::atomic::AtomicUsize::new(0);
        let stopped = std::sync::atomic::AtomicBool::new(false);
        let runs = x_bucket_items(
            parts,
            parts,
            || (),
            |(), part| {
                // The arena mapped whichever thread fills it (`x_mapped_record_entries`).
                let mut appends =
                    XBucketChunks::new(buckets, Vec::with_capacity(x_mapped_record_entries()));
                let mut range_appends = Vec::new();
                let mut stop = None;
                // The open count range, the appends at its start, and those not yet added to
                // `appended`.
                let (mut open, mut range_start, mut reported) = (usize::MAX, 0usize, 0usize);
                for &(index, segment) in &cuts[part] {
                    if index != open {
                        if open != usize::MAX {
                            range_appends.push((open, appends.total - range_start));
                        }
                        (open, range_start) = (index, appends.total);
                    }
                    let is_table = ranges[index].1;
                    let whole = x_for_each_offset(word_size, query, &[segment], |offset, word| {
                        appends.high = (offset as u32) << 16;
                        words.for_each_word(word_of(word), &mut appends);
                        if appends.total > max_appends {
                            stop = Some(XBucketStop::TooLarge);
                            return false;
                        }
                        // An append takes at most one overflow cell: within the capacity the
                        // table's test passes.
                        if is_table && appends.total - range_start > overflow_cell_capacity {
                            stop = Some(XBucketStop::Declined);
                            return false;
                        }
                        if parts > 1 && appends.total - reported >= X_COUNT_CHUNK {
                            let more =
                                appends.total - std::mem::replace(&mut reported, appends.total);
                            if appended.fetch_add(more, Ordering::Relaxed) + more > max_appends {
                                stop = Some(XBucketStop::TooLarge);
                            }
                            if stop.is_some() || stopped.load(Ordering::Relaxed) {
                                return false;
                            }
                        }
                        true
                    });
                    if !whole {
                        stopped.store(true, Ordering::Relaxed);
                        // Another run's stop: the serial fill decides (below).
                        stop.get_or_insert(XBucketStop::TooLarge);
                        break;
                    }
                }
                if open != usize::MAX {
                    range_appends.push((open, appends.total - range_start));
                }
                XBucketRun {
                    appends,
                    range_appends,
                    stop,
                }
            },
        );
        if let Some(stop) = runs.iter().find_map(|run| run.stop) {
            return Err(stop);
        }
        if parts > 1 {
            // The runs' appends together: the serial fill's at its end.
            let total: usize = runs.iter().map(|run| run.appends.total).sum();
            let mut table_appends = vec![0usize; ranges.len()];
            for run in &runs {
                for &(index, appends) in &run.range_appends {
                    table_appends[index] += appends;
                }
            }
            let past_capacity = table_appends
                .iter()
                .zip(&ranges)
                .any(|(&appends, (_, is_table))| *is_table && appends > overflow_cell_capacity);
            if total > max_appends || past_capacity {
                return Err(XBucketStop::TooLarge);
            }
        }
        Ok(runs)
    };
    let runs = match runs(parts) {
        Ok(runs) => runs,
        // The serial fill stops at the first of the two (more than one part: the runs only show
        // that it stops).
        Err(_) if parts > 1 => match runs(1) {
            Ok(runs) => runs,
            Err(stop) => return x_bucket_stop(stop),
        },
        Err(stop) => return x_bucket_stop(stop),
    };

    // Per bucket: its first slot of `hits`, after those of the buckets before.
    let mut bucket_starts = vec![0u32; buckets + 1];
    for b in 0..buckets {
        let appends: u32 = runs.iter().map(|run| run.appends.count(b)).sum();
        bucket_starts[b + 1] = bucket_starts[b] + appends;
    }
    let total = bucket_starts[buckets] as usize;

    let kept = if pooled {
        x_flat_pool::<H>().lend()
    } else {
        XFlatTableVectors::default()
    };
    // One more entry: the number of hits.
    let mut starts = x_sized_on::<S, _>(kept.starts, backbone_size + 1, false);
    let mut hits = x_sized_on::<S, _>(kept.hits, total, pooled);
    let mut pv = x_zeroed_on(kept.pv, (backbone_size >> PV_ARRAY_BTS) + 1);

    // Ranges of whole buckets of about as many appends each, with their slices of the vectors.
    const { assert!(X_BUCKET_SHIFT as usize >= PV_ARRAY_BTS) };
    let pieces = if parts > 1 {
        (parts * X_BUCKET_RANGES_PER_PART).min(buckets)
    } else {
        1
    };
    let mut bounds = vec![0usize; pieces + 1];
    for (piece, bound) in bounds.iter_mut().enumerate().skip(1) {
        let target = (piece * total).div_ceil(pieces) as u32;
        *bound = bucket_starts
            .partition_point(|&start| start < target)
            .min(buckets);
    }
    bounds[pieces] = buckets;
    let mut slices = Vec::with_capacity(pieces);
    {
        let (mut starts_left, mut hits_left, mut pv_left) =
            (&mut starts[..backbone_size], &mut hits[..], &mut pv[..]);
        for piece in 0..pieces {
            let (b0, b1) = (bounds[piece], bounds[piece + 1].max(bounds[piece]));
            let cells = (b1 << X_BUCKET_SHIFT).min(backbone_size)
                - (b0 << X_BUCKET_SHIFT).min(backbone_size);
            let slots = (bucket_starts[b1] - bucket_starts[b0]) as usize;
            // The last range's words end with the array's (one past the last cell's).
            let words = if piece + 1 == pieces {
                pv_left.len()
            } else {
                ((b1 - b0) << (X_BUCKET_SHIFT as usize - PV_ARRAY_BTS)).min(pv_left.len())
            };
            let (s, rest) = std::mem::take(&mut starts_left).split_at_mut(cells);
            starts_left = rest;
            let (h, rest) = std::mem::take(&mut hits_left).split_at_mut(slots);
            hits_left = rest;
            let (p, rest) = std::mem::take(&mut pv_left).split_at_mut(words);
            pv_left = rest;
            slices.push(Mutex::new(Some((b0..b1, s, h, p))));
        }
    }

    // Per bucket, as `x_buckets_on`: the counts of its cells (its runs' chunks in order), their
    // starts, the PV bits and `finalize`'s statistics, the fill from the chunks, then the scan order
    // of the cells of more than 5 offsets.
    let statistics = x_bucket_items(
        parts,
        pieces,
        || (vec![0u32; 1 << X_BUCKET_SHIFT], Vec::<u32>::new()),
        |(cell_slots, overflowed), piece| {
            let (piece_buckets, starts, hits, pv) = slices[piece]
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .take()
                .expect("each bucket range taken once");
            let first_cell = piece_buckets.start << X_BUCKET_SHIFT;
            let first_slot = bucket_starts[piece_buckets.start];
            let mut nonempty_count = 0usize;
            let mut longest_chain = 0i32;
            for b in piece_buckets {
                let first = (b << X_BUCKET_SHIFT) - first_cell;
                let cells = (backbone_size - (b << X_BUCKET_SHIFT)).min(1 << X_BUCKET_SHIFT);
                // Each cell's count, then its cursor (in the range's slots).
                let slots = &mut cell_slots[..cells];
                slots.fill(0);
                for run in &runs {
                    run.appends
                        .for_each(b, |entry| slots[(entry & mask) as usize] += 1);
                }
                overflowed.clear();
                let mut slot = bucket_starts[b];
                for (cell, (start, count)) in starts[first..first + cells]
                    .iter_mut()
                    .zip(slots.iter_mut())
                    .enumerate()
                {
                    let n = *count;
                    *start = slot;
                    if n > 0 {
                        nonempty_count += 1;
                        longest_chain = longest_chain.max(n as i32);
                        compressed_pv_set(pv, first + cell, PV_ARRAY_BTS);
                        if n > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1 {
                            overflowed.push(cell as u32);
                        }
                    }
                    *count = slot - first_slot;
                    slot += n;
                }
                for run in &runs {
                    run.appends.for_each(b, |entry| {
                        let at = &mut slots[(entry & mask) as usize];
                        hits[*at as usize] = H::from_offset(entry >> 16);
                        *at += 1;
                    });
                }
                // A cursor ends at the cell's end.
                for &cell in overflowed.iter() {
                    let cell = cell as usize;
                    let begin = (starts[first + cell] - first_slot) as usize;
                    let end = slots[cell] as usize;
                    x_to_scan_order(&mut hits[begin..end]);
                }
            }
            (nonempty_count, longest_chain)
        },
    );
    drop(slices);
    starts[backbone_size] = total as u32;
    drop(runs);
    let nonempty_count: usize = statistics.iter().map(|&(count, _)| count).sum();
    let longest_chain = statistics
        .iter()
        .map(|&(_, chain)| chain)
        .max()
        .unwrap_or(0);
    let mut pv_array_bts = PV_ARRAY_BTS;
    if (nonempty_count as f64) <= 0.01 * (backbone_size as f64) {
        pv_array_bts += ilog2(backbone_size / (8 * X_TARGET_PV_BYTES));
    }

    XBucketBuild::Built(Some(XFlatBuild::Table(XFlatCompressedLookup {
        word_length: header.word_length as i32,
        compressed_alphabet_size: header.alphabet as i32,
        reciprocal_alphabet_size: header.reciprocal_alphabet_size,
        longest_chain,
        pv,
        pv_array_bts: pv_array_bts as i32,
        scaled_compress_table_u8: build_scaled_compress_table_u8(&header.scaled_compress_table),
        starts,
        hits,
    })))
}

// No NCBI counterpart: the outcome of a bucket fill that stopped (`XBucketStop`).
fn x_bucket_stop<H: XQueryOffset>(stop: XBucketStop) -> XBucketBuild<H> {
    match stop {
        XBucketStop::TooLarge => XBucketBuild::Built(Some(XFlatBuild::TooLarge)),
        XBucketStop::Declined => XBucketBuild::Declined,
    }
}

// No NCBI counterpart: the count-and-fill builder's table (where the bucket fill declines at one
// part) with the `u16` offsets of BLASTX's (the query has at most `X_U16_LETTERS` letters).
fn x_narrow_flat(table: XFlatCompressedLookup) -> XFlatCompressedLookup<u16> {
    XFlatCompressedLookup {
        word_length: table.word_length,
        compressed_alphabet_size: table.compressed_alphabet_size,
        reciprocal_alphabet_size: table.reciprocal_alphabet_size,
        longest_chain: table.longest_chain,
        pv: table.pv,
        pv_array_bts: table.pv_array_bts,
        scaled_compress_table_u8: table.scaled_compress_table_u8,
        starts: table.starts,
        hits: table
            .hits
            .iter()
            .map(|&q| {
                debug_assert!(q as u32 <= u32::from(u16::MAX));
                u16::from_offset(q as u32)
            })
            .collect(),
    }
}

impl<H: XQueryOffset> XFlatCompressedLookup<H> {
    /// The hits of cell `index`, in NCBI's scan order.
    #[inline(always)]
    pub(crate) fn cell(&self, index: usize) -> &[H] {
        &self.hits[self.starts[index] as usize..self.starts[index + 1] as usize]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:264-409
    /// `BlastCompressedAaLookupTable::scan_subject` on the flat table: the same statements, a cell's
    /// hits copied from its slice. `CALLER` gives each caller a copy of its own (`X_BLASTP`,
    /// `X_TBLASTN`, `X_BLASTX`, `X_BLASTP_SOFT`, `X_BLASTP_MATRIX`), so that the compiler inlines
    /// BLASTP's as it did when that was the only caller.
    pub(crate) fn scan_subject<const CALLER: u8>(
        &self,
        subject: &[u8],
        seq_ranges: &[(i32, i32)],
        offset_pairs: &mut [OffsetPair],
        array_size: i32,
        s_range: &mut [i32; 3],
    ) -> i32 {
        let mut totalhits = 0i32;
        let word_length = usize::try_from(self.word_length)
            .expect("NCBI BLAST compressed word length must fit usize");
        let compressed_alphabet_size = usize::try_from(self.compressed_alphabet_size)
            .expect("NCBI BLAST compressed alphabet size must fit usize");
        let recip = i64::from(self.reciprocal_alphabet_size);
        let pv_array_bts = usize::try_from(self.pv_array_bts)
            .expect("NCBI BLAST compressed PV shift must fit usize");

        while determine_compressed_scanning_offsets(
            seq_ranges,
            self.word_length,
            self.word_length,
            s_range,
        ) {
            let s_first = usize::try_from(s_range[1])
                .expect("NCBI BLAST compressed scan start must be non-negative");
            let s_last = usize::try_from(s_range[2])
                .expect("NCBI BLAST compressed scan end must be non-negative");
            let mut s = s_first;
            let mut index = 0i32;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:496-500
            // As in `scan_subject`: letters past the subject are read as NULLB.
            let mut padded_word = [0u8; 8];
            while s <= s_last {
                let end = s + word_length - 1;
                let word = match subject.get(s..end) {
                    Some(word) => word,
                    None => {
                        let available = subject.get(s..).unwrap_or(&[]);
                        let word = &mut padded_word[..word_length - 1];
                        word.fill(0);
                        word[..available.len()].copy_from_slice(available);
                        &padded_word[..word_length - 1]
                    }
                };
                let computed = compute_compressed_index_u8(
                    word_length - 1,
                    word,
                    compressed_alphabet_size,
                    &self.scaled_compress_table_u8,
                );
                index = computed.0;
                let skip = computed.1;
                if skip == 0 {
                    break;
                }
                s += 1;
            }

            let mut next_char = if s <= s_last {
                subject.get(s + word_length - 1).copied().unwrap_or(0)
            } else {
                0
            };
            let mut preshift = compressed_reciprocal_divide(index, recip);
            let scaled_compress_table = &self.scaled_compress_table_u8;

            while s <= s_last {
                let mut compressed_char = scaled_compress_table[usize::from(next_char)];
                next_char = subject.get(s + word_length).copied().unwrap_or(0);

                if compressed_char < 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:293-312
                    preshift = 0;
                    s += 1;
                    let mut skip = self.word_length - 1;
                    while skip != 0 && s <= s_last {
                        compressed_char = scaled_compress_table[usize::from(next_char)];
                        next_char = subject.get(s + word_length).copied().unwrap_or(0);
                        if compressed_char < 0 {
                            skip = self.word_length - 1;
                            preshift = 0;
                            s += 1;
                            continue;
                        }
                        index = preshift + compressed_char;
                        preshift = compressed_reciprocal_divide(index, recip);
                        skip -= 1;
                        s += 1;
                    }
                    continue;
                }

                index = preshift + compressed_char;
                preshift = compressed_reciprocal_divide(index, recip);

                let index_usize = usize::try_from(index)
                    .expect("NCBI BLAST compressed scan index must be non-negative");
                if compressed_pv_test(&self.pv, index_usize, pv_array_bts) {
                    let cell = self.cell(index_usize);
                    // `num_used` of the cell.
                    let numhits = cell.len() as i32;
                    if numhits != 0 {
                        if numhits <= array_size - totalhits {
                            let dest_base = usize::try_from(totalhits)
                                .expect("NCBI BLAST compressed total hits fit usize");
                            let s_off = s as u32;
                            for (pair, &q_off) in offset_pairs[dest_base..dest_base + cell.len()]
                                .iter_mut()
                                .zip(cell)
                            {
                                *pair = OffsetPair {
                                    q_off: q_off.offset(),
                                    s_off,
                                };
                            }
                            totalhits += numhits;
                        } else {
                            s_range[1] = i32::try_from(s)
                                .expect("NCBI BLAST compressed subject offset fits Int4");
                            return totalhits;
                        }
                    }
                }
                s += 1;
            }

            s_range[1] = i32::try_from(s).expect("NCBI BLAST compressed subject offset fits Int4");
        }

        totalhits
    }
}

impl XFlatCompressedLookup {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-503
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341
    /// The scan of `x_seed_fused_compressed::x_scan_cells` (the rounds of the two-hit scan of a
    /// subject's one range, without the offset-pair array) on the flat table: calls
    /// `visit(subject_offset, hits)` for every cell found, with its slice. `cells` is scratch.
    ///
    /// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): with `BRANCH_FREE` the offsets are listed without a branch on
    /// the presence bit, as `x_seed_fused::x_scan_cells_bf` lists them: every offset that reaches
    /// `PV_TEST` is written at the list's end, and the list's length advances by the bit. The list
    /// holds the same offsets and indices in the same order, so the calls are the same; `cells`
    /// then reserves a place for every offset of the subject.
    pub(crate) fn for_each_cell<const BRANCH_FREE: bool, V>(
        &self,
        subject: &[u8],
        subject_length: usize,
        cells: &mut Vec<(i32, u32)>,
        mut visit: V,
    ) where
        V: FnMut(i32, &[i32]),
    {
        let word_length = usize::try_from(self.word_length)
            .expect("NCBI BLAST compressed word length must fit usize");
        let compressed_alphabet_size = usize::try_from(self.compressed_alphabet_size)
            .expect("NCBI BLAST compressed alphabet size must fit usize");
        let recip = i64::from(self.reciprocal_alphabet_size);
        let pv_array_bts = usize::try_from(self.pv_array_bts)
            .expect("NCBI BLAST compressed PV shift must fit usize");
        let scaled_compress_table = &self.scaled_compress_table_u8;
        let pv = &self.pv;
        cells.clear();

        let mut s_last =
            i32::try_from(subject_length).expect("subject length must fit Int4") - self.word_length;
        if s_last < 0 {
            s_last = 0;
        }
        let s_last = s_last as usize;
        // A letter past the subject's buffer (a subject shorter than a word) is read as NULLB.
        let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0);

        let mut s = 0usize;
        let mut index = 0i32;
        let mut padded_word = [0u8; 8];
        while s <= s_last {
            let end = s + word_length - 1;
            let word = match subject.get(s..end) {
                Some(word) => word,
                None => {
                    let available = subject.get(s..).unwrap_or(&[]);
                    let word = &mut padded_word[..word_length - 1];
                    word.fill(0);
                    word[..available.len()].copy_from_slice(available);
                    &padded_word[..word_length - 1]
                }
            };
            let (computed, skip) = compute_compressed_index_u8(
                word_length - 1,
                word,
                compressed_alphabet_size,
                scaled_compress_table,
            );
            index = computed;
            if skip == 0 {
                break;
            }
            s += 1;
        }

        let mut next_char = if s <= s_last {
            byte_at(s + word_length - 1)
        } else {
            0
        };
        let mut preshift = compressed_reciprocal_divide(index, recip);

        // EXPERIMENT (LOSAT_X_BPCLUTSCANBF): a place for every offset from `s` to `s_last` (`s <=
        // s_last + 1` here), written through `list`; `listed` of them are the list.
        let list = if BRANCH_FREE {
            cells.reserve(s_last + 1 - s);
            cells.as_mut_ptr()
        } else {
            std::ptr::null_mut()
        };
        let mut listed = 0usize;
        while s <= s_last {
            let mut compressed_char = scaled_compress_table[usize::from(next_char)];
            next_char = byte_at(s + word_length);

            if compressed_char < 0 {
                preshift = 0;
                s += 1;
                let mut skip = self.word_length - 1;
                while skip != 0 && s <= s_last {
                    compressed_char = scaled_compress_table[usize::from(next_char)];
                    next_char = byte_at(s + word_length);
                    if compressed_char < 0 {
                        skip = self.word_length - 1;
                        preshift = 0;
                        s += 1;
                        continue;
                    }
                    index = preshift + compressed_char;
                    preshift = compressed_reciprocal_divide(index, recip);
                    skip -= 1;
                    s += 1;
                }
                continue;
            }

            index = preshift + compressed_char;
            preshift = compressed_reciprocal_divide(index, recip);

            // index = preshift + compressed_char >= 0: both terms are (from) non-negative letters.
            if BRANCH_FREE {
                // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:55-57 (PV_TEST)
                let ix = index as usize;
                let bit = (pv[ix >> pv_array_bts] >> (ix & PV_ARRAY_MASK)) & 1;
                // SAFETY: each offset tested adds at most one to `listed`, so with `s0` the offset
                // at the reserve, `listed <= s - s0 < s_last + 1 - s0`: inside the reserve of the
                // empty `cells`.
                unsafe { list.add(listed).write((s as i32, index as u32)) };
                listed += bit as usize;
            } else if compressed_pv_test(pv, index as usize, pv_array_bts) {
                cells.push((s as i32, index as u32));
            }
            s += 1;
        }
        if BRANCH_FREE {
            // SAFETY: the first `listed` places were written above, within the reserve.
            unsafe { cells.set_len(listed) };
        }

        let starts = &self.starts;
        for (n, &(s, index)) in cells.iter().enumerate() {
            // Cache hints only: the slice bounds of a cell X_STARTS_AHEAD cells ahead, and the hits of
            // the cell X_HITS_AHEAD ahead (whose bounds were read ahead before).
            if let Some(&(_, later)) = cells.get(n + X_STARTS_AHEAD) {
                x_prefetch(starts.as_ptr().wrapping_add(later as usize));
            }
            if let Some(&(_, later)) = cells.get(n + X_HITS_AHEAD) {
                let start = starts[later as usize] as usize;
                x_prefetch(self.hits.as_ptr().wrapping_add(start));
            }
            let cell = self.cell(index as usize);
            if !cell.is_empty() {
                visit(s, cell);
            }
        }
    }

    // No NCBI counterpart: the range of LOSAT's two-offsets-per-step scan; it does not change any
    // value NCBI computes.
    /// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): whether [`Self::x_list_pairs_lean`] lists the offsets of
    /// `subject` (`subject_length` residues): the compressed alphabet has 15 letters (NCBI's
    /// reciprocal for 15), the word length is 2 to 6 (the indices stay below 15^6, where the
    /// divisions by 15 and by 225 are exact), the presence bits are one per cell (shift
    /// `PV_ARRAY_BTS`: only this two-offset step needs one bit per cell, `for_each_cell::<true, _>`
    /// handles any shift; the tables of words of 5 letters, over a backbone of 15^5 cells, below
    /// NCBI's 8 x 262,144 for a coarser shift, always have one bit per cell, while those of 6 letters
    /// (15^6 cells) have the shift 7, so this scan declines them) and they and the
    /// slice bounds cover every index, and the subject keeps a byte
    /// past its last word (the sentinel, which a search's subjects of a word or more keep), so that
    /// the scan reads no letter out of bounds.
    pub(crate) fn x_pairs_fit(&self, subject: &[u8], subject_length: usize) -> bool {
        let Ok(word_length) = usize::try_from(self.word_length) else {
            return false;
        };
        if self.compressed_alphabet_size != 15
            || self.reciprocal_alphabet_size != 286_331_154
            || !(2..=6).contains(&word_length)
            || subject_length > subject.len()
        {
            return false;
        }
        let s_last = subject_length.saturating_sub(word_length);
        let indices = 15u64.pow(word_length as u32);
        subject.len() > s_last + word_length
            && self.starts.len() as u64 > indices
            && self.pv_array_bts == PV_ARRAY_BTS as i32
            && (self.pv.len() as u64) << PV_ARRAY_BTS >= indices
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341 (cited at
    // `for_each_cell`)
    /// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): the list of `for_each_cell::<true, _>` (each offset whose
    /// index has its presence bit, with the index, in order) for a subject where
    /// [`Self::x_pairs_fit`] holds, made two offsets per step. With both new letters `c0` and `c1`
    /// in the alphabet (scaled values `C0`, `C1`, multiples of 15^(word length - 1)) and `p` the
    /// preshift before offset `s`, the index of `s` is `i0 = p + C0`, the index of `s + 1` is
    /// `i0 / 15 + C1` (NCBI's step), and the preshift after `s + 1` is
    /// `(i0 / 15 + C1) / 15 = i0 / 225 + C1 / 15` (as `C1` is a multiple of 15), so the chain from
    /// preshift to preshift takes one multiply per two offsets; `C1 / 15` is off the chain. The
    /// divisions are NCBI's reciprocal multiply and, for 225, `x * 19,088,744 >> 32`, exact below
    /// 15^6 (lib test); `C1 / 15`, an exact division, is `C1` times the inverse of 15 modulo 2^32.
    /// The two-offset steps are a loop of their own, out of line ([`x_pairs_loop`]), so that their
    /// state (the offset, the list's length, the preshift and the next letter) stays in registers:
    /// in one loop with the other steps, the list's length was kept on the stack, loaded and stored
    /// at every pair. That loop takes NCBI's reciprocal for 15 as the constant `x_pairs_fit` asks
    /// for, computes in `u32` (every index is below 15^6) and reads the subject through a pointer a
    /// word length in. The first word, a letter outside the alphabet (NCBI's step and its skip, as
    /// in `for_each_cell`) and the last offset take one offset per step. `cells` reserves
    /// [`X_STARTS_AHEAD`] places past the list.
    ///
    /// # Safety
    ///
    /// [`Self::x_pairs_fit`]`(subject, subject_length)` holds: the scan reads the subject's bytes
    /// and the presence bits without bounds checks. `CALLER` gives each caller a copy of its own,
    /// as in [`Self::scan_subject`] (`X_BLASTP`, `X_BLASTP_MATRIX`).
    #[inline(never)]
    pub(crate) unsafe fn x_list_pairs_lean<const CALLER: u8>(
        &self,
        subject: &[u8],
        subject_length: usize,
        cells: &mut Vec<(i32, u32)>,
    ) {
        debug_assert!(self.x_pairs_fit(subject, subject_length));
        // The scaled letters are at most 14 * 15^(word length - 1), so every index and preshift
        // below stays under 15^word length.
        let word_length = self.word_length as usize;
        let recip = i64::from(self.reciprocal_alphabet_size);
        let table = &self.scaled_compress_table_u8;
        let pv = &self.pv;
        cells.clear();
        let s_last = subject_length.saturating_sub(word_length);
        // SAFETY: the loops read up to s_last + word length, below subject.len() (x_pairs_fit, the
        // caller's contract).
        let byte_at = |i: usize| unsafe { *subject.get_unchecked(i) };
        let bit_of = |index: i32| {
            let ix = index as usize;
            // SAFETY: an index is below 15^word length (the scaled letters, above), which the bits
            // cover (x_pairs_fit, the caller's contract).
            let word = unsafe { *pv.get_unchecked(ix >> PV_ARRAY_BTS) };
            (word >> (ix & PV_ARRAY_MASK)) & 1
        };

        let mut s = 0usize;
        let mut index = 0i32;
        while s <= s_last {
            let (computed, skip) = compute_compressed_index_u8(
                word_length - 1,
                &subject[s..s + word_length - 1],
                15,
                table,
            );
            index = computed;
            if skip == 0 {
                break;
            }
            s += 1;
        }
        let mut next_char = if s <= s_last {
            byte_at(s + word_length - 1)
        } else {
            0
        };
        let mut preshift = compressed_reciprocal_divide(index, recip);
        cells.reserve(s_last + 1 - s + X_STARTS_AHEAD);
        let list = cells.as_mut_ptr();
        let mut listed = 0usize;
        // SAFETY: word length < subject.len() (x_pairs_fit), so the pointer is inside the subject.
        let sub_w = unsafe { subject.as_ptr().add(word_length) };
        loop {
            let mut at = XPairsAt {
                s,
                listed,
                preshift,
                next_char,
            };
            // SAFETY: x_pairs_fit holds (the caller's contract), `next_char` is the letter at
            // s + word length - 1 and `list` has the reserve above, `listed` places of it written.
            unsafe { x_pairs_loop(sub_w, table, pv.as_ptr(), list, &mut at, s_last) };
            XPairsAt {
                s,
                listed,
                preshift,
                next_char,
            } = at;
            if s > s_last {
                break;
            }
            // Offset s_last, or a pair with a letter outside the alphabet: one offset's step.
            let mut compressed_char = table[usize::from(next_char)];
            next_char = byte_at(s + word_length);
            if compressed_char < 0 {
                preshift = 0;
                s += 1;
                let mut skip = self.word_length - 1;
                while skip != 0 && s <= s_last {
                    compressed_char = table[usize::from(next_char)];
                    next_char = byte_at(s + word_length);
                    if compressed_char < 0 {
                        skip = self.word_length - 1;
                        preshift = 0;
                        s += 1;
                        continue;
                    }
                    index = preshift + compressed_char;
                    preshift = compressed_reciprocal_divide(index, recip);
                    skip -= 1;
                    s += 1;
                }
                continue;
            }
            index = preshift + compressed_char;
            preshift = compressed_reciprocal_divide(index, recip);
            let bit = bit_of(index);
            // SAFETY: each offset adds at most one to `listed`, so `listed` stays below the
            // reserve (as in `for_each_cell`).
            unsafe { list.add(listed).write((s as i32, index as u32)) };
            listed += bit as usize;
            s += 1;
        }
        // SAFETY: the first `listed` places were written above, within the reserve.
        unsafe { cells.set_len(listed) };
    }
}

// No NCBI counterpart: the state of LOSAT's list around its two-offset steps.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): where [`x_pairs_loop`] starts and stops: the next offset `s`,
/// the list's length `listed`, the preshift before `s` and the letter at `s + word length - 1`.
struct XPairsAt {
    s: usize,
    listed: usize,
    preshift: i32,
    next_char: u8,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341 (cited at
// `XFlatCompressedLookup::for_each_cell`)
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): the two-offset steps of
/// [`XFlatCompressedLookup::x_list_pairs_lean`] from `at`, until fewer than two offsets up to
/// `s_last` are left or a pair has a letter outside the alphabet; `at` then holds that pair's first
/// offset, for the caller's step. `sub_w` points a word length into the subject, `table` holds the
/// scaled letters and `pv` the presence bits. The divisions are those of `x_list_pairs_lean`'s
/// header with NCBI's reciprocal for 15, 286,331,154, as a constant.
///
/// # Safety
///
/// `sub_w`, `table` and `pv` are those of a table and subject where
/// [`XFlatCompressedLookup::x_pairs_fit`] holds (the loop reads up to `sub_w[s_last]` and the
/// bits of indices below 15^word length without bounds checks), `at.next_char` is the letter at
/// `at.s + word length - 1`, and `list` has room for an entry per offset from `at.listed`.
#[inline(never)]
unsafe fn x_pairs_loop(
    sub_w: *const u8,
    table: &[i32; 256],
    pv: *const u32,
    list: *mut (i32, u32),
    at: &mut XPairsAt,
    s_last: usize,
) {
    let (mut s, mut listed, mut preshift, mut next_char) =
        (at.s, at.listed, at.preshift, at.next_char);
    let bit_of = |index: u32| {
        let ix = index as usize;
        // SAFETY: an index is below 15^word length, which the bits cover (this function's
        // contract).
        let word = unsafe { *pv.add(ix >> PV_ARRAY_BTS) };
        (word >> (ix & PV_ARRAY_MASK)) & 1
    };
    while s < s_last {
        // SAFETY: s + word length < s_last + word length < subject.len() (x_pairs_fit).
        let c1 = unsafe { *sub_w.add(s) };
        let scaled0 = table[usize::from(next_char)];
        let scaled1 = table[usize::from(c1)];
        if (scaled0 | scaled1) < 0 {
            break;
        }
        let i0 = (preshift + scaled0) as u32;
        let i1 = ((u64::from(i0) * 286_331_154) >> 32) as u32 + scaled1 as u32;
        // i0 / 225 + C1 / 15, as in x_list_pairs_lean's header.
        preshift = (((u64::from(i0) * 19_088_744) >> 32) as u32)
            .wrapping_add((scaled1 as u32).wrapping_mul(0xEEEE_EEEF)) as i32;
        // SAFETY: s + 1 + word length <= s_last + word length < subject.len() (x_pairs_fit).
        next_char = unsafe { *sub_w.add(s + 1) };
        let bit0 = bit_of(i0);
        // SAFETY: each offset adds at most one to `listed`, within `list`'s room (this function's
        // contract).
        unsafe { list.add(listed).write((s as i32, i0)) };
        listed += bit0 as usize;
        let bit1 = bit_of(i1);
        // SAFETY: as above.
        unsafe { list.add(listed).write((s as i32 + 1, i1)) };
        listed += bit1 as usize;
        s += 2;
    }
    *at = XPairsAt {
        s,
        listed,
        preshift,
        next_char,
    };
}

// No NCBI counterpart: how many cells ahead `for_each_cell` reads the slice bounds and the hits.
// Tuning only. blastp-fast's lean loop after the branch-free list (`x_scan_flat_lean` in
// `blastp/x_seed_fused_clut.rs`) reads ahead as far, over X_STARTS_AHEAD entries past its list, so
// the hits must not be read further ahead than the bounds.
pub(crate) const X_STARTS_AHEAD: usize = 16;
pub(crate) const X_HITS_AHEAD: usize = 8;
const _: () = assert!(X_HITS_AHEAD <= X_STARTS_AHEAD);

// No NCBI counterpart: a cache hint; it reads nothing and changes no value. `ptr` may point anywhere
// (a prefetch does not fault).
#[inline(always)]
pub(crate) fn x_prefetch<T>(ptr: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch has no effect on program state and does not fault.
    unsafe {
        std::arch::x86_64::_mm_prefetch(ptr.cast::<i8>(), std::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = ptr;
}

/// EXPERIMENT (LOSAT_X_BPCLUT): blastp-fast's compressed lookup table as the reference builds it or
/// flat.
pub(crate) enum XCompressedLookup {
    Reference(BlastCompressedAaLookupTable),
    Flat(XFlatCompressedLookup),
}

/// BLASTX's compressed lookup table: the reference table, or flat with `u16` query offsets, or flat
/// with `i32` offsets for a query batch of more than `X_U16_LETTERS` letters (several queries in
/// one batch are not split, so a batch that reaches 10002 letters with a long last query keeps it
/// whole). Matched only out of line (`blastx/seed.rs`).
pub(crate) enum XBxCompressedLookup {
    Reference(BlastCompressedAaLookupTable),
    Flat(XFlatCompressedLookup<u16>),
    Wide(XFlatCompressedLookup),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
/// EXPERIMENT (LOSAT_X_BPCLUT): `try_build_blosum62_compressed_lookup` for BLASTP. With the switch
/// off it calls it; on, it builds the flat table (the reference builder where it declines), on
/// `pool`'s threads when it has more than one (`x_clut_par.rs`). `x_flat_args` goes to the pool's
/// build (LOSAT_X_BPLKPFLAT). A build ahead of its block (`XFlatArgs::x_ahead`, LOSAT_X_BPLKPAHEAD)
/// is only the pool's build with the switch on; elsewhere it returns `None` at once.
pub(crate) fn x_bp_build_compressed_lookup(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<Result<XCompressedLookup, ()>> {
    if x_flat_args.x_ahead && (!x_bpclut_on() || !x_pool_builds(pool)) {
        return None;
    }
    let reference = || {
        try_build_blosum62_compressed_lookup(word_size, scaled_threshold, query, locations, tables)
    };
    if !x_bpclut_on() {
        return reference().map(|built| built.map(XCompressedLookup::Reference));
    }
    match x_build_flat_in_pool(
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        pool,
        x_flat_args,
    )? {
        XFlatBuild::Table(table) => Some(Ok(XCompressedLookup::Flat(table))),
        XFlatBuild::Refused => Some(Err(())),
        // LOSAT_X_BPLKPAHEAD: a build ahead of its block runs no reference builder; the block
        // builds its table itself (`XFlatArgs::x_ahead`).
        XFlatBuild::TooLarge if x_flat_args.x_ahead => None,
        XFlatBuild::TooLarge => reference().map(|built| built.map(XCompressedLookup::Reference)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
/// Item SC2d (LOSAT_X_BPLKPAHEAD): the first part of `x_bp_build_compressed_lookup` (of
/// `x_bp_build_compressed_lookup_of` for another search matrix `matrix`, item SC2d-2) for a build
/// ahead of its block (`XFlatArgs::x_ahead`): the appends pass of the pool's flat table, with the
/// extra tasks (`x_clut_par.rs`). `None` where the build would give no table ahead (the switch off, a
/// pool of one thread, `TooLarge`) or the build ends at the appends pass; the second part is
/// `x_bp_flat_ahead_table`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_bp_flat_ahead_appends(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<XFlatAppends> {
    if !x_bpclut_on() || !x_pool_builds(pool) {
        return None;
    }
    x_clut_par::x_flat_ahead_appends_in_pool(
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        pool,
        x_flat_args,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_bp_flat_ahead_appends`)
/// Item SC2d (LOSAT_X_BPLKPAHEAD): the second part of a build ahead of its block: the flat table of
/// `appends` (`x_bp_flat_ahead_appends` of `matrix`, `query` and `locations`) on the vectors that
/// `spare` holds (`XFlatSpare::x_make_room`). `None` where they lack the room (the spare keeps them)
/// or a table refuses (which the block's own build then raises).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_bp_flat_ahead_table(
    matrix: ScoringMatrix,
    appends: XFlatAppends,
    query: &[u8],
    locations: &[(i32, i32)],
    pool: &SearchPool<'_>,
    spare: &XFlatSpare,
    segments: Option<&XSegmentPool>,
) -> Option<XFlatCompressedLookup> {
    match x_clut_par::x_flat_ahead_stripes_in_pool(
        matrix, appends, query, locations, pool, spare, segments,
    )? {
        XFlatBuild::Table(table) => Some(table),
        XFlatBuild::Refused | XFlatBuild::TooLarge => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
/// TBLASTN's compressed lookup table (word sizes 5 to 7, `tblastn/search_seed.rs`): what
/// `x_bp_build_compressed_lookup` builds on one thread, the reference table with the switch off and
/// otherwise the flat one (the bucket fill, the count-and-fill builder where it declines, the
/// reference builder where that declines), with the same refusal. A function of its own, which
/// calls only the builders that blastp-fast's build calls out of line, so that blastp-fast's build
/// keeps its code.
#[inline(never)]
pub(crate) fn x_tbn_build_compressed_lookup(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
) -> Option<Result<XCompressedLookup, ()>> {
    x_build_compressed_on(
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        false,
        1,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310,1352-1354
/// `x_tbn_build_compressed_lookup` with the words of the search matrix `matrix`, any but BLOSUM62
/// (`build_compressed_score_matrix`): the same builders in the same order, the reference table with
/// the switch off. Out of line, so that BLOSUM62's build keeps its code.
#[inline(never)]
pub(crate) fn x_tbn_build_compressed_lookup_of(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    parts: usize,
) -> Option<Result<XCompressedLookup, ()>> {
    let reference = || {
        try_build_compressed_lookup(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
        )
        .map(|built| built.map(XCompressedLookup::Reference))
    };
    if !x_bpclut_on() {
        return reference();
    }
    // The bucket fill on this thread, and where it declines the count-and-fill builder's table in
    // `parts` (`x_flat_in_parts`, BLASTX's past `X_U16_LETTERS` letters; 1: on this thread).
    let bucket = if x_bpffill_on() {
        x_buckets_of(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
        )
    } else {
        XBucketBuild::Declined
    };
    let built = match bucket {
        XBucketBuild::Built(built) => built,
        XBucketBuild::Declined => x_flat_in_parts::<i32>(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            parts,
            false,
        )
        .unwrap_or_else(|| {
            x_try_build_flat_of(
                matrix,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                NCBI_OVERFLOW_CELL_CAPACITY,
                X_COUNT_CHUNK,
            )
        }),
    };
    match built? {
        XFlatBuild::Table(table) => Some(Ok(XCompressedLookup::Flat(table))),
        XFlatBuild::Refused => Some(Err(())),
        XFlatBuild::TooLarge => reference(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310,1352-1354 (see
// `x_tbn_build_compressed_lookup_of`)
/// BLASTX's compressed table of a query batch with the words of the search matrix `matrix`, any but
/// BLOSUM62 (items SC2b, SC9 of line NP): as `x_bx_build_compressed_lookup`, its bucket fill
/// (`x_bx_buckets_with`, word size 5) in `parts` and where it declines the count-and-fill builder's
/// table in `parts` (`x_flat_in_parts`; 1: on this thread), with `Int4` offsets, on fresh vectors as
/// in dev (line NP SC9e: lent by the pool of BLASTX's tables, BLOSUM45 word size 6's tables of 50
/// cDNAs x scer_pep at 8 threads took about 5 % more maximum RSS, 2,756 against 2,632-2,646 MB over
/// 20 runs each, in the same time); `x_tbn_build_compressed_lookup_of` past `X_U16_LETTERS`
/// letters. Out of line, so that BLOSUM62's build keeps its code.
#[inline(never)]
pub(crate) fn x_bx_build_compressed_lookup_of(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    parts: usize,
) -> Option<Result<XCompressedLookup, ()>> {
    if query.len() > X_U16_LETTERS || !x_bpclut_on() {
        return x_tbn_build_compressed_lookup_of(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            parts,
        );
    }
    let bucket = if x_bpffill_on() {
        x_bx_buckets_with::<_, i32, [u8]>(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            false,
            parts,
            |word| word,
        )
    } else {
        XBucketBuild::Declined
    };
    let built = match bucket {
        XBucketBuild::Built(built) => built,
        XBucketBuild::Declined => x_flat_in_parts::<i32>(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            parts,
            false,
        )
        .unwrap_or_else(|| {
            x_try_build_flat_of(
                matrix,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                NCBI_OVERFLOW_CELL_CAPACITY,
                X_COUNT_CHUNK,
            )
        }),
    };
    match built? {
        XFlatBuild::Table(table) => Some(Ok(XCompressedLookup::Flat(table))),
        XFlatBuild::Refused => Some(Err(())),
        XFlatBuild::TooLarge => try_build_compressed_lookup(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
        )
        .map(|built| built.map(XCompressedLookup::Reference)),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310,1352-1354 (see
// `x_tbn_build_compressed_lookup_of`)
/// blastp-fast's compressed lookup table with the words of the search matrix `matrix`, any but
/// BLOSUM62 (items SC2b, SC2d-2): `x_bp_build_compressed_lookup` with the matrix's words. The flat
/// table by `pool`'s threads where `x_pool_builds` (`x_try_build_flat_in_pool_of`, with
/// `x_flat_args`: the spare vectors and the build ahead, as BLOSUM62's; the reference builder where
/// it declines), else `x_tbn_build_compressed_lookup_of` on this thread; the reference table with
/// the switch off. A build ahead of its block (`XFlatArgs::x_ahead`) is only the pool's build with
/// the switch on; elsewhere it returns `None` at once. Out of line, so that BLOSUM62's build keeps
/// its code.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_bp_build_compressed_lookup_of(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<Result<XCompressedLookup, ()>> {
    if x_flat_args.x_ahead && (!x_bpclut_on() || !x_pool_builds(pool)) {
        return None;
    }
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_bpclut_on() && x_pool_builds(pool) {
        return match x_clut_par::x_try_build_flat_in_pool_of(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            pool,
            x_flat_args,
        )? {
            XFlatBuild::Table(table) => Some(Ok(XCompressedLookup::Flat(table))),
            XFlatBuild::Refused => Some(Err(())),
            // LOSAT_X_BPLKPAHEAD: a build ahead of its block runs no reference builder; the block
            // builds its table itself (`XFlatArgs::x_ahead`).
            XFlatBuild::TooLarge if x_flat_args.x_ahead => None,
            XFlatBuild::TooLarge => try_build_compressed_lookup(
                matrix,
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
            )
            .map(|built| built.map(XCompressedLookup::Reference)),
        };
    }
    let _ = (pool, x_flat_args);
    x_tbn_build_compressed_lookup_of(
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        1,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_tbn_build_compressed_lookup`)
/// BLASTX's compressed lookup table: `x_tbn_build_compressed_lookup`'s with `u16` offsets
/// (`x_bx_buckets_on`; the count-and-fill builder's table narrowed where the bucket fill declines),
/// on vectors kept between builds with `LOSAT_X_BXCLUTPOOL` (`blastx::seed::Lookup`'s drop gives
/// them back). A query of more than `X_U16_LETTERS` letters gets `x_tbn_build_compressed_lookup`'s
/// table, also on kept vectors. One instance per word length `WORD` (`-word_size` 5 to 7, chosen
/// once per query batch by `blastx::seed`), so that the bucket fill's walk of the neighbouring
/// words takes its length as a constant, as it did when 5 was the only one (item SC9 of line NP;
/// 6 and 7 decline at its PV test and take the count-and-fill builder, with the length at run
/// time). A unit under the cap (`parts` above 1, `LOSAT_X_BXBUILDPAR`) builds its table in that
/// many parts: the bucket fill in parts (`x_bx_buckets_on`, line NP SC9e) and where it declines the
/// count-and-fill builder's table in parts, the stripes writing the `u16` offsets into vectors the
/// pool lends (`x_flat_in_parts`, items SC9b (2) and SC9d of line NP).
#[inline(never)]
pub(crate) fn x_bx_build_compressed_lookup<const WORD: usize>(
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    parts: usize,
) -> Option<Result<XBxCompressedLookup, ()>> {
    let word_size = WORD;
    if query.len() > X_U16_LETTERS {
        return x_build_compressed_on(
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            x_bx_clut_pool(),
            parts,
        )
        .map(|built| {
            built.map(|table| match table {
                XCompressedLookup::Reference(table) => XBxCompressedLookup::Reference(table),
                XCompressedLookup::Flat(table) => XBxCompressedLookup::Wide(table),
            })
        });
    }
    let reference = || {
        try_build_blosum62_compressed_lookup(word_size, scaled_threshold, query, locations, tables)
            .map(|built| built.map(XBxCompressedLookup::Reference))
    };
    if !x_bpclut_on() {
        return reference();
    }
    // The bucket fill (word size 5) in `parts`, and where it declines the count-and-fill builder's
    // table in `parts` (`x_flat_in_parts`; 1: on this thread), on the pool's vectors.
    let pooled = x_bx_clut_pool();
    let bucket = if x_bpffill_on() {
        x_bx_buckets_on::<WORD>(
            scaled_threshold,
            query,
            locations,
            tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            pooled,
            parts,
        )
    } else {
        XBucketBuild::Declined
    };
    let built = match bucket {
        XBucketBuild::Built(built) => built,
        XBucketBuild::Declined => x_flat_in_parts::<u16>(
            ScoringMatrix::Blosum62,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            parts,
            pooled,
        )
        .unwrap_or_else(|| {
            x_try_build_flat_with_capacity(
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                NCBI_OVERFLOW_CELL_CAPACITY,
                X_COUNT_CHUNK,
                0,
            )
            .map(XFlatBuild::x_offsets)
        }),
    };
    match built? {
        XFlatBuild::Table(table) => Some(Ok(XBxCompressedLookup::Flat(table))),
        XFlatBuild::Refused => Some(Err(())),
        XFlatBuild::TooLarge => reference(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// NCBI builds a query chunk's table on one thread before the search threads start
// (c++/src/algo/blast/api/prelim_stage.cpp:145-188 starts them for the search).
/// Item SC9b (2) of line NP: the count-and-fill builder's table of `matrix`'s words built by `parts`
/// tasks of the pool the build runs on (`x_clut_par::x_try_build_flat_in_parts`, the same table slot
/// for slot), with offsets of type `H` written by the stripes and, with `pooled`
/// (`LOSAT_X_BXCLUTPOOL`), on vectors the pool of BLASTX's tables lends: a BLASTX unit under the cap
/// (`LOSAT_X_BXBUILDPAR`, `blastx/x_bx_mem.rs`) where the bucket fill declines (word sizes 6 and 7).
/// `None` where the builder runs on this thread: `parts` 1, or a build without threads.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn x_flat_in_parts<H: XQueryOffset + XPooledOffset>(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    parts: usize,
    pooled: bool,
) -> Option<Option<XFlatBuild<H>>> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parts > 1 {
        // LOSAT_X_BXCLUTPOOL: the table on a set of vectors the pool lends (allocated by the
        // build's caller once, not by the pool's threads for each table: freed there, glibc's
        // raised mmap threshold left them in those threads' arenas); a build that makes no table
        // drops them and gives an empty set back against the one lent.
        let lent = pooled.then(|| x_flat_pool::<H>().lend());
        let built = x_clut_par::x_try_build_flat_in_parts(
            matrix,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            parts,
            lent,
        );
        if pooled && !matches!(built, Some(XFlatBuild::Table(_))) {
            drop(x_flat_pool::<H>().take_back(XFlatTableVectors::default()));
        }
        return Some(built);
    }
    let _ = (
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        parts,
        pooled,
    );
    None
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354 (see
// `x_tbn_build_compressed_lookup`)
/// `x_tbn_build_compressed_lookup`, with `pooled` the builds on vectors kept between builds
/// (`LOSAT_X_BXCLUTPOOL`) and, at more than one part, the count-and-fill builder's table in `parts`
/// (`x_flat_in_parts`; 1: the bucket fill, or the count-and-fill builder where it declines, on this
/// thread).
fn x_build_compressed_on(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pooled: bool,
    parts: usize,
) -> Option<Result<XCompressedLookup, ()>> {
    let reference = || {
        try_build_blosum62_compressed_lookup(word_size, scaled_threshold, query, locations, tables)
            .map(|built| built.map(XCompressedLookup::Reference))
    };
    if !x_bpclut_on() {
        return reference();
    }
    // The bucket fill on this thread, and where it declines the count-and-fill builder's table in
    // `parts` on the pool's vectors (`x_flat_in_parts`; 1: on this thread).
    let bucket = if x_bpffill_on() {
        x_buckets_on(
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
            pooled,
            0,
        )
    } else {
        XBucketBuild::Declined
    };
    let built = match bucket {
        XBucketBuild::Built(built) => built,
        XBucketBuild::Declined => x_flat_in_parts::<i32>(
            ScoringMatrix::Blosum62,
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            parts,
            pooled,
        )
        .unwrap_or_else(|| {
            x_try_build_flat_with_capacity(
                word_size,
                scaled_threshold,
                query,
                locations,
                tables,
                NCBI_OVERFLOW_CELL_CAPACITY,
                X_COUNT_CHUNK,
                0,
            )
        }),
    };
    match built? {
        XFlatBuild::Table(table) => Some(Ok(XCompressedLookup::Flat(table))),
        XFlatBuild::Refused => Some(Err(())),
        XFlatBuild::TooLarge => reference(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
/// Whether the flat table is built by `pool`'s threads: a pool of more than one thread, in a build
/// with threads.
pub(crate) fn x_pool_builds(pool: &SearchPool<'_>) -> bool {
    cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )) && pool.enabled()
        && pool.threads() > 1
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
/// EXPERIMENT (LOSAT_X_BPCLUT): the flat table, by `pool`'s threads when `x_pool_builds` (the same
/// table, `x_clut_par.rs`; `x_flat_args` as LOSAT_X_BPLKPFLAT gives it), else by the serial builder
/// (which uses only `x_flat_args.appends`: it runs none of the extra tasks and takes no vectors).
pub(crate) fn x_build_flat_in_pool(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    pool: &SearchPool<'_>,
    x_flat_args: XFlatArgs<'_>,
) -> Option<XFlatBuild> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_pool_builds(pool) {
        return x_clut_par::x_try_build_flat_in_pool(
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            pool,
            x_flat_args,
        );
    }
    let _ = pool;
    x_try_build_flat_reserving(
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        x_flat_args.appends,
    )
}

// No NCBI counterpart: compares two of LOSAT's flat tables (the tests' check).
/// Panics unless the two flat tables are equal, field by field and slot by slot.
#[cfg(test)]
pub(crate) fn x_assert_same_flat_tables(got: &XFlatCompressedLookup, want: &XFlatCompressedLookup) {
    x_assert_same_flat_tables_as("flat tables:", got, want);
}

// No NCBI counterpart: compares two of LOSAT's flat tables (the tests' check).
/// `x_assert_same_flat_tables` with the messages starting with `what`.
#[cfg(test)]
pub(crate) fn x_assert_same_flat_tables_as(
    what: &str,
    got: &XFlatCompressedLookup,
    want: &XFlatCompressedLookup,
) {
    assert!(
        got.word_length == want.word_length
            && got.compressed_alphabet_size == want.compressed_alphabet_size
            && got.reciprocal_alphabet_size == want.reciprocal_alphabet_size
            && got.longest_chain == want.longest_chain
            && got.pv_array_bts == want.pv_array_bts
            && got.scaled_compress_table_u8 == want.scaled_compress_table_u8,
        "{what} table header differs"
    );
    assert!(got.pv == want.pv, "{what} pv differs");
    if let Some(index) = got
        .starts
        .iter()
        .zip(&want.starts)
        .position(|(got, want)| got != want)
        .or((got.starts.len() != want.starts.len()).then_some(0))
    {
        panic!("{what} starts differ at cell {index}");
    }
    if let Some(slot) = got
        .hits
        .iter()
        .zip(&want.hits)
        .position(|(got, want)| got != want)
        .or((got.hits.len() != want.hits.len()).then_some(0))
    {
        panic!("{what} hits differ at slot {slot}");
    }
}

/// Panics unless `flat` is `reference` (every cell's hits in the reference's scan order, the PV
/// array, the scan's parameters).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:337-348
#[cfg(test)]
pub(crate) fn x_assert_same_compressed_table(
    flat: &XFlatCompressedLookup,
    reference: &BlastCompressedAaLookupTable,
) {
    assert!(
        flat.word_length == reference.word_length
            && flat.compressed_alphabet_size == reference.compressed_alphabet_size
            && flat.reciprocal_alphabet_size == reference.reciprocal_alphabet_size
            && flat.longest_chain == reference.longest_chain
            && flat.pv_array_bts == reference.pv_array_bts
            && flat.scaled_compress_table_u8 == reference.scaled_compress_table_u8,
        "compressed table: table header differs"
    );
    assert!(flat.pv == reference.pv, "compressed table: pv differs");
    assert_eq!(
        flat.starts.len(),
        reference.backbone.len() + 1,
        "compressed table: backbone size differs"
    );
    for (index, cell) in reference.backbone.iter().enumerate() {
        let hits = flat.cell(index);
        assert!(
            hits.len() == cell.num_used as usize,
            "compressed table: cell {index} count differs ({} vs {})",
            hits.len(),
            cell.num_used
        );
        if cell.num_used > 0 {
            let order = reference.query_offsets_in_scan_order(index);
            assert!(
                hits == order.as_slice(),
                "compressed table: cell {index} differs ({hits:?} vs {order:?})"
            );
        }
    }
}

#[cfg(test)]
#[path = "x_clut_direct_tests.rs"]
mod tests;
