//! EXPERIMENT (LOSAT_X_BPCLUT; BPCLUT is a strict-set switch): the compressed lookup table of
//! blastp-fast (`CompressedAaLookupTable`) built flat, in three passes (count, layout, fill): each
//! cell's hits lie in one contiguous slice, already in the order NCBI's compressed scan copies them,
//! and the scans read a cell as that slice.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849
//! ```c
//!     CompressedLookupBackboneCell *backbone_cell = lookup->backbone + index;
//!     Int4 num_entries = backbone_cell->num_used;
//!
//!     switch (num_entries) {
//!     case 0:
//!     	backbone_cell->query_offset = query_offset;
//!     	break;
//!     case 1:
//!     case 2:
//!     case 3:
//!     case 4:
//!         backbone_cell->payload.query_offsets[num_entries-1] = query_offset;
//!         break;
//!     case 5:
//!     {
//!         Int4 tmp_offset_0 = backbone_cell->payload.query_offsets[0];
//!         Int4 tmp_offset_1 = backbone_cell->payload.query_offsets[1];
//!         ...
//!         CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
//!         new_cell->next = NULL;
//!         new_cell->query_offsets[0] = backbone_cell->payload.query_offsets[2];
//!         new_cell->query_offsets[1] = backbone_cell->payload.query_offsets[3];
//!         new_cell->query_offsets[2] = query_offset;
//!         backbone_cell->payload.overflow_list.query_offsets[0] = tmp_offset_0;
//!         backbone_cell->payload.overflow_list.query_offsets[1] = tmp_offset_1;
//!         backbone_cell->payload.overflow_list.head = new_cell;
//!     }
//!     break;
//!     default:
//!     {
//!         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
//!         if (cell_index == 0 ) { /* can't be empty => it's full  */
//!         	CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
//!             new_cell->next = backbone_cell->payload.overflow_list.head;
//!             backbone_cell->payload.overflow_list.head = new_cell;
//!         }
//!         backbone_cell->payload.overflow_list.head->query_offsets[cell_index] = query_offset;
//!     }
//!     break;
//!     }
//!     backbone_cell->num_used++;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:348-396
//! ```c
//!              dest[0].qs_offsets.q_off = backbone_cell->query_offset;
//!              ...
//!              if (numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL+1) {
//!             	 dest ++;
//!                 for (i = 0; i < numhits-1; i++) {
//!                    dest[i].qs_offsets.q_off = backbone_cell->payload.query_offsets[i];
//!                 ...
//!              else {
//!                 CompressedOverflowCell* curr_cell =
//!                                     backbone_cell->payload.overflow_list.head;
//!                 Int4 first_cell_entries = (numhits - 3) & COMPRESSED_HITS_CELL_MASK;
//!                 dest[1].qs_offsets.q_off = backbone_cell->payload.overflow_list.query_offsets[0];
//!                 dest[2].qs_offsets.q_off = backbone_cell->payload.overflow_list.query_offsets[1];
//!                 dest += 3;
//!                 for (i = 0; i < first_cell_entries; i++) {
//!                    dest[i].qs_offsets.q_off =  curr_cell->query_offsets[i];
//!                 ...
//!                 if (first_cell_entries) {
//!                    curr_cell = curr_cell->next;
//!                 }
//!                    while (curr_cell != NULL) {
//!                       dest += i;
//!                       for (i = 0; i < COMPRESSED_HITS_PER_OVERFLOW_CELL; i++) {
//!                          dest[i].qs_offsets.q_off =  curr_cell->query_offsets[i];
//!                       ...
//!                       curr_cell = curr_cell->next;    /* prefetch */
//! ```
//!
//! The reference (`try_build_blosum62_compressed_lookup`) runs `add_word_hits` for every query offset
//! of every lookup segment, in order: each call appends the offset to the cells of its words (the word
//! itself, then its neighbors in the order of the `add_word_hits_core` recursion) with `add_word_hit`,
//! which keeps a cell's first offsets in the backbone and the rest in a list of overflow cells. Here:
//! - count pass: the same appends in the same order (`XWords::for_each_word`: the statements of
//!   `add_word_hits` and `add_word_hits_core`, the index of a neighbor formed as its letters are
//!   chosen, `sum c[i] * A^i` as `compute_encoded_compressed_index` forms it), recorded as the cell of
//!   each append and the number of appends of each query offset; then each cell's count goes up by
//!   one per recorded append, in order.
//! - layout: the cells' slices one after the other in index order (`starts`, the prefix sums of the
//!   counts), and the PV array and `longest_chain` from the counts by `finalize`'s statements.
//! - fill: the recorded appends again in order, each written at its cell's cursor (the cell's start,
//!   moved on by one per append), so that a cell's slice holds its offsets in the order of their
//!   appends; then the offset of rank k (0-based) of a cell whose count n is above 5 is moved to slot
//!   `x_emission_position(k, n)` (for n <= 5 that slot is k).
//!
//! The permutation. Let a cell take the offsets o_0, ..., o_{n-1} in this order (o_k is appended when
//! `num_used` is k). By `add_word_hit` (above): o_0 is `query_offset` and o_1..o_4 are
//! `payload.query_offsets[0..4]`; when n >= 6, o_5 makes the first overflow cell, which takes o_3, o_4
//! and o_5 in slots 0, 1, 2 while o_1, o_2 stay as `overflow_list.query_offsets[0..2]`; each o_k with
//! k >= 6 goes to slot (k - 3) & 3 of the head, and a new head is put in front when that slot is 0.
//! So for n >= 6 the offset o_k with k >= 3 lies in the g-th overflow cell of the cell, g = (k - 3) >> 2,
//! at slot s = (k - 3) & 3, the list runs from the newest cell (g = G - 1) to the first (g = 0, whose
//! `next` is NULL), G = ceil((n - 3) / 4) = floor(n / 4) cells, and the newest holds
//! m = n - 3 - 4 (G - 1) offsets, 1 <= m <= 4. The scan (above) writes o_0 at dest[0]; for n <= 5 it
//! writes o_k at dest[k]. For n >= 6 it writes o_1, o_2 at dest[1], dest[2], then
//! f = (n - 3) & 3 slots of the head; n - 3 = m (mod 4), so f = m when f != 0, and the walk then moves
//! to the next cell, while for f = 0 (m = 4) the walk starts at the head: either way the head's m
//! slots come first, then cells G - 2, ..., 0 with 4 slots each. Hence the emission position of o_k:
//! - k when n <= 5 or k < 3;
//! - 3 + s when g = G - 1 (the head);
//! - 3 + m + 4 (G - 2 - g) + s = n - 4 (g + 1) + s when g < G - 1.
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
//! The refusal. The reference tests NCBI's tables (`tables`, ranges of the segments) while it adds the
//! words: after each query offset of a table, while the table's appends have taken no more cells of the
//! table of all batches than NCBI's banks hold, the table passes (a batch's n offsets of a cell take
//! cells(c + n) - cells(c) >= cells(n) of them, cells(n) = 0 for n <= 5 and floor(n / 4) from 6 on);
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
//! The table is the count-and-fill builder's slot for slot; only the order of the memory accesses
//! differs.
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
use std::sync::{Mutex, PoisonError};

// EXPERIMENT (LOSAT_X_BPCLUT at more than one thread): the same table built by the search's pool
// (exactness argument in the file's header).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_clut_par.rs"]
mod x_clut_par;
/// PROBE (v4-BPF s11): the table-set pool's entry points for the block search.
#[cfg(feature = "parallel")]
pub(crate) use x_clut_par::{x_s11_give_table, x_s11_t};
/// PROBE (v4-BPF s11): no pool without threads.
#[cfg(not(feature = "parallel"))]
pub(crate) fn x_s11_t() -> Option<usize> {
    None
}
/// PROBE (v4-BPF s11): no pool without threads.
#[cfg(not(feature = "parallel"))]
pub(crate) fn x_s11_give_table(_starts: Vec<u32>, _hits: Vec<i32>) {}

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
// ```c
//     lookup->backbone = (CompressedLookupBackboneCell* )calloc(
//                                     lookup->backbone_size,
//                                     sizeof(CompressedLookupBackboneCell));
// ...
// BlastCompressedAaLookupTable *BlastCompressedAaLookupTableDestruct(
//                                  BlastCompressedAaLookupTable * lookup)
// {
//     ...
//     sfree(lookup->backbone);
// ```
// NCBI allocates a zeroed table for each search and frees it with the search; a block search makes
// a table per block. Buffer reuse only: the parallel build writes every entry of the vectors it keeps
// before any entry is read (`x_clut_par.rs`), so the values left in them do not matter.
/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the `starts` and `hits` of a flat table that its search no longer
/// reads, kept by the block search (`blastp/x_query_blocks.rs` `XBlockSubjects`) for the next
/// parallel build: the table's vectors are allocated by the caller of the builds, once, and then
/// recycled instead of being zeroed again for each block. The second field (LOSAT_X_BPLKPAHEAD) holds
/// a PV array's vector for a build ahead of its block (`x_make_room`); it is empty otherwise.
#[derive(Default)]
pub(crate) struct XFlatSpare(Mutex<(Vec<u32>, Vec<i32>)>, Mutex<Vec<u32>>);

impl XFlatSpare {
    /// PROBE (v4-BPF s11): holds `starts` and `hits` (a set of the pool) for the build that takes
    /// this spare.
    pub(crate) fn x_s11_put(&self, starts: Vec<u32>, hits: Vec<i32>) {
        *self.0.lock().unwrap_or_else(PoisonError::into_inner) = (starts, hits);
    }

    // No NCBI counterpart: buffer reuse only (above).
    /// Keeps `table`'s `starts` and `hits` (in place of any vectors kept before); the rest of the
    /// table is dropped.
    pub(crate) fn give(&self, table: XFlatCompressedLookup) {
        let old = std::mem::replace(
            &mut *self.0.lock().unwrap_or_else(PoisonError::into_inner),
            (table.starts, table.hits),
        );
        // PROBE (v4-BPF s8, LOSAT_BPF_S8STASH=1): the vectors held before are kept for later builds.
        if x_clut_par::x_s8_stash() && (old.0.capacity() > 0 || old.1.capacity() > 0) {
            let mut stash = x_clut_par::X_S8_TABLES
                .lock()
                .unwrap_or_else(PoisonError::into_inner);
            if stash.len() < 4 {
                stash.push(old);
            }
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
    /// the block before gave back), and each one short of its room is replaced here by an empty one
    /// of `backbone_size + 1`, `hits_grown` (at least `hits`) and `(backbone_size >> PV_ARRAY_BTS) +
    /// 1` entries of capacity. Nothing is written into them. Returns the number of vectors allocated.
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
            kept.1 = Vec::with_capacity(hits_grown.max(hits));
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
    pub(crate) spare: Option<&'a XFlatSpare>,
    pub(crate) extra: Option<&'a dyn XExtraTasks>,
    /// LOSAT_X_BPLKPAHEAD: the build of a block's table ahead of that block, on a pool task
    /// (`blastp/blast_engine/x_lookup_ahead.rs`): it allocates none of the table's vectors, takes
    /// them all from `spare` (made ready by `XFlatSpare::x_make_room`), and makes no table where they
    /// lack the room, the appends pass the bound (`TooLarge`) or the switch's mode is `Off`; the
    /// block then builds its table itself.
    pub(crate) x_ahead: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1189
// ```c
//     const Int4 kTargetPVBytes = 262144;
// ```
const X_TARGET_PV_BYTES: usize = 262_144;

/// EXPERIMENT (LOSAT_X_BPCLUT): the compressed lookup table with each cell's hits in one slice of
/// `hits` (`starts[index]..starts[index + 1]`), in NCBI's scan order. The other fields are those of
/// `BlastCompressedAaLookupTable` that the scans read.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:236-267
// ```c
// typedef struct BlastCompressedAaLookupTable {
//     Int4 threshold;
//     Int4 word_length;
//     Int4 alphabet_size;
//     Int4 compressed_alphabet_size;
//     Int4 reciprocal_alphabet_size;
//     ...
//     PV_ARRAY_TYPE *pv;
//     Int4 pv_array_bts;
//     Uint1* compress_table;
//     Int4* scaled_compress_table;
// } BlastCompressedAaLookupTable;
// ```
pub(crate) struct XFlatCompressedLookup {
    pub(crate) word_length: i32,
    pub(crate) compressed_alphabet_size: i32,
    pub(crate) reciprocal_alphabet_size: i32,
    pub(crate) longest_chain: i32,
    /// The cells with hits (`finalize`'s `count`). EXPERIMENT (LOSAT_X_BPCLUTSCANBF): chooses how
    /// `for_each_cell` lists the offsets ([`XFlatCompressedLookup::x_scan_bf_pays`]).
    pub(crate) used_cells: usize,
    pub(crate) pv: Vec<u32>,
    pub(crate) pv_array_bts: i32,
    pub(crate) scaled_compress_table_u8: [i32; 256],
    /// `backbone_size + 1` entries: cell `index`'s hits are `hits[starts[index]..starts[index + 1]]`.
    pub(crate) starts: Vec<u32>,
    pub(crate) hits: Vec<i32>,
}

/// What the flat builder made.
pub(crate) enum XFlatBuild {
    Table(XFlatCompressedLookup),
    /// One of NCBI's tables fills its overflow banks (the reference's `Err(())`).
    Refused,
    /// More appends than the slices' `u32` offsets take: the reference builder runs instead.
    TooLarge,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 and
// c++/src/algo/blast/core/blast_aascan.c:348-396 (module header)
// ```c
//         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
//         ...
//                 Int4 first_cell_entries = (numhits - 3) & COMPRESSED_HITS_CELL_MASK;
// ```
/// The position in NCBI's scan order of the offset appended with rank `k` (0-based) to a cell that
/// ends with `n` offsets (module header).
#[inline(always)]
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:800-806,825-836
// ```c
//     case 5:
//     {
//         ...
//         CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
//     ...
//     default:
//     {
//         ...
//         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
//
//         if (cell_index == 0 ) { /* can't be empty => it's full  */
//
//             /* fetch next vacant cell */
//         	CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
// ```
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
}

impl<'a> XWords<'a> {
    fn new(
        header: &'a XHeader,
        threshold: i32,
        matrix: &'a CompressedScoreMatrix,
        info: &'a CompressedNeighborInfo,
    ) -> Self {
        let mut power = [0u32; 8];
        let mut scale = 1u32;
        for p in power.iter_mut().take(header.word_length) {
            *p = scale;
            scale *= header.alphabet as u32;
        }
        let mut last = [[0u32; BLASTAA_SIZE]; BLASTAA_SIZE];
        for (c, row) in last.iter_mut().enumerate() {
            for (i, term) in row.iter_mut().enumerate().take(header.alphabet) {
                *term = u32::from(info.matrix_sorted_char[c][i]) * power[header.word_length - 1];
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
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1075-1121
    // ```c
    //     score = 0;
    //     for (i = 0; i < lookup->word_length; i++) {
    //         int c = lookup->compress_table[w[i]];
    //
    //         if (c >= lookup->compressed_alphabet_size) /* "non-20 aa": skip it*/
    //             return;
    //
    //         score += info->matrix[w[i]][c];
    //     }
    //     ...
    //     if (lookup->threshold == 0 || score < lookup->threshold) {
    //         s_CompressedLookupAddUnencoded(lookup, w, query_offset);
    //     }
    //     ...
    //     if (lookup->threshold == 0)
    //         return;
    //     ...
    //     score = info->row_max[w[0]];
    //     for (i = 1; i < lookup->word_length; i++)
    //         score += info->row_max[w[i]];
    //
    //     s_CompressedAddWordHitsCore(info, score, 0);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:919-929
    // ```c
    //     index = s_ComputeCompressedIndex(lookup->word_length, w,
    //                                      lookup->compressed_alphabet_size,
    //                                      &skip, lookup);
    //     if (skip == 0)
    //         s_CompressedLookupAddWordHit(lookup, index, query_offset);
    // ```
    /// Gives `sink` every cell that `add_word_hits` appends the offset of `word` to, in its order.
    #[inline]
    fn for_each_word(&self, word: &[u8], sink: &mut impl XSink) {
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
                self.word_length,
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
        self.neighbors(word, score, 0, 0, sink);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053
    // ```c
    //     score -= info->row_max[currQueryChar];
    //     rowSorted = info->matrixSorted[currQueryChar];
    //     charSorted = info->matrixSortedChar[currQueryChar];
    //
    //     if (current_pos == info->wordsize - 1) {
    //         ...
    //         for (i = 0; i < info->compressed_alphabet_size &&
    //                (score + rowSorted[i] >= info->threshold); i++) {
    //             subject_word[current_pos] = charSorted[i];
    //             s_CompressedLookupAddEncoded(lookup, subject_word,
    //                                          query_offset);
    //         ...
    //         return;
    //     }
    //     ...
    //     for (i = 0; i < info->compressed_alphabet_size &&
    //            (score + rowSorted[i] >= info->threshold); i++) {
    //         subject_word[current_pos] = charSorted[i];
    //         s_CompressedAddWordHitsCore(info, score + rowSorted[i],
    //                                     current_pos + 1);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:894-896
    // ```c
    //     case 5:
    //         index = w[0] + W6p1[w[1]] + W6p2[w[2]] +
    //                 W6p3[w[3]] + W6p4[w[4]];
    //         break;
    // ```
    // `partial` is the index of the letters chosen so far, sum c[i] * A^i over i < current_pos (the
    // terms of `compute_encoded_compressed_index`). At the last position the loop's letters are the
    // prefix of the sorted row up to the first score below the threshold (`partition_point`; the row
    // does not increase, checked in `new`).
    fn neighbors(
        &self,
        query_word: &[u8],
        mut score: i32,
        current_pos: usize,
        partial: u32,
        sink: &mut impl XSink,
    ) {
        let curr_query_char = query_word[current_pos] as usize;
        score -= self.info.row_max[curr_query_char];
        let row_sorted = &self.info.matrix_sorted[curr_query_char][..self.alphabet];
        let threshold = self.threshold;

        if current_pos == self.word_length - 1 {
            let count =
                row_sorted.partition_point(|&letter_score| score + letter_score >= threshold);
            sink.many(partial, &self.last[curr_query_char][..count]);
            return;
        }

        let char_sorted = &self.info.matrix_sorted_char[curr_query_char];
        let scale = self.power[current_pos];
        if current_pos + 2 == self.word_length {
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
                query_word,
                score + letter_score,
                current_pos + 1,
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
    // ```c
    //     if (word_size == 6 || word_size == 5) {
    //         lookup->compressed_alphabet_size = 15;
    //         lookup->reciprocal_alphabet_size = 286331154;
    //     }
    //     else {
    //         lookup->compressed_alphabet_size = 10;
    //         lookup->reciprocal_alphabet_size = 429496730;
    //     }
    //
    //     lookup->backbone_size = (Int4)pow(lookup->compressed_alphabet_size,
    //                                       word_size) + 1;
    // ```
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
// ```c
//     for (loc = location; loc; loc = loc->next){
//         Int4 from = loc->ssr->left;
//         Int4 to = loc->ssr->right - lookup->word_length + 1;
//
//         for (offset = from; offset <= to; offset++){
//             s_CompressedAddWordHits(&info, query->sequence, offset);
//         }
//     }
// ```
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
// ```c
//     if (lookup->curr_overflow_cell ==
//                         COMPRESSED_OVERFLOW_CELLS_IN_BANK) {
//         /* need a new bank */
//         Int4 bank_idx = lookup->curr_overflow_bank + 1;
//         lookup->overflow_banks[bank_idx] = (CompressedOverflowCell*) malloc(
//                                            COMPRESSED_OVERFLOW_CELLS_IN_BANK *
//                                            sizeof(CompressedOverflowCell));
//         ASSERT(bank_idx < COMPRESSED_OVERFLOW_MAX_BANKS);
// ```
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

// No NCBI counterpart: the most appends the flat table takes; past it the reference builder runs. A
// query offset appends to each cell at most twice (its word and a neighbor), so the counts and the
// slices' `u32` offsets stay below `i32::MAX` (`num_used` is an Int4).
fn x_max_appends(backbone_size: usize) -> usize {
    (i32::MAX as usize).saturating_sub(2 * backbone_size)
}

/// EXPERIMENT (LOSAT_X_BPCLUT): `try_build_blosum62_compressed_lookup` as a flat table. With
/// LOSAT_X_BPFFILL the bucket fill builds it (the count-and-fill builder where it declines).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
pub(crate) fn x_try_build_flat_compressed_lookup(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
) -> Option<XFlatBuild> {
    if x_bpffill_on() {
        match x_try_build_flat_by_buckets(
            word_size,
            scaled_threshold,
            query,
            locations,
            tables,
            NCBI_OVERFLOW_CELL_CAPACITY,
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
    )
}

// No NCBI counterpart: how many appends the count pass records before it counts them (and tests NCBI's
// table, as the reference does after each query offset); it bounds the record kept before a refusal.
const X_COUNT_CHUNK: usize = 1 << 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop over the
// segments, see `x_for_each_offset`)
// ```c
//     for (loc = location; loc; loc = loc->next){
// ```
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
    // ```c
    //     Int4 num_entries = backbone_cell->num_used;
    //     ...
    //     backbone_cell->num_used++;
    // ```
    /// Counts the appends recorded since the last call, in order.
    fn count_new(&mut self) {
        for &index in &self.targets[self.counted..] {
            let c = self.counts[index as usize];
            self.cells += usize::from(x_takes_overflow_cell(c));
            self.counts[index as usize] = c + 1;
        }
        self.counted = self.targets.len();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop, see
// `x_for_each_offset`), :787-849 (the appends, module header), :1182-1242 (`finalize`, below)
// ```c
//         for (offset = from; offset <= to; offset++){
//             s_CompressedAddWordHits(&info, query->sequence, offset);
//         }
// ```
// The passes of the module header; the tests of NCBI's tables are those of
// `try_build_blosum62_compressed_lookup_with_capacity`, made after each chunk of appends and at the end
// of a table instead of after each query offset (the reference tests a table once, at the first query
// offset where its appends have taken more cells of the table of all batches than the capacity, and its
// outcome, the batch's own table, does not depend on that offset; the cells taken only grow, so a table
// that passes that bound passes it at the end of the table).
fn x_try_build_flat_with_capacity(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
    count_chunk: usize,
) -> Option<XFlatBuild> {
    let header = XHeader::new(word_size)?;
    let matrix = build_blosum62_compressed_score_matrix(header.alphabet)?;
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let backbone_size = header.backbone_size;
    let max_appends = x_max_appends(backbone_size);

    // Count pass, over the segments before, inside and between NCBI's tables (the words between
    // tables, a split query batch, are added without a test).
    let ranges = x_count_ranges(tables, locations.len());
    let mut st = XCount {
        targets: Vec::new(),
        per_offset: Vec::new(),
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
            refused = !(words_before
                && !x_table_fills_banks(
                    &words,
                    backbone_size,
                    query,
                    range_locations,
                    overflow_cell_capacity,
                ));
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
    //             longest_chain = MAX(count, longest_chain);
    // ```
    // The statements of `finalize` on the counts (`num_used`); then each count becomes its cell's start
    // (the sum of the counts before it), the last entry the number of hits.
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
    for (index, entry) in starts[..backbone_size].iter_mut().enumerate() {
        let count = *entry;
        *entry = total;
        if count > 0 {
            compressed_pv_set(&mut pv, index, pv_array_bts);
            longest_chain = longest_chain.max(count as i32);
            if count > COMPRESSED_HITS_PER_BACKBONE_CELL as u32 + 1 {
                overflowed.push(index as u32);
            }
        }
        total += count;
    }
    starts[backbone_size] = total;

    // Fill pass: the appends again in order, each at the next slot of its cell (the starts serve as the
    // cells' cursors), so each cell holds o_0..o_{n-1}; a cursor ends at the next cell's start, so the
    // starts are then those of the cells after them, moved back by one entry. Then each cell of more
    // than 5 offsets puts o_k at `x_emission_position(k, n)` (cells of 5 or fewer: position k).
    let mut hits = vec![0i32; total as usize];
    let mut at = 0usize;
    for &(query_offset, appends) in &per_offset {
        let end = at + appends as usize;
        for &index in &targets[at..end] {
            let cursor = &mut starts[index as usize];
            hits[*cursor as usize] = query_offset;
            *cursor += 1;
        }
        at = end;
    }
    drop(targets);
    drop(per_offset);
    // The last cell (index `backbone_size - 1` = A^w) is never a word's index, so its cursor stayed at
    // `total`, the last entry.
    starts.copy_within(0..backbone_size, 1);
    starts[0] = 0;
    let mut appended: Vec<i32> = Vec::new();
    for &index in &overflowed {
        let index = index as usize;
        let cell = &mut hits[starts[index] as usize..starts[index + 1] as usize];
        let n = cell.len() as u32;
        appended.clear();
        appended.extend_from_slice(cell);
        for (k, &query_offset) in appended.iter().enumerate() {
            cell[x_emission_position(k as u32, n) as usize] = query_offset;
        }
    }

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
enum XBucketBuild {
    /// `x_try_build_flat_with_capacity`'s outcome: its table, its `TooLarge`, or `None` where the
    /// table's header cannot be made.
    Built(Option<XFlatBuild>),
    /// The count-and-fill builder decides (`x_try_build_flat_by_buckets`).
    Declined,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175 (the loop, see
// `x_for_each_offset`), :787-849 (the appends, module header), :1182-1242 (`finalize`, below)
// ```c
//         for (offset = from; offset <= to; offset++){
//             s_CompressedAddWordHits(&info, query->sequence, offset);
//         }
// ```
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
    let Some(matrix) = build_blosum62_compressed_score_matrix(header.alphabet) else {
        return XBucketBuild::Built(None);
    };
    let info = CompressedNeighborInfo::new(&matrix, header.alphabet);
    let words = XWords::new(&header, scaled_threshold, &matrix, &info);
    let max_appends = x_max_appends(backbone_size);

    // The vectors in the count-and-fill builder's order: the starts (its counts) first, the hits
    // after the appends; the record of the appends is dropped at the end, the offsets' appends
    // after the partition (they bound the list of a cell's offsets moved to the scan order: a query
    // offset appends to a cell at most twice).
    // One more entry: the number of hits.
    let mut starts = vec![0u32; backbone_size + 1];

    // Appends, over the segments before, inside and between NCBI's tables, recorded as the
    // count-and-fill builder records them; not counted.
    let ranges = x_count_ranges(tables, locations.len());
    let mut targets: Vec<u32> = Vec::new();
    let mut per_offset: Vec<(i32, u32)> = Vec::new();
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
        if too_large {
            return XBucketBuild::Built(Some(XFlatBuild::TooLarge));
        }
        if declined {
            return XBucketBuild::Declined;
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
    let mut pv = vec![0u32; (backbone_size >> PV_ARRAY_BTS) + 1];
    let mut hits = vec![0i32; total];
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
    //             longest_chain = MAX(count, longest_chain);
    // ```
    // `pv_array_bts` is `PV_ARRAY_BTS` whatever `count` (the larger shift is the same, tested
    // above), so the bits are set as the counts are read.
    let mut cell_slots = vec![0u32; 1 << X_BUCKET_SHIFT];
    let mut overflowed: Vec<u32> = Vec::new();
    let mut appended: Vec<i32> = Vec::new();
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
            let n = (end - begin) as u32;
            let slice = &mut hits[begin..end];
            appended.clear();
            appended.extend_from_slice(slice);
            for (k, &query_offset) in appended.iter().enumerate() {
                slice[x_emission_position(k as u32, n) as usize] = query_offset;
            }
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
        used_cells: nonempty_count,
        pv,
        pv_array_bts: pv_array_bts as i32,
        scaled_compress_table_u8: build_scaled_compress_table_u8(&header.scaled_compress_table),
        starts,
        hits,
    })))
}

impl XFlatCompressedLookup {
    /// The hits of cell `index`, in NCBI's scan order.
    #[inline(always)]
    pub(crate) fn cell(&self, index: usize) -> &[i32] {
        &self.hits[self.starts[index] as usize..self.starts[index + 1] as usize]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:264-409
    // ```c
    // while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
    //     for(s = s_first; s <= s_last; s++){
    //         index = s_ComputeCompressedIndex(word_length - 1, s, ...);
    //         if(!skip) break;
    //     }
    //     next_char = ((s <= s_last)? s[word_length-1] : 0);
    //     preshift = (Int4)((((Int8)index) * recip) >> 32);
    //     for (; s <= s_last; s++) {
    //         compressed_char = scaled_compress_table[next_char];
    //         next_char = s[word_length];
    //         ...
    //         if (PV_TEST(pv, index, pv_array_bts)) {
    //            ...
    //            numhits = backbone_cell->num_used;
    //            if (numhits) {
    //            if (numhits <= (array_size - totalhits)) {
    //               /* copy the hits to the destination */
    //               ...
    //            else
    //            {
    //               s_range[1] = (Int4)(s - subject->sequence);
    //               return totalhits;
    //            }
    //            }
    //         }
    //     }
    //     s_range[1] = (Int4)(s - subject->sequence);
    // }
    // ```
    /// `BlastCompressedAaLookupTable::scan_subject` on the flat table: the same statements, a cell's
    /// hits copied from its slice.
    pub(crate) fn scan_subject(
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
            // ```c
            //     scan_range[1] = subject->seq_ranges[0].left;
            //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
            //
            //     if (scan_range[2] < scan_range[1])
            //         scan_range[2] = scan_range[1];
            // ```
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
                    // ```c
                    //        if(compressed_char < 0){ /* flush (rare) "bad" character(s) */
                    //          preshift = 0;
                    //          s++;
                    //          for(skip = word_length-1; skip && (s <= s_last) ; s++){
                    //            ...
                    //          }
                    //          s--; /*undo the following increment*/
                    //          continue;
                    //        }
                    // ```
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
                                    q_off: q_off as u32,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-503
    // ```c
    //     scan_range[0] = 0;
    //     scan_range[1] = subject->seq_ranges[0].left;
    //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
    //
    //     if (scan_range[2] < scan_range[1])
    //         scan_range[2] = scan_range[1];
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341
    // ```c
    //     /* prime the index */
    //     for(s = s_first; s <= s_last; s++){
    //         index = s_ComputeCompressedIndex(word_length - 1, s,
    //                                          compressed_alphabet_size,
    //                                          &skip, lookup);
    //         if(!skip)
    //           break;
    //     }
    //     ...
    //        index = preshift + compressed_char;
    //        preshift = (Int4)((((Int8)( index )) * recip) >> 32);
    //
    //        /* if there are hits */
    //        if (PV_TEST(pv, index, pv_array_bts)) {
    //           ...
    //           numhits = backbone_cell->num_used;
    // ```
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
                // NCBI reference: c++/include/algo/blast/core/blast_lookup.h:55-57 (PV_TEST)
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
    /// EXPERIMENT (LOSAT_X_BPFSTAGE): whether [`Self::x_list_pairs`] lists the offsets of `subject`
    /// (`subject_length` residues): the compressed alphabet has 15 letters (NCBI's reciprocal for
    /// 15), the word length is 2 to 6 (the indices stay below 15^6, where the divisions by 15 and
    /// by 225 are exact), the presence bits are one per cell (shift `PV_ARRAY_BTS`, as on every
    /// table with 1 % of its cells used or more, which the branch-free list needs) and they and the
    /// slice bounds cover every index, and the subject keeps a byte past its last word (the
    /// sentinel, which a search's subjects of a word or more keep), so that the scan reads no letter
    /// out of bounds.
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341 (quoted at
    // `for_each_cell`)
    /// EXPERIMENT (LOSAT_X_BPFSTAGE): the list of `for_each_cell::<true, _>` (each offset whose
    /// index has its presence bit, with the index, in order) for a subject where
    /// [`Self::x_pairs_fit`] holds, made two offsets per step and out of line (the scan inlined in
    /// a larger loop kept its state on the stack). With both new letters `c0` and `c1` in the
    /// alphabet (scaled values `C0`, `C1`, multiples of 15^(word length - 1)) and `p` the preshift
    /// before offset `s`, the index of `s` is `i0 = p + C0`, the index of `s + 1` is `i0 / 15 + C1`
    /// (NCBI's step), and the preshift after `s + 1` is `(i0 / 15 + C1) / 15 = i0 / 225 + C1 / 15`
    /// (as `C1` is a multiple of 15), so the chain from preshift to preshift takes one multiply per
    /// two offsets; `C1 / 15` is off the chain. The divisions are NCBI's reciprocal multiply and,
    /// for 225, `x * 19,088,744 >> 32`, exact below 15^6 (lib test); `C1 / 15`, an exact division,
    /// is `C1` times the inverse of 15 modulo 2^32. A letter outside the alphabet
    /// takes NCBI's step and its skip, as `for_each_cell` does. `cells` reserves
    /// [`X_STARTS_AHEAD`] places past the list.
    ///
    /// # Safety
    ///
    /// [`Self::x_pairs_fit`]`(subject, subject_length)` holds: the scan reads the subject's bytes
    /// and the presence bits without bounds checks.
    #[inline(never)]
    pub(crate) unsafe fn x_list_pairs(
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
        let div225 = |x: i32| ((i64::from(x) * 19_088_744) >> 32) as i32;
        // 15 * 0xEEEE_EEEF = 1 modulo 2^32.
        let div15_exact = |x: i32| (x as u32).wrapping_mul(0xEEEE_EEEF) as i32;

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
        while s <= s_last {
            if s < s_last {
                let c1 = byte_at(s + word_length);
                let scaled0 = table[usize::from(next_char)];
                let scaled1 = table[usize::from(c1)];
                if (scaled0 | scaled1) >= 0 {
                    let i0 = preshift + scaled0;
                    let i1 = compressed_reciprocal_divide(i0, recip) + scaled1;
                    preshift = div225(i0) + div15_exact(scaled1);
                    next_char = byte_at(s + word_length + 1);
                    let bit0 = bit_of(i0);
                    // SAFETY: each offset adds at most one to `listed`, so `listed` stays below
                    // the reserve (as in `for_each_cell`).
                    unsafe { list.add(listed).write((s as i32, i0 as u32)) };
                    listed += bit0 as usize;
                    let bit1 = bit_of(i1);
                    // SAFETY: as above.
                    unsafe { list.add(listed).write((s as i32 + 1, i1 as u32)) };
                    listed += bit1 as usize;
                    s += 2;
                    continue;
                }
            }
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
            // SAFETY: as above.
            unsafe { list.add(listed).write((s as i32, index as u32)) };
            listed += bit as usize;
            s += 1;
        }
        // SAFETY: the first `listed` places were written above, within the reserve.
        unsafe { cells.set_len(listed) };
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341 (quoted at
    // `for_each_cell`)
    /// EXPERIMENT (LOSAT_X_BPFLEAN): [`Self::x_list_pairs`] with its two-offset steps in a loop of
    /// their own, out of line ([`x_pairs_loop`]), so that their state (the offset, the list's
    /// length, the preshift and the next letter) stays in registers: in one loop with the other
    /// steps, the list's length was kept on the stack, loaded and stored at every pair. That loop
    /// takes NCBI's reciprocal for 15 as the constant `x_pairs_fit` asks for, computes in `u32`
    /// (every index is below 15^6) and reads the subject through a pointer a word length in. The
    /// first word, a letter outside the alphabet and the last offset take the steps of
    /// `x_list_pairs`, so the list, and the reserve past it, are the same.
    ///
    /// # Safety
    ///
    /// As [`Self::x_list_pairs`]: [`Self::x_pairs_fit`]`(subject, subject_length)` holds.
    #[inline(never)]
    pub(crate) unsafe fn x_list_pairs_lean(
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
            // Offset s_last, or a pair with a letter outside the alphabet: x_list_pairs's step.
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

    // No NCBI counterpart: chooses LOSAT's list form; both list the same offsets.
    /// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): whether blastp-fast lists this table's offsets without a
    /// branch (`for_each_cell::<true, _>`): at least [`x_scan_bf_min`] of every 10,000 backbone
    /// cells have hits. With fewer, few offsets pass `PV_TEST`, the branch is mostly predicted and
    /// the list's store at every offset costs more than the mispredicts it saves.
    pub(crate) fn x_scan_bf_pays(&self) -> bool {
        x_dense_enough(self.used_cells, self.starts.len() - 1, x_scan_bf_min())
    }
}

// No NCBI counterpart: the state of LOSAT's list around its two-offset steps.
/// EXPERIMENT (LOSAT_X_BPFLEAN): where [`x_pairs_loop`] starts and stops: the next offset `s`,
/// the list's length `listed`, the preshift before `s` and the letter at `s + word length - 1`.
struct XPairsAt {
    s: usize,
    listed: usize,
    preshift: i32,
    next_char: u8,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-341 (quoted at
// `XFlatCompressedLookup::for_each_cell`)
/// EXPERIMENT (LOSAT_X_BPFLEAN): the two-offset steps of [`XFlatCompressedLookup::x_list_pairs`]
/// from `at`, until fewer than two offsets up to `s_last` are left or a pair has a letter outside
/// the alphabet; `at` then holds that pair's first offset, for the caller's step. `sub_w` points a
/// word length into the subject, `table` holds the scaled letters and `pv` the presence bits. The
/// divisions are those of `x_list_pairs` with NCBI's reciprocal for 15, 286,331,154, as a constant.
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
        // i0 / 225 + C1 / 15, as in x_list_pairs.
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

// No NCBI counterpart: the rule of LOSAT_X_BPCLUTSCANBFMIN; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBF): `used` of `cells` backbone cells reach `min` per 10,000 (a
/// `min` above 10,000 takes no table).
fn x_dense_enough(used: usize, cells: usize, min: u64) -> bool {
    u128::from(used as u64) * 10_000 >= u128::from(min) * cells as u128
}

// No NCBI counterpart: a tuning value of LOSAT's list form; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBFMIN): the built-in value of [`x_scan_bf_min`]. Natively 500: on
/// the workstation the branch-free list cost x1.07 the cycles at 234 used cells per 10,000, x1.00 at
/// 442, x0.96 at 728 and x0.80 at 2,260. On wasm32 0 (every table): in Node the loop with the branch
/// is the slower one at every density measured, the branch-free list x0.81 the time at 119 and x0.80
/// at 234.
const X_SCAN_BF_MIN: u64 = if cfg!(target_arch = "wasm32") { 0 } else { 500 };

// No NCBI counterpart: reads LOSAT_X_BPCLUTSCANBFMIN once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPCLUTSCANBFMIN): the fewest backbone cells with hits, per 10,000, for which
/// blastp-fast's scan lists the offsets without a branch (built-in [`X_SCAN_BF_MIN`]); `0` every
/// table; a value that is not a number takes the built-in one.
fn x_scan_bf_min() -> u64 {
    static MIN: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BPCLUTSCANBFMIN")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(X_SCAN_BF_MIN)
    })
}

// No NCBI counterpart: how many cells ahead `for_each_cell` reads the slice bounds and the hits. Tuning
// only. LOSAT_X_BPFSTAGE's and LOSAT_X_BPFLEAN's loops (`blastp/x_seed_fused_clut.rs`) read ahead as
// far, over X_STARTS_AHEAD entries past their list, so the hits must not be read further ahead than the
// bounds.
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-149
// ```c
//     TBlastThreads the_threads(GetNumberOfThreads());
// ```
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
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
/// EXPERIMENT (LOSAT_X_BPCLUT): the flat table, by `pool`'s threads when `x_pool_builds` (the same
/// table, `x_clut_par.rs`; `x_flat_args` as LOSAT_X_BPLKPFLAT gives it), else by the serial builder
/// (which does not use `x_flat_args`: it runs none of the extra tasks and takes no vectors).
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
    let _ = (pool, x_flat_args);
    x_try_build_flat_compressed_lookup(word_size, scaled_threshold, query, locations, tables)
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
            && got.used_cells == want.used_cells
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
// ```c
//           numhits = backbone_cell->num_used;
//           ...
//              dest[0].qs_offsets.q_off = backbone_cell->query_offset;
// ```
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
mod tests {
    use super::*;
    use crate::utils::matrix::ncbistdaa;

    // A small deterministic generator (xorshift64).
    pub(super) struct XorShift(pub(super) u64);

    impl XorShift {
        pub(super) fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x % n as u64) as usize
        }
    }

    // NCBISTDAA letters inside the 15-letter compressed alphabet and outside it (NULLB, X, U, *, O).
    pub(super) const LETTERS: [u8; 23] = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 27,
    ];
    const OUTSIDE: [u8; 5] = [0, 21, 24, 25, 26];
    pub(super) const FEW_LETTERS: [u8; 4] =
        [ncbistdaa::A, ncbistdaa::L, ncbistdaa::S, ncbistdaa::K];

    fn flat_table(build: Option<XFlatBuild>) -> XFlatCompressedLookup {
        match build.expect("word size 5") {
            XFlatBuild::Table(table) => table,
            XFlatBuild::Refused => panic!("refused"),
            XFlatBuild::TooLarge => panic!("too large"),
        }
    }

    // The reference's rounds of `scan_subject` over a subject's one range (`blastp_for_each_offset_pair`
    // in two-hit mode), with an array of `array_size` pairs: the hits in order and the rounds.
    fn rounds<F>(
        word_length: i32,
        subject_length: usize,
        array_size: usize,
        mut scan: F,
    ) -> (Vec<(u32, u32)>, usize)
    where
        F: FnMut(&[(i32, i32)], &mut [OffsetPair], i32, &mut [i32; 3]) -> i32,
    {
        let seq_ranges = [(0, subject_length as i32)];
        let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - word_length];
        if scan_range[2] < scan_range[1] {
            scan_range[2] = scan_range[1];
        }
        let mut offset_pairs = vec![OffsetPair::default(); array_size];
        let mut hits = Vec::new();
        let mut count = 0;
        while scan_range[1] <= scan_range[2] {
            let n = scan(
                &seq_ranges,
                &mut offset_pairs,
                array_size as i32,
                &mut scan_range,
            );
            hits.extend(
                offset_pairs[..n as usize]
                    .iter()
                    .map(|pair| (pair.s_off, pair.q_off)),
            );
            count += 1;
        }
        (hits, count)
    }

    fn random_subject(rng: &mut XorShift, length: usize, outside_one_in: usize) -> Vec<u8> {
        let mut subject: Vec<u8> = (0..length)
            .map(|_| {
                if outside_one_in != 0 && rng.below(outside_one_in) == 0 {
                    OUTSIDE[rng.below(OUTSIDE.len())]
                } else {
                    LETTERS[rng.below(LETTERS.len())]
                }
            })
            .collect();
        subject.push(0);
        subject
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:825-841
    // ```c
    //     default:
    //     {
    //         ...
    //         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
    //
    //         if (cell_index == 0 ) { /* can't be empty => it's full  */
    // ```
    // The permutation for cells of 1 to 40 offsets: the reference's scan order of a cell
    // (`query_offsets_in_scan_order`, and the hits of the reference's `scan_subject` over a subject
    // with every word once) puts the offset of rank k at `x_emission_position(k, n)`, which is a
    // permutation of 0..n; the reference's overflow cells number cells(n).
    #[test]
    fn emission_positions_are_the_reference_scan_order_for_one_to_forty_offsets() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let mut lookup = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
        let mut subject = Vec::new();
        let mut expected_hits = Vec::new();
        for n in 1..=40u32 {
            let (word, index) = loop {
                let word: Vec<u8> = (0..5).map(|_| LETTERS[rng.below(LETTERS.len())]).collect();
                let encoded: Vec<u8> = word
                    .iter()
                    .map(|&letter| lookup.compress_table[letter as usize])
                    .collect();
                let index = compute_encoded_compressed_index(5, &encoded, 15);
                if lookup.backbone[index].num_used == 0 {
                    break (word, index);
                }
            };
            let cells_before = lookup.overflow_cells.len();
            let mut cells = 0usize;
            for k in 0..n {
                cells += usize::from(x_takes_overflow_cell(k));
                lookup.add_word_hit(index, (1000 * n + k) as i32);
            }
            assert_eq!(lookup.overflow_cells.len() - cells_before, cells, "n {n}");
            assert_eq!(cells, if n <= 5 { 0 } else { n as usize / 4 }, "n {n}");
            let mut expected = vec![-1i32; n as usize];
            for k in 0..n {
                let position = x_emission_position(k, n) as usize;
                assert_eq!(expected[position], -1, "n {n}: position {position} twice");
                expected[position] = (1000 * n + k) as i32;
            }
            assert_eq!(lookup.query_offsets_in_scan_order(index), expected, "n {n}");
            let s = subject.len() as u32;
            expected_hits.extend(expected.iter().map(|&q| (s, q as u32)));
            subject.extend_from_slice(&word);
            subject.push(OUTSIDE[rng.below(OUTSIDE.len())]);
        }
        lookup.finalize();
        let length = subject.len();
        subject.push(0);
        for array_size in [40, 41, 77, 4096 + 40] {
            let (hits, _) = rounds(5, length, array_size, |ranges, pairs, size, range| {
                lookup.scan_subject(&subject, ranges, pairs, size, range)
            });
            assert_eq!(hits, expected_hits, "array {array_size}");
        }
    }

    // Queries of `lengths` letters from `letters`, concatenated as BLASTP's lookup query (each query
    // followed by a sentinel). Each query has one lookup segment, or `masked` gaps (SEG-masked runs)
    // that cut it into several.
    pub(super) fn queries(
        rng: &mut XorShift,
        lengths: &[usize],
        letters: &[u8],
        masked: usize,
    ) -> (Vec<u8>, Vec<(i32, i32)>, Vec<usize>) {
        let mut query = Vec::new();
        let mut locations = Vec::new();
        // The first segment of each query.
        let mut first = Vec::new();
        for &length in lengths {
            let left = query.len() as i32;
            for _ in 0..length {
                query.push(letters[rng.below(letters.len())]);
            }
            first.push(locations.len());
            let mut from = left;
            let end = left + length as i32 - 1;
            for _ in 0..masked {
                if end - from < 30 {
                    break;
                }
                let gap = from + 5 + rng.below((end - from - 20) as usize) as i32;
                let gap_end = (gap + 1 + rng.below(12) as i32).min(end);
                if gap > from {
                    locations.push((from, gap - 1));
                }
                from = gap_end + 1;
            }
            if from <= end {
                locations.push((from, end));
            }
            query.push(0);
        }
        first.push(locations.len());
        (query, locations, first)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
    // ```c
    //     for (loc = location; loc; loc = loc->next){
    //         Int4 from = loc->ssr->left;
    //         Int4 to = loc->ssr->right - lookup->word_length + 1;
    //
    //         for (offset = from; offset <= to; offset++){
    //             s_CompressedAddWordHits(&info, query->sequence, offset);
    //         }
    //     }
    // ```
    // Random queries (few letters and low thresholds, so that cells overflow; all letters with letters
    // outside the compressed alphabet), SEG-masked runs (several segments per query), several NCBI
    // tables and a split batch (segments outside every table), threshold 0 (the words themselves
    // only): the flat table is the reference's, cell by cell in scan order, with the same PV array and
    // header.
    #[test]
    fn flat_table_equals_reference_on_random_queries() {
        let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
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
            // NCBI's tables: consecutive queries' segments; some batches are left out (split).
            let mut tables = Vec::new();
            let mut q = 0;
            while q < lengths.len() {
                let end = (q + 1 + rng.below(2)).min(lengths.len());
                if rng.below(3) != 0 {
                    tables.push(first[q]..first[end]);
                }
                q = end;
            }
            let reference =
                try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                    .expect("word size 5")
                    .expect("not refused");
            let flat = flat_table(x_try_build_flat_compressed_lookup(
                5, threshold, &query, &locations, &tables,
            ));
            x_assert_same_compressed_table(&flat, &reference);
            overflow_seen |= reference.longest_chain > 5;
        }
        assert!(overflow_seen, "cells with overflow lists");
        // No words at all (every segment shorter than a word).
        let (query, locations, _) = queries(&mut rng, &[3, 4, 2], &FEW_LETTERS, 0);
        let reference = try_build_blosum62_compressed_lookup(5, 1100, &query, &locations, &[])
            .expect("word size 5")
            .expect("not refused");
        let flat = flat_table(x_try_build_flat_compressed_lookup(
            5,
            1100,
            &query,
            &locations,
            &[],
        ));
        x_assert_same_compressed_table(&flat, &reference);
        assert!(flat.hits.is_empty() && flat.longest_chain == 0);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
    // ```c
    //     if (lookup->curr_overflow_cell ==
    //                         COMPRESSED_OVERFLOW_CELLS_IN_BANK) {
    //         /* need a new bank */
    //         Int4 bank_idx = lookup->curr_overflow_bank + 1;
    //         lookup->overflow_banks[bank_idx] = (CompressedOverflowCell*) malloc(
    // ```
    // The refusal at small capacities: on random queries and tables (with split batches), at
    // capacities around each table's own cells and its cells of the table of all batches, the flat
    // builder refuses exactly when the reference does, and otherwise builds the reference's table;
    // `x_table_fills_banks` is `compressed_table_fills_banks`.
    #[test]
    fn flat_refusal_is_the_references_at_small_capacities() {
        let mut rng = XorShift(0x1234_5678_9abc_def1);
        let mut refused_seen = 0usize;
        let mut built_seen = 0usize;
        let header = XHeader::new(5).expect("word size 5");
        let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
        let info = CompressedNeighborInfo::new(&matrix, 15);
        for round in 0..8 {
            let threshold = [1100, 1600, 2000][round % 3];
            let words = XWords::new(&header, threshold, &matrix, &info);
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
            // Each table's own cells (the reference's table of its segments).
            let mut capacities = vec![1usize];
            for table in &tables {
                let mut own = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
                own.threshold = threshold;
                own.overflow_cell_capacity = usize::MAX;
                own.add_neighboring_words(&matrix, &query, &locations[table.clone()]);
                let cells = own.overflow_cells.len();
                capacities.extend([cells.saturating_sub(1), cells, cells + 1, 2 * cells]);
                for capacity in [cells.saturating_sub(1), cells]
                    .into_iter()
                    .filter(|&c| c > 0)
                {
                    assert_eq!(
                        x_table_fills_banks(
                            &words,
                            header.backbone_size,
                            &query,
                            &locations[table.clone()],
                            capacity
                        ),
                        compressed_table_fills_banks(
                            5,
                            threshold,
                            &matrix,
                            &query,
                            &locations[table.clone()],
                            capacity
                        ),
                        "round {round}, capacity {capacity}"
                    );
                }
            }
            capacities.retain(|&capacity| capacity > 0);
            for capacity in capacities {
                let reference = try_build_blosum62_compressed_lookup_with_capacity(
                    5, threshold, &query, &locations, &tables, capacity,
                )
                .expect("word size 5");
                // The bucket fill (LOSAT_X_BPFFILL) builds the same table or declines.
                match x_try_build_flat_by_buckets(
                    5, threshold, &query, &locations, &tables, capacity,
                ) {
                    XBucketBuild::Declined => {}
                    XBucketBuild::Built(Some(XFlatBuild::Table(flat))) => {
                        let reference = reference.as_ref().unwrap_or_else(|_| {
                            panic!("round {round}, capacity {capacity}: buckets built, reference refused")
                        });
                        x_assert_same_compressed_table(&flat, reference);
                    }
                    XBucketBuild::Built(_) => {
                        panic!("round {round}, capacity {capacity}: buckets refused")
                    }
                }
                // Chunks of one append (a test after each query offset with appends, as the
                // reference tests), of a few, and the default (a test at each table's end).
                for chunk in [1, 7, X_COUNT_CHUNK] {
                    let built = x_try_build_flat_with_capacity(
                        5, threshold, &query, &locations, &tables, capacity, chunk,
                    )
                    .expect("word size 5");
                    match (&reference, built) {
                        (Err(()), XFlatBuild::Refused) => refused_seen += 1,
                        (Ok(reference), XFlatBuild::Table(flat)) => {
                            x_assert_same_compressed_table(&flat, reference);
                            built_seen += 1;
                        }
                        (reference, _) => panic!(
                            "round {round}, capacity {capacity}, chunk {chunk}: reference refused {}",
                            reference.is_err()
                        ),
                    }
                }
            }
        }
        assert!(refused_seen > 0 && built_seen > 0);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 (module header)
    // ```c
    //     Int4 num_entries = backbone_cell->num_used;
    //     ...
    //     backbone_cell->num_used++;
    // ```
    // The bucket fill (LOSAT_X_BPFFILL): on random queries (few letters and low thresholds, so that
    // cells overflow; letters outside the compressed alphabet), SEG-masked runs, NCBI's tables with
    // split batches, threshold 0, and capacities from 1 to NCBI's, it builds the count-and-fill
    // builder's table field by field and slot by slot, or declines where a table's appends pass the
    // capacity (never at NCBI's); a query of 2^22 letters or more does not fit its packed appends.
    #[test]
    fn x_bucket_fill_is_the_count_and_fill_build() {
        let mut rng = XorShift(0x5be0_cd19_137e_2179);
        let (mut built, mut declined, mut overflow_seen) = (0usize, 0usize, false);
        for round in 0..24 {
            let (letters, threshold): (&[u8], i32) = match round % 5 {
                0 => (&FEW_LETTERS[..], 1100),
                1 => (&FEW_LETTERS[..], 1600),
                2 => (&LETTERS[..], 1100),
                3 => (&[1u8, 11, 17, 21, 24, 10, 3][..], 2000),
                _ => (&FEW_LETTERS[..], 0),
            };
            let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(400)).collect();
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
            for capacity in [1usize, 50, 400, NCBI_OVERFLOW_CELL_CAPACITY] {
                let want = x_try_build_flat_with_capacity(
                    5,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    capacity,
                    X_COUNT_CHUNK,
                );
                let got = x_try_build_flat_by_buckets(
                    5, threshold, &query, &locations, &tables, capacity,
                );
                match (got, want) {
                    (XBucketBuild::Declined, _) => {
                        assert!(
                            capacity < NCBI_OVERFLOW_CELL_CAPACITY,
                            "round {round}: declined"
                        );
                        declined += 1;
                    }
                    (
                        XBucketBuild::Built(Some(XFlatBuild::Table(got))),
                        Some(XFlatBuild::Table(want)),
                    ) => {
                        x_assert_same_flat_tables(&got, &want);
                        overflow_seen |= want.longest_chain > 5;
                        built += 1;
                    }
                    (XBucketBuild::Built(_), _) => {
                        panic!("round {round}, capacity {capacity}: another outcome")
                    }
                }
            }
        }
        assert!(built > 0 && declined > 0 && overflow_seen);
        // Offsets up to 2^22 - 1 (bit 31 of the packed append set): a query of 2^22 - 1 letters
        // with one segment at its end; a query of 2^22 letters is declined, and so is word size 6
        // (the PV shift depends on the counts).
        let length = (1usize << 22) - 1;
        let mut query: Vec<u8> = (0..length)
            .map(|_| FEW_LETTERS[rng.below(FEW_LETTERS.len())])
            .collect();
        let locations = [(length as i32 - 400, length as i32 - 1)];
        let capacity = NCBI_OVERFLOW_CELL_CAPACITY;
        let want = x_try_build_flat_with_capacity(
            5,
            1100,
            &query,
            &locations,
            &[0..1],
            capacity,
            X_COUNT_CHUNK,
        );
        let got = x_try_build_flat_by_buckets(5, 1100, &query, &locations, &[0..1], capacity);
        let (XBucketBuild::Built(Some(XFlatBuild::Table(got))), Some(XFlatBuild::Table(want))) =
            (got, want)
        else {
            panic!("2^22 - 1 letters: not two tables");
        };
        assert!(got.hits.iter().any(|&offset| offset >= 1 << 21));
        x_assert_same_flat_tables(&got, &want);
        query.push(FEW_LETTERS[0]);
        let got = x_try_build_flat_by_buckets(5, 1100, &query, &locations, &[0..1], capacity);
        assert!(matches!(got, XBucketBuild::Declined), "2^22 letters");
        let got =
            x_try_build_flat_by_buckets(6, 1100, &query[..1000], &[(0, 999)], &[0..1], capacity);
        assert!(matches!(got, XBucketBuild::Declined), "word size 6");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:340-404
    // ```c
    //           if (numhits) {
    //           if (numhits <= (array_size - totalhits)) {
    //              ...
    //           else
    //               /* not enough space in the destination array */
    //           {
    //               s_range[1] = (Int4)(s - subject->sequence);
    //               return totalhits;
    //           }
    // ```
    // The scans: on random subjects (letters outside the alphabet none, rare, frequent; short
    // subjects), the flat `scan_subject` gives the reference's hits and rounds for arrays from the
    // longest cell up, and `for_each_cell` gives the same hits, with the branch-free list too
    // (LOSAT_X_BPCLUTSCANBF); `used_cells` counts the cells with hits.
    #[test]
    fn flat_scans_give_the_reference_scans_hits() {
        let mut rng = XorShift(0x3c6e_f372_fe94_f82b);
        for (threshold, letters) in [
            (1100, &FEW_LETTERS[..]),
            (900, &[1u8, 11, 17][..]),
            (1300, &LETTERS[..]),
        ] {
            let lengths = [120usize, 75, 200, 33];
            let (query, locations, _) = queries(&mut rng, &lengths, letters, 1);
            let tables = [0..locations.len()];
            let reference =
                try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                    .expect("word size 5")
                    .expect("not refused");
            let flat = flat_table(x_try_build_flat_compressed_lookup(
                5, threshold, &query, &locations, &tables,
            ));
            assert!(reference.longest_chain > 5, "cells with overflow lists");
            let used = flat.starts.windows(2).filter(|w| w[1] > w[0]).count();
            assert!(used > 0 && flat.used_cells == used, "used cells");
            let longest = reference.longest_chain as usize;
            let mut subjects = Vec::new();
            for outside_one_in in [0usize, 50, 6, 2] {
                for _ in 0..5 {
                    let length = 1 + rng.below(300);
                    let mut subject = random_subject(&mut rng, length, outside_one_in);
                    if length > 60 && rng.below(2) == 0 {
                        let from = rng.below(query.len() - 60);
                        let to = rng.below(length - 50);
                        subject[to..to + 50].copy_from_slice(&query[from..from + 50]);
                    }
                    subjects.push(subject);
                }
            }
            for length in 1..=6 {
                subjects.push(random_subject(&mut rng, length, 0));
            }
            let mut cells = Vec::new();
            for (n, subject) in subjects.iter().enumerate() {
                let length = subject.len() - 1;
                let mut expected = None;
                for array_size in [longest, longest + 1, 2 * longest + 3, 4096 + longest] {
                    let (want, want_rounds) =
                        rounds(5, length, array_size, |ranges, pairs, size, range| {
                            reference.scan_subject(subject, ranges, pairs, size, range)
                        });
                    let (got, got_rounds) =
                        rounds(5, length, array_size, |ranges, pairs, size, range| {
                            flat.scan_subject(subject, ranges, pairs, size, range)
                        });
                    assert_eq!(got, want, "subject {n}, array {array_size}");
                    assert_eq!(got_rounds, want_rounds, "subject {n}, array {array_size}");
                    expected.get_or_insert(want);
                }
                let mut fused = Vec::new();
                flat.for_each_cell::<false, _>(subject, length, &mut cells, |s, hits| {
                    fused.extend(hits.iter().map(|&q| (s as u32, q as u32)));
                });
                assert_eq!(Some(fused), expected, "subject {n}: for_each_cell");
                // LOSAT_X_BPCLUTSCANBF: the branch-free list, on the used scratch and on a new one.
                for mut scratch in [std::mem::take(&mut cells), Vec::new()] {
                    let mut fused = Vec::new();
                    flat.for_each_cell::<true, _>(subject, length, &mut scratch, |s, hits| {
                        fused.extend(hits.iter().map(|&q| (s as u32, q as u32)));
                    });
                    assert_eq!(
                        Some(fused),
                        expected,
                        "subject {n}: branch-free for_each_cell"
                    );
                    cells = scratch;
                }
                // LOSAT_X_BPFSTAGE: the list two offsets per step, on a subject of a word or more
                // with its sentinel; a shorter one, or one without the sentinel, is declined.
                assert_eq!(
                    flat.x_pairs_fit(subject, length),
                    length >= 5,
                    "subject {n}: fits"
                );
                assert!(
                    !flat.x_pairs_fit(&subject[..length], length),
                    "subject {n}: no sentinel"
                );
                if length < 5 {
                    continue;
                }
                let want_list = cells.clone();
                // SAFETY: x_pairs_fit holds (asserted above).
                unsafe { flat.x_list_pairs(subject, length, &mut cells) };
                assert_eq!(cells, want_list, "subject {n}: x_list_pairs's list");
                assert!(
                    cells.capacity() >= cells.len() + X_STARTS_AHEAD,
                    "subject {n}: reserve"
                );
                // LOSAT_X_BPFLEAN: the list with the two-offset steps out of line.
                let mut lean = Vec::new();
                // SAFETY: x_pairs_fit holds (asserted above).
                unsafe { flat.x_list_pairs_lean(subject, length, &mut lean) };
                assert_eq!(lean, want_list, "subject {n}: x_list_pairs_lean's list");
                assert!(
                    lean.capacity() >= lean.len() + X_STARTS_AHEAD,
                    "subject {n}: lean reserve"
                );
                let listed: Vec<(u32, u32)> = cells
                    .iter()
                    .flat_map(|&(s, index)| {
                        flat.cell(index as usize)
                            .iter()
                            .map(move |&q| (s as u32, q as u32))
                    })
                    .collect();
                assert_eq!(Some(listed), expected, "subject {n}: x_list_pairs");
            }
        }
    }

    // LOSAT_X_BPFSTAGE: `x * 19,088,744 >> 32` is `x / 225` below 15^6, NCBI's reciprocal
    // multiply is `x / 15` there (the indices of word lengths up to 6), and the multiply by the
    // inverse of 15 is `x / 15` on its multiples.
    #[test]
    fn the_pairs_scan_divisions_are_exact_below_fifteen_to_the_sixth() {
        for x in 0..15i64.pow(6) {
            assert_eq!((x * 19_088_744) >> 32, x / 225, "x {x}");
            assert_eq!(
                i64::from(compressed_reciprocal_divide(x as i32, 286_331_154)),
                x / 15,
                "x {x}"
            );
            if x % 15 == 0 {
                assert_eq!(
                    i64::from((x as u32).wrapping_mul(0xEEEE_EEEF)),
                    x / 15,
                    "x {x}"
                );
            }
        }
    }

    // LOSAT_X_BPFSTAGE: the list two offsets per step on word length 6 too (and LOSAT_X_BPFLEAN's),
    // and the tables and subjects it declines.
    #[test]
    fn the_pairs_scan_lists_word_length_six_and_declines_outside_its_range() {
        let mut rng = XorShift(0x6a09_e667_f3bc_c909);
        let (query, locations, _) = queries(&mut rng, &[200usize, 150], &LETTERS[..], 1);
        let tables = [0..locations.len()];
        let mut flat6 = flat_table(x_try_build_flat_compressed_lookup(
            6, 1100, &query, &locations, &tables,
        ));
        // A table with under 1 % of its cells used takes a coarser shift of its bits (NCBI's
        // `pv_array_bts`), which the scan declines; this one has one bit per cell.
        let probe = random_subject(&mut rng, 60, 0);
        assert_eq!(
            flat6.pv_array_bts, PV_ARRAY_BTS as i32,
            "a denser word-6 table"
        );
        assert!(flat6.x_pairs_fit(&probe, 60));
        flat6.pv_array_bts += 2;
        assert!(!flat6.x_pairs_fit(&probe, 60), "coarse bits");
        flat6.pv_array_bts -= 2;
        let mut cells = Vec::new();
        let mut found = 0;
        for outside_one_in in [0usize, 6, 2] {
            for length in [1usize, 5, 6, 7, 8, 60, 301] {
                let subject = random_subject(&mut rng, length, outside_one_in);
                assert_eq!(flat6.x_pairs_fit(&subject, length), length >= 6);
                if length < 6 {
                    continue;
                }
                let mut want = Vec::new();
                flat6.for_each_cell::<true, _>(&subject, length, &mut cells, |s, hits| {
                    want.extend(hits.iter().map(|&q| (s, q)));
                });
                let want_list = cells.clone();
                // SAFETY: x_pairs_fit holds (asserted above).
                unsafe { flat6.x_list_pairs(&subject, length, &mut cells) };
                assert_eq!(cells, want_list, "length {length}: list");
                let mut lean = Vec::new();
                // SAFETY: x_pairs_fit holds (asserted above).
                unsafe { flat6.x_list_pairs_lean(&subject, length, &mut lean) };
                assert_eq!(lean, want_list, "length {length}: lean list");
                found += want.len();
                let got: Vec<(i32, i32)> = cells
                    .iter()
                    .flat_map(|&(s, index)| flat6.cell(index as usize).iter().map(move |&q| (s, q)))
                    .collect();
                assert_eq!(got, want, "length {length}");
            }
        }
        assert!(found > 0, "word length 6 finds hits");
        let subject = random_subject(&mut rng, 50, 0);
        assert!(
            !flat6.x_pairs_fit(&subject, 51),
            "no byte past the last word"
        );
        assert!(!flat6.x_pairs_fit(&subject, 52), "length past the buffer");
        flat6.reciprocal_alphabet_size += 1;
        assert!(!flat6.x_pairs_fit(&subject, 50), "another reciprocal");
        flat6.reciprocal_alphabet_size -= 1;
        flat6.pv.pop();
        assert!(!flat6.x_pairs_fit(&subject, 50), "bits short of 15^6");
    }

    // LOSAT_X_BPCLUTSCANBFMIN: a table takes the branch-free list from `min` used cells per 10,000 up;
    // `0` takes it for every table, an empty one included.
    #[test]
    fn the_branch_free_list_starts_at_min_used_cells_per_ten_thousand() {
        assert!(x_dense_enough(500, 10_000, 500));
        assert!(!x_dense_enough(499, 10_000, 500));
        assert!(x_dense_enough(37_969, 759_376, 500));
        assert!(!x_dense_enough(37_968, 759_376, 500));
        assert!(x_dense_enough(0, 759_376, 0));
        assert!(!x_dense_enough(759_376, 759_376, u64::MAX));
    }
}
