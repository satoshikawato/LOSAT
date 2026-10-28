//! The BLASTX protein lookup table filled in place, with the neighbour-word sets memoised per
//! thread (LOSAT_X_BPLKPAHEAD: or kept by the blastp block search for its builds ahead,
//! `XAaSpares`).
//!
//! NCBI builds the table in two steps: `s_AddNeighboringWords` appends, for every distinct
//! query word in ascending backbone index, the word's query offsets to the chain of every
//! word it neighbours (`BlastLookupAddWordHit`), and `BlastAaLookupFinalize` copies the chains
//! into the thick backbone and the overflow array.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
//! ```c
//!     for (i = 0; i < lookup->backbone_size; i++) {
//!         if (exact_backbone[i] != NULL) {
//!             s_AddWordHits(lookup, matrix, query->sequence,
//!                           exact_backbone[i], query_bias, row_max);
//!             sfree(exact_backbone[i]);
//!         }
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590
//! ```c
//!         for (i = 0; i < alphabet_size; i++) {
//!             if (score + row[i] >= threshold) {
//!                 subject_word[current_pos] = i;
//!                 for (j = 0; j < offset_list[1]; j++) {
//!                     BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
//!                                           charsize, subject_word,
//!                                           query_bias + offset_list[j + 2]);
//!                 }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:74-76
//! ```c
//!     /* add the hit */
//!     chain[chain[1] + 2] = query_offset;
//!     chain[1]++;
//! ```
//!
//! This module computes the same table without the chains. For a source word `s` (a cell of
//! the exact backbone) let T(s) be the list of target cells it appends to, in NCBI's order:
//! `s` itself when `threshold == 0 || self_score < threshold` (blast_aalookup.c:504-509), then
//! the neighbour words in the order of the `s_AddWordHitsCore` recursion (546-606). The
//! reference appends, for s ascending and w in T(s), the offsets of s (in chain order) to w.
//! Here a count pass adds |offsets(s)| to count[w] for every such (s, w); the cells are laid
//! out from the counts by `BlastAaLookupFinalize`'s rule; and a fill pass writes, for s
//! ascending and w in T(s), the offsets of s at w's fill cursor.
//!
//! Identity argument. The fill pass performs the same sequence of (cell, query offset)
//! appends as the reference, so the entries of every cell are the reference chain's entries
//! in the same order. `num_used`, the pv bits, the overflow placement (cursor assigned in
//! ascending cell order to cells with more than AA_HITS_PER_CELL entries), the overflow size
//! and `longest_chain` are functions of the per-cell counts and follow the same rule as
//! `BlastAaLookupFinalize`; unused cell slots stay zero as after the C `calloc`. T(s) reads
//! only the residues of the word (the exact backbone index identifies them), the matrix, the
//! threshold, the word length and the alphabet, so a list memoised under the same key is the
//! list the recursion would produce again. The exact backbone is built by the reference
//! function (a partitioned build, LOSAT_X_BXBUILDPAR, takes the same cells and offsets from a
//! counting sort of that function's loop, `x_counted_sources`), and every other field of the table
//! is computed by the reference expressions. In diagnostics mode the reference builder runs (its
//! prints read the chains).

use super::*;
use std::cell::RefCell;

// EXPERIMENT (PAIRLUT): the per-group tables of the TBLASTX unit
// path built by this file's count, layout and fill passes on the lookup locations of one group; a
// child module so that it uses `x_target_list` and the memo below unchanged (see x_pair_lut.rs).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
// ```c
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
//             sfree(exact_backbone[i]);
//         }
//     }
// ```
// The module runs this loop's appends for the exact backbone of a subset of the lookup locations;
// it is used only when one of its switches is set.
#[path = "x_pair_lut.rs"]
pub(crate) mod x_pair_lut;

/// Neighbour lists of one (matrix, threshold, word length, alphabet) setting, by source cell.
struct XNeighborMemo {
    key: (ScoringMatrix, i32, usize, usize),
    lists: Vec<Option<Box<[u16]>>>,
}

thread_local! {
    static X_MEMO: RefCell<Option<XNeighborMemo>> = const { RefCell::new(None) };
}

// No NCBI counterpart (buffer reuse and one memo for all threads; the table is the one
// `x_build_direct` builds): LOSAT_X_BXLUTPOOL, in the strict set (on unless `=0` or
// `LOSAT_X_DEFAULT=0`). NCBI builds a BLASTX table with calloc'd vectors and frees them with the
// search (blast_aalookup.c BlastAaLookupTableNew / BlastAaLookupFinalize, cited above). LOSAT builds
// one per query batch or chunk searched side by side; a thread's glibc arena keeps what the thread
// allocated after it is freed, so every building thread kept a table's worth (3-6 MiB with the
// build's counts, exact backbone and its own neighbour memo). Pooled builds read one set of
// neighbour lists shared by all threads (`XSharedMemo`, filled as the per-thread memo is), take the
// table's vectors from the tables dropped before (`x_give_back_pooled`) and the build's counts and
// exact backbone from the builds before (`XBuildTemps`), each cleared and zero-filled to the sizes
// of the new table as calloc gives them; the pool holds as many sets as builds and tables were alive
// at once. The values written are those of `x_build_direct`.
/// Whether BLASTX builds its tables pooled (`LOSAT_X_BXLUTPOOL`).
pub(crate) fn x_bx_lut_pool() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXLUTPOOL"))
}

// No NCBI counterpart: buffer reuse only (above).
/// Whether a BLASTX table is built pooled, so that its vectors go back to the pool when it is
/// dropped: the switch is on and the in-place builder runs (diagnostics mode runs the reference
/// builder, whose tables are freed as before).
pub(crate) fn x_bx_pooled_builds() -> bool {
    x_bx_lut_pool() && !diagnostics_enabled()
}

// No NCBI counterpart (buffer reuse only): LOSAT_X_BXPOOLN1, in the strict set (on unless `=0`).
// The pool moves a table's vectors and its build's temporaries from the thread that drops them to
// the thread that builds next, so that no thread's glibc arena keeps a table's worth after the
// table is freed (`x_bx_lut_pool`). A search on one thread has no other thread: with the switch its
// tables are built by `x_build_direct` (its own neighbour memo, `vec!` vectors) and freed with the
// search's lookup (`blastx/seed.rs` `Lookup`), as NCBI builds and frees them and as `=0` of
// LOSAT_X_BXLUTPOOL does. Kept, the pool held a table and its build's temporaries through each
// query batch's end (the composition-based statistics and the traceback), which then took its
// memory on top: a 660 kb genome (AP027131) against 543 proteins (AP027133) at one thread had a
// peak `RssAnon` of 26.0 MiB, 22.6 without the pool. Searches of two threads or more build pooled as before.
/// Whether a BLASTX search of `threads` threads builds its tables pooled: the builds are pooled
/// ([`x_bx_pooled_builds`]) and, with `LOSAT_X_BXPOOLN1`, the search has two threads or more.
pub(crate) fn x_bx_pooled_for(threads: usize) -> bool {
    static N1: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if !x_bx_pooled_builds() {
        return false;
    }
    if !*N1.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXPOOLN1")) {
        return true;
    }
    let own = threads <= 1;
    crate::utils::xreach::x_reach!("LOSAT_X_BXPOOLN1", own);
    !own
}

/// No NCBI counterpart (memoised neighbour lists shared by threads, `XNeighborMemo` otherwise).
struct XSharedMemo {
    key: (ScoringMatrix, i32, usize, usize),
    lists: Box<[std::sync::OnceLock<Box<[u16]>>]>,
}

// No NCBI counterpart: the shared neighbour lists of the last key built pooled.
static X_SHARED_MEMO: std::sync::Mutex<Option<std::sync::Arc<XSharedMemo>>> =
    std::sync::Mutex::new(None);

// No NCBI counterpart: the shared lists for `key` (a new set where the kept one has another key).
fn x_shared_memo(
    key: (ScoringMatrix, i32, usize, usize),
    backbone_size: usize,
) -> std::sync::Arc<XSharedMemo> {
    let mut kept = x_pool_lock(&X_SHARED_MEMO);
    if let Some(memo) = kept
        .as_ref()
        .filter(|m| m.key == key && m.lists.len() == backbone_size)
    {
        return std::sync::Arc::clone(memo);
    }
    let memo = std::sync::Arc::new(XSharedMemo {
        key,
        lists: (0..backbone_size)
            .map(|_| std::sync::OnceLock::new())
            .collect(),
    });
    *kept = Some(std::sync::Arc::clone(&memo));
    memo
}

// No NCBI counterpart (scheduling only): the lists of `x_target_list`, made before the build that
// reads them; each is a function of its word alone (module comment), so a list made here is the
// one the build would make.
/// EXPERIMENT (LOSAT_X_BXLISTAHEAD): makes the shared neighbour lists (the pooled builds' memo,
/// [`x_shared_memo`]) of BLASTX's standard table (BLOSUM62, words of 3) with `threshold` for the
/// words of `query` at `locations` (as `x_for_each_exact_word` reads them), on the pool this runs
/// on, so that a pooled build of a table of these words finds them made.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(crate) fn x_bx_lists_ahead(query: &[u8], locations: &[(i32, i32)], threshold: i32) {
    use rayon::prelude::*;
    let (matrix, word_length) = (ScoringMatrix::Blosum62, 3);
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    if word_length * charsize > 16 {
        return;
    }
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);
    // Same expression as build_lookup_from_prepared (blast_aalookup.c:438-444).
    let row_max: Vec<i32> = (0..alphabet_size)
        .map(|i| {
            (0..alphabet_size)
                .map(|j| lookup_matrix_score(matrix, i, j))
                .max()
                .unwrap_or(i32::MIN)
        })
        .collect();
    let shared = x_shared_memo((matrix, threshold, word_length, alphabet_size), backbone_size);
    // The first offset of each word not listed yet.
    let mut first = vec![u32::MAX; backbone_size];
    x_for_each_exact_word(query, locations, word_length, charsize, |idx, offset| {
        if first[idx] == u32::MAX && shared.lists[idx].get().is_none() {
            first[idx] = offset as u32;
        }
    });
    let words: Vec<(usize, usize)> = first
        .iter()
        .enumerate()
        .filter(|(_, &offset)| offset != u32::MAX)
        .map(|(idx, &offset)| (idx, offset as usize))
        .collect();
    words.par_iter().with_min_len(64).for_each(|&(idx, offset)| {
        shared.lists[idx].get_or_init(|| {
            x_target_list(
                &query[offset..],
                alphabet_size,
                word_length,
                charsize,
                &row_max,
                threshold,
                matrix,
            )
        });
    });
}

/// No NCBI counterpart: the neighbour lists a build reads, its thread's or the shared ones.
enum XLists<'a> {
    Own(&'a mut Vec<Option<Box<[u16]>>>),
    Shared(&'a [std::sync::OnceLock<Box<[u16]>>]),
}

impl XLists<'_> {
    // No NCBI counterpart: the memo's lookup (the lists are those s_AddWordHitsCore would produce).
    fn get_or_make(&mut self, idx: usize, make: impl FnOnce() -> Box<[u16]>) -> &[u16] {
        match self {
            Self::Own(lists) => lists[idx].get_or_insert_with(make),
            Self::Shared(lists) => lists[idx].get_or_init(make),
        }
    }

    // No NCBI counterpart: as `get_or_make`, for a cell the count pass listed.
    fn get(&self, idx: usize) -> &[u16] {
        match self {
            Self::Own(lists) => lists[idx].as_deref(),
            Self::Shared(lists) => lists[idx].get().map(|list| &list[..]),
        }
        .expect("listed in the count pass")
    }
}

/// No NCBI counterpart: a build's temporary vectors, kept between pooled builds (`cells`, `starts`
/// and `offsets`: the present source cells of a partitioned build, `x_counted_sources`, which keeps
/// its cursors in `base`; the build's `counts` holds its ranges' counts, `x_partitioned_passes`).
#[derive(Default)]
struct XBuildTemps {
    counts: Vec<u32>,
    base: Vec<u32>,
    fill: Vec<u32>,
    chains: Option<LookupBackboneChains>,
    cells: Vec<u32>,
    starts: Vec<u32>,
    offsets: Vec<i32>,
}

// SAFETY: the exact backbone's block pointers point into slabs it owns (boxes, which do not move
// when the backbone moves) and are only read through it; a kept backbone is reset (every cell
// empty, every slab spare) before the next build, so no pointer outlives its build. The kept
// temporaries are moved between threads under the pool's lock, never shared.
unsafe impl Send for XBuildTemps {}

// No NCBI counterpart: the kept build temporaries and table vectors (buffer reuse only).
static X_BX_TEMPS: std::sync::Mutex<Vec<XBuildTemps>> = std::sync::Mutex::new(Vec::new());
static X_BX_TABLES: std::sync::Mutex<Vec<XAaVectors>> = std::sync::Mutex::new(Vec::new());

// No NCBI counterpart: a lock held only around a move (a poisoned lock's value is still valid).
fn x_pool_lock<T>(mutex: &std::sync::Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

// No NCBI counterpart: `vector` cleared and zero-filled to `n` entries (as calloc gives them); one
// short of room is replaced by one of `n` entries, or with `slack` (the overflow, whose length
// changes from table to table) `n + n / 4`, not the doubling of `resize`.
fn x_zeroed<T: Default + Clone>(mut vector: Vec<T>, n: usize, slack: bool) -> Vec<T> {
    if vector.capacity() < n {
        vector = Vec::with_capacity(if slack { n + n / 4 } else { n });
    }
    vector.clear();
    vector.resize(n, T::default());
    vector
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467 (see
// `x_build_direct`)
// ```c
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
// ```
/// [`x_build_direct`] with the shared neighbour lists and kept vectors (`LOSAT_X_BXLUTPOOL`); with
/// `parts` above 1, its passes on that many tasks of the pool it runs on (`LOSAT_X_BXBUILDPAR`,
/// [`x_partitioned_passes`]).
pub(super) fn x_build_direct_pooled(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    parts: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
    let mut unused = None;
    x_build_direct_in::<BackboneCell>(
        prepared,
        threshold,
        matrix,
        word_length,
        &mut unused,
        None,
        true,
        parts,
        false,
    )
    .map(|table| table.expect("pooled vectors grow to any table"))
}

// No NCBI counterpart: buffer reuse (`LOSAT_X_BXLUTPOOL`); NCBI frees the vectors.
/// Keeps `table`'s backbone, overflow and pv vectors for the next pooled build (the table is left
/// empty).
pub(crate) fn x_give_back_pooled(table: &mut BlastAaLookupTable) {
    let vectors = XAaVectors {
        backbone: std::mem::take(&mut table.backbone),
        overflow: std::mem::take(&mut table.overflow),
        pv: std::mem::take(&mut table.pv),
        grow: false,
        ..XAaVectors::default()
    };
    if !vectors.is_empty() {
        x_pool_lock(&X_BX_TABLES).push(vectors);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
// ```c
//       lookup->thick_backbone =
//             calloc(lookup->backbone_size, sizeof(AaLookupBackboneCell));
// ...
//       pv = lookup->pv = (PV_ARRAY_TYPE *) calloc(
//             (lookup->backbone_size >> PV_ARRAY_BTS) + 1,
//              sizeof(PV_ARRAY_TYPE));
// ...
//       if (overflow_cells_needed > 0) {
//           lookup->overflow =
//                  calloc(overflow_cells_needed, sizeof(Int4));
// ```
// NCBI allocates a zeroed table for each search and frees it with the search. Buffer reuse only:
// the vectors taken are cleared and filled with zeros to the table's sizes before the build writes
// them, as the C `calloc` gives them.
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the backbone, overflow and pv vectors of a BLASTP table, kept
/// for another table's build ([`XAaSpares`]).
#[derive(Default)]
pub(crate) struct XAaVectors {
    backbone: Vec<BackboneCell>,
    overflow: Vec<i32>,
    pv: Vec<u32>,
    /// EXPERIMENT (LOSAT_X_BPLK16): the backbone and overflow vectors of a table of `u16` hits;
    /// the pair of the other width is empty (`x_take_zeroed` and `XAaSpares::x_make_room` drop
    /// the pair a build or a room does not use, so that a spare keeps one width).
    backbone16: Vec<BackboneCell16>,
    overflow16: Vec<u16>,
    /// EXPERIMENT (LOSAT_X_BPLKREUSE): a vector short of the table's size grows in place
    /// (`x_take_zeroed`) instead of refusing the build.
    grow: bool,
}

// No NCBI counterpart: buffer reuse and the width of the table's hits (LOSAT_X_BPLK16); it
// changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLK16): a cell width the in-place builder ([`x_build_direct_in`]) builds
/// tables of, and the vectors of [`XAaVectors`] it builds on.
trait XVec: XCell {
    /// The width's backbone and overflow vectors of `vectors`.
    fn vectors(vectors: &mut XAaVectors) -> (&mut Vec<Self>, &mut Vec<Self::Hit>);

    /// Drops the other width's backbone and overflow vectors.
    fn drop_other(vectors: &mut XAaVectors);

    /// The table of a partitioned build (LOSAT_X_BXBUILDPAR: `Int4` hits only, BLASTX's pooled
    /// build) as this width's table. The narrow width builds no partitioned table.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn partitioned(table: BlastAaLookupTable) -> BlastAaLookupTable<Self>;
}

impl XVec for BackboneCell {
    #[inline(always)]
    fn vectors(vectors: &mut XAaVectors) -> (&mut Vec<Self>, &mut Vec<i32>) {
        (&mut vectors.backbone, &mut vectors.overflow)
    }

    #[inline(always)]
    fn drop_other(vectors: &mut XAaVectors) {
        vectors.backbone16 = Vec::new();
        vectors.overflow16 = Vec::new();
    }

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[inline(always)]
    fn partitioned(table: BlastAaLookupTable) -> BlastAaLookupTable {
        table
    }
}

impl XVec for BackboneCell16 {
    fn vectors(vectors: &mut XAaVectors) -> (&mut Vec<Self>, &mut Vec<u16>) {
        (&mut vectors.backbone16, &mut vectors.overflow16)
    }

    fn drop_other(vectors: &mut XAaVectors) {
        vectors.backbone = Vec::new();
        vectors.overflow = Vec::new();
    }

    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn partitioned(_table: BlastAaLookupTable) -> BlastAaLookupTable<Self> {
        unreachable!("a partitioned build (pooled, parts above 1) makes tables of Int4 hits only")
    }
}

impl XAaVectors {
    // No NCBI counterpart: buffer reuse only (above).
    /// The three vectors, cleared and zero-filled to `backbone_size`, `pv_size` and `overflow_size`
    /// entries, where each has the capacity (nothing is allocated); else `None`, nothing taken.
    /// With `grow` (LOSAT_X_BPLKREUSE) a vector short of its size grows first, by a reallocation
    /// of its own allocation to that size (`reserve_exact`; the first allocation where it has
    /// none), never a new allocation and a free of the old one: glibc
    /// moves a mapped allocation with `mremap`, so the freed table never raises its dynamic mmap
    /// threshold, which a free of a mapped allocation does (the next table then lands in the main
    /// heap, whose free space stays resident below the allocations made after it).
    ///
    /// EXPERIMENT (LOSAT_X_BPLK16): of the width `C`; the sizes are those of the table without the
    /// zero padding of its width (`XCell::PAD_CELLS`, `XCell::PAD_HITS`: none for `Int4` hits),
    /// which the vectors are given as well. The vectors of the other width are dropped.
    fn x_take_zeroed<C: XVec>(
        &mut self,
        backbone_size: usize,
        pv_size: usize,
        overflow_size: usize,
    ) -> Option<(Vec<C>, Vec<u32>, Vec<C::Hit>)> {
        let backbone_size = backbone_size + C::PAD_CELLS;
        let overflow_size = overflow_size + C::PAD_HITS;
        let grow = self.grow;
        let pv_capacity = self.pv.capacity();
        let (own_backbone, own_overflow) = C::vectors(self);
        if !grow
            && (own_backbone.capacity() < backbone_size
                || pv_capacity < pv_size
                || own_overflow.capacity() < overflow_size)
        {
            return None;
        }
        let (mut backbone, mut overflow) =
            (std::mem::take(own_backbone), std::mem::take(own_overflow));
        // The other width's vectors are freed before this width's grow.
        C::drop_other(self);
        // Grown and filled in the order of before LOSAT_X_BPLK16: backbone, pv, overflow.
        backbone.clear();
        backbone.reserve_exact(backbone_size);
        backbone.resize(backbone_size, C::default());
        let mut pv = std::mem::take(&mut self.pv);
        pv.clear();
        pv.reserve_exact(pv_size);
        pv.resize(pv_size, 0);
        overflow.clear();
        overflow.reserve_exact(overflow_size);
        overflow.resize(overflow_size, C::Hit::default());
        Some((backbone, pv, overflow))
    }

    // No NCBI counterpart: buffer reuse only.
    fn is_empty(&self) -> bool {
        self.backbone.capacity() == 0
            && self.overflow.capacity() == 0
            && self.backbone16.capacity() == 0
            && self.overflow16.capacity() == 0
            && self.pv.capacity() == 0
    }
}

// No NCBI counterpart: buffer reuse, and the thread that allocates a table's vectors and the
// neighbour lists (LOSAT_X_BPLKPAHEAD, `blastp/blast_engine/x_lookup_ahead.rs`); it changes no
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): what the blastp block search (`blastp/x_query_blocks.rs`
/// `XBlockSubjects`) keeps for the build of the next block's exact table ahead of that block, on a
/// task of the search's pool (`x_build_ahead`): the vectors of the last table its search no longer
/// reads (`give`), the vectors the next build takes (`x_make_room`), and the neighbour lists that
/// build uses. Every vector of a table built ahead is allocated by the block search's thread
/// (`x_make_room`), never by the pool's task (a pool thread's allocations stay in its own glibc
/// arena, which keeps them resident after they are freed), and one set of neighbour lists serves
/// every build ahead (the block search's thread's, moved here), not one set per pool thread.
/// EXPERIMENT (LOSAT_X_BPLKREUSE): where no table is built ahead (one thread), each block's own
/// build takes the vectors of the block before from here (`x_build_on_block_spare`).
#[derive(Default)]
pub(crate) struct XAaSpares {
    block: std::sync::Mutex<XAaVectors>,
    ahead: std::sync::Mutex<XAaVectors>,
    memo: std::sync::Mutex<Option<XNeighborMemo>>,
}

// No NCBI counterpart: a lock held only around a move; a poisoned lock's value is still a valid
// spare (buffer reuse only).
fn x_locked<T>(mutex: &std::sync::Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

impl XAaSpares {
    // No NCBI counterpart: where the build ahead applies; it reads switches and changes no value
    // NCBI computes.
    /// Whether exact blastp's table is built by [`x_build_direct`] (outside diagnostics mode,
    /// whose reference builder prints its phases), the build [`Self::x_build_ahead`] makes.
    pub(crate) fn x_applies() -> bool {
        !diagnostics_enabled()
    }

    // No NCBI counterpart: buffer reuse only (`XAaSpares`).
    /// Keeps `table`'s backbone, overflow and pv vectors (in place of any kept before); the rest of
    /// the table is dropped. The next build ahead (`x_make_room`) or, with LOSAT_X_BPLKREUSE, the
    /// next block's own build (`x_build_on_block_spare`) takes them.
    pub(crate) fn give(&self, table: BlastAaLookupTable) {
        *x_locked(&self.block) = XAaVectors {
            backbone: table.backbone,
            overflow: table.overflow,
            pv: table.pv,
            grow: false,
            ..XAaVectors::default()
        };
    }

    // No NCBI counterpart: buffer reuse only (`XAaSpares`).
    /// EXPERIMENT (LOSAT_X_BPLK16): [`Self::give`] for a table of `u16` hits (its padding stays in
    /// the vectors; the next build of that width clears and fills them to its own sizes).
    pub(crate) fn give16(&self, table: BlastAaLookupTable<BackboneCell16>) {
        *x_locked(&self.block) = XAaVectors {
            backbone16: table.backbone,
            overflow16: table.overflow,
            pv: table.pv,
            grow: false,
            ..XAaVectors::default()
        };
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (see
    // `build_ncbi_lookup_from_prepared`)
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    // ```
    /// EXPERIMENT (LOSAT_X_BPLKREUSE): `x_bp_build_lookup(true, queries, threshold, karlin_params,
    /// check_ideal)` (where [`Self::x_applies`]) on the calling thread, with its neighbour lists,
    /// for a block search whose tables are not built ahead: the table built on the vectors the
    /// block before gave back (`give`), each grown in place where it is short
    /// (`XAaVectors::x_take_zeroed`), or on new ones for the first block. Where the in-place
    /// builder declines (it does not where `x_applies`), the reference builder builds it and the
    /// vectors stay kept.
    pub(crate) fn x_build_on_block_spare(
        &self,
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        karlin_params: &KarlinParams,
        check_ideal: bool,
    ) -> (BlastAaLookupTable, Vec<QueryContext>) {
        let _ = karlin_params;
        let prepared = x_bp_prepared_query(queries, check_ideal);
        // build_lookup_from_prepared returns the prepared contexts unchanged.
        let contexts = prepared.contexts.clone();
        match self.x_on_block_spare::<BackboneCell>(prepared, threshold) {
            Ok(table) => (table, contexts),
            Err(prepared) => build_lookup_from_prepared(
                prepared,
                threshold,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
            ),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (see
    // `build_ncbi_lookup_from_prepared`)
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    // ```
    /// The in-place build of [`Self::x_build_on_block_spare`] with the cells of width `C`: the
    /// table built on the vectors the block before gave back (`give`, `give16`), each grown in
    /// place where it is short, or `Err(prepared)` where the in-place builder declines (the vectors
    /// stay kept).
    fn x_on_block_spare<C: XVec>(
        &self,
        prepared: PreparedLookupQuery,
        threshold: i32,
    ) -> std::result::Result<BlastAaLookupTable<C>, PreparedLookupQuery> {
        let mut room = std::mem::take(&mut *x_locked(&self.block));
        room.grow = true;
        let built = X_MEMO.with(|memo| {
            x_build_direct_in::<C>(
                prepared,
                threshold,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
                &mut memo.borrow_mut(),
                Some(&mut room),
                false,
                1,
                crate::seqdb::x_bp_db_keep(),
            )
        });
        match built {
            Ok(table) => Ok(table.expect("grown vectors have the room for any table")),
            Err(prepared) => {
                room.grow = false;
                *x_locked(&self.block) = room;
                Err(prepared)
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (see
    // `build_ncbi_lookup_from_prepared`)
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    // ```
    /// EXPERIMENT (LOSAT_X_BPLK16): [`Self::x_build_on_block_spare`] with the width of the table's
    /// hits chosen by the build (`x_bplk16_width`): `u16` where the lookup query is at most 65,535
    /// letters and the in-place builder applies, else `Int4` as before.
    pub(crate) fn x_build_on_block_spare16(
        &self,
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        check_ideal: bool,
    ) -> (XBpLookup, Vec<QueryContext>) {
        let prepared = x_bp_prepared_query(queries, check_ideal);
        // build_lookup_from_prepared returns the prepared contexts unchanged.
        let contexts = prepared.contexts.clone();
        let prepared = if x_bplk16_fits(&prepared) {
            match self.x_on_block_spare::<BackboneCell16>(prepared, threshold) {
                Ok(table) => {
                    x_bplk16_shadow(queries, threshold, check_ideal, &table);
                    return (XBpLookup::Narrow(table), contexts);
                }
                Err(prepared) => prepared,
            }
        } else {
            match self.x_on_block_spare::<BackboneCell>(prepared, threshold) {
                Ok(table) => return (XBpLookup::Wide(table), contexts),
                Err(prepared) => prepared,
            }
        };
        let (table, contexts) = build_lookup_from_prepared(
            prepared,
            threshold,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        );
        (XBpLookup::Wide(table), contexts)
    }

    // No NCBI counterpart: buffer reuse, and the thread that allocates the next table's vectors
    // (`XAaSpares`).
    /// On the block search's thread, before the build of the next block's table ahead of that
    /// block: keeps for that build vectors with the room for `backbone_size` cells, `pv_size` pv
    /// words and `overflow` overflow entries: its own (left by a build that made no table), else
    /// those the block before gave back (`give`; dropped here otherwise), and each one short of its
    /// room is replaced here by an empty one of that capacity (`overflow_grown`, at least
    /// `overflow`, for the overflow). Nothing is written into them. Also keeps this thread's
    /// neighbour lists for the builds ahead, if none are kept. Returns the number of vectors
    /// allocated.
    pub(crate) fn x_make_room(
        &self,
        backbone_size: usize,
        pv_size: usize,
        overflow: usize,
        overflow_grown: usize,
    ) -> usize {
        self.x_make_room_of(backbone_size, pv_size, overflow, overflow_grown, false)
    }

    // No NCBI counterpart: buffer reuse only (`XAaSpares`).
    /// EXPERIMENT (LOSAT_X_BPLK16): [`Self::x_make_room`] for the vectors of a table of `u16` hits
    /// (`backbone_size` and `overflow` count their padding, as the vectors of the table built
    /// before have it); the vectors of the other width are dropped. A block whose table turns out
    /// to have the other width builds on its own vectors.
    pub(crate) fn x_make_room16(
        &self,
        backbone_size: usize,
        pv_size: usize,
        overflow: usize,
        overflow_grown: usize,
    ) -> usize {
        self.x_make_room_of(backbone_size, pv_size, overflow, overflow_grown, true)
    }

    /// [`Self::x_make_room`] for the vectors of the width `narrow` names.
    fn x_make_room_of(
        &self,
        backbone_size: usize,
        pv_size: usize,
        overflow: usize,
        overflow_grown: usize,
        narrow: bool,
    ) -> usize {
        let mut kept = std::mem::take(&mut *x_locked(&self.ahead));
        let given_back = std::mem::take(&mut *x_locked(&self.block));
        if kept.is_empty() {
            kept = given_back;
        }
        let mut allocated = 0;
        if narrow {
            BackboneCell16::drop_other(&mut kept);
        } else {
            BackboneCell::drop_other(&mut kept);
        }
        // The allocations in the order of the vectors' use before LOSAT_X_BPLK16: backbone, pv,
        // overflow.
        if narrow {
            if kept.backbone16.capacity() < backbone_size {
                kept.backbone16 = Vec::with_capacity(backbone_size);
                allocated += 1;
            }
        } else if kept.backbone.capacity() < backbone_size {
            kept.backbone = Vec::with_capacity(backbone_size);
            allocated += 1;
        }
        if kept.pv.capacity() < pv_size {
            kept.pv = Vec::with_capacity(pv_size);
            allocated += 1;
        }
        if narrow {
            if kept.overflow16.capacity() < overflow {
                kept.overflow16 = Vec::with_capacity(overflow_grown.max(overflow));
                allocated += 1;
            }
        } else if kept.overflow.capacity() < overflow {
            kept.overflow = Vec::with_capacity(overflow_grown.max(overflow));
            allocated += 1;
        }
        *x_locked(&self.ahead) = kept;
        let mut memo = x_locked(&self.memo);
        if memo.is_none() {
            *memo = X_MEMO.with(|own| own.borrow_mut().take());
        }
        allocated
    }

    // No NCBI counterpart: buffer reuse and the thread that keeps the neighbour lists
    // (`XAaSpares`).
    /// On the block search's thread, where it may build a table itself (a block whose build ahead
    /// made no table, the last block, after which the search's thread may search again): the
    /// neighbour lists kept here go back to this thread, if it has none, and with `release` the
    /// vectors kept here are dropped (the last block: no build ahead follows).
    pub(crate) fn x_home(&self, release: bool) {
        let memo = x_locked(&self.memo).take();
        if memo.is_some() {
            X_MEMO.with(|own| {
                let mut own = own.borrow_mut();
                if own.is_none() {
                    *own = memo;
                }
            });
        }
        if release {
            *x_locked(&self.block) = XAaVectors::default();
            *x_locked(&self.ahead) = XAaVectors::default();
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (see
    // `build_ncbi_lookup_from_prepared`)
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    // ```
    /// `x_bp_build_lookup(queries, threshold, check_ideal)` (where
    /// [`Self::x_applies`]), the table built on the vectors `x_make_room` kept and with the
    /// neighbour lists kept here (both taken while it builds and given back after, the vectors only
    /// where no table was made). `None` where no table was made: the vectors lack the room, or the
    /// in-place builder declines (the reference builder would build it).
    pub(crate) fn x_build_ahead(
        &self,
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        karlin_params: &KarlinParams,
        check_ideal: bool,
    ) -> Option<(BlastAaLookupTable, Vec<QueryContext>)> {
        let _ = karlin_params;
        let prepared = x_bp_prepared_query(queries, check_ideal);
        // build_lookup_from_prepared returns the prepared contexts unchanged.
        let contexts = prepared.contexts.clone();
        let mut room = std::mem::take(&mut *x_locked(&self.ahead));
        let mut memo = x_locked(&self.memo).take();
        let built = x_build_direct_in::<BackboneCell>(
            prepared,
            threshold,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
            &mut memo,
            Some(&mut room),
            false,
            1,
            crate::seqdb::x_bp_db_keep(),
        );
        *x_locked(&self.memo) = memo;
        *x_locked(&self.ahead) = room;
        match built {
            Ok(Some(table)) => Some((table, contexts)),
            Ok(None) | Err(_) => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (see
    // `build_ncbi_lookup_from_prepared`)
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    // ```
    /// EXPERIMENT (LOSAT_X_BPLK16): [`Self::x_build_ahead`] with the width of the table's hits
    /// chosen by the build (`x_bplk16_fits`): `u16` where the lookup query is at most 65,535
    /// letters, else `Int4` as before; `None` where no table was made (the kept vectors of the
    /// width lack the room, or the in-place builder declines).
    pub(crate) fn x_build_ahead16(
        &self,
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        check_ideal: bool,
    ) -> Option<(XBpLookup, Vec<QueryContext>)> {
        let prepared = x_bp_prepared_query(queries, check_ideal);
        // build_lookup_from_prepared returns the prepared contexts unchanged.
        let contexts = prepared.contexts.clone();
        let narrow = x_bplk16_fits(&prepared);
        let mut room = std::mem::take(&mut *x_locked(&self.ahead));
        let mut memo = x_locked(&self.memo).take();
        let built = if narrow {
            x_build_direct_in::<BackboneCell16>(
                prepared,
                threshold,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
                &mut memo,
                Some(&mut room),
                false,
                1,
                crate::seqdb::x_bp_db_keep(),
            )
            .map(|table| table.map(XBpLookup::Narrow))
        } else {
            x_build_direct_in::<BackboneCell>(
                prepared,
                threshold,
                ScoringMatrix::Blosum62,
                LOOKUP_WORD_LENGTH,
                &mut memo,
                Some(&mut room),
                false,
                1,
                crate::seqdb::x_bp_db_keep(),
            )
            .map(|table| table.map(XBpLookup::Wide))
        };
        *x_locked(&self.memo) = memo;
        *x_locked(&self.ahead) = room;
        match built {
            Ok(Some(table)) => {
                if let XBpLookup::Narrow(narrow) = &table {
                    x_bplk16_shadow(queries, threshold, check_ideal, narrow);
                }
                Some((table, contexts))
            }
            Ok(None) | Err(_) => None,
        }
    }
}

/// EXPERIMENT (LOSAT_X_BPLK16): the table of a blastp block: `Int4` hits (every table before the
/// switch, and the blocks whose lookup query is above 65,535 letters) or `u16` hits.
pub(crate) enum XBpLookup {
    Wide(BlastAaLookupTable),
    Narrow(BlastAaLookupTable<BackboneCell16>),
}

// No NCBI counterpart: chooses the width of a table's hits (LOSAT_X_BPLK16); it changes no value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLK16): whether the table built from `prepared` takes `u16` hits: every
/// query offset (below the lookup query's length, the concatenation of the contexts with their
/// sentinels) and every cell's hit count (at most one hit per offset) fits 16 bits.
fn x_bplk16_fits(prepared: &PreparedLookupQuery) -> bool {
    // The expression of x_build_direct_in (BLAST_SequenceBlk length).
    let query_length = prepared
        .contexts
        .last()
        .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
        .unwrap_or(0);
    (0..=i32::from(u16::MAX)).contains(&query_length)
}

// No NCBI counterpart: a diagnostic of LOSAT_X_BPLK16; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLK16SHADOW): with the variable set, the table of `u16` hits just built for
/// `queries` is compared with the table of `Int4` hits of the reference builder
/// (`build_lookup_from_prepared`) cell by cell ([`x_assert_narrow_same_as_wide`]); a difference
/// panics with the cell's index.
fn x_bplk16_shadow(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
    narrow: &BlastAaLookupTable<BackboneCell16>,
) {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if !*ON.get_or_init(|| std::env::var_os("LOSAT_X_BPLK16SHADOW").is_some()) {
        return;
    }
    let (wide, _) = build_lookup_from_prepared(
        x_bp_prepared_query(queries, check_ideal),
        threshold,
        ScoringMatrix::Blosum62,
        LOOKUP_WORD_LENGTH,
    );
    x_assert_narrow_same_as_wide("LOSAT_X_BPLK16SHADOW:", narrow, &wide);
}

// No NCBI counterpart: the check of LOSAT_X_BPLK16SHADOW and of the unit tests; it changes no value
// NCBI computes.
/// Asserts that a table of `u16` hits holds what the table of `Int4` hits of the same query holds:
/// for every cell the same `num_used` and the same hit list (in order), the same pv and
/// `longest_chain` and scalars, and the padding (`XCell::PAD_CELLS` zero cells, `XCell::PAD_HITS`
/// zero hits) after the cells and the overflow. The first difference panics with a message that
/// starts with `prefix` and names the cell.
pub(crate) fn x_assert_narrow_same_as_wide(
    prefix: &str,
    narrow: &BlastAaLookupTable<BackboneCell16>,
    wide: &BlastAaLookupTable,
) {
    assert_eq!(
        narrow.backbone.len(),
        wide.backbone.len() + BackboneCell16::PAD_CELLS,
        "{prefix} backbone size differs"
    );
    assert_eq!(
        narrow.overflow.len(),
        wide.overflow.len() + BackboneCell16::PAD_HITS,
        "{prefix} overflow size differs"
    );
    for (idx, (a, b)) in narrow.backbone.iter().zip(&wide.backbone).enumerate() {
        let (a_hits, b_hits) = (a.hits(&narrow.overflow), b.hits(&wide.overflow));
        assert!(
            a.count() == b.count()
                && a_hits.len() == b_hits.len()
                && a_hits.iter().zip(b_hits).all(|(a, b)| a.q() == b.q()),
            "{prefix} backbone cell {idx} differs ({} {:?} vs {} {:?})",
            a.count(),
            a_hits,
            b.count(),
            b_hits
        );
    }
    for (idx, cell) in narrow.backbone[wide.backbone.len()..].iter().enumerate() {
        assert!(
            cell.num_used == 0 && cell.entries == [0; AA_HITS_PER_CELL],
            "{prefix} padding cell {idx} is not zero"
        );
    }
    for (idx, (a, b)) in narrow.overflow.iter().zip(&wide.overflow).enumerate() {
        assert!(
            a.q() == b.q(),
            "{prefix} overflow entry {idx} differs ({a} vs {b})"
        );
    }
    assert!(
        narrow.overflow[wide.overflow.len()..]
            .iter()
            .all(|&h| h == 0),
        "{prefix} overflow padding is not zero"
    );
    assert!(narrow.pv == wide.pv, "{prefix} pv differs");
    assert!(
        narrow.frame_bases == wide.frame_bases
            && narrow.num_contexts == wide.num_contexts
            && narrow.query_length == wide.query_length
            && narrow.word_length == wide.word_length
            && narrow.alphabet_size == wide.alphabet_size
            && narrow.charsize == wide.charsize
            && narrow.mask == wide.mask
            && narrow.longest_chain == wide.longest_chain
            && narrow.threshold == wide.threshold
            && narrow.row_max == wide.row_max,
        "{prefix} table fields differ"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:562-563,580-582,600-604
// ```c
//     score -= info->row_max[query_word[current_pos]];
//     row = info->matrix[query_word[current_pos]];
// ...
//         for (i = 0; i < alphabet_size; i++) {
//             if (score + row[i] >= threshold) {
//                 subject_word[current_pos] = i;
// ...
//     for (i = 0; i < alphabet_size; i++) {
//         if (score + row[i] >= threshold) {
//             subject_word[current_pos] = i;
//             s_AddWordHitsCore(info, score + row[i], current_pos + 1);
// ```
// The recursion of `blast_aa_add_word_hits_core` (the port of s_AddWordHitsCore) with the same
// tests in the same order; where the C code calls BlastLookupAddWordHit for a complete word it
// records the word's backbone index (blast_lookup.c:46 ComputeTableIndex, the same shifts as
// `blast_lookup_add_word_hit`).
#[allow(clippy::too_many_arguments)]
fn x_record_neighbors(
    query_word: &[u8],
    subject_word: &mut [u8; 32],
    alphabet_size: usize,
    wordsize: usize,
    charsize: usize,
    row_max: &[i32],
    threshold: i32,
    matrix: ScoringMatrix,
    score: i32,
    current_pos: usize,
    out: &mut Vec<u16>,
) {
    let query_residue = query_word[current_pos] as usize;
    let score = score - row_max[query_residue];
    if current_pos == wordsize - 1 {
        for residue in 0..alphabet_size {
            let residue_score = lookup_matrix_score(matrix, query_residue, residue);
            if score + residue_score >= threshold {
                subject_word[current_pos] = residue as u8;
                out.push(x_word_index(subject_word, wordsize, charsize));
            }
        }
        return;
    }
    for residue in 0..alphabet_size {
        let residue_score = lookup_matrix_score(matrix, query_residue, residue);
        if score + residue_score >= threshold {
            subject_word[current_pos] = residue as u8;
            x_record_neighbors(
                query_word,
                subject_word,
                alphabet_size,
                wordsize,
                charsize,
                row_max,
                threshold,
                matrix,
                score + residue_score,
                current_pos + 1,
                out,
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:44-46
// ```c
//     /* compute the backbone cell to update */
//
//     index = ComputeTableIndex(wordsize, charsize, seq);
// ```
// Same shifts as `blast_lookup_add_word_hit`; the caller checks that the index fits 16 bits.
#[inline(always)]
fn x_word_index(word: &[u8], wordsize: usize, charsize: usize) -> u16 {
    let mut index = 0usize;
    for &residue in word.iter().take(wordsize) {
        index = (index << charsize) | residue as usize;
    }
    index as u16
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:490-496,504-509,539-543
// ```c
//     w = query + offset_list[2];
//
//     /* Compute the self-score of this word */
//
//     score = matrix[w[0]][w[0]];
//     for (i = 1; i < lookup->word_length; i++)
//         score += matrix[w[i]][w[i]];
// ...
//     if (lookup->threshold == 0 || score < lookup->threshold) {
//         for (i = 0; i < offset_list[1]; i++) {
//             BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
//                                   lookup->charsize, w,
//                                   query_bias + offset_list[i + 2]);
//         }
// ...
//     score = row_max[w[0]];
//     for (i = 1; i < lookup->word_length; i++)
//         score += row_max[w[i]];
//
//     s_AddWordHitsCore(&info, score, 0);
// ```
// T(s) for the word at `query_word`: the exact cell when s_AddWordHits adds the exact offsets,
// then the neighbour cells in recursion order (none when threshold == 0, lines 518-519).
#[allow(clippy::too_many_arguments)]
fn x_target_list(
    query_word: &[u8],
    alphabet_size: usize,
    wordsize: usize,
    charsize: usize,
    row_max: &[i32],
    threshold: i32,
    matrix: ScoringMatrix,
) -> Box<[u16]> {
    let mut out = Vec::new();
    let mut self_score =
        lookup_matrix_score(matrix, query_word[0] as usize, query_word[0] as usize);
    for residue in query_word.iter().take(wordsize).skip(1) {
        let residue = *residue as usize;
        self_score += lookup_matrix_score(matrix, residue, residue);
    }
    if threshold == 0 || self_score < threshold {
        out.push(x_word_index(query_word, wordsize, charsize));
    }
    if threshold != 0 {
        let mut max_score = row_max[query_word[0] as usize];
        for residue in query_word.iter().take(wordsize).skip(1) {
            max_score += row_max[*residue as usize];
        }
        let mut subject_word = [0u8; 32];
        x_record_neighbors(
            query_word,
            &mut subject_word,
            alphabet_size,
            wordsize,
            charsize,
            row_max,
            threshold,
            matrix,
            max_score,
            0,
            &mut out,
        );
    }
    out.into_boxed_slice()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469
// ```c
//     exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
// ...
//     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
//                                       lookup->charsize, lookup->word_length,
//                                       query, location);
// ...
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-297,338-358
// ```c
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (lookup->thin_backbone[i]) {
// ...
//             if (lookup->thin_backbone[i][1] > AA_HITS_PER_CELL){
// ...
//                 overflow_cells_needed += lookup->thin_backbone[i][1];
//             }
//             if (lookup->thin_backbone[i][1] > longest_chain)
//                 longest_chain = lookup->thin_backbone[i][1];
// ...
//       for (i = 0; i < lookup->backbone_size; i++) {
//         /* if there are hits there, */
//         if (lookup->thin_backbone[i] ) {
//             Int4 * dest = NULL;
//             /* set the corresponding bit in the pv_array */
//             PV_SET(pv, i, PV_ARRAY_BTS);
//             bbc[i].num_used = lookup->thin_backbone[i][1];
//             /* if there are three or fewer hits, */
//             if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
//                 /* copy them into the thick_backbone cell */
//                 dest = bbc[i].payload.entries;
//             else /* more than three hits; copy to overflow array */
//             {
//                 bbc[i].payload.overflow_cursor = overflow_cursor;
//                 dest = (Int4 *) lookup->overflow;
//                 dest += overflow_cursor;
//                 overflow_cursor += lookup->thin_backbone[i][1];
//             }
//             for (j=0; j <lookup->thin_backbone[i][1]; j++)
//                 dest[j] = lookup->thin_backbone[i][j + 2];
// ```
/// The table `build_lookup_from_prepared` returns, built without chains (see the module
/// comment). `Err(prepared)` when the reference builder must run (diagnostics mode, or a
/// backbone index wider than 16 bits).
pub(super) fn x_build_direct(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
    X_MEMO.with(|memo| {
        x_build_direct_in::<BackboneCell>(
            prepared,
            threshold,
            matrix,
            word_length,
            &mut memo.borrow_mut(),
            None,
            false,
            1,
            false,
        )
        .map(|table| table.expect("fresh vectors have the room for any table"))
    })
}

/// EXPERIMENT (LOSAT_X_BPDBKEEP): [`x_build_direct`] for BLASTP, whose chains glibc maps and
/// unmaps without raising its mmap threshold where the switch is on
/// (`LookupBackboneChains::x_new_unthresholded`).
fn x_bp_build_direct(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
    X_MEMO.with(|memo| {
        x_build_direct_in::<BackboneCell>(
            prepared,
            threshold,
            matrix,
            word_length,
            &mut memo.borrow_mut(),
            None,
            false,
            1,
            crate::seqdb::x_bp_db_keep(),
        )
        .map(|table| table.expect("fresh vectors have the room for any table"))
    })
}

/// [`x_build_direct`] with the neighbour lists of `memo` (the calling thread's, or those an
/// `XAaSpares` keeps), and, with `room` (LOSAT_X_BPLKPAHEAD), the table's backbone, overflow and pv
/// vectors taken from `room` instead of allocated: `Ok(None)` (nothing taken) where one of them
/// lacks the capacity. A pooled build with `parts` above 1 runs its passes on that many tasks of
/// the pool it runs on ([`x_partitioned_passes`]).
///
/// EXPERIMENT (LOSAT_X_BPLK16): the table has the cells of width `C`: `BackboneCell` (the table of
/// every program before the switch) or `BackboneCell16`, whose caller has checked that the lookup
/// query is at most 65,535 letters (`x_bplk16_fits`) and passes no `pooled` or `parts` above 1. A
/// table of `u16` hits has `XCell::PAD_CELLS` zero cells after its `backbone_size` cells and
/// `XCell::PAD_HITS` zero hits after its overflow; neither is in any count.
#[allow(clippy::too_many_arguments)]
fn x_build_direct_in<C: XVec>(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    memo: &mut Option<XNeighborMemo>,
    room: Option<&mut XAaVectors>,
    pooled: bool,
    parts: usize,
    x_unthresholded: bool,
) -> std::result::Result<Option<BlastAaLookupTable<C>>, PreparedLookupQuery> {
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    if diagnostics_enabled() || word_length * charsize > 16 || word_length > 32 {
        return Err(prepared);
    }
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);
    let PreparedLookupQuery {
        concat_query,
        lookup_locations,
        frame_bases,
        contexts,
        skipped_seg_mask: _,
    } = prepared;
    // Same expression as build_lookup_from_prepared (BLAST_SequenceBlk length).
    let query_length = contexts
        .last()
        .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
        .unwrap_or(0);
    // Same expression as build_lookup_from_prepared (blast_aalookup.c:438-444).
    let row_max: Vec<i32> = (0..alphabet_size)
        .map(|i| {
            (0..alphabet_size)
                .map(|j| lookup_matrix_score(matrix, i, j))
                .max()
                .unwrap_or(i32::MIN)
        })
        .collect();
    let mut temps = if pooled {
        x_pool_lock(&X_BX_TEMPS).pop().unwrap_or_default()
    } else {
        XBuildTemps::default()
    };
    let key = (matrix, threshold, word_length, alphabet_size);
    // Dispatch point (LOSAT_X_BXBUILDPAR, in the strict set): a pooled build (the shared lists)
    // with parts asked for runs its passes by source ranges on the pool, on the exact backbone's
    // cells and offsets as a counting sort gives them.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pooled && parts > 1 && room.is_none() {
        x_counted_sources(
            &concat_query,
            &lookup_locations,
            word_length,
            charsize,
            backbone_size,
            &mut temps,
        );
        let shared = x_shared_memo(key, backbone_size);
        let (backbone, pv, overflow, longest_chain) = x_partitioned_passes(
            parts,
            backbone_size,
            &shared.lists,
            &|first: usize| {
                x_target_list(
                    &concat_query[first..],
                    alphabet_size,
                    word_length,
                    charsize,
                    &row_max,
                    threshold,
                    matrix,
                )
            },
            &mut temps,
        );
        x_pool_lock(&X_BX_TEMPS).push(temps);
        return Ok(Some(C::partitioned(BlastAaLookupTable {
            backbone,
            overflow,
            pv,
            frame_bases,
            num_contexts: contexts.len(),
            query_length,
            word_length: word_length as i32,
            alphabet_size: alphabet_size as i32,
            charsize: charsize as i32,
            mask: mask as i32,
            longest_chain,
            threshold,
            row_max,
        })));
    }
    // The reference exact backbone (blast_aalookup.c:448-456).
    let mut exact_backbone = match temps.chains.take() {
        Some(mut chains) => {
            chains.x_reset(backbone_size);
            chains
        }
        // EXPERIMENT (LOSAT_X_BPDBKEEP): the chains' memory freed without raising glibc's mmap
        // threshold.
        None if x_unthresholded => LookupBackboneChains::x_new_unthresholded(backbone_size),
        None => LookupBackboneChains::new(backbone_size),
    };
    let mut skipped_invalid_residue = 0usize;
    blast_lookup_index_query_exact_matches(
        &mut exact_backbone,
        word_length as i32,
        charsize as i32,
        word_length as i32,
        &concat_query,
        &lookup_locations,
        &mut skipped_invalid_residue,
    );

    let shared = pooled.then(|| x_shared_memo(key, backbone_size));
    if shared.is_none()
        && memo
            .as_ref()
            .is_none_or(|m| m.key != key || m.lists.len() != backbone_size)
    {
        *memo = Some(XNeighborMemo {
            key,
            lists: std::iter::repeat_with(|| None)
                .take(backbone_size)
                .collect(),
        });
    }
    let mut lists = match &shared {
        Some(shared) => XLists::Shared(&shared.lists),
        None => XLists::Own(&mut memo.as_mut().expect("memo set above").lists),
    };

    // Count pass: count[w] = number of appends the reference makes to cell w.
    let mut counts: Vec<u32> = if pooled {
        x_zeroed(std::mem::take(&mut temps.counts), backbone_size, false)
    } else {
        vec![0; backbone_size]
    };
    for idx in 0..backbone_size {
        let Some(offset_list) = exact_backbone[idx].as_deref() else {
            continue;
        };
        let k = lookup_chain_num_used(offset_list) as u32;
        let list = lists.get_or_make(idx, || {
            // blast_aalookup.c:490: w = query + offset_list[2];
            let first = usize::try_from(offset_list[2]).expect("NCBI BLAST query offset must fit");
            x_target_list(
                &concat_query[first..],
                alphabet_size,
                word_length,
                charsize,
                &row_max,
                threshold,
                matrix,
            )
        });
        for &w in list.iter() {
            counts[w as usize] += k;
        }
    }

    // Layout: BlastAaLookupFinalize's rule applied to the counts.
    let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
    let mut overflow_size = 0usize;
    let mut longest_chain: i32 = 0;
    for &count in &counts {
        if count as usize > AA_HITS_PER_CELL {
            overflow_size += count as usize;
        }
        if (count as i32) > longest_chain {
            longest_chain = count as i32;
        }
    }
    let (mut backbone, mut pv, mut overflow): (Vec<C>, Vec<u32>, Vec<C::Hit>) = match room {
        None if pooled => {
            let mut kept = x_pool_lock(&X_BX_TABLES).pop().unwrap_or_default();
            let (kept_backbone, kept_overflow) = C::vectors(&mut kept);
            let (kept_backbone, kept_overflow) =
                (std::mem::take(kept_backbone), std::mem::take(kept_overflow));
            (
                x_zeroed(kept_backbone, backbone_size + C::PAD_CELLS, false),
                x_zeroed(kept.pv, pv_size, false),
                x_zeroed(kept_overflow, overflow_size + C::PAD_HITS, true),
            )
        }
        None => (
            vec![C::default(); backbone_size + C::PAD_CELLS],
            vec![0u32; pv_size],
            vec![C::Hit::default(); overflow_size + C::PAD_HITS],
        ),
        Some(room) => match room.x_take_zeroed::<C>(backbone_size, pv_size, overflow_size) {
            Some(vectors) => vectors,
            None => return Ok(None),
        },
    };
    let mut base: Vec<u32> = if pooled {
        x_zeroed(std::mem::take(&mut temps.base), backbone_size, false)
    } else {
        vec![0; backbone_size]
    };
    let mut overflow_cursor = 0usize;
    for (idx, &count) in counts.iter().enumerate() {
        if count == 0 {
            continue;
        }
        pv_set(&mut pv, idx);
        backbone[idx].set_count(count as usize);
        if count as usize > AA_HITS_PER_CELL {
            backbone[idx].set_cursor(overflow_cursor);
            base[idx] = overflow_cursor as u32;
            overflow_cursor += count as usize;
        }
    }

    // Fill pass: the reference's appends, in the reference's order.
    let mut fill: Vec<u32> = if pooled {
        x_zeroed(std::mem::take(&mut temps.fill), backbone_size, false)
    } else {
        vec![0; backbone_size]
    };
    for idx in 0..backbone_size {
        let Some(offset_list) = exact_backbone[idx].as_deref() else {
            continue;
        };
        // query_bias is 0 here (build_lookup_from_prepared passes 0).
        let offsets = lookup_chain_entries(offset_list);
        let list = lists.get(idx);
        for &w in list {
            let w = w as usize;
            let at = fill[w] as usize;
            if counts[w] as usize <= AA_HITS_PER_CELL {
                C::Hit::x_store(
                    &mut backbone[w].entries_mut()[at..at + offsets.len()],
                    offsets,
                );
            } else {
                let start = base[w] as usize + at;
                C::Hit::x_store(&mut overflow[start..start + offsets.len()], offsets);
            }
            fill[w] += offsets.len() as u32;
        }
    }

    if pooled {
        temps.counts = counts;
        temps.base = base;
        temps.fill = fill;
        temps.chains = Some(exact_backbone);
        x_pool_lock(&X_BX_TEMPS).push(temps);
    }
    Ok(Some(BlastAaLookupTable {
        backbone,
        overflow,
        pv,
        frame_bases,
        num_contexts: contexts.len(),
        query_length,
        word_length: word_length as i32,
        alphabet_size: alphabet_size as i32,
        charsize: charsize as i32,
        mask: mask as i32,
        longest_chain,
        threshold,
        row_max,
    }))
}

// No NCBI counterpart: spins of `x_claimed`'s caller before it yields while it waits for parts.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const X_CLAIM_SPINS: u32 = 1 << 12;

// No NCBI counterpart (scheduling only): `f` on each of `parts` on tasks of the pool the caller runs
// on. The caller and the helpers it spawns claim the parts in turn; once every part is claimed, the
// caller waits for the parts others still run without taking other work of the pool (a part is
// short, and a job taken while waiting, such as another batch's search, could hold the caller far
// longer than the part; the scope's own wait at its end, for a helper that claimed nothing, can
// still run one). A part that panics still counts as ended, and the scope raises its panic.
/// EXPERIMENT (LOSAT_X_BXBUILDPAR): `f(part)` for each part, its results in part order.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_claimed<T: Send, R: Send>(parts: Vec<T>, f: &(dyn Fn(T) -> R + Sync)) -> Vec<R> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct Ended<'a>(&'a AtomicUsize);
    impl Drop for Ended<'_> {
        fn drop(&mut self) {
            self.0.fetch_add(1, Ordering::Release);
        }
    }
    let n = parts.len();
    let slots: Vec<std::sync::Mutex<Option<T>>> = parts
        .into_iter()
        .map(|part| std::sync::Mutex::new(Some(part)))
        .collect();
    let results: Vec<std::sync::Mutex<Option<R>>> =
        (0..n).map(|_| std::sync::Mutex::new(None)).collect();
    let (next, ended) = (AtomicUsize::new(0), AtomicUsize::new(0));
    let work = || loop {
        let i = next.fetch_add(1, Ordering::Relaxed);
        if i >= n {
            return;
        }
        let _ended = Ended(&ended);
        let part = x_pool_lock(&slots[i])
            .take()
            .expect("each part claimed once");
        let result = f(part);
        *x_pool_lock(&results[i]) = Some(result);
    };
    rayon::in_place_scope(|scope| {
        for _ in 1..n.min(rayon::current_num_threads()) {
            scope.spawn(|_| work());
        }
        work();
        let mut spins = 0u32;
        while ended.load(Ordering::Acquire) < n {
            // A part is short; past a bound the helper running it was likely preempted (more
            // threads than CPUs), so give the CPU away instead of spinning through its timeslice.
            if spins < X_CLAIM_SPINS {
                spins += 1;
                std::hint::spin_loop();
            } else {
                std::thread::yield_now();
            }
        }
    });
    results
        .into_iter()
        .map(|result| {
            result
                .into_inner()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .expect("each part ended")
        })
        .collect()
}

// No NCBI counterpart: `vector` with `n` entries for a build that writes every entry itself: the
// entries kept are not cleared (new ones are zero), and one short of room is replaced as by
// `x_zeroed`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_sized<T: Default + Clone>(mut vector: Vec<T>, n: usize, slack: bool) -> Vec<T> {
    if vector.capacity() < n {
        vector = Vec::with_capacity(if slack { n + n / 4 } else { n });
    }
    vector.truncate(n);
    vector.resize(n, T::default());
    vector
}

/// No NCBI counterpart: the most source ranges of a partitioned build (each holds a count per
/// backbone cell, 112 KiB for words of three letters, kept with the build's temporaries: about
/// threads x 112 KiB for the K builds at a time, up to 8 parts each).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const X_MOST_PARTS: usize = 8;

/// No NCBI counterpart: the weight of a source whose neighbour list is not made yet, per offset
/// (the ranges are cut by offsets times list length).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const X_UNLISTED_LENGTH: usize = 24;

/// No NCBI counterpart: the bit of a range's position that says the position is in the overflow
/// (a cell of more than AA_HITS_PER_CELL entries), not in the cell's own entries.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const X_IN_OVERFLOW: u32 = 1 << 31;

/// No NCBI counterpart: where the source ranges of a partitioned build write their appends (the
/// table's backbone and overflow, written by the ranges side by side at disjoint positions).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XDisjointOut {
    backbone: *mut BackboneCell,
    overflow: *mut i32,
    #[cfg(debug_assertions)]
    lens: (usize, usize),
}

// SAFETY: the tasks of the fill write through the pointers only, each at positions no other task
// writes (`XDisjointOut::put`), while nothing else reads or writes the vectors.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
unsafe impl Sync for XDisjointOut {}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl XDisjointOut {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:74-76 (an append at the
    // chain's end, cited in the module comment); the copy of BlastAaLookupFinalize (blast_aalookup.c
    // 746-747, cited at `x_build_direct`) to the cell's entries or the overflow.
    /// Writes `src` at `at` of cell `w`: the cell's entries from `at`, or with `X_IN_OVERFLOW` set
    /// the overflow from `at`'s other bits. Copies of one to four offsets are made element by element
    /// (no call of the C library's memmove for the short lists of most cells and sources).
    ///
    /// # Safety
    /// `w` is a cell of the backbone, the positions `at..at + src.len()` lie in the cell's entries or
    /// in the overflow, and no other thread reads or writes them while this runs.
    #[inline(always)]
    unsafe fn put(&self, w: usize, at: u32, src: &[i32]) {
        let dst = if at & X_IN_OVERFLOW == 0 {
            #[cfg(debug_assertions)]
            debug_assert!(w < self.lens.0 && at as usize + src.len() <= AA_HITS_PER_CELL);
            // SAFETY: (caller) `w` is a cell of the backbone; no reference to the cell is made.
            unsafe { std::ptr::addr_of_mut!((*self.backbone.add(w)).entries) }
                .cast::<i32>()
                .wrapping_add(at as usize)
        } else {
            let at = (at & !X_IN_OVERFLOW) as usize;
            #[cfg(debug_assertions)]
            debug_assert!(at + src.len() <= self.lens.1);
            self.overflow.wrapping_add(at)
        };
        // SAFETY: (caller) `dst..dst + src.len()` lies in the cell's entries or in the overflow and
        // is written by this call alone.
        unsafe {
            match src {
                [s0] => dst.write(*s0),
                [s0, s1] => {
                    dst.write(*s0);
                    dst.add(1).write(*s1);
                }
                [s0, s1, s2] => {
                    dst.write(*s0);
                    dst.add(1).write(*s1);
                    dst.add(2).write(*s2);
                }
                [s0, s1, s2, s3] => {
                    dst.write(*s0);
                    dst.add(1).write(*s1);
                    dst.add(2).write(*s2);
                    dst.add(3).write(*s3);
                }
                src => std::ptr::copy_nonoverlapping(src.as_ptr(), dst, src.len()),
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:79-132 (the loop of
// `blast_lookup_index_query_exact_matches`, backbone.rs, which carries it)
// ```c
// for (loc = locations; loc; loc = loc->next) {
//     Int4 from = loc->ssr->left;
//     Int4 to = loc->ssr->right;
//     if (word_length > to - from + 1)
//         continue;
//     seq = query->sequence + from;
//     word_target = seq + lut_word_length;
//     for (offset = from; offset <= to; offset++, seq++) {
//         if (seq >= word_target) {
//             BlastLookupAddWordHit(backbone, lut_word_length, charsize,
//                                   seq - lut_word_length, offset - lut_word_length);
//         }
//         if (*seq & invalid_mask)
//             word_target = seq + lut_word_length + 1;
//     }
//     if (seq >= word_target) {
//         BlastLookupAddWordHit(backbone, lut_word_length, charsize,
//                               seq - lut_word_length, offset - lut_word_length);
//     }
// }
// ```
// The same loop with the same tests; where it calls BlastLookupAddWordHit, `visit` gets the
// word's backbone index (blast_lookup.c:46 ComputeTableIndex, the shifts of
// `blast_lookup_add_word_hit`) and the query offset the call appends.
/// EXPERIMENT (LOSAT_X_BXBUILDPAR): `visit(index, offset)` for each word the exact backbone's
/// build adds, in the order it adds them.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_for_each_exact_word(
    query: &[u8],
    locations: &[(i32, i32)],
    word_length: usize,
    charsize: usize,
    mut visit: impl FnMut(usize, usize),
) {
    let invalid_mask: u8 = 0xffu8 << (charsize as u32);
    for &(from, to) in locations {
        if word_length as i32 > to - from + 1 {
            continue;
        }
        let from = usize::try_from(from).expect("NCBI BLAST query offset must fit usize");
        let to = usize::try_from(to).expect("NCBI BLAST query offset must fit usize");
        let mut word_target = from + word_length;
        for seq in from..=to {
            if seq >= word_target {
                let start = seq - word_length;
                visit(
                    x_word_index(&query[start..], word_length, charsize) as usize,
                    start,
                );
            }
            if query[seq] & invalid_mask != 0 {
                word_target = seq + word_length + 1;
            }
        }
        if to + 1 >= word_target {
            let start = to + 1 - word_length;
            visit(
                x_word_index(&query[start..], word_length, charsize) as usize,
                start,
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-456 (cited at
// `x_build_direct`): the exact backbone, a chain per backbone cell holding the query offsets of the
// cell's words in the order BlastLookupIndexQueryExactMatches adds them. No NCBI counterpart for the
// form: a stable counting sort of the words of `x_for_each_exact_word` by cell (count, then each
// cell's start in ascending cell order, then each word's offset at its cell's cursor in the loop's
// order) gives every cell its chain's offsets in chain order, the cells that have a chain in
// ascending order, with no chain allocated.
/// EXPERIMENT (LOSAT_X_BXBUILDPAR): the exact backbone of a partitioned build as plain slices:
/// `temps.cells` the cells that have a chain, ascending, and cell i's offsets (chain order)
/// `temps.offsets[temps.starts[i]..temps.starts[i + 1]]`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_counted_sources(
    query: &[u8],
    locations: &[(i32, i32)],
    word_length: usize,
    charsize: usize,
    backbone_size: usize,
    temps: &mut XBuildTemps,
) {
    let mut at = x_zeroed(std::mem::take(&mut temps.base), backbone_size, false);
    let mut words = 0usize;
    x_for_each_exact_word(query, locations, word_length, charsize, |idx, _| {
        at[idx] += 1;
        words += 1;
    });
    let (cells, starts, offsets) = (&mut temps.cells, &mut temps.starts, &mut temps.offsets);
    cells.clear();
    starts.clear();
    let mut running = 0u32;
    for (idx, slot) in at.iter_mut().enumerate() {
        if *slot != 0 {
            cells.push(idx as u32);
            starts.push(running);
            running += std::mem::replace(slot, running);
        }
    }
    starts.push(running);
    offsets.clear();
    offsets.resize(words, 0);
    x_for_each_exact_word(query, locations, word_length, charsize, |idx, offset| {
        offsets[at[idx] as usize] =
            i32::try_from(offset).expect("NCBI BLAST query offset must fit Int4");
        at[idx] += 1;
    });
    temps.base = at;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467 and 283-358 (cited
// at `x_build_direct`): the appends of s_AddWordHits for s ascending, then BlastAaLookupFinalize's
// layout. No NCBI counterpart for the parts (scheduling only): NCBI builds a table on one thread.
// Here the passes of `x_build_direct_in` run on tasks of the pool. The present source cells, in
// ascending order, are cut into P contiguous ranges (of about equal appends).
// (a) Each range counts its own appends per target cell: c_p[w], the offsets its sources append to
//     w (making the neighbour lists not made yet; each is a function of its word alone, module
//     comment).
// (b) Each stripe of target cells (a multiple of 32 cells, so that no pv word is shared) sums them:
//     count[w] = sum over p of c_p[w], the reference's count, from which it sets num_used and the
//     pv bit of each of its cells, and its cells' overflow size and longest count.
// (c) Each stripe, after the overflow of the stripes before it (their sizes summed in stripe order:
//     BlastAaLookupFinalize assigns the overflow cursors in ascending cell order), gives its
//     overflow cells their cursors and turns c_p[w] into the position of range p's first append to
//     w: the cell's start (0 in its own entries, or its cursor in the overflow) plus the sum of
//     c_q[w] over q < p.
// (d) Each range makes its appends, s ascending and w in T(s) in list order, at its position for w,
//     advanced by the offsets written.
// Every append to a cell from a range before p comes before range p's appends to it, and within a
// range the appends are made in the reference's order, so each cell holds the reference's entries
// in the reference's order; the ranges write disjoint positions, every position of a non-empty cell
// once. Each source is walked once in (a) and once in (d). The backbone and the overflow are not
// zero-filled first: (b) writes every cell whole (its unused entries zero, as after the C `calloc`)
// and (d) every overflow position.
/// EXPERIMENT (LOSAT_X_BXBUILDPAR): the count pass, the layout and the fill pass of a pooled build
/// (`x_build_direct_in`) on `parts` tasks of the pool it runs on (at most `X_MOST_PARTS`): the
/// table's backbone, pv and overflow vectors (taken from the pool) and its `longest_chain`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[allow(clippy::type_complexity)]
fn x_partitioned_passes(
    parts: usize,
    backbone_size: usize,
    lists: &[std::sync::OnceLock<Box<[u16]>>],
    target_list: &(dyn Fn(usize) -> Box<[u16]> + Sync),
    temps: &mut XBuildTemps,
) -> (Vec<BackboneCell>, Vec<u32>, Vec<i32>, i32) {
    use crate::algorithm::tblastx::lookup::PV_ARRAY_MASK;
    let parts = parts.clamp(1, X_MOST_PARTS);
    // The present source cells (`x_counted_sources`) and the appends' weights.
    let (cells, starts, offsets) = (
        std::mem::take(&mut temps.cells),
        std::mem::take(&mut temps.starts),
        std::mem::take(&mut temps.offsets),
    );
    let weight = |idx: usize, k: usize| {
        (k * lists[idx]
            .get()
            .map_or(X_UNLISTED_LENGTH, |list| list.len())) as u64
    };
    let total: u64 = (0..cells.len())
        .map(|i| weight(cells[i] as usize, (starts[i + 1] - starts[i]) as usize))
        .sum();
    let source = |i: usize| &offsets[starts[i] as usize..starts[i + 1] as usize];
    // The ranges' bounds: range p ends at the first source whose weights so far reach
    // total * (p + 1) / P.
    let mut bounds = [cells.len(); X_MOST_PARTS + 1];
    bounds[0] = 0;
    let (mut cut, mut so_far) = (1usize, 0u64);
    for (i, &cell) in cells.iter().enumerate() {
        if cut == parts {
            break;
        }
        so_far += weight(cell as usize, source(i).len());
        while cut < parts && so_far * parts as u64 >= total * cut as u64 {
            bounds[cut] = i + 1;
            cut += 1;
        }
    }
    let ranges = || bounds[..=parts].windows(2).map(|bound| bound[0]..bound[1]);
    // (a) The ranges' counts per target cell.
    let mut per_range = x_sized(
        std::mem::take(&mut temps.counts),
        parts * backbone_size,
        false,
    );
    x_claimed(
        per_range.chunks_mut(backbone_size).zip(ranges()).collect(),
        &|(counts, range): (&mut [u32], std::ops::Range<usize>)| {
            counts.fill(0);
            for i in range {
                let k = source(i).len() as u32;
                // blast_aalookup.c:490: w = query + offset_list[2];
                let list =
                    lists[cells[i] as usize].get_or_init(|| target_list(source(i)[0] as usize));
                for &w in list.iter() {
                    counts[w as usize] += k;
                }
            }
        },
    );
    let list = |i: usize| -> &[u16] { lists[cells[i] as usize].get().expect("made in (a)") };
    // (b) The stripes' counts, num_used, pv bits, overflow sizes and longest counts.
    let stripe = backbone_size.div_ceil(parts).next_multiple_of(32);
    let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
    let kept = x_pool_lock(&X_BX_TABLES).pop().unwrap_or_default();
    let mut backbone = x_sized(kept.backbone, backbone_size, false);
    let mut pv = x_zeroed(kept.pv, pv_size, false);
    let mut columns: Vec<Vec<&mut [u32]>> = (0..backbone_size.div_ceil(stripe))
        .map(|_| Vec::with_capacity(parts))
        .collect();
    for counts in per_range.chunks_mut(backbone_size) {
        for (column, part) in columns.iter_mut().zip(counts.chunks_mut(stripe)) {
            column.push(part);
        }
    }
    let laid: Vec<_> = x_claimed(
        columns
            .into_iter()
            .zip(backbone.chunks_mut(stripe))
            .zip(pv.chunks_mut(stripe >> PV_ARRAY_BTS))
            .collect(),
        &|((column, backbone), pv): ((Vec<&mut [u32]>, &mut [BackboneCell]), &mut [u32])| {
            let (mut overflow_size, mut longest_chain) = (0usize, 0i32);
            for (local, cell) in backbone.iter_mut().enumerate() {
                let count: u32 = column.iter().map(|counts| counts[local]).sum();
                *cell = BackboneCell {
                    num_used: count as i32,
                    entries: [0; AA_HITS_PER_CELL],
                };
                if count == 0 {
                    continue;
                }
                // pv_set at the stripe's first cell + local, a multiple of 32.
                pv[local >> PV_ARRAY_BTS] |= 1u32 << (local & PV_ARRAY_MASK);
                if count as usize > AA_HITS_PER_CELL {
                    overflow_size += count as usize;
                }
                longest_chain = longest_chain.max(count as i32);
            }
            (overflow_size, longest_chain, column, backbone)
        },
    );
    let overflow_size: usize = laid.iter().map(|part| part.0).sum();
    let longest_chain = laid.iter().map(|part| part.1).max().unwrap_or(0);
    assert!(
        overflow_size < X_IN_OVERFLOW as usize,
        "BLASTX lookup overflow cursor must fit in Int4"
    );
    // (c) The overflow cursors and each range's first position per cell.
    let mut cursor = 0usize;
    let mut placed = Vec::with_capacity(laid.len());
    for (size, _, column, backbone) in laid {
        placed.push((cursor, column, backbone));
        cursor += size;
    }
    x_claimed(placed, &|(mut cursor, mut column, backbone): (
        usize,
        Vec<&mut [u32]>,
        &mut [BackboneCell],
    )| {
        for (local, cell) in backbone.iter_mut().enumerate() {
            let count = cell.num_used as usize;
            if count == 0 {
                continue;
            }
            let mut at = if count > AA_HITS_PER_CELL {
                cell.entries[0] = cursor as i32;
                cursor += count;
                (cursor - count) as u32 | X_IN_OVERFLOW
            } else {
                0
            };
            for counts in column.iter_mut() {
                let n = counts[local];
                counts[local] = at;
                at += n;
            }
        }
    });
    // (d) The appends.
    let mut overflow = x_sized(kept.overflow, overflow_size, true);
    let out = XDisjointOut {
        backbone: backbone.as_mut_ptr(),
        overflow: overflow.as_mut_ptr(),
        #[cfg(debug_assertions)]
        lens: (backbone.len(), overflow.len()),
    };
    x_claimed(
        per_range.chunks_mut(backbone_size).zip(ranges()).collect(),
        &|(at, range): (&mut [u32], std::ops::Range<usize>)| {
            for i in range {
                let offsets = source(i);
                for &w in list(i) {
                    let w = w as usize;
                    // SAFETY: by (c), range p's appends to cell w take the positions from at[w]
                    // (its first, set in (c)) for c_p[w] offsets, inside the cell's count, after the
                    // ranges before p and before the ranges after p; an overflow cell's count lies
                    // in the overflow from its cursor and the cursors give the cells disjoint parts.
                    // So no two tasks write a position, and the backbone and the overflow are read
                    // or written by nothing else until the tasks end.
                    unsafe { out.put(w, at[w], offsets) };
                    at[w] += offsets.len() as u32;
                }
            }
        },
    );
    temps.counts = per_range;
    temps.cells = cells;
    temps.starts = starts;
    temps.offsets = offsets;
    (backbone, pv, overflow, longest_chain)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469
// ```c
//     /* create an empty backbone */
//
//     exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
// ...
//     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
//                                       lookup->charsize, lookup->word_length,
//                                       query, location);
// ...
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
//             sfree(exact_backbone[i]);
//         }
//     }
// ```
// No NCBI counterpart: the BLASTP choice between the reference builder and the in-place builder
// (`build_ncbi_lookup_from_prepared` in backbone.rs carries the NCBI code); it changes no value
// NCBI computes.
/// `build_ncbi_lookup` for BLASTP: prepares the query as `build_ncbi_lookup_for_profile(..,
/// Blosum62, LOOKUP_WORD_LENGTH)` does and builds the table with [`x_build_direct`] (the same
/// table without the chains, so without the chains' memory while it is built), falling back to
/// the reference builder where that declines.
pub(super) fn x_bp_build_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    x_build_direct_or_reference(queries, threshold, check_ideal, |_| {}, x_bp_build_direct)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (see
// `x_build_direct_or_reference`)
// ```c
//     for (i = 0; i < lookup->backbone_size; i++) {
//         if (exact_backbone[i] != NULL) {
//             s_AddWordHits(lookup, matrix, query->sequence,
//                           exact_backbone[i], query_bias, row_max);
// ```
/// EXPERIMENT (LOSAT_X_BPLK16): [`x_bp_build_lookup`] with the width of the table's hits chosen by
/// the build: with `narrow`, `u16` hits where the lookup query is at most 65,535 letters
/// ([`x_bplk16_fits`]) and the in-place builder applies; otherwise `Int4` hits, the table
/// [`x_bp_build_lookup`] builds.
pub(super) fn x_bp_build_lookup_as(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
    narrow: bool,
) -> (XBpLookup, Vec<QueryContext>) {
    let matrix = ScoringMatrix::Blosum62;
    let prepared = x_bp_prepared_query(queries, check_ideal);
    // build_lookup_from_prepared returns the prepared contexts unchanged.
    let contexts = prepared.contexts.clone();
    let prepared = if narrow && x_bplk16_fits(&prepared) {
        let built = X_MEMO.with(|memo| {
            x_build_direct_in::<BackboneCell16>(
                prepared,
                threshold,
                matrix,
                LOOKUP_WORD_LENGTH,
                &mut memo.borrow_mut(),
                None,
                false,
                1,
                crate::seqdb::x_bp_db_keep(),
            )
            .map(|table| table.expect("fresh vectors have the room for any table"))
        });
        match built {
            Ok(table) => {
                x_bplk16_shadow(queries, threshold, check_ideal, &table);
                return (XBpLookup::Narrow(table), contexts);
            }
            Err(prepared) => prepared,
        }
    } else {
        prepared
    };
    // As x_build_direct_or_reference with BLASTP's direct builder.
    match x_bp_build_direct(prepared, threshold, matrix, LOOKUP_WORD_LENGTH) {
        Ok(table) => (XBpLookup::Wide(table), contexts),
        Err(prepared) => {
            let (table, contexts) =
                build_lookup_from_prepared(prepared, threshold, matrix, LOOKUP_WORD_LENGTH);
            (XBpLookup::Wide(table), contexts)
        }
    }
}

/// EXPERIMENT (LOSAT_X_TBXLUT): `true` = TBLASTX's full table built in place by [`x_build_direct`],
/// `false` = off (`build_ncbi_lookup`, the chain builder).
// No NCBI counterpart: chooses between two LOSAT code paths with the same result. A strict-set
// switch (on unless it is `0`; `utils/xswitch.rs`).
pub(super) fn x_tbxlut_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXLUT"))
}

// No NCBI counterpart: records the side the switch took (`LOSAT_X_REACH`).
/// EXPERIMENT (LOSAT_X_TBXLUT): a batch whose unit path builds its group tables directly
/// (PAIRLUT) builds no full table, so the switch is not used for it.
pub(crate) fn x_tbxlut_reach_stub() {
    if x_tbxlut_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXLUT", false);
    }
}

/// EXPERIMENT (LOSAT_X_TBXLUT): `build_ncbi_lookup` for the full table of
/// a TBLASTX query batch (the table of the reference frame loop), with TBLASTX's switch; the same
/// builder as BLASTP's [`x_bp_build_lookup`].
pub(super) fn x_tbx_build_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    karlin_params: &KarlinParams,
    check_ideal: bool,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    x_build_ncbi_lookup_with(
        x_tbxlut_on(),
        queries,
        threshold,
        karlin_params,
        check_ideal,
        |direct| crate::utils::xreach::x_reach!("LOSAT_X_TBXLUT", direct),
        x_build_direct,
    )
}

/// The query of `queries` prepared as `build_ncbi_lookup_for_profile(.., Blosum62,
/// LOOKUP_WORD_LENGTH)` prepares it ([`x_build_direct_or_reference`] and
/// [`XAaSpares::x_build_ahead`]).
fn x_bp_prepared_query(queries: &[Vec<QueryFrame>], check_ideal: bool) -> PreparedLookupQuery {
    let matrix = ScoringMatrix::Blosum62;
    // The arguments of build_ncbi_lookup_for_profile for BLOSUM62.
    let bounds = (-4, 11);
    prepare_lookup_query(
        queries,
        ideal_karlin_params_for_matrix(matrix, bounds),
        &compute_std_aa_composition(),
        LOOKUP_WORD_LENGTH,
        check_ideal,
        matrix,
        bounds,
    )
}

/// The direct builder of [`x_build_direct_or_reference`]: [`x_build_direct`] or BLASTP's
/// [`x_bp_build_direct`].
type XBuildDirect = fn(
    PreparedLookupQuery,
    i32,
    ScoringMatrix,
    usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery>;

/// The body of [`x_tbx_build_lookup`]: with `direct` false it calls `build_ncbi_lookup`; with
/// `direct` true it is [`x_build_direct_or_reference`] with `build_direct`. `reach` is called with
/// `direct` true with whether the direct builder built the table (`LOSAT_X_REACH`).
fn x_build_ncbi_lookup_with(
    direct: bool,
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    karlin_params: &KarlinParams,
    check_ideal: bool,
    reach: fn(bool),
    build_direct: XBuildDirect,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    if !direct {
        return build_ncbi_lookup(queries, threshold, karlin_params, check_ideal);
    }
    x_build_direct_or_reference(queries, threshold, check_ideal, reach, build_direct)
}

/// The table of [`x_bp_build_lookup`] and of [`x_build_ncbi_lookup_with`] with `direct` true: the
/// query is prepared as `build_ncbi_lookup_for_profile(.., Blosum62, LOOKUP_WORD_LENGTH)` prepares
/// it and the table is built with `build_direct` ([`x_build_direct`], or [`x_bp_build_direct`] for
/// BLASTP), falling back to the reference builder where that declines. `reach` is called with
/// whether the direct builder built the table.
fn x_build_direct_or_reference(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
    reach: fn(bool),
    build_direct: XBuildDirect,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let matrix = ScoringMatrix::Blosum62;
    let prepared = x_bp_prepared_query(queries, check_ideal);
    // build_lookup_from_prepared returns the prepared contexts unchanged.
    let contexts = prepared.contexts.clone();
    match build_direct(prepared, threshold, matrix, LOOKUP_WORD_LENGTH) {
        Ok(table) => {
            reach(true);
            (table, contexts)
        }
        Err(prepared) => {
            reach(false);
            build_lookup_from_prepared(prepared, threshold, matrix, LOOKUP_WORD_LENGTH)
        }
    }
}

/// Assert that two tables are equal field by field (the unit tests).
#[cfg(test)]
pub(super) fn x_assert_same_table(direct: &BlastAaLookupTable, reference: &BlastAaLookupTable) {
    x_assert_same_table_as("lookup table:", direct, reference);
}

// No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
/// Whether two tables are equal field by field ([`x_assert_same_table_as`] names the first
/// difference).
#[cfg(test)]
pub(crate) fn x_same_table(a: &BlastAaLookupTable, b: &BlastAaLookupTable) -> bool {
    a.backbone.len() == b.backbone.len()
        && a.backbone
            .iter()
            .zip(&b.backbone)
            .all(|(a, b)| a.num_used == b.num_used && a.entries == b.entries)
        && a.overflow == b.overflow
        && a.pv == b.pv
        && a.frame_bases == b.frame_bases
        && a.num_contexts == b.num_contexts
        && a.query_length == b.query_length
        && a.word_length == b.word_length
        && a.alphabet_size == b.alphabet_size
        && a.charsize == b.charsize
        && a.mask == b.mask
        && a.longest_chain == b.longest_chain
        && a.threshold == b.threshold
        && a.row_max == b.row_max
}

// No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
/// Asserts that two tables are equal field by field; the first difference panics with a message
/// that starts with `prefix`.
#[cfg(test)]
pub(crate) fn x_assert_same_table_as(
    prefix: &str,
    direct: &BlastAaLookupTable,
    reference: &BlastAaLookupTable,
) {
    assert_eq!(
        direct.backbone.len(),
        reference.backbone.len(),
        "{prefix} backbone size differs"
    );
    for (idx, (a, b)) in direct.backbone.iter().zip(&reference.backbone).enumerate() {
        assert!(
            a.num_used == b.num_used && a.entries == b.entries,
            "{prefix} backbone cell {idx} differs ({} {:?} vs {} {:?})",
            a.num_used,
            a.entries,
            b.num_used,
            b.entries
        );
    }
    assert!(
        direct.overflow == reference.overflow,
        "{prefix} overflow array differs"
    );
    assert!(direct.pv == reference.pv, "{prefix} pv differs");
    assert!(
        direct.frame_bases == reference.frame_bases
            && direct.num_contexts == reference.num_contexts
            && direct.query_length == reference.query_length
            && direct.word_length == reference.word_length
            && direct.alphabet_size == reference.alphabet_size
            && direct.charsize == reference.charsize
            && direct.mask == reference.mask
            && direct.longest_chain == reference.longest_chain
            && direct.threshold == reference.threshold
            && direct.row_max == reference.row_max,
        "{prefix} table header differs"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(super) fn frame(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>, frame: i8) -> QueryFrame {
        let mut aa_seq = Vec::with_capacity(seq.len() + 2);
        aa_seq.push(0);
        aa_seq.extend_from_slice(&seq);
        aa_seq.push(0);
        QueryFrame {
            frame,
            aa_len: seq.len(),
            orig_len: seq.len(),
            aa_seq,
            aa_seq_nomask: None,
            seg_masks,
        }
    }

    pub(super) fn residues(seed: &mut u64, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                // NCBISTDAA 1..=27 (every letter of the alphabet but the gap)
                1 + ((*seed >> 33) % 27) as u8
            })
            .collect()
    }

    fn prepared(queries: &[Vec<QueryFrame>], matrix: ScoringMatrix) -> PreparedLookupQuery {
        let bounds = if matrix == ScoringMatrix::Blosum62 {
            (-4, 11)
        } else {
            matrix_score_bounds(matrix)
        };
        prepare_lookup_query(
            queries,
            ideal_karlin_params_for_matrix(matrix, bounds),
            &compute_std_aa_composition(),
            LOOKUP_WORD_LENGTH,
            true,
            matrix,
            bounds,
        )
    }

    fn same(queries: &[Vec<QueryFrame>], threshold: i32, matrix: ScoringMatrix) -> usize {
        let reference =
            build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
        let direct = match x_build_direct(prepared(queries, matrix), threshold, matrix, 3) {
            Ok(table) => table,
            Err(_) => panic!("in-place builder declined"),
        };
        x_assert_same_table(&direct, &reference);
        reference.overflow.len()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-509,580-590
    // ```c
    //     if (lookup->threshold == 0 || score < lookup->threshold) {
    //         for (i = 0; i < offset_list[1]; i++) {
    //             BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
    //                                   lookup->charsize, w,
    //                                   query_bias + offset_list[i + 2]);
    //         }
    // ...
    //         for (i = 0; i < alphabet_size; i++) {
    //             if (score + row[i] >= threshold) {
    //                 subject_word[current_pos] = i;
    //                 for (j = 0; j < offset_list[1]; j++) {
    //                     BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
    //                                           charsize, subject_word,
    //                                           query_bias + offset_list[j + 2]);
    //                 }
    // ```
    #[test]
    fn in_place_table_equals_reference_for_thresholds_masks_and_overflow() {
        let mut seed = 7u64;
        let mut low_complexity = vec![1u8; 300]; // poly-A: cells far above AA_HITS_PER_CELL
        low_complexity.extend(residues(&mut seed, 40));
        let queries = vec![
            vec![
                frame(residues(&mut seed, 700), vec![(10, 30), (400, 401)], 1),
                frame(residues(&mut seed, 2), vec![], 2),
                frame(residues(&mut seed, 3), vec![], 3),
                frame(low_complexity, vec![(250, 299)], -1),
                frame(residues(&mut seed, 650), vec![(0, 0)], -2),
                frame(Vec::new(), vec![], -3),
            ],
            vec![frame(residues(&mut seed, 1200), vec![], 1)],
        ];
        let mut overflow_seen = false;
        for threshold in [0, 1, 11, 12, 13, 19, 30, 1000] {
            overflow_seen |= same(&queries, threshold, ScoringMatrix::Blosum62) > 0;
        }
        assert!(overflow_seen, "the fixture must reach the overflow array");
        // A second query set on the same thread reads the neighbour lists memoised above.
        let other = vec![vec![frame(residues(&mut seed, 900), vec![(5, 50)], 1)]];
        same(&other, 12, ScoringMatrix::Blosum62);
        // Another matrix: the memo key changes and the lists are recomputed.
        same(&other, 12, ScoringMatrix::Blosum45);
        same(&queries, 11, ScoringMatrix::Pam30);
        same(&queries, 12, ScoringMatrix::Blosum62);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:338-343
    // ```c
    //       for (i = 0; i < lookup->backbone_size; i++) {
    //         /* if there are hits there, */
    //         if (lookup->thin_backbone[i] ) {
    //             Int4 * dest = NULL;
    //             /* set the corresponding bit in the pv_array */
    //             PV_SET(pv, i, PV_ARRAY_BTS);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
    // ```c
    //       lookup->thick_backbone =
    //             calloc(lookup->backbone_size, sizeof(AaLookupBackboneCell));
    // ```
    // LOSAT_X_BXLUTPOOL: a table built on vectors another table and another build left (cleared
    // and zero-filled as calloc gives them) and with the neighbour lists shared by all threads is
    // the reference table: after a smaller and a larger table, with two tables alive at once, and
    // after the shared lists change key.
    #[test]
    fn pooled_tables_equal_the_reference_on_kept_vectors() {
        let mut seed = 23u64;
        let mut low_complexity = vec![1u8; 400];
        low_complexity.extend(residues(&mut seed, 60));
        let small = vec![vec![frame(residues(&mut seed, 300), vec![], 1)]];
        let mut invalid = residues(&mut seed, 700);
        for at in [0, 5, 6, 300, 699] {
            invalid[at] = 0x40;
        }
        let large = vec![
            vec![
                frame(residues(&mut seed, 2_500), vec![(100, 140), (143, 200)], 1),
                frame(low_complexity, vec![(350, 399)], -1),
            ],
            vec![frame(residues(&mut seed, 1_800), vec![], 2)],
            vec![frame(invalid, vec![], 3)],
        ];
        let pooled = |queries: &[Vec<QueryFrame>], threshold: i32, matrix: ScoringMatrix| {
            let reference =
                build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
            let table =
                match x_build_direct_pooled(prepared(queries, matrix), threshold, matrix, 3, 1) {
                    Ok(table) => table,
                    Err(_) => panic!("in-place builder declined"),
                };
            x_assert_same_table(&table, &reference);
            table
        };
        for (queries, threshold) in [(&small, 12), (&large, 11), (&small, 13), (&large, 12)] {
            let mut table = pooled(queries, threshold, ScoringMatrix::Blosum62);
            x_give_back_pooled(&mut table);
            assert!(table.backbone.is_empty() && table.overflow.is_empty() && table.pv.is_empty());
        }
        let mut first = pooled(&large, 12, ScoringMatrix::Blosum62);
        let mut second = pooled(&small, 12, ScoringMatrix::Blosum62);
        x_give_back_pooled(&mut first);
        x_give_back_pooled(&mut second);
        for (queries, threshold, matrix) in [
            (&small, 12, ScoringMatrix::Blosum45),
            (&large, 11, ScoringMatrix::Pam30),
            (&large, 12, ScoringMatrix::Blosum62),
        ] {
            let mut table = pooled(queries, threshold, matrix);
            x_give_back_pooled(&mut table);
        }
    }

    // No NCBI counterpart (scheduling only): LOSAT_X_BXBUILDPAR's passes by source ranges give the
    // reference table, with 2 to 8 ranges (the most; more asked for are clamped), with overflow
    // cells (the low-complexity frame), a location shorter than a word (between the first frame's
    // masks) and residues outside the alphabet's bits (the third query: the exact backbone skips
    // their words), a query of fewer source cells than ranges and an all-masked one, thresholds
    // with and without neighbours and other matrices, on vectors other builds left.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn partitioned_tables_equal_the_reference() {
        let mut seed = 31u64;
        let mut low_complexity = vec![1u8; 400];
        low_complexity.extend(residues(&mut seed, 60));
        let small = vec![vec![frame(residues(&mut seed, 300), vec![], 1)]];
        let mut invalid = residues(&mut seed, 700);
        for at in [0, 5, 6, 300, 699] {
            invalid[at] = 0x40;
        }
        let large = vec![
            vec![
                frame(residues(&mut seed, 2_500), vec![(100, 140), (143, 200)], 1),
                frame(low_complexity, vec![(350, 399)], -1),
            ],
            vec![frame(residues(&mut seed, 1_800), vec![], 2)],
            vec![frame(invalid, vec![], 3)],
        ];
        // Fewer source cells than parts, and none at all.
        let tiny = vec![vec![frame(residues(&mut seed, 6), vec![], 1)]];
        let masked = vec![vec![frame(residues(&mut seed, 50), vec![(0, 49)], 1)]];
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            for (queries, threshold, matrix) in [
                (&large, 11, ScoringMatrix::Blosum62),
                (&small, 12, ScoringMatrix::Blosum62),
                (&large, 13, ScoringMatrix::Blosum45),
                (&large, 0, ScoringMatrix::Blosum62),
                (&small, 12, ScoringMatrix::Pam30),
                (&tiny, 11, ScoringMatrix::Blosum62),
                (&masked, 11, ScoringMatrix::Blosum62),
            ] {
                let reference =
                    build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
                let edge = std::ptr::eq(queries, &tiny) || std::ptr::eq(queries, &masked);
                assert!(edge || threshold == 0 || !reference.overflow.is_empty());
                for parts in [2, 3, 4, 7, 8, 1_000] {
                    let mut table = match x_build_direct_pooled(
                        prepared(queries, matrix),
                        threshold,
                        matrix,
                        3,
                        parts,
                    ) {
                        Ok(table) => table,
                        Err(_) => panic!("in-place builder declined"),
                    };
                    x_assert_same_table(&table, &reference);
                    x_give_back_pooled(&mut table);
                }
            }
        });
    }

    #[test]
    fn in_place_table_of_an_all_masked_query_is_empty() {
        let mut seed = 11u64;
        let queries = vec![vec![frame(residues(&mut seed, 50), vec![(0, 49)], 1)]];
        assert_eq!(same(&queries, 12, ScoringMatrix::Blosum62), 0);
        let table = match x_build_direct(
            prepared(&queries, ScoringMatrix::Blosum62),
            12,
            ScoringMatrix::Blosum62,
            3,
        ) {
            Ok(table) => table,
            Err(_) => panic!("in-place builder declined"),
        };
        assert!(table.pv.iter().all(|&w| w == 0));
        assert_eq!(table.longest_chain, 0);
    }

    /// BLASTP's builder against `build_ncbi_lookup`: every table field and every context field.
    /// Returns the overflow length.
    fn same_bp(queries: &[Vec<QueryFrame>], threshold: i32, check_ideal: bool) -> usize {
        let params = blosum62_ideal_karlin_params();
        let (reference, ref_contexts) = build_ncbi_lookup(queries, threshold, &params, check_ideal);
        let (got, contexts) = x_bp_build_lookup(queries, threshold, check_ideal);
        x_assert_same_table(&got, &reference);
        assert_eq!(contexts.len(), ref_contexts.len());
        for (a, b) in contexts.iter().zip(&ref_contexts) {
            assert!(
                a.q_idx == b.q_idx
                    && a.f_idx == b.f_idx
                    && a.frame == b.frame
                    && a.aa_seq == b.aa_seq
                    && a.aa_seq_nomask == b.aa_seq_nomask
                    && a.aa_len == b.aa_len
                    && a.orig_len == b.orig_len
                    && a.frame_base == b.frame_base
                    && a.is_valid == b.is_valid
                    && format!("{:?}", a.karlin_params) == format!("{:?}", b.karlin_params),
                "context {} differs",
                a.q_idx
            );
        }
        reference.overflow.len()
    }

    /// One BLASTP query: a single frame (as blastp builds them), optionally masked.
    fn bp_query(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>, nomask: bool) -> Vec<QueryFrame> {
        let mut f = frame(seq, seg_masks, 1);
        if nomask {
            f.aa_seq_nomask = Some(f.aa_seq.clone());
        }
        vec![f]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
    // ```c
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    //             sfree(exact_backbone[i]);
    //         }
    //     }
    // ```
    #[test]
    fn bplut_table_and_contexts_equal_build_ncbi_lookup() {
        let mut seed = 23u64;
        let mut low_complexity = vec![1u8; 200];
        low_complexity.extend(residues(&mut seed, 30));
        let queries = vec![
            bp_query(residues(&mut seed, 900), vec![(40, 90)], true),
            bp_query(residues(&mut seed, 2), vec![], false),
            bp_query(low_complexity, vec![(150, 199)], true),
            bp_query(residues(&mut seed, 1), vec![], false),
            bp_query(Vec::new(), vec![], false),
            bp_query(residues(&mut seed, 1500), vec![], false),
            bp_query(residues(&mut seed, 60), vec![(0, 59)], true),
        ];
        let mut overflow_seen = false;
        for check_ideal in [false, true] {
            for threshold in [0, 11, 12, 13, 21] {
                overflow_seen |= same_bp(&queries, threshold, check_ideal) > 0;
            }
        }
        assert!(overflow_seen, "the fixture must reach the overflow array");
        // A threshold beyond Int4 (`INT_MIN`, every word a neighbour) on short queries.
        let short: Vec<Vec<QueryFrame>> = (0..3)
            .map(|_| bp_query(residues(&mut seed, 15), vec![], false))
            .collect();
        same_bp(&short, i32::MIN, false);
    }

    /// The query offsets of every cell of `table`, as the scan reads them.
    fn cell_lists(table: &BlastAaLookupTable) -> impl Iterator<Item = &[i32]> {
        table.backbone.iter().map(|cell| {
            let n = cell.num_used as usize;
            if n <= AA_HITS_PER_CELL {
                &cell.entries[..n]
            } else {
                let cursor = cell.entries[0] as usize;
                &table.overflow[cursor..cursor + n]
            }
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-518
    // ```c
    //     if (lookup->threshold == 0 || score < lookup->threshold) {
    //         for (i = 0; i < offset_list[1]; i++) {
    //             BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
    //                                   lookup->charsize, w,
    //                                   query_bias + offset_list[i + 2]);
    //         }
    //     }
    //     ...
    //     if (lookup->threshold == 0)
    //         return;
    // ```
    /// LOSAT_X_BPDIAG (`blastp/x_seed_fused.rs`) tests the hits of a cell in batches, which is
    /// exact because a cell lists a query offset at most once (the word itself is added only when
    /// the neighbour enumeration does not return it). Both BLASTP builders, thresholds without
    /// neighbours (0), with few and with every word as a neighbour, masks and the overflow array.
    #[test]
    fn every_cell_lists_a_query_offset_once() {
        let mut seed = 31u64;
        let mut low_complexity = vec![1u8; 200];
        low_complexity.extend(residues(&mut seed, 30));
        let queries = vec![
            bp_query(residues(&mut seed, 900), vec![(40, 90)], true),
            bp_query(residues(&mut seed, 2), vec![], false),
            bp_query(low_complexity, vec![(150, 199)], true),
            bp_query(residues(&mut seed, 1500), vec![], false),
            bp_query(residues(&mut seed, 60), vec![(0, 59)], true),
        ];
        let short: Vec<Vec<QueryFrame>> = (0..3)
            .map(|_| bp_query(residues(&mut seed, 15), vec![], false))
            .collect();
        let params = blosum62_ideal_karlin_params();
        let mut long_cells = 0;
        for direct in [false, true] {
            for (set, threshold) in [
                (&queries, 0),
                (&queries, 11),
                (&queries, 13),
                (&short, 1),
                (&short, i32::MIN),
            ] {
                let (table, _) = if direct {
                    x_bp_build_lookup(set, threshold, false)
                } else {
                    build_ncbi_lookup(set, threshold, &params, false)
                };
                for hits in cell_lists(&table) {
                    let mut sorted = hits.to_vec();
                    sorted.sort_unstable();
                    assert!(
                        sorted.windows(2).all(|pair| pair[0] != pair[1]),
                        "direct {direct} threshold {threshold}: a cell lists an offset twice: {hits:?}"
                    );
                    long_cells += usize::from(hits.len() > 16);
                }
            }
        }
        assert!(long_cells > 0, "the fixture must give long cells");
    }

    /// BLASTP queries whose table reaches the overflow array (as in
    /// `bplut_table_and_contexts_equal_build_ncbi_lookup`).
    fn bp_queries(seed: &mut u64) -> Vec<Vec<QueryFrame>> {
        let mut low_complexity = vec![1u8; 200];
        low_complexity.extend(residues(seed, 30));
        vec![
            bp_query(residues(seed, 900), vec![(40, 90)], true),
            bp_query(residues(seed, 2), vec![], false),
            bp_query(low_complexity, vec![(150, 199)], true),
            bp_query(Vec::new(), vec![], false),
            bp_query(residues(seed, 1500), vec![], false),
        ]
    }

    /// `x_bp_build_lookup(..)`'s table and contexts are those of `got`: every table field,
    /// and every context field (the Karlin blocks by their bits).
    fn assert_same_bp(got: &(BlastAaLookupTable, Vec<QueryContext>), queries: &[Vec<QueryFrame>]) {
        let params = blosum62_ideal_karlin_params();
        let (want, want_contexts) = x_bp_build_lookup(queries, 11, false);
        x_assert_same_table_as("ahead:", &got.0, &want);
        assert!(x_same_table(&got.0, &want));
        assert_eq!(got.1.len(), want_contexts.len());
        for (a, b) in got.1.iter().zip(&want_contexts) {
            assert!(
                a.q_idx == b.q_idx
                    && a.aa_seq == b.aa_seq
                    && a.frame_base == b.frame_base
                    && a.is_valid == b.is_valid
                    && format!("{:?}", a.karlin_params) == format!("{:?}", b.karlin_params),
                "context {} differs",
                a.q_idx
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
    // ```c
    //       lookup->thick_backbone =
    //             calloc(lookup->backbone_size, sizeof(AaLookupBackboneCell));
    // ...
    //           lookup->overflow =
    //                  calloc(overflow_cells_needed, sizeof(Int4));
    // ```
    /// LOSAT_X_BPLKPAHEAD: the table built ahead on vectors kept by `XAaSpares` (a table given back
    /// with other values, then moved to the build by `x_make_room`) is `x_bp_build_lookup(..)`'s;
    /// `x_make_room` allocates only the vectors short of their room, and keeps the calling
    /// thread's neighbour lists for the builds, which give them back.
    #[test]
    fn aa_spares_build_the_bplut_table_on_kept_vectors() {
        let mut seed = 41u64;
        let queries = bp_queries(&mut seed);
        let params = blosum62_ideal_karlin_params();
        let (first, _) = x_bp_build_lookup(&queries, 11, false);
        let (backbone, pv, overflow) = (first.backbone.len(), first.pv.len(), first.overflow.len());
        assert!(overflow > 0, "the fixture must reach the overflow array");
        assert!(X_MEMO.with(|memo| memo.borrow().is_some()));
        let spares = XAaSpares::default();
        // Fresh spares: all three allocated here, and this thread's lists move to the spares.
        assert_eq!(spares.x_make_room(backbone, pv, overflow, overflow), 3);
        assert!(X_MEMO.with(|memo| memo.borrow().is_none()));
        assert!(x_locked(&spares.memo).is_some());
        let addresses = |v: &XAaVectors| {
            [
                v.backbone.as_ptr() as usize,
                v.overflow.as_ptr() as usize,
                v.pv.as_ptr() as usize,
            ]
        };
        let table_addresses = |t: &BlastAaLookupTable| {
            [
                t.backbone.as_ptr() as usize,
                t.overflow.as_ptr() as usize,
                t.pv.as_ptr() as usize,
            ]
        };
        let kept = addresses(&x_locked(&spares.ahead));
        let built = spares
            .x_build_ahead(&queries, 11, &params, false)
            .expect("room for the table");
        // Built on the kept vectors (the build allocates none of them) with the kept lists, which
        // it gives back (this thread still has none).
        assert_eq!(table_addresses(&built.0), kept);
        assert!(X_MEMO.with(|memo| memo.borrow().is_none()));
        assert!(x_locked(&spares.memo).is_some(), "the lists are given back");
        assert!(
            x_locked(&spares.ahead).is_empty(),
            "the build took the vectors"
        );
        assert_same_bp(&built, &queries);
        // A table given back with other values: its vectors serve the next build.
        let (mut table, _) = built;
        for cell in table.backbone.iter_mut() {
            cell.num_used = 99;
            cell.entries = [7; AA_HITS_PER_CELL];
        }
        table.overflow.fill(-5);
        table.pv.fill(u32::MAX);
        let given = table_addresses(&table);
        spares.give(table);
        assert_eq!(spares.x_make_room(backbone, pv, overflow, overflow), 0);
        assert_eq!(addresses(&x_locked(&spares.ahead)), given);
        assert!(x_locked(&spares.block).is_empty());
        // Its own vectors are kept (those given back meanwhile are dropped); only a vector short of
        // its room is replaced.
        let (other_table, _) = x_bp_build_lookup(&queries, 11, false);
        spares.give(other_table);
        assert_eq!(spares.x_make_room(backbone, pv, overflow, overflow), 0);
        assert_eq!(addresses(&x_locked(&spares.ahead)), given);
        assert!(x_locked(&spares.block).is_empty());
        assert_eq!(
            spares.x_make_room(backbone, pv, overflow * 4, overflow * 8),
            1
        );
        assert!(x_locked(&spares.ahead).overflow.capacity() >= overflow * 8);
        let other = bp_queries(&mut seed);
        let built = spares
            .x_build_ahead(&other, 11, &params, false)
            .expect("room for the table");
        assert_same_bp(&built, &other);
        // Home: the lists go back to a thread without any (not over its own), and with `release`
        // the kept vectors are dropped.
        spares.give(built.0);
        X_MEMO.with(|memo| *memo.borrow_mut() = None);
        spares.x_home(false);
        assert!(X_MEMO.with(|memo| memo.borrow().is_some()));
        assert!(x_locked(&spares.memo).is_none());
        assert!(!x_locked(&spares.block).is_empty());
        spares.x_home(true);
        assert!(x_locked(&spares.block).is_empty() && x_locked(&spares.ahead).is_empty());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:321-336
    // ```c
    //       lookup->thick_backbone =
    //             calloc(lookup->backbone_size, sizeof(AaLookupBackboneCell));
    // ...
    //           lookup->overflow =
    //                  calloc(overflow_cells_needed, sizeof(Int4));
    // ```
    /// LOSAT_X_BPLKREUSE: each block's own build on the block spare (`x_build_on_block_spare`) is
    /// `x_bp_build_lookup(true, ..)`'s table: the first on new vectors; the next ones on the
    /// vectors the block before gave back with other values, the same allocations where they have
    /// the room (a table as large or smaller), grown in place for a larger one. The spare is empty
    /// while a table built on it is in use.
    #[test]
    fn aa_block_spare_builds_each_block_on_the_vectors_given_back() {
        let mut seed = 47u64;
        let queries = bp_queries(&mut seed);
        let small = vec![bp_query(residues(&mut seed, 300), vec![], false)];
        let large = vec![
            bp_query(residues(&mut seed, 4000), vec![(100, 180)], true),
            bp_query(residues(&mut seed, 3000), vec![], false),
        ];
        let params = blosum62_ideal_karlin_params();
        let table_addresses = |t: &BlastAaLookupTable| {
            [
                t.backbone.as_ptr() as usize,
                t.overflow.as_ptr() as usize,
                t.pv.as_ptr() as usize,
            ]
        };
        let spare = XAaSpares::default();
        let first = spare.x_build_on_block_spare(&queries, 11, &params, false);
        assert!(
            !first.0.overflow.is_empty(),
            "the fixture must reach the overflow array"
        );
        assert!(x_locked(&spare.block).is_empty());
        assert_same_bp(&first, &queries);
        // Given back with other values: the next table is built on the same allocations.
        let (mut table, _) = first;
        let overflow = table.overflow.len();
        for cell in table.backbone.iter_mut() {
            cell.num_used = 99;
            cell.entries = [7; AA_HITS_PER_CELL];
        }
        table.overflow.fill(-5);
        table.pv.fill(u32::MAX);
        let given = table_addresses(&table);
        spare.give(table);
        let again = spare.x_build_on_block_spare(&queries, 11, &params, false);
        assert_eq!(table_addresses(&again.0), given);
        assert!(x_locked(&spare.block).is_empty());
        assert_same_bp(&again, &queries);
        spare.give(again.0);
        let smaller = spare.x_build_on_block_spare(&small, 11, &params, false);
        assert_eq!(table_addresses(&smaller.0), given);
        assert_same_bp(&smaller, &small);
        // A larger table: the overflow vector grows to its size.
        spare.give(smaller.0);
        let grown = spare.x_build_on_block_spare(&large, 11, &params, false);
        assert!(grown.0.overflow.len() > overflow);
        assert_same_bp(&grown, &large);
    }

    /// LOSAT_X_BPLKPAHEAD: vectors without the room for the table give no table and are left as
    /// they were; the neighbour lists are kept.
    #[test]
    fn aa_spares_without_the_room_give_no_table() {
        let mut seed = 43u64;
        let queries = bp_queries(&mut seed);
        let params = blosum62_ideal_karlin_params();
        let (table, _) = x_bp_build_lookup(&queries, 11, false);
        let (backbone, pv, overflow) = (table.backbone.len(), table.pv.len(), table.overflow.len());
        for (b, p, o) in [
            (backbone, pv, overflow - 1),
            (backbone - 1, pv, overflow),
            (backbone, pv - 1, overflow),
        ] {
            let spares = XAaSpares::default();
            spares.x_make_room(b, p, o, o);
            let capacities = |spares: &XAaSpares| {
                let ahead = x_locked(&spares.ahead);
                (
                    ahead.backbone.capacity(),
                    ahead.pv.capacity(),
                    ahead.overflow.capacity(),
                )
            };
            let before = capacities(&spares);
            assert!(
                before.0 < backbone || before.1 < pv || before.2 < overflow,
                "{before:?}"
            );
            assert!(spares.x_build_ahead(&queries, 11, &params, false).is_none());
            assert_eq!(capacities(&spares), before);
            assert!(x_locked(&spares.memo).is_some());
        }
    }

    /// The comparison finds one backbone cell, one overflow entry, one pv bit and one header field
    /// that differ, and names the first difference after its prefix.
    #[test]
    fn same_table_finds_one_cell_entry_bit_or_field() {
        let mut seed = 47u64;
        let queries = bp_queries(&mut seed);
        let params = blosum62_ideal_karlin_params();
        let build = || x_bp_build_lookup(&queries, 11, false).0;
        let want = build();
        assert!(x_same_table(&build(), &want));
        let cell = want
            .backbone
            .iter()
            .position(|cell| cell.num_used > 0 && cell.num_used as usize <= AA_HITS_PER_CELL)
            .expect("a cell with entries");
        let changes: [(&str, &dyn Fn(&mut BlastAaLookupTable)); 4] = [
            ("backbone cell", &|t| t.backbone[cell].entries[0] += 1),
            ("overflow array", &|t| t.overflow[0] += 1),
            ("pv", &|t| t.pv[0] ^= 1 << 31),
            ("table header", &|t| t.longest_chain += 1),
        ];
        for (what, change) in changes {
            let mut got = build();
            change(&mut got);
            assert!(!x_same_table(&got, &want), "{what}");
            let message = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                x_assert_same_table_as("PREFIX:", &got, &want)
            }))
            .expect_err(what);
            let message = message
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| message.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default();
            assert!(message.starts_with(&format!("PREFIX: {what}")), "{message}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-358 (quoted at
    // `x_build_direct`)
    // LOSAT_X_BPLK16: the table of 16-bit hits holds, cell by cell, the hits of the table of
    // 32-bit hits of the same block in the same order (the comparison of LOSAT_X_BPLK16SHADOW,
    // against the reference builder), for short queries, a low-complexity query whose cells are far
    // above AA_HITS_PER_CELL (the overflow array, with its cursor split over two entries), several
    // thresholds, and the block that is just at the 65,535 letters of the choice (16-bit) and one
    // letter above it (32-bit).
    #[test]
    fn table_of_u16_hits_equals_the_table_of_i32_hits() {
        let mut seed = 16u64;
        let mut low_complexity = vec![1u8; 300];
        low_complexity.extend(residues(&mut seed, 40));
        let queries = vec![
            vec![
                frame(residues(&mut seed, 700), vec![(10, 30), (400, 401)], 1),
                frame(residues(&mut seed, 2), vec![], 2),
                frame(low_complexity, vec![(250, 299)], -1),
                frame(Vec::new(), vec![], -3),
            ],
            vec![frame(residues(&mut seed, 1200), vec![], 1)],
        ];
        let mut overflow_seen = false;
        for threshold in [0, 11, 13, 1000] {
            let (built, contexts) = x_bp_build_lookup_as(&queries, threshold, false, true);
            let XBpLookup::Narrow(narrow) = built else {
                panic!("a small block takes 16-bit hits (threshold {threshold})");
            };
            let (wide, wide_contexts) = x_bp_build_lookup(&queries, threshold, false);
            assert_eq!(contexts.len(), wide_contexts.len());
            x_assert_narrow_same_as_wide(&format!("threshold {threshold}:"), &narrow, &wide);
            overflow_seen |= !wide.overflow.is_empty();
            assert!(wide.longest_chain as usize > AA_HITS_PER_CELL);
        }
        assert!(overflow_seen, "the fixture must reach the overflow array");
        // Not asked for: the table of 32-bit hits.
        assert!(matches!(
            x_bp_build_lookup_as(&queries, 11, false, false).0,
            XBpLookup::Wide(_)
        ));
        // The choice at 65,535 letters: the lookup query's length is the last context's frame_base
        // plus its aa_len.
        let base = x_bp_prepared_query(&[vec![frame(residues(&mut seed, 50), vec![], 1)]], false)
            .contexts[0]
            .frame_base as usize;
        for (length, narrow) in [(65_535 - base, true), (65_536 - base, false)] {
            let block = vec![vec![frame(residues(&mut seed, length), vec![], 1)]];
            assert_eq!(
                x_bplk16_fits(&x_bp_prepared_query(&block, false)),
                narrow,
                "{length} letters"
            );
            match x_bp_build_lookup_as(&block, 11, false, true).0 {
                XBpLookup::Narrow(table) => {
                    assert!(narrow, "{length} letters took 16-bit hits");
                    assert_eq!(table.query_length, 65_535);
                    let reference = build_lookup_from_prepared(
                        x_bp_prepared_query(&block, false),
                        11,
                        ScoringMatrix::Blosum62,
                        LOOKUP_WORD_LENGTH,
                    )
                    .0;
                    x_assert_narrow_same_as_wide("65,535 letters:", &table, &reference);
                }
                XBpLookup::Wide(table) => {
                    assert!(!narrow, "{length} letters took 32-bit hits");
                    assert_eq!(table.query_length, 65_536);
                }
            }
        }
    }

    // LOSAT_X_BPLK16: the vectors of a table of 16-bit hits given back (`give16`) serve the next
    // block's build on the same spares (`x_build_on_block_spare16`, grown in place where short),
    // and the vectors made ready for a build ahead (`x_make_room` with the narrow width, the
    // padding counted in the sizes the table before reports) serve `x_build_ahead16`; a spare of the
    // other width is not reused (the room is made for the width asked for) and the tables are the
    // tables of 32-bit hits.
    #[test]
    fn spares_serve_the_builds_of_both_widths() {
        let mut seed = 160u64;
        let mut low_complexity = vec![1u8; 200];
        low_complexity.extend(residues(&mut seed, 30));
        let blocks: Vec<Vec<Vec<QueryFrame>>> = vec![
            vec![vec![frame(residues(&mut seed, 900), vec![], 1)]],
            vec![
                vec![frame(low_complexity, vec![], 1)],
                vec![frame(residues(&mut seed, 1500), vec![(7, 70)], 1)],
            ],
            vec![vec![frame(residues(&mut seed, 300), vec![], 1)]],
        ];
        let params = blosum62_ideal_karlin_params();
        let spares = XAaSpares::default();
        let narrow_of = |built: XBpLookup| match built {
            XBpLookup::Narrow(table) => table,
            XBpLookup::Wide(_) => panic!("a small block takes 16-bit hits"),
        };
        // The builds on the block spare, each on the vectors the block before gave back.
        for block in &blocks {
            let table = narrow_of(spares.x_build_on_block_spare16(block, 11, false).0);
            let wide = x_bp_build_lookup(block, 11, false).0;
            x_assert_narrow_same_as_wide("on the block spare:", &table, &wide);
            spares.give16(table);
        }
        // The build ahead on the room made for the next table, from the table before it.
        let first = narrow_of(spares.x_build_on_block_spare16(&blocks[0], 11, false).0);
        spares.x_make_room16(first.backbone.len(), first.pv.len(), 1 << 16, 1 << 17);
        spares.give16(first);
        for block in &blocks {
            let (built, _) = spares
                .x_build_ahead16(block, 11, false)
                .expect("the room made holds the table");
            let table = narrow_of(built);
            let wide = x_bp_build_lookup(block, 11, false).0;
            x_assert_narrow_same_as_wide("built ahead:", &table, &wide);
            spares.x_make_room16(table.backbone.len(), table.pv.len(), 1 << 16, 1 << 17);
        }
        // A table of 32-bit hits built ahead takes room of its own width, not the 16-bit vectors.
        assert!(spares
            .x_build_ahead(&blocks[0], 11, &params, false)
            .is_none());
        let wide = x_bp_build_lookup(&blocks[0], 11, false).0;
        spares.x_make_room(wide.backbone.len(), wide.pv.len(), 1 << 16, 1 << 17);
        let (table, _) = spares
            .x_build_ahead(&blocks[0], 11, &params, false)
            .expect("the room made holds the table");
        x_assert_same_table(&table, &wide);
    }
}
