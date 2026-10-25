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
//! function, and every other field of the table is computed by the reference expressions.
//! In diagnostics mode the reference builder runs (its prints read the chains).

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

/// No NCBI counterpart: a build's temporary vectors, kept between pooled builds.
#[derive(Default)]
struct XBuildTemps {
    counts: Vec<u32>,
    base: Vec<u32>,
    fill: Vec<u32>,
    chains: Option<LookupBackboneChains>,
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
/// [`x_build_direct`] with the shared neighbour lists and kept vectors (`LOSAT_X_BXLUTPOOL`).
pub(super) fn x_build_direct_pooled(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
    let mut unused = None;
    x_build_direct_in(
        prepared,
        threshold,
        matrix,
        word_length,
        &mut unused,
        None,
        true,
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
}

impl XAaVectors {
    // No NCBI counterpart: buffer reuse only (above).
    /// The three vectors, cleared and zero-filled to `backbone_size`, `pv_size` and `overflow_size`
    /// entries, where each has the capacity (nothing is allocated); else `None`, nothing taken.
    fn x_take_zeroed(
        &mut self,
        backbone_size: usize,
        pv_size: usize,
        overflow_size: usize,
    ) -> Option<(Vec<BackboneCell>, Vec<u32>, Vec<i32>)> {
        if self.backbone.capacity() < backbone_size
            || self.pv.capacity() < pv_size
            || self.overflow.capacity() < overflow_size
        {
            return None;
        }
        let mut backbone = std::mem::take(&mut self.backbone);
        backbone.clear();
        backbone.resize(backbone_size, BackboneCell::default());
        let mut pv = std::mem::take(&mut self.pv);
        pv.clear();
        pv.resize(pv_size, 0);
        let mut overflow = std::mem::take(&mut self.overflow);
        overflow.clear();
        overflow.resize(overflow_size, 0);
        Some((backbone, pv, overflow))
    }

    // No NCBI counterpart: buffer reuse only.
    fn is_empty(&self) -> bool {
        self.backbone.capacity() == 0 && self.overflow.capacity() == 0 && self.pv.capacity() == 0
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
    /// the table is dropped.
    pub(crate) fn give(&self, table: BlastAaLookupTable) {
        *x_locked(&self.block) = XAaVectors {
            backbone: table.backbone,
            overflow: table.overflow,
            pv: table.pv,
        };
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
        let mut kept = std::mem::take(&mut *x_locked(&self.ahead));
        let given_back = std::mem::take(&mut *x_locked(&self.block));
        if kept.is_empty() {
            kept = given_back;
        }
        let mut allocated = 0;
        if kept.backbone.capacity() < backbone_size {
            kept.backbone = Vec::with_capacity(backbone_size);
            allocated += 1;
        }
        if kept.pv.capacity() < pv_size {
            kept.pv = Vec::with_capacity(pv_size);
            allocated += 1;
        }
        if kept.overflow.capacity() < overflow {
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
        let built = x_build_direct_in(
            prepared,
            threshold,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
            &mut memo,
            Some(&mut room),
            false,
        );
        *x_locked(&self.memo) = memo;
        *x_locked(&self.ahead) = room;
        match built {
            Ok(Some(table)) => Some((table, contexts)),
            Ok(None) | Err(_) => None,
        }
    }
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
        x_build_direct_in(
            prepared,
            threshold,
            matrix,
            word_length,
            &mut memo.borrow_mut(),
            None,
            false,
        )
        .map(|table| table.expect("fresh vectors have the room for any table"))
    })
}

/// [`x_build_direct`] with the neighbour lists of `memo` (the calling thread's, or those an
/// `XAaSpares` keeps), and, with `room` (LOSAT_X_BPLKPAHEAD), the table's backbone, overflow and pv
/// vectors taken from `room` instead of allocated: `Ok(None)` (nothing taken) where one of them
/// lacks the capacity.
fn x_build_direct_in(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    memo: &mut Option<XNeighborMemo>,
    room: Option<&mut XAaVectors>,
    pooled: bool,
) -> std::result::Result<Option<BlastAaLookupTable>, PreparedLookupQuery> {
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
    // The reference exact backbone (blast_aalookup.c:448-456).
    let mut temps = if pooled {
        x_pool_lock(&X_BX_TEMPS).pop().unwrap_or_default()
    } else {
        XBuildTemps::default()
    };
    let mut exact_backbone = match temps.chains.take() {
        Some(mut chains) => {
            chains.x_reset(backbone_size);
            chains
        }
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

    let key = (matrix, threshold, word_length, alphabet_size);
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
    let (mut backbone, mut pv, mut overflow) = match room {
        None if pooled => {
            let kept = x_pool_lock(&X_BX_TABLES).pop().unwrap_or_default();
            (
                x_zeroed(kept.backbone, backbone_size, false),
                x_zeroed(kept.pv, pv_size, false),
                x_zeroed(kept.overflow, overflow_size, true),
            )
        }
        None => (
            vec![BackboneCell::default(); backbone_size],
            vec![0u32; pv_size],
            vec![0i32; overflow_size],
        ),
        Some(room) => match room.x_take_zeroed(backbone_size, pv_size, overflow_size) {
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
        backbone[idx].num_used = count as i32;
        if count as usize > AA_HITS_PER_CELL {
            backbone[idx].entries[0] = overflow_cursor as i32;
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
                backbone[w].entries[at..at + offsets.len()].copy_from_slice(offsets);
            } else {
                let start = base[w] as usize + at;
                overflow[start..start + offsets.len()].copy_from_slice(offsets);
            }
            fill[w] += offsets.len() as u32;
        }
    }

    if pooled {
        x_pool_lock(&X_BX_TEMPS).push(XBuildTemps {
            counts,
            base,
            fill,
            chains: Some(exact_backbone),
        });
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
    x_build_direct_or_reference(queries, threshold, check_ideal, |_| {})
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

/// The body of [`x_tbx_build_lookup`]: with `direct` false it calls `build_ncbi_lookup`; with
/// `direct` true it is [`x_build_direct_or_reference`]. `reach` is called with `direct` true with
/// whether the direct builder built the table (`LOSAT_X_REACH`).
fn x_build_ncbi_lookup_with(
    direct: bool,
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    karlin_params: &KarlinParams,
    check_ideal: bool,
    reach: fn(bool),
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    if !direct {
        return build_ncbi_lookup(queries, threshold, karlin_params, check_ideal);
    }
    x_build_direct_or_reference(queries, threshold, check_ideal, reach)
}

/// The table of [`x_bp_build_lookup`] and of [`x_build_ncbi_lookup_with`] with `direct` true: the
/// query is prepared as `build_ncbi_lookup_for_profile(.., Blosum62, LOOKUP_WORD_LENGTH)` prepares
/// it and the table is built with [`x_build_direct`], falling back to the reference builder where
/// that declines. `reach` is called with whether the direct builder built the table.
fn x_build_direct_or_reference(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    check_ideal: bool,
    reach: fn(bool),
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let matrix = ScoringMatrix::Blosum62;
    let prepared = x_bp_prepared_query(queries, check_ideal);
    // build_lookup_from_prepared returns the prepared contexts unchanged.
    let contexts = prepared.contexts.clone();
    match x_build_direct(prepared, threshold, matrix, LOOKUP_WORD_LENGTH) {
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
        let large = vec![
            vec![
                frame(residues(&mut seed, 2_500), vec![(100, 140)], 1),
                frame(low_complexity, vec![(350, 399)], -1),
            ],
            vec![frame(residues(&mut seed, 1_800), vec![], 2)],
        ];
        let pooled = |queries: &[Vec<QueryFrame>], threshold: i32, matrix: ScoringMatrix| {
            let reference =
                build_lookup_from_prepared(prepared(queries, matrix), threshold, matrix, 3).0;
            let table = match x_build_direct_pooled(prepared(queries, matrix), threshold, matrix, 3)
            {
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
}
