//! SC10b PR-3: the group tables of the TBLASTX unit path (`x_pair_par.rs`) for words whose backbone
//! index does not fit 16 bits (`-word_size 4`: 20 bits), built for all the batch's tables at once,
//! with 32-bit indices, by counting and filling on the search pool. A child module of
//! `x_pair_lut.rs`, whose `x_build_group` (16-bit neighbour lists) keeps building the tables of
//! words of 2 and 3 letters.
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
//! # What is built
//!
//! For lists of lookup locations L_0, .., L_{n-1} (the groups of `x_pair_par.rs::group_starts`, as
//! `XPairLutQuery::x_group_locations` cuts the batch's locations), table t is the table
//! `XPairLutQuery::x_build_group` builds from L_t: the full table's chains filtered to the query
//! offsets of L_t, order kept, in `BlastAaLookupFinalize`'s layout, compacted where the branch-free
//! scan takes the table (x_pair_lut.rs, items 1-6). With E_t(s) the offsets of the source word s in
//! L_t in the order `BlastLookupIndexQueryExactMatches` adds them, and T(s) the cells
//! `s_AddWordHits` appends E(s) to, the chain of cell w in table t is: for s ascending with w in
//! T(s), E_t(s) in order (x_pair_lut.rs, items 1-5).
//!
//! # The steps
//!
//! 1. Sources. A stable counting sort of the words of L_0, then L_1, .., by source cell s (the loop
//!    of `blast_lookup_index_query_exact_matches`, `x_wide_exact_words`): the sources with a word,
//!    ascending, and each source's offsets ordered by (table, order of the table's loop). Table
//!    t's part of source s is therefore E_t(s), in its order (x_pair_lut.rs item 2: each location
//!    is indexed on its own, so the loop over L_t adds exactly the words of L_t, in L_t's order).
//! 2. Lists. T(s) of each source, once for all tables, on the pool over ranges of sources: the
//!    statements of `x_target_list` (the exact cell first when s_AddWordHits adds it, then the
//!    neighbours in the order of the `s_AddWordHitsCore` recursion) with 32-bit indices. T(s) reads
//!    only the residues of s, the matrix, the threshold and the word length (x_lut_direct.rs header),
//!    so it is read at s's first offset, as the reference reads it (blast_aalookup.c:490).
//! 3. Each table on its own, its sources cut into P contiguous ranges (of about equal work):
//!    (a) each range counts its appends per cell, c_p[w] = sum of |E_t(s)| over its s with w in
//!    T(s); (b) each stripe of cells (a multiple of 32 cells) sums them, count[w] = sum over p of
//!    c_p[w], the reference's chain length, with its cells' overflow size, longest count and used
//!    cells; (c) each stripe, after the overflow and the kept cells of the stripes before it (cells
//!    ascending: `BlastAaLookupFinalize` assigns the overflow cursors in ascending cell order, and a
//!    compacted table keeps its cells in index order), writes its cells (`num_used`, the pv bit,
//!    the cursor of an overflow cell; a cell without entries all zero, as after the C `calloc`) and
//!    turns c_p[w] into range p's first position in cell w: the cell's start plus the counts of
//!    the ranges before p; (d) each range appends, for s ascending and w in T(s) in list order,
//!    E_t(s) at its position for w, advanced by the offsets written.
//!
//! # Why the tables are those of `x_build_group`
//!
//! In each cell the appends of range p come after those of every range before it and before those
//! of every range after it, and inside a range they come in source order, so cell w of table t
//! holds, for s ascending with w in T(s), E_t(s) in order: the chain of x_pair_lut.rs item 5. The
//! ranges write disjoint positions, every position of a cell with entries once. `num_used`, the pv
//! bits, the inline entries or the overflow cursor, the overflow array and `longest_chain` are
//! `BlastAaLookupFinalize`'s function of these chains; the compaction keeps the cells with entries
//! in index order (`x_scan_bf::x_compact`) with the rank of the presence vector (`x_scan_bf::x_rank`);
//! `q_min` and `q_max` are the least and greatest offset of a source word with a target, as in
//! `x_build_group`; the header fields are the expressions of `x_build_group`. Only the scheduling
//! (ranges, stripes, the pool's threads) differs, and it changes no value.
//!
//! Memory: the tables, the lists T(s) of the batch's sources (freed when the tables are built; no
//! memo is kept across batches: at 4 letters the lists of a batch take hundreds of megabytes) and P
//! counts per backbone cell (P at most `X_WIDE_MOST_PARTS`, 3.6 MB each at 4 letters). No chain is
//! made: the reference builder holds the chains and the table at once.

use super::*;

/// No NCBI counterpart: the letter width of the lookup alphabet (`charsize`, ilog2(28) + 1), as a
/// constant.
const X_WIDE_CHARSIZE: usize = LOOKUP_ALPHABET_SIZE.ilog2() as usize + 1;

/// No NCBI counterpart: the widest backbone index the joint builder takes, in bits. Each source
/// range counts per backbone cell (`X_WIDE_MOST_PARTS` arrays): 3.6 MB each at 20 bits (words of
/// 4 letters, the longest TBLASTX takes: `Word-size must be less than 6`, checked at 5).
const X_WIDE_INDEX_BITS: usize = 20;

/// No NCBI counterpart: the most source ranges of a table (each holds a count per backbone cell).
const X_WIDE_MOST_PARTS: usize = 16;

/// No NCBI counterpart: the least work of a table's source range, per backbone cell (`x_wide_table`).
/// A range zeroes, sums and lays out a count per backbone cell (passes of about a nanosecond a
/// cell), against several nanoseconds an append, so the ranges of a small table (short queries,
/// many groups) cost more than they split. Measured (60 queries of 1.5-3 kb against a 0.6 Mb
/// subject, 3 rounds): with a range per thread (up to 16) 0.80, 1.02 and 1.02 s at 8, 16 and 32
/// threads (user 4.1, 8.0, 9.5 s); with one append a cell 0.74, 0.85 and 0.89 s (user 3.9, 5.1,
/// 6.4 s); with two 0.82, 1.23 and 1.30 s (too few ranges). Two 0.6 Mb genomes at 8 threads: 1.98
/// s against 2.15 (their tables take every range either way).
const X_WIDE_RANGE_WORK_PER_CELL: u64 = 1;

/// No NCBI counterpart: ranges of sources per pool thread for the lists T(s) (they differ in
/// length; the pool's threads claim the ranges in turn).
const X_WIDE_LIST_PARTS_PER_THREAD: usize = 4;

/// No NCBI counterpart: the bit of a range's position that says the position is in the overflow
/// array; without it the position is `(cell slot << 2) | entry`.
const X_WIDE_IN_OVERFLOW: u32 = 1 << 31;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:103-105 (ComputeTableIndex,
// cited at `x_direct_covers`)
/// Whether the joint builder takes words of `word_length` letters: their backbone index, 5 bits a
/// letter, fits `X_WIDE_INDEX_BITS` (words of up to 4 letters).
#[inline(always)]
pub(crate) const fn x_wide_covers(word_length: usize) -> bool {
    word_length >= 1 && word_length * X_WIDE_CHARSIZE <= X_WIDE_INDEX_BITS
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:44-46
// ```c
//     /* compute the backbone cell to update */
//
//     index = ComputeTableIndex(wordsize, charsize, seq);
// ```
/// The backbone index of the `word_length`-letter word at `word` (32 bits; `x_wide_covers`).
#[inline(always)]
fn x_wide_index(word: &[u8], word_length: usize) -> u32 {
    let mut index = 0u32;
    for &residue in &word[..word_length] {
        index = (index << X_WIDE_CHARSIZE) | u32::from(residue);
    }
    index
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:79-132 (the loop of
// `blast_lookup_index_query_exact_matches`, backbone.rs)
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
// The statements of `x_for_each_exact_word` (x_lut_direct.rs) with 32-bit indices: where the loop
// calls BlastLookupAddWordHit, `visit` gets the word's backbone index and the query offset.
/// `visit(index, offset)` for each word the exact backbone's build adds, in the order it adds them.
fn x_wide_exact_words(
    query: &[u8],
    locations: &[(i32, i32)],
    word_length: usize,
    mut visit: impl FnMut(u32, i32),
) {
    let invalid_mask: u8 = 0xffu8 << (X_WIDE_CHARSIZE as u32);
    let offset =
        |start: usize| i32::try_from(start).expect("NCBI BLAST query offset must fit Int4");
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
                visit(x_wide_index(&query[start..], word_length), offset(start));
            }
            if query[seq] & invalid_mask != 0 {
                word_target = seq + word_length + 1;
            }
        }
        if to + 1 >= word_target {
            let start = to + 1 - word_length;
            visit(x_wide_index(&query[start..], word_length), offset(start));
        }
    }
}

/// No NCBI counterpart: the matrix rows and row maxima the lists read (`info->matrix`,
/// `info->row_max` of s_AddWordHits), as arrays of the alphabet's size.
struct XWideScores {
    rows: [[i32; LOOKUP_ALPHABET_SIZE]; LOOKUP_ALPHABET_SIZE],
    row_max: [i32; LOOKUP_ALPHABET_SIZE],
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
// The recursion of `x_record_neighbors` with the same tests in the same order; the letters chosen
// so far are carried as their index (`prefix`, ComputeTableIndex's shifts over `subject_word`), so a
// complete word's backbone index is `(prefix << charsize) | i`.
#[allow(clippy::too_many_arguments)]
fn x_wide_neighbors(
    query_word: &[u8],
    word_length: usize,
    scores: &XWideScores,
    threshold: i32,
    score: i32,
    current_pos: usize,
    prefix: u32,
    out: &mut Vec<u32>,
) {
    let query_residue = query_word[current_pos] as usize;
    let score = score - scores.row_max[query_residue];
    let row = &scores.rows[query_residue];
    if current_pos == word_length - 1 {
        for (residue, &residue_score) in row.iter().enumerate() {
            if score + residue_score >= threshold {
                out.push((prefix << X_WIDE_CHARSIZE) | residue as u32);
            }
        }
        return;
    }
    for (residue, &residue_score) in row.iter().enumerate() {
        if score + residue_score >= threshold {
            x_wide_neighbors(
                query_word,
                word_length,
                scores,
                threshold,
                score + residue_score,
                current_pos + 1,
                (prefix << X_WIDE_CHARSIZE) | residue as u32,
                out,
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:490-496,504-509,518-519,539-543
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
//     if (lookup->threshold == 0)
//         return;
// ...
//     score = row_max[w[0]];
//     for (i = 1; i < lookup->word_length; i++)
//         score += row_max[w[i]];
//
//     s_AddWordHitsCore(&info, score, 0);
// ```
// The statements of `x_target_list`, appending T(s) to `out`.
/// Appends T(s) of the word at `query_word` to `out`.
fn x_wide_targets(
    query_word: &[u8],
    word_length: usize,
    scores: &XWideScores,
    threshold: i32,
    out: &mut Vec<u32>,
) {
    let first = query_word[0] as usize;
    let mut self_score = scores.rows[first][first];
    for &residue in &query_word[1..word_length] {
        let residue = residue as usize;
        self_score += scores.rows[residue][residue];
    }
    if threshold == 0 || self_score < threshold {
        out.push(x_wide_index(query_word, word_length));
    }
    if threshold != 0 {
        let mut max_score = scores.row_max[first];
        for &residue in &query_word[1..word_length] {
            max_score += scores.row_max[residue as usize];
        }
        x_wide_neighbors(
            query_word,
            word_length,
            scores,
            threshold,
            max_score,
            0,
            0,
            out,
        );
    }
}

/// No NCBI counterpart: table t's part of one source (step 1): source `i` (its index among the
/// sources) and its offsets E_t(s), `offsets[from..to]`.
#[derive(Clone, Copy)]
struct XWideSource {
    i: u32,
    from: u32,
    to: u32,
}

/// No NCBI counterpart: the sources of the batch's tables (step 1).
struct XWideSources {
    /// The cell of each source, ascending.
    cells: Vec<u32>,
    /// The offsets of the sources, by source, then by table, then in the order of the table's
    /// loop.
    offsets: Vec<i32>,
    /// Source i's first offset (`offsets[first[i]]`), where its word is read.
    first: Vec<u32>,
    /// Each table's sources, in source order.
    by_table: Vec<Vec<XWideSource>>,
}

impl XWideSources {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-456 (cited at
    // `XPairLutQuery::x_build_group`): the exact backbone of each table's locations. No NCBI
    // counterpart for the form (a counting sort, as `x_build_lean_in` makes it).
    /// The sources of the tables of `tables` (step 1).
    fn new(
        query: &[u8],
        tables: &[Vec<(i32, i32)>],
        word_length: usize,
        backbone_size: usize,
    ) -> Self {
        assert!(
            tables.len() <= usize::from(u16::MAX) + 1,
            "joint lookup build: too many tables"
        );
        let mut cursor: Vec<u32> = vec![0; backbone_size];
        let mut words = 0usize;
        for locations in tables {
            x_wide_exact_words(query, locations, word_length, |idx, _| {
                cursor[idx as usize] += 1;
                words += 1;
            });
        }
        assert!(
            words < u32::MAX as usize,
            "joint lookup build: query words must fit 32 bits"
        );
        let mut cells: Vec<u32> = Vec::new();
        let mut starts: Vec<u32> = Vec::new();
        let mut running = 0u32;
        for (idx, slot) in cursor.iter_mut().enumerate() {
            if *slot != 0 {
                cells.push(idx as u32);
                starts.push(running);
                running += std::mem::replace(slot, running);
            }
        }
        starts.push(running);
        let mut offsets: Vec<i32> = vec![0; words];
        let mut table_of: Vec<u16> = vec![0; words];
        for (t, locations) in tables.iter().enumerate() {
            x_wide_exact_words(query, locations, word_length, |idx, offset| {
                let at = cursor[idx as usize] as usize;
                offsets[at] = offset;
                table_of[at] = t as u16;
                cursor[idx as usize] += 1;
            });
        }
        drop(cursor);
        // Each source's offsets come table by table (the tables were walked in order), so table
        // t's part of source i is one run.
        let mut by_table: Vec<Vec<XWideSource>> = vec![Vec::new(); tables.len()];
        for i in 0..cells.len() {
            let (start, end) = (starts[i], starts[i + 1]);
            let mut from = start;
            while from < end {
                let t = table_of[from as usize];
                let mut to = from + 1;
                while to < end && table_of[to as usize] == t {
                    to += 1;
                }
                by_table[usize::from(t)].push(XWideSource {
                    i: i as u32,
                    from,
                    to,
                });
                from = to;
            }
        }
        starts.pop();
        Self {
            cells,
            offsets,
            first: starts,
            by_table,
        }
    }

    /// The offsets of `source` (E_t(s) of its table).
    #[inline(always)]
    fn offsets(&self, source: XWideSource) -> &[i32] {
        &self.offsets[source.from as usize..source.to as usize]
    }
}

/// No NCBI counterpart: the lists T(s) of the sources (step 2), made by ranges of sources: range
/// r holds the sources `bounds[r]..bounds[r + 1]`, source i's list
/// `targets[ends[j - 1]..ends[j]]` of its range (j = i - bounds[r], `ends[-1]` = 0).
struct XWideLists {
    bounds: Vec<usize>,
    parts: Vec<(Vec<u32>, Vec<usize>)>,
}

impl XWideLists {
    /// T(s) of source `i`.
    #[inline(always)]
    fn list(&self, i: u32) -> &[u32] {
        let i = i as usize;
        let r = self.bounds.partition_point(|&b| b <= i) - 1;
        let (targets, ends) = &self.parts[r];
        let j = i - self.bounds[r];
        let from = if j == 0 { 0 } else { ends[j - 1] };
        &targets[from..ends[j]]
    }
}

// No NCBI counterpart (scheduling only): NCBI builds a table on one thread. `f` on each of `parts`,
// the results in part order: on tasks of `pool` (`x_claimed`) when it has two threads or more, else
// in order on this thread.
fn x_wide_parts<T: Send, R: Send>(
    pool: &crate::utils::threading::SearchPool<'_>,
    parts: Vec<T>,
    f: &(dyn Fn(T) -> R + Sync),
) -> Vec<R> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parts.len() > 1 && x_wide_threads(pool) > 1 {
        return pool.install(|| x_claimed(parts, f));
    }
    let _ = pool;
    parts.into_iter().map(f).collect()
}

// No NCBI counterpart (scheduling only): the threads the build's tasks run on.
fn x_wide_threads(pool: &crate::utils::threading::SearchPool<'_>) -> usize {
    if cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )) && pool.enabled()
    {
        pool.threads().max(1)
    } else {
        1
    }
}

/// No NCBI counterpart: where the ranges of step (d) write their appends (a table's backbone and
/// overflow, written by the ranges side by side at disjoint positions).
struct XWideOut {
    backbone: *mut BackboneCell,
    overflow: *mut i32,
    #[cfg(debug_assertions)]
    lens: (usize, usize),
}

// SAFETY: the ranges write through the pointers only, each at positions no other range writes
// (`XWideOut::put`), while nothing else reads or writes the vectors.
unsafe impl Sync for XWideOut {}

impl XWideOut {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:74-76 (an append at the
    // chain's end, cited in the module comment); the copy of BlastAaLookupFinalize
    // (blast_aalookup.c:351-359, cited at `x_build_group`) to the cell's entries or the overflow.
    /// Writes `src` at position `at`: with `X_WIDE_IN_OVERFLOW` the overflow from `at`'s other
    /// bits, else the entries of the cell in slot `at >> 2` from entry `at & 3`.
    ///
    /// # Safety
    /// The positions `at..at + src.len()` lie in one cell's entries or in the overflow, and no other
    /// thread reads or writes them while this runs.
    #[inline(always)]
    unsafe fn put(&self, at: u32, src: &[i32]) {
        let dst: *mut i32 = if at & X_WIDE_IN_OVERFLOW == 0 {
            let (slot, entry) = ((at >> 2) as usize, (at & 3) as usize);
            #[cfg(debug_assertions)]
            debug_assert!(slot < self.lens.0 && entry + src.len() <= AA_HITS_PER_CELL);
            // SAFETY: (caller) the slot is a cell of the backbone; no reference to it is made.
            unsafe { std::ptr::addr_of_mut!((*self.backbone.add(slot)).entries) }
                .cast::<i32>()
                .wrapping_add(entry)
        } else {
            let at = (at & !X_WIDE_IN_OVERFLOW) as usize;
            #[cfg(debug_assertions)]
            debug_assert!(at + src.len() <= self.lens.1);
            self.overflow.wrapping_add(at)
        };
        for (k, &offset) in src.iter().enumerate() {
            // SAFETY: (caller) `dst..dst + src.len()` lies in the cell's entries or the overflow
            // and is written by this call alone.
            unsafe { dst.add(k).write(offset) };
        }
    }
}

/// No NCBI counterpart: the totals of a stripe of cells (step (b)): overflow entries, longest
/// count, cells with entries.
#[derive(Clone, Copy, Default)]
struct XWideStripe {
    overflow: usize,
    longest: u32,
    used: usize,
}

impl XPairLutQuery {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 and 283-358 (cited
    // at `x_build_group`): the exact backbone, the appends of s_AddWordHits for s ascending, then
    // BlastAaLookupFinalize's layout, for the locations of each table.
    /// The tables `x_build_group(&tables[t], mapped)` builds, for every t (module comment), built
    /// together with 32-bit indices on `pool` (on this thread when it has one thread).
    pub(crate) fn x_build_tables_wide(
        &self,
        tables: &[Vec<(i32, i32)>],
        pool: &crate::utils::threading::SearchPool<'_>,
        mapped: bool,
    ) -> Vec<XPairLutGroup> {
        self.x_build_tables_wide_with(tables, pool, mapped, X_WIDE_RANGE_WORK_PER_CELL)
    }

    // No NCBI counterpart (scheduling only): `x_build_tables_wide` with the least work of a source
    // range per backbone cell as an argument (the unit tests take 0: as many ranges as threads, up to
    // `X_WIDE_MOST_PARTS`, whatever a table's size).
    /// `x_build_tables_wide` with `range_work` in place of `X_WIDE_RANGE_WORK_PER_CELL`.
    fn x_build_tables_wide_with(
        &self,
        tables: &[Vec<(i32, i32)>],
        pool: &crate::utils::threading::SearchPool<'_>,
        mapped: bool,
        range_work: u64,
    ) -> Vec<XPairLutGroup> {
        let word_length = self.word_length;
        assert!(
            x_wide_covers(word_length),
            "joint lookup build: word length {word_length}"
        );
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        debug_assert_eq!(charsize, X_WIDE_CHARSIZE);
        let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);
        assert!(backbone_size < 1 << 29, "joint lookup build: backbone size");
        let query = &self.prepared.concat_query;
        let threads = x_wide_threads(pool);

        // 1. The sources.
        let sources = XWideSources::new(query, tables, word_length, backbone_size);

        // 2. The lists T(s), once for every table.
        let mut scores = XWideScores {
            rows: [[0; LOOKUP_ALPHABET_SIZE]; LOOKUP_ALPHABET_SIZE],
            row_max: [0; LOOKUP_ALPHABET_SIZE],
        };
        for (a, row) in scores.rows.iter_mut().enumerate() {
            for (b, score) in row.iter_mut().enumerate() {
                *score = lookup_matrix_score(self.matrix, a, b);
            }
        }
        scores.row_max.copy_from_slice(&self.row_max);
        let n = sources.cells.len();
        let parts = (threads * X_WIDE_LIST_PARTS_PER_THREAD).min(n).max(1);
        let bounds: Vec<usize> = (0..=parts).map(|r| n * r / parts).collect();
        let threshold = self.threshold;
        let made = x_wide_parts(
            pool,
            bounds.windows(2).map(|b| b[0]..b[1]).collect(),
            &|range: std::ops::Range<usize>| {
                let mut targets: Vec<u32> = Vec::new();
                let mut ends: Vec<usize> = Vec::with_capacity(range.len());
                for i in range {
                    // blast_aalookup.c:490: w = query + offset_list[2];
                    let first = usize::try_from(sources.offsets[sources.first[i] as usize])
                        .expect("NCBI BLAST query offset must fit");
                    x_wide_targets(
                        &query[first..],
                        word_length,
                        &scores,
                        threshold,
                        &mut targets,
                    );
                    ends.push(targets.len());
                }
                targets.shrink_to_fit();
                (targets, ends)
            },
        );
        let lists = XWideLists {
            bounds,
            parts: made,
        };

        // 3. Each table, its counts kept from one table to the next.
        let mut counts: Vec<u32> = Vec::new();
        sources
            .by_table
            .iter()
            .map(|table_sources| {
                self.x_wide_table(
                    table_sources,
                    &sources,
                    &lists,
                    &mut counts,
                    (pool, threads),
                    (backbone_size, mapped, range_work),
                )
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-299,338-358 (cited at
    // `x_build_group`): the counts, the cells' layout and the copy of the chains.
    /// One table of `x_build_tables_wide` from its sources `table_sources` (steps (a)-(d)).
    fn x_wide_table(
        &self,
        table_sources: &[XWideSource],
        sources: &XWideSources,
        lists: &XWideLists,
        counts: &mut Vec<u32>,
        (pool, threads): (&crate::utils::threading::SearchPool<'_>, usize),
        (backbone_size, mapped, range_work): (usize, bool, u64),
    ) -> XPairLutGroup {
        let word_length = self.word_length;
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = X_WIDE_CHARSIZE;
        // The ranges: range p ends at the first source whose work so far reaches total * (p + 1) /
        // P (the work of a source: its list, read in (a) and (d), and its appends in (d)). Each
        // range costs passes over a count per backbone cell ((a)-(c)), so a table takes no more
        // ranges than its work pays for (`X_WIDE_RANGE_WORK_PER_CELL`).
        let work =
            |s: &XWideSource| (lists.list(s.i).len() as u64 + 1) * (u64::from(s.to - s.from) + 1);
        let total: u64 = table_sources.iter().map(work).sum();
        let paid = total
            .checked_div(range_work * backbone_size as u64)
            .unwrap_or(u64::MAX);
        let parts = threads
            .min(X_WIDE_MOST_PARTS)
            .min(usize::try_from(paid).unwrap_or(usize::MAX))
            .min(table_sources.len())
            .max(1);
        let mut bounds: Vec<usize> = vec![table_sources.len(); parts + 1];
        bounds[0] = 0;
        let (mut cut, mut so_far) = (1usize, 0u64);
        for (j, s) in table_sources.iter().enumerate() {
            if cut == parts {
                break;
            }
            so_far += work(s);
            while cut < parts && so_far * parts as u64 >= total * cut as u64 {
                bounds[cut] = j + 1;
                cut += 1;
            }
        }
        let ranges = || bounds.windows(2).map(|b| &table_sources[b[0]..b[1]]);
        // Each range zeroes its own counts in (a).
        counts.truncate(parts * backbone_size);
        counts.resize(parts * backbone_size, 0);

        // (a) Each range's appends per cell, and the bounds of the entries of its sources with a
        // target.
        let bounds_of: Vec<(u32, u32)> = x_wide_parts(
            pool,
            counts.chunks_mut(backbone_size).zip(ranges()).collect(),
            &|(c, range): (&mut [u32], &[XWideSource])| {
                c.fill(0);
                let (mut q_min, mut q_max) = (u32::MAX, 0u32);
                for &s in range {
                    let list = lists.list(s.i);
                    let k = s.to - s.from;
                    for &w in list {
                        c[w as usize] += k;
                    }
                    if !list.is_empty() {
                        for &q in sources.offsets(s) {
                            q_min = q_min.min(q as u32);
                            q_max = q_max.max(q as u32);
                        }
                    }
                }
                (q_min, q_max)
            },
        );
        let q_min = bounds_of.iter().map(|b| b.0).min().unwrap_or(u32::MAX);
        let q_max = bounds_of.iter().map(|b| b.1).max().unwrap_or(0);

        // (b) The stripes' totals.
        let counted: &[u32] = counts;
        let stripe = backbone_size.div_ceil(parts).next_multiple_of(32);
        let stripes = backbone_size.div_ceil(stripe);
        let summed: Vec<XWideStripe> = x_wide_parts(pool, (0..stripes).collect(), &|j: usize| {
            let cells = j * stripe..((j + 1) * stripe).min(backbone_size);
            let mut totals = XWideStripe::default();
            for w in cells {
                let count: u32 = (0..parts).map(|p| counted[p * backbone_size + w]).sum();
                if count == 0 {
                    continue;
                }
                totals.used += 1;
                if count as usize > AA_HITS_PER_CELL {
                    totals.overflow += count as usize;
                }
                totals.longest = totals.longest.max(count);
            }
            totals
        });
        let overflow_size: usize = summed.iter().map(|s| s.overflow).sum();
        let used: usize = summed.iter().map(|s| s.used).sum();
        let longest_chain = i32::try_from(summed.iter().map(|s| s.longest).max().unwrap_or(0))
            .expect("NCBI BLAST lookup chain hit count must fit Int4");
        assert!(
            overflow_size < X_WIDE_IN_OVERFLOW as usize,
            "NCBI BLAST lookup overflow cursor must fit Int4"
        );

        // (c) The layout: BlastAaLookupFinalize's rule, compacted where the branch-free scan takes
        // the table (`x_scan_bf_takes`, as `x_build_group` decides).
        let compact = crate::algorithm::tblastx::x_scan_bf::x_scan_bf_takes(used, backbone_size);
        let backbone_len = if compact { used } else { backbone_size };
        let mut backbone: Vec<BackboneCell> =
            Vec::with_capacity(x_group_capacity::<BackboneCell>(backbone_len, mapped));
        let mut overflow: Vec<i32> =
            Vec::with_capacity(x_group_capacity::<i32>(overflow_size, mapped));
        let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
        let mut pv: Vec<u32> = vec![0u32; pv_size];
        {
            // Each stripe's columns of counts, its cells (kept ones only when compacted), its pv
            // words (a stripe is a multiple of 32 cells) and its part of the overflow.
            let mut columns: Vec<Vec<&mut [u32]>> =
                (0..stripes).map(|_| Vec::with_capacity(parts)).collect();
            for c in counts.chunks_mut(backbone_size) {
                for (column, part) in columns.iter_mut().zip(c.chunks_mut(stripe)) {
                    column.push(part);
                }
            }
            let mut cells_left = &mut backbone.spare_capacity_mut()[..backbone_len];
            let mut overflow_left = &mut overflow.spare_capacity_mut()[..overflow_size];
            let mut pv_left = &mut pv[..];
            let mut placed = Vec::with_capacity(stripes);
            let (mut slot_base, mut cursor_base) = (0usize, 0usize);
            for (j, (column, totals)) in columns.into_iter().zip(&summed).enumerate() {
                let first = j * stripe;
                let len = (stripe).min(backbone_size - first);
                let kept = if compact { totals.used } else { len };
                let (cells, rest) = std::mem::take(&mut cells_left).split_at_mut(kept);
                cells_left = rest;
                let (part, rest) = std::mem::take(&mut overflow_left).split_at_mut(totals.overflow);
                overflow_left = rest;
                let words = (stripe >> PV_ARRAY_BTS).min(pv_left.len());
                let (pv_part, rest) = std::mem::take(&mut pv_left).split_at_mut(words);
                pv_left = rest;
                placed.push((
                    (first, slot_base, cursor_base),
                    column,
                    cells,
                    part,
                    pv_part,
                ));
                slot_base += if compact { totals.used } else { len };
                cursor_base += totals.overflow;
            }
            debug_assert!(cells_left.is_empty() && overflow_left.is_empty());
            x_wide_parts(pool, placed, &|(
                (first, slot_base, mut cursor),
                mut column,
                cells,
                part,
                pv_part,
            ): (
                (usize, usize, usize),
                Vec<&mut [u32]>,
                &mut [std::mem::MaybeUninit<BackboneCell>],
                &mut [std::mem::MaybeUninit<i32>],
                &mut [u32],
            )| {
                let len = column[0].len();
                let mut kept = 0usize;
                for local in 0..len {
                    let count: u32 = column.iter().map(|c| c[local]).sum();
                    if count == 0 {
                        if !compact {
                            cells[local].write(BackboneCell::default());
                        }
                        continue;
                    }
                    // pv_set at the stripe's first cell + local (a multiple of 32).
                    pv_part[local >> PV_ARRAY_BTS] |=
                        1u32 << (local & crate::algorithm::tblastx::lookup::PV_ARRAY_MASK);
                    let (at_cell, slot) = if compact {
                        kept += 1;
                        (kept - 1, slot_base + kept - 1)
                    } else {
                        (local, first + local)
                    };
                    let mut cell = BackboneCell {
                        num_used: count as i32,
                        ..BackboneCell::default()
                    };
                    let mut at = if count as usize > AA_HITS_PER_CELL {
                        cell.entries[0] = cursor as i32;
                        cursor += count as usize;
                        (cursor - count as usize) as u32 | X_WIDE_IN_OVERFLOW
                    } else {
                        (slot as u32) << 2
                    };
                    cells[at_cell].write(cell);
                    for c in column.iter_mut() {
                        let n = c[local];
                        c[local] = at;
                        at += n;
                    }
                }
                debug_assert!(!compact || kept == cells.len());
                // The overflow's entries are all written in (d); zeroed here first, by the
                // stripes side by side.
                for entry in part.iter_mut() {
                    entry.write(0);
                }
            });
        }
        // SAFETY: (c) wrote every cell of the backbone (each stripe its `kept` cells: all of its
        // cells, or, compacted, its cells with entries, which `summed` counted) and every entry of
        // the overflow (each stripe its part, the parts' sizes summing to `overflow_size`).
        unsafe {
            backbone.set_len(backbone_len);
            overflow.set_len(overflow_size);
        }

        // (d) The appends.
        let out = XWideOut {
            backbone: backbone.as_mut_ptr(),
            overflow: overflow.as_mut_ptr(),
            #[cfg(debug_assertions)]
            lens: (backbone_len, overflow_size),
        };
        x_wide_parts(
            pool,
            counts.chunks_mut(backbone_size).zip(ranges()).collect(),
            &|(at, range): (&mut [u32], &[XWideSource])| {
                for &s in range {
                    let src = sources.offsets(s);
                    for &w in lists.list(s.i) {
                        let w = w as usize;
                        // SAFETY: by (c), range p's appends to cell w take the positions from
                        // at[w] (its first) for c_p[w] offsets, inside the cell's entries or its
                        // part of the overflow, after the ranges before p and before the ranges
                        // after p; the cells' parts of the overflow are disjoint. So no two tasks
                        // write a position, and the vectors are read or written by nothing else
                        // until the tasks end.
                        unsafe { out.put(at[w], src) };
                        at[w] += src.len() as u32;
                    }
                }
            },
        );

        let rank = if compact {
            crate::algorithm::tblastx::x_scan_bf::x_rank(&pv)
        } else {
            Vec::new()
        };
        XPairLutGroup {
            q_min: if longest_chain > 0 { q_min } else { 0 },
            q_max,
            rank,
            table: BlastAaLookupTable {
                backbone,
                overflow,
                pv,
                frame_bases: self.prepared.frame_bases.clone(),
                num_contexts: self.prepared.contexts.len(),
                query_length: self.query_length,
                word_length: word_length as i32,
                alphabet_size: alphabet_size as i32,
                charsize: charsize as i32,
                mask: compute_mask(word_length, charsize) as i32,
                longest_chain,
                threshold: self.threshold,
                row_max: self.row_max.clone(),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::tests::{frame, residues};
    use super::*;

    const FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];

    /// The table of the query offsets `bounds[g] <= q < bounds[g + 1]` filtered out of the full
    /// table, order kept, in `BlastAaLookupFinalize`'s layout (the rule of `x_pair_par.rs`
    /// `build_groups`), with the bounds of its entries.
    fn filtered(full: &BlastAaLookupTable, bounds: &[i32], g: usize) -> XPairLutGroup {
        let lo = bounds[g];
        let hi = bounds.get(g + 1).copied().unwrap_or(i32::MAX);
        let mut table = BlastAaLookupTable {
            backbone: vec![BackboneCell::default(); full.backbone.len()],
            overflow: Vec::new(),
            pv: vec![0u32; full.pv.len()],
            frame_bases: full.frame_bases.clone(),
            num_contexts: full.num_contexts,
            query_length: full.query_length,
            word_length: full.word_length,
            alphabet_size: full.alphabet_size,
            charsize: full.charsize,
            mask: full.mask,
            longest_chain: 0,
            threshold: full.threshold,
            row_max: full.row_max.clone(),
        };
        let (mut q_min, mut q_max) = (u32::MAX, 0u32);
        for idx in 0..full.backbone.len() {
            let chain: Vec<i32> = full
                .get_hits(idx)
                .iter()
                .copied()
                .filter(|&q| lo <= q && q < hi)
                .collect();
            if chain.is_empty() {
                continue;
            }
            pv_set(&mut table.pv, idx);
            table.backbone[idx].num_used = chain.len() as i32;
            if chain.len() <= AA_HITS_PER_CELL {
                table.backbone[idx].entries[..chain.len()].copy_from_slice(&chain);
            } else {
                table.backbone[idx].entries[0] = table.overflow.len() as i32;
                table.overflow.extend_from_slice(&chain);
            }
            table.longest_chain = table.longest_chain.max(chain.len() as i32);
            for &q in &chain {
                q_min = q_min.min(q as u32);
                q_max = q_max.max(q as u32);
            }
        }
        XPairLutGroup {
            q_min: if table.longest_chain > 0 { q_min } else { 0 },
            q_max,
            table,
            rank: Vec::new(),
        }
    }

    /// Asserts that `wide` equals `reference` (a table of the full backbone), compacting the
    /// reference where the wide table is compacted.
    fn assert_same(g: usize, wide: &XPairLutGroup, mut reference: XPairLutGroup) {
        if !wide.rank.is_empty() {
            let rank = crate::algorithm::tblastx::x_scan_bf::x_compact(&mut reference.table);
            assert_eq!(wide.rank, rank, "table {g}: rank differs");
        } else {
            assert!(
                !crate::algorithm::tblastx::x_scan_bf::x_scan_bf_choice(&reference.table),
                "table {g}: the wide table is not compacted where the scan takes it"
            );
        }
        x_assert_same_group(
            g,
            &wide.table,
            (wide.q_min, wide.q_max),
            &reference.table,
            (reference.q_min, reference.q_max),
        );
    }

    /// The counts of `check`: the full table's entries, the tables with entries, the entries in
    /// the tables' overflow arrays and the compacted tables.
    #[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
    struct Seen {
        entries: usize,
        with_entries: usize,
        in_overflow: usize,
        compacted: usize,
    }

    /// Builds the wide tables of `queries` (`matrix`, `word`-letter words, `threshold`) for the
    /// group starts `starts` (None: one table per context) on `pool`, and compares each, field by
    /// field, with the filter of the chain builder's full table (`build_lookup_from_prepared`) and,
    /// for words of 2 and 3 letters, with `x_build_group`'s table; with one group, the table is the
    /// full table itself. The stub (`x_into_stub`) carries the full table's header and
    /// `longest_chain`.
    fn check(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        matrix: ScoringMatrix,
        word: usize,
        starts: Option<&[usize]>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Seen {
        let params = crate::stats::lookup_protein_params_ungapped(matrix);
        let (full, full_contexts) =
            crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
                queries, threshold, &params, true, matrix, word,
            );
        assert_eq!(full.word_length, word as i32);
        let q = XPairLutQuery::x_prepare_for(queries, threshold, true, matrix, word)
            .expect("a builder covers the word length");
        x_assert_same_contexts(q.contexts(), &full_contexts);
        let starts: Vec<usize> = match starts {
            Some(starts) => starts.to_vec(),
            None => (0..q.contexts().len()).collect(),
        };
        let bounds: Vec<i32> = starts.iter().map(|&c| full.frame_bases[c]).collect();
        let locations = q
            .x_group_locations(&starts)
            .expect("every location inside one group");
        let wide = q.x_build_tables_wide(&locations, pool, true);
        assert_eq!(wide.len(), starts.len());
        // As many source ranges as the pool has threads (up to `X_WIDE_MOST_PARTS`), whatever the
        // tables' size: the same tables.
        let ranged = q.x_build_tables_wide_with(&locations, pool, true, 0);
        for (g, (a, b)) in wide.iter().zip(&ranged).enumerate() {
            assert_eq!(a.rank, b.rank, "table {g}: rank differs with every range");
            x_assert_same_group(
                g,
                &b.table,
                (b.q_min, b.q_max),
                &a.table,
                (a.q_min, a.q_max),
            );
        }
        let mut seen = Seen {
            entries: full.backbone.iter().map(|c| c.num_used as usize).sum(),
            ..Seen::default()
        };
        for (g, table) in wide.iter().enumerate() {
            assert_same(g, table, filtered(&full, &bounds, g));
            if x_direct_covers(word) {
                let direct = q.x_build_group(&locations[g], true);
                assert_eq!(
                    direct.rank, table.rank,
                    "table {g}: rank differs from x_build_group"
                );
                x_assert_same_group(
                    g,
                    &table.table,
                    (table.q_min, table.q_max),
                    &direct.table,
                    (direct.q_min, direct.q_max),
                );
            }
            seen.with_entries += usize::from(table.table.longest_chain > 0);
            seen.in_overflow += table.table.overflow.len();
            seen.compacted += usize::from(!table.rank.is_empty());
        }
        if starts == [0] {
            // One group: the full table itself.
            let entries: Vec<i32> = (0..full.backbone.len())
                .flat_map(|idx| full.get_hits(idx).iter().copied())
                .collect();
            let reference = XPairLutGroup {
                q_min: entries.iter().map(|&q| q as u32).min().unwrap_or(0),
                q_max: entries.iter().map(|&q| q as u32).max().unwrap_or(0),
                rank: Vec::new(),
                table: BlastAaLookupTable {
                    backbone: full.backbone.clone(),
                    overflow: full.overflow.clone(),
                    pv: full.pv.clone(),
                    frame_bases: full.frame_bases.clone(),
                    row_max: full.row_max.clone(),
                    ..full
                },
            };
            assert_same(0, &wide[0], reference);
        }
        let (stub, contexts) = q.x_into_stub(wide.iter().map(|g| &g.table));
        x_assert_stub_scalars(&stub, &full);
        x_assert_same_contexts(&contexts, &full_contexts);
        seen
    }

    /// Frames of `len` random letters (1 to 27: every letter but the gap) with a SEG-like mask
    /// now and then.
    fn frames(seed: &mut u64, n: usize, len: usize) -> Vec<QueryFrame> {
        (0..n)
            .map(|c| {
                let len = len + 13 * c;
                let masks = match c % 3 {
                    0 => vec![(5, 12)],
                    1 => vec![],
                    _ => vec![(0, 3), (len / 2, len / 2 + 9)],
                };
                frame(residues(seed, len), masks, FRAMES[c])
            })
            .collect()
    }

    /// `f` on a one-thread pool and on pools of 4 and 8 threads.
    fn on_pools(f: impl Fn(&crate::utils::threading::SearchPool<'_>) + Sync) {
        f(&crate::utils::threading::SearchPool::x_serial());
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        for threads in [4usize, 8] {
            crate::utils::threading::with_search_pool(threads, "tblastx", |pool| {
                assert_eq!(x_wide_threads(pool), threads);
                f(pool);
                Ok(())
            })
            .expect("search pool");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467,504-509,580-590
    // (cited in the module comment): the chain of every cell, then BlastAaLookupFinalize's layout.
    // Words of 2, 3 and 4 letters, four matrices, thresholds 0 (no neighbours), 11 and 13, one
    // query of six contexts with SEG masks, a context shorter than a word, an empty one, a context
    // of a repeated letter (cells far above AA_HITS_PER_CELL), and three queries of six contexts in
    // groups of one, three and all contexts, on pools of 1, 4 and 8 threads.
    #[test]
    fn wide_tables_equal_the_filtered_full_table() {
        let mut seed = 0x3b17_5eed_0000_0004u64;
        let mut repeated = vec![1u8; 120];
        repeated.extend(residues(&mut seed, 40));
        let one_query = vec![vec![
            frame(residues(&mut seed, 260), vec![(10, 30), (200, 201)], 1),
            frame(residues(&mut seed, 3), vec![], 2),
            frame(repeated, vec![(70, 99)], 3),
            frame(residues(&mut seed, 230), vec![(0, 5)], -1),
            frame(Vec::new(), vec![], -2),
            frame(residues(&mut seed, 240), vec![], -3),
        ]];
        let three: Vec<Vec<QueryFrame>> = (0..3).map(|_| frames(&mut seed, 6, 60)).collect();
        let matrices = [
            ScoringMatrix::Blosum62,
            ScoringMatrix::Pam30,
            ScoringMatrix::Blosum45,
            ScoringMatrix::Pam250,
        ];
        let totals = std::sync::Mutex::new(Seen::default());
        on_pools(|pool| {
            for word in [2usize, 3, 4] {
                for (m, &matrix) in matrices.iter().enumerate() {
                    for threshold in [0, 11, 13] {
                        // The larger thresholds at 4 letters only with the matrices of few
                        // neighbours there (BLOSUM45 and PAM250 at 11 make millions of entries).
                        if word == 4 && threshold == 11 && m >= 2 {
                            continue;
                        }
                        let seen = check(&one_query, threshold, matrix, word, None, pool);
                        // The 3-letter context has words of 2 and 3 letters, none of 4.
                        let filled = if word <= 3 { 5 } else { 4 };
                        assert_eq!(
                            seen.with_entries, filled,
                            "word {word} threshold {threshold}"
                        );
                        let mut t = totals.lock().expect("totals");
                        t.entries += seen.entries;
                        t.in_overflow += seen.in_overflow;
                        t.compacted += seen.compacted;
                        drop(t);
                        for starts in [None, Some(&[0, 3, 6, 9, 12, 15][..]), Some(&[0, 1, 7][..])]
                        {
                            check(&three, threshold, matrix, word, starts, pool);
                        }
                        let seen = check(&three, threshold, matrix, word, Some(&[0]), pool);
                        assert_eq!(seen.with_entries, 1);
                    }
                }
            }
        });
        let t = totals.into_inner().expect("totals");
        assert!(
            t.entries > 0 && t.in_overflow > 0 && t.compacted > 0,
            "{t:?}"
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590 (cited in the
    // module comment). Threshold 1 (most words are neighbours) and INT_MIN (every word, NCBI's
    // conversion of a threshold beyond Int4): tables of the full backbone (more than half of the
    // cells used), cells with long chains, at each word length.
    #[test]
    fn wide_tables_with_most_or_every_word_a_neighbour() {
        let mut seed = 0x3b17_5eed_0000_0005u64;
        let query = vec![vec![
            frame(residues(&mut seed, 9), vec![], 1),
            frame(residues(&mut seed, 8), vec![(2, 4)], 2),
            frame(residues(&mut seed, 3), vec![], 3),
            frame(residues(&mut seed, 7), vec![], -1),
        ]];
        on_pools(|pool| {
            for word in [2usize, 3, 4] {
                for matrix in [ScoringMatrix::Blosum62, ScoringMatrix::Pam30] {
                    let seen = check(&query, 1, matrix, word, None, pool);
                    assert!(seen.with_entries >= 3);
                    let seen = check(&query, i32::MIN, matrix, word, None, pool);
                    // Every valid word (28 letters) of every query word.
                    let cells = 28usize.pow(word as u32);
                    let words: usize = [9usize, 8, 3, 7]
                        .iter()
                        .enumerate()
                        .map(|(c, &len)| {
                            if c == 1 {
                                // The mask (2, 4) (letters 2 and 3) leaves [0, 1] and [4, 7].
                                [2usize, 4]
                                    .iter()
                                    .map(|&l| (l + 1).saturating_sub(word))
                                    .sum()
                            } else {
                                (len + 1).saturating_sub(word)
                            }
                        })
                        .sum();
                    assert_eq!(seen.entries, words * cells, "word {word}");
                    // Every table with entries holds every valid word: of the full backbone; the
                    // context of 3 letters has no 4-letter word (an empty table, compacted).
                    let filled = if word <= 3 { 4 } else { 3 };
                    assert_eq!(
                        (seen.with_entries, seen.compacted),
                        (filled, 4 - filled),
                        "word {word}: every table with entries of the full backbone"
                    );
                    check(&query, i32::MIN, matrix, word, Some(&[0]), pool);
                }
            }
        });
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1047 (cited at
    // `x_group_locations`). Contexts masked from end to end, groups without a location, and a
    // query without any word: the tables without entries are empty in both builders.
    #[test]
    fn wide_tables_of_masked_and_empty_groups() {
        let mut seed = 0x3b17_5eed_0000_0006u64;
        let masked: Vec<Vec<QueryFrame>> = vec![(0..6)
            .map(|c| {
                let masks = if c % 2 == 0 {
                    vec![(0, 50)]
                } else {
                    vec![(3, 9)]
                };
                frame(residues(&mut seed, 50), masks, FRAMES[c])
            })
            .collect()];
        let none: Vec<Vec<QueryFrame>> = vec![(0..6)
            .map(|c| frame(residues(&mut seed, 30), vec![(0, 30)], FRAMES[c]))
            .collect()];
        on_pools(|pool| {
            for word in [2usize, 3, 4] {
                let seen = check(&masked, 13, ScoringMatrix::Blosum62, word, None, pool);
                assert_eq!(seen.with_entries, 3);
                let seen = check(&none, 13, ScoringMatrix::Blosum62, word, None, pool);
                assert_eq!(
                    (seen.entries, seen.with_entries, seen.in_overflow),
                    (0, 0, 0)
                );
            }
        });
    }

    // No NCBI counterpart: the word lengths each builder takes.
    #[test]
    fn the_builders_cover_words_of_two_to_four_letters() {
        assert!(!x_wide_covers(0));
        assert!((1..=4).all(x_wide_covers));
        assert!(!x_wide_covers(5));
        assert!((2..=4).all(x_tables_cover));
        assert!(x_direct_covers(3) && !x_direct_covers(4));
    }
}
