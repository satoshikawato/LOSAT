//! EXPERIMENT (PAIRLUT, graduated on 2026-10-22: no variable): the lookup tables of the TBLASTX
//! unit path (`tblastx/x_pair_par.rs`), one per query-context group, built directly
//! from the query instead of being filtered out of the batch's full table.
//!
//! This file is a child module of `x_lut_direct.rs` (declared there with `#[path]`, as that file
//! is for `backbone.rs`), so it uses the neighbour lists of the in-place builder (`x_target_list`,
//! the memos of TBLASTX's builds, `X_TBX_MEMOS`) and the reference builder pieces of `backbone.rs`
//! (`prepare_lookup_query`, `blast_lookup_index_query_exact_matches`) without changing them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:92-94,105-115
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:338-358
//!
//! # What is built
//!
//! `x_pair_par.rs::build_groups` makes the table of group g from the batch's full table: every
//! backbone chain filtered to the query offsets of the group (`group_of(q) == g`, the group being
//! a run of consecutive contexts), order kept, laid out by the rule of `BlastAaLookupFinalize`
//! above. `x_build_group` here makes the same table from the query: the reference exact backbone
//! (`blast_lookup_index_query_exact_matches`) of **the group's lookup locations only**, then the
//! count, layout and fill passes of `x_lut_direct.rs`. The full table is then not needed; the
//! reference code that still holds a `BlastAaLookupTable` gets `x_into_stub` (the full table's
//! header fields, no entries).
//!
//! # Why the direct table equals the filtered one
//!
//! Notation: E(s) is the exact chain of the source word s (cell s of the exact backbone, the
//! query offsets of the word in the order `BlastLookupIndexQueryExactMatches` adds them); T(s)
//! is the list of cells `s_AddWordHits` appends E(s) to, in its order (`x_target_list`).
//!
//! 1. The full chain of cell w is, for s ascending with w in T(s), E(s) in order
//!    (`x_lut_direct.rs` header; blast_aalookup.c:461-467 and 504-509, 580-590). Its filter to
//!    group g is therefore, for s ascending with w in T(s), E(s) filtered to g, in order.
//! 2. `BlastLookupIndexQueryExactMatches` (blast_lookup.c:92-131; `backbone.rs`
//!    `blast_lookup_index_query_exact_matches`) walks the locations in list order, and every
//!    location on its own (`seq`, `word_target` start again at `from`). So the exact backbone of
//!    a sub-list of the locations is the full one with each chain filtered to the offsets of
//!    those locations, order kept.
//! 3. A location [from, to] adds the word starts from..=to-word_length+1 at most. When `from`
//!    and `to` belong to group g (`x_group_locations` checks both ends with the `group_of` of
//!    `build_groups`; the groups are intervals of query offsets), all its word starts do. So
//!    "the offsets of the locations of group g" are exactly "the offsets q with group_of(q) = g",
//!    and by item 2 the exact backbone of the group's locations is E_g(s) = E(s) filtered to g.
//!    `prepare_lookup_query` makes the locations per context, inside
//!    [frame_base, frame_base + aa_len - 1] (blast_filter.c:1045-1056,1099,1111-1114: the
//!    complement of the context's masks), so the check holds for SEG masks; a location that fails
//!    it sends the batch to the reference builder (`x_group_locations` returns None).
//! 4. T(s) depends only on the residues of s, the matrix, the threshold, the word length and the
//!    alphabet (`x_lut_direct.rs` header). The backbone index of s fixes its residues (5 bits
//!    each; a word with a residue of 32 or more is not indexed, blast_lookup.c:117-120), so the
//!    list read at the first offset of E_g(s) is the list the reference reads at the first offset
//!    of E(s) (blast_aalookup.c:490, where `w` is `query` plus `offset_list[2]`), and a memoised list is that
//!    list again.
//! 5. The fill pass of group g appends, for s ascending with E_g(s) not empty and w in T(s),
//!    E_g(s) to cell w: by items 1-4 the entries of the filtered full chain in the same order (a
//!    source word with an empty E_g(s) adds nothing in either). The counts are the lengths of
//!    these chains, and `num_used`, the pv bits, the inline entries or the overflow cursor
//!    (assigned in ascending cell order to the cells with more than AA_HITS_PER_CELL entries),
//!    the overflow array and `longest_chain` are the same function of the counts and the chains
//!    in `build_groups` and in the layout pass here (`BlastAaLookupFinalize`'s rule); unused
//!    slots are zero in both. `q_min` / `q_max` are the minimum and maximum of the table's
//!    entries (0 / 0 for a table without entries), as in `build_groups`.
//! 6. The other fields (`frame_bases`, `num_contexts`, `query_length`, `word_length`,
//!    `alphabet_size`, `charsize`, `mask`, `threshold`, `row_max`) are computed by the expressions
//!    of `build_lookup_from_prepared` from the same prepared query; `build_groups` copies them
//!    from the full table. The stub's `longest_chain` is the maximum over the cells of the summed
//!    group counts: the groups partition every full chain, so this is the full table's value.
//!
//! Corner cases: threshold 0 (T(s) = \[s\], no neighbours) and a threshold of INT_MIN (every word
//! is a neighbour) are `x_target_list`'s cases, unchanged; a group without a location (contexts
//! shorter than a word, or fully masked) has an empty exact backbone and an all-zero table, as the
//! filter of the full table has; with more than `X_MAX_GROUPS` contexts a group holds several
//! contexts and item 3 is checked per group, not per context.
//!
//! The unit tests of `x_pair_par.rs` build the full table and `build_groups` as well and stop at
//! the first difference (`x_assert_same_group`, `x_assert_stub_scalars`,
//! `x_assert_same_contexts`).

use super::*;

// SC10b PR-3: the group tables of words whose backbone index does not fit 16 bits (`-word_size 4`),
// built for all groups at once with 32-bit indices on the search pool (`x_build_tables_wide`).
#[path = "x_pair_lut_wide.rs"]
pub(crate) mod x_pair_lut_wide;

/// The prepared query of one batch (the reference `prepare_lookup_query`) with the settings of
/// its lookup table.
pub(crate) struct XPairLutQuery {
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    /// `BLAST_SequenceBlk->length` of the concatenated query (expression of
    /// `build_lookup_from_prepared`).
    query_length: i32,
    /// Row maxima of the matrix (expression of `build_lookup_from_prepared`).
    row_max: Vec<i32>,
}

/// One direct group table with the bounds of its entries.
pub(crate) struct XPairLutGroup {
    /// With `rank` not empty, a compacted table (`x_scan_bf::x_compact`): its backbone holds only
    /// the cells with entries.
    pub table: BlastAaLookupTable,
    pub q_min: u32,
    pub q_max: u32,
    /// The presence-vector rank of a compacted table (`x_scan_bf::x_rank`), else empty.
    pub rank: Vec<u32>,
}

/// No NCBI counterpart: a block of a group's build from which glibc maps it (`x_group_capacity`).
const X_GROUP_MAPPED_FROM: usize = 64 << 10;

/// No NCBI counterpart: the capacity of a block of `elements` of a group's build. Batches are built
/// and searched on the pool's threads at several threads, and blocks of 100 KiB to 1.5 MiB that
/// stay in those threads' arenas between the kept tables hold their heaps' high water: 50 cDNAs x
/// 1 Mb and 40 x 40 Begomovirus genomes at 8 threads, the arenas 17-23 MiB, with one arena (the
/// same live data) a max RSS below that of the reference loop. So a block of at least
/// `X_GROUP_MAPPED_FROM` bytes takes the capacity that glibc maps whatever its size and whose free
/// leaves its mmap threshold as it is (`x_unthresholded_capacity`, as TBLASTN's tables,
/// x_lut_direct.rs); only what is written is touched.
fn x_group_capacity<T>(elements: usize, mapped: bool) -> usize {
    let size = std::mem::size_of::<T>();
    if mapped && elements.saturating_mul(size) >= X_GROUP_MAPPED_FROM {
        crate::utils::x_mapped::x_unthresholded_capacity(elements, size)
    } else {
        elements
    }
}

/// No NCBI counterpart: `vec![0; elements]` with `x_group_capacity`; `calloc` touches none of a
/// mapped block.
fn x_group_zeroed<T: Copy + Default>(elements: usize, mapped: bool) -> Vec<T> {
    let mut zeroed = vec![T::default(); x_group_capacity::<T>(elements, mapped)];
    zeroed.truncate(elements);
    zeroed
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:103-105
// (ComputeTableIndex: `wordsize` letters of `charsize` bits)
// The precondition of the direct builders (`x_build_direct`, `x_build_group`): their neighbour
// lists hold backbone indices in 16 bits (`x_target_list`), so the index of a word, `word_length`
// letters of `charsize` bits (5 for the 28 letters), must fit 16 bits: words of 2 and 3 letters
// (10 and 15 bits), not of 4 (20 bits). Words of 4 letters take the joint builder of 32-bit
// indices instead (`x_pair_lut_wide.rs`, `x_tables_cover`), so that the unit path's group tables
// are built from the query and not filtered out of a full table (`build_groups`).
/// Whether the direct builder takes a table of `word_length`-letter words: their backbone index
/// fits 16 bits.
#[inline(always)]
pub(crate) const fn x_direct_covers(word_length: usize) -> bool {
    let charsize = LOOKUP_ALPHABET_SIZE.ilog2() as usize + 1;
    word_length * charsize <= 16 && word_length <= 32
}

// No NCBI counterpart: which builder makes the unit path's group tables (the tables are the same
// either way, module comment).
/// Whether the unit path builds the group tables of `word_length`-letter words from the query:
/// the 16-bit direct builder (`x_direct_covers`: 2 and 3 letters) or the joint builder of 32-bit
/// indices (`x_pair_lut_wide::x_wide_covers`: 4 letters).
#[inline(always)]
pub(crate) const fn x_tables_cover(word_length: usize) -> bool {
    x_direct_covers(word_length) || x_pair_lut_wide::x_wide_covers(word_length)
}

impl XPairLutQuery {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
    // The query preparation of `build_ncbi_lookup` -> `build_ncbi_lookup_for_profile` (BLOSUM62,
    // LOOKUP_WORD_LENGTH, the same bounds and ideal parameters), by the reference function. The
    // table itself is not built here.
    /// The prepared query that `build_ncbi_lookup` (given `queries`, `threshold` and `check_ideal`) builds its
    /// table from. `None` when the direct builder does not apply (the preconditions of
    /// `x_build_direct`: diagnostics mode, or a backbone index wider than 16 bits).
    pub(crate) fn x_prepare(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        check_ideal: bool,
    ) -> Option<Self> {
        Self::x_prepare_for(
            queries,
            threshold,
            check_ideal,
            ScoringMatrix::Blosum62,
            LOOKUP_WORD_LENGTH,
        )
    }

    /// `x_prepare` with the scores of `matrix` (`-matrix`) and words of `word_length` letters
    /// (`-word_size`): the query preparation of `build_ncbi_lookup_for_profile(queries, threshold,
    /// _, check_ideal, matrix, word_length)`. `None` also for a word length whose backbone index
    /// neither builder covers (`x_tables_cover`).
    #[inline(always)]
    pub(crate) fn x_prepare_for(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        check_ideal: bool,
        matrix: ScoringMatrix,
        word_length: usize,
    ) -> Option<Self> {
        if diagnostics_enabled() || !x_tables_cover(word_length) {
            return None;
        }
        // Same expressions as build_ncbi_lookup_for_profile.
        let bounds = if matrix == ScoringMatrix::Blosum62 {
            (-4, 11)
        } else {
            matrix_score_bounds(matrix)
        };
        let prepared = prepare_lookup_query(
            queries,
            ideal_karlin_params_for_matrix(matrix, bounds),
            &compute_std_aa_composition(),
            word_length,
            check_ideal,
            matrix,
            bounds,
        );
        Some(Self::x_from_prepared(
            prepared,
            threshold,
            matrix,
            word_length,
        ))
    }

    /// The query of `build_lookup_from_prepared(prepared, threshold, matrix, word_length)`: the
    /// lookup query `prepared` with the scores of `matrix` and words of `word_length` letters, as
    /// `x_prepare_for` makes it for TBLASTX and `blastx::seed::Lookup` hands it over for BLASTX's
    /// table of 4-letter words (`x_build_full_wide`, item SC9 of line NP).
    pub(in super::super) fn x_from_prepared(
        prepared: PreparedLookupQuery,
        threshold: i32,
        matrix: ScoringMatrix,
        word_length: usize,
    ) -> Self {
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        // The expressions of `build_lookup_from_prepared`.
        let query_length = prepared
            .contexts
            .last()
            .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
            .unwrap_or(0);
        let row_max: Vec<i32> = (0..alphabet_size)
            .map(|i| {
                (0..alphabet_size)
                    .map(|j| lookup_matrix_score(matrix, i, j))
                    .max()
                    .unwrap_or(i32::MIN)
            })
            .collect();
        Self {
            prepared,
            threshold,
            matrix,
            word_length,
            query_length,
            row_max,
        }
    }

    /// The query contexts (those `build_ncbi_lookup` returns with its table).
    pub(crate) fn contexts(&self) -> &[QueryContext] {
        &self.prepared.contexts
    }

    /// The word length of the tables built from this query (`lookup->word_length`).
    pub(crate) fn word_length(&self) -> usize {
        self.word_length
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1056
    // The lookup locations are made per context, inside the context. Here they are split by
    // group (`starts[g]` = first context of group g, the groups of `x_pair_par.rs::group_starts`)
    // with the `group_of` of `build_groups`, and both ends of a location must fall into the same
    // group (module comment, item 3).
    /// The lookup locations of each group, in list order. `None` when a location does not lie
    /// inside one group (the reference builder must run).
    pub(crate) fn x_group_locations(&self, starts: &[usize]) -> Option<Vec<Vec<(i32, i32)>>> {
        let bounds: Vec<i32> = starts
            .iter()
            .map(|&c| self.prepared.frame_bases[c])
            .collect();
        let mut out: Vec<Vec<(i32, i32)>> = vec![Vec::new(); bounds.len()];
        if bounds.is_empty() {
            return self.prepared.lookup_locations.is_empty().then_some(out);
        }
        // The group of a query offset: the last group whose first context starts at or before it.
        let group_of = |q: i32| -> usize {
            let mut g = bounds.len() - 1;
            while g > 0 && q < bounds[g] {
                g -= 1;
            }
            g
        };
        for &(from, to) in &self.prepared.lookup_locations {
            let g = group_of(from);
            if from > to || group_of(to) != g {
                return None;
            }
            out[g].push((from, to));
        }
        Some(out)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-297,338-358
    // The count, layout and fill passes of `x_lut_direct.rs::x_build_direct` on the reference
    // exact backbone of `locations` (one group's lookup locations), reading the prepared query in
    // place. The neighbour lists come from a memo of TBLASTX's builds (`x_with_tbx_memo`).
    /// The lookup table of one group: the table `build_groups` filters out of the full table for
    /// the group whose lookup locations are `locations` (module comment).
    /// `mapped`: the build's large blocks take `x_group_capacity`'s mapped capacity.
    pub(crate) fn x_build_group(&self, locations: &[(i32, i32)], mapped: bool) -> XPairLutGroup {
        let threshold = self.threshold;
        let matrix = self.matrix;
        let word_length = self.word_length;
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        let mask = compute_mask(word_length, charsize);
        let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);
        let concat_query = &self.prepared.concat_query;
        let row_max = &self.row_max;
        // The reference exact backbone (blast_aalookup.c:448-456) of this group's locations.
        let mut exact_backbone = LookupBackboneChains::new(backbone_size);
        let mut skipped_invalid_residue = 0usize;
        blast_lookup_index_query_exact_matches(
            &mut exact_backbone,
            word_length as i32,
            charsize as i32,
            word_length as i32,
            concat_query,
            locations,
            &mut skipped_invalid_residue,
        );

        x_with_tbx_memo(|memo| {
            let key = (matrix, threshold, word_length, alphabet_size);
            if memo
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
            let lists = &mut memo.as_mut().expect("memo set above").lists;

            // Count pass: count[w] = number of appends the reference makes to cell w from the
            // source words of this group. The bounds of the entries are taken on the way: every
            // offset of a source word with a target is an entry of the table.
            let mut counts: Vec<u32> = x_group_zeroed(backbone_size, mapped);
            let mut q_min = u32::MAX;
            let mut q_max = 0u32;
            for idx in 0..backbone_size {
                let Some(offset_list) = exact_backbone[idx].as_deref() else {
                    continue;
                };
                let k = lookup_chain_num_used(offset_list) as u32;
                let list = lists[idx].get_or_insert_with(|| {
                    // blast_aalookup.c:490: w = query + offset_list[2];
                    let first =
                        usize::try_from(offset_list[2]).expect("NCBI BLAST query offset must fit");
                    x_target_list(
                        &concat_query[first..],
                        alphabet_size,
                        word_length,
                        charsize,
                        row_max,
                        threshold,
                        matrix,
                    )
                });
                for &w in list.iter() {
                    counts[w as usize] += k;
                }
                if !list.is_empty() {
                    for &q in lookup_chain_entries(offset_list) {
                        q_min = q_min.min(q as u32);
                        q_max = q_max.max(q as u32);
                    }
                }
            }

            // Layout: BlastAaLookupFinalize's rule applied to the counts. A table the branch-free
            // scan takes (`x_scan_bf_takes`, as `XGroup::new` decides for a table of the full
            // backbone) is laid out compacted (`x_scan_bf::x_compact`): its backbone holds the cells
            // with entries in index order, cell w at `slot[w]`; the entries, the overflow array and
            // the presence vector are those of the full layout.
            let used = counts.iter().filter(|&&count| count > 0).count();
            let compact =
                crate::algorithm::tblastx::x_scan_bf::x_scan_bf_takes(used, backbone_size);
            let backbone_len = if compact { used } else { backbone_size };
            let mut backbone: Vec<BackboneCell> =
                Vec::with_capacity(x_group_capacity::<BackboneCell>(backbone_len, mapped));
            backbone.resize(backbone_len, BackboneCell::default());
            let mut slot: Vec<u32> = if compact {
                x_group_zeroed(backbone_size, mapped)
            } else {
                Vec::new()
            };
            let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
            let mut pv: Vec<u32> = vec![0u32; pv_size];
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
            let mut overflow: Vec<i32> =
                Vec::with_capacity(x_group_capacity::<i32>(overflow_size, mapped));
            overflow.resize(overflow_size, 0);
            // Where the fill pass writes each cell's next offset: 0 in the cell, the start of its
            // run in the overflow array (a `usize`: a position of an array of 2^32 entries or
            // more, which a `u32` wrapped).
            let mut next: Vec<usize> = x_group_zeroed(backbone_size, mapped);
            let mut overflow_cursor = 0usize;
            let mut kept = 0usize;
            for (idx, &count) in counts.iter().enumerate() {
                if count == 0 {
                    continue;
                }
                pv_set(&mut pv, idx);
                let at = if compact {
                    slot[idx] = kept as u32;
                    kept += 1;
                    kept - 1
                } else {
                    idx
                };
                backbone[at].num_used = count as i32;
                if count as usize > AA_HITS_PER_CELL {
                    backbone[at].set_cursor(overflow_cursor);
                    next[idx] = overflow_cursor;
                    overflow_cursor += count as usize;
                }
            }

            // Fill pass: the reference's appends of this group, in the reference's order.
            for idx in 0..backbone_size {
                let Some(offset_list) = exact_backbone[idx].as_deref() else {
                    continue;
                };
                // query_bias is 0 here (build_lookup_from_prepared passes 0).
                let offsets = lookup_chain_entries(offset_list);
                let list = lists[idx].as_deref().expect("listed in the count pass");
                for &w in list {
                    let w = w as usize;
                    let at = next[w];
                    if counts[w] as usize <= AA_HITS_PER_CELL {
                        let cell = if compact { slot[w] as usize } else { w };
                        backbone[cell].entries[at..at + offsets.len()].copy_from_slice(offsets);
                    } else {
                        overflow[at..at + offsets.len()].copy_from_slice(offsets);
                    }
                    next[w] += offsets.len();
                }
            }

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
                    mask: mask as i32,
                    longest_chain,
                    threshold,
                    row_max: row_max.clone(),
                },
            }
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-266
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:294-299
    // The table the reference code holds when the unit path scans with the group tables only: the
    // header fields of the full table (expressions of `build_lookup_from_prepared`; the full
    // chain of a cell is the union of its group chains, so `longest_chain` is the largest summed
    // count) over an all-zero presence vector of the normal size, no backbone cells and no
    // overflow array. A scan of it tests presence bits only, finds nothing and reads nothing out
    // of bounds; `x_pair_par.rs` and the per-subject dispatch make sure it is never scanned.
    /// The stub of the full table (no entries) and the contexts `build_ncbi_lookup` returns.
    /// `groups` are the tables of every group of this query, of the full backbone or compacted
    /// (`x_scan_bf::x_compact`).
    pub(crate) fn x_into_stub<'a>(
        self,
        groups: impl Iterator<Item = &'a BlastAaLookupTable>,
    ) -> (BlastAaLookupTable, Vec<QueryContext>) {
        let word_length = self.word_length;
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        let mask = compute_mask(word_length, charsize);
        let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);
        let mut totals: Vec<i64> = vec![0; backbone_size];
        for table in groups {
            for (index, cell) in x_cells_with_entries(table, backbone_size) {
                totals[index] += i64::from(cell.num_used);
            }
        }
        let longest_chain = i32::try_from(totals.iter().copied().max().unwrap_or(0))
            .expect("NCBI BLAST lookup chain hit count must fit Int4");
        let PreparedLookupQuery {
            concat_query: _,
            lookup_locations: _,
            frame_bases,
            contexts,
            skipped_seg_mask: _,
        } = self.prepared;
        (
            BlastAaLookupTable {
                backbone: Vec::new(),
                overflow: Vec::new(),
                pv: vec![0u32; (backbone_size >> PV_ARRAY_BTS) + 1],
                frame_bases,
                num_contexts: contexts.len(),
                query_length: self.query_length,
                word_length: word_length as i32,
                alphabet_size: alphabet_size as i32,
                charsize: charsize as i32,
                mask: mask as i32,
                longest_chain,
                threshold: self.threshold,
                row_max: self.row_max,
            },
            contexts,
        )
    }
}

// No NCBI counterpart: reads a group table of the unit path in either layout (`x_into_stub`).
/// The cells with entries of `table`, with their indices, in index order: cell `i` of a table of
/// `backbone_size` cells, or, for a compacted table (`x_scan_bf::x_compact`), its `k`-th kept
/// cell for its `k`-th set presence bit.
pub(crate) fn x_cells_with_entries(
    table: &BlastAaLookupTable,
    backbone_size: usize,
) -> impl Iterator<Item = (usize, &BackboneCell)> + '_ {
    let full = table.backbone.len() == backbone_size;
    assert!(
        full || table.backbone.len() < backbone_size,
        "group backbone size"
    );
    table
        .pv
        .iter()
        .enumerate()
        .flat_map(|(word, &bits)| {
            (0..32usize)
                .filter(move |bit| (bits >> bit) & 1 != 0)
                .map(move |bit| (word << PV_ARRAY_BTS) | bit)
        })
        .enumerate()
        .map(move |(k, index)| (index, &table.backbone[if full { index } else { k }]))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:343-346,351-354
// No computation: the unit tests compare every field of the two tables of a group (the cells, the
// overflow array, the presence vector and the header) and the bounds of the entries.
/// Assert that the direct table of group `g` equals the one `build_groups` filtered out of the
/// full table, field by field.
#[cfg(test)]
pub(crate) fn x_assert_same_group(
    g: usize,
    direct: &BlastAaLookupTable,
    direct_bounds: (u32, u32),
    reference: &BlastAaLookupTable,
    reference_bounds: (u32, u32),
) {
    assert_eq!(
        direct.backbone.len(),
        reference.backbone.len(),
        "group {g}: backbone size differs"
    );
    for (idx, (a, b)) in direct.backbone.iter().zip(&reference.backbone).enumerate() {
        assert!(
            a.num_used == b.num_used && a.entries == b.entries,
            "group {g}: backbone cell {idx} differs (direct {} {:?} vs reference {} {:?})",
            a.num_used,
            a.entries,
            b.num_used,
            b.entries
        );
    }
    assert_eq!(
        direct.overflow.len(),
        reference.overflow.len(),
        "group {g}: overflow size differs"
    );
    if let Some(at) =
        (0..direct.overflow.len()).find(|&i| direct.overflow[i] != reference.overflow[i])
    {
        panic!(
            "group {g}: overflow entry {at} differs (direct {} vs reference {})",
            direct.overflow[at], reference.overflow[at]
        );
    }
    assert!(direct.pv == reference.pv, "group {g}: pv differs");
    assert_eq!(
        direct.longest_chain, reference.longest_chain,
        "group {g}: longest_chain differs"
    );
    x_assert_same_header("group", g, direct, reference);
    assert_eq!(
        direct_bounds, reference_bounds,
        "group {g}: q_min / q_max differ (direct vs reference)"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:236-245
// No computation: the header fields of two tables compared (all but `longest_chain`, which the
// callers compare with the value that applies).
#[cfg(test)]
fn x_assert_same_header(
    what: &str,
    g: usize,
    direct: &BlastAaLookupTable,
    reference: &BlastAaLookupTable,
) {
    assert!(
        direct.frame_bases == reference.frame_bases
            && direct.num_contexts == reference.num_contexts
            && direct.query_length == reference.query_length
            && direct.word_length == reference.word_length
            && direct.alphabet_size == reference.alphabet_size
            && direct.charsize == reference.charsize
            && direct.mask == reference.mask
            && direct.threshold == reference.threshold
            && direct.row_max == reference.row_max,
        "{what} {g}: table header differs"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-266
// No computation: the stub must carry the values the reference code reads from the full table
// (`longest_chain` for the offset-pair array, the header for the context search), the full
// table's presence-vector size and no entry.
/// Assert that the stub has the full table's header fields and `longest_chain`, the full table's
/// pv size with no bit set, no backbone cell and no overflow entry.
#[cfg(test)]
pub(crate) fn x_assert_stub_scalars(stub: &BlastAaLookupTable, full: &BlastAaLookupTable) {
    x_assert_same_header("stub", 0, stub, full);
    assert_eq!(
        stub.longest_chain, full.longest_chain,
        "stub: longest_chain differs from the full table"
    );
    assert_eq!(
        stub.pv.len(),
        full.pv.len(),
        "stub: pv size differs from the full table"
    );
    assert!(
        stub.overflow.is_empty()
            && stub.backbone.is_empty()
            && stub.pv.iter().all(|&word| word == 0),
        "stub: the stub table must have no entry"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
// No computation: the contexts handed out with the stub come from the same `prepare_lookup_query`
// as those `build_ncbi_lookup` returns; the unit tests check that they are the same, field by field
// (the Karlin parameters bitwise).
/// Assert that two context lists are equal.
#[cfg(test)]
pub(crate) fn x_assert_same_contexts(direct: &[QueryContext], reference: &[QueryContext]) {
    x_assert_same_contexts_as("", direct, reference);
}

// No NCBI counterpart: the check of `x_assert_same_contexts` for other unit tests.
/// `x_assert_same_contexts` with the messages starting with `what`.
#[cfg(test)]
pub(crate) fn x_assert_same_contexts_as(
    what: &str,
    direct: &[QueryContext],
    reference: &[QueryContext],
) {
    assert_eq!(
        direct.len(),
        reference.len(),
        "{what} number of contexts differs"
    );
    for (c, (a, b)) in direct.iter().zip(reference).enumerate() {
        let QueryContext {
            q_idx,
            frame,
            aa_seq,
            aa_seq_nomask,
            aa_len,
            orig_len,
            frame_base,
            is_valid,
            karlin_params,
        } = a;
        let bits = |p: &KarlinParams| {
            [
                p.lambda.to_bits(),
                p.k.to_bits(),
                p.h.to_bits(),
                p.alpha.to_bits(),
                p.beta.to_bits(),
            ]
        };
        assert!(
            *q_idx == b.q_idx
                && *frame == b.frame
                && *aa_seq == b.aa_seq
                && *aa_seq_nomask == b.aa_seq_nomask
                && *aa_len == b.aa_len
                && *orig_len == b.orig_len
                && *frame_base == b.frame_base
                && *is_valid == b.is_valid
                && bits(karlin_params) == bits(&b.karlin_params),
            "{what} context {c} differs"
        );
    }
}

#[cfg(test)]
#[path = "x_pair_lut_tests.rs"]
mod tests;
