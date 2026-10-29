//! EXPERIMENT (PAIRLUT, graduated on 2026-10-22: no variable): the lookup tables of the TBLASTX
//! unit path (`tblastx/x_pair_par.rs`), one per query-context group, built directly
//! from the query instead of being filtered out of the batch's full table.
//!
//! This file is a child module of `x_lut_direct.rs` (declared there with `#[path]`, as that file
//! is for `backbone.rs`), so it uses the neighbour lists of the in-place builder (`x_target_list`,
//! the per-thread memo) and the reference builder pieces of `backbone.rs` (`prepare_lookup_query`,
//! `blast_lookup_index_query_exact_matches`) without changing them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:92-94,105-115
//! ```c
//!     for (loc = locations; loc; loc = loc->next) {
//!         Int4 from = loc->ssr->left;
//!         Int4 to = loc->ssr->right;
//! ...
//!         seq = query->sequence + from;
//!         word_target = seq + lut_word_length;
//!
//!         for (offset = from; offset <= to; offset++, seq++) {
//!
//!             if (seq >= word_target) {
//!                 BlastLookupAddWordHit(backbone,
//!                                       lut_word_length, charsize,
//!                                       seq - lut_word_length,
//!                                       offset - lut_word_length);
//!             }
//! ```
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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:338-358
//! ```c
//!       for (i = 0; i < lookup->backbone_size; i++) {
//!         /* if there are hits there, */
//!         if (lookup->thin_backbone[i] ) {
//!             Int4 * dest = NULL;
//!             /* set the corresponding bit in the pv_array */
//!             PV_SET(pv, i, PV_ARRAY_BTS);
//!             bbc[i].num_used = lookup->thin_backbone[i][1];
//!             /* if there are three or fewer hits, */
//!             if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
//!                 /* copy them into the thick_backbone cell */
//!                 dest = bbc[i].payload.entries;
//!             else /* more than three hits; copy to overflow array */
//!             {
//!                 bbc[i].payload.overflow_cursor = overflow_cursor;
//!                 dest = (Int4 *) lookup->overflow;
//!                 dest += overflow_cursor;
//!                 overflow_cursor += lookup->thin_backbone[i][1];
//!             }
//!             for (j=0; j <lookup->thin_backbone[i][1]; j++)
//!                 dest[j] = lookup->thin_backbone[i][j + 2];
//! ```
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
//!    of E(s) (blast_aalookup.c:490 `w = query + offset_list[2]`), and a memoised list is that
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
//! Corner cases: threshold 0 (T(s) = [s], no neighbours) and a threshold of INT_MIN (every word
//! is a neighbour) are `x_target_list`'s cases, unchanged; a group without a location (contexts
//! shorter than a word, or fully masked) has an empty exact backbone and an all-zero table, as the
//! filter of the full table has; with more than `X_MAX_GROUPS` contexts a group holds several
//! contexts and item 3 is checked per group, not per context.
//!
//! The unit tests of `x_pair_par.rs` build the full table and `build_groups` as well and stop at
//! the first difference (`x_assert_same_group`, `x_assert_stub_scalars`,
//! `x_assert_same_contexts`).

use super::*;

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
    pub table: BlastAaLookupTable,
    pub q_min: u32,
    pub q_max: u32,
}

impl XPairLutQuery {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
    // ```c
    // BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
    //                       &lookup_wrap->lut);
    // ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
    // BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
    //                           query, lookup_segments, 0);
    // ```
    // The query preparation of `build_ncbi_lookup` -> `build_ncbi_lookup_for_profile` (BLOSUM62,
    // LOOKUP_WORD_LENGTH, the same bounds and ideal parameters), by the reference function. The
    // table itself is not built here.
    /// The prepared query `build_ncbi_lookup(queries, threshold, _, check_ideal)` builds its
    /// table from. `None` when the direct builder does not apply (the preconditions of
    /// `x_build_direct`: diagnostics mode, or a backbone index wider than 16 bits).
    pub(crate) fn x_prepare(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        check_ideal: bool,
    ) -> Option<Self> {
        Self::x_prepare_for(queries, threshold, check_ideal, ScoringMatrix::Blosum62)
    }

    /// `x_prepare` with the scores of `matrix` (another matrix's search, `-matrix`): the query
    /// preparation of `build_ncbi_lookup_for_profile(queries, threshold, _, check_ideal, matrix,
    /// LOOKUP_WORD_LENGTH)`.
    #[inline(always)]
    pub(crate) fn x_prepare_for(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        check_ideal: bool,
        matrix: ScoringMatrix,
    ) -> Option<Self> {
        let word_length = LOOKUP_WORD_LENGTH;
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        if diagnostics_enabled() || word_length * charsize > 16 || word_length > 32 {
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
        // Same expression as build_lookup_from_prepared (BLAST_SequenceBlk length).
        let query_length = prepared
            .contexts
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
        Some(Self {
            prepared,
            threshold,
            matrix,
            word_length,
            query_length,
            row_max,
        })
    }

    /// The query contexts (those `build_ncbi_lookup` returns with its table).
    pub(crate) fn contexts(&self) -> &[QueryContext] {
        &self.prepared.contexts
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1056
    // ```c
    //       start_offset = query_info->contexts[context].query_offset;
    //       end_offset = query_info->contexts[context].query_length
    //           + start_offset - 1;
    // ...
    //       if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
    // ...
    //          tail = BlastSeqLocNew(tail ? &tail : complement_mask,
    //                                start_offset, end_offset);
    // ```
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
    // The count, layout and fill passes of `x_lut_direct.rs::x_build_direct` on the reference
    // exact backbone of `locations` (one group's lookup locations), reading the prepared query in
    // place. The neighbour lists come from the calling thread's memo, as there.
    /// The lookup table of one group: the table `build_groups` filters out of the full table for
    /// the group whose lookup locations are `locations` (module comment).
    pub(crate) fn x_build_group(&self, locations: &[(i32, i32)]) -> XPairLutGroup {
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

        X_MEMO.with(|memo| {
            let mut memo = memo.borrow_mut();
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
            let mut counts: Vec<u32> = vec![0; backbone_size];
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

            // Layout: BlastAaLookupFinalize's rule applied to the counts.
            let mut backbone: Vec<BackboneCell> = vec![BackboneCell::default(); backbone_size];
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
            let mut overflow: Vec<i32> = vec![0; overflow_size];
            let mut base: Vec<u32> = vec![0; backbone_size];
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

            // Fill pass: the reference's appends of this group, in the reference's order.
            let mut fill: Vec<u32> = vec![0; backbone_size];
            for idx in 0..backbone_size {
                let Some(offset_list) = exact_backbone[idx].as_deref() else {
                    continue;
                };
                // query_bias is 0 here (build_lookup_from_prepared passes 0).
                let offsets = lookup_chain_entries(offset_list);
                let list = lists[idx].as_deref().expect("listed in the count pass");
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

            XPairLutGroup {
                q_min: if longest_chain > 0 { q_min } else { 0 },
                q_max,
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
    // ```c
    //    case eAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:294-299
    // ```c
    //             if (lookup->thin_backbone[i][1] > longest_chain)
    //                 longest_chain = lookup->thin_backbone[i][1];
    //         }
    //     }
    //     lookup->overflow_size = overflow_cells_needed;
    //     lookup->longest_chain = longest_chain;
    // ```
    // The table the reference code holds when the unit path scans with the group tables only: the
    // header fields of the full table (expressions of `build_lookup_from_prepared`; the full
    // chain of a cell is the union of its group chains, so `longest_chain` is the largest summed
    // count) over an all-zero backbone and presence vector of the normal size and no overflow
    // array. A scan of it finds nothing and reads nothing out of bounds; `x_pair_par.rs` and the
    // per-subject dispatch make sure it is never scanned.
    /// The stub of the full table (no entries) and the contexts `build_ncbi_lookup` returns.
    /// `groups` are the tables of every group of this query.
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
            assert_eq!(table.backbone.len(), backbone_size, "group backbone size");
            for (total, cell) in totals.iter_mut().zip(&table.backbone) {
                *total += i64::from(cell.num_used);
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
                backbone: vec![BackboneCell::default(); backbone_size],
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:343-346,351-354
// ```c
//             PV_SET(pv, i, PV_ARRAY_BTS);
//             bbc[i].num_used = lookup->thin_backbone[i][1];
//             /* if there are three or fewer hits, */
//             if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
// ...
//                 bbc[i].payload.overflow_cursor = overflow_cursor;
//                 dest = (Int4 *) lookup->overflow;
//                 dest += overflow_cursor;
//                 overflow_cursor += lookup->thin_backbone[i][1];
// ```
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
// ```c
//     lookup->charsize = ilog2(BLASTAA_SIZE) + 1;
//     lookup->word_length = opt->word_size;
//
//     for (i = 0; i < lookup->word_length; i++)
//         lookup->backbone_size |= (BLASTAA_SIZE - 1) << (i * lookup->charsize);
//     lookup->backbone_size++;
//
//     lookup->mask = (1 << (opt->word_size * lookup->charsize)) - 1;
//     lookup->alphabet_size = BLASTAA_SIZE;
//     lookup->threshold = (Int4)opt->threshold;
// ```
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
// ```c
//    case eAaLookupTable:
//       offset_array_size = OFFSET_ARRAY_SIZE +
//          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
// ```
// No computation: the stub must carry the values the reference code reads from the full table
// (`longest_chain` for the offset-pair array, the header for the context search) and the full
// table's array sizes with no entry.
/// Assert that the stub has the full table's header fields and `longest_chain`, the full table's
/// backbone and pv sizes, and no entry.
#[cfg(test)]
pub(crate) fn x_assert_stub_scalars(stub: &BlastAaLookupTable, full: &BlastAaLookupTable) {
    x_assert_same_header("stub", 0, stub, full);
    assert_eq!(
        stub.longest_chain, full.longest_chain,
        "stub: longest_chain differs from the full table"
    );
    assert!(
        stub.backbone.len() == full.backbone.len() && stub.pv.len() == full.pv.len(),
        "stub: backbone or pv size differs from the full table"
    );
    assert!(
        stub.overflow.is_empty()
            && stub.pv.iter().all(|&word| word == 0)
            && stub
                .backbone
                .iter()
                .all(|cell| cell.num_used == 0 && cell.entries == [0; AA_HITS_PER_CELL]),
        "stub: the stub table must have no entry"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2797
// ```c
// Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
// ...
// loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ...
// if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
//    Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
// ```
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
            f_idx,
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
                && *f_idx == b.f_idx
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
mod tests {
    use super::super::tests::{frame, residues};
    use super::*;

    /// One context per group: the first context of each group.
    fn one_group_per_context(q: &XPairLutQuery) -> Vec<usize> {
        (0..q.contexts().len()).collect()
    }

    /// The table of group g as the filter of the full table: every chain filtered to the query
    /// offsets of the group (`bounds[g] <= q < bounds[g + 1]`), order kept, in the layout of
    /// `BlastAaLookupFinalize` (the rule `x_pair_par.rs::build_groups` applies).
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
        let mut q_min = u32::MAX;
        let mut q_max = 0u32;
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
        }
    }

    /// Builds the direct groups of `queries` for the group starts `starts` (None: one group per
    /// context) and compares each with the filter of the reference full table, and the stub with
    /// the full table. Returns (entries of the full table, groups with entries, entries in the
    /// overflow arrays of the groups).
    fn check(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        starts: Option<&[usize]>,
    ) -> (usize, usize, usize) {
        let params = crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (full, full_contexts) = build_ncbi_lookup(queries, threshold, &params, true);
        let q = XPairLutQuery::x_prepare(queries, threshold, true).expect("direct builder applies");
        x_assert_same_contexts(q.contexts(), &full_contexts);
        let starts: Vec<usize> = match starts {
            Some(starts) => starts.to_vec(),
            None => one_group_per_context(&q),
        };
        let bounds: Vec<i32> = starts.iter().map(|&c| full.frame_bases[c]).collect();
        let locations = q
            .x_group_locations(&starts)
            .expect("every location inside one group");
        assert_eq!(locations.len(), starts.len());
        let groups: Vec<XPairLutGroup> = locations.iter().map(|l| q.x_build_group(l)).collect();
        let mut with_entries = 0usize;
        let mut in_overflow = 0usize;
        for (g, direct) in groups.iter().enumerate() {
            let reference = filtered(&full, &bounds, g);
            x_assert_same_group(
                g,
                &direct.table,
                (direct.q_min, direct.q_max),
                &reference.table,
                (reference.q_min, reference.q_max),
            );
            with_entries += usize::from(direct.table.longest_chain > 0);
            in_overflow += direct.table.overflow.len();
        }
        let (stub, contexts) = q.x_into_stub(groups.iter().map(|g| &g.table));
        x_assert_stub_scalars(&stub, &full);
        x_assert_same_contexts(&contexts, &full_contexts);
        let entries = full.backbone.iter().map(|c| c.num_used as usize).sum();
        (entries, with_entries, in_overflow)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:92-94,105-106
    // ```c
    //     for (loc = locations; loc; loc = loc->next) {
    //         Int4 from = loc->ssr->left;
    //         Int4 to = loc->ssr->right;
    // ...
    //         seq = query->sequence + from;
    //         word_target = seq + lut_word_length;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:504-509,518-519
    // ```c
    //     if (lookup->threshold == 0 || score < lookup->threshold) {
    //         for (i = 0; i < offset_list[1]; i++) {
    //             BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
    //                                   lookup->charsize, w,
    //                                   query_bias + offset_list[i + 2]);
    //         }
    // ...
    //     if (lookup->threshold == 0)
    //         return;
    // ```
    #[test]
    fn direct_groups_equal_the_filtered_full_table() {
        let mut seed = 0x5eed_0001u64;
        let mut low_complexity = vec![1u8; 200]; // poly-A: cells far above AA_HITS_PER_CELL
        low_complexity.extend(residues(&mut seed, 60));
        // One query: six contexts with masks, a frame shorter than a word, an empty frame.
        let one_query = vec![vec![
            frame(
                residues(&mut seed, 500),
                vec![(10, 30), (400, 401), (498, 500)],
                1,
            ),
            frame(residues(&mut seed, 2), vec![], 2),
            frame(low_complexity, vec![(150, 199)], 3),
            frame(residues(&mut seed, 450), vec![(0, 5)], -1),
            frame(Vec::new(), vec![], -2),
            frame(residues(&mut seed, 480), vec![], -3),
        ]];
        for threshold in [0, 11, 13] {
            let (entries, with_entries, in_overflow) = check(&one_query, threshold, None);
            assert!(entries > 0);
            assert_eq!(with_entries, 4, "threshold {threshold}");
            assert!(in_overflow > 0, "the fixture must reach the overflow array");
        }
        // Every number of contexts from 1 to 6.
        for n in 1..=6usize {
            let frames: Vec<QueryFrame> = (0..n)
                .map(|c| {
                    let len = 40 + 37 * c;
                    frame(
                        residues(&mut seed, len),
                        if c % 2 == 0 { vec![(5, 12)] } else { vec![] },
                        [1i8, 2, 3, -1, -2, -3][c],
                    )
                })
                .collect();
            for threshold in [0, 11, 13] {
                let (_, with_entries, _) = check(std::slice::from_ref(&frames), threshold, None);
                assert_eq!(with_entries, n);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590
    // ```c
    //         for (i = 0; i < alphabet_size; i++) {
    //             if (score + row[i] >= threshold) {
    //                 subject_word[current_pos] = i;
    //                 for (j = 0; j < offset_list[1]; j++) {
    //                     BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
    //                                           charsize, subject_word,
    //                                           query_bias + offset_list[j + 2]);
    //                 }
    // ```
    // A threshold of INT_MIN (NCBI's conversion of a threshold beyond Int4) makes every word a
    // neighbour of every query word.
    #[test]
    fn direct_groups_with_every_word_a_neighbour() {
        let mut seed = 0x5eed_0002u64;
        let query = vec![vec![
            frame(residues(&mut seed, 9), vec![], 1),
            frame(residues(&mut seed, 7), vec![(2, 4)], 2),
            frame(residues(&mut seed, 2), vec![], 3),
            frame(residues(&mut seed, 8), vec![], -1),
        ]];
        let (entries, with_entries, in_overflow) = check(&query, i32::MIN, None);
        assert_eq!(with_entries, 3);
        // 7 + 1 + 6 word starts (the mask (2, 4) of the 7-letter frame leaves the locations [0, 1],
        // too short for a word, and [4, 6], one word), each in every one of the 28^3 cells.
        assert_eq!(entries, 14 * 28 * 28 * 28);
        // The cells of the groups with 7 and 6 words are in their overflow arrays.
        assert_eq!(in_overflow, 13 * 28 * 28 * 28);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1045-1047
    // ```c
    //       start_offset = query_info->contexts[context].query_offset;
    //       end_offset = query_info->contexts[context].query_length
    //           + start_offset - 1;
    // ```
    // Several queries: more than six contexts, taken as groups of several contexts, and a context
    // that is masked from end to end.
    #[test]
    fn direct_groups_of_several_contexts_and_an_all_masked_context() {
        let mut seed = 0x5eed_0003u64;
        let queries: Vec<Vec<QueryFrame>> = (0..3)
            .map(|q| {
                (0..6)
                    .map(|c| {
                        let len = 30 + 11 * ((q * 6 + c) % 7);
                        let masks = if q == 1 && c == 2 {
                            vec![(0, len)] // all masked
                        } else if c == 4 {
                            vec![(3, 9), (20, 22)]
                        } else {
                            vec![]
                        };
                        frame(residues(&mut seed, len), masks, [1i8, 2, 3, -1, -2, -3][c])
                    })
                    .collect()
            })
            .collect();
        for threshold in [0, 11, 13] {
            // 18 groups of one context: the all-masked context gives an empty group.
            let (_, with_entries, _) = check(&queries, threshold, None);
            assert_eq!(with_entries, 17);
            // Six groups of three contexts, and uneven groups.
            let (_, with_entries, _) = check(&queries, threshold, Some(&[0, 3, 6, 9, 12, 15]));
            assert_eq!(with_entries, 6);
            check(&queries, threshold, Some(&[0, 1, 7, 8, 9, 17]));
            // One group: the direct table is the full table.
            let (_, with_entries, _) = check(&queries, threshold, Some(&[0]));
            assert_eq!(with_entries, 1);
        }
        // A query whose every context is masked: no group has an entry.
        let masked: Vec<Vec<QueryFrame>> = vec![(0..6)
            .map(|c| {
                frame(
                    residues(&mut seed, 40),
                    vec![(0, 40)],
                    [1i8, 2, 3, -1, -2, -3][c],
                )
            })
            .collect()];
        assert_eq!(check(&masked, 13, None), (0, 0, 0));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1095-1099
    // ```c
    //          right = filter_start - 1;
    //
    //          /* Cache the tail of the list to avoid the overhead of traversing the
    //           * list when appending to it */
    //          tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
    // ```
    // A mask interval that starts past the end of its context (SEG never makes one) gives a
    // location that runs into the next context: such a batch is left to the reference builder.
    #[test]
    fn a_location_outside_one_group_declines() {
        let mut seed = 0x5eed_0004u64;
        let queries = vec![vec![
            frame(residues(&mut seed, 20), vec![(25, 28)], 1),
            frame(residues(&mut seed, 20), vec![], 2),
        ]];
        let q = XPairLutQuery::x_prepare(&queries, 13, true).expect("direct builder applies");
        assert!(q.x_group_locations(&[0, 1]).is_none());
        // Inside one group of both contexts the location is fine, and the tables agree.
        assert!(q.x_group_locations(&[0]).is_some());
        check(&queries, 13, Some(&[0]));
    }
}
