//! EXPERIMENT (LOSAT_X_BXLUT / LOSAT_X_BXLUTSHADOW): the BLASTX protein lookup table
//! filled in place, with the neighbour-word sets memoised per thread.
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

/// EXPERIMENT (LOSAT_X_BXLUT): 0 = off (reference builder), 1 = in-place builder,
/// 2 = shadow (`LOSAT_X_BXLUTSHADOW`: both builders, compared).
pub(super) fn x_bxlut_mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BXLUTSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_BXLUT").is_some() {
            1
        } else {
            0
        }
    })
}

/// Neighbour lists of one (matrix, threshold, word length, alphabet) setting, by source cell.
struct XNeighborMemo {
    key: (ScoringMatrix, i32, usize, usize),
    lists: Vec<Option<Box<[u16]>>>,
}

thread_local! {
    static X_MEMO: RefCell<Option<XNeighborMemo>> = const { RefCell::new(None) };
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
/// EXPERIMENT (LOSAT_X_BXLUT): the table `build_lookup_from_prepared` returns, built without
/// chains (see the module comment). `None` when the reference builder must run (diagnostics
/// mode, or a backbone index wider than 16 bits).
pub(super) fn x_build_direct(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
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
    let mut exact_backbone = LookupBackboneChains::new(backbone_size);
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

        // Count pass: count[w] = number of appends the reference makes to cell w.
        let mut counts: Vec<u32> = vec![0; backbone_size];
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

        // Fill pass: the reference's appends, in the reference's order.
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

        Ok(BlastAaLookupTable {
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
        })
    })
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
/// EXPERIMENT (LOSAT_X_BXLUT / LOSAT_X_BXLUTSHADOW): the switch's dispatch for
/// `build_ncbi_lookup_from_prepared`. `Ok(table)` when this module produced the table (with
/// the shadow switch: the reference table, after comparing it with the in-place one);
/// `Err(prepared)` hands the input back to the reference builder (switch off, or not
/// applicable).
pub(super) fn x_dispatch(
    prepared: PreparedLookupQuery,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
) -> std::result::Result<BlastAaLookupTable, PreparedLookupQuery> {
    match x_bxlut_mode() {
        0 => Err(prepared),
        1 => x_build_direct(prepared, threshold, matrix, word_length),
        _ => {
            let copy = PreparedLookupQuery {
                concat_query: prepared.concat_query.clone(),
                lookup_locations: prepared.lookup_locations.clone(),
                frame_bases: prepared.frame_bases.clone(),
                contexts: prepared.contexts.clone(),
                skipped_seg_mask: prepared.skipped_seg_mask,
            };
            let reference = build_lookup_from_prepared(prepared, threshold, matrix, word_length).0;
            if let Ok(direct) = x_build_direct(copy, threshold, matrix, word_length) {
                x_assert_same_table(&direct, &reference);
                X_SHADOW_TABLES.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(reference)
        }
    }
}

/// EXPERIMENT (LOSAT_X_BXLUTSHADOW): tables compared so far.
pub(crate) static X_SHADOW_TABLES: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BXLUTSHADOW): `Some(tables compared)` in shadow mode, for the exit summary.
pub(crate) fn x_shadow_tables() -> Option<u64> {
    (x_bxlut_mode() == 2).then(|| X_SHADOW_TABLES.load(std::sync::atomic::Ordering::Relaxed))
}

/// EXPERIMENT (LOSAT_X_BXLUTSHADOW): assert that two tables are equal field by field.
pub(super) fn x_assert_same_table(direct: &BlastAaLookupTable, reference: &BlastAaLookupTable) {
    assert_eq!(
        direct.backbone.len(),
        reference.backbone.len(),
        "LOSAT_X_BXLUTSHADOW: backbone size differs"
    );
    for (idx, (a, b)) in direct.backbone.iter().zip(&reference.backbone).enumerate() {
        assert!(
            a.num_used == b.num_used && a.entries == b.entries,
            "LOSAT_X_BXLUTSHADOW: backbone cell {idx} differs ({} {:?} vs {} {:?})",
            a.num_used,
            a.entries,
            b.num_used,
            b.entries
        );
    }
    assert!(
        direct.overflow == reference.overflow,
        "LOSAT_X_BXLUTSHADOW: overflow array differs"
    );
    assert!(direct.pv == reference.pv, "LOSAT_X_BXLUTSHADOW: pv differs");
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
        "LOSAT_X_BXLUTSHADOW: table header differs"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>, frame: i8) -> QueryFrame {
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

    fn residues(seed: &mut u64, n: usize) -> Vec<u8> {
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
}
