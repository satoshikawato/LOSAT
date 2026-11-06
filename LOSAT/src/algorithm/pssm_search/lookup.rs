//! The lookup table of a PSSM search: the neighbouring words of each query position scored
//! against the position's PSSM rows (no exact-word rule), finalized into the table that blastp's
//! subject scan reads (`BlastAaLookupTable<BackboneCell>`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:84-99
//! ```c
//!    case eAaLookupTable:
//!        {
//!        Int4** matrix = NULL;
//!        Boolean has_pssm = FALSE;
//!        if (sbp->psi_matrix && sbp->psi_matrix->pssm) {
//!            matrix = sbp->psi_matrix->pssm->data;
//!            has_pssm = TRUE;
//!        } else {
//!            matrix = sbp->matrix->data;
//!        }
//!        BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
//!                              &lookup_wrap->lut);
//!        ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
//!        BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
//!                                  query, lookup_segments, 0);
//!        /* if query length less than 64k, we can save cache by using small bone */
//!        bone_type = ( query->length >= INT2_MAX*2) ? eBackbone: eSmallbone;
//!        BlastAaLookupFinalize((BlastAaLookupTable*) lookup_wrap->lut, bone_type);
//! ```
//! NCBI's small backbone (`Uint2` cells, for queries below 65,534 letters) holds the same hits in
//! the same order as the large one; LOSAT keeps `Int4` cells (`BackboneCell`) for every length,
//! the table blastp's scan reads (`x_seed_fused::x_scan_two_hit`). Storage width only.

use super::rows::PssmRows;
use crate::algorithm::tblastx::lookup::{BackboneCell, BlastAaLookupTable, AA_HITS_PER_CELL};
use crate::core::blast_util::BLASTAA_SIZE;

/// NCBI's `PV_ARRAY_BTS` (blast_lookup.h:41-44): 32 presence bits per word.
const PV_ARRAY_BTS: usize = 5;
const PV_ARRAY_MASK: usize = 31;

/// The thin backbone (`lookup->thin_backbone`): for each word index the query offsets added to
/// it, in the order they were added.
struct ThinBackbone {
    chains: Vec<Vec<i32>>,
}

/// What the recursion of `s_AddPSSMWordHitsCore` reads (`NeighborInfo`).
struct PssmNeighbors<'a, R: PssmRows> {
    rows: &'a R,
    /// The query offset of the word (`info->query_bias`).
    offset: usize,
    wordsize: usize,
    charsize: usize,
    threshold: i32,
    /// `row_max` of the word's `wordsize` positions.
    row_max: &'a [i32],
    subject_word: [u8; 32],
}

/// The lookup table of `rows` for a query of one context (offset 0) indexed whole: NCBI's
/// `BlastAaLookupTableNew`, `s_AddPSSMNeighboringWords` over the lookup segment `[0, len - 1]`
/// and `BlastAaLookupFinalize`. `threshold` is `(Int4)opt->threshold`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:236-247
/// ```c
///     lookup->charsize = ilog2(BLASTAA_SIZE) + 1;
///     lookup->word_length = opt->word_size;
///
///     for (i = 0; i < lookup->word_length; i++)
///         lookup->backbone_size |= (BLASTAA_SIZE - 1) << (i * lookup->charsize);
///     lookup->backbone_size++;
///
///     lookup->mask = (1 << (opt->word_size * lookup->charsize)) - 1;
///     lookup->alphabet_size = BLASTAA_SIZE;
///     lookup->threshold = (Int4)opt->threshold;
/// ```
pub(crate) fn build_pssm_lookup<R: PssmRows>(
    rows: &R,
    word_length: usize,
    threshold: i32,
) -> BlastAaLookupTable<BackboneCell> {
    // ilog2(28) + 1
    let charsize = 5usize;
    let mut backbone_size = 0usize;
    for i in 0..word_length {
        backbone_size |= (BLASTAA_SIZE - 1) << (i * charsize);
    }
    backbone_size += 1;
    let mask = (1usize << (word_length * charsize)) - 1;
    let mut thin = ThinBackbone {
        // Built through an iterator of this module rather than `vec![...; n]`, so that the PSSM
        // table shares no `SpecFromElem` instance with blastp's and tblastx's builders (design
        // §4.3: no generic helper shared with the default path).
        chains: (0..backbone_size).map(|_| Vec::new()).collect(),
    };
    add_pssm_neighboring_words(&mut thin, rows, word_length, charsize, threshold);
    finalize(thin, rows.len(), word_length, charsize, mask, threshold)
}

/// `s_AddPSSMNeighboringWords` for the one location `[0, len - 1]`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:623-658
/// ```c
///     for (loc = location; loc; loc = loc->next) {
///         Int4 from = loc->ssr->left;
///         Int4 to = loc->ssr->right - wordsize + 1;
///         Int4 **row = matrix + from;
///
///         /* prepare to start another run of adjacent query words. Find the
///            maximum possible score for the first wordsize-1 rows of the PSSM */
///
///         if (to >= from) {
///             for (i = 0; i < wordsize - 1; i++) {
///                 row_max[i] = row[i][0];
///                 for (j = 1; j < lookup->alphabet_size; j++)
///                     row_max[i] = MAX(row_max[i], row[i][j]);
///             }
///         }
///
///         for (offset = from; offset <= to; offset++, row++) {
///             /* find the maximum score of the next PSSM row */
///
///             row_max[wordsize - 1] = row[wordsize - 1][0];
///             for (i = 1; i < lookup->alphabet_size; i++)
///                 row_max[wordsize - 1] = MAX(row_max[wordsize - 1],
///                                             row[wordsize - 1][i]);
///
///             /* find all neighboring words */
///
///             s_AddPSSMWordHits(lookup, row, offset + query_bias, row_max);
///
///             /* shift the list of maximum scores over by one, to make room for
///                the next maximum in the next loop iteration */
///
///             for (i = 0; i < wordsize - 1; i++)
///                 row_max[i] = row_max[i + 1];
///         }
///     }
/// ```
fn add_pssm_neighboring_words<R: PssmRows>(
    thin: &mut ThinBackbone,
    rows: &R,
    wordsize: usize,
    charsize: usize,
    threshold: i32,
) {
    let len = rows.len();
    if len < wordsize {
        return;
    }
    let row_max_of = |pos: usize| -> i32 {
        let row = rows.row(pos);
        let mut max = row[0];
        for &score in &row[1..BLASTAA_SIZE] {
            max = max.max(score);
        }
        max
    };
    let mut row_max = vec![0i32; wordsize];
    for (i, slot) in row_max.iter_mut().enumerate().take(wordsize - 1) {
        *slot = row_max_of(i);
    }
    for offset in 0..=len - wordsize {
        row_max[wordsize - 1] = row_max_of(offset + wordsize - 1);
        add_pssm_word_hits(thin, rows, offset, wordsize, charsize, threshold, &row_max);
        for i in 0..wordsize - 1 {
            row_max[i] = row_max[i + 1];
        }
    }
}

/// `s_AddPSSMWordHits`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:681-692
/// ```c
///     /* compute the largest possible score that any neighboring word can have;
///        this maximum will gradually be replaced by exact scores as subject
///        words are built up */
///
///     score = row_max[0];
///     for (i = 1; i < lookup->word_length; i++)
///         score += row_max[i];
///
///     s_AddPSSMWordHitsCore(&info, score, 0);
/// ```
fn add_pssm_word_hits<R: PssmRows>(
    thin: &mut ThinBackbone,
    rows: &R,
    offset: usize,
    wordsize: usize,
    charsize: usize,
    threshold: i32,
    row_max: &[i32],
) {
    let mut info = PssmNeighbors {
        rows,
        offset,
        wordsize,
        charsize,
        threshold,
        row_max,
        subject_word: [0; 32],
    };
    let mut score = row_max[0];
    for &max in &row_max[1..wordsize] {
        score += max;
    }
    add_pssm_word_hits_core(thin, &mut info, score, 0);
}

/// `s_AddPSSMWordHitsCore`: the recursion over the alphabet at each word position, keeping the
/// words whose exact score reaches the threshold.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:695-748
/// ```c
///     score -= info->row_max[current_pos];
///     row = info->matrix[current_pos];
///
///     if (current_pos == info->wordsize - 1) {
///         ...
///         for (i = 0; i < alphabet_size; i++) {
///             if (score + row[i] >= threshold) {
///                 subject_word[current_pos] = i;
///                 BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
///                                       charsize, subject_word, offset);
///             }
///         }
///         return;
///     }
///
///     for (i = 0; i < alphabet_size; i++) {
///         if (score + row[i] >= threshold) {
///             subject_word[current_pos] = i;
///             s_AddPSSMWordHitsCore(info, score + row[i], current_pos + 1);
///         }
///     }
/// ```
/// `info->matrix` is the PSSM from the word's first position (`row = matrix + from` advanced to
/// the offset), so `matrix[current_pos]` is the row of query position `offset + current_pos`.
fn add_pssm_word_hits_core<R: PssmRows>(
    thin: &mut ThinBackbone,
    info: &mut PssmNeighbors<'_, R>,
    score: i32,
    current_pos: usize,
) {
    let score = score - info.row_max[current_pos];
    let row = info.rows.row(info.offset + current_pos);
    if current_pos == info.wordsize - 1 {
        for (i, &cell) in row.iter().enumerate() {
            if score + cell >= info.threshold {
                info.subject_word[current_pos] = i as u8;
                blast_lookup_add_word_hit(
                    thin,
                    info.wordsize,
                    info.charsize,
                    &info.subject_word,
                    info.offset as i32,
                );
            }
        }
        return;
    }
    for (i, &cell) in row.iter().enumerate() {
        if score + cell >= info.threshold {
            info.subject_word[current_pos] = i as u8;
            add_pssm_word_hits_core(thin, info, score + cell, current_pos + 1);
        }
    }
}

/// `BlastLookupAddWordHit`: the offset appended to the chain of the word's index.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:33-77
/// ```c
///     index = ComputeTableIndex(wordsize, charsize, seq);
///     ...
///     chain[chain[1] + 2] = query_offset;
///     chain[1]++;
/// ```
/// NCBI's chain header (its allocated size and count) is the `Vec`'s; the offsets are the same.
fn blast_lookup_add_word_hit(
    thin: &mut ThinBackbone,
    wordsize: usize,
    charsize: usize,
    seq: &[u8],
    query_offset: i32,
) {
    let mut index = 0usize;
    for &residue in &seq[..wordsize] {
        index = (index << charsize) | residue as usize;
    }
    thin.chains[index].push(query_offset);
}

/// `BlastAaLookupFinalize` with the large backbone (`Int4` cells): the presence vector, the cells
/// with up to `AA_HITS_PER_CELL` hits in place, the others in the overflow array in index order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:285-366
/// ```c
///     /* find out how many cells need the overflow array */
///     for (i = 0; i < lookup->backbone_size; i++) {
///         if (lookup->thin_backbone[i]) {
///             if (lookup->thin_backbone[i][1] > AA_HITS_PER_CELL){
///                 overflow_cells_needed += lookup->thin_backbone[i][1];
///             }
///             if (lookup->thin_backbone[i][1] > longest_chain)
///                 longest_chain = lookup->thin_backbone[i][1];
///         }
///     }
///     ...
///       for (i = 0; i < lookup->backbone_size; i++) {
///         /* if there are hits there, */
///         if (lookup->thin_backbone[i] ) {
///             Int4 * dest = NULL;
///             /* set the corresponding bit in the pv_array */
///             PV_SET(pv, i, PV_ARRAY_BTS);
///             bbc[i].num_used = lookup->thin_backbone[i][1];
///             /* if there are three or fewer hits, */
///             if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
///                 /* copy them into the thick_backbone cell */
///                 dest = bbc[i].payload.entries;
///             else /* more than three hits; copy to overflow array */
///             {
///                 bbc[i].payload.overflow_cursor = overflow_cursor;
///                 dest = (Int4 *) lookup->overflow;
///                 dest += overflow_cursor;
///                 overflow_cursor += lookup->thin_backbone[i][1];
///             }
///             for (j=0; j <lookup->thin_backbone[i][1]; j++)
///                 dest[j] = lookup->thin_backbone[i][j + 2];
/// ```
fn finalize(
    thin: ThinBackbone,
    query_length: usize,
    word_length: usize,
    charsize: usize,
    mask: usize,
    threshold: i32,
) -> BlastAaLookupTable<BackboneCell> {
    let backbone_size = thin.chains.len();
    let mut overflow_cells_needed = 0usize;
    let mut longest_chain = 0usize;
    for chain in &thin.chains {
        if chain.len() > AA_HITS_PER_CELL {
            overflow_cells_needed += chain.len();
        }
        longest_chain = longest_chain.max(chain.len());
    }
    // Built through an iterator of this module (see `ThinBackbone`'s chains).
    let mut backbone: Vec<BackboneCell> = (0..backbone_size)
        .map(|_| BackboneCell::default())
        .collect();
    let mut pv = vec![0u32; (backbone_size >> PV_ARRAY_BTS) + 1];
    let mut overflow = vec![0i32; overflow_cells_needed];
    let mut overflow_cursor = 0usize;
    for (i, chain) in thin.chains.into_iter().enumerate() {
        if chain.is_empty() {
            continue;
        }
        pv[i >> PV_ARRAY_BTS] |= 1u32 << (i & PV_ARRAY_MASK);
        backbone[i].num_used = chain.len() as i32;
        if chain.len() <= AA_HITS_PER_CELL {
            backbone[i].entries[..chain.len()].copy_from_slice(&chain);
        } else {
            backbone[i].entries[0] = overflow_cursor as i32;
            overflow[overflow_cursor..overflow_cursor + chain.len()].copy_from_slice(&chain);
            overflow_cursor += chain.len();
        }
    }
    BlastAaLookupTable {
        backbone,
        overflow,
        pv,
        // One context, at offset 0 of the query.
        frame_bases: vec![0],
        num_contexts: 1,
        query_length: query_length as i32,
        word_length: word_length as i32,
        alphabet_size: BLASTAA_SIZE as i32,
        charsize: charsize as i32,
        mask: mask as i32,
        longest_chain: longest_chain as i32,
        threshold,
        // NCBI's PSSM table keeps no per-letter row maximum (`row_max` is the recursion's).
        row_max: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::pssm_search::rows::QueryPssm;

    fn pssm(rows: Vec<[i32; BLASTAA_SIZE]>) -> QueryPssm {
        QueryPssm { rows }
    }

    /// The hits of a word (NCBIstdaa letters) in a finalized table, through its presence bit.
    fn hits(table: &BlastAaLookupTable<BackboneCell>, word: &[u8]) -> Vec<i32> {
        let mut index = 0usize;
        for &letter in word {
            index = (index << 5) | letter as usize;
        }
        if table.pv[index >> PV_ARRAY_BTS] & (1 << (index & PV_ARRAY_MASK)) == 0 {
            return Vec::new();
        }
        table.get_hits(index).to_vec()
    }

    // A PSSM of rows that score one letter per position high: each position's word of three
    // letters reaches the threshold only with those letters, so the table holds exactly the
    // query's words (the offsets in position order), and a word whose sum misses the threshold
    // is not added (no exact-word rule).
    #[test]
    fn words_reaching_the_threshold_are_added_in_position_order() {
        let mut rows = vec![[-4i32; BLASTAA_SIZE]; 6];
        // A C D A C D (NCBIstdaa 1, 3, 4)
        for (pos, letter) in [1usize, 3, 4, 1, 3, 4].into_iter().enumerate() {
            rows[pos][letter] = 5;
        }
        let rows = pssm(rows);
        let table = build_pssm_lookup(&rows, 3, 11);
        assert_eq!(hits(&table, &[1, 3, 4]), vec![0, 3]);
        assert_eq!(hits(&table, &[3, 4, 1]), vec![1]);
        assert_eq!(hits(&table, &[4, 1, 3]), vec![2]);
        assert_eq!(table.longest_chain, 2);
        // 15 >= 11, but one other letter gives 5 + 5 - 4 = 6 < 11.
        assert!(hits(&table, &[1, 3, 3]).is_empty());
        // Threshold 16: no word reaches it (the best is 15), not even the query's own.
        let table = build_pssm_lookup(&rows, 3, 16);
        assert!(table.backbone.iter().all(|cell| cell.num_used == 0));
    }

    // Chains of more than three offsets go to the overflow array in index order, with the cell's
    // cursor in its first entry.
    #[test]
    fn long_chains_go_to_the_overflow_array() {
        let mut rows = vec![[-4i32; BLASTAA_SIZE]; 7];
        for row in rows.iter_mut() {
            row[1] = 5;
        }
        let rows = pssm(rows);
        let table = build_pssm_lookup(&rows, 3, 15);
        assert_eq!(hits(&table, &[1, 1, 1]), vec![0, 1, 2, 3, 4]);
        assert_eq!(table.overflow, vec![0, 1, 2, 3, 4]);
        assert_eq!(table.longest_chain, 5);
        assert_eq!(table.query_length, 7);
    }

    /// The de Bruijn sequence of order 3 over the 28 letters, unwrapped (`debruijn` of the unit
    /// test below): every word of 3 letters appears once.
    fn debruijn3() -> Vec<u8> {
        fn db(t: usize, p: usize, a: &mut [usize; 7], out: &mut Vec<u8>) {
            if t > 3 {
                if 3 % p == 0 {
                    out.extend(a[1..=p].iter().map(|&letter| letter as u8));
                }
                return;
            }
            a[t] = a[t - p];
            db(t + 1, p, a, out);
            for letter in a[t - p] + 1..BLASTAA_SIZE {
                a[t] = letter;
                db(t + 1, t, a, out);
            }
        }
        let mut out = Vec::new();
        db(1, 1, &mut [0; 7], &mut out);
        let wrap: Vec<u8> = out[..2].to_vec();
        out.extend(wrap);
        out
    }

    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:310-421 testDebruijnPSSM
    // NCBI compiles the case out (`#if 0`, "Needs to be fixed to actually use a PSSM"); this port
    // runs it on the PSSM table: the query is the de Bruijn sequence of order 3 (every word once),
    // the PSSM row of each position BLOSUM62's column of its letter (`posMatrix[i][j] =
    // matrix[j][sequence[i+1]]`), threshold 11, word size 3. Every word's cell holds exactly the
    // query words whose score with it reaches 11 (the brute-force count of the test): the PSSM
    // table has no exact-word rule.
    #[test]
    fn debruijn_pssm_cells_hold_the_brute_force_neighbours() {
        let sequence = debruijn3();
        assert_eq!(sequence.len(), BLASTAA_SIZE.pow(3) + 2);
        let matrix = |q: usize, s: usize| crate::core::blast_util::blosum62_score(q as u8, s as u8);
        let rows = pssm(
            sequence
                .iter()
                .map(|&letter| {
                    let mut row = [0i32; BLASTAA_SIZE];
                    for (j, cell) in row.iter_mut().enumerate() {
                        *cell = matrix(j, letter as usize);
                    }
                    row
                })
                .collect(),
        );
        let table = build_pssm_lookup(&rows, 3, 11);
        // count(i, j, k) = #{(x, y, z): m[i][x] + m[j][y] + m[k][z] >= 11}, by the number of z of
        // each row k reaching each threshold.
        let sorted: Vec<Vec<i32>> = (0..BLASTAA_SIZE)
            .map(|k| {
                let mut row: Vec<i32> = (0..BLASTAA_SIZE).map(|z| matrix(k, z)).collect();
                row.sort_unstable();
                row
            })
            .collect();
        let reaching = |k: usize, t: i32| BLASTAA_SIZE - sorted[k].partition_point(|&v| v < t);
        for i in 0..BLASTAA_SIZE {
            for j in 0..BLASTAA_SIZE {
                let mut pairs = Vec::with_capacity(BLASTAA_SIZE * BLASTAA_SIZE);
                for x in 0..BLASTAA_SIZE {
                    for y in 0..BLASTAA_SIZE {
                        pairs.push(matrix(i, x).saturating_add(matrix(j, y)));
                    }
                }
                for k in 0..BLASTAA_SIZE {
                    let count: usize = pairs
                        .iter()
                        .map(|&pair| reaching(k, 11i32.saturating_sub(pair)))
                        .sum();
                    let index = (i << 10) | (j << 5) | k;
                    assert_eq!(
                        table.backbone[index].num_used as usize, count,
                        "word ({i}, {j}, {k})"
                    );
                }
            }
        }
    }
}
