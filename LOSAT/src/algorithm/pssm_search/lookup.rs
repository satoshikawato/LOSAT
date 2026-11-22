//! The lookup table of a PSSM search: the neighbouring words of each query position scored
//! against the position's PSSM rows (no exact-word rule), finalized into the table that blastp's
//! subject scan reads (`BlastAaLookupTable<BackboneCell>`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:83-100
//! NCBI's small backbone (`Uint2` cells, for queries below 65,534 letters) holds the same hits in
//! the same order as the large one; LOSAT keeps `Int4` cells (`BackboneCell`) for every length,
//! the table blastp's scan reads (`x_seed_fused::x_scan_two_hit`). Storage width only.

use super::rows::PssmRows;
use crate::algorithm::tblastx::lookup::{BackboneCell, BlastAaLookupTable, AA_HITS_PER_CELL};
use crate::core::blast_util::BLASTAA_SIZE;

/// NCBI's `PV_ARRAY_BTS` (c++/include/algo/blast/core/blast_lookup.h:41-44): 32 presence bits per word.
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
/// and `BlastAaLookupFinalize`. `threshold` is `opt->threshold` cast to `Int4`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:236-247
#[cfg(test)]
pub(crate) fn build_pssm_lookup<R: PssmRows>(
    rows: &R,
    word_length: usize,
    threshold: i32,
) -> BlastAaLookupTable<BackboneCell> {
    build_pssm_lookup_in(rows, word_length, threshold, &[(0, rows.len() as i64 - 1)])
}

/// `build_pssm_lookup` over the lookup segments `locations` (NCBI's `lookup_segments`, inclusive
/// `[left, right]` ranges of query offsets; one with `left > right` adds no word): the
/// complement of the query's filtering locations, whose masks (`-seg`, `-lcase_masking`, also
/// with `-soft_masking true`) keep words out of the table: psiblast's PSSM searches
/// (`blastp::blast_engine::pssm_hooks::pssm_lookup`, blastp's lookup segments) and tblastn
/// `-in_pssm` (`algorithm::tblastn::in_pssm`).
pub(crate) fn build_pssm_lookup_in<R: PssmRows>(
    rows: &R,
    word_length: usize,
    threshold: i32,
    locations: &[(i64, i64)],
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
    add_pssm_neighboring_words(&mut thin, rows, word_length, charsize, threshold, locations);
    finalize(thin, rows.len(), word_length, charsize, mask, threshold)
}

/// `s_AddPSSMNeighboringWords` over the lookup segments `locations`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:623-658
fn add_pssm_neighboring_words<R: PssmRows>(
    thin: &mut ThinBackbone,
    rows: &R,
    wordsize: usize,
    charsize: usize,
    threshold: i32,
    locations: &[(i64, i64)],
) {
    let row_max_of = |pos: usize| -> i32 {
        let row = rows.row(pos);
        let mut max = row[0];
        for &score in &row[1..BLASTAA_SIZE] {
            max = max.max(score);
        }
        max
    };
    let mut row_max = vec![0i32; wordsize];
    for &(from, right) in locations {
        let to = right - wordsize as i64 + 1;
        if to < from {
            continue;
        }
        let (from, to) = (from as usize, to as usize);
        for (i, slot) in row_max.iter_mut().enumerate().take(wordsize - 1) {
            *slot = row_max_of(from + i);
        }
        for offset in from..=to {
            row_max[wordsize - 1] = row_max_of(offset + wordsize - 1);
            add_pssm_word_hits(thin, rows, offset, wordsize, charsize, threshold, &row_max);
            for i in 0..wordsize - 1 {
                row_max[i] = row_max[i + 1];
            }
        }
    }
}

/// `s_AddPSSMWordHits`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:681-693
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
/// `info->matrix` is the PSSM from the word's first position (`matrix` advanced by `from`, then to
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:282-366
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
            backbone[i].set_cursor(overflow_cursor);
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
#[path = "lookup_tests.rs"]
mod tests;
