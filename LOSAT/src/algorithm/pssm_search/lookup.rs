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
/// it, in the order they were added. A table is indexed once per query (`build_pssm_lookup_in`);
/// makeprofiledb indexes every profile of a volume into one backbone, each at its offset in the
/// concatenated profiles (`query_bias`) and with its own threshold
/// (`makeprofiledb::rps_files`).
pub(crate) struct ThinBackbone {
    chains: Vec<Vec<i32>>,
    word_length: usize,
}

/// `ilog2(BLASTAA_SIZE) + 1`: the bits of a letter in a word's index.
const CHARSIZE: usize = 5;

impl ThinBackbone {
    /// `BlastAaLookupTableNew`'s backbone: one empty chain per word index of `word_length`
    /// letters of the 28-letter alphabet.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:225-247
    pub(crate) fn new(word_length: usize) -> Self {
        let mut backbone_size = 0usize;
        for i in 0..word_length {
            backbone_size |= (BLASTAA_SIZE - 1) << (i * CHARSIZE);
        }
        backbone_size += 1;
        ThinBackbone {
            // Built through an iterator of this module rather than `vec![...; n]`, so that the PSSM
            // table shares no `SpecFromElem` instance with blastp's and tblastx's builders (design
            // §4.3: no generic helper shared with the default path).
            chains: (0..backbone_size).map(|_| Vec::new()).collect(),
            word_length,
        }
    }

    /// `BlastAaLookupIndexQuery` of a PSSM (`s_AddPSSMNeighboringWords`): the words of `rows`
    /// over the lookup segments `locations`, each stored at its offset plus `query_bias`, kept
    /// when their score reaches `threshold` (`lookup->threshold`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:413-426
    /// NCBI's threshold and the recursion's sums are `Int4`. A search's threshold is an `Int4` and
    /// its PSSM scores stay within `BLAST_SCORE_MIN..=BLAST_SCORE_MAX` (`in_pssm`, the PSSM engine),
    /// so their sums never leave `Int4`; makeprofiledb's profiles may hold any `Int4` score and
    /// its scaled threshold any value, where NCBI's sums overflow and its threshold's conversion
    /// is undefined: LOSAT adds and compares in 64 bits, the mathematical values (a proposed
    /// exception of `PD-LOSAT-NCBI-DEFECTS`, psiblast line DECISIONS 80).
    pub(crate) fn add_pssm_words<R: PssmRows>(
        &mut self,
        rows: &R,
        query_bias: i32,
        threshold: i64,
        locations: &[(i64, i64)],
    ) {
        let word_length = self.word_length;
        add_pssm_neighboring_words(
            self,
            rows,
            word_length,
            CHARSIZE,
            threshold,
            query_bias,
            locations,
        );
    }

    /// `BlastAaLookupFinalize` with `eBackbone` of the backbone, for a table of `query_length`
    /// letters searched with `threshold` (the values the table keeps).
    pub(crate) fn finalize(
        self,
        query_length: usize,
        threshold: i32,
    ) -> BlastAaLookupTable<BackboneCell> {
        let word_length = self.word_length;
        let mask = (1usize << (word_length * CHARSIZE)) - 1;
        finalize(self, query_length, word_length, CHARSIZE, mask, threshold)
    }
}

/// What the recursion of `s_AddPSSMWordHitsCore` reads (`NeighborInfo`).
struct PssmNeighbors<'a, R: PssmRows> {
    rows: &'a R,
    /// The row of the word's first letter in `rows`.
    offset: usize,
    /// The offset stored for the word (`info->query_bias`: the offset plus the query's bias).
    stored_offset: i32,
    wordsize: usize,
    charsize: usize,
    threshold: i64,
    /// `row_max` of the word's `wordsize` positions.
    row_max: &'a [i64],
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
    let mut thin = ThinBackbone::new(word_length);
    thin.add_pssm_words(rows, 0, i64::from(threshold), locations);
    thin.finalize(rows.len(), threshold)
}

/// `s_AddPSSMNeighboringWords` over the lookup segments `locations`, the offsets stored plus
/// `query_bias`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:623-658
fn add_pssm_neighboring_words<R: PssmRows>(
    thin: &mut ThinBackbone,
    rows: &R,
    wordsize: usize,
    charsize: usize,
    threshold: i64,
    query_bias: i32,
    locations: &[(i64, i64)],
) {
    let row_max_of = |pos: usize| -> i64 {
        let row = rows.row(pos);
        let mut max = row[0];
        for &score in &row[1..BLASTAA_SIZE] {
            max = max.max(score);
        }
        i64::from(max)
    };
    let mut row_max = vec![0i64; wordsize];
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
            add_pssm_word_hits(
                thin,
                rows,
                (offset, query_bias),
                wordsize,
                charsize,
                threshold,
                &row_max,
            );
            for i in 0..wordsize - 1 {
                row_max[i] = row_max[i + 1];
            }
        }
    }
}

/// `s_AddPSSMWordHits` of the word at `offset` of `rows`, stored at `offset + query_bias`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:681-693
fn add_pssm_word_hits<R: PssmRows>(
    thin: &mut ThinBackbone,
    rows: &R,
    (offset, query_bias): (usize, i32),
    wordsize: usize,
    charsize: usize,
    threshold: i64,
    row_max: &[i64],
) {
    let mut info = PssmNeighbors {
        rows,
        offset,
        stored_offset: offset as i32 + query_bias,
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
    score: i64,
    current_pos: usize,
) {
    let score = score - info.row_max[current_pos];
    let row = info.rows.row(info.offset + current_pos);
    if current_pos == info.wordsize - 1 {
        for (i, &cell) in row.iter().enumerate() {
            if score + i64::from(cell) >= info.threshold {
                info.subject_word[current_pos] = i as u8;
                blast_lookup_add_word_hit(
                    thin,
                    info.wordsize,
                    info.charsize,
                    &info.subject_word,
                    info.stored_offset,
                );
            }
        }
        return;
    }
    for (i, &cell) in row.iter().enumerate() {
        let cell = i64::from(cell);
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
