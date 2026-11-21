//! SC10b PR-3: the group tables of the TBLASTX unit path (`x_pair_par.rs`) for words whose backbone
//! index does not fit 16 bits (`-word_size 4`: 20 bits), built for all the batch's tables at once,
//! with 32-bit indices, by counting and filling on the search pool. BLASTX's table of 4-letter words
//! is its table of one list of all the query's lookup locations (`x_build_full_wide`, item SC9 of
//! line NP), and so are TBLASTX's full table of a batch and TBLASTN's table of 4-letter words
//! (`backbone.rs` `x_full_wide_table`, item X2e of line RG). A child module of
//! `x_pair_lut.rs`, whose `x_build_group` (16-bit neighbour lists) keeps building the tables of
//! words of 2 and 3 letters.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:461-467
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-590
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:74-76
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
//! 1. Sources. A stable sort of the words of L_0, then L_1, .., by source cell s (the loop of
//!    `blast_lookup_index_query_exact_matches`, `x_wide_exact_words`; a radix sort of the words'
//!    order, `x_wide_radix_order`): the sources with a word, ascending, and each source's offsets
//!    ordered by (table, order of the table's loop). Table t's part of source s is therefore one
//!    run, E_t(s), in its order (x_pair_lut.rs item 2: each location is indexed on its own, so the
//!    loop over L_t adds exactly the words of L_t, in L_t's order).
//! 2. Lists. T(s) of each source, once for all tables, on the pool over ranges of sources: the
//!    statements of `x_target_list` (the neighbours of the `s_AddWordHitsCore` recursion, in its
//!    order, which is ascending: `x_wide_neighbors`; and the exact cell when s_AddWordHits adds it,
//!    where its index falls among them) with 32-bit indices. T(s) reads only the residues of s, the
//!    matrix, the threshold and the word length (x_lut_direct.rs header), so it is read at s's first
//!    offset, as the reference reads it (blast_aalookup.c:490). The cells of T(s) are distinct, so
//!    the place of the exact cell changes no chain. A range keeps its lists by stripe of the backbone
//!    (2^b cells, 32 to 2^16; step 3): for stripe j, its sources with cells in the stripe, ascending,
//!    each with those cells (their index in the stripe, a `u16`), one block of T(s) since T(s) is
//!    ascending. A range also takes, by table, the least and greatest offset of its sources with a
//!    target.
//! 3. The tables, all of them at once, by stripe (the stripes side by side): (a) each stripe counts
//!    the appends to each of its cells in each table, `count_t[w]` = sum of |E_t(s)| over the s with
//!    w in T(s) (the reference's chain length), from every range's lists for the stripe, and sums
//!    each table's overflow size, longest count and used cells in the stripe; (b) each table's
//!    cells and overflow go to the stripes in cell order (`BlastAaLookupFinalize` assigns the
//!    overflow cursors in ascending cell order, and a compacted table keeps its cells in index
//!    order); (c) each stripe, from its counts of (a) (zeroed by the stripe there, kept since),
//!    writes its cells (`num_used`, the pv bit, the cursor of an overflow cell; a cell without
//!    entries all zero, as after the C `calloc`), turns `count_t[w]` into the first position of cell
//!    w, then walks the lists as (a) did and appends, for each s with w in T(s), E_t(s) at the
//!    position of w in table t, advanced by the offsets written.
//!
//! # Why the tables are those of `x_build_group`
//!
//! A cell is one stripe's, and the stripe appends to it in the order of the ranges, each range's
//! sources ascending, so cell w of table t holds, for s ascending with w in T(s), E_t(s) in order:
//! the chain of x_pair_lut.rs item 5. The positions of a cell are written once each (a cell's run
//! starts at its first position and takes `count_t[w]` offsets), and the cells' runs of the overflow
//! are disjoint. `num_used`, the pv bits, the inline entries or the overflow cursor, the overflow
//! array and `longest_chain` are `BlastAaLookupFinalize`'s function of these chains; the compaction
//! keeps the cells with entries in index order (`x_scan_bf::x_compact`) with the rank of the
//! presence vector (`x_scan_bf::x_rank`); `q_min` and `q_max` are the least and greatest offset of a
//! source word with a target, as in `x_build_group`; the header fields are the expressions of
//! `x_build_group`. Only the scheduling (ranges, stripes, the pool's threads) differs, and it changes
//! no value.
//!
//! Memory: the tables, the lists T(s) of the batch's sources by stripe (2 bytes a cell of a list
//! and 8 a block, freed when the tables are built; no memo is kept across batches: at 4 letters the
//! lists of a batch take tens to hundreds of megabytes), the counts of step 3 (4 bytes a backbone cell and
//! table: TBLASTX's unit path builds at most 6 tables at once, x_pair_par.rs `X_MAX_GROUPS`; a full
//! table, of BLASTX, TBLASTN or TBLASTX's reference frame loop, is one), and
//! for a while the sources' words in step 1 (20 bytes a query word) and a range's lists as made
//! before they are put by stripe (step 2). No chain is made: the reference builder holds the chains
//! and the table at once.

use super::*;

/// No NCBI counterpart: the letter width of the lookup alphabet (`charsize`, ilog2(28) + 1), as a
/// constant.
const X_WIDE_CHARSIZE: usize = LOOKUP_ALPHABET_SIZE.ilog2() as usize + 1;

/// No NCBI counterpart: the widest backbone index the joint builder takes, in bits (words of 4
/// letters, the longest TBLASTX takes: `Word-size must be less than 6`, checked at 5).
const X_WIDE_INDEX_BITS: usize = 20;

/// No NCBI counterpart: the bits of a digit of the radix sort of the sources (`x_wide_radix_order`):
/// two passes over 1,024 counts for 4-letter words.
const X_WIDE_RADIX_BITS: usize = 10;

/// No NCBI counterpart: stripes of the backbone per pool thread (step 3; the pool's threads claim
/// the stripes in turn, and the stripes' work differs with the cells' counts).
const X_WIDE_STRIPES_PER_THREAD: usize = 2;

/// No NCBI counterpart: the widest stripe, in bits of its cells (a cell's index in its stripe is a
/// `u16`).
const X_WIDE_STRIPE_BITS_MOST: u32 = 16;

/// No NCBI counterpart: ranges of sources per pool thread for the lists T(s) (they differ in
/// length; the pool's threads claim the ranges in turn).
const X_WIDE_LIST_PARTS_PER_THREAD: usize = 4;

/// No NCBI counterpart: the bit of a position in a stripe that says the position is in the stripe's
/// part of the overflow array; without it the position is `(cell's slot in the stripe << 2) | entry`.
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
/// `info->row_max` of s_AddWordHits), as arrays of the alphabet's size, and the letters of each row
/// that score at least a floor, one bit each (`x_wide_neighbors`): every letter for a floor up to
/// the matrix's least score `least` (the gap letter's `BLAST_SCORE_MIN` in NCBI's matrices, far below
/// the others), else for the floors `low + v`, v from 0 to `ladder - 1`, with `low` the least score
/// above `least` and `low + ladder - 1` the greatest.
struct XWideScores {
    rows: [[i32; LOOKUP_ALPHABET_SIZE]; LOOKUP_ALPHABET_SIZE],
    row_max: [i32; LOOKUP_ALPHABET_SIZE],
    least: i32,
    low: i32,
    ladder: usize,
    /// Row a's letters scoring at least `low + v`: `at_least[a * ladder + v]`.
    at_least: Vec<u32>,
}

const _: () = assert!(
    LOOKUP_ALPHABET_SIZE <= 32,
    "a row's letters are bits of a u32"
);

/// No NCBI counterpart: every letter of the lookup alphabet, one bit each.
const X_WIDE_ALL_LETTERS: u32 = ((1u64 << LOOKUP_ALPHABET_SIZE) - 1) as u32;

impl XWideScores {
    /// The rows of `matrix`, the row maxima `row_max`, and each row's letters by floor.
    fn new(matrix: ScoringMatrix, row_max: &[i32]) -> Self {
        let mut rows = [[0; LOOKUP_ALPHABET_SIZE]; LOOKUP_ALPHABET_SIZE];
        for (a, row) in rows.iter_mut().enumerate() {
            for (b, score) in row.iter_mut().enumerate() {
                *score = lookup_matrix_score(matrix, a, b);
            }
        }
        let scores = || rows.iter().flatten().copied();
        let least = scores().min().unwrap_or(0);
        let most = scores().max().unwrap_or(0);
        let low = scores().filter(|&s| s > least).min().unwrap_or(least);
        let ladder = usize::try_from(i64::from(most) - i64::from(low) + 1)
            .expect("a matrix's scores span a ladder");
        let mut at_least = vec![0u32; LOOKUP_ALPHABET_SIZE * ladder];
        for (a, row) in rows.iter().enumerate() {
            for (v, letters) in at_least[a * ladder..(a + 1) * ladder]
                .iter_mut()
                .enumerate()
            {
                let floor = low + v as i32;
                for (b, &score) in row.iter().enumerate() {
                    if score >= floor {
                        *letters |= 1 << b;
                    }
                }
            }
        }
        let mut scores = XWideScores {
            rows,
            row_max: [0; LOOKUP_ALPHABET_SIZE],
            least,
            low,
            ladder,
            at_least,
        };
        scores.row_max.copy_from_slice(row_max);
        scores
    }

    /// No NCBI counterpart: row a's letters i with `row[i] >= floor`, one bit each.
    #[inline(always)]
    fn letters_at_least(&self, a: usize, floor: i64) -> u32 {
        if floor <= i64::from(self.least) {
            return X_WIDE_ALL_LETTERS;
        }
        // Up to `low` every letter but those of the least score passes (the first rung), above the
        // greatest score none.
        let v = (floor - i64::from(self.low)).max(0);
        if v >= self.ladder as i64 {
            return 0;
        }
        self.at_least[a * self.ladder + v as usize]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:562-563,580-582,600-604
// The recursion of `x_record_neighbors` with the same tests in the same order; the letters chosen
// so far are carried as their index (`prefix`, ComputeTableIndex's shifts over `subject_word`), so a
// complete word's backbone index is `(prefix << charsize) | i`. Each level takes the letters that
// pass `score + row[i] >= threshold` from the row's bits of letters scoring at least
// `threshold - score` (`XWideScores::letters_at_least`, in 64 bits), lowest first: the C loop's
// letters, in its order. The level before the last runs the last level's loop itself (the same
// statements, without a call per letter). The words come out in ascending order of their index (the
// letters are its digits, highest first).
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
    let floor = |score: i32| i64::from(threshold) - i64::from(score);
    let query_residue = query_word[current_pos] as usize;
    let score = score - scores.row_max[query_residue];
    let mut letters = scores.letters_at_least(query_residue, floor(score));
    if current_pos == word_length - 1 {
        while letters != 0 {
            let residue = letters.trailing_zeros();
            letters &= letters - 1;
            out.push((prefix << X_WIDE_CHARSIZE) | residue);
        }
        return;
    }
    let row = &scores.rows[query_residue];
    if current_pos + 2 == word_length {
        // The last level for each letter here: its `score -= row_max[...]` and loop.
        let last = query_word[current_pos + 1] as usize;
        let last_max = scores.row_max[last];
        while letters != 0 {
            let residue = letters.trailing_zeros();
            letters &= letters - 1;
            let at_last = score + row[residue as usize] - last_max;
            let mut last_letters = scores.letters_at_least(last, floor(at_last));
            let last_prefix = ((prefix << X_WIDE_CHARSIZE) | residue) << X_WIDE_CHARSIZE;
            while last_letters != 0 {
                out.push(last_prefix | last_letters.trailing_zeros());
                last_letters &= last_letters - 1;
            }
        }
        return;
    }
    while letters != 0 {
        let residue = letters.trailing_zeros();
        letters &= letters - 1;
        x_wide_neighbors(
            query_word,
            word_length,
            scores,
            threshold,
            score + row[residue as usize],
            current_pos + 1,
            (prefix << X_WIDE_CHARSIZE) | residue,
            out,
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:490-496,504-509,518-519,539-543
// The statements of `x_target_list`, appending T(s) to `out` in ascending order of the cells: the
// exact cell, which s_AddWordHits adds first, goes where its index falls among the neighbours (it
// is none of them: they score at least the threshold, the word itself less), which changes no table
// (the cells of T(s) are distinct; module comment, step 2).
/// Appends T(s) of the word at `query_word` to `out`, its cells ascending.
fn x_wide_targets(
    query_word: &[u8],
    word_length: usize,
    scores: &XWideScores,
    threshold: i32,
    out: &mut Vec<u32>,
) {
    let first = query_word[0] as usize;
    let start = out.len();
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
    let mut self_score = scores.rows[first][first];
    for &residue in &query_word[1..word_length] {
        let residue = residue as usize;
        self_score += scores.rows[residue][residue];
    }
    if threshold == 0 || self_score < threshold {
        let exact = x_wide_index(query_word, word_length);
        let at = start + out[start..].partition_point(|&w| w < exact);
        out.insert(at, exact);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:580-582,600-604 (cited at
// `x_wide_neighbors`)
// No NCBI counterpart for the count itself: `x_wide_neighbors`' words, counted (the same tests,
// the last level's letters by their number), for LOSAT's choice of the blocks' ends.
/// The number of words `x_wide_neighbors` appends for `query_word` from `current_pos` on.
fn x_wide_neighbor_count(
    query_word: &[u8],
    word_length: usize,
    scores: &XWideScores,
    threshold: i32,
    score: i32,
    current_pos: usize,
) -> u64 {
    let query_residue = query_word[current_pos] as usize;
    let score = score - scores.row_max[query_residue];
    let mut letters =
        scores.letters_at_least(query_residue, i64::from(threshold) - i64::from(score));
    if current_pos == word_length - 1 {
        return u64::from(letters.count_ones());
    }
    let row = &scores.rows[query_residue];
    let mut words = 0u64;
    while letters != 0 {
        let residue = letters.trailing_zeros() as usize;
        letters &= letters - 1;
        words += x_wide_neighbor_count(
            query_word,
            word_length,
            scores,
            threshold,
            score + row[residue],
            current_pos + 1,
        );
    }
    words
}

// No NCBI counterpart: the longest word `XWideWordCount` counts (the joint builder's, 32-bit
// indices of 5-bit letters).
const X_WIDE_COUNT_MOST_LETTERS: usize = 6;

// No NCBI counterpart: counts the cells of `x_wide_targets` without a table, for LOSAT's choice of
// the blocks' ends (EXPERIMENT LOSAT_X_BPQBLOCK, `blastp/blast_engine.rs` `x_word_blocks`); it
// changes no value NCBI computes.
/// The number of cells T(s) of a query word, the appends of `s_AddWordHits` for one query offset
/// to the lookup table (`x_wide_targets`, counted: the neighbours, and the word itself when there
/// is no threshold or its own score is below it), for every word length of the lookup table, the
/// direct builders' (2 and 3) and the joint builder's (4). Neither the number nor the word's own
/// score depends on the order of its letters, so a count is kept per set of letters (`counts`, at
/// the set's rank among the sets of `word_length` letters), and the letters of a proteome take a
/// few tens of thousands of enumerations; as `CompressedWordCount` (compressed.rs) for the
/// compressed table. The kept counts are shared: threads that count with one counter keep them for
/// each other (two threads may count a set both; they keep the same count).
pub(crate) struct XWideWordCount {
    scores: XWideScores,
    threshold: i32,
    word_length: usize,
    /// `binomial[n][k]` is C(n, k), for n below the alphabet's size plus `word_length` and k up
    /// to `word_length`.
    binomial: Vec<[u32; X_WIDE_COUNT_MOST_LETTERS + 1]>,
    /// Per rank of a set of letters, its count plus 1; 0 while it is not counted.
    counts: Box<[std::sync::atomic::AtomicU32]>,
}

impl XWideWordCount {
    /// The counts of the table of `matrix`'s words of `word_length` letters with the threshold
    /// `threshold` (the table's `Int4` threshold, as `blastp_aa_lookup` builds it); `None` for a
    /// word length the joint builder's indices do not cover.
    pub(crate) fn new(matrix: ScoringMatrix, word_length: usize, threshold: i32) -> Option<Self> {
        if !(1..=X_WIDE_COUNT_MOST_LETTERS).contains(&word_length) || !x_wide_covers(word_length) {
            return None;
        }
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        // The row maxima of `build_lookup_from_prepared` (`XPairLutQuery::x_from_prepared`).
        let row_max: Vec<i32> = (0..alphabet_size)
            .map(|i| {
                (0..alphabet_size)
                    .map(|j| lookup_matrix_score(matrix, i, j))
                    .max()
                    .unwrap_or(i32::MIN)
            })
            .collect();
        // A set of `word_length` of the letters, sorted, a_0 <= ... <= a_{w-1}, has the rank sum
        // C(a_i + i, i + 1), one of C(letters + w - 1, w) (as `CompressedWordCount`).
        let mut binomial = vec![[0u32; X_WIDE_COUNT_MOST_LETTERS + 1]; alphabet_size + word_length];
        for n in 0..binomial.len() {
            binomial[n][0] = 1;
            for k in 1..=word_length.min(n) {
                binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
            }
        }
        let sets = binomial[alphabet_size + word_length - 1][word_length] as usize;
        Some(Self {
            scores: XWideScores::new(matrix, &row_max),
            threshold,
            word_length,
            binomial,
            counts: (0..sets)
                .map(|_| std::sync::atomic::AtomicU32::new(0))
                .collect(),
        })
    }

    /// The cells of the table's backbone (`compute_backbone_size`, as `build_lookup_from_prepared`).
    pub(crate) fn backbone_cells(&self) -> u64 {
        compute_backbone_size(
            self.word_length,
            LOOKUP_ALPHABET_SIZE,
            ilog2(LOOKUP_ALPHABET_SIZE) + 1,
        ) as u64
    }

    /// The rank of the set of `word`'s letters; `None` for a letter outside the alphabet.
    fn rank(&self, word: &[u8]) -> Option<usize> {
        let mut sorted = [0u8; X_WIDE_COUNT_MOST_LETTERS];
        for (to, &residue) in sorted.iter_mut().zip(word) {
            if usize::from(residue) >= LOOKUP_ALPHABET_SIZE {
                return None;
            }
            *to = residue;
        }
        let sorted = &mut sorted[..self.word_length];
        sorted.sort_unstable();
        Some(
            sorted
                .iter()
                .enumerate()
                .map(|(i, &letter)| self.binomial[usize::from(letter) + i][i + 1] as usize)
                .sum(),
        )
    }

    /// The cells T(s) of `word` (NCBISTDAA letters, `word_length` of them); 0 for a word with a
    /// letter outside the alphabet.
    fn words(&self, word: &[u8]) -> u64 {
        use std::sync::atomic::Ordering;
        let Some(rank) = self.rank(word) else {
            return 0;
        };
        let kept = &self.counts[rank];
        match kept.load(Ordering::Relaxed) {
            0 => {}
            words_and_one => return u64::from(words_and_one - 1),
        }
        let words = self.count(word);
        if let Some(words_and_one) = u32::try_from(words).ok().and_then(|w| w.checked_add(1)) {
            kept.store(words_and_one, Ordering::Relaxed);
        }
        words
    }

    /// `words` without the kept counts: the statements of `x_wide_targets`, counted.
    fn count(&self, word: &[u8]) -> u64 {
        let scores = &self.scores;
        let first = usize::from(word[0]);
        let mut self_score = scores.rows[first][first];
        for &residue in &word[1..self.word_length] {
            let residue = usize::from(residue);
            self_score += scores.rows[residue][residue];
        }
        let own = u64::from(self.threshold == 0 || self_score < self.threshold);
        if self.threshold == 0 {
            return own;
        }
        let mut max_score = scores.row_max[first];
        for &residue in &word[1..self.word_length] {
            max_score += scores.row_max[usize::from(residue)];
        }
        own + x_wide_neighbor_count(word, self.word_length, scores, self.threshold, max_score, 0)
    }

    /// The cells of `queries`' words (each query's NCBISTDAA letters): the most that a table whose
    /// words are some of them appends (SEG and `-lcase_masking` only leave words out).
    pub(crate) fn batch_words(&self, queries: &mut dyn Iterator<Item = &[u8]>) -> u64 {
        let mut words = 0u64;
        for query in queries {
            for word in query.windows(self.word_length) {
                words = words.saturating_add(self.words(word));
            }
        }
        words
    }
}

/// No NCBI counterpart: table t's part of one source (step 1): its offsets E_t(s),
/// `offsets[from..to]`.
#[derive(Clone, Copy)]
struct XWideRun {
    t: u32,
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
    /// The runs of source i, `runs[run_starts[i]..run_starts[i + 1]]`: one for each table with an
    /// offset of the source, tables ascending.
    runs: Vec<XWideRun>,
    run_starts: Vec<u32>,
}

// No NCBI counterpart (the form of a stable sort): the order of `keys`, each below `2^bits`,
// ascending, equal keys in their order: a least-significant-digit radix sort, digits of at most
// `X_WIDE_RADIX_BITS` bits.
/// The positions of `keys` sorted by key, stably.
fn x_wide_radix_order(keys: &[u32], bits: usize) -> Vec<u32> {
    debug_assert!(bits >= 32 || keys.iter().all(|&key| key >> bits == 0));
    let passes = bits.div_ceil(X_WIDE_RADIX_BITS).max(1);
    let digit = bits.div_ceil(passes);
    let mask = (1u32 << digit) - 1;
    let mut order: Vec<u32> = (0..keys.len() as u32).collect();
    let mut next: Vec<u32> = vec![0; keys.len()];
    let mut at: Vec<usize> = vec![0; 1 << digit];
    for pass in 0..passes {
        let shift = pass * digit;
        at.fill(0);
        for &key in keys {
            at[((key >> shift) & mask) as usize] += 1;
        }
        let mut sum = 0usize;
        for slot in at.iter_mut() {
            sum += std::mem::replace(slot, sum);
        }
        for &k in &order {
            let d = ((keys[k as usize] >> shift) & mask) as usize;
            next[at[d]] = k;
            at[d] += 1;
        }
        std::mem::swap(&mut order, &mut next);
    }
    order
}

impl XWideSources {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-456 (cited at
    // `XPairLutQuery::x_build_group`): the exact backbone of each table's locations. No NCBI
    // counterpart for the form (a stable sort of the words by cell).
    /// The sources of the tables of `tables` (step 1).
    fn new(query: &[u8], tables: &[Vec<(i32, i32)>], word_length: usize) -> Self {
        assert!(
            u32::try_from(tables.len()).is_ok(),
            "joint lookup build: too many tables"
        );
        // The words in the order of the loops: table by table, each in its loop's order.
        let (mut keys, mut words, mut table_of) = (Vec::new(), Vec::new(), Vec::new());
        for (t, locations) in tables.iter().enumerate() {
            x_wide_exact_words(query, locations, word_length, |idx, offset| {
                keys.push(idx);
                words.push(offset);
                table_of.push(t as u32);
            });
        }
        assert!(
            keys.len() < u32::MAX as usize,
            "joint lookup build: query words must fit 32 bits"
        );
        let order = x_wide_radix_order(&keys, word_length * X_WIDE_CHARSIZE);
        let mut cells: Vec<u32> = Vec::new();
        let mut offsets: Vec<i32> = Vec::with_capacity(order.len());
        let mut runs: Vec<XWideRun> = Vec::new();
        let mut run_starts: Vec<u32> = Vec::new();
        for (at, &k) in order.iter().enumerate() {
            let (k, at) = (k as usize, at as u32);
            let (cell, t) = (keys[k], table_of[k]);
            let new_source = cells.last() != Some(&cell);
            if new_source {
                cells.push(cell);
                run_starts.push(runs.len() as u32);
            }
            // Equal cells keep the loops' order, so a source's offsets come table by table.
            match runs.last_mut() {
                Some(run) if !new_source && run.t == t => run.to += 1,
                _ => runs.push(XWideRun {
                    t,
                    from: at,
                    to: at + 1,
                }),
            }
            offsets.push(words[k]);
        }
        run_starts.push(runs.len() as u32);
        Self {
            cells,
            offsets,
            runs,
            run_starts,
        }
    }

    /// No NCBI counterpart: the runs of source `i`.
    #[inline(always)]
    fn runs_of(&self, i: usize) -> &[XWideRun] {
        &self.runs[self.run_starts[i] as usize..self.run_starts[i + 1] as usize]
    }

    /// No NCBI counterpart: the offsets of `run` (E_t(s) of its table).
    #[inline(always)]
    fn offsets(&self, run: XWideRun) -> &[i32] {
        &self.offsets[run.from as usize..run.to as usize]
    }
}

/// No NCBI counterpart: the stripes of the backbone (step 2): `1 << bits` cells each (`bits` from 5
/// to 16: a stripe's presence-vector words are its own, and a cell's index in it is a `u16`), `count`
/// of them.
#[derive(Clone, Copy)]
struct XWideStripes {
    bits: u32,
    count: usize,
}

/// No NCBI counterpart: a source of a range with cells in a stripe (step 2): its index in the
/// range, and the end of its cells in the range's `cells`.
#[derive(Clone, Copy, Default)]
struct XWideSpan {
    i: u32,
    end: u32,
}

/// No NCBI counterpart: a range's lists T(s) by stripe (step 2). For stripe j, the cells (their
/// index in the stripe) of the range's lists, source by source ascending, are
/// `cells[cell_starts[j]..cell_starts[j + 1]]`, and the sources with cells in the stripe, ascending,
/// with where their cells end, `spans[span_starts[j]..span_starts[j + 1]]`. Also, by table, the
/// least and greatest offset of the range's sources with a target.
struct XWideListPart {
    cells: Vec<u16>,
    cell_starts: Vec<u32>,
    spans: Vec<XWideSpan>,
    span_starts: Vec<u32>,
    q_bounds: Vec<(u32, u32)>,
}

/// No NCBI counterpart: the lists T(s) of the sources (step 2), by ranges of sources: range r holds
/// the sources `bounds[r]..bounds[r + 1]`.
struct XWideLists {
    bounds: Vec<usize>,
    parts: Vec<XWideListPart>,
}

impl XWideLists {
    /// No NCBI counterpart: `visit(runs of s, cells of T(s) in stripe j)` for each source s with
    /// cells in stripe `j`, ascending (the ranges in order).
    #[inline(always)]
    fn for_each_in(
        &self,
        sources: &XWideSources,
        j: usize,
        mut visit: impl FnMut(&[XWideRun], &[u16]),
    ) {
        for (part, &base) in self.parts.iter().zip(&self.bounds) {
            let spans = &part.spans[part.span_starts[j] as usize..part.span_starts[j + 1] as usize];
            let mut from = part.cell_starts[j] as usize;
            for span in spans {
                let end = span.end as usize;
                visit(
                    sources.runs_of(base + span.i as usize),
                    &part.cells[from..end],
                );
                from = end;
            }
        }
    }
}

/// No NCBI counterpart: `visit(k, j, block)` for each list k of `targets` (list k ends at
/// `ends[k]`) and each stripe j it reaches, with its cells in stripe j: a list is ascending
/// (`x_wide_targets`), so they are one block.
#[inline(always)]
fn x_wide_blocks(
    targets: &[u32],
    ends: &[usize],
    stripes: XWideStripes,
    mut visit: impl FnMut(usize, usize, &[u32]),
) {
    let mut from = 0usize;
    for (k, &end) in ends.iter().enumerate() {
        let mut list = &targets[from..end];
        while let Some(&w) = list.first() {
            let j = (w >> stripes.bits) as usize;
            let n = list.partition_point(|&w| (w >> stripes.bits) as usize == j);
            visit(k, j, &list[..n]);
            list = &list[n..];
        }
        from = end;
    }
}

/// No NCBI counterpart (the form of step 2): the lists of the sources `range` (each T(s)
/// ascending, `x_wide_targets`) by stripe, and their offsets' bounds by table (of `groups`).
#[allow(clippy::too_many_arguments)]
fn x_wide_list_part(
    range: std::ops::Range<usize>,
    sources: &XWideSources,
    query: &[u8],
    word_length: usize,
    scores: &XWideScores,
    threshold: i32,
    stripes: XWideStripes,
    groups: usize,
) -> XWideListPart {
    let mut targets: Vec<u32> = Vec::new();
    let mut ends: Vec<usize> = Vec::with_capacity(range.len());
    let mut q_bounds = vec![(u32::MAX, 0u32); groups];
    for i in range {
        let runs = sources.runs_of(i);
        // blast_aalookup.c:490: w = query + offset_list[2];
        let first = usize::try_from(sources.offsets[runs[0].from as usize])
            .expect("NCBI BLAST query offset must fit");
        let before = targets.len();
        x_wide_targets(
            &query[first..],
            word_length,
            scores,
            threshold,
            &mut targets,
        );
        ends.push(targets.len());
        if targets.len() > before {
            for &run in runs {
                let (q_min, q_max) = &mut q_bounds[run.t as usize];
                for &q in sources.offsets(run) {
                    *q_min = (*q_min).min(q as u32);
                    *q_max = (*q_max).max(q as u32);
                }
            }
        }
    }
    assert!(
        u32::try_from(targets.len()).is_ok(),
        "joint lookup build: a range's lists must fit 32 bits"
    );
    // The sizes by stripe, then the cells and spans in place.
    let mut cell_starts: Vec<u32> = vec![0; stripes.count + 1];
    let mut span_starts: Vec<u32> = vec![0; stripes.count + 1];
    x_wide_blocks(&targets, &ends, stripes, |_, j, block| {
        cell_starts[j + 1] += block.len() as u32;
        span_starts[j + 1] += 1;
    });
    for j in 0..stripes.count {
        cell_starts[j + 1] += cell_starts[j];
        span_starts[j + 1] += span_starts[j];
    }
    let mut cells: Vec<u16> = vec![0; targets.len()];
    let mut spans: Vec<XWideSpan> = vec![XWideSpan::default(); span_starts[stripes.count] as usize];
    let mut cell_at: Vec<u32> = cell_starts[..stripes.count].to_vec();
    let mut span_at: Vec<u32> = span_starts[..stripes.count].to_vec();
    let mask = (1u32 << stripes.bits) - 1;
    x_wide_blocks(&targets, &ends, stripes, |k, j, block| {
        let at = cell_at[j] as usize;
        for (cell, &w) in cells[at..at + block.len()].iter_mut().zip(block) {
            *cell = (w & mask) as u16;
        }
        cell_at[j] += block.len() as u32;
        spans[span_at[j] as usize] = XWideSpan {
            i: k as u32,
            end: cell_at[j],
        };
        span_at[j] += 1;
    });
    XWideListPart {
        cells,
        cell_starts,
        spans,
        span_starts,
        q_bounds,
    }
}

// No NCBI counterpart (scheduling only): where a build's tasks run: on TBLASTX's search pool (the
// unit path's group tables), or, for BLASTX's table of 4-letter words (item SC9 of line NP), on the
// pool the build already runs on, in as many tasks as BLASTX's batch gives its build
// (`LOSAT_X_BXBUILDPAR`, `blastx/x_bx_mem.rs`).
#[derive(Clone, Copy)]
enum XWideOn<'p, 'a> {
    Pool(&'p crate::utils::threading::SearchPool<'a>),
    Here(usize),
}

// No NCBI counterpart (scheduling only): NCBI builds a table on one thread. `f` on each of `parts`,
// the results in part order: on tasks of the pool (`x_claimed`) when the build has two threads or
// more, else in order on this thread.
fn x_wide_parts<T: Send, R: Send>(
    on: XWideOn<'_, '_>,
    parts: Vec<T>,
    f: &(dyn Fn(T) -> R + Sync),
) -> Vec<R> {
    // Without threads every part runs on this thread, wherever the build was to run.
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let _ = on;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parts.len() > 1 && x_wide_threads(on) > 1 {
        return match on {
            XWideOn::Pool(pool) => pool.install(|| x_claimed(parts, f)),
            XWideOn::Here(_) => x_claimed(parts, f),
        };
    }
    parts.into_iter().map(f).collect()
}

// No NCBI counterpart (scheduling only): the threads the build's tasks run on.
fn x_wide_threads(on: XWideOn<'_, '_>) -> usize {
    if !cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )) {
        return 1;
    }
    match on {
        XWideOn::Pool(pool) if pool.enabled() => pool.threads().max(1),
        XWideOn::Pool(_) => 1,
        XWideOn::Here(parts) => parts.max(1),
    }
}

/// No NCBI counterpart: `written` as initialized elements.
///
/// # Safety
/// Every element of `written` was written.
unsafe fn x_wide_written<T>(written: &mut [std::mem::MaybeUninit<T>]) -> &mut [T] {
    // SAFETY: (caller) every element is initialized; `MaybeUninit<T>` has `T`'s layout.
    unsafe { std::slice::from_raw_parts_mut(written.as_mut_ptr().cast::<T>(), written.len()) }
}

/// No NCBI counterpart: a table's totals in a stripe of cells (step 3 (a)): overflow entries,
/// longest count, cells with entries.
#[derive(Clone, Copy, Default)]
struct XWideStripe {
    overflow: usize,
    longest: u32,
    used: usize,
}

/// No NCBI counterpart: a table being built (step 3 (b)): its vectors, the lengths they take, and
/// its totals.
struct XWideTable {
    backbone: Vec<BackboneCell>,
    overflow: Vec<i32>,
    pv: Vec<u32>,
    backbone_len: usize,
    overflow_size: usize,
    longest: u32,
    compact: bool,
}

/// No NCBI counterpart: a table's part of a stripe of cells (step 3 (b)): its cells (all the
/// stripe's, or, compacted, the stripe's cells with entries), its part of the overflow (cells
/// ascending) and where that starts in the overflow, and the stripe's presence-vector words.
struct XWidePlace<'a> {
    cells: &'a mut [std::mem::MaybeUninit<BackboneCell>],
    overflow: &'a mut [std::mem::MaybeUninit<i32>],
    cursor_base: usize,
    pv: &'a mut [u32],
    compact: bool,
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
        self.x_build_tables_wide_with(tables, XWideOn::Pool(pool), mapped, None, true)
    }

    /// The lookup table of 4-letter words (`-word_size 4`) of BLASTX (item SC9 of line NP), TBLASTN
    /// and TBLASTX's reference frame loop (item X2e of line RG), whose backbone index the 16-bit
    /// direct builder does not cover: the table of one list of lookup locations, all of the query's,
    /// which is `build_lookup_from_prepared`'s table (module comment: the full table filtered to its
    /// own locations), in `BlastAaLookupFinalize`'s layout and never compacted (their scans index the
    /// whole backbone), built in tasks for `parts` threads of the pool the build runs on.
    pub(crate) fn x_build_full_wide(&self, parts: usize) -> BlastAaLookupTable {
        self.x_build_full_wide_with(XWideOn::Here(parts))
    }

    /// BLASTP's lookup table of 4-letter words (`-word_size 4`, item SC8b of line NP): the table
    /// [`Self::x_build_full_wide`] builds, in tasks on `pool`'s threads.
    pub(crate) fn x_build_full_wide_on(
        &self,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> BlastAaLookupTable {
        self.x_build_full_wide_with(XWideOn::Pool(pool))
    }

    // No NCBI counterpart (scheduling only): the one table of the query's locations, its tasks on
    // `on`.
    /// The body of [`Self::x_build_full_wide`] and [`Self::x_build_full_wide_on`].
    fn x_build_full_wide_with(&self, on: XWideOn<'_, '_>) -> BlastAaLookupTable {
        let tables = [self.prepared.lookup_locations.clone()];
        let mut built = self.x_build_tables_wide_with(&tables, on, false, None, false);
        debug_assert!(built.len() == 1 && built[0].rank.is_empty());
        built.pop().expect("one table").table
    }

    // No NCBI counterpart (scheduling only): `x_build_tables_wide` with the stripes' size as an
    // argument (the unit tests take stripes of 32 cells: as many stripes as the backbone holds).
    /// `x_build_tables_wide` with stripes of `2^stripe_bits` cells where given, the tasks on `on`,
    /// and tables compacted only where `compactable`.
    fn x_build_tables_wide_with(
        &self,
        tables: &[Vec<(i32, i32)>],
        on: XWideOn<'_, '_>,
        mapped: bool,
        stripe_bits: Option<u32>,
        compactable: bool,
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
        let threads = x_wide_threads(on);
        let groups = tables.len();

        // 1. The sources.
        let sources = XWideSources::new(query, tables, word_length);

        // The stripes: about `X_WIDE_STRIPES_PER_THREAD` a thread, of 32 to 2^16 cells.
        let stripes = {
            let cells = (backbone_size / (threads * X_WIDE_STRIPES_PER_THREAD))
                .clamp(1 << PV_ARRAY_BTS, 1 << X_WIDE_STRIPE_BITS_MOST);
            let bits = stripe_bits.unwrap_or(cells.ilog2());
            assert!(
                (PV_ARRAY_BTS as u32..=X_WIDE_STRIPE_BITS_MOST).contains(&bits),
                "joint lookup build: stripes of 2^{bits} cells"
            );
            XWideStripes {
                bits,
                count: backbone_size.div_ceil(1 << bits),
            }
        };
        let stripe_len =
            |j: usize| (1usize << stripes.bits).min(backbone_size - (j << stripes.bits));

        // 2. The lists T(s) by stripe, once for every table.
        let scores = XWideScores::new(self.matrix, &self.row_max);
        let n = sources.cells.len();
        let ranges = (threads * X_WIDE_LIST_PARTS_PER_THREAD).min(n).max(1);
        let bounds: Vec<usize> = (0..=ranges).map(|r| n * r / ranges).collect();
        let threshold = self.threshold;
        let made = x_wide_parts(
            on,
            bounds.windows(2).map(|b| b[0]..b[1]).collect(),
            &|range: std::ops::Range<usize>| {
                x_wide_list_part(
                    range,
                    &sources,
                    query,
                    word_length,
                    &scores,
                    threshold,
                    stripes,
                    groups,
                )
            },
        );
        let lists = XWideLists {
            bounds,
            parts: made,
        };
        let mut q_bounds = vec![(u32::MAX, 0u32); groups];
        for part in &lists.parts {
            for (bound, &(q_min, q_max)) in q_bounds.iter_mut().zip(&part.q_bounds) {
                bound.0 = bound.0.min(q_min);
                bound.1 = bound.1.max(q_max);
            }
        }

        // 3 (a). Each stripe's counts: its cells' in each table (`groups` blocks of
        // `stripe_len(j)`), the stripes' blocks in stripe order; and each table's totals in it.
        let counted_len = groups * backbone_size;
        let mut counts: Vec<u32> = Vec::with_capacity(counted_len);
        let mut blocks: Vec<(usize, &mut [std::mem::MaybeUninit<u32>])> =
            Vec::with_capacity(stripes.count);
        {
            let mut left = &mut counts.spare_capacity_mut()[..counted_len];
            for j in 0..stripes.count {
                let (block, rest) = std::mem::take(&mut left).split_at_mut(groups * stripe_len(j));
                blocks.push((j, block));
                left = rest;
            }
            debug_assert!(left.is_empty());
        }
        let totals: Vec<Vec<XWideStripe>> = x_wide_parts(on, blocks, &|(j, block): (
            usize,
            &mut [std::mem::MaybeUninit<u32>],
        )| {
            let len = stripe_len(j);
            block.fill(std::mem::MaybeUninit::new(0));
            // SAFETY: every element of `block` was written just above.
            let block = unsafe { x_wide_written(block) };
            lists.for_each_in(&sources, j, |runs, cells| {
                for run in runs {
                    let k = run.to - run.from;
                    let counts = &mut block[run.t as usize * len..][..len];
                    for &w in cells {
                        counts[usize::from(w)] += k;
                    }
                }
            });
            block
                .chunks(len)
                .map(|counts| {
                    let mut totals = XWideStripe::default();
                    for &count in counts {
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
                })
                .collect()
        });
        // SAFETY: the stripes' blocks cover the `counted_len` counts (each stripe `groups` times its
        // cells), and each stripe wrote all of its block.
        unsafe { counts.set_len(counted_len) };

        // 3 (b). Each table: its size and layout (BlastAaLookupFinalize's rule, compacted where
        // the branch-free scan takes the table: `x_scan_bf_takes`, as `x_build_group` decides),
        // and its parts of the stripes, cells ascending.
        let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
        let mut built: Vec<XWideTable> = (0..groups)
            .map(|t| {
                let (mut overflow_size, mut used, mut longest) = (0usize, 0usize, 0u32);
                for stripe in &totals {
                    overflow_size += stripe[t].overflow;
                    used += stripe[t].used;
                    longest = longest.max(stripe[t].longest);
                }
                // The positions in a stripe's part of the overflow are `u32`s beside the flag (2^31
                // entries of a stripe of at most 2^16 cells are 8 GiB); a cell's
                // cursor takes 64 bits (`set_cursor`), so the table's overflow passes 2^31 entries,
                // where NCBI's `Int4` cursor overflows, as the other builders' do.
                assert!(
                    totals.iter().all(|stripe| stripe[t].overflow < X_WIDE_IN_OVERFLOW as usize),
                    "a stripe of a lookup table of 4-letter words holds at most 2^31 - 1 overflow entries"
                );
                let compact = compactable
                    && crate::algorithm::tblastx::x_scan_bf::x_scan_bf_takes(used, backbone_size);
                let backbone_len = if compact { used } else { backbone_size };
                XWideTable {
                    backbone: Vec::with_capacity(x_group_capacity::<BackboneCell>(
                        backbone_len,
                        mapped,
                    )),
                    overflow: Vec::with_capacity(x_group_capacity::<i32>(overflow_size, mapped)),
                    pv: vec![0u32; pv_size],
                    backbone_len,
                    overflow_size,
                    longest,
                    compact,
                }
            })
            .collect();
        let mut places: Vec<Vec<XWidePlace<'_>>> = (0..stripes.count)
            .map(|_| Vec::with_capacity(groups))
            .collect();
        for (t, table) in built.iter_mut().enumerate() {
            let mut cells_left = &mut table.backbone.spare_capacity_mut()[..table.backbone_len];
            let mut overflow_left = &mut table.overflow.spare_capacity_mut()[..table.overflow_size];
            let mut pv_left = &mut table.pv[..];
            let mut cursor_base = 0usize;
            for (j, (stripe, place)) in totals.iter().zip(places.iter_mut()).enumerate() {
                let len = stripe_len(j);
                let kept = if table.compact { stripe[t].used } else { len };
                let (cells, rest) = std::mem::take(&mut cells_left).split_at_mut(kept);
                cells_left = rest;
                let (overflow, rest) =
                    std::mem::take(&mut overflow_left).split_at_mut(stripe[t].overflow);
                overflow_left = rest;
                // A stripe starts at a multiple of 32 cells; the last takes the words left.
                let words = len.div_ceil(1 << PV_ARRAY_BTS).min(pv_left.len());
                let (pv, rest) = std::mem::take(&mut pv_left).split_at_mut(words);
                pv_left = rest;
                place.push(XWidePlace {
                    cells,
                    overflow,
                    cursor_base,
                    pv,
                    compact: table.compact,
                });
                cursor_base += stripe[t].overflow;
            }
            debug_assert!(cells_left.is_empty() && overflow_left.is_empty());
            // The backbone's cells are 28 more than a multiple of 32 at 1 to 4 letters, so the last
            // stripe takes the extra word (`pv_size`); none is left past 32 cells either way.
            debug_assert!(pv_left.len() <= 1);
        }

        // 3 (c). Each stripe: its cells and each cell's first position (from its counts), then the
        // appends.
        let mut stripe_counts: Vec<&mut [u32]> = Vec::with_capacity(stripes.count);
        {
            let mut left = &mut counts[..];
            for j in 0..stripes.count {
                let (block, rest) = std::mem::take(&mut left).split_at_mut(groups * stripe_len(j));
                stripe_counts.push(block);
                left = rest;
            }
        }
        x_wide_parts(
            on,
            places.into_iter().zip(stripe_counts).enumerate().collect(),
            &|(j, (place, counts)): (usize, (Vec<XWidePlace<'_>>, &mut [u32]))| {
                let len = stripe_len(j);
                let mut cells: Vec<&mut [BackboneCell]> = Vec::with_capacity(groups);
                let mut overflows: Vec<&mut [std::mem::MaybeUninit<i32>]> =
                    Vec::with_capacity(groups);
                let mut compacts: Vec<bool> = Vec::with_capacity(groups);
                for (place, counts) in place.into_iter().zip(counts.chunks_mut(len)) {
                    let mut kept = 0usize;
                    let mut cursor = 0usize;
                    for (local, at) in counts.iter_mut().enumerate() {
                        let count = *at;
                        if count == 0 {
                            if !place.compact {
                                place.cells[local].write(BackboneCell::default());
                            }
                            continue;
                        }
                        // pv_set at the stripe's first cell (a multiple of 32) + local.
                        place.pv[local >> PV_ARRAY_BTS] |=
                            1u32 << (local & crate::algorithm::tblastx::lookup::PV_ARRAY_MASK);
                        let slot = if place.compact {
                            kept += 1;
                            kept - 1
                        } else {
                            local
                        };
                        let mut cell = BackboneCell {
                            num_used: count as i32,
                            ..BackboneCell::default()
                        };
                        *at = if count as usize > AA_HITS_PER_CELL {
                            cell.set_cursor(place.cursor_base + cursor);
                            cursor += count as usize;
                            (cursor - count as usize) as u32 | X_WIDE_IN_OVERFLOW
                        } else {
                            (slot as u32) << 2
                        };
                        place.cells[slot].write(cell);
                    }
                    // The cells written and the overflow taken are the stripe's parts (the
                    // casts below and `set_len` rely on them).
                    assert!(place.compact || place.cells.len() == len);
                    assert!(!place.compact || kept == place.cells.len());
                    assert_eq!(cursor, place.overflow.len());
                    // SAFETY: the loop wrote every cell of `place.cells` (all the stripe's, or,
                    // compacted, its cells with entries, `kept` of them, as asserted).
                    cells.push(unsafe { x_wide_written(place.cells) });
                    overflows.push(place.overflow);
                    compacts.push(place.compact);
                }
                // The appends: a cell's position advances by the offsets written, from the
                // cell's first position, over the walk the counts took.
                lists.for_each_in(&sources, j, |runs, targets| {
                    for &run in runs {
                        let t = run.t as usize;
                        let src = sources.offsets(run);
                        let (cells, overflow) = (&mut *cells[t], &mut *overflows[t]);
                        let positions = &mut counts[t * len..][..len];
                        // One offset (most sources in a table) apart from several.
                        if let [q] = *src {
                            for &w in targets {
                                let at = &mut positions[usize::from(w)];
                                if *at & X_WIDE_IN_OVERFLOW == 0 {
                                    cells[(*at >> 2) as usize].entries[(*at & 3) as usize] = q;
                                } else {
                                    overflow[(*at & !X_WIDE_IN_OVERFLOW) as usize].write(q);
                                }
                                *at += 1;
                            }
                            continue;
                        }
                        for &w in targets {
                            let at = &mut positions[usize::from(w)];
                            if *at & X_WIDE_IN_OVERFLOW == 0 {
                                let (slot, entry) = ((*at >> 2) as usize, (*at & 3) as usize);
                                for (k, &q) in src.iter().enumerate() {
                                    cells[slot].entries[entry + k] = q;
                                }
                            } else {
                                let from = (*at & !X_WIDE_IN_OVERFLOW) as usize;
                                for (k, &q) in src.iter().enumerate() {
                                    overflow[from + k].write(q);
                                }
                            }
                            *at += src.len() as u32;
                        }
                    }
                });
                // Each cell with entries took its count of offsets from its first position (its
                // slot's entries, or its run of the overflow from the cursor the stripe's part
                // gave it, cells ascending), so the overflow entries were written, once each.
                for (t, ends) in counts.chunks(len).enumerate() {
                    let (mut kept, mut overflow) = (0usize, 0usize);
                    for (local, &end) in ends.iter().enumerate() {
                        if end == 0 {
                            continue;
                        }
                        let slot = if compacts[t] {
                            kept += 1;
                            kept - 1
                        } else {
                            local
                        };
                        let count = cells[t][slot].num_used as usize;
                        if count > AA_HITS_PER_CELL {
                            overflow += count;
                            assert_eq!(end, overflow as u32 | X_WIDE_IN_OVERFLOW);
                        } else {
                            assert_eq!(end, ((slot as u32) << 2) + count as u32);
                        }
                    }
                    assert_eq!(overflow, overflows[t].len());
                }
            },
        );

        built
            .into_iter()
            .zip(q_bounds)
            .map(|(mut table, (q_min, q_max))| {
                // SAFETY: 3 (c) wrote every cell of the backbone (each stripe its part: all its
                // cells, or, compacted, its cells with entries, which the totals counted) and
                // every entry of the overflow (each stripe's part, as it asserted after its
                // appends; the parts cover the overflow).
                unsafe {
                    table.backbone.set_len(table.backbone_len);
                    table.overflow.set_len(table.overflow_size);
                }
                let longest_chain = i32::try_from(table.longest)
                    .expect("NCBI BLAST lookup chain hit count must fit Int4");
                let rank = if table.compact {
                    crate::algorithm::tblastx::x_scan_bf::x_rank(&table.pv)
                } else {
                    Vec::new()
                };
                XPairLutGroup {
                    q_min: if longest_chain > 0 { q_min } else { 0 },
                    q_max,
                    rank,
                    table: BlastAaLookupTable {
                        backbone: table.backbone,
                        overflow: table.overflow,
                        pv: table.pv,
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
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "x_pair_lut_wide_tests.rs"]
mod tests;
